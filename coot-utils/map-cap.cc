/*
 * coot-utils/map-cap.cc
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include <glm/glm.hpp>

#include "utils/colour-holder.hh"
#include "coot-map-utils.hh"
#include "map-cap.hh"

namespace {

   // 2D grid for storing density values and coordinates on the cap plane
   class coord_array_2d {
   public:
      std::vector<std::pair<clipper::Coord_orth, float>> v;
      unsigned int nx;
      coord_array_2d(unsigned int nx_in, unsigned int ny) : nx(nx_in) {
         v.resize(nx * ny);
      }
      void set(const unsigned int &i, const unsigned int &j, const clipper::Coord_orth &co, const float &f) {
         v[j * nx + i] = std::pair<clipper::Coord_orth, float>(co, f);
      }
      clipper::Coord_orth get_co(const unsigned int &i, const unsigned int &j) const {
         return v[j * nx + i].first;
      }
      float get_f(const unsigned int &i, const unsigned int &j) const {
         return v[j * nx + i].second;
      }
   };

   enum { MS_NO_CROSSING = -2,
          MS_NO_SQUARE = -1,
          MS_ALL_ABOVE,
          MS_UP_0_0,
          MS_UP_0_1,
          MS_UP_1_0,
          MS_UP_1_1,
          MS_UP_0_0_and_0_1,
          MS_UP_0_0_and_1_0,
          MS_UP_0_0_and_1_1, // hideous valley
          MS_UP_0_1_and_1_0, // hideous valley
          MS_UP_0_1_and_1_1,
          MS_UP_1_0_and_1_1,
          MS_UP_0_0_and_0_1_and_1_0,
          MS_UP_0_0_and_0_1_and_1_1,
          MS_UP_0_0_and_1_0_and_1_1,
          MS_UP_0_1_and_1_0_and_1_1,
   };

   int get_square_type(const unsigned int &i, const unsigned int &j,
                       const coord_array_2d &arr, const float &contour_level) {

      int square_type = MS_NO_SQUARE;
      const float &v00 = arr.get_f(i,   j);
      const float &v10 = arr.get_f(i+1, j);
      const float &v01 = arr.get_f(i,   j+1);
      const float &v11 = arr.get_f(i+1, j+1);
      if (v00 < contour_level)
         if (v01 < contour_level)
            if (v10 < contour_level)
               if (v11 < contour_level)
                  return MS_NO_CROSSING;
      if (v00 > contour_level)
         if (v01 > contour_level)
            if (v10 > contour_level)
               if (v11 > contour_level)
                  return MS_ALL_ABOVE;

      if (v00 < contour_level) {
         if (v01 < contour_level) {
            if (v10 < contour_level) {
               if (v11 < contour_level) {
                  return MS_NO_CROSSING;
               } else {
                  return MS_UP_1_1;
               }
            } else {
               if (v11 < contour_level) {
                  return MS_UP_1_0;
               } else {
                  return MS_UP_1_0_and_1_1;
               }
            }
         } else {
            if (v10 < contour_level) {
               if (v11 < contour_level) {
                  return MS_UP_0_1;
               } else {
                  return MS_UP_0_1_and_1_1;
               }
            } else {
               if (v11 < contour_level) {
                  return MS_UP_0_1_and_1_0;
               } else {
                  return MS_UP_0_1_and_1_0_and_1_1;
               }
            }
         }
      } else {
         if (v01 < contour_level) {
            if (v10 < contour_level) {
               if (v11 < contour_level) {
                  return MS_UP_0_0;
               } else {
                  return MS_UP_0_0_and_1_1;
               }
            } else {
               if (v11 < contour_level) {
                  return MS_UP_0_0_and_1_0;
               } else {
                  return MS_UP_0_0_and_1_0_and_1_1;
               }
            }
         } else {
            if (v10 < contour_level) {
               if (v11 < contour_level) {
                  return MS_UP_0_0_and_0_1;
               } else {
                  return MS_UP_0_0_and_0_1_and_1_1;
               }
            } else {
               if (v11 < contour_level) {
                  return MS_UP_0_0_and_0_1_and_1_0;
               }
            }
         }
      }
      return square_type;
   }

   // Make the 2D marching-squares cap mesh on the plane
   std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle>>
   make_2d_cap(const clipper::Xmap<float> &xmap,
               float contour_level,
               const clipper::Coord_orth &base_point,
               const clipper::Coord_orth &x_axis_uv,
               const clipper::Coord_orth &y_axis_uv,
               double x_axis_step_size,
               double y_axis_step_size,
               unsigned int n_x_axis_points,
               unsigned int n_y_axis_points,
               const clipper::Xmap<float> *other_map_for_colouring_p,
               float other_map_for_colouring_min_value,
               float other_map_for_colouring_max_value,
               float radial_map_colour_saturation) {

      std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle>> p;
      if (n_x_axis_points < 2) return p;
      if (n_y_axis_points < 2) return p;

      clipper::Coord_orth one_step_along_x = x_axis_step_size * x_axis_uv;
      clipper::Coord_orth one_step_along_y = y_axis_step_size * y_axis_uv;

      clipper::Coord_orth n(clipper::Coord_orth::cross(x_axis_uv, y_axis_uv));
      glm::vec3 normal = glm::normalize(glm::vec3(n.x(), n.y(), n.z()));

      // Fill the 2D grid with density values
      coord_array_2d arr(n_x_axis_points, n_y_axis_points);
      for (unsigned int i = 0; i < n_x_axis_points; i++) {
         for (unsigned int j = 0; j < n_y_axis_points; j++) {
            clipper::Coord_orth co = base_point +
               static_cast<double>(i) * one_step_along_x +
               static_cast<double>(j) * one_step_along_y;
            const float &f_m = coot::util::density_at_point(xmap, co);
            arr.set(i, j, co, f_m);
         }
      }

      const unsigned int OFFS = 10000;
      std::map<unsigned int, int> grid_coords_index_map;
      auto get_vertex_index = [](const unsigned int &i, const unsigned int &j,
                                 const coord_array_2d &arr_l,
                                 const std::map<unsigned int, int> &gcim) {
         const unsigned int OFFS_L = 10000;
         unsigned int idx = j * OFFS_L + i;
         auto it = gcim.find(idx);
         if (it == gcim.end()) return static_cast<int>(-1);
         return static_cast<int>(it->second);
      };

      float minv = other_map_for_colouring_min_value;
      float maxv = other_map_for_colouring_max_value;
      float cl = contour_level;
      float sat = radial_map_colour_saturation;

      if (!other_map_for_colouring_p) {
         minv = 0.9f * cl;
         maxv = 6.0f * cl;
      }

      const clipper::Xmap<float> *xmap_p = &xmap;
      const clipper::Xmap<float> *om = other_map_for_colouring_p;
      std::map<unsigned int, float> map_f;

      auto map_col_func_1 = [sat, minv, maxv](const float &dv) {
         coot::colour_holder cc(0.6f + 0.4f * sat, 0.6f - 0.6f * sat, 0.6f - 0.6f * sat);
         float fraction = 0.0f;
         if (dv > minv) {
            if (dv > maxv) {
               fraction = 1.0f;
            } else {
               float range = maxv - minv;
               float m = dv - minv;
               fraction = m / range;
            }
            cc.rotate_by(0.66f * fraction);
         }
         return glm::vec4(cc.red, cc.green, cc.blue, 1.0f);
      };

      auto map_col_func_2 = [sat, minv, maxv, &arr, xmap_p, om, &map_col_func_1](
                                const unsigned int &i, const unsigned int &j,
                                std::map<unsigned int, float> &mf) {
         const unsigned int OFFS_L = 10000;
         unsigned int idx = j * OFFS_L + i;
         auto it = mf.find(idx);
         if (it != mf.end()) {
            return map_col_func_1(it->second);
         } else {
            const clipper::Coord_orth &co = arr.get_co(i, j);
            float f = -1;
            if (om)
               f = coot::util::density_at_point(*om, co);
            else
               f = coot::util::density_at_point(*xmap_p, co);
            mf[idx] = f;
            return map_col_func_1(f);
         }
      };

      auto map_col_func_3 = [xmap_p, om, &map_col_func_1](const clipper::Coord_orth &co) {
         float f = -1;
         if (om)
            f = coot::util::density_at_point(*om, co);
         else
            f = coot::util::density_at_point(*xmap_p, co);
         return map_col_func_1(f);
      };

      auto clipper_to_glm = [](const clipper::Coord_orth &co) {
         return glm::vec3(co.x(), co.y(), co.z());
      };

      for (unsigned int i = 0; i < (n_x_axis_points - 1); i++) {
         for (unsigned int j = 0; j < (n_y_axis_points - 1); j++) {
            int ms_type = get_square_type(i, j, arr, contour_level);
            if ((ms_type != MS_NO_CROSSING) && (ms_type != MS_NO_SQUARE)) {
               const float &v00 = arr.get_f(i, j);
               const float &v10 = arr.get_f(i+1, j);
               const float &v01 = arr.get_f(i, j+1);
               const float &v11 = arr.get_f(i+1, j+1);
               clipper::Coord_orth co00 = arr.get_co(i, j);
               clipper::Coord_orth co10 = arr.get_co(i+1, j);
               clipper::Coord_orth co01 = arr.get_co(i, j+1);
               clipper::Coord_orth co11 = arr.get_co(i+1, j+1);

               int idx_00 = get_vertex_index(i,   j,   arr, grid_coords_index_map);
               int idx_10 = get_vertex_index(i+1, j,   arr, grid_coords_index_map);
               int idx_01 = get_vertex_index(i,   j+1, arr, grid_coords_index_map);
               int idx_11 = get_vertex_index(i+1, j+1, arr, grid_coords_index_map);
               if (idx_00 == -1) {
                  p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(co00), normal, map_col_func_2(i, j, map_f)));
                  idx_00 = p.first.size() - 1;
                  grid_coords_index_map[j * OFFS + i] = idx_00;
               }
               if (idx_10 == -1) {
                  p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(co10), normal, map_col_func_2(i+1, j, map_f)));
                  idx_10 = p.first.size() - 1;
                  grid_coords_index_map[j * OFFS + i+1] = idx_10;
               }
               if (idx_01 == -1) {
                  p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(co01), normal, map_col_func_2(i, j+1, map_f)));
                  idx_01 = p.first.size() - 1;
                  grid_coords_index_map[(j+1) * OFFS + i] = idx_01;
               }
               if (idx_11 == -1) {
                  p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(co11), normal, map_col_func_2(i+1, j+1, map_f)));
                  idx_11 = p.first.size() - 1;
                  grid_coords_index_map[(j+1) * OFFS + i + 1] = idx_11;
               }

               switch (ms_type) {

               case MS_ALL_ABOVE:
                  {
                     p.second.push_back(g_triangle(idx_00, idx_10, idx_11));
                     p.second.push_back(g_triangle(idx_00, idx_11, idx_01));
                  }
                  break;

               case MS_UP_0_0:
                  {
                     float frac_x1 = (v00 - contour_level) / (v00 - v10);
                     float frac_y1 = (v00 - contour_level) / (v00 - v01);
                     clipper::Coord_orth pxc(co00 + clipper::Coord_orth(frac_x1 * one_step_along_x));
                     clipper::Coord_orth pyc(co00 + clipper::Coord_orth(frac_y1 * one_step_along_y));
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pyc), normal, map_col_func_3(pyc)));
                     p.second.push_back(g_triangle(idx_00, idx_base, idx_base + 1));
                  }
                  break;

               case MS_UP_0_1_and_1_0_and_1_1:
                  {
                     float frac_x1 = (v00 - contour_level) / (v00 - v10);
                     float frac_y1 = (v00 - contour_level) / (v00 - v01);
                     clipper::Coord_orth pxc(co00 + frac_x1 * one_step_along_x);
                     clipper::Coord_orth pyc(co00 + frac_y1 * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pyc), normal, map_col_func_3(pyc)));
                     p.second.push_back(g_triangle(idx_11, idx_01, idx_base + 1));
                     p.second.push_back(g_triangle(idx_11, idx_base + 1, idx_base));
                     p.second.push_back(g_triangle(idx_11, idx_base, idx_10));
                  }
                  break;

               case MS_UP_0_1:
                  {
                     float frac_x = (v01 - contour_level) / (v01 - v11);
                     float frac_y = (contour_level - v00) / (v01 - v00);
                     clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                     clipper::Coord_orth pxy(co00 + frac_y * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxy), normal, map_col_func_3(pxy)));
                     p.second.push_back(g_triangle(idx_01, idx_base, idx_base + 1));
                  }
                  break;

               case MS_UP_0_0_and_1_0_and_1_1:
                  {
                     float frac_x = (v01 - contour_level) / (v01 - v11);
                     float frac_y = (contour_level - v00) / (v01 - v00);
                     clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                     clipper::Coord_orth pxy(co00 + frac_y * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxy), normal, map_col_func_3(pxy)));
                     p.second.push_back(g_triangle(idx_00, idx_10, idx_base + 1));
                     p.second.push_back(g_triangle(idx_10, idx_base, idx_base + 1));
                     p.second.push_back(g_triangle(idx_10, idx_11, idx_base));
                  }
                  break;

               case MS_UP_1_0:
                  {
                     float frac_x = (contour_level - v00) / (v10 - v00);
                     float frac_y = (v10 - contour_level) / (v10 - v11);
                     clipper::Coord_orth pxc(co00 + frac_x * one_step_along_x);
                     clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxy), normal, map_col_func_3(pxy)));
                     p.second.push_back(g_triangle(idx_base, idx_10, idx_base + 1));
                  }
                  break;

               case MS_UP_0_0_and_0_1_and_1_1:
                  {
                     float frac_x = (contour_level - v00) / (v10 - v00);
                     float frac_y = (v10 - contour_level) / (v10 - v11);
                     clipper::Coord_orth pxc(co00 + frac_x * one_step_along_x);
                     clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxy), normal, map_col_func_3(pxy)));
                     p.second.push_back(g_triangle(idx_01, idx_00, idx_base));
                     p.second.push_back(g_triangle(idx_01, idx_base, idx_base + 1));
                     p.second.push_back(g_triangle(idx_01, idx_base + 1, idx_11));
                  }
                  break;

               case MS_UP_1_1:
                  {
                     float frac_x = (contour_level - v01) / (v11 - v01);
                     float frac_y = (contour_level - v10) / (v11 - v10);
                     clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                     clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxy), normal, map_col_func_3(pxy)));
                     p.second.push_back(g_triangle(idx_11, idx_base, idx_base + 1));
                  }
                  break;

               case MS_UP_0_0_and_0_1_and_1_0:
                  {
                     float frac_x = (contour_level - v01) / (v11 - v01);
                     float frac_y = (contour_level - v10) / (v11 - v10);
                     clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                     clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxc), normal, map_col_func_3(pxc)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(pxy), normal, map_col_func_3(pxy)));
                     p.second.push_back(g_triangle(idx_00, idx_10, idx_base + 1));
                     p.second.push_back(g_triangle(idx_00, idx_base + 1, idx_base));
                     p.second.push_back(g_triangle(idx_00, idx_base, idx_01));
                  }
                  break;

               // 2 up cases

               case MS_UP_0_0_and_0_1:
                  {
                     float frac_x1 = -(contour_level - v01) / (v01 - v11);
                     float frac_x2 = -(contour_level - v00) / (v00 - v10);
                     clipper::Coord_orth px1c(co01 + frac_x1 * one_step_along_x);
                     clipper::Coord_orth px2c(co00 + frac_x2 * one_step_along_x);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(px1c), normal, map_col_func_3(px1c)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(px2c), normal, map_col_func_3(px2c)));
                     p.second.push_back(g_triangle(idx_00, idx_base, idx_01));
                     p.second.push_back(g_triangle(idx_00, idx_base + 1, idx_base));
                  }
                  break;

               case MS_UP_1_0_and_1_1:
                  {
                     float frac_x1 = -(contour_level - v01) / (v01 - v11);
                     float frac_x2 = -(contour_level - v00) / (v00 - v10);
                     clipper::Coord_orth px1c(co01 + frac_x1 * one_step_along_x);
                     clipper::Coord_orth px2c(co00 + frac_x2 * one_step_along_x);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(px1c), normal, map_col_func_3(px1c)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(px2c), normal, map_col_func_3(px2c)));
                     p.second.push_back(g_triangle(idx_10, idx_base, idx_base + 1));
                     p.second.push_back(g_triangle(idx_10, idx_11, idx_base));
                  }
                  break;

               case MS_UP_0_0_and_1_0:
                  {
                     float frac_y1 = (v00 - contour_level) / (v00 - v01);
                     float frac_y2 = (v10 - contour_level) / (v10 - v11);
                     clipper::Coord_orth py1c(co00 + frac_y1 * one_step_along_y);
                     clipper::Coord_orth py2c(co10 + frac_y2 * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(py1c), normal, map_col_func_3(py1c)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(py2c), normal, map_col_func_3(py2c)));
                     p.second.push_back(g_triangle(idx_00, idx_10, idx_base + 1));
                     p.second.push_back(g_triangle(idx_00, idx_base + 1, idx_base));
                  }
                  break;

               case MS_UP_0_1_and_1_1:
                  {
                     float frac_y1 = (v00 - contour_level) / (v00 - v01);
                     float frac_y2 = (v10 - contour_level) / (v10 - v11);
                     clipper::Coord_orth py1c(co00 + frac_y1 * one_step_along_y);
                     clipper::Coord_orth py2c(co10 + frac_y2 * one_step_along_y);
                     unsigned int idx_base = p.first.size();
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(py1c), normal, map_col_func_3(py1c)));
                     p.first.push_back(coot::api::vnc_vertex(clipper_to_glm(py2c), normal, map_col_func_3(py2c)));
                     p.second.push_back(g_triangle(idx_11, idx_base + 1, idx_base));
                     p.second.push_back(g_triangle(idx_11, idx_01, idx_base));
                  }
                  break;

               case MS_UP_0_0_and_1_1:
                  break;

               case MS_UP_0_1_and_1_0:
                  break;
               }
            }
         }
      }
      return p;
   }

   // Trim the 3D isosurface mesh, keeping only triangles fully behind the plane,
   // find boundary edges, and stitch to the 2D cap mesh boundaries.
   void trim_and_stitch(const std::vector<coot::api::vnc_vertex> &iso_verts,
                        const std::vector<g_triangle> &iso_tris,
                        const glm::vec3 &plane_normal,
                        const glm::vec3 &plane_point,
                        const std::vector<coot::api::vnc_vertex> &cap_verts,
                        const std::vector<g_triangle> &cap_tris,
                        coot::simple_mesh_t &mesh) {

      // Trim: keep triangles fully behind the plane
      std::vector<g_triangle> trimmed_triangles;
      trimmed_triangles.reserve(iso_tris.size());
      for (const auto &tri : iso_tris) {
         const glm::vec3 &p0 = iso_verts[tri.point_id[0]].pos;
         const glm::vec3 &p1 = iso_verts[tri.point_id[1]].pos;
         const glm::vec3 &p2 = iso_verts[tri.point_id[2]].pos;
         float d0 = glm::dot(p0 - plane_point, plane_normal);
         float d1 = glm::dot(p1 - plane_point, plane_normal);
         float d2 = glm::dot(p2 - plane_point, plane_normal);
         if (d0 < 0.0f && d1 < 0.0f && d2 < 0.0f)
            trimmed_triangles.push_back(tri);
      }

      // Find boundary edges of the cap mesh
      std::map<std::pair<unsigned int, unsigned int>, int> cap_edge_count;
      for (const auto &tri : cap_tris) {
         for (int e = 0; e < 3; e++) {
            unsigned int a = tri.point_id[e];
            unsigned int b = tri.point_id[(e+1) % 3];
            auto edge = std::make_pair(std::min(a, b), std::max(a, b));
            cap_edge_count[edge]++;
         }
      }
      std::vector<std::pair<unsigned int, unsigned int>> cap_boundary_edges;
      for (const auto &ec : cap_edge_count)
         if (ec.second == 1)
            cap_boundary_edges.push_back(ec.first);

      // Find boundary edges of trimmed 3D mesh near the plane
      std::map<std::pair<unsigned int, unsigned int>, int> mesh_edge_count;
      for (const auto &tri : trimmed_triangles) {
         for (int e = 0; e < 3; e++) {
            unsigned int a = tri.point_id[e];
            unsigned int b = tri.point_id[(e+1) % 3];
            auto edge = std::make_pair(std::min(a, b), std::max(a, b));
            mesh_edge_count[edge]++;
         }
      }
      float plane_threshold = 2.0f;
      std::vector<std::pair<unsigned int, unsigned int>> mesh_boundary_edges_near_plane;
      for (const auto &ec : mesh_edge_count) {
         if (ec.second == 1) {
            const glm::vec3 &pa = iso_verts[ec.first.first].pos;
            const glm::vec3 &pb = iso_verts[ec.first.second].pos;
            float da = std::abs(glm::dot(pa - plane_point, plane_normal));
            float db = std::abs(glm::dot(pb - plane_point, plane_normal));
            if (da < plane_threshold && db < plane_threshold)
               mesh_boundary_edges_near_plane.push_back(ec.first);
         }
      }

      // The cap mesh vertices are already in the output mesh (indices 0..cap_vert_count-1).
      // The trimmed 3D mesh vertices will be appended, offset by cap_vert_count.
      unsigned int cap_vert_count = cap_verts.size();

      // Add trimmed 3D mesh vertices and triangles to the output mesh
      for (const auto &v : iso_verts)
         mesh.vertices.push_back(v);
      for (const auto &tri : trimmed_triangles) {
         g_triangle t(tri.point_id[0] + cap_vert_count,
                      tri.point_id[1] + cap_vert_count,
                      tri.point_id[2] + cap_vert_count);
         mesh.triangles.push_back(t);
      }

      // Project a point onto the cap plane
      auto project_onto_plane = [&plane_point, &plane_normal](const glm::vec3 &p) -> glm::vec3 {
         float d = glm::dot(p - plane_point, plane_normal);
         return p - d * plane_normal;
      };

      // Precompute projected midpoints
      std::vector<glm::vec3> cap_edge_midpoints(cap_boundary_edges.size());
      for (std::size_t i = 0; i < cap_boundary_edges.size(); i++) {
         const glm::vec3 &a = cap_verts[cap_boundary_edges[i].first].pos;
         const glm::vec3 &b = cap_verts[cap_boundary_edges[i].second].pos;
         cap_edge_midpoints[i] = project_onto_plane(0.5f * (a + b));
      }

      std::vector<glm::vec3> mesh_edge_midpoints(mesh_boundary_edges_near_plane.size());
      for (std::size_t i = 0; i < mesh_boundary_edges_near_plane.size(); i++) {
         const glm::vec3 &a = iso_verts[mesh_boundary_edges_near_plane[i].first].pos;
         const glm::vec3 &b = iso_verts[mesh_boundary_edges_near_plane[i].second].pos;
         mesh_edge_midpoints[i] = project_onto_plane(0.5f * (a + b));
      }

      // Stitch: match edges bidirectionally and create quads
      std::set<std::pair<int, int>> stitched_pairs;
      float max_stitch_dist_sq = 9.0f; // 3.0 Angstroms squared

      auto add_stitch_quad = [&](std::size_t mi, std::size_t ci) {
         auto key = std::make_pair(static_cast<int>(mi), static_cast<int>(ci));
         if (stitched_pairs.count(key)) return;
         stitched_pairs.insert(key);

         unsigned int m0 = mesh_boundary_edges_near_plane[mi].first  + cap_vert_count;
         unsigned int m1 = mesh_boundary_edges_near_plane[mi].second + cap_vert_count;
         unsigned int c0 = cap_boundary_edges[ci].first;
         unsigned int c1 = cap_boundary_edges[ci].second;

         // Choose vertex pairing to avoid twisted quad
         const glm::vec3 &ma = iso_verts[mesh_boundary_edges_near_plane[mi].first].pos;
         const glm::vec3 &mb = iso_verts[mesh_boundary_edges_near_plane[mi].second].pos;
         float dist_straight = glm::length(ma - cap_verts[c0].pos) + glm::length(mb - cap_verts[c1].pos);
         float dist_crossed  = glm::length(ma - cap_verts[c1].pos) + glm::length(mb - cap_verts[c0].pos);
         if (dist_crossed < dist_straight)
            std::swap(c0, c1);

         mesh.triangles.push_back(g_triangle(m0, m1, c0));
         mesh.triangles.push_back(g_triangle(m1, c1, c0));
      };

      // Forward: 3D mesh edges -> nearest cap edges
      for (std::size_t mi = 0; mi < mesh_boundary_edges_near_plane.size(); mi++) {
         float best_dist_sq = max_stitch_dist_sq;
         int best_ci = -1;
         for (std::size_t ci = 0; ci < cap_boundary_edges.size(); ci++) {
            glm::vec3 diff = cap_edge_midpoints[ci] - mesh_edge_midpoints[mi];
            float d2 = glm::dot(diff, diff);
            if (d2 < best_dist_sq) { best_dist_sq = d2; best_ci = ci; }
         }
         if (best_ci >= 0) add_stitch_quad(mi, best_ci);
      }

      // Reverse: cap edges -> nearest 3D mesh edges
      for (std::size_t ci = 0; ci < cap_boundary_edges.size(); ci++) {
         float best_dist_sq = max_stitch_dist_sq;
         int best_mi = -1;
         for (std::size_t mi = 0; mi < mesh_boundary_edges_near_plane.size(); mi++) {
            glm::vec3 diff = mesh_edge_midpoints[mi] - cap_edge_midpoints[ci];
            float d2 = glm::dot(diff, diff);
            if (d2 < best_dist_sq) { best_dist_sq = d2; best_mi = mi; }
         }
         if (best_mi >= 0) add_stitch_quad(best_mi, ci);
      }
   }

} // anonymous namespace


coot::simple_mesh_t
coot::make_map_cap_mesh(const clipper::Xmap<float> &xmap,
                        float contour_level,
                        const clipper::Coord_orth &base_point,
                        const clipper::Coord_orth &x_axis_uv,
                        const clipper::Coord_orth &y_axis_uv,
                        double x_axis_step_size,
                        double y_axis_step_size,
                        unsigned int n_x_axis_points,
                        unsigned int n_y_axis_points,
                        const std::vector<api::vnc_vertex> &isosurface_vertices,
                        const std::vector<g_triangle> &isosurface_triangles,
                        const clipper::Xmap<float> *other_map_for_colouring_p,
                        float other_map_for_colouring_min_value,
                        float other_map_for_colouring_max_value,
                        float radial_map_colour_saturation) {

   simple_mesh_t mesh;

   // 1. Make the 2D cap mesh (marching squares)
   auto cap = make_2d_cap(xmap, contour_level, base_point, x_axis_uv, y_axis_uv,
                          x_axis_step_size, y_axis_step_size,
                          n_x_axis_points, n_y_axis_points,
                          other_map_for_colouring_p,
                          other_map_for_colouring_min_value,
                          other_map_for_colouring_max_value,
                          radial_map_colour_saturation);

   // Start the output mesh with the cap vertices and triangles
   mesh.vertices = cap.first;
   mesh.triangles = cap.second;

   // 2. Trim the 3D isosurface, find boundaries, and stitch
   if (!isosurface_vertices.empty() && !isosurface_triangles.empty()) {
      clipper::Coord_orth z_axis_uv(clipper::Coord_orth::cross(x_axis_uv, y_axis_uv));
      glm::vec3 plane_normal(z_axis_uv.x(), z_axis_uv.y(), z_axis_uv.z());
      glm::vec3 plane_point(base_point.x(), base_point.y(), base_point.z());

      trim_and_stitch(isosurface_vertices, isosurface_triangles,
                      plane_normal, plane_point,
                      cap.first, cap.second, mesh);
   }

   return mesh;
}
