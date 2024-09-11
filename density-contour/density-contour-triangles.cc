/*
 * density-contour/density-contour-triangles.cc
 * 
 * Copyright 2007 by The University of York
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include <algorithm>

#include "CIsoSurface.h"

void
coot::density_contour_triangles_container_t::depth_sort(const clipper::Coord_orth &back_plane_point,
							const clipper::Coord_orth &front_plane_point) {

   clipper::Coord_orth back_front = front_plane_point - back_plane_point;
   double bf_squared = back_front.lengthsq();
   if (bf_squared < 0.000001)
      bf_squared = 0.000001;
   for (unsigned int i=0; i<point_indices.size(); i++) {
      clipper::Coord_orth h_back = point_indices[i].mid_point - back_plane_point;
      double dot = clipper::Coord_orth::dot(back_front, h_back);
      point_indices[i].back_front_projection_distance = dot * dot / bf_squared;
   }
   std::sort(point_indices.begin(), point_indices.end());

}

void
coot::density_contour_triangles_container_t::remove_small_triangles() {

   float min_length = 0.1; // just for now
   float min_length_sqrd = min_length * min_length;

   for (unsigned int i=0; i<point_indices.size(); i++) {
      TRIANGLE &tri = point_indices[i];
      double dd_01 = (points[tri.pointID[0]] - points[tri.pointID[1]]).lengthsq();
      double dd_12 = (points[tri.pointID[1]] - points[tri.pointID[2]]).lengthsq();
      double dd_02 = (points[tri.pointID[0]] - points[tri.pointID[2]]).lengthsq();
      if (false) { // 0.1 is quite small for 1.8A resolution data
         std::cout << "tri-side " << dd_01 << "\n";
         std::cout << "tri-side " << dd_12 << "\n";
         std::cout << "tri-side " << dd_02 << "\n";
      }
      if (dd_01 < min_length_sqrd) {
         // merge points
         tri.reject_this = true;
      }
      if (dd_12 < min_length_sqrd) {
         // merge points
         tri.reject_this = true;
      }
      if (dd_02 < min_length_sqrd) {
         // merge points
         tri.reject_this = true;
      }
   }

   auto local_points = points;

   auto remover = [min_length_sqrd, local_points] (const TRIANGLE &tri) {
                     double dd_01 = (local_points[tri.pointID[0]] - local_points[tri.pointID[1]]).lengthsq();
                     double dd_12 = (local_points[tri.pointID[1]] - local_points[tri.pointID[2]]).lengthsq();
                     double dd_02 = (local_points[tri.pointID[0]] - local_points[tri.pointID[2]]).lengthsq();
                     if (dd_01 < min_length_sqrd) return true;
                     if (dd_12 < min_length_sqrd) return true;
                     if (dd_02 < min_length_sqrd) return true;
                     return false;
                  };

   point_indices.erase(std::remove_if(point_indices.begin(), point_indices.end(), remover),
                       point_indices.end());
   

}

void
coot::density_contour_triangles_container_t::calculate_normals() {

   bool do_normal_averaging = true;

   // for each triangle
   //    for each of the 3 vertices
   //       add the flat shading normal to the sum of normals of that vertex
   // then normalise the normals

   // I think that some normals are not set because there is no triangle that uses them.
   // (sphere filtering, I guess).
   // So, just set that normal to (0,0,1) to remove the nans generated when calculating
   // the unit() of (0,0,0);

   // validate
   if (false) {
      for (unsigned int i=0; i<point_indices.size(); i++) {
         std::cout << "normal for point index " << i << " " << point_indices[i].normal_for_flat_shading.format()
                   << std::endl;
      }
   }

   if (do_normal_averaging) {

      std::vector<clipper::Coord_orth> sum_normals(normals.size());
      std::vector<unsigned int> n_contribs(normals.size(), 0);
      clipper::Coord_orth zero(0,0,0);
      for (unsigned int i=0; i<sum_normals.size(); i++)
         sum_normals[i] = zero;

      for (unsigned int i=0; i<point_indices.size(); i++) {
         for (int j=0; j<3; j++) {
            sum_normals[point_indices[i].pointID[j]] += point_indices[i].normal_for_flat_shading;
            n_contribs[point_indices[i].pointID[j]]++;
         }
      }
      for (unsigned int i=0; i<points.size(); i++) {
         if (n_contribs[i] > 0)
            normals[i] = clipper::Coord_orth(sum_normals[i].unit());
         else
            normals[i] = clipper::Coord_orth(0,0,1);
      }

      // more validation that I don't want to delete at the moment.
      if (false) {
         std::cout << "in calculate_normals() points.size()  " << points.size() << std::endl;
         std::cout << "in calculate_normals() normals.size() " << normals.size() << std::endl;
         for (unsigned int i=0; i<normals.size(); i++) {
            std::cout << "in calculate_normals() normal " << i << " " << normals[i].format() << std::endl;
         }
      }
   } else {

      // something here.
   }
}

#include <chrono>
#include "coot-utils/coot-map-utils.hh"

// use the function gradient
void
coot::density_contour_triangles_container_t::calculate_normals_for_vertices(const clipper::Xmap<float> &xmap) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   float delta = 0.03;
   for (unsigned int i=0; i<points.size(); i++) {
      const auto &pos = points[i];
      clipper::Coord_orth p_x_1(pos.x() - delta, pos.y(), pos.z());
      clipper::Coord_orth p_x_2(pos.x() + delta, pos.y(), pos.z());
      clipper::Coord_orth p_y_1(pos.x(), pos.y() - delta, pos.z());
      clipper::Coord_orth p_y_2(pos.x(), pos.y() + delta, pos.z());
      clipper::Coord_orth p_z_1(pos.x(), pos.y(), pos.z() - delta);
      clipper::Coord_orth p_z_2(pos.x(), pos.y(), pos.z() + delta);
      float f_x_1 = util::density_at_point(xmap, p_x_1);
      float f_x_2 = util::density_at_point(xmap, p_x_2);
      float f_y_1 = util::density_at_point(xmap, p_y_1);
      float f_y_2 = util::density_at_point(xmap, p_y_2);
      float f_z_1 = util::density_at_point(xmap, p_z_1);
      float f_z_2 = util::density_at_point(xmap, p_z_2);
      clipper::Coord_orth grr(f_x_1 - f_x_2, f_y_1 - f_y_2, f_z_1 - f_z_2);
      clipper::Coord_orth gr(grr.unit());
      if (false) {
         float f = util::density_at_point(xmap, pos);
         const auto &n = normals[i];
         std::cout << "pos " << pos.x() << " " << pos.y() << " " << pos.z() << " "
                   << "grr " << grr.x() << " " << grr.y() << " " << grr.z() << " "
                   << "normal " << n.x() << " " << n.y() << " " << n.z() << " "
                   << "gr  " <<  gr.x() << " " <<  gr.y() << " " <<  gr.z() << " "
                   << " f " << f << " "
                   << "fx  " << f_x_1 << " " << f_x_2 << " "
                   << "fy  " << f_y_1 << " " << f_y_2 << " "
                   << "fz  " << f_z_1 << " " << f_z_2 << " "
                   << std::endl;
      }
      normals[i] = gr;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();

   if (false)
      std::cout << "normals_from_function_gradient(): time " << d10 << " ms " << std::endl;
}
