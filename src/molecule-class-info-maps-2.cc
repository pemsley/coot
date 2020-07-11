/* src/molecule-class-info-maps.cc
 * 
 * Copyright 2015 by Medical Research Council
 * 
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#define GLM_ENABLE_EXPERIMENTAL // # for to_string
#include <glm/ext.hpp>

#include "molecule-class-info.h"

// shift "bottom left" to the origin and make sure that it's on a grid that is
// "nice" (acceptable?) for molrep.
//
// Now the resulting map has gompertz scaling of density (drops off
// towards the edge) and border added (8% of radius).
// 
int
molecule_class_info_t::export_map_fragment_with_origin_shift(float radius,
							     clipper::Coord_orth centre,
							     const std::string &file_name) const {

   // centre is the centre of the map of this molecule (that we are extracting from)

   int r = 0;
   if (has_xmap()) {
      clipper::Cell          xmap_cell = xmap.cell();
      clipper::Grid_sampling xmap_grid_sampling = xmap.grid_sampling();
      clipper::Coord_orth centre_moved = centre - clipper::Coord_orth(0.1, 0.1, 0.1);

      clipper::Grid_range gr0(xmap_cell, xmap_grid_sampling, radius);
      clipper::Grid_range gr1(gr0.min() + centre.coord_frac(xmap_cell).coord_grid(xmap_grid_sampling),
			      gr0.max() + centre.coord_frac(xmap_cell).coord_grid(xmap_grid_sampling));
      
      int nx_u = 2*gr0.max().u();
      int nx_v = 2*gr0.max().v();
      int nx_w = 2*gr0.max().w();
      
      clipper::Coord_grid nxmap_grid(nx_u, nx_v, nx_w);
      clipper::Coord_grid nxmap_grid_max(nx_u-1, nx_v-1, nx_w-1);
      clipper::Coord_grid nxmap_origin(0,0,0);
      clipper::Grid_range nxmap_grid_range(nxmap_origin, nxmap_grid_max);
      clipper::Grid_sampling nxmap_grid_sampling(nx_u, nx_v, nx_w);
      std::cout << "--------------- nxmap grid_sampling init with " << nx_u << " " << nx_v << " " << nx_w
		<< std::endl;
      
      std::cout << "--------------- gr0.min() " << gr0.min().format() << std::endl;
      std::cout << "--------------- gr0.max() " << gr0.max().format() << std::endl;
      std::cout << "--------------- nxmap_grid_sampling "   << nxmap_grid_sampling.format() << std::endl;
      std::cout << "--------------- nxmap grid range min: " << nxmap_grid_range.min().format() << std::endl;
      std::cout << "--------------- nxmap grid range max: " << nxmap_grid_range.max().format() << std::endl;
      

      double nx_alpha = xmap.cell().descr().alpha();
      double nx_beta  = xmap.cell().descr().beta();
      double nx_gamma = xmap.cell().descr().gamma();
      double nx_a = xmap.cell().descr().a() * double(nxmap_grid_sampling.nu())/double(xmap_grid_sampling.nu());
      double nx_b = xmap.cell().descr().b() * double(nxmap_grid_sampling.nv())/double(xmap_grid_sampling.nv());
      double nx_c = xmap.cell().descr().c() * double(nxmap_grid_sampling.nw())/double(xmap_grid_sampling.nw());
      clipper::Cell_descr nxmap_cell_descr(nx_a, nx_b, nx_c, nx_alpha, nx_beta, nx_gamma);
      clipper::Cell nxmap_cell(nxmap_cell_descr);

      // init nxmap
      clipper::NXmap<float> nxmap_local(nxmap_cell, nxmap_grid_sampling, nxmap_grid_range);
      
      clipper::Xmap<float>::Map_reference_coord ix(xmap);
      clipper::Coord_orth centre_radius(centre.x() - radius,
					centre.y() - radius,
					centre.z() - radius);
      clipper::Coord_orth nxmap_centre(radius, radius, radius);
      clipper::Coord_grid offset = xmap.coord_map(centre_radius).coord_grid();
      
      std::cout << "--------------- xmap offset to centre        " << centre_radius.format() << std::endl;
      std::cout << "--------------- xmap offset to centre (grid) " << offset.format() << std::endl;

      typedef clipper::NXmap<float>::Map_reference_index NRI;
      double limited_radius = radius * 0.92;
      for (NRI inx = nxmap_local.first(); !inx.last(); inx.next()) {
	 clipper::Coord_orth p = inx.coord().coord_frac(nxmap_grid_sampling).coord_orth(nxmap_cell);
	 double d_to_c_sq = clipper::Coord_orth(p-nxmap_centre).lengthsq();
	 if (d_to_c_sq > limited_radius*limited_radius) {
	    nxmap_local[inx] = 0.0;
	    if (0) 
	       std::cout << " inx " << inx.coord().format() << " " << d_to_c_sq << "  " << p.format() << " "
			 << centre.format() << " vs " << limited_radius*limited_radius << " is outside "
			 << std::endl;
	 } else {
	    // ix indexes the xmap
	    clipper::Coord_grid gp = p.coord_frac(xmap_cell).coord_grid(xmap_grid_sampling);
	    ix.set_coord(gp + offset);

	    // make a function that is y=1 around x=0 and y=1 around x=1 and
	    // falls of to 0.5 around x=0.8 or so.
	    //
	    double x = sqrt(d_to_c_sq)/limited_radius;
	    double gompertz_a = 0.14;
	    double gompertz_b = 0.1; 
	    double gompertz_c = 3;
	    double gompertz_scale = 1 - (-gompertz_a*1.1 + gompertz_a * exp (gompertz_b * exp(gompertz_c * x)));
	    nxmap_local[inx] = xmap[ix] * gompertz_scale;
	    if (0)
	       std::cout << " inx " << inx.coord().format() << " " << d_to_c_sq << "  " << p.format() << " " << centre.format()
			 << " vs " << limited_radius*limited_radius << " " << gompertz_scale << std::endl;
	 } 
	    
      }
      clipper::CCP4MAPfile mapout;
      mapout.open_write(file_name);
      mapout.set_cell(nxmap_cell);
      mapout.export_nxmap(nxmap_local);
      mapout.close_write();
      std::cout << "Exported map " << file_name << std::endl;
   }

   return r;
}


mean_and_variance<float>
molecule_class_info_t::set_and_get_histogram_values(unsigned int n_bins) {

   // fill map_histogram_values and return it.

   if (map_histogram_values.size() == 1) { // failed before
      //
   } else {
      if (map_histogram_values.size() > 0) {
          // use the cache
      } else {
         // uncached
         bool ignore_pseudo_zeros = false;
         mean_and_variance<float> mv =
            map_density_distribution(xmap, n_bins, false, ignore_pseudo_zeros);
         if (mv.size() == 0) {
            // add fake result
            mv.bins.push_back(0);
         }
         map_histogram_values = mv;
      }
   }
   return map_histogram_values;
}


// radial colouring
void
molecule_class_info_t::set_radial_map_colouring_centre(float x, float y, float z) {
   clipper::Coord_orth c(x,y,z);
   radial_map_colour_centre = c;
}

void
molecule_class_info_t::set_radial_map_colouring_min_radius(float r) {
   radial_map_colour_radius_min = r;
}

void
molecule_class_info_t::set_radial_map_colouring_max_radius(float r) {
   radial_map_colour_radius_max = r;

}

void
molecule_class_info_t::set_radial_map_colouring_invert(bool invert_state) {
   radial_map_colour_invert_flag = invert_state;
}

void
molecule_class_info_t::set_radial_map_colouring_saturation(float saturation) {
   radial_map_colour_saturation = saturation;
}

// we don't want to have triangles for MS_NO_CROSSING, but we want 2 triangles (a quad that fills
// the grid_ for MS_ALL_ABOVE
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

#include "coot-utils/coot-map-utils.hh"

// using current contour level,
// return world coordinates and normals
std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecule_class_info_t::make_map_cap(const clipper::Coord_orth &base_point,
                                    const clipper::Coord_orth &x_axis_uv,
                                    const clipper::Coord_orth &y_axis_uv,
                                    double x_axis_step_size,
                                    double y_axis_step_size,
                                    unsigned int n_x_axis_points,
                                    unsigned int n_y_axis_points) const {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > p;
   if (n_x_axis_points < 2) return p;
   if (n_y_axis_points < 2) return p;

   clipper::Coord_orth one_step_along_x = x_axis_step_size * x_axis_uv;
   clipper::Coord_orth one_step_along_y = y_axis_step_size * y_axis_uv;

   std::cout << "debug:: one_step_along_x " << one_step_along_x.format() << std::endl;
   std::cout << "debug:: one_step_along_y " << one_step_along_y.format() << std::endl;

   clipper::Coord_orth n(clipper::Coord_orth::cross(x_axis_uv, y_axis_uv));
   glm::vec3 normal = glm::normalize(glm::vec3(n.x(), n.y(), n.z()));
   glm::vec4 col(0.7, 0.7, 0.7, 1.0);

   //    [01]             [11]
   //      o ----------- o 
   //      |             |
   //      |             |
   //      |             |
   //      |             |
   //      |             |
   //      o ----------- o 
   //     [00]            [10]
   //

   // first fill the 2D grid with density values
   std::cout << "........... fill the grid other_map_for_colouring_p " << other_map_for_colouring_p
             << std::endl;

   float f = -999.9;
   coord_array_2d arr(n_x_axis_points, n_y_axis_points);
   for (unsigned int i=0; i<n_x_axis_points; i++) {
      for (unsigned int j=0; j<n_y_axis_points; j++) {
         clipper::Coord_orth co = base_point +
            static_cast<double>(i) * one_step_along_x +
            static_cast<double>(j) * one_step_along_y;
         const float &f_m = density_at_point(co);
         arr.set(i,j,co,f_m);
      }
   }


   // grid coordinates index map:
   // We want to store the index in p.second of coordinates of the grid points
   // because they will be use for multiple (4) triangles
   // No need to store (reuse) the "fraction" points
   const unsigned int OFFS = 10000;
   std::map<unsigned int, int> grid_coords_index_map;
   auto get_vertex_index = [] (const unsigned int &i,
                               const unsigned int &j,
                               const coord_array_2d &arr,
                               const std::map<unsigned int, int> &grid_coords_index_map) {
                              const unsigned int OFFS = 10000;
                              unsigned int idx = j * OFFS + i;
                              std::map<unsigned int, int>::const_iterator it;
                              it = grid_coords_index_map.find(idx);
                              if (it == grid_coords_index_map.end()) {
                                 return static_cast<int>(-1);
                              } else {
                                 return static_cast<int>(it->second);
                              }
                           };

   float minv = other_map_for_colouring_min_value;
   float maxv = other_map_for_colouring_max_value;
   float cl = contour_level;

   if (! other_map_for_colouring_p) {
      float s = map_sigma();
      minv = 0.9 * cl;
      maxv = 6.0 * s;
   }

   float sat = radial_map_colour_saturation;

   const clipper::Xmap<float> *xmap_this = &xmap;
   std::map<unsigned int, float> map_f;
   const clipper::Xmap<float> *om = 0;
   if (colour_map_using_other_map_flag)
      om = other_map_for_colouring_p;

   auto map_col_func_1 = [sat, minv, maxv] (const float &dv) {
                            coot::colour_t cc(0.6+0.4*sat, 0.6-0.6*sat, 0.6-0.6*sat);
                            float fraction = 0.0;
                            if (dv > minv) {
                               if (dv > maxv) {
                                  fraction = 1.0;
                               } else {
                                  float range = maxv - minv;
                                  float m = dv - minv;
                                  fraction = m/range;
                               }
                               cc.rotate(0.66 * fraction);
                            }
                            return glm::vec4(cc.col[0], cc.col[1], cc.col[2], 1.0f);
                         };

   auto map_col_func_2 = [sat, minv, maxv, arr, xmap_this, om, map_col_func_1] (const unsigned int &i, const unsigned int &j,
                                                                     std::map<unsigned int, float> &map_f) {
                            const unsigned int OFFS = 10000;
                            unsigned int idx = j * OFFS + i;
                            std::map<unsigned int, float>::const_iterator it;
                            it = map_f.find(idx);
                            if (it != map_f.end()) {
                               return map_col_func_1(it->second);
                            } else {
                               const clipper::Coord_orth &co = arr.get_co(i,j);
                               float f = -1;
                               if (om) {
                                  f = coot::util::density_at_point(*om, co);
                               } else {
                                  f = coot::util::density_at_point(*xmap_this, co);
                               }
                               map_f[idx] = f;
                               return map_col_func_1(f);
                            }
                         };

   auto map_col_func_3 = [xmap_this, om, map_col_func_1] (const clipper::Coord_orth &co) {
                            float f = -1;
                            if (om) {
                               f = coot::util::density_at_point(*om, co);
                            } else {
                               f = coot::util::density_at_point(*xmap_this, co);
                            }
                            return map_col_func_1(f);
                         };

   auto map_colour_function = [map_col_func_1] (const float &dv) { return map_col_func_1(dv); };


   auto clipper_to_glm = [] (const clipper::Coord_orth &co) { return glm::vec4(co.x(), co.y(), co.z(), 1.0); };

   for (unsigned int i=0; i<(n_x_axis_points-1); i++) {
      for (unsigned int j=0; j<(n_y_axis_points-1); j++) {
         int ms_type = get_square_type(i, j, arr, contour_level);
	 if ((ms_type != MS_NO_CROSSING) && (ms_type != MS_NO_SQUARE)) {
            const float &v00 = arr.get_f(i, j);
            const float &v10 = arr.get_f(i+1, j);
            const float &v01 = arr.get_f(i, j+1);
            const float &v11 = arr.get_f(i+1, j+1);
            clipper::Coord_orth co00 = arr.get_co(i,j);
            clipper::Coord_orth co10 = arr.get_co(i+1,j);
            clipper::Coord_orth co01 = arr.get_co(i,j+1);
            clipper::Coord_orth co11 = arr.get_co(i+1,j+1);

            int idx_00 = get_vertex_index(i,   j,   arr, grid_coords_index_map);
            int idx_10 = get_vertex_index(i+1, j,   arr, grid_coords_index_map);
            int idx_01 = get_vertex_index(i,   j+1, arr, grid_coords_index_map);
            int idx_11 = get_vertex_index(i+1, j+1, arr, grid_coords_index_map);
            if (idx_00 == -1) {
               p.first.push_back(s_generic_vertex(clipper_to_glm(co00), normal, map_col_func_2(i, j, map_f)));
               idx_00 = p.first.size() - 1;
               grid_coords_index_map[j * OFFS + i] = idx_00;
            }
            if (idx_10 == -1) {
               p.first.push_back(s_generic_vertex(clipper_to_glm(co10), normal, map_col_func_2(i+1, j, map_f)));
               idx_10 = p.first.size() - 1;
               grid_coords_index_map[j * OFFS + i+1] = idx_10;
            }
            if (idx_01 == -1) {
               p.first.push_back(s_generic_vertex(clipper_to_glm(co01), normal, map_col_func_2(i, j+1, map_f)));
               idx_01 = p.first.size() - 1;
               grid_coords_index_map[(j+1) * OFFS + i] = idx_01;
            }
            if (idx_11 == -1) {
               p.first.push_back(s_generic_vertex(clipper_to_glm(co11), normal, map_col_func_2(i, j+1, map_f)));
               idx_11 = p.first.size() - 1;
               grid_coords_index_map[(j+1) * OFFS + i + 1] = idx_11;
            }

            switch (ms_type) {

            case MS_ALL_ABOVE:

               {
                  g_triangle gt_0(idx_00, idx_10, idx_11);
                  g_triangle gt_1(idx_00, idx_11, idx_01);
                  p.second.push_back(gt_0);
                  p.second.push_back(gt_1);
               }

               break;

            case MS_UP_0_0:
               {
                  float frac_x1 = (v00-contour_level)/(v00-v10);
                  float frac_y1 = (v00-contour_level)/(v00-v01);  // tick
                  clipper::Coord_orth co00 = arr.get_co(i,j);
                  clipper::Coord_orth co10 = arr.get_co(i+1,j);
                  clipper::Coord_orth co01 = arr.get_co(i,j+1);
                  clipper::Coord_orth pxc(co00 + clipper::Coord_orth(frac_x1 * one_step_along_x));
                  clipper::Coord_orth pyc(co00 + clipper::Coord_orth(frac_y1 * one_step_along_y));
                  glm::vec3 px_glm(pxc.x(), pxc.y(), pxc.z());
                  glm::vec3 py_glm(pyc.x(), pyc.y(), pyc.z());
                  unsigned int idx_base = p.first.size();
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pyc);
                  s_generic_vertex vx(px_glm, normal, col_1);
                  s_generic_vertex vy(py_glm, normal, col_2);
                  p.first.push_back(vx);
                  p.first.push_back(vy);
                  g_triangle tri(idx_00, idx_base, idx_base + 1);
                  p.second.push_back(tri);
               }
               break;

            case MS_UP_0_1_and_1_0_and_1_1:  // reverse of above
               {
                  float frac_x1 = (v00-contour_level)/(v00-v10);
                  float frac_y1 = (v00-contour_level)/(v00-v01); // tick
                  clipper::Coord_orth pxc(co00 + frac_x1 * one_step_along_x);
                  clipper::Coord_orth pyc(co00 + frac_y1 * one_step_along_y);
                  glm::vec3 px_glm = clipper_to_glm(pxc);
                  glm::vec3 py_glm = clipper_to_glm(pyc);
                  // 3 triangles, 2 extra vertices
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pyc);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned int idx_base = p.first.size();
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  g_triangle gt0(idx_11, idx_01, idx_base + 1);
                  g_triangle gt1(idx_11, idx_base + 1, idx_base);
                  g_triangle gt2(idx_11, idx_base, idx_10);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
                  p.second.push_back(gt2);
               }
               break;

            case MS_UP_0_1:
               {
                  float frac_x = (v01-contour_level)/(v01-v11); 
                  float frac_y = (contour_level-v00)/(v01-v00);  // tick
                  if (frac_x < 0.0) std::cout << "ERROR MS_UP_0_1 frac_x " << frac_x << std::endl;
                  if (frac_y < 0.0) std::cout << "ERROR MS_UP_0_1 frac_y " << frac_y << std::endl;
                  if (frac_x > 1.0) std::cout << "ERROR MS_UP_0_1 frac_x " << frac_x << std::endl;
                  if (frac_y > 1.0) std::cout << "ERROR MS_UP_0_1 frac_y " << frac_y << std::endl;
                  clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                  clipper::Coord_orth pxy(co00 + frac_y * one_step_along_y);
                  glm::vec4 px_glm = clipper_to_glm(pxc);
                  glm::vec4 py_glm = clipper_to_glm(pxy);
                  // one triangle, 2 vertices
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pxy);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  g_triangle gt(idx_01, idx_base, idx_base + 1);
                  p.second.push_back(gt);

               }
               break;

            case MS_UP_0_0_and_1_0_and_1_1:
               {
                  float frac_x = (v01-contour_level)/(v01-v11);
                  float frac_y = (contour_level-v00)/(v01-v00);  // tick
                  if (frac_x < 0.0) std::cout << "ERROR MS_UP_0_0_and_1_0_and_1_1 frac_x " << frac_x << std::endl;
                  if (frac_y < 0.0) std::cout << "ERROR MS_UP_0_0_and_1_0_and_1_1 frac_y " << frac_y << std::endl;
                  if (frac_x > 1.0) std::cout << "ERROR MS_UP_0_0_and_1_0_and_1_1 frac_x " << frac_x << std::endl;
                  if (frac_y > 1.0) std::cout << "ERROR MS_UP_0_0_and_1_0_and_1_1 frac_y " << frac_y << std::endl;
                  clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                  clipper::Coord_orth pxy(co00 + frac_y * one_step_along_y);
                  glm::vec4 px_glm = clipper_to_glm(pxc);
                  glm::vec4 py_glm = clipper_to_glm(pxy);
                  // 3 triangle, same 2 vertices (as above)
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pxy);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  g_triangle gt0(idx_00, idx_10, idx_base + 1);
                  g_triangle gt1(idx_10, idx_base, idx_base + 1);
                  g_triangle gt2(idx_10, idx_11, idx_base);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
                  p.second.push_back(gt2);
               }
               break;

            case MS_UP_1_0:
               {
#if 1
                  float frac_x = (contour_level-v00)/(v10-v00);
                  float frac_y = (v10-contour_level)/(v10-v11); // tick
                  if (frac_x < 0.0) std::cout << "ERROR MS_UP_1_0 frac_x " << frac_x << std::endl;
                  if (frac_y < 0.0) std::cout << "ERROR MS_UP_1_0 frac_y " << frac_y << std::endl;
                  if (frac_x > 1.0) std::cout << "ERROR MS_UP_1_0 frac_x " << frac_x << std::endl;
                  if (frac_y > 1.0) std::cout << "ERROR MS_UP_1_0 frac_y " << frac_y << std::endl;
                  clipper::Coord_orth pxc(co00 + frac_x * one_step_along_x);
                  clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                  glm::vec4 px_glm = clipper_to_glm(pxc);
                  glm::vec4 py_glm = clipper_to_glm(pxy);
                  // 1 triangle,  2 new vertices
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pxy);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  g_triangle gt0(idx_base, idx_10, idx_base + 1);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
#endif
               }
               break;

            case MS_UP_0_0_and_0_1_and_1_1:  // reverse of above
               {
#if 1
                  float frac_x = (contour_level-v00)/(v10-v00);
                  float frac_y = (v10-contour_level)/(v10-v11); // tick
                  if (frac_x < 0.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_1 frac_x " << frac_x << std::endl;
                  if (frac_y < 0.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_1 frac_y " << frac_y << std::endl;
                  if (frac_x > 1.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_1 frac_x " << frac_x << std::endl;
                  if (frac_y > 1.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_1 frac_y " << frac_y << std::endl;
                  clipper::Coord_orth pxc(co00 + frac_x * one_step_along_x);
                  clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                  glm::vec4 px_glm = clipper_to_glm(pxc);
                  glm::vec4 py_glm = clipper_to_glm(pxy);
                  // 3 triangles, 2 vertices as above
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pxy);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  g_triangle gt0(idx_01, idx_00, idx_base);
                  g_triangle gt1(idx_01, idx_base, idx_base + 1);
                  g_triangle gt2(idx_01, idx_base + 1, idx_11);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
                  p.second.push_back(gt2);
#endif
               }
               break;
      
            case MS_UP_1_1:
               {
#if 1
                  float frac_x = (contour_level-v01)/(v11-v01);
                  float frac_y = (contour_level-v10)/(v11-v10);
                  if (frac_x < 0.0) std::cout << "ERROR MS_UP_1_1 frac_x " << frac_x << std::endl;
                  if (frac_y < 0.0) std::cout << "ERROR MS_UP_1_1 frac_y " << frac_y << std::endl;
                  if (frac_x > 1.0) std::cout << "ERROR MS_UP_1_1 frac_x " << frac_x << std::endl;
                  if (frac_y > 1.0) std::cout << "ERROR MS_UP_1_1 frac_y " << frac_y << std::endl;
                  clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                  clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                  glm::vec4 px_glm = clipper_to_glm(pxc);
                  glm::vec4 py_glm = clipper_to_glm(pxy);
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pxy);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned int idx_base = p.first.size();
                  g_triangle gt0(idx_11, idx_base, idx_base + 1);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
#endif
               }
               break;

            case MS_UP_0_0_and_0_1_and_1_0: // reverse of above
               {
#if 1
                  float frac_x = (contour_level-v01)/(v11-v01);
                  float frac_y = (contour_level-v10)/(v11-v10);
                  if (frac_x < 0.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_0 frac_x " << frac_x << std::endl;
                  if (frac_y < 0.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_0 frac_y " << frac_y << std::endl;
                  if (frac_x > 1.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_0 frac_x " << frac_x << std::endl;
                  if (frac_y > 1.0) std::cout << "ERROR MS_UP_0_0_and_0_1_and_1_0 frac_y " << frac_y << std::endl;
                  clipper::Coord_orth pxc(co01 + frac_x * one_step_along_x);
                  clipper::Coord_orth pxy(co10 + frac_y * one_step_along_y);
                  glm::vec4 px_glm = clipper_to_glm(pxc);
                  glm::vec4 py_glm = clipper_to_glm(pxy);
                  // 3 triangles, 2 vertices as above
                  glm::vec4 col_1 = map_col_func_3(pxc);
                  glm::vec4 col_2 = map_col_func_3(pxy);
                  s_generic_vertex v0(px_glm, normal, col_1);
                  s_generic_vertex v1(py_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  g_triangle gt0(idx_00, idx_10, idx_base + 1);
                  g_triangle gt1(idx_00, idx_base + 1, idx_base);
                  g_triangle gt2(idx_00, idx_base, idx_01);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
                  p.second.push_back(gt2);
#endif
               }
               break;

            // ------------------------------- 2 up --------------------

            case MS_UP_0_0_and_0_1:
               {
#if 1
                  float frac_x1 = -(contour_level-v01)/(v01-v11);
                  float frac_x2 = -(contour_level-v00)/(v00-v10);
                  if (frac_x1 < 0.0) std::cout << "ERROR MS_UP_0_0_and_0_1 frac_x1 " << frac_x1 << std::endl;
                  if (frac_x1 > 1.0) std::cout << "ERROR MS_UP_0_0_and_0_1 frac_x1 " << frac_x1 << std::endl;
                  if (frac_x2 < 0.0) std::cout << "ERROR MS_UP_0_0_and_0_1 frac_x2 " << frac_x2 << std::endl;
                  if (frac_x2 > 1.0) std::cout << "ERROR MS_UP_0_0_and_0_1 frac_x2 " << frac_x2 << std::endl;
                  clipper::Coord_orth px1c(co01 + frac_x1 * one_step_along_x);
                  clipper::Coord_orth px2c(co00 + frac_x2 * one_step_along_x);
                  glm::vec4 px1_glm = clipper_to_glm(px1c);
                  glm::vec4 px2_glm = clipper_to_glm(px2c);
                  // 2 triangles, 2 vertices
                  glm::vec4 col_1 = map_col_func_3(px1c);
                  glm::vec4 col_2 = map_col_func_3(px2c);
                  s_generic_vertex v0(px1_glm, normal, col_1);
                  s_generic_vertex v1(px2_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  g_triangle gt0(idx_00, idx_base, idx_01);
                  g_triangle gt1(idx_00, idx_base + 1, idx_base);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
#endif
               }
               break;

            case MS_UP_1_0_and_1_1:  // reverse of above
               {
#if 1
                  float frac_x1 = -(contour_level-v01)/(v01-v11);
                  float frac_x2 = -(contour_level-v00)/(v00-v10);
                  if (frac_x1 < 0.0) std::cout << "ERROR MS_UP_1_0_and_1_1 frac_x1 " << frac_x1 << std::endl;
                  if (frac_x1 > 1.0) std::cout << "ERROR MS_UP_1_0_and_1_1 frac_x1 " << frac_x1 << std::endl;
                  if (frac_x2 < 0.0) std::cout << "ERROR MS_UP_1_0_and_1_1 frac_x2 " << frac_x2 << std::endl;
                  if (frac_x2 > 1.0) std::cout << "ERROR MS_UP_1_0_and_1_1 frac_x2 " << frac_x2 << std::endl;
                  clipper::Coord_orth px1c(co01 + frac_x1 * one_step_along_x);
                  clipper::Coord_orth px2c(co00 + frac_x2 * one_step_along_x);
                  glm::vec4 px1_glm = clipper_to_glm(px1c);
                  glm::vec4 px2_glm = clipper_to_glm(px2c);
                  // 2 triangles, 2 vertices as above
                  glm::vec4 col_1 = map_col_func_3(px1c);
                  glm::vec4 col_2 = map_col_func_3(px2c);
                  s_generic_vertex v0(px1_glm, normal, col_1);
                  s_generic_vertex v1(px2_glm, normal, col_2);
                  unsigned  int idx_base = p.first.size();
                  g_triangle gt0(idx_10, idx_base, idx_base + 1);
                  g_triangle gt1(idx_10, idx_11, idx_base);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
#endif
               }
               break;

            case MS_UP_0_0_and_1_0:
               {
                  float frac_y1 = (v00-contour_level)/(v00-v01);
                  float frac_y2 = (v10-contour_level)/(v10-v11);
                  if (frac_y1 < 0.0) std::cout << "ERROR MS_UP_0_0_and_1_0 frac_y1 " << frac_y1 << std::endl;
                  if (frac_y1 > 1.0) std::cout << "ERROR MS_UP_0_0_and_1_0 frac_y1 " << frac_y1 << std::endl;
                  if (frac_y2 < 0.0) std::cout << "ERROR MS_UP_0_0_and_1_0 frac_y2 " << frac_y2 << std::endl;
                  if (frac_y2 > 1.0) std::cout << "ERROR MS_UP_0_0_and_1_0 frac_y2 " << frac_y2 << std::endl;
                  clipper::Coord_orth py1c(co00 + frac_y1 * one_step_along_y);
                  clipper::Coord_orth py2c(co10 + frac_y2 * one_step_along_y);
                  glm::vec4 py1_glm = clipper_to_glm(py1c);
                  glm::vec4 py2_glm = clipper_to_glm(py2c);
                  glm::vec4 col(0.3, 0.0, 0.3, 1.0);
                  glm::vec4 col_1 = map_col_func_3(py1c);
                  glm::vec4 col_2 = map_col_func_3(py2c);
                  s_generic_vertex v0(py1_glm, normal, col_1);
                  s_generic_vertex v1(py2_glm, normal, col_2);
                  unsigned int idx_base = p.first.size();
                  g_triangle gt0(idx_00, idx_10, idx_base + 1);
                  g_triangle gt1(idx_00, idx_base + 1, idx_base);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
               }
               break;

            case MS_UP_0_1_and_1_1: // reverse of above
               {
                  float frac_y1 = (v00-contour_level)/(v00-v01);
                  float frac_y2 = (v10-contour_level)/(v10-v11);
                  if (frac_y1 < 0.0) std::cout << "ERROR MS_UP_0_1_and_1_1 frac_y1 " << frac_y1 << std::endl;
                  if (frac_y1 > 1.0) std::cout << "ERROR MS_UP_0_1_and_1_1 frac_y1 " << frac_y1 << std::endl;
                  if (frac_y2 < 0.0) std::cout << "ERROR MS_UP_0_1_and_1_1 frac_y2 " << frac_y2 << std::endl;
                  if (frac_y2 > 1.0) std::cout << "ERROR MS_UP_0_1_and_1_1 frac_y2 " << frac_y2 << std::endl;
                  clipper::Coord_orth py1c(co00 + frac_y1 * one_step_along_y);
                  clipper::Coord_orth py2c(co10 + frac_y2 * one_step_along_y);
                  glm::vec4 py1_glm = clipper_to_glm(py1c);
                  glm::vec4 py2_glm = clipper_to_glm(py2c);
                  glm::vec4 col_1 = map_col_func_3(py1c);
                  glm::vec4 col_2 = map_col_func_3(py2c);
                  s_generic_vertex v0(py1_glm, normal, col_1);
                  s_generic_vertex v1(py2_glm, normal, col_2);
                  unsigned int idx_base = p.first.size();
                  g_triangle gt0(idx_11, idx_base + 1, idx_base);
                  g_triangle gt1(idx_11, idx_01, idx_base);
                  p.first.push_back(v0);
                  p.first.push_back(v1);
                  p.second.push_back(gt0);
                  p.second.push_back(gt1);
               }
               break;

            case MS_UP_0_0_and_1_1:
               std::cout << "Fill MS_UP_0_0_and_1_1" << std::endl;
               break;

            case MS_UP_0_1_and_1_0:
               std::cout << "Fill MS_UP_1_0_and_0_1" << std::endl;
               break;
            }
         }
      }
   }

   return p;
}

int
molecule_class_info_t::get_square_type(const unsigned int &i,
                                       const unsigned int &j,
                                       const coord_array_2d &arr,
                                       const float & contour_level) const {

   // maybe this could be a member function of the array?

   // c.f lbg-flev.cc:square_type(). We don't include it because that
   // depends on MAKE_ENHANCED_LIGAND_TOOLS

   // this may have the other convention re: i, j indexing! (c.f. ligand_grid contouring)

   int square_type = MS_NO_SQUARE;
   const float &v00 = arr.get_f(i,   j);
   const float &v10 = arr.get_f(i+1, j);
   const float &v01 = arr.get_f(i,   j+1);
   const float &v11 = arr.get_f(i+1, j+1);
   if (v00 < contour_level) {
      if (v01 < contour_level) {
         if (v10 < contour_level) {
            if (v11 < contour_level) {
               return MS_NO_CROSSING;
            }
         }
      }
   }
   if (v00 > contour_level) {
      if (v01 > contour_level) {
         if (v10 > contour_level) {
            if (v11 > contour_level) {
               return MS_ALL_ABOVE;
            }
         }
      }
   }

   // OK, so it is not either of the trivial cases (no
   // crossing), there are 14 (?) other variants.
   // 
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

         // 0,1 is up
	       
         if (v10 < contour_level) { 
            if (v11 < contour_level) {
               return MS_UP_0_1;
            } else {
               return MS_UP_0_1_and_1_1;
            }
         } else {
            if (v11 < contour_level) {
               return MS_UP_0_1_and_1_0;      // hideous valley
            } else {
               return MS_UP_0_1_and_1_0_and_1_1; // (only 0,0 down)
            }
         }
      }
   } else {

      // 0,0 is up
	    
      if (v01 < contour_level) { 
         if (v10 < contour_level) { 
            if (v11 < contour_level) {
               return MS_UP_0_0;
            } else {
               return MS_UP_0_0_and_1_1; // another hideous valley
            }
         } else {
            // 1,0 is up
            if (v11 < contour_level) {
               return MS_UP_0_0_and_1_0;
            } else {
               return MS_UP_0_0_and_1_0_and_1_1; // 0,1 is down
            }
         }
      } else {

         // 0,1 is up
	       
         if (v10 < contour_level) { 
            if (v11 < contour_level) {
               return MS_UP_0_0_and_0_1;
            } else {
               return MS_UP_0_0_and_0_1_and_1_1; // 1,0 is down
            }
         } else {
            // if we get here, this test must pass.
            if (v11 < contour_level) {
               return MS_UP_0_0_and_0_1_and_1_0; // only 1,1 is down
            }
         }
      }
   }
   return square_type;
}
