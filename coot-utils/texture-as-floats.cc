/*
 * coot-utils/texture-as-floats.cc
 *
 * Copyright 2024 by Medical Research Council
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

#include "coot-map-utils.hh"
#include "texture-as-floats.hh"

texture_as_floats_t::texture_as_floats_t(const clipper::Xmap<float> &xmap, int section_index, int axis) {

   // Things we need to fill:
   // int width;
   // int height;
   // float x_size;
   // float y_size;
   // float z_position;
   // std::vector<float> image_data;

   std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
   float mean = mv.first;
   float std_dev = std::sqrt(mv.second);

   // colour map info needs to be passed
   // for now we will do grey-scale
   float data_value_for_top    =  4.0;  // full white
   float data_value_for_bottom = -1.0;  // full black

   data_value_for_top    = mean + 2.5f * std_dev;
   data_value_for_bottom = mean - std_dev;
   float f_range = data_value_for_top - data_value_for_bottom;

   clipper::Cell cell = xmap.cell();
   clipper::Grid_sampling gs = xmap.grid_sampling();
   int new_size = gs.nu() * gs.nv();
   std::cout << "texture_as_floats_t constructor: "  << gs.format()
             << " image data new size " << new_size << std::endl;
   image_data.resize(new_size);
   int image_data_size = image_data.size();

   x_size = cell.a();
   y_size = cell.b();
   if (section_index >= gs.nw()) section_index = gs.nw() -1;
   if (section_index < 0) section_index = 0;
   float z_frac = static_cast<float>(section_index) / static_cast<float>(gs.nw());
   z_position = z_frac * cell.c();
   
   clipper::Coord_grid cg_0(0,0,section_index);
   clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_index);
 
   if (axis == 0) cg_1 = clipper::Coord_grid(section_index, gs.nv()-1, gs.nw()-1); // X
   if (axis == 1) cg_1 = clipper::Coord_grid(gs.nu()-1, section_index, gs.nw()-1); // Y

   clipper::Grid_map grid(cg_0, cg_1);
   clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;
   int nv = gs.nv();

   int c_u = 0;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      int c_v = 0;
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
         int c_w = 0;
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            const float &f = xmap[iw];
            float f_in_range = (f-data_value_for_bottom)/f_range;
            if (f_in_range < 0.0) f_in_range = 0.0;
            if (f_in_range > 1.0) f_in_range = 1.0;
            int idx = c_v + nv * c_u;
            if (idx >= image_data_size) {
               std::cout << "ERROR:: image data index out of range " << idx << " " << image_data_size << std::endl;
            } else {
               image_data[idx] = f_in_range;
            }
            c_w++;
         }
         c_v++;
      }
      c_u++;
   }

   width  = gs.nv();
   height = gs.nu();
}
