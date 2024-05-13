/*
 * coot-utils/mini-texture.cc
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

#include <chrono>
#include "coot-utils/coot-map-utils.hh"
#include "mini-texture.hh"

mini_texture_t::mini_texture_t(const clipper::Xmap<float> &xmap,
                               int section_index, int axis,
                               float data_value_for_bottom, float data_value_for_top) {

   auto tp_start = std::chrono::high_resolution_clock::now();

   data_value_for_top_of_range    = data_value_for_top;
   data_value_for_bottom_of_range = data_value_for_bottom;

   float f_range = data_value_for_top - data_value_for_bottom;
   const float inv_f_range = 1.0/f_range; // because we divide by the f_range;

   clipper::Grid_sampling gs = xmap.grid_sampling();
   // std::cout << "mini_texture_t  constructor: "  << gs.format() << std::endl;
   int image_data_size = -1;
   if (axis == 0) image_data_size = gs.nv() * gs.nw() * 4;
   if (axis == 1) image_data_size = gs.nu() * gs.nw() * 4;
   if (axis == 2) image_data_size = gs.nu() * gs.nv() * 4;

   if (false) {
      if (axis == 0) std::cout << "axis 0: image_data_size 4 * " << gs.nv() << " * " << gs.nw() << std::endl;
      if (axis == 1) std::cout << "axis 1: image_data_size 4 * " << gs.nu() << " * " << gs.nw() << std::endl;
      if (axis == 2) std::cout << "axis 2: image_data_size 4 * " << gs.nu() << " * " << gs.nv() << std::endl;
   }

   image_data = new unsigned char[image_data_size];
   if (true) { // pre-purple
      for (int ii=0; ii<image_data_size; ii+=4) {
         image_data[ii]   = 66;
         image_data[ii+1] = 0;
         image_data[ii+2] = 66;
         image_data[ii+3] = 255;
      }
   }
   // auto tp_t_id2 = std::chrono::high_resolution_clock::now();
   // auto deltaid = std::chrono::duration_cast<std::chrono::milliseconds>(tp_t_id2 - tp_t_id1);

   clipper::Cell cell = xmap.cell();
   x_size = -1;
   y_size = -1;

   float z_frac = static_cast<float>(section_index) / static_cast<float>(gs.nw());
   z_position = z_frac * cell.c();

   clipper::Coord_grid cg_0;
   clipper::Coord_grid cg_1;
   if (axis == 0) {
      int section_index_X = z_frac * gs.nu();
      cg_0 = clipper::Coord_grid(section_index_X, 0, 0);
      cg_1 = clipper::Coord_grid(section_index_X, gs.nv()-1, gs.nw()-1); // X
      x_size = cell.c();
      y_size = cell.b();
   }
   if (axis == 1) {
      int section_index_Y = z_frac * gs.nv();
      cg_0 = clipper::Coord_grid(0, section_index_Y, 0);
      cg_1 = clipper::Coord_grid(gs.nu()-1, section_index_Y, gs.nw()-1); // Y
      x_size = cell.a();
      y_size = cell.c();
   }
   if (axis == 2) {
      int section_index_Z = z_frac * gs.nw();
      cg_0 = clipper::Coord_grid(0, 0, section_index_Z);
      cg_1 = clipper::Coord_grid(gs.nu()-1, gs.nv()-1, section_index); // Z
      x_size = cell.a();
      y_size = cell.b();
   }

   clipper::Grid_map grid(cg_0, cg_1);
   int nu = gs.nu();
   int nv = gs.nv();
   int nw = gs.nw();

   if (false)
      std::cout << "here with axis " << axis << " "
                << " grid " << grid.format() << " from " << cg_0.format() << " " << cg_1.format() << std::endl;

   auto tp_pre_grid = std::chrono::high_resolution_clock::now();
   clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;
   int c_u = 0;
   for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u()) {
      int c_v = 0;
      for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v()) {
         int c_w = 0;
         for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w()) {
            const float &f = xmap[iw];
            float f_in_range = (f-data_value_for_bottom) * inv_f_range;
            if (f_in_range < 0.0) f_in_range = 0.0;
            if (f_in_range > 1.0) f_in_range = 1.0;
            int img_x_coord = -1;
            int img_y_coord = -1;
            int img_n_rows  = -1;
            if (axis == 0) { img_x_coord = c_v; img_y_coord = c_w; img_n_rows = nw; }
            if (axis == 1) { img_x_coord = c_u; img_y_coord = c_w; img_n_rows = nw; }
            if (axis == 2) { img_x_coord = c_u; img_y_coord = c_v; img_n_rows = nv; }
            int idx = 4 * (img_y_coord + img_n_rows * img_x_coord);
            if (axis == 1111110) {
               std::cout << "debug:: " << c_v << " " << c_w << " " << idx/4 << std::endl;
            }
            if (idx >= image_data_size) {
               std::cout << "mini_texture_t(): out of index " << idx << " " << image_data_size
                         << " for axis " << axis << " coords " << img_x_coord << " " << img_y_coord << " row " << img_n_rows << "\n";
            } else {
               if (idx < 0) continue;
               image_data[idx]   = static_cast<unsigned char>(255 * f_in_range);
               image_data[idx+1] = static_cast<unsigned char>(255 * f_in_range);
               image_data[idx+2] = static_cast<unsigned char>(255 * f_in_range);
               image_data[idx+3] = 255;
               // if (image_data[idx] > 200) { std::cout << idx << std::endl;}
               // std::cout << idx << " " << iw.coord().format() << " " << f << " " << f_in_range << "\n";
               if (axis == 110) {
                  // for (int ii=0; ii<4; ii++) image_data[idx+ii] = 255;
                  float f_1 = static_cast<float>(c_u)/static_cast<float>(gs.nu());
                  float f_2 = static_cast<float>(c_v)/static_cast<float>(gs.nv());
                  float f_3 = static_cast<float>(c_w)/static_cast<float>(gs.nw());
                  f_1 = z_frac;
                  image_data[idx]   = static_cast<unsigned char>(255 * f_1);
                  image_data[idx+1] = static_cast<unsigned char>(255 * f_2);
                  image_data[idx+2] = static_cast<unsigned char>(255 * f_3);
                  image_data[idx+3] = static_cast<unsigned char>(255);
               }
            }
            c_w++;
         }
         c_v++;
      }
      c_u++;
   }

   if (axis == 0) { width = gs.nw(); height = gs.nv(); } // correct for index colour
   if (axis == 1) { width = gs.nw(); height = gs.nv(); }
   if (axis == 2) { width = gs.nv(); height = gs.nu(); }

   if (false) {
      auto tp_post_grid = std::chrono::high_resolution_clock::now();
      auto delta2  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_post_grid - tp_pre_grid);
      auto deltaa  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_post_grid - tp_start);
      // auto deltatr = std::chrono::duration_cast<std::chrono::milliseconds>(tp_t_1 - tp_start);
      std::cout << "constructor all:             " << deltaa.count()   << " ms" << std::endl;
      std::cout << "constructor grid time:       " << delta2.count()   << " ms" << std::endl;
      // std::cout << "constructor image-data time: " << deltaid.count()  << " ms" << std::endl;
      // std::cout << "constructor test_range:      " << deltatr.count()  << " ms" << std::endl;
   }

}

mini_texture_t::mini_texture_t(const texture_as_floats_t &taf, float data_value_for_bot, float data_value_for_top) {
   width = taf.width;
   height = taf.height;
   x_size = taf.x_size;
   y_size = taf.y_size;
   z_position = 0.0;
   data_value_for_bottom_of_range = data_value_for_bot;
   data_value_for_top_of_range = data_value_for_top;

   data_value_for_top_of_range    = data_value_for_top;
   data_value_for_bottom_of_range = data_value_for_bot;

   float f_range = data_value_for_top - data_value_for_bot;
   const float inv_f_range = 1.0/f_range; // because we divide by the f_range;

   size_t image_data_size = 4 * width * height;
   image_data = new unsigned char[image_data_size];
   for (size_t i = 0; i < image_data_size; i += 4) {
      const float &f = taf.image_data[i/4];
      float f_in_range = (f-data_value_for_bot) * inv_f_range;
      if (f_in_range < 0.0) f_in_range = 0.0;
      if (f_in_range > 1.0) f_in_range = 1.0;
      image_data[i]   = static_cast<unsigned char>(255 * f_in_range);
      image_data[i+1] = static_cast<unsigned char>(255 * f_in_range);
      image_data[i+2] = static_cast<unsigned char>(255 * f_in_range);
      image_data[i+3] = 255;
   }
}

mini_texture_t::~mini_texture_t() {
   // delete [] image_data;
}

void
mini_texture_t::clear() {
   delete [] image_data;
}
