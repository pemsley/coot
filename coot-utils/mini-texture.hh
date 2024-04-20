/*
 * coot-utils/mini-texture.hh
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

#ifndef COOT_SRC_MINI_TEXTURE_HH
#define COOT_SRC_MINI_TEXTURE_HH

#include <clipper/core/xmap.h>

class mini_texture_t {
   public:
   mini_texture_t(int w, int h, unsigned char *d) :
      width(w), height(h), x_size(0), y_size(0), z_position(0), image_data(d),
      data_value_for_top_of_range(0), data_value_for_bottom_of_range(0) {}
   // 0 is for X-sections, 1 is for Y-sections and 2 is for Z-sections
   mini_texture_t(const clipper::Xmap<float> &xmap, int section_index, int axis,
                  float data_value_for_bottom, float data_value_for_top);
   ~mini_texture_t();
   void clear(); // delets image_data
   int width;
   int height;
   float x_size;
   float y_size;
   float z_position;
   unsigned char *image_data;
   float data_value_for_top_of_range;
   float data_value_for_bottom_of_range;
   bool empty() { return image_data == nullptr; }
};

#endif // COOT_SRC_MINI_TEXTURE_HH
