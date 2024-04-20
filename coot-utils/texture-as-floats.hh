/*
 * coot-utils/texture-as-floats.hh
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

#ifndef COOT_TEXTURE_AS_FLOATS
#define COOT_TEXTURE_AS_FLOATS

#include <clipper/core/xmap.h>

class texture_as_floats_t {
   public:
   texture_as_floats_t() : width(0), height(0), x_size(0), y_size(0), z_position(0) {}
   //! axis_id is 0 for X-sections, 1 for Y-sections and 2 for Z-sections
   texture_as_floats_t(const clipper::Xmap<float> &xmap, int section_index, int axis_id);
   int width;
   int height;
   float x_size;
   float y_size;
   float z_position;
   std::vector<float> image_data;
   bool empty() { return image_data.empty(); }
};

#endif // COOT_TEXTURE_AS_FLOATS
