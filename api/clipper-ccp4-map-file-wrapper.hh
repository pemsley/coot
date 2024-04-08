/*
 * api/clipper-ccp4-map-file-wrapper.hh
 * 
 * Copyright 2020 by Medical Research Council
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


#include <clipper/ccp4/ccp4_map_io.h>

class clipper_map_file_wrapper : public clipper::CCP4MAPfile {
public:
   clipper_map_file_wrapper() : clipper::CCP4MAPfile() { }
   void wrap_open_read(const clipper::String &filename_in) {
      open_read(filename_in);
   }
   bool starts_at_zero() const {
      if (grid_map_.min() == clipper::Coord_grid(0,0,0)) {
	 return true;
      } else {
	 return false;
      }
   }
};

