/*
 * src/updating-map-params.hh
 *
 * Copyright 2019 by Medical Research Council
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

#ifndef UPDATING_MAP_PARAMS_T_HH
#define UPDATING_MAP_PARAMS_T_HH

#include <string>

class updating_map_params_t {
public:
   updating_map_params_t() {
      imol = -1;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
      is_difference_map = false;
      use_weights = false;
   }
   updating_map_params_t(int imol_in,
			 const std::string &mtz_file_name_in,
			 const std::string &f_col_in,
			 const std::string &phi_col_in,
			 const std::string &w_col_in,
			 bool use_weights_in, bool is_difference_map_in) :
      mtz_file_name(mtz_file_name_in), f_col(f_col_in), phi_col(phi_col_in), weight_col(w_col_in)
   {
      imol = imol_in;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
      use_weights = use_weights_in;
      is_difference_map = is_difference_map_in;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
   }
   int imol;
   std::string mtz_file_name;
   std::string f_col;
   std::string phi_col;
   std::string weight_col;
   bool use_weights;
   bool is_difference_map;
   timespec ctime;
};

#endif // UPDATING_MAP_PARAMS_T_HH
