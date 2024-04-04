/*
 * src/mtz-column-auto-read.hh
 *
 * Copyright 2014 by Medical Research Council
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

#include <string>

namespace coot { 
   class mtz_column_trials_info_t {
   public:
      std::string f_col;
      std::string phi_col;
      std::string w_col;
      bool use_weights;
      bool is_diff_map;
      mtz_column_trials_info_t(const std::string &f_col_in,
			       const std::string &phi_col_in,
			       const std::string &w_col_in,
			       bool w_in, bool d_in) {
	 f_col = f_col_in;
	 phi_col = phi_col_in;
	 w_col = w_col_in;
	 use_weights = w_in;
	 is_diff_map = d_in;
      }
      mtz_column_trials_info_t(const std::string &f_col_in,
			       const std::string &phi_col_in,
			       bool d_in) {
	 f_col = f_col_in;
	 phi_col = phi_col_in;
	 w_col = "";
	 use_weights = false;
	 is_diff_map = d_in;
      }
   };
}
