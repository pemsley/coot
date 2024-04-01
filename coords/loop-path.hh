/*
 * coords/loop-path.hh
 * 
 * Copyright 2018 by Medical Research Council
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
 */


#ifndef LOOP_PATH_HH
#define LOOP_PATH_HH
#include <mmdb2/mmdb_manager.h>

#include "Cartesian.h"

namespace coot {

   // first is: these points need have bad CA-CA distance spots added
   // second is the vector of line segments
   std::pair<bool, std::vector<coot::CartesianPair> >
   loop_path(mmdb::Atom *start_back_2,
	     mmdb::Atom *start,
	     mmdb::Atom *end,
	     mmdb::Atom *end_plus_2,
	     unsigned int n_line_segments);

   bool is_sane_inter_residue_distance(double dist_between_residues, int res_no_delta, bool is_NA);

}

#endif // LOOP_PATH_HH
