/*
 * geometry/match-results.hh
 *
 * Copyright 2016 by Medical Research Council
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


namespace coot {

   // a container for the results of the comparison vs CCP4SRS graph matching.
   //
   class match_results_t {
   public:
      bool success;
      std::string name;
      std::string comp_id;
      mmdb::Residue *res;
      std::vector<std::pair<int, int> > graph_match_atom_indices;
      // clipper::RTop_orth
      match_results_t(const std::string &comp_id_in, const std::string &name_in, mmdb::Residue *res_in) {
	 name = name_in;
	 comp_id = comp_id_in;
	 res = res_in;
	 if (res_in)
	    success = true;
	 else
	    success = false;
      }
   };

}
