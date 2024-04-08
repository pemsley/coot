/*
 * coot-utils/helix-like.hh
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
 *
 */

#include <mmdb2/mmdb_manager.h>
#include <vector>
#include <string>

namespace coot {

   // select the residues in a chain and call this for every chain
   // in the molecule.
   // return a vector of helical residues in that chain
   //
   std::vector<mmdb::Residue *> like_a_helix(mmdb::Manager *mol, int selection_handle);

   class helical_results_t {
   public:
      bool is_alpha_helix_like;
      bool is_pi_helix_like;
      bool is_3_10_helix_like;
      float sum_delta;
      helical_results_t() {
	 is_alpha_helix_like = false;
	 is_pi_helix_like    = false;
	 is_3_10_helix_like  = false;
	 sum_delta = 0;
      }
   };

   // test_helical_residues should be in the right order. 
   helical_results_t compare_to_helix(const std::vector<mmdb::Residue *> &test_helical_residues,
                                      const std::vector<clipper::Coord_orth> &alpha_helix_ref_positions);

   helical_results_t compare_to_helix(const std::vector<mmdb::Residue *> &helical_residues);

   std::vector<clipper::Coord_orth> alpha_helical_reference_positions();
}
