/*
 * coot-utils/strand-like.hh
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

#ifndef STRAND_LIKE_HH
#define STRAND_LIKE_HH

#include <mmdb2/mmdb_manager.h>
#include <vector>
#include <string>

namespace coot {

   // select the residues in a chain and call this for every chain
   // in the molecule.
   void like_a_strand(mmdb::Manager *mol, int selection_handle);

   class strand_results_t {
   public:
      bool is_beta_strand_like;
      bool is_beta_bulge_like;
      float sum_delta;
      helical_results_t() {
	 is_beta_bulge_like  = false;
	 is_beta_strand_like = false;
	 sum_delta = 0;
      }
   };

   // the atoms need to be in the correct order
   strand_results_t compare_to_strand(const std::vector<mmdb::Residue *> &test_residues,
                                      const std::vector<clipper::Coord_orth> &beta_strand_ref_positions);

   strand_results_t compare_to_strand(const std::vector<mmdb::Residue *> &strand_residues);

   std::vector<clipper::Coord_orth> beta_strand_reference_positions();
}

#endif
