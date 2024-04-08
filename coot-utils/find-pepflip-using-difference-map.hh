/*
 * coot-utils/find-pepflip-using-difference-map.hh
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

#include <vector>
include "geometry/atom-and-residue-specs.hh"

namespace coot {

   class flip_atom_triplet_t {
   public:
      flip_atom_triplet_t(mmdb::Atom *at_1, mmdb::Atom *at_2, mmdb::Atom *at_3) {
	 CA_this = at_1;
	 O_this  = at_2;
	 CA_next = at_3;
      }
      mmdb::Atom *CA_this;
      mmdb::Atom *O_this;
      mmdb::Atom *CA_next;
      clipper::Coord_orth current_O_position() const;
      clipper::Coord_orth flipped_O_position() const;
   };

   class pepflip_using_difference_map {
      mmdb::Manager *mol;
      const clipper::Xmap<float> &diff_xmap;

      std::vector<flip_atom_triplet_t> get_peptide_atom_triplets() const;

      std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > make_random_other_pairs(unsigned int n_others);

   public:
      pepflip_using_difference_map(mmdb::Manager *mol, const clipper::Xmap<float> &xmap) const;

      // return a std::vector of residue specs
      std::vector<residue_spec_t> get_suggested_flips(float n_sigma = 3.0f) const;
   };
}

