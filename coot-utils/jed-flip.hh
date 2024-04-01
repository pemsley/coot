/*
 * coot-utils/jed-flip.hh
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


#ifndef JED_FLIP_HH
#define JED_FLIP_HH

#include "geometry/protein-geometry.hh"
#include "coot-coord-extras.hh"
#include "atom-tree.hh"

namespace coot {

   namespace util {

      // return a diagnostic message (set if needed)
      std::string jed_flip(int imol,
			   mmdb::Residue *residue_p, mmdb::Atom *clicked_atom,
			   bool invert_selection, 
			   protein_geometry *geom);

      std::string
      jed_flip_internal(coot::atom_tree_t &tree,
			const std::vector<dict_torsion_restraint_t> &interesting_torsions,
			const std::string &atom_name,
			int atom_idx,
			bool invert_selection);

      // return a non-null string on a problem
      //
      std::string
      jed_flip_internal(atom_tree_t &tree,
			const dict_torsion_restraint_t &torsion,
			const std::string &atom_name,
			int clicked_atom_idx,
			bool invert_selection);

   }
}

#endif // JED_FLIP_HH
