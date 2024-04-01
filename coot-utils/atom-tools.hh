/*
 * coot-utils/atom-tools.hh
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

// the atom_selection_container_t is not part of coot-coord-utils.hh

#ifndef ATOM_TOOLS_HH
#define ATOM_TOOLS_HH

#include "atom-selection-container.hh"

namespace coot {

   
   // So that we don't draw the atoms in the "static" molecules that have intermediate atoms.
   // - should make the view more clear.
   //
   std::set<int> atom_indices_in_other_molecule(atom_selection_container_t mol_atom_sel,
						atom_selection_container_t moving_atom_sel);

   void find_out_of_register_errors(mmdb::Manager *post_mutations_mol, mmdb::Manager *ref_mol);

}


#endif // ATOM_TOOLS_HH


