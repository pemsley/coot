/*
 * coot-utils/bonded-atoms.hh
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

#include <vector>
#include <mmdb2/mmdb_manager.h>

// I want a fast way of knowing if atoms are bonded
// (for non-bonded contacts)
// so instead of using a dictionary, let's try this
// by distance within the same residue and by
// hard-coded atom name comparison between residues.
//
// RNA, DNA and polypeptides

namespace coot {

   // this should only be called for atoms in the same residue
   // or atoms that are next to each other in sequence
   bool are_polymer_bonded(mmdb::Atom *at_1, mmdb::Atom *at_2);

   std::vector<std::vector<unsigned int> > make_bonds(mmdb::Manager *mol, int n_selected_atoms, int mol_atom_index_handle);

   // do I want a version that works with an atom selection?

   std::vector<std::vector<unsigned int> >
   find_1_4_connections(const std::vector<std::vector<unsigned int> > &bonds_vec);
   
}
