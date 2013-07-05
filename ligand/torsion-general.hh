/* ligand/torsion-general.hh
 * 
 * Copyright 2007 by The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms ofn the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */


#include <vector>
#include <string>
#include "coot-utils/coot-coord-utils.hh"

#ifndef __MMDB_Manager__
#include <mmdb_manager.h>
#endif
#include "ccp4mg-utils/mgtree.h"

namespace coot {

   class torsion_general {
      enum setup_params { NO, YES, FAIL=-1 };
      bool setup_correctly;
      CResidue *residue_p;
      CMMDBManager *mol;
      std::vector<atom_spec_t> user_defined_torsion_atoms;
      std::vector<int> clicked_atom_indices;
      std::vector<std::vector<int> > contact_indices;
      // return the atom index in residue of the atom with the given
      // spec.  Return -1 on failure to find the atom.
      int atom_index(const coot::atom_spec_t &spec) const;
      std::vector<std::vector<int> > get_contact_indices() const;
   public:
      torsion_general(CResidue *res, CMMDBManager *residue_mol_in,
		      const std::vector<atom_spec_t> &user_defined_torsion_atoms_in);
      int change_by(double diff, Tree *tree); // tree is modified
      Tree GetTree() const;
   };

} 
