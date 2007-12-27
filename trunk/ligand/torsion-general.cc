/* ligand/torsion-general.cc
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


#include "torsion-general.hh"
#include "mgtree.h"

coot::torsion_general::torsion_general(CResidue *res, CMMDBManager *residue_mol_in,
				       const std::vector<atom_spec_t> &user_defined_torsion_atoms_in) {
   setup_correctly = NO;
   mol = residue_mol_in;
   residue_p = res;
   user_defined_torsion_atoms = user_defined_torsion_atoms_in;

   int index_0 = atom_index(user_defined_torsion_atoms[0]);
   int index_1 = atom_index(user_defined_torsion_atoms[1]);
   int index_2 = atom_index(user_defined_torsion_atoms[2]);
   int index_3 = atom_index(user_defined_torsion_atoms[3]);

   if (index_0 != FAIL) { 
      if (index_1 != FAIL) { 
	 if (index_2 != FAIL) { 
	    if (index_3 != FAIL) {
	       contact_indices = get_contact_indices();
	       clicked_atom_indices.clear();
	       clicked_atom_indices.push_back(index_0);
	       clicked_atom_indices.push_back(index_1);
	       clicked_atom_indices.push_back(index_2);
	       clicked_atom_indices.push_back(index_3);
	       setup_correctly = 1;
	    } else {
	       std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms[3] << std::endl;
	    }
	 } else {
	    std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms[2] << std::endl;
	 }
      } else {
	 std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms[1] << std::endl;
      }
   } else {
      std::cout << "ERROR:: failed to find " << user_defined_torsion_atoms[0] << std::endl;
   }
}


int
coot::torsion_general::change_by(double diff) {

   int r=1;
//    std::cout << "user_defined_torsion_atoms.size(): " << user_defined_torsion_atoms.size()
// 	     << std::endl;
   if (setup_correctly) {

      PPCAtom residue_atoms;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector< ::Cartesian > coords;
      for(int i=0; i<n_residue_atoms; i++) {
	 ::Cartesian c(residue_atoms[i]->x,
		       residue_atoms[i]->y,
		       residue_atoms[i]->z);
	 coords.push_back(c);
      }
      int base_index = clicked_atom_indices[0];
      Tree tree;
      tree.SetCoords(coords, base_index, contact_indices);
      TreeVertex *tv = tree.GetCoord(clicked_atom_indices[1]);
      if (tv->GetNumberOfChildren() > 0) {
	 TreeVertex *tvc0 = tv->GetChild(0);
	 float tors = clipper::Util::d2rad(diff);
	 tree.RotateAboutBond(clicked_atom_indices[1],
			      clicked_atom_indices[2], tors);
      
	 std::vector< ::Cartesian > coords_rotatated =
	    tree.GetAllCartesians();
	 if (int(coords_rotatated.size()) != n_residue_atoms) {
	    std::cout << "disaster in atom selection, tors_general\n";
	 } else {
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       residue_atoms[iat]->x = coords_rotatated[iat].get_x();
	       residue_atoms[iat]->y = coords_rotatated[iat].get_y();
	       residue_atoms[iat]->z = coords_rotatated[iat].get_z();
	    }
	    r = 0; // return good status
	 }
      } else {
	 std::cout << "WARNING: this vertex " << clicked_atom_indices[2]
		   << " has no children (strangely)\n";
      } 
   } else {
      std::cout << "Sorry torsion_general not setup correctly" << std::endl;
   }
   return r;
}

// return the atom index in residue of the atom with the given
// spec.  Return -1 on failure to find the atom.
int
coot::torsion_general::atom_index(const coot::atom_spec_t &spec) const {

   int r = FAIL;

   if (residue_p) {
      PPCAtom residue_atoms;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 if (spec.matches_spec(residue_atoms[i])) {
	    // 	    std::cout << "Found a spec match for " << spec << std::endl;
	    return i;
// 	 } else {
// 	    std::cout << spec << " does not match " << coot::atom_spec_t(residue_atoms[i]) << std::endl;
	 }
      }
   }
   return r;
}

std::vector<std::vector<int> >
coot::torsion_general::get_contact_indices() const {

   std::vector<std::vector<int> > v;
   PSContact pscontact = NULL;
   int n_contacts;
   float min_dist = 0.1;
   float max_dist = 1.9; // CB->SG CYS 1.8A
   if (std::string(residue_p->GetResName()) == "MSE")
      max_dist = 2.0;
   long i_contact_group = 1;
   mat44 my_matt;
   CSymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   PPCAtom residue_atoms;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   
   mol->SeekContacts(residue_atoms, n_residue_atoms,
		     residue_atoms, n_residue_atoms,
		     min_dist, max_dist, // min, max distances
		     0,        // seqDist 0 -> in same res also
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   v.resize(n_residue_atoms);
   for (int ic=0; ic< n_contacts; ic++) {
      // std::cout << pscontact[ic].id1 << " " << pscontact[ic].id2 << std::endl;
      v[pscontact[ic].id1].push_back(pscontact[ic].id2);
      v[pscontact[ic].id2].push_back(pscontact[ic].id1);
   }
   return v;
}
