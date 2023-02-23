// -*-c++-*-
/* coords/Bond_lines_ext.cc
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
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
 * 02110-1301, USA
 */

#include <string>
#include <vector>

#include "Cartesian.h"

#include <mmdb2/mmdb_manager.h> 

#include "mmdb-extras.h"
#include "mmdb-crystal.h"


#include "Bond_lines.h"

#include "Bond_lines_ext.h"


// This is like construct_from_atom_sel, except we don't care about
// dual conformations and atom types.
// 
void 
Bond_lines_ext::find_skel_atom_bonds(atom_selection_container_t SelAtom) {

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
   
   // initialize each colour in the Bond_lines_container
   //
   int col = 0; // atom colour
   
   Bond_lines a(col);
   bonds.push_back(a);

   if (SelAtom.n_selected_atoms <= 0) {
      std::cout << "No skeleton atoms found" << std::endl;
      return;
   }

   mmdb::Contact *contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   // matrix stuff
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;

   // update my_matt;  You can't do this if you haven't set the space group.
   // 
   // symm.GetTMatrix(my_matt, 0);

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   int model_number = 1;

   SelAtom.mol->SeekContacts(SelAtom.atom_selection, SelAtom.n_selected_atoms,
 			     SelAtom.atom_selection, SelAtom.n_selected_atoms,
 			     0.01, 0.7, // min, max distances
 			     0,          // seqDist 0 -> in same res also
 			     contact, ncontacts,
 			     0, &my_matt, i_contact_group);
   
   std::cout << "found " << ncontacts << " bone contacts from "
	 << SelAtom.n_selected_atoms << " selected bone atoms. " <<  std::endl;

   mmdb::PPAtom atom_sel = SelAtom.atom_selection; // shorter form
 
   if (ncontacts > 0) {
      
      for (int i=0; i< ncontacts; i++) {

	 if ( contact[i].id2 >  contact[i].id1 ) {

	    int iat_1 = contact[i].id1;
	    int iat_2 = contact[i].id2;

	    coot::Cartesian atom_1(atom_sel[ contact[i].id1 ]->x,
			     atom_sel[ contact[i].id1 ]->y,
			     atom_sel[ contact[i].id1 ]->z);

	    coot::Cartesian atom_2(atom_sel[ contact[i].id2 ]->x,
			     atom_sel[ contact[i].id2 ]->y,
			     atom_sel[ contact[i].id2 ]->z);

	    addBond(col, atom_1, atom_2, cc, model_number, iat_1, iat_2, true, true);

	 } // contact atom is higher up the list check.
      } // i over ncontacts

      delete [] contact;

   } else {

      std::cout << "There were no skeleton bonds!?" << std::endl;
   }

}

coot::Cartesian 
Bond_lines_ext::find_molecule_middle(atom_selection_container_t SelAtom,
				     float max_neighbour_dist) {

   // find the coordinate that has most bone neighbours.

   coot::Cartesian centre; 
   
   if (SelAtom.n_selected_atoms <= 0) {
      std::cout << "No skeleton atoms found" << std::endl;
      return centre;
   }

   mmdb::Contact *contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   // matrix stuff
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;

   // update my_matt;  You can't do this if you haven't set the space group.
   // 
   // symm.GetTMatrix(my_matt, 0);

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   SelAtom.mol->SeekContacts(SelAtom.atom_selection, SelAtom.n_selected_atoms,
 			     SelAtom.atom_selection, SelAtom.n_selected_atoms,
 			     0.01, max_neighbour_dist, // min, max distances
 			     0,          // seqDist 0 -> in same res also
 			     contact, ncontacts,
 			     0, &my_matt, i_contact_group);
   
   std::cout << "found " << ncontacts << " bone contacts from "
	 << SelAtom.n_selected_atoms << " selected bone atoms. " <<  std::endl;

   mmdb::PPAtom atom_sel = SelAtom.atom_selection; // shorter form
 
   if (ncontacts > 0) {

      std::vector<int> contact_counter(SelAtom.n_selected_atoms, 0);
      
      for (int i=0; i< ncontacts; i++) {

	 contact_counter[contact[i].id1]++; 

      } // i over ncontacts

      delete [] contact;

      // now which atom has the most in contact_counter?
      //
      int max_neighbour_atom_index = -1; // initially a wrong atom index.
      int max_neighbours = 0; 
      
      for(int i=0; i<SelAtom.n_selected_atoms; i++) {
	 if (contact_counter[i] > max_neighbours) {
	    max_neighbour_atom_index = i;
	    max_neighbours = contact_counter[i];
	 }
      }

      if (max_neighbour_atom_index != -1) {

	 //
	 centre.set_them(atom_sel[max_neighbour_atom_index]->x,
			 atom_sel[max_neighbour_atom_index]->y,
			 atom_sel[max_neighbour_atom_index]->z);
      } else {

	 //
	 std::cout << "Pathalogical case in find_molecule_middle.\n"
	      << "WARNING! BAD CENTRE" << std::endl;
      }

   } else {

      std::cout << "There were no skeleton bonds!?" << std::endl;
   }

   return centre;

}

