/* ideal/extra-restraints-kk.cc
 * 
 * Copyright 2011, 2012 by Kevin Keating
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

#include "coords/Cartesian.h"
#include "coords/mmdb-extras.h"
#include "coords/mmdb-crystal.h"

#include "molecule-class-info.h"

void molecule_class_info_t::remove_extra_start_pos_restraint(coot::atom_spec_t atom_1) {

   std::vector<coot::extra_restraints_t::extra_start_pos_restraint_t>::iterator it;
   for (it=extra_restraints.start_pos_restraints.begin(); it != extra_restraints.start_pos_restraints.end(); it++) { 
      if (it->atom_1 == atom_1) {
	 extra_restraints.start_pos_restraints.erase(it);
	 std::cout << "deleted extra start_pos restraint " << atom_1 << std::endl;
	 break;
      }
   }
   
   //there is currenctly no representation for start pos restraints
   //update_extra_restraints_representation();
}

void molecule_class_info_t::remove_extra_torsion_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2,
                                                           coot::atom_spec_t atom_3, coot::atom_spec_t atom_4) {

   std::vector<coot::extra_restraints_t::extra_torsion_restraint_t>::iterator it;
   for (it=extra_restraints.torsion_restraints.begin(); it != extra_restraints.torsion_restraints.end(); it++) { 
      if ((it->atom_1 == atom_1) &&
	  (it->atom_2 == atom_2) &&
	  (it->atom_3 == atom_3) &&
	  (it->atom_4 == atom_4)) {
	 extra_restraints.torsion_restraints.erase(it);
	 std::cout << "deleted extra torsion restraint " << atom_1 << " to " << atom_2 << " to " << atom_3 << " to " << atom_4 << std::endl;
	 break;
      }
   }
   
   //there is currenctly no representation for torsion restraints
   //update_extra_restraints_representation();
}

// return an index of the new restraint
int
molecule_class_info_t::add_extra_angle_restraint(coot::atom_spec_t atom_1,
						   coot::atom_spec_t atom_2,
						   coot::atom_spec_t atom_3,
						   double angle, double esd) {

   CAtom *at_1 = get_atom(atom_1);
   CAtom *at_2 = get_atom(atom_2);
   CAtom *at_3 = get_atom(atom_3);
   if (at_1) {
      int atom_index = -1;
      at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_1.int_user_data = atom_index;
   }
   if (at_2) {
      int atom_index = -1;
      at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_2.int_user_data = atom_index;
   }
   if (at_3) {
      int atom_index = -1;
      at_3->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_3.int_user_data = atom_index;
   }
   coot::extra_restraints_t::extra_angle_restraint_t ang(atom_1, atom_2,
							    atom_3,
							    angle, esd);
   extra_restraints.angle_restraints.push_back(ang);
   update_extra_restraints_representation();
   return extra_restraints.angle_restraints.size() -1;
}

void molecule_class_info_t::remove_extra_angle_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2, coot::atom_spec_t atom_3) {

   std::vector<coot::extra_restraints_t::extra_angle_restraint_t>::iterator it;
   for (it=extra_restraints.angle_restraints.begin(); it != extra_restraints.angle_restraints.end(); it++) { 
      if (((it->atom_1 == atom_1) &&
	   (it->atom_2 == atom_2) &&
           (it->atom_3 == atom_3)) ||
	  ((it->atom_3 == atom_1) &&
	   (it->atom_2 == atom_2) &&
           (it->atom_1 == atom_3))) {
	 extra_restraints.angle_restraints.erase(it);
	 std::cout << "deleted extra angle restraint " << atom_1 << " to " << atom_2 << " to " << atom_3 << std::endl;
	 break;
      }
   }
   //there is currenctly no representation for angle restraints
   //update_extra_restraints_representation();
}

// return an index of the new restraint
int
molecule_class_info_t::add_extra_start_pos_restraint(coot::atom_spec_t atom_1,
						double esd) {
   int r = -1; // unset
   CAtom *at_1 = get_atom(atom_1);
   if (at_1) {
      int atom_index = -1;
      at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_1.int_user_data = atom_index;
   
      // updates restraint on atom if it can, else adds
      extra_restraints.add_start_pos_restraint(coot::atom_spec_t(atom_1), esd); 
      update_extra_restraints_representation();
      r = extra_restraints.start_pos_restraints.size() -1;
   }
   return r;
}
