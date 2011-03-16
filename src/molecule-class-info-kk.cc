/* ideal/extra-restraints-kk.cc
 * 
 * Copyright 2011 by Kevin Keating
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

#include "Cartesian.h"
#include "mmdb-extras.h"
#include "mmdb-crystal.h"

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
molecule_class_info_t::add_extra_start_pos_restraint(coot::atom_spec_t atom_1,
						double esd) {
   int r = -1; // unset
   CAtom *at_1 = get_atom(atom_1);
   if (at_1) {
      int atom_index = -1;
      at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_1.int_user_data = atom_index;
   
      coot::extra_restraints_t::extra_start_pos_restraint_t start_pos(atom_1, esd);
      extra_restraints.start_pos_restraints.push_back(start_pos);
      update_extra_restraints_representation();
      r = extra_restraints.start_pos_restraints.size() -1;
   }
   return r;
}
