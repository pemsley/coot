
// these are needed to compile molecule-compile-info.h:
// 
#include "Cartesian.h"
#include "mmdb-extras.h"
#include "mmdb-crystal.h"

#include "molecule-class-info.h"

// return an index of the new restraint
int
molecule_class_info_t::add_extra_bond_restraint(coot::atom_spec_t atom_1,
						coot::atom_spec_t atom_2,
						double bond_dist, double esd) {
   CAtom *at_1 = get_atom(atom_1);
   CAtom *at_2 = get_atom(atom_2);
   if (at_1) {
      int atom_index = -1;
      at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_1.int_user_data = atom_index;
   }
   if (at_2) {
      int atom_index = -1;
      at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_2.int_user_data = atom_index;
   }
   coot::extra_restraints_t::extra_bond_restraint_t bond(atom_1, atom_2, bond_dist, esd);
   extra_restraints.bond_restraints.push_back(bond);
   update_extra_restraints_representation();
   return extra_restraints.bond_restraints.size() -1;
}


void molecule_class_info_t::remove_extra_bond_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2) {

   std::vector<coot::extra_restraints_t::extra_bond_restraint_t>::iterator it;
   for (it=extra_restraints.bond_restraints.begin(); it != extra_restraints.bond_restraints.end(); it++) { 
      if (((it->atom_1 == atom_1) &&
	   (it->atom_2 == atom_2)) ||
	  ((it->atom_2 == atom_1) &&
	   (it->atom_1 == atom_2))) {
	 extra_restraints.bond_restraints.erase(it);
	 std::cout << "deleted extra bond restraint " << atom_1 << " to " << atom_2 << std::endl;
	 break;
      }
   }
   update_extra_restraints_representation();
}
