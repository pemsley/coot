
#include <map>
#include <algorithm>

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
   int r = -1; // unset
   CAtom *at_1 = get_atom(atom_1);
   CAtom *at_2 = get_atom(atom_2);
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
   if (at_1 && at_2) { 
      coot::extra_restraints_t::extra_bond_restraint_t bond(atom_1, atom_2, bond_dist, esd);
      extra_restraints.bond_restraints.push_back(bond);
      update_extra_restraints_representation();
      r = extra_restraints.bond_restraints.size() -1;
   }
   return r;
}


void molecule_class_info_t::remove_extra_bond_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2) {

   int n_bonds_pre = extra_restraints.bond_restraints.size();
   std::vector<coot::extra_restraints_t::extra_bond_restraint_t>::iterator it;
   for (it=extra_restraints.bond_restraints.begin(); it != extra_restraints.bond_restraints.end(); it++) { 
      if (((it->atom_1 == atom_1) &&
	   (it->atom_2 == atom_2)) ||
	  ((it->atom_2 == atom_1) &&
	   (it->atom_1 == atom_2))) {
	 extra_restraints.bond_restraints.erase(it);
	 std::cout << "deleted extra bond restraint " << atom_1 << " to " << atom_2 << std::endl;
	 // break;
      }
   }
   int n_bonds_post = extra_restraints.bond_restraints.size();
   std::cout << "DEBUG:: pre: " << n_bonds_pre << " post " << n_bonds_post << std::endl;
   update_extra_restraints_representation();
}

void
molecule_class_info_t::add_refmac_extra_restraints(const std::string &file_name) {

   coot::extra_restraints_t r;
   r.read_refmac_distance_restraints(file_name);
   extra_restraints.add_restraints(r);
   std::cout << "in add_refmac_extra_restraints we have " << extra_restraints.bond_restraints.size()
	     << " bond restraints " << std::endl;
   update_extra_restraints_representation();
}

void
molecule_class_info_t::delete_extra_restraints_for_residue(const coot::residue_spec_t &rs) {

   unsigned int pre_n = extra_restraints.bond_restraints.size(); 
   extra_restraints.delete_restraints_for_residue(rs);
   unsigned int post_n = extra_restraints.bond_restraints.size();
   if (post_n != pre_n)
      update_extra_restraints_representation();
}

void
molecule_class_info_t::delete_extra_restraints_worse_than(const double &n_sigma) {

   unsigned int pre_n = extra_restraints.bond_restraints.size();

   // the real dist, with atom specs for keys
   // (to be used by the erasor)
   // 
   std::map<std::pair<coot::atom_spec_t, coot::atom_spec_t>, double> dist_map;
   std::map<coot::atom_spec_t, CAtom *> atom_map;
   std::map<coot::atom_spec_t, CAtom *>::const_iterator it_1;
   std::map<coot::atom_spec_t, CAtom *>::const_iterator it_2;
   // first fill the dist_map and fill the atom_map as you do so.
   for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {
      coot::extra_restraints_t::extra_bond_restraint_t &br = extra_restraints.bond_restraints[i];
      CAtom *at_1 = NULL;
      CAtom *at_2 = NULL;
      it_1 = atom_map.find(br.atom_1);
      it_2 = atom_map.find(br.atom_2);
      if (it_1 == atom_map.end()) {
	 at_1 = get_atom(br.atom_1);
	 atom_map[br.atom_1] = at_1;
      } else {
	 at_1 = it_1->second; // most of the hits, I hope
      }
      // and the same for 2:
      if (it_2 == atom_map.end()) {
	 at_2 = get_atom(br.atom_2);
	 atom_map[br.atom_2] = at_2;
      } else {
	 at_2 = it_2->second; 
      }
      if (at_1 && at_2) {
	 double dx = at_1->x - at_2->x;
	 double dy = at_1->y - at_2->y;
	 double dz = at_1->z - at_2->z;
	 double dist_sq = dx*dx + dy*dy + dz*dz;
	 if (dist_sq < 0) dist_sq = 0;
	 dist_map[std::pair<coot::atom_spec_t, coot::atom_spec_t>(br.atom_1, br.atom_2)] = sqrt(dist_sq);
      }
   }

   extra_restraints.bond_restraints.erase(std::remove_if(extra_restraints.bond_restraints.begin(), extra_restraints.bond_restraints.end(), coot::extra_restraints_t::bond_erasor(dist_map, n_sigma)), extra_restraints.bond_restraints.end());

   unsigned int post_n = extra_restraints.bond_restraints.size();
   if (post_n != pre_n)
      update_extra_restraints_representation();
   std::cout << "INFO deleted : " << pre_n - post_n << " extra bond restraints" << std::endl;
} 


void
molecule_class_info_t::clear_extra_restraints() {
   extra_restraints.clear();
} 

// return an index of the new restraint
int
molecule_class_info_t::add_extra_torsion_restraint(coot::atom_spec_t atom_1,
						   coot::atom_spec_t atom_2,
						   coot::atom_spec_t atom_3,
						   coot::atom_spec_t atom_4,
						   double torsion_angle, double esd, int period) {

   CAtom *at_1 = get_atom(atom_1);
   CAtom *at_2 = get_atom(atom_2);
   CAtom *at_3 = get_atom(atom_3);
   CAtom *at_4 = get_atom(atom_4);
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
   if (at_4) {
      int atom_index = -1;
      at_4->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_4.int_user_data = atom_index;
   }
   coot::extra_restraints_t::extra_torsion_restraint_t tors(atom_1, atom_2,
							    atom_3,atom_4,
							    torsion_angle, esd, period);
   extra_restraints.torsion_restraints.push_back(tors);
   update_extra_restraints_representation();
   return extra_restraints.torsion_restraints.size() -1;
}
