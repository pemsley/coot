/* src/molecule-class-info-refine.cc
 *
 * Copyright 2010, 2011, 2012 by the University of Oxford
 * Copyright 2013 by Medical Research Council
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include "compat/coot-sysdep.h"

#include <map>
#include <algorithm>

// these are needed to compile molecule-compile-info.h:
//
#include "coords/Cartesian.hh"
#include "coords/mmdb-extras.hh"
#include "coords/mmdb-crystal.hh"

// morphing
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/xmap-stats.hh"
#include "ligand/rigid-body.hh"

#include "molecule-class-info.h"

class atom_morph_info_t {
public:
   double weight;
   clipper::Coord_orth residue_centre;
   clipper::RTop_orth residue_rtop;
   atom_morph_info_t(const clipper::RTop_orth &rtop_in,
		     const clipper::Coord_orth &co_in,
		     double w_in) {
      weight = w_in;
      residue_rtop = rtop_in;
      residue_centre = co_in;
   }
   std::pair<clipper::RTop_orth, float> rtop_and_weight() const {
      return std::pair<clipper::RTop_orth, float> (residue_rtop, weight);
   }
};

// return an index of the new restraint
int
molecule_class_info_t::add_extra_bond_restraint(coot::atom_spec_t atom_1,
						coot::atom_spec_t atom_2,
						double bond_dist, double esd) {
   int r = -1; // unset
   mmdb::Atom *at_1 = get_atom(atom_1);
   mmdb::Atom *at_2 = get_atom(atom_2);
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
   } else {
      std::cout << "WARNING:: add_extra_bond_restraint() oops: " << at_1 << " " << atom_1 << " "
		<< at_2 << " " << atom_2 << std::endl;
   }
   return r;
}

//! arguments are modified, so they are not const.
int
molecule_class_info_t::add_extra_geman_mcclure_restraint(coot::atom_spec_t atom_1,
                                                         coot::atom_spec_t atom_2,
                                                         double bond_dist, double esd) {
   int r = -1; // unset
   mmdb::Atom *at_1 = get_atom(atom_1);
   mmdb::Atom *at_2 = get_atom(atom_2);
   int atom_index_1 = -1;
   int atom_index_2 = -1;
   if (at_1) {
      at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index_1); // set atom_index
      atom_1.int_user_data = atom_index_1;
   }
   if (at_2) {
      at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index_2); // set atom_index
      atom_2.int_user_data = atom_index_2;
   }
   if (at_1 && at_2) {
      coot::extra_restraints_t::extra_geman_mcclure_restraint_t bond(atom_1, atom_2, bond_dist, esd);
      extra_restraints.geman_mcclure_restraints.push_back(bond);
      update_extra_restraints_representation();
      r = extra_restraints.geman_mcclure_restraints.size() -1;
   } else {
      std::cout << "WARNING:: add_extra_geman_mcclure_restraint() oops: " << at_1 << " " << atom_1 << " "
		<< at_2 << " " << atom_2 << std::endl;
   }
   return r;
}


int
molecule_class_info_t::add_extra_bond_restraints(const std::vector<coot::extra_restraints_t::extra_bond_restraint_t> &bond_specs) {

   int r = -1; // unset
   for (unsigned int i=0; i<bond_specs.size(); i++) {
      coot::extra_restraints_t::extra_bond_restraint_t bond_spec = bond_specs[i]; // gets modified by addition of atom indices
      mmdb::Atom *at_1 = get_atom(bond_spec.atom_1);
      mmdb::Atom *at_2 = get_atom(bond_spec.atom_2);
      if (at_1) {
	 int atom_index = -1;
	 at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
	 bond_spec.atom_1.int_user_data = atom_index;
      }
      if (at_2) {
	 int atom_index = -1;
	 at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
	 bond_spec.atom_2.int_user_data = atom_index;
      }
      if (at_1 && at_2) {
	 // coot::extra_restraints_t::extra_bond_restraint_t bond(bond_spec.atom_1_spec, bond_spec.atom_2_spec, bond_spec.dist, bond_spec.esd);
	 extra_restraints.bond_restraints.push_back(bond_spec);
	 r = extra_restraints.bond_restraints.size() -1;
      } else {
	 std::cout << "WARNING:: add_extra_bond_restraint() oops: " << at_1 << " " << bond_spec.atom_1 << " "
		   << at_2 << " " << bond_spec.atom_2 << std::endl;
      }
   }
   update_extra_restraints_representation();
   return r;
}

int
molecule_class_info_t::add_extra_geman_mcclure_restraints(const std::vector<coot::extra_restraints_t::extra_geman_mcclure_restraint_t> &bond_specs) {

   int r = -1; // unset
   for (unsigned int i=0; i<bond_specs.size(); i++) {
      coot::extra_restraints_t::extra_geman_mcclure_restraint_t bond_spec = bond_specs[i]; // gets modified by addition of atom indices
      mmdb::Atom *at_1 = get_atom(bond_spec.atom_1);
      mmdb::Atom *at_2 = get_atom(bond_spec.atom_2);
      int atom_index_1 = -1;
      int atom_index_2 = -1;
      if (at_1) {
	 at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index_1); // set atom_index
	 bond_spec.atom_1.int_user_data = atom_index_1;
      }
      if (at_2) {
	 at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index_2); // set atom_index
	 bond_spec.atom_2.int_user_data = atom_index_2;
      }
      if (at_1 && at_2) {
	 extra_restraints.geman_mcclure_restraints.push_back(bond_spec);
	 r = extra_restraints.geman_mcclure_restraints.size() -1;
      } else {
	 std::cout << "WARNING:: add_extra_bond_restraint() oops: " << at_1 << " "
		   << bond_spec.atom_1 << " " << at_2 << " " << bond_spec.atom_2 << std::endl;
      }
   }
   update_extra_restraints_representation();
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
	 break; // this break was commented out, but it must not be, because the iterator
	        // is no longer valid after we erase from a vector.
      }
   }
   int n_bonds_post = extra_restraints.bond_restraints.size();
   std::cout << "DEBUG:: pre: " << n_bonds_pre << " post " << n_bonds_post << std::endl;
   update_extra_restraints_representation();
}


void
molecule_class_info_t::remove_extra_geman_mcclure_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2) {

   int n_bonds_pre = extra_restraints.geman_mcclure_restraints.size();
   std::vector<coot::extra_restraints_t::extra_geman_mcclure_restraint_t>::iterator it;
   for (it=extra_restraints.geman_mcclure_restraints.begin(); it != extra_restraints.geman_mcclure_restraints.end(); it++) {
      if (((it->atom_1 == atom_1) &&
	   (it->atom_2 == atom_2)) ||
	  ((it->atom_2 == atom_1) &&
	   (it->atom_1 == atom_2))) {
	 extra_restraints.geman_mcclure_restraints.erase(it);
	 std::cout << "deleted extra bond restraint " << atom_1 << " to " << atom_2 << std::endl;
	 break; // this break was commented out, but it must not be, because the iterator
	        // is no longer valid after we erase from a vector.
      }
   }
   int n_bonds_post = extra_restraints.geman_mcclure_restraints.size();
   std::cout << "DEBUG:: pre: GM bonds " << n_bonds_pre << " post " << n_bonds_post << std::endl;
   update_extra_restraints_representation();
}

void
molecule_class_info_t::add_refmac_extra_restraints(const std::string &file_name) {

   coot::extra_restraints_t r;
   r.read_refmac_extra_restraints(file_name);
   extra_restraints.add_restraints(r);
   std::cout << "INFO:: add_refmac_extra_restraints(): have "
	     << extra_restraints.bond_restraints.size() << " extra bond restraints " << std::endl;
   std::cout << "INFO:: add_refmac_extra_restraints(): have "
	     << extra_restraints.angle_restraints.size() << " extra angle restraints " << std::endl;
   std::cout << "INFO:: add_refmac_extra_restraints(): have "
	     << extra_restraints.torsion_restraints.size() << " extraa torsion restraints " << std::endl;
   update_extra_restraints_representation();
}

// extra target position restraints refine like pull atom restraints
// but are stored differently (like other extra restraints)
int
molecule_class_info_t::add_extra_target_position_restraint(coot::atom_spec_t &spec,
							   const clipper::Coord_orth &pos,
							   float weight) {
   int r = -1;
   mmdb::Atom *at = get_atom(spec);
   if (at) {
      int atom_index = -1;
      at->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      spec.int_user_data = atom_index;
      coot::extra_restraints_t::extra_target_position_restraint_t tpr(spec, pos, weight);
      std::cout << "debug:: adding target position restraint for " << spec << std::endl;
      extra_restraints.target_position_restraints.push_back(tpr);
      r = 1;
   }
   if (r == -1)
      std::cout << "WARNING:: Failure to add_extra_target_position_restraint for " << spec << std::endl;
   return r;
}

#ifdef HAVE_CXX11
// extra target position restraints refine like pull atom restraints
// but are stored differently (like other extra restraints)
int
molecule_class_info_t::add_extra_target_position_restraints(const std::vector<std::tuple<coot::atom_spec_t, const clipper::Coord_orth , float > > &etprs) {

   int r = -1;
   for (std::size_t i=0; i<etprs.size(); i++) {

      coot::atom_spec_t spec         = std::get<0>(etprs[i]); // copy
      const clipper::Coord_orth &pos = std::get<1>(etprs[i]);
      float  weight                  = std::get<2>(etprs[i]);
      mmdb::Atom *at = get_atom(spec);
      if (at) {
	 int atom_index = -1;
	 at->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
	 spec.int_user_data = atom_index;
	 coot::extra_restraints_t::extra_target_position_restraint_t tpr(spec, pos, weight);
	 extra_restraints.target_position_restraints.push_back(tpr);
	 r = 1;
      }
   }
   return r;
}
#endif

void
molecule_class_info_t::add_parallel_plane_restraint(coot::residue_spec_t spec_1,
						    coot::residue_spec_t spec_2) {

   std::vector<std::string> ap_1_names;
   std::vector<std::string> ap_2_names;
   std::string alt_conf_1; // only a nutter would want to add parallel plane restraints
   std::string alt_conf_2; // to a residue with an alt conf.

   mmdb::Residue *r_1 = get_residue(spec_1);
   mmdb::Residue *r_2 = get_residue(spec_2);

   if (r_1) {
      if (r_2) {

	 std::string rn_1 = r_1->GetResName();
	 std::string rn_2 = r_2->GetResName();

	 ap_1_names = nucelotide_residue_name_to_base_atom_names(rn_1);
	 ap_2_names = nucelotide_residue_name_to_base_atom_names(rn_2);

	 if (ap_1_names.empty()) ap_1_names = residue_name_to_plane_atom_names(rn_1);
	 if (ap_2_names.empty()) ap_2_names = residue_name_to_plane_atom_names(rn_2);

	 std::cout << "ap_2_names ";
	 for (auto i: ap_2_names)
	    std::cout << i << " ";
	 std::cout << "" << std::endl;

	 std::cout << "Adding parallel plane restraint " << spec_1 << " " << spec_2 << std::endl;
	 coot::parallel_planes_t pp(spec_1, spec_2, ap_1_names, ap_2_names,
				    alt_conf_1, alt_conf_2);

	 extra_restraints.parallel_plane_restraints.push_back(pp);

      } else {
	 std::cout << "INFO:: missing residue 2 " << spec_2 << std::endl;
      }
   } else {
	 std::cout << "INFO:: missing residue 1 " << spec_1 << std::endl;
   }

   update_extra_restraints_representation_parallel_planes();
}

std::vector<std::string>
molecule_class_info_t::nucelotide_residue_name_to_base_atom_names(const std::string &rn) const {

   std::vector<std::string> names;

   if (rn == "A" || rn == "G" || rn == "DA" || rn == "DG") {
      names.push_back("N9"); names.push_back("C8");
      names.push_back("N7"); names.push_back("C8");
      names.push_back("C5"); names.push_back("C6");
      names.push_back("C4"); names.push_back("C6");
      names.push_back("N3"); names.push_back("C2");
      names.push_back("N1");
   }

   if (rn == "C"  || rn == "T" || rn == "U" ||
       rn == "DC" || rn == "DT") {
      names.push_back("N1"); names.push_back("C2");
      names.push_back("C3"); names.push_back("C4");
      names.push_back("N3"); names.push_back("C5");
      names.push_back("O2");
   }

   return names;
}


std::vector<std::string>
molecule_class_info_t::residue_name_to_plane_atom_names(const std::string &rn) const {

   std::vector<std::string> names;

   if (rn == "PHE" || rn == "TYR") {
      names.push_back("CG");  names.push_back("CZ");
      names.push_back("CD1"); names.push_back("CD2");
      names.push_back("CE1"); names.push_back("CE2");
   }
   if (rn == "ARG") {
      names.push_back("CD");  names.push_back("NE");
      names.push_back("CZ");
      names.push_back("NH1"); names.push_back("NH2");
   }
   if (rn == "TRP") {
      names.push_back("CG");   names.push_back("CD1");
      names.push_back("NE1");  names.push_back("CE2");
      names.push_back("CD2");  names.push_back("CE3");
      names.push_back("CZ2");  names.push_back("CH2");
      names.push_back("CZ3");
   }
   return names;
}




void
molecule_class_info_t::delete_extra_restraints_for_residue(const coot::residue_spec_t &rs) {

   unsigned int pre_n = extra_restraints.bond_restraints.size();
   extra_restraints.delete_restraints_for_residue(rs);
   unsigned int post_n = extra_restraints.bond_restraints.size();
   if (post_n != pre_n)
      update_extra_restraints_representation();
}


// 20131014 unused currently.
bool
spec_pair_comparer(const std::pair<coot::atom_spec_t, coot::atom_spec_t> &p1,
		   const std::pair<coot::atom_spec_t, coot::atom_spec_t> &p2) {

//    if (p1.first < p2.first) {
//       return true;
//    } else {
//       if (p2.first < p1.first) {
// 	 return false;
//       } else {
// 	 return (p1.second < p2.second);
//       }
//    }

   if (p1.first < p2.first) {
      std::cout << "spec_pair_comparer A " << "[" << p1.first << " , " << p1.second << "]" << " < " "[" << p2.first << " , " << p2.second << "]" << std::endl;
      return true;
   } else {
      if (p2.first < p1.first) {
	 std::cout << "spec_pair_comparer B " << "[" << p2.first << " , " << p2.second << "]" " < " "[" << p1.first << " , " << p1.second << "]" << std::endl;
	 return false;
      } else {
	 bool v = (p1.second < p2.second);
	 if (v) {
	    std::cout << "spec_pair_comparer C " << "[" << p1.first << " , " << p1.second << "]" " < " "[" << p2.first << " , " << p2.second << "]" << std::endl;
	    return true;
	 } else {
	    std::cout << "spec_pair_comparer D " << "[" << p2.first << " , " << p2.second << "]" " < " "[" << p1.first << " , " << p1.second << "]" << std::endl;
	    return false;
	 }
      }
   }
}

void
molecule_class_info_t::delete_extra_restraints_worse_than(const double &n_sigma) {


   unsigned int pre_n = extra_restraints.bond_restraints.size();

   std::vector<coot::extra_restraints_t::extra_bond_restraint_t> ebrv_l;
   ebrv_l.reserve(extra_restraints.bond_restraints.size());

   std::map<coot::atom_spec_t, mmdb::Atom *> atom_map;
   std::map<coot::atom_spec_t, mmdb::Atom *>::const_iterator it_1;
   std::map<coot::atom_spec_t, mmdb::Atom *>::const_iterator it_2;
   // first fill the dist_map and fill the atom_map as you do so.
   for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {
      coot::extra_restraints_t::extra_bond_restraint_t &br = extra_restraints.bond_restraints[i];
      mmdb::Atom *at_1 = NULL;
      mmdb::Atom *at_2 = NULL;
      int ifast_index_1 = br.atom_1.int_user_data;
      int ifast_index_2 = br.atom_2.int_user_data;
      // std::cout << "fast indices: " << ifast_index_1 << " " << ifast_index_2 << std::endl;

      if (ifast_index_1 != -1) {
	 if (ifast_index_1 < atom_sel.n_selected_atoms) {
	    mmdb::Atom *at = atom_sel.atom_selection[ifast_index_1];
	    if (br.atom_1.matches_spec(at)) {
	       at_1 = at;
	    }
	 }
      }
      if (ifast_index_2 != -1) {
	 if (ifast_index_2 < atom_sel.n_selected_atoms) {
	    mmdb::Atom *at = atom_sel.atom_selection[ifast_index_2];
	    if (br.atom_2.matches_spec(at)) {
	       at_2 = at;
	    }
	 }
      }

      // If fast indexing doesn't do the job, let's try the map

      if (! at_1) {
	 it_1 = atom_map.find(br.atom_1);
	 if (it_1 == atom_map.end()) {
	    at_1 = get_atom(br.atom_1);
	    atom_map[br.atom_1] = at_1;
	 } else {
	    at_1 = it_1->second;
	 }
      }

      // and the same for 2:

      if (! at_1) {
	 it_2 = atom_map.find(br.atom_2);
	 if (it_2 == atom_map.end()) {
	    at_2 = get_atom(br.atom_2);
	    atom_map[br.atom_2] = at_2;
	 } else {
	    at_2 = it_2->second;
	 }
      }

      if (at_1 && at_2) {
	 double dx = at_1->x - at_2->x;
	 double dy = at_1->y - at_2->y;
	 double dz = at_1->z - at_2->z;
	 double dist_sq = dx*dx + dy*dy + dz*dz;
	 double dist = sqrt(dist_sq);

	 double this_n_sigma = fabs((br.bond_dist -dist)/br.esd);
	 if (this_n_sigma < n_sigma) {
	    ebrv_l.push_back(br);
	 }
      } else {
	 if (! at_1)
	    std::cout << "WARNING: missing atom_1 " << br.atom_1 << " when constructing dist_map" << std::endl;
	 if (! at_2)
	    std::cout << "WARNING: missing atom_2 " << br.atom_2 << " when constructing dist_map" << std::endl;
      }
   }
   extra_restraints.bond_restraints = ebrv_l;

   // remove_if and erase formulation.  Should work.  Crashes for some reason.
   //
   // extra_restraints.bond_restraints.erase(std::remove_if(extra_restraints.bond_restraints.begin(), extra_restraints.bond_restraints.end(), coot::extra_restraints_t::bond_eraser(dist_map, n_sigma)), extra_restraints.bond_restraints.end());

   unsigned int post_n = extra_restraints.bond_restraints.size();
   if (post_n != pre_n)
      update_extra_restraints_representation();
   std::cout << "INFO deleted : " << pre_n - post_n << " of " << pre_n << " extra bond restraints" << std::endl;
}

void
molecule_class_info_t::set_extra_restraints_prosmart_sigma_limits(double limit_low, double limit_high) {

   extra_restraints_representation.prosmart_restraint_display_limit_low  = limit_low;
   extra_restraints_representation.prosmart_restraint_display_limit_high = limit_high;

   // and redraw

   update_extra_restraints_representation();
}

void
molecule_class_info_t::generate_self_restraints(float local_dist_max,
						const coot::protein_geometry &geom) {

   // Find all the contacts in chain_id that are less than or equal to local_dist_max
   // that are not bonded or related by an angle.

   int selHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(selHnd, 0, "*",
			     mmdb::ANY_RES, "*", // start, insertion code
			     mmdb::ANY_RES, "*", // end, insertion code
			     "*", // residue name
			     "*",
			     "*", // elements
			     "*"); // alt locs

   generate_local_self_restraints(selHnd, local_dist_max, geom);
   atom_sel.mol->DeleteSelection(selHnd);

}


// make them yourself - easy as pie.
void
molecule_class_info_t::generate_local_self_restraints(float local_dist_max,
						      const std::string &chain_id,
						      const coot::protein_geometry &geom) {


   // Find all the contacts in chain_id that are less than or equal to local_dist_max
   // that are not bonded or related by an angle.

   int selHnd = atom_sel.mol->NewSelection(); // - check the deletion

   atom_sel.mol->SelectAtoms(selHnd, 0, chain_id.c_str(),
			     mmdb::ANY_RES, "*", // start, insertion code
			     mmdb::ANY_RES, "*", // end, insertion code
			     "*", // residue name
			     "*",
			     "*", // elements
			     "*"); // alt locs

   generate_local_self_restraints(selHnd, local_dist_max, geom);

   atom_sel.mol->DeleteSelection(selHnd);
}

void
molecule_class_info_t::generate_local_self_restraints(int selHnd, float local_dist_max,
						      const coot::protein_geometry &geom) {

   if (false)
      std::cout << "DEBUG:: here we are in mci::generate_local_self_restraints()! " << local_dist_max
                << std::endl;

   // clear what's already there - if anything
   extra_restraints.bond_restraints.clear();

   int nSelAtoms;
   mmdb::PPAtom SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, SelAtom, nSelAtoms);

   // std::cout << "here we are in mci::generate_local_self_restraints()! nSelAtoms " << nSelAtoms << std::endl;
   // bonded_neighbours in this case, means bonded or angle-related
   // bonded_neighbours["ALA"] -> all bond pairs and 1-3 angles
   std::map<std::string, std::vector<std::pair<std::string, std::string> > > bonded_neighbours;

   // now find contacts:
   //
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
	 my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   //
   atom_sel.mol->SeekContacts(SelAtom,   nSelAtoms,
			      SelAtom,   nSelAtoms,
			      0.1, local_dist_max, // min, max distances
			      0,        // seqDist 0 -> in same res also
			      pscontact, n_contacts,
			      0, &my_matt, i_contact_group);

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i=0; i<n_contacts; i++) {

            // 20221223-PE don't go both ways:
            if (pscontact[i].id1 > pscontact[i].id2) continue;

	    mmdb::Atom *at_1 = SelAtom[pscontact[i].id1];
	    mmdb::Atom *at_2 = SelAtom[pscontact[i].id2];
	    std::string ele_1 = at_1->element;
	    std::string ele_2 = at_2->element;
	    if (ele_1 != " H" && ele_2 != " H") {
	       bool ignore_this = false; // set for bonded and angled atoms
	       bool in_same_res = false;
	       if (at_1->residue == at_2->residue)
		  in_same_res = true;

	       if (in_same_res) {
		  std::string comp_id = at_1->residue->GetResName();
		  std::string at_name_1 = at_1->GetAtomName();
		  std::string at_name_2 = at_2->GetAtomName();

		  // This is slow
		  // if (geom.are_bonded_or_angle_related(comp_id, at_name_1, at_name_2))
		  // ignore_this = true;
		  //

		  std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator it;
		  it = bonded_neighbours.find(comp_id);
		  std::vector<std::pair<std::string, std::string> > bps;
		  if (it == bonded_neighbours.end()) {
		     bps = geom.get_bonded_and_1_3_angles(comp_id, imol_no);
		     bonded_neighbours[comp_id] = bps;
		  } else {
		     bps = it->second;
		  }

		  for (unsigned int ipr=0; ipr<bps.size(); ipr++) {
		     if (at_name_1 == bps[ipr].first) {
			if (at_name_2 == bps[ipr].second) {
			   ignore_this = true;
			   break;
			}
		     }
		     if (at_name_2 == bps[ipr].first) {
			if (at_name_1 == bps[ipr].second) {
			   ignore_this = true;
			   break;
			}
		     }
		  }
	       } else {
		  std::string at_name_1 = at_1->GetAtomName();
		  std::string at_name_2 = at_2->GetAtomName();
		  // by hand (hmmm)
		  // PDBv3 FIXME
		  if ((at_name_1 == " N  " && at_name_2 == " C  ") ||
		      (at_name_1 == " C  " && at_name_2 == " N  ") ||
		      (at_name_1 == " N  " && at_name_2 == " O  ") ||
		      (at_name_1 == " N  " && at_name_2 == " CA ") ||
		      (at_name_1 == " O  " && at_name_2 == " N  ") ||
		      (at_name_1 == " CA " && at_name_2 == " N  ")) {
		     ignore_this = true;
		  }
	       }

	       if (! ignore_this) {
		  clipper::Coord_orth p1 = coot::co(at_1);
		  clipper::Coord_orth p2 = coot::co(at_2);
		  double dist = sqrt((p1-p2).lengthsq());
		  double esd  = 0.05;
		  coot::atom_spec_t atom_spec_1(at_1);
		  coot::atom_spec_t atom_spec_2(at_2);
		  int idx_1 = -1;
		  int idx_2 = -1;
		  at_1->GetUDData(atom_sel.UDDAtomIndexHandle, idx_1);
		  at_2->GetUDData(atom_sel.UDDAtomIndexHandle, idx_2);
		  atom_spec_1.int_user_data = idx_1;
		  atom_spec_2.int_user_data = idx_2;
		  coot::extra_restraints_t::extra_geman_mcclure_restraint_t gmr(atom_spec_1, atom_spec_2,
										dist, esd);

		  // 20191120-PE self restraints are GM restraints, not bond
		  // restraints (previously bond restraints were GM only)
		  //
		  // extra_restraints.bond_restraints.push_back(br);

		  extra_restraints.geman_mcclure_restraints.push_back(gmr);

	       }
	    }
	 }
         delete [] pscontact;
      }
   }

   if (true) {
      std::cout << "----------- extra restraints ---------" << std::endl;
      std::cout << "     " << extra_restraints.bond_restraints.size() << std::endl;
      std::cout << "     " << extra_restraints.geman_mcclure_restraints.size() << std::endl;
   }

   // 20180510-PE surely I want to update the representation even if there are no
   //             extra bond restraints?
   //
   // if (extra_restraints.bond_restraints.size())
   // update_extra_restraints_representation();
   //
   update_extra_restraints_representation();

   // delete the selection from the calling function, not here.
   // atom_sel.mol->DeleteSelection(selHnd);
}

void
molecule_class_info_t::generate_local_self_restraints(float local_dist_max,
						      const std::vector<coot::residue_spec_t> &residue_specs,
						      const coot::protein_geometry &geom) {

   // Find all the contacts in chain_id that are less than or equal to local_dist_max
   // that are not bonded or related by an angle.

   int selHnd = coot::specs_to_atom_selection(residue_specs, atom_sel.mol, 0);
   if (selHnd >= 0) {
      generate_local_self_restraints(selHnd, local_dist_max, geom);
   }
   atom_sel.mol->DeleteSelection(selHnd);
}




void
molecule_class_info_t::clear_extra_restraints() {
   extra_restraints.clear();
   update_extra_restraints_representation();
}

// return an index of the new restraint
int
molecule_class_info_t::add_extra_torsion_restraint(coot::atom_spec_t atom_1,
						   coot::atom_spec_t atom_2,
						   coot::atom_spec_t atom_3,
						   coot::atom_spec_t atom_4,
						   double torsion_angle, double esd, int period) {

   mmdb::Atom *at_1 = get_atom(atom_1);
   mmdb::Atom *at_2 = get_atom(atom_2);
   mmdb::Atom *at_3 = get_atom(atom_3);
   mmdb::Atom *at_4 = get_atom(atom_4);
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


int
molecule_class_info_t::morph_fit_all(const clipper::Xmap<float> &xmap_in, float transformation_average_radius) {

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_neighb=2; // either side of central residue
   int n_chains = model_p->GetNumberOfChains();

   // the central residue and it's upstream and downstream neighbours (if it has them)
   std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > > moving_residues;

   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();

      for (int ires=0; ires<nres; ires++) { // residue-in-chain loop
	 std::vector<mmdb::Residue *> v; // up and downstream neighbours

	 for (int ifragres=-n_neighb; ifragres<=n_neighb; ifragres++) {
	    if (ifragres != 0) {
	       int idx = ires+ifragres;
	       if ((idx >=0) && (idx<nres)) {
		  mmdb::Residue *r = chain_p->GetResidue(idx);
		  if (r)
		     v.push_back(r);
	       }
	    }
	 }
	 std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > p(chain_p->GetResidue(ires), v);
	 moving_residues.push_back(p);
      }
   }
   return morph_fit_residues(moving_residues, xmap_in, transformation_average_radius);
}

int
molecule_class_info_t::morph_fit_chain(const std::string &chain_id,
				       const clipper::Xmap<float> &xmap_in, float transformation_average_radius) {

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_neighb=2; // either side of central residue
   int n_chains = model_p->GetNumberOfChains();

   // the central residue and it's upstream and downstream neighbours (if it has them)
   std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > > moving_residues;

   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (std::string(chain_p->GetChainID()) == chain_id) {
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) { // residue-in-chain loop
	    std::vector<mmdb::Residue *> v; // up and downstream neighbours

	    for (int ifragres=-n_neighb; ifragres<=n_neighb; ifragres++) {
	       if (ifragres != 0) {
		  int idx = ires+ifragres;
		  if ((idx >=0) && (idx<nres)) {
		     mmdb::Residue *r = chain_p->GetResidue(idx);
		     if (r)
			v.push_back(r);
		  }
	       }
	    }
	    std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > p(chain_p->GetResidue(ires), v);
	    moving_residues.push_back(p);
	 }
      }
   }
   return morph_fit_residues(moving_residues, xmap_in, transformation_average_radius);

}

int
molecule_class_info_t::morph_fit_residues(const std::vector<coot::residue_spec_t> &residue_specs,
					  const clipper::Xmap<float> &xmap_in, float transformation_average_radius) {

   // fill this from specs:
   std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > > moving_residues;
   for (unsigned int i=0; i<residue_specs.size(); i++) {
      mmdb::Residue *r = get_residue(residue_specs[i]);
      if (r) {
	 std::vector<mmdb::Residue *> env_residues =
	    coot::residues_near_residue(r, atom_sel.mol, transformation_average_radius);
	 std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > p(r, env_residues);
	 moving_residues.push_back(p);
      }
   }
   int r = 0;
   if (! moving_residues.empty())
      r = morph_fit_residues(moving_residues, xmap_in, transformation_average_radius);
   return r;
}

std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
molecule_class_info_t::peptide_C_N_pairs(const std::vector<mmdb::Residue *> &residues) const {

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v;
   std::vector<mmdb::Chain *> chains;
   for (unsigned int ires=0; ires<residues.size(); ires++) {
      mmdb::Residue *r = residues[ires];
      mmdb::Chain *this_chain = r->GetChain();
      if (std::find(chains.begin(), chains.end(), this_chain) == chains.end()) {
	 chains.push_back(this_chain);
      }
   }
   for (std::size_t ic=0; ic<chains.size(); ic++) {
      mmdb::Chain *chain_p = chains[ic];
      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > chain_pairs = coot::util::peptide_C_N_pairs(chain_p);
      v.insert(v.end(), chain_pairs.begin(), chain_pairs.end());
   }

   return v;
}


int
molecule_class_info_t::morph_fit_residues(std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Residue *> > > moving_residues,
					  const clipper::Xmap<float> &xmap_in,
					  float transformation_average_radius) {

   int success = 0;

   // construct minimol fragments
   bool ignore_pseudo_zeros = coot::util::is_EM_map(xmap_in);
   int n_bins = 40;
   mean_and_variance<float> mv = map_density_distribution(xmap_in, n_bins, false, ignore_pseudo_zeros);
   float map_rmsd = sqrt(mv.variance);

   // store the local origin too.
   std::map<mmdb::Residue *, morph_rtop_triple> rtop_map;

   for (unsigned int ires=0; ires<moving_residues.size(); ires++) {

      std::vector<mmdb::Residue *> fragment_residues;
      mmdb::Residue *residue_p = moving_residues[ires].first;
      std::pair<bool, clipper::Coord_orth> local_centre = residue_centre(residue_p);
      if (local_centre.first) {
	 std::cout << "\rINFO:: Getting RTops for " << coot::residue_spec_t(residue_p);
	 std::cout.flush();

	 fragment_residues.push_back(moving_residues[ires].first);
	 for (unsigned int ires_l=0; ires_l<moving_residues[ires].second.size(); ires_l++) {
	    mmdb::Residue *r = moving_residues[ires].second[ires_l];
	    fragment_residues.push_back(r);
	 }
	 coot::minimol::fragment f;
	 for (unsigned int ifr=0; ifr<fragment_residues.size(); ifr++) {
	    coot::minimol::residue fr(fragment_residues[ifr]->GetSeqNum());

	    bool add_all_residue = true; // change this if the residue is an amino acid
	    if (fragment_residues[ifr]->isAminoacid())
	       add_all_residue = false; // because we want only the main_chain

	    mmdb::PAtom *residue_atoms = 0;
	    int n_residue_atoms;
	    fragment_residues[ifr]->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {

	       if (! add_all_residue) {
		  if (coot::is_main_chain_p(residue_atoms[iat])) {
		     fr.addatom(coot::minimol::atom(residue_atoms[iat]));
		  }
	       } else {
		  // add all atoms of this residue
		  fr.addatom(coot::minimol::atom(residue_atoms[iat]));
	       }
	    }
	    f.addresidue(fr, false);
	 }
	 coot::minimol::molecule m(f);

	 coot::minimol::molecule m_copy = m; // for debugging

	 // returns the local rtop (relative to local centre) to move m into map.
	 std::pair<bool, clipper::RTop_orth> rtop = coot::get_rigid_body_fit_rtop(&m, local_centre.second, xmap_in, map_rmsd);
	 if (rtop.first) {
	    morph_rtop_triple rt(local_centre.second, rtop);
	    rtop_map[residue_p] = rt;

	    // debugging
	    //
	    if (false) {
	       coot::rigid_body_fit(&m_copy, xmap_in, map_rmsd);
	       std::string file_name = "morph-" + coot::util::int_to_string(ires);
	       file_name += ".pdb";
	       m_copy.write_file(file_name, 10);
	       std::cout << "    rtop for " << residue_p->GetSeqNum() << " " << residue_p->GetResName()
			 << " local centre " << local_centre.second.format() << " is " << std::endl;
	       std::cout << rt.rtop.format() << std::endl;
	    }
	 }
      }
      std::cout << std::endl; // for \r RTop specs
   }

   std::map<mmdb::Residue *, morph_rtop_triple>::const_iterator it;
   // std::cout << "rtop_map.size(): " << rtop_map.size() << std::endl;
   if (rtop_map.size()) {
      success = 1;
      make_backup(__FUNCTION__);

      // get the "key" residues, for use to find the peptide_C_N_pairs
      std::vector<mmdb::Residue *> residues_vec;
      for (unsigned int ires=0; ires<moving_residues.size(); ires++)
	 residues_vec.push_back(moving_residues[ires].first);
      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > c_n_pairs = peptide_C_N_pairs(residues_vec);

      std::map<mmdb::Residue *, morph_rtop_triple> simple_shifts;
      std::map<mmdb::Residue *, morph_rtop_triple> smooth_shifts;

      for (it=rtop_map.begin(); it!=rtop_map.end(); it++) {
	 mmdb::Residue *this_residue = it->first;
	 if (it->second.valid) {

	    // Morphing step is super-fast
	    // std::cout << "\rINFO:: Morphing " << coot::residue_spec_t(this_residue);
	    // std::cout.flush();

	    std::vector<std::pair<clipper::RTop_orth, float> > rtops;
	    // std::cout << "this residue:\n" << it->second.second.format() << std::endl;
	    rtops.push_back(std::pair<clipper::RTop_orth,float>(it->second.rtop, 1));
 	    std::vector<mmdb::Residue *> neighb_residues =
 	       coot::residues_near_residue(this_residue, atom_sel.mol, transformation_average_radius);

  	    for (unsigned int i_n_res=0; i_n_res<neighb_residues.size(); i_n_res++) {
	       std::map<mmdb::Residue *, morph_rtop_triple>::const_iterator it_for_neighb =
		  rtop_map.find(neighb_residues[i_n_res]);
	       if (it_for_neighb != rtop_map.end()) {
		  if (it_for_neighb->second.valid) {
		     float weight = 0.1;
		     float d_r = distance_between_residues(this_residue, neighb_residues[i_n_res]);
		     if (d_r > 0) {
			weight = 3.8/d_r;
			if (weight > 1.0)
			   weight = 1.0; // weight of central residue, shouldn't be more than that.
			// std::cout << "distance " << d_r << " weight " << weight << std::endl;
		     }
		     std::pair<clipper::RTop_orth, float> p(it_for_neighb->second.rtop, weight);
		     rtops.push_back(p);
		     if (0)
			std::cout << "adding rtop for "
				  << coot::residue_spec_t(neighb_residues[i_n_res]) << "\n"
				  << it_for_neighb->second.rtop.format() << std::endl;
		  }
	       }
	    }

	    // pre-local shifts and quaternion-based rtop averaging
	    // morph_residue_atoms_by_average_rtops(this_residue, rtops);

	    coot::util::quaternion q(0,0,0,0);
	    // clipper::RTop_orth rtop = q.centroid_rtop(rtops);
	    bool robust_filter = true;
	    clipper::RTop_orth rtop = q.centroid_rtop(rtops, robust_filter);

	    // debugging: save just to view them
	    simple_shifts[this_residue] = it->second;
	    smooth_shifts[this_residue] =
                  morph_rtop_triple(it->second.co, std::pair<bool, clipper::RTop_orth>(true, rtop));

	 } else {
	    std::cout << "no RTop for " << coot::residue_spec_t(it->first) << std::endl;
	 }
      }

      crunch_model_t crunch_model = morph_fit_crunch_analysis(smooth_shifts);
      morph_fit_uncrunch(&smooth_shifts, crunch_model); // alter smooth shifts

      for (it=smooth_shifts.begin(); it!=smooth_shifts.end(); it++) {

	 mmdb::Residue *this_residue = it->first;
	 if (it->second.valid) {
 	    translate_by_internal(-it->second.co,  it->first);
 	    transform_by_internal(it->second.rtop, it->first);
 	    translate_by_internal(it->second.co,   it->first);
	 }
      }

      std::cout << std::endl;
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();


      morph_show_shifts(simple_shifts, smooth_shifts);
      coot::util::standardize_peptide_C_N_distances(c_n_pairs);

   }
   return success;
}


crunch_model_t
molecule_class_info_t::morph_fit_crunch_analysis(const std::map<mmdb::Residue *, morph_rtop_triple> &residue_shifts) const {

   // Fill centres
   //
   std::vector<clipper::Coord_orth> centres;
   std::map<mmdb::Residue *, morph_rtop_triple>::const_iterator it;
   for (it=residue_shifts.begin(); it!=residue_shifts.end(); it++) {
      mmdb::Residue *residue_p = it->first;
      if (residue_p) {
	 if (it->second.valid) {
	    std::pair<bool, clipper::Coord_orth> rc = residue_centre(residue_p);
	    if (rc.first) {
	       centres.push_back(rc.second);
	    }
	 }
      }
   }

   // get centre_pos from centres
   //
   clipper::Coord_orth sum_pos(0,0,0);
   for (unsigned int i=0; i<centres.size(); i++)
      sum_pos += centres[i];
   double d = 1.0/double(centres.size());
   clipper::Coord_orth centre_pos(sum_pos[0]*d, sum_pos[1]*d, sum_pos[2]*d);

   // fill data
   //
   std::vector<std::pair<double, double> > data;
   for (it=residue_shifts.begin(); it!=residue_shifts.end(); it++) {
      mmdb::Residue *residue_p = it->first;
      if (residue_p) {
	 if (it->second.valid) {
	    std::pair<bool, clipper::Coord_orth> rc = residue_centre(residue_p);
	    if (rc.first) {
	       clipper::Coord_orth local_vec = rc.second-centre_pos;
	       clipper::Coord_orth local_vec_unit(local_vec.unit());
	       double d_c = sqrt(local_vec.clipper::Coord_orth::lengthsq());
	       clipper::Coord_orth trn(it->second.rtop.trn());
	       double dp = clipper::Coord_orth::dot(local_vec_unit, trn);
	       // data.push_back(std::pair<double, double> (d, sqrt(trn.lengthsq())));
	       data.push_back(std::pair<double, double> (d_c, dp));
	    }
	 }
      }
   }


   coot::least_squares_fit lsq(data);

   if (true) { // debug
      std::cout << "data: m " << lsq.m() <<  "  c: " << lsq.c() << std::endl;
      std::ofstream f("morph.tab");
      for (unsigned int idata=0; idata<data.size(); idata++)
	 f << idata << "   " << data[idata].first << " " << data[idata].second << "\n";
      std::ofstream fm("model.tab");
      for (unsigned int r=0; r<40; r++)
	 fm << r << " " << lsq.m() * r + lsq.c() <<  "\n";
   }

   crunch_model_t crunch_model(lsq, centre_pos);
   return crunch_model;

}

// alter shifts
void
molecule_class_info_t::morph_fit_uncrunch(std::map<mmdb::Residue *, morph_rtop_triple> *shifts,
					  crunch_model_t crunch_model) {

   std::map<mmdb::Residue *, morph_rtop_triple>::iterator it;
   for (it=shifts->begin(); it!=shifts->end(); it++) {
      mmdb::Residue *residue_p = it->first;
      if (residue_p) {
	 if (it->second.valid) {
	    std::pair<bool, clipper::Coord_orth> rc = residue_centre(residue_p);
	    if (rc.first) {
	       clipper::Coord_orth local_vec = rc.second-crunch_model.centre;
	       clipper::Coord_orth local_vec_unit(local_vec.unit());
	       double d_c = sqrt(local_vec.clipper::Coord_orth::lengthsq()); // distance from centre
	       double uncrunch_ampl = crunch_model.lsq.m() * d_c + crunch_model.lsq.c();
	       clipper::Coord_orth uncrunch_vec(- uncrunch_ampl * local_vec_unit);
	       clipper::Coord_orth new_vec(it->second.rtop.trn() + uncrunch_vec);
	       clipper::RTop_orth new_rtop(it->second.rtop.rot(), new_vec);
	       it->second.rtop = new_rtop;
	    }
	 }
      }
   }
}




// I fail to make a function that does a good "average" of RTops,
// so do it long-hand by generating sets of coordinates by applying
// each rtop to each atom - weights are transfered in the second part of the pair.
//
// This doesn't do backups or unsaved changes marking of course.
void
molecule_class_info_t::morph_residue_atoms_by_average_rtops(mmdb::Residue *residue_p,
							    const std::vector<std::pair<clipper::RTop_orth, float> > &rtops) {

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   if (rtops.size()) {
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 clipper::Coord_orth pt(at->x, at->y, at->z);
	 clipper::Coord_orth sum_transformed_pts(0,0,0);
	 double sum_weights = 0.0;
	 for (unsigned int i_rtop=0; i_rtop<rtops.size(); i_rtop++) {
	    clipper::Coord_orth t_pt = pt.transform(rtops[i_rtop].first);
	    double weight = rtops[i_rtop].second;
	    sum_weights += weight;
	    sum_transformed_pts += t_pt*weight;
	 }
	 if (sum_weights > 0.0) {
	    double inv_weight = 1.0/sum_weights;
	    clipper::Coord_orth new_pt(sum_transformed_pts.x() * inv_weight,
				       sum_transformed_pts.y() * inv_weight,
				       sum_transformed_pts.z() * inv_weight);
	    at->x = new_pt.x();
	    at->y = new_pt.y();
	    at->z = new_pt.z();
	 }
      }
   }
}



void
molecule_class_info_t::morph_show_shifts(const std::map<mmdb::Residue *, morph_rtop_triple> &simple_shifts,
					 const std::map<mmdb::Residue *, morph_rtop_triple> &smooth_shifts) const {

   // write a file

   std::map<mmdb::Residue *, morph_rtop_triple>::const_iterator it;
   std::ofstream f("morph-shifts.scm");

   std::string ss;
   ss = "(define simple-shifts (new-generic-object-number \"simple-shifts\"))";
   f << ss << "\n";
   ss = "(define smooth-shifts (new-generic-object-number \"smooth-shifts\"))";
   f << ss << "\n";
   ss = "(set-display-generic-object simple-shifts 1)";
   f << ss << "\n";
   ss = "(set-display-generic-object smooth-shifts 1)";
   f << ss << "\n";

   for (it=simple_shifts.begin(); it!=simple_shifts.end(); it++) {
      mmdb::Residue *r = it->first;
      std::pair<bool, clipper::Coord_orth> rc = residue_centre(r);
      mmdb::Atom *C_alpha = r->GetAtom(" CA ");
      if (! C_alpha)
	 C_alpha = r->GetAtom(" P  ");
      if (C_alpha) {
	 clipper::Coord_orth ca_pos(C_alpha->x, C_alpha->y, C_alpha->z);
	 if (rc.first) {
	    std::string s;
	    std::string line_colour = "yellow";
	    std::string ball_colour = "yellow";

	    clipper::RTop_orth rtop_for_centre(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), it->second.co);
	    clipper::RTop_orth rtop_for_centre_i(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), -it->second.co);
	    clipper::Coord_orth tp_1   = ca_pos.transform(rtop_for_centre_i);
	    clipper::Coord_orth tp_2   =   tp_1.transform(it->second.rtop);
	    clipper::Coord_orth to_pos =   tp_2.transform(rtop_for_centre);

	    s += "(to-generic-object-add-line  simple-shifts ";
	    s += "\"";
	    s += line_colour;
	    s += "\"";
	    s += "  2 ";
	    s += coot::util::float_to_string(ca_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.z());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	    s = "";
	    s += "(to-generic-object-add-point simple-shifts ";
	    s += "\"";
	    s += ball_colour;
	    s += "\"";
	    s += " 12                   ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	 }
      }
   }
   for (it=smooth_shifts.begin(); it!=smooth_shifts.end(); it++) {
      mmdb::Residue *r = it->first;
      std::pair<bool, clipper::Coord_orth> rc = residue_centre(r);
      mmdb::Atom *C_alpha = r->GetAtom(" CA ");
      if (! C_alpha)
	 C_alpha = r->GetAtom(" P  ");
      if (C_alpha) {
	 clipper::Coord_orth ca_pos(C_alpha->x, C_alpha->y, C_alpha->z);
	 if (rc.first) {
	    std::string s;
	    std::string line_colour = "red";
	    std::string ball_colour = "red";

	    clipper::RTop_orth rtop_for_centre(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), it->second.co);
	    clipper::RTop_orth rtop_for_centre_i(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), -it->second.co);
	    clipper::Coord_orth tp_1   = ca_pos.transform(rtop_for_centre_i);
	    clipper::Coord_orth tp_2   =   tp_1.transform(it->second.rtop);
	    clipper::Coord_orth to_pos =   tp_2.transform(rtop_for_centre);

	    s += "(to-generic-object-add-line  smooth-shifts ";
	    s += "\"";
	    s += line_colour;
	    s += "\"";
	    s += " 2 ";
	    s += coot::util::float_to_string(ca_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.z());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	    s = "";
	    s += "(to-generic-object-add-point smooth-shifts ";
	    s += "\"";
	    s += ball_colour;
	    s += "\"";
	    s += " 14                 ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	 }
      }
   }
   f.close();
}

int
molecule_class_info_t::morph_fit_by_secondary_structure_elements(const std::string &chain_id,
								 const clipper::Xmap<float> &xmap_in,
								 float map_rmsd) {

   int status = 0;
   float local_radius = 16;

   int imodel = 1;
   bool simple_move = false;

   if (atom_sel.mol) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imodel);
      if (model_p) {

	 make_backup(__FUNCTION__);
	 bool model_changed = true;

	 int nhelix = model_p->GetNumberOfHelices();
	 int nsheet = model_p->GetNumberOfSheets();

	 if (nhelix == 0 && nsheet == 0)
	    add_secondary_structure_header_records(true);
	 mmdb::Helix * helix_p;
	 mmdb::Sheet * sheet_p;
	 mmdb::Strand * strand_p;

	 mmdb::Chain *chain_p = model_p->GetChain(chain_id.c_str());
	 if (chain_p) {

	    std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > c_n_pairs =
	       coot::util::peptide_C_N_pairs(chain_p);

	    // store the RTops for some residues (we also need the
	    // local around which the rtop_orth is performed)
	    //
	    std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> > rtop_map;

	    std::cout << "INFO:: " << nhelix << " helices" << std::endl;
	    std::cout << "INFO:: " << nsheet << " sheets"  << std::endl;

	    for (int ih=1; ih<=nhelix; ih++) {
	       helix_p = model_p->GetHelix(ih);
	       if (helix_p) {

		  std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> > rtops_fragment =
		     morph_fit_by_secondary_structure_fragment(chain_p, chain_id,
							       helix_p->initSeqNum, helix_p->endSeqNum,
							       xmap_in, map_rmsd,
							       simple_move);
		  // add rtops_fragment bits to overall rtops_map;
		  std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> >::const_iterator it;
		  for (it=rtops_fragment.begin(); it!=rtops_fragment.end(); it++)
		     rtop_map[it->first] = it->second;

	       } else {
		  std::cout << "ERROR: no helix!?" << std::endl;
	       }
	    }

	    for (int is=1; is<=nsheet; is++) {
	       sheet_p = model_p->GetSheet(is);

	       int nstrand = sheet_p->nStrands;
	       for (int istrand=0; istrand<nstrand; istrand++) {
		  strand_p = sheet_p->strand[istrand];
		  if (strand_p) {
		     std::cout << "---- handle strand ------ id: " << strand_p->sheetID << " # "
			       << strand_p->strandNo << " " << strand_p->initChainID << " "
			       << strand_p->initSeqNum << " "
			       << strand_p->endChainID << " "
			       << strand_p->endSeqNum
			       << std::endl;

		     if (std::string(strand_p->initChainID) == chain_id) {
			if (std::string(strand_p->endChainID) == chain_id) {
			   std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> > rtops_fragment =
			      morph_fit_by_secondary_structure_fragment(chain_p, chain_id,
									strand_p->initSeqNum, strand_p->endSeqNum,
									xmap_in, map_rmsd,
									simple_move);

			   // add rtops_fragment bits to overall rtops_map;
			   std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> >::const_iterator it;
			   for (it=rtops_fragment.begin(); it!=rtops_fragment.end(); it++)
			      rtop_map[it->first] = it->second;
			}
		     }
		  }
	       }
	    }

	    // OK, so now some residues (those in SSE) have rtops
	    //
	    // Now run over the residues and atoms of the chain and apply
	    // the local weighted average of the RTops to the coordinates
	    //
	    // The residues in a given SSE all have the same RTop_orth.
	    //
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::Residue *residue_p;
	    mmdb::Atom *at;
	    std::map<mmdb::Residue *, clipper::Coord_orth> residue_centres;
	    for (int ires=0; ires<nres; ires++) {

	       residue_p = chain_p->GetResidue(ires);
	       std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> >::const_iterator it_ss =
		  rtop_map.find(residue_p);
	       if (it_ss != rtop_map.end()) {

		  // OK this was a residue in a SSE.  We know how to move these atoms (i.e. use their own
		  // RTop, not morphing). This block is what Israel Sanchez-Fernandez wanted.
		  //
		  mmdb::PPAtom residue_atoms = 0;
		  int n_residue_atoms;
		  residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
		  clipper::Coord_orth centre = it_ss->second.first;
		  for (int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at_l = residue_atoms[iat];
		     clipper::Coord_orth pt_1 = coot::co(at_l);
		     clipper::Coord_orth pt_2 = pt_1 - centre;
		     clipper::Coord_orth pt_3 = pt_2.transform(it_ss->second.second);
		     clipper::Coord_orth pt_4 = pt_3 + centre;
		     coot::update_position(at_l, pt_4);
		  }

	       } else {

		  // ------------------------------------------------------------------
		  //    morph:  move the atoms in these residues based on the RTops of
		  //            the residues in their environments.
		  // ------------------------------------------------------------------

		  int n_atoms = residue_p->GetNumberOfAtoms();

		  // get the centre of this residue (residue_p) from
		  // the centre cache (if not in the cache, add it to
		  // the cache).
		  clipper::Coord_orth this_residue_centre(0,0,0);
		  if (residue_centres.find(residue_p) == residue_centres.end()) {
		     std::pair<bool, clipper::Coord_orth> pp = coot::util::get_residue_centre(residue_p);
		     if (pp.first) {
			this_residue_centre = pp.second;
			residue_centres[residue_p] = pp.second;
		     }
		  } else {
		     this_residue_centre = residue_centres[residue_p];
		  }

		  // get the centre of the environment residues from
		  // the cache (if not in the cache, add them to the
		  // cache).
		  std::vector<mmdb::Residue *> env_residues =
		     coot::residues_near_residue(residue_p, atom_sel.mol, local_radius);
		  for (unsigned int jres=0; jres<env_residues.size(); jres++) {
		     if (residue_centres.find(env_residues[jres]) == residue_centres.end()) {
			std::pair<bool, clipper::Coord_orth> pp =
			   coot::util::get_residue_centre(env_residues[jres]);
			if (pp.first)
			   residue_centres[env_residues[jres]] = pp.second;
		     }
		  }

		  for (int iat=0; iat<n_atoms; iat++) {
		     at = residue_p->GetAtom(iat);
		     clipper::Coord_orth pt_atom = coot::co(at);
		     // std::tuple<clipper::RTop_orth, clipper::Coord_orth, float>
		     std::vector<atom_morph_info_t> rtops_for_atom;
		     for (unsigned int ier=0; ier<env_residues.size(); ier++) {
			std::map<mmdb::Residue *, clipper::Coord_orth>::const_iterator it;
			it = residue_centres.find(env_residues[ier]);
			if (it != residue_centres.end()) {
			   const clipper::Coord_orth &pt_e_r = it->second;
			   double d_sqrd = (pt_e_r - pt_atom).lengthsq();
			   if (d_sqrd < 1.0) d_sqrd = 1.0;
			   double d = sqrt(d_sqrd);
			   double w = 1.0/d;
			   std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> >::const_iterator it_rtop =
			      rtop_map.find(env_residues[ier]);
			   if (it_rtop != rtop_map.end()) {
			      atom_morph_info_t t(it_rtop->second.second, pt_e_r, w);
			      rtops_for_atom.push_back(t);
			   }
			}
		     }

		     if (rtops_for_atom.size()) {

			std::vector<std::pair<clipper::RTop_orth,float> > rtop_pairs_for_atom(rtops_for_atom.size());
			for (unsigned int i=0; i<rtops_for_atom.size(); i++)
			   rtop_pairs_for_atom[i] = rtops_for_atom[i].rtop_and_weight();

			coot::util::quaternion q(0,0,0,0);
			bool robust_filter = true;
			clipper::RTop_orth rtop_for_atom = q.centroid_rtop(rtop_pairs_for_atom, robust_filter);
			clipper::Coord_orth p_1 = coot::co(at);
			clipper::Coord_orth p_2 = p_1 - this_residue_centre;
			clipper::Coord_orth p_3 = p_2.transform(rtop_for_atom);
			clipper::Coord_orth p_4 = p_3 + this_residue_centre;
			coot::update_position(at, p_4);
		     }
		  }
	       }
	    }

	    if (model_changed) {
	       coot::util::standardize_peptide_C_N_distances(c_n_pairs); // move the C and Ns closer if needed
	       atom_sel.mol->FinishStructEdit();
	       have_unsaved_changes_flag = 1;
	       make_bonds_type_checked();
	    }
	 }
      }
   }
   return status;
}

std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> >
molecule_class_info_t::morph_fit_by_secondary_structure_fragment(mmdb::Chain *chain_p,
								 const std::string &chain_id,
								 int initSeqNum,
								 int endSeqNum,
								 const clipper::Xmap<float> &xmap_in,
								 float map_rmsd,
								 bool simple_move) {


   std::map<mmdb::Residue *, std::pair<clipper::Coord_orth, clipper::RTop_orth> > rtop_map;

   coot::minimol::fragment f(chain_id);
   std::vector<mmdb::Residue *> added_residues;
   for (int res_no=initSeqNum; res_no<=endSeqNum; res_no++) {
      mmdb::Residue *residue_p = chain_p->GetResidue(res_no, "");
      if (residue_p) {
	 f.addresidue(coot::minimol::residue(residue_p), false);
	 added_residues.push_back(residue_p);
      } else {
	 std::cout << "Null residue for " << chain_id << " " << res_no << std::endl;
      }
   }

   if (!added_residues.size()) {
      std::cout << "no added residues for helix " << chain_id << " "
		<< initSeqNum << " " << endSeqNum << std::endl;
   } else {
      coot::minimol::molecule m(f);
      std::pair<bool, clipper::Coord_orth> sse_centre = coot::centre_of_residues(added_residues);
      if (sse_centre.first) {

	 // returns the local rtop (relative to local centre) to move m into map.
	 //
	 std::pair<bool, clipper::RTop_orth> rtop =
	    coot::get_rigid_body_fit_rtop(&m, sse_centre.second, xmap_in, map_rmsd);

	 if (rtop.first) {
	    if (0)
	       std::cout << "Got and RTop for SSE " << chain_id << " "
			 << initSeqNum << " -- " << endSeqNum
			 << std::endl;
	    for (unsigned int ires=0; ires<added_residues.size(); ires++) {
	       std::pair<clipper::Coord_orth, clipper::RTop_orth> p(sse_centre.second, rtop.second);
	       rtop_map[added_residues[ires]] = p;
	    }

	    if (simple_move) {
	       // simple move the coordinates
	       m.transform(rtop.second);
	       clipper::Mat33<double> mat(1,0,0,0,1,0,0,0,1);
	       clipper::RTop_orth rtop_synth(mat, sse_centre.second);
	       m.transform(rtop_synth);
	       atom_selection_container_t asc = make_asc(m.pcmmdbmanager());
	       replace_fragment(asc);
	    }
	 }
      }
   }
   return rtop_map;
}




// maybe extra parameters are needed here (e.g. for colouring later, perhaps).
void
coot::extra_restraints_representation_t::add_parallel_plane(const lsq_plane_info_t &pi_1,
							    const lsq_plane_info_t &pi_2) {

      double offset_d = 0.4;
      clipper::Coord_orth n_1 = pi_1.normal();
      clipper::Coord_orth n_2 = pi_2.normal();
      clipper::Coord_orth offset_pt_1_a = pi_1.centre() + offset_d * n_1;
      clipper::Coord_orth offset_pt_1_b = pi_1.centre() - offset_d * n_1;
      clipper::Coord_orth offset_pt_2_a = pi_2.centre() + offset_d * n_2;
      clipper::Coord_orth offset_pt_2_b = pi_2.centre() - offset_d * n_2;

      clipper::Coord_orth op_1 = offset_pt_1_a;
      if (clipper::Coord_orth(offset_pt_1_b - pi_2.centre()).lengthsq() <
	  clipper::Coord_orth(offset_pt_1_a - pi_2.centre()).lengthsq())
	 op_1 = offset_pt_1_b;
      clipper::Coord_orth op_2 = offset_pt_2_a;
      if (clipper::Coord_orth(offset_pt_2_b - pi_1.centre()).lengthsq() <
	  clipper::Coord_orth(offset_pt_2_a - pi_1.centre()).lengthsq())
	 op_2 = offset_pt_2_b;

      double f_top_1 = -(pi_2.a()*op_1.x() + pi_2.b()*op_1.y() + pi_2.c()*op_1.z() - pi_2.d());
      double f_top_2 = -(pi_1.a()*op_2.x() + pi_1.b()*op_2.y() + pi_1.c()*op_2.z() - pi_1.d());
      double f_bot   = pi_2.a() * pi_1.a() + pi_2.b() * pi_1.b() + pi_2.c() * pi_1.c();
      double s_1 = f_top_1/f_bot;
      double s_2 = f_top_2/f_bot;

      // projected plane points
      clipper::Coord_orth ppp_1 = op_1 + s_1 * pi_1.normal();
      clipper::Coord_orth ppp_2 = op_2 + s_2 * pi_2.normal();

      extra_parallel_planes_restraints_representation_t eppr_1(op_1, ppp_1, pi_1.normal(), 1.3, 0.2);
      extra_parallel_planes_restraints_representation_t eppr_2(op_2, ppp_2, pi_2.normal(), 1.3, 0.2);
      parallel_planes.push_back(eppr_1);
      parallel_planes.push_back(eppr_2);

}


#include "coot-utils/coot_shiftfield.h"

void
molecule_class_info_t::shiftfield_b_factor_refinement(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                                      const clipper::HKL_data<clipper::data32::Flag> &free) {
   if (atom_sel.mol) {
      make_backup(__FUNCTION__);
      int n_cycles = 3;
      coot::shift_field_b_factor_refinement(fobs, free, atom_sel.mol, n_cycles);
   }

}

void
molecule_class_info_t::shiftfield_xyz_factor_refinement(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                                        const clipper::HKL_data<clipper::data32::Flag> &free) {

   if (atom_sel.mol) {
      make_backup(__FUNCTION__);
      coot::shift_field_xyz_refinement(fobs, free, atom_sel.mol, 2.0);
   }
}
