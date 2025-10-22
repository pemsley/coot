/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013, 2015, 2016 by Medical Research Council
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

#include <optional>

#include "compat/coot-sysdep.h"

#include "coot-utils/bonded-pairs.hh"
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "simple-restraint.hh"

#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

// if they were not passed in the constructor.
void
coot::restraints_container_t::fill_links(mmdb::Manager *mol) {

   // fill std::vector<mmdb::Link> links

   links.clear(); // hmmm!

   if (mol) {
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
	 unsigned int n_links = model_p->GetNumberOfLinks();
	 for (unsigned int i=1; i<=n_links; i++) {
	    mmdb::Link *ref_link = model_p->GetLink(i);
	    if (ref_link) {
	       mmdb::Link l(*ref_link);
	       links.push_back(l);
	    }
	 }
      }
   }
   if (false)
      std::cout << "INFO:: refinement transfered " << links.size() << " links" << std::endl;
}


// Need to add test that residues are linked with trans
//
std::vector<coot::rama_triple_t>
coot::restraints_container_t::make_rama_triples(int SelResHnd,
						const coot::protein_geometry &geom) const {
   std::vector<coot::rama_triple_t> v;
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);

   for (int i=0; i<(nSelResidues-2); i++) {
      if (SelResidue[i] && SelResidue[i+1] && SelResidue[i+2]) {
	 std::string link_type = "TRANS";
	 rama_triple_t t(SelResidue[i], SelResidue[i+1], SelResidue[i+2], link_type);
	 v.push_back(t);
      }
   }
   return v;
}


coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_by_linear(int SelResHnd,
							const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t c;
   mmdb::PPResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {

      std::string link_type("TRANS");
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage

      for (int i=0; i<(nSelResidues-1); i++) {
	 if (SelResidue[i] && SelResidue[i+1]) {


	    // There are a couple of ways we allow links.  First, that
	    // the residue numbers are consecutive.
	    //
	    // If there is a gap in residue numbering, the C and N
	    // have to be within 3A.

	    // Are these residues neighbours?  We can add some sort of
	    // ins code test here in the future.
	    if ( abs(SelResidue[i]->GetSeqNum() - SelResidue[i+1]->GetSeqNum()) <= 1) {
	       link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       // std::cout << "DEBUG ------------ in bonded_residues_by_linear() link_type is :"
	       // << link_type << ":" << std::endl;
	       if (link_type != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					whole_first_residue_is_fixed,
					whole_second_residue_is_fixed, link_type);
		  c.try_add(p);
	       }
	    }

	    // distance check this one.... it could be the opposite of
	    // an insertion code, or simply a gap - and we don't want
	    // to make a bond for a gap.
	    //
	    if (abs(SelResidue[i]->index - SelResidue[i+1]->index) <= 1) {
	       // link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
               std::cout << "####################### find_link_type_compli() called from bonded_residues_by_linear()"
                         << std::endl;
	       std::pair<std::string, bool> link_info =
		  find_link_type_complicado(SelResidue[i], SelResidue[i+1], geom);
	       if (false)
		  std::cout << "DEBUG:: ---------- in bonded_residues_by_linear() link_info is :"
			    << link_info.first << " " << link_info.second << ":" << std::endl;

	       if (link_info.first != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  if (link_info.second == 0) {
		     coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_info.first);
		     c.try_add(p);
		  } else {
		     // order switch
		     coot::bonded_pair_t p(SelResidue[i+1], SelResidue[i],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_info.first);
		     c.try_add(p);
		  }
	       }
	    }
	 }
      }
   }
   return c;
}


coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_from_res_vec(const coot::protein_geometry &geom) const {

   bool debug = false; // Are your residues in the same chain?  If not, filter() will not bond them.
                       // 20221120-PE Hmm... that may not be true any more.

   coot::bonded_pair_container_t bpc;

   if (verbose_geometry_reporting == VERBOSE)
      debug = true;

   int nres = residues_vec.size();

   if (debug) {
      std::cout << "debug:: ------------------- bonded_residues_from_res_vec() residues_vec.size() "
		<< residues_vec.size() << std::endl;
      for (unsigned int i=0; i<residues_vec.size(); i++) {
	 std::cout << "   fixed: " << residues_vec[i].first << " spec: "
		   << residue_spec_t(residues_vec[i].second) << std::endl;
      }
      for (unsigned int ii=0; ii<residues_vec.size(); ii++) {
	 mmdb::Residue *res_f = residues_vec[ii].second;
	 for (unsigned int jj=ii+1; jj<residues_vec.size(); jj++) {
	    mmdb::Residue *res_s = residues_vec[jj].second;

	    std::cout << "debug:: ------------ test here with res_f and res_s "
		      << residue_spec_t(res_f) << " " << residue_spec_t(res_s) << std::endl;

	 }
      }
   }

   for (unsigned int ii=0; ii<residues_vec.size(); ii++) {
      mmdb::Residue *res_f = residues_vec[ii].second;
      for (unsigned int jj=ii+1; jj<residues_vec.size(); jj++) {
	 mmdb::Residue *res_s = residues_vec[jj].second;

	 if (debug)
	    std::cout << "debug:: ----- in bonded_resdues_from_res_vec " << residue_spec_t(res_f) << " "
		      << residue_spec_t(res_s) << "\n";

	 if (res_f == res_s) continue;

         // std::cout << "####################### find_link_type_compli() called from bonded_resdues_from_res_vec()"
         // << std::endl;
	 // 20180911 I now no longer want to evaluate closest approach here.
	 //
	 // std::pair<bool, float> d = closest_approach(res_f, res_s);
	 // Linking should be resolved by find_link_type_complicado(), not
	 // here by distance between residues.

         // Return the link type and a residue order switch flag.
         // Return link_type as "" if not found.
	 std::pair<std::string, bool> l = find_link_type_complicado(res_f, res_s, geom);

	 // too verbose?
	 if (false)
	    std::cout << "   INFO:: find_link_type_complicado() for: "
		      << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
		      << " returns link_type -> \"" << l.first << "\""
                      << " order-switch-flag: " << l.second << std::endl;

	 std::string link_type = l.first;
	 if (!link_type.empty()) {

	    bool whole_first_residue_is_fixed = 0;
	    bool whole_second_residue_is_fixed = 0;
	    bool order_switch_flag = l.second;

	    if (!order_switch_flag) {
	       coot::bonded_pair_t p(res_f, res_s,
				     whole_first_residue_is_fixed,
				     whole_second_residue_is_fixed, link_type);
	       bool previously_added_flag = bpc.try_add(p);
	    } else {
	       coot::bonded_pair_t p(res_s, res_f,
				     whole_first_residue_is_fixed,
				     whole_second_residue_is_fixed,
				     link_type);
	       bool previously_added_flag = bpc.try_add(p);
	       // std::cout << "              previously_added_flag " << previously_added_flag
	       // << std::endl;
	    }

	    // if the link type is a straight-forward TRANS link of 2 residue next to each other
	    // in residue numbers and serial numbers, then we don't need to find any other type
	    // of link for this residue (so break out of the inner for-loop).
	    bool was_straight_forward_trans_link = false;
	    int resno_1 = res_f->GetSeqNum();
	    int resno_2 = res_s->GetSeqNum();
	    int ser_num_1 = res_f->index;
	    int ser_num_2 = res_s->index;
	    if (resno_2 == (resno_1 + 1)) {
	       if (ser_num_2 == (ser_num_1 + 1)) {
		  std::string rn_1 = res_f->GetResName();
		  if (rn_1 != "ASN" && rn_1 != "CYS" && rn_1 != "SER" && rn_1 != "TYR") {
		     if (link_type == "TRANS" || link_type == "PTRANS")
			was_straight_forward_trans_link = true;
		  }
	       }
	    }
	    if (was_straight_forward_trans_link) {
	       // std::cout << "------------ was straight_forward TRANS link! - breaking"  << std::endl;
	       break;
	    }

	 } else {
	    if (debug)
	       std::cout << "DEBUG:: find_link_type_complicado() blank result: "
			 << "link_type find_link_type_complicado() for "
			 << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
			 << " returns \"" << l.first << "\" " << l.second << std::endl;
	 }
      }
   }

   bpc.filter(); // removes 1-3 bond items and if 1-2 and 1-3 bonds exist

   // std::cout << "---------------- done bonded_residues_from_res_vec()" << std::endl;

   return bpc;
}


// Add a trans bond linkage
//
// Residue 1 of the link is the first atom of the link
int
coot::restraints_container_t::add_link_bond(std::string link_type,
					    mmdb::PResidue first, mmdb::PResidue second,
					    short int is_fixed_first,
					    short int is_fixed_second,
					    const coot::protein_geometry &geom) {

   bool debug = false;

   mmdb::PPAtom first_sel;
   mmdb::PPAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   std::vector<bool> fixed_atom_flags(2);  // 2 atoms in this restraint.

   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);
   short int found_link_type = 0;
   int index1, index2;

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   if (debug) {
      std::cout << "INFO:: geom.link_size() is " << geom.link_size() << std::endl;
      std::cout << "first residue:\n";
      for (int i=0; i<n_first_res_atoms; i++)
	 std::cout << "    " << first_sel[i]->name  << " " << first_sel[i]->GetSeqNum() << "\n";
      std::cout << "second residue:\n";
      for (int i=0; i<n_second_res_atoms; i++)
	 std::cout << "    " << second_sel[i]->name  << " " << second_sel[i]->GetSeqNum() << "\n";
   }

   int nbond = 0;
   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically "TRANS"
	 found_link_type = 1;
	 for (unsigned int j=0; j<geom.link(i).link_bond_restraint.size(); j++) {
	    if (geom.link(i).link_bond_restraint[j].atom_1_comp_id == 1 && 
		geom.link(i).link_bond_restraint[j].atom_2_comp_id == 2) {
	       // as expected.
	    } else {
	       std::cout << "PROGRAMMER ERROR (shortsighted fool)" << std::endl;
	       std::cout << "bad things will now happen..." << std::endl; 
	    }
	    for (int ifat=0; ifat<n_first_res_atoms; ifat++) { 
	       std::string pdb_atom_name_1(first_sel[ifat]->name);

	       if (pdb_atom_name_1 == geom.link(i).link_bond_restraint[j].atom_id_1_4c()) {
		  for (int isat=0; isat<n_second_res_atoms; isat++) { 
		     std::string pdb_atom_name_2(second_sel[isat]->name);

		     if (pdb_atom_name_2 == geom.link(i).link_bond_restraint[j].atom_id_2_4c()) {

			if (debug)
			   std::cout << "DEBUG::  adding " << link_type << " bond for "
				     << first->seqNum
				     << " -> " << second->seqNum << " atoms "
				     << first_sel [ifat]->GetAtomName() << " to "
				     << second_sel[isat]->GetAtomName()
				     << std::endl;

			// Now, do the alt confs match?
			//
			std::string alt_conf_1 =  first_sel[ifat]->altLoc;
			std::string alt_conf_2 = second_sel[isat]->altLoc;
			if ((alt_conf_1 == alt_conf_2) || (alt_conf_1 == "") || (alt_conf_2 == "")) {

                           // 20230110-PE are you here again? Check that udd_atom_index_handle
                           // is set (by calling init_shared_post()).

			   first_sel [ifat]->GetUDData(udd_atom_index_handle, index1);
			   second_sel[isat]->GetUDData(udd_atom_index_handle, index2);

			   // set the UDD flag for this residue being bonded/angle with 
			   // the other

			   bonded_atom_indices[index1].insert(index2);
			   bonded_atom_indices[index2].insert(index1);

			   fixed_atom_flags[0] = is_fixed_first;
			   fixed_atom_flags[1] = is_fixed_second;

			   std::vector<bool> other_fixed_flags = make_fixed_flags(index1, index2);
			   for (int ii=0; ii<2; ii++)
			      if (other_fixed_flags[ii])
				 fixed_atom_flags[ii] = 1;

			   add(BOND_RESTRAINT, index1, index2,
			       fixed_atom_flags,
			       geom.link(i).link_bond_restraint[j].dist(),
			       geom.link(i).link_bond_restraint[j].esd(),
			       1.2); // junk value
			   nbond++;
			}
		     }
		  }
	       }
	    }
	 }
      } 
   } 

   if (found_link_type == 0)
      std::cout << "link type \"" << link_type << "\" not found in dictionary!!\n";

   // std::cout << "add link bond returns " <<  nbond << std::endl;
   return nbond;
}


// Note that we have to convert between fixed residues and fixed
// atoms.  It was easier with bonds of course where there was a 1:1
// relationship.  Residue 1 of the link was the first atom of the link
int
coot::restraints_container_t::add_link_angle(std::string link_type,
					     mmdb::PResidue first, mmdb::PResidue second,
					     short int is_fixed_first,  // residue
					     short int is_fixed_second, // residue
					     const coot::protein_geometry &geom) {

   int nangle = 0;

   mmdb::PPAtom first_sel;
   mmdb::PPAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   int n_atom_1, n_atom_2, n_atom_3;
   int index1, index2, index3;

   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);

   mmdb::PPAtom atom_1_sel, atom_2_sel, atom_3_sel; // assigned to either
					       // first_sel or
					       // second_sel when
					       // atom_1_comp_id
					       // (etc.) are known.

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   std::vector<bool> fixed_flag(3);
   fixed_flag[0] = 0;  // not fixed
   fixed_flag[1] = 0;
   fixed_flag[2] = 0;

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically TRANS
	 for (unsigned int j=0; j<geom.link(i).link_angle_restraint.size(); j++) {

	    if (geom.link(i).link_angle_restraint[j].atom_1_comp_id == 1) {
	       atom_1_sel = first_sel;
	       n_atom_1   = n_first_res_atoms;
	       fixed_flag[0] = is_fixed_first;
	    } else {
	       atom_1_sel = second_sel;
	       n_atom_1   = n_second_res_atoms;
	       fixed_flag[0] = is_fixed_second;
	    }
	    if (geom.link(i).link_angle_restraint[j].atom_2_comp_id == 1) {
	       atom_2_sel = first_sel;
	       n_atom_2   = n_first_res_atoms;
	       fixed_flag[1] = is_fixed_first;
	    } else {
	       atom_2_sel = second_sel;
	       n_atom_2   = n_second_res_atoms;
	       fixed_flag[1] = is_fixed_second;
	    }
	    if (geom.link(i).link_angle_restraint[j].atom_3_comp_id == 1) {
	       atom_3_sel = first_sel;
	       n_atom_3   = n_first_res_atoms;
	       fixed_flag[2] = is_fixed_first;
	    } else {
	       atom_3_sel = second_sel;
	       n_atom_3   = n_second_res_atoms;
	       fixed_flag[2] = is_fixed_second;
	    }
	    
	    for (int ifat=0; ifat<n_atom_1; ifat++) { 
	       std::string pdb_atom_name_1(atom_1_sel[ifat]->name);

	       if (pdb_atom_name_1 == geom.link(i).link_angle_restraint[j].atom_id_1_4c()) {
		  for (int isat=0; isat<n_atom_2; isat++) { 
		     std::string pdb_atom_name_2(atom_2_sel[isat]->name);

		     if (pdb_atom_name_2 == geom.link(i).link_angle_restraint[j].atom_id_2_4c()) {
			for (int itat=0; itat<n_atom_3; itat++) { 
			   std::string pdb_atom_name_3(atom_3_sel[itat]->name);
			   
			   if (pdb_atom_name_3 == geom.link(i).link_angle_restraint[j].atom_id_3_4c()) {

// 			      std::cout << "INFO: res "
// 					<< atom_1_sel[ifat]->residue->seqNum << " "
// 					<< atom_1_sel[ifat]->name << " to "
// 					<< atom_2_sel[isat]->residue->seqNum << " "
// 					<< atom_2_sel[isat]->name << " to "
// 					<< atom_3_sel[itat]->residue->seqNum << " "
// 					<< atom_3_sel[itat]->name << std::endl;
				 
//  			      int index1_old = get_asc_index(pdb_atom_name_1,
//  							 atom_1_sel[ifat]->residue->seqNum,
//  							 atom_1_sel[ifat]->residue->GetChainID());
			
//  			      int index2_old = get_asc_index(pdb_atom_name_2,
//  							 atom_2_sel[isat]->residue->seqNum,
//  							 atom_2_sel[isat]->residue->GetChainID());

//  			      int index3_old = get_asc_index(pdb_atom_name_3,
//  							 atom_3_sel[itat]->residue->seqNum,
//  							 atom_3_sel[itat]->residue->GetChainID());

			      std::string alt_conf_1 = atom_1_sel[ifat]->altLoc;
			      std::string alt_conf_2 = atom_2_sel[isat]->altLoc;
			      std::string alt_conf_3 = atom_3_sel[itat]->altLoc;

			      // either they are all the same (including the ususal case of all "")
			      // or at_1 and at_2 are the same and at_3 is blank
			      // or at_1 is non-blank and at_2 and at_2 are blank
			      // 
			      if (((alt_conf_1 == alt_conf_2) && (alt_conf_1 == alt_conf_3)) ||
				  ((alt_conf_1 == alt_conf_2) && (alt_conf_3 == "")) ||
				  ((alt_conf_2 == alt_conf_3) && (alt_conf_1 == "")) ||
				  ((alt_conf_1 == "") && (alt_conf_2 == "")) ||
				  ((alt_conf_2 == "") && (alt_conf_3 == ""))

				  ) {

				 // set the UDD flag for this residue being bonded/angle with 
				 // the other
				    
				 atom_1_sel[ifat]->GetUDData(udd_atom_index_handle, index1);
				 atom_2_sel[isat]->GetUDData(udd_atom_index_handle, index2);
				 atom_3_sel[itat]->GetUDData(udd_atom_index_handle, index3);
			      
				 if (false) {  // debug
				    std::cout << "bonded_atom_indices.size(): "
					      <<  bonded_atom_indices.size() << std::endl;
				    std::cout << "   add_link_angle: "  << " " << index1 << std::endl;
				    std::cout << "   add_link_angle: "  << " " << index2 << std::endl;
				    std::cout << "   add_link_angle: "  << " " << index3 << std::endl;

				    std::cout << "adding link angle: "
					      << atom_spec_t(atom[index1]) << " "
					      << atom_spec_t(atom[index2]) << " "
					      << atom_spec_t(atom[index3]) << std::endl;
				    
				    std::cout << "same as? test    : "
					      << atom_spec_t(atom_1_sel[ifat]) << " "
					      << atom_spec_t(atom_2_sel[isat]) << " "
					      << atom_spec_t(atom_3_sel[itat]) << std::endl;
				 }

			     
				 bonded_atom_indices[index1].insert(index3);
				 bonded_atom_indices[index3].insert(index1);

				 std::vector<bool> other_fixed_flags = make_fixed_flags(index1,
											index2,
											index3);
				 for (unsigned int ii=0; ii<other_fixed_flags.size(); ii++)
				    if (other_fixed_flags[ii])
				       fixed_flag[ii] = 1;

				 bool is_single_H_atom_angle_restraint = false;
				 unsigned int nH = 0;
				 if (is_hydrogen(atom_1_sel[ifat])) nH++;
				 if (is_hydrogen(atom_3_sel[itat])) nH++;
				 if (nH == 1) is_single_H_atom_angle_restraint = true;

				 add(ANGLE_RESTRAINT, index1, index2, index3,
				     fixed_flag,
				     geom.link(i).link_angle_restraint[j].angle(),
				     geom.link(i).link_angle_restraint[j].angle_esd(),
				     is_single_H_atom_angle_restraint);
				 nangle++;
			      } 
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return nangle;
}


int
coot::restraints_container_t::add_link_torsion(std::string link_type,
                                               mmdb::Residue *first,
                                               mmdb::Residue *second,
                                               short int is_fixed_first,
                                               short int is_fixed_second,
                                               const coot::protein_geometry &geom) {
   int n_torsions = 0;

   mmdb::PAtom *first_sel = 0;
   mmdb::PAtom *second_sel = 0;
   int n_first_res_atoms, n_second_res_atoms;
   int n_atom_1, n_atom_2, n_atom_3, n_atom_4;
   int index1, index2, index3, index4;

   first->GetAtomTable(first_sel,   n_first_res_atoms);
   second->GetAtomTable(second_sel, n_second_res_atoms);

   mmdb::PPAtom atom_1_sel, atom_2_sel, atom_3_sel, atom_4_sel; // assigned to either
                                                                   // first_sel or
                                                                // second_sel when
                                                                // atom_1_comp_id
                                                                // (etc.) are known.

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   std::vector<bool> fixed_flag(4);
   fixed_flag[0] = 0;  // not fixed
   fixed_flag[1] = 0;
   fixed_flag[2] = 0;
   fixed_flag[3] = 0;

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically TRANS
         for (unsigned int j=0; j<geom.link(i).link_torsion_restraint.size(); j++) {

            const auto &ltr = geom.link(i).link_torsion_restraint[j];

            if (ltr.atom_1_comp_id == 1) {
               atom_1_sel = first_sel;
               n_atom_1 = n_first_res_atoms;
               fixed_flag[0] = is_fixed_first;
            } else {
               atom_1_sel = second_sel;
               n_atom_1 = n_second_res_atoms;
               fixed_flag[0] = is_fixed_second;
            }

            if (ltr.atom_2_comp_id == 1) {
               atom_2_sel = first_sel;
               n_atom_2 = n_first_res_atoms;
               fixed_flag[1] = is_fixed_first;
            } else {
               atom_2_sel = second_sel;
               n_atom_2 = n_second_res_atoms;
               fixed_flag[1] = is_fixed_second;
            }

            if (ltr.atom_3_comp_id == 1) {
               atom_3_sel = first_sel;
               n_atom_3 = n_first_res_atoms;
               fixed_flag[2] = is_fixed_first;
            } else {
               atom_3_sel = second_sel;
               n_atom_3 = n_second_res_atoms;
               fixed_flag[2] = is_fixed_second;
            }

            if (ltr.atom_4_comp_id == 1) {
               atom_4_sel = first_sel;
               n_atom_4 = n_first_res_atoms;
               fixed_flag[3] = is_fixed_first;
            } else {
               atom_4_sel = second_sel;
               n_atom_4 = n_second_res_atoms;
               fixed_flag[3] = is_fixed_second;
            }

            for (int ifat=0; ifat<n_atom_1; ifat++) {
               std::string pdb_atom_name_1(atom_1_sel[ifat]->GetAtomName());

               if (pdb_atom_name_1 == ltr.atom_id_1_4c()) {

                  for (int isat=0; isat<n_atom_2; isat++) {
                     std::string pdb_atom_name_2(atom_2_sel[isat]->GetAtomName());

                     if (pdb_atom_name_2 == ltr.atom_id_2_4c()) {

                        for (int itat=0; itat<n_atom_3; itat++) {
                           std::string pdb_atom_name_3(atom_3_sel[itat]->GetAtomName());

                           if (pdb_atom_name_3 == ltr.atom_id_3_4c()) {

                              for (int iffat=0; iffat<n_atom_4; iffat++) {
                                 std::string pdb_atom_name_4(atom_4_sel[iffat]->GetAtomName());

                                 if (pdb_atom_name_4 == ltr.atom_id_4_4c()) {

                                    mmdb::Atom *atom_1 = atom_1_sel[ifat];
                                    mmdb::Atom *atom_2 = atom_2_sel[isat];
                                    mmdb::Atom *atom_3 = atom_3_sel[itat];
                                    mmdb::Atom *atom_4 = atom_4_sel[iffat];

                                    int index_1 = -1, index_2 = -1, index_3 = -1, index_4 = -1;
                                    atom_1->GetUDData(udd_atom_index_handle, index_1);
                                    atom_2->GetUDData(udd_atom_index_handle, index_2);
                                    atom_3->GetUDData(udd_atom_index_handle, index_3);
                                    atom_4->GetUDData(udd_atom_index_handle, index_4);

                                    // skip dictionary mainchain torsions
                                    if (pdb_atom_name_1 == " N  " && pdb_atom_name_4 == " N  ") continue;
                                    if (pdb_atom_name_1 == " CA " && pdb_atom_name_4 == " CA ") continue;
                                    if (pdb_atom_name_1 == " C  " && pdb_atom_name_4 == " C  ") continue;

				    if (false)
				       std::cout << "adding link torsion "
						 << coot::atom_spec_t(atom_1) << " "
						 << coot::atom_spec_t(atom_2) << " "
						 << coot::atom_spec_t(atom_3) << " "
						 << coot::atom_spec_t(atom_4) << " "
						 << ltr.angle() << " " << ltr.period()
						 << std::endl;

                                    add(TORSION_RESTRAINT, index_1, index_2, index_3, index_4,
                                        fixed_flag, ltr.angle(), ltr.angle_esd(), 1.2, ltr.period());
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   return n_torsions;
}

int
coot::restraints_container_t::add_link_torsion_for_phi_psi(std::string link_type,
                                                           int phi_psi_restraints_type,
                                                           mmdb::Residue *first,
                                                           mmdb::Residue *second,
                                                           short int is_fixed_first,
                                                           short int is_fixed_second,
                                                           const coot::protein_geometry &geom) {
   // link_type is "p", "TRANS" etc.

//    std::cout << "--------- :: Adding link torsion, link_type: " << link_type << " phi_psi_restraints_type: "
//          << phi_psi_restraints_type << std::endl;

   int n_torsion = 0;

   mmdb::PPAtom first_sel;
   mmdb::PPAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   int n_atom_1, n_atom_2, n_atom_3, n_atom_4;

   first->GetAtomTable(first_sel,   n_first_res_atoms);
   second->GetAtomTable(second_sel, n_second_res_atoms);

   // assigned to either first_sel or second_sel when atom_1_comp_id
   // (etc.) are known.
   mmdb::PPAtom atom_1_sel, atom_2_sel, atom_3_sel, atom_4_sel;

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   std::vector<bool> fixed_flag(4);
   fixed_flag[0] = 0;  // not fixed
   fixed_flag[1] = 0;
   fixed_flag[2] = 0;
   fixed_flag[3] = 0;

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) {
	 for (unsigned int j=0; j<geom.link(i).link_torsion_restraint.size(); j++) {

	    // This could have been more compact if we were using a
	    // vector of atom ids... (or lisp)... heyho.
	    // 
	    if (geom.link(i).link_torsion_restraint[j].atom_1_comp_id == 1) {
	       atom_1_sel = first_sel;
	       n_atom_1 = n_first_res_atoms;
	       fixed_flag[0] = is_fixed_first;
	    } else {
	       atom_1_sel = second_sel;
	       n_atom_1 = n_second_res_atoms; 
	       fixed_flag[0] = is_fixed_second;
	    }
	    if (geom.link(i).link_torsion_restraint[j].atom_2_comp_id == 1) {
	       atom_2_sel = first_sel;
	       n_atom_2 = n_first_res_atoms; 
	       fixed_flag[1] = is_fixed_first;
	    } else {
	       atom_2_sel = second_sel;
	       n_atom_2 = n_second_res_atoms; 
	       fixed_flag[1] = is_fixed_second;
	    }
	    if (geom.link(i).link_torsion_restraint[j].atom_3_comp_id == 1) {
	       atom_3_sel = first_sel;
	       n_atom_3 = n_first_res_atoms; 
	       fixed_flag[2] = is_fixed_first;
	    } else {
	       atom_3_sel = second_sel;
	       n_atom_3 = n_second_res_atoms; 
	       fixed_flag[2] = is_fixed_second;
	    }
	    if (geom.link(i).link_torsion_restraint[j].atom_4_comp_id == 1) {
	       n_atom_4 = n_first_res_atoms; 
	       atom_4_sel = first_sel;
	       fixed_flag[3] = is_fixed_first;
	    } else {
	       atom_4_sel = second_sel;
	       n_atom_4 = n_second_res_atoms; 
	       fixed_flag[3] = is_fixed_second;
	    }
	    for (int ifat=0; ifat<n_atom_1; ifat++) { 
	       std::string pdb_atom_name_1(atom_1_sel[ifat]->name);
	       
	       if (pdb_atom_name_1 == geom.link(i).link_torsion_restraint[j].atom_id_1_4c()) {
		  for (int isat=0; isat<n_atom_2; isat++) { 
		     std::string pdb_atom_name_2(atom_2_sel[isat]->name);
		     
		     if (pdb_atom_name_2 == geom.link(i).link_torsion_restraint[j].atom_id_2_4c()) {
			for (int itat=0; itat<n_atom_3; itat++) { 
			   std::string pdb_atom_name_3(atom_3_sel[itat]->name);
			   
			   if (pdb_atom_name_3 == geom.link(i).link_torsion_restraint[j].atom_id_3_4c()) {
			      for (int iffat=0; iffat<n_atom_4; iffat++) {
				 std::string pdb_atom_name_4(atom_4_sel[iffat]->name);
				 			   
				 if (pdb_atom_name_4 == geom.link(i).link_torsion_restraint[j].atom_id_4_4c()) {
				    
				    int index1 = get_asc_index(atom_1_sel[ifat]->name,
							       atom_1_sel[ifat]->altLoc,
							       atom_1_sel[ifat]->residue->seqNum,
							       atom_1_sel[ifat]->GetInsCode(),
							       atom_1_sel[ifat]->GetChainID());
			
				    int index2 = get_asc_index(atom_2_sel[isat]->name,
							       atom_2_sel[isat]->altLoc,
							       atom_2_sel[isat]->residue->seqNum,
							       atom_2_sel[isat]->GetInsCode(),
							       atom_2_sel[isat]->GetChainID());
				    
				    int index3 = get_asc_index(atom_3_sel[itat]->name,
							       atom_3_sel[itat]->altLoc,
							       atom_3_sel[itat]->residue->seqNum,
							       atom_3_sel[itat]->GetInsCode(),
							       atom_3_sel[itat]->GetChainID());

				    int index4 = get_asc_index(atom_4_sel[iffat]->name,
							       atom_4_sel[iffat]->altLoc,
							       atom_4_sel[iffat]->residue->seqNum,
							       atom_4_sel[iffat]->GetInsCode(),
							       atom_4_sel[iffat]->residue->GetChainID());

//  				    std::cout << "torsion restraint.... " << geom.link(i).link_torsion_restraint[j].id()
//  					      << " from atoms \n    "
// 					      << atom_1_sel[ifat]->name << " " 
//  					      << atom_1_sel[ifat]->GetSeqNum() << "\n    " 
//  					      << atom_2_sel[isat]->name << " " 
//  					      << atom_2_sel[isat]->GetSeqNum() << "\n    " 
//  					      << atom_3_sel[itat]->name << " " 
//  					      << atom_3_sel[itat]->GetSeqNum() << "\n    " 
//  					      << atom_4_sel[iffat]->name << " " 
//  					      << atom_4_sel[iffat]->GetSeqNum() << "\n";
				       

				    double target_phi = -57.82 + 360.0;
				    double target_psi = -47.0  + 360.0;
				    double esd = 5.01;
				    if (phi_psi_restraints_type == coot::restraints_container_t::LINK_TORSION_ALPHA_HELIX) {
				       // theortical helix CA-N bond: -57.82 phi
				       //                  CA-C bond: -47    psi
				       target_phi = -57.82 + 360.0;
				       target_psi = -47.00 + 360.0;
				    }
				       
				    if (phi_psi_restraints_type == coot::restraints_container_t::LINK_TORSION_BETA_STRAND) {
				       // beta strand
				       target_phi = -110.0 + 360.0;
				       target_psi =  120.0; // approx values
				       esd = 15.0; // guess
				    }
				    
				    if (geom.link(i).link_torsion_restraint[j].id() == "phi") { 
				       add(TORSION_RESTRAINT, index1, index2, index3, index4,
					   fixed_flag,
					   target_phi,
					   esd,
					   1.2, // junk value
					   1);
				       n_torsion++;
				       // std::cout << "!!!!!!!!!!!!!!!! added link torsion restraint phi" << std::endl;
				    }
				    if (geom.link(i).link_torsion_restraint[j].id() == "psi") { 
				       add(TORSION_RESTRAINT, index1, index2, index3, index4,
					   fixed_flag,
					   target_psi,
					   esd, 
					   1.2, // junk value (obs)
					   1);
				       //   std::cout << "!!!!!!!!!!!!!!!! added link torsion restraint psi" << std::endl;
				       n_torsion++;
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_torsion; 
}


// Using bonded pairs internal copy, modify residues as needed
// by deleting atoms in chem mods.
//
void
coot::restraints_container_t::apply_link_chem_mods(const coot::protein_geometry &geom) {

   bonded_pairs_container.apply_chem_mods(geom);
} 

int
coot::restraints_container_t::make_link_restraints(const coot::protein_geometry &geom,
						   bool do_rama_plot_restraints,
						   bool do_trans_peptide_restraints) {

   std::cout << "debug:: make_link_restraints() ***********************************************************"
             << " I don't think that this function is called. "<< std::endl;

   if (from_residue_vector) {
      // coot::bonded_pair_container_t bonded_pairs_container;
      bonded_pairs_container = make_link_restraints_from_res_vec(geom,
								 do_rama_plot_restraints,
								 do_trans_peptide_restraints);
      return bonded_pairs_container.size();
   } else { 
      return make_link_restraints_by_linear(geom,
					    do_rama_plot_restraints,
					    do_trans_peptide_restraints); // conventional
   }
   
   int ir = 0;
   return ir; 
}

int
coot::restraints_container_t::make_link_restraints_by_linear(const coot::protein_geometry &geom,
							     bool do_rama_plot_restraints,
							     bool do_trans_peptide_restraints) {


   // Last time (for monomer geometry), we got a residue and added
   // bonds, angles and torsions by checking the atoms of that single
   // residue.
   //
   // This time, we need 2 consecutive residues, checking the atom
   // types of each - note that the atoms have to correspond to the
   // correct comp_id for that atom.
   //
   int selHnd = mol->NewSelection();
   mmdb::PPResidue     SelResidue;
   int nSelResidues;


//    timeval start_time;
//    timeval current_time;
//    double td;


   mol->Select ( selHnd,mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
		 chain_id_save.c_str(), // Chain(s)
		 istart_res, "*",  // starting res
		 iend_res,   "*",  // ending res
		 "*",  // residue name
		 "*",  // Residue must contain this atom name?
		 "*",  // Residue must contain this Element?
		 "*",  // altLocs
		 mmdb::SKEY_NEW // selection key
		 );
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   std::cout << "INFO:: GetSelIndex (make_link_restraints) returned " << nSelResidues
	     << " residues (for link restraints) between (and including) residues "
	     << istart_res << " and " << iend_res << " of chain " << chain_id_save
	     << std::endl;
   
   coot::bonded_pair_container_t bonded_residue_pairs =
      bonded_residues_conventional(selHnd, geom);

   if (0)
      std::cout << "debug -------- in make_link_restraints_by_linear() "
		<< bonded_residue_pairs.size() << " bonded residue pairs "
		<< std::endl;
   // std::cout << "  " << bonded_residue_pairs << std::endl;

   int iv = make_link_restraints_by_pairs(geom, bonded_residue_pairs, do_trans_peptide_restraints, "Link");

   if (do_rama_plot_restraints) {
      add_rama_links(selHnd, geom); // uses TRANS links
   }
    

   mol->DeleteSelection(selHnd);
   return iv;
}


coot::bonded_pair_container_t
coot::restraints_container_t::make_link_restraints_from_res_vec(const coot::protein_geometry &geom,
                                                                bool do_rama_plot_restraints,
                                                                bool do_trans_peptide_restraints) {

   // std::cout << "DEBUG:: starting make_link_restraints_from_res_vec() ::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;

   // 2025-10-14-PE do these arguments make sense?
   // return make_link_restraints_by_distance(geom, do_rama_plot_restraints, do_trans_peptide_restraints);

   // return make_link_restraints_from_links(geom);

   return coot::bonded_pair_container_t();
}

coot::bonded_pair_container_t
coot::restraints_container_t::make_link_restraints_from_links(const coot::protein_geometry &geom) {

   auto choose_a_link = [] (const std::string &comp_id_1, const std::string &comp_id_2,
                            const std::string &group_1, const std::string &group_2,
                            const std::vector<chem_link> &link_infos_f, const std::vector<chem_link> &link_infos_b,
                            const std::pair<atom_spec_t, atom_spec_t> &l_atoms,
                            const coot::protein_geometry &geom) -> std::optional<dictionary_residue_link_restraints_t> {

      // std::cout << "choose_a_link(): " << l_atoms.first << " " << l_atoms.second << std::endl;
      std::vector<dictionary_residue_link_restraints_t> matched_links; // by Id and atom names
      // look in geom.dict_link_res_restraints for a matching link restraint.
      for (unsigned int i=0; i<link_infos_f.size(); i++) {
         std::string id = link_infos_f[i].Id();
         dictionary_residue_link_restraints_t lr = geom.link(id);
         if (! lr.empty()) {
            // std::cout << "   non-empty link found " << lr.link_id << std::endl;
            for (unsigned int j=0; j<lr.link_bond_restraint.size(); j++) {
               const auto &lbr = lr.link_bond_restraint[j];
               if (false)
                  std::cout << "   compare " << id << " " << l_atoms.first << " " << l_atoms.second << " "
                            << lbr.atom_id_1() << " " << lbr.atom_id_2() << " " << std::endl;
               std::string a1 = util::remove_whitespace(l_atoms.first.atom_name);
               std::string a2 = util::remove_whitespace(l_atoms.second.atom_name);
               bool m = lbr.matches(a1, a2);
               if (m) {
                  matched_links.push_back(lr);
               }
            }
         }
      }

      // std::cout << "DEBUG:: found " << matched_links.size() << " matched links" << std::endl;

      // this is pretty hacky!
      if (matched_links.empty()) {
         return std::nullopt;
      } else {
         unsigned int idx_pref = 0;
         if (matched_links.size() > 1) {
            if (matched_links[1].link_id == "BETA1-4")
               idx_pref = 1;
         }
         return matched_links[idx_pref];
      }
   };

   // std::cout << "DEBUG:: starting make_link_restraints_from_links() ::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;

   coot::bonded_pair_container_t v;
   // go through the list of links in the manager and see if there are links that match
   // residues in our selection.
   int imod = 1;
   if (mol) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         mmdb::LinkContainer *links = model_p->GetLinks();
         int n_links = model_p->GetNumberOfLinks();
         for (int ilink=0; ilink<=n_links; ilink++) {
            mmdb::Link *link_p = model_p->GetLink(ilink);
            if (link_p) {
               std::pair<atom_spec_t, atom_spec_t> l_atoms = link_atoms(link_p, model_p);
               residue_spec_t res_spec_1(l_atoms.first);
               residue_spec_t res_spec_2(l_atoms.second);
               mmdb::Residue *first_residue  = util::get_residue(res_spec_1, mol);
               mmdb::Residue *second_residue = util::get_residue(res_spec_2, mol);
               if (first_residue) {
                  if (second_residue) {
                     std::string comp_id_1 =  first_residue->GetResName();
                     std::string comp_id_2 = second_residue->GetResName();
                     std::string group_1 = geom.get_group( first_residue);
                     std::string group_2 = geom.get_group(second_residue);
                     if (group_1 == "DNA") group_1 = "DNA/RNA";
                     if (group_1 == "RNA") group_1 = "DNA/RNA";
                     if (group_2 == "DNA") group_2 = "DNA/RNA";
                     if (group_2 == "RNA") group_2 = "DNA/RNA";
                     // handle share/coot/data/cho-acedrg/NAG-acedrg.cif and friends (20230720-PE they should be replaced now).
                     if (group_1 == "D-SACCHARIDE") group_1 = "pyranose";
                     if (group_2 == "D-SACCHARIDE") group_2 = "pyranose";
                     std::vector<coot::chem_link> link_infos_f = geom.matching_chem_links(comp_id_1, group_1, comp_id_2, group_2);
                     std::vector<coot::chem_link> link_infos_b = geom.matching_chem_links(comp_id_2, group_2, comp_id_1, group_1);

                     if (false) {
                        std::cout << "DEBUG:: make_link_restraints_from_links(): found " << link_infos_f.size() << " forward links for "
                                  << comp_id_1 << " " << comp_id_2 << std::endl;
                        std::cout << "DEBUG:: make_link_restraints_from_links(): found " << link_infos_b.size() << " backward links for "
                                  << comp_id_1 << " " << comp_id_2 << std::endl;
                        for (unsigned int i=0; i<link_infos_f.size(); i++)
                           std::cout << "DEBUG:: make_link_restraints_from_links(): link_infos_f " << link_infos_f[i] << std::endl;
                        for (unsigned int i=0; i<link_infos_b.size(); i++)
                           std::cout << "DEBUG:: make_link_restraints_from_links(): link_infos_b " << link_infos_b[i] << std::endl;
                     }

                     std::optional<dictionary_residue_link_restraints_t> lr =
                        choose_a_link(comp_id_1, comp_id_2, group_1, group_2, link_infos_f, link_infos_b, l_atoms, geom);
                     if (lr) {

                        bool fixed_1 = false;
                        bool fixed_2 = false;
                        std::string link_id = lr.value().link_id;
                        // 2025-10-22-PE now find if the residues of this link were fixed. Is there a better way than this
                        // code block!?
                        for (unsigned int i=0; i<residues_vec.size(); i++) {
                           if (residues_vec[i].second ==  first_residue) if (residues_vec[i].first) fixed_1 = true;
                           if (residues_vec[i].second == second_residue) if (residues_vec[i].first) fixed_2 = true;
                        }

                        if (false)
                           std::cout << "OK, eventually calliing make_link_restraints_for_link_ng() for "
                                     << coot::residue_spec_t(first_residue) << " and " << coot::residue_spec_t(second_residue)
                                     << " " << link_id << std::endl;

                        make_link_restraints_for_link_ng(link_id, first_residue, second_residue, fixed_1, fixed_2, false, geom);
                     }
                  }
               }
            }
         }
      }
   }
   // std::cout << "DEBUG:: done make_link_restraints_from_links() ::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
   return v;
}

coot::bonded_pair_container_t
coot::restraints_container_t::make_link_restraints_by_distance(const coot::protein_geometry &geom,
                                                               bool do_rama_plot_restraints,
                                                               bool do_trans_peptide_restraints) {

   // this determines the link type

   auto tp_0 = std::chrono::high_resolution_clock::now();
   bonded_pair_container_t bonded_residue_pairs = bonded_residues_from_res_vec(geom); // calc or get?
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   std::cout << "INFO:: Timing for bonded_residues_from_res_vec " << d10 << " milliseconds" << std::endl;

   if (false) {
      std::cout << "   DEBUG:: in make_link_restraints_from_res_vec() found "
		<< bonded_residue_pairs.size() << " bonded residues " << std::endl;
      for (std::size_t i=0; i<bonded_residue_pairs.size(); i++) {
	 const bonded_pair_t &bp = bonded_residue_pairs.bonded_residues[i];
	 std::cout << "   " << i << " " << residue_spec_t(bp.res_1) << " " << residue_spec_t(bp.res_2)
		   << std::endl;
      }
   }

   auto tp_2 = std::chrono::high_resolution_clock::now();
   int iv = make_link_restraints_by_pairs(geom, bonded_residue_pairs, do_trans_peptide_restraints, "Link");
   auto tp_3 = std::chrono::high_resolution_clock::now();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   std::cout << "INFO:: Timing for make_link_restraints_by_pairs " << d32 << " milliseconds" << std::endl;

   if (do_rama_plot_restraints)
      add_rama_links_from_res_vec(bonded_residue_pairs, geom);

   return bonded_residue_pairs;
}

#include <iomanip>

int
coot::restraints_container_t::make_link_restraints_by_pairs(const coot::protein_geometry &geom,
							    const coot::bonded_pair_container_t &bonded_residue_pairs,
							    bool do_trans_peptide_restraints,
							    std::string link_flank_link_string) {

   bool debug = false;
   int iret = 0;
   int n_link_bond_restr = 0;
   int n_link_angle_restr = 0;
   int n_link_torsion_restr = 0;
   int n_link_trans_peptide = 0;
   int n_link_plane_restr = 0;
   int n_link_parallel_plane_restr = 0;

   if (false) {
      std::cout << "Here are the bonded pairs: " << std::endl;
      for (unsigned int ibonded_residue=0;
	   ibonded_residue<bonded_residue_pairs.size();
	   ibonded_residue++) {

	 std::string link_type = bonded_residue_pairs[ibonded_residue].link_type;
	 mmdb::Residue *sel_res_1 = bonded_residue_pairs[ibonded_residue].res_1;
	 mmdb::Residue *sel_res_2 = bonded_residue_pairs[ibonded_residue].res_2;

	 std::cout << " ------- looking for link :" << link_type
		   << ": restraints etc. between residues "
		   << residue_spec_t(sel_res_1) << " " << sel_res_1->GetResName()
		   << " - "
		   << residue_spec_t(sel_res_2) << " " << sel_res_2->GetResName()
		   << std::endl;
      }
   }
   // std::cout << "---------- Done bonded pairs" << std::endl;

   for (unsigned int ibonded_residue=0;
	ibonded_residue<bonded_residue_pairs.size();
	ibonded_residue++) {
      
      // now select pairs of residues: but not for the last one, which
      // doesn't have a peptide link on its C atom.
      //
      // It has a modification there, I suppose, but I am ignoring
      // that for now.
      //
      //
      std::string link_type = bonded_residue_pairs[ibonded_residue].link_type;

      mmdb::Residue *sel_res_1 = bonded_residue_pairs[ibonded_residue].res_1;
      mmdb::Residue *sel_res_2 = bonded_residue_pairs[ibonded_residue].res_2;

      if (verbose_geometry_reporting == VERBOSE) {
	 std::cout << " ------- looking for link :" << link_type
		   << ": restraints etc. between residues "
		   << residue_spec_t(sel_res_1) << " " << sel_res_1->GetResName()
		   << " - "
		   << residue_spec_t(sel_res_2) << " " << sel_res_2->GetResName()
		   << std::endl;
      }
	 

      // Are these residues neighbours?  We can add some sort of
      // ins code test here in the future.
      
      //if ( abs(sel_res_1->GetSeqNum() - sel_res_2->GetSeqNum()) <= 1
      // || abs(sel_res_1->index - sel_res_2->index) <= 1) {

      {

	 bool is_fixed_first_residue  = bonded_residue_pairs[ibonded_residue].is_fixed_first;
	 bool is_fixed_second_residue = bonded_residue_pairs[ibonded_residue].is_fixed_second;

	 // link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);

	 if (restraints_usage_flag & BONDS_MASK)
	    n_link_bond_restr += add_link_bond(link_type,
					       sel_res_1, sel_res_2,
					       is_fixed_first_residue,
					       is_fixed_second_residue,
					       geom);

	 if (restraints_usage_flag & ANGLES_MASK)
	    n_link_angle_restr += add_link_angle(link_type,
						 sel_res_1, sel_res_2,
						 is_fixed_first_residue,
						 is_fixed_second_residue,
						 geom);

	 if (restraints_usage_flag & TRANS_PEPTIDE_MASK)
	    if (do_trans_peptide_restraints)
	       n_link_trans_peptide += add_link_trans_peptide(sel_res_1, sel_res_2,
							      is_fixed_first_residue,
							      is_fixed_second_residue, false); // don't add if cis
	 
	 // 	    gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t1 = td;

	 if (restraints_usage_flag & PLANES_MASK)
	    n_link_plane_restr += add_link_plane(link_type,
						 sel_res_1, sel_res_2,
						 is_fixed_first_residue,
						 is_fixed_second_residue,
						 geom);

	 // 	    gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t2 = td;

	 // link_torsions_type_flag is the type of phi/psi restraints
	 // std::cout << "---------------link_torsions_type_flag: "
	 // << link_torsions_type_flag
	 // << std::endl;

	 // gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t3 = td;

	 // 	    std::cout << "after bonds    " << t0 << std::endl;
	 // 	    std::cout << "after angles   " << t1 << std::endl;
	 // 	    std::cout << "after planes   " << t2 << std::endl;
	 // 	    std::cout << "after torsions " << t3 << std::endl;
      }
   }

   if (verbose_geometry_reporting != QUIET) { 
      std::cout << link_flank_link_string << " restraints: " << std::endl;
      std::cout << "   " << n_link_bond_restr    << " bond    links" << std::endl;
      std::cout << "   " << n_link_angle_restr   << " angle   links" << std::endl;
      std::cout << "   " << n_link_plane_restr   << " plane   links" << std::endl;
      std::cout << "   " << n_link_trans_peptide << " trans-peptide links";
      if (! do_trans_peptide_restraints)
	 std::cout << " (not requested)";
      std::cout << std::endl;
      std::cout << "   " << n_link_parallel_plane_restr   << " parallel plane restraints"
		<< std::endl;
   }
   return iret; 
}

// Uses TRANS links only
void
coot::restraints_container_t::add_rama_links(int selHnd, const coot::protein_geometry &geom) {
   // 
   // Note do_rama_plot_restraints is true before we come here.
   
   int n_link_torsion_restr = 0;
   std::vector<coot::rama_triple_t> v = make_rama_triples(selHnd, geom);
   for (unsigned int ir=0; ir<v.size(); ir++) {
      std::string link_type = "TRANS";
      add_rama(link_type,
	       v[ir].r_1, v[ir].r_2, v[ir].r_3,
		  0,0,0, geom); // no residues are fixed.
      n_link_torsion_restr++;
   }
   std::cout << "   " << n_link_torsion_restr << " torsion/rama links" << std::endl;
}

// Uses TRANS and PTRANS links
void
coot::restraints_container_t::add_rama_links_from_res_vec(const coot::bonded_pair_container_t &bonded_residue_pairs,
							  const coot::protein_geometry &geom) {


   std::vector<rama_triple_t> rama_triples;
   //
   // This relies on the bonded_residues in bonded_residue_pairs being
   // ordered (i.e. res_1 is the first residue in the comp_id of the
   // restraints)
   // 
   for (unsigned int i=0; i<bonded_residue_pairs.bonded_residues.size(); i++) { 
      for (unsigned int j=0; j<bonded_residue_pairs.bonded_residues.size(); j++) { 
	 if (i != j) {
	    if (bonded_residue_pairs.bonded_residues[i].res_2 ==
		bonded_residue_pairs.bonded_residues[j].res_1) {
	       const std::string &lt_i = bonded_residue_pairs.bonded_residues[i].link_type;
	       const std::string &lt_j = bonded_residue_pairs.bonded_residues[j].link_type;
	       if (lt_i == "TRANS" || lt_i == "PTRANS") {
		  if (lt_j == "TRANS" || lt_j == "PTRANS") {
		     coot::rama_triple_t rt(bonded_residue_pairs.bonded_residues[i].res_1,
					    bonded_residue_pairs.bonded_residues[i].res_2,
					    bonded_residue_pairs.bonded_residues[j].res_2,
					    lt_j, // link between residues 2 and 3
					    bonded_residue_pairs.bonded_residues[i].is_fixed_first,
					    bonded_residue_pairs.bonded_residues[i].is_fixed_second,
					    bonded_residue_pairs.bonded_residues[j].is_fixed_second);
		     rama_triples.push_back(rt);
		  }
	       }
	    }
	 }
      }
   }

   for (unsigned int ir=0; ir<rama_triples.size(); ir++) {
      const rama_triple_t &rt = rama_triples[ir];
      add_rama(rt.link_type,
	       rt.r_1, rt.r_2, rt.r_3,
	       rt.fixed_1, rt.fixed_2, rt.fixed_3,
	       geom);
   }
   
}


// Simple, residue-name and restraints group type link finder
// 
std::string
coot::restraints_container_t::find_link_type(mmdb::Residue *first,
					     mmdb::Residue *second,
					     const coot::protein_geometry &geom) const {

   // Should return TRANS, PTRANS (PRO-TRANS), CIS, PCIS, p, BETA1-2, BETA1-4 etc.
   std::string link_type(""); // unset

   std::string residue_type_1 = first->name;
   std::string residue_type_2 = second->name;
   if (residue_type_1 == "UNK") residue_type_1 = "ALA"; // hack for KDC.
   if (residue_type_2 == "UNK") residue_type_2 = "ALA";

   std::string t1="";
   std::string t2="";

   for (unsigned int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_1)) {
	 t1 = geom[idr].second.residue_info.group;
	 break;
      }
   }
   for (unsigned int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_2)) {
	 t2 = geom[idr].second.residue_info.group;
	 break;
      }
   }

   if (t1 == "L-peptide" || t1 == "D-peptide" || t1 == "M-peptide" || t1 == "P-peptide" || t1 == "peptide") {
      if (t2 == "L-peptide" || t2 == "D-peptide" || t2 == "M-peptide" || t2 == "P-peptide" || t2 == "peptide") {
	 if (residue_type_2 == "PRO" || residue_type_2 == "HYP") {
	    link_type = "PTRANS";
	 } else {
	    link_type = "TRANS";
	 }
      }
   }

   if (coot::util::is_nucleotide_by_dict(first, geom))
      link_type = "p"; // phosphodiester linkage

   if (t1 == "D-pyranose" || t1 == "D-furanose" || t1 == "L-pyranose" || t1 == "L-furanose" ||
       t1 == "pyranose" || t1 == "furanose" /* new-style refmac dictionary */) { 
      if (t2 == "D-pyranose" || t2 == "D-furanose" || t2 == "L-pyranose" || t2 == "L-furanose" ||
	  t2 == "pyranose" || t2 == "furanose"  /* new-style refmac dictionary */) {
	 bool use_links_in_molecule = false; // because this is the simple version
	 link_type = find_glycosidic_linkage_type(first, second, geom, use_links_in_molecule);
	 std::cout << "INFO:: glycosidic_linkage type :" << link_type << ":\n";
      }
   }

   return link_type;
}

// Return "" on failure to find link.
//
// If this returns "", consider calling this function again, with
// reversed arguments.
//
std::string
coot::restraints_container_t::find_glycosidic_linkage_type(mmdb::Residue *first, mmdb::Residue *second,
							   const protein_geometry &geom,
							   bool use_links_in_molecule) const {

   std::string r;
   if (use_links_in_molecule)
      r = geom.find_glycosidic_linkage_type_by_distance(first, second); // 20150714 FIXME
   else
      r = geom.find_glycosidic_linkage_type(first, second, mol);
   return r;
}

bool
coot::restraints_container_t::link_infos_are_glycosidic_by_name_p(const std::vector<coot::chem_link> &link_infos) const {

   bool is_sweet = false;
   for (unsigned int i=0; i<link_infos.size(); i++) {
      std::string id = link_infos[i].Id();
      if (id.length() > 4) {
         if ((id.substr(0,5) == "ALPHA") ||
             (id.substr(0,4) == "BETA")) {
            is_sweet = 1;
            break;
         }
      }
   }
   return is_sweet;
}

// Return the link type and a residue order switch flag.
// Return link_type as "" if not found.
//
// to help our correct judgement of making the right links, we need
// the the neighbouring residues too (flanking residues) so that
// have_intermediate_residue_by_seqnum() can use them (to check for
// spurious bonding from highly distorted structure).
//
std::pair<std::string, bool>
coot::restraints_container_t::find_link_type_complicado(mmdb::Residue *first,
							mmdb::Residue *second,
							const coot::protein_geometry &geom) const {

   // std::cout << "####################### find_link_type_compli() called with first_residue "
   // << coot::residue_spec_t(first)  << " " <<  first->GetResName() << " and second residue "
   // << coot::residue_spec_t(second) << " " << second->GetResName() << std::endl;

   return find_link_type_2022(first, second, geom);
}


std::pair<std::string, bool>
coot::restraints_container_t::find_link_type_2022(mmdb::Residue *first_residue,
                                                  mmdb::Residue *second_residue,
                                                  const coot::protein_geometry &geom) const {

   auto get_consecutive = [] (mmdb::Residue *first_residue, mmdb::Residue *second_residue) {
      bool state = false;
      int idx_1 =  first_residue->index;
      int idx_2 = second_residue->index;
      int d = idx_2 - idx_1;
      if (d <= 1)
         if (d >= -1)
            state = true;
      return state;
   };

   auto SS_filter = [] (mmdb::Residue *first, mmdb::Residue *second) {

      // return either "SS" or "".

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      bool found = false;
      first->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string name(at->name);
            if (name == " SG ") {
               found = true;
               break;
            }
         }
      }
      if (found) {
         found = false;
         n_residue_atoms = 0;
         residue_atoms = 0;
         second->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string name(at->name);
               if (name == " SG ") {
                  found = true;
                  break; // micro-optimiziation!
               }
            }
         }
      }
      return found ? "SS" : "";
   };

   auto AA_RNA_filter = [] (mmdb::Residue *first, mmdb::Residue *second, bool order_switch_flag) {

      // return either "AA-RNA" or "".

      if (order_switch_flag)
         std::swap(first, second);

      bool found_1 = false;
      bool found_2 = false;
      clipper::Coord_orth pt_1(0,0,0);
      clipper::Coord_orth pt_2(0,0,0);
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      first->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string name(at->name);
            if (name == " C  ") {
               found_1 = true;
               pt_1 = clipper::Coord_orth(at->x, at->y, at->z);
            }
         }
      }
      residue_atoms = 0;
      n_residue_atoms = 0;
      second->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string name(at->name);
            if (name == " O3'") {
               found_2 = true;
               pt_2 = clipper::Coord_orth(at->x, at->y, at->z);
            }
         }
      }
      bool close = false;
      if (found_1 && found_2) {
         double dd = (pt_1 - pt_2).lengthsq();
         double d = std::sqrt(dd);
         if (d < 3.0) close = true;
      }
      return close ? "AA-RNA" : "";
   };

   auto pyr_SER_filter = [] (mmdb::Residue *first, mmdb::Residue *second, bool order_switch_flag) {

      // return either "SER-pyr" or "".

      if (order_switch_flag)
         std::swap(first, second);

      bool found_1 = false;
      bool found_2 = false;
      clipper::Coord_orth pt_1(0,0,0);
      clipper::Coord_orth pt_2(0,0,0);
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      first->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string name(at->name);
            if (name == " C1 ") {
               found_1 = true;
               pt_1 = clipper::Coord_orth(at->x, at->y, at->z);
            }
         }
      }
      residue_atoms = 0;
      n_residue_atoms = 0;
      second->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string name(at->name);
            if (name == " OG ") {
               found_2 = true;
               pt_2 = clipper::Coord_orth(at->x, at->y, at->z);
            }
         }
      }
      bool close = false;
      if (found_1 && found_2) {
         double dd = (pt_1 - pt_2).lengthsq();
         double d = std::sqrt(dd);
         if (d < 3.0) close = true;
      }
      return close ? "pyr-SER" : "";
   };

   // return the link_id if the linking atoms are close (3A)
   // otherwise return empty string
   auto link_type_filter_general = [&geom] (mmdb::Residue *first, mmdb::Residue *second,
                                            const std::string &link_id,
                                            bool order_switch_flag) {

      bool debug = false;
      if (debug)
         std::cout << "link_type_filter_general() starting with link_id " << link_id << std::endl;

      double dist_crit = 3.0; // A
      std::string found_link; // fail initially.

      dictionary_residue_link_restraints_t link = geom.link(link_id);
      if (link.link_id.empty()) {
      } else {
         if (order_switch_flag)
            std::swap(first, second);
         if (link.link_bond_restraint.empty()) {
         } else {
            // very bizarre to have more than one of these
            const dict_link_bond_restraint_t &lbr = link.link_bond_restraint[0];
            bool found_1 = false;
            bool found_2 = false;
            std::string link_atom_1_name = lbr.atom_id_1_4c();
            std::string link_atom_2_name = lbr.atom_id_2_4c();
            clipper::Coord_orth pt_1(0,0,0);
            clipper::Coord_orth pt_2(0,0,0);
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            first->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (! at->isTer()) {
                  std::string name(at->name);
                  if (name == link_atom_1_name) {
                     found_1 = true;
                     pt_1 = clipper::Coord_orth(at->x, at->y, at->z);
                  }
               }
            }
            residue_atoms = 0;
            n_residue_atoms = 0;
            second->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (! at->isTer()) {
                  std::string name(at->name);
                  if (name == link_atom_2_name) {
                     found_2 = true;
                     pt_2 = clipper::Coord_orth(at->x, at->y, at->z);
                  }
               }
            }
            if (found_1 && found_2) {
               double dd = (pt_1 - pt_2).lengthsq();
               double d = std::sqrt(dd);
               if (d < dist_crit) {
                  found_link = link_id;
               }
            }
         }
      }
      if (debug)
         std::cout << "link_type_filter_general() checking type " << link_id << " and returns \""
                   << found_link << "\"" << std::endl;
      return found_link;
   };

   bool debug_links = false;

   if (debug_links) {
      std::cout << "\n####################### find_link_type_2022() called with first_residue       "
                << coot::residue_spec_t(first_residue)  << " " <<  first_residue->GetResName() << " and second residue "
                << coot::residue_spec_t(second_residue) << " " << second_residue->GetResName() << std::endl;
   }

   std::string link_type; // set this
   bool order_switch_was_needed = false; // and this

   try {
      std::string comp_id_1 =  first_residue->GetResName();
      std::string comp_id_2 = second_residue->GetResName();
      std::string group_1 = geom.get_group( first_residue);
      std::string group_2 = geom.get_group(second_residue);
      if (group_1 == "DNA") group_1 = "DNA/RNA";
      if (group_1 == "RNA") group_1 = "DNA/RNA";
      if (group_2 == "DNA") group_2 = "DNA/RNA";
      if (group_2 == "RNA") group_2 = "DNA/RNA";
      // handle share/coot/data/cho-acedrg/NAG-acedrg.cif and friends (20230720-PE they should be replaced now).
      if (group_1 == "D-SACCHARIDE") group_1 = "pyranose";
      if (group_2 == "D-SACCHARIDE") group_2 = "pyranose";

      if (debug_links)
         std::cout << "  "
                   << " comp_id_1 " << comp_id_1 << " group_1 " << group_1
                   << " comp_id_2 " << comp_id_2 << " group_2 " << group_2
                   << std::endl;

      if (group_1 == "pyranose" && group_2 == "pyranose") { // does this link O-linked carbohydrates?
         std::string link_type_glyco;
         bool use_links_in_molecule = true;
         link_type_glyco = find_glycosidic_linkage_type(first_residue, second_residue, geom, use_links_in_molecule);
         if (link_type_glyco.empty()) {
            link_type_glyco = find_glycosidic_linkage_type(second_residue, first_residue, geom, use_links_in_molecule);
            if (! link_type_glyco.empty()) {
               link_type = link_type_glyco;
               order_switch_was_needed =  true;
            }
         } else {
            link_type = link_type_glyco;
         }

      } else {

         std::vector<coot::chem_link> link_infos_f = geom.matching_chem_links(comp_id_1, group_1, comp_id_2, group_2);
         std::vector<coot::chem_link> link_infos_b = geom.matching_chem_links(comp_id_2, group_2, comp_id_1, group_1);

         if (debug_links) {
            std::cout << "   ###### found n-forward  links: " << link_infos_f.size() << std::endl;
            std::cout << "   ###### found n-backward links: " << link_infos_b.size() << std::endl;
         }

         std::vector<std::pair<coot::chem_link, bool> > chem_links;
         for (const auto &link : link_infos_f) chem_links.push_back(std::make_pair(link, false));
         for (const auto &link : link_infos_b) chem_links.push_back(std::make_pair(link,  true));

         if (debug_links) {
            for (unsigned int ilink=0; ilink<chem_links.size(); ilink++) {
               const coot::chem_link &link = chem_links[ilink].first;
               bool order_switch_is_needed = chem_links[ilink].second;
               std::cout << "............... matching links " << ilink << " of " << chem_links.size() << " "
                         << link << " order_switch_is_needed: " << order_switch_is_needed << std::endl;
            }
         }

         if (! chem_links.empty()) {

            bool residues_are_consecutive = get_consecutive(first_residue, second_residue);

            if (residues_are_consecutive) {
               // try to choose TRANS if it is there:
               for (unsigned int ilink=0; ilink<chem_links.size(); ilink++) {
                  const coot::chem_link &link = chem_links[ilink].first;
                  bool order_switch_is_needed = chem_links[ilink].second;
                  if (link.Id() == "TRANS") {
                     link_type = "TRANS";
                     order_switch_was_needed = order_switch_is_needed;
                  }
               }
            } else {

               // It cannot be a peptide bond of any kind then.
               // (20240923-PE unless it is a peptide macrocycle...)

               if (debug_links)
                  std::cout << "   .... find_link_type_2022() here 000 with chem_links size() " << chem_links.size()
                            << " and current link-type :" << link_type << ":" << std::endl;

               // let's try to find a link where the comp_id1 and comp_id2 are not "any"
               for (unsigned int ilink=0; ilink<chem_links.size(); ilink++) {
                  const coot::chem_link &link = chem_links[ilink].first;
                  bool order_switch_is_needed = chem_links[ilink].second;
                  if (! link.chem_link_comp_id_1.empty()) {
                     if (! link.chem_link_comp_id_2.empty()) {
                        // this is our boy (typically an bespoke acedrg link)
                        link_type = link.Id();
                        order_switch_was_needed = order_switch_is_needed;
                     }
                  }
               }

               if (link_type.empty()) {
                  // It might be a "SS" then, choose it if it is in the list.
                  // The options at this stage are SS, TRANS, or CIS.
                  for (unsigned int ilink=0; ilink<chem_links.size(); ilink++) {
                     const coot::chem_link &link = chem_links[ilink].first;
                     bool order_switch_is_needed = chem_links[ilink].second;
                     if (link.Id() == "SS") {
                        link_type = "SS";
                        order_switch_was_needed = order_switch_is_needed;
                     }
                  }
               }
            }

            // 20221120-PE just choose the top one. I can be more clever if/when needed.
            if (link_type.empty()) {
               link_type               = chem_links[0].first.Id();
               order_switch_was_needed = chem_links[0].second;
            }
         }
      }

      if (debug_links)
         std::cout << "   .... find_link_type_2022() here A with link type: " << link_type
                   << " and order_switch_was_needed: " << order_switch_was_needed << std::endl;

   }

   catch (const std::runtime_error &e) {
      std::cout << "WARNING:: no group found in find_link_type_2022() " << std::endl;
      std::cout << e.what() << std::endl;
      return std::pair<std::string, bool> ("", order_switch_was_needed);
   }

   // 20230605-PE
   // the match for SS is:
   // SS . CYS-SS peptide . CYS-SS peptide SS-bridge
   // i.e. comp_id_1 and comp_id_2 are "." i.e. the SS link matches all peptides.
   // We don't want that.
   // Let's filter them out
   if (link_type == "SS")
      link_type = SS_filter(first_residue, second_residue); // return "" if these residues don't contain SG atoms

   if (link_type == "AA-RNA")
      link_type = AA_RNA_filter(first_residue, second_residue, order_switch_was_needed);

   if (link_type == "pyr-SER")
      link_type = pyr_SER_filter(first_residue, second_residue, order_switch_was_needed);

   if (link_type == "p") // 20241126-PE
      if (! get_consecutive(first_residue, second_residue))
         link_type = ""; // nope

   // now check other links (but no need to check the links that we have already checked)
   if (! link_type.empty()) {
      if (link_type != "SS") {
         if (link_type != "AA-RNA") {
            if (link_type != "pyr-SER") {
               if (link_type != "TRANS") {
                  if (link_type != "CIS") {
                     if (link_type != "p") {
                        // if the link is real, this should come back with the link id that it was sent.
                        link_type = link_type_filter_general(first_residue, second_residue, link_type, order_switch_was_needed);
                     }
                  }
               }
            }
         }
      }
   }

   if (debug_links)
      std::cout << "####################### find_link_type_2022() returns \"" << link_type << "\""
                << " order_switch: " << order_switch_was_needed << std::endl;

   return std::pair<std::string, bool> (link_type, order_switch_was_needed);
}



// this test might not be in the right place.
//
// the residues passed to this function are (should be, if you're going to
// reuse it) sorted.
//
// Ideally we should check neighbours of neighbours too (flanking residues)
// 
bool
coot::restraints_container_t::have_intermediate_residue_by_seqnum(mmdb::Residue *first,
								  mmdb::Residue *second) const {

   bool r = false;
   mmdb::Chain *c_1 =  first->GetChain();
   mmdb::Chain *c_2 = second->GetChain();
   if (c_1 == c_2) {
      if (c_1) {
	 int res_no_1 =  first->GetSeqNum();
	 int res_no_2 = second->GetSeqNum();
	 int res_no_diff = res_no_2 - res_no_1;

	 if (res_no_diff != 1) {

	    // try to find a residue that has resno more than res_no_1 and
	    // less than res_no_2
	    //
	    for (unsigned int ii=0; ii<residues_vec.size(); ii++) {
	       int resno_this = residues_vec[ii].second->GetSeqNum();
	       // std::cout << "    res_no_this: " << resno_this << " "
	       //  	    << residue_spec_t(residues_vec[ii].second) << std::endl;
	       if (resno_this > res_no_1) {
		  if (resno_this < res_no_2) {
		     mmdb::Chain *c_this = residues_vec[ii].second->GetChain();
		     if (c_this == c_1) {
			r = true;
			break;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return r;
} 


// a pair, first is if C and N (or whatever the bonding atoms were,
// given the chem_link) are close and second is if an order switch is
// needed to make it so.
//
// return "" as first if no close link found.
//
// alt confs are ignored when finding the atoms for the potential
// link.  Perhaps they should not be (but the calling function does
// not have them).
// 
std::pair<std::string, bool>
coot::restraints_container_t::general_link_find_close_link(const std::vector<coot::chem_link> &li,
							   mmdb::Residue *r1, mmdb::Residue *r2,
							   bool order_switch_flag,
							   const coot::protein_geometry &geom) const {

   
   std::pair<std::string, bool> r("", order_switch_flag);
   std::string rs = general_link_find_close_link_inner(li, r1, r2, order_switch_flag, geom);
   if (rs != "") {
      r.first = rs;
   } else { 
      rs = general_link_find_close_link_inner(li, r2, r1, order_switch_flag, geom);
      if (rs  != "") {
	 r.first = rs;
	 r.second = 1;
      }
   }
   return r;
}

std::string
coot::restraints_container_t::general_link_find_close_link_inner(const std::vector<coot::chem_link> &li,
								 mmdb::Residue *r1, mmdb::Residue *r2,
								 bool order_switch_flag,
								 const coot::protein_geometry &geom) const {

   float dist_crit = 3.0; // Angstroms.
   bool debug = false;

   if (order_switch_flag)
      std::swap(r1, r2);

   if (debug) {
      std::cout << "DEBUG:: general_link_find_close_link_inner() 1: " << coot::residue_spec_t(r1) << std::endl;
      std::cout << "DEBUG:: general_link_find_close_link_inner() 2: " << coot::residue_spec_t(r2) << std::endl;
   }
   
   std::string r("");
   std::pair<bool,float> close = closest_approach(r1, r2);
   if (close.first) {
      if (close.second < dist_crit) {
	 for (unsigned int ilink=0; ilink<li.size(); ilink++) {
	    coot::chem_link link = li[ilink];
	    // now find link with id in dict_link_res_restraints;
	    dictionary_residue_link_restraints_t lr = geom.link(link.Id());
	    if (lr.link_id != "") { // found link with lind_id link.Id() { 
	       for (unsigned int ib=0; ib<lr.link_bond_restraint.size(); ib++) {
		  std::string atom_id_1 = lr.link_bond_restraint[ib].atom_id_1_4c();
		  std::string atom_id_2 = lr.link_bond_restraint[ib].atom_id_2_4c();
		  mmdb::Atom *at_1 = r1->GetAtom(atom_id_1.c_str());
		  mmdb::Atom *at_2 = r2->GetAtom(atom_id_2.c_str());
		  if (at_1 && at_2) {
		     clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		     double d = clipper::Coord_orth::length(p1,p2);
		     if (debug) 
			std::cout << "  dist check " << " link-bond-number: "
				  << ib << " " << coot::atom_spec_t(at_1) << " -to- "
				  <<  coot::atom_spec_t(at_2) << "   " << d << " for link type "
				  <<  lr.link_id << std::endl;
		     if (d < dist_crit) {
			r = link.Id();
			break;
		     } 
		  } else {

		     if (debug) {
			// It's fine that we don't find SG to SG for
			// L-peptide residues in an SS link
			std::cout << "WARNING:: dictionary vs model problem? "
				  << "We didn't find the bonded linking atoms :"
				  << atom_id_1 << ": :"
				  << atom_id_2 << ":"
				  << " for link-id " << link.Id() << " linking "
				  << coot::residue_spec_t(r1) << " and "
				  << coot::residue_spec_t(r2)
				  << std::endl;
			std::cout << "       :: at_1: ";
			if (at_1)
			   std::cout << coot::atom_spec_t(at_1);
			else
			   std::cout << "NULL";
			std::cout << " at_2: ";
			if (at_2)
			   std::cout << coot::atom_spec_t(at_2);
			else
			   std::cout << "NULL";
			std::cout << std::endl;
			mmdb::PPAtom residue_atoms = 0;
			int n_residue_atoms;
			r1->GetAtomTable(residue_atoms, n_residue_atoms);
			for (int i=0; i<n_residue_atoms; i++) { 
			   std::cout << "   " << r1->GetResName() << " "
				     << i << " :"  << residue_atoms[i]->name
				     << ":" << std::endl;
			}
			residue_atoms = 0;
			r2->GetAtomTable(residue_atoms, n_residue_atoms);
			for (int i=0; i<n_residue_atoms; i++) { 
			   std::cout << "   " << r2->GetResName() << " "
				     << i << " :"  << residue_atoms[i]->name
				     << ":" << std::endl;
			}
		     }
		  } 
	       }
	    }
	    if (r != "")
	       break;
	 }
      } else {
         if (debug)
	    std::cout << "FAIL close enough with closer than dist_crit " << close.second << std::endl;
      }
   }
   if (debug)
      std::cout << "DEBUG:: general_link_find_close_link_inner() return \"" << r << "\""
		<< " when called with order_switch_flag " << order_switch_flag<< std::endl;
   return r;
}

int coot::restraints_container_t::add_link_plane(std::string link_type,
						 mmdb::PResidue first, mmdb::PResidue second,
						 short int is_fixed_first_res,
						 short int is_fixed_second_res,
						 const coot::protein_geometry &geom) {

   bool debug = false;
   if (debug)
      std::cout << "DEBUG:: add_link_plane() ::::::::  for type " << link_type << " "
		<< first->GetChainID() << " " << first->GetSeqNum()
		<< " :" << first->GetInsCode() << ":"
		<< " -> " << second->GetChainID() << " " << second->GetSeqNum()
		<< " :" << second->GetInsCode() << ":" << std::endl;

   int n_plane = 0;

   mmdb::PPAtom first_sel = 0;
   mmdb::PPAtom second_sel = 0;
   mmdb::PPAtom atom_sel;  // gets assigned to either first or second
   int n_first_res_atoms, n_second_res_atoms;
   int link_res_n_atoms; // gets assigned to either n_first_res_atoms, n_second_res_atoms
   mmdb::PResidue res; 

   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   // 20100424 there can be any number of plane restraints made here (not just one for each
   // link_plane_restraint) because there can be any number of alt confs passed in the
   // residues.  So we need to find all the alt conf in the 2 residues passed and make a
   // mapping between the alt conf and the plane restraint's vector of atom indices.
   //
   std::map<std::string, std::vector<int> > atom_indices_map;
   for (int iat=0; iat<n_first_res_atoms; iat++) {
      std::string alt_loc(first_sel[iat]->altLoc);
      std::map<std::string, std::vector<int> >::const_iterator it = atom_indices_map.find(alt_loc);
      if (it == atom_indices_map.end()) {
	 std::vector<int> v;
	 atom_indices_map[alt_loc] = v;
      }
   }
   for (int iat=0; iat<n_second_res_atoms; iat++) {
      std::string alt_loc(second_sel[iat]->altLoc);
      std::map<std::string, std::vector<int> >::const_iterator it = atom_indices_map.find(alt_loc);
      if (it == atom_indices_map.end()) {
	 std::vector<int> v;
	 atom_indices_map[alt_loc] = v;
      }
   }

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically "TRANS"
	 for (unsigned int ip=0; ip<geom.link(i).link_plane_restraint.size(); ip++) {
	    std::vector<bool> fixed_flags(geom.link(i).link_plane_restraint[ip].n_atoms(), false);

	    std::map<std::string, std::vector<int> >::iterator it;
	    for (it = atom_indices_map.begin(); it != atom_indices_map.end(); it++)
	       it->second.clear();

	    for (unsigned int irest_at=0; irest_at<geom.link(i).link_plane_restraint[ip].n_atoms(); irest_at++) {
	       if (geom.link(i).link_plane_restraint[ip].atom_comp_ids[irest_at] == 1) {
		  // std::cout << "comp_id: 1 " << std::endl;
		  link_res_n_atoms = n_first_res_atoms;
		  atom_sel = first_sel;
		  res = first;
		  fixed_flags[irest_at] = is_fixed_first_res;
	       } else {
		  // std::cout << "comp_id: 2 " << std::endl;
		  link_res_n_atoms = n_second_res_atoms;
		  atom_sel = second_sel;
		  res = second; 
		  fixed_flags[irest_at] = is_fixed_second_res;
	       }
	       for (int iat=0; iat<link_res_n_atoms; iat++) { 
		  std::string pdb_atom_name(atom_sel[iat]->name);
		  if (geom.link(i).link_plane_restraint[ip].atom_id(irest_at) == pdb_atom_name) {
		     if (debug)
			std::cout << "     pushing back to :" << atom_sel[iat]->altLoc << ": vector "
				  << res->GetChainID() << " "
				  << res->seqNum << " :"
				  << atom_sel[iat]->name << ": :"
				  << atom_sel[iat]->altLoc << ":" << std::endl;


		     // Too slow for ribosomes
		     // atom_indices_map[atom_sel[iat]->altLoc].push_back(get_asc_index(atom_sel[iat]->name,
		     // atom_sel[iat]->altLoc,
		     // res->seqNum,
		     // res->GetInsCode(),
		     // res->GetChainID()));

		     std::string key(atom_sel[iat]->altLoc);
		     int idx_t_2 = -1;
		     atom_sel[iat]->GetUDData(udd_atom_index_handle, idx_t_2);
		     atom_indices_map[key].push_back(idx_t_2);

		  }
	       }
	    }


	    for (it = atom_indices_map.begin(); it != atom_indices_map.end(); ++it) {
	       
	       if (it->second.size() > 3) {
		  
		  if (debug) {
		     std::cout << "DEBUG:: add_link_plane() atom indices: for plane restraint for"
			       << " atoms with alt-conf \"" << it->first << "\"" << std::endl;
		     for (unsigned int indx=0; indx<it->second.size(); indx++) {
			std::cout << "           " << it->second[indx] << std::endl;
		     }
		     for (unsigned int ind=0; ind<it->second.size(); ind++) {
			std::cout << ind << " " << atom[it->second[ind]]->GetChainID() << " "
				  << atom[it->second[ind]]->GetSeqNum() << " :"
				  << atom[it->second[ind]]->name << ": :"
				  << atom[it->second[ind]]->altLoc << ":\n";
		     }
		     std::cout << "DEBUG:: add_link_plane() with pos indexes ";
		     for (unsigned int ipos=0; ipos<it->second.size(); ipos++)
			std::cout << " " << it->second[ipos];
		     std::cout << "\n";
		  }
		  
		  std::vector<bool> other_fixed_flags = make_fixed_flags(it->second);

		  if (debug) {
		     for (unsigned int ii=0; ii<other_fixed_flags.size(); ii++) {
			std::cout << "other fixed flags " << other_fixed_flags[ii] << std::endl;
		     }
		  }
		  for (unsigned int ii=0; ii<other_fixed_flags.size(); ii++)
		     if (other_fixed_flags[ii])
			fixed_flags[ii] = true;

		  // position and sigma
		  std::vector<std::pair<int, double> > pos_sigma;
		  for (unsigned int ii=0; ii<it->second.size(); ii++) { 
		     double sigma_esd = geom.link(i).link_plane_restraint[ip].dist_esd();
		     pos_sigma.push_back(std::pair<int, double> (it->second[ii], sigma_esd));
		  }

		  add_plane(pos_sigma, fixed_flags);
		  n_plane++;
	       }
	    }
	 }
      }
   }

   return n_plane;
}
