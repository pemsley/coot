/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013, 2015, 2016 by Medical Research Council
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

// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL

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

   bool debug = false; // Are your residues in the same chain?  If not filter() will not bond them.

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

	 // 20180911 I now no longer want to evaluate closest approach here.
	 //
	 // std::pair<bool, float> d = closest_approach(res_f, res_s);
	 // Linking should be resolved by find_link_type_complicado(), not
	 // here by distance between residues.

	 std::pair<std::string, bool> l = find_link_type_complicado(res_f, res_s, geom);

	 // too verbose?
	 if (false)
	    std::cout << "   INFO:: find_link_type_complicado() for: "
		      << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
		      << " returns link_type -> \"" << l.first << "\"" << std::endl;

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

                                    std::cout << "----------------------- adding link torsion! "
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

   bool debug = false;
   // this determines the link type

   auto tp_0 = std::chrono::high_resolution_clock::now();
   bonded_pair_container_t bonded_residue_pairs = bonded_residues_from_res_vec(geom);
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
      r = geom.find_glycosidic_linkage_type(first, second); // 20150714 FIXME
   else 
      r = geom.find_glycosidic_linkage_type(first, second, mol);
   return r;
}

bool
coot::restraints_container_t::link_infos_are_glycosidic_p(const std::vector<std::pair<coot::chem_link, bool> > &link_infos) const {

   bool is_sweet = 0;
   for (unsigned int i=0; i<link_infos.size(); i++) {
      std::string id = link_infos[i].first.Id();
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

   bool by_distance = false; // 20150714 new style looks for LINKs and
			     // SSBONDs in the molecule before
			     // assigning pyranose or disulphides etc.

   bool debug = false; // TMI.

   std::string link_type = "";
   bool order_switch_flag = false;
   std::string comp_id_1 = first->GetResName();
   std::string comp_id_2 = second->GetResName();

   try {
      std::string group_1 = geom.get_group(first);
      std::string group_2 = geom.get_group(second);
      if (group_1 == "DNA") group_1 = "DNA/RNA";
      if (group_1 == "RNA") group_1 = "DNA/RNA";
      if (group_2 == "DNA") group_2 = "DNA/RNA";
      if (group_2 == "RNA") group_2 = "DNA/RNA";
      if (debug)
	 std::cout << "====== DEBUG:: find_link_type_complicado() from "
		   << first->GetChainID() << " " << first->GetSeqNum() << " " << first->GetResName()
		   << " <--> "
		   << second->GetChainID() << " " << second->GetSeqNum() << " " << second->GetResName()
		   << "  comp-id-1 \"" << comp_id_1 << "\" group-1 \"" << group_1 << "\" comp-id-2 \""
                   << comp_id_2 << "\" group-2 \"" << group_2 << "\"" << std::endl;

      try {
	 std::vector<std::pair<coot::chem_link, bool> > link_infos =
	    geom.matching_chem_link(comp_id_1, group_1, comp_id_2, group_2);

	 if (false) { // what are the possible matching chem links?
                     // This is useful when trying to find why a cif file doesn't link the
                     // atoms that you want it to link.

	    std::cout << "  ====== DEBUG:: find_link_type_complicado() found " << link_infos.size()
                      << " link infos for residue compids/groups "
		      << comp_id_1 << " " << group_1 << " and "
		      << comp_id_2 << " " << group_2 << std::endl;
	    std::cout << "   ===== DEBUG:: find_link_type_complicado() first  " << coot::residue_spec_t(first)
		      << std::endl;
	    std::cout << "   ===== DEBUG:: find_link_type_complicado() second " << coot::residue_spec_t(second)
		      << std::endl;
	    std::cout << "   ===== DEBUG:: find_link_type_complicado() " << link_infos.size()
		      << " possible links: (link_infos):\n";
	    for (unsigned int il=0; il<link_infos.size(); il++)
	       std::cout << "   ===== DEBUG:: find_link_type_complicado() link_info["
			 << il << "]: " << link_infos[il].first << "    order swap: "
			 << link_infos[il].second << " "
			 << std::endl;
	 }

	 // Now, if link is a TRANS (default-peptide-link), then
	 // make sure that the C and N (or N and C) atoms of the
	 // first and second residue are within dist_crit (2.0A (is this right?)) of
	 // each other.  If not, then we should fail to find a link
	 // between these 2 residues.
	 // 
	 for (unsigned int ilink=0; ilink<link_infos.size(); ilink++) {

	    if (debug)
	       std::cout << "   ===== DEBUG:: find_link_type_complicado() LINKS:: for-loop testing ilink " << ilink << " "
			 << link_infos[ilink].first << " bool: " << link_infos[ilink].second
			 << std::endl;

	    if (link_infos[ilink].first.is_peptide_link_p()) {

	       // ===================== is a peptide ========================================
	       
	       // note to self: use an enum for the first, don't use magic numbers
	       std::pair<int, bool> peptide_info = peptide_C_and_N_are_in_order_p(first, second);

	       if (false)
		  std::cout << "   peptide_C_and_N_are_in_order_p returns " <<  peptide_info.first
			    << " and order-switch: " << peptide_info.second << std::endl;


	       //
	       if (peptide_info.first == IS_PEPTIDE) {
		  order_switch_flag = peptide_info.second;

		  mmdb::Residue *r_1 = first;
		  mmdb::Residue *r_2 = second;

		  if (order_switch_flag)
		     std::swap(r_1, r_2);

		  if (! have_intermediate_residue_by_seqnum(r_1, r_2)) {

		     std::string rn2 = r_2->name;
		     if (rn2 == "PRO") {
			link_type = "PTRANS";
		     } else {

			link_type = "TRANS"; // 200100415 for now, we force all peptide links to be
			                     // TRANS.  (We don't yet (as of today) know if this link
				             // was CIS or TRANS). TRANS has 5-atom (plane3) plane
				             // restraints, CIS does not.
		     }
		  }
	       } else {

		  if (peptide_info.first == UNKNOWN) {

		     // maybe it was a link between residues in different molecules
		     // i.e. a moving residue and a non moving one (things are in a mess)
		     // or there was an insertion code.
		     // so let's try by distance and residue details.

		     peptide_info = peptide_C_and_N_are_close_p(first, second);

                     if (false) {
		        std::cout << "   ===== DEBUG:: find_link_type_complicado() unknown, so moving-to-fixed? for "
				  << residue_spec_t(first) << " " << residue_spec_t(second) << " "
                                  << peptide_info.first << " " << peptide_info.second
                                  << std::endl;
		     }
		  }
	       }

	       if (peptide_info.first != IS_PEPTIDE) {

		  // is not a peptide link

		  std::vector<std::pair<coot::chem_link, bool> > link_infos_non_peptide =
		     geom.matching_chem_link_non_peptide(comp_id_1, group_1, comp_id_2, group_2);

		  if (debug)
		     for (unsigned int il=0; il<link_infos_non_peptide.size(); il++)
			std::cout << "  ===== DEBUG:: find_link_type_complicado() non-peptide link: "
				  << link_infos_non_peptide[il].first.Id() << std::endl;

		  // 20100330 eh?  is something missing here?  What
		  // shall we do with link_infos_non_peptide?  did I
		  // mean that, now the peptide test has failed,
		  // perhaps we should distance check potential other
		  // (non-peptide) links?
		  // 
		  // Yes, I think I did.
		  // 
		  // returns "" if no link found
		  // second is the order switch flag;
		  // 
		  std::pair<std::string, bool> non_peptide_close_link_info =
		     general_link_find_close_link(link_infos_non_peptide, first, second, order_switch_flag, geom);
                  if (debug)
		     std::cout << "---- DEBUG:: find_link_type_complicado() general_link_find_close_link() "
			       << non_peptide_close_link_info.first << " "
			       << non_peptide_close_link_info.second << std::endl;

		  if (!non_peptide_close_link_info.first.empty()) {

                     bool is_linked = false;
                     mmdb::Model *model_p = first->GetModel();
                     if (model_p) {
                        int n_links = model_p->GetNumberOfLinks();
                        // std::cout << "-------------------------------------- n_links " << n_links << std::endl;
                        if (n_links > 0) {
                           mmdb::Residue *r_1 = first;
                           mmdb::Residue *r_2 = second;
                           residue_spec_t SS_residue_1_spec(r_1);
                           residue_spec_t SS_residue_2_spec(r_2);
                           for (int i_link=1; i_link<=n_links; i_link++) {
                              mmdb::PLink link = model_p->GetLink(i_link);
                              std::pair<coot::atom_spec_t, coot::atom_spec_t> link_atom_specs = coot::link_atoms(link, model_p);
                              residue_spec_t link_residue_spec_1(link_atom_specs.first);
                              residue_spec_t link_residue_spec_2(link_atom_specs.second);
                              if (link_residue_spec_1  == SS_residue_1_spec) is_linked = true;
                              if (link_residue_spec_2  == SS_residue_1_spec) is_linked = true;
                              if (link_residue_spec_1  == SS_residue_2_spec) is_linked = true;
                              if (link_residue_spec_2  == SS_residue_2_spec) is_linked = true;
                              if (is_linked) break;
                           }
                        }
                     }

                     if (! is_linked) {
                        link_type = non_peptide_close_link_info.first;
                        order_switch_flag = non_peptide_close_link_info.second;
                     }
		  }
	       }
	       
	    } else {

	       // ===================== is not a peptide ========================================

	       if (debug)
		  std::cout << "   ====== DEBUG:: find_link_type_complicado() this "
			    << link_infos[ilink].first
			    << " is not a peptide link... "
			    << std::endl;
	       //
	       order_switch_flag = link_infos[ilink].second;
	       bool use_links_in_molecule = true; // 20150714

	       // we should not be testing link_infos inside a loop that
	       // selects a particular link info. Hmmmm.
	       // The test should be link_info.is_glycosidic(link_infos[ilink])
	       //
	       if (debug)
		  std::cout << "   ===== DEBUG:: find_link_type_complicado() link_infos are glycosidic: "
			    << link_infos_are_glycosidic_p(link_infos) << std::endl;
	    
	       if (link_infos_are_glycosidic_p(link_infos)) {

		  // ===================== is  glycosidic ========================================

		  link_type = find_glycosidic_linkage_type(first, second, geom, use_links_in_molecule);
		  if (debug)
		     std::cout << " ===== DEBUG:: find_link_type_complicado() find_glycosidic_linkage_type() returns \""
			       << link_type << "\"" << std::endl;
		  if (link_type == "") {
		     link_type = find_glycosidic_linkage_type(second, first, geom, use_links_in_molecule);
		     if (link_type != "")
			order_switch_flag = 1;
		  }

		  // 20150714 This bit of code here is indicative of sloppy thinking previously.
		  // We have a list of matching link_infos.  We should know before we run the
		  // for(ilink) loop to run find_glycosidic_linkage_type().  That loops over
		  // links independently of link_infos - it should be done outside the for(ilink)
		  // loop.
		  // 
		  // For now, we will break out of the for(ilink) loop - we don't need to
		  // rerun find_glycosidic_linkage_type() for every link in link_infos.
		  //
		  if (! link_type.empty())
		      break;

	       } else {

		  // ===================== is not glycosidic ========================================

		  if (debug)
		     std::cout << "    ===== DEBUG::find_link_type_complicado() not a glycosidic_linkage"
			       << std::endl;

		  std::vector<std::pair<coot::chem_link, bool> > link_infos_non_peptide;

		  if (! use_links_in_molecule) { 
		     link_infos_non_peptide = 
			geom.matching_chem_link_non_peptide(comp_id_1, group_1, comp_id_2, group_2);
		  } else {
		     link_infos_non_peptide = 
			geom.matching_chem_link_non_peptide(comp_id_1, group_1, comp_id_2, group_2, mol);
		  }

		  if (false)
		     for (unsigned int il=0; il<link_infos_non_peptide.size(); il++)
			std::cout << "   ===== DEBUG:: find_link_type_complicado() pre-calling general_link_find_close_link() "
				  << il << " of " << link_infos_non_peptide.size()
				  << " possible non-peptide link: " << link_infos_non_peptide[il].first.Id()
				  << std::endl;

		  // returns "" if no link found
		  // second is the order switch flag;
		  //
		  // std::cout << "debug:: calling general_link_find_close_link() here" << std::endl;
		  std::pair<std::string, bool> non_peptide_close_link_info = 
		     general_link_find_close_link(link_infos_non_peptide, first, second, order_switch_flag, geom);
		  if (non_peptide_close_link_info.first != "") {

                     // Filter out SS bonds that are to SGs that have links - The Joosten Filter
                     if (non_peptide_close_link_info.first != "SS") {
                        link_type = non_peptide_close_link_info.first;
                        order_switch_flag = non_peptide_close_link_info.second;
                     } else {
                        // SS Bond. Do we have a link in the header?

                        bool is_linked = false;
                        mmdb::Model *model_p = first->GetModel();
                        if (model_p) {
                           int n_links = model_p->GetNumberOfLinks();
                           // std::cout << "-------------------------------------- n_links "
                           //           << n_links << std::endl;
                           if (n_links > 0) {
                              mmdb::Residue *r_1 = first;
                              mmdb::Residue *r_2 = second;
                              residue_spec_t SS_residue_1_spec(r_1);
                              residue_spec_t SS_residue_2_spec(r_2);
                              for (int i_link=1; i_link<=n_links; i_link++) {
                                 mmdb::PLink link = model_p->GetLink(i_link);
                                 std::pair<atom_spec_t, atom_spec_t> link_atom_specs = link_atoms(link, model_p);
                                 residue_spec_t link_residue_spec_1(link_atom_specs.first);
                                 residue_spec_t link_residue_spec_2(link_atom_specs.second);
                                 if (link_residue_spec_1 == SS_residue_1_spec) is_linked = true;
                                 if (link_residue_spec_2 == SS_residue_1_spec) is_linked = true;
                                 if (link_residue_spec_1 == SS_residue_2_spec) is_linked = true;
                                 if (link_residue_spec_2 == SS_residue_2_spec) is_linked = true;
                                 if (is_linked) break;
                              }
                           }
                        }
                        // std::cout << "debug:: here with is_linked " << is_linked << " link_type \"" << link_type << "\"" << std::endl;

                        if (! is_linked) {
                           link_type = non_peptide_close_link_info.first;
                           order_switch_flag = non_peptide_close_link_info.second;
                        }
                     }
		  }
	       }
	    }
	 }
      }
      catch (const std::runtime_error &mess_in) {
	 if (debug) {
	    // didn't find a chem_link for this pair, that's OK sometimes.
	    std::cout << "CAUGHT exception in find_link_type_complicado(): "
		      << mess_in.what() << std::endl;
	    geom.print_chem_links();
	 }
      }
   }
   catch (const std::runtime_error &mess) {
      // Failed to get group.  We don't want to hear about not getting
      // the group of thousands of HOHs.
      // std::cout << mess.what() << std::endl;
   }

   if (debug)
      std::cout << "   ===== DEBUG:: find_link_type_complicado() given "
		<< coot::residue_spec_t(first) << " and "
		<< coot::residue_spec_t(second)
		<< " find_link_type_complicado returns link type :"
		<< link_type << ": and order-switch flag " << order_switch_flag << std::endl;
   return std::pair<std::string, bool> (link_type, order_switch_flag);
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
coot::restraints_container_t::general_link_find_close_link(const std::vector<std::pair<coot::chem_link, bool> > &li,
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
coot::restraints_container_t::general_link_find_close_link_inner(const std::vector<std::pair<coot::chem_link, bool> > &li,
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
	    coot::chem_link link = li[ilink].first;
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


	    for (it = atom_indices_map.begin(); it != atom_indices_map.end(); it++) {
	       
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

		  if (debug) { // debug
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

#endif // HAVE_GSL

