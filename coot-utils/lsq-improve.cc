/* coot-utils/lsq-improve.cc
 * 
 * Copyright 2011 by The University of Oxford
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

// c.f. with O-function lsq_improve

#include <map>
#include <stdexcept>
#include <algorithm>

#include <string.h> // strncpy

#include "lsq-improve.hh"

coot::lsq_improve::lsq_improve(mmdb::Manager *mol_ref, const std::string &ref_selection_string,
			       mmdb::Manager *mol_moving, const std::string &moving_selection_string) {

   mol = new mmdb::Manager;
   mol_initial_copy = NULL; // reset later hopefully.
   crit_close = 6; // A
   n_res_for_frag = 6;
   n_rounds_max = 10;
   
   n_ref_atoms = CAs_to_model(mol_ref,    1);
   n_mov_atoms = CAs_to_model(mol_moving, 2);

   if (! n_ref_atoms) {
      std::cout << "no CA atoms from ref mol " << std::endl;
   } else { 
      if (! n_mov_atoms) { 
	 std::cout << "no CA atoms from moving mol " << std::endl;
      } else {
	 // delete on deconstruction.
	 sel_hnd_1 = mol->NewSelection();
	 sel_hnd_2 = mol->NewSelection();

	 // now apply the passed selections
	 // 
	 mol->Select(sel_hnd_1, mmdb::STYPE_ATOM,    ref_selection_string.c_str(), mmdb::SKEY_OR);
	 mol->Select(sel_hnd_2, mmdb::STYPE_ATOM, moving_selection_string.c_str(), mmdb::SKEY_OR);

	 mmdb::PPAtom atom_sel_1=NULL;
	 mmdb::PPAtom atom_sel_2=NULL;
	 int n_sel_1;
	 int n_sel_2;
	 mol->SelectAtoms (sel_hnd_1, 1, "*",
			   mmdb::ANY_RES, // starting resno, an int
			   "*", // any insertion code
			   mmdb::ANY_RES, // ending resno
			   "*", // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*",  // alt loc.
			   mmdb::SKEY_AND
			   );
   
	 mol->SelectAtoms (sel_hnd_2, 2, "*",
			   mmdb::ANY_RES, // starting resno, an int
			   "*", // any insertion code
			   mmdb::ANY_RES, // ending resno
			   "*", // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*",  // alt loc.
			   mmdb::SKEY_AND
			   );
	 
	 mol->GetSelIndex(sel_hnd_1, atom_sel_1, n_sel_1);
	 mol->GetSelIndex(sel_hnd_2, atom_sel_2, n_sel_2);

	 // throw an exception here if n_sel_1 or n_sel_1 are 2 or less.

	 mol_initial_copy = new mmdb::Manager;
	 mol_initial_copy->Copy(mol, mmdb::MMDBFCM_All);
      }
   }
}

// return the number of Atoms
int 
coot::lsq_improve::CAs_to_model(mmdb::Manager *mol_in, int model_number) {

   int n_atoms = 0; // increment and return
   if (mol_in) { 
      int imod = 1;
      mmdb::Model *model_p = mol_in->GetModel(imod);
      if (!model_p) {
	 std::cout << "Oops no MODEL 1 in input molecule for synthmol model-no "
		   << model_number << std::endl;
      } else {

	 mmdb::Model *model_new = new mmdb::Model;
	 mol->AddModel(model_new);
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    mmdb::Chain *chain_new = new mmdb::Chain(model_new, chain_p->GetChainID());
	    model_new->AddChain(chain_new);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) { 
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       mmdb::Atom *at = residue_p->GetAtom(" CA ");
	       if (at) {
		  mmdb::Residue *residue_new = new mmdb::Residue(chain_new);
		  chain_new->AddResidue(residue_new);
		  residue_new->SetResName(residue_p->GetResName());
		  residue_new->seqNum = residue_p->GetSeqNum();
		  strncpy(residue_new->insCode, residue_p->GetInsCode(), 3);
		  mmdb::Atom *atom_new = new mmdb::Atom(residue_new);
		  residue_new->AddAtom(atom_new);
		  atom_new->Copy(at);
		  n_atoms++;
	       }
	    }
	 }
	 mol->FinishStructEdit();
      }
   }
   return n_atoms;
} 

void coot::lsq_improve::improve() {

   if (! n_ref_atoms) {
      std::cout << "no CA atoms from ref mol " << std::endl;
   } else { 
      if (! n_mov_atoms) { 
	 std::cout << "no CA atoms from moving mol " << std::endl;
      } else { 
	 std::vector<std::vector<coot::lsq_range_match_info_t> > previous_matches;
	 bool found_match_before = 0;
	 int i_round = 0;

	 while (i_round <= n_rounds_max) {

	    bool summary_to_screen_flag = 0;
	    if (i_round == 0 || i_round == n_rounds_max)
	       summary_to_screen_flag = 1;
	    std::vector<coot::lsq_range_match_info_t> new_matches =
	       get_new_matches(i_round, n_rounds_max, summary_to_screen_flag);
	    previous_matches.push_back(new_matches);
	    apply_matches(new_matches);
	    
	    i_round++;
	 }
      }
   }
}



std::vector<coot::lsq_range_match_info_t>
coot::lsq_improve::get_new_matches(int round_number, int round_max,
				   bool summary_to_screen_flag) const {

   std::vector<coot::lsq_range_match_info_t> r;

   mmdb::PPAtom atom_sel_1;
   mmdb::PPAtom atom_sel_2;
   int n_sel_1;
   int n_sel_2;
   mol->GetSelIndex(sel_hnd_1, atom_sel_1, n_sel_1);
   mol->GetSelIndex(sel_hnd_2, atom_sel_2, n_sel_2);

   int ncontacts = 0;
   mmdb::Contact *contact = NULL;

   // round_number  round_max  -> multiplier
   //    0             10         1
   //    5             10         1 - (1-0.3)*round_number/round_max
   //   10             10         0.3
   mmdb::realtype min_dist = 0.3;
   mmdb::realtype multiplier = 1.0 - (1.0 - min_dist)*mmdb::realtype(round_number)/mmdb::realtype(round_max);
   mmdb::realtype max_dist = crit_close * multiplier;

   std::cout << "   round " << round_number << " multiplier  " << multiplier
	     << "  max_dist " << max_dist << std::endl;

   mol->SeekContacts(atom_sel_1, n_sel_1,
		     atom_sel_2, n_sel_2,
		     0.0, max_dist,
		     0,
		     contact, ncontacts);

   // find residues in mol2 (moving) that are close to residues of mol 1 (reference), key
   // 
   std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> > contact_residues;

//    std::cout << "n_sel_1:: " << n_sel_1 << std::endl;
//    std::cout << "n_sel_2:: " << n_sel_2 << std::endl;
//    std::cout << "SeekContacts() found ncontacts: " << ncontacts << std::endl;
   
   
   if (ncontacts) {
      for (int icon=0; icon<ncontacts; icon++) {
	 mmdb::Atom *atom_ref = atom_sel_1[contact[icon].id1];
	 mmdb::Atom *atom_mov = atom_sel_2[contact[icon].id2];
	 coot::residue_spec_t r_1(atom_ref->GetResidue());
	 coot::residue_spec_t r_2(atom_mov->GetResidue());
	 contact_residues[r_1].push_back(r_2);
	 // include an "already matched" flag with the same key
      }
   }
   delete [] contact; // we are done with contact now.

   if (0) { 
      std::cout << "residue specs after contacts " << std::endl;
      std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> >::iterator iter;
      for (iter=contact_residues.begin(); iter!=contact_residues.end(); ++iter) {
	 std::cout << " " << iter->first << std::endl;
	 for (unsigned int ispec=0; ispec<iter->second.size(); ispec++) { 
	    std::cout << "      " << iter->second[ispec] << std::endl;
	 }
      }
   }
   
   // the heart of the function
   // 
   r = get_new_matches(contact_residues);

   if (summary_to_screen_flag) { 
      std::cout << "INFO:: --- round " << round_number
		<< " found " << r.size() << " matches ---- " << std::endl;
      for (unsigned int i=0; i<r.size(); i++)
	 std::cout << "    " << r[i] << std::endl;
   }
   
   return r;
}

// core function
// 
std::vector<coot::lsq_range_match_info_t>
coot::lsq_improve::get_new_matches(const std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> > &contact_residues) const {

   std::vector<coot::lsq_range_match_info_t> r;
   
   std::map<coot::residue_spec_t, bool> residue_done; // markers for contact_residues
   std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> >::const_iterator iter;
   std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> >::const_iterator iter_running;
   std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> >::const_iterator it;
   for (iter=contact_residues.begin(); iter!=contact_residues.end(); iter++)
      residue_done[iter->first] = 0;

   
   for (iter=contact_residues.begin(); iter!=contact_residues.end(); iter++) {
      if (! residue_done[iter->first]) {

	 // std::vector<coot::lsq_range_match_info_t> contiguous_frag_spec;
	 std::vector<std::pair<coot::residue_spec_t, coot::residue_spec_t> > contiguous_frag_spec;

	 // std::cout << "==== outer loop ref res: " << iter->first << std::endl;
	       
	 bool found_next_pair = 1; // init
	 iter_running = iter;
	 int depth = 0;
	 coot::residue_spec_t mov_running;
	 
	 while (found_next_pair) {

	    depth++;
	    
	    // Is there a next (ref) which has a next (moving)? (There is
	    // a vector of moving that we need to run through)
	    
	    it = contact_residues.find(iter_running->first.next());
	    if (it == contact_residues.end()) {
	       found_next_pair = 0;
	       // std::cout << "fail to find next ref res ("<< iter_running->first.next()
	       // << ")" << std::endl;
	    } else {

	       bool ifound = 0; // init
	       for (unsigned int imov=0; imov<iter_running->second.size(); imov++) {

		  if (depth == 1 || iter_running->second[imov] == mov_running) { 
		     mov_running = iter_running->second[imov];
		     // std::cout << "     pushing pair " << iter_running->first << "  "
		     // << iter_running->second[imov] << std::endl;
		     std::pair<coot::residue_spec_t, coot::residue_spec_t> pair(iter_running->first, iter_running->second[imov]);
		     contiguous_frag_spec.push_back(pair);
		     ifound = 1;
		     break;
		  }
	       }
	       if (! ifound)
		  found_next_pair = 0;
	    }

	    if (found_next_pair) {
	       iter_running++;
	       mov_running = mov_running.next();
	    } else {
	       // we will fail out of the while next time it is tested
	       if (int(contiguous_frag_spec.size()) >= n_res_for_frag) {
		  for (unsigned int ipair=0; ipair<contiguous_frag_spec.size(); ipair++) { 
		     residue_done[contiguous_frag_spec[ipair].first] = 1;
		  }
		  coot::lsq_range_match_info_t range(contiguous_frag_spec.begin()->first.res_no,
						     contiguous_frag_spec.back().first.res_no,
						     contiguous_frag_spec.begin()->first.chain_id,
						     contiguous_frag_spec.begin()->second.res_no,
						     contiguous_frag_spec.back().second.res_no,
						     contiguous_frag_spec.begin()->second.chain_id,
						     COOT_LSQ_CA);
		  range.set_model_number_reference(1);
		  range.set_model_number_matcher(2);
		  r.push_back(range);
	       } 
	    } 
	 } 
      }
   }
   return r;
}


// move the moving model (with model number 2) in mol.
void
coot::lsq_improve::apply_matches(const std::vector<coot::lsq_range_match_info_t> &matches) {

   std::pair<short int, clipper::RTop_orth> mat = 
      coot::util::get_lsq_matrix(mol, mol, matches, 1, 0);
   if (mat.first) { 
      coot::util::transform_selection(mol, sel_hnd_2, mat.second);
      if (0) { // debug
	 std::cout << "apply_matches() moving model by matrix:\n"
		   << mat.second.format()
		   << std::endl;
	 std::string file_name = "applied.pdb";
	 mol->WritePDBASCII(file_name.c_str());
      } 
   } else {
      std::cout << "OOOpps!  bad matrix in apply_matches() "
		<< " - this should not happen" << std::endl;
   }

}


// can throw an exception
//
clipper::RTop_orth
coot::lsq_improve::rtop_of_moving() const {

   std::vector<coot::lsq_range_match_info_t> matches =
      get_new_matches(n_rounds_max, n_rounds_max);
   clipper::RTop_orth rtop = rtop_of_moving(matches);
   return rtop;
} 

// can throw an exception
//
clipper::RTop_orth
coot::lsq_improve::rtop_of_moving(const std::vector<coot::lsq_range_match_info_t> &matches) const {

   // std::cout << "rtop_of_moving() using ranges: " << std::endl;
   // for (unsigned int i=0; i<matches.size(); i++)
   // std::cout << "    " << matches[i] << std::endl;
   
   clipper::RTop_orth rtop;
   if (! n_ref_atoms) {
      std::string message = "no CA atoms from ref mol ";
      throw(std::runtime_error(message));
   } else { 
      if (! n_mov_atoms) { 
	 std::string message = "no CA atoms from moving mol ";
	 throw(std::runtime_error(message));
      } else {
	 if (!mol_initial_copy) {
	    std::string message = "Null copy of initial! ";
	    throw(std::runtime_error(message));
	 } else { 
	    std::pair<short int, clipper::RTop_orth> mat = 
	       coot::util::get_lsq_matrix(mol, mol_initial_copy, matches, 1, 0);
	    if (! mat.first) {
	       std::string message = "Bad matrix ";
	       throw(std::runtime_error(message));
	    } else {
	       rtop = mat.second;
	    }
	 }
      }
   }
   return rtop;
} 
