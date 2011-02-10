
// c.f. with O-function lsq_improve

#include <map>
#include <stdexcept>
#include <algorithm>

#include <string.h> // strncpy

#include "lsq-improve.hh"

coot::lsq_improve::lsq_improve(CMMDBManager *mol_ref, CMMDBManager *mol_moving) {
   mol = new CMMDBManager;
   crit_close = 6; // A
   n_res_for_frag = 6;
   
   n_ref_atoms = CAs_to_model(mol_ref,    1);
   n_mov_atoms = CAs_to_model(mol_moving, 2);

   std::cout << "addded n_ref_atoms: " << n_ref_atoms << std::endl;
   std::cout << "addded n_mov_atoms: " << n_mov_atoms << std::endl;

   if (! n_ref_atoms) {
      std::cout << "no CA atoms from ref mol " << std::endl;
   } else { 
      if (! n_mov_atoms) { 
	 std::cout << "no CA atoms from moving mol " << std::endl;
      } else { 
	 sel_hnd_1 = mol->NewSelection();
	 sel_hnd_2 = mol->NewSelection();

	 mol->SelectAtoms (sel_hnd_1, 1, "*",
			   ANY_RES, // starting resno, an int
			   "*", // any insertion code
			   ANY_RES, // ending resno
			   "*", // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*"  // alt loc.
			   );
   
	 mol->SelectAtoms (sel_hnd_2, 2, "*",
			   ANY_RES, // starting resno, an int
			   "*", // any insertion code
			   ANY_RES, // ending resno
			   "*", // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*"  // alt loc.
			   );

	 
	 PPCAtom atom_sel_1 = 0;
	 PPCAtom atom_sel_2 = 0;
	 int n_sel_1;
	 int n_sel_2;
	 mol->GetSelIndex(sel_hnd_1, atom_sel_1, n_sel_1);
	 mol->GetSelIndex(sel_hnd_2, atom_sel_2, n_sel_2);
	 
	 for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	    CModel *model_p = mol->GetModel(imod);
	    CChain *chain_p;
	    int nchains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       CResidue *residue_p;
	       CAtom *at;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  for (int iat=0; iat<n_atoms; iat++) {
		     at = residue_p->GetAtom(iat);
		  }
	       }
	    }
	 }
      }
   }
}

// return the number of Atoms
int 
coot::lsq_improve::CAs_to_model(CMMDBManager *mol_in, int model_number) {

   int n_atoms = 0; // increment and return
   if (mol_in) { 
      int imod = 1;
      CModel *model_p = mol_in->GetModel(imod);
      if (!model_p) {
	 std::cout << "Oops no MODEL 1 in input molecule for synthmol model-no " << model_number
		   << std::endl;
      } else {

	 CModel *model_new = new CModel;
	 mol->AddModel(model_new);
	 int nchains = model_p->GetNumberOfChains();
	 std::cout << "refmol has " << nchains << " chains " << std::endl;
	 for (int ichain=0; ichain<nchains; ichain++) {
	    CChain *chain_p = model_p->GetChain(ichain);
	    CChain *chain_new = new CChain(model_new, chain_p->GetChainID());
	    model_new->AddChain(chain_new);
	    int nres = chain_p->GetNumberOfResidues();
	    std::cout << "   refmol chain " << chain_p->GetChainID()
		      << " has " << nres << " residues " << std::endl;
	    for (int ires=0; ires<nres; ires++) { 
	       CResidue *residue_p = chain_p->GetResidue(ires);
	       CAtom *at = residue_p->GetAtom(" CA ");
	       if (at) {
		  CResidue *residue_new = new CResidue(chain_new);
		  chain_new->AddResidue(residue_new);
		  residue_new->SetResName(residue_p->GetResName());
		  residue_new->seqNum = residue_p->GetSeqNum();
		  strncpy(residue_new->insCode, residue_p->GetInsCode(), 3);
		  CAtom *atom_new = new CAtom(residue_new);
		  residue_new->AddAtom(atom_new);
		  atom_new->Copy(at);
		  n_atoms++;
	       }
	    }
	 }
	 mol->FinishStructEdit();

	 std::cout << "model " << model_number << " has " << model_new->GetNumberOfAtoms(1)
		   << " atoms "<< std::endl;
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
	 int n_rounds_max = 1;

	 while (! found_match_before) {

	    std::cout << "DEBUG:: improve round " << i_round << std::endl;
	    std::vector<coot::lsq_range_match_info_t> new_matches = get_new_matches();

	    for (unsigned int imatch=0; imatch<previous_matches.size(); imatch++) { 
	       found_match_before = same_matches(new_matches, previous_matches[imatch]);
	       if (found_match_before)
		  break;
	    }
	    if (! found_match_before) { 
	       previous_matches.push_back(new_matches);
	       apply_matches(new_matches);
	    }

	    // just for sanity sake
	    i_round++;
	    if (i_round >= n_rounds_max)
	       break;
	 }
      }
   }
}


bool
coot::lsq_improve::same_matches(const std::vector<coot::lsq_range_match_info_t> &ranges_1,
				const std::vector<coot::lsq_range_match_info_t> &ranges_2) const {
   return 0;
} 


std::vector<coot::lsq_range_match_info_t>
coot::lsq_improve::get_new_matches() const {

   std::vector<coot::lsq_range_match_info_t> r;

   PPCAtom atom_sel_1;
   PPCAtom atom_sel_2;
   int n_sel_1;
   int n_sel_2;
   mol->GetSelIndex(sel_hnd_1, atom_sel_1, n_sel_1);
   mol->GetSelIndex(sel_hnd_2, atom_sel_2, n_sel_2);

   int ncontacts = 0;
   PSContact contact = NULL;

   mol->SeekContacts(atom_sel_1, n_sel_1,
		     atom_sel_2, n_sel_2,
		     0.0, crit_close,
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
	 CAtom *atom_ref = atom_sel_1[contact[icon].id1];
	 CAtom *atom_mov = atom_sel_2[contact[icon].id2];
	 coot::residue_spec_t r_1(atom_ref);
	 coot::residue_spec_t r_2(atom_mov);
	 contact_residues[r_1].push_back(r_2);
	 // include an "already matched" flag with the same key
      }
   }
   delete [] contact; // we are done with contact now.

   if (0) { 
      std::cout << "residue specs after contacts " << std::endl;
      std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> >::iterator iter;
      for (iter=contact_residues.begin(); iter!=contact_residues.end(); iter++) {
	 std::cout << " " << iter->first << std::endl;
	 for (unsigned int ispec=0; ispec<iter->second.size(); ispec++) { 
	    std::cout << "      " << iter->second[ispec] << std::endl;
	 }
      }
   }
   
   // the heart of the function
   // 
   r = get_new_matches(contact_residues);
   std::cout << "--------------- found " << r.size() << " matches ------------ " << std::endl;


   
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

		  if (0)
		     std::cout << "   depth " << depth << " ref-res:   " << iter_running->first
			       << " considering match to mov spec "
			       << imov << "/" << iter_running->second.size() << "  "
			       << iter_running->second[imov] << std::endl;

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
	       // std::cout << ".... updating mov_running from " << mov_running << " to "
	       // << mov_running.next() << std::endl;
	       mov_running = mov_running.next();
	    } else {
	       // we will fail out of the while next time it is tested
	       // std::cout << "Found a frag of size " << contiguous_frag_spec.size() << std::endl;
	       if (contiguous_frag_spec.size() >= n_res_for_frag) {
		  for (unsigned int ipair=0; ipair<contiguous_frag_spec.size(); ipair++) { 
		     residue_done[contiguous_frag_spec[ipair].first] = 1;
		  }
		  coot::lsq_range_match_info_t range(contiguous_frag_spec.begin()->first.resno,
						     contiguous_frag_spec.back().first.resno,
						     contiguous_frag_spec.begin()->first.chain,
						     contiguous_frag_spec.begin()->second.resno,
						     contiguous_frag_spec.back().second.resno,
						     contiguous_frag_spec.begin()->second.chain,
						     COOT_LSQ_CA);
						     
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
   std::pair<short int, clipper::RTop_orth> mat = coot::util::get_lsq_matrix(mol, mol, matches, 1);
}


// can throw an exception
//
clipper::RTop_orth
coot::lsq_improve::rtop_of_moving() const {

   std::vector<coot::lsq_range_match_info_t> matches = get_new_matches();
   clipper::RTop_orth rtop = rtop_of_moving(matches);
   return rtop;
} 

// can throw an exception
//
clipper::RTop_orth
coot::lsq_improve::rtop_of_moving(const std::vector<coot::lsq_range_match_info_t> &matches) const {

   clipper::RTop_orth rtop;
   if (! n_ref_atoms) {
      std::string message = "no CA atoms from ref mol ";
      throw(std::runtime_error(message));
   } else { 
      if (! n_mov_atoms) { 
	 std::string message = "no CA atoms from moving mol ";
	 throw(std::runtime_error(message));
      } else {
	 // do something
      }
   }
   return rtop;
} 
