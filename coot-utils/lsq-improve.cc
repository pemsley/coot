
// c.f. with O-function lsq_improve

#include <map>
#include <stdexcept>

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
	 std::cout << ".... in constructor selected " << n_sel_1 << " atoms " << std::endl;
	 std::cout << ".... in constructor selected " << n_sel_2 << " atoms " << std::endl;

	 
	 std::cout << "... n_models " << mol->GetNumberOfModels() << std::endl;
	 for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	    CModel *model_p = mol->GetModel(imod);
	    std::cout << "... model_p " << model_p->GetSerNum() << " " << model_p << std::endl;
	    CChain *chain_p;
	    int nchains = model_p->GetNumberOfChains();
	    std::cout << "... n_chains: " << nchains << std::endl;
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       std::cout << "... chain_p: " << chain_p << std::endl;
	       int nres = chain_p->GetNumberOfResidues();
	       std::cout << "... ... nres " << nres << std::endl;
	       CResidue *residue_p;
	       CAtom *at;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  std::cout << "... ... residue_p " << residue_p << std::endl;
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
	 int n_rounds_max = 5;

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
   std::map<coot::residue_spec_t, bool> residue_done; // markers for contact_residues

   std::cout << "n_sel_1:: " << n_sel_1 << std::endl;
   std::cout << "n_sel_2:: " << n_sel_2 << std::endl;
   std::cout << "SeekContacts() found ncontacts: " << ncontacts << std::endl;
   
   if (ncontacts) {
      for (int icon=0; icon<ncontacts; icon++) {
	 CAtom *atom_ref = atom_sel_1[contact[icon].id1];
	 CAtom *atom_mov = atom_sel_2[contact[icon].id2];
	 coot::residue_spec_t r_1(atom_ref);
	 coot::residue_spec_t r_2(atom_mov);
	 contact_residues[r_1].push_back(r_2);
	 // include an "already matched" flag with the same key
	 residue_done[r_1] = 0;
      }
   }

   std::map<coot::residue_spec_t, std::vector<coot::residue_spec_t> >::iterator iter;
   std::map<coot::residue_spec_t, bool>::iterator done_iter; // markers for contact_residues

   std::cout << "residue specs after contacts " << std::endl;
   for (iter=contact_residues.begin(); iter!=contact_residues.end(); iter++) {

      if (! residue_done[iter->first]) {

	 std::cout << "     " << iter->first << std::endl;

	 std::vector<coot::residue_spec_t> contiguous_frag_spec;

	 

	 // now set as already matched
	 residue_done[iter->first] = 1;
      }
   } 
   
   delete [] contact;
   return r;
} 

// move the moving model (with model number 2) in mol.
void
coot::lsq_improve::apply_matches(const std::vector<coot::lsq_range_match_info_t> &matches) {

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
