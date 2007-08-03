/* ligand/ideal-rna.cc
 * 
 * Copyright 2006 The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include <iostream>

#include "coot-coord-utils.hh"
#include "ideal-rna.hh"

coot::ideal_rna::ideal_rna(const std::string &RNA_or_DNA, const std::string &form,
			   short int single_stranged_flag_in,
			   const std::string &sequence, CMMDBManager *standard_residues_in) {

   RNA_or_DNA_ = RNA_or_DNA;
   form_ = form;
   seq = sequence;
   standard_residues = standard_residues_in;
   single_stranged_flag = single_stranged_flag_in;

} 

// return a null pointer on something bad happened.
CMMDBManager *
coot::ideal_rna::make_molecule() {

   CMMDBManager *mol = 0; 
   CResidue *ur;
   short int is_dna_flag;
   coot::ideal_rna::form_t form_flag = A_FORM;

   if (RNA_or_DNA_ == "RNA") { 
      // in a special place
      is_dna_flag = 0;
   } else {
      is_dna_flag = 1;
   }

   if (form_ == "B") { 
      ur = get_standard_residue_instance("Cd", standard_residues);
      form_flag = B_FORM;
   } else {
      form_flag = A_FORM;
      ur = get_standard_residue_instance("Gr", standard_residues);
      if (is_dna_flag) {
	 delete_o2_prime(ur);
      }
   }

   if (ur) {
      mol = new CMMDBManager;
      CModel *model_p = new CModel;
      CChain *sense_chain_p = new CChain;
      sense_chain_p->SetChainID("A");
      model_p->AddChain(sense_chain_p);

      clipper::Mat33<double> antisense_base_mat(1, 0, 0, 0, -1, 0, 0, 0, -1);
      clipper::RTop_orth antisense_base_rtop(antisense_base_mat,
					     clipper::Coord_orth(0,0,0));
      CResidue *antisense_ref = coot::util::deep_copy_this_residue(ur, "", 1);
      // now transform antisense base to the right place:
      coot::util::transform_atoms(antisense_ref, antisense_base_rtop);

      for(int iseq=0; iseq<seq.length(); iseq++) {
	 if (is_valid_base(seq[iseq])) { 

	    // sense residue
	    CResidue *res = coot::util::deep_copy_this_residue(ur, "", 1);
	    res->seqNum = 1 + iseq ;
	    sense_chain_p->AddResidue(res);
	    clipper::RTop_orth o = n_turns(iseq, seq.length(), form_flag);
	    coot::util::transform_atoms(res, o);
	    mutate_res(res, seq[iseq], is_dna_flag);
	 }
      }

      if (! single_stranged_flag) { 
	 CChain *antisense_chain_p = new CChain;
	 // People complain that the resdiues of the antisense chain
	 // appear in the pdb file in the wrong order (high residue
	 // numbers first). So lets just make a vector of them, push
	 // them back and add them backwards at the end.
	 //
	 std::vector<CResidue *> residues_v;
	 antisense_chain_p->SetChainID("B");
	 model_p->AddChain(antisense_chain_p);
	 for(int iseq=0; iseq<seq.length(); iseq++) {
	    if (is_valid_base(seq[iseq])) { 

	       // antisense residue
	       CResidue *res = coot::util::deep_copy_this_residue(antisense_ref, "", 1);
	       res->seqNum = seq.length() - iseq;
	       // antisense_chain_p->AddResidue(res);  "backwards in pdb"
	       residues_v.push_back(res);  // modern
	       clipper::RTop_orth o = n_turns(iseq, seq.length(), form_flag);
	       coot::util::transform_atoms(res, o);
	       char complimentary_base = antisense_base(seq[iseq], is_dna_flag);
	       if (complimentary_base) 
		  mutate_res(res, complimentary_base, is_dna_flag);
	    }
	 }
 	 if (residues_v.size() > 0) {
	    // test must fail, unsigned is bad then
 	    for (int i=int(residues_v.size()-1); i>=0; i--) {
 	       antisense_chain_p->AddResidue(residues_v[i]);
 	    }
 	 } 
      }
      mol->AddModel(model_p);
      mol->FinishStructEdit();
      mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      
   } else {
      std::cout << "WARNING:: Bad standard residue Ur/Td" << std::endl;
   }
   return mol;
}


clipper::RTop_orth
coot::ideal_rna::n_turns(int nbase, int n_in_chain, coot::ideal_rna::form_t form_flag) const {

   clipper::Mat33<double> m_identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
   clipper::Coord_orth t_identity(0, 0, 0);
   clipper::RTop_orth o(m_identity, t_identity);

   // A form
   clipper::Mat33<double> base_mat(0.8415, -0.5402, 0,
				   0.5402, 0.8415, 0,
				   0, 0, 1);
   clipper::Coord_orth base_trans(0, 0, 2.81);
   clipper::RTop_orth one_base(base_mat, base_trans);

   if (form_flag == coot::ideal_rna::B_FORM) {
      // change one_base
      clipper::Mat33<double> b_base_mat(0.809,   -0.5878, 0,
					0.5878,   0.809,  0,
					0, 0, 1);
      clipper::Coord_orth b_base_trans(0, 0, 3.38);
      one_base = clipper::RTop_orth(b_base_mat, b_base_trans);
   } 

   for(int ibase=0; ibase<nbase; ibase++) {
      clipper::RTop_orth o1(one_base * o);
      o = o1;
   }

   return o;

}

short int
coot::ideal_rna::is_valid_base(char base) const {

   if (base == 'g') {
      return 1;
   } else {
      if (base == 'c') {
	 return 1;
      } else {
	 if (base == 'a') {
	    return 1;
	 } else {
	    if (base == 't') {
	       return 1;
	    } else {
	       if (base == 'u') {
		  return 1;
	       }
	    }
	 }
      }
   }
   return 0;
}

char
coot::ideal_rna::antisense_base(char base, short int is_dna_flag) const {

   if (base == 'g') {
      return 'c';
   } else {
      if (base == 'c') {
	 return 'g';
      } else {
	 if (base == 'a') {
	    if (is_dna_flag) 
	       return 't';
	    else 
	       return 'u';
	 } else {
	    if (base == 't') {
	       return 'a';
	    } else {
	       if (base == 'u') {
		  return 'a';
	       }
	    }
	 }
      }
   }
   return 0;
}


//
// Get a deep copy:
// return NULL on failure
// 
CResidue *
coot::ideal_rna::get_standard_residue_instance(const std::string &residue_type,
					       CMMDBManager *standard_residues) const {

   CResidue *std_residue = 0;
   int selHnd = standard_residues->NewSelection();
   standard_residues->Select (selHnd, STYPE_RESIDUE, 1, // .. TYPE, iModel
			       "*", // Chain(s) it's "A" in this case.
			       ANY_RES,"*",  // starting res
			       ANY_RES,"*",  // ending res
			       (char *) residue_type.c_str(),  // residue name
			       "*",  // Residue must contain this atom name?
			       "*",  // Residue must contain this Element?
			       "*",  // altLocs
			       SKEY_NEW // selection key
			       );
   PPCResidue SelResidue;
   int nSelResidues;
   standard_residues->GetSelIndex(selHnd, SelResidue, nSelResidues);
   
   if (nSelResidues != 1) {
      std::cout << "This should never happen - ";
      std::cout << "badness in get_standard_residue_instance, we selected "
		<< nSelResidues
		<< " residues looking for residues of type :"
		<< residue_type << ":\n";
   } else {
      std_residue = coot::util::deep_copy_this_residue(SelResidue[0], "", 1);
   }
   standard_residues->DeleteSelection(selHnd);
   return std_residue;
}

// return a status:
int 
coot::ideal_rna::mutate_res(CResidue *res, char base, short int is_dna_flag) const {

   // we need to get instances of bases from the standard residues
   int status = 0;

   std::string residue_type = "None";
   if (is_dna_flag) { 
      if (base == 'a')
	 residue_type = "Ad";
      if (base == 'g')
	 residue_type = "Gd";
      if (base == 't')
	 residue_type = "Td";
      if (base == 'c')
	 residue_type = "Cd";
   } else {
      if (base == 'a')
	 residue_type = "Ar";
      if (base == 'g')
	 residue_type = "Gr";
      if (base == 'u')
	 residue_type = "Ur";
      if (base == 'c')
	 residue_type = "Cr";
   }

   CResidue *std_res = get_standard_residue_instance(residue_type, standard_residues);
   if (std_res) {
      coot::util::mutate_base(res, std_res);
      status = 1;
   }
   return status;
}

void
coot::ideal_rna::delete_o2_prime(CResidue *res) const {

   int natoms;
   PPCAtom residue_atoms;
   short int deleted = 0;
   
   res->GetAtomTable(residue_atoms, natoms);
   for (int i=0; i<natoms; i++) {
      std::string atname(residue_atoms[i]->name);
      if (atname == " O2*") {
	 res->DeleteAtom(i);
	 deleted=1;
      } 
      if (atname == " O2'") {
	 res->DeleteAtom(i);
	 deleted=1;
      }
   }
   if (deleted)
      res->TrimAtomTable();
}
