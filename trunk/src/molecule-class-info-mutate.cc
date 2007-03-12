
/* src/molecule-class-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef _MSC_VER
#include <windows.h>
#undef AddAtom
#define AddAtomA AddAtom
#endif

#include <stdlib.h>

#include "clipper/core/xmap.h"
#include "CIsoSurface.h"

#include "clipper/core/hkl_compute.h"
#include "clipper/clipper-phs.h"
#include "clipper/core/map_utils.h" // Map_stats


#include "mmdb_manager.h"
// #include "mmdb-extras.h"
// #include "mmdb.h"
// #include "mmdb-crystal.h"

#include "graphics-info.h"

// #include "xmap-utils.h"
// #include "coot-coord-utils.hh"
// #include "coot-utils.hh"

#include "molecule-class-info.h"
#include "mmdb_align.h"
#include "mmdb_tables.h"

void
molecule_class_info_t::mutate_chain(const std::string &chain_id,
				    const coot::chain_mutation_info_container_t &mutation_info,
				    PCResidue *SelResidues,
				    int nSelResidues) {

   if (mutation_info.insertions.size() > 0 ||
       mutation_info.deletions.size() > 0 ||
       mutation_info.mutations.size() > 0) {
      make_backup();

      // Don't backup each mutation, insertion etc - just do it before
      // and after.
      short int save_backup_state = backup_this_molecule;
      backup_this_molecule = 0;

      std::cout << "mutate chain " << mutation_info.insertions.size()
		<< " insertions" << std::endl;
      std::cout << "mutate chain " << mutation_info.deletions.size()
		<< " deletions" << std::endl;
      std::cout << "mutate chain " << mutation_info.mutations.size()
		<< " mutations" << std::endl;

      // do the operations in this order:
      // mutations
      // deletions
      // insertions
      int n_mutations  = 0;
      int n_deletions  = 0;
      int n_insertions = 0;

      // But first let's save the original ordering
      std::vector<int> original_seqnums(nSelResidues);
      for (int i=0; i<nSelResidues; i++)
	 original_seqnums[i] = SelResidues[i]->GetSeqNum();

      // --------------- mutations ----------------
      for (unsigned int i=0; i<mutation_info.mutations.size(); i++) {
	 std::string target_type = mutation_info.mutations[i].second;
	 coot::residue_spec_t rs = mutation_info.mutations[i].first;
	 // std::cout << "DEBUG:: chains: " << chain_id << " " << rs.chain << std::endl;
	 int SelectionHandle = atom_sel.mol->NewSelection();
	 PCResidue *local_residues;
	 int local_n_selected_residues;
	 atom_sel.mol->Select(SelectionHandle, STYPE_RESIDUE, 0,
			      (char *) chain_id.c_str(),
			      rs.resno, (char *) rs.insertion_code.c_str(),
			      rs.resno, (char *) rs.insertion_code.c_str(),
			      "*", "*", "*", "*",
			      SKEY_NEW
			      );
	 atom_sel.mol->GetSelIndex(SelectionHandle, local_residues,
				   local_n_selected_residues);
	 std::vector<CResidue *> residues_vec(local_n_selected_residues);
	 for (int i=0; i<local_n_selected_residues; i++)
	    residues_vec[i] = local_residues[i];
	 atom_sel.mol->DeleteSelection(SelectionHandle);
	 if (local_n_selected_residues > 0) {
	    mutate(residues_vec[0], target_type);
	    n_mutations++;
	 } else {
	    std::cout << "ERROR:: bad select in mutations " << chain_id
		      << " " << rs.resno << " " << rs.insertion_code << std::endl;
	 }

	 // Nope.... Can't DeleteSelection after mods.
	 // atom_sel.mol->DeleteSelection(SelectionHandle);
      }

      // --------------- deletions ----------------
      std::vector<std::pair<CResidue *, int> > residues_for_deletion;
      for (unsigned int i=0; i<mutation_info.deletions.size(); i++) {
	 coot::residue_spec_t rs = mutation_info.deletions[i];
	 int SelectionHandle = atom_sel.mol->NewSelection();
	 PCResidue *local_residues;
	 int local_n_selected_residues;
	 atom_sel.mol->Select(SelectionHandle, STYPE_RESIDUE, 0,
			      (char *) chain_id.c_str(),
			      rs.resno, (char *) rs.insertion_code.c_str(),
			      rs.resno, (char *) rs.insertion_code.c_str(),
			      "*", "*", "*", "*",
			      SKEY_NEW
			      );
	 atom_sel.mol->GetSelIndex(SelectionHandle, local_residues,
				   local_n_selected_residues);
// 	 std::cout << "DEBUG:: delete section: residue spec :" << chain_id << ": "
// 		   << rs.resno << " :" << rs.insertion_code
// 		   << ": selected " << local_n_selected_residues
// 		   << " residues\n";
	 if (local_n_selected_residues > 0) {
	    std::cout << "DEBUG:: marking for deleting "
		      << local_residues[0]->GetChainID() << " "
		      << local_residues[0]->GetSeqNum()  << " "
		      << local_residues[0]->GetResName() << "\n";
	    n_deletions++;
	    residues_for_deletion.push_back(std::pair<CResidue *, int> (local_residues[0], rs.resno));
	 }
	 atom_sel.mol->DeleteSelection(SelectionHandle);

      }

      for (unsigned int ird=0; ird<residues_for_deletion.size(); ird++) {
	 delete residues_for_deletion[ird].first;
	 residues_for_deletion[ird].first = NULL;
	 // now renumber:
	 for (int ires=0; ires<nSelResidues; ires++) {
	    if (SelResidues[ires]) {
	       if (original_seqnums[ires] > residues_for_deletion[ird].second) {
		  SelResidues[ires]->seqNum--;
	       }
	    } else {
	       std::cout << "INFO:: found a null residue at " << ires
			 << std::endl;
	    }
	 }
      }
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      // Is this dangerous? Yes, here (after mods) it certainly is.

      // --------------- insertions ----------------
      for (unsigned int i=0; i<mutation_info.insertions.size(); i++) {
	 coot::mutate_insertion_range_info_t r = mutation_info.insertions[i];
	 int offset = r.types.size();
	 n_insertions++;
	 for (int ires=0; ires<nSelResidues; ires++) {
	    if (SelResidues[ires]) {
	       if (original_seqnums[ires] > r.start_resno) {
		  SelResidues[ires]->seqNum += offset;
	       }
	    }
	 }
      }

      std::cout << "Applied " << n_mutations << " mutations " << std::endl;
      std::cout << "Applied " << n_deletions << " deletions " << std::endl;
      std::cout << "Applied " << n_insertions << " insertions " << std::endl;
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
      backup_this_molecule = save_backup_state;
   }

}


void
molecule_class_info_t::align_and_mutate(const std::string chain_id,
					const coot::fasta &fasta_seq) {

   std::string target = fasta_seq.sequence;

   CMMDBManager *mol = atom_sel.mol;
   if (mol) { 
      int selHnd = mol->NewSelection();
      PCResidue *SelResidues = NULL;
      int nSelResidues;
   
      mol->Select(selHnd, STYPE_RESIDUE, 0,
		  (char *) chain_id.c_str(),
		  ANY_RES, "*",
		  ANY_RES, "*",
		  "*",  // residue name
		  "*",  // Residue must contain this atom name?
		  "*",  // Residue must contain this Element?
		  "*",  // altLocs
		  SKEY_NEW // selection key
		  );
      mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      if (nSelResidues > 0) {

	 // I don't know if we can do this here, but I do know we
	 // mol->DeleteSelection(selHnd); // can't DeleteSelection after mods.
	 mutate_chain(chain_id,
		      align_on_chain(chain_id, SelResidues,
				     nSelResidues, target),
		      SelResidues, nSelResidues);
      }
   } else {
      std::cout << "ERROR:: null mol in align_and_mutate" << std::endl;
   }
}

coot::chain_mutation_info_container_t
molecule_class_info_t::align_on_chain(const std::string &chain_id,
				      PCResidue *SelResidues,
				      int nSelResidues,
				      const std::string &target) const {

   coot::chain_mutation_info_container_t ch_info(chain_id);

   std::cout << "\n\n-------------------------------------------------" << std::endl;
   std::cout << "        chain " << chain_id << std::endl;
   std::cout << "-------------------------------------------------\n\n" << std::endl;

   std::vector<std::pair<CResidue *, int> > vseq =
      coot::util::sort_residues_by_seqno(SelResidues, nSelResidues);

   // old style
   // std::string model = make_model_string_for_alignment(SelResidues, nSelResidues);
   // new style
   std::string model = coot::util::model_sequence(vseq);
   // std::cout << "INFO:: " << model << std::endl;

   CAlignment align;
   align.SetScores(0.5, -0.2);; // 2.0, -1 are the defaults
   align.Align((char *)model.c_str(), (char *)target.c_str());

   ch_info.alignedS = align.GetAlignedS();
   ch_info.alignedT = align.GetAlignedT();

   std::cout << ">  ";
   std::cout << name_ << std::endl;
   std::cout << align.GetAlignedS() << std::endl;
   std::cout << "> target seq: \n" << align.GetAlignedT() << std::endl;
   std::cout << "INFO:: alignment score " << align.GetScore() << std::endl;

   // get something like:
   // DVSGTVCLSALPPEATDTLNLIASDGPFPYSQD----F--------Q-------NRESVLPTQSYGYYHEYTVITP
   // DTSGTVCLS-LPPEA------IASDGPFPYSQDGTERFDSCVNCAWQRTGVVFQNRESVLPTQSYGYYHEYTVITP

   // Consider
   // pdb seq: "AVSGTVCLSA"
   // target: "TAVSGTVCLSAT"
   // s ->    "-AVSGTVCLSA-"
   // indexing 001234567899

   // So GetAlignedS() and GetAlignedT() return strings of the same length
   //
   // We run across the length of the returned string, looking for differences:

   // model  has residue             :  Simple mutate.
   // target has different residue   :  

   // model  has residue             :  There is an insertion in the model.
   // target has "-"                 :  Delete that residue, offset residues numbers 
   //                                :  to the right by -1.
   // 
   // model  has "-"                 :  The is a deletion in the model.
   // target has residue             :  Insert a gap for the residue by offsetting 
   //                                :  residue numbers to the right by +1.
   // 
   // s is from source (our pdb file)
   // t is the target sequence
   
   std::string s=align.GetAlignedS();
   std::string t=align.GetAlignedT();

   // we need to match the poistion in SelResidues to the position
   // after alignment (it has had "-"s inserted into it).
   std::vector<int> selindex(s.length());
   int sel_offset = 0;
   for (unsigned int iseq_indx=0; iseq_indx<s.length(); iseq_indx++) {
      selindex[iseq_indx] = iseq_indx - sel_offset;
//       std::cout << "assigned: selindex[" << iseq_indx << "]=" << selindex[iseq_indx]
// 		<< std::endl;
      if (s[iseq_indx] == '-')
	 sel_offset++;
   }

   if (s.length() == t.length()) {

      for (unsigned int iseq_indx=0; iseq_indx<s.length(); iseq_indx++) {
      }
      
      std::vector<int> resno_offsets(s.length(), 0);
//       std::cout << "DEBUG:: s.length() " << s.length() << std::endl;
//       std::cout << "DEBUG:: nSelResidues " << nSelResidues << std::endl;
      std::string inscode("");
      int ires = 0;
      
      for (unsigned int iseq_indx=0; iseq_indx<s.length(); iseq_indx++) {
	 if (s[iseq_indx] != t[iseq_indx]) {

// 	    std::cout << "DEBUG:: iseq_indx: " << iseq_indx 
// 		      << " selindex[" << iseq_indx << "]="
// 		      << selindex[iseq_indx] << std::endl;
	    
	    // These only make sense when the aligned residue (in s) was not "-"
	    if (s[iseq_indx] != '-') { 
	       ires            = SelResidues[selindex[iseq_indx]]->GetSeqNum();
	       inscode = SelResidues[selindex[iseq_indx]]->GetInsCode();
	    }

	    //	    std::cout << "DEBUG:: ires: " << ires << std::endl;
	    
	    // Case 1: (simple mutate)
	    if ((s[iseq_indx] != '-') && t[iseq_indx] != '-') {
// 	       std::cout << "mutate res number " << ires << " "
// 			 << s[iseq_indx] << " to " << t[iseq_indx] << std::endl;
	       std::string target_type =
		  coot::util::single_letter_to_3_letter_code(t[iseq_indx]);
	       coot::residue_spec_t res_spec(ires);
	       ch_info.add_mutation(res_spec, target_type);
	    }

	    // Case 2: model had insertion
	    if ((s[iseq_indx] != '-') && t[iseq_indx] == '-') {
	       for (unsigned int i=iseq_indx+1; i<s.length(); i++)
		  resno_offsets[i] -= 1;
// 	       std::cout << "Delete residue number " << iseq_indx << " "
// 			 << s[iseq_indx] << std::endl;
	       coot::residue_spec_t res_spec(ires);
	       ch_info.add_deletion(res_spec);
	    }

	    // Case 3: model has a deletion
	    if ((s[iseq_indx] == '-') && t[iseq_indx] != '-') {
	       for (unsigned int i=iseq_indx+1; i<s.length(); i++)
		  resno_offsets[i] += 1;
	       // ires will be for the previous residue.  It was not
	       // set for this one.
	       coot::residue_spec_t res_spec(ires+1);
	       std::string target_type =
		  coot::util::single_letter_to_3_letter_code(t[iseq_indx]);
	       ch_info.add_insertion(res_spec, target_type);
// 	       std::cout << "Insert residue number " << iseq_indx << " "
// 			 << t[iseq_indx] << std::endl;
	    }
	 }
      }
   }
   ch_info.rationalize_insertions();
   return ch_info;
}


// redundant now that we have coot-util functions.
//
std::string
molecule_class_info_t::make_model_string_for_alignment(PCResidue *SelResidues,
						       int nSelResidues) const {

   std::vector<std::pair<CResidue *, int> > vseq =
      coot::util::sort_residues_by_seqno(SelResidues, nSelResidues);
   return coot::util::model_sequence(vseq);

}



// Here is something that does DNA/RNA
int
molecule_class_info_t::mutate_base(const coot::residue_spec_t &res_spec, std::string type) {

   int istat=0;
   std::string refmac_nuc_type = type;
   // we match the requested residue type to the residue type that is
   // in the standard residues file.
   if (refmac_nuc_type.length() == 1) {
      if (refmac_nuc_type == "A")
	 refmac_nuc_type = "Ad";
      if (refmac_nuc_type == "G")
	 refmac_nuc_type = "Gr";
      if (refmac_nuc_type == "C")
	 refmac_nuc_type = "Cd";
      if (refmac_nuc_type == "T")
	 refmac_nuc_type = "Td";
      if (refmac_nuc_type == "U")
	 refmac_nuc_type = "Ur";
   } 
   
   if (atom_sel.n_selected_atoms > 0) { 

      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 CChain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p != NULL) {  
		  std::string mol_chain_id = chain_p->GetChainID();
		  if (mol_chain_id == res_spec.chain) { 
		     int nres = chain_p->GetNumberOfResidues();
		     PCResidue residue_p;
		     for (int ires=0; ires<nres; ires++) { 
			residue_p = chain_p->GetResidue(ires);
			if (residue_p->GetSeqNum() == res_spec.resno) {
			   if (res_spec.insertion_code == residue_p->GetInsCode()) {

			      // Found the residue (nucleotide in this case):

			      CResidue *std_base =
				 get_standard_residue_instance(refmac_nuc_type);
			      if (std_base) { 
				 mutate_base_internal(residue_p, std_base);
				 istat = 1;
			      } else {
				 std::cout << "Oops - can't find standard residue for type "
					   << type << std::endl;
			      }
			   }
			}
			if (istat)
			   break;
		     }
		  }
	       }
	       if (istat)
		  break;
	    }
	 }
      }
   }
   if (istat) { 
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked();
   }
   return istat; 
}

			

// Here std_base is at some arbitary position when passed.
// 
void
molecule_class_info_t::mutate_base_internal(CResidue *residue, CResidue *std_base) {

   make_backup();
   std::cout << "mutate_base_internal: pre:  "
	     << residue->GetNumberOfAtoms() << std::endl;
   coot::util::mutate_base(residue, std_base);
   std::cout << "mutate_base_internal: post: "
	     << residue->GetNumberOfAtoms() << std::endl;
   have_unsaved_changes_flag = 1;
}

int
molecule_class_info_t::exchange_chain_ids_for_seg_ids() {

   int istat = 0;
   short int changed_flag = 0;
   if (atom_sel.n_selected_atoms > 0) {

      // Make a list/vector of Chains so that we can delete them after
      // we have created the new chains.
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 CChain *chain_p = 0;
	 std::vector<int> chain_vec;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_vec.push_back(ichain);
	 }
      
	 //
	 CAtom *at = 0;
	 std::vector<std::pair<std::vector<CAtom *>, std::string> > atom_chain_vec;
	 std::vector<CAtom *> running_atom_vec;
	 std::string current_seg_id = "unassigned-camel-burger";
	 std::string current_chain_id = "unassigned-camel-burger";
	 for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
	    at = atom_sel.atom_selection[iat];
	    std::string this_atom_segid = at->segID;
	    std::string this_atom_chain_id = at->GetChainID();
	    if (current_seg_id != this_atom_segid) {
	       // add the running chain to the list then.  We construct a chain-id
	       // for it.
	       if (running_atom_vec.size() > 0) {
		  std::string constructed_seg_id = current_chain_id;  // (that of the previous atom)
		  constructed_seg_id += ":";
		  constructed_seg_id += current_seg_id; // also previous atom
		  atom_chain_vec.push_back(std::pair<std::vector<CAtom *>, std::string>(running_atom_vec, constructed_seg_id));
	       }
	       // start a new chain
	       running_atom_vec.clear();
	       running_atom_vec.push_back(at);
	       // set for next push
	       current_seg_id = this_atom_segid;
	       current_chain_id = this_atom_chain_id;
	    } else {
	       // all in all it's just a-nother atom in the chain...
	       running_atom_vec.push_back(at);
	    } 
	 }
	 if (running_atom_vec.size() > 0) {
	    atom_chain_vec.push_back(std::pair<std::vector<CAtom *>, std::string>(running_atom_vec, current_seg_id));
	 }
	 
	 // OK, so we have vector of vectors of atoms.  We need to make
	 // new atoms and residues to put them in.
	 std::cout << "INFO:: Creating " << atom_chain_vec.size() << " new chains\n";
	 for (unsigned int inch=0; inch<atom_chain_vec.size(); inch++) {
	    CChain *chain_p = new CChain;
	    const char *chid = atom_chain_vec[inch].second.c_str();
	    chain_p->SetChainID(chid);
	    // Add the chain to the model:
	    model_p->AddChain(chain_p);
	    CResidue *residue_p = 0;
	    int prev_resno = -999999;
	    char *prev_ins = "";
	    for (unsigned int iat=0; iat<atom_chain_vec[inch].first.size(); iat++) {
	       at = atom_chain_vec[inch].first[iat];
	       char *resname = at->GetResName();
	       char *ins     = at->GetInsCode();
	       int seqnum      = at->GetSeqNum();
	       if ((seqnum != prev_resno) || strcmp(ins,prev_ins)) {
		  // we need to make a new residue attached to chain_p
// 		  std::cout << "debug:: Making new residue " << resname << " "
// 			    << seqnum << " :" << ins << ": prev_resno " 
// 			    << prev_resno << " prev_ins :" << prev_ins << ":"
// 			    << std::endl;
		  residue_p = new CResidue(chain_p, resname, seqnum, ins);
		  prev_resno = seqnum;
		  prev_ins = ins;
	       }
	       CAtom *atom_p = new CAtom(residue_p); // does an AddAtom
	       atom_p->Copy(at); // doesn't touch atom's res
	    }
	 }
	 // now (finally) delete the chains of the model:
	 for (int ich=0; ich<chain_vec.size(); ich++) {
	    model_p->DeleteChain(chain_vec[ich]);
	 }
      }
   }

   // Let's just imagine that it always happens...
   have_unsaved_changes_flag = 1;
   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   make_bonds_type_checked();
   return changed_flag;
}
