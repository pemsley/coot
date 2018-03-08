
#include "molecule-class-info.h"

// 20180302 add PIR alignment parsing
void
molecule_class_info_t::associate_pir_alignment(const std::string &chain_id, const std::string &alignment) {

   if (alignment.size()) {
      pir_alignments[chain_id] = alignment;
   }
}

std::pair<bool,int>
molecule_class_info_t::max_res_no_in_chain(mmdb::Chain *chain_p) const {

   bool status = false;
   int res_no_max = -9999;
   if (chain_p) {
      int n = chain_p->GetNumberOfResidues();
      for (int i=0; i<n; i++) { // do this backwards?
	 mmdb::Residue *residue_p = chain_p->GetResidue(i);
	 int resno_this = residue_p->GetSeqNum();
	 if (resno_this > res_no_max) {
	    res_no_max = resno_this;
	    status = true;
	 }
      }
   }
   return std::pair<bool, int> (status, res_no_max);
}


// apply the alignment
void
molecule_class_info_t::apply_pir_alignment(const std::string &chain_id) {

   std::map<std::string, coot::pir_alignment_t>::const_iterator it;

   it = pir_alignments.find(chain_id);
   if (it != pir_alignments.end()) {
      const coot::pir_alignment_t &a = it->second;

      if (a.matches.size() > 0) {

	 mmdb::Chain *chain_p = 0;
	 int imod = 1;
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
 	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p_this = model_p->GetChain(ichain);
	       std::string this_chain_id = chain_p_this->GetChainID();
	       if (this_chain_id == chain_id) {
		  chain_p = chain_p_this;
		  break;
	       }
	    }
	 }

	 if (chain_p) {

	    // so now we have the alignment and the chain to which it should be applied

	    int bs = backups_state();
	    turn_off_backup();

	    mmdb::PResidue *residues = 0;
	    int n_residues;
	    chain_p->GetResidueTable(residues, n_residues);
	    int i_res = 0; // to start with
	    int i_res_offset_counter = 0;

	    if (a.size() > 0) {
	       if (a.size(0) > 0) {
		  std::vector<coot::pir_alignment_t::matched_residue_t> matches = a.get_matches(0);

		  // here set i_res_offset_counter to be the serial number
		  // of the residue of the first residue number in the model
		  
		  for (std::size_t i_align_pair=0; i_align_pair<matches.size(); i_align_pair++) {

		     const coot::pir_alignment_t::matched_residue_t &mr = matches[i_align_pair];
		     // std::cout << "Here 2 " << i_align_pair << " " << mr << std::endl;

		     if (mr.aligned == '-') {

			// That means there is no residue in the model for this pair.
			// We will need to build this loop. We can't mutate it. Move along
			// the sequence

		     } else {

			std::string pir_res_type = coot::util::single_letter_to_3_letter_code(mr.aligned);

			int i_res_index = i_res + i_res_offset_counter;
			if (i_res_index < n_residues) {
			   mmdb::Residue *residue_p = residues[i_res_index];

			   bool found_matching_residue = false;
			   while (!found_matching_residue) {
			      if (residue_p) {
				 std::string mol_res_type = residue_p->GetResName();
				 if (false)
				    std::cout << "looking for \"" << mr.aligned << "\" \""
					      << pir_res_type << "\" found "
					      << mol_res_type << std::endl;
				 if (pir_res_type == mol_res_type) {
				    found_matching_residue = true;

				    // now actually mutate (if needed)
				    std::string to = "to";
				    if (mr.aligned == mr.target) to = "..";
				    std::cout << "INFO:: mutate " << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName()
					      << " from " << mr.aligned << " " << to << " " << mr.target << std::endl;
				    if (mr.aligned != mr.target) {
				       std::string new_residue_type = coot::util::single_letter_to_3_letter_code(mr.target);
				       mutate(residue_p, new_residue_type);
				    }
				 }
			      }

			      // try something new
			      i_res++;
			      i_res_index = i_res + i_res_offset_counter;
			      if (i_res_index<n_residues) {
				 residue_p = residues[i_res_index];
			      } else {
				 break; // the while loop
			      }
			   }

			} else {
			   break; // we have gone past the last residue of the chain
			}
		     }
		  }
	       }
	    }
	    
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	    if (bs)
	       turn_on_backup();
	 }
      }
   }
}
