
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
   std::cout << "apply_pir_alignment " << chain_id << std::endl;

   it = pir_alignments.find(chain_id);
   if (it != pir_alignments.end()) {
      const coot::pir_alignment_t &a = it->second;

      if (a.matches.size() > 0) {

	 std::cout << "Here 1" << std::endl;

	 int bs = backups_state();
	 turn_off_backup();

	 const std::map<int, coot::pir_alignment_t::matched_residue_t> matches = a.matches[0];
	 std::map<int, coot::pir_alignment_t::matched_residue_t>::const_iterator it;
	 int res_no_offset = 0; // initially we have found no gaps in our model

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

	    std::pair<bool, int> res_no_max = max_res_no_in_chain(chain_p);
	    if (res_no_max.first) {
	       for(it=matches.begin(); it!=matches.end(); it++) {
		  int res_no = it->first;
		  const coot::pir_alignment_t::matched_residue_t &mr = it->second;

		  bool match = false;

		  while (match == false) {

		     int res_no_for_checking = res_no + res_no_offset;
		     if (res_no_for_checking > res_no_max.second)
			break;

		     mmdb::Residue *residue_p = get_residue(chain_id, res_no_for_checking, "");

		     if (residue_p) {

			std::string pir_res_type = coot::util::single_letter_to_3_letter_code(mr.aligned);
			std::string mol_res_type = residue_p->GetResName();

			if (pir_res_type == mol_res_type) {
			   match = true;

			   std::string to = "to";
			   if (mr.aligned == mr.target) to = "..";
			   std::cout << "INFO:: mutate " << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName()
				     << " from " << mr.aligned << " " << to << " " << mr.target << std::endl;
			   if (mr.aligned != mr.target) {
			      std::string new_residue_type = coot::util::single_letter_to_3_letter_code(mr.target);
			      mutate(res_no, "", chain_id, new_residue_type);
			   }
			} else {
			   res_no_offset++;
			}
		     } else {
			std::cout << "No residue for residue number " << res_no << std::endl;
			// res_no_offset++;
			break;
		     }
		  }
	       }
	    }
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	 }

	 if (bs)
	    turn_on_backup();
      }
   }
}
