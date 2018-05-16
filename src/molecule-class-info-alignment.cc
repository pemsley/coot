
#include "molecule-class-info.h"

// 20180302 add PIR alignment parsing
void
molecule_class_info_t::associate_pir_alignment(const std::string &chain_id, const std::string &alignment) {

   if (alignment.size()) {
      pir_alignments[chain_id] = alignment;
   }
}

// apply the alignment
void
molecule_class_info_t::apply_pir_alignment(const std::string &chain_id) {

   std::map<std::string, coot::pir_alignment_t>::const_iterator it;
   std::cout << "INFO:: apply_pir_alignment " << chain_id << std::endl;

   it = pir_alignments.find(chain_id);
   if (it != pir_alignments.end()) {
      const coot::pir_alignment_t &a = it->second;

      if (a.matches.size() > 0) {

	 int bs = backups_state();
	 turn_off_backup();

	 const std::map<int, coot::pir_alignment_t::matched_residue_t> matches = a.matches[0];
	 std::map<int, coot::pir_alignment_t::matched_residue_t>::const_iterator it;
	 for(it=matches.begin(); it!=matches.end(); it++) {
	    int res_no = it->first;
	    const coot::pir_alignment_t::matched_residue_t &mr = it->second;
	    mmdb::Residue *residue_p = get_residue(chain_id, res_no, "");

	    if (residue_p) {
	       std::string to = "to";
	       if (mr.aligned == mr.target) to = "..";
	       std::cout << "INFO:: mutate " << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName()
			 << " from " << mr.aligned << " " << to << " " << mr.target << std::endl;
	       if (mr.aligned != mr.target) {
		  std::string new_residue_type = coot::util::single_letter_to_3_letter_code(mr.target);
		  mutate(res_no, "", chain_id, new_residue_type);
	       }
	    } else {
	       std::cout << "No residue for residue number " << res_no << std::endl;
	    }
	 }

	 if (bs)
	    turn_on_backup();
      }
   }
}
