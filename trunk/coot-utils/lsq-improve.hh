
#include "coot-coord-utils.hh"

namespace coot {

   class lsq_improve{
      int n_rounds_max;
      CMMDBManager *mol;
      CMMDBManager *mol_initial_copy;
      int sel_hnd_1;
      int sel_hnd_2;
      int n_ref_atoms;
      int n_mov_atoms;
      int CAs_to_model(CMMDBManager *mol_in, int model_number);
      
      // the crit_close has a multipler that is dependent on
      // round_number so that at high round number the criterion for
      // close atoms is reduced (i.e. made more tight).
      // 
      std::vector<lsq_range_match_info_t> get_new_matches(int round_number, int rounds_max, bool summary_to_screen_flag=0) const;
      std::vector<lsq_range_match_info_t> get_new_matches(const std::map<residue_spec_t, std::vector<residue_spec_t> > &contact_residues) const;
      // move the moving model (with model number 2) in mol.
      void apply_matches(const std::vector<lsq_range_match_info_t> &matches);
      realtype crit_close;
      int n_res_for_frag;
      clipper::RTop_orth rtop_of_moving(const std::vector<lsq_range_match_info_t> &matches) const;
   public:
      lsq_improve(CMMDBManager *mol_ref, CMMDBManager *mol_moving);
      void set_crit_close(realtype val) { crit_close = val; }
      void set_n_res_for_frag(int n_res_in) { n_res_for_frag = n_res_in; }
      void improve();
      
      // so, what is the RTop after we've done all that improving?
      // 
      // this can throw an exception.
      // 
      clipper::RTop_orth rtop_of_moving() const;
      
      ~lsq_improve() {
	 delete mol;
	 mol->DeleteSelection(sel_hnd_1);
	 mol->DeleteSelection(sel_hnd_2);
      }
   };
}


