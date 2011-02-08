
#include "coot-coord-utils.hh"

namespace coot {

   class lsq_improve{
      CMMDBManager *mol;
      int sel_hnd_1;
      int sel_hnd_2;
      int n_ref_atoms;
      int n_mov_atoms;
      int CAs_to_model(CMMDBManager *mol_in, int model_number);
      std::vector<coot::lsq_range_match_info_t> get_new_matches() const;
      // move the moving model (with model number 2) in mol.
      void apply_matches(const std::vector<coot::lsq_range_match_info_t> &matches);
      // have we reached convergence?
      bool same_matches(const std::vector<coot::lsq_range_match_info_t> &ranges_1,
			const std::vector<coot::lsq_range_match_info_t> &ranges_2) const;
      realtype crit_close;
      int n_res_for_frag;
      clipper::RTop_orth rtop_of_moving(const std::vector<coot::lsq_range_match_info_t> &matches) const;
   public:
      lsq_improve(CMMDBManager *mol_ref, CMMDBManager *mol_moving);
      void set_crit_close(realtype val) { crit_close = val; }
      void set_n_res_for_frag(int n_res_in) { n_res_for_frag = n_res_in; }
      void improve();
      
      // so, what is the RTop after we've done all that improving?
      // question to self: should this throw an exception?
      // 
      clipper::RTop_orth rtop_of_moving() const;
      
      ~lsq_improve() {
	 delete mol;
	 mol->DeleteSelection(sel_hnd_1);
	 mol->DeleteSelection(sel_hnd_2);
      }
   };
}


