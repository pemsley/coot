#ifndef NEW_RESIDUE_BY_3_PHI_PSI_HH
#define NEW_RESIDUE_BY_3_PHI_PSI_HH

#include <string>
#include <clipper/core/xmap.h>
#include <clipper/core/ramachandran.h>
#include <mmdb2/mmdb_manager.h>

#include "utils/ctpl.h"
#include "mini-mol/mini-mol.hh"
#include "ideal/phi-psi.hh"

#include "ligand.hh"

namespace coot {

   // probably there is a better hhome for this function
   //
   float get_random_float(); // betwween 0 and 1
   
   class new_residue_by_3_phi_psi {

      // this is copied from residue_by_phi_psi. Perhaps extract it.
      class connecting_atoms_t {
	 bool filled_flag;
      public:
	 connecting_atoms_t() { filled_flag = false; }
	 connecting_atoms_t(const clipper::Coord_orth &N_pos_in,
			    const clipper::Coord_orth &CA_pos_in,
			    const clipper::Coord_orth &C_pos_in) :
	    N_pos(N_pos_in), CA_pos(CA_pos_in), C_pos(C_pos_in) { filled_flag = true; }
	 void set_upstream_C(const clipper::Coord_orth &C_pos_in) {
	    upstream_C.first = true;
	    upstream_C.second = C_pos_in;
	 }
	 void set_downstream_N(const clipper::Coord_orth &N_pos_in) {
	    downstream_N.first = true;
	    downstream_N.second = N_pos_in;
	 }
	 clipper::Coord_orth N_pos;
	 clipper::Coord_orth CA_pos;
	 clipper::Coord_orth C_pos;
	 // so that we can calculate phi of the existing residue (from residue N-1)
	 std::pair<bool, clipper::Coord_orth> upstream_C; // for building forwards
	 // so that we can calculate psi of the existing residue (from residue N+1)
	 std::pair<bool, clipper::Coord_orth> downstream_N; // for building backwards
	 bool empty() const { return ! filled_flag; }
	 std::pair<bool, double> get_phi() const {
	    if (upstream_C.first) {
	       double t = clipper::Coord_orth::torsion(upstream_C.second, N_pos, CA_pos, C_pos);
	       return std::pair<bool, double> (true, t);
	    } else {
	       return std::pair<bool, double> (false, -999.9);
	    }
	 }
	 std::pair<bool, double> get_psi() const {
	    if (downstream_N.first) {
	       double t = clipper::Coord_orth::torsion(N_pos, CA_pos, C_pos, downstream_N.second);
	       return std::pair<bool, double> (true, t);
	    } else {
	       return std::pair<bool, double> (false, -999.9);
	    }
	 }
      };

      ctpl::thread_pool *thread_pool_p;
      unsigned int n_threads;
      std::string chain_id;
      std::string terminus_type;
      mmdb::Residue *upstream_neighbour_residue_p;
      mmdb::Residue *downstream_neighbour_residue_p;
      mmdb::Residue *residue_p;
      float rama_max;
      float rama_pro_max;
      clipper::Ramachandran rama;
      clipper::Ramachandran rama_pro;

      void init_phi_psi_plot();
      connecting_atoms_t get_connecting_residue_atoms() const;

      // When we are building backwards, we have psi for the selected residue, but not phi.
      // To be clear, we only have psi, when we have also have the downstream (higher residue number) residue
      // as well.
      static double get_phi_by_random_given_psi(double psi, const clipper::Ramachandran &rama); // phi and psi in radians
      // When we are building forward, we have phi for the selected residue, but not psi.
      // To be clear, we only have phi, when we have also have the upstream (lower residue number) residue
      // as well.
      static double get_psi_by_random_given_phi(double phi, const clipper::Ramachandran &rama); // phi and psi in radians

      static phi_psi_t get_phi_psi_by_random(const clipper::Ramachandran &rama,
				      const float &rama_max,
				      bool is_pro_rama);

      static
      minimol::fragment make_3_res_joining_frag_forward(const std::string &chain_id, const connecting_atoms_t &current_res_pos,
                                                        const double &psi_conditional_deg,
                                                        const phi_psi_t &pp_1, const phi_psi_t &tpp_2, const phi_psi_t &pp_3,
                                                        int seq_num);

      static
      minimol::fragment make_3_res_joining_frag_backward(const std::string &chain_id, const connecting_atoms_t &current_res_pos,
                                                         const double &phi_conditional_deg,
                                                         const phi_psi_t &pp_1, const phi_psi_t &tpp_2, const phi_psi_t &pp_3,
                                                         int seq_num);

      static minimol::residue                          
      construct_next_res_from_rama_angles(float phi, float psi, float tau, int seqno, const connecting_atoms_t &current_res_pos, float occupancy);
      
      static minimol::residue                          
      construct_prev_res_from_rama_angles(float phi, float psi, float tau, int seqno, const connecting_atoms_t &current_res_pos, float occupancy);
      
      static float score_fragment_basic(const minimol::fragment &frag,
                                        const connecting_atoms_t &current_res_pos,
                                        const clipper::Xmap<float> &xmap);
      static float score_fragment_using_peptide_fingerprint(const minimol::fragment &frag,
                                                            const connecting_atoms_t &current_res_pos,
                                                            const clipper::Xmap<float> &xmap,
                                                            int res_no_base, int i_trial); // pass i_trial for debugging

   public:
      new_residue_by_3_phi_psi(const std::string &terminus_type, mmdb::Residue *residue_p, const std::string &chain_id);
      void add_thread_pool(ctpl::thread_pool  *thread_pool_p, unsigned int n_threads);
      minimol::fragment best_fit_phi_psi(unsigned int n_trials, const clipper::Xmap<float> &xmap) const;
      void set_downstream_neighbour(mmdb::Residue *r) { downstream_neighbour_residue_p = r; }
      void set_upstream_neighbour(mmdb::Residue *r) { upstream_neighbour_residue_p = r; }
   };
}

#endif // NEW_RESIDUE_BY_3_PHI_PSI_HH
