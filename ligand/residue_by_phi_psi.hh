/* ligand/residue_by_phi_psi.cc
 * 
 * Copyright 2005 by Paul Emsley, The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

#ifndef RESIDUE_BY_PHI_PSI_HH
#define RESIDUE_BY_PHI_PSI_HH

// 20180101-PE should I put this in ligand.hh?
#ifdef HAVE_BOOST
#ifdef HAVE_CXX_THREAD
#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#include "utils/ctpl.h"
#endif // HAVE_CXX_THREAD
#endif // HAVE_BOOST

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#include <thread>
#include <chrono>
#include "utils/ctpl.h" // match that included in simple-restraint.hh
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

#include "ligand.hh"
#include "clipper/core/ramachandran.h"

#include "ideal/phi-psi.hh"

namespace coot { 

   class residue_by_phi_psi : public ligand {

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

      int ires_terminus;
      std::string chain_id;
      std::string residue_type;
      std::string terminus_type;
      float rama_max, rama_max_pro;
      clipper::Ramachandran rama;
      clipper::Ramachandran rama_pro;
      float b_factor;


      // This is not const because GetAtom() fo mmdb::Residue is not a const
      // function.  Argh!
      //
      mmdb::Residue *residue_p; // the residue of the last atom (we
                                // clicked on an atom of it).
      mmdb::Residue *upstream_neighbour_residue_p;
      mmdb::Residue *downstream_neighbour_residue_p;

      bool debug_trials_flag;
      phi_psi_t get_phi_psi_by_random(const clipper::Ramachandran &rama_local,
				      const float &rama_max_local,
				      bool is_pro_rama) const;

      // When we are building backwards, we have psi for the selected residue, but not phi.
      // To be clear, we only have psi, when we have also have the downstream (higher residue number) residue
      // as well.
      static double get_phi_by_random_given_psi(double psi, const clipper::Ramachandran &rama); // phi and psi in radians
      // When we are building forward, we have phi for the selected residue, but not psi.
      // To be clear, we only have phi, when we have also have the upstream (lower residue number) residue
      // as well.
      static double get_psi_by_random_given_phi(double phi, const clipper::Ramachandran &rama); // phi and psi in radians
     void init_phi_psi_plot(); 

      // get rid of offset - we know it's the NEXT                                 
     minimol::residue                          
     construct_next_res_from_rama_angles(float phi, float psi, float tau,
					 int seqno, const connecting_atoms_t &current_res_pos) const;
     minimol::residue
     construct_prev_res_from_rama_angles(float phi, float psi, float tau,
					 int seqno, const connecting_atoms_t &current_res_pos) const;

//      minimol::residue construct_joining_res(const phi_psi_t &pp,
// 						  int seqno,
// 						  const clipper::Coord_orth &next_n,
// 						  const clipper::Coord_orth &next_ca,
// 						  const clipper::Coord_orth &next_c) const;

     connecting_atoms_t get_connecting_residue_atoms() const; 
     minimol::fragment fit_terminal_residue_generic(int n_trials, int offset,
						    bool do_rigid_body_refinement,
						    const clipper::Xmap<float> &xmap_in);

     // should this be static?
     std::pair<ligand_score_card, coot::minimol::fragment>
     fit_terminal_residue_generic_trial_inner(int itrial,
					      int offset,
					      int next_residue_seq_num,
					      const connecting_atoms_t &pos,
					      bool two_residues_flag,
					      bool do_rigid_body_refinement,
					      const clipper::Xmap<float> &xmap_in) const;

     static
     void fit_terminal_residue_generic_trial_inner_multithread(int ithread_idx,
							       int itrial_start,
							       int itrial_end,
							       int offset,
							       mmdb::Residue *res_p,
							       int next_residue_seq_num,
							       const std::string &terminus_type,
							       const std::string &residue_type,
							       float b_factor_in,
							       const connecting_atoms_t &pos,
							       const clipper::Xmap<float> &xmap_in,
							       float map_rms_in,
							       const clipper::Ramachandran &rama, float rama_max,
							       const clipper::Ramachandran &rama_pro, float rama_max_pro,
							       bool debug_solutions,
							       std::vector<std::pair<ligand_score_card, minimol::fragment> > *results,
                                                               std::atomic<unsigned int> &thread_count);

      // retire this at some stage                                                       
//      minimol::fragment
//      make_2_res_joining_frag(const std::string &chain_id,
// 			     const phi_psi_t &pp1,
// 			     const phi_psi_t &pp2,
// 			     int seqnum,
// 			     int offset, // + or - 1
// 			     const connecting_atoms_t &current_atom_positions) const;

     minimol::fragment
     make_2_res_joining_frag_new(const std::string &chain_id,
				 const connecting_atoms_t &current_atom_positions,
				 const double &phi_conditional,
				 const double &psi_conditional,
				 const phi_psi_t &pp1,
				 const phi_psi_t &pp2,
				 int seqnum,
				 int offset // + or - 1
				 ) const;
     minimol::fragment
     make_2_res_joining_frag_new_building_forwards(const std::string &chain_id,
						   const connecting_atoms_t &current_atom_positions,
						   const double &psi_conditional,
						   const phi_psi_t &pp1,
						   const phi_psi_t &pp2,
						   int seqnum) const;
     minimol::fragment
     make_2_res_joining_frag_new_building_backwards(const std::string &chain_id,
						    const connecting_atoms_t &current_atom_positions,
						    const double &phi_conditional,
						    const phi_psi_t &pp1,
						    const phi_psi_t &pp2,
						    int seqnum) const;
     void
     add_characteristic_low_points(coot::ligand_score_card *s,
				   int itrial,
				   const connecting_atoms_t &current_res_pos,
				   const coot::phi_psi_t &p1,
				   const coot::phi_psi_t &p2,
				   mmdb::Residue *residue_p,
				   int offset,
				   const clipper::Coord_orth &next_n,
				   const clipper::Coord_orth &next_ca,
				   const clipper::Coord_orth &next_c, // for forward building low density points
				   const coot::minimol::fragment &frag,
				   const clipper::Xmap<float> &xmap_in) const;

     void debug_compare_check(const coot::minimol::residue &mres,
			      std::vector<minimol::atom *> atoms_p);

     
   public:
      residue_by_phi_psi(const std::string &terminus_type, // "N", or "C"
			 mmdb::Residue *res_p,
			 const std::string &chain_id, 
			 const std::string &res_type,
			 float b_factor_in);
      void set_downstream_neighbour(mmdb::Residue *r) { downstream_neighbour_residue_p = r; }
      void set_upstream_neighbour(mmdb::Residue *r) { upstream_neighbour_residue_p = r; }

      void write_trial_pdbs() { debug_trials_flag = true; }

      static void debug_trials(const coot::minimol::fragment &frag, int itrial,
			       int offset,
			       int current_seqnum,
			       const connecting_atoms_t &current_res_pos,
			       const double &phi_conditional,
			       const double &psi_conditional,
			       const phi_psi_t &pp1,
			       const phi_psi_t &pp2,
			       const coot::ligand_score_card &score_card);

      minimol::molecule best_fit_phi_psi(int n_trials,
					 const clipper::Xmap<float> &xmap_in);

      // This is the externally called function (from graphics-info-modelling)
      //
      minimol::molecule best_fit_phi_psi(int n_trials,
					 bool do_rigid_body_refinement,
					 bool add_other_residue_flag,
					 const clipper::Xmap<float> &xmap_in);

      static clipper::Coord_orth best_fit_phi_psi_attaching_oxygen_position_update(const minimol::molecule &mm,
										   mmdb::Residue *residue_p);

      // offset: N or C addition (-1 or 1).
      minimol::fragment best_fit_phi_psi(int n_trials, int offset); 

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
     ctpl::thread_pool *thread_pool_p;
     unsigned int n_threads;
     void thread_pool(ctpl::thread_pool *tp_in, int n_threads_in) {
	thread_pool_p = tp_in;
	n_threads = n_threads_in;
     }
   
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

  };

  minimol::residue
  build_N_terminal_ALA(float phi, float psi, int seqno,
		       const clipper::Coord_orth &previous_n,
		       const clipper::Coord_orth &previous_ca,
		       const clipper::Coord_orth &previous_c,
		       float b_factor); 
   minimol::residue
   build_C_terminal_ALA(float phi, float psi, int seqno,
			const clipper::Coord_orth &next_n,
			const clipper::Coord_orth &next_ca,
			const clipper::Coord_orth &next_c,
			float b_factor);

   minimol::fragment
   multi_build_N_terminal_ALA(const mmdb::Manager *mol_in,
			      const std::string &terminus_type, // "N", or "C"
			      mmdb::Residue *res_p,
			      const std::string &chain_id, 
			      const std::string &res_type,
			      float b_factor_in,
			      int n_trials);
} // namespace coot

#endif // RESIDUE_BY_PHI_PSI_HH
