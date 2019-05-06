/* ligand/residue_by_phi_psi.cc
 *
 * Copyright 2005 by Paul Emsley, The University of York
 * Copyright 2015 by Medical Research Council
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

#include <fstream>

#include "utils/coot-utils.hh"
// no dependency on coords files
// #include "coords/mmdb-extras.h"
// #include "coords/mmdb.h"

#include "coot-utils/coot-coord-utils.hh" // for co()
#include "residue_by_phi_psi.hh"
#include "mini-mol/mini-mol-utils.hh"

coot::residue_by_phi_psi::residue_by_phi_psi(const std::string &terminus,
					     mmdb::Residue *res_p,
					     const std::string &chain_id_in, 
					     const std::string &res_type,
					     float b_factor_in) {
   
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   thread_pool_p = 0;
   n_threads = 0;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   set_acceptable_fit_fraction(0.0); // allow all solutions to be scored.

   chain_id      = chain_id_in;
   residue_type  = res_type;
   terminus_type = terminus;
   b_factor      = b_factor_in;
   residue_p     = res_p;
   upstream_neighbour_residue_p = 0;
   downstream_neighbour_residue_p = 0;
   init_phi_psi_plot(); // uses residue_type
   set_dont_test_rotations();  // for ligand fitting
   set_dont_write_solutions(); // for ligand fitting
   debug_trials_flag = false; // (don't) write as PDBs, the various phi/phi/phi/psi solutions
}

// This is the new externally called function
// 
coot::minimol::molecule
coot::residue_by_phi_psi::best_fit_phi_psi(int n_trials,
					   const clipper::Xmap<float> &xmap_in) {

   return best_fit_phi_psi(n_trials, false, false, xmap_in);

}

// This is the externally called function
// 
coot::minimol::molecule
coot::residue_by_phi_psi::best_fit_phi_psi(int n_trials,
					   bool do_rigid_body_refinement,
					   bool add_other_residue_flag,
					   const clipper::Xmap<float> &xmap_in) {

   coot::minimol::molecule m;
   int offset = 0;
   if (terminus_type == "C") 
      offset = 1;
   else 
      if (terminus_type == "N") 
	 offset = -1;
      else
	 if (terminus_type == "MN")
	    offset = -1; // act like a N terminal addition (we have
			 // both Ca to define the added residue's
			 // carbonyl O position).
	 else 
	    if (terminus_type == "MC") 
	       offset = 1;
            else
	       if (terminus_type == "singleton")
		  offset = 1;

   // std::cout << "--- in best_fit_phi_psi() with offset " << offset << std::endl;

   if (offset == 0) {
      std::cout <<  "not a terminal residue\n";
   } else {
      if (false) { // debug
	 std::cout << "INFO:: Fitting terminal residue ";
	 if (do_rigid_body_refinement)
	    std::cout << " with individual rigid body fitting.\n";
	 else
	    std::cout << " without individual rigid body fitting.\n";

	 std::cout << "--- in best_fit_phi_psi() calling fit_terminal_residue_generic() with n_trials "
		   << n_trials << std::endl;
      }
      minimol::fragment frag = fit_terminal_residue_generic(n_trials, offset,
							    do_rigid_body_refinement, xmap_in);

      if (false) {
	 std::string file_name = "frag-best-fit.pdb";
	 minimol::molecule mmm(frag);
	 mmm.write_file(file_name, 10);
      }

      if (add_other_residue_flag) {
	 m.fragments.push_back(frag);
      } else {
	 // Get rid of the other residue by finding the first residues
	 // with atoms, then jumping out.

	 int ifrag = m.fragment_for_chain(chain_id); // chain_id is class variable

	 if (terminus_type == "C" || terminus_type == "MC" || terminus_type == "singleton") {
	    for (int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
	       if (frag[ires].atoms.size() > 0) {
		  try { 
		     m.fragments[ifrag].addresidue(frag[ires],0);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "ERROR:: best_fit_phi_psi() " << rte.what() << std::endl;
		  } 
		  break;
	       }
	    }
	 } else { // terminus_type == "N" or "MN"
	    for (int ires=frag.max_residue_number(); ires>=frag.min_res_no(); ires--) {
	       if (frag[ires].atoms.size() > 0) {
		  try {
		     m.fragments[ifrag].addresidue(frag[ires],0);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "ERROR:: best_fit_phi_psi() " << rte.what() << std::endl;
		  }
		  break;
	       }
	    }
	 } 
      }
   }

   return m;
}


// ---------------------------- Modern (2018) interface ----------------
//                              used by multi_peptide
// ---------------------------------------------------------------------
//
coot::minimol::fragment
coot::residue_by_phi_psi::best_fit_phi_psi(int n_trials, int offset) {


   minimol::fragment f = fit_terminal_residue_generic(n_trials, offset, false, Xmap());
   return f;
}

coot::minimol::fragment
coot::residue_by_phi_psi::fit_terminal_residue_generic(int n_trials, int offset, 
						       bool do_rigid_body_refinement,
						       const clipper::Xmap<float> &xmap_in) {

   coot::minimol::fragment best_fragment; // the returned thing
   bool debug = false; // shall we write out hypothesis di-peptides?

   coot::minimol::residue rres(residue_p->GetSeqNum() + offset);

   if (false) {
      std::cout << "                     called fit_terminal_residue_generic() residue-seqnum: "
		<< residue_p->GetSeqNum() << " offset: " << offset << std::endl;
   }

   // the positions of the current residue (and the residue connecting to that)
   connecting_atoms_t current_res_pos = get_connecting_residue_atoms();

   if (current_res_pos.empty()) {
      std::cout << "WARNING:: Failed to find atoms of terminal residue." << std::endl;
      std::cout << "WARNING:: Something strange in coordinates!? " << std::endl;
   } else {

      float best_score = 0.0;
      bool two_residues_flag = true;
      if (terminus_type == "MC")
	 two_residues_flag = 0;
      if (terminus_type == "singleton")
	 two_residues_flag = 0;
      if (terminus_type == "MN")
	 two_residues_flag = 0;
      int next_residue_seq_num;
      if (terminus_type == "C" || terminus_type == "MC" || terminus_type == "singleton") {
	 next_residue_seq_num = residue_p->GetSeqNum() + 1;
      } else {
	 next_residue_seq_num = residue_p->GetSeqNum() - 1;
      }

      std::vector<std::pair<ligand_score_card, minimol::fragment> > results(n_trials);

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

      if (n_threads >= 1) {

	 auto tp_1 = std::chrono::high_resolution_clock::now();
	 int n_per_thread = n_trials/n_threads;
	 for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
	    int trial_idx_start = i_thread * n_per_thread;
	    int trial_idx_end   = trial_idx_start + n_per_thread;
	    if (i_thread == (n_threads - 1))
	       trial_idx_end = n_trials; // for loop uses iat_start and tests for < iat_end

	    // do_rigid_body_refinement is ignored
	    //
	    if (false)
	       std::cout << "DEBUG:: fit_terminal_residue_generic() dispatching trial set "
			 << trial_idx_start << " to " << trial_idx_end << "\n";
	    thread_pool_p->push(fit_terminal_residue_generic_trial_inner_multithread,
				trial_idx_start, trial_idx_end, offset, residue_p, next_residue_seq_num,
				terminus_type, residue_type, b_factor,
				std::ref(current_res_pos), std::ref(xmap_in), map_rms,
				std::ref(rama), rama_max,
				std::ref(rama_pro), rama_max_pro,
				debug_trials_flag, &results);
	 }

	 auto tp_2 = std::chrono::high_resolution_clock::now();
	 // wait for thread pool to finish jobs.
	 bool wait_continue = true;
	 while (wait_continue) {
	    std::this_thread::sleep_for(std::chrono::milliseconds(20));
	    if (thread_pool_p->n_idle() == thread_pool_p->size())
	       wait_continue = false;
	 }
	 auto tp_3 = std::chrono::high_resolution_clock::now();

	 // std::cout << "fit_terminal_residue_generic(): multithread results: " << std::endl;
	 for (int itrial=0; itrial<n_trials; itrial++) {
	    // std::cout << "   comparing scores: " << results[itrial].first.get_score() << " "
	    // << best_score << std::endl;
	    if (results[itrial].first.get_score() > best_score) {
	       best_score = results[itrial].first.get_score();
	       best_fragment = results[itrial].second;
	       if (false)
		  std::cout << "debug ... best fragment has first residue and n residues: "
			    << best_fragment.first_residue()
			    << " " << best_fragment.n_filled_residues() << std::endl;
	    }
	 }

	 // std::cout << "DEBUG:: multi-thread: best_score " << best_score << std::endl;
	 auto tp_4 = std::chrono::high_resolution_clock::now();
	 auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
	 auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
	 auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
	 // std::cout << "Timings: " << d21 << " " << d32 << " " << d43 << " ms" << std::endl;

      } else {

	 // standard non-threaded
	 auto tp_1 = std::chrono::high_resolution_clock::now();

	 fit_terminal_residue_generic_trial_inner_multithread(0, // dummy thread index
							      0, n_trials,
							      offset, residue_p,
							      next_residue_seq_num,
							      terminus_type, residue_type,
							      b_factor, current_res_pos, xmap_in,
							      map_rms,
							      rama, rama_max,
							      rama_pro, rama_max_pro,
							      debug_trials_flag,
							      &results);

	 // find the best
	 for (int itrial=0; itrial<n_trials; itrial++) {
	    const ligand_score_card &s    = results[itrial].first;
	    const minimol::fragment &frag = results[itrial].second;
	    if (s.get_score() > best_score) {
	       best_score = s.get_score();
	       best_fragment = frag;
	    }
	 }
	 auto tp_2 = std::chrono::high_resolution_clock::now();
	 auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
	 std::cout << "INFO:: Single Thread Trials Timings: " << d21 << " ms " << std::endl;
      }

#else

      {
	 fit_terminal_residue_generic_trial_inner_multithread(0, // dummy thread index
							      0, n_trials,
							      offset, residue_p,
							      next_residue_seq_num,
							      terminus_type, residue_type,
							      b_factor, current_res_pos, xmap_in,
							      map_rms,
							      rama, rama_max,
							      rama_pro, rama_max_pro,
							      debug_trials_flag,
							      &results);

	 // find the best
	 for (int itrial=0; itrial<n_trials; itrial++) {
	    const ligand_score_card &s    = results[itrial].first;
	    const minimol::fragment &frag = results[itrial].second;
	    if (s.get_score() > best_score) {
	       best_score = s.get_score();
	       best_fragment = frag;
	    }
	 }
      }

#endif

   }
   if (best_fragment.residues.size() == 0) { 
      std::cout << "WARNING! fit_terminal_residue_generic:"
		<< " best_fragment is empty" << std::endl;
   }
   return best_fragment;
}

// The function will do no rigid body refinement and only to the 2 residue addition version
//
// For the second time it seems that passing a pointer to a xmap doesn't do a copy (when used
// with *xmap) but having a const ref &xmap *does* do a copy.  i.e. pushing threads is
// realllly slow with const ref xmaps!
// Why?
//
// [Later] Oh... it's because threads copy their const ref arguments by default. To force
// a reference, use std::ref.  Done.  Super-fast thread pushing now.
//
// static
void
coot::residue_by_phi_psi::fit_terminal_residue_generic_trial_inner_multithread(int ithread_idx,
									       int itrial_start,
									       int itrial_end,
									       int offset,
									       mmdb::Residue *res_p,
									       int next_residue_seq_num,
									       const std::string &terminus_type,
									       const std::string &residue_type,
									       float b_factor_in,
									       const connecting_atoms_t &current_res_pos,
									       const clipper::Xmap<float> &xmap_in,
									       float map_rms_in,
									       const clipper::Ramachandran &rama, float rama_max,
									       const clipper::Ramachandran &rama_pro, float rama_max_pro,
									       bool debug_solutions,
									       std::vector<std::pair<ligand_score_card, minimol::fragment> > *results) {

   bool do_low_density_fingerprint_positions = true;

   // I copied the contents of generic_trial here - the idea is that
   // score_orientation() doesn't need a map that has been imported into a residue_by_phi_psi.
   //

   std::string chain_id = res_p->GetChainID();
   residue_by_phi_psi rphipsi(terminus_type, res_p, chain_id, residue_type, b_factor_in);
   bool done_unexpected_missing_phi_message = false;
   bool done_unexpected_missing_psi_message = false;

   coot::ligand_score_card current_best; // zero score

   for (int itrial=itrial_start; itrial<itrial_end; itrial++) {

      coot::minimol::fragment frag;

      const clipper::Coord_orth &next_n  = current_res_pos.N_pos;
      const clipper::Coord_orth &next_c  = current_res_pos.C_pos;
      const clipper::Coord_orth &next_ca = current_res_pos.CA_pos;

      std::pair<bool,double> phi_current = current_res_pos.get_phi();
      std::pair<bool,double> psi_current = current_res_pos.get_psi();

      double phi_conditional = 0;
      double psi_conditional = 0;

      if (offset == 1) {
	 // forwards, we have a psi, need to generate a psi to place the N
	 if (phi_current.first) {
	    psi_conditional = get_psi_by_random_given_phi(phi_current.second, rama); // in radians
	 } else {
	    if (! done_unexpected_missing_phi_message) {
	       // std::cout << "unexpected missing phi_current" << residue_spec_t(res_p) << std::endl;
	       done_unexpected_missing_phi_message = true;
	    }
	 }
      } else {
	 // backwards, we have a phi, need a psi
	 if (psi_current.first) {
	    phi_conditional = get_phi_by_random_given_psi(psi_current.second, rama); // in radians
	 } else {
	    if (! done_unexpected_missing_psi_message) {
	       // std::cout << "unexpected missing psi_current" << residue_spec_t(res_p) << std::endl;
	       done_unexpected_missing_psi_message = true;
	    }
	 }
      }

      // how can this work? get_phi_psi_by_random() is private.
      phi_psi_t pp1 = rphipsi.get_phi_psi_by_random(rama, rama_max, false);
      phi_psi_t pp2 = rphipsi.get_phi_psi_by_random(rama, rama_max, false);

      frag = rphipsi.make_2_res_joining_frag_new(chain_id, current_res_pos,
						 clipper::Util::rad2d(phi_conditional),
						 clipper::Util::rad2d(psi_conditional),
						 pp1, pp2,
						 res_p->GetSeqNum(), offset);  // + or - 1

      std::vector<minimol::atom *> atoms_p = frag.select_atoms_serial();
      coot::ligand_score_card s = rphipsi.score_orientation(atoms_p, xmap_in, true);

      if (do_low_density_fingerprint_positions) {
	 rphipsi.add_characteristic_low_points(&s, itrial, current_res_pos,
					       pp1, pp2, res_p, offset, next_n, next_ca, next_c, frag, xmap_in);
      }

      if (s.get_score() > current_best.get_score()) {
	 current_best = s;
	 std::pair<ligand_score_card, minimol::fragment> result(s, frag);

	 results->at(itrial) = result;

	 if (debug_solutions)
	    debug_trials(frag, itrial, offset, res_p->GetSeqNum(),
			 current_res_pos,
			 phi_conditional,
			 psi_conditional,
			 pp1, pp2, s);
      }
   }
}

// static
void
coot::residue_by_phi_psi::debug_trials(const coot::minimol::fragment &frag, int itrial,
				       int offset,
				       int current_seqnum, // of the residue to which we are adding other residues
				       const connecting_atoms_t &current_res_pos,
				       const double &phi_conditional,
				       const double &psi_conditional,
				       const phi_psi_t &pp1,
				       const phi_psi_t &pp2,
				       const coot::ligand_score_card &score_card) {

   // DEBUGGING:  Let's write a pdb file for this fragment
   // Then look at them all.  Are they sensibly placed?
   coot::minimol::molecule m_tmp;
   m_tmp.fragments.push_back(frag);

   if (true) { // I want statistics, not pdb trial files for the moment
      std::string tmp_filename = "phi-psi-";
      tmp_filename += util::int_to_string(itrial);
      tmp_filename += ".pdb";
      m_tmp.write_file(tmp_filename, 10);
   }

   if (offset == 1) {
      std::string values_fn = "rama-values/values.tab";
      std::ofstream f(values_fn.c_str(), std::ios::app);
      if (f) {
	 // we have res_current (N), res_new_first (N+1), res_new_second, (N+2)
	 //
	 clipper::Coord_orth res_new_first_n   = frag[current_seqnum+1][" N  "].pos;
	 clipper::Coord_orth res_new_first_ca  = frag[current_seqnum+1][" CA "].pos;
	 clipper::Coord_orth res_new_first_c   = frag[current_seqnum+1][" C  "].pos;
	 clipper::Coord_orth res_new_second_n  = frag[current_seqnum+2][" N  "].pos;
	 clipper::Coord_orth res_new_second_ca = frag[current_seqnum+2][" CA "].pos;
	 clipper::Coord_orth res_new_second_c  = frag[current_seqnum+2][" C  "].pos;
	 double psi_current_rad = clipper::Coord_orth::torsion(current_res_pos.N_pos,
							       current_res_pos.CA_pos,
							       current_res_pos.C_pos,
							       res_new_first_n);
	 double phi_r1_rad = clipper::Coord_orth::torsion(current_res_pos.C_pos,
							  res_new_first_n,
							  res_new_first_ca,
							  res_new_first_c);
	 double psi_r1_rad = clipper::Coord_orth::torsion(res_new_first_n,
							  res_new_first_ca,
							  res_new_first_c,
							  res_new_second_n);

	 double psi_current = clipper::Util::rad2d(psi_current_rad);
	 double phi_r_1 = clipper::Util::rad2d(phi_r1_rad);
	 double psi_r_1 = clipper::Util::rad2d(psi_r1_rad);

	 f << "itrial " << itrial
	   << " seqnum " << current_seqnum
	   << " score " << score_card.get_score()
	   << " input-psi " << psi_conditional << " output-psi " << psi_current
	   << " new phi input " << pp1.phi << " output " << phi_r_1
	   << " new psi input " << pp1.psi << " output " << psi_r_1
	   << "\n";

	 f.close();
      }
   }
}


// change the score card
// which N do we want?  It depends on which direction we are building - it seems?
// next_n when building forwards is wrong as it stands.
//
// This works fine when building on the N-terminus.
//
void
coot::residue_by_phi_psi::add_characteristic_low_points(coot::ligand_score_card *s,
							int itrial,
							const connecting_atoms_t &current_res_pos,
							const coot::phi_psi_t &pp1,
							const coot::phi_psi_t &pp2,
							mmdb::Residue *res_p,
							int offset,
							const clipper::Coord_orth &next_n,
							const clipper::Coord_orth &next_ca,
							const clipper::Coord_orth &next_c,
							const coot::minimol::fragment &frag,
							const clipper::Xmap<float> &xmap_in) const {

   bool debug = false;
   if (debug) {
      std::string file_name = "frag-" + util::int_to_string(itrial) + ".pdb";
      minimol::molecule m(frag);
      m.write_file(file_name, 10);
   }
   // the ~HA point (more extended than an H bond distance of course).
   // Recall that psi is the rotation of the Ns around the C-CA bond.
   // and the N position that we need for ld_pos_1 is the N of the
   // residue with resno resno_this+1.
   //
   // res_p is the residue that we are adding to
   // subject_res_num is the residue number of the residue that we
   // are about to add.
   //
   // If we are adding to the N-terminus, the N is in res_p.
   // If we are adding to the C-terminus, the N is in the second residue of the fragment
   //
   double a = clipper::Util::d2rad(pp1.tau-5); // maybe more like 120?
   double t = clipper::Util::d2rad(pp1.psi-125);
   // this is for N-terminal addition
   int neighb_seqnum = res_p->GetSeqNum();
   int subject_res_num = neighb_seqnum + offset;
   if (frag[subject_res_num].atoms.size() > 0) {
      clipper::Coord_orth c_pos  = frag[subject_res_num][" C  "].pos;
      clipper::Coord_orth ca_pos = frag[subject_res_num][" CA "].pos;
      clipper::Coord_orth n_pos  = frag[subject_res_num][" N  "].pos;

      // ca_pos is the bonded neighbor of ld_pos_1 (not next_n).
      clipper::Coord_orth ld_pos_1(next_n, c_pos, ca_pos, 1.8, a, t);

      // need to change if we are building on the C-terminus,
      // the N comes from the other residue, not the N of the
      // residue to which we are adding subject_res
      //
      if (offset == 1) {
	 // set the ~HA low density position as you would set a C,
	 // but rotate phi by ~ 109 degrees
	 t = clipper::Util::d2rad(pp1.phi+129); // checked
	 a = clipper::Util::d2rad(120);
	 const clipper::Coord_orth &c_pos = current_res_pos.C_pos;
	 ld_pos_1 = clipper::Coord_orth(c_pos, n_pos, ca_pos, 1.6, a, t);

	 if (false) {
	    std::cout << "debug:: reset ld_pos_1 to " << ld_pos_1.format() << std::endl;
	    std::cout << "        using CA pos " << ca_pos.format() << std::endl;
	    std::cout << "        using C  pos " << c_pos.format() << std::endl;
	 }
      }

      // ld_pos_2 ca-c-n 120 1.8 (pseudo torsion) in the direction of NH
      // The N of the NH is in res_p if we are building on the N terminus (obviously)
      //    next_n and next_c are passed in this case.
      // The N of the NH is in subject_res_num, as is the CA. The C is in res_p.
      //
      // setup as if building on N-terminus first
      double a_2 = clipper::Util::d2rad(120.0);
      double t_2 = clipper::Util::d2rad(180.0); // maybe we can have 2 of these with different torsions
                                                // say 155 and 205.
      clipper::Coord_orth ld_pos_2(c_pos, next_ca, next_n, 1.5, a_2, t_2);
      // or just 2 more...
      clipper::Coord_orth ld_pos_3(c_pos, next_ca, next_n, 1.8, a_2, t_2-clipper::Util::d2rad(50));
      clipper::Coord_orth ld_pos_4(c_pos, next_ca, next_n, 1.8, a_2, t_2+clipper::Util::d2rad(50));

      if (offset == 1) {
	 // building on C-terminus
	 clipper::Coord_orth n_pos = frag[subject_res_num][" N  "].pos;
	 ld_pos_2 = clipper::Coord_orth(next_c, ca_pos, n_pos, 1.5, a_2, t_2);
      }

      float dv_ld_pos_1 = score_position(ld_pos_1, xmap_in);
      float dv_ld_pos_2 = score_position(ld_pos_2, xmap_in);
      float dv_ld_pos_3 = score_position(ld_pos_3, xmap_in);
      float dv_ld_pos_4 = score_position(ld_pos_4, xmap_in);

      if (false) { // debugging
	 std::cout << "ld_pos_1: " << ld_pos_1.format() << " scores " << dv_ld_pos_1 << " c.f. " << s->score_per_atom
		   << " psi " << pp1.psi << " t: " << clipper::Util::rad2d(t) << std::endl;
	 std::cout << "ld_pos_2: " << ld_pos_2.format() << " score  " << dv_ld_pos_2 << std::endl;
	 std::ofstream f("debug-ld-pos.py");
	 f << "set_rotation_centre( " << next_n.x() << ", " << next_n.y() << ", " << next_n.z() << ")\n";
	 f << "set_rotation_centre( " << c_pos.x() << ", " << c_pos.y() << ", " << c_pos.z() << ")\n";
	 f << "set_rotation_centre( " << ca_pos.x() << ", " << ca_pos.y() << ", " << ca_pos.z() << ")\n";
	 f << "set_rotation_centre( " << ld_pos_1.x() << ", " << ld_pos_1.y() << ", " << ld_pos_1.z() << ")\n";
      }

      std::pair<clipper::Coord_orth, float> sp(ld_pos_1, dv_ld_pos_1);
      //       s->scored_characteristic_low_density_points.push_back(sp);
      //       s->atom_point_score -= 2.0 * dv_ld_pos_1; // scale needs optimizing
      //       s->atom_point_score -= 2.0 * dv_ld_pos_2; //  likewise
      //       s->atom_point_score -= 1.0 * dv_ld_pos_3; //  likewise
      //       s->atom_point_score -= 1.0 * dv_ld_pos_4; //  likewise
   } else {
      std::cout << "DEBUG:: oops in add_characteristic_low_points() "
		<< " residue subject_res_num "
		<< subject_res_num << " - No atoms" << std::endl;
   }
}



// If there are 3 positions on return their order is N, C, CA
// 
coot::residue_by_phi_psi::connecting_atoms_t
coot::residue_by_phi_psi::get_connecting_residue_atoms() const {

   connecting_atoms_t atoms_in_residue;

   // Note, it seems that if there are 2 atoms in a residue with the
   // selected name (which happens if there is an altconf), then
   // GetAtom() returns NULL.  Ugh.
   // 20180403-PE But is that something we care about in this class?

   // std::vector<clipper::Coord_orth> pos; old

   mmdb::Atom *N_at = 0;
   mmdb::Atom *C_at = 0;
   mmdb::Atom *CA_at = 0;

   mmdb::PPAtom residue_atoms = 0;
   int nResidueAtoms;
   residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      mmdb::Atom *at = residue_atoms[i];
      std::string atom_name(at->name);
      if (atom_name == " N  ") N_at = at;
      if (atom_name == " C  ") C_at = at;
      if (atom_name == " CA ") CA_at = at;
   }

   if (N_at) {
      if (C_at) {
	 if (CA_at) {
	    clipper::Coord_orth  N_at_pos = co(N_at);
	    clipper::Coord_orth CA_at_pos = co(CA_at);
	    clipper::Coord_orth  C_at_pos = co(C_at);
	    atoms_in_residue = connecting_atoms_t(N_at_pos, CA_at_pos, C_at_pos);
	 }
      }
   }

   if (upstream_neighbour_residue_p) {
      mmdb::PPAtom residue_atoms = 0;
      int nResidueAtoms;
      upstream_neighbour_residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
	 mmdb::Atom *at = residue_atoms[i];
	 std::string atom_name(at->GetAtomName());
	 if (atom_name == " C  ") {  // PDBv3 FIXME
	    clipper::Coord_orth pos = co(at);
	    atoms_in_residue.set_upstream_C(pos);
	    break;
	 }
      }
   } else {
      // std::cout << "............ no upstream_neighbour_residue_p" << std::endl;
   }

   if (downstream_neighbour_residue_p) {
      mmdb::PPAtom residue_atoms = 0;
      int nResidueAtoms;
      downstream_neighbour_residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
	 mmdb::Atom *at = residue_atoms[i];
	 std::string atom_name(at->GetAtomName());
	 if (atom_name == " N  ") {  // PDBv3 FIXME
	    clipper::Coord_orth pos = co(at);
	    atoms_in_residue.set_downstream_N(pos);
	    break;
	 }
      }
   }
   return atoms_in_residue;
}


// a forward-built residue is made from psi of the previous residue and a phi from "this" one.
// phi_this, psi_prev are in degrees (what a mess)
// 
coot::minimol::residue
coot::residue_by_phi_psi::construct_next_res_from_rama_angles(float phi_this, float psi_prev,
							      float tau, int seqno,
							      const connecting_atoms_t &current_res_pos) const {

   const clipper::Coord_orth &previous_n  = current_res_pos.N_pos;
   const clipper::Coord_orth &previous_ca = current_res_pos.CA_pos;
   const clipper::Coord_orth &previous_c  = current_res_pos.C_pos;

   // +/- 10 degrees (but in radians)
   double omega_jitter =  20.0 * (M_PI/180) * (util::random()/static_cast<double>(RAND_MAX) - 0.5);
   double O_torsion    = 720.0 * (M_PI/180) * (util::random()/static_cast<double>(RAND_MAX) - 0.5);

   double angle, torsion;

   coot::minimol::residue mres(seqno);
   mres.name = residue_type;

   // N  of next is placed by psi_previous
   // CA of next is placed by omega
   // C  of next is placed by phi_next
   //
   // So, for clarity, say we are building on residue 54 forwards
   // 55 N  is placed by psi of residue 54
   // 55 CA is places by omega
   // 55 C  is placed by phi of residue 55
   // 55 O  is placed by psi of residue 55 (because it places the N of 56)

   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(psi_prev);
   clipper::Coord_orth n_pos(previous_n, previous_ca, previous_c, 
			     1.329, angle, torsion); // C-N bond

   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(180.0) + omega_jitter;
   clipper::Coord_orth ca_pos(previous_ca, previous_c, n_pos,
			      1.458, angle, torsion); // N-CA bond

   angle = clipper::Util::d2rad(tau); // N-CA-C  // was 111.200
   torsion = clipper::Util::d2rad(phi_this);
   clipper::Coord_orth c_pos(previous_c, n_pos, ca_pos,
			     1.525, angle, torsion); // CA-C bond

   angle = clipper::Util::d2rad(120.800);  // CA-C-O
   torsion = O_torsion; // Random.  It depends on the *next* residue, not the prevous one.
                        // i.e. should be: ca_next, n_next, c_pos with torsion 0.0
   clipper::Coord_orth o_pos(n_pos, ca_pos, c_pos,
			     1.231, angle, torsion); // C-O bond

   mres.addatom(coot::minimol::atom(" N  ", " N", n_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" C  ", " C", c_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" CA ", " C", ca_pos, "", b_factor));
   mres.addatom(coot::minimol::atom(" O  ", " O", o_pos,  "", b_factor));

   return mres;
}

// phi is phi of the residue to which we are joining.
coot::minimol::residue
coot::residue_by_phi_psi::construct_prev_res_from_rama_angles(float phi, float psi, float tau,
							      int seqno, const connecting_atoms_t &current_res_pos) const {

   coot::minimol::residue mres(seqno);
   mres.name = residue_type;

   double angle, torsion;

   // C
   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(phi);
   clipper::Coord_orth c_pos(current_res_pos.C_pos,
			     current_res_pos.CA_pos,
			     current_res_pos.N_pos,
			     1.329, angle, torsion); // C-N bond

   // Ca
   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(180.0);
   clipper::Coord_orth ca_pos(current_res_pos.CA_pos, current_res_pos.N_pos, c_pos,
			      1.525, angle, torsion); // Ca-C bond

   // N
   angle = clipper::Util::d2rad(tau); // N-Ca-C
   torsion = clipper::Util::d2rad(psi);
   clipper::Coord_orth n_pos(current_res_pos.N_pos, c_pos, ca_pos,
			     1.458, angle, torsion); // Ca-N

   // O
   angle = clipper::Util::d2rad(120.800);  // Ca-C-O
   torsion = clipper::Util::d2rad(0.0);
//    clipper::Coord_orth o_pos(n_pos, ca_pos, c_pos,
// 			     1.231, angle, torsion); // C=0
   clipper::Coord_orth o_pos(current_res_pos.CA_pos, current_res_pos.N_pos, c_pos,
			     1.231, angle, torsion); // C=0


   // 20170925 there is no CBs. Why? The get added in execute_add_terminal_residue()
   // add_side_chain_to_terminal_res().
   // I've added it here now.
   //
   // 201880517-PE: Or have I?

   // CB
   angle = clipper::Util::d2rad(tau);
   torsion = clipper::Util::d2rad(psi+123.4);
   clipper::Coord_orth cb_pos(current_res_pos.N_pos, c_pos, ca_pos,
			     1.52, angle, torsion); // Ca-Cb


   // ---------------------------

   mres.addatom(coot::minimol::atom(" N  ", " N", n_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" C  ", " C", c_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" CA ", " C", ca_pos, "", b_factor));
   mres.addatom(coot::minimol::atom(" O  ", " O", o_pos,  "", b_factor));

   // I guess that we should be consistent with construct_next_res_from_rama_angles().
   // i.e. here is not the place to add the CB.
   //
   // mres.addatom(coot::minimol::atom(" CB ", " C", cb_pos, "", b_factor));

   return mres;
}


// retire this function. Building forwards and backwards treats the
// argments differently
//
// coot::minimol::residue
// coot::residue_by_phi_psi::construct_joining_res(const phi_psi_t &pp,
// 						int seqno,
// 						const connecting_atoms_t &connecting_atoms) const {

//    if (terminus_type == "C" || terminus_type == "MC" || terminus_type == "singleton")
//       return construct_next_res_from_rama_angles(pp.phi, pp.psi, pp.tau,
// 						 seqno, next_n, next_ca, next_c);
//    else // terminus_type == "N" or "MN"
//       return construct_prev_res_from_rama_angles(pp.phi, pp.psi, pp.tau,
// 						 seqno, next_n, next_ca, next_c);
// }

void
coot::residue_by_phi_psi::init_phi_psi_plot() {

   rama.init(clipper::Ramachandran::Gly5);

   rama_pro.init(clipper::Ramachandran::Pro5);

   rama_max = 0.0;
   rama_max_pro = 0.0;
   
   for (float phi=0.0; phi<360.0; phi+=3.0) {
      for (float psi=0.0; psi<360.0; psi+=3.0) {
	 float v = rama.probability(clipper::Util::d2rad(phi),
				    clipper::Util::d2rad(psi));
	 if (v > rama_max)
	    rama_max = v;
      }
   }
   for (float phi=0.0; phi<360.0; phi+=3.0) {
      for (float psi=0.0; psi<360.0; psi+=3.0) {
	 float v = rama_pro.probability(clipper::Util::d2rad(phi),
					clipper::Util::d2rad(psi));
	 if (v > rama_max_pro)
	    rama_max_pro = v;
      }
   }
}

// if is_pro_rama is not set correctly, this can fail to terminate
//
coot::phi_psi_t
coot::residue_by_phi_psi::get_phi_psi_by_random(const clipper::Ramachandran &rama_local,
						const float &rama_max_local,
						bool is_pro_rama) const {
   
   float phi, psi, r;

   // some fraction of the time we want (to force) phi values between 0 and 180
   // (most of the time we don't)
   bool force_positive_phi = false;
   if (! is_pro_rama) {
      float phi_test_random = float (coot::util::random())/float (RAND_MAX);
      if (phi_test_random < 0.2) // say
	 force_positive_phi = true;
   }

   for (;;) {
      phi = 360.0 * float (coot::util::random())/float (RAND_MAX);
      psi = 360.0 * float (coot::util::random())/float (RAND_MAX);

      r = rama_max_local * coot::util::random()/ float (RAND_MAX);
      float prob = rama_local.probability(clipper::Util::d2rad(phi),
					  clipper::Util::d2rad(psi));

      // std::cout << "            compare " << prob << " " << r << std::endl;

      if (prob > r) {
	 if (force_positive_phi) {
	    if (phi < 180.0) { // i.e. on the right hand side of a Rama plot (centred on 0,0)
	       break;
	    }
	 } else {
	    break;
	 }
      }
   }
   float u1 = static_cast<float>(coot::util::random())/static_cast<float>(RAND_MAX);
   float u2 = 2.0 * u1 - 1.0; // -1 -> 1.
   float tau = 111.2 + u2 * 6.2; // 4.0 needs investigation
   return phi_psi_t(phi, psi, tau);
}

// phi and psi in radians
// When we are building backwards, we have psi for the selected residue, but not phi.
// To be clear, we only have psi, when we have also have the downstream (higher residue number) residue
// as well.
// static
double
coot::residue_by_phi_psi::get_phi_by_random_given_psi(double psi, const clipper::Ramachandran &rama) {

   double phi;
   double one_over_rand_max = 1.0 / static_cast<double> (RAND_MAX);

   // find rand_max_local (for probability of phi given psi
   //
   std::vector<double> pr_phi(72);
   double two_pi = 2.0 * M_PI;
   double step = M_PI / 36.0;
   double conditional_pr_rama_max = 0.0;
   for (unsigned int i=0; i<72; i++) {
      double phi = (i+0.5)*step;
      double pr = rama.probability(phi, psi);
      if (pr > conditional_pr_rama_max)
	 conditional_pr_rama_max = pr;
   }
   
   for (;;) {
      phi = 2.0 * M_PI * util::random() * one_over_rand_max;
      double r = conditional_pr_rama_max * util::random() * one_over_rand_max;
      double prob = rama.probability(phi, psi);
      if (prob > r) {
	 break;
      }
   }
   return phi;
}

// phi and psi in radians
// When we are building forward, we have phi for the selected residue, but not psi.
// To be clear, we only have phi, when we have also have the upstream (lower residue number) residue
// as well.
// static
double
coot::residue_by_phi_psi::get_psi_by_random_given_phi(double phi, const clipper::Ramachandran &rama) {

   double psi;
   double one_over_rand_max = 1.0 / static_cast<double> (RAND_MAX);

   // find rand_max_local (for probability of phi given psi
   //
   std::vector<double> pr_psi(72);
   double two_pi = 2.0 * M_PI;
   double step = M_PI / 36.0;
   double conditional_pr_rama_max = 0.0;
   for (unsigned int i=0; i<72; i++) {
      double psi = (i+0.5)*step;
      double pr = rama.probability(phi, psi);
      // std::cout << "pr for i " << i << " " << pr << std::endl;
      if (pr > conditional_pr_rama_max)
	 conditional_pr_rama_max = pr;
   }

   if (conditional_pr_rama_max < 0.0001) {
      // something went wrong, hack a return value
      // so that we don't stay in the below loop forever
      //
      psi = 2.0 * M_PI * util::random() * one_over_rand_max;
   } else{

      for (;;) {
	 psi = 2.0 * M_PI * util::random() * one_over_rand_max;
	 double r = conditional_pr_rama_max * util::random() * one_over_rand_max;
	 double prob = rama.probability(phi, psi);
	 if (false) // debug
	    std::cout << "prob " << prob << " from phi: " << phi << " psi: " << psi << " cf r " << r
		      << " phi is " << phi << " conditional_pr_rama_max " << conditional_pr_rama_max
		      << std::endl;
	 if (prob > r) {
	    break;
	 }
      }
   }
   return psi;
}



void
coot::residue_by_phi_psi::debug_compare_check(const coot::minimol::residue &mres,
					      std::vector<minimol::atom *> atoms)  {

   int icount = 0;

   std::cout << "mres has " << mres.atoms.size() << " atoms, "
	     << "atoms has " << atoms.size() << " atoms." << std::endl;
   for (unsigned int iat=0; iat<mres.atoms.size(); iat++) {
      std::cout << "check " <<  mres[iat].pos.format()  << " vs. "
		<< atoms[icount]->pos.format() << std::endl;
      icount++;
   }

} 



// now the non-classed version of the next residue building functions:

coot::minimol::residue
coot::build_C_terminal_ALA(float phi, float psi, 
			   int seqno,
			   const clipper::Coord_orth &previous_n,
			   const clipper::Coord_orth &previous_ca,
			   const clipper::Coord_orth &previous_c,
			   float b_factor) {


   coot::minimol::residue mres(seqno);
   mres.name = "ALA"; // can be changed by caller

   double angle, torsion;
	 
   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(psi);
   clipper::Coord_orth n_pos(previous_n, previous_ca, previous_c, 
			     1.329, angle, torsion); // C-N bond

   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(180.0); 
   clipper::Coord_orth ca_pos(previous_ca, previous_c, n_pos,
			      1.458, angle, torsion); // N-CA bond

   angle = clipper::Util::d2rad(111.200); // N-CA-C
   torsion = clipper::Util::d2rad(phi);
   clipper::Coord_orth c_pos(previous_c, n_pos, ca_pos,
			     1.525, angle, torsion); // CA-C bond

   angle = clipper::Util::d2rad(120.800);  // CA-C-O
   torsion = clipper::Util::d2rad(90.0); // just a guess.  It depends
					  // on the *next* residue,
					  // not the prevous one.  i.e. should be:
                                          // ca_next, n_next, c_pos with torsion 0.0
   clipper::Coord_orth o_pos(n_pos, ca_pos, c_pos,
			     1.231, angle, torsion); // C-O bond

   mres.addatom(coot::minimol::atom(" N  ", " N", n_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" C  ", " C", c_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" CA ", " C", ca_pos, "", b_factor));
   mres.addatom(coot::minimol::atom(" O  ", " O", o_pos,  "", b_factor));

   return mres;
}   

coot::minimol::residue
coot::build_N_terminal_ALA(float phi, float psi, 
			   int seqno,
			   const clipper::Coord_orth &next_n,
			   const clipper::Coord_orth &next_ca,
			   const clipper::Coord_orth &next_c,
			   float b_factor) {
   
   coot::minimol::residue mres(seqno);
   mres.name = "ALA"; // can be change to UNK by caller

   double angle, torsion;

   // C
   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(phi);
   clipper::Coord_orth c_pos(next_c, next_ca, next_n,
			     1.329, angle, torsion); // C-N bond

   // Ca
   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(180.0);
   clipper::Coord_orth ca_pos(next_ca, next_n, c_pos,
			      1.525, angle, torsion); // Ca-C bond

   // N
   angle = clipper::Util::d2rad(111.200); // N-Ca-C
   torsion = clipper::Util::d2rad(psi);
   clipper::Coord_orth n_pos(next_n, c_pos, ca_pos,
			     1.458, angle, torsion); // Ca-N

   // O
   angle = clipper::Util::d2rad(120.800);  // Ca-C-O
   torsion = clipper::Util::d2rad(0.0);
//    clipper::Coord_orth o_pos(n_pos, ca_pos, c_pos,
// 			     1.231, angle, torsion); // C=0
   clipper::Coord_orth o_pos(next_ca, next_n, c_pos,
			     1.231, angle, torsion); // C=0

   // ---------------------------

   mres.addatom(coot::minimol::atom(" N  ", " N", n_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" C  ", " C", c_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" CA ", " C", ca_pos, "", b_factor));
   mres.addatom(coot::minimol::atom(" O  ", " O", o_pos,  "", b_factor));

   return mres;
}

// static
clipper::Coord_orth
coot::residue_by_phi_psi::best_fit_phi_psi_attaching_oxygen_position_update(const minimol::molecule &mm,
									    mmdb::Residue *residue_p) {

   clipper::Coord_orth pos(0,0,0);

   // These are my 3 atoms:
   // find the C atom in residue_p and the N and CA of the first residue in mm

   minimol::residue res_with_CA_C(residue_p);
   if (mm.fragments.size() == 1) {
      const minimol::fragment &frag = mm.fragments[0];
      int idx_first = frag.first_residue();
      const minimol::residue &res_with_N = frag[idx_first];
      for (std::size_t iat=0; iat<res_with_N.atoms.size(); iat++) {
	 std::cout << "    " << iat << " " << res_with_N[iat] << std::endl;
      }

      std::pair<bool, clipper::Coord_orth> o_pos = o_position(res_with_CA_C, res_with_N);
      if (o_pos.first)
	 pos = o_pos.second;
   }
   return pos;
}
