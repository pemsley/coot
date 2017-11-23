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

#include "utils/coot-utils.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#include "residue_by_phi_psi.hh"


coot::residue_by_phi_psi::residue_by_phi_psi(const std::string &terminus,
					     mmdb::Residue *res_p,
					     const std::string &chain_id_in, 
					     const std::string &res_type,
					     float b_factor_in) {
   
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   thread_pool_p = 0;
   n_threads = 0;
#endif
   chain_id      = chain_id_in;
   residue_type  = res_type;
   terminus_type = terminus;
   b_factor      = b_factor_in;
   residue_p     = res_p;
   init_phi_psi_plot();
   set_dont_test_rotations();
   set_dont_write_solutions();
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

coot::minimol::fragment
coot::residue_by_phi_psi::best_fit_phi_psi(int n_trials, int offset) {

   std::cout << "                      called (fragment) best_fit_phi_psi() with offset "
	     << offset << " -- is this even used? " << std::endl;
   minimol::fragment f = fit_terminal_residue_generic(n_trials, offset, false, Xmap());
   return f;
}

coot::minimol::fragment
coot::residue_by_phi_psi::fit_terminal_residue_generic(int n_trials, int offset, 
						       bool do_rigid_body_refinement,
						       const clipper::Xmap<float> &xmap_in) {

   coot::minimol::fragment best_fragment; // the returned thing
   bool debug = false; // shall we write out hypothesis di-peptides?
   
   // float angle;
   // float torsion;
   coot::minimol::residue rres(residue_p->GetSeqNum() + offset);
   
   std::cout << "                     called fit_terminal_residue_generic() residue-seqnum: "
	     << residue_p->GetSeqNum() << " offset: " << offset << std::endl;

   std::vector<clipper::Coord_orth> pos = get_connecting_residue_atoms();

   if (pos.size() != 3) { 
      std::cout << "WARNING:: Failed to find atoms of terminal residue, found "
		<< pos.size() << " atoms" << std::endl;
      std::cout << "WARNING:: Something strange in coordinates!? " << std::endl;
   } else {

      clipper::Coord_orth next_n  = pos[0];
      clipper::Coord_orth next_c  = pos[1];
      clipper::Coord_orth next_ca = pos[2];

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

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

      if (n_threads > 1) {
	 auto tp_1 = std::chrono::high_resolution_clock::now();
	 std::vector<std::pair<ligand_score_card, minimol::fragment> > results(n_trials);
	 int n_per_thread = n_trials/n_threads;
	 for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
	    int trial_idx_start = i_thread * n_per_thread;
	    int trial_idx_end   = trial_idx_start + n_per_thread;
	    if (i_thread == (n_threads - 1))
	       trial_idx_end = n_trials; // for loop uses iat_start and tests for < iat_end

	    // do_rigid_body_refinement is ignored
	    //
	    std::cout << "dispatching trial set " << trial_idx_start << " to " << trial_idx_end << "\n";
	    thread_pool_p->push(fit_terminal_residue_generic_trial_inner_multithread,
				trial_idx_start, trial_idx_end, offset, residue_p, next_residue_seq_num,
				terminus_type, residue_type, b_factor,
				pos, std::ref(xmap_in), map_rms, &results);
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
	 // std::cout << "multithread results: " << std::endl;
	 for (int itrial=0; itrial<n_trials; itrial++) {
	    // std::cout << "mt   " << itrial << " " << results[itrial].first.atom_point_score << std::endl;
	    if (results[itrial].first.atom_point_score > best_score) {
	       best_score = results[itrial].first.atom_point_score;
	       best_fragment = results[itrial].second;
	    }
	 }
	 auto tp_4 = std::chrono::high_resolution_clock::now();
	 auto d21 = chrono::duration_cast<chrono::milliseconds>(tp_2 - tp_1).count();
	 auto d32 = chrono::duration_cast<chrono::milliseconds>(tp_3 - tp_2).count();
	 auto d43 = chrono::duration_cast<chrono::milliseconds>(tp_4 - tp_3).count();
	 std::cout << "Timings: " << d21 << " " << d32 << " " << d43 << " ms" << std::endl;
      } else {
	 // standard non-threaded
	 auto tp_1 = std::chrono::high_resolution_clock::now();
	 for (int itrial=0; itrial<n_trials; itrial++) {
	    std::pair<ligand_score_card, minimol::fragment> fitted_thing =
	       fit_terminal_residue_generic_trial_inner(itrial, offset, next_residue_seq_num,
							pos, two_residues_flag, do_rigid_body_refinement,
							xmap_in);
	    const ligand_score_card &s = fitted_thing.first;
	    const minimol::fragment &frag = fitted_thing.second;
	    if (s.atom_point_score > best_score) {
	       best_score = s.atom_point_score;
	       best_fragment = frag;
	    }
	 }
	 auto tp_2 = std::chrono::high_resolution_clock::now();
	 auto d21 = chrono::duration_cast<chrono::milliseconds>(tp_2 - tp_1).count();
	 std::cout << "Timings: " << d21 << " ms " << std::endl;
      }

#else

      for (int itrial=0; itrial<n_trials; itrial++) {
	 std::pair<ligand_score_card, minimol::fragment> fitted_thing = 
	    fit_terminal_residue_generic_trial_inner(itrial, offset, next_residue_seq_num, pos, two_residues_flag, do_rigid_body_refinement, xmap_in);

	 const ligand_score_card &s    = fitted_thing.first;
	 const minimol::fragment &frag = fitted_thing.second;
	 if (s.atom_point_score > best_score) {
	    best_score = s.atom_point_score;
	    best_fragment = frag;
	 }
      } // ends the n_trials for-loop
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
									       const std::vector<clipper::Coord_orth> &pos,
									       const clipper::Xmap<float> &xmap_in,
									       float map_rms_in,
									       std::vector<std::pair<ligand_score_card, minimol::fragment> > *results) {

   // I copied the contents of generic_trial here - the idea is that
   // score_orientation() doesn't need a map that has been imported into a residue_by_phi_psi.
   //

   std::string chain_id = res_p->GetChainID();
   residue_by_phi_psi rphipsi(terminus_type, res_p, chain_id, residue_type, b_factor_in);

   for (int itrial=itrial_start; itrial<itrial_end; itrial++) {
      // std::pair<coot::ligand_score_card, coot::minimol::fragment> result;

      // replaced this with explicit function details
      // rphipsi.fit_terminal_residue_generic_trial_inner(itrial, offset, next_residue_seq_num, pos, true, true);
      // to use the reference xmap_in.

      coot::minimol::fragment frag;
      clipper::Coord_orth next_n  = pos[0];
      clipper::Coord_orth next_c  = pos[1];
      clipper::Coord_orth next_ca = pos[2];

      phi_psi_pair p1 = rphipsi.get_phi_psi_by_random();
      phi_psi_pair p2 = rphipsi.get_phi_psi_by_random();
      frag = rphipsi.make_2_res_joining_frag(chain_id, p1, p2,
					     res_p->GetSeqNum(),
					     offset,  // + or - 1
					     next_n,
					     next_ca,
					     next_c);

      std::vector<minimol::atom *> atoms_p = frag.select_atoms_serial();
      coot::ligand_score_card s = rphipsi.score_orientation(atoms_p, xmap_in);
      std::pair<coot::ligand_score_card, coot::minimol::fragment> result(s, frag);

      bool do_low_density_fingerprint_positions = false;
      if (do_low_density_fingerprint_positions) {
	 rphipsi.add_characteristic_low_points(&s, itrial, p1, p2, res_p, next_n, frag, xmap_in);
      }

      (*results)[itrial] = result;
   }

}


// change the score card
void
coot::residue_by_phi_psi::add_characteristic_low_points(coot::ligand_score_card *s,
							int itrial,
							const coot::phi_psi_pair &p1,
							const coot::phi_psi_pair &p2,
							mmdb::Residue *res_p,
							const clipper::Coord_orth &next_n,
							const coot::minimol::fragment &frag,
							const clipper::Xmap<float> &xmap_in) const {

   bool debug = false;
   if (debug) {
      std::string file_name = "frag-" + util::int_to_string(itrial) + ".pdb";
      minimol::molecule m(frag);
      m.write_file(file_name, 10);
   }
   // the ~HA point (more extended than an H bond distance of course)
   double a = clipper::Util::d2rad(p1.tau-5); // maybe more like 120?
   double t = clipper::Util::d2rad(p1.psi-125);
   // this is for N-terminal addition
   int neighb_seqnum = res_p->GetSeqNum();
   int subject_res_num = neighb_seqnum-1;
   if (frag[subject_res_num].atoms.size() > 0) {
      clipper::Coord_orth c_pos  = frag[subject_res_num][" C  "].pos;
      clipper::Coord_orth ca_pos = frag[subject_res_num][" CA "].pos;
      clipper::Coord_orth ld_pos(next_n, c_pos, ca_pos, 1.8, a, t);
      float db = score_position(ld_pos, xmap_in);

      if (false) {
	 std::cout << "ld_pos: " << ld_pos << " scores " << db << " c.f. " << s->score_per_atom
		   << std::endl;
	 std::ofstream f("debug-pos.txt");
	 f << "set_rotation_centre( " << ld_pos.x() << ", " << ld_pos.y() << ", " << ld_pos.z()
	   << ")\n";
      }

      std::pair<clipper::Coord_orth, float> sp(ld_pos, db);
      s->scored_characteristic_low_density_points.push_back(sp);
      s->atom_point_score -= db;
   } else {
      std::cout << "DEBUG:: oops residue subject_res_num " << subject_res_num << " - No atoms"
		<< std::endl;
   }
}


std::pair<coot::ligand_score_card, coot::minimol::fragment>
coot::residue_by_phi_psi::fit_terminal_residue_generic_trial_inner(int itrial,
								   int offset,
								   int next_residue_seq_num,
								   const std::vector<clipper::Coord_orth> &pos,
								   bool two_residues_flag,
								   bool do_rigid_body_refinement,
								   const clipper::Xmap<float> &xmap_in) const {

   bool debug = false;
   coot::minimol::fragment frag;
   clipper::Coord_orth next_n  = pos[0];
   clipper::Coord_orth next_c  = pos[1];
   clipper::Coord_orth next_ca = pos[2];

   phi_psi_pair p1 = get_phi_psi_by_random();
   phi_psi_pair p2 = get_phi_psi_by_random();

   if (two_residues_flag) { // the majority of times

      if (false) {
	 std::ofstream out("phi-psis.tab", std::ios::app);
	 out << "phi-psi 1 " << p1.phi << " " << p1.psi << " " << p1.tau << std::endl;
	 // out << "phi-psi 2 " << p2.phi << " " << p2.psi << " " << p2.tau << std::endl;
      }

      frag = make_2_res_joining_frag(chain_id, p1, p2,
				     residue_p->GetSeqNum(),
				     offset,  // + or - 1
				     next_n,
				     next_ca,
				     next_c);

      // 	    std::cout << "DEBUG:: back from make_2: checking " << std::endl;
      // 	    frag.check();
      // 	    std::cout << "DEBUG:: (back from make_2) make_2_res_joining_frag returns "
      // 		      << "fragment of size " <<  frag.residues.size()
      // 		      << " residues " << std::endl;

      // 	    std::cout << "Here is some residue information about frag (residues_offset: " <<
      // 		      frag.min_res_no() -1 << "):\n";
      // 	    for (int ii=0; ii<frag.residues.size(); ii++)
      // 	       std::cout << "   residue index " << ii << " has seqnum: "
      // 			 << frag.residues[ii].seqnum << " and "
      // 			 << frag.residues[ii].atoms.size()
      // 			 << " atoms" << std::endl;

      if (debug) {
	 int neighb_seqnum = residue_p->GetSeqNum();

	 clipper::Coord_orth frag_n;
	 clipper::Coord_orth frag_ca;
	 clipper::Coord_orth frag_c;
	 clipper::Coord_orth frag_n_plus_1;
	 if (offset == 1) {
	    frag_n        = frag[neighb_seqnum+1][" N  "].pos;
	    frag_ca       = frag[neighb_seqnum+1][" CA "].pos;
	    frag_c        = frag[neighb_seqnum+1][" C  "].pos;
	    frag_n_plus_1 = frag[neighb_seqnum+2][" N  "].pos;
	 } else {
	    // std::cout << "CHECKME! getting atoms from residues "
	    // << neighb_seqnum-1 << " and " << neighb_seqnum - 2 << std::endl;
	    frag_n        = frag[neighb_seqnum-1][" N  "].pos;
	    frag_ca       = frag[neighb_seqnum-1][" CA "].pos;
	    frag_c        = frag[neighb_seqnum-1][" C  "].pos;
	    frag_n_plus_1 = frag[neighb_seqnum-2][" N  "].pos;
	 }

	 float phi_real = clipper::Util::rad2d(clipper::Coord_orth::torsion(next_c, frag_n, frag_ca, frag_c));
	 float psi_real = clipper::Util::rad2d(clipper::Coord_orth::torsion(frag_n, frag_ca, frag_c, frag_n_plus_1));
	 std::cout <<  "Rama-values:: " << phi_real << " " << psi_real << std::endl;
      } // end debug

   } else {
      frag.addresidue(construct_joining_res(get_phi_psi_by_random(),
					    next_residue_seq_num,
					    next_n,
					    next_ca,
					    next_c), 0);
   }

   std::vector<minimol::atom *> atoms_p = frag.select_atoms_serial();

   // we can't call this because rigid_body_refine_ligand() is not cont
   //
   // if (do_rigid_body_refinement)
   //     rigid_body_refine_ligand(&atoms_p, Xmap()); // tinker with atoms_p

   coot::ligand_score_card s = score_orientation(atoms_p, xmap_in);

   bool do_low_density_fingerprint_positions = true;
   if (do_low_density_fingerprint_positions) {

      if (debug) {
	 std::string file_name = "frag-" + util::int_to_string(itrial) + ".pdb";
	 minimol::molecule m(frag);
	 m.write_file(file_name, 10);
      }
      // the ~HA point (more extended than an H bond distance of course)
      double a = clipper::Util::d2rad(p1.tau-5); // maybe more like 120?
      double t = clipper::Util::d2rad(p1.psi-125);
      // this is for N-terminal addition
      int neighb_seqnum = residue_p->GetSeqNum();
      int subject_res_num = neighb_seqnum-1;
      if (frag[subject_res_num].atoms.size() > 0) {
	 clipper::Coord_orth c_pos  = frag[subject_res_num][" C  "].pos;
	 clipper::Coord_orth ca_pos = frag[subject_res_num][" CA "].pos;
	 clipper::Coord_orth ld_pos(next_n, c_pos, ca_pos, 1.8, a, t);
	 float db = score_position(ld_pos, xmap_in);

	 if (false) {
	    std::cout << "ld_pos: " << ld_pos << " scores " << db << " c.f. " << s.score_per_atom
		      << std::endl;
	    std::ofstream f("debug-pos.txt");
	    f << "set_rotation_centre( " << ld_pos.x() << ", " << ld_pos.y() << ", " << ld_pos.z()
	      << ")\n";
	 }

	 std::pair<clipper::Coord_orth, float> sp(ld_pos, db);
	 s.scored_characteristic_low_density_points.push_back(sp);
	 s.atom_point_score -= db;
      } else {
	 std::cout << "DEBUG:: oops residue subject_res_num " << subject_res_num << " - No atoms"
		   << std::endl;
      }
   }


   if (debug) {
      // DEBUGGING:  Let's write a pdb file for this fragment
      // Then look at them all.  Are they sensibly placed?
      coot::minimol::molecule m_tmp;
      m_tmp.fragments.push_back(frag);
      std::string tmp_filename = "phi-psi-";
      tmp_filename += util::int_to_string(itrial);
      tmp_filename += ".pdb";
      m_tmp.write_file(tmp_filename, 10);
   }

   return std::pair<ligand_score_card, minimol::fragment> (s, frag);

}

// Return a 2 residue fragment that attaches to the starting residue.
// The fragment is either upstream or downstream of the starting residue 
// depending on offset.
// 
coot::minimol::fragment
coot::residue_by_phi_psi::make_2_res_joining_frag(const std::string &chain_id_in,
						  const phi_psi_pair &pp1,
						  const phi_psi_pair &pp2,
						  int seqnum,
						  int offset, // + or - 1
						  const clipper::Coord_orth &next_n_in,
						  const clipper::Coord_orth &next_ca_in,
						  const clipper::Coord_orth &next_c) const {
   bool debug = false;
   minimol::fragment frag(chain_id_in);
   minimol::residue res1;
   minimol::residue res2;

   clipper::Coord_orth next_n   = next_n_in;
   clipper::Coord_orth next_ca = next_ca_in;

   double rand_lim = 0.0;
   next_n += clipper::Coord_orth(rand_lim * (util::random()/float (RAND_MAX) - 0.5),
				 rand_lim * (util::random()/float (RAND_MAX) - 0.5),
				 rand_lim * (util::random()/float (RAND_MAX) - 0.5));
   next_ca += clipper::Coord_orth(0.1 * (util::random()/float (RAND_MAX) - 0.5),
				  0.1 * (util::random()/float (RAND_MAX) - 0.5),
				  0.1 * (util::random()/float (RAND_MAX) - 0.5));
   
   // std::cout << "seed N " << next_n.x() << " " << next_n.y() << " " << next_n.z() << std::endl;

   if (terminus_type == "C" || terminus_type == "MC" || terminus_type == "singleton") { 
      res1 = construct_next_res_from_rama_angles(pp1.phi, pp1.psi, pp1.tau,
						 seqnum + offset,
						 next_n, next_ca, next_c);
      res2 = construct_next_res_from_rama_angles(pp2.phi, pp2.psi, pp2.tau,
						 seqnum + offset + offset,
						 res1[" N  "].pos,
						 res1[" CA "].pos,
						 res1[" C  "].pos);
      
      // now let's tinker with (correct the position of) the O atom of residue 1,
      // now that we have residue 2.
      //
      // double angle =  clipper::Util::d2rad(120.800);  // CA-C-O no 
      double angle = clipper::Util::d2rad(123.0); // N-C-O
      double torsion = clipper::Util::d2rad(0.0);
      clipper::Coord_orth o_pos(res2[" CA "].pos, res2[" N  "].pos, res1[" C  "].pos,
				1.231, angle, torsion);
      for (unsigned int iat=0; iat<res1.atoms.size(); iat++)
	 if (res1[iat].name == " O  ")
	    res1[iat].pos = o_pos;

      // Here we need to tinker with the O position of the *second*
      // residue so that it is not necesarility in the plane of N(n)
      // Ca(n) C(n).

      torsion = 360.0 * float (util::random())/float (RAND_MAX); 
      angle = clipper::Util::d2rad(120.8); // CA-C-O
      clipper::Coord_orth o_pos2(res2[" N  "].pos, res2[" CA "].pos, res2[" C  "].pos,
				 1.231, angle, torsion);
      for (unsigned int iat=0; iat<res2.atoms.size(); iat++)
	 if (res2[iat].name == " O  ") { 
// 	    std::cout << "INFO:: adjusting O of second res from " 
// 		      << res2[iat].pos.format() << " to " << o_pos2.format() << std::endl;;
	    res2[iat].pos = o_pos2;
	 }

   } else { // terminus_type == "N"
      res1 = construct_prev_res_from_rama_angles(pp1.phi, pp1.psi, pp1.tau,
						 seqnum + offset,
						 next_n, next_ca, next_c);
      res2 = construct_prev_res_from_rama_angles(pp2.phi, pp2.psi, pp2.tau,
						 seqnum + offset + offset,
						 res1[" N  "].pos,
						 res1[" CA "].pos,
						 res1[" C  "].pos);
   }

   try { 
      frag.addresidue(res1, 0);
      frag.addresidue(res2, 0);
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR:: make_2_res_joining_frag() " << rte.what() << std::endl;
   }


   int n_atoms1, iseqnum1;
   n_atoms1 = frag[res1.seqnum].atoms.size();
   iseqnum1 = frag[res1.seqnum].seqnum;

   int n_atoms2, iseqnum2;
   n_atoms2 = frag[res2.seqnum].atoms.size();
   iseqnum2 = frag[res2.seqnum].seqnum;
   
   if (debug) {
      std::cout << "DEBUG:: terminus_type == N: frag[" << res2.seqnum << "] has "
		<< n_atoms2 << " atoms and seqnum " << iseqnum2 << "\n";

   
      std::cout << "END of make_2_res_joining_frag frag has " << frag.residues.size()
		<< " residues - it should have 3" << std::endl;

      std::cout << "END of make_2_res_joining_frag: check..." << std::endl;
      frag.check();
   }

   return frag;
} 



// If there are 3 positions on return their order is N, C, CA
// 
std::vector<clipper::Coord_orth>
coot::residue_by_phi_psi::get_connecting_residue_atoms() const {

   // Note, it seems that if there are 2 atoms in a residue with the
   // selected name (which happens if there is an altconf), then
   // GetAtom() returns NULL.  Ugh.
   // 
   
   std::vector<clipper::Coord_orth> pos;

   // Lets look through the atoms of the residue, looking first for
   // altconf "", if not that then "A", if not that, then anything.
   // We do this by assigning backwards, the "A" gets overwritten ""
   // if it exists (actually, this is a pathological case, "" and "A"
   // should not exist in the same residue for the same atom name.
   // 
   // But to do this we need flag to see if each of the
   // search_atom_name were found.

   std::string search_atom_name[3];
   search_atom_name[0] = " N  ";
   search_atom_name[1] = " C  ";
   search_atom_name[2] = " CA ";

   short int search_atom_flag[3];
   for(int i=0; i<3; i++)
      search_atom_flag[i] = 0; // initially not found.


   std::string search_altconf[3];
   search_altconf[0] = "A";
   search_altconf[1] = "";

   clipper::Coord_orth found_atom[3];

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int isearch=0; isearch<3; isearch++) {
      for(int ialtconf=0; ialtconf<3; ialtconf++) {
	 for (int i=0; i<nResidueAtoms; i++) { 
	    if (std::string(residue_atoms[i]->name) == search_atom_name[isearch]) { 
	       if (std::string(residue_atoms[i]->altLoc) == search_altconf[ialtconf]) {
		  found_atom[isearch] = clipper::Coord_orth(residue_atoms[i]->x,
							    residue_atoms[i]->y,
							    residue_atoms[i]->z);
		  // std::cout << "found a " << residue_atoms[i]->name << std::endl;
		  search_atom_flag[isearch] = 1;
	       }
	    }
	 } 
      }
   } 
   


   //    if (pos.size() != 3) { 
   if (! search_atom_flag[0] || ! search_atom_flag[1] || ! search_atom_flag[2]) { 
      std::cout << "ERROR: missing atoms in get_connecting_residue_atoms: we found "
		<< pos.size() << " atoms\n";
      if (residue_p) {
	 mmdb::PPAtom residue_atoms_l;
	 int nResidueAtoms_l;
	 residue_p->GetAtomTable(residue_atoms_l, nResidueAtoms_l);
	 std::cout << "residue has " << nResidueAtoms_l << " atoms " << std::endl;
	 for (int i=0; i<nResidueAtoms_l; i++) { 
	    std::cout << i << " " << residue_atoms_l[i] << std::endl;
	 } 
      } else { 
	 std::cout << "ERROR: missing atoms in get_connecting_residue_atoms: we have "
		   << "NULL residue_p\n";
      }
   } else { 
      // All hunkey dorey
      pos.push_back(found_atom[0]);
      pos.push_back(found_atom[1]);
      pos.push_back(found_atom[2]);
   }
   return pos;
}


coot::minimol::residue
coot::residue_by_phi_psi::construct_next_res_from_rama_angles(float phi, float psi, float tau,
						     int seqno,
						     const clipper::Coord_orth &previous_n,
						     const clipper::Coord_orth &previous_ca,
						     const clipper::Coord_orth &previous_c) const {


   coot::minimol::residue mres(seqno);
   mres.name = residue_type;

   double angle, torsion;
	 
   angle = clipper::Util::d2rad(116.200); // Ca-C-N
   torsion = clipper::Util::d2rad(psi);
   clipper::Coord_orth n_pos(previous_n, previous_ca, previous_c, 
			     1.329, angle, torsion); // C-N bond

   angle = clipper::Util::d2rad(121.700); // C-N-Ca
   torsion = clipper::Util::d2rad(180.0); 
   clipper::Coord_orth ca_pos(previous_ca, previous_c, n_pos,
			      1.458, angle, torsion); // N-CA bond

   angle = clipper::Util::d2rad(tau); // N-CA-C  // was 111.200
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
coot::residue_by_phi_psi::construct_prev_res_from_rama_angles(float phi, float psi, float tau,
							      int seqno,
							      const clipper::Coord_orth &next_n,
							      const clipper::Coord_orth &next_ca,
							      const clipper::Coord_orth &next_c) const {

   coot::minimol::residue mres(seqno);
   mres.name = residue_type;

   double angle, torsion;

   // std::cout << "constructing previous residue trial:  next n:" << next_n.format() << std::endl;

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
   angle = clipper::Util::d2rad(tau); // N-Ca-C
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


   // 20170925 there is no CBs. Why? The get added in execute_add_terminal_residue()
   // add_side_chain_to_terminal_res().
   // I've added it here now.

   // CB
   angle = clipper::Util::d2rad(tau);
   torsion = clipper::Util::d2rad(psi+123.4);
   clipper::Coord_orth cb_pos(next_n, c_pos, ca_pos,
			     1.52, angle, torsion); // Ca-Cb


   // ---------------------------

   mres.addatom(coot::minimol::atom(" N  ", " N", n_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" C  ", " C", c_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" CA ", " C", ca_pos, "", b_factor));
   mres.addatom(coot::minimol::atom(" O  ", " O", o_pos,  "", b_factor));
   mres.addatom(coot::minimol::atom(" CB ", " C", cb_pos, "", b_factor));

   return mres;
}

coot::minimol::residue
coot::residue_by_phi_psi::construct_joining_res(const phi_psi_pair &pp,
						int seqno,
						const clipper::Coord_orth &next_n,
						const clipper::Coord_orth &next_ca,
						const clipper::Coord_orth &next_c) const { 

   if (terminus_type == "C" || terminus_type == "MC" || terminus_type == "singleton") 
      return construct_next_res_from_rama_angles(pp.phi, pp.psi, pp.tau,
						 seqno, next_n, next_ca, next_c);
   else // terminus_type == "N" or "MN"
      return construct_prev_res_from_rama_angles(pp.phi, pp.psi, pp.tau,
						 seqno, next_n, next_ca, next_c);
}
						
void 
coot::residue_by_phi_psi::init_phi_psi_plot() { 

   if (residue_type == "GLY") {
      rama.init(clipper::Ramachandran::Gly5);
   } else { 
      if (residue_type == "PRO") {
	 rama.init(clipper::Ramachandran::Pro5);
      } else { 
	 rama.init(clipper::Ramachandran::NonGlyPro5);
      }
   }
   
   rama_max = 0.0;
   
   for (float phi=0.0; phi<360.0; phi+=3.0) { 
      for (float psi=0.0; psi<360.0; psi+=3.0) { 
	 float v = rama.probability(clipper::Util::d2rad(phi), 
				    clipper::Util::d2rad(psi));
	 if (v > rama_max) 
	    rama_max = v;
      }
   }
}

coot::phi_psi_pair
coot::residue_by_phi_psi::get_phi_psi_by_random() const { 
   
   float phi, psi, prob, r;
   
   for (;;) {
      phi = 360.0 * float (coot::util::random())/float (RAND_MAX);
      psi = 360.0 * float (coot::util::random())/float (RAND_MAX);

      r = rama_max * coot::util::random()/ float (RAND_MAX);
      prob = rama.probability(clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi));
      
      if (prob > r) break;
   }
   float u1 = static_cast<float>(coot::util::random())/static_cast<float>(RAND_MAX);
   float u2 = 2.0 * u1 - 1.0; // -1 -> 1.
   float tau = 111.2 + u2 * 2.0; // 2.0 needs investigation
   return phi_psi_pair(phi, psi, tau);
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

   // std::cout << "constructing previous residue trial:  next n:" << next_n.format() << std::endl;

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
