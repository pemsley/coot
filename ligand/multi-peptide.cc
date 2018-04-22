/* ligand/multi-peptide.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
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
 * 02110-1301, USA
 */

#include "geometry/residue-and-atom-specs.hh"
#include "mini-mol/mini-mol-utils.hh"
#include "ideal/simple-restraint.hh"
#include "multi-peptide.hh"
#include "residue_by_phi_psi.hh"
#include "coot-utils/coot-map-utils.hh"
#include "analysis/stats.hh"


coot::minimol::fragment
coot::multi_build_N_terminal_ALA(mmdb::Residue *res_p,
				 mmdb::Residue *up_or_down_stream_neighbour_p, // depends on offset
				 const std::string &chain_id,
				 float b_factor_in,
				 int n_trials,
				 const coot::protein_geometry &geom,
				 const clipper::Xmap<float> &xmap,
				 std::pair<float, float> mv,
				 bool debug_trials) {

   return multi_build_terminal_ALA(-1, res_p, up_or_down_stream_neighbour_p,
				   chain_id, b_factor_in, n_trials, geom, xmap, mv, debug_trials);

}

coot::minimol::fragment
coot::multi_build_C_terminal_ALA(mmdb::Residue *res_p,
				 mmdb::Residue *upstream_neighb_p,
				 const std::string &chain_id,
				 float b_factor_in,
				 int n_trials,
				 const coot::protein_geometry &geom,
				 const clipper::Xmap<float> &xmap,
				 std::pair<float, float> mv,
				 bool debug_trials) {

   return multi_build_terminal_ALA(1, res_p, upstream_neighb_p, chain_id, b_factor_in, n_trials, geom, xmap, mv, debug_trials);

}

// The new algorithm is quite a bit different from the previous one - rather than edit that,
// I'll just start from scratch (there are common elements).
//
coot::minimol::fragment
coot::multi_build_terminal_residue_addition_forwards_2018(mmdb::Residue *res_p,
							  mmdb::Residue *upstream_neighbour_p,
							  const std::string &chain_id,
							  float b_factor_in,
							  int n_trials,
							  const coot::protein_geometry &geom,
							  const clipper::Xmap<float> &xmap,
							  std::pair<float, float> mv,
							  bool debug_trials) {

   ctpl::thread_pool thread_pool(coot::get_max_number_of_threads());

   minimol::fragment many_residues; // add residue to this, and return it
   many_residues.addresidue(res_p, false);

   if (true) { // debug
      std::string file_name = "frag-with-start-residue.pdb";
      minimol::molecule mmm(many_residues);
      mmm.write_file(file_name, 10);
   }

   // setup up initial conditions for the while loop
   bool prev_happy_fit = true;
   bool happy_fit = true; // starting
   int icount = 0;
   mmdb::Residue *res_prev_p = 0; // whatever was res_p the prevous round
   minimol::residue res_prev; // what was res on the previous round.

   while (happy_fit) {

      if (false) { // debug
	 std::string file_name = "while-new-";
	 file_name += coot::util::int_to_string(res_p->GetSeqNum()) + std::string(".pdb");
	 minimol::molecule mmm(many_residues);
	 mmm.write_file(file_name, 10);
      }

      // fill res, the first of the 2 residues fitted by residue_by_phi_psi
      minimol::residue res;

      std::string residue_type = "ALA"; // for the future, if residue is a bad fit or CB is a bad fit
                                        // then try fitting a GLY and compare statistics, choosing
                                        // a GLY if suitable.  Similar for PRO - where we deliberately
                                        // downweight solutions that have density where PRO CD is.

      residue_by_phi_psi addres("C", res_p, chain_id, residue_type, 20);
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      unsigned int n_threads = coot::get_max_number_of_threads();
      if (n_threads >= 1)
	    addres.thread_pool(&thread_pool, n_threads);
#endif
      if (debug_trials)
	 addres.write_trial_pdbs();
      addres.set_upstream_neighbour(res_prev_p);
      addres.import_map_from(xmap);
      minimol::fragment f = addres.best_fit_phi_psi(n_trials, 1); // offset = 1
      // add frag to many_residues, and set res
      bool found_res = false;
      for (int ires=f.min_res_no(); ires<=f.max_residue_number(); ires++) {
	 if (f[ires].atoms.size() > 0) {
	    if (! found_res) {
	       res = f[ires];
	       found_res = true;
	    }
	    many_residues.addresidue(f[ires], false);
	 }
      }

      std::cout << "== res.seqnum " << res.seqnum << std::endl;

      // now we have 2 residues in fragment f. We need to move the O of our current res_p
      // to be consistent with those 2 new residues
      {
	 minimol::atom C_res_prev = many_residues[res_p->GetSeqNum()][" C  "];
	 // std::cout << "DEBUG:: here with C_res_prev " << C_res_prev << std::endl;
	 double angle   = clipper::Util::d2rad(123.0); // N-C-O
	 double torsion = clipper::Util::d2rad(0.0);
	 clipper::Coord_orth new_o_pos_for_current_res_new(res[" CA "].pos,
							   res[" N  "].pos,
							   C_res_prev.pos, 1.231, angle, torsion);
	 minimol::residue &res_current = many_residues[res_p->GetSeqNum()];
	 if (false) { // debug res_current
	    std::cout << "debug:: res_current " << res_current << " has "
		      << res_current.atoms.size() << " atoms " << std::endl;
	    for (std::size_t iat=0; iat<res_current.atoms.size(); iat++) {
	       std::cout << "    " << iat << " " << res_current.atoms[iat] << std::endl;
	    }
	 }
	 minimol::atom &o_at = many_residues[res_p->GetSeqNum()].at(" O  ");
	 o_at.pos = new_o_pos_for_current_res_new;
	 if (false)
	    std::cout << "    O of residue " << res_p->GetSeqNum() << "\n    "
		      << " using CA " << res[" CA "].pos.format()  << "\n    "
		      << " using N "  << res[" N  "].pos.format()  << "\n    "
		      << " using C "  << C_res_prev.pos.format()   << "\n    "
		      << " updated to "
		      << many_residues[res_p->GetSeqNum()].at(" O  ").pos.format()
		      << " cf " << new_o_pos_for_current_res_new.format() << std::endl;
      }

      std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(res);
      if (cbeta_info.first) {
	 res.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 20.0);
	 many_residues[res.seqnum].addatom(" CB ", " C", cbeta_info.second, "", 1.0, 30.0);
      }

      // res is a copy of the residue in many_residues

      // refine res and its neighbour and one more, update the atoms of many_residues
      //
      int offset = 1;
      refine_end(&many_residues, res.seqnum, offset, geom, xmap);

      // copy res after refinement.  Maybe I can use a reference?

      res = many_residues[res.seqnum];

      // after refinement we need to ask if res is well-fitting.
      // So res needs to be properly filled.
      bool this_happy_fit = does_residue_fit(res, xmap, mv);

      if (! this_happy_fit) {
	 // try stripping off the CB and making the residue a GLY?
	 // and refining again
      }

      if (! this_happy_fit) {
	 // exit on just one failure for the moment (testing).
	 // Seems reasonable.
	 happy_fit = false;
      }

      if (false) { // debug
	 std::string file_name = "post-refined-";
	 file_name += util::int_to_string(res.seqnum) + std::string(".pdb");
	 minimol::molecule mmm(many_residues);
	 mmm.write_file(file_name, 10);
      }

      many_residues.residues.pop_back();

      if (crashing_into_self(many_residues, res.seqnum, offset)) {
	 happy_fit = false;
      }

      if (false) { // debug
	 std::string file_name = "post-chop-last-residue-";
	 file_name += util::int_to_string(res.seqnum) + std::string(".pdb");
	 minimol::molecule mmm(many_residues);
	 mmm.write_file(file_name, 10);
      }

      // res is the first new residue

      if (! this_happy_fit) {
	 std::pair<bool, minimol::residue> recover =
	    try_to_recover_from_bad_fit_forwards(res_p, res_prev_p, chain_id, n_trials, geom, xmap, mv);
	 if (recover.first) {
	    this_happy_fit = true;
	    happy_fit = true;
	    res = recover.second;
	    many_residues[res.seqnum] = res;

	    // this is necessary if I refine before adding res to many_residues
	    // many_residues.residues.pop_back();

	    // update the O pos of the previous residue
	    //
	    // make this a function update_O_position_in_prev_residue(res_p, res);
	    update_O_position_in_prev_residue(res_p, &many_residues, res);

	    // what does this do? Optimize the O pos?
	    // refine_end(&many_residues, res.seqnum, offset, geom, xmap);

	    std::cout << "post recovery" << std::endl;
	 }
      }

      if (! this_happy_fit && ! prev_happy_fit) {
	 happy_fit = false;
	 std::cout << "INFO:: 2 unhappy fits - stopping now " << std::endl;
      }

      if (res.is_empty()) {
	 happy_fit = false;
      } else {
	 if (happy_fit) {
	    icount++;

	    // set res_p to be the newly constructed residue
	    res_prev   = res;
	    res_prev_p = res_p; // save to set upstream_C
	    res_p = res.make_residue();
	    prev_happy_fit = this_happy_fit;
	 }
      }
   }

   return many_residues;
}

void
coot::update_O_position_in_prev_residue(mmdb::Residue *res_p,
					coot::minimol::fragment *many_residues,
					const coot::minimol::residue &res) {

   minimol::atom C_res_prev  = many_residues->at(res_p->GetSeqNum())[" C  "];
   minimol::atom CA_res_prev = many_residues->at(res_p->GetSeqNum())[" CA "];
   double angle   = clipper::Util::d2rad(123.0); // N-C-O
   double tors_deg = 0.0; // O is trans to the CA of the next residue
   // unless peptide is cis.
   double tors_peptide = clipper::Coord_orth::torsion(CA_res_prev.pos, C_res_prev.pos,
						      res[" N  "].pos, res[" CA "].pos);
   // cis or trans
   if (std::abs(tors_peptide) < M_PI_2) tors_deg = 180.0;
   double torsion = clipper::Util::d2rad(tors_deg);
   clipper::Coord_orth new_o_pos_for_current_res_new(res[" CA "].pos,
						     res[" N  "].pos,
						     C_res_prev.pos, 1.231, angle, torsion);
   minimol::atom &o_at = many_residues->at(res_p->GetSeqNum()).at(" O  ");
   o_at.pos = new_o_pos_for_current_res_new;
}


// try again with more trials
// try again with PRO
// try again with GLY
std::pair<bool, coot::minimol::residue>
coot::try_to_recover_from_bad_fit_forwards(mmdb::Residue *res_p, mmdb::Residue *res_prev_p,
					   const std::string &chain_id, int n_trials,
					   const coot::protein_geometry &geom,
					   const clipper::Xmap<float> &xmap,
					   std::pair<float, float> mv) {

   std::cout << "try to recover.. " << std::endl;

   std::pair<bool, coot::minimol::residue> r;
   r.first = false;

   std::string residue_type = "ALA";
   std::string dir = "C";
   residue_by_phi_psi rpp(dir, res_p, chain_id, residue_type, 20);
   rpp.set_upstream_neighbour(res_prev_p);
   rpp.import_map_from(xmap);

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   // Without adding the thread pool, we go through a different code path in fit_terminal_residue_generic()
   // and that can be useful for debugging.
   //
   unsigned int n_threads_max = coot::get_max_number_of_threads();
   ctpl::thread_pool thread_pool(n_threads_max); // pass the tread pool
   if (n_threads_max >= 1)
      rpp.thread_pool(&thread_pool, n_threads_max);
#endif

   minimol::fragment f = rpp.best_fit_phi_psi(n_trials * 8, 1);

   // we don't refine here - is that a good idea?
   //
   // No, it's not. Add refinement.

   int new_res_seqnum = res_p->GetSeqNum()+1;
   refine_end(&f, new_res_seqnum, 1, geom, xmap);

   const minimol::residue &res = f[new_res_seqnum];
   bool liberal_fit = true;
   bool this_happy_fit = does_residue_fit(res, xmap, mv);
   if (this_happy_fit) {

      std::cout << "... recovered with more trials " << std::endl;
      r.first = true;
      r.second = res;
   } else {

      std::cout << "try to recover as a PRO..." << std::endl;
      // try as a PRO
      // note to self: only the first residue of the pair should be a PRO
      residue_by_phi_psi rpp(dir, res_p, chain_id, "PRO", 30.0);
      rpp.set_upstream_neighbour(res_prev_p);
      rpp.import_map_from(xmap);
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      if (n_threads_max >= 1)
	 rpp.thread_pool(&thread_pool, n_threads_max);
#endif

      minimol::fragment f = rpp.best_fit_phi_psi(800, 1);

      // we don't refine here - is that a good idea?

      minimol::residue &res = f[res_p->GetSeqNum()+1];
      bool this_happy_fit = does_residue_fit(res, xmap, mv);
      if (this_happy_fit) {
	 std::cout << "... recovered as a PRO " << std::endl;
	 // add CG and CD to res
	 double bl_cg   = 1.515;
	 double bl_cd   = 1.508;
	 double ang_cg  = 103.5;
	 double ang_cd  = 104.5;
	 double tors_cg =  21.5;
	 double tors_cd = -30.5;
	 double b_factor = 37.7; // Hmm
	 std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(res);
	 if (cbeta_info.first) {
	    res.addatom(" CB ", " C", cbeta_info.second, "", 1.0, b_factor);
	    clipper::Coord_orth cg_pos = clipper::Coord_orth(res[" N  "].pos,
							     res[" CA "].pos,
							     res[" CB "].pos, bl_cg,
							     clipper::Util::d2rad(ang_cg),
							     clipper::Util::d2rad(tors_cg));
	    clipper::Coord_orth cd_pos = clipper::Coord_orth(res[" CA "].pos,
							     res[" CB "].pos,
							     cg_pos, bl_cd,
							     clipper::Util::d2rad(ang_cd),
							     clipper::Util::d2rad(tors_cd));
	    res.addatom(coot::minimol::atom(" CG ", " C", cg_pos, "", b_factor));
	    res.addatom(coot::minimol::atom(" CD ", " C", cd_pos, "", b_factor));
	    res.name = "PRO";
	 }

	 r.first = true;
	 r.second = res;
      } else {

	 // Try as a GLY?

      }
   }

   std::cout << "debug:: recover status: " << r.first << std::endl;

   return r;
}


// up_or_down_stream_neighbour_p is allowed to be null
//
coot::minimol::fragment
coot::multi_build_terminal_ALA(int offset, // direction
			       mmdb::Residue *res_p,
			       mmdb::Residue *up_or_down_stream_neighbour_p, // depends on offset
			       const std::string &chain_id,
			       float b_factor_in,
			       int n_trials,
			       const coot::protein_geometry &geom,
			       const clipper::Xmap<float> &xmap,
			       std::pair<float, float> mv,
			       bool debug_trials) {

   
   // res_p might not be a member of hierarchy, so pass the chain id too.

   minimol::fragment many_residues;

   std::string dir = "N";
   if (offset == 1)
      dir = "C";

   bool prev_happy_fit = true;
   bool happy_fit = true; // starting
   int icount = 0;

   mmdb::Residue *res_prev = 0; // whatever was res_p the prevous round

   while (happy_fit) {

      // res from fragment f is what we really want
      minimol::residue res;

      std::string residue_type = "ALA";

      residue_by_phi_psi rpp(dir, res_p, chain_id, residue_type, 20);
      if (debug_trials)
	 rpp.write_trial_pdbs();
      if (offset ==  1) rpp.set_upstream_neighbour(res_prev);
      if (offset == -1) rpp.set_downstream_neighbour(res_prev);

      rpp.import_map_from(xmap);

      minimol::fragment f = rpp.best_fit_phi_psi(n_trials, offset);
      if (offset == -1) {
	 for (int ires=f.max_residue_number(); ires>=f.min_res_no(); ires--) {
	    if (f[ires].atoms.size() > 0) {
	       res = f[ires];
	       break;
	    }
	 }
      } else {
	 // for C-terminal builds, get the first residue.
	 for (int ires=f.min_res_no(); ires<=f.max_residue_number(); ires++) {
	    if (f[ires].atoms.size() > 0) {
	       res = f[ires];
	       break;
	    }
	 }
      }

      if (res.is_empty()) {
	 happy_fit = false;
      } else {

	 std::cout << "   debug::multi_build_terminal_ALA() with offset "
		   << offset << " was passed residue seqnum " << res_p->GetSeqNum()
		   << " and found a newly-built residue with seqnum " << res.seqnum << std::endl;

	 std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(res);
	 if (cbeta_info.first) { 
	    res.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 20.0);
	    res.name = "ALA";
	 }

	 if (offset == -1)
	    res.seqnum = res_p->GetSeqNum() - 1;
	 if (offset ==  1)
	    res.seqnum = res_p->GetSeqNum() + 1;
	    
	 bool this_happy_fit = does_residue_fit(res, xmap, mv);

	 // std::cout << "Here with this_happy_fit " << this_happy_fit << " for residue " << res << std::endl;

 	 if (! this_happy_fit && ! prev_happy_fit) {
 	    happy_fit = false;
	    std::cout << "   info:: 2 unhappy fits - stopping now " << std::endl;
	 }

	 if (happy_fit) { 
	    icount++;
	    // add res to mol and refine
	    many_residues.addresidue(res, 20);

	    if (true) {
	       // debug:  Let's write a pdb file for this fragment
	       coot::minimol::molecule m_tmp;
	       m_tmp.fragments.push_back(many_residues);
	       std::string tmp_filename = "frag-from-phi-psi-";
	       tmp_filename += util::int_to_string(icount); 
	       tmp_filename += ".pdb";
	       m_tmp.write_file(tmp_filename, 10);
	    }

	    // refine res and its neighbour, update the atoms of many_residues
	    //
	    refine_end(&many_residues, res.seqnum, offset, geom, xmap);

	    // ------------- For next round ----------------
	    //
	    // set res_p to be the newly constructed residue
	    res_prev = res_p; // save to set upstream_C or downstream_N for a residue_by_phi_psi
	    res_p = res.make_residue();
	    prev_happy_fit = this_happy_fit;
	 }
      }
   }
   return many_residues;
}

// move the atoms (the end residue and its neighbour) of many_residue
// 
// offset tells us which direction the new residue was added (-1 means
// new N terminus)
//
// seqnum here is the residue number of the first addd residue
//
void
coot::refine_end(coot::minimol::fragment *many_residues,
		 int seqnum, int offset,
		 const coot::protein_geometry &geom,
		 const clipper::Xmap<float> &xmap) {

   mmdb::Manager *mol = new mmdb::Manager; // d
   mmdb::Model *model_p = new mmdb::Model;
   int ires_new = -1;

   if (offset == -1) {
      ires_new = many_residues->first_residue(); // with atoms
   } else {
      // the resno of the last residue in the fragment
      ires_new = many_residues->max_residue_number();
   }

   mmdb::Chain *chain_p = many_residues->make_chain();

   model_p->AddChain(chain_p);
   mol->AddModel(model_p);
 
   mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   mol->FinishStructEdit();

   if (false) {
      std::string pdb_refine_file_name = "refine-in-";
      pdb_refine_file_name += coot::util::int_to_string(seqnum);
      pdb_refine_file_name += ".pdb";
      mol->WritePDBASCII(pdb_refine_file_name.c_str());
   }

   std::vector<std::pair<bool, mmdb::Residue *> > residues;
   std::vector<mmdb::Residue *> moving_residues; 

   int nres = chain_p->GetNumberOfResidues();

   for (int ires=0; ires<nres; ires++) {
      mmdb::Residue *residue_p = chain_p->GetResidue(ires);
      std::pair<bool, mmdb::Residue *> p(false, residue_p);
      int res_seq_num = residue_p->GetSeqNum();

      if (res_seq_num == ires_new) { 
	 moving_residues.push_back(residue_p);
	 residues.push_back(p);
      }
      if (offset == -1) {
	 if (res_seq_num == (ires_new + 1)) { // the next residue after the N-terminal residue
	    moving_residues.push_back(residue_p);
	    residues.push_back(p);
	 }
	 if (res_seq_num == (ires_new + 2)) {
	    std::pair<bool, mmdb::Residue *> fixed_p(false, residue_p);
	    residues.push_back(fixed_p);
	 }
      }
      if (offset == 1) {
	 if (res_seq_num == (ires_new - 1)) { // the residue before the C-terminal residue
	    moving_residues.push_back(residue_p);
	    residues.push_back(p);
	 }
	 if (res_seq_num == (ires_new - 2)) {
	    std::pair<bool, mmdb::Residue *> fixed_p(false, residue_p);
	    // 20180414-PE I think I want the residue that the new residue
	    // attaches to to be moving too!
	    moving_residues.push_back(residue_p);
	    residues.push_back(fixed_p);
	 }
      }
   }

   std::vector<mmdb::Link> links;
   std::vector<atom_spec_t> fixed_atom_specs;
   restraint_usage_Flags flags = TYPICAL_RESTRAINTS;

   if (false) { // debug
      std::cout << "refining with " << residues.size() << " residues " << std::endl;
      for (std::size_t ires=0; ires<residues.size(); ires++) {
	 std::cout << "    " << coot::residue_spec_t(residues[ires].second) << std::endl;
      }
   }
   restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, xmap);

   // Does this make things slower? (seems so, try passing the thread pool
   // or test how long it takes to create and add).
   //
// #ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
//    ctpl::thread_pool thread_pool(coot::get_max_number_of_threads());
//    restraints.thread_pool(&thread_pool, get_max_number_of_threads());
// #endif

   pseudo_restraint_bond_type pseudos = NO_PSEUDO_BONDS;
   bool do_internal_torsions = false;
   float weight = 60;

   restraints.set_quiet_reporting();
   restraints.add_map(weight);
   bool do_trans_peptide_restraints = true;
   do_trans_peptide_restraints = false;
   int imol = 0;
   restraints.make_restraints(imol, geom, flags, do_internal_torsions, do_internal_torsions, 0, 0, pseudos);
   restraints.minimize(flags);

   if (false) {
      std::string pdb_refine_file_name = "refine-out-";
      pdb_refine_file_name += coot::util::int_to_string(seqnum);
      pdb_refine_file_name += ".pdb";
      mol->WritePDBASCII(pdb_refine_file_name.c_str());
   }

   for (unsigned int ii=0; ii<moving_residues.size(); ii++) {
      int seqnum_l = moving_residues[ii]->GetSeqNum();
      (*many_residues)[seqnum_l].update_positions_from(moving_residues[ii]);
   }

   delete mol;
} 

bool
coot::does_residue_fit(const coot::minimol::residue &res, const clipper::Xmap<float> &xmap,
		       std::pair<float, float> mv) {

   bool r = true;
   double z_crit = 1.3;

   float mean = mv.first;
   float rmsd = sqrtf(mv.second);

   std::vector<double> rho(res.atoms.size());
   for (unsigned int iat=0; iat<res.atoms.size(); iat++) { 
      float d = util::density_at_point(xmap, res.atoms[iat].pos);
      rho[iat] = d;
   }

   stats::single st(rho);
   for (unsigned int iat=0; iat<res.atoms.size(); iat++) {
      if (rho[iat] < (mean + rmsd * z_crit)) {

	 if (res.atoms[iat].name != " CB ") {  // PDBv3 FIXME
	    std::cout << "   low density for atom residue: " << res.seqnum
		      << " atom: " << res.atoms[iat].name
		      << rho[iat] << " vs " << mean <<  " + " << rmsd << " * " << z_crit << " at "
		      << res.atoms[iat].pos.format()
		      << std::endl;
	    r = false;
	    break;
	 }
      } 
   }

   return r;
} 

bool
coot::crashing_into_self(const minimol::fragment &many_residues, int seqnum, int offset) {

   bool status = false;

   std::vector<clipper::Coord_orth> positions_on_trial;
   if (offset == 1) {
      // const minimol::residue &last_res = many_residues.residues.back(); // interesting but not useful here
      const minimol::residue &last_res = many_residues[seqnum];
      for (std::size_t i=0; i<last_res.atoms.size(); i++) {
	 positions_on_trial.push_back(last_res.atoms[i].pos);
      }
   }
   std::set<int> dont_check_these;
   if (offset == 1) {
      dont_check_these.insert(seqnum);
      dont_check_these.insert(seqnum+1);
   } else {
      dont_check_these.insert(seqnum);
      dont_check_these.insert(seqnum-1);
   }

   unsigned int n_close = 0;
   for (int i=many_residues.min_res_no(); i<=many_residues.max_residue_number(); i++) {
      if (dont_check_these.find(i) == dont_check_these.end()) {
	 for (std::size_t iat=0; iat<many_residues[i].atoms.size(); iat++) {
	    for (std::size_t j=0; j<positions_on_trial.size(); j++) {
	       double l = clipper::Coord_orth::length(positions_on_trial[j], many_residues[i].atoms[iat].pos);
	       if (l < 1.8) { // good value?
		  n_close++;
	       }
	    }
	 }
      }
   }

   if (n_close > 1)
      status = true; // is clashing with self.

   return status;
}
