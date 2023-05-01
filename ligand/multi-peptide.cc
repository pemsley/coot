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

#ifdef HAVE_CXX_THREAD
#include <future>
#endif // HAVE_CXX_THREAD

#include "ideal/simple-restraint.hh"
#include "multi-peptide.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "mini-mol/mini-mol-utils.hh"
#include "residue_by_phi_psi.hh"
#include "coot-utils/coot-map-utils.hh"
#include "analysis/stats.hh"
#include "cootaneer/cootaneer-sequence.h"

#include "trace.hh"

void
coot::multi_build_terminal_residue_addition::setup_standard_residues_mol() {

   standard_residues_mol = 0;
   std::string standard_env_dir = "COOT_STANDARD_RESIDUES";
   const char *filename_str = getenv(standard_env_dir.c_str());
   std::string filename;
   if (! filename_str) {
      std::string standard_residues_file_name = package_data_dir(); // check COOT_PREFIX
      standard_residues_file_name += "/";
      standard_residues_file_name += "standard-residues.pdb";
      filename = standard_residues_file_name;
   } else {
      filename = std::string(filename_str);
   }

   if (file_exists(filename)) {

      mmdb::Manager *mol = new mmdb::Manager;
      mmdb::ERROR_CODE err = mol->ReadCoorFile(filename.c_str());
      if (err) {
	 std::cout << "There was an error reading " << filename.c_str() << ". \n";
	 std::cout << "ERROR " << err << " READ: "
		   << mmdb::GetErrorDescription(err) << std::endl;
	 delete mol;
      } else {
	 standard_residues_mol = mol;
      }
   }
}

void
coot::multi_build_terminal_residue_addition::setup_symms() {

   int symm_shift_range = 2; // 1

   const clipper::Spacegroup &spg = xmap.spacegroup();
   const clipper::Cell &cell = xmap.cell();
   unsigned int ns = spg.num_symops();
   for (int x_shift = -symm_shift_range; x_shift<(1+symm_shift_range); x_shift++) {
      for (int y_shift = -symm_shift_range; y_shift<(1+symm_shift_range); y_shift++) {
	 for (int z_shift = -symm_shift_range; z_shift<(1+symm_shift_range); z_shift++) {
	    clipper::Coord_frac cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift);
	    for (unsigned int is = 0; is < ns; is++) {
	       const clipper::Symop &symop = spg.symop(is);
	       clipper::RTop_frac rtf(symop.rot(), symop.trn() + cell_shift);
	       clipper::RTop_orth rtop = symop.rtop_orth(cell);
	       symms.push_back(rtop);
	    }
	 }
      }
   }
}

void
coot::multi_build_terminal_residue_addition::init_no_go() {

   // we can't use a bool as the type, don't try.
   no_go.init(xmap.spacegroup(),
	      xmap.cell(),
	      xmap.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = no_go.first(); !ix.last(); ix.next() )
      no_go[ix] = 0; // not already traced.
}

void
coot::multi_build_terminal_residue_addition::start_from_map(const coot::protein_geometry &geom) {

   bool debugging = false;
   trace t(xmap);
   std::vector<minimol::fragment> seeds = t.make_seeds();
   std::cout << "INFO:: Made " << seeds.size() << " seeds " << std::endl;

   //kludge for testing backwards building

   // if (seeds.size() > 0)
   // seeds.resize(4);

   // This (and the use of set_cell() and set_spacegroup() later) is inelegant.
   float acell[6];
   acell[0] = xmap.cell().descr().a();
   acell[1] = xmap.cell().descr().b();
   acell[2] = xmap.cell().descr().c();
   acell[3] = clipper::Util::rad2d(xmap.cell().descr().alpha());
   acell[4] = clipper::Util::rad2d(xmap.cell().descr().beta());
   acell[5] = clipper::Util::rad2d(xmap.cell().descr().gamma());
   std::string spacegroup_str_hm = xmap.spacegroup().symbol_hm();

   // We need to worry about fragment store going out of scope.
   // How about: wait for this async to terminate before finishing this function
   //
   // std::vector<std::pair<std::string, std::string> > sequences;
   // std::string k = "SMNPPPPET SNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKI IKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMA EALEKLFLQKINELPTEE";
   // sequences.push_back(std::pair<std::string, std::string> ("A", k));

   // This can be a thread (with later thread.join()), it doesn't need
   // to be a future and async.
   // 201806016-PE
   // In fact a thread is better/necessary - a thread will run, but an async need not.

   //    std::future<double> store_score(std::async(store_manager, std::ref(fragment_store),
   // 					      std::ref(store_lock), std::cref(xmap), std::cref(sequences)));

   std::thread store_thread(store_manager, std::ref(fragment_store),
			    std::ref(store_lock), std::cref(xmap), std::cref(sequences));

   std::cout << "INFO:: sequences: " << sequences.size() << " " << std::endl;
   for (std::size_t ii=0; ii<sequences.size(); ii++) {
      std::cout << sequences[ii].first << "  " << sequences[ii].second << std::endl;
   }

   unsigned int n_seeds_max = 10; // seeds.size();

   for (std::size_t iseed=0; iseed<n_seeds_max; iseed++) {
      minimol::molecule mmm(seeds[iseed]);
      mmdb::Manager *mol = mmm.pcmmdbmanager();
      mmdb::Residue *r = util::get_residue(residue_spec_t("A", 99, ""), mol);
      mmdb::Residue *r_prev = util::get_previous_residue(residue_spec_t(r), mol);
      float b_fact = 30.0;
      int n_trials = 20010; // 20000 is a reasonable minimal number
      n_trials = 100;

      minimol::fragment fC =
         forwards_2018(iseed, r, r_prev, "A", b_fact, n_trials, geom, xmap, mv, debugging);
      if (fC.n_filled_residues() > 2) { // needs tweaking
         int build_dir = 1;
         add_to_fragment_store(fC, build_dir);
      }

      if (true) { // debug
         std::string file_name = "trace-frag-forward-build-from-seed-";
         file_name += util::int_to_string(iseed);
         file_name += ".pdb";
         minimol::molecule mmm_out(fC);
         // add cell from map
         mmm_out.set_cell(acell);
         mmm_out.set_spacegroup(spacegroup_str_hm);
         mmm_out.write_file(file_name, 10);
      }

      if (true) {
	 // The variable names are confusing here, when building backwards
	 // the "next" residue from current (r, with res-no N) has res-no N+1.
	 // That is a residue that has been built (if we can find it).
	 //
	 // I need the full/partial residues in the seeds to be the other way round
	 // if I want to set "r_next" - which get's called r_prev
	 // (previous value of r) in backwards_2018().
	 //

	 // does this crash when r is null?
	 mmdb::Residue *r_next = 0;
	 minimol::fragment fN =
	    backwards_2018(iseed, r, r_next, "A", b_fact, n_trials, geom, xmap, mv, debugging);

	 if (true) { // debug
	    std::string file_name = "trace-frag-backwards-build-from-seed-";
	    file_name += util::int_to_string(iseed);
	    file_name += ".pdb";
	    minimol::molecule mmm_out(fN);
	    // add cell from map
	    mmm_out.set_cell(acell);
	    mmm_out.set_spacegroup(spacegroup_str_hm);
	    mmm_out.write_file(file_name, 10);
	 }

	 if (fN.n_filled_residues() > 3) { // needs tweaking
	    int build_dir = -1;
	    add_to_fragment_store(fN, build_dir);
	 }
      }
   }

   fragment_store.all_fragments_stored = true;

#ifdef HAVE_CXX_THREAD
   store_thread.join();
#endif // HAVE_CXX_THREAD

   std::cout << "The fragment store has " << fragment_store.stored_fragments.size()
	     << " fragments, of which " << fragment_store.n_sequenced() << " have sidechains"
	     << std::endl;

   minimol::molecule mmm;
   mmm.set_cell(acell);
   mmm.set_spacegroup(spacegroup_str_hm);
   for (std::size_t i=0; i<fragment_store.stored_fragments.size(); i++) {
      minimol::fragment frag = fragment_store.stored_fragments[i].frag;
      std::string chid = t.frag_idx_to_chain_id(i);
      frag.fragment_id = chid;
      mmm.fragments.push_back(frag);
   }

   mmm.write_file("built.pdb", 20.0);
   mmdb::Manager *mol = mmm.pcmmdbmanager();
   clipper::MMDBfile* mmdbfile = static_cast<clipper::MMDBfile*>(mol);
   clipper::MiniMol mm;
   mmdbfile->import_minimol(mm);
   bool r = ProteinTools::globularise(mm);
   mmdbfile->export_minimol(mm);
   mol->WritePDBASCII("built-2.pdb");

}


// with_sidechains is an optional arg, default false
void
coot::multi_build_terminal_residue_addition::add_to_fragment_store(const coot::minimol::fragment &new_fragment,
								   int build_dir,
								   bool with_sidechains) {

   stored_fragment_t f(new_fragment, build_dir, with_sidechains, standard_residues_mol);
   fragment_store.add(f, store_lock);
   mask_no_go_map(new_fragment);

}

void
coot::multi_build_terminal_residue_addition::mask_no_go_map(const coot::minimol::fragment &frag) {

   float atom_radius = 2.0;
   for (int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      const minimol::residue &r = frag[ires];
      if (r.n_atoms() > 0) {
	 for (unsigned int iat=0; iat<r.atoms.size(); iat++) {
	    const minimol::atom &at = r.atoms[iat];
	    clipper::Coord_frac cf = at.pos.coord_frac(no_go.cell());
	    clipper::Coord_frac box0(
				     cf.u() - atom_radius/no_go.cell().descr().a(),
				     cf.v() - atom_radius/no_go.cell().descr().b(),
				     cf.w() - atom_radius/no_go.cell().descr().c());

	    clipper::Coord_frac box1(
				     cf.u() + atom_radius/no_go.cell().descr().a(),
				     cf.v() + atom_radius/no_go.cell().descr().b(),
				     cf.w() + atom_radius/no_go.cell().descr().c());

	    clipper::Grid_map grid(box0.coord_grid(no_go.grid_sampling()),
				   box1.coord_grid(no_go.grid_sampling()));

	    float atom_radius_sq = atom_radius * atom_radius;
	    int nhit = 0;
	    int nmiss = 0;

	    clipper::Xmap_base::Map_reference_coord ix( no_go, grid.min() ), iu, iv, iw;
	    for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
	       for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
		  for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
		     if ( (iw.coord().coord_frac(no_go.grid_sampling()).coord_orth(no_go.cell()) - at.pos).lengthsq() < atom_radius_sq) {

			if (false)
			   std::cout << "masked point at "
				     << iw.coord().coord_frac(no_go.grid_sampling()).coord_orth(no_go.cell()).format()
				     << " centre point: " << at.pos.format() << " "
				     << (iw.coord().coord_frac(no_go.grid_sampling()).coord_orth(no_go.cell()) - at.pos).lengthsq()
				     << std::endl;
			no_go[iw] = 1;
			nhit++;
		     } else {
			nmiss++;
		     }
		  }
	       }
	    }
	 }
      }
   }
}

void
coot::stored_fragment_container_t::add(const stored_fragment_t &sf, std::atomic<unsigned int> &store_lock) {

   unsigned int unlocked = 0;
   while (! store_lock.compare_exchange_weak(unlocked, 1)) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
   }
   stored_fragments.push_back(sf);
   store_lock = 0;
}

// Do I need the lock? I guess that only this thread will try to alter frag, thre is no race.
//
// static
bool
coot::stored_fragment_t::try_assign_sidechains(coot::stored_fragment_t &stored_frag,
					       std::atomic<unsigned int> &locked,
					       const clipper::Xmap<float> &xmap,
					       const std::vector<std::pair<std::string, std::string> > &sequences,
					       mmdb::Manager *standard_residues_mol) {

   std::cout << "################## try_assign_sidechains() " << std::endl;

   if (stored_frag.sidechains_tried == false) {

      // ugly!
      //
      std::string pkg_data_dir = package_data_dir();
      std::string llkdfile = pkg_data_dir + "/cootaneer-llk-2.40.dat";
      // was that over-ridden?

      const char *cp = getenv("COOT_PREFIX");
      if (cp) {
	 llkdfile = cp;
	 llkdfile += "/share/coot/cootaneer-llk-2.40.dat";
      }
      if (!file_exists(llkdfile)) {
	 std::cout << "Ooops! Can't find cootaneer likelihood data! - failure"
		   << std::endl;
      } else {
	 Coot_sequence sequencer(llkdfile);
	 std::string chain_id = "A";
	 stored_frag.frag.fragment_id = "A"; // is this really what I want to do?
	 minimol::molecule mmm(stored_frag.frag);
	 mmdb::Manager *mol = mmm.pcmmdbmanager();
         if (! sequences.empty()) {
            std::cout << "---------- calling sequencer.sequence_chain "
                      << sequences.size() << " " << mol << " " << chain_id << std::endl;
            sequencer.sequence_chain(xmap, sequences, *mol, chain_id);
            std::cout << "---------- done sequencer.sequence_chain" << std::endl;
            std::string best_seq = sequencer.best_sequence();
            std::string full_seq = sequencer.full_sequence();
            double conf = sequencer.confidence();
            int chnnum = sequencer.chain_number();
            int chnoff = sequencer.chain_offset();
            std::cout << "Sequence: " << best_seq << "\nConfidence: " << conf << "\n";
            if (chnnum >= 0) {
               std::cout << "\nFrom    : " << full_seq << "\nChain id: "
                         << chnnum << "\tOffset: " << chnoff+1 << "\n";
               if (conf > 0.9) {
                  std::cout << "----------------------------- sequenced --------------------" << std::endl;
	       // modify stored_frag
                  apply_sequence(stored_frag, mol, best_seq, chnoff, standard_residues_mol, locked);
               }
            }
	 }
	 delete mol;
      }
      stored_frag.sidechains_tried = true;
   }
   return false;
}

// modify stored flag.  Use the store_lock for mutex access.
//
// static
void
coot::stored_fragment_t::apply_sequence(coot::stored_fragment_t &stored_frag,
					mmdb::Manager *mol, // frag as a mmdb::Manager *
					const std::string &best_seq,
					int resno_offset,
					mmdb::Manager *standard_residues_mol,
					std::atomic<unsigned int> &store_lock) {

   std::vector<mmdb::Residue *> residues;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) {
	    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    residues.push_back(residue_p);
	 }
      }
   }

   if (residues.size() != best_seq.length()) {
      std::cout << "ERROR:: oops residue mismatch " << best_seq.length() << " " << residues.size()
		<< std::endl;
   } else {
      // selected_chain_p->GetChain()->SetChainID(chain_id.c_str());
      std::vector<mmdb::Residue *> res_deletion; // after mutation, delete these
      for (unsigned int ichar=0; ichar<best_seq.length(); ichar++) {
	 const char &seq_char = best_seq[ichar];
	 if (seq_char == '?') {
	    std::cout << "DEBUG:: bypassing ? at " << ichar << std::endl;
	    // but we still need to set the sequence offset
	    mmdb::Residue *poly_ala_res = residues[ichar];
	    poly_ala_res->seqNum = resno_offset + ichar;
	 } else {
	    std::string res_type = coot::util:: single_letter_to_3_letter_code(best_seq[ichar]);
	    if (res_type != "") {
	       std::cout << "Mutating to " << res_type << " at " << ichar << std::endl;
	       mmdb::Residue *std_res = get_standard_residue_instance(res_type, standard_residues_mol);
	       if (std_res) {
		  mmdb::Residue *poly_ala_res = residues[ichar];
		  std::cout << "Mutating poly_ala residue number " << poly_ala_res->GetSeqNum() << std::endl;
		  std::string alt_conf = "";
		  util::mutate(poly_ala_res, std_res, alt_conf, 0); // not shelx
		  poly_ala_res->seqNum = resno_offset + ichar;
		  if (ichar < residues.size()) {
		     res_deletion.push_back(poly_ala_res);
		  }
	       }
	    }
	 }
      }

      if (residues.size() > 0) {
	 // now put the residues in residues into frag
	 unsigned int unlocked = 0;
	 while (! store_lock.compare_exchange_weak(unlocked, 1)) {
	    std::this_thread::sleep_for(std::chrono::microseconds(10));
	 }
	 for (std::size_t ires=0; ires<residues.size(); ires++) {
	    mmdb::Residue *r = residues[ires];
	    // std::cout << "  res-to-frag: " << residue_spec_t(r) << " " << r->GetResName() << std::endl;
	    minimol::residue mr(r);
	    stored_frag.frag[ires] = mr;
	 }
	 stored_frag.with_sidechains = true;
	 store_lock = 0;
      }

   }

   if (false)
      mol->WritePDBASCII("apply-sequence-results.pdb");
}

// static
mmdb::Residue *
coot::stored_fragment_t::get_standard_residue_instance(const std::string &residue_type_in,
						       mmdb::Manager *standard_residues_mol) {

   if (standard_residues_mol) {
      int imod = 1;
      mmdb::Model *model_p = standard_residues_mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       std::string res_name(residue_p->GetResName());
	       if (res_name == residue_type_in) {
		  mmdb::Residue *std_residue = util::deep_copy_this_residue(residue_p);
		  // std::cout << "get_standard_residue_instance returning " << std_residue << std::endl;
		  return std_residue;
	       }
	    }
	 }
      }
   }
   return 0;
}




// static
void
coot::multi_build_terminal_residue_addition::store_manager(stored_fragment_container_t &fragment_store,
							   std::atomic<unsigned int> &locked,
							   const clipper::Xmap<float> &xmap,
							   const std::vector<std::pair<std::string, std::string> > &sequences) {

   double score = 0; // nothing useful, but we wait for this value

#ifdef HAVE_CXX_THREAD

   bool all_done = false;

   unsigned int unlocked = 0;
   std::vector<std::thread> threads;

   while (! all_done) {

      while (! locked.compare_exchange_weak(unlocked, 1)) {
	 std::this_thread::sleep_for(std::chrono::microseconds(10));
      }

      // do something

      // first, the exit condition
      std::size_t n_done = 0;
      for (std::size_t i=0; i<fragment_store.stored_fragments.size(); i++) {
	 if (fragment_store.stored_fragments[i].sidechains_tried == false) {
	    if (fragment_store.stored_fragments[i].with_sidechains == true) {
	       n_done++;
	    }
	 } else {
	    n_done++;
	 }
      }
      if (false)
	 std::cout << "Here in store_manager: " << n_done << " vs "
		   << fragment_store.stored_fragments.size() << " "
		   << fragment_store.all_fragments_stored << std::endl;

      if (fragment_store.all_fragments_stored)
	 if (n_done == fragment_store.stored_fragments.size())
	    all_done = true;

      for (std::size_t i=0; i<fragment_store.stored_fragments.size(); i++) {
	 if (fragment_store.stored_fragments[i].sidechains_tried == false) {
	    if (fragment_store.stored_fragments[i].with_sidechains == false) {

	       threads.push_back(std::thread(stored_fragment_t::try_assign_sidechains,
					     std::ref(fragment_store.stored_fragments[i]),
					     std::ref(locked),
					     std::cref(xmap),
					     std::cref(sequences),
					     fragment_store.stored_fragments[i].standard_residues_mol));
	    }
	 }
      }

      locked = 0;
      std::this_thread::sleep_for(std::chrono::milliseconds(800)); // wait before checking the store again
   }
   unsigned int n_threads = threads.size();
   for (unsigned int i_thread=0; i_thread<n_threads; i_thread++)
      threads.at(i_thread).join();

#endif // HAVE_CXX_THREAD

}

// static
bool
coot::stored_fragment_t::matches_position(const position_triple_t &p1,
					  const position_triple_t &p2,
					  const std::vector<clipper::RTop_orth> &symms,
					  double d_crit_sqrd) {

   if (false) {
      std::cout << "Here are the symms: -------------- " << std::endl;
      for (std::size_t is=0; is<symms.size(); is++) {
	 std::cout << "     -------- " << is << std::endl;
	 std::cout << symms[is].format() << std::endl;
      }
   }

   unsigned int n_match = 0;
   for (std::size_t i=0; i<3; i++) {
      for (std::size_t is=0; is<symms.size(); is++) {
	 clipper::Coord_orth sp2(symms[is] * p2.positions[i]);
	 double dd = (p1.positions[i]-sp2).lengthsq();
	 if (dd < d_crit_sqrd) {
	    n_match++;
	    std::cout << "n_match " << n_match << " " << sqrt(dd) << " " << sqrt(d_crit_sqrd)
		      << " for isymm " << is << std::endl;
	    std::cout << symms[i].format() << std::endl;
	    break;
	 }
      }
   }

   if (n_match == 3)
      return true;

   return false;
}

bool
coot::stored_fragment_t::matches_position_in_fragment(const position_triple_t &res_triple_for_testing,
						      const std::vector<clipper::RTop_orth> &symms) const {

   bool match = false;
   double d_crit = 0.05;
   double dd = d_crit * d_crit;

   // std::cout << "in matches_position_in_fragment() here is res_triple_for_testing: " << std::endl;
   // for (std::size_t i=0; i<3; i++)
   // std::cout << "    " << res_triple_for_testing.positions[i].format() << std::endl;

   // std::cout << "in matches_position_in_fragment() there are " << residue_atom_positions.size()
   // << " residue atom triples " << std::endl;

   for (std::size_t i=0; i<residue_atom_positions.size(); i++) {
      if (false) {
	 std::cout << " residue_atom_positions " << i << std::endl;
	 std::cout << "c.f. " << res_triple_for_testing.positions[0].format() << " "
		   << residue_atom_positions[i].second.positions[0].format() << std::endl;
      }
      bool r = matches_position(res_triple_for_testing, residue_atom_positions[i].second, symms, dd);
      if (r) {
	 match = true;
	 break;
      }
   }

   return match;
}



void
coot::stored_fragment_t::fill_residue_atom_positions() {

   // fragment starts at 1 - hmm
   //    std::cout << "-------------------- this fragment range "
   // 	     << frag.min_res_no() << " " << frag.max_residue_number()
   // 	     << std::endl;

   for (int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
      const minimol::residue &r = frag[ires];
      if (false)
	 std::cout << "constructing a position_triple_t " << r << " with "
		   << " ires " << ires << " and "<< r.atoms.size() << " atoms " << std::endl;
      if (r.n_atoms() > 0) {
	 try {
	    position_triple_t pt(r);
	    std::pair<int, position_triple_t> rap(ires, pt);
	    residue_atom_positions.push_back(rap);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: " << rte.what() << std::endl;
	 }
      }
   }
}

void
coot::stored_fragment_t::position_triple_t::fill_residue_atom_positions(const minimol::residue &r) {

   unsigned int n_found = 0;

   // std::cout << "in fill_residue_atom_positions() " << r.n_atoms() << std::endl;
   if (r.atoms.size() > 0) {
      std::vector<clipper::Coord_orth> p(3);
      for (unsigned int iat=0; iat<r.atoms.size(); iat++) {
	 const minimol::atom &at = r.atoms[iat];
	 // std::cout << "setting position from atom " << at << std::endl;
	 if (at.name == " N  ") {
	    p[0] = at.pos;
	    n_found++;
	 }
	 if (at.name == " CA ") {
	    p[1] = at.pos;
	    n_found++;
	 }
	 if (at.name == " C  ") {
	    p[2] = at.pos;
	    n_found++;
	 }
      }
      if (n_found == 3) {
	 for (std::size_t i=0; i<3; i++)
	    positions[i] = p[i];
      }
   }


   if (n_found != 3) {
      std::cout << "ERROR in fill_residue_atom_positions() n_found " << n_found << " in residue with "
		<< r.atoms.size() << " atoms " << std::endl;
      std::string m("in fill_residue_atom_positions(): missing atoms: ");
      m += util::int_to_string(n_found);
      throw(std::runtime_error(m));
   }

   // std::cout << "done fill_residue_atom_positions() " << n_found << std::endl;

}


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
coot::multi_build_terminal_residue_addition::forwards_2018(unsigned int iseed,
							   mmdb::Residue *res_p,
							   mmdb::Residue *upstream_neighbour_p,
							   const std::string &chain_id,
							   float b_factor_in,
							   int n_trials,
							   const coot::protein_geometry &geom,
							   const clipper::Xmap<float> &xmap,
							   std::pair<float, float> mv,
							   bool debug_trials) {


   unsigned int n_threads = coot::get_max_number_of_threads();
   if (n_threads < 1) n_threads = 1;
   ctpl::thread_pool thread_pool(n_threads);

   minimol::fragment many_residues; // add residue to this, and return it
   many_residues.addresidue(minimol::residue(res_p), false);

   // setup up initial conditions for the while loop
   bool prev_happy_fit = true;
   bool happy_fit = true; // starting
   int icount = 0;
   mmdb::Residue *res_prev_p = 0; // whatever was res_p the prevous round
   minimol::residue res_prev; // what was res on the previous round.

   while (happy_fit) {

      if (true) { // debug
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

      addres.thread_pool(&thread_pool, n_threads);

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
      update_O_position_in_prev_residue(res_p, &many_residues, res);

      std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(res);
      if (cbeta_info.first) {
	 res.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 20.0);
	 many_residues[res.seqnum].addatom(" CB ", " C", cbeta_info.second, "", 1.0, 30.0);
      }

      // res is a copy of the residue in many_residues

      // refine res and its neighbour and one more, update the atoms of many_residues
      //
      int offset = 1;
      refine_end(&many_residues, res.seqnum, offset, geom, xmap, &thread_pool, n_threads);

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

	    int build_dir = 1; // forwards
	    bool ab_flag = was_this_already_built_p(res, iseed, build_dir, store_lock);
	    if (ab_flag) {
	       happy_fit = false; // stop building
	       std::cout << "Forward-building with ab_flag: " << ab_flag << " for seed " << iseed << std::endl;
	    }

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
   }

   return many_residues;
}

coot::minimol::fragment
coot::multi_build_terminal_residue_addition::backwards_2018(unsigned int iseed,
							    mmdb::Residue *res_p,
							    mmdb::Residue *upstream_neighbour_p,
							    const std::string &chain_id,
							    float b_factor_in,
							    int n_trials,
							    const coot::protein_geometry &geom,
							    const clipper::Xmap<float> &xmap,
							    std::pair<float, float> mv,
							    bool debug_trials) {

   unsigned int n_threads = get_max_number_of_threads();
   if (n_threads < 1) n_threads = 1;
   ctpl::thread_pool thread_pool(n_threads);

   minimol::fragment many_residues;
   many_residues.addresidue(minimol::residue(res_p), false);

   // setup initial conditions for the while loop
   bool happy_fit = true;
   int icount = 0;
   mmdb::Residue *res_prev_p = 0; // this will become whatever res_p was the previous round
   minimol::residue res_prev; // whatever res was on the previous round

   while (happy_fit) {

      if (false) { // debugging
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 res_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    std::cout << "Atom " << iat << " " << coot::atom_spec_t(residue_atoms[iat]) << std::endl;
	 }
      }

      // fill res, the first of the 2 residues fitted by residue_by_phi_psi
      minimol::residue res;
      std::string residue_type = "ALA"; // for the first try. If it's ALA we will add CB.
      residue_by_phi_psi addres("MN", res_p, chain_id, residue_type, 20);

      addres.thread_pool(&thread_pool, n_threads);

      addres.set_downstream_neighbour(res_prev_p);
      addres.import_map_from(xmap);
      minimol::fragment f = addres.best_fit_phi_psi(n_trials, -1); // offset = -1

      int ires = res_p->GetSeqNum() - 1;

      if (false) { // debug
	 std::string file_name = "N-terminal-best-fit-phi-psi-";
	 file_name += util::int_to_string(ires) + std::string(".pdb");
	 minimol::molecule mmm(f);
	 mmm.write_file(file_name, 10);
      }

      // add frag to many_residues, and set res
      bool found_res = false;
      res = f[ires];

      many_residues.addresidue(f[ires],   false);
      many_residues.addresidue(f[ires-1], false); // the "next" residue, Needed so that we
                                                  // can well refine res.
      std::cout << "=== res.seqnum " << res.seqnum << " ===== ires " << ires << std::endl;
      std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(res);
      if (cbeta_info.first) {
	 res.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 20.0);
	 many_residues[res.seqnum].addatom(" CB ", " C", cbeta_info.second, "", 1.0, 30.0);
      }
      int offset = -1;
      refine_end(&many_residues, res.seqnum, offset, geom, xmap, &thread_pool, n_threads);
      // now extract res from many residues now that it has been updated.
      res = many_residues[res.seqnum];
      happy_fit = does_residue_fit(res, xmap, mv);

      if (false) { // debug
	 std::string file_name = "N-terminal-build-post-refined-";
	 file_name += util::int_to_string(res.seqnum) + std::string(".pdb");
	 minimol::molecule mmm(many_residues);
	 mmm.write_file(file_name, 10);
      }

      many_residues.delete_first_residue(); // remove the "other" residue

      if (res.is_empty())
	 happy_fit = false;

      if (happy_fit) {

	 bool ab_flag = was_this_already_built_p(res, iseed, -1, store_lock);
	 // std::cout << "Here with ab_flag: " << ab_flag << std::endl;
	 if (ab_flag) {
	    happy_fit = false; // stop building this fragment
	    std::cout << "Backwards building with ab_flag: " << ab_flag << " for seed " << iseed << std::endl;
	 }
      }

      if (happy_fit) {

	 icount++;

	 // setup for next round

	 // set res_p to be the newly constructed residue
	 res_prev = res;
	 res_prev_p = res_p; // save to set downstream_N
	 res_p = res.make_residue();
	 if (false) { // debug
	    mmdb::Atom **residue_atoms = 0;
	    int n_residue_atoms;
	    res_p->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       std::cout << "Atom of res_p from make_residue() " << iat << " "
			 << atom_spec_t(residue_atoms[iat]) << std::endl;
	    }
	 }
      }
   }
   return many_residues;
}

// for now make this in the main thread.  We should pass the thread pool
// if we want to speed it up.
//
bool
coot::multi_build_terminal_residue_addition::was_this_already_built_p(coot::minimol::residue &res,
								      unsigned int seed_number,
								      int build_dir,
								      std::atomic<unsigned int> &store_lock) const {

   bool matches_previous = false;
   unsigned int unlocked = 0;
   stored_fragment_t::position_triple_t res_triple_pos(res);

   if (is_in_no_go_map(res)) {
      matches_previous = true;
   } else {

      if (false) {
	 std::cout << "in was_this_already_built_p, this is: " << std::endl;
	 for (std::size_t i=0; i<3; i++)
	    std::cout << "    " << res_triple_pos.positions[i].format() << std::endl;
      }

      while (! store_lock.compare_exchange_weak(unlocked, 1)) {
	 std::this_thread::sleep_for(std::chrono::microseconds(10));
      }

      for (std::size_t i=0; i<fragment_store.stored_fragments.size(); i++) {
	 if (fragment_store.stored_fragments[i].build_dir == build_dir) {
	    if (fragment_store.stored_fragments[i].matches_position_in_fragment(res_triple_pos, symms)) {
	       matches_previous = true;
	       std::cout << "|||||||||||||| seed number " << seed_number << " build-dir " << build_dir
			 << " matched by stored fragment number " << i << std::endl;
	       break;
	    }
	 }
      }
      store_lock = 0;
   }
   return matches_previous;
}


void
coot::multi_build_terminal_residue_addition::update_O_position_in_prev_residue(mmdb::Residue *res_p,
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
coot::multi_build_terminal_residue_addition::try_to_recover_from_bad_fit_forwards(mmdb::Residue *res_p, mmdb::Residue *res_prev_p,
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

   // Without adding the thread pool, we go through a different code path in fit_terminal_residue_generic()
   // and that can be useful for debugging.
   //
   unsigned int n_threads_max = coot::get_max_number_of_threads();
   if (n_threads_max < 1) n_threads_max = 1;

   ctpl::thread_pool thread_pool(n_threads_max); // pass the tread pool
   rpp.thread_pool(&thread_pool, n_threads_max);

   minimol::fragment f = rpp.best_fit_phi_psi(n_trials * 8, 1);

   // we don't refine here - is that a good idea?
   //
   // No, it's not. Add refinement.

   int new_res_seqnum = res_p->GetSeqNum()+1;
   refine_end(&f, new_res_seqnum, 1, geom, xmap, &thread_pool, n_threads_max);

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

   unsigned int n_threads = coot::get_max_number_of_threads();
   if (n_threads < 1) n_threads = 1;
   ctpl::thread_pool thread_pool(n_threads);

   // kludge, for backwards compatibility. Maybe I should just delete this function
   //
   std::vector<std::pair<std::string, std::string> > seqs;
   multi_build_terminal_residue_addition mbtra(geom, xmap, mv, seqs);

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
	    
	 bool this_happy_fit = mbtra.does_residue_fit(res, xmap, mv);

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
	    mbtra.refine_end(&many_residues, res.seqnum, offset, geom, xmap, &thread_pool, n_threads);

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
coot::multi_build_terminal_residue_addition::refine_end(coot::minimol::fragment *many_residues,
							int seqnum, int offset,
							const coot::protein_geometry &geom,
							const clipper::Xmap<float> &xmap,
							ctpl::thread_pool *thread_pool_p, int n_threads) {

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
   restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, &xmap);
   restraints.thread_pool(thread_pool_p, n_threads);

   // Does this make things slower? (seems so, try passing the thread pool
   // or test how long it takes to create and add).
   //
   ctpl::thread_pool thread_pool(coot::get_max_number_of_threads());
   restraints.thread_pool(&thread_pool, get_max_number_of_threads());

   pseudo_restraint_bond_type pseudos = NO_PSEUDO_BONDS;
   bool do_internal_torsions = false;
   float weight = 60.0f;

   restraints.set_quiet_reporting();
   restraints.add_map(weight);
   bool do_trans_peptide_restraints = true;
   do_trans_peptide_restraints = false;
   int imol = 0;
   restraints.make_restraints(imol, geom, flags, do_internal_torsions, do_internal_torsions, 0, 0, true, true, false, pseudos);
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
coot::multi_build_terminal_residue_addition::does_residue_fit(const coot::minimol::residue &res,
							      const clipper::Xmap<float> &xmap,
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
coot::multi_build_terminal_residue_addition::crashing_into_self(const minimol::fragment &many_residues, int seqnum, int offset) {

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


// return true if the atoms of the residue are in the no-go map (i.e. it
// had been masked because it has been built before).
//
bool
coot::multi_build_terminal_residue_addition::is_in_no_go_map(minimol::residue &r) const {

   bool status = false;
   float atom_radius = 1.0;
   bool big_enough_atom_radius = false; // initial value to start the while

   int n_hit = 0;
   int n_miss = 0;
   while (! big_enough_atom_radius) {
      if (r.n_atoms() > 0) {
	 n_hit = 0;
	 n_miss = 0;
	 float atom_radius_sq = atom_radius * atom_radius;
	 for (unsigned int iat=0; iat<r.atoms.size(); iat++) {
	    const minimol::atom &at = r.atoms[iat];
	    clipper::Coord_frac cf = at.pos.coord_frac(no_go.cell());
	    clipper::Coord_frac box0(
				     cf.u() - atom_radius/no_go.cell().descr().a(),
				     cf.v() - atom_radius/no_go.cell().descr().b(),
				     cf.w() - atom_radius/no_go.cell().descr().c());

	    clipper::Coord_frac box1(
				     cf.u() + atom_radius/no_go.cell().descr().a(),
				     cf.v() + atom_radius/no_go.cell().descr().b(),
				     cf.w() + atom_radius/no_go.cell().descr().c());

	    clipper::Grid_map grid(box0.coord_grid(no_go.grid_sampling()),
				   box1.coord_grid(no_go.grid_sampling()));

	    clipper::Xmap_base::Map_reference_coord ix( no_go, grid.min() ), iu, iv, iw;
	    for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
	       for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
		  for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
		     if ( (iw.coord().coord_frac(no_go.grid_sampling()).coord_orth(no_go.cell()) - at.pos).lengthsq() < atom_radius_sq) {
			if (no_go[iw])
			   ++n_hit;
			else
			   ++n_miss;
		     }
		  }
	       }
	    }
	 }
      }

      // std::cout << "debug:: nhit " << n_hit << " nmiss: " << n_miss << std::endl;

      if ((n_hit + n_miss)  > 12) {
	 big_enough_atom_radius = true;
	 if (n_hit > 6) {
	    status = true;
	 }
      } else {
	 atom_radius *= 1.5;
      }
   }
   return status;
}
