/* ligand/trace.cc
 * 
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

#include <iomanip>
#include <algorithm>
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "trace.hh"
#include "ligand.hh"
#include "analysis/stats.hh"
#include "analysis/kolmogorov.hh" // for testing
#include "multi-peptide.hh"
#include "mini-mol/mini-mol-utils.hh"

// First we will flood the map, then look for atom pair that are 2.8
// to 4.8A apart and spin a search the density

coot::trace::trace(const clipper::Xmap<float>& xmap_in) {

   xmap = xmap_in;
   flood_atom_mask_radius = 1.5;
   flood_atom_mask_radius = 1.1; 
   rmsd_cut_off = 1.5;

   mol = 0;
   atom_selection = 0;
   n_selected_atoms = 0;

   using_test_model = false; 

   frag_score_crit = 0;

   std::pair<float, float> mv = util::mean_and_variance(xmap);

   map_mean = mv.first;
   map_rmsd = sqrt(mv.second);

   // set refinement weight to 20/rmsd

   set_inital_density_scales();
   
}

std::ostream&
coot::operator<<(std::ostream &s, scored_node_t sn) {

   s << "[" << sn.atom_idx << " " << sn.spin_score << " "
     << sn.alpha << " " << sn.reversed_flag << "]";
   return s;
}

std::vector<coot::minimol::fragment>
coot::trace::make_seeds() {

   // seeds are a 2-residue fragment with the "interesting" residue as 2 and
   // with a few atoms (CA,C,O) of the paired (upstream) residue.

   // if not cryo-EM, then move seeds to around origin before returning.

   // note to future self: this function does not need a backwards version.

   std::vector<coot::minimol::fragment> seeds;
   bool debug = true;
   minimol::molecule flood_mol = get_flood_molecule();

   mmdb::Manager *action_mol = flood_mol.pcmmdbmanager();

   if (action_mol) {

      double variation = 0.25; // make bigger at lower resolutions (maybe up to 0.5?)
      std::vector<std::pair<unsigned int, unsigned int> > apwd =
	 atom_pairs_within_distance(action_mol, 3.81, variation);
      std::vector<std::pair<unsigned int, scored_node_t> > scores = spin_score_pairs(apwd);

      unsigned int n_top_spin_pairs = 1000;

      std::sort(scores.begin(), scores.end(), scored_node_t::sort_pair_scores);
      if (scores.size() < n_top_spin_pairs)
	 n_top_spin_pairs = scores.size();
      else
	 scores.resize(n_top_spin_pairs); // this is the right way

      if (false) { // debugging

	 for (std::size_t i=0; i<scores.size(); i++)
	    std::cout << "   scores: " << i << " " << scores[i].first << " "
		      << scores[i].second.spin_score << std::endl;

	 for (std::size_t i=0; i<scores.size(); i++) {
	    const std::pair<unsigned int, scored_node_t> &node = scores[i];
	    minimol::fragment frag = make_fragment(node, 1, "A");
	    if (true) { // debug
	       std::string file_name = "frag-seeds-";
	       file_name += coot::util::int_to_string(i) + std::string(".pdb");
	       minimol::molecule mmm(frag);
	       mmm.write_file(file_name, 10);
	    }
	 }
      }

      double good_enough_score = 0.5; // overwritten

      if (scores.size() > 0)
	 good_enough_score = scores.back().second.spin_score;

      make_connection_map(scores); // fills fwd_connection_map

      dir_t dir = FORWARDS;
      std::string chain_id("A");
      for (std::size_t i=0; i<scores.size(); i++) {
	 const std::pair<unsigned int, scored_node_t> &node = scores[i];
	 if (node.second.spin_score > good_enough_score) {

	    std::map<unsigned int, std::vector<scored_node_t> >::const_iterator it;
	    it = fwd_connection_map.find(node.second.atom_idx);
	    if (it != fwd_connection_map.end()) {

	       const std::vector<scored_node_t> &v = it->second;
	       if (v.size() > 0) {
		  // maybe check _all_ connections?
		  unsigned int idx_next = v[0].atom_idx;
		  std::pair<unsigned int, scored_node_t> next_node(node.second.atom_idx, v[0]);
		  if (next_node.second.spin_score > good_enough_score) {
		     minimol::fragment frag = make_residue(node, next_node, 98, "A");
		     seeds.push_back(frag);
		  }
	       }
	    }
	 }
      }
   }

   if (seeds.size() > 0) {
      if (! util::is_EM_map(xmap)) {
	 if (true) {
	    for (std::size_t i=0; i<seeds.size(); i++) {
	       // std::string fn = "seed-residue-pre-moved-" + util::int_to_string(i) + ".pdb";
	       // frag_to_pdb(seeds[i], fn); // adds cell and symmetry for output
	    }
	 }
	 move_seeds_close_to_origin(&seeds); // by translation
	 if (true) {
	    for (std::size_t i=0; i<seeds.size(); i++) {
	       // std::string fn = "seed-residue-" + util::int_to_string(i) + ".pdb";
	       // frag_to_pdb(seeds[i], fn); // adds cell and symmetry for output
	    }
	 }
      }
   }
   return seeds;
}

void
coot::trace::move_seeds_close_to_origin(std::vector<coot::minimol::fragment> *seeds_in) const {

   for (std::size_t i=0; i<seeds_in->size(); i++) {
      minimol::fragment &seed = seeds_in->at(i);

      std::vector<clipper::Coord_orth> coords;
      unsigned int n = 0;
      for (int ires=seed.min_res_no(); ires<=seed.max_residue_number(); ires++) {
	 for (unsigned int iatom=0; iatom<seed[ires].atoms.size(); iatom++) {
	    const clipper::Coord_orth &p = seed[ires][iatom].pos;
	    coords.push_back(p);
	    n++;
	 }
      }
      if (n > 0) {
	 const clipper::Spacegroup &spg = xmap.spacegroup();
	 const clipper::Cell &cell = xmap.cell();
	 clipper::Coord_frac cf = util::shift_to_origin(coords, xmap.cell(), xmap.spacegroup());
	 if (! util::is_000_shift(cf)) {
	    // shift 'em
	    clipper::Coord_orth co = cf.coord_orth(xmap.cell());
	    for (int ires=seed.min_res_no(); ires<=seed.max_residue_number(); ires++) {
	       for (unsigned int iatom=0; iatom<seed[ires].atoms.size(); iatom++) {
		  if (false) { // debug
		     double d1 = sqrt(seed[ires][iatom].pos.lengthsq());
		     double d2 = sqrt((seed[ires][iatom].pos+co).lengthsq());
		     std::cout << "seed origin-moved " << ires << " " << iatom
			       << " from " << seed[ires][iatom].pos.format() << " to "
			       << (seed[ires][iatom].pos + co).format()
			       << " by " << co.format() << " d1 " << d1 << " d2 " << d2
			       << std::endl;
		  }
		  seed[ires][iatom].pos += co;
	       }
	    }
	 }

	 clipper::Coord_orth sum(0,0,0);
	 for (int ires=seed.min_res_no(); ires<=seed.max_residue_number(); ires++) {
	    for (unsigned int iatom=0; iatom<seed[ires].atoms.size(); iatom++) {
	       const clipper::Coord_orth &p = seed[ires][iatom].pos;
	       sum += p;
	    }
	 }
	 double oon = 1.0/static_cast<double>(n);
	 clipper::Coord_orth frag_av_pos(sum.x() * oon, sum.y() * oon, sum.z() * oon);
	 // OK, can I rotate them with symmetry so that they are closer?
	 double d2min = frag_av_pos.lengthsq() -0.01; // small delta is useful
	 int ibest = -1;
	 clipper::RTop_orth rtop_best;
	 unsigned int ns = spg.num_symops();
	 for (int x_shift = -1; x_shift<2; x_shift++) {
	    for (int y_shift = -1; y_shift<2; y_shift++) {
	       for (int z_shift = -1; z_shift<2; z_shift++) {
		  clipper::Coord_frac cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift);
		  for (unsigned int is = 0; is < ns; is++) {
		     const clipper::Symop &symop = spg.symop(is);
		     clipper::RTop_frac rtf(symop.rot(), symop.trn() + cell_shift);
		     clipper::RTop_orth rtop = symop.rtop_orth(cell);
		     clipper::Coord_orth np = rtop * frag_av_pos;
		     double d2 = np.lengthsq();
		     if (d2 < d2min) {
			d2min = d2;
			ibest = is;
			rtop_best = rtop;
		     }
		  }
	       }
	    }
	 }
	 // rotate 'em
	 if (ibest >= 0) {
	    // number 0 would be surprising
	    if (false) // debug
	       std::cout << "DEBUG:: rotating fragment seed  " << i << " using symop idx "
			 << ibest << std::endl;
	    for (int ires=seed.min_res_no(); ires<=seed.max_residue_number(); ires++) {
	       for (unsigned int iatom=0; iatom<seed[ires].atoms.size(); iatom++) {
		  if (false) { // debug
		     double d1 = sqrt(seed[ires][iatom].pos.lengthsq());
		     double d2 = sqrt((rtop_best * seed[ires][iatom].pos).lengthsq());
		     std::cout << "seed symm-application " << ires << " " << iatom
			       << " from " << seed[ires][iatom].pos.format() << " to "
			       << (rtop_best * seed[ires][iatom].pos).format()
			       << " d1 " << d1 << " d2 " << d2
			       << std::endl;
		  }
		  seed[ires][iatom].pos = rtop_best * seed[ires][iatom].pos;
	       }
	    }
	 }
      }
   }
}

void
coot::trace::action() {

   bool do_trace_by_spin_pairs = true;
   minimol::molecule flood_mol = get_flood_molecule();

   mmdb::Manager *action_mol = flood_mol.pcmmdbmanager();

   if (action_mol) {
      std::vector<std::pair<unsigned int, unsigned int> > apwd =
	 atom_pairs_within_distance(action_mol, 3.81, 1.0);
      std::vector<std::pair<unsigned int, scored_node_t> > scores = spin_score_pairs(apwd);

      // Now start from the best peptide and try to put coordinate there

      unsigned int n_top_spin_pairs_for_tracing_starts = 200;
      unsigned int n_top_spin_pairs = 300; // use only the top 1000 for
                                            // potential connections

      std::sort(scores.begin(), scores.end(), scored_node_t::sort_pair_scores);
      if (scores.size() < n_top_spin_pairs)
	 n_top_spin_pairs = scores.size();
      else
	 scores.resize(n_top_spin_pairs);

      make_connection_map(scores);
      set_frag_score_crit(scores);
      if (using_test_model)
	 ks_test(scores);

      std::vector<std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> >
	 frag_store;

      if (do_trace_by_spin_pairs) {  // fill the frag_store
	 for (unsigned int i=0; i<n_top_spin_pairs_for_tracing_starts; i++) {
	    std::string chain_id = frag_idx_to_chain_id(i);

	    std::vector<scored_node_t> start_path;
	    scored_node_t dummy_first;
	    dummy_first.atom_idx = scores[i].first;
	    start_path.push_back(dummy_first);
      
	    std::cout << "----------- test_model() starting point number  " << i
		      << " atom_idx: " << scores[i].first << " node: " << scores[i].second 
		      << " chain_id " << chain_id << std::endl;

	    int res_no_base = scores[i].first;

	    res_no_base = i;

	    // actually, we *have* the first spin pair but it gets thrown
	    // away and "rediscovered" when build_2_choose_1() does just
	    // that (and hopefully picks the idx scores[i].second.atom_ids).
	    // 

	    dir_t dir = FORWARDS;
      
	    std::cout << "------------------------------- follow forwards ---------"
		      << std::endl;
	    std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> ff = 
	       follow_fragment(scores[i].first, start_path, res_no_base, chain_id, dir);

	    std::cout << "------------------------------- follow backwards --------"
		      << std::endl;
	    dir = BACKWARDS;
	    std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> bf = 
	       follow_fragment(scores[1].second.atom_idx, start_path, 1000+res_no_base, chain_id, dir);
	 
	    add_replace_reject(frag_store, ff);
	    add_replace_reject(frag_store, bf);

	    std::cout << "For spin-pair " << i << " we have " << frag_store.size()
		      << " fragment " << std::endl;
	 }

	 if (false) { 
	    for (unsigned int ifrag=0; ifrag<frag_store.size(); ifrag++) {
	       std::string  node_string = util::int_to_string(ifrag);
	       std::string fn_ = "frag-store-" + node_string + ".pdb";
	       frag_to_pdb(frag_store[ifrag].second, fn_);
	    }
	 }
      }


      // Rama terminal addition/refine trace.
      //
      std::pair<float, float> mv = util::mean_and_variance(xmap);
      protein_geometry geom;
      geom.init_standard();
      geom.remove_planar_peptide_restraint();
      multi_peptide(frag_store, geom, mv);
      
   }
}

void
coot::trace::frag_to_pdb(const minimol::fragment &frag, const std::string &fn) const {

   minimol::molecule m(frag);
   if (! m.is_empty()) {
      m.set_cell(xmap.cell());
      m.set_spacegroup(xmap.spacegroup().symbol_hm());
      m.write_file(fn, 10);
   }
}

std::string
coot::trace::frag_idx_to_chain_id(unsigned int idx) const {

   std::string s = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";

   if (idx < s.length()) {
      char c = s[idx];
      std::string ss(1,c);
      return ss;
   } else {
      return "Z";
   }
}

double
coot::trace::ks_test(const std::vector<std::pair<unsigned int, scored_node_t> > &scores) {

   double ks = 0;

   std::vector<double> t1;
   std::vector<double> t2;
   for (unsigned int i=0; i<scores.size(); i++) {
      if (scores[i].second.udd_flag)
	 t2.push_back(scores[i].second.spin_score);
      else 
	 t1.push_back(scores[i].second.spin_score);
   }
   ks = nicholls::get_KS(t1, t2);
   std::cout << "debug:: ks: t1.size() " << t1.size() << " t2.size() " << t2.size()
	     << " ks: " << ks << std::endl;
   return ks;
}

double
coot::trace::ks_test() {

   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atom_pairs_within_distance(mol, 3.81, 1.0);
   std::vector<std::pair<unsigned int, scored_node_t> > scores = spin_score_pairs(apwd);
   return ks_test(scores);

} 



coot::minimol::molecule
coot::trace::get_flood_molecule() const {

   bool debug = true;
   coot::ligand lig;

   lig.set_cluster_size_check_off();
   lig.set_chemically_sensible_check_off();
   lig.set_sphericity_test_off();
	       
   lig.set_map_atom_mask_radius(flood_atom_mask_radius);
   lig.set_water_to_protein_distance_limits(10.0, 1.5);
   
   lig.import_map_from(xmap);
   
   lig.flood2(rmsd_cut_off);
   coot::minimol::molecule water_mol = lig.water_mol();

   if (debug) {
      std::string output_pdb = "flood-mol.pdb";
      water_mol.write_file(output_pdb, 30.0);
      lig.output_map("find-waters-masked-flooded.map");
   }
   
   return water_mol;
   
}

std::vector<std::pair<unsigned int, unsigned int> >
coot::trace::atom_pairs_within_distance(mmdb::Manager *mol_in,
					double trans_dist,
					double trans_dist_variation) {

   // set class members mol, atom_selection and n_sel_atoms.
   mol = mol_in;  // the peaks in the map - some of which are CAs hopefully.
   

   std::vector<std::pair<unsigned int, unsigned int> > v;
   if (mol) {
      selhnd = mol->NewSelection(); // d (in destructor)
      mol->SelectAtoms(selhnd, 0, "*",
		       mmdb::ANY_RES, // starting resno, an int
		       "*", // any insertion code
		       mmdb::ANY_RES, // ending resno
		       "*", // ending insertion code
		       "*", // any residue name
		       "*", // atom name
		       "*", // elements
		       "");

      atom_selection = 0; // member data - cleared on destruction
      n_selected_atoms = 0;
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

      std::cout << "selected " << n_selected_atoms << " for distance pair check"
		<< std::endl;

      int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM, "index");
      if (uddHnd<0)  {
	 std::cout << " atom bonding registration failed.\n";
      } else {
	 for (int i=0; i< n_selected_atoms; i++) {
	    mmdb::Atom *at = atom_selection[i];
	    at->PutUDData(uddHnd, i); // is this needed any more?
	 }
 
	 mmdb::Contact *pscontact = NULL;
	 int n_contacts;
	 long i_contact_group = 1;
	 mmdb::mat44 my_matt;
	 for (int i=0; i<4; i++) 
	    for (int j=0; j<4; j++) 
	       my_matt[i][j] = 0.0;      
	 for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
	 //
	 mmdb::realtype local_dist_max = trans_dist + trans_dist_variation;
	 mmdb::realtype local_dist_min = trans_dist - trans_dist_variation;

	 std::cout << "debug:: SeekContacts with distance limits "
		   << local_dist_min << " " << local_dist_max << std::endl;
	 
	 mol->SeekContacts(atom_selection, n_selected_atoms,
			   atom_selection, n_selected_atoms,
			   local_dist_min, local_dist_max,
			   0,        // seqDist 0 -> in same res also
			   pscontact, n_contacts,
			   0, &my_matt, i_contact_group);

	 if (n_contacts > 0) {
	    if (pscontact) {
	       for (int i=0; i<n_contacts; i++) {
		  if (pscontact[i].id1 < pscontact[i].id2) { 
		     mmdb::Atom *at_1 = atom_selection[pscontact[i].id1];
		     mmdb::Atom *at_2 = atom_selection[pscontact[i].id2];
		     int idx_1, idx_2;
		     at_1->GetUDData(uddHnd, idx_1);
		     at_2->GetUDData(uddHnd, idx_2);
		     std::pair<unsigned int, unsigned int> p(idx_1, idx_2);

		     if (1) {  // debug
			// check that saa atoms are at the same place
			// as sel_atoms atoms.
			clipper::Coord_orth sel_atom_1_co = co(at_1);
			clipper::Coord_orth sel_atom_2_co = co(at_2);
		     } 
			
		     v.push_back(p);
		  }
	       }
	    }
	 }

	 std::cout << "found " << n_contacts << " distance pairs " << std::endl;
	 std::cout << "made  " << v.size() << " distance pairs " << std::endl;
      }
   }
   return v;

}
   

// return sorted scores
// 
std::vector<std::pair<unsigned int, coot::scored_node_t> > 
coot::trace::spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd) {

   // apwd : atom (index) pairs within distance

   unsigned int n_top = 5000; // top-scoring spin-score pairs for tracing
   bool debug = false;
   
   std::vector<std::pair<unsigned int, scored_node_t> > scores(apwd.size()*2);

#pragma omp parallel for shared(scores)   
   for (unsigned int i=0; i<apwd.size(); i++) {
      const unsigned &at_idx_1 = apwd[i].first;
      const unsigned &at_idx_2 = apwd[i].second;
      scores[2*i ]  = spin_score(at_idx_1, at_idx_2);
      if (debug)
	 output_spin_score(scores[2*i  ], at_idx_1, at_idx_2);
   }
   for (unsigned int i=0; i<apwd.size(); i++) {
      const unsigned &at_idx_1 = apwd[i].first;
      const unsigned &at_idx_2 = apwd[i].second;
      scores[2*i+1] = spin_score(at_idx_2, at_idx_1);
      if (debug)
	 output_spin_score(scores[2*i+1], at_idx_2, at_idx_1);
   }

   if (scores.size() > n_top) {
      std::sort(scores.begin(), scores.end(), scored_node_t::sort_pair_scores);
      scores.resize(n_top);
   } else {
      n_top = scores.size();
   }

   if (debug) { 
      std::cout << "---- sorted scores ----- " << n_top <<  std::endl;
      for (unsigned int i=0; i<n_top; i++) {
	 const std::string &at_name_1 = index_to_name(scores[i].first);
	 const std::string &at_name_2 = index_to_name(scores[i].second.atom_idx);
	 int res_no_1 = atom_selection[scores[i].first          ]->GetSeqNum();
	 int res_no_2 = atom_selection[scores[i].second.atom_idx]->GetSeqNum();
	 std::string chain_id_1 = atom_selection[scores[i].first          ]->GetChainID();
	 std::string chain_id_2 = atom_selection[scores[i].second.atom_idx]->GetChainID();
      
	 std::cout << "sorted spin scores " << " " << std::setw(4) << scores[i].first << " "
		   << " to " << std::setw(4) << scores[i].second.atom_idx << " "
		   << at_name_1 << " " << std::setw(3) << res_no_1 << " " << chain_id_1 << " " 
		   << at_name_2 << " " << std::setw(3) << res_no_2 << " " << chain_id_2 << " " 
		   << index_to_pos(scores[i].first).format() << " "
		   << index_to_pos(scores[i].second.atom_idx).format() << " "
		   << scores[i].second.spin_score << std::endl;
      }
   }
   return scores;
   
}

void
coot::trace::output_spin_score(const std::pair<unsigned int, scored_node_t> &score,
			       unsigned int at_idx_1,
			       unsigned int at_idx_2) const {

   // debugging output
   //
   //
   bool ca_1_flag = false;
   bool ca_2_flag = false;
   bool consecutive_flag = false;
   // 
   // indexing looks OK.
   if (index_to_name(at_idx_1) == " CA ") ca_1_flag = true;
   if (index_to_name(at_idx_2) == " CA ") ca_2_flag = true;
   // 
   if (ca_1_flag && ca_2_flag) {
      int resno_delta =
	 atom_selection[at_idx_2]->GetSeqNum() -
	 atom_selection[at_idx_1]->GetSeqNum();
      if (resno_delta == 1)
	 consecutive_flag = true;
   } 
   // 
   std::string at_name_1 = index_to_name(at_idx_1);
   std::string at_name_2 = index_to_name(at_idx_2);
   int res_no_1 = atom_selection[at_idx_1]->GetSeqNum();
   int res_no_2 = atom_selection[at_idx_2]->GetSeqNum();
   std::string chain_id_1 = atom_selection[at_idx_1]->GetChainID();
   std::string chain_id_2 = atom_selection[at_idx_2]->GetChainID();

   clipper::Coord_orth co_1 = index_to_pos(at_idx_1);
   clipper::Coord_orth co_2 = index_to_pos(at_idx_2);
   double dist = clipper::Coord_orth::length(co_1, co_2);

   std::cout << "spin-scores " << std::setw(4) << at_idx_1 << " ";
   if (using_test_model)
      std::cout << at_name_1 << " " << res_no_1 << " " << chain_id_1 << " ";
   std::cout << std::setw(4) << at_idx_2 << " ";
   if (using_test_model)
      std::cout << at_name_2 << " " << res_no_2 << " " << chain_id_2 << " ";
   std::cout // << " dist: " << dist << "    "
      << " score: "
      << std::right << std::setw(8) << std::setprecision(3) << std::fixed
      << score.second.spin_score << "  " 
      << co_1.x() << " " << co_1.y() << " " << co_1.z() << "    "
      << co_2.x() << " " << co_2.y() << " " << co_2.z() << "    "
      << ca_1_flag << " " << ca_2_flag << " " << consecutive_flag
      << std::endl;
}

void
coot::trace::make_connection_map(const std::vector<std::pair<unsigned int, scored_node_t> > &scores) {
   
   double good_enough_score = 0.5; // overwritten

   if (scores.size() > 0)
      good_enough_score = scores.back().second.spin_score;
   std::cout << "DEBUG:: set good_enough_score to make connection map: "
	     << good_enough_score << std::endl;

   // add them to the connection map:
   // 
   std::vector<scored_node_t>::const_iterator it;

   for (unsigned int i=0; i<scores.size(); i++) {
      if (scores[i].second.spin_score  > good_enough_score) {
	 it = std::find(fwd_connection_map[scores[i].first].begin(),
			fwd_connection_map[scores[i].first].end(),
			scores[i].second);
	 if (it == fwd_connection_map[scores[i].first].end())
	    fwd_connection_map[scores[i].first].push_back(scores[i].second);
      }
   }

   // make the reverse connection map from that: "a node was (fwd)
   // connected by these nodes"
   //
   std::map<unsigned int, std::vector<scored_node_t> >::const_iterator itm;

   if (false) { 
      for (itm = fwd_connection_map.begin(); itm != fwd_connection_map.end(); itm++) {
	 const std::vector<scored_node_t> &v = itm->second;
	 for (it=v.begin(); it!=v.end(); it++) {
	    scored_node_t reverse_node(itm->first, it->spin_score, it->alpha, true);
	    bck_connection_map[it->atom_idx].push_back(reverse_node);
	 }
      }
   }

   for (unsigned int i=0; i<scores.size(); i++) {
      if (scores[i].second.spin_score > good_enough_score) {
	 scored_node_t node(scores[i].first,
			    scores[i].second.spin_score,
			    scores[i].second.alpha,
			    true);
	 if (false) // debug
	    std::cout << "bck_connection_map["
		      << std::setw(4) << scores[i].second.atom_idx
		      << "] push back node of " << std::setw(4)
		      << node.atom_idx << std::endl;
	 
	 bck_connection_map[scores[i].second.atom_idx].push_back(node);
      }
   }

   if (false) { // debugging output
      for (itm=fwd_connection_map.begin(); itm!=fwd_connection_map.end(); itm++) {
	 std::cout << "fwd-map: " << itm->first << " ";
	 if (using_test_model)
	    std::cout << atom_selection[itm->first]->name << " " 
		      << atom_selection[itm->first]->GetSeqNum() << ":  ";
	 for (unsigned int jj=0; jj<itm->second.size(); jj++) { 
	    std::cout << " " << itm->second[jj].atom_idx;
	    if (using_test_model)
	       std::cout << " "
			 << atom_selection[itm->second[jj].atom_idx]->name << " "
			 << atom_selection[itm->second[jj].atom_idx]->GetSeqNum() << " ";
	 }
	 std::cout << std::endl;

	 const std::vector<scored_node_t> &v = bck_connection_map[itm->first];
	 std::cout << "bck-map: " << itm->first << " ";
	 if (using_test_model)
	    std::cout << atom_selection[itm->first]->name << " "
		      << atom_selection[itm->first]->GetSeqNum() << ":  ";
	 
	 for (unsigned int jj=0; jj<v.size(); jj++) {
	    std::cout << " " << v[jj].atom_idx;
	    if (using_test_model)
	       std::cout << " " << atom_selection[v[jj].atom_idx]->name
			 << " " << atom_selection[v[jj].atom_idx]->GetSeqNum() << " ";
	 }
	 std::cout << std::endl;
      }
   }

   bool debug = false;
   
   // debug::
   //
   if (debug) { 
      for (itm=fwd_connection_map.begin(); itm!=fwd_connection_map.end(); itm++) {
	 int res_no = atom_selection[itm->first]->GetSeqNum();
	 std::cout << "fwd-map " << itm->first <<  "  ["
		   << std::setw(2) << itm->second.size() << "] ";
	 if (using_test_model)
	    std::cout << atom_selection[itm->first]->name << " " << res_no << " ";
	 std::cout << index_to_pos(itm->first).format() << "  ";
	 for (unsigned int jj=0; jj<itm->second.size(); jj++) {
	    int res_no_n = atom_selection[itm->second[jj].atom_idx]->GetSeqNum();
	    std::cout << "  " << itm->second[jj].atom_idx << " ";
	    if (using_test_model)
	       std::cout << atom_selection[itm->second[jj].atom_idx]->name << " "
			 << res_no_n << " ";
	 }
	 std::cout << std::endl;
      }

      for (itm=bck_connection_map.begin(); itm!=bck_connection_map.end(); itm++) {
	 int res_no = atom_selection[itm->first]->GetSeqNum();
	 std::cout << "bck-map " << itm->first <<  "  ["
		   << std::setw(2) << itm->second.size() << "] ";
	 if (using_test_model)
	    std::cout << atom_selection[itm->first]->name << " " << res_no << " ";
	 std::cout << index_to_pos(itm->first).format() << "  ";
	 for (unsigned int jj=0; jj<itm->second.size(); jj++) {
	    int res_no_n = atom_selection[itm->second[jj].atom_idx]->GetSeqNum();
	    std::cout << "  " << itm->second[jj].atom_idx << " ";
	    if (using_test_model)
	       std::cout << atom_selection[itm->second[jj].atom_idx]->name << " "
			 << res_no_n << " ";
	 }
	 std::cout << std::endl;
      }
   }
}

void
coot::trace::set_frag_score_crit(const std::vector<std::pair<unsigned int, scored_node_t> > &scores) {

   if (scores.size()) { 
      double sum = 0;
      for (unsigned int i=0; i<scores.size(); i++) { 
	 sum += scores[i].second.spin_score;
      }

      frag_score_crit = 2 * sum/double(scores.size());
   }
}


coot::minimol::fragment
coot::trace::make_fragment(const std::pair<unsigned int, coot::scored_node_t> &scored_node,
			   int res_no_base,
			   std::string chain_id) {

   if (false)
      std::cout << "   make_fragment() called with  node: " << scored_node.first
		<< " scored_node.idx: " << scored_node.second.atom_idx << std::endl;
   
   clipper::Coord_orth pos_1 = index_to_pos(scored_node.first);
   clipper::Coord_orth pos_2 = index_to_pos(scored_node.second.atom_idx);

   if (scored_node.second.reversed_flag)
      std::swap(pos_1, pos_2);

   // as in spin_score():
   // 
   clipper::Coord_orth arb(0,0,1);
   clipper::Coord_orth diff_p(pos_2 - pos_1);
   clipper::Coord_orth diff_p_unit(diff_p.unit());

   clipper::Coord_orth perp(clipper::Coord_orth::cross(arb, diff_p));
   clipper::Coord_orth perp_unit(perp.unit());

   clipper::Coord_orth double_perp(clipper::Coord_orth::cross(diff_p, perp));
   clipper::Coord_orth double_perp_unit(double_perp.unit());

   double diff_p_len = sqrt(diff_p.lengthsq());

   double alpha = scored_node.second.alpha;

   // << " using alpha " << alpha << std::endl;

   double along_CA_CA_pt_O = 1.73; // the C is lower down than the O.
   double along_CA_CA_pt_C = 1.46;
   double along_CA_CA_pt_N = 2.44;

   double ideal_peptide_length = 3.81;

   // we don't want the peptide to be scrunched up on one side of a
   // "long" peptide... let the atom positions expand along a long peptide.

   double f_ca_ca_o = along_CA_CA_pt_O * diff_p_len/ideal_peptide_length;
   double f_ca_ca_c = along_CA_CA_pt_C * diff_p_len/ideal_peptide_length;
   double f_ca_ca_n = along_CA_CA_pt_N * diff_p_len/ideal_peptide_length;

   clipper::Coord_orth rel_line_pt_O(diff_p_unit * f_ca_ca_o + perp_unit * 1.66);
   clipper::Coord_orth rel_line_pt_C(diff_p_unit * f_ca_ca_c + perp_unit * 0.48);
   clipper::Coord_orth rel_line_pt_N(diff_p_unit * f_ca_ca_n - perp_unit * 0.47);

   clipper::Coord_orth p_N = util::rotate_around_vector(diff_p_unit,
							pos_1 + rel_line_pt_N,
							pos_1, alpha);
   
   clipper::Coord_orth p_O = util::rotate_around_vector(diff_p_unit,
							pos_1 + rel_line_pt_O,
							pos_1, alpha);
   
   clipper::Coord_orth p_C = util::rotate_around_vector(diff_p_unit,
							pos_1 + rel_line_pt_C,
							pos_1, alpha);

   minimol::residue r1(res_no_base,   "ALA");
   minimol::residue r2(res_no_base+1, "ALA");
   minimol::atom at_O   (" O  ", "O", p_O,   "", 10);
   minimol::atom at_C   (" C  ", "C", p_C,   "", 10);
   minimol::atom at_N   (" N  ", "N", p_N,   "", 10);
   minimol::atom at_CA_1(" CA ", "C", pos_1, "", 10);
   minimol::atom at_CA_2(" CA ", "C", pos_2, "", 10);

   r1.addatom(at_CA_1);
   r1.addatom(at_C);
   r1.addatom(at_O);
   
   r2.addatom(at_N);
   r2.addatom(at_CA_2);

   // minimol::fragment f(chain_id, false); // use non-expanding constructor
   minimol::fragment f(chain_id);
   
   try { 
      f.addresidue(r1, false);
      f.addresidue(r2, false);
      if (false) { // debug
	 // is f in a good state now?
	 minimol::molecule m_debug(f);
	 std::string  node_string_1 = util::int_to_string(scored_node.first);
	 std::string  node_string_2 = util::int_to_string(scored_node.second.atom_idx);
	 std::string fn_ = "just-a-pair-" + node_string_1 + "-" + node_string_2 + ".pdb";
	 m_debug.write_file(fn_, 10);
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR:: rte on addresidue() in make_fragment() "
		<< rte.what() << std::endl;
   }
   return f;
}

coot::minimol::fragment
// coot::minimol::residue
coot::trace::make_residue(const std::pair<unsigned int, coot::scored_node_t> &scored_node_a,
			  const std::pair<unsigned int, coot::scored_node_t> &scored_node_b,
			  int res_no_base,
			  std::string chain_id) const {

   // The returned residue will be residue b not residue a.
   // Residue a makes the N of the returned residue. Residue b make
   //

   if (false)
      std::cout << " make_residue() called with scored_node_a.first: " << scored_node_a.first
		<< " scored_node_a.atom_idx: " << scored_node_a.second.atom_idx << " b: "
		<< scored_node_b.first << " " << scored_node_b.second.atom_idx
		<< std::endl;

   clipper::Coord_orth pos_1 = index_to_pos(scored_node_a.first);
   clipper::Coord_orth pos_2 = index_to_pos(scored_node_a.second.atom_idx);
   clipper::Coord_orth pos_3 = index_to_pos(scored_node_b.second.atom_idx);

   if (scored_node_a.second.reversed_flag) {
      std::swap(pos_1, pos_2);
      std::cout << "WARNING:: FIXME - can't deal with reversed peptides ATM " << std::endl;
      minimol::fragment dummy_frag;
      return dummy_frag;
   }

   // as in spin_score():
   //
   clipper::Coord_orth arb(0,0,1);
   clipper::Coord_orth diff_p_1(pos_2 - pos_1);
   clipper::Coord_orth diff_p_1_unit(diff_p_1.unit());
   clipper::Coord_orth diff_p_2(pos_3 - pos_2);
   clipper::Coord_orth diff_p_2_unit(diff_p_2.unit());

   clipper::Coord_orth perp_1(clipper::Coord_orth::cross(arb, diff_p_1));
   clipper::Coord_orth perp_1_unit(perp_1.unit());
   clipper::Coord_orth perp_2(clipper::Coord_orth::cross(arb, diff_p_2));
   clipper::Coord_orth perp_2_unit(perp_2.unit());

   clipper::Coord_orth double_perp_1(clipper::Coord_orth::cross(diff_p_1, perp_1));
   clipper::Coord_orth double_perp_1_unit(double_perp_1.unit());
   clipper::Coord_orth double_perp_2(clipper::Coord_orth::cross(diff_p_2, perp_2));
   clipper::Coord_orth double_perp_2_unit(double_perp_2.unit());

   double diff_p_1_len = sqrt(diff_p_1.lengthsq());
   double diff_p_2_len = sqrt(diff_p_2.lengthsq());

   double alpha_1 = scored_node_a.second.alpha;
   double alpha_2 = scored_node_b.second.alpha;

   // << " using alpha " << alpha << std::endl;

   double along_CA_CA_pt_O = 1.73; // the C is lower down than the O.
   double along_CA_CA_pt_C = 1.46;
   double along_CA_CA_pt_N = 2.44;

   double ideal_peptide_length = 3.81;

   // we don't want the peptide to be scrunched up on one side of a
   // "long" peptide... let the atom positions expand along a long peptide.

   double f_ca_ca_o_1 = along_CA_CA_pt_O * diff_p_1_len/ideal_peptide_length;
   double f_ca_ca_c_1 = along_CA_CA_pt_C * diff_p_1_len/ideal_peptide_length;
   double f_ca_ca_n_1 = along_CA_CA_pt_N * diff_p_1_len/ideal_peptide_length;
   double f_ca_ca_o_2 = along_CA_CA_pt_O * diff_p_2_len/ideal_peptide_length;
   double f_ca_ca_c_2 = along_CA_CA_pt_C * diff_p_2_len/ideal_peptide_length;
   double f_ca_ca_n_2 = along_CA_CA_pt_N * diff_p_2_len/ideal_peptide_length;

   clipper::Coord_orth rel_line_pt_O(diff_p_1_unit * f_ca_ca_o_1 + perp_1_unit * 1.66);
   clipper::Coord_orth rel_line_pt_C(diff_p_1_unit * f_ca_ca_c_1 + perp_1_unit * 0.48);
   clipper::Coord_orth rel_line_pt_N(diff_p_1_unit * f_ca_ca_n_1 - perp_1_unit * 0.47);

   clipper::Coord_orth rel_line_pt_C_2(diff_p_2_unit * f_ca_ca_c_2 + perp_2_unit * 0.48);
   clipper::Coord_orth rel_line_pt_O_2(diff_p_2_unit * f_ca_ca_o_2 + perp_2_unit * 1.66);

   clipper::Coord_orth p_N_2 = util::rotate_around_vector(diff_p_1_unit,
							  pos_1 + rel_line_pt_N,
							  pos_1, alpha_1);

   clipper::Coord_orth p_O_1 = util::rotate_around_vector(diff_p_1_unit,
							  pos_1 + rel_line_pt_O,
							  pos_1, alpha_1);

   clipper::Coord_orth p_O_2 = util::rotate_around_vector(diff_p_2_unit,
							  pos_2 + rel_line_pt_O_2,
							  pos_2, alpha_2);

   clipper::Coord_orth p_C_1 = util::rotate_around_vector(diff_p_1_unit,
							  pos_1 + rel_line_pt_C,
							  pos_1, alpha_1);

   clipper::Coord_orth p_C_2 = util::rotate_around_vector(diff_p_2_unit,
							  pos_2 + rel_line_pt_C_2,
							  pos_2, alpha_2);

   minimol::residue r1(res_no_base,   "ALA");
   minimol::residue r2(res_no_base+1, "ALA");
   minimol::atom at_O_1 (" O  ", "O", p_O_1,   "", 10);
   minimol::atom at_C_1 (" C  ", "C", p_C_1,   "", 10);
   minimol::atom at_N   (" N  ", "N", p_N_2,   "", 10);
   minimol::atom at_CA_1(" CA ", "C", pos_1, "", 10);
   minimol::atom at_CA_2(" CA ", "C", pos_2, "", 10);
   minimol::atom at_O_2 (" O  ", "O", p_O_2, "", 10);
   minimol::atom at_C_2 (" C  ", "C", p_C_2, "", 10);

   r1.addatom(at_CA_1);
   r1.addatom(at_C_1);
   r1.addatom(at_O_1);

   r2.addatom(at_N);
   r2.addatom(at_CA_2);
   r2.addatom(at_C_2);
   r2.addatom(at_O_2);

   if (false) { // debugging
      minimol::atom at_Cnext(" Cnext", "C", pos_3, "", 30);
      r2.addatom(at_Cnext);
   }

   std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position(r2);
   if (cbeta_info.first) {
      r2.addatom(" CB ", " C", cbeta_info.second, "", 1.0, 30);
   }

   minimol::fragment f("A");
   f.addresidue(r1, false);
   f.addresidue(r2, false);
   return f;

   // return r2;
}

bool
coot::trace::nice_fit(const minimol::fragment &f) const {

   return true;
}

bool
coot::trace::nice_fit(const minimol::residue &r1, const minimol::residue &r2) const {

   double s = get_fit_score(r1, r2);
   return true;

}

double
coot::trace::get_fit_score(const minimol::residue &r1, const minimol::residue &r2) const {

   return 0;
}

coot::minimol::fragment
coot::trace::merge_fragments(const coot::minimol::fragment &f1,
			     const coot::minimol::fragment &f2) const {

   bool debug = false;
   minimol::fragment merged;

   if (debug) { 
      std::cout << "---- in merge_fragments() f1: " << f1 << std::endl;
      std::cout << "---- in merge_fragments() f2: " << f2 << std::endl;

      for (int ires=f1.min_res_no(); ires<=f1.max_residue_number(); ires++)
	 for (unsigned int iat=0; iat<f1[ires].atoms.size(); iat++)
	    std::cout << "   mol1: res: " << ires << " atom " << f1[ires].atoms[iat]
		      << std::endl;
      for (int ires=f2.min_res_no(); ires<=f2.max_residue_number(); ires++)
	 for (unsigned int iat=0; iat<f2[ires].atoms.size(); iat++)
	    std::cout << "   mol2: res: " << ires << " atom " << f2[ires].atoms[iat]
		      << std::endl;
   }

   merged = f1;
   for (int ires=f2.min_res_no(); ires<=f2.max_residue_number(); ires++) {

      // Does this residue exist in merged already?
      // If so, the we want to replace or add to the atoms that are already there
      // If not, add this residue

      bool have_residue = true;
      try {
	 const minimol::residue &r = merged[ires];
      }
      catch (const std::runtime_error &rte) {
	 have_residue = false;
	 std::cout << "caught no residue for " << ires << std::endl;
      } 

      if (! have_residue) {
	 minimol::residue r = f2[ires];
	 r.seqnum = ires;
	 std::cout << "merge_residues() calling addresidue() for ires " << ires
		   << " and name " << r.name << std::endl;
	 merged.addresidue(r, false);
      } else {

	 if (f2[ires].atoms.size()) { 
	    // atoms to be added or replace existing atoms

	    if (debug)
	       std::cout << "atoms to be added or replace existing atoms for ires "
			 << ires << " and internal seqnum " << f2[ires].seqnum
			 << std::endl;

	    // OK, residue merges[ires] does exist but it's empty

	    merged[ires].seqnum = ires;
	    merged[ires].name = f2[ires].name;

	    for (unsigned int jat=0; jat<f2[ires].atoms.size(); jat++) {

	       bool found = false;
	       // current atoms
	       for (unsigned int iat=0; iat<merged[ires].atoms.size(); iat++) {
	       
		  if (merged[ires][iat].name == f2[ires][jat].name) {
		     merged[ires][iat] = f2[ires][jat];
		     if (debug)
			std::cout << "replaced atom in ires  " << ires << ": "
				  << merged[ires][iat] << std::endl;
		     found = true;
		     break;
		  }
	       }
	       if (! found) {
		  if (debug)
		     std::cout << "add atom " << f2[ires][jat] << std::endl;
		  merged[ires].addatom(f2[ires][jat]);
	       }
	    }
	 }
      }
   }
   return merged;
} 


void
coot::trace::trace_graph() {

   // fill interesting_trees

   std::map<unsigned int, std::vector<scored_node_t> >::const_iterator it;

   std::cout << "in ---- trace_graph() --- tr is of size "
	     << fwd_connection_map.size() << std::endl;

   for(it=fwd_connection_map.begin(); it!=fwd_connection_map.end(); it++) {

      std::vector<scored_node_t> path;
      unsigned int iv=it->first;

      // check here if iv is a leaf

      if (it->second.size() == 1) { 
	 // std::cout << "---- trace start with node " << iv << std::endl;
	 scored_node_t leaf(iv, 0, 0);
	 next_vertex(path, 0, leaf);
      }
   }
   sort_filter_interesting_trees();
}

void
coot::trace::sort_filter_interesting_trees() {

   std::sort(interesting_trees.begin(),
	     interesting_trees.end(),
	     sort_trees_by_length);

   if (false)  // debuggging. 
      for (unsigned int itree=0; itree<interesting_trees.size(); itree++)
	 std::cout << "   itree " << itree << " " << interesting_trees[itree].size()
		   << std::endl;
}



std::pair<unsigned int, coot::scored_node_t> 
coot::trace::spin_score(unsigned int idx_1, unsigned int idx_2) const {

   double score = 0;

   mmdb::Atom *at_1 = atom_selection[idx_1]; 
   mmdb::Atom *at_2 = atom_selection[idx_2];

   const clipper::Coord_orth pos_1 = index_to_pos(idx_1);

   const clipper::Coord_orth pos_2 = index_to_pos(idx_2);

   // draw a line between pos_1 and pos_2
   
   // find a point A that is 1.56 down the stick and 1.57 away from the line
   // find a point B that is 1.9  down the stick and 1.91 away from the line
   // find a point C that is 1.9  down the stick and 1.91 away from the line
   //      in the opposite direction to B.

   // spin this points A, B and C around the pos_1 - pos_2 line and the score is
   // rho(A) - rho(B) - rho(C)

   // Note to self: also "above" and "below" the peptide plane we expect to have
   // little to no density - so that should be added to the scoring system too.

   clipper::Coord_orth arb(0,0,1);
   clipper::Coord_orth diff_p(pos_2 - pos_1);
   clipper::Coord_orth diff_p_unit(diff_p.unit());

   clipper::Coord_orth perp(clipper::Coord_orth::cross(arb, diff_p));
   clipper::Coord_orth perp_unit(perp.unit());

   clipper::Coord_orth double_perp(clipper::Coord_orth::cross(diff_p, perp));
   clipper::Coord_orth double_perp_unit(double_perp.unit());

   double along_CA_CA_pt_O = 1.53; // the C is lower down than the O.
   double along_CA_CA_pt_for_perp = 2.33;
                                         
                                         
   double along_CA_CA_pt_N = 2.5;
   double ideal_peptide_length = 3.81;

   // we don't want the peptide to be scrunched up on one side of a
   // "long" peptide... let the atom positions expand along a long peptide. 
   double diff_p_len = sqrt(diff_p.lengthsq());
   double f_ca_ca_o = along_CA_CA_pt_O * diff_p_len/ideal_peptide_length;
   // double f_ca_ca_c = along_CA_CA_pt_C * diff_p_len/ideal_peptide_length;
   double f_ca_ca_n = along_CA_CA_pt_N * diff_p_len/ideal_peptide_length;
   double f_ca_ca_pt_for_perp = along_CA_CA_pt_for_perp * diff_p_len/ideal_peptide_length;

   // clipper::Coord_orth rel_line_pt_C(diff_p_unit * f_ca_ca_c + perp_unit * 0.7);
   // clipper::Coord_orth rel_line_pt_N(diff_p_unit * f_ca_ca_n - perp_unit * 0.5);

   // there is good density 1.9A away from the mid-line in the direction of the CO
   // (at the O). 
   // there is little density 3.7A away from the mid-line in the direction of the CO
   clipper::Coord_orth rel_line_pt_O(      diff_p_unit * f_ca_ca_o + perp_unit * 1.89);
   clipper::Coord_orth rel_line_pt_O_low(  diff_p_unit * f_ca_ca_o + perp_unit * 3.9);
   clipper::Coord_orth rel_line_pt_CO_anti(diff_p_unit * f_ca_ca_o - perp_unit * 0.5);
   clipper::Coord_orth rel_line_pt_N(      diff_p_unit * f_ca_ca_n - perp_unit * 0.3);
   clipper::Coord_orth rel_line_pt_perp1(diff_p_unit * f_ca_ca_pt_for_perp  + double_perp_unit * 1.85);
   clipper::Coord_orth rel_line_pt_perp2(diff_p_unit * f_ca_ca_pt_for_perp  - double_perp_unit * 1.72);


   // Idea from Andrea T (originally from George, I understand)
   // 
   // There should be density for a hydrogen bond acceptor in the
   // extension of the direction to the N from the CA-CA line.
   // However, we need to spin-search this a bit because it's often somewhat off axis.
   // 
   // Currently doesn't work.
   clipper::Coord_orth rel_line_pt_N_accpt(diff_p_unit * f_ca_ca_n - perp_unit * 3.0);
   clipper::Coord_orth rel_line_pt_N_accpt_off(diff_p_unit * (f_ca_ca_n + 0.5) - perp_unit * 3.0);
   
   float rho_at_1 = util::density_at_point(xmap, pos_1);
   float rho_at_2 = util::density_at_point(xmap, pos_2);

   // the mid-point between CAs should have density too.
   clipper::Coord_orth pt_mid(pos_1 * 0.50 + pos_2 * 0.5);
   float rho_mid = util::density_at_point(xmap, pt_mid);

   int n_steps = 36;
   float best_score = -999;

   float rho_CO_best      = -999; // for testing scoring
   float rho_CO_low_best  = -999;
   float rho_CO_anti_best = -999;
   float rho_N_best       = -999;
   float rho_perp_1_best  = -999;
   float rho_perp_2_best  = -999;
   double alpha_best = -1; // negative indicates this was not set - which is a bad thing.

   // // these can be optimized with machine learning?
   // float scale_CO       =  0.5;
   // float scale_CO_low   = -0.6;
   // float scale_CO_anti  = -0.1;
   // float scale_perp     = -0.7;
   // float scale_mid      =  1.6;
   // float scale_non_line =  1.0;
   // float scale_N        = -0.0;

   // unsigned int idx_test_1 = 1228;
   // unsigned int idx_test_2 = 1238;
   // unsigned int idx_test_1 = 101;
   // unsigned int idx_test_2 = 109;
   
   unsigned int idx_test_1 = 131;
   unsigned int idx_test_2 = 139;
   
   for (int i=0; i< int(n_steps); i++) { 
      double alpha = 2 * M_PI * double(i)/double(n_steps);

      // direction position orig-shift angle
      // 
      clipper::Coord_orth p_CO = util::rotate_around_vector(diff_p_unit,
							    pos_1 + rel_line_pt_O,
							    pos_1, alpha);
      
      clipper::Coord_orth p_CO_low = util::rotate_around_vector(diff_p_unit,
								pos_1 + rel_line_pt_O_low,
								pos_1, alpha);
      
      clipper::Coord_orth p_CO_anti = util::rotate_around_vector(diff_p_unit,
								 pos_1 + rel_line_pt_CO_anti,
								 pos_1, alpha);
      
      clipper::Coord_orth p_N = util::rotate_around_vector(diff_p_unit,
							   pos_1 + rel_line_pt_N,
							   pos_1, alpha);
      
      clipper::Coord_orth p_2 = util::rotate_around_vector(diff_p_unit,
							   pos_1 + rel_line_pt_perp1,
							   pos_1, alpha);
      clipper::Coord_orth p_3 = util::rotate_around_vector(diff_p_unit,
							   pos_1 + rel_line_pt_perp2,
							   pos_1, alpha);

      // clipper::Coord_orth p_N_acceptor_1 = util::rotate_round_vector(diff_p_unit,
      // 								     pos_1 + rel_line_pt_N_accpt,
      // 								     pos_1, alpha);
      // clipper::Coord_orth p_N_acceptor_2 = util::rotate_round_vector(diff_p_unit,
      // 								     pos_1 + rel_line_pt_N_accpt,
      // 								     pos_1, (alpha + 0.26));
      // clipper::Coord_orth p_N_acceptor_3 = util::rotate_round_vector(diff_p_unit,
      // 								     pos_1 + rel_line_pt_N_accpt,
      // 								     pos_1, (alpha - 0.26));
      // clipper::Coord_orth p_N_acceptor_4 = util::rotate_round_vector(diff_p_unit,
      // 								     pos_1 + rel_line_pt_N_accpt_off,
      // 								     pos_1, alpha);
      // clipper::Coord_orth p_N_acceptor_5 = util::rotate_round_vector(diff_p_unit,
      // 								     pos_1 + rel_line_pt_N_accpt_off,
      // 								     pos_1, (alpha + 0.26));
      // clipper::Coord_orth p_N_acceptor_6 = util::rotate_round_vector(diff_p_unit,
      // 								     pos_1 + rel_line_pt_N_accpt_off,
      // 								     pos_1, (alpha - 0.26));
      
      
      float rho_CO      = util::density_at_point(xmap, p_CO);
      float rho_CO_low  = util::density_at_point(xmap, p_CO_low);
      float rho_CO_anti = util::density_at_point(xmap, p_CO_anti);
      float rho_N       = util::density_at_point(xmap, p_N);
      float rho_perp_1  = util::density_at_point(xmap, p_2);
      float rho_perp_2  = util::density_at_point(xmap, p_3);

      // float rho_acceptor_1 = util::density_at_point(xmap, p_N_acceptor_1);
      // float rho_acceptor_2 = util::density_at_point(xmap, p_N_acceptor_2);
      // float rho_acceptor_3 = util::density_at_point(xmap, p_N_acceptor_3);
      // float rho_acceptor_4 = util::density_at_point(xmap, p_N_acceptor_4);
      // float rho_acceptor_5 = util::density_at_point(xmap, p_N_acceptor_5);
      // float rho_acceptor_6 = util::density_at_point(xmap, p_N_acceptor_6);

      // float rho_acceptor_best = rho_acceptor_1;
      // if (rho_acceptor_2 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_2;
      // if (rho_acceptor_3 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_3;
      // if (rho_acceptor_4 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_4;
      // if (rho_acceptor_5 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_5;
      // if (rho_acceptor_6 > rho_acceptor_best) rho_acceptor_best = rho_acceptor_6;
      
      float this_score =
	 scale_CO      * rho_CO      +
	 scale_CO_low  * rho_CO_low  + 
	 scale_CO_anti * rho_CO_anti +
	 scale_N       * rho_N       +
	 scale_perp    * rho_perp_1  +
	 scale_perp    * rho_perp_2;       
      // scale_N_accpt * rho_acceptor_best


      if (idx_1 == idx_test_1 && idx_2 == idx_test_2) {
	 std::cout << "debug_pos:: CO     " << p_CO.x() << " " << p_CO.y() << " " << p_CO.z()
		   << " " << rho_CO << std::endl;
	 std::cout << "debug_pos:: CO_low " << p_CO_low.x() << " " << p_CO_low.y()
		   << " " << p_CO_low.z() << " " << rho_CO_low << std::endl;
	 std::cout << "debug_pos:: CO_anti " << p_CO_anti.x() << " " << p_CO_anti.y()
		   << " " << p_CO_anti.z() << " " << rho_CO_anti << std::endl;
	 std::cout << "debug_pos:: perp-1 " << p_2.x() << " " << p_2.y() << " " << p_2.z()
		   << " " << rho_perp_1 << std::endl;
	 std::cout << "debug_pos:: perp-2 " << p_3.x() << " " << p_3.y() << " " << p_3.z()
		   << " " << rho_perp_2 << std::endl;
	 std::cout << "debug_pos:: N " << p_N.x() << " " << p_N.y() << " " << p_N.z()
		   << " " << rho_N << std::endl;
      } 
      
      if (this_score > best_score) { 
	 best_score       = this_score;
	 rho_CO_best      = rho_CO;
	 rho_CO_low_best  = rho_CO_low;
	 rho_CO_anti_best = rho_CO_anti;
	 rho_N_best       = rho_N;
	 rho_perp_1_best  = rho_perp_1;
	 rho_perp_2_best  = rho_perp_2;
	 alpha_best       = alpha;
      }
   }

   bool output_density_values = false;
   if (output_density_values) {
      std::cout << "debug-rho:: CO     "  << rho_CO_best      << " "
		<< "debug-rho:: CO_low "  << rho_CO_low_best  << " "
		<< "debug-rho:: CO_anti " << rho_CO_anti_best << " "
		<< "debug-rho:: perp-1 "  << rho_perp_1_best  << " "
		<< "debug-rho:: perp-2 "  << rho_perp_2_best  << " "
		<< "debug-rho:: N "       << rho_N_best << std::endl;
   }

   
   float non_line_equal_density_penalty_1 = rho_at_1 + rho_at_2 - 2 * rho_mid;
   float non_line_equal_density_penalty = 
      - non_line_equal_density_penalty_1 * non_line_equal_density_penalty_1 /
      map_rmsd;

   best_score += scale_mid * rho_mid;
   best_score += scale_non_line * non_line_equal_density_penalty;

   scored_node_t best_node(idx_2, best_score, alpha_best);
   // for testing
   if (using_test_model) {
      std::string atom_name_1(at_1->name);
      std::string atom_name_2(at_2->name);
      if (atom_name_1 == " CA ") {
	 if (atom_name_2 == " CA ") {
	    if ((at_1->GetSeqNum() + 1) == at_2->GetSeqNum())
	       best_node.udd_flag = true;
	 }
      }
   }
   return std::pair<unsigned int, scored_node_t> (idx_1, best_node);
}

void
coot::trace::print_interesting_trees() const {

   for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
      std::cout << "interesting tree " << itree << ": ";
      for (unsigned int j=0; j<interesting_trees[itree].size(); j++) { 
	 const unsigned int &idx = interesting_trees[itree][j].atom_idx;
	 int res_no = atom_selection[idx]->GetSeqNum();
	 std::cout << "  " << idx;
	 if (using_test_model)
	    std::cout << " (" << atom_selection[idx]->name << " " << res_no << ")";
      }
      std::cout << std::endl;
   }

   minimol::molecule m;
   int offset = 0;

   if (false) { 
      for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
	 minimol::fragment f(util::int_to_string(itree));
	 for (unsigned int j=0; j<interesting_trees[itree].size(); j++) {
	    minimol::residue r(j + offset, "ALA");
	    minimol::atom at(" CA ", "C", index_to_pos(interesting_trees[itree][j].atom_idx), "", 10);
	    r.addatom(at);
	    f.addresidue(r, false);
	 }
	 // offset += interesting_trees[itree].size();
	 m.fragments.push_back(f);
      }
      std::cout << "writing interesting.pdb " << std::endl;
      m.write_file("interesting.pdb", 20);
   }
} 

void
coot::trace::optimize_weights(mmdb::Manager *mol_in) {

   using_test_model = true; // makes sense only when use use template pdb
                            // for atom seeding
   
   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atom_pairs_within_distance(mol_in, 3.81, 1);

   double s[8], s_orig[8];
   s[0] = scale_CO;
   s[1] = scale_CO_low;
   s[2] = scale_CO_anti;
   s[3] = scale_perp;
   s[4] = scale_mid;
   s[5] = scale_non_line;
   s[6] = scale_N;
   s[7] = scale_N_accpt;
   for (unsigned int i=0; i<8; i++)
      s_orig[i] = s[i];

   std::vector<std::pair<unsigned int, scored_node_t> > scores =
      spin_score_pairs(apwd);
   set_scales(s);
   double ks = ks_test(scores);
   // The unchanged set
   std::cout << "ks: " << ks << "  ";
   for (unsigned int i=0; i<8; i++) std::cout << " " << s[i];
   std::cout << std::endl;
   

   unsigned int n_factors = 21;

   unsigned int i_scale = 7; // don't forget to set/reset the "moving" scale
   
   for (unsigned int ifact=0; ifact<n_factors; ifact++) { 
      double f_1 = (double(ifact) - double(n_factors-1) * 0.5); // -10 to 10
      double f_2 = 1 + 0.2 * f_1; // -1 to 3 // -3 to 5
      std::cout << "f_2: " << f_2 << std::endl;
      s[i_scale] = s_orig[i_scale] * f_2;
	 
      set_scales(s);

      std::vector<std::pair<unsigned int, scored_node_t> > spin_scores = spin_score_pairs(apwd);
      double ks = ks_test(spin_scores);

      std::cout << "ks: " << ks << "  ";
      for (unsigned int i=0; i<8; i++) std::cout << " " << s[i];
      std::cout << std::endl;
   }
}

void
coot::trace::test_model(mmdb::Manager *mol_in) {

   using_test_model = true; // makes sense only when use use template pdb
                                        // for atom seeding
   
   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atom_pairs_within_distance(mol_in, 3.81, 1);

   std::vector<std::pair<unsigned int, scored_node_t> > scores = spin_score_pairs(apwd);

   unsigned int n_top = 1000; // top-scoring spin-score pairs for tracing
   
   std::sort(scores.begin(), scores.end(), scored_node_t::sort_pair_scores);
   if (scores.size() > n_top) {
      scores.resize(n_top);
   } else {
      n_top = scores.size();
   }
   
   make_connection_map(scores);

   set_frag_score_crit(scores);

   int n_top_fragments = 100;
   //

   std::vector<std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> >
      frag_store;


#pragma omp parallel for // currently we write out the fragments, not store them
   
   for (int i=0; i<n_top_fragments; i++) {

      // int i = 2;

      // int tid = omp_get_thread_num();
      // std::cout << "OMP::: " << tid << std::endl;

      std::vector<scored_node_t> start_path;
      scored_node_t dummy_first;
      dummy_first.atom_idx = scores[i].first;
      start_path.push_back(dummy_first);
      
      std::string chain_id = frag_idx_to_chain_id(i);
      std::cout << "----------- test_model() starting point number  " << i
		<< " atom_idx: " << scores[i].first << " node: " << scores[i].second 
		<< " chain_id " << chain_id << std::endl;

      int res_no_base = scores[i].first;

      res_no_base = i;

      // actually, we *have* the first spin pair but it gets thrown
      // away and "rediscovered" when build_2_choose_1() does just
      // that (and hopefully picks the idx scores[i].second.atom_ids).
      // 

      dir_t dir = FORWARDS;
      
      std::cout << "------------------------------- follow forwards ---------"
       		<< std::endl;


      std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> ff = 
	 follow_fragment(scores[i].first, start_path, res_no_base, chain_id, dir);

      std::cout << "------------------------------- follow backwards --------"
       		<< std::endl;
      dir = BACKWARDS;
      std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> bf = 
	 follow_fragment(scores[1].second.atom_idx, start_path, 1000+res_no_base, chain_id, dir);

      add_replace_reject(frag_store, ff);
      add_replace_reject(frag_store, bf);

      std::cout << "For i_top_frag " << i << " we have " << frag_store.size()
		<< " fragment " << std::endl;
      
   }


   if (false) { 
      for (unsigned int ifrag=0; ifrag<frag_store.size(); ifrag++) {
	 std::string  fn = "frag-store-" + util::int_to_string(ifrag) + ".pdb";
	 frag_to_pdb(frag_store[ifrag].second, fn);
      }
   }


   // Rama terminal addition/refine trace.
   //
   std::pair<float, float> mv = util::mean_and_variance(xmap);
   protein_geometry geom;
   geom.init_standard();
   geom.remove_planar_peptide_restraint();
   multi_peptide(frag_store, geom, mv);
   
}

// if return.first is more than 1, then we have something interesting
// 
std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment>
coot::trace::follow_fragment(unsigned int atom_idx, const std::vector<scored_node_t> &path_in,
			     int res_no_base,
			     const std::string &chain_id,
			     dir_t dir) {

   std::vector<scored_node_t> local_path = path_in;
   minimol::fragment running_fragment;
   unsigned int atom_idx_in = atom_idx;
   
   bool status = true;
   
   while (status) {
      // this function should return the atom_idx of the node that
      // is a best neighbour of atom_idx
      //

      int res_no = res_no_base + local_path.size();
      if (dir == BACKWARDS)
	 res_no = res_no_base - local_path.size();
      
      std::pair<bool, scored_node_t> sn =
	 build_2_choose_1(atom_idx, local_path, res_no, chain_id, dir);
      status = sn.first;
      if (false) // debug
	 std::cout << "From node " << atom_idx << " got next: status: "
		   << sn.first << " and index: " << sn.second.atom_idx << std::endl;
      if (sn.first) {
	 std::pair<unsigned int, scored_node_t> p(atom_idx, sn.second);

	 // f needs to get the "right" residue numbers.
	 minimol::fragment f = make_fragment(p, res_no, chain_id);
	 
	 running_fragment = merge_fragments(running_fragment, f);

	 // setup for next round
	 atom_idx = sn.second.atom_idx;
	 local_path.push_back(sn.second);
	 
      } else {
	 // what's in running fragment?

	 if (false) { 
	    minimol::molecule m_debug(running_fragment);
	    std::string  node_string = util::int_to_string(res_no_base);
	    std::string fn_ = "running-frag-node-" + node_string + ".pdb";
	    std::cout << ":::::: fragment output: " << fn_ << std::endl;
	    if (! m_debug.is_empty())
	       m_debug.write_file(fn_, 10);
	 }
      }
   }

   // debug: count the residues of running_fragment
   unsigned int n_fragment_residues = 0;
   for (int ires=running_fragment.min_res_no();
	   ires<=running_fragment.max_residue_number();
	   ires++) { 
	 if (running_fragment[ires].atoms.size())
	    n_fragment_residues++;
   }
   
   std::cout << "debug:: follow_fragment(): returning running_fragment of size: "
	     <<  local_path.size()  << " " << n_fragment_residues << std::endl;

   add_cbetas(&running_fragment);
   std::pair<std::vector<scored_node_t>, minimol::fragment> p(local_path, running_fragment);
   
   return p;
}

// modify frag
void
coot::trace::add_cbetas(minimol::fragment *frag) {

   for (int ires=frag->first_residue(); ires<=frag->max_residue_number(); ires++) {
      if ((*frag)[ires].atoms.size() > 2) {
	 std::pair<bool, clipper::Coord_orth> cbeta_info = cbeta_position((*frag)[ires]);
	 if (cbeta_info.first) {
	    (*frag)[ires].addatom(" CB ", " C", cbeta_info.second, "", 1.0, 20);
	 }
      }
   }
}



void
coot::trace::add_replace_reject(std::vector<std::pair<std::vector<scored_node_t>, minimol::fragment> > &frag_store, const std::pair<std::vector<scored_node_t>, minimol::fragment> &trial) const {

   // is this like anything else?
   unsigned int n_for_match = 4;
   unsigned int n_match_for_replace = 4; 
   bool add_it = true;

   unsigned int n_match = 0;
   for (unsigned int i=0; i<frag_store.size(); i++) { 
      for (unsigned int j=0; j<frag_store[i].first.size(); j++) { 
	 for (unsigned int k=0; k<trial.first.size(); k++) { 
	    if (trial.first[k] == frag_store[i].first[j]) {
	       n_match ++;
	    } 
	 }
      }
      if (n_match >= n_for_match) {
	 if (trial.first.size() < frag_store[i].first.size()) {
	    // there is already something like and longer than trial in frag_store
	    add_it = false;
	    break;
	 } 
      }
   }



   if (add_it) {

      // delete things already in frag_store that are like trial

      frag_store.erase(std::remove_if(frag_store.begin(),
				      frag_store.end(),
				      frag_store_eraser(trial, &xmap, n_match_for_replace)),
		       frag_store.end());
      frag_store.push_back(trial);
   } 

} 


std::pair<bool, coot::scored_node_t>
coot::trace::build_2_choose_1(unsigned int atom_idx,
			      const std::vector<scored_node_t> &start_path,
			      int res_no_base,
			      const std::string &chain_id,
			      dir_t dir) {

   bool status = false;

   if (false) { // debug

      std::cout << "build_2_choose_1(): called with atom_idx " << atom_idx;
      if (using_test_model)
	 std::cout << " ( "
		   << atom_selection[atom_idx]->name << " "
		   << atom_selection[atom_idx]->GetSeqNum() << " "
		   << atom_selection[atom_idx]->GetChainID() << " ) ";
      std::cout << " start_path: ";
      for (unsigned int jj=0; jj<start_path.size(); jj++) { 
	 std::cout << " " << start_path[jj];
	 if (using_test_model)
	    std::cout << "( "
		      << atom_selection[start_path[jj].atom_idx]->name << " "
		      << atom_selection[start_path[jj].atom_idx]->GetSeqNum() << " "
		      << atom_selection[start_path[jj].atom_idx]->GetChainID() << " ) ";
      }
      std::cout << " and dir " << dir << std::endl;


      std::vector<scored_node_t> nf =
	 get_neighbours_of_vertex_excluding_path(atom_idx, start_path, FORWARDS);
      std::cout << "   INFO:: build_2_choose_1(): node-lev-1 " << atom_idx << " has "
		<< nf.size() << " forward neighbours:" << std::endl;
      for (unsigned int j=0; j<nf.size(); j++)
	 std::cout << "   " << nf[j].atom_idx;
      std::cout << std::endl;

      std::vector<scored_node_t> nb =
	 get_neighbours_of_vertex_excluding_path(atom_idx, start_path, BACKWARDS);

      std::cout << "   INFO:: build_2_choose_1(): node-lev-1 " << atom_idx << " has "
		<< nb.size() << " backward neighbours:" << std::endl;
      for (unsigned int j=0; j<nb.size(); j++)
	 std::cout << "   " << nb[j].atom_idx;
      std::cout << std::endl;
   }  // end debug

   // this functions checks that the neighbors don't have to same atom_idx as
   // those in path 
   //
   std::vector<scored_node_t> neighbs_1 = get_neighbours_of_vertex_excluding_path(atom_idx, start_path, dir);

   if (true) { 
      std::cout << "   INFO:: node-lev-1 " << atom_idx << " has "
		<< neighbs_1.size() << " neighbours:" << std::endl;
      for (unsigned int j=0; j<neighbs_1.size(); j++)
	 std::cout << "   " << neighbs_1[j].atom_idx;
      std::cout << std::endl;
   }
   
   std::vector<indexed_frag_t> frag_store; // chose the best one of these

   for (unsigned int j=0; j<neighbs_1.size(); j++) { 
      std::pair<unsigned int, scored_node_t> p_1(atom_idx, neighbs_1[j]);

      bool ang_1_OK = true;
      if (start_path.size() > 1) {
	 double ang_1 = path_candidate_angle(neighbs_1[j].atom_idx, start_path);

	 if (ang_1 <= 0.8 * M_PI * 0.5)
	    ang_1_OK = false;
	
	 if (false) { // debug
	    unsigned int l = start_path.size();
	    if (using_test_model) 
	       std::cout << "   angle_1 from "
			 << neighbs_1[j].atom_idx << " ( "
			 << atom_selection[neighbs_1[j].atom_idx]->name         << " " 
			 << atom_selection[neighbs_1[j].atom_idx]->GetSeqNum()  << " " 
			 << atom_selection[neighbs_1[j].atom_idx]->GetChainID() << " ) "
			 << start_path[l-1] << " ( "
			 << atom_selection[start_path[l-1].atom_idx]->name         << " " 
			 << atom_selection[start_path[l-1].atom_idx]->GetSeqNum()  << " " 
			 << atom_selection[start_path[l-1].atom_idx]->GetChainID() << " ) "
			 << start_path[l-2] << " ( "
			 << atom_selection[start_path[l-2].atom_idx]->name         << " " 
			 << atom_selection[start_path[l-2].atom_idx]->GetSeqNum()  << " " 
			 << atom_selection[start_path[l-2].atom_idx]->GetChainID() << " ) "
			 << " ang: " << ang_1 << " " << " (should be more than "
			 << 0.8 * M_PI * 0.5 << ") " << ang_1_OK << std::endl;
	 } // end debug
      }

      if (ang_1_OK) { 
	 minimol::fragment f_1 = make_fragment(p_1, res_no_base, chain_id);

	 if (nice_fit(f_1)) {

	    std::vector<scored_node_t> local_path = start_path;
	    local_path.push_back(neighbs_1[j]);

	    std::vector<scored_node_t> neighbs_2 =
	       get_neighbours_of_vertex_excluding_path(neighbs_1[j].atom_idx, local_path, dir);

	    if (false) {
	       std::cout << "   INFO:: node-lev-2 " << neighbs_1[j].atom_idx << " has "
			 << neighbs_2.size() << " neighbours:" << std::endl;
	       for (unsigned int jj=0; jj<neighbs_2.size(); jj++)
		  std::cout << "   " << neighbs_2[jj].atom_idx;
	       std::cout << std::endl;
	    }
	 
	    for (unsigned int jj=0; jj<neighbs_2.size(); jj++) {

	       double ang_2 = path_candidate_angle(neighbs_2[jj].atom_idx, local_path);

	       if (false) { 
		  unsigned int l = local_path.size();
		  if (using_test_model) 
		     std::cout << "   angle_2 from "
			       << neighbs_2[jj].atom_idx << " ( "
			       << atom_selection[neighbs_2[jj].atom_idx]->name         << " " 
			       << atom_selection[neighbs_2[jj].atom_idx]->GetSeqNum()  << " " 
			       << atom_selection[neighbs_2[jj].atom_idx]->GetChainID() << " ) "
			       << local_path[l-1] << " ( "
			       << atom_selection[local_path[l-1].atom_idx]->name         << " " 
			       << atom_selection[local_path[l-1].atom_idx]->GetSeqNum()  << " " 
			       << atom_selection[local_path[l-1].atom_idx]->GetChainID() << " ) "
			       << local_path[l-2] << " ( "
			       << atom_selection[local_path[l-2].atom_idx]->name         << " " 
			       << atom_selection[local_path[l-2].atom_idx]->GetSeqNum()  << " " 
			       << atom_selection[local_path[l-2].atom_idx]->GetChainID() << " ) "
			       << " ang: " << ang_2 << " " << " (should be more than "
			       << 0.8 * M_PI * 0.5 << ")" << std::endl;
		  else
		     std::cout << "   angle from "
			       << neighbs_2[jj].atom_idx << " " 
			       << local_path[l-1] << " " << local_path[l-2] << " ang: "
			       << ang_2 << " " << " (should be more than "
			       << 0.8 * M_PI * 0.5 << ")" << std::endl;
		  
	       }
		     
	       if (ang_2 > 0.8 * M_PI * 0.5) {

		  std::pair<unsigned int, scored_node_t> p_2(neighbs_1[j].atom_idx, neighbs_2[jj]);
		  minimol::fragment f_2 = make_fragment(p_2, res_no_base, chain_id);

		  // Now store f1, f2 indexed on neighbs_1[j].atom_idx, neighbs_2[jj].atom_idx
		  // so that we can pick choose the best fit and continue the fragment starting
		  // from neighbs_1[j]

		  indexed_frag_t idx_frag(neighbs_1[j], neighbs_2[jj], f_1, f_2);

		  frag_store.push_back(idx_frag);
	       } else  {
		  // std::cout << "   reject by angle " << neighbs_2[jj] << std::endl;
	       }
	    }
	 } else {
	    std::cout << "   reject by fit " << std::endl;
	 }
      }
   }

   // now choose the best frag in frag_store

   std::cout << ":: building fragment starting at res " << res_no_base
	     << " I have to choose the best of " << frag_store.size()
	     << " frags " << std::endl;

   double best_score = -9999;
   scored_node_t best_scored_node;

   if (frag_store.size() > 0) {
      for (unsigned int i=0; i<frag_store.size(); i++) { 
	 if (frag_store[i].get_score() > best_score) {
	    best_score       = frag_store[i].get_score();
	    best_scored_node = frag_store[i].node_1;
	 }
      }
      // double frag_score_crit = -33.0;  // hack value
      if (best_score > frag_score_crit) {
	 status = true;
      } else {
	 std::cout << "   None of the " << frag_store.size() << " fragments "
		   << "have score better than " << frag_score_crit << std::endl;
      } 
   } else {
      std::cout << "No fragments from which to choose in build_2_choose_1()" << std::endl;
   } 
   
   return std::pair<bool, scored_node_t> (status, best_scored_node);

} 


void
coot::trace::print_tree(const std::vector<unsigned int> &path) const {

   std::cout << "path: ";
   for (unsigned int i=0; i<path.size(); i++) {
      int res_no = atom_selection[path[i]]->GetSeqNum();
      std::cout << "  " << path[i] << " (" << index_to_name(path[i]) << " "
		<< res_no << ")";
   }
   std::cout << std::endl;

   if (false) 
      for (unsigned int i=0; i<path.size(); i++) {
	 const clipper::Coord_orth &pt = index_to_pos(path[i]);
	 std::cout << "long path " << i
		   << " " << pt.x()
		   << " " << pt.y()
		   << " " << pt.z()
		   << std::endl;
   }
}



double
coot::trace::path_candidate_angle(unsigned int candidate_vertex,
				  const std::vector<scored_node_t> &path) const {

   unsigned int l = path.size();
   const clipper::Coord_orth &pt_1 = index_to_pos(candidate_vertex);
   const clipper::Coord_orth &pt_2 = index_to_pos(path[l-1].atom_idx);
   const clipper::Coord_orth &pt_3 = index_to_pos(path[l-2].atom_idx);
   double angle = clipper::Coord_orth::angle(pt_1, pt_2, pt_3);

   return angle;
}


void
coot::trace::add_tree_maybe(const std::vector<scored_node_t> &path) {


   // add this path if there is not already a path that is longer that
   // contains at least n atoms of the same path points.

   bool add_this = true;
   unsigned int n_match_crit = 2;
   unsigned int n_match_for_replace = 2; // at least this number of matches to replace an existing tree

   for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
      unsigned int n_match = 0;

      // do I already have something that's longer than path and
      // similar to path already in interesting_trees? (if so, set
      // add_this to false).

      // the rejection should take account of the length of the
      // overlapping paths.
      //
      // if they are in the same direction and the length of the paths
      // is much longer than the overlaps, then they should be merged.
      // How can that happen though? The path explorer (next_vertex)
      // should find the merged tree.
      
      if (path.size() < interesting_trees[itree].size()) {
	 for (unsigned int i=0; i<interesting_trees[itree].size(); i++) {
	    for (unsigned int j=0; j<path.size(); j++) {
	       if (path[j] == interesting_trees[itree][i]) {
		  n_match ++;
	       }
	    }

	    if (n_match >= n_match_crit) {
	       add_this = false;
	       break;
	    }
	 }
      }
   }

   if (0) { 
      std::cout << "add status " << add_this << " for tree of length " << path.size()
		<< " because " << interesting_trees.size() << " trees of lengths ";
      for (unsigned int ii=0; ii<interesting_trees.size(); ii++) { 
	 std::cout << "  " << interesting_trees[ii].size();
      }
      std::cout << std::endl;
   }

   if (add_this) {
      interesting_trees.erase(std::remove_if(interesting_trees.begin(),
					     interesting_trees.end(),
					     trace_path_eraser(path, n_match_for_replace)),
			      interesting_trees.end());
      interesting_trees.push_back(path);
   } 
}



// Rama terminal addition/refine trace.
void
coot::trace::multi_peptide(const std::vector<std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> > &frag_store, const coot::protein_geometry &geom, std::pair<float, float> &mv) {

   unsigned int n_top = 20;
   if (frag_store.size() < n_top) n_top = frag_store.size();

   float b_factor = 20;
   int n_trials = 3000; // for terminal residue addition

   std::cout << "multi_peptide(): we have " << frag_store.size() << " fragments in the store " << std::endl;

   for (unsigned int i=0; i<n_top; i++) {

      { 
	 std::string node_string = util::int_to_string(i);
	 std::string fn_ = "from-multi-peptide:frag-store-" + node_string + ".pdb";
	 frag_to_pdb(frag_store[i].second, fn_);
      }

      int min_res_no     = frag_store[i].second.first_residue();
      int max_res_no     = frag_store[i].second.max_residue_number();
      int n_terminal_res = frag_store[i].second.first_residue() + 1;
      int c_terminal_res = frag_store[i].second.max_residue_number() -1; // this should contain C,O, CA and N.

      bool in_range = false;
      if (n_terminal_res <= max_res_no)
	 if (n_terminal_res >= min_res_no)
	    if (c_terminal_res <= max_res_no)
	       if (c_terminal_res >= min_res_no)
		  if (n_terminal_res < c_terminal_res) 
		     in_range = true;

      if (in_range) { 

	 try {
      
	    unsigned int n_atoms_in_N_res = frag_store[i].second[n_terminal_res].atoms.size();
	    unsigned int n_atoms_in_C_res = frag_store[i].second[c_terminal_res].atoms.size();

	    std::cout << "   multi_peptide(): fragstore frag[" << i << "] N-terminal residue with seqnum "
		      <<  n_terminal_res << " has " << n_atoms_in_N_res << " atoms " << std::endl;
	    std::cout << "   multi_peptide(): fragstore frag[" << i << "] C-terminal residue with seqnum "
		      <<  c_terminal_res << " has " << n_atoms_in_C_res << " atoms " << std::endl;

	    if (n_atoms_in_N_res > 2) {
	       bool debugging = false;
	       mmdb::Residue *res_p            = frag_store[i].second[n_terminal_res].make_residue();
	       // 20180406-PE does this work? not tested
	       mmdb::Residue *res_downstream_p = frag_store[i].second[n_terminal_res-1].make_residue();
	       minimol::fragment f = multi_build_N_terminal_ALA(res_p,
								res_downstream_p,
								frag_store[i].second.fragment_id,
								b_factor,
								n_trials,
								geom,
								xmap, mv, debugging);
      
	       std::cout << "multi-build on N on frag_store fragment index " << i
			 << " made a fragment of size " << f.n_filled_residues() << std::endl;

	       { 
		  std::string  node_string = util::int_to_string(i);
		  std::string fn_ = "from-multi-peptide:multi-build-from-N-" + node_string + ".pdb";
		  frag_to_pdb(f, fn_);
	       }
	    
	    }
	    if (n_atoms_in_C_res > 2) {
	       bool debugging = false;
	       mmdb::Residue *res_p = frag_store[i].second[c_terminal_res].make_residue();
	       // 20180406-PE does this work? not tested
	       mmdb::Residue *res_upstream_p = frag_store[i].second[c_terminal_res-1].make_residue();
	       minimol::fragment f = multi_build_C_terminal_ALA(res_p,
								res_upstream_p,
								frag_store[i].second.fragment_id,
								b_factor,
								n_trials,
								geom,
								xmap, mv, debugging);
      
	       std::cout << "multi-build on C on frag_store fragment index " << i
			 << " made a fragment of size " << f.n_filled_residues() << std::endl;

	       {
		  std::string  node_string = util::int_to_string(i);
		  std::string fn_ = "from-multi-peptide:multi-build-from-C-" + node_string + ".pdb";
		  frag_to_pdb(f, fn_);
	       }
	    }
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "skip " << rte.what() << std::endl;
	 }
      }
   }
}
