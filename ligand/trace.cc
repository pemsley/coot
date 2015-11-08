
#include <iomanip>
#include <algorithm>
#include "coot-utils/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "trace.hh"
#include "ligand.hh"
// #include "coot-utils/tree.hh"
#include "analysis/stats.hh"

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

   add_atom_names_in_map_output = false; 

}

void
coot::trace::action() {

   minimol::molecule flood_mol = get_flood_molecule();

   mmdb::Manager *mol = flood_mol.pcmmdbmanager();

   if (mol) { 
      std::vector<std::pair<unsigned int, unsigned int> > apwd =
	 atom_pairs_within_distance(mol, 3.81, 1.0);
      spin_score_pairs(apwd);
      trace_graph();
      print_interesting_trees();
   }

}

coot::minimol::molecule
coot::trace::get_flood_molecule() const {

   coot::ligand lig;
   lig.set_cluster_size_check_off();
   lig.set_chemically_sensible_check_off();
   lig.set_sphericity_test_off();
	       
   lig.set_map_atom_mask_radius(flood_atom_mask_radius);
   lig.set_water_to_protein_distance_limits(10.0, 1.5);
   
   lig.import_map_from(xmap);
   
   lig.flood2(rmsd_cut_off);
   coot::minimol::molecule water_mol = lig.water_mol();
   // debug
   std::string output_pdb = "flood-mol.pdb";
   water_mol.write_file(output_pdb, 30.0);
   lig.output_map("find-waters-masked-flooded.map");
   
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

      atom_selection = 0; // member data - cleared on descruction
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
   


void
coot::trace::spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd) {

   // apwd : atom (index) pairs within distance

   unsigned int n_top = 500; // top-scoring spin-score pairs for tracing
   
   double good_enough_quartile = 0.9;  // only allow the top 90% of scores to be in a trace/tree
   double good_enough_score = 0.5; // overwritten

   std::vector<std::pair<unsigned int, scored_node_t> > scores(apwd.size()*2);
   
   for (unsigned int i=0; i<apwd.size(); i++) {
      const unsigned &at_idx_1 = apwd[i].first;
      const unsigned &at_idx_2 = apwd[i].second;

      scores[2*i ]  = spin_score(at_idx_1, at_idx_2);
      scores[2*i+1] = spin_score(at_idx_2, at_idx_1);

      // debugging output
      //
      //
      bool ca_1_flag = false;
      bool ca_2_flag = false;
      bool consecutive_flag = false;

      // indexing looks OK.
      if (index_to_name(at_idx_1) == " CA ") ca_1_flag = true;
      if (index_to_name(at_idx_2) == " CA ") ca_2_flag = true;

      if (ca_1_flag && ca_2_flag) {
	 int resno_delta =
	    atom_selection[at_idx_2]->GetSeqNum() -
	    atom_selection[at_idx_1]->GetSeqNum();
	 if (resno_delta == 1)
	    consecutive_flag = true;
      } 

      if (true) { // debugging
	 const std::string &at_name_1 = index_to_name(at_idx_1);
	 const std::string &at_name_2 = index_to_name(at_idx_2);
	 int res_no_1 = atom_selection[at_idx_1]->GetSeqNum();
	 int res_no_2 = atom_selection[at_idx_2]->GetSeqNum();

	 clipper::Coord_orth co_1 = index_to_pos(at_idx_1);
	 clipper::Coord_orth co_2 = index_to_pos(at_idx_2);
	 double dist = clipper::Coord_orth::length(co_1, co_2);

	 std::cout << "spin-scores "
		   << at_idx_1 << " ";
	 if (add_atom_names_in_map_output)
	    std::cout << at_name_1 << " " << res_no_1 << " "
	    << at_idx_2 << " ";
	 if (add_atom_names_in_map_output)
	    std::cout << at_name_2 << " " << res_no_2 << " ";
	 std::cout << dist << " " 
		   << ca_1_flag << " " << ca_2_flag << " " << consecutive_flag << " "
		   << std::endl;
      }
   }

   if (scores.size() > n_top) {
      std::sort(scores.begin(), scores.end(), scored_node_t::sort_pair_scores);
      scores.resize(n_top);
   } else {
      n_top = scores.size();
   }

   if (scores.size() > 0)
      good_enough_score = scores[n_top-1].second.spin_score;
   
   std::cout << "DEBUG:: set good_enough_score to " << good_enough_score << std::endl;
   

   // add them to the connection map:
   // 
   std::vector<scored_node_t>::const_iterator it;

   for (unsigned int i=0; i<n_top; i++) {
      
      if (scores[i].second.spin_score  > good_enough_score) {

	 it = std::find(tr[scores[i].first].begin(),
			tr[scores[i].first].end(),
			scores[i].second);
	 if (it == tr[scores[i].first].end())
	    tr[scores[i].first].push_back(scores[i].second);
      }
   }

   std::map<unsigned int, std::vector<scored_node_t> >::const_iterator itm;
   for (itm=tr.begin(); itm!=tr.end(); itm++) {
      int res_no = atom_selection[itm->first]->GetSeqNum();
      std::cout << "map " << itm->first <<  "  ["
		<< std::setw(2) << itm->second.size() << "] ";
      if (add_atom_names_in_map_output)
	 std::cout << atom_selection[itm->first]->name << " " << res_no << " ";
      std::cout << index_to_pos(itm->first).format() << "  ";
      for (unsigned int jj=0; jj<itm->second.size(); jj++) {
	 int res_no_n = atom_selection[itm->second[jj].atom_idx]->GetSeqNum();
	 std::cout << "  " << itm->second[jj].atom_idx << " ";
	 if (add_atom_names_in_map_output)
	    std::cout << atom_selection[itm->second[jj].atom_idx]->name << " " <<res_no_n << " ";
      }
      std::cout << std::endl;
   }
   
}

void
coot::trace::trace_graph() {

   // fill interesting_trees

   std::map<unsigned int, std::vector<scored_node_t> >::const_iterator it;

   std::cout << "in ---- trace_graph() --- tr is of size "
	     << tr.size() << std::endl;

   for(it=tr.begin(); it!=tr.end(); it++) {

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

   clipper::Coord_orth arb(0,0,1);
   clipper::Coord_orth diff_p(pos_2 - pos_1);
   clipper::Coord_orth diff_p_unit(diff_p.unit());

   clipper::Coord_orth perp(clipper::Coord_orth::cross(arb, diff_p));
   clipper::Coord_orth perp_unit(perp.unit());

   clipper::Coord_orth double_perp(clipper::Coord_orth::cross(diff_p, perp));
   clipper::Coord_orth double_perp_unit(double_perp.unit());

   // there is good   density 1.9A away from the mid-line in the direction of the CO (at the O)
   // there is little density 3.7A away from the mid-line in the direction of the CO
   clipper::Coord_orth rel_line_pt_A(    diff_p_unit * 1.5  +        perp_unit * 1.9);
   clipper::Coord_orth rel_line_pt_A_low(diff_p_unit * 1.5  +        perp_unit * 3.7);
   clipper::Coord_orth rel_line_pt_B(    diff_p_unit * 2.3  + double_perp_unit * 1.8);
   clipper::Coord_orth rel_line_pt_C(    diff_p_unit * 2.3  - double_perp_unit * 1.8);

   float rho_at_1 = util::density_at_point(xmap, pos_1);
   float rho_at_2 = util::density_at_point(xmap, pos_2);

   // the mid-point between CAs should have density too.
   clipper::Coord_orth pt_mid(pos_1 * 0.50 + pos_2 * 0.5);
   float rho_mid = util::density_at_point(xmap, pt_mid);

   int n_steps = 36;
   float best_score = -999;
   float rho_CO_best = -999; // for testing scoring
   float rho_CO_low_best = -999;
   float rho_2_best  = -999;
   float rho_3_best  = -999;
   double alpha_best = 0;

   // these can be optimized with machine learning?
   float scale_CO       =  0.6;
   float scale_CO_low   = -0.7;
   float scale_perp     = -0.7;
   float scale_mid      =  1.4;
   float scale_non_line =  1.0;

   // unsigned int idx_test_1 = 1228;
   // unsigned int idx_test_2 = 1238;
   // unsigned int idx_test_1 = 101;
   // unsigned int idx_test_2 = 109;
   
   unsigned int idx_test_1 = 0; // stop output
   unsigned int idx_test_2 = 0;
   
   for (int i=0; i< int(n_steps); i++) { 
      double alpha = 2 * M_PI * double(i)/double(n_steps);

                                                         // direction position orig-shift angle
      clipper::Coord_orth p_CO = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_A,
							  pos_1, alpha);
      
      clipper::Coord_orth p_CO_low = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_A_low,
							  pos_1, alpha);
      
      clipper::Coord_orth p_2 = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_B,
							  pos_1, alpha);
      clipper::Coord_orth p_3 = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_C,
							  pos_1, alpha);
      
      float rho_CO     = util::density_at_point(xmap, p_CO);
      float rho_CO_low = util::density_at_point(xmap, p_CO_low);
      float rho_2      = util::density_at_point(xmap, p_2);
      float rho_3      = util::density_at_point(xmap, p_3);

      float this_score =
	 scale_CO *     rho_CO     +
	 scale_CO_low * rho_CO_low + 
	 scale_perp *   rho_2      +
	 scale_perp *   rho_3;

      if (idx_1 == idx_test_1 && idx_2 == idx_test_2) {
	 std::cout << "debug_pos:: CO     " << p_CO.x() << " " << p_CO.y() << " " << p_CO.z()
		   << " " << rho_CO << std::endl;
	 std::cout << "debug_pos:: CO_low " << p_CO_low.x() << " " << p_CO_low.y()
		   << " " << p_CO_low.z() << " " << rho_CO_low << std::endl;
	 std::cout << "debug_pos:: perp-1 " << p_2.x() << " " << p_2.y() << " " << p_2.z()
		   << " " << rho_2 << std::endl;
	 std::cout << "debug_pos:: perp-2 " << p_3.x() << " " << p_3.y() << " " << p_3.z()
		   << " " << rho_3 << std::endl;
      } 
      
      if (this_score > best_score) { 
	 best_score = this_score;
	 rho_CO_best = rho_CO;
	 rho_CO_low_best = rho_CO_low;
	 rho_2_best = rho_2;
	 rho_3_best = rho_3;
	 alpha_best = alpha;
      }
   }

   float non_line_equal_density_penalty_1 = rho_at_1 + rho_at_2 - 2 * rho_mid;
   float non_line_equal_density_penalty = 
      - non_line_equal_density_penalty_1 * non_line_equal_density_penalty_1;

   best_score += scale_mid * rho_mid;
   best_score += scale_non_line * non_line_equal_density_penalty;

   if (idx_1 == idx_test_1 && idx_2 == idx_test_2) {
      std::cout << "debug score-parts: " << rho_at_1 << " " << rho_at_2 << " mid: "
		<< rho_mid << " CO " << rho_CO_best
		<< " perp-1 " << rho_2_best << " perp-2 " << rho_3_best
		<< " = " << best_score << std::endl;
   }

   scored_node_t best_node(idx_2, best_score, alpha_best);
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
	 if (add_atom_names_in_map_output)
	    std::cout << " (" << atom_selection[idx]->name << " " << res_no << ")";
      }
      std::cout << std::endl;
   }


   minimol::molecule m;
   int offset = 0;
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

void
coot::trace::test_model(mmdb::Manager *mol) {

   add_atom_names_in_map_output = true; // makes sense only when use use template pdb
                                        // for atom seeding
   
   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atom_pairs_within_distance(mol, 3.81, 1);

   spin_score_pairs(apwd);
   trace_graph();
   print_interesting_trees();

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
coot::trace::path_candidate_angle(const std::vector<scored_node_t> &path,
				  unsigned int candidate_vertex) const {

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

