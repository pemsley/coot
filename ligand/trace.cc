
#include <iomanip>
#include "coot-utils/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "trace.hh"
#include "ligand.hh"
#include "coot-utils/tree.hh"

// First we will flood the map, then look for atom pair that are 2.8
// to 4.8A apart and spin a search the density

coot::trace::trace(const clipper::Xmap<float>& xmap_in) {

   xmap = xmap_in;
   flood_atom_mask_radius = 1.5;
   rmsd_cut_off = 1.0;
   
}

void
coot::trace::action() {

   minimol::molecule flood_mol = get_flood_molecule();

   sas = flood_mol.select_atoms_serial();
   std::cout << "sas length " << sas.size() << std::endl;
   
   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atoms_pairs_within_distance(flood_mol, 3.81, 1.0);
   spin_score_pairs(apwd);
   trace_graph();
   print_interesting_trees();

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
   std::string output_pdb = "out.pdb";
   water_mol.write_file(output_pdb, 30.0);
   lig.output_map("find-waters-masked-flooded.map");
   
   return water_mol;
   
}

std::vector<std::pair<unsigned int, unsigned int> >
coot::trace::atoms_pairs_within_distance(mmdb::Manager *mol,
					 double trans_dist,
					 double trans_dist_variation) {

   std::vector<std::pair<unsigned int, unsigned int> > v;
   if (mol) {
      int SelHnd = mol->NewSelection(); // d
      mol->SelectAtoms(SelHnd, 0, "*",
		       mmdb::ANY_RES, // starting resno, an int
		       "*", // any insertion code
		       mmdb::ANY_RES, // ending resno
		       "*", // ending insertion code
		       "*", // any residue name
		       "*", // atom name
		       "*", // elements
		       "");

      mmdb::PAtom *sel_atoms = 0;
      int n_sel_atoms;
      mol->GetSelIndex(SelHnd, sel_atoms, n_sel_atoms);

      std::cout << "selected " << n_sel_atoms << " for distance pair check"
		<< " (this should be the same as sas length)" << std::endl;

      int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM, "index");
      if (uddHnd<0)  {
	 std::cout << " atom bonding registration failed.\n";
      } else {
	 for (int i=0; i< n_sel_atoms; i++) { 
	    sel_atoms[i]->PutUDData(uddHnd, i);
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
	 mol->SeekContacts(sel_atoms, n_sel_atoms,
			   sel_atoms, n_sel_atoms,
			   local_dist_min, local_dist_max,
			   0,        // seqDist 0 -> in same res also
			   pscontact, n_contacts,
			   0, &my_matt, i_contact_group);

	 if (n_contacts > 0) {
	    if (pscontact) {
	       for (int i=0; i<n_contacts; i++) {
		  if (pscontact[i].id1 < pscontact[i].id2) { 
		     mmdb::Atom *at_1 = sel_atoms[pscontact[i].id1];
		     mmdb::Atom *at_2 = sel_atoms[pscontact[i].id2];
		     int idx_1, idx_2;
		     at_1->GetUDData(uddHnd, idx_1);
		     at_2->GetUDData(uddHnd, idx_2);
		     std::pair<unsigned int, unsigned int> p(idx_1, idx_2);
		     v.push_back(p);
		  }
	       }
	    }
	 }

	 std::cout << "found " << n_contacts << " distance pairs " << std::endl;
	 std::cout << "made  " << v.size() << " distance pairs " << std::endl;
      }
      
      mol->DeleteSelection(SelHnd);
   }
   return v;

}
   


std::vector<std::pair<unsigned int, unsigned int> >
coot::trace::atoms_pairs_within_distance(const coot::minimol::molecule &flood_mol,
					 double trans_dist,
					 double trans_dist_variation) {
   
   std::vector<std::pair<unsigned int, unsigned int> > v;
   mmdb::Manager *mol = flood_mol.pcmmdbmanager();
   return atoms_pairs_within_distance(mol, trans_dist, trans_dist_variation);
}

void
coot::trace::spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd) {

   // apwd : atom (index) pairs within distance
   
   double good_enough_score = 0.7;

   std::vector<unsigned int>::const_iterator it; 
   
   for (unsigned int i=0; i<apwd.size(); i++) {
      const unsigned &at_idx_1 = apwd[i].first;
      const unsigned &at_idx_2 = apwd[i].second;
      const double &d_f = spin_score(at_idx_1, at_idx_2);
      const double &d_b = spin_score(at_idx_2, at_idx_1);

      // debugging output
      //
      //
      bool ca_1_flag = false;
      bool ca_2_flag = false;
      bool consecutive_flag = false;

      // indexing looks OK.
      if (sas[at_idx_1]->name == " CA ") ca_1_flag = true;
      if (sas[at_idx_2]->name == " CA ") ca_2_flag = true;

      if (ca_1_flag && ca_2_flag) {
	 std::cout << "  occs " << sas[at_idx_2]->occupancy
		   << "       " << sas[at_idx_1]->occupancy << std::endl;
	 float occ_delta = sas[at_idx_1]->occupancy - sas[at_idx_2]->occupancy;
	 if (fabsf(occ_delta) < 1.1) {
	    if (occ_delta < 0)
	       consecutive_flag = true;
	 } 
      } 

      std::cout << "spin-scores "
		<< apwd[i].first << " " << apwd[i].second << " "
		<< ca_1_flag << " " << ca_2_flag << " " << consecutive_flag << " "
		<< d_f << " " << d_b << std::endl;

      
      if ((d_f > good_enough_score) || (d_b > good_enough_score)) {

	 it = std::find(tr[apwd[i].first].begin(),
			tr[apwd[i].first].end(),
			apwd[i].second);
	 if (it == tr[apwd[i].first].end())
	    tr[apwd[i].first].push_back(apwd[i].second);

	 // and reverse indexing
	       
	 it = std::find(tr[apwd[i].second].begin(),
			tr[apwd[i].second].end(),
			apwd[i].first);
	 if (it == tr[apwd[i].second].end())
	    tr[apwd[i].second].push_back(apwd[i].first);

      }
   }

   std::map<unsigned int, std::vector<unsigned int> >::const_iterator itm;
   for (itm=tr.begin(); itm!=tr.end(); itm++) {
      std::cout << "map " << itm->first <<  "  ["
		<< std::setw(2) << itm->second.size() << "] ";
      std::cout << sas[itm->first]->name << " ";
      std::cout << sas[itm->first]->pos.format() << "  ";
      for (unsigned int jj=0; jj<itm->second.size(); jj++) { 
	 std::cout << "  " << itm->second[jj];
      }
      std::cout << std::endl;
   }
   
}

void
coot::trace::trace_graph() {

   std::map<unsigned int, std::vector<unsigned int> >::const_iterator it;

   std::cout << "in trace_graph tr is of size " << tr.size() << std::endl;

   for(it=tr.begin(); it!=tr.end(); it++) {

      std::vector<unsigned int> path;
      unsigned int iv=it->first;
      path = next_vertex(path, 0, iv);

      // print the path
      std::cout << iv << " path size " << path.size() << " ";
      for (unsigned int ii=0; ii<path.size(); ii++) { 
	 std::cout << path[ii] << "  ";
      }
      std::cout << std::endl;
   }
}



// We presume that idx_1 and idx_2 are less than the size of sas
// 
double
coot::trace::spin_score(unsigned int idx_1, unsigned int idx_2) const {

   double score = 0;

   minimol::atom *at_1 = sas[idx_1];
   minimol::atom *at_2 = sas[idx_2];

   const clipper::Coord_orth &pos_1 = at_1->pos;
   const clipper::Coord_orth &pos_2 = at_2->pos;

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

   // these can be optimized with machine learning?
   float scale_CO       =  0.6;
   float scale_CO_low   = -0.7;
   float scale_perp     = -0.7;
   float scale_mid      =  1.4;
   float scale_non_line =  1.0;

   // unsigned int idx_test_1 = 1228;
   // unsigned int idx_test_2 = 1238;
   
   unsigned int idx_test_1 = 101;
   unsigned int idx_test_2 = 109;
   
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

   return best_score;
   
}

void
coot::trace::print_interesting_trees() const {

   for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
      std::cout << "interesting tree " << itree << ": ";
      for (unsigned int j=0; j<interesting_trees[itree].size(); j++) { 
	 const unsigned int &idx = interesting_trees[itree][j];
	 int res_no = nearbyint(sas[idx]->occupancy);
	 std::cout << "  " << idx << " " << " (" << sas[idx]->name << " " << res_no << ")";
      }
      std::cout << std::endl;
   }


   minimol::molecule m;
   minimol::fragment f("A");
   int offset = 0;
   for (unsigned int itree=0; itree<interesting_trees.size(); itree++) {
      for (unsigned int j=0; j<interesting_trees[itree].size(); j++) {
	 minimol::residue r(j + offset, "HOH");
	 minimol::atom at(" O  ", "O", sas[interesting_trees[itree][j]]->pos, "", 10);
	 std::cout << "added atom " << at << std::endl;
	 r.addatom(at);
	 f.addresidue(r, false);
      }
      offset += interesting_trees[itree].size();
   }
   m.fragments.push_back(f);
   std::cout << "writing interesting.pdb " << std::endl;
   m.write_file("interesting.pdb", 20);
   

} 

void
coot::trace::test_model(mmdb::Manager *mol) {

   mol_for_sas = minimol::molecule(mol);

   for (unsigned int ich=0; ich<mol_for_sas.fragments.size(); ich++) {
      for (int ires=mol_for_sas.fragments[ich].min_res_no();
	   ires<=mol_for_sas.fragments[ich].max_residue_number(); ires++) {
	 for (unsigned int iat=0; iat<mol_for_sas.fragments[ich][ires].atoms.size(); iat++) { 
	    mol_for_sas.fragments[ich][ires][iat].occupancy = ires + 100*ich;
	 }
      }
   }
   
   sas = mol_for_sas.select_atoms_serial();
   
   std::cout << "sas length " << sas.size() << std::endl;

   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atoms_pairs_within_distance(mol, 3.81, 1);

   spin_score_pairs(apwd);
   trace_graph();
   print_interesting_trees();

} 
