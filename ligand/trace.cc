
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
   std::vector<std::pair<unsigned int, unsigned int> > apwd =
      atoms_pairs_within_distance(flood_mol, 3.81, 0.2);
   spin_score_pairs(apwd);

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
coot::trace::atoms_pairs_within_distance(const coot::minimol::molecule &flood_mol,
					 double trans_dist,
					 double trans_dist_variation) {
   
   std::vector<std::pair<unsigned int, unsigned int> > v;
   mmdb::Manager *mol = flood_mol.pcmmdbmanager();
   if (mol) {
      int SelHnd = mol->NewSelection();
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
		<< std::endl;

      sas = flood_mol.select_atoms_serial();
      std::cout << "sas length " << sas.size() << std::endl;

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

	 std::cout << "found " << n_contacts << " distance pairs " << std::endl;
	 std::cout << "made  " << v.size() << " distance pairs " << std::endl;
      }
      
      mol->DeleteSelection(SelHnd);
   }
   return v;
} 

void
coot::trace::spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd) const {

   // for (unsigned int i=0; i<apwd.size(); i++) { 
   //    std::cout << "   " << apwd[i].first << " " << apwd[i].second << std::endl;
   // }

   tree<trace_node> tr;
   tree<trace_node>::iterator it;
   tree<trace_node>::iterator top = tr.begin();
   double good_enough_score = 1.499;
   
   for (unsigned int i=0; i<apwd.size(); i++) {
      const double &d_f = spin_score(apwd[i].first,  apwd[i].second);
      const double &d_b = spin_score(apwd[i].second, apwd[i].first);
      std::cout << apwd[i].first << " " << apwd[i].second << " "
		<< d_f << " " << d_b << std::endl;
      if (d_f > good_enough_score || d_b > good_enough_score) {

	 trace_node n(apwd[i].second, d_f, d_b);

	 // is apwd[i].second already in the tree?
	 //
	 bool found = false;
	 unsigned int n_tree = 0;
	 for (it=tr.begin(); it!=tr.end(); it++) {
	    n_tree++;
	    if (it->atom_idx == apwd[i].first) {
	       // yes it was
	       found = true;
	       tr.append_child(it, n);
	       std::cout << "found" << std::endl;
	       break;
	    }
	 }
	 if (! found) {
	    std::cout << "new start for node with idx " << apwd[i].first
		      << " and " << n_tree << " in tree for index "
		      << apwd[i].second << std::endl;
	    tr.insert(tr.begin(), n); 
	 } 
      }
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

   clipper::Coord_orth rel_line_pt_A(diff_p_unit * 1.56 +        perp_unit * 1.57);
   clipper::Coord_orth rel_line_pt_B(diff_p_unit * 1.9  + double_perp_unit * 1.91);
   clipper::Coord_orth rel_line_pt_C(diff_p_unit * 1.9  - double_perp_unit * 1.91);

   int n_steps = 180;
   float best_score = -999;
   for (unsigned int i=0; i<n_steps; i++) { 
      double alpha = 2 * M_PI * double(i)/double(n_steps);
      clipper::Coord_orth p_1 = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_A,
							  pos_1, alpha);
      clipper::Coord_orth p_2 = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_B,
							  pos_1, alpha);
      clipper::Coord_orth p_3 = util::rotate_round_vector(diff_p_unit,
							  pos_1 + rel_line_pt_C,
							  pos_1, alpha);
      float this_rho_1 = util::density_at_point(xmap, p_1);
      float this_rho_2 = util::density_at_point(xmap, p_2);
      float this_rho_3 = util::density_at_point(xmap, p_3);

      float this_score = this_rho_1 - this_rho_2 - this_rho_3;

      if (false) 
	 std::cout << "debug_pos:: " << p_1.x() << " " << p_1.y() << " " << p_1.z()
		   << " " << this_rho_1 << " " << this_rho_2 << " " << this_rho_3
		   << std::endl;
      
      if (this_score > best_score) 
	 best_score = this_score;
   }

   return best_score;
   
}
