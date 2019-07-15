
#include <chrono>

#include "coot-coord-utils.hh"
#include "stack-and-pair.hh"

coot::stack_and_pair::stack_and_pair(mmdb::Manager *mol, int atom_selection_handle) {

   init();
   mmdb::PAtom *sel_atoms = 0;
   int n_sel_atoms;
   mol->GetSelIndex(atom_selection_handle, sel_atoms, n_sel_atoms);
   normal_map = calculate_residue_normals(sel_atoms, n_sel_atoms);

}

coot::stack_and_pair::stack_and_pair(mmdb::Manager *mol,
				     const std::vector<std::pair<bool,mmdb::Residue *> > &residues_vec) {

   init();
   normal_map = calculate_residue_normals(residues_vec);
}


void
coot::stack_and_pair::init() {

   // PDBv3 FIXME
   std::vector<std::string> base_atom_names = {" N1 ", " C2 ", " N3 ", " C4 ",
					       " C5 ", " C6 ", " N7 ", " C8 ", " N9 "};
   for (std::size_t i=0; i<base_atom_names.size(); i++)
      base_atom_name_set.insert(base_atom_names[i]);

   angle_crit = clipper::Util::d2rad(30.0);

}

bool
coot::stack_and_pair::contains_nucleic_acid(mmdb::Atom **SelAtom, int n_sel_atoms) {

   return true;
}

std::pair<bool, clipper::Coord_orth>
coot::stack_and_pair::get_base_normal(mmdb::Residue *residue_p) const {

   std::pair<bool, clipper::Coord_orth> r(false, clipper::Coord_orth(0,0,0));

   std::vector<clipper::Coord_orth> v;
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string atom_name(at->name);
      if (base_atom_name_set.find(atom_name) != base_atom_name_set.end()) {
	 v.push_back(co(at));
      }
   }
   if (v.size() > 2) {
      lsq_plane_info_t lsq(v);
      r.first = true;
      r.second = clipper::Coord_orth(lsq.normal().unit());
   }
   return r;
}

std::vector<std::string>
coot::stack_and_pair::get_base_atom_names(mmdb::Residue *residue_p) const {

   std::vector<std::string> v;
   v.reserve(6);
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string atom_name(at->name);
      if (base_atom_name_set.find(atom_name) != base_atom_name_set.end()) {
	 v.push_back(atom_name);
      }
   }
   return v;
}

std::pair<bool,clipper::Coord_orth>
coot::stack_and_pair::get_base_centre(mmdb::Residue *residue_p) const {

   std::pair<bool,clipper::Coord_orth> p(false, clipper::Coord_orth(0,0,0));

   unsigned int n_centres = 0;
   clipper::Coord_orth centre_sum(0.0, 0.0, 0.0);
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      std::string atom_name(at->name);
      if (base_atom_name_set.find(atom_name) != base_atom_name_set.end()) {
	 centre_sum += co(at);
	 n_centres++;
      }
   }

   if (n_centres > 3) {
      p.first = true;
      double n_d = static_cast<double>(n_centres);
      double r_n_d = 1.0/n_d;
      p.second = clipper::Coord_orth(centre_sum.x() * r_n_d,
				     centre_sum.y() * r_n_d,
				     centre_sum.z() * r_n_d);
   }
   return p;
}

std::map<mmdb::Residue *, clipper::Coord_orth>
coot::stack_and_pair::calculate_residue_normals(const std::vector<std::pair<bool,mmdb::Residue *> > &residues_vec) const {

   std::map<mmdb::Residue *, clipper::Coord_orth> m;
   for (std::size_t i=0; i<residues_vec.size(); i++) {
      mmdb::Residue *r = residues_vec[i].second;
      if (r) {
	 std::pair<bool, clipper::Coord_orth> normal = get_base_normal(r);
	 if (normal.first) {
	    m[r] = normal.second;
	 }
      }
   }
   return m;
}

std::map<mmdb::Residue *, clipper::Coord_orth>
coot::stack_and_pair::calculate_residue_normals(mmdb::Atom **SelAtom, int n_sel_atoms) {

   std::set<mmdb::Residue *> done_res;
   std::map<mmdb::Residue *, clipper::Coord_orth> m;
   for (int i=0; i<n_sel_atoms; i++) {
      mmdb::Residue *r = SelAtom[i]->residue;
      if (done_res.find(r) == done_res.end()) {
	 std::pair<bool, clipper::Coord_orth> bn = get_base_normal(r);
	 if (bn.first) {
	    m[r] = bn.second;
	 }
	 done_res.insert(r);
      }
   }
   return m;
}

int
coot::stack_and_pair::mark_donors_and_acceptors(mmdb::Manager *mol, int selection_handle,
						const protein_geometry &geom) {

   mmdb::PAtom *sel_atoms = 0;
   int n_sel_atoms;
   mol->GetSelIndex(selection_handle, sel_atoms, n_sel_atoms);

   int udd_h_bond_type_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "hb_type");
   for (int i=0; i<n_sel_atoms; i++) {
      mmdb::Atom *at = sel_atoms[i];
      std::string name = at->name;
      std::string res_name = at->GetResName();
      int h_bond_type = geom.get_h_bond_type(name, res_name, protein_geometry::IMOL_ENC_ANY);
      at->PutUDData(udd_h_bond_type_handle, h_bond_type);
   }

   return udd_h_bond_type_handle;

}

bool
coot::stack_and_pair::similar_normals(mmdb::Residue *res_1, mmdb::Residue *res_2,
				      const std::map<mmdb::Residue *, clipper::Coord_orth> &normal_map) const {

   // normals are presumed to be normalized.

   bool status = false;
   std::map<mmdb::Residue *, clipper::Coord_orth>::const_iterator it_1;
   std::map<mmdb::Residue *, clipper::Coord_orth>::const_iterator it_2;

   it_1 = normal_map.find(res_1);
   it_2 = normal_map.find(res_2);

   if (it_1 != normal_map.end()) {
      if (it_2 != normal_map.end()) {
	 const clipper::Coord_orth &n1 = it_1->second;
	 const clipper::Coord_orth &n2 = it_2->second;
	 double dp = clipper::Coord_orth::dot(n1, n2);
	 double cos_angle_crit = cos(angle_crit);
	 if ((dp > cos_angle_crit) || (dp < -cos_angle_crit)) {
	    status = true;
	 }
      }
   }
   return status;
}


std::vector<coot::stack_and_pair::paired_residues_info_t>
coot::stack_and_pair::paired_residues(mmdb::Manager *mol,
				      const std::vector<std::pair<bool, mmdb::Residue *> > &residues_vec,
				      bool residues_are_all_moving_flag,
				      const coot::protein_geometry &geom) {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   std::vector<coot::stack_and_pair::paired_residues_info_t> v;
   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > atom_vec;

   float dist_crit = 3.8; // less than this to form a pairing bond
   float dist_crit_sqrt = dist_crit * dist_crit;

   int selection_handle_moving;
   int n_selected_atoms_moving = 0;
   int selection_handle_all = mol->NewSelection(); // d
   mol->SelectAtoms(selection_handle_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");

   mmdb::Atom **selected_atoms_all = 0;
   mmdb::Atom **selected_atoms_moving = 0;
   int n_selected_atoms_all;
   mol->GetSelIndex(selection_handle_all, selected_atoms_all, n_selected_atoms_all);

   auto tp_1 = std::chrono::high_resolution_clock::now();

   if (residues_are_all_moving_flag) {
      selection_handle_moving = selection_handle_all;
      n_selected_atoms_moving = n_selected_atoms_all;
      selected_atoms_moving = selected_atoms_all;
   } else {
      auto tp_2 = std::chrono::high_resolution_clock::now();
      selection_handle_moving = mol->NewSelection(); // d
      for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	 if (residues_vec[ires].first == false) {
	    mol->SelectAtoms(selection_handle_moving, 1,
			     residues_vec[ires].second->GetChainID(),
			     residues_vec[ires].second->GetSeqNum(),
			     residues_vec[ires].second->GetInsCode(),
			     residues_vec[ires].second->GetSeqNum(),
			     residues_vec[ires].second->GetInsCode(),
			     "*", // any residue name
			     "*", // atom name
			     "*", // elements
			     "*",  // alt loc.
			     mmdb::SKEY_OR
			     );
	 }
      }
      auto tp_3 = std::chrono::high_resolution_clock::now();
      auto d32 = std::chrono::duration_cast<std::chrono::microseconds>(tp_3 - tp_2).count();
      std::cout << "------------------ timings: for residues SelectAtoms() "
		<< d32 << " ms" << std::endl;
      mol->GetSelIndex(selection_handle_moving, selected_atoms_moving, n_selected_atoms_moving);
   }
   
   auto tp_4 = std::chrono::high_resolution_clock::now();
   auto tp_5 = std::chrono::high_resolution_clock::now();
   if (n_selected_atoms_moving > 1) {

      // Does this selection contain nucleic acid?
      if (contains_nucleic_acid(selected_atoms_moving, n_selected_atoms_moving)) {
	 mmdb::mat44 my_matt;
	 for (int i=0; i<4; i++) 
	    for (int j=0; j<4; j++) 
	       my_matt[i][j] = 0.0;      
	 for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   
	 mmdb::Contact *pscontact = NULL;
	 int n_contacts;
	 long i_contact_group = 1;

	 tp_4 = std::chrono::high_resolution_clock::now();
	 mol->SeekContacts(selected_atoms_moving, n_selected_atoms_moving,
			   selected_atoms_all, n_selected_atoms_all,
			   0.01, dist_crit,
			   0, // seqDist 0 -> also in same res.
			   pscontact, n_contacts,
			   0, &my_matt, i_contact_group);

	 if (false)
	    std::cout << "pairing:: found n_contacts in atom selection: " << n_contacts << std::endl;

	 tp_5 = std::chrono::high_resolution_clock::now();
	 if (n_contacts > 0) {
	    if (pscontact) {

	       // interesting Hydrogen bonds - but not base pairing
	       std::set<std::string> excluded_oxygens;
	       excluded_oxygens.insert(" OP1"); excluded_oxygens.insert(" OP2"); excluded_oxygens.insert(" O2'");
	       excluded_oxygens.insert(" O3'"); excluded_oxygens.insert(" O5'");

	       int hb_type_udd_handle = mark_donors_and_acceptors(mol, selection_handle_all, geom); // using UDD data

	       for (int i_contact=0; i_contact<n_contacts; i_contact++) {
		  mmdb::Atom *at_1 = selected_atoms_moving[pscontact[i_contact].id1];
		  mmdb::Atom *at_2 = selected_atoms_all[pscontact[i_contact].id2];

		  // in different residues and both are are O or N, and are close
		  // enough together and both are nucleic acids

		  if (at_1->residue != at_2->residue) {
		     std::string ele_1(at_1->element);
		     std::string ele_2(at_2->element);

		     // PDBv3 FIXME
		     if (ele_1 == " O" || ele_1 == " N") {
			if (ele_2 == " O" || ele_2 == " N") {

			   float dx = at_1->x - at_2->x;
			   float dy = at_1->y - at_2->y;
			   float dz = at_1->z - at_2->z;
			   float dd = dx * dx + dy * dy + dz * dz;
			   if (dd < dist_crit_sqrt) {
			      int hb_type_1 = coot::HB_UNASSIGNED;
			      int hb_type_2 = coot::HB_UNASSIGNED;

			      at_1->GetUDData(hb_type_udd_handle, hb_type_1);
			      at_2->GetUDData(hb_type_udd_handle, hb_type_2);

			      if (hb_type_1 == coot::HB_ACCEPTOR || hb_type_1 == coot::HB_BOTH) {
				 if (hb_type_2 == coot::HB_DONOR || hb_type_2 == coot::HB_BOTH) {

				    if (at_1->GetChain() == at_2->GetChain()) {
				       int residue_index_1 = at_1->residue->index;
				       int residue_index_2 = at_2->residue->index;
				       int residue_index_delta = residue_index_2 - residue_index_1;
				       if (abs(residue_index_delta) < 2)
					  continue;
				    }

				    if (util::is_nucleotide(at_1->residue)) {
				       if (util::is_nucleotide(at_2->residue)) {

					  if (similar_normals(at_1->residue, at_2->residue, normal_map)) {

					     // also, to stop base pairing the residue above
					     // or below, we need to check the dot product
					     // of the base and the line between the H-bonded atoms.
					     // They should be almost pependicular - and won't be
					     // when this base is paired with the (wrong) base
					     // above or below on the other side.

					     clipper::Coord_orth pt_1 = co(at_1);
					     clipper::Coord_orth pt_2 = co(at_2);
					     clipper::Coord_orth atom_atom_unit_vector((pt_2 - pt_1).unit());

					     double dp_1 = clipper::Coord_orth::dot(atom_atom_unit_vector, normal_map[at_1->residue]);
					     double dp_2 = clipper::Coord_orth::dot(atom_atom_unit_vector, normal_map[at_2->residue]);

					     if (false) {
						std::cout << " dot product 1 " << dp_1 << " " << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << "\n";
						std::cout << " dot product 2 " << dp_2 << " " << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << "\n";
					     }

					     if (std::abs(dp_1) < 0.5) {
						if (std::abs(dp_2) < 0.5) {

						   // no ribose or phosphate atoms:
						   std::string name_1(at_1->name);
						   std::string name_2(at_2->name);
						   if (excluded_oxygens.find(name_1) == excluded_oxygens.end()) {
						      if (excluded_oxygens.find(name_2) == excluded_oxygens.end()) {

							 std::pair<mmdb::Atom *, mmdb::Atom *> p(at_1, at_2);
							 atom_vec.push_back(p);
						      }
						   }
						}
					     }
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   auto tp_6 = std::chrono::high_resolution_clock::now();


   // If there are 2 or more bonds to any residue from the same (other) residue, then we should add
   // long distance restraints (I think)
   //
   // Another time.

   // now convert atom_vec into a vector of paired_residues_info_t
   for (std::size_t i=0; i<atom_vec.size(); i++) {
      mmdb::Atom *at_1 = atom_vec[i].first;
      mmdb::Atom *at_2 = atom_vec[i].second;
      mmdb::Residue *res_1 = at_1->residue;
      mmdb::Residue *res_2 = at_2->residue;
      bool found = false;

      for (std::size_t j=0; j<v.size(); j++) {
	 if (v[j].res_1 == res_1) {
	    if (v[j].res_2 == res_2) {
	       v[j].atom_pair_vec.push_back(atom_vec[i]);
	       found = true;
	       break;
	    }
	 }
	 if (v[j].res_1 == res_2) {
	    if (v[j].res_2 == res_1) {
	       std::swap(atom_vec[i].first, atom_vec[i].second);
	       v[j].atom_pair_vec.push_back(atom_vec[i]);
	       found = true;
	       break;
	    }
	 }
      }
      if (! found) {
	 std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > apv;
	 apv.push_back(std::pair<mmdb::Atom *, mmdb::Atom *>(at_1, at_2));
	 paired_residues_info_t pri(res_1, res_2, apv);
	 v.push_back(pri);
      }
   }
   auto tp_7 = std::chrono::high_resolution_clock::now();

   mol->DeleteSelection(selection_handle_all);
   mol->DeleteSelection(selection_handle_moving);

   auto tp_8 = std::chrono::high_resolution_clock::now();

   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
   auto d65 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_6 - tp_5).count();
   auto d76 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_7 - tp_6).count();
   auto d87 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_8 - tp_7).count();
   std::cout << "------------------ timings: for paired_residues(): "
	     << d10 << " " << d54 << " " << d65 << " "
	     << d76 << "  " << d87 << " " << std::endl;
   return v;
}

std::vector<coot::stack_and_pair::stacked_planes_info_t>
coot::stack_and_pair::stacked_residues(mmdb::Manager *mol) {

   // Dealing with alt-confs here (alt-confs are used in parallel plane restraitns)
   // is messy. Add it another time.

   std::vector<stacked_planes_info_t> spc;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 mmdb::Chain *chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<(nres-1); ires++) {
	    mmdb::Residue *residue_this = chain_p->GetResidue(ires);
	    mmdb::Residue *residue_next = chain_p->GetResidue(ires+1);
	    if (residue_this) {
	       if (residue_next) {
		  std::pair<bool,clipper::Coord_orth> base_centre_1 = get_base_centre(residue_this);
		  std::pair<bool,clipper::Coord_orth> base_centre_2 = get_base_centre(residue_next);
		  if (base_centre_1.first) {
		     if (base_centre_2.first) {
			double dd = (base_centre_2.second-base_centre_1.second).lengthsq();
			double d = sqrt(dd);
			if (d < 5.0) { // generous?
			   if (similar_normals(residue_this, residue_next, normal_map)) {
			      std::vector<std::string> an1 = get_base_atom_names(residue_this);
			      std::vector<std::string> an2 = get_base_atom_names(residue_next);
			      stacked_planes_info_t sp(residue_this, residue_next, an1, an2);
			      spc.push_back(sp);
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return spc;
}


