
#include "coot_molecule.hh"


void
coot::molecule_t::generate_self_restraints(float local_dist_max, const coot::protein_geometry &geom) {

   // Find all the contacts in chain_id that are less than or equal to local_dist_max
   // that are not bonded or related by an angle.

   int selHnd = atom_sel.mol->NewSelection(); // deleted in generate_local_self_restraints()
                                              // Confusing! bleugh! FIXME

   atom_sel.mol->SelectAtoms(selHnd, 0, "*",
			     mmdb::ANY_RES, "*", // start, insertion code
			     mmdb::ANY_RES, "*", // end, insertion code
			     "*", // residue name
			     "*",
			     "*", // elements
			     "*"); // alt locs

   generate_local_self_restraints(selHnd, local_dist_max, geom);
}

void
coot::molecule_t::generate_chain_self_restraints(float local_dist_max,
                                                 const std::string &chain_id,
                                                 const coot::protein_geometry &geom) {

   // Find all the contacts in chain_id that are less than or equal to local_dist_max
   // that are not bonded or related by an angle.

   int selHnd = atom_sel.mol->NewSelection(); // - check the deletion

   atom_sel.mol->SelectAtoms(selHnd, 0, chain_id.c_str(),
			     mmdb::ANY_RES, "*", // start, insertion code
			     mmdb::ANY_RES, "*", // end, insertion code
			     "*", // residue name
			     "*",
			     "*", // elements
			     "*"); // alt locs

   generate_local_self_restraints(selHnd, local_dist_max, geom);

   // atom_sel.mol->DeleteSelection(selHnd);
}

void
coot::molecule_t::generate_local_self_restraints(float local_dist_max,
                                                 const std::vector<coot::residue_spec_t> &residue_specs,
                                                 const coot::protein_geometry &geom) {

   // Find all the contacts in chain_id that are less than or equal to local_dist_max
   // that are not bonded or related by an angle.

   int selHnd = coot::specs_to_atom_selection(residue_specs, atom_sel.mol, 0);
   if (selHnd >= 0) {
      generate_local_self_restraints(selHnd, local_dist_max, geom);
   }
}

void
coot::molecule_t::generate_local_self_restraints(int selHnd, float local_dist_max,
                                                 const coot::protein_geometry &geom) {

   // clear what's already there - if anything
   extra_restraints.bond_restraints.clear();

   int nSelAtoms;
   mmdb::PPAtom SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, SelAtom, nSelAtoms);

   // bonded_neighbours in this case, means bonded or angle-related
   // bonded_neighbours["ALA"] -> all bond pairs and 1-3 angles
   std::map<std::string, std::vector<std::pair<std::string, std::string> > > bonded_neighbours;

   // now find contacts:
   //
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
	 my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   //
   atom_sel.mol->SeekContacts(SelAtom,   nSelAtoms,
			      SelAtom,   nSelAtoms,
			      0.1, local_dist_max, // min, max distances
			      0,        // seqDist 0 -> in same res also
			      pscontact, n_contacts,
			      0, &my_matt, i_contact_group);

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i=0; i<n_contacts; i++) {

            // 20221223-PE don't go both ways:
            if (pscontact[i].id1 > pscontact[i].id2) continue;

	    mmdb::Atom *at_1 = SelAtom[pscontact[i].id1];
	    mmdb::Atom *at_2 = SelAtom[pscontact[i].id2];
	    std::string ele_1 = at_1->element;
	    std::string ele_2 = at_2->element;
	    if (ele_1 != " H" && ele_2 != " H") {
	       bool ignore_this = false; // set for bonded and angled atoms
	       bool in_same_res = false;
	       if (at_1->residue == at_2->residue)
		  in_same_res = true;

	       if (in_same_res) {
		  std::string comp_id = at_1->residue->GetResName();
		  std::string at_name_1 = at_1->GetAtomName();
		  std::string at_name_2 = at_2->GetAtomName();

		  // This is slow
		  // if (geom.are_bonded_or_angle_related(comp_id, at_name_1, at_name_2))
		  // ignore_this = true;
		  //

		  std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator it;
		  it = bonded_neighbours.find(comp_id);
		  std::vector<std::pair<std::string, std::string> > bps;
		  if (it == bonded_neighbours.end()) {
		     bps = geom.get_bonded_and_1_3_angles(comp_id, imol_no);
		     bonded_neighbours[comp_id] = bps;
		  } else {
		     bps = it->second;
		  }

		  for (unsigned int ipr=0; ipr<bps.size(); ipr++) {
		     if (at_name_1 == bps[ipr].first) {
			if (at_name_2 == bps[ipr].second) {
			   ignore_this = true;
			   break;
			}
		     }
		     if (at_name_2 == bps[ipr].first) {
			if (at_name_1 == bps[ipr].second) {
			   ignore_this = true;
			   break;
			}
		     }
		  }
	       } else {
		  std::string at_name_1 = at_1->GetAtomName();
		  std::string at_name_2 = at_2->GetAtomName();
		  // by hand (hmmm)
		  // PDBv3 FIXME
		  if ((at_name_1 == " N  " && at_name_2 == " C  ") ||
		      (at_name_1 == " C  " && at_name_2 == " N  ") ||
		      (at_name_1 == " N  " && at_name_2 == " O  ") ||
		      (at_name_1 == " N  " && at_name_2 == " CA ") ||
		      (at_name_1 == " O  " && at_name_2 == " N  ") ||
		      (at_name_1 == " CA " && at_name_2 == " N  ")) {
		     ignore_this = true;
		  }
	       }

	       if (! ignore_this) {
		  clipper::Coord_orth p1 = coot::co(at_1);
		  clipper::Coord_orth p2 = coot::co(at_2);
		  double dist = sqrt((p1-p2).lengthsq());
		  double esd  = 0.05;
		  coot::atom_spec_t atom_spec_1(at_1);
		  coot::atom_spec_t atom_spec_2(at_2);
		  int idx_1 = -1;
		  int idx_2 = -1;
		  at_1->GetUDData(atom_sel.UDDAtomIndexHandle, idx_1);
		  at_2->GetUDData(atom_sel.UDDAtomIndexHandle, idx_2);
		  atom_spec_1.int_user_data = idx_1;
		  atom_spec_2.int_user_data = idx_2;
		  coot::extra_restraints_t::extra_geman_mcclure_restraint_t gmr(atom_spec_1, atom_spec_2,
										dist, esd);

		  // 20191120-PE self restraints are GM restraints, not bond
		  // restraints (previously bond restraints were GM only)
		  //
		  // extra_restraints.bond_restraints.push_back(br);

		  extra_restraints.geman_mcclure_restraints.push_back(gmr);

	       }
	    }
	 }
      }
   }

   // 20180510-PE surely I want to update the representation even if there are no
   //             extra bond restraints?
   //
   // if (extra_restraints.bond_restraints.size())
   // update_extra_restraints_representation();
   //
   // update_extra_restraints_representation(); libcootapi doesn't care about this

   atom_sel.mol->DeleteSelection(selHnd);
}

void
coot::molecule_t::add_parallel_plane_restraint(coot::residue_spec_t spec_1,
                                               coot::residue_spec_t spec_2) {

   std::vector<std::string> ap_1_names;
   std::vector<std::string> ap_2_names;
   std::string alt_conf_1; // only a nutter would want to add parallel plane restraints
   std::string alt_conf_2; // to a residue with an alt conf.

   mmdb::Residue *r_1 = get_residue(spec_1);
   mmdb::Residue *r_2 = get_residue(spec_2);

   if (r_1) {
      if (r_2) {

	 std::string rn_1 = r_1->GetResName();
	 std::string rn_2 = r_2->GetResName();

	 ap_1_names = nucelotide_residue_name_to_base_atom_names(rn_1);
	 ap_2_names = nucelotide_residue_name_to_base_atom_names(rn_2);

	 if (ap_1_names.empty()) ap_1_names = residue_name_to_plane_atom_names(rn_1);
	 if (ap_2_names.empty()) ap_2_names = residue_name_to_plane_atom_names(rn_2);

	 std::cout << "ap_2_names ";
	 for (auto i: ap_2_names)
	    std::cout << i << " ";
	 std::cout << "" << std::endl;

	 std::cout << "Adding parallel plane restraint " << spec_1 << " " << spec_2 << std::endl;
	 coot::parallel_planes_t pp(spec_1, spec_2, ap_1_names, ap_2_names,
				    alt_conf_1, alt_conf_2);

	 extra_restraints.parallel_plane_restraints.push_back(pp);

      } else {
	 std::cout << "INFO:: missing residue 2 " << spec_2 << std::endl;
      }
   } else {
	 std::cout << "INFO:: missing residue 1 " << spec_1 << std::endl;
   }

}



// which uses:
std::vector<std::string>
coot::molecule_t::nucelotide_residue_name_to_base_atom_names(const std::string &rn) const {

   std::vector<std::string> v;

   return v;
      
}


// for non-bases, normal amino acids (simple-minded, currently).
std::vector<std::string>
coot::molecule_t::residue_name_to_plane_atom_names(const std::string &rn) const {

   std::vector<std::string> v;

   return v;
}


void
coot::molecule_t::clear_extra_restraints() {
   extra_restraints.clear();
}
