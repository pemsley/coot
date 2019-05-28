
#include "simple-restraint.hh"
#include "coot-utils/contacts-by-bricks.hh"

int
coot::restraints_container_t::make_restraints_ng(int imol,
						 const coot::protein_geometry &geom,
						 coot::restraint_usage_Flags flags_in,
						 bool do_residue_internal_torsions,
						 bool do_trans_peptide_restraints,
						 float rama_plot_target_weight,
						 bool do_rama_plot_restraints,
						 bool do_auto_helix_restraints,
						 bool do_auto_strand_restraints,
						 coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds,
						 bool do_link_restraints,
						 bool do_flank_restraints) {

   std::cout << "-------------------- in make_restraints_ng() mol: " << mol << std::endl;
   debug_sets();

   auto tp_0 = std::chrono::high_resolution_clock::now();
   restraints_usage_flag = flags_in;
   if (n_atoms) {
      mark_OXT(geom);
      make_monomer_restraints(imol, geom, do_residue_internal_torsions);
      auto tp_1 = std::chrono::high_resolution_clock::now();

      std::map<mmdb::Residue *, unsigned int> residue_link_count_map =
	 make_link_restraints_ng(geom, do_rama_plot_restraints, do_trans_peptide_restraints);

      auto tp_2 = std::chrono::high_resolution_clock::now();
      // now, what is linked to what?
      // that is in bonded_atom_indices. That is a vector of what is bonded or angle-bonded
      // to each atom.

      // don't forget about rama to flanking residues here:

      make_flanking_atoms_restraints_ng(residue_link_count_map);

      auto tp_3 = std::chrono::high_resolution_clock::now();
      reduced_angle_info_container_t raic(restraints_vec);

      auto tp_4 = std::chrono::high_resolution_clock::now();
      make_non_bonded_contact_restraints_ng(imol, raic, geom);
      auto tp_5 = std::chrono::high_resolution_clock::now();

      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      std::cout << "------------------- timing for make_restraints_ng(): monomers: "
		<< d10 << " links: " << d21 << " " << " milliseconds " << std::endl;

   }
   make_df_restraints_indices();
   make_distortion_electron_density_ranges();

   // info();  - are the NBCs correct?

   return restraints_vec.size();
}

void
coot::restraints_container_t::make_flanking_atoms_restraints_ng(const std::map<mmdb::Residue *, unsigned int>
 &residue_link_count_map) const {

   // residue 20 is being refined, we need to make link restraints to residues 19 and 21 using
   // fixed_neighbours_set
   //
   // residues_vec contains residue 20
   //
   // fixed_neighbours_set are the (fixed) neighbours of residue 20.
   // 
   // mol contains residue 20 and the residues around it (made by the function that
   // creates the restraints object)

   if (true) {
      std::cout << "residue_link_count_map size " << residue_link_count_map.size() << std::endl;
      std::map<mmdb::Residue *, unsigned int>::const_iterator it;
      for (it=residue_link_count_map.begin();
	   it!=residue_link_count_map.end(); it++) {
	 std::cout << "        residue_link_count_map: "
		   << residue_spec_t(it->first) << " " << it->second << std::endl;
      }
   }

   std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
   for (it=fixed_neighbours_set.begin(); it!=fixed_neighbours_set.end(); it++) {
      mmdb::Residue *residue_p = it->first;
      std::cout << "Considering residue " << residue_spec_t(residue_p) << " for link restraints"
		<< std::endl;

      unsigned int n_link_for_residue = 0;
      std::map<mmdb::Residue *, unsigned int>::const_iterator itm = residue_link_count_map.find(residue_p);
      if (itm != residue_link_count_map.end()) {
	 n_link_for_residue = itm->second;
      } else {
	 if (true)
	    std::cout << "OOps " << residue_spec_t(residue_p) << " not found in link map" << std::endl;
      }

      // was it polymer linked at both ends?
      //
      if (n_link_for_residue != 2) {

	 int index_for_residue = residue_p->index;
	 const std::set<mmdb::Residue *> &s = it->second;
	 std::set<mmdb::Residue *>::const_iterator its;
	 for (its=s.begin(); its!=s.end(); its++) {
	    mmdb::Residue *neighb = *its;
	    int index_for_neighb = neighb->index;
	    int index_delta = index_for_neighb - index_for_residue;
	    if (index_delta == -1 || index_delta == 1) {
	       std::cout << "     fixed neigb: " << residue_spec_t(*its) << std::endl;
	    }
	 }
      }
   }
   
   
}


// Use raic to know what is 1-4 related
void
coot::restraints_container_t::make_non_bonded_contact_restraints_ng(int imol,
								    const coot::restraints_container_t::reduced_angle_info_container_t &raic,
								    const coot::protein_geometry &geom) {

   // potentially multithreadable.

   // raic can be a member of the class
   // to make it: reduced_angle_info_container_t raic(restraints_vec);

   float dist_max = 8.0;

   // cache the energy types:
   std::vector<std::string> energy_type_for_atom(n_atoms);

   // needs timing test - might be slow
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atom[i];
      energy_type_for_atom[i] = get_type_energy(imol, at, geom);
   }

   std::map<std::string, std::pair<bool, std::vector<std::list<std::string> > > > residue_ring_map_cache;

   std::set<unsigned int> fixed_atom_flags_set; // signed to unsigned conversion - bleugh.
   std::set<int>::const_iterator it;
   for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
      fixed_atom_flags_set.insert(*it);

   contacts_by_bricks cb(atom, n_atoms, fixed_atom_flags_set);
   cb.set_dist_max(dist_max);
   std::vector<std::set<unsigned int> > vcontacts;
   cb.find_the_contacts(&vcontacts);

   for (std::size_t i=0; i<vcontacts.size(); i++) {
      const std::set<unsigned int> &n_set = vcontacts[i];
      mmdb::Atom *at_1 = atom[i];
      std::set<unsigned int>::const_iterator it;
      for (it=n_set.begin(); it!=n_set.end(); it++) {

	 // bonded_atom_indices should be a vector of *sets*!
	 // use bonded_atom_indices[i].find(*it)
	 if (std::find(bonded_atom_indices[i].begin(),
		       bonded_atom_indices[i].end(), *it) != bonded_atom_indices[i].end()) {
	    // std::cout << "found a bonded atom pair " << i << " " << *it << std::endl;
	    continue;
	 }
	 
	 mmdb::Atom *at_2 = atom[*it];

	 if (fixed_atom_indices.find(i) != fixed_atom_indices.end())
	    if (fixed_atom_indices.find(*it) != fixed_atom_indices.end())
	       continue;

	 if (*it < i) /* only add NBC one way round */
	    continue;

	 std::string res_name_1 = at_1->GetResName();
	 std::string res_name_2 = at_2->GetResName();
	 int res_no_1 = at_1->GetSeqNum();
	 int res_no_2 = at_2->GetSeqNum();

	 bool in_same_residue = (at_1->residue == at_2->residue);

	 // Hacketty-hack!
	 bool add_it = true;
	 if (res_name_1 == "PRO" || res_name_1 == "HYP") {
	    int res_no_pro   = res_no_1;
	    int res_no_other = res_no_2;
	    if (res_no_pro == (res_no_other + 1)) {
	       std::string atom_name = at_1->name;
	       if (atom_name == " CD ")  // PDBv3 FIXME
		  add_it = false;
	    }
	 }
	 if (res_name_1 == "ASN" || res_name_2 == "NAG") {
	    std::string atom_name_1(at_1->name);
	    std::string atom_name_2(at_2->name);
	    if (atom_name_1 == " OD1")
	       if (atom_name_2 == " C1 ")
		  add_it = false;
	    if (atom_name_1 == "HD21")
	       if (atom_name_2 == " C1 ")
		  add_it = false;
	 }
	 if (! add_it) continue;

	 bool in_same_ring_flag = false;
	 std::string atom_name_1(at_1->name);
	 std::string atom_name_2(at_2->name);
	 if (in_same_residue)
	    in_same_ring_flag = is_in_same_ring(imol, at_2->residue,
						residue_ring_map_cache, // add to this
						atom_name_1, atom_name_2, geom);
	 bool is_1_4_related = raic.is_1_4(i, *it);
	 bool is_O_C_1_5_related = check_for_O_C_1_5_relation(at_1, at_2);

	 double dist_min = 3.4;

	 if (is_1_4_related) {
	    dist_min = 2.64; // was 2.7 but c.f. guanine ring distances
	 }

	 if (is_hydrogen(at_1)) dist_min -= 0.7;
	 if (is_hydrogen(at_2)) dist_min -= 0.7;

	 if (is_O_C_1_5_related)
	    dist_min = 2.84;


	 bool strange_exception = false;
	 int rn_diff = abs(res_no_2 - res_no_1);
	 std::vector<bool> fixed_atom_flags = make_fixed_flags(i, *it);
	 if (rn_diff == 1) { // or 0 for insertion codes?
	    std::string atom_name_1 = at_1->GetAtomName();
	    std::string atom_name_2 = at_2->GetAtomName();
	    if (fixed_atom_flags.size()) {
	       if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
		  if (atom_name_1 == " O  ")
		     if (atom_name_2 == " CA ")
			strange_exception = true;
		  if (atom_name_1 == " CA ")
		     if (atom_name_2 == " O  ")
			strange_exception = true;
		  if (atom_name_1 == " N  ")
		     if (atom_name_2 == " CB ")
			strange_exception = true;
		  if (atom_name_1 == " CB ")
		     if (atom_name_2 == " N  ")
			strange_exception = true;
		  if (atom_name_1 == " C  ")
		     if (atom_name_2 == " CB ")
			strange_exception = true;
	       }
	    }
	 }
	 if (rn_diff == 2) {
	    if (fixed_atom_flags.size()) {
	       if (fixed_atom_flags[0] || fixed_atom_flags[1]) {
		  std::string atom_name_1 = at_1->GetAtomName();
		  std::string atom_name_2 = at_2->GetAtomName();
		  if (atom_name_1 == " C  ")
		     if (atom_name_2 == " N  ")
			strange_exception = true;
		  if (atom_name_1 == " N  ")
		     if (atom_name_2 == " C  ")
			strange_exception = true; // 3.1 would be enough

	       }
	    }
	 }

	 if (strange_exception)
	    dist_min = 2.7;

	 bool is_H_non_bonded_contact = false;

	 if (is_hydrogen(at_1)) {
	    is_H_non_bonded_contact = true;
	    if (H_parent_atom_is_donor(at_1))
	       if (is_acceptor(energy_type_for_atom[*it], geom))
		  dist_min -= 0.7;
	 }
	 if (is_hydrogen(at_2)) {
	    is_H_non_bonded_contact = true;
	    if (H_parent_atom_is_donor(at_2))
	       if (is_acceptor(energy_type_for_atom[i], geom))
		  dist_min -= 0.7;
	 }

	 simple_restraint::nbc_function_t nbcf = simple_restraint::LENNARD_JONES;
	 simple_restraint r(NON_BONDED_CONTACT_RESTRAINT,
			    nbcf, i, *it,
			    energy_type_for_atom[i],
			    energy_type_for_atom[*it],
			    is_H_non_bonded_contact,
			    fixed_atom_flags, dist_min);
	 restraints_vec.push_back(r);

      }
   }
}


coot::restraints_container_t::link_restraints_counts
coot::restraints_container_t::make_link_restraints_for_link_ng(const std::string &link_type,
							       mmdb::Residue *res_1,
							       mmdb::Residue *res_2,
							       bool is_fixed_first_residue,
							       bool is_fixed_second_residue,
							       bool do_trans_peptide_restraints,
							       const coot::protein_geometry &geom) {

   link_restraints_counts lrc;

   if (restraints_usage_flag & BONDS_MASK)
      lrc.n_link_bond_restr += add_link_bond(link_type,
					     res_1, res_2,
					     is_fixed_first_residue,
					     is_fixed_second_residue,
					     geom);

   if (restraints_usage_flag & ANGLES_MASK)
      lrc.n_link_angle_restr += add_link_angle(link_type,
					   res_1, res_2,
					   is_fixed_first_residue,
					   is_fixed_second_residue,
					   geom);

   if (restraints_usage_flag & TRANS_PEPTIDE_MASK)
      if (do_trans_peptide_restraints)
	 lrc.n_link_trans_peptide += add_link_trans_peptide(res_1, res_2,
							is_fixed_first_residue,
							is_fixed_second_residue,
							geom);

   if (restraints_usage_flag & PLANES_MASK)
      lrc.n_link_plane_restr += add_link_plane(link_type,
					   res_1, res_2,
					   is_fixed_first_residue,
					   is_fixed_second_residue,
					   geom);
   return lrc;
}

std::string
coot::restraints_container_t::find_peptide_link_type_ng(mmdb::Residue *res_1,
							mmdb::Residue *res_2,
							const coot::protein_geometry &geom) const {

   std::string link_type;

   std::string t1;
   std::string t2;

   std::string residue_type_1 = res_1->name;
   std::string residue_type_2 = res_2->name;

   for (unsigned int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_1)) {
	 t1 = geom[idr].second.residue_info.group;
	 break;
      }
   }
   for (unsigned int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_2)) {
	 t2 = geom[idr].second.residue_info.group;
	 break;
      }
   }

   // This should be true because we do a test on amino acid type at the start of try_make_peptide_link_ng
   //
   if (t1 == "L-peptide" || t1 == "D-peptide" || t1 == "M-peptide" || t1 == "P-peptide" || t1 == "peptide") {
      if (t2 == "L-peptide" || t2 == "D-peptide" || t2 == "M-peptide" || t2 == "P-peptide" || t2 == "peptide") {
	 if (residue_type_2 == "PRO" || residue_type_2 == "HYP") {
	    link_type = "PTRANS";
	 } else {
	    link_type = "TRANS";
	 }
      }
   }

   // CIS and PCIS will need a torsion check around
   // r1::CA - r1::C - r2::N - r2::CA

   return link_type;
}



std::pair<bool, coot::restraints_container_t::link_restraints_counts>
coot::restraints_container_t::try_make_peptide_link_ng(const coot::protein_geometry &geom,
						       std::pair<bool, mmdb::Residue *> res_1_pair,
						       std::pair<bool, mmdb::Residue *> res_2_pair,
						       bool do_rama_plot_restraints,
						       bool do_trans_peptide_restraints) {

   mmdb::Residue *res_1 = res_1_pair.second;
   mmdb::Residue *res_2 = res_2_pair.second;
   std::string res_name_1(res_1->GetResName());
   std::string res_name_2(res_2->GetResName());
   bool status = false;
   link_restraints_counts lrc;

   if (util::is_standard_amino_acid_name(res_name_1)) {
      if (util::is_standard_amino_acid_name(res_name_2)) {
	 mmdb::Atom **residue_1_atoms = 0;
	 mmdb::Atom **residue_2_atoms = 0;
	 int n_residue_1_atoms;
	 int n_residue_2_atoms;
	 res_1->GetAtomTable(residue_1_atoms, n_residue_1_atoms);
	 res_2->GetAtomTable(residue_2_atoms, n_residue_2_atoms);
	 for (int iat_1=0; iat_1<n_residue_1_atoms; iat_1++) {
	    mmdb::Atom *at_1 = residue_1_atoms[iat_1];
	    std::string at_name_1(at_1->GetAtomName());
	    if (at_name_1 == " C  ") { // PDBv3 FIXE
	       std::string alt_conf_1(at_1->altLoc);
	       for (int iat_2=0; iat_2<n_residue_2_atoms; iat_2++) {
		  mmdb::Atom *at_2 = residue_2_atoms[iat_2];
		  std::string at_name_2(at_2->GetAtomName());
		  if (at_name_2 == " N  ") { // PDBv3 FIXE
		     std::string alt_conf_2(at_2->altLoc);
		     if (alt_conf_1 == alt_conf_2 || alt_conf_1.empty() || alt_conf_2.empty()) {

			// find_peptide_link_type_ng doesn't test the geometry -
			// just the residues types
			std::string link_type = find_peptide_link_type_ng(res_1, res_2, geom);
			if (! link_type.empty()) {
			   bool is_fixed_first_residue  = res_1_pair.first;
			   bool is_fixed_second_residue = res_2_pair.first;
			   lrc = make_link_restraints_for_link_ng(link_type, res_1, res_2,
								  is_fixed_first_residue,
								  is_fixed_second_residue,
								  do_trans_peptide_restraints,
								  geom);

			   // Rama restraints here!

			   status = true;

			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   return std::pair<bool, coot::restraints_container_t::link_restraints_counts>(status, lrc);

}

std::pair<bool, coot::restraints_container_t::link_restraints_counts>
coot::restraints_container_t::try_make_phosphodiester_link_ng(const coot::protein_geometry &geom,
						       std::pair<bool, mmdb::Residue *> res_1_pair,
						       std::pair<bool, mmdb::Residue *> res_2_pair,
						       bool do_trans_peptide_restraints) {

   mmdb::Residue *res_1 = res_1_pair.second;
   mmdb::Residue *res_2 = res_2_pair.second;
   std::string res_name_1(res_1->GetResName());
   std::string res_name_2(res_2->GetResName());
   bool status = false;
   link_restraints_counts lrc;

   if (util::is_standard_amino_acid_name(res_name_1)) {
      if (util::is_standard_amino_acid_name(res_name_2)) {
	 mmdb::Atom **residue_1_atoms = 0;
	 mmdb::Atom **residue_2_atoms = 0;
	 int n_residue_1_atoms;
	 int n_residue_2_atoms;
	 res_1->GetAtomTable(residue_1_atoms, n_residue_1_atoms);
	 res_2->GetAtomTable(residue_2_atoms, n_residue_2_atoms);
	 for (int iat_1=0; iat_1<n_residue_1_atoms; iat_1++) {
	    mmdb::Atom *at_1 = residue_1_atoms[iat_1];
	    std::string at_name_1(at_1->GetAtomName());
	    if (at_name_1 == " C  ") { // PDBv3 FIXE
	       std::string alt_conf_1(at_1->altLoc);
	       for (int iat_2=0; iat_2<n_residue_2_atoms; iat_2++) {
		  mmdb::Atom *at_2 = residue_2_atoms[iat_2];
		  std::string at_name_2(at_2->GetAtomName());
		  if (at_name_2 == " N  ") { // PDBv3 FIXE
		     std::string alt_conf_2(at_2->altLoc);
		     if (alt_conf_1 == alt_conf_2 || alt_conf_1.empty() || alt_conf_2.empty()) {

			// find_peptide_link_type_ng doesn't test the geometry -
			// just the residues types
			std::string link_type = find_peptide_link_type_ng(res_1, res_2, geom);
			if (! link_type.empty()) {
			   bool is_fixed_first_residue  = res_1_pair.first;
			   bool is_fixed_second_residue = res_2_pair.first;
			   lrc = make_link_restraints_for_link_ng(link_type, res_1, res_2,
								  is_fixed_first_residue,
								  is_fixed_second_residue,
								  do_trans_peptide_restraints,
								  geom);
			   status = true;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   return std::pair<bool, coot::restraints_container_t::link_restraints_counts>(status, lrc);

}

void
coot::restraints_container_t::make_polymer_links_ng(const coot::protein_geometry &geom,
						    bool do_rama_plot_restraints,
						    bool do_trans_peptide_restraints,
						    std::map<mmdb::Residue *, unsigned int> *residue_link_count_map_p,
						    std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > *residue_pair_link_set_p) {

   // make_polymer_links_ng uses the residues_vec (which is sorted in initialization)
   // so, only the passed active atoms are consider for polymer links

   std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > &residue_pair_link_set = *residue_pair_link_set_p;

   link_restraints_counts accum_links;

   // the residues are sorted by now, we can't have a forward link on the last residue of residues_vec
   //
   int n=residues_vec.size()-1;
   for (int i1=0; i1<n; i1++) {
      int i2=i1+1;

      if (residues_vec[i1].first && residues_vec[i2].first)
	 continue;

      mmdb::Residue *res_1 = residues_vec[i1].second;
      mmdb::Residue *res_2 = residues_vec[i2].second;

      std::string res_name_1(res_1->GetResName());
      std::string res_name_2(res_2->GetResName());
      if (res_name_1 == "HOH") continue; // break?
      if (res_name_2 == "HOH") continue;

      if (res_1->chain == res_2->chain) {
	 int serial_delta = res_2->index - res_1->index;

	 if (serial_delta == 1) {
	    // possibly polymer

	    int res_no_delta = res_2->GetSeqNum() - res_1->GetSeqNum();

	    std::pair<bool, link_restraints_counts> results;
	    // first test for normal peptide residue
	    //
	    if (res_no_delta == 1) {

	       // this can fail if there are missing link (C,N) atoms
	       results = try_make_peptide_link_ng(geom, residues_vec[i1], residues_vec[i2],
						  do_rama_plot_restraints,
						  do_trans_peptide_restraints);
	    } else {
	       if (res_no_delta == 0) { // insertion code?
		  std::string ins_code_1(res_1->GetInsCode());
		  std::string ins_code_2(res_2->GetInsCode());
		  if (ins_code_1 != ins_code_2) {
		     results = try_make_peptide_link_ng(geom, residues_vec[i1], residues_vec[i2],
							do_rama_plot_restraints,
							do_trans_peptide_restraints);
		  }
	       }
	    }

	    if (results.first) {
	       accum_links.add(results.second);
	    } else {
	       results = try_make_phosphodiester_link_ng(geom, residues_vec[i1], residues_vec[i2],
							 do_trans_peptide_restraints);
	    }

	    if (results.first) {

	       // make a note that those residues were linked (needed?)
	       //
	       std::map<mmdb::Residue *, unsigned int>::const_iterator it_1 =
		  residue_link_count_map_p->find(res_1);
	       if (it_1 == residue_link_count_map_p->end()) {
		  (*residue_link_count_map_p)[res_1] = 1; // a bit ugly.  Better syntax today?
	       } else {
		  residue_link_count_map_p->at(res_1)++;
	       }
	       std::map<mmdb::Residue *, unsigned int>::const_iterator it_2 =
		  residue_link_count_map_p->find(res_2);
	       if (it_2 == residue_link_count_map_p->end()) {
		  (*residue_link_count_map_p)[res_2] = 1;
	       } else {
		  residue_link_count_map_p->at(res_2)++;
	       }

	       std::pair<mmdb::Residue *, mmdb::Residue *> residue_pair_f(res_1, res_2);
	       std::pair<mmdb::Residue *, mmdb::Residue *> residue_pair_b(res_2, res_1);
	       residue_pair_link_set.insert(residue_pair_f);
	       residue_pair_link_set.insert(residue_pair_b);
	    }
	 }
      }
   }

   accum_links.report();

}

void
coot::restraints_container_t::fill_residue_link_count_map_ng(std::map<mmdb::Residue *, unsigned int> *residue_link_count_map_p) const {

   // not needed. Not yet at least.
}

// test for N and C linked peptide residue (or RNA/DNA)
bool
coot::restraints_container_t::is_fully_linked_ng(mmdb::Residue *r,
						 const std::map<mmdb::Residue *, unsigned int> &residue_link_count_map) const {

   std::map<mmdb::Residue *, unsigned int>::const_iterator it = residue_link_count_map.find(r);
   if (it == residue_link_count_map.end()) {
      return false;
   } else {
      unsigned int n_link_for_residue = it->second;
      if (n_link_for_residue == 2)
	 return true;
   }
   return false;

}

// non polymer links and carbohydrate links.
// And metal links
coot::restraints_container_t::link_restraints_counts
coot::restraints_container_t::make_other_types_of_link(const coot::protein_geometry &geom,
						       const std::map<mmdb::Residue *, unsigned int> &residue_link_count_map,
						       std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > residue_pair_link_set) {

   // I now doubt that residue_link_count_map is useful. Remove it?

   // ------------- Generic/non-polymer and carbohydrate Links -------------------

   link_restraints_counts lrc;

   if (false) {
      std::cout << "debug----------- fixed_atom_indices " << std::endl;
      std::set<int>::const_iterator it;
      for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
	 std::cout << "   fixed atoms: " << *it << std::endl;
   }

   class new_linked_residue {
   public:
      new_linked_residue(mmdb::Atom *at_1, mmdb::Atom *at_2,
			 mmdb::Residue *res_1, mmdb::Residue *res_2,
			 const std::string &link_type,
			 bool order_switch_flag) : at_1(at_1), at_2(at_2), res_1(res_1), res_2(res_2),
						   link_type(link_type), order_switch_flag(order_switch_flag) {}
      mmdb::Atom *at_1, *at_2;
      mmdb::Residue *res_1, *res_2;
      std::string link_type;
      bool order_switch_flag;
      bool match_p(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 if (r1 == res_1 && r2 == res_2) {
	    return true;
	 } else {
	    return (r2 == res_1 && r1 == res_2);
	 }
      }
   };

   class new_linked_residue_list {
   public:
      new_linked_residue_list() {}
      std::vector<new_linked_residue> nlr_vec;
      void insert(mmdb::Atom *at_1, mmdb::Atom *at_2,
		  mmdb::Residue *res_1, mmdb::Residue *res_2,
		  const std::string &link_type,
		  bool order_switch_flag) {
	 bool found = already_added(res_1, res_2);
	 if (! found) {
	    new_linked_residue nlr(at_1, at_2, res_1, res_2, link_type, order_switch_flag);
	    nlr_vec.push_back(nlr);
	 }
      }
      bool already_added(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 bool found = false;
	 for (std::size_t i=0; i<nlr_vec.size(); i++) {
	    if (nlr_vec[i].match_p(r1, r2)) {
	       found = true;
	       break;
	    }
	 }
	 return found;
      }
   };

   new_linked_residue_list nlrs;

   std::set<unsigned int> fixed_atom_flags_set;
   std::set<int>::const_iterator it;
   for (it=fixed_atom_indices.begin(); it!=fixed_atom_indices.end(); it++)
      fixed_atom_flags_set.insert(*it);
   contacts_by_bricks cb(atom, n_atoms, fixed_atom_flags_set);
   cb.set_dist_max(3.0);
   std::vector<std::set<unsigned int> > vcontacts;
   bool only_between_different_residues_flag = true;
   cb.find_the_contacts(&vcontacts, only_between_different_residues_flag);

   for (std::size_t i=0; i<vcontacts.size(); i++) {
      const std::set<unsigned int> &n_set = vcontacts[i];
      mmdb::Atom *at_1 = atom[i];
      // if (! is_fully_linked_ng(at_1->residue, residue_link_count_map)) {
      if (true) {
	 std::set<unsigned int>::const_iterator it;
	 for (it=n_set.begin(); it!=n_set.end(); it++) {
	    mmdb::Atom *at_2 = atom[*it];

	    mmdb::Residue *res_1 = at_1->residue;
	    mmdb::Residue *res_2 = at_2->residue;

	    if (res_1 == res_2) {
	       // these should not be here - they should be filtered out in contacts_by_bricks.
	       // At some state test if this is still needed.
	       continue;
	    }

	    if (at_2->residue->chain == at_1->residue->chain) {

	       // no links for non-moving bumping residues
	       //
	       if (fixed_atom_indices.find(i) != fixed_atom_indices.end())
		  if (fixed_atom_indices.find(*it) != fixed_atom_indices.end())
		     continue;

	       // at_2 might have a lower residue number than at_1.

	       // At this stage most of the contacts will be peptide contacts (generally speaking)
	       // How can I know if the residues of these atoms have already been linked
	       // (quite likely they will have been)

	       std::pair<mmdb::Residue *, mmdb::Residue *> test_residue_pair(res_1, res_2);

	       if (residue_pair_link_set.find(test_residue_pair) == residue_pair_link_set.end())  {
		  // it was a new pair.

		  if (false) {
		     std::cout << "failed to find these residues in the polyer-linked set: "
			       << residue_spec_t(res_1) << " "
			       << residue_spec_t(res_2) << std::endl;
		     std::cout << "Here are the pairs in the linked set " << std::endl;
		     std::set<std::pair<mmdb::Residue *, mmdb::Residue *> >::const_iterator it;
		     for (it=residue_pair_link_set.begin(); it!=residue_pair_link_set.end(); it++)
			std::cout << "   " << residue_spec_t(it->first) << " " << residue_spec_t(it->second) << std::endl;
		  }

		  mmdb::Residue *res_1 = at_1->residue;
		  mmdb::Residue *res_2 = at_2->residue;

		  if (nlrs.already_added(res_1, res_2))
		     continue;

		  // not sure that this is what I want now, really
		  std::pair<std::string, bool> lt = find_link_type_complicado(res_1, res_2, geom);
		  // Returns first (link_type) as "" if not found, second is order switch flag
		  if (! lt.first.empty()) {
		     std::cout << "Make a link restraints " << residue_spec_t(res_1) << " " << residue_spec_t(res_2)
			       << " with type " << lt.first << " and order switch " << lt.second
			       << std::endl;
		     nlrs.insert(at_1, at_2, res_1, res_2, lt.first, lt.second);
		  }
	       }
	    }
	 }
      }
   }

   return lrc;
}

std::map<mmdb::Residue *, unsigned int>
coot::restraints_container_t::make_link_restraints_ng(const coot::protein_geometry &geom,
						      bool do_rama_plot_restraints,
						      bool do_trans_peptide_restraints) {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   std::map<mmdb::Residue *, unsigned int> residue_link_count_map;
   std::set<std::pair<mmdb::Residue *, mmdb::Residue *> > residue_pair_link_set;
   auto tp_1 = std::chrono::high_resolution_clock::now();

   // residue 2 will be linked to residue 1 no matter how far apart they are
   //
   make_polymer_links_ng(geom, do_rama_plot_restraints, do_trans_peptide_restraints,
			 &residue_link_count_map,
			 &residue_pair_link_set);

   auto tp_2 = std::chrono::high_resolution_clock::now();
   make_other_types_of_link(geom, residue_link_count_map, residue_pair_link_set);


   auto tp_3 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   std::cout << "------------------- timing for make_link_restraints_ng(): : "
	     << d10 << " " << d21 << " " << d32 << " milliseconds " << std::endl;

   return residue_link_count_map;
}
