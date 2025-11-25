/* pli/protein-ligand-interactions.cc
 * 
 * Copyright 2010, 2011, 2012 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#include "coot-utils/coot-h-bonds.hh"
#include "protein-ligand-interactions.hh"

// Use coot::h_bonds class to generate ligands.  We do that by creating a synthetic
// temporary  molecule and atom selections.
//
std::vector<pli::fle_ligand_bond_t>
pli::get_fle_ligand_bonds(mmdb::Residue *ligand_res,
                          const std::vector<mmdb::Residue *> &residues,
                          mmdb::Manager *mol,
                          const std::map<std::string, std::string> &name_map,
                          const coot::protein_geometry &geom, int imol,
                          float water_dist_max,
                          float h_bond_dist_max) {

   std::vector<fle_ligand_bond_t> v; // returned value
   bool debug = true;

   if (debug) {
      std::cout << "::::::::::::::::::::: get_fle_ligand_bonds() inputs: " << std::endl;
      std::cout << "::::::::: ligand_res: " << ligand_res << " " << coot::residue_spec_t(ligand_res) <<
	 std::endl;
      std::cout << "::::::::: n residues: " << residues.size() << std::endl;
      for (unsigned int ires=0; ires<residues.size(); ires++)
	 std::cout << ":::::::::      residue: " << ires << " " << coot::residue_spec_t(residues[ires])
		   << " " << residues[ires]->GetResName() << std::endl;
      std::cout << "::::::::: mol: " << mol << std::endl;
      std::cout << "::::::::: name map size: " << name_map.size() << std::endl;
      std::map<std::string, std::string>::const_iterator it;
      for (it=name_map.begin(); it!=name_map.end(); it++)
	 std::cout << ":::::::::    name map: " << it->first << "->" << it->second << std::endl;
      std::cout << "::::::::: water_dist_max: " << water_dist_max << std::endl;
      std::cout << "::::::::: h_bond_dist_max: " << h_bond_dist_max << std::endl;
   }

   std::vector<mmdb::Residue *> rv = residues;
   rv.push_back(ligand_res);

   std::pair<bool, mmdb::Manager *> m = coot::util::create_mmdbmanager_from_residue_vector(rv, mol);
   coot::residue_spec_t ligand_spec(ligand_res);

   if (m.first) {
      int SelHnd_all = m.second->NewSelection(); // d
      int SelHnd_lig = m.second->NewSelection(); // d
      m.second->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      m.second->SelectAtoms(SelHnd_lig, 0, ligand_spec.chain_id.c_str(),
			    ligand_spec.res_no, ligand_spec.ins_code.c_str(),
			    ligand_spec.res_no, ligand_spec.ins_code.c_str(),
			    "*", "*", "*", "*");

      // -----------------------
      //   hydrogen bonds
      // -----------------------

      coot::h_bonds hb;
      // std::vector<coot::h_bond> hbonds = hb.get(SelHnd_lig, SelHnd_all, m.second, geom);
      std::pair<bool, int> status = hb.check_hb_status(SelHnd_lig, m.second, geom, imol);
      if (! status.first)
	 std::cout << "WARNING:: ===================== no HB status on atoms of ligand! ======="
		   << "=========" << std::endl;

      std::vector<coot::h_bond> hbonds = hb.get_mcdonald_and_thornton(SelHnd_lig, SelHnd_all, m.second,
                                                                      geom, imol, h_bond_dist_max);

      if (debug) {
         std::cout << "DEBUG:: get_mcdonald_and_thornton() returned " << hbonds.size() << " H-bonds" << std::endl;
         for (int i=0; i<hbonds.size(); i++) {
            std::cout << "      " << i << " " << hbonds[i] << std::endl;
         }
      }

      if (debug)
	 std::cout << "DEBUG:: get_fle_ligand_bonds from h_bonds class found "
		   << hbonds.size() << " H bonds." << std::endl;

      for (unsigned int i=0; i<hbonds.size(); i++) {
	 if (debug)
	    std::cout << "DEBUG:: in get_fle_ligand_bonds() hbond [" << i << "] "
		      << coot::atom_spec_t(hbonds[i].donor) << "...to... "
		      << coot::atom_spec_t(hbonds[i].acceptor) << " with ligand donor flag "
		      << hbonds[i].ligand_atom_is_donor << std::endl;

	 int bond_type = fle_ligand_bond_t::get_bond_type(hbonds[i].donor,
							  hbonds[i].acceptor,
							  hbonds[i].ligand_atom_is_donor);
	 // override these 2 if ligand atom is donor
	 //
	 mmdb::Atom      *ligand_atom = hbonds[i].acceptor;
	 mmdb::Atom *env_residue_atom = hbonds[i].donor;
	 double explict_H_bond_fudge_factor = 0.0; // for H-bonds with no Hs.
	 //
	 if (hbonds[i].ligand_atom_is_donor) {
	    ligand_atom = hbonds[i].donor;
	    env_residue_atom = hbonds[i].acceptor;
	 }

	 // OK, 20110511 new style, where is the hydrogen?
	 //
	 if (hbonds[i].has_hydrogen()) {
	    if (hbonds[i].ligand_atom_is_H()) {
	       ligand_atom      = hbonds[i].hb_hydrogen;
	       env_residue_atom = hbonds[i].acceptor;
	    } else {
	       ligand_atom      = hbonds[i].acceptor;
	       env_residue_atom = hbonds[i].hb_hydrogen;
	    }
	    explict_H_bond_fudge_factor = 1.2;
	 }

	 // This map no longer works because we don't pass ligand atom
	 // name any more (we pass a spec).
	 //

// 	 // Now, in 3D (pre-prodrgification) we don't have (polar) Hs on the ligand
// 	 // (but we do in 2D), the map allows transfer from the ligand O or N to the
// 	 // polar H in FLEV.
// 	 // 
// 	 std::map<std::string, std::string>::const_iterator it = name_map.find(ligand_atom->name);
// 	 if (it != name_map.end()) {
// 	    // If the map happens, that's presumably because we found a H
// 	    // attached to an N (or an H attached to an O), either way, we
// 	    // are sitting now on an H.
// 	    ligand_atom_name = it->second;
// 	 }

	 if (debug)
	    std::cout << "constructing fle ligand bond " << ligand_atom->name
		      << " " << bond_type << " " << hbonds[i].dist << " "
		      << coot::atom_spec_t(env_residue_atom) << " "
		      << env_residue_atom->GetResName()
		      << std::endl;

	 bool is_bond_to_water = false;
	 bool ok_to_add = true; // reset to 0 for certain waters
	 if (env_residue_atom) {
	    std::string env_residue_name(env_residue_atom->GetResName());
	    if (env_residue_name == "HOH") {
	       is_bond_to_water = true;
	       if (hbonds[i].dist > water_dist_max)
		  ok_to_add = false;
	    }
	 }
	 if (ok_to_add) {
	    // we want to pass the atom specifics (not just the environ residue spec)
	    //
	    // coot::fle_ligand_bond_t bond(ligand_atom_name, bond_type, hbonds[i].dist, res_spec);
	    //
	    fle_ligand_bond_t bond(coot::atom_spec_t(ligand_atom),
					 coot::atom_spec_t(env_residue_atom),
					 bond_type,
					 hbonds[i].dist+explict_H_bond_fudge_factor, is_bond_to_water);

	    std::string residue_name = ligand_atom->GetResName();
	    if (residue_name == "HOH")
	       bond.water_protein_length = find_water_protein_length(ligand_atom->residue, mol);

	    v.push_back(bond);
	 }
      }

      if (debug)
	 std::cout << ".... get_fle_ligand_bonds(): after h-bonds v.size() is " << v.size() << std::endl;

      // -----------------------
      //   covalent bonds
      // -----------------------

      // by distance and by LINK

      std::vector<fle_ligand_bond_t> covalent_bonds_d =
	 get_covalent_bonds_by_distance(m.second, SelHnd_lig, SelHnd_all, ligand_spec, geom);
      std::vector<fle_ligand_bond_t> covalent_bonds_l =
	 get_covalent_bonds_by_links(ligand_res, mol);
      std::vector<fle_ligand_bond_t> covalent_bonds = covalent_bonds_d;
      // now add in bonds if they are not in there already
      for (unsigned int i=0; i<covalent_bonds_l.size(); i++)
	 if (std::find(covalent_bonds_d.begin(), covalent_bonds_d.end(), covalent_bonds_l[i])
	     == covalent_bonds_d.end())
	    covalent_bonds.push_back(covalent_bonds_l[i]);

      // finally add the covent bonds to all-bonds
      for (unsigned int i=0; i<covalent_bonds.size(); i++)
	 v.push_back(covalent_bonds[i]);

      // -----------------------
      //   metal bonds
      // -----------------------

      std::vector<fle_ligand_bond_t> metal_bonds = get_metal_bonds(ligand_res, residues);
      for (unsigned int i=0; i<metal_bonds.size(); i++)
	 v.push_back(metal_bonds[i]);

      if (debug)
	 std::cout << ".... get_fle_ligand_bonds(): after metal bonds v.size() is " << v.size()
		   << std::endl;

      
      // -----------------------
      //   clean up 
      // -----------------------
      
      m.second->DeleteSelection(SelHnd_lig);
      m.second->DeleteSelection(SelHnd_all);
      delete m.second;

   }

   if (debug) {
      std::cout << ":::: get_fle_ligand_bonds returns these " << v.size()
		<< " bonds: " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) { 
	 std::cout << "   " << i << " :  " << v[i] << std::endl;
      }
   } 
   return v;
}




std::vector<pli::fle_ligand_bond_t>
pli::get_covalent_bonds_by_distance(mmdb::Manager *mol,
                                    int SelHnd_lig,
                                    int SelHnd_all,
                                    const coot::residue_spec_t &ligand_spec,
                                    const coot::protein_geometry &geom) {

   // 20101016, hydrogens don't make covalent bonds between ligands
   // and protein (or residues to residues in general) - so these get
   // filtered out.
   
   std::vector<fle_ligand_bond_t> v;
   int SelHnd_local = mol->NewSelection();
   mol->SelectAtoms(SelHnd_local, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
   mol->Select(SelHnd_local, mmdb::STYPE_ATOM, 0, ligand_spec.chain_id.c_str(),
	       ligand_spec.res_no, ligand_spec.ins_code.c_str(),
	       ligand_spec.res_no, ligand_spec.ins_code.c_str(),
	       "*", "*", "*", "*", mmdb::SKEY_XOR);

   // now find contacts:
   // 
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   mmdb::PPAtom lig_atom_selection = 0;
   int n_lig_atoms;
   mol->GetSelIndex(SelHnd_lig, lig_atom_selection, n_lig_atoms);

   mmdb::PPAtom other_atom_selection = 0;
   int n_other_atoms;
   mol->GetSelIndex(SelHnd_local, other_atom_selection, n_other_atoms);

   mmdb::realtype min_dist = 0.1;
   mmdb::realtype max_dist = 2.3; //  even S-S is shorter than this, I think
   mmdb::realtype cno_max_dist = 1.8;

   mol->SeekContacts(lig_atom_selection,   n_lig_atoms,
		     other_atom_selection, n_other_atoms,
		     min_dist, max_dist, // min, max distances
		     0,        // seqDist 0 -> in same res also
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > contacting_pairs_vec;

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i=0; i<n_contacts; i++) {
	    mmdb::Atom *at_1 =   lig_atom_selection[pscontact[i].id1];
	    mmdb::Atom *at_2 = other_atom_selection[pscontact[i].id2];

	    // move on if these are interacting atoms
	    std::string alt_conf_1 = at_1->altLoc;
	    std::string alt_conf_2 = at_2->altLoc;
	    if (!alt_conf_1.empty() && ! alt_conf_2.empty())
	       if (alt_conf_1 != alt_conf_2)
		  continue;

	    // 	    std::cout << "DEUBG:: Covalent test "
	    // 		      << coot::atom_spec_t(at_1) << "..."
	    // 		      << coot::atom_spec_t(at_2) << std::endl;
	    std::pair<mmdb::Residue *, mmdb::Residue *> pair(at_1->GetResidue(),
						   at_2->GetResidue());

	    std::string ele_1 = at_1->element;
	    std::string ele_2 = at_2->element;
	    if (ele_1 != " H") { 
	       if (ele_2 != " H") { 
		  clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
		  clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
		  double d = (pt_1-pt_2).lengthsq();

		  double dist_for_bond = max_dist;
		  if (((ele_1 == " C") || (ele_1 == " N") || (ele_1 == "O")) &&
		      ((ele_2 == " C") || (ele_2 == " N") || (ele_2 == "O")))
		     dist_for_bond = cno_max_dist;

		  if (d < dist_for_bond) { 

		     // only add this pair if it is not already in the list:
		     if (std::find(contacting_pairs_vec.begin(), contacting_pairs_vec.end(), pair) ==
			 contacting_pairs_vec.end()) {
			contacting_pairs_vec.push_back(pair);
			int bond_type = fle_ligand_bond_t::BOND_COVALENT;
			fle_ligand_bond_t bond(coot::atom_spec_t(at_1), // ligand
                                               coot::atom_spec_t(at_2), // env residue
                                               bond_type, d, false);
			v.push_back(bond);
		     }
		  }
	       }
	    }
	 }
      }
   }

   mol->DeleteSelection(SelHnd_local);
   return v;
}

std::vector<pli::fle_ligand_bond_t>
pli::get_covalent_bonds_by_links(mmdb::Residue *residue_ligand_p,
				  mmdb::Manager *mol) {

   std::vector<fle_ligand_bond_t> v;

   std::string residue_ligand_chain_id;
   int residue_ligand_res_no;
   std::string residue_ligand_ins_code;
   if (residue_ligand_p) { 
      residue_ligand_chain_id = residue_ligand_p->GetChainID();
      residue_ligand_res_no   = residue_ligand_p->GetSeqNum();
      residue_ligand_ins_code = residue_ligand_p->GetInsCode();
   }
   
   if (residue_ligand_p) { 
      if (mol) {
	 for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	    mmdb::Model *model_p = mol->GetModel(imod);
	    int n_links = model_p->GetNumberOfLinks();

	    for (int i_link=1; i_link<=n_links; i_link++) {
	       mmdb::PLink link = model_p->GetLink(i_link);

	       // Ligand then Protein
	       if (residue_ligand_chain_id == link->chainID1) {
		  if (residue_ligand_res_no == link->seqNum1) {
		     if (residue_ligand_ins_code == link->insCode1) {
			std::pair<coot::atom_spec_t, coot::atom_spec_t> linked_atoms = coot::link_atoms(link, model_p);
			mmdb::Atom *at_1 = coot::util::get_atom(linked_atoms.first,  mol);
			mmdb::Atom *at_2 = coot::util::get_atom(linked_atoms.second, mol);
			if (at_1 && at_2) {

			   // move on if these are interacting atoms
			   std::string alt_conf_1 = at_1->altLoc;
			   std::string alt_conf_2 = at_2->altLoc;
			   if (!alt_conf_1.empty() && ! alt_conf_2.empty())
			      if (alt_conf_1 != alt_conf_2)
				 continue;

			   double dist = coot::distance(at_1, at_2);
			   fle_ligand_bond_t b(linked_atoms.first, linked_atoms.second,
					       fle_ligand_bond_t::BOND_COVALENT,
					       dist, false);
			   v.push_back(b);
			}
		     }
		  }
	       }

	       // Protein then Ligand
	       if (residue_ligand_chain_id == link->chainID2) {
		  if (residue_ligand_res_no == link->seqNum2) {
		     if (residue_ligand_ins_code == link->insCode2) {
			std::pair<coot::atom_spec_t, coot::atom_spec_t> linked_atoms = coot::link_atoms(link, model_p);
			mmdb::Atom *at_1 = coot::util::get_atom(linked_atoms.first,  mol);
			mmdb::Atom *at_2 = coot::util::get_atom(linked_atoms.second, mol);
			if (at_1 && at_2) { 

			   // move on if these are interacting atoms
			   std::string alt_conf_1 = at_1->altLoc;
			   std::string alt_conf_2 = at_2->altLoc;
			   if (!alt_conf_1.empty() && ! alt_conf_2.empty())
			      if (alt_conf_1 != alt_conf_2)
				 continue;

			   double dist = coot::distance(at_1, at_2);
			   fle_ligand_bond_t b(linked_atoms.second, linked_atoms.first,
					       fle_ligand_bond_t::BOND_COVALENT,
					       dist, false);
			   v.push_back(b);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return v;
}


std::vector<pli::fle_ligand_bond_t>
pli::get_metal_bonds(mmdb::Residue *ligand_residue, const std::vector<mmdb::Residue *> &residues) {

   // a non-Hydrogen, non-Carbon ligand atom, that is.
   double max_dist_metal_to_ligand_atom = 3.5; // pass this parameter?
   
   std::vector<fle_ligand_bond_t> v;

   double best_dist_sqrd = max_dist_metal_to_ligand_atom * max_dist_metal_to_ligand_atom;
   mmdb::Atom *ligand_atom = NULL;
   mmdb::Atom *env_residue_atom = NULL;
   
   mmdb::PPAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_residue->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);
   for (unsigned int i=0; i<residues.size(); i++) { 
      if (is_a_metal(residues[i])) {
	 mmdb::PPAtom residue_atoms = 0;
	 int n_residue_atoms;
	 residues[i]->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int irat=0; irat<n_residue_atoms; irat++) {
	    for (int ilat=0; ilat<n_residue_atoms; ilat++) {

	       // move on if these are interacting atoms
	       mmdb::Atom *at_1 = ligand_residue_atoms[ilat];
	       mmdb::Atom *at_2 = residue_atoms[irat];
	       std::string alt_conf_1 = at_1->altLoc;
	       std::string alt_conf_2 = at_2->altLoc;
	       if (!alt_conf_1.empty() && ! alt_conf_2.empty())
		  if (alt_conf_1 != alt_conf_2)
		     continue;

	       std::string ele(residue_atoms[irat]->element);
	       if ((ele == " H") || (ele == " C")) { 
		  clipper::Coord_orth pt_1(ligand_residue_atoms[ilat]->x,
					   ligand_residue_atoms[ilat]->y,
					   ligand_residue_atoms[ilat]->z);
		  clipper::Coord_orth pt_2(residue_atoms[irat]->x,
					   residue_atoms[irat]->y,
					   residue_atoms[irat]->z);
		  double d2 = (pt_1-pt_2).clipper::Coord_orth::lengthsq();
		  if (d2 < best_dist_sqrd) {
		     best_dist_sqrd = d2;
		     ligand_atom = ligand_residue_atoms[ilat];
		     env_residue_atom = residue_atoms[irat];
		  }
	       }
	    }
	 }
	 if (best_dist_sqrd < max_dist_metal_to_ligand_atom * max_dist_metal_to_ligand_atom) {
	    fle_ligand_bond_t bond(coot::atom_spec_t(ligand_atom),
                                   coot::atom_spec_t(env_residue_atom),
                                   fle_ligand_bond_t::METAL_CONTACT_BOND,
                                   sqrt(best_dist_sqrd), false);
	    v.push_back(bond);
	 }
      }
   }

   return v;
}


// should be in coot-utils perhaps?
//
bool
pli::is_a_metal(mmdb::Residue *res) {

   bool r = 0;
   std::string res_name = res->GetResName();
   if (res_name == "MG")
      return 1;
   if (res_name == "CA")
      return 1;
   if (res_name == "MN")
      return 1;
   if (res_name == "FE")
      return 1;
   if (res_name == "K")
      return 1;
   if (res_name == "NA")
      return 1;
   if (res_name == "CO")
      return 1;
   if (res_name == "NI")
      return 1;
   if (res_name == "CU")
      return 1;
   if (res_name == "ZN")
      return 1;
   if (res_name == "RU")
      return 1;
   if (res_name == "PT")
      return 1;
   if (res_name == "AU")
      return 1;
   if (res_name == "AG")
      return 1;

   return 0;

}



// consider where a peptide is the ligand
std::vector<pli::fle_ligand_bond_t>
pli::protein_ligand_interactions(mmdb::Residue *ligand_residue_p, mmdb::Manager *mol,
                                 coot::protein_geometry *geom_p, int imol,
                                 float h_bond_dist_max) {

   float water_dist_max = 3.6; // pass this
   float residues_near_radius = 5.0; // pass this

   coot::residue_spec_t spec(ligand_residue_p);

   int SelHnd_all = mol->NewSelection(); // d
   int SelHnd_lig = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES,
		    "*", "*", "*", "*", "*");
   mol->SelectAtoms(SelHnd_lig, 0, spec.chain_id.c_str(),
		    spec.res_no, spec.ins_code.c_str(),
		    spec.res_no, spec.ins_code.c_str(),
		    "*", "*", "*", "*");

   std::vector<mmdb::Residue *> residues =
      coot::residues_near_residue(ligand_residue_p, mol, residues_near_radius);

   std::map<std::string, std::string> dummy_name_map;
   std::vector<fle_ligand_bond_t> bonds = get_fle_ligand_bonds(ligand_residue_p, residues, mol,
							       dummy_name_map, *geom_p, imol,
							       water_dist_max, h_bond_dist_max);

   coot::h_bonds hb;
   std::pair<bool, int> status = hb.check_hb_status(SelHnd_lig, mol, *geom_p, imol);
   if (! status.first)
      std::cout << "WARNING:: no HB status on atoms of ligand\n";
   std::vector<coot::h_bond> hbonds = hb.get_mcdonald_and_thornton(SelHnd_lig,
                                                                   SelHnd_all,
                                                                   mol, *geom_p, imol, h_bond_dist_max);

   for (unsigned int i=0; i<hbonds.size(); i++) {
      if (true)
	 std::cout << "DEBUG:: in process_ligand() hbond [" << i << "] donor "
		   << coot::atom_spec_t(hbonds[i].donor) << "...to... "
		   << coot::atom_spec_t(hbonds[i].acceptor) << " with ligand donor flag "
		   << hbonds[i].ligand_atom_is_donor << std::endl;

      // override these 2 if ligand atom is donor
      //
      mmdb::Atom      *ligand_atom = hbonds[i].acceptor;
      mmdb::Atom *env_residue_atom = hbonds[i].donor;
      if (hbonds[i].ligand_atom_is_donor) {
	 ligand_atom = hbonds[i].donor;
	 env_residue_atom = hbonds[i].acceptor;
      }
   }

   // Other stuff

   mol->DeleteSelection(SelHnd_all);
   mol->DeleteSelection(SelHnd_lig);

   return bonds;
}


// return 100 if no other contact found (strange!)
// 
double
pli::find_water_protein_length(mmdb::Residue *ligand_residue, mmdb::Manager *mol) {

   double dist = 100;

   double dist_sqrd = dist * dist;
   double dist_sqrd_init = dist_sqrd;
   
   mmdb::PPAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_residue->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 if (ligand_residue != residue_p) {
	    std::string residue_name(residue_p->GetResName());
	    if (residue_name != "HOH") { 
	       mmdb::PPAtom residue_atoms = 0;
	       int n_residue_atoms;
	       residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (int il=0; il<n_ligand_residue_atoms; il++) {
		  for (int irat=0; irat<n_residue_atoms; irat++) {
		     std::string ele(residue_atoms[irat]->element);
		     if ((ele == " O") || (ele == " N")) { 
			clipper::Coord_orth pt_1(ligand_residue_atoms[il]->x,
						 ligand_residue_atoms[il]->y,
						 ligand_residue_atoms[il]->z);
			clipper::Coord_orth pt_2(residue_atoms[irat]->x,
						 residue_atoms[irat]->y,
						 residue_atoms[irat]->z);
			double d2 = (pt_1-pt_2).clipper::Coord_orth::lengthsq();
			if (d2 < dist_sqrd) {
			   dist_sqrd = d2;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (dist_sqrd < dist_sqrd_init)  // usually is.
      dist = sqrt(dist_sqrd);

   return dist;
}
