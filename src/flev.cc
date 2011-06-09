/* src/flev.cc
 * 
 * Copyright 2010 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

/*  ----------------------------------------------------------------------- */
/*               Flattened Ligand Environment View  Interface               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

// For reasons I don't understand, this should come near the top of
// includes, otherwise we get RDKit dcgettext() include file problems.
//
#ifdef MAKE_ENTERPRISE_TOOLS
#include "lbg.hh"
#endif

#include "c-interface-ligands.hh"
#include "mmdb-extras.h"
#include "mmdb.h"

#include "graphics-info.h"
#include "c-interface.h"
#include "lbg-graph.hh"

#include "coot-h-bonds.hh"

#include "flev.hh"

#include "cc-interface.hh" // for add_animated_ligand_interactions


std::ostream&
coot::operator<< (std::ostream& s, const pi_stacking_instance_t &stack) {

   s << "[" << stack.type << " " << coot::residue_spec_t(stack.res) << " "
     << stack.overlap_score << " ligand-atom-name :"
     <<  stack.ligand_cationic_atom_name
     << ": ";
   for (unsigned int i=0; i<stack.ligand_ring_atom_names.size(); i++)
      s << "  :" << stack.ligand_ring_atom_names[i] << ":   " ;
   s << "]";
   return s;
}

int
sprout_hydrogens(int imol,
		 const char *chain_id,
		 int res_no,
		 const char *ins_code) {
   
   int r = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<bool, std::string> r_add =
	 g.molecules[imol].sprout_hydrogens(chain_id, res_no, ins_code, *(g.Geom_p()));
      r = r_add.first;
      if (r)
	 graphics_draw();
      else
	 info_dialog(r_add.second.c_str());
   }
   return r;
}


void fle_view_internal(int imol, const char *chain_id, int res_no, const char *ins_code, 
		       int imol_ligand_fragment, 
		       const char *prodrg_output_flat_mol_file_name,
		       const char *prodrg_output_flat_pdb_file_name,
		       const char *prodrg_output_3d_pdb_file_name,
		       const char *prodrg_output_dict_cif_file_name) {
   
   float residues_near_radius = 4.5;
   float water_dist_max = 3.25;
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();
   
   // Plan:
   //
   // read in the prodrg output flat pdb file.
   // 
   // lsq match from the residue at imol/chain_id/resno/ins_code to that the first
   // residue in the prodrg flat molecule.  Get the RTop.
   //
   // Get the residues in imol that are within 4A (say) of the residue
   // at imol/chain_id/resno/ins_code.  Find the centres of these
   // residue (construct a helper class)
   // 
   // class fle_residues_helper_t { 
   //    centre
   //    residue spec
   //    residue name
   // }
   //
   //
   // Rotate centres by RTop
   //
   // Write a file of lines like:
   //
   // RES flat_pos_x flat_pos_y TYPE spec
   // 
   // (flat_pos_x, flat_pos_y are in the coordinate system of the
   // flattened ligand)
   //
   // Additionally, a RES card can have the immediately following record:
   //
   // BOND to_atom_x to_atom_y colour_name
   //

   atom_selection_container_t flat = get_atom_selection(prodrg_output_flat_pdb_file_name, 1);
   if (flat.read_success) {
      if (is_valid_model_molecule(imol_ligand_fragment)) {

	 CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 CResidue  *res_ref = coot::util::get_residue(chain_id, res_no, ins_code, mol);
	 CResidue *flat_res = coot::util::get_residue("", 1, "", flat.mol); // prodrg attribs

	 std::string res_name_res_ref; // just in case res_ref is NULL
	 if (res_ref)
	    res_name_res_ref = res_ref->GetResName();

	 // We only want to read a cif dictionary if is is not one of
	 // the standard ones (this is because "standard" residues
	 // have corrrect group for linking (typically "L-peptide")
	 // but the file from PRODRG is always "monomer" - and that
	 // means bad news if we subsequently try to look for covalent
	 // bonds.
	 //
	 if (! (coot::standard_residue_name_p(res_name_res_ref)))
	    handle_cif_dictionary(prodrg_output_dict_cif_file_name);
	 
	 
	 if (res_ref) { 
	    
	    std::string ligand_res_name(res_ref->GetResName());
	    CMMDBManager *ligand_mol =
	       graphics_info_t::molecules[imol_ligand_fragment].atom_sel.mol;
	    int every_nth = 1;
	    std::vector<coot::lsq_range_match_info_t> matches;
	    coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
					       COOT_LSQ_ALL);
	    matches.push_back(match);
	    std::pair<short int, clipper::RTop_orth> lsq_mat = 
	       coot::util::get_lsq_matrix(flat.mol, ligand_mol, matches, every_nth);
	    
	    if (lsq_mat.first) { 
	       std::vector<CResidue *> residues =
		  coot::residues_near_residue(res_ref, mol, residues_near_radius);
	       
	       // residues needs to be filtered to remove waters that
	       // are not connected to a protein atom.

	       // 20101228 water_max_dist was 3.6, not we tweak it to
	       // 3.25.  Should be a user-setable param.
	       // 
	       std::vector<CResidue *> filtered_residues =
		  coot::filter_residues_by_solvent_contact(res_ref, mol, residues, water_dist_max);
	       
	       // for the atoms in the ligand only, for the moment.
	       // More would require a class containing with and
	       // without solvent accessibilites for each residue.
	       // Maybe that can be another function.  And that would
	       // need to consider neighbours of neighbours, perhaps
	       // done with a larger radius.
	       //
	       coot::dots_representation_info_t dots;
	       std::vector<std::pair<coot::atom_spec_t, float> > s_a_v = 
		  dots.solvent_accessibilities(res_ref, filtered_residues);
	       
	       // for the ligand environment residues:
	       std::vector<coot::solvent_exposure_difference_helper_t> sed = 
		  dots.solvent_exposure_differences(res_ref, filtered_residues);
	       
	       if (0)
		  for (unsigned int i=0; i<s_a_v.size(); i++)
		     std::cout << "   " << i << " " << s_a_v[i].first << " "
			       << s_a_v[i].second << std::endl;

	       //
	       std::map<std::string, std::string> name_map =
		  coot::make_flat_ligand_name_map(flat_res);

	       // Now what are the bonds between the ligand and residues?
	       // 
	       // a vector of fle-ligand-bond which contain
	       // ligand-atom-name residue-spec and bond type
	       // (acceptor/donor).
	       //
	       // std::vector<coot::fle_ligand_bond_t> bonds_to_ligand =
	       // coot::get_fle_ligand_bonds(res_ref, filtered_residues, name_map);

	       // using new (20100522) h_bond class (based on geometry)
	       //
 	       std::vector<coot::fle_ligand_bond_t> bonds_to_ligand = 
 		  coot::get_fle_ligand_bonds(res_ref, filtered_residues, mol,
 					     name_map, *geom_p, water_dist_max);

	       add_animated_ligand_interactions(imol, bonds_to_ligand);

	       if (1) 
		  std::cout << "Found ================== " << bonds_to_ligand.size()
			    << " ==================== bonds to ligand " << std::endl;

	       if (filtered_residues.size()) {
		  std::vector<coot::fle_residues_helper_t> centres(filtered_residues.size());
		  for (unsigned int ires=0; ires<filtered_residues.size(); ires++) { 
		     CResidue *res_copy = coot::util::deep_copy_this_residue(filtered_residues[ires]);
		     std::string res_name = filtered_residues[ires]->GetResName();
		     coot::util::transform_atoms(res_copy, lsq_mat.second);
		     std::pair<bool, clipper::Coord_orth> c =
			coot::util::get_residue_centre(res_copy);
		     if (c.first) {
			if (0) 
			   std::cout << "DEBUG:: creating fle_centre with centre "
				     << c.second.format() << std::endl;
			coot::fle_residues_helper_t fle_centre(c.second,
							       coot::residue_spec_t(filtered_residues[ires]),
							       res_name);
			centres[ires] = fle_centre;
		     } else {
			std::cout << "WARNING:: failed to get residue centre for "
				  << coot::residue_spec_t(res_copy) << std::endl;
		     }
		     delete res_copy;
		  }
		  
		  std::pair<bool, coot::dictionary_residue_restraints_t> p = 
		     geom_p->get_monomer_restraints(ligand_res_name);
		  
		  if (! p.first) {
		     std::cout << "WARNING:: fle_view_internal(): "
			       << "Failed to get monomer_restraints for PRODRG residue"
			       << std::endl;
		  } else {

		     // ----------- residue infos ----------
		     // 
		     coot::pi_stacking_container_t pi_stack_info(p.second, filtered_residues, res_ref);
		     write_fle_centres(centres, bonds_to_ligand, sed, pi_stack_info, flat_res);
		     
		     // ----------- ligand atom infos ------
		     // 
		     coot::flev_attached_hydrogens_t ah(p.second);
		     ah.cannonballs(res_ref, prodrg_output_3d_pdb_file_name, p.second);
		     ah.distances_to_protein(res_ref, mol);
		     write_ligand_atom_accessibilities(s_a_v, ah, flat_res);

		  }
	       }
	    }
	 }
      }
   } 
}

void fle_view_with_rdkit(int imol, const char *chain_id, int res_no,
			 const char *ins_code, float residues_near_radius) { 

#ifndef MAKE_ENTERPRISE_TOOLS
# else

   double weight_for_3d_distances = 0.4; // for 3d distances
   double water_dist_max = 3.25;
   
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();

   if (is_valid_model_molecule(imol)) { 
      CResidue  *res_ref = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      CMMDBManager *mol_for_res_ref = g.molecules[imol].atom_sel.mol;
      if (res_ref) {
	 std::string ligand_res_name(res_ref->GetResName());

	 std::pair<bool, coot::dictionary_residue_restraints_t> p = 
	    geom_p->get_monomer_restraints_at_least_minimal(ligand_res_name);
	 
	 if (! p.first) {
	    std::cout << "WARNING:: fle_view_with_rdkit(): "
		      << "Failed to get monomer_restraints for ligand of type "
		      << ligand_res_name << std::endl;
	 } else {
	    std::vector<CResidue *> residues =
	       coot::residues_near_residue(res_ref, mol_for_res_ref, residues_near_radius);
	 
	    // residues needs to be filtered to remove waters that
	    // are not connected to a protein atom.
	 
	    std::vector<CResidue *> filtered_residues =
	       coot::filter_residues_by_solvent_contact(res_ref, mol_for_res_ref, residues, 3.25);

	    // for the atoms in the ligand only, for the moment.
	    // More would require a class containing with and
	    // without solvent accessibilites for each residue.
	    // Maybe that can be another function.  And that would
	    // need to consider neighbours of neighbours, perhaps
	    // done with a larger radius.
	    //
	    coot::dots_representation_info_t dots;
	    std::vector<std::pair<coot::atom_spec_t, float> > s_a_v = 
	       dots.solvent_accessibilities(res_ref, filtered_residues);
	 
	    // for the ligand environment residues:
	    std::vector<coot::solvent_exposure_difference_helper_t> sed = 
	       dots.solvent_exposure_differences(res_ref, filtered_residues);

	    try {
	       // can throw a std::runtime_error
	       RDKit::RWMol rdkm = coot::rdkit_mol(res_ref, *g.Geom_p());

	       // assign atom names
	       if (rdkm.getNumAtoms() < res_ref->GetNumberOfAtoms()) {
		  std::cout << "oops failure to construct rdkit molecule " << std::endl;
	       } else {
		  PPCAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_ref->GetAtomTable(residue_atoms, n_residue_atoms);

		  // polar Hs only, that is - need new function here.
		  // (can throw a std::exception)
		  std::cout << "DEBUG:: calling remove_non_polar_Hs() " << std::endl;
		  coot::undelocalise(&rdkm);
		  coot::assign_formal_charges(&rdkm);
		  coot::remove_non_polar_Hs(&rdkm);
		  std::cout << "DEBUG::    done remove_non_polar_Hs() " << std::endl;

		  int mol_2d_depict_conformer =
		     coot::add_2d_conformer(&rdkm, weight_for_3d_distances);
		  lig_build::molfile_molecule_t m =
		     coot::make_molfile_molecule(rdkm, mol_2d_depict_conformer);

		  CResidue *residue_flat = coot::make_residue(rdkm, mol_2d_depict_conformer);
		  CMMDBManager *mol_for_flat_residue =
		     coot::util::create_mmdbmanager_from_residue(NULL, residue_flat);

		  if (0) 
		     for (unsigned int iat=0; iat<m.atoms.size(); iat++)
			std::cout << iat << "  " << m.atoms[iat] << std::endl;
	       
		  std::string view_name = "Molecule ";
		  view_name += coot::util::int_to_string(imol);
		  view_name += " ";
		  view_name += chain_id;
		  view_name += coot::util::int_to_string(res_no);
		  view_name += ins_code;
		  view_name += "    code: ";
		  view_name += ligand_res_name;

		  std::pair<bool, coot::residue_spec_t> ligand_spec_pair(1, coot::residue_spec_t(res_ref));

		  bool use_graphics_flag = graphics_info_t::use_graphics_interface_flag;
		  bool stand_alone_flag = 0; // no, it isn't from here.
		  lbg_info_t *lbg_local_p = lbg(m, ligand_spec_pair,
						NULL, view_name, ligand_res_name, imol,
						use_graphics_flag, stand_alone_flag);

 		  std::map<std::string, std::string> name_map =
 		     coot::make_flat_ligand_name_map(res_ref);
	       
 		  std::vector<coot::fle_ligand_bond_t> bonds_to_ligand =
 		     coot::get_fle_ligand_bonds(res_ref, filtered_residues,
						mol_for_res_ref, name_map, *geom_p, water_dist_max);

		  add_animated_ligand_interactions(imol, bonds_to_ligand);

		  std::vector<coot::fle_residues_helper_t> res_centres =
		     coot::get_flev_residue_centres(res_ref,
						    mol_for_res_ref,
						    filtered_residues,
						    mol_for_flat_residue);
		  std::vector<int> add_reps_vec =
		     coot::make_add_reps_for_near_residues(filtered_residues, imol);

		  if (0) { 
		     std::cout << "------------- in flev: centres.size() is "
			       << res_centres.size() << std::endl;
		     for (unsigned int ic=0; ic<res_centres.size(); ic++)
			std::cout << "   " << ic << "  " << res_centres[ic]
				  << std::endl;
		  }


		  // ----------- residue infos ----------
		  // 
		  coot::pi_stacking_container_t
		     pi_stack_info(p.second, filtered_residues, res_ref);
		     
		  // ----------- ligand atom infos ------
		  // 
 		  coot::flev_attached_hydrogens_t ah(p.second);
		  // ah.cannonballs(res_ref, mol_for_res_ref, p.second);
 		  ah.distances_to_protein_using_correct_Hs(res_ref, mol_for_res_ref, *geom_p);

		  // ------------ show it -----------------
		  // 
 		  lbg_local_p->annotate(s_a_v, res_centres, add_reps_vec, bonds_to_ligand, sed,
					ah, pi_stack_info, p.second);
		  delete mol_for_flat_residue;
	       }
	    }
	    catch (std::runtime_error rte) {
	       std::cout << "ERROR:: (runtime error) in fle_view_with_rdkit(): "
			 << rte.what() << std::endl;
	    } 
	    catch (std::exception e) {
	       std::cout << "ERROR (exception) in fle_view_with_rdkit(): " << e.what() << std::endl;
	    } 
	 }
      }
   }
#endif // MAKE_ENTERPRISE_TOOLS   
}

std::vector<int>
coot::make_add_reps_for_near_residues(std::vector<CResidue *> filtered_residues,
				      int imol) {

   std::vector<int> v(filtered_residues.size());

   int rep_type = coot::SIMPLE_LINES;
   int bonds_box_type = coot::NORMAL_BONDS;
   int bond_width = 8;
   int draw_hydrogens_flag = 1;
   for (unsigned int i=0; i<filtered_residues.size(); i++) {
      v[i] = additional_representation_by_attributes(imol,
						     filtered_residues[i]->GetChainID(),
						     filtered_residues[i]->GetSeqNum(),
						     filtered_residues[i]->GetSeqNum(),
						     filtered_residues[i]->GetInsCode(),
						     rep_type, bonds_box_type, bond_width,
						     draw_hydrogens_flag);
      set_show_additional_representation(imol, v[i], 0); // undisplay it
   }
   return v;
}

void
coot::add_animated_ligand_interactions(int imol,
				       const std::vector<coot::fle_ligand_bond_t> &ligand_bonds) {

   for (unsigned int i=0; i<ligand_bonds.size(); i++) {
      std::cout << "Here....  adding animated ligand interaction " << i << std::endl;
      add_animated_ligand_interaction(imol, ligand_bonds[i]);
   }
}



std::ostream &
coot::operator<<(std::ostream &s, fle_residues_helper_t fler) {

   s << fler.is_set;
   if (fler.is_set) {
      s << " " << fler.centre.format() << " " << fler.spec << " " << fler.residue_name;
   }
   return s;
}


std::vector<coot::fle_residues_helper_t>
coot::get_flev_residue_centres(CResidue *residue_ligand_3d,
			       CMMDBManager *mol_containing_residue_ligand, 
			       std::vector<CResidue *> residues,
			       CMMDBManager *flat_mol) {

   std::vector<coot::fle_residues_helper_t> centres;

   if (flat_mol) { 

      // get the lsq matrix that maps the ligand in 3D onto the flat ligand
      int res_no = residue_ligand_3d->GetSeqNum();
      std::string chain_id = residue_ligand_3d->GetChainID();
      int every_nth = 1;
      std::vector<coot::lsq_range_match_info_t> matches;
      coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
 					 COOT_LSQ_ALL);
       matches.push_back(match);
       std::pair<short int, clipper::RTop_orth> lsq_mat = 
	  coot::util::get_lsq_matrix(flat_mol, mol_containing_residue_ligand, matches, every_nth);
      // Now make the residues
      
      // std::vector<coot::fle_residues_helper_t> centres(residues.size());
      centres.resize(residues.size());
      for (unsigned int ires=0; ires<residues.size(); ires++) { 
	 CResidue *res_copy = coot::util::deep_copy_this_residue(residues[ires]);
	 std::string res_name = residues[ires]->GetResName();
	 coot::util::transform_atoms(res_copy, lsq_mat.second);
	 std::pair<bool, clipper::Coord_orth> c =
	    coot::util::get_residue_centre(res_copy);
	 if (c.first) {
	    coot::fle_residues_helper_t fle_centre(c.second,
						   coot::residue_spec_t(residues[ires]),
						   res_name);
	    centres[ires] = fle_centre;
	 } else {
	    std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
		      << coot::residue_spec_t(res_copy) << std::endl;
	 }
	 delete res_copy;
      }
   }
   return centres;
} 


std::map<std::string, std::string>
coot::make_flat_ligand_name_map(CResidue *flat_res) {

   double bond_to_H_dist = 1.1;

   double b2Hd2 = bond_to_H_dist * bond_to_H_dist;
   std::map<std::string, std::string> map;
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   flat_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      CAtom *at_i = residue_atoms[iat];
      std::string ele_i = at_i->element;
      clipper::Coord_orth pt_i(at_i->x, at_i->y, at_i->z);
      if (ele_i == " H") { 
	 for (int jat=0; jat<n_residue_atoms; jat++) {
	    if (iat != jat) { 
	       CAtom *at_j = residue_atoms[jat];
	       std::string ele_j = at_j->element;
	       if (ele_j != " H") { 
		  clipper::Coord_orth pt_j(at_j->x, at_j->y, at_j->z);
		  if ((pt_i - pt_j).lengthsq() < b2Hd2) {
		     map[at_j->name] = at_i->name;
		     break;
		  } 
	       }
	    }
	 }
      }
   }

//    std::map<std::string, std::string>::const_iterator it;
//    std::cout << "=== name map: === " << std::endl;
//    for (it=map.begin(); it!=map.end(); it++) {
//       std::cout << "  :" << it->first << ": :" << it->second << ":" << std::endl;
//    }
//    std::cout << "=== === " << std::endl;

   return map;
}


coot::pi_stacking_container_t::pi_stacking_container_t(const coot::dictionary_residue_restraints_t &monomer_restraints,
						       
						       const std::vector<CResidue *> &residues,
						       CResidue *res_ref) {

   bool debug = 1;

   // --------------------------------------------------------------------
   //          ligand ring systems
   // --------------------------------------------------------------------
   
   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   std::vector<std::vector<std::string> > ring_list = get_ligand_aromatic_ring_list(monomer_restraints);

   // float pi_overlap_thresh = 0.0015; // play value
   float pi_pi_overlap_thresh = 0.1;    // HLZ in 3LTW
   float pi_cation_overlap_thresh = 30;  // ZZG in 2wot is 27, close but spurious interaction,
                                         // a bit more than that then.

   
   for (unsigned int iring=0; iring<ring_list.size(); iring++) {
      try {
	 std::pair<clipper::Coord_orth, clipper::Coord_orth> ligand_ring_pi_pts = 
	    get_ring_pi_centre_points(ring_list[iring], res_ref);

	 if (debug) {
	    std::cout << "========= ligand ring ";
	    for (unsigned int iat=0; iat<ring_list[iring].size(); iat++)
	       std::cout << ring_list[iring][iat] << "  ";
	    
	    std::cout << " ====== points " << ligand_ring_pi_pts.first.format() << " "
		      << ligand_ring_pi_pts.second.format() << std::endl;
	 }
	 
	 for (unsigned int ires=0; ires<residues.size(); ires++) {

	    if (debug) {
	       std::string res_name(residues[ires]->GetResName());
	       std::cout << "==== Environment residue " << coot::residue_spec_t(residues[ires])
			 << " " << res_name << std::endl;
	    }

	    // return a pair that is the score and the stacking type
	    std::pair<float, pi_stacking_instance_t::stacking_t> pi_overlap_1 =
	       get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.first);
	    std::pair<float, pi_stacking_instance_t::stacking_t> pi_overlap_2 =
	       get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.second);

	    if (debug) 
	       std::cout << "   protein cation:ligand ring: Overlaps:  score "
			 << pi_overlap_1.first << " type: " << pi_overlap_1.second << "  score: "
			 << pi_overlap_2.first << " type: " << pi_overlap_2.second << std::endl;

	    float thresh = -1;
	    if (pi_overlap_1.second == coot::pi_stacking_instance_t::PI_PI_STACKING)
	       thresh = pi_pi_overlap_thresh;
	    if (pi_overlap_1.second == coot::pi_stacking_instance_t::PI_CATION_STACKING)
	       thresh = pi_cation_overlap_thresh;
	    
	    if (pi_overlap_1.first > thresh) {
	       coot::pi_stacking_instance_t st(residues[ires],
					       pi_overlap_1.second,
					       ring_list[iring]);
	       st.overlap_score = pi_overlap_1.first;
	       std::cout << "adding a stacking " << st << std::endl;
	       stackings.push_back(st);
	    }
	    if (pi_overlap_2.first > thresh) {
	       coot::pi_stacking_instance_t st(residues[ires],
					       pi_overlap_2.second,
					       ring_list[iring]);
	       st.overlap_score = pi_overlap_2.first;
	       std::cout << "adding a stacking " << st << std::endl;
	       stackings.push_back(st);
	    }
	 }
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   // --------------------------------------------------------------------
   //          ligand cations (-> protein ring systems)
   // --------------------------------------------------------------------
   
   // Now, we need to test the ligand cations against the environment
   // residues that have pi system (not ligand ring systems).
   //
   // OK, so what are the ligand cations?

   std::vector<std::pair<std::string, clipper::Coord_orth> > cation_points =
      get_ligand_cations(res_ref, monomer_restraints);

   for (unsigned int icat=0; icat<cation_points.size(); icat++) {
      for (unsigned int ires=0; ires<residues.size(); ires++) {
	 // get_ligand_cation_residue_pi_overlap
	 float score = get_pi_overlap_to_ligand_cation(residues[ires], cation_points[icat].second);

	 // std::cout << "debug:: in pi_stacking_container_t constructor cation point " << icat << " of "
	 // << cation_points.size() << "  score: " << score << " c.f. " << pi_overlap_thresh
	 // << std::endl;

	 if (score > pi_cation_overlap_thresh) { 
	    // add a stacking to stackings.
	    coot::pi_stacking_instance_t stacking(residues[ires], cation_points[icat].first);
	    stacking.overlap_score = score;
	    stackings.push_back(stacking);
	 }
      }
   }

}

std::vector<std::vector<std::string> >
coot::pi_stacking_container_t::get_ligand_aromatic_ring_list(const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::vector<std::string> > ring_list = monomer_restraints.get_ligand_aromatic_ring_list();
   return ring_list;
}


// by search through res_ref
std::vector<std::pair<std::string, clipper::Coord_orth> >
coot::pi_stacking_container_t::get_ligand_cations(CResidue *res_ref,
						 const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   std::vector<std::pair<std::string, clipper::Coord_orth> > v;
   int n_residue_atoms;
   PPCAtom residue_atoms = NULL;
   res_ref->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) { 
      std::string ele(residue_atoms[iat]->element);
      if (ele == " N") {
	 // how many bonds does this N have?
	 int n_bonds = 0;
	 std::string atom_name(residue_atoms[iat]->name);
	 for (unsigned int ibond=0; ibond<monomer_restraints.bond_restraint.size(); ibond++) { 
	    if ((monomer_restraints.bond_restraint[ibond].atom_id_1_4c() == atom_name) || (monomer_restraints.bond_restraint[ibond].atom_id_2_4c() == atom_name)) {
	       if (monomer_restraints.bond_restraint[ibond].type() == "single")
		  n_bonds++;
	       if (monomer_restraints.bond_restraint[ibond].type() == "double")
		  n_bonds += 2;
	       if (monomer_restraints.bond_restraint[ibond].type() == "triple")
		  n_bonds += 3;
	    } 
	 }

	 if (n_bonds > 3) { // i.e. 4
	    clipper::Coord_orth pt(residue_atoms[iat]->x,
				   residue_atoms[iat]->y,
				   residue_atoms[iat]->z);
	    std::pair<std::string, clipper::Coord_orth> p(atom_name, pt);
	    v.push_back(p);
	 }
      }
   }
   if (0) { // debug
      std::cout << "debug:: in coot::pi_stacking_container_t::get_ligand_cations() found " << v.size()
		<< " ligand cation points " << std::endl;
      for (unsigned int i=0; i<v.size(); i++)
	 std::cout << "   :" << v[i].first << ":  " << v[i].second.format() << std::endl;
   }
   return v;
}


// can throw an exception because it calls get_ring_pi_centre_points()
//
// should return the stacking type, e.g. PI_CATION_STACKING.
// 
std::pair<float, coot::pi_stacking_instance_t::stacking_t>
coot::pi_stacking_container_t::get_pi_overlap_to_ligand_ring(CResidue *res,
							     const clipper::Coord_orth &ligand_pi_point) const {

   float pi_pi_score = 0;
   float pi_cation_score = 0;
   
   std::string res_name(res->GetResName());
   coot::pi_stacking_instance_t::stacking_t stacking_type = coot::pi_stacking_instance_t::PI_PI_STACKING;

   // First test all the ring systems in the residue
   // 
   std::vector<std::vector<std::string> > atom_names = ring_atom_names(res_name);
   for (unsigned int iring=0; iring<atom_names.size(); iring++) {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> residue_pi_points =
	 get_ring_pi_centre_points(atom_names[iring], res);
      float score_1 = overlap_of_pi_spheres(ligand_pi_point, residue_pi_points.first,
					    0.78, -1, 0.78, -1);
      float score_2 = overlap_of_pi_spheres(ligand_pi_point, residue_pi_points.second,
					    0.78, -1, 0.78, -1);
      if (score_1 > pi_pi_score)
	 pi_pi_score = score_1;
      if (score_2 > pi_pi_score)
	 pi_pi_score = score_2;
   }

   // Now test all the cation-pi interactions
   // 
   std::vector<clipper::Coord_orth> cation_atom_point = get_cation_atom_positions(res);

   // std::cout << "DEBUG:: there are " << cation_atom_point.size()
   // << " cation atom points" << std::endl;
   
   for (unsigned int icat=0; icat<cation_atom_point.size(); icat++) {
      pi_cation_score += overlap_of_cation_pi(ligand_pi_point, cation_atom_point[icat]);
   }

   float score = pi_pi_score;
   if (pi_cation_score > pi_pi_score) {
      score = pi_cation_score;
      stacking_type = coot::pi_stacking_instance_t::PI_CATION_STACKING;
   }
   
   return std::pair<float, coot::pi_stacking_instance_t::stacking_t> (score, stacking_type);
}

// Return the best score of the ligand cation overlap to any of the
// ring systems in res (if any of course, typically this function just
// falls through returning 0.0).
// 
float
coot::pi_stacking_container_t::get_pi_overlap_to_ligand_cation(CResidue *res,
							       const clipper::Coord_orth &pt) const {

   float score = 0.0;
   std::string res_name(res->GetResName());
   std::vector<std::vector<std::string> > atom_names = ring_atom_names(res_name);
   for (unsigned int iring=0; iring<atom_names.size(); iring++) {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> residue_pi_points =
	 get_ring_pi_centre_points(atom_names[iring], res);
      float pi_cation_score_1 = overlap_of_cation_pi(pt, residue_pi_points.first);
      float pi_cation_score_2 = overlap_of_cation_pi(pt, residue_pi_points.second );
      if (pi_cation_score_1 > score)
	 score = pi_cation_score_1;
      if (pi_cation_score_2 > score)
	 score = pi_cation_score_2;
   }
   return score;
}


// the sum of the product of the function: m1 exp [ m2 d^2 ]. Where d
// is the is distance of this particular point from either the control
// point (the case of a ring system) or the atom holding the charge (in
// the case of cation in cation-pi interactions).
// 
// m1, m2 pairs from Clarke & Labute:
// 0.78 -1.0
// 7.8  -0.05
// 
float
coot::pi_stacking_container_t::overlap_of_pi_spheres(const clipper::Coord_orth &pt_1,
						     const clipper::Coord_orth &pt_2,
						     const double &m1_pt_1, const double &m2_pt_1,
						     const double &m1_pt_2, const double &m2_pt_2) const {

   // 0.78 exp(-|r-V|^2) where r is a point in space and V is the
   // normally extended ring centre point.  Return the "integral of
   // the function over space".
   
   double score = 0;
   double grid_size = 3.0;
   double grid_step = 0.2;
   double fudge = 5; // so that ideal overlap gives score of 1.0;
   
   for (double delta_x = -grid_size; delta_x<=grid_size; delta_x += grid_step) {
      for (double delta_y = -grid_size; delta_y<=grid_size; delta_y += grid_step) {
	 for (double delta_z = -grid_size; delta_z<=grid_size; delta_z += grid_step) {
	    clipper::Coord_orth t(pt_1.x() + delta_x,
				  pt_1.y() + delta_y,
				  pt_1.z() + delta_z);
	    double d_1_sqd = (pt_1 - t).lengthsq();
	    double d_2_sqd = (pt_2 - t).lengthsq();
	    double score_1 = 0.0;
	    double score_2 = 0.0;
	    if (d_1_sqd < 12.0)
	       score_1 = m1_pt_1 * exp(m2_pt_1 * d_1_sqd);
	    if (d_2_sqd < 12.0)
	       score_2 = m1_pt_2 * exp(m2_pt_2 * d_2_sqd);
	    // std::cout << "scores: " << d_1_sqd << " -> " << score_1 << "     "
	    // << d_2_sqd << " -> " << score_2 << "\n";
	    score += score_1 * score_2;
	 }
      }
   }
   score *= fudge * (grid_step * grid_step * grid_step);
   return score; // double to float
}


// Can throw an exception if not all the ring_atom_names were found in
// res_ref, or lsq plane normal goes bad.
// 
std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::pi_stacking_container_t::get_ring_pi_centre_points(const std::vector<std::string> &ring_atom_names,
							 CResidue *res_ref) const {
   // dummy points, overwritten.
   clipper::Coord_orth pt_1(0,0,0);
   clipper::Coord_orth pt_2(0,0,0);

   int n_residue_atoms;
   PPCAtom residue_atoms = NULL;
   res_ref->GetAtomTable(residue_atoms, n_residue_atoms);

   // fill this vector
   std::vector<clipper::Coord_orth> aromatic_plane_points;
   for (unsigned int iring_at=0; iring_at<ring_atom_names.size(); iring_at++) {
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == ring_atom_names[iring_at]) {
	    clipper::Coord_orth at_pt(residue_atoms[iat]->x,
				      residue_atoms[iat]->y,
				      residue_atoms[iat]->z);
	    aromatic_plane_points.push_back(at_pt);
	    break;
	 }
      }
   }
   if (aromatic_plane_points.size() != ring_atom_names.size()) {
      coot::residue_spec_t spec(res_ref);
      std::string mess = "Not all aromatic atoms were found in residue ";
      mess += spec.chain; 
      mess += " "; 
      mess += spec.resno; 
      throw std::runtime_error(mess);
   } else {
      // can throw an exception
      std::pair<clipper::Coord_orth, clipper::Coord_orth> p =
	 ring_centre_and_normal(aromatic_plane_points);
      double plane_offset = 1; // Angstroms plane normal offset, Clark & Labute 2007
      pt_1 = p.first + plane_offset * p.second;
      pt_2 = p.first - plane_offset * p.second;
      if (0) { 
	 std::cout << "   Centre: " << p.first.x() << " " << p.first.y() << " " << p.first.z() << std::endl;
	 std::cout << "   Normal: " << p.second.x() << " " << p.second.y() << " " << p.second.z() << std::endl;
	 std::cout << " Pi point: " << pt_1.x() << " " << pt_1.y() << " " << pt_1.z() << std::endl;
      }
   }
   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (pt_1, pt_2);
}

// can throw an exception if not enough points found in pts.
//
std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::pi_stacking_container_t::ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const {

   clipper::Coord_orth centre(0,0,0);
   clipper::Coord_orth normal(0,0,0);

   if (pts.size() < 3) {
      std::string mess = "in ring_centre_and_normal() not enough point to calculate lsq plane";
      throw runtime_error(mess);
   }

   double sum_x = 0, sum_y = 0, sum_z = 0;
   for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
      sum_x += pts[ipl].x();
      sum_y += pts[ipl].y();
      sum_z += pts[ipl].z();
   }

   double divisor = 1.0/double(pts.size());
   centre = clipper::Coord_orth(sum_x*divisor, sum_y*divisor, sum_z*divisor);
   clipper::Matrix<double> mat(3,3);
   for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
      mat(0,0) += (pts[ipl].x()-centre.x()) * (pts[ipl].x()-centre.x());
      mat(1,1) += (pts[ipl].y()-centre.y()) * (pts[ipl].y()-centre.y());
      mat(2,2) += (pts[ipl].z()-centre.z()) * (pts[ipl].z()-centre.z());
      mat(0,1) += (pts[ipl].x()-centre.x()) * (pts[ipl].y()-centre.y());
      mat(0,2) += (pts[ipl].x()-centre.x()) * (pts[ipl].z()-centre.z());
      mat(1,2) += (pts[ipl].y()-centre.y()) * (pts[ipl].z()-centre.z());
   }
   mat(1,0) = mat(0,1);
   mat(2,0) = mat(0,2);
   mat(2,1) = mat(1,2);
   std::vector<double> eigens = mat.eigen(false);
   if (0) 
      std::cout << "   we get eigen values: "
		<< eigens[0] << "  "
		<< eigens[1] << "  "
		<< eigens[2] << std::endl;

   int eigen_index = 0;
   if (eigens[1] < eigens[eigen_index])
      eigen_index = 1;
   if (eigens[2] < eigens[eigen_index])
      eigen_index = 2;

   clipper::Coord_orth eigen_vec(mat(0, eigen_index),
				 mat(1, eigen_index),
				 mat(2, eigen_index));

   double sum_sq = 1e-20;
   sum_sq += eigen_vec.x() * eigen_vec.x();
   sum_sq += eigen_vec.y() * eigen_vec.y();
   sum_sq += eigen_vec.z() * eigen_vec.z();

   normal = clipper::Coord_orth(eigen_vec.x()/sum_sq,
				eigen_vec.y()/sum_sq,
				eigen_vec.z()/sum_sq);

   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (centre, normal);
}

// TRP (for example) has 2 rings, so we have to return a vector
//
// Don't forget RNA and DNA! FIXME
//
std::vector<std::vector<std::string> >
coot::pi_stacking_container_t::ring_atom_names(const std::string &residue_name) const {

   std::vector<std::vector<std::string> > v_outer;
   if ((residue_name == "PHE") || (residue_name == "TYR") || (residue_name == "PTY")) {
      std::vector<std::string> v;
      v.push_back(" CG ");
      v.push_back(" CZ ");
      v.push_back(" CD1");
      v.push_back(" CD2");
      v.push_back(" CE1");
      v.push_back(" CE2");
      v_outer.push_back(v);
   }
   if (residue_name == "TRP") {
      std::vector<std::string> v;
      v.push_back(" CG ");
      v.push_back(" CD1");
      v.push_back(" CD2");
      v.push_back(" NE1");
      v.push_back(" CE2");
      v_outer.push_back(v);
      v.clear();
      v.push_back(" CD2");
      v.push_back(" CE2");
      v.push_back(" CE3");
      v.push_back(" CZ3");
      v.push_back(" CZ2");
      v.push_back(" CH2");
      v_outer.push_back(v);
   }

   if ((residue_name == "DA") ||   // or RNA equivalent
       (residue_name == "Ad") ||
       (residue_name == "Ar")) { 
      std::vector<std::string> v;
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N9 ");
      v.push_back(" C8 ");
      v.push_back(" N7 ");
      v_outer.push_back(v);
      v.clear();
      v.push_back(" C6 ");
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N3 ");
      v.push_back(" C2 ");
      v.push_back(" N1 ");
   }
       
   if ((residue_name == "DG") ||   // or RNA equivalent
       (residue_name == "Gd") ||
       (residue_name == "Gr")) { 
      std::vector<std::string> v;
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N9 ");
      v.push_back(" C8 ");
      v.push_back(" N7 ");
      v_outer.push_back(v);
      v.clear();
      v.push_back(" C6 ");
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N3 ");
      v.push_back(" C2 ");
      v.push_back(" N1 ");
   }

   if ((residue_name == "DC") ||   // or RNA equivalent
       (residue_name == "Cd") ||
       (residue_name == "Cr")) { 
      std::vector<std::string> v;
      v.push_back(" N1 ");
      v.push_back(" C2 ");
      v.push_back(" N3 ");
      v.push_back(" C4 ");
      v.push_back(" C5 ");
      v.push_back(" C6 ");
      v_outer.push_back(v);
   }
   
   if ((residue_name == "DT") ||   // or RNA equivalent
       (residue_name == "Td") ||
       (residue_name == "Tr")) { 
      std::vector<std::string> v;
      v.push_back(" N1 ");
      v.push_back(" C2 ");
      v.push_back(" N3 ");
      v.push_back(" C4 ");
      v.push_back(" C5 ");
      v.push_back(" C6 ");
      v_outer.push_back(v);
   }
       
   return v_outer;
}


std::vector<clipper::Coord_orth>
coot::pi_stacking_container_t::get_cation_atom_positions(CResidue *res) const {
   
   std::vector<clipper::Coord_orth> v;

   std::string res_name(res->GetResName());

   if (res_name == "LYS") {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 if (atom_name == " NZ ") {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    v.push_back(pt);
	 }
      }
   }

   if (res_name == "ARG") {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 if ((atom_name == " NH1") ||
	     (atom_name == " NH2")) {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    v.push_back(pt);
	 }
      }
   }
   

   return v;
}

float
coot::pi_stacking_container_t::overlap_of_cation_pi(const clipper::Coord_orth &ligand_pi_point,
						    const clipper::Coord_orth &cation_atom_point) const {

   // multipliers follow the order of the atom/control points
   // 
   float score = overlap_of_pi_spheres(ligand_pi_point, cation_atom_point, 0.78, -1, 7.8, -0.05);

   // std::cout << "overlap_of_cation_pi() calling overlap_of_pi_spheres() returns "
   // << score << std::endl;
   
   return score;
}
 

void
coot::write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			const coot::pi_stacking_container_t &stack_info,
			CResidue *res_flat) {

   if (0) {
      std::cout << "DEBUG:: in write_fle_centres() " << stack_info.stackings.size()
		<< " stackings:\n";
      for (unsigned int i=0; i<stack_info.stackings.size(); i++)
	 std::cout << "    " << i << ":  " << stack_info.stackings[i] << std::endl;
      std::cout << std::endl;
   }
   
   std::string file_name = "coot-tmp-fle-view-residue-info.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      
      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].is_set) { 
	    of << "RES " << v[i].centre.x() << " " << v[i].centre.y() << " "
	       << v[i].centre.z() << " "
	       << v[i].residue_name << " " << v[i].spec.chain << v[i].spec.resno;

	    // Add the qualifying distance to protein if this is a water atom.
	    // 
	    if (v[i].residue_name == "HOH") { 
	       for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
		  if (bonds_to_ligand[ib].res_spec == v[i].spec) {
		     of << " water->protein-length " << bonds_to_ligand[ib].water_protein_length;
		     break;
		  }
	       }
	       of << "\n";
	    } else { 
	       of << "\n";
	    }


	    // now, was there a bond with this res spec?
	    // 
	    // if so find the x y coordinates of the atom with this
	    // name in the res_flat and write them out.
	    // 
	    for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	       if (bonds_to_ligand[ib].res_spec == v[i].spec) {
		  PPCAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_flat->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
		     CAtom *at = residue_atoms[iat];
		     std::string atom_name(at->name);
		     if (atom_name == bonds_to_ligand[ib].ligand_atom_spec.atom_name) {
			// the coordinates of the matching flat ligand
			// atom (done like solvent accessibilites).
			of << "BOND " << at->name
			   << " length: " << bonds_to_ligand[ib].bond_length
			   << " type: " << bonds_to_ligand[ib].bond_type
			   << "\n";
		     }
		  }
	       }
	    }


	    // Is there a solvent accessibility associated with this residue?
	    // 
	    for (unsigned int ised=0; ised<sed.size(); ised++) { 
	       if (sed[ised].res_spec == v[i].spec) {
		  of << "SOLVENT_ACCESSIBILITY_DIFFERENCE "
		     << sed[ised].exposure_fraction_holo << " " 
		     << sed[ised].exposure_fraction_apo << "\n";
	       }
	    }

	    // Is there a stacking interaction associated with this residue?
	    //
	    // There could, in fact, be several stackings associated
	    // with this residue.  We want to pick the best one.  So
	    // we need to run through the list of stacking first,
	    // looking for the best.
	    //
	    // 
	    double best_stacking_score = -1.0;
	    coot::pi_stacking_instance_t best_stacking(NULL, "");
	    
	    for (unsigned int istack=0; istack<stack_info.stackings.size(); istack++) {
	       coot::residue_spec_t spec(stack_info.stackings[istack].res);
	       if (spec == v[i].spec) {
		  if (stack_info.stackings[istack].overlap_score > best_stacking_score) {
		     best_stacking_score = stack_info.stackings[istack].overlap_score;
		     best_stacking = stack_info.stackings[istack];
		  }
	       }
	    }

	    if (best_stacking_score > 0.0) {
	       std::string type = "    pi-pi"; // leading spaces
					       // important in format
					       // of atom names

	       // Recall that CATION_PI_STACKING is for ligand ring
	       // systems, PI_CATION_STACKING is for protein residue
	       // ring systems.
	       // 
	       if (best_stacking.type == coot::pi_stacking_instance_t::CATION_PI_STACKING)
		  type = "cation-pi";
	       if (best_stacking.type == coot::pi_stacking_instance_t::PI_CATION_STACKING)
		  type = "pi-cation";
	       if (0) // debug
		  of << "# Best stacking for RES " 
		     << v[i].residue_name << " " << v[i].spec.chain
		     << v[i].spec.resno << " is " << best_stacking.type << " "
		     << best_stacking_score << std::endl;
	       of << "STACKING " << type << " ";
	       switch (best_stacking.type) {
	       case coot::pi_stacking_instance_t::NO_STACKING:
		  break;
	       case coot::pi_stacking_instance_t::PI_PI_STACKING: 
		  for (unsigned int jat=0;
		       jat<best_stacking.ligand_ring_atom_names.size();
		       jat++) {
		     of << best_stacking.ligand_ring_atom_names[jat] << "  ";
		  }
		  break;
		  
	       case coot::pi_stacking_instance_t::CATION_PI_STACKING:
		  of << best_stacking.ligand_cationic_atom_name;
		  break;
		  
	       case coot::pi_stacking_instance_t::PI_CATION_STACKING:
		  for (unsigned int jat=0;
		       jat<best_stacking.ligand_ring_atom_names.size();
		       jat++) {
		     of << best_stacking.ligand_ring_atom_names[jat] << "  ";
		  }
		  break;
	       }
	       of << "\n";
	    }
	    
	 } else {
	    std::cout << "WARNING centre " << i << " was not set " << std::endl;
	 } 
      }
   }
}

std::ostream&
coot::operator<< (std::ostream& s, const coot::bash_distance_t &bd) {

   if (bd.unlimited()) { 
      s << "unlimited"; // C&L.
   } else {
      s << bd.dist;
   } 
   return s;
} 


// the reference_residue is the flat residue
void
coot::write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					const coot::flev_attached_hydrogens_t &attached_hydrogens,
					CResidue *reference_residue) {

   // for each of the atoms in reference_residue, find the spec in sav
   // and write out the coordinates and the atom name and the
   // solvent_accessibility.

   std::string file_name = "coot-tmp-fle-view-solvent-accessibilites.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      reference_residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int i=0; i<n_residue_atoms; i++) {
	 CAtom *at = residue_atoms[i];
	 std::string atom_name(at->GetAtomName());
	 for (unsigned int j=0; j<sav.size(); j++) {
	    if (atom_name == sav[j].first.atom_name) {
	       of << "ATOM:" << atom_name << " " << at->x << " " <<  at->y << " " << at->z
		  << " " << sav[j].second << "\n";
	       std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it =
		  attached_hydrogens.atom_bashes.find(atom_name);
	       if (it != attached_hydrogens.atom_bashes.end()) {
		  for (unsigned int ibash=0; ibash<it->second.size(); ibash++) { 
		     of << "   BASH: " << it->second[ibash] << "\n";
		  }
	       } 
	       break;
	    }
	 }
      }
   }
}


// Use coot::h_bonds class to generate ligands.  We do that by creating a synthetic
// temporary  molecule and atom selections.
// 
std::vector<coot::fle_ligand_bond_t>
coot::get_fle_ligand_bonds(CResidue *ligand_res,
			   const std::vector<CResidue *> &residues,
			   CMMDBManager *mol, 
			   const std::map<std::string, std::string> &name_map,
			   const coot::protein_geometry &geom,
			   float water_dist_max) {
   
   std::vector<coot::fle_ligand_bond_t> v; // returned value

   bool debug = 0;
   std::vector<CResidue *> rv = residues;
   rv.push_back(ligand_res);

   std::pair<bool, CMMDBManager *> m = coot::util::create_mmdbmanager_from_residue_vector(rv);
   coot::residue_spec_t ligand_spec(ligand_res);

   if (m.first) { 
      int SelHnd_all = m.second->NewSelection();
      int SelHnd_lig = m.second->NewSelection();
      m.second->SelectAtoms(SelHnd_all, 0, "*", ANY_RES, "*", ANY_RES, "*", "*", "*", "*", "*");
      m.second->SelectAtoms(SelHnd_lig, 0, ligand_spec.chain.c_str(),
			    ligand_spec.resno, ligand_spec.insertion_code.c_str(),
			    ligand_spec.resno, ligand_spec.insertion_code.c_str(),
			    "*", "*", "*", "*");

      // -----------------------
      //   hydrogen bonds 
      // -----------------------
      
      coot::h_bonds hb;
      // std::vector<coot::h_bond> hbonds = hb.get(SelHnd_lig, SelHnd_all, m.second, geom);
      std::vector<coot::h_bond> hbonds = hb.get_mcdonald_and_thornton(SelHnd_lig, SelHnd_all, m.second, geom);

      if (debug)
	 std::cout << "DEBUG:: get_fle_ligand_bonds from h_bonds class found "
		   << hbonds.size() << " H bonds." << std::endl;

      for (unsigned int i=0; i<hbonds.size(); i++) { 
	 std::cout << coot::atom_spec_t(hbonds[i].donor) << "..."
		   << coot::atom_spec_t(hbonds[i].acceptor) << " with ligand donor flag "
		   << hbonds[i].ligand_atom_is_donor << std::endl;

	 int bond_type = fle_ligand_bond_t::get_bond_type(hbonds[i].donor,
							  hbonds[i].acceptor,
							  hbonds[i].ligand_atom_is_donor);
	 // override these 2 if ligand atom is donor
	 // 
	 CAtom      *ligand_atom = hbonds[i].acceptor;
	 CAtom *env_residue_atom = hbonds[i].donor;
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
		      << env_residue_atom << std::endl;

	 // we want to pass the atom specifics (not just the environ residue spec)
	 // 
	 // coot::fle_ligand_bond_t bond(ligand_atom_name, bond_type, hbonds[i].dist, res_spec);
	 // 
	 coot::fle_ligand_bond_t bond(coot::atom_spec_t(ligand_atom),
				      coot::atom_spec_t(env_residue_atom), 
				      bond_type,
				      hbonds[i].dist+explict_H_bond_fudge_factor);
	 
	 std::string residue_name = ligand_atom->GetResName();
	 if (residue_name == "HOH")
	    bond.water_protein_length = find_water_protein_length(ligand_atom->residue, mol);
	 
	 bool ok_to_add = 1; // reset to 0 for certain waters
	 if (ligand_atom) {
	    std::string ligand_residue_name(ligand_atom->GetResName());
	    if (ligand_residue_name == "HOH") {
	       if (hbonds[i].dist > water_dist_max)
		  ok_to_add = 0;
	    }
	 }
	 if (ok_to_add)
	    v.push_back(bond);
      }

      std::cout << ".... get_fle_ligand_bonds(): after h-bonds v.size() is " << v.size() << std::endl;

      // -----------------------
      //   covalent bonds 
      // -----------------------

      std::vector<coot::fle_ligand_bond_t> covalent_bonds =
	 get_covalent_bonds(m.second, SelHnd_lig, SelHnd_all, ligand_spec, geom);
      for (unsigned int i=0; i<covalent_bonds.size(); i++)
	 v.push_back(covalent_bonds[i]);

      std::cout << ".... get_fle_ligand_bonds(): after covalent bonds v.size() is " << v.size()
		<< std::endl;

      // -----------------------
      //   metal bonds 
      // -----------------------

      std::vector<coot::fle_ligand_bond_t> metal_bonds = get_metal_bonds(ligand_res, residues);
      for (unsigned int i=0; i<metal_bonds.size(); i++)
	 v.push_back(metal_bonds[i]);

      std::cout << ".... get_fle_ligand_bonds(): after metal bonds v.size() is " << v.size()
		<< std::endl;

      
      // -----------------------
      //   clean up 
      // -----------------------
      
      m.second->DeleteSelection(SelHnd_lig);
      m.second->DeleteSelection(SelHnd_all);
      delete m.second;

   }

   if (1) {
      std::cout << ":::: get_fle_ligand_bonds returns these " << v.size()
		<< " bonds: " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) { 
	 std::cout << "   " << i << "  " << v[i] << std::endl;
      }
   } 
   return v;
}

std::vector<coot::fle_ligand_bond_t>
coot::get_covalent_bonds(CMMDBManager *mol,
			 int SelHnd_lig,
			 int SelHnd_all,
			 const residue_spec_t &ligand_spec,
			 const protein_geometry &geom) {

   // 20101016, hydrogens don't make covalent bonds between ligands
   // and protein (or residues to residues in general).
   
   std::vector<coot::fle_ligand_bond_t> v;
   int SelHnd_local = mol->NewSelection();
   mol->SelectAtoms(SelHnd_local, 0, "*", ANY_RES, "*", ANY_RES, "*", "*", "*", "*", "*");
   mol->Select(SelHnd_local, STYPE_ATOM, 0, ligand_spec.chain.c_str(),
	       ligand_spec.resno, ligand_spec.insertion_code.c_str(),
	       ligand_spec.resno, ligand_spec.insertion_code.c_str(),
	       "*", "*", "*", "*", SKEY_XOR);

   // now find contacts:
   // 
   PSContact pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mat44 my_matt;
   CSymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   PPCAtom lig_atom_selection = 0;
   int n_lig_atoms;
   mol->GetSelIndex(SelHnd_lig, lig_atom_selection, n_lig_atoms);

   PPCAtom other_atom_selection = 0;
   int n_other_atoms;
   mol->GetSelIndex(SelHnd_local, other_atom_selection, n_other_atoms);

   realtype min_dist = 0.1;
   realtype max_dist = 2.3; //  even S-S is shorter than this, I think
   realtype cno_max_dist = 1.8;

   mol->SeekContacts(  lig_atom_selection,   n_lig_atoms,
		     other_atom_selection, n_other_atoms,
		     min_dist, max_dist, // min, max distances
		     0,        // seqDist 0 -> in same res also
		     pscontact, n_contacts,
		       0, &my_matt, i_contact_group);

   std::vector<std::pair<CResidue *, CResidue *> > contacting_pairs_vec;

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i=0; i<n_contacts; i++) {
	    CAtom *at_1 =   lig_atom_selection[pscontact[i].id1];
	    CAtom *at_2 = other_atom_selection[pscontact[i].id2];
// 	    std::cout << "DEUBG:: Covalent test "
// 		      << coot::atom_spec_t(at_1) << "..."
// 		      << coot::atom_spec_t(at_2) << std::endl;
	    std::pair<CResidue *, CResidue *> pair(at_1->GetResidue(),
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
			int bond_type = coot::fle_ligand_bond_t::BOND_COVALENT;
			coot::fle_ligand_bond_t bond(coot::atom_spec_t(at_1), // ligand
						     coot::atom_spec_t(at_2), // env residue
						     bond_type, d);
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

std::vector<coot::fle_ligand_bond_t>
coot::get_metal_bonds(CResidue *ligand_residue, const std::vector<CResidue *> &residues) {

   // a non-Hydrogen, non-Carbon ligand atom, that is.
   double max_dist_metal_to_ligand_atom = 3.5; // pass this parameter?
   
   std::vector<coot::fle_ligand_bond_t> v;

   double best_dist_sqrd = max_dist_metal_to_ligand_atom * max_dist_metal_to_ligand_atom;
   std::string ligand_atom_name; // goes with best_dist_sqrd
   CAtom *ligand_atom = NULL;
   CAtom *env_residue_atom = NULL;
   
   PPCAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_residue->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);
   for (unsigned int i=0; i<residues.size(); i++) { 
      if (coot::is_a_metal(residues[i])) {
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 residues[i]->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (unsigned int irat=0; irat<n_residue_atoms; irat++) {
	    for (unsigned int ilat=0; ilat<n_residue_atoms; ilat++) {
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
	    coot::fle_ligand_bond_t bond(ligand_atom,
					 env_residue_atom,
 					 coot::fle_ligand_bond_t::METAL_CONTACT_BOND,
 					 sqrt(best_dist_sqrd));
	    v.push_back(bond);
	 }
      }
   }

   return v;
}


// should be in coot-utils perhaps?
//
bool
coot::is_a_metal(CResidue *res) {

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


// return 100 if no other contact found (strange!)
// 
double
coot::find_water_protein_length(CResidue *ligand_residue, CMMDBManager *mol) {

   double dist = 100;

   double dist_sqrd = dist * dist;
   double dist_sqrd_init = dist_sqrd;
   
   PPCAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_residue->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);


   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 if (ligand_residue != residue_p) {
	    std::string residue_name(residue_p->GetResName());
	    if (residue_name != "HOH") { 
	       PPCAtom residue_atoms = 0;
	       int n_residue_atoms;
	       residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (unsigned int il=0; il<n_ligand_residue_atoms; il++) {
		  for (unsigned int irat=0; irat<n_residue_atoms; irat++) {
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

   

bool
coot::standard_residue_name_p(const std::string &rn) {

   std::vector<std::string> s;

   s.push_back("ALA");
   s.push_back("ASP");
   s.push_back("ASN");
   s.push_back("CYS");
   s.push_back("GLN");
   s.push_back("GLY");
   s.push_back("GLU");
   s.push_back("PHE");
   s.push_back("HIS");
   s.push_back("ILE");
   s.push_back("LYS");
   s.push_back("LEU");
   s.push_back("MET");
   s.push_back("MSE");
   s.push_back("PRO");
   s.push_back("ARG");
   s.push_back("SER");
   s.push_back("THR");
   s.push_back("VAL");
   s.push_back("TRP");
   s.push_back("TYR");

   s.push_back("AR");
   s.push_back("AD");
   s.push_back("CR");
   s.push_back("CD");
   s.push_back("GR");
   s.push_back("GD");
   s.push_back("TD");

   if (std::find(s.begin(), s.end(), rn) == s.end())
      return 0;
   else
      return 1;
}

// The substitution contour:
// 
// The substitution contour depends on the hydridisation of the atoms
// of the ligands.  We need to "fire cannonballs" along the direction
// of the hydrogens and see when they bash into the atoms of the
// protein.
//
// Most (I guess) hydrogens of ligands are "riding", by that, I mean
// that the position of the hydrogen is determined by the geometry of
// the atom to which it is attached (e.g. the H on an N that has 2
// single bonds to 2 carbons).
//
// So, PRODRG will give us a ligand that has hydrogens on it - good
// show (we will need to run it in "regular" mode ("MINI PREP", say)
// to generate such hydrogens on the model).
//
// We need to determine which are the riding hydrogens and which are
// the rotatable ones. The riding hydrogens give us the cannonball
// vector trivially. For the rotable ones, we should draw a number
// from a probability distribution defined by the parameters in the
// cif dictionary.
// 
// How do we know which are riding and which are not?
//
// Here is my guess: riding hydrogens do not have torsions for them in
// the dictionary.  Rotable hydrogens do.  Note though, that PRODRG
// hydrogen names (in the cif file) are problematic.  You can't rely
// on name matching of hydrogens between the PDB and the dictionary.
//
// So here is what we have to do:
//
// i) From the dictionary, look through the list of non-Hydrogen atoms
// to see if there are hydrogens bonded to the atom.
// 
//    ii) If so, does this non-Hydrogen atom appear as atom 2 in a 
//    torsion restraint with a hydrogen as atom 1 or appear as atom 3
//    in a torsion restraint with a hydrogen as atom 4?
//
//        iii) If so, then this is a torsionable hydrogen, we need to
//        generate vectors by random sampling from the probability
//        distribution of this torsion.
//
//        However, if it is not, then this is a riding hydrogen.  Add
//        the vector from the ligand atom to the hydrogen as a
//        cannonball vector associated with that ligand atom.
//
// OK, at the end of the hybridisation analysis, we want a vector of
// ligand "heavy" (in this case, non-Hydrogen) atoms and associated
// with each atom is a vector (i.e. list) of cannonball vectors
// (i.e. directions).  It is not yet clear to me if we need anything
// else, so maybe a vector of pairs of ligand atoms and a vector of
// Coord_orths will do.
//
// Fine, but note that we get the hydrogen directions from the prodrg
// ligand, but we want the hydrogen position on the atoms of our
// (reference) ligand.  So, we'll have to do a name mapping of the
// non-hydrogen atoms using mmdb graph match.  Complicated.
//
// This vector of pairs will in turn be processed into a "bash"
// distance for each atom (how far can we go from the atom before
// bashing in the atoms of the protein?)  Note that, in doing so, we
// need to be able to specify "infinite bash distance". Let's make a
// trivial class for bash distance.
//
// Then these distances get overlayed onto the ligand grid.
//
// The the ligand grid gets contoured.
//
// And plotted to the canvas.

std::vector<std::pair<std::string, int> >
coot::get_prodrg_hybridizations(const std::string &prodrg_log_file_name) {

   // This is currently not used. Now we get the cannonball directions
   // from the positions of the hydrogens directly, not infered from
   // the hybridisation of the atoms.

   std::cout << "in get_prodrg_hybridizations() "  << std::endl;
   std::vector<std::pair<std::string, int> > v;
   std::ifstream f(prodrg_log_file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: failed to open " << prodrg_log_file_name << std::endl;
   } else {
      std::string line;
      bool hybridizations_flag = 0;
      while (std::getline(f, line)) {
	 std::vector<std::string> words =
	    coot::util::split_string_no_blanks(line, " ");
	 if (words.size() > 1) {
	    if (words[0] == "DETERMINING") {
	       if (words[1] == "HYBRIDISATION") {
		  hybridizations_flag = 1;
	       }
	    }
	    if (words[0] == "PRODRG>") {
	       hybridizations_flag = 0;
	    }
	 }
	 if (hybridizations_flag) {
	    if (line.find("qual = ") != std::string::npos ) {
	       std::string atom_name=line.substr(1,4);
	       std::string hybridization_str = line.substr(13,3);
	       std::cout << ":" << atom_name <<": :" << hybridization_str
			 << ":" << std::endl;
	       int hybridization = coot::SP_HYBRIDIZATION;
	       if (hybridization_str == "sp3")
		  hybridization = coot::SP3_HYBRIDIZATION;
	       if (hybridization_str == "sp2")
		  hybridization = coot::SP2_HYBRIDIZATION;
	       std::pair<std::string, int> p(atom_name, hybridization);
	       v.push_back(p);
	    }
	 }
      }
   }
   return v;
}

std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > >
coot::get_cannonball_vectors(CResidue *ligand_res_3d,
			     const coot::dictionary_residue_restraints_t &monomer_restraints) {

   std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > > v;

   return v;
}

// examine the dictionary and find the atoms to which the hydrogens
// are attached.  Is the hydrogen riding or rotatable?
// 
//            ... riding hydrogens do not have torsions for them in
// the dictionary.  Rotable hydrogens do.  Note though, that PRODRG
// hydrogen names (in the cif file) are problematic.  You can't rely
// on name matching of hydrogens between the PDB and the dictionary.
//
// So here is what we have to do:
//
// i) From the dictionary, look through the list of non-Hydrogen atoms
// to see if there are hydrogens bonded to the atom.
// 
//    ii) If so, does this non-Hydrogen atom appear as atom 2 in a 
//    torsion restraint with a hydrogen as atom 1 or appear as atom 3
//    in a torsion restraint with a hydrogen as atom 4?
//
//        iii) If so, then this is a torsionable hydrogen, we need to
//        generate vectors by random sampling from the probability
//        distribution of this torsion.
//
//        However, if it is not, then this is a riding hydrogen.  Add
//        the vector from the ligand atom to the hydrogen as a
//        cannonball vector associated with that ligand atom.
//
coot::flev_attached_hydrogens_t::flev_attached_hydrogens_t(const coot::dictionary_residue_restraints_t &restraints) {

   for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
      std::string atom_name_1 = restraints.bond_restraint[ibond].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ibond].atom_id_2_4c();
      if ((restraints.is_hydrogen(atom_name_1)) && (! restraints.is_hydrogen(atom_name_2))) {
	 std::swap(atom_name_1, atom_name_2);
      }
      
      if ((! restraints.is_hydrogen(atom_name_1)) && (restraints.is_hydrogen(atom_name_2))) {
	 // a heavy atom connected to a hydrogen.
	 // Does it exist in the torsions?

	 bool found = 0;
	 std::pair<std::string, std::string> p(atom_name_1, atom_name_2);
	 for (unsigned int itor=0; itor<restraints.torsion_restraint.size(); itor++) { 
	    if (((restraints.torsion_restraint[itor].atom_id_1_4c() == atom_name_2) &&
		 (restraints.torsion_restraint[itor].atom_id_2_4c() == atom_name_1)) ||
		((restraints.torsion_restraint[itor].atom_id_4_4c() == atom_name_2) &&
		 (restraints.torsion_restraint[itor].atom_id_3_4c() == atom_name_1))) { 
	       if (restraints.torsion_restraint[itor].is_const()) {
		  atoms_with_riding_hydrogens.push_back(p);
	       } else {
		  atoms_with_rotating_hydrogens.push_back(p);
	       }
	       found = 1;
	       break;
	    }
	 }
	 if (! found) {
	    atoms_with_riding_hydrogens.push_back(p);
	 }
      }
   }
} 

void
coot::flev_attached_hydrogens_t::cannonballs(CResidue *ligand_residue_3d,
					     const std::string &prodrg_3d_ligand_file_name,
					     const coot::dictionary_residue_restraints_t &restraints) {

   atom_selection_container_t asc = get_atom_selection(prodrg_3d_ligand_file_name, 1);
   if (asc.read_success) {
      cannonballs(ligand_residue_3d, asc.mol, restraints);
   }
}

void
coot::flev_attached_hydrogens_t::cannonballs(CResidue *ligand_residue_3d,
					     CMMDBManager *mol, 
					     const coot::dictionary_residue_restraints_t &restraints) {

   if (! mol)
      return;
   
   PSContact pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mat44 my_matt;
   CSymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;


   int SelHnd_H = mol->NewSelection();
   int SelHnd_non_H = mol->NewSelection();

   PPCAtom hydrogen_selection = 0;
   PPCAtom non_hydrogen_selection = 0;
   int n_hydrogen_atoms;
   int n_non_hydrogen_atoms;
      
      
   mol->SelectAtoms(SelHnd_H,     0, "*", ANY_RES, "*", ANY_RES, "*", "*", "*", " H", "*");
   mol->SelectAtoms(SelHnd_non_H, 0, "*", ANY_RES, "*", ANY_RES, "*", "*", "*", "!H", "*");
      
   mol->GetSelIndex(SelHnd_H, hydrogen_selection, n_hydrogen_atoms);
   mol->GetSelIndex(SelHnd_non_H, non_hydrogen_selection, n_non_hydrogen_atoms);
      
   std::cout << "Found " << n_hydrogen_atoms << " Hydrogens " << std::endl;
   std::cout << "Found " << n_non_hydrogen_atoms << " non Hydrogens " << std::endl;

   if (n_hydrogen_atoms == 0) {
      std::cout << "WARNING:: Oops found no hydrogens for cannonballs" << std::endl;
      return;
   }
   if (n_non_hydrogen_atoms == 0) {
      std::cout << "WARNING:: Oops found no non-hydrogens for cannonballs" << std::endl;
      return;
   }

   mol->SeekContacts(hydrogen_selection, n_hydrogen_atoms,
		     non_hydrogen_selection, n_non_hydrogen_atoms,
		     0.1, 1.5,
		     0, // in same res also
		     pscontact, n_contacts,
		     0, &my_matt, i_contact_group);

   std::cout << "Found " << n_contacts << " contacts to Hydrogens " << std::endl;

   // We need to find the torsion of the hydrogen,
   // A torsion (that can be mapped to the reference ligand) is:
   //
   // At_name_base At_name_2 At_name_bond_to_H bond_length bond_angle torsion_angle
   //
   // that is, we work from "inside" ligand atoms out to the hydrogen
   // 
   // 
   if (n_contacts > 0) {
      for (int i=0; i< n_contacts; i++) {
	 CAtom *at = non_hydrogen_selection[pscontact[i].id2];
	 std::string atom_name_bonded_to_H(at->name);

	 bool found_torsion_for_this_H = 0;

	 // riding hydrogens:
	 // 
	 for (unsigned int iat=0; iat<atoms_with_riding_hydrogens.size(); iat++) { 
	    if (atom_name_bonded_to_H == atoms_with_riding_hydrogens[iat].first) {
	       CAtom *h_at = hydrogen_selection[pscontact[i].id1];
	       found_torsion_for_this_H = add_named_torsion(h_at, at, restraints, mol, coot::H_IS_RIDING);
	    }
	    if (found_torsion_for_this_H)
	       break;
	 }

	 // rotating hydrogens:
	 // 
	 for (unsigned int iat=0; iat<atoms_with_rotating_hydrogens.size(); iat++) { 
	    if (atom_name_bonded_to_H == atoms_with_rotating_hydrogens[iat].first) {
	       CAtom *h_at = hydrogen_selection[pscontact[i].id1];
	       found_torsion_for_this_H = add_named_torsion(h_at, at, restraints, mol, coot::H_IS_ROTATABLE);
	    }
	    if (found_torsion_for_this_H)
	       break;
	 }
      }
   }

   mol->DeleteSelection(SelHnd_H);
   mol->DeleteSelection(SelHnd_non_H);

   named_hydrogens_to_reference_ligand(ligand_residue_3d, restraints);

}

// hydrogen_type is either H_IS_RIDING or H_IS_ROTATABLE
//
bool
coot::flev_attached_hydrogens_t::add_named_torsion(CAtom *h_at, CAtom *at,
						   const coot::dictionary_residue_restraints_t &restraints,
						   CMMDBManager *mol, // 3d prodrg ligand mol
						   int hydrogen_type)  {

   bool found_torsion_for_this_H = 0;
   std::string atom_name_bonded_to_H(at->name);
   clipper::Coord_orth p_h(h_at->x, h_at->y, h_at->z);
   clipper::Coord_orth p_1(at->x, at->y, at->z);

   // now we work back through the restraints, finding
   // an atom that bonds to at/atom_name_bonded_to_H,
   // and one that bonds to that (by name).
   for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
      std::string atom_name_1 = restraints.bond_restraint[ibond].atom_id_1_4c();
      std::string atom_name_2 = restraints.bond_restraint[ibond].atom_id_2_4c();

      if (atom_name_bonded_to_H == atom_name_2)
	 std::swap(atom_name_1, atom_name_2);

      if (atom_name_bonded_to_H == atom_name_1) {

	 if (! restraints.is_hydrogen(atom_name_2)) { 
	    std::string At_name_2 = atom_name_2;

	    // now for the base...
	    //
	    for (unsigned int jbond=0; jbond<restraints.bond_restraint.size(); jbond++) {
	       std::string atom_name_b_1 = restraints.bond_restraint[jbond].atom_id_1_4c();
	       std::string atom_name_b_2 = restraints.bond_restraint[jbond].atom_id_2_4c();

	       if (At_name_2 == atom_name_b_2)
		  std::swap(atom_name_b_1, atom_name_b_2); // same trick

	       if (At_name_2 == atom_name_b_1) {

		  // atom_name_b_2 is the base then (maybe)

		  if (atom_name_b_2 != atom_name_bonded_to_H) {
		     if (! restraints.is_hydrogen(atom_name_b_1)) {
			std::string base_atom_name = atom_name_b_2;

			// now, where are those atoms (At_name_2 and base_atom_name)?
			CAtom *At_2 = NULL;
			CAtom *base_atom = NULL;

			int imod = 1;
			CModel *model_p = mol->GetModel(imod);
			CChain *chain_p;
			int nchains = model_p->GetNumberOfChains();
			for (int ichain=0; ichain<nchains; ichain++) {
			   chain_p = model_p->GetChain(ichain);
			   int nres = chain_p->GetNumberOfResidues();
			   CResidue *residue_p;
			   CAtom *residue_at;
			   for (int ires=0; ires<nres; ires++) { 
			      residue_p = chain_p->GetResidue(ires);
			      int n_atoms = residue_p->GetNumberOfAtoms();
			      for (int iat=0; iat<n_atoms; iat++) {
				 residue_at = residue_p->GetAtom(iat);
				 std::string res_atom_name(residue_at->name);
				 if (res_atom_name == At_name_2)
				    At_2 = residue_at;
				 if (res_atom_name == base_atom_name)
				    base_atom = residue_at;
			      }
			   }
			}

			if (!base_atom || !At_2) {

			   if (!base_atom)
			      std::cout << "Failed to find base in 3d prodrg residue "
					<< base_atom_name << std::endl;
			   if (!At_2)
			      std::cout << "Failed to find base or At_2 in 3d prodrg residue "
					<< At_name_2 << std::endl;
			} else {
			   try { 
			      clipper::Coord_orth p_2(At_2->x, At_2->y, At_2->z);
			      clipper::Coord_orth p_base(base_atom->x, base_atom->y, base_atom->z);
					     
			      double tors_r = clipper::Coord_orth::torsion(p_base, p_2, p_1, p_h);
			      double tors = clipper::Util::rad2d(tors_r);
			      double angle = coot::angle(h_at, at, At_2);
			      double dist = clipper::Coord_orth::length(p_h, p_1);

			      coot::named_torsion_t torsion(base_atom_name,
							    At_name_2,
							    atom_name_bonded_to_H,
							    dist, angle, tors, hydrogen_type);

			      // std::cout << "  Yeah!!! adding named torsion " << std::endl;
			      named_torsions.push_back(torsion);
			      found_torsion_for_this_H = 1;
			   }
			   catch (std::runtime_error rte) {
			      std::cout << "WARNING:: " << rte.what() << std::endl;
			   } 
			} 
		     }
		  }
	       }
	       if (found_torsion_for_this_H)
		  break;
	    }
	 }
      }
      if (found_torsion_for_this_H)
	 break;
   }

   return found_torsion_for_this_H;
}

// For (each?) of the atoms in our real reference residue
// ligand_residue_3d (that should have hydrogens attached) give us a
// unit vector from the bonding atom in the direction of (each of, if
// there are more than one) the hydrogen(s).
// 
std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > >
coot::flev_attached_hydrogens_t::named_hydrogens_to_reference_ligand(CResidue *ligand_residue_3d,
								     const coot::dictionary_residue_restraints_t &restraints) const {

   std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > > v;

   for (unsigned int i=0; i<named_torsions.size(); i++) { 
      if (named_torsions[i].hydrogen_type == coot::H_IS_RIDING) {
	 CAtom *atom_base = NULL;
	 CAtom *atom_2 = NULL;
	 CAtom *atom_bonded_to_H = NULL;
	 
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 ligand_residue_3d->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) { 
	    std::string atom_name(residue_atoms[iat]->name);
	    if (atom_name == named_torsions[i].base_atom_name) {
	       atom_base = residue_atoms[iat];
	    }
	    if (atom_name == named_torsions[i].atom_name_2) {
	       atom_2 = residue_atoms[iat];
	    }
	    if (atom_name == named_torsions[i].atom_name_bonded_to_H) {
	       atom_bonded_to_H = residue_atoms[iat];
	    }
	 }

	 if (atom_base && atom_2 && atom_bonded_to_H) {
	    clipper::Coord_orth pos_atom_base(atom_base->x, atom_base->y, atom_base->z);
	    clipper::Coord_orth pos_atom_2(atom_2->x, atom_2->y, atom_2->z);
	    clipper::Coord_orth pos_atom_bonded_to_H(atom_bonded_to_H->x, atom_bonded_to_H->y, atom_bonded_to_H->z);

	    clipper::Coord_orth new_pt(pos_atom_base, pos_atom_2, pos_atom_bonded_to_H,
				       1.0, // unit vector
				       clipper::Util::d2rad(named_torsions[i].angle),
				       clipper::Util::d2rad(named_torsions[i].torsion));
	    clipper::Coord_orth vect = new_pt = pos_atom_bonded_to_H;


	    // add that to the pile
	    bool found_atom = 0;
	    for (unsigned int irv=0; irv<v.size(); irv++) { 
	       if (v[irv].first == atom_bonded_to_H) {
		  v[irv].second.push_back(vect);
		  found_atom = 1;
	       }
	    }
	    if (! found_atom) {
	       std::vector<clipper::Coord_orth> cov;
	       cov.push_back(vect);
	       std::pair<CAtom *, std::vector<clipper::Coord_orth> > p(atom_bonded_to_H, cov);
	       v.push_back(p);
	    } 
	 }
      }
   }

   return v;
}
								     

// apply those cannonball directions onto the real reference ligand:
void
coot::flev_attached_hydrogens_t::distances_to_protein(CResidue *residue_reference, 
						      CMMDBManager *mol_reference) {

   float radius = 6.0;
   std::vector<CResidue *> env_residues =
      coot::residues_near_residue(residue_reference, mol_reference, radius);
   bool debug = 0;

   if (debug) { 
      std::cout << "named torsion: " << named_torsions.size() << std::endl;
      for (unsigned int i=0; i<named_torsions.size(); i++) {
	 std::cout << "  " << i
		   << named_torsions[i].base_atom_name << "  "
		   << named_torsions[i].atom_name_2 << "  "
		   << named_torsions[i].atom_name_bonded_to_H << "  "
		   << named_torsions[i].dist << "  "
		   << named_torsions[i].angle << "  "
		   << named_torsions[i].torsion << "  "
		   << std::endl;
      }
   }
   
   for (unsigned int i=0; i<named_torsions.size(); i++) {
      try {
	 if (named_torsions[i].hydrogen_type == coot::H_IS_RIDING) {

	    std::cout << "hydrogen on " << named_torsions[i].atom_name_bonded_to_H 
		      << " is riding " << std::endl;
	    
	    std::pair<clipper::Coord_orth, clipper::Coord_orth> pt_base_and_H =
	       hydrogen_pos(named_torsions[i], residue_reference);

	    std::cout << "pt_base_and_H: ligand atom at: " << pt_base_and_H.first.format()
		      << " H atom at: " << pt_base_and_H.second.format() << std::endl;
	    
	    
	    std::vector<CAtom *> atoms = close_atoms(pt_base_and_H.second, env_residues);

	    // bash is one of a (potentially) number of bash distances
	    // for a given ligand (non-hydrogen) atom.
	    // (passing the H coord, the lig-atom coord, a vector of CAtoms *s.)
	    // 
	    coot::bash_distance_t bash = find_bash_distance(pt_base_and_H.first,
							    pt_base_and_H.second,
							    atoms);
	    std::cout << "   found bash: " << bash << std::endl;
	    atom_bashes[named_torsions[i].atom_name_bonded_to_H].push_back(bash);
	 }
	 
	 if (named_torsions[i].hydrogen_type == coot::H_IS_ROTATABLE) {

	    std::cout << "hydrogen on " << named_torsions[i].atom_name_bonded_to_H 
		      << " is rotatable " << std::endl;
	    
	    int n_torsion_samples = 8;
	    std::pair<clipper::Coord_orth, clipper::Coord_orth> pt_base_and_H =
	       hydrogen_pos(named_torsions[i], residue_reference);
	    std::vector<CAtom *> atoms = close_atoms(pt_base_and_H.second, env_residues);
	    for (unsigned int itor=0; itor<n_torsion_samples; itor++) {

	       coot::named_torsion_t tmp_tor = named_torsions[i];
	       tmp_tor.torsion += double(itor) * 360.0/double(n_torsion_samples);
	       if (tmp_tor.torsion > 360)
		  tmp_tor.torsion -= 360;
	       pt_base_and_H = hydrogen_pos(tmp_tor, residue_reference);
	       
	       coot::bash_distance_t bash = find_bash_distance(pt_base_and_H.first,
							       pt_base_and_H.second,
							       atoms);
	       if (debug)
		  std::cout << "Adding bash distance " << bash << " to atom "
			    << named_torsions[i].atom_name_bonded_to_H
			    << std::endl;
	       atom_bashes[named_torsions[i].atom_name_bonded_to_H].push_back(bash);
	    }
	 }
      }
      catch (std::runtime_error rte) {
	 std::cout << rte.what() << std::endl;
      } 
   }
}

// RDKit version (well, the version that is used when the hydrogens
// are correctly named (according to the dictionary) and placed on the
// ligand of interest.
//
// This fills the atom_bashes vector.
// 
void
coot::flev_attached_hydrogens_t::distances_to_protein_using_correct_Hs(CResidue *ligand_residue,
								       CMMDBManager *mol,
								       const coot::protein_geometry &geom) {

   // the constructor (called just before this) should fill
   // atoms_with_rotating_hydrogens and atoms_with_riding_hydrogens
   // vectors (using the restraints).
   
   float radius = 6.0;
   std::vector<CResidue *> env_residues =
      coot::residues_near_residue(ligand_residue, mol, radius);

   // -------------------------------------------------------------------
   //                    riding hydrogens
   // -------------------------------------------------------------------
   // 
   PPCAtom residue_atoms = 0;
   int n_ligand_atoms;
   ligand_residue->GetAtomTable(residue_atoms, n_ligand_atoms);
   for (unsigned int irh=0; irh<atoms_with_riding_hydrogens.size(); irh++) { 
      CAtom *lig_at = NULL;
      CAtom *H_at = NULL;
      for (unsigned int iat=0; iat<n_ligand_atoms; iat++) { 
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == atoms_with_riding_hydrogens[irh].first)
	    lig_at = residue_atoms[iat];
	 if (atom_name == atoms_with_riding_hydrogens[irh].second)
	    H_at = residue_atoms[iat];
	 if (lig_at && H_at)
	    break;
      }
      if (lig_at && H_at) {
	 clipper::Coord_orth H_pt(H_at->x, H_at->y, H_at->z);
	 clipper::Coord_orth lig_atom_pt(lig_at->x, lig_at->y, lig_at->z);
	 
	 std::vector<CAtom *> atoms = close_atoms(H_pt, env_residues);
	 coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, H_pt, atoms);
	 atom_bashes[atoms_with_riding_hydrogens[irh].first].push_back(bash);
	 if (0) 
	    std::cout << " adding bash distance " << bash << " to atom " 
		      << atoms_with_riding_hydrogens[irh].first << std::endl;
      } 
   }
   
   // -------------------------------------------------------------------
   //                 rotatable hydrogens (more complex)
   // -------------------------------------------------------------------
   // 
   
   for (unsigned int irh=0; irh<atoms_with_rotating_hydrogens.size(); irh++) { 
      CAtom *lig_at = NULL;
      CAtom *H_at = NULL;
      for (unsigned int iat=0; iat<n_ligand_atoms; iat++) { 
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == atoms_with_rotating_hydrogens[irh].first)
	    lig_at = residue_atoms[iat];
	 if (atom_name == atoms_with_rotating_hydrogens[irh].second)
	    H_at = residue_atoms[iat];
	 if (lig_at && H_at)
	    break;
      }
      if (lig_at && H_at) {
	 clipper::Coord_orth H_pt(H_at->x, H_at->y, H_at->z);
	 clipper::Coord_orth lig_atom_pt(lig_at->x, lig_at->y, lig_at->z);
	 
	 std::vector<CAtom *> atoms = close_atoms(H_pt, env_residues);

	 try { 
	    clipper::Coord_orth vector_pt = get_atom_pos_bonded_to_atom(lig_at, H_at, // not H_at
									ligand_residue, geom);
	    clipper::Coord_orth base_ref_pt(0,0,0);
	    double tors = clipper::Coord_orth::torsion(base_ref_pt, vector_pt, lig_atom_pt, H_pt);
	    double dist = sqrt((lig_atom_pt - H_pt).lengthsq());
	    double angle = clipper::Coord_orth::angle(vector_pt, lig_atom_pt, H_pt);

	    int n_torsion_samples = 8;
	    for (unsigned int itor=0; itor<n_torsion_samples; itor++) {

	       double tmp_tor_d =  double(itor) * 360.0/double(n_torsion_samples);
	       double tmp_tor = clipper::Util::d2rad(tmp_tor_d);
	       tmp_tor += tors;
	       clipper::Coord_orth new_pt =
		  clipper::Coord_orth(base_ref_pt, vector_pt, lig_atom_pt, dist, angle, tmp_tor);
	       coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, new_pt, atoms);
	       atom_bashes[atoms_with_rotating_hydrogens[irh].first].push_back(bash);
	    }
	 }
	 catch (std::runtime_error rte) {
	    std::cout << rte.what() << std::endl;
	 } 
      }
   }

}



// Can throw an exception
// 
// Return the position of the H-ligand atom (the atom to which the H
// is attached) and the hydrogen position - in that order.
// 
std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::flev_attached_hydrogens_t::hydrogen_pos(const coot::named_torsion_t &named_tor,
					      CResidue *residue_p) const {
   clipper:: Coord_orth pt(0,0,0);
   clipper:: Coord_orth pt_ligand_atom(0,0,0);
   CAtom *at_1 = NULL;
   CAtom *at_2 = NULL;
   CAtom *at_3 = NULL;

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == named_tor.base_atom_name)
	 at_1 = residue_atoms[i];
      if (atom_name == named_tor.atom_name_2)
	 at_2 = residue_atoms[i];
      if (atom_name == named_tor.atom_name_bonded_to_H)
	 at_3 = residue_atoms[i];
   }

   if (! (at_1 && at_2 && at_3)) {
      throw(std::runtime_error("missing atoms in residue"));
   } else {
      clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
      clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
      clipper::Coord_orth pt_3(at_3->x, at_3->y, at_3->z);
      clipper::Coord_orth p4_h(pt_1, pt_2, pt_3,
			       named_tor.dist,
			       clipper::Util::d2rad(named_tor.angle),
			       clipper::Util::d2rad(named_tor.torsion));
      pt = p4_h;
      pt_ligand_atom = pt_3;
      
//       std::cout << "in hydrogen_pos() constructed H pos " << pt.format()
// 		<< " from atom at_1: " << at_1 << " "
// 		<< "atom at_2: " << at_2 << " "
// 		<< "atom at_3: " << at_3 << " "
// 		<< " dist: " << named_tor.dist << " " 
// 		<< " angle: " << named_tor.angle << " " 
// 		<< " torsion: " << named_tor.torsion << std::endl;
   } 
   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (pt_ligand_atom, pt);
}

// find an atom (the atom, perhaps) bonded to lig_at that is not H_at.
// Return its position.
// 
// Can throw a std::runtime_error if not found.
// 
clipper::Coord_orth
coot::flev_attached_hydrogens_t::get_atom_pos_bonded_to_atom(CAtom *lig_at, CAtom *H_at, // not H_at
							     CResidue *ligand_residue,
							     const coot::protein_geometry &geom) const {
   std::string res_name(lig_at->residue->GetResName());
   std::pair<bool, dictionary_residue_restraints_t> p = 
      geom.get_monomer_restraints_at_least_minimal(res_name);

   if (! p.first) {
      std::string m = "No monomer type ";
      m += res_name;
      m += " found in dictionary";
      throw std::runtime_error(m);
   } else {
      CAtom *bonded_atom = NULL;
      std::string bonded_atom_name;
      std::string lig_at_name = lig_at->name;
      std::string H_at_name = H_at->name;
      for (unsigned int ibond=0; ibond<p.second.bond_restraint.size(); ibond++) { 
	 std::string atom_name_1 = p.second.bond_restraint[ibond].atom_id_1_4c();
	 std::string atom_name_2 = p.second.bond_restraint[ibond].atom_id_2_4c();
	 if (atom_name_1 == lig_at_name) {
	    if (atom_name_2 != H_at_name) {
	       bonded_atom_name = atom_name_2;
	       break;
	    }
	 }
	 if (atom_name_2 == lig_at_name) {
	    if (atom_name_1 != H_at_name) {
	       bonded_atom_name = atom_name_1;
	       break;
	    }
	 }
      }
      if (bonded_atom_name != "") {
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 ligand_residue->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (unsigned int iat=0; iat<n_residue_atoms ; iat++) {
	    std::string atom_name = residue_atoms[iat]->name;
	    if (atom_name == bonded_atom_name) {
	       bonded_atom = residue_atoms[iat];
	       break;
	    } 
	 }
      }

      if (! bonded_atom) {
	 std::string m = "No atom bonded to ";
	 m += lig_at_name;
	 m += " found in dictionary for ";
	 m += res_name;
	 throw std::runtime_error(m);
      } else {
	 // good
	 return clipper::Coord_orth(bonded_atom->x, bonded_atom->y, bonded_atom->z);
      }
   }

}


// What are the atoms that are close (distance < 6A) to pt?
// 
// waters are not counted as close atoms.
// 
std::vector<CAtom *>
coot::flev_attached_hydrogens_t::close_atoms(const clipper::Coord_orth &pt,
					     const std::vector<CResidue *> &env_residues) const {

   std::vector<CAtom *> v;
   double dist_crit = 6.0;
   double dist_crit_squared = dist_crit * dist_crit;
   
   for (unsigned int i=0; i<env_residues.size(); i++) {
      CResidue *residue_p = env_residues[i];
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 clipper::Coord_orth atom_pos(residue_atoms[iat]->x, residue_atoms[iat]->y, residue_atoms[iat]->z);
	 double d_squared = (pt - atom_pos).lengthsq();
	 if (d_squared < dist_crit_squared) { 
	    std::string rn(residue_atoms[iat]->GetResName());
	    if (rn != "HOH")
	       v.push_back(residue_atoms[iat]);
	 }
      }
   }
   return v;
}


coot::bash_distance_t
coot::flev_attached_hydrogens_t::find_bash_distance(const clipper::Coord_orth &ligand_atom_pos,
						    const clipper::Coord_orth &hydrogen_pos,
						    const std::vector<CAtom *> &close_residue_atoms) const {

   double cannonball_radius = 0.8; // radius of the cannonball, c.f. at least a hydrogen.
   
   double max_dist = 4.05; // if we have travelled 4A without bashing
			   // into a protein atom then this has
			   // practically unlimited substitution distance.

   double max_dist_squared = max_dist * max_dist;
   clipper::Coord_orth h_vector((hydrogen_pos - ligand_atom_pos).unit());

   if (0)
      std::cout << "h_vector: " << h_vector.format() << " from hydrogen pos: "
		<< hydrogen_pos.format() << "and ligand atom pos: " << ligand_atom_pos.format()
		<< std::endl;

   // set the atomic radii:
   // 
   std::vector<double> radius(close_residue_atoms.size());
   for (unsigned int iat=0; iat<close_residue_atoms.size(); iat++) {
      std::string ele(close_residue_atoms[iat]->element);
      radius[iat] = get_radius(ele);
   }
   
   coot::bash_distance_t dd;
   
   std::vector<clipper::Coord_orth> atom_positions(close_residue_atoms.size());
   // likewise set the atom positions so that we don't have to keep doing it.
   for (unsigned int i=0; i<close_residue_atoms.size(); i++)
      atom_positions[i] = clipper::Coord_orth(close_residue_atoms[i]->x,
					      close_residue_atoms[i]->y,
					      close_residue_atoms[i]->z);
   
   for (double slide=0; slide<=max_dist; slide+=0.04) {
      clipper::Coord_orth test_pt = ligand_atom_pos + slide * h_vector;
      if (0)
	 std::cout << "   bash distance for ligand atom at " << ligand_atom_pos.format() << " "
		   << "determined from " << atom_positions.size() << " atom positions"
		   << std::endl;
      for (unsigned int iat=0; iat<atom_positions.size(); iat++) {
	 double atom_radius_plus_cbr = radius[iat] + cannonball_radius;
	 double d_squared = (test_pt - atom_positions[iat]).lengthsq();
	 if (0)
	    std::cout << "   atom " << iat << " "
		      << close_residue_atoms[iat]->GetChainID() << " "
		      << close_residue_atoms[iat]->GetSeqNum() << " "
		      << close_residue_atoms[iat]->GetAtomName() << " "
		      << " slide: " << slide
		      << " comparing " << sqrt(d_squared) << "^2  and "
 		   << atom_radius_plus_cbr << "^2" << std::endl;
	 if (d_squared < atom_radius_plus_cbr*atom_radius_plus_cbr) {
	    dd = coot::bash_distance_t(slide);
	    break;
	 }
      }
      if (dd.limited) {
	 break;
      }
   }
   return dd;
}

// c.f. 
// double
// coot::dots_representation_info_t::get_radius(const std::string &ele)

double
coot::flev_attached_hydrogens_t::get_radius(const std::string &ele) const { 

   double radius = 1.70;
   if (ele == " H")
      radius = 1.20;
   if (ele == " N")
      radius = 1.55;
   if (ele == " O")
      radius = 1.52;
   if (ele == " S")
      radius = 1.8;
   return radius;
} 
