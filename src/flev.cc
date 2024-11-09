/* src/flev.cc
 * 
 * Copyright 2010, 2011, 2012 The University of Oxford
 * Copyright 2014 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */

/*  ----------------------------------------------------------------------- */
/*               Flattened Ligand Environment View  Interface               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


// For reasons I don't understand, this should come near the top of
// includes, otherwise we get RDKit dcgettext() include file problems.
//
// #include "lbg/lbg.hh"

#include "c-interface-ligands.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "geometry/lbg-graph.hh"

#include "coot-utils/coot-h-bonds.hh"
#include "pli/protein-ligand-interactions.hh"

// 20140311: also, we want to include pi-stacking explicitly here
// because we need it to compile flev if we don't HAVE_GOOCANVAS
// (pi-stacking.hh is included from lbg.hh if HAVE_GOOCANVAS).
// 
#include "pli/pi-stacking.hh"

#include "flev.hh"

#include "cc-interface.hh" // for add_animated_ligand_interactions


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

//  return true if both ligand_res and env_res have hydrogens in the
//  dictionary and (at least one) hydrogen in the model
// (if ligand_res is a SO4 (no hydrogens) that should not cause a fail).
// 
bool flev_check_for_hydrogens(int imol, mmdb::Residue *ligand_res,
			      const std::vector<mmdb::Residue *> &env_residues,
			      coot::protein_geometry *geom_p) {

   bool status = false;
   int n_hydrogens_ligand = 0;
   int n_hydrogens_others = 0;

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   ligand_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) { 
      std::string ele = residue_atoms[i]->element;
      if (ele == " H") { // needs PDBv3 update
	 n_hydrogens_ligand++;
      } 
   }
   for (unsigned int ires=0; ires<env_residues.size(); ires++) { 
      residue_atoms = 0;
      n_residue_atoms = 0;
      env_residues[ires]->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) { 
	 std::string ele = residue_atoms[iat]->element;
	 if (ele == " H") { // needs PDBv3 update
	    n_hydrogens_others++;
	 }
      }
   }
   if (n_hydrogens_others > 0)
      if (n_hydrogens_ligand > 0)
	 return true;


   if (n_hydrogens_ligand == 0) {
      // perhaps the ligand was not supposed to have any ligands
      std::string residue_type = ligand_res->GetResName();
      if (geom_p->have_dictionary_for_residue_type(residue_type, imol,
						   graphics_info_t::cif_dictionary_read_number++)) {
	 int n_hydrogens_ligand_dict = geom_p->n_hydrogens(residue_type);
	 if (n_hydrogens_ligand_dict > 0)
	    status = false;
	 else
	    if (n_hydrogens_others > 0)
	       // environment has hydrogens, ligand does not (and is
	       // not supposed to have)
	       status = true;
	    else
	       // environment does not have hydrogens (and is supposed
	       // to have (we presume))
	       status = false;
      }
   }
   return status;
} 


void fle_view_internal(int imol, const char *chain_id, int res_no, const char *ins_code, 
		       int imol_ligand_fragment, 
		       const char *prodrg_output_flat_mol_file_name,
		       const char *prodrg_output_flat_pdb_file_name,
		       const char *prodrg_output_3d_pdb_file_name,
		       const char *prodrg_output_dict_cif_file_name) {

   fle_view_internal_to_png(imol, chain_id, res_no, ins_code, imol_ligand_fragment,
			    prodrg_output_flat_mol_file_name,
			    prodrg_output_flat_pdb_file_name,
			    prodrg_output_3d_pdb_file_name,
			    prodrg_output_dict_cif_file_name, 0, "");

}

/* for command-line operation */
void fle_view_internal_to_png(int imol, const char *chain_id, int res_no, 
			      const char *ins_code, 
			      int imol_ligand_fragment, 
			      const char *prodrg_output_flat_mol_file_name,
			      const char *prodrg_output_flat_pdb_file_name,
			      const char *prodrg_output_3d_pdb_file_name,
			      const char *prodrg_output_dict_cif_file_name,
			      int output_to_png_file_flag,
			      const char *png_file_name) {

   if (0) 
      std::cout << "debug:: fle_view_internal() called with " 
		<< imol << " " 
		<< chain_id << " " 
		<< res_no << " " 
		<< ins_code << " " 
		<< imol_ligand_fragment << " " 
		<< prodrg_output_flat_mol_file_name << " " 
		<< prodrg_output_flat_pdb_file_name << " " 
		<< prodrg_output_3d_pdb_file_name   << " " 
		<< prodrg_output_dict_cif_file_name << " " 
		<< std::endl;
   
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

   atom_selection_container_t flat = get_atom_selection(prodrg_output_flat_pdb_file_name, false, true, false);
   if (flat.read_success) {
      if (is_valid_model_molecule(imol_ligand_fragment)) {

	 mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 mmdb::Residue  *res_ref = coot::util::get_residue(chain_id, res_no, ins_code, mol);
	 mmdb::Residue *flat_res = coot::util::get_residue("", 1, "", flat.mol); // prodrg attribs

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
	    mmdb::Manager *ligand_mol =
	       graphics_info_t::molecules[imol_ligand_fragment].atom_sel.mol;
	    int every_nth = 1;
	    std::vector<coot::lsq_range_match_info_t> matches;
	    coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id, coot::lsq_t::ALL);
	    matches.push_back(match);
	    std::pair<short int, clipper::RTop_orth> lsq_mat = 
	       coot::util::get_lsq_matrix(flat.mol, ligand_mol, matches, every_nth);
	    
	    if (lsq_mat.first) { 
	       std::vector<mmdb::Residue *> residues =
		  coot::residues_near_residue(res_ref, mol, residues_near_radius);
	       
	       // residues needs to be filtered to remove waters that
	       // are not connected to a protein atom.

	       // 20101228 water_max_dist was 3.6, not we tweak it to
	       // 3.25.  Should be a user-setable param.
	       // 
	       std::vector<mmdb::Residue *> filtered_residues =
		  coot::filter_residues_by_solvent_contact(res_ref, mol, residues, water_dist_max);
	       
	       // for the atoms in the ligand only, for the moment.
	       // More would require a class containing with and
	       // without solvent accessibilites for each residue.
	       // Maybe that can be another function.  And that would
	       // need to consider neighbours of neighbours, perhaps
	       // done with a larger radius.
	       //
	       pli::dots_representation_info_t dots;
	       std::vector<std::pair<coot::atom_spec_t, float> > s_a_v = 
		  dots.solvent_accessibilities(res_ref, filtered_residues);

#if 0 // no longer do we use lbg. Put this lbg function in pli I guess
	       // for the ligand environment residues:
	       std::vector<coot::solvent_exposure_difference_helper_t> sed = 
		  dots.solvent_exposure_differences(res_ref, filtered_residues);
#endif

	       if (false)
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
	       std::vector<pli::fle_ligand_bond_t> bonds_to_ligand =
		  pli::get_fle_ligand_bonds(res_ref, filtered_residues, mol,
					     name_map, *geom_p, imol,
					     graphics_info_t::fle_water_dist_max,
					     graphics_info_t::fle_h_bond_dist_max);

	       // add_animated_ligand_interactions(imol, bonds_to_ligand);

	       if (true)
		  std::cout << "Found ================== " << bonds_to_ligand.size()
			    << " ==================== bonds to ligand " << std::endl;

	       if (filtered_residues.size()) {
		  std::vector<pli::fle_residues_helper_t> centres(filtered_residues.size());
		  for (unsigned int ires=0; ires<filtered_residues.size(); ires++) { 
		     mmdb::Residue *res_copy = coot::util::deep_copy_this_residue(filtered_residues[ires]);
		     std::string res_name = filtered_residues[ires]->GetResName();
		     coot::util::transform_atoms(res_copy, lsq_mat.second);
		     std::pair<bool, clipper::Coord_orth> c =
			coot::util::get_residue_centre(res_copy);
		     if (c.first) {
			if (0) 
			   std::cout << "DEBUG:: creating fle_centre with centre "
				     << c.second.format() << std::endl;
			pli::fle_residues_helper_t fle_centre(c.second,
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
		     geom_p->get_monomer_restraints(ligand_res_name, imol);
		  
		  if (! p.first) {
		     std::cout << "WARNING:: fle_view_internal(): "
			       << "Failed to get monomer_restraints for PRODRG residue"
			       << std::endl;
		  } else {

		     // ----------- residue infos ----------
		     // 
		     pli::pi_stacking_container_t pi_stack_info(p.second, filtered_residues, res_ref);

#ifdef VERY_OLD_FLEV_FUNCTIONS		     
		     write_fle_centres(centres, bonds_to_ligand, sed, pi_stack_info, flat_res);
#endif		     
		     // ----------- ligand atom infos ------
		     // 
		     pli::flev_attached_hydrogens_t ah(p.second);
		     ah.cannonballs(res_ref, prodrg_output_3d_pdb_file_name, p.second);
		     ah.distances_to_protein(res_ref, mol);
		     write_ligand_atom_accessibilities(s_a_v, ah, flat_res);

		     std::pair<bool, coot::residue_spec_t>
			ligand_spec_pair(1, coot::residue_spec_t(res_ref));
		     
		     std::string view_name;
		     
		     bool use_graphics_flag = graphics_info_t::use_graphics_interface_flag;
		     if (output_to_png_file_flag)
			use_graphics_flag = false;
		     
		     bool stand_alone_flag = 0; // no, it isn't from here.

		     lig_build::molfile_molecule_t m;
		     m.read(prodrg_output_flat_mol_file_name);

		     std::vector<pli::fle_residues_helper_t> res_centres =
			coot::get_flev_residue_centres(res_ref,
						       mol, // mol_for_res_ref,
						       filtered_residues,
						       flat.mol // mol_for_flat_residue
						       );
		     std::vector<int> add_reps_vec;
		     if (graphics_info_t::use_graphics_interface_flag)
		     add_reps_vec = coot::make_add_reps_for_near_residues(filtered_residues, imol);
		     
		     if (0) { 
			std::cout << "------------- in flev: centres.size() is "
				  << res_centres.size() << std::endl;
			for (unsigned int ic=0; ic<res_centres.size(); ic++)
			   std::cout << "   " << ic << "  " << res_centres[ic]
				     << std::endl;
		     }
		     

#ifdef HAVE_GOOCANVAS		     
		     lbg_info_t *lbg_local_p = lbg(m, ligand_spec_pair,
						   flat.mol, view_name, ligand_res_name, imol,
						   graphics_info_t::Geom_p(),
						   use_graphics_flag, stand_alone_flag,
						   coot_get_url,
						   prodrg_import_function,
						   sbase_import_function,
						   get_drug_mdl_via_wikipedia_and_drugbank);
		     lbg_local_p->annotate(s_a_v, res_centres, add_reps_vec, bonds_to_ligand, sed,
					   ah, pi_stack_info, p.second);

		     lbg_local_p->set_orient_view_func(orient_view);
		     lbg_local_p->set_set_rotation_centre_func(set_rotation_centre);
		     lbg_local_p->set_set_show_additional_representation_func(set_show_additional_representation);
		     lbg_local_p->set_all_additional_representations_off_except_func(all_additional_representations_off_except);
		     
		     if (! use_graphics_flag) { 
			lbg_local_p->set_draw_flev_annotations(true);
			lbg_local_p->draw_all_flev_annotations();
			lbg_local_p->write_png(png_file_name);
		     } 

#endif // HAVE_GOOCANVAS		     
		     

		  }
	       }
	    }
	 }
      }
   } 
}

void fle_view_with_rdkit(int imol, const char *chain_id, int res_no,
			 const char *ins_code, float residues_near_radius) {

   fle_view_with_rdkit_internal(imol, chain_id, res_no, ins_code, residues_near_radius, "", "");
}

void fle_view_with_rdkit_to_png(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *png_file_name) {

   fle_view_with_rdkit_internal(imol, chain_id, res_no, ins_code, residues_near_radius, "png", png_file_name);
}

void fle_view_with_rdkit_to_svg(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *png_file_name) {

   fle_view_with_rdkit_internal(imol, chain_id, res_no, ins_code, residues_near_radius, "svg", png_file_name);
} 


void fle_view_with_rdkit_internal(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *file_format, const char *output_image_file_name) {

#if COMPILE_WITH_LBG // don't configure this test when the time comes, just delete this line

#ifndef MAKE_ENHANCED_LIGAND_TOOLS

   std::cout << "WARNING:: fle_view_with_rdkit_internal() not enhanced ligand " << std::endl;

# else

   double weight_for_3d_distances = 0.4; // for 3d distances
   double water_dist_max = 3.25;
   
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();
   std::string output_format = file_format;

   if (is_valid_model_molecule(imol)) { 
      mmdb::Residue  *res_ref = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      mmdb::Manager *mol_for_res_ref = g.molecules[imol].atom_sel.mol;
      if (res_ref) {
	 std::string ligand_res_name(res_ref->GetResName());

	 std::pair<bool, coot::dictionary_residue_restraints_t> p = 
	    geom_p->get_monomer_restraints_at_least_minimal(ligand_res_name, imol);
	 
	 if (! p.first) {
	    std::string s1 = "WARNING:: fle_view_with_rdkit(): ";
	    std::string s2 = "WARNING:: ";
	    s1 += "Failed to get \nmonomer_restraints for ligand of type ";
	    s2 += "Failed to get \nmonomer_restraints for ligand of type ";
	    s1 += ligand_res_name;
	    s2 += ligand_res_name;
	    std::cout << s1 << std::endl;
	    info_dialog(s2.c_str());
	 } else {
	    std::vector<mmdb::Residue *> residues =
	       coot::residues_near_residue(res_ref, mol_for_res_ref, residues_near_radius);
	 
	    // residues needs to be filtered to remove waters that
	    // are not connected to a protein atom.
	 
	    std::vector<mmdb::Residue *> filtered_residues =
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
	       RDKit::RWMol rdkm = coot::rdkit_mol(res_ref, imol, *g.Geom_p());

	       // assign atom names
	       if (int(rdkm.getNumAtoms()) < res_ref->GetNumberOfAtoms()) {
		  std::cout << "WARNING:: failure to construct rdkit molecule " << rdkm.getNumAtoms()
                            << " vs " << res_ref->GetNumberOfAtoms() << std::endl;
	       } else {
		  mmdb::PPAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_ref->GetAtomTable(residue_atoms, n_residue_atoms);

		  // polar Hs only, that is - need new function here.
		  // (can throw a std::exception)
		  coot::undelocalise(&rdkm);
		  coot::assign_formal_charges(&rdkm);
		  coot::remove_non_polar_Hs(&rdkm);

		  // we need to sanitizeMol() after remove_non_polar_Hs, and have
		  // a kekulized representation.
		  // we need to sanitize to get ring info,
		  // then we need to kekulize because we are making a 2D chemical diagram
		  //
		  // failed_op sets set by sanitizeMol() - perhaps we shouldn't ignore it.
		  //
		  unsigned int failed_op_1 = 0;
		  unsigned int failed_op_2 = 0;
		  RDKit::MolOps::sanitizeMol(rdkm, failed_op_1, RDKit::MolOps::SANITIZE_ALL);
		  RDKit::MolOps::sanitizeMol(rdkm, failed_op_2, RDKit::MolOps::SANITIZE_KEKULIZE);

		  std::cout << "DEBUG:: sanitizeMol() returned with failed_op: "
			    << failed_op_1 << " " << failed_op_2
			    << " (note 'no-failure' is value 0)." << std::endl;

		  try {
		     RDKit::RingInfo *ri = rdkm.getRingInfo();
		     unsigned int n_rings = ri->numRings();
		  }
		  catch (const std::runtime_error &rte) {
		     std::vector<std::vector<int> > ring_info;
		     RDKit::MolOps::findSSSR(rdkm, ring_info);
		  }


		  int mol_2d_depict_conformer =
		     coot::add_2d_conformer(&rdkm, weight_for_3d_distances);
		  lig_build::molfile_molecule_t m =
		     coot::make_molfile_molecule(rdkm, mol_2d_depict_conformer);

		  mmdb::Residue *residue_flat =
		     coot::make_residue(rdkm, mol_2d_depict_conformer, "XXX");
		  mmdb::Manager *mol_for_flat_residue =
		     coot::util::create_mmdbmanager_from_residue(residue_flat);

		  if (false)
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
						graphics_info_t::Geom_p(),
						use_graphics_flag, stand_alone_flag,
						coot_get_url,
						prodrg_import_function,
						sbase_import_function,
						get_drug_mdl_via_wikipedia_and_drugbank);

		  lbg_local_p->set_orient_view_func(orient_view);
		  lbg_local_p->set_set_rotation_centre_func(set_rotation_centre);
		  lbg_local_p->set_set_show_additional_representation_func(set_show_additional_representation);
		  lbg_local_p->set_all_additional_representations_off_except_func(all_additional_representations_off_except);

 		  std::map<std::string, std::string> name_map =
 		     coot::make_flat_ligand_name_map(res_ref);
	       
 		  std::vector<coot::fle_ligand_bond_t> bonds_to_ligand =
 		     coot::get_fle_ligand_bonds(res_ref, filtered_residues,
						mol_for_res_ref, name_map, *geom_p, 
						graphics_info_t::fle_water_dist_max,
						graphics_info_t::fle_h_bond_dist_max);

		  if (graphics_info_t::use_graphics_interface_flag)
		     add_animated_ligand_interactions(imol, bonds_to_ligand);

		  std::vector<coot::fle_residues_helper_t> res_centres =
		     coot::get_flev_residue_centres(res_ref,
						    mol_for_res_ref,
						    filtered_residues,
						    mol_for_flat_residue);

		  std::vector<int> add_reps_vec;
		  if (graphics_info_t::use_graphics_interface_flag)
		     add_reps_vec = coot::make_add_reps_for_near_residues(filtered_residues, imol);

		  if (0) { 
		     for (unsigned int ic=0; ic<res_centres.size(); ic++)
			std::cout << "   " << ic << "  " << res_centres[ic]
				  << std::endl;
		  }

		  // ----------- residue infos ----------
		  //
		  coot::pi_stacking_container_t
		     pi_stack_info(p.second, filtered_residues, res_ref, rdkm);

		     
		  // ----------- ligand atom infos ------
		  // 
 		  coot::flev_attached_hydrogens_t ah(p.second);
		  // ah.cannonballs(res_ref, mol_for_res_ref, p.second);
 		  ah.distances_to_protein_using_correct_Hs(res_ref, mol_for_res_ref, *geom_p);

		  // ------------ show it -----------------
		  //
 		  lbg_local_p->annotate(s_a_v, res_centres, add_reps_vec, bonds_to_ligand, sed,
					ah, pi_stack_info, p.second);
		  if (! use_graphics_flag) { 
		     lbg_local_p->set_draw_flev_annotations(true);
		     lbg_local_p->draw_all_flev_annotations();
		     if (output_format == "png")
			lbg_local_p->write_png(output_image_file_name);
		     if (output_format == "svg")
			lbg_local_p->write_svg(output_image_file_name);
		  } 
		  delete mol_for_flat_residue;
	       }
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "ERROR:: (runtime error) in fle_view_with_rdkit(): "
			 << rte.what() << std::endl;
	    } 
	    catch (const std::exception &e) {
	       std::cout << "ERROR (exception) in fle_view_with_rdkit(): " << e.what() << std::endl;
	    } 
	 }
      }
   }
#endif // MAKE_ENHANCED_LIGAND_TOOLS   
#endif // COMPILE_WITH_LBG
}

void fle_view_set_water_dist_max(float dist_max) { 
   graphics_info_t::fle_water_dist_max = dist_max;   // default 3.25
} 

void fle_view_set_h_bond_dist_max(float h_bond_dist_max) {
   graphics_info_t::fle_h_bond_dist_max = h_bond_dist_max;   // default 3.9
} 


std::vector<int>
coot::make_add_reps_for_near_residues(std::vector<mmdb::Residue *> filtered_residues,
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
				       const std::vector<pli::fle_ligand_bond_t> &ligand_bonds) {

   for (unsigned int i=0; i<ligand_bonds.size(); i++) {
      // std::cout << "Here....  adding animated ligand interaction " << i << std::endl;
      add_animated_ligand_interaction(imol, ligand_bonds[i]);
   }
}


std::vector<pli::fle_residues_helper_t>
coot::get_flev_residue_centres(mmdb::Residue *residue_ligand_3d,
			       mmdb::Manager *mol_containing_residue_ligand, 
			       std::vector<mmdb::Residue *> residues,
			       mmdb::Manager *flat_mol) {

   std::vector<pli::fle_residues_helper_t> centres;

   if (flat_mol) { 

      // get the lsq matrix that maps the ligand in 3D onto the flat ligand
      int res_no = residue_ligand_3d->GetSeqNum();
      std::string chain_id = residue_ligand_3d->GetChainID();
      int every_nth = 1;
      std::vector<coot::lsq_range_match_info_t> matches;
      coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
 					 coot::lsq_t::ALL);
       matches.push_back(match);
       std::pair<short int, clipper::RTop_orth> lsq_mat =
	  coot::util::get_lsq_matrix(flat_mol, mol_containing_residue_ligand, matches, every_nth);
      // Now make the residues

      // std::vector<coot::fle_residues_helper_t> centres(residues.size());
      centres.resize(residues.size());
      for (unsigned int ires=0; ires<residues.size(); ires++) { 
	 mmdb::Residue *res_copy = coot::util::deep_copy_this_residue(residues[ires]);
	 std::string res_name = residues[ires]->GetResName();
	 std::pair<bool, clipper::Coord_orth> absolute_centre =
	    coot::util::get_residue_centre(res_copy);
	 if (absolute_centre.first) { 
	    coot::util::transform_atoms(res_copy, lsq_mat.second);
	    std::pair<bool, clipper::Coord_orth> c =
	       coot::util::get_residue_centre(res_copy);
	    if (c.first) {
	       pli::fle_residues_helper_t fle_centre(c.second,
                                                     coot::residue_spec_t(residues[ires]),
                                                     res_name);

	       // Setting the interaction position to the residue
	       // centre is a hack.  What we need to be is the middle
	       // of the hydrogen bond or the middle of the pi-pi
	       // stacking (for instance).  To do that we need to know
	       // what positions of the interacting in the ligand (and
	       // of the residue).  Tricky from here?
	       // 
	       fle_centre.set_interaction_position(absolute_centre.second);
	       centres[ires] = fle_centre;
	    } else {
	       std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
			 << coot::residue_spec_t(res_copy) << std::endl;
	    }
	 } else {
	    std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
		      << coot::residue_spec_t(res_copy) << std::endl;
	 }
	 delete res_copy;
      }
   }
   if (0) 
      for (unsigned int ic=0; ic<centres.size(); ic++) { 
	 std::cout << "centre " << ic << " has TR_centre "
		   << centres[ic].transformed_relative_centre.format() << std::endl;
      }
   return centres;
} 


std::map<std::string, std::string>
coot::make_flat_ligand_name_map(mmdb::Residue *flat_res) {

   double bond_to_H_dist = 1.1;

   double b2Hd2 = bond_to_H_dist * bond_to_H_dist;
   std::map<std::string, std::string> map;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   flat_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at_i = residue_atoms[iat];
      std::string ele_i = at_i->element;
      clipper::Coord_orth pt_i(at_i->x, at_i->y, at_i->z);
      if (ele_i == " H") { 
	 for (int jat=0; jat<n_residue_atoms; jat++) {
	    if (iat != jat) { 
	       mmdb::Atom *at_j = residue_atoms[jat];
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


 
#ifdef VERY_OLD_FLEV_FUNCTIONS		     

void
coot::write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			const coot::pi_stacking_container_t &stack_info,
			mmdb::Residue *res_flat) {

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
	       << v[i].residue_name << " " << v[i].spec.chain << v[i].spec.res_no;

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
		  mmdb::PPAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_flat->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at = residue_atoms[iat];
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
		     << v[i].residue_name << " " << v[i].spec.chain_id
		     << v[i].spec.res_no << " is " << best_stacking.type << " "
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

#endif //  VERY_OLD_FLEV_FUNCTIONS		     


// the reference_residue is the flat residue
void
coot::write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					const pli::flev_attached_hydrogens_t &attached_hydrogens,
					mmdb::Residue *reference_residue) {

   // for each of the atoms in reference_residue, find the spec in sav
   // and write out the coordinates and the atom name and the
   // solvent_accessibility.

   std::string file_name = "coot-tmp-fle-view-solvent-accessibilites.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      reference_residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 mmdb::Atom *at = residue_atoms[i];
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

std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
coot::get_cannonball_vectors(mmdb::Residue *ligand_res_3d,
			     const coot::dictionary_residue_restraints_t &monomer_restraints) {

   std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > > v;

   return v;
}

