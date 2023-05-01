/* src/flev.cc
 * 
 * Copyright 2010, 2011, 2012 The University of Oxford
 * Copyright 2014 by Medical Research Council
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

#include "compat/coot-sysdep.h"


// For reasons I don't understand, this should come near the top of
// includes, otherwise we get RDKit dcgettext() include file problems.
//
#include "lbg/lbg.hh"

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

   atom_selection_container_t flat = get_atom_selection(prodrg_output_flat_pdb_file_name, true, false, false);
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
	    coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
					       COOT_LSQ_ALL);
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
 					     name_map, *geom_p,
					     graphics_info_t::fle_water_dist_max,
					     graphics_info_t::fle_h_bond_dist_max);

	       add_animated_ligand_interactions(imol, bonds_to_ligand);

	       if (1) 
		  std::cout << "Found ================== " << bonds_to_ligand.size()
			    << " ==================== bonds to ligand " << std::endl;

	       if (filtered_residues.size()) {
		  std::vector<coot::fle_residues_helper_t> centres(filtered_residues.size());
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
		     geom_p->get_monomer_restraints(ligand_res_name, imol);
		  
		  if (! p.first) {
		     std::cout << "WARNING:: fle_view_internal(): "
			       << "Failed to get monomer_restraints for PRODRG residue"
			       << std::endl;
		  } else {

		     // ----------- residue infos ----------
		     // 
		     coot::pi_stacking_container_t pi_stack_info(p.second, filtered_residues, res_ref);

#ifdef VERY_OLD_FLEV_FUNCTIONS		     
		     write_fle_centres(centres, bonds_to_ligand, sed, pi_stack_info, flat_res);
#endif		     
		     // ----------- ligand atom infos ------
		     // 
		     coot::flev_attached_hydrogens_t ah(p.second);
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

		     std::vector<coot::fle_residues_helper_t> res_centres =
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
		  std::cout << "WARNING:: failure to construct rdkit molecule " << rdkm.getNumAtoms() << " vs " << res_ref->GetNumberOfAtoms()
	                    << std::endl;
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
				       const std::vector<coot::fle_ligand_bond_t> &ligand_bonds) {

   for (unsigned int i=0; i<ligand_bonds.size(); i++) {
      // std::cout << "Here....  adding animated ligand interaction " << i << std::endl;
      add_animated_ligand_interaction(imol, ligand_bonds[i]);
   }
}



std::ostream &
coot::operator<<(std::ostream &s, fle_residues_helper_t fler) {

   s << fler.is_set;
   if (fler.is_set) {
      s << " " << fler.transformed_relative_centre.format() << " "
	<< fler.spec << " " << fler.residue_name;
   }
   return s;
}


std::vector<coot::fle_residues_helper_t>
coot::get_flev_residue_centres(mmdb::Residue *residue_ligand_3d,
			       mmdb::Manager *mol_containing_residue_ligand, 
			       std::vector<mmdb::Residue *> residues,
			       mmdb::Manager *flat_mol) {

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
	 mmdb::Residue *res_copy = coot::util::deep_copy_this_residue(residues[ires]);
	 std::string res_name = residues[ires]->GetResName();
	 std::pair<bool, clipper::Coord_orth> absolute_centre =
	    coot::util::get_residue_centre(res_copy);
	 if (absolute_centre.first) { 
	    coot::util::transform_atoms(res_copy, lsq_mat.second);
	    std::pair<bool, clipper::Coord_orth> c =
	       coot::util::get_residue_centre(res_copy);
	    if (c.first) {
	       coot::fle_residues_helper_t fle_centre(c.second,
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
					const coot::flev_attached_hydrogens_t &attached_hydrogens,
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
coot::flev_attached_hydrogens_t::cannonballs(mmdb::Residue *ligand_residue_3d,
					     const std::string &prodrg_3d_ligand_file_name,
					     const coot::dictionary_residue_restraints_t &restraints) {

   atom_selection_container_t asc = get_atom_selection(prodrg_3d_ligand_file_name, true, false, false);
   if (asc.read_success) {
      cannonballs(ligand_residue_3d, asc.mol, restraints);
   }
}

void
coot::flev_attached_hydrogens_t::cannonballs(mmdb::Residue *ligand_residue_3d,
					     mmdb::Manager *mol, 
					     const coot::dictionary_residue_restraints_t &restraints) {

   if (! mol)
      return;
   
   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;


   int SelHnd_H = mol->NewSelection();
   int SelHnd_non_H = mol->NewSelection();

   mmdb::PPAtom hydrogen_selection = 0;
   mmdb::PPAtom non_hydrogen_selection = 0;
   int n_hydrogen_atoms;
   int n_non_hydrogen_atoms;
      
      
   mol->SelectAtoms(SelHnd_H,     0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", " H", "*");
   mol->SelectAtoms(SelHnd_non_H, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "!H", "*");
      
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
	 mmdb::Atom *at = non_hydrogen_selection[pscontact[i].id2];
	 std::string atom_name_bonded_to_H(at->name);

	 bool found_torsion_for_this_H = 0;

	 // riding hydrogens:
	 // 
	 for (unsigned int iat=0; iat<atoms_with_riding_hydrogens.size(); iat++) { 
	    if (atom_name_bonded_to_H == atoms_with_riding_hydrogens[iat].first) {
	       mmdb::Atom *h_at = hydrogen_selection[pscontact[i].id1];
	       found_torsion_for_this_H = add_named_torsion(h_at, at, restraints, mol, coot::H_IS_RIDING);
	    }
	    if (found_torsion_for_this_H)
	       break;
	 }

	 // rotating hydrogens:
	 // 
	 for (unsigned int iat=0; iat<atoms_with_rotating_hydrogens.size(); iat++) { 
	    if (atom_name_bonded_to_H == atoms_with_rotating_hydrogens[iat].first) {
	       mmdb::Atom *h_at = hydrogen_selection[pscontact[i].id1];
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
coot::flev_attached_hydrogens_t::add_named_torsion(mmdb::Atom *h_at, mmdb::Atom *at,
						   const coot::dictionary_residue_restraints_t &restraints,
						   mmdb::Manager *mol, // 3d prodrg ligand mol
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
			mmdb::Atom *At_2 = NULL;
			mmdb::Atom *base_atom = NULL;

			int imod = 1;
			mmdb::Model *model_p = mol->GetModel(imod);
			mmdb::Chain *chain_p;
			int nchains = model_p->GetNumberOfChains();
			for (int ichain=0; ichain<nchains; ichain++) {
			   chain_p = model_p->GetChain(ichain);
			   int nres = chain_p->GetNumberOfResidues();
			   mmdb::Residue *residue_p;
			   mmdb::Atom *residue_at;
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
			   catch (const std::runtime_error &rte) {
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
std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
coot::flev_attached_hydrogens_t::named_hydrogens_to_reference_ligand(mmdb::Residue *ligand_residue_3d,
								     const coot::dictionary_residue_restraints_t &restraints) const {

   std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > > v;

   for (unsigned int i=0; i<named_torsions.size(); i++) { 
      if (named_torsions[i].hydrogen_type == coot::H_IS_RIDING) {
	 mmdb::Atom *atom_base = NULL;
	 mmdb::Atom *atom_2 = NULL;
	 mmdb::Atom *atom_bonded_to_H = NULL;
	 
	 mmdb::PPAtom residue_atoms = 0;
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
	       std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > p(atom_bonded_to_H, cov);
	       v.push_back(p);
	    } 
	 }
      }
   }

   return v;
}
								     

// apply those cannonball directions onto the real reference ligand:
void
coot::flev_attached_hydrogens_t::distances_to_protein(mmdb::Residue *residue_reference, 
						      mmdb::Manager *mol_reference) {

   float radius = 6.0;
   std::vector<mmdb::Residue *> env_residues =
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

	    if (debug)
	       std::cout << "hydrogen on " << named_torsions[i].atom_name_bonded_to_H 
			 << " is riding " << std::endl;
	    
	    std::pair<clipper::Coord_orth, clipper::Coord_orth> pt_base_and_H =
	       hydrogen_pos(named_torsions[i], residue_reference);

	    if (debug)
	       std::cout << "pt_base_and_H: ligand atom at: " << pt_base_and_H.first.format()
			 << " H atom at: " << pt_base_and_H.second.format() << std::endl;
	    
	    std::vector<mmdb::Atom *> atoms = close_atoms(pt_base_and_H.second, env_residues);

	    // bash is one of a (potentially) number of bash distances
	    // for a given ligand (non-hydrogen) atom.
	    // (passing the H coord, the lig-atom coord, a vector of mmdb::Atoms *s.)
	    // 
	    coot::bash_distance_t bash = find_bash_distance(pt_base_and_H.first,
							    pt_base_and_H.second,
							    atoms);
	    if (debug)
	       std::cout << "   found bash: " << bash << std::endl;
	    atom_bashes[named_torsions[i].atom_name_bonded_to_H].push_back(bash);
	 }
	 
	 if (named_torsions[i].hydrogen_type == coot::H_IS_ROTATABLE) {

	    if (debug)
	       std::cout << "hydrogen on " << named_torsions[i].atom_name_bonded_to_H 
			 << " is rotatable " << std::endl;
	    
	    int n_torsion_samples = 8;
	    std::pair<clipper::Coord_orth, clipper::Coord_orth> pt_base_and_H =
	       hydrogen_pos(named_torsions[i], residue_reference);
	    std::vector<mmdb::Atom *> atoms = close_atoms(pt_base_and_H.second, env_residues);
	    for (int itor=0; itor<n_torsion_samples; itor++) {

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
      catch (const std::runtime_error &rte) {
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
coot::flev_attached_hydrogens_t::distances_to_protein_using_correct_Hs(mmdb::Residue *ligand_residue,
								       mmdb::Manager *mol,
								       const coot::protein_geometry &geom) {

   // the constructor (called just before this) should fill
   // atoms_with_rotating_hydrogens and atoms_with_riding_hydrogens
   // vectors (using the restraints).
   
   float radius = 6.0;
   std::vector<mmdb::Residue *> env_residues =
      coot::residues_near_residue(ligand_residue, mol, radius);

   // -------------------------------------------------------------------
   //                    riding hydrogens
   // -------------------------------------------------------------------
   // 
   mmdb::PPAtom residue_atoms = 0;
   int n_ligand_atoms;
   ligand_residue->GetAtomTable(residue_atoms, n_ligand_atoms);
   for (unsigned int irh=0; irh<atoms_with_riding_hydrogens.size(); irh++) { 
      mmdb::Atom *lig_at = NULL;
      mmdb::Atom *H_at = NULL;
      for (int iat=0; iat<n_ligand_atoms; iat++) { 
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
	 
	 std::vector<mmdb::Atom *> atoms = close_atoms(H_pt, env_residues);
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
      mmdb::Atom *lig_at = NULL;
      mmdb::Atom *H_at = NULL;
      for (int iat=0; iat<n_ligand_atoms; iat++) { 
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
	 
	 std::vector<mmdb::Atom *> atoms = close_atoms(H_pt, env_residues);

	 try { 
	    clipper::Coord_orth vector_pt = get_atom_pos_bonded_to_atom(lig_at, H_at, // not H_at
									ligand_residue, geom);
	    clipper::Coord_orth base_ref_pt(0,0,0);
	    double tors = clipper::Coord_orth::torsion(base_ref_pt, vector_pt, lig_atom_pt, H_pt);
	    double dist = sqrt((lig_atom_pt - H_pt).lengthsq());
	    double angle = clipper::Coord_orth::angle(vector_pt, lig_atom_pt, H_pt);

	    int n_torsion_samples = 8;
	    for (int itor=0; itor<n_torsion_samples; itor++) {

	       double tmp_tor_d =  double(itor) * 360.0/double(n_torsion_samples);
	       double tmp_tor = clipper::Util::d2rad(tmp_tor_d);
	       tmp_tor += tors;
	       clipper::Coord_orth new_pt =
		  clipper::Coord_orth(base_ref_pt, vector_pt, lig_atom_pt, dist, angle, tmp_tor);
	       coot::bash_distance_t bash = find_bash_distance(lig_atom_pt, new_pt, atoms);
	       atom_bashes[atoms_with_rotating_hydrogens[irh].first].push_back(bash);
	    }
	 }
	 catch (const std::runtime_error &rte) {
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
					      mmdb::Residue *residue_p) const {
   clipper:: Coord_orth pt(0,0,0);
   clipper:: Coord_orth pt_ligand_atom(0,0,0);
   mmdb::Atom *at_1 = NULL;
   mmdb::Atom *at_2 = NULL;
   mmdb::Atom *at_3 = NULL;

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
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
coot::flev_attached_hydrogens_t::get_atom_pos_bonded_to_atom(mmdb::Atom *lig_at, mmdb::Atom *H_at, // not H_at
							     mmdb::Residue *ligand_residue,
							     const coot::protein_geometry &geom) const {
   int imol = 0; // FIXME needs checking
   std::string res_name(lig_at->residue->GetResName());
   std::pair<bool, dictionary_residue_restraints_t> p = 
      geom.get_monomer_restraints_at_least_minimal(res_name, imol);

   if (! p.first) {
      std::string m = "No monomer type ";
      m += res_name;
      m += " found in dictionary";
      throw std::runtime_error(m);
   } else {
      mmdb::Atom *bonded_atom = NULL;
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
	 mmdb::PPAtom residue_atoms = 0;
	 int n_residue_atoms;
	 ligand_residue->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms ; iat++) {
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
std::vector<mmdb::Atom *>
coot::flev_attached_hydrogens_t::close_atoms(const clipper::Coord_orth &pt,
					     const std::vector<mmdb::Residue *> &env_residues) const {

   std::vector<mmdb::Atom *> v;
   double dist_crit = 6.0;
   double dist_crit_squared = dist_crit * dist_crit;
   
   for (unsigned int i=0; i<env_residues.size(); i++) {
      mmdb::Residue *residue_p = env_residues[i];
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
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
						    const std::vector<mmdb::Atom *> &close_residue_atoms) const {

   // find the residue from the close residue atoms and cache the
   // dictionaries here so that we can then call
   // dictionary_map[residue_p].type_energy(atom_name) and use that
   // type energy (checking for not "") to find if the atom is a
   // hydrogen bond accetor (or both) to use
   // energy_lib_t::some_new_accessor_function_hb_type(type_energy).
   // If hb_type is acceptor, then decrease bash distance,
   // atom_radius_plus_cbr by 0.8A or so.
   // 
   std::map<mmdb::Residue *, coot::dictionary_residue_restraints_t> dictionary_map;

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
