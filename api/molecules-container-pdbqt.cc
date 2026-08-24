/*
 * api/molecules-container-pdbqt.cc
 *
 * Thin molecules_container_t shims over the PDBQT writers. The receptor writer
 * lives in coot-utils/pdbqt.hh (RDKit-free) and the flexible-ligand writer in
 * lidia-core/pdbqt-ligand.hh, so the same functionality is reusable by the GUI
 * (src/) without depending on api/.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#include <iostream>

#include "molecules-container.hh"
#include "coot-utils/pdbqt.hh"
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/pdbqt-ligand.hh"
#endif

int
molecules_container_t::export_molecule_as_pdbqt(int imol, const std::string &file_name) {

   if (! is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: export_molecule_as_pdbqt(): invalid model molecule " << imol << std::endl;
      return 0;
   }
   return coot::pdbqt::write_receptor(molecules[imol].atom_sel.mol, geom, imol, file_name);
}

int
molecules_container_t::export_flexible_receptor_as_pdbqt(int imol, const std::string &flex_residues_cid,
                                                         const std::string &rigid_file_name,
                                                         const std::string &flex_file_name) {

   if (! is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: export_flexible_receptor_as_pdbqt(): invalid model molecule "
                << imol << std::endl;
      return 0;
   }
   std::vector<mmdb::Residue *> residues = molecules[imol].cid_to_residues(flex_residues_cid);
   std::vector<coot::residue_spec_t> flex_specs;
   for (unsigned int i=0; i<residues.size(); i++)
      if (residues[i]) flex_specs.push_back(coot::residue_spec_t(residues[i]));
   if (flex_specs.empty()) {
      std::cout << "WARNING:: export_flexible_receptor_as_pdbqt(): no residues for cid "
                << flex_residues_cid << std::endl;
      return 0;
   }
   return coot::pdbqt::write_flexible_receptor(molecules[imol].atom_sel.mol, geom, imol,
                                               flex_specs, rigid_file_name, flex_file_name);
}

int
molecules_container_t::export_flexible_receptor_near_point_as_pdbqt(int imol,
                                                                    float x, float y, float z,
                                                                    float radius,
                                                                    const std::string &rigid_file_name,
                                                                    const std::string &flex_file_name) {

   if (! is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: export_flexible_receptor_near_point_as_pdbqt(): invalid model molecule "
                << imol << std::endl;
      return 0;
   }
   std::vector<coot::residue_spec_t> flex_specs =
      coot::pdbqt::flexible_side_chain_residues_near(molecules[imol].atom_sel.mol, x, y, z, radius);
   if (flex_specs.empty()) {
      std::cout << "WARNING:: export_flexible_receptor_near_point_as_pdbqt(): no flexible side chains within "
                << radius << " A of point" << std::endl;
      return 0;
   }
   return coot::pdbqt::write_flexible_receptor(molecules[imol].atom_sel.mol, geom, imol,
                                               flex_specs, rigid_file_name, flex_file_name);
}

int
molecules_container_t::read_pdbqt(const std::string &file_name) {

   int imol = -1;
   mmdb::Manager *mol = coot::pdbqt::read(file_name);
   if (mol) {
      imol = molecules.size();
      atom_selection_container_t asc = make_asc(mol);
      molecules.push_back(coot::molecule_t(asc, imol, file_name));
   }
   return imol;
}

std::vector<coot::pdbqt::pose_score_t>
molecules_container_t::get_vina_scores(int imol) const {

   if (! is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: get_vina_scores(): invalid model molecule " << imol << std::endl;
      return std::vector<coot::pdbqt::pose_score_t>();
   }
   return coot::pdbqt::get_scores(molecules[imol].atom_sel.mol);
}

int
molecules_container_t::export_ligand_as_pdbqt(int imol, const std::string &cid,
                                              const std::string &file_name) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   if (! is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: export_ligand_as_pdbqt(): invalid model molecule " << imol << std::endl;
      return 0;
   }
   mmdb::Residue *residue_p = get_residue_using_cid(imol, cid);
   if (! residue_p) {
      std::cout << "WARNING:: export_ligand_as_pdbqt(): no residue for cid " << cid << std::endl;
      return 0;
   }
   std::string res_name(residue_p->GetResName());
   std::pair<bool, coot::dictionary_residue_restraints_t> rp =
      geom.get_monomer_restraints(res_name, imol);
   const coot::dictionary_residue_restraints_t *rest_p = rp.first ? &rp.second : nullptr;
   return coot::pdbqt::write_ligand(residue_p, rest_p, res_name, cid, file_name);
#else
   std::cout << "WARNING:: export_ligand_as_pdbqt(): needs MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
   return 0;
#endif
}
