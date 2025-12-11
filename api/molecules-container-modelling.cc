
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/atom-selection-container.hh"
#include "molecules-container.hh"

#include "utils/logging.hh"
extern logging logger;


//! Copy the molecule
//!
//! @param imol the specified molecule
//! @return the new molecule number
int
molecules_container_t::copy_molecule(int imol) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      imol_new = molecules.size();
      mmdb::Manager *mol = coot::util::copy_molecule(molecules[imol].atom_sel.mol);
      atom_selection_container_t asc = make_asc(mol);
      std::string new_name = "copy-of-molecule-" + std::to_string(imol);
      molecules.push_back(coot::molecule_t(asc, imol_new, new_name));
   }

   if (is_valid_map_molecule(imol)) {
      imol_new = molecules.size();
      std::string new_name = "copy-of-molecule-" + std::to_string(imol);
      const clipper::Xmap<float> &xmap = molecules[imol].xmap;
      bool is_em = molecules[imol].is_EM_map();
      molecules.push_back(coot::molecule_t(new_name, imol_new, xmap, is_em));
   }
   return imol_new;
}



//! return the new molecule number (or -1 on no atoms selected)
int
molecules_container_t::copy_fragment_using_cid(int imol, const std::string &multi_cids) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = get_mol(imol);
      int selHnd = mol->NewSelection(); // d
      std::vector<std::string> v = coot::util::split_string(multi_cids, "||");
      if (! v.empty())
         for (const auto &cid : v)
            mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_OR);
      // replace_fragment() uses UDDOldAtomIndexHandle for fast indexing.
      mmdb::Manager *new_manager = coot::util::create_mmdbmanager_from_atom_selection(mol, selHnd);
      if (new_manager) {
         int transfer_atom_index_handle = new_manager->GetUDDHandle(mmdb::UDR_ATOM, "transfer atom index");

         // std::cout << "..... transfer_atom_index_handle A " << transfer_atom_index_handle << std::endl;

         imol_new = molecules.size();
         atom_selection_container_t asc = make_asc(new_manager);
         asc.UDDOldAtomIndexHandle = transfer_atom_index_handle;
         std::string new_name = "copy-fragment-from-molecule-" + std::to_string(imol);
         molecules.push_back(coot::molecule_t(asc, imol_new, new_name));
         if (false)
            std::cout << "debug:: in mc::copy_fragment_using_cid(): the UDDOldAtomIndexHandle for molecule " << imol_new
                      << " is " << molecules[imol_new].atom_sel.UDDOldAtomIndexHandle << std::endl;

      }
      mol->DeleteSelection(selHnd);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return imol_new;
}

//! copy a fragment - use this in preference to `copy_fragment_using_cid()` when copying
//! a molecule fragment to make a molten zone for refinement.
//! That is because this version quietly also copies the residues near the residues of the selection.
//! so that those residues can be used for links and non-bonded contact restraints.
//! @return the new molecule number (or -1 on no atoms selected)
int
molecules_container_t::copy_fragment_for_refinement_using_cid(int imol, const std::string &multi_cids) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = get_mol(imol);
      int selHnd = mol->NewSelection(); // d

      // 20230606-PE old -  simple
      // mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);

      std::vector<std::string> v = coot::util::split_string(multi_cids, "||");
      if (! v.empty())
         for (const auto &cid : v)
            mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_OR);

      // replace_fragment() uses UDDOldAtomIndexHandle for fast indexing.
      mmdb::Manager *new_manager = coot::util::create_mmdbmanager_from_atom_selection(mol, selHnd);
      if (new_manager) {

         int transfer_atom_index_handle = new_manager->GetUDDHandle(mmdb::UDR_ATOM, "transfer atom index");
         // create a new molecule
         imol_new = molecules.size();
         atom_selection_container_t asc = make_asc(new_manager);
         asc.UDDOldAtomIndexHandle = transfer_atom_index_handle;
         std::string new_name = "copy-fragment-from-molecule-" + std::to_string(imol);
         molecules.push_back(coot::molecule_t(asc, imol_new, new_name));
         auto &coot_mol = molecules[imol_new];
         // add a new molecule of the neighboring residues internal to the new molecule
         coot_mol.add_neighbor_residues_for_refinement_help(mol);
      } else {
         std::cout << "WARNING:: copy_fragment_for_refinement_using_cid() new_manager was null" << std::endl;
      }
      mol->DeleteSelection(selHnd);
   }
   return imol_new;
 }



//! return the new molecule number (or -1 on no atoms selected)

int
molecules_container_t::copy_fragment_using_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = get_mol(imol);
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         std::vector<mmdb::Residue *> selected_residues;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id_this(chain_p->GetChainID());
            if (chain_id_this == chain_id) {
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int res_no_this = residue_p->GetSeqNum();
                     if (res_no_this <= res_no_end) {
                        if (res_no_this >= res_no_start) {
                           selected_residues.push_back(residue_p);
                        }
                     }
                  }
               }
            }
         }
         if (! selected_residues.empty()) {
            std::pair<bool,std::string> use_alt_conf = std::pair<bool, std::string>(false, "");
            std::pair<bool, mmdb::Manager *> new_manager =
               coot::util::create_mmdbmanager_from_residue_vector(selected_residues, mol, use_alt_conf);
            if (new_manager.first) {
               imol_new = molecules.size();
               atom_selection_container_t asc = make_asc(new_manager.second);
               std::string new_name = "atom-selection-from-molecule-" + std::to_string(imol);
               molecules.push_back(coot::molecule_t(asc, imol_new, new_name));
            }
         }
      }
   }
   return imol_new;
}


void
molecules_container_t::eigen_flip_ligand(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      coot::minimol::molecule mm = molecules[imol].eigen_flip_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

}


void
molecules_container_t::eigen_flip_ligand_using_cid(int imol, const std::string &residue_cid) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_cid_to_residue_spec(imol, residue_cid);
      coot::minimol::molecule mm = molecules[imol].eigen_flip_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
}


std::pair<int, std::vector<merge_molecule_results_info_t> >
molecules_container_t::merge_molecules(int imol, const std::string &list_of_other_molecules) {

   int istat = 0;
   std::vector<merge_molecule_results_info_t> resulting_merge_info;
   if (is_valid_model_molecule(imol)) {

      std::vector<atom_selection_container_t> atom_selections;
      std::vector<std::string> number_strings = coot::util::split_string(list_of_other_molecules, ":");
      for (const auto &item : number_strings) {
         int idx = coot::util::string_to_int(item);
         if (is_valid_model_molecule(idx))
            atom_selections.push_back(molecules[idx].atom_sel);
      }
      auto r = molecules[imol].merge_molecules(atom_selections);
      istat = r.first;
      resulting_merge_info = r.second;
      set_updating_maps_need_an_update(imol);
   }
   return std::pair<int, std::vector<merge_molecule_results_info_t> > (istat, resulting_merge_info);

}

std::pair<int, std::vector<merge_molecule_results_info_t> >
molecules_container_t::merge_molecules(int imol, std::vector<mmdb::Manager *> mols) {

   int istat = 0;
   std::vector<merge_molecule_results_info_t> resulting_merge_info;
   if (is_valid_model_molecule(imol)) {

      std::vector<atom_selection_container_t> atom_selections;
      for (const auto &mol : mols) {
         atom_selection_container_t atom_sel = make_asc(mol);
         atom_selections.push_back(atom_sel);
      }
      auto r = molecules[imol].merge_molecules(atom_selections);
      istat = r.first;
      resulting_merge_info = r.second;
   }
   return std::pair<int, std::vector<merge_molecule_results_info_t> > (istat, resulting_merge_info);

}


int
molecules_container_t::cis_trans_convert(int imol, const std::string &atom_cid) {

   int status = 0;
   mmdb::Manager *standard_residues_mol = standard_residues_asc.mol;
   if (standard_residues_mol) {
      if (is_valid_model_molecule(imol)) {
         status = molecules[imol].cis_trans_conversion(atom_cid, standard_residues_mol);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
      }
   } else {
      logger.log(log_t::ERROR, logging::function_name_t("mc::cis_trans_convert"),
                 "Null standard_residues_asc.mol");
   }
   return status;

}


//! This function is for adding compounds/molecules like buffer agents and precipitants or anions and cations.
//! i.e. those ligands that can be positioned without need for internal torsion angle manipulation.
//! @return the success status, 1 or good, 0 for not goo.
int
molecules_container_t::add_compound(int imol, const std::string &tlc, int imol_dict, int imol_map, float x, float y, float z) {

   auto move_residue = [] (mmdb::Residue *residue_p, const coot::Cartesian &position) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      coot::Cartesian position_sum;
      unsigned int n_atoms = 0;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            n_atoms++;
            position_sum + coot::Cartesian(at->x, at->y, at->z);
         }
      }
      if (n_atoms > 0) {
         float inv = 1.0 / static_cast<float>(n_atoms);
         coot::Cartesian current_position(position_sum * inv);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               at->x += position.x() - current_position.x();
               at->y += position.y() - current_position.y();
               at->z += position.z() - current_position.z();
            }
         }
      }
   };

   auto make_residue_spec_from_first_residue = [] (mmdb::Manager *mol, const std::string &chain_id) {
      coot::residue_spec_t spec;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (std::string(chain_p->GetChainID()) == chain_id) {
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     spec = coot::residue_spec_t(residue_p);
                  }
               }
            }
         }
      }
      return spec;
   };

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         coot::Cartesian position(x,y,z);
         get_monomer(tlc); // just load the dictionary
         auto mr = geom.get_monomer_restraints(tlc, imol_dict);
         std::cout << "in add_compound() mr.first is " << mr.first << std::endl;
         if (mr.first) {
            const auto &monomer_restraints = mr.second;
            bool idealised_flag = true;
            float b_factor = coot::util::median_temperature_factor(molecules[imol].atom_sel.atom_selection,
                                                                   molecules[imol].atom_sel.n_selected_atoms,
                                                                   0, 100, false, false);
            mmdb::Residue *residue_p = monomer_restraints.GetResidue(idealised_flag, b_factor);
            if (residue_p) {
               mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
               if (mol) {
                  move_residue(residue_p, position);
                  const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
                  float map_rmsd = molecules[imol_map].get_map_rmsd_approx();
                  std::vector<mmdb::Manager *> mols = { mol };
                  auto merge_results = merge_molecules(imol, mols);

                  std::cout << "debug:: in add_compound() with merge_results size "
                            << merge_results.first << " "
                            << merge_results.second.size() << std::endl;

                  if (merge_results.first == 1) {
                     if (merge_results.second.size() == 1) {
                        coot::residue_spec_t res_spec = merge_results.second[0].spec;

                        std::cout << "in add_compound() is_chain is " << merge_results.second[0].is_chain << std::endl;
                        std::cout << "in add_compound() chain_id is " << merge_results.second[0].chain_id << std::endl;
                        std::cout << "in add_compound() res_spec is " << res_spec << std::endl;

                        if (merge_results.second[0].is_chain) {
                           // set the res-spec to be the first residue of the chain
                           res_spec = make_residue_spec_from_first_residue(get_mol(imol), merge_results.second[0].chain_id);
                           std::cout << "in add_compound() updated res_spec is " << res_spec << std::endl;
                        }

                        int n_trials = 100;
                        float translation_scale_factor = 1.0;
                        float d = molecules[imol].fit_to_map_by_random_jiggle(res_spec, xmap, map_rmsd, n_trials, translation_scale_factor);
                        std::cout << "score from fit_to_map_by_random_jiggle() " << d << std::endl;
                        set_updating_maps_need_an_update(imol);
                        status = 1;
                     }
                  }
               }
            }
         } else {
            std::cout << "WARNING:: Not compound " << tlc << " found in dictionary store" << std::endl;
         }
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


//! replace a fragment
//!
//! _i.e._ replace the atoms of ``imol_base`` by those of the atom selection ``atom_selection`` in ``imol_reference``
//! (``imol_base`` is the molecule that is modified).
//!
//! @return the success status
int
molecules_container_t::replace_fragment(int imol_base, int imol_reference, const std::string &atom_selection) {

   int status = 0;
   if (is_valid_model_molecule(imol_base)) {
      if (is_valid_model_molecule(imol_reference)) {

         std::string multi_cids = atom_selection;
         mmdb::Manager *mol_ref = molecules[imol_reference].atom_sel.mol;
         int old_atom_index_handle = molecules[imol_reference].atom_sel.UDDOldAtomIndexHandle;
         int SelHnd = mol_ref->NewSelection(); // d
         // 20231113-PE replace with multi-selection
         // mol_ref->Select(SelHnd, mmdb::STYPE_ATOM, atom_selection.c_str(), mmdb::SKEY_NEW);
         std::vector<std::string> v = coot::util::split_string(multi_cids, "||");
         if (! v.empty())
            for (const auto &cid : v)
               mol_ref->Select(SelHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_OR);

#if 0 // 20240201-PE replace by creating a new molecule
         mmdb::Manager *mol_select = coot::util::create_mmdbmanager_from_atom_selection(mol_ref, SelHnd);
         atom_selection_container_t asc_moving = make_asc(mol_select);
         status = molecules[imol_base].replace_fragment(asc_moving);
#endif

         status = molecules[imol_base].replace_fragment(mol_ref, old_atom_index_handle, SelHnd);

         mol_ref->DeleteSelection(SelHnd);

#if 0 // 20240201-PE now replaced
         asc_moving.clear_up(); // deletes mol_select
#endif
         set_updating_maps_need_an_update(imol_base);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_reference << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_base << std::endl;
   }
   return status;
}


//! Rigid-body fitting
//!
//! `multi_cids" is a "||"-separated list of residues CIDs, e.g. "//A/12-52||//A/14-15||/B/56-66"
int
molecules_container_t::rigid_body_fit(int imol, const std::string &multi_cid, int imol_map) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         status = molecules[imol].rigid_body_fit(multi_cid, xmap);
      } else {
         std::cout << "ERROR:: in rigid_body_fit() bad map index " << imol_map << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;


}


// minimize/optimize the geometry of the specified residue(s)
// @return the success status 1 if the minimization was performed and 0 if it was not.
std::pair<int, coot::instanced_mesh_t>
molecules_container_t::minimize_energy(int imol, const std::string &atom_selection_cid,
                                       int n_cycles,
                                       bool do_rama_plot_restraints, float rama_plot_weight,
                                       bool do_torsion_restraints, float torsion_weight, bool refinement_is_quiet) {
   int status = 0;
   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].minimize(atom_selection_cid, n_cycles,
                                        do_rama_plot_restraints, rama_plot_weight,
                                        do_torsion_restraints, torsion_weight,
                                        refinement_is_quiet, &geom);
      std::string mode = "COLOUR-BY-CHAIN-AND-DICTIONARY";
      bool draw_hydrogen_atoms_flag = true; // use data member as we do for draw_missing_residue_loops_flag?
      unsigned int smoothness_factor = 1;
      bool show_atoms_as_aniso_flag = false;
      bool show_aniso_atoms_as_ortep = false;
      float aniso_probability = 0.5f;
      im = molecules[imol].get_bonds_mesh_instanced(mode, &geom, true, 0.12, 1.4,
                                                    show_atoms_as_aniso_flag,
                                                    aniso_probability,
                                                    show_aniso_atoms_as_ortep,
                                                    smoothness_factor,
                                                    draw_hydrogen_atoms_flag, draw_missing_residue_loops_flag);

   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   // use std::move() here?
   return std::make_pair(status, im);;

}

//! @param imol is the model molecule index
//! @param atom_selection_cid is the selection CID e.g. "//A/15" (residue 15 of chain A)
//! @param n_cycles is the number of refinement cycles. If you pass n_cycles = 100 (or some such) then you can
//!         get the mesh for the partially optimized ligand/residues
//! @param do_rama_plot_restraints is the flag for the usage of Ramachandran plot restraints
//! @param rama_plot_weight is the flag to set the Ramachandran plot restraints weight
//! @param do_torsion_restraints is the flag for the usage of torsion restraints
//! @param torsion_weight is the flag to set the torsion restraints weight
//! @param refinement_is_quiet is used to reduce the amount of diagnostic text written to the output
//!
//! @return the function value at termination
float
molecules_container_t::minimize(int imol, const std::string &atom_selection_cid,
                                int n_cycles,
                                bool do_rama_plot_restraints, float rama_plot_weight,
                                bool do_torsion_restraints, float torsion_weight, bool refinement_is_quiet) {

   int status = 0;
   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].minimize(atom_selection_cid, n_cycles,
                                        do_rama_plot_restraints, rama_plot_weight,
                                        do_torsion_restraints, torsion_weight,
                                        refinement_is_quiet, &geom);

   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   return status; // this is the wrong return value.
}



coot::graph_match_info_t
molecules_container_t::overlap_ligands_internal(int imol_ligand, int imol_ref, const std::string &chain_id_ref,
                                                int resno_ref, bool apply_rtop_flag) {

   coot::graph_match_info_t graph_info;

   mmdb::Residue *residue_moving = 0;
   mmdb::Residue *residue_reference = 0;

   // running best ligands:
   mmdb::Residue *best_residue_moving = nullptr;
   double best_score = 99999999.9; // low score good.
   clipper::RTop_orth best_rtop;


   if (! is_valid_model_molecule(imol_ligand))
      return graph_info;

   if (! is_valid_model_molecule(imol_ref))
      return graph_info;

   mmdb::Manager *mol_moving = molecules[imol_ligand].atom_sel.mol;
   mmdb::Manager *mol_ref    = molecules[imol_ref].atom_sel.mol;

   for (int imod=1; imod<=mol_moving->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol_moving->GetModel(imod);
      int nchains_moving = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains_moving; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::PResidue residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               int n_atoms = residue_p->GetNumberOfAtoms();
               if (n_atoms > 0) {
                  residue_moving = residue_p;
                  break;
               }
            }
         }
         if (residue_moving)
            break;
      }

      if (! residue_moving) {
         std::cout << "Oops.  overlap_ligands_internal()): Failed to find moving residue" << std::endl;
      } else {
         int imodel_ref = 1;
         mmdb::Model *model_ref_p = mol_ref->GetModel(imodel_ref);
         mmdb::Chain *chain_p;
         int nchains = model_ref_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_ref_p->GetChain(ichain);
            if (std::string(chain_p->GetChainID()) == std::string(chain_id_ref)) {
               int nres = chain_p->GetNumberOfResidues();
               mmdb::PResidue residue_p;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int seqnum = residue_p->GetSeqNum();
                     if (seqnum == resno_ref) {
                        residue_reference = residue_p;
                        break;
                     }
                  }
               }
               if (residue_reference)
                  break;
            }
         }

         if (!residue_reference) {
            std::cout << "Oops.  Failed to find reference residue" << std::endl;
         } else {
            bool match_hydrogens_also = 0;
            coot::graph_match_info_t rtop_info =
               coot::graph_match(residue_moving, residue_reference, apply_rtop_flag, match_hydrogens_also);
            if (rtop_info.success) {
               if (rtop_info.dist_score < best_score) { // low score good.
                  best_score = rtop_info.dist_score;
                  best_residue_moving = residue_moving;
                  best_rtop = rtop_info.rtop;
                  graph_info = rtop_info;
               }
            } else {
               // std::cout << "Oops.  Match failed somehow" << std::endl;
            }
         }
      }
   }
   if (apply_rtop_flag) {
      if (best_residue_moving) {
         // move just the best ligand:
         molecules[imol_ligand].transform_by(best_rtop, best_residue_moving);
         // delete everything except the best_residue_moving
         // molecules[imol_ligand].delete_all_except_res(best_residue_moving);
      }
   }
   return graph_info;
}


//! match ligand torsions and positions, different api
bool
molecules_container_t::match_ligand_torsions_and_position_using_cid(int imol_ligand, int imol_ref, const std::string &cid) {

   bool status = false;
   if (is_valid_model_molecule(imol_ligand)) {
      if (is_valid_model_molecule(imol_ref)) {
         std::pair<bool, coot::residue_spec_t> rs_pair = molecules[imol_ref].cid_to_residue_spec(cid);
         if (rs_pair.first) {
            status = match_ligand_torsions_and_position(imol_ligand, imol_ref, rs_pair.second.chain_id, rs_pair.second.res_no);
         }
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_ref << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_ligand << std::endl;
   }
   return status;

}

//! match ligand torsions - return the success status
bool
molecules_container_t::match_ligand_torsions(int imol_ligand, int imol_ref, const std::string &chain_id_ref, int resno_ref) {

   bool status = false;
   if (is_valid_model_molecule(imol_ligand)) {
      if (is_valid_model_molecule(imol_ref)) {
         coot::residue_spec_t res_spec(chain_id_ref, resno_ref, "");
         mmdb::Residue *res_reference = molecules[imol_ref].get_residue(res_spec);
         if (res_reference) {
            std::string res_name_ref_res(res_reference->GetResName());
            std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
               geom.get_monomer_restraints(res_name_ref_res, imol_ref);
            if (restraints_info.first) {
               std::vector <coot::dict_torsion_restraint_t> tr_ref_res =
                  geom.get_monomer_torsions_from_geometry(res_name_ref_res, 0);
               int n_torsions_moved = molecules[imol_ligand].match_torsions(res_reference, tr_ref_res, geom);

               if (n_torsions_moved > 0) status = true;

               // currently, this doesn't make sense, but it will when imol_ligand contains a whole protein
               set_updating_maps_need_an_update(imol_ligand);
            }
         }
      }
   }
   return status;
}

//! match ligand positions - return the success status
bool
molecules_container_t::match_ligand_position(int imol_ligand, int imol_ref, const std::string &chain_id_ref, int resno_ref) {
   bool status = false;
   if (is_valid_model_molecule(imol_ligand)) {
      if (is_valid_model_molecule(imol_ref)) {
         bool apply_rtop_flag = true;
         overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, apply_rtop_flag);
         status = true;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_ligand << std::endl;
   }
   return status;
}

//! match ligand torsions and positions
bool
molecules_container_t::match_ligand_torsions_and_position(int imol_ligand, int imol_ref, const std::string &chain_id_ref, int resno_ref) {
   bool status = false;
   if (is_valid_model_molecule(imol_ligand)) {
      if (is_valid_model_molecule(imol_ref)) {
         coot::residue_spec_t res_spec(chain_id_ref, resno_ref, "");
         status = match_ligand_torsions(imol_ligand, imol_ref, chain_id_ref, resno_ref);
         match_ligand_position(imol_ligand, imol_ref, chain_id_ref, resno_ref);
         set_updating_maps_need_an_update(imol_ligand);
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_ligand << std::endl;
   }
   return status;
}

void 
molecules_container_t::associate_sequence(int imol, const std::string &chain_id, const std::string &sequence) {

   if (is_valid_model_molecule(imol)) { // but also maps?
      molecules[imol].associate_sequence_with_molecule(chain_id, sequence);
   }

}

void
molecules_container_t::assign_sequence(int imol_model, int imol_map) {

   if (is_valid_model_molecule(imol_model)) {
       if (is_valid_map_molecule(imol_map)) {
           const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
           molecules[imol_model].assign_sequence(xmap, geom);
       }
   }
}


//! Get the conformers that can be generated by variation around rotatable bonds as described in the dictionary.
//! Torsions that are marked as "const" are excluded from the variation, as are pyranose ring torsions
//! and torsions that rotate hydrogen atoms.
//! @return a vector of indices of the new molecules
std::vector<int>
molecules_container_t::get_dictionary_conformers(const std::string &comp_id, int imol_enc, bool remove_internal_clash_conformers) {

   std::vector<int> mol_indices;
   std::pair<bool, coot::dictionary_residue_restraints_t> r = geom.get_monomer_restraints(comp_id, imol_enc);
   if (r.first) {
      bool delete_clash_confs = true;
      std::vector<mmdb::Residue *> confs = coot::util::get_dictionary_conformers(r.second, delete_clash_confs);
      for (unsigned int i=0; i<confs.size(); i++) {
         mmdb::Residue *res = confs[i];
         mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(res);
         // std::cout << "conf " << i  << " mol: " << mol << std::endl;
         std::string name = comp_id + "-conf-" + std::to_string(i);
         atom_selection_container_t asc = make_asc(mol);
         int imol_no = molecules.size();
         molecules.push_back(coot::molecule_t(asc, imol_no, name));
         mol_indices.push_back(imol_no);
      }
      for (unsigned int i=0; i<confs.size(); i++)
         delete confs[i];
   }
   if (false) {
      std::cout << "mol_indices:" << std::endl;
      for (unsigned int i=0; i<mol_indices.size(); i++) {
         std::cout << "       " << mol_indices[i] << std::endl;
      }
   }
   return mol_indices;
}

//! replace a residue
//!
//! Change the type of a residue (for example, "TYR" to "PTY").
//! The algorithm will superpose the mainchain CA, C and N and try to set matching torsion
//! to the angle that they were in the reference structure.
void
molecules_container_t::replace_residue(int imol, const std::string &residue_cid,
                                       const std::string &new_residue_type, int imol_enc) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].replace_residue(residue_cid, new_residue_type, imol_enc, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   // return status
}


//! Rotate atoms around torsion
//!
//! the bond is presumed to be between atom-2 and atom-3. Atom-1 and atom-4 are
//! used to define the absolute torsion angle.
//!
//! @return status 1 if successful, 0 if not.
int
molecules_container_t::rotate_around_bond(int imol, const std::string &residue_cid,
                                          const std::string &atom_name_1,
                                          const std::string &atom_name_2,
                                          const std::string &atom_name_3,
                                          const std::string &atom_name_4,
                                          double torsion_angle) {

   std::string alt_conf = "";
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_name_quad quad(atom_name_1, atom_name_2, atom_name_3, atom_name_4);
      status = molecules[imol].rotate_around_bond(residue_cid, alt_conf, quad, torsion_angle, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol
                << std::endl;
   }
   return status;
}




//! change the alt confs
//!
//! @param change_mode is either "residue", "side-chain" or a comma-separated atom-name
//! pairs (e.g "N,CA") - you can (of course) specify just one atom: "N".
// @return the success status (1 is done, 0 means failed to do)
int
molecules_container_t::change_alt_locs(int imol, const std::string &cid, const std::string &change_mode) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].change_alt_locs(cid, change_mode);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


coot::instanced_mesh_t
molecules_container_t::get_HOLE(int imol,
                                float start_pos_x, float start_pos_y, float start_pos_z,
                                float end_pos_x, float end_pos_y, float end_pos_z) const {

   coot::instanced_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      clipper::Coord_orth start_pos(start_pos_x, start_pos_y, start_pos_z);
      clipper::Coord_orth end_pos(end_pos_x, end_pos_y, end_pos_z);
      m = molecules[imol].get_HOLE(start_pos, end_pos, geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return m;

}

//! Get SVG for 2d ligand environment view (FLEV)
//!
//! The caller should make sure that the dictionary for the ligand has been loaded - this
//! function won't do that. It will add hydrogen atoms if needed.
//!
//! @param imol is the model molecule index
//! @param residue_cid is the cid for the residue
std::string
molecules_container_t::get_svg_for_2d_ligand_environment_view(int imol,
                                                              const std::string &residue_cid,
                                                              bool add_key) {
   float radius = 4.2; // pass this
   std::string s;
   if (is_valid_model_molecule(imol)) {
      s = molecules[imol].get_svg_for_2d_ligand_environment_view(residue_cid, &geom, add_key);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return s;
}



//! delete all carbohydrate
//!
//! @param imol is the model molecule index
//!
//! @return true on successful deletion, return false on no deletion.
bool
molecules_container_t::delete_all_carbohydrate(int imol) {

   bool status = false;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].delete_all_carbohydrate();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}
