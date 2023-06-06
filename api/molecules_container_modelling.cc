
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/atom-selection-container.hh"
#include "molecules_container.hh"

//! return the new molecule number (or -1 on no atoms selected)
int
molecules_container_t::copy_fragment_using_cid(int imol, const std::string &cid) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = get_mol(imol);
      int selHnd = mol->NewSelection(); // d
      mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      mmdb::Manager *new_manager = coot::util::create_mmdbmanager_from_atom_selection(mol, selHnd);
      if (new_manager) {
         imol_new = molecules.size();
         atom_selection_container_t asc = make_asc(new_manager);
         std::string new_name = "copy-fragment-from-molecule-" + std::to_string(imol);
         molecules.push_back(coot::molecule_t(asc, imol_new, new_name));
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

      mmdb::Manager *new_manager = coot::util::create_mmdbmanager_from_atom_selection(mol, selHnd);
      if (new_manager) {

         // create a new molecule
         imol_new = molecules.size();
         atom_selection_container_t asc = make_asc(new_manager);
         std::string new_name = "copy-fragment-from-molecule-" + std::to_string(imol);
         molecules.push_back(coot::molecule_t(asc, imol_new, new_name));
         auto &coot_mol = molecules[imol_new];
         // add a new molecule of the neigbhoring residues internal to the new molecule
         coot_mol.add_neighbor_residues_for_refinement_help(mol);
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
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].cis_trans_conversion(atom_cid, standard_residues_mol);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
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

                  std::cout << "llllllllllllllllllllllllll in add_compound() with merge_results size "
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

         mmdb::Manager *mol_ref = molecules[imol_reference].atom_sel.mol;
         int SelHnd = mol_ref->NewSelection(); // d
         mol_ref->Select(SelHnd, mmdb::STYPE_ATOM, atom_selection.c_str(), mmdb::SKEY_NEW);
         mmdb::Manager *mol_select = coot::util::create_mmdbmanager_from_atom_selection(mol_ref, SelHnd);
         atom_selection_container_t asc_moving = make_asc(mol_select);
         status = molecules[imol_base].replace_fragment(asc_moving);
         mol_ref->DeleteSelection(SelHnd);
         asc_moving.clear_up(); // deletes mol_select
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

