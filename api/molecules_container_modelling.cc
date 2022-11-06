
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
