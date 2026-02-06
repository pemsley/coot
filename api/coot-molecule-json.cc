
#include "coot-molecule.hh"
#include "coot-utils/json.hpp"
#include "utils/coot-utils.hh"

std::string coot::molecule_t::get_molecule_selection_as_json(const std::string &cid) const {

   auto atom_to_json = [] (mmdb::Atom *at) {
      nlohmann::json j;
      // 2025-10-10-PE add more attributes later
      std::string se = util::remove_whitespace(std::string(at->element));
      std::string sn = util::remove_whitespace(std::string(at->name));
      j["x"] = at->x;
      j["y"] = at->y;
      j["z"] = at->z;
      j["tempFactor"] = at->tempFactor;
      j["occupancy"] = at->occupancy;
      j["name"] = sn;
      j["element"] = se;
      return j;
   };

   std::string s;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Atom **SelAtoms = nullptr;
      int nSelAtoms = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         nlohmann::json j_models = nlohmann::json::array();
         for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               nlohmann::json j_chains = nlohmann::json::array();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_res = chain_p->GetNumberOfResidues();
                  nlohmann::json j_residues = nlohmann::json::array();
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        int n_atoms = residue_p->GetNumberOfAtoms();
                        nlohmann::json j_atoms = nlohmann::json::array();
                        int atom_count = 0;
                        for (int iat=0; iat<n_atoms; iat++) {
                           mmdb::Atom *at = residue_p->GetAtom(iat);
                           if (! at->isTer()) {
                              if (at->isInSelection(selHnd)) {
                                 nlohmann::json j = atom_to_json(at);
                                 j_atoms.push_back(j);
                                 atom_count++;
                              }
                           }
                        }
                        if (atom_count > 0) {
                           nlohmann::json j_residue;
                           j_residue["name"]     = residue_p->GetResName();
                           j_residue["ins_code"] = residue_p->GetInsCode();
                           j_residue["seqnum"]   = residue_p->GetSeqNum();
                           j_residue["atoms"] = j_atoms;
                           j_residues.push_back(j_residue);
                        }
                     }
                  }
                  if (j_residues.size() > 0) {
                     nlohmann::json j_chain;
                     j_chain["chain_id"] = chain_p->GetChainID();
                     j_chain["residues"] = j_residues;
                     j_chains.push_back(j_chain);
                  }
               }
               if (j_chains.size() > 0) {
                  nlohmann::json j_model;
                  j_model["chains"] = j_chains;
                  j_models.push_back(j_model);
               }
            }
         }
         nlohmann::json j_manager;
         j_manager["models"] = j_models;
         s = j_manager.dump(2);
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return s;
}

