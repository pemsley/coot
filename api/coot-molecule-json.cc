
#include <clipper/core/coords.h>

#include "coot-molecule.hh"
#include "coot-utils/json.hpp"
#include "utils/coot-utils.hh"
#include "rama-plot-phi-psi.hh"
#include "ligand/primitive-chi-angles.hh"

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

std::string
coot::molecule_t::get_torsions_for_residues_in_chain_as_json(const std::string &chain_id) const {

   auto residue_has_alt_confs = [](mmdb::Residue *r) -> bool {
      int n_atoms = r->GetNumberOfAtoms();
      for (int i=0; i<n_atoms; i++) {
         mmdb::Atom *at = r->GetAtom(i);
         if (at->isTer()) continue;
         if (std::string(at->altLoc) != "")
            return true;
      }
      return false;
   };

   auto get_backbone_coord = [](mmdb::Residue *r, const char *name) -> std::pair<bool, clipper::Coord_orth> {
      int n_atoms = r->GetNumberOfAtoms();
      for (int i=0; i<n_atoms; i++) {
         mmdb::Atom *at = r->GetAtom(i);
         if (!at->isTer()) {
            if (std::string(at->name) == name)
               return {true, clipper::Coord_orth(at->x, at->y, at->z)};
         }
      }
      return {false, clipper::Coord_orth()};
   };

   nlohmann::json j_residues = nlohmann::json::array();
   if (!atom_sel.mol) return j_residues.dump();

   for (int imod=1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (!model_p) continue;
      int n_chains = model_p->GetNumberOfChains();
      for (int ic=0; ic<n_chains; ic++) {
         mmdb::Chain *chain_p = model_p->GetChain(ic);
         if (std::string(chain_p->GetChainID()) != chain_id) continue;
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=1; ires<(nres-1); ires++) {
            mmdb::Residue *res_prev = chain_p->GetResidue(ires-1);
            mmdb::Residue *res_this = chain_p->GetResidue(ires);
            mmdb::Residue *res_next = chain_p->GetResidue(ires+1);
            if (!res_prev || !res_this || !res_next) continue;

            if (residue_has_alt_confs(res_this)) continue;

            std::pair<bool, rama_plot::phi_psi_t> pp =
               rama_plot::util::get_phi_psi(res_prev, res_this, res_next);
            if (!pp.first) continue;

            auto [ok_n,  n_pos]  = get_backbone_coord(res_this, " N  ");
            auto [ok_ca, ca_pos] = get_backbone_coord(res_this, " CA ");
            auto [ok_c,  c_pos]  = get_backbone_coord(res_this, " C  ");
            if (!ok_n || !ok_ca || !ok_c) continue;
            double tau = clipper::Util::rad2d(clipper::Coord_orth::angle(n_pos, ca_pos, c_pos));

            nlohmann::json j_chis = nlohmann::json::array();
            try {
               coot::primitive_chi_angles pca(res_this);
               std::vector<coot::alt_confed_chi_angles> chi_info = pca.get_chi_angles();
               if (chi_info.size() == 1) {
                  for (const auto &chi_pair : chi_info[0].chi_angles) {
                     nlohmann::json j_chi;
                     j_chi["chi"] = chi_pair.first;
                     j_chi["value"] = chi_pair.second;
                     j_chis.push_back(j_chi);
                  }
               }
            } catch (const std::runtime_error &) {
               // GLY, ALA etc. have no chi angles
            }

            nlohmann::json j_res;
            j_res["residue_type"] = std::string(res_this->GetResName());
            j_res["res_no"] = res_this->GetSeqNum();
            j_res["chain_id"] = chain_id;
            j_res["ins_code"] = std::string(res_this->GetInsCode());
            j_res["phi"] = pp.second.phi;
            j_res["psi"] = pp.second.psi;
            j_res["tau"] = tau;
            j_res["chi_angles"] = j_chis;

            j_residues.push_back(j_res);
         }
      }
   }
   return j_residues.dump();
}

