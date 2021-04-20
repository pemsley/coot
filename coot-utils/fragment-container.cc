
#include <iostream>
#include "fragment-container.hh"

coot::fragment_container_t
coot::make_fragments(mmdb::Manager *mol) {

   coot::fragment_container_t fc;

   if (! mol) return fc;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (! model_p) return fc;

   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      std::string chain_id(chain_p->GetChainID());
      int n_res = chain_p->GetNumberOfResidues();
      if (n_res > 0) {
         mmdb::Residue *residue_start = 0;
         mmdb::Residue *residue_prev = 0;
         std::vector<mmdb::Residue *> residues_running;
         for (int i_res=0; i_res<n_res; i_res++) {
            mmdb::Residue *res = chain_p->GetResidue(i_res);
            if (! res) continue;
            if (residue_start == 0) {
               residue_start = res;
               residues_running.push_back(res);
            } else {
               if (residue_prev) {
                  int resno_prev = residue_prev->GetSeqNum();
                  int resno_this = res->GetSeqNum();
                  if (resno_prev != (resno_this-1)) {
                     fragment_container_t::fragment_range_t fr(chain_id, residue_spec_t(residue_start), residue_spec_t(residue_prev));
                     std::cout << "transferring  residues running of size " << residues_running.size() << std::endl;
                     fr.residues = residues_running;
                     residue_start = res;
                     fc.add(fr);
                     residues_running.clear();
                  } else {
                     residues_running.push_back(res);
                  }
               }
            }
            residue_prev = res; // for next round
         }

         if (! residues_running.empty()) {
            if (residue_start) {
               if (residue_prev) {
                  if (residue_start != residue_prev) {
                     fragment_container_t::fragment_range_t fr(chain_id, residue_spec_t(residue_start), residue_spec_t(residue_prev));
                     fr.residues = residues_running;
                     std::cout << "transferring  residues running of size " << residues_running.size() << std::endl;
                     fc.add(fr);
                  }
               }
            }
         }
      }
   }

   return fc;
}

void
coot::fragment_container_t::print_fragments() const {

   for (const auto &r : ranges) {
      std::cout << "Fragment: in Chain " << r.chain_id << "  " << r.start_res << " " << r.end_res << std::endl;
   }

}
