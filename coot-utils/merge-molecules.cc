
#include <string>
#include "merge-molecules.hh"


void
coot::merge_molecules(mmdb::Manager *mol_first, std::vector<mmdb::Manager *> mol_others) {

   // only copy the first chain of the mol_others;

   auto index_to_chain_id = [] (unsigned int idx) {
                               unsigned int tens = idx/10;
                               unsigned int units = idx - 10 * tens;
                               std::string abc("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
                               std::string r = "A";
                               char t = abc[tens];
                               std::string rr(1, t);
                               r = rr + std::to_string(units);
                               return rr;
                            };

   mmdb::Model *mol_first_model_p = mol_first->GetModel(1);
   if (! mol_first_model_p) return; //   quietly, is that sensible?

   unsigned int chain_idx = 0;
   for (unsigned int i=0; i<mol_others.size(); i++) {
      mmdb::Manager *mol = mol_others[i];
      mmdb::Model *model_p = mol->GetModel(1);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ich=0; ich<n_chains; ich++) {
            mmdb::Chain *chain_p = model_p->GetChain(ich);
            if (chain_p) {
               mmdb::Chain *copy_chain_p = new mmdb::Chain;
               copy_chain_p->Copy(chain_p);
               std::string ch_id = index_to_chain_id(chain_idx);
               copy_chain_p->SetChainID(ch_id.c_str());
               mol_first_model_p->AddChain(copy_chain_p);
               chain_idx++;
            }
         }
      }
   }
}

