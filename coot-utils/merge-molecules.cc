/*
 * coot-utils/merge-molecules.cc
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

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

