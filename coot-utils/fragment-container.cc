/*
 * coot-utils/fragment-container.cc
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
                     // std::cout << "transferring residues running of size " << residues_running.size() << std::endl;
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


// "fragment" in this case means a run of residues of the given size. There will be lots of fragments
// typically and they will be overlapping.
//
coot::fragment_container_t
coot::make_overlapping_fragments(mmdb::Manager *mol, const std::string &chain_id, unsigned int fragment_length) {

   coot::fragment_container_t fc;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (! model_p) return fc;

   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      std::string chain_id_local(chain_p->GetChainID());
      if (chain_id_local == chain_id) {
         int n_res = chain_p->GetNumberOfResidues();
         if (n_res > 0) {
            for (int i_res=0; i_res<n_res; i_res++) {
               mmdb::Residue *res = chain_p->GetResidue(i_res);
               if (! res) continue;
               mmdb::Residue *residue_start = res;
               int resno_start = res->GetSeqNum();
               int this_run_max = fragment_length;
               if ((i_res + this_run_max) >= n_res) this_run_max = n_res - i_res - 1;
               if (false)
                  std::cout << "debug:: combi: " << i_res + this_run_max << " this_run_max " << this_run_max
                            << " for n_res " << n_res << std::endl;
               std::vector<mmdb::Residue *> residues_running;
               residues_running.push_back(residue_start);
               for (int j_res=1; j_res<this_run_max; j_res++) {
                  mmdb::Residue *r = chain_p->GetResidue(i_res + j_res);
                  if (r) {
                     int resno_this = r->GetSeqNum();
                     if ((resno_this - resno_start) == j_res) {
                        residues_running.push_back(r);
                     }
                  }
               }
               if (false)
                  std::cout << "debug:: residues_running size " << residues_running.size() << std::endl;
               if (residues_running.size() > 10) { // or something
                  fragment_container_t::fragment_range_t fr(chain_id,
                                                            residue_spec_t(*residues_running.begin()),
                                                            residue_spec_t( residues_running.back()));
                  fc.add(fr);
               }
            }
         }
      }
   }

   return fc;
}
