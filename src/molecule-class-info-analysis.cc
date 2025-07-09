/*
 * src/molecule-class-info-analysis.cc
 *
 * Copyright 2020 by Medical Research Council
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


#include "molecule-class-info.h"


coot::model_composition_stats_t
molecule_class_info_t::get_model_composition_statistics() const {

   coot::model_composition_stats_t m;

   return m;

}

std::vector<std::string>
molecule_class_info_t::get_types_in_molecule() const {

   // 20250513-PE move this into coot-utils

   std::vector<std::string> v;
   std::set<std::string> s;
   mmdb::Manager *mol = atom_sel.mol;
   if (mol) {
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	 mmdb::Model *model_p = mol->GetModel(imod);
	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       int n_res = chain_p->GetNumberOfResidues();
	       for (int ires=0; ires<n_res; ires++) {
		  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  if (residue_p) {
		     std::string type = residue_p->GetResName();
		     s.insert(type);
		  }
	       }
	    }
	 }
      }
   }
   for (const auto &item : s) {
      v.push_back(item);
   }
   return v;
}
