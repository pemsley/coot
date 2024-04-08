/*
 * geometry/mol-utils-2.cc
 *
 * Copyright 2017 by Medical Research Council
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

#include "mol-utils-2.hh"

// put this in mol-utils.cc when that arrives in this branch

std::map<std::string, std::pair<int, int> >
coot::get_residue_number_limits(mmdb::Manager *mol) {

   std::map<std::string, std::pair<int, int> > limits;
   if (mol) {
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    std::string chain_id = chain_p->GetChainID();
	    int nres = chain_p->GetNumberOfResidues();
	    int min_res_no =  9999999;
	    int max_res_no = -9999999;
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       int res_no = residue_p->GetSeqNum();
	       if (res_no < min_res_no) min_res_no = res_no;
	       if (res_no > max_res_no) max_res_no = res_no;
	    }
	    if (min_res_no < 9999999) {
	       if (max_res_no > -9999999) {
		  std::pair<int, int> res_lims(min_res_no, max_res_no);
		  limits[chain_id] = res_lims;
	       }
	    }
	 }
      }
   }
   return limits;
}


