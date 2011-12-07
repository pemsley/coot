/* ideal/molmantest.cc
 * 
 * Copyright 2004 The University of York
 * Author Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <iostream>
#include <string>
#include <vector>

#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "mmdb.h"


void molman(MyCMMDBManager *mol) {

   int nmodels = mol->GetNumberOfModels();
   std::cout << "models: " << nmodels << std::endl;
   for (int imodel=1; imodel<=nmodels; imodel++) {
      int nchains = mol->GetNumberOfChains(imodel);
      std::cout << "model " << imodel << " has "
		<< nchains << " chains"  << std::endl;
      for (int ichain=0; ichain<nchains; ichain++) {
	 int nres = mol->GetNumberOfResidues(imodel,ichain);
	 PCChain chn = mol->GetChain(imodel, ichain);
	 std::string chain_name = chn->GetChainID();
	 std::cout << "   chain: " << chain_name << " " << ichain
		   << " has " << nres << " residues" << std::endl;
	 for (int ires=0; ires<nres; ires++) {
	    PCResidue res = mol->GetResidue(imodel,ichain,ires);
	    int natoms = res->GetNumberOfAtoms();
	    int seqno  = res->GetSeqNum();
	    std::cout << "      residue " << ires << ", seqno "
		      << seqno << " has " << natoms << " atoms "
		      << std::endl;
	 }
      }
   }
}


main(int argc, char **argv) {
   
   if (argc < 2) {
      std::cout << "usage: " << argv[0] << " pdb_filename " << std::endl;
   } else {
      string pdb_file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name); 
      molman(asc.mol);
   }

   return 0; 

} 
