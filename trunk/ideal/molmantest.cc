
#include <iostream>
#include <string>
#include <vector>

#include "mmdb_manager.h"
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
