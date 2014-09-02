
#include <iostream>
#include <string>
#include "elastic.hh"

void
coot::test_elastic() {

   std::string file_name = "bit.pdb";
    file_name = "1x8b.pdb";
   CMMDBManager *mol = new CMMDBManager;
   mol->ReadCoorFile(file_name.c_str());

   int selhnd = mol->NewSelection();

   mol->SelectAtoms(selhnd, 0,
		    "*",
		    ANY_RES, "*",
		    ANY_RES, "*",
		    "*",  // residue name
		    "*",  // Residue must contain this atom name?
		    "*",  // Residue must contain this Element?
		    "*"  // altLocs
		    );

   PPCAtom atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

   elastic_network_model_t enm(mol, selhnd, 1, 20, 2000000);
   mol->DeleteSelection(selhnd);
   
} 


int main(int argc, char **argv) {

   coot::test_elastic();

   if (0) {
      int n = 100000000; // ~100 million sqrt()s/second
      double v = 0.9; 
      for (unsigned int i=0; i<n; i++) { 
	 // v = sqrt(v + 65.7);
	 v = acos(v * 0.01);
      }
      std::cout << "   " << v << std::endl;
   }

   return 0;
} 
