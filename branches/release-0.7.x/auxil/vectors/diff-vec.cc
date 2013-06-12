
#include <iostream>

#include <mmdb/mmdb_manager.h>

int main(int argc, char **argv) {

   if (argc > 1) {
      char *filename = argv[1];
      CMMDBManager *mol = new CMMDBManager;
      mol->ReadCoorFile(filename);

      std::cout << "(define obj (new-generic-object-number \"EJD's vectors\"))"
		<< std::endl;

      if (mol) {
	 int imod = 1;
	 int ncount = 11; 
	 
	 CModel *model_p = mol->GetModel(imod);
	 CChain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    CAtom *at;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);

		  // inner atom loop
		  for (int jchain=0; jchain<nchains; jchain++) {
		     CChain *chain2_p = model_p->GetChain(jchain);
		     int nres2 = chain2_p->GetNumberOfResidues();
		     PCResidue residue2_p;
		     CAtom *at2;
		     for (int jres=0; jres<nres2; jres++) { 
			residue2_p = chain2_p->GetResidue(jres);
			int n_atoms2 = residue2_p->GetNumberOfAtoms();
			
			for (int jat=0; jat<n_atoms2; jat++) {
			   at2 = residue2_p->GetAtom(jat);

			   if (at != at2) {
			      ncount++;

			      if (ncount > 10) {
				 ncount = 0;
				 float diff_x = at->x - at2->x;
				 float diff_y = at->y - at2->y;
				 float diff_z = at->z - at2->z;
				 
				 std::string colour = "grey";
				 int width = 2;
				 float m = at->occupancy * at2->occupancy;
				 if (std::string(at->element) != " C")
				    colour = "green";
				 if (std::string(at2->element) != " C")
				    colour = "green";
				 if (std::string(at->element) == " S") {
				    colour = "yellow";
				    width = 4;
				 }
				 if (std::string(at2->element) == " S") {
				    colour = "yellow";
				    width = 4;
				 }
				 if (std::string(at->element) == " N")
				    colour = "blue";
				 if (std::string(at2->element) == " N")
				    colour = "blue";
				 if (std::string(at->element) == " O")
				    colour = "red";
				 if (std::string(at2->element) == " O")
				    colour = "red";
				 
// 				 std::cout << "(to-generic-object-add-line obj "
// 					   << "\"" << colour << "\" "
// 					   << width << " "
// 					   << 0.0 << " " << 0.0 << " " << 0.0 << " "
// 					   << diff_x << " " << diff_y << " " << diff_z
// 					   << ")" << std::endl;
				 std::cout << "(to-generic-object-add-point obj "
					   << "\"" << colour << "\" "
					   << width << " "
					   << diff_x << " " << diff_y << " " << diff_z
					   << ")" << std::endl;
			      }
			   } 
			}
		     }
		  }
	       }
	    }
	 }
      }
      std::cout << "(set-display-generic-object obj 1)" << std::endl;
   } 
} 

