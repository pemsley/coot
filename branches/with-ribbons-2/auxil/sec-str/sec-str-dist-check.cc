
#include <unistd.h>
#include <iostream>
#include "sec-str-dist-check.hh"
#include "clipper/core/coords.h"

int main(int argc, char **argv) {

   if (argc > 1) { 
      std::string filename = argv[1];
      CMMDBManager *mol = get_mol(filename);
      if (mol) { 
	 CModel *model_p = mol->GetModel(1);
	 int status = model_p->CalcSecStructure(1);
	 if (status == SSERC_Ok) {
	    std::cout << "INFO:: SSE status was OK\n";
            CSSContainer helices;
            
	    distance_checks(model_p);
	 } else {
	    std::cout << "INFO:: SSE status was bad\n" << status << "\n";
	 }
      }
   } else {
      std::cout << "Usage: " << argv[0] << " <pdb-file-name>"
		<< std::endl;
   } 
   return 0;
}

// Running over all of Top500:
// oxygen distance
//   Strand stats: 4.5171  +/- 0.356533
//   Helix  stats: 3.56436 +/- 0.333756
// 
void
distance_checks(CModel *model_p) {

   int imod = 1;
      
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p_1;
      PCResidue residue_p_2;
      for (int ires=0; ires<(nres-1); ires++) {
	 residue_p_1 = chain_p->GetResidue(ires);
	 residue_p_2 = chain_p->GetResidue(ires+1);
	 int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
	 int n_atoms_2 = residue_p_2->GetNumberOfAtoms();
	 if (n_atoms_1 && n_atoms_2) {
	    if (residue_p_1->GetSeqNum() == (residue_p_2->GetSeqNum()-1)) {
	       int sse_1 = SSE(residue_p_1);
	       int sse_2 = SSE(residue_p_2);
	       if (((sse_1 == SSE_Strand) && (sse_2 == SSE_Strand)) ||
		   ((sse_1 == SSE_Helix)  && (sse_2 == SSE_Helix))) {
		  float d = oxygen_check(residue_p_1, residue_p_2);
		  if (d > 0) {
		     std::string ss("Helix ");
		     if (sse_1 == SSE_Strand)
			ss = "Strand";
		     std::cout << ss << " " << d << "  " 
			       << residue_p_1->GetChainID()
			       << " "
			       << residue_p_1->GetSeqNum()
			       << std::endl;
		  }
	       }
	    }
	 }
      }
   }
}

int SSE(CResidue *res) {
   return res->SSE;
} 

CMMDBManager *get_mol(const std::string &filename) {

   CMMDBManager *MMDBManager = new CMMDBManager();
   int err = MMDBManager->ReadCoorFile((char *)filename.c_str());
   if (err) {
      std::cout << "Error reading " << filename << std::endl;
      delete MMDBManager;
      MMDBManager = 0;
   }
   return MMDBManager;
} 



// return -1 on failure
// 
float oxygen_check(CResidue *residue_p_1, CResidue *residue_p_2) {

   int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
   int n_atoms_2 = residue_p_2->GetNumberOfAtoms();
   CAtom *at, *at2;
   float d = -1;
   short int done_it = 0; 
   
   for (int iat=0; iat<n_atoms_1; iat++) {
      at = residue_p_1->GetAtom(iat);
      std::string atom_name (at->name);
      if (atom_name == " O  ") {
	 for (int iat2=0; iat2<n_atoms_2; iat2++) {
	    at2 = residue_p_2->GetAtom(iat2);
	    std::string atom_name_2 (at2->name);
	    if (atom_name_2 == " O  ") {
	       std::string altconf1 = at->altLoc;
	       std::string altconf2 = at2->altLoc;
	       if ((altconf1 == "") && (altconf2 == "")) { 
		  clipper::Coord_orth c1(at->x,  at->y,  at->z);
		  clipper::Coord_orth c2(at2->x, at2->y, at2->z);
		  d = clipper::Coord_orth::length(c1, c2);
		  done_it = 1;
		  break;
	       }
	    }
	 }
      }
      if (done_it)
	 break;
   }
   return d;
}
