/* mini-mol/testminimol.cc
 * 
 * Copyright  2004 The University of York
 * Author: Paul Emsley
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

#include "mini-mol.hh"

int
main(int argc, char **argv) {

//     std::vector<int> vec(5);
//     vec[3] = 12;
//     vec.resize(10);
//     std::vector<int> vec2 = vec;
    
//     for(int i=0; i<vec.size(); i++)
//        std::cout << i << " " << vec2[i] << std::endl;

   if (argc > 1) {
      std::string pdb_name(argv[1]);
      CMMDBManager mol;
      coot::minimol::molecule myminimol;
      myminimol.read_file(pdb_name);
      
      std::vector<coot::minimol::atom *> atoms =
	 myminimol.select_atoms_serial();

      atoms[2]->pos = clipper::Coord_orth(1000, 1000, 1000);

      myminimol.write_file("minimoled-mol.pdb", 20.0);

      coot::minimol::molecule m2 = myminimol.molecule_of_atom_types(std::string(" CA "));

      m2.write_file("minimoled-cas.pdb", 20.0);
      for(int ifrag=0; ifrag<m2.fragments.size(); ifrag++)
	 std::cout << "m2 fragment " << ifrag << " has "
		   << m2[ifrag].residues.size()
		   << " residues " << std::endl; 
      std::string unused_chain = m2.unused_chain_id("W");
      std::cout << "Found unused chain " << unused_chain << std::endl;

   }
   
   coot::minimol::molecule a;
   int ifrag = a.fragment_for_chain("A");
   coot::minimol::residue res(1,"ALA");
   res.addatom(" CA ", " C", clipper::Coord_orth(0,0,0), "", 10.0, 1.0);
   a[ifrag].addresidue(res, 0);

   a.write_file("atest.pdb", 20.0);

   // double v = 0.12;
   // good - returns negative when given negative 0.119429  -0.119429
   // std::cout << atan2(v,1) << "  " << atan2(-v,1) << std::endl;

   // Does mmdb unmangle hydrogen names?  No, it turns out.
   // 
//    if (0) { 
//       CMMDBManager *mol = new CMMDBManager;
//       mol->ReadCoorFile("prodrg.pdb");
//       int imod = 1;
//       CModel *model_p = mol->GetModel(imod);
//       CChain *chain_p;
//       // run over chains of the existing mol
//       int nchains = model_p->GetNumberOfChains();
//       for (int ichain=0; ichain<nchains; ichain++) {
// 	 chain_p = model_p->GetChain(ichain);
// 	 int nres = chain_p->GetNumberOfResidues();
// 	 CResidue *residue_p;
// 	 CAtom *at;
// 	 for (int ires=0; ires<nres; ires++) { 
// 	    residue_p = chain_p->GetResidue(ires);
// 	    int n_atoms = residue_p->GetNumberOfAtoms();
	 
// 	    for (int iat=0; iat<n_atoms; iat++) {
// 	       at = residue_p->GetAtom(iat);
// 	       std::cout << chain_p->GetChainID() << " " << residue_p->GetSeqNum() << " "
// 			 << at->GetAtomName() << std::endl;
// 	    }
// 	 }
//       }
//    }
   
   
   return 0;
}
