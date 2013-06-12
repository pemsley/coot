/* lbg/test_lbg_functions.cc
 * 
 * Copyright 2012 by The University of Oxford
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

#include <goocanvas.h>

#include <iostream>

#include "lbg.hh"

#include "lig-build.hh"
#include "lbg-molfile.hh"

// #include "graphics-c-interface-functions-blanks.cc"

int test_molfile() {

   lig_build::molfile_molecule_t m;
   m.read("piper-3.mol");

   return 0;
}

int test_split_molecule() {
   int status = 0;
#ifdef MAKE_ENTERPRISE_TOOLS

   lig_build::molfile_molecule_t mol;
   // mol.read("FPX.mdl");
   mol.read("../src/recap-test.mol");

   std::cout << "------------- just after reading recap-test: " << std::endl;
   std::cout << "   scale correction " << mol.get_scale_correction().first << " "
	     << mol.get_scale_correction().second << std::endl;
   
   widgeted_molecule_t wmol(mol, NULL);
   lbg_info_t l;
   l.no_graphics_mode();
   RDKit::RWMol rdkm = l.rdkit_mol(wmol);
   int bond_index = 11;
   int atom_index = 6;
   RDKit::ROMol *main_sans_R = coot::split_molecule(rdkm, bond_index, atom_index);
   if (main_sans_R) {
      RDKit::MolToMolFile(*main_sans_R, "main-sans-R.mol");
      RDKit::MolToMolFile(rdkm, "starting.mol");


      // std::string smiles_string = "[*]n1cncn1"; // [*]c1ccnn1[*]
      // [*]c1ccccn1 [*]n1cnnc1 [*]n1cncn1 [*]c1ccncc1 [*]c1nc(C)cs1
      // [*]c1ncccn1 [*]n1cccn1

      std::vector<std::string> fragment_smiles;
      fragment_smiles.push_back("[*]n1cncn1");
      fragment_smiles.push_back("[*]c1ccnn1[*]");
      fragment_smiles.push_back("[*]c1ccccn1");
      fragment_smiles.push_back("[*]n1cnnc1");
      fragment_smiles.push_back("[*]n1cncn1");
      fragment_smiles.push_back("[*]c1ccncc1");
      fragment_smiles.push_back("[*]c1nc(C)cs1");
      fragment_smiles.push_back("[*]c1ncccn1");
      fragment_smiles.push_back("[*]n1cccn1");
      
      for (unsigned int is=0; is<fragment_smiles.size(); is++) { 
	 std::string smiles_string = fragment_smiles[is];
	 RDKit::ROMol *test_fragment = RDKit::SmilesToMol(smiles_string);
	 std::vector<RDKit::ROMol *> joined_mols = 
	    coot::join_molecules(*main_sans_R, atom_index, *test_fragment);
	 for (unsigned int imol=0; imol<joined_mols.size(); imol++) {
	    std::string file_name = "joined-";
	    file_name += coot::util::int_to_string(is);
	    file_name += ".mol";
	    RDKit::MolToMolFile(*joined_mols[imol], file_name);
	 }
      }
      
   } else {
      std::cout << "main_sans_R was null :-(" << std::endl;
   }
      
#endif
   return status;
}

int test_topological_equivalence(const std::string &file_name) {

   int r = 0;

   std::cout << "reading mol file: " << file_name << std::endl;
   lig_build::molfile_molecule_t m;
   // m.read("test.mdl"); // caffeine
   // m.read("coot-non-chiral.mol"); // has equivalent atoms
   // m.read("propane.mol"); // has equivalent atoms
   m.read(file_name);
   CMMDBManager *pdb_mol = NULL;
   widgeted_molecule_t wm(m, pdb_mol);

   topological_equivalence_t top(wm.atoms, wm.bonds);

   return r;
}

int main(int argc, char **argv) {

   gtk_init(&argc, &argv);
   Py_Initialize();
   PySys_SetArgv(argc, argv);

   int r = 0;

   std::string file_name = "phe.mol";
   if (argc > 1)
      file_name = argv[1];
      
   // r = test_topological_equivalence(file_name);

   // r = test_molfile();

   r = test_split_molecule();

   if (r == 0)
      std::cout << "test failed" << std::endl;

   return 0;

}
