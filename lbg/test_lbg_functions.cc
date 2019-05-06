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

#include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"

// #include "graphics-c-interface-functions-blanks.cc"

int test_molfile() {

   lig_build::molfile_molecule_t m;
   m.read("piper-3.mol");

   return 0;
}

int test_split_molecule() {

   int status = 0;

   // ifdef this out because:
   // 
   // Undefined symbols for architecture x86_64:
   //   "lig_build::molecule_t<widgeted_atom_t, widgeted_bond_t>::~molecule_t()", referenced from:
   //       widgeted_molecule_t::widgeted_molecule_t() in test_lbg_functions.o
   // ld: symbol(s) not found for architecture x86_64
   // clang: error: linker command failed with exit code 1 (use -v to see invocation)

#ifdef XYZ_ABC

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
   mmdb::Manager *pdb_mol = NULL;
   widgeted_molecule_t wm(m, pdb_mol);

   topological_equivalence_t top(wm.atoms, wm.bonds);

   return r;
}

int test_ccp4srs_graph_search() {

   int r = 0;

#ifdef HAVE_CCP4SRS

   double search_similarity = 0.9;
   bool inc_Hs = false;
   std::string monomer_type = "TRP";
   unsigned int n_atoms = 12;

   coot::protein_geometry *geom_p = new coot::protein_geometry;
   const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
   std::string srs_dir = PKGDATADIR;
   if (d1)
      srs_dir = d1;
   std::cout << "---- geom_p init_ccp4srs with srs_dir " << srs_dir
	     << std::endl;
   geom_p->init_ccp4srs(srs_dir);
   geom_p->try_dynamic_add(monomer_type, true);
   std::pair<bool, coot::dictionary_residue_restraints_t> rest =
      geom_p->get_monomer_restraints(monomer_type);
   mmdb::math::Graph *graph = rest.second.make_graph(inc_Hs);
   graph->SetName ("Coot-LBG-Query");
   graph->MakeVertexIDs();

   int build_result = graph->Build(false);

   if (build_result != 0) {

      std::cout << "Bad graph build result" << std::endl;

   } else {

      graph->MakeSymmetryRelief(false);
      graph->Print();
      std::cout << "graph search using similarity  " << search_similarity << std::endl;
      std::cout << "graph build returns: " << build_result << std::endl;
      std::vector<coot::match_results_t> v =
	 geom_p->compare_vs_ccp4srs(graph, search_similarity, n_atoms);
      delete graph;
      std::cout << "found " << v.size() << " close matches" << std::endl;
   }
#endif // HAVE_CCP4SRS
   return r;
}


// #include "cairo-molecule.hh"

// this needs to be here (also).  It is in wmolecule.cc and hence the library also.
// But if this is not here I get unresovled symbol for this destructor when compiling
// this exectuable on the mac (clang).
template<class cairo_atom_t, class cairo_bond_t> lig_build::molecule_t<cairo_atom_t, cairo_bond_t>::~molecule_t() {}

int test_mol_to_cairo(int argc, char **argv) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#endif // MAKE_ENHANCED_LIGAND_TOOLS
   return 1;
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

   // r = test_split_molecule();

   // r = test_ccp4srs_graph_search();

   r = test_mol_to_cairo(argc, argv);

   if (r == 0)
      std::cout << "test failed" << std::endl;

   return 0;

}
