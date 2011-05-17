
#include <goocanvas.h>

#include <iostream>

#include "lbg.hh"

#include "lig-build.hh"
#include "lbg-molfile.hh"

#include "graphics-c-interface-functions-blanks.cc"

int test_molfile() {

   lig_build::molfile_molecule_t m;
   m.read("piper-3.mol");

   return 0;
} 

int test_topological_equivalence(const std::string &file_name) {

   int r = 0;

   std::cout << "reading mol2 file: " << file_name << std::endl;
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

   int r = 0;

   std::string file_name = "phe.mol";
   if (argc > 1)
      file_name = argv[1];
      
   // r = test_topological_equivalence(file_name);

   r = test_molfile();

   if (r == 0)
      std::cout << "test failed" << std::endl;

   return 0;

}
