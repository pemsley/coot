
#include <goocanvas.h>

#include <iostream>

#include "lbg.hh"

#include "lig-build.hh"
#include "lbg-molfile.hh"

#include "graphics-c-interface-functions-blanks.cc"

int test_topological_equivalence() {

   int r = 0;

   lig_build::molfile_molecule_t m;
   // m.read("test.mdl"); // caffeine
   // m.read("coot-non-chiral.mol"); // has equivalent atoms
   // m.read("propane.mol"); // has equivalent atoms
   m.read("phe.mol");
   CMMDBManager *pdb_mol = NULL;
   widgeted_molecule_t wm(m, pdb_mol);

   topological_eqivalence_t top(wm.atoms, wm.bonds);

   return r;
}

int main(int argc, char **argv) {

   gtk_init(&argc, &argv);

   int r = 0;

   r = test_topological_equivalence();

   if (r == 0)
      std::cout << "test failed" << std::endl;

   return 0;

}
