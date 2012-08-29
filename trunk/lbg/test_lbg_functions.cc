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
