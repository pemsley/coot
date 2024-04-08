/*
 * coot-utils/test-elastic.cc
 *
 * Copyright 2012 by University of York
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
//
#include <iostream>
#include <string>
#include "elastic.hh"

void
coot::test_elastic() {

   std::string file_name = "bit.pdb";
    file_name = "1x8b.pdb";
   mmdb::Manager *mol = new mmdb::Manager;
   mol->ReadCoorFile(file_name.c_str());

   int selhnd = mol->NewSelection();

   mol->SelectAtoms(selhnd, 0,
                    "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*",  // residue name
                    "*",  // Residue must contain this atom name?
                    "*",  // Residue must contain this Element?
                    "*"  // altLocs
                    );

   mmdb::PPAtom atom_selection = NULL;
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
      for (int i=0; i<n; i++) {
         // v = sqrt(v + 65.7);
         v = acos(v * 0.01);
      }
      std::cout << "   " << v << std::endl;
   }

   return 0;
}
