/*
 * src/fix-nomenclature.cc
 *
 * Copyright 2007 by University of York
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

#include <stdlib.h>
#include <iostream> // fixes undefined strchr, strchrr problems

#include "geometry/protein-geometry.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"

#include "coot-nomenclature.hh"

int
main(int argc, char **argv) {

//    int a = 5;
//    std::cout << a << "\n";
//    a = a | 4; 
//    std::cout << a << "\n";
//    a = a | 8; 
//    std::cout << a << "\n";
//    a = a | 16; 
//    std::cout << a << "\n";
//    a = a | 16; 
//    std::cout << a << "\n";

//    a = 1;
//    a = a << 1;
//    std::cout << a << "\n";
//    a = a << 1;
//    std::cout << a << "\n";
//    a = a <<= 1;
//    std::cout << a << "\n";


   if (argc < 3) {
      std::cout << "Usage: " << argv[0] << " in-filename out-filename\n";
      exit(1);
   } else {
      coot::protein_geometry geom;
      std::string filename = argv[1];
      atom_selection_container_t asc = get_atom_selection(filename, true, false, false);
      coot::nomenclature n(asc.mol);
      asc.mol->WritePDBASCII(argv[2]);
      n.fix(&geom);
   }
   return 0;
}
