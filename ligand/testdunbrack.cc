/* ligand/test-dunbrack.cc
 * 
 * Copyright 2005 by Paul Emsley, The University of York
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
 * 02110-1301, USA.
 */

#include <iostream>
#include "dunbrack.hh"
#include "coot-utils.hh"

int main(int argc, char **argv) {

   CResidue *r = 0; 
   coot::dunbrack d(r, "FRE");

//    std::string test_string = "45%";
//    std::vector<std::string> parts = coot::util::split_string(test_string, "%");
//    std::cout << test_string << " has " << parts.size() << " parts" << std::endl;
//    test_string = "45";
//    parts = coot::util::split_string(test_string, "%");
//    std::cout << test_string << " has " << parts.size() << " parts" << std::endl;

   if (argc > 1) {
      std::string n(argv[1]);
      std::cout << "Reading Library file: " << n << std::endl;
      d.read_penultimate_library(n);
   } else {
      std::cout << "Usage: " << argv[0] << " dictionary-file-name " << std::endl;
   }


   return 0;
}
