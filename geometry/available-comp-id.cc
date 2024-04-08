/* geometry/available-comp-id.cc
 * 
 * Copyright 2015 Paul Emsley
 * Copyright 2015 by Medical Research Council
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include <iostream>
#include "utils/coot-utils.hh"
#include "protein-geometry.hh"

int
main(int argc, char **argv) {

#ifndef HAVE_CCP4SRS   
   std::cout << "This build does not incorporate the necessary CCP4 SRS library."
	     << std::endl;
#else
   
   std::string code_start;

   if (argc > 1) {
      std::string t = argv[1];
      if (t.length() <= 3)
	 code_start = t;
   }
   
   int n_top = 10;


   if (argc > 2) {
      std::string t = argv[2];
      try {
	 n_top = coot::util::string_to_int(t);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   } 
   
   coot::protein_geometry pg;
   int status = pg.init_ccp4srs("");
   std::vector<std::string> new_codes = pg.get_available_ligand_comp_id(code_start, n_top);

   if (new_codes.size() == 0) {
      if (code_start.length() < 3) 
	 std::cout << "INFO:: No matching codes starting with \"" << code_start
		   << "\"" << std::endl;
      if (code_start.length() >= 3) 
	 std::cout << "INFO:: comp-id " << code_start << " is not available" << std::endl;

   } else {

      if (code_start.length() == 3) {
	 std::cout << "INFO:: comp-id " << code_start << " is available" << std::endl;
      } else { 
	 for (unsigned int i=0; i<new_codes.size(); i++) {
	    if (i>0)
	       if (i%10 == 0) 
		  std::cout << std::endl;
	    std::cout << "    " << new_codes[i];
	 }
	 std::cout << std::endl;
      }
   }

#endif   
   
   return 0; 
}

