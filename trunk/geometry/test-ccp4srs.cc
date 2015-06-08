/* geometry/test-geometry.cc
 * 
 * Copyright 2004  The University of York
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


#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>

#include "protein-geometry.hh"

std::vector<std::string>
get_new_code(const coot::protein_geometry &pg, const std::string &code_start, unsigned int n_top) {

   std::vector<std::string> r = pg.get_available_ligand_comp_id(code_start, n_top);
   return r;
} 

int
main(int argc, char **argv) {

#ifdef HAVE_CCP4SRS   
   std::string filename;
   int read_number = 1;

   // if srs-dir is not given, then the default should be $CCP4_MASTER/share/ccp4srs

   coot::protein_geometry pg;
   int status = pg.init_ccp4srs("nothing");
   std::cout << "INFO:: test-ccp4srs: status: " << status << std::endl;
   std::string comp_id = "LYS";
   pg.try_dynamic_add(comp_id, read_number++);
   std::pair<bool, coot::dictionary_residue_restraints_t> r = pg.get_monomer_restraints(comp_id);

   if (r.first) {
      coot::dictionary_residue_restraints_t &restraints = r.second;
      double local_search_similarity = 0.9;
      int n_atoms = restraints.residue_info.number_atoms_nh;
      mmdb::math::Graph *graph = restraints.make_graph(true);
      if (false) { 
	 std::vector<coot::match_results_t> v =
	    pg.compare_vs_ccp4srs(graph, local_search_similarity, n_atoms);
	 std::cout << "INFO:: test-ccp4srs: got " << v.size() << " SRS matches " << std::endl;
      }
      delete graph;
   } else {
      std::cout << "No comp_id " << comp_id << " in dictionary." << std::endl;
   }

   std::string code_start;

   if (argc > 1) {
      std::string t = argv[1];
      if (t.length() <= 3) 
	 code_start = t;
   }
   int n_top = 10;
   
   std::vector<std::string> new_codes = get_new_code(pg, code_start, n_top);
   for (unsigned int i=0; i<new_codes.size(); i++) {
      if (i>0)
	 if (i%10 == 0) 
	    std::cout << std::endl;
      std::cout << "    " << new_codes[i];
   }
   std::cout << std::endl;
   
#endif // HAVE_CCP4SRS   
   return 0; 
}

