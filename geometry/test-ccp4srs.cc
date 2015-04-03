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

int
main(int argc, char **argv) {

#ifdef HAVE_CCP4SRS   
   std::string filename;
   int read_number = 1;

   coot::protein_geometry pg;
   int status = pg.init_ccp4srs("nothing");
   std::cout << "status: " << status << std::endl;
   std::string comp_id = "LYS";
   pg.try_dynamic_add(comp_id, read_number++);
   std::pair<bool, coot::dictionary_residue_restraints_t> r = pg.get_monomer_restraints(comp_id);

   if (r.first) {
      coot::dictionary_residue_restraints_t &restraints = r.second;
      double local_search_similarity = 0.9;
      int n_atoms = restraints.residue_info.number_atoms_nh;
      mmdb::math::Graph *graph = restraints.make_graph(true);
      std::vector<coot::match_results_t> v =
	 pg.compare_vs_ccp4srs(graph, local_search_similarity, n_atoms);
      std::cout << "got " << v.size() << " SRS matches " << std::endl;
   }
#endif // HAVE_CCP4SRS   
   return 0; 
}

