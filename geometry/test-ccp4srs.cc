/* geometry/test-geometry.cc
 * 
 * Copyright 2004  The University of York
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


#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>

#include <iostream>

#include "utils/coot-utils.hh"
#include "protein-geometry.hh"

void
list_monomers() { 
         
#ifdef HAVE_CCP4SRS   
   int RC = ccp4srs::CCP4SRS_FileNotFound; // initial status.
   std::string dir;
   const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
   const char *d2 = getenv("CCP4");
  
   if (d1) {
      if (coot::file_exists(d1))
         dir = d1;
   } else {
      if (d2) {
         std::string dir_a = coot::util::append_dir_dir(d2, "share");
         std::string dir_b = coot::util::append_dir_dir(dir_a, "ccp4srs");
         if (coot::file_exists(dir_b))
            dir = dir_b;
      }
   }
         
   if (dir.length()) {
      std::cout << "INFO:: CCP4SRS::loadIndex from dir: " << dir << std::endl;
      ccp4srs::Manager *ccp4srs = new ccp4srs::Manager;
      RC = ccp4srs->loadIndex(dir.c_str());
      if (true)
         std::cout << "init_ccprsrs() ... loadIndex() returned " << RC << std::endl;
      if (RC != ccp4srs::CCP4SRS_Ok) {
         std::cout << "CCP4SRS init problem." << std::endl;
         delete ccp4srs;
         ccp4srs = NULL;
      } else {
         int l = ccp4srs->n_entries();
         for (int i=0; i<l ;i++)  {
         ccp4srs::Monomer *Monomer = ccp4srs->getMonomer(i, NULL);
         if (Monomer) { 
             std::string group = coot::util::downcase(Monomer->chem_type());
             if (group == "non-polymer") { 
                std::cout << Monomer->ID() << "\n";
             } else { 
                // std::cout << i << " " << group << std::endl;
             }
           }
         }
      }
   } else {
     std::cout << "WARNING:: init_ccp4srs() no dir" << std::endl;
   } 
#endif
}

void
make_dict(const std::string &comp_id, const std::string &file_name) {

#ifdef HAVE_CCP4SRS   
   coot::protein_geometry pg;
   int status = pg.init_ccp4srs("nothing");

   pg.fill_using_ccp4srs(comp_id);

   std::pair<bool, coot::dictionary_residue_restraints_t> r = pg.get_monomer_restraints(comp_id);
   if (r.first) { 
      r.second.write_cif(file_name);
   }

#endif

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
 
   if (argc > 1) {
      std::string a1(argv[1]);
      if (a1 == "list") {  // non-polymers only
         list_monomers();
      } else { 
        if (argc == 3) { 
           if (a1 == "make-dict") { 
              std::string comp_id = argv[2];
              std::string out_file_name(comp_id + "-SRS.cif");
              make_dict(comp_id, out_file_name);
           } 
        }
     }
   } else { 

      std::string comp_id = "LYS";
      pg.try_dynamic_add(comp_id, read_number++);
      std::pair<bool, coot::dictionary_residue_restraints_t> r = pg.get_monomer_restraints(comp_id);

      if (r.first) {
         coot::dictionary_residue_restraints_t &restraints = r.second;
         double local_search_similarity = 0.9;
         int n_atoms = restraints.residue_info.number_atoms_nh;
         mmdb::math::Graph *graph = restraints.make_graph(true);

	 // this crashes in mmdb graph matching code.
         if (false) {
	    std::vector<coot::match_results_t> v =
	       pg.compare_vs_ccp4srs(graph, local_search_similarity, n_atoms);
	    std::cout << "INFO:: test-ccp4srs: got " << v.size() << " SRS matches " << std::endl;
         }

	 if (true) {
	    pg.fill_using_ccp4srs("TRP");
	    std::string match_str = "trypto";
	    std::vector<std::pair<std::string, std::string> > v = pg.matching_ccp4srs_residues_names(match_str);
	    std::cout << " Found " << v.size() << " monomers with names matching " << match_str << std::endl;
	    for (unsigned int i=0; i<v.size(); i++)
	       std::cout << v[i].first << " " << v[i].second << std::endl;
	 }
	 
         delete graph;
      } else {
         std::cout << "No comp_id " << comp_id << " in dictionary." << std::endl;
      }
   }

#else
   std::cout << "This build does not include CCP4 SRS" << std::endl;
#endif // HAVE_CCP4SRS   
   return 0; 
}

