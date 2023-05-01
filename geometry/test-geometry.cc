/* geometry/test-geometry.cc
 * 
 * Copyright 2004  The University of York
 * Copyright 2015 by Medical Research Council
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
#include "protein-donor-acceptors.hh"

std::vector <std::string> protein_monomers();

std::vector <std::string>
protein_monomers() {

   std::vector <std::string> s;

   s.push_back("/lib/data/monomers/a/ALA.cif");
   s.push_back("/lib/data/monomers/c/CYS.cif");
   s.push_back("/lib/data/monomers/a/ASP.cif");
   s.push_back("/lib/data/monomers/a/ASN.cif");
   s.push_back("/lib/data/monomers/g/GLN.cif");
   s.push_back("/lib/data/monomers/g/GLY.cif");
   s.push_back("/lib/data/monomers/g/GLU.cif");
   s.push_back("/lib/data/monomers/p/PHE.cif");
   s.push_back("/lib/data/monomers/h/HIS.cif");
   s.push_back("/lib/data/monomers/i/ILE.cif");
   s.push_back("/lib/data/monomers/l/LYS.cif");
   s.push_back("/lib/data/monomers/l/LEU.cif");
   s.push_back("/lib/data/monomers/m/MET.cif");
   s.push_back("/lib/data/monomers/p/PRO.cif");
   s.push_back("/lib/data/monomers/a/ARG.cif");
   s.push_back("/lib/data/monomers/s/SER.cif");
   s.push_back("/lib/data/monomers/t/THR.cif");
   s.push_back("/lib/data/monomers/v/VAL.cif");
   s.push_back("/lib/data/monomers/t/TRP.cif");
   s.push_back("/lib/data/monomers/t/TYR.cif");   

   return s;
}

int
ccp4_setup(int argc, char **argv) {
 
   std::string filename;
   int read_number = 1;

   if (argc <= 1) {
      
      char *s = getenv("CCP4");
      if (s) {
	 filename = s;
	 // contains the linkages:
	 filename += "/lib/data/monomers/list/mon_lib_list.cif";
	 // now check that that file is there:
	 struct stat buf;
	 stat(filename.c_str(), &buf);
	 if (S_ISREG(buf.st_mode)) {
	    std::cout << "reading mon_lib_list..." << std::endl; 
	    coot::protein_geometry prot;
	    prot.init_refmac_mon_lib(filename, -1);

	    // now the protein monomers:
	    std::vector <std::string> protein_mono = protein_monomers();
	    for (unsigned int i=0; i<protein_mono.size(); i++) {
	       
	       filename = s; 
	       filename += protein_mono[i];

	       // filename = "ALA.cif";
	       stat(filename.c_str(), &buf);
	       if (S_ISREG(buf.st_mode)) {
		  prot.init_refmac_mon_lib(filename, read_number++);
	       } else {
		  std::cout << "ERROR: dictionary " << filename
			    << " is not a regular file" << std::endl;
	       }
	    }
	    prot.info();
	 } else {
	    std::cout << filename << " is not a regular file (for links)"
		      << std::endl;
	 }
      } else { 
	 std::cout << argv[0] << " <path-to-mon_lib_prot.cif>" << std::endl;
	 std::cout << ".... or setup ccp4. " << std::endl; 
      }
      
   } else { 

      for (int i=1; i<argc; i++) {
	 filename = argv[i];
	 coot::protein_geometry prot;
	 prot.init_refmac_mon_lib(filename, read_number++);
      }
   }
   return 0; 
}

int test_improper_dihdedrals(int argc, char **argv) {

   coot::protein_geometry prot;
   prot.init_standard();
   int n = prot.size();

   for (int i=0; i<n; i++) {
      const std::pair<int, coot::dictionary_residue_restraints_t> rest = prot[i];
      int n_pl_rest = rest.second.plane_restraint.size();
      std::cout << "# type " << rest.second.residue_info.comp_id << " " << n_pl_rest << " plane restraints"
                << std::endl;
      for (int j=0; j<n_pl_rest; j++) {
	 std::vector<coot::atom_name_quad> qv = rest.second.plane_restraint_to_improper_dihedrals(j);
	 for (std::size_t k=0; k<qv.size(); k++) {
	    std::cout << "    Here they are: " << k << " " << qv[k] << std::endl;
	 }
	 std::cout << "   " << rest.second.residue_info.comp_id << " " << j  << " got "
		   << qv.size() << " atom name quads " << std::endl;
      }
   }
   return 0;
}

void
test_quick_hbs() {

   coot::quick_protein_donor_acceptors qpda;
   qpda.test();

}


int main(int argc, char **argv) {

   // test_quick_hbs();

   test_improper_dihdedrals(argc, argv);


   return 0;
}
