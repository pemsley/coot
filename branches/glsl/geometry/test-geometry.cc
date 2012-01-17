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

std::vector <std::string> protein_monomers();

std::vector <std::string>
protein_monomers() {

   std::vector <std::string> s;

   s.push_back("/data/monomers/a/ALA.cif");
   s.push_back("/data/monomers/c/CYS.cif");
   s.push_back("/data/monomers/a/ASP.cif");
   s.push_back("/data/monomers/a/ASN.cif");
   s.push_back("/data/monomers/g/GLN.cif");
   s.push_back("/data/monomers/g/GLY.cif");
   s.push_back("/data/monomers/g/GLU.cif");
   s.push_back("/data/monomers/p/PHE.cif");
   s.push_back("/data/monomers/h/HIS.cif");
   s.push_back("/data/monomers/i/ILE.cif");
   s.push_back("/data/monomers/l/LYS.cif");
   s.push_back("/data/monomers/l/LEU.cif");
   s.push_back("/data/monomers/m/MET.cif");
   s.push_back("/data/monomers/p/PRO.cif");
   s.push_back("/data/monomers/a/ARG.cif");
   s.push_back("/data/monomers/s/SER.cif");
   s.push_back("/data/monomers/t/THR.cif");
   s.push_back("/data/monomers/v/VAL.cif");
   s.push_back("/data/monomers/t/TRP.cif");
   s.push_back("/data/monomers/t/TYR.cif");   

   return s;
}

int
main(int argc, char **argv) {
 
   std::string filename;
   int read_number = 1;

   if (argc <= 1) {
      
      char *s = getenv("CCP4_LIB");
      // s = ".";
      if (s) {
	 filename = s;
	 // contains the linkages:
	 filename += "/data/monomers/list/mon_lib_list.cif";
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

