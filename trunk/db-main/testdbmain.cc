/* db-main/test-dbmain.cc
 * 
 * Copyright 2005 The University of York
 * Author: Paul Emsley
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


#include <stdlib.h>
#include "db-main.hh"

int 
main(int argc, char **argv) { 

   int ilength = 6;
   std::string filename; 
   int iresno_start = 1;
   int iresno_end   = 46;
   std::string chain_id;

   if (argc < 5) {
      std::cout << "Usage: " << argv[0] << " [pdb_name ChainID StartResNo LastResno]"
		<< std::endl;
      exit(1); 
   } else {
      filename = argv[1];
      chain_id = argv[2];
      char **endptr = NULL; 
      iresno_start = strtol(argv[3],endptr,10);
      endptr = NULL;
      iresno_end   = strtol(argv[4],endptr,10);
   } 


   // First get a set of Ca that the library will fit.  Typically, we
   // use 1crn.pdb to provide the Ca.  Other pdbs can be used if their
   // filnename is provided.

   coot::minimol::molecule target_all_coords;
   std::cout << "Reading target coordinates: " << filename << std::endl;
   target_all_coords.read_file(filename);
   coot::minimol::molecule target_ca_coords = target_all_coords.molecule_of_atom_types(std::string(" CA "));

   if (target_all_coords.fragments.size() == 0) {
      std::cout << "Problem reading testfragment/Ca-selection thereof"
		<< std::endl;
   } else { 

      iresno_end = target_ca_coords[0].max_residue_number();
      // iresno_end = 40; 
	 
      // Now setup the library, fill it with data (coordinates)
      coot::db_main main_chain;
      int idbfrags = main_chain.fill_with_fragments(ilength);

      // could also use main_chain.is_empty() here.
      if (idbfrags > 0) { 

	 // now fill big results.
	 //
	 for(unsigned int ifrag=0; ifrag<target_ca_coords.fragments.size(); ifrag++) {
//  	    std::cout << "target_ca_coords fragment " << ifrag << " has "
//  		      << target_ca_coords[ifrag].n_residues()
//  		      << " residues " << std::endl; 
	 }

	 main_chain.match_target_fragment(target_ca_coords,
					  iresno_start,
					  iresno_end,
					  ilength);

	 main_chain.merge_fragments();

	 coot::minimol::molecule mol;

	 mol.fragments.push_back(main_chain.mainchain_fragment());

	 mol.write_file("db-mainchain.pdb", 20.0);
      }
   }
   return 0; 
}
