/* db-main/test-dbmain.cc
 * 
 * Copyright 2005 The University of York
 * Copyright 2008 The University of Oxford
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

class generate_fragments_t {
public:
   coot::minimol::fragment five_residues_as_cas;
   clipper::Coord_orth ox_pos;
   bool have_oxt;
};

generate_fragments_t
make_target_5_res_frag(const coot::minimol::fragment &tf, int ires_start) {

   coot::minimol::fragment f;
   generate_fragments_t g;
   g.have_oxt = 0;
   for (int ii=0; ii<5; ii++) {
      int ires = ii+ires_start;
      coot::minimol::residue r(ires);
      int natoms = tf[ires].atoms.size();
      // std::cout << tf[ires] << std::endl;
      for (int iat=0; iat<natoms; iat++) {
	 // std::cout << tf[ires][iat] << std::endl;
	 if (tf[ires][iat].name == " CA ") {
	    r.addatom(tf[ires][iat]);
	 }
	 if (ii==2) { 
	    if (tf[ires][iat].name == " O  ") {
	       g.ox_pos = tf[ires][iat].pos;
	       g.have_oxt = 1;
	    }
	 }
      }
//       std::cout << " in make_target_5_res_frag adding residue with seqnum: "
// 		<< r.seqnum << std::endl;
      try { 
	 f.addresidue(r, 0);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: make_target_5_res_frag() " << rte.what() << std::endl;
      } 
   }
   g.five_residues_as_cas = f;
   return g;
}

int 
main(int argc, char **argv) { 

   int ilength = 5;
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
      
   coot::db_main main_chain;
   main_chain.fill_with_fragments(ilength);

   for (unsigned int ifrag=0; ifrag<target_all_coords.fragments.size(); ifrag++) {

      if (target_all_coords[ifrag].fragment_id == chain_id) { 

	 for (int ires=iresno_start; (ires+5)<=iresno_end; ires++) {
	    std::cout << "---------------------------------------" << std::endl;
	    generate_fragments_t frag_info =
	       make_target_5_res_frag(target_all_coords[ifrag], ires);
	    
	    if (frag_info.five_residues_as_cas.n_filled_residues() > 0) {

	       if (frag_info.have_oxt) { 

		  // load up an internal vector of fragments that have been
		  // rtop onto this target_ca_coords_5_res_frag
		  //
// 		  std::cout << "Inner Matching to target starting at " << ires
// 			    << std::endl;
		  main_chain.match_targets_for_pepflip(frag_info.five_residues_as_cas);

		  float cutoff = 2.5;
		  // check that internal vector against the oxygen position
		  float frac = main_chain.mid_oxt_outliers(frag_info.ox_pos, ires+2, cutoff);
// 		  std::cout << "Fraction better: " << 100.0*frac
// 			    << "% have oxygen dist closer than " << cutoff << std::endl;
		  main_chain.clear_results();
	       }
	    }
	 }
      }
   }
   return 0; 
}
