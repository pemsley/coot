/* high-res/testseqassignment.cc
 * 
 * Copyright 2003, 2004  by Paul Emsley, The University of York
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

#include "mmdb-extras.h" // for atom_selection_container_t
#include "mmdb.h"
#include "coot-map-utils.hh"
#include "sequence-assignment.hh"

int
main(int argc, char **argv) {


   if (argc > 1) {
      coot::sequence_assignment::side_chain_score_t scs;
      std::string sequence_chain_id("A");
      // std::string fasta_seq("> sequence\nVADGFLKIVCSEQVFSAGVT");
      std::string fasta_seq;
      fasta_seq += ">chain_A   ";
      fasta_seq += "\n";
      fasta_seq += "DVSGTVCLSA LPPEATDTLN LIASDGPFPY SQDGVVFQNR ESVLPTQSYG YYHEYTVITP";
      fasta_seq += "\n";
      fasta_seq += "GARTRGTRRI ICGEATQEDY YTGDHYATFS LID";
      fasta_seq += "\n";

      scs.add_fasta_sequence(sequence_chain_id, fasta_seq);
      std::string pdb_filename = argv[1];

      atom_selection_container_t asc = get_atom_selection(pdb_filename, 1);
      if (asc.read_success) {

	 if (argc > 4) {

	    scs.add_fasta_sequence(sequence_chain_id, fasta_seq);

	    std::string mtz_file_name = argv[2];
	    std::string f_col = argv[3];
	    std::string phi_col = argv[4];
	    std::string w_col = "dummy";
	    short int weight_flag = 0;
	    short int is_diff_map_flag = 0;
	    clipper::Xmap<float> xmap;
	    coot::util::map_fill_from_mtz(&xmap, mtz_file_name, f_col, phi_col,
					  w_col, weight_flag, is_diff_map_flag);
	    scs.generate_scores(asc.mol, xmap);
	    
	 }
      }

   } else {
      std::cout << "Usage: " << argv[0] << " pdb-file-name mtz-file-name F-col Phi-col\n";
   }

   
   return 0;

}
