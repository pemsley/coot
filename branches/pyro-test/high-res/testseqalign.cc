/* high-res/testseqalign.cc
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

#include "sequence-assignment.hh"

int
main(int argc, char **argv) { 
   
   int r=0;
   coot::sequence_assignment::side_chain_score_t scs;

   std::string sequence_chain_id("A");
   std::string fasta_seq("> sequence\nVADGFLKIVCSEQVFSAGVT");

   scs.add_fasta_sequence(sequence_chain_id, fasta_seq);
   
   int max_resno_in_chain = 100;
   std::string fragment_id("C");
   double score;

   for (int resno=1; resno<100; resno++) { 
      for (int residue_idx = 0; residue_idx<20; residue_idx++) { 

	 score = (float(resno)/100.0) * (float(residue_idx)/20.0);
	 scs.add_score(0, fragment_id, resno, max_resno_in_chain, residue_idx, score);
      }
   }

   scs.debug();

   scs.slider();

   return r;
}
