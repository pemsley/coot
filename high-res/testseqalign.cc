
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
