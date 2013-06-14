
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <gtk/gtk.h>

#include "c-interface-sequence.hh"
#include "coot-align.hh"

#include "c-interface.h"
#include "graphics-info.h"

std::vector<coot::chain_mutation_info_container_t>
sequence_comparison_to_chains(int imol, std::string sequence) {

   std::vector<coot::chain_mutation_info_container_t> r;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].sequence_comparison_to_chains(sequence);
   }
   return r;
}

bool
assign_sequence_to_best_matching_chain(std::string sequence) {

   bool done = false;
   graphics_info_t g;
   std::vector<std::pair<int, std::vector<coot::chain_mutation_info_container_t> > > r;

   int best_match_mol = -1;
   std::string best_match_chain;
   // the lower the dissimilarity_score the better
   double best_match_score = -1;
   
   for (unsigned int imol=0; imol<graphics_n_molecules(); imol++) { 
      if (is_valid_model_molecule(imol)) {
	 std::vector<coot::chain_mutation_info_container_t> matches =
	    g.molecules[imol].sequence_comparison_to_chains(sequence);

	 if (0) // debug
	    for (unsigned int im=0; im<matches.size(); im++)
	       std::cout << "    \"" << matches[im].chain_id << "\" dissimilarity_score: "
			 << matches[im].dissimilarity_score() << " alignment score: "
			 << matches[im].alignment_score.first << " "
			 << matches[im].alignment_score.second << " "
			 << std::endl;

	 for (unsigned int im=0; im<matches.size(); im++) {
	    if (matches[im].alignment_score.first) { 
	       double score = matches[im].alignment_score.second;
	       if (score > best_match_score) {
		  best_match_score = score;
		  best_match_mol = imol;
		  best_match_chain = matches[im].chain_id;
	       }
	    }
	 }
      }
   }

   if (is_valid_model_molecule(best_match_mol)) {
      // assign to NCS-related molecules too if it finds a match
      g.molecules[best_match_mol].assign_sequence_to_NCS_related_chains_from_string(best_match_chain, sequence);
      done = true;
   }

   return done;
}

bool
assign_sequences_to_best_matching_chain_from_fasta(std::string fasta_file_name) {

   bool done = false;
   coot::fasta_multi f(fasta_file_name);
   
   for (unsigned int i=0; i<f.size(); i++) { 
      bool d = assign_sequence_to_best_matching_chain(f[i].sequence);
      if (d)
	 done = true;
   }
   return done;
}



   
