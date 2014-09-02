
#ifndef C_INTERFACE_SEQUENCE_HH
#define C_INTERFACE_SEQUENCE_HH

#include<vector>
#include<string>
#include "coot-align.hh"

class sequence_to_chain_results_t {
public:
   sequence_to_chain_results_t() {
      matches = false;
      match_fraction = -1;
   }
   bool matches;
   double match_fraction;
};

std::vector<coot::chain_mutation_info_container_t>
sequence_comparison_to_chains(int imol, std::string sequence);
   
bool assign_sequence_to_best_matching_chain(std::string sequence);

bool assign_sequences_to_best_matching_chain_from_fasta(std::string fasta_file_name);

#endif // C_INTERFACE_SEQUENCE_HH
