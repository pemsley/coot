
#ifndef C_INTERFACE_SEQUENCE_HH
#define C_INTERFACE_SEQUENCE_HH

#include<vector>
#include<string>
#include "coot-align.hh"
#include "utils/coot-fasta.hh"

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

void apply_fasta_multi_to_fragment(int imol, const std::string &chain_id, int resno_start, int resno_end, int imol_map,
                                   const coot::fasta_multi &fam);

void assign_sequence_to_active_fragment();

#endif // C_INTERFACE_SEQUENCE_HH
