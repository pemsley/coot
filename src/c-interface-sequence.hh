/*
 * src/c-interface-sequence.hh
 *
 * Copyright 2013 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

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
