/* coot-utils/smith-waterman.cc
 *
 * Copyright 2021 by Medical Research Council
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

#include <iostream>
#include <iomanip>
#include <algorithm>

#include "smith-waterman.hh"

#include <mmdb2/mmdb_tables.h>  // for mmdb::Get1LetterCode()
#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>

#include "utils/coot-utils.hh"
#include "analysis/stats.hh"


// get rid of this on merge.
std::string
coot_util_single_letter_to_3_letter_code(char code) {

   if (code == 'G') return std::string("GLY");
   if (code == 'A') return std::string("ALA");
   if (code == 'V') return std::string("VAL");
   if (code == 'S') return std::string("SER");
   if (code == 'N') return std::string("ASN");
   if (code == 'P') return std::string("PRO");
   if (code == 'D') return std::string("ASP");
   if (code == 'C') return std::string("CYS");
   if (code == 'Q') return std::string("GLN");
   if (code == 'E') return std::string("GLU");
   if (code == 'H') return std::string("HIS");
   if (code == 'I') return std::string("ILE");
   if (code == 'L') return std::string("LEU");
   if (code == 'K') return std::string("LYS");
   if (code == 'M') return std::string("MET");
   if (code == 'F') return std::string("PHE");
   if (code == 'T') return std::string("THR");
   if (code == 'W') return std::string("TRP");
   if (code == 'Y') return std::string("TYR");
   if (code == 'R') return std::string("ARG");

   return std::string("");
}



float sm_wat::s(char a, const std::map<std::string, std::pair<std::string, double> > &b) {

   float ss = 0.0;   // think about how bad this needs to be when accumulating log likelihoods.
   std::string tlc = coot_util_single_letter_to_3_letter_code(a);
   std::map<std::string, std::pair<std::string, double> >::const_iterator it = b.find(tlc);

   if (it != b.end())
      ss = it->second.second;

   return ss;
}

// k is at least 1.
// this is the penalty for a gap of size k
float sm_wat::W_gap_sequence(int k) {
   float gp = 3.5 + (k-1) * 0.5;
   return gp;
}

// k is at least 1.
// this is the penalty for a gap of size k
float sm_wat::W_gap_residues(int k) {
   float gp = 6.5 + (k-1) * 0.5;
   return gp;
}


float sm_wat::score_with_method_1(int seq_idx,
                                  int types_idx,
                                  const std::vector<std::vector<std::pair<bool, float> > > &H,
                                  const std::string &target_sequence,
                                  const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   char a = target_sequence[seq_idx-1];
   const std::map<std::string, std::pair<std::string, double> > &b = scored_residues[types_idx-1].second;
   float s_a_b = s(a,b);

   if (false) {
      std::string tlc = coot_util_single_letter_to_3_letter_code(a);
      mmdb::Residue *residue_p = scored_residues[types_idx-1].first;
      if (tlc == std::string(residue_p->GetResName())) {
         s_a_b = 1.0;
      } else {
         s_a_b = -0.5;
      }
   }
   float sv = H_i_j(seq_idx-1, types_idx-1, H, target_sequence, scored_residues) + s_a_b;
   // std::cout << "----    score method_1 " << sv << std::endl;
   return sv;
}

float sm_wat::score_with_method_2(int seq_idx,
                                  int types_idx,
                                  const std::vector<std::vector<std::pair<bool, float> > > &H,
                                  const std::string &target_sequence,
                                  const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   const int max_gap_size = 13;
   std::vector<float> scores(max_gap_size+1, -1000.0);
   for(int k=1; k<=max_gap_size; k++) {
      int i_minus_k = seq_idx - k; // i - k;
      int j = types_idx;
      if (i_minus_k > 0) {
         float H_score = H_i_j(i_minus_k, j, H, target_sequence, scored_residues);
         float W_score = W_gap_sequence(k);
         scores[k] = H_score - W_score;
         // std::cout << "debug score[k] k = " << k << " " << H_score << " - " << W_score << " =  " << scores[k] << std::endl;
      }
   }
   // pick the best score
   float score_best = -10000.0;
   for(int k=1; k<=max_gap_size; k++)
      if (scores[k] > score_best) score_best = scores[k];

   if (false) {
      std::cout << "method_2 scores ";
      for (int k=1; k<=max_gap_size; k++)
         std::cout << "k:" << k << " " << std::setprecision(1) << scores[k] << "  ";
      std::cout << std::endl;
      std::cout << "---- best_score method_2 " << score_best << std::endl;
   }
   return score_best;
}

float sm_wat::score_with_method_3(int seq_idx, int types_idx,
                                  const std::vector<std::vector<std::pair<bool, float> > > &H,
                                  const std::string &target_sequence,
                                  const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   int max_gap_size = 13; // one of these need only be 2
   std::vector<float> scores(max_gap_size+1, -1000.0);
   for(int k=1; k<=max_gap_size; k++) {
      int j_minus_k = types_idx - k; // j - k;
      int i = seq_idx;
      if (j_minus_k > 0) {
         float H_score = H_i_j(i, j_minus_k, H, target_sequence, scored_residues);
         float W_score = W_gap_residues(k);
         scores[k] = H_score - W_score;
      }
   }
   // pick the best score
   float score_best = -10000.0;
   for(int k=1; k<=max_gap_size; k++)
      if (scores[k] > score_best) score_best = scores[k];

   if (false) {
      std::cout << "method_3 scores ";
      for (int k=1; k<=max_gap_size; k++)
         std::cout << "k:" << k << " " << std::setprecision(1) << scores[k] << "  ";
      std::cout << std::endl;
      std::cout << "---- best_score method_3 " << score_best << std::endl;
   }
   return score_best;
}

std::pair<float, std::vector<sm_wat::cell_t> >
sm_wat::backtrack(const std::vector<std::vector<std::pair<bool, float> > > &H) {

   std::vector<cell_t> blank;
   return backtrack_others(H, blank);

}

std::pair<float, std::vector<sm_wat::cell_t> >
sm_wat::backtrack_others(const std::vector<std::vector<std::pair<bool, float> > > &H,
                         const std::vector<cell_t> &cells_from_better_traces) {

   std::vector<sm_wat::cell_t> seq_index_for_matches;

   auto as_yet_undiscovered = [cells_from_better_traces] (const cell_t &cell) {
                                 // maybe we should check for gaps also e,g  (8, -1)
                                 return (! std::any_of(cells_from_better_traces.cbegin(), cells_from_better_traces.cend(),
                                                     [cell] (const cell_t &o) { return (o==cell);}));
                              };

   auto selected_move_i_was_in_already_accepted_cells = [] (const cell_t &selected_move,
                                                            const std::vector<cell_t> others) {
                                                           return (std::any_of(others.cbegin(), others.cend(),
                                                                               [selected_move] (const cell_t &c) { return (selected_move.i == c.i); }));};
   auto selected_move_j_was_in_already_accepted_cells = [] (const cell_t &selected_move,
                                                            const std::vector<cell_t> others) {
                                                           return (std::any_of(others.cbegin(), others.cend(),
                                                                               [selected_move] (const cell_t &c) { return (selected_move.j == c.j); }));};
   // First find where the max is and store it in idx_max_value

   cell_t idx_max_value(0,0);
   float max_value = -1000.0;
   float max_value_first = max_value;
   unsigned int lim_1 = H.size();
   for (unsigned int i=0; i<lim_1; i++) {
      unsigned int lim_2 = H[i].size();
      for (unsigned int j=0; j<lim_2; j++) {
         if (H[i][j].second > max_value) {
            if (as_yet_undiscovered(cell_t(i,j))) {
               max_value = H[i][j].second;
               max_value_first = max_value;
               idx_max_value = cell_t(i,j);
            }
         }
      }
   }

   if (true)
      std::cout << "debug:: grid search found max value " << max_value << " at "
                << idx_max_value.i << " " << idx_max_value.j << std::endl;

   seq_index_for_matches.push_back(idx_max_value);

   // backtrack until we find a 0
   bool zero_found = false;

   auto sort_moves_fn = [] (const std::pair<cell_t, float> &m1,
                            const std::pair<cell_t, float> &m2) {
                           return m1.second > m2.second;
                        };

   unsigned int round_number = 0;
   bool debug = false;
   while (! zero_found) {
      round_number++;
      if (debug)
         std::cout << "debug:: loop start idx_max_value " << idx_max_value.i << " " << idx_max_value.j
                   << " value " << max_value << std::endl;

      if (idx_max_value.i > 1) {
         if (idx_max_value.j > 1) {
            const unsigned int &i = idx_max_value.i;
            const unsigned int &j = idx_max_value.j;
            cell_t cell_up(i, j-1);
            cell_t cell_left(i-1, j);
            cell_t cell_diag(i-1, j-1);
            const float &diag = H[i-1][j-1].second;
            const float &left = H[i-1][j  ].second;
            const float &up   = H[i  ][j-1].second;
            std::pair<bool, float> available_scored_up(  as_yet_undiscovered(cell_up),     up);
            std::pair<bool, float> available_scored_left(as_yet_undiscovered(cell_left), left);
            std::pair<bool, float> available_scored_diag(as_yet_undiscovered(cell_diag), diag);
            bool move_selected = false;
            cell_t selected_move;

            std::vector<std::pair<cell_t, float> > available_moves;
            if (cell_up.non_zero_indexed())
               if (available_scored_up.first)
                  if (up > 0)
                     available_moves.push_back(std::pair<cell_t, float>(cell_up, up));
            if (cell_left.non_zero_indexed())
               if (available_scored_left.first)
                  if (left > 0)
                     available_moves.push_back(std::pair<cell_t, float>(cell_left, left));
            if (cell_diag.non_zero_indexed())
               if (available_scored_diag.first)
                  if (diag > 0)
                     available_moves.push_back(std::pair<cell_t, float>(cell_diag, diag));

            if (available_moves.size() == 3) {
               std::sort(available_moves.begin(), available_moves.end(), sort_moves_fn);
               selected_move = available_moves[0].first;
               max_value = available_moves[0].second;
               move_selected = true;
            }
            
            if (available_moves.size() == 2) {
               if (available_moves[0].second > available_moves[1].second) {
                  selected_move = available_moves[0].first;
                  max_value = available_moves[0].second;
                  move_selected = true;
               } else {
                  selected_move = available_moves[1].first;
                  max_value = available_moves[1].second;
                  move_selected = true;
               }

               std::cout << "in n_availible_moves == 2 and move_selected is " << move_selected << std::endl;
            }

            if (available_moves.size() == 1) {
               selected_move = available_moves[0].first;
               max_value = available_moves[0].second;
               move_selected = true;
            }

            if (move_selected) {

               idx_max_value = selected_move;
               if (seq_index_for_matches.empty()) {
                  seq_index_for_matches.push_back(selected_move);
               } else {
                  if (selected_move_i_was_in_already_accepted_cells(selected_move, seq_index_for_matches))
                     selected_move.i = -1;
                  if (selected_move_j_was_in_already_accepted_cells(selected_move, seq_index_for_matches))
                     selected_move.j = -1;
                  seq_index_for_matches.push_back(selected_move);
               }
                  
            } else {

               if (false) {
                  std::cout << "------------ oops - no move selected - nowhere to go "
                            << "n_availible_moves " << available_moves.size() << std::endl;
                  std::cout << "up:   " << cell_up.i   << " " << cell_up.j   << "  value  " << up   << std::endl;
                  std::cout << "left: " << cell_left.i << " " << cell_left.j << "  value  " << left << std::endl;
                  std::cout << "diag: " << cell_diag.i << " " << cell_diag.j << "  value  " << diag << std::endl;
                  if (! cells_from_better_traces.empty()) {
                     std::cout << "previous" << std::endl;
                     for (const auto &c : cells_from_better_traces)
                        std::cout << "   " << c.i << " " << c.j << std::endl;
                  }
               }
               
               zero_found = true;
            }
         } else {
            zero_found = true;
         }
      } else {
         zero_found = true;
      }
   }
   
   return std::pair<float, std::vector<sm_wat::cell_t> >(max_value_first, seq_index_for_matches);
}




std::vector<std::vector<std::pair<bool, float> > >
sm_wat::construct_H(const std::string &target_sequence,
                    const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   // an array of [target_sequence.size() + 1] x [n_residues + 1]

   unsigned int n_model_residues = scored_residues.size();
   unsigned int n = target_sequence.size() + 1;
   unsigned int m = n_model_residues + 1;

   std::pair<bool, float> initial_value(false, -1.0);
   std::pair<bool, float> true_0(true, 0.0);
   std::vector<std::vector<std::pair<bool, float> > > H(n);
   for (unsigned int i=0; i<H.size(); i++)
      H[i].resize(m, initial_value);
   for (unsigned int i=0; i<H.size(); i++)
      H[i][0] = true_0;
   for (unsigned int j=0; j<m; j++)
      H[0][j] = true_0;
   std::cout << "H is constructed " << n << " " << m << std::endl;
   return H;
}

void
sm_wat::fill_scoring_matrix(std::vector<std::vector<std::pair<bool, float> > > &H,
                            const std::string &target_sequence,
                            const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   int n = target_sequence.length();
   int m = scored_residues.size();
   int lim_1 = n + 1;
   int lim_2 = m + 1;

   // std::cout << "::::::::::::::::: fill_scoring_matrix lims: " << lim_1 << " " << lim_2 << std::endl;
   for (int i=0; i<lim_1; i++) {
      for (int j=0; j<lim_2; j++) {
         // std::cout << "calc H for " << i << " " << j << std::endl;
         float v = H_i_j(i, j, H, target_sequence, scored_residues);
         // std::cout << "filling H " << i << " " << j << " with " << std::setprecision(4) << v << std::endl;
         H[i][j] = std::pair<bool, float> (true, v);
      }
   }
   // std::cout << "H is filled" << std::endl;
}

float
sm_wat::H_i_j(int seq_idx, int residue_idx,
              const std::vector<std::vector<std::pair<bool, float> > > &H,
              const std::string &target_sequence,
              const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   // std::cout << "H_i_j called with " << seq_idx << " " << types_idx << std::endl;

   const int &i =     seq_idx;
   const int &j = residue_idx;
   if (seq_idx     == 0) return 0.0;
   if (residue_idx == 0) return 0.0;

   if (H[i][j].first) return H[i][j].second;

   float offset = 0.0;

   float v_1 = score_with_method_1(i, j, H, target_sequence, scored_residues) + offset;
   float v_2 = score_with_method_2(i, j, H, target_sequence, scored_residues) + offset;
   float v_3 = score_with_method_3(i, j, H, target_sequence, scored_residues) + offset;
   float v_4 = 0.0;
   float v = v_1;
   if (v_2 > v) v = v_2;
   if (v_3 > v) v = v_3;
   if (v_4 > v) v = v_4;

   // std::cout << "debug:: H " << i << " " << j << " choose between " << v_1 << " " << v_2 << " " << v_3 << " " << v_4 << std::endl;

   return v;
}

void sm_wat::print_H(const std::vector<std::vector<std::pair<bool, float> > > &H) {

   unsigned int lim_1 = H.size();

   if (true) {
      std::cout << "booleans" << std::endl;
      for (unsigned int i=0; i<lim_1; i++) {
         unsigned int lim_2 = H[i].size();
         for (unsigned int j=0; j<lim_2; j++) {
            std::cout << H[i][j].first << "  ";
         }
         std::cout << std::endl;
      }
   }

   std::cout << "values" << std::fixed << std::endl;
   for (unsigned int i=0; i<lim_1; i++) {
      unsigned int lim_2 = H[i].size();
      for (unsigned int j=0; j<lim_2; j++) {
         std::cout << std::setw(5) << std::setprecision(5) << std::left << H[i][j].second << " ";
      }
      std::cout << std::endl;
   }
}


void
sm_wat::print_alignment(const std::vector<cell_t> &indexed_sequences_in, // currently reverse order
                        const std::string &sequence,
                        const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   // std::cout << "in print_alignment() sequence length      " << sequence.length() << std::endl;
   // std::cout << "in print_alignment() scored residues size " << scored_residues.size() << std::endl;

   auto indexed_sequences = indexed_sequences_in;
   std::reverse(indexed_sequences.begin(), indexed_sequences.end());
   std::string seq_align_string;
   std::string res_align_string;
   std::string res_align_string_slc;
   for (const auto &p : indexed_sequences) {
      int idx_for_seq = p.i - 1;
      int idx_for_res = p.j - 1;
      if (p.i == -1) {
         seq_align_string += "-";
      } else {
         // happy path
         char c = sequence[idx_for_seq];
         seq_align_string += c;
      }
      if (p.j == -1) {
         res_align_string += "--- ";
         res_align_string_slc += "-";
      } else {
         mmdb::Residue *res = scored_residues[idx_for_res].first;;
         std::string res_name = res->GetResName();
         res_align_string += res_name;
         res_align_string += " ";
         char r[10];
         mmdb::Get1LetterCode(res->GetResName(), r);
         res_align_string_slc += r[0];
      }
   }
   std::cout << "   input sequence: " << sequence << std::endl;
   std::cout << "      aligned target-seq: " << seq_align_string << std::endl;
   std::cout << "   aligned current-model: " << res_align_string_slc << std::endl;

   if (false) // not very useful
      std::cout << "   model: " << res_align_string << std::endl;
}


std::vector<sm_wat::cell_t> // return the alignment
sm_wat::smith_waterman(const std::string &sequence,
                       const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues_in) {

   // because I want to fiddle with them for debugging
   std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > scored_residues = scored_residues_in;

   coot::stats::single data_pre;
   coot::stats::single data_post;

   for (const auto &scored_residue : scored_residues) {
      mmdb::Residue *residue_p = scored_residue.first;
      std::string actual_residue_name(residue_p->GetResName());
      for (const auto &scored_type : scored_residue.second) {
         data_pre.add(scored_type.second.second);
      }
   }

   float squidge_offset = -data_pre.mean() + 1.2 * sqrt(data_pre.variance()); // 1.1 may need some adjustment, related to gap penalty
   float squidge_scale = 1.0/sqrt(data_pre.variance());

   for (auto &scored_residue : scored_residues) {
      mmdb::Residue *residue_p = scored_residue.first;
      std::string actual_residue_name(residue_p->GetResName());
      // std::cout << "Residue " << residue_p->GetChainID() << " " << residue_p->GetSeqNum() << " "
      // << residue_p->GetResName() << std::endl;
      for (auto &scored_type : scored_residue.second) {
         scored_type.second.second += squidge_offset;
         scored_type.second.second *= squidge_scale;
         data_post.add(scored_type.second.second);
         if (scored_type.first == actual_residue_name)
            scored_type.second.second += 0.0;
      }
   }

   std::cout << "--------- statistics ---------" << std::endl;
   std::cout << "   data_pre:  mean: " << data_pre.mean()  << " sd: " << sqrtf(data_pre.variance()) << std::endl;
   std::cout << "   data_post: mean: " << data_post.mean() << " sd: " << sqrt(data_post.variance()) << std::endl;

   bool print_the_residue_scores = false;
   
   if (print_the_residue_scores) { // print the residue scores
      for (const auto &scored_residue : scored_residues) {
         mmdb::Residue *residue_p = scored_residue.first;
         std::cout << "Residue " << residue_p->GetChainID() << " " << residue_p->GetSeqNum()
                   << " " << residue_p->GetResName() << std::endl;
         for (const auto &scored_type : scored_residue.second) {
            std::cout << "   " << scored_type.first << " " << std::setw(7) << std::setprecision(4) << std::fixed << std::right << scored_type.second.second
                      << std::endl;
         }
      }
   }

   std::vector<std::vector<std::pair<bool, float> > > H = construct_H(sequence, scored_residues);
   // print_H(H);
   fill_scoring_matrix(H, sequence, scored_residues);
   // print_H(H);

   // these indices are for residue numbers, not indexing into the string.
   // So we need to offset by -1 get the correct single-letter-code from the sequence
   // and the residue from the scored_residues (or run-of-residues)
   //
   std::pair<float, std::vector<cell_t> > scored_trace   = sm_wat::backtrack(H);

   // std::cout << "Top Score " << scored_trace.first << std::endl;

   bool print_sequence_indices = true; // debugging really.
   if (print_sequence_indices) {
      std::cout << "Best indexes:" << std::endl;
      std::cout << "sequence indices (rev): ";
      for (unsigned int i=0; i<scored_trace.second.size(); i++)
         std::cout << std::setw(2) << std::right << scored_trace.second[i].i << " ";
      std::cout << std::endl;
      std::cout << "residue  indices (rev): ";
      for (unsigned int j=0; j<scored_trace.second.size(); j++)
         std::cout << std::setw(2) << std::right << scored_trace.second[j].j << " ";
      std::cout << std::endl;
   }

   print_alignment(scored_trace.second, sequence, scored_residues);

#if 0
   std::cout << " " << std::endl;

   std::pair<float, std::vector<cell_t> > scored_trace_2 = sm_wat::backtrack_others(H, scored_trace.second);
   std::cout << "Next Best " << scored_trace_2.first << std::endl;
   std::cout << "Second Best indexes:" << std::endl;
   std::cout << "sequence indices (rev): ";
   for (unsigned int i=0; i<scored_trace_2.second.size(); i++)
      std::cout << std::setw(2) << std::right << scored_trace_2.second[i].i << " ";
   std::cout << std::endl;
   std::cout << "residue  indices (rev): ";
   for (unsigned int j=0; j<scored_trace_2.second.size(); j++)
      std::cout << std::setw(2) << std::right << scored_trace_2.second[j].j << " ";
   std::cout << std::endl;

   print_alignment(scored_trace_2.second, sequence, scored_residues);
#endif

   return scored_trace.second;

}

#include <ligand/side-chain-densities.hh>


std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > >
sm_wat::get_side_chain_density_scores_for_residues(const std::vector<mmdb::Residue *> &a_run_of_residues,
                                           const clipper::Xmap<float> &xmap) {
   
   std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > v;

   std::cout << "------------------- a_run_of_residues has " << a_run_of_residues.size() << " residues" << std::endl;

   if (! a_run_of_residues.empty()) {

      coot::side_chain_densities scd;

      scd.fill_residue_blocks(a_run_of_residues, xmap); // return fast if already filled, uses atomic locking.
      int n_residues = a_run_of_residues.size();
      for (int i=0; i<n_residues; i++) {
         mmdb::Residue *residue_p = a_run_of_residues[i];
         std::map<std::string, std::pair<std::string, double> > likelihood_map =
            scd.likelihood_of_each_rotamer_at_this_residue(residue_p, xmap);
         std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > p(residue_p, likelihood_map);
         v.push_back(p);
      }
   }
   return v;
}

void
sm_wat::apply_alignment_to_model(const std::vector<sm_wat::cell_t> &alignment,
                                 const std::string &target_sequence,
                                 const std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > &scored_residues) {

   // mutate away!
   for (unsigned int i=0; i<scored_residues.size(); i++) {
      std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > scored_residue = scored_residues[i];
      mmdb::Residue *residue_p = scored_residue.first;
      std::string current_res_type(residue_p->GetResName());
      bool mutated_this = false;
      int idx_residue_seq = i + 1;
      for (unsigned int j=0; j<alignment.size(); j++) {
         if (alignment[j].j == idx_residue_seq) {
            int raw_index = alignment[j].i-1;
            if (raw_index >= 0) {
               char c = target_sequence[raw_index];
               // std::cout << "for j " << j << " raw index is " << alignment[j].i-1 << " c is '" << c << "'" << std::endl;
               if (c != '=') {
                  std::string new_res_name = coot_util_single_letter_to_3_letter_code(c);
                  if (false)
                     std::cout << "deubg:: mutate " << coot::residue_spec_t(residue_p) << " " << current_res_type << " to "
                               << "\"" << new_res_name << "\" from char c " << c << std::endl;
                  coot::util::mutate(residue_p, new_res_name);
                  mutated_this = true;
                  break;
               }
            }
         }
      }
      if (mutated_this == false) {
         std::cout << "Missed mutation for residue " << coot::residue_spec_t(residue_p) << std::endl;
      }
   }
}

#include <coot-utils/fragment-container.hh>
#include <ligand/backrub-rotamer.hh>

void
sm_wat::align_and_mutate_and_backrub(mmdb::Manager *mol, const std::string &seq, const clipper::Xmap<float> &xmap,
                                     const coot::protein_geometry &pg) {

   coot::fragment_container_t fragments = coot::make_fragments(mol);
   fragments.print_fragments();

   for (const auto &fragment : fragments.ranges) {

      std::cout << "----------------- fragment has " << fragment.residues.size() << " residues " << std::endl;
      std::cout << "----------------- fragment: " << std::endl;
      std::vector<std::pair<mmdb::Residue *, std::map<std::string, std::pair<std::string, double> > > > scored_residues =
         get_side_chain_density_scores_for_residues(fragment.residues, xmap);

      std::cout << "-------------------- we got scored_residues of size " << scored_residues.size() << std::endl;

      if (! scored_residues.empty()) {
         const std::vector<sm_wat::cell_t> &alignment = sm_wat::smith_waterman(seq, scored_residues);
         apply_alignment_to_model(alignment, seq, scored_residues);
         coot::backrub_molecule(mol, &xmap, pg);
      }
   }
}
