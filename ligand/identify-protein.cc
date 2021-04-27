/* ligand/identify-protein.cc
 * 
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

#include <string>
#include <iostream>
#include <iomanip>
#include <thread>
#include <numeric>
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include "coot-utils/atom-selection-container.hh"
#include "utils/coot-utils.hh"
#include "utils/coot-fasta.hh"
#include "side-chain-densities.hh"


#include "utils/coot-fasta.hh"
#include "coot-utils/fragment-container.hh"


void test_sequence(const std::string &pdb_file_name,
		   const std::string &map_file_name,
		   const std::string &multi_sequence_file_name) {

   coot::fasta_multi fam(multi_sequence_file_name);
   if (fam.size() > 0) {

      try {
         clipper::CCP4MAPfile file;
         file.open_read(map_file_name);
         clipper::Xmap<float> xmap;
         file.import_xmap(xmap);
         atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false, false);
         if (asc.read_success) {
            std::vector<coot::side_chain_densities::results_t> scores = get_fragment_sequence_scores(asc.mol, fam, xmap);

            bool print_scores_for_alignments = false;
            if (print_scores_for_alignments) {
               std::cout << "----- score size ----- " << scores.size() << std::endl;
               for (const auto &score : scores) {
                  std::cout << "   " << score.sequence_name << " " << score.sequence << " " << score.sum_score << std::endl;
               }
            }


            // search results - actually, capture the top 10 or so, so that we can see how good
            // the best is compared to the rest.
            double best_score = -9999999999999999.9;
            coot::side_chain_densities::results_t best_result;
            bool found_a_better = false;
            for (const auto &r : scores) {
               if (r.sum_score > best_score) {
                  best_score = r.sum_score;
                  best_result = r;
                  found_a_better = true;
               }
            }
            if (found_a_better) {
               std::cout << "\nBest fit and Top 30:\n"
                         << "   " << best_result.sequence << " "
                            << std::setw(9) << best_result.sum_score + 0.0 << " "
                            << best_result.sequence_name << "\n";
               std::cout << "   --------------------------------------------------------------------\n";
            }
            

            // print out the top 10

            unsigned int n_top = 30;
            std::vector<std::pair<double, int> > top_10_scores_indices(n_top, std::pair<double, int>(-999999.9, -1));
            float worst_score_in_top_10 = -999999.9;
            for (unsigned int i=0; i<scores.size(); i++) {
               auto &score(scores[i]);
               if (score.sum_score > worst_score_in_top_10) {
                  // add this somewhere in the top_10_scores_indices
                  for (auto &s : top_10_scores_indices) {
                     if (s.first < score.sum_score) {
                        s.first  = score.sum_score;
                        s.second = i;
                        // now find and set the (new) worse score in top_10_scores_indices
                        double worst_score_in_top_10_local = 999999.9;
                        auto update_worst = [&worst_score_in_top_10_local] (const std::pair<double, int> &s_inner) mutable {
                                               if (s_inner.first < worst_score_in_top_10_local) worst_score_in_top_10_local = s_inner.first; };
                        std::for_each(top_10_scores_indices.cbegin(), top_10_scores_indices.cend(), update_worst);
                        worst_score_in_top_10 = worst_score_in_top_10_local;
                        break;
                     }
                  }
               }
            }

            for (const auto &s : top_10_scores_indices)
               if (s.second != -1)
                  std::cout << "   " // << s.second << " "
                            << scores[s.second].sequence << " "
                            << std::setw(9) << scores[s.second].sum_score + 0.0 << " "
                            << scores[s.second].sequence_name << "\n";

            
         } else {
            std::cout << "WARNING:: failed to open " << pdb_file_name << std::endl;
         }
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
      }
   }

}


int main(int argc, char **argv) {

   int status = 0;

   if (argc > 3) {
      std::string map_file_name(argv[1]);
      std::string pdb_file_name(argv[2]); // poly-ALA model
      std::string multi_sequence_file_name(argv[3]);
      test_sequence(pdb_file_name, map_file_name, multi_sequence_file_name);
   } else {
      std::cout << "Usage: coot-identify-protein <map-file-name> <coordinates-file-name> <fasta-file-name>\n";
   }
   return status;
}
