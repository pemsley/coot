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
#include <thread>
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/fragment-container.hh"
#include "utils/coot-utils.hh"
#include "utils/coot-fasta.hh"
#include "side-chain-densities.hh"
#include "utils/split-indices.hh"

void test_sequence(const std::string &pdb_file_name,
		   const std::string &map_file_name,
		   const std::string &multi_sequence_file_name) {

   auto proc_threads = [] (const std::pair<unsigned int, unsigned int> &start_stop_pair,
                           const coot::fasta_multi &fam,
                           const std::vector<mmdb::Residue *> &a_run_of_residues,
                           const clipper::Xmap<float> &xmap,
                           coot::side_chain_densities &scd) { // fill scd

                          for(unsigned int idx=start_stop_pair.first; idx!=start_stop_pair.second; ++idx) {
                             scd.test_sequence(a_run_of_residues, xmap, fam[idx].name, fam[idx].sequence);
                          }
                       };

   coot::fasta_multi fam(multi_sequence_file_name);
   if (fam.size() > 0) {

      unsigned int n_sequences = fam.size();
      try {
         clipper::CCP4MAPfile file;
         file.open_read(map_file_name);
         clipper::Xmap<float> xmap;
         file.import_xmap(xmap);
         atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false, false);
         if (asc.read_success) {
            // "analysis" constructor
            // coot::side_chain_densities scd(n_steps, grid_box_radius, useable_grid_points_file_name);
            // scd.set_data_dir("side-chain-data");
            coot::side_chain_densities scd;

            coot::fragment_container_t fc = coot::make_fragments(asc.mol);

            std::cout << "INFO:: n_sequences: " << n_sequences << std::endl;

            for (const auto &range : fc.ranges) {
               std::vector<mmdb::Residue *> a_run_of_residues =
                  scd.setup_test_sequence(asc.mol, range.chain_id, range.start_res.res_no, range.end_res.res_no, xmap);

#if 1 // threaded version - reinstate when multisequence comparison is fast
               unsigned int n_threads = coot::get_max_number_of_threads();
               std::vector<std::pair<unsigned int, unsigned int> > seq_index_vector =
                  coot::atom_index_ranges(n_sequences, n_threads);
               std::vector<std::thread> threads;

               for (unsigned int i=0; i<seq_index_vector.size(); i++) {
                  std::pair<unsigned int, unsigned int> index_pair = seq_index_vector[i];
                  threads.push_back(std::thread(proc_threads, index_pair, fam, a_run_of_residues, xmap, std::ref(scd)));
               }

               for (unsigned int i=0; i<seq_index_vector.size(); i++)
                  threads[i].join();
#endif

#if 0 // the single threaded way

               for (unsigned int idx=0; idx<n_sequences; idx++) {
                  std::string sequence = fam[idx].sequence;
                  // std::cout << "Input Sequence:\n" << sequence << std::endl;
                  const std::string &name = fam[idx].name;
                  scd.test_sequence(a_run_of_residues, xmap, name, sequence);
               }
#endif

               bool  print_results = false;
               std::map<std::string, std::vector<coot::side_chain_densities::results_t> >::const_iterator it;
               if (print_results) {
                  for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
                     const std::string &sequence = it->first;
                     std::cout << "sequence: " << sequence << std::endl;
                     for (const auto &r : it->second) {
                        std::cout << "   " << r.offset << " " << r.sequence << " " << r.sum_score << std::endl;
                     }
                  }
               }
               // search results - actually, capture the top 10 or so, so that we can see how good
               // the best is compared to the rest.
               double best_score = -9999999999999999.9;
               coot::side_chain_densities::results_t best_result;
               bool found_a_better = false;
               for (it=scd.results_container.begin(); it!=scd.results_container.end(); ++it) {
                  for (const auto &r : it->second) {
                     if (r.sum_score > best_score) {
                        best_score = r.sum_score;
                        best_result = r;
                        found_a_better = true;
                     }
                  }
               }
               if (found_a_better) {
                  std::cout << "Best fit: " << best_result.sequence_name << " "
                            << best_result.offset << " " << best_result.sequence << " " << best_result.sum_score
                            << std::endl;
               }
            }
         }
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << pdb_file_name << std::endl;
         std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
      }
   }

   // Now, what's best in fam?
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
