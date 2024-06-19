/* coot-utils/residue-and-atom-specs.hh
  * 
 * Copyright 2013 by Medical Research Council
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

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <clipper/ccp4/ccp4_map_io.h>
#include "atom-selection-container.hh"
#include "xmap-stats.hh"
#include "q-score.hh"

int main(int argc, char **argv) {

   if (argc > 2) {
      std::string pdb_file_name(argv[1]);
      std::string map_file_name(argv[2]);
      clipper::CCP4MAPfile file;
      try {
         file.open_read(map_file_name);
         clipper::Grid_sampling fgs = file.grid_sampling();
         clipper::Xmap<float> xmap;
         file.import_xmap(xmap);
         bool use_gemmi = false;
         atom_selection_container_t asc = get_atom_selection(pdb_file_name, use_gemmi);
         if (asc.read_success) {
            mean_and_variance<float> mv = map_density_distribution(xmap, 1000, false, true);

            auto tp_0 = std::chrono::high_resolution_clock::now();
            coot::q_score_t q_score(asc.mol);
            q_score.calc(xmap, mv.mean, std::sqrt(mv.variance));
            q_score.close();
            auto tp_1 = std::chrono::high_resolution_clock::now();
            auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
            std::cout << "# Timing " << d10 << " milliseconds" << std::endl;
         }
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
      }
   }
   return 0;
}
