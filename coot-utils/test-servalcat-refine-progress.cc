/* coot-utils/test-servalcat-refine-progress.cc
 *
 * Copyright 2026 by Medical Research Council
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
 * 02110-1301, USA.
 */

#include <iostream>
#include <iomanip>
#include "servalcat-refine-progress.hh"

int main(int argc, char **argv) {

   // default to the x-ray stats file written by a servalcat test run
   std::string file_name = "coot-servalcat/thing_stats.json";
   if (argc > 1) file_name = argv[1];

   coot::servalcat_refine_progress_t progress;
   bool ok = progress.parse_json_file(file_name);
   if (! ok) {
      std::cout << "FAIL: could not parse " << file_name << std::endl;
      std::cout << "Usage: test-servalcat-refine-progress <stats-json-file>" << std::endl;
      return 1;
   }

   std::cout << "parsed " << progress.size() << " cycles from " << file_name << std::endl;

   std::cout << "summary keys:";
   for (const auto &k : progress.summary_keys()) std::cout << " [" << k << "]";
   std::cout << std::endl;

   std::cout << "geom outer keys:";
   for (const auto &k : progress.geom_outer_keys()) std::cout << " [" << k << "]";
   std::cout << std::endl;

   // R-factors across cycles
   std::vector<std::pair<int, double> > rwork = progress.summary_series("Rwork");
   std::vector<std::pair<int, double> > rfree = progress.summary_series("Rfree");
   std::cout << std::fixed << std::setprecision(4);
   std::cout << "cycle    Rwork    Rfree" << std::endl;
   for (std::size_t i=0; i<rwork.size(); i++) {
      double rf = (i < rfree.size()) ? rfree[i].second : 0.0;
      std::cout << "  " << std::setw(3) << rwork[i].first
                << "  " << std::setw(7) << rwork[i].second
                << "  " << std::setw(7) << rf << std::endl;
   }

   // binned R-free at the final cycle, against 1/d^2
   int li = progress.last_cycle_index();
   std::vector<std::pair<double, double> > binned_rfree = progress.binned_series(li, "Rfree");
   std::cout << "binned Rfree at final cycle: " << binned_rfree.size() << " shells" << std::endl;

   // a geometry trace
   std::vector<std::pair<int, double> > gz =
      progress.geom_series("r.m.s.Z", "Bond distances, non H");
   if (! gz.empty())
      std::cout << "r.m.s.Z (bonds, non H): " << gz.front().second
                << " -> " << gz.back().second << std::endl;

   // a couple of sanity checks
   int status = 0;
   if (progress.empty())                                 { std::cout << "FAIL: no cycles" << std::endl; status = 1; }
   if (rwork.empty())                                    { std::cout << "FAIL: no Rwork series" << std::endl; status = 1; }
   if (progress.binned_series(li, "Rfree").empty())      { std::cout << "FAIL: no binned Rfree" << std::endl; status = 1; }
   if (status == 0) std::cout << "PASS" << std::endl;
   return status;
}
