/* lidia-core/test-cod-types-similarity.cc
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
 * 02110-1301, USA
 */

// Tests the RDKit-free parts of the COD-types similarity machinery:
//   - normalize_atom_type() parity with compute-cod-types.py --level normalized
//   - binary DB loader + Jaccard search against the shipped database
//
// Usage: test-cod-types-similarity [path-to-cod-types-db-normalized.bin]

#include <iostream>
#include <iomanip>
#include <cmath>

#include "cod-types-similarity.hh"

static int n_fail = 0;

static void check_eq(const std::string &got, const std::string &expected, const std::string &what) {
   if (got != expected) {
      std::cout << "FAIL: " << what << "\n   got:      " << got
                << "\n   expected: " << expected << std::endl;
      n_fail++;
   }
}

int main(int argc, char **argv) {

   // 1. normalize_atom_type() parity. Ground-truth pairs were produced by
   //    acedrg-similarity.py/compute-cod-types.py simplify_atom_type(level='normalized').
   const std::vector<std::pair<std::string, std::string> > pairs = {
      {"C[6a](C[6a]C[6a]H)2(COO){1|C<3>,2|H<1>}",                  "C[6a](CCH)2(COO)"},
      {"C[6a](C[6a]C[6a]2)(C[6a]C[6a]H)(H){1|Cl<1>,1|H<1>,3|C<3>}","C[6a](CC)(CCH)(H)"},
      {"Cl(C[6a]C[6a]2)",                                          "Cl(CC)"},
      {"O(CC[6a]O)",                                               "O(CCO)"},
      {"C(C[6a]C[6a]2)(O)2",                                       "C(CC)(O)2"},
      {"N(C[5a]C[5a])(CH2)(H){1|C<3>}",                            "N(CC)(CH)(H)"},
      {"H",                                                        "H"},
      {"C[6a](C[6a]C[6a]H)2(Cl){1|C<3>,2|H<1>}",                   "C[6a](CCH)2(Cl)"},
      {"S(CC)(=O)(=O)",                                            "S(CC)(O)(O)"},
      {"P(OO)(OO)(=O)",                                            "P(O)(OO)(OO)"},
   };
   for (const auto &p : pairs)
      check_eq(cod::normalize_atom_type(p.first), p.second, "normalize " + p.first);

   std::cout << "normalize_atom_type: " << pairs.size() << " cases checked" << std::endl;

   // 2. DB load + search (optional - only if a DB path is given).
   if (argc > 1) {
      std::string db_path = argv[1];
      try {
         cod::cod_types_db_t db = cod::load_cod_types_db_binary(db_path);
         std::cout << "Loaded DB: " << db.size() << " molecules from " << db_path << std::endl;

         // The A1E4R query (normalized, 9 non-H types) - see acedrg-similarity.py
         // search A1E4R-rest.cif ... --level normalized, whose top hit is 174 @ 0.7778.
         std::set<std::string> query = {
            "C(CC)(O)2",
            "C[6a](CC)(CCH)(H)",
            "C[6a](CC)(CCH)2",
            "C[6a](CCC)(CCH)(H)",
            "C[6a](CCCl)(CCH)(H)",
            "C[6a](CCH)2(COO)",
            "C[6a](CCH)2(Cl)",
            "Cl(CC)",
            "O(CCO)",
         };
         // Tversky alpha=1, beta=0.2 (the default). For this query the top hit
         // 174 is fully contained in the query (|B-A|=0), so the beta term is 0
         // and the score equals 7/9 = 0.7778.
         std::vector<cod::cod_match_t> hits = cod::search(query, db, 0.35);
         if (hits.empty()) {
            std::cout << "FAIL: search returned no hits" << std::endl;
            n_fail++;
         } else {
            const cod::cod_match_t &top = hits[0];
            std::cout << "Top hit: " << top.comp_id << "  score=" << std::fixed
                      << std::setprecision(4) << top.score
                      << "  common=" << top.n_common << std::endl;
            if (top.comp_id != "174") {
               std::cout << "FAIL: expected top hit 174, got " << top.comp_id << std::endl;
               n_fail++;
            }
            if (std::fabs(top.score - 0.7778) > 0.0005) {
               std::cout << "FAIL: expected top score ~0.7778, got " << top.score << std::endl;
               n_fail++;
            }
         }
      }
      catch (const std::exception &e) {
         std::cout << "FAIL: DB load/search threw: " << e.what() << std::endl;
         n_fail++;
      }
   } else {
      std::cout << "(no DB path given - skipping load/search test)" << std::endl;
   }

   if (n_fail == 0) {
      std::cout << "ALL TESTS PASSED" << std::endl;
      return 0;
   }
   std::cout << n_fail << " TEST(S) FAILED" << std::endl;
   return 1;
}
