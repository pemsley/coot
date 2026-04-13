/*
 * analysis/test-daca.cc
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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
#include <cmath>
#include <map>
#include <set>
#include "daca.hh"

int test_box_index_round_trip() {

   int n_fail = 0;

   // A point at (2.3, -1.7, 4.9) with box_width 1.0 should land in box (2, -2, 4)
   clipper::Coord_orth pos(2.3, -1.7, 4.9);
   coot::daca::box_index_t bi(pos);

   if (bi.idx_x != 2)  { std::cout << "FAIL: idx_x " << bi.idx_x << " expected 2"  << std::endl; n_fail++; }
   if (bi.idx_y != -2) { std::cout << "FAIL: idx_y " << bi.idx_y << " expected -2" << std::endl; n_fail++; }
   if (bi.idx_z != 4)  { std::cout << "FAIL: idx_z " << bi.idx_z << " expected 4"  << std::endl; n_fail++; }

   // coord_orth() should return the centre of that box: (2.5, -1.5, 4.5)
   clipper::Coord_orth centre = bi.coord_orth();
   double tol = 1.0e-6;
   if (std::abs(centre.x() - 2.5) > tol) { std::cout << "FAIL: centre x " << centre.x() << std::endl; n_fail++; }
   if (std::abs(centre.y() - (-1.5)) > tol) { std::cout << "FAIL: centre y " << centre.y() << std::endl; n_fail++; }
   if (std::abs(centre.z() - 4.5) > tol) { std::cout << "FAIL: centre z " << centre.z() << std::endl; n_fail++; }

   // Test the origin box
   clipper::Coord_orth origin(0.1, 0.2, 0.3);
   coot::daca::box_index_t bi_origin(origin);
   if (bi_origin.idx_x != 0) { std::cout << "FAIL: origin idx_x" << std::endl; n_fail++; }
   if (bi_origin.idx_y != 0) { std::cout << "FAIL: origin idx_y" << std::endl; n_fail++; }
   if (bi_origin.idx_z != 0) { std::cout << "FAIL: origin idx_z" << std::endl; n_fail++; }

   // Negative coordinates: floor(-0.1) = -1
   clipper::Coord_orth neg(-0.1, -0.1, -0.1);
   coot::daca::box_index_t bi_neg(neg);
   if (bi_neg.idx_x != -1) { std::cout << "FAIL: neg idx_x " << bi_neg.idx_x << std::endl; n_fail++; }
   if (bi_neg.idx_y != -1) { std::cout << "FAIL: neg idx_y " << bi_neg.idx_y << std::endl; n_fail++; }
   if (bi_neg.idx_z != -1) { std::cout << "FAIL: neg idx_z " << bi_neg.idx_z << std::endl; n_fail++; }

   if (n_fail == 0)
      std::cout << "PASS test_box_index_round_trip" << std::endl;
   return n_fail;
}

int test_box_index_d_squared() {

   int n_fail = 0;
   double tol = 0.01;

   // Box (1, 2, 3) has centre (1.5, 2.5, 3.5), d_squared = 1.5^2 + 2.5^2 + 3.5^2 = 20.75
   coot::daca::box_index_t bi(1, 2, 3);
   float ds = bi.d_squared();
   float expected = 1.5f*1.5f + 2.5f*2.5f + 3.5f*3.5f;
   if (std::abs(ds - expected) > tol) {
      std::cout << "FAIL: d_squared " << ds << " expected " << expected << std::endl;
      n_fail++;
   }

   float d = bi.d();
   if (std::abs(d - std::sqrt(expected)) > tol) {
      std::cout << "FAIL: d " << d << " expected " << std::sqrt(expected) << std::endl;
      n_fail++;
   }

   if (n_fail == 0)
      std::cout << "PASS test_box_index_d_squared" << std::endl;
   return n_fail;
}

int test_box_index_ordering() {

   int n_fail = 0;

   coot::daca::box_index_t a(1, 2, 3);
   coot::daca::box_index_t b(1, 2, 4);
   coot::daca::box_index_t c(1, 3, 0);
   coot::daca::box_index_t d(2, 0, 0);
   coot::daca::box_index_t a2(1, 2, 3);

   // irreflexive: !(a < a)
   if (a < a) { std::cout << "FAIL: a < a should be false" << std::endl; n_fail++; }

   // consistent ordering: a < b < c < d
   if (!(a < b)) { std::cout << "FAIL: a < b" << std::endl; n_fail++; }
   if (!(b < c)) { std::cout << "FAIL: b < c" << std::endl; n_fail++; }
   if (!(c < d)) { std::cout << "FAIL: c < d" << std::endl; n_fail++; }
   if (!(a < d)) { std::cout << "FAIL: a < d (transitivity)" << std::endl; n_fail++; }

   // antisymmetric: a < b implies !(b < a)
   if (b < a) { std::cout << "FAIL: b < a should be false" << std::endl; n_fail++; }

   // equivalence: !(a < a2) && !(a2 < a)
   if (a < a2) { std::cout << "FAIL: a < a2 should be false" << std::endl; n_fail++; }
   if (a2 < a) { std::cout << "FAIL: a2 < a should be false" << std::endl; n_fail++; }

   // usable as std::map key
   std::map<coot::daca::box_index_t, int> m;
   m[a] = 10;
   m[b] = 20;
   m[a2] = 30; // should overwrite a
   if (m.size() != 2) { std::cout << "FAIL: map size " << m.size() << " expected 2" << std::endl; n_fail++; }
   if (m[a] != 30) { std::cout << "FAIL: map[a] " << m[a] << " expected 30" << std::endl; n_fail++; }

   if (n_fail == 0)
      std::cout << "PASS test_box_index_ordering" << std::endl;
   return n_fail;
}

int test_atom_names_for_fragments() {

   int n_fail = 0;
   coot::daca daca;

   // Expected fragment counts per residue type (including the common CA/C/O fragment)
   std::map<std::string, unsigned int> expected_counts = {
      {"GLY", 2}, {"ALA", 2}, {"CYS", 3}, {"ASP", 4}, {"GLU", 5},
      {"PHE", 4}, {"HIS", 4}, {"ILE", 4}, {"LYS", 6}, {"LEU", 4},
      {"MET", 5}, {"MSE", 5}, {"ASN", 4}, {"PRO", 2}, {"GLN", 5},
      {"ARG", 6}, {"SER", 3}, {"THR", 3}, {"VAL", 3}, {"TRP", 4},
      {"TYR", 4}
   };

   for (const auto &p : expected_counts) {
      std::vector<std::vector<std::string> > frags = daca.atom_names_for_fragments(p.first);
      if (frags.size() != p.second) {
         std::cout << "FAIL: " << p.first << " fragment count " << frags.size()
                   << " expected " << p.second << std::endl;
         n_fail++;
      }

      // Every residue should have CA/C/O as the first fragment
      if (frags.empty()) {
         std::cout << "FAIL: " << p.first << " has no fragments" << std::endl;
         n_fail++;
         continue;
      }
      const std::vector<std::string> &first_frag = frags[0];
      bool has_ca = false, has_c = false, has_o = false;
      for (const auto &name : first_frag) {
         if (name == " CA ") has_ca = true;
         if (name == " C  ") has_c = true;
         if (name == " O  ") has_o = true;
      }
      if (!has_ca || !has_c || !has_o) {
         std::cout << "FAIL: " << p.first << " first fragment missing CA/C/O" << std::endl;
         n_fail++;
      }

      // Every fragment should have at least 3 atoms (minimum for an RTop fit)
      for (unsigned int i=0; i<frags.size(); i++) {
         if (frags[i].size() < 3) {
            std::cout << "FAIL: " << p.first << " fragment " << i
                      << " has only " << frags[i].size() << " atoms (need >= 3)" << std::endl;
            n_fail++;
         }
      }

      // No duplicate atom names within a fragment
      for (unsigned int i=0; i<frags.size(); i++) {
         std::set<std::string> names(frags[i].begin(), frags[i].end());
         if (names.size() != frags[i].size()) {
            std::cout << "FAIL: " << p.first << " fragment " << i
                      << " has duplicate atom names" << std::endl;
            n_fail++;
         }
      }
   }

   // Unknown residue should return just the common CA/C/O fragment
   std::vector<std::vector<std::string> > unk = daca.atom_names_for_fragments("UNK");
   if (unk.size() != 1) {
      std::cout << "FAIL: UNK fragment count " << unk.size() << " expected 1" << std::endl;
      n_fail++;
   }

   if (n_fail == 0)
      std::cout << "PASS test_atom_names_for_fragments" << std::endl;
   return n_fail;
}

int test_gompertz_scale() {

   int n_fail = 0;
   coot::daca daca;

   // gompertz_scale takes squared distance.
   // At small distances (near origin), scale should be close to 1.0
   // At large distances, scale should approach 0.0
   // Should be monotonically decreasing with distance

   float g_near = daca.gompertz_scale(0.0f);
   if (g_near < 0.9f) {
      std::cout << "FAIL: gompertz at origin " << g_near << " expected close to 1.0" << std::endl;
      n_fail++;
   }

   float g_far = daca.gompertz_scale(100.0f); // dist = 10 Angstroms
   if (g_far > 0.5f) {
      std::cout << "FAIL: gompertz at 10A " << g_far << " expected < 0.5" << std::endl;
      n_fail++;
   }

   // Monotonicity check: increasing distance should give decreasing (or equal) scale
   float prev = daca.gompertz_scale(0.0f);
   bool monotonic = true;
   for (int i=1; i<=15; i++) {
      float dist = static_cast<float>(i);
      float dist_sq = dist * dist;
      float g = daca.gompertz_scale(dist_sq);
      if (g > prev + 1.0e-6f) {
         std::cout << "FAIL: gompertz not monotonic at dist=" << dist
                   << " g=" << g << " prev=" << prev << std::endl;
         monotonic = false;
         n_fail++;
      }
      prev = g;
   }

   // All values should be in [0, 1]
   for (int i=0; i<=20; i++) {
      float dist = static_cast<float>(i) * 0.5f;
      float dist_sq = dist * dist;
      float g = daca.gompertz_scale(dist_sq);
      if (g < 0.0f || g > 1.0f) {
         std::cout << "FAIL: gompertz out of [0,1] at dist=" << dist
                   << " g=" << g << std::endl;
         n_fail++;
      }
   }

   if (n_fail == 0)
      std::cout << "PASS test_gompertz_scale" << std::endl;
   return n_fail;
}


int main(int argc, char **argv) {

   int n_fail = 0;

   n_fail += test_box_index_round_trip();
   n_fail += test_box_index_d_squared();
   n_fail += test_box_index_ordering();
   n_fail += test_atom_names_for_fragments();
   n_fail += test_gompertz_scale();

   if (n_fail == 0)
      std::cout << "All DACA tests passed" << std::endl;
   else
      std::cout << n_fail << " DACA test(s) FAILED" << std::endl;

   return n_fail;
}
