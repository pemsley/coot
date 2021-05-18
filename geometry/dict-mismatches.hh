/* geometry/dict-mismatches.cc
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

#ifndef DICT_MISMATCHES_HH
#define DICT_MISMATCHES_HH

#include<string>
#include<math.h>

namespace coot {
   class bond_mismatch_t {
   public:
      std::string atom_id_1;
      std::string atom_id_2;
      double dist_1;
      double dist_2;
      double abs_diff;
      double diff;
      bond_mismatch_t(const std::string &a1,
                      const std::string &a2,
                      double d1, double d2) : atom_id_1(a1), atom_id_2(a2) {
         dist_1 = d1;
         dist_2 = d2;
         diff = dist_2 - dist_1;
         abs_diff = fabs(diff);
      }
      bool operator<(const bond_mismatch_t &m) const {
         return m.abs_diff < abs_diff;
      }
   };

   class angle_mismatch_t {
   public:
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      double angle_1, angle_2;
      double abs_diff;
      double diff;
      angle_mismatch_t(const std::string &a1,
                       const std::string &a2,
                       const std::string &a3,
                       double a_1, double a_2)  : atom_id_1(a1), atom_id_2(a2), atom_id_3(a3) {
         angle_1 = a_1;
         angle_2 = a_2;
         diff = angle_2 - angle_1;
         abs_diff = fabs(diff);
      }
      bool operator<(const angle_mismatch_t &m) const {
         return m.abs_diff < abs_diff;
      }
   };
}

#endif // DICT_MISMATCHES_HH
