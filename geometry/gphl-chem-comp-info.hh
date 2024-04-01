/*
 * geometry/gphl-chem-comp-info.hh
 *
 * Copyright 2023 by Medical Research Council
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
#ifndef GPHL_CHEM_COMP_INFO_HH
#define GPHL_CHEM_COMP_INFO_HH

#include <string>
#include <vector>

namespace coot {

   class gphl_chem_comp_info_t {
   public:
      gphl_chem_comp_info_t() {}
      std::vector<std::pair<std::string, std::string> > info;
      void add(const std::string &key, const std::string &data_string) {
         auto p = std::make_pair(key, data_string);
         info.push_back(p);
      }
      int get_index(const std::string &key) const {
         int idx = -1;
         for (unsigned int i=0; i<info.size(); i++) {
            if (info[i].first == key) {
               idx = i;
               break;
            }
         }
         return idx;
      }
      // caller ensures that the idx is in range (using the above function)
      std::pair<std::string, std::string> & operator[](int idx) {
         return info[idx];
      }
   };

}

#endif // GPHL_CHEM_COMP_INFO_HH
