/*
 * coot-utils/fragment-container.hh
 *
 * Copyright 2021 by Medical Research Council
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
#ifndef FRAGMENT_INFO_HH
#define FRAGMENT_INFO_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class fragment_container_t {
   public:

      class fragment_range_t {
      public:
         std::vector<mmdb::Residue *> residues; // ordered residues in range
         std::string chain_id;
         residue_spec_t start_res;
         residue_spec_t end_res;
         fragment_range_t(const std::string &chain_id_in,
                          const residue_spec_t &r1, const residue_spec_t &r2) : chain_id(chain_id_in), start_res(r1), end_res(r2) {}
         // friend std::ostream &operator<<(std::ostream &s, const fragment_range_t &fc); // see below
      };
      // fix this another time!
      // std::ostream& operator<<(std::ostream &s, const fragment_range_t &fc);

      std::vector<fragment_range_t> ranges;
      fragment_container_t() {}
      void add(const fragment_range_t &r) {
         ranges.push_back(r);
      }
      void print_fragments() const;
   };

   // "fragment" here means whole chains, or parts of chains if there are chain breaks.
   fragment_container_t make_fragments(mmdb::Manager *mol);

   // "fragment" in this case means a run of residues of the given size. There will be lots of fragments
   // typically and they will be overlapping.
   //
   fragment_container_t make_overlapping_fragments(mmdb::Manager *mol, const std::string &chain_id, unsigned int fragment_length);

}



#endif // FRAGMENT_INFO_HH
