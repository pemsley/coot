/*
 * coot-utils/map-index.hh
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

#ifndef COOT_MAP_INDEX_HH
#define COOT_MAP_INDEX_HH

namespace coot {
   class map_index_t {
      int index_;
   public:
      enum index_type { UNASSIGNED = -1 };
      map_index_t() { index_ = UNASSIGNED; }
      explicit map_index_t(int i) { index_ = i; }
      int index() const { return index_; }
      bool is_assigned() const { return (index_ != UNASSIGNED); }
      bool operator==(const map_index_t &ti) const {
	 return (ti.index() == index_);
      }
   };

}


#endif // MAP_INDEX_HH
