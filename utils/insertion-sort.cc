/*
 * utils/insertion-sort.cc
 *
 * Copyright 2009 by University of Oxford
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


#include <algorithm>

#include "insertion-sort.hh"

void insertionSort(std::vector<int> &vec) {
   for (auto it = vec.begin(); it != vec.end(); it++) {
      // Searching the upper bound, i.e., first element greater than *it from beginning
      auto const insertion_point = std::upper_bound(vec.begin(), it, *it);
      std::rotate(insertion_point, it, it+1);
   }
}
