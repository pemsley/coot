/*
 * utils/split-indices.hh
 *
 * Copyright 2018 by Medical Research Council
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


#include <vector>

namespace coot {

   // split and fill indices
   void split_indices(std::vector<std::vector<unsigned int> > *indices,
		      unsigned int n_items,
		      unsigned int n_threads);

   std::vector<std::pair<unsigned int, unsigned int> >
   atom_index_ranges(unsigned int n_atoms, unsigned int n_threads);

}
