/*
 * coot-utils/atom-vertex.hh
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
#ifndef ATOM_VERTEX_HH
#define ATOM_VERTEX_HH

#include "mini-mol/atom-quads.hh"

namespace coot {
   class atom_vertex {

   public:
      enum connection_type_t { START, END, STANDARD, NONE };
      connection_type_t connection_type;
      std::vector<int> forward;
      std::vector<int> backward;
      std::pair<bool,atom_index_quad> torsion_quad;
      atom_vertex() {
	 connection_type = NONE;
	 torsion_quad.first = 0;
      }
   };
}


#endif // ATOM_VERTEX_HH
