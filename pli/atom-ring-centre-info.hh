/* pli/atom-ring-centre-info.hh
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

#ifndef ATOM_RING_CENTRE_INFO_HH
#define ATOM_RING_CENTRE_INFO_HH

#include "lidia-core/lig-build.hh"

// trivial container for a (copy of an) atom an its ring centre (if
// it has one)
class atom_ring_centre_info_t {
public:
   lig_build::atom_t atom;
   bool has_ring_centre_flag;
   lig_build::pos_t ring_centre;
   explicit atom_ring_centre_info_t(const lig_build::atom_t &at) : atom(at) {
      has_ring_centre_flag = false;
   }
   void add_ring_centre(const lig_build::pos_t &pos) {
      ring_centre = pos;
      has_ring_centre_flag = true;
   }
};
std::ostream& operator<<(std::ostream &s, atom_ring_centre_info_t wa);



#endif // ATOM_RING_CENTRE_INFO_HH
