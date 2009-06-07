/* mini-mol/min-mol-utils.hh
 * 
 * Copyright 2003  The University of York
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

#include "mini-mol.hh"

namespace coot { 

   // C-beta position return a pair, the first of which indicates if
   // the cbeta was properly positioned.
   std::pair<short int, clipper::Coord_orth> cbeta_position(const minimol::residue &res);

   // Carbonyl O position return a pair, the first of which indicates
   // if the O was properly positioned.  There are distance sanity
   // checks applied also.
   std::pair<short int, clipper::Coord_orth> o_position(const minimol::residue &res_with_CA_C, const minimol::residue &res_with_N);

}
