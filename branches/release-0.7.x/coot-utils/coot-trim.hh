/* coot-utils/coot-trim.hh
 * 
 * Copyright 2006, The University of York
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

// #include "mmdb-extras.h"

namespace coot {

   namespace util { 

      enum {TRIM_BY_MAP_DELETE, TRIM_BY_MAP_ZERO_OCC};

      // return the number of trimmed atoms
      int
      trim_molecule_by_map(CMMDBManager *mol,
			   const clipper::Xmap<float> &xmap,
			   float map_level,
			   short int remove_or_zero_occ_flag,
			   short int waters_only_flag);

   }

}   
