/* ideal/pepflip.hh
 * 
 * Copyright 2002, 2003 The University of York
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


#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <mmdb/mmdb_manager.h>
#include "clipper/core/coords.h"

namespace coot {

   // mmdb-style interface.
   //
   // Given a mol and a residue number and chain of the first residue
   // in the peptide (i.e. the residue with the C and O atoms), flip
   // the C and O atoms of this peptide and the N of the next one (if
   // exists) round a line joining this Ca to the next one.
   //
   // Return status is 0 if the flip did not happen (because, for
   // example, either or both of the Ca's could not be found).
   //
   // Typically, one would copy one's mol (and save it) before calling
   // this.
   // 
   int pepflip(CMMDBManager *mol,
	       const std::string &chain_id,
	       int resno, 
	       const std::string &inscode,
	       const std::string &altconf);

   // You are advised against using this externally.
   // 
   std::vector<clipper::Coord_orth> 
   flip_internal(const std::vector<clipper::Coord_orth> &ca,
		 const std::vector<CAtom *> &atoms); 
   

}
