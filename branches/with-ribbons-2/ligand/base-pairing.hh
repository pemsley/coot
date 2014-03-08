/* ligand/base-pairing.hh
 * 
 * Copyright 2006 by The University of York
 * Copyright 2009 by The University of Oxford
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
 * 02110-1301, USA.
 */

#include <mmdb/mmdb_manager.h>
#include "clipper/core/coords.h"
#include "coot-utils/coot-lsq-types.h"

namespace coot {
   // we need standard residues to make the ideal RNA
   CResidue * watson_crick_partner(CResidue *res_ref, CMMDBManager *standard_residues);
   std::pair<bool, clipper::RTop_orth> 
      base_pair_match_matix(CResidue *res_ref, CResidue *res_moving);

}

