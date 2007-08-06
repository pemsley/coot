/* src/richardson-rotamer.hh
 * 
 * Copyright 2007 by The Unversity of Oxford
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

#ifndef RICHARDSON_ROTAMER_HH
#define RICHARDSON_ROTAMER_HH

#include "rotamer.hh"

namespace coot {

   class richardson_rotamer : public rotamer {
   public:
      richardson_rotamer(CResidue *res) : 
	 rotamer(residue, 0) {}
   };

}


#endif // RICHARDSON_ROTAMER_HH
