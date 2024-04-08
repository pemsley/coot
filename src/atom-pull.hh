/*
 * src/atom-pull.hh
 *
 * Copyright 2015 by Medical Research Council
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


#ifndef ATOM_PULL_HH
#define ATOM_PULL_HH

#include <clipper/core/coords.h>
#include "geometry/residue-and-atom-specs.hh"

class atom_pull_info_t {

   bool status;
public:
   coot::atom_spec_t spec;
   clipper::Coord_orth pos;


   atom_pull_info_t() { status = false; }
   atom_pull_info_t(const coot::atom_spec_t &spec_in, const clipper::Coord_orth &pos_in) :
      spec(spec_in), pos(pos_in) {
      status = true;
   }

   void off() { status = false; }

   void on() { status = true; }

   bool get_status() const { return status; }

   // return first false if number on not found
   std::pair<bool, int> find_spec(mmdb::PAtom *atoms, int n_atoms) const {

      int idx = -1;
      bool local_status = false;
      if (status) { 
	 for (int iat=0; iat<n_atoms; iat++) {
	    if (coot::atom_spec_t(atoms[iat]).is_same(spec)) {
	       idx = iat;
	       local_status = true;
	       break;
	    } 
	 }
      } 
      return std::pair<bool, int> (local_status, idx);
   }
};

#endif // ATOM_PULL_HH

