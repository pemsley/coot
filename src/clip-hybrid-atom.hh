/*
 * src/clip-hybrid-atom.hh
 *
 * Copyright 2022 by Medical Research Council
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
#ifndef CLIP_HYBRID_ATOM_HH
#define CLIP_HYBRID_ATOM_HH

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-crystal.hh"

namespace coot {
   class clip_hybrid_atom {
   public:
      mmdb::Atom *atom;
      // clipper::Coord_orth pos;
      coot::Cartesian pos;
      clip_hybrid_atom() { atom = NULL; }
      clip_hybrid_atom(mmdb::Atom *mmdb_atom_p, const coot::Cartesian &p) : atom(mmdb_atom_p), pos(p) {}
   };
}

#endif // CLIP_HYBRID_ATOM_HH
