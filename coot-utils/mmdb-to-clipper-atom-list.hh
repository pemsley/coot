/*
 * coot-utils/mmdb-to-clipper-atom-list.hh
 *
 * Copyright 2026 Jordan Dialpuri, Medical Research Council Laboratory of Molecular Biology
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth
 * Floor, Boston, MA 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */

#ifndef COOT_UTILS_MMDB_TO_CLIPPER_ATOM_LIST_HH
#define COOT_UTILS_MMDB_TO_CLIPPER_ATOM_LIST_HH

#include <clipper/core/coords.h>
#include <clipper/core/clipper_util.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   //! A drop-in replacement for clipper::MMDBAtom_list.
   //
   // clipper::MMDBAtom (clipper/mmdb/clipper_mmdb.cpp) derives from mmdb::CAtom and
   // reads the atom's coordinate/element/occupancy/B as *public data members*
   // (`return Coord_orth(x,y,z);`, `String(mmdb::CAtom::element)`, ...). Those method
   // bodies are compiled into libclipper-mmdb against real MMDB's memory layout, so
   // reinterpret-casting an mmdb-shim atom (which is gemmi-backed, with a completely
   // different layout, and exposes x()/element via accessor methods) to MMDBAtom*
   // reads garbage — every structure-factor calculation then sees junk coordinates.
   //
   // clipper::Atom / Atom_list are pure clipper value types with no MMDB dependency,
   // so we build the Atom_list ourselves through the shim's accessor methods. The
   // result feeds clipper's SFcalc unchanged. Matches clipper::MMDBAtom's field
   // semantics: element as the PDB-aligned 2-char string, isotropic U from B, and
   // anisotropic U only when ANISOU was present.
   class MMDBAtom_list : public clipper::Atom_list {
   public:
      MMDBAtom_list(const mmdb::PPAtom ppcatom, const int natom) {
         reserve(natom);
         for (int i = 0; i < natom; i++) {
            mmdb::Atom *at = ppcatom[i];
            if (!at) continue;
            clipper::Atom a;
            a.set_element(clipper::String(at->GetElementName()));
            if (!at->isTer())
               a.set_coord_orth(clipper::Coord_orth(at->x(), at->y(), at->z()));
            else
               a.set_coord_orth(clipper::Coord_orth(clipper::Coord_orth::null()));
            a.set_occupancy(at->occupancy());
            a.set_u_iso(clipper::Util::b2u(at->tempFactor()));
            if (at->WhatIsSet & ::mmdb::ASET_Anis_tFac)
               a.set_u_aniso_orth(clipper::U_aniso_orth(at->u11(), at->u22(), at->u33(),
                                                        at->u12(), at->u13(), at->u23()));
            else
               a.set_u_aniso_orth(clipper::U_aniso_orth(clipper::U_aniso_orth::null()));
            push_back(a);
         }
      }
   };

}  // namespace coot

#endif  // COOT_UTILS_MMDB_TO_CLIPPER_ATOM_LIST_HH
