/*
 * api/plain-atom-overlap.nn
 *
 * Copyright 2024 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#ifndef PLAIN_ATOM_OVERLAP_HH
#define PLAIN_ATOM_OVERLAP_HH

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class plain_atom_overlap_t {
   public:
      int ligand_atom_index;
      atom_spec_t atom_spec_1;
      atom_spec_t atom_spec_2;
      double overlap_volume;
      double r_1;
      double r_2;
      bool is_h_bond;
      plain_atom_overlap_t() : ligand_atom_index(-1), overlap_volume(-1), r_1(-1), r_2(-1), is_h_bond(false) {}
      plain_atom_overlap_t(int ligand_atom_index_in,
                           const atom_spec_t &as_1,
                           const atom_spec_t &as_2,
                           double ov,
                           double r1,
                           double r2,
                           bool isHb) : ligand_atom_index(ligand_atom_index_in),
                                        atom_spec_1(as_1),
                                        atom_spec_2(as_2),
                                        overlap_volume(ov),
                                        r_1(r1),
                                        r_2(r2),
                                        is_h_bond(isHb) {}
   };

}

#endif // PLAIN_ATOM_OVERLAP_HH
