/*
 * coords/graphics-line.hh
 * 
 * Copyright 2020 by Medical Research Council
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
 */

#ifndef GRAPHICS_LINE_HH
#define GRAPHICS_LINE_HH

#include "Cartesian.h"

class graphics_line_t {
public:
   enum cylinder_class_t { UNK, SINGLE, DOUBLE, TRIPLE, KEK_DOUBLE_BOND_INNER_BOND }; // so that double bonds can be drawn thinner
                                                          // than single bonds (likewise triple)
   cylinder_class_t cylinder_class;
   coot::CartesianPair positions;
   bool has_begin_cap;
   bool has_end_cap;
   // int residue_index;
   // restore this when finished
   // mmdb::Residue *residue_p; // the residue for the bond (maybe there should be 2 residues_ps? because
                             // sometimes there will be 2 residues for the same graphics_line_t.
                             // Hmm.
   int model_number; // -1 is unset
   int atom_index_1;
   int atom_index_2;
#if 0
   // default single bond constructor
   graphics_line_t(const coot::CartesianPair &p, bool b, bool e) {
      positions = p;
      has_begin_cap = b;
      has_end_cap = e;
      cylinder_class = SINGLE;
      // residue_index = -1; // unset
      //residue_p = 0;
      atom_index_1 = -1;
      atom_index_2 = -1;
      model_number = -1;
   }
#endif
   // we want atom indices now, not just the residue
   // graphics_line_t(const coot::CartesianPair &p, cylinder_class_t cc, bool b, bool e, mmdb::Residue *residue_p_in);

   graphics_line_t(const coot::CartesianPair &p, cylinder_class_t cc, bool b, bool e,
		   int model_no_in,
		   int atom_index_1_in, int atom_index_2_in) : positions(p) {
      has_begin_cap = b;
      has_end_cap = e;
      cylinder_class = cc;
      atom_index_1 = atom_index_1_in;
      atom_index_2 = atom_index_2_in;
      model_number = model_no_in;
   }
   graphics_line_t() { }
};

#endif // GRAPHICS_LINE_HH

