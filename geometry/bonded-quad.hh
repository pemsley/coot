/*
 * geometry/bonded-quad.hh
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
 *
 */
//
#ifndef BONDED_QUAD_HH
#define BONDED_QUAD_HH

#include "mini-mol/atom-quads.hh"

// Bond between atoms 2 and 3, with 1 and 4 used to find the positions of the bond
//
class bonded_quad_atoms : public coot::atom_quad {
public:
   enum bond_t { NONE, SINGLE, DOUBLE, TRIPLE, DELOC };
   bonded_quad_atoms() { bond_type = NONE; }
   bond_t bond_type;
};

class bonded_quad_atom_names : public coot::atom_name_quad {
public:
   enum bond_t { UNASSIGNED, SINGLE, DOUBLE, TRIPLE, DELOC };
   bonded_quad_atom_names() {}
   bonded_quad_atom_names(const std::string &atom_name_0,
                          const std::string &atom_name_1,
                          const std::string &atom_name_2,
                          const std::string &atom_name_3) : atom_name_quad(atom_name_0,
                                                                           atom_name_1,
                                                                           atom_name_2,
                                                                           atom_name_3) {
      bond_type = UNASSIGNED;
   }
   bond_t bond_type;
};


#endif // BONDED_QUAD_HH

//
