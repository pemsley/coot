/*
 * api/bond-colour.hh
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef API_BOND_COLOUR_HH
#define API_BOND_COLOUR_HH

namespace coot {
   enum class api_bond_colour_t { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
      COLOUR_BY_CHAIN_BONDS=3,
      CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
      BONDS_NO_HYDROGENS=15,
      CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
      CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
      CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
      COLOUR_BY_MOLECULE_BONDS=8,
      COLOUR_BY_RAINBOW_BONDS=9,
      COLOUR_BY_B_FACTOR_BONDS=10,
      COLOUR_BY_OCCUPANCY_BONDS=11,
      COLOUR_BY_USER_DEFINED_COLOURS____BONDS=12,
      COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS=13,
      COLOUR_BY_CHAIN_GOODSELL=21 // 20230826-PE from coords/Bond_lines.h
   };
}

#endif // API_BOND_COLOUR_HH
