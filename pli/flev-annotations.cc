/*
 * pli/flev-annotations.cc
 *
 * Copyright 2017 by Medical Research Council
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

#include <string>
#include <vector>

#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"
#include "lidia-core/lbg-shared.hh"
#include "flev-annotations.hh"

//
std::ostream& coot::operator<<(std::ostream &s, coot::fle_ligand_bond_t flb) {

   s << "Ligand-H-bond: " << flb.bond_type << " lig-at: " << flb.ligand_atom_spec
     << " " << flb.interacting_residue_atom_spec << " length: " << flb.bond_length;
   if (flb.is_H_bond_to_water)
      s << " (water)";
   return s;
}
