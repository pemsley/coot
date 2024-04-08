/*
 * coot-utils/plane-utils.hh
 *
 * Copyright 2018 by Medical Research Council
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

#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

namespace coot {

   // return a true in the first if the second is a valid angle (in degress)
   //
   // vector will be typically the difference in positoin between a ring atom and
   // an atom attached to that ring atom.  This function is used to distinugish
   // alpha and beta linking in carbohydrates
   //
   std::pair<bool, double> angle_betwen_plane_and_vector(mmdb::Residue *residue_p,
							 const std::vector<std::string> &ring_atom_names,
							 const std::string &altconf,
							 const clipper::Coord_orth &vector);

   std::pair<bool, double> angle_betwen_plane_and_vector(mmdb::Residue *residue_p,
							 mmdb::Atom *atom_in_ring,
							 mmdb::Atom *bonding_atom);

}
