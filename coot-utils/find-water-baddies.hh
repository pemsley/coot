/*
 * coot-utils/find-water-baddies.hh
 *
 * Copyright 2023 by Medical Research Council
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
#ifndef FIND_WATER_BADDIES_HH
#define FIND_WATER_BADDIES_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "atom-selection-container.hh"
#include <clipper/core/xmap.h>

namespace coot {
   std::vector <atom_spec_t>
   find_water_baddies_OR(atom_selection_container_t asc,
                         float b_factor_lim, const clipper::Xmap<float> &xmap_in,
                         float map_in_sigma,
                         float outlier_sigma_level,
                         float min_dist, float max_dist,
                         short int ignore_part_occ_contact_flag,
                         short int ignore_zero_occ_flag);

}


#endif // FIND_WATER_BADDIES_HH
