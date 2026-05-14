/*
 * coot-utils/edcalc.hh
 *
 * Copyright 2026 by Medical Research Council
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

#ifndef COOT_EDCALC_HH
#define COOT_EDCALC_HH

#include <clipper/core/xmap.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   // Isotropic electron density calculation — Coot's own version of
   // clipper::EDcalc_iso, with the same scattering factor tables and
   // Gaussian summation. This exists so that we can later parallelise
   // the calculation using a block-colouring scheme.
   //
   // The map is zeroed, then for each atom a 5-Gaussian (+constant)
   // density is accumulated within a sphere of the given radius.
   // Finally the crystallographic multiplicity correction is applied.

   clipper::Xmap<float> calc_atom_map_edcalc(mmdb::Manager *mol,
                                              int atom_selection_handle,
                                              const clipper::Cell &cell,
                                              const clipper::Spacegroup &space_group,
                                              const clipper::Grid_sampling &sampling,
                                              float radius = 2.5f);

} // namespace coot

#endif // COOT_EDCALC_HH
