/*
 * docking/haddock-utils.hh
 *
 * Copyright 2025 by Medical Research Council
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef HADDOCK_UTILS_HH
#define HADDOCK_UTILS_HH

#include <vector>
#include <mmdb2/mmdb_manager.h>

#include "haddock-types.hh"

namespace coot {
namespace haddock {

   // Filter a list of candidate active residues by solvent accessibility.
   // Returns only those residues with relative SASA > threshold (default 50%).
   //
   std::vector<coot::residue_spec_t>
   filter_by_sasa(mmdb::Manager *mol,
                  const std::vector<coot::residue_spec_t> &candidates,
                  double sasa_threshold = 0.5);

   // Find passive residues: surface neighbours of the active residues
   // that have relative SASA > threshold and are within neighbour_dist
   // of any active residue atom.
   //
   std::vector<coot::residue_spec_t>
   find_passive_residues(mmdb::Manager *mol,
                         const std::vector<coot::residue_spec_t> &active_residues,
                         double sasa_threshold = 0.5,
                         double neighbour_dist = 6.5);

   // Build the complete active_passive_residues_t from a list of
   // candidate active residues.  Filters by SASA, then finds passive
   // residues automatically.
   //
   active_passive_residues_t
   define_residues(mmdb::Manager *mol,
                   const std::vector<coot::residue_spec_t> &active_candidates,
                   double sasa_threshold = 0.5,
                   double neighbour_dist = 6.5);

} // namespace haddock
} // namespace coot

#endif // HADDOCK_UTILS_HH
