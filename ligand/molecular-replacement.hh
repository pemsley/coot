/*
 * ligand/molecular-replacement.hh
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

#ifndef MOLECULAR_REPLACEMENT_HH
#define MOLECULAR_REPLACEMENT_HH

#include <vector>
#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include <mmdb2/mmdb_manager.h>
#include <glm/gtc/quaternion.hpp>

namespace coot {

   /// A single molecular replacement solution from the core search.
   /// The caller takes ownership of placed_mol and must delete it when done.
   struct mr_solution_t {
      glm::quat rotation;
      float rotation_score;
      clipper::Coord_orth translation;
      float translation_score;
      mmdb::Manager *placed_mol;   ///< placed model (caller owns); nullptr if not created
      mr_solution_t() : rotation(glm::quat(1,0,0,0)), rotation_score(0.0f),
                        translation(0,0,0), translation_score(0.0f),
                        placed_mol(nullptr) {}
   };

   /// Core molecular replacement pipeline.
   ///
   /// Given an observed map, a search model, and an approximate centre,
   /// find the orientation and position that best fit the model into the
   /// density. Returns solutions sorted by translation score (descending).
   ///
   /// The pipeline:
   ///  1. Iterative density-weighted recentring (cf. helix_placement)
   ///  2. Crowther fast rotation function + 2-pass refinement
   ///  3. Phased translation search for each rotation, filtered by distance
   ///
   /// @param xmap_obs the observed (e.g. cryo-EM) map
   /// @param mol_model the search model (not modified; deep-copied internally)
   /// @param target_centre approximate centre to search around
   /// @param n_rotation_solutions number of rotation peaks to refine
   /// @param n_translation_solutions number of translation peaks per rotation
   std::vector<mr_solution_t>
   molecular_replacement_search(const clipper::Xmap<float> &xmap_obs,
                                mmdb::Manager *mol_model,
                                const clipper::Coord_orth &target_centre,
                                int n_rotation_solutions = 10,
                                int n_translation_solutions = 10);

} // namespace coot

#endif // MOLECULAR_REPLACEMENT_HH
