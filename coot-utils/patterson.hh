/*
 * coot-utils/patterson.hh
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

#ifndef PATTERSON_HH
#define PATTERSON_HH

#include <vector>
#include <clipper/core/xmap.h>
#include <clipper/core/hkl_info.h>
#include <clipper/core/hkl_data.h>
#include <clipper/core/hkl_datatypes.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

namespace coot {

   // Make a Patterson map from an existing map (P1 cryo-EM use case).
   // The resolution controls how many reflections are generated for the FFT.
   clipper::Xmap<float> make_patterson_from_map(const clipper::Xmap<float> &xmap,
                                                const clipper::Resolution &reso);

   // Generate a quasi-uniform grid of rotations on SO(3) using the Hopf fibration.
   // angular_step_degrees controls the sampling density.
   std::vector<glm::quat> make_quaternion_grid(float angular_step_degrees);

   struct rotation_function_result_t {
      glm::quat rotation;
      float score;
      rotation_function_result_t() : rotation(glm::quat(1,0,0,0)), score(0.0f) {}
      rotation_function_result_t(const glm::quat &q, float s) : rotation(q), score(s) {}
   };

   // Rotate a P1 map about the centre of the cell.
   // Points that rotate outside the cell are set to zero (no wrapping).
   clipper::Xmap<float> rotate_p1_map(const clipper::Xmap<float> &xmap,
                                       const clipper::Mat33<> &rotation_matrix);

   struct translation_search_result_t {
      clipper::Coord_orth position;
      float score;
      translation_search_result_t() : position(0,0,0), score(0.0f) {}
      translation_search_result_t(const clipper::Coord_orth &p, float s) : position(p), score(s) {}
   };

   // Phased translation function.
   //
   // Given an observed map and a model map (computed from the oriented model
   // placed at the origin of the same cell), compute the cross-correlation
   // map via FFT:
   //
   //   PTF(t) = IFFT[ F_obs(h) * F_model*(h) ]
   //
   // The peak in this map gives the translation that best overlaps the model
   // with the observed density.
   //
   // observed_map:  the experimental electron density map
   // model_map:     density calculated from the oriented model at the origin,
   //                in the same cell/grid as observed_map
   // n_results:     number of top peaks to return
   // tf_map_p:      if non-null, the translation function map is copied here
   //
   // Returns results sorted by descending score.
   std::vector<translation_search_result_t>
   phased_translation_search(const clipper::Xmap<float> &observed_map,
                             const clipper::Xmap<float> &model_map,
                             unsigned int n_results = 10,
                             clipper::Xmap<float> *tf_map_p = nullptr);

   // Overload that accepts pre-computed observed-map structure factors,
   // avoiding the expensive observed_map.fft_to() on every call.
   //
   // fphi_obs:  pre-computed F/phi from the observed map
   // hkls:      the HKL_info used to compute fphi_obs
   // model_map: density calculated from the oriented model at the origin,
   //            in the same cell/grid as observed_map
   // observed_map: needed only for cell/grid/spacegroup (not FFT'd again)
   // n_results: number of top peaks to return
   // tf_map_p:  if non-null, the translation function map is copied here
   std::vector<translation_search_result_t>
   phased_translation_search(const clipper::HKL_data<clipper::datatypes::F_phi<float>> &fphi_obs,
                             const clipper::HKL_info &hkls,
                             const clipper::Xmap<float> &observed_map,
                             const clipper::Xmap<float> &model_map,
                             unsigned int n_results = 10,
                             clipper::Xmap<float> *tf_map_p = nullptr);

   // Rotation function search: compare two Patterson maps in a spherical shell,
   // scanning over SO(3). Returns results sorted by descending score.
   //
   // target_patterson: Patterson of the cryo-EM map
   // search_patterson: Patterson of the search model map
   // shell_inner_radius: inner radius of spherical shell (excludes origin peak)
   // shell_outer_radius: outer radius of spherical shell
   // angular_step_degrees: rotation sampling step
   std::vector<rotation_function_result_t>
   rotation_function_search(const clipper::Xmap<float> &target_patterson,
                            const clipper::Xmap<float> &search_patterson,
                            float shell_inner_radius,
                            float shell_outer_radius,
                            float angular_step_degrees);

} // namespace coot

#endif // PATTERSON_HH
