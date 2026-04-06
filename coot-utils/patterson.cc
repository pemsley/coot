/*
 * coot-utils/patterson.cc
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

#include <thread>
#include <algorithm>
#include <cmath>
#include <iostream>

#include <clipper/clipper-contrib.h>
#include <clipper/core/map_interp.h>

#include "utils/split-indices.hh"
#include "utils/coot-utils.hh"
#include "coot-map-utils.hh"
#include "patterson.hh"

clipper::Xmap<float>
coot::make_patterson_from_map(const clipper::Xmap<float> &xmap,
                              const clipper::Resolution &reso) {

   float mg = coot::util::max_gridding(xmap);
   clipper::Resolution map_reso(2.0 * mg);
   // Use the coarser of the requested and map Nyquist resolution
   clipper::Resolution use_reso = (reso.limit() > map_reso.limit()) ? reso : map_reso;
   clipper::HKL_info hkls(xmap.spacegroup(), xmap.cell(), use_reso, true);
   clipper::HKL_data<clipper::datatypes::F_phi<float>> fphi(hkls);
   xmap.fft_to(fphi);

   // Normalise to E-values in resolution shells, then use E^2 as Patterson coefficients.
   // This removes the fall-off with resolution and the molecular envelope contribution,
   // giving a sharper Patterson with better discrimination in the rotation function.
   typedef clipper::HKL_data_base::HKL_reference_index HRI;

   // Compute <F^2> in resolution shells
   int n_bins = 20;
   std::vector<double> shell_sum_f2(n_bins, 0.0);
   std::vector<int> shell_count(n_bins, 0);
   double inv_reso_sq_max = 1.0 / (use_reso.limit() * use_reso.limit());

   for (HRI ih = fphi.first(); !ih.last(); ih.next()) {
      if (!fphi[ih].missing()) {
         float f = fphi[ih].f();
         double inv_reso_sq = ih.invresolsq();
         int bin = static_cast<int>(inv_reso_sq / inv_reso_sq_max * n_bins);
         if (bin >= n_bins) bin = n_bins - 1;
         if (bin < 0) bin = 0;
         shell_sum_f2[bin] += static_cast<double>(f) * f;
         shell_count[bin]++;
      }
   }

   std::vector<double> shell_mean_f2(n_bins, 1.0);
   for (int i=0; i<n_bins; i++) {
      if (shell_count[i] > 0)
         shell_mean_f2[i] = shell_sum_f2[i] / shell_count[i];
   }

   // Compute E^2 Patterson coefficients
   for (HRI ih = fphi.first(); !ih.last(); ih.next()) {
      if (!fphi[ih].missing()) {
         float f = fphi[ih].f();
         double inv_reso_sq = ih.invresolsq();
         int bin = static_cast<int>(inv_reso_sq / inv_reso_sq_max * n_bins);
         if (bin >= n_bins) bin = n_bins - 1;
         if (bin < 0) bin = 0;
         double e = f / std::sqrt(shell_mean_f2[bin]);
         fphi[ih].f() = static_cast<float>(e * e);
      }
      fphi[ih].phi() = 0.0;
   }

   clipper::Grid_sampling grid(xmap.spacegroup(), xmap.cell(), reso);
   clipper::Xmap<float> patterson(xmap.spacegroup(), xmap.cell(), grid);
   patterson.fft_from(fphi);
   return patterson;
}

// Original interface: computes fphi_obs internally (one-off use)
std::vector<coot::translation_search_result_t>
coot::phased_translation_search(const clipper::Xmap<float> &observed_map,
                                const clipper::Xmap<float> &model_map,
                                unsigned int n_results,
                                clipper::Xmap<float> *tf_map_p) {

   clipper::HKL_info hkls(observed_map.spacegroup(), observed_map.cell(),
                           clipper::Resolution(2.0 * coot::util::max_gridding(observed_map)), true);
   clipper::HKL_data<clipper::datatypes::F_phi<float>> fphi_obs(hkls);
   observed_map.fft_to(fphi_obs);

   return phased_translation_search(fphi_obs, hkls, observed_map, model_map, n_results, tf_map_p);
}

#include <clipper/clipper-ccp4.h>

// Overload that accepts pre-computed observed-map structure factors
std::vector<coot::translation_search_result_t>
coot::phased_translation_search(const clipper::HKL_data<clipper::datatypes::F_phi<float>> &fphi_obs,
                                const clipper::HKL_info &hkls,
                                const clipper::Xmap<float> &observed_map,
                                const clipper::Xmap<float> &model_map,
                                unsigned int n_results,
                                clipper::Xmap<float> *tf_map_p) {

   clipper::HKL_data<clipper::datatypes::F_phi<float>> fphi_model(hkls);
   model_map.fft_to(fphi_model);

   if (true) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write("atom-map.map");
      mapout.export_xmap(model_map);
      mapout.close_write();
   }

   typedef clipper::HKL_data_base::HKL_reference_index HRI;

   // Cross-correlation coefficients: F_obs * conj(F_model)
   clipper::HKL_data<clipper::datatypes::F_phi<float>> fphi_cc(hkls);

   unsigned int n_total = 0;
   unsigned int n_missing_obs = 0;
   unsigned int n_missing_model = 0;
   unsigned int n_filled = 0;
   for (HRI ih = fphi_obs.first(); !ih.last(); ih.next()) {
      n_total += 1;
      if (!fphi_obs[ih].missing() && !fphi_model[ih].missing()) {
         fphi_cc[ih].f() = fphi_obs[ih].f() * fphi_model[ih].f();
         fphi_cc[ih].phi() = fphi_obs[ih].phi() - fphi_model[ih].phi();
         n_filled++;
      } else {
         fphi_cc[ih].f() = 0.0f;
         fphi_cc[ih].phi() = 0.0f;
         if (fphi_obs[ih].missing()) n_missing_obs++;
         if (fphi_model[ih].missing()) n_missing_model++;
      }
   }
   std::cout << "DEBUG:: reflections: n_total " << n_total << " n_missing_obs " << n_missing_obs << " n_missing_model " << n_missing_model << " n_filled " << n_filled << std::endl;

   // IFFT to get the translation function map
   clipper::Xmap<float> tf_map(observed_map.spacegroup(), observed_map.cell(),
                                observed_map.grid_sampling());
   tf_map.fft_from(fphi_cc);

   // Optionally return the TF map
   if (tf_map_p)
      *tf_map_p = tf_map;

   // Find the top peaks by scanning the map
   // Collect all grid points and their values, then partial-sort
   struct grid_peak_t {
      clipper::Coord_grid coord;
      float value;
   };
   std::vector<grid_peak_t> peaks;
   peaks.reserve(1000);

   // First pass: find the global statistics
   double sum = 0.0;
   double sum_sq = 0.0;
   unsigned int count = 0;
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = tf_map.first(); !ix.last(); ix.next()) {
      float v = tf_map[ix];
      sum += v;
      sum_sq += static_cast<double>(v) * v;
      count++;
   }
   double mean = sum / count;
   double variance = sum_sq / count - mean * mean;
   double sd = std::sqrt(std::max(variance, 0.0));

   std::cout << "INFO:: phased_translation_search() TF map: mean=" << mean
             << " sd=" << sd << " n_grid=" << count << std::endl;

   // Second pass: collect points above a threshold (mean + 2*sigma),
   // or just collect all if the map is small
   float threshold = static_cast<float>(mean + 2.0 * sd);
   for (ix = tf_map.first(); !ix.last(); ix.next()) {
      float v = tf_map[ix];
      if (v > threshold) {
         peaks.push_back({ix.coord(), v});
      }
   }

   // Sort by descending value
   std::sort(peaks.begin(), peaks.end(),
             [](const grid_peak_t &a, const grid_peak_t &b) {
                return a.value > b.value;
             });

   // Convert to results with orthogonal coordinates
   unsigned int n_out = std::min(n_results, static_cast<unsigned int>(peaks.size()));
   std::vector<translation_search_result_t> results(n_out);
   for (unsigned int i=0; i<n_out; i++) {
      clipper::Coord_orth pos = peaks[i].coord.coord_frac(tf_map.grid_sampling())
                                                .coord_orth(tf_map.cell());
      float sigma_score = static_cast<float>((peaks[i].value - mean) / sd);
      results[i] = translation_search_result_t(pos, sigma_score);
   }

   return results;
}

std::vector<glm::quat>
coot::make_quaternion_grid(float angular_step_degrees) {

   std::vector<glm::quat> quats;

   // Hopf fibration parameterization of SO(3):
   //   q = (cos(phi)*cos(theta1), sin(phi)*cos(theta2), sin(phi)*sin(theta2), cos(phi)*sin(theta1))
   //   theta1 in [0, 2*pi), theta2 in [0, 2*pi), phi in [0, pi/2]
   //   Uniform measure requires sampling uniformly in cos(2*phi)

   float step_rad = glm::radians(angular_step_degrees);

   // Number of steps for each angular variable.
   // theta1 and theta2 range over 2*pi, phi is parameterized by cos(2*phi) in [-1, 1].
   int n_theta = static_cast<int>(std::ceil(2.0f * M_PI / step_rad));
   int n_phi   = static_cast<int>(std::ceil(2.0f / step_rad));  // range of cos(2*phi) is 2

   if (n_theta < 1) n_theta = 1;
   if (n_phi   < 1) n_phi   = 1;

   quats.reserve(n_theta * n_theta * n_phi);

   for (int i=0; i<n_theta; i++) {
      float theta1 = 2.0f * M_PI * static_cast<float>(i) / static_cast<float>(n_theta);
      for (int j=0; j<n_theta; j++) {
         float theta2 = 2.0f * M_PI * static_cast<float>(j) / static_cast<float>(n_theta);
         for (int k=0; k<n_phi; k++) {
            float cos2phi = 1.0f - (2.0f * k + 1.0f) / static_cast<float>(n_phi);
            float phi = 0.5f * std::acos(cos2phi);
            float sp = std::sin(phi);
            float cp = std::cos(phi);
            glm::quat q(cp * std::cos(theta1),   // w
                        sp * std::cos(theta2),   // x
                        sp * std::sin(theta2),   // y
                        cp * std::sin(theta1));  // z
            quats.push_back(q);
         }
      }
   }

   return quats;
}

clipper::Xmap<float>
coot::rotate_p1_map(const clipper::Xmap<float> &xmap,
                    const clipper::Mat33<> &rotation_matrix) {

   clipper::Xmap<float> rotated(xmap.spacegroup(), xmap.cell(), xmap.grid_sampling());

   clipper::Mat33<> rot_inv = rotation_matrix.inverse();
   clipper::Cell cell = xmap.cell();

   // Cell centre in orthogonal coordinates
   clipper::Coord_frac centre_frac(0.5, 0.5, 0.5);
   clipper::Coord_orth centre = centre_frac.coord_orth(cell);

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = rotated.first(); !ix.last(); ix.next()) {
      // Orthogonal position of this grid point
      clipper::Coord_orth pos = ix.coord().coord_frac(rotated.grid_sampling()).coord_orth(cell);

      // Rotate about cell centre: p' = R_inv * (p - c) + c
      clipper::Vec3<> shifted(pos.x() - centre.x(),
                              pos.y() - centre.y(),
                              pos.z() - centre.z());
      clipper::Vec3<> rotated_shifted = rot_inv * shifted;
      clipper::Coord_orth source_pos(rotated_shifted[0] + centre.x(),
                                     rotated_shifted[1] + centre.y(),
                                     rotated_shifted[2] + centre.z());

      // Check if the source point is inside the cell
      clipper::Coord_frac source_frac = source_pos.coord_frac(cell);
      if (source_frac.u() >= 0.0 && source_frac.u() < 1.0 &&
          source_frac.v() >= 0.0 && source_frac.v() < 1.0 &&
          source_frac.w() >= 0.0 && source_frac.w() < 1.0) {
         float val;
         clipper::Interp_linear::interp(xmap, xmap.coord_map(source_pos), val);
         rotated[ix] = val;
      } else {
         rotated[ix] = 0.0f;
      }
   }

   return rotated;
}

// Worker function for threaded rotation function search.
// Evaluates rotations in the range [range_start, range_end) and writes
// the results into results_for_thread.
static void
rotation_function_search_thread_fn(unsigned int range_start,
                                   unsigned int range_end,
                                   const std::vector<glm::quat> &quats,
                                   const clipper::Xmap<float> &target_patterson,
                                   const clipper::Xmap<float> &search_patterson,
                                   float shell_inner_radius,
                                   float shell_outer_radius,
                                   float grid_step,
                                   std::vector<coot::rotation_function_result_t> *results_p) {

   float r_inner_sq = shell_inner_radius * shell_inner_radius;
   float r_outer_sq = shell_outer_radius * shell_outer_radius;

   // Sample points within the spherical shell on a Cartesian grid
   // centred on the Patterson origin (0,0,0).
   // The rotation function rotates Patterson vectors around the origin.
   struct sample_point_t {
      clipper::Coord_orth pos;
      float target_val;
   };
   std::vector<sample_point_t> sample_points;

   // Sample over the full spherical shell around the Patterson origin.
   // The Patterson is periodic (computed via FFT), so negative coordinates
   // wrap correctly via the Xmap.
   float r_max = shell_outer_radius;
   for (float x=-r_max; x<=r_max; x+=grid_step) {
      for (float y=-r_max; y<=r_max; y+=grid_step) {
         for (float z=-r_max; z<=r_max; z+=grid_step) {
            float r_sq = x*x + y*y + z*z;
            if (r_sq >= r_inner_sq && r_sq <= r_outer_sq) {
               clipper::Coord_orth pos(x, y, z);
               float target_val;
               clipper::Interp_linear::interp(target_patterson, target_patterson.coord_map(pos), target_val);
               sample_points.push_back({pos, target_val});
            }
         }
      }
   }

   // Pre-compute target mean and sum-of-squares for correlation coefficient
   double target_sum = 0.0;
   for (const auto &sp : sample_points)
      target_sum += sp.target_val;
   double target_mean = target_sum / static_cast<double>(sample_points.size());

   double target_var = 0.0;
   for (const auto &sp : sample_points) {
      double d = sp.target_val - target_mean;
      target_var += d * d;
   }
   double target_sd = std::sqrt(target_var);

   for (unsigned int iq=range_start; iq<range_end; iq++) {
      const glm::quat &q = quats[iq];
      glm::mat3 rot_mat = glm::mat3_cast(q);

      // Collect rotated search Patterson values and compute correlation
      double search_sum = 0.0;
      double cross_sum = 0.0;
      double search_sq_sum = 0.0;
      unsigned int n = sample_points.size();

      for (const auto &sp : sample_points) {
         glm::vec3 p(sp.pos.x(), sp.pos.y(), sp.pos.z());
         glm::vec3 rp = rot_mat * p;
         clipper::Coord_orth rotated_pos(rp.x, rp.y, rp.z);

         float search_val;
         clipper::Interp_linear::interp(search_patterson, search_patterson.coord_map(rotated_pos), search_val);
         search_sum += search_val;
         search_sq_sum += static_cast<double>(search_val) * search_val;
         cross_sum += sp.target_val * static_cast<double>(search_val);
      }

      double search_mean = search_sum / static_cast<double>(n);
      double search_var = search_sq_sum - n * search_mean * search_mean;
      double search_sd = std::sqrt(std::max(search_var, 0.0));

      // Pearson correlation coefficient
      double cc = 0.0;
      if (target_sd > 0.0 && search_sd > 0.0) {
         double cov = cross_sum - n * target_mean * search_mean;
         cc = cov / (target_sd * search_sd);
      }
      (*results_p)[iq - range_start] = coot::rotation_function_result_t(q, static_cast<float>(cc));
   }
}

std::vector<coot::rotation_function_result_t>
coot::rotation_function_search(const clipper::Xmap<float> &target_patterson,
                               const clipper::Xmap<float> &search_patterson,
                               float shell_inner_radius,
                               float shell_outer_radius,
                               float angular_step_degrees) {

   std::vector<glm::quat> quats = make_quaternion_grid(angular_step_degrees);
   unsigned int n_quats = quats.size();

   std::cout << "INFO:: rotation_function_search(): " << n_quats << " rotations to evaluate" << std::endl;

   // Grid step for sampling within the spherical shell.
   // Use ~2 Angstroms — finer than the typical Patterson feature width.
   float grid_step = 2.0f;

   unsigned int n_threads = coot::get_max_number_of_threads();
   std::vector<std::pair<unsigned int, unsigned int> > ranges = atom_index_ranges(n_quats, n_threads);

   // Pre-allocate per-thread result vectors
   std::vector<std::vector<rotation_function_result_t>> results_per_thread(ranges.size());
   for (unsigned int i=0; i<ranges.size(); i++)
      results_per_thread[i].resize(ranges[i].second - ranges[i].first);

   std::vector<std::thread> threads;
   for (unsigned int i_thread=0; i_thread<ranges.size(); i_thread++) {
      threads.push_back(std::thread(rotation_function_search_thread_fn,
                                    ranges[i_thread].first,
                                    ranges[i_thread].second,
                                    std::cref(quats),
                                    std::cref(target_patterson),
                                    std::cref(search_patterson),
                                    shell_inner_radius,
                                    shell_outer_radius,
                                    grid_step,
                                    &results_per_thread[i_thread]));
   }

   for (unsigned int i_thread=0; i_thread<threads.size(); i_thread++)
      threads[i_thread].join();

   // Merge results
   std::vector<rotation_function_result_t> results;
   results.reserve(n_quats);
   for (const auto &thread_results : results_per_thread)
      for (const auto &r : thread_results)
         results.push_back(r);

   // Sort by descending score
   std::sort(results.begin(), results.end(),
             [](const rotation_function_result_t &a, const rotation_function_result_t &b) {
                return a.score > b.score;
             });

   return results;
}
