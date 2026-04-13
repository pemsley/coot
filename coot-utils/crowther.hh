/*
 * coot-utils/crowther.hh
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

#ifndef CROWTHER_HH
#define CROWTHER_HH

#include <vector>
#include <complex>
#include <clipper/core/xmap.h>
#include <clipper/core/hkl_info.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "patterson.hh"

namespace coot {

   // Crowther fast rotation function.
   //
   // Expands two Pattersons in spherical harmonics and evaluates the rotation
   // function R(alpha, beta, gamma) as a Fourier series in alpha and gamma
   // for each sampled value of beta.
   //
   // Reference: Crowther (1972), Navaza (1987), Foadi (2011).
   //
   class crowther_t {
   public:
      // ell_max: maximum spherical harmonic order (controls angular resolution)
      // radius:  integration radius in Angstroms
      // reso:    resolution limit for structure factors
      crowther_t(int ell_max, float radius, const clipper::Resolution &reso);

      // Set target and search maps (fragments centred at origin, P1)
      void set_target(const clipper::Xmap<float> &xmap);
      void set_search(const clipper::Xmap<float> &xmap);

      // Compute the rotation function.
      // n_beta: number of beta samples in (0, pi).
      // Returns results sorted by descending score.
      std::vector<rotation_function_result_t> compute(int n_beta);

      // Evaluate the rotation function at a single set of ZYZ Euler angles (radians).
      // Must call compute() or precompute_c_coefficients() first.
      double evaluate_at_euler(double alpha, double beta, double gamma) const;

      // Refine the top n_refine orientations from a coarse search using a
      // 3x3x3 grid search in ZYZ Euler angles with the given angular step (degrees).
      // Returns the refined results sorted by descending score.
      std::vector<rotation_function_result_t>
      refine_orientations(const std::vector<rotation_function_result_t> &coarse,
                          int n_refine, float step_degrees) const;

      // Print diagnostic information about the expansion coefficients.
      // Useful for validating each step of the computation.
      void print_diagnostics() const;

   private:
      int ell_max_;
      float radius_;
      clipper::Resolution reso_;
      int n_radial_;

      // Per-reflection data extracted from map FFT
      struct reflection_data_t {
         double e_sq;    // normalised |E|^2
         double s;       // reciprocal space distance |s|
         double theta;   // polar angle of s vector
         double phi;     // azimuthal angle of s vector
      };

      std::vector<reflection_data_t> target_refs_;
      std::vector<reflection_data_t> search_refs_;

      // Stored c_{l,m,m'} coefficients from radial integration.
      // c_coeffs_[l] is a (2l+1) x (2l+1) matrix indexed by [m+l][m'+l].
      std::vector<std::vector<std::vector<std::complex<double>>>> c_coeffs_;
      bool c_coeffs_valid_;

      // Compute and store c_{l,m,m'} coefficients.
      void precompute_c_coefficients();

      // Extract reflection data from an xmap, with E-value normalisation
      std::vector<reflection_data_t> extract_reflections(const clipper::Xmap<float> &xmap) const;

      // Precompute Y_lm(theta_s, phi_s) for all reflections and all (l,m).
      // Indexed as [ref_index][lm_index] where lm_index = l*(l+1) + m.
      std::vector<std::vector<std::complex<double>>>
      precompute_ylm(const std::vector<reflection_data_t> &refs) const;

      // Compute expansion coefficients a'_lm(r) at a single radial value.
      // (Omits the constant 4*pi*i^l/V factor which cancels in the product.)
      // Returns vector indexed by lm_index = l*(l+1) + m.
      std::vector<std::complex<double>>
      compute_alm_at_r(const std::vector<reflection_data_t> &refs,
                       const std::vector<std::vector<std::complex<double>>> &ylm_cache,
                       double r) const;

      // Mathematical helper functions
      static double spherical_bessel_j(int l, double x);
      static std::complex<double> spherical_harmonic(int l, int m, double theta, double phi);
      static double associated_legendre(int l, int m, double x);
      static double wigner_d(int l, int mp, int m, double beta);
      static double log_factorial(int n);
   };

} // namespace coot

#endif // CROWTHER_HH
