/*
 * coot-utils/zernike.hh
 *
 * Copyright 2025 by Medical Research Council
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

#ifndef COOT_ZERNIKE_HH
#define COOT_ZERNIKE_HH

#include <vector>
#include <complex>
#include <clipper/core/xmap.h>
#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   // 3D Zernike moments for rotation-invariant map comparison
   //
   // Use case: Place AlphaFold domain models into experimental cryo-EM maps.
   // Zernike moments are rotation-invariant, so only translation needs to be searched.
   //
   // The rotation-invariant descriptor F_nl = sqrt(sum_m |Omega_nl^m|^2)
   // allows matching a probe to a target regardless of orientation.

   class ZernikeDescriptor {
   public:
      // Compute Zernike descriptor from an Xmap around a centre point
      //
      // @param xmap The electron density map
      // @param centre The centre point for the spherical sample region
      // @param radius The radius of the sampling sphere (Angstroms)
      // @param n_grid Number of grid points along each dimension for sampling
      // @param order Maximum Zernike order n (typically 15-20)
      ZernikeDescriptor(const clipper::Xmap<float> &xmap,
                        const clipper::Coord_orth &centre,
                        float radius,
                        int n_grid = 64,
                        int order = 20);

      // Compute Zernike descriptor with protein-shaped mask
      // Only samples points within mask_radius of CA atoms from the model
      //
      // @param xmap The electron density map
      // @param centre The centre point for the spherical sample region
      // @param radius The radius of the sampling sphere (Angstroms)
      // @param mol The model used to define the protein mask (CA atoms)
      // @param mask_radius Distance threshold from CA atoms (Angstroms)
      // @param n_grid Number of grid points along each dimension for sampling
      // @param order Maximum Zernike order n (typically 15-20)
      ZernikeDescriptor(const clipper::Xmap<float> &xmap,
                        const clipper::Coord_orth &centre,
                        float radius,
                        mmdb::Manager *mol,
                        float mask_radius,
                        int n_grid = 64,
                        int order = 20);

      // Compute Zernike descriptor with pre-computed relative mask
      // Uses relative CA positions for consistent masking at any centre position
      //
      // @param xmap The electron density map
      // @param centre The centre point for the spherical sample region
      // @param radius The radius of the sampling sphere (Angstroms)
      // @param relative_ca_positions CA positions relative to reference centre
      // @param mask_radius Distance threshold from CA atoms (Angstroms)
      // @param n_grid Number of grid points along each dimension for sampling
      // @param order Maximum Zernike order n (typically 15-20)
      ZernikeDescriptor(const clipper::Xmap<float> &xmap,
                        const clipper::Coord_orth &centre,
                        float radius,
                        const std::vector<clipper::Coord_orth> &relative_ca_positions,
                        float mask_radius,
                        int n_grid = 64,
                        int order = 20);

      // Default constructor for empty descriptor
      ZernikeDescriptor() : order_(0), radius_(0.0f), n_grid_(0) {}

      // Get the rotation-invariant descriptor (vector of F_nl norms)
      // The vector contains F_nl for all valid (n,l) pairs where:
      //   - n ranges from 0 to order
      //   - l = n, n-2, n-4, ... down to 0 or 1
      std::vector<float> get_invariants() const;

      // Get the number of invariant values
      // For order N, this is sum over n=0..N of floor((n+2)/2)
      size_t n_invariants() const;

      // Compare two descriptors using Euclidean distance
      // Returns the L2 distance between the invariant vectors
      static float distance(const ZernikeDescriptor &d1, const ZernikeDescriptor &d2);

      // Compare two descriptors using normalized correlation
      // Returns value in [-1, 1], where 1 means identical
      static float correlation(const ZernikeDescriptor &d1, const ZernikeDescriptor &d2);

      // Compare two descriptors using cosine similarity
      // Returns value in [-1, 1], where 1 means identical direction
      // Scale-invariant: doesn't depend on magnitude of density values
      static float cosine_similarity(const ZernikeDescriptor &d1, const ZernikeDescriptor &d2);

      // Get order
      int order() const { return order_; }

      // Get radius
      float radius() const { return radius_; }

      // Check if descriptor is valid (has been computed)
      bool is_valid() const { return !invariants_.empty(); }

      // Diagnostic information
      struct Diagnostics {
         size_t n_grid_points;      // Points sampled (inside sphere AND mask if used)
         size_t n_points_in_sphere; // Points inside sphere (before mask)
         size_t n_points_in_mask;   // Points that passed mask test
         size_t n_points_outside_mask; // Points rejected by mask
         float density_min;
         float density_max;
         float density_mean;
         float density_sum;
      };

      // Get diagnostic information from the computation
      Diagnostics get_diagnostics() const { return diagnostics_; }

      // Print diagnostic information
      void print_diagnostics() const;

      // Print first few raw moments (complex values)
      void print_moments(int max_n = 3) const;

   private:
      int order_;
      float radius_;
      int n_grid_;
      std::vector<float> invariants_;  // F_nl values
      Diagnostics diagnostics_;

      // Full complex moments Omega_nl^m stored for potential alignment recovery
      // Index: moments_[index(n,l,m)] where index maps (n,l,m) to linear array
      std::vector<std::complex<double>> moments_;

      // Compute moments from grid sample values
      // grid_values: vector of (x, y, z, density) tuples in normalized coords [-1,1]
      void compute_moments(const std::vector<std::tuple<float,float,float,float>> &grid_values);

      // Normalize moments by density sum (shape, not mass)
      void normalize_moments();

      // Convert moments to rotation-invariant descriptors
      void compute_invariants();

      // Map (n,l,m) to linear index in moments_ array
      size_t moment_index(int n, int l, int m) const;

      // Compute 3D Zernike basis function Z_nl^m at a point
      // Point (r, theta, phi) in spherical coordinates, r in [0,1]
      std::complex<double> zernike_basis(int n, int l, int m,
                                          double r, double theta, double phi) const;

      // Zernike radial polynomial R_nl(r)
      double radial_polynomial(int n, int l, double r) const;

      // Spherical harmonic Y_l^m(theta, phi)
      std::complex<double> spherical_harmonic(int l, int m, double theta, double phi) const;

      // Associated Legendre polynomial P_l^m(x) (unnormalized)
      double assoc_legendre(int l, int m, double x) const;

      // Factorial and double factorial helpers
      static double factorial(int n);
      static double double_factorial(int n);

      // Precomputed Q_nlk coefficients for radial polynomials
      // Q_nlk = (-1)^k * sqrt(2n+3) * binomial(k, n-l)/2) * ... (Novotni formula)
      mutable std::vector<std::vector<std::vector<double>>> Q_coeffs_;
      void precompute_Q_coefficients() const;
      double Q_nlk(int n, int l, int k) const;
   };

   // Helper: count (n,l) pairs for given max order
   // Returns the dimension of the invariant vector
   inline size_t zernike_invariant_count(int order) {
      size_t count = 0;
      for (int n = 0; n <= order; ++n) {
         // l = n, n-2, n-4, ... down to 0 or 1
         // number of valid l values = floor((n+2)/2)
         count += (n + 2) / 2;
      }
      return count;
   }

} // namespace coot

#endif // COOT_ZERNIKE_HH
