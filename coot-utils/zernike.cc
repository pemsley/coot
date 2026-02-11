/*
 * coot-utils/zernike.cc
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

#include <cmath>
#include <iostream>
#include <numeric>
#include <limits>
#include <iomanip>
#include "zernike.hh"
#include "coot-map-utils.hh"

namespace coot {

// Constructor: compute Zernike descriptor from map
ZernikeDescriptor::ZernikeDescriptor(const clipper::Xmap<float> &xmap,
                                     const clipper::Coord_orth &centre,
                                     float radius,
                                     int n_grid,
                                     int order)
   : order_(order), radius_(radius), n_grid_(n_grid) {

   // Precompute Q coefficients for radial polynomials
   precompute_Q_coefficients();

   // Allocate moments array
   // Total number of (n,l,m) combinations
   size_t n_moments = 0;
   for (int n = 0; n <= order_; ++n) {
      for (int l = n; l >= 0; l -= 2) {
         n_moments += 2 * l + 1;  // m from -l to +l
      }
   }
   moments_.resize(n_moments, std::complex<double>(0.0, 0.0));

   // Sample the map on a grid within the sphere
   std::vector<std::tuple<float,float,float,float>> grid_values;
   grid_values.reserve(n_grid * n_grid * n_grid);

   float grid_step = 2.0f * radius / (n_grid - 1);

   // Initialize diagnostics (no mask for this constructor)
   diagnostics_.n_grid_points = 0;
   diagnostics_.n_points_in_sphere = 0;
   diagnostics_.n_points_in_mask = 0;
   diagnostics_.n_points_outside_mask = 0;
   diagnostics_.density_min = std::numeric_limits<float>::max();
   diagnostics_.density_max = std::numeric_limits<float>::lowest();
   diagnostics_.density_sum = 0.0f;

   for (int ix = 0; ix < n_grid; ++ix) {
      float x = -radius + ix * grid_step;
      for (int iy = 0; iy < n_grid; ++iy) {
         float y = -radius + iy * grid_step;
         for (int iz = 0; iz < n_grid; ++iz) {
            float z = -radius + iz * grid_step;

            // Check if point is inside the unit ball (normalized coordinates)
            float r_norm = std::sqrt(x*x + y*y + z*z) / radius;
            if (r_norm > 1.0f) continue;

            diagnostics_.n_points_in_sphere++;

            // Get density at this point
            clipper::Coord_orth pt(centre.x() + x, centre.y() + y, centre.z() + z);
            float density = util::density_at_point(xmap, pt);

            // Update diagnostics
            diagnostics_.n_grid_points++;
            diagnostics_.density_sum += density;
            if (density < diagnostics_.density_min) diagnostics_.density_min = density;
            if (density > diagnostics_.density_max) diagnostics_.density_max = density;

            // Store normalized coordinates and density
            grid_values.push_back(std::make_tuple(x/radius, y/radius, z/radius, density));
         }
      }
   }

   // Finalize diagnostics
   diagnostics_.density_mean = (diagnostics_.n_grid_points > 0) ?
      diagnostics_.density_sum / diagnostics_.n_grid_points : 0.0f;

   // Compute moments from grid values
   compute_moments(grid_values);

   // Normalize moments by density sum (shape, not mass)
   normalize_moments();

   // Convert moments to rotation-invariant descriptors
   compute_invariants();
}

// Constructor with protein-shaped mask (CA atoms)
ZernikeDescriptor::ZernikeDescriptor(const clipper::Xmap<float> &xmap,
                                     const clipper::Coord_orth &centre,
                                     float radius,
                                     mmdb::Manager *mol,
                                     float mask_radius,
                                     int n_grid,
                                     int order)
   : order_(order), radius_(radius), n_grid_(n_grid) {

   // Precompute Q coefficients for radial polynomials
   precompute_Q_coefficients();

   // Allocate moments array
   size_t n_moments = 0;
   for (int n = 0; n <= order_; ++n) {
      for (int l = n; l >= 0; l -= 2) {
         n_moments += 2 * l + 1;
      }
   }
   moments_.resize(n_moments, std::complex<double>(0.0, 0.0));

   // Select CA atoms from the model
   int sel_hnd = mol->NewSelection();
   mol->SelectAtoms(sel_hnd, 0, "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*",  // residue name
                    " CA ",  // atom name (with spaces for PDB format)
                    "*",  // element
                    "*"); // alt loc

   mmdb::PPAtom ca_atoms = nullptr;
   int n_ca_atoms = 0;
   mol->GetSelIndex(sel_hnd, ca_atoms, n_ca_atoms);

   // Build vector of CA positions for distance checking
   std::vector<clipper::Coord_orth> ca_positions;
   ca_positions.reserve(n_ca_atoms);
   double ca_min_x = 1e9, ca_max_x = -1e9;
   double ca_min_y = 1e9, ca_max_y = -1e9;
   double ca_min_z = 1e9, ca_max_z = -1e9;
   for (int i = 0; i < n_ca_atoms; ++i) {
      double x = ca_atoms[i]->x;
      double y = ca_atoms[i]->y;
      double z = ca_atoms[i]->z;
      ca_positions.push_back(clipper::Coord_orth(x, y, z));
      if (x < ca_min_x) ca_min_x = x;
      if (x > ca_max_x) ca_max_x = x;
      if (y < ca_min_y) ca_min_y = y;
      if (y > ca_max_y) ca_max_y = y;
      if (z < ca_min_z) ca_min_z = z;
      if (z > ca_max_z) ca_max_z = z;
   }

   std::cout << "CA atoms: " << n_ca_atoms << std::endl;
   std::cout << "CA range X: " << ca_min_x << " to " << ca_max_x << std::endl;
   std::cout << "CA range Y: " << ca_min_y << " to " << ca_max_y << std::endl;
   std::cout << "CA range Z: " << ca_min_z << " to " << ca_max_z << std::endl;

   mol->DeleteSelection(sel_hnd);

   // Lambda to check if a point is within mask_radius of any CA atom
   auto is_near_ca = [&ca_positions, mask_radius](const clipper::Coord_orth &pt) -> bool {
      float mask_radius_sq = mask_radius * mask_radius;
      for (const auto &ca : ca_positions) {
         float dx = pt.x() - ca.x();
         float dy = pt.y() - ca.y();
         float dz = pt.z() - ca.z();
         if (dx*dx + dy*dy + dz*dz < mask_radius_sq) {
            return true;
         }
      }
      return false;
   };

   // Sample the map on a grid within the sphere, masked by CA proximity
   std::vector<std::tuple<float,float,float,float>> grid_values;
   grid_values.reserve(n_grid * n_grid * n_grid);

   float grid_step = 2.0f * radius / (n_grid - 1);

   // Initialize diagnostics
   diagnostics_.n_grid_points = 0;
   diagnostics_.n_points_in_sphere = 0;
   diagnostics_.n_points_in_mask = 0;
   diagnostics_.n_points_outside_mask = 0;
   diagnostics_.density_min = std::numeric_limits<float>::max();
   diagnostics_.density_max = std::numeric_limits<float>::lowest();
   diagnostics_.density_sum = 0.0f;

   for (int ix = 0; ix < n_grid; ++ix) {
      float x = -radius + ix * grid_step;
      for (int iy = 0; iy < n_grid; ++iy) {
         float y = -radius + iy * grid_step;
         for (int iz = 0; iz < n_grid; ++iz) {
            float z = -radius + iz * grid_step;

            // Check if point is inside the unit ball
            float r_norm = std::sqrt(x*x + y*y + z*z) / radius;
            if (r_norm > 1.0f) continue;

            diagnostics_.n_points_in_sphere++;

            // Get world coordinates
            clipper::Coord_orth pt(centre.x() + x, centre.y() + y, centre.z() + z);

            // Check if point is within mask_radius of any CA atom
            if (!is_near_ca(pt)) {
               diagnostics_.n_points_outside_mask++;
               continue;
            }

            diagnostics_.n_points_in_mask++;

            // Get density at this point
            float density = util::density_at_point(xmap, pt);

            // Update diagnostics
            diagnostics_.n_grid_points++;
            diagnostics_.density_sum += density;
            if (density < diagnostics_.density_min) diagnostics_.density_min = density;
            if (density > diagnostics_.density_max) diagnostics_.density_max = density;

            // Store normalized coordinates and density
            grid_values.push_back(std::make_tuple(x/radius, y/radius, z/radius, density));
         }
      }
   }

   // Finalize diagnostics
   diagnostics_.density_mean = (diagnostics_.n_grid_points > 0) ?
      diagnostics_.density_sum / diagnostics_.n_grid_points : 0.0f;

   // Compute moments from grid values
   compute_moments(grid_values);

   // Normalize moments by density sum (shape, not mass)
   normalize_moments();

   // Convert moments to rotation-invariant descriptors
   compute_invariants();
}

// Constructor with pre-computed relative CA positions for consistent masking
ZernikeDescriptor::ZernikeDescriptor(const clipper::Xmap<float> &xmap,
                                     const clipper::Coord_orth &centre,
                                     float radius,
                                     const std::vector<clipper::Coord_orth> &relative_ca_positions,
                                     float mask_radius,
                                     int n_grid,
                                     int order)
   : order_(order), radius_(radius), n_grid_(n_grid) {

   // Precompute Q coefficients for radial polynomials
   precompute_Q_coefficients();

   // Allocate moments array
   size_t n_moments = 0;
   for (int n = 0; n <= order_; ++n) {
      for (int l = n; l >= 0; l -= 2) {
         n_moments += 2 * l + 1;
      }
   }
   moments_.resize(n_moments, std::complex<double>(0.0, 0.0));

   float mask_radius_sq = mask_radius * mask_radius;

   // Lambda to check if a point (relative to centre) is near any relative CA position
   auto is_near_ca = [&relative_ca_positions, mask_radius_sq](float x, float y, float z) -> bool {
      for (const auto &ca_rel : relative_ca_positions) {
         float dx = x - ca_rel.x();
         float dy = y - ca_rel.y();
         float dz = z - ca_rel.z();
         if (dx*dx + dy*dy + dz*dz < mask_radius_sq) {
            return true;
         }
      }
      return false;
   };

   // Sample the map on a grid within the sphere, masked by relative CA proximity
   std::vector<std::tuple<float,float,float,float>> grid_values;
   grid_values.reserve(n_grid * n_grid * n_grid);

   float grid_step = 2.0f * radius / (n_grid - 1);

   // Initialize diagnostics
   diagnostics_.n_grid_points = 0;
   diagnostics_.n_points_in_sphere = 0;
   diagnostics_.n_points_in_mask = 0;
   diagnostics_.n_points_outside_mask = 0;
   diagnostics_.density_min = std::numeric_limits<float>::max();
   diagnostics_.density_max = std::numeric_limits<float>::lowest();
   diagnostics_.density_sum = 0.0f;

   for (int ix = 0; ix < n_grid; ++ix) {
      float x = -radius + ix * grid_step;
      for (int iy = 0; iy < n_grid; ++iy) {
         float y = -radius + iy * grid_step;
         for (int iz = 0; iz < n_grid; ++iz) {
            float z = -radius + iz * grid_step;

            // Check if point is inside the unit ball
            float r_norm = std::sqrt(x*x + y*y + z*z) / radius;
            if (r_norm > 1.0f) continue;

            diagnostics_.n_points_in_sphere++;

            // Check if point is within mask (using relative coordinates)
            if (!is_near_ca(x, y, z)) {
               diagnostics_.n_points_outside_mask++;
               continue;
            }

            diagnostics_.n_points_in_mask++;

            // Get density at this point (absolute position in map)
            clipper::Coord_orth pt(centre.x() + x, centre.y() + y, centre.z() + z);
            float density = util::density_at_point(xmap, pt);

            // Update diagnostics
            diagnostics_.n_grid_points++;
            diagnostics_.density_sum += density;
            if (density < diagnostics_.density_min) diagnostics_.density_min = density;
            if (density > diagnostics_.density_max) diagnostics_.density_max = density;

            // Store normalized coordinates and density
            grid_values.push_back(std::make_tuple(x/radius, y/radius, z/radius, density));
         }
      }
   }

   // Finalize diagnostics
   diagnostics_.density_mean = (diagnostics_.n_grid_points > 0) ?
      diagnostics_.density_sum / diagnostics_.n_grid_points : 0.0f;

   // Compute moments from grid values
   compute_moments(grid_values);

   // Normalize moments by density sum (shape, not mass)
   normalize_moments();

   // Convert moments to rotation-invariant descriptors
   compute_invariants();
}

// Normalize moments by density sum to make descriptor shape-based, not mass-based
void ZernikeDescriptor::normalize_moments() {
   if (std::abs(diagnostics_.density_sum) < 1e-10) return;  // avoid division by zero

   double scale = 1.0 / diagnostics_.density_sum;
   for (auto &m : moments_) {
      m *= scale;
   }
}

// Compute moments from grid sample values
void ZernikeDescriptor::compute_moments(
   const std::vector<std::tuple<float,float,float,float>> &grid_values) {

   // Volume element for integration (approximate)
   double dV = 8.0 / (n_grid_ * n_grid_ * n_grid_);  // normalized to unit ball

   // For each grid point, accumulate contribution to all moments
   for (const auto &gv : grid_values) {
      float x = std::get<0>(gv);
      float y = std::get<1>(gv);
      float z = std::get<2>(gv);
      float density = std::get<3>(gv);

      // Convert to spherical coordinates
      double r = std::sqrt(x*x + y*y + z*z);
      if (r < 1e-10) {
         // At origin, only l=0 contributes
         r = 0.0;
      }

      double theta = (r > 1e-10) ? std::acos(z / r) : 0.0;
      double phi = std::atan2(y, x);

      // Accumulate contribution to each moment
      for (int n = 0; n <= order_; ++n) {
         for (int l = n; l >= 0; l -= 2) {
            for (int m = -l; m <= l; ++m) {
               std::complex<double> Z = zernike_basis(n, l, m, r, theta, phi);
               size_t idx = moment_index(n, l, m);
               // Moment = integral of f(r) * conj(Z_nl^m) dV
               moments_[idx] += static_cast<double>(density) * std::conj(Z) * dV;
            }
         }
      }
   }
}

// Convert moments to rotation-invariant descriptors
void ZernikeDescriptor::compute_invariants() {
   invariants_.clear();
   invariants_.reserve(zernike_invariant_count(order_));

   for (int n = 0; n <= order_; ++n) {
      for (int l = n; l >= 0; l -= 2) {
         // F_nl = sqrt(sum_m |Omega_nl^m|^2)
         double sum_sq = 0.0;
         for (int m = -l; m <= l; ++m) {
            size_t idx = moment_index(n, l, m);
            sum_sq += std::norm(moments_[idx]);  // |z|^2
         }
         invariants_.push_back(static_cast<float>(std::sqrt(sum_sq)));
      }
   }
}

// Get invariants
std::vector<float> ZernikeDescriptor::get_invariants() const {
   return invariants_;
}

// Get number of invariants
size_t ZernikeDescriptor::n_invariants() const {
   return invariants_.size();
}

// Map (n,l,m) to linear index
size_t ZernikeDescriptor::moment_index(int n, int l, int m) const {
   size_t idx = 0;
   // Sum over all previous n values
   for (int nn = 0; nn < n; ++nn) {
      for (int ll = nn; ll >= 0; ll -= 2) {
         idx += 2 * ll + 1;
      }
   }
   // Add contribution from current n, previous l values
   for (int ll = n; ll > l; ll -= 2) {
      idx += 2 * ll + 1;
   }
   // Add m offset (m goes from -l to +l)
   idx += m + l;
   return idx;
}

// Compute 3D Zernike basis function Z_nl^m
std::complex<double> ZernikeDescriptor::zernike_basis(int n, int l, int m,
                                                       double r, double theta, double phi) const {
   double R = radial_polynomial(n, l, r);
   std::complex<double> Y = spherical_harmonic(l, m, theta, phi);
   return R * Y;
}

// Zernike radial polynomial R_nl(r)
// R_nl(r) = sum_k Q_nlk * r^(2k+l)
double ZernikeDescriptor::radial_polynomial(int n, int l, double r) const {
   double result = 0.0;
   int k_max = (n - l) / 2;
   double r_power = std::pow(r, l);  // r^l

   for (int k = 0; k <= k_max; ++k) {
      result += Q_nlk(n, l, k) * r_power;
      r_power *= r * r;  // r^(l+2k)
   }
   return result;
}

// Precompute Q coefficients
// Q_nlk from Novotni & Klein 2003:
// Q_nlk = (-1)^k * sqrt(2n+3) * C(k, (n-l)/2) * C(2k+2l+1, k) * C((n+l+2k+3)/2, n) / C(k+l+1, k)
// where C(n,k) is binomial coefficient
void ZernikeDescriptor::precompute_Q_coefficients() const {
   Q_coeffs_.resize(order_ + 1);
   for (int n = 0; n <= order_; ++n) {
      Q_coeffs_[n].resize(n + 1);
      for (int l = n; l >= 0; l -= 2) {
         int k_max = (n - l) / 2;
         Q_coeffs_[n][l].resize(k_max + 1);
         for (int k = 0; k <= k_max; ++k) {
            // Using the formula from the literature
            double sign = (k % 2 == 0) ? 1.0 : -1.0;
            double norm = std::sqrt(2.0 * n + 3.0);

            // Compute binomial coefficients
            // C(n, k) = n! / (k! * (n-k)!)
            auto binomial = [](int n, int k) -> double {
               if (k < 0 || k > n) return 0.0;
               if (k == 0 || k == n) return 1.0;
               double result = 1.0;
               for (int i = 0; i < k; ++i) {
                  result *= (n - i);
                  result /= (i + 1);
               }
               return result;
            };

            // Calculate coefficients using recurrence or direct formula
            // Using Novotni & Klein equation (10)
            double c1 = binomial((n - l) / 2, k);
            double c2 = binomial(n, (n - l) / 2 - k);

            // Normalization factor for orthonormality
            double c3 = factorial(2 * k + l);
            double c4 = std::pow(2.0, 2 * k + l) * factorial(k) * factorial(k + l);

            Q_coeffs_[n][l][k] = sign * norm * c1 * c2 * c3 / c4;
         }
      }
   }
}

double ZernikeDescriptor::Q_nlk(int n, int l, int k) const {
   if (n < 0 || n > order_ || l < 0 || l > n || (n - l) % 2 != 0) return 0.0;
   if (k < 0 || k > (n - l) / 2) return 0.0;
   return Q_coeffs_[n][l][k];
}

// Spherical harmonic Y_l^m(theta, phi)
// Using the standard physics convention
std::complex<double> ZernikeDescriptor::spherical_harmonic(int l, int m,
                                                            double theta, double phi) const {
   // Normalization constant
   double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * M_PI) *
                           factorial(l - std::abs(m)) / factorial(l + std::abs(m)));

   // Associated Legendre polynomial
   double P = assoc_legendre(l, std::abs(m), std::cos(theta));

   // Complex exponential
   std::complex<double> exp_im_phi(std::cos(m * phi), std::sin(m * phi));

   // Condon-Shortley phase for m < 0
   double phase = 1.0;
   if (m < 0) {
      phase = (std::abs(m) % 2 == 0) ? 1.0 : -1.0;
   }

   return phase * norm * P * exp_im_phi;
}

// Associated Legendre polynomial P_l^m(x)
// Using recurrence relation for stability
double ZernikeDescriptor::assoc_legendre(int l, int m, double x) const {
   if (m < 0 || m > l) return 0.0;
   if (l < 0) return 0.0;

   // P_m^m(x) = (-1)^m * (2m-1)!! * (1-x^2)^(m/2)
   double pmm = 1.0;
   if (m > 0) {
      double somx2 = std::sqrt((1.0 - x) * (1.0 + x));
      double fact = 1.0;
      for (int i = 1; i <= m; ++i) {
         pmm *= -fact * somx2;
         fact += 2.0;
      }
   }

   if (l == m) return pmm;

   // P_{m+1}^m(x) = x * (2m + 1) * P_m^m(x)
   double pmmp1 = x * (2.0 * m + 1.0) * pmm;
   if (l == m + 1) return pmmp1;

   // Use recurrence: (l-m) * P_l^m = x * (2l-1) * P_{l-1}^m - (l+m-1) * P_{l-2}^m
   double pll = 0.0;
   for (int ll = m + 2; ll <= l; ++ll) {
      pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
      pmm = pmmp1;
      pmmp1 = pll;
   }
   return pll;
}

// Factorial
double ZernikeDescriptor::factorial(int n) {
   if (n <= 1) return 1.0;
   double result = 1.0;
   for (int i = 2; i <= n; ++i) {
      result *= i;
   }
   return result;
}

// Double factorial: n!! = n * (n-2) * (n-4) * ... * 1 or 2
double ZernikeDescriptor::double_factorial(int n) {
   if (n <= 1) return 1.0;
   double result = 1.0;
   for (int i = n; i >= 1; i -= 2) {
      result *= i;
   }
   return result;
}

// Euclidean distance between descriptors
float ZernikeDescriptor::distance(const ZernikeDescriptor &d1, const ZernikeDescriptor &d2) {
   const auto &v1 = d1.invariants_;
   const auto &v2 = d2.invariants_;

   if (v1.size() != v2.size()) {
      std::cerr << "ZernikeDescriptor::distance: descriptor size mismatch" << std::endl;
      return -1.0f;
   }

   double sum_sq = 0.0;
   for (size_t i = 0; i < v1.size(); ++i) {
      double diff = v1[i] - v2[i];
      sum_sq += diff * diff;
   }
   return static_cast<float>(std::sqrt(sum_sq));
}

// Correlation between descriptors
float ZernikeDescriptor::correlation(const ZernikeDescriptor &d1, const ZernikeDescriptor &d2) {
   const auto &v1 = d1.invariants_;
   const auto &v2 = d2.invariants_;

   if (v1.size() != v2.size() || v1.empty()) {
      return 0.0f;
   }

   // Compute means
   double mean1 = std::accumulate(v1.begin(), v1.end(), 0.0) / v1.size();
   double mean2 = std::accumulate(v2.begin(), v2.end(), 0.0) / v2.size();

   // Compute correlation
   double num = 0.0, den1 = 0.0, den2 = 0.0;
   for (size_t i = 0; i < v1.size(); ++i) {
      double d1i = v1[i] - mean1;
      double d2i = v2[i] - mean2;
      num += d1i * d2i;
      den1 += d1i * d1i;
      den2 += d2i * d2i;
   }

   double denom = std::sqrt(den1 * den2);
   if (denom < 1e-10) return 0.0f;

   return static_cast<float>(num / denom);
}

// Cosine similarity between descriptors
// cos_sim = (A · B) / (||A|| × ||B||)
// Scale-invariant: 1 = identical direction, 0 = orthogonal, -1 = opposite
float ZernikeDescriptor::cosine_similarity(const ZernikeDescriptor &d1, const ZernikeDescriptor &d2) {
   const auto &v1 = d1.invariants_;
   const auto &v2 = d2.invariants_;

   if (v1.size() != v2.size() || v1.empty()) {
      return 0.0f;
   }

   double dot = 0.0, norm1_sq = 0.0, norm2_sq = 0.0;
   for (size_t i = 0; i < v1.size(); ++i) {
      dot += v1[i] * v2[i];
      norm1_sq += v1[i] * v1[i];
      norm2_sq += v2[i] * v2[i];
   }

   double denom = std::sqrt(norm1_sq * norm2_sq);
   if (denom < 1e-10) return 0.0f;

   return static_cast<float>(dot / denom);
}

// Print diagnostic information
void ZernikeDescriptor::print_diagnostics() const {
   std::cout << "Zernike Diagnostics:" << std::endl;
   std::cout << "  Points in sphere:      " << diagnostics_.n_points_in_sphere << std::endl;
   std::cout << "  Points in mask:        " << diagnostics_.n_points_in_mask << std::endl;
   std::cout << "  Points outside mask:   " << diagnostics_.n_points_outside_mask << std::endl;
   std::cout << "  Grid points used:      " << diagnostics_.n_grid_points << std::endl;
   if (diagnostics_.n_points_in_sphere > 0) {
      float pct_in_mask = 100.0f * diagnostics_.n_points_in_mask / diagnostics_.n_points_in_sphere;
      std::cout << "  Percent in mask:       " << pct_in_mask << "%" << std::endl;
   }
   std::cout << "  Density min:  " << diagnostics_.density_min << std::endl;
   std::cout << "  Density max:  " << diagnostics_.density_max << std::endl;
   std::cout << "  Density mean: " << diagnostics_.density_mean << std::endl;
   std::cout << "  Density sum:  " << diagnostics_.density_sum << std::endl;
}

// Print first few raw moments
void ZernikeDescriptor::print_moments(int max_n) const {
   std::cout << "Raw moments (n,l,m) -> Omega_nl^m:" << std::endl;
   for (int n = 0; n <= std::min(max_n, order_); ++n) {
      for (int l = n; l >= 0; l -= 2) {
         for (int m = -l; m <= l; ++m) {
            size_t idx = moment_index(n, l, m);
            const auto &z = moments_[idx];
            std::cout << "  (" << n << "," << l << "," << std::setw(2) << m << "): "
                      << std::setw(12) << z.real() << " + "
                      << std::setw(12) << z.imag() << "i"
                      << "  |z|=" << std::abs(z) << std::endl;
         }
      }
   }
}

} // namespace coot
