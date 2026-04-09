/*
 * coot-utils/crowther.cc
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

#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <thread>

#include <clipper/clipper.h>
#include <clipper/core/map_interp.h>

#include "crowther.hh"
#include "utils/split-indices.hh"

using namespace coot;

// ---------------------------------------------------------------------------
// Mathematical helper functions
// ---------------------------------------------------------------------------

double crowther_t::log_factorial(int n) {
   if (n <= 1) return 0.0;
   double result = 0.0;
   for (int i=2; i<=n; i++)
      result += std::log(static_cast<double>(i));
   return result;
}

// Associated Legendre polynomial P_l^m(x).
// Works for all m in [-l, l].
double crowther_t::associated_legendre(int l, int m, double x) {

   int am = std::abs(m);
   if (am > l) return 0.0;

   // Compute P_l^{|m|}(x) via upward recurrence in l starting from P_{|m|}^{|m|}.
   double pmm = 1.0;
   if (am > 0) {
      double somx2 = std::sqrt(1.0 - x * x);
      double fact = 1.0;
      for (int i=1; i<=am; i++) {
         pmm *= -fact * somx2;
         fact += 2.0;
      }
   }

   double result = pmm;
   if (l > am) {
      double pmmp1 = x * (2.0 * am + 1.0) * pmm;
      result = pmmp1;
      for (int ll=am+2; ll<=l; ll++) {
         double pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + am - 1.0) * pmm) / (ll - am);
         pmm = pmmp1;
         pmmp1 = pll;
         result = pll;
      }
   }

   // For negative m: P_l^{-m} = (-1)^m (l-m)!/(l+m)! P_l^m
   if (m < 0) {
      double sign = (am % 2 == 0) ? 1.0 : -1.0;
      result *= sign * std::exp(log_factorial(l - am) - log_factorial(l + am));
   }
   return result;
}

// Spherical harmonic Y_l^m(theta, phi) using the convention of equation (19)
// in Foadi (2011), without Condon-Shortley phase prepended.
std::complex<double> crowther_t::spherical_harmonic(int l, int m, double theta, double phi) {

   int am = std::abs(m);
   double cos_theta = std::cos(theta);
   double plm = associated_legendre(l, am, cos_theta);
   double norm = std::sqrt((2.0 * l + 1.0) / (4.0 * M_PI) *
                           std::exp(log_factorial(l - am) - log_factorial(l + am)));

   if (m >= 0) {
      return norm * plm * std::polar(1.0, m * phi);
   } else {
      // Y_l^{-|m|} = (-1)^|m| conj(Y_l^{|m|})
      double sign = (am % 2 == 0) ? 1.0 : -1.0;
      return sign * norm * plm * std::polar(1.0, m * phi);
   }
}

// Spherical Bessel function j_l(x) via forward recurrence.
double crowther_t::spherical_bessel_j(int l, double x) {

   if (x < 1e-15) return (l == 0) ? 1.0 : 0.0;
   if (l == 0) return std::sin(x) / x;
   if (l == 1) return std::sin(x) / (x * x) - std::cos(x) / x;

   // Forward recurrence: j_{n+1}(x) = (2n+1)/x j_n(x) - j_{n-1}(x)
   // Stable when x >= l.
   if (x >= static_cast<double>(l)) {
      double jm1 = std::sin(x) / x;
      double jcurr = std::sin(x) / (x * x) - std::cos(x) / x;
      for (int n=1; n<l; n++) {
         double jnext = (2.0 * n + 1.0) / x * jcurr - jm1;
         jm1 = jcurr;
         jcurr = jnext;
      }
      return jcurr;
   }

   // Backward recurrence (Miller's algorithm) for x < l.
   int l_start = l + 15 + static_cast<int>(std::sqrt(40.0 * l));
   std::vector<double> f(l_start + 2, 0.0);
   f[l_start] = 1.0;
   for (int n=l_start-1; n>=0; n--)
      f[n] = (2.0 * (n + 1) + 1.0) / x * f[n + 1] - f[n + 2];

   if (std::abs(f[0]) < 1e-30) return 0.0;
   double scale = (std::sin(x) / x) / f[0];
   return f[l] * scale;
}

// Wigner (small) d-matrix element d^l_{m',m}(beta).
// Uses the explicit sum formula with log factorials for stability.
double crowther_t::wigner_d(int l, int mp, int m, double beta) {

   double cb2 = std::cos(beta * 0.5);
   double sb2 = std::sin(beta * 0.5);

   int s_min = std::max(0, m - mp);
   int s_max = std::min(l + m, l - mp);
   if (s_max < s_min) return 0.0;

   // Precompute the square root of the four factorials (same for all s)
   double log_prefactor = 0.5 * (log_factorial(l + mp) + log_factorial(l - mp) +
                                  log_factorial(l + m) + log_factorial(l - m));

   double sum = 0.0;
   for (int s=s_min; s<=s_max; s++) {
      int sign_exp = mp - m + s;
      double sign = ((sign_exp % 2 == 0) || (sign_exp < 0 && (-sign_exp) % 2 == 0)) ? 1.0 : -1.0;
      // Cleaner: (-1)^n for any integer n
      sign = 1.0;
      int ae = std::abs(sign_exp);
      if (ae % 2 != 0) sign = -1.0;

      double log_denom = log_factorial(s) +
                         log_factorial(l + m - s) +
                         log_factorial(l - mp - s) +
                         log_factorial(mp - m + s);

      int pow_cos = 2 * l + m - mp - 2 * s;
      int pow_sin = mp - m + 2 * s;

      double log_term = log_prefactor - log_denom;
      double term = sign * std::exp(log_term);

      if (pow_cos > 0) term *= std::pow(cb2, pow_cos);
      else if (pow_cos == 0) { /* multiply by 1 */ }
      else term *= std::pow(cb2, pow_cos);  // cb2^0 = 1, but pow_cos < 0 shouldn't happen

      if (pow_sin > 0) term *= std::pow(sb2, pow_sin);
      else if (pow_sin == 0) { /* multiply by 1 */ }
      else term *= std::pow(sb2, pow_sin);

      sum += term;
   }
   return sum;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

crowther_t::crowther_t(int ell_max, float radius, const clipper::Resolution &reso)
   : ell_max_(ell_max), radius_(radius), reso_(reso), n_radial_(50), c_coeffs_valid_(false) {
}

// ---------------------------------------------------------------------------
// Extract reflections with E-value normalisation
// ---------------------------------------------------------------------------

std::vector<crowther_t::reflection_data_t>
crowther_t::extract_reflections(const clipper::Xmap<float> &xmap) const {

   clipper::HKL_info hkls(xmap.spacegroup(), xmap.cell(), reso_, true);
   clipper::HKL_data<clipper::datatypes::F_phi<float>> fphi(hkls);
   xmap.fft_to(fphi);

   clipper::Cell cell = xmap.cell();

   // First pass: collect |F|^2 and reciprocal space coordinates
   struct raw_ref_t {
      double f_sq;
      double s;
      double theta;
      double phi;
      double inv_reso_sq;
   };
   std::vector<raw_ref_t> raw_refs;

   typedef clipper::HKL_data_base::HKL_reference_index HRI;
   for (HRI ih = fphi.first(); !ih.last(); ih.next()) {
      if (fphi[ih].missing()) continue;
      float f = fphi[ih].f();
      if (f < 1e-10f) continue;

      clipper::HKL hkl = ih.hkl();
      clipper::Coord_reci_orth s_orth = hkl.coord_reci_orth(cell);
      double sx = s_orth[0];
      double sy = s_orth[1];
      double sz = s_orth[2];
      double s = std::sqrt(sx*sx + sy*sy + sz*sz);
      if (s < 1e-10) continue;

      raw_ref_t rr;
      rr.f_sq = static_cast<double>(f) * f;
      rr.s = s;
      rr.theta = std::acos(sz / s);
      rr.phi = std::atan2(sy, sx);
      rr.inv_reso_sq = ih.invresolsq();
      raw_refs.push_back(rr);
   }

   // E-value normalisation: compute <|F|^2> in resolution shells
   int n_bins = 20;
   double inv_reso_sq_max = 1.0 / (reso_.limit() * reso_.limit());
   std::vector<double> shell_sum(n_bins, 0.0);
   std::vector<int> shell_count(n_bins, 0);

   for (const auto &rr : raw_refs) {
      int bin = static_cast<int>(rr.inv_reso_sq / inv_reso_sq_max * n_bins);
      if (bin >= n_bins) bin = n_bins - 1;
      if (bin < 0) bin = 0;
      shell_sum[bin] += rr.f_sq;
      shell_count[bin]++;
   }

   std::vector<double> shell_mean(n_bins, 1.0);
   for (int i=0; i<n_bins; i++) {
      if (shell_count[i] > 0)
         shell_mean[i] = shell_sum[i] / shell_count[i];
   }

   // Build output with |E|^2 = |F|^2 / <|F|^2>
   std::vector<reflection_data_t> refs;
   refs.reserve(raw_refs.size());
   for (const auto &rr : raw_refs) {
      int bin = static_cast<int>(rr.inv_reso_sq / inv_reso_sq_max * n_bins);
      if (bin >= n_bins) bin = n_bins - 1;
      if (bin < 0) bin = 0;
      reflection_data_t rd;
      rd.e_sq = rr.f_sq / shell_mean[bin];
      rd.s = rr.s;
      rd.theta = rr.theta;
      rd.phi = rr.phi;
      refs.push_back(rd);
   }

   std::cout << "INFO:: crowther_t::extract_reflections() " << refs.size() << " reflections" << std::endl;
   return refs;
}

// ---------------------------------------------------------------------------
// Precompute Y_lm for all reflections
// ---------------------------------------------------------------------------

std::vector<std::vector<std::complex<double>>>
crowther_t::precompute_ylm(const std::vector<reflection_data_t> &refs) const {

   int n_lm = (ell_max_ + 1) * (ell_max_ + 1);
   std::vector<std::vector<std::complex<double>>> ylm(refs.size(),
      std::vector<std::complex<double>>(n_lm, 0.0));

   for (size_t i=0; i<refs.size(); i++) {
      double th = refs[i].theta;
      double ph = refs[i].phi;
      for (int l=0; l<=ell_max_; l++) {
         for (int m=-l; m<=l; m++) {
            int idx = l * (l + 1) + m;
            ylm[i][idx] = spherical_harmonic(l, m, th, ph);
         }
      }
   }
   return ylm;
}

// ---------------------------------------------------------------------------
// Compute a'_lm(r) = sum_s E^2(s) j_l(2*pi*s*r) Y_l^{m*}(theta_s, phi_s)
// ---------------------------------------------------------------------------

std::vector<std::complex<double>>
crowther_t::compute_alm_at_r(const std::vector<reflection_data_t> &refs,
                              const std::vector<std::vector<std::complex<double>>> &ylm_cache,
                              double r) const {

   int n_lm = (ell_max_ + 1) * (ell_max_ + 1);
   std::vector<std::complex<double>> alm(n_lm, 0.0);

   for (size_t i=0; i<refs.size(); i++) {
      double e_sq = refs[i].e_sq;
      double x = 2.0 * M_PI * refs[i].s * r;

      // Compute j_l(x) for all l
      std::vector<double> jl(ell_max_ + 1);
      for (int l=0; l<=ell_max_; l++)
         jl[l] = spherical_bessel_j(l, x);

      for (int l=0; l<=ell_max_; l++) {
         double coeff = e_sq * jl[l];
         for (int m=-l; m<=l; m++) {
            int idx = l * (l + 1) + m;
            alm[idx] += coeff * std::conj(ylm_cache[i][idx]);
         }
      }
   }

   return alm;
}

// ---------------------------------------------------------------------------
// Set target / search maps
// ---------------------------------------------------------------------------

void crowther_t::set_target(const clipper::Xmap<float> &xmap) {
   target_refs_ = extract_reflections(xmap);
}

void crowther_t::set_search(const clipper::Xmap<float> &xmap) {
   search_refs_ = extract_reflections(xmap);
}

// ---------------------------------------------------------------------------
// Diagnostics: print expansion coefficients for validation
// ---------------------------------------------------------------------------

void crowther_t::print_diagnostics() const {

   std::cout << "=== Crowther rotation function diagnostics ===" << std::endl;
   std::cout << "ell_max=" << ell_max_ << " radius=" << radius_
             << " reso=" << reso_.limit() << " n_radial=" << n_radial_ << std::endl;
   std::cout << "target reflections: " << target_refs_.size()
             << "  search reflections: " << search_refs_.size() << std::endl;

   // Check E^2 statistics for both sets
   for (int which=0; which<2; which++) {
      const auto &refs = (which == 0) ? target_refs_ : search_refs_;
      const char *label = (which == 0) ? "target" : "search";

      double e2_sum = 0.0, e2_min = 1e30, e2_max = -1e30;
      double s_min = 1e30, s_max = -1e30;
      for (const auto &r : refs) {
         e2_sum += r.e_sq;
         if (r.e_sq < e2_min) e2_min = r.e_sq;
         if (r.e_sq > e2_max) e2_max = r.e_sq;
         if (r.s < s_min) s_min = r.s;
         if (r.s > s_max) s_max = r.s;
      }
      double e2_mean = e2_sum / refs.size();
      std::cout << label << " E^2: mean=" << std::fixed << std::setprecision(4) << e2_mean
                << " min=" << e2_min << " max=" << e2_max << std::endl;
      std::cout << label << " |s| range: [" << s_min << ", " << s_max << "]" << std::endl;
   }

   // Compute a_lm at a few radial values and print magnitudes by l
   auto ylm_target = precompute_ylm(target_refs_);
   auto ylm_search = precompute_ylm(search_refs_);

   std::cout << "\n--- a_lm magnitudes by ell (summed |a_lm|^2 over m) ---" << std::endl;
   std::cout << std::setw(6) << "r(A)";
   for (int l=0; l<=std::min(ell_max_, 10); l++)
      std::cout << std::setw(14) << ("l=" + std::to_string(l));
   std::cout << std::endl;

   double dr = radius_ / n_radial_;
   for (int ir=0; ir<n_radial_; ir+=10) {
      double r = (ir + 0.5) * dr;
      auto alm_t = compute_alm_at_r(target_refs_, ylm_target, r);
      auto alm_s = compute_alm_at_r(search_refs_, ylm_search, r);

      std::cout << "target r=" << std::fixed << std::setprecision(2) << std::setw(6) << r;
      for (int l=0; l<=std::min(ell_max_, 10); l++) {
         double sum_sq = 0.0;
         for (int m=-l; m<=l; m++) {
            int idx = l * (l + 1) + m;
            sum_sq += std::norm(alm_t[idx]);
         }
         std::cout << std::scientific << std::setprecision(3) << std::setw(14) << std::sqrt(sum_sq);
      }
      std::cout << std::endl;

      std::cout << "search r=" << std::fixed << std::setprecision(2) << std::setw(6) << r;
      for (int l=0; l<=std::min(ell_max_, 10); l++) {
         double sum_sq = 0.0;
         for (int m=-l; m<=l; m++) {
            int idx = l * (l + 1) + m;
            sum_sq += std::norm(alm_s[idx]);
         }
         std::cout << std::scientific << std::setprecision(3) << std::setw(14) << std::sqrt(sum_sq);
      }
      std::cout << std::endl;
   }

   // Compute c_lmm' and report |c_l| = sqrt(sum_{m,m'} |c_{lmm'}|^2) for each l
   std::vector<std::vector<std::vector<std::complex<double>>>> c(ell_max_ + 1);
   for (int l=0; l<=ell_max_; l++) {
      int dim = 2 * l + 1;
      c[l].resize(dim, std::vector<std::complex<double>>(dim, 0.0));
   }

   for (int ir=0; ir<n_radial_; ir++) {
      double r = (ir + 0.5) * dr;
      double r2dr = r * r * dr;
      auto alm = compute_alm_at_r(target_refs_, ylm_target, r);
      auto blm = compute_alm_at_r(search_refs_, ylm_search, r);
      for (int l=0; l<=ell_max_; l++) {
         for (int m=-l; m<=l; m++) {
            int idx_m = l * (l + 1) + m;
            for (int mp=-l; mp<=l; mp++) {
               int idx_mp = l * (l + 1) + mp;
               c[l][m + l][mp + l] += std::conj(alm[idx_m]) * blm[idx_mp] * r2dr;
            }
         }
      }
   }

   std::cout << "\n--- |c_l| = sqrt(sum_{m,m'} |c_{lmm'}|^2) ---" << std::endl;
   for (int l=0; l<=ell_max_; l++) {
      double sum_sq = 0.0;
      for (int m=-l; m<=l; m++) {
         for (int mp=-l; mp<=l; mp++) {
            sum_sq += std::norm(c[l][m + l][mp + l]);
         }
      }
      std::cout << "l=" << std::setw(3) << l
                << " |c_l|=" << std::scientific << std::setprecision(6) << std::sqrt(sum_sq)
                << std::endl;
   }

   std::cout << "=== end diagnostics ===" << std::endl;
}

// ---------------------------------------------------------------------------
// Precompute c_{l,m,m'} coefficients from radial integration
// ---------------------------------------------------------------------------

void crowther_t::precompute_c_coefficients() {

   auto t0 = std::chrono::high_resolution_clock::now();
   auto ylm_target = precompute_ylm(target_refs_);
   auto ylm_search = precompute_ylm(search_refs_);

   // c_{l,m,m'} = integral_0^R conj(a_{lm}(r)) * b_{lm'}(r) * r^2 dr
   // Parallelise over radial shells, each thread accumulates its own c_coeffs,
   // then merge.

   unsigned int n_threads = std::thread::hardware_concurrency();
   if (n_threads == 0) n_threads = 4;
   if (n_threads > static_cast<unsigned int>(n_radial_)) n_threads = n_radial_;

   auto ranges = coot::atom_index_ranges(n_radial_, n_threads);

   // Per-thread local c_coeffs
   using c_type = std::vector<std::vector<std::vector<std::complex<double>>>>;
   std::vector<c_type> thread_c(ranges.size());
   for (unsigned int t=0; t<ranges.size(); t++) {
      thread_c[t].resize(ell_max_ + 1);
      for (int l=0; l<=ell_max_; l++) {
         int dim = 2 * l + 1;
         thread_c[t][l].assign(dim, std::vector<std::complex<double>>(dim, 0.0));
      }
   }

   double dr = radius_ / n_radial_;

   auto worker = [&](unsigned int thread_idx, unsigned int ir_start, unsigned int ir_end) {
      auto &local_c = thread_c[thread_idx];
      for (unsigned int ir=ir_start; ir<ir_end; ir++) {
         double r = (ir + 0.5) * dr;
         double r2dr = r * r * dr;

         auto alm = compute_alm_at_r(target_refs_, ylm_target, r);
         auto blm = compute_alm_at_r(search_refs_, ylm_search, r);

         for (int l=0; l<=ell_max_; l++) {
            for (int m=-l; m<=l; m++) {
               int idx_m = l * (l + 1) + m;
               for (int mp=-l; mp<=l; mp++) {
                  int idx_mp = l * (l + 1) + mp;
                  local_c[l][m + l][mp + l] += std::conj(alm[idx_m]) * blm[idx_mp] * r2dr;
               }
            }
         }
      }
   };

   std::vector<std::thread> threads;
   for (unsigned int i=0; i<ranges.size(); i++)
      threads.emplace_back(worker, i, ranges[i].first, ranges[i].second);
   for (auto &t : threads)
      t.join();

   // Merge thread results
   c_coeffs_.resize(ell_max_ + 1);
   for (int l=0; l<=ell_max_; l++) {
      int dim = 2 * l + 1;
      c_coeffs_[l].assign(dim, std::vector<std::complex<double>>(dim, 0.0));
      for (unsigned int t=0; t<ranges.size(); t++) {
         for (int m=0; m<dim; m++) {
            for (int mp=0; mp<dim; mp++) {
               c_coeffs_[l][m][mp] += thread_c[t][l][m][mp];
            }
         }
      }
   }

   c_coeffs_valid_ = true;
   auto t1 = std::chrono::high_resolution_clock::now();
   double c_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
   std::cout << "INFO:: crowther_t c_{lmm'} coefficients computed in "
             << std::fixed << std::setprecision(1) << c_ms << " ms ("
             << ranges.size() << " threads)" << std::endl;
}

// ---------------------------------------------------------------------------
// Evaluate the rotation function at a single orientation (ZYZ Euler, radians)
// ---------------------------------------------------------------------------

double crowther_t::evaluate_at_euler(double alpha, double beta, double gamma) const {

   if (!c_coeffs_valid_) {
      std::cerr << "ERROR:: evaluate_at_euler() called before c coefficients computed" << std::endl;
      return 0.0;
   }

   // C_{m,m'}(beta) = sum_l c_{lmm'} d^l_{mm'}(beta)
   // R(alpha,beta,gamma) = sum_{m,m'} C_{mm'} exp(-i*m*alpha) exp(-i*m'*gamma)

   std::complex<double> val(0.0, 0.0);
   for (int l=2; l<=ell_max_; l++) {
      for (int m=-l; m<=l; m++) {
         std::complex<double> exp_ma = std::polar(1.0, -m * alpha);
         for (int mp=-l; mp<=l; mp++) {
            double d = wigner_d(l, m, mp, beta);
            std::complex<double> exp_mpg = std::polar(1.0, -mp * gamma);
            val += c_coeffs_[l][m + l][mp + l] * d * exp_ma * exp_mpg;
         }
      }
   }
   return val.real();
}

// ---------------------------------------------------------------------------
// Refine orientations with a 3x3x3 grid search in ZYZ Euler angles
// ---------------------------------------------------------------------------

std::vector<rotation_function_result_t>
crowther_t::refine_orientations(const std::vector<rotation_function_result_t> &coarse,
                                int n_refine, float step_degrees) const {

   if (!c_coeffs_valid_) {
      std::cerr << "ERROR:: refine_orientations() called before c coefficients computed" << std::endl;
      return coarse;
   }

   double step = glm::radians(static_cast<double>(step_degrees));
   int n = std::min(n_refine, static_cast<int>(coarse.size()));

   std::vector<rotation_function_result_t> refined(n);

   unsigned int n_threads = std::thread::hardware_concurrency();
   if (n_threads == 0) n_threads = 4;
   if (n_threads > static_cast<unsigned int>(n)) n_threads = n;

   auto ranges = coot::atom_index_ranges(n, n_threads);
   std::vector<std::thread> threads;
   threads.reserve(ranges.size());

   for (unsigned int t=0; t<ranges.size(); t++) {
      threads.emplace_back([&, t]() {
         for (unsigned int i=ranges[t].first; i<ranges[t].second; i++) {
            const glm::quat &q = coarse[i].rotation;

            // Convert quaternion back to ZYZ Euler angles.
            glm::mat3 R = glm::mat3_cast(q);
            double cos_beta = static_cast<double>(R[2][2]);
            cos_beta = std::clamp(cos_beta, -1.0, 1.0);
            double beta_0 = std::acos(cos_beta);

            double alpha_0, gamma_0;
            if (std::abs(std::sin(beta_0)) > 1e-6) {
               alpha_0 = std::atan2(static_cast<double>(R[2][1]),
                                    static_cast<double>(R[2][0]));
               gamma_0 = std::atan2(static_cast<double>(R[1][2]),
                                    static_cast<double>(-R[0][2]));
            } else {
               alpha_0 = std::atan2(static_cast<double>(R[1][0]),
                                    static_cast<double>(R[0][0]));
               gamma_0 = 0.0;
            }

            // Search the 3x3x3 grid
            double best_score = -1e30;
            glm::quat best_quat = q;

            for (int da=-1; da<=1; da++) {
               double alpha = alpha_0 + da * step;
               for (int db=-1; db<=1; db++) {
                  double beta = beta_0 + db * step;
                  for (int dg=-1; dg<=1; dg++) {
                     double gamma = gamma_0 + dg * step;

                     double score = evaluate_at_euler(alpha, beta, gamma);
                     if (score > best_score) {
                        best_score = score;
                        float fa = static_cast<float>(alpha);
                        float fb = static_cast<float>(beta);
                        float fg = static_cast<float>(gamma);
                        glm::quat qza = glm::angleAxis(fa, glm::vec3(0, 0, 1));
                        glm::quat qyb = glm::angleAxis(fb, glm::vec3(0, 1, 0));
                        glm::quat qzg = glm::angleAxis(fg, glm::vec3(0, 0, 1));
                        best_quat = qza * qyb * qzg;
                     }
                  }
               }
            }

            refined[i] = rotation_function_result_t(best_quat, static_cast<float>(best_score));
         }
      });
   }
   for (auto &th : threads) th.join();

   // Sort refined results by descending score
   std::sort(refined.begin(), refined.end(),
             [](const rotation_function_result_t &a, const rotation_function_result_t &b) {
                return a.score > b.score;
             });

   return refined;
}

// ---------------------------------------------------------------------------
// Compute the rotation function
// ---------------------------------------------------------------------------

std::vector<rotation_function_result_t>
crowther_t::compute(int n_beta) {

   int N = 2 * ell_max_ + 1;  // grid size for alpha and gamma

   std::cout << "INFO:: crowther_t::compute() ell_max=" << ell_max_
             << " n_beta=" << n_beta << " N_alpha_gamma=" << N
             << " n_reflections: target=" << target_refs_.size()
             << " search=" << search_refs_.size() << std::endl;

   // --- Step 1+2: Precompute c_{l,m,m'} coefficients ---
   if (!c_coeffs_valid_)
      precompute_c_coefficients();

   // --- Step 3: For each beta, compute C_{m,m'}(beta) = sum_l c_{lmm'} d^l_{mm'}(beta) ---
   // --- Step 4: Evaluate R(alpha, beta, gamma) by 2D DFT ---
   // Parallelise over beta values using std::thread.

   unsigned int n_threads = std::thread::hardware_concurrency();
   if (n_threads == 0) n_threads = 4;
   if (n_threads > static_cast<unsigned int>(n_beta)) n_threads = n_beta;

   auto ranges = coot::atom_index_ranges(n_beta, n_threads);
   std::vector<std::vector<rotation_function_result_t>> thread_results(ranges.size());

   // Capture what each thread needs by value/reference
   auto worker = [&](unsigned int thread_idx, unsigned int ib_start, unsigned int ib_end) {
      auto &local_results = thread_results[thread_idx];
      local_results.reserve((ib_end - ib_start) * N * N);

      for (unsigned int ib=ib_start; ib<ib_end; ib++) {
         double beta = M_PI * (ib + 0.5) / n_beta;

         // C_{m,m'}(beta) indexed by [m + ell_max_][m' + ell_max_]
         std::vector<std::vector<std::complex<double>>> C(N,
            std::vector<std::complex<double>>(N, 0.0));

         // Start from l=2: l=0 is constant, l=1 is dipole
         for (int l=2; l<=ell_max_; l++) {
            for (int m=-l; m<=l; m++) {
               for (int mp=-l; mp<=l; mp++) {
                  double d = wigner_d(l, m, mp, beta);
                  C[m + ell_max_][mp + ell_max_] += c_coeffs_[l][m + l][mp + l] * d;
               }
            }
         }

         // 2D DFT: R(alpha_j, beta, gamma_k) = sum_{m,m'} C_{mm'} exp(-i*m*alpha_j - i*m'*gamma_k)
         for (int ja=0; ja<N; ja++) {
            double alpha = 2.0 * M_PI * ja / N;
            for (int jg=0; jg<N; jg++) {
               double gamma = 2.0 * M_PI * jg / N;

               std::complex<double> val(0.0, 0.0);
               for (int m=-ell_max_; m<=ell_max_; m++) {
                  std::complex<double> exp_ma = std::polar(1.0, -m * alpha);
                  for (int mp=-ell_max_; mp<=ell_max_; mp++) {
                     std::complex<double> exp_mpg = std::polar(1.0, -mp * gamma);
                     val += C[m + ell_max_][mp + ell_max_] * exp_ma * exp_mpg;
                  }
               }

               float score = static_cast<float>(val.real());

               // Convert ZYZ Euler angles to quaternion: Omega = Rz(alpha) Ry(beta) Rz(gamma)
               float fa = static_cast<float>(alpha);
               float fb = static_cast<float>(beta);
               float fg = static_cast<float>(gamma);
               glm::quat qza = glm::angleAxis(fa, glm::vec3(0, 0, 1));
               glm::quat qyb = glm::angleAxis(fb, glm::vec3(0, 1, 0));
               glm::quat qzg = glm::angleAxis(fg, glm::vec3(0, 0, 1));
               glm::quat q = qza * qyb * qzg;

               local_results.push_back(rotation_function_result_t(q, score));
            }
         }
      }
   };

   std::cout << "INFO:: crowther_t::compute() using " << ranges.size() << " threads" << std::endl;

   std::vector<std::thread> threads;
   for (unsigned int i=0; i<ranges.size(); i++)
      threads.emplace_back(worker, i, ranges[i].first, ranges[i].second);
   for (auto &t : threads)
      t.join();

   // Merge thread results
   std::vector<rotation_function_result_t> results;
   results.reserve(n_beta * N * N);
   for (auto &tr : thread_results)
      results.insert(results.end(), tr.begin(), tr.end());

   // Sort by descending score
   std::sort(results.begin(), results.end(),
             [](const rotation_function_result_t &a, const rotation_function_result_t &b) {
                return a.score > b.score;
             });

   std::cout << "INFO:: crowther_t::compute() " << results.size() << " results" << std::endl;
   return results;
}
