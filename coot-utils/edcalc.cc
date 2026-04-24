/*
 * coot-utils/edcalc.cc
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

#include <cmath>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>

#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include <clipper/core/atomsf.h>

#include "edcalc.hh"

namespace coot {

   // Per-atom precomputed Gaussian coefficients for real-space density.
   // rho(r^2) = occ * sum_i aw[i] * exp(bw[i] * r^2)
   struct atom_edcalc_data_t {
      clipper::Coord_orth pos;
      float occ;
      float aw[6];
      float bw[6];
   };

} // namespace coot

clipper::Xmap<float>
coot::calc_atom_map_edcalc(mmdb::Manager *mol,
                            int atom_selection_handle,
                            const clipper::Cell &cell,
                            const clipper::Spacegroup &space_group,
                            const clipper::Grid_sampling &sampling,
                            float radius) {

   int nu = sampling.nu(), nv = sampling.nv(), nw = sampling.nw();

   mmdb::PPAtom sel_atoms = nullptr;
   int n_atoms;
   mol->GetSelIndex(atom_selection_handle, sel_atoms, n_atoms);

   const float rescale_b_u = 1.0f / (8.0f * M_PI * M_PI);
   const float fourpi2 = 4.0f * M_PI * M_PI;

   // Precompute per-atom Gaussian coefficients (same as clipper::AtomShapeFn::init)
   std::vector<atom_edcalc_data_t> atom_data;
   atom_data.reserve(n_atoms);

   for (int iat = 0; iat < n_atoms; iat++) {
      mmdb::Atom *at = sel_atoms[iat];
      if (at->isTer()) continue;

      atom_edcalc_data_t ad;
      ad.pos = clipper::Coord_orth(at->x, at->y, at->z);
      ad.occ = at->occupancy;

      float u_iso = at->tempFactor * rescale_b_u;
      if (u_iso < 0.1f) u_iso = 0.1f; // avoid NaN from very low B-factors

      std::string ele(at->element);
      const clipper::ScatteringFactorsData &sf = clipper::ScatteringFactors::instance()[ele];

      for (int i = 0; i < 6; i++) {
         float b_total = sf.b[i] + 8.0f * M_PI * M_PI * u_iso; // b[i] + U2B(u_iso)
         ad.bw[i] = -fourpi2 / b_total;
         ad.aw[i] = sf.a[i] * std::pow(-ad.bw[i] / M_PI, 1.5f);
      }

      atom_data.push_back(ad);
   }

   // Grid range offset for the radius (handles non-orthogonal cells correctly)
   clipper::Grid_range gd(cell, sampling, radius);
   float radius_sq = radius * radius;

   // Block sizes in grid points — at least as large as the radius extent in each direction
   int block_u = std::max(1, std::max(-gd.min().u(), gd.max().u()));
   int block_v = std::max(1, std::max(-gd.min().v(), gd.max().v()));
   int block_w = std::max(1, std::max(-gd.min().w(), gd.max().w()));

   // Assign each atom to a colour (0-26) based on its block position in the 3x3x3 tiling.
   // All atoms with the same colour are in blocks that are non-adjacent, so their
   // density contributions (which extend at most 1 block in each direction) cannot
   // write to the same grid points. This allows lock-free parallel accumulation.
   std::vector<std::vector<unsigned int>> atoms_by_colour(27);
   for (unsigned int i = 0; i < atom_data.size(); i++) {
      clipper::Coord_frac cf = atom_data[i].pos.coord_frac(cell);
      int gu = ((int)std::floor(cf.u() * nu) % nu + nu) % nu;
      int gv = ((int)std::floor(cf.v() * nv) % nv + nv) % nv;
      int gw = ((int)std::floor(cf.w() * nw) % nw + nw) % nw;
      int colour = ((gu / block_u) % 3) * 9 + ((gv / block_v) % 3) * 3 + ((gw / block_w) % 3);
      atoms_by_colour[colour].push_back(i);
   }

   // Allocate raw grid (full unit cell, no symmetry folding) for thread-safe accumulation
   size_t grid_size = static_cast<size_t>(nu) * nv * nw;
   std::vector<float> raw_grid(grid_size, 0.0f);

   // Orthogonalization matrix for inline frac->orth conversion
   clipper::Mat33<> m_orth = cell.matrix_orth();
   float m00 = m_orth(0,0), m01 = m_orth(0,1), m02 = m_orth(0,2);
   float m10 = m_orth(1,0), m11 = m_orth(1,1), m12 = m_orth(1,2);
   float m20 = m_orth(2,0), m21 = m_orth(2,1), m22 = m_orth(2,2);
   float inv_nu = 1.0f / nu, inv_nv = 1.0f / nv, inv_nw = 1.0f / nw;
   size_t stride_u = static_cast<size_t>(nv) * nw;

   // Lambda: accumulate density for a range of atoms into the raw grid
   auto process_atoms = [&](const std::vector<unsigned int> &indices,
                            unsigned int start, unsigned int end) {
      for (unsigned int ai = start; ai < end; ai++) {
         const auto &ad = atom_data[indices[ai]];

         clipper::Coord_frac cf = ad.pos.coord_frac(cell);
         clipper::Coord_grid cg = cf.coord_grid(sampling);

         int g0u = cg.u() + gd.min().u(), g1u = cg.u() + gd.max().u();
         int g0v = cg.v() + gd.min().v(), g1v = cg.v() + gd.max().v();
         int g0w = cg.w() + gd.min().w(), g1w = cg.w() + gd.max().w();

         float ax = ad.pos.x(), ay = ad.pos.y(), az = ad.pos.z();

         for (int u = g0u; u <= g1u; u++) {
            int uw = ((u % nu) + nu) % nu;
            float fu = u * inv_nu;
            for (int v = g0v; v <= g1v; v++) {
               int vw = ((v % nv) + nv) % nv;
               float fv = v * inv_nv;
               size_t uv_off = static_cast<size_t>(uw) * stride_u + static_cast<size_t>(vw) * nw;
               for (int w = g0w; w <= g1w; w++) {
                  int ww = ((w % nw) + nw) % nw;
                  float fw = w * inv_nw;
                  // grid point in orthogonal Angstroms
                  float gx = m00*fu + m01*fv + m02*fw;
                  float gy = m10*fu + m11*fv + m12*fw;
                  float gz = m20*fu + m21*fv + m22*fw;
                  float dx = gx - ax, dy = gy - ay, dz = gz - az;
                  float rsq = dx*dx + dy*dy + dz*dz;
                  if (rsq < radius_sq) {
                     float rho = ad.occ * (
                        ad.aw[0] * std::exp(ad.bw[0] * rsq) +
                        ad.aw[1] * std::exp(ad.bw[1] * rsq) +
                        ad.aw[2] * std::exp(ad.bw[2] * rsq) +
                        ad.aw[3] * std::exp(ad.bw[3] * rsq) +
                        ad.aw[4] * std::exp(ad.bw[4] * rsq) +
                        ad.aw[5] * std::exp(ad.bw[5] * rsq));
                     raw_grid[uv_off + ww] += rho;
                  }
               }
            }
         }
      }
   };

   auto tp_accum_start = std::chrono::high_resolution_clock::now();

   // Threading setup
   unsigned int n_threads = std::thread::hardware_concurrency();
   if (n_threads == 0) n_threads = 4;

   // For each of the 27 colours: all same-colour blocks are non-adjacent,
   // so their atoms can be processed in parallel without write conflicts.
   for (int colour = 0; colour < 27; colour++) {
      const auto &indices = atoms_by_colour[colour];
      unsigned int n_in_colour = indices.size();
      if (n_in_colour == 0) continue;

      if (n_in_colour < n_threads * 2 || n_threads <= 1) {
         // not enough work to justify threading overhead
         process_atoms(indices, 0, n_in_colour);
      } else {
         std::vector<std::thread> threads;
         threads.reserve(n_threads);
         unsigned int chunk = (n_in_colour + n_threads - 1) / n_threads;
         for (unsigned int t = 0; t < n_threads; t++) {
            unsigned int start = t * chunk;
            unsigned int end = std::min(start + chunk, n_in_colour);
            if (start >= end) break;
            threads.emplace_back(process_atoms, std::cref(indices), start, end);
         }
         for (auto &th : threads)
            th.join();
      }
   }

   auto tp_accum_end = std::chrono::high_resolution_clock::now();

   // Fold the raw grid into the Xmap (handles symmetry folding to ASU)
   clipper::Xmap<float> xmap;
   xmap.init(space_group, cell, sampling);
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap.first(); !ix.last(); ix.next())
      xmap[ix] = 0.0f;

   for (int u = 0; u < nu; u++) {
      for (int v = 0; v < nv; v++) {
         size_t uv_off = static_cast<size_t>(u) * stride_u + static_cast<size_t>(v) * nw;
         for (int w = 0; w < nw; w++) {
            float val = raw_grid[uv_off + w];
            if (val != 0.0f) {
               clipper::Xmap<float>::Map_reference_coord rc(xmap, clipper::Coord_grid(u, v, w));
               xmap[rc] += val;
            }
         }
      }
   }

   auto tp_mult_start = std::chrono::high_resolution_clock::now();

   // Crystallographic multiplicity correction — not needed for correlation
   // (same multiplicity applied to both calc and reference maps cancels out)
   if (false)
      for (ix = xmap.first(); !ix.last(); ix.next())
         xmap[ix] *= xmap.multiplicity(ix.coord());

   auto tp_end = std::chrono::high_resolution_clock::now();

   auto d_accum = std::chrono::duration_cast<std::chrono::milliseconds>(tp_accum_end - tp_accum_start).count();
   auto d_fold = std::chrono::duration_cast<std::chrono::milliseconds>(tp_mult_start - tp_accum_end).count();
   auto d_mult = std::chrono::duration_cast<std::chrono::milliseconds>(tp_end - tp_mult_start).count();
   std::cout << "TIMINGS:: calc_atom_map_edcalc: accum " << d_accum << " ms, fold " << d_fold
             << " ms, multiplicity " << d_mult << " ms, n_threads " << n_threads << std::endl;

   return xmap;
}
