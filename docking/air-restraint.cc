/*
 * docking/air-restraint.cc
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

#include "air-restraint.hh"

// ---------------------------------------------------------------------------
// AIR effective distance
// ---------------------------------------------------------------------------

// d_eff = ( sum_m sum_k sum_n  1/d_{m,n}^6 )^(-1/6)
//
// where m indexes atoms of the active residue on protein A,
// k indexes partner residues (active+passive) on protein B,
// n indexes atoms of partner residue k.
//
// The 1/r^6 averaging mimics the attractive part of a Lennard-Jones
// potential and ensures that the AIR is satisfied as soon as any two
// atoms of the two proteins are in contact.

double
coot::haddock::air_effective_distance(const residue_atoms_t &active_residue_A,
                                      const std::vector<residue_atoms_t> &partner_residues_B) {

   double sum_inv_r6 = 0.0;

   for (unsigned int m=0; m<active_residue_A.positions.size(); m++) {
      const clipper::Coord_orth &pos_m = active_residue_A.positions[m];
      for (unsigned int k=0; k<partner_residues_B.size(); k++) {
         const residue_atoms_t &partner_res = partner_residues_B[k];
         for (unsigned int n=0; n<partner_res.positions.size(); n++) {
            const clipper::Coord_orth &pos_n = partner_res.positions[n];
            double dx = pos_m.x() - pos_n.x();
            double dy = pos_m.y() - pos_n.y();
            double dz = pos_m.z() - pos_n.z();
            double dist_sq = dx*dx + dy*dy + dz*dz;
            // Stabilise: clamp minimum distance to 0.5 A
            if (dist_sq < 0.25) dist_sq = 0.25;
            double dist_sq_cubed = dist_sq * dist_sq * dist_sq; // r^6
            sum_inv_r6 += 1.0 / dist_sq_cubed;
         }
      }
   }

   if (sum_inv_r6 < 1.0e-30) return 1.0e5; // effectively infinite distance

   // d_eff = sum^(-1/6)
   double d_eff = std::pow(sum_inv_r6, -1.0/6.0);
   return d_eff;
}

// ---------------------------------------------------------------------------
// AIR penalty function
// ---------------------------------------------------------------------------

// Flat-bottom harmonic: zero when d_eff <= distance_cutoff,
// (d_eff - distance_cutoff)^2 otherwise.

double
coot::haddock::air_penalty(double d_eff, double distance_cutoff) {

   if (d_eff <= distance_cutoff) return 0.0;
   double delta = d_eff - distance_cutoff;
   return delta * delta;
}

// ---------------------------------------------------------------------------
// AIR energy and gradients
// ---------------------------------------------------------------------------

// For each active residue i on protein A, we have one AIR to the
// partner residues on protein B.  The total AIR energy is:
//
//   E_AIR = sum_i  penalty(d_eff_i)
//
// where penalty is the flat-bottom harmonic above.
//
// The gradient of the penalty with respect to atom position x_m
// (on protein A, in active residue i) is:
//
//   dE/dx_m = 2 * (d_eff - d0) * dd_eff/dx_m    (when d_eff > d0)
//
// where dd_eff/dx_m = d_eff^7 * sum_n  6*(x_m - x_n) / d_{m,n}^8
//
// Similarly, the gradient with respect to x_n (on protein B) is:
//
//   dE/dx_n = 2 * (d_eff - d0) * dd_eff/dx_n
//
// where dd_eff/dx_n = d_eff^7 * sum_m  6*(x_n - x_m) / d_{m,n}^8

double
coot::haddock::air_energy_and_gradients(const std::vector<residue_atoms_t> &active_residues_A,
                                        const std::vector<residue_atoms_t> &partner_residues_B,
                                        double distance_cutoff,
                                        std::vector<residue_atoms_t> *grad_A,
                                        std::vector<residue_atoms_t> *grad_B) {

   double total_energy = 0.0;

   for (unsigned int i=0; i<active_residues_A.size(); i++) {
      const residue_atoms_t &res_A = active_residues_A[i];

      // First pass: compute sum_inv_r6 and d_eff
      double sum_inv_r6 = 0.0;
      for (unsigned int m=0; m<res_A.positions.size(); m++) {
         const clipper::Coord_orth &pos_m = res_A.positions[m];
         for (unsigned int k=0; k<partner_residues_B.size(); k++) {
            const residue_atoms_t &res_B = partner_residues_B[k];
            for (unsigned int n=0; n<res_B.positions.size(); n++) {
               const clipper::Coord_orth &pos_n = res_B.positions[n];
               double dx = pos_m.x() - pos_n.x();
               double dy = pos_m.y() - pos_n.y();
               double dz = pos_m.z() - pos_n.z();
               double dist_sq = dx*dx + dy*dy + dz*dz;
               if (dist_sq < 0.25) dist_sq = 0.25;
               double r6 = dist_sq * dist_sq * dist_sq;
               sum_inv_r6 += 1.0 / r6;
            }
         }
      }

      if (sum_inv_r6 < 1.0e-30) continue;

      double d_eff = std::pow(sum_inv_r6, -1.0/6.0);
      double penalty = air_penalty(d_eff, distance_cutoff);
      total_energy += penalty;

      // Gradients only needed when the restraint is violated
      if (d_eff <= distance_cutoff) continue;
      if (!grad_A && !grad_B) continue;

      // Prefactor: 2 * (d_eff - d0) * d_eff^7
      // The d_eff^7 comes from differentiating d_eff = S^(-1/6):
      //   dd_eff/dS = (-1/6) * S^(-7/6) = (-1/6) * d_eff^7 / d_eff^0
      // Actually: d(S^(-1/6))/dS = (-1/6) * S^(-7/6)
      // and S = sum_inv_r6, so dd_eff/dS = (-1/6) * sum_inv_r6^(-7/6)
      //
      // dd_eff/d(d_{mn}) = dd_eff/dS * dS/d(d_{mn})
      // dS/d(d_{mn}) = -6 / d_{mn}^7
      //
      // So dd_eff/d(d_{mn}) = (-1/6) * S^(-7/6) * (-6 / d_{mn}^7)
      //                     = S^(-7/6) / d_{mn}^7
      //                     = d_eff^7 / (d_eff^0 ... no)
      //
      // Let me redo this carefully.
      //
      // S = sum_{m,k,n} d_{mn}^{-6}
      // d_eff = S^{-1/6}
      //
      // dd_eff/dx_m = dd_eff/dS * dS/dx_m
      //
      // dd_eff/dS = (-1/6) * S^{-7/6}
      //
      // dS/dx_m = sum_{k,n} d/dx_m [ d_{mn}^{-6} ]
      //         = sum_{k,n} (-6) * d_{mn}^{-7} * d(d_{mn})/dx_m
      //         = sum_{k,n} (-6) * d_{mn}^{-7} * (x_m - x_n)/d_{mn}
      //         = sum_{k,n} (-6) * (x_m - x_n) / d_{mn}^8
      //
      // dd_eff/dx_m = (-1/6) * S^{-7/6} * sum_{k,n} (-6) * (x_m - x_n) / d_{mn}^8
      //             = S^{-7/6} * sum_{k,n} (x_m - x_n) / d_{mn}^8
      //
      // And S^{-7/6} = (S^{-1/6})^7 = d_eff^7
      //
      // So dd_eff/dx_m = d_eff^7 * sum_{k,n} (x_m - x_n) / d_{mn}^8
      //
      // The full gradient of the penalty:
      // dE/dx_m = 2*(d_eff - d0) * d_eff^7 * sum_{k,n} (x_m - x_n) / d_{mn}^8

      double delta = d_eff - distance_cutoff;
      double d_eff_7 = d_eff * d_eff * d_eff * d_eff * d_eff * d_eff * d_eff;
      double prefactor = 2.0 * delta * d_eff_7;

      // Second pass: accumulate gradients
      for (unsigned int m=0; m<res_A.positions.size(); m++) {
         const clipper::Coord_orth &pos_m = res_A.positions[m];
         double gx_m = 0.0, gy_m = 0.0, gz_m = 0.0;

         for (unsigned int k=0; k<partner_residues_B.size(); k++) {
            const residue_atoms_t &res_B = partner_residues_B[k];
            for (unsigned int n=0; n<res_B.positions.size(); n++) {
               const clipper::Coord_orth &pos_n = res_B.positions[n];
               double dx = pos_m.x() - pos_n.x();
               double dy = pos_m.y() - pos_n.y();
               double dz = pos_m.z() - pos_n.z();
               double dist_sq = dx*dx + dy*dy + dz*dz;
               if (dist_sq < 0.25) dist_sq = 0.25;
               // 1/d^8 = 1/(d^2)^4
               double dist_sq_sq = dist_sq * dist_sq;
               double inv_r8 = 1.0 / (dist_sq_sq * dist_sq_sq);

               double contrib_x = dx * inv_r8;
               double contrib_y = dy * inv_r8;
               double contrib_z = dz * inv_r8;

               gx_m += contrib_x;
               gy_m += contrib_y;
               gz_m += contrib_z;

               // Gradient on atom n of protein B is the negative
               if (grad_B) {
                  (*grad_B)[k].positions[n] = clipper::Coord_orth(
                     (*grad_B)[k].positions[n].x() + prefactor * (-contrib_x),
                     (*grad_B)[k].positions[n].y() + prefactor * (-contrib_y),
                     (*grad_B)[k].positions[n].z() + prefactor * (-contrib_z));
               }
            }
         }

         if (grad_A) {
            (*grad_A)[i].positions[m] = clipper::Coord_orth(
               (*grad_A)[i].positions[m].x() + prefactor * gx_m,
               (*grad_A)[i].positions[m].y() + prefactor * gy_m,
               (*grad_A)[i].positions[m].z() + prefactor * gz_m);
         }
      }
   }

   return total_energy;
}

// ---------------------------------------------------------------------------
// Numerical gradient check
// ---------------------------------------------------------------------------

void
coot::haddock::air_numerical_gradient_check(const std::vector<residue_atoms_t> &active_residues_A,
                                            const std::vector<residue_atoms_t> &partner_residues_B,
                                            double distance_cutoff,
                                            double micro_step) {

   // Compute analytical gradients
   std::vector<residue_atoms_t> grad_A(active_residues_A.size());
   for (unsigned int i=0; i<active_residues_A.size(); i++) {
      grad_A[i].positions.resize(active_residues_A[i].positions.size(),
                                 clipper::Coord_orth(0,0,0));
   }
   std::vector<residue_atoms_t> grad_B(partner_residues_B.size());
   for (unsigned int k=0; k<partner_residues_B.size(); k++) {
      grad_B[k].positions.resize(partner_residues_B[k].positions.size(),
                                 clipper::Coord_orth(0,0,0));
   }

   air_energy_and_gradients(active_residues_A, partner_residues_B,
                            distance_cutoff, &grad_A, &grad_B);

   std::cout << "AIR numerical gradient check (micro_step=" << micro_step << ")\n";
   std::cout << std::setw(12) << "component"
             << std::setw(14) << "analytical"
             << std::setw(14) << "numerical"
             << std::setw(14) << "ratio"
             << std::endl;

   // Check gradients on protein A atoms
   for (unsigned int i=0; i<active_residues_A.size(); i++) {
      for (unsigned int m=0; m<active_residues_A[i].positions.size(); m++) {
         for (int coord=0; coord<3; coord++) {
            // Perturb +
            std::vector<residue_atoms_t> perturbed_A = active_residues_A;
            clipper::Coord_orth &p = perturbed_A[i].positions[m];
            double orig = (coord == 0) ? p.x() : (coord == 1) ? p.y() : p.z();
            if (coord == 0) p = clipper::Coord_orth(orig + micro_step, p.y(), p.z());
            if (coord == 1) p = clipper::Coord_orth(p.x(), orig + micro_step, p.z());
            if (coord == 2) p = clipper::Coord_orth(p.x(), p.y(), orig + micro_step);
            double e_plus = air_energy_and_gradients(perturbed_A, partner_residues_B,
                                                     distance_cutoff, nullptr, nullptr);

            // Perturb -
            perturbed_A = active_residues_A;
            clipper::Coord_orth &p2 = perturbed_A[i].positions[m];
            if (coord == 0) p2 = clipper::Coord_orth(orig - micro_step, p2.y(), p2.z());
            if (coord == 1) p2 = clipper::Coord_orth(p2.x(), orig - micro_step, p2.z());
            if (coord == 2) p2 = clipper::Coord_orth(p2.x(), p2.y(), orig - micro_step);
            double e_minus = air_energy_and_gradients(perturbed_A, partner_residues_B,
                                                      distance_cutoff, nullptr, nullptr);

            double numerical = (e_plus - e_minus) / (2.0 * micro_step);
            double analytical = (coord == 0) ? grad_A[i].positions[m].x() :
                                (coord == 1) ? grad_A[i].positions[m].y() :
                                               grad_A[i].positions[m].z();

            double ratio = (std::abs(numerical) > 1.0e-10) ? analytical / numerical : 0.0;
            std::string label = "A[" + std::to_string(i) + "][" + std::to_string(m) + "]."
                                + std::string(1, "xyz"[coord]);

            std::cout << std::setw(12) << label
                      << std::setw(14) << std::setprecision(6) << std::fixed << analytical
                      << std::setw(14) << std::setprecision(6) << std::fixed << numerical
                      << std::setw(14) << std::setprecision(4) << std::fixed << ratio
                      << std::endl;
         }
      }
   }

   // Check gradients on protein B atoms
   for (unsigned int k=0; k<partner_residues_B.size(); k++) {
      for (unsigned int n=0; n<partner_residues_B[k].positions.size(); n++) {
         for (int coord=0; coord<3; coord++) {
            std::vector<residue_atoms_t> perturbed_B = partner_residues_B;
            clipper::Coord_orth &p = perturbed_B[k].positions[n];
            double orig = (coord == 0) ? p.x() : (coord == 1) ? p.y() : p.z();
            if (coord == 0) p = clipper::Coord_orth(orig + micro_step, p.y(), p.z());
            if (coord == 1) p = clipper::Coord_orth(p.x(), orig + micro_step, p.z());
            if (coord == 2) p = clipper::Coord_orth(p.x(), p.y(), orig + micro_step);
            double e_plus = air_energy_and_gradients(active_residues_A, perturbed_B,
                                                     distance_cutoff, nullptr, nullptr);

            perturbed_B = partner_residues_B;
            clipper::Coord_orth &p2 = perturbed_B[k].positions[n];
            if (coord == 0) p2 = clipper::Coord_orth(orig - micro_step, p2.y(), p2.z());
            if (coord == 1) p2 = clipper::Coord_orth(p2.x(), orig - micro_step, p2.z());
            if (coord == 2) p2 = clipper::Coord_orth(p2.x(), p2.y(), orig - micro_step);
            double e_minus = air_energy_and_gradients(active_residues_A, perturbed_B,
                                                      distance_cutoff, nullptr, nullptr);

            double numerical = (e_plus - e_minus) / (2.0 * micro_step);
            double analytical = (coord == 0) ? grad_B[k].positions[n].x() :
                                (coord == 1) ? grad_B[k].positions[n].y() :
                                               grad_B[k].positions[n].z();

            double ratio = (std::abs(numerical) > 1.0e-10) ? analytical / numerical : 0.0;
            std::string label = "B[" + std::to_string(k) + "][" + std::to_string(n) + "]."
                                + std::string(1, "xyz"[coord]);

            std::cout << std::setw(12) << label
                      << std::setw(14) << std::setprecision(6) << std::fixed << analytical
                      << std::setw(14) << std::setprecision(6) << std::fixed << numerical
                      << std::setw(14) << std::setprecision(4) << std::fixed << ratio
                      << std::endl;
         }
      }
   }
}
