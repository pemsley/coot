/*
 * docking/air-restraint.hh
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

#ifndef AIR_RESTRAINT_HH
#define AIR_RESTRAINT_HH

#include <vector>
#include <clipper/core/coords.h>

namespace coot {
namespace haddock {

   // Atom positions for one residue, extracted from mmdb for efficiency.
   // We pre-extract coordinates so the inner loops don't touch mmdb.
   //
   struct residue_atoms_t {
      std::vector<clipper::Coord_orth> positions;
   };

   // Compute the AIR effective distance between one active residue
   // (on protein A) and the set of active+passive residues (on protein B).
   //
   // d_eff = ( sum_m sum_k sum_n  1/d_{m,n}^6 )^(-1/6)
   //
   // Returns d_eff in Angstroms.
   //
   double air_effective_distance(const residue_atoms_t &active_residue_A,
                                 const std::vector<residue_atoms_t> &partner_residues_B);

   // Compute the AIR penalty for a single restraint.
   //
   // Flat-bottom harmonic: zero when d_eff <= distance_cutoff,
   // (d_eff - distance_cutoff)^2 when d_eff > distance_cutoff.
   //
   double air_penalty(double d_eff, double distance_cutoff);

   // Compute the total AIR energy and gradients on all atoms.
   //
   // active_residues_A: pre-extracted atom positions for each active residue on protein A.
   // partner_residues_B: pre-extracted atom positions for all active+passive residues on protein B.
   // distance_cutoff: the AIR distance cutoff (default 3.0 A).
   //
   // Returns the total AIR energy.
   //
   // If grad_A is non-null, accumulates gradients on the atoms of the
   // active residues of protein A (same indexing as active_residues_A).
   // If grad_B is non-null, accumulates gradients on the atoms of the
   // partner residues of protein B.
   //
   double air_energy_and_gradients(const std::vector<residue_atoms_t> &active_residues_A,
                                   const std::vector<residue_atoms_t> &partner_residues_B,
                                   double distance_cutoff,
                                   std::vector<residue_atoms_t> *grad_A,
                                   std::vector<residue_atoms_t> *grad_B);

   // Numerical gradient check for debugging.
   // Perturbs each coordinate by micro_step and compares with analytical gradient.
   // Prints comparison to stdout.
   //
   void air_numerical_gradient_check(const std::vector<residue_atoms_t> &active_residues_A,
                                     const std::vector<residue_atoms_t> &partner_residues_B,
                                     double distance_cutoff,
                                     double micro_step = 0.0001);

} // namespace haddock
} // namespace coot

#endif // AIR_RESTRAINT_HH
