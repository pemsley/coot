/*
 * docking/haddock-types.hh
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

#ifndef HADDOCK_TYPES_HH
#define HADDOCK_TYPES_HH

#include <vector>

#include <clipper/core/coords.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {
namespace haddock {

   // Active and passive residues for one protein.
   //
   // Active residues are those experimentally identified as being involved
   // in the interaction (e.g. significant chemical shift perturbation upon
   // complex formation, or mutagenesis data) and having >50% relative
   // solvent accessibility.
   //
   // Passive residues are surface neighbours of the active residues that
   // also have >50% relative solvent accessibility.  They do not have
   // direct AIRs to the partner protein but can satisfy the partner's
   // active restraints.
   //
   struct active_passive_residues_t {
      std::vector<coot::residue_spec_t> active;
      std::vector<coot::residue_spec_t> passive;
      // convenience: all active + passive
      std::vector<coot::residue_spec_t> all() const {
         std::vector<coot::residue_spec_t> r = active;
         r.insert(r.end(), passive.begin(), passive.end());
         return r;
      }
   };

   // An Ambiguous Interaction Restraint (AIR).
   //
   // Defined as an ambiguous intermolecular distance with a maximum value
   // of distance_cutoff (default 3.0 A) between any atom of an active
   // residue on one protein and any atom of the active+passive residues
   // on the partner protein.
   //
   // The effective distance uses 1/r^6 sum averaging:
   //
   //   d_eff = ( sum_m sum_k sum_n  1/d_{m,n}^6 )^(-1/6)
   //
   // where m indexes atoms of the active residue, k indexes partner
   // residues (active+passive), and n indexes atoms of residue k.
   //
   struct air_restraint_t {
      coot::residue_spec_t active_residue;  // on this protein
      // partner residues (active + passive) are stored externally
      double distance_cutoff;               // 3.0 A default
      air_restraint_t() : distance_cutoff(3.0) {}
      air_restraint_t(const coot::residue_spec_t &res, double d)
         : active_residue(res), distance_cutoff(d) {}
   };

   // Parameters for the docking protocol.
   //
   struct docking_parameters_t {
      double separation_distance = 150.0;   // initial placement distance (A)
      int n_trials = 1000;                  // number of random starting orientations
      int n_keep = 200;                     // best solutions to keep
      int n_orient_cycles = 4;              // alternating orientation optimisation cycles
      double vdw_cutoff = 8.5;             // A
      double elec_cutoff = 8.5;            // A
      double lj_epsilon = 1.0;             // kcal/mol
      double dielectric = 10.0;            // effective dielectric constant
      double air_weight = 1.0;
      double elec_weight = 1.0;
      double vdw_weight = 1.0;
      unsigned int n_threads = 0;           // 0 = auto-detect
   };

   // Result of a single rigid body docking trial.
   //
   struct docking_result_t {
      glm::quat rotation;                   // rotation applied to protein B
      clipper::Coord_orth translation;       // translation applied to protein B
      double e_vdw;
      double e_elec;
      double e_air;
      double e_total;                       // weighted sum
      docking_result_t() :
         rotation(glm::quat(1,0,0,0)),
         translation(0,0,0),
         e_vdw(0), e_elec(0), e_air(0), e_total(0) {}
   };

   // Intermolecular energy breakdown.
   //
   struct intermolecular_energy_t {
      double e_vdw;
      double e_elec;
      double e_air;
      double e_total;
      intermolecular_energy_t() : e_vdw(0), e_elec(0), e_air(0), e_total(0) {}
      intermolecular_energy_t(double vdw, double elec, double air, double total)
         : e_vdw(vdw), e_elec(elec), e_air(air), e_total(total) {}
   };

   // ---------------------------------------------------------------
   // Semi-flexible refinement types (Stage 2)
   // ---------------------------------------------------------------

   // Which atoms rotate when a torsion angle is changed.
   // All indices refer to the molecule's atom_data_t vector.
   //
   struct torsion_move_t {
      int axis_atom_1_idx;   // bond origin (closer to backbone)
      int axis_atom_2_idx;   // bond end (farther from backbone)
      std::vector<int> moving_atom_indices; // atoms beyond axis_atom_2
      double current_angle;  // radians
      torsion_move_t() : axis_atom_1_idx(-1), axis_atom_2_idx(-1),
                          current_angle(0) {}
   };

   // Per-residue flexibility information, built once per pose.
   //
   struct residue_flex_info_t {
      coot::residue_spec_t spec;
      std::string residue_name;
      std::vector<int> all_atom_indices;    // into molecule atom vector
      std::vector<torsion_move_t> chi_moves; // chi1..chi4
      torsion_move_t phi_move;              // only used in phase 3
      torsion_move_t psi_move;              // only used in phase 3
   };

   // Parameters for semi-flexible simulated annealing refinement.
   //
   struct semiflex_parameters_t {
      // Phase 1: rigid body SA
      double T_start_rigid = 2000.0;
      double T_end_rigid = 500.0;
      int n_steps_rigid = 500;
      // Phase 2: side-chain flexibility
      double T_start_sc = 1000.0;
      double T_end_sc = 50.0;
      int n_steps_sc = 1000;
      // Phase 3: backbone + side-chain flexibility
      double T_start_both = 1000.0;
      double T_end_both = 50.0;
      int n_steps_both = 1000;
      // Interface detection
      double interface_cutoff = 5.0;     // A
      // Move sizes
      double max_chi_step = 20.0;        // degrees
      double max_phi_psi_step = 10.0;    // degrees
      double max_rb_rot_step = 5.0;      // degrees
      double max_rb_trans_step = 0.5;    // A
      unsigned int n_threads = 0;        // 0 = auto-detect
   };

   // Result of semi-flexible refinement.
   //
   struct semiflex_result_t {
      docking_result_t pose;
      std::vector<residue_flex_info_t> flex_residues;
      double e_total;
      semiflex_result_t() : e_total(0) {}
   };

} // namespace haddock
} // namespace coot

#endif // HADDOCK_TYPES_HH
