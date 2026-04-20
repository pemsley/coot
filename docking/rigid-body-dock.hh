/*
 * docking/rigid-body-dock.hh
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

#ifndef RIGID_BODY_DOCK_HH
#define RIGID_BODY_DOCK_HH

#include <vector>
#include <random>

#include <gsl/gsl_multimin.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "haddock-types.hh"
#include "intermolecular-energy.hh"

namespace coot {
namespace haddock {

   // The rigid body docking engine.
   //
   // Performs HADDOCK Stage 1: random orientation generation followed by
   // rigid body energy minimisation driven by AIRs plus VDW/electrostatic
   // intermolecular energy.
   //
   class rigid_body_docker_t {
   public:
      rigid_body_docker_t(const molecule_data_t &mol_A,
                          const molecule_data_t &mol_B,
                          const docking_parameters_t &params);

      // Run the full rigid body docking protocol.
      // Returns the best n_keep results sorted by total energy.
      std::vector<docking_result_t> dock();

   private:

      molecule_data_t mol_A_;
      molecule_data_t mol_B_;
      docking_parameters_t params_;

      // Spatial grid for protein A (fixed) — built once, reused per evaluation
      spatial_grid_t grid_A_;

      // Centre of mass of each molecule (original coordinates)
      clipper::Coord_orth com_A_;
      clipper::Coord_orth com_B_;

      // Atoms of protein B centred at origin (for rotation)
      std::vector<atom_data_t> centred_B_atoms_;

      // AIR residue atoms of protein B centred at origin
      std::vector<residue_atoms_t> centred_B_partner_atoms_;

      // Compute centre of mass of atom positions
      static clipper::Coord_orth centre_of_mass(const std::vector<atom_data_t> &atoms);

      // Generate a uniform random quaternion
      static glm::quat random_quaternion(std::mt19937 &rng);

      // Apply a rotation (quaternion) and translation to a set of atom positions
      static std::vector<clipper::Coord_orth>
      transform_positions(const std::vector<atom_data_t> &atoms,
                          const glm::quat &rotation,
                          const clipper::Coord_orth &translation);

      // Apply transform to grouped residue atoms (for AIRs)
      static std::vector<residue_atoms_t>
      transform_residue_atoms(const std::vector<residue_atoms_t> &res_atoms,
                              const glm::quat &rotation,
                              const clipper::Coord_orth &translation);

      // Run a single docking trial (one random start)
      docking_result_t run_single_trial(std::mt19937 &rng);

      // Rigid body energy minimisation.
      // Optimises the rotation and translation of protein B to minimise
      // intermolecular energy against fixed protein A.
      //
      // rotation_p and translation_p are modified in place.
      //
      void rigid_body_minimise(glm::quat *rotation_p,
                               clipper::Coord_orth *translation_p,
                               bool rotation_only);

      // --- GSL minimiser interface ---
      //
      // The GSL minimiser operates on a gsl_vector of 6 parameters:
      //   [0..2] = ZYZ Euler angles (alpha, beta, gamma)
      //   [3..5] = translation (tx, ty, tz)
      //
      // When rotation_only is true, only [0..2] are varied.
      //
      // The params pointer is cast to a minimiser_context_t.

      struct minimiser_context_t {
         const rigid_body_docker_t *docker;
         bool rotation_only;
      };

      static double gsl_f(const gsl_vector *v, void *params);
      static void gsl_df(const gsl_vector *v, void *params, gsl_vector *df);
      static void gsl_fdf(const gsl_vector *v, void *params,
                          double *f, gsl_vector *df);

      // Convert Euler angles to quaternion (ZYZ convention)
      static glm::quat euler_to_quat(double alpha, double beta, double gamma);

      // Convert quaternion to Euler angles (ZYZ convention)
      static void quat_to_euler(const glm::quat &q,
                                double *alpha_p, double *beta_p, double *gamma_p);

      // Evaluate energy for a given pose (quaternion + translation)
      // of protein B against fixed protein A.
      intermolecular_energy_t evaluate_energy(const glm::quat &rotation,
                                              const clipper::Coord_orth &translation) const;

      // Evaluate energy and optionally compute analytical rigid body gradient.
      // If gradient_out is non-null, fills it with the 6-DOF gradient
      // [dE/dalpha, dE/dbeta, dE/dgamma, dE/dtx, dE/dty, dE/dtz].
      // alpha and beta are needed for the rotation axis computation.
      // If rotation_only, only the first 3 components are filled.
      intermolecular_energy_t
      evaluate_with_rb_gradient(const glm::quat &rotation,
                                const clipper::Coord_orth &translation,
                                double alpha, double beta,
                                bool rotation_only,
                                double *gradient_out) const;
   };

} // namespace haddock
} // namespace coot

#endif // RIGID_BODY_DOCK_HH
