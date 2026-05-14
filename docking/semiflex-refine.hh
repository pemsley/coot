/*
 * docking/semiflex-refine.hh
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

#ifndef SEMIFLEX_REFINE_HH
#define SEMIFLEX_REFINE_HH

#include <vector>
#include <random>
#include <set>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "haddock-types.hh"
#include "intermolecular-energy.hh"

namespace coot {
namespace haddock {

   // Semi-flexible simulated annealing refinement engine.
   //
   // Takes the output of rigid body docking (Stage 1) and refines
   // through three sequential SA phases:
   //
   //   Phase 1: rigid body SA (2000 -> 500 K)
   //   Phase 2: interface side-chain chi angles (1000 -> 50 K)
   //   Phase 3: interface backbone phi/psi + chi angles (1000 -> 50 K)
   //
   // Only protein B (the mobile protein) is made flexible.
   // Protein A remains rigid throughout.
   //
   class semiflex_refiner_t {
   public:
      semiflex_refiner_t(const molecule_data_t &data_A,
                         const molecule_data_t &data_B,
                         const docking_parameters_t &energy_params,
                         const semiflex_parameters_t &sf_params);

      // Refine a batch of rigid body results.
      // Thread-parallel over the input vector.
      // Returns refined results sorted by total energy.
      std::vector<semiflex_result_t>
      refine(const std::vector<docking_result_t> &rigid_results);

   private:

      molecule_data_t data_A_;
      molecule_data_t data_B_;
      docking_parameters_t energy_params_;
      semiflex_parameters_t sf_params_;

      // Spatial grid for protein A (fixed)
      spatial_grid_t grid_A_;

      // Centre of mass and centred atoms for protein B
      clipper::Coord_orth com_B_;
      std::vector<atom_data_t> centred_B_atoms_;

      // AIR residue atoms of protein B centred at origin
      std::vector<residue_atoms_t> centred_B_partner_atoms_;

      // --- Downstream atom tables ---
      // For each standard amino acid, for each chi angle, the list of
      // atom names that rotate (i.e. are "downstream" of the chi bond).
      // Built once in the constructor.
      //
      // downstream_atoms_["ARG"][0] = {"CG", "CD", "NE", "CZ", "NH1", "NH2"}
      //                                (chi1 rotates these)
      //
      std::map<std::string, std::vector<std::vector<std::string>>>
         downstream_atoms_;

      void build_downstream_atom_tables();

      // --- Interface detection ---
      // Given a rigid body pose, identify which residues of B are
      // within interface_cutoff of any atom of A.
      std::set<int> detect_interface_residue_indices(
         const std::vector<atom_data_t> &transformed_B) const;

      // --- Flexibility setup ---
      // Build residue_flex_info_t for each interface residue.
      // Maps chi atom names to atom indices in the B atom vector.
      std::vector<residue_flex_info_t>
      build_flex_info(const std::set<int> &interface_residue_indices,
                      const std::vector<atom_data_t> &B_atoms,
                      bool backbone_flexible) const;

      // --- Coordinate generation ---
      // Apply all torsion angle changes to a set of B atom positions.
      // Modifies atoms_p in place. Applies chi1 then chi2 etc.
      static void apply_torsion_state(
         std::vector<atom_data_t> *atoms_p,
         const std::vector<residue_flex_info_t> &flex);

      // Transform centred B atoms by a rigid body pose
      std::vector<atom_data_t>
      transform_B_atoms(const glm::quat &rotation,
                        const clipper::Coord_orth &translation) const;

      // Transform AIR residue atoms by a rigid body pose
      std::vector<residue_atoms_t>
      transform_B_air_atoms(const glm::quat &rotation,
                            const clipper::Coord_orth &translation) const;

      // --- Euler angle conversion (reuse from rigid_body_docker_t) ---
      static glm::quat euler_to_quat(double alpha, double beta, double gamma);
      static void quat_to_euler(const glm::quat &q,
                                double *alpha_p, double *beta_p,
                                double *gamma_p);

      // --- SA state ---
      // The Metropolis state is a flat double vector:
      //   [0..2] = ZYZ Euler angles (alpha, beta, gamma)
      //   [3..5] = translation (tx, ty, tz)
      //   [6..] = torsion angles for flexible residues
      //
      // The phase determines which DOF are perturbed.

      enum sa_phase_t { RIGID_SA, SIDECHAIN_SA, BACKBONE_SC_SA };

      // Run a single SA phase with Metropolis acceptance.
      // Updates state_vec in place.
      void run_sa_phase(std::vector<double> *state_vec_p,
                        const std::vector<residue_flex_info_t> &flex,
                        sa_phase_t phase,
                        double T_start, double T_end, int n_steps,
                        std::mt19937 &rng) const;

      // Evaluate energy for a given state vector
      intermolecular_energy_t
      evaluate_state(const std::vector<double> &state_vec,
                     const std::vector<residue_flex_info_t> &flex) const;

      // Pack/unpack between state vector and structured types
      std::vector<double>
      pack_state(const glm::quat &rotation,
                 const clipper::Coord_orth &translation,
                 const std::vector<residue_flex_info_t> &flex) const;

      void unpack_state(const std::vector<double> &state_vec,
                        glm::quat *rotation_p,
                        clipper::Coord_orth *translation_p,
                        std::vector<residue_flex_info_t> *flex_p) const;

      // Refine one result through all three phases
      semiflex_result_t refine_single(const docking_result_t &pose,
                                      std::mt19937 &rng) const;

      // Centre of mass computation
      static clipper::Coord_orth
      centre_of_mass(const std::vector<atom_data_t> &atoms);
   };

} // namespace haddock
} // namespace coot

#endif // SEMIFLEX_REFINE_HH
