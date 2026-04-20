/*
 * docking/intermolecular-energy.hh
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

#ifndef INTERMOLECULAR_ENERGY_HH
#define INTERMOLECULAR_ENERGY_HH

#include <vector>
#include <string>
#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>

#include "haddock-types.hh"
#include "air-restraint.hh"

namespace coot {
namespace haddock {

   // Pre-extracted atom data for efficient energy evaluation.
   // Extracted once from mmdb, then used repeatedly as the rigid body
   // pose is varied.
   //
   struct atom_data_t {
      clipper::Coord_orth position;
      double vdw_radius;      // van der Waals radius (A)
      double charge;          // partial charge (e)
      std::string atom_name;
      std::string residue_name;
      int residue_index;      // index into a parallel residue spec vector (-1 = unset)
      atom_data_t() : position(0,0,0), vdw_radius(1.7), charge(0),
                       residue_index(-1) {}
   };

   // Spatial grid (cell list) for fast neighbour lookup on fixed protein A.
   // Built once, reused for every energy evaluation against moving protein B.
   // Cell size = cutoff distance, so neighbours are found by checking the
   // 3×3×3 block of cells around a query point.
   //
   class spatial_grid_t {
   public:
      spatial_grid_t() : cell_size_(0), nx_(0), ny_(0), nz_(0) {}
      spatial_grid_t(const std::vector<atom_data_t> &atoms, double cutoff);

      // Return indices of atoms within cutoff of the query point.
      // (Candidates only — caller should still check exact distance.)
      void get_neighbours(const clipper::Coord_orth &query,
                          std::vector<unsigned int> *indices_out) const;
   private:
      double cell_size_;
      double origin_x_, origin_y_, origin_z_;
      int nx_, ny_, nz_;
      // cells_[ix + iy*nx_ + iz*nx_*ny_] = list of atom indices
      std::vector<std::vector<unsigned int>> cells_;
      int cell_index(int ix, int iy, int iz) const {
         return ix + iy * nx_ + iz * nx_ * ny_;
      }
   };

   // All the data needed to evaluate intermolecular energy.
   // Extracted from mmdb::Manager molecules before docking begins.
   //
   struct molecule_data_t {
      std::vector<atom_data_t> atoms;
      // Parallel vector: residue_specs[i] is the spec for residue_index i
      std::vector<coot::residue_spec_t> residue_specs;
      // AIR-related: atoms grouped by residue for active/passive residues
      std::vector<residue_atoms_t> active_residue_atoms;
      std::vector<residue_atoms_t> partner_residue_atoms; // active + passive
   };

   // Extract atom data from an mmdb::Manager for docking.
   // The selection should cover the atoms to include (typically all
   // non-hydrogen atoms, or all atoms).
   //
   molecule_data_t extract_molecule_data(mmdb::Manager *mol,
                                         const active_passive_residues_t &ap_residues,
                                         bool is_partner);

   // Compute the intermolecular energy between two sets of atoms.
   //
   // atoms_A: fixed protein atoms
   // atoms_B: moving protein atoms (already transformed)
   // active_A_atoms: active residue atom groups for protein A (for AIRs)
   // partner_B_atoms: active+passive residue atom groups for protein B
   // params: docking parameters (cutoffs, weights, etc.)
   //
   // Returns the energy breakdown.
   //
   // If grad_B is non-null, accumulates the gradient on each atom of
   // protein B (for projecting onto rigid body DOF).
   //
   intermolecular_energy_t
   compute_intermolecular_energy(const std::vector<atom_data_t> &atoms_A,
                                 const std::vector<atom_data_t> &atoms_B,
                                 const std::vector<residue_atoms_t> &active_A_atoms,
                                 const std::vector<residue_atoms_t> &partner_B_atoms,
                                 const docking_parameters_t &params,
                                 std::vector<clipper::Coord_orth> *grad_B,
                                 const spatial_grid_t *grid_A = nullptr);

} // namespace haddock
} // namespace coot

#endif // INTERMOLECULAR_ENERGY_HH
