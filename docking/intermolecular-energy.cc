/*
 * docking/intermolecular-energy.cc
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>

#include <mmdb2/mmdb_tables.h>

#include "intermolecular-energy.hh"
#include "MoleculesToTriangles/CXXSurface/CXXChargeTable.h"

// ---------------------------------------------------------------------------
// Extract atom data from mmdb
// ---------------------------------------------------------------------------

coot::haddock::molecule_data_t
coot::haddock::extract_molecule_data(mmdb::Manager *mol,
                                     const active_passive_residues_t &ap_residues,
                                     bool is_partner) {

   molecule_data_t data;

   // Extract all non-hydrogen atoms
   int sel_hnd = mol->NewSelection();
   mol->SelectAtoms(sel_hnd, 0, "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*", "*", "!H", "*");

   mmdb::Atom **atoms = nullptr;
   int n_atoms = 0;
   mol->GetSelIndex(sel_hnd, atoms, n_atoms);

   CXXChargeTable charge_table;

   // Build residue index: map each unique (chain, resno, ins_code)
   // to a sequential index, storing the spec in residue_specs.
   std::map<std::string, int> res_key_to_index;

   data.atoms.resize(n_atoms);
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = atoms[i];
      data.atoms[i].position = clipper::Coord_orth(at->x, at->y, at->z);
      data.atoms[i].vdw_radius = mmdb::getVdWaalsRadius(at->element);
      data.atoms[i].charge = charge_table.getCharge(
         std::string(at->GetResName()),
         std::string(at->name));
      data.atoms[i].atom_name = at->name;
      data.atoms[i].residue_name = at->GetResName();

      // Assign residue index
      mmdb::Residue *res = at->GetResidue();
      if (res) {
         std::string key = std::string(res->GetChainID()) + "/"
                         + std::to_string(res->GetSeqNum()) + "/"
                         + std::string(res->GetInsCode());
         auto it = res_key_to_index.find(key);
         if (it == res_key_to_index.end()) {
            int idx = static_cast<int>(data.residue_specs.size());
            res_key_to_index[key] = idx;
            data.residue_specs.push_back(coot::residue_spec_t(res));
            data.atoms[i].residue_index = idx;
         } else {
            data.atoms[i].residue_index = it->second;
         }
      }
   }

   mol->DeleteSelection(sel_hnd);

   // Extract atom positions grouped by AIR residues.
   //
   // For active residues (used as the "source" of AIRs)
   for (unsigned int i=0; i<ap_residues.active.size(); i++) {
      const coot::residue_spec_t &spec = ap_residues.active[i];
      int res_sel = mol->NewSelection();
      mol->SelectAtoms(res_sel, 0, spec.chain_id.c_str(),
                       spec.res_no, spec.ins_code.c_str(),
                       spec.res_no, spec.ins_code.c_str(),
                       "*", "*", "!H", "*");
      mmdb::Atom **res_atoms = nullptr;
      int n_res_atoms = 0;
      mol->GetSelIndex(res_sel, res_atoms, n_res_atoms);

      residue_atoms_t ra;
      ra.positions.resize(n_res_atoms);
      for (int j=0; j<n_res_atoms; j++) {
         ra.positions[j] = clipper::Coord_orth(res_atoms[j]->x,
                                               res_atoms[j]->y,
                                               res_atoms[j]->z);
      }
      data.active_residue_atoms.push_back(ra);
      mol->DeleteSelection(res_sel);
   }

   // For partner residues (active + passive, used as the "target" of AIRs)
   if (is_partner) {
      std::vector<coot::residue_spec_t> all = ap_residues.all();
      for (unsigned int i=0; i<all.size(); i++) {
         const coot::residue_spec_t &spec = all[i];
         int res_sel = mol->NewSelection();
         mol->SelectAtoms(res_sel, 0, spec.chain_id.c_str(),
                          spec.res_no, spec.ins_code.c_str(),
                          spec.res_no, spec.ins_code.c_str(),
                          "*", "*", "!H", "*");
         mmdb::Atom **res_atoms = nullptr;
         int n_res_atoms = 0;
         mol->GetSelIndex(res_sel, res_atoms, n_res_atoms);

         residue_atoms_t ra;
         ra.positions.resize(n_res_atoms);
         for (int j=0; j<n_res_atoms; j++) {
            ra.positions[j] = clipper::Coord_orth(res_atoms[j]->x,
                                                  res_atoms[j]->y,
                                                  res_atoms[j]->z);
         }
         data.partner_residue_atoms.push_back(ra);
         mol->DeleteSelection(res_sel);
      }
   }

   return data;
}

// ---------------------------------------------------------------------------
// Spatial grid (cell list) for fast neighbour lookup
// ---------------------------------------------------------------------------

coot::haddock::spatial_grid_t::spatial_grid_t(const std::vector<atom_data_t> &atoms,
                                              double cutoff) {

   cell_size_ = cutoff;
   if (atoms.empty()) { nx_ = ny_ = nz_ = 0; return; }

   // Find bounding box with padding
   double xmin = atoms[0].position.x(), xmax = xmin;
   double ymin = atoms[0].position.y(), ymax = ymin;
   double zmin = atoms[0].position.z(), zmax = zmin;
   for (const auto &at : atoms) {
      double x = at.position.x(), y = at.position.y(), z = at.position.z();
      if (x < xmin) xmin = x; if (x > xmax) xmax = x;
      if (y < ymin) ymin = y; if (y > ymax) ymax = y;
      if (z < zmin) zmin = z; if (z > zmax) zmax = z;
   }

   // Pad by one cell on each side so query points slightly outside still work
   origin_x_ = xmin - cell_size_;
   origin_y_ = ymin - cell_size_;
   origin_z_ = zmin - cell_size_;
   nx_ = static_cast<int>((xmax - xmin) / cell_size_) + 3;
   ny_ = static_cast<int>((ymax - ymin) / cell_size_) + 3;
   nz_ = static_cast<int>((zmax - zmin) / cell_size_) + 3;

   cells_.resize(nx_ * ny_ * nz_);

   for (unsigned int i=0; i<atoms.size(); i++) {
      int ix = static_cast<int>((atoms[i].position.x() - origin_x_) / cell_size_);
      int iy = static_cast<int>((atoms[i].position.y() - origin_y_) / cell_size_);
      int iz = static_cast<int>((atoms[i].position.z() - origin_z_) / cell_size_);
      ix = std::clamp(ix, 0, nx_ - 1);
      iy = std::clamp(iy, 0, ny_ - 1);
      iz = std::clamp(iz, 0, nz_ - 1);
      cells_[cell_index(ix, iy, iz)].push_back(i);
   }
}

void
coot::haddock::spatial_grid_t::get_neighbours(const clipper::Coord_orth &query,
                                              std::vector<unsigned int> *indices_out) const {

   indices_out->clear();
   if (nx_ == 0) return;

   int cx = static_cast<int>((query.x() - origin_x_) / cell_size_);
   int cy = static_cast<int>((query.y() - origin_y_) / cell_size_);
   int cz = static_cast<int>((query.z() - origin_z_) / cell_size_);

   // Check 3x3x3 block of cells
   for (int dz=-1; dz<=1; dz++) {
      int iz = cz + dz;
      if (iz < 0 || iz >= nz_) continue;
      for (int dy=-1; dy<=1; dy++) {
         int iy = cy + dy;
         if (iy < 0 || iy >= ny_) continue;
         for (int dx=-1; dx<=1; dx++) {
            int ix = cx + dx;
            if (ix < 0 || ix >= nx_) continue;
            const auto &cell = cells_[cell_index(ix, iy, iz)];
            indices_out->insert(indices_out->end(), cell.begin(), cell.end());
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Intermolecular energy
// ---------------------------------------------------------------------------

// The intermolecular energy has three terms:
//
// E = w_vdw * E_vdw + w_elec * E_elec + w_air * E_AIR
//
// E_vdw:  Lennard-Jones 12-6 between all inter-protein atom pairs
//         within vdw_cutoff distance.
//
// E_elec: Coulomb q_i*q_j / (epsilon * r_ij) between all inter-protein
//         atom pairs within elec_cutoff distance.  Only atoms with
//         non-zero partial charges contribute.
//
// E_AIR:  Sum of ambiguous interaction restraint penalties.
//
// If grid_A is provided, it is used for fast neighbour lookup on protein A.
// Otherwise a brute-force O(N_A * N_B) loop is used.

// Inline helper to evaluate VDW + Coulomb for one (i, j) pair
static inline void
eval_pair(const coot::haddock::atom_data_t &ai,
          const coot::haddock::atom_data_t &aj,
          unsigned int j,
          double vdw_cutoff_sq, double elec_cutoff_sq,
          double coulomb_constant,
          const coot::haddock::docking_parameters_t &params,
          double *e_vdw_p, double *e_elec_p,
          std::vector<clipper::Coord_orth> *grad_B) {

   double dx = ai.position.x() - aj.position.x();
   double dy = ai.position.y() - aj.position.y();
   double dz = ai.position.z() - aj.position.z();
   double dist_sq = dx*dx + dy*dy + dz*dz;

   double cutoff_sq = std::max(vdw_cutoff_sq, elec_cutoff_sq);
   if (dist_sq >= cutoff_sq) return;

   // Van der Waals (Lennard-Jones 12-6)
   if (dist_sq < vdw_cutoff_sq) {
      double r_min = ai.vdw_radius + aj.vdw_radius;
      double d2 = dist_sq;
      if (d2 < 1.0) d2 = 1.0;
      double r_min_sq = r_min * r_min;
      double alpha_sq = r_min_sq / d2;
      double alpha_6 = alpha_sq * alpha_sq * alpha_sq;
      double alpha_12 = alpha_6 * alpha_6;
      *e_vdw_p += params.lj_epsilon * (alpha_12 - 2.0 * alpha_6);

      if (grad_B) {
         double dist = std::sqrt(d2);
         double dv_dr = 12.0 * params.lj_epsilon * alpha_6 * (1.0 - alpha_6) / dist;
         double inv_r = 1.0 / dist;
         double gx = dv_dr * (-dx) * inv_r;
         double gy = dv_dr * (-dy) * inv_r;
         double gz = dv_dr * (-dz) * inv_r;
         (*grad_B)[j] = clipper::Coord_orth(
            (*grad_B)[j].x() + params.vdw_weight * gx,
            (*grad_B)[j].y() + params.vdw_weight * gy,
            (*grad_B)[j].z() + params.vdw_weight * gz);
      }
   }

   // Electrostatics (Coulomb)
   if (dist_sq < elec_cutoff_sq) {
      double qi = ai.charge;
      double qj = aj.charge;
      if (std::abs(qi) > 1.0e-6 && std::abs(qj) > 1.0e-6) {
         double d2 = dist_sq;
         if (d2 < 1.0) d2 = 1.0;
         double dist = std::sqrt(d2);
         *e_elec_p += coulomb_constant * qi * qj / (params.dielectric * dist);

         if (grad_B) {
            double dv_dr = -coulomb_constant * qi * qj / (params.dielectric * d2);
            double inv_r = 1.0 / dist;
            double gx = dv_dr * (-dx) * inv_r;
            double gy = dv_dr * (-dy) * inv_r;
            double gz = dv_dr * (-dz) * inv_r;
            (*grad_B)[j] = clipper::Coord_orth(
               (*grad_B)[j].x() + params.elec_weight * gx,
               (*grad_B)[j].y() + params.elec_weight * gy,
               (*grad_B)[j].z() + params.elec_weight * gz);
         }
      }
   }
}

coot::haddock::intermolecular_energy_t
coot::haddock::compute_intermolecular_energy(const std::vector<atom_data_t> &atoms_A,
                                             const std::vector<atom_data_t> &atoms_B,
                                             const std::vector<residue_atoms_t> &active_A_atoms,
                                             const std::vector<residue_atoms_t> &partner_B_atoms,
                                             const docking_parameters_t &params,
                                             std::vector<clipper::Coord_orth> *grad_B,
                                             const spatial_grid_t *grid_A) {

   double e_vdw = 0.0;
   double e_elec = 0.0;

   double vdw_cutoff_sq = params.vdw_cutoff * params.vdw_cutoff;
   double elec_cutoff_sq = params.elec_cutoff * params.elec_cutoff;
   const double coulomb_constant = 332.0637;

   if (grid_A) {
      // Grid-accelerated path: for each B atom, find nearby A atoms
      std::vector<unsigned int> nbr_indices;
      nbr_indices.reserve(64);

      for (unsigned int j=0; j<atoms_B.size(); j++) {
         grid_A->get_neighbours(atoms_B[j].position, &nbr_indices);
         for (unsigned int idx : nbr_indices) {
            eval_pair(atoms_A[idx], atoms_B[j], j,
                      vdw_cutoff_sq, elec_cutoff_sq, coulomb_constant,
                      params, &e_vdw, &e_elec, grad_B);
         }
      }
   } else {
      // Brute-force O(N_A * N_B) path
      for (unsigned int i=0; i<atoms_A.size(); i++) {
         for (unsigned int j=0; j<atoms_B.size(); j++) {
            eval_pair(atoms_A[i], atoms_B[j], j,
                      vdw_cutoff_sq, elec_cutoff_sq, coulomb_constant,
                      params, &e_vdw, &e_elec, grad_B);
         }
      }
   }

   // AIR energy
   double e_air = air_energy_and_gradients(active_A_atoms, partner_B_atoms,
                                           3.0, nullptr, nullptr);

   double e_total = params.vdw_weight * e_vdw
                  + params.elec_weight * e_elec
                  + params.air_weight * e_air;

   return intermolecular_energy_t(e_vdw, e_elec, e_air, e_total);
}
