/*
 * docking/semiflex-refine.cc
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
#include <thread>
#include <chrono>
#include <map>

#include <glm/gtc/matrix_transform.hpp>

#include "semiflex-refine.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "utils/split-indices.hh"

using namespace coot::haddock;

// ---------------------------------------------------------------------------
// Downstream atom tables for standard amino acids.
//
// For each residue type, for each chi angle (in order chi1..chiN),
// the list of atom names that rotate when that chi is changed.
// These are the atoms "downstream" of the chi bond axis
// (i.e. further from the backbone than atom 3 of the quad).
//
// Chi axis is atom2 -> atom3 of the quad (e.g. CA->CB for chi1).
// Moving atoms = everything beyond atom3 in the residue topology.
// ---------------------------------------------------------------------------

static std::map<std::string, std::vector<std::vector<std::string>>>
build_downstream_tables() {

   std::map<std::string, std::vector<std::vector<std::string>>> t;

   // ARG: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD, chi3 CB-CG-CD-NE, chi4 CG-CD-NE-CZ
   t["ARG"] = {
      {"CG", "CD", "NE", "CZ", "NH1", "NH2"},   // chi1: rotate CG and beyond
      {"CD", "NE", "CZ", "NH1", "NH2"},           // chi2: rotate CD and beyond
      {"NE", "CZ", "NH1", "NH2"},                  // chi3
      {"CZ", "NH1", "NH2"}                         // chi4
   };

   // ASP: chi1 N-CA-CB-CG, chi2 CA-CB-CG-OD1
   t["ASP"] = {
      {"CG", "OD1", "OD2"},
      {"OD1", "OD2"}
   };

   // ASN: chi1 N-CA-CB-CG, chi2 CA-CB-CG-OD1
   t["ASN"] = {
      {"CG", "OD1", "ND2"},
      {"OD1", "ND2"}
   };

   // CYS: chi1 N-CA-CB-SG
   t["CYS"] = {
      {"SG"}
   };

   // GLU: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD, chi3 CB-CG-CD-OE1
   t["GLU"] = {
      {"CG", "CD", "OE1", "OE2"},
      {"CD", "OE1", "OE2"},
      {"OE1", "OE2"}
   };

   // GLN: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD, chi3 CB-CG-CD-OE1
   t["GLN"] = {
      {"CG", "CD", "OE1", "NE2"},
      {"CD", "OE1", "NE2"},
      {"OE1", "NE2"}
   };

   // HIS: chi1 N-CA-CB-CG, chi2 CA-CB-CG-ND1
   t["HIS"] = {
      {"CG", "ND1", "CD2", "CE1", "NE2"},
      {"ND1", "CD2", "CE1", "NE2"}
   };

   // ILE: chi1 N-CA-CB-CG1, chi2 CA-CB-CG1-CD1
   // Note: CG2 is on the other branch, does NOT rotate with chi1
   // For chi1 of ILE, the axis is CA->CB, and CG1 + CD1 rotate,
   // but CG2 also rotates (it's attached to CB on the other side?).
   // Actually no — CG2 is directly bonded to CB. When we rotate
   // around CA-CB, everything beyond CB rotates, which includes
   // CG1, CG2, and CD1. But chi1 is defined as N-CA-CB-CG1, so
   // for the standard chi rotation, only the CG1 branch moves.
   // Wait — that's not right either. In standard chi definitions,
   // chi1 rotates around the CA-CB bond, and ALL atoms beyond CB
   // move together. The chi angle is measured to CG1 but CG2 also
   // rotates. Let me be precise: chi1 rotates around the CA-CB bond
   // axis. Everything beyond CB (CG1, CG2, CD1) all rotate.
   t["ILE"] = {
      {"CG1", "CG2", "CD1"},  // chi1: everything beyond CB
      {"CD1"}                   // chi2: only CD1 beyond CG1
   };

   // LEU: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD1
   t["LEU"] = {
      {"CG", "CD1", "CD2"},
      {"CD1", "CD2"}           // both CD1 and CD2 are beyond CG
   };

   // LYS: chi1..chi4
   t["LYS"] = {
      {"CG", "CD", "CE", "NZ"},
      {"CD", "CE", "NZ"},
      {"CE", "NZ"},
      {"NZ"}
   };

   // MET: chi1 N-CA-CB-CG, chi2 CA-CB-CG-SD, chi3 CB-CG-SD-CE
   t["MET"] = {
      {"CG", "SD", "CE"},
      {"SD", "CE"},
      {"CE"}
   };

   // PHE: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD1
   t["PHE"] = {
      {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
      {"CD1", "CD2", "CE1", "CE2", "CZ"}
   };

   // PRO: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD
   // (Pro is unusual — ring constraint; skip for flexibility)
   t["PRO"] = {
      {"CG", "CD"},
      {"CD"}
   };

   // SER: chi1 N-CA-CB-OG
   t["SER"] = {
      {"OG"}
   };

   // THR: chi1 N-CA-CB-OG1
   // CG2 is on the other branch from CB
   t["THR"] = {
      {"OG1", "CG2"}   // both beyond CB
   };

   // TRP: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD1
   t["TRP"] = {
      {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
      {"CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}
   };

   // TYR: chi1 N-CA-CB-CG, chi2 CA-CB-CG-CD1
   t["TYR"] = {
      {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
      {"CD1", "CD2", "CE1", "CE2", "CZ", "OH"}
   };

   // VAL: chi1 N-CA-CB-CG1
   // CG2 is attached to CB; both branch atoms rotate with chi1
   t["VAL"] = {
      {"CG1", "CG2"}
   };

   return t;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

semiflex_refiner_t::semiflex_refiner_t(const molecule_data_t &data_A,
                                       const molecule_data_t &data_B,
                                       const docking_parameters_t &energy_params,
                                       const semiflex_parameters_t &sf_params)
   : data_A_(data_A), data_B_(data_B),
     energy_params_(energy_params), sf_params_(sf_params) {

   // Build spatial grid for fixed protein A
   double cutoff = std::max(energy_params_.vdw_cutoff, energy_params_.elec_cutoff);
   grid_A_ = spatial_grid_t(data_A_.atoms, cutoff);

   // Centre of mass of protein B
   com_B_ = centre_of_mass(data_B_.atoms);

   // Centre B atoms at origin
   centred_B_atoms_ = data_B_.atoms;
   for (auto &at : centred_B_atoms_) {
      at.position = clipper::Coord_orth(at.position.x() - com_B_.x(),
                                        at.position.y() - com_B_.y(),
                                        at.position.z() - com_B_.z());
   }

   // Centre B partner (AIR) atoms at origin
   centred_B_partner_atoms_ = data_B_.partner_residue_atoms;
   for (auto &ra : centred_B_partner_atoms_) {
      for (auto &pos : ra.positions) {
         pos = clipper::Coord_orth(pos.x() - com_B_.x(),
                                   pos.y() - com_B_.y(),
                                   pos.z() - com_B_.z());
      }
   }

   // Build downstream atom tables
   downstream_atoms_ = build_downstream_tables();
}

// ---------------------------------------------------------------------------
// Centre of mass
// ---------------------------------------------------------------------------

clipper::Coord_orth
semiflex_refiner_t::centre_of_mass(const std::vector<atom_data_t> &atoms) {

   double sx = 0, sy = 0, sz = 0;
   for (const auto &at : atoms) {
      sx += at.position.x();
      sy += at.position.y();
      sz += at.position.z();
   }
   double n = static_cast<double>(atoms.size());
   return clipper::Coord_orth(sx/n, sy/n, sz/n);
}

// ---------------------------------------------------------------------------
// Euler angle / quaternion conversion (ZYZ convention)
// (Duplicated from rigid_body_docker_t to keep semiflex self-contained)
// ---------------------------------------------------------------------------

glm::quat
semiflex_refiner_t::euler_to_quat(double alpha, double beta, double gamma) {

   float fa = static_cast<float>(alpha);
   float fb = static_cast<float>(beta);
   float fg = static_cast<float>(gamma);
   glm::quat qza = glm::angleAxis(fa, glm::vec3(0, 0, 1));
   glm::quat qyb = glm::angleAxis(fb, glm::vec3(0, 1, 0));
   glm::quat qzg = glm::angleAxis(fg, glm::vec3(0, 0, 1));
   return qza * qyb * qzg;
}

void
semiflex_refiner_t::quat_to_euler(const glm::quat &q,
                                  double *alpha_p, double *beta_p,
                                  double *gamma_p) {

   glm::mat3 R = glm::mat3_cast(q);
   double cos_beta = static_cast<double>(R[2][2]);
   cos_beta = std::clamp(cos_beta, -1.0, 1.0);
   double beta = std::acos(cos_beta);

   double alpha, gamma;
   if (std::abs(std::sin(beta)) > 1e-6) {
      alpha = std::atan2(static_cast<double>(R[2][1]),
                         static_cast<double>(R[2][0]));
      gamma = std::atan2(static_cast<double>(R[1][2]),
                         static_cast<double>(-R[0][2]));
   } else {
      alpha = std::atan2(static_cast<double>(R[1][0]),
                         static_cast<double>(R[0][0]));
      gamma = 0.0;
   }

   *alpha_p = alpha;
   *beta_p = beta;
   *gamma_p = gamma;
}

// ---------------------------------------------------------------------------
// Transform centred B atoms by a rigid body pose
// ---------------------------------------------------------------------------

std::vector<atom_data_t>
semiflex_refiner_t::transform_B_atoms(const glm::quat &rotation,
                                      const clipper::Coord_orth &translation) const {

   glm::mat3 R = glm::mat3_cast(rotation);
   std::vector<atom_data_t> result = centred_B_atoms_;
   for (auto &at : result) {
      glm::vec3 v(static_cast<float>(at.position.x()),
                   static_cast<float>(at.position.y()),
                   static_cast<float>(at.position.z()));
      glm::vec3 rv = R * v;
      at.position = clipper::Coord_orth(rv.x + translation.x(),
                                         rv.y + translation.y(),
                                         rv.z + translation.z());
   }
   return result;
}

// ---------------------------------------------------------------------------
// Transform AIR residue atoms by a rigid body pose
// ---------------------------------------------------------------------------

std::vector<residue_atoms_t>
semiflex_refiner_t::transform_B_air_atoms(const glm::quat &rotation,
                                          const clipper::Coord_orth &translation) const {

   glm::mat3 R = glm::mat3_cast(rotation);
   std::vector<residue_atoms_t> result = centred_B_partner_atoms_;
   for (auto &ra : result) {
      for (auto &pos : ra.positions) {
         glm::vec3 v(static_cast<float>(pos.x()),
                      static_cast<float>(pos.y()),
                      static_cast<float>(pos.z()));
         glm::vec3 rv = R * v;
         pos = clipper::Coord_orth(rv.x + translation.x(),
                                    rv.y + translation.y(),
                                    rv.z + translation.z());
      }
   }
   return result;
}

// ---------------------------------------------------------------------------
// Interface detection
// ---------------------------------------------------------------------------
//
// Given transformed B atoms, find which residue indices of B have any
// atom within interface_cutoff of any atom of A (using the spatial grid).

std::set<int>
semiflex_refiner_t::detect_interface_residue_indices(
   const std::vector<atom_data_t> &transformed_B) const {

   std::set<int> interface_res_indices;
   double cutoff_sq = sf_params_.interface_cutoff * sf_params_.interface_cutoff;
   std::vector<unsigned int> nbr_indices;
   nbr_indices.reserve(64);

   for (unsigned int j=0; j<transformed_B.size(); j++) {
      const atom_data_t &bj = transformed_B[j];
      if (bj.residue_index < 0) continue;

      // Already found this residue?
      if (interface_res_indices.count(bj.residue_index)) continue;

      grid_A_.get_neighbours(bj.position, &nbr_indices);
      for (unsigned int idx : nbr_indices) {
         double dx = bj.position.x() - data_A_.atoms[idx].position.x();
         double dy = bj.position.y() - data_A_.atoms[idx].position.y();
         double dz = bj.position.z() - data_A_.atoms[idx].position.z();
         double d_sq = dx*dx + dy*dy + dz*dz;
         if (d_sq < cutoff_sq) {
            interface_res_indices.insert(bj.residue_index);
            break;
         }
      }
   }

   return interface_res_indices;
}

// ---------------------------------------------------------------------------
// Build flexibility info for interface residues
// ---------------------------------------------------------------------------
//
// For each interface residue, look up chi angle quads in the downstream
// atom tables, map atom names to indices in the B atom vector, and build
// torsion_move_t entries.

std::vector<residue_flex_info_t>
semiflex_refiner_t::build_flex_info(const std::set<int> &interface_residue_indices,
                                    const std::vector<atom_data_t> &B_atoms,
                                    bool backbone_flexible) const {

   std::vector<residue_flex_info_t> flex_vec;

   for (int res_idx : interface_residue_indices) {
      if (res_idx < 0 || res_idx >= static_cast<int>(data_B_.residue_specs.size()))
         continue;

      const coot::residue_spec_t &spec = data_B_.residue_specs[res_idx];

      // Collect all atom indices for this residue and build name->index map
      std::map<std::string, int> name_to_idx;
      std::vector<int> all_indices;

      for (unsigned int i=0; i<B_atoms.size(); i++) {
         if (B_atoms[i].residue_index == res_idx) {
            all_indices.push_back(static_cast<int>(i));
            // Strip whitespace from atom name for lookup
            std::string aname = B_atoms[i].atom_name;
            // Remove leading/trailing spaces
            size_t start = aname.find_first_not_of(' ');
            size_t end = aname.find_last_not_of(' ');
            if (start != std::string::npos)
               aname = aname.substr(start, end - start + 1);
            name_to_idx[aname] = static_cast<int>(i);
         }
      }

      if (all_indices.empty()) continue;

      const std::string &resname = B_atoms[all_indices[0]].residue_name;

      // Look up downstream atom table for this residue type
      auto dt_it = downstream_atoms_.find(resname);
      if (dt_it == downstream_atoms_.end()) continue; // Not a standard residue or GLY/ALA

      const auto &downstream_for_res = dt_it->second;

      residue_flex_info_t flex;
      flex.spec = spec;
      flex.residue_name = resname;
      flex.all_atom_indices = all_indices;

      // Chi angle quads: N-CA-CB-CG pattern
      // The chi axis is atoms[1]->atoms[2] of the quad (0-indexed).
      // For chi1: axis = CA -> CB
      // For chi2: axis = CB -> CG
      // etc.
      //
      // Chi atom quads (standard):
      //   chi1: N  CA CB  X   (axis CA->CB, X = CG/SG/OG/OG1/CG1)
      //   chi2: CA CB CG  Y   (axis CB->CG, Y varies)
      //   chi3: CB CG CD  Z   (axis CG->CD)
      //   chi4: CG CD NE/CE Z (axis CD->NE/CE)
      //
      // Axis atom 1 = quad[1], axis atom 2 = quad[2]
      // Moving atoms = downstream table entry for this chi

      // The chi quads are well-known. We use fixed patterns:
      // chi1 axis: CA->CB (or equivalent)
      // chi2 axis: CB->CG (or CB->CG1 for ILE)
      // etc.
      // Rather than pulling in chi_angles class, use the quad patterns directly.

      struct chi_quad_t {
         std::string a1, a2; // axis atoms (a1->a2 is axis)
      };

      // Build axis atoms from the downstream table structure.
      // chi1 axis is always CA -> CB (except for weird cases).
      // For subsequent chis, axis is: quad[1] -> quad[2].
      // Standard patterns per number of chis:
      //
      //   chi1: CA->CB
      //   chi2: CB->{CG, CG1} (ILE uses CG1)
      //   chi3: CG->CD  or  CG->SD (MET)
      //   chi4: CD->NE  or  CD->CE (LYS)

      // Hard-coded chi axis atoms per residue type per chi index
      // This matches the chi_angles setup_chi_atom_quads() definitions
      struct chi_axes_table_t {
         std::string axis1;
         std::string axis2;
      };

      // Build axes from quad definitions:
      // Quad is: atom0-atom1-atom2-atom3, axis = atom1->atom2
      static const std::map<std::string, std::vector<chi_axes_table_t>> chi_axes = {
         {"ARG", {{"CA","CB"}, {"CB","CG"}, {"CG","CD"}, {"CD","NE"}}},
         {"ASP", {{"CA","CB"}, {"CB","CG"}}},
         {"ASN", {{"CA","CB"}, {"CB","CG"}}},
         {"CYS", {{"CA","CB"}}},
         {"GLU", {{"CA","CB"}, {"CB","CG"}, {"CG","CD"}}},
         {"GLN", {{"CA","CB"}, {"CB","CG"}, {"CG","CD"}}},
         {"HIS", {{"CA","CB"}, {"CB","CG"}}},
         {"ILE", {{"CA","CB"}, {"CB","CG1"}}},
         {"LEU", {{"CA","CB"}, {"CB","CG"}}},
         {"LYS", {{"CA","CB"}, {"CB","CG"}, {"CG","CD"}, {"CD","CE"}}},
         {"MET", {{"CA","CB"}, {"CB","CG"}, {"CG","SD"}}},
         {"PHE", {{"CA","CB"}, {"CB","CG"}}},
         {"PRO", {{"CA","CB"}, {"CB","CG"}}},
         {"SER", {{"CA","CB"}}},
         {"THR", {{"CA","CB"}}},
         {"TRP", {{"CA","CB"}, {"CB","CG"}}},
         {"TYR", {{"CA","CB"}, {"CB","CG"}}},
         {"VAL", {{"CA","CB"}}}
      };

      auto axes_it = chi_axes.find(resname);
      if (axes_it == chi_axes.end()) continue;

      const auto &axes_for_res = axes_it->second;
      unsigned int n_chi = std::min(downstream_for_res.size(), axes_for_res.size());

      for (unsigned int c=0; c<n_chi; c++) {
         const std::string &ax1_name = axes_for_res[c].axis1;
         const std::string &ax2_name = axes_for_res[c].axis2;

         auto it1 = name_to_idx.find(ax1_name);
         auto it2 = name_to_idx.find(ax2_name);
         if (it1 == name_to_idx.end() || it2 == name_to_idx.end())
            continue; // Missing axis atoms

         torsion_move_t move;
         move.axis_atom_1_idx = it1->second;
         move.axis_atom_2_idx = it2->second;
         move.current_angle = 0.0; // deltas from initial conformation

         // Map downstream atom names to indices
         for (const std::string &aname : downstream_for_res[c]) {
            auto ait = name_to_idx.find(aname);
            if (ait != name_to_idx.end()) {
               move.moving_atom_indices.push_back(ait->second);
            }
         }

         if (!move.moving_atom_indices.empty()) {
            flex.chi_moves.push_back(move);
         }
      }

      // Backbone flexibility (phase 3 only)
      if (backbone_flexible) {
         // Phi: C(i-1) - N - CA - C     axis = N -> CA
         // Psi: N - CA - C - N(i+1)     axis = CA -> C
         //
         // For phi, the moving atoms are everything on the side-chain side
         // of the N-CA bond (i.e. CB and beyond + C, O).
         // For psi, the moving atoms are C and O of this residue.
         //
         // However, in the semi-flexible context, we only make small
         // perturbations and the "moving atoms" are the side-chain atoms
         // (which move with backbone changes anyway via the chi rotations).
         // We approximate: phi rotates everything beyond CA (CB side),
         // psi rotates C, O.

         auto n_it = name_to_idx.find("N");
         auto ca_it = name_to_idx.find("CA");
         auto c_it = name_to_idx.find("C");
         auto o_it = name_to_idx.find("O");

         if (n_it != name_to_idx.end() && ca_it != name_to_idx.end()) {
            // Phi: axis N -> CA
            torsion_move_t phi;
            phi.axis_atom_1_idx = n_it->second;
            phi.axis_atom_2_idx = ca_it->second;
            phi.current_angle = 0.0;
            // Moving atoms: everything except N (backbone C, O, side chain)
            for (int ai : all_indices) {
               std::string aname = B_atoms[ai].atom_name;
               size_t s = aname.find_first_not_of(' ');
               size_t e = aname.find_last_not_of(' ');
               if (s != std::string::npos) aname = aname.substr(s, e - s + 1);
               if (aname != "N" && aname != "CA") {
                  phi.moving_atom_indices.push_back(ai);
               }
            }
            if (!phi.moving_atom_indices.empty())
               flex.phi_move = phi;
         }

         if (ca_it != name_to_idx.end() && c_it != name_to_idx.end()) {
            // Psi: axis CA -> C
            torsion_move_t psi;
            psi.axis_atom_1_idx = ca_it->second;
            psi.axis_atom_2_idx = c_it->second;
            psi.current_angle = 0.0;
            // Moving atoms: C and O only (next residue's N would move
            // but we don't modify cross-residue atoms)
            if (o_it != name_to_idx.end()) {
               psi.moving_atom_indices.push_back(c_it->second);
               psi.moving_atom_indices.push_back(o_it->second);
            }
            if (!psi.moving_atom_indices.empty())
               flex.psi_move = psi;
         }
      }

      if (!flex.chi_moves.empty() ||
          (backbone_flexible && flex.phi_move.axis_atom_1_idx >= 0))
         flex_vec.push_back(flex);
   }

   return flex_vec;
}

// ---------------------------------------------------------------------------
// Apply torsion angle state to atom positions
// ---------------------------------------------------------------------------
//
// For each flexible residue, apply chi rotations sequentially
// (chi1 first, then chi2, etc.) using rotate_around_vector.
// Each chi rotation moves its downstream atoms relative to the
// current positions (which may already be modified by earlier chis
// in the same residue).

void
semiflex_refiner_t::apply_torsion_state(
   std::vector<atom_data_t> *atoms_p,
   const std::vector<residue_flex_info_t> &flex) {

   for (const auto &fi : flex) {

      // Apply phi if set (phase 3)
      if (fi.phi_move.axis_atom_1_idx >= 0 &&
          std::abs(fi.phi_move.current_angle) > 1e-10) {

         const auto &move = fi.phi_move;
         const clipper::Coord_orth &origin =
            (*atoms_p)[move.axis_atom_1_idx].position;
         clipper::Coord_orth axis(
            (*atoms_p)[move.axis_atom_2_idx].position.x() - origin.x(),
            (*atoms_p)[move.axis_atom_2_idx].position.y() - origin.y(),
            (*atoms_p)[move.axis_atom_2_idx].position.z() - origin.z());

         for (int mi : move.moving_atom_indices) {
            (*atoms_p)[mi].position =
               coot::util::rotate_around_vector(axis,
                  (*atoms_p)[mi].position, origin, move.current_angle);
         }
      }

      // Apply psi if set (phase 3)
      if (fi.psi_move.axis_atom_1_idx >= 0 &&
          std::abs(fi.psi_move.current_angle) > 1e-10) {

         const auto &move = fi.psi_move;
         const clipper::Coord_orth &origin =
            (*atoms_p)[move.axis_atom_1_idx].position;
         clipper::Coord_orth axis(
            (*atoms_p)[move.axis_atom_2_idx].position.x() - origin.x(),
            (*atoms_p)[move.axis_atom_2_idx].position.y() - origin.y(),
            (*atoms_p)[move.axis_atom_2_idx].position.z() - origin.z());

         for (int mi : move.moving_atom_indices) {
            (*atoms_p)[mi].position =
               coot::util::rotate_around_vector(axis,
                  (*atoms_p)[mi].position, origin, move.current_angle);
         }
      }

      // Apply chi angles in order (chi1 first, then chi2, etc.)
      for (const auto &chi : fi.chi_moves) {
         if (std::abs(chi.current_angle) < 1e-10) continue;

         const clipper::Coord_orth &origin =
            (*atoms_p)[chi.axis_atom_1_idx].position;
         clipper::Coord_orth axis(
            (*atoms_p)[chi.axis_atom_2_idx].position.x() - origin.x(),
            (*atoms_p)[chi.axis_atom_2_idx].position.y() - origin.y(),
            (*atoms_p)[chi.axis_atom_2_idx].position.z() - origin.z());

         for (int mi : chi.moving_atom_indices) {
            (*atoms_p)[mi].position =
               coot::util::rotate_around_vector(axis,
                  (*atoms_p)[mi].position, origin, chi.current_angle);
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Pack / unpack state vector
// ---------------------------------------------------------------------------
//
// State vector layout:
//   [0..2] = ZYZ Euler angles (alpha, beta, gamma) in radians
//   [3..5] = translation (tx, ty, tz) in Angstroms
//   [6..] = torsion angles: for each flex residue, first the chi angles
//           (chi1, chi2, ...), then phi, psi (if backbone_flexible)

std::vector<double>
semiflex_refiner_t::pack_state(const glm::quat &rotation,
                               const clipper::Coord_orth &translation,
                               const std::vector<residue_flex_info_t> &flex) const {

   double alpha, beta, gamma;
   quat_to_euler(rotation, &alpha, &beta, &gamma);

   std::vector<double> state;
   state.reserve(6 + flex.size() * 6); // generous estimate
   state.push_back(alpha);
   state.push_back(beta);
   state.push_back(gamma);
   state.push_back(translation.x());
   state.push_back(translation.y());
   state.push_back(translation.z());

   for (const auto &fi : flex) {
      for (const auto &chi : fi.chi_moves)
         state.push_back(chi.current_angle);
      state.push_back(fi.phi_move.current_angle);
      state.push_back(fi.psi_move.current_angle);
   }

   return state;
}

void
semiflex_refiner_t::unpack_state(const std::vector<double> &state_vec,
                                 glm::quat *rotation_p,
                                 clipper::Coord_orth *translation_p,
                                 std::vector<residue_flex_info_t> *flex_p) const {

   *rotation_p = euler_to_quat(state_vec[0], state_vec[1], state_vec[2]);
   *translation_p = clipper::Coord_orth(state_vec[3], state_vec[4], state_vec[5]);

   unsigned int offset = 6;
   for (auto &fi : *flex_p) {
      for (auto &chi : fi.chi_moves) {
         chi.current_angle = state_vec[offset++];
      }
      fi.phi_move.current_angle = state_vec[offset++];
      fi.psi_move.current_angle = state_vec[offset++];
   }
}

// ---------------------------------------------------------------------------
// Evaluate energy for a given state vector
// ---------------------------------------------------------------------------

intermolecular_energy_t
semiflex_refiner_t::evaluate_state(const std::vector<double> &state_vec,
                                   const std::vector<residue_flex_info_t> &flex) const {

   // Unpack rigid body DOF
   glm::quat rotation = euler_to_quat(state_vec[0], state_vec[1], state_vec[2]);
   clipper::Coord_orth translation(state_vec[3], state_vec[4], state_vec[5]);

   // Transform centred B atoms by rigid body pose
   std::vector<atom_data_t> B_atoms = transform_B_atoms(rotation, translation);

   // Apply torsion angles (need a mutable copy of flex to unpack into)
   std::vector<residue_flex_info_t> flex_copy = flex;
   unsigned int offset = 6;
   for (auto &fi : flex_copy) {
      for (auto &chi : fi.chi_moves)
         chi.current_angle = state_vec[offset++];
      fi.phi_move.current_angle = state_vec[offset++];
      fi.psi_move.current_angle = state_vec[offset++];
   }
   apply_torsion_state(&B_atoms, flex_copy);

   // Transform AIR atoms too
   std::vector<residue_atoms_t> B_air = transform_B_air_atoms(rotation, translation);
   // Note: AIR atoms are not flexed (they are a separate copy).
   // This is an approximation — in a full implementation we'd update
   // the AIR atoms with the torsion changes too.

   return compute_intermolecular_energy(data_A_.atoms, B_atoms,
                                        data_A_.active_residue_atoms,
                                        B_air,
                                        energy_params_,
                                        nullptr,  // no gradient
                                        &grid_A_);
}

// ---------------------------------------------------------------------------
// Simulated annealing phase
// ---------------------------------------------------------------------------
//
// Metropolis Monte Carlo with geometric temperature schedule:
//   T(step) = T_start * (T_end/T_start)^(step/(n_steps-1))
//
// In each step, one DOF is randomly perturbed within its step size.
// The perturbation is accepted with probability min(1, exp(-dE/T)).

void
semiflex_refiner_t::run_sa_phase(std::vector<double> *state_vec_p,
                                 const std::vector<residue_flex_info_t> &flex,
                                 sa_phase_t phase,
                                 double T_start, double T_end, int n_steps,
                                 std::mt19937 &rng) const {

   const double deg_to_rad = M_PI / 180.0;
   const double kB = 0.001987; // Boltzmann constant in kcal/(mol·K)

   // Determine which DOF indices are active for this phase
   // state_vec layout: [0..2]=euler, [3..5]=translation, [6..]=torsions
   std::vector<int> active_dof;

   if (phase == RIGID_SA) {
      // Rigid body: euler angles and translation
      for (int i=0; i<6; i++) active_dof.push_back(i);
   } else {
      // Torsion phases: chi angles, and optionally phi/psi
      unsigned int offset = 6;
      for (const auto &fi : flex) {
         // Chi angles
         for (unsigned int c=0; c<fi.chi_moves.size(); c++) {
            active_dof.push_back(static_cast<int>(offset));
            offset++;
         }
         // phi
         int phi_offset = static_cast<int>(offset);
         offset++;
         // psi
         int psi_offset = static_cast<int>(offset);
         offset++;

         if (phase == BACKBONE_SC_SA) {
            if (fi.phi_move.axis_atom_1_idx >= 0)
               active_dof.push_back(phi_offset);
            if (fi.psi_move.axis_atom_1_idx >= 0)
               active_dof.push_back(psi_offset);
         }
      }
   }

   if (active_dof.empty()) return;

   // Step sizes for each DOF
   // [0..2] = rotation (radians), [3..5] = translation (A), [6..] = torsions
   std::vector<double> step_sizes(state_vec_p->size(), 0.0);
   step_sizes[0] = sf_params_.max_rb_rot_step * deg_to_rad;
   step_sizes[1] = sf_params_.max_rb_rot_step * deg_to_rad;
   step_sizes[2] = sf_params_.max_rb_rot_step * deg_to_rad;
   step_sizes[3] = sf_params_.max_rb_trans_step;
   step_sizes[4] = sf_params_.max_rb_trans_step;
   step_sizes[5] = sf_params_.max_rb_trans_step;

   // Torsion step sizes
   unsigned int offset = 6;
   for (const auto &fi : flex) {
      for (unsigned int c=0; c<fi.chi_moves.size(); c++) {
         step_sizes[offset] = sf_params_.max_chi_step * deg_to_rad;
         offset++;
      }
      // phi
      step_sizes[offset] = sf_params_.max_phi_psi_step * deg_to_rad;
      offset++;
      // psi
      step_sizes[offset] = sf_params_.max_phi_psi_step * deg_to_rad;
      offset++;
   }

   // Current energy
   double current_E = evaluate_state(*state_vec_p, flex).e_total;

   std::uniform_int_distribution<int> dof_dist(0, static_cast<int>(active_dof.size()) - 1);
   std::uniform_real_distribution<double> uniform_01(0.0, 1.0);
   std::uniform_real_distribution<double> uniform_pm(-1.0, 1.0);

   double log_ratio = (T_end > 0 && T_start > 0) ? std::log(T_end / T_start) : 0.0;

   int n_accepted = 0;

   for (int step=0; step<n_steps; step++) {

      // Temperature: geometric schedule
      double frac = (n_steps > 1) ? static_cast<double>(step) / (n_steps - 1) : 1.0;
      double T = T_start * std::exp(log_ratio * frac);

      // Pick a random active DOF to perturb
      int dof_idx = active_dof[dof_dist(rng)];
      double delta = uniform_pm(rng) * step_sizes[dof_idx];

      // Make trial move
      double old_val = (*state_vec_p)[dof_idx];
      (*state_vec_p)[dof_idx] = old_val + delta;

      // Evaluate trial energy
      double trial_E = evaluate_state(*state_vec_p, flex).e_total;

      // Metropolis criterion
      double dE = trial_E - current_E;
      bool accept = false;
      if (dE <= 0) {
         accept = true;
      } else if (T > 1e-10) {
         double prob = std::exp(-dE / (kB * T));
         accept = (uniform_01(rng) < prob);
      }

      if (accept) {
         current_E = trial_E;
         n_accepted++;
      } else {
         // Reject: restore old value
         (*state_vec_p)[dof_idx] = old_val;
      }
   }

   if (false) { // debug output
      std::string phase_name;
      if (phase == RIGID_SA) phase_name = "rigid";
      else if (phase == SIDECHAIN_SA) phase_name = "sidechain";
      else phase_name = "backbone+sc";
      std::cout << "    SA " << phase_name << ": accepted "
                << n_accepted << "/" << n_steps
                << " (" << std::setprecision(1) << std::fixed
                << 100.0 * n_accepted / n_steps << "%)"
                << " E=" << std::setprecision(2) << current_E << std::endl;
   }
}

// ---------------------------------------------------------------------------
// Refine a single docking result through all three SA phases
// ---------------------------------------------------------------------------

semiflex_result_t
semiflex_refiner_t::refine_single(const docking_result_t &pose,
                                  std::mt19937 &rng) const {

   // Step 1: Transform B atoms by the input rigid body pose
   std::vector<atom_data_t> transformed_B = transform_B_atoms(pose.rotation,
                                                               pose.translation);

   // Step 2: Detect interface residues
   std::set<int> interface_indices = detect_interface_residue_indices(transformed_B);

   // Step 3: Build flex info (side-chain only initially)
   std::vector<residue_flex_info_t> flex_sc =
      build_flex_info(interface_indices, transformed_B, false);

   // Step 4: Build flex info with backbone (for phase 3)
   std::vector<residue_flex_info_t> flex_both =
      build_flex_info(interface_indices, transformed_B, true);

   // Pack initial state
   std::vector<double> state_vec = pack_state(pose.rotation, pose.translation, flex_sc);

   // Phase 1: Rigid body SA
   run_sa_phase(&state_vec, flex_sc, RIGID_SA,
                sf_params_.T_start_rigid, sf_params_.T_end_rigid,
                sf_params_.n_steps_rigid, rng);

   // Phase 2: Side-chain SA
   run_sa_phase(&state_vec, flex_sc, SIDECHAIN_SA,
                sf_params_.T_start_sc, sf_params_.T_end_sc,
                sf_params_.n_steps_sc, rng);

   // Phase 3: Backbone + side-chain SA
   // Need to repack state to include phi/psi DOF.
   // Unpack current state into rigid body + sc torsions
   glm::quat rotation;
   clipper::Coord_orth translation;
   unpack_state(state_vec, &rotation, &translation, &flex_sc);

   // Transfer chi angles from flex_sc to flex_both
   // (they share the same residue order and chi count)
   for (unsigned int r=0; r<flex_both.size() && r<flex_sc.size(); r++) {
      for (unsigned int c=0; c<flex_both[r].chi_moves.size() &&
           c<flex_sc[r].chi_moves.size(); c++) {
         flex_both[r].chi_moves[c].current_angle =
            flex_sc[r].chi_moves[c].current_angle;
      }
   }

   // Repack with backbone DOF included
   state_vec = pack_state(rotation, translation, flex_both);

   run_sa_phase(&state_vec, flex_both, BACKBONE_SC_SA,
                sf_params_.T_start_both, sf_params_.T_end_both,
                sf_params_.n_steps_both, rng);

   // Unpack final state
   unpack_state(state_vec, &rotation, &translation, &flex_both);

   // Evaluate final energy
   intermolecular_energy_t final_E = evaluate_state(state_vec, flex_both);

   // Build result
   semiflex_result_t result;
   result.pose.rotation = rotation;
   result.pose.translation = translation;
   result.pose.e_vdw = final_E.e_vdw;
   result.pose.e_elec = final_E.e_elec;
   result.pose.e_air = final_E.e_air;
   result.pose.e_total = final_E.e_total;
   result.flex_residues = flex_both;
   result.e_total = final_E.e_total;

   return result;
}

// ---------------------------------------------------------------------------
// Main entry point: refine a batch of rigid body results
// ---------------------------------------------------------------------------

std::vector<semiflex_result_t>
semiflex_refiner_t::refine(const std::vector<docking_result_t> &rigid_results) {

   unsigned int n_threads = sf_params_.n_threads;
   if (n_threads == 0) {
      n_threads = std::thread::hardware_concurrency();
      if (n_threads == 0) n_threads = 4;
   }
   unsigned int n_input = static_cast<unsigned int>(rigid_results.size());
   if (n_threads > n_input)
      n_threads = n_input;

   auto t_start = std::chrono::high_resolution_clock::now();

   std::cout << "HADDOCK semi-flexible refinement: "
             << n_input << " structures on "
             << n_threads << " threads" << std::endl;

   auto ranges = coot::atom_index_ranges(n_input, n_threads);

   std::vector<std::vector<semiflex_result_t>> thread_results(ranges.size());
   std::vector<std::thread> threads;

   for (unsigned int t=0; t<ranges.size(); t++) {
      threads.emplace_back([&, t]() {
         std::mt19937 rng(123 + t * 777 + ranges[t].first);
         auto &local_results = thread_results[t];
         local_results.reserve(ranges[t].second - ranges[t].first);

         for (unsigned int i=ranges[t].first; i<ranges[t].second; i++) {
            semiflex_result_t result = refine_single(rigid_results[i], rng);
            local_results.push_back(result);

            std::cout << "  Refined " << i+1 << "/" << n_input
                      << "  E_in=" << std::setprecision(1) << std::fixed
                      << rigid_results[i].e_total
                      << " -> E_out=" << result.e_total
                      << " (" << result.flex_residues.size()
                      << " flex res)" << std::endl;
         }
      });
   }

   for (auto &th : threads) th.join();

   // Merge results
   std::vector<semiflex_result_t> all_results;
   all_results.reserve(n_input);
   for (auto &tr : thread_results)
      all_results.insert(all_results.end(), tr.begin(), tr.end());

   // Sort by total energy
   std::sort(all_results.begin(), all_results.end(),
             [](const semiflex_result_t &a, const semiflex_result_t &b) {
                return a.e_total < b.e_total;
             });

   auto t_end = std::chrono::high_resolution_clock::now();
   double elapsed = std::chrono::duration<double>(t_end - t_start).count();
   std::cout << "HADDOCK semi-flexible refinement complete ("
             << std::setprecision(1) << std::fixed << elapsed << " s, "
             << std::setprecision(1) << elapsed / n_input << " s/structure)."
             << std::endl;
   if (!all_results.empty()) {
      std::cout << "  Best energy: " << std::setprecision(2) << std::fixed
                << all_results.front().e_total << std::endl;
      std::cout << "  Worst energy: " << all_results.back().e_total << std::endl;
   }

   return all_results;
}
