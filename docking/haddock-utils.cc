/*
 * docking/haddock-utils.cc
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

#include <iostream>
#include <set>
#include <cmath>

#include <mmdb2/mmdb_tables.h>

#include "haddock-utils.hh"

// ---------------------------------------------------------------------------
// SASA filtering
// ---------------------------------------------------------------------------

// Per-residue relative SASA is computed by summing absolute SASA over
// all atoms of the residue and dividing by the maximum SASA for that
// residue type (from extended tripeptide Ala-X-Ala).
//
// Maximum SASA values for standard amino acids (from Tien et al. 2013,
// "Maximum allowed solvent accessibilities of residues in proteins",
// PLOS ONE 8:e80635).  Units: A^2.

static double max_sasa_for_residue(const std::string &resname) {

   // Empirical maximum SASA values for standard amino acids
   // (from extended Ala-X-Ala tripeptide, all-atom)
   if (resname == "ALA") return 129.0;
   if (resname == "ARG") return 274.0;
   if (resname == "ASN") return 195.0;
   if (resname == "ASP") return 193.0;
   if (resname == "CYS") return 167.0;
   if (resname == "GLN") return 225.0;
   if (resname == "GLU") return 223.0;
   if (resname == "GLY") return 104.0;
   if (resname == "HIS") return 224.0;
   if (resname == "ILE") return 197.0;
   if (resname == "LEU") return 201.0;
   if (resname == "LYS") return 236.0;
   if (resname == "MET") return 224.0;
   if (resname == "PHE") return 240.0;
   if (resname == "PRO") return 159.0;
   if (resname == "SER") return 155.0;
   if (resname == "THR") return 172.0;
   if (resname == "TRP") return 285.0;
   if (resname == "TYR") return 263.0;
   if (resname == "VAL") return 174.0;
   return 200.0; // default
}

// Simple per-atom SASA using the Shrake-Rupley dot method.
// Returns the total SASA per residue for the specified residues.
//
// This is a simplified implementation.  For production use, the
// CSASArea class from ccp4mg-utils could be used instead.

static std::vector<double>
compute_residue_sasa(mmdb::Manager *mol,
                     const std::vector<coot::residue_spec_t> &residues) {

   // Select all atoms for the neighbour search
   int all_sel = mol->NewSelection();
   mol->SelectAtoms(all_sel, 0, "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*", "*", "!H", "*");
   mmdb::Atom **all_atoms = nullptr;
   int n_all = 0;
   mol->GetSelIndex(all_sel, all_atoms, n_all);

   double probe_radius = 1.4; // water probe radius (A)
   int n_points = 92;         // number of test points per atom (Shrake-Rupley)

   // Generate Fibonacci sphere points
   std::vector<clipper::Coord_orth> sphere_points(n_points);
   double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
   for (int i=0; i<n_points; i++) {
      double theta = std::acos(1.0 - 2.0*(i + 0.5)/n_points);
      double phi = 2.0 * M_PI * i / golden_ratio;
      sphere_points[i] = clipper::Coord_orth(std::sin(theta)*std::cos(phi),
                                              std::sin(theta)*std::sin(phi),
                                              std::cos(theta));
   }

   std::vector<double> result(residues.size(), 0.0);

   for (unsigned int r=0; r<residues.size(); r++) {
      const coot::residue_spec_t &spec = residues[r];

      // Select atoms of this residue
      int res_sel = mol->NewSelection();
      mol->SelectAtoms(res_sel, 0, spec.chain_id.c_str(),
                       spec.res_no, spec.ins_code.c_str(),
                       spec.res_no, spec.ins_code.c_str(),
                       "*", "*", "!H", "*");
      mmdb::Atom **res_atoms = nullptr;
      int n_res = 0;
      mol->GetSelIndex(res_sel, res_atoms, n_res);

      double total_sasa = 0.0;

      for (int ia=0; ia<n_res; ia++) {
         mmdb::Atom *atom = res_atoms[ia];
         double atom_radius = mmdb::getVdWaalsRadius(atom->element) + probe_radius;
         double point_area = 4.0 * M_PI * atom_radius * atom_radius / n_points;

         int accessible = 0;
         for (int ip=0; ip<n_points; ip++) {
            double px = atom->x + atom_radius * sphere_points[ip].x();
            double py = atom->y + atom_radius * sphere_points[ip].y();
            double pz = atom->z + atom_radius * sphere_points[ip].z();

            bool buried = false;
            for (int j=0; j<n_all; j++) {
               if (all_atoms[j] == atom) continue;
               double nr = mmdb::getVdWaalsRadius(all_atoms[j]->element) + probe_radius;
               double dx = px - all_atoms[j]->x;
               double dy = py - all_atoms[j]->y;
               double dz = pz - all_atoms[j]->z;
               if (dx*dx + dy*dy + dz*dz < nr*nr) {
                  buried = true;
                  break;
               }
            }
            if (!buried) accessible++;
         }
         total_sasa += accessible * point_area;
      }

      result[r] = total_sasa;
      mol->DeleteSelection(res_sel);
   }

   mol->DeleteSelection(all_sel);
   return result;
}

std::vector<coot::residue_spec_t>
coot::haddock::filter_by_sasa(mmdb::Manager *mol,
                              const std::vector<coot::residue_spec_t> &candidates,
                              double sasa_threshold) {

   std::vector<double> sasa_values = compute_residue_sasa(mol, candidates);
   std::vector<coot::residue_spec_t> filtered;

   for (unsigned int i=0; i<candidates.size(); i++) {
      // Get residue name for max SASA lookup
      int res_sel = mol->NewSelection();
      mol->SelectAtoms(res_sel, 0, candidates[i].chain_id.c_str(),
                       candidates[i].res_no, candidates[i].ins_code.c_str(),
                       candidates[i].res_no, candidates[i].ins_code.c_str(),
                       "*", "*", "*", "*");
      mmdb::Atom **res_atoms = nullptr;
      int n_res = 0;
      mol->GetSelIndex(res_sel, res_atoms, n_res);

      std::string resname;
      if (n_res > 0) resname = res_atoms[0]->GetResName();
      mol->DeleteSelection(res_sel);

      double max_sasa = max_sasa_for_residue(resname);
      double relative_sasa = sasa_values[i] / max_sasa;

      if (relative_sasa > sasa_threshold) {
         filtered.push_back(candidates[i]);
      } else {
         std::cout << "HADDOCK: filtering out residue " << candidates[i]
                   << " (" << resname << ")"
                   << " relative SASA=" << relative_sasa
                   << " (threshold=" << sasa_threshold << ")" << std::endl;
      }
   }

   return filtered;
}

// ---------------------------------------------------------------------------
// Find passive residues
// ---------------------------------------------------------------------------

std::vector<coot::residue_spec_t>
coot::haddock::find_passive_residues(mmdb::Manager *mol,
                                     const std::vector<coot::residue_spec_t> &active_residues,
                                     double sasa_threshold,
                                     double neighbour_dist) {

   // Build a set of active residue specs for exclusion
   std::set<coot::residue_spec_t> active_set(active_residues.begin(),
                                              active_residues.end());

   // For each active residue, find neighbouring residues within
   // neighbour_dist of any of its atoms
   std::set<coot::residue_spec_t> candidate_passive;

   for (const auto &active_spec : active_residues) {
      // Select atoms of this active residue
      int act_sel = mol->NewSelection();
      mol->SelectAtoms(act_sel, 0, active_spec.chain_id.c_str(),
                       active_spec.res_no, active_spec.ins_code.c_str(),
                       active_spec.res_no, active_spec.ins_code.c_str(),
                       "*", "*", "!H", "*");
      mmdb::Atom **act_atoms = nullptr;
      int n_act = 0;
      mol->GetSelIndex(act_sel, act_atoms, n_act);

      // Find neighbouring atoms within neighbour_dist
      // Select all non-hydrogen atoms in the same chain
      int nbr_sel = mol->NewSelection();
      mol->SelectNeighbours(nbr_sel, mmdb::STYPE_ATOM,
                            act_atoms, n_act,
                            0.0, neighbour_dist,
                            mmdb::SKEY_NEW);
      mmdb::Atom **nbr_atoms = nullptr;
      int n_nbr = 0;
      mol->GetSelIndex(nbr_sel, nbr_atoms, n_nbr);

      for (int j=0; j<n_nbr; j++) {
         mmdb::Residue *res = nbr_atoms[j]->GetResidue();
         if (!res) continue;
         coot::residue_spec_t spec(res);
         // Exclude active residues themselves
         if (active_set.find(spec) == active_set.end()) {
            candidate_passive.insert(spec);
         }
      }

      mol->DeleteSelection(nbr_sel);
      mol->DeleteSelection(act_sel);
   }

   // Filter candidates by SASA
   std::vector<coot::residue_spec_t> candidates(candidate_passive.begin(),
                                                  candidate_passive.end());
   return filter_by_sasa(mol, candidates, sasa_threshold);
}

// ---------------------------------------------------------------------------
// Complete residue definition
// ---------------------------------------------------------------------------

coot::haddock::active_passive_residues_t
coot::haddock::define_residues(mmdb::Manager *mol,
                               const std::vector<coot::residue_spec_t> &active_candidates,
                               double sasa_threshold,
                               double neighbour_dist) {

   active_passive_residues_t result;

   // Filter active candidates by SASA
   result.active = filter_by_sasa(mol, active_candidates, sasa_threshold);

   std::cout << "HADDOCK: " << result.active.size() << " active residues (from "
             << active_candidates.size() << " candidates)" << std::endl;

   // Find passive residues
   result.passive = find_passive_residues(mol, result.active,
                                          sasa_threshold, neighbour_dist);

   std::cout << "HADDOCK: " << result.passive.size()
             << " passive residues identified" << std::endl;

   return result;
}
