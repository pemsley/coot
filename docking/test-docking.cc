/*
 * docking/test-docking.cc
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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sys/stat.h>

#include <clipper/core/coords.h>
#include <mmdb2/mmdb_manager.h>
#include <mmdb2/mmdb_tables.h>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "air-restraint.hh"
#include "intermolecular-energy.hh"
#include "rigid-body-dock.hh"
#include "haddock-types.hh"
#include "haddock-utils.hh"
#include "semiflex-refine.hh"
#include "coot-utils/coot-coord-utils.hh"

using namespace coot::haddock;

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

static int n_tests_passed = 0;
static int n_tests_failed = 0;

static void report(const std::string &test_name, bool passed) {
   if (passed) {
      std::cout << "  PASS: " << test_name << std::endl;
      n_tests_passed++;
   } else {
      std::cout << "  FAIL: " << test_name << std::endl;
      n_tests_failed++;
   }
}

// ---------------------------------------------------------------------------
// Test 1: AIR effective distance
// ---------------------------------------------------------------------------

static void test_air_effective_distance() {

   std::cout << "\n--- Test: AIR effective distance ---\n";

   residue_atoms_t res_A;
   res_A.positions.push_back(clipper::Coord_orth(0, 0, 0));

   residue_atoms_t res_B;
   res_B.positions.push_back(clipper::Coord_orth(5, 0, 0));

   std::vector<residue_atoms_t> partner_B;
   partner_B.push_back(res_B);

   double d_eff = air_effective_distance(res_A, partner_B);
   std::cout << "  d_eff = " << d_eff << " (expected 5.0)" << std::endl;
   report("single pair d_eff", std::abs(d_eff - 5.0) < 0.01);

   // Add a second atom at 3 A
   double expected = std::pow(1.0/std::pow(5.0,6) + 1.0/std::pow(3.0,6), -1.0/6.0);
   partner_B[0].positions.push_back(clipper::Coord_orth(3, 0, 0));

   d_eff = air_effective_distance(res_A, partner_B);
   std::cout << "  d_eff = " << d_eff << " (expected " << expected << ")" << std::endl;
   report("two-atom d_eff", std::abs(d_eff - expected) < 0.01);
   report("d_eff dominated by closer atom", d_eff < 5.0 && d_eff > 2.5);
}

// ---------------------------------------------------------------------------
// Test 2: AIR penalty
// ---------------------------------------------------------------------------

static void test_air_penalty() {

   std::cout << "\n--- Test: AIR penalty ---\n";

   report("no penalty within cutoff", air_penalty(2.0, 3.0) == 0.0);
   report("no penalty at cutoff", air_penalty(3.0, 3.0) == 0.0);
   report("correct penalty at 5 A", std::abs(air_penalty(5.0, 3.0) - 4.0) < 0.01);
   report("correct penalty at 10 A", std::abs(air_penalty(10.0, 3.0) - 49.0) < 0.01);
}

// ---------------------------------------------------------------------------
// Test 3: AIR gradient numerical check
// ---------------------------------------------------------------------------

static void test_air_gradient() {

   std::cout << "\n--- Test: AIR gradient (numerical check) ---\n";

   residue_atoms_t res_A;
   res_A.positions.push_back(clipper::Coord_orth(0, 0, 0));
   res_A.positions.push_back(clipper::Coord_orth(1.5, 0, 0));
   res_A.positions.push_back(clipper::Coord_orth(0, 1.5, 0));

   residue_atoms_t res_B1;
   res_B1.positions.push_back(clipper::Coord_orth(6, 0, 0));
   res_B1.positions.push_back(clipper::Coord_orth(7, 1, 0));

   residue_atoms_t res_B2;
   res_B2.positions.push_back(clipper::Coord_orth(6, 2, 1));
   res_B2.positions.push_back(clipper::Coord_orth(7, 0, 1));

   std::vector<residue_atoms_t> active_A;
   active_A.push_back(res_A);

   std::vector<residue_atoms_t> partner_B;
   partner_B.push_back(res_B1);
   partner_B.push_back(res_B2);

   double d_eff = air_effective_distance(res_A, partner_B);
   std::cout << "  d_eff = " << d_eff << " A (cutoff 3.0 A)\n";
   report("AIR is violated (d_eff > 3)", d_eff > 3.0);

   std::cout << "\n  Numerical vs analytical gradient comparison:\n";
   air_numerical_gradient_check(active_A, partner_B, 3.0, 0.0001);

   // Check gradient direction
   std::vector<residue_atoms_t> grad_A(1);
   grad_A[0].positions.resize(3, clipper::Coord_orth(0,0,0));
   std::vector<residue_atoms_t> grad_B(2);
   grad_B[0].positions.resize(2, clipper::Coord_orth(0,0,0));
   grad_B[1].positions.resize(2, clipper::Coord_orth(0,0,0));

   air_energy_and_gradients(active_A, partner_B, 3.0, &grad_A, &grad_B);

   // A is at x~0, B is at x~6.  Moving A toward B (+x) decreases penalty,
   // so dE/dx_A should be negative.
   report("grad_A x-component is negative (force toward B)",
          grad_A[0].positions[0].x() < 0);
   report("grad_B x-component is positive (force toward A)",
          grad_B[0].positions[0].x() > 0);
}

// ---------------------------------------------------------------------------
// Test 4: Dock chain E of 9V3F back to the ABC trimer
// ---------------------------------------------------------------------------
//
// Strategy:
// 1. Read pdb9v3f.ent
// 2. Identify interface residues between the receptor (A+B+C) and the
//    ligand (E) — in practice the binding cleft is formed by A and C
//    but we don't know that beforehand
// 3. Extract the receptor chains and chain E into separate mmdb::Managers
// 4. Use interface residues as active residues
// 5. Run rigid body docking
// 6. Check that the lowest energy solution places E near its original position

// Extract a set of chains from a molecule, deleting all others
static mmdb::Manager *extract_chains(mmdb::Manager *mol,
                                     const std::vector<std::string> &keep_chains) {

   mmdb::Manager *new_mol = new mmdb::Manager();
   new_mol->Copy(mol, mmdb::MMDBFCM_All);

   std::set<std::string> keep_set(keep_chains.begin(), keep_chains.end());

   mmdb::Model *model = new_mol->GetModel(1);
   if (!model) return new_mol;

   int n_chains = model->GetNumberOfChains();
   for (int i=n_chains-1; i>=0; i--) {
      mmdb::Chain *chain = model->GetChain(i);
      if (chain) {
         std::string cid(chain->GetChainID());
         if (keep_set.find(cid) == keep_set.end()) {
            model->DeleteChain(i);
         }
      }
   }
   new_mol->FinishStructEdit();
   return new_mol;
}

// Find interface residues between a set of receptor chains and a ligand
// chain by contact distance.  Returns (receptor_interface, ligand_interface).
static std::pair<std::vector<coot::residue_spec_t>, std::vector<coot::residue_spec_t>>
find_interface(mmdb::Manager *mol,
               const std::vector<std::string> &receptor_chains,
               const std::string &ligand_chain,
               double contact_dist) {

   std::set<coot::residue_spec_t> iface_rec, iface_lig;

   // Select atoms from all receptor chains
   int sel_rec = mol->NewSelection();
   for (const auto &cid : receptor_chains) {
      mol->SelectAtoms(sel_rec, 0, cid.c_str(),
                       mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                       "*", "*", "!H", "*",
                       mmdb::SKEY_OR);
   }

   int sel_lig = mol->NewSelection();
   mol->SelectAtoms(sel_lig, 0, ligand_chain.c_str(),
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "!H", "*");

   mmdb::Atom **atoms_rec = nullptr;
   mmdb::Atom **atoms_lig = nullptr;
   int n_rec = 0, n_lig = 0;
   mol->GetSelIndex(sel_rec, atoms_rec, n_rec);
   mol->GetSelIndex(sel_lig, atoms_lig, n_lig);

   double contact_dist_sq = contact_dist * contact_dist;

   for (int i=0; i<n_rec; i++) {
      for (int j=0; j<n_lig; j++) {
         double dx = atoms_rec[i]->x - atoms_lig[j]->x;
         double dy = atoms_rec[i]->y - atoms_lig[j]->y;
         double dz = atoms_rec[i]->z - atoms_lig[j]->z;
         if (dx*dx + dy*dy + dz*dz < contact_dist_sq) {
            mmdb::Residue *r_rec = atoms_rec[i]->GetResidue();
            mmdb::Residue *r_lig = atoms_lig[j]->GetResidue();
            if (r_rec) iface_rec.insert(coot::residue_spec_t(r_rec));
            if (r_lig) iface_lig.insert(coot::residue_spec_t(r_lig));
         }
      }
   }

   mol->DeleteSelection(sel_rec);
   mol->DeleteSelection(sel_lig);

   return std::make_pair(
      std::vector<coot::residue_spec_t>(iface_rec.begin(), iface_rec.end()),
      std::vector<coot::residue_spec_t>(iface_lig.begin(), iface_lig.end()));
}

// Apply a docking result's transform to an mmdb::Manager.
// The docker centres protein B by subtracting com_B, then applies
// R * centred_pos + translation.  So the overall transform on an
// original atom position is:
//   new_pos = R * (old_pos - com_B) + translation
//
static void apply_dock_transform(mmdb::Manager *mol,
                                 const docking_result_t &result,
                                 const clipper::Coord_orth &com_B) {

   glm::mat3 R = glm::mat3_cast(result.rotation);

   int sel = mol->NewSelection();
   mol->SelectAtoms(sel, 0, "*",
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "*", "*");
   mmdb::Atom **atoms = nullptr;
   int n_atoms = 0;
   mol->GetSelIndex(sel, atoms, n_atoms);

   for (int i=0; i<n_atoms; i++) {
      glm::vec3 p(static_cast<float>(atoms[i]->x - com_B.x()),
                   static_cast<float>(atoms[i]->y - com_B.y()),
                   static_cast<float>(atoms[i]->z - com_B.z()));
      glm::vec3 rp = R * p;
      atoms[i]->x = rp.x + result.translation.x();
      atoms[i]->y = rp.y + result.translation.y();
      atoms[i]->z = rp.z + result.translation.z();
   }

   mol->DeleteSelection(sel);
}

// Write a complex PDB: chain A (fixed) + docked chain E
static void write_docked_pdb(mmdb::Manager *mol_A,
                             mmdb::Manager *mol_E_template,
                             const docking_result_t &result,
                             const clipper::Coord_orth &com_E,
                             const std::string &filename) {

   // Make a copy of chain E and apply the transform
   mmdb::Manager *docked_E = new mmdb::Manager();
   docked_E->Copy(mol_E_template, mmdb::MMDBFCM_All);
   apply_dock_transform(docked_E, result, com_E);

   // Build a combined molecule: copy A, then add chain E
   mmdb::Manager *combined = new mmdb::Manager();
   combined->Copy(mol_A, mmdb::MMDBFCM_All);

   mmdb::Model *model = combined->GetModel(1);
   mmdb::Model *e_model = docked_E->GetModel(1);
   if (model && e_model) {
      int n_chains_e = e_model->GetNumberOfChains();
      for (int i=0; i<n_chains_e; i++) {
         mmdb::Chain *chain = e_model->GetChain(i);
         if (chain) {
            mmdb::Chain *new_chain = new mmdb::Chain();
            new_chain->Copy(chain);
            model->AddChain(new_chain);
         }
      }
   }
   combined->FinishStructEdit();
   combined->WritePDBASCII(filename.c_str());

   std::cout << "  Wrote " << filename << std::endl;

   delete docked_E;
   delete combined;
}

// Write CSV with all docking results for plotting
static void write_results_csv(const std::vector<docking_result_t> &results,
                              const clipper::Coord_orth &original_com_E,
                              const std::string &filename) {

   std::ofstream csv(filename);
   csv << "rank,e_total,e_vdw,e_elec,e_air,com_dist" << std::endl;

   for (unsigned int i=0; i<results.size(); i++) {
      const docking_result_t &r = results[i];
      double d = clipper::Coord_orth::length(r.translation, original_com_E);
      csv << i << ","
          << std::setprecision(2) << std::fixed
          << r.e_total << ","
          << r.e_vdw << ","
          << r.e_elec << ","
          << r.e_air << ","
          << d << std::endl;
   }
   std::cout << "  Wrote " << filename << " (" << results.size()
             << " results)" << std::endl;
}

static void test_dock_9v3f_chain_E(const std::string &pdb_file) {

   std::cout << "\n--- Test: Dock chain E of 9V3F onto ABC trimer ---\n";

   // Read the PDB
   mmdb::Manager *mol = new mmdb::Manager();
   int rc = mol->ReadCoorFile(pdb_file.c_str());
   if (rc != 0) {
      std::cout << "  ERROR: could not read " << pdb_file << std::endl;
      report("read PDB file", false);
      delete mol;
      return;
   }
   report("read PDB file", true);

   // The receptor is the ABC trimer.  The ligand is chain E (nanobody).
   // The binding site is formed by chains A and C, but we treat the
   // whole trimer as the receptor since in a real scenario we wouldn't
   // know which chains contribute.
   std::vector<std::string> receptor_chains = {"A", "B", "C"};
   std::string ligand_chain = "E";

   // Find interface residues between receptor (A+B+C) and ligand (E)
   auto [iface_rec, iface_lig] = find_interface(mol, receptor_chains,
                                                 ligand_chain, 5.0);
   std::cout << "  Interface residues: receptor=" << iface_rec.size()
             << " ligand=" << iface_lig.size() << std::endl;
   report("found interface residues",
          iface_rec.size() > 5 && iface_lig.size() > 5);

   // Show which receptor chains contribute to the interface
   std::set<std::string> interface_chains;
   for (const auto &spec : iface_rec)
      interface_chains.insert(spec.chain_id);
   std::cout << "  Receptor interface chains:";
   for (const auto &cid : interface_chains)
      std::cout << " " << cid;
   std::cout << std::endl;

   // Record original centre of mass of chain E for later comparison
   int sel_E_orig = mol->NewSelection();
   mol->SelectAtoms(sel_E_orig, 0, ligand_chain.c_str(),
                    mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "!H", "*");
   mmdb::Atom **orig_E_atoms = nullptr;
   int n_orig_E = 0;
   mol->GetSelIndex(sel_E_orig, orig_E_atoms, n_orig_E);
   double ox = 0, oy = 0, oz = 0;
   for (int i=0; i<n_orig_E; i++) {
      ox += orig_E_atoms[i]->x;
      oy += orig_E_atoms[i]->y;
      oz += orig_E_atoms[i]->z;
   }
   clipper::Coord_orth original_com_E(ox/n_orig_E, oy/n_orig_E, oz/n_orig_E);
   mol->DeleteSelection(sel_E_orig);

   std::cout << "  Original chain E COM: ("
             << original_com_E.x() << ", "
             << original_com_E.y() << ", "
             << original_com_E.z() << ")\n";

   // Extract the receptor (A+B+C) and ligand (E) into separate molecules
   mmdb::Manager *mol_rec = extract_chains(mol, receptor_chains);
   mmdb::Manager *mol_lig = extract_chains(mol, {ligand_chain});

   // Build active_passive_residues for each.
   // Use interface residues directly as active (skip SASA filtering for now
   // since these are known interface residues from the complex structure).
   active_passive_residues_t ap_rec;
   ap_rec.active = iface_rec;
   active_passive_residues_t ap_lig;
   ap_lig.active = iface_lig;

   std::cout << "  Active residues: receptor=" << ap_rec.active.size()
             << " ligand=" << ap_lig.active.size() << std::endl;

   // Extract molecule data for docking.
   // Receptor is fixed, ligand is moved.
   // For AIRs: receptor's active residues are the "source", ligand's
   // active+passive are the "partner".
   molecule_data_t data_rec = extract_molecule_data(mol_rec, ap_rec, false);
   molecule_data_t data_lig = extract_molecule_data(mol_lig, ap_lig, true);

   std::cout << "  Atom counts: receptor=" << data_rec.atoms.size()
             << " ligand=" << data_lig.atoms.size() << std::endl;
   std::cout << "  AIR residues: receptor active="
             << data_rec.active_residue_atoms.size()
             << " ligand partner="
             << data_lig.partner_residue_atoms.size() << std::endl;

   // Run rigid body docking with modest parameters for testing
   docking_parameters_t params;
   params.n_trials = 200;
   params.n_keep = 20;
   params.n_orient_cycles = 4;
   params.separation_distance = -1.0;  // auto-compute from protein sizes
   params.n_threads = 4;

   rigid_body_docker_t docker(data_rec, data_lig, params);
   std::vector<docking_result_t> results = docker.dock();

   report("got docking results", !results.empty());

   if (!results.empty()) {
      // Print top 5 results
      int n_show = std::min(5, static_cast<int>(results.size()));
      for (int i=0; i<n_show; i++) {
         const docking_result_t &r = results[i];
         clipper::Coord_orth dcom = r.translation;
         double d = clipper::Coord_orth::length(dcom, original_com_E);
         std::cout << "  Result " << i << ": E=" << std::setprecision(1) << std::fixed
                   << r.e_total
                   << " (VDW=" << r.e_vdw
                   << " Elec=" << r.e_elec
                   << " AIR=" << r.e_air
                   << ") dist=" << d << " A" << std::endl;
      }

      const docking_result_t &best = results[0];
      std::cout << "  Best energy: " << std::setprecision(2) << std::fixed
                << best.e_total
                << " (VDW=" << best.e_vdw
                << " Elec=" << best.e_elec
                << " AIR=" << best.e_air << ")\n";

      // Compute where the docked chain E COM ends up
      // The docker centres B at origin, rotates, then translates.
      // So the docked COM = R * (original_com_E - com_E_centred) + translation
      // But we need to figure out what com_E the docker used...
      // Actually, the docker centres B by subtracting com_B_, so:
      //   docked_pos = R * (original_pos - com_B) + translation
      // and docked_com = R * (0, 0, 0) + translation = translation
      // (since centring moves COM to origin)

      clipper::Coord_orth docked_com = best.translation;
      double dist = clipper::Coord_orth::length(docked_com, original_com_E);
      std::cout << "  Docked chain E COM: ("
                << docked_com.x() << ", "
                << docked_com.y() << ", "
                << docked_com.z() << ")\n";
      std::cout << "  Distance from original COM: " << dist << " A\n";

      // The docked position should be reasonably close to the original.
      // With 200 trials and rigid body only, within ~20 A of the
      // correct position is good.  Within 15 A is encouraging —
      // subsequent semi-flexible refinement (stage 2) would improve this.
      report("docked COM within 20 A of original", dist < 20.0);
      report("docked COM within 15 A of original", dist < 15.0);

      // The best energy should be negative (attractive interactions)
      report("best total energy is negative", best.e_total < 0);

      // Write output files for visualisation
      std::cout << "\n  Writing output files...\n";

      // Write the receptor (A+B+C) alone
      mol_rec->WritePDBASCII("dock-receptor.pdb");
      std::cout << "  Wrote dock-receptor.pdb (A+B+C, original coords)"
                << std::endl;

      // Write the native complex (all chains) for comparison
      mol->WritePDBASCII("dock-native-complex.pdb");
      std::cout << "  Wrote dock-native-complex.pdb (original)" << std::endl;

      // Write the top 5 docking results: receptor + docked ligand
      int n_write = std::min(5, static_cast<int>(results.size()));
      for (int i=0; i<n_write; i++) {
         std::string fname = "dock-result-" + std::to_string(i) + ".pdb";
         write_docked_pdb(mol_rec, mol_lig, results[i],
                          original_com_E, fname);
      }

      // Write CSV of all kept results
      write_results_csv(results, original_com_E, "dock-results.csv");
   }

   delete mol_rec;
   delete mol_lig;
   delete mol;
}

// ---------------------------------------------------------------------------
// Test 5: E2A-HPr docking with real NMR chemical shift perturbation data
// ---------------------------------------------------------------------------
//
// This is the classic HADDOCK benchmark case from Dominguez et al.,
// JACS 125, 1731-1737 (2003).  The active residues were determined
// experimentally by NMR titration:
//
//   E2A (EIN domain):  van Nuland et al., J. Mol. Biol. 246, 180-193 (1995)
//   HPr:               Zhou & Bhatt, Biochemistry 37, 6269-6277 (1998)
//
// Structures:
//   E2A unbound: PDB 1F3G (X-ray, chain A, residues 19-168)
//   HPr unbound: PDB 1HDN (NMR ensemble, chain A, residues 1-85)
//   Reference:   PDB 1GGR (E2A=chain A, HPr=chain B, residues 301-385)
//
// The AIR restraint definitions come from the HADDOCK3 tutorial:
//   haddocking/haddock3/examples/docking-protein-protein/data/e2a-hpr_air.tbl
//

// Compute centre of mass from an mmdb::Manager (all atoms or selected chain)
static clipper::Coord_orth compute_com(mmdb::Manager *mol,
                                       const std::string &chain_id = "") {

   int sel = mol->NewSelection();
   if (chain_id.empty()) {
      mol->SelectAtoms(sel, 1, "*",
                       mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                       "*", "*", "!H", "*");
   } else {
      mol->SelectAtoms(sel, 1, chain_id.c_str(),
                       mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                       "*", "*", "!H", "*");
   }

   mmdb::Atom **atoms = nullptr;
   int n_atoms = 0;
   mol->GetSelIndex(sel, atoms, n_atoms);

   double sx = 0, sy = 0, sz = 0;
   for (int i=0; i<n_atoms; i++) {
      sx += atoms[i]->x;
      sy += atoms[i]->y;
      sz += atoms[i]->z;
   }
   mol->DeleteSelection(sel);

   if (n_atoms == 0) return clipper::Coord_orth(0, 0, 0);
   return clipper::Coord_orth(sx/n_atoms, sy/n_atoms, sz/n_atoms);
}

// Keep only model 1 of an NMR ensemble
static void keep_model_1(mmdb::Manager *mol) {

   int n_models = mol->GetNumberOfModels();
   if (n_models <= 1) return;

   for (int i=n_models; i>=2; i--) {
      mol->DeleteModel(i);
   }
   mol->FinishStructEdit();
}

static void test_dock_e2a_hpr(const std::string &pdb_e2a,
                               const std::string &pdb_hpr,
                               const std::string &pdb_ref) {

   std::cout << "\n--- Test: E2A-HPr docking (real NMR CSP data) ---\n";

   // Read the three structures
   mmdb::Manager *mol_e2a = new mmdb::Manager();
   mmdb::Manager *mol_hpr = new mmdb::Manager();
   mmdb::Manager *mol_ref = new mmdb::Manager();

   int rc_e2a = mol_e2a->ReadCoorFile(pdb_e2a.c_str());
   int rc_hpr = mol_hpr->ReadCoorFile(pdb_hpr.c_str());
   int rc_ref = mol_ref->ReadCoorFile(pdb_ref.c_str());

   if (rc_e2a != 0 || rc_hpr != 0 || rc_ref != 0) {
      std::cout << "  ERROR: could not read PDB files"
                << " (e2a=" << rc_e2a
                << " hpr=" << rc_hpr
                << " ref=" << rc_ref << ")" << std::endl;
      report("read E2A-HPr PDB files", false);
      delete mol_e2a; delete mol_hpr; delete mol_ref;
      return;
   }
   report("read E2A-HPr PDB files", true);

   // 1HDN is a 30-model NMR ensemble — use only model 1
   keep_model_1(mol_hpr);

   // 1GGR also has 3 models — use model 1
   keep_model_1(mol_ref);

   std::cout << "  E2A (1F3G): chain A" << std::endl;
   std::cout << "  HPr (1HDN): chain A, model 1 of NMR ensemble" << std::endl;

   // Superimpose 1GGR chain A onto 1F3G chain A so that the native HPr
   // COM (1GGR chain B) is in the same coordinate frame as our receptor.
   // Collect matching CA atoms by residue number.
   std::vector<clipper::Coord_orth> cas_ref;  // 1GGR chain A
   std::vector<clipper::Coord_orth> cas_tgt;  // 1F3G chain A

   {
      // Build a map of resno -> CA position for 1F3G
      std::map<int, clipper::Coord_orth> f3g_ca_map;
      int sel_f = mol_e2a->NewSelection();
      mol_e2a->SelectAtoms(sel_f, 1, "A",
                           mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                           "*", " CA ", "*", "*");
      mmdb::Atom **atoms_f = nullptr;
      int n_f = 0;
      mol_e2a->GetSelIndex(sel_f, atoms_f, n_f);
      for (int i=0; i<n_f; i++)
         f3g_ca_map[atoms_f[i]->GetSeqNum()] =
            clipper::Coord_orth(atoms_f[i]->x, atoms_f[i]->y, atoms_f[i]->z);
      mol_e2a->DeleteSelection(sel_f);

      // Collect matching pairs from 1GGR chain A
      int sel_r = mol_ref->NewSelection();
      mol_ref->SelectAtoms(sel_r, 1, "A",
                           mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                           "*", " CA ", "*", "*");
      mmdb::Atom **atoms_r = nullptr;
      int n_r = 0;
      mol_ref->GetSelIndex(sel_r, atoms_r, n_r);
      for (int i=0; i<n_r; i++) {
         int resno = atoms_r[i]->GetSeqNum();
         auto it = f3g_ca_map.find(resno);
         if (it != f3g_ca_map.end()) {
            cas_ref.push_back(clipper::Coord_orth(atoms_r[i]->x,
                                                  atoms_r[i]->y,
                                                  atoms_r[i]->z));
            cas_tgt.push_back(it->second);
         }
      }
      mol_ref->DeleteSelection(sel_r);
   }

   std::cout << "  Superposition: matched " << cas_ref.size()
             << " CA atoms between 1GGR and 1F3G chain A" << std::endl;

   // LSQ superposition: transform that maps 1GGR onto 1F3G frame
   clipper::RTop_orth rtop(cas_ref, cas_tgt);

   // Apply the transform to the native HPr COM (chain B of 1GGR)
   clipper::Coord_orth raw_com_hpr = compute_com(mol_ref, "B");
   clipper::Coord_orth native_com_hpr = rtop * raw_com_hpr;

   std::cout << "  Native HPr COM (1GGR chain B, in 1F3G frame): ("
             << std::setprecision(1) << std::fixed
             << native_com_hpr.x() << ", "
             << native_com_hpr.y() << ", "
             << native_com_hpr.z() << ")\n";

   // ---------------------------------------------------------------
   // Define active/passive residues from real NMR CSP experiments
   // ---------------------------------------------------------------
   //
   // E2A active residues: from NMR titration with HPr
   //   (van Nuland et al., J. Mol. Biol. 246, 180-193, 1995)
   //
   std::vector<int> e2a_active_resnos = {38, 40, 45, 46, 69, 71, 78, 80, 94, 96, 141};
   std::vector<int> e2a_passive_resnos = {37, 39, 43, 68, 72, 97, 109, 132};

   // HPr active residues: from NMR titration with E2A
   //   (Zhou & Bhatt, Biochemistry 37, 6269-6277, 1998)
   //
   std::vector<int> hpr_active_resnos = {15, 16, 17, 20, 48, 51, 52, 54, 56};
   std::vector<int> hpr_passive_resnos = {12, 21, 24, 47, 49, 57, 85};

   // Build active_passive_residues_t for each protein.
   // E2A is chain A in 1F3G, HPr is chain A in 1HDN.
   active_passive_residues_t ap_e2a;
   for (int r : e2a_active_resnos)
      ap_e2a.active.push_back(coot::residue_spec_t("A", r, ""));
   for (int r : e2a_passive_resnos)
      ap_e2a.passive.push_back(coot::residue_spec_t("A", r, ""));

   active_passive_residues_t ap_hpr;
   for (int r : hpr_active_resnos)
      ap_hpr.active.push_back(coot::residue_spec_t("A", r, ""));
   for (int r : hpr_passive_resnos)
      ap_hpr.passive.push_back(coot::residue_spec_t("A", r, ""));

   std::cout << "  E2A: " << ap_e2a.active.size() << " active, "
             << ap_e2a.passive.size() << " passive residues (NMR CSP)\n";
   std::cout << "  HPr: " << ap_hpr.active.size() << " active, "
             << ap_hpr.passive.size() << " passive residues (NMR CSP)\n";

   // Extract molecule data.
   // E2A is the fixed receptor; HPr is the mobile ligand.
   // For AIRs: E2A active residues generate restraints to HPr's
   // active+passive set (and vice versa).
   // is_partner=false means this molecule's active residues are the
   // "source" of AIRs; is_partner=true means this molecule provides
   // the "partner" (active+passive) targets.
   molecule_data_t data_e2a = extract_molecule_data(mol_e2a, ap_e2a, false);
   molecule_data_t data_hpr = extract_molecule_data(mol_hpr, ap_hpr, true);

   std::cout << "  Atom counts: E2A=" << data_e2a.atoms.size()
             << " HPr=" << data_hpr.atoms.size() << std::endl;
   std::cout << "  AIR residues: E2A active="
             << data_e2a.active_residue_atoms.size()
             << " HPr partner="
             << data_hpr.partner_residue_atoms.size() << std::endl;

   report("E2A has atoms", data_e2a.atoms.size() > 100);
   report("HPr has atoms", data_hpr.atoms.size() > 100);
   report("E2A has active AIR residues", data_e2a.active_residue_atoms.size() > 0);
   report("HPr has partner AIR residues", data_hpr.partner_residue_atoms.size() > 0);

   // Run rigid body docking
   docking_parameters_t params;
   params.n_trials = 200;
   params.n_keep = 50;
   params.n_orient_cycles = 4;
   params.separation_distance = -1.0;  // auto from protein sizes
   params.n_threads = 4;

   rigid_body_docker_t docker(data_e2a, data_hpr, params);
   std::vector<docking_result_t> results = docker.dock();

   report("got E2A-HPr docking results", !results.empty());

   if (!results.empty()) {
      // Print top results
      int n_show = std::min(5, static_cast<int>(results.size()));
      for (int i=0; i<n_show; i++) {
         const docking_result_t &r = results[i];
         double d = clipper::Coord_orth::length(r.translation, native_com_hpr);
         std::cout << "  Result " << i << ": E=" << std::setprecision(1) << std::fixed
                   << r.e_total
                   << " (VDW=" << r.e_vdw
                   << " Elec=" << r.e_elec
                   << " AIR=" << r.e_air
                   << ") dist=" << d << " A" << std::endl;
      }

      const docking_result_t &best = results[0];
      double best_dist = clipper::Coord_orth::length(best.translation, native_com_hpr);
      std::cout << "\n  Best energy: " << std::setprecision(2) << std::fixed
                << best.e_total << " kcal/mol\n";
      std::cout << "  Distance from native COM: " << best_dist << " A\n";

      // With real NMR data (not oracle AIRs), we expect the rigid body
      // stage to place the ligand in the correct region.  Within 25 A
      // is a reasonable expectation for rigid body + imperfect restraints.
      report("best docked COM within 30 A of native", best_dist < 30.0);
      report("best docked COM within 20 A of native", best_dist < 20.0);
      report("best total energy is negative", best.e_total < 0);

      // Write output files
      std::cout << "\n  Writing output files...\n";

      // Write 1GGR superimposed into the 1F3G coordinate frame as
      // the ground truth reference.  This lets you open e2a-hpr-receptor.pdb
      // + e2a-hpr-native.pdb + e2a-hpr-result-*.pdb together and see
      // how the docked HPr compares to the crystallographic answer.
      {
         mmdb::Manager *ref_copy = new mmdb::Manager();
         ref_copy->Copy(mol_ref, mmdb::MMDBFCM_All);
         // Apply the superposition transform to ALL atoms
         int sel_all = ref_copy->NewSelection();
         ref_copy->SelectAtoms(sel_all, 1, "*",
                               mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                               "*", "*", "*", "*");
         mmdb::Atom **ref_atoms = nullptr;
         int n_ref_atoms = 0;
         ref_copy->GetSelIndex(sel_all, ref_atoms, n_ref_atoms);
         for (int i=0; i<n_ref_atoms; i++) {
            clipper::Coord_orth orig(ref_atoms[i]->x,
                                     ref_atoms[i]->y,
                                     ref_atoms[i]->z);
            clipper::Coord_orth transformed = rtop * orig;
            ref_atoms[i]->x = transformed.x();
            ref_atoms[i]->y = transformed.y();
            ref_atoms[i]->z = transformed.z();
         }
         ref_copy->DeleteSelection(sel_all);
         ref_copy->WritePDBASCII("e2a-hpr-native.pdb");
         std::cout << "  Wrote e2a-hpr-native.pdb (1GGR superimposed onto"
                   << " 1F3G frame)" << std::endl;
         delete ref_copy;
      }

      mol_e2a->WritePDBASCII("e2a-hpr-receptor.pdb");
      std::cout << "  Wrote e2a-hpr-receptor.pdb (E2A, original coords)" << std::endl;

      int n_write = std::min(5, static_cast<int>(results.size()));
      for (int i=0; i<n_write; i++) {
         std::string fname = "e2a-hpr-result-" + std::to_string(i) + ".pdb";
         // HPr COM for apply_dock_transform: the docker centres HPr at
         // its COM, so we compute that from the unbound structure
         clipper::Coord_orth com_hpr = compute_com(mol_hpr);
         write_docked_pdb(mol_e2a, mol_hpr, results[i], com_hpr, fname);
      }

      write_results_csv(results, native_com_hpr, "e2a-hpr-results.csv");

      // ---------------------------------------------------------------
      // Build native HPr CA map for l-RMSD style CA-CA distances.
      // 1GGR chain B residues are numbered 301-385; HPr (1HDN) uses 1-85.
      // We transform 1GGR by rtop into the 1F3G frame, then store
      // CA positions keyed by (resno - 300) to match the docked HPr.
      // ---------------------------------------------------------------
      std::map<int, clipper::Coord_orth> native_hpr_ca_map;
      {
         int sel_b = mol_ref->NewSelection();
         mol_ref->SelectAtoms(sel_b, 1, "B",
                              mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                              "*", " CA ", "*", "*");
         mmdb::Atom **b_atoms = nullptr;
         int n_b = 0;
         mol_ref->GetSelIndex(sel_b, b_atoms, n_b);
         for (int i=0; i<n_b; i++) {
            clipper::Coord_orth orig(b_atoms[i]->x, b_atoms[i]->y,
                                     b_atoms[i]->z);
            clipper::Coord_orth tfm = rtop * orig;
            int mapped_resno = b_atoms[i]->GetSeqNum() - 300;
            native_hpr_ca_map[mapped_resno] = tfm;
         }
         mol_ref->DeleteSelection(sel_b);
         std::cout << "  Native HPr CA map: " << native_hpr_ca_map.size()
                   << " CA atoms (chain B of 1GGR, resno-300)\n";
      }

      // Build HPr (unbound) CA positions keyed by resno, for
      // transforming by each docked pose.
      struct ca_entry_t {
         int resno;
         clipper::Coord_orth position; // original (unbound) position
      };
      std::vector<ca_entry_t> hpr_ca_entries;
      {
         int sel_ca = mol_hpr->NewSelection();
         mol_hpr->SelectAtoms(sel_ca, 1, "A",
                              mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                              "*", " CA ", "*", "*");
         mmdb::Atom **ca_atoms = nullptr;
         int n_ca = 0;
         mol_hpr->GetSelIndex(sel_ca, ca_atoms, n_ca);
         for (int i=0; i<n_ca; i++) {
            ca_entry_t e;
            e.resno = ca_atoms[i]->GetSeqNum();
            e.position = clipper::Coord_orth(ca_atoms[i]->x,
                                              ca_atoms[i]->y,
                                              ca_atoms[i]->z);
            hpr_ca_entries.push_back(e);
         }
         mol_hpr->DeleteSelection(sel_ca);
      }

      // Compute HPr COM for transforming CA positions by docked pose
      clipper::Coord_orth com_hpr_dock = compute_com(mol_hpr);

      // Lambda: compute mean CA-CA distance for a docked pose
      auto mean_ca_ca_dist = [&](const docking_result_t &pose) -> double {
         glm::mat3 R = glm::mat3_cast(pose.rotation);
         double sum = 0;
         int count = 0;
         for (const auto &entry : hpr_ca_entries) {
            auto it = native_hpr_ca_map.find(entry.resno);
            if (it == native_hpr_ca_map.end()) continue;
            // Centre at HPr COM, rotate, translate
            clipper::Coord_orth centred(
               entry.position.x() - com_hpr_dock.x(),
               entry.position.y() - com_hpr_dock.y(),
               entry.position.z() - com_hpr_dock.z());
            glm::vec3 v(static_cast<float>(centred.x()),
                         static_cast<float>(centred.y()),
                         static_cast<float>(centred.z()));
            glm::vec3 rv = R * v;
            clipper::Coord_orth docked(rv.x + pose.translation.x(),
                                        rv.y + pose.translation.y(),
                                        rv.z + pose.translation.z());
            double d = clipper::Coord_orth::length(docked, it->second);
            sum += d;
            count++;
         }
         return (count > 0) ? sum / count : 999.0;
      };

      // ---------------------------------------------------------------
      // Stage 2: Semi-flexible simulated annealing refinement
      // ---------------------------------------------------------------
      //
      // Take the top 5 rigid body results and refine them through
      // three SA phases (rigid body, side-chain, backbone+SC).

      std::cout << "\n--- Semi-flexible refinement (Stage 2) ---\n";

      int n_refine = std::min(50, static_cast<int>(results.size()));
      std::vector<docking_result_t> top_results(results.begin(),
                                                 results.begin() + n_refine);

      semiflex_parameters_t sf_params;
      // Use fewer steps for testing speed
      sf_params.n_steps_rigid = 200;
      sf_params.n_steps_sc = 400;
      sf_params.n_steps_both = 400;
      sf_params.n_threads = 4;

      semiflex_refiner_t refiner(data_e2a, data_hpr, params, sf_params);
      std::vector<semiflex_result_t> sf_results = refiner.refine(top_results);

      report("got semiflex results", !sf_results.empty());

      if (!sf_results.empty()) {
         std::cout << "\n  Semi-flexible results (top 10 by energy):\n";
         int n_show_sf = std::min(10, static_cast<int>(sf_results.size()));
         for (int i=0; i<n_show_sf; i++) {
            const semiflex_result_t &sr = sf_results[i];
            double d = clipper::Coord_orth::length(sr.pose.translation,
                                                    native_com_hpr);
            std::cout << "  SF " << std::setw(2) << i
                      << ": E=" << std::setprecision(1) << std::fixed
                      << std::setw(8) << sr.e_total
                      << " (VDW=" << std::setw(7) << sr.pose.e_vdw
                      << " Elec=" << std::setw(7) << sr.pose.e_elec
                      << " AIR=" << std::setw(6) << sr.pose.e_air
                      << ") dist=" << std::setw(5) << d << " A"
                      << " flex=" << sr.flex_residues.size() << " res"
                      << std::endl;
         }

         // Check that refinement improved energy for at least some structures
         int n_improved = 0;
         for (unsigned int i=0; i<sf_results.size() && i<top_results.size(); i++) {
            if (sf_results[i].e_total < top_results[i].e_total)
               n_improved++;
         }
         std::cout << "  Energy improved for " << n_improved << "/"
                   << sf_results.size() << " structures" << std::endl;
         report("at least one structure improved by semiflex",
                n_improved > 0);

         // Find the structure closest to native
         double best_dist = 1e30;
         int best_dist_idx = -1;
         for (unsigned int i=0; i<sf_results.size(); i++) {
            double d = clipper::Coord_orth::length(
               sf_results[i].pose.translation, native_com_hpr);
            if (d < best_dist) { best_dist = d; best_dist_idx = i; }
         }
         std::cout << "  Closest to native: SF " << best_dist_idx
                   << " at " << std::setprecision(1) << best_dist
                   << " A (E=" << sf_results[best_dist_idx].e_total
                   << ")" << std::endl;

         // Write top 10 semiflex results
         int n_write_sf = std::min(10, static_cast<int>(sf_results.size()));
         for (int i=0; i<n_write_sf; i++) {
            std::string fname = "e2a-hpr-semiflex-" + std::to_string(i) + ".pdb";
            clipper::Coord_orth com_hpr_sf = compute_com(mol_hpr);
            write_docked_pdb(mol_e2a, mol_hpr, sf_results[i].pose,
                             com_hpr_sf, fname);
         }
         std::cout << "  Wrote " << n_write_sf
                   << " semiflex PDB files" << std::endl;

         // Write CSV with energy, COM distance, and mean CA-CA distance
         {
            std::ofstream csv("e2a-hpr-semiflex.csv");
            csv << "rank,e_total,e_vdw,e_elec,e_air,com_dist,mean_ca_ca,n_flex\n";
            for (unsigned int i=0; i<sf_results.size(); i++) {
               const semiflex_result_t &sr = sf_results[i];
               double com_d = clipper::Coord_orth::length(
                  sr.pose.translation, native_com_hpr);
               double ca_d = mean_ca_ca_dist(sr.pose);
               csv << i << ","
                   << std::setprecision(2) << std::fixed
                   << sr.e_total << ","
                   << sr.pose.e_vdw << ","
                   << sr.pose.e_elec << ","
                   << sr.pose.e_air << ","
                   << com_d << ","
                   << ca_d << ","
                   << sr.flex_residues.size() << "\n";
            }
            std::cout << "  Wrote e2a-hpr-semiflex.csv ("
                      << sf_results.size() << " results)" << std::endl;
         }
      }
   }

   delete mol_e2a;
   delete mol_hpr;
   delete mol_ref;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

static bool file_exists(const std::string &path) {
   struct stat st;
   return stat(path.c_str(), &st) == 0;
}

int main(int argc, char **argv) {

   std::cout << "=== HADDOCK Docking Tests ===" << std::endl;

   // Pure math tests (no PDB needed)
   test_air_effective_distance();
   test_air_penalty();
   test_air_gradient();

   // Real docking test 1: 9V3F chain E (oracle AIRs from known interface)
   std::string pdb_file = "pdb9v3f.ent";
   if (argc > 1) {
      pdb_file = argv[1];
   }

   if (file_exists(pdb_file)) {
      test_dock_9v3f_chain_E(pdb_file);
   } else {
      std::cout << "\n  SKIP: " << pdb_file << " not found"
                << " (pass path as argument)" << std::endl;
   }

   // Real docking test 2: E2A-HPr (real NMR CSP data)
   std::string pdb_e2a = "pdb1f3g.ent";
   std::string pdb_hpr = "pdb1hdn.ent";
   std::string pdb_ref = "pdb1ggr.ent";

   if (file_exists(pdb_e2a) && file_exists(pdb_hpr) && file_exists(pdb_ref)) {
      test_dock_e2a_hpr(pdb_e2a, pdb_hpr, pdb_ref);
   } else {
      std::cout << "\n  SKIP: E2A-HPr PDB files not found"
                << " (need pdb1f3g.ent, pdb1hdn.ent, pdb1ggr.ent)"
                << std::endl;
   }

   std::cout << "\n=== Summary ===" << std::endl;
   std::cout << "  Passed: " << n_tests_passed << std::endl;
   std::cout << "  Failed: " << n_tests_failed << std::endl;

   return (n_tests_failed > 0) ? 1 : 0;
}
