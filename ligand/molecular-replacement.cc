/*
 * ligand/molecular-replacement.cc
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
#include <chrono>
#include <cmath>
#include <algorithm>
#include <thread>
#include <map>
#include <cstdlib>

#include <clipper/core/xmap.h>
#include <clipper/core/map_utils.h>
#include <clipper/core/map_interp.h>
#include <clipper/contrib/edcalc.h>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "molecular-replacement.hh"
#include "rigid-body.hh"
#include "utils/coot-utils.hh"
#include "utils/split-indices.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/crowther.hh"
#include "coot-utils/patterson.hh"

namespace {

   // Diagnostics for testing the MR pipeline against a known answer.
   //
   // If the environment variable COOT_MR_DIAGNOSTIC_TRUE_MODEL names a coordinates file
   // that contains the search model in its correct position (same chain ids, residue
   // numbers and atom names), then each stage of the search is compared with the true
   // rotation and translation.
   //
   // The MR placement convention is: placed = R * (x - model_centre) + t
   //
   class mr_truth_t {
   public:
      bool valid;
      glm::quat rotation;
      clipper::Coord_orth translation;
      std::map<std::string, clipper::Coord_orth> true_positions;
      mr_truth_t() : valid(false), rotation(1,0,0,0), translation(0,0,0) {}

      static std::string atom_key(mmdb::Atom *at) {
         return std::string(at->GetChainID()) + "/" + std::to_string(at->GetSeqNum()) +
            std::string(at->GetInsCode()) + "/" + std::string(at->GetAtomName()) + "," +
            std::string(at->altLoc);
      }

      static std::vector<mmdb::Atom *> model_atoms(mmdb::Manager *mol) {
         std::vector<mmdb::Atom *> v;
         mmdb::Model *model_p = mol->GetModel(1);
         if (! model_p) return v;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  if (! at->isTer()) v.push_back(at);
               }
            }
         }
         return v;
      }

      void init(mmdb::Manager *mol_model, const clipper::Coord_orth &model_centre) {
         const char *fn = std::getenv("COOT_MR_DIAGNOSTIC_TRUE_MODEL");
         if (! fn) return;
         mmdb::Manager *mol_true = new mmdb::Manager;
         mmdb::ERROR_CODE err = mol_true->ReadCoorFile(fn);
         if (err) {
            std::cout << "MR-DIAG:: failed to read true model " << fn << std::endl;
            delete mol_true;
            return;
         }
         for (mmdb::Atom *at : model_atoms(mol_true))
            true_positions[atom_key(at)] = clipper::Coord_orth(at->x, at->y, at->z);
         delete mol_true;

         std::vector<clipper::Coord_orth> src, tgt;
         for (mmdb::Atom *at : model_atoms(mol_model)) {
            auto it = true_positions.find(atom_key(at));
            if (it != true_positions.end()) {
               src.push_back(clipper::Coord_orth(at->x, at->y, at->z));
               tgt.push_back(it->second);
            }
         }
         std::cout << "MR-DIAG:: true model " << fn << " " << true_positions.size()
                   << " atoms, " << src.size() << " matched to search model" << std::endl;
         if (src.size() < 3) return;

         clipper::RTop_orth rtop(src, tgt);
         double sum_sq = 0.0;
         for (unsigned int i=0; i<src.size(); i++)
            sum_sq += (src[i].transform(rtop) - tgt[i]).lengthsq();
         clipper::Mat33<> r = rtop.rot();
         glm::mat3 m;
         for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
               m[j][i] = static_cast<float>(r(i,j)); // glm is column-major
         rotation = glm::quat_cast(m);
         translation = model_centre.transform(rtop);
         valid = true;
         std::cout << "MR-DIAG:: superposition rmsd " << std::fixed << std::setprecision(3)
                   << std::sqrt(sum_sq / src.size()) << " A; true rotation is "
                   << std::setprecision(1) << angle_to(glm::quat(1,0,0,0))
                   << " deg from identity; true translation (" << translation.x() << ", "
                   << translation.y() << ", " << translation.z() << ")" << std::endl;
      }

      // degrees from the nearest of identity and the 2-folds about x, y and z - i.e. how close
      // is this rotation to one that maps an axis-aligned (orthorhombic) lattice onto itself?
      static double angle_to_lattice_ops(const glm::quat &q) {
         const glm::quat ops[4] = { glm::quat(1,0,0,0), glm::quat(0,1,0,0), glm::quat(0,0,1,0), glm::quat(0,0,0,1) };
         double best = 1000.0;
         for (const auto &op : ops) {
            double dot = std::abs(static_cast<double>(glm::dot(glm::normalize(q), op)));
            if (dot > 1.0) dot = 1.0;
            double a = glm::degrees(2.0 * std::acos(dot));
            if (a < best) best = a;
         }
         return best;
      }

      // degrees
      double angle_to(const glm::quat &q) const {
         double dot = std::abs(static_cast<double>(glm::dot(glm::normalize(q), rotation)));
         if (dot > 1.0) dot = 1.0;
         return glm::degrees(2.0 * std::acos(dot));
      }

      double distance_to(const clipper::Coord_orth &t) const {
         return std::sqrt((t - translation).lengthsq());
      }

      // rmsd from the true positions (the search model atoms only)
      double rmsd(mmdb::Manager *mol) const {
         double sum_sq = 0.0;
         unsigned int n = 0;
         for (mmdb::Atom *at : model_atoms(mol)) {
            auto it = true_positions.find(atom_key(at));
            if (it != true_positions.end()) {
               sum_sq += (clipper::Coord_orth(at->x, at->y, at->z) - it->second).lengthsq();
               n++;
            }
         }
         return (n > 0) ? std::sqrt(sum_sq / n) : -1.0;
      }

      double fraction_within(const clipper::Coord_orth &c, float radius) const {
         if (true_positions.empty()) return 0.0;
         double r_sq = radius * radius;
         unsigned int n = 0;
         for (const auto &tp : true_positions)
            if ((tp.second - c).lengthsq() < r_sq) n++;
         return static_cast<double>(n) / static_cast<double>(true_positions.size());
      }
   };

   // ZYZ Euler angles (radians) of a rotation, matching crowther_t (Omega = Rz(a) Ry(b) Rz(g))
   void diagnostic_quat_to_zyz(const glm::quat &q, double *alpha_p, double *beta_p, double *gamma_p) {
      glm::mat3 R = glm::mat3_cast(q);
      double cos_beta = std::clamp(static_cast<double>(R[2][2]), -1.0, 1.0);
      *beta_p = std::acos(cos_beta);
      if (std::abs(std::sin(*beta_p)) > 1e-6) {
         *alpha_p = std::atan2(static_cast<double>(R[2][1]), static_cast<double>(R[2][0]));
         *gamma_p = std::atan2(static_cast<double>(R[1][2]), static_cast<double>(-R[0][2]));
      } else {
         *alpha_p = std::atan2(static_cast<double>(R[1][0]), static_cast<double>(R[0][0]));
         *gamma_p = 0.0;
      }
   }

   // Run the phased translation function for the given rotation (as the main pipeline does)
   // and report where the true translation lies in it.
   void diagnostic_translation_at_rotation(const glm::quat &rot_q, const std::string &label,
                                           mmdb::Manager *mol_model,
                                           const clipper::Coord_orth &mol_centre,
                                           const clipper::Cell &local_cell,
                                           const clipper::Grid_sampling &local_gs,
                                           const clipper::HKL_data<clipper::datatypes::F_phi<float>> &fphi_obs_local,
                                           const clipper::HKL_info &local_hkls,
                                           const clipper::Xmap<float> &local_obs_map,
                                           const clipper::Coord_orth &local_origin_in_obs,
                                           const mr_truth_t &truth) {

      mmdb::Manager *mol_rot = new mmdb::Manager();
      mol_rot->Copy(mol_model, mmdb::MMDBFCM_All);
      int sel_rot = mol_rot->NewSelection();
      mol_rot->SelectAtoms(sel_rot, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mmdb::PPAtom atoms_rot = 0;
      int na_rot;
      mol_rot->GetSelIndex(sel_rot, atoms_rot, na_rot);
      glm::mat3 rot_mat = glm::mat3_cast(rot_q);
      for (int i=0; i<na_rot; i++) {
         glm::vec3 pos(atoms_rot[i]->x - mol_centre.x(), atoms_rot[i]->y - mol_centre.y(),
                       atoms_rot[i]->z - mol_centre.z());
         glm::vec3 rotated = rot_mat * pos;
         atoms_rot[i]->x = rotated.x;
         atoms_rot[i]->y = rotated.y;
         atoms_rot[i]->z = rotated.z;
      }
      clipper::Xmap<float> rotated_model_map =
         coot::util::calc_atom_map(mol_rot, sel_rot, local_cell, clipper::Spacegroup::p1(), local_gs);
      mol_rot->DeleteSelection(sel_rot);
      delete mol_rot;

      clipper::Xmap<float> tf_map;
      std::vector<coot::translation_search_result_t> tr =
         coot::phased_translation_search(fphi_obs_local, local_hkls, local_obs_map,
                                         rotated_model_map, 100000, &tf_map);

      double sum = 0.0, sum_sq = 0.0;
      unsigned int count = 0;
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = tf_map.first(); !ix.last(); ix.next()) {
         float v = tf_map[ix];
         sum += v; sum_sq += static_cast<double>(v) * v; count++;
      }
      double mean = sum / count;
      double sd = std::sqrt(std::max(sum_sq / count - mean * mean, 1e-30));

      clipper::Coord_orth true_local = truth.translation - local_origin_in_obs;
      float v_true = 0.0f;
      clipper::Interp_linear::interp(tf_map, true_local.coord_frac(local_cell).coord_map(local_gs), v_true);

      std::cout << "MR-DIAG:: TF " << label << " rotation (" << std::fixed << std::setprecision(1)
                << truth.angle_to(rot_q) << " deg from true): TF at true translation "
                << std::setprecision(2) << (v_true - mean) / sd << " sigma" << std::endl;

      // distinct peaks: a grid point at least 5 A from all higher-scoring retained peaks
      std::vector<coot::translation_search_result_t> distinct;
      for (auto &r : tr) {
         clipper::Coord_orth obs_pos = r.position + local_origin_in_obs;
         bool close = false;
         for (const auto &d : distinct)
            if ((d.position - obs_pos).lengthsq() < 25.0) { close = true; break; }
         if (! close) distinct.push_back(coot::translation_search_result_t(obs_pos, r.score));
         if (distinct.size() >= 8) break;
      }
      for (unsigned int i=0; i<distinct.size(); i++)
         std::cout << "MR-DIAG::    TF peak " << i << " " << std::fixed << std::setprecision(2)
                   << distinct[i].score << " sigma at (" << std::setprecision(1)
                   << distinct[i].position.x() << ", " << distinct[i].position.y() << ", "
                   << distinct[i].position.z() << ") " << truth.distance_to(distinct[i].position)
                   << " A from true translation" << std::endl;
   }
}

namespace {

   // Diagnostic experiment: if both COOT_MR_DIAGNOSTIC_CALC_TARGET and COOT_MR_DIAGNOSTIC_TRUE_MODEL
   // are set, replace the observed map by a map calculated from the true model (on the same
   // cell and grid) - i.e. test the search with perfect data.
   bool diagnostic_calc_target_map(const clipper::Xmap<float> &xmap_ref, clipper::Xmap<float> *xmap_calc_p) {
      const char *fn = std::getenv("COOT_MR_DIAGNOSTIC_TRUE_MODEL");
      if (! fn) return false;
      if (! std::getenv("COOT_MR_DIAGNOSTIC_CALC_TARGET")) return false;
      mmdb::Manager *mol_true = new mmdb::Manager;
      if (mol_true->ReadCoorFile(fn)) { delete mol_true; return false; }
      int sel = mol_true->NewSelection();
      mol_true->SelectAtoms(sel, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      *xmap_calc_p = coot::util::calc_atom_map(mol_true, sel, xmap_ref.cell(), xmap_ref.spacegroup(),
                                               xmap_ref.grid_sampling());
      mol_true->DeleteSelection(sel);
      delete mol_true;
      std::cout << "MR-DIAG:: using a map calculated from the true model as the target map" << std::endl;
      return true;
   }
}

std::vector<coot::mr_solution_t>
coot::molecular_replacement_search(const clipper::Xmap<float> &xmap_obs_in,
                                   mmdb::Manager *mol_model,
                                   const clipper::Coord_orth &target_centre,
                                   int n_rotation_solutions,
                                   int n_translation_solutions) {

   std::vector<mr_solution_t> all_solutions;

   clipper::Xmap<float> xmap_calc_diag;
   bool use_calc = diagnostic_calc_target_map(xmap_obs_in, &xmap_calc_diag);
   const clipper::Xmap<float> &xmap_obs = use_calc ? xmap_calc_diag : xmap_obs_in;

   auto t_total_start = std::chrono::high_resolution_clock::now();

   // --- Determine radii from the model's radius of gyration ---

   std::pair<bool, double> rg_pair = coot::radius_of_gyration(mol_model);
   if (! rg_pair.first) {
      std::cout << "ERROR:: molecular_replacement_search: could not compute radius of gyration"
                << std::endl;
      return all_solutions;
   }
   double rg = rg_pair.second;
   float fragment_radius = static_cast<float>(rg * 1.5);
   if (fragment_radius < 15.0f) fragment_radius = 15.0f;
   float crowther_radius = fragment_radius * 0.9f;

   std::cout << "INFO:: molecular_replacement_search"
             << " radius of gyration: " << std::fixed << std::setprecision(1)
             << rg << " A, fragment radius: " << fragment_radius
             << " A, Crowther radius: " << crowther_radius << " A" << std::endl;

   // --- Iterative density-weighted recentring ---
   //
   // The user's click may be offset from the true molecule centre. The Tukey
   // window must be symmetric about the molecule for the Patterson to match the
   // model's. Iteratively step toward the centroid of high density, limiting
   // each step to 3 A (cf. helix_placement::move_helix_centre_point_guess).

   clipper::Coord_orth recentred = target_centre;
   auto t_recentre_start = std::chrono::high_resolution_clock::now();
   {
      float max_step = 3.0f;
      float max_total_sq = static_cast<float>(rg * rg);
      float centroid_radius = static_cast<float>(rg * 0.75);
      if (centroid_radius < 7.0f) centroid_radius = 7.0f;
      int max_iterations = 10;
      for (int iter=0; iter<max_iterations; iter++) {
         coot::util::map_fragment_info_t mf_rough(xmap_obs, recentred, centroid_radius, true, 0.0f);
         if (mf_rough.xmap.is_null()) break;

         double sum = 0.0, sum_sq = 0.0;
         unsigned int n = 0;
         clipper::Xmap<float>::Map_reference_index inx;
         for (inx = mf_rough.xmap.first(); !inx.last(); inx.next()) {
            float v = mf_rough.xmap[inx];
            if (v != 0.0f) {
               sum += v;
               sum_sq += v * v;
               n++;
            }
         }
         if (n == 0) break;
         double mean = sum / n;
         double var = sum_sq / n - mean * mean;
         double sigma = (var > 0.0) ? std::sqrt(var) : 0.0;
         float threshold = static_cast<float>(mean + 3.0 * sigma);

         double wx = 0.0, wy = 0.0, wz = 0.0, wt = 0.0;
         clipper::Cell rough_cell = mf_rough.xmap.cell();
         clipper::Grid_sampling rough_gs = mf_rough.xmap.grid_sampling();
         for (inx = mf_rough.xmap.first(); !inx.last(); inx.next()) {
            float v = mf_rough.xmap[inx];
            if (v > threshold) {
               clipper::Coord_orth p = inx.coord().coord_frac(rough_gs).coord_orth(rough_cell);
               double px = p.x(); if (px > 0.5*rough_cell.a()) px -= rough_cell.a();
               double py = p.y(); if (py > 0.5*rough_cell.b()) py -= rough_cell.b();
               double pz = p.z(); if (pz > 0.5*rough_cell.c()) pz -= rough_cell.c();
               wx += v * px;
               wy += v * py;
               wz += v * pz;
               wt += v;
            }
         }
         if (wt <= 0.0) break;
         wx /= wt; wy /= wt; wz /= wt;

         double step_dist = std::sqrt(wx*wx + wy*wy + wz*wz);
         if (step_dist < 0.1) break;

         if (step_dist > max_step) {
            double scale = max_step / step_dist;
            wx *= scale; wy *= scale; wz *= scale;
         }

         recentred = clipper::Coord_orth(recentred.x() + wx,
                                         recentred.y() + wy,
                                         recentred.z() + wz);

         if ((recentred - target_centre).lengthsq() > max_total_sq) {
            std::cout << "INFO:: molecular_replacement_search"
                      << " recentring movement limit reached" << std::endl;
            break;
         }
      }
      double total_shift = std::sqrt((recentred - target_centre).lengthsq());
      if (total_shift > 0.1) {
         std::cout << "INFO:: molecular_replacement_search recentred from ("
                   << std::fixed << std::setprecision(1)
                   << target_centre.x() << ", " << target_centre.y() << ", " << target_centre.z()
                   << ") to ("
                   << recentred.x() << ", " << recentred.y() << ", " << recentred.z()
                   << ") shift " << total_shift << " A" << std::endl;
      }
   }
   auto t_recentre_end = std::chrono::high_resolution_clock::now();
   double t_recentre_ms = std::chrono::duration<double, std::milli>(t_recentre_end - t_recentre_start).count();
   std::cout << "INFO:: molecular_replacement_search recentring: " << std::fixed
             << std::setprecision(1) << t_recentre_ms << " ms" << std::endl;

   // --- Extract target map fragment ---

   coot::util::map_fragment_info_t mf(xmap_obs, recentred, fragment_radius, true, 0.5f);
   if (mf.xmap.is_null()) {
      std::cout << "ERROR:: molecular_replacement_search: target map fragment is null" << std::endl;
      return all_solutions;
   }

   // --- Prepare the model: centre it and compute an atom map ---

   std::pair<bool, clipper::Coord_orth> mol_centre = coot::centre_of_molecule(mol_model);
   if (! mol_centre.first) {
      std::cout << "ERROR:: molecular_replacement_search: could not find centre of molecule"
                << std::endl;
      return all_solutions;
   }

   mr_truth_t truth;
   truth.init(mol_model, mol_centre.second);
   if (truth.valid) {
      std::cout << "MR-DIAG:: target_centre is " << std::fixed << std::setprecision(1)
                << truth.distance_to(target_centre) << " A from the true centre; recentred is "
                << truth.distance_to(recentred) << " A from the true centre" << std::endl;
      std::cout << "MR-DIAG:: fraction of true-model atoms within the fragment radius "
                << fragment_radius << " A of: true centre " << std::setprecision(2)
                << truth.fraction_within(truth.translation, fragment_radius)
                << ", recentred " << truth.fraction_within(recentred, fragment_radius)
                << ", within local-box half-side " << fragment_radius + 10.0f
                << " A of recentred " << truth.fraction_within(recentred, fragment_radius + 10.0f)
                << std::endl;
   }

   mmdb::Manager *mol = new mmdb::Manager();
   mol->Copy(mol_model, mmdb::MMDBFCM_All);

   int SelHnd = mol->NewSelection();
   mol->SelectAtoms(SelHnd, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                    "*", "*", "*", "*");

   std::pair<clipper::Coord_orth, clipper::Coord_orth> ext = coot::util::extents(mol, SelHnd);
   double x_range = ext.second.x() - ext.first.x();
   double y_range = ext.second.y() - ext.first.y();
   double z_range = ext.second.z() - ext.first.z();
   double border = 10.0;
   double r_90 = clipper::Util::d2rad(90.0);
   clipper::Cell_descr cd_model(x_range + 2.0*border, y_range + 2.0*border,
                                z_range + 2.0*border, r_90, r_90, r_90);
   clipper::Cell cell_model(cd_model);
   clipper::Spacegroup spacegroup = clipper::Spacegroup::p1();
   clipper::Resolution reso(3.0);
   clipper::Grid_sampling gs_model(spacegroup, cell_model, reso);

   clipper::Coord_orth cell_centre_model(0.5 * cd_model.a(), 0.5 * cd_model.b(), 0.5 * cd_model.c());
   clipper::Coord_orth shift(cell_centre_model.x() - mol_centre.second.x(),
                             cell_centre_model.y() - mol_centre.second.y(),
                             cell_centre_model.z() - mol_centre.second.z());

   mmdb::PPAtom sel_atoms = 0;
   int n_atoms;
   mol->GetSelIndex(SelHnd, sel_atoms, n_atoms);

   std::cout << "INFO:: molecular_replacement_search n_atoms in copy: " << n_atoms
             << " cell: " << std::fixed << std::setprecision(1)
             << cd_model.a() << " x " << cd_model.b() << " x " << cd_model.c()
             << " grid: " << gs_model.nu() << "x" << gs_model.nv() << "x" << gs_model.nw()
             << std::endl;

   for (int i=0; i<n_atoms; i++) {
      sel_atoms[i]->x += shift.x();
      sel_atoms[i]->y += shift.y();
      sel_atoms[i]->z += shift.z();
   }

   clipper::Xmap<float> model_map = coot::util::calc_atom_map(mol, SelHnd, cell_model, spacegroup, gs_model);

   // Clean NaN values from the model map. These can arise from TER records
   // or atoms with invalid B-factors being included in the atom selection.
   // NaN values corrupt the FFT, making all structure factors NaN.
   {
      unsigned int nan_count = 0;
      clipper::Xmap<float>::Map_reference_index inx;
      for (inx = model_map.first(); !inx.last(); inx.next()) {
         if (std::isnan(model_map[inx]) || std::isinf(model_map[inx])) {
            model_map[inx] = 0.0f;
            nan_count++;
         }
      }
      if (nan_count > 0)
         std::cout << "INFO:: molecular_replacement_search cleaned " << nan_count
                   << " NaN/inf values from model map" << std::endl;
   }

   coot::util::map_fragment_info_t mf_model(model_map, cell_centre_model, fragment_radius, true, 0.5f);

   // --- Crowther rotation function ---

   auto t_rot_start = std::chrono::high_resolution_clock::now();

   int ell_max = 15;
   int n_beta = 36;
   coot::crowther_t crowther(ell_max, crowther_radius, reso);
   crowther.set_target(mf.xmap);
   crowther.set_search(mf_model.xmap);

   std::vector<coot::rotation_function_result_t> rot_results = crowther.compute(n_beta);

   auto t_rot_end = std::chrono::high_resolution_clock::now();
   double t_rot_ms = std::chrono::duration<double, std::milli>(t_rot_end - t_rot_start).count();
   std::cout << "INFO:: molecular_replacement_search rotation coarse search: "
             << rot_results.size() << " results in " << std::fixed << std::setprecision(1)
             << t_rot_ms << " ms" << std::endl;

   if (truth.valid && ! rot_results.empty()) {
      auto rank_of_score = [&rot_results](double s) {
         unsigned int n = 0;
         for (const auto &r : rot_results) if (r.score > s) n++;
         return n;
      };
      double a, b, g;
      diagnostic_quat_to_zyz(truth.rotation, &a, &b, &g);
      double s_true = crowther.evaluate_at_euler(a, b, g);
      diagnostic_quat_to_zyz(glm::inverse(truth.rotation), &a, &b, &g);
      double s_inv = crowther.evaluate_at_euler(a, b, g);
      std::cout << "MR-DIAG:: RF top coarse score " << std::scientific << std::setprecision(3)
                << rot_results[0].score << "; RF at true rotation " << s_true
                << " (would rank " << rank_of_score(s_true) << " of " << rot_results.size()
                << "); RF at inverse of true rotation " << s_inv << " (rank " << rank_of_score(s_inv)
                << ")" << std::endl;
      std::vector<coot::rotation_function_result_t> tv(1, coot::rotation_function_result_t(truth.rotation, s_true));
      auto tv_r = crowther.refine_orientations(tv, 1, 2.0f);
      tv_r = crowther.refine_orientations(tv_r, 1, 0.7f);
      std::cout << "MR-DIAG:: RF refined from true rotation: " << tv_r[0].score << " at "
                << std::fixed << std::setprecision(1) << truth.angle_to(tv_r[0].rotation)
                << " deg from true" << std::endl;
      unsigned int i_nearest = 0;
      double best_angle = 1000.0;
      for (unsigned int i=0; i<rot_results.size(); i++) {
         double ang = truth.angle_to(rot_results[i].rotation);
         if (ang < best_angle) { best_angle = ang; i_nearest = i; }
      }
      std::cout << "MR-DIAG:: nearest coarse grid point to true rotation: " << best_angle
                << " deg, rank " << i_nearest << " score " << std::scientific << std::setprecision(3)
                << rot_results[i_nearest].score << std::endl;
      // first coarse peak within 15 degrees of the truth
      for (unsigned int i=0; i<rot_results.size(); i++) {
         if (truth.angle_to(rot_results[i].rotation) < 15.0) {
            std::cout << "MR-DIAG:: first coarse result within 15 deg of true: rank " << i
                      << " score " << rot_results[i].score << " ("
                      << 100.0 * rot_results[i].score / rot_results[0].score << "% of top)" << std::endl;
            break;
         }
      }
      for (unsigned int i=0; i<std::min(static_cast<std::size_t>(10), rot_results.size()); i++)
         std::cout << "MR-DIAG::    coarse " << i << " score " << std::scientific << std::setprecision(3)
                   << rot_results[i].score << " " << std::fixed << std::setprecision(1)
                   << truth.angle_to(rot_results[i].rotation) << " deg from true, "
                   << mr_truth_t::angle_to_lattice_ops(rot_results[i].rotation)
                   << " deg from identity/2-fold-about-axis" << std::endl;
      auto print_q = [](const std::string &label, const glm::quat &q) {
         std::cout << "MR-DIAG-Q:: " << label << " " << std::fixed << std::setprecision(6)
                   << q.w << " " << q.x << " " << q.y << " " << q.z << std::endl;
      };
      print_q("true", truth.rotation);
      for (unsigned int i=0; i<std::min(static_cast<std::size_t>(30), rot_results.size()); i++)
         print_q("coarse-" + std::to_string(i) + " " + std::to_string(rot_results[i].score), rot_results[i].rotation);
      std::cout << "MR-DIAG:: true rotation is " << std::fixed << std::setprecision(1)
                << mr_truth_t::angle_to_lattice_ops(truth.rotation)
                << " deg from identity/2-fold-about-axis" << std::endl;
   }

   // Filter: keep only rotation peaks >= 70% of the top score
   if (! rot_results.empty()) {
      float top_rot = rot_results[0].score;
      float rot_cutoff = top_rot * 0.7f;
      auto it = std::remove_if(rot_results.begin(), rot_results.end(),
                                [rot_cutoff](const coot::rotation_function_result_t &r) {
                                   return r.score < rot_cutoff;
                                });
      unsigned int n_before = rot_results.size();
      rot_results.erase(it, rot_results.end());
      std::cout << "INFO:: molecular_replacement_search rotation pre-filter: "
                << n_before << " -> " << rot_results.size()
                << " (cutoff 70% of " << std::scientific << std::setprecision(2)
                << top_rot << ")" << std::endl;
   }

   // Two-pass refinement
   auto t_refine_start = std::chrono::high_resolution_clock::now();
   auto refined = crowther.refine_orientations(rot_results, n_rotation_solutions, 2.0f);
   refined = crowther.refine_orientations(refined, n_rotation_solutions, 0.7f);

   // Filter: keep only refined rotations >= 80% of the top refined score
   if (! refined.empty()) {
      float top_refined = refined[0].score;
      float refined_cutoff = top_refined * 0.8f;
      auto it = std::remove_if(refined.begin(), refined.end(),
                                [refined_cutoff](const coot::rotation_function_result_t &r) {
                                   return r.score < refined_cutoff;
                                });
      unsigned int n_before = refined.size();
      refined.erase(it, refined.end());
      std::cout << "INFO:: molecular_replacement_search rotation post-refinement filter: "
                << n_before << " -> " << refined.size()
                << " (cutoff 80% of " << std::scientific << std::setprecision(2)
                << top_refined << ")" << std::endl;
   }

   // Filter: remove rotation solutions within 4 degrees of a higher-ranked one
   {
      float min_angle_rad = glm::radians(4.0f);
      std::vector<coot::rotation_function_result_t> unique;
      unique.reserve(refined.size());
      for (const auto &r : refined) {
         bool too_close = false;
         for (const auto &u : unique) {
            float dot = std::abs(glm::dot(r.rotation, u.rotation));
            if (dot > 1.0f) dot = 1.0f;
            float angle = 2.0f * std::acos(dot);
            if (angle < min_angle_rad) {
               too_close = true;
               break;
            }
         }
         if (! too_close)
            unique.push_back(r);
      }
      unsigned int n_before = refined.size();
      refined = unique;
      if (refined.size() < n_before) {
         std::cout << "INFO:: molecular_replacement_search rotation proximity filter: "
                   << n_before << " -> " << refined.size()
                   << " (minimum separation 4 deg)" << std::endl;
      }
   }
   auto t_refine_end = std::chrono::high_resolution_clock::now();
   double t_refine_ms = std::chrono::duration<double, std::milli>(t_refine_end - t_refine_start).count();
   std::cout << "INFO:: molecular_replacement_search rotation refinement: " << std::fixed
             << std::setprecision(1) << t_refine_ms << " ms (" << refined.size()
             << " solutions)" << std::endl;

   if (truth.valid)
      for (unsigned int i=0; i<refined.size(); i++)
         std::cout << "MR-DIAG:: refined rotation " << i << " score " << std::scientific
                   << std::setprecision(3) << refined[i].score << " " << std::fixed << std::setprecision(1)
                   << truth.angle_to(refined[i].rotation) << " deg from true, "
                   << mr_truth_t::angle_to_lattice_ops(refined[i].rotation)
                   << " deg from identity/2-fold-about-axis" << std::endl;

   // --- Translation search for each refined rotation (parallelised) ---
   //
   // Two key optimisations:
   //
   // 1. Extract a local P1 box from the observed map around the recentred position
   //    (side ~ 2*fragment_radius + border). Both the observed and model
   //    density are computed on this smaller grid, so the FFTs are much
   //    cheaper than using the full EM cell.
   //
   // 2. Pre-compute the FFT of the observed local box once, rather than
   //    recomputing it for every rotation.

   auto t_trans_start = std::chrono::high_resolution_clock::now();

   // Build a local P1 box centred on the recentred position (not the user's target_centre)
   float local_box_radius = fragment_radius + 10.0f;
   double r_90_local = clipper::Util::d2rad(90.0);
   double local_side = 2.0 * static_cast<double>(local_box_radius);
   clipper::Cell_descr cd_local(local_side, local_side, local_side, r_90_local, r_90_local, r_90_local);
   clipper::Cell local_cell(cd_local);
   clipper::Resolution local_reso(3.0);
   clipper::Grid_sampling local_gs(spacegroup, local_cell, local_reso);
   clipper::Coord_orth local_cell_centre(0.5 * local_side, 0.5 * local_side, 0.5 * local_side);

   // Interpolate the observed map into the local box
   clipper::Xmap<float> local_obs_map(spacegroup, local_cell, local_gs);
   {
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = local_obs_map.first(); !ix.last(); ix.next()) {
         clipper::Coord_orth local_pos = ix.coord().coord_frac(local_gs).coord_orth(local_cell);
         clipper::Coord_orth obs_pos(local_pos.x() - local_cell_centre.x() + recentred.x(),
                                    local_pos.y() - local_cell_centre.y() + recentred.y(),
                                    local_pos.z() - local_cell_centre.z() + recentred.z());
         float val;
         clipper::Interp_linear::interp(xmap_obs, xmap_obs.coord_map(obs_pos), val);
         local_obs_map[ix] = val;
      }
   }

   std::cout << "INFO:: molecular_replacement_search translation local box: "
             << std::fixed << std::setprecision(1) << local_side << " A side, grid "
             << local_gs.nu() << "x" << local_gs.nv() << "x" << local_gs.nw() << std::endl;

   // Pre-compute the FFT of the local observed map (optimisation 2)
   clipper::HKL_info local_hkls(spacegroup, local_cell,
                                 clipper::Resolution(2.0 * coot::util::max_gridding(local_obs_map)), true);
   clipper::HKL_data<clipper::datatypes::F_phi<float>> fphi_obs_local(local_hkls);
   local_obs_map.fft_to(fphi_obs_local);

   if (truth.valid) {
      clipper::Coord_orth local_origin(recentred.x() - local_cell_centre.x(),
                                       recentred.y() - local_cell_centre.y(),
                                       recentred.z() - local_cell_centre.z());
      diagnostic_translation_at_rotation(truth.rotation, "true", mol_model, mol_centre.second,
                                         local_cell, local_gs, fphi_obs_local, local_hkls,
                                         local_obs_map, local_origin, truth);
      if (! refined.empty())
         diagnostic_translation_at_rotation(refined[0].rotation, "top-refined", mol_model, mol_centre.second,
                                            local_cell, local_gs, fphi_obs_local, local_hkls,
                                            local_obs_map, local_origin, truth);
   }

   unsigned int n_rot = refined.size();
   unsigned int n_trans_threads = coot::get_max_number_of_threads();
   if (n_trans_threads == 0) n_trans_threads = 4;
   if (n_trans_threads > n_rot) n_trans_threads = n_rot;

   std::vector<std::vector<mr_solution_t>> per_rotation_solutions(n_rot);
   auto trans_ranges = coot::atom_index_ranges(n_rot, n_trans_threads);
   std::vector<std::thread> trans_threads;
   trans_threads.reserve(trans_ranges.size());

   for (unsigned int t=0; t<trans_ranges.size(); t++) {
      trans_threads.emplace_back([&, t]() {
         for (unsigned int irot=trans_ranges[t].first; irot<trans_ranges[t].second; irot++) {

            const glm::quat &rot_q = refined[irot].rotation;
            float rot_score = refined[irot].score;
            glm::mat3 rot_mat = glm::mat3_cast(rot_q);

            // Fresh copy of the model, centred at origin, rotated
            mmdb::Manager *mol_rot = new mmdb::Manager();
            mol_rot->Copy(mol_model, mmdb::MMDBFCM_All);

            int sel_rot = mol_rot->NewSelection();
            mol_rot->SelectAtoms(sel_rot, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                                 "*", "*", "*", "*");
            mmdb::PPAtom atoms_rot = 0;
            int na_rot;
            mol_rot->GetSelIndex(sel_rot, atoms_rot, na_rot);

            for (int i=0; i<na_rot; i++) {
               atoms_rot[i]->x -= mol_centre.second.x();
               atoms_rot[i]->y -= mol_centre.second.y();
               atoms_rot[i]->z -= mol_centre.second.z();
            }
            for (int i=0; i<na_rot; i++) {
               glm::vec3 pos(atoms_rot[i]->x, atoms_rot[i]->y, atoms_rot[i]->z);
               glm::vec3 rotated = rot_mat * pos;
               atoms_rot[i]->x = rotated.x;
               atoms_rot[i]->y = rotated.y;
               atoms_rot[i]->z = rotated.z;
            }

            // Compute atom map on the local box cell/grid (optimisation 1)
            clipper::Xmap<float> rotated_model_map = coot::util::calc_atom_map(mol_rot, sel_rot,
                                                                                local_cell, spacegroup, local_gs);

            unsigned int n_trans_raw = n_translation_solutions * 20;
            std::vector<coot::translation_search_result_t> trans_results_raw =
               coot::phased_translation_search(fphi_obs_local, local_hkls, local_obs_map,
                                               rotated_model_map, n_trans_raw);

            // Convert local-box translations back to observed-map coordinates.
            // The local box has its origin at (recentred - local_cell_centre),
            // so a peak at local position p corresponds to observed-map position
            // p + (recentred - local_cell_centre).
            clipper::Coord_orth local_origin_in_obs(recentred.x() - local_cell_centre.x(),
                                                    recentred.y() - local_cell_centre.y(),
                                                    recentred.z() - local_cell_centre.z());

            // Filter: keep only peaks within fragment_radius of the recentred position
            float max_dist = fragment_radius;
            std::vector<coot::translation_search_result_t> trans_results;
            for (auto &tr : trans_results_raw) {
               // Convert to observed-map coordinates
               clipper::Coord_orth obs_pos(tr.position.x() + local_origin_in_obs.x(),
                                           tr.position.y() + local_origin_in_obs.y(),
                                           tr.position.z() + local_origin_in_obs.z());
               tr.position = obs_pos;

               double dx = obs_pos.x() - recentred.x();
               double dy = obs_pos.y() - recentred.y();
               double dz = obs_pos.z() - recentred.z();
               double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
               if (dist < max_dist) {
                  trans_results.push_back(tr);
                  if (static_cast<int>(trans_results.size()) >= n_translation_solutions)
                     break;
               }
            }

            // Create placed models for each translation peak
            auto &local_solutions = per_rotation_solutions[irot];
            for (unsigned int itrans=0; itrans<trans_results.size(); itrans++) {

               const clipper::Coord_orth &trans_pos = trans_results[itrans].position;

               mmdb::Manager *mol_placed = new mmdb::Manager();
               mol_placed->Copy(mol_rot, mmdb::MMDBFCM_All);

               int sel_placed = mol_placed->NewSelection();
               mol_placed->SelectAtoms(sel_placed, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                                       "*", "*", "*", "*");
               mmdb::PPAtom atoms_placed = 0;
               int na_placed;
               mol_placed->GetSelIndex(sel_placed, atoms_placed, na_placed);

               for (int i=0; i<na_placed; i++) {
                  atoms_placed[i]->x += trans_pos.x();
                  atoms_placed[i]->y += trans_pos.y();
                  atoms_placed[i]->z += trans_pos.z();
               }
               mol_placed->DeleteSelection(sel_placed);

               mr_solution_t sol;
               sol.rotation = rot_q;
               sol.rotation_score = rot_score;
               sol.translation = trans_pos;
               sol.translation_score = trans_results[itrans].score;
               sol.placed_mol = mol_placed;
               local_solutions.push_back(sol);
            }

            mol_rot->DeleteSelection(sel_rot);
            delete mol_rot;
         }
      });
   }
   for (auto &th : trans_threads) th.join();

   // Merge per-rotation results
   for (auto &rv : per_rotation_solutions)
      all_solutions.insert(all_solutions.end(), rv.begin(), rv.end());

   auto t_trans_end = std::chrono::high_resolution_clock::now();
   double t_trans_ms = std::chrono::duration<double, std::milli>(t_trans_end - t_trans_start).count();

   auto t_total_end = std::chrono::high_resolution_clock::now();
   double t_total_ms = std::chrono::duration<double, std::milli>(t_total_end - t_total_start).count();

   // Sort by translation score (descending)
   std::sort(all_solutions.begin(), all_solutions.end(),
             [](const mr_solution_t &a, const mr_solution_t &b) {
                return a.translation_score > b.translation_score;
             });

   // Filter: keep only translation solutions >= 70% of the top TF score
   if (! all_solutions.empty()) {
      float top_tf = all_solutions[0].translation_score;
      float tf_cutoff = top_tf * 0.7f;
      auto it = std::remove_if(all_solutions.begin(), all_solutions.end(),
                                [tf_cutoff](const mr_solution_t &s) {
                                   return s.translation_score < tf_cutoff;
                                });
      unsigned int n_before = all_solutions.size();
      // Delete placed_mol for filtered-out solutions
      for (auto jt = it; jt != all_solutions.end(); ++jt) {
         if (jt->placed_mol) {
            delete jt->placed_mol;
            jt->placed_mol = nullptr;
         }
      }
      all_solutions.erase(it, all_solutions.end());
      std::cout << "INFO:: molecular_replacement_search TF filter: "
                << n_before << " -> " << all_solutions.size()
                << " (cutoff 70% of " << std::fixed << std::setprecision(1)
                << top_tf << " sigma)" << std::endl;
   }

   if (truth.valid)
      for (unsigned int i=0; i<all_solutions.size(); i++)
         std::cout << "MR-DIAG:: pre-rigid-body solution " << i << " TF " << std::fixed << std::setprecision(2)
                   << all_solutions[i].translation_score << " rot " << std::setprecision(1)
                   << truth.angle_to(all_solutions[i].rotation) << " deg, trans "
                   << truth.distance_to(all_solutions[i].translation) << " A from true, rmsd "
                   << truth.rmsd(all_solutions[i].placed_mol) << " A" << std::endl;

   // --- Rigid-body refinement of surviving solutions (parallelised) ---

   auto t_rigid_start = std::chrono::high_resolution_clock::now();

   if (! all_solutions.empty()) {
      unsigned int n_sol = all_solutions.size();
      unsigned int n_rigid_threads = coot::get_max_number_of_threads();
      if (n_rigid_threads == 0) n_rigid_threads = 4;
      if (n_rigid_threads > n_sol) n_rigid_threads = n_sol;

      auto rigid_ranges = coot::atom_index_ranges(n_sol, n_rigid_threads);
      std::vector<std::thread> rigid_threads;
      rigid_threads.reserve(rigid_ranges.size());

      for (unsigned int t=0; t<rigid_ranges.size(); t++) {
         rigid_threads.emplace_back([&, t]() {
            for (unsigned int isol=rigid_ranges[t].first; isol<rigid_ranges[t].second; isol++) {
               mr_solution_t &sol = all_solutions[isol];
               if (! sol.placed_mol) continue;

               // For MR all atoms move; the masking molecule is empty.
               coot::minimol::molecule moving_mol(sol.placed_mol);
               coot::minimol::molecule empty_mol;
               coot::minimol::molecule moved_mol =
                  coot::rigid_body_fit_with_masking(empty_mol, moving_mol, xmap_obs);

               // Transfer refined coordinates back into placed_mol
               for (unsigned int ifrag=0; ifrag<moved_mol.fragments.size(); ifrag++) {
                  for (int ires=moved_mol[ifrag].min_res_no();
                       ires<=moved_mol[ifrag].max_residue_number(); ires++) {
                     for (unsigned int iat=0; iat<moved_mol[ifrag][ires].atoms.size(); iat++) {
                        const auto &moved_at = moved_mol[ifrag][ires][iat];
                        // Find the corresponding residue in placed_mol
                        coot::residue_spec_t spec(moved_mol[ifrag].fragment_id, ires, "");
                        mmdb::Residue *res = coot::util::get_residue(spec, sol.placed_mol);
                        if (res) {
                           mmdb::PPAtom res_atoms = nullptr;
                           int n_res_atoms = 0;
                           res->GetAtomTable(res_atoms, n_res_atoms);
                           for (int j=0; j<n_res_atoms; j++) {
                              if (! res_atoms[j]->isTer()) {
                                 std::string name(res_atoms[j]->GetAtomName());
                                 std::string altloc(res_atoms[j]->altLoc);
                                 if (name == moved_at.name && altloc == moved_at.altLoc) {
                                    res_atoms[j]->x = moved_at.pos.x();
                                    res_atoms[j]->y = moved_at.pos.y();
                                    res_atoms[j]->z = moved_at.pos.z();
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         });
      }
      for (auto &th : rigid_threads) th.join();

      std::cout << "INFO:: molecular_replacement_search rigid-body refined "
                << n_sol << " solutions" << std::endl;

      // Remove duplicate solutions that have converged to the same position.
      // Compare centres of mass of the placed models; keep only solutions
      // whose centre is at least 3 A from all higher-ranked (earlier) ones.
      float min_sep = 3.0f;
      float min_sep_sq = min_sep * min_sep;
      std::vector<mr_solution_t> unique_solutions;
      std::vector<clipper::Coord_orth> unique_centres;
      unique_solutions.reserve(all_solutions.size());
      unique_centres.reserve(all_solutions.size());

      for (auto &sol : all_solutions) {
         if (! sol.placed_mol) continue;
         std::pair<bool, clipper::Coord_orth> com = coot::centre_of_molecule(sol.placed_mol);
         if (! com.first) continue;

         bool too_close = false;
         for (const auto &uc : unique_centres) {
            double dx = com.second.x() - uc.x();
            double dy = com.second.y() - uc.y();
            double dz = com.second.z() - uc.z();
            if (dx*dx + dy*dy + dz*dz < min_sep_sq) {
               too_close = true;
               break;
            }
         }
         if (! too_close) {
            unique_centres.push_back(com.second);
            unique_solutions.push_back(sol);
         } else {
            // Clean up the duplicate
            if (sol.placed_mol) {
               delete sol.placed_mol;
               sol.placed_mol = nullptr;
            }
         }
      }

      if (unique_solutions.size() < all_solutions.size()) {
         std::cout << "INFO:: molecular_replacement_search post-rigid-body proximity filter: "
                   << all_solutions.size() << " -> " << unique_solutions.size()
                   << " (minimum separation " << min_sep << " A)" << std::endl;
      }
      all_solutions = unique_solutions;

      if (truth.valid)
         for (unsigned int i=0; i<all_solutions.size(); i++)
            std::cout << "MR-DIAG:: post-rigid-body solution " << i << " TF " << std::fixed
                      << std::setprecision(2) << all_solutions[i].translation_score << " rmsd "
                      << truth.rmsd(all_solutions[i].placed_mol) << " A from true" << std::endl;
   }

   auto t_rigid_end = std::chrono::high_resolution_clock::now();
   double t_rigid_ms = std::chrono::duration<double, std::milli>(t_rigid_end - t_rigid_start).count();

   // Recompute t_total_end after rigid-body refinement
   t_total_end = std::chrono::high_resolution_clock::now();
   t_total_ms = std::chrono::duration<double, std::milli>(t_total_end - t_total_start).count();

   // Print summary table
   std::cout << std::endl;
   std::cout << "INFO:: =================================================================" << std::endl;
   std::cout << "INFO:: Molecular Replacement Results" << std::endl;
   std::cout << "INFO:: Target centre: (" << std::fixed << std::setprecision(1)
             << target_centre.x() << ", " << target_centre.y() << ", " << target_centre.z()
             << ") search centre (recentred): ("
             << recentred.x() << ", " << recentred.y() << ", " << recentred.z()
             << ")" << std::endl;
   std::cout << "INFO:: =================================================================" << std::endl;
   std::cout << "INFO:: " << std::setw(4) << "Rank"
             << std::setw(10) << "Rot_scr"
             << std::setw(10) << "TF_sigma"
             << "   Position" << std::endl;
   std::cout << "INFO:: " << std::string(60, '-') << std::endl;

   unsigned int n_print = std::min(static_cast<unsigned int>(all_solutions.size()),
                                   static_cast<unsigned int>(n_rotation_solutions * n_translation_solutions));
   for (unsigned int i=0; i<n_print; i++) {
      const auto &s = all_solutions[i];
      std::cout << "INFO:: " << std::setw(4) << i
                << std::scientific << std::setprecision(2) << std::setw(10) << s.rotation_score
                << std::fixed << std::setprecision(1) << std::setw(10) << s.translation_score
                << "   (" << std::setprecision(1)
                << s.translation.x() << ", " << s.translation.y() << ", " << s.translation.z()
                << ")" << std::endl;
   }

   std::cout << "INFO:: =================================================================" << std::endl;
   std::cout << "INFO:: Timings: recentring " << std::fixed << std::setprecision(1) << t_recentre_ms << " ms"
             << ", rotation " << t_rot_ms << " ms"
             << ", refinement " << t_refine_ms << " ms"
             << ", translation " << t_trans_ms << " ms"
             << ", rigid-body " << t_rigid_ms << " ms"
             << ", total " << t_total_ms << " ms (" << t_total_ms / 1000.0 << " s)" << std::endl;
   std::cout << "INFO:: =================================================================" << std::endl;

   // Clean up the model working copy
   mol->DeleteSelection(SelHnd);
   delete mol;

   return all_solutions;
}
