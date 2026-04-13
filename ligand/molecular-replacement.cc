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

#include <clipper/core/xmap.h>
#include <clipper/core/map_utils.h>
#include <clipper/core/map_interp.h>
#include <clipper/contrib/edcalc.h>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#include "molecular-replacement.hh"
#include "rigid-body.hh"
#include "utils/split-indices.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/crowther.hh"
#include "coot-utils/patterson.hh"

std::vector<coot::mr_solution_t>
coot::molecular_replacement_search(const clipper::Xmap<float> &xmap_obs,
                                   mmdb::Manager *mol_model,
                                   const clipper::Coord_orth &target_centre,
                                   int n_rotation_solutions,
                                   int n_translation_solutions) {

   std::vector<mr_solution_t> all_solutions;

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

   // --- Translation search for each refined rotation (parallelised) ---
   //
   // Two key optimisations:
   //
   // 1. Extract a local P1 box from the observed map around target_centre
   //    (side ~ 2*fragment_radius + border). Both the observed and model
   //    density are computed on this smaller grid, so the FFTs are much
   //    cheaper than using the full EM cell.
   //
   // 2. Pre-compute the FFT of the observed local box once, rather than
   //    recomputing it for every rotation.

   auto t_trans_start = std::chrono::high_resolution_clock::now();

   // Build a local P1 box centred on target_centre
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
         clipper::Coord_orth obs_pos(local_pos.x() - local_cell_centre.x() + target_centre.x(),
                                    local_pos.y() - local_cell_centre.y() + target_centre.y(),
                                    local_pos.z() - local_cell_centre.z() + target_centre.z());
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

   unsigned int n_rot = refined.size();
   unsigned int n_trans_threads = std::thread::hardware_concurrency();
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
            // The local box has its origin at (target_centre - local_cell_centre),
            // so a peak at local position p corresponds to observed-map position
            // p + (target_centre - local_cell_centre).
            clipper::Coord_orth local_origin_in_obs(target_centre.x() - local_cell_centre.x(),
                                                    target_centre.y() - local_cell_centre.y(),
                                                    target_centre.z() - local_cell_centre.z());

            // Filter: keep only peaks within fragment_radius of the target centre
            float max_dist = fragment_radius;
            std::vector<coot::translation_search_result_t> trans_results;
            for (auto &tr : trans_results_raw) {
               // Convert to observed-map coordinates
               clipper::Coord_orth obs_pos(tr.position.x() + local_origin_in_obs.x(),
                                           tr.position.y() + local_origin_in_obs.y(),
                                           tr.position.z() + local_origin_in_obs.z());
               tr.position = obs_pos;

               double dx = obs_pos.x() - target_centre.x();
               double dy = obs_pos.y() - target_centre.y();
               double dz = obs_pos.z() - target_centre.z();
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

   // --- Rigid-body refinement of surviving solutions (parallelised) ---

   auto t_rigid_start = std::chrono::high_resolution_clock::now();

   if (! all_solutions.empty()) {
      unsigned int n_sol = all_solutions.size();
      unsigned int n_rigid_threads = std::thread::hardware_concurrency();
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
