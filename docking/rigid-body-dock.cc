/*
 * docking/rigid-body-dock.cc
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

#include <glm/gtc/matrix_transform.hpp>

#include "rigid-body-dock.hh"
#include "utils/split-indices.hh"

using namespace coot::haddock;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

rigid_body_docker_t::rigid_body_docker_t(const molecule_data_t &mol_A,
                                         const molecule_data_t &mol_B,
                                         const docking_parameters_t &params)
   : mol_A_(mol_A), mol_B_(mol_B), params_(params) {

   com_A_ = centre_of_mass(mol_A_.atoms);
   com_B_ = centre_of_mass(mol_B_.atoms);

   // Build spatial grid for protein A (fixed throughout docking).
   // Use the larger of VDW and electrostatic cutoffs as the cell size.
   double cutoff = std::max(params_.vdw_cutoff, params_.elec_cutoff);
   grid_A_ = spatial_grid_t(mol_A_.atoms, cutoff);

   // Auto-compute separation distance if not explicitly set:
   // 2 × max(radius_A, radius_B) + 15 A (following HADDOCK convention)
   if (params_.separation_distance <= 0) {
      double max_r_A = 0;
      for (const auto &at : mol_A_.atoms) {
         double r = clipper::Coord_orth::length(at.position, com_A_);
         if (r > max_r_A) max_r_A = r;
      }
      double max_r_B = 0;
      for (const auto &at : mol_B_.atoms) {
         double r = clipper::Coord_orth::length(at.position, com_B_);
         if (r > max_r_B) max_r_B = r;
      }
      params_.separation_distance = max_r_A + max_r_B + 15.0;
      std::cout << "HADDOCK: auto separation = " << params_.separation_distance
                << " A (r_A=" << max_r_A << " r_B=" << max_r_B << ")"
                << std::endl;
   }

   // Centre protein B atoms at origin for rotation
   centred_B_atoms_ = mol_B_.atoms;
   for (auto &at : centred_B_atoms_) {
      at.position = clipper::Coord_orth(at.position.x() - com_B_.x(),
                                        at.position.y() - com_B_.y(),
                                        at.position.z() - com_B_.z());
   }

   // Centre the AIR partner residue atoms of B at origin too
   centred_B_partner_atoms_ = mol_B_.partner_residue_atoms;
   for (auto &ra : centred_B_partner_atoms_) {
      for (auto &pos : ra.positions) {
         pos = clipper::Coord_orth(pos.x() - com_B_.x(),
                                   pos.y() - com_B_.y(),
                                   pos.z() - com_B_.z());
      }
   }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

clipper::Coord_orth
rigid_body_docker_t::centre_of_mass(const std::vector<atom_data_t> &atoms) {

   double sx = 0, sy = 0, sz = 0;
   for (const auto &at : atoms) {
      sx += at.position.x();
      sy += at.position.y();
      sz += at.position.z();
   }
   double n = static_cast<double>(atoms.size());
   return clipper::Coord_orth(sx/n, sy/n, sz/n);
}

glm::quat
rigid_body_docker_t::random_quaternion(std::mt19937 &rng) {

   // Uniform random quaternion using Marsaglia's method:
   // generate two pairs of uniform random numbers on the unit disk,
   // then combine.
   std::uniform_real_distribution<double> dist(-1.0, 1.0);

   double s1, s2, x1, y1, x2, y2;
   do {
      x1 = dist(rng);
      y1 = dist(rng);
      s1 = x1*x1 + y1*y1;
   } while (s1 >= 1.0);

   do {
      x2 = dist(rng);
      y2 = dist(rng);
      s2 = x2*x2 + y2*y2;
   } while (s2 >= 1.0);

   double root = std::sqrt((1.0 - s1) / s2);
   return glm::quat(static_cast<float>(x1),
                     static_cast<float>(y1),
                     static_cast<float>(x2 * root),
                     static_cast<float>(y2 * root));
}

std::vector<clipper::Coord_orth>
rigid_body_docker_t::transform_positions(const std::vector<atom_data_t> &atoms,
                                         const glm::quat &rotation,
                                         const clipper::Coord_orth &translation) {

   glm::mat3 R = glm::mat3_cast(rotation);
   std::vector<clipper::Coord_orth> result(atoms.size());
   for (unsigned int i=0; i<atoms.size(); i++) {
      const clipper::Coord_orth &p = atoms[i].position;
      glm::vec3 v(static_cast<float>(p.x()),
                   static_cast<float>(p.y()),
                   static_cast<float>(p.z()));
      glm::vec3 rv = R * v;
      result[i] = clipper::Coord_orth(rv.x + translation.x(),
                                       rv.y + translation.y(),
                                       rv.z + translation.z());
   }
   return result;
}

std::vector<residue_atoms_t>
rigid_body_docker_t::transform_residue_atoms(const std::vector<residue_atoms_t> &res_atoms,
                                             const glm::quat &rotation,
                                             const clipper::Coord_orth &translation) {

   glm::mat3 R = glm::mat3_cast(rotation);
   std::vector<residue_atoms_t> result(res_atoms.size());
   for (unsigned int k=0; k<res_atoms.size(); k++) {
      result[k].positions.resize(res_atoms[k].positions.size());
      for (unsigned int n=0; n<res_atoms[k].positions.size(); n++) {
         const clipper::Coord_orth &p = res_atoms[k].positions[n];
         glm::vec3 v(static_cast<float>(p.x()),
                      static_cast<float>(p.y()),
                      static_cast<float>(p.z()));
         glm::vec3 rv = R * v;
         result[k].positions[n] = clipper::Coord_orth(rv.x + translation.x(),
                                                       rv.y + translation.y(),
                                                       rv.z + translation.z());
      }
   }
   return result;
}

// ---------------------------------------------------------------------------
// Euler angle / quaternion conversion (ZYZ convention)
// ---------------------------------------------------------------------------

glm::quat
rigid_body_docker_t::euler_to_quat(double alpha, double beta, double gamma) {

   float fa = static_cast<float>(alpha);
   float fb = static_cast<float>(beta);
   float fg = static_cast<float>(gamma);
   glm::quat qza = glm::angleAxis(fa, glm::vec3(0, 0, 1));
   glm::quat qyb = glm::angleAxis(fb, glm::vec3(0, 1, 0));
   glm::quat qzg = glm::angleAxis(fg, glm::vec3(0, 0, 1));
   return qza * qyb * qzg;
}

void
rigid_body_docker_t::quat_to_euler(const glm::quat &q,
                                   double *alpha_p, double *beta_p, double *gamma_p) {

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
// Energy evaluation for a given pose
// ---------------------------------------------------------------------------

intermolecular_energy_t
rigid_body_docker_t::evaluate_energy(const glm::quat &rotation,
                                     const clipper::Coord_orth &translation) const {

   return evaluate_with_rb_gradient(rotation, translation,
                                    0.0, 0.0, false, nullptr);
}

// ---------------------------------------------------------------------------
// Energy + analytical rigid body gradient
// ---------------------------------------------------------------------------
//
// The position of atom j of protein B (centred at origin) is:
//
//   x_j = R(alpha, beta, gamma) * x_j^0  +  t
//
// where R = R_z(alpha) * R_y(beta) * R_z(gamma) is the ZYZ Euler rotation.
//
// Translation gradient:
//   dE/dt_x = sum_j  dE/dx_j
//   (and similarly for y, z)
//
// Rotational gradients use the fact that differentiating a rotation by
// angle theta about axis n gives  dx_j/dtheta = n x (x_j - t).
// For ZYZ Euler angles the three axes are:
//   alpha:  z-axis                     = (0, 0, 1)
//   beta:   R_z(alpha) * y-axis        = (-sin(alpha), cos(alpha), 0)
//   gamma:  R_z(alpha)*R_y(beta)*z-axis = (sin(beta)*cos(alpha), sin(beta)*sin(alpha), cos(beta))
//
// So dE/d(angle_k) = sum_j g_j . (axis_k x r_j)
//                   = sum_j axis_k . (r_j x g_j)    [scalar triple product]
//
// where g_j = dE/dx_j is the per-atom force and r_j = x_j - t.
//

intermolecular_energy_t
rigid_body_docker_t::evaluate_with_rb_gradient(const glm::quat &rotation,
                                               const clipper::Coord_orth &translation,
                                               double alpha, double beta,
                                               bool rotation_only,
                                               double *gradient_out) const {

   // Transform protein B atoms
   glm::mat3 R = glm::mat3_cast(rotation);
   std::vector<atom_data_t> transformed_B = centred_B_atoms_;
   for (auto &at : transformed_B) {
      glm::vec3 p(static_cast<float>(at.position.x()),
                   static_cast<float>(at.position.y()),
                   static_cast<float>(at.position.z()));
      glm::vec3 rp = R * p;
      at.position = clipper::Coord_orth(rp.x + translation.x(),
                                         rp.y + translation.y(),
                                         rp.z + translation.z());
   }

   // Transform AIR partner B residue atoms
   std::vector<residue_atoms_t> transformed_partner_B =
      transform_residue_atoms(centred_B_partner_atoms_, rotation, translation);

   // Compute VDW + Coulomb energy (and per-atom gradients on B if needed)
   std::vector<clipper::Coord_orth> grad_B;
   std::vector<clipper::Coord_orth> *grad_B_p = nullptr;
   if (gradient_out) {
      grad_B.resize(transformed_B.size(), clipper::Coord_orth(0, 0, 0));
      grad_B_p = &grad_B;
   }

   intermolecular_energy_t e =
      compute_intermolecular_energy(mol_A_.atoms, transformed_B,
                                    mol_A_.active_residue_atoms,
                                    transformed_partner_B,
                                    params_, grad_B_p, &grid_A_);

   if (!gradient_out) return e;

   // Compute AIR gradients on partner B residue atoms
   std::vector<residue_atoms_t> air_grad_B(transformed_partner_B.size());
   for (unsigned int k=0; k<transformed_partner_B.size(); k++) {
      air_grad_B[k].positions.resize(transformed_partner_B[k].positions.size(),
                                      clipper::Coord_orth(0, 0, 0));
   }
   air_energy_and_gradients(mol_A_.active_residue_atoms,
                            transformed_partner_B,
                            3.0, nullptr, &air_grad_B);

   // Compute the three ZYZ rotation axes
   double sa = std::sin(alpha), ca = std::cos(alpha);
   double sb = std::sin(beta),  cb = std::cos(beta);

   // axis_alpha = z = (0, 0, 1)
   // axis_beta  = R_z(alpha) * y = (-sin(alpha), cos(alpha), 0)
   // axis_gamma = R_z(alpha) * R_y(beta) * z = (sb*ca, sb*sa, cb)
   double ax_a[3] = {0.0, 0.0, 1.0};
   double ax_b[3] = {-sa, ca, 0.0};
   double ax_g[3] = {sb*ca, sb*sa, cb};

   double dE_dalpha = 0, dE_dbeta = 0, dE_dgamma = 0;
   double dE_dtx = 0, dE_dty = 0, dE_dtz = 0;

   // Project VDW + Coulomb per-atom gradients onto 6 DOF
   for (unsigned int j=0; j<transformed_B.size(); j++) {
      double gx = grad_B[j].x();
      double gy = grad_B[j].y();
      double gz = grad_B[j].z();

      // r_j = transformed position - translation = R * x_j^0
      double rx = transformed_B[j].position.x() - translation.x();
      double ry = transformed_B[j].position.y() - translation.y();
      double rz = transformed_B[j].position.z() - translation.z();

      // torque_j = r_j x g_j
      double tqx = ry*gz - rz*gy;
      double tqy = rz*gx - rx*gz;
      double tqz = rx*gy - ry*gx;

      dE_dalpha += ax_a[0]*tqx + ax_a[1]*tqy + ax_a[2]*tqz;
      dE_dbeta  += ax_b[0]*tqx + ax_b[1]*tqy + ax_b[2]*tqz;
      dE_dgamma += ax_g[0]*tqx + ax_g[1]*tqy + ax_g[2]*tqz;

      dE_dtx += gx;
      dE_dty += gy;
      dE_dtz += gz;
   }

   // Project AIR gradients (on partner B residue atoms) onto 6 DOF
   double air_w = params_.air_weight;
   for (unsigned int k=0; k<air_grad_B.size(); k++) {
      for (unsigned int n=0; n<air_grad_B[k].positions.size(); n++) {
         double gx = air_grad_B[k].positions[n].x();
         double gy = air_grad_B[k].positions[n].y();
         double gz = air_grad_B[k].positions[n].z();

         double rx = transformed_partner_B[k].positions[n].x() - translation.x();
         double ry = transformed_partner_B[k].positions[n].y() - translation.y();
         double rz = transformed_partner_B[k].positions[n].z() - translation.z();

         double tqx = ry*gz - rz*gy;
         double tqy = rz*gx - rx*gz;
         double tqz = rx*gy - ry*gx;

         dE_dalpha += air_w * (ax_a[0]*tqx + ax_a[1]*tqy + ax_a[2]*tqz);
         dE_dbeta  += air_w * (ax_b[0]*tqx + ax_b[1]*tqy + ax_b[2]*tqz);
         dE_dgamma += air_w * (ax_g[0]*tqx + ax_g[1]*tqy + ax_g[2]*tqz);

         dE_dtx += air_w * gx;
         dE_dty += air_w * gy;
         dE_dtz += air_w * gz;
      }
   }

   gradient_out[0] = dE_dalpha;
   gradient_out[1] = dE_dbeta;
   gradient_out[2] = dE_dgamma;
   if (!rotation_only) {
      gradient_out[3] = dE_dtx;
      gradient_out[4] = dE_dty;
      gradient_out[5] = dE_dtz;
   }

   return e;
}

// ---------------------------------------------------------------------------
// GSL minimiser callbacks
// ---------------------------------------------------------------------------

double
rigid_body_docker_t::gsl_f(const gsl_vector *v, void *params) {

   auto *ctx = static_cast<minimiser_context_t *>(params);
   const rigid_body_docker_t *docker = ctx->docker;

   double alpha = gsl_vector_get(v, 0);
   double beta  = gsl_vector_get(v, 1);
   double gamma = gsl_vector_get(v, 2);
   glm::quat rotation = euler_to_quat(alpha, beta, gamma);

   clipper::Coord_orth translation(0, 0, 0);
   if (!ctx->rotation_only) {
      translation = clipper::Coord_orth(gsl_vector_get(v, 3),
                                         gsl_vector_get(v, 4),
                                         gsl_vector_get(v, 5));
   }

   intermolecular_energy_t e = docker->evaluate_with_rb_gradient(
      rotation, translation, alpha, beta, ctx->rotation_only, nullptr);
   return e.e_total;
}

void
rigid_body_docker_t::gsl_df(const gsl_vector *v, void *params, gsl_vector *df) {

   auto *ctx = static_cast<minimiser_context_t *>(params);
   const rigid_body_docker_t *docker = ctx->docker;

   double alpha = gsl_vector_get(v, 0);
   double beta  = gsl_vector_get(v, 1);
   double gamma = gsl_vector_get(v, 2);
   glm::quat rotation = euler_to_quat(alpha, beta, gamma);

   clipper::Coord_orth translation(0, 0, 0);
   if (!ctx->rotation_only) {
      translation = clipper::Coord_orth(gsl_vector_get(v, 3),
                                         gsl_vector_get(v, 4),
                                         gsl_vector_get(v, 5));
   }

   int n_var = ctx->rotation_only ? 3 : 6;
   double gradient[6] = {0};
   docker->evaluate_with_rb_gradient(rotation, translation,
                                     alpha, beta,
                                     ctx->rotation_only, gradient);

   for (int i=0; i<n_var; i++)
      gsl_vector_set(df, i, gradient[i]);
}

void
rigid_body_docker_t::gsl_fdf(const gsl_vector *v, void *params,
                              double *f, gsl_vector *df) {

   auto *ctx = static_cast<minimiser_context_t *>(params);
   const rigid_body_docker_t *docker = ctx->docker;

   double alpha = gsl_vector_get(v, 0);
   double beta  = gsl_vector_get(v, 1);
   double gamma = gsl_vector_get(v, 2);
   glm::quat rotation = euler_to_quat(alpha, beta, gamma);

   clipper::Coord_orth translation(0, 0, 0);
   if (!ctx->rotation_only) {
      translation = clipper::Coord_orth(gsl_vector_get(v, 3),
                                         gsl_vector_get(v, 4),
                                         gsl_vector_get(v, 5));
   }

   int n_var = ctx->rotation_only ? 3 : 6;
   double gradient[6] = {0};
   intermolecular_energy_t e = docker->evaluate_with_rb_gradient(
      rotation, translation, alpha, beta, ctx->rotation_only, gradient);

   *f = e.e_total;
   for (int i=0; i<n_var; i++)
      gsl_vector_set(df, i, gradient[i]);
}

// ---------------------------------------------------------------------------
// Rigid body minimisation
// ---------------------------------------------------------------------------

void
rigid_body_docker_t::rigid_body_minimise(glm::quat *rotation_p,
                                         clipper::Coord_orth *translation_p,
                                         bool rotation_only) {

   int n_var = rotation_only ? 3 : 6;

   gsl_multimin_function_fdf my_func;
   my_func.n = n_var;
   my_func.f = &gsl_f;
   my_func.df = &gsl_df;
   my_func.fdf = &gsl_fdf;

   minimiser_context_t ctx;
   ctx.docker = this;
   ctx.rotation_only = rotation_only;
   my_func.params = static_cast<void *>(&ctx);

   // Convert current quaternion to Euler angles
   double alpha, beta, gamma;
   quat_to_euler(*rotation_p, &alpha, &beta, &gamma);

   gsl_vector *x = gsl_vector_alloc(n_var);
   gsl_vector_set(x, 0, alpha);
   gsl_vector_set(x, 1, beta);
   gsl_vector_set(x, 2, gamma);
   if (!rotation_only) {
      gsl_vector_set(x, 3, translation_p->x());
      gsl_vector_set(x, 4, translation_p->y());
      gsl_vector_set(x, 5, translation_p->z());
   }

   const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;
   gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, n_var);

   // For rotation-only, 0.1 radians is a good step; for full 6-DOF,
   // the translation needs larger steps to cover the separation distance.
   double step_size = rotation_only ? 0.1 : 1.0;
   gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, 1e-3);

   unsigned int iter = 0;
   unsigned int max_iter = rotation_only ? 100 : 500;
   int status;

   do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);
      if (status) break;
      status = gsl_multimin_test_gradient(s->gradient, 1e-2);
   } while (status == GSL_CONTINUE && iter < max_iter);

   // Read back the optimised parameters
   alpha = gsl_vector_get(s->x, 0);
   beta  = gsl_vector_get(s->x, 1);
   gamma = gsl_vector_get(s->x, 2);
   *rotation_p = euler_to_quat(alpha, beta, gamma);

   if (!rotation_only) {
      *translation_p = clipper::Coord_orth(gsl_vector_get(s->x, 3),
                                            gsl_vector_get(s->x, 4),
                                            gsl_vector_get(s->x, 5));
   }

   gsl_multimin_fdfminimizer_free(s);
   gsl_vector_free(x);
}

// ---------------------------------------------------------------------------
// Single trial
// ---------------------------------------------------------------------------

docking_result_t
rigid_body_docker_t::run_single_trial(std::mt19937 &rng) {

   // Step 1: Random rotation for protein A (applied to A's atoms
   // temporarily) and random rotation for protein B.
   //
   // Actually, for the rigid body stage we fix protein A and only
   // vary protein B's pose.  But we randomise the initial orientation
   // of B relative to A.  Randomising A would be equivalent to
   // randomising B with a different random rotation, so we only need
   // to randomise B.

   glm::quat rotation = random_quaternion(rng);

   // Step 2: Initial translation — place protein B at separation_distance
   // from protein A's centre of mass in a random direction.
   // Use a uniform random point on the sphere (Marsaglia's method)
   // so we sample all approach directions equally.
   std::uniform_real_distribution<double> uniform(-1.0, 1.0);
   double u1, u2, s;
   do {
      u1 = uniform(rng);
      u2 = uniform(rng);
      s = u1*u1 + u2*u2;
   } while (s >= 1.0);
   double factor = 2.0 * std::sqrt(1.0 - s);
   double nx = u1 * factor;
   double ny = u2 * factor;
   double nz = 1.0 - 2.0 * s;

   clipper::Coord_orth translation(
      com_A_.x() + params_.separation_distance * nx,
      com_A_.y() + params_.separation_distance * ny,
      com_A_.z() + params_.separation_distance * nz);

   // Step 3: Orientational optimisation (rotation only)
   for (int cycle=0; cycle<params_.n_orient_cycles; cycle++) {
      rigid_body_minimise(&rotation, &translation, true); // rotation only
   }

   // Step 4: Full rigid body minimisation (rotation + translation)
   rigid_body_minimise(&rotation, &translation, false);

   // Step 5: Second round of full minimisation to refine further
   rigid_body_minimise(&rotation, &translation, false);

   // Step 6: Evaluate final energy
   intermolecular_energy_t e = evaluate_energy(rotation, translation);

   docking_result_t result;
   result.rotation = rotation;
   result.translation = translation;
   result.e_vdw = e.e_vdw;
   result.e_elec = e.e_elec;
   result.e_air = e.e_air;
   result.e_total = e.e_total;

   return result;
}

// ---------------------------------------------------------------------------
// Main docking protocol
// ---------------------------------------------------------------------------

std::vector<docking_result_t>
rigid_body_docker_t::dock() {

   unsigned int n_threads = params_.n_threads;
   if (n_threads == 0) {
      n_threads = std::thread::hardware_concurrency();
      if (n_threads == 0) n_threads = 4;
   }
   if (n_threads > static_cast<unsigned int>(params_.n_trials))
      n_threads = params_.n_trials;

   std::cout << "HADDOCK rigid body docking: "
             << params_.n_trials << " trials on "
             << n_threads << " threads" << std::endl;
   std::cout << "  Protein A: " << mol_A_.atoms.size() << " atoms, "
             << mol_A_.active_residue_atoms.size() << " active residues" << std::endl;
   std::cout << "  Protein B: " << mol_B_.atoms.size() << " atoms, "
             << centred_B_partner_atoms_.size() << " partner residues" << std::endl;

   auto ranges = coot::atom_index_ranges(params_.n_trials, n_threads);

   std::vector<std::vector<docking_result_t>> thread_results(ranges.size());
   std::vector<std::thread> threads;

   for (unsigned int t=0; t<ranges.size(); t++) {
      threads.emplace_back([&, t]() {
         // Each thread gets its own RNG seeded differently
         std::mt19937 rng(42 + t * 1000 + ranges[t].first);
         auto &local_results = thread_results[t];
         local_results.reserve(ranges[t].second - ranges[t].first);

         for (unsigned int trial=ranges[t].first; trial<ranges[t].second; trial++) {
            docking_result_t result = run_single_trial(rng);
            local_results.push_back(result);

            if (trial % 100 == 0) {
               std::cout << "  Trial " << trial << "/" << params_.n_trials
                         << " E=" << std::setprecision(1) << std::fixed
                         << result.e_total << std::endl;
            }
         }
      });
   }

   for (auto &th : threads) th.join();

   // Merge results
   std::vector<docking_result_t> all_results;
   all_results.reserve(params_.n_trials);
   for (auto &tr : thread_results)
      all_results.insert(all_results.end(), tr.begin(), tr.end());

   // Sort by total energy (ascending — lowest energy first)
   std::sort(all_results.begin(), all_results.end(),
             [](const docking_result_t &a, const docking_result_t &b) {
                return a.e_total < b.e_total;
             });

   // Keep best n_keep
   if (static_cast<int>(all_results.size()) > params_.n_keep)
      all_results.resize(params_.n_keep);

   std::cout << "HADDOCK rigid body docking complete." << std::endl;
   if (!all_results.empty()) {
      std::cout << "  Best energy: " << std::setprecision(2) << std::fixed
                << all_results.front().e_total
                << " (VDW=" << all_results.front().e_vdw
                << " Elec=" << all_results.front().e_elec
                << " AIR=" << all_results.front().e_air << ")" << std::endl;
      std::cout << "  Worst kept:  " << all_results.back().e_total << std::endl;
   }

   return all_results;
}
