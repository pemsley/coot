/*
 * coot-utils/surface-on-torus.cc
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <cmath>
#include "surface-on-torus.hh"

glm::vec4
coot::height_to_colour(float h, float height_scale) {

   // map h into [0,1] using the scale factor as the expected peak value
   float t = h * height_scale;
   if (t < 0.0f) t = 0.0f;
   if (t > 1.0f) t = 1.0f;

   // dark blue (0) -> muted green (0.5) -> red (1)
   float r, g, b;
   if (t < 0.5f) {
      float s = t * 2.0f; // 0 to 1 over the first half
      r = 0.0f;
      g = 0.55f * s;      // up to muted green
      b = 0.3f * (1.0f - s) + 0.0f;  // dark blue fades out
   } else {
      float s = (t - 0.5f) * 2.0f; // 0 to 1 over the second half
      r = 0.8f * s;        // red rises
      g = 0.55f * (1.0f - s); // green fades
      b = 0.0f;
   }

   return glm::vec4(r, g, b, 1.0f);
}

coot::simple_mesh_t
coot::make_surface_on_torus(const std::vector<std::vector<float>> &height_data,
                            float R, float r, float height_scale) {

   const float pi = 3.14159265358979323846f;

   unsigned int n_psi = height_data.size();
   if (n_psi == 0) return simple_mesh_t();
   unsigned int n_phi = height_data[0].size();
   if (n_phi == 0) return simple_mesh_t();

   float d_phi = 2.0f * pi / static_cast<float>(n_phi);
   float d_psi = 2.0f * pi / static_cast<float>(n_psi);

   std::vector<api::vnc_vertex> vertices(n_psi * n_phi);
   std::vector<g_triangle> triangles;
   triangles.reserve(2 * n_psi * n_phi);

   // vertex positions and normals
   for (unsigned int ip=0; ip<n_psi; ip++) {
      float psi = -pi + (static_cast<float>(ip) + 0.5f) * d_psi;
      float cos_psi = cosf(psi);
      float sin_psi = sinf(psi);

      for (unsigned int jp=0; jp<n_phi; jp++) {
         float phi = -pi + (static_cast<float>(jp) + 0.5f) * d_phi;
         float cos_phi = cosf(phi);
         float sin_phi = sinf(phi);

         float h = height_data[ip][jp] * height_scale;

         // central differences with periodic wrapping
         unsigned int ip_next = (ip + 1) % n_psi;
         unsigned int ip_prev = (ip + n_psi - 1) % n_psi;
         unsigned int jp_next = (jp + 1) % n_phi;
         unsigned int jp_prev = (jp + n_phi - 1) % n_phi;

         float dh_dphi = height_scale * (height_data[ip][jp_next] - height_data[ip][jp_prev]) / (2.0f * d_phi);
         float dh_dpsi = height_scale * (height_data[ip_next][jp] - height_data[ip_prev][jp]) / (2.0f * d_psi);

         // displaced position
         float r_eff = r + h;
         float R_eff = R + r_eff * cos_psi;

         glm::vec3 pos(R_eff * cos_phi,
                       R_eff * sin_phi,
                       r_eff * sin_psi);

         // tangent vectors
         float Tx_phi = -R_eff * sin_phi + dh_dphi * cos_psi * cos_phi;
         float Ty_phi =  R_eff * cos_phi + dh_dphi * cos_psi * sin_phi;
         float Tz_phi =  dh_dphi * sin_psi;

         float Tx_psi = -r_eff * sin_psi * cos_phi + dh_dpsi * cos_psi * cos_phi;
         float Ty_psi = -r_eff * sin_psi * sin_phi + dh_dpsi * cos_psi * sin_phi;
         float Tz_psi =  r_eff * cos_psi            + dh_dpsi * sin_psi;

         // cross product for normal
         float Nx = Ty_phi * Tz_psi - Tz_phi * Ty_psi;
         float Ny = Tz_phi * Tx_psi - Tx_phi * Tz_psi;
         float Nz = Tx_phi * Ty_psi - Ty_phi * Tx_psi;

         float mag = sqrtf(Nx * Nx + Ny * Ny + Nz * Nz);
         if (mag < 1e-12f) mag = 1.0f;
         glm::vec3 normal(Nx / mag, Ny / mag, Nz / mag);

         unsigned int idx = ip * n_phi + jp;
         glm::vec4 col = height_to_colour(height_data[ip][jp], height_scale);
         vertices[idx] = api::vnc_vertex(pos, normal, col);
      }
   }

   // triangulate the grid with periodic wrapping on both axes
   for (unsigned int ip=0; ip<n_psi; ip++) {
      unsigned int ip_next = (ip + 1) % n_psi;
      unsigned int row_this = ip      * n_phi;
      unsigned int row_next = ip_next * n_phi;
      for (unsigned int jp=0; jp<n_phi; jp++) {
         unsigned int jp_next = (jp + 1) % n_phi;

         unsigned int idx_00 = row_this + jp;
         unsigned int idx_01 = row_this + jp_next;
         unsigned int idx_10 = row_next + jp;
         unsigned int idx_11 = row_next + jp_next;

         triangles.push_back(g_triangle(idx_00, idx_10, idx_11));
         triangles.push_back(g_triangle(idx_00, idx_11, idx_01));
      }
   }

   return simple_mesh_t(vertices, triangles);
}
