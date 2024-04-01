/*
 * coot-utils/cylinder.hh
 *
 * Copyright 2023 by Medical Research Council
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

#ifndef COOT_CYLINDER_HH
#define COOT_CYLINDER_HH

#include <vector>
#include "vertex.hh"
#include "coot-utils/g_triangle.hh"

class cylinder {
   void add_flat_cap(float z);
   float height;
   float base_radius;
   float top_radius;
   unsigned int n_slices;
   glm::mat4 ori;
   glm::vec3 start;
   glm::vec4 basic_colour;
   float unstubby_rounded_cap_factor;
   // add these new triangles to the mesh
   void add_vertices_and_triangles(const std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > &vt);
   void init(const std::pair<glm::vec3, glm::vec3> &cart_pair,
             float base_radius, float top_radius, float height,
             const glm::vec4 &base_colour,
             unsigned int n_slices=8, unsigned int n_stacks=2);
   void add_flat_cap(int end_type); // 0 for base, 1 for top
public:
   std::vector<coot::api::vnc_vertex> vertices;
   std::vector<g_triangle> triangles;
   cylinder() {
      height = 1.0; base_radius = 1.0; top_radius = 1.0; n_slices = 16;
      unstubby_rounded_cap_factor = 1.0; } // needs to be filled
   cylinder(const std::pair<glm::vec3, glm::vec3> &pos_pair,
            float base_radius, float top_radius, float height,
            unsigned int n_slices=8, unsigned int n_stacks=2);
   cylinder(const std::pair<glm::vec3, glm::vec3> &pos_pair,
            float base_radius, float top_radius, float height,
            const glm::vec4 &base_colour,
            unsigned int n_slices=8, unsigned int n_stacks=2);
   void add_flat_end_cap();
   void add_flat_start_cap();
   void add_octahemisphere_end_cap();
   void add_octahemisphere_start_cap();
   void set_unstubby_rounded_cap_factor(float f) { unstubby_rounded_cap_factor = f; }
   void z_translate(float zt);
   void add_sad_face(); // to move the eyesballs. the eyes would have to be in a different mesh.
   void crenulations();
};

#endif
