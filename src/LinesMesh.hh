/*
 * src/LinesMesh.hh
 *
 * Copyright 2020 by Medical Research Council
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

#ifndef LINES_MESH_HH
#define LINES_MESH_HH

#include <vector>
#include <clipper/core/coords.h>

#include "generic-vertex.hh"
#include "Shader.hh"

#define COMPILED_WITH_CLIPPER

// Currently only ambient-lit - e.g. unit cell

class LinesMesh {
   enum { VAO_NOT_SET = 99999999 };
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   bool first_time;
   void init();
   void make_vertices_for_pulse(const glm::vec4 &colour, float radius,
                                unsigned int n_rings,
                                float theta_offset, bool broken_mode);
   glm::vec3 central_position;
   std::string name;
   bool offset_positions_have_been_set; // because we need to know if we should send the position_offset to
                                 // the shader as a unifrom
   bool scales_have_been_set; // because we need to know if we should send the scales to
                                 // the shader as a uniform
   glm::vec2 offset_positions; // zero by default
   glm::vec2 scales; // 1.0 by default

public:
   LinesMesh() { init(); }
   // e.g. a box will have 8 vertices and 12 * 2 indices
   LinesMesh(const std::vector<s_generic_vertex> &vertices_in,
             const std::vector<unsigned int> &indices_in) : vertices(vertices_in), indices(indices_in) {
      init();
   }
#ifdef COMPILED_WITH_CLIPPER
   explicit LinesMesh(const clipper::Cell &cell);
#endif
   std::vector<s_generic_vertex> vertices;
   std::vector<unsigned int> indices;
   void set_name(const std::string &n) { name = n; }
   void setup();
   void setup_red_pulse(float radius_overall, unsigned int n_rings, bool broken_line_mode, const glm::vec4 &col=glm::vec4(0.9, 0.2, 0.2, 1.0));
   void setup_green_pulse(bool broken_line_mode);
   void set_scales(const glm::vec2 &s) { scales = s; scales_have_been_set = true; }
   void set_offset_positions(const glm::vec2 &p) { offset_positions = p; offset_positions_have_been_set = true; }
   void update_buffers_for_pulse(float delta_time, int direction=1); // delta time in ms.
   void update_buffers_for_invalid_residue_pulse(unsigned int n_times_called);
   void update_buffers_by_resize(float resize_factor); // say 1.1 for growing rings
   void setup_vertices_and_indices(const std::vector<s_generic_vertex> &vertices,
                                   const std::vector<unsigned int> &indices); // calls setup().
   void update_vertices_and_indices(const std::vector<s_generic_vertex> &vertices,
                                    const std::vector<unsigned int> &indices); // no call to setup(). Just update sub buffer data
   void update_radius_ring_vertices(float new_radius);
   void draw(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation, bool use_view_rotation=false);
   void draw(Shader *shader_p, const glm::vec3 &atom_position, const glm::mat4 &mvp,
             const glm::mat4 &view_rotation, bool use_view_rotation=false);
   void clear();
   bool empty() const { return (vertices.size() == 0); }
};


#endif // LINES_MESH_HH
