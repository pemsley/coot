/*
 * src/LinesMesh.cc
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

#include <iostream>
#include <epoxy/gl.h>

#define GLM_ENABLE_EXPERIMENTAL
// #include <glm/ext.hpp>
#include "LinesMesh.hh"

// should this have its own header?
glm::vec3 coord_orth_to_glm(const clipper::Coord_orth &co);


void
LinesMesh::init() {

   index_buffer_id = 999999;
   buffer_id = 999999;
   vao = VAO_NOT_SET;
   scales_have_been_set = false;
   offset_positions_have_been_set = false;
   scales = glm::vec2(1.0, 1.0);
   offset_positions = glm::vec2(0,0);
   first_time = true;
}

#ifdef COMPILED_WITH_CLIPPER
LinesMesh::LinesMesh(const clipper::Cell &cell) {

   float corners[8][3] = {
                          {0,0,0}, //0
                          {0,0,1}, //1
                          {0,1,0}, //2
                          {0,1,1}, //3
                          {1,0,0}, //4
                          {1,0,1}, //5
                          {1,1,0}, //6
                          {1,1,1}};//7

   init();

   vertices.resize(8);
   for (int ii=0; ii<8; ii++) {
      clipper::Coord_frac c_f(corners[ii][0],corners[ii][1],corners[ii][2]);
      clipper::Coord_orth c_o = c_f.coord_orth(cell);
      glm::vec3 pos = coord_orth_to_glm(c_o);
      // std::cout << ii << " " << glm::to_string(pos) << std::endl;
      vertices[ii].pos = pos;
      vertices[ii].normal = glm::vec3(0,0,1); // not used
      vertices[ii].color  = glm::vec4(0.6, 0.6, 0.1, 1.0);
   }

   indices.push_back(0); indices.push_back(1);
   indices.push_back(1); indices.push_back(3);
   indices.push_back(3); indices.push_back(2);
   indices.push_back(2); indices.push_back(0);

   indices.push_back(4); indices.push_back(5);
   indices.push_back(5); indices.push_back(7);
   indices.push_back(7); indices.push_back(6);
   indices.push_back(6); indices.push_back(4);

   indices.push_back(0); indices.push_back(4);
   indices.push_back(1); indices.push_back(5);
   indices.push_back(2); indices.push_back(6);
   indices.push_back(3); indices.push_back(7);

}
#endif

void
LinesMesh::setup_vertices_and_indices(const std::vector<s_generic_vertex> &vertices_in,
                                      const std::vector<unsigned int> &indices_in) {

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: --- update_vertices_and_indices() start" << std::endl;

   vertices = vertices_in;
   indices = indices_in;

   setup();

   if (false) {
      std::cout << "debug::::::::: setup_vertices_and_indices() vertices.size " << vertices.size() << std::endl;
      std::cout << "debug::::::::: setup_vertices_and_indices() indices.size " << indices.size() << std::endl;
      std::cout << "debug::::::::: post setup_vertices_and_indices() vao is " << vao << std::endl;
   }
}

void
LinesMesh::update_vertices_and_indices(const std::vector<s_generic_vertex> &vertices_in,
                                       const std::vector<unsigned int> &indices_in) {

   // 20211004-PE
   // When then object is initialized, this function is called then setup() is called.
   //
   // When it is actually displayed, then this function is called to update the buffer data

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: --- update_vertices_and_indices() start" << std::endl;

   // Here check that the new vertices and indices vectors are smaller or equal to the
   // starting sizes (500, 1500).

   vertices = vertices_in;
   indices = indices_in; // the size of indices is used in the glDrawElements() function in draw()
                         // we don't actually need the indices or the vertices as their
                         // contents gets shoved into the buffer data below and that's the end
                         // of it. This can be reworked if needed. 20211004-PE

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: update_vertices_and_indices() You forgot to setup this LinesMesh "
                << name << " " << std::endl;
   glBindVertexArray(vao);

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: A LinesMesh::update_vertices_and_indices() " << err << "\n";

   unsigned int n_vertices = vertices.size();
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: LinesMesh::update_vertices_and_indices() B1 " << err << "\n";

   // std::cout << "debug:: update buffersubdata with vertices of size " << vertices.size() << std::endl;
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(s_generic_vertex), &(vertices[0]));
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: LinesMesh::update_vertices_and_indices() B2 " << err << "\n";

   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: C LinesMesh::update_vertices_and_indices() " << err << "\n";
   unsigned int n_bytes = indices.size() * sizeof(unsigned int);
   glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, n_bytes, &(indices[0]));

}

void
LinesMesh::update_radius_ring_vertices(float new_radius) {

   float r = new_radius;
   unsigned int n_points = vertices.size();
   for (unsigned int i=0; i<n_points; i++) {
      double theta = 2.0 * M_PI * static_cast<double>(i) / 100.0;
      glm::vec3 pt(r * cos(theta), r * sin(theta), 0.0);
      vertices[i].pos = pt;
   }
   update_vertices_and_indices(vertices, indices);
}



void
LinesMesh::clear() {

   vertices.clear();
   indices.clear();
}

void
LinesMesh::draw(Shader *shader_p,
                const glm::mat4 &mvp, const glm::mat4 &view_rotation, bool use_view_rotation) {

   if (vertices.empty()) return;
   if (indices.empty()) return;

   GLenum err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw() -- start --\n";
   shader_p->Use();
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw A()\n";
   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: LinesMesh::draw() You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;
   glBindVertexArray(vao);
   err = glGetError(); if (err) std::cout << "GL ERROR:: LinesMesh::draw() B binding vao " << vao << "\n";
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   err = glGetError(); if (err) std::cout << "GL ERROR:: LinesMesh::draw C()\n";

   // no atom_position uniform

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError(); if (err) std::cout << "error:: " << shader_p->name
                                          << " LinesMesh.draw() post mvp uniform "
                                          << err << std::endl;

   if (scales_have_been_set)
      shader_p->set_vec2_for_uniform("scales", scales);
   if (offset_positions_have_been_set)
      shader_p->set_vec2_for_uniform("offset_positions", offset_positions);

   GLuint n_indices = indices.size();

   if (false)
      std::cout << "DEBUG:: LinesMesh draw() drawing n_indices " << n_indices
                << " vertices.size() " << vertices.size() << std::endl;

   glDrawElements(GL_LINES, n_indices, GL_UNSIGNED_INT, nullptr);
   err = glGetError(); if (err) std::cout << "error LinesMesh::draw() glDrawElements()"
                                          << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glUseProgram(0);
}

// identification pulse
void
LinesMesh::draw(Shader *shader_p, const glm::vec3 &atom_position,
                const glm::mat4 &mvp, const glm::mat4 &view_rotation, bool use_view_rotation) {

   if (vertices.empty()) return;
   if (indices.empty()) return;

   GLenum err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw() -- start --\n";
   shader_p->Use();
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw A()\n";
   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: LinesMesh::draw() (identification pulse) You forgot to setup this mesh "
                << name << " " << shader_p->name << std::endl;
   glBindVertexArray(vao);
   err = glGetError(); if (err) std::cout << "GL ERROR:: LinesMesh::draw() B vao\n";
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   err = glGetError(); if (err) std::cout << "GL ERROR:: LinesMesh::draw C()\n";

   // we are using 2 (at the moment) different shaders for this class.
   // Hmmm... lines-pulse.shader uses atom_position, and that's the only shader
   // for this mesh. so use_view_rotation should always be true?
   if (use_view_rotation) {
      // std::cout << "sending atom_centre " << glm::to_string(atom_position) << std::endl;
      glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
      shader_p->set_vec3_for_uniform("atom_centre", atom_position);
   }

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError(); if (err) std::cout << "GL ERROR:: " << shader_p->name
                                          << " LinesMesh.draw() post mvp uniform "
                                          << err << std::endl;

   GLuint n_vertices = indices.size();

   if (false)
      std::cout << "debug:: LinesMesh draw() drawing n_vertices " << n_vertices << std::endl;

   glDrawElements(GL_LINES, n_vertices, GL_UNSIGNED_INT, nullptr);
   err = glGetError(); if (err) std::cout << "error LinesMesh::draw() glDrawElements()"
                                          << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glUseProgram(0);
}

void
LinesMesh::setup() {

   if (vertices.empty()) std::cout << "error:: LinesMesh::setup() called before vertices filled " << std::endl;
   if (indices.empty())  std::cout << "error:: LinesMesh::setup() called before indices filled " << std::endl;

   if (vertices.empty()) return;
   if (indices.empty()) return;

   if (vao == VAO_NOT_SET)
      glGenVertexArrays (1, &vao);

   glBindVertexArray (vao);

   if (! first_time)
      glDeleteBuffers(GL_ARRAY_BUFFER, &buffer_id);

   glGenBuffers(1, &buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   unsigned int n_vertices = vertices.size();
   glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_STATIC_DRAW);

   // position
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex), 0);
   // normal
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                         reinterpret_cast<void *>(sizeof(glm::vec3)));
   // colour
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                         reinterpret_cast<void *>(2 * sizeof(glm::vec3)));

   if (! first_time)
      glDeleteBuffers(1, &index_buffer_id);

   glGenBuffers(1, &index_buffer_id);
   GLenum err = glGetError(); if (err) std::cout << "GL error A LinesMesh::setup()\n";
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error B LinesMesh::setup()\n";
   unsigned int n_bytes = indices.size() * sizeof(unsigned int);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &indices[0], GL_STATIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error B LinesMesh::setup() glBufferData()\n";

}

void
LinesMesh::make_vertices_for_pulse(const glm::vec4 &colour, float radius_overall,
                                   unsigned int n_rings, float theta_offset, bool broken_mode) {

   // we have just simply constructed this object. Now add some vertices.

   // we don't need to regenerate the indices if they have already been generated

#if 0 // quick fix for glm compilation problems
   std::cout << "debug::  make_vertices_for_pulse() --- start --- "
             << glm::to_string(colour) << " r: "
             << radius_overall << " n-rings: " << n_rings << " "
             << std::endl;
#endif

   vertices.clear();
   indices.clear();

   glm::vec3 normal(0,0,1);
   const unsigned int n_line_segments = 30;
   for (unsigned int j=0; j<n_rings; j++) {
      unsigned int idx_base = vertices.size();
      float r = 0.06f * radius_overall * static_cast<float>(j+1);
      for (unsigned int i=0; i<n_line_segments; i++) {
         float theta = 2.0 * M_PI * static_cast<float>(i) / static_cast<float>(n_line_segments);
         theta += theta_offset;
         glm::vec3 p(r * sinf(theta), r * cosf(theta), 0.0);
         s_generic_vertex v(p, normal, colour);
         // std::cout << "vertex position " << glm::to_string(v.pos) << " " << glm::to_string(colour) << std::endl;
         vertices.push_back(v);
      }
      for (unsigned int i=0; i<n_line_segments; i++) {
         if (broken_mode)
            if ((j+i) %2 == 0)
               continue;
         unsigned int i_this = i;
         unsigned int i_next = i+1;
         if (i_next == n_line_segments)
            i_next = 0;
         indices.push_back(i_this + idx_base);
         indices.push_back(i_next + idx_base);
      }
   }
}

void
LinesMesh::setup_green_pulse(bool broken_line_mode) {

   glm::vec4 colour(0.2, 0.8, 0.2, 1.0); // green(!?)
   unsigned int n_rings = 13;

#if 0 // quick fix for glm compilation problems
   std::cout << "DEBUG:: calling make_vertices_for_pulse() with colour " << glm::to_string(colour) << std::endl;
#endif
   make_vertices_for_pulse(colour, 13.0, n_rings, 0.0, broken_line_mode);
#if 0 // quick fix for glm compilation problems
   std::cout << "DEBUG:: done make_vertices_for_pulse() " << vertices.size() << " " << indices.size() << std::endl;
#endif
   setup();
}

void
LinesMesh::setup_red_pulse(float radius_overall, unsigned int n_rings, bool broken_line_mode, const glm::vec4 &col) {

   // unsigned int n_rings = 3;
   // float radius_overall = 6.0;
   make_vertices_for_pulse(col, radius_overall, n_rings, 0.0, broken_line_mode);
   setup();
}


void
LinesMesh::update_buffers_for_pulse(float n_steps, int direction) { // delta time in ms.

   float r = 0.4 * n_steps;
   glm::vec4 colour(0.2, 0.8, 0.2, 1.0);
   if (direction == -1) {
      colour = glm::vec4(0.8, 0.4, 0.5, 1.0);
      r = 6.0 - 6.0/20.0 * static_cast<float>(n_steps);
   }

   // nice swirl effect

   float theta_offset = 0.0;
   if (direction == -1)
      theta_offset = static_cast<float>(n_steps) * - 0.05;

   bool broken_line_mode =  true;
   unsigned int n_vertices = vertices.size();
   unsigned int n_rings = 3;
   glBindVertexArray(vao);
   make_vertices_for_pulse(colour, r, n_rings, theta_offset, broken_line_mode);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(s_generic_vertex), &(vertices[0]));
}

void
LinesMesh::update_buffers_by_resize(float resize_factor) { // say 1.1 for growing rings

   unsigned int n_vertices = vertices.size();
   glBindVertexArray(vao);
   for (auto &v : vertices)
      v.pos *= resize_factor;
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(s_generic_vertex), &(vertices[0]));
}



void
LinesMesh::update_buffers_for_invalid_residue_pulse(unsigned int n_times_called) {

   unsigned int n_vertices = vertices.size();
   float r_1 = static_cast<float>(n_times_called)/48.0;
   float r = 2.5 + 2.5 * sin(r_1 * 9.0);
   // std::cout << "update_buffers_for_invalid_residue_pulse: n " << n_times_called << " r " << r << std::endl;
   float theta_offset = 0.0;
   bool broken_line_mode = false;
   glBindVertexArray(vao);
   unsigned int n_rings = 1;
   glm::vec4 colour(0.8, 0.2, 0.2, 1.0);
   make_vertices_for_pulse(colour, r, n_rings, theta_offset, broken_line_mode);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(s_generic_vertex), &(vertices[0]));

}
