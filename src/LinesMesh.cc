
#include <iostream>
#include <epoxy/gl.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
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

LinesMesh::LinesMesh(const clipper::Cell &cell) {

   index_buffer_id = 999999;

   float corners[8][3] = {
                          {0,0,0}, //0
                          {0,0,1}, //1
                          {0,1,0}, //2
                          {0,1,1}, //3
                          {1,0,0}, //4
                          {1,0,1}, //5
                          {1,1,0}, //6
                          {1,1,1}};//7

   vertices.resize(8);
   for (int ii=0; ii<8; ii++) {
      clipper::Coord_frac c_f(corners[ii][0],corners[ii][1],corners[ii][2]);
      clipper::Coord_orth c_o = c_f.coord_orth(cell);
      glm::vec3 pos = coord_orth_to_glm(c_o);
      std::cout << ii << " " << glm::to_string(pos) << std::endl;
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

void
LinesMesh::update_vertices_and_indices(const std::vector<s_generic_vertex> &vertices_in,
                                       const std::vector<unsigned int> &indices_in) {

   vertices = vertices_in;
   indices = indices_in;
   setup();
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
      std::cout << "You forgot to setup this mesh " << name << " " << shader_p->name << std::endl;
   glBindVertexArray(vao);
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw() B binding vao\n";
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw C()\n";

   // no atom_position uniform

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError(); if (err) std::cout << "error:: " << shader_p->name
                                          << " LinesMesh.draw() post mvp uniform "
                                          << err << std::endl;

   if (scales_have_been_set)
      shader_p->set_vec2_for_uniform("scales", scales);
   if (offset_positions_have_been_set)
      shader_p->set_vec2_for_uniform("offset_positions", offset_positions);

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
      std::cout << "You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;
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

   if (first_time)
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

   vertices.clear();
   indices.clear();

   glm::vec3 normal(0,0,1);
   const unsigned int n_line_segments = 30;
   for (unsigned int j=0; j<n_rings; j++) {
      unsigned int idx_base = vertices.size();
      float r = 0.06 * radius_overall * static_cast<float>(j+1);
      for (unsigned int i=0; i<n_line_segments; i++) {
         float theta = 2.0 * M_PI * static_cast<float>(i) / static_cast<float>(n_line_segments);
         theta += theta_offset;
         glm::vec3 p(r * sinf(theta), r * cosf(theta), 0.0);
         s_generic_vertex v(p, normal, colour);
         // std::cout << "vertex position " << glm::to_string(v.pos) << std::endl;
         vertices.push_back(v);
      }
      for (unsigned int i=0; i<n_line_segments; i++) {
         if (broken_mode)
            if ((j+i) %2 == 0) continue;
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
LinesMesh::setup_pulse(bool broken_line_mode) {

   glm::vec4 colour(0.2, 0.8, 0.2, 1.0);
   unsigned int n_rings = 3;
   make_vertices_for_pulse(colour, 2.0, n_rings, 0.0, broken_line_mode);
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
