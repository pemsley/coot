
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
   vao = 999999;
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
LinesMesh::clear() {

   vertices.clear();
   indices.clear();
}

void
LinesMesh::draw(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation, bool use_view_rotation) {

   if (vertices.empty()) return;
   if (indices.empty()) return;

   GLenum err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw() -- start --\n";
   shader_p->Use();
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw A()\n";
   if (vao == 999999)
      std::cout << "You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;
   glBindVertexArray(vao);
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw() B vao\n";
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   err = glGetError(); if (err) std::cout << "error:: LinesMesh::draw C()\n";

   if (use_view_rotation) {
      glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
      shader_p->set_vec3_for_uniform("atom_centre", central_position);
   }

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError(); if (err) std::cout << "error:: " << shader_p->name
                                          << " LinesMesh.draw() post mvp uniform "
                                          << err << std::endl;

   GLuint n_vertices = indices.size();
   // std::cout << "debug:: LinesMesh draw() drawing n_vertices " << n_vertices << std::endl;
   glDrawElements(GL_LINES, n_vertices, GL_UNSIGNED_INT, nullptr);
   err = glGetError(); if (err) std::cout << "error LinesMesh::draw() glDrawElements()"
                                          << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glUseProgram(0);
}

void
LinesMesh::setup(Shader *shader_p) {

   shader_p->Use();

   if (vertices.empty()) std::cout << "error:: LinesMesh::setup() called before vertices filled " << std::endl;
   if (indices.empty())  std::cout << "error:: LinesMesh::setup() called before indices filled " << std::endl;

   if (vertices.empty()) return;
   if (indices.empty()) return;

   glGenVertexArrays (1, &vao);
   glBindVertexArray (vao);

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

   glGenBuffers(1, &index_buffer_id);
   GLenum err = glGetError(); if (err) std::cout << "GL error A LinesMesh::setup()\n";
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error B LinesMesh::setup()\n";
   unsigned int n_bytes = indices.size() * sizeof(unsigned int);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &indices[0], GL_STATIC_DRAW);

}

void
LinesMesh::make_vertices_for_pulse(float radius_overall) {

   // we have just simply constructed this object. Now add some vertices.

   // we don't need to regenerate the indices if they have already been generated
   
   vertices.clear();
   indices.clear();

   glm::vec3 normal(0,0,1);
   glm::vec4 colour(0.2, 0.8, 0.2, 1.0);
   const unsigned int n_line_segments = 30;
   unsigned int n_rings = 3;
   for (unsigned int j=0; j<n_rings; j++) {
      unsigned int idx_base = vertices.size();
      float r = 0.06 * radius_overall * static_cast<float>(j+1);
      for (unsigned int i=0; i<n_line_segments; i++) {
         float theta = 2.0 * M_PI * static_cast<float>(i) / static_cast<float>(n_line_segments);
         glm::vec3 p(r * sinf(theta), r * cosf(theta), 0.0);
         s_generic_vertex v(p + central_position, normal, colour);
         // std::cout << "vertex position " << glm::to_string(v.pos) << std::endl;
         vertices.push_back(v);
      }
      for (unsigned int i=0; i<n_line_segments; i++) {
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
LinesMesh::setup_pulse(const glm::vec3 &central_position_in, Shader *shader_p) {

   central_position = central_position_in;
   make_vertices_for_pulse(2.01);
   setup(shader_p);
}



void
LinesMesh::update_buffers_for_pulse(float n_steps) { // delta time in ms.

   float r = 0.4 * n_steps;
   unsigned int n_vertices = vertices.size();
   glBindVertexArray(vao);
   make_vertices_for_pulse(r);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(s_generic_vertex), &(vertices[0]));
}
