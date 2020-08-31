
#include "HUDTextureMesh.hh"

#include <iostream>

void
HUDTextureMesh::init() {

   draw_this_mesh = true;
   vao = 99999999; // unset
}

void
HUDTextureMesh::setup_quad() {

   vertices.clear();
   triangles.clear();

   // invert Y as we do the map
   vertices.push_back(HUDTextureMesh_attribs_t(glm::vec2(-1.0f,  1.0f), glm::vec2(0.0f, 0.0f)));
   vertices.push_back(HUDTextureMesh_attribs_t(glm::vec2( 1.0f,  1.0f), glm::vec2(1.0f, 0.0f)));
   vertices.push_back(HUDTextureMesh_attribs_t(glm::vec2( 1.0f, -1.0f), glm::vec2(1.0f, 1.0f)));
   vertices.push_back(HUDTextureMesh_attribs_t(glm::vec2(-1.0f, -1.0f), glm::vec2(0.0f, 1.0f)));

   triangles.push_back(g_triangle(0, 1, 2));
   triangles.push_back(g_triangle(2, 3, 0));

   set_position_and_scale(glm::vec2(0,0), 1.0);

   setup_buffers();

}

void
HUDTextureMesh::setup_texture_coords_for_nbcs_only() {

   vertices[2].texture_coords.y = 0.5;
   vertices[3].texture_coords.y = 0.5;
   vertices[2].position.y = 0.0;
   vertices[3].position.y = 0.0;
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, 4 * 2 * sizeof(glm::vec2), &(vertices[0]));

}


void
HUDTextureMesh::setup_texture_coords_for_nbcs_and_rama() {

   vertices[2].texture_coords.y = 1.0;
   vertices[3].texture_coords.y = 1.0;
   vertices[2].position.y = -1.0;
   vertices[3].position.y = -1.0;
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, 4 * 2 * sizeof(glm::vec2), &(vertices[0]));

}

void
HUDTextureMesh::set_position_and_scale(const glm::vec2 &pos, float scale_in) {

   position =  pos;
   scale = scale_in;
}


void
HUDTextureMesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   bool first_time = true; // make this member data if needed.

   if (first_time)
      glGenVertexArrays(1, &vao);

   glBindVertexArray(vao);

   unsigned int n_vertices = vertices.size(); // 4

   if (first_time) {
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * 2 * sizeof(glm::vec2), &(vertices[0]), GL_STATIC_DRAW);
   } else {
      glDeleteBuffers(1, &buffer_id);
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * 2 * sizeof(glm::vec2), &(vertices[0]), GL_STATIC_DRAW);
   }

   // position (of the quad)
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(glm::vec2), 0);
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2 *  sizeof(glm::vec2), reinterpret_cast<void *>(sizeof(glm::vec2)));

   // 2 triangles - it's a quad - instanced.
   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);

   if (first_time) {
      glGenBuffers(1, &index_buffer_id);
      GLenum err = glGetError(); if (err) std::cout << "GL error HUDTextureMesh setup_buffers()\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
   } else {
      glDeleteBuffers(1, &index_buffer_id);
      glGenBuffers(1, &index_buffer_id);
      GLenum err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
   }

   // std::cout << "HUDMesh::setup_buffers() indices " << n_bytes << " bytes" << std::endl;
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_DYNAMIC_DRAW);
   GLenum err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_simple_triangles()\n";

   glDisableVertexAttribArray(0);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

   first_time = false;

}

void
HUDTextureMesh::draw(Shader *shader_p) {

   if (! draw_this_mesh) return;

   shader_p->Use();

   if (vao == 99999999)
      std::cout << "error:: You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;

   glBindVertexArray(vao);

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id); // needed?

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);

   shader_p->set_vec2_for_uniform("position", position);
   shader_p->set_float_for_uniform("scale", scale);

   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id); // needed?
   glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
   GLenum err = glGetError();
   if (err) std::cout << "error HUDMesh::draw() glDrawElementsInstanced()"
                      << " of HUDMesh \"" << name << "\""
                      << " with shader" << shader_p->name
                      << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glUseProgram(0);

}

