
#include <iostream>
#include <epoxy/gl.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext.hpp>
#include "LinesMesh.hh"

// should this have its own header?
glm::vec3 coord_orth_to_glm(const clipper::Coord_orth &co);


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
LinesMesh::draw(Shader *shader_p, const glm::mat4 &mvp) {

   shader_p->Use();
   glBindVertexArray(vao);
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   GLenum err = glGetError(); if (err) std::cout << "GL error A LinesMesh::draw()\n";
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "   GL error LinesMesh::draw() glBindBuffer() v "
                                          << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "   GL error LinesMesh::draw() glBindBuffer() i "
                                          << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError(); if (err) std::cout << "   error:: " << shader_p->name
                                          << " LinesMesh.draw() post mvp uniform "
                                          << err << std::endl;
   GLuint n_vertices = indices.size();
   // std::cout << "debug:: LinesMesh draw() drawing n_vertices " << n_vertices << std::endl;
   glDrawElements(GL_LINES, n_vertices, GL_UNSIGNED_INT, nullptr);
   err = glGetError(); if (err) std::cout << "   GL error LinesMesh::draw() glDrawElements()"
                                          << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glUseProgram(0);
}

void
LinesMesh::setup(Shader *shader_p) {

   shader_p->Use();

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
