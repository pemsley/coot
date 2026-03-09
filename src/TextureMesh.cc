/*
 * src/TextureMesh.cc
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

#include "stereo-eye.hh"
#ifdef USE_PYTHON
#include "Python.h"
#endif

#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
// #include <glm/ext.hpp> // 20240326-PE
#include <glm/gtx/string_cast.hpp> // for to_string()
#include <glm/gtc/type_ptr.hpp>  // for value_ptr() 20240326-PE


#include "ft-character.hh"
#include "TextureMesh.hh"

// #define THIS_IS_HMT

#ifdef THIS_IS_HMT
#include "display-info.hh"
#else
#include "graphics-info.h"
#endif

#ifdef THIS_IS_HMT
#include "coot-utils.hh" // added to compile
#endif

#include <Texture.hh>

// Define these only in *one* .cc file.
// #define TINYGLTF_IMPLEMENTATION
// #define STB_IMAGE_IMPLEMENTATION
// #define STB_IMAGE_WRITE_IMPLEMENTATION

#include "coot-utils/tiny_gltf.h"


// for the moment make the scales explict, when fixed make the scales default
void
TextureMesh::setup_camera_facing_quad(float scale_x, float scale_y, float offset_x, float offset_y) {

   draw_this_mesh = true;

   glm::vec3 n(0,0,1);
   glm::vec4 col(1.0, 1.0, 1.0, 1.0);

   vertices.clear();
   triangles.clear();

   // the indexing might well be wrong here - I'm sort of guessing
   vertices.push_back(TextureMeshVertex(glm::vec3(-scale_x,  scale_y, 0.0f), n, col, glm::vec2(0,0)));
   vertices.push_back(TextureMeshVertex(glm::vec3( scale_x,  scale_y, 0.0f), n, col, glm::vec2(1,0)));
   vertices.push_back(TextureMeshVertex(glm::vec3( scale_x, -scale_y, 0.0f), n, col, glm::vec2(1,1)));
   vertices.push_back(TextureMeshVertex(glm::vec3(-scale_x, -scale_y, 0.0f), n, col, glm::vec2(0,1)));

   // angry_diego has y = 0 at the bottom of the image
   //
   for (unsigned int i=0; i<vertices.size(); i++)
      vertices[i].position += glm::vec3(offset_x, offset_y, 0.0f);

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));

   setup_buffers();

}

void
TextureMesh::setup_tomo_quad(float scale_x, float scale_y, float x_offset, float y_offset, float z_pos,
                             bool texture_x_y_swap_flag) {

  draw_this_mesh = true;

   glm::vec3 n(0,0,1);
   glm::vec4 col(1.0, 1.0, 1.0, 1.0);

   vertices.clear();
   triangles.clear();

   // the logic/indexing might well be wrong here - fix in future

   if (texture_x_y_swap_flag) {
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset,           y_offset + scale_y, z_pos), n, col, glm::vec2(0,1)));
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset + scale_x, y_offset + scale_y, z_pos), n, col, glm::vec2(1,1)));
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset + scale_x, y_offset, z_pos), n, col, glm::vec2(1,0)));
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset,           y_offset, z_pos), n, col, glm::vec2(0,0)));
   } else {
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset,           y_offset,           z_pos), n, col, glm::vec2(0,0)));
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset + scale_x, y_offset,           z_pos), n, col, glm::vec2(0,1)));
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset + scale_x, y_offset + scale_y, z_pos), n, col, glm::vec2(1,1)));
      vertices.push_back(TextureMeshVertex(glm::vec3(x_offset,           y_offset + scale_y, z_pos), n, col, glm::vec2(1,0)));
   }

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));

   setup_buffers();

}


void
TextureMesh::set_colour(const glm::vec4 &col_in) {

   for (unsigned int i=0; i<vertices.size(); i++) {
      vertices[i].color = col_in;
   }
}

// static
std::string
TextureMesh::_(int err) {

   std::string s = std::to_string(err);
   if (err == GL_INVALID_ENUM)      s = "GL_INVALID_ENUM";
   if (err == GL_INVALID_OPERATION) s = "GL_INVALID_OPERATION";
   if (err == GL_INVALID_VALUE)     s = "GL_INVALID_VALUE";
   return s;
}

void
TextureMesh::setup_buffers() {

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers() --- start --- " << _(err) << "\n";
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers() --- start --- " << _(err) << "\n";
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers() --- start --- " << _(err) << "\n";

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   glGenVertexArrays(1, &vao);
   glBindVertexArray(vao);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers() A" << _(err) << std::endl;

   setup_tbn(vertices.size());

   glGenBuffers(1, &buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   unsigned int n_vertices = vertices.size();
   if (false)
      std::cout << "DEBUG:: in TextureMesh::setup_buffers() " << name << " n_vertices is " << n_vertices
                << " buffer_id " << buffer_id << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(TextureMeshVertex), &(vertices[0]), GL_STATIC_DRAW);
   // std::cout << "in TextureMesh::setup_buffers() " << name << " done glBufferData() " << std::endl;

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers() B\n";

   // position
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex), 0);

   // normal
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                          reinterpret_cast<void *>(sizeof(glm::vec3)));

   // tangent
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                          reinterpret_cast<void *>(2 * sizeof(glm::vec3)));

   // bitangent
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                          reinterpret_cast<void *>(3 * sizeof(glm::vec3)));

   // colour
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                         reinterpret_cast<void *>(4 * sizeof(glm::vec3)));

   // texture coordinates
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 2, GL_FLOAT, GL_FALSE, sizeof(TextureMeshVertex),
                         reinterpret_cast<void *>(4 * sizeof(glm::vec3) + sizeof(glm::vec4)));


   glGenBuffers(1, &index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers()" << _(err) << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL ERROR:: TextureMesh::setup_buffers()" << _(err) << std::endl;
   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);
   if (false)
      std::cout << "debug:: in TextureMesh::setup_buffers(): " << name
                << " glBufferData for index buffer_id " << index_buffer_id
                << " n_triangles: " << n_triangles
                << " allocating with size: " << n_bytes << " bytes" << std::endl;
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_STATIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL ERROR TextureMesh::setup_buffers()" << _(err) << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

}

void
TextureMesh::setup_tbn(unsigned int n_vertices) {

   auto get_tangent_bitangent = [] (const TextureMeshVertex &vertex_1,
                                    const TextureMeshVertex &vertex_2,
                                    const TextureMeshVertex &vertex_3) {

                                   const glm::vec3 &pos_1 = vertex_1.position;
                                   const glm::vec3 &pos_2 = vertex_2.position;
                                   const glm::vec3 &pos_3 = vertex_3.position;

                                   const glm::vec2 &uv_1 = vertex_1.texCoord;
                                   const glm::vec2 &uv_2 = vertex_2.texCoord;
                                   const glm::vec2 &uv_3 = vertex_3.texCoord;

                                   glm::vec3 tangent, bitangent;
                                   glm::vec3 edge_1 = pos_2 - pos_1;
                                   glm::vec3 edge_2 = pos_3 - pos_1;
                                   glm::vec2 deltaUV_1 = uv_2 - uv_1;
                                   glm::vec2 deltaUV_2 = uv_3 - uv_1;
                                   GLfloat f = 1.0f/(deltaUV_1.x * deltaUV_2.y - deltaUV_2.x * deltaUV_1.y);

                                   tangent.x = f * (deltaUV_2.y * edge_1.x - deltaUV_1.y * edge_2.x);
                                   tangent.y = f * (deltaUV_2.y * edge_1.y - deltaUV_1.y * edge_2.y);
                                   tangent.z = f * (deltaUV_2.y * edge_1.z - deltaUV_1.y * edge_2.z);
                                   tangent =  glm::normalize(tangent);

                                   bitangent.x = f * (-deltaUV_2.x * edge_1.x - deltaUV_1.x * edge_2.x);
                                   bitangent.y = f * (-deltaUV_2.x * edge_1.y - deltaUV_1.x * edge_2.y);
                                   bitangent.z = f * (-deltaUV_2.x * edge_1.z - deltaUV_1.x * edge_2.z);
                                   bitangent =  glm::normalize(bitangent);

                                   return std::make_pair(tangent, bitangent);
                                };

   for (const auto &tri : triangles) {
      unsigned int idx_0 = tri[0];
      unsigned int idx_1 = tri[1];
      unsigned int idx_2 = tri[2];
      // std::cout << "indices " << idx_0 << " " << idx_1 << " " << idx_2 << " vs "<< vertices.size() << std::endl;
      if (idx_0 < n_vertices) {
         if (idx_1 < n_vertices) {
            if (idx_2 < n_vertices) {
               TextureMeshVertex &v0 = vertices[idx_0];
               TextureMeshVertex &v1 = vertices[idx_1];
               TextureMeshVertex &v2 = vertices[idx_2];
               std::pair<glm::vec3, glm::vec3> tangent_bitangent = get_tangent_bitangent(v0, v1, v2);

               v0.tangent   = tangent_bitangent.first;
               v0.bitangent = tangent_bitangent.second;
               v1.tangent   = tangent_bitangent.first;
               v1.bitangent = tangent_bitangent.second;
               v2.tangent   = tangent_bitangent.first;
               v2.bitangent = tangent_bitangent.second;
            }
         }
      }
   }
}


void
TextureMesh::draw_atom_label(const std::string &atom_label,
                             const glm::vec3 &atom_label_position,
                             const glm::vec4 &text_colour, // set using glBufferSubData
                             Shader *shader_p,
                             stereo_eye_t eye,
                             const glm::mat4 &mvp,
                             const glm::mat4 &view_rotation_matrix,
                             const glm::vec4 &background_colour,
                             bool do_depth_fog,
                             bool is_perspective_projection) {

   auto get_stereo_x_scale_and_offset = [] (stereo_eye_t eye) {

      // for side by side stereo, of course
      float stereo_x_scale  = 1.0;
      float stereo_x_offset = 0.0;
      if (eye == stereo_eye_t::LEFT_EYE) {
         stereo_x_scale = 2.0f;
      }
      if (eye == stereo_eye_t::RIGHT_EYE) {
         stereo_x_offset = -0.5f;
      }
      return std::pair<float, float>(stereo_x_scale, stereo_x_offset);
   };


   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   // std::cout << "debug:: in draw_atom_label() n-trianges: " << n_triangles << " n-verts: " << n_verts  << std::endl;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() " << shader_p->name << " -- start -- error "
                      << err << std::endl;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   shader_p->set_vec3_for_uniform("label_position", atom_label_position);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post set label_position " << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post set mvp " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post set view rotation " << err << std::endl;

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() :" << name << ": " << shader_p->name
                      << " post background_colour " << err << std::endl;
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() " << name << " " << shader_p->name
                      << " post do_depth_fog " << err << std::endl;
   shader_p->set_bool_for_uniform("is_perspective_projection", is_perspective_projection);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() " << name << " " << shader_p->name
                      << " post is_perspective_projection " << err << std::endl;

   std::pair<float, float> stereo_x_scale_and_offset = get_stereo_x_scale_and_offset(eye);
   const float &stereo_x_scale  = stereo_x_scale_and_offset.first;
   const float &stereo_x_offset = stereo_x_scale_and_offset.second;

   // std::cout << "DEBUG:: draw_atom_label() stereo_x_scale and offset " << stereo_x_scale << " " << stereo_x_offset << std::endl;

   shader_p->set_float_for_uniform("stereo_x_scale",  stereo_x_scale);
   shader_p->set_float_for_uniform("stereo_x_offset", stereo_x_offset);


   if (vao == VAO_NOT_SET)
      std::cout << "You forget to setup this TextureMesh " << name << " " << shader_p->name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label()) " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw_atom_label() A3 " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id); // not needed? test by removing
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() glBindBuffer() v " << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id); // needed? test by removing
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() glBindBuffer() i " << err << std::endl;

   glEnableVertexAttribArray(0); // position
   glEnableVertexAttribArray(1); // normal    // not used for atom labels
   glEnableVertexAttribArray(2); // tangent   // ditto
   glEnableVertexAttribArray(3); // bitangent // ditto
   glEnableVertexAttribArray(4); // colour
   glEnableVertexAttribArray(5); // texCoord

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_atom_label() " << name << " pre-draw " << err << std::endl;

   // scale of 20 make letter separation good
   GLfloat scale = 100.1;

   // ------------------------------- text texture code here -----------------------

   // consider refactoring when working

   // display_info_t di;

   float x = 0;
   float y = 0;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   err = glGetError(); if (err) std::cout << "GL ERROR:: draw_atom_label() A0 " << err << std::endl;

#ifdef THIS_IS_HMT
   std::map<GLchar, FT_character> &ft_characters = display_info_t::ft_characters;
#else
   std::map<GLchar, FT_character> &ft_characters = graphics_info_t::ft_characters;
#endif

   for (std::string::const_iterator it_c = atom_label.begin(); it_c != atom_label.end(); ++it_c) {
      // bitsy spider...
      std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*it_c);
      if (it == ft_characters.end()) {
         std::cout << "ERROR:: TextureMesh::draw_atom_label():: Failed to lookup glyph for " << *it_c
                   << " in " << atom_label << std::endl;
         continue;
      };
      const FT_character &ch = it->second;
      GLfloat xpos = x + ch.Bearing.x * scale;
      GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;
      GLfloat w = ch.Size.x * scale;
      GLfloat h = ch.Size.y * scale;

      // std::cout << "---------- char " << *c << " x " << x << std::endl;

      // Update the vertices for each character
      //
      std::vector<TextureMeshVertex> texture_mesh_vertices = vertices;

      // 0,0 -> add h
      // 1,0 -> add w,h
      // 1,1 -> add w

      bool debug = false;

      if (debug) {
         std::cout << "texture_mesh_vertices 0 " << glm::to_string(texture_mesh_vertices[0].position) << std::endl;
         std::cout << "texture_mesh_vertices 1 " << glm::to_string(texture_mesh_vertices[1].position) << std::endl;
         std::cout << "texture_mesh_vertices 2 " << glm::to_string(texture_mesh_vertices[2].position) << std::endl;
         std::cout << "texture_mesh_vertices 3 " << glm::to_string(texture_mesh_vertices[3].position) << std::endl;
         std::cout << "here with w " << w << " and h " << h << std::endl;
      }

      for (unsigned int i=0; i<4; i++)
         texture_mesh_vertices[i].color = text_colour;

      for (unsigned int i=0; i<4; i++)
         texture_mesh_vertices[i].position += glm::vec3(xpos, ypos, 0.0f);
      
      texture_mesh_vertices[0].position.y += h;
      texture_mesh_vertices[1].position.x += w;
      texture_mesh_vertices[1].position.y += h;
      texture_mesh_vertices[2].position.x += w;

      if (debug) {
         std::cout << "post texture_mesh_vertices 0 " << glm::to_string(texture_mesh_vertices[0].position) << std::endl;
         std::cout << "post texture_mesh_vertices 1 " << glm::to_string(texture_mesh_vertices[1].position) << std::endl;
         std::cout << "post texture_mesh_vertices 2 " << glm::to_string(texture_mesh_vertices[2].position) << std::endl;
         std::cout << "post texture_mesh_vertices 3 " << glm::to_string(texture_mesh_vertices[3].position) << std::endl;
      }

      glBindTexture(GL_TEXTURE_2D, ch.TextureID);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      unsigned int n_vertices = vertices.size();
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(TextureMeshVertex), &texture_mesh_vertices[0]);

      unsigned int n_draw_verts = 6;
      glDrawElements(GL_TRIANGLES, n_draw_verts, GL_UNSIGNED_INT, nullptr);

      err = glGetError();
      if (err) std::cout << "TextureMesh::draw_atom_label() glDrawArrays() " << err << std::endl;

       // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
      x += (ch.Advance >> 6) * scale * 1.0;

   }


   // ------------------------------- done text texture code  ----------------------

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_atom_label() glDrawElements()"
                      << " of Mesh \"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);

   glUseProgram (0);
}

void
TextureMesh::draw(Shader *shader_p,
                  const glm::mat4 &mvp,
                  const glm::mat4 &view_rotation_matrix,
                  const std::map<unsigned int, lights_info_t> &lights,
                  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                  const glm::vec4 &background_colour,
                  bool do_depth_fog) {

   // std::cout << "Here in TextureMesh::draw() " << name << " with shader " << shader_p->name  << std::endl;

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw() " << shader_p->name << " -- start -- "
                      << err << std::endl;

   shader_p->Use();

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw::draw() "
                      << shader_p->name << " shader::Use() " << err << std::endl;

   const std::string &shader_name = shader_p->name;

   // Bind the textures
   GLuint unit = 0;
   for (auto &texture : textures) {
      if (false)
         std::cout << "Mesh::draw() " << name << " binding texture " << texture.name << " to unit " << unit << std::endl;
      texture.texture.Bind(unit);
      unit++;
   }

   // std::cout << "mvp: " << glm::to_string(mvp) << std::endl;

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw::draw() "
                      << shader_p->name << " pre mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, glm::value_ptr(mvp));
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw::draw() "
                      << shader_p->name << " post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, glm::value_ptr(view_rotation_matrix));
   err = glGetError();
   if (err) std::cout << "GL ERROR: TextureMesh::draw(): "
                      << shader_p->name << " post view rotation uniform " << err << std::endl;

   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw() " << shader_name << " pre-set eye position"
                      << " with GL err " << err << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw() " << shader_name << " post-set eye position"
                      << " with GL err " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;


   // ------------------------ more texture mesh uniforms here ------------------

   shader_p->set_bool_for_uniform("do_animation", do_animation);
   shader_p->set_float_for_uniform("animation_A", animation_A); // overall
   shader_p->set_float_for_uniform("animation_k", animation_k); // wave number
   shader_p->set_float_for_uniform("animation_w", animation_w); // angular freq (per second)

   std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
   float time = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_constructed).count();

   // std::cout << "sending time " << time << std::endl;
   shader_p->set_float_for_uniform("time", time); // thousands of milliseconds


   // ------------------------ texture mesh uniforms done ------------------

   // this lights block can be in it's own function (same as Mesh)
   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   if (vao == VAO_NOT_SET)
      std::cout << "You forgot to setup this mesh (or setup with empty vertices or triangles) "
                << name << " " << shader_p->name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "   error draw() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   shader_p->set_int_for_uniform("base_texture", 0);

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw() A3 " << err << std::endl;
   glActiveTexture(GL_TEXTURE1);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw() A4 " << err << std::endl;

   glEnableVertexAttribArray(0);  // position
   glEnableVertexAttribArray(1);  // normal
   glEnableVertexAttribArray(2);  // colour (not used)
   glEnableVertexAttribArray(3);  // texture coordinates
   glEnableVertexAttribArray(4);
   glEnableVertexAttribArray(5);

   err = glGetError();
   if (err) std::cout << "   error draw() " << name << " pre-draw " << err << std::endl;

   // if (use_blending) {
   // glEnable(GL_BLEND);
   // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   //}

   // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
   // before making a new VAO?

   if (false)
      std::cout << "debug:: TextureMesh::draw() " << name << " shader " << shader_p->name
                << " vao " << vao
                << " drawing " << n_verts << " triangle vertices"  << std::endl;

   glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "   error TextureMesh::draw() glDrawElements()"
                      << " of Mesh \"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   // if (use_blending) {
   // glDisable(GL_BLEND);
   // }

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);

   glUseProgram (0);

}

// Holy Zarquon swimming fish
void
TextureMesh::set_animation_paramaters(float amplitide_overall, float wave_number, float freq) {

   animation_A = amplitide_overall;
   animation_k = wave_number;
   animation_w = freq;

}



void
TextureMesh::draw_with_shadows(Shader *shader_p,
                               const glm::mat4 &mvp,
                               const glm::mat4 &view_rotation_matrix,
                               const std::map<unsigned int, lights_info_t> &lights,
                               const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                               const glm::vec4 &background_colour,
                               bool do_depth_fog,
                               const glm::mat4 &light_space_mvp,
                               unsigned int shadow_depthMap,
                               float shadow_strength,
                               unsigned int shadow_softness,
                               bool show_just_shadows) {

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_with_shadows() " << shader_p->name << " -- start -- "
                      << err << std::endl;

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   // Bind the textures
   bool specular_map_included = false;
   bool normal_map_included   = false;
   bool reversed_normals      = false;
   GLuint unit = 0;
   GLuint specular_map_unit = 99999;
   GLuint   normal_map_unit = 99999;
   for (auto &texture : textures) {
      texture.texture.Bind(unit);
      if (texture.texture.type == Texture::SPECULAR) {
         specular_map_included = true;
         specular_map_unit = unit;
      }
      if (texture.texture.type == Texture::NORMAL) {
         normal_map_included = true;
         normal_map_unit = unit;
         if (texture.texture.reversed_normals)
            reversed_normals = true;
      }

      if (false)
         std::cout << "Mesh::draw_with_shadows() " << name << " binding texture " << texture.name
                   << " type " << texture.texture.type
                   << " to unit " << unit << " texture-handle: " << texture.texture.m_texture_handle << std::endl;

      unit++;
   }

   // std::cout << "mvp: " << glm::to_string(mvp) << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, glm::value_ptr(mvp));
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw::draw_with_shadows() "
                      << shader_p->name << " post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, glm::value_ptr(view_rotation_matrix));
   err = glGetError();
   if (err) std::cout << "GL ERROR: TextureMesh::draw_with_shadows(): "
                      << shader_p->name << " post view rotation uniform " << err << std::endl;

   shader_p->set_mat4_for_uniform("light_space_mvp", light_space_mvp);
   err = glGetError();
   if (err) std::cout << "GL ERROR: TextureMesh::draw_with_shadows(): "
                      << shader_p->name << " post light-space-mvp " << err << std::endl;

   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   shader_p->set_bool_for_uniform("reversed_normals", reversed_normals);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_with_shadows() " << shader_name << " pre-set eye position"
                      << " with GL err " << err << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_with_shadows() " << shader_name << " post-set eye position"
                      << " with GL err " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_with_shadows() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   // this lights block can be in its own function (same as Mesh)
   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: TextureMess::draw_with_shadows() You forgot to setup this mesh"
                << " (or setup with empty vertices or triangles) " << name << " " << shader_p->name
                << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_with_shadows() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   shader_p->set_bool_for_uniform("specular_map_included", specular_map_included);
   shader_p->set_bool_for_uniform("normal_map_included",     normal_map_included);
   shader_p->set_int_for_uniform("base_texture", 0);
   if (specular_map_included)
      shader_p->set_int_for_uniform("specular_map", specular_map_unit);
   if (normal_map_included)
      shader_p->set_int_for_uniform("normal_map",     normal_map_unit);
   shader_p->set_int_for_uniform("shadow_map",   4);
   shader_p->set_float_for_uniform("shadow_strength", shadow_strength);
   shader_p->set_int_for_uniform("shadow_softness", shadow_softness); // maybe unsigned int?

   shader_p->set_bool_for_uniform("show_shadows", show_just_shadows);


   // other Textures have already been bound in display_info_t::draw_models()
   glActiveTexture(GL_TEXTURE4);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_with_shadows() A4 " << err << std::endl;

   glBindTexture(GL_TEXTURE_2D, shadow_depthMap);

   glEnableVertexAttribArray(0);  // position
   glEnableVertexAttribArray(1);  // normal
   glEnableVertexAttribArray(2);  // tangent
   glEnableVertexAttribArray(3);  // bitangent
   glEnableVertexAttribArray(4);  // colour
   glEnableVertexAttribArray(5);  // texCoord

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_with_shadows() " << name << " pre-draw " << err << std::endl;

   // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
   // before making a new VAO?

   if (false)
      std::cout << "DEBUG:: TextureMesh::draw_with_shadows() " << name << " shader " << shader_p->name
                << " vao " << vao
                << " drawing " << n_verts << " triangle vertices"  << std::endl;

   glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_with_shadows() glDrawElements()"
                      << " of Mesh \"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   // if (use_blending) {
   // glDisable(GL_BLEND);
   // }

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);

   glUseProgram (0);

}

void
TextureMesh::draw_for_ssao(Shader *shader_p,
                           const glm::mat4 &model,
                           const glm::mat4 &view,
                           const glm::mat4 &projection) {

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;
   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_for_ssao() " << shader_p->name << " -- start -- "
                      << err << std::endl;

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   // vertex Shader:
   //
   // layout(location = 0) in vec3 position;
   // layout(location = 1) in vec3 normal;
   // layout(location = 2) in vec3 tangent;
   // layout(location = 3) in vec3 bitangent;
   // layout(location = 4) in vec4 colour;
   // layout(location = 5) in vec2 texCoord;

   // uniform mat4 model; // include the view rotation and item translation
   // uniform mat4 view; // the lookat matrix
   // uniform mat4 projection; // projection matrix

   shader_p->set_mat4_for_uniform("model",      model);
   shader_p->set_mat4_for_uniform("view",       view);
   shader_p->set_mat4_for_uniform("projection", projection);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_for_ssao() " << shader_name << " post uniforms" << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "TextureMesh::draw_for_ssao() You forgot to setup this mesh (or setup with empty vertices or triangles) "
                << "\"" << name << "\" \"" << shader_p->name << "\"" << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_for_ssao() \"" << shader_name << "\" \""
		      << name << "\"" << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   glEnableVertexAttribArray(0); // position
   glEnableVertexAttribArray(1); // normal
   glEnableVertexAttribArray(2); // tangent
   glEnableVertexAttribArray(3); // bitangent
   glEnableVertexAttribArray(4); // colour (not used)
   glEnableVertexAttribArray(5); // texture coordinates

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_ao() " << name << " pre-draw " << err << std::endl;

   glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: TextureMesh::draw_ao() glDrawElements() of Mesh "
                      << "\"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);

   glUseProgram(0);

}



void
TextureMesh::import(const std::vector<TextureMeshVertex> &verts_in, const std::vector<g_triangle> &triangles_in) {

   vertices = verts_in;
   triangles = triangles_in;
   draw_this_mesh = true;
}


void
TextureMesh::import(const IndexedModel &ind_model, float scale) {

   // std::cout << "TextureMesh::import(const IndexedModel &ind_model)" << std::endl;

   bool sane_input = false;
   if (ind_model.positions.size() == ind_model.texCoords.size())
      if (ind_model.positions.size() == ind_model.normals.size())
         sane_input = true;

   if (ind_model.positions.size() == 0)
      sane_input = false;

   std::cout << "TextureMesh::import() indices.size() " << ind_model.indices.size() << std::endl;

   if (sane_input) {
      for (unsigned int i=0; i<ind_model.positions.size(); i++) {
         glm::vec4 col(0.5, 0.5, 0.5, 1.0);
         // std::cout << "debug normal " << i << " " << glm::to_string(ind_model.normals[i]) << std::endl;
         // std::cout << "debug texCoords " << i << " " << glm::to_string(ind_model.texCoords[i]) << std::endl;
         TextureMeshVertex v(scale * ind_model.positions[i],
                             ind_model.normals[i],
                             col,
                             ind_model.texCoords[i]);
         vertices.push_back(v);
      }

      for (unsigned int i=0; i<ind_model.indices.size(); i += 3) {
         g_triangle gt(ind_model.indices[i],
                       ind_model.indices[i+1],
                       ind_model.indices[i+2]);
         triangles.push_back(gt);
      }
   }

   setup_buffers();

}


// 20220129-PE why did I delete this for crows?
// 20250310-PE Ah, possibly because it used glBufferData and instancing should
// update using glBufferSubData() - which this function now does.
void
TextureMesh::update_instancing_buffer_data(const std::vector<glm::vec3> &positions) {

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: in update_instancing_buffer_data(): You forgot to setup this TextureMesh "
                << name << std::endl;

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: TextureMesh::update_instancing_buffers() --- start --- " << _(err) << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: TextureMesh::setup_instancing_buffers() post binding vao " << _(err) << std::endl;

   // glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
   // err = glGetError();
   // if (err)
   //    std::cout << "GL ERROR:: TextureMesh::setup_instancing_buffers() post bind buffer " << _(err) << std::endl;

   int n_positions = positions.size();
   n_instances = n_positions;
   if (n_positions > n_instances_allocated)
      n_positions = n_instances_allocated;

   if (positions.empty()) return;

   // glBufferData(GL_ARRAY_BUFFER, n_positions * sizeof(glm::vec3), &(positions[0]), GL_STATIC_DRAW);

   glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id); // 20250312-PE needed.
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_positions * sizeof(glm::vec3), &(positions[0]));

}


// for happy faces that drift up the screen
void
TextureMesh::update_instancing_buffer_data_for_happy_faces(const std::vector<glm::vec3> &positions_in, // original positions
                                           unsigned int draw_count_in,
                                           unsigned int draw_count_max,
                                           const glm::vec3 &screen_y_uv) {

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: in update_instancing_buffer_data_for_happy_faces(): You forgot to setup this TextureMesh "
                << name << std::endl;

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: TextureMesh::update_instancing_buffer_data_for_happy_faces() --- start --- " << _(err) << std::endl;
   glBindVertexArray(vao);
   draw_count = draw_count_in; // do I need this to be a class data item?
   std::vector<glm::vec3> positions(positions_in);
   int n_positions = positions.size();
   if (n_positions > n_instances_allocated) {
      std::cout << "Too many TextureMesh instances " << n_positions << " " << n_instances_allocated
                << std::endl;
   } else {

      // Use the screen centre to generate tp.
      // Use the index of the position to change the phaase of the wiggle.

      auto get_position_delta = [draw_count_max] (const glm::vec3 &screen_y_uv,
                                                  unsigned int draw_count_in,
                                                  unsigned int index) {
                                   float f1 = static_cast<float>(draw_count_in)/static_cast<float>(draw_count_max);
                                   float f2 = f1 * f1 * 2.5f;
                                   glm::vec3 f_uv = f2 * screen_y_uv;
                                   glm::vec3 tp = glm::normalize(glm::vec3(0.1f, 0.2f, 0.3f));
                                   glm::vec3 cp_1 = glm::cross(screen_y_uv, tp);
                                   glm::vec3 cp_2 = glm::cross(screen_y_uv, cp_1);
                                   float phase = 0.1 * static_cast<float>(index);
                                   f_uv += 0.9f * sinf(9.0f * f1 + phase) * cp_2;
                                   return f_uv;
                                };

      n_instances = positions.size();

      // now update the positions
      for (unsigned int i=0; i<positions.size(); i++)
         positions[i] += get_position_delta(screen_y_uv, draw_count_in, i);

      glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_positions * sizeof(glm::vec3), &(positions[0]));
   }
}

void
TextureMesh::setup_instancing_buffers(unsigned int n_happy_faces_max) {

   n_instances = 0;
   n_instances_allocated = n_happy_faces_max;
   is_instanced = true;

   // in the layout this is 6 called "instance_translation"

   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "GL error ####"
                      << " TextureMesh::setup_instancing_buffers() A " << err << std::endl;
   unsigned int n_bytes = n_happy_faces_max * sizeof(glm::vec3);

   glGenBuffers(1, &inst_positions_id);
   // std::cout << "setup_instancing_buffers: glGenBuffers inst_positions_id " << inst_positions_id << std::endl;
   glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
   glBufferData(GL_ARRAY_BUFFER, n_bytes, nullptr, GL_DYNAMIC_DRAW);
   // prevous attributes are position, normal, colour, texCoords
   glEnableVertexAttribArray(6);

   void *step_over_previous_attrib_bytes = 0; // for instanced attributes, we don't need to step
                                              // over the standard vertex attributes.

   glVertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), step_over_previous_attrib_bytes);
   glVertexAttribDivisor(6, 1);
   err = glGetError();
   if (err) std::cout << "GL error #####"
                      << " TextureMesh::setup_instancing_buffers() B " << err << std::endl;
}


// draw_count and draw_count_max are used to set the opactity (it counts the number of times drawn)
//
void
TextureMesh::draw_instances(Shader *shader_p,
                            const glm::mat4 &mvp,
                            const glm::mat4 &view_rotation,
                            const glm::vec4 &background_colour,
                            bool is_perspective_projection) {

   if (false)
      std::cout << "TextureMesh::draw_instances() A " << name << " n_instances: " << n_instances
                << " n_triangles: " << triangles.size() <<std::endl;

   if (! draw_this_mesh) return;
   // this can happen when all the particles have life 0 - and have been removed.
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_instances() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // tangent // not used for camera-facing textures
   glEnableVertexAttribArray(3); // bitangent // not used
   glEnableVertexAttribArray(4); // colour
   glEnableVertexAttribArray(5); // texCoord
   glEnableVertexAttribArray(6); // instanced position

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_instances() " << shader_p->name
                      << " draw_instances() post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_instances() " << shader_p->name
                      << " draw_instances() post view_rotation uniform " << err << std::endl;

   shader_p->set_bool_for_uniform("is_perspective_projection", is_perspective_projection);
   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_float_for_uniform("opacity", 1.0);

   float scale = 1.0; // the shader says:  // 0.8 for happy faces, 0.2 for anchored/fixed atoms
   // the scale should be set in the arguments to setup_camera_facing_quad() function.
   shader_p->set_float_for_uniform("canvas_scale", scale);

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw_instances() activetexture "
                                          << err << std::endl;

   glEnable(GL_DEPTH_TEST);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   unsigned int n_verts = 6;

   if (false)
      std::cout << "TextureMesh::draw_instances() C " << name << " shader: " << shader_p->name << " "
                << "n_verts " <<  n_verts << " n_instances " << n_instances << std::endl;

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);

   err = glGetError();
   if (err) std::cout << "error draw_instances() on glDrawElementsInstanced() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);

}

void
TextureMesh::draw_instances_for_ssao(Shader *shader_p,
                                     const glm::mat4 &model,
                                     const glm::mat4 &view,
                                     const glm::mat4 &projection) {

   if (! draw_this_mesh) return;
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_instances() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // tangent // not used for camera-facing textures
   glEnableVertexAttribArray(3); // bitangent // not used
   glEnableVertexAttribArray(4); // colour
   glEnableVertexAttribArray(5); // texCoord
   glEnableVertexAttribArray(6); // instanced position


   shader_p->set_mat4_for_uniform("model",      model);
   shader_p->set_mat4_for_uniform("view",       view);
   shader_p->set_mat4_for_uniform("projection", projection);

   unsigned int n_verts = 6;
   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);


   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);

}



// draw_count and draw_count_max are used to set the opactity (it counts the number of times drawn)
//
void
TextureMesh::draw_fading_instances(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation,
                                   unsigned int draw_count, unsigned int draw_count_max) {

   // std::cout << "TextureMesh::draw_instances() A " << n_instances << " " << triangles.size()
   // <<std::endl;

   if (! draw_this_mesh) return;
   // this can happen when all the particles have life 0 - and have been removed.
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   const float pi = 3.1415926f;
   float draw_count_frac = static_cast<float>(draw_count)/static_cast<float>(draw_count_max);
   const float &f = draw_count_frac;  // shorthand
   float opacity = sinf(sqrt(f) * pi); // maybe this goes negative?
   // std::cout << "opacity: f " << f << " opacity " << opacity << std::endl;

   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_instances() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // tangent // not used for camera-facing textures
   glEnableVertexAttribArray(3); // bitangent // not used
   glEnableVertexAttribArray(4); // colour
   glEnableVertexAttribArray(5); // texCoord
   glEnableVertexAttribArray(6); // instanced position

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_instances() " << shader_p->name
                      << " draw_instances() post mvp uniform " << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: TextureMesh::draw_instances() " << shader_p->name
                      << " draw_instances() post view_rotation uniform " << err << std::endl;

   shader_p->set_float_for_uniform("opacity", opacity);

   float scale = 0.3; // the shader says:  // 0.8 for happy faces, 0.2 for anchored/fixed atoms
   shader_p->set_float_for_uniform("canvas_scale", scale);

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "error:: TextureMesh::draw_instances() activetexture "
                                          << err << std::endl;

   glEnable(GL_DEPTH_TEST); // not the problem
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   unsigned int n_verts = 6;

   if (false)
      std::cout << "TextureMesh::draw_instances() C " << name << " shader: " << shader_p->name << " "
                << "n_verts " <<  n_verts << " n_instances " << n_instances << std::endl;

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);

   err = glGetError();
   if (err) std::cout << "error draw_instances() on glDrawElementsInstanced() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);
   
}

void
TextureMesh::apply_scale(const float &s) {

   glm::vec3 sf(s, s, s);
   for (unsigned int ii=0; ii<vertices.size(); ii++)
      vertices[ii].position *= s;
   setup_buffers(); // transfer the coordinates
}


void
TextureMesh::translate(const glm::vec3 &t) {

   for (auto &v : vertices) {
      v.position += t;
   }
   setup_buffers();
}


void
TextureMesh::apply_transformation(const glm::mat4 &m) {

   for (unsigned int ii=0; ii<vertices.size(); ii++) {
      glm::vec4 p_1(vertices[ii].position, 1.0);
      glm::vec4 p_2 = p_1 * m;
      glm::vec3 p_3(p_2);
      vertices[ii].position = p_3;
   }
   setup_buffers(); // transfer the coordinates
}


// include_call_to_setup_buffers is by default true
bool
TextureMesh::load_from_glTF(const std::string &file_name_in, bool include_call_to_setup_buffers) {

   // is this what's in a node?
   class extracted_buffer_info_t {
   public:
      int buffer_view_index;
      std::vector<glm::vec3> positions;
      std::vector<glm::vec3> normals;
      std::vector<glm::vec2> texture_coords;
      std::vector<g_triangle> triangles;
      glm::vec4 base_colour;
   };


   auto indent = [] (unsigned int il) {
                   std::string s;
                   for (unsigned int i=0; i<il; i++) s += "   ";
                   return s;
                   };

   auto proc_indices = [indent] (tinygltf::Model &model, const tinygltf::Accessor &indices_accessor) {

                          std::vector<g_triangle> triangles;
                          std::cout << indent(3) << "proc_indices()" << std::endl;

                          std::vector<unsigned int> mesh_indices;
                          const tinygltf::BufferView &buffer_view = model.bufferViews[indices_accessor.bufferView];
                          const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                          const uint8_t *base = &buffer.data.at(buffer_view.byteOffset + indices_accessor.byteOffset);

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                             const uint32_t *p = reinterpret_cast<const uint32_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                             const uint16_t *p = reinterpret_cast<const uint16_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_BYTE) {
                             const uint8_t *p = reinterpret_cast<const uint8_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          unsigned int n_triangles = mesh_indices.size()/3;
                          if (n_triangles * 3 == mesh_indices.size()) {
                             for (unsigned int ii=0; ii<mesh_indices.size(); ii+=3) {
                                g_triangle t(mesh_indices[ii], mesh_indices[ii+1], mesh_indices[ii+2]);
                                triangles.push_back(t);
                             }
                          }
                          return triangles;
                       };

   auto proc_material = [] (tinygltf::Model &model, int material_index) {

                           std::cout << "debug:: proc_material(): material_index " << material_index
                                     << " materials.size() " << model.materials.size() << std::endl;

                           const tinygltf::Material &material = model.materials[material_index];

                           std::cout << "debug:: proc_material(): material.pbrMetallicRoughness() basecolour size"
                                     << material.pbrMetallicRoughness.baseColorFactor.size() << std::endl;

                           if (false)
                              std::cout << "Material " << material_index << " name: \"" << material.name << "\""
                                        << " alphaMode" << material.alphaMode
                                        << " pbr-metallicroughness colour "
                                        << material.pbrMetallicRoughness.baseColorFactor[0] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[1] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[2] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[3] << " "
                                        << " metalicFactor " << material.pbrMetallicRoughness.metallicFactor
                                        << " roughnessFactor " << material.pbrMetallicRoughness.roughnessFactor
                                        << std::endl;
                           glm::vec4 colour(material.pbrMetallicRoughness.baseColorFactor[0],
                                            material.pbrMetallicRoughness.baseColorFactor[1],
                                            material.pbrMetallicRoughness.baseColorFactor[2],
                                            material.pbrMetallicRoughness.baseColorFactor[3]);
                           return colour;
                   };

   auto proc_primitive = [proc_indices, proc_material, indent] (tinygltf::Model &model, const tinygltf::Primitive &primitive) {

                            std::cout << indent(2) << "proc_primitive()" << std::endl;

                            extracted_buffer_info_t ebi; // returned object

                            // ------------------- primitive material -------------------

                            std::cout << indent(3) << "primitive material " << primitive.material << std::endl;
                            if (primitive.material >= 0) {
                               glm::vec4 colour = proc_material(model, primitive.material);
                               ebi.base_colour = colour;
                            } else {
                               std::cout << "Ooops skipping proc_material()" << std::endl;
                            }

                            // ------------------- primitive attributes -------------------

                            std::map<std::string, int>::const_iterator it;
                            unsigned int icount = 0;
                            for (it=primitive.attributes.begin(); it!=primitive.attributes.end(); ++it) {
                               const std::string &key = it->first;
                               std::cout << "::: key " << std::setw(8) << key << " for attribute index "
                                         << icount << " of " << primitive.attributes.size() << std::endl;
                               icount++;
                            }

                            for (it=primitive.attributes.begin(); it!=primitive.attributes.end(); ++it) {
                               const std::string &key = it->first;
                               int index = it->second;
                               std::cout << indent(3) << "primvitive attribute " << key << " " << index << std::endl;

                               const tinygltf::Accessor &accessor = model.accessors[index];
                               const tinygltf::BufferView &buffer_view = model.bufferViews[accessor.bufferView];
                               const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER)
                                  std::cout << indent(3) << "an array buffer" << std::endl;
                               if (buffer_view.target == TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER)
                                  std::cout << indent(3) << "an index buffer" << std::endl;

                               if (false)
                                  std::cout << indent(2) << "accessor componentType " << accessor.componentType
                                            << " accessor type " << accessor.type << std::endl;

                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_BYTE)
                                  std::cout << "      component type byte " << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE)
                                  std::cout << "      component type unsigned byte " << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_SHORT)
                                  std::cout << "      component type short" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT)
                                  std::cout << "      component type unsigned short" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_INT)
                                  std::cout << "      component type int" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT)
                                  std::cout << "      component type unsigned int" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT)
                                  std::cout << "      component type float" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE)
                                  std::cout << "      component type double" << std::endl;

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "POSITION") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                           // std::cout << "      ---- float array for positions " << std::endl;
                                           const float *positions = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = positions[ii*3 + 0];
                                              float y = positions[ii*3 + 1];
                                              float z = positions[ii*3 + 2];
                                              ebi.positions.push_back(glm::vec3(x,y,z));
                                           }
                                        }
                                     }
                                  }
                               }

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "NORMAL") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                           // std::cout << "      ---- float array for normals " << std::endl;
                                           const float *normals = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = normals[ii*3 + 0];
                                              float y = normals[ii*3 + 1];
                                              float z = normals[ii*3 + 2];
                                              ebi.normals.push_back(glm::vec3(x,y,z));
                                           }
                                        }
                                     }
                                  }
                               }

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "TEXCOORD_0") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC2) {
                                           std::cout << "      ---- float array for texture coordinates " << std::endl;
                                           const float *tcs = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = tcs[ii*2 + 0];
                                              float y = tcs[ii*2 + 1];
                                              ebi.texture_coords.push_back(glm::vec2(x,y));
                                           }
                                        }
                                     }
                                  }
                               }
                            }

                            // ------------------- primitive indices -------------------

                            const tinygltf::Accessor &indices_accessor = model.accessors[primitive.indices];
                            std::vector<g_triangle> triangles = proc_indices(model, indices_accessor);
                            ebi.triangles = triangles; // std::move here?

                            std::cout << indent(2) << "debug:: proc_primitive() returns ebi with " << ebi.positions.size() << " positions "
                                      << ebi.normals.size() << " normals " << ebi.texture_coords.size() << " texture-coords and "
                                      << triangles.size() << " triangles" << std::endl;

                            return ebi;
                         };

   auto proc_mesh = [proc_primitive, indent] (tinygltf::Model &model, const tinygltf::Mesh &mesh) {

                       std::cout << indent(1) << "proc_mesh(): " << mesh.name << std::endl;
                       std::vector<extracted_buffer_info_t> ebi_vec;
                       for (unsigned int i=0; i<mesh.primitives.size(); i++) {
                          std::cout << indent(1) << "proc_mesh(): primitive number " << i << " of " << mesh.primitives.size() << std::endl;
                          const tinygltf::Primitive &primitive = mesh.primitives[i];
                          if (primitive.indices >= 0) {
                             auto new_mesh_primitive = proc_primitive(model, primitive);
                             ebi_vec.push_back(new_mesh_primitive);
                          }
                       }
                       return ebi_vec;
                   };


   auto proc_node = [proc_mesh] (tinygltf::Model &model, size_t i_node_index, const tinygltf::Node &node) {

                       std::vector<extracted_buffer_info_t> r;

                       std::cout << "Node " << i_node_index << " info:" << std::endl;
                       std::cout << "   Node name: " << node.name << std::endl;
                       std::cout << "   Node mesh index: " << node.mesh << std::endl;
                       std::cout << "   Node rotation vec elements count: " << node.rotation.size() << std::endl;
                       std::cout << "   Node scale vec elements count:    " << node.scale.size() << std::endl;
                       std::cout << "   Node scale vec translation count: " << node.translation.size() << std::endl;
                       if (node.mesh > -1) { // not a light or a camera
                          r = proc_mesh(model, model.meshes[node.mesh]);
                       } else {
                       }
                       return r;
                   };

   auto proc_images = [indent] (const tinygltf::Image &image) {
                         std::cout << indent(1) << "Image name: "   << image.name   << std::endl;
                         std::cout << indent(1) << "Image bits: "   << image.bits   << std::endl;
                         std::cout << indent(1) << "Image width: "  << image.width  << std::endl;
                         std::cout << indent(1) << "Image height: " << image.height << std::endl;
                         std::cout << indent(1) << "Image component: " << image.component << std::endl;
                         std::cout << indent(1) << "Image data buffer size:" << image.image.size() << std::endl;

                         Texture texture;
                         texture.handle_raw_image_data(image.name, image.image, image.width, image.height);

                         return texture;
                   };

   auto proc_model_v2 = [proc_node, proc_images] (tinygltf::Model &model) {

                           std::vector<extracted_buffer_info_t> r;
                           std::vector<Texture> textures;

                           unsigned int n_accessors = model.accessors.size();
                           std::cout << "debug:: model has " << n_accessors << " accessors" << std::endl;

                           // --- Materials ---

                           std::cout << "INFO:: this model contains " << model.materials.size() << " materials" << std::endl;
                           for (unsigned int imat=0; imat<model.materials.size(); imat++) {
                              const tinygltf::Material &material = model.materials[imat];
                              std::cout << "Material " << imat << " name: " << material.name << " alphaMode" << material.alphaMode
                                        << " pbr-metallicroughness colour "
                                        << material.pbrMetallicRoughness.baseColorFactor[0] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[1] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[2] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[3] << " "
                                        << " metalicFactor " << material.pbrMetallicRoughness.metallicFactor
                                        << " roughnessFactor " << material.pbrMetallicRoughness.roughnessFactor
                                        << std::endl;
                           }

                           // --- Vertices and Indices ---
#if 0
                           // scenes

                           for (unsigned int iscene=0; iscene<model.scenes.size(); iscene++) {
                              const tinygltf::Scene &scene = model.scenes[iscene];
                              std::cout << ":::: This scene contains " << scene.nodes.size() << " nodes "<< std::endl;
                              for (size_t i = 0; i < scene.nodes.size(); i++) {
                                 auto nodes = proc_node(model, i, model.nodes[scene.nodes[i]]);
                                 r.insert(r.end(), nodes.begin(), nodes.end());
                              }
                           }
#endif

                           // model

                           std::cout << "This model contains " << model.nodes.size() << " nodes" << std::endl;

                           for (size_t i = 0; i < model.nodes.size(); i++) {
                              auto nodes = proc_node(model, i, model.nodes[i]);
                              r.insert(r.end(), nodes.begin(), nodes.end());
                           }

                           // images

                           std::cout << "This model contains " << model.images.size() << " images and "
                                     << model.textures.size() << " textures" << std::endl;
                           for (unsigned int i=0; i<model.textures.size(); i++) {
                              int image_index = model.textures[i].source;
                              if (image_index >= 0 && image_index < static_cast<int>(model.images.size())) {
                                 Texture texture = proc_images(model.images[image_index]);
                                 textures.push_back(texture);
                              }
                           }

                           return std::make_pair(r, textures);
                        };

   bool status = true;
   tinygltf::Model model;
   tinygltf::TinyGLTF loader;
   std::string err;
   std::string warn;

   std::string file_name = file_name_in;
   if (coot::file_exists(file_name)) {
      // do nothing
   } else {
      std::string dir = coot::package_data_dir();
      std::string dir_2 = coot::util::append_dir_dir(dir, "glTF");
      file_name = coot::util::append_dir_file(dir_2, file_name);
   }

   bool use_binary = false;
   std::string ext = coot::util::file_name_extension(file_name_in);
   if (ext == ".glb") use_binary = true;

   // std::cout << "debug:: TextureMesh::load_from_glTF(): " << file_name_in << " use_binary: " << use_binary << std::endl;

   // use the extension to check which function to use
   //
   // bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
   bool read_success = false;

   if (use_binary) {
      read_success = loader.LoadBinaryFromFile(&model, &err, &warn, file_name); // for binary glTF(.glb)
   } else {
      read_success = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
   }

   if (!warn.empty()) {
      std::cout << "WARNING:: load_from_glTF():" << warn << std::endl;
   }

   if (!err.empty()) {
      std::cout << "ERROR:: load_from_glTF(): " << err << std::endl;
   }

   if (read_success == false) {
      std::cout << "WARNING:: failed to parse glTF from " << file_name << std::endl;
   }

   if (err.empty() && read_success) { // OK, good, let's go!

      auto r_pair = proc_model_v2(model);
      std::vector<extracted_buffer_info_t> &r = r_pair.first;
      std::vector<Texture> &extracted_textures = r_pair.second;

      for (const auto &texture : extracted_textures) {
         std::string texture_name("texture-name-here");
         if (!texture.file_name.empty())
            texture_name = texture.file_name;
         TextureInfoType ti(texture, texture_name);
         textures.push_back(ti);
      }

      // std::cout << "load_from_glTF(): found " << r.size() << " mesh primitives" << std::endl;
      for (unsigned int i=0; i<r.size(); i++) {
         const extracted_buffer_info_t &ebi = r[i];
         if (ebi.normals.size() > 0) {
            if (ebi.normals.size() == ebi.positions.size()) {
               unsigned int max_vertex_index = 0; // set this
               if (max_vertex_index < ebi.normals.size()) {
                  unsigned int idx_vert_base = vertices.size();
                  unsigned int idx_tri_base = triangles.size();
                  if (false)
                     std::cout << "debug ebi sizes " << ebi.positions.size() << " "
                               << ebi.normals.size() << " " << ebi.texture_coords.size() << std::endl;
                  if (ebi.texture_coords.size() == ebi.positions.size()) {
                     for (unsigned int j=0; j<ebi.normals.size(); j++) {
                        // s_generic_vertex g(ebi.positions[j], ebi.normals[j], ebi.base_colour);
                        TextureMeshVertex tmv(ebi.positions[j], ebi.normals[j], ebi.base_colour, ebi.texture_coords[j]);
                        vertices.push_back(tmv);
                     }
                  }

                  // 20211011-PE If ebi.texture_coords.size() is 0 - as we would expect for a untextured model, then perhaps we can
                  // just add a dummy texture_coords and set a flag to say to use the colours and not the textures.
                  // Another time.

                  triangles.insert(triangles.end(), ebi.triangles.begin(), ebi.triangles.end());
                  if (idx_vert_base != 0) {
                     for (unsigned int i=idx_tri_base; i<triangles.size(); i++) {
                        triangles[i].rebase(idx_vert_base);
                     }
                  }
               }
            }
         }
      }

      if (include_call_to_setup_buffers) {
         std::cout << "pre-setup_buffer() " << vertices.size() << " vertices "  << triangles.size() << " triangles "  << std::endl;
         setup_buffers();
      }
   } else {
      std::cout << "::::::::::::: non-success path" << std::endl;
      status = false; // boo
   }

   std::cout << "debug:: load_from_glTF() " << file_name_in << " " << vertices.size() << " vertices "
             << triangles.size() << " triangles "  << std::endl;
   std::cout << "debug:: load_from_glTF() returns status " << status << std::endl;
   return status;
}

