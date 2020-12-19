
#include "HUDTextureMesh.hh"

#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()


void
HUDTextureMesh::init() {

   draw_this_mesh = true;
   position = glm::vec2(0,0);
   scales = glm::vec2(1,1);
   vao = VAO_NOT_SET; // unset
   first_time = true;
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

   set_position_and_scales(glm::vec2(0,0), glm::vec2(1.0, 1.0));

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
HUDTextureMesh::set_position_and_scales(const glm::vec2 &pos, const glm::vec2 &scales_in) {

   position =  pos;
   scales = scales_in;
}

void
HUDTextureMesh::set_position(const glm::vec2 &pos) {

   position =  pos;
}

void
HUDTextureMesh::set_scales(const glm::vec2 &scales_in) {

   scales = scales_in;
}

void
HUDTextureMesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

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
   glDisableVertexAttribArray(1);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

   first_time = false;

}

void
HUDTextureMesh::draw(Shader *shader_p) {

   if (! draw_this_mesh) return;

   shader_p->Use();

   if (vao == VAO_NOT_SET)
      std::cout << "error:: You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;

   glBindVertexArray(vao);

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id); // needed?
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id); // needed?

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);

   if (false) {
      std::cout << "HUDTextureMesh::draw() " << name << " sending position " << glm::to_string(position) << std::endl;
      std::cout << "HUDTextureMesh::draw() " << name << " sending scales "   << glm::to_string(scales) << std::endl;
   }
   shader_p->set_vec2_for_uniform("position", position);
   shader_p->set_vec2_for_uniform("scales", scales);

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

void
HUDTextureMesh::draw_label(const std::string &label, glm::vec4 &text_colour, Shader *shader_p,
                           const std::map<GLchar, FT_character> &ft_characters) {

   if (! draw_this_mesh) return;
   unsigned int n_triangles = triangles.size();
   unsigned int n_vertices = vertices.size();
   if (n_triangles == 0) return;

   // std::cout << "we didn't return early" << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "error draw_label() " << shader_p->name << " -- start -- " << err << std::endl;

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   if (vao == VAO_NOT_SET)
      std::cout << "error:: draw_label(): You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error HUDTextureMesh::draw_label()) " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err "
                      << err << std::endl;

   shader_p->set_vec2_for_uniform("position", position);
   shader_p->set_vec2_for_uniform("scales", scales);
   shader_p->set_vec4_for_uniform("text_colour", text_colour);

   glActiveTexture(GL_TEXTURE0);
   err = glGetError();
   if (err) std::cout << "error:: HUDTextureMesh::draw_atom_label() A3 " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id); // needed?
   err = glGetError(); if (err) std::cout << "error HUDTextureMesh::draw_label() glBindBuffer() v "
                                          << err << std::endl;

   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id); // needed?
   err = glGetError(); if (err) std::cout << "error HUDTextureMesh::draw_label() glBindBuffer() i "
                                          << err << std::endl;

   glEnableVertexAttribArray(0);  // vertex
   glEnableVertexAttribArray(1);  // texCoord

   // extract this and put it in its own file and function? Hmm tricky, because of the use of vertices.
   // More trouble than it's worth?

   float x = 0;
   float y = 0;
   GLfloat scale = 10.1; // guess
   std::string::const_iterator it_c;
   for (it_c = label.begin(); it_c != label.end(); ++it_c) {
      std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*it_c);
      if (it == ft_characters.end()) {
         std::cout << "ERROR:: HUDTextureMesh::draw_label() Failed to lookup glyph for " << *it_c << std::endl;
         continue;
      };
      const FT_character &ch = it->second;
      GLfloat xpos = x + ch.Bearing.x * scale;
      GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;
      GLfloat w = ch.Size.x * scale;
      GLfloat h = ch.Size.y * scale;

      // Update the vertices for each character
      //
      std::vector<HUDTextureMesh_attribs_t> texture_mesh_vertices = vertices;

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
         texture_mesh_vertices[i].position += glm::vec2(xpos, ypos);

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
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(HUDTextureMesh_attribs_t), &texture_mesh_vertices[0]);

      unsigned int n_draw_verts = 6;
      glDrawElements(GL_TRIANGLES, n_draw_verts, GL_UNSIGNED_INT, nullptr);

      err = glGetError(); if (err) std::cout << "error HUDTextureMesh::draw_label() glDrawElements() " << err << std::endl;

       // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
      x += (ch.Advance >> 6) * scale * 1.0;

   }


   // ------------------------------- done text texture code  ----------------------

   err = glGetError();
   if (err) std::cout << "   error HUDTextureMesh::draw_label() glDrawElements()"
                      << " of \"" << name << "\""
                      << " shader: " << shader_name
                      << " vao " << vao
                      << " n_vertices " << n_vertices
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glUseProgram (0);

}


void
HUDTextureMesh::draw_label(const std::string &label, bool highlight_label_flag, Shader *shader_p,
                           const std::map<GLchar, FT_character> &ft_characters) {
   glm::vec4 text_colour(0.8, 0.8, 0.8, 1.0);
   if (highlight_label_flag) text_colour = glm::vec4(1.0, 1.0, 0.6, 1.0);
   draw_label(label, text_colour, shader_p, ft_characters);

}
