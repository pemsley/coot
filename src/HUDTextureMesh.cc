
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
   is_instanced = false;
   window_resize_scales_correction_set = false;
   window_resize_position_correction_set = false;
   window_resize_position_correction = glm::vec2(0,0);
   window_resize_scales_correction   = glm::vec2(1,1);
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

   vertices[2].texture_coords.y = 1.0; // was 0.5;
   vertices[3].texture_coords.y = 1.0; // was 0.5;
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
HUDTextureMesh::update_instancing_buffer_data(const std::vector<glm::vec2> &new_positions) {

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: HUDTextureMesh::update_instancing_buffer_data() --start-- err "
                << err << std::endl;

   unsigned int n_phi_psis = new_positions.size();
   if (n_phi_psis > n_instances_max)
      n_phi_psis = n_instances_max;

   n_instances = n_phi_psis;

   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: HUDTextureMesh::update_instancing_buffer_data() You forgot to setup this Mesh "
                << name << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: HUDTextureMesh::update_instancing_buffer_data() binding vao err "
                << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_phi_psis * sizeof(glm::vec2), &(new_positions[0]));
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: HUDTextureMesh::update_instancing_buffer_data() binding buffersubdata err "
                << err << std::endl;

}


void
HUDTextureMesh::draw(Shader *shader_p, screen_position_origins_t screen_position_origin) {

   // 20220308-PE if you are looking at rama data points, then they are drawn using draw_instances()

   if (false)
      std::cout << "HUDTextureMesh::draw() " << name << " " << shader_p->name
                << " draw_this_mesh " << draw_this_mesh << std::endl;

   if (! draw_this_mesh) return;

   if (is_instanced) {
      std::cout << "GL ERROR:: wrong draw call in HUDTextureMesh::draw()" << std::endl;
      return;
   }

   shader_p->Use();

   if (vao == VAO_NOT_SET)
      std::cout << "error:: You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl;

   glBindVertexArray(vao);

   glEnableVertexAttribArray(0);   // vec2 vertices
   glEnableVertexAttribArray(1);   // vec2 texCoords;

   glDisable(GL_DEPTH_TEST);

   bool rel_top   = false;
   bool rel_right = false;
   if (screen_position_origin == TOP_RIGHT) {
      rel_top   = true;
      rel_right = true;
   }
   if (screen_position_origin == BOTTOM_RIGHT) {
      rel_right = true;
   }
   if (screen_position_origin == TOP_LEFT) {
      rel_top   = true;
   }
   shader_p->set_bool_for_uniform("relative_to_top",   rel_top);
   shader_p->set_bool_for_uniform("relative_to_right", rel_right);

   glm::vec4 text_colour(0.8, 0.7, 0.5, 1.0);                 // what's this used for?
   shader_p->set_vec2_for_uniform("position", position);
   shader_p->set_vec2_for_uniform("scales",   scales);
   shader_p->set_vec4_for_uniform("text_colour", text_colour);

   shader_p->set_int_for_uniform("image_texture", 0); // sampler2D - we don't *need* to set this uniform for the sampler
                                                      // (as there are only 1 sampler2Ds) but it's good practice to do
                                                      // so, I suppose

   if (window_resize_position_correction_set)
      shader_p->set_vec2_for_uniform("window_resize_position_correction", window_resize_position_correction);
   if (window_resize_scales_correction_set)
      shader_p->set_vec2_for_uniform("window_resize_scales_correction", window_resize_scales_correction);

   if (false) {
      std::cout << "HUDTextureMesh::draw() " << name << " sending"
                << " scales " << glm::to_string(scales)
                << " position " << glm::to_string(position)
                << " window_resize_scales_correction "   << glm::to_string(window_resize_scales_correction)
                << " window_resize_position_correction " << glm::to_string(window_resize_position_correction)
                << std::endl;
   }

   // std::cout << "debug:: HUDTextureMesh::draw() glDrawElements()" << std::endl;
   glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: HUDMesh::draw() glDrawElementsInstanced()"
                      << " of HUDMesh \"" << name << "\""
                      << " with shader" << shader_p->name
                      << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glUseProgram(0);

}

void
HUDTextureMesh::draw_instances(Shader *shader_p) {

   if (! draw_this_mesh) return;

   if (! is_instanced) { // set by setup_instancing_buffers()
      std::cout << "GL ERROR:: wrong draw call in HUDTextureMesh::draw_instances()" << std::endl;
      return;
   }

   shader_p->Use();
   if (vao == VAO_NOT_SET)
      std::cout << "error:: You forgot to setup this mesh " << name << " "
                << shader_p->name << std::endl; // or maybe you didn't attach_buffers() before
                                                // creating this mesh
   glBindVertexArray(vao);

   glEnableVertexAttribArray(0);  // vec2 vertices   - standard attribute
   glEnableVertexAttribArray(1);  // vec2 texCoords  - standard attribute
   glEnableVertexAttribArray(2); // position of this (instanced) phi_psi

   shader_p->set_vec2_for_uniform("position", position); // these are for position global positioning and
   shader_p->set_vec2_for_uniform("scales", scales);     // scaling - used by all rama points

   if (window_resize_position_correction_set)
      shader_p->set_vec2_for_uniform("window_resize_position_correction", window_resize_position_correction);
   if (window_resize_scales_correction_set)
      shader_p->set_vec2_for_uniform("window_resize_scales_correction", window_resize_scales_correction);

   if (false) {
      std::cout << "HUDTextureMesh::draw_instances() " << name << " sending"
                << " scales " << glm::to_string(scales)
                << " position " << glm::to_string(position)
                << " window_resize_scales_correction "   << glm::to_string(window_resize_scales_correction)
                << " window_resize_position_correction " << glm::to_string(window_resize_position_correction)
                << std::endl;
   }

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERORR:: in HUDTextureMesh::draw_instances() err " << err << std::endl;

   unsigned int n_verts = 6; // 2 triangles

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);

   err = glGetError();
   if (err) std::cout << "error HUDMesh::draw() glDrawElementsInstanced()"
                      << " of HUDMesh \"" << name << "\""
                      << " with shader" << shader_p->name
                      << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glUseProgram(0);
}

void
HUDTextureMesh::setup_instancing_buffers(unsigned int n_phi_psi_max) {

   n_instances = 0;
   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: HUDTextureMesh::setup_instancing_buffers() You forgot to setup this mesh "
                << name << std::endl;
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERORR:: in HUDTextureMesh::setup_instancing_buffers() err  " << err
                << " on binding vao " << vao << std::endl;
   is_instanced = true;
   n_instances_max = n_phi_psi_max;
   unsigned int n_bytes = n_phi_psi_max * sizeof(glm::vec2);
   glGenBuffers(1, &inst_positions_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_positions_id);
   glBufferData(GL_ARRAY_BUFFER, n_bytes, nullptr, GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), nullptr);
   glVertexAttribDivisor(2, 1);
   err = glGetError();
   if (err)
      std::cout << "GL ERORR:: in HUDTextureMesh::setup_instancing_buffers() err " << err << std::endl;

}

float
HUDTextureMesh::get_sum_x_advance(const std::string &label, const std::map<GLchar, FT_character> &ft_characters) const {

   float x_advance_sum = 0.0;
   std::string::const_iterator it_c;
   GLfloat scale = 10.1; // Hmmmmm!
   for (it_c = label.begin(); it_c != label.end(); ++it_c) {
      std::map<GLchar, FT_character>::const_iterator it = ft_characters.find(*it_c);
      if (it == ft_characters.end()) {
         std::cout << "ERROR:: HUDTextureMesh::draw_label() Failed to lookup glyph for " << *it_c << std::endl;
         continue;
      };
      const FT_character &ch = it->second;
      x_advance_sum += (ch.Advance >> 6) * scale;
   }
   return x_advance_sum;
}

void
HUDTextureMesh::draw_label(const std::string &label, glm::vec4 &text_colour, Shader *shader_p,
                           const std::map<GLchar, FT_character> &ft_characters) {

   if (! draw_this_mesh) return;
   unsigned int n_triangles = triangles.size();
   unsigned int n_vertices = vertices.size();
   if (n_triangles == 0) return;

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

   if (false) {
      std::cout << "sending scales " << glm::to_string(scales) << std::endl;
      std::cout << "sending position " << glm::to_string(position) << std::endl;
   }
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
   scale = 9.0; // testing
   std::string::const_iterator it_c;
   for (it_c = label.begin(); it_c != label.end(); ++it_c) {

      char cc = *it_c;

      err = glGetError();
      if (err)
         std::cout << "error HUDTextureMesh::draw_label() glDrawElements() 0 loop start " << err <<
            " when drawing char for " << cc << std::endl;
      
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

      err = glGetError();
      if (err)
         std::cout << "error HUDTextureMesh::draw_label() glDrawElements() --pre-bind texture " << err << " when drawing char for " << cc << std::endl;
      glBindTexture(GL_TEXTURE_2D, ch.TextureID);
      err = glGetError();
      if (err)
         std::cout << "error HUDTextureMesh::draw_label() glDrawElements() --post-bind texture " << err << " when drawing char for " << cc
                   << " with ch.TextureID " << ch.TextureID << std::endl;
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      err = glGetError();
      if (err)
         std::cout << "error HUDTextureMesh::draw_label() glDrawElements() B " << err << " when drawing char for " << cc << std::endl;
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_vertices * sizeof(HUDTextureMesh_attribs_t), &texture_mesh_vertices[0]);
      err = glGetError();
      if (err)
         std::cout << "error HUDTextureMesh::draw_label() glDrawElements() C " << err << " when drawing char for " << cc << std::endl;
      unsigned int n_draw_verts = 6;
      glDrawElements(GL_TRIANGLES, n_draw_verts, GL_UNSIGNED_INT, nullptr);
      err = glGetError();
      if (err) {
         std::cout << "error HUDTextureMesh::draw_label() glDrawElements() D " << err << " when drawing char for " << cc << std::endl;
      }

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
   glDisable(GL_BLEND);
   glUseProgram (0);

}


void
HUDTextureMesh::draw_label(const std::string &label, bool highlight_label_flag, Shader *shader_p,
                           const std::map<GLchar, FT_character> &ft_characters) {
   glm::vec4 text_colour(0.8, 0.8, 0.8, 1.0);
   if (highlight_label_flag) text_colour = glm::vec4(1.0, 1.0, 0.6, 1.0);
   draw_label(label, text_colour, shader_p, ft_characters);

}
