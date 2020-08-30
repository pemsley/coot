
#include <iostream>
#include "HUDMesh.hh"

void
HUDMesh::init() {

   n_instances = 0;
   first_time = true;
   use_blending = false;
   inst_hud_bar_attribs_buffer_id = 0;
}

void
HUDMesh::setup_camera_facing_quad() {

   float scale_x = 0.4; // pass?
   float scale_y = 0.2;

   glm::vec3 n(0,0,1);
   glm::vec4 col(1.0, 1.0, 1.0, 1.0);

   vertices.clear();
   triangles.clear();

   vertices.push_back(glm::vec2(0.0f, 0.0f ));
   vertices.push_back(glm::vec2(1.0f, 0.0f ));
   vertices.push_back(glm::vec2(1.0f, 0.03f));
   vertices.push_back(glm::vec2(0.0f, 0.03f));

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));

   setup_buffers();

}

void
HUDMesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   if (first_time)
      glGenVertexArrays(1, &vao);

   glBindVertexArray(vao);

   unsigned int n_vertices = vertices.size(); // 4

   if (first_time) {
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(glm::vec2), &(vertices[0]), GL_DYNAMIC_DRAW);
   } else {
      glDeleteBuffers(1, &buffer_id);
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(glm::vec2), &(vertices[0]), GL_DYNAMIC_DRAW);
   }

   // position (of the quad)
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), 0);

   // 2 triangles - it's a quad - instanced.
   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes = n_triangles * 3 * sizeof(unsigned int);

   if (first_time) {
      glGenBuffers(1, &index_buffer_id);
      GLenum err = glGetError(); if (err) std::cout << "GL error HUDMesh setup_buffers()\n";
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
HUDMesh::setup_instancing_buffer(unsigned int n_boxes) {

   // make *space* for the instancing values, but don't fill them with data (here)

   n_instances = n_boxes; // this is how many we want to draw now.

   glBindVertexArray (vao);

   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);

   glGenBuffers(1, &inst_hud_bar_attribs_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_hud_bar_attribs_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_boxes * sizeof(HUD_bar_attribs_t), nullptr, GL_DYNAMIC_DRAW);

   // layout
   // 0 vec2 position        (quad)
   // 1 vec4 colour          (instanced)
   // 2 vec2 position_offset (instanced)
   // 3 float scale (width)  (instanced)

   // colour - instanced
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(HUD_bar_attribs_t), 0);
   glVertexAttribDivisor(1, 1);

   // position_offset - instanced
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(HUD_bar_attribs_t),
                         reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(2, 1);

   // scale (width of the bar - along the x axis) - instanced
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(HUD_bar_attribs_t),
                         reinterpret_cast<void *>(sizeof(glm::vec4) + sizeof(glm::vec2)));
   glVertexAttribDivisor(3, 1);

   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);

   glBindBuffer(GL_ARRAY_BUFFER, 0);

}

void
HUDMesh::update_instancing_buffer_data(const std::vector<HUD_bar_attribs_t> &new_bars) {

   if (false)
      std::cout << "debug:: HUDMesh::update_instancing_buffer_data() "
                << new_bars.size() << " " << n_instances << std::endl;
   if (!new_bars.empty()) {
      unsigned int s = new_bars.size();
      if (s > n_instances)
         s = n_instances; // only draw as many as we have allocated
      else
         n_instances = s;
      glBindBuffer(GL_ARRAY_BUFFER, inst_hud_bar_attribs_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, s * sizeof(HUD_bar_attribs_t), &(new_bars[0]));
   }
}

void
HUDMesh::draw(Shader *shader_p) {

   // std::cout << "debug:: HUDMesh::draw() --- start --- " << n_instances << std::endl;

   if (this_mesh_is_closed) return;

   if (n_instances == 0) return;

   shader_p->Use();

   glBindVertexArray(vao);
   // glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glEnableVertexAttribArray(0);
   // glBindBuffer(GL_ARRAY_BUFFER, inst_hud_bar_attribs_buffer_id);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);

   glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr, n_instances);
   GLenum err = glGetError();
   if (err) std::cout << "   error HUDMesh::draw() glDrawElementsInstanced()"
                      << " of HUDMesh \"" << name << "\""
                      << " with shader" << shader_p->name
                      << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glUseProgram(0);

}
