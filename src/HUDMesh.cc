
#include <iostream>
#include "HUDMesh.hh"

void
HUDMesh::init() {

   max_n_instances = 0;
   n_instances = 0;
   first_time = true;
   use_blending = false;
   inst_hud_bar_attribs_buffer_id = 0;
}

void
HUDMesh::setup_camera_facing_quad_for_bar() {

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
HUDMesh::setup_vertices_and_triangles_for_button() {

   vertices.clear();
   triangles.clear();

   vertices.push_back(glm::vec2(0.0f, 0.0f));
   vertices.push_back(glm::vec2(1.0f, 0.0f));
   vertices.push_back(glm::vec2(1.0f, 1.0f));
   vertices.push_back(glm::vec2(0.0f, 1.0f));

   // Add more vertices (and triangles) here for rounded edge look

   float delta = 0.03;
   float r = 3.0; // because buttons are wider than the height
   vertices.push_back(glm::vec2(   -delta,     delta*r));
   vertices.push_back(glm::vec2(   -delta, 1.0-delta*r));
   vertices.push_back(glm::vec2(1.0+delta, 1.0-delta*r));
   vertices.push_back(glm::vec2(1.0+delta,     delta*r));

   // and the corresponding shades:
   shades.push_back(-1);
   shades.push_back(-1);
   shades.push_back(1);
   shades.push_back(1);
   shades.push_back(-1.0 + 2.0 * delta*r);
   shades.push_back( 1.0 - 2.0 * delta*r);
   shades.push_back( 1.0 - 2.0 * delta*r);
   shades.push_back(-1.0 + 2.0 * delta*r);

   // triangles:

   triangles.push_back(g_triangle(0,4,5));
   triangles.push_back(g_triangle(5,3,0));
   triangles.push_back(g_triangle(1,6,7));
   triangles.push_back(g_triangle(1,2,6));

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));
   //

   setup_buffers();

   // now caller of this should now call setup_instancing_buffer(n_buttons_max)
}

void
HUDMesh::setup_buffers() {

   if (triangles.empty()) return;
   if (vertices.empty()) return;

   if (first_time)
      glGenVertexArrays(1, &vao);

   glBindVertexArray(vao);

   // ------------------------- vertices -----------------------------------------

   unsigned int n_vertices = vertices.size(); // 4 for hud bars

   // if (shades.size() != vertices.size())
      shades.resize(vertices.size(), 0.0);

   if (first_time) {
      glGenBuffers(1, &vertex_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(glm::vec2), &(vertices[0]), GL_DYNAMIC_DRAW);
      glGenBuffers(1, &shades_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, shades_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(float), &(shades[0]), GL_DYNAMIC_DRAW);
   } else {
      glDeleteBuffers(1, &vertex_buffer_id);
      glGenBuffers(1, &vertex_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(glm::vec2), &(vertices[0]), GL_DYNAMIC_DRAW);
      glDeleteBuffers(1, &shades_buffer_id);
      glGenBuffers(1, &shades_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, shades_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(float), &(shades[0]), GL_DYNAMIC_DRAW);
   }

   // vertices (of the quad)
   glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_id);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), 0);
   // shades of the quad vertices
   glBindBuffer(GL_ARRAY_BUFFER, shades_buffer_id);
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(float), 0);

   // ------------------------- triangles -----------------------------------------

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
   glDisableVertexAttribArray(1);
   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

   first_time = false;

}

void
HUDMesh::setup_instancing_buffer(unsigned int n_boxes, unsigned int size_of_bar) {

   // make *space* for the instancing values, but don't fill them with data (here)

   max_n_instances = n_boxes; // as much space as we have allocated.
   n_instances = 0; // this is how many we want to draw now.

   glBindVertexArray(vao);

   glGenBuffers(1, &inst_hud_bar_attribs_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_hud_bar_attribs_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_boxes * size_of_bar, nullptr, GL_DYNAMIC_DRAW);

   // layout
   // 0 vec2 position        (quad)
   // 1 float shade          (quad)
   // 2 vec4 colour          (instanced)
   // 3 vec2 position_offset (instanced)
   // 4 float scale (width)  (instanced)
   // 5 float scale (height) (instanced)

   // colour - instanced
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, size_of_bar, 0);
   glVertexAttribDivisor(2, 1);

   // position_offset - instanced
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, size_of_bar, reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(3, 1);

   // scale_x (width of the bar - along the x axis) - instanced
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, size_of_bar, reinterpret_cast<void *>(sizeof(glm::vec4) + sizeof(glm::vec2)));
   glVertexAttribDivisor(4, 1);

   // scale_y (height of the bar) - instanced
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, size_of_bar,
                         reinterpret_cast<void *>(sizeof(glm::vec4) + sizeof(glm::vec2) + sizeof(float)));
   glVertexAttribDivisor(5, 1);

   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);

   glBindBuffer(GL_ARRAY_BUFFER, 0);

   // unbind the vao here I think.

}

void
HUDMesh::update_instancing_buffer_data(const std::vector<HUD_bar_attribs_t> &new_bars) {

   unsigned int s = new_bars.size();
   n_instances = s;
   if (s > max_n_instances)
         n_instances = max_n_instances; // only draw as many as we have allocated
   glBindBuffer(GL_ARRAY_BUFFER, inst_hud_bar_attribs_buffer_id);
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(HUD_bar_attribs_t), &(new_bars[0]));

}

void
HUDMesh::update_instancing_buffer_data(const std::vector<HUD_button_info_t> &buttons_info) {

   unsigned int s = buttons_info.size();
   n_instances = s;
   if (s> max_n_instances)
      n_instances = max_n_instances;
   glBindBuffer(GL_ARRAY_BUFFER, inst_hud_bar_attribs_buffer_id);
   // argh! - 20210828-PE OK! solved with inheritance.
   glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(HUD_button_info_t), &(buttons_info[0]));
}

// static
HUD_button_limits_t
HUD_button_info_t::get_button_limits(unsigned int button_index, int width, int height) {

   // matches below function (needed to test for hit: mouse-over/click)

   // 1.3 is used here as a multiplier to give the buttons some vertical room between them

   float left = 0.6;
   float right = left + button_width;
   float bottom = -0.9 + 1.3 * button_height * static_cast<float>(button_index);
   float top = bottom + HUD_button_info_t::button_height;
   HUD_button_limits_t lims(top, bottom, left, right);
   return lims;
   
}


#include <glm/gtx/string_cast.hpp> // to_string

// counting from the bottom! (at the moment)
//
// get the right position for the button
void
HUD_button_info_t::set_position_offset(unsigned int button_index, int width, int height) {

   // matches above function

   // 1.3 is used here as a multiplier to give the buttons some vertical room between them

   glm::vec2 po(0.6, -0.9 + 1.3 * button_height * static_cast<float>(button_index));
   std::cout << "set_position_offset(): button_index " << button_index << " " << glm::to_string(po) << std::endl;
   set_position_offset(po);
}


void
HUDMesh::draw(Shader *shader_p) {

   if (false)
      std::cout << "debug:: HUDMesh::draw() --- start --- n_instances: " << n_instances
                << " vao " << vao << std::endl;

   if (this_mesh_is_closed) return;

   if (n_instances == 0) return;

   shader_p->Use();

   glBindVertexArray(vao);

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glEnableVertexAttribArray(3);
   glEnableVertexAttribArray(4);
   glEnableVertexAttribArray(5);

   if (false)
      std::cout << "debug:: HUDMesh::draw() glDrawElementsInstanced()"
                << " of HUDMesh \"" << name << "\""
                << " with shader " << shader_p->name
                << " and " << n_instances << " instances" << std::endl;

   unsigned int n_triangle_vertices = triangles.size() * 3;

   glDrawElementsInstanced(GL_TRIANGLES, n_triangle_vertices, GL_UNSIGNED_INT, nullptr, n_instances);
   GLenum err = glGetError();
   if (err) std::cout << "error HUDMesh::draw() glDrawElementsInstanced()"
                      << " of HUDMesh \"" << name << "\""
                      << " with shader " << shader_p->name << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glUseProgram(0);

}
