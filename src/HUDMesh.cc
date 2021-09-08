
#include <iostream>
#include "HUDMesh.hh"

void
HUDMesh::init() {

   max_n_instances = 0;
   n_instances = 0;
   first_time = true;
   use_blending = false;
   inst_hud_bar_attribs_buffer_id = 0;
   scales_have_been_set = false;
   offset_position_has_been_set = false;
   scales = glm::vec2(1,1);
   offset_position = glm::vec2(0,0);
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
HUDMesh::setup_simple_camera_facing_quad() {

   vertices.clear();
   triangles.clear();

   vertices.push_back(glm::vec2(0.0f, 0.0f));
   vertices.push_back(glm::vec2(1.0f, 0.0f));
   vertices.push_back(glm::vec2(1.0f, 1.0f));
   vertices.push_back(glm::vec2(0.0f, 1.0f));

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
HUDMesh::setup_vertices_and_triangles_for_tooltip_background() {

   vertices.clear();
   triangles.clear();

   vertices.push_back(glm::vec2(0.0f, 0.0f));
   vertices.push_back(glm::vec2(1.0f, 0.0f));
   vertices.push_back(glm::vec2(1.0f, 1.0f));
   vertices.push_back(glm::vec2(0.0f, 1.0f));
   vertices.push_back(glm::vec2(0.0f, 1.3f)); // guess at geometry
   vertices.push_back(glm::vec2(0.2f, 1.0f));

   // now move the whole mesh down... good idea?
   for (auto &vertex : vertices)
      vertex += glm::vec2(0.0, -0.3);

   shades.push_back(-1);
   shades.push_back(-1);
   shades.push_back(1);
   shades.push_back(1);
   shades.push_back(1);

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));
   triangles.push_back(g_triangle(2,4,5));

   setup_buffers();
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

#include <glm/gtx/string_cast.hpp> // to_string


HUD_button_limits_t
HUD_button_info_t::get_button_limits(int width, int height) const {

   // matches below function (needed to test for hit: mouse-over/click)

   // 1.3 is used here as a multiplier to give the buttons some vertical room between them

   float ww = static_cast<float>(width);
   float wh = static_cast<float>(height);

   glm::vec2 po = calculate_position_offset(position_offset_index, width, height);
   float left   = po.x;
   float bottom = po.y;
   float adjusted_button_width_opengl_coords  = button_width  * static_cast<float>(900)/ww;
   float adjusted_button_height_opengl_coords = button_height * static_cast<float>(900)/wh;
   float right = left + adjusted_button_width_opengl_coords;
   float top  = bottom + adjusted_button_height_opengl_coords;

   HUD_button_limits_t lims(top, bottom, left, right); // opengl coords
   return lims;

}

// static
glm::vec2
HUD_button_info_t::calculate_position_offset(unsigned int button_index, int width, int height) {

   // OK so with a GL widget of 900x900, the buttons are (more or less) the right size and the right width
   // and the right height.  That is with a default button_width of 0.3 and a default button_height of 0.06
   //
   // Let's scale things to match
   //
   float ww = static_cast<float>(width);
   float wh = static_cast<float>(height);
   //
   // So, in GL coords, the left and right of the button is 0.6 and 0.9 - in pixels, that's 540 (ww-360) and 810 (ww-90)
   //                   the bottom of the bottom button is -0.9         - in pixels, that's 810 (www-90)
   //                   thb button height is 0.06                       - in pixels, that's  54

   const float standard_offset_right_margin_n_pixels  = 90.0;
   const float standard_offset_bottom_margin_n_pixels = 90.0;

   const float button_right_margin_opengl_coords  = standard_offset_right_margin_n_pixels  / ww;
   const float button_bottom_margin_opengl_coords = standard_offset_bottom_margin_n_pixels / wh;

   float adjusted_button_width_opengl_coords  = button_width  * static_cast<float>(900)/ww;
   float adjusted_button_height_opengl_coords = button_height * static_cast<float>(900)/wh;
   float button_left_opengl_coords = (1.0-button_right_margin_opengl_coords) - adjusted_button_width_opengl_coords;

   glm::vec2 po(button_left_opengl_coords,
                (-1+button_bottom_margin_opengl_coords) + 1.3 * adjusted_button_height_opengl_coords * static_cast<float>(button_index));

   if (false)
      std::cout << "debug:: button_right_margin_opengl_coords " << button_right_margin_opengl_coords
                << " adjusted_button_width_opengl_coords " << adjusted_button_width_opengl_coords << " "
                << "set_position_offset(): button_index " << button_index << " " << glm::to_string(po) << std::endl;

   return po;

}

void
HUD_button_info_t::set_scales_and_position_offset(unsigned int button_index, int glarea_width, int glarea_height) {

   glm::vec2 po = calculate_position_offset(button_index, glarea_width, glarea_height);

   set_position_offset(button_index, po);

   float ww = static_cast<float>(glarea_width);
   float wh = static_cast<float>(glarea_height);
   scale_x *= static_cast<float>(900)/ww;
   scale_y *= static_cast<float>(900)/wh;
}



void
HUDMesh::draw(Shader *shader_p) { // in this case draw() is draw_instanced() --- maybe rename it later
                                  // to be less confusing/more consistent.

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

   if (scales_have_been_set)
      shader_p->set_vec2_for_uniform("scales", scales);
   if (offset_position_has_been_set)
      shader_p->set_vec2_for_uniform("offset_position", offset_position);

   if (false) {
      if (scales_have_been_set)
         std::cout << "HUDMesh sending scales          " << glm::to_string(scales) << std::endl;
      if (offset_position_has_been_set)
         std::cout << "HUDMesh sending offset_position " << glm::to_string(offset_position) << std::endl;
   }

   unsigned int n_triangle_vertices = triangles.size() * 3;

   if (false)
      std::cout << "debug:: HUDMesh::draw() glDrawElementsInstanced()"
                << " of HUDMesh \"" << name << "\""
                << " with shader " << shader_p->name
                << " with " << n_triangle_vertices << " triangle verices "
                << " and " << n_instances << " instances" << std::endl;

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
