#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "globjects.h"
#include "graphics-info.h"

#include "draw.hh"
#include "draw-2.hh"
#include "framebuffer.hh"

#include "text-rendering-utils.hh"
#include "cc-interface-scripting.hh"

// header
glm::mat4 get_view_rotation();
glm::vec3 get_eye_position();

enum {VIEW_CENTRAL_CUBE, ORIGIN_CUBE};

gint idle_contour_function(gpointer data);

// maybe this can go in the draw-2.hh header
glm::vec4 new_unproject(float z);


float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
      // positions   // texCoords
      -1.0f,  1.0f,  0.0f, 1.0f,
      -1.0f, -1.0f,  0.0f, 0.0f,
       1.0f, -1.0f,  1.0f, 0.0f,

      -1.0f,  1.0f,  0.0f, 1.0f,
       1.0f, -1.0f,  1.0f, 0.0f,
       1.0f,  1.0f,  1.0f, 1.0f
   };

void init_shaders() {

   graphics_info_t g;
   g.init_shaders();

}

void
graphics_info_t::init_screen_quads() {

   graphics_info_t::shader_for_screen.Use();
   // screen quad VAO
   unsigned int quadVBO;
   glGenVertexArrays(1, &graphics_info_t::screen_quad_vertex_array_id);
   glBindVertexArray(graphics_info_t::screen_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   GLenum err = glGetError();
   if (true) std::cout << "init_screen_quads() err is " << err << std::endl;

}

void
graphics_info_t::init_blur_quads() {

   graphics_info_t::shader_for_blur.Use();
   // screen quad VAO
   unsigned int quadVBO;
   glGenVertexArrays(1, &graphics_info_t::blur_quad_vertex_array_id);
   glBindVertexArray(graphics_info_t::blur_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   GLenum err = glGetError();
   if (true) std::cout << "init_blur_quads() err is " << err << std::endl;

}

// void graphics_info_t::init_central_cube();

void
graphics_info_t::init_buffers() {
   init_central_cube();
   init_screen_quads();
   init_blur_quads();
}

void
graphics_info_t::init_central_cube() {

   float positions[24] = {
                          -0.5,  -0.5, -0.5,
                          -0.5,  -0.5,  0.5,
                          -0.5,   0.5, -0.5,
                          -0.5,   0.5,  0.5,
                           0.5,  -0.5, -0.5,
                           0.5,  -0.5,  0.5,
                           0.5,   0.5, -0.5,
                           0.5,   0.5,  0.5
   };

   glUseProgram(graphics_info_t::shader_for_central_cube.get_program_id());
   GLenum err = glGetError();
   if (err) std::cout << "init_central_cube() glUseProgram() err is " << err << std::endl;

   // number of lines * 2:
   unsigned int indices[24] { 0,1, 1,5, 5,4, 4,0, 2,3, 3,7, 7,6, 6,2, 0,2, 1,3, 5,7, 4,6 };

   // GLuint VertexArrayID;
   glGenVertexArrays(1, &graphics_info_t::central_cube_vertexarray_id);
   glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);

   // GLuint vertexbuffer;
   glGenBuffers(1, &graphics_info_t::central_cube_array_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::central_cube_array_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 24, &positions[0], GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

   // unsigned int ibo;
   glGenBuffers(1, &graphics_info_t::central_cube_index_buffer_id);
   err = glGetError();
   if (err) std::cout << "init_central_cube() index glGenBuffers() err is " << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::central_cube_index_buffer_id);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 24, &indices[0], GL_STATIC_DRAW);
   err = glGetError();
   if (err) std::cout << "init_central_cube() glBufferData() err is " << err << std::endl;
   glBindVertexArray(0);

}

void
graphics_info_t::init_hud_text() {

   std::cout << "------------------ init_hud_text() ---------------------\n";
   graphics_info_t g;
   g.load_freetype_font_textures();
   glUseProgram(g.shader_for_hud_text.get_program_id());
   GLenum err = glGetError();
   if (err) std::cout << "init_hud_text() glUseProgram() err is " << err << std::endl;
   glGenVertexArrays(1, &graphics_info_t::hud_text_vertexarray_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glGenVertexArrays() err is " << err << std::endl;
   glBindVertexArray(graphics_info_t::hud_text_vertexarray_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glBindVertexArray() err is " << err << std::endl;
   glGenBuffers(1, &graphics_info_t::hud_text_array_buffer_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glGenBuffers() err is " << err << std::endl;
   glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::hud_text_array_buffer_id);
   err = glGetError(); if (err) std::cout << "init_hud_text() glBindBuffer() err is " << err << std::endl;
   glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
   err = glGetError(); if (err) std::cout << "init_hud_text() glBufferData() err is " << err << std::endl;
   glEnableVertexAttribArray(0);
   err = glGetError(); if (err) std::cout << "init_hud_text() glEnableVertexAttribArray() err is " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glBindVertexArray(0);
   std::cout << "------------------ done init_hud_text() ---------------------\n";
}

// static
glm::mat4
graphics_info_t::get_model_view_matrix() {

   // where is this used?

   if (! perspective_projection_flag) {
      return glm::toMat4(glm_quat);
   } else {
      return glm::mat4(1.0);
   }
}

// static
glm::mat4
graphics_info_t::get_molecule_mvp() {

   // presumes that we are in the correct programID

   float w = static_cast<float>(graphics_info_t::graphics_x_size);
   float h = static_cast<float>(graphics_info_t::graphics_y_size);
   float screen_ratio = w/h;

   // I don't think that the quaternion belongs to the model matrix, it should be
   // part of the view matrix I think.
   // Yes. That's right.
   glm::mat4 model_matrix = glm::mat4(1.0);

   float z = graphics_info_t::zoom * 0.04;
   glm::vec3 sc(z,z,z);

   GLfloat near = -0.1 * zoom * clipping_front;
   GLfloat far  =  0.3 * zoom * clipping_back;

   if (false)
      std::cout << "near " << near << " far " << far << " clipping front "
                << clipping_front << " back " << clipping_back << std::endl;

   float sr = screen_ratio;
   glm::mat4 projection_matrix = glm::ortho(-0.3f*zoom*sr, 0.3f*zoom*sr,
                                            -0.3f*zoom,    0.3f*zoom,
                                            near, far);
   

   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   // std::cout << "rotation centre " << glm::to_string(rc) << std::endl;
   glm::mat4 view_matrix = glm::toMat4(graphics_info_t::glm_quat);

   view_matrix = glm::translate(view_matrix, -rc);
   // view_matrix = glm::scale(view_matrix, reverse_z); causes weirdness - not sure about handedness
   glm::mat4 mvp = projection_matrix * view_matrix * model_matrix;

   if (false) {
      std::cout << "get_molecule_mvp: " << glm::to_string(projection_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(view_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(model_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(mvp) << std::endl;
   }

   if (graphics_info_t::perspective_projection_flag) {

      // for fun/testing
      // turn off view scaling when tinkering with this?
      // there should not be a concept of "zoom" with perspective view, just translation
      // along screen-Z.

      float fov = 40.0;

      glm::vec4 up_1(0,1,0,1);
      glm::mat4 trackball_matrix = glm::toMat4(graphics_info_t::glm_quat);
      glm::vec4 up_2 = trackball_matrix * up_1;
      glm::vec3 up = glm::vec3(up_2);
      
      glm::vec3 ep = eye_position;

      if (false)
         std::cout << "eye position " << glm::to_string(ep) << " rc " << glm::to_string(rc)
                   << " up " << glm::to_string(up) << std::endl;

      view_matrix = glm::lookAt(ep, rc, up);
      fov = 30.0; // degrees

      glm::mat4 projection_matrix_persp = glm::perspective(glm::radians(fov),
                                                           screen_ratio,
                                                           screen_z_near_perspective,
                                                           screen_z_far_perspective);
      mvp = projection_matrix_persp * view_matrix * model_matrix;
   }


   return mvp;
}

// can we work out the eye position without needing to unproject? (because that depends
// on get_molecule_mvp()...
//
glm::vec3
graphics_info_t::get_eye_position() {

   if (! graphics_info_t::perspective_projection_flag) {

      // orthograph eye position is inferred from centre psotion
      // and zoom and view rotation.

      glm::vec3 test_vector_1(0.0, 0.0, 1.0);
      glm::vec3 test_vector_2(1.0, 1.0, 0.0);

      glm::mat4 vr = get_view_rotation();
      glm::vec4 rot_test_vector_1 = glm::vec4(test_vector_1, 1.0) * vr;
      glm::vec4 rot_test_vector_2 = glm::vec4(test_vector_2, 1.0) * vr;

      glm::vec3 ep = graphics_info_t::zoom * glm::vec3(rot_test_vector_1);
      glm::vec3 rc = graphics_info_t::get_rotation_centre();
      ep += rc;

      return ep;
   } else {

      // Perspective projection we manipulate the eye position directly

      return graphics_info_t::eye_position;

   }

}

glm::vec4
graphics_info_t::new_unproject(float z) {

   std::cout << "This is not used (is it?) A " << std::endl;
   return glm::vec4(0,1,0,1);
#if 0
   // z is 1 and -1 for front and back (or vice verse).
   GtkAllocation allocation;
   gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
   float w = allocation.width;
   float h = allocation.height;
   graphics_info_t g;
   float mouseX = 2.0 *    g.GetMouseBeginX()/w  - 1.0f;
   float mouseY = 2.0 * (h-g.GetMouseBeginY())/h - 1.0f;
   std::cout << "debug in new_unproject widget w and h " << w << " " << h << std::endl;
   std::cout << "debug in new_unproject mouse x and y widget  " << g.GetMouseBeginX() << " " << g.GetMouseBeginY() << std::endl;
   std::cout << "debug in new_unproject mouse x and y GL      " << mouseX << " " << mouseY << std::endl;
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   float real_y = - mouseY; // in range -1 -> 1
   glm::vec4 screenPos_f = glm::vec4(mouseX, real_y, z, 1.0f);
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   // std::cout << "screen_pos " << glm::to_string(screenPos_f) << std::endl;
   // std::cout << "world_pos " << glm::to_string(worldPos_f) << std::endl;
   return worldPos_f;
#endif
}


glm::vec4
graphics_info_t::new_unproject(float x, float y, float z) {

   std::cout << "This is not used (is it?) B" << std::endl;
   return glm::vec4(0,1,0,1);

#if 0
   // z is 1 and -1 for front and back (or vice verse).
   // GtkAllocation allocation;
   // gtk_widget_get_allocation(graphics_info_t::glarea, &allocation);
   // int w = allocation.width;
   // int h = allocation.height;
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos_f = glm::vec4(x, y, z, 1.0f); // maybe +1
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   if (false) {
      std::cout << "   " << glm::to_string(mvp) << std::endl;
      std::cout << "   " << glm::to_string(vp_inv) << std::endl;
      std::cout << "   " << glm::to_string(screenPos_f) << std::endl;
      std::cout << "   " << glm::to_string(worldPos_f) << std::endl;
   }
   return worldPos_f;
#endif
}

glm::mat4
graphics_info_t::get_view_rotation() {

   // need to be in the correct program (well, the model-drawing part)

   if (! perspective_projection_flag)
      return glm::toMat4(graphics_info_t::glm_quat);
   else
      return glm::mat4(1.0);
}


void
graphics_info_t::setup_map_uniforms(const Shader &shader,
                                    const glm::mat4 &mvp,
                                    const glm::mat4 &view_rotation,
                                    float density_surface_opacity) {

   glUniformMatrix4fv(shader.mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   GLenum err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniformMatrix4fv() mvp " << err << std::endl;
   glUniformMatrix4fv(shader.view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniformMatrix4fv() vr  " << err << std::endl;

   GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
   glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
   glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
   err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniform4fv() for bg  " << err << std::endl;

   // opacity: (I can't get this to work for lines)
   GLuint opacity_uniform_location = shader.map_opacity_uniform_location;
   float opacity = density_surface_opacity;
   glUniform1f(opacity_uniform_location, opacity);
   err = glGetError(); if (err) std::cout << "   setup_map_uniforms() glUniformf() for opacity "
                                          << err << std::endl;

   GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
   glm::vec4 ep = glm::vec4(get_eye_position(), 1.0);
   glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
   err = glGetError(); if (err) std::cout << "   setup_map_uniforms() glUniform4fv() for eye position "
                                          << err << std::endl;

   // lights
   std::map<unsigned int, gl_lights_info_t>::const_iterator it;
   it = graphics_info_t::lights.find(0);
   if (it != graphics_info_t::lights.end()) {
      const gl_lights_info_t &light = it->second;
      glUniform1i(shader.light_0_is_on_uniform_location, light.is_on);
      glUniform4fv(shader.light_0_position_uniform_location, 1, glm::value_ptr(light.position));
   }
   it = graphics_info_t::lights.find(1);
   if (it != graphics_info_t::lights.end()) {
      const gl_lights_info_t &light = it->second;
      glUniform1i(shader.light_1_is_on_uniform_location, light.is_on);
      glUniform4fv(shader.light_1_position_uniform_location, 1, glm::value_ptr(light.position));
   }
}


void
graphics_info_t::draw_map_molecules(bool draw_transparent_maps) {

   // run throgh this molecule loop twice - for opaque then transparent maps
   // first, a block that decides if we need to do anything.

   bool needs_blend_reset = false;

   unsigned int n_transparent_maps = 0;
   if (draw_transparent_maps) {
      for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         const molecule_class_info_t &m = graphics_info_t::molecules[ii];
         if (! m.draw_it_for_map) continue;
         if (! m.is_an_opaque_map())
            n_transparent_maps++;
      }
      if (n_transparent_maps > 0) {
         needs_blend_reset = true;
         glEnable(GL_BLEND);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }
   }

   if (!draw_transparent_maps || n_transparent_maps > 0) {

      glLineWidth(map_line_width);
      GLenum err = glGetError();
      if (err) std::cout << "gtk3_draw_map_molecules() glLineWidth " << err << std::endl;

      GLuint pid = graphics_info_t::shader_for_maps.get_program_id();
      glUseProgram(pid);
      err = glGetError();
      if (err) std::cout << "gtk3_draw_map_molecules() glUseProgram with GL err "
                         << err << std::endl;

      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation = get_view_rotation();
      if (perspective_projection_flag)
         view_rotation = glm::mat4(1.0);

      glEnable(GL_DEPTH_TEST); // this needs to be in the draw loop!?
      glDepthFunc(GL_LESS);

      Shader &shader = graphics_info_t::shader_for_maps;

      glm::vec4 ep(get_eye_position(), 1.0);

      for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
         const molecule_class_info_t &m = graphics_info_t::molecules[ii];
         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         if (! m.draw_it_for_map) continue;
         if (draw_transparent_maps)
            if (m.is_an_opaque_map())
               continue; // not this round

         if (m.n_vertices_for_map_VertexArray > 0) {

            bool draw_with_lines = true;
            if (!m.draw_it_for_map_standard_lines) draw_with_lines = false;

            glUniform1i(shader.is_perspective_projection_uniform_location,
                        graphics_info_t::perspective_projection_flag);
            err = glGetError(); if (err) std::cout << "   draw_map_molecules() error B " << std::endl;

            if (draw_with_lines) {
               // I don't see why this is needed - but it is.
               if (! m.is_an_opaque_map())
                  glEnable(GL_BLEND);

               glBindVertexArray(m.m_VertexArrayID_for_map);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glBindVertexArray() "
                                  << graphics_info_t::molecules[ii].m_VertexArrayID_for_map
                                  << " with GL err " << err << std::endl;

               glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_map_lines_ID);

               graphics_info_t::setup_map_uniforms(shader, mvp, view_rotation, m.density_surface_opacity);
               glDrawElements(GL_LINES, m.n_vertices_for_map_VertexArray,
                              GL_UNSIGNED_INT, nullptr);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glDrawElements() n_vertices: "
                                  << m.n_vertices_for_map_VertexArray
                                  << " with GL err " << err << std::endl;
            }

            if (!draw_with_lines) { // draw as a solid object
               if (false)
                  std::cout << "   draw_map_molecules(): imol " << ii
                            << " array_id and n_vertices_for_VertexArray: "
                            << m.m_VertexArrayID_for_map << " "
                            << m.n_indices_for_triangles
                            << std::endl;

               if (! m.is_an_opaque_map()) {
                  // sort the triangles
                  clipper::Coord_orth eye_pos_co(ep.x, ep.y, ep.z);
                  graphics_info_t::molecules[ii].sort_map_triangles(eye_pos_co);
               }

               glBindVertexArray(m.m_VertexArrayID_for_map);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glBindVertexArray() "
                                  << m.m_VertexArrayID_for_map << " with GL err " << err << std::endl;
               glEnable(GL_BLEND);
               glBindBuffer(GL_ARRAY_BUFFER,         m.m_VertexBufferID);
               glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_map_triangles_ID);

               glUniformMatrix4fv(graphics_info_t::shader_for_maps.mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() " << err << std::endl;
               glUniformMatrix4fv(shader.view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() " << err << std::endl;

               GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
               glm::vec4 bgc(background_colour, 1.0);
               glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glUniform4fv() for bg  " << err << std::endl;

               // opacity:
               GLuint opacity_uniform_location = graphics_info_t::shader_for_maps.map_opacity_uniform_location;
               float opacity = m.density_surface_opacity;
               glUniform1f(opacity_uniform_location, opacity);
               err = glGetError(); if (err) std::cout << "   draw_map_molecules() glUniformf() for opacity "
                                                      << err << std::endl;

               GLuint eye_position_uniform_location = graphics_info_t::shader_for_maps.eye_position_uniform_location;
               glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));

               glDrawElements(GL_TRIANGLES, m.n_indices_for_triangles, GL_UNSIGNED_INT, nullptr);

               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glDrawElements() n_indices_for_triangles "
                                  << graphics_info_t::molecules[ii].n_indices_for_triangles
                                  << " with GL err " << err << std::endl;
            }
         }
      }
   }

   // to be clean we should use
   // glDisableVertexAttribArray(0);
   // here.
   // that would mean adding glEnableVertexAttribArray() for the attributes (position, normal, colour).
   // in the above block.

   if (needs_blend_reset) {
      glDisable(GL_BLEND);
   }
}

void
graphics_info_t::draw_model_molecules() {

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();

   Shader &shader = graphics_info_t::shader_for_models;
   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {

      const molecule_class_info_t &m = graphics_info_t::molecules[ii];
      if (! graphics_info_t::is_valid_model_molecule(ii)) continue;
      if (! m.draw_it) continue;

      if (false)
         std::cout << "imol " << ii << " n_vertices_for_model_VertexArray "
                   << m.n_vertices_for_model_VertexArray << std::endl;
      if (m.n_vertices_for_model_VertexArray > 0) {

         glDisable(GL_BLEND); // stop semi-transparent bonds - but why do we have them?
         gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));

         shader.Use();
         GLuint err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glUseProgram() "
                                                       << err << std::endl;

         glUniform1i(shader.is_perspective_projection_uniform_location,
                     graphics_info_t::perspective_projection_flag);

         glBindVertexArray(m.m_VertexArray_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glBindVertexArray() "
                            << m.m_VertexArray_for_model_ID
                            << " with GL err " << err << std::endl;

         glBindBuffer(GL_ARRAY_BUFFER, m.m_VertexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() v " << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() i " << err << std::endl;

         GLuint mvp_location           = shader.mvp_uniform_location;
         GLuint view_rotation_location = shader.view_rotation_uniform_location;

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() pre mvp " << err << std::endl;
         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for mvp " << err << std::endl;
         glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for view_rotation " << err << std::endl;
         // std::cout << glm::to_string(mvp) << std::endl;
         // std::cout << glm::to_string(view_rotation) << std::endl;

         GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
         glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec4 ep(get_eye_position(), 1.0);
         glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for eye position " << err << std::endl;

         // draw with the vertex count, not the index count.
         GLuint n_verts = graphics_info_t::molecules[ii].n_indices_for_model_triangles;
         // std::cout << "   Drawing " << n_verts << " model vertices" << std::endl;
         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glDrawElements() "
                            << n_verts << " with GL err " << err << std::endl;

         draw_molecule_atom_labels(m, mvp, view_rotation);

      }
   }
}

void
graphics_info_t::draw_molecule_atom_labels(const molecule_class_info_t &m,
                                           const glm::mat4 &mvp,
                                           const glm::mat4 &view_rotation) {

   // pass the glarea widget width and height.

   // put a triangle or square where the atom label should be, facing the camera
   // "billboarding"

   glm::vec3 label_colour(font_colour.red, font_colour.green, font_colour.blue);

   if (false) { // test label
      // Put atom label test at 42, 9, 13
      glm::vec3 point(42, 9, 13);
      // point = glm::vec3(0,0,0);

      glm::vec4 projected_point_2 = mvp * glm::vec4(point, 1.0);
      if (true)
         std::cout << "projected point " << glm::to_string(projected_point_2) << std::endl;

      projected_point_2.x = 0.5 * (projected_point_2.x + 1.0f);
      projected_point_2.y = 0.5 * (projected_point_2.y + 1.0f);

      projected_point_2.x *= 900.0;
      projected_point_2.y *= 900.0;

      glm::vec3 pp(projected_point_2);

      glEnable(GL_DEPTH_TEST); // or we don't see the label. Either a blurred label or nothing.
      render_atom_label(shader_for_atom_labels, ". Test Label", pp, 1.0, label_colour);
   }

   int n_atoms_to_label = m.labelled_atom_index_list.size();
   if (n_atoms_to_label == 0) return;

   // maybe pass these?
   GtkAllocation allocation;
   GtkWidget *widget = graphics_info_t::glareas[0];
   if (! widget) return;
   gtk_widget_get_allocation(widget, &allocation);

   for (int ii=0; ii<n_atoms_to_label ; ii++) {
      std::pair<std::string, clipper::Coord_orth> lab_pos =
         m.make_atom_label_string(ii, brief_atom_labels_flag, seg_ids_in_atom_labels_flag);
      const clipper::Coord_orth &co = lab_pos.second;
      const std::string &label = lab_pos.first;
      if (false)
         std::cout << "Atom label at " << co.format() << " "
                   << coot::util::single_quote(lab_pos.first) << std::endl;

      glm::vec3 point(co.x(), co.y(), co.z());
      glm::vec4 projected_point = mvp * glm::vec4(point, 1.0);
      // convert axes from -1 -> 1 to 0 -> 1
      projected_point.x = 0.5 * (projected_point.x + 1.0f);
      projected_point.y = 0.5 * (projected_point.y + 1.0f);
      projected_point.x *= allocation.width;
      projected_point.y *= allocation.height;
      glm::vec3 pp(projected_point);
      render_atom_label(shader_for_atom_labels, label, pp, 1.0, label_colour);
   }

}

void
graphics_info_t::draw_intermediate_atoms() { // draw_moving_atoms()

   // all these draw functions should be moved int graphics_info_t.

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();

   Shader &shader = graphics_info_t::shader_for_models;
   if (true) {
      const molecule_class_info_t &m = graphics_info_t::moving_atoms_molecule;

      if (false)
         std::cout << " moving atoms molecule: n_vertices_for_model_VertexArray "
                   << m.n_vertices_for_model_VertexArray << std::endl;

      if (m.n_vertices_for_model_VertexArray > 0) {

         glDisable(GL_BLEND); // stop semi-transparent bonds - but why do we have them?
         gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));


         shader.Use();
         GLuint err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glUseProgram() "
                                                       << err << std::endl;

         glBindVertexArray(m.m_VertexArray_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glBindVertexArray() "
                            << m.m_VertexArray_for_model_ID
                            << " with GL err " << err << std::endl;

         // should not be needed?
         glBindBuffer(GL_ARRAY_BUFFER, m.m_VertexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() v " << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() i " << err << std::endl;

         GLuint mvp_location           = shader.mvp_uniform_location;
         GLuint view_rotation_location = shader.view_rotation_uniform_location;

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() pre mvp " << err << std::endl;
         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for mvp " << err << std::endl;
         glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for view_rotation " << err << std::endl;
         // std::cout << glm::to_string(mvp) << std::endl;
         // std::cout << glm::to_string(view_rotation) << std::endl;

         GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
         glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec4 ep = glm::vec4(get_eye_position(), 1.0);
         glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for eye position " << err << std::endl;

         // draw with the vertex count, not the index count.
         GLuint n_verts = m.n_indices_for_model_triangles;
         // std::cout << "   Drawing " << n_verts << " model vertices" << std::endl;
         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glDrawElements() "
                            << n_verts << " with GL err " << err << std::endl;

      }
   }
}

void
graphics_info_t::draw_molecular_triangles() {
#ifdef USE_MOLECULES_TO_TRIANGLES
   // Martin's triangular molecules
   //
   // centre of the screen
   FCXXCoord pos(graphics_info_t::RotationCentre_x(),
   graphics_info_t::RotationCentre_y(),
   graphics_info_t::RotationCentre_z());
   // where is the eye?  That's what we want.
   // front plane is at z=0;
   GtkAllocation allocation;
   GtkWidget *widget = graphics_info_t::glareas[0];
   if (! widget) return;
   gtk_widget_get_allocation(widget, &allocation);
   coot::Cartesian tp_1_cart = unproject_xyz(allocation.width/2,
                                             allocation.height/2, 1);
   FCXXCoord tp_1(tp_1_cart.x(), tp_1_cart.y(), tp_1_cart.z());
   FCXXCoord diff = tp_1 - pos;
   FCXXCoord eye_pos = pos + diff * 5.0;
   // std::cout << "eye_pos: " << eye_pos << "\n";
   // coot::Cartesian eye_cart = pos + 20 * diff;
   // FCXXCoord eye_pos(eye_cart.x(), eye_cart.y(), eye_cart.z());
   if (graphics_info_t::mol_tri_scene_setup) {
      if (graphics_info_t::mol_tri_renderer) {
         //Can retrieve reference to the light if so preferred
         // This doesn't move the lights
         // FCXXCoord random_trans(50.0 * coot::util::random()/float(RAND_MAX),
         // 		           50.0 * coot::util::random()/float(RAND_MAX),
         //                        50.0 * coot::util::random()/float(RAND_MAX));
         FCXXCoord light_pos = pos + diff * 10; //  + random_trans;
         FCXXCoord neg_light_pos = pos + diff * 10; // - random_trans;

         graphics_info_t::mol_tri_scene_setup->getLight(0)->setTranslation(light_pos);
         graphics_info_t::mol_tri_scene_setup->getLight(1)->setTranslation(neg_light_pos);

         for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
            if (graphics_info_t::is_valid_model_molecule(ii)) {
               if (graphics_info_t::molecules[ii].draw_it) {
                  if (graphics_info_t::molecules[ii].molrepinsts.size()) {
                     std::cout << "----------------------- in draw_molecular_triangles() calling Martin code now... \n";
                     // molrepinsts get added to mol_tri_scene_setup when then are made
                     GLenum err = glGetError();
                     if (err) std::cout << "gl error pre-renderer in draw_molecular_triangles() " << err << std::endl;
                     // turns on glLighting.
                     graphics_info_t::mol_tri_scene_setup->renderWithRendererFromViewpoint(graphics_info_t::mol_tri_renderer,
                                                                                           eye_pos);
                     err = glGetError();
                     if (err) std::cout << "gl error in draw_molecular_triangles() " << err << std::endl;
                  }
               }
            }
         }
      }
   }
#endif
}

void
graphics_info_t::draw_molecules() {

   // opaque things

   draw_model_molecules();
   draw_intermediate_atoms();
   // draw_molecular_triangles(); // Martin's renderings

   // transparent things... (maybe this function (draw_molecules()) should not be split out).

   draw_map_molecules(false); // transparency
   draw_map_molecules(true);
}

void
graphics_info_t::draw_cube(GtkGLArea *glarea, unsigned int cube_type) {

   gtk_gl_area_make_current(glarea);
   glLineWidth(2.0);  // GLv4 antialiasing - OpenGL implementations are not required to support this
   GLenum err = glGetError();
   if (err) std::cout << "   error draw_central_cube() A err " << err << std::endl;

   // To see the possible values of the line width in aliased mode:
   // GLfloat line_width_max_min[2] = {0.0f, 0.0f};
   // glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, lineWidthRange);
   // This may not be possible in GL_LINE_SMOOTH mode.

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation(); // hhmm... naming

   glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);
   err = glGetError(); if (err) std::cout << "   error draw_central_cube() B err " << err << std::endl;
   glUseProgram(graphics_info_t::shader_for_central_cube.get_program_id());
   err = glGetError(); if (err) std::cout << "   error draw_central_cube() C err " << err << std::endl;
   glm::mat4 view_orientation = glm::toMat4(graphics_info_t::glm_quat);
   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   if (cube_type == VIEW_CENTRAL_CUBE) {
      mvp = glm::translate(mvp, rc);
      glm::vec3 sc(0.2f, 0.2f, 0.2f);
      mvp = glm::scale(mvp, sc);
   }
   if (cube_type == ORIGIN_CUBE) {
      glm::vec3 sc(0.3f, 0.3f, 0.3f);
      mvp = glm::scale(mvp, sc);
   }

   // we don't diverge here on the cube tye. Maybe change the name of the shader
   // because it does both
   Shader &shader = graphics_info_t::shader_for_central_cube;

   // we do this for all the shaders - Hmm.
   {
      GLuint mvp_location           = shader.mvp_uniform_location;
      GLuint view_rotation_location = shader.view_rotation_uniform_location;

      glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
      err = glGetError();
      if (err) std::cout << "   error draw_central_cube() glUniformMatrix4fv() for mvp " << err << std::endl;
      glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
      err = glGetError();
      if (err) std::cout << "   error draw_central_cube() glUniformMatrix4fv() for view_rotation " << err << std::endl;

      GLuint line_colour_uniform_location = shader.line_colour_uniform_location;
      glm::vec4 lc(0.5, 0.4, 0.4, 1.0);
      if (cube_type == ORIGIN_CUBE)
         lc = glm::vec4(0.6, 0.6, 0.4, 1.0);
      glUniform4fv(line_colour_uniform_location, 1, glm::value_ptr(lc));
      err = glGetError();
      if (err) std::cout << "   error draw_central_cube() glUniform4fv() for line colour " << err << std::endl;

      GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
      glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
      glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
      err = glGetError();
      if (err) std::cout << "   error draw_central_cube() glUniform4fv() for background " << err << std::endl;

   }

   glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "draw_central_cube() F glDrawElements() err " << err << std::endl;

   glBindVertexArray(0); // unbind
   glUseProgram(0);
}


void
graphics_info_t::draw_central_cube(GtkGLArea *glarea) {
   draw_cube(glarea, VIEW_CENTRAL_CUBE);
}

void
graphics_info_t::draw_origin_cube(GtkGLArea *glarea) {
   draw_cube(glarea, ORIGIN_CUBE);
}

GtkWidget *my_gtkglarea(GtkWidget *vbox) {

   GtkWidget *w = gtk_gl_area_new();
   gtk_widget_set_size_request(w, 900, 900);
   gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 2);
   return w;
}

void
graphics_info_t::setup_lights() {

   // not your old style lights

   gl_lights_info_t light;
   light.position = glm::vec4(-2.0f, -1.0f, 5.0f, 1.0f);
   graphics_info_t::lights[0] = light;
   light.position = glm::vec4( 3.0f, -2.0f, 4.0f, 1.0f);
   graphics_info_t::lights[1] = light;
}

void
on_glarea_realize(GtkGLArea *glarea) {

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   gtk_gl_area_make_current(glarea);
   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);
   GLenum err = glGetError();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() A err " << err << std::endl;

   // GLX_SAMPLE_BUFFERS_ARB
   // https://www.khronos.org/registry/OpenGL/extensions/ARB/ARB_multisample.txt
   // gdk/x11/gdkglcontext-x11.c
   // glXGetConfig(dpy, &visual_list[0], GLX_SAMPLE_BUFFERS_ARB, &gl_info[i].num_multisample);


   // glEnable(GL_MULTISAMPLE); // seems not to work at the moment. Needs work on the GTK->OpenGL interface 
   err = glGetError();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() A1 err " << err << std::endl;

   graphics_info_t g;
   g.init_shaders();
   g.init_buffers();
   err = glGetError();
   std::cout << "on_glarea_realize() post init_shaders() err is " << err << std::endl;

   graphics_info_t::shader_for_screen.Use(); // needed?

   err = glGetError(); std::cout << "start on_glarea_realize() err is " << err << std::endl;

   unsigned int index_offset = 0;
   graphics_info_t::screen_framebuffer.init(w,h, index_offset, "screen/occlusion");
   err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                          << err << std::endl;
   index_offset = 1;
   graphics_info_t::blur_framebuffer.init(w,h, index_offset, "blur");
   err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                          << err << std::endl;

   setup_hud_text(w, h, graphics_info_t::shader_for_hud_text, false);
   setup_hud_text(w, h, graphics_info_t::shader_for_atom_labels, true);

   graphics_info_t::shader_for_screen.Use();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() B screen framebuffer err " << err << std::endl;
   graphics_info_t::shader_for_screen.set_int_for_uniform("screenTexture", 0);
   err = glGetError(); if (err) std::cout << "on_glarea_realize() C screen framebuffer err " << err << std::endl;
   graphics_info_t::shader_for_screen.set_int_for_uniform("screenDepth", 1);
   err = glGetError(); if (err) std::cout << "on_glarea_realize() D screen framebuffer err " << err << std::endl;

   graphics_info_t::shader_for_blur.Use();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() blur shader-framebuffer B err " << err << std::endl;
   graphics_info_t::shader_for_blur.set_int_for_uniform("screenTexture", 0);
   err = glGetError(); if (err) std::cout << "on_glarea_realize() blur C shader-framebuffer err " << err << std::endl;
   graphics_info_t::shader_for_blur.set_int_for_uniform("screenDepth", 1);
   err = glGetError(); if (err) std::cout << "on_glarea_realize() blur D shader-framebuffer err " << err << std::endl;

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

   glEnable(GL_DEPTH_TEST);
   // glEnable(GL_CULL_FACE); // if I enable this, then I get to see th back side
                              // of the bonds (and atoms, possibly) It's a weird look

   // glEnable(GL_BLEND);
   // glEnable(GL_LINE_SMOOTH);

   // Make antialised lines
   if (false) {
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_LINE_SMOOTH);
   }

   g.setup_lights();

   GLfloat light0pos[4];
   glGetLightfv(GL_LIGHT0, GL_POSITION, light0pos);
   err = glGetError(); if (err) std::cout << "realsize() " << err << std::endl;
   std::cout << "... light 0: " << light0pos[0] << " " << light0pos[1] << " " << light0pos[2] << " " << std::endl;

   // Martin's Molecular triangles
   // setup_molecular_triangles();

#if !defined(USE_GUILE) && !defined(USE_PYTHON)
   // handle_command_line_data(cld);
#endif

   err = glGetError(); if (true) std::cout << "on_glarea_realize() --end-- with err " << err << std::endl;

   g.setup_key_bindings();
}

gboolean
on_glarea_render(GtkGLArea *glarea) {

   return graphics_info_t::render(glarea);

}

gboolean
graphics_info_t::render(GtkGLArea *glarea) {

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   GLenum err = glGetError();
   if (err) std::cout << "render() start " << err << std::endl;

   // is this needed?
   gtk_gl_area_make_current(glarea);
   err = glGetError(); if (err) std::cout << "render() post gtk_gl_area_make_current() " << err << std::endl;

   screen_framebuffer.bind();
   err = glGetError(); if (err) std::cout << "render() post screen_buffer bind() " << err << std::endl;

   glEnable(GL_DEPTH_TEST);

   {
      const glm::vec3 &bg = graphics_info_t::background_colour;
      glClearColor (bg[0], bg[1], bg[2], 1.0);
      err = glGetError(); if (err) std::cout << "render() B err " << err << std::endl;
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      err = glGetError(); if (err) std::cout << "render() C err " << err << std::endl;

      draw_central_cube(glarea);
      draw_origin_cube(glarea);
      err = glGetError(); if (err) std::cout << "render()  pre-draw-text err " << err << std::endl;

      draw_molecules();

      // Put atom label test at 42, 9, 13
      // draw_hud_text(w, h, shader_for_hud_text);
      // err = glGetError(); if (err) std::cout << "render() post-draw-text err " << err << std::endl;

      glBindVertexArray(0);
   }


   graphics_info_t::blur_framebuffer.bind();
   // glDisable(GL_DEPTH_TEST);
   glEnable(GL_DEPTH_TEST);

   // Screen shader (ambient occlusion)

   {
      graphics_info_t::shader_for_screen.Use();
      glBindVertexArray(graphics_info_t::screen_quad_vertex_array_id);

      // glClearColor(0.5, 0.2, 0.2, 1.0);
      const glm::vec3 &bg = graphics_info_t::background_colour;
      glClearColor (bg[0], bg[1], bg[2], 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      GLuint pid = graphics_info_t::shader_for_screen.get_program_id();
      glActiveTexture(GL_TEXTURE0 + 1);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::screen_framebuffer.get_texture_colour());
      glUniform1i(glGetUniformLocation(pid, "screenTexture"), 1);
      glActiveTexture(GL_TEXTURE0 + 2);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::screen_framebuffer.get_texture_depth());
      glUniform1i(glGetUniformLocation(pid, "screenDepth"), 2);
      err = glGetError(); if (err) std::cout << "on_glarea_render() D err " << err << std::endl;

      glDrawArrays(GL_TRIANGLES, 0, 6);
      err = glGetError(); if (err) std::cout << "on_glarea_render() E err " << err << std::endl;
   }

   // use this, rather than glBindFramebuffer(GL_FRAMEBUFFER, 0); ... just Gtk things.
   gtk_gl_area_attach_buffers(glarea);
   glEnable(GL_DEPTH_TEST);

   // z-blur shader
   {
      graphics_info_t::shader_for_blur.Use();
      glBindVertexArray(graphics_info_t::blur_quad_vertex_array_id);

      glClearColor(0.5, 0.2, 0.2, 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      Shader &shader = graphics_info_t::shader_for_blur;

      glUniform1f(shader.zoom_uniform_location, graphics_info_t::zoom);
      err = glGetError(); if (err) std::cout << "on_glarea_render() blur-A err " << err << std::endl;
      glUniform1i(shader.is_perspective_projection_uniform_location,
                  perspective_projection_flag);
      err = glGetError(); if (err) std::cout << "   on_glarea_render() blur-A2 error " << std::endl;

      GLuint pid = graphics_info_t::shader_for_blur.get_program_id();
      glActiveTexture(GL_TEXTURE0 + 1);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::blur_framebuffer.get_texture_colour());
      glUniform1i(glGetUniformLocation(pid, "screenTexture"), 1); // was 1
      glActiveTexture(GL_TEXTURE0 + 2);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::blur_framebuffer.get_texture_depth());
      glUniform1i(glGetUniformLocation(pid, "screenDepth"), 2); // was 2
      err = glGetError(); if (err) std::cout << "on_glarea_render() blur-B err " << err << std::endl;

      glDrawArrays(GL_TRIANGLES, 0, 6);
      err = glGetError(); if (err) std::cout << "on_glarea_render() blur-C err " << err << std::endl;
   }


   if (true) { //test to check if there is HUD text to draw, don't enter here if
               // it is not needed.

      // True HUD text (not atom labels) should be added *after* blurring (and here we are).
      // Currently, we don't see it.
      glDisable(GL_DEPTH_TEST); // needed because HUD text gets put at z=1.0 (it seems).
      draw_hud_text(w, h, shader_for_hud_text);
      err = glGetError(); if (err) std::cout << "render() post-draw-text err " << err << std::endl;
   }

   graphics_info_t::frame_counter++;
   if (graphics_info_t::frame_draw_queue.size() > 0) {
      std::chrono::time_point<std::chrono::system_clock> tp_now = std::chrono::high_resolution_clock::now();
      std::chrono::time_point<std::chrono::system_clock> queue_time = graphics_info_t::frame_draw_queue.front();
      graphics_info_t::frame_draw_queue.pop();
      auto delta_time = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - queue_time).count();
      if (false)
         std::cout << "INFO:: ---------- Timing check frame " << delta_time << " milliseconds" << " queue size "
                   << graphics_info_t::frame_draw_queue.size() << std::endl;
   }

   glFlush();
   return FALSE;
}

void
graphics_info_t::reset_frame_buffers(int width, int height) {
   graphics_info_t g;
   unsigned int index_offset = 0;
   g.screen_framebuffer.init(width, height, index_offset, "screen");
   index_offset = 1;
   g.blur_framebuffer.init(width, height, index_offset, "blur");
}

void
on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   graphics_info_t g;
   g.graphics_x_size = width;
   g.graphics_y_size = height;
   setup_hud_text(width, height, g.shader_for_hud_text, false);
   setup_hud_text(width, height, g.shader_for_atom_labels, true); // change the function name
   g.reset_frame_buffers(width, height);
}

gboolean
on_glarea_scroll(GtkWidget *widget, GdkEventScroll *event) {

   int direction = 1;
   if (event->direction == GDK_SCROLL_UP)
      direction = -1;

   graphics_info_t g;
   int imol_scroll = g.scroll_wheel_map;

   if (! g.is_valid_map_molecule(imol_scroll)) {

      std::vector<int> dm = g.displayed_map_imols();
      if (std::find(dm.begin(), dm.end(), imol_scroll) == dm.end()) {
         if (dm.size() > 0)
            imol_scroll = dm[0];
      }
   }

   if (! g.molecules[imol_scroll].is_displayed_p()) {
      // don't scroll the map if the map is not displayed. Scroll the
      // map that *is* displayed
      std::vector<int> dm = g.displayed_map_imols();
      if (dm.size() > 0)
         imol_scroll = dm[0];
   }
 
   if (g.is_valid_map_molecule(imol_scroll)) {
      // use direction
      if (direction == 1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
      if (direction == -1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;
      int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);
      std::cout << "INFO:: contour level for map " << imol_scroll << " is "
                << g.molecules[imol_scroll].contour_level << std::endl;
      g.set_density_level_string(imol_scroll, g.molecules[imol_scroll].contour_level);
      g.display_density_level_this_image = 1;
      // g.update_maps();
      g.graphics_draw(); // queue
   }
   return TRUE;
}

void
graphics_info_t::try_label_unlabel_active_atom() {

   std::pair<int, mmdb::Atom *> aa = get_active_atom();
   int im = aa.first;
   if (im >= 0) {
      mmdb::Atom *at = aa.second;
      if (at) {
         int atom_index;
         // this is a bit convoluteed :-)
         int ierr = at->GetUDData(molecules[im].atom_sel.UDDAtomIndexHandle, atom_index);
	 if (ierr == mmdb::UDDATA_Ok) {
            molecules[im].add_to_labelled_atom_list(atom_index);
            graphics_draw();
         } else {
            std::cout << "WARNING:: Bad UDData for atom_index for atom " << std::endl;
         }
      }
   }
}


gboolean
on_glarea_button_press(GtkWidget *widget, GdkEventButton *event) {

   // std::cout << "button press!" << std::endl;
   graphics_info_t g;
   g.SetMouseBegin(event->x,event->y);
   g.SetMouseClicked(event->x, event->y);
   int x_as_int, y_as_int;
   GdkModifierType mask;
   // gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state); Old-style - keep for grepping
   GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());
   GdkDevice *mouse = gdk_seat_get_pointer(seat);
   gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);

   bool was_a_double_click = false;
   if (event->type==GDK_2BUTTON_PRESS)
      was_a_double_click = true;

   GdkModifierType state;

   if (true) { // check here for left-mouse click
      if (was_a_double_click) {
         pick_info nearest_atom_index_info = atom_pick(event);
         if (nearest_atom_index_info.success == GL_TRUE) {
            int im = nearest_atom_index_info.imol;
            g.molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
            g.graphics_draw();
         }
      }
   }


   g.check_if_in_range_defines(event, mask);
   return TRUE;
}

gboolean
on_glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

   if (event->state & GDK_BUTTON2_MASK) {
      graphics_info_t g;
      pick_info nearest_atom_index_info = atom_pick_gtk3();
      double delta_x = g.GetMouseClickedX() - event->x;
      double delta_y = g.GetMouseClickedY() - event->y;
      if (std::abs(delta_x) < 10.0) {
         if (std::abs(delta_y) < 10.0) {
            if (nearest_atom_index_info.success == GL_TRUE) {
               g.setRotationCentre(nearest_atom_index_info.atom_index,
				                       nearest_atom_index_info.imol);
            }
         }
      }
   }
   return TRUE;
}

void
do_drag_pan_gtk3(GtkWidget *widget) {

   // This should be a graphics_info_t function?

   GtkAllocation allocation;
   gtk_widget_get_allocation(widget, &allocation);
   int w = allocation.width;
   int h = allocation.height;

   graphics_info_t g;
   glm::mat4 mvp = g.get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion

   float mouseX_1 = g.GetMouseBeginX() / (w * 0.5f) - 1.0f;
   float mouseY_1 = g.GetMouseBeginY() / (h * 0.5f) - 1.0f;
   float mouseX_2 = g.mouse_current_x  / (w * 0.5f) - 1.0f;
   float mouseY_2 = g.mouse_current_y  / (h * 0.5f) - 1.0f;

   glm::mat4 vp_inv = glm::inverse(mvp);

   glm::vec4 screenPos_1 = glm::vec4(mouseX_1, -mouseY_1, 1.0f, 1.0f);
   glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2, 1.0f, 1.0f);
   glm::vec4 worldPos_1 = vp_inv * screenPos_1;
   glm::vec4 worldPos_2 = vp_inv * screenPos_2;

   glm::vec4 delta(worldPos_1 - worldPos_2);
   glm::vec3 delta_v3(delta);

   float delta_scale_factor = 1.0;
   if (graphics_info_t::perspective_projection_flag)
      delta_scale_factor = 20.0; // move the front a lot more
   g.add_to_rotation_centre(delta_scale_factor * delta_v3);
   g.update_maps();
   if (graphics_info_t::glareas.size() > 0)
      int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);
}

gboolean
on_glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   int r = 0;
   graphics_info_t g;

   // split this function up before it gets too big.

   bool control_is_pressed = false;
   bool   shift_is_pressed = false;
   if (event->state & GDK_CONTROL_MASK) control_is_pressed = true;
   if (event->state & GDK_SHIFT_MASK) shift_is_pressed = true;

   g.mouse_current_x = event->x;
   g.mouse_current_y = event->y;

   if (event->state & GDK_BUTTON1_MASK) {
      if (control_is_pressed) {
         do_drag_pan_gtk3(widget);
      } else {
         GtkAllocation allocation;
         gtk_widget_get_allocation(widget, &allocation);
         int w = allocation.width;
         int h = allocation.height;
         graphics_info_t::update_view_quaternion(w, h);
      }
   }

   if (event->state & GDK_BUTTON2_MASK) {
      // View Panning
      do_drag_pan_gtk3(widget);
   }

   if (event->state & GDK_BUTTON3_MASK) {
      // Zooming
      double delta_x = event->x - g.GetMouseBeginX();
      double delta_y = event->y - g.GetMouseBeginY();
      double fx = 1.0 + delta_x/300.0;
      double fy = 1.0 + delta_y/300.0;
      if (fx > 0.0) g.zoom /= fx;
      if (fy > 0.0) g.zoom /= fy;
      if (false)
         std::cout << "zooming with perspective_projection_flag "
                   << graphics_info_t::perspective_projection_flag
                   << " " << g.zoom << std::endl;
      if (! graphics_info_t::perspective_projection_flag) {
         // std::cout << "now zoom: " << g.zoom << std::endl;
      } else {
         // Move the eye towards the rotation centre (don't move the rotation centre)

         // this is the same as translate_in_screen_z. Excpet for moving the rotation centre
         glm::vec3 ep = graphics_info_t::get_eye_position();
         glm::vec3 rc = graphics_info_t::get_rotation_centre();
         glm::vec3 delta = rc - ep;
         glm::vec3 delta_uv = normalize(delta);

         // more zoomed needs to have smaller step than when zoomed out.
         // Zoomed out, zoom is ~100. Zoomed in is ~25
         float step_size = 0.005;
         glm::vec3 step = step_size * graphics_info_t::zoom * delta_x * delta_uv;
         // try again
         step_size = 0.05;
         step = step_size * delta_x * delta_uv;

         graphics_info_t::eye_position += step;
         double l = glm::distance(graphics_info_t::eye_position, rc);

         if (false) // debug
            std::cout << "motion: ep " << glm::to_string(ep) << " rc " << glm::to_string(rc)
                      << " zoom " << graphics_info_t::zoom << " step "
                      << glm::to_string(step) << std::endl;

         // if the distance to the centre of rotation changes, but the clipping planes do not, then
         // as we zoom in, then the object get clipped. Bleugh. The perspective near and far
         // need to be dependent on zoom

         float screen_z_near_perspective_limit = l * 0.99;
         float screen_z_far_perspective_limit  = l * 1.01;
         if (graphics_info_t::screen_z_near_perspective > screen_z_near_perspective_limit)
            graphics_info_t::screen_z_near_perspective = screen_z_near_perspective_limit;
         if (graphics_info_t::screen_z_far_perspective < screen_z_far_perspective_limit)
            graphics_info_t::screen_z_far_perspective = screen_z_far_perspective_limit;

      }
   }

   // for next motion
   g.SetMouseBegin(event->x,event->y);
   // gtk_widget_queue_draw(widget);
   g.graphics_draw(); // queue
   return TRUE;
}

gint
view_spin_func(gpointer data) {

   float delta = 0.002;
   glm::vec3 EulerAngles(0, delta, 0);
   glm::quat quat_delta(EulerAngles);
   glm::quat normalized_quat_delta(glm::normalize(quat_delta));
   glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
   graphics_info_t::glm_quat = glm::normalize(product);
   graphics_info_t::graphics_draw(); // queue

   std::chrono::time_point<std::chrono::system_clock> tp_now = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed_seconds = tp_now - graphics_info_t::previous_frame_time_for_per_second_counter;
   if (elapsed_seconds.count() > 1.0) {
      float nf = graphics_info_t::frame_counter - graphics_info_t::frame_counter_at_last_display;
      std::cout << "INFO:: Time/frame: " << 1000 * elapsed_seconds.count()/nf << " milliseconds "
                << nf/elapsed_seconds.count() << " frames/second\n";
      graphics_info_t::previous_frame_time_for_per_second_counter = tp_now;
      graphics_info_t::frame_counter_at_last_display = graphics_info_t::frame_counter;
   }

   // now the stutter checker:
   std::chrono::duration<double> elapsed_seconds_fast = tp_now - graphics_info_t::previous_frame_time;
   if (elapsed_seconds_fast.count() > 0.03)
      std::cout << "INFO:: " << 1000 * elapsed_seconds_fast.count() << " milliseconds for that frame\n";
   graphics_info_t::previous_frame_time = tp_now;

   // kludge/race condition?
   if (graphics_info_t::idle_function_spin_rock_token == -1)
      return FALSE;
   else
      return TRUE;
}

#include "c-interface.h" // for update_go_to_atom_from_current_position()

void
graphics_info_t::translate_in_screen_z(float step_size) {

   // this looks a bit weird without perspective view

   glm::vec3 ep = get_eye_position();
   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   glm::vec3 delta = rc - ep;
   glm::vec3 delta_uv = normalize(delta);

   // more zoomed in has smaller zoom than zoomed out. Zoomed out is ~100. Zoomed in is ~25
   glm::vec3 step = 0.01 * step_size * zoom * delta_uv;

   if (true) // debug
      std::cout << "ep " << glm::to_string(ep) << " rc " << glm::to_string(rc)
                << " zoom " << zoom << " step " << glm::to_string(step) << std::endl;
   graphics_info_t::add_to_rotation_centre(step);

   eye_position += step;

}

void
graphics_info_t::move_forwards() {
   // these are the other way round in perspective - that's interesting.
   translate_in_screen_z(1.0);
}

void
graphics_info_t::move_backwards() {
   translate_in_screen_z(-1.0);
}

void
graphics_info_t::setup_key_bindings() {

   graphics_info_t g;

   // if we are serious about user-defined key-bindings all of these functions should be thunks in the user API
   // (and returning gboolean).

   auto l1 = []() { graphics_info_t g; g.adjust_clipping(-0.1); return gboolean(TRUE); };
   auto l2 = []() { graphics_info_t g; g.adjust_clipping( 0.1); return gboolean(TRUE); };
   auto l5 = []() { graphics_info_t g; g.blob_under_pointer_to_screen_centre(); return gboolean(TRUE); };

   auto l6 = []() {
                if (graphics_info_t::idle_function_spin_rock_token != -1) {
                   std::cout << "Removing the idle function\n";
                   g_idle_remove_by_data(GINT_TO_POINTER(66)); // just a kludge for the moment
                   graphics_info_t::idle_function_spin_rock_token = -1;
                } else {
                   int toi = g_timeout_add(5, view_spin_func, GINT_TO_POINTER(66));
                   graphics_info_t::idle_function_spin_rock_token = toi;
                }
                return gboolean(TRUE);
             };

   auto l7 = []() {
                int imol_scroll = graphics_info_t::scroll_wheel_map;
                if (graphics_info_t::is_valid_map_molecule(imol_scroll))
                   graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
                if (graphics_info_t::glareas.size() > 0)
                   int contour_idle_token = g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
                graphics_info_t g;
                g.set_density_level_string(imol_scroll, graphics_info_t::molecules[imol_scroll].contour_level);
                graphics_info_t::display_density_level_this_image = 1;
                return gboolean(TRUE);
             };

   auto l8 = []() {
                int imol_scroll = graphics_info_t::scroll_wheel_map;
                if (graphics_info_t::is_valid_map_molecule(imol_scroll))
                   graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;
                if (graphics_info_t::glareas.size() > 0)
                   int contour_idle_token = g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
                graphics_info_t g;
                g.set_density_level_string(imol_scroll, graphics_info_t::molecules[imol_scroll].contour_level);
                graphics_info_t::display_density_level_this_image = 1;
                return gboolean(TRUE);
             };

   auto l9 = []() { update_go_to_atom_from_current_position(); return gboolean(TRUE); };

   auto l10 = []() { graphics_info_t::zoom *= 0.9; return gboolean(TRUE); };

   auto l11 = []() { graphics_info_t::zoom *= 1.1; return gboolean(TRUE); };

   auto l12 = []() { graphics_info_t g; g.move_forwards(); return gboolean(TRUE); };

   auto l13 = []() { graphics_info_t g; g.move_backwards(); return gboolean(TRUE); };

   auto l14 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('forward')"); return gboolean(TRUE); };

   auto l15 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('backward')"); return gboolean(TRUE); };

   auto l16 = []() { graphics_info_t g; g.undo_last_move(); return gboolean(TRUE); };

   auto l18 = []() { graphics_info_t g; g.accept_moving_atoms(); return gboolean(TRUE); };

   auto l19 = []() { graphics_info_t g; g.clear_up_moving_atoms_wrapper(); return gboolean(TRUE); };

   auto l20 = []() { graphics_info_t g; g.eigen_flip_active_residue(); return gboolean(TRUE); };

   auto l21 = []() { graphics_info_t g; g.try_label_unlabel_active_atom(); return gboolean(TRUE); };

   // Note to self, Space and Shift Space are key *Release* functions

   std::vector<std::pair<keyboard_key_t, key_bindings_t> > kb_vec;
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_d,      key_bindings_t(l1, "increase clipping")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_f,      key_bindings_t(l2, "decrease clipping")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_g,      key_bindings_t(l5, "go to blob")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_i,      key_bindings_t(l6, "spin")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_plus,   key_bindings_t(l7, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_equal,  key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_minus,  key_bindings_t(l8, "decrease contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_p,      key_bindings_t(l9, "update go-to atom by position")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_n,      key_bindings_t(l10, "Zoom in")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_m,      key_bindings_t(l11, "Zoom out")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_w,      key_bindings_t(l12, "Move forward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_s,      key_bindings_t(l13, "Move backward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l14, "NCS Skip forward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_O,      key_bindings_t(l15, "NCS Skip backward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_u,      key_bindings_t(l16, "Undo Move")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Return, key_bindings_t(l18, "Accept Moving Atoms")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Escape, key_bindings_t(l19, "Reject Moving Atoms")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_e,      key_bindings_t(l20, "EigenFlip Active Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_l,      key_bindings_t(l21, "Label/Unlabel Active Atom")));

   // control keys

   auto lc1 = []() { show_go_to_residue_keyboarding_mode_window(); return gboolean(TRUE); };
   key_bindings_t go_to_residue_key_binding(lc1, "Show Go To Residue Keyboarding Window");
   std::pair<keyboard_key_t, key_bindings_t> p1(keyboard_key_t(GDK_KEY_g, true), go_to_residue_key_binding);
   kb_vec.push_back(p1);

   auto lc2 = []() { graphics_info_t g; g.apply_undo(); return gboolean(TRUE); };
   key_bindings_t undo_key_binding(lc2, "Undo");
   std::pair<keyboard_key_t, key_bindings_t> p2(keyboard_key_t(GDK_KEY_z, true), undo_key_binding);
   kb_vec.push_back(p2);

   auto lc3 = []() { graphics_info_t g; g.apply_redo(); return gboolean(TRUE);};
   key_bindings_t redo_key_binding(lc3, "Redo");
   std::pair<keyboard_key_t, key_bindings_t> p3(keyboard_key_t(GDK_KEY_y, true), redo_key_binding);
   kb_vec.push_back(p3);

   // left
   auto lc4 = []() {
                 if (graphics_info_t::control_is_pressed) {
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Left);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Left);
                 } else {
                    keypad_translate_xyz(1, 1);
                 }
                 return gboolean(TRUE);
              };

   // right
   auto lc5 = []() {
                 if (graphics_info_t::control_is_pressed) {
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Right);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Right);
                 } else {
                    keypad_translate_xyz(1, -1);
                 }
                 return gboolean(TRUE);
              };

   // up
   auto lc6 = []() {
                 if (graphics_info_t::control_is_pressed) {
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Up);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Up);
                 } else {
                    keypad_translate_xyz(2, 1);
                 }
                 return gboolean(TRUE);
              };
   // down
   auto lc7 = []() {
                 if (graphics_info_t::control_is_pressed) {
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Down);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Down);
                 } else {
                    keypad_translate_xyz(2, -1);
                 }
                 return gboolean(TRUE);
              };

   key_bindings_t ctrl_arrow_left_key_binding(lc4, "R/T Left");
   key_bindings_t ctrl_arrow_right_key_binding(lc5, "R/T Right");
   key_bindings_t ctrl_arrow_up_key_binding(lc6, "R/T Up");
   key_bindings_t ctrl_arrow_down_key_binding(lc7, "R/T Down");
   std::pair<keyboard_key_t, key_bindings_t> p4(keyboard_key_t(GDK_KEY_Left,  true), ctrl_arrow_left_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p5(keyboard_key_t(GDK_KEY_Right, true), ctrl_arrow_right_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p6(keyboard_key_t(GDK_KEY_Up,    true), ctrl_arrow_up_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p7(keyboard_key_t(GDK_KEY_Down,  true), ctrl_arrow_down_key_binding);
   kb_vec.push_back(p4);
   kb_vec.push_back(p5);
   kb_vec.push_back(p6);
   kb_vec.push_back(p7);

   std::vector<std::pair<keyboard_key_t, key_bindings_t> >::const_iterator it;
   for (it=kb_vec.begin(); it!=kb_vec.end(); it++)
     g.key_bindings_map[it->first] = it->second;

}



gboolean
on_glarea_key_press_notify(GtkWidget *widget, GdkEventKey *event) {

   // move this function into graphics_info_t?

   graphics_info_t g;
   gboolean handled = false;

   // "space" and "shift space" have the same keyval. So ctrl and shift handling are different.
   bool control_is_pressed_flag = false;
   g.shift_is_pressed = false;
   if (event->state & GDK_CONTROL_MASK) control_is_pressed_flag = true;
   if (event->state & GDK_SHIFT_MASK) g.shift_is_pressed = true;
   if (event->keyval == GDK_KEY_Shift_L) g.shift_is_pressed = true;


   keyboard_key_t kbk(event->keyval, control_is_pressed_flag);

   std::map<keyboard_key_t, key_bindings_t>::const_iterator it = g.key_bindings_map.find(kbk);

   if (it != g.key_bindings_map.end()) {
     const key_bindings_t &kb = it->second;
     if (true)
        std::cout << "key-binding for key " << it->first.gdk_key << " "
                  << it->first.ctrl_is_pressed << " " << kb.description << std::endl;
     kb.run();
   } else {
      std::cout << "on_glarea_key_press_notify() key not found in map: " << event->keyval << std::endl;
   }


   // fix the type here
   if (int(event->keyval) == graphics_info_t::update_go_to_atom_from_current_residue_key) {
      update_go_to_atom_from_current_position();
      handled = TRUE;
   }

   graphics_info_t::graphics_draw(); // queue
   // gtk_widget_queue_draw(widget);

   return handled;

}

gboolean
on_glarea_key_release_notify(GtkWidget *widget, GdkEventKey *event) {

   graphics_info_t g;

   // We need to check the GDK_KEY_Shift_R also, I guess. Not clear
   // to me how to do that now. Fix later.
   g.shift_is_pressed = false;
   if (event->state & GDK_SHIFT_MASK) g.shift_is_pressed = true;
   if (event->keyval == GDK_KEY_Shift_L) g.shift_is_pressed = true;

   if (event->keyval == GDK_KEY_space) {
      g.reorienting_next_residue_mode = false; // hack
      bool reorienting = graphics_info_t::reorienting_next_residue_mode;
      if (reorienting) {
         if (graphics_info_t::shift_is_pressed) {
            g.reorienting_next_residue(false); // backwards
         } else {
            std::cout << "Do forward reorienting" << std::endl;
            g.reorienting_next_residue(true); // forwards
         }
      } else {
         // old/standard simple translation
         if (graphics_info_t::shift_is_pressed) {
            g.intelligent_previous_atom_centring(g.go_to_atom_window);
         } else {
            g.intelligent_next_atom_centring(g.go_to_atom_window);
         }
      }
   }
   return TRUE;
}

void
my_glarea_add_signals_and_events(GtkWidget *glarea) {

   gtk_widget_add_events(glarea, GDK_SCROLL_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_PRESS_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_RELEASE_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_POINTER_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MASK);
   gtk_widget_add_events(glarea, GDK_KEY_PRESS_MASK);

   // key presses for the glarea:

   gtk_widget_set_can_focus(glarea, TRUE);
   gtk_widget_grab_focus(glarea);

   g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
   g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
   g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);
   g_signal_connect(glarea, "scroll-event",          G_CALLBACK(on_glarea_scroll),             NULL);
   g_signal_connect(glarea, "button-press-event",    G_CALLBACK(on_glarea_button_press),       NULL);
   g_signal_connect(glarea, "button-release-event",  G_CALLBACK(on_glarea_button_release),     NULL);
   g_signal_connect(glarea, "motion-notify-event",   G_CALLBACK(on_glarea_motion_notify),      NULL);
   g_signal_connect(glarea, "key-press-event",       G_CALLBACK(on_glarea_key_press_notify),   NULL);
   g_signal_connect(glarea, "key-release-event",     G_CALLBACK(on_glarea_key_release_notify), NULL);

}
