/*
 * src/graphics-info-draw.cc
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
#include <Python.h>
#endif // USE_PYTHON

#include "compat/coot-sysdep.h"

#include "glm/matrix.hpp"
#define GLM_ENABLE_EXPERIMENTAL // # for norm things
// #include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/type_ptr.hpp>  // for value_ptr() 20240326-PE

#include <iostream>
#include <string>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "c-interface.h" // for update_go_to_atom_from_current_position()
#include "graphics-info.h"

#include "draw-2.hh"
#include "framebuffer.hh"

#include "text-rendering-utils.hh"
#include "cc-interface-scripting.hh"
#include "coot-utils/cylinder-with-rotation-translation.hh"

#include "screendump-tga.hh"
#include "widget-from-builder.hh"

#include "utils/logging.hh"
extern logging logger;


enum {VIEW_CENTRAL_CUBE, ORIGIN_CUBE};


glm::vec3
get_camera_up_direction(const glm::mat4 &mouse_quat_mat) {

   glm::vec4 z_p(0.0f, 1.0f, 0.0f, 1.0f);
   glm::vec4 r = z_p * mouse_quat_mat;
   glm::vec3 r3(r);
   return r3;
}

float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
      // positions   // texCoords
      -1.0f,  1.0f,  0.0f, 1.0f,
      -1.0f, -1.0f,  0.0f, 0.0f,
       1.0f, -1.0f,  1.0f, 0.0f,

      -1.0f,  1.0f,  0.0f, 1.0f,
       1.0f, -1.0f,  1.0f, 0.0f,
       1.0f,  1.0f,  1.0f, 1.0f
};


void
graphics_info_t::init_screen_quads() {

   // screen quad VAO
   unsigned int quadVBO;
   glGenVertexArrays(1, &screen_quad_vertex_array_id);
   glBindVertexArray(screen_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   GLenum err = glGetError();
   if (err) std::cout << "init_screen_quads() A err is " << err << std::endl;

   glGenVertexArrays(1, &blur_y_quad_vertex_array_id);
   glBindVertexArray(blur_y_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   err = glGetError();
   if (err) std::cout << "init_screen_quads() B err is " << err << std::endl;

   glGenVertexArrays(1, &blur_x_quad_vertex_array_id);
   glBindVertexArray(blur_x_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   err = glGetError();
   if (err) std::cout << "init_screen_quads() C err is " << err << std::endl;

   glGenVertexArrays(1, &combine_textures_using_depth_quad_vertex_array_id);
   glBindVertexArray(combine_textures_using_depth_quad_vertex_array_id);
   glGenBuffers(1, &quadVBO);
   glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
   glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
   err = glGetError();
   if (err) std::cout << "init_screen_quads() D err is " << err << std::endl;

}
void
graphics_info_t::init_blur_quads() {

   // This is the 2020 version - this function can be deleted when the new version works.

   // graphics_info_t::shader_for_blur.Use(); // setting up buffers doesn't need a shader
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
   if (err) std::cout << "init_blur_quads() err is " << err << std::endl;

}

// void graphics_info_t::init_central_cube();

void
graphics_info_t::init_buffers() {
   // std::cout << "debug:: init_buffers() init_central_cube" << std::endl;
   init_central_cube();
   // std::cout << "debug:: init_buffers() init_screen_quads" << std::endl;
   init_screen_quads();
   // std::cout << "debug:: init_buffers() init_blur_quads" << std::endl;
   init_blur_quads();
}

void
graphics_info_t::init_central_cube() {

   float cube_positions[24] = {
                          -0.5,  -0.5, -0.5,
                          -0.5,  -0.5,  0.5,
                          -0.5,   0.5, -0.5,
                          -0.5,   0.5,  0.5,
                           0.5,  -0.5, -0.5,
                           0.5,  -0.5,  0.5,
                           0.5,   0.5, -0.5,
                           0.5,   0.5,  0.5
   };

   float crosshair_positions[18] = {
                                    -0.5f,  0.0f,  0.0,
                                     0.5f,  0.0f,  0.0,
                                     0.0f, -0.5f,  0.0,
                                     0.0f,  0.5f,  0.0,
                                     0.0f,  0.0f, -0.5,
                                     0.0f,  0.0f,  0.5
   };

   // graphics_info_t::shader_for_central_cube.Use(); setting up buffers doesn't need a shader
   GLenum err = glGetError();
   if (err) std::cout << "init_central_cube() glUseProgram() err is " << err << std::endl;

   // number of lines * 2:
   unsigned int cube_indices[24] { 0,1, 1,5, 5,4, 4,0, 2,3, 3,7, 7,6, 6,2, 0,2, 1,3, 5,7, 4,6 };

   unsigned int crosshair_indices[6] = {0,1,2,3,4,5};

   // GLuint VertexArrayID;
   glGenVertexArrays(1, &graphics_info_t::central_cube_vertexarray_id);
   glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);

   // GLuint vertexbuffer;
   glGenBuffers(1, &graphics_info_t::central_cube_array_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::central_cube_array_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 24, &cube_positions[0], GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

   // unsigned int ibo;
   glGenBuffers(1, &graphics_info_t::central_cube_index_buffer_id);
   err = glGetError();
   if (err) std::cout << "init_central_cube() index glGenBuffers() err is " << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::central_cube_index_buffer_id);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 24, &cube_indices[0], GL_STATIC_DRAW);
   err = glGetError();
   if (err) std::cout << "init_central_cube() glBufferData() err is " << err << std::endl;
   glBindVertexArray(0);

   // now the crosshairs

   glGenVertexArrays(1, &graphics_info_t::rotation_centre_crosshairs_vertexarray_id);
   glBindVertexArray(graphics_info_t::rotation_centre_crosshairs_vertexarray_id);
   // positions
   glGenBuffers(1, &graphics_info_t::rotation_centre_crosshairs_vertex_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::rotation_centre_crosshairs_vertex_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 18, &crosshair_positions[0], GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

   // indices
   glGenBuffers(1, &graphics_info_t::rotation_centre_crosshairs_index_buffer_id);
   err = glGetError();
   if (err) std::cout << "init_central_cube() index buffer glGenBuffers() for crosshairs A err is "
                      << err << std::endl;
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::rotation_centre_crosshairs_index_buffer_id);
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 6, &crosshair_indices[0], GL_STATIC_DRAW);
   if (err) std::cout << "init_central_cube() index buffer glGenBuffers() for crosshairs B err is "
                      << err << std::endl;
   glBindVertexArray(0);

}

void
graphics_info_t::init_hud_text() {

   // std::cout << ":::::::::::: init_hud_text() " << std::endl;

   graphics_info_t g;
   // g.load_freetype_font_textures(); 20220226-PE already done by now.
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
}


void
graphics_info_t::handle_delete_item_curor_change(GtkWidget *widget) {

   if (delete_item_widget) {
      if (delete_item_water) {
         graphics_info_t g;
	 pick_info naii = g.atom_pick_gtk3(false);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME 
#else
         GdkDisplay *display = gdk_display_get_default();
         GdkWindow *window = 0;
         window = gtk_widget_get_window(GTK_WIDGET(widget));
         if (window) {
            // GdkCursor *current_cursor = gdk_window_get_cursor(window);
            // std::cout << "current cursor " << gdk_cursor_get_cursor_type(current_cursor) << std::endl;
            if (naii.success == GL_TRUE) {
               int imol = naii.imol;
               molecule_class_info_t &m = graphics_info_t::molecules[imol];
               std::string res_name = m.atom_sel.atom_selection[naii.atom_index]->GetResName();
               if (res_name == "HOH") {
                  GdkCursor *c = gdk_cursor_new_from_name (display, "crosshair");
                  // std::cout << "crosshair type " << gdk_cursor_get_cursor_type(c) << std::endl;
                  gdk_window_set_cursor(window, c);
               } else {
                  GdkCursor *c = gdk_cursor_new_from_name (display, "not-allowed");
                  // std::cout << "not-allowed type " << gdk_cursor_get_cursor_type(c) << std::endl;
                  gdk_window_set_cursor(window, c);
               }
            } else {
               GdkCursor *c = gdk_cursor_new_from_name (display, "not-allowed");
               // std::cout << "not-allowed type " << gdk_cursor_get_cursor_type(c) << std::endl;
               gdk_window_set_cursor(window, c);
            }
         }
#endif
      }
   }
}

// Called by pinch zoom gesture
// static
void
graphics_info_t::mouse_zoom_by_scale_factor(double sf) {

   // with a "long" pinch gesture, sf can be 0.4
   // so using sf directly is not what we want.

   // So try adding a fixed zoom depending on which size
   // of 1.0 sf is
   float zf = 1.0;
   if (sf > 1.0) zf = 1.05;
   if (sf < 1.0) zf = 0.95;

   zoom /= zf;
   // sensible limits for looking at proteins
   if (zoom <    0.2) zoom = 0.2;
   if (zoom > 2000.0) zoom = 2000.0;
   // std::cout << "debug:: mouse_zoom_by_scale_factor() sf " << sf << " zoom " << zoom << std::endl;
   //  mouse_zoom_by_scale_factor_inner(sf);
   graphics_draw();

}


// static
void
graphics_info_t::mouse_zoom_by_scale_factor_inner(double sf) {

   if (perspective_projection_flag) {

      { // own graphics_info_t function - c.f. adjust clipping

         eye_position.z *= sf;

         double  l = graphics_info_t::eye_position.z;
         double zf = graphics_info_t::screen_z_far_perspective;
         double zn = graphics_info_t::screen_z_near_perspective;

         graphics_info_t::screen_z_near_perspective *= sf;
         graphics_info_t::screen_z_far_perspective  *= sf;

         float screen_z_near_perspective_limit = l * 0.95;
         float screen_z_far_perspective_limit  = l * 1.05;
         if (graphics_info_t::screen_z_near_perspective < 2.0)
            graphics_info_t::screen_z_near_perspective = 2.0;
         if (graphics_info_t::screen_z_far_perspective > 1000.0)
            graphics_info_t::screen_z_far_perspective = 1000.0;

         if (graphics_info_t::screen_z_near_perspective > screen_z_near_perspective_limit)
            graphics_info_t::screen_z_near_perspective = screen_z_near_perspective_limit;
         if (graphics_info_t::screen_z_far_perspective < screen_z_far_perspective_limit)
            graphics_info_t::screen_z_far_perspective = screen_z_far_perspective_limit;
         if (false)
            std::cout << "mouse_zoom_by_scale_factor_inner(): debug l: " << l << " post-manip: "
                      << graphics_info_t::screen_z_near_perspective << " "
                      << graphics_info_t::screen_z_far_perspective << std::endl;
      }

   } else {

      // stabilize the scale factor
      if (sf < 0.5) sf = 0.5;
      if (sf > 2.0) sf = 2.0;
      graphics_info_t::eye_position.z *= sf;

   }
}

// static
void
graphics_info_t::mouse_zoom(double delta_x_drag, double delta_y_drag) {

   double current_mouse_x = drag_begin_x + delta_x_drag;
   double current_mouse_y = drag_begin_y + delta_y_drag;

   double delta_x = current_mouse_x - get_mouse_previous_position_x();
   double delta_y = current_mouse_y - get_mouse_previous_position_y();

   // Zooming
   double fx = 1.0 + delta_x/300.0;
   double fy = 1.0 + delta_y/300.0;
   if (fx > 0.0) graphics_info_t::zoom /= fx;
   if (fy > 0.0) graphics_info_t::zoom /= fy;
   if (false)
      std::cout << "zooming with perspective_projection_flag "
                << graphics_info_t::perspective_projection_flag
                << " " << graphics_info_t::zoom << std::endl;

   if (perspective_projection_flag) {

      // std::cout << "now zoom: " << screen_z_near_perspective << " " << screen_z_far_perspective << std::endl;
      if (fabs(delta_y) > fabs(delta_x))
         delta_x = delta_y;
      float sf = 1.0 - delta_x * 0.003;
      mouse_zoom_by_scale_factor_inner(sf);

   } else {

      // Move the eye towards the rotation centre (don't move the rotation centre)
      if (fabs(delta_y) > fabs(delta_x))
         delta_x = delta_y;
      float sf = 1.0 - delta_x * 0.003;
      mouse_zoom_by_scale_factor_inner(sf);
   }
   graphics_draw(); // or should this be called by the function that calls this function?
}

// static
void
graphics_info_t::scroll_zoom(int direction) {

   // c.f. mouse_zoom() it was copied from there. I am not sure that I like this yet.
   // and don't want to refactor if I'm going to dump either this or that later.

   // scroll up mean direction -1

   // Zooming
   double delta_x = 15.0;
   if (direction == 1) delta_x = -delta_x;
   double fx = 1.0 + delta_x/300.0;
   if (fx > 0.0) graphics_info_t::zoom /= fx;
   if (false)
      std::cout << "zooming with perspective_projection_flag "
                << graphics_info_t::perspective_projection_flag
                << " " << graphics_info_t::zoom << std::endl;
   if (! graphics_info_t::perspective_projection_flag) {
      // std::cout << "now zoom: " << g.zoom << std::endl;
   } else {
      // Move the eye towards the rotation centre (don't move the rotation centre)
      float sf = 1.0 - delta_x * 0.003;
      graphics_info_t::eye_position.z *= sf;

      { // own graphics_info_t function - c.f. adjust clipping
         double  l = graphics_info_t::eye_position.z;

         graphics_info_t::screen_z_near_perspective *= sf;
         graphics_info_t::screen_z_far_perspective  *= sf;

         float screen_z_near_perspective_limit = l * 0.95;
         float screen_z_far_perspective_limit  = l * 1.05;
         if (graphics_info_t::screen_z_near_perspective < 2.0)
            graphics_info_t::screen_z_near_perspective = 2.0;
         if (graphics_info_t::screen_z_far_perspective > 1000.0)
            graphics_info_t::screen_z_far_perspective = 1000.0;

         if (graphics_info_t::screen_z_near_perspective > screen_z_near_perspective_limit)
            graphics_info_t::screen_z_near_perspective = screen_z_near_perspective_limit;
         if (graphics_info_t::screen_z_far_perspective < screen_z_far_perspective_limit)
            graphics_info_t::screen_z_far_perspective = screen_z_far_perspective_limit;
         if (false)
            std::cout << "on_glarea_motion_notify(): debug l: " << l << " post-manip: "
                      << graphics_info_t::screen_z_near_perspective << " "
                      << graphics_info_t::screen_z_far_perspective << std::endl;
      }
   }
   graphics_draw(); // or should this be called by the function that calls this function?
}

// static
glm::mat4
graphics_info_t::get_model_matrix() {

   // rework this so that the view rotation matrix becomes part of the *model* matrix
   // for both orthographic and perspective

   glm::mat4 m(1.0f);
   glm::vec3 rc = get_rotation_centre();
   glm::mat4 model_translation_matrix = glm::translate(m, -rc);
   glm::mat4 model_rotation_matrix = get_model_rotation();
   glm::mat4 model_matrix = model_rotation_matrix * model_translation_matrix;
   return model_matrix;
}

// static
glm::mat4
graphics_info_t::get_view_matrix() { // the lookAt() matrix

   // make this static

   // this function reworked this so that the view rotation matrix became part of the *model* matrix
   // for both orthograph and perspective

      glm::vec3 ep = eye_position; // in view space i.e. (0,0,z) (z = 40, say)

      glm::vec3 ep_ws = get_world_space_eye_position();
      glm::vec3 rc = get_rotation_centre();
      glm::vec3 t_ep_ws = ep_ws - rc;
      // glm::vec3 rc = glm::vec3(0,0,0);
      rc = glm::vec3(0,0,0);
      glm::vec3 up(0,1,0);
      glm::mat4 view_matrix = glm::lookAt(ep, rc, up);
      // std::cout << "debug:: get_view_matrix() t_ep_ws " << glm::to_string(t_ep_ws) << " rc " << glm::to_string(rc) << std::endl;
      return view_matrix;
}

// static
glm::mat4
graphics_info_t::get_projection_matrix(bool do_orthographic_projection,
                                       int graphics_x_size, int graphics_y_size) {

   float w = static_cast<float>(graphics_x_size);
   float h = static_cast<float>(graphics_y_size);
   float screen_ratio = static_cast<float>(w)/static_cast<float>(h);
   if (do_orthographic_projection) {

      // 20240814-PE if clipping front and back are somehow zero (I don't know
      // how that happened, but those were the values in the state script)
      // the we get nan and infs in matrices
      //
      if (clipping_front < 0.00001) clipping_front = 0.00001;
      if (clipping_back  < 0.00001) clipping_back  = 0.00001;

      float sr = screen_ratio;
      GLfloat near =  -0.1 * zoom * clipping_front + eye_position.z;
      GLfloat far  =   0.3 * zoom * clipping_back  + eye_position.z;

      glm::mat4 projection_matrix = glm::ortho(-0.3f*zoom*sr, 0.3f*zoom*sr,
                                               -0.3f*zoom,    0.3f*zoom,
                                               near, far);
      if (false) {
         std::cout << "debug:: get_projection_matrix() near " << near << " far " << far
                   << " clipping-front: " << clipping_front << " clipping_back: " << clipping_back << " "
                   << "eye_position " << glm::to_string(eye_position)
                   << " zoom " << zoom << std::endl;
         std::cout << "projection matrix ortho " << glm::to_string(projection_matrix) << std::endl;
      }
      return projection_matrix;
   } else {
      // perspective_fov is in degrees
      glm::mat4 projection_matrix_persp = glm::perspective(glm::radians(perspective_fov),
                                                           screen_ratio,
                                                           screen_z_near_perspective,
                                                           screen_z_far_perspective);
      // std::cout << "projection matrix persp " << glm::to_string(projection_matrix_persp)
      // << " " << graphics_x_size << " " << graphics_y_size << std::endl;
      return projection_matrix_persp;
   }
}

glm::mat4
graphics_info_t::get_mvp_for_shadow_map(const glm::vec3 &light_direction_eye_space) const {

   // In draw(), the light direction is converted to molecular space in the shader (setup_light())
   //    glm::vec4 p4   = glm::vec4(light.direction,1.0) * view_rotation_matrix;
   //    glm::vec3 direction_in_molecule_coordinates_space = glm::vec3(p4);
   // 20211116-PE I think this is wrong. which is why, in texture-meshes.shader we (need to) do:
   //          mat4 ivr = transpose(view_rotation);
   //          vec3 light_dir = (vec4(light_sources[i].direction_in_molecule_coordinates_space, 1.0) * ivr).xyz;
   // i.e. the rotation matrix inversion could and should happen in setup_light().
   //      Let's come back to that later.

   // 20220217-PE how long does this function take to run? It is called many times per frame.

   // conversion from crows
   glm::vec3 rotation_centre = get_rotation_centre();

   glm::mat4 model_rotation = get_model_rotation();
   glm::mat3 model_rotation_mat3 = glm::mat3(model_rotation);
   glm::mat3 inv_model_rotation_mat3 = glm::transpose(model_rotation_mat3);

   // glm::vec4 light_position_world_space = glm::vec4(light_position_eye_space, 1.0) * view_rotation;

   glm::vec3 light_direction_world_space = inv_model_rotation_mat3 * (light_direction_eye_space);

   glm::mat4 model_matrix = glm::mat4(1.0);

   // 20231119-PE Calculate this from the extents of the displayed molecules? suggest_shadow_box_size() ?
   //             Not every frame though. Hmm.
   float box_size = shadow_box_size; // user setable, default 66.
   if (box_size < 0.0) box_size = 120.0;
   glm::mat4 projection_matrix = glm::ortho(-box_size, box_size, -box_size, box_size, -box_size, box_size);

   glm::vec3 rc = get_rotation_centre();
   glm::vec3 light_position = 1.0f * light_direction_world_space + rc;

   glm::vec3 y_screen = get_screen_y_uv();
   glm::mat4 light_position_matrix = glm::lookAt(light_position, rotation_centre, y_screen);

   glm::vec3 eye_position_ws = get_world_space_eye_position();

   if (false)
      std::cout << "debug:: get_mvp_for_shadow_map() rotation-centre: " << glm::to_string(rc)
                << " light_position: " << glm::to_string(light_position)
                << " eye_position: "   << glm::to_string(eye_position_ws)
                << std::endl;

   glm::mat4 mvp = projection_matrix * light_position_matrix * model_matrix;
   return mvp;

}

// static
glm::mat4
graphics_info_t::get_light_space_mvp(int light_index) {

   glm::mat4 m(1.0f);
   std::map<unsigned int, lights_info_t>::const_iterator it;
   it = lights.find(light_index);
   if (it != lights.end()) {
      graphics_info_t g;
      glm::vec3 dir = it->second.direction;
      m = g.get_mvp_for_shadow_map(dir); // make this static?
   } else {
      std::cout << "ERROR:: get_light_space_mvp() bad light index " << light_index << std::endl;
   }
   return m;
}




// static
glm::mat4
graphics_info_t::get_molecule_mvp(stereo_eye_t eye, bool debug_matrices) {

   int w = graphics_x_size;
   int h = graphics_y_size;

   if (false) {  // debug problematic matrices - get rid of this, make sure that it doesn't do anything
      GtkAllocation allocation;
      gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
      w = allocation.width;
      h = allocation.height;
   }

   if (scale_up_graphics != 1) {
      w *= scale_up_graphics;
      h *= scale_up_graphics;
   }
   if (scale_down_graphics != 1) {
      w /= scale_down_graphics;
      h /= scale_down_graphics;
   }
   // std::cout << scale_up_graphics << " " << scale_down_graphics << " " << w << " " << h << std::endl;

   bool do_orthographic_projection = ! perspective_projection_flag; // weird
   glm::mat4  view_matrix =       get_view_matrix();
   glm::mat4 model_matrix =      get_model_matrix();
   glm::mat4  proj_matrix = get_projection_matrix(do_orthographic_projection, w, h);

   // this setup is for cross-eye : i.e left image is on the right.
   // we use the stereo mode (wall vs cross) to get the correct sign
   //
   float angle = -3.0f; // degrees, default cross-eye
   if (display_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) angle = -angle;
   if (eye == stereo_eye_t::LEFT_EYE) {
      glm::mat4 rot_z = glm::rotate(glm::mat4(1.0f), glm::radians(angle), glm::vec3(0.0f, 1.0f, 0.0f));
      view_matrix = rot_z * view_matrix;
   }
   if (eye == stereo_eye_t::RIGHT_EYE) {
      glm::mat4 rot_z = glm::rotate(glm::mat4(1.0f), glm::radians(-angle), glm::vec3(0.0f, 1.0f, 0.0f));
      view_matrix = rot_z * view_matrix;
   }

   glm::mat4 mvp = proj_matrix * view_matrix * model_matrix;

   if (false) {
      std::cout << "debug:: in get_molecule_mvp() model " << glm::to_string(model_matrix) << std::endl;
      std::cout << "debug:: in get_molecule_mvp() view  " << glm::to_string(view_matrix)  << std::endl;
      std::cout << "debug:: in get_molecule_mvp() proj  " << glm::to_string(proj_matrix)  << std::endl;
      std::cout << "debug:: in get_molecule_mvp() mvp   " << glm::to_string(mvp)          << std::endl;
   }
   return mvp;

}

// can we work out the eye position without needing to unproject? (because that depends
// on get_molecule_mvp()...
//
glm::vec3
graphics_info_t::get_world_space_eye_position() {

   if (! graphics_info_t::perspective_projection_flag) {

      // orthograph eye position is inferred from centre position
      // and zoom and view rotation.

      // does this work? How can I tell?

      glm::vec3 test_vector_1(0.0, 0.0, 1.0);
      glm::vec3 test_vector_2(1.0, 1.0, 0.0);

      glm::mat4 mr = get_model_rotation();
      glm::vec4 rot_test_vector_1 = glm::vec4(test_vector_1, 1.0) * mr;
      glm::vec4 rot_test_vector_2 = glm::vec4(test_vector_2, 1.0) * mr;

      glm::vec3 ep = graphics_info_t::zoom * glm::vec3(rot_test_vector_1);
      glm::vec3 rc = graphics_info_t::get_rotation_centre();
      ep += rc;

      return ep;

   } else {

      // the eye_position is in view-coordinates is stored directly and is
      // by default and often (always?) (0,0,40).

      // I need to convert that to world coordinates and then rotate
      // and translate the world according to rotation centre and mouse-based
      // quaternion

      glm::vec3 ep = eye_position;
      glm::vec4 ep_4(ep, 1.0);
      glm::vec3 up(0,1,0);
      glm::vec3 origin(0,0,0);

      glm::mat4 trackball_matrix = glm::toMat4(graphics_info_t::view_quaternion);
      glm::vec3 rc = graphics_info_t::get_rotation_centre();
      glm::mat4 model_matrix = glm::mat4(1.0);
      model_matrix = glm::translate(model_matrix, -rc);
      model_matrix = trackball_matrix * model_matrix;
      glm::mat4 model_inv = glm::inverse(model_matrix);

      glm::mat4 view_matrix = glm::lookAt(ep, origin, up);
      glm::mat4 view_inv = glm::inverse(view_matrix);
      glm::vec4 ep_world_4 = model_inv * ep_4; // * view_inv;

      if (false) {
         std::cout << "model_inv " << glm::to_string(model_inv) << std::endl;
         std::cout << "view_inv  " << glm::to_string( view_inv) << std::endl;
      }

      glm::vec3 ep_world(ep_world_4);

      return ep_world;
   }

}

glm::vec4
graphics_info_t::unproject(float z) {

   // z is 1 and -1 for front and back (or vice verse).
 
   stereo_eye_t eye = stereo_eye_t::MONO;

   GtkAllocation allocation;
   gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
   float w = allocation.width;
   float h = allocation.height;
   graphics_info_t g;
   float mouseX = 2.0 *    g.GetMouseBeginX()/w  - 1.0f;
   float mouseY = 2.0 * (h-g.GetMouseBeginY())/h - 1.0f;
   std::cout << "debug in new_unproject widget w and h " << w << " " << h << std::endl;
   std::cout << "debug in new_unproject mouse x and y widget  "
             << g.GetMouseBeginX() << " "
             << g.GetMouseBeginY() << std::endl;
   std::cout << "debug in new_unproject mouse x and y GL      " << mouseX << " " << mouseY << std::endl;
   glm::mat4 mvp = get_molecule_mvp(eye);
   glm::mat4 vp_inv = glm::inverse(mvp);
   float real_y = - mouseY; // in range -1 -> 1
   glm::vec4 screenPos_f = glm::vec4(mouseX, real_y, z, 1.0f);
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   std::cout << "debug in new_unproject() screen_pos " << glm::to_string(screenPos_f) << std::endl;
   std::cout << "debug in new_unproject() world_pos " << glm::to_string(worldPos_f) << std::endl;
   return worldPos_f;

}


// the mouse-based quaternion now rotates the model, not the view!
glm::mat4
graphics_info_t::get_model_rotation() {

   // need to be in the correct program (well, the model-drawing part)

   return glm::toMat4(graphics_info_t::view_quaternion);

}


// If the next time you want to use this, but don't have access to graphics_info_t, then move
// this function outside the graphics_info_t class - mabye it's own file/header!
//
void
graphics_info_t::myglLineWidth(int n_pixels) {

#ifdef __APPLE__

   GLint range[2];
   glGetIntegerv(GL_ALIASED_LINE_WIDTH_RANGE, range);
   if (n_pixels < range[1])
      glLineWidth(n_pixels);
   else
      glLineWidth(range[1]);
#else
   glLineWidth(n_pixels);
#endif
   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: in myglLineWidth()  " << n_pixels << " " << err << std::endl;
}

void
graphics_info_t::draw_map_molecules(stereo_eye_t eye, bool draw_transparent_maps) {

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: g.draw_map_molecules() -- start -- " << err << std::endl;

   // run through this molecule loop twice - for opaque then transparent maps
   // first, a block that decides if we need to do anything.

   bool needs_blend_reset = false;

   //

   unsigned int n_transparent_maps = 0;
   unsigned int n_maps_to_draw = 0;
   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
      const molecule_class_info_t &m = graphics_info_t::molecules[ii];
      if (draw_transparent_maps) {
         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         if (! m.draw_it_for_map) continue;
         if (! m.is_an_opaque_map()) {
            n_transparent_maps++;
            n_maps_to_draw += 1;
         }
      } else {
         if (m.is_an_opaque_map()) {
            if (m.draw_it_for_map)
               n_maps_to_draw += 1;
         }
      }
   }

   if (n_maps_to_draw == 0) return;

   if (n_transparent_maps > 0) {
      needs_blend_reset = true;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   bool cosine_dependent_map_opacity = true;

   if (cosine_dependent_map_opacity) {
      needs_blend_reset = true;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   err = glGetError();
   if (err) std::cout << "gtk3_draw_map_molecules() A " << err << std::endl;

   if (!draw_transparent_maps || n_transparent_maps > 0) {

      myglLineWidth(map_line_width * framebuffer_scale);
      err = glGetError();
      if (err) std::cout << "gtk3_draw_map_molecules() glLineWidth " << err << std::endl;


      Shader &shader = shader_for_meshes;
      shader.Use(); // needed? I think not.

      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation = get_model_rotation();

      glEnable(GL_DEPTH_TEST); // this needs to be in the draw loop!?
      glDepthFunc(GL_LESS);
      glDisable(GL_BLEND); // 20220211-PE testing (where has the previous map gone?)
                           // Hmm.. - seems not to be it.
      glm::vec4 ep(get_world_space_eye_position(), 1.0);
      glm::vec3 ep3 = ep/ep.w;

      for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {

         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         molecule_class_info_t &m = graphics_info_t::molecules[ii]; // not const because shader changes

         if (false)
            std::cout << "--------- in draw_map_molecules() calling m.draw_map_molecule() imol " << ii
                      << " shader " << shader.name << std::endl;

         m.map_as_mesh_gl_lines_version.set_material(m.material_for_maps); // how/why is this needed? (seems that it is)
         m.draw_map_molecule(eye, draw_transparent_maps, shader, mvp, model_rotation, eye_position, ep,
                             lights, background_colour, perspective_projection_flag);
      }
   }

   if (needs_blend_reset) {
      glDisable(GL_BLEND);
   }
}

void
graphics_info_t::draw_model_molecules(stereo_eye_t eye) {

   // std::cout << "draw_model_molecules() --- start ---" << std::endl;

   // This is only called in "Plain" mode - i.e. it is not used in "Fancy" mode.
   // This function is called by draw_molecules(), which in turn is called by render_3d_scene()

   glm::mat4 mvp = get_molecule_mvp(eye);

   // std::cout << "debug:: mvp in draw_model_molecules() is     " << glm::to_string(mvp) << std::endl;
   glm::mat4 model_rotation = get_model_rotation();

   glm::vec4 bgc(background_colour, 1.0);

   for (int ii=n_molecules()-1; ii>=0; ii--) {
      if (! is_valid_model_molecule(ii)) continue;

      molecule_class_info_t &m = molecules[ii];
      // std::cout << "draw_model_mqolecules() A " << ii << " m.draw_it " << m.draw_it << std::endl;

      if (! m.draw_it) continue;

      // 20230827-PE this is for the new/consolidated api-based instanced meshes.
      //
      Shader &shader_instances_p = shader_for_instanced_objects;
      float opacity = 1.0f;
      bool gl_lines_mode = false;
      bool show_just_shadows = false;
      m.model_molecule_meshes.draw(&shader_for_meshes, &shader_instances_p, eye, mvp, model_rotation, lights, eye_position,
                                   opacity, bgc, gl_lines_mode, shader_do_depth_fog_flag, show_just_shadows);

      if (show_symmetry) {
         Shader &symm_shader_p = shader_for_symmetry_atoms_bond_lines;
         m.draw_symmetry(&symm_shader_p, mvp, model_rotation, lights, eye_position, bgc, shader_do_depth_fog_flag);
      }
   }

   // this block of code should be a member function of molecule_class_info_t
   // (20220208-PE it mostly is now - but still it should be a one-liner here in draw_model_molecules()

   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {

      molecule_class_info_t &m = graphics_info_t::molecules[ii];
      if (! graphics_info_t::is_valid_model_molecule(ii)) continue;
      if (! m.draw_it) continue;

      if (m.draw_model_molecule_as_lines) {
         float lw = m.get_bond_thickness(); // returns an int.
         m.model_molecule_meshes.draw_simple_bond_lines(&shader_for_symmetry_atoms_bond_lines, mvp, bgc, lw, shader_do_depth_fog_flag);
      } else {
#if 0 // the molecule_as_mesh is not filled at the moment, because the bond generation is now on the instanced path.
         bool show_just_shadows = false;
         bool wireframe_mode = false;
         float opacity = 1.0f;
         Shader *shader_p = &shader_for_meshes;
         m.molecule_as_mesh.draw(shader_p, mvp, model_rotation, lights, eye_position, opacity, bgc,
                                  wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);
#endif
      }
      m.draw_dots(&shader_for_rama_balls, mvp, model_rotation, lights, eye_position,
                     bgc, shader_do_depth_fog_flag);

      m.draw_ncs_ghosts(&shader_for_meshes, eye, mvp, model_rotation, lights, eye_position, bgc);

      glEnable(GL_BLEND);
      draw_molecule_atom_labels(m, eye, mvp, model_rotation);

   }
}

void
graphics_info_t::draw_model_molecules_symmetry_with_shadows(stereo_eye_t eye) {

   if (show_symmetry) {
      for (int ii=n_molecules()-1; ii>=0; ii--) {
         if (! is_valid_model_molecule(ii)) continue;
         molecule_class_info_t &m = molecules[ii];
         if (! m.draw_it) continue;
         Shader &symm_shader_p = shader_for_symmetry_atoms_bond_lines;
         glm::mat4 model_rotation = get_model_rotation();
         glm::vec4 bgc(background_colour, 1.0);
         glm::mat4 mvp = get_molecule_mvp(eye);
         m.draw_symmetry(&symm_shader_p, mvp, model_rotation, lights, eye_position, bgc, shader_do_depth_fog_flag);
      }
   }
}

void
graphics_info_t::draw_map_molecules_with_shadows() {

   // called by draw_molecules_with_shadows()
}

void
graphics_info_t::draw_molecule_atom_labels(molecule_class_info_t &m,
                                           stereo_eye_t eye,
                                           const glm::mat4 &mvp,
                                           const glm::mat4 &view_rotation) {

   // pass the glarea widget width and height.

   // put a triangle or square where the atom label should be, facing the camera
   // "billboarding"

   glm::vec4 label_colour(font_colour.red, font_colour.green, font_colour.blue, 1.0);

   if (false) { // test label

      // Put atom label test at 42, 9, 13
      glm::vec3 point(42, 9, 13);
      // point = glm::vec3(0,0,0);

      glm::vec4 projected_point_2 = mvp * glm::vec4(point, 1.0);
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
   int n_symm_atoms_to_label = m.labelled_symm_atom_index_list.size();

   // std::cout << "draw_molecule_atom_labels " << n_atoms_to_label << " " << n_symm_atoms_to_label << std::endl;

   if (n_atoms_to_label == 0 && n_symm_atoms_to_label == 0) {
   } else {
      m.draw_atom_labels(brief_atom_labels_flag, seg_ids_in_atom_labels_flag,
                         label_colour, eye, mvp, view_rotation);
   }

   // this is draw_generic_texts() - but not in its own function

   if (! generic_texts.empty()) {
      auto atom_label_colour = label_colour;
      for (unsigned int i=0; i<generic_texts.size(); i++) {
         const coot::generic_text_object_t gto = generic_texts[i];
         const std::string &label = gto.s;
         glm::vec3 position(gto.x, gto.y, gto.z);
         tmesh_for_labels.draw_atom_label(label, position, atom_label_colour,
                                          &shader_for_atom_labels, eye, mvp, view_rotation,
                                          glm::vec4(background_colour, 1.0),
                                          shader_do_depth_fog_flag,
                                          perspective_projection_flag);
      }
   }


   glDisable(GL_BLEND);

}

void
graphics_info_t::draw_intermediate_atoms(stereo_eye_t eye, unsigned int pass_type) { // draw_moving_atoms()

   // std::cout << "draw_intermediate_atoms() --- start --- " << std::endl;

   // ----------------------------------------
   // move this function into graphics-info-draw-model-molecules.cc
   // ----------------------------------------

   // this function gets called from draw_with_shadows() - but doesn't yet
   // use the shodows meshes shader.

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glm::mat4 mvp = get_molecule_mvp(eye);
   glm::mat4 model_rotation = get_model_rotation();

   molecule_class_info_t &m = graphics_info_t::moving_atoms_molecule;
   glm::vec4 bgc(background_colour, 1.0);
   bool show_just_shadows = false; // make this a member data item.

   float opacity = 1.0f;

   if (pass_type == PASS_TYPE_STANDARD) {

      // instanced:
      Shader &shader_p = shader_for_instanced_objects;
      // m.model_molecule_meshes.set_debug_mode(true);
      m.draw_molecule_as_meshes(&shader_p, mvp, model_rotation, lights, eye_position, bgc, shader_do_depth_fog_flag);
   }

   if (pass_type == PASS_TYPE_SSAO) {
      bool do_orthographic_projection = ! perspective_projection_flag;
      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
      int w = allocation.width;
      int h = allocation.height;
      auto model_matrix = get_model_matrix();
      auto view_matrix = get_view_matrix();
      auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
      m.model_molecule_meshes.draw_for_ssao(&shader_for_meshes_for_ssao,
                                            &shader_for_instanced_meshes_for_ssao,
                                            model_matrix, view_matrix, projection_matrix);
   }

   if (pass_type == PASS_TYPE_GEN_SHADOW_MAP) {

      // 20231011-PE have I used the right shader here?
      Shader &shader = shader_for_meshes_shadow_map;
      glm::vec3 dummy_eye_position;
      bool gl_lines_mode = false;
      glm::mat4 mvp_orthogonal = glm::mat4(1.0f); // placeholder
      unsigned int light_index = 0;
      std::map<unsigned int, lights_info_t>::const_iterator it = lights.find(light_index);
      if (it != lights.end()) {
         graphics_info_t g;
         const auto &light = it->second;
         mvp_orthogonal = g.get_mvp_for_shadow_map(light.direction);
      }
      bool do_depth_fog = true;
      glm::vec4 bg_col(background_colour, 1.0);
      m.model_molecule_meshes.draw(&shader_for_models, &shader_for_instanced_objects,
                                   eye, mvp_orthogonal, model_rotation, lights, dummy_eye_position,
                                   opacity, bg_col, gl_lines_mode, do_depth_fog, show_just_shadows);
   }

}


#include "Instanced-Markup-Mesh.hh"

void
graphics_info_t::setup_rama_balls() {

   rama_balls_mesh.setup_octasphere(2);
   rama_balls_mesh.setup_instancing_buffers(3000);
}

void
graphics_info_t::update_rama_balls(std::vector<Instanced_Markup_Mesh_attrib_t> *balls) {

   // so this function should be called set_rama_balls_new positions_and_colours

   // the calling function calls
   // rama_balls_mesh.update_instancing_buffers(balls) after this function

   const auto &rr = saved_dragged_refinement_results;

   balls->clear();

  glm::vec3 screen_up_dir(0.2, 0.3, 0.3);

   // std::cout << "update rama ball for " << rr.all_ramas.size() << " balls " << std::endl;
   for (unsigned int i=0; i<rr.all_ramas.size(); i++) {

      const coot::atom_spec_t &spec_CA = rr.all_ramas[i].atom_spec_CA;
      mmdb::Atom *at = spec_CA.get_atom(moving_atoms_asc->mol);
      if (at) {

         float d = rr.all_ramas[i].distortion;
         glm::vec3 atom_position(at->x, at->y, at->z);
         glm::vec3 ball_position(rr.all_ramas[i].ball_pos_x,
                                 rr.all_ramas[i].ball_pos_y,
                                 rr.all_ramas[i].ball_pos_z);
         float size = 0.38;
         // std::cout << "debug d " << d << std::endl;
         float ra = hud_geometry_distortion_to_rotation_amount_rama(d);
         coot::colour_t cc(0.1, 0.9, 0.2);
         cc.rotate(ra);
         glm::vec4 col = cc.to_glm();
         Instanced_Markup_Mesh_attrib_t ball(col, ball_position, size);
         // float d1 = d + 85.0; 20210902-PE Hmm.
         float d1 = d + 16.0;
         float d2 = - d1 * 0.4; // 20210902-PE was 0.016;
         // std::cout << "d2: " << d2 << std::endl;
         if (d2 < 0.0) d2 = 0.0;
         if (d2 > 1.0) d2 = 1.0;
         ball.specular_strength = 0.01 + d2;
         ball.shininess = 0.9 + 155.0 * d2;
         balls->push_back(ball);
      }
   }
}



void
graphics_info_t::draw_intermediate_atoms_rama_balls(stereo_eye_t eye, unsigned int pass_type) {

   // 20220302-PE Currently I don't draw rama balls with instancing
   //             It would be nice to have.
   return;

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   Shader &shader = graphics_info_t::shader_for_rama_balls;

   glm::mat4 mvp = get_molecule_mvp(eye);
   glm::vec3 eye_position = get_world_space_eye_position();
   glm::mat4 model_rotation = get_model_rotation();
   glm::vec4 bg_col(background_colour, 1.0);
   bool do_depth_fog = shader_do_depth_fog_flag;
   // this is a bit ugly

   // the balls are updated after a refinement cycle has finished -
   // no need to do it here
   // graphics_info_t g;
   // std::vector<Instanced_Markup_Mesh_attrib_t> balls;
   // update_rama_balls(&balls);
   // rama_balls_mesh.update_instancing_buffers(balls);

   // note: from graphics-info.h: static Instanced_Markup_Mesh rama_balls_mesh;
   rama_balls_mesh.draw(&shader, mvp, model_rotation, lights, eye_position, bg_col, do_depth_fog);

}

void
graphics_info_t::setup_atom_pull_restraints_glsl() {

   // build the triangles for the cylinder and
   // set m_VertexArray_for_pull_restraints_ID
   //     m_VertexBuffer_for_pull_restraints_ID
   //     n_indices_for_atom_pull_triangles

   n_atom_pulls = 0; // now a class variable.
   for (std::size_t i=0; i<atom_pulls.size(); i++) {
      const atom_pull_info_t &atom_pull = atom_pulls[i];
      if (atom_pull.get_status()) {
         std::pair<bool, int> spec = atom_pull.find_spec(moving_atoms_asc->atom_selection,
                                                         moving_atoms_asc->n_selected_atoms);
         if (spec.first)
            n_atom_pulls++;
      }
   }

   if (n_atom_pulls > 0) {
      unsigned int n_slices = 10;
      unsigned int n_stacks = 2;
      n_vertices_for_atom_pull_restraints  = n_slices * (n_stacks +1) * n_atom_pulls;
      n_triangles_for_atom_pull_restraints = n_slices * n_stacks * 2  * n_atom_pulls;
      // add in the triangles for the arrow-head (lots of vertices at the arrow tip because from cylinder)
      unsigned int n_stacks_for_arrow_tip = 6;
      n_vertices_for_atom_pull_restraints  += n_stacks_for_arrow_tip * n_slices * 2 * n_atom_pulls;
      n_triangles_for_atom_pull_restraints += n_stacks_for_arrow_tip * n_slices * 2 * n_atom_pulls;
      // the indices of the vertices in the triangles (3 indices per triangle)
      unsigned int *flat_indices = new unsigned int[n_triangles_for_atom_pull_restraints * 3];
      unsigned int *flat_indices_start = flat_indices;
      unsigned int ifi = 0; // index into flat indices - running
      coot::api::vertex_with_rotation_translation *vertices = new coot::api::vertex_with_rotation_translation[n_vertices_for_atom_pull_restraints];
      coot::api::vertex_with_rotation_translation *vertices_start = vertices;
      unsigned int iv = 0; // index into vertices - running

      // auto vertex_with_rotation_translation_to_generic_vertex = [] (const coot::api::vertex_with_rotation_translation &v) {
      // return vertex_with_rotation_translation(v.pos, v.normal, v.colour);
      // };

      for (std::size_t i=0; i<atom_pulls.size(); i++) {
         const atom_pull_info_t &atom_pull = atom_pulls[i];
         if (atom_pull.get_status()) {
            std::pair<bool, int> spec = atom_pull.find_spec(moving_atoms_asc->atom_selection,
                                                            moving_atoms_asc->n_selected_atoms);
            if (spec.first) {
               float arrow_head_length = 0.2;

               mmdb::Atom *at = moving_atoms_asc->atom_selection[spec.second];

               // coot::Cartesian pt_start_c(at->x, at->y, at->z);
               // coot::Cartesian pt_end_c(atom_pull.pos.x(), atom_pull.pos.y(), atom_pull.pos.z());
               // coot::Cartesian b = pt_end_c - pt_start_c;

               glm::vec3 pt_start_g(at->x, at->y, at->z);
               glm::vec3 pt_end_g(atom_pull.pos.x(), atom_pull.pos.y(), atom_pull.pos.z());
               glm::vec3 b = pt_end_g - pt_start_g;

               float bl_pull = glm::distance(b, glm::vec3(0,0,0));
               float bl = bl_pull - arrow_head_length;
               if (arrow_head_length > bl_pull)
                  arrow_head_length = bl_pull;
               glm::vec3 b_uv = glm::normalize(b);
               float bl_stick = bl;
               if (bl_stick < 0.0) bl_stick = 0.0;

               // coot::Cartesian meeting_point = pt_start_c + b_uv * bl_stick;
               // coot::CartesianPair pos_pair(pt_start_c, meeting_point);
               glm::vec3 meeting_point = pt_start_g + b_uv * bl_stick;
               std::pair<glm::vec3, glm::vec3> pos_pair(pt_start_g, meeting_point);
               float radius = 0.1;

               cylinder_with_rotation_translation c(pos_pair, radius, radius, bl, n_slices, n_stacks);
               for (std::size_t j=0; j<c.triangle_indices_vec.size(); j++) {
                  flat_indices[ifi  ] = c.triangle_indices_vec[j].point_id[0]+iv;
                  flat_indices[ifi+1] = c.triangle_indices_vec[j].point_id[1]+iv;
                  flat_indices[ifi+2] = c.triangle_indices_vec[j].point_id[2]+iv;
                  ifi += 3;
               }
               for (std::size_t j=0; j<c.vertices.size(); j++) {
                  // Use a constructor here when hmt code has been correctly integrated
                  vertices[iv] = c.vertices[j];
                  vertices[iv].colour = glm::vec4(0.8, 0.5, 0.3, 1.0);
                  iv++;
               }

               // coot::CartesianPair pp(pt_end_c, meeting_point);
               std::pair<glm::vec3, glm::vec3> pp(pt_end_g, meeting_point);
               // std::cout << "arrow-head: " << radius << " " << arrow_head_length << std::endl;
               cylinder_with_rotation_translation c_arrow_head(pp, 2.0 * radius, 0.0, arrow_head_length,
                                                               n_slices, n_stacks_for_arrow_tip);
               for (std::size_t j=0; j<c_arrow_head.triangle_indices_vec.size(); j++) {
                  if (ifi < n_triangles_for_atom_pull_restraints * 3) {
                     flat_indices[ifi  ] = c_arrow_head.triangle_indices_vec[j].point_id[0]+iv;
                     flat_indices[ifi+1] = c_arrow_head.triangle_indices_vec[j].point_id[1]+iv;
                     flat_indices[ifi+2] = c_arrow_head.triangle_indices_vec[j].point_id[2]+iv;
                  } else {
                     std::cout << "ERROR:: indexing for c_arrow_head "
                               << ifi << " " << n_triangles_for_atom_pull_restraints << std::endl;
                  }
                  ifi += 3;
               }
               for (std::size_t j=0; j<c_arrow_head.vertices.size(); j++) {
                  vertices[iv] = c_arrow_head.vertices[j];
                  vertices[iv].colour = glm::vec4(0.8,0.5,0.3,1.0);
                  iv++;
               }
            }
         }
      }

      // does this need to be done every time? I doubt it. Needs check.
      //
      glGenVertexArrays(1, &m_VertexArray_for_pull_restraints_ID);
      GLenum err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() A"
                          << " with GL err " << err << std::endl;
      glBindVertexArray(m_VertexArray_for_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() B"
                          << " with GL err " << err << std::endl;
      glGenBuffers(1, &m_VertexBuffer_for_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() C"
                          << " with GL err " << err << std::endl;
      glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer_for_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() D"
                          << " with GL err " << err << std::endl;
      GLuint n_bytes = sizeof(coot::api::vertex_with_rotation_translation) * n_vertices_for_atom_pull_restraints;
      // maybe STATIC_DRAW, maybe not
      glBufferData(GL_ARRAY_BUFFER, n_bytes, vertices, GL_DYNAMIC_DRAW);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() E"
                          << " with GL err " << err << std::endl;


      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);
      glEnableVertexAttribArray(2);
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17c\n";
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(0 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17d\n";
      glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(1 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17e\n";
      glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17f\n";

      // translate position, 3, size 3 floats
      glEnableVertexAttribArray(3);

      // surely this (annd below) has been set-up already? -- CheckMe.
      glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(3 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17aa\n";

      // positions, 4, size 3 floats
      glEnableVertexAttribArray(4);
      err = glGetError(); if (err) std::cout << "GL error bonds 6\n";
      glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(4 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 7\n";

      //  normals, 5, size 3 floats
      glEnableVertexAttribArray(5);
      err = glGetError(); if (err) std::cout << "GL error bonds 11\n";
      glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(5 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 12\n";

      //  colours, 6, size 4 floats
      glEnableVertexAttribArray(6);
      err = glGetError(); if (err) std::cout << "GL error bonds 16\n";
      glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation),
                            reinterpret_cast<void *>(6 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17\n";

      glGenBuffers(1, &m_IndexBuffer_for_atom_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() G"
                         << " with GL err " << err << std::endl;
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_atom_pull_restraints_ID);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() H"
                         << " with GL err " << err << std::endl;
      n_bytes = n_triangles_for_atom_pull_restraints * 3 * sizeof(unsigned int);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, flat_indices, GL_STATIC_DRAW);
      err = glGetError();
      if (err) std::cout << "setup_atom_pull_restraints_glsl() --end-- err " << err << std::endl;

      delete [] flat_indices;
      delete [] vertices;

   }
}


// static
void
graphics_info_t::draw_atom_pull_restraints(stereo_eye_t eye) {

   // Note to self: do this first with standard (modern) OpenGL.
   //
   // Then do it again with instances. It will be faster to draw bonds and atoms that way.
   // Maybe density lines too.

   // don't draw this if intermediate atoms are not shown
   //
   if (! regularize_object_bonds_box.empty()) {
      if (!moving_atoms_asc) return;
      if (moving_atoms_asc->n_selected_atoms > 0) {

         // std::cout << "drawing atom pull restraints with n_atom_pulls " << n_atom_pulls << std::endl;

         if (n_atom_pulls > 0) { // class variable now.

            if (false)
               std::cout << "drawing atom pull restraints with n_atom_pulls " << n_atom_pulls
                         << " and n_triangles_for_atom_pull_restraints " << n_triangles_for_atom_pull_restraints
                         << std::endl;

            Shader &shader = shader_for_models;
            shader.Use();
            GLuint err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glUseProgram() "
                               << err << std::endl;

            glBindVertexArray(m_VertexArray_for_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindVertexArray()"
                               << " with GL err " << err << std::endl;

            glm::mat4 mvp = get_molecule_mvp(eye);
            glm::mat4 model_rotation = get_model_rotation();
            GLuint mvp_location = shader.mvp_uniform_location;
            GLuint view_rotation_location = shader.view_rotation_uniform_location;
            glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
            glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &model_rotation[0][0]);

            std::map<unsigned int, lights_info_t>::const_iterator it;
            unsigned int light_idx = 0;
            it = lights.find(light_idx);
            if (it != lights.end())
               shader.setup_light(light_idx, it->second, model_rotation);
            light_idx = 1;
            it = lights.find(light_idx);
            if (it != lights.end())
               shader.setup_light(light_idx, it->second, model_rotation);

            glm::vec4 bg_col(background_colour, 1.0f);
            shader.set_vec4_for_uniform("background_colour", bg_col);
            shader.set_bool_for_uniform("do_depth_fog", shader_do_depth_fog_flag);

            glEnableVertexAttribArray(0);
            glEnableVertexAttribArray(1);
            glEnableVertexAttribArray(2);
            glEnableVertexAttribArray(3);
            glEnableVertexAttribArray(4);
            glEnableVertexAttribArray(5);
            glEnableVertexAttribArray(6);

            GLuint n_verts = 3 * n_triangles_for_atom_pull_restraints;
            err = glGetError();
            if (err) std::cout << "      error draw_atom_pull_restraints() pre-glDrawElements() "
                               << n_verts << " with GL err " << err << std::endl;
            glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
            err = glGetError();
            if (err) std::cout << "   error in draw_atom_pull_restraints() glDrawElements() n_verts: "
                               << n_verts << " with GL err " << err << std::endl;
         }
      }
   }
}

// static
void // draw_proportional_editing_neighbour_displacement_max_radius
graphics_info_t::draw_intermediate_atoms_pull_restraint_neighbour_displacement_max_radius_ring() {

   if (! regularize_object_bonds_box.empty()) {
      if (!moving_atoms_asc) return;
      if (moving_atoms_asc->n_selected_atoms > 0) {
         auto rcc = RotationCentre();
         glm::vec3 rc(rcc.x(), rcc.y(), rcc.z());
         rc = glm::vec3(0,0,0);
         glm::mat4 unit(1.0);
         glm::mat4 trans = glm::translate(unit, -rc);
         glm::mat4 view = get_view_matrix();
         int w = graphics_x_size;
         int h = graphics_y_size;
         bool ortho_flag = true;
         glm::mat4 proj = get_projection_matrix(ortho_flag, w, h);
         glm::mat4 mvp = proj * view * trans;
         glm::mat4 model_rotation = get_model_rotation();
         bool use_model_rotation = false;
         lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring.draw(&shader_for_lines,
                                                                                   rc, mvp,
                                                                                   model_rotation,
                                                                                   use_model_rotation);
      }
   }
}


void
graphics_info_t::draw_molecular_triangles() {

   // goodby innards
}

// static
void
graphics_info_t::draw_particles(stereo_eye_t eye) {

   if (curmudgeon_mode) return;

   if (! particles.empty()) {
      if (mesh_for_particles.have_instances()) {
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation = get_model_rotation();
         mesh_for_particles.draw_particles(&shader_for_particles, mvp, model_rotation);
      }
   }

   // std::cout << "debug:: draw_particles(): gone_diego_particles size " << meshed_particles_for_gone_diegos.size() << std::endl;
   if (! meshed_particles_for_gone_diegos.empty()) {
      for (unsigned int i=0; i<meshed_particles_for_gone_diegos.size(); i++) {
         Mesh &mesh(meshed_particles_for_gone_diegos[i].mesh);
         if (mesh.have_instances()) {
            glm::mat4 mvp = get_molecule_mvp(eye);
            glm::mat4 model_rotation = get_model_rotation();
            // std::cout << "debug:: draw_particles(): drawing gone diego particles! imesh: " << i << std::endl;
            mesh.draw_particles(&shader_for_particles, mvp, model_rotation);
         } else {
            std::cout << "draw_particles(): Ooops gone-diego imesh: " << i << " " << mesh.name << " has no instances" << std::endl;
         }
      }
   }

   { // gone difference map peaks.
      Mesh &mesh = meshed_particles_for_gone_diff_map_peaks.mesh;
      if (mesh.have_instances()) {
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation = get_model_rotation();
         mesh.draw_particles(&shader_for_particles, mvp, model_rotation);
      }
   }
}

// static
void
graphics_info_t::draw_happy_face_residue_markers(stereo_eye_t eye) {

   if (curmudgeon_mode) return;

   // make it work (somewhat) like particles, but we are using a screen-facing
   // texture, not a bespoke screen-facing n-triangle polygon.
   // So it's somewhat like the atom labels (a texture in 3D), in perspective
   // happy faces at the back are smaller.
   // But unlike labels, happy faces should be scaled by distance eye to rotation
   // centre, so that they always appear at a constant size (constant n-pixels wide
   // on the screen).

   if (tmesh_for_happy_face_residues_markers.draw_this_mesh) {

      if (tmesh_for_happy_face_residues_markers.have_instances()) {

         // the update of the instanced positions is done in the tick function

         graphics_info_t g; // needed for draw_count_max_for_happy_face_residue_markers. Use a better way?
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation = get_model_rotation();
         texture_for_happy_face_residue_marker.Bind(0);
         unsigned int draw_count = draw_count_for_happy_face_residue_markers;
         unsigned int draw_count_max = g.draw_count_max_for_happy_face_residue_markers;
         tmesh_for_happy_face_residues_markers.draw_fading_instances(&shader_for_happy_face_residue_markers,
                                                                     mvp, model_rotation,
                                                                     draw_count, draw_count_max);
      }
   }
}

// static
void
graphics_info_t::draw_anchored_atom_markers(stereo_eye_t eye) {

   if (tmesh_for_anchored_atom_markers.draw_this_mesh) {
      if (tmesh_for_anchored_atom_markers.have_instances()) {
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 view_rotation = get_model_rotation();
         glm::vec4 bg_col(background_colour, 1.0);
         texture_for_anchored_atom_markers.Bind(0);
         tmesh_for_anchored_atom_markers.draw_instances(&shader_for_happy_face_residue_markers,
                                                        mvp, view_rotation, bg_col, perspective_projection_flag);
      }
   }
}

void
graphics_info_t::draw_texture_meshes(stereo_eye_t eye) {

   // std::cout << "draw_texture_meshes() --- start --- " << std::endl;

   if (! texture_meshes.empty()) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      Shader &shader = shader_for_texture_meshes;
      for (unsigned int i=0; i<texture_meshes.size(); i++) {
         TextureMesh &tm = texture_meshes[i];
         // std::cout << "debug:: in draw_textures_meshes textures size " << tm.textures.size() << std::endl;
         if (! tm.textures.empty()) {
            // std::cout << "Binding and drawing the texture mesh" << std::endl;

            // std::cout << "............ get crow texture drawing code" << std::endl;

            bool do_depth_fog = true;
            int idx_start = tm.textures.size() - 1;
            for (int idx_texture=idx_start; idx_texture>=0; idx_texture--) {
               const auto &texture = tm.textures[idx_texture];
               if (false)
                  std::cout << "binding texture " << idx_texture << " to unit " << texture.unit << std::endl;
               tm.textures[idx_texture].texture.Bind(texture.unit);
            }
            glEnable(GL_BLEND);
            // we need some user control over the map section opacity
            tm.draw(&shader, mvp, model_rotation, lights, eye_position, bg_col, do_depth_fog);
            glDisable(GL_BLEND);

#if 0
            //
            // 20211018-PE it matters that these get called in the right order!
            //
            // BASE_TEXTURE before NORMAL_MAP
            // I've reversed the indexing of idx_texture for now - there's probably a better
            // way of sorting this out.
            // texture.unit for BASE_TEXTURE is 0
            // texture.unit for NORMAL_MAP is 1
            // as it should be (AFAICS)
            //
            // for (unsigned int idx_texture=0; idx_texture<tm.textures.size(); idx_texture++) {
            int idx_start = tm.textures.size() - 1;
            for (int idx_texture=idx_start; idx_texture>=0; idx_texture--) {
               const auto &texture = tm.textures[idx_texture];
               if (texture.texture_type == TextureInfoType::BASE_TEXTURE) {
                  // std::cout << "   binding BASE_TEXTURE texture_unit " << texture.unit << std::endl;
                  tm.textures[idx_texture].texture.Bind(texture.unit);
               }
               if (texture.texture_type == TextureInfoType::NORMAL_MAP) {
                  // std::cout << "   binding NORMAL_MAP texture_unit " << texture.unit << std::endl;
                  tm.textures[idx_texture].texture.Bind(texture.unit);
               }
               tm.draw(&shader, mvp, view_rotation, lights, eye_position, bg_col, do_depth_fog);
            }
#endif
         }
      }
   }
}

void
graphics_info_t::draw_hud_refinement_dialog_arrow_tab() {

   if (showing_intermediate_atoms_from_refinement()) {
      // show a (clickable/highlighting) HUD texture - indicating that refinement parameters dialog
      // can be shown (as an overlay)

      // std::cout << "here in draw_hud_refinement_dialog_arrow_tab() B " << std::endl;

      auto get_munged_offset_and_scale =  [] (HUDTextureMesh::screen_position_origins_t spo,
                                              const glm::vec2 &offset_natural,
                                              float scale_x_natural, float scale_y_natural,
                                              int glarea_width, int glarea_height) {

                                             glm::vec2 offset_rel = glm::vec2(0,0);

                                             // we don't need to be clever now that the shader is passed
                                             // the relative origin.
                                             // So this code may not be needed.

                                             float w = static_cast<float>(glarea_width);
                                             float h = static_cast<float>(glarea_height);

                                             float wr = static_cast<float>(900)/static_cast<float>(glarea_width);
                                             float hr = static_cast<float>(900)/static_cast<float>(glarea_height);

                                             if (spo == HUDTextureMesh::TOP_LEFT)
                                                offset_rel = glm::vec2(-1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr) - offset_natural;
                                             if (spo == HUDTextureMesh::BOTTOM_LEFT)
                                                offset_rel = glm::vec2(wr - 1.0, hr - 1.0) * offset_natural;
                                             if (spo == HUDTextureMesh::BOTTOM_RIGHT)
                                                offset_rel = glm::vec2(1.0 + offset_natural.x/wr, -1.0 + offset_natural.y/hr);

                                             if (spo == HUDTextureMesh::TOP_RIGHT) {
                                             }

                                             glm::vec2 scales_new(scale_x_natural * wr, scale_y_natural * hr);

                                             return std::pair<glm::vec2, glm::vec2>(offset_rel, scales_new);
                                          };

      glDisable(GL_DEPTH_TEST);
      if (hud_refinement_dialog_arrow_is_moused_over) {
         // std::cout << "hud_refinement_dialog_arrow_is_moused_over " << std::endl;
         texture_for_hud_refinement_dialog_arrow_highlighted.Bind(0);
      } else {
         // std::cout << "hud_refinement_dialog_arrow_is_moused_over not " << std::endl;
         texture_for_hud_refinement_dialog_arrow.Bind(0);
      }

      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
      int w = allocation.width;
      int h = allocation.height;
      float wf = static_cast<float>(w);

      tmesh_for_hud_refinement_dialog_arrow.set_scales(glm::vec2(0.04, 0.04));
      glm::vec2 position_natural(-0.04f, -0.1f); // relative to top right
      tmesh_for_hud_refinement_dialog_arrow.set_position(position_natural);
      auto p_s = get_munged_offset_and_scale(HUDTextureMesh::TOP_RIGHT, position_natural, 1.0, 1.0, w, h);
      glm::vec2 munged_position_offset = p_s.first;
      glm::vec2 munged_scales = p_s.second;
      tmesh_for_hud_refinement_dialog_arrow.set_window_resize_position_correction(munged_position_offset);
      tmesh_for_hud_refinement_dialog_arrow.set_window_resize_scales_correction(munged_scales);

      Shader &shader = shader_for_hud_image_texture;
      tmesh_for_hud_refinement_dialog_arrow.draw(&shader, HUDTextureMesh::TOP_RIGHT);
   }
}

void
graphics_info_t::draw_hud_colour_bar() {

   // 20230919-PE draw_hud_colour_bar_flag needs to be turned on somewhere when alphafold is used. It isn't.
   // Previously I had used user_defined_colours.empty(), which is the wrong test when we have
   // user defined colours and not alphafold (probably most of the time this is the case).
   //
   if (! draw_hud_colour_bar_flag) return;

   if (user_defined_colours.empty()) return;

   // this is the colour bar for Alphafold pLDDTs and the like

   // I think that all the draw_hud_*() functions should be passed h, w.
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   float aspect_ratio = static_cast<float>(w)/static_cast<float>(h);

   // ---------------- draw HUD colour texture ---------------------------------------

   auto get_munged_offset_and_scale =  [] (HUDTextureMesh::screen_position_origins_t spo,
                                           const glm::vec2 &offset_natural,
                                           float scale_x_natural, float scale_y_natural,
                                           int glarea_width, int glarea_height) {

      glm::vec2 offset_rel = glm::vec2(0,0);

      // we don't need to be clever now that the shader is passed
      // the relative origin.
      // So this code may not be needed.

      float w = static_cast<float>(glarea_width);
      float h = static_cast<float>(glarea_height);

      float wr = static_cast<float>(900)/static_cast<float>(glarea_width);
      float hr = static_cast<float>(900)/static_cast<float>(glarea_height);

      if (spo == HUDTextureMesh::TOP_LEFT)
         offset_rel = glm::vec2(-1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr) - offset_natural;
      if (spo == HUDTextureMesh::BOTTOM_LEFT)
         offset_rel = glm::vec2(wr - 1.0, hr - 1.0) * offset_natural;
      if (spo == HUDTextureMesh::BOTTOM_RIGHT)
         offset_rel = glm::vec2(1.0 + offset_natural.x/wr, -1.0 + offset_natural.y/hr);

      if (spo == HUDTextureMesh::TOP_RIGHT) {
      }

      glm::vec2 scales_new(scale_x_natural * wr, scale_y_natural * hr);

      return std::pair<glm::vec2, glm::vec2>(offset_rel, scales_new);
   };


   glDisable(GL_DEPTH_TEST);
   texture_for_hud_colour_bar.Bind(0);

   tmesh_for_hud_colour_bar.set_scales(glm::vec2(0.5, 0.02));

   glm::vec2 position_natural(-1.0f, -0.07f); // relative to top right // was -0.06
   tmesh_for_hud_colour_bar.set_position(position_natural);

   // I think that we need a HUDTextureMesh::TOP_MIDDLE position "target"
   //
   auto p_s = get_munged_offset_and_scale(HUDTextureMesh::TOP_RIGHT, position_natural, 1.0, 1.0, w, h);
   glm::vec2 munged_position_offset = p_s.first;
   glm::vec2 munged_scales = p_s.second;
   tmesh_for_hud_colour_bar.set_window_resize_position_correction(munged_position_offset);
   tmesh_for_hud_colour_bar.set_window_resize_scales_correction(munged_scales);

   Shader &shader = shader_for_hud_image_texture;
   tmesh_for_hud_colour_bar.draw(&shader, HUDTextureMesh::TOP_RIGHT);

   // ---------------- draw HUD colour text ---------------------------------------

   // label and tick marks

   // these text positions need to be relative to the colour bar on window-resize. They are not.

   std::vector<std::pair<std::string, glm::vec2> > positioned_texts;
   positioned_texts.push_back(std::make_pair("pLDDT", glm::vec2(-0.63, 0.896)));
   positioned_texts.push_back(std::make_pair(  "0.0", glm::vec2(-0.53, 0.95)));
   positioned_texts.push_back(std::make_pair( "25.0", glm::vec2(-0.28, 0.95)));
   positioned_texts.push_back(std::make_pair( "50.0", glm::vec2(-0.03, 0.95)));
   positioned_texts.push_back(std::make_pair( "75.0", glm::vec2( 0.22, 0.95)));
   positioned_texts.push_back(std::make_pair("100.0", glm::vec2( 0.45, 0.95)));

   for (unsigned int i=0; i<positioned_texts.size(); i++) {

      std::string label = positioned_texts[i].first;
      const auto &pos = positioned_texts[i].second;
      bool use_label_highlight = false;
      glm::vec2 label_scale(0.00008, 0.00008 * aspect_ratio);
      tmesh_for_hud_geometry_tooltip_label.set_scales(label_scale);
      tmesh_for_hud_geometry_tooltip_label.set_position(pos);
      tmesh_for_hud_geometry_tooltip_label.set_window_resize_position_correction(munged_position_offset);
      tmesh_for_hud_geometry_tooltip_label.set_window_resize_scales_correction(munged_scales);
      tmesh_for_hud_geometry_tooltip_label.draw_label(label, use_label_highlight,
                                                      &shader_for_hud_geometry_tooltip_text,
                                                      ft_characters);
   }

}



void
graphics_info_t::draw_molecules(stereo_eye_t eye) {

   // this is not called in fancy mode.

   // opaque things

   draw_outlined_active_residue(eye);

   draw_intermediate_atoms(eye, PASS_TYPE_STANDARD);

   draw_intermediate_atoms_rama_balls(eye, PASS_TYPE_STANDARD);

   draw_intermediate_atoms_pull_restraint_neighbour_displacement_max_radius_ring(); // proportional editing

   draw_atom_pull_restraints(eye);

   // return; // no draw

   draw_molecules_other_meshes(eye, PASS_TYPE_STANDARD);

   draw_instanced_meshes(eye);

   draw_map_molecules(eye, false); // transparency

   draw_unit_cells(eye);

   draw_environment_graphics_object(eye);

   draw_hydrogen_bonds_mesh(eye); // like boids

   draw_boids(eye);

   draw_particles(eye);

   draw_happy_face_residue_markers(eye);

   draw_bad_nbc_atom_pair_markers(eye, PASS_TYPE_STANDARD);

   draw_bad_nbc_atom_pair_dashed_lines(eye, PASS_TYPE_STANDARD);

   draw_chiral_volume_outlier_markers(eye, PASS_TYPE_STANDARD);

   draw_unhappy_atom_markers(eye, PASS_TYPE_STANDARD);

   draw_anchored_atom_markers(eye);

   // this is the last opaque thing to be drawn because the atom labels are blended.
   // It should be easy to break out the atom label code into its own function. That
   // might be better.
   //
   draw_model_molecules(eye);

   // transparent things...

   draw_map_molecules(eye, true); // transparent

   draw_generic_objects(PASS_TYPE_STANDARD, eye);

   // moved down
   draw_meshed_generic_display_object_meshes(eye, PASS_TYPE_STANDARD);

}

void
graphics_info_t::draw_molecules_with_shadows(stereo_eye_t eye) {

   int n_mols = n_molecules();
   bool show_just_shadows = false;
   glm::mat4 mvp = get_molecule_mvp(eye);
   auto model_rotation_matrix = get_model_rotation();

   int light_index = 0; // 20220215-PE for now.
   glm::mat4 light_view_mvp = get_light_space_mvp(light_index);
   glm::vec4 bg_col_v4(background_colour, 1.0f);

   // models

   for (int i=0; i<n_mols; i++) {
      if (is_valid_model_molecule(i)) {
         molecule_class_info_t &m = molecules[i];
         if (m.draw_it) {

            if (m.draw_model_molecule_as_lines) {

               float lw = m.get_bond_thickness(); // returns an int.
               m.model_molecule_meshes.draw_simple_bond_lines(&shader_for_symmetry_atoms_bond_lines, mvp, bg_col_v4, lw, shader_do_depth_fog_flag);

            } else {

               float opacity = 1.0;
               shader_for_instanced_meshes_with_shadows.Use();
#if 0 // before model molecules were instanced
               m.molecule_as_mesh.draw_with_shadows(&shader_for_meshes_with_shadows, mvp, model_rotation_matrix, lights,
                                                    eye_position, opacity, bg_col_v4, shader_do_depth_fog_flag, light_view_mvp,
                                                    shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
               m.draw_molecule_as_meshes_with_shadows(&shader_for_instanced_meshes_with_shadows, mvp, model_rotation_matrix, lights,
                                                      eye_position, opacity, bg_col_v4, shader_do_depth_fog_flag, light_view_mvp,
                                                      shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
#endif

               // shader_for_instanced_meshes_with_shadows.set_bool_for_uniform("do_fresnel", false); // models should not fresnel
               // 20230904-PE Let's try using the new model_molecule_as_meshes class
               m.model_molecule_meshes.draw_molecule_with_shadows(&shader_for_instanced_meshes_with_shadows, mvp, model_rotation_matrix, lights,
                                                                  eye_position, opacity, bg_col_v4, shader_do_depth_fog_flag, light_view_mvp,
                                                                  shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
            }

            // this has not been shadowified:
            m.draw_dots(&shader_for_rama_balls, mvp, model_rotation_matrix, lights, eye_position,
                        bg_col_v4, shader_do_depth_fog_flag);

            m.draw_ncs_ghosts(&shader_for_meshes, eye, mvp, model_rotation_matrix, lights, eye_position, bg_col_v4);

            glEnable(GL_BLEND);
            // good idea to not use shadows on atom labels?
            // 20220226-PE not here.
            // draw_molecule_atom_labels(m, mvp, model_rotation_matrix);
         }
      }
   }

   // maps


   // put this function inside molecule_class_info_t
   //
   auto draw_map_with_shadow = [] (molecule_class_info_t &m,
                                   const glm::mat4 &mvp,
                                   const glm::mat4 &model_rotation_matrix,
                                   const glm::vec4 &bg_col_v4,
                                   const glm::mat4 &light_view_mvp) {

                                  // draw_map_molecule() (the non-shadow version) is inside molecule_class_info_t - so should this function be too
                                  if (m.draw_it_for_map) {

                                     bool show_just_shadows = false;

                                     // this block should be inside draw_with_shadows(), indeed, fresnel_settings should be passed
                                     // to Mesh::draw_with_shadows()
                                     //
                                     shader_for_meshes_with_shadows.Use();
                                     shader_for_meshes_with_shadows.set_bool_for_uniform("do_fresnel",     m.fresnel_settings.state);
                                     shader_for_meshes_with_shadows.set_float_for_uniform("fresnel_bias",  m.fresnel_settings.bias);
                                     shader_for_meshes_with_shadows.set_float_for_uniform("fresnel_scale", m.fresnel_settings.scale);
                                     shader_for_meshes_with_shadows.set_float_for_uniform("fresnel_power", m.fresnel_settings.power);
                                     shader_for_meshes_with_shadows.set_vec4_for_uniform("fresnel_colour", m.fresnel_settings.colour);

                                     float opacity = m.density_surface_opacity;

                                     if (m.draw_it_for_map_standard_lines) {
                                        myglLineWidth(map_line_width);
                                        if (opacity < 1.0) m.map_as_mesh_gl_lines_version.use_blending = true;
                                        m.map_as_mesh_gl_lines_version.set_material(m.material_for_maps);
                                        m.map_as_mesh_gl_lines_version.draw_with_shadows(&shader_for_meshes_with_shadows, mvp, model_rotation_matrix, lights,
                                                                                         eye_position, opacity, bg_col_v4, shader_do_depth_fog_flag, light_view_mvp,
                                                                                         shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
                                     } else {
                                        if (opacity < 1.0) {
                                           m.map_as_mesh.use_blending = true;
                                           glm::vec3 eye_position_ws = get_world_space_eye_position();
                                           m.map_as_mesh.sort_map_triangles(eye_position_ws);
                                        }
                                        m.map_as_mesh.set_material(m.material_for_maps);
                                        m.map_as_mesh.draw_with_shadows(&shader_for_meshes_with_shadows, mvp, model_rotation_matrix, lights,
                                                                        eye_position, opacity, bg_col_v4, shader_do_depth_fog_flag, light_view_mvp,
                                                                        shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
                                     }
                                  }
                               };



   for (int i=0; i<n_mols; i++) {
      if (is_valid_map_molecule(i)) {
         molecule_class_info_t &m = molecules[i];
         draw_map_with_shadow(m, mvp, model_rotation_matrix, bg_col_v4, light_view_mvp);
      }
   }


   // convert these to read the shadow texture

   draw_model_molecules_symmetry_with_shadows(eye); // does symmetry

   draw_outlined_active_residue(eye);

   draw_intermediate_atoms(eye, PASS_TYPE_STANDARD);

   draw_intermediate_atoms_rama_balls(eye, PASS_TYPE_STANDARD);

   draw_atom_pull_restraints(eye);

   draw_meshed_generic_display_object_meshes(eye, PASS_TYPE_WITH_SHADOWS);

   draw_molecules_other_meshes(eye, PASS_TYPE_STANDARD);

   draw_instanced_meshes(eye);

   // draw_map_molecules(false); // transparency

   draw_unit_cells(eye);

   draw_environment_graphics_object(eye);

   draw_generic_objects(PASS_TYPE_STANDARD, eye);

   draw_hydrogen_bonds_mesh(eye); // like boids

   draw_anchored_atom_markers(eye);

   draw_boids(eye);

   draw_particles(eye);

   draw_happy_face_residue_markers(eye);

   // now drawn later, like atom labels
   // draw_bad_nbc_atom_pair_markers(PASS_TYPE_STANDARD);

   // this is the last opaque thing to be drawn because the atom labels are blended.
   // It should be easy to break out the atom label code into its own function. That
   // might be better.
   //
   // draw_model_molecules();

   // transparent things...

   // draw_map_molecules(true);


}


void
graphics_info_t::draw_molecules_atom_labels(stereo_eye_t eye) {

   // calls draw_molecule_atom_labels(m, mvp, model_rotation_matrix);

   int n_mols = n_molecules();
   glm::mat4 mvp = get_molecule_mvp(eye);
   auto model_rotation_matrix = get_model_rotation();

   for (int i=0; i<n_mols; i++) {
      if (is_valid_model_molecule(i)) {
         // glEnable(GL_BLEND); // surely this is called inside draw_labels()?
         molecule_class_info_t &m = molecules[i];
         if (m.draw_it) {
            draw_molecule_atom_labels(m, eye, mvp, model_rotation_matrix);
         }
      }
   }
}


// This does (draws) symmetry too.
//
// static
void
graphics_info_t::draw_environment_graphics_object(stereo_eye_t eye) {

#if 0   // old... keep for reference (for a while)
   graphics_info_t g;
   if (is_valid_model_molecule(mol_no_for_environment_distances)) {
      if (g.molecules[mol_no_for_environment_distances].is_displayed_p()) {
      g.environment_graphics_object_internal(environment_object_bonds_box);
      if (g.show_symmetry)
         g.environment_graphics_object_internal(symmetry_environment_object_bonds_box);
      }
   }
#endif

   if (is_valid_model_molecule(mol_no_for_environment_distances)) {
      molecule_class_info_t &m = molecules[mol_no_for_environment_distances];
      if (m.is_displayed_p()) {
         if (environment_show_distances) {
            glm::mat4 mvp = get_molecule_mvp(eye);
            glm::vec3 eye_position = get_world_space_eye_position();
            glm::mat4 model_rotation = get_model_rotation();
            glm::vec4 bg_col(background_colour, 1.0);

            bool do_depth_fog = shader_do_depth_fog_flag;
            bool show_just_shadows = false;
            bool wireframe_mode = false;
            float opacity = 1.0f;
            auto ccrc = RotationCentre();
            glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
            mesh_for_environment_distances.mesh.draw(&shader_for_moleculestotriangles,
                                                     eye, mvp, model_rotation,
                                                     lights, eye_position, rc, opacity, bg_col,
                                                     wireframe_mode, do_depth_fog, show_just_shadows);

            Shader *shader_p = &shader_for_atom_labels;

            GLenum err = glGetError();
            if (err) std::cout << "error draw_environment_graphics_object() before labela err "
                               << err << std::endl;

            if (! labels.empty()) {
               for (unsigned int i=0; i<labels.size(); i++) {
                  const std::string &label  = labels[i].label;
                  const glm::vec3 &position = labels[i].position;
                  const glm::vec4 &colour   = labels[i].colour;
                  // caches these textures in a map std::map<std::string, thing> where
                  // the key is the label.
                  tmesh_for_labels.draw_atom_label(label, position, colour, shader_p,
                                                   eye, mvp, model_rotation, bg_col, do_depth_fog,
                                                   perspective_projection_flag);
               }
            }

            if (show_symmetry) {

               // Fill me.

            }
         }
      }
   }
}

#include "molecular-mesh-generator.hh"
#include "pulse-data.hh"

void
graphics_info_t::update_mesh_for_outline_of_active_residue(int imol, const coot::atom_spec_t &spec, int n_press) {

   // if the range was not valid, but n_press was 5, we want some visual feedback - use
   // the invalid residue pulse (maybe it needs a new name then?). It uses delete_item_pulse_centres
   // for the centres of the pulses.
   auto setup_invalid_residue_pulse_for_invalid_range_outine = [] () {
      bool broken_line_mode = false;
      float radius_overall = 6.0;
      unsigned int n_rings = 3;
      lines_mesh_for_generic_pulse.setup_red_pulse(radius_overall, n_rings, broken_line_mode);
      pulse_data_t *pulse_data = new pulse_data_t(0, 12);
      gpointer user_data = reinterpret_cast<void *>(pulse_data);
      std::vector<glm::vec3> positions = { get_rotation_centre() };
      generic_pulse_centres = positions;
      gtk_widget_add_tick_callback(glareas[0], generic_pulse_function, user_data, NULL);
   };


   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      if (mol) {
         coot::residue_spec_t res_spec(spec);
         mmdb::Residue *residue_p = molecules[imol].get_residue(res_spec);
         if (residue_p) {
            attach_buffers();
            int bond_width = 9;
            int model_number = residue_p->GetModelNum();
            molecular_mesh_generator_t mmg;

            // setup the range, if the user had defined one in this molecule
            molecular_mesh_generator_t::range_t range;
            int imol_1 = in_range_first_picked_atom.int_user_data;
            int imol_2 = in_range_second_picked_atom.int_user_data;
            if (imol_1 == imol_2) {
               if (imol_1 == imol) {
                  coot::atom_spec_t spec_1 = in_range_first_picked_atom;
                  coot::atom_spec_t spec_2 = in_range_second_picked_atom;
                  if (spec_1.chain_id == spec_2.chain_id) {
                     int rn1 = spec_1.res_no;
                     int rn2 = spec_2.res_no;
                     if (rn2 < rn1) std::swap(rn1, rn2);
                     range = molecular_mesh_generator_t::range_t(spec_1.chain_id, rn1, rn2);
                  }
               }
            }

            std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > p =
               mmg.get_molecular_triangles_mesh_for_active_residue(imol, mol, model_number, residue_p, Geom_p(),
                                                                   bond_width, range, n_press);

            // if the range was not valid, but n_press was 5, we want some visual feedback - use
            // the invalid residue pulse (maybe it needs a new name then?) It uses delete_item_pulse_centres
            // for the centres of the pulses
            if (n_press == 5) {
               if (! range.is_valid) {
                  setup_invalid_residue_pulse_for_invalid_range_outine();
               }
            }

            mesh_for_outline_of_active_residue.clear();
            if (!p.first.empty()) {
               mesh_for_outline_of_active_residue.import(p);
               Material mat;
               mesh_for_outline_of_active_residue.setup(mat);
            }
         }
      }
   }
}


void
graphics_info_t::draw_outlined_active_residue(stereo_eye_t eye) {

   if (outline_for_active_residue_frame_count > 0) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      std::map<unsigned int, lights_info_t> dummy_lights;
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      Shader &shader = shader_for_outline_of_active_residue;
      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      auto ccrc = RotationCentre();
      glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
      mesh_for_outline_of_active_residue.draw(&shader, eye, mvp, model_rotation, dummy_lights, eye_position, rc,
                                              opacity, bg_col, wireframe_mode, false, show_just_shadows);
   }
};


void
graphics_info_t::draw_unit_cells(stereo_eye_t eye) {

   glm::mat4 mvp = get_molecule_mvp(eye);
   for (int ii=n_molecules()-1; ii>=0; ii--) {
      molecule_class_info_t &m = molecules[ii];
      m.draw_unit_cell(&shader_for_lines, mvp);
   }

}

void
graphics_info_t::draw_meshed_generic_display_object_meshes(stereo_eye_t eye, unsigned int pass_type) {

   // non-instanced.

   // std::cout << "draw_meshed_generic_display_object_meshes() with pass_type " << pass_type << std::endl;

   auto have_generic_display_objects_to_draw = [] () {
      bool generic_display_objects_to_draw = false;
      if (!generic_display_objects.empty()) {
         for (unsigned int i=0; i<generic_display_objects.size(); i++) {
            if (generic_display_objects[i].mesh.get_draw_this_mesh()) {
               generic_display_objects_to_draw = true;
               break;
            }
         }
      }
      return generic_display_objects_to_draw;
   };

   if (pass_type == PASS_TYPE_STANDARD) {

      glEnable(GL_BLEND); // 20250714-PE
      if (have_generic_display_objects_to_draw()) {
         glm::mat4 model_rotation = get_model_rotation();
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::vec4 bg_col(background_colour, 1.0);
         bool wireframe_mode = false;
         float opacity = 0.5;
         auto ccrc = RotationCentre();
         glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
         for (unsigned int i=0; i<generic_display_objects.size(); i++) {
            if (false)
               std::cout << "drawing i " << i << std::endl;
            generic_display_objects[i].mesh.draw(&shader_for_moleculestotriangles,
                                                 eye, mvp, model_rotation, lights, eye_position, rc, opacity,
                                                 bg_col, wireframe_mode, false, show_just_shadows);
         }
      }
   }

   if (pass_type == PASS_TYPE_SSAO) {
      if (have_generic_display_objects_to_draw()) {
         glm::vec4 bg_col(background_colour, 1.0);
         auto ccrc = RotationCentre();
         glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
         bool do_orthographic_projection = ! perspective_projection_flag;
         GtkAllocation allocation;
         gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
         int w = allocation.width;
         int h = allocation.height;
         auto model_matrix = get_model_matrix();
         auto view_matrix = get_view_matrix();
         auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
         for (unsigned int i=0; i<generic_display_objects.size(); i++) {
            generic_display_objects[i].mesh.draw_for_ssao(&shader_for_meshes_for_ssao,
                                                          model_matrix, view_matrix, projection_matrix);
         }
      }
   }

   if (pass_type == PASS_TYPE_GEN_SHADOW_MAP) {
      if (have_generic_display_objects_to_draw()) {
         int light_index = 0;
         std::map<unsigned int, lights_info_t>::const_iterator it;
         it = lights.find(light_index);
         if (it != lights.end()) {
            const auto &light = it->second;
            graphics_info_t g;
            glm::mat4 mvp_orthogonal = g.get_mvp_for_shadow_map(light.direction); // make this static?
            glm::mat4 model_rotation = get_model_rotation();
            glm::vec4 bg_col_v4(background_colour, 1.0f);
            auto ccrc = RotationCentre();
            glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
            glm::vec3 dummy_eye_position;
            float opacity = 1.0;
            bool do_depth_fog = false;
            bool gl_lines_mode = false;
            for (unsigned int i=0; i<generic_display_objects.size(); i++) {
               generic_display_objects[i].mesh.draw(&shader_for_meshes_shadow_map,
                                                    eye, mvp_orthogonal, model_rotation, lights, dummy_eye_position,
                                                    rc, opacity, bg_col_v4, gl_lines_mode,
                                                    do_depth_fog, show_just_shadows);
            }
         }
      }
   }

   if (pass_type == PASS_TYPE_WITH_SHADOWS) {

      // std::cout << "--------------------------- pass_type WITH SHADOWS!!!!!!!!!!!!!!!!!" << std::endl;
      if (have_generic_display_objects_to_draw()) {
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation = get_model_rotation();
         glm::vec4 bg_col_v4(background_colour, 1.0f);
         auto ccrc = RotationCentre();
         glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
         glm::vec3 eye_position;
         float opacity = 1.0;
         bool do_depth_fog = false;
         int light_index =  0;
         glm::mat4 light_view_mvp = get_light_space_mvp(light_index);
         bool show_just_shadows = false;

         for (unsigned int i=0; i<generic_display_objects.size(); i++) {
            generic_display_objects[i].mesh.draw_with_shadows(&shader_for_meshes_with_shadows,
                                                              mvp, model_rotation, lights, eye_position, opacity,
                                                              bg_col_v4, do_depth_fog, light_view_mvp,
                                                              shadow_depthMap_texture, shadow_strength, shadow_softness,
                                                              show_just_shadows);
         }
      }
   }
}



void
graphics_info_t::draw_molecules_other_meshes(stereo_eye_t eye, unsigned int pass_type) {

   // std::cout << "debug:: draw_molecules_other_meshes() ---start--- " << pass_type << std::endl;

   // This function doesn't draw these
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots H-bond
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots wide-contact
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots close-contact
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots small-overlap
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots vdw-surface
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots big-overlap

   bool draw_meshes = true;
   bool draw_mesh_normals = false;

   glm::vec3 eye_position = get_world_space_eye_position();
   glm::mat4 mvp = get_molecule_mvp(eye);
   glm::mat4 mvp_orthogonal = glm::mat4(1.0f); // placeholder
   glm::mat4 model_rotation = get_model_rotation();
   glm::vec4 bg_col(background_colour, 1.0);
   auto ccrc = RotationCentre();
   glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
   bool do_depth_fog = shader_do_depth_fog_flag;

   unsigned int light_index = 0;
   std::map<unsigned int, lights_info_t>::const_iterator it = lights.find(light_index);
   if (it != lights.end()) {
      graphics_info_t g;
      const auto &light = it->second;
      mvp_orthogonal = g.get_mvp_for_shadow_map(light.direction);
   }

   //std::cout << "mvp diag "
   // << mvp[0][0] << " " << mvp[1][1] << " " << mvp[2][2] << std::endl;

   glm::mat3 vrm(glm::toMat4(graphics_info_t::view_quaternion));

   if (draw_meshes) {
      bool have_meshes_to_draw = false;
      for (int i=n_molecules()-1; i>=0; i--) {
         if (! molecules[i].meshes.empty()) {
            have_meshes_to_draw = true;
            break;
         }
      }

      // std::cout << "in draw_meshed_generic_display_object_meshes() with have_meshes_to_draw " << have_meshes_to_draw << std::endl;

      if (have_meshes_to_draw) {

         glDisable(GL_BLEND);
         for (int ii=n_molecules()-1; ii>=0; ii--) {

            molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
            if (! is_valid_model_molecule(ii)) continue;
            for (unsigned int jj=0; jj<m.meshes.size(); jj++) {

               Mesh &mesh = m.meshes[jj];

               if (false)
                  std::cout << "debug:: mesh jj " << jj << " of " << m.meshes.size()
                            << " instanced: " << m.meshes[jj].is_instanced << std::endl;

               if (mesh.is_instanced) {

                  bool transferred_colour_is_instanced = false;
                  mesh.draw_instanced(pass_type,
                                      &shader_for_moleculestotriangles, eye, mvp,
                                      model_rotation, lights, eye_position,
                                      bg_col, do_depth_fog, transferred_colour_is_instanced);
               } else {

                  if (pass_type == PASS_TYPE_STANDARD) {
                     bool show_just_shadows = false;
                     bool wireframe_mode = false;
                     float opacity = 1.0f;
                     if (true)
                        m.meshes[jj].draw(&shader_for_moleculestotriangles, eye, mvp,
                                          model_rotation, lights, eye_position, rc, opacity, bg_col,
                                          wireframe_mode, do_depth_fog, show_just_shadows);

                  }
                  if (pass_type == PASS_TYPE_WITH_SHADOWS) {
                     if (false)
                        std::cout << "draw-molecule-other-meshes() " << m.meshes[jj].name << " "
                                  << shader_for_moleculestotriangles_with_shadows.name << std::endl;
                     bool show_just_shadows = false;
                     float opacity = 1.0f;
                     glm::mat4 light_view_mvp = get_light_space_mvp(light_index);
                     show_just_shadows = false;

                     // we don't do eye for shadows yet.
                     m.meshes[jj].draw_with_shadows(&shader_for_moleculestotriangles_with_shadows, mvp,
                                                    model_rotation, lights, eye_position, opacity, bg_col,
                                                    do_depth_fog, light_view_mvp,
                                                    shadow_depthMap_texture, shadow_strength,
                                                    shadow_softness, show_just_shadows);

                  }
                  if (pass_type == PASS_TYPE_SSAO) {
                     bool do_orthographic_projection = ! perspective_projection_flag;
                     GtkAllocation allocation;
                     gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
                     int w = allocation.width;
                     int h = allocation.height;
                     auto model_matrix = get_model_matrix();
                     auto view_matrix = get_view_matrix();
                     auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
                     m.meshes[jj].draw_for_ssao(&shader_for_meshes_for_ssao,
                                                model_matrix,
                                                view_matrix,
                                                projection_matrix);
                  }
                  if (pass_type == PASS_TYPE_GEN_SHADOW_MAP) { // i.e. generating the shadow map, not using it.

                     glm::vec3 dummy_eye_position;
                     bool gl_lines_mode = false;
                     bool show_just_shadows = false;
                     bool opacity = 1.0;

                     mesh.draw(&shader_for_meshes_shadow_map,
                               eye,
                               mvp_orthogonal,
                               model_rotation,
                               lights,
                               dummy_eye_position, rc,
                               opacity,
                               bg_col,
                               gl_lines_mode,
                               do_depth_fog,
                               show_just_shadows);

                  }
               }
            }
            glUseProgram(0);
         }

         if (draw_mesh_normals) {
            if (draw_normals_flag) {
               for (int ii=n_molecules()-1; ii>=0; ii--) {
                  molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
                  m.mesh_draw_normals(mvp);
               }
            }
         }
      }
   }
}

void
graphics_info_t::draw_instanced_meshes(stereo_eye_t eye) {

   // presumes opaque-only

   // ---------------------- Draw molecule instanced_meshes

   bool have_mol_meshes_to_draw = false;
   for (int i=n_molecules()-1; i>=0; i--) {
      if (! molecules[i].instanced_meshes.empty()) {
         if (! is_valid_model_molecule(i)) continue;
         if (molecules[i].draw_it) {
            have_mol_meshes_to_draw = true;
            break;
         }
      }
   }

   if (have_mol_meshes_to_draw) {
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool do_depth_fog = shader_do_depth_fog_flag;
      glDisable(GL_BLEND);
      for (int ii=n_molecules()-1; ii>=0; ii--) {
         if (! is_valid_model_molecule(ii)) continue;
         molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
         if (m.draw_it) {
            for (unsigned int jj=0; jj<m.instanced_meshes.size(); jj++) {
               // std::cout << "   graphics_info_t::draw_instanced_meshes() A " << m.instanced_meshes[jj].get_name() << std::endl;
               m.instanced_meshes[jj].draw(&shader_for_rama_balls, mvp,
                                           model_rotation, lights, eye_position, bg_col, do_depth_fog);
            }
         }
      }
   }

   // ---------------------- And draw our own

   if (! instanced_meshes.empty()) {
      glm::mat4 model_rotation = get_model_rotation();
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::vec4 bg_col(background_colour, 1.0);
      bool do_depth_fog = shader_do_depth_fog_flag;
      for (unsigned int jj=0; jj<instanced_meshes.size(); jj++) {
         // std::cout << "   graphics_info_t::draw_instanced_meshes() our own " << jj << " "
         // << instanced_meshes[jj].get_name() << std::endl;
         instanced_meshes[jj].draw(&shader_for_rama_balls, mvp,
                                   model_rotation, lights, eye_position, bg_col, do_depth_fog);
      }
   }
}

void
graphics_info_t::draw_meshes() {

   // I don't think that this function is called - which is good because draw_meshed_generic_display_object_meshes()
   // and draw_instanced_meshes() are called from elsewhere.

   // presumes only opaques

   // std::cout << "------------------- draw_meshes() " << std::endl;

   // draw_meshed_generic_display_object_meshes();
   // draw_instanced_meshes();
}

void
graphics_info_t::draw_cube(stereo_eye_t eye, GtkGLArea *glarea, unsigned int cube_type) {

   // std::cout << "draw_cube() with cube_type " << cube_type << std::endl;

   gtk_gl_area_make_current(glarea);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_central_cube() A0 err " << err << std::endl;

// wrap this if you use it again. myglLineWidth()
#ifdef __APPLE__
   // Shut up Mac OS X. This should not give an error
#else
   glLineWidth(2.0);  // GLv4 antialiasing - OpenGL implementations are not required to support this
                      // but they should not create an error. Mac does.
   err = glGetError();
   if (err) std::cout << "error draw_central_cube() A1 glLineWidth() err " << err << std::endl;
#endif


   // To see the possible values of the line width in aliased mode:
   // GLfloat line_width_max_min[2] = {0.0f, 0.0f};
   // glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, lineWidthRange);
   // This may not be possible in GL_LINE_SMOOTH mode.

   glm::mat4 mvp = get_molecule_mvp(eye);
   glm::mat4 model_rotation = get_model_rotation(); // hhmm... naming ... 20220212-PE  fixed now.

   glBindVertexArray(central_cube_vertexarray_id);
   err = glGetError(); if (err) std::cout << "   error::draw_central_cube() B err " << err << std::endl;
   glUseProgram(shader_for_central_cube.get_program_id());
   err = glGetError(); if (err) std::cout << "   error::draw_central_cube() C err " << err << std::endl;
   glm::vec3 rc = get_rotation_centre();
   if (cube_type == VIEW_CENTRAL_CUBE) {
      mvp = glm::translate(mvp, rc);
      float s = user_defined_rotation_centre_crosshairs_size_scale_factor;
      glm::vec3 sc(s,s,s);
      mvp = glm::scale(mvp, sc);
   }
   if (cube_type == ORIGIN_CUBE) {
      glm::vec3 sc(0.3f, 0.3f, 0.3f);
      mvp = glm::scale(mvp, sc);
   }

   // we don't diverge here on the cube type. Maybe change the name of the shader
   // because it does both
   Shader &shader = shader_for_central_cube;

   // we do this for all the shaders - Hmm.
   {
      GLuint mvp_location           = shader.mvp_uniform_location;
      GLuint view_rotation_location = shader.view_rotation_uniform_location;

      glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniformMatrix4fv() for mvp " << err << std::endl;
      glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &model_rotation[0][0]);
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniformMatrix4fv() for view_rotation " << err
                         << std::endl;

      GLuint line_colour_uniform_location = shader.line_colour_uniform_location;
      glm::vec4 lc(0.5, 0.4, 0.4, 1.0);
      if (cube_type == ORIGIN_CUBE)
         lc = glm::vec4(0.6, 0.6, 0.4, 1.0);
      glUniform4fv(line_colour_uniform_location, 1, glm::value_ptr(lc));
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniform4fv() for line colour " << err << std::endl;

      GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
      glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
      glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
      err = glGetError();
      if (err) std::cout << "error::draw_central_cube() glUniform4fv() for background " << err << std::endl;

   }

   glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "error::draw_central_cube() F glDrawElements() err " << err << std::endl;

   glBindVertexArray(0); // unbind
   glUseProgram(0);
}


void
graphics_info_t::draw_central_cube(stereo_eye_t eye, GtkGLArea *glarea) {
   draw_cube(eye, glarea, VIEW_CENTRAL_CUBE);
}

void
graphics_info_t::draw_origin_cube(stereo_eye_t eye, GtkGLArea *glarea) {

   draw_cube(eye, glarea, ORIGIN_CUBE);
}

void
graphics_info_t::draw_rotation_centre_crosshairs(stereo_eye_t eye, GtkGLArea *glarea, unsigned int pass_type) {

   // gtk_gl_area_make_current(glarea); // needed?, no it isn't.
   GLenum err = glGetError();
   if (err) std::cout << "error draw_rotation_centre_crosshairs() A0 err " << err << std::endl;

   glLineWidth(1.0);
   err = glGetError();
   if (err) std::cout << "error draw_rotation_centre_crosshairs() A1 err " << err << std::endl;

   glm::mat4 mvp = get_molecule_mvp(eye);
   glm::mat4 model_rotation = get_model_rotation();

   glBindVertexArray(rotation_centre_crosshairs_vertexarray_id);
   if (err) std::cout << "error draw_rotation_centre_crosshairs() B err " << err << std::endl;

   if (pass_type == PASS_TYPE_STANDARD) {
      shader_for_central_cube.Use(); // (it's drawing the crosshairs though - same shader)
   }

   if (pass_type == PASS_TYPE_SSAO) {
      shader_for_rotation_centre_cross_hairs_for_ssao.Use();
   }

   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   mvp = glm::translate(mvp, rc);
   // 20241105-PE is this a good idea?
   if (user_defined_rotation_centre_crosshairs_size_scale_factor < 0.02)
      user_defined_rotation_centre_crosshairs_size_scale_factor = 0.02;
   float s = 0.2f * user_defined_rotation_centre_crosshairs_size_scale_factor * zoom;
   glm::vec3 sc(s,s,s);
   mvp = glm::scale(mvp, sc);

   GLuint mvp_location           = shader_for_central_cube.mvp_uniform_location;
   GLuint view_rotation_location = shader_for_central_cube.view_rotation_uniform_location;

   glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error::draw_rotation_centre_crosshairs() glUniformMatrix4fv() for mvp " << err << std::endl;
   glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &model_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "error::draw_rotation_centre_crosshairs() glUniformMatrix4fv() for view_rotation " << err
                      << std::endl;

   if (pass_type == PASS_TYPE_STANDARD) {

      bool is_bb = graphics_info_t::background_is_black_p();
      glm::vec4 line_colour = rotation_centre_cross_hairs_colour;
      if (! is_bb)
         line_colour = glm::vec4(1.0f - rotation_centre_cross_hairs_colour[0],
                                 1.0f - rotation_centre_cross_hairs_colour[1],
                                 1.0f - rotation_centre_cross_hairs_colour[0],
                                 1.0f);

      GLuint line_colour_uniform_location = shader_for_central_cube.line_colour_uniform_location;
      glUniform4fv(line_colour_uniform_location, 1, glm::value_ptr(line_colour));

      GLuint background_colour_uniform_location = shader_for_central_cube.background_colour_uniform_location;
      glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
      glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));

      err = glGetError();
      if (err) std::cout << "error::draw_rotation_centre_crosshairs() glUniformMatrix4fv() for background " << err
                         << std::endl;
   }

   if (pass_type == PASS_TYPE_SSAO) {

      // need to be sent uniforms model, view. projection

      int w = graphics_x_size;
      int h = graphics_y_size;
      bool do_orthographic_projection = ! perspective_projection_flag;
      glm::mat4 model_matrix = get_model_matrix();
      glm::mat4  view_matrix = get_view_matrix();
      glm::mat4  proj_matrix = get_projection_matrix(do_orthographic_projection, w, h);

      shader_for_rotation_centre_cross_hairs_for_ssao.set_mat4_for_uniform("model",     model_matrix);
      shader_for_rotation_centre_cross_hairs_for_ssao.set_mat4_for_uniform("view",       view_matrix);
      shader_for_rotation_centre_cross_hairs_for_ssao.set_mat4_for_uniform("projection", proj_matrix);

   }

   glDrawElements(GL_LINES, 6, GL_UNSIGNED_INT, nullptr);
   if (err) std::cout << "error::draw_rotation_centre_crosshairs() glDrawElements " << err << std::endl;
   glBindVertexArray(0); // unbind
   glUseProgram(0);

}


// #include "event-controller-callbacks.hh"

void print_opengl_info() {

   // std::cout << "----------------------- print_opengl_info() ----------" << std::endl;

   const char *s1 = reinterpret_cast<const char *>(glGetString(GL_VERSION));
   const char *s2 = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
   const char *s3 = reinterpret_cast<const char *>(glGetString(GL_RENDERER));
   const char *s4 = reinterpret_cast<const char *>(glGetString(GL_VENDOR));
   if (s1 && s2 && s3 && s4) {
      std::string ss1(s1);
      std::string ss2(s2);
      std::string ss3(s3);
      std::string ss4(s4);
      // std::cout << "INFO:: GL Version:                  " << ss1 << std::endl;
      // std::cout << "INFO:: GL Shading Language Version: " << ss2 << std::endl;
      // std::cout << "INFO:: GL Renderer:                 " << ss3 << std::endl;
      // std::cout << "INFO:: GL Vendor:                   " << ss4 << std::endl;
      logger.log(log_t::INFO, "GL Version:",                  ss1);
      logger.log(log_t::INFO, "GL Shading Language Version:", ss2);
      logger.log(log_t::INFO, "GL Renderer:",                 ss3);
      logger.log(log_t::INFO, "GL Vendor:",                   ss4);
   } else {
      std::cout << "ERROR:: on_glarea_realize() null from glGetString()" << std::endl;
   }

}



gboolean
on_glarea_render(GtkGLArea *glarea) {

   // std::cout << "INFO:: in on_glarea_render() " << std::endl;

   bool screen_dump_frame_buffer = false;
   return graphics_info_t::render(screen_dump_frame_buffer);

}


void
on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   graphics_info_t g;

   // for the GL widget, not the window.
   g.graphics_x_size = width;
   g.graphics_y_size = height;

   std::cout << "INFO:: in on_glarea_resize() GtkGLArea widget dimensions " << width << " " << height << std::endl;

   // why do I need to do this?
   // setup_hud_text(width, height, g.shader_for_hud_text, false);
   // setup_hud_text(width, height, g.shader_for_atom_labels, true); // change the function name

   // g.setup_hud_geometry_bars(); // because they depend on the aspect ratio - but can't that be
                                   // passed as a uniform?

   // std::cout << "INFO:: Reset frame buffers " << width << "x" << height << std::endl;
   g.reset_frame_buffers(width, height);

   g.resize_framebuffers_textures_renderbuffers(width, height); // 20220131-PE added from crows merge

   g.reset_hud_buttons_size_and_position();
}


void
graphics_info_t::draw_measure_distance_and_angles(stereo_eye_t eye) {

   if (mesh_for_measure_distance_object_vec.get_draw_this_mesh()) {
      Shader &shader = shader_for_moleculestotriangles;
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation_matrix = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      auto ccrc = RotationCentre();
      glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
      mesh_for_measure_distance_object_vec.draw(&shader, eye, mvp, model_rotation_matrix, lights, eye_position, rc,
                                                opacity, bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

      mesh_for_measure_angle_object_vec.draw(&shader, eye, mvp, model_rotation_matrix, lights, eye_position, rc,
                                             opacity, bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

      if (! labels_for_measure_distances_and_angles.empty()) {
         for (unsigned int i=0; i<labels_for_measure_distances_and_angles.size(); i++) {
            const auto &label = labels_for_measure_distances_and_angles[i];
            tmesh_for_labels.draw_atom_label(label.label, label.position, label.colour, &shader_for_atom_labels,
                                             eye, mvp, model_rotation_matrix, bg_col,
                                             shader_do_depth_fog_flag, perspective_projection_flag);
         }
      }
   }

}


void
graphics_info_t::setup_lights() {

   lights_info_t light;
   light.position = glm::vec4(-2.0f, 2.0f, 5.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(0.5, 0.0, 1.0));
   // light.direction = glm::normalize(glm::vec3(0.0, 0.0, -1.0));
   light.diffuse *= 1.4;
   light.ambient *= 1.4;
   graphics_info_t::lights[0] = light;

   // dim(0.99) hardly changes anything (the argument is a multiplier)
   light.dim(0.25);
   light.position = glm::vec4(3.0f, -2.0f, 4.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(-1.0, 0.5, 1.0));
   // light.is_on = false;

   if (false)
      std::cout << "light 1 light: "
                << "ambient " << glm::to_string(light.ambient) << " "
                << "diffuse " << glm::to_string(light.diffuse) << std::endl;
   graphics_info_t::lights[1] = light;
}

// this is called from realize()
void
graphics_info_t::setup_hud_geometry_bars() {

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: setup_hud_geometry_bars() --start-- error " << err << std::endl;

   if (! glareas[0]) return;

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   float aspect_ratio = static_cast<float>(w)/static_cast<float>(h);

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: setup_hud_geometry_bars() A error " << err << std::endl;

   // gtk_gl_area_attach_buffers(gl_area); // needed? I think not (because we are in the base framebuffer when
                                           // this function is called) and it causes an glError to be set.

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: setup_hud_geometry_bars() B error " << err << std::endl;

   mesh_for_hud_geometry.setup_camera_facing_quad_for_bar();
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: setup_hud_geometry_bars() C error " << err << std::endl;

   mesh_for_hud_geometry.setup_instancing_buffer(500, sizeof(HUD_bar_attribs_t));

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: setup_hud_geometry_bars() C error " << err << std::endl;

   // If not found in this directory, then try default directory.
   texture_for_hud_geometry_labels_map["Rama"].init("hud-label-rama-small.png");
   // texture_for_hud_geometry_labels_map["Rama"].init("rama-plot-other-normal.png");
   texture_for_hud_geometry_labels_map["NBC"].init("hud-label-NBC-small.png");
   //texture_for_hud_geometry_labels_map["NBC"].init("rama-plot-other-normal.png");
   texture_for_hud_geometry_labels_map["Rota"].init("hud-label-rota-small.png");
   //texture_for_hud_geometry_labels_map["Rota"].init("rama-plot-other-normal.png");
   texture_for_hud_geometry_labels_map["Pull"].init("hud-label-pull-small.png");
   // texture_for_hud_geometry_labels_map["Pull"].init("rama-plot-other-normal.png");

   texture_for_hud_geometry_labels_map["Chiral"].init("hud-label-chiral-small.png");

   // texture_for_hud_tooltip_background.set_default_directory(coot::package_data_dir());
   texture_for_hud_tooltip_background.init("hud-tooltip.png"); // 94x47

   // Do I need to Use() the shader_for_hud_geometry_labels here?
   shader_for_hud_geometry_labels.Use();
   mesh_for_hud_geometry_labels.setup_quad();

   // glm::vec2 position(-0.98, 0.903);
   // glm::vec2 position(-0.0, 0.0);
   // glm::vec2 scales(0.56/aspect_ratio, 0.56);
   // mesh_for_hud_geometry_labels.set_position_and_scales(position, scales); // ""NBC, Pull"" texture

   // 20220319-PE  this is no longer drawn (at the moment)
   // mesh_for_hud_tooltip_background.setup_quad(); // does setup_buffers()
   // mesh_for_hud_tooltip_background.set_scales(glm::vec2(sc_x, sc_y));

   tmesh_for_hud_geometry_tooltip_label.setup_quad();
   glm::vec2 label_scale(0.000095, 0.000095/aspect_ratio);
   tmesh_for_hud_geometry_tooltip_label.set_scales(label_scale);

   // std::cout << "---------- done setup_hud_geometry_bars()" << std::endl;

}

void
graphics_info_t::setup_hud_buttons() {

   if (! glareas[0]) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: setup_hud_buttons() --start-- error " << err << std::endl;

   // std::cout << "debug:: in setup_hud_buttons() use_graphics_interface_flag " << use_graphics_interface_flag
   //          << " glareas[0] " << glareas[0] << std::endl;

   // attach_buffers(__FUNCTION__); // 20221005-PE not need, we are in the right framebuffer already.
                                    // And it causes a gl error to be set if it is called

   GError* error = gtk_gl_area_get_error(GTK_GL_AREA(glareas[0]));
   if (error)
      std::cout << "debug:: in setup_hud_buttons() current GError on glarea " << error->message << std::endl;

   err = glGetError();
   if (err) std::cout << "GL ERROR:: setup_hud_buttons() post attach_buffers() error " << err << std::endl;
   error = gtk_gl_area_get_error(GTK_GL_AREA(glareas[0]));
   if (error)
      std::cout << "debug:: in setup_hud_buttons() 2 current GError on glarea " << error->message << std::endl;

   mesh_for_hud_buttons.setup_vertices_and_triangles_for_button(); // instanced button

   unsigned int n_buttons_max = 20; // surely 6 is enough?
   mesh_for_hud_buttons.setup_instancing_buffer(n_buttons_max, sizeof(HUD_button_info_t));

   err = glGetError();
   if (err)
      std::cout << "debug:: in setup_hud_buttons() finish " << std::endl;

   // std::cout << "---------- done setup_hud_buttons()" << std::endl;
}

void
graphics_info_t::clear_hud_buttons() {

   attach_buffers();
   hud_button_info.clear();
   mesh_for_hud_buttons.update_instancing_buffer_data(hud_button_info); // empty
}

float
graphics_info_t::hud_geometry_distortion_to_bar_size_nbc(float distortion) {
   return distortion * 0.002;
}


float
graphics_info_t::hud_geometry_distortion_to_bar_size_atom_pull(float distortion) {
   return distortion * 0.0008;
}


float
graphics_info_t::hud_geometry_distortion_to_bar_size_rama(float distortion) {

#if 0 // 20210902-PE this is how it was
   float d1 = distortion + 200.0;
   float d2 = d1 * 0.0003;
   if (d2 < 0.0) d2 = 0.0;
   float d3 = 100.0 * d2 * d2;
   return d3;
#endif

#if 0 // 20220328-PE this is how it was
   float d1 = distortion + 16.0;
   float d2 = d1 / 6.0;
   if (d2 < 0.0) d2 = 0.0;
   float d3 = 0.1 * d2 * d2;
#endif

   float d1 = distortion + 18.0; // 18.0 is carefully chosen (16 too small, 22 too big)
   float d2 = d1 / 6.0f;
   if (d2 < 0.0f) d2 = 0.0f;
   float d3 = 0.055f * d2 * d2; // carefully chosen.
   if (d3 > 0.08) d3 = 0.08; // carefully chosen

   return d3;
}

// this function is used to colour the rama balla and colour the HUD geometry bars for rama
// (a good idea to use the same function, it turns out).
//
float
graphics_info_t::hud_geometry_distortion_to_rotation_amount_rama(float distortion) {

   // 20210902-PE note to self - the numbers coming here (distortion) need to be unscaled
   // by the rama restraints weight.
   // But ignoring that for now... Good values are less than -15 (all the way down to ~ -24)
   // Bad numbers are -9
   //
   // Final rotation amounts: 1.0 is pure green
   // 0.68 is pure red

   // When we don't have ramachandran restraints, then the rama balls are calculated
   // by make_generic_vertices_for_rama_balls() called by make_glsl_bonds_type_checked()
   // (if graphics_info_t::do_rama_restraints is false)
   //
   // Note also that cis peptides don't have ramachandran restraints.
   //
   float d2 = distortion + 16.0;
   float rotation_amount = 1.0 - 0.05 * d2; // 20220329-PE less "bad" colour than they used to be
   if (rotation_amount < 0.68) rotation_amount = 0.68; // red cap
   if (rotation_amount > 1.00) rotation_amount = 1.0;

   // std::cout << "debug:: distortion " << distortion << " rotation_amount " << rotation_amount << std::endl;

   return rotation_amount;
}

void
graphics_info_t::draw_hud_buttons() {

   // If you can't see them, then there are no buttons to draw.

   auto get_munged_offset_and_scale =  [] (HUDTextureMesh::screen_position_origins_t spo,
                                           const glm::vec2 &offset_natural,
                                           float scale_x_natural, float scale_y_natural,
                                           int glarea_width, int glarea_height) {

                                          glm::vec2 offset_rel = glm::vec2(0,0);

                                          // we don't need to be clever now that the shader is passed
                                          // the relative origin.
                                          // So this code may not be needed.

                                          float w = static_cast<float>(glarea_width);
                                          float h = static_cast<float>(glarea_height);

                                          float wr = static_cast<float>(900)/static_cast<float>(glarea_width);
                                          float hr = static_cast<float>(900)/static_cast<float>(glarea_height);

                                          if (spo == HUDTextureMesh::TOP_LEFT)
                                             offset_rel = glm::vec2(-1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr) - offset_natural;
                                          if (spo == HUDTextureMesh::BOTTOM_LEFT)
                                             offset_rel = glm::vec2(wr - 1.0, hr - 1.0) * offset_natural;
                                          if (spo == HUDTextureMesh::BOTTOM_RIGHT)
                                             offset_rel = glm::vec2(1.0 + offset_natural.x/wr, -1.0 + offset_natural.y/hr);

                                          if (spo == HUDTextureMesh::TOP_RIGHT) {
                                          }

                                          glm::vec2 scales_new(scale_x_natural * wr, scale_y_natural * hr);

                                          return std::pair<glm::vec2, glm::vec2>(offset_rel, scales_new);
                                       };

   if (hud_button_info.empty()) return;

   glEnable(GL_DEPTH_TEST); // needed?
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   float aspect_ratio = static_cast<float>(w)/static_cast<float>(h);

   float height_adjust = static_cast<float>(900)/static_cast<float>(h);
   float button_width  = HUD_button_info_t::button_width  * static_cast<float>(900)/static_cast<float>(w);
   float button_height = HUD_button_info_t::button_height * static_cast<float>(900)/static_cast<float>(h);

   glm::vec2 position_natural(-0.1f, .1f); // relative to bottom right
   auto p_s = get_munged_offset_and_scale(HUDTextureMesh::BOTTOM_RIGHT, position_natural, 1.0, 1.0, w, h);
   glm::vec2 munged_position_offset = p_s.first;
   glm::vec2 munged_scales = p_s.second;
   // std::cout << "window corrections scales: " << glm::to_string(munged_scales) << " pos-off: "
   //           << glm::to_string(munged_position_offset) << std::endl;
   mesh_for_hud_buttons.set_window_resize_scales_correction(munged_scales);
   mesh_for_hud_buttons.set_window_resize_position_correction(munged_position_offset);

   mesh_for_hud_buttons.draw(&shader_for_hud_buttons); // we have added the button instances before now.
                                                       // (actually in show_accept_reject_hud_buttons()).

   // do the texture for the labels all on the fly - is that sound?
   //
   glm::vec4 text_colour_white(0.95f, 0.95f, 0.95f, 1.0f);
   Shader &shader = shader_for_hud_geometry_tooltip_text;
   shader.Use();
   for (unsigned int i=0; i<hud_button_info.size(); i++) {
      const auto &button = hud_button_info[i];
      const std::string &label = button.button_label;
      if (! label.empty()) {
         std::string mesh_name = "HUDTexturemesh for button with label " + label;
         HUDTextureMesh htm(mesh_name);
         htm.setup_quad();
         float text_scale_raw = 0.4 * 0.00018;
         text_scale_raw *= 1.2; // 20220324-PE
         // text_scale_raw *= 100.0;
         float text_scale = text_scale_raw * height_adjust;
         glm::vec2 label_scale(text_scale / aspect_ratio, text_scale);
         htm.set_scales(label_scale);
         unsigned int n_chars = label.size(); // really I want the sum of x_advance for the letters. Can I get that?
         float x_advance = htm.get_sum_x_advance(label, ft_characters);
         float width_adjust = static_cast<float>(900)/static_cast<float>(w);
         float tl_adjust = - static_cast<float>(n_chars-1) * text_scale_raw * 2.2 * 50.0 * width_adjust;
         glm::vec2 pos = button.position_offset;
         pos += glm::vec2(0.0, 0.3 * button_height); // vertical adjustment for label
         pos += glm::vec2(0.5 * button_width, 0.00); // horizontal adjustment for label (lefttext is middle of button)
         pos += glm::vec2(tl_adjust, 0.00); // horizontal adjustment for text length
         htm.set_position(pos);
         htm.draw_label(label, text_colour_white, &shader, ft_characters);
      }
   }
}

void
graphics_info_t::setup_draw_for_translation_gizmo() {

   GLenum err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "--start--", stringify_error_code(err));

   attach_buffers(); // this causes a GL error - why?

   err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "A", stringify_error_code(err));

   size_t s = translation_gizmo.mesh.vertices.size();
   std::vector<s_generic_vertex> cv(s); //  conveted vertices
   for (unsigned int i=0; i<cv.size(); i++) {
      cv[i].pos    = translation_gizmo.mesh.vertices[i].pos;
      cv[i].normal = translation_gizmo.mesh.vertices[i].normal;
      cv[i].color  = translation_gizmo.mesh.vertices[i].color;
   }
   err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "B", stringify_error_code(err));
   translation_gizmo_mesh.clear(); // so that we don't add to the mesh!
   err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "C", stringify_error_code(err));
   translation_gizmo_mesh.import(cv, translation_gizmo.mesh.triangles);
   err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "D", stringify_error_code(err));
   translation_gizmo_mesh.setup_buffers();
   err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "E",  stringify_error_code(err));
   translation_gizmo_mesh.set_draw_this_mesh(false);
   err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("setup_draw_for_translation_gizmo"),
                 "F", stringify_error_code(err));

}


void
graphics_info_t::clear_gl_rama_plot() {

   gl_rama_plot.clear();
}

void
graphics_info_t::draw_hud_ramachandran_plot() {

   // 20240719-PE we do both the setup and draw here! I am not sure that is the right way.
   //             It might be vey slow. Ah, no, I was a bit clever, there is a position hash.

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   if (draw_gl_ramachandran_plot_flag) {
      if (draw_gl_ramachandran_plot_user_control_flag) {
         if (moving_atoms_asc) {
            if (moving_atoms_asc->n_selected_atoms > 0) {
               std::string residue_selection = "//";
               gl_rama_plot_t::draw_mode_t draw_mode = gl_rama_plot_t::draw_mode_t::DRAW_MODE;
               gl_rama_plot.setup_from(imol_moving_atoms, moving_atoms_asc->mol, residue_selection, draw_mode); // checks to see if an update is acutally needed.
               // no context switch needed for the HUD Rama plot
               bool clear_needed_flag = false;
               gl_rama_plot.draw(&shader_for_rama_plot_axes_and_ticks,
                                 &shader_for_rama_plot_phi_phis_markers, // instanced
                                 &shader_for_hud_image_texture, w, h, w, h,
                                 clear_needed_flag); // background texture (not text!), uses window_resize_position_correction
            }
         }
      }
   }

}

void
graphics_info_t::draw_hud_fps() {

   // these are "relative to"
   enum screen_position_origins_t { TOP_LEFT, TOP_RIGHT, BOTTOM_LEFT, BOTTOM_RIGHT};

   auto get_munged_offset_and_scale = [] (screen_position_origins_t spo,
                                          const glm::vec2 &offset_natural,
                                          float scale_x_natural, float scale_y_natural) {

                                         glm::vec2 offset_new = offset_natural;

                                         // glm::vec2 scales_new(scale_x_natural, scale_y_natural);

                                         // calculating the aspect_ratio like this takes 0.22 microseconds

                                         GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
                                         GtkAllocation allocation;
                                         gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
                                         int w = allocation.width;
                                         int h = allocation.height;

                                         // for 900 pixels and offset if 0.1 is 90 pixels.
                                         // 90 pixels in a 1000 pixels widths is 0.1/wr

                                         float wr = static_cast<float>(w)/static_cast<float>(900);
                                         float hr = static_cast<float>(h)/static_cast<float>(900);

                                         if (spo == TOP_LEFT)
                                            offset_new = glm::vec2(-1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr);
                                         if (spo == BOTTOM_LEFT)
                                            offset_new = glm::vec2(-1.0 + offset_natural.x/wr, -1.0 + offset_natural.y/hr);
                                         if (spo == BOTTOM_RIGHT)
                                            offset_new = glm::vec2(1.0 + offset_natural.x/wr, -1.0 + offset_natural.y/hr);
                                         if (spo == TOP_RIGHT)
                                            offset_new = glm::vec2(1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr);

                                         glm::vec2 scales_new(scale_x_natural/wr, scale_y_natural/hr);

                                         return std::pair<glm::vec2, glm::vec2>(offset_new, scales_new);
                                      };

   if (GetFPSFlag()) {

      // ----------------- HUD FPS string ---------------------------------

      std::string s = "FPS: " + coot::util::float_to_string_using_dec_pl(fps, 2);
      if (fps > 0) {
         float ms_per_frame = 1000.0 / fps;
         s += "  " + coot::util::float_to_string_using_dec_pl(ms_per_frame, 2) + " ms/frame";
      }

      if (fps_std_dev >= 0.0) {
         s += "  std.dev.: ";
         s += coot::util::float_to_string_using_dec_pl(fps_std_dev, 2);
         s += " ms/frame";
      }
      // attach_buffers(); // needed?
      HUDTextureMesh htm("mesh for FPS");
      htm.setup_quad(); // oops! Does this use the right framebuffer?
      Shader &shader = shader_for_hud_geometry_tooltip_text;  // change the name of this - it's for general (real) HUD text
      glm::vec4 col(0.7, 0.7, 0.4, 1.0);
      glm::vec4 grey(0.5, 0.5, 0.5, 0.4);
      glm::vec4 full_grey(0.5, 0.5, 0.5, 1.0);
      auto p_s = get_munged_offset_and_scale(TOP_LEFT, glm::vec2(0.1, -0.1), 0.0001, 0.0001);
      const glm::vec2 &munged_position_offset = p_s.first;
      const glm::vec2 &munged_scales = p_s.second;
      htm.set_scales(munged_scales);
      htm.set_position(munged_position_offset);
      htm.draw_label(s, col, &shader, ft_characters);

      // ----------------- HUD graph (in ms/frame) ---------------------------------

      if (frame_time_history_list.size() > 2) {

         myglLineWidth(1.0); // even Apple should have no problem with this.

         std::vector<glm::vec2> data;
         data.reserve(frame_time_history_list.size()+2); // it would be better if this was outside the hot path

         // base line
         float x_o = munged_position_offset.x;
         float y_o = munged_position_offset.y - 0.3;

         //LinesMesh lines_mesh; // 3D! (because I don't have a HUDLines class)

         // now convert those data to 3D vertices indices to be used by lines_mesh...
         // (we'll just use a unit matrix for the mvp when drawing them)
         std::vector<s_generic_vertex> vertices;
         vertices.reserve(data.size() + 6);
         std::vector<unsigned int> indices;
         glm::vec3 norm(0,0,1);  // not used

         // make glm::vec2 data and then convert that to OpenGL screen coordinates
         //
         float ms_to_opengl_y = fps_times_scale_factor; // default 0.0025
         // ms_to_opengl_y = 0.001;
         unsigned int time_count = 0;
         std::list<std::chrono::time_point<std::chrono::high_resolution_clock> >::const_iterator it;
         for (it = frame_time_history_list.begin(); it != frame_time_history_list.end(); ++it) {
            if (it != frame_time_history_list.begin()) {
               float x = static_cast<float>(time_count);
               const std::chrono::time_point<std::chrono::high_resolution_clock> &tp_this = *it;
               const std::chrono::time_point<std::chrono::high_resolution_clock> &tp_prev = *std::prev(it);
               auto delta_t = std::chrono::duration_cast<std::chrono::milliseconds>(tp_this - tp_prev).count();
               data.push_back(glm::vec2(x_o + 0.001 * x, y_o + ms_to_opengl_y * delta_t));
               time_count++;
            }
         }

         // base line and grid lines into vertices first
         //
         float y_tick_mark = 20.0 * ms_to_opengl_y; // 20ms converted to OpenGL y coord
         vertices.push_back(s_generic_vertex(glm::vec3(x_o,       y_o,                   -1), norm, full_grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o + 0.5, y_o,                   -1), norm, full_grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o,       y_o + y_tick_mark,     -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o + 0.5, y_o + y_tick_mark,     -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o,       y_o + 2 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o + 0.5, y_o + 2 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o,       y_o + 3 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o + 0.5, y_o + 3 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o,       y_o + 4 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o + 0.5, y_o + 4 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o,       y_o + 5 * y_tick_mark, -1), norm, grey));
         vertices.push_back(s_generic_vertex(glm::vec3(x_o + 0.5, y_o + 5 * y_tick_mark, -1), norm, grey));

         for (unsigned int i=0; i<data.size(); i++)
            vertices.push_back(s_generic_vertex(glm::vec3(data[i], -1), norm, col));

         for (unsigned int i=0; i<(data.size()-2+10); i++) {
            if (i == 1 || i == 3 || i == 5 || i == 7 || i == 9 || i == 11) {
               // no line betwween base line and grid lines and start of real data
            } else {
               indices.push_back(i);
               indices.push_back(i+1);
            }
         }

         // this looks like it can, from time to time, draw to the wrong framebuffer. Hmm.

         lines_mesh_for_hud_lines.update_vertices_and_indices(vertices, indices);
         glm::mat4 dummy_mat4(1.0);
         lines_mesh_for_hud_lines.draw(&shader_for_hud_lines, dummy_mat4, dummy_mat4);
      }
   }
}

void
graphics_info_t::show_atom_pull_toolbar_buttons() {

#if 0 // this is old-school graphics/gui, isn't it?
   if (use_graphics_interface_flag) {
      GtkWidget *button_1 = get_widget_from_builder("clear_atom_pull_restraints_toolbutton");
      GtkWidget *button_2 = get_widget_from_builder("auto_clear_atom_pull_restraints_togglebutton");

      if (button_1)
         gtk_widget_set_visible(button_1, TRUE);
      else
         std::cout << "in show_atom_pull_toolbar_buttons() missing button1" << std::endl;
      if (button_2)
         gtk_widget_set_visible(button_2, TRUE);
      else
         std::cout << "in show_atom_pull_toolbar_buttons() missing button2" << std::endl;
   }
#endif
}


void
graphics_info_t::hide_atom_pull_toolbar_buttons() {

   if (use_graphics_interface_flag) {
      GtkWidget *button_1 = get_widget_from_builder("clear_atom_pull_restraints_toolbutton");
      GtkWidget *button_2 = get_widget_from_builder("auto_clear_atom_pull_restraints_togglebutton");

      if (button_1)
         gtk_widget_set_visible(button_1, FALSE);
      if (button_2)
         gtk_widget_set_visible(button_2, FALSE);
   }
}

void
graphics_info_t::show_accept_reject_hud_buttons() {


   if (false)
      std::cout << "--------------------- show_accept_reject_hud_buttons() " << std::endl;

   // add some HUD buttons

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   HUD_button_info_t button_1("Accept");
   HUD_button_info_t button_2("Reject");
   HUD_button_info_t button_3("Sidechain 180"); // failure to lookup glyph for degree symbol ° :-(
   HUD_button_info_t button_4("Pep-Flip This");
   HUD_button_info_t button_5("Pep-Flip Next");
   HUD_button_info_t button_6("Backrub Rotamer");
   HUD_button_info_t button_7("JED Flip");
   HUD_button_info_t button_8("Cis/Trans");

   button_1.set_colour(glm::vec4(0.4, 0.7, 0.4, 0.5));
   button_2.set_colour(glm::vec4(0.7, 0.4, 0.4, 0.5));

   button_1.set_scales_and_position_offset(0, w, h);
   button_2.set_scales_and_position_offset(1, w, h);
   button_3.set_scales_and_position_offset(2, w, h);
   button_7.set_scales_and_position_offset(3, w, h);
   button_8.set_scales_and_position_offset(4, w, h);
   button_6.set_scales_and_position_offset(5, w, h);
   button_5.set_scales_and_position_offset(6, w, h);
   button_4.set_scales_and_position_offset(7, w, h);

   auto button_1_func = [] () {

                           graphics_info_t g;
                           g.stop_refinement_internal();
                           g.accept_moving_atoms();
                           // g.hud_button_info.clear();
                           g.hide_atom_pull_toolbar_buttons();
                           // g.draw_bad_nbc_atom_pair_markers = false;
                           g.clear_gl_rama_plot();
                           g.graphics_draw();
                           return true;
                   };

   auto button_2_func = [] () {
                           graphics_info_t g;
                           g.stop_refinement_internal();
                           g.clear_up_moving_atoms();
                           // g.hud_button_info.clear();
                           // g.draw_bad_nbc_atom_pair_markers = false;
                           g.rebond_molecule_corresponding_to_moving_atoms();
                           g.graphics_draw();
                           g.hide_atom_pull_toolbar_buttons();
                           g.clear_gl_rama_plot();
                           return true;
                        };
   auto button_3_func = [] () {
                           graphics_info_t g;
                           g.side_chain_flip_180_intermediate_atoms();
                           return true;
                        };

   auto button_4_func = [] () {
                           graphics_info_t g;
                           g.pepflip_intermediate_atoms();
                           return true;
                        };

   auto button_5_func = [] () {
                           graphics_info_t g;
                           g.pepflip_intermediate_atoms_other_peptide();
                           return true;
                        };

   auto button_6_func = [] () {
                           graphics_info_t g;
                           return g.backrub_rotamer_intermediate_atoms();
                        };

   auto button_7_func = [] () {
                           graphics_info_t g;
                           g.jed_flip_intermediate_atoms(false);
                           return true;
                        };

   auto button_8_func = [] () {
                           graphics_info_t g;
                           g.cis_trans_conversion_intermediate_atoms();
                           return true;
                        };

   button_1.connect(button_1_func);
   button_2.connect(button_2_func);
   button_3.connect(button_3_func);
   button_4.connect(button_4_func);
   button_5.connect(button_5_func);
   button_6.connect(button_6_func);
   button_7.connect(button_7_func);
   button_8.connect(button_8_func);

   hud_button_info.push_back(button_1);
   hud_button_info.push_back(button_2);
   hud_button_info.push_back(button_3);
   hud_button_info.push_back(button_4);
   hud_button_info.push_back(button_5);
   hud_button_info.push_back(button_6);
   hud_button_info.push_back(button_7);
   hud_button_info.push_back(button_8);

   gtk_gl_area_attach_buffers(gl_area);
   mesh_for_hud_buttons.update_instancing_buffer_data(hud_button_info);

}

void
graphics_info_t::reset_hud_buttons_size_and_position() {

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   for (unsigned int i=0; i<hud_button_info.size(); i++) {
      auto &button = hud_button_info[i];
      button.set_scales_and_position_offset(i, w, h);
   }
}

// static
float
graphics_info_t::get_x_base_for_hud_geometry_bars() {

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;

   float w_adjust = static_cast<float>(w)/static_cast<float>(900);

   // shift to more negative x when the window is wider
   return -0.83 - 0.02 * w_adjust;

}


void
graphics_info_t::draw_hud_geometry_bars() {

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glEnable(GL_DEPTH_TEST); // needed?
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   GLenum err = glGetError(); if (err) std::cout << "GL ERROR:: draw_hud_geometry_bars() A error " << err << std::endl;

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   coot::refinement_results_t &rr = saved_dragged_refinement_results;

   // --------------------- first draw the text (labels) texture -----------------------

   class hud_label_info_t {
   public:
      std::string name;
      unsigned int bar_index;
      float label_relative_width;
      hud_label_info_t(const std::string &n, unsigned int i, float rw) : name(n), bar_index(i), label_relative_width(rw) {}
      hud_label_info_t(const std::string &n, unsigned int i) : name(n), bar_index(i) { label_relative_width = 1.0; }
   };
   std::vector<hud_label_info_t> hud_label_info;
   hud_label_info.push_back(hud_label_info_t("Pull",   0, 0.7));
   hud_label_info.push_back(hud_label_info_t("Rama",   2, 1.0));
   hud_label_info.push_back(hud_label_info_t("Rota",   3, 0.9));
   hud_label_info.push_back(hud_label_info_t("NBC",    1, 0.9));
   hud_label_info.push_back(hud_label_info_t("Chiral", 4, 1.0));
   float x_base = get_x_base_for_hud_geometry_bars();

   // Note that the x-positions are are not the left-most edge of the label (hmm)

   // Don't forget these are *images* not actual text.

   for (const auto &hud_label : hud_label_info) {
      texture_for_hud_geometry_labels_map[hud_label.name].Bind(0);
      unsigned int bar_index = hud_label.bar_index;
      float text_y_offset = 0.017; // relative to the the bars in add_bars()
      float y = 0.943 + text_y_offset - 0.05 * static_cast<float>(bar_index); // c.f. add_bars()
      float width_adjust = static_cast<float>(900)/static_cast<float>(w);
      glm::vec2 scales(0.046 * hud_label.label_relative_width * width_adjust, 0.015);
      glm::vec2 position(x_base - 0.05 * width_adjust, y);
      mesh_for_hud_geometry_labels.set_position_and_scales(position, scales);
      mesh_for_hud_geometry_labels.draw(&shader_for_hud_geometry_labels);
      err = glGetError();
      if (err) std::cout << "GL ERROR:: draw_hud_geometry_bars() error Textures "
                         << hud_label.name << " " << err << std::endl;
   }

      // ----------------------- now draw the bars -----------------------

   auto probability_to_rotation_amount = [] (float probability) {
                                            // high probability should have low rotation
                                            float q_1 = 0.01 * (100.0 - probability);
                                            // if (q < 0) q = 0;
                                            // if (q > 1) q = 1;
                                            float q_2 = 0.68 * q_1;
                                            return q_2;
                                         };

   auto distortion_to_rotation_amount_nbc = [] (float distortion) {
                                               // we want to rotate to red (which is negative direction) but
                                               // rotate() doesn't work with negative rotations, so make it
                                               // 1.0 - amount (1.0 being a full rotation).
                                               float fac = 0.02; // 20220503-PE was 0.012; - I want the colours to be less green now
                                               float rotation_amount = 1.0 - fac * distortion;
                                               if (rotation_amount < 0.68) rotation_amount = 0.68;
                                               return rotation_amount;
                                            };

   auto chiral_volume_distortion_to_rotation_amount = [] (float distortion) {
                                               float fac = 0.05;
                                               float rotation_amount = 1.0 - fac * distortion;
                                               if (rotation_amount < 0.68) rotation_amount = 0.68;
                                               return rotation_amount;
   };


   // Other hud_geometry_distortion_to_x function are class members
   // So make this function match the function of the same name in check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked()
   //
   auto hud_geometry_distortion_to_bar_size_chiral = +[] (float distortion) {
      return distortion * 0.01f; // the f is needed
   };

   auto add_bars = [] (const std::vector<std::pair<coot::atom_spec_t, float> > &baddies,
                       unsigned int bar_index,
                       std::vector<HUD_bar_attribs_t> *new_bars_p,
                       const coot::residue_spec_t &active_residue_spec,
                       float x_base_for_hud_geometry_bars,
                       float (*distortion_to_rotation_amount)(float),
                       float (*distortion_to_bar_size)(float)) {

                      // to_top_left() needs to be the same as check_bars()
                      glm::vec2 to_top_left(x_base_for_hud_geometry_bars, 0.943 - 0.05 * static_cast<float>(bar_index));
                      float sum_l = 0;
                      int n = baddies.size();
                      glm::vec4 col_white(0.8,0.8, 0.8, 0.7);
                      for (int i=(n-1); i>=0; i--) {
                         coot::colour_t cc(0.1, 0.9, 0.2);
                         float d = baddies[i].second;
                         float rotation_amount = distortion_to_rotation_amount(d);
                         cc.rotate(rotation_amount);
                         glm::vec4 col = cc.to_glm();
                         col.w = 0.7;
                         glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                         float bar_length = distortion_to_bar_size(d);
                         bool this_atom_is_in_a_moving_atoms_residue = baddies[i].first.int_user_data;

                         // the vector of HUD_bar_attribs_t is fed directly to a opengl buffer.
                         // So I need "expand" to 2 bars right here - one with a "thin bar" attribute
                         //
                         if (! this_atom_is_in_a_moving_atoms_residue) {
                            float bar_height = 0.03; // universal
                            float bar_slither_y_scale = 0.3; // looks nice
                            float y_offset_main = bar_height * bar_slither_y_scale;
                            glm::vec2 position_offset_for_main = position_offset + glm::vec2(0, y_offset_main);
                            HUD_bar_attribs_t bar_main(col, position_offset_for_main, bar_length);
                            bar_main.scale_y = 1.0 - bar_slither_y_scale;
                            new_bars_p->push_back(bar_main);
                            // slither bar
                            glm::vec2 position_offset_for_slither = position_offset + glm::vec2(0,0);
                            HUD_bar_attribs_t bar_slither(col_white, position_offset_for_slither, bar_length);
                            bar_slither.scale_y = bar_slither_y_scale;
                            new_bars_p->push_back(bar_slither);
                         } else {
                            HUD_bar_attribs_t bar(col, position_offset, bar_length);
                            new_bars_p->push_back(bar);
                         }

                         // active residue purple bar (overdraws)
                         coot::residue_spec_t residue_for_bar(baddies[i].first);
                         if (residue_for_bar == active_residue_spec) {
                            glm::vec4 col_pink(1.0f, 0.2f, 1.0f, 0.8f);
                            HUD_bar_attribs_t bar(col_pink, position_offset + glm::vec2(0.0f, -0.01f), bar_length);
                            // new_bars_p->push_back(bar);
                            new_bars_p->insert(new_bars_p->begin(), bar);
                         }

                         sum_l += bar_length + 0.005; // with a gap between bars
                      }
                   };

   auto rota_sorter = [] (const rotamer_markup_container_t &rmc_1,
                          const rotamer_markup_container_t &rmc_2) {
                         if (rmc_2.rpi.probability < rmc_1.rpi.probability)
                            return true;
                         else
                            return false;
                      };

   auto add_rotamer_bars = [rota_sorter] (std::vector<HUD_bar_attribs_t> *new_bars_p,
                                          unsigned int bar_index,
                                          const coot::residue_spec_t &active_residue_spec,
                                          float x_base_for_hud_geometry_bars,
                                          rotamer_markup_container_t *rotamer_markups,
                                          int n_rotamer_markups) {

                              // std::cout << "debug:: add_rotamer_bars() n_rotamer_markups: " << n_rotamer_markups
                              // << " " << std::endl;
                              // this code has to be the same as the check_if_hud_bar_clicked code

                              // needs to be consitent with above and check_bars()
                              glm::vec2 to_top_left(x_base_for_hud_geometry_bars, 0.943 - 0.05 * static_cast<float>(bar_index));
                              glm::vec4 col_white(0.8,0.8, 0.8, 0.7);
                              std::vector<rotamer_markup_container_t> v;
                              // filter out the goodies
                              for (int i=0; i<n_rotamer_markups; i++)
                                 if (rotamer_markups[i].rpi.probability < 40) // 40 %
                                    if (rotamer_markups[i].rpi.probability >= 0)
                                       v.push_back(rotamer_markups[i]);
                              // sort the baddies
                              std::sort(v.begin(), v.end(), rota_sorter);
                              unsigned int n_rota_max = 20;
                              if (v.size() > n_rota_max) {
                                 unsigned int n_for_deletion = v.size() - n_rota_max;
                                 std::vector<rotamer_markup_container_t>::iterator v_begin = v.begin();
                                 std::vector<rotamer_markup_container_t>::iterator v_last;
                                 v_last = v_begin + n_for_deletion; // (line length)
                                 v.erase(v_begin, v_last);
                              }

                              float sum_l = 0;
                              for (unsigned int i=0; i<v.size(); i++) {

                                 rotamer_markup_container_t &rm = v[i];
                                 bool this_atom_is_in_a_moving_atoms_residue = rm.spec.int_user_data;

                                 float pr = rm.rpi.probability;
                                 float q = 0.01 * (48.0f - v[i].rpi.probability);
                                 if (q > 1.0) q = 1.0;
                                 if (q < 0.0) q = 0.0;
                                 float bar_length = std::pow(q, 6.0) * 4.0;

                                 if (false)
                                    std::cout << "bar i " << i << " " << v[i].spec << " " << v[i].col
                                              << " pr " << pr << " length " << bar_length << std::endl;

                                 const coot::colour_holder &ch = v[i].col;
                                 glm::vec4 col(ch.red, ch.green, ch.blue, 0.7);

                                 if (! this_atom_is_in_a_moving_atoms_residue) {
                                    float bar_height = 0.03; // universal
                                    float bar_slither_y_scale = 0.3; // looks nice
                                    float y_offset_main = bar_height * bar_slither_y_scale;
                                    glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                                    glm::vec2 position_offset_for_main = position_offset + glm::vec2(0, y_offset_main);
                                    HUD_bar_attribs_t bar_main(col, position_offset_for_main, bar_length);
                                    bar_main.scale_y = 1.0 - bar_slither_y_scale;
                                    new_bars_p->push_back(bar_main);
                                    // slither bar
                                    glm::vec2 position_offset_for_slither = position_offset + glm::vec2(0,0);
                                    HUD_bar_attribs_t bar_slither(col_white, position_offset_for_slither, bar_length);
                                    bar_slither.scale_y = bar_slither_y_scale;
                                    new_bars_p->push_back(bar_slither);


                                 } else {
                                    glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                                    HUD_bar_attribs_t bar(col, position_offset, bar_length);
                                    new_bars_p->push_back(bar);
                                 }

                                 // active residue purple bar (overdraws)
                                 coot::residue_spec_t residue_for_bar(v[i].spec);
                                 if (residue_for_bar == active_residue_spec) {
                                    glm::vec4 col_pink(1.0f, 0.2f, 1.0f, 0.8f);
                                    glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0) + glm::vec2(0.0f, -0.01f);
                                    HUD_bar_attribs_t bar(col_pink, position_offset, bar_length);
                                    // new_bars_p->push_back(bar);
                                    new_bars_p->insert(new_bars_p->begin(), bar);
                                 }

                                 sum_l += bar_length + 0.005; // with a gap between bars
                              }
                           };

   std::vector<HUD_bar_attribs_t> new_bars;

   // set the residue spec for the moving molecule.
   coot::Cartesian rc = get_rotation_centre_cart();
   float within_radius_limit = 9.0;
   coot::residue_spec_t moving_atoms_active_residue;
   mmdb::Atom *mv_at = get_moving_atoms_active_atom(rc, within_radius_limit);
   if (mv_at) {
      coot::atom_spec_t mv_at_spec(mv_at);
      moving_atoms_active_residue = coot::residue_spec_t(mv_at_spec);
   }

   float x_base_for_hud_geometry_bars = get_x_base_for_hud_geometry_bars();
   // add to new_bars
   add_bars(rr.sorted_atom_pulls, 0, &new_bars,
            moving_atoms_active_residue, x_base_for_hud_geometry_bars,
            distortion_to_rotation_amount_nbc, hud_geometry_distortion_to_bar_size_atom_pull);


   if (rr.refinement_results_contain_overall_nbc_score) {

      // 20220503-PE add_bars() take the argument std::vector<std::pair<coot::atom_spec_t, float> > &baddies
      // but now rr.sorted_nbc_baddies is std::vector<refinement_results_nbc_baddie_t>
      // so now I need to convert

      std::vector<std::pair<coot::atom_spec_t, float> > converted_baddies(rr.sorted_nbc_baddies.size() * 2); // 20230519-PE both ways
      for (unsigned int i=0; i<rr.sorted_nbc_baddies.size(); i++) {
         const auto &bip = rr.sorted_nbc_baddies[i];
         std::pair<coot::atom_spec_t, float> p_1(bip.atom_spec_1, bip.score);
         std::pair<coot::atom_spec_t, float> p_2(bip.atom_spec_2, bip.score);
         // 20230813-PE fixes a crash, I hope.
         if ((2*i+1) < converted_baddies.size()) {
            converted_baddies[2*i  ] = p_1;
            converted_baddies[2*i+1] = p_2;
         } else {
            std::cout << "ERROR:: out of range in converted_baddies  " << 2*i << " " << converted_baddies.size() << std::endl;
         }
      }
      add_bars(converted_baddies, 1, &new_bars, moving_atoms_active_residue,
               x_base_for_hud_geometry_bars, distortion_to_rotation_amount_nbc,
               hud_geometry_distortion_to_bar_size_nbc);
   }

   if (rr.refinement_results_contain_overall_rama_plot_score) {
      // std::cout << "add_bars() for rama with " << rr.sorted_rama_baddies.size() << " sorted baddies" << std::endl;
      add_bars(rr.sorted_rama_baddies, 2, &new_bars, moving_atoms_active_residue, x_base_for_hud_geometry_bars,
               hud_geometry_distortion_to_rotation_amount_rama, hud_geometry_distortion_to_bar_size_rama);
   }

   if (! rr.sorted_chiral_volume_baddies.empty()) {
      std::vector<std::pair<coot::atom_spec_t, float> > converted_baddies(rr.sorted_chiral_volume_baddies.size());
      for (unsigned int ii=0; ii<rr.sorted_chiral_volume_baddies.size(); ii++) {
         converted_baddies[ii] = std::make_pair(rr.sorted_chiral_volume_baddies[ii].atom_spec,
                                                rr.sorted_chiral_volume_baddies[ii].distortion);
      }
      add_bars(converted_baddies, 4, &new_bars, moving_atoms_active_residue, x_base_for_hud_geometry_bars,
               chiral_volume_distortion_to_rotation_amount, hud_geometry_distortion_to_bar_size_chiral);
   }

   // add rotas to new_bars

   // note to self - for rotation/translation, it's regularize_object_bonds_box that gets updated.
   // moving_atoms_molecule is not used until make_moving_atoms_graphics_object().
   // execute_rotate_translate_ready() calls make_moving_atoms_graphics_object()

   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         int nrms = moving_atoms_molecule.bonds_box.n_rotamer_markups;
         if (false) {
            std::cout << "Here are the rotamer markups: " <<  nrms << std::endl;
            for (int ii=0; ii<nrms; ii++) {
               auto &rm = moving_atoms_molecule.bonds_box.rotamer_markups[ii];
               std::cout << "   " << rm.spec << " " << rm.pos.format() << " " << rm.col << std::endl;
            }
         }
         if (nrms > 0) {
            add_rotamer_bars(&new_bars, 3, moving_atoms_active_residue,
                             x_base_for_hud_geometry_bars,
                             moving_atoms_molecule.bonds_box.rotamer_markups, nrms);
         }
      }
   }

   // std::cout << "in draw_hud_geometry_bars() " << new_bars.size() << std::endl;
   if (! new_bars.empty()) {
      // std::cout << "new bar size " << new_bars.size() << std::endl;
      mesh_for_hud_geometry.update_instancing_buffer_data(new_bars);
      mesh_for_hud_geometry.draw(&shader_for_hud_geometry_bars);
   }
   glDisable(GL_BLEND);

}

std::pair<bool, mmdb::Atom *>
graphics_info_t::check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(double mouse_x, double mouse_y, bool act_on_hit) {

   std::pair<bool, mmdb::Atom *> status_pair(false, 0);
   if (! moving_atoms_asc) return std::pair<bool, mmdb::Atom *>(false, 0);
   if (! moving_atoms_asc->mol) return std::pair<bool, mmdb::Atom *>(false, 0);

   coot::refinement_results_t &rr = saved_dragged_refinement_results;

   // this values in this loop must match those in the loop above
   // (draw_hud_geometry_bars())

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   double frac_x = mouse_x/static_cast<double>(w);
   double frac_y = 1.0 - mouse_y/static_cast<double>(h);
   glm::vec2 mouse_in_opengl_coords(2.0 * frac_x - 1.0, 2.0 * frac_y - 1.0);

   // these functions are copies of those in the above function. If you edit them again, make them
   // member functions.

   auto rota_sorter = [] (const rotamer_markup_container_t &rmc_1,
                          const rotamer_markup_container_t &rmc_2) {
                         if (rmc_2.rpi.probability < rmc_1.rpi.probability)
                            return true;
                         else
                            return false;
                      };

   // the act_on_hit flag is check to see if the move should be made, or that we merely return
   // a success status (we want to act when clicked, but return a status when moused-over)
   //
   auto check_blocks = [mouse_in_opengl_coords] (const std::vector<std::pair<coot::atom_spec_t, float> > &baddies,
                                                 unsigned int bar_index,
                                                 float x_base_for_hud_geometry_bars,
                                                 float (*distortion_to_bar_size)(float),
                                                 bool act_on_hit) {

                          bool status = false;
                          mmdb::Atom *at_out = 0;
                          glm::vec2 to_top_left(x_base_for_hud_geometry_bars, 0.943 - 0.05 * static_cast<float>(bar_index));
                          float sum_l = 0;
                          int n = baddies.size();
                          for (int i=(n-1); i>=0; i--) {
                             float d = baddies[i].second;
                             glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                             float bar_length = distortion_to_bar_size(d);
                             sum_l += bar_length + 0.005; // with a gap between bars

                             if (false) {
                                glm::vec2 position_offset_far_point = position_offset;
                                position_offset_far_point.x += bar_length;
                                std::cout << "checking "
                                          << glm::to_string(position_offset) << " "
                                          << glm::to_string(position_offset_far_point) << " "
                                          << " " << bar_length
                                          << " vs mouse " << glm::to_string(mouse_in_opengl_coords)
                                          << std::endl;
                             }

                             if (mouse_in_opengl_coords.x >= position_offset.x) {
                                if (mouse_in_opengl_coords.x <= (position_offset.x + bar_length)) {

                                   // std::cout << ":::::::::: x hit bar_index " << bar_index
                                   // << " i " << i << " " << baddies[i].first << std::endl;

                                   float tiny_y_offset = -0.01; // not sure why I need this
                                   if (mouse_in_opengl_coords.y >= (to_top_left.y + tiny_y_offset)) {
                                      // 0.03 is the bar height in setup_camera_facing_quad()
                                      float bar_height = 0.03;
                                      if (mouse_in_opengl_coords.y <= (to_top_left.y+tiny_y_offset+bar_height)) {
                                         coot::atom_spec_t spec(baddies[i].first);
                                         if (moving_atoms_asc->mol) {
                                            mmdb::Atom *at = spec.get_atom(moving_atoms_asc->mol);
                                            if (at) {
                                               at_out = at;
                                               status = true;
                                               if (act_on_hit) {
                                                  clipper::Coord_orth pt = coot::co(at);
                                                  std::cout << "INFO: geom bar atom: " << coot::atom_spec_t(at)
                                                            << std::endl;
                                                  set_rotation_centre(pt);
                                               }
                                            }
                                         } else {
                                            std::cout << "ERROR:: no moving atoms mol" << std::endl;
                                         }
                                      }
                                   }
                                }
                             }
                          }
                          return std::pair<bool, mmdb::Atom *>(status, at_out);
                       };

   auto check_rota_blocks = [mouse_in_opengl_coords,
                             rota_sorter] (unsigned int bar_index,
                                           float x_base_for_hud_geometry_bars,
                                           rotamer_markup_container_t *rotamer_markups,
                                           int n_rotamer_markups,
                                           bool act_on_hit) {

                               bool status = false;
                               mmdb::Atom *at_out = 0;
                               coot::residue_spec_t spec_for_at_out;
                               glm::vec2 to_top_left(x_base_for_hud_geometry_bars, 0.943 - 0.05 * static_cast<float>(bar_index));
                               std::vector<rotamer_markup_container_t> v;
                               // filter out the goodies
                               for (int i=0; i<n_rotamer_markups; i++)
                                  if (rotamer_markups[i].rpi.probability < 40) // 40 %
                                    if (rotamer_markups[i].rpi.probability >= 0)
                                       v.push_back(rotamer_markups[i]);
                               // sort the baddies
                               std::sort(v.begin(), v.end(), rota_sorter);
                               unsigned int n_rota_max = 20;
                               if (v.size() > n_rota_max) {
                                  unsigned int n_for_deletion = v.size() - n_rota_max;
                                  std::vector<rotamer_markup_container_t>::iterator v_begin = v.begin();
                                  std::vector<rotamer_markup_container_t>::iterator v_last  = v_begin + n_for_deletion;
                                  v.erase(v_begin, v_last);
                               }

                               float sum_l = 0;
                               for (unsigned int i=0; i<v.size(); i++) {
                                  // float pr = v[i].rpi.probability;
                                  float q = 0.01 * (48.0f - v[i].rpi.probability);
                                  if (q > 1.0) q = 1.0;
                                  if (q < 0.0) q = 0.0;
                                  float bar_length = std::pow(q, 6.0) * 4.0;
                                  glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);

                                  if (mouse_in_opengl_coords.x >= position_offset.x) {
                                     if (mouse_in_opengl_coords.x <= (position_offset.x + bar_length)) {
                                        // std::cout << ":::::::::: x hit bar_index " << bar_index
                                        //           << " i " << i << " " << baddies[i].first << std::endl;
                                        float tiny_y_offset = -0.01; // not sure why I need this
                                        if (mouse_in_opengl_coords.y >= (to_top_left.y + tiny_y_offset)) {
                                           // 0.03 is the bar height in setup_camera_facing_quad()
                                           float bar_height = 0.03;
                                           if (mouse_in_opengl_coords.y <= (to_top_left.y+tiny_y_offset+bar_height)) {

                                              if (false)
                                                 std::cout << "rama bar hit! " << i << " "
                                                           << v[i].spec << " "
                                                           << v[i].col << " "
                                                           << "probability " << v[i].rpi.probability << std::endl;

                                              if (moving_atoms_asc->mol) {
                                                 status = true;
                                                 spec_for_at_out = v[i].spec;
                                                 if (act_on_hit) {
                                                    clipper::Coord_orth pos = v[i].pos;
                                                    graphics_info_t::set_rotation_centre(pos);
                                                 }
                                              }
                                           }
                                        }
                                     }
                                  }
                                  sum_l += bar_length + 0.005; // with a gap between bars
                               }
                               // I need this function to return an atom (so that it's like the other geometry bars)
                               // - because that atom spec gets turned into an tooltip label.
                               if (status) {
                                  // I first need to find the first residue
                                  mmdb::Residue *residue_p = spec_for_at_out.get_residue(moving_atoms_asc->mol);
                                  if (residue_p) {
                                     int n_atoms = residue_p->GetNumberOfAtoms();
                                     if (n_atoms > 0) at_out = residue_p->GetAtom(0);
                                     if (n_atoms > 1) at_out = residue_p->GetAtom(1); // CA, usually
                                  }
                               }
                               return std::pair<bool, mmdb::Atom *>(status, at_out);
                            };

   float x_base_for_hud_geometry_bars = get_x_base_for_hud_geometry_bars();

   status_pair = check_blocks(rr.sorted_atom_pulls, 0, x_base_for_hud_geometry_bars,
                              hud_geometry_distortion_to_bar_size_atom_pull, act_on_hit);

   if (!status_pair.first) {

      if (rr.refinement_results_contain_overall_nbc_score) {

         // 20220503-PE add_bars() take the argument std::vector<std::pair<coot::atom_spec_t, float> > &baddies
         // but now rr.sorted_nbc_baddies is std::vector<refinement_results_nbc_baddie_t>
         // so now I need to convert

         std::vector<std::pair<coot::atom_spec_t, float> > converted_baddies(rr.sorted_nbc_baddies.size() * 2); // 20230519-PE both ways
         for (unsigned int i=0; i<rr.sorted_nbc_baddies.size(); i++) {
            const auto &bip = rr.sorted_nbc_baddies[i];
            std::pair<coot::atom_spec_t, float> p_1(bip.atom_spec_1, bip.score);
            std::pair<coot::atom_spec_t, float> p_2(bip.atom_spec_2, bip.score);
            // 20230813-PE fixes a crash, I hope.
            if ((2*i+1) < converted_baddies.size()) {
               converted_baddies[2*i  ] = p_1;
               converted_baddies[2*i+1] = p_2;
            } else {
               std::cout << "ERROR:: out of range in converted_baddies  " << 2*i << " " << converted_baddies.size() << std::endl;
            }
         }

         status_pair = check_blocks(converted_baddies, 1, x_base_for_hud_geometry_bars,
                                    hud_geometry_distortion_to_bar_size_nbc, act_on_hit);
      }
   }

   if (!status_pair.first)
      if (rr.refinement_results_contain_overall_rama_plot_score)
         status_pair = check_blocks(rr.sorted_rama_baddies, 2, x_base_for_hud_geometry_bars,
                                    hud_geometry_distortion_to_bar_size_rama, act_on_hit);

   if (!status_pair.first) {
      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            int nrms = moving_atoms_molecule.bonds_box.n_rotamer_markups;
            if (nrms > 0) {
               status_pair = check_rota_blocks(3, x_base_for_hud_geometry_bars,
                                               moving_atoms_molecule.bonds_box.rotamer_markups, nrms, act_on_hit);
            }
         }
      }
   }

   if (!status_pair.first) {
      if (! rr.sorted_chiral_volume_baddies.empty()) {
         if (moving_atoms_asc) {
            if (moving_atoms_asc->mol) {

               auto hud_geometry_distortion_to_bar_size_chiral = +[] (float distortion) {
                  return distortion * 0.01f; // the f is needed - match lambda function of the same name
                                             // in draw_hud_geometry_bars().
               };

               std::vector<std::pair<coot::atom_spec_t, float> > converted_baddies(rr.sorted_chiral_volume_baddies.size());
               for (unsigned int ii=0; ii<rr.sorted_chiral_volume_baddies.size(); ii++) {
                  converted_baddies[ii] = std::make_pair(rr.sorted_chiral_volume_baddies[ii].atom_spec,
                                                         rr.sorted_chiral_volume_baddies[ii].distortion);
               }
               status_pair = check_blocks(converted_baddies, 4, x_base_for_hud_geometry_bars,
                                          hud_geometry_distortion_to_bar_size_chiral, act_on_hit);
            }
         }
      }
   }

   if (act_on_hit) {
      if (status_pair.first) {
         mmdb::Atom *at = status_pair.second;
         if (at) {
            mmdb::Residue *residue_p = at->residue;
            moving_atoms_visited_residues.insert(residue_p);
            active_atom_for_hud_geometry_bar = at;
         }
      }
   }

   return status_pair;
}

std::pair<bool, mmdb::Atom *>
graphics_info_t::check_if_moused_over_hud_bar(double mouse_x, double mouse_y) {

   // copied from check_if_hud_bar_clicked().
   // Now that we have act_on_hit, we can extract most of this code to a common function
   // called check_if_hud_bar_mouse_over_or_act_on_hurd_bar_click()

   bool act_on_hit = false;
   auto r = check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(mouse_x, mouse_y, act_on_hit);
   if (false)
      std::cout << ":::::::::: debug:: check_if_moused_over_hud_bar() returns "
                << r.first << " " << r.second << std::endl;
   return r;
}

bool
graphics_info_t::check_if_hud_bar_clicked(double mouse_x, double mouse_y) {

   if (! moving_atoms_asc) return false;
   if (! moving_atoms_asc->mol) return false;
   bool act_on_hit = true;
   std::pair<bool, mmdb::Atom *> r = check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(mouse_x, mouse_y, act_on_hit);
   if (false)
      std::cout << ":::::::::: debug:: check_if_hud_bar_clicked() returns "
                << r.first << " " << r.second << std::endl;
   return  r.first;
}

bool
graphics_info_t::check_if_hud_button_moused_over(double mouse_x, double mouse_y, bool button_1_is_down) {

   // std::cout << "Here in check_if_hud_button_moused_over() with button_1_is_down " << button_1_is_down << std::endl;

   bool act_on_hit = false;
   check_if_hud_button_moused_over_or_act_on_hit(mouse_x, mouse_y, act_on_hit, button_1_is_down);
   return false;
}

bool
graphics_info_t::check_if_hud_button_clicked(double mouse_x, double mouse_y) {

   bool act_on_hit = true;
   bool status = check_if_hud_button_moused_over_or_act_on_hit(mouse_x, mouse_y, act_on_hit, false);
   return status;
}

// this function needs to be passed mouse press or mouse release button info
// so that it can do the button highlighting correctly.
//
// button_1_is_down is used for the highlighting.
//
bool
graphics_info_t::check_if_hud_button_moused_over_or_act_on_hit(double x, double y, bool act_on_hit, bool button_1_is_down) {

   auto highlight_just_button_with_index = [button_1_is_down] (unsigned int idx_active) {
                                              for (unsigned int i=0; i<hud_button_info.size(); i++) {
                                                 auto &button = hud_button_info[i];
                                                 if (i == idx_active) {
                                                    if (button_1_is_down) {
                                                       button.set_button_colour_for_mode(HUD_button_info_t::PRESSED);
                                                    } else {
                                                       button.set_button_colour_for_mode(HUD_button_info_t::HIGHLIGHTED);
                                                    }
                                                 } else {
                                                    button.set_button_colour_for_mode(HUD_button_info_t::BASIC);
                                                 }
                                              }
                                              GLenum err = glGetError();
                                              if (err) std::cout << "GL ERROR:: highlight_just_button_with_index pos-B "
                                                                 << err << std::endl;
                                              attach_buffers();
                                              err = glGetError();
                                              if (err) std::cout << "GL ERROR:: highlight_just_button_with_index pos-C "
                                                                 << err << std::endl;
                                              mesh_for_hud_buttons.update_instancing_buffer_data(hud_button_info);
                                              if (err) std::cout << "GL ERROR:: highlight_just_button_with_index pos-D "
                                                                 << err << std::endl;
                                              graphics_draw(); // let's see the changes then
                                           };
   auto unhighlight_all_buttons = [] () {
                                              for (unsigned int i=0; i<hud_button_info.size(); i++) {
                                                 auto &button = hud_button_info[i];
                                                 button.set_button_colour_for_mode(HUD_button_info_t::BASIC);
                                              }
                                              GLenum err = glGetError();
                                              if (err) std::cout << "GL ERROR:: unhighlight_all_buttons pos-B "
                                                                 << err << std::endl;
                                              attach_buffers();
                                              err = glGetError();
                                              if (err) std::cout << "GL ERROR:: unhighlight_all_buttons pos-C "
                                                                 << err << std::endl;
                                              mesh_for_hud_buttons.update_instancing_buffer_data(hud_button_info);
                                              err = glGetError();
                                              if (err) std::cout << "GL ERROR:: unhighlight_all_buttons pos-D "
                                                                 << err << std::endl;
                                  };

   bool status = false;
   bool clear_HUD_buttons_flag = false;
   if (! hud_button_info.empty()) {
      GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
      int w = allocation.width;
      int h = allocation.height;

      double x_gl_coords =  2.0 * x/static_cast<double>(w) - 1.0;
      double y_gl_coords = -2.0 * y/static_cast<double>(h) + 1.0;

      for (unsigned int i=0; i<hud_button_info.size(); i++) {
         const auto &button = hud_button_info[i];
         // are we on that button?
         HUD_button_limits_t lims = button.get_button_limits(w, h);
         if (lims.is_hit(x_gl_coords,y_gl_coords)) {
            if (act_on_hit) {
               // std::cout << "Act on button " << i << " callback" << std::endl;
               if (button.callback_function) {
                  button.callback_function();
               }
               // don't clear the hud buttons in the callback.
               if (button.button_label == "Accept") clear_HUD_buttons_flag = true;
               if (button.button_label == "Reject") clear_HUD_buttons_flag = true;
            }
            status = true;
            highlight_just_button_with_index(i);
         }
      }
      if (!status) {
         unhighlight_all_buttons();
      }
   }
   if (clear_HUD_buttons_flag) clear_hud_buttons();

   return status;
}


void
graphics_info_t::draw_hud_geometry_tooltip() {

   // this flag is set when the user mouses over a HUD bar
   // and removed when they move from a hud geometry bar.

   if (draw_hud_tooltip_flag) {

      glEnable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);

      bool draw_background = false; // 20211124-PE backgrounds are not correctly scaled at the moment

      if (draw_background) {
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

         texture_for_hud_tooltip_background.Bind(0);
         mesh_for_hud_tooltip_background.set_scales(glm::vec2(0.163, 0.05)); // hud-tooltip.png is 103x50
         mesh_for_hud_tooltip_background.draw(&shader_for_hud_geometry_labels);
      }

      // now the text that goes into (on top of) the background

      std::string label = "W 356 CA"; // checking the HUD bars (on mouse-over) should
                                     // return this - and the position, which means
                                     // that they need to return
                                     // something other than a bool. Store it in
                                     // graphics_info_t somewhere.
                                     // HUD_geometry_tooltip_text_position
                                     // HUD_geometry_tooltip_text_label
      label = label_for_hud_geometry_tooltip;

      // do this elsewhere (on_glarea_motion_notify())
      // glm::vec2 label_position(-0.64, 0.72);
      // tmesh_for_hud_geometry_tooltip_label.set_position(label_position);

      // this is now done in setup_hud_geometry_bars() which is called by resize()
      // glm::vec2 label_scale(0.00015, 0.00015); // fixed.
      // tmesh_for_hud_geometry_tooltip_label.set_scales(label_scale);

      bool use_label_highlight = true;
      mmdb::Residue *residue_p = 0;
      if (active_atom_for_hud_geometry_bar)
         residue_p = active_atom_for_hud_geometry_bar->residue;
      if (residue_p)
         if (moving_atoms_visited_residues.find(residue_p) != moving_atoms_visited_residues.end())
            use_label_highlight = false;

      // we don't want the residue label text to stretch when the window is wide
      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
      int w = allocation.width;
      int h = allocation.height;
      float aspect_ratio = static_cast<float>(w)/static_cast<float>(h);
      // 20220215-PE Hmmm about 0.00006 will do.
      glm::vec2 label_scale(0.00006, 0.00006 * aspect_ratio);
      label_scale *= 1.5f;
      tmesh_for_hud_geometry_tooltip_label.set_scales(label_scale);

      tmesh_for_hud_geometry_tooltip_label.draw_label(label, use_label_highlight,
                                                      &shader_for_hud_geometry_tooltip_text,
                                                      ft_characters);
   }
}

#include "analysis/stats.hh"

void
graphics_info_t::draw_hud_elements() {

   draw_hud_ligand_view();

   draw_hud_geometry_bars();

   draw_hud_geometry_tooltip(); // background and text

   draw_hud_ramachandran_plot();

   draw_hud_buttons();

   draw_hud_fps();

   draw_hud_refinement_dialog_arrow_tab();

   draw_hud_colour_bar();

}

bool
graphics_info_t::check_if_hud_rama_plot_clicked(double mouse_x, double mouse_y) {

   auto get_rotation_centre_from_intermediate_atoms_residue_spec = [] (const coot::residue_spec_t &rs) {
      bool status;
      coot::Cartesian pos;
      if (moving_atoms_asc) {
         mmdb::Manager *mol = moving_atoms_asc->mol;
         mmdb::Residue *residue_p = coot::util::get_residue(rs, mol);
         if (residue_p) {
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (! at->isTer()) {
                  std::string atom_name(at->GetAtomName());
                  if (atom_name == " CA ") {
                     pos = coot::Cartesian(at->x, at->y, at->z);
                     status = true;
                  }
               }
            }
         }
      }
      return std::make_pair(status, pos);
   };

   // c.f. draw_hud_ramachandran_plot()

   bool status = false;

   if (! moving_atoms_asc) return false;
   if (! moving_atoms_asc->mol) return false;

   if (draw_gl_ramachandran_plot_flag) {
      if (draw_gl_ramachandran_plot_user_control_flag) {
         if (moving_atoms_asc) {
            if (moving_atoms_asc->n_selected_atoms > 0) {
               int imol = imol_moving_atoms;
               //gl_rama_plot_t rama;
               // std::string residue_selection = "//";     // Hmm!
               // gl_rama_plot_t::draw_mode_t draw_mode = gl_rama_plot_t::draw_mode_t::CHECK_IF_PICKED;
               // rama.setup_from(imol, moving_atoms_asc->mol, residue_selection, draw_mode);
               GtkAllocation allocation;
               GtkWidget *gl_area = glareas[0];
               gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
               int w = allocation.width;
               int h = allocation.height;
               mouse_over_hit_t hit = gl_rama_plot.get_mouse_over_hit(mouse_x, mouse_y, w, h);
               if (false)
                  std::cout << "hit: plot clicked: " << hit.plot_was_clicked
                            << " residue_was_clicked: " << hit.residue_was_clicked
                            << " spec " << hit.residue_spec << std::endl;
               if (hit.plot_was_clicked) status = true;
               if (hit.residue_was_clicked) {
                  std::pair<bool, coot::Cartesian> rc = get_rotation_centre_from_intermediate_atoms_residue_spec(hit.residue_spec);
                  if (rc.first) {
                     setRotationCentre(rc.second, false);
                  }
               }
            }
         }
      }
   }
   return status;
}


void
graphics_info_t::render_3d_scene(GtkGLArea *gl_area, stereo_eye_t eye) {

   // note: this function is called from render_scene_sans_depth_blur()
   // 20230814-PE Is it?
   //             It is not used by the "Fancy" frame-buffer path

   //  ------------------- render scene ----------------------------

   // std::cout << "render_3d_scene() start" << std::endl;

   glEnable(GL_DEPTH_TEST);

   // const glm::vec3 &bg = graphics_info_t::background_colour;
   // glClearColor (bg[0], bg[1], bg[2], 1.0); // what difference does this make?

   GLenum err = glGetError(); if (err) std::cout << "render_3d_scene lambda B err " << err << std::endl;
   // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   err = glGetError(); if (err) std::cout << "render_3d_scene lambda C err " << err << std::endl;

   draw_origin_cube(eye, gl_area);
   err = glGetError(); if (err) std::cout << "render scene lambda post cubes err " << err << std::endl;

   draw_molecules(eye); // includes particles, happy-faces and boids (should they be there (maybe not))
                     // so rename this function? Or just bring everything here?  Put this render() function
                     // into new file graphics-info-opengl-render.cc

   draw_at_screen_centre_pulse(eye);

   draw_invalid_residue_pulse(eye);

   draw_delete_item_pulse(eye);

   draw_generic_pulses(eye);

   draw_measure_distance_and_angles(eye); // maybe in draw_molecules()?

   draw_extra_distance_restraints(eye, PASS_TYPE_STANDARD); // GM_restraints

   draw_pointer_distances_objects(eye);

   draw_translation_gizmo(eye); // maybe rotation gizmo too, later.

   draw_texture_meshes(eye);

}


void
graphics_info_t::render_3d_scene_with_shadows(stereo_eye_t eye) {

   // std::cout << "render_3d_scene_with_shadows() --- start ---" << std::endl;

   // note: this function is called from render_scene_sans_depth_blur()

   //  ------------------- render scene ----------------------------

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);

   // Try not using the meshes-with-shadows.shader
   // render_3d_scene(gl_area);
   // return;

   glEnable(GL_DEPTH_TEST);
   GLenum err = glGetError();
   if (err) std::cout << "render_3d_scene_with_shadows B err " << err << std::endl;

   draw_origin_cube(eye, gl_area);
   err = glGetError(); if (err) std::cout << "render scene lambda post cubes err " << err << std::endl;

   // draw_rotation_centre_crosshairs(gl_area, PASS_TYPE_STANDARD);

   draw_molecules_with_shadows(eye); // includes particles, happy-faces and boids (should they be there (maybe not))
                                  // so rename this function? Or just bring everything here?  Put this render() function
                                  // into new file graphics-info-opengl-render.cc


   draw_molecules_other_meshes(eye, PASS_TYPE_WITH_SHADOWS);

   draw_at_screen_centre_pulse(eye);


   draw_invalid_residue_pulse(eye);

   draw_generic_pulses(eye);

   draw_delete_item_pulse(eye);

   draw_measure_distance_and_angles(eye); // maybe in draw_molecules()?

   draw_pointer_distances_objects(eye);

   draw_extra_distance_restraints(eye, PASS_TYPE_WITH_SHADOWS); // GM_restraints. 20231121-PE is this the right pass type?

   draw_texture_meshes(eye);

}


void
graphics_info_t::render_3d_scene_for_ssao() {

}


// these are both optional arguments.
gboolean
graphics_info_t::render(bool to_screendump_framebuffer_flag, const std::string &output_file_name) {

   // auto tp_0 = std::chrono::high_resolution_clock::now();

   // auto draw_hud_elements = [] () { };
   // auto render_3d_scene = [] (GtkGLArea *gl_area) { };

   auto do_fps_std_dev_stuff = [] {
                              if (GetFPSFlag()) {
                                 unsigned int n_fps_history = frame_time_history_list.size();
                                 unsigned int n_history_max = 60;
                                 if (n_fps_history > 5) {
                                    coot::stats::single data;
                                    int n_history_count = n_fps_history - n_history_max;
                                    int count = 0;
                                    std::list<std::chrono::time_point<std::chrono::high_resolution_clock> >::const_iterator it;
                                    for (it = frame_time_history_list.begin(); it != frame_time_history_list.end(); ++it) {
                                       if (it != frame_time_history_list.begin()) {
                                          if (count > n_history_count) {
                                             const std::chrono::time_point<std::chrono::high_resolution_clock> &tp_this = *it;
                                             const std::chrono::time_point<std::chrono::high_resolution_clock> &tp_prev = *std::prev(it);
                                             auto delta_t = std::chrono::duration_cast<std::chrono::milliseconds>(tp_this - tp_prev).count();
                                             data.add(delta_t);
                                          }
                                          count++;
                                       }
                                    }
                                    if (data.size() > 5) {
                                       auto v = data.variance();
                                       fps_std_dev = sqrt(v);
                                    }
                                 }
                              }
                           };

   auto update_fps_statistics = [do_fps_std_dev_stuff] () {
                                   if (GetFPSFlag()) {
                                      frame_counter++;
                                      std::chrono::time_point<std::chrono::high_resolution_clock> tp_now = std::chrono::high_resolution_clock::now();
                                      auto delta_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - previous_frame_time_for_per_second_counter);
                                      auto elapsed_seconds = 0.001 * delta_time_ms;
                                      if (elapsed_seconds.count() >= 1.0) {
                                         float num_frames_delta = frame_counter - frame_counter_at_last_display;
                                         previous_frame_time_for_per_second_counter = tp_now;
                                         frame_counter_at_last_display = frame_counter;
                                         fps = num_frames_delta/elapsed_seconds.count();
                                         do_fps_std_dev_stuff();
                                      }
                                   }
                                };

   auto screendump_image = [update_fps_statistics] (const std::string &output_file_name) {
      // this works! Nice framebuffer scaling with screendump_tga().

      std::cout << "debug:: in screendump_image() with use_framebuffers " << use_framebuffers << std::endl;

      GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
      int w = allocation.width;
      int h = allocation.height;

// #ifdef __APPLE__
//      use_framebuffers = false;
// #endif

      if (use_framebuffers) { // static class variable

         glViewport(0, 0, framebuffer_scale * w, framebuffer_scale * h);
         GLenum err = glGetError();
         if (err) std::cout << "GL ERROR:: render() post glViewport() err " << err << std::endl;
         screen_framebuffer.bind(); // screen_ao, that is
         err = glGetError();
         if (err) std::cout << "GL ERROR:: render() post screen_framebuffer bind() err " << err << std::endl;

         render_3d_scene(gl_area, stereo_eye_t::MONO);

         // screendump
         glDisable(GL_DEPTH_TEST);
         const unsigned int &sf = framebuffer_scale;
         glViewport(0, 0, sf * w, sf * h);
         framebuffer screendump_framebuffer;
         unsigned int index_offset = 0;
         screendump_framebuffer.init(sf * w, sf * h, index_offset, "screendump");
         screendump_framebuffer.bind();

         // render_3d_scene(gl_area);
         // render_scene_with_screen_ao_shader();

         render_scene();

         // gtk_gl_area_attach_buffers(gl_area);
         screendump_tga_internal(output_file_name, w, h, sf, screendump_framebuffer.get_fbo());

      } else {

         gtk_gl_area_attach_buffers(gl_area);
         render_3d_scene(gl_area, stereo_eye_t::MONO);
         draw_hud_elements();

      }

      // 20211112-PE
      // This seems to do bad things to the frame-rate on the PC (although fullscreen mode seems
      // unaffected and looks to be *faster* than windowed mode (could be a gtk thing)).
      // This is vital to see anything sane on the Mac.
      glFlush();

      // auto tp_1 = std::chrono::high_resolution_clock::now();
      // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      // std::cout << "INFO:: render() " << d10 << " microseconds" << std::endl;

      // std::cout << "calling update_fps_statistics() " << std::endl;
      update_fps_statistics();

      return FALSE;

   };

   auto update_frame_time_history = [] () {
      unsigned int frame_time_history_list_max_n_elements = 500;
      GtkWidget *glarea = glareas[0];
      if (glarea) {
         auto tp_now = std::chrono::high_resolution_clock::now();
         frame_time_history_list.push_back(tp_now);
         if (frame_time_history_list.size() >= (frame_time_history_list_max_n_elements+1))
            frame_time_history_list.pop_front();
      }
   };


   // ################################## main line #############################################

   update_frame_time_history();

   if (to_screendump_framebuffer_flag) {
      return screendump_image(output_file_name);
   } else {
      gboolean state = render_scene();
      draw_hud_elements();
#ifdef __APPLE__
      glFinish();
#else
      glFlush();
#endif
      update_fps_statistics();
      return state;
   }
}

void
graphics_info_t::render_scene_with_x_blur() {

   shader_for_x_blur.Use();
   glBindVertexArray(blur_x_quad_vertex_array_id);

   const glm::vec3 &bg = background_colour;
   glClearColor(bg[0], bg[1], bg[2], 1.0); // needed?
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glActiveTexture(GL_TEXTURE0 + 0);
   glBindTexture(GL_TEXTURE_2D, blur_x_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, blur_x_framebuffer.get_texture_depth());
   shader_for_x_blur.set_int_for_uniform("screenTexture", 0);
   GLenum err = glGetError(); if (err) std::cout << "GL ERROR:: render_scene_with_x_blur() D err " << err << std::endl;

   glDrawArrays(GL_TRIANGLES, 0, 6);
   err = glGetError(); if (err) std::cout << "GL ERROR:: render_scene_with_x_blur() E err " << err << std::endl;

}

void
graphics_info_t::render_scene_with_y_blur() {

   shader_for_y_blur.Use();
   glBindVertexArray(blur_y_quad_vertex_array_id);

   const glm::vec3 &bg = background_colour;
   glClearColor(bg[0], bg[1], bg[2], 1.0); // needed?
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glActiveTexture(GL_TEXTURE0 + 0);
   glBindTexture(GL_TEXTURE_2D, blur_y_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, blur_y_framebuffer.get_texture_depth());
   shader_for_y_blur.set_int_for_uniform("screenTexture", 0);
   GLenum err = glGetError(); if (err) std::cout << "GL ERROR:: render_scene_with_x_blur() D err " << err << std::endl;

   glDrawArrays(GL_TRIANGLES, 0, 6);
   err = glGetError(); if (err) std::cout << "GL ERROR:: render_scene_with_x_blur() E err " << err << std::endl;
}

void
graphics_info_t::render_scene_with_texture_combination_for_depth_blur() {

   shader_for_dof_blur_by_texture_combination.Use();
   glBindVertexArray(combine_textures_using_depth_quad_vertex_array_id);

   const glm::vec3 &bg = background_colour;
   glClearColor(bg[0], bg[1], bg[2], 1.0); // needed?
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // shader:
   // uniform sampler2D screenTexture1;
   // uniform sampler2D screenTexture2;
   // uniform sampler2D screenDepth;

   shader_for_dof_blur_by_texture_combination.set_bool_for_uniform("do_outline_mode", shader_do_outline_flag);

   shader_for_dof_blur_by_texture_combination.set_float_for_uniform("focus_blur_z_depth",  focus_blur_z_depth);
   shader_for_dof_blur_by_texture_combination.set_float_for_uniform("focus_blur_strength", focus_blur_strength);

   glActiveTexture(GL_TEXTURE0 + 0);
   glBindTexture(GL_TEXTURE_2D, combine_textures_using_depth_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, blur_y_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 2);
   // glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_depth()); // 20220222-PE pre-crow code
   glBindTexture(GL_TEXTURE_2D, framebuffer_for_effects.get_texture_depth());

   shader_for_dof_blur_by_texture_combination.set_int_for_uniform("screenTexture1", 0);
   shader_for_dof_blur_by_texture_combination.set_int_for_uniform("screenTexture2", 1);
   shader_for_dof_blur_by_texture_combination.set_int_for_uniform("screenDepth",    2);
   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: render_scene_with_texture_combination_for_depth_blur() D err " << err << std::endl;

   glDrawArrays(GL_TRIANGLES, 0, 6);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: render_scene_with_texture_combination_for_depth_blur() E err " << err << std::endl;
}



void
graphics_info_t::reset_frame_buffers(int window_width, int window_height) {

   if (false)
      std::cout << "DEBUG:: reset_frame_buffers() " << window_width << " " << window_height
                << " use_framebuffers: " << use_framebuffers << std::endl;

   if (use_framebuffers) {

      // 20220108-PE note to self. Try using the framebuffer::reset() function instead

      unsigned int sf = framebuffer_scale; // this is set by set_framebuffer_scale_factor()
      unsigned int index_offset = 0;

      // width  = width;
      // height = height;
      if (false)
         std::cout << "debug:: reset_frame_buffers() with sf " << sf << " "
                   << window_width << " x " << window_height << std::endl;

      screen_framebuffer.init(sf * window_width, sf * window_height, index_offset, "screen");
      GLenum err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      blur_x_framebuffer.init(sf * window_width, sf * window_height, index_offset, "blur-x");
      err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      blur_y_framebuffer.init(sf * window_width, sf * window_height, index_offset, "blur-y");
      err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      combine_textures_using_depth_framebuffer.init(sf * window_width, sf * window_height, index_offset, "combine");
      err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      // std::cout << "debug:: reset_frame_buffers() sf " << sf << " width " << width << " height " << height << std::endl;

      // index_offset = 0;
      // g.blur_framebuffer.init(width, height, index_offset, "blur");

      // ------------------ crows code --------------------------

      // note to self:
      // the shadow texture doesn't need to change - it's under user control, not
      // dependent on the window size

      framebuffer_for_ssao_gbuffer.reset_test(window_width, window_height);

      gint w = window_width;
      gint h = window_height;

      // cut and paste from init_joey_ssao_stuff() for now - do better later.

      {
         glBindFramebuffer(GL_FRAMEBUFFER, ssaoFBO);
         glBindTexture(GL_TEXTURE_2D, ssaoColorBuffer);
         glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_FLOAT, NULL);
         glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaoColorBuffer, 0);

         glBindFramebuffer(GL_FRAMEBUFFER, ssaoBlurFBO);
         glBindTexture(GL_TEXTURE_2D, ssaoColorBufferBlur);
         glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_FLOAT, NULL);
         glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaoColorBufferBlur, 0);
         glBindFramebuffer(GL_FRAMEBUFFER, 0);
      }

      // the render bufffer rboDepth does something related to the SSAO
      {
         glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
         glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h);
      }
   }

}

void
graphics_info_t::try_label_unlabel_active_atom() {

   std::pair<int, mmdb::Atom *> aa = get_active_atom();
   int im = aa.first;
   if (im >= 0) {
      mmdb::Atom *at = aa.second;
      if (at) {
         int atom_index;
         // this is a bit convoluted :-)
         int ierr = at->GetUDData(molecules[im].atom_sel.UDDAtomIndexHandle, atom_index);
	 if (ierr == mmdb::UDDATA_Ok) {
            molecules[im].add_to_labelled_atom_list(atom_index);
	    add_picked_atom_info_to_status_bar(im, atom_index);
            graphics_draw();
         } else {
            std::cout << "WARNING:: Bad UDData for atom_index for atom " << std::endl;
         }
      }
   }
}


// static
glm::vec3
graphics_info_t::get_screen_y_uv() {

   glm::vec3 minus_y = graphics_info_t::unproject_to_world_coordinates(glm::vec3(0.0f, -1.0f, 0.0f));
   glm::vec3  plus_y = graphics_info_t::unproject_to_world_coordinates(glm::vec3(0.0f,  1.0f, 0.0f));
   glm::vec3 delta = plus_y - minus_y;
   glm::vec3 d_uv = glm::normalize(delta);
   return d_uv;
}

// static
glm::vec3
graphics_info_t::get_screen_x_uv() {

   glm::vec3 minus_x = graphics_info_t::unproject_to_world_coordinates(glm::vec3(-1.0f, 0.0f, 0.0f));
   glm::vec3  plus_x = graphics_info_t::unproject_to_world_coordinates(glm::vec3( 1.0f, 0.0f, 0.0f));
   glm::vec3 delta = plus_x - minus_x;
   glm::vec3 d_uv = glm::normalize(delta);
   return d_uv;
}


void
graphics_info_t::translate_in_screen_z(float step_size) {

   // The step size is good when were zoomed in but too big when we are zoomed out.

   // this looks a bit weird without perspective view

   glm::vec3 ep = get_world_space_eye_position();
   glm::vec3 rc = get_rotation_centre();
   glm::vec3 delta = rc - ep;
   glm::vec3 delta_uv = normalize(delta);

   // more zoomed in has smaller zoom than zoomed out. Zoomed out is ~100. Zoomed in is ~25
   glm::vec3 step = 0.005f * step_size * zoom * delta_uv;

   if (false) // debug
      std::cout << "ep " << glm::to_string(ep) << " rc " << glm::to_string(rc)
                << " zoom " << zoom << " step " << glm::to_string(step) << std::endl;

   add_to_rotation_centre(step);

}

void
graphics_info_t::translate_in_screen_x(float step_size) {

   // The step size is good when were zoomed in but too big when we are zoomed out.

   glm::vec3 screen_x_uv = get_screen_x_uv();
   glm::vec3 step = 0.005f * step_size * zoom * screen_x_uv;
   add_to_rotation_centre(step);
}



// static
std::vector<glm::vec3>
graphics_info_t::get_particle_centre_positions() {

   auto mmdb_to_glm = [] (mmdb::Atom *at) { return glm::vec3(at->x, at->y, at->z); };

   get_moving_atoms_lock(__FUNCTION__);

   std::vector<glm::vec3> v;
   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
            mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
            if (! at->isTer()) {
               std::string atom_name(at->GetAtomName());
               if (atom_name == " CA " || atom_name == " N1 " || atom_name == " N9 ") {
                  glm::vec3 p = mmdb_to_glm(at);
                  v.push_back(p);
               }
            }
         }
      }
   }
   release_moving_atoms_lock(__FUNCTION__);

   if (v.empty()) {
      glm::vec3 rc = get_rotation_centre();
      v.push_back(rc);
   }

   return v;
}

void
graphics_info_t::setup_draw_for_particles() {

   if (false) // from the days when particle drawing was a problem!
      std::cout << "setup_draw_for_particles(): -- start -- n_particles " << particles.size()
                <<  std::endl;

   if (particles.empty()) {
      std::cout << "setup_draw_for_particles(): let's make new particles " << std::endl;

      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
      GLenum err = glGetError();
      if (err) std::cout << "Error:: setup_draw_for_particles() Post attach buffers err is "
                         << err << std::endl;

      // 20231202-PE we don't need to use the shader, do we?
      shader_for_particles.Use();

      err = glGetError();
      if (err) std::cout << "GL ERROR:: setup_draw_for_particles() Post Use() err is "
                         << err << std::endl;

      std::vector<glm::vec3> positions = get_particle_centre_positions();
      particles.make_particles(n_particles, positions);
      std::cout << "setup_draw_for_particles(): done making " << n_particles << " particles "
                << " for " << positions.size() << " positions" << std::endl;

      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed? 20240610-PE, yes, I think so.
      mesh_for_particles.setup_vertex_and_instancing_buffers_for_particles(particles.size());
      mesh_for_particles.update_instancing_buffer_data_for_particles(particles);
      glUseProgram(0);
   }
   // passing user_data and Notify function at the end
   if (! do_tick_particles) {
      if (! tick_function_is_active()) {
         int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
         idle_function_spin_rock_token = new_tick_id;
      }
      do_tick_particles = true;
   }

   // std::cout << "setup_draw_for_particles(): -- done -- " << std::endl;
}

void
graphics_info_t::setup_draw_for_happy_face_residue_markers_init() {

   // run this once - call from realize()

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: setup_draw_for_happy_face_residue_markers_init() -- start -- "
                      << std::endl;

   const unsigned int max_happy_faces = 200; // surely enough?

   // 20221005-PE this causes an error on startup (like hud buttons and hud geometry bars)
   // gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   // err = glGetError();
   // if (err) std::cout << "GL ERROR:: setup_draw_for_happy_face_residue_markers_init() "
   //                    << "Post attach buffers err is " << err << std::endl;

   // If not found in this directory, then try default directory.
   // texture_for_happy_face_residue_marker.set_default_directory(coot::package_data_dir());
   texture_for_happy_face_residue_marker.init("happy-face-marker.png");

   // shader_for_happy_face_residue_markers.Use(); // needed?
   tmesh_for_happy_face_residues_markers.setup_camera_facing_quad(0.8, 0.8, 0.0, 0.0);
   tmesh_for_happy_face_residues_markers.setup_instancing_buffers(max_happy_faces);
   tmesh_for_happy_face_residues_markers.draw_this_mesh = false;

   err = glGetError();
   if (err) std::cout << "GL ERROR::- setup_draw_for_happy_face_residue_markers_init() "
                      << "--- end --- err is " << err << std::endl;
}

void
graphics_info_t::setup_draw_for_anchored_atom_markers_init() {

   // run this once - called from realize()


   // 20221005-PE this causes an error on startup (like hud buttons and hud geometry bars)
   //             and setup_draw_for_happy_face_residue_markers_init()
   // attach_buffers();
   // GLenum err = glGetError();
   // if (err) std::cout << "Error::- setup_draw_for_anchored_atom_markers_init() "
   //                    << "--- start --- err is " << err << std::endl;

   const unsigned int max_anchored_atoms = 200;

   // attach_buffers();

   GLenum err = glGetError();
   if (err) std::cout << "Error::- setup_draw_for_anchored_atom_markers_init() "
                      << "Post attach_buffers() err is " << err << std::endl;


   texture_for_anchored_atom_markers.init("anchor-for-fixed-atoms.png");
   // texture_for_anchored_atom_markers.Bind(0); // why is this needed?
   tmesh_for_anchored_atom_markers.setup_camera_facing_quad(0.3, 0.3, 0.0, 0.0);
   tmesh_for_anchored_atom_markers.setup_instancing_buffers(max_anchored_atoms);
   tmesh_for_anchored_atom_markers.draw_this_mesh = false;

}

void
graphics_info_t::setup_draw_for_happy_face_residue_markers() {

   // run this at the start of a "show animated happy faces"

   std::vector<glm::vec3> positions = get_happy_face_residue_marker_positions();
   happy_face_residue_marker_starting_positions = positions;

   glm::vec3 up_uv = get_screen_y_uv();
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?

   if (false)
      std::cout << "setup_draw_for_happy_face_residue_markers() calling update_instancing_buffer_data()"
                << " with draw_count_for_happy_face_residue_markers "
                << draw_count_for_happy_face_residue_markers << std::endl;
   unsigned int n_max = draw_count_max_for_happy_face_residue_markers;
   tmesh_for_happy_face_residues_markers.update_instancing_buffer_data_for_happy_faces(positions, 0, n_max, up_uv);
   tmesh_for_happy_face_residues_markers.draw_this_mesh = true;
   draw_count_for_happy_face_residue_markers = 0;
   if (! tick_function_is_active()) {
      tick_function_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
   }
   do_tick_happy_face_residue_markers = true;

}

void
graphics_info_t::setup_draw_for_anchored_atom_markers() {

   // these move frame by frame (on moving atoms coordinates updates)
   // 
   auto get_intermediate_atoms_anchored_atoms_positions = [] () {
                                          std::vector<glm::vec3> positions;
                                          if (moving_atoms_asc) {
                                             if (moving_atoms_asc->n_selected_atoms > 0) {
                                                for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
                                                   std::vector<coot::atom_spec_t> fixed = molecules[imol_moving_atoms].get_fixed_atoms();
                                                   for (unsigned int ifixed=0; ifixed<fixed.size(); ifixed++) {
                                                   }
                                                }
                                             }
                                          }
                                          return positions;
                                       };

   // these are fixed on the molecule
   //
   auto get_anchored_atoms_positions = [] (int imol, const glm::vec3 up_uv) {
                                          std::vector<glm::vec3> positions;
                                          const auto &m = molecules[imol];
                                          std::vector<coot::Cartesian> fap = m.fixed_atom_positions;
                                          for (unsigned int i=0; i<fap.size(); i++) {
                                             glm::vec3 p(fap[i].x(), fap[i].y(), fap[i].z());
                                             p += up_uv * 0.3f;
                                             positions.push_back(p);
                                          }
                                          return positions;
                                       };

   std::pair<int, mmdb::Atom *> aa = get_active_atom();
   if (aa.second) {
      int imol = aa.first;
      glm::vec3 up_uv = get_screen_y_uv();
      std::vector<glm::vec3> positions = get_anchored_atoms_positions(imol, up_uv);

      if (positions.empty()) {
         tmesh_for_anchored_atom_markers.draw_this_mesh = false;
      } else {
         attach_buffers();
         tmesh_for_anchored_atom_markers.draw_this_mesh = true;
         tmesh_for_anchored_atom_markers.update_instancing_buffer_data(positions);
      }
   }
}

#include "coot-utils/fib-sphere.hh"

std::vector<glm::vec3>
graphics_info_t::get_happy_face_residue_marker_positions() {

   const unsigned int max_happy_faces = 200; // surely enough? - If changed, change above
   std::vector<glm::vec3> v;
   bool make_fake_points = false;

   auto clipper_to_glm = [] (const clipper::Coord_orth &co) {
                            return glm::vec3(co.x(), co.y(), co.z());
                         };

   if (make_fake_points) {
      glm::vec3 sc = get_rotation_centre();
      std::vector<clipper::Coord_orth> cv = coot::fibonacci_sphere(80);
      for (auto p : cv)
         v.push_back(sc + 5.0f * clipper_to_glm(p));
   } else {

      // This is just a bit of fun... actually, I will need to ask something like
      // last_restraints->get_improved_residues();
      // last_restraints->get_damaged_residues();

      // How about the set-fixed-during-refinement udd?
      // see set_fixed_during_refinement_udd()

      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            int uddHnd = moving_atoms_asc->mol->GetUDDHandle(mmdb::UDR_ATOM , "FixedDuringRefinement");
            std::vector<mmdb::Residue *> residues;
            int imod = 1;
            mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        // I need to not add this if it's a fixed residue. How can I know that?
                        bool is_fixed = false;
                        mmdb::Atom **residue_atoms = 0;
                        int n_residue_atoms = 0;
                        residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                        for(int iat=0; iat<n_residue_atoms; iat++) {
                           mmdb::Atom *at = residue_atoms[iat];
                           int is_fixed_status = 0;
                           int ierr = at->GetUDData(uddHnd, is_fixed_status);
                           if (ierr == mmdb::Error_Ok) {
                              if (is_fixed_status == 1) {
                                 is_fixed = true;
                                 break;
                              }
                           }
                        }
                        if (! is_fixed)
                           residues.push_back(residue_p);
                     }
                  }
               }
            }

            for (auto r : residues) {
               std::pair<bool, clipper::Coord_orth> rc = coot::util::get_residue_centre(r);
               if (rc.first) {
                  glm::vec3 p = clipper_to_glm(rc.second);
                  v.push_back(p);
               }
            }
         }
      }
   }


   if (v.size() > max_happy_faces)
      std::cout << "error:: ------------------ too many happy faces" << std::endl;

   return v;
}

#include "sound.hh"

// static
void
graphics_info_t::update_bad_nbc_atom_pair_marker_positions() {

   auto gone_contacts_from_nbc_baddies = [] (const coot::refinement_results_t &rr_c) {

                                            // 20220505-PE note to self: previous_round_nbc_baddies_atom_index_map needs to be updated
                                            // every round - before we do a last_restraints->minimize(). Worry about these vectors
                                            // going out of scope. This may need locking? This function needs to be called from the main
                                            // thread of course.

                                            std::vector<std::pair<int, int> > gone_atom_pairs;

                                            const std::map<int, std::vector<int> > &current_round_nbc_baddies_atom_index_map = rr_c.nbc_baddies_atom_index_map;

                                            std::map<int, std::vector<int> >::const_iterator it_1;
                                            for (it_1=previous_round_nbc_baddies_atom_index_map.begin();
                                                 it_1!=previous_round_nbc_baddies_atom_index_map.end(); ++it_1) {
                                               const int &index_1 = it_1->first;
                                               const std::vector<int> &v_1(it_1->second);
                                               std::vector<int>::const_iterator it_2;
                                               for (it_2=v_1.begin(); it_2!=v_1.end(); ++it_2) {
                                                  const int &index_2(*it_2);
                                                  // can I find that index_1, index_2 pair in the current set?
                                                  std::map<int, std::vector<int> >::const_iterator it_3 = current_round_nbc_baddies_atom_index_map.find(index_1);
                                                  if (it_3 == current_round_nbc_baddies_atom_index_map.end()) {
                                                     std::pair<int, int> p(index_1, index_2);
                                                     gone_atom_pairs.push_back(p);
                                                  } else {
                                                     const std::vector<int> &v_2(it_3->second);
                                                     // so the first atom was there - what about the second atom?
                                                     std::vector<int>::const_iterator it_4 = std::find(v_2.begin(), v_2.end(), index_2);
                                                     if (it_4 == v_2.end()) {
                                                        std::pair<int, int> p(index_1, index_2);
                                                        gone_atom_pairs.push_back(p);
                                                     }
                                                  }
                                               }
                                            }
                                            return gone_atom_pairs;
                                         };

   auto get_gone_nbc_baddie_positions = [] (const std::vector<std::pair<int, int> > &gone_atom_pairs,
                                            mmdb::Atom **atom_selection, int n_selected_atoms) {
      glm::vec3 y_screen = get_screen_y_uv();
      std::vector<glm::vec3> positions;
      std::vector<std::pair<int, int> >::const_iterator it;
      for (it=gone_atom_pairs.begin(); it!=gone_atom_pairs.end(); ++it) {
         const auto &pair(*it);
         if (pair.first < n_selected_atoms) {
            if (pair.second < n_selected_atoms) {
               mmdb::Atom *at_1 = atom_selection[pair.first];
               mmdb::Atom *at_2 = atom_selection[pair.second];
               float x = 0.5 * (at_1->x + at_2->x);
               float y = 0.5 * (at_1->y + at_2->y);
               float z = 0.5 * (at_1->z + at_2->z);
               positions.push_back(glm::vec3(x,y,z) + 0.84f * y_screen);
            }
         }
      }
      return positions;
   };

   auto get_gone_count = [] (const std::map<int, std::vector<int> > &nbc_baddies_atom_index_map) {
      unsigned int n = 0;
      std::map<int, std::vector<int> >::const_iterator it_1;
      for (it_1=nbc_baddies_atom_index_map.begin(); it_1!=nbc_baddies_atom_index_map.end(); ++it_1) {
         const std::vector<int> &v_1(it_1->second);
         n += v_1.size();
      }
      return n;
   };

   auto nbc_baddies_count_delta = [get_gone_count] (const std::map<int, std::vector<int> > &nbc_baddies_atom_index_map_prev,
                                                    const std::map<int, std::vector<int> > &nbc_baddies_atom_index_map_new) {
      int n1 = get_gone_count(nbc_baddies_atom_index_map_prev); // it's just a count, not a gone count.
      int n2 = get_gone_count(nbc_baddies_atom_index_map_new);
      std::cout << "n1: " << n1 << " n2: " << n2 << std::endl;
      return n2-n1;
   };

   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         coot::refinement_results_t &rr = saved_dragged_refinement_results;

         if (false) {
            unsigned int gone_count_prev = get_gone_count(previous_round_nbc_baddies_atom_index_map);
            unsigned int gone_count_this = get_gone_count(rr.nbc_baddies_atom_index_map);
            std::cout << "compare nbc sizes: " << gone_count_prev << " " << gone_count_this << std::endl;
         }

         // if (nbc_baddies_count_delta(previous_round_nbc_baddies_atom_index_map, rr.nbc_baddies_atom_index_map) > 1)
         // play_sound("diego-arrives");

         int bad_nbc_atom_pair_marker_positions_size_pre = bad_nbc_atom_pair_marker_positions.size();
         bad_nbc_atom_pair_marker_positions.clear();
         std::vector<coot::refinement_results_nbc_baddie_t> &baddies(rr.sorted_nbc_baddies);
         for (unsigned int i=0; i<baddies.size(); i++) {
            bad_nbc_atom_pair_marker_positions.push_back(coord_orth_to_glm(baddies[i].mid_point));
         }
         int bad_nbc_atom_pair_marker_positions_size_post = bad_nbc_atom_pair_marker_positions.size();

         int bad_nbc_size_delta = bad_nbc_atom_pair_marker_positions_size_post - bad_nbc_atom_pair_marker_positions_size_pre;
         if (bad_nbc_size_delta > 0)
            play_sound("diego-arrives");  // maybe new-bump

         GLenum err = glGetError();
         if (err)
            std::cout << "GL ERROR:: update_bad_nbc_atom_pair_marker_positions() pos-B " << stringify_error_code(err) << std::endl;

         attach_buffers();
         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: update_bad_nbc_atom_pair_marker_positions() pos-C - post attach_buffers "
                      << stringify_error_code(err) << std::endl;
         tmesh_for_bad_nbc_atom_pair_markers.draw_this_mesh = true;
         tmesh_for_bad_nbc_atom_pair_markers.update_instancing_buffer_data(bad_nbc_atom_pair_marker_positions);
         err = glGetError();
         if (err)
            std::cout << "GL ERROR:: update_bad_nbc_atom_pair_marker_positions() pos-C - post update_instancing_buffer_data "
                      << stringify_error_code(err) << std::endl;
         if (! bad_nbc_atom_pair_marker_positions.empty())
            draw_bad_nbc_atom_pair_markers_flag = true;

         if (true) { // gone diego particles
            std::vector<std::pair<int, int> > gone_atom_pairs = gone_contacts_from_nbc_baddies(rr);
            // std::cout << "debug:: gone_atoms_pair size " << gone_atom_pairs.size() << std::endl;
            // for (unsigned int i=0; i<gone_atom_pairs.size(); i++)
            // std::cout << "       gone " << gone_atom_pairs[i].first << " " << gone_atom_pairs[i].second << std::endl;

            if (! gone_atom_pairs.empty()) {

               std::vector<glm::vec3> gone_diego_positions = get_gone_nbc_baddie_positions(gone_atom_pairs,
                                                                                           moving_atoms_asc->atom_selection,
                                                                                           moving_atoms_asc->n_selected_atoms);
               if (!gone_diego_positions.empty()) {
                  setup_draw_for_particles_for_new_gone_diegos(gone_diego_positions);
               }
            }
         }
         // for next round
         previous_round_nbc_baddies_atom_index_map = rr.nbc_baddies_atom_index_map;
      } else {
         bad_nbc_atom_pair_marker_positions.clear();
      }
   } else {
      bad_nbc_atom_pair_marker_positions.clear();
   }
}

#include "coot-utils/cylinder-utils.hh"

// static
void graphics_info_t::update_bad_nbc_atom_pair_dashed_lines() {

   // look at above

   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &vertices) {
      std::vector<s_generic_vertex> v_out(vertices.size());
      for (unsigned int i=0; i<vertices.size(); i++) {
         const auto &v = vertices[i];
         v_out[i].pos    = v.pos;
         v_out[i].normal = v.normal;
         v_out[i].color  = v.color;
      }
      return v_out;
   };

   auto baddies_to_dashed_lines = [] (const std::vector<coot::refinement_results_nbc_baddie_t> &baddies,
                                      unsigned int n_dashes) {
      dashed_cylinders_info_t dci;
      std::vector<std::pair<glm::vec3, glm::vec3> > positions;
      for (unsigned int i=0; i<baddies.size(); i++) {
         const auto &baddie = baddies[i];
         glm::vec3 p1(baddie.atom_1_pos.x(), baddie.atom_1_pos.y(), baddie.atom_1_pos.z());
         glm::vec3 p2(baddie.atom_2_pos.x(), baddie.atom_2_pos.y(), baddie.atom_2_pos.z());
         positions.push_back(std::make_pair(p1, p2));
      }
      dci = get_dashed_cylinders(positions, n_dashes);
      return dci;
   };

   // Note:
   //    class dashed_cylinders_info_t {
   // public:
   //    dashed_cylinders_info_t() {}
   //    cylinder c;
   //    glm::mat3 rot;
   //    glm::vec3 scales;
   //    std::vector<glm::vec3> offsets;
   // };

   GLenum err = glGetError();
   if (err)
      logger.log(log_t::GL_ERROR, logging::function_name_t("update_bad_nbc_atom_pair_dashed_lines"),
                 "--start--", stringify_error_code(err));

   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         coot::refinement_results_t &rr = saved_dragged_refinement_results;

         unsigned int n_dashes = 23;

         attach_buffers();
         err = glGetError();
         if (err)
            logger.log(log_t::GL_ERROR, logging::function_name_t("update_bad_nbc_atom_pair_dashed_lines"),
                       "A", stringify_error_code(err));

         // make instances
         std::vector<coot::refinement_results_nbc_baddie_t> &baddies(rr.sorted_nbc_baddies);
         if (! baddies.empty()) {
            dashed_cylinders_info_t dl = baddies_to_dashed_lines(baddies, n_dashes);
            std::vector<glm::mat4> mats;
            std::vector<glm::vec4> colours;
            for (unsigned int i=0; i<dl.oris_and_offsets.size(); i++) {
               const auto &ori = dl.oris_and_offsets[i].first;
               for (unsigned int j=0; j<dl.oris_and_offsets[i].second.size(); j++) {
                  const auto &offset = dl.oris_and_offsets[i].second[j];
                  glm::mat4 u(1.0f);
                  glm::mat4 s = glm::scale(u, dl.scales);
                  glm::mat4 r = dl.oris_and_offsets[i].first;
                  glm::mat4 t = glm::translate(u, offset);
                  glm::mat4 m1 = t * r * s;
                  mats.push_back(m1);
                  colours.push_back(glm::vec4(0.6, 0.4, 0.2, 1.0));
               }
            }
            Shader *shader_p = nullptr;
            unsigned int n_instances = mats.size();
            Material material;
            if (n_instances > 0) {
               for (unsigned int i=0; i<mats.size(); i++) {
                  // std::cout << "mat " << i << " " << glm::to_string(mats[i]) << std::endl;
               }
            }
            bad_nbc_atom_pair_dashed_line.setup_rtsc_instancing(shader_p, mats, colours, n_instances, material);
         }
      }
   }
}


void
graphics_info_t::setup_draw_for_particles_for_new_gone_diegos(const std::vector<glm::vec3> &positions) {

    // usually only one

    // gone_diego_particles and meshes_for_gone_diego_particles live and die together.
    // formalise that.

    // std::cout << "setup_draw_for_particles_for_gone_diegos() of " << positions.size() << std::endl;

    if (! positions.empty()) {

       play_sound("diego-gone-pop");

       glm::vec3 screen_x_uv = get_screen_x_uv();
       glm::vec3 screen_y_uv = get_screen_y_uv();

       meshed_particle_container_t mp(Mesh("gone-diego"), particle_container_t());
       meshed_particles_for_gone_diegos.push_back(mp);
       particle_container_t &last_particles = meshed_particles_for_gone_diegos.back().particle_container;
       Mesh &last_mesh                      = meshed_particles_for_gone_diegos.back().mesh;

       attach_buffers();
       unsigned int n_particles_per_burst = 10;
       int n_instances = n_particles_per_burst * positions.size();
       last_particles.make_gone_diego_particles(n_particles_per_burst, positions, screen_x_uv, screen_y_uv);
       // last_mesh.setup_vertex_and_instancing_buffers_for_particles(n_instances, 8, 0.2);
       last_mesh.setup_vertex_and_instancing_buffers_for_particles(n_instances);
       last_mesh.update_instancing_buffer_data_for_particles(last_particles);

       if (! do_tick_gone_diegos) {
          if (! tick_function_is_active()) {
             int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
             idle_function_spin_rock_token = new_tick_id;
          }
          do_tick_gone_diegos = true;
       }
    }
 }

// static
void
graphics_info_t::setup_draw_for_particles_for_gone_diff_map_peaks(const std::vector<std::pair<glm::vec3, float> > &positions) {

   // std::cout << "********* setup_draw_for_particles_for_gone_diff_map_peaks() " << positions.size() << std::endl;
   play_sound("diff-map-peak-gone-pop");

   glm::vec3 screen_x_uv = get_screen_x_uv();
   glm::vec3 screen_y_uv = get_screen_y_uv();
   particle_container_t &particles = meshed_particles_for_gone_diff_map_peaks.particle_container;
   Mesh &mesh                      = meshed_particles_for_gone_diff_map_peaks.mesh;
   attach_buffers();
   unsigned int n_particles_per_burst = 5;
   int n_instances = n_particles_per_burst * positions.size();
   particles.make_gone_diff_map_peaks_particles(n_particles_per_burst, positions, screen_x_uv, screen_y_uv);
   mesh.setup_vertex_and_instancing_buffers_for_particles(n_instances); // wrong polygon.
   mesh.clear();
   mesh.setup_camera_facing_polygon(8, 0.1, false, 0);
   mesh.update_instancing_buffer_data_for_particles(particles);

   if (! do_tick_gone_diff_map_peaks) {
      if (! tick_function_is_active()) {
         int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
         idle_function_spin_rock_token = new_tick_id;
      }
      do_tick_gone_diff_map_peaks = true;
   }
}


// static
void
graphics_info_t::setup_draw_for_bad_nbc_atom_pair_markers() {

   // 20221005-PE this causes an error on startup (like hud buttons and hud geometry bars)
   //             and setup_draw_for_happy_face_residue_markers_init()
   //             We are already in the correct framebuffer.
   //
   // attach_buffers();
   // GLenum err = glGetError();
   // if (err)
   //    std::cout << "GL ERROR:: start of setup_draw_bad_nbc_atom_pair_markers() "  << err << std::endl;

   texture_for_bad_nbc_atom_pair_markers.init("angry-diego.png");
   float ts = 0.7; // relative texture size
   tmesh_for_bad_nbc_atom_pair_markers.setup_camera_facing_quad(ts, ts, 0.0, 0.7);
   tmesh_for_bad_nbc_atom_pair_markers.setup_instancing_buffers(200);
   tmesh_for_bad_nbc_atom_pair_markers.draw_this_mesh = true;

}

void graphics_info_t::setup_draw_for_bad_nbc_atom_pair_dashed_line() {

   // the mesh gets setup on update
   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &vertices) {
      std::vector<s_generic_vertex> v_out(vertices.size());
      for (unsigned int i=0; i<vertices.size(); i++) {
         const auto &v = vertices[i];
         v_out[i].pos    = v.pos;
         v_out[i].normal = v.normal;
         v_out[i].color  = v.color;
      }
      return v_out;
   };

   GLenum err = glGetError();
   if (err)
      logger.log(log_t::WARNING, logging::function_name_t("setup_draw_for_bad_nbc_atom_pair_dashed_line"),
                 "---start---", stringify_error_code(err));

   attach_buffers();
   cylinder c;
   c.init_unit(20);
   c.add_flat_end_cap();
   c.add_flat_start_cap();
   bad_nbc_atom_pair_dashed_line = Mesh(convert_vertices(c.vertices), c.triangles);
   bad_nbc_atom_pair_dashed_line.set_name("bad_nbc_atom_pair_dashed_line Mesh");
   bad_nbc_atom_pair_dashed_line.setup_buffers();

   if (err)
      logger.log(log_t::WARNING, logging::function_name_t("setup_draw_for_bad_nbc_atom_pair_dashed_line"),
                 "---end---", stringify_error_code(err));

}


// static
void
graphics_info_t::draw_bad_nbc_atom_pair_markers(stereo_eye_t eye, unsigned int pass_type) {

   if (curmudgeon_mode) return;

   if (draw_bad_nbc_atom_pair_markers_flag) {
      if (! bad_nbc_atom_pair_marker_positions.empty()) {

         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation = get_model_rotation();
         glm::vec4 bg_col(background_colour, 1.0);
         texture_for_bad_nbc_atom_pair_markers.Bind(0);

         if (pass_type == PASS_TYPE_STANDARD)
            tmesh_for_bad_nbc_atom_pair_markers.draw_instances(&shader_for_happy_face_residue_markers,
                                                               mvp, model_rotation, bg_col, perspective_projection_flag);

         if (pass_type == PASS_TYPE_SSAO) {
            GtkAllocation allocation;
            gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
            int w = allocation.width;
            int h = allocation.height;
            bool do_orthographic_projection = ! perspective_projection_flag;
            auto model_matrix = get_model_matrix();
            auto view_matrix = get_view_matrix();
            auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
            tmesh_for_bad_nbc_atom_pair_markers.draw_instances_for_ssao(&shader_for_happy_face_residue_markers_for_ssao,
                                                                        model_matrix, view_matrix, projection_matrix);
         }
      }
   }
}

void graphics_info_t::draw_bad_nbc_atom_pair_dashed_lines(stereo_eye_t eye, unsigned int pass_type) {

   if (curmudgeon_mode) return;

   if (draw_bad_nbc_atom_pair_markers_flag) {

      if (! bad_nbc_atom_pair_marker_positions.empty()) {
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation = get_model_rotation();
         glm::vec4 bg_col(background_colour, 1.0);

         bool transfered_colour_is_instanced = true;
         if (pass_type == PASS_TYPE_STANDARD)
            bad_nbc_atom_pair_dashed_line.draw_instanced(pass_type,
                                                         &shader_for_instanced_objects,
                                                         eye,
                                                         mvp,
                                                         model_rotation,
                                                         lights,
                                                         eye_position,
                                                         bg_col,
                                                         shader_do_depth_fog_flag,
                                                         transfered_colour_is_instanced);

      }
   }

}


void graphics_info_t::setup_draw_for_chiral_volume_outlier_markers() {

    texture_for_chiral_volume_outlier_markers.init("chiral-volume-outlier-marker.png");
    float ts = 0.7; // relative texture size
    tmesh_for_chiral_volume_outlier_markers.setup_camera_facing_quad(ts, ts, 0.0, 0.7);
    tmesh_for_chiral_volume_outlier_markers.setup_instancing_buffers(200);
    tmesh_for_chiral_volume_outlier_markers.draw_this_mesh = true;

}

// static
void graphics_info_t::draw_chiral_volume_outlier_markers(stereo_eye_t eye, unsigned int pass_type) {

   if (curmudgeon_mode) return;

   // unlike NBC markers, each molecule can have it's own chiral volume outlier markers
   for (unsigned int imol=0; imol<molecules.size(); imol++) {
       if (is_valid_model_molecule(imol)) {
          if (molecules[imol].draw_it) {
             if (molecules[imol].draw_chiral_volume_outlier_markers_flag) {
                if (! molecules[imol].chiral_volume_outlier_marker_positions.empty()) {

                   unsigned int n = graphics_info_t::molecules[imol].chiral_volume_outlier_marker_positions.size();

                   glm::mat4 mvp = get_molecule_mvp(eye);
                   glm::mat4 model_rotation = get_model_rotation();
                   glm::vec4 bg_col(background_colour, 1.0);
                   texture_for_chiral_volume_outlier_markers.Bind(0);

                   if (pass_type == PASS_TYPE_STANDARD) {
                      tmesh_for_chiral_volume_outlier_markers.draw_instances(&shader_for_happy_face_residue_markers,
                                                                             mvp, model_rotation, bg_col, perspective_projection_flag);
                   }

                   if (pass_type == PASS_TYPE_SSAO) {
                      GtkAllocation allocation;
                      gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
                      int w = allocation.width;
                      int h = allocation.height;
                      bool do_orthographic_projection = ! perspective_projection_flag;
                      auto model_matrix = get_model_matrix();
                      auto view_matrix = get_view_matrix();
                      auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
                      tmesh_for_chiral_volume_outlier_markers.draw_instances_for_ssao(&shader_for_happy_face_residue_markers_for_ssao,
                                                                                      model_matrix, view_matrix, projection_matrix);
                   }
                }
             }
          }
       }
    }
 }

 //static
 void
    graphics_info_t::update_chiral_volume_outlier_marker_positions() {

    for (unsigned int imol=0; imol<molecules.size(); imol++) {
       if (is_valid_model_molecule(imol)) {
          if (molecules[imol].draw_chiral_volume_outlier_markers_flag) {
             unsigned int n_prev = molecules[imol].chiral_volume_outlier_marker_positions.size();
             molecules[imol].fill_chiral_volume_outlier_marker_positions(1);
             const auto &positions = molecules[imol].chiral_volume_outlier_marker_positions;
             if (positions.size() < n_prev) {
                play_sound("STARS");
             }
             if (! positions.empty()) {
                // update the instancing mesh
                attach_buffers();
                tmesh_for_chiral_volume_outlier_markers.draw_this_mesh = true;
                tmesh_for_chiral_volume_outlier_markers.update_instancing_buffer_data(positions);
                molecules[imol].draw_chiral_volume_outlier_markers_flag = true;
             }
          }
       }
    }
 }

void graphics_info_t::add_unhappy_atom_marker(int imol, const coot::atom_spec_t &atom_spec) {

   if (curmudgeon_mode) return;

   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *at = molecules[imol].get_atom(atom_spec);
      if (at) {
         glm::vec3 p(at->x, at->y, at->z);
         auto &positions = molecules[imol].unhappy_atom_marker_positions;
         positions.push_back(p);
         attach_buffers();
         tmesh_for_unhappy_atom_markers.draw_this_mesh = true;
         tmesh_for_unhappy_atom_markers.update_instancing_buffer_data(positions);
         unsigned int n_instances = tmesh_for_unhappy_atom_markers.get_n_instances();
         if (false)
            std::cout << "debug:: :::::::::::::::::::::::::::::::: add position " << glm::to_string(p)
                      << "  " << n_instances << std::endl;
      }
   }
}

void graphics_info_t::remove_all_unhappy_atom_markers() {

   for (unsigned int imol=0; imol<molecules.size(); imol++) {
      if (is_valid_model_molecule(imol)) {
         if (! molecules[imol].unhappy_atom_marker_positions.empty()) {
            molecules[imol].unhappy_atom_marker_positions.clear();
         }
      }
   }
   std::vector<glm::vec3> empty;
   tmesh_for_unhappy_atom_markers.draw_this_mesh = false;
   tmesh_for_unhappy_atom_markers.update_instancing_buffer_data(empty);
   graphics_draw();
}

void graphics_info_t::setup_draw_for_unhappy_atom_markers() {

   texture_for_unhappy_atom_markers.init("sad-santiago.png");
   float ts = 0.7; // relative texture size
   tmesh_for_unhappy_atom_markers.setup_camera_facing_quad(ts, ts, 0.0, 0.7);
   tmesh_for_unhappy_atom_markers.setup_instancing_buffers(200);
   tmesh_for_unhappy_atom_markers.draw_this_mesh = true;
}

// static
void graphics_info_t::draw_unhappy_atom_markers(stereo_eye_t eye, unsigned int pass_type) {

   if (curmudgeon_mode) return;

   for (unsigned int imol=0; imol<molecules.size(); imol++) {
      if (is_valid_model_molecule(imol)) {
         if (molecules[imol].draw_it) {
                   glm::mat4 mvp = get_molecule_mvp(eye);
                   glm::mat4 model_rotation = get_model_rotation();
                   glm::vec4 bg_col(background_colour, 1.0);
                   texture_for_unhappy_atom_markers.Bind(0);

                   if (pass_type == PASS_TYPE_STANDARD) {
                      tmesh_for_unhappy_atom_markers.draw_instances(&shader_for_happy_face_residue_markers,
                                                                    mvp, model_rotation, bg_col, perspective_projection_flag);
                   }

                   if (pass_type == PASS_TYPE_SSAO) {
                      GtkAllocation allocation;
                      gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
                      int w = allocation.width;
                      int h = allocation.height;
                      bool do_orthographic_projection = ! perspective_projection_flag;
                      auto model_matrix = get_model_matrix();
                      auto view_matrix = get_view_matrix();
                      auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
                      tmesh_for_unhappy_atom_markers.draw_instances_for_ssao(&shader_for_happy_face_residue_markers_for_ssao,
                                                                             model_matrix, view_matrix, projection_matrix);
                   }
         }
      }
   }

}


// static
void
graphics_info_t::update_hydrogen_bond_positions() {

   auto atom_to_glm = [] (mmdb::Atom *at) { return glm::vec3(at->x, at->y, at->z); };

   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         coot::refinement_results_t &rr = saved_dragged_refinement_results;
         if (! rr.hydrogen_bond_atom_index_vec.empty()) {

            // fill hydrogen_bonds_atom_position_pairs
            unsigned int n = rr.hydrogen_bond_atom_index_vec.size();
            hydrogen_bonds_atom_position_pairs.clear();
            hydrogen_bonds_atom_position_pairs.reserve(n);

            for (unsigned int i=0; i<rr.hydrogen_bond_atom_index_vec.size(); i++) {
               const int &idx_1 = rr.hydrogen_bond_atom_index_vec[i].first;
               const int &idx_2 = rr.hydrogen_bond_atom_index_vec[i].second;
               if (idx_1 < moving_atoms_asc->n_selected_atoms) {
                  if (idx_2 < moving_atoms_asc->n_selected_atoms) {
                     mmdb::Atom *at_1 = moving_atoms_asc->atom_selection[idx_1];
                     mmdb::Atom *at_2 = moving_atoms_asc->atom_selection[idx_2];
                     glm::vec3 p_1 = atom_to_glm(at_1);
                     glm::vec3 p_2 = atom_to_glm(at_2);
                     hydrogen_bonds_atom_position_pairs.push_back(std::make_pair(p_1, p_2));
                  }
               }
            }

            attach_buffers();
            std::string label = "Hydrogen Bonds";
            update_hydrogen_bond_mesh(label);
         }
      }
   }
}


//static
gboolean
graphics_info_t::wait_for_hooray_refinement_tick_func(GtkWidget *widget,
                                                      GdkFrameClock *frame_clock,
                                                      gpointer data) {
   gboolean continue_status = 1;

   if (setup_draw_for_particles_semaphore) {
      if (! particles_have_been_shown_already_for_this_round_flag) {
         graphics_info_t g;
         g.setup_draw_for_particles();
         setup_draw_for_particles_semaphore = false; // it's done it's job
         particles_have_been_shown_already_for_this_round_flag = true; // only once per round
         continue_status = 0; // job done.
      }
   }
   return continue_status;
}


void
graphics_info_t::setup_draw_for_boids() {

   if (boids.size() == 0) {
      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?

      unsigned int n_boids = 30;
      boids.make_boids(n_boids);

      meshed_generic_display_object m;
      coot::colour_holder col(0.4, 0.5, 0.6);
      std::pair<glm::vec3, glm::vec3> start_end(glm::vec3(1.95,0,0), glm::vec3(-1.95,0,0));
      m.add_cone(start_end, col, 3.0, 0.0, 24, false, true,
                 meshed_generic_display_object::FLAT_CAP,
                 meshed_generic_display_object::FLAT_CAP);
      mesh_for_boids = m.mesh;

      std::vector<glm::mat4>    mats(n_boids);
      std::vector<glm::vec4> colours(n_boids);
      for (unsigned int i=0; i<n_boids; i++) {
         mats[i] = glm::mat4(1.0f);
         colours[i] = glm::vec4(0.2, 0.6, 0.4, 1.0);
      }
      Material material;
      mesh_for_boids.setup_rtsc_instancing(&shader_for_instanced_objects,
                                           mats, colours, n_boids, material);

      if (! tick_function_is_active()) {
         int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
      }
      do_tick_boids = true;

      // boids box

      std::vector<s_generic_vertex> vertices;
      std::vector<unsigned int> indices;
      glm::vec3 n(0,0,1);
      glm::vec4 c(0.3f, 0.3f, 0.3f, 1.0f);
      float boids_box_lim = boids.boids_box_limit;

      float corners[8][3] = {
                          {0,0,0}, //0
                          {0,0,1}, //1
                          {0,1,0}, //2
                          {0,1,1}, //3
                          {1,0,0}, //4
                          {1,0,1}, //5
                          {1,1,0}, //6
                          {1,1,1}};//7
      for (unsigned int i=0; i<8; i++) {
         for (unsigned int j=0; j<3; j++) {
            corners[i][j] *= 2.0f;
            corners[i][j] -= 1.0;
            corners[i][j] *= boids_box_lim;
         }
      }
      for (unsigned int ii=0; ii<8; ii++)
         vertices.push_back(s_generic_vertex(glm::vec3(corners[ii][0],corners[ii][1],corners[ii][2]), n, c));

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

      lines_mesh_for_boids_box = LinesMesh(vertices, indices);
      lines_mesh_for_boids_box.setup();
   }
}

void
graphics_info_t::draw_hud_ligand_view() {

   // 20230618-PE we don't want to do anything here if graphics_ligand_view_flag is false.
   // graphics_ligand_view_flag is false when we are not centred on a ligand.

   if (graphics_ligand_view_flag) {
      if  (is_valid_model_molecule(graphics_ligand_view_imol)) {
         if (molecules[graphics_ligand_view_imol].is_displayed_p()) {

            GtkAllocation allocation;
            gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
            float w = allocation.width;
            float h = allocation.height;
            GLenum err = glGetError();
            if (err)
               std::cout << "draw_ligand_view() --- start --- " << err << std::endl;

            graphics_ligand_mesh_molecule.draw(&shader_for_ligand_view,
                                               &shader_for_hud_geometry_tooltip_text,
                                               w, h, ft_characters);
            err = glGetError();
            if (err)
               std::cout << "GL ERROR:: draw_ligand_view() --- end --- " << err << std::endl;
         }
      }
   }
}



void
graphics_info_t::draw_boids(stereo_eye_t eye) {

   if (boids.size() > 0) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      auto ccrc = RotationCentre();
      glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
      mesh_for_boids.draw(&shader_for_instanced_objects,
                          eye, mvp, model_rotation_matrix, lights, eye_position, rc, opacity, bg_col,
                          wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

      lines_mesh_for_boids_box.draw(&shader_for_lines, mvp, model_rotation_matrix);
   }
}

void
graphics_info_t::update_hydrogen_bond_mesh(const std::string &label) {

   return; // 20230823-PE for now we don't want hydrogen-bond mesh, I don't want
           // to debug while the tick function is active.

   // caller fills static std::vector<std::pair<glm::vec3, glm::vec3> > hydrogen_bonds_atom_position_pairs
   // before this function

   Material material;
   material.shininess = 10.0;
   material.specular_strength = 0.02;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   Mesh mesh(label);
   mesh_for_hydrogen_bonds = mesh;
   Shader &shader = shader_for_instanced_objects;
   mesh_for_hydrogen_bonds.setup_hydrogen_bond_cyclinders(&shader, material);

   std::chrono::time_point<std::chrono::high_resolution_clock> tp_now = std::chrono::high_resolution_clock::now();
   std::chrono::time_point<std::chrono::high_resolution_clock> tp_prev = tick_hydrogen_bond_mesh_t_previous;
   auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - tp_prev);
   float theta = 0.002 * delta.count();
   // std::cout << "delta from time " << delta.count() << " theta " << theta << std::endl;
   std::vector<glm::mat4> mats;
   for (unsigned int i=0; i<hydrogen_bonds_atom_position_pairs.size(); i++) {
      const std::pair<glm::vec3, glm::vec3> &p = hydrogen_bonds_atom_position_pairs[i];
      mats.push_back(Mesh::make_hydrogen_bond_cylinder_orientation(p.first, p.second, theta));
   }
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // Needed? Yes! Vital
   mesh_for_hydrogen_bonds.update_instancing_buffer_data_standard(mats);
   add_a_tick();
   do_tick_hydrogen_bonds_mesh = true;

}

void
graphics_info_t::draw_hydrogen_bonds_mesh(stereo_eye_t eye) {

   // 20210827-PE  each molecule should have its own hydrogen bond mesh. Not just one of them.
   // Fix that later.

   if (mesh_for_hydrogen_bonds.get_draw_this_mesh()) {

      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);

      int pass_type = PASS_TYPE_STANDARD;
      mesh_for_hydrogen_bonds.draw_instanced(pass_type,
                                             &shader_for_instanced_objects, eye,
                                             mvp, model_rotation_matrix, lights, eye_position, bg_col,
                                             shader_do_depth_fog_flag, false, false, true, 0, 0, 0, 0.2);
   }

}


std::vector<glm::vec3>
graphics_info_t::residue_to_positions(mmdb::Residue *residue_p) const {
   std::vector<glm::vec3> v;
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for(int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         glm::vec3 p(at->x, at->y, at->z);
         v.push_back(p);
      }
   }
   return v;
};

#include "pulse-data.hh"

void
graphics_info_t::setup_delete_item_pulse(mmdb::Residue *residue_p) {

   // next you use this functionn make it a member of graphics_info_t
   // gboolean delete_item_pulse_func(GtkWidget *widget,
   //                                 GdkFrameClock *frame_clock,
   //                                 gpointer data)
   //
   auto delete_item_pulse_func = [] (GtkWidget *widget,
                                     GdkFrameClock *frame_clock,
                                     gpointer data) {

                                    gboolean continue_status = 1;
                                    pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
                                    pulse_data->n_pulse_steps += 1;
                                    if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
                                       continue_status = 0;
                                       lines_mesh_for_generic_pulse.clear();
                                       generic_pulse_centres.clear();
                                    } else {
                                       float ns = pulse_data->n_pulse_steps;
                                       lines_mesh_for_generic_pulse.update_buffers_for_pulse(ns, -1);
                                    }
                                    graphics_draw();
                                    return gboolean(continue_status);
                                 };

   pulse_data_t *pulse_data = new pulse_data_t(0, 20); // 20 matches the number in update_buffers_for_pulse()
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> positions = residue_to_positions(residue_p);
   generic_pulse_centres = positions;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = true;
   float radius_overall = 6.0;
   unsigned int n_rings = 3;
   lines_mesh_for_generic_pulse.setup_red_pulse(radius_overall, n_rings, broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], delete_item_pulse_func, user_data, NULL);

};

void
graphics_info_t::setup_delete_residues_pulse(const std::vector<mmdb::Residue *> &residues) {

   // next you use this function make it a member of graphics_info_t
   // gboolean delete_item_pulse_func(GtkWidget *widget,
   //                                 GdkFrameClock *frame_clock,
   //                                 gpointer data)
   //
   auto delete_item_pulse_func = [] (GtkWidget *widget,
                                     GdkFrameClock *frame_clock,
                                     gpointer data) {

                                    gboolean continue_status = 1;
                                    pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
                                    pulse_data->n_pulse_steps += 1;
                                    if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
                                       continue_status = 0;
                                       lines_mesh_for_generic_pulse.clear();
                                       generic_pulse_centres.clear();
                                    } else {
                                       float ns = pulse_data->n_pulse_steps;
                                       lines_mesh_for_generic_pulse.update_buffers_for_pulse(ns, -1);
                                    }
                                    graphics_draw();
                                    return gboolean(continue_status);
                                 };

   pulse_data_t *pulse_data = new pulse_data_t(0, 20); // 20 matches the number in update_buffers_for_pulse()
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> all_positions;
   for (unsigned int i=0; i<residues.size(); i++) {
      mmdb::Residue *residue_p = residues[i];
      std::vector<glm::vec3> residue_positions = residue_to_positions(residue_p);
      all_positions.insert(all_positions.end(), residue_positions.begin(), residue_positions.end());
   }
   generic_pulse_centres = all_positions;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = true;
   unsigned int n_rings = 3;
   float radius_overall = 6.0;
   lines_mesh_for_generic_pulse.setup_red_pulse(radius_overall, n_rings, broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], delete_item_pulse_func, user_data, NULL);

};


// static
gboolean
graphics_info_t::invalid_residue_pulse_function(GtkWidget *widget,
                                                GdkFrameClock *frame_clock,
                                                gpointer data) {

   gboolean continue_status = 1;
   pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
   pulse_data->n_pulse_steps += 1;
   if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
      continue_status = 0;
      lines_mesh_for_identification_pulse.clear();
      generic_pulse_centres.clear(); // we sneakily use this vector (but no longer)
   } else {
      float ns = pulse_data->n_pulse_steps;
      lines_mesh_for_identification_pulse.update_buffers_for_invalid_residue_pulse(ns);
   }
   graphics_draw();
   return gboolean(continue_status);
}

// static
gboolean
graphics_info_t::screen_centre_pulse_function(GtkWidget *widget,
                                              GdkFrameClock *frame_clock,
                                              gpointer data) {

   gboolean continue_status = 1;
   pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
   pulse_data->n_pulse_steps += 1;
   if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
      continue_status = 0;
      lines_mesh_for_identification_pulse.clear();
      generic_pulse_centres.clear();
   }
   return gboolean(continue_status);
}

// static
gboolean
graphics_info_t::generic_pulse_function(GtkWidget *widget,
                                        GdkFrameClock *frame_clock,
                                        gpointer data) {

   gboolean continue_status = 1;
   pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
   pulse_data->n_pulse_steps += 1;
   if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
      continue_status = 0;
      lines_mesh_for_generic_pulse.clear();
      generic_pulse_centres.clear();
   } else {
      lines_mesh_for_generic_pulse.update_buffers_by_resize(pulse_data->resize_factor);
   }
   graphics_draw();
   return gboolean(continue_status);
}


void
graphics_info_t::setup_invalid_residue_pulse(mmdb::Residue *residue_p) {

   pulse_data_t *pulse_data = new pulse_data_t(0, 24);
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> residue_positions = residue_to_positions(residue_p);
   generic_pulse_centres = residue_positions; // sneakily use a wrongly named function
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   unsigned int n_rings = 3;
   float radius_overall = 6.0;
   bool broken_line_mode = false;
   lines_mesh_for_identification_pulse.setup_red_pulse(radius_overall, n_rings, broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], invalid_residue_pulse_function, user_data, NULL);

}


void graphics_info_t::draw_at_screen_centre_pulse(stereo_eye_t eye) {

   // 2025-12-06-PE identification and screen-centre are the same thing. "identification" should be renamed.

   if (! lines_mesh_for_identification_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation_matrix = get_model_rotation();
      myglLineWidth(2.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_at_screen_centre_pulse() post glLineWidth " << err << std::endl;
      lines_mesh_for_identification_pulse.draw(&shader_for_lines_pulse,
                                               identification_pulse_centre,
                                               mvp, model_rotation_matrix, true);
   }
}

void graphics_info_t::draw_generic_pulses(stereo_eye_t eye) {

   if (false)
      std::cout << "draw_generic_pulses()  -- start -- "
                << lines_mesh_for_generic_pulse.empty() << " " << generic_pulse_centres.size() << std::endl;

   if (! lines_mesh_for_generic_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation_matrix = get_model_rotation();
      for (auto pulse_centre : generic_pulse_centres)
         lines_mesh_for_generic_pulse.draw(&shader_for_lines_pulse,
                                           pulse_centre, mvp,
                                           model_rotation_matrix, true);
   }
}

void graphics_info_t::draw_invalid_residue_pulse(stereo_eye_t eye) {

   return; // because we do it in draw_generic_pulses()

   if (! lines_mesh_for_generic_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation_matrix = get_model_rotation();
      myglLineWidth(3.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_invalid_residue_pulse() glLineWidth " << err << std::endl;
      for (auto pulse_centre : generic_pulse_centres)
         lines_mesh_for_generic_pulse.draw(&shader_for_lines_pulse,
                                           pulse_centre, mvp,
                                           model_rotation_matrix, true);
   }
}

void
graphics_info_t::draw_delete_item_pulse(stereo_eye_t eye) {

   if (! lines_mesh_for_generic_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp(eye);
      glm::mat4 model_rotation_matrix = get_model_rotation();
      myglLineWidth(2.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_delete_item_pulse() glLineWidth " << err << std::endl;
      for (unsigned int i=0; i<generic_pulse_centres.size(); i++) {
         lines_mesh_for_generic_pulse.draw(&shader_for_lines_pulse,
                                           generic_pulse_centres[i],
                                           mvp, model_rotation_matrix, true);
      }
   }
}

void
graphics_info_t::draw_pointer_distances_objects(stereo_eye_t eye) {

   if (show_pointer_distances_flag) {
      if (! pointer_distances_object_vec.empty()) {
         Shader &shader = shader_for_moleculestotriangles;
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation_matrix = get_model_rotation();
         glm::vec4 bg_col(background_colour, 1.0);
         bool show_just_shadows = false;
         bool wireframe_mode = false;
         float opacity = 1.0f;
         auto ccrc = RotationCentre();
         glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
         mesh_for_pointer_distances.mesh.draw(&shader, eye, mvp, model_rotation_matrix, lights, eye_position, rc, opacity,
                                              bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

         if (! labels_for_pointer_distances.empty()) {
            Shader &shader_labels = shader_for_atom_labels;
            for (unsigned int i=0; i<labels_for_pointer_distances.size(); i++) {
               const auto &label = labels_for_pointer_distances[i];
               tmesh_for_labels.draw_atom_label(label.label, label.position, label.colour, &shader_labels,
                                                eye, mvp, model_rotation_matrix, bg_col,
                                                shader_do_depth_fog_flag, perspective_projection_flag);
            }
         }
      }
   }
}

void
graphics_info_t::draw_translation_gizmo(stereo_eye_t eye) { // maybe rotation gizmo too, later.

   if (translation_gizmo_mesh.get_draw_this_mesh()) {
      bool do_it = false;
      int tgagdo_num = translation_gizmo.attached_to_generic_display_object_number;
      int tgam_num   = translation_gizmo.attached_to_molecule_number;
      if (tgagdo_num != translation_gizmo_t::UNATTACHED)
         if (tgagdo_num >= 0)
            if (tgagdo_num < int(generic_display_objects.size()))
               if (generic_display_objects[tgagdo_num].mesh.get_draw_this_mesh())
                  do_it = true;
      if (is_valid_model_molecule(tgam_num))
         if (molecules[tgagdo_num].get_mol_is_displayed())
            do_it = true;
      if (is_valid_map_molecule(tgam_num))
         if (molecules[tgagdo_num].is_displayed_p())
            do_it = true;
      if (do_it) {
         Shader &shader = shader_for_moleculestotriangles;
         glm::mat4 mvp = get_molecule_mvp(eye);
         glm::mat4 model_rotation_matrix = get_model_rotation();
         glm::vec4 bg_col(background_colour, 1.0);
         bool show_just_shadows = false;
         bool wireframe_mode = false;
         float opacity = 1.0f;
         auto ccrc = RotationCentre();
         glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());
         translation_gizmo_mesh.draw(&shader, eye, mvp, model_rotation_matrix, lights, eye_position, rc, opacity,
                                     bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);
      }
   }

}



void
graphics_info_t::make_extra_distance_restraints_objects() {

   // c.f. update_hydrogen_bond_mesh().

   double penalty_min = 0.05; // only restraints that have more than this "distortion" are considered for drawing.
                               // Make this user-setable.

   // the model has been updated, we need to update the positions and orientations using in the instancing

   auto clipper_to_glm = [] (const clipper::Coord_orth &co) {
                            return glm::vec3(co.x(), co.y(), co.z());
                         };

   // How frequently should this function be called? Every frame? Hmm..
   if (moving_atoms_extra_restraints_representation.bonds.empty()) return;

   unsigned int maerrb_size = moving_atoms_extra_restraints_representation.bonds.size();
   attach_buffers();
   mesh_for_extra_distance_restraints.setup_instancing_buffer_data_for_extra_distance_restraints(maerrb_size);
   // now fill extra_distance_restraints_markup_data

   extra_distance_restraints_markup_data.clear();
   extra_distance_restraints_markup_data.reserve(moving_atoms_extra_restraints_representation.bonds.size());
   for (unsigned int i=0; i<moving_atoms_extra_restraints_representation.bonds.size(); i++) {
      const coot::extra_restraints_representation_t::extra_bond_restraints_respresentation_t &ebrr =
         moving_atoms_extra_restraints_representation.bonds[i];
      double dd = clipper::Coord_orth(ebrr.first - ebrr.second).lengthsq();
      double d = std::sqrt(dd);
      extra_distance_restraint_markup_instancing_data_t edrmid;
      // the width should represent the pulling power (i.e. the size of the penalty/distortion)
      // make a function extra_bond_restraints_respresentation_t::get_penalty(alpha, sigma);
      double sigma = 0.1; // what is this actually?
      double penalty = ebrr.distortion_score_GM(sigma, geman_mcclure_alpha);
      if (penalty < penalty_min) continue;
      double width = 0.2 * penalty; // 20230823-PE was 0.3
      if (width < 0.01) width = 0.01;
      if (width > 0.10) width = 0.10;
      edrmid.width = width;
      edrmid.length = static_cast<float>(d);
      edrmid.position = clipper_to_glm(ebrr.second);

      clipper::Coord_orth delta = ebrr.second - ebrr.first;
      clipper::Coord_orth delta_uv = clipper::Coord_orth(delta.unit());
      glm::vec3 delta_uv_glm = clipper_to_glm(delta_uv);

      glm::mat4 ori44 = glm::orientation(delta_uv_glm, glm::vec3(0.0, 0.0, 1.0));
      glm::mat3 ori33 = glm::mat3(ori44);
      edrmid.orientation = ori33;

      // std::cout << "edrmid " << i << " position " << glm::to_string(edrmid.position) << " length " << d
      // << "ori " << glm::to_string(edrmid.orientation) << std::endl;

      double delta_length = ebrr.length_delta();
      glm::vec4 colour_base = glm::vec4(0.5f, 0.5f, 0.5f, 1.0f);
      // std::cout << "delta length " << delta_length << std::endl;

      // for colouring, limit the delta_length)
      if (delta_length >  1.0) delta_length =  1.0;
      if (delta_length < -1.0) delta_length = -1.0;
      glm::vec4 colour = colour_base + static_cast<float>(delta_length) * glm::vec4(-0.8f, 0.8f, -0.8f, 0.0f);
      edrmid.colour = 0.8f * colour;
      extra_distance_restraints_markup_data.push_back(edrmid);
   }

   if (false) { // 20230519-PE come back to this.
      std::cout << "in make_extra_distance_restraints_objects() bond size "
                << moving_atoms_extra_restraints_representation.bonds.size() << std::endl;
      std::cout << "in make_extra_distance_restraints_objects() extra_distance_restraints_markup_data size "
                << extra_distance_restraints_markup_data.size() << std::endl;
   }
   mesh_for_extra_distance_restraints.update_instancing_buffer_data_for_extra_distance_restraints(extra_distance_restraints_markup_data);

}

// static
void
graphics_info_t::draw_extra_distance_restraints(stereo_eye_t eye, int pass_type) {

   // 20230825-PE we don't want to see these if there are no intermediate atoms being displayed
   // Maybe they should be cleared up on "clear_moving_atoms()" (or whatever the function is called).

   if (!moving_atoms_asc)
      return;
   if (!moving_atoms_asc->mol)
      return;

   // 20231121-PE HACK for now:
   if (pass_type == PASS_TYPE_WITH_SHADOWS) pass_type = PASS_TYPE_STANDARD;

   if (! draw_it_for_moving_atoms_restraints_graphics_object_user_control) return;

   // std::cout << "draw_extra_distance_restraints() pass_type: " << pass_type << std::endl;
   // std::cout << "draw_extra_distance_restraints() mesh_for_extra_distance_restraints " << std::endl;;
   // mesh_for_extra_distance_restraints.debug();

   // it used to be called draw_it_for_moving_atoms_restraints_graphics_object - why not use that varible?
   //
   // what about draw_it_for_moving_atoms_restraints_graphics_object?
   //
   if (pass_type == PASS_TYPE_STANDARD) {
      if (false)
         std::cout << "in draw_extra_distance_restraints() with show_extra_distance_restraints_flag "
                   << show_extra_distance_restraints_flag << " " << extra_distance_restraints_markup_data.size()
                   << std::endl;
      if (show_extra_distance_restraints_flag) {
         if (! extra_distance_restraints_markup_data.empty()) {
            glm::mat4 mvp = get_molecule_mvp(eye);
            glm::mat4 model_rotation_matrix = get_model_rotation();
            glm::vec4 bg_col(background_colour, 1.0f);
            glDisable(GL_BLEND);
            Shader &shader = shader_for_extra_distance_restraints;
            mesh_for_extra_distance_restraints.draw_extra_distance_restraint_instances(&shader, mvp, model_rotation_matrix, lights,
                                                                                       eye_position, bg_col, shader_do_depth_fog_flag);
         }
      }
   }

   if (pass_type == PASS_TYPE_SSAO) {
      if (show_extra_distance_restraints_flag) {
         if (! extra_distance_restraints_markup_data.empty()) {
            Shader &shader = shader_for_extra_distance_restraints; // wrong shader - needs a new one.
            GtkAllocation allocation;
            gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
            int w = allocation.width;
            int h = allocation.height;
            bool do_orthographic_projection = ! perspective_projection_flag;
            auto model_matrix = get_model_matrix();
            auto view_matrix = get_view_matrix();
            auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
            mesh_for_extra_distance_restraints.draw_instances_for_ssao(&shader, model_matrix, view_matrix, projection_matrix);
         }
      }
   }

}

// called from gl widget realize function
// static
 void
    graphics_info_t::setup_lines_mesh_for_proportional_editing() {

    unsigned int n_points = 100;
    std::vector<s_generic_vertex> vertices(n_points);
    glm::vec3 n(0,0,1);
    glm::vec4 c(0.7,0.7,0.7,1.0);
    double r = 0.001;
    for (unsigned int i=0; i<n_points; i++) {
       double theta = 2.0 * M_PI * static_cast<double>(i) / 100.0;
       glm::vec3 pt(r * cos(theta), r * sin(theta), 0.0);
       vertices[i] = s_generic_vertex(pt, n, c);
    }

    std::vector<unsigned int> indices;
    for (unsigned int i=0; i<n_points; i++) {
       unsigned int i_next = i+1;
       if (i_next == n_points) i_next = 0;
       indices.push_back(i);
       indices.push_back(i_next);
    }

    lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring = LinesMesh(vertices, indices);
    std::string name = "lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring";
    lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring.set_name(name);
    lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring.setup();
 }



void
graphics_info_t::move_forwards() {
   // these are the other way round in perspective - that's interesting.
   translate_in_screen_z(3.0);
}

void
graphics_info_t::move_backwards() {
   translate_in_screen_z(-3.0);
}

void
graphics_info_t::step_screen_left() {
   translate_in_screen_x(-1.0);  // function uses zoom
}

void
graphics_info_t::step_screen_right() {
   translate_in_screen_x(1.0);
}

#include <glm/gtx/rotate_vector.hpp>
#include "matrix-utils.hh"

// widget is the glarea.
//
gint
graphics_info_t::idle_contour_function(gpointer data) {

   gint continue_status = 0;
   bool something_changed = false;

   bool is_from_contour_level_change(GPOINTER_TO_INT(data));

   // when there's nothing else to do, update the contour levels
   //
   // then update maps

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_xmap()) { // FIXME or nxmap : needs test for being a map molecule
         int &cc = graphics_info_t::molecules[imol].pending_contour_level_change_count;

         if (cc != 0) {

	          if (cc < 0) {
	             while (cc != 0) {
	                cc++;
	                graphics_info_t::molecules[imol].change_contour(-1);
	             }
	          }

	          if (cc > 0) {
	              while (cc != 0) {
	                 cc--;
	                 graphics_info_t::molecules[imol].change_contour(1);
	              }
	          }

           graphics_info_t g;
           bool really_change_the_map_contours = true;
           if (! is_from_contour_level_change) really_change_the_map_contours = false;
	   g.molecules[imol].update_map(really_change_the_map_contours);
           float map_rmsd = g.molecules[imol].map_sigma();
	   continue_status = 0;
           float cl = g.molecules[imol].contour_level;
           float r = cl/map_rmsd;
           // std::cout << "DEBUG:: idle_contour_function() imol: " << imol << " contour level: "
	   // << g.molecules[imol].contour_level << " n-rmsd: " << r << std::endl;
	   logger.log(log_t::DEBUG, logging::function_name_t("idle_contour_function"),
		      {imol, "contour_level", g.molecules[imol].contour_level, "n-rmsd:",
		       r});
           g.set_density_level_string(imol, g.molecules[imol].contour_level);
           std::string s = "Map " + std::to_string(imol) + "  contour_level " +
              coot::util::float_to_string_using_dec_pl(cl, 3) + "  n-rmsd: " +
              coot::util::float_to_string_using_dec_pl(r, 3);
           add_status_bar_text(s.c_str());
           g.display_density_level_this_image = 1;
           something_changed = true;
         }
      }
   }

   if (something_changed)
      graphics_draw();

   // std::cout << "--- debug:: idle_contour_function() done " << continue_status << std::endl;
   return continue_status;
}


void
graphics_info_t::contour_level_scroll_scrollable_map(int direction) {

   int imol_scroll = scroll_wheel_map;
   if (! is_valid_map_molecule(imol_scroll)) {

      std::vector<int> dm = displayed_map_imols();
      if (std::find(dm.begin(), dm.end(), imol_scroll) == dm.end()) {
         if (dm.size() > 0)
            imol_scroll = dm[0];
      }
   }

   if (is_valid_map_molecule(imol_scroll)) {
      if (! molecules[imol_scroll].is_displayed_p()) {
         // don't scroll the map if the map is not displayed. Scroll the
         // map that *is* displayed
         std::vector<int> dm = displayed_map_imols();
         if (dm.size() > 0)
            imol_scroll = dm[0];
      }
   }

   if (is_valid_map_molecule(imol_scroll)) {
      // use direction
      if (direction ==  1) molecules[imol_scroll].pending_contour_level_change_count--;
      if (direction == -1) molecules[imol_scroll].pending_contour_level_change_count++;

      // std::cout << "INFO:: contour level for map " << imol_scroll << " is "
      //           << molecules[imol_scroll].contour_level
      //           << " pending: " << molecules[imol_scroll].pending_contour_level_change_count
      //           << std::endl;
      logger.log(log_t::INFO, "contour level for map", imol_scroll, "is",
		 molecules[imol_scroll].contour_level,
                 "step-size iso:", iso_level_increment,
                 "diff: ", diff_map_iso_level_increment,
                 "pending",
		 molecules[imol_scroll].pending_contour_level_change_count);

      set_density_level_string(imol_scroll, molecules[imol_scroll].contour_level);
      display_density_level_this_image = 1;

      graphics_draw(); // queue
   }
}

//static
void
graphics_info_t::fullscreen() {

   GtkWidget *window = graphics_info_t::get_main_window();

   if (GTK_IS_WINDOW(window)) {

      // GtkWidget *overlay    = widget_from_builder("main_window_graphics_overlay");
      GtkWidget *status_bar = widget_from_builder("main_window_statusbar");
      GtkWidget *tool_bar   = widget_from_builder("main_window_toolbar_hbox_outer");

      GtkWidget* sidebar = widget_from_builder("main_window_vbox_inner");

      gtk_widget_set_visible(tool_bar, FALSE);
      gtk_widget_set_visible(sidebar, FALSE);
      gtk_widget_set_visible(status_bar, FALSE);

      gtk_window_fullscreen(GTK_WINDOW(window));
      gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(window), FALSE);

      // gtk_box_remove(GTK_BOX(vbox), status_bar);
      // gtk_box_remove(GTK_BOX(vbox), tool_bar);

      gtk_widget_set_visible(status_bar, FALSE);
      gtk_widget_set_visible(tool_bar,   FALSE);

      graphics_info_t::add_status_bar_text(""); // clear it

      graphics_grab_focus();

   } else {
      g_error("%p is not a Gtk.Window !", window);
   }
}

//static
void
graphics_info_t::unfullscreen() {

   GtkWidget *window = graphics_info_t::get_main_window();
   if (GTK_IS_WINDOW(window)) {
      gtk_window_unfullscreen(GTK_WINDOW(window));
      gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(window),TRUE);

      GtkWidget* sidebar    = widget_from_builder("main_window_vbox_inner");
      GtkWidget* tool_bar   = widget_from_builder("main_window_toolbar_hbox_outer");
      GtkWidget* status_bar = widget_from_builder("main_window_statusbar");

      gtk_widget_set_visible(status_bar, TRUE);
      gtk_widget_set_visible(tool_bar,   TRUE);

      if(false) {
         GtkWidget *vbox       = widget_from_builder("main_window_vbox");
         GtkWidget *overlay    = widget_from_builder("main_window_graphics_overlay");

         gtk_overlay_remove_overlay(GTK_OVERLAY(overlay), tool_bar);
         gtk_box_append(GTK_BOX(vbox), tool_bar);

         gtk_overlay_remove_overlay(GTK_OVERLAY(overlay), status_bar);
         gtk_box_append(GTK_BOX(vbox), status_bar);
      }

      //gtk_widget_set_visible(tool_bar_frame, TRUE);
      gtk_widget_set_visible(tool_bar, TRUE);
      gtk_widget_set_visible(sidebar, TRUE);
      gtk_widget_set_visible(status_bar, TRUE);

   } else {
      g_error("%p is not a Gtk.Window !", window);
   }
}

#include "pumpkin.hh"
 void
    graphics_info_t::pumpkin() {

    std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > p1 = ::pumpkin();
    std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > p2 = ::pumpkin_stalk();

    glm::vec4 col_1(0.85, 0.45, 0.19, 1.0);
    glm::vec4 col_2(0.35, 0.45, 0.19, 1.0);
    attach_buffers();
    std::vector<s_generic_vertex> v1(p1.first.size());
    std::vector<s_generic_vertex> v2(p2.first.size());

    coot::Cartesian scc = get_rotation_centre_cart();
    glm::vec3 sc(scc.x(), scc.y(), scc.z());
    for (unsigned int i=0; i<p1.first.size(); i++)
       v1[i] = s_generic_vertex(2.0f * p1.first[i].pos + sc, p1.first[i].normal, col_1);
    for (unsigned int i=0; i<p2.first.size(); i++)
       v2[i] = s_generic_vertex(2.0f * p2.first[i].pos + sc, p2.first[i].normal, col_2);

    // add v2 to v1
    unsigned int idx_base = p1.first.size();
    unsigned int idx_tri_base = p1.second.size();
    v1.insert(v1.end(), v2.begin(), v2.end());
    std::vector<g_triangle> triangles = p1.second;
    triangles.insert(triangles.end(), p2.second.begin(), p2.second.end());
    for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
       triangles[i].rebase(idx_base);

    Mesh m(v1, triangles);
    m.set_name("Pumpkin");
    Material material;
    m.setup(material);
    meshed_generic_display_object obj(m);

    generic_display_objects.push_back(obj);

    graphics_draw();

 }
