
#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#include "compat/coot-sysdep.h"

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "c-interface.h" // for update_go_to_atom_from_current_position()
#include "globjects.h"
#include "graphics-info.h"

#include "draw.hh"
#include "draw-2.hh"
#include "framebuffer.hh"

#include "text-rendering-utils.hh"
#include "cc-interface-scripting.hh"
#include "cylinder-with-rotation-translation.hh"

#include "screendump-tga.hh"

enum {VIEW_CENTRAL_CUBE, ORIGIN_CUBE};


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
   if (err) std::cout << "init_blur_quads() err is " << err << std::endl;

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

   graphics_info_t::shader_for_central_cube.Use();
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
         GdkDisplay *display = gdk_display_get_default();
         GdkWindow *window = 0;
#if (GTK_MAJOR_VERSION < 4)
         window = gtk_widget_get_window(GTK_WIDGET(widget));
#endif
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
      }
   }
}


// static
void
graphics_info_t::mouse_zoom(double delta_x, double delta_y) {

   // Zooming
   double fx = 1.0 + delta_x/300.0;
   double fy = 1.0 + delta_y/300.0;
   if (fx > 0.0) graphics_info_t::zoom /= fx;
   if (fy > 0.0) graphics_info_t::zoom /= fy;
   if (false)
      std::cout << "zooming with perspective_projection_flag "
                << graphics_info_t::perspective_projection_flag
                << " " << graphics_info_t::zoom << std::endl;
   if (! graphics_info_t::perspective_projection_flag) {
      // std::cout << "now zoom: " << g.zoom << std::endl;
   } else {
      // Move the eye towards the rotation centre (don't move the rotation centre)
      if (fabs(delta_y) > fabs(delta_x))
         delta_x = delta_y;
      float sf = 1.0 - delta_x * 0.003;
      graphics_info_t::eye_position.z *= sf;

      { // own graphics_info_t function - c.f. adjust clipping
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
            std::cout << "on_glarea_motion_notify(): debug l: " << l << " post-manip: "
                      << graphics_info_t::screen_z_near_perspective << " "
                      << graphics_info_t::screen_z_far_perspective << std::endl;
      }
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
   // for both orthograph and perspective

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
      float sr = screen_ratio;
      GLfloat near =  -0.1 * zoom * clipping_front + eye_position.z;
      GLfloat far  =   0.3 * zoom * clipping_back  + eye_position.z;

      if (false)
         std::cout << "debug:: get_projection_matrix() near " << near << " far " << far
                   << " clipping-front: " << clipping_front << " clipping_back: " << clipping_back << " "
                   << "eye_position " << glm::to_string(eye_position)
                   << " zoom " << zoom << std::endl;

      glm::mat4 projection_matrix = glm::ortho(-0.3f*zoom*sr, 0.3f*zoom*sr,
                                               -0.3f*zoom,    0.3f*zoom,
                                               near, far);
      // std::cout << "projection matrix ortho " << glm::to_string(projection_matrix) << std::endl;
      return projection_matrix;
   } else {
      float fov = 30.0; // put this in graphics_info_t and add setters and getters to the API.
      glm::mat4 projection_matrix_persp = glm::perspective(glm::radians(fov),
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

   float box_size = shadow_box_size; // user setable, default 66.
   if (box_size < 0.0) box_size = 66.0;
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
graphics_info_t::get_molecule_mvp(bool debug_matrices) {

   int w = graphics_x_size;
   int h = graphics_y_size;

   if (true) {  // debug problematic matrices - get rid of this, make sure that it doesn't do anything                    
      GtkAllocation allocation;
      gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
      w = allocation.width;
      h = allocation.height;
   }

   bool do_orthographic_projection = ! perspective_projection_flag; // weird
   glm::mat4  view_matrix =       get_view_matrix();
   glm::mat4 model_matrix =      get_model_matrix();
   glm::mat4  proj_matrix = get_projection_matrix(do_orthographic_projection, w, h);
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

      glm::mat4 trackball_matrix = glm::toMat4(graphics_info_t::glm_quat);
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
   glm::mat4 mvp = get_molecule_mvp();
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

   return glm::toMat4(graphics_info_t::glm_quat);

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
graphics_info_t::draw_map_molecules(bool draw_transparent_maps) {

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: gtk3_draw_map_molecules() -- start -- " << err << std::endl;

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

      glm::mat4 mvp = get_molecule_mvp();
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
         m.draw_map_molecule(draw_transparent_maps, shader, mvp, model_rotation, eye_position, ep,
                             lights, background_colour, perspective_projection_flag);
      }
   }

   if (needs_blend_reset) {
      glDisable(GL_BLEND);
   }
}

void
graphics_info_t::draw_model_molecules() {

   glm::mat4 mvp = get_molecule_mvp();

   // std::cout << "debug:: mvp in draw_model_molecules() is     " << glm::to_string(mvp) << std::endl;
   glm::mat4 model_rotation = get_model_rotation();

   glm::vec4 bgc(background_colour, 1.0);

   for (int ii=n_molecules()-1; ii>=0; ii--) {
      if (! is_valid_model_molecule(ii)) continue;
      molecule_class_info_t &m = molecules[ii];
      if (! m.draw_it) continue;

      // I think that this is for the instanced meshes.
      // Shader &shader_p = shader_for_model_as_meshes; // is this actually instanced meshes?
      // m.draw_molecule_as_meshes(&shader_p, mvp, model_rotation, lights, eye_position, bgc, shader_do_depth_fog_flag);

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

      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      Shader *shader_p = &shader_for_meshes;
      if (m.draw_model_molecule_as_lines) {
         float lw = m.get_bond_thickness(); // returns an int.
         m.molecule_as_mesh.draw_simple_bond_lines(&shader_for_symmetry_atoms_bond_lines, mvp, bgc, lw, shader_do_depth_fog_flag);
      } else {
         m.molecule_as_mesh.draw(shader_p, mvp, model_rotation, lights, eye_position, opacity, bgc,
                                 wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);
      }

      m.draw_dots(&shader_for_rama_balls, mvp, model_rotation, lights, eye_position,
                     bgc, shader_do_depth_fog_flag);

      glEnable(GL_BLEND);
      draw_molecule_atom_labels(m, mvp, model_rotation);

   }
}

void
graphics_info_t::draw_model_molecules_with_shadows() {

   if (show_symmetry) {
      for (int ii=n_molecules()-1; ii>=0; ii--) {
         if (! is_valid_model_molecule(ii)) continue;
         molecule_class_info_t &m = molecules[ii];
         if (! m.draw_it) continue;
         Shader &symm_shader_p = shader_for_symmetry_atoms_bond_lines;
         glm::mat4 model_rotation = get_model_rotation();
         glm::vec4 bgc(background_colour, 1.0);
         glm::mat4 mvp = get_molecule_mvp();
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
   if (n_atoms_to_label == 0 && n_symm_atoms_to_label == 0) return;

   m.draw_atom_labels(brief_atom_labels_flag, seg_ids_in_atom_labels_flag,
                      label_colour, mvp, view_rotation);

   glDisable(GL_BLEND);

}

void
graphics_info_t::draw_intermediate_atoms(unsigned int pass_type) { // draw_moving_atoms()

   // this function gets called from draw_with_shadows() - but doesn't yet
   // use the shodows meshes shader.

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 model_rotation = get_model_rotation();


   molecule_class_info_t &m = graphics_info_t::moving_atoms_molecule;
   glm::vec4 bgc(background_colour, 1.0);
   bool show_just_shadows = false; // make this a member data item.
   if (false)
      std::cout << "----------------- draw_intermediate_atoms() intermediate atoms mesh has "
                << m.molecule_as_mesh.vertices.size() << " vertices and " << m.molecule_as_mesh.triangles.size()
                << " triangles" << std::endl;
   float opacity = 1.0f;

   if (pass_type == PASS_TYPE_STANDARD) {
      Shader &shader = shader_for_meshes_with_shadows;
      bool wireframe_mode = false;
      m.molecule_as_mesh.draw(&shader, mvp, model_rotation, lights, eye_position, opacity, bgc,
                              wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);
   }

   if (pass_type == PASS_TYPE_SSAO) {
      Shader &shader = shader_for_meshes_for_ssao;
      bool do_orthographic_projection = ! perspective_projection_flag;
      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(glareas[0]), &allocation);
      int w = allocation.width;
      int h = allocation.height;
      auto model_matrix = get_model_matrix();
      auto view_matrix = get_view_matrix();
      auto projection_matrix = get_projection_matrix(do_orthographic_projection, w, h);
      m.molecule_as_mesh.draw_for_ssao(&shader, model_matrix, view_matrix, projection_matrix);
   }

   if (pass_type == PASS_TYPE_FOR_SHADOWS) { // generating, not using - PASS_TYPE_GEN_SHADOW_MAP is a clearer name.
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
      m.molecule_as_mesh.draw(&shader, mvp_orthogonal, model_rotation, lights, dummy_eye_position,
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

   auto rr = saved_dragged_refinement_results;

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
graphics_info_t::draw_intermediate_atoms_rama_balls(unsigned int pass_type) {

   // 20220302-PE Currently I don't draw rama balls with instancing
   //             It would be nice to have.
   return;

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   Shader &shader = graphics_info_t::shader_for_rama_balls;

   glm::mat4 mvp = get_molecule_mvp();
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
      vertex_with_rotation_translation *vertices =
         new vertex_with_rotation_translation[n_vertices_for_atom_pull_restraints];
      vertex_with_rotation_translation *vertices_start = vertices;
      unsigned int iv = 0; // index into vertices - running

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
      GLuint n_bytes = sizeof(vertex_with_rotation_translation) * n_vertices_for_atom_pull_restraints;
      // maybe STATIC_DRAW, maybe not
      glBufferData(GL_ARRAY_BUFFER, n_bytes, vertices, GL_DYNAMIC_DRAW);
      err = glGetError();
      if (err) std::cout << "   error setup_atom_pull_restraints_glsl() E"
                          << " with GL err " << err << std::endl;


      glEnableVertexAttribArray(0);
      glEnableVertexAttribArray(1);
      glEnableVertexAttribArray(2);
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17c\n";
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(0 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17d\n";
      glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(1 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17e\n";
      glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error setup_atom_pull_restraints_glsl() 17f\n";

      // translate position, 3, size 3 floats
      glEnableVertexAttribArray(3);

      // surely this (annd below) has been set-up already? -- CheckMe.
      glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(3 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17aa\n";

      // positions, 4, size 3 floats
      glEnableVertexAttribArray(4);
      err = glGetError(); if (err) std::cout << "GL error bonds 6\n";
      glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(4 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 7\n";

      //  normals, 5, size 3 floats
      glEnableVertexAttribArray(5);
      err = glGetError(); if (err) std::cout << "GL error bonds 11\n";
      glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(5 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 12\n";

      //  colours, 6, size 4 floats
      glEnableVertexAttribArray(6);
      err = glGetError(); if (err) std::cout << "GL error bonds 16\n";
      glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
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
graphics_info_t::draw_atom_pull_restraints() {

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

            glm::mat4 mvp = get_molecule_mvp();
            glm::mat4 model_rotation = get_model_rotation();
            GLuint mvp_location = shader.mvp_uniform_location;
            GLuint view_rotation_location = shader.view_rotation_uniform_location;
            glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
            glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &model_rotation[0][0]);

            std::map<unsigned int, lights_info_t>::const_iterator it;
            unsigned int light_idx = 0;
            it = lights.find(light_idx);
            if (it != lights.end())
               shader.setup_light(light_idx, it->second, model_rotation, eye_position);
            light_idx = 1;
            it = lights.find(light_idx);
            if (it != lights.end())
               shader.setup_light(light_idx, it->second, model_rotation, eye_position);

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

void
graphics_info_t::draw_molecular_triangles() {
#ifdef USE_MOLECULES_TO_TRIANGLES
   // Martin's triangular molecules
   //
   // centre of the screen
   FCXXCoord pos(graphics_info_t::RotationCentre_x(),
                 graphics_info_t::RotationCentre_y(),
                 graphics_info_t::RotationCentre_z());

   glm::vec3 eye_position = get_world_space_eye_position();
   FCXXCoord eye_pos(eye_position.x, eye_position.y, eye_position.z);

   // std::cout << "eye_pos: " << eye_pos << "\n";
   // coot::Cartesian eye_cart = pos + 20 * diff;
   // FCXXCoord eye_pos(eye_cart.x(), eye_cart.y(), eye_cart.z());
   if (graphics_info_t::mol_tri_scene_setup) {
      if (graphics_info_t::mol_tri_renderer) {
         //Can retrieve reference to the light if so preferred
         FCXXCoord light_pos = pos;
         FCXXCoord neg_light_pos = pos;

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

// static
void
graphics_info_t::draw_particles() {

   if (! particles.empty()) {
      if (mesh_for_particles.have_instances()) {
         glm::mat4 mvp = get_molecule_mvp();
         glm::mat4 model_rotation = get_model_rotation();
         mesh_for_particles.draw_particles(&shader_for_particles, mvp, model_rotation);
      }
   }
}

// static
void
graphics_info_t::draw_happy_face_residue_markers() {

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
         glm::mat4 mvp = get_molecule_mvp();
         glm::mat4 model_rotation = get_model_rotation();
         texture_for_happy_face_residue_marker.Bind(0);
         unsigned int draw_count = draw_count_for_happy_face_residue_markers;
         unsigned int draw_count_max = g.draw_count_max_for_happy_face_residue_markers;
         tmesh_for_happy_face_residues_markers.draw_instances(&shader_for_happy_face_residue_markers,
                                                              mvp, model_rotation,
                                                              draw_count, draw_count_max);
      }
   }
}

// static
void
graphics_info_t::draw_anchored_atom_markers()  {

   if (tmesh_for_anchored_atom_markers.draw_this_mesh) {
      if (tmesh_for_anchored_atom_markers.have_instances()) {
         glm::mat4 mvp = get_molecule_mvp();
         glm::mat4 model_rotation = get_model_rotation();
         texture_for_anchored_atom_markers.Bind(0);
         // tmesh_for_anchored_atom_markers.draw_instances(&shader_for_happy_face_residue_markers, mvp, view_rotation);
         std::cout << "FIXME draw_anchored_atom_markers() merged code " << std::endl;
      }
   }
}

void
graphics_info_t::draw_texture_meshes() {

   // std::cout << "draw_texture_meshes() --- start --- " << std::endl;

   if (! texture_meshes.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      Shader &shader = shader_for_texture_meshes;
      for (unsigned int i=0; i<texture_meshes.size(); i++) {
         TextureMesh &tm = texture_meshes[i];
         // std::cout << "debug:: in draw_textures_meshes textures size " << tm.textures.size() << std::endl;
         if (! tm.textures.empty()) {
            // std::cout << "Binding and drawing the texture mesh" << std::endl;

            std::cout << "............ get crow texture drawing code" << std::endl;
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
         texture_for_hud_refinement_dialog_arrow_highlighted.Bind(0);
      } else {
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
graphics_info_t::draw_molecules() {

   // opaque things

   draw_outlined_active_residue();

   draw_intermediate_atoms(PASS_TYPE_STANDARD);

   draw_intermediate_atoms_rama_balls(PASS_TYPE_STANDARD);

   draw_atom_pull_restraints();

   draw_meshed_generic_display_object_meshes(PASS_TYPE_STANDARD);

   draw_instanced_meshes();

   draw_map_molecules(false); // transparency

   draw_unit_cells();

   draw_environment_graphics_object();

   draw_generic_objects();

   draw_hydrogen_bonds_mesh(); // like boids

   draw_boids();

   draw_particles();

   draw_happy_face_residue_markers();

   draw_anchored_atom_markers();

   // this is the last opaque thing to be drawn because the atom labels are blended.
   // It should be easy to break out the atom label code into its own function. That
   // might be better.
   //
   draw_model_molecules();

   // transparent things...

   draw_map_molecules(true);

}

void
graphics_info_t::draw_molecules_with_shadows() {

   int n_mols = n_molecules();
   bool show_just_shadows = false;
   glm::mat4 mvp = get_molecule_mvp();
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
               m.molecule_as_mesh.draw_simple_bond_lines(&shader_for_symmetry_atoms_bond_lines, mvp, bg_col_v4, lw, shader_do_depth_fog_flag);

            } else {

               float opacity = 1.0;
               m.molecule_as_mesh.draw_with_shadows(&shader_for_meshes_with_shadows, mvp, model_rotation_matrix, lights,
                                                    eye_position, opacity, bg_col_v4, shader_do_depth_fog_flag, light_view_mvp,
                                                    shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
            }

            // this has not been shadowified:
            m.draw_dots(&shader_for_rama_balls, mvp, model_rotation_matrix, lights, eye_position,
                        bg_col_v4, shader_do_depth_fog_flag);

            m.draw_ncs_ghosts(&shader_for_meshes, mvp, model_rotation_matrix, lights, eye_position, bg_col_v4);

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


   // convert these to shadow versions

   draw_model_molecules_with_shadows(); // does symmetry

   draw_outlined_active_residue();

   draw_intermediate_atoms(PASS_TYPE_STANDARD);

   draw_intermediate_atoms_rama_balls(PASS_TYPE_STANDARD);

   draw_atom_pull_restraints();

   draw_meshed_generic_display_object_meshes(PASS_TYPE_STANDARD);

   draw_instanced_meshes();

   draw_map_molecules(false); // transparency

   draw_unit_cells();

   draw_environment_graphics_object();

   draw_generic_objects();

   draw_hydrogen_bonds_mesh(); // like boids

   draw_boids();

   draw_particles();

   draw_happy_face_residue_markers();

   draw_anchored_atom_markers();

   // this is the last opaque thing to be drawn because the atom labels are blended.
   // It should be easy to break out the atom label code into its own function. That
   // might be better.
   //
   // draw_model_molecules();

   // transparent things...

   // draw_map_molecules(true);


}


void
graphics_info_t::draw_molecules_atom_labels() {


   // calls draw_molecule_atom_labels(m, mvp, model_rotation_matrix);

   int n_mols = n_molecules();
   glm::mat4 mvp = get_molecule_mvp();
   auto model_rotation_matrix = get_model_rotation();

   for (int i=0; i<n_mols; i++) {
      if (is_valid_model_molecule(i)) {
         // glEnable(GL_BLEND); // surely this is called inside draw_labels()?
         molecule_class_info_t &m = molecules[i];
         if (m.draw_it) {
            draw_molecule_atom_labels(m, mvp, model_rotation_matrix);
         }
      }
   }


}


// This does (draws) symmetry too.
//
// static
void
graphics_info_t::draw_environment_graphics_object() {

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
            glm::mat4 mvp = get_molecule_mvp();
            glm::vec3 eye_position = get_world_space_eye_position();
            glm::mat4 model_rotation = get_model_rotation();
            glm::vec4 bg_col(background_colour, 1.0);

            bool do_depth_fog = shader_do_depth_fog_flag;
            bool show_just_shadows = false;
            bool wireframe_mode = false;
            float opacity = 1.0f;
            mesh_for_environment_distances.mesh.draw(&shader_for_moleculestotriangles,
                                                     mvp, model_rotation,
                                                     lights, eye_position, opacity, bg_col,
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
                                                   mvp, model_rotation, bg_col, do_depth_fog,
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

void
graphics_info_t::update_mesh_for_outline_of_active_residue(int imol, const coot::atom_spec_t &spec) {

   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      if (mol) {
         coot::residue_spec_t res_spec(spec);
         mmdb::Residue *residue_p = molecules[imol].get_residue(res_spec);
         if (residue_p) {
            // what about Mesh's make_graphical_bonds_bonds_bonds adn make_graphical_bonds_hemispherical_atoms
            molecular_mesh_generator_t mmg;
            std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > p =
               mmg.get_molecular_triangles_mesh_for_active_residue(imol, mol, residue_p, Geom_p());
            mesh_for_outline_of_active_residue.clear();
            mesh_for_outline_of_active_residue.import(p);
            Material mat;
            mesh_for_outline_of_active_residue.setup(mat);
         }
      }
   }
}


void
graphics_info_t::draw_outlined_active_residue() {

   if (outline_for_active_residue_frame_count > 0) {
      glm::mat4 mvp = get_molecule_mvp();
      std::map<unsigned int, lights_info_t> dummy_lights;
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      Shader &shader = shader_for_outline_of_active_residue;
      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      mesh_for_outline_of_active_residue.draw(&shader, mvp, model_rotation, dummy_lights, eye_position, opacity,
                                              bg_col, wireframe_mode, false, show_just_shadows);
   }
};


void
graphics_info_t::draw_unit_cells() {

   glm::mat4 mvp = get_molecule_mvp();
   for (int ii=n_molecules()-1; ii>=0; ii--) {
      molecule_class_info_t &m = molecules[ii];
      m.draw_unit_cell(&shader_for_lines, mvp);
   }

}

void
graphics_info_t::draw_meshed_generic_display_object_meshes(unsigned int pass_type) {

   // This function doesn't draw these
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots H-bond
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots wide-contact
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots close-contact
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots small-overlap
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots vdw-surface
   // graphics_info_t::draw_instanced_meshes() A Molecule 2: Ligand Contact Dots big-overlap

   // std::cout << "------------- draw_meshed_generic_display_object_meshes() " << std::endl;

   bool draw_meshes = true;
   bool draw_mesh_normals = false;

   glm::vec3 eye_position = get_world_space_eye_position();
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 mvp_orthogonal = glm::mat4(1.0f); // placeholder
   glm::mat4 model_rotation = get_model_rotation();
   glm::vec4 bg_col(background_colour, 1.0);
   bool do_depth_fog = true;

   unsigned int light_index = 0;
   std::map<unsigned int, lights_info_t>::const_iterator it = lights.find(light_index);
   if (it != lights.end()) {
      graphics_info_t g;
      const auto &light = it->second;
      mvp_orthogonal = g.get_mvp_for_shadow_map(light.direction);
   }

   //std::cout << "mvp diag "
   // << mvp[0][0] << " " << mvp[1][1] << " " << mvp[2][2] << std::endl;

   glm::mat3 vrm(glm::toMat4(graphics_info_t::glm_quat));

   // Yes, identity matrix
   // std::cout << "p: " << glm::to_string(p) << std::endl;

   if (draw_meshes) { //local, debugging
      bool have_meshes_to_draw = false;
      for (int i=n_molecules()-1; i>=0; i--) {
         if (! molecules[i].meshes.empty()) {
            have_meshes_to_draw = true;
            break;
         }
      }

      if (have_meshes_to_draw) {
         // std::cout << "   Here A in draw_meshed_generic_display_object_meshes() " << std::endl;
         glDisable(GL_BLEND);
         for (int ii=n_molecules()-1; ii>=0; ii--) {
            // std::cout << "Here B in draw_meshed_generic_display_object_meshes() " << ii  << std::endl;
            molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
            for (unsigned int jj=0; jj<m.meshes.size(); jj++) {
               Mesh &mesh = m.meshes[jj];

               if (false)
                  std::cout << "mesh jj " << jj << " of " << m.meshes.size()
                            << " instanced: " << m.meshes[jj].is_instanced << std::endl;

               if (mesh.is_instanced) {

                  mesh.draw_instanced(&shader_for_moleculestotriangles, mvp,
                                              model_rotation, lights, eye_position,
                                              bg_col, do_depth_fog);
               } else {
                  if (pass_type == PASS_TYPE_STANDARD) {
                     bool show_just_shadows = false;
                     bool wireframe_mode = false;
                     float opacity = 1.0f;
                     m.meshes[jj].draw(&shader_for_meshes_with_shadows, mvp,
                                       model_rotation, lights, eye_position, opacity, bg_col,
                                       wireframe_mode, do_depth_fog, show_just_shadows);
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
                  if (pass_type == PASS_TYPE_FOR_SHADOWS) { // i.e. generating the shadow map, not using it.

                     glm::vec3 dummy_eye_position;
                     bool gl_lines_mode = false;
                     bool show_just_shadows = false;
                     bool opacity = 1.0;

                     mesh.draw(&shader_for_meshes_shadow_map,
                               mvp_orthogonal,
                               model_rotation,
                               lights,
                               dummy_eye_position,
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
graphics_info_t::draw_instanced_meshes() {

   // presumes opaque-only

   bool have_meshes_to_draw = false;
   for (int i=n_molecules()-1; i>=0; i--) {
      if (! molecules[i].instanced_meshes.empty()) {
         if (molecules[i].draw_it) {
            have_meshes_to_draw = true;
            break;
         }
      }
   }
   if (! instanced_meshes.empty())
      have_meshes_to_draw = true;

   if (have_meshes_to_draw) {
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool do_depth_fog = shader_do_depth_fog_flag;
      glDisable(GL_BLEND);
      for (int ii=n_molecules()-1; ii>=0; ii--) {
         molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
         if (molecules[ii].draw_it) {
            for (unsigned int jj=0; jj<m.instanced_meshes.size(); jj++) {
               // std::cout << "   graphics_info_t::draw_instanced_meshes() A " << m.instanced_meshes[jj].get_name() << std::endl;
               m.instanced_meshes[jj].draw(&shader_for_rama_balls, mvp,
                                           model_rotation, lights, eye_position, bg_col, do_depth_fog);
            }
         }
      }

      // And draw our own

      if (! instanced_meshes.empty()) {
         for (unsigned int jj=0; jj<instanced_meshes.size(); jj++) {
            // std::cout << "   graphics_info_t::draw_instanced_meshes() our own " << jj << " "
            // << instanced_meshes[jj].get_name() << std::endl;
            instanced_meshes[jj].draw(&shader_for_rama_balls, mvp,
                                      model_rotation, lights, eye_position, bg_col, do_depth_fog);
         }
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
graphics_info_t::draw_cube(GtkGLArea *glarea, unsigned int cube_type) {

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

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 model_rotation = get_model_rotation(); // hhmm... naming ... 20220212-PE  fixed now.

   glBindVertexArray(central_cube_vertexarray_id);
   err = glGetError(); if (err) std::cout << "   error::draw_central_cube() B err " << err << std::endl;
   glUseProgram(shader_for_central_cube.get_program_id());
   err = glGetError(); if (err) std::cout << "   error::draw_central_cube() C err " << err << std::endl;
   glm::vec3 rc = get_rotation_centre();
   if (cube_type == VIEW_CENTRAL_CUBE) {
      mvp = glm::translate(mvp, rc);
      float s = rotation_centre_cube_size;
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
graphics_info_t::draw_central_cube(GtkGLArea *glarea) {
   draw_cube(glarea, VIEW_CENTRAL_CUBE);
}

void
graphics_info_t::draw_origin_cube(GtkGLArea *glarea) {

   draw_cube(glarea, ORIGIN_CUBE);
}

void
graphics_info_t::draw_rotation_centre_crosshairs(GtkGLArea *glarea, unsigned int pass_type) {

   // gtk_gl_area_make_current(glarea); // needed?, no it isn't.
   GLenum err = glGetError();
   if (err) std::cout << "error draw_rotation_centre_crosshairs() A0 err " << err << std::endl;

   glLineWidth(1.0);
   err = glGetError();
   if (err) std::cout << "error draw_rotation_centre_crosshairs() A1 err " << err << std::endl;

   glm::mat4 mvp = get_molecule_mvp();
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
   float s = 6.0f * rotation_centre_cube_size;
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
      glm::vec4 line_colour(0.8f, 0.8f, 0.8f, 1.0f);
      if (! is_bb) 
         line_colour = glm::vec4(0.2f, 0.2f, 0.2f, 1.0f);

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



GtkWidget *create_and_pack_gtkglarea(GtkWidget *vbox, bool use_gtk_builder) {

   // the use_gtk_builder flag really means "was invoked from the path that..."

   GtkWidget *w = gtk_gl_area_new();

   auto get_gl_widget_dimension_scale_factor  = [] () {
                                                   int sf = 1;
                                                   char *e = getenv("COOT_OPENGL_WIDGET_SCALE_FACTOR");
                                                   if (e) {
                                                      std::string ee(e);
                                                      sf = std::stoi(ee);
                                                   }
                                                   return sf;
                                                };

   // allow the user to set the major and minor version (for debugging)

   int opengl_major_version = 3;
   int opengl_minor_version = 3;
   char *e1 = getenv("COOT_OPENGL_MAJOR_VERSION");
   char *e2 = getenv("COOT_OPENGL_MINOR_VERSION");
   if (e1) {
      std::string e1s(e1);
      opengl_major_version = std::stoi(e1s);
   }
   if (e2) {
      std::string e2s(e2);
      opengl_minor_version = std::stoi(e2s);
   }
   std::cout << "DEBUG:: setting OpenGL required version to "
             << opengl_major_version << " " << opengl_minor_version << std::endl;

   gtk_gl_area_set_required_version(GTK_GL_AREA(w), opengl_major_version, opengl_minor_version);

   unsigned int dimensions = 700;
   if (! use_gtk_builder) dimensions = 900;
   dimensions = 900;
   int gl_widget_dimension_scale_factor = get_gl_widget_dimension_scale_factor();
   gtk_widget_set_size_request(w,
                               gl_widget_dimension_scale_factor * dimensions,
                               gl_widget_dimension_scale_factor * dimensions);
   gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 0);
   return w;
}

void
graphics_info_t::draw_measure_distance_and_angles() {

   if (mesh_for_measure_distance_object_vec.get_draw_this_mesh()) {
      Shader &shader = shader_for_moleculestotriangles;
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      mesh_for_measure_distance_object_vec.draw(&shader, mvp, model_rotation_matrix, lights, eye_position,
                                                opacity, bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

      mesh_for_measure_angle_object_vec.draw(&shader, mvp, model_rotation_matrix, lights, eye_position,
                                             opacity, bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

      if (! labels_for_measure_distances_and_angles.empty()) {
         Shader &shader = shader_for_atom_labels;
         for (unsigned int i=0; i<labels_for_measure_distances_and_angles.size(); i++) {
            const auto &label = labels_for_measure_distances_and_angles[i];
            tmesh_for_labels.draw_atom_label(label.label, label.position, label.colour, &shader,
                                             mvp, model_rotation_matrix, bg_col,
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
   graphics_info_t::lights[0] = light;

   light.position = glm::vec4(3.0f, -2.0f, 4.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(-1.0, 0.5, 1.0));
   // light.is_on = false;
   graphics_info_t::lights[1] = light;
}

// this is called from realize()
void
graphics_info_t::setup_hud_geometry_bars() {

   if (! glareas[0]) return;

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;
   float aspect_ratio = static_cast<float>(w)/static_cast<float>(h);

   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   shader_for_hud_geometry_bars.Use();

   mesh_for_hud_geometry.setup_camera_facing_quad_for_bar();
   mesh_for_hud_geometry.setup_instancing_buffer(500, sizeof(HUD_bar_attribs_t));

   // If not found in this directory, then try default directory.
   texture_for_hud_geometry_labels_map["Rama"].set_default_directory(coot::package_data_dir());
   texture_for_hud_geometry_labels_map["Rama"].init("hud-label-rama-small.png");
   // texture_for_hud_geometry_labels_map["Rama"].init("rama-plot-other-normal.png");
   texture_for_hud_geometry_labels_map["NBC"].set_default_directory(coot::package_data_dir());
   texture_for_hud_geometry_labels_map["NBC"].init("hud-label-NBC-small.png");
   //texture_for_hud_geometry_labels_map["NBC"].init("rama-plot-other-normal.png");
   texture_for_hud_geometry_labels_map["Rota"].set_default_directory(coot::package_data_dir());
   texture_for_hud_geometry_labels_map["Rota"].init("hud-label-rota-small.png");
   //texture_for_hud_geometry_labels_map["Rota"].init("rama-plot-other-normal.png");
   texture_for_hud_geometry_labels_map["Pull"].set_default_directory(coot::package_data_dir());
   texture_for_hud_geometry_labels_map["Pull"].init("hud-label-pull-small.png");
   // texture_for_hud_geometry_labels_map["Pull"].init("rama-plot-other-normal.png");

   texture_for_hud_tooltip_background.set_default_directory(coot::package_data_dir());
   texture_for_hud_tooltip_background.init("hud-tooltip.png"); // 94x47
   float sc_x = 0.1 * static_cast<float>(103) / aspect_ratio;
   float sc_y = 0.01 * static_cast<float>(50);

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

}

void
graphics_info_t::setup_hud_buttons() {

   if (! glareas[0]) return;

   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   shader_for_hud_buttons.Use();
   mesh_for_hud_buttons.setup_vertices_and_triangles_for_button(); // instanced button
   unsigned int n_buttons_max = 20; // surely 6 is enough?
   mesh_for_hud_buttons.setup_instancing_buffer(n_buttons_max, sizeof(HUD_button_info_t));
   // maybe mesh_for_hud_buttons.close() ?
}

void
graphics_info_t::clear_hud_buttons() {

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
graphics_info_t::clear_gl_rama_plot() {

   gl_rama_plot.clear();
}

void
graphics_info_t::draw_hud_ramachandran_plot() {

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   // auto tp_0 = std::chrono::high_resolution_clock::now();
   bool draw_gl_ramachandran_plot = true;
   // 20220403-PE Do I want to draw this only if there are Rama restraints? No.
   // This test was added in 5cba65a245693929a276961928124b112e10291b
   // if (! saved_dragged_refinement_results.refinement_results_contain_overall_rama_plot_score)
   //       draw_gl_ramachandran_plot = false;p
   if (draw_gl_ramachandran_plot) { // make this a member of graphics_info_t
      if (moving_atoms_asc) {
         if (moving_atoms_asc->n_selected_atoms > 0) {
            gl_rama_plot.setup_from(imol_moving_atoms, moving_atoms_asc->mol); // checks to see if an update is acutally needed.
            gl_rama_plot.draw(&shader_for_rama_plot_axes_and_ticks,
                              &shader_for_rama_plot_phi_phis_markers, // instanced
                              &shader_for_hud_image_texture, w, h, w, h); // background texture (not text!), uses window_resize_position_correction
         }
      }
   }

   // auto tp_1 = std::chrono::high_resolution_clock::now();
   // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: draw_ramachandran_plot() " << d10 << " microseconds" << std::endl;
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

   if (use_graphics_interface_flag) {
      GtkWidget *button_1 = get_widget_from_builder("clear_atom_pull_restraints_toolbutton");
      GtkWidget *button_2 = get_widget_from_builder("auto_clear_atom_pull_restraints_togglebutton");

      gtk_widget_show(button_1);
      gtk_widget_show(button_2);
   }
}


void
graphics_info_t::hide_atom_pull_toolbar_buttons() {

   if (use_graphics_interface_flag) {
      GtkWidget *button_1 = get_widget_from_builder("clear_atom_pull_restraints_toolbutton");
      GtkWidget *button_2 = get_widget_from_builder("auto_clear_atom_pull_restraints_togglebutton");
      
      gtk_widget_hide(button_1);
      gtk_widget_hide(button_2);
   }
}

void
graphics_info_t::show_accept_reject_hud_buttons() {

   // add some HUD buttons

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   HUD_button_info_t button_1("OK");
   HUD_button_info_t button_2("Cancel");
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
                           g.hud_button_info.clear();
                           g.hide_atom_pull_toolbar_buttons();
                           g.clear_gl_rama_plot();
                           g.graphics_draw();
                           return true;
                   };

   auto button_2_func = [] () {
                           std::cout << "------------------- HUD button Cancel callback start " << std::endl;
                           graphics_info_t g;
                           g.stop_refinement_internal();
                           g.clear_up_moving_atoms();
                           g.hud_button_info.clear();
                           g.rebond_molecule_corresponding_to_moving_atoms();
                           g.graphics_draw();
                           g.hide_atom_pull_toolbar_buttons();
                           g.clear_gl_rama_plot();
                           std::cout << "------------------- HUD button Cancel callback done " << std::endl;
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
   hud_label_info.push_back(hud_label_info_t("Pull", 0, 0.7));
   hud_label_info.push_back(hud_label_info_t("Rama", 2));
   hud_label_info.push_back(hud_label_info_t("Rota", 3, 0.9));
   hud_label_info.push_back(hud_label_info_t("NBC",  1, 0.9));
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
                                               float rotation_amount = 1.0 - 0.012 * distortion;
                                               if (rotation_amount < 0.68) rotation_amount = 0.68;
                                               return rotation_amount;
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

   if (rr.refinement_results_contain_overall_nbc_score)
      add_bars(rr.sorted_nbc_baddies, 1, &new_bars, moving_atoms_active_residue,
               x_base_for_hud_geometry_bars, distortion_to_rotation_amount_nbc,
               hud_geometry_distortion_to_bar_size_nbc);

   if (rr.refinement_results_contain_overall_rama_plot_score) {
      // std::cout << "add_bars() for rama with " << rr.sorted_rama_baddies.size() << " sorted baddies" << std::endl;
      add_bars(rr.sorted_rama_baddies, 2, &new_bars, moving_atoms_active_residue, x_base_for_hud_geometry_bars,
               hud_geometry_distortion_to_rotation_amount_rama, hud_geometry_distortion_to_bar_size_rama);
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

   if (! new_bars.empty()) {
      // std::cout << "new bar size " << new_bars.size() << std::endl;
      mesh_for_hud_geometry.update_instancing_buffer_data(new_bars);
      Shader &shader = shader_for_hud_geometry_bars;
      mesh_for_hud_geometry.draw(&shader);
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

   if (!status_pair.first)
      if (rr.refinement_results_contain_overall_nbc_score)
         status_pair = check_blocks(rr.sorted_nbc_baddies, 1, x_base_for_hud_geometry_bars,
                                    hud_geometry_distortion_to_bar_size_nbc, act_on_hit);

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
                                              mesh_for_hud_buttons.update_instancing_buffer_data(hud_button_info);
                                              graphics_draw(); // let's see the changes then
                                           };
   auto unhighlight_all_buttons = [] () {
                                              for (unsigned int i=0; i<hud_button_info.size(); i++) {
                                                 auto &button = hud_button_info[i];
                                                 button.set_button_colour_for_mode(HUD_button_info_t::BASIC);
                                              }
                                              mesh_for_hud_buttons.update_instancing_buffer_data(hud_button_info);
                                  };

   bool status = false;
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
               std::cout << "Act on button " << i << " callback" << std::endl;
               if (button.callback_function) {
                  button.callback_function();
               }
            }
            status = true;
            highlight_just_button_with_index(i);
         }
      }
      if (!status) {
         unhighlight_all_buttons();
      }
   }

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

}

void
graphics_info_t::render_3d_scene(GtkGLArea *gl_area) {

   // note: this function is called from render_scene_sans_depth_blur()

   //  ------------------- render scene ----------------------------

   glEnable(GL_DEPTH_TEST);

   // const glm::vec3 &bg = graphics_info_t::background_colour;
   // glClearColor (bg[0], bg[1], bg[2], 1.0); // what difference does this make?

   GLenum err = glGetError(); if (err) std::cout << "render_3d_scene lambda B err " << err << std::endl;
   // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   err = glGetError(); if (err) std::cout << "render_3d_scene lambda C err " << err << std::endl;

   draw_origin_cube(gl_area);
   err = glGetError(); if (err) std::cout << "render scene lambda post cubes err " << err << std::endl;

   draw_rotation_centre_crosshairs(gl_area, PASS_TYPE_STANDARD);

   draw_molecules(); // includes particles, happy-faces and boids (should they be there (maybe not))
                     // so rename this function? Or just bring everything here?  Put this render() function
                     // into new file graphics-info-opengl-render.cc

   draw_invalid_residue_pulse();

   draw_identification_pulse();

   draw_delete_item_pulse();

   draw_measure_distance_and_angles(); // maybe in draw_molecules()?

   draw_pointer_distances_objects();

   draw_texture_meshes();

}

void
graphics_info_t::render_3d_scene_with_shadows() {

   // note: this function is called from render_scene_sans_depth_blur()

   //  ------------------- render scene ----------------------------

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);

   // Try not using the meshes-with-shadows.shader
   // render_3d_scene(gl_area);
   // return;

   glEnable(GL_DEPTH_TEST);
   GLenum err = glGetError();
   if (err) std::cout << "render_3d_scene_with_shadows B err " << err << std::endl;

   draw_origin_cube(gl_area);
   err = glGetError(); if (err) std::cout << "render scene lambda post cubes err " << err << std::endl;

   draw_rotation_centre_crosshairs(gl_area, PASS_TYPE_STANDARD);

   draw_molecules_with_shadows(); // includes particles, happy-faces and boids (should they be there (maybe not))
                                  // so rename this function? Or just bring everything here?  Put this render() function
                                  // into new file graphics-info-opengl-render.cc

   draw_invalid_residue_pulse();

   draw_identification_pulse();

   draw_delete_item_pulse();

   draw_measure_distance_and_angles(); // maybe in draw_molecules()?

   draw_pointer_distances_objects();

   draw_texture_meshes();

}


void
graphics_info_t::render_3d_scene_for_ssao() {

}



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

   unsigned int frame_time_history_list_max_n_elements = 500;
   GtkWidget *glarea = glareas[0];
   if (glarea) {
      auto tp_now = std::chrono::high_resolution_clock::now();
      frame_time_history_list.push_back(tp_now);
      if (frame_time_history_list.size() >= (frame_time_history_list_max_n_elements+1))
         frame_time_history_list.pop_front();
   }

   // 20220405-PE the control flow is confusing/unconentional here - I should make it more clear

   bool use_crow_code = true;
   if (to_screendump_framebuffer_flag) use_crow_code = false;

   if (use_crow_code) {
      gboolean state = render_scene();

      draw_hud_elements();

      glFlush();
      update_fps_statistics();
      return state;
   }

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   bool make_image_for_screen = ! to_screendump_framebuffer_flag;

#ifdef __APPLE__
   use_framebuffers = false;
#endif

   if (use_framebuffers) { // static class variable

      glViewport(0, 0, framebuffer_scale * w, framebuffer_scale * h);
      GLenum err = glGetError();
      if (err) std::cout << "GL ERROR:: render() post glViewport() err " << err << std::endl;
      screen_framebuffer.bind(); // screen_ao, that is
      err = glGetError();
      if (err) std::cout << "GL ERROR:: render() post screen_framebuffer bind() err " << err << std::endl;

      render_3d_scene(gl_area);

      if (make_image_for_screen) { // a normal draw

         glViewport(0, 0, w, h);

         // use this, rather than glBindFramebuffer(GL_FRAMEBUFFER, 0); ... just Gtk things.
         // gtk_gl_area_attach_buffers(gl_area);

         if (shader_do_depth_of_field_blur_flag) {

            blur_y_framebuffer.bind();
            render_scene_with_screen_ao_shader();
            blur_x_framebuffer.bind();
            render_scene_with_y_blur();

            combine_textures_using_depth_framebuffer.bind();

            render_scene_with_x_blur();

            gtk_gl_area_attach_buffers(gl_area);

            render_scene_with_texture_combination_for_depth_blur();

            // And finally draw the HUD elements to the GTK framebuffer
            draw_hud_elements();

         } else {

            gtk_gl_area_attach_buffers(gl_area);
            render_scene_with_screen_ao_shader();
            // And finally draw the HUD elements to the GTK framebuffer
            draw_hud_elements();

         }

      } else {

         // screendump
         glDisable(GL_DEPTH_TEST);
         unsigned int sf = framebuffer_scale;
         glViewport(0, 0, sf * w, sf * h);
         framebuffer screendump_framebuffer;
         unsigned int index_offset = 0;
         screendump_framebuffer.init(sf * w, sf * h, index_offset, "screendump");
         screendump_framebuffer.bind();
         render_3d_scene(gl_area);
         // render_scene_with_screen_ao_shader();
         gtk_gl_area_attach_buffers(gl_area);
         screendump_tga_internal(output_file_name, w, h, sf, screendump_framebuffer.get_fbo());

      }
   } else {

      // simple/direct - for debugging framebuffers
      gtk_gl_area_attach_buffers(gl_area);
      render_3d_scene(gl_area);
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
}

void
graphics_info_t::render_scene_with_screen_ao_shader() {

   glEnable(GL_DEPTH_TEST);
   shader_for_screen.Use();
   glBindVertexArray(screen_quad_vertex_array_id);

   const glm::vec3 &bg = background_colour;
   glClearColor(bg[0], bg[1], bg[2], 1.0); // this can be seen
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glActiveTexture(GL_TEXTURE0 + 0);
   glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_depth());
   shader_for_screen.set_int_for_uniform("screenTexture", 0);
   shader_for_screen.set_int_for_uniform("screenDepth", 1);
   GLenum err = glGetError(); if (err) std::cout << "render() D err " << err << std::endl;
   shader_for_screen.set_bool_for_uniform("do_ambient_occlusion", shader_do_ambient_occlusion_flag);
   shader_for_screen.set_bool_for_uniform("do_outline", shader_do_outline_flag);
   glm::vec4 background_col(background_colour, 1.0);
   bool background_is_black = background_is_black_p();
   shader_for_screen.set_bool_for_uniform("background_is_dark", background_is_black);

   glDrawArrays(GL_TRIANGLES, 0, 6);
   err = glGetError(); if (err) std::cout << "render() E err " << err << std::endl;

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
graphics_info_t::reset_frame_buffers(int width, int height) {

   if (use_framebuffers) {


      // 20220108-PE note to self. Try using the framebuffer::reset() function instead

      unsigned int sf = framebuffer_scale;
      unsigned int index_offset = 0;

      // width  = width;
      // height = height;
      if (false)
         std::cout << "debug:: reset_frame_buffers() with sf " << sf << " "
                   << width << " x " << height << std::endl;
      screen_framebuffer.init(sf * width, sf * height, index_offset, "screen");
      GLenum err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      blur_x_framebuffer.init(sf * width, sf * height, index_offset, "blur-x");
      err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      blur_y_framebuffer.init(sf * width, sf * height, index_offset, "blur-y");
      err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      combine_textures_using_depth_framebuffer.init(sf * width, sf * height, index_offset, "combine");
      err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

      // std::cout << "debug:: reset_frame_buffers() sf " << sf << " width " << width << " height " << height << std::endl;

      // index_offset = 0;
      // g.blur_framebuffer.init(width, height, index_offset, "blur");

      // ------------------ crows code --------------------------

      // note to self:
      // the shadow texture doesn't need to change - it's under user control, not
      // dependent on the window size

      framebuffer_for_ssao_gbuffer.reset_test(width, height);

      gint w = width;
      gint h = height;

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
   glm::vec3 step = 0.005 * step_size * zoom * screen_x_uv;
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

      shader_for_particles.Use();

      err = glGetError();
      if (err) std::cout << "Error::- setup_draw_for_particles() Post Use() err is "
                         << err << std::endl;

      std::vector<glm::vec3> positions = get_particle_centre_positions();
      particles.make_particles(n_particles, positions);
      std::cout << "setup_draw_for_particles(): done making " << n_particles << " particles "
                << " for " << positions.size() << " positions" << std::endl;

      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
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

   const unsigned int max_happy_faces = 200; // surely enough?
   
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   GLenum err = glGetError();
   if (err) std::cout << "Error::- setup_draw_for_happy_face_residue_markers_init() "
                      << "Post attach buffers err is " << err << std::endl;

   // If not found in this directory, then try default directory.
   texture_for_happy_face_residue_marker.set_default_directory(coot::package_data_dir());
   texture_for_happy_face_residue_marker.init("happy-face-marker.png");

   shader_for_happy_face_residue_markers.Use();
   tmesh_for_happy_face_residues_markers.setup_camera_facing_quad(&shader_for_happy_face_residue_markers, 1.0, 1.0);
   tmesh_for_happy_face_residues_markers.setup_instancing_buffers(max_happy_faces);
   tmesh_for_happy_face_residues_markers.draw_this_mesh = false;

}

void
graphics_info_t::setup_draw_for_anchored_atom_markers_init() {

   // run this once - call from realize()

   const unsigned int max_anchored_atoms = 200; // surely enough?
   
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   GLenum err = glGetError();
   if (err) std::cout << "Error::- setup_draw_for_anchored_atom_markers_init() "
                      << "Post attach buffers err is " << err << std::endl;

   // If not found in this directory, then try default directory.
   texture_for_anchored_atom_markers.set_default_directory(coot::package_data_dir());
   texture_for_anchored_atom_markers.init("anchor-for-fixed-atoms.png");

   shader_for_anchored_atom_markers.Use();
   tmesh_for_anchored_atom_markers.setup_camera_facing_quad(&shader_for_happy_face_residue_markers, 1.0, 1.0);
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
         gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
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
         v.push_back(sc + 5.0 * clipper_to_glm(p));
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
      std::pair<glm::vec3, glm::vec3> start_end(glm::vec3(0.95,0,0), glm::vec3(-0.95,0,0));
      m.add_cone(start_end, col, 1.0, 0.0, 12, false, true,
                 meshed_generic_display_object::FLAT_CAP,
                 meshed_generic_display_object::FLAT_CAP);
      mesh_for_boids = m.mesh;

      std::vector<glm::mat4>    mats(n_boids);
      std::vector<glm::vec4> colours(n_boids);
      for (unsigned int i=0; i<n_boids; i++) {
         const fun::boid boid = boids[i];
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

   GtkAllocation allocation;
   gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
   float w = allocation.width;
   float h = allocation.height;
   GLenum err = glGetError();
   if (err)
      std::cout << "draw_ligand_view() --- start --- " << err << std::endl;

   myglLineWidth(2.0);
   graphics_ligand_mesh_molecule.draw(&shader_for_ligand_view,
                                      &shader_for_hud_geometry_tooltip_text,
                                      w, h, ft_characters);
   err = glGetError();
   if (err)
      std::cout << "draw_ligand_view() --- end --- " << err << std::endl;
}



void
graphics_info_t::draw_boids() {

   if (boids.size() > 0) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool show_just_shadows = false;
      bool wireframe_mode = false;
      float opacity = 1.0f;
      mesh_for_boids.draw(&shader_for_instanced_objects,
                          mvp, model_rotation_matrix, lights, eye_position, opacity, bg_col,
                          wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

      lines_mesh_for_boids_box.draw(&shader_for_lines, mvp, model_rotation_matrix);
   }
}

void
graphics_info_t::draw_hydrogen_bonds_mesh() {

   // 20210827-PE  each molecule should have its own hydrogen bond mesh. Not just one of them.
   // Fix that later.

   if (mesh_for_hydrogen_bonds.get_draw_this_mesh()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      glm::vec4 bg_col(background_colour, 1.0);

      // 20211210-PE  note that we are not calling the draw_instanced() - that seems perverse to me.
      // Hmm.
      mesh_for_hydrogen_bonds.draw_instanced(&shader_for_instanced_objects,
                                             mvp, model_rotation_matrix, lights, eye_position, bg_col,
                                             shader_do_depth_fog_flag, false, true, 0, 0, 0, 0.2);
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
                                       lines_mesh_for_delete_item_pulse.clear();
                                       delete_item_pulse_centres.clear();
                                    } else {
                                       float ns = pulse_data->n_pulse_steps;
                                       lines_mesh_for_delete_item_pulse.update_buffers_for_pulse(ns, -1);
                                    }
                                    graphics_draw();
                                    return gboolean(continue_status);
                                 };

   pulse_data_t *pulse_data = new pulse_data_t(0, 20); // 20 matches the number in update_buffers_for_pulse()
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> positions = residue_to_positions(residue_p);
   delete_item_pulse_centres = positions;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = true;
   lines_mesh_for_delete_item_pulse.setup_pulse(broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], delete_item_pulse_func, user_data, NULL);

};

void
graphics_info_t::setup_delete_residues_pulse(const std::vector<mmdb::Residue *> &residues) {

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
                                       lines_mesh_for_delete_item_pulse.clear();
                                       delete_item_pulse_centres.clear();
                                    } else {
                                       float ns = pulse_data->n_pulse_steps;
                                       lines_mesh_for_delete_item_pulse.update_buffers_for_pulse(ns, -1);
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
   delete_item_pulse_centres = all_positions;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = true;
   lines_mesh_for_delete_item_pulse.setup_pulse(broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], delete_item_pulse_func, user_data, NULL);

};


// I think that this function should be part of glarea_tick_func().
// (bring glarea_tick_func() into graphics_info_t.)
//
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
      delete_item_pulse_centres.clear(); // we sneakily use this vector (but no longer)
   } else {
      float ns = pulse_data->n_pulse_steps;
      lines_mesh_for_identification_pulse.update_buffers_for_invalid_residue_pulse(ns);
   }
   graphics_draw();
   return gboolean(continue_status);
}

void
graphics_info_t::setup_invalid_residue_pulse(mmdb::Residue *residue_p) {

   pulse_data_t *pulse_data = new pulse_data_t(0, 24);
   gpointer user_data = reinterpret_cast<void *>(pulse_data);
   std::vector<glm::vec3> residue_positions = residue_to_positions(residue_p);
   delete_item_pulse_centres = residue_positions; // sneakily use a wrongly named function
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
   bool broken_line_mode = false;
   lines_mesh_for_identification_pulse.setup_pulse(broken_line_mode);
   gtk_widget_add_tick_callback(glareas[0], invalid_residue_pulse_function, user_data, NULL);

}


void
graphics_info_t::draw_invalid_residue_pulse() {

   if (! lines_mesh_for_identification_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      myglLineWidth(3.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_invalid_residue_pulse() glLineWidth " << err << std::endl;
      for (auto pulse_centre : delete_item_pulse_centres)
         lines_mesh_for_identification_pulse.draw(&shader_for_lines_pulse,
                                                  pulse_centre, mvp,
                                                  model_rotation_matrix, true);
   }
}


void
graphics_info_t::draw_identification_pulse() {

   if (! lines_mesh_for_identification_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      myglLineWidth(2.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_identification_pulse() glLineWidth " << err << std::endl;
      lines_mesh_for_identification_pulse.draw(&shader_for_lines_pulse,
                                               identification_pulse_centre,
                                               mvp, model_rotation_matrix, true);
   }
}

void
graphics_info_t::draw_delete_item_pulse() {

   if (! lines_mesh_for_delete_item_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 model_rotation_matrix = get_model_rotation();
      myglLineWidth(2.0);
      GLenum err = glGetError();
      if (err) std::cout << "draw_delete_item_pulse() glLineWidth " << err << std::endl;
      for (unsigned int i=0; i<delete_item_pulse_centres.size(); i++) {
         lines_mesh_for_delete_item_pulse.draw(&shader_for_lines_pulse,
                                               delete_item_pulse_centres[i],
                                               mvp, model_rotation_matrix, true);
      }
   }
}

void
graphics_info_t::draw_pointer_distances_objects() {

   if (show_pointer_distances_flag) {
      if (! pointer_distances_object_vec.empty()) {
         Shader &shader = shader_for_moleculestotriangles;
         glm::mat4 mvp = get_molecule_mvp();
         glm::mat4 model_rotation_matrix = get_model_rotation();
         glm::vec4 bg_col(background_colour, 1.0);
         bool show_just_shadows = false;
         bool wireframe_mode = false;
         float opacity = 1.0f;
         mesh_for_pointer_distances.mesh.draw(&shader, mvp, model_rotation_matrix, lights, eye_position, opacity,
                                              bg_col, wireframe_mode, shader_do_depth_fog_flag, show_just_shadows);

         if (! labels_for_pointer_distances.empty()) {
            Shader &shader = shader_for_atom_labels;
            for (unsigned int i=0; i<labels_for_pointer_distances.size(); i++) {
               const auto &label = labels_for_pointer_distances[i];
               tmesh_for_labels.draw_atom_label(label.label, label.position, label.colour, &shader,
                                                mvp, model_rotation_matrix, bg_col,
                                                shader_do_depth_fog_flag, perspective_projection_flag);
            }
         }
      }
   }

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

void
graphics_info_t::setup_key_bindings() {

   graphics_info_t g;

   // if we are serious about user-defined key-bindings all of these functions should be thunks in the user API
   // (and returning gboolean).

   auto l1 = []() { graphics_info_t g; g.adjust_clipping(-0.1); return gboolean(TRUE); };
   auto l2 = []() { graphics_info_t g; g.adjust_clipping( 0.1); return gboolean(TRUE); };
   auto l5 = []() { graphics_info_t g; g.blob_under_pointer_to_screen_centre(); return gboolean(TRUE); };

   auto l6 = []() {

                if (do_tick_spin) {
                   do_tick_spin = false;
                } else {
                   if (! tick_function_is_active()) {
                      int spin_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
                      // this is not a good name if we are storing a generic tick function id.
                      idle_function_spin_rock_token = spin_tick_id;
                   }
                   do_tick_spin = true;

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

   auto l9 = [] () {
                update_go_to_atom_from_current_position();
                return gboolean(TRUE);
             };

   auto l10 = []() { graphics_info_t::zoom *= 0.9; return gboolean(TRUE); };

   auto l11 = []() { graphics_info_t::zoom *= 1.1; return gboolean(TRUE); };

   auto l12 = []() { graphics_info_t g; g.move_forwards(); return gboolean(TRUE); };

   auto l13 = []() { graphics_info_t g; g.move_backwards(); return gboolean(TRUE); };

   auto l13l = []() { graphics_info_t g; g.step_screen_left();   return gboolean(TRUE); };
   auto l13r = []() { graphics_info_t g; g.step_screen_right(); return gboolean(TRUE); };

   auto l14 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('forward')"); return gboolean(TRUE); };

   auto l15 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('backward')"); return gboolean(TRUE); };

   auto l16 = []() { graphics_info_t g; g.undo_last_move(); return gboolean(TRUE); };

   auto l18 = []() { graphics_info_t g; g.clear_hud_buttons(); g.accept_moving_atoms(); return gboolean(TRUE); };

   auto l18_space = []() {
                       graphics_info_t g;
                       if (g.hud_button_info.size()) {
                          g.clear_hud_buttons(); g.accept_moving_atoms();
                       } else {
                          
                          // Move the view - don't click the button
                          
                          // g.reorienting_next_residue_mode = false; // hack
                          bool reorienting = graphics_info_t::reorienting_next_residue_mode;
                          if (reorienting) {
                             if (graphics_info_t::shift_is_pressed) {
                                g.reorienting_next_residue(false); // backwards
                             } else {
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
                       return gboolean(TRUE);
                    };

   auto l19 = []() {
                 graphics_info_t g;
                 if (g.moving_atoms_asc) {
                    g.clear_up_moving_atoms_wrapper(); g.clear_gl_rama_plot();
                 } else {
                    unfullscreen();
                 }
                 return gboolean(TRUE);
              };

   auto l20 = []() { graphics_info_t g; g.eigen_flip_active_residue(); return gboolean(TRUE); };

   auto l21 = []() { graphics_info_t g; g.try_label_unlabel_active_atom(); return gboolean(TRUE); };

   auto l22 = []() {
                  graphics_info_t g;
                  g.setup_draw_for_particles();
                  return gboolean(TRUE);
             };

   // boids
   auto l23 = [] () {
      graphics_info_t g;

      if (false) {
         if (! graphics_info_t::do_tick_boids)
            graphics_info_t::do_tick_boids = true;
         else
            graphics_info_t::do_tick_boids = false;
         g.setup_draw_for_boids();
         if (! graphics_info_t::do_tick_boids)
            std::cout << "--------- key press ----------- do_tick_boids "
                      << graphics_info_t::do_tick_boids << std::endl;
      }
      return gboolean(TRUE);
   };

   auto l24 = [] () {
                 // using the C API
                 // do_add_terminal_residue(1); // waits for user click :-)
                 graphics_info_t g;
                 g.add_terminal_residue_using_active_atom();
                 return gboolean(TRUE);
      };

   auto l25 = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    int imol_map = g.imol_refinement_map;
                    if (residue_p) {
                       mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].fill_partial_residue(residue_spec, g.Geom_p(), imol_map);

                       // now refine that
                       int saved_state = g.refinement_immediate_replacement_flag;
                       g.refinement_immediate_replacement_flag = 1;
                       std::string alt_conf("");
                       std::vector<mmdb::Residue *> rs = { residue_p };
                       g.refine_residues_vec(imol, rs, alt_conf, mol);
                       g.conditionally_wait_for_refinement_to_finish();
                       g.accept_moving_atoms();
                       g.refinement_immediate_replacement_flag = saved_state;
                    }
                 }
                 return gboolean(TRUE);
              };


   auto l26 = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].delete_residue_sidechain(residue_spec);
                    }
                 }
                 return gboolean(TRUE);
              };

   auto l28 = [] () {

                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       std::string this_chain_id = residue_p->GetChainID();
                       coot::residue_spec_t residue_spec(residue_p);
                       std::vector<std::vector<std::string> > ghost_chains_sets = molecules[imol].ncs_ghost_chains();
                       unsigned int n_ghost_chain_sets = ghost_chains_sets.size();
                       for (unsigned int i=0; i<n_ghost_chain_sets; i++) {
                          const std::vector<std::string> &chain_ids = ghost_chains_sets[i];
                          if (std::find(chain_ids.begin(), chain_ids.end(), this_chain_id) != chain_ids.end()) {
                             unsigned int idx_next = 0;
                             for (unsigned int j=0; j<chain_ids.size(); j++) {
                                if (chain_ids[j] == this_chain_id) {
                                   idx_next = j + 1;
                                   if (idx_next == chain_ids.size())
                                      idx_next = 0;
                                   break;
                                }
                             }
                             std::string chain_id_next = chain_ids[idx_next];
                             clipper::Coord_orth current_position = coot::co(at);
                             bool forward_flag = true;
                             glm::mat4 quat_mat = glm::toMat4(glm_quat);
                             clipper::Mat33<double> current_view_mat = glm_to_mat33(quat_mat);

                             if (molecules[imol].ncs_ghosts_have_rtops_p() == 0)
                                molecules[imol].fill_ghost_info(1, ncs_homology_level);

                             std::pair<bool, clipper::RTop_orth> new_ori =
                                molecules[imol].apply_ncs_to_view_orientation(current_view_mat,
                                                                              current_position,
                                                                              this_chain_id, chain_id_next,
                                                                              forward_flag);
                             if (new_ori.first) {
                                coot::util::quaternion q(new_ori.second.rot());
                                glm::quat q_ncs = coot_quaternion_to_glm(q);
                                glm_quat = glm::normalize(glm_quat * q_ncs); // wrong
                                clipper::Coord_orth t(new_ori.second.trn());
                                set_rotation_centre(t);

                                glm::quat q_ncs_1 = glm::rotate(q_ncs,   3.1415926f, glm::vec3(1,0,0));
                                glm::quat q_ncs_2 = glm::rotate(q_ncs_1, 3.1415926f, glm::vec3(1,0,0));
                                glm::quat q_ncs_3 = glm::inverse(q_ncs_2);

                                coot::util::quaternion cq = glm_to_coot_quaternion(q_ncs_3);

                                std::cout << "debug q_ncs  : " << glm::to_string(q_ncs)   << std::endl;
                                std::cout << "debug q_ncs_1: " << glm::to_string(q_ncs_1) << std::endl;
                                std::cout << "debug q_ncs_2: " << glm::to_string(q_ncs_2) << std::endl;
                                std::cout << "debug q_ncs_3: " << glm::to_string(q_ncs_3) << std::endl;
                                std::cout << "before: " << q << " after " << cq << std::endl;

                                graphics_info_t g;
                                g.update_things_on_move();

                             }
                             break;
                          }
                       }
                    } else {
                       std::cout << "ERROR:: no residue" << std::endl;
                    }
                 }
                 graphics_draw();
                 return gboolean(TRUE);
              };

   auto l29 = [] () {

                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
                 if (pp.first) {
                    g.update_mesh_for_outline_of_active_residue(pp.second.first, pp.second.second);
                    if (! tick_function_is_active()) {
                       int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
                    }
                    outline_for_active_residue_frame_count = 30;
                    do_tick_outline_for_active_residue = true;
                 }
                 return gboolean(TRUE);
              };

   auto l31 = [] () {
                 graphics_info_t g;
                 g.decrease_clipping_front();
                 return gboolean(TRUE);
              };

   auto l32 = [] () {
                 graphics_info_t g;
                 g.increase_clipping_front();
                 return gboolean(TRUE);
              };

   auto l33 = [] () {
                 graphics_info_t g;
                 g.decrease_clipping_back();
                 return gboolean(TRUE);
              };

   auto l34 = [] () {
                 graphics_info_t g;
                 g.increase_clipping_back();
                 return gboolean(TRUE);
              };

   // Note to self, Space and Shift Space are key *Release* functions

   std::vector<std::pair<keyboard_key_t, key_bindings_t> > kb_vec;
   // kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_d,      key_bindings_t(l1, "increase clipping")));
   kb_vec.push_back(std::make_pair(GDK_KEY_d, key_bindings_t(l13r, "step right")));
   kb_vec.push_back(std::make_pair(GDK_KEY_a, key_bindings_t(l13l, "step left")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_f,      key_bindings_t(l2, "decrease clipping")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_g,      key_bindings_t(l5, "go to blob")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_i,      key_bindings_t(l6, "spin")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_plus,   key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_equal,  key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_minus,  key_bindings_t(l7, "decrease contour level")));
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
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_l,      key_bindings_t(l21, "Label/Unlabel Active Atom")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_q,      key_bindings_t(l22, "Particles")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_b,      key_bindings_t(l23, "Murmuration")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_y,      key_bindings_t(l24, "Add Terminal Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_k,      key_bindings_t(l25, "Fill Partial Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_K,      key_bindings_t(l26, "Delete Sidechain")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l28, "NCS Other Chain")));

   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_space,  key_bindings_t(l18_space, "Accept Moving Atoms")));

   // clipping
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_1,      key_bindings_t(l31, "Clipping Front Expand")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_2,      key_bindings_t(l32, "Clipping Front Reduce")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_3,      key_bindings_t(l33, "Clipping Back Expand")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_4,      key_bindings_t(l34, "Clipping Back Reduce")));

   // control
   // meh - ugly and almost useless. Try again.
   // kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Control_L, key_bindings_t(l29, "Highlight Active Residue")));

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

   auto ldr = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       // for this to work I need to move setup_delete_item_pulse() into
                       // graphics_info_t. Not today.
                       g.setup_delete_item_pulse(residue_p);
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].delete_residue(residue_spec);
                    }
                 }
                 return gboolean(TRUE);
              };
   key_bindings_t delete_residue_key_binding(ldr, "Delete Residue");
   std::pair<keyboard_key_t, key_bindings_t> pdel(keyboard_key_t(GDK_KEY_d, true), delete_residue_key_binding);
   kb_vec.push_back(pdel);


   // ctrl left
   auto lc4 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed. No need to test it again here.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Left);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Left);
                 } else {
                    keypad_translate_xyz(1, 1);
                 }
                 return gboolean(TRUE);
              };

   // ctrl right
   auto lc5 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Right);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Right);
                 } else {
                    keypad_translate_xyz(1, -1);
                 }
                 return gboolean(TRUE);
              };

   // ctrl up
   auto lc6 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Up);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Up);
                 } else {
                    keypad_translate_xyz(2, 1);
                 }
                 return gboolean(TRUE);
              };
   // ctrl down
   auto lc7 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Down);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Down);
                 } else {
                    keypad_translate_xyz(2, -1);
                 }
                 return gboolean(TRUE);
              };

   auto lc_qsa = [] () {
                    graphics_info_t g;
                    g.quick_save();
                    return gboolean(TRUE);
                 };

   key_bindings_t ctrl_arrow_left_key_binding(lc4, "R/T Left");
   key_bindings_t ctrl_arrow_right_key_binding(lc5, "R/T Right");
   key_bindings_t ctrl_arrow_up_key_binding(lc6, "R/T Up");
   key_bindings_t ctrl_arrow_down_key_binding(lc7, "R/T Down");
   key_bindings_t ctrl_eigen_flip(l20, "Eigen-Flip");
   key_bindings_t ctrl_quick_save(lc_qsa, "Quick Save");

   std::pair<keyboard_key_t, key_bindings_t> p4(keyboard_key_t(GDK_KEY_Left,  true), ctrl_arrow_left_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p5(keyboard_key_t(GDK_KEY_Right, true), ctrl_arrow_right_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p6(keyboard_key_t(GDK_KEY_Up,    true), ctrl_arrow_up_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p7(keyboard_key_t(GDK_KEY_Down,  true), ctrl_arrow_down_key_binding);
   std::pair<keyboard_key_t, key_bindings_t> p8(keyboard_key_t(GDK_KEY_e,     true), ctrl_eigen_flip);
   std::pair<keyboard_key_t, key_bindings_t> p9(keyboard_key_t(GDK_KEY_s,     true), ctrl_quick_save);

   kb_vec.push_back(p4);
   kb_vec.push_back(p5);
   kb_vec.push_back(p6);
   kb_vec.push_back(p7);
   kb_vec.push_back(p8);
   kb_vec.push_back(p9);

   std::vector<std::pair<keyboard_key_t, key_bindings_t> >::const_iterator it;
   for (it=kb_vec.begin(); it!=kb_vec.end(); it++)
     g.key_bindings_map[it->first] = it->second;

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
      if (direction == 1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
      if (direction == -1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;

      // std::cout << "INFO:: contour level for map " << imol_scroll << " is "
      // << molecules[imol_scroll].contour_level << std::endl;
      set_density_level_string(imol_scroll, molecules[imol_scroll].contour_level);
      display_density_level_this_image = 1;

      graphics_draw(); // queue
   }
}

#include "widget-from-builder.hh"

//static
void
graphics_info_t::fullscreen() {

   GtkWidget *window = widget_from_builder("main_window");

   if (GTK_IS_WINDOW(window)) {

      GtkWidget *vbox       = widget_from_builder("main_window_vbox");
      GtkWidget *overlay    = widget_from_builder("main_window_graphics_overlay");
      GtkWidget *status_bar = widget_from_builder("main_window_statusbar");
      GtkWidget *tool_bar   = widget_from_builder("main_window_toolbar");
      // GtkWidget *tool_bar_frame   = widget_from_builder("main_window_model_fit_dialog_frame");
      GtkWidget *menu_bar_frame   = widget_from_builder("main_window_menubar");

      // std::cout << "hiding menu_bar_frame " << menu_bar_frame << std::endl;
      gtk_widget_hide(menu_bar_frame);
      // std::cout << "hiding tool_bar " << tool_bar << std::endl;
      gtk_widget_hide(tool_bar);
      gtk_window_fullscreen(GTK_WINDOW(window));

      gtk_container_remove(GTK_CONTAINER(vbox), status_bar);
      gtk_overlay_add_overlay(GTK_OVERLAY(overlay), status_bar);
      // gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), status_bar, TRUE);
      gtk_widget_set_halign(status_bar, GTK_ALIGN_START);
      gtk_widget_set_valign(status_bar, GTK_ALIGN_END);
      gtk_widget_grab_focus(glareas[0]);

      if (false) {
         gtk_container_remove(GTK_CONTAINER(vbox), tool_bar);
         gtk_overlay_add_overlay(GTK_OVERLAY(overlay), tool_bar);
         gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), tool_bar, TRUE);
         gtk_widget_set_halign(tool_bar, GTK_ALIGN_START);
         gtk_widget_set_valign(tool_bar, GTK_ALIGN_START);
      }

      if (false) {
         gtk_container_remove(GTK_CONTAINER(vbox), menu_bar_frame);
         gtk_overlay_add_overlay(GTK_OVERLAY(overlay), menu_bar_frame);
         gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), menu_bar_frame, TRUE);
         gtk_widget_set_halign(menu_bar_frame, GTK_ALIGN_START);
         gtk_widget_set_valign(menu_bar_frame, GTK_ALIGN_START);
      }

      graphics_info_t g;
      g.add_status_bar_text(""); // clear it
   }
}

//static
void
graphics_info_t::unfullscreen() {

   GtkWidget *window = widget_from_builder("main_window");
   if (GTK_IS_WINDOW(window)) {
      gtk_window_unfullscreen(GTK_WINDOW(window));

      // unoverlay and reparent here.

      GtkWidget *vbox       = widget_from_builder("main_window_vbox");
      GtkWidget *overlay    = widget_from_builder("main_window_graphics_overlay");
      GtkWidget *status_bar = widget_from_builder("main_window_statusbar");
      GtkWidget *tool_bar   = widget_from_builder("main_window_toolbar");
      GtkWidget *menu_bar   = widget_from_builder("main_window_menubar");
      GtkWidget *tool_bar_frame   = widget_from_builder("main_window_model_fit_dialog_frame");

#if (GTK_MAJOR_VERSION > 3)
      gtk_overlay_remove_overlay(GTK_OVERLAY(overlay), tool_bar);
      gtk_container_add(GTK_CONTAINER(vbox), tool_bar);

      gtk_overlay_remove_overlay(GTK_OVERLAY(overlay), status_bar);
      gtk_container_add(GTK_CONTAINER(vbox), status_bar);
#endif

      std::cout << "show widget " << menu_bar << std::endl;
      std::cout << "show widget " << tool_bar << std::endl;

      gtk_widget_show(menu_bar);
      gtk_widget_show(tool_bar);
      gtk_widget_show(tool_bar_frame);
   }
}
