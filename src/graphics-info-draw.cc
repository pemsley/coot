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
   if (err) std::cout << "init_screen_quads() err is " << err << std::endl;

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
graphics_info_t::get_molecule_mvp(bool debug_matrices) {

   // presumes that we are in the correct programID

   float w = static_cast<float>(graphics_info_t::graphics_x_size);
   float h = static_cast<float>(graphics_info_t::graphics_y_size);
   float screen_ratio = static_cast<float>(w)/static_cast<float>(h);

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

   if (graphics_info_t::perspective_projection_flag) {

      // for fun/testing
      // turn off view scaling when tinkering with this?
      // there should not be a concept of "zoom" with perspective view, just translation
      // along screen-Z.

      glm::mat4 trackball_matrix = glm::toMat4(graphics_info_t::glm_quat);

      glm::vec3 ep = eye_position; // in view space i.e. (0,0,z) (z = 40, say)
      glm::vec3 up(0,1,0);
      glm::vec3 origin(0,0,0);

      model_matrix = glm::mat4(1.0);
      model_matrix = glm::translate(model_matrix, -rc);
      model_matrix = trackball_matrix * model_matrix;

      view_matrix = glm::lookAt(ep, origin, up);

      float fov = 50.0/zoom; // degrees
      if (fov > 50.0) fov = 50.0;
      fov = 35.0;
      fov = 30.0;

      glm::mat4 projection_matrix_persp = glm::perspective(glm::radians(fov),
                                                           screen_ratio,
                                                           screen_z_near_perspective,
                                                           screen_z_far_perspective);
      projection_matrix = projection_matrix_persp; // for debugging below
      mvp = projection_matrix_persp * view_matrix * model_matrix;
   }

   // bool debug_matrices = false;
   if (debug_matrices) {
      std::cout << "model, view, projection, mvp" << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(model_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(view_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(projection_matrix) << std::endl;
      std::cout << "get_molecule_mvp: " << glm::to_string(mvp) << std::endl;
   }


   return mvp;
}

// static
glm::mat4
graphics_info_t::get_particle_mvp() {

   float w = static_cast<float>(graphics_x_size);
   float h = static_cast<float>(graphics_y_size);
   float screen_ratio = static_cast<float>(w)/static_cast<float>(h);
   float sr = screen_ratio;

   glm::mat4 model_matrix(1.0);

   float z = zoom * 0.04;
   GLfloat near = -0.1 * zoom * clipping_front;
   GLfloat far  =  0.3 * zoom * clipping_back;

   glm::mat4 proj = glm::ortho(-0.3f*zoom*sr, 0.3f*zoom*sr,
                               -0.3f*zoom,    0.3f*zoom,
                               near, far);

   glm::vec3 rc = get_rotation_centre();
   glm::mat4 view_matrix = glm::toMat4(glm_quat);
   // view_matrix = glm::translate(view_matrix, -rc);
   if (true) {
      glm::mat3 rot_mat(view_matrix);
      if (false)
         std::cout << "   view " << glm::to_string(view_matrix)
                   << "\nrot_mat " << glm::to_string(rot_mat) << std::endl;
      glm::mat3 tt = glm::transpose(rot_mat);
      glm::vec4 trans = view_matrix[3];
      view_matrix = glm::mat4(1.0f);
      // view_matrix = glm::mat4(tt);
      // view_matrix = glm::translate(view_matrix, rc);
   }

   glm::mat4 mvp = proj * view_matrix * model_matrix;

   if (perspective_projection_flag) {

      float fov = 35.0; // degrees, the smaller this value, the more we seem to be
                        // "zoomed in."
      float screen_ratio = 1.0;
      float z_front = 0.1;
      float z_back = 160;
      // glm::mat4 quat_mat(mouse_button_info.quat);
      glm::mat4 quat_mat = glm::toMat4(glm_quat);

      glm::vec3 ep_ori = get_world_space_eye_position();
      glm::vec3 rc = get_rotation_centre();
      glm::vec3 ep = ep_ori + rc; // is this right?
      glm::vec3 up = get_camera_up_direction(quat_mat);
      view_matrix = glm::lookAt(ep, rc, up);
      for (unsigned int i=0; i<3; i++)
         for (unsigned int j=0; j<3; j++)
            view_matrix[i][j] = 0.0f;
      for (unsigned int j=0; j<3; j++)
         view_matrix[j][j] = 1.0f;
      proj = glm::perspective(glm::radians(fov), screen_ratio, z_front, z_back);
      mvp = proj * view_matrix * model_matrix;
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

      glm::mat4 vr = get_view_rotation();
      glm::vec4 rot_test_vector_1 = glm::vec4(test_vector_1, 1.0) * vr;
      glm::vec4 rot_test_vector_2 = glm::vec4(test_vector_2, 1.0) * vr;

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


glm::mat4
graphics_info_t::get_view_rotation() {

   // need to be in the correct program (well, the model-drawing part)

   return glm::toMat4(graphics_info_t::glm_quat);

}


void
graphics_info_t::setup_map_uniforms(Shader *shader_p,
                                    const glm::mat4 &mvp,
                                    const glm::mat4 &view_rotation,
                                    float density_surface_opacity) {

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   GLenum err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniformMatrix4fv() mvp " << err << std::endl;
   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
   err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniformMatrix4fv() vr  " << err << std::endl;

   GLuint background_colour_uniform_location = shader_p->background_colour_uniform_location;
   glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
   glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
   err = glGetError();
   if (err) std::cout << "   setup_map_uniforms() glUniform4fv() for bg  " << err << std::endl;

   // GLuint opacity_uniform_location = shader.map_opacity_uniform_location;
   shader_p->set_float_for_uniform("map_opacity", density_surface_opacity);
   err = glGetError(); if (err) std::cout << "   setup_map_uniforms() glUniformf() for opacity "
                                          << err << std::endl;

   GLuint eye_position_uniform_location = shader_p->eye_position_uniform_location;
   glm::vec4 ep = glm::vec4(get_world_space_eye_position(), 1.0);
   glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
   err = glGetError(); if (err) std::cout << "   setup_map_uniforms() glUniform4fv() for eye position "
                                          << err << std::endl;

}


void
graphics_info_t::draw_map_molecules(bool draw_transparent_maps) {

   // run through this molecule loop twice - for opaque then transparent maps
   // first, a block that decides if we need to do anything.

   bool needs_blend_reset = false;

   //

   bool cosine_dependent_map_opacity = true;

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

   if (cosine_dependent_map_opacity) {
      needs_blend_reset = true;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   if (!draw_transparent_maps || n_transparent_maps > 0) {

      glLineWidth(map_line_width * framebuffer_scale);
      GLenum err = glGetError();
      if (err) std::cout << "gtk3_draw_map_molecules() glLineWidth " << err << std::endl;

      shader_for_maps.Use();

      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation = get_view_rotation();

      glEnable(GL_DEPTH_TEST); // this needs to be in the draw loop!?
      glDepthFunc(GL_LESS);

      Shader &shader = shader_for_maps;

      glm::vec4 ep(get_world_space_eye_position(), 1.0);
      glm::vec3 ep3 = ep/ep.w;

      for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {

         molecule_class_info_t &m = graphics_info_t::molecules[ii]; // not const because shader changes
         if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
         if (! m.draw_it_for_map) continue;
         if (draw_transparent_maps) {
            if (m.is_an_opaque_map())
               continue; // not this round
         } else {
            // only draw (completely) opaque (that's what the question means)
            if (! m.is_an_opaque_map())
               continue;
         }

         if (m.n_vertices_for_map_VertexArray > 0) {

            err = glGetError();
            if (err) std::cout << "draw_map_molecules() --- draw map loop start --- error "
                               << std::endl;

            bool draw_with_lines = true;
            if (!m.draw_it_for_map_standard_lines) draw_with_lines = false;

            //glUniform1i(shader.is_perspective_projection_uniform_location,
            // graphics_info_t::perspective_projection_flag);
            shader.set_bool_for_uniform("is_perspective_projection", perspective_projection_flag);
            err = glGetError(); if (err) std::cout << "   draw_map_molecules() error B " << std::endl;

            shader.set_bool_for_uniform("do_depth_fog", graphics_info_t::shader_do_depth_fog_flag);
            shader.set_bool_for_uniform("do_diffuse_lighting", true);
            shader.set_float_for_uniform("shininess", m.shader_shininess);
            shader.set_float_for_uniform("specular_strength", m.shader_specular_strength);

            // --- lights ----

            std::map<unsigned int, lights_info_t>::const_iterator it; // iterate over the lights map
            for (it=lights.begin(); it!=lights.end(); it++) {
               unsigned int light_idx = it->first;
               const lights_info_t &light = it->second;
               shader.setup_light(light_idx, light, view_rotation);
            }

            // --- material ---

            Material &material = m.material_for_maps;
            shader.set_vec4_for_uniform( "material.ambient",   material.ambient);
            shader.set_vec4_for_uniform( "material.diffuse",   material.diffuse);
            shader.set_vec4_for_uniform( "material.specular",  material.specular);
            shader.set_float_for_uniform("material.shininess", material.shininess);
            shader.set_float_for_uniform("material.specular_strength", material.specular_strength);

            // --- background ---

            GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
            glm::vec4 bgc(background_colour, 1.0);
            glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniform4fv() for bg  " << err << std::endl;

            // --- fresnel ---

            shader.set_bool_for_uniform("do_fresnel",     m.fresnel_settings.state);
            shader.set_float_for_uniform("fresnel_bias",  m.fresnel_settings.bias);
            shader.set_float_for_uniform("fresnel_scale", m.fresnel_settings.scale);
            shader.set_float_for_uniform("fresnel_power", m.fresnel_settings.power);

            // --- draw ---

            if (draw_with_lines) {
               // I don't see why this is needed - but it is.
               if (! m.is_an_opaque_map())
                  glEnable(GL_BLEND);

               glBindVertexArray(m.m_VertexArrayID_for_map);
               err = glGetError();
               if (err) std::cout << "   draw_map_molecules() glBindVertexArray() "
                                  << m.m_VertexArrayID_for_map
                                  << " with GL err " << err << std::endl;

               glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_map_lines_ID);

               setup_map_uniforms(&shader, mvp, view_rotation, m.density_surface_opacity);
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
                  m.sort_map_triangles(eye_pos_co);
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

               // opacity:
               // GLuint opacity_uniform_location = graphics_info_t::shader_for_maps.map_opacity_uniform_location;
               // float opacity = m.density_surface_opacity;
               // glUniform1f(opacity_uniform_location, opacity);
               shader.set_float_for_uniform("map_opacity", m.density_surface_opacity);
               err = glGetError(); if (err) std::cout << "   draw_map_molecules() glUniformf() for opacity "
                                                      << err << std::endl;
               // cosine_dependent_map_opacity
               shader.set_bool_for_uniform("cosine_dependent_map_opacity", cosine_dependent_map_opacity);

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

      molecule_class_info_t &m = graphics_info_t::molecules[ii];
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

         // glUniform1i(shader.is_perspective_projection_uniform_location,
         // graphics_info_t::perspective_projection_flag);
         shader.set_bool_for_uniform("is_perspective_projection", perspective_projection_flag);

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
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for view_rotation "
                            << err << std::endl;

         // std::cout << glm::to_string(mvp) << std::endl;
         // std::cout << glm::to_string(view_rotation) << std::endl;

         GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
         glm::vec4 bgc(background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec3 ep(get_world_space_eye_position());
         glUniform3fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform3fv() for eye position " << err << std::endl;

         shader.set_bool_for_uniform("do_depth_fog", shader_do_depth_fog_flag);
         shader.set_bool_for_uniform("do_diffuse_lighting", true); // false for demo c.f. old style graphics

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() pre-lights glDrawElements() "
                            << shader.name << " with err " << err << std::endl;

         // --- lights ----

         std::map<unsigned int, lights_info_t>::const_iterator it; // iterate over the lights map
         for (it=lights.begin(); it!=lights.end(); it++) {
            unsigned int light_idx = it->first;
            const lights_info_t &light = it->second;
            shader.setup_light(light_idx, light, view_rotation);
         }
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() post-lights "
                            << shader.name << " with err " << err << std::endl;


         // --- material ---

         const Material &material = m.material_for_models;
         shader.set_vec4_for_uniform( "material.ambient",   material.ambient);
         shader.set_vec4_for_uniform( "material.diffuse",   material.diffuse);
         shader.set_vec4_for_uniform( "material.specular",  material.specular);
         shader.set_float_for_uniform("material.shininess", material.shininess);
         shader.set_float_for_uniform("material.specular_strength", material.specular_strength);

         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() post-material "
                            << shader.name << " with err " << err << std::endl;

         // draw with the vertex count, not the index count.
         GLuint n_verts = graphics_info_t::molecules[ii].n_indices_for_model_triangles;

         // std::cout << "   Drawing " << n_verts << " model vertices" << std::endl;
         err = glGetError();
         if (err) std::cout << "   error pre draw_model_molecules() glDrawElements() " << shader.name
                            << " n_vertices " << n_verts << " with GL err " << err << std::endl;
         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glDrawElements() " << shader.name
                            << " n_vertices " << n_verts << " with GL err " << err << std::endl;

         draw_molecule_atom_labels(m, mvp, view_rotation);

      }
   }
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

   // this doesn't seem sensibly arranged.
   glm::vec3 unused(0,0,0);
   m.draw_atom_labels(brief_atom_labels_flag,
                      seg_ids_in_atom_labels_flag,
                      label_colour,
                      mvp, view_rotation, unused);

   glDisable(GL_BLEND);

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

      if (m.n_vertices_for_model_VertexArray > 0) {

         glDisable(GL_BLEND); // stop semi-transparent bonds - but why do we have them?
         gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));

         shader.Use();
         GLuint err = glGetError();
         if (err) std::cout << "   error draw_intermediate_atoms() glUseProgram() "
                            << err << std::endl;

         glBindVertexArray(m.m_VertexArray_for_model_ID);
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glBindVertexArray() "
                            << m.m_VertexArray_for_model_ID
                            << " with GL err " << err << std::endl;

         // should not be needed?
         glBindBuffer(GL_ARRAY_BUFFER, m.m_VertexBuffer_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_intermediate_atoms() glBindBuffer() v "
                            << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.m_IndexBuffer_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_intermediate_atoms() glBindBuffer() i "
                            << err << std::endl;

         GLuint mvp_location           = shader.mvp_uniform_location;
         GLuint view_rotation_location = shader.view_rotation_uniform_location;

         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniformMatrix4fv() pre mvp "
                            << err << std::endl;
         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniformMatrix4fv() for mvp "
                            << err << std::endl;
         glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniformMatrix4fv() "
                            << "for view_rotation " << err << std::endl;

         GLuint background_colour_uniform_location = shader.background_colour_uniform_location;
         glm::vec4 bgc(background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "error draw_model_molecules() glUniform4fv() for background "
                            << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec3 ep = get_world_space_eye_position();
         glUniform3fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "error draw_intermediate_atoms() glUniform4fv() for eye position "
                            << err << std::endl;

         shader.set_bool_for_uniform("do_depth_fog", shader_do_depth_fog_flag);

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


#include "Instanced-Markup-Mesh.hh"

void
graphics_info_t::setup_rama_balls() {

   rama_balls_mesh.setup_octasphere(2);
   rama_balls_mesh.setup_instancing_buffers(3000);
}

void
graphics_info_t::update_rama_balls(std::vector<Instanced_Markup_Mesh_attrib_t> *balls) {

   std::pair<std::vector<vertex_with_rotation_translation>, std::vector<g_triangle> > v;

   auto rr = saved_dragged_refinement_results;

   balls->clear();

   float rama_ball_pos_offset_scale = 0.6;
   glm::vec3 screen_up_dir(0.2, 0.3, 0.3);

   unsigned int n_atoms_found = 0;

   // std::cout << "update rama ball for " << rr.all_ramas.size() << " balls " << std::endl;
   for (unsigned int i=0; i<rr.all_ramas.size(); i++) {

      float rama_score = rr.all_ramas[i].distortion;

      const coot::atom_spec_t &spec_CA = rr.all_ramas[i].atom_spec_CA;
      mmdb::Atom *at = spec_CA.get_atom(moving_atoms_asc->mol);
      if (at) {

         float d = rr.all_ramas[i].distortion;
         glm::vec3 atom_position(at->x, at->y, at->z);
         glm::vec3 ball_position(rr.all_ramas[i].ball_pos_x,
                                 rr.all_ramas[i].ball_pos_y,
                                 rr.all_ramas[i].ball_pos_z);
         // std::cout << "update_rama_balls() " << i << " " << spec_CA << " " << d << std::endl;
         float size = 0.38;
         // std::cout << "debug d " << d << std::endl;
         float ra = hud_geometry_distortion_to_rotation_amount_rama(d);
         coot::colour_t cc(0.1, 0.9, 0.2);
         cc.rotate(ra);
         glm::vec4 col = cc.to_glm();
         Instanced_Markup_Mesh_attrib_t ball(col, ball_position, size);
         float d1 = d + 70.0;
         float d2 = - d1 * 0.01;
         if (d2 < 0.0) d2 = 0.0;
         if (d2 > 1.0) d2 = 1.0;
         ball.specular_strength = 0.04 + d2;
         ball.shininess = 2.0 + 65.0 * d2;
         balls->push_back(ball);
      }
   }
}



void
graphics_info_t::draw_intermediate_atoms_rama_balls() {

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();

   Shader &shader = graphics_info_t::shader_for_rama_balls;
   if (true) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 view_rotation = get_view_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      bool do_depth_fog = true;
      // this is a bit ugly

      // the balls are updated after a refinement cycle has finished -
      // no need to do it here
      // graphics_info_t g;
      // std::vector<Instanced_Markup_Mesh_attrib_t> balls;
      // update_rama_balls(&balls);
      // rama_balls_mesh.update_instancing_buffers(balls);
      rama_balls_mesh.draw(&shader, mvp, view_rotation, lights, eye_position, bg_col, do_depth_fog);
   }

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

         if (n_atom_pulls > 0) { // class variable now.

            Shader &shader = shader_for_models;
            shader.Use();
            GLuint err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glUseProgram() "
                               << err << std::endl;

            glBindVertexArray(m_VertexArray_for_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindVertexArray()"
                               << " with GL err " << err << std::endl;
            glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer_for_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindBuffer()"
                               << " with GL err " << err << std::endl;

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_atom_pull_restraints_ID);
            err = glGetError();
            if (err) std::cout << "   error draw_atom_pull_restraints() glBindBuffer() for index"
                               << " with GL err " << err << std::endl;

            // Uniforms
            glm::mat4 mvp = get_molecule_mvp();
            glm::mat4 view_rotation = get_view_rotation();
            GLuint mvp_location = shader.mvp_uniform_location;

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
         glm::mat4 mvp_particle = get_particle_mvp();
         glm::mat4 mvp = get_molecule_mvp();
         mesh_for_particles.draw_particles(&shader_for_particles, mvp_particle);
      }
   }

}

void
graphics_info_t::draw_molecules() {

   // opaque things

   draw_intermediate_atoms();

   draw_intermediate_atoms_rama_balls();

   draw_atom_pull_restraints();

   draw_meshes(); // get a better name

   draw_map_molecules(false); // transparency

   draw_unit_cells();

   draw_environment_graphics_object();

   draw_generic_objects();

   draw_boids();

   // this is the last opaque thing to be drawn because the atom labels are blended.
   // It should be easy to break out the atom label code into its own function. That
   // might be better.
   //
   draw_model_molecules();

   // transparent things...

   draw_particles();

   draw_map_molecules(true);

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
            glm::mat4 view_rotation = get_view_rotation();
            glm::vec4 bg_col(background_colour, 1.0);

            bool do_depth_fog = shader_do_depth_fog_flag;
            mesh_for_environment_distances.mesh.draw(&shader_for_moleculestotriangles,
                                                     mvp, view_rotation,
                                                     lights, eye_position, bg_col,
                                                     do_depth_fog);

            Shader *shader_p = &shader_for_atom_labels;

            GLenum err = glGetError();
            if (err) std::cout << "error draw_environment_graphics_object() before labela err "
                               << err << std::endl;

            if (! labels.empty()) {
               for (unsigned int i=0; i<labels.size(); i++) {
                  const std::string &label  = labels[i].label;
                  const glm::vec3 &position = labels[i].position;
                  const glm::vec4 &colour   = labels[i].colour;
                   tmesh_for_labels.draw_atom_label(label, position, colour, shader_p,
                                                    mvp, view_rotation, lights, eye_position,
                                                    bg_col, do_depth_fog,
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


void
graphics_info_t::draw_unit_cells() {

   glm::mat4 mvp = get_molecule_mvp();
   for (int ii=n_molecules()-1; ii>=0; ii--) {
      molecule_class_info_t &m = molecules[ii];
      m.draw_unit_cell(&shader_for_lines, mvp);
   }

}

void
graphics_info_t::draw_meshes() {

   bool draw_meshes = true;
   bool draw_mesh_normals = false;

   glm::vec3 eye_position = get_world_space_eye_position();
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation();
   glm::vec4 bg_col(background_colour, 1.0);

   bool do_depth_fog = true;

   //std::cout << "mvp diag "
   // << mvp[0][0] << " " << mvp[1][1] << " " << mvp[2][2] << std::endl;

   glm::mat3 vrm(glm::toMat4(graphics_info_t::glm_quat));
   glm::mat3 vrmt = glm::transpose(vrm);
   glm::mat3 p = vrmt * vrm;

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
         glDisable(GL_BLEND);
         for (int ii=n_molecules()-1; ii>=0; ii--) {
            molecule_class_info_t &m = molecules[ii]; // not const because the shader changes
            for (unsigned int jj=0; jj<m.meshes.size(); jj++) {
               m.meshes[jj].draw(&shader_for_moleculestotriangles, mvp,
                                 view_rotation, lights, eye_position, bg_col, do_depth_fog);
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
      float s = rotation_centre_cube_size;
      glm::vec3 sc(s,s,s);
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
      if (err) std::cout << "   error draw_central_cube() glUniformMatrix4fv() for view_rotation " << err
                         << std::endl;

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

   lights_info_t light;
   light.position = glm::vec4(-2.0f, 2.0f, 5.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(0.5, 0, 1.0));
   graphics_info_t::lights[0] = light;

   light.position = glm::vec4(3.0f, -2.0f, 4.0f, 1.0f);
   light.direction = glm::normalize(glm::vec3(-1.0, 0.5, 1.0));
   graphics_info_t::lights[1] = light;
}

void
graphics_info_t::setup_hud_geometry_bars() {

   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   shader_for_hud_geometry_bars.Use();

   mesh_for_hud_geometry.setup_camera_facing_quad_for_bar();
   mesh_for_hud_geometry.setup_instancing_buffer(100);

   // If not found in this directory, then try default directory.
   texture_for_hud_geometry_labels.set_default_directory(coot::package_data_dir());
   texture_for_hud_geometry_labels.init("hud-label-nbc-rama.png");

   // Do I need to Use() the shader_for_hud_geometry_labels here?
   shader_for_hud_geometry_labels.Use();
   mesh_for_hud_geometry_labels.setup_quad();
   mesh_for_hud_geometry_labels.set_position_and_scale(glm::vec2(-0.96, 0.89), 0.037); // was 0.35
}

float
graphics_info_t::hud_geometry_distortion_to_bar_size_rama(float distortion) {
   float d1 = distortion + 200.0;
   float d2 = d1 * 0.0003;
   if (d2 < 0.0) d2 = 0.0;
   float d3 = 100.0 * d2 * d2;
   return d3;
}

float
graphics_info_t::hud_geometry_distortion_to_bar_size_nbc(float distortion) {
   return distortion * 0.002;
}


float
graphics_info_t::hud_geometry_distortion_to_rotation_amount_rama(float distortion) {
   distortion += 200.0;
   float rotation_amount = 1.0 - 0.0022 * distortion;
   if (rotation_amount < 0.68) rotation_amount = 0.68; // red cap
   if (rotation_amount > 1.0) rotation_amount = 1.0;
   return rotation_amount;
}


void
graphics_info_t::draw_hud_geometry_bars() {

   if (! moving_atoms_asc) return;
   if (! moving_atoms_asc->mol) return;

   glEnable(GL_DEPTH_TEST); // needed?
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   // first draw the text (labels) texture

   if (! saved_dragged_refinement_results.refinement_results_contain_overall_rama_plot_score) {
      // std::cout << "chopping the texture " << std::endl;
      mesh_for_hud_geometry_labels.setup_texture_coords_for_nbcs_only();
   } else {
      mesh_for_hud_geometry_labels.setup_texture_coords_for_nbcs_and_rama();
   }

   texture_for_hud_geometry_labels.Bind(0);
   mesh_for_hud_geometry_labels.draw(&shader_for_hud_geometry_labels);

   // now draw the bars

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
                       auto distortion_to_rotation_amount,
                       auto distortion_to_bar_size) {

                         glm::vec2 to_top_left(-0.91, 0.9 - 0.05 * static_cast<float>(bar_index));
                         float sum_l = 0;
                         int n = baddies.size();
                         for (int i=(n-1); i>=0; i--) {
                            coot::colour_t cc(0.1, 0.9, 0.2);
                            float d = baddies[i].second;
                            float rotation_amount = distortion_to_rotation_amount(d);
                            cc.rotate(rotation_amount);
                            glm::vec4 col = cc.to_glm();
                            col.w = 0.7;
                            glm::vec2 position_offset = to_top_left + glm::vec2(sum_l, 0.0);
                            float bar_length = distortion_to_bar_size(d);
                            if (bar_index == 111)
                               std::cout << "index i " << i << " distortion " << d
                                         << " bar length " << bar_length
                                         << " atom " << baddies[i].first
                                         << " position_offset " << glm::to_string(position_offset)  << std::endl;

                            HUD_bar_attribs_t bar(col, position_offset, bar_length);
                            new_bars_p->push_back(bar);
                            sum_l += bar_length + 0.005; // with a gap between bars
                          }
                       };

   std::vector<HUD_bar_attribs_t> new_bars;

   if (saved_dragged_refinement_results.refinement_results_contain_overall_rama_plot_score)
      add_bars(saved_dragged_refinement_results.sorted_rama_baddies, 1, &new_bars,
               hud_geometry_distortion_to_rotation_amount_rama, hud_geometry_distortion_to_bar_size_rama);

   if (saved_dragged_refinement_results.refinement_results_contain_overall_nbc_score)
      add_bars(saved_dragged_refinement_results.sorted_nbc_baddies, 0, &new_bars,
               distortion_to_rotation_amount_nbc, hud_geometry_distortion_to_bar_size_nbc);

   if (! new_bars.empty()) {
      // std::cout << "new bar size " << new_bars.size() << std::endl;
      mesh_for_hud_geometry.update_instancing_buffer_data(new_bars);
      Shader &shader = shader_for_hud_geometry_bars;
      mesh_for_hud_geometry.draw(&shader);
   }
   glDisable(GL_BLEND);

}

bool
graphics_info_t::check_if_hud_bar_clicked(double mouse_x, double mouse_y) {

   bool status = false;
   if (! moving_atoms_asc) return false;
   if (! moving_atoms_asc->mol) return false;

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

   auto distortion_to_bar_size_rama = [] (float distortion) {
                                         distortion += 200.0;
                                         float d2 = distortion * 0.0003;
                                         if (d2 < 0.0) d2 = 0.0;
                                         float d3 = 100.0 * d2 * d2;
                                         return d3;
                                      };

   auto distortion_to_bar_size_nbc = [] (float distortion) {
                                        return distortion * 0.002;
                                     };

   auto check_blocks = [mouse_in_opengl_coords] (const std::vector<std::pair<coot::atom_spec_t, float> > &baddies,
                                                 unsigned int bar_index,
                                                 auto distortion_to_bar_size) {

                          bool status = false;
                          glm::vec2 to_top_left(-0.91, 0.9 - 0.05 * static_cast<float>(bar_index));
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
                                   //           << " i " << i << " " << baddies[i].first << std::endl;
                                   float tiny_y_offset = -0.01; // not sure why I need this
                                   if (mouse_in_opengl_coords.y >= (to_top_left.y + tiny_y_offset)) {
                                      // 0.03 is the bar height in setup_camera_facing_quad()
                                      float bar_height = 0.03;
                                      if (mouse_in_opengl_coords.y <= (to_top_left.y+tiny_y_offset+bar_height)) {
                                         coot::atom_spec_t spec(baddies[i].first);
                                         if (moving_atoms_asc->mol) {
                                            mmdb::Atom *at = spec.get_atom(moving_atoms_asc->mol);
                                            if (at) {
                                               clipper::Coord_orth pt = coot::co(at);
                                               std::cout << "INFO: geom bar atom: " << coot::atom_spec_t(at)
                                                         << std::endl;
                                               set_rotation_centre(pt);
                                               status = true;
                                            }
                                         } else {
                                            std::cout << "ERROR:: no moving atoms mol" << std::endl;
                                         }
                                      }
                                   }
                                }
                             }
                          }
                          return status;
                       };

   if (saved_dragged_refinement_results.refinement_results_contain_overall_nbc_score)
      status = check_blocks(saved_dragged_refinement_results.sorted_nbc_baddies, 0, distortion_to_bar_size_nbc);

   if (saved_dragged_refinement_results.refinement_results_contain_overall_rama_plot_score)
      status = check_blocks(saved_dragged_refinement_results.sorted_rama_baddies, 1, distortion_to_bar_size_rama);

   return status;
}


gboolean
graphics_info_t::render(bool to_screendump_framebuffer, const std::string &output_file_name) {

   GtkGLArea *gl_area = GTK_GL_AREA(glareas[0]);
   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   GLenum err = glGetError();
   if (err) std::cout << "render() start " << err << std::endl;

   // is this needed? - does the context ever change?
   gtk_gl_area_make_current(gl_area);
   err = glGetError(); if (err) std::cout << "render() post gtk_gl_area_make_current() err " << err << std::endl;

   glViewport(0, 0, framebuffer_scale * w, framebuffer_scale * h);
   err = glGetError(); if (err) std::cout << "render() post glViewport() err " << err << std::endl;
   screen_framebuffer.bind();
   err = glGetError(); if (err) std::cout << "render() post screen_framebuffer bind() err " << err << std::endl;

   glEnable(GL_DEPTH_TEST);

   { //  ------------------- render scene ----------------------------
      const glm::vec3 &bg = graphics_info_t::background_colour;
      glClearColor (bg[0], bg[1], bg[2], 1.0); // what difference does this make?
      err = glGetError(); if (err) std::cout << "render() B err " << err << std::endl;
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      err = glGetError(); if (err) std::cout << "render() C err " << err << std::endl;

      draw_central_cube(gl_area);
      draw_origin_cube(gl_area);
      err = glGetError(); if (err) std::cout << "render()  pre-draw-text err " << err << std::endl;

      draw_molecules();

      draw_hud_geometry_bars();

      draw_identification_pulse();

      glBindVertexArray(0); // here is not the place to call this.
   }

   if (to_screendump_framebuffer) {

      glDisable(GL_DEPTH_TEST);
      unsigned int sf = framebuffer_scale;
      glViewport(0, 0, sf * w, sf *h);
      framebuffer screendump_framebuffer;
      unsigned int index_offset = 0;
      screendump_framebuffer.init(sf * w, sf * h, index_offset, "screendump");
      screendump_framebuffer.bind();
      render_scene_to_base_framebuffer();
      gtk_gl_area_attach_buffers(gl_area);
      screendump_tga_internal(output_file_name, w, h, sf, screendump_framebuffer.get_fbo());

   } else {
      glViewport(0, 0, w, h);
      // use this, rather than glBindFramebuffer(GL_FRAMEBUFFER, 0); ... just Gtk things.
      gtk_gl_area_attach_buffers(gl_area);
      render_scene_to_base_framebuffer();
   }

   glEnable(GL_DEPTH_TEST);


   return FALSE;
}

void
graphics_info_t::render_scene_to_base_framebuffer() {

   glEnable(GL_DEPTH_TEST);
   shader_for_screen.Use();
   glBindVertexArray(screen_quad_vertex_array_id);

   const glm::vec3 &bg = background_colour;
   glClearColor(bg[0], bg[1], bg[2], 1.0); // this can be seen
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   GLuint pid = shader_for_screen.get_program_id();
   glActiveTexture(GL_TEXTURE0 + 0);
   glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_colour());
   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, screen_framebuffer.get_texture_depth());
   shader_for_screen.set_int_for_uniform("screenTexture", 0);
   shader_for_screen.set_int_for_uniform("screenDepth", 1);
   GLenum err = glGetError(); if (err) std::cout << "render() D err " << err << std::endl;
   shader_for_screen.set_bool_for_uniform("do_ambient_occlusion", shader_do_ambient_occlusion_flag);
   shader_for_screen.set_bool_for_uniform("do_outline", shader_do_outline_flag);

   glDrawArrays(GL_TRIANGLES, 0, 6);
   err = glGetError(); if (err) std::cout << "render() E err " << err << std::endl;

}


void
graphics_info_t::reset_frame_buffers(int width, int height) {

   graphics_info_t g;
   unsigned int sf = framebuffer_scale;
   unsigned int index_offset = 0;
   // std::cout << "debug:: reset_frame_buffers() with sf " << sf << " "
   // << width << " x " << height << std::endl;
   g.screen_framebuffer.init(sf * width, sf * height, index_offset, "screen");
   GLenum err = glGetError(); if (err) std::cout << "reset_frame_buffers() err " << err << std::endl;

   // index_offset = 0;
   // g.blur_framebuffer.init(width, height, index_offset, "blur");
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

void
graphics_info_t::translate_in_screen_z(float step_size) {

   // this looks a bit weird without perspective view

   glm::vec3 ep = get_world_space_eye_position();
   glm::vec3 rc = get_rotation_centre();
   glm::vec3 delta = rc - ep;
   glm::vec3 delta_uv = normalize(delta);

   // more zoomed in has smaller zoom than zoomed out. Zoomed out is ~100. Zoomed in is ~25
   glm::vec3 step = 0.005 * step_size * zoom * delta_uv;

   if (true) // debug
      std::cout << "ep " << glm::to_string(ep) << " rc " << glm::to_string(rc)
                << " zoom " << zoom << " step " << glm::to_string(step) << std::endl;

   add_to_rotation_centre(step);

}

void
graphics_info_t::setup_draw_for_particles() {
   graphics_info_t g;
   if (particles.empty()) {
      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
      glm::vec3 rc = g.get_rotation_centre();
      // std::cout << "making " << n_particles << " around " << glm::to_string(rc) << std::endl;
      particles.make_particles(n_particles, rc);
   }
   // passing user_data and Notify function at the end
   if (! do_tick_particles) {
      do_tick_particles = true;
      int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
      graphics_info_t::idle_function_spin_rock_token = new_tick_id;
   }
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

      do_tick_boids = true;
      int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);

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
      lines_mesh_for_boids_box.setup(&shader_for_lines);
   }
}


void
graphics_info_t::draw_boids() {

   if (boids.size() > 0) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::vec3 eye_position = get_world_space_eye_position();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      glm::vec4 bg_col(background_colour, 1.0);
      mesh_for_boids.draw(&shader_for_instanced_objects,
                          mvp, view_rotation_matrix, lights, eye_position, bg_col,
                          shader_do_depth_fog_flag);

      lines_mesh_for_boids_box.draw(&shader_for_lines, mvp, view_rotation_matrix);
   }
}

//static
void
graphics_info_t::setup_pulse_identification() {

   // we don't set it up here - we set it up in the function that sets up the tick timeout.
   //
   // lines_mesh_for_identification_pulse.setup(&shader_for_lines);
}


void
graphics_info_t::draw_identification_pulse() {

   if (! lines_mesh_for_identification_pulse.empty()) {
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 view_rotation_matrix = get_view_rotation();
      lines_mesh_for_identification_pulse.draw(&shader_for_lines_pulse, mvp, view_rotation_matrix, true);
   }
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
                if (idle_function_spin_rock_token != -1) {
                   std::cout << "Removing the idle function\n";
                   gtk_widget_remove_tick_callback(glareas[0], idle_function_spin_rock_token);
                   idle_function_spin_rock_token = -1;
                } else {

                   do_tick_spin = true;
                   int spin_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);

                   // this is not a good name if we are storing a generic tick function id.
                   idle_function_spin_rock_token = spin_tick_id;
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

   auto l22 = []() {
                  graphics_info_t g;
                  g.setup_draw_for_particles();
                  return gboolean(TRUE);
             };

   // boids
   auto l23 = [] () {
                 graphics_info_t g;
                 if (! graphics_info_t::do_tick_boids)
                    graphics_info_t::do_tick_boids = true;
                 else
                    graphics_info_t::do_tick_boids = false;

                 g.setup_draw_for_boids();

                 if (! graphics_info_t::do_tick_boids)
                    std::cout << "--------- key press ----------- do_tick_boids "
                           << graphics_info_t::do_tick_boids << std::endl;
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
   //   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_e,      key_bindings_t(l20, "EigenFlip Active Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_l,      key_bindings_t(l21, "Label/Unlabel Active Atom")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_q,      key_bindings_t(l22, "Particles")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_b,      key_bindings_t(l23, "Murmuration")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_y,      key_bindings_t(l24, "Add Terminal Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_k,      key_bindings_t(l25, "Fill Partial Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_K,      key_bindings_t(l26, "Delete Sidechain")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l28, "NCS Other Chain")));

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
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].delete_residue(residue_spec);
                    }
                 }
                 return gboolean(TRUE);
              };
   key_bindings_t delete_residue_key_binding(ldr, "Delete Residue");
   std::pair<keyboard_key_t, key_bindings_t> pdel(keyboard_key_t(GDK_KEY_d, true), delete_residue_key_binding);
   kb_vec.push_back(pdel);



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
      int contour_idle_token = g_idle_add(idle_contour_function, glareas[0]);
      std::cout << "INFO:: contour level for map " << imol_scroll << " is "
                << molecules[imol_scroll].contour_level << std::endl;
      set_density_level_string(imol_scroll, molecules[imol_scroll].contour_level);
      display_density_level_this_image = 1;

      graphics_draw(); // queue
   }
}
