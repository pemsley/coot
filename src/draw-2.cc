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

enum {VIEW_CENTRAL_CUBE, ORIGIN_CUBE};

gint idle_contour_function(gpointer data);


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


glm::vec3
get_camera_up_direction(const glm::mat4 &mouse_quat_mat) {

   glm::vec4 z_p(0.0f, 1.0f, 0.0f, 1.0f);
   glm::vec4 r = z_p * mouse_quat_mat;
   glm::vec3 r3(r);
   return r3;
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
      // quaternion (ther order of operations is not yet clear to me).

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


// Make these be part of graphics_info_t

Mesh mesh_for_particles("mesh-for-particles");
int n_particles = 100;
particle_container_t particles;

bool do_tick_particles = false;
bool do_tick_spin = false;

gboolean
glarea_tick_func(GtkWidget *widget,
                 GdkFrameClock *frame_clock,
                 gpointer data) {

   if (do_tick_particles) {
      if (particles.empty()) {
         do_tick_particles = false;
         return FALSE;
      } else {
         particles.update_particles();
         mesh_for_particles.update_instancing_buffer_data_for_particles(particles);
      }
   }

   if (do_tick_spin) {
      float delta = 0.002;
      glm::vec3 EulerAngles(0, delta, 0);
      glm::quat quat_delta(EulerAngles);
      glm::quat normalized_quat_delta(glm::normalize(quat_delta));
      glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
      graphics_info_t::glm_quat = glm::normalize(product);
   }

   gtk_widget_queue_draw(widget); // needed?             

   return TRUE;
}



void
graphics_info_t::draw_map_molecules(bool draw_transparent_maps) {

   // run through this molecule loop twice - for opaque then transparent maps
   // first, a block that decides if we need to do anything.

   bool needs_blend_reset = false;

   //

   bool cosine_dependent_map_opacity = false;

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

   if (cosine_dependent_map_opacity) {
      needs_blend_reset = true;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   if (!draw_transparent_maps || n_transparent_maps > 0) {

      glLineWidth(map_line_width);
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
         if (draw_transparent_maps)
            if (m.is_an_opaque_map())
               continue; // not this round

         if (m.n_vertices_for_map_VertexArray > 0) {

            err = glGetError(); if (err) std::cout << "   draw_map_molecules() --- map start --- error " << std::endl;

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
         glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec4 ep(get_world_space_eye_position(), 1.0);
         glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for eye position " << err << std::endl;

         shader.set_bool_for_uniform("do_depth_fog", graphics_info_t::shader_do_depth_fog_flag);
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
         glm::vec4 bgc(background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         GLuint eye_position_uniform_location = shader.eye_position_uniform_location;
         glm::vec4 ep = glm::vec4(get_world_space_eye_position(), 1.0);
         glUniform4fv(eye_position_uniform_location, 1, glm::value_ptr(ep));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for eye position " << err << std::endl;

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
      unsigned int     *flat_indices = new unsigned int[n_triangles_for_atom_pull_restraints * 3];
      unsigned int     *flat_indices_start = flat_indices;
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
               float bl_stick = bl_pull - arrow_head_length;
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
      err = glGetError(); if (err) std::cout << "GL error bonds 17c\n";
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(0 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17c\n";
      glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(1 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17c\n";
      glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(vertex_with_rotation_translation),
                            reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
      err = glGetError(); if (err) std::cout << "GL error bonds 17c\n";

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

   draw_model_molecules();
   draw_intermediate_atoms();
   draw_atom_pull_restraints();

   draw_meshes(); // get a better name

   draw_map_molecules(false); // transparency

   draw_unit_cells();

   draw_environment_graphics_object();

   draw_generic_objects();

   // transparent things...

   draw_particles();

   draw_map_molecules(true);

}


// This does (draws) symmetry too.
//
// static
void
graphics_info_t::draw_environment_graphics_object() {

#if 0   
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

            bool do_depth_fog = true;
            mesh_for_environment_distances.mesh.draw(&shader_for_moleculestotriangles,
                                                     mvp, view_rotation,
                                                     lights, eye_position, bg_col,
                                                     do_depth_fog);

            if (show_symmetry) {
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

   // At some stage, enable this.  Currently (I think) the winding on the atoms is the wrong
   // way around - replace with octaspheres and octahemispheres.
   // I have just now changed the winding on the "solid" map triangles - and now it looks
   // fine.
   // 
   // glEnable(GL_CULL_FACE); // if I enable this, then I get to see the back side
                              // of the atoms. It's a weird look.

   // Make antialised lines
   if (false) {
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_LINE_SMOOTH);
   }

   g.setup_lights();

   gtk_gl_area_attach_buffers(GTK_GL_AREA(g.glareas[0])); // needed?   
   particles.make_particles(n_particles, g.get_rotation_centre());
   mesh_for_particles.setup_instancing_buffers_for_particles(particles.size());

   err = glGetError();
   if (err) std::cout << "on_glarea_realize() --end-- with err " << err << std::endl;

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

   // is this needed? - does the context ever change?
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

      glBindVertexArray(0); // here is not the place to call this.
   }


   graphics_info_t::blur_framebuffer.bind();
   // glDisable(GL_DEPTH_TEST);
   glEnable(GL_DEPTH_TEST);

   // Screen shader (ambient occlusion)

   {

      bool do_ambient_occlusion = shader_do_ambient_occlusion_flag;
      graphics_info_t::shader_for_screen.Use();
      glBindVertexArray(screen_quad_vertex_array_id);

      // glClearColor(0.5, 0.2, 0.2, 1.0);
      const glm::vec3 &bg = graphics_info_t::background_colour;
      glClearColor(bg[0], bg[1], bg[2], 1.0);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      GLuint pid = graphics_info_t::shader_for_screen.get_program_id();
      glActiveTexture(GL_TEXTURE0 + 1);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::screen_framebuffer.get_texture_colour());
      glUniform1i(glGetUniformLocation(pid, "screenTexture"), 1);
      glActiveTexture(GL_TEXTURE0 + 2);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::screen_framebuffer.get_texture_depth());
      glUniform1i(glGetUniformLocation(pid, "screenDepth"), 2);
      err = glGetError(); if (err) std::cout << "on_glarea_render() D err " << err << std::endl;
      graphics_info_t::shader_for_screen.set_bool_for_uniform("do_ambient_occlusion", do_ambient_occlusion);

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

      shader.set_float_for_uniform("zoom", zoom);
      // glUniform1f(shader.zoom_uniform_location, graphics_info_t::zoom);
      err = glGetError(); if (err) std::cout << "on_glarea_render() blur-A err " << err << std::endl;
      // glUniform1i(shader.is_perspective_projection_uniform_location, perspective_projection_flag);
      shader.set_bool_for_uniform("is_perspective_projection", perspective_projection_flag);
      err = glGetError(); if (err) std::cout << "   on_glarea_render() blur-A2 error " << std::endl;

      GLuint pid = graphics_info_t::shader_for_blur.get_program_id();
      glActiveTexture(GL_TEXTURE0 + 1);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::blur_framebuffer.get_texture_colour());
      glUniform1i(glGetUniformLocation(pid, "screenTexture"), 1); // was 1
      glActiveTexture(GL_TEXTURE0 + 2);
      glBindTexture(GL_TEXTURE_2D, graphics_info_t::blur_framebuffer.get_texture_depth());
      glUniform1i(glGetUniformLocation(pid, "screenDepth"), 2); // was 2
      err = glGetError(); if (err) std::cout << "on_glarea_render() blur-B err " << err << std::endl;
      shader.set_bool_for_uniform("do_depth_blur", graphics_info_t::shader_do_depth_blur_flag);
      shader.set_bool_for_uniform("do_outline", graphics_info_t::shader_do_outline_flag);

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

      // implicit type cast
      gboolean handled = g.check_if_moving_atom_pull(was_a_double_click);

      if (! handled) {
         if (was_a_double_click) {
            pick_info nearest_atom_index_info = g.atom_pick_gtk3(false);
            if (nearest_atom_index_info.success == GL_TRUE) {
               int im = nearest_atom_index_info.imol;
               g.molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
               g.graphics_draw();
            }
         }
      }
   }


   g.check_if_in_range_defines(event, mask);
   return TRUE;
}

gboolean
on_glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

   if (graphics_info_t::in_moving_atoms_drag_atom_mode_flag) {
      graphics_info_t g;
      g.unset_moving_atoms_currently_dragged_atom_index();
      g.do_post_drag_refinement_maybe();
      graphics_info_t::in_moving_atoms_drag_atom_mode_flag = 0;
   }

   if (event->state & GDK_BUTTON2_MASK) {
      graphics_info_t g;
      pick_info nearest_atom_index_info = g.atom_pick_gtk3(false);
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

   // This should be a graphics_info_t function

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

   glm::vec4 delta(worldPos_1 / worldPos_1.w - worldPos_2 / worldPos_2.w);
   glm::vec3 delta_v3(delta);

   g.add_to_rotation_centre(delta_v3);
   g.update_maps();
   if (graphics_info_t::glareas.size() > 0)
      int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);
}

gboolean
on_glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   int r = 0;
   graphics_info_t g;

   // split this function up before it gets too big.

   int x_as_int, y_as_int;
   GdkModifierType mask;
   GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());
   GdkDevice *mouse = gdk_seat_get_pointer(seat);
   gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);

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

         bool handled = false;
         if (g.in_moving_atoms_drag_atom_mode_flag) {
            if (g.last_restraints_size() > 0) {
               // move an already picked atom
               g.move_atom_pull_target_position(x_as_int, y_as_int);
               handled = true;
            } else {
               // don't allow translation drag of the
               // intermediate atoms when they are a rotamer:
               //
               if (! g.rotamer_dialog) {
                  // e.g. translate an added peptide fragment.
                  g.move_moving_atoms_by_simple_translation(x_as_int, y_as_int);
               }
            }
         }
         if (! handled) {
            GtkAllocation allocation;
            gtk_widget_get_allocation(widget, &allocation);
            int w = allocation.width;
            int h = allocation.height;
            graphics_info_t::update_view_quaternion(w, h);
         }
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
      std::cout << "making " << n_particles << " around " << glm::to_string(rc) << std::endl;
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
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_q,      key_bindings_t(l22, "Particles")));

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
        std::cout << "key-binding for key: " << it->first.gdk_key << " : "
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
