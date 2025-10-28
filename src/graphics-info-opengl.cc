/*
 * src/graphics-info-opengl.cc
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

#ifdef USE_PYTHON
#include <Python.h>
#endif

#include <random>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
// #include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>

#define ENABLE_NLS // 20220606-PE fixes weird dcgettext() compiler errors
#include "graphics-info.h"

std::vector<std::reference_wrapper<Shader> > get_shader_refs() {

   std::vector<std::reference_wrapper<Shader> > rs = { graphics_info_t::shader_for_maps,
                                                       graphics_info_t::shader_for_map_caps,
                                                       graphics_info_t::shader_for_models,
                                                       graphics_info_t::shader_for_outline_of_active_residue,
                                                       graphics_info_t::shader_for_model_as_meshes,
                                                       graphics_info_t::shader_for_symmetry_atoms_bond_lines,
                                                       graphics_info_t::shader_for_central_cube,
                                                       graphics_info_t::shader_for_origin_cube,
                                                       graphics_info_t::shader_for_hud_text,
                                                       graphics_info_t::shader_for_hud_geometry_bars,
                                                       graphics_info_t::shader_for_hud_geometry_labels,
                                                       graphics_info_t::shader_for_hud_geometry_tooltip_text,
                                                       graphics_info_t::shader_for_hud_buttons,
                                                       graphics_info_t::shader_for_hud_image_texture,
                                                       graphics_info_t::shader_for_atom_labels,
                                                       graphics_info_t::shader_for_moleculestotriangles,
                                                       graphics_info_t::shader_for_hud_lines,
                                                       graphics_info_t::shader_for_lines,
                                                       graphics_info_t::shader_for_lines_pulse,
                                                       graphics_info_t::shader_for_rama_balls,
                                                       graphics_info_t::shader_for_particles,
                                                       graphics_info_t::shader_for_instanced_objects,
                                                       graphics_info_t::shader_for_extra_distance_restraints,
                                                       graphics_info_t::shader_for_happy_face_residue_markers,
                                                       graphics_info_t::shader_for_rama_plot_phi_phis_markers,
                                                       graphics_info_t::shader_for_rama_plot_axes_and_ticks,
                                                       graphics_info_t::shader_for_ligand_view,
                                                       graphics_info_t::shader_for_texture_meshes,
                                                       graphics_info_t::shader_for_meshes,
                                                       graphics_info_t::shader_for_background_image,
                                                       graphics_info_t::shader_for_tmeshes};

   if (! graphics_info_t::graphics_is_gl_es) {
      rs.push_back(graphics_info_t::shader_for_effects);
      rs.push_back(graphics_info_t::shader_for_dof_blur_by_texture_combination);
      rs.push_back(graphics_info_t::shader_for_meshes_with_shadows);
      rs.push_back(graphics_info_t::shader_for_meshes_shadow_map);
      rs.push_back(graphics_info_t::shader_for_instanced_meshes_shadow_map);
      rs.push_back(graphics_info_t::shader_for_meshes_for_ssao);
      rs.push_back(graphics_info_t::shader_for_instanced_meshes_for_ssao);
      rs.push_back(graphics_info_t::shader_for_instanced_meshes_with_shadows);
      rs.push_back(graphics_info_t::shader_for_tmeshes_for_ssao);
      rs.push_back(graphics_info_t::shader_for_tmeshes_with_shadows);
      rs.push_back(graphics_info_t::shader_for_texture_meshes_shadow_map);
      rs.push_back(graphics_info_t::shader_for_rotation_centre_cross_hairs_for_ssao);
      rs.push_back(graphics_info_t::shader_for_tmeshes_for_ssao);
      rs.push_back(graphics_info_t::shader_for_shadow_map_image_texture_mesh);
      rs.push_back(graphics_info_t::shaderGeometryPass);
      rs.push_back(graphics_info_t::shader_for_happy_face_residue_markers_for_ssao);
      rs.push_back(graphics_info_t::shader_for_x_blur);
      rs.push_back(graphics_info_t::shader_for_y_blur);
      rs.push_back(graphics_info_t::shaderSSAO);
      rs.push_back(graphics_info_t::shaderSSAOBlur);
   }

   return rs;
}


bool
graphics_info_t::init_shader(const std::string &shader_file_name) {

   bool status = false;
   std::vector<std::reference_wrapper<Shader> > shader_refs = get_shader_refs();
   std::vector<std::reference_wrapper<Shader> >::iterator it;
   for (it=shader_refs.begin(); it!=shader_refs.end(); ++it) {
      if (it->get().name == shader_file_name) {
         Shader &shader(it->get());
         std::cout << "init_shader(): found the shader " << shader.name << std::endl;
         shader.init(shader_file_name, Shader::Entity_t::NONE);
         status = true;
      }
   }
   std::cout << "--- done init_shader() ---" << std::endl;
   return status;
}



bool
graphics_info_t::init_shaders() {

   std::vector<std::reference_wrapper<Shader> > shaders = get_shader_refs();

   bool status = true;  // success

   std::string p = coot::package_data_dir();
   std::string d = coot::util::append_dir_dir(p, "shaders");
   if (false)
      std::cout << "INFO:: shader default dir: " << d << std::endl;
   std::vector<std::reference_wrapper<Shader> >::iterator it;
   for (it=shaders.begin(); it!=shaders.end(); ++it)
      it->get().set_default_directory(d);

   shader_for_tmeshes.init("texture-meshes.shader",                                 Shader::Entity_t::MAP);  // Hmm! where is this used? (duplicate)
   shader_for_outline_of_active_residue.init("outline-of-active-residue.shader", Shader::Entity_t::MODEL);
   shader_for_maps.init("map.shader", Shader::Entity_t::MAP);
   shader_for_map_caps.init("draw-map-cap.shader", Shader::Entity_t::MAP);
   shader_for_models.init("model.shader", Shader::Entity_t::MODEL);
   shader_for_central_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   shader_for_origin_cube.init("central-cube.shader", Shader::Entity_t::INFRASTRUCTURE);
   shader_for_hud_text.init("hud-text.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_geometry_bars.init("hud-bars.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_geometry_labels.init("hud-labels.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_image_texture.init("hud-image-texture.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_atom_labels.init("atom-label.shader", Shader::Entity_t::MODEL);
   shader_for_moleculestotriangles.init("moleculestotriangles.shader", Shader::Entity_t::MAP);
   shader_for_lines.init("lines.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_lines_pulse.init("lines-pulse.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_rama_balls.init("rama-balls.shader", Shader::Entity_t::MODEL);
   shader_for_particles.init("particles.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_instanced_objects.init("instanced-objects.shader", Shader::Entity_t::INSTANCED_DISPLAY_OBJECT);
   shader_for_extra_distance_restraints.init("extra-distance-restraints.shader", Shader::Entity_t::INSTANCED_DISPLAY_OBJECT);
   shader_for_hud_geometry_tooltip_text.init("hud-geometry-tooltip-text.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_happy_face_residue_markers.init("residue-markers.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   shader_for_ligand_view.init("ligand-view.shader", Shader::Entity_t::NONE);
   shader_for_model_as_meshes.init("model-as-mesh.shader", Shader::Entity_t::MODEL);
   shader_for_symmetry_atoms_bond_lines.init("symmetry-atoms-lines.shader", Shader::Entity_t::MAP);
   shader_for_hud_buttons.init("hud-bars.shader", Shader::Entity_t::HUD_TEXT); // ! needs a better name, c.f. shader_for_hud_geometry_bars
   shader_for_rama_plot_axes_and_ticks.init("rama-plot-axes-and-ticks.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_rama_plot_phi_phis_markers.init("rama-plot-phi-psi-markers.shader", Shader::Entity_t::HUD_TEXT);
   shader_for_hud_lines.init("hud-lines.shader", Shader::Entity_t::MODEL);
   shader_for_background_image.init("background-image.shader", Shader::Entity_t::NONE);
   shader_for_meshes.init("meshes.shader", Shader::Entity_t::MAP); // 20220208-PE temporay while crow code is merged.
   shader_for_texture_meshes.init("texture-meshes.shader", Shader::Entity_t::MAP);

   if (graphics_is_gl_es) {
   } else {
      // crows
      shader_for_happy_face_residue_markers_for_ssao.init("residue-markers-for-ssao.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
      shader_for_meshes_with_shadows.init("meshes-with-shadows.shader",                Shader::Entity_t::MAP);
      shader_for_meshes_shadow_map.init("meshes-for-shadow-map.shader",                Shader::Entity_t::MAP);
      shader_for_instanced_meshes_shadow_map.init("instanced-meshes-for-shadow-map.shader", Shader::Entity_t::MAP);
      shader_for_meshes_for_ssao.init("meshes-for-ssao.shader",                        Shader::Entity_t::MAP);
      shader_for_instanced_meshes_for_ssao.init("instanced-meshes-for-ssao.shader",    Shader::Entity_t::MAP);
      shader_for_tmeshes_for_ssao.init("texture-meshes-for-ssao.shader",               Shader::Entity_t::MAP);
      shader_for_tmeshes_with_shadows.init("texture-meshes-with-shadows.shader",       Shader::Entity_t::MAP);
      shader_for_texture_meshes_shadow_map.init("texture-meshes-shadow-map.shader",    Shader::Entity_t::MAP);
      shader_for_shadow_map_image_texture_mesh.init("shadow-map-image-texture.shader", Shader::Entity_t::MAP);
      shaderGeometryPass.init("9.ssao_geometry.shader", Shader::Entity_t::NONE);
      shaderSSAO.init(        "9.ssao.shader",          Shader::Entity_t::NONE);
      shaderSSAOBlur.init(    "9.ssao_blur.shader",     Shader::Entity_t::NONE);
      shader_for_instanced_meshes_with_shadows.init("instanced-meshes-with-shadows.shader", Shader::Entity_t::MAP);

      shader_for_effects.init("effects.shader", Shader::Entity_t::NONE);

      // testing image textures
      // camera_facing_quad_shader.init("camera-facing-quad-shader-for-testing.shader", Shader::Entity_t::MODEL);

      // we use the above to make an image/texture in the framebuffer and use then
      // shader_for_screen to convert that framebuffer to the screen buffer.
      shader_for_x_blur.init("blur-x.shader", Shader::Entity_t::SCREEN);
      shader_for_y_blur.init("blur-y.shader", Shader::Entity_t::SCREEN);
      shader_for_dof_blur_by_texture_combination.init("depth-of-field.shader", Shader::Entity_t::SCREEN);

      // long name at the bottom
      shader_for_rotation_centre_cross_hairs_for_ssao.init("rotation-centre-cross-hairs-for-ssao.shader", Shader::Entity_t::NONE);
   }

   for (it=shaders.begin(); it!=shaders.end(); ++it) {
      if (! it->get().get_success_status()) {
         std::cout << "ERROR:: shader \"" <<it->get().name << "\" failed" << std::endl;
         status = false;
      }
   }

   shaders_have_been_compiled = true; // tested in the window resize callback.
   return status;
}


void
graphics_info_t::init_framebuffers(unsigned int width, unsigned int height) { // 20220129-PE a crows thing

   // width and height are passed becasue the window has not been realised yet so it has zero width and height.

   unsigned int index_offset = 0;
   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: init_framebuffers start () err is " << err << std::endl;

   float w = width; // allocation.width;
   float h = height; // allocation.height;

   auto make_generic_framebuffer = [] (const std::string &framebuffer_name,
                                       unsigned int &depthMap_framebuffer,
                                       unsigned int &depthMap_texture,
                                       unsigned int texture_width,
                                       unsigned int texture_height) {

      glGenFramebuffers(1, &depthMap_framebuffer);
      glGenTextures(1, &depthMap_texture);
      glBindTexture(GL_TEXTURE_2D, depthMap_texture);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, texture_width, texture_height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
      float clampColour[] = {1.0f, 1.0f, 1.0f, 1.0f};
      glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, clampColour);

      glBindFramebuffer(GL_FRAMEBUFFER, depthMap_framebuffer);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap_texture, 0);

      // for shadow framebuffer we do this
      glDrawBuffer(GL_NONE);
      glReadBuffer(GL_NONE);

      if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
         std::cout << "Framebuffer for " << framebuffer_name << " not complete!" << std::endl;
      // else
      //    std::cout << "Framebuffer for " << framebuffer_name << " was complete!" << std::endl;
   
      GLenum err = glGetError();
      if (err)
         std::cout << "GL ERROR:: init_framebuffers() post shadow depthmap, error is " << err << std::endl;

      // glBindFramebuffer(GL_FRAMEBUFFER, 0); // standard OpenGL
      // gtk_gl_area_attach_buffers(GTK_GL_AREA(gl_area));  // GTK OpenGL comment out - testing
   };

   make_generic_framebuffer("shadow-depth-framebuffer", shadow_depthMap_framebuffer, shadow_depthMap_texture,
                            shadow_texture_width, shadow_texture_height);

   // ------------------ AO framebuffer ------------------------

   unsigned int attachment_index_colour_texture = 0;
   framebuffer_for_effects.init(w, h, attachment_index_colour_texture, "effects-framebuffer");

   // unsigned int effects_fbo = 0;
   // unsigned int effects_depth_texture = 0;
   // make_generic_framebuffer("effects-framebuffer", effects_fbo, effects_depth_texture, w, h);
   // make_generic_framebuffer("effects-framebuffer", effects_fbo, effects_depth_texture, w, h);

   // ------------------ DOF blur framebuffers ------------------------

   // index_offset is added to GL_COLOR_ATTACHMENT0 in the call too glFramebufferTexture()
   //
   index_offset = 0;
   blur_y_framebuffer.init(w, h, index_offset, "blur-y");
   err = glGetError();
   if (err) std::cout << "GL ERROR:: post blur_y_framebuffer init() err is " << err << std::endl;
   index_offset = 0;
   blur_x_framebuffer.init(w, h, index_offset, "blur-x");
   err = glGetError();
   if (err) std::cout << "GL ERROR::post blur_x_framebuffer init() err is " << err << std::endl;
   index_offset = 0;
   combine_textures_using_depth_framebuffer.init(w, h, index_offset, "new-blur");
   err = glGetError();
   if (err)
      std::cout << "GL ERR:: init_framebuffers() post blur_combine framebuffer init() err is "
                << err << std::endl;

}


glm::vec4
graphics_info_t::unproject(float x, float y, float z) {

   // z is 1 and -1 for front and back (or vice verse).

   if (! glareas[0]) return glm::vec4(0,0,0,0);

   GtkAllocation allocation;
   gtk_widget_get_allocation(glareas[0], &allocation);
   int w = allocation.width;
   int h = allocation.height;

   float mouseX = x / (w * 0.5f) - 1.0f;
   float mouseY = (h - y) / (h * 0.5f) - 1.0f;

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos_f = glm::vec4(mouseX, mouseY, z, 1.0f); // maybe +1
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   if (false) {
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(mvp) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(vp_inv) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(screenPos_f) << std::endl;
      std::cout << "unproject(" << x << "," << y << "," << z << ")   " << glm::to_string(worldPos_f) << std::endl;
   }

   // to turn these into points in the world space, don't forget to divide by w.
   return worldPos_f;

}

// projected_coords are in clip space, so mouse position will have to be converted.
//
// static
glm::vec3
graphics_info_t::unproject_to_world_coordinates(const glm::vec3 &projected_coords) {

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos = glm::vec4(projected_coords, 1.0f);
   glm::vec4 c = vp_inv * screenPos;

   // std::cout << "debug:: unproject_to_world_coordinates() c " << glm::to_string(c) << std::endl;
   double oow = 1.0/c.w;
   return glm::vec3(c.x * oow, c.y * oow, c.z * oow);

}


int
graphics_info_t::blob_under_pointer_to_screen_centre() {

   graphics_info_t g; // needed?
   int r = 0;
   if (use_graphics_interface_flag) {
      int imol_map = Imol_Refinement_Map();
      if (imol_map != -1) {
         if (molecules[imol_map].is_displayed_p()) {
            // OK we have a map to search.
            // coot::Cartesian front = unproject(0.0);
            // coot::Cartesian back  = unproject(1.0);
            // glm::vec4 glm_front = new_unproject(-0.3);
            // glm::vec4 glm_back  = new_unproject( 1.0);

            GtkAllocation allocation = graphics_info_t::get_glarea_allocation();
            int w = allocation.width;
            int h = allocation.height;

            glm::mat4 mvp = graphics_info_t::get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion
            glm::mat4 vp_inv = glm::inverse(mvp);

            // 20220811-PE mouse_current_x and mouse_current_y are set by the motion callback.
            float mouseX_2 = mouse_current_x  / (w * 0.5f) - 1.0f; // mouse_x and mouse_y are updated in the motion callback.
            float mouseY_2 = mouse_current_y  / (h * 0.5f) - 1.0f;
            // I revered the sign here - it does the right thing now.
            glm::vec4 screenPos_1 = glm::vec4(mouseX_2, -mouseY_2, -1.0f, 1.0f);
            glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2,  1.0f, 1.0f);
            glm::vec4 worldPos_1 = vp_inv * screenPos_1;
            glm::vec4 worldPos_2 = vp_inv * screenPos_2;

            double oowp_1 = 1.0/worldPos_1.w;
            double oowp_2 = 1.0/worldPos_2.w;
            coot::Cartesian front(worldPos_1.x * oowp_1, worldPos_1.y * oowp_1, worldPos_1.z * oowp_1);
            coot::Cartesian  back(worldPos_2.x * oowp_2, worldPos_2.y * oowp_2, worldPos_2.z * oowp_2);

            clipper::Coord_orth p1(front.x(), front.y(), front.z());
            clipper::Coord_orth p2( back.x(),  back.y(),  back.z());
            if (true) {
               std::cout << "debug:: blob_under_pointer_to_screen_centre() " << mouse_x << " " << mouse_y << std::endl;
               std::cout << "debug:: blob_under_pointer_to_screen_centre() " << mouseX_2 << " " << mouseY_2 << std::endl;
               std::cout << "debug:: blob_under_pointer_to_screen_centre() " << glm::to_string(screenPos_1) << " "
                         << glm::to_string(screenPos_2) << std::endl;
               std::cout << "debug:: blob_under_pointer_to_screen_centre() " << front << " " << back << std::endl;
               // std::cout << "blob_under_pointer_to_screen_centre() " << p1.format() << " "
               // << p2.format() << std::endl;
            }
            coot::Cartesian rc = g.RotationCentre();

            try {
               clipper::Coord_orth blob =
                  molecules[imol_refinement_map].find_peak_along_line_favour_front(p1, p2);
               coot::Cartesian cc(blob.x(), blob.y(), blob.z());
               // coot::Cartesian cc = front.mid_point(back);
               coot::Cartesian delta = rc - cc;
               // std::cout << "Delta: " << delta << std::endl;
               g.setRotationCentre(cc);
               for(int ii=0; ii<n_molecules(); ii++) {
                  molecules[ii].update_map(auto_recontour_map_flag);
                  molecules[ii].update_symmetry();
               }
               g.make_pointer_distance_objects();
               graphics_draw();
            }
            catch (const std::runtime_error &mess) {
               // 20220202-PE deprecated copy constexpr coot::Cartesian... I wonder what that means.
               std::cout << "debug:: given front " << front << " and back " << back << std::endl;
               std::cout << mess.what() << std::endl;
            }
         }
      } else {
         // 2025-10-01-PE I don't like this popping up when there are no molecules
         // maybe a visual effect would be better - like the red rings.
         if (! molecules.empty()) {
	    std::string s = "WARNING:: Refinement map not selected - no action";
	    std::cout << s << std::endl;
	    info_dialog(s.c_str());
         }
      }
   }
   return r;
}


void
graphics_info_t::set_clipping_front(float v) {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_near_perspective_limit = l * 0.99;
      if (v < screen_z_near_perspective_limit)
         if (v > 2.0)
            screen_z_near_perspective = v;
   } else {
      clipping_front = v;
   }

   // std::cout << "DEBUG:: in set_clipping_front() now planes: front: " << clipping_front
   //           << " back: " << clipping_back
   //           << " eye-position " << glm::to_string(eye_position) << std::endl;

   graphics_draw();
}


void
graphics_info_t::set_clipping_back(float v) {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_far_perspective_limit = l * 1.01;
      if (v > screen_z_far_perspective_limit)
         if (v < 1000.0)
            screen_z_far_perspective = v;
   } else {
      clipping_back = v;
   }
   graphics_draw();
}


void
graphics_info_t::adjust_clipping(float d) {

   if (! perspective_projection_flag) {

      clipping_front = clipping_front * (1.0 + d);
      clipping_back  = clipping_back  * (1.0 + d);

      // std::cout << "now clipping_front" << clipping_front << " clipping_back " << clipping_back << std::endl;

   } else {

      // --- perspective ---

      double l = eye_position.z;
      double zf = screen_z_far_perspective;
      double zn = screen_z_near_perspective;

      // we should (and now do) concern ourselves with the distance to
      // the rotation centre so that the clipping planes are not
      // changed so that rotation centre is clipped.

      if (d < 0) {

         // close down (narrow)

         screen_z_near_perspective = l - (l-zn) * 0.97;
         screen_z_far_perspective  = l + (zf-l) * 0.95;

      } else {

         // expand

         screen_z_far_perspective  = l + (zf-l) * 1.05;
         screen_z_near_perspective = l - (l-zn) * 1.03;

      }

      float screen_z_near_perspective_limit = l * 0.99;
      float screen_z_far_perspective_limit  = l * 1.01;
      if (screen_z_near_perspective > screen_z_near_perspective_limit)
         screen_z_near_perspective = screen_z_near_perspective_limit;
      if (screen_z_far_perspective < screen_z_far_perspective_limit)
         screen_z_far_perspective = screen_z_far_perspective_limit;

      if (screen_z_near_perspective <    2.0) screen_z_near_perspective =    2.0;
      if (screen_z_far_perspective  > 1000.0) screen_z_far_perspective  = 1000.0;

      std::cout << "adjust_clipping(): debug l " << l
                << "    post-manip: " << screen_z_near_perspective << " "
                << screen_z_far_perspective << std::endl;
   }
}

void
graphics_info_t::increase_clipping_front() {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_near_perspective_limit = l * 0.99;
      float v = screen_z_near_perspective * 1.05263;
      if (v < screen_z_near_perspective_limit) {
         if (v > 2.0)
            screen_z_near_perspective = v;
      } else {
         std::cout << "Not moving screen_z_near_perspective to " << v << " eye_position.z " << eye_position.z << std::endl;
      }
   } else {
      float d = 0.05;
      clipping_front = clipping_front * (1.0 - d);
   }
   // std::cout << "in increase_clipping_front() clipping_fron now " << clipping_front << std::endl;

   graphics_draw();
}

void
graphics_info_t::increase_clipping_back() {

   if (perspective_projection_flag) {
      float szfp_start = screen_z_far_perspective;
      screen_z_far_perspective *= 1.02;
      // std::cout << "increase_clipping_back() was " << szfp_start << " now " << screen_z_far_perspective << std::endl;
   } else {
      float d = 0.05;
      clipping_back = clipping_back * (1.0 + d);
   }
   // std::cout << "in increase_clipping_back() clipping_back now " << clipping_back << std::endl;
   graphics_draw();
}

void
graphics_info_t::decrease_clipping_front() {

   if (perspective_projection_flag) {
      double l = eye_position.z;
      float screen_z_near_perspective_limit = l * 0.99;
      float v = screen_z_near_perspective * 0.95;
      if (v < screen_z_near_perspective_limit) {
         if (v > 2.0) {
            screen_z_near_perspective = v;
         }
      } else {
         std::cout << "Not moving screen_z_near_perspective to " << v << " eye_position.z " << eye_position.z << std::endl;
      }
   } else {
      float d = 0.05;
      clipping_front = clipping_front * (1.0 + d);
      // std::cout << "in decrease_clipping_front() clipping_front now " << clipping_front << std::endl;
   }
   graphics_draw();
}

void
graphics_info_t::decrease_clipping_back() {

   if (perspective_projection_flag) {
      float szfp = screen_z_far_perspective;
      float v = screen_z_far_perspective * 0.98;
      if (v > eye_position.z) {
         screen_z_far_perspective = v;
         std::cout << "decrease_clipping_back() was " << szfp << " now " << screen_z_near_perspective << std::endl;
      } else {
         std::cout << "Not moving screen_z_far_perspective_limit " << std::endl;
      }
   } else {
      float d = 0.05;
      clipping_back = clipping_back * (1.0 - d);
   }
   // std::cout << "in decrease_clipping_ack() clipping_back now " << clipping_back << std::endl;
   graphics_draw();
}


void
graphics_info_t::set_view_quaternion(float i, float j, float k, float l) {

   // currently sets glm_quat (that's not a good name)
   // change it it view_quaternion

   // maybe the order will be wrong
   glm::quat q(i, j, k, l);
   view_quaternion = q;
}


void
graphics_info_t::update_view_quaternion(int glarea_width, int glarea_height,
                                        double delta_x_drag, double delta_y_drag) {

   // deltas from when the drag started

   const float &tbs = trackball_size;
   float w = static_cast<float>(glarea_width);
   float h = static_cast<float>(glarea_height);
   bool do_it = true;
   if (mouse_x == 0.0 && mouse_y == 0) {
      do_it = false;
   }
   double current_mouse_x = drag_begin_x + delta_x_drag;
   double current_mouse_y = drag_begin_y + delta_y_drag;
   if (false)
      std::cout << "debug:: update_view_quaterion(): "
                << "new-current-pos: " << std::setw(8) << current_mouse_x << " " << std::setw(8) << current_mouse_y
                << " stored-pos " << std::setw(8) << mouse_x << " " << std::setw(8) << mouse_y
                << " delta " << std::setw(8) << current_mouse_x - mouse_x << " "
                << std::setw(8) << current_mouse_y - mouse_y << std::endl;

   // if (abs(current_mouse_y - mouse_y) > 50) do_it = false;
   // if (abs(current_mouse_x - mouse_x) > 50) do_it = false;

   if (do_it) {
      glm::quat tb_quat = trackball_to_quaternion((2.0 * mouse_x - w)/w, (h - 2.0 * mouse_y)/h,
                                                  (2.0 * current_mouse_x - w)/w,
                                                  (h - 2.0 * current_mouse_y)/h, tbs);
      tb_quat = glm::conjugate(tb_quat);
      auto prod = tb_quat * view_quaternion;
      view_quaternion = glm::normalize(prod);
      // starting down z: view quaternion quat( 0.999986, { 0.003957,  0.002008, 0.002724})
      // rotate above around screen z 180: view quaternion quat(0.000280, {0.003231, 0.006105, 0.999976})
      // rotate above around screen y 180: view quaternion quat(-0.001343, {-0.999951, -0.006849, 0.006969})
      // rotate above aournd screen z 180: view quaternion quat(0.008224, {-0.003389, 0.999938, 0.006631})

      // down x: view quaternion quat(0.499870, {-0.501352, 0.497861, 0.500910})
      // rotate screen z 180: view quaternion quat(-0.502290, {-0.498534, -0.498968, 0.500200})
      // rotate screen y 180: view quaternion quat(0.497855, {0.505183, -0.492305, 0.504545})
      // rotate screen z 180: view quaternion quat(0.503984, {-0.497813, -0.503251, -0.494895})

      // down y: view quaternion quat(0.001501, {0.007544, -0.707126, -0.707045})
      // rotate screen z: view quaternion quat(0.707459, {0.706673, 0.010734, 0.000401})
      // rotate screen y: view quaternion quat(-0.005150, {0.003813, -0.709497, 0.704679})
      // rotate screen z: view quaternion quat(-0.705722, {0.708459, 0.005714, -0.002985})

      if (false)
         std::cout << "view quaternion " << glm::to_string(view_quaternion) << std::endl;
   }
   mouse_x = current_mouse_x;
   mouse_y = current_mouse_y;

}



#include "glarea_tick_function.hh"

// static
// extra_annotation is a default argument, default false;
void
graphics_info_t::setup_cylinder_clashes(const coot::atom_overlaps_dots_container_t &c,
                                        int imol, float tube_radius, bool extra_annotation) {

   auto get_clashes_object_name = [] (int imol) {
      std::string clashes_name = "  Molecule " + coot::util::int_to_string(imol) + ":";
      clashes_name += " clashes insta-mesh";
      return clashes_name;
   };

   // clashes - 20211008-PE this should be in a lambda for clarity
   //
   //             We can't do cylinders with this shader! So make a ball instead.
   // 20210910-PE Let's try to make another separate instancing mesh for clashes.
   //
   graphics_info_t g;
   if (c.clashes.size() == 0) {
      if (false)
         std::cout << "zero clashes" << std::endl;
      std::string clashes_name = get_clashes_object_name(imol);
      int clashes_obj_index = generic_object_index(clashes_name);
      // std::cout << "debug:: setup_cylinder_clashes() here with clashes_obj_index " << clashes_obj_index << std::endl;
      if (clashes_obj_index == -1) {
         clashes_obj_index = g.new_generic_object_number_for_molecule(clashes_name, imol); // make static?
         if (imol == -1) {
	    // std::cout << "........ attaching to intermediate atoms!" << std::endl;
	    g.generic_display_objects[clashes_obj_index].attach_to_intermediate_atoms();
	 }
      } else {
	 std::cout << "clearing clashes..." << std::endl;
         g.generic_display_objects[clashes_obj_index].clear();
         if (imol == -1)
            g.generic_display_objects[clashes_obj_index].attach_to_intermediate_atoms();
      }

   } else {
      std::string clashes_name = get_clashes_object_name(imol);
      int clashes_obj_index = generic_object_index(clashes_name);
      if (clashes_obj_index == -1) {
         clashes_obj_index = g.new_generic_object_number_for_molecule(clashes_name, imol); // make static?
         if (imol == -1)
            g.generic_display_objects[clashes_obj_index].attach_to_intermediate_atoms();
      } else {
	 // std::cout << "clearing (2) clashes..." << std::endl;
         g.generic_display_objects[clashes_obj_index].clear();
         if (imol == -1)
            g.generic_display_objects[clashes_obj_index].attach_to_intermediate_atoms();
      }

      float dimmer = 0.66;
      float z_scale = 0.37;
      float unstubby_cap_factor = 1.1/z_scale; // see below
      if (extra_annotation) {
         dimmer = 0.95;
         unstubby_cap_factor = 2.1;
         z_scale = 0.7;
      }
      coot::colour_holder clash_col = colour_values_from_colour_name("#ff59c9");
      clash_col.scale_intensity(dimmer);
      glm::vec4 clash_col_glm(clash_col.red, clash_col.green, clash_col.blue, 1.0);

      // instancing for capped cylinders
      meshed_generic_display_object &obj = g.generic_display_objects[clashes_obj_index];
      // std::cout << ":::::::::: in setup_cylinder_clashes() obj.get_imol() " << obj.get_imol() << std::endl;
      float line_radius = 0.062f;
      line_radius = tube_radius;
      const unsigned int n_slices = 16;
      std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(0,0,0), glm::vec3(0,0,1));
      obj.add_cylinder(pos_pair, clash_col, line_radius, n_slices, true, true,
                       meshed_generic_display_object::ROUNDED_CAP,
                       meshed_generic_display_object::ROUNDED_CAP, extra_annotation, unstubby_cap_factor); // does obj.mesh.import()
      // now I need to init the buffers of obj.mesh.
      Material material;
      if (extra_annotation) {
         material.shininess = 199.9;
         material.specular_strength = 0.9;
      }
      obj.mesh.setup(material); // calls setup_buffers()
      //
      // now accumulate the instancing matrices (the colours will stay the same)
      std::vector<glm::mat4> mats;
      mats.reserve(c.clashes.size());
      for (unsigned int i=0; i<c.clashes.size(); i++) {
         std::pair<glm::vec3, glm::vec3> pos_pair_clash(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                        glm::vec3(coord_orth_to_glm(c.clashes[i].second)));
         const glm::vec3 &start  = pos_pair_clash.first;
         const glm::vec3 &finish = pos_pair_clash.second;
         glm::vec3 b = finish - start;
         glm::vec3 normalized_bond_orientation(glm::normalize(b));
         glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0));
         glm::vec3 sc(1.1, 1.1, z_scale);
         glm::mat4 unit(1.0);
         glm::mat4 mt_1 = glm::translate(unit, start);
         glm::mat4 mt_2 = mt_1 * ori;
         glm::mat4 mt_3 = glm::scale(mt_2, sc);
         mats.push_back(mt_3);
      }

      auto set_display_generic_object_simple = [] (int object_number, short int istate) {
                                                  if (object_number >=0  && object_number < int(generic_display_objects.size())) {
                                                     generic_display_objects[object_number].mesh.set_draw_this_mesh(istate);
                                                  } else {
                                                     std::cout << "ERROR:: object_number in to_generic_object_add_point: "
                                                               << object_number << std::endl;
                                                  }
                                               };

      std::vector<glm::vec4> cols(c.clashes.size(), clash_col_glm);
      unsigned int n_instances = mats.size();
      obj.mesh.setup_rtsc_instancing(nullptr, mats, cols, n_instances, material); // also does setup_buffers()
      obj.mesh.update_instancing_buffer_data(mats, cols); // is this needed?
      set_display_generic_object_simple(clashes_obj_index, 1);

      // add continuous updating
      if (! tick_function_is_active()) {
         tick_function_id = gtk_widget_add_tick_callback(g.glareas[0], glarea_tick_func, 0, 0); // turn off by turn off FPS monitor
      }
      g.do_tick_constant_draw = true;
   }
};


//static
bool
graphics_info_t::get_exta_annotation_state() {

   bool extra_annotation = false;
   time_t times = time(NULL);
   struct tm result;
   localtime_r(&times, &result);

   if (result.tm_mday == 1)
      if (result.tm_mon == 3)
         if (result.tm_sec%5==0)
            extra_annotation = true;
   if (result.tm_mon == 9)
      if (result.tm_mday > 15)
         if (result.tm_sec%5==0)
            extra_annotation = true;
   return extra_annotation;
}


#include "coot-utils/oct.hh"

glm::vec4 colour_holder_to_glm(const coot::colour_holder &ch); // in meshed-generic-display-objects.cc


bool
graphics_info_t::coot_all_atom_contact_dots_are_begin_displayed_for(int imol) const {

   bool status = false;
   for (unsigned int i=0; i<generic_display_objects.size(); i++) {
      const auto &gdo = generic_display_objects[i];
      if (gdo.imol == imol) {
         const std::string &mesh_name = gdo.mesh.name;
         unsigned int n_instances = generic_display_objects[i].mesh.get_n_instances();
         std::cout << "debug mesh " << i << " has name " << mesh_name
                   << " and " << n_instances << " instances" << std::endl;
         if (mesh_name.find("Contact Dots for Molecule") != std::string::npos) {
            status = true;
            break;
         }
         if (mesh_name.find("insta-mesh") != std::string::npos) {
            status = true;
            break;
         }
      }
   }
   return status;
}


// probably not the right place for this function
//
// This should be called with imol = -1 for intermediate atoms.
void
graphics_info_t::coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol) {

   // 20230521-PE why are the dots part of a real molecule, but the clashes
   // are graphics_info_t generic display objects?

   auto get_generic_object_index = [this] (const std::string &name) {
      int idx_new = -1;
      int size = generic_display_objects.size();
      for (int i=0; i<size; i++) {
         if (generic_display_objects[i].mesh.name == name)
            return i;
      }
      idx_new = this->new_generic_object_number(name);
      return idx_new;
   };

   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &v_in) {
      std::vector<s_generic_vertex> v_out(v_in.size());
      for (unsigned int i=0; i<v_in.size(); i++) {
         v_out[i].pos    = v_in[i].pos;
         v_out[i].normal = v_in[i].normal;
         v_out[i].color  = glm::vec4(0.5f, 0.5f, 0.95f, 1.0f);
      }
      return v_out;
   };

   unsigned int octasphere_subdivisions = 1; // make a member of graphics_info_t with an API
   octasphere_subdivisions = contact_dot_sphere_subdivisions; // 20211129-PE done
   // octasphere_subdivisions = 2;

   bool ignore_waters = true;
   // I don't like this part
   std::map<std::string, coot::colour_holder> colour_map;
   colour_map["blue"      ] = colour_values_from_colour_name("blue");
   colour_map["sky"       ] = colour_values_from_colour_name("sky");
   colour_map["sea"       ] = colour_values_from_colour_name("sea");
   colour_map["greentint" ] = colour_values_from_colour_name("greentint");
   colour_map["darkpurple"] = colour_values_from_colour_name("darkpurple");
   colour_map["green"     ] = colour_values_from_colour_name("green");
   colour_map["orange"    ] = colour_values_from_colour_name("orange");
   colour_map["orangered" ] = colour_values_from_colour_name("orangered");
   colour_map["yellow"    ] = colour_values_from_colour_name("yellow");
   colour_map["yellowtint"] = colour_values_from_colour_name("yellowtint");
   colour_map["red"       ] = colour_values_from_colour_name("red");
   colour_map["#55dd55"   ] = colour_values_from_colour_name("#55dd55");
   colour_map["hotpink"   ] = colour_values_from_colour_name("hotpink");
   colour_map["grey"      ] = colour_values_from_colour_name("grey");
   colour_map["magenta"   ] = colour_values_from_colour_name("magenta");
   colour_map["royalblue" ] = colour_values_from_colour_name("royalblue");

   auto colour_string_to_colour_holder = [colour_map] (const std::string &c) {
      return colour_map.find(c)->second;
   };

   coot::atom_overlaps_dots_container_t c;

#if 0 // why did I want to add contact dots for a (static) molecule when moving atoms were being displayed?
      // Weird.

      // get_moving_atoms_lock(__FUNCTION__);
   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         if (moving_atoms_asc->n_selected_atoms > 0) {
            coot::atom_overlaps_container_t overlaps(mol, graphics_info_t::Geom_p(), ignore_waters, 0.5, 0.25);
            c = overlaps.all_atom_contact_dots(contact_dots_density, true);
         }
      }
   }
   // release_moving_atoms_lock(__FUNCTION__);
#endif

   auto colour_string_to_col4 = [colour_map] (const std::string &colour_string) {
      std::map<std::string, coot::colour_holder>::const_iterator it;
      it = colour_map.find(colour_string);
      if (it == colour_map.end()) {
         return glm::vec4(0.95f, 0.15f, 0.95f, 1.0f);
      } else {
         return colour_holder_to_glm(it->second);
      }
   };

   // more sensible

   coot::atom_overlaps_container_t overlaps(mol, graphics_info_t::Geom_p(), ignore_waters, 0.5, 0.25);
   bool do_vdw_surface = all_atom_contact_dots_do_vdw_surface; // static class variable
   c = overlaps.all_atom_contact_dots(contact_dots_density, do_vdw_surface);

   gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   std::string molecule_name_stub = "Contact Dots for Molecule ";
   molecule_name_stub += coot::util::int_to_string(imol);
   molecule_name_stub += ": ";

   glm::vec4 colour_4(0.5f, 0.5f, 0.5f, 1.0f); // tmp
   std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > oct =
      make_octasphere(octasphere_subdivisions, glm::vec3(0,0,0),
                      1.0f, colour_4, true);

   const std::vector<coot::api::vnc_vertex> &vertices = oct.first;
   const std::vector<g_triangle>           &triangles = oct.second;

   float point_size = 0.08f;
   std::unordered_map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
   for (it=c.dots.begin(); it!=c.dots.end(); ++it) {
      float point_size_inner = point_size;
      const std::string &type = it->first;
      const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
      // std::cout << "dots type " << type << " count " << v.size() << std::endl;
      Material material;
      material.specular_strength = 0.5; // default
      if (type == "vdw-surface") {
         material.specular_strength = 0.1; // dull, reduces zoomed out speckles
         point_size_inner = 0.03;
      }
      std::string mesh_name = molecule_name_stub + type;
      int object_index = get_generic_object_index(mesh_name);
      meshed_generic_display_object &obj = generic_display_objects[object_index];
      obj.imol = imol;
      if (imol == -1)
         obj.attach_to_intermediate_atoms();

      obj.mesh.import(convert_vertices(vertices), triangles);
      obj.mesh.setup(material); // calls setup_buffers()
      std::vector<glm::mat4> mats(v.size());
      glm::vec4 col_glm(0.95f, 0.6f, 0.97f, 1.0f);
      std::vector<glm::vec4> cols(v.size(), col_glm);

      for (unsigned int i=0; i<v.size(); i++) {
         // 20230522-PE this seems a bit wierd, but it puts the balls int the right place and the right size.
         glm::vec3 position(v[i].pos.x()/point_size_inner, v[i].pos.y()/point_size_inner, v[i].pos.z()/point_size_inner);
         glm::mat4 ori(1.0f);
         glm::vec3 sc3(point_size_inner, point_size_inner, point_size_inner);
         // mats[i] = glm::scale(ori, sc3) + glm::translate(position);
         mats[i] = glm::translate(glm::scale(ori, sc3), position);

         // colour
         const std::string &colour_string = v[i].col;
         glm::vec4 col4 = colour_string_to_col4(colour_string);
         cols[i] = col4;
      }

      unsigned int n_instances = mats.size();
      obj.mesh.setup_rtsc_instancing(nullptr, mats, cols, n_instances, material); // also does setup_buffers()
      // these vectors need not be of the same size
      obj.mesh.update_instancing_buffer_data(mats, cols);
      obj.mesh.set_draw_this_mesh(1);

   }


   // and now the spikes (which are not spherical)
   //
   float tube_size = point_size * 1.1f;
   setup_cylinder_clashes(c, imol, tube_size, get_exta_annotation_state());

}

void
graphics_info_t::generate_ssao_kernel_samples() {

   auto lerp = [] (float a, float b, float f) {
                  return a + f * (b - a);
               };

   // generate sample kernel
   // ----------------------
   std::uniform_real_distribution<GLfloat> randomFloats(0.0, 1.0); // generates random floats between 0.0 and 1.0
   std::default_random_engine generator;
   ssaoKernel.clear();

   // std::cout << "debug:: generating " << n_ssao_kernel_samples << " SSAO kernel samples" << std::endl;

   for (unsigned int i = 0; i < n_ssao_kernel_samples; ++i) {
      glm::vec3 sample(randomFloats(generator) * 2.0 - 1.0, randomFloats(generator) * 2.0 - 1.0, randomFloats(generator));
      sample = glm::normalize(sample);
      sample *= randomFloats(generator);
      float scale = static_cast<float>(i) / static_cast<float>(n_ssao_kernel_samples); // was 64
      // scale = static_cast<float>(i)/64.0;
      // scale samples s.t. they're more aligned to center of kernel
      scale = lerp(0.1f, 1.0f, scale * scale);
      sample *= scale;
      ssaoKernel.push_back(sample);
   }

}


void
graphics_info_t::init_joey_ssao_stuff(int w, int h) {

   // need to init this stuff:
   //
   // static unsigned int gBuffer;
   // static unsigned int ssaoFBO;
   // static unsigned int ssaoBlurFBO;
   // static Shader shaderGeometryPass;         done
   // static Shader shaderSSAO;                 done
   // static Shader shaderSSAOBlur;             done
   // static Shader shaderLightingPass;         done
   // static unsigned int gPosition;
   // static unsigned int gNormal;
   // static unsigned int gAlbedo;
   // static unsigned int noiseTexture;
   // static unsigned int ssaoColorBuffer;
   // static unsigned int ssaoColorBufferBlur;

   // static void renderQuad();
   // static void renderCube();
   // static std::vector<glm::vec3> ssaoKernel;
   // // Camera camera(glm::vec3(0.0f, 0.0f, 5.0f));
   // static Camera camera;

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: init_joey_ssao_stuff() --- start --- err is " << err << std::endl;

   // attach_buffers();
   // err = glGetError();
   // if (err) std::cout << "GL ERROR:: init_joey_ssao_stuff() post attach_buffers() err is " << err << std::endl;

   // thw window is not realized yet, so the w and h are passed.
   // GtkAllocation allocation;
   // auto gl_area = glareas[0];
   // gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   // int w = allocation.width;
   // int h = allocation.height;


   {

         // the quad for the screen AO texture
   //
      {
         float quadVertices[] = { // vertex attributes for a quad
                                 // positions   // texCoords
                                 -1.0f,  1.0f,  0.0f, 1.0f,
                                 -1.0f, -1.0f,  0.0f, 0.0f,
                                  1.0f, -1.0f,  1.0f, 0.0f,

                                 -1.0f,  1.0f,  0.0f, 1.0f,
                                  1.0f, -1.0f,  1.0f, 0.0f,
                                  1.0f,  1.0f,  1.0f, 1.0f
         };
         glGenVertexArrays(1, &screen_AO_quad_vertex_array_id);  // Use a HUDTextureMesh when this is working
         glBindVertexArray(screen_AO_quad_vertex_array_id);
         glGenBuffers(1, &screen_AO_quad_VBO);
         glBindBuffer(GL_ARRAY_BUFFER, screen_AO_quad_VBO);
         glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
         glEnableVertexAttribArray(0);
         glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), static_cast<void *>(0));
         glEnableVertexAttribArray(1);
         glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void *>(2 * sizeof(float)));
         err = glGetError();
         if (err) std::cout << "init_screen_quads() err is " << err << std::endl;

         unsigned int quadVBO;  // this separate names and put it into the display_info_t as data members.
         // (ideally). We can get away with it for the moment, if we don't use them explicitly
         // if we bind each of the VAOs at render time.

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



   }

   // std::cout << "debug:: in init_joey_ssao_stuff() window w h " << w << " " << h << std::endl;

   glEnable(GL_DEPTH_TEST);

   // the SSAO shaders are done with the others
   // shaderGeometryPass.init("9.ssao_geometry.shader", Shader::Entity_t::NONE);
   // shaderSSAO.init(        "9.ssao.shader",          Shader::Entity_t::NONE);
   // shaderSSAOBlur.init(    "9.ssao_blur.shader",     Shader::Entity_t::NONE);

   // shader_for_effects is done in regular init_shaders()

   // glGenFramebuffers(1, &gBufferFBO);
   // glBindFramebuffer(GL_FRAMEBUFFER, gBufferFBO);

   // when attachment_index_colour_texture is set to 10 there is an OpenGL error on setting
   // up the framebuffer, but it doesn't seem to have a visual effect.
   unsigned int attachment_index_colour_texture = 0; // CHECKME

   // std::cout << "---------------------------------- init_joey_ssao_stuff() here 1 " << std::endl;
   framebuffer_for_ssao_gbuffer.init(w, h, attachment_index_colour_texture, "SSAO-gBuffer-framebuffer");
   // std::cout << "---------------------------------- init_joey_ssao_stuff() here 2 " << std::endl;

   framebuffer_for_ssao_gbuffer.do_gbuffer_stuff(w, h);

   if (true) {
      // tell OpenGL which color attachments we'll use (of this framebuffer) for rendering
      unsigned int attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
      glDrawBuffers(3, attachments);
      // create and attach depth buffer (renderbuffer)
      // unsigned int rboDepth; now a memeber of display_info_t
      glGenRenderbuffers(1, &rboDepth);
      glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
      glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboDepth);
      // finally check if framebuffer is complete
      if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
         std::cout << "Framebuffer for SSAO GBuffer not complete!" << w << " " << h << std::endl;
      // else
      //    std::cout << "Framebuffer for SSAO GBuffer was complete!" << w << " " << h << std::endl;

      // you are not allowed to bind this in opengl.
      // glBindFramebuffer(GL_FRAMEBUFFER, 0);
   }

   // also create framebuffer to hold SSAO processing stage
   // -----------------------------------------------------

   glGenFramebuffers(1, &ssaoFBO);
   glGenFramebuffers(1, &ssaoBlurFBO);
   glBindFramebuffer(GL_FRAMEBUFFER, ssaoFBO);

   // SSAO color buffer
   glGenTextures(1, &ssaoColorBuffer);
   glBindTexture(GL_TEXTURE_2D, ssaoColorBuffer);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_FLOAT, NULL);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaoColorBuffer, 0);
   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout << "SSAO Framebuffer not complete! " << w << " " << h << std::endl;
   // else
   //    std::cout << "SSAO Framebuffer was complete!" << w << " " << h << std::endl;

   // and blur stage
   glBindFramebuffer(GL_FRAMEBUFFER, ssaoBlurFBO);
   glGenTextures(1, &ssaoColorBufferBlur);
   glBindTexture(GL_TEXTURE_2D, ssaoColorBufferBlur);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_FLOAT, NULL);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ssaoColorBufferBlur, 0);
   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout << "SSAO Blur Framebuffer not complete!" << w << " " << h << std::endl;
   // else
   //    std::cout << "SSAO Blur Framebuffer was complete!" << w << " " << h << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, 0);

   generate_ssao_kernel_samples();

   // generate noise texture
   // ----------------------
   std::uniform_real_distribution<GLfloat> randomFloats(0.0, 1.0); // generates random floats between 0.0 and 1.0
   std::default_random_engine generator;
   std::vector<glm::vec3> ssaoNoise;
   for (unsigned int i = 0; i < 16; i++) {
      glm::vec3 noise(randomFloats(generator) * 2.0 - 1.0, randomFloats(generator) * 2.0 - 1.0, 0.0f); // rotate around z-axis (in tangent space)
      ssaoNoise.push_back(noise);
   }

   glGenTextures(1, &noiseTexture);
   glBindTexture(GL_TEXTURE_2D, noiseTexture);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 4, 4, 0, GL_RGB, GL_FLOAT, &ssaoNoise[0]);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   // std::cout << "DEBUG:: in init_joey_ssao_stuff() gPosition gNormal gAlbedo noiseTexture "
   // << gPosition << " " << gNormal << " " << gAlbedo << " " << noiseTexture << std::endl;

   err = glGetError();
   if (err)
      std::cout << "ERROR init_joey_ssao_stuff() end err is " << err << std::endl;

}

void
graphics_info_t::read_test_gltf_models() {

   read_some_test_models();
}


void
graphics_info_t::read_some_test_models() {

   auto setup_tmesh_for_floor = [] () {
                                   std::vector<TextureMeshVertex> vertices;
                                   std::vector<g_triangle> triangles;
                                   glm::vec4 c(0.5, 0.5, 0.5, 1.0); // colour
                                   // top
                                   TextureMeshVertex v0(glm::vec3(0, 0, 0), glm::vec3(0,0,1), c, glm::vec2(0,0)); vertices.push_back(v0);
                                   TextureMeshVertex v1(glm::vec3(1, 0, 0), glm::vec3(0,0,1), c, glm::vec2(3,0)); vertices.push_back(v1);
                                   TextureMeshVertex v2(glm::vec3(0, 1, 0), glm::vec3(0,0,1), c, glm::vec2(0,3)); vertices.push_back(v2);
                                   TextureMeshVertex v3(glm::vec3(1, 1, 0), glm::vec3(0,0,1), c, glm::vec2(3,3)); vertices.push_back(v3);
                                   for (auto &vert : vertices)
                                      vert.position -= glm::vec3(0.5, 0.5, 0.0);
                                   for (auto &vert : vertices)
                                      vert.position *= 34.0f;
                                   triangles.push_back(g_triangle(0,1,2));
                                   triangles.push_back(g_triangle(1,2,3));
                                   TextureMesh tmesh("floor-tmesh");
                                   tmesh.import(vertices, triangles);
                                   tmesh.setup_buffers();
                                   Texture texture_a("wooden-floor.png", Texture::DIFFUSE);
                                   Texture texture_n("wooden-floor-normal-mild.png", Texture::NORMAL);
                                   tmesh.add_texture(TextureInfoType(texture_a, "wooden-floor-diffuse"));
                                   // tmesh.add_texture(TextureInfoType(texture_n, "wooden-floor-normal"));
                                   Model model;
                                   model.add_tmesh(tmesh);
                                   add_model(model);
                                };

   auto setup_tmesh_for_wall = [] () {
                                   std::vector<TextureMeshVertex> vertices;
                                   std::vector<g_triangle> triangles;
                                   glm::vec4 c(0.5, 0.5, 0.5, 1.0); // colour
                                   // top
                                   TextureMeshVertex v0(glm::vec3(0, 0, 0), glm::vec3(0,0,1), c, glm::vec2(0,0)); vertices.push_back(v0);
                                   TextureMeshVertex v1(glm::vec3(1, 0, 0), glm::vec3(0,0,1), c, glm::vec2(3,0)); vertices.push_back(v1);
                                   TextureMeshVertex v2(glm::vec3(0, 1, 0), glm::vec3(0,0,1), c, glm::vec2(0,3)); vertices.push_back(v2);
                                   TextureMeshVertex v3(glm::vec3(1, 1, 0), glm::vec3(0,0,1), c, glm::vec2(3,3)); vertices.push_back(v3);
                                   float angle = M_PI * 0.5f;
                                   for (auto &vert : vertices)
                                      vert.position = glm::rotate(vert.position, angle, glm::vec3(1,0,0));
                                   for (auto &vert : vertices)
                                      vert.position *= 33.0f; // matches floor length
                                   auto save_vertices = vertices;
                                   for (auto &vert : vertices)
                                      vert.position -= glm::vec3(16.0f, -17.0f, 0.0f);
                                   triangles.push_back(g_triangle(0,1,2));
                                   triangles.push_back(g_triangle(1,2,3));
                                   TextureMesh tmesh_1("wall-1-tmesh");
                                   // TextureMesh tmesh_3("wall-3-tmesh");
                                   tmesh_1.import(vertices, triangles);
                                   tmesh_1.setup_buffers();
                                   //Texture texture("brick-wall-gareth-david-KREQ7J7nsNw-unsplash.jpg", Texture::DIFFUSE);
                                   // ORM is
                                   // red:   ambient occlusion
                                   // green: roughness
                                   // blue:  metalicity
                                   bool reversed_normals = true;
                                   Texture texture_a("wall-6_albedo.png", Texture::DIFFUSE);
                                   Texture texture_n("wall-6_normal.png", Texture::NORMAL, reversed_normals);
                                   Texture texture_o("wall-6_orm.png",    Texture::AMBIENT_OCCLUSION_ROUGHNESS_METALICITY);
                                   tmesh_1.add_texture(TextureInfoType(texture_a, "wall-diffuse"));
                                   tmesh_1.add_texture(TextureInfoType(texture_n, "wall-normal"));
                                   tmesh_1.add_texture(TextureInfoType(texture_o, "wall-ambient-occlusion"));
                                   Model model;
                                   model.add_tmesh(tmesh_1);
                                   add_model(model);
                                };

   auto setup_tmesh_for_grass = [] () {
                                   std::vector<TextureMeshVertex> vertices;
                                   std::vector<g_triangle> triangles;
                                   glm::vec4 c(0.5, 0.5, 0.5, 1.0); // colour
                                   // top
                                   TextureMeshVertex v0(glm::vec3(0, 0, 0), glm::vec3(0,0,1), c, glm::vec2(0,0)); vertices.push_back(v0);
                                   TextureMeshVertex v1(glm::vec3(1, 0, 0), glm::vec3(0,0,1), c, glm::vec2(3,0)); vertices.push_back(v1);
                                   TextureMeshVertex v2(glm::vec3(0, 1, 0), glm::vec3(0,0,1), c, glm::vec2(0,3)); vertices.push_back(v2);
                                   TextureMeshVertex v3(glm::vec3(1, 1, 0), glm::vec3(0,0,1), c, glm::vec2(3,3)); vertices.push_back(v3);

                                   for (auto &vert : vertices)
                                      vert.position -= glm::vec3(0.5, 0.5, 0.0);
                                   for (auto &vert : vertices)
                                      vert.position *= 35.0f;

                                   triangles.push_back(g_triangle(0,1,2));
                                   triangles.push_back(g_triangle(1,2,3));
                                   TextureMesh tmesh_1("grass-tmesh");
                                   tmesh_1.import(vertices, triangles);
                                   tmesh_1.setup_buffers();
                                   //Texture texture("brick-wall-gareth-david-KREQ7J7nsNw-unsplash.jpg", Texture::DIFFUSE);
                                   // ORM is
                                   // red:   ambient occlusion
                                   // green: roughness
                                   // blue:  metalicity
                                   // (but looking at the orm images, I find this hard to believe)
                                   //
                                   Texture texture_a("grass_albedo.png", Texture::DIFFUSE);
                                   Texture texture_n("grass_normal.png", Texture::NORMAL);
                                   Texture texture_o("grass_orm.png",    Texture::AMBIENT_OCCLUSION_ROUGHNESS_METALICITY);
                                   tmesh_1.add_texture(TextureInfoType(texture_a, "grass-diffuse"));
                                   tmesh_1.add_texture(TextureInfoType(texture_n, "grass-normal"));
                                   tmesh_1.add_texture(TextureInfoType(texture_o, "grass-orm"));
                                   Model model;
                                   model.add_tmesh(tmesh_1);
                                   add_model(model);
                                };

   attach_buffers();

   std::map<std::string, bool> show_mesh;
   show_mesh["grass"] = true;
   show_mesh["wall"]  = true;
   show_mesh["boxes"] = true;
   show_mesh["crow"]  = true;
   show_mesh["tim"]   = true;
   show_mesh["spike"] = false;
   show_mesh["ribo"]  = false;
   show_mesh["little_chestnut"] = false;
   show_mesh["blacksmith"] = true;
   show_mesh["dwarf_blacksmith"] = true;
   show_mesh["teapots"] = false;
   show_mesh["sponza"] = false;
   show_mesh["port"] = false;
   show_mesh["vila"] = false;

   if (false) { // just one
      show_mesh["grass"] = false;
      show_mesh["wall"]  = false;
      show_mesh["crow"]  = false;
      show_mesh["tim"]   = false;
   }

   if (show_mesh["grass"])
      setup_tmesh_for_grass();

   if (show_mesh["wall"])
      setup_tmesh_for_wall();

   // --- Crow ---

   if (show_mesh["crow"]) {
      TextureMesh crow_tmesh("crow");
      // crow_tmesh.load_from_glTF("crow-17-with-grey-surface.glb");
      crow_tmesh.load_from_glTF("crow-21.glb");
      // crow_tmesh.load_from_glTF("spike-protein-with-ace2-light-green-v8.glb");

      Model crow_model;
      crow_model.add_tmesh(crow_tmesh);
      add_model(crow_model);
   }


   // --- 1TIM Protein ---

   if (show_mesh["tim"]) {
      Mesh tim_mesh("tim");
      tim_mesh.load_from_glTF("1tim-A.glb");
      // tim_mesh.debug_to_file();
      Material mat;
      mat.shininess = 512.0;
      mat.specular_strength = 1.0;
      mat.ambient  = glm::vec4(0.9, 0.9, 0.9, 1.0);
      mat.diffuse  = glm::vec4(0.9, 0.9, 0.9, 1.0);
      mat.turn_specularity_on(true);
      tim_mesh.set_material(mat); // override the material extracted from the gltf
      Model tim_model;
      tim_model.add_mesh(tim_mesh);
      tim_model.scale(0.05f);
      tim_model.translate(glm::vec3(0, 0, 3));
      add_model(tim_model);
   }

   // --- Blacksmith ---

   if (show_mesh["blacksmith"]) {
      TextureMesh blacksmith_tmesh("blacksmith");
      blacksmith_tmesh.load_from_glTF("blacksmith.glb", false);
      Model blacksmith_model;
      blacksmith_model.add_tmesh(blacksmith_tmesh);
      blacksmith_model.translate(glm::vec3(-13.0f, 0.0f, 0.0f));
      add_model(blacksmith_model);
   }

   if (show_mesh["dwarf_blacksmith"]) {
      TextureMesh dwarf_blacksmith_tmesh("blacksmith");
      dwarf_blacksmith_tmesh.load_from_glTF("dwarf_blacksmith.glb", false);
      Model dwarf_blacksmith_model;
      dwarf_blacksmith_model.add_tmesh(dwarf_blacksmith_tmesh);
      dwarf_blacksmith_model.translate(glm::vec3(-13.0f, 0.0f, 0.0f));
      add_model(dwarf_blacksmith_model);
   }
}

int
graphics_info_t::load_gltf_model(const std::string &gltf_file_name) {

   attach_buffers();

   TextureMesh tm("some name"); // extract/replace this from the gltf data
   tm.load_from_glTF(gltf_file_name);
   // e.invert_normals(); // it is shiny on the inside either way around - hmm.

   // why do this?
   if (false) {
      Material mat;
      mat.shininess = 64.0;
      mat.specular_strength = 1.0;
      mat.ambient  = glm::vec4(0.7, 0.7, 0.7, 1.0);
      mat.diffuse  = glm::vec4(0.7, 0.7, 0.7, 1.0);
      mat.turn_specularity_on(true);
      // tm.set_material(mat); // override the material extracted from the gltf
   }
   Model e_model;
   e_model.add_tmesh(tm);
   add_model(e_model);

   // add continuous updating
   if (! tick_function_is_active()) {
      tick_function_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
   }
   do_tick_constant_draw = true;

   return models.size() - 1;

}


void
graphics_info_t::resize_framebuffers_textures_renderbuffers(int width, int height) {

   // std::cout << "DEBUG:: resize_framebuffers_textures_renderbuffers() " << width << " " << height << std::endl;

   framebuffer_for_effects.reset(width, height);
   blur_x_framebuffer.reset(width, height);
   blur_y_framebuffer.reset(width, height);
   combine_textures_using_depth_framebuffer.reset(width, height);

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


// static
Mesh &
graphics_info_t::get_mesh_for_eyelashes() {

   if (mesh_for_eyelashes.empty()) {
      std::string fn = "grey-eyelashes-many-lashes.glb";
      mesh_for_eyelashes.load_from_glTF(fn);
   }
   return mesh_for_eyelashes;
}
