/*
 * src/graphics-info-draw-models.cc
 *
 * Copyright 2022 by Medical Research Council
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

#include "graphics-info.h"
#include "stereo-eye.hh"

void
graphics_info_t::draw_models_for_ssao() {

   Shader *shader_for_tmeshes_p = &shader_for_tmeshes_for_ssao;
   Shader *shader_for_meshes_p  = &shader_for_meshes_for_ssao;
   bool do_orthographic_projection = ! perspective_projection_flag;
   GtkAllocation allocation;
   auto gl_area = glareas[0];
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   graphics_info_t g;
   glm::mat4 model_mat      =      g.get_model_matrix();
   glm::mat4 view_mat       =       g.get_view_matrix();
   glm::mat4 projection_mat = g.get_projection_matrix(do_orthographic_projection, w, h);

   for (unsigned int ii=0; ii<models.size(); ii++) {
      auto &model = models[ii];
      model.draw_for_ssao(shader_for_tmeshes_p, shader_for_meshes_p, model_mat, view_mat, projection_mat);
   }
}

void
graphics_info_t::draw_molecules_for_ssao() {

   bool do_orthographic_projection = ! perspective_projection_flag;
   GtkAllocation allocation;
   auto gl_area = glareas[0];
   gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   glm::mat4 model_mat      =      get_model_matrix();
   glm::mat4 view_mat       =       get_view_matrix();
   glm::mat4 projection_mat = get_projection_matrix(do_orthographic_projection, w, h);

   glDisable(GL_BLEND);
   for (int i=0; i<n_molecules(); i++) {
      if (is_valid_model_molecule(i)) {
         if (molecules[i].draw_it)
            if (! molecules[i].draw_model_molecule_as_lines) {
               // molecules[i].molecule_as_mesh.draw_for_ssao(&shader_for_meshes_for_ssao, model_mat, view_mat, projection_mat);
               molecules[i].draw_molecule_as_meshes_for_ssao(&shader_for_meshes_for_ssao,
                                                             &shader_for_instanced_meshes_for_ssao,
                                                             model_mat, view_mat, projection_mat);
            }
      }
      if (is_valid_map_molecule(i)) {
         if (molecules[i].draw_it_for_map) {
            molecules[i].draw_map_molecule_for_ssao(&shader_for_meshes_for_ssao, model_mat, view_mat, projection_mat);
         }
      }
   }

   glEnable(GL_DEPTH_TEST);
   draw_intermediate_atoms(PASS_TYPE_SSAO);
   draw_intermediate_atoms_rama_balls(PASS_TYPE_SSAO); // currently rama balls are part of intermediate atoms
   draw_rotation_centre_crosshairs(GTK_GL_AREA(gl_area), PASS_TYPE_SSAO);
   draw_bad_nbc_atom_pair_markers(PASS_TYPE_SSAO);

#if 0 // Nice things to have, but they need to work with shader_for_meshes_for_ssao.
   draw_atom_pull_restraints(PASS_TYPE_SSAO);
   draw_instanced_meshes(PASS_TYPE_SSAO);
   draw_environment_graphics_object(PASS_TYPE_SSAO);
   draw_generic_objects(PASS_TYPE_SSAO);
#endif

   draw_meshed_generic_display_object_meshes(PASS_TYPE_SSAO);
   draw_molecules_other_meshes(PASS_TYPE_SSAO); // draws the meshes in a molecules std::vector<Messh> meshes;
   draw_generic_objects(PASS_TYPE_SSAO);
   glDisable(GL_BLEND);

}

// Models with a captial M.
void
graphics_info_t::draw_models(Shader *shader_for_tmeshes_p,
                            Shader *shader_for_meshes_p,
                            Shader *shader_for_tmeshes_with_shadows_p,
                            Shader *shader_for_meshes_with_shadows_p,
                            int graphics_x_size,
                            int graphics_y_size,
                            bool draw_shadows, // draw_shadows default false
                            float shadow_strength, // default 0.4
                            bool show_just_shadows) {

   // 20211118-PE shader_for_meshes_with_shadows_p is not used yet.

   // std::cout << "------ draw_models() with draw_shadows " << draw_shadows << std::endl;

   auto get_Model_rotation_and_translation = [] (const Model &model, unsigned int model_index) {
      glm::mat4 mr(1.0f);
      glm::mat4 mt(1.0f);
      if (model.do_animation) {
         float time_ms = model.duration();
         float time_s = time_ms * 0.001f;
         float sf = 0.3f;
         float fidx = static_cast<float>(model_index);
         float x = -5.0 * sinf(2.0f * fidx + sf * time_s);
         float y =  5.0 * cosf(2.0f * fidx + sf * time_s);
         glm::vec3 t(x,y,0.0f);
         mr = glm::translate(mr, t);
         float theta = sf * time_s + 2.0f * fidx + 0.05; // head toward to inside of the circle a small bit
         glm::vec3 up(0.0f, 0.0f, 1.0f);
         mr = glm::rotate(mr, theta, up);
      }
      return std::make_pair(mr, mt);
   };

   auto get_mvp = [get_Model_rotation_and_translation] (const Model &model, unsigned int model_index, int w, int h) {
      bool do_orthographic_projection = ! perspective_projection_flag; // weird
      glm::mat4 view_matrix = get_view_matrix();

      glm::mat4 m(1.0f);
      glm::vec3 rc = get_rotation_centre();
      rc = glm::vec3(0,0,0);
      // we have to do the translations before the view rotation (not view_matrix - that's a different thing)
      // c.f. get_model_matrix()
      glm::mat4 translation_matrix = glm::translate(m, rc);
      std::pair<glm::mat4, glm::mat4> p = get_Model_rotation_and_translation(model, model_index);
      glm::mat4 Model_translation_matrix = p.second;
      glm::mat4 Model_rotation_matrix = p.first;
      glm::mat4 rotation_matrix = get_model_rotation();
      glm::mat4 model_at_screen_centre = rotation_matrix * Model_rotation_matrix * Model_translation_matrix * translation_matrix;

      glm::mat4  proj_matrix = get_projection_matrix(do_orthographic_projection, w, h);
      glm::mat4 mvp = proj_matrix * view_matrix * model_at_screen_centre;
      return mvp;
   };

   glm::vec3 bg_col = background_colour;
   bool do_depth_fog = shader_do_depth_fog_flag;

   glm::mat4 model_rotation = get_model_rotation();
   glm::vec4 bg_col_v4(bg_col, 1.0f);
   auto ccrc = RotationCentre();
   glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());

   if (draw_shadows) {
      // done elsewhere
   } else {
      for (unsigned int ii=0; ii<models.size(); ii++) {
         auto &model = models[ii];
         glm::mat4 mvp = get_mvp(model, ii, graphics_x_size, graphics_y_size);
	 if (shader_for_tmeshes_p) {
            if (false)
               std::cout << "DEBUG:: draw_models(): drawing texture mesh model " << ii << " "
                         << shader_for_tmeshes_p->name << std::endl;
	    shader_for_tmeshes_p->Use();
	    model.draw_tmeshes(shader_for_tmeshes_p, mvp, model_rotation, lights, eye_position,
                               bg_col_v4, do_depth_fog);
         }

         // now the coloured vertices mesh (for molecule things, not textured things)
         if (shader_for_meshes_p) {
            float opacity = 1.0f;
            model.draw_meshes(shader_for_meshes_p, mvp, model_rotation, lights, eye_position, rc, opacity, bg_col_v4, do_depth_fog);
	 }
      }
   }
}


// Models with a capital M.
void
graphics_info_t::draw_models_with_shadows(Shader *shader_for_tmeshes_with_shadows_p,
					 Shader *shader_for_meshes_with_shadows_p,
					 int graphics_x_size,
					 int graphics_y_size,
					 bool draw_shadows, // maybe pass this flag to the model drawer
					 float shadow_strength,
					 bool show_just_shadows) {

   // crow convertion vars
   auto get_mvp = [] (int x, int y) {
                     return get_molecule_mvp(); // no args passed.... Hmm.
                   };


   glm::vec3 bg_col = background_colour;
   bool do_depth_fog = shader_do_depth_fog_flag;

   glm::mat4 mvp = get_mvp(graphics_x_size, graphics_y_size);
   // std::cout << "debug:: mvp draw_models_with_shadows() is " << glm::to_string(mvp) << std::endl;
   glm::mat4 model_rotation = get_model_rotation();
   glm::vec4 bg_col_v4(bg_col, 1.0f);
   for (unsigned int ii=0; ii<models.size(); ii++) {
      auto &model = models[ii];

      int light_index = 0; // 20211121-PE for now
      glm::mat4 light_view_mvp = get_light_space_mvp(light_index);
      float opacity = 1.0;

      model.draw_with_shadows(shader_for_tmeshes_with_shadows_p,
			      shader_for_meshes_with_shadows_p,
			      mvp, model_rotation, lights, eye_position, opacity,
			      bg_col_v4, do_depth_fog, light_view_mvp,
                              shadow_depthMap_texture, shadow_strength, shadow_softness,
			      show_just_shadows);
   }
}



// Models as in the class Model
void
graphics_info_t::draw_Models_for_shadow_map(unsigned int light_index) {

   // we need to adjust the view_rotation and hence the mvp so that they match that light's point of view

   // non-textured meshes have been removed for the moment.

   std::map<unsigned int, lights_info_t>::const_iterator it;
   it = lights.find(light_index);
   if (it != lights.end()) {
      const auto &light = it->second;
      glm::mat4 mvp_orthogonal = get_mvp_for_shadow_map(light.direction);
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col_v4(background_colour, 1.0f);
      auto ccrc = RotationCentre();
      glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());

      for (unsigned int ii=0; ii<models.size(); ii++) {
         auto &m = models[ii];
         shader_for_texture_meshes_shadow_map.Use(); // needed?


	 // I don't like iterating through tmeshes and meshes here in graphics_info_t

         for (unsigned int j=0; j<m.tmeshes.size(); j++) {
            std::map<unsigned int, lights_info_t> dummy_lights;
            glm::vec3 dummy_eye_position;
            m.draw_tmesh(j, &shader_for_texture_meshes_shadow_map,
			 mvp_orthogonal, model_rotation, dummy_lights, dummy_eye_position,
			 bg_col_v4, false);
         }

	 for (unsigned int j=0; j<m.meshes.size(); j++) {
            std::map<unsigned int, lights_info_t> dummy_lights;
            glm::vec3 dummy_eye_position;
            float opacity = 1.0f;
            m.draw_mesh(j, &shader_for_meshes_shadow_map, mvp_orthogonal, model_rotation,
			dummy_lights, dummy_eye_position, rc, opacity, bg_col_v4, false, false);
	 }
      }
   }
}

void
graphics_info_t::draw_molecules_for_shadow_map(unsigned int light_index) {

   stereo_eye_t eye = stereo_eye_t::MONO; // pass this

   // std::cout << "draw_molecules_for_shadow_map() " << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_molecules_for_shadow_map() -- start -- " << err << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   it = lights.find(light_index);
   if (it != lights.end()) {
      const auto &light = it->second;
      glm::mat4 mvp_orthogonal = get_mvp_for_shadow_map(light.direction);
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col_v4(background_colour, 1.0f);
      auto ccrc = RotationCentre();
      glm::vec3 rc(ccrc.x(), ccrc.y(), ccrc.z());

      // --------- maps ----------

      int n_mols = n_molecules();
      for (int i=0; i<n_mols; i++) {
         if (is_valid_map_molecule(i)) {
            molecule_class_info_t &m = molecules[i];
            if (m.draw_it_for_map) {
               glm::vec3 dummy_eye_position;
               bool gl_lines_mode = false;
               bool show_just_shadows = false;
               bool opacity = 1.0;
               bool do_depth_fog = false;
               // draw the surface - not the mesh - not sure if this is a good idea yet.
               m.map_as_mesh.draw(&shader_for_meshes_shadow_map,
                                  mvp_orthogonal, model_rotation, lights, dummy_eye_position, rc,
                                  opacity, bg_col_v4, gl_lines_mode, do_depth_fog, show_just_shadows);
            }
         }
      }

      // --------- models ----------

      for (int i=0; i<n_mols; i++) {
         if (is_valid_model_molecule(i)) {
            molecule_class_info_t &m = molecules[i];
            if (m.draw_it) {
               glm::vec3 dummy_eye_position;
               bool gl_lines_mode = false;
               bool show_just_shadows = false;
               bool opacity = 1.0;
               bool do_depth_fog = false;
               // model molecule, that is of course.
#if 0 // 20230818-PE pre-model_molecule_meshes
               m.molecule_as_mesh.draw(&shader_for_meshes_shadow_map,
                                       mvp_orthogonal, model_rotation, lights, dummy_eye_position,
                                       opacity, bg_col_v4, gl_lines_mode, do_depth_fog, show_just_shadows);
#endif

               m.model_molecule_meshes.draw(&shader_for_meshes_shadow_map,
                                            &shader_for_instanced_meshes_shadow_map,
                                            eye, mvp_orthogonal, model_rotation, lights, dummy_eye_position,
                                            opacity, bg_col_v4, gl_lines_mode, do_depth_fog, show_just_shadows);
            }
         }
      }
   }

   // ribbons

   draw_meshed_generic_display_object_meshes(PASS_TYPE_GEN_SHADOW_MAP); // draws generic objects

   draw_molecules_other_meshes(PASS_TYPE_GEN_SHADOW_MAP); // draws the meshes in a molecules std::vector<Messh> meshes;

   // more stuff needs to be added here - see draw_molecules_for_ssao() (where they haven't been added either)
}
