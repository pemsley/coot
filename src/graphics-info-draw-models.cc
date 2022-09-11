
#include "graphics-info.h"

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
            if (! molecules[i].draw_model_molecule_as_lines)
               molecules[i].molecule_as_mesh.draw_for_ssao(&shader_for_meshes_for_ssao, model_mat, view_mat, projection_mat);
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

   draw_meshed_generic_display_object_meshes(PASS_TYPE_SSAO); // draws the meshes in a molecules std::vector<Messh> meshes;
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

   // crow convertion vars
   auto get_mvp = [] (int x, int y) {
                     return get_molecule_mvp(); // no args passed.... Hmm.
                   };
   glm::vec3 bg_col = background_colour;
   bool do_depth_fog = shader_do_depth_fog_flag;

   glm::mat4 mvp = get_mvp(graphics_x_size, graphics_y_size);
   glm::mat4 model_rotation = get_model_rotation();
   glm::vec4 bg_col_v4(bg_col, 1.0f);
   if (draw_shadows) {
      // done elsewhere
   } else {
      for (unsigned int ii=0; ii<models.size(); ii++) {
         auto &model = models[ii];
	 if (shader_for_tmeshes_p) {
	    shader_for_tmeshes_p->Use();
	    model.draw_tmeshes(shader_for_tmeshes_p, mvp, model_rotation, lights, eye_position, bg_col_v4,
			       do_depth_fog);
         }

         // now the coloured vertices mesh (for molecule things, not textured things)
         if (shader_for_meshes_p) {
            float opacity = 1.0f;
            model.draw_meshes(shader_for_meshes_p, mvp, model_rotation, lights, eye_position, opacity, bg_col_v4, do_depth_fog);
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
			dummy_lights, dummy_eye_position, opacity, bg_col_v4, false, false);
	 }
      }
   }
}

void
graphics_info_t::draw_molecules_for_shadow_map(unsigned int light_index) {

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_molecules_for_shadow_map() -- start -- " << err << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   it = lights.find(light_index);
   if (it != lights.end()) {
      const auto &light = it->second;
      glm::mat4 mvp_orthogonal = get_mvp_for_shadow_map(light.direction);
      glm::mat4 model_rotation = get_model_rotation();
      glm::vec4 bg_col_v4(background_colour, 1.0f);

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
                                  mvp_orthogonal, model_rotation, lights, dummy_eye_position,
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
               m.molecule_as_mesh.draw(&shader_for_meshes_shadow_map,
                                       mvp_orthogonal, model_rotation, lights, dummy_eye_position,
                                       opacity, bg_col_v4, gl_lines_mode, do_depth_fog, show_just_shadows);
            }
         }
      }
   }

   // ribbons

   draw_meshed_generic_display_object_meshes(PASS_TYPE_FOR_SHADOWS); // draws the meshes in a molecules std::vector<Messh> meshes;

   // more stuff needs to be added here - see draw_molecules_for_ssao() (where they haven't been added either)
}
