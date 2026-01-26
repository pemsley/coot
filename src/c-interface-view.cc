/* src/c-interface.cc
 * 
 * Copyright 2010 by The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */


#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#ifdef USE_GUILE
#include <cstddef> // define std::ptrdiff_t for clang
#include <libguile.h>
#endif

#include "cc-interface.hh"
#include "graphics-info.h"

#include "c-interface.h"

/*  ----------------------------------------------------------------------- */
/*                  Functions for FLEV layout callbacks                     */
/*  ----------------------------------------------------------------------- */
// orient the graphics somehow so that the interaction between
// central_residue and neighbour_residue is perpendicular to screen z.
void orient_view(int imol,
		 const coot::residue_spec_t &central_residue_spec, // ligand typically
		 const coot::residue_spec_t &neighbour_residue_spec) {

   graphics_info_t g;
   if (g.use_graphics_interface_flag) {
      if (g.is_valid_model_molecule(imol)) {
	 try { 
	    clipper::Coord_orth vec = g.molecules[imol].get_vector(central_residue_spec,
								   neighbour_residue_spec);

            // Use coot::ScreenVectors here (remove duplication)

            glm::vec4 glm_back  = g.unproject(1.0);
            glm::vec4 glm_front = g.unproject(0.0);
	    coot::Cartesian b(glm_back.x, glm_back.y, glm_back.z);
	    coot::Cartesian f(glm_front.x, glm_front.y, glm_front.z);
	    coot::Cartesian vec_cart(vec);
	    coot::Cartesian b_to_f_cart = f - b;

	    glm::vec4 glm_centre = g.unproject(0, 0, 0.5);
	    glm::vec4 glm_right  = g.unproject(1, 0, 0.5);
	    glm::vec4 glm_top    = g.unproject(0, 1, 0.5);

	    coot::Cartesian centre(glm_centre.x, glm_centre.y, glm_centre.z);
	    coot::Cartesian front(glm_front.x, glm_front.y, glm_front.z);
	    coot::Cartesian right(glm_right.x, glm_right.y, glm_right.z);
	    coot::Cartesian top(glm_top.x, glm_top.y, glm_top.z);

	    coot::Cartesian screen_x = (right - centre);
	    coot::Cartesian screen_y = (top   - centre);
	    coot::Cartesian screen_z = (front - centre);

	    screen_x.unit_vector_yourself();
	    screen_y.unit_vector_yourself();
	    screen_z.unit_vector_yourself();

	    double a_s_x = dot_product(screen_x, vec_cart);
	    double a_s_z = dot_product(screen_z, vec_cart);

	    double theta = atan2(a_s_z, a_s_x) * 0.5;

	    if (false)
	       std::cout << "theta: " << clipper::Util::rad2d(theta) << " degrees with asx "
			 << a_s_x << " and asz " << a_s_z << std::endl;

	    rotate_y_scene(100, 0.01 * clipper::Util::rad2d(theta));


	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << rte.what() << std::endl;
	 }
      }
   }
}

#include "screendump-tga.hh"

/*  ----------------------------------------------------------------------- */
/*                  Mew Screendump                                          */
/*  ----------------------------------------------------------------------- */
void screendump_tga(const std::string &file_name) {

   graphics_info_t g;
   bool do_screendump = true;
   g.render(do_screendump, file_name);

}


/*  ----------------------------------------------------------------------- */
/*                         perspective,blur,AO on/off */
/*  ----------------------------------------------------------------------- */

// maybe these functions need their own file? Yes. shader-settings.cc

void set_use_perspective_projection(short int state) {

   graphics_info_t::perspective_projection_flag = state;
   graphics_draw();
}

int use_perspective_projection_state() {
   return graphics_info_t::perspective_projection_flag;
}

//! \brief set the perspective fov. Default 20 degrees.
void set_perspective_fov(float degrees) {
   graphics_info_t::perspective_fov = degrees;
   graphics_draw();
}


//! \brief set use ambient occlusion
void set_use_ambient_occlusion(short int state) {
   // user interface is set_use_xxx
   graphics_info_t g;
   g.set_do_ambient_occlusion(state);
   graphics_draw();
}

//! \brief query use ambient occlusion
int use_ambient_occlusion_state() {
   return graphics_info_t::shader_do_ambient_occlusion_flag;
}

//! \brief set use depth blur
void set_use_depth_blur(short int state) {
   // graphics_info_t::shader_do_depth_blur_flag = state; // 2020
   graphics_info_t::shader_do_depth_of_field_blur_flag = state;
   graphics_draw();
}
//! \brief query use depth blur
int use_depth_blur_state() {

   // return graphics_info_t::shader_do_depth_blur_flag; 2020
   return graphics_info_t::shader_do_depth_of_field_blur_flag;

}


//! \brief set use fog
void set_use_fog(short int state) {
   graphics_info_t::shader_do_depth_fog_flag = state;
   graphics_draw();
}

//! \brief query use fog
int use_fog_state() {
   return graphics_info_t::shader_do_depth_fog_flag;
}

void set_use_outline(short int state) {
   graphics_info_t g;
   g.shader_do_outline_flag = state;
   graphics_draw();
}

int use_outline_state() {
   graphics_info_t g;
   return g.shader_do_outline_flag;
}

void set_map_shininess(int imol, float shininess) {
   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].shader_shininess = shininess;
      graphics_draw();
   }
}

void set_map_specular_strength(int imol, float specular_strength) {
   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].shader_specular_strength = specular_strength;
      graphics_draw();
   }
}

void set_map_fresnel_settings(int imol, short int state, float bias, float scale, float power) {

   if (is_valid_map_molecule(imol)) {
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      m.fresnel_settings.update_settings(state, bias, scale, power);
      graphics_draw();
   }

}

void set_draw_normals(short int state) {

   graphics_info_t::draw_normals_flag = state;
   graphics_draw();

}

int  draw_normals_state() {
   return graphics_info_t::draw_normals_flag;
}


void set_draw_mesh(int imol, int mesh_index, short int state) {
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      int size = graphics_info_t::molecules[imol].meshes.size();
      if (mesh_index >= 0 && mesh_index < size) {
         graphics_info_t::molecules[imol].meshes[mesh_index].set_draw_this_mesh(state);
         graphics_info_t::graphics_draw();
      }
   }
}

int draw_mesh_state(int imol, int mesh_index) {
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      int size = graphics_info_t::molecules[imol].meshes.size();
      if (mesh_index >= 0 && mesh_index < size) {
         return graphics_info_t::molecules[imol].meshes[mesh_index].get_draw_this_mesh();
      }
   }
   return -1;
}

void set_map_material_specular(int imol, float specular_strength, float shininess) {

   if (is_valid_map_molecule(imol)) {
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      m.material_for_maps.turn_specularity_on(true);
      m.material_for_maps.specular_strength = specular_strength;
      m.material_for_maps.shininess         = shininess;
      graphics_draw();
   }

}

void set_model_material_specular(int imol, float specular_strength, float shininess) {

   if (is_valid_model_molecule(imol)) {
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      m.material_for_models.specular_strength = specular_strength;
      m.material_for_models.shininess = shininess;

      // m.molecule_as_mesh.set_material_specularity(specular_strength, shininess);
      m.model_molecule_meshes.set_material_specularity(specular_strength, shininess);

      // how about doing this instead of above? (not tested)
      // m.set_material(m.material_for_models);
      graphics_draw();
   }
}

void set_model_material_diffuse(int imol, float r, float g, float b, float a) {

   if (is_valid_model_molecule(imol)) {
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      glm::vec4 d(r,g,b,a);
      m.material_for_models.diffuse = d;
      m.model_molecule_meshes.set_material_diffuse(d);
      graphics_draw();
   }
}

//! \brief set the ambient material multipler - default is 0.2
void set_model_material_ambient(int imol, float r, float g, float b, float a) {

   if (is_valid_model_molecule(imol)) {
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      glm::vec4 ambient(r,g,b,a);
      m.material_for_models.ambient = ambient;
      m.model_molecule_meshes.set_material_ambient(ambient);
   }
   graphics_draw();
}

//! \brief
void set_model_goodselliness(float pastelization_factor) {

   graphics_info_t::goodselliness = pastelization_factor;

   // if the molecule is drawn in goodsell mode, then force a redraw of it
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
         short int f = graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag;
         std::set<int> s; // dummy
         bool g = false; // goodsell_mode
         bool force_rebonding = true;
         // graphics_info_t::molecules[imol].make_colour_by_chain_bonds(s,f,g, force_rebonding);
         set_colour_by_chain_goodsell_mode(imol); // you wanted goodsell mode, right?
      }
   }
   graphics_draw();
}






void reload_map_shader() {

   graphics_info_t g;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(g.glareas[0]));
   std::cout << "reload map shader" << std::endl;
   g.shader_for_maps.init("map.shader", Shader::Entity_t::MAP);
   g.shader_for_meshes.init("meshes.shader", Shader::Entity_t::MAP);
   graphics_draw();
}

void reload_model_shader() {

   graphics_info_t g;
   g.shader_for_models.init("model.shader", Shader::Entity_t::MODEL);
   graphics_draw();
}

void set_atom_radius_scale_factor(int imol, float scale_factor) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_atom_radius_scale_factor(scale_factor);
   }
   graphics_draw();
}

void set_fresnel_colour(int imol, float red, float green, float blue, float opacity) {

   if (is_valid_map_molecule(imol)) {
      glm::vec4 col(red, green, blue, opacity);
      graphics_info_t::molecules[imol].set_fresnel_colour(col);
      graphics_draw();
   }
}

//! \brief
void set_focus_blur_z_depth(float z) {
   graphics_info_t::focus_blur_z_depth = z;
   graphics_draw();
}

//! \brief
void set_focus_blur_strength(float st) {
   graphics_info_t::focus_blur_strength = st;
   graphics_draw();
}

//! \brief set the shadow strength (0.0 to 1.0)
void set_shadow_strength(float strength) {
   graphics_info_t::shadow_strength = strength;
   graphics_draw();
}

//! \brief set the shadow softness (1, 2 or 3)
void set_shadow_texture_resolution_multiplier(unsigned int m) {
   graphics_info_t g;
   g.set_shadow_texture_resolution_multiplier(m);

   graphics_draw();
}


//! \brief
void set_ssao_kernel_n_samples(unsigned int n_samples) {

   graphics_info_t::n_ssao_kernel_samples = n_samples;
   graphics_info_t::generate_ssao_kernel_samples();
   graphics_draw();
}

//! \brief set SSAO strength
void set_ssao_strength(float strength) {
   graphics_info_t::ssao_strength = strength;
   graphics_draw();
}

//! \brief set SSAO strength
void set_ssao_radius(float radius) {
   graphics_info_t::SSAO_radius = radius;
   graphics_draw();
}

//! \brief set SSAO bias
void set_ssao_bias(float bias) {
   graphics_info_t::SSAO_bias = bias;
   graphics_draw();
}

//! \brief set SSAO blur size (0, 1, or 2)
void set_ssao_blur_size(unsigned int blur_size) {
   graphics_info_t::ssao_blur_size = blur_size;
   graphics_draw();
}

//! \brief adjust the effects shader output type (for debugging effects)
void set_effects_shader_output_type(unsigned int type) {
   graphics_info_t::effects_shader_output_type = type;
   graphics_draw();
}

//! \brief adjust the effects shader brightness
void set_effects_shader_brightness(float f) {
   graphics_info_t::effects_brightness = f;
   graphics_draw();
}

//! \brief adjust the effects shader gamma
void set_effects_shader_gamma(float f) {
   graphics_info_t::effects_gamma = f;
   graphics_draw();
}

void set_fps_timing_scale_factor(float f) {

   graphics_info_t::fps_times_scale_factor = f;
   graphics_draw();

}

//! \brief draw background image
void set_draw_background_image(bool state) {
   graphics_info_t::draw_background_image_flag = state;
   graphics_draw();
}

//! \brief set the shadow softness (1, 2 or 3)
void set_shadow_softness(unsigned int softness) {
   graphics_info_t::shadow_softness = softness;
   graphics_draw();
}


//! \brief set the shadow resolution (1,2,3,4)
void set_shadow_resolution(int reso_multiplier) {
   graphics_info_t g;
   g.set_shadow_texture_resolution_multiplier(reso_multiplier);
   graphics_draw();

}
//! \brief set shadow box size - default 66;
void set_shadow_box_size(float size) {
   graphics_info_t g;
   g.shadow_box_size = size;
   graphics_draw();
}



//! \brief set use fancy lighting (default 1 = true);
void set_use_fancy_lighting(short int state) {

   graphics_info_t g;
   if (state) {
      // default (this means - and shadow also)
      g.displayed_image_type = graphics_info_t::SHOW_AO_SCENE;
   } else {
      g.displayed_image_type = graphics_info_t::SHOW_BASIC_SCENE;
   }
   graphics_draw();

}

//! \brief set bond smoothness (default 1 (not smooth))
void set_bond_smoothness_factor(unsigned int fac) {
   graphics_info_t::bond_smoothness_factor = fac;

   // rebonding of the molecules needed here.
   //
   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
         graphics_info_t::molecules[imol].make_glsl_bonds_type_checked(__FUNCTION__);
      }
   }
   graphics_draw();
}


//! \brief set the draw state of the Ramachandran plot display during Real Space Refinement
void set_draw_gl_ramachandran_plot_during_refinement(short int state) {

   graphics_info_t::draw_gl_ramachandran_plot_user_control_flag = state;
   graphics_draw();
}

//! \brief reset the frame buffers
void reset_framebuffers() {

   graphics_info_t g;
   GtkAllocation allocation = g.get_glarea_allocation();
   g.reset_frame_buffers(allocation.width, allocation.height);
   g.graphics_draw();

}






//! \brief set use simple lines for model molecule
void set_use_simple_lines_for_model_molecules(short int state) {

   for (int i=0; i<graphics_n_molecules(); i++) {
      if (is_valid_model_molecule(i)) {
         graphics_info_t::molecules[i].set_draw_model_molecule_as_lines(state);
      }
   }

   // and we need to have a setting in graphics_info_t so that new molecules
   // are draw as lines also. Another time.

   graphics_draw();
}



// testing function
void read_test_gltf_models() {
   graphics_info_t g;
   g.read_test_gltf_models();
   g.graphics_draw();
}

//! \brief load a gltf model
int load_gltf_model(const std::string &gltf_file_name) {
   graphics_info_t g;
   int idx = g.load_gltf_model(gltf_file_name);
   g.graphics_draw();
   return idx;
}

//! \brief set the model animation parameters
void set_model_animation_parameters(unsigned int model_index, float amplitude, float wave_numer, float freq) {

   graphics_info_t g;
   g.set_model_animation_parameters(model_index, amplitude, wave_numer, freq);
}

//! \brief enable/disable the model animation (on or off)
void set_model_animation_state(unsigned int model_index, bool state) {

   graphics_info_t g;
   g.set_model_animation_state(model_index, state);
}


//! \brief load a gltf model
void scale_model(unsigned int model_index, float scale_factor) {

   graphics_info_t g;
   g.scale_model(model_index, scale_factor);
   g.graphics_draw();
}





/*  ----------------------------------------------------------------------- */
/*                         single-model view */
/*  ----------------------------------------------------------------------- */
/*! \name single-model view */
/* \{ */
/*! \brief put molecule number imol to display only model number imodel */
void single_model_view_model_number(int imol, int imodel) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].single_model_view_model_number(imodel);
      graphics_draw();
      std::string s = "Model number ";
      s += coot::util::int_to_string(imodel);
      add_status_bar_text(s.c_str());
   } 
}

/*! \brief the current model number being displayed */
int single_model_view_this_model_number(int imol) {
   int r = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].single_model_view_this_model_number();
      std::string s = "Model number ";
      s += coot::util::int_to_string(r);
      add_status_bar_text(s.c_str());
      graphics_draw();
   }
   return r;
}

/*! \brief the next model number to be displayed */
int single_model_view_next_model_number(int imol) {
   int r = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].single_model_view_next_model_number();
      std::string s = "Model number ";
      s += coot::util::int_to_string(r);
      add_status_bar_text(s.c_str());
      graphics_draw();
   }
   return r;
}

/*! \brief the previous model number to be displayed */
int single_model_view_prev_model_number(int imol) {
   int r = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].single_model_view_prev_model_number();
      std::string s = "Model number ";
      s += coot::util::int_to_string(r);
      add_status_bar_text(s.c_str());
      graphics_draw();
   }
   return r;
}

