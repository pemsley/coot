
#ifdef USE_PYTHON
#include "Python.h"
#endif
#include "graphics-info.h"

bool graphics_info_t::residue_type_selection_was_user_picked_residue_range = false;

bool graphics_info_t::make_auto_h_bond_restraints_flag = false;

bool graphics_info_t::do_rotamer_restraints = false;

bool graphics_info_t::do_debug_refinement = false;

std::atomic<bool> graphics_info_t::on_going_updating_map_lock(false);

std::map<keyboard_key_t, key_bindings_t> graphics_info_t::key_bindings_map;

std::string graphics_info_t::mtz_file_for_refmac;

bool graphics_info_t::convert_dictionary_planes_to_improper_dihedrals_flag = false;

float graphics_info_t::screen_z_near_perspective =  15.0;
float graphics_info_t::screen_z_far_perspective  = 400.0;

GtkWidget *graphics_info_t::main_window = NULL;

// GLuint graphics_info_t::programID_for_maps = 0; in a shader now  - as
//programID_for_central_cube should be
Shader graphics_info_t::shader_for_maps;
Shader graphics_info_t::shader_for_models;
Shader graphics_info_t::shader_for_central_cube;
Shader graphics_info_t::shader_for_origin_cube;
Shader graphics_info_t::shader_for_hud_text;
Shader graphics_info_t::shader_for_screen;
Shader graphics_info_t::shader_for_blur;
std::chrono::time_point<std::chrono::system_clock> graphics_info_t::previous_frame_time = std::chrono::high_resolution_clock::now();
std::chrono::time_point<std::chrono::system_clock> graphics_info_t::previous_frame_time_for_per_second_counter = std::chrono::high_resolution_clock::now();
long graphics_info_t::frame_counter = 0;
long graphics_info_t::frame_counter_at_last_display = 0;
std::queue<std::chrono::time_point<std::chrono::system_clock> > graphics_info_t::frame_draw_queue;

glm::vec3 graphics_info_t::eye_position = glm::vec3(0,0,0);
std::map<unsigned int, gl_lights_info_t> graphics_info_t::lights;

molecule_class_info_t graphics_info_t::moving_atoms_molecule;

bool graphics_info_t::vera_font_loaded = false;

// static
void
graphics_info_t::make_gl_context_current(bool gl_context_current_request_index) {
   if (glareas.empty()) return;
   if (display_mode_use_secondary_p()) {
      if (gl_context_current_request_index == GL_CONTEXT_SECONDARY) {
         if (glareas.size() > 1) {
            GtkWidget *glarea = glareas[1];
            if (glarea) {
               make_current_gl_context(glarea);
            }
         }
      }
      if (gl_context_current_request_index == GL_CONTEXT_MAIN) {
         GtkWidget *glarea = glareas[0];
	 if (glarea) {
            make_current_gl_context(glarea);
	 }
      }
   } else {
      if (gl_context_current_request_index == GL_CONTEXT_MAIN) {
         GtkWidget *glarea = glareas[0];
	 if (glarea) {
            make_current_gl_context(glarea);
	 }
      }
   }
}
