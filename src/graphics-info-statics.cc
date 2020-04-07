
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

GtkWidget *graphics_info_t::main_window = NULL;


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
