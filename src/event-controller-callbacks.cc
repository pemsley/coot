/*
 * src/event-controller-callbacks.cc
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

#include "geometry/residue-and-atom-specs.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp> // for to_string()

#include "event-controller-callbacks.hh"
#include "graphics-info.h"
#include "sound.hh"

#include "utils/logging.hh"
extern logging logger;


void play_sound_left_click() {

   // play_sound_file("538554_3725923-lq-Sjonas88-success.ogg");
   // play_sound_file("325112_3246658-lq-fisch12345-success.ogg");
   // play_sound_file("538546_3725923-lq-Sjonas_Rising.ogg");
   // play_sound_file("538548_3725923-lq-Sjonas-Select-3.ogg"); // nice soft click
   // play_sound_file("538549_3725923-lq-Sjonas-Select-2.ogg"); // "tink"
   // play_sound_file("538550_3725923-lq-Sjonas88-Deep-tone.ogg"); // marimba?
   // play_sound_file("538553_3725923-lq-Sjonas88-Stars.ogg"); // high pitch couple of notes

}

// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
//                              primary
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

void
graphics_info_t::on_glarea_drag_begin_primary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   SetMouseBegin(x,y);
   SetMouseClicked(x, y);
   mouse_x = x;
   mouse_y = y;
   drag_begin_x = x;
   drag_begin_y = y;
   // 2022 code
   set_mouse_previous_position(x,y);

   // pick_info nearest_atom_index_info = atom_pick_gtk3(false);

   bool handled = false;
   bool was_a_double_click = false;

   if (false) {       // check here translation_gizmo_picked()?

   } else {
      handled = check_if_moving_atom_pull(was_a_double_click);

      if (! handled) {
         check_if_in_range_defines();
      }
   }

   play_sound_left_click();

}

translation_gizmo_t::pick_info_t
graphics_info_t::translation_gizmo_picked() {

   // 20250719-PE c.f. atom_pick_gtk3()

   translation_gizmo_t::pick_info_t pick_info = translation_gizmo_t::pick_info_t::NONE;
   if (translation_gizmo_mesh.get_draw_this_mesh()) {
      std::cout << "translation gizmo is being drawn" << std::endl;
      graphics_info_t g;
      GtkAllocation allocation = get_glarea_allocation();
      int w = allocation.width;
      int h = allocation.height;
      float mouseX = g.GetMouseBeginX() / (w * 0.5f) - 1.0f;  // should be static?
      float mouseY = g.GetMouseBeginY() / (h * 0.5f) - 1.0f;
      glm::mat4 mvp = get_molecule_mvp();
      glm::mat4 vp_inv = glm::inverse(mvp);
      float real_y = - mouseY; // in range -1 -> 1
      glm::vec4 screenPos_f = glm::vec4(mouseX, real_y, -1.0f, 1.0f);
      glm::vec4 screenPos_b = glm::vec4(mouseX, real_y,  1.0f, 1.0f);
      glm::vec4 worldPos_f = vp_inv * screenPos_f;
      glm::vec4 worldPos_b = vp_inv * screenPos_b;
      float w_scale_f = 1.0/worldPos_f.w;
      float w_scale_b = 1.0/worldPos_b.w;
      coot::Cartesian front(worldPos_f.x * w_scale_f, worldPos_f.y * w_scale_f, worldPos_f.z * w_scale_f);
      coot::Cartesian  back(worldPos_b.x * w_scale_b, worldPos_b.y * w_scale_b, worldPos_b.z * w_scale_b);

      pick_info = translation_gizmo.pick(front, back);
      // translation_gizmo_axis_dragged = pick_info;

      if (pick_info != translation_gizmo_t::pick_info_t::NONE) {
         std::cout << "translation gizmo picked! " << pick_info << std::endl;
      }
   }
   return pick_info;
};


// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_primary(GtkGestureDrag *gesture,
                                               double drag_delta_x, double drag_delta_y,
                                               GtkWidget *gl_area) {

   auto do_view_rotation = [gl_area] (double delta_x, double delta_y) {
      GtkAllocation allocation;
      gtk_widget_get_allocation(gl_area, &allocation);
      int w = allocation.width;
      int h = allocation.height;
      update_view_quaternion(w, h, delta_x, delta_y);
   };

   auto do_view_zoom = [] (double drag_delta_x, double drag_delta_y) {
      mouse_zoom(drag_delta_x, drag_delta_y);
   };

   auto move_translation_gizmo = [] (double screen_delta_x, double screen_delta_y) {

      if (false)
         std::cout << "move_translation_gizmo() " << screen_delta_x << " "  << screen_delta_y << " axis: "
                   << translation_gizmo_axis_dragged << std::endl;

      if (translation_gizmo_axis_dragged != translation_gizmo_t::pick_info_t::NONE) {

         double x = drag_begin_x + screen_delta_x;
         double y = drag_begin_y + screen_delta_y;
         double delta_delta_x = x - get_mouse_previous_position_x();
         double delta_delta_y = y - get_mouse_previous_position_y();

         glm::mat4 model_rotation_matrix = get_model_rotation();
         double sf = 0.00075 * zoom;
         glm::vec4 screen_vec(sf * delta_delta_x, -sf * delta_delta_y, 0.0, 1.0);
         glm::vec4 mol_space_vec = glm::transpose(model_rotation_matrix) * screen_vec;
         coot::Cartesian t(0,0,0);
         if (translation_gizmo_axis_dragged == translation_gizmo_t::pick_info_t::X_AXIS) t = coot::Cartesian(mol_space_vec.x, 0, 0);
         if (translation_gizmo_axis_dragged == translation_gizmo_t::pick_info_t::Y_AXIS) t = coot::Cartesian(0, mol_space_vec.y, 0);
         if (translation_gizmo_axis_dragged == translation_gizmo_t::pick_info_t::Z_AXIS) t = coot::Cartesian(0, 0, mol_space_vec.z);

         // 2025-10-01-PE I need to call setup_draw_for_translation_gizmo() here?
         // That doesn't seem like a good design.
         translation_gizmo.translate(t);
         setup_draw_for_translation_gizmo();

         int tgagdo = translation_gizmo.attached_to_generic_display_object_number;
         if (tgagdo != translation_gizmo_t::UNATTACHED) {
            if (is_valid_generic_display_object_number(tgagdo)) {
               generic_display_objects[tgagdo].translate(t);
            }
         }
      }

   };

   if (false)
      std::cout << "debug:: use_primary_mouse_for_view_rotation_flag "
                << use_primary_mouse_for_view_rotation_flag << std::endl;

   // Ctrl left-mouse means pan
   GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
   bool control_is_pressed = (modifier & GDK_CONTROL_MASK);
   bool   shift_is_pressed = (modifier & GDK_SHIFT_MASK);
   double x = drag_begin_x + drag_delta_x;
   double y = drag_begin_y + drag_delta_y;
   double delta_delta_x = x - get_mouse_previous_position_x();
   double delta_delta_y = y - get_mouse_previous_position_y();

   bool handled = false;

   if (translation_gizmo_mesh.get_draw_this_mesh()) {
      if (translation_gizmo_is_being_dragged) {
         move_translation_gizmo(drag_delta_x, drag_delta_y);
         handled = true;
      }
   }

   if (! handled) {
      if (in_moving_atoms_drag_atom_mode_flag) {
         if (last_restraints_size() > 0) {
            // move an already picked atom
            bool this_atom_is_anchored = false;
            mmdb::Atom *dragged_anchored_atom = nullptr;
            // use molecule-class-info's fixed atom specs to see if this is a fixed atom.

            if (moving_atoms_asc) {
               mmdb::Atom *at = moving_atoms_asc->atom_selection[moving_atoms_currently_dragged_atom_index];
               coot::atom_spec_t at_spec(at);
               for (unsigned int ispec=0; ispec<molecules[imol_moving_atoms].fixed_atom_specs.size(); ispec++) {
                  if (at_spec == molecules[imol_moving_atoms].fixed_atom_specs[ispec]) {
                     this_atom_is_anchored = true;
                     dragged_anchored_atom =  at;
                     break;
                  }
               }
            }

            if (this_atom_is_anchored) {
               std::cout << "debug:: update primary: this atom is anchored! "
                         << coot::atom_spec_t(dragged_anchored_atom) << std::endl;
               // move_dragged_anchored_atom(dragged_anchored_atom)
            } else {
               move_atom_pull_target_position(x, y, control_is_pressed);
            }
            handled = true;
         } else {
         }
      } else {
         if (control_is_pressed) {
            do_drag_pan_gtk3(gl_area, drag_delta_x, drag_delta_y); // 20220613-PE no redraw here currently
            handled = true;
            graphics_draw();
         } else {
            if (shift_is_pressed) {
               do_view_zoom(drag_delta_x, drag_delta_y);
            } else {
               if (use_primary_mouse_for_view_rotation_flag) {
                  do_view_rotation(drag_delta_x, drag_delta_y);
                  graphics_draw();
               } else {
                  // is this logic correct?
                  rotate_chi(delta_delta_x, delta_delta_y); // does its own graphics_draw()
               }
            }
         }
      }
   }

   graphics_draw();
   mouse_current_x = mouse_clicked_begin.first  + drag_delta_x;
   mouse_current_y = mouse_clicked_begin.second + drag_delta_y;

   // we do this in the update of the secondary. Let's do it here too.
   SetMouseBegin(mouse_current_x, mouse_current_y); // not really "begin", but "previous position"
   set_mouse_previous_position(x, y);
}


void
graphics_info_t::on_glarea_drag_end_primary(G_GNUC_UNUSED GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   double xx = drag_begin_x + x;
   double yy = drag_begin_y + y;
   bool hud_clicked = check_if_hud_button_clicked(xx, yy);

   if (!hud_clicked) {
      if (last_restraints_size() > 0) {
         moving_atoms_currently_dragged_atom_index = -1; // because we have dropped it now.
         poke_the_refinement(); // this will remove the pull restraint if the pulled atom position was close to its target.
      }
   }

   translation_gizmo_is_being_dragged = false;

}


// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
//                              secondary
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

void
graphics_info_t::on_glarea_drag_begin_secondary(G_GNUC_UNUSED GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   SetMouseBegin(x,y);
   SetMouseClicked(x, y);
   mouse_x = x;
   mouse_y = y;
   drag_begin_x = x;
   drag_begin_y = y;
   // 2022 code
   set_mouse_previous_position(x,y);

   bool trackpad_drag = false;
#if 0 // try to use "click" event for this
   if (using_trackpad) {
      trackpad_drag = true;
      check_if_in_range_defines(); // this does a pick and looks for distances and angle defines.
   }
#endif
   if (use_primary_mouse_for_view_rotation_flag) {
      bool was_a_double_click = false; // maybe set this correctly?
      bool handled = check_if_moving_atom_pull(was_a_double_click);
   }

}

// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_secondary(GtkGestureDrag *gesture,
                                                 double drag_delta_x, double drag_delta_y,
                                                 GtkWidget *gl_area) {

   if (false)
      std::cout << "graphics_info_t::on_glarea_drag_update_secondary() " << std::endl;

   auto do_view_zoom = [] (double drag_delta_x, double drag_delta_y) {
      mouse_zoom(drag_delta_x, drag_delta_y);
   };

   auto do_view_rotation = [gl_area] (double delta_x, double delta_y) {
      GtkAllocation allocation;
      gtk_widget_get_allocation(gl_area, &allocation);
      int w = allocation.width;
      int h = allocation.height;
      update_view_quaternion(w, h, delta_x, delta_y);
   };

   double x = drag_begin_x + drag_delta_x;
   double y = drag_begin_y + drag_delta_y;

   // std::cout << "drag update_secondary: " << x << " " << y << std::endl;

   GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
   bool control_is_pressed = (modifier & GDK_CONTROL_MASK);
   bool   shift_is_pressed = (modifier & GDK_SHIFT_MASK);

   if (false)
      std::cout << "on_glarea_drag_update_secondary shift is pressed " << shift_is_pressed
                << " control_is_pressed " << control_is_pressed << " "
                << drag_delta_x << " " << drag_delta_y
                << std::endl;

   if (shift_is_pressed) {
      do_view_zoom(drag_delta_x, drag_delta_y);
   } else {
      if (control_is_pressed) {
         do_drag_pan_gtk4(gl_area, drag_delta_x, drag_delta_y);
      } else {
         // zoom with chording. Check both because currently
         // APPLE has primary swapped.
         bool do_chorded_view_zoom = false;
         if (modifier & GDK_BUTTON1_MASK)
            if (modifier & GDK_BUTTON3_MASK)
               do_chorded_view_zoom = true;

         if (do_chorded_view_zoom) {
            do_view_zoom(drag_delta_x, drag_delta_y);
         } else {

            bool old_style_mouse = false;
            if (use_primary_mouse_for_view_rotation_flag)
               old_style_mouse = true;
            bool handled = false;
            if (old_style_mouse) {
               do_view_zoom(drag_delta_x, drag_delta_y);
               handled = true;
            }

            if (! handled) {
               do_view_rotation(view_rotation_per_pixel_scale_factor * drag_delta_x,
                                view_rotation_per_pixel_scale_factor * drag_delta_y);
            }
         }
      }
   }

   if (false)
      std::cout << "in on_glarea_drag_update() mouse_clicked_begin,xy "
                << mouse_clicked_begin.first << " " << mouse_clicked_begin.second
                << " drag_begin " << drag_begin_x << " " << drag_begin_y << std::endl;

   graphics_draw();

   mouse_current_x = mouse_clicked_begin.first  + drag_delta_x;
   mouse_current_y = mouse_clicked_begin.second + drag_delta_y;

   SetMouseBegin(mouse_current_x, mouse_current_y); // not really "begin", but "previous position"

   set_mouse_previous_position(x, y);

}

void
graphics_info_t::on_glarea_drag_end_secondary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

}


// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
//                              middle
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

void
graphics_info_t::on_glarea_drag_begin_middle(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   SetMouseBegin(x,y);
   SetMouseClicked(x, y);
   mouse_x = x;
   mouse_y = y;
   drag_begin_x = x;
   drag_begin_y = y;
   // 2022 code
   set_mouse_previous_position(x,y);

   // note to self middle mouse pick happens on button press, but the
   // recentre happens on button release

   if (false)
      std::cout << "in on_glarea_drag_begin_middle() set previous position and drag_begin to "
                << x << " " << y << std::endl;

}

// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_middle(GtkGestureDrag *gesture,
                                              double drag_delta_x, double drag_delta_y,
                                              GtkWidget *gl_area) {

   do_drag_pan_gtk3(gl_area, drag_delta_x, drag_delta_y); // 20220613-PE no redraw here currently
   graphics_draw();
   double x = drag_begin_x + drag_delta_x;
   double y = drag_begin_y + drag_delta_y;
   set_mouse_previous_position(x, y);
   // std::cout << "drag update_middle: " << x << " " << y << std::endl;
}

void
graphics_info_t::on_glarea_drag_end_middle(GtkGestureDrag *gesture, double drag_delta_x, double drag_delta_y, GtkWidget *gl_area) {

   if (fabs(drag_delta_x) < 5.0) {
      if (fabs(drag_delta_y) < 5.0) {
         pick_info nearest_atom_index_info = atom_pick_gtk3(false);
         if (nearest_atom_index_info.success == GL_TRUE) {
            setRotationCentre(nearest_atom_index_info.atom_index,
                              nearest_atom_index_info.imol);
            add_picked_atom_info_to_status_bar(nearest_atom_index_info.imol,
                                               nearest_atom_index_info.atom_index);
            graphics_grab_focus();
         } else {
            // std::cout << "debug:: on_glarea_drag_end_middle() calling symmetry_atom_pick()" << std::endl;
            coot::Symm_Atom_Pick_Info_t sap = symmetry_atom_pick();
            if (sap.success == GL_TRUE) {
               if (is_valid_model_molecule(sap.imol)) {

                  std::pair<symm_trans_t, Cell_Translation> symtransshiftinfo(sap.symm_trans, sap.pre_shift_to_origin);
                  setRotationCentre(translate_atom_with_pre_shift(molecules[sap.imol].atom_sel,
                                                                  sap.atom_index, symtransshiftinfo));
                  graphics_draw();
                  graphics_grab_focus();
               }
            }
         }
      }
   }
}



void
graphics_info_t::on_glarea_click(GtkGestureClick *controller,
                                 gint n_press,
                                 gdouble x,
                                 gdouble y,
                                 G_GNUC_UNUSED gpointer user_data) {

   auto check_if_refinement_dialog_arrow_tab_was_clicked = [] () {
      graphics_info_t g;
      gboolean handled = FALSE;
      if (g.hud_refinement_dialog_arrow_is_moused_over) {
         g.show_refinement_and_regularization_parameters_frame();
         g.hud_refinement_dialog_arrow_is_moused_over = false; // job done
         handled = TRUE;
         g.graphics_draw(); // unhighlight the arrow
      }
      if (false)
         std::cout << "debug:: check_if_refinement_dialog_arrow_tab_was_clicked() returns " << handled << std::endl;
      return gboolean(handled);
   };

   // no longer useful
   // std::cout << "----------(mouse) click!" << std::endl;

   SetMouseBegin(x,y);

   bool clicked = check_if_hud_bar_clicked(x,y);

   // std::cout << "status for HUD bar clicked: " << clicked << " x " << x << " y " << y << std::endl;

   if (! clicked)
      clicked = check_if_hud_rama_plot_clicked(x,y);

   // std::cout << "status for HUD Rama clicked: " << clicked << std::endl;

   if (!clicked) {

      translation_gizmo_t::pick_info_t pi = translation_gizmo_picked(); // typically NONE
      if (pi != translation_gizmo_t::pick_info_t::NONE) {

         translation_gizmo_is_being_dragged = true;
         translation_gizmo_axis_dragged = pi;

      } else {

         // std::cout << "n_press " << n_press << std::endl;

         // n_press can go up to 20, 30...
         //
         if (n_press == 2) { // otherwise triple clicking would toggle the label off, we don't want that.

            bool handled = false;

            std::cout << "########## double-click!" << std::endl;

            if (in_moving_atoms_drag_atom_mode_flag) {
               if (last_restraints_size() > 0) {
                  handled = check_if_moving_atom_pull(true); // passing was-a-double-click
               }
            }

            if (! handled) {
               bool intermediate_atoms_only_flag = false;
               pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
               if (naii.success) {
                  int imol = naii.imol;
                  molecules[imol].add_to_labelled_atom_list(naii.atom_index);
                  add_picked_atom_info_to_status_bar(imol, naii.atom_index);
                  handled = true;
                  graphics_draw();

               } else {
                  coot::Symm_Atom_Pick_Info_t sap = symmetry_atom_pick();
                  if (sap.success == GL_TRUE) {
                     if (is_valid_model_molecule(sap.imol)) {
                        if (graphics_info_t::molecules[sap.imol].show_symmetry) {
                           int imol = sap.imol;
                           std::pair<symm_trans_t, Cell_Translation> symtransshiftinfo(sap.symm_trans, sap.pre_shift_to_origin);
                           molecules[imol].add_atom_to_labelled_symm_atom_list(sap.atom_index, sap.symm_trans,
                                                                               sap.pre_shift_to_origin);
                           handled = true;
                           graphics_draw();
                        }
                     }
                  }
               }
            }

            if (! handled) {
               bool was_on_a_hud_button = check_if_hud_button_moused_over_or_act_on_hit(x, y, false, true);
               if (! was_on_a_hud_button)
                  blob_under_pointer_to_screen_centre();
            }
         }

         if (n_press == 1) {

            // std::cout << "##################### on_glarea_click() 1 click " << std::endl;

            bool handled = check_if_refinement_dialog_arrow_tab_was_clicked();

            if (! handled) {
               // test for user-defined click here
               if (in_user_defined_define > 0) {
                  bool intermediate_atoms_only_flag = false;
                  std::cout << "DEBUG:: in user-defined " << in_user_defined_define << std::endl;
                  pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
                  if (naii.success) {
                     mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
                     coot::atom_spec_t spec(at);
                     user_defined_atom_pick_specs.push_back(spec);
                     in_user_defined_define -= 1;
                     if (in_user_defined_define == 0) {
                        // run the function then
                        run_user_defined_click_func();
                     }
                  }
                  handled = true;
               }
            }

            if (! handled) {
               GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(controller));
               // std::cout << "debug:: on_glarea_click(); modifier: " << modifier << std::endl;
               if (modifier == 8) { // "option" key on Mac (ALT on PC is 24)
                  bool intermediate_atoms_only_flag = false;
                  pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
                  if (naii.success) {
                     setRotationCentre(naii.atom_index, naii.imol);
                     add_picked_atom_info_to_status_bar(naii.imol, naii.atom_index);
                  }

               } else { // not "option" modifier

                  GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(controller));
                  // std::cout << "debug:: on_glarea_click(); modifier: " << modifier << std::endl;

                  if (tomo_picker_flag) {

                     bool shift_is_pressed = (modifier & GDK_SHIFT_MASK);
                     handled = tomo_pick(x,y, n_press, shift_is_pressed);

                  } else {

                     if (modifier & GDK_SHIFT_MASK) { // shift

                        bool intermediate_atoms_only_flag = false;
                        pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
                        if (naii.success) {
                           int imol = naii.imol;
                           mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[naii.atom_index];
                           molecules[imol].add_to_labelled_atom_list(naii.atom_index);
                           graphics_draw();
                           handled = true;
                        }
                        if (! handled) {
                           coot::Symm_Atom_Pick_Info_t sapi = symmetry_atom_pick();
                           if (sapi.success == GL_TRUE) {
                              int imol = sapi.imol;
                              molecules[imol].add_atom_to_labelled_symm_atom_list(sapi.atom_index, sapi.symm_trans,
                                                                                  sapi.pre_shift_to_origin);
                              graphics_draw();
                           }
                        }

                     } else {

                        if (delete_item_atom == 1) {
                           std::cout << "here A " << std::endl;
                           bool intermediate_atoms_only_flag = false;
                           pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
                           if (naii.success) {
                              std::cout << "here C " << std::endl;
                              // this is convoluted!
                              mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
                              coot::atom_spec_t at_spec(at);
                              std::cout << "here D " << at_spec << std::endl;
                              graphics_info_t::molecules[naii.imol].delete_atom(at_spec);
                              graphics_info_t::graphics_draw();
                              if (modifier & GDK_CONTROL_MASK) {
                                 std::cout << "mulit-pick" << std::endl;
                              } else {
                                 delete_item_atom = 0; // unset
                              }
                           } else {
                              std::cout << "Missed" << std::endl;
                              // do red ring ping here.
                           }
                           handled = true;
                        }

                        // std::cout << "Here with in_range_define " << in_range_define << std::endl;
                        if (! handled) {
                           if (in_range_define == 1 || in_range_define == 2) {
                              bool intermediate_atoms_only_flag = false;
                              pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
                              if (naii.success) {
                                 int imol = naii.imol;
                                 mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[naii.atom_index];
                                 if (in_range_define == 1) {
                                    in_range_first_picked_atom  = coot::atom_spec_t(at);
                                    in_range_first_picked_atom.int_user_data = imol;
                                    molecules[imol].add_to_labelled_atom_list(naii.atom_index);
                                 }
                                 if (in_range_define == 2) {
                                    in_range_second_picked_atom = coot::atom_spec_t(at);
                                    in_range_second_picked_atom.int_user_data = imol;
                                    molecules[imol].add_to_labelled_atom_list(naii.atom_index);
                                 }
                                 in_range_define = 2;
                                 graphics_draw(); // make the label appear
                                 handled =  true;
                              }
                           }
                        }
                     }
                  }

                  if (! handled) {
                     bool intermediate_atoms_only_flag = true;
                     pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
                     if (naii.success) {
                        mmdb::Atom *at = moving_atoms_asc->atom_selection[naii.atom_index];
                        moving_atoms_currently_dragged_atom_index = naii.atom_index;
                        // std::cout << "debug:: in on_glarea_click() picked an intermediate atom " << coot::atom_spec_t(at) << std::endl;
                     }
                  }

                  if (! handled) {

                     // 20240902-PE maybe it should run (and act on the symmtry atom pick) if this is a middle-mouse click?

                     // does this ever run?
                     // 20240902-PE yes it does - maybe it shouldn't.
                     // std::cout << "debug:: click handler: Symmetry atom pick here B - does this run? When? " << std::endl;
                     // coot::Symm_Atom_Pick_Info_t sap = symmetry_atom_pick();
                  }
               }
            }
         }
      }
   }

   graphics_grab_focus(); // 20250615-PE is this a good idea?
}

void
graphics_info_t::do_drag_pan_gtk4(GtkWidget *widget, double drag_delta_x, double drag_delta_y) {

   GtkAllocation allocation;
   gtk_widget_get_allocation(widget, &allocation);
   int w = allocation.width;
   int h = allocation.height;

   graphics_info_t g;
   glm::mat4 mvp = g.get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion

   mouse_current_x = mouse_clicked_begin.first  + drag_delta_x;
   mouse_current_y = mouse_clicked_begin.second + drag_delta_y;

   float mouseX_1 = mouse_current_x / (w * 0.5f) - 1.0f;
   float mouseY_1 = mouse_current_y / (h * 0.5f) - 1.0f;
   float mouseX_2 = get_mouse_previous_position_x() / (w * 0.5f) - 1.0f;
   float mouseY_2 = get_mouse_previous_position_y() / (h * 0.5f) - 1.0f;

   glm::mat4 vp_inv = glm::inverse(mvp);
   glm::vec4 screenPos_1 = glm::vec4(mouseX_1, -mouseY_1, 1.0f, 1.0f);
   glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2, 1.0f, 1.0f);
   glm::vec4 worldPos_1 = vp_inv * screenPos_1;
   glm::vec4 worldPos_2 = vp_inv * screenPos_2;

   glm::vec4 delta(worldPos_1 / worldPos_1.w - worldPos_2 / worldPos_2.w);
   glm::vec3 delta_v3(-delta);

   // it was a typo in the caller of this function, not the static variables!

   if (false)
      std::cout << "mouse-clicked-begin " << mouse_clicked_begin.first << " " << mouse_clicked_begin.second
                << std::endl;
 
   if (false)
      std::cout << "drag_delta_x " << drag_delta_x << " drag_delta_y " << drag_delta_y
                << std::endl;
 
   if (false)
      std::cout << "mouse_current_x " << mouse_current_x << " mouse_current_y " << mouse_current_y
                << std::endl;

   if (false) {
      std::cout << "screen pos delta "
                << mouse_current_x - get_mouse_previous_position_x() << " "
                << mouse_current_y - get_mouse_previous_position_y() << std::endl;
   }
   if (false) {
      std::cout << "in do_drag_pan_gtk4() mouse-current: "
                << mouse_current_x << " " << mouse_current_y << " "
                << "mouse-prev: " << get_mouse_previous_position_x() << " " << get_mouse_previous_position_y()
                << " mouse-1 " << mouseX_1 << " " << mouseY_1
                << " mouse-2 " << mouseX_2 << " " << mouseY_2
                << " delta " << glm::to_string(delta) << " delta_v3 " << glm::to_string(delta_v3) << std::endl;
   }

   g.add_to_rotation_centre(delta_v3);

   // g.update_maps();
   // if (graphics_info_t::glareas.size() > 0)
   // int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);

   g.update_things_on_move(); // 20211013-PE do I need the _and_redraw() version of this function?

   set_mouse_previous_position(mouse_current_x, mouse_current_y); // for next round
}


void
graphics_info_t::do_drag_pan_gtk3(GtkWidget *widget, double drag_delta_x, double drag_delta_y) {

   // who calls this function now?

   // std::cout << "do_drag_pan_gtk3() " << std::endl;

   // This should be a graphics_info_t function

   GtkAllocation allocation;
   gtk_widget_get_allocation(widget, &allocation);
   int w = allocation.width;
   int h = allocation.height;

   graphics_info_t g;
   glm::mat4 mvp = g.get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion

   mouse_current_x = mouse_clicked_begin.first  + drag_delta_x;
   mouse_current_y = mouse_clicked_begin.second + drag_delta_y;

   float mouseX_1 = mouse_current_x / (w * 0.5f) - 1.0f;
   float mouseY_1 = mouse_current_y / (h * 0.5f) - 1.0f;
   float mouseX_2 = get_mouse_previous_position_x() / (w * 0.5f) - 1.0f;
   float mouseY_2 = get_mouse_previous_position_y() / (h * 0.5f) - 1.0f;

   glm::mat4 vp_inv = glm::inverse(mvp);

   glm::vec4 screenPos_1 = glm::vec4(mouseX_1, -mouseY_1, 1.0f, 1.0f);
   glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2, 1.0f, 1.0f);
   glm::vec4 worldPos_1 = vp_inv * screenPos_1;
   glm::vec4 worldPos_2 = vp_inv * screenPos_2;

   glm::vec4 delta(worldPos_1 / worldPos_1.w - worldPos_2 / worldPos_2.w);
   glm::vec3 delta_v3(-delta);

   if (false)
      std::cout << "in do_drag_pan_gtk3() mouse-current: "
                << mouse_current_x << " " << mouse_current_y << " "
                << "mouse-prev: " << get_mouse_previous_position_x() << " " << get_mouse_previous_position_y()
                << " mouse-1 " << mouseX_1 << " " << mouseY_1
                << " mouse-2 " << mouseX_2 << " " << mouseY_2
                << " delta " << glm::to_string(delta) << " delta_v3 " << glm::to_string(delta_v3) << std::endl;

   g.add_to_rotation_centre(delta_v3);

   // g.update_maps();
   // if (graphics_info_t::glareas.size() > 0)
   // int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);

   g.update_things_on_move(); // 20211013-PE do I need the _and_redraw() version of this function?

}

#include "c-interface.h"

void load_tutorial_model_and_data_ec() {

   std::string p = coot::package_data_dir();
   std::string d = coot::util::append_dir_dir(p, "data");

   std::string pdb_fn = coot::util::append_dir_file(d, "tutorial-modern.pdb");
   std::string mtz_fn = coot::util::append_dir_file(d, "rnasa-1.8-all_refmac1.mtz");

   // int imol = handle_read_draw_molecule_with_recentre(pdb_fn.c_str(), true);

   graphics_info_t g;
   int imol = g.create_molecule();
   float bw = graphics_info_t::default_bond_width;
   int bonds_box_type = graphics_info_t::default_bonds_box_type;
   bool recentre_on_read_pdb_flag = true;
   int istat = g.molecules[imol].handle_read_draw_molecule(imol, pdb_fn,
							  coot::util::current_working_dir(),
							  graphics_info_t::Geom_p(),
							  recentre_on_read_pdb_flag, 0,
							  g.allow_duplseqnum,
							  g.convert_to_v2_atom_names_flag,
							  bw, bonds_box_type, true);

   int imol_map = make_and_draw_map_with_refmac_params(mtz_fn.c_str(), "FWT", "PHWT", "",
                                                       0, 0, 1, "FGMP18", "SIGFGMP18", "FreeR_flag", 1);
   int imol_diff_map = make_and_draw_map(mtz_fn.c_str(), "DELFWT", "PHDELWT", "", 0, 1);

}

gboolean
graphics_info_t::on_glarea_key_controller_key_pressed(GtkEventControllerKey *controller,
                                                      guint                  keyval,
                                                      guint                  keycode,
                                                      guint                  modifiers) {

   auto coot_points_frame_callback = +[] (gpointer user_data) {
      GtkWidget *frame = get_widget_from_builder("coot-points-frame");
      if (frame) {
         gtk_widget_set_visible(frame, FALSE);
      }
      return FALSE;
   };

   // I like this function. It should be the callback of a button that gets added
   // to the toolbar when you add Updating Maps.
   auto test_function = [coot_points_frame_callback] () {
      GtkWidget *frame = get_widget_from_builder("coot-points-frame");
      if (frame) {
         gtk_widget_set_visible(frame, TRUE);
         GSourceFunc cb = G_SOURCE_FUNC(coot_points_frame_callback);
         g_timeout_add(4000, cb, nullptr);
      }
   };


   gboolean handled = false;

   control_is_pressed = (modifiers & GDK_CONTROL_MASK);
   shift_is_pressed   = (modifiers & GDK_SHIFT_MASK);

   if (false)
      std::cout << "DEBUG:: on_glarea_key_controller_key_pressed() keyval: " << keyval << std::endl;

   if (false)
      std::cout << "DEBUG:: on_glarea_key_controller_key_pressed() control_is_pressed " << control_is_pressed
                << " shift_is_pressed " << shift_is_pressed << std::endl;

   // key-bindings for testing
   //
   // if (keyval == 101)  E
   //    test_function();

   // if (keyval == 113)  Q
   //    load_tutorial_model_and_data_ec(); // ec: event-controller

   keyboard_key_t kbk(keyval, control_is_pressed);
   add_key_to_history(kbk);

   bool found = false;
   std::map<keyboard_key_t, key_bindings_t>::const_iterator it = key_bindings_map.find(kbk);
   if (it != key_bindings_map.end()) {
     const key_bindings_t &kb = it->second;
     if (true) {
        // std::cout << "INFO:: key-binding for key: " << it->first.gdk_key << " : "
        // << it->first.ctrl_is_pressed << " " << kb.description << std::endl;
        keyboard_key_t key = it->first;
        logger.log(log_t::INFO, {std::string("key-binding for key:"), key.gdk_key, it->first.ctrl_is_pressed,
                                 kb.description});
     }
     handled = kb.run();
     found = true;
   }

   if (! found)
      std::cout << "on_glarea_key_controller_key_pressed() key not found in map: " << keyval << std::endl;

   graphics_draw();
   return gboolean(handled);
}

void
graphics_info_t::on_glarea_key_controller_key_released(GtkEventControllerKey *controller,
                                                       guint                  keyval,
                                                       guint                  keycode,
                                                       guint                  modifiers) {

   control_is_pressed = (modifiers & GDK_CONTROL_MASK);
   shift_is_pressed   = (modifiers & GDK_SHIFT_MASK);

   // 20251122-PE
   delete_item_atom = 0; // turn off pick-delete, otherwise if Ctrl is released, delete_item_atom
                         // is still active, leading to delete of atom on next atom pick - which
                         // is probably not what is wanted.

}


void
graphics_info_t::on_glarea_scrolled(GtkEventControllerScroll *controller,
                                    double                    dx,
                                    double                    dy,
                                    gpointer                  user_data) {

   auto do_mouse_zoom = [] (double dy) {
      int dir = 1;
      if (dy > 0) dir = -1;
      // mouse_zoom(zz, 0.0); // 20250519-PE don't call mouse zoom - it uses drag argument
      // call scroll_zoom
      scroll_zoom(dir);
   };

   // std::cout << "debug:: on_glarea_scrolled() --- start --- dy: " << dy << std::endl;

   GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(controller));
   control_is_pressed = (modifier & GDK_CONTROL_MASK);
   shift_is_pressed = (modifier & GDK_SHIFT_MASK);

   bool handled = false;
   if (false)
      std::cout << "on_glarea_scrolled() control_is_pressed " << control_is_pressed
                << " shift_is_pressed " << shift_is_pressed << std::endl;

   if (control_is_pressed) {
      if (shift_is_pressed) {
         if (dy > 0)
            change_model_molecule_representation_mode(-1);
         else
            change_model_molecule_representation_mode(1);
         graphics_draw();
         handled = true;
      } else {

         // 20250519-PE this is how it used to be! Seems esoteric - I am not
         // sure what it actually does
         if (false) {
            // dy is either 1.0 or -1.0
            // std::cout << "change the proportional editing " << dx << " " << dy << std::endl;
            bool dir = false;
            if (dy < 0.0) dir = true;
            pull_restraint_neighbour_displacement_change_max_radius(dir);
            graphics_draw();
            handled = true;
         }

         // ctrl-scroll zoom, same as shift-scroll zoom
         if (true) {
            do_mouse_zoom(dy);
            handled = true;
         }
      }
   }

   if (! handled) {
      if (shift_is_pressed) {

            do_mouse_zoom(dy);
            handled = true;

      } else {
         // scroll density

         // start the idle function - why is this needed? The contouring used to
         // work (i.e. the idle function was added somewhere (else)).
         if (graphics_info_t::glareas.size() > 0) {
            g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
         }
         contour_level_scroll_scrollable_map(dy);
      }
   }
   // std::cout << "debug:: on_glarea_scrolled() done " << std::endl;
}


void
graphics_info_t::on_glarea_motion(G_GNUC_UNUSED GtkEventControllerMotion* controller,
                                  gdouble x,
                                  gdouble y,
                                  G_GNUC_UNUSED gpointer user_data) {


   // The widget here is the glarea. Pass the height and width of the glarea instead.
   //
   auto check_for_hud_refinemement_dialog_arrow_mouse_over = [] (double mouse_x, double mouse_y, int w, int h) {

      // set hud_refinement_dialog_arrow_is_moused_over as needed.

      bool state = false;
      if (showing_intermediate_atoms_from_refinement()) {
         float xx =    2.0 * mouse_x/static_cast<float>(w) - 1.0f;
         float yy = - (2.0 * mouse_y/static_cast<float>(h) - 1.0f);
         // std::cout << "xx " << xx << " yy " << yy << std::endl;
         float arrow_size = 0.04;
         if (xx > (1.0 - 2.0 * arrow_size)) {
            if (yy > (0.9-arrow_size)) {
               if (yy < (0.9+arrow_size)) {
                  state = true;
               }
            }
         }
      }
      if (state != hud_refinement_dialog_arrow_is_moused_over) {
         hud_refinement_dialog_arrow_is_moused_over = state;
         graphics_draw();
      }
   };

   // So that I can change the highlighting for the moused-over HUD buttons.

   // We can't easily use mouse_x_m mouse_y because they are used by update_view_quaternion().

   // mouse_x = x;
   // mouse_y = y;

   mouse_current_x = x;
   mouse_current_y = y;

   // set_mouse_previous_position(x, y);

   GtkAllocation allocation;
   gtk_widget_get_allocation(glareas[0], &allocation);
   int w = allocation.width;
   int h = allocation.height;

   check_if_hud_button_moused_over(x, y, false);
   check_for_hud_refinemement_dialog_arrow_mouse_over(x, y, w, h);

}



void
graphics_info_t::change_model_molecule_representation_mode(int up_or_down) {

   auto debug_line = [] (const std::string &m) {
      if (false)
         std::cout << m << std::endl;
   };

   // enum { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
   //        COLOUR_BY_CHAIN_BONDS=3,
   //        CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
   //        BONDS_NO_HYDROGENS=15,
   //        CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
   //        CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
   //        CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
   //        COLOUR_BY_MOLECULE_BONDS=8,
   //        COLOUR_BY_RAINBOW_BONDS=9,
   //        COLOUR_BY_B_FACTOR_BONDS=10,
   //        COLOUR_BY_OCCUPANCY_BONDS=11,
   //        COLOUR_BY_USER_DEFINED_COLOURS_BONDS=12 };

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      int bond_type = molecules[imol].Bonds_box_type();

      // up_or_down  1 means scroll up
      // up_or_down -1 means scroll down

      if (up_or_down == 1) {
         if (bond_type == coot::NORMAL_BONDS) { /* atom type bonds */
            molecules[imol].occupancy_representation();
         }
         if (bond_type == coot::COLOUR_BY_CHAIN_BONDS) {
            bool force_rebonding = false;
            molecules[imol].bond_representation(Geom_p(), force_rebonding);
         }
         if (bond_type == coot::COLOUR_BY_MOLECULE_BONDS) {
            // graphics_to_colour_by_chain(imol);
            bool force_rebonding = false;
            molecules[imol].make_colour_by_chain_bonds(force_rebonding);
         }
         if (bond_type == coot::CA_BONDS) {
            // graphics_to_colour_by_molecule(imol);
            bool force_rebonding = false;
            molecules[imol].make_colour_by_molecule_bonds(force_rebonding);
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS) {
            // graphics_to_ca_representation(imol);
            bool force_rebonding = false;
            molecules[imol].ca_representation(force_rebonding);
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR) {
            // graphics_to_ca_plus_ligands_representation(imol);
            bool force_rebonding = false;
            molecules[imol].ca_plus_ligands_representation(Geom_p(), force_rebonding);
         }
         if (bond_type == coot::COLOUR_BY_RAINBOW_BONDS) {
            // graphics_to_ca_plus_ligands_sec_struct_representation(imol);
            molecules[imol].ca_plus_ligands_sec_struct_representation(Geom_p());
         }
         if (bond_type == coot::BONDS_SEC_STRUCT_COLOUR) {
            // graphics_to_rainbow_representation(imol);
            molecules[imol].ca_plus_ligands_rainbow_representation(Geom_p());
         }
         if (bond_type == coot::BONDS_NO_WATERS) {
            // graphics_to_sec_struct_bonds_representation(imol);
            molecules[imol].bonds_sec_struct_representation();
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR) {
            // graphics_to_bonds_no_waters_representation(imol);
            molecules[imol].bonds_no_waters_representation();
         }
         if (bond_type == coot::COLOUR_BY_B_FACTOR_BONDS) {
            // graphics_to_b_factor_cas_representation(imol);
            molecules[imol].b_factor_representation_as_cas();
         }
         if (bond_type == coot::COLOUR_BY_OCCUPANCY_BONDS) {
            // graphics_to_b_factor_representation(imol);
            molecules[imol].b_factor_representation();
         }
      }

      if (up_or_down == -1) {
         if (bond_type == coot::NORMAL_BONDS) {
            // graphics_to_colour_by_chain(imol);
            debug_line("A");
            bool force_rebonding = false;
            molecules[imol].make_colour_by_chain_bonds(force_rebonding);
         }
         if (bond_type == coot::COLOUR_BY_CHAIN_BONDS) {
            // graphics_to_colour_by_molecule(imol);
            bool force_rebonding = false;
            debug_line("B");
            molecules[imol].make_colour_by_molecule_bonds(force_rebonding);
         }
         if (bond_type == coot::COLOUR_BY_MOLECULE_BONDS) {
            // graphics_to_ca_representation(imol);
            bool force_rebonding = false;
            debug_line("C");
            molecules[imol].ca_representation(force_rebonding);
         }
         if (bond_type == coot::CA_BONDS) {
            // graphics_to_ca_plus_ligands_representation(imol);
            debug_line("D");
            bool force_rebonding = false;
            molecules[imol].ca_plus_ligands_representation(Geom_p(), force_rebonding);
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS) {
            // graphics_to_ca_plus_ligands_sec_struct_representation(imol);
            debug_line("E");
            molecules[imol].ca_plus_ligands_sec_struct_representation(Geom_p());
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR) {
            // graphics_to_rainbow_representation(imol);
            debug_line("F");
            molecules[imol].ca_plus_ligands_rainbow_representation(Geom_p());
         }
         if (bond_type == coot::COLOUR_BY_RAINBOW_BONDS) {
            // graphics_to_sec_struct_bonds_representation(imol);
            debug_line("G");
            molecules[imol].bonds_sec_struct_representation();
         }
         if (bond_type == coot::BONDS_SEC_STRUCT_COLOUR) {
            // graphics_to_bonds_no_waters_representation(imol);
            debug_line("H");
            molecules[imol].bonds_no_waters_representation();
         }
         if (bond_type == coot::BONDS_NO_WATERS) {
            //graphics_to_b_factor_cas_representation(imol);
            debug_line("I");
            molecules[imol].b_factor_representation_as_cas();
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR) {
            // graphics_to_b_factor_representation(imol);
            debug_line("J");
            molecules[imol].b_factor_representation();
         }
         if (bond_type == coot::COLOUR_BY_B_FACTOR_BONDS) {
            // graphics_to_occupancy_representation(imol);
            debug_line("K");
            molecules[imol].occupancy_representation();
         }
         if (bond_type == coot::COLOUR_BY_OCCUPANCY_BONDS) {
            // graphics_to_bonds_representation(imol);
            debug_line("L");
            bool force_rebonding = false;
            molecules[imol].bond_representation(Geom_p(), force_rebonding);
         }
      }
   }
}
