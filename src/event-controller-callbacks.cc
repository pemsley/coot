
#include "event-controller-callbacks.hh"
#include "graphics-info.h"

void
graphics_info_t::on_glarea_drag_begin(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   // 20220429-PE is this controller for left-mouse or right-mouse?

   auto check_if_refinement_dialog_arrow_tab_was_clicked = [] () {
                                                              graphics_info_t g;
                                                              gboolean handled = FALSE;
                                                              if (g.hud_refinement_dialog_arrow_is_moused_over) {
                                                                 g.show_refinement_and_regularization_parameters_dialog();
                                                                 g.hud_refinement_dialog_arrow_is_moused_over = false; // job done
                                                                 handled = TRUE;
                                                                 g.graphics_draw(); // unhighlight the arrow
                                                              }
                                                              return gboolean(handled);
                                                           };

   graphics_info_t g;
   g.SetMouseBegin(x,y);
   g.SetMouseClicked(x, y);

   // GdkModifierType mask;
   // gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state); Old-style - keep for grepping
   // GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());

   bool was_a_double_click = false;

   // 20220429-PE we don't do this now
   // if (event->type == GDK_2BUTTON_PRESS)
   // was_a_double_click = true;

   // 20220416-PE  I need the state because that tells me which button was pressed.
   //              How do I do that if I don't use gdk_window_get_pointer()?

   // 20220429-PE we don't do this now
   // GdkModifierType state;
   // gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state); // deprecated.

   // if (true) { // check here for left-mouse click
   // if (event->state & GDK_BUTTON1_PRESS) {
   if (true) { // check here for left-mouse click

      bool handled = false;

      // first thing to test is the HUD bar
      handled = g.check_if_hud_bar_clicked(x, y);

      if (! handled)
         handled = check_if_refinement_dialog_arrow_tab_was_clicked();

      if (! handled) {

         // OK...
         { // rama plot click
            GtkAllocation allocation;
            gtk_widget_get_allocation(gl_area, &allocation);
            int w = allocation.width;
            int h = allocation.height;
            auto rama_plot_hit = g.gl_rama_plot.get_mouse_over_hit(x, y, w, h);
            if (rama_plot_hit.plot_was_clicked) {
               if (rama_plot_hit.residue_was_clicked) {
                  std::cout << "::::::::::::::::: click " << rama_plot_hit.residue_was_clicked << std::endl;
                  std::string message = "Ramachandran plot clicked residue: ";
                  message += rama_plot_hit.residue_spec.chain_id;
                  message += " ";
                  message += std::to_string(rama_plot_hit.residue_spec.res_no);
                  if (! rama_plot_hit.residue_spec.ins_code.empty()) {
                     message += " ";
                     message += rama_plot_hit.residue_spec.ins_code;
                  }
                  g.add_status_bar_text(message);

                  g.set_go_to_residue_intelligent(rama_plot_hit.residue_spec.chain_id,
                                                  rama_plot_hit.residue_spec.res_no,
                                                  rama_plot_hit.residue_spec.ins_code);
                  int success = g.try_centre_from_new_go_to_atom();
                  if (success) {
                     g.update_things_on_move_and_redraw();
                  }
               }
               handled = true;
            }
         }

         // 20210829-PE This should be in *button-release* I think.
         // Here we could check for button-down (to give a "button pressed but not activatetd" look)
         // Also I need to check that right-mouse is not being used before calling this.
         //
         // 20210830-PE OK, let's comment out the button_clicked then, and merely act as if
         // the mouse had been moved when the button is down
         // handled = g.check_if_hud_button_clicked(event->x, event->y);
         //

         // std::cout << "::::::::::::::::::: Here A event type " << event->type << std::endl;
         // std::cout << "::::::::::::::::::: Here A event button " << event->button << std::endl;
         // std::cout << "::::::::::::::::::: Here A debug " << event->state << " " << GDK_BUTTON1_MASK  << std::endl;
         // std::cout << "::::::::::::::::::: Here A debug " << event->state << " " << GDK_BUTTON2_MASK  << std::endl;
         // std::cout << "::::::::::::::::::: Here A debug " << event->state << " " << GDK_BUTTON3_MASK  << std::endl;

         if (! handled) {

            // this is not the place for checking the mouse button
            // a different mouse button gets a different controller.

            //             if (event->button == 1) // event->state & GDK_BUTTON1_MASK didn't work because event->state
                                    // was 16 GDK_MOD2_MASK (I don't know why)

            handled = g.check_if_hud_button_moused_over(x, y, true);
         }
      }

   GdkModifierType mouse_pick_button_mask        = GDK_BUTTON1_MASK;
   GdkModifierType mouse_view_rotate_button_mask = GDK_BUTTON3_MASK;

#ifdef __APPLE__
   mouse_view_rotate_button_mask = GDK_BUTTON1_MASK;
   mouse_pick_button_mask        = GDK_BUTTON1_MASK;
#endif

      if (! handled) {
         // implicit type cast
         // std::cout << "debug event->state " << event->state << " mouse_pick_button_mask " << mouse_pick_button_mask << std::endl;

         // 20220429-PE was: if (state & mouse_pick_button_mask) {
         // 
         if (true) {
            // std::cout << "yes, was a mouse pick button" << std::endl;
            handled = g.check_if_moving_atom_pull(was_a_double_click);
         } else {
            // std::cout << "no, was not a mouse pick button" << std::endl;
         }

         if (! handled) {
            if (was_a_double_click) {
               bool intermediate_atoms_only_flag = false;
               pick_info nearest_atom_index_info = g.atom_pick_gtk3(intermediate_atoms_only_flag);
               if (nearest_atom_index_info.success == GL_TRUE) {
                  handled = true;
                  int im = nearest_atom_index_info.imol;
                  g.molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
                  g.add_picked_atom_info_to_status_bar(im, nearest_atom_index_info.atom_index);
                  g.graphics_draw();
               } else {

                  // try symmetry atom click (c.f. middle button release)
                  //
                  if (g.show_symmetry) {
                     coot::Symm_Atom_Pick_Info_t sap = g.symmetry_atom_pick();
                     if (sap.success) {
                        g.add_picked_atom_info_to_status_bar(sap.imol, sap.atom_index);
                        g.molecules[sap.imol].add_atom_to_labelled_symm_atom_list(sap.atom_index,
                                                                                  sap.symm_trans,
                                                                                  sap.pre_shift_to_origin);
                        g.graphics_draw();
                     }
                  }
               }
            }
         }
      }
   }

   // 20220429-PE Oh - more work check_if_in_range_defines() will need to be rewritten.
   //
   // g.check_if_in_range_defines(event, mask);

}

void
graphics_info_t::on_glarea_drag_update(GtkGestureDrag *gesture, double delta_x, double delta_y, GtkWidget *gl_area) {

   
   auto check_for_hud_bar_tooltip = [gl_area] (double event_x, double event_y) {
                               graphics_info_t g;
                               std::pair<bool, mmdb::Atom *> handled_pair = g.check_if_moused_over_hud_bar(event_x, event_y);

                               if (handled_pair.first) {
                                  g.draw_hud_tooltip_flag = true;

                                  // gtk mouse position to OpenGL (clip?) coordinates
                                  GtkAllocation allocation;
                                  gtk_widget_get_allocation(gl_area, &allocation);
                                  int w = allocation.width;
                                  int h = allocation.height;
                                  float xx =    2.0 * g.mouse_current_x/static_cast<float>(w) - 1.0f;
                                  float yy = - (2.0 * g.mouse_current_y/static_cast<float>(h) - 1.0f);
                                  glm::vec2 pos(xx, yy);
                                  // this makes the top-left of the tooltip bubble point at the hud geometry bar box (mouse position)
                                  // without it, the tooltip middle is at the cursor position
                                  // 0.1  too much to the left
                                  // 0.07 too much to the left
                                  // 0.05 too much to the left (not much)
                                  // 0.0  too much to the right
                                  float ww = 0.04f * (static_cast<float>(w)/900.0 - 1.0); //  hard-coded inital width - hmmm.
                                  glm::vec2 background_texture_offset(0.08f - ww, -0.058f);
                                  glm::vec2 label_texture_offset(0.0f, -0.086f);
                                  glm::vec2 background_texture_pos = pos + background_texture_offset;
                                  glm::vec2 atom_label_position = pos + label_texture_offset;
                                  g.mesh_for_hud_tooltip_background.set_position(background_texture_pos); // used in uniforms
                                  g.tmesh_for_hud_geometry_tooltip_label.set_position(atom_label_position);

                                  mmdb::Atom *at = handled_pair.second;
                                  coot::atom_spec_t at_spec(at);
                                  g.label_for_hud_geometry_tooltip = at_spec.simple_label(at->residue->GetResName()); // e.g. A 65 CA
                                  g.active_atom_for_hud_geometry_bar = at;
                                  graphics_draw();
                                  // return TRUE;
                               } else {
                                  g.draw_hud_tooltip_flag = false;
                               }
                            };

   auto check_for_hud_refinemement_dialog_arrow_mouse_over = [gl_area] (double mouse_x, double mouse_y) {
                                                                 graphics_info_t g;
                                                                 // set hud_refinement_dialog_arrow_is_moused_over as needed.
                                                                 g.hud_refinement_dialog_arrow_is_moused_over = false; // initially
                                                                 if (g.showing_intermediate_atoms_from_refinement()) {
                                                                    GtkAllocation allocation;
                                                                    gtk_widget_get_allocation(gl_area, &allocation);
                                                                    int w = allocation.width;
                                                                    int h = allocation.height;
                                                                    float xx =    2.0 * mouse_x/static_cast<float>(w) - 1.0f;
                                                                    float yy = - (2.0 * mouse_y/static_cast<float>(h) - 1.0f);
                                                                    // std::cout << "xx " << xx << " yy " << yy << std::endl;
                                                                    float arrow_size = 0.04;
                                                                    if (xx > (1.0 - 2.0 * arrow_size)) {
                                                                       if (yy > (0.9-arrow_size)) {
                                                                          if (yy < (0.9+arrow_size)) {
                                                                             g.hud_refinement_dialog_arrow_is_moused_over = true;
                                                                          }
                                                                       }
                                                                    }
                                                                 }
                                                             };


   graphics_info_t g;

   // split this function up before it gets too big.

   int x_as_int, y_as_int;
   GdkModifierType mask;
   GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());
   GdkDevice *mouse = gdk_seat_get_pointer(seat);
   // gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);

   bool control_is_pressed = false;
   bool   shift_is_pressed = false;

   // old style
   // if (event->state & GDK_CONTROL_MASK) control_is_pressed = true;
   // if (event->state & GDK_SHIFT_MASK) shift_is_pressed = true;

   // Hmmm! 
   // if (modifiers & GDK_SHIFT_MASK) shift_is_pressed = true;
   // if (modifiers & GDK_CONTROL_MASK) control_is_pressed = true;

   mouse_current_x = mouse_clicked_begin.first  + delta_x;
   mouse_current_y = mouse_clicked_begin.second + delta_y;

   check_for_hud_bar_tooltip(mouse_current_x, mouse_current_y);

   check_for_hud_refinemement_dialog_arrow_mouse_over(mouse_current_x, mouse_current_y);

   // 20220429-PE this sort of test is done with a different conntroller... right?
   //
   //
   // if not right mouse pressed:
   // if (event->state & GDK_BUTTON3_MASK) {
   if (false) {
   
   } else {
      bool button_1_is_down = false;

      // 20220429-PE how do we test for this now?
      // if (event->state & GDK_BUTTON1_MASK) button_1_is_down = true;
      g.check_if_hud_button_moused_over(mouse_current_x, mouse_current_y, button_1_is_down);
   }

   auto mouse_view_rotate = [control_is_pressed] (GtkWidget *w) {
                               if (control_is_pressed) {
                                  graphics_info_t g;
                                  g.do_drag_pan_gtk3(w);
                               } else {
                                  GtkAllocation allocation;
                                  gtk_widget_get_allocation(w, &allocation);
                                  int w = allocation.width;
                                  int h = allocation.height;
                                  graphics_info_t::update_view_quaternion(w, h);
                               }
                            };

   // atom pulls with left mouse (but not with right mouse also (that's zoom)
   //
   // if (event->state & GDK_BUTTON1_MASK) {

   // 20220429-PE Hmmm this needs more thinking
   if (false) {

      // if (! (event->state & GDK_BUTTON3_MASK)) {
      if (true) {

         if (g.in_moving_atoms_drag_atom_mode_flag) {
            if (g.last_restraints_size() > 0) {
               // move an already picked atom
               g.move_atom_pull_target_position(x_as_int, y_as_int);
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
      }
   }

   { // rama plot mouse-over
      GtkAllocation allocation;
      gtk_widget_get_allocation(gl_area, &allocation);
      int w = allocation.width;
      int h = allocation.height;
      auto rama_plot_hit = g.gl_rama_plot.get_mouse_over_hit(mouse_current_x, mouse_current_y, w, h);
      if (rama_plot_hit.plot_was_clicked) {
         if (rama_plot_hit.residue_was_clicked) {
            // std::cout << "::::::::::::::::: hit " << rama_plot_hit.second << std::endl;
            std::string message = "Rama plot residue: ";
            message += rama_plot_hit.residue_spec.chain_id;
            message += " ";
            message += std::to_string(rama_plot_hit.residue_spec.res_no);
            if (! rama_plot_hit.residue_spec.ins_code.empty()) {
               message += " ";
               message += rama_plot_hit.residue_spec.ins_code;
            }
            add_status_bar_text(message.c_str());
         }
      }
   }

   //if (event->state & GDK_BUTTON2_MASK) {
   if (false) {
      if (shift_is_pressed) {
         // mouse_view_rotate(widget);
         std::cout << "shift middle mouse - what to do here?" << std::endl;
      } else {
         do_drag_pan_gtk3(gl_area);          // View Panning
      }
   }

   int mouse_action_button_mask = GDK_BUTTON3_MASK;
   int mouse_other_button       = GDK_BUTTON1_MASK;

   // test for being a mac laptop? - or a user setting
#ifdef __APPLE__  // this needs improvement
   if (true) {
      mouse_action_button_mask = GDK_BUTTON1_MASK;
      mouse_other_button       = GDK_BUTTON3_MASK;
   }
#endif

   // if (event->state & mouse_action_button_mask) {
   if (false) {
      // if (event->state & mouse_other_button) {
      if (false) {
         // chording
         g.mouse_zoom(delta_x, delta_y);
      } else {
         if (! shift_is_pressed) {
            // don't rotate the view if we are in atom drag mode
            if (!g.in_moving_atoms_drag_atom_mode_flag) {
               mouse_view_rotate(gl_area);
            }
         } else {
            g.mouse_zoom(delta_x, delta_y);
         }
      }
   }

   g.handle_delete_item_curor_change(gl_area);

   // for next motion
   g.SetMouseBegin(mouse_current_x, mouse_current_y);
   // gtk_widget_queue_draw(widget);
   g.graphics_draw(); // queue


}


void
graphics_info_t::on_glarea_drag_end(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

}

void on_glarea_drag_begin(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   graphics_info_t g;
   g.on_glarea_drag_begin(gesture, x, y, gl_area);

}

void on_glarea_drag_update(GtkGestureDrag *gesture, double delta_x, double delta_y, GtkWidget *gl_area) {

   // was mouse motion callback

   graphics_info_t g;
   g.on_glarea_drag_update(gesture, delta_x, delta_y, gl_area);

}

void on_glarea_drag_end(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   graphics_info_t g;
   g.on_glarea_drag_end(gesture, x, y, gl_area);
}

gboolean
on_key_controller_key_pressed(GtkEventControllerKey *controller,
                              guint                  keyval,
                              guint                  keycode,
                              guint                  modifiers,
                              GtkButton             *button) {

   return gboolean(TRUE);

}

void
on_key_controller_key_released(GtkEventControllerKey *controller,
                               guint                  keyval,
                               guint                  keycode,
                               guint                  modifiers,
                               GtkButton             *button) {

}
