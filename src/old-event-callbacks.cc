
#include <gtk/gtk.h>

gboolean
on_glarea_scroll(GtkWidget *widget, GdkEventScroll *event) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
   int direction = 1;
   if (event->direction == GDK_SCROLL_UP)
      direction = -1;

   graphics_info_t g;
   bool handled = false;
   bool control_is_pressed = false;
   bool   shift_is_pressed = false;
   if (event->state & GDK_CONTROL_MASK) control_is_pressed = true;
   if (event->state & GDK_SHIFT_MASK) shift_is_pressed = true;

   if (control_is_pressed) {
      if (shift_is_pressed){
         if (direction == 1)
            change_model_molecule_representation_mode(-1);
         else
            change_model_molecule_representation_mode(1);
         handled = true;
      }
   }

   if (! handled) {

      if (shift_is_pressed) {
         graphics_info_t::scroll_zoom(direction);
      } else {
         // scroll density

         // start the idle function - why is this needed? The contouring used to
         // work (i.e. the idle function was added somewhere (else)).
         if (graphics_info_t::glareas.size() > 0) {
            g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
         }
         g.contour_level_scroll_scrollable_map(direction);
      }
   }
#endif
   return TRUE;
}


gboolean
on_glarea_button_press(GtkWidget *widget, GdkEventButton *event) {

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

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE FIXME mouse
#else
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
   if (event->type == GDK_2BUTTON_PRESS)
      was_a_double_click = true;

   // 20220416-PE  I need the state because that tells me which button was pressed.
   //              How do I do that if I don't use gdk_window_get_pointer()?
   GdkModifierType state;
   gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state); // deprecated.

   bool handled = false;

   // if (true) { // check here for left-mouse click
   // if (event->state & GDK_BUTTON1_PRESS) {
   if (true) { // check here for left-mouse click

      if (false)
         std::cout << "click event: " << event->x << " " << event->y << " "
                   << x_as_int << " " << y_as_int << std::endl;

      // first thing to test is the HUD bar
      handled = g.check_if_hud_bar_clicked(event->x, event->y);

      if (! handled)
         handled = check_if_refinement_dialog_arrow_tab_was_clicked();

      if (! handled) {

         // OK...
         { // rama plot click
            GtkAllocation allocation;
            gtk_widget_get_allocation(widget, &allocation);
            int w = allocation.width;
            int h = allocation.height;
            auto rama_plot_hit = g.gl_rama_plot.get_mouse_over_hit(event->x, event->y, w, h);
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
                  add_status_bar_text(message.c_str());

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

            if (event->button == 1) // event->state & GDK_BUTTON1_MASK didn't work because event->state
                                    // was 16 GDK_MOD2_MASK (I don't know why)
               handled = g.check_if_hud_button_moused_over(event->x, event->y, true);
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

         if (state & mouse_pick_button_mask) {
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

   if (! handled)
      g.check_if_in_range_defines(event, mask);
#endif
   return TRUE;
}

gboolean
on_glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
   // 20220528-PE FIXME mouse
   graphics_info_t g;
   if (graphics_info_t::in_moving_atoms_drag_atom_mode_flag) {
      g.unset_moving_atoms_currently_dragged_atom_index();
      g.do_post_drag_refinement_maybe();
      graphics_info_t::in_moving_atoms_drag_atom_mode_flag = 0;
   }

   if (event->state & GDK_BUTTON1_MASK)
      g.check_if_hud_button_clicked(event->x, event->y);

   if (event->state & GDK_BUTTON2_MASK) {
      double delta_x = g.GetMouseClickedX() - event->x;
      double delta_y = g.GetMouseClickedY() - event->y;
      if (std::abs(delta_x) < 10.0) {
         if (std::abs(delta_y) < 10.0) {
            pick_info nearest_atom_index_info = g.atom_pick_gtk3(false);
            if (nearest_atom_index_info.success == GL_TRUE) {
               g.setRotationCentre(nearest_atom_index_info.atom_index,
                                   nearest_atom_index_info.imol);
               g.add_picked_atom_info_to_status_bar(nearest_atom_index_info.imol,
                                                    nearest_atom_index_info.atom_index);
            } else {
               if (g.show_symmetry) {
                  coot::Symm_Atom_Pick_Info_t sap = g.symmetry_atom_pick();
                  if (sap.success) {
                     coot::Cartesian pos = sap.hybrid_atom.pos;
                     g.setRotationCentre(pos);
                     g.add_picked_atom_info_to_status_bar(sap.imol, sap.atom_index);
                     g.molecules[sap.imol].add_atom_to_labelled_symm_atom_list(sap.atom_index,
                                                                               sap.symm_trans,
                                                                               sap.pre_shift_to_origin);
                  }
               }
            }
         }
      }
   }
#endif
   return TRUE;
}

gboolean
on_glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME mouse
#else

   auto check_for_hud_bar_tooltip = [widget] (double event_x, double event_y) {
                               graphics_info_t g;
                               std::pair<bool, mmdb::Atom *> handled_pair = g.check_if_moused_over_hud_bar(event_x, event_y);

                               if (handled_pair.first) {
                                  g.draw_hud_tooltip_flag = true;

                                  // gtk mouse position to OpenGL (clip?) coordinates
                                  GtkAllocation allocation;
                                  gtk_widget_get_allocation(widget, &allocation);
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

   auto check_for_hud_refinemement_dialog_arrow_mouse_over = [widget] (double mouse_x, double mouse_y) {
                                                                 graphics_info_t g;
                                                                 // set hud_refinement_dialog_arrow_is_moused_over as needed.
                                                                 g.hud_refinement_dialog_arrow_is_moused_over = false; // initially
                                                                 if (g.showing_intermediate_atoms_from_refinement()) {
                                                                    GtkAllocation allocation;
                                                                    gtk_widget_get_allocation(widget, &allocation);
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
   gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);

   bool control_is_pressed = false;
   bool   shift_is_pressed = false;
   if (event->state & GDK_CONTROL_MASK) control_is_pressed = true;
   if (event->state & GDK_SHIFT_MASK) shift_is_pressed = true;

   g.mouse_current_x = event->x;
   g.mouse_current_y = event->y;

   check_for_hud_bar_tooltip(event->x, event->y);

   check_for_hud_refinemement_dialog_arrow_mouse_over(event->x, event->y);

   // if not right mouse pressed:
   if (event->state & GDK_BUTTON3_MASK) {
   } else {
      bool button_1_is_down = false;
      if (event->state & GDK_BUTTON1_MASK) button_1_is_down = true;
      g.check_if_hud_button_moused_over(event->x, event->y, button_1_is_down);
   }

   auto mouse_view_rotate = [control_is_pressed] (GtkWidget *w) {
                               if (control_is_pressed) {
                                  do_drag_pan_gtk3(w);
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
   if (event->state & GDK_BUTTON1_MASK) {

      if (! (event->state & GDK_BUTTON3_MASK)) {

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
      gtk_widget_get_allocation(widget, &allocation);
      int w = allocation.width;
      int h = allocation.height;
      auto rama_plot_hit = g.gl_rama_plot.get_mouse_over_hit(event->x, event->y, w, h);
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

   if (event->state & GDK_BUTTON2_MASK) {
      if (shift_is_pressed) {
         // mouse_view_rotate(widget);
         std::cout << "shift middle mouse - what to do here?" << std::endl;
      } else {
         do_drag_pan_gtk3(widget);          // View Panning
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

   if (event->state & mouse_action_button_mask) {
      double delta_x = event->x - g.GetMouseBeginX();
      double delta_y = event->y - g.GetMouseBeginY();
      if (event->state & mouse_other_button) {
         // chording
         g.mouse_zoom(delta_x, delta_y);
      } else {
         if (! shift_is_pressed) {
            // don't rotate the view if we are in atom drag mode
            if (!g.in_moving_atoms_drag_atom_mode_flag) {
               mouse_view_rotate(widget);
            }
         } else {
            g.mouse_zoom(delta_x, delta_y);
         }
      }
   }

   g.handle_delete_item_curor_change(widget);

   // for next motion
   g.SetMouseBegin(event->x,event->y);
   // gtk_widget_queue_draw(widget);
   g.graphics_draw(); // queue
#endif
   return TRUE;
}


gboolean
on_glarea_key_press_notify(GtkWidget *widget, GdkEventKey *event) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE FIXME mouse
   return 1;
#else

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
   g.add_key_to_history(kbk);

   std::map<keyboard_key_t, key_bindings_t>::const_iterator it = g.key_bindings_map.find(kbk);

   bool found = false;
   if (it != g.key_bindings_map.end()) {
     const key_bindings_t &kb = it->second;
     if (true)
        std::cout << "INFO:: key-binding for key: " << it->first.gdk_key << " : "
                  << it->first.ctrl_is_pressed << " " << kb.description << std::endl;
     handled = kb.run();
     found = true;
   }

   // Don't make a special case for P now.
   // int kv = event->keyval;
   // if (kv == graphics_info_t::update_go_to_atom_from_current_residue_key) {
   // update_go_to_atom_from_current_position();
   // handled = TRUE;
   // }

   if (! found)
      if (! handled)
         std::cout << "on_glarea_key_press_notify() key not found in map: " << event->keyval << std::endl;

   g.check_keyboard_history_for_easter_egg_codes();
   g.graphics_draw(); // queue

   return handled;
#endif

}

gboolean
on_glarea_key_release_notify(GtkWidget *widget, GdkEventKey *event) {

#if (GTK_MAJOR_VERSION == 4)
   return 1;
#else
   graphics_info_t g;

   // We need to check the GDK_KEY_Shift_R also, I guess. Not clear
   // to me how to do that now. Fix later.
   g.shift_is_pressed = false;
   if (event->state & GDK_SHIFT_MASK) g.shift_is_pressed = true;
   if (event->keyval == GDK_KEY_Shift_L) g.shift_is_pressed = true;

   // key release is a very special event  - normally we act on key-press.

#endif
   return TRUE;
}
