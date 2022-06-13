
#include "event-controller-callbacks.hh"
#include "graphics-info.h"


// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
//                              primary
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

void
graphics_info_t::on_glarea_drag_begin_primary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

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

   SetMouseBegin(x,y);
   SetMouseClicked(x, y);
   mouse_x = x;
   mouse_y = y;
   drag_begin_x = x;
   drag_begin_y = y;
   // 2022 code
   set_mouse_previous_position(x,y);

}

// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_primary(GtkGestureDrag *gesture, double drag_delta_x, double drag_delta_y, GtkWidget *gl_area) {

#if 0 // 20220613-PE why is this here now?
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
#endif


   // Ctrl left-mouse means pan
   GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
   bool control_is_pressed = (modifier & GDK_CONTROL_MASK);

   if (control_is_pressed) {
      do_drag_pan_gtk3(gl_area, drag_delta_x, drag_delta_y); // 20220613-PE no redraw here currently
      graphics_draw();
   }
   double x = drag_begin_x + drag_delta_x;
   double y = drag_begin_y + drag_delta_y;
   set_mouse_previous_position(x, y);

}


void
graphics_info_t::on_glarea_drag_end_primary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

}


// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
//                              secondary
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

void
graphics_info_t::on_glarea_drag_begin_secondary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   // 20220429-PE is this controller for left-mouse or right-mouse?

   SetMouseBegin(x,y);
   SetMouseClicked(x, y);
   mouse_x = x;
   mouse_y = y;
   drag_begin_x = x;
   drag_begin_y = y;
   // 2022 code
   set_mouse_previous_position(x,y);
}

// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_secondary(GtkGestureDrag *gesture,
                                                 double drag_delta_x, double drag_delta_y,
                                                 GtkWidget *gl_area) {

   auto do_view_zoom = [] (double drag_delta_x, double drag_delta_y) {
      // std::cout << "calling mouse_zoom with " << drag_delta_x << " " << drag_delta_y << std::endl;
      mouse_zoom(drag_delta_x, drag_delta_y);
   };

   auto do_view_rotation = [gl_area] (double delta_x, double delta_y) {
      GtkAllocation allocation;
      gtk_widget_get_allocation(gl_area, &allocation);
      int w = allocation.width;
      int h = allocation.height;
      update_view_quaternion(w, h, delta_x, delta_y);
   };

   GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(gesture));
   bool control_is_pressed = (modifier & GDK_CONTROL_MASK);
   bool   shift_is_pressed = (modifier & GDK_SHIFT_MASK);

   if (shift_is_pressed) {
      do_view_zoom(drag_delta_x, drag_delta_y);
   } else {
      if (control_is_pressed) {
         do_drag_pan_gtk3(gl_area, drag_delta_x, drag_delta_x);
      } else {
         do_view_rotation(drag_delta_x, drag_delta_y);
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

   double x = drag_begin_x + drag_delta_x;
   double y = drag_begin_y + drag_delta_y;
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

   // 20220429-PE is this controller for left-mouse or right-mouse?

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

   std::cout << "in on_glarea_drag_begin_middle() set previous position and drag_begin to "
             << x << " " << y << std::endl;

}

// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_middle(GtkGestureDrag *gesture,
                                              double drag_delta_x, double drag_delta_y,
                                              GtkWidget *gl_area) {

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
         }
      }
   }

}



void
graphics_info_t::on_glarea_click(GtkGestureClick* self,
                                 gint n_press,
                                 gdouble x,
                                 gdouble y,
                                 gpointer user_data) {

   // n_press can go up to 20, 30...
   //
   if (n_press == 2) { // otherwise triple clicking would toggle the label off, we don't want that.
      bool intermediate_atoms_only_flag = false;
      pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
      int im = naii.imol;
      molecules[im].add_to_labelled_atom_list(naii.atom_index);
      add_picked_atom_info_to_status_bar(im, naii.atom_index);
      graphics_draw();
   }

}



void
graphics_info_t::do_drag_pan_gtk3(GtkWidget *widget, double drag_delta_x, double drag_delta_y) {

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



gboolean
graphics_info_t::on_glarea_key_controller_key_pressed(GtkEventControllerKey *controller,
                                                      guint                  keyval,
                                                      guint                  keycode,
                                                      guint                  modifiers) {


   gboolean handled = false;
   control_is_pressed = (modifiers & GDK_CONTROL_MASK);
   shift_is_pressed   = (modifiers & GDK_SHIFT_MASK);

   keyboard_key_t kbk(keyval, control_is_pressed);
   add_key_to_history(kbk);

   bool found = false;
   std::map<keyboard_key_t, key_bindings_t>::const_iterator it = key_bindings_map.find(kbk);
   if (it != key_bindings_map.end()) {
     const key_bindings_t &kb = it->second;
     if (true)
        std::cout << "INFO:: key-binding for key: " << it->first.gdk_key << " : "
                  << it->first.ctrl_is_pressed << " " << kb.description << std::endl;
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

}
