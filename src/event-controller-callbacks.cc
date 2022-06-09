
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
graphics_info_t::on_glarea_drag_update(GtkGestureDrag *gesture, double drag_delta_x, double drag_delta_y, GtkWidget *gl_area) {


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

   // std::cout << "mods control_is_pressed " << control_is_pressed << " shift_is_pressed " << shift_is_pressed << std::endl;

   if (shift_is_pressed) {
      do_view_zoom(drag_delta_x, drag_delta_y);
   } else {
      if (control_is_pressed) {
         do_drag_pan_gtk3(gl_area);
      } else {
         do_view_rotation(drag_delta_x, drag_delta_y);
      }
   }

   if (false)
      std::cout << "in on_glarea_drag_update() mouse_clicked_begin,xy " << mouse_clicked_begin.first << " " << mouse_clicked_begin.second
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
graphics_info_t::on_glarea_drag_end(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

}


void
graphics_info_t::do_drag_pan_gtk3(GtkWidget *widget) {

   // This should be a graphics_info_t function

   GtkAllocation allocation;
   gtk_widget_get_allocation(widget, &allocation);
   int w = allocation.width;
   int h = allocation.height;

   graphics_info_t g;
   glm::mat4 mvp = g.get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion

   float mouseX_1 = g.GetMouseBeginX() / (w * 0.5f) - 1.0f;
   float mouseY_1 = g.GetMouseBeginY() / (h * 0.5f) - 1.0f;
   float mouseX_2 = g.mouse_current_x  / (w * 0.5f) - 1.0f;
   float mouseY_2 = g.mouse_current_y  / (h * 0.5f) - 1.0f;

   glm::mat4 vp_inv = glm::inverse(mvp);

   glm::vec4 screenPos_1 = glm::vec4(mouseX_1, -mouseY_1, 1.0f, 1.0f);
   glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2, 1.0f, 1.0f);
   glm::vec4 worldPos_1 = vp_inv * screenPos_1;
   glm::vec4 worldPos_2 = vp_inv * screenPos_2;

   glm::vec4 delta(worldPos_1 / worldPos_1.w - worldPos_2 / worldPos_2.w);
   glm::vec3 delta_v3(delta);

   g.add_to_rotation_centre(delta_v3);

   // g.update_maps();
   // if (graphics_info_t::glareas.size() > 0)
   // int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);

   g.update_things_on_move(); // 20211013-PE do I need the _and_redraw() version of this function?

}
