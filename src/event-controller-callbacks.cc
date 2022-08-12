
#include "event-controller-callbacks.hh"
#include "graphics-info.h"


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
   handled = check_if_moving_atom_pull(was_a_double_click);

   if (! handled) {
      check_if_in_range_defines();
   }

}

// drag_delta_x and drag_delta_y are the delta coordinates relative to where the drag began.
void
graphics_info_t::on_glarea_drag_update_primary(GtkGestureDrag *gesture, double drag_delta_x, double drag_delta_y, GtkWidget *gl_area) {

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

   if (in_moving_atoms_drag_atom_mode_flag) {
      if (last_restraints_size() > 0) {
         // move an already picked atom
         move_atom_pull_target_position(x, y);
      }
   }

}


void
graphics_info_t::on_glarea_drag_end_primary(G_GNUC_UNUSED GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

   double xx = drag_begin_x + x;
   double yy = drag_begin_y + y;
   bool clicked = check_if_hud_button_clicked(xx, yy);

   std::cout << "hud_button was clicked: " << clicked << std::endl;

}


// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
//                              secondary
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------

void
graphics_info_t::on_glarea_drag_begin_secondary(G_GNUC_UNUSED GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area) {

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

   // std::cout << "shift is pressed " << shift_is_pressed << std::endl;

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

   do_drag_pan_gtk3(gl_area, drag_delta_x, drag_delta_y); // 20220613-PE no redraw here currently
   graphics_draw();
   double x = drag_begin_x + drag_delta_x;
   double y = drag_begin_y + drag_delta_y;
   set_mouse_previous_position(x, y);
   std::cout << "update_middle: " << x << " " << y << std::endl;
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
graphics_info_t::on_glarea_click(GtkGestureClick *controller,
                                 gint n_press,
                                 G_GNUC_UNUSED gdouble x,
                                 G_GNUC_UNUSED gdouble y,
                                 G_GNUC_UNUSED gpointer user_data) {

   SetMouseBegin(x,y);

   bool clicked = check_if_hud_bar_clicked(x,y);

   std::cout << "status for HUD bar clicked: " << clicked << " x " << x << " y " << y << std::endl;

   if (clicked) {
      // action occurs in above function
   } else {

      // std::cout << "n_press " << n_press << std::endl;
      // n_press can go up to 20, 30...
      //
      if (n_press == 2) { // otherwise triple clicking would toggle the label off, we don't want that.

         bool handled = false;

         if (in_moving_atoms_drag_atom_mode_flag) {
            if (last_restraints_size() > 0) {
               handled = check_if_moving_atom_pull(true); // passing was-a-double-click
            }
         }

         if (! handled) {
            bool intermediate_atoms_only_flag = false;
            pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
            int imol = naii.imol;
            if (naii.success) {
               molecules[imol].add_to_labelled_atom_list(naii.atom_index);
               add_picked_atom_info_to_status_bar(imol, naii.atom_index);
               graphics_draw();
            }
         }
      }

      if (n_press == 1) {

         GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(controller));
         std::cout << "debug:: on_glarea_click(); modifier: " << modifier << std::endl;
         if (modifier == 8) { // "option" key on Mac (ALT on PC is 24)
            bool intermediate_atoms_only_flag = false;
            pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
            if (naii.success) {
               setRotationCentre(naii.atom_index, naii.imol);
               add_picked_atom_info_to_status_bar(naii.imol, naii.atom_index);
            }

         } else { // not "option" modifier

            bool intermediate_atoms_only_flag = true;
            pick_info naii = atom_pick_gtk3(intermediate_atoms_only_flag);
            if (naii.success) {
               mmdb::Atom *at = moving_atoms_asc->atom_selection[naii.atom_index];
               std::cout << "debug:: in on_glarea_click() picked an intermediate atom " << coot::atom_spec_t(at) << std::endl;
            }
         }
      }
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

   // return; // for now

   int imol_map = make_and_draw_map_with_refmac_params(mtz_fn.c_str(), "FWT", "PHWT", "", 0, 0, 1, "FGMP18", "SIGFGMP18", "FreeR_flag", 1);
   int imol_diff_map = make_and_draw_map(mtz_fn.c_str(), "DELFWT", "PHDELWT", "", 0, 1);


}

gboolean
graphics_info_t::on_glarea_key_controller_key_pressed(GtkEventControllerKey *controller,
                                                      guint                  keyval,
                                                      guint                  keycode,
                                                      guint                  modifiers) {


   gboolean handled = false;

   control_is_pressed = (modifiers & GDK_CONTROL_MASK);
   shift_is_pressed   = (modifiers & GDK_SHIFT_MASK);

   if (false)
      std::cout << "on_glarea_key_controller_key_pressed() control_is_pressed " << control_is_pressed
                << " shift_is_pressed " << shift_is_pressed << std::endl;

   if (false)
      std::cout << "keyval: " << keyval << std::endl;

   if (keyval == 113)
      load_tutorial_model_and_data_ec();

   keyboard_key_t kbk(keyval, control_is_pressed);
   add_key_to_history(kbk);

   bool found = false;
   std::map<keyboard_key_t, key_bindings_t>::const_iterator it = key_bindings_map.find(kbk);
   if (it != key_bindings_map.end()) {
     const key_bindings_t &kb = it->second;
     if (false)
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


void
graphics_info_t::on_glarea_scrolled(GtkEventControllerScroll *controller,
                                    double                    dx,
                                    double                    dy,
                                    gpointer                  user_data) {

   GdkModifierType modifier = gtk_event_controller_get_current_event_state(GTK_EVENT_CONTROLLER(controller));
   control_is_pressed = (modifier & GDK_CONTROL_MASK);
   shift_is_pressed = (modifier & GDK_SHIFT_MASK);

   bool handled = false;
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
      }
   }

   if (! handled) {
      if (shift_is_pressed) {
         graphics_info_t::scroll_zoom(dy);
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
}


void
graphics_info_t::on_glarea_motion(G_GNUC_UNUSED GtkEventControllerMotion* controller,
                                  gdouble x,
                                  gdouble y,
                                  G_GNUC_UNUSED gpointer user_data) {


   // So that I can change the highlighting for the moused-over HUD buttons.

   // We can't easily use mouse_x_m mouse_y because they are used by update_view_quaternion().

   // mouse_x = x;
   // mouse_y = y;

   mouse_current_x = x;
   mouse_current_y = y;

   // set_mouse_previous_position(x, y);

   check_if_hud_button_moused_over(x, y, false);

}



void
graphics_info_t::change_model_molecule_representation_mode(int up_or_down) {

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
            bool force_rebonding = false;
            molecules[imol].make_colour_by_chain_bonds(force_rebonding);
         }
         if (bond_type == coot::COLOUR_BY_CHAIN_BONDS) {
            // graphics_to_colour_by_molecule(imol);
            bool force_rebonding = false;
            molecules[imol].make_colour_by_molecule_bonds(force_rebonding);
         }
         if (bond_type == coot::COLOUR_BY_MOLECULE_BONDS) {
            // graphics_to_ca_representation(imol);
            bool force_rebonding = false;
            molecules[imol].ca_representation(force_rebonding);
         }
         if (bond_type == coot::CA_BONDS) {
            // graphics_to_ca_plus_ligands_representation(imol);
            bool force_rebonding = false;
            molecules[imol].ca_plus_ligands_representation(Geom_p(), force_rebonding);
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS) {
            // graphics_to_ca_plus_ligands_sec_struct_representation(imol);
            molecules[imol].ca_plus_ligands_sec_struct_representation(Geom_p());
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR) {
            // graphics_to_rainbow_representation(imol);
            molecules[imol].ca_plus_ligands_rainbow_representation(Geom_p());
         }
         if (bond_type == coot::COLOUR_BY_RAINBOW_BONDS) {
            // graphics_to_sec_struct_bonds_representation(imol);
            molecules[imol].bonds_sec_struct_representation();
         }
         if (bond_type == coot::BONDS_SEC_STRUCT_COLOUR) {
            // graphics_to_bonds_no_waters_representation(imol);
            molecules[imol].bonds_no_waters_representation();
         }
         if (bond_type == coot::BONDS_NO_WATERS) {
            //graphics_to_b_factor_cas_representation(imol);
            molecules[imol].b_factor_representation_as_cas();
         }
         if (bond_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR) {
            // graphics_to_b_factor_representation(imol);
            molecules[imol].b_factor_representation();
         }
         if (bond_type == coot::COLOUR_BY_B_FACTOR_BONDS) {
            // graphics_to_occupancy_representation(imol);
            molecules[imol].occupancy_representation();
         }
         if (bond_type == coot::COLOUR_BY_OCCUPANCY_BONDS) {
            // graphics_to_bonds_representation(imol);
            bool force_rebonding = false;
            molecules[imol].bond_representation(Geom_p(), force_rebonding);
         }
      }
   }
}
