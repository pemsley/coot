
// there are the only two that we need from c-interface.h
#include "c-interface.h"
#include "cc-interface.hh"
#include "glib.h"

// fromm cc-interface.hh
void set_bond_smoothness_factor(unsigned int fac);

#include "matrix-utils.hh"
#include "rsr-functions.hh"
#include "graphics-info.h"

// can't be a lambda funtion because of capture issues

void keypad_translate_xyz(short int axis, short int direction) {

   graphics_info_t g;
   if (axis == 3) {
      coot::Cartesian v = screen_z_to_real_space_vector(graphics_info_t::glareas[0]);
      v *= 0.05 * float(direction);
      g.add_vector_to_RotationCentre(v);
   } else {
      gdouble x_diff, y_diff;
      x_diff = y_diff = 0;
      coot::CartesianPair vec_x_y = screen_x_to_real_space_vector(graphics_info_t::glareas[0]);
      if (axis == 1) x_diff = 1;
      if (axis == 2) y_diff = 1;
      g.add_to_RotationCentre(vec_x_y, x_diff * 0.1 * float(direction),
                              y_diff * 0.1 * float(direction));
      if (g.GetActiveMapDrag() == 1) {
         for (int ii=0; ii<g.n_molecules(); ii++) {
            g.molecules[ii].update_map(true); // to take account
            // of new rotation centre.
         }
      }
      for (int ii=0; ii<g.n_molecules(); ii++) {
         g.molecules[ii].update_symmetry();
      }
      g.graphics_draw();
   }
}



// static
void
graphics_info_t::print_key_bindings() {

   std::vector<std::pair<std::string, int>> v = {
      std::make_pair("0", 0x030),
      std::make_pair("1", 0x031),
      std::make_pair("2", 0x032),
      std::make_pair("3", 0x033),
      std::make_pair("4", 0x034),
      std::make_pair("5", 0x035),
      std::make_pair("6", 0x036),
      std::make_pair("7", 0x037),
      std::make_pair("8", 0x038),
      std::make_pair("9", 0x039),
      std::make_pair("colon",     0x03a),
      std::make_pair("semicolon", 0x03b),
      std::make_pair("less",      0x03c),
      std::make_pair("equal",     0x03d),
      std::make_pair("greater",   0x03e),
      std::make_pair("question",  0x03f),
      std::make_pair("minus",    GDK_KEY_minus),
      std::make_pair("plus",    GDK_KEY_plus),
      std::make_pair("at", 0x040),
      std::make_pair("A", 0x041),
      std::make_pair("B", 0x042),
      std::make_pair("C", 0x043),
      std::make_pair("D", 0x044),
      std::make_pair("E", 0x045),
      std::make_pair("F", 0x046),
      std::make_pair("G", 0x047),
      std::make_pair("H", 0x048),
      std::make_pair("I", 0x049),
      std::make_pair("J", 0x04a),
      std::make_pair("K", 0x04b),
      std::make_pair("L", 0x04c),
      std::make_pair("M", 0x04d),
      std::make_pair("N", 0x04e),
      std::make_pair("O", 0x04f),
      std::make_pair("P", 0x050),
      std::make_pair("Q", 0x051),
      std::make_pair("R", 0x052),
      std::make_pair("S", 0x053),
      std::make_pair("T", 0x054),
      std::make_pair("U", 0x055),
      std::make_pair("V", 0x056),
      std::make_pair("W", 0x057),
      std::make_pair("X", 0x058),
      std::make_pair("Y", 0x059),
      std::make_pair("Z", 0x05a),
      std::make_pair("bracketleft",  0x05b),
      std::make_pair("backslash",    0x05c),
      std::make_pair("bracketright", 0x05d),
      std::make_pair("asciicircum",  0x05e),
      std::make_pair("underscore",   0x05f),
      std::make_pair("grave",        0x060),
      std::make_pair("quoteleft",    0x060),
      std::make_pair("a", 0x061),
      std::make_pair("b", 0x062),
      std::make_pair("c", 0x063),
      std::make_pair("d", 0x064),
      std::make_pair("e", 0x065),
      std::make_pair("f", 0x066),
      std::make_pair("g", 0x067),
      std::make_pair("h", 0x068),
      std::make_pair("i", 0x069),
      std::make_pair("j", 0x06a),
      std::make_pair("k", 0x06b),
      std::make_pair("l", 0x06c),
      std::make_pair("m", 0x06d),
      std::make_pair("n", 0x06e),
      std::make_pair("o", 0x06f),
      std::make_pair("p", 0x070),
      std::make_pair("q", 0x071),
      std::make_pair("r", 0x072),
      std::make_pair("s", 0x073),
      std::make_pair("t", 0x074),
      std::make_pair("u", 0x075),
      std::make_pair("v", 0x076),
      std::make_pair("w", 0x077),
      std::make_pair("x", 0x078),
      std::make_pair("y", 0x079),
      std::make_pair("z", 0x07a),
      std::make_pair("KP_plus",  0xffab),
      std::make_pair("KP_minus", 0xffad),
      std::make_pair("KP_0", 0xffb0),
      std::make_pair("KP_1", 0xffb1),
      std::make_pair("KP_2", 0xffb2),
      std::make_pair("KP_3", 0xffb3),
      std::make_pair("KP_4", 0xffb4),
      std::make_pair("KP_5", 0xffb5),
      std::make_pair("KP_6", 0xffb6),
      std::make_pair("KP_7", 0xffb7),
      std::make_pair("KP_8", 0xffb8),
      std::make_pair("KP_9", 0xffb9),
      std::make_pair("F1",  0xffbe),
      std::make_pair("F2",  0xffbf),
      std::make_pair("F3",  0xffc0),
      std::make_pair("F4",  0xffc1),
      std::make_pair("F5",  0xffc2),
      std::make_pair("F6",  0xffc3),
      std::make_pair("F7" , 0xffc4),
      std::make_pair("F8",  0xffc5),
      std::make_pair("F9",  0xffc6),
      std::make_pair("F10", 0xffc7),
      std::make_pair("F11", 0xffc8),
      std::make_pair("L1",  0xffc8),
      std::make_pair("F12", 0xffc9),
      std::make_pair("L2",  0xffc9),
      std::make_pair("Escape",  GDK_KEY_Escape),
      std::make_pair("Return",  GDK_KEY_Return),
      std::make_pair("Left",    GDK_KEY_Left),
      std::make_pair("Right",   GDK_KEY_Right),
      std::make_pair("Up",      GDK_KEY_Up),
      std::make_pair("Down",    GDK_KEY_Down),
   };

   // std::map<keyboard_key_t, key_bindings_t> key_bindings_map;

   for (const auto &kb : key_bindings_map) {
      std::string key_string = "gdk-symbol-" + std::to_string(kb.first.gdk_key);
      for (unsigned int ii=0; ii<v.size(); ii++) {
         if (v[ii].second == kb.first.gdk_key) {
            key_string = v[ii].first;
            break;
         }
      }
      std::string ctrl_string = "    ";
      if (kb.first.ctrl_is_pressed) ctrl_string = "Ctrl";
      std::cout << "binding: "
                << ctrl_string << " " << std::setw(5) << key_string
                << "  ->  " << std::setw(20) << std::left << kb.second.description
                << " " << key_bindings_t::type_to_string(kb.second.type)
                << std::endl;
   }
}

void
graphics_info_t::setup_key_bindings() {

   graphics_info_t g;

   // if we are serious about user-defined key-bindings all of these functions should be thunks in the user API
   // (and returning gboolean).

   auto l1 = []() { graphics_info_t g; g.adjust_clipping(-0.1); return gboolean(TRUE); };
   auto l2 = []() { graphics_info_t g; g.adjust_clipping( 0.1); return gboolean(TRUE); };
   auto l5 = []() { graphics_info_t g; g.blob_under_pointer_to_screen_centre(); return gboolean(TRUE); };

   auto l6 = []() {

                if (do_tick_spin) {
                   std::cout << "removing tick spin flag" << std::endl;
                   do_tick_spin = false;
                } else {
                   std::cout << "adding tick spin flag A" << std::endl;
                   if (! tick_function_is_active()) {
                      std::cout << "adding tick spin flag B" << std::endl;
                      int spin_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
                      // this is not a good name if we are storing a generic tick function id.
                      idle_function_spin_rock_token = spin_tick_id;
                   }
                   do_tick_spin = true;
                }
                return gboolean(TRUE);
             };

   auto l7 = []() {
                graphics_info_t g;
                int imol_scroll = g.intelligent_get_scroll_wheel_map();
                if (graphics_info_t::is_valid_map_molecule(imol_scroll)) {
                   graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
                }
                if (graphics_info_t::glareas.size() > 0)
                   int contour_idle_token = g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
                g.set_density_level_string(imol_scroll, graphics_info_t::molecules[imol_scroll].contour_level);
                graphics_info_t::display_density_level_this_image = 1;
                return gboolean(TRUE);
             };

   auto l8 = []() {
                graphics_info_t g;
                int imol_scroll = g.intelligent_get_scroll_wheel_map();
                if (graphics_info_t::is_valid_map_molecule(imol_scroll)) {
                   graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;
                }
                if (graphics_info_t::glareas.size() > 0)
                   int contour_idle_token = g_idle_add(idle_contour_function, graphics_info_t::glareas[0]);
                g.set_density_level_string(imol_scroll, graphics_info_t::molecules[imol_scroll].contour_level);
                graphics_info_t::display_density_level_this_image = 1;
                return gboolean(TRUE);
             };

   auto l9 = [] () {
                update_go_to_atom_from_current_position();
                return gboolean(TRUE);
             };

   auto l10 = []() { graphics_info_t::zoom *= 0.9; return gboolean(TRUE); };

   auto l11 = []() { graphics_info_t::zoom *= 1.1; return gboolean(TRUE); };

   auto l12 = []() { graphics_info_t g; g.move_forwards(); return gboolean(TRUE); };

   auto l13 = []() { graphics_info_t g; g.move_backwards(); return gboolean(TRUE); };

   auto l13l = []() { graphics_info_t g; g.step_screen_left();   return gboolean(TRUE); };
   auto l13r = []() { graphics_info_t g; g.step_screen_right(); return gboolean(TRUE); };

   //    auto l14 = []() { safe_python_command("import ncs; ncs.skip_to_next_ncs_chain('forward')"); return gboolean(TRUE); };

   // auto l14 = []() { /* use l28 */ return gboolean(TRUE); };

   // auto l15 = []() { /* use l28 */ return gboolean(TRUE); };

   auto l16 = []() { graphics_info_t g; g.undo_last_move(); return gboolean(TRUE); };

   auto l18 = []() { graphics_info_t g; g.clear_hud_buttons(); g.accept_moving_atoms(); return gboolean(TRUE); };

   auto l18_space = []() {
                       graphics_info_t g;
                       if (g.hud_button_info.size()) {
                          g.clear_hud_buttons(); g.accept_moving_atoms();
                       } else {

                          // Move the view - don't click the button

                          // g.reorienting_next_residue_mode = false; // hack
                          bool reorienting = graphics_info_t::reorienting_next_residue_mode;
                          if (reorienting) {
                             if (graphics_info_t::shift_is_pressed) {
                                g.reorienting_next_residue(false); // backwards
                             } else {
                                g.reorienting_next_residue(true); // forwards
                             }
                          } else {
                             // old/standard simple translation
                             if (graphics_info_t::shift_is_pressed) {
                                g.intelligent_previous_atom_centring(g.go_to_atom_window);
                             } else {
                                g.intelligent_next_atom_centring(g.go_to_atom_window);
                             }
                          }
                       }
                       return gboolean(TRUE);
                    };

   auto l19 = []() {
                 graphics_info_t g;
                 if (g.moving_atoms_asc) {
                    g.clear_up_moving_atoms_wrapper(); g.clear_gl_rama_plot();
                 } else {
                    unfullscreen();
                 }
                 return gboolean(TRUE);
              };

   auto l20 = []() { graphics_info_t g; g.eigen_flip_active_residue(); return gboolean(TRUE); };

   auto l21 = []() { graphics_info_t g; g.try_label_unlabel_active_atom(); return gboolean(TRUE); };

   auto l22 = []() {
                  graphics_info_t g;
                  g.setup_draw_for_particles();
                  return gboolean(TRUE);
             };

   // boids
   auto l23 = [] () {
      graphics_info_t g;

      if (true) {
         if (! graphics_info_t::do_tick_boids)
            graphics_info_t::do_tick_boids = true;
         else
            graphics_info_t::do_tick_boids = false;

	 // add_a_tick();

         g.setup_draw_for_boids();

	 std::cout << "----- key press ------ do_tick_boids "
		   << graphics_info_t::do_tick_boids << std::endl;
	 bool state = tick_function_is_active();
	 std::cout << "l23: tick_function_is_active() returns " << state << std::endl;
	 if (true) // use state in future?
	    gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
      }
      return gboolean(TRUE);
   };

   auto l24 = [] () {
                 // using the C API
                 // do_add_terminal_residue(1); // waits for user click :-)
                 graphics_info_t g;
                 g.add_terminal_residue_using_active_atom();
                 return gboolean(TRUE);
      };

   auto l25 = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    int imol_map = g.imol_refinement_map;
                    if (residue_p) {
                       mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].fill_partial_residue(residue_spec, g.Geom_p(), imol_map);

                       // now refine that
                       int saved_state = g.refinement_immediate_replacement_flag;
                       g.refinement_immediate_replacement_flag = 1;
                       std::string alt_conf("");
                       std::vector<mmdb::Residue *> rs = { residue_p };
                       g.refine_residues_vec(imol, rs, alt_conf, mol);
                       g.conditionally_wait_for_refinement_to_finish();
                       g.accept_moving_atoms();
                       g.refinement_immediate_replacement_flag = saved_state;
                    }
                 }
                 return gboolean(TRUE);
              };


   auto l26 = [] () {
                 graphics_info_t g;
                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       coot::residue_spec_t residue_spec(residue_p);
                       g.molecules[imol].delete_residue_sidechain(residue_spec);
                    }
                 }
                 return gboolean(TRUE);
              };

   auto l28 = [] () {

                 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
                 if (aa_spec_pair.first) {
                    int imol = aa_spec_pair.second.first;
                    mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
                    mmdb::Residue *residue_p = at->GetResidue();
                    if (residue_p) {
                       std::string this_chain_id = residue_p->GetChainID();
                       coot::residue_spec_t residue_spec(residue_p);
                       std::vector<std::vector<std::string> > ghost_chains_sets = molecules[imol].ncs_ghost_chains();
                       unsigned int n_ghost_chain_sets = ghost_chains_sets.size();
                       for (unsigned int i=0; i<n_ghost_chain_sets; i++) {
                          const std::vector<std::string> &chain_ids = ghost_chains_sets[i];
                          if (std::find(chain_ids.begin(), chain_ids.end(), this_chain_id) != chain_ids.end()) {
                             unsigned int idx_next = 0;
                             for (unsigned int j=0; j<chain_ids.size(); j++) {
                                if (chain_ids[j] == this_chain_id) {
                                   idx_next = j + 1;
                                   if (idx_next == chain_ids.size())
                                      idx_next = 0;
                                   break;
                                }
                             }
                             std::string chain_id_next = chain_ids[idx_next];
                             clipper::Coord_orth current_position = coot::co(at);
                             bool forward_flag = true;
                             glm::mat4 quat_mat = glm::toMat4(view_quaternion);
                             clipper::Mat33<double> current_view_mat = glm_to_mat33(quat_mat);

                             if (molecules[imol].ncs_ghosts_have_rtops_p() == 0)
                                molecules[imol].fill_ghost_info(1, ncs_homology_level);

                             std::pair<bool, clipper::RTop_orth> new_ori =
                                molecules[imol].apply_ncs_to_view_orientation(current_view_mat,
                                                                              current_position,
                                                                              this_chain_id, chain_id_next,
                                                                              forward_flag);
                             if (new_ori.first) {

                                clipper::Coord_orth t(new_ori.second.trn());
                                set_rotation_centre(t);

				view_quaternion = matrix_to_quaternion(new_ori.second.rot());

                                graphics_info_t g;
                                g.update_things_on_move(); // not static
                             }
                             break;
                          }
                       }
                    } else {
                       std::cout << "ERROR:: no residue" << std::endl;
                    }
                 }
                 graphics_draw();
                 return gboolean(TRUE);
              };

   auto l29 = [] () {

      graphics_info_t g;
      auto tp_now = std::chrono::high_resolution_clock::now();
      int n_press = g.get_n_pressed_for_leftquote_tap(tp_now);
      // std::cout << "highlighting active residue " << n_press << std::endl;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
         int imol = pp.second.first;
         g.update_mesh_for_outline_of_active_residue(imol, pp.second.second, n_press);
         if (! tick_function_is_active()) {
            int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
         }
         outline_for_active_residue_frame_count = 40;
         do_tick_outline_for_active_residue = true;
      }
      return gboolean(TRUE);
   };

   auto l31 = [] () {
                 graphics_info_t g;
                 g.decrease_clipping_front();
                 return gboolean(TRUE);
              };

   auto l32 = [] () {
                 graphics_info_t g;
                 g.increase_clipping_front();
                 return gboolean(TRUE);
              };

   auto l33 = [] () {
                 graphics_info_t g;
                 g.decrease_clipping_back();
                 return gboolean(TRUE);
              };

   auto l34 = [] () {
                 graphics_info_t g;
                 g.increase_clipping_back();
                 return gboolean(TRUE);
              };

   auto l35 = [] () {
                 graphics_info_t g;
                 g.wrapped_create_display_control_window();
                 return gboolean(TRUE);
              };

   auto l36 = [] () {
                 graphics_info_t g;
                 g.triple_refine_auto_accept();
                 return gboolean(TRUE);
              };

   auto l37 = [] {
      graphics_info_t g;
      g.display_next_map(); // one at a time, all, none.
      return gboolean(TRUE);
   };

   auto l38 = [] () {
      graphics_info_t g;
      g.toggle_display_of_last_model();
      return gboolean(TRUE);
   };

   auto l40 = [] () {
      rsr_sphere_refine_plus();
      return gboolean(TRUE);
   };

   auto l40c = [] () {
      rsr_refine_chain();
      return gboolean(TRUE);
   };

   auto l41 = [] () {
      bool is_all_em   = true;
      bool is_all_xray = true;
      for (int ii=0; ii<n_molecules(); ii++) {
         if (is_valid_map_molecule(ii)) {
            if (molecules[ii].is_EM_map()) {
               is_all_xray = false;
            } else {
               is_all_em = false;
            }
         }
      }
      if (is_all_xray)
         box_radius_xray *= (1.0/1.15);
      if (is_all_em)
         box_radius_em   *= (1.0/1.15);

      // 20250531-PE as it used to be:
      if ((! is_all_em) && (! is_all_em)) {
         box_radius_xray *= (1.0/1.15);
         box_radius_em   *= (1.0/1.15);
      }

      // is there an "update maps" function?
      for (int ii=0; ii<n_molecules(); ii++) {
         if (is_valid_map_molecule(ii))
            molecules[ii].update_map(true);
      }
      return gboolean(TRUE);
   };

   auto l42 = [] () {
      box_radius_xray *= 1.15;
      box_radius_em *= 1.15;
      for (int ii=0; ii<n_molecules(); ii++) {
         if (is_valid_map_molecule(ii))
            molecules[ii].update_map(true);
      }
      return gboolean(TRUE);
   };

   auto l43 = [] () {
      bool done = false;
      int scroll_wheel_map_prev = scroll_wheel_map;
      for (int ii=0; ii<n_molecules(); ii++) {
         if (is_valid_map_molecule(ii)) {
            if (ii > scroll_wheel_map) {
               scroll_wheel_map = ii;
               done = true;
               break;
            }
         }
      }
      if (! done) {
         for (int ii=0; ii<n_molecules(); ii++) {
            if (is_valid_map_molecule(ii)) {
               scroll_wheel_map = ii;
               break;
            }
         }
      }
      if (scroll_wheel_map != scroll_wheel_map_prev) {
         // we need to update the Display Manager
         graphics_info_t g;
         g.set_scrollable_map(scroll_wheel_map); // calls activate_scroll_radio_button_in_display_manager()
      }
      return gboolean(TRUE);
   };

   auto l44 = [] () {
      graphics_info_t g;
      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            g.backrub_rotamer_intermediate_atoms();
         }
      } else {
         std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
         int imol = aa.first;
         if (is_valid_model_molecule(imol)) {
            std::string alt_conf = aa.second->altLoc;
            coot::residue_spec_t res_spec(coot::atom_spec_t(aa.second));
            g.auto_fit_rotamer_ng(imol, res_spec, alt_conf);
         }
      }
      return gboolean(TRUE);
   };

   auto l45 = [] () {

      graphics_info_t g;
      bool done = false;
      // I need to be consistent about checking for moving_atoms_asc or moving_atoms_asc->mol
      // being null to mean if moving atoms are being displayed.
      // init() does a `new` for moving_atoms_asc.
      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            g.pepflip_intermediate_atoms();
            done = true;
         }
      }

      if (! done) {
         std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
         int imol = pp.second.first;
         if (is_valid_model_molecule(imol)) {
            coot::atom_spec_t as(pp.second.second);
            g.pepflip(imol, as);
         }
      }
      return gboolean(TRUE);
   };

   auto l46 = [] {
      show_keyboard_mutate_frame();
      return gboolean(TRUE);
   };

   auto l47 = [] {
      undo_symmetry_view();
      return gboolean(TRUE);
   };

   // Note to self, Space and Shift Space are key *Release* functions

   std::vector<std::pair<keyboard_key_t, key_bindings_t> > kb_vec;
   // kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_d,      key_bindings_t(l1, "increase clipping")));
   kb_vec.push_back(std::make_pair(GDK_KEY_d, key_bindings_t(l13r, "step right")));
   kb_vec.push_back(std::make_pair(GDK_KEY_a, key_bindings_t(l13l, "step left")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_e,      key_bindings_t(l44, "Auto-fit Rotamer")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_f,      key_bindings_t(l2, "decrease clipping")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_g,      key_bindings_t(l5, "go to blob")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_h,      key_bindings_t(l36, "Triple Refine with Auto-accept")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_i,      key_bindings_t(l6, "spin")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_j,      key_bindings_t(l44, "Auto-fit Rotamer"))); // where it used to be
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_plus,   key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_equal,  key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_minus,  key_bindings_t(l7, "decrease contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_KP_Add,      key_bindings_t(l8, "increase contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_KP_Subtract, key_bindings_t(l7, "decrease contour level")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_p,      key_bindings_t(l9, "update go-to atom by position")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_q,      key_bindings_t(l45, "Pep-flip")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_n,      key_bindings_t(l10, "Zoom in")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_m,      key_bindings_t(l11, "Zoom out")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_w,      key_bindings_t(l12, "Move forward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_s,      key_bindings_t(l13, "Move backward")));
   // kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l14, "NCS Skip forward")));
   // kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_O,      key_bindings_t(l15, "NCS Skip backward")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_u,      key_bindings_t(l16, "Undo Move")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Return, key_bindings_t(l18, "Accept Moving Atoms")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Escape, key_bindings_t(l19, "Reject Moving Atoms")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_l,      key_bindings_t(l21, "Label/Unlabel Active Atom")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_b,      key_bindings_t(l23, "Murmuration")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_y,      key_bindings_t(l24, "Add Terminal Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_k,      key_bindings_t(l25, "Fill Partial Residue")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_K,      key_bindings_t(l26, "Delete Sidechain")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_o,      key_bindings_t(l28, "NCS Other Chain")));

   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_A,      key_bindings_t(l38, "Toggle Display of Last Model")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_E,      key_bindings_t(l40c, "Chain Refine")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_M,      key_bindings_t(l46, "Keyboard Mutate")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_R,      key_bindings_t(l40, "Sphere Refine")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_Q,      key_bindings_t(l37, "Display Next Map")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_V,      key_bindings_t(l47, "Undo Symmetry View")));

   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_space,  key_bindings_t(l18_space, "Accept Moving Atoms")));

   // clipping
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_1,      key_bindings_t(l31, "Clipping Front Expand")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_2,      key_bindings_t(l32, "Clipping Front Reduce")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_3,      key_bindings_t(l33, "Clipping Back Reduce")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_4,      key_bindings_t(l34, "Clipping Back Expand")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_F8,     key_bindings_t(l35, "Show Display Manager")));

   // map radius
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_bracketleft,  key_bindings_t(l41, "Decrease Map Radius")));
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_bracketright, key_bindings_t(l42, "Increase Map Radius")));

   // scroll-wheel map change
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_W, key_bindings_t(l43, "Change Scroll-wheel map")));

   // control
   // meh - ugly and almost useless. Try again.
   // kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_asciitilde, key_bindings_t(l29, "Highlight Active Residue")));
   // try backtick:
   kb_vec.push_back(std::pair<keyboard_key_t, key_bindings_t>(GDK_KEY_quoteleft, key_bindings_t(l29, "Highlight Active Residue")));

   // control keys

   auto lc_copy = [] () { graphics_info_t g; g.copy_active_atom_molecule(); return gboolean(TRUE); };
   key_bindings_t copy_mol_key_binding(lc_copy, "Copy Model Molecule");
   std::pair<keyboard_key_t, key_bindings_t> p_copy(keyboard_key_t(GDK_KEY_c, true), copy_mol_key_binding);
   kb_vec.push_back(p_copy);

   auto lc1 = []() { show_go_to_residue_keyboarding_mode_window(); return gboolean(TRUE); };
   key_bindings_t go_to_residue_key_binding(lc1, "Show Go To Residue Keyboarding Window");
   std::pair<keyboard_key_t, key_bindings_t> p1(keyboard_key_t(GDK_KEY_g, true), go_to_residue_key_binding);
   kb_vec.push_back(p1);

   auto lc2 = []() { graphics_info_t g; g.apply_undo(); return gboolean(TRUE); };
   key_bindings_t undo_key_binding(lc2, "Undo");
   std::pair<keyboard_key_t, key_bindings_t> p2(keyboard_key_t(GDK_KEY_z, true), undo_key_binding);
   kb_vec.push_back(p2);

   auto lc3 = []() { graphics_info_t g; g.apply_redo(); return gboolean(TRUE);};
   key_bindings_t redo_key_binding(lc3, "Redo");
   std::pair<keyboard_key_t, key_bindings_t> p3(keyboard_key_t(GDK_KEY_y, true), redo_key_binding);
   kb_vec.push_back(p3);

   auto lc_res_info = []() {
      // this blob was copied from residue_info_action() - it could be refactored
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
         int imol = pp.second.first;
         coot::residue_spec_t res_spec(pp.second.second);
         output_residue_info_dialog(imol, res_spec);
      }
      return gboolean(TRUE);
   };
   key_bindings_t residue_info_key_binding(lc_res_info, "Residue Info");
   std::pair<keyboard_key_t, key_bindings_t> p_res_info(keyboard_key_t(GDK_KEY_i, true), residue_info_key_binding);
   kb_vec.push_back(p_res_info);

   auto ldr = [] () {

      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
      if (aa_spec_pair.first) {
         int imol = aa_spec_pair.second.first;
         mmdb::Atom *at = molecules[imol].get_atom(aa_spec_pair.second.second);
         mmdb::Residue *residue_p = at->GetResidue();
         if (residue_p) {
            // for this to work I need to move setup_delete_item_pulse() into
            // graphics_info_t. Not today.
            g.setup_delete_item_pulse(residue_p);
            coot::residue_spec_t residue_spec(residue_p);
            g.molecules[imol].delete_residue(residue_spec);
         }
      }
      return gboolean(TRUE);
   };
   key_bindings_t delete_residue_key_binding(ldr, "Delete Residue");
   std::pair<keyboard_key_t, key_bindings_t> pdel(keyboard_key_t(GDK_KEY_d, true), delete_residue_key_binding);
   kb_vec.push_back(pdel);

   auto law = [] () {
      graphics_info_t g;
      g.place_typed_atom_at_pointer("Water");
      return gboolean(TRUE);
   };
   key_bindings_t add_water_key_binding(law, "Add Water");
   std::pair<keyboard_key_t, key_bindings_t> paw(keyboard_key_t(GDK_KEY_w, true), add_water_key_binding);
   kb_vec.push_back(paw);

   // 2025-10-03-PE Thanks for the reminder AAAAdragon.
   auto l_go_to_lig = [] {
      go_to_ligand();
      return gboolean(TRUE);
   };
   key_bindings_t go_to_ligand_binding(l_go_to_lig, "Go To Ligand");
   std::pair<keyboard_key_t, key_bindings_t> pgl(keyboard_key_t(GDK_KEY_l, true), go_to_ligand_binding);
   kb_vec.push_back(pgl);

   // Direction is either +1 or -1 (in or out)
   //

   // ctrl left
   auto lc4 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed. No need to test it again here.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Left);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Left);
                 } else {
                    keypad_translate_xyz(1, 1);
                 }
                 return gboolean(TRUE);
              };

   // ctrl right
   auto lc5 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Right);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Right);
                 } else {
                    keypad_translate_xyz(1, -1);
                 }
                 return gboolean(TRUE);
              };

   // ctrl up
   auto lc6 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Up);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Up);
                 } else {
                    keypad_translate_xyz(2, 1);
                 }
                 return gboolean(TRUE);
              };
   // ctrl down
   auto lc7 = []() {
                 if (true) { // we don't get here unless Ctrl is pressed.
                    if (graphics_info_t::shift_is_pressed)
                       graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Down);
                    else
                       graphics_info_t::nudge_active_residue(GDK_KEY_Down);
                 } else {
                    keypad_translate_xyz(2, -1);
                 }
                 return gboolean(TRUE);
              };

   auto lc_qsa = [] () {
                    graphics_info_t g;
                    g.quick_save();
                    g.graphics_grab_focus();
                    return gboolean(TRUE);
                 };

   auto lc_smooth_bonds = [] () {
      graphics_info_t g;
      set_bond_smoothness_factor(2); // updates all molecules
      return gboolean(TRUE);
   };

   auto lc_toggle_validation_side_panel = [] () {
      GtkWidget* pane = widget_from_builder("main_window_ramchandran_and_validation_pane");
      if (pane) {
         if (gtk_widget_get_visible(pane) == TRUE) {
            gtk_widget_set_visible(pane, FALSE);
         } else {
            gtk_widget_set_visible(pane, TRUE);
         }
      }
      return gboolean(TRUE);
   };

   auto lc_toggle_alt_conf_view = [] () {
      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      const std::string &current_alt_conf = pp.second.second.alt_conf;
      if (pp.first) {
         int imol = pp.second.first;
         g.molecules[imol].alt_conf_view_next_alt_conf(current_alt_conf);
      }
      return gboolean(TRUE);
   };

   key_bindings_t ctrl_arrow_left_key_binding(lc4, "R/T Left");
   key_bindings_t ctrl_arrow_right_key_binding(lc5, "R/T Right");
   key_bindings_t ctrl_arrow_up_key_binding(lc6, "R/T Up");
   key_bindings_t ctrl_arrow_down_key_binding(lc7, "R/T Down");
   key_bindings_t ctrl_eigen_flip(l20, "Eigen-Flip");
   key_bindings_t ctrl_quick_save(lc_qsa, "Quick Save");
   key_bindings_t ctrl_toggle_panel(lc_toggle_validation_side_panel, "Toggle Validation Panel");
   key_bindings_t ctrl_toggle_alt_conf_view(lc_toggle_alt_conf_view, "Toggle Alt Conf View");
   key_bindings_t ctrl_bond_smoothness(lc_smooth_bonds, "Smooth bonds");

   std::pair<keyboard_key_t, key_bindings_t>  p4(keyboard_key_t(GDK_KEY_Left,  true), ctrl_arrow_left_key_binding);
   std::pair<keyboard_key_t, key_bindings_t>  p5(keyboard_key_t(GDK_KEY_Right, true), ctrl_arrow_right_key_binding);
   std::pair<keyboard_key_t, key_bindings_t>  p6(keyboard_key_t(GDK_KEY_Up,    true), ctrl_arrow_up_key_binding);
   std::pair<keyboard_key_t, key_bindings_t>  p7(keyboard_key_t(GDK_KEY_Down,  true), ctrl_arrow_down_key_binding);
   std::pair<keyboard_key_t, key_bindings_t>  p8(keyboard_key_t(GDK_KEY_e,    true), ctrl_eigen_flip);
   std::pair<keyboard_key_t, key_bindings_t>  p9(keyboard_key_t(GDK_KEY_s,    true), ctrl_quick_save);
   std::pair<keyboard_key_t, key_bindings_t> p10(keyboard_key_t(GDK_KEY_b,    true), ctrl_toggle_panel);
   std::pair<keyboard_key_t, key_bindings_t> p11(keyboard_key_t(GDK_KEY_a,    true), ctrl_toggle_alt_conf_view);
   std::pair<keyboard_key_t, key_bindings_t> p12(keyboard_key_t(GDK_KEY_f,    true), ctrl_bond_smoothness);

   kb_vec.push_back(p4);
   kb_vec.push_back(p5);
   kb_vec.push_back(p6);
   kb_vec.push_back(p7);
   kb_vec.push_back(p8);
   kb_vec.push_back(p9);
   kb_vec.push_back(p10);
   kb_vec.push_back(p11);
   kb_vec.push_back(p12);

   std::vector<std::pair<keyboard_key_t, key_bindings_t> >::const_iterator it;
   for (it=kb_vec.begin(); it!=kb_vec.end(); ++it)
     g.key_bindings_map[it->first] = it->second;

}

