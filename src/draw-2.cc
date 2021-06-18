#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "c-interface.h" // for update_go_to_atom_from_current_position()
#include "globjects.h"
#include "graphics-info.h"

#include "draw.hh"
#include "draw-2.hh"
#include "framebuffer.hh"

#include "text-rendering-utils.hh"
#include "cc-interface-scripting.hh"
#include "cylinder-with-rotation-translation.hh"

gint idle_contour_function(gpointer data);


glm::vec3
get_camera_up_direction(const glm::mat4 &mouse_quat_mat) {

   glm::vec4 z_p(0.0f, 1.0f, 0.0f, 1.0f);
   glm::vec4 r = z_p * mouse_quat_mat;
   glm::vec3 r3(r);
   return r3;
}


gboolean
glarea_tick_func(GtkWidget *widget,
                 GdkFrameClock *frame_clock,
                 gpointer data) {

   if (graphics_info_t::do_tick_particles) {
      if (graphics_info_t::particles.empty()) {
         graphics_info_t::do_tick_particles = false;
      } else {
         gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0])); // needed?
         // std::cout << "glarea_tick_func() calls update_particles() " << std::endl;
         graphics_info_t::particles.update_particles();
         graphics_info_t::mesh_for_particles.update_instancing_buffer_data_for_particles(graphics_info_t::particles);
      }
   }

   if (graphics_info_t::do_tick_spin) {
      float delta = 0.002;
      glm::vec3 EulerAngles(0, delta, 0);
      glm::quat quat_delta(EulerAngles);
      glm::quat normalized_quat_delta(glm::normalize(quat_delta));
      glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
      graphics_info_t::glm_quat = glm::normalize(product);
   }
   
   if (graphics_info_t::do_tick_boids) {
      graphics_info_t::boids.update();
      std::vector<glm::mat4> mats(graphics_info_t::boids.size());
      for (unsigned int ii=0; ii<graphics_info_t::boids.size(); ii++)
         mats[ii] = graphics_info_t::boids[ii].make_mat();
      graphics_info_t::mesh_for_boids.update_instancing_buffer_data(mats);
   }

   if (graphics_info_t::do_tick_hydrogen_bonds_mesh) {
      if (graphics_info_t::mesh_for_hydrogen_bonds.draw_this_mesh) {
         std::chrono::time_point<std::chrono::high_resolution_clock> tp_now =
            std::chrono::high_resolution_clock::now();
         std::chrono::time_point<std::chrono::high_resolution_clock> tp_prev =
            graphics_info_t::tick_hydrogen_bond_mesh_t_previous;
         auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - tp_prev);
         float theta = 0.002 * delta.count();
         // std::cout << "delta from time " << delta.count() << " theta " << theta << std::endl;
         std::vector<glm::mat4> mats;
         for (unsigned int i=0; i<graphics_info_t::hydrogen_bonds_atom_position_pairs.size(); i++) {
            const std::pair<glm::vec3, glm::vec3> &p = graphics_info_t::hydrogen_bonds_atom_position_pairs[i];
            mats.push_back(Mesh::make_hydrogen_bond_cylinder_orientation(p.first, p.second, theta));
         }
         // std::cout << "mesh_for_hydrogen_bonds.update_instancing_buffer_data(mats) " << mats.size() << std::endl;
         graphics_info_t::mesh_for_hydrogen_bonds.update_instancing_buffer_data(mats);
      }
   }

   if (false)
      std::cout << "### in the glarea_tick_func() "
                << graphics_info_t::do_tick_happy_face_residue_markers << "  "
                << graphics_info_t::draw_count_for_happy_face_residue_markers << std::endl;
   if (graphics_info_t::do_tick_happy_face_residue_markers) {
      if (graphics_info_t::tmesh_for_happy_face_residues_markers.draw_this_mesh) {
         graphics_info_t::draw_count_for_happy_face_residue_markers += 1;
         graphics_info_t g;
         if (graphics_info_t::draw_count_for_happy_face_residue_markers >= g.draw_count_max_for_happy_face_residue_markers) {
            graphics_info_t::do_tick_happy_face_residue_markers = false;
            graphics_info_t::draw_count_for_happy_face_residue_markers = 0;

            graphics_info_t::tmesh_for_happy_face_residues_markers.draw_this_mesh =  false;
         }

         // repeating code in setup_draw_for_happy_face_residue_markers()

         const std::vector<glm::vec3> &positions = graphics_info_t::happy_face_residue_marker_starting_positions;
         glm::vec3 up_uv = g.get_screen_y_uv();
         unsigned int draw_count = g.draw_count_for_happy_face_residue_markers;
         unsigned int draw_count_max = g.draw_count_max_for_happy_face_residue_markers;
         g.tmesh_for_happy_face_residues_markers.update_instancing_buffer_data_for_happy_faces(positions, draw_count, draw_count_max, up_uv);
      }
   }

   gtk_widget_queue_draw(widget); // needed?

   if (graphics_info_t::do_tick_particles ||
       graphics_info_t::do_tick_spin      ||
       graphics_info_t::do_tick_boids     ||
       graphics_info_t::do_tick_hydrogen_bonds_mesh ||
       graphics_info_t::do_tick_happy_face_residue_markers)
      return gboolean(TRUE);
   else
      return gboolean(FALSE);
}


#include <gdk/gdk.h>
#include "screendump-tga.hh"

void
on_glarea_realize(GtkGLArea *glarea) {

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   std::cout << "debug:: on_glarea_realize() about to make_current()" << std::endl;
   gtk_gl_area_make_current(glarea);
   GLenum err = glGetError();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() A err " << err << std::endl;
   if (gtk_gl_area_get_error(glarea) != NULL) {
      std::cout << "OOPS:: on_glarea_realize() error on gtk_gl_area_make_current()" << std::endl;
      return;
   }

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

   // GLX_SAMPLE_BUFFERS_ARB
   // https://www.khronos.org/registry/OpenGL/extensions/ARB/ARB_multisample.txt
   // gdk/x11/gdkglcontext-x11.c
   // glXGetConfig(dpy, &visual_list[0], GLX_SAMPLE_BUFFERS_ARB, &gl_info[i].num_multisample);

   // glEnable(GL_MULTISAMPLE); // seems not to work at the moment. Needs work on the GTK->OpenGL interface

   const char *s1 = reinterpret_cast<const char *>(glGetString(GL_VERSION));
   const char *s2 = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
   const char *s3 = reinterpret_cast<const char *>(glGetString(GL_RENDERER));
   const char *s4 = reinterpret_cast<const char *>(glGetString(GL_VENDOR));
   if (s1 && s2 && s3 && s4) {
      std::string ss1(s1);
      std::string ss2(s2);
      std::string ss3(s3);
      std::string ss4(s4);
      std::cout << "INFO:: GL Version:                  " << ss1 << std::endl;
      std::cout << "INFO:: GL Shading Language Version: " << ss2 << std::endl;
      std::cout << "INFO:: GL Renderer:                 " << ss3 << std::endl;
      std::cout << "INFO:: GL Vendor:                   " << ss4 << std::endl;
   } else {
      std::cout << "error:: on_glarea_realize() null from glGetString()" << std::endl;
   }

   graphics_info_t g;
   bool status = g.init_shaders();

   if (status == true) {
      // happy path

      g.init_buffers();
      err = glGetError();
      if (err) std::cout << "error:: on_glarea_realize() post init_shaders() err is " << err << std::endl;

      graphics_info_t::shader_for_screen.Use(); // needed?

      err = glGetError();
      if (err) std::cout << "error:: start on_glarea_realize() err is " << err << std::endl;

      if (graphics_info_t::use_framebuffers) {
         unsigned int index_offset = 0;
         graphics_info_t::screen_framebuffer.init(w, h, index_offset, "screen/occlusion");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 1;
         graphics_info_t::blur_framebuffer.init(w,h, index_offset, "blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                                << err << std::endl;

         graphics_info_t::shader_for_screen.Use();
         err = glGetError(); if (err) std::cout << "on_glarea_realize() B screen framebuffer err " << err << std::endl;
         graphics_info_t::shader_for_screen.set_int_for_uniform("screenTexture", 0);
         err = glGetError(); if (err) std::cout << "on_glarea_realize() C screen framebuffer err " << err << std::endl;
         graphics_info_t::shader_for_screen.set_int_for_uniform("screenDepth", 1);
         err = glGetError(); if (err) std::cout << "on_glarea_realize() D screen framebuffer err " << err << std::endl;

         graphics_info_t::shader_for_blur.Use();
         err = glGetError(); if (err) std::cout << "on_glarea_realize() blur shader-framebuffer B err " << err << std::endl;
         graphics_info_t::shader_for_blur.set_int_for_uniform("screenTexture", 0);
         err = glGetError(); if (err) std::cout << "on_glarea_realize() blur C shader-framebuffer err " << err << std::endl;
         graphics_info_t::shader_for_blur.set_int_for_uniform("screenDepth", 1);
         err = glGetError(); if (err) std::cout << "on_glarea_realize() blur D shader-framebuffer err " << err << std::endl;
      }

      setup_hud_text(w, h, graphics_info_t::shader_for_hud_text, false);
      setup_hud_text(w, h, graphics_info_t::shader_for_atom_labels, true);

      gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

      glEnable(GL_DEPTH_TEST);

      // At some stage, enable this.  Currently (I think) the winding on the atoms is the wrong
      // way around - replace with octaspheres and octahemispheres.
      // I have just now changed the winding on the "solid" map triangles - and now it looks
      // fine.
      // 
      // glEnable(GL_CULL_FACE); // if I enable this, then I get to see the back side
      // of the atoms. It's a weird look.

      // Make antialised lines - not in this
      if (false) {
         glEnable (GL_BLEND);
         glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
         glEnable(GL_LINE_SMOOTH);
      }

      g.setup_lights();

      float x_scale = 4.4;  // what are these numbers!?
      float y_scale = 1.2;
      g.tmesh_for_labels.setup_camera_facing_quad(&g.shader_for_atom_labels, x_scale, y_scale);

      g.setup_hud_geometry_bars();

      g.setup_rama_balls();

      g.setup_key_bindings();

      g.setup_draw_for_happy_face_residue_markers_init();

      err = glGetError();
      if (err) std::cout << "################ GL ERROR on_glarea_realize() --end-- with err "
                         << err << std::endl;

      GdkGLContext *context = gtk_gl_area_get_context(GTK_GL_AREA(glarea));
      gboolean legacy_flag = gdk_gl_context_is_legacy(context);
      std::cout << "INFO:: gdk_gl_context_is_legacy() returns " << legacy_flag << std::endl;
      std::cout << "------------------------ realize() done " << std::endl;

   } else {
      std::cout << "ERROR:: Shader compilation failed " << std::endl;
      exit(1);
   }

}

gboolean
on_glarea_render(GtkGLArea *glarea) {

   return graphics_info_t::render(false); // not to screendump framebuffer.

}


void
on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   graphics_info_t g;
   g.graphics_x_size = width;
   g.graphics_y_size = height;

   // why do I need to do this?
   setup_hud_text(width, height, g.shader_for_hud_text, false);
   setup_hud_text(width, height, g.shader_for_atom_labels, true); // change the function name

   g.setup_hud_geometry_bars(); // because they depend on the aspect ratio - but can't that be
                                // passed as a uniform?

   // std::cout << "INFO:: Reset frame buffers " << width << "x" << height << std::endl;
   g.reset_frame_buffers(width, height);
}


gboolean
on_glarea_scroll(GtkWidget *widget, GdkEventScroll *event) {


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
   return TRUE;
}


gboolean
on_glarea_button_press(GtkWidget *widget, GdkEventButton *event) {

   // std::cout << "button press!" << std::endl;
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

   GdkModifierType state;

   // if (true) { // check here for left-mouse click
   // if (event->state & GDK_BUTTON1_PRESS) {
   if (true) { // check here for left-mouse click

      bool handled = false;

      if (false)
         std::cout << "click event: " << event->x << " " << event->y << " "
                   << x_as_int << " " << y_as_int << std::endl;

      // first thing to test is the HUD bar
      handled = g.check_if_hud_bar_clicked(event->x, event->y);

      if (! handled) {
         // implicit type cast
         handled = g.check_if_moving_atom_pull(was_a_double_click);

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
               }
            }
         }
      }
   }


   g.check_if_in_range_defines(event, mask);
   return TRUE;
}

gboolean
on_glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

   if (graphics_info_t::in_moving_atoms_drag_atom_mode_flag) {
      graphics_info_t g;
      g.unset_moving_atoms_currently_dragged_atom_index();
      g.do_post_drag_refinement_maybe();
      graphics_info_t::in_moving_atoms_drag_atom_mode_flag = 0;
   }

   if (event->state & GDK_BUTTON2_MASK) {
      graphics_info_t g;
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
                     std::cout << "sap success" << std::endl;
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
   return TRUE;
}

void
do_drag_pan_gtk3(GtkWidget *widget) {

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
   g.update_maps();
   if (graphics_info_t::glareas.size() > 0)
      int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);
}

gboolean
on_glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   int r = 0;
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

   std::pair<bool, mmdb::Atom *> handled_pair = g.check_if_moused_over_hud_bar(event->x, event->y);

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
      return TRUE;
   } else {
      g.draw_hud_tooltip_flag = false;
   }


   auto mouse_view_rotate = [control_is_pressed] (GtkWidget *widget) {
                               if (control_is_pressed) {
                                  do_drag_pan_gtk3(widget);
                               } else {
                                  GtkAllocation allocation;
                                  gtk_widget_get_allocation(widget, &allocation);
                                  int w = allocation.width;
                                  int h = allocation.height;
                                  graphics_info_t::update_view_quaternion(w, h);
                               }
                            };

   if (event->state & GDK_BUTTON1_MASK) {

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

   if (event->state & GDK_BUTTON2_MASK) {
      if (shift_is_pressed) {
         // mouse_view_rotate(widget);
         std::cout << "shift middle mouse - what to do here?" << std::endl;
      } else {
         do_drag_pan_gtk3(widget);          // View Panning
      }
   }

   if (event->state & GDK_BUTTON3_MASK) {
      double delta_x = event->x - g.GetMouseBeginX();
      double delta_y = event->y - g.GetMouseBeginY();
      if (event->state & GDK_BUTTON1_MASK) {
         // chording
         g.mouse_zoom(delta_x, delta_y);
      } else {
         if (! shift_is_pressed) {
            mouse_view_rotate(widget);
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
   return TRUE;
}

gint
view_spin_func(gpointer data) {

   float delta = 0.002;
   glm::vec3 EulerAngles(0, delta, 0);
   glm::quat quat_delta(EulerAngles);
   glm::quat normalized_quat_delta(glm::normalize(quat_delta));
   glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
   graphics_info_t::glm_quat = glm::normalize(product);
   graphics_info_t::graphics_draw(); // queue

   std::chrono::time_point<std::chrono::high_resolution_clock> tp_now = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed_seconds = tp_now - graphics_info_t::previous_frame_time_for_per_second_counter;
   if (elapsed_seconds.count() > 1.0) {
      float nf = graphics_info_t::frame_counter - graphics_info_t::frame_counter_at_last_display;
      std::cout << "INFO:: Time/frame: " << 1000 * elapsed_seconds.count()/nf << " milliseconds "
                << nf/elapsed_seconds.count() << " frames/second\n";
      graphics_info_t::previous_frame_time_for_per_second_counter = tp_now;
      graphics_info_t::frame_counter_at_last_display = graphics_info_t::frame_counter;
   }

   // now the stutter checker:
   std::chrono::duration<double> elapsed_seconds_fast = tp_now - graphics_info_t::previous_frame_time;
   if (elapsed_seconds_fast.count() > 0.03)
      std::cout << "INFO:: " << 1000 * elapsed_seconds_fast.count() << " milliseconds for that frame\n";
   graphics_info_t::previous_frame_time = tp_now;

   // kludge/race condition?
   if (graphics_info_t::idle_function_spin_rock_token == -1)
      return FALSE;
   else
      return TRUE;
}



gboolean
on_glarea_key_press_notify(GtkWidget *widget, GdkEventKey *event) {

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

   int kv = event->keyval;
   if (kv == graphics_info_t::update_go_to_atom_from_current_residue_key) {
      update_go_to_atom_from_current_position();
      handled = TRUE;
   }

   if (! found)
      if (! handled)
         std::cout << "on_glarea_key_press_notify() key not found in map: " << event->keyval << std::endl;

   graphics_info_t::graphics_draw(); // queue

   return handled;

}

gboolean
on_glarea_key_release_notify(GtkWidget *widget, GdkEventKey *event) {

   graphics_info_t g;

   // We need to check the GDK_KEY_Shift_R also, I guess. Not clear
   // to me how to do that now. Fix later.
   g.shift_is_pressed = false;
   if (event->state & GDK_SHIFT_MASK) g.shift_is_pressed = true;
   if (event->keyval == GDK_KEY_Shift_L) g.shift_is_pressed = true;

   if (event->keyval == GDK_KEY_space) {
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
   return TRUE;
}

void
my_glarea_add_signals_and_events(GtkWidget *glarea) {

   gtk_widget_add_events(glarea, GDK_SCROLL_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_PRESS_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_RELEASE_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_POINTER_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MASK);
   gtk_widget_add_events(glarea, GDK_KEY_PRESS_MASK);

   // key presses for the glarea:

   gtk_widget_set_can_focus(glarea, TRUE);
   gtk_widget_grab_focus(glarea);

   g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
   g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
   g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);
   g_signal_connect(glarea, "scroll-event",          G_CALLBACK(on_glarea_scroll),             NULL);
   g_signal_connect(glarea, "button-press-event",    G_CALLBACK(on_glarea_button_press),       NULL);
   g_signal_connect(glarea, "button-release-event",  G_CALLBACK(on_glarea_button_release),     NULL);
   g_signal_connect(glarea, "motion-notify-event",   G_CALLBACK(on_glarea_motion_notify),      NULL);
   g_signal_connect(glarea, "key-press-event",       G_CALLBACK(on_glarea_key_press_notify),   NULL);
   g_signal_connect(glarea, "key-release-event",     G_CALLBACK(on_glarea_key_release_notify), NULL);

}
