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
         return FALSE;
      } else {
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

   gtk_widget_queue_draw(widget); // needed?             

   return TRUE;
}



#include "screendump-tga.hh"

void
on_glarea_realize(GtkGLArea *glarea) {

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   gtk_gl_area_make_current(glarea);
   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);
   GLenum err = glGetError();
   err = glGetError(); if (err) std::cout << "on_glarea_realize() A err " << err << std::endl;

   // GLX_SAMPLE_BUFFERS_ARB
   // https://www.khronos.org/registry/OpenGL/extensions/ARB/ARB_multisample.txt
   // gdk/x11/gdkglcontext-x11.c
   // glXGetConfig(dpy, &visual_list[0], GLX_SAMPLE_BUFFERS_ARB, &gl_info[i].num_multisample);

   // glEnable(GL_MULTISAMPLE); // seems not to work at the moment. Needs work on the GTK->OpenGL interface 


   graphics_info_t g;
   g.init_shaders();
   g.init_buffers();
   err = glGetError();
   std::cout << "on_glarea_realize() post init_shaders() err is " << err << std::endl;

   graphics_info_t::shader_for_screen.Use(); // needed?

   err = glGetError(); std::cout << "start on_glarea_realize() err is " << err << std::endl;

   unsigned int index_offset = 0;
   graphics_info_t::screen_framebuffer.init(w, h, index_offset, "screen/occlusion");
   err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                          << err << std::endl;
   index_offset = 1;
   graphics_info_t::blur_framebuffer.init(w,h, index_offset, "blur");
   err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                          << err << std::endl;

   setup_hud_text(w, h, graphics_info_t::shader_for_hud_text, false);
   setup_hud_text(w, h, graphics_info_t::shader_for_atom_labels, true);

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

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

   glEnable(GL_DEPTH_TEST);

   // At some stage, enable this.  Currently (I think) the winding on the atoms is the wrong
   // way around - replace with octaspheres and octahemispheres.
   // I have just now changed the winding on the "solid" map triangles - and now it looks
   // fine.
   // 
   // glEnable(GL_CULL_FACE); // if I enable this, then I get to see the back side
                              // of the atoms. It's a weird look.

   // Make antialised lines
   if (false) {
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_LINE_SMOOTH);
   }

   g.setup_lights();

   gtk_gl_area_attach_buffers(GTK_GL_AREA(g.glareas[0])); // needed?   
   g.particles.make_particles(g.n_particles, g.get_rotation_centre());
   g.mesh_for_particles.setup_instancing_buffers_for_particles(g.particles.size());

   err = glGetError();
   if (err) std::cout << "on_glarea_realize() --end-- with err " << err << std::endl;

   g.mesh_for_particles.set_name("mesh for particles");

   g.tmesh_for_labels.setup_camera_facing_quad(&g.shader_for_atom_labels);

   g.setup_pulse_identification(); // not needed I think

   g.setup_hud_geometry_bars();

   g.setup_rama_balls();

   g.setup_key_bindings();
   
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
   setup_hud_text(width, height, g.shader_for_hud_text, false);
   setup_hud_text(width, height, g.shader_for_atom_labels, true); // change the function name

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

   if (! handled)
      g.contour_level_scroll_scrollable_map(direction);
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
      pick_info nearest_atom_index_info = g.atom_pick_gtk3(false);
      double delta_x = g.GetMouseClickedX() - event->x;
      double delta_y = g.GetMouseClickedY() - event->y;
      if (std::abs(delta_x) < 10.0) {
         if (std::abs(delta_y) < 10.0) {
            if (nearest_atom_index_info.success == GL_TRUE) {
               g.setRotationCentre(nearest_atom_index_info.atom_index,
                                   nearest_atom_index_info.imol);
               g.add_picked_atom_info_to_status_bar(nearest_atom_index_info.imol,
                                                    nearest_atom_index_info.atom_index);
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

   auto mouse_view_rotate = [control_is_pressed] (GtkWidget *widget, int x_as_int, int y_as_int) {
                               graphics_info_t g;
                               if (control_is_pressed) {
                                  do_drag_pan_gtk3(widget);
                               } else {

                                  bool handled = false;
                                  if (! handled) {
                                     GtkAllocation allocation;
                                     gtk_widget_get_allocation(widget, &allocation);
                                     int w = allocation.width;
                                     int h = allocation.height;
                                     graphics_info_t::update_view_quaternion(w, h);
                                  }
                               }
                            };

   auto mouse_zoom = [] (double delta_x, double delta_y) {
                        // Zooming
                        double fx = 1.0 + delta_x/300.0;
                        double fy = 1.0 + delta_y/300.0;
                        if (fx > 0.0) graphics_info_t::zoom /= fx;
                        if (fy > 0.0) graphics_info_t::zoom /= fy;
                        if (false)
                           std::cout << "zooming with perspective_projection_flag "
                                     << graphics_info_t::perspective_projection_flag
                                     << " " << graphics_info_t::zoom << std::endl;
                        if (! graphics_info_t::perspective_projection_flag) {
                           // std::cout << "now zoom: " << g.zoom << std::endl;
                        } else {
                           // Move the eye towards the rotation centre (don't move the rotation centre)
                           if (fabs(delta_y) > fabs(delta_x))
                              delta_x = delta_y;
                           float sf = 1.0 - delta_x * 0.003;
                           graphics_info_t::eye_position.z *= sf;

                           { // own graphics_info_t function - c.f. adjust clipping
                              double  l = graphics_info_t::eye_position.z;
                              double zf = graphics_info_t::screen_z_far_perspective;
                              double zn = graphics_info_t::screen_z_near_perspective;

                              graphics_info_t::screen_z_near_perspective *= sf;
                              graphics_info_t::screen_z_far_perspective  *= sf;

                              float screen_z_near_perspective_limit = l * 0.95;
                              float screen_z_far_perspective_limit  = l * 1.05;
                              if (graphics_info_t::screen_z_near_perspective < 2.0)
                                 graphics_info_t::screen_z_near_perspective = 2.0;
                              if (graphics_info_t::screen_z_far_perspective > 1000.0)
                                 graphics_info_t::screen_z_far_perspective = 1000.0;

                              if (graphics_info_t::screen_z_near_perspective > screen_z_near_perspective_limit)
                                 graphics_info_t::screen_z_near_perspective = screen_z_near_perspective_limit;
                              if (graphics_info_t::screen_z_far_perspective < screen_z_far_perspective_limit)
                                 graphics_info_t::screen_z_far_perspective = screen_z_far_perspective_limit;
                              if (false)
                                 std::cout << "on_glarea_motion_notify(): debug l: " << l << " post-manip: "
                                           << graphics_info_t::screen_z_near_perspective << " "
                                           << graphics_info_t::screen_z_far_perspective << std::endl;
                           }
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
      // View Panning
      do_drag_pan_gtk3(widget);
   }

   if (event->state & GDK_BUTTON3_MASK) {

      if (event->state & GDK_BUTTON1_MASK) {
         double delta_x = event->x - g.GetMouseBeginX();
         double delta_y = event->y - g.GetMouseBeginY();
         mouse_zoom(delta_x, delta_y);
      } else {
         if (!shift_is_pressed) {
            mouse_view_rotate(widget, x_as_int, y_as_int);
         }
         if (shift_is_pressed) {
            double delta_x = event->x - g.GetMouseBeginX();
            double delta_y = event->y - g.GetMouseBeginY();
            mouse_zoom(delta_x, delta_y);
         }
      }
   }

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

   std::chrono::time_point<std::chrono::system_clock> tp_now = std::chrono::high_resolution_clock::now();
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
        std::cout << "key-binding for key: " << it->first.gdk_key << " : "
                  << it->first.ctrl_is_pressed << " " << kb.description << std::endl;
     kb.run();
     found =  true;
   }

   // fix the type here
   if (int(event->keyval) == graphics_info_t::update_go_to_atom_from_current_residue_key) {
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
