/*
 * src/new-startup.cc
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

#include <cstddef>
#include <iostream>
#include <string>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include <clipper/core/test_core.h>
#include <clipper/contrib/test_contrib.h>

#include "glib-object.h"
#include "glib.h"
#include "glibconfig.h"
#include "utils/xdg-base.hh"

#include "graphics-info.h"
#include "create-menu-item-actions.hh"
#include "setup-gui-components.hh"
#include "coot-setup-python.hh"
#include "utils/coot-utils.hh"
#include "command-line.hh"
#include "c-interface-preferences.h"
#include "boot-python.hh"
#include "layla/layla_embedded.hpp"

#include "testing.hh" // for test_internal();

#include "utils/logging.hh"
#include "widget-from-builder.hh"
std::string git_commit(); // use a header?

extern logging logger;

void print_opengl_info();

void init_framebuffers(GtkWidget *glarea) {

   // put this into graphics-info I suppose.

   // std::cout << "DEBUG:: use_framebuffers: " << graphics_info_t::use_framebuffers << std::endl;

   // std::cout << "----- start init_framebuffers() ----" << std::endl;

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

   GLenum err = glGetError();
   if (err)
      std::cout << "ERROR:: init_framebuffers() --- start --- err is " << err << std::endl;

   if (graphics_info_t::use_framebuffers) {
      unsigned int index_offset = 0;
      graphics_info_t::screen_framebuffer.init(w, h, index_offset, "screen/occlusion");
      err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                             << err << std::endl;
      // index_offset = 1;
      graphics_info_t::blur_y_framebuffer.init(w, h, index_offset, "blur-y");
      err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_y_framebuffer init() err is "
                                             << err << std::endl;
      // index_offset = 2;
      graphics_info_t::blur_x_framebuffer.init(w, h, index_offset, "blur-x");
      err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_x_framebuffer init() err is "
                                             << err << std::endl;
      // index_offset = 3;
      graphics_info_t::combine_textures_using_depth_framebuffer.init(w, h, index_offset, "new-blur");
      err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_combine framebuffer init() err is "
                                             << err << std::endl;
      // index_offset = 4;
      graphics_info_t::blur_framebuffer.init(w, h, index_offset, "blur");
      err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                             << err << std::endl;
      err = glGetError();
      if (err)
         std::cout << "ERROR:: init_framebuffers() --- done --- err is " << err << std::endl;
   }

   // std::cout << "----- done init_framebuffers() ----" << std::endl;
}


#include "text-rendering-utils.hh"
#include "stringify-error-code.hh"
// from c-inteerface.cc
extern "C" void run_command_line_scripts();


void
new_startup_realize(GtkWidget *gl_area) {

   GdkDisplay *display = gdk_display_get_default();
   GListModel *lm = gdk_display_get_monitors(display);

   guint n_items = g_list_model_get_n_items(lm);
   if (n_items > 0) {
      for (unsigned int i=0; i<n_items; i++) {
         gpointer item = g_list_model_get_item(lm, i);
         int imon = i + 1;
         GdkMonitor *monitor = GDK_MONITOR(item);
         const char *monitor_description = gdk_monitor_get_description(monitor);
         const char *monitor_connection  = gdk_monitor_get_connector(monitor);
         if (monitor_description)
            // std::cout << "INFO:: monitor " << imon << " description " << monitor_description << std::endl;
            logger.log(log_t::INFO, "monitor", std::to_string(imon), "description", monitor_description);
         else
            // std::cout << "INFO:: monitor " << imon << " no description " << std::endl;
            logger.log(log_t::INFO, "monitor", std::to_string(imon), "no description");
         if (monitor_connection)
            // std::cout << "INFO:: monitor " << imon << " connection "  << monitor_connection  << std::endl;
            logger.log(log_t::INFO, "monitor ", std::to_string(imon), " connection ", monitor_connection);
         int monitor_refresh_rate = gdk_monitor_get_refresh_rate(monitor);
         // std::cout << "INFO:: monitor " << imon << " refresh rate " << monitor_refresh_rate << " mHz"  << std::endl;
         logger.log(log_t::INFO, "monitor", imon, "refresh rate", monitor_refresh_rate, "mHz");
         int monitor_scale_factor = gdk_monitor_get_scale_factor(monitor);
         // std::cout << "INFO:: monitor " << imon << " scale_factor " << monitor_scale_factor << std::endl;
         logger.log(log_t::INFO, "monitor", std::to_string(imon), "scale_factor", monitor_scale_factor);

#if GTK_MINOR_VERSION >= 14
         double monitor_scale = gdk_monitor_get_scale(monitor);
         // std::cout << "INFO:: monitor " << imon << " scale " << monitor_scale << std::endl;
         logger.log(log_t::INFO, "monitor", std::to_string(imon), "scale", monitor_scale);
#endif
      }
   }

   gtk_gl_area_make_current(GTK_GL_AREA (gl_area));

   GError* error = gtk_gl_area_get_error(GTK_GL_AREA (gl_area));
   if (error != NULL) {
      std::cout << "WARNING:: new_startup_realize() gtk_gl_area_get_error() returned an error: " << std::endl;
      std::cout << error->message << std::endl;
      return;
   }

   GdkGLContext *context = gtk_gl_area_get_context(GTK_GL_AREA(gl_area));

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(gl_area), TRUE);

   print_opengl_info();


   int w = 500;
   int h = 500;

   graphics_info_t g;
   // init_framebuffers(gl_area); // Hmm - I don't know what this does compared to below.
   g.init_buffers();
   g.init_shaders();
   g.setup_lights();
   // 20241001-PE glDrawBuffer() in init_framebuffers() barfs with GL ES.
   if (!g.graphics_is_gl_es) {
      g.init_framebuffers(w, h);
      g.init_joey_ssao_stuff(w, h);
   }

   float x_scale = 4.4;  // what are these numbers!?
   float y_scale = 1.2;
   x_scale = 1.002;
   y_scale = 1.002;
   g.tmesh_for_labels.setup_camera_facing_quad(x_scale, y_scale, 0.0, 0.0);
   g.setup_hud_geometry_bars();
   g.setup_hud_buttons();
   g.setup_rama_balls();
   g.setup_key_bindings();
   float double_rama_size = 0.8; // scaled by 0.5 in the gl-rama draw call.
   g.gl_rama_plot.setup_buffers(double_rama_size); // rama relative size, put it into graphics_info_t
   // and allow it to be set in the API
   g.setup_draw_for_happy_face_residue_markers_init();
   g.setup_draw_for_bad_nbc_atom_pair_markers(); // angry diego
   g.setup_draw_for_bad_nbc_atom_pair_dashed_line();
   g.setup_draw_for_chiral_volume_outlier_markers();
   g.setup_draw_for_anchored_atom_markers_init();
   g.setup_draw_for_unhappy_atom_markers();
   g.setup_lines_mesh_for_proportional_editing();
   g.lines_mesh_for_hud_lines.set_name("lines mesh for fps graph");
   unsigned int frame_time_history_list_max_n_elements = 500;
   // +40 for base and grid lines
   std::vector<s_generic_vertex> empty_vertices(frame_time_history_list_max_n_elements + 40);
   std::vector<unsigned int> empty_indices(1500, 0); // or some number
   g.lines_mesh_for_hud_lines.setup_vertices_and_indices(empty_vertices, empty_indices);

   setup_hud_text(w, h, graphics_info_t::shader_for_hud_text, false);
   setup_hud_text(w, h, graphics_info_t::shader_for_atom_labels, true);
   g.tmesh_for_hud_refinement_dialog_arrow = HUDTextureMesh("HUD tmesh for refinement dialog arrow");
   g.tmesh_for_hud_refinement_dialog_arrow.setup_quad();
   g.texture_for_hud_refinement_dialog_arrow             = Texture("refinement-dialog-arrrow.png", Texture::DIFFUSE);
   g.texture_for_hud_refinement_dialog_arrow_highlighted = Texture("refinement-dialog-arrrow-highlighted.png", Texture::DIFFUSE);

   g.tmesh_for_shadow_map.setup_quad();

   Material material;
   GLenum err = glGetError();
   if (err)
      std::cout << "ERROR:: new_startup_realize() pos-D err is " << stringify_error_code(err)
                << std::endl;
   // g.attach_buffers();
   err = glGetError();
   if (err)
      std::cout << "ERROR:: new_startup_realize() pos-E post attach_buffers() err is "
                << stringify_error_code(err) << std::endl;
   g.mesh_for_extra_distance_restraints.setup_extra_distance_restraint_cylinder(material); // init

   // scale the gizmo to the object being translated
   float scale_factor = 22.2;
   g.translation_gizmo.scale(scale_factor);
   g.setup_draw_for_translation_gizmo();

   g.setup_key_bindings();

   err = glGetError();
   if (err)
      std::cout << "ERROR:: new_startup_realize() --end-- err is " << stringify_error_code(err)
                << std::endl;

   auto run_command_line_scripts_callback = +[] (gpointer user_data) {
      run_command_line_scripts();
      return G_SOURCE_REMOVE;
   };
   g_idle_add(run_command_line_scripts_callback, nullptr);

   // Hmm! - causes weird graphics problems
   // setup_python(0, NULL); // needs to called after GTK has started - because it depends on gtk.
                             // 20220629-PE not at the moment though - I removed the gobject parts from the code path


}


void
new_startup_unrealize(GtkWidget *widget) {

   gtk_gl_area_make_current (GTK_GL_AREA (widget));
   if (gtk_gl_area_get_error (GTK_GL_AREA (widget)) != NULL)
      return;

}


gboolean
new_startup_on_glarea_render(GtkGLArea *glarea) {

   // std::cout << "DEBUG: new_startup_on_glarea_render()!" << std::endl;
   bool screen_dump_frame_buffer = false;
   return graphics_info_t::render(screen_dump_frame_buffer);
}


#include "c-interface.h" // for run_script()
void
new_startup_on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   graphics_info_t g;
   // for the GL widget, not the window.
   g.graphics_x_size = width;
   g.graphics_y_size = height;
   if (g.graphics_is_gl_es) {
      // don't touch the framebuffers
   } else {
       g.reset_frame_buffers(width, height); // currently makes the widget blank (not drawn)
       g.resize_framebuffers_textures_renderbuffers(width, height); // 20220131-PE added from crows merge
   }
   g.reset_hud_buttons_size_and_position();
   g.mouse_speed = static_cast<double>(width) / 900.0;

   if (false) {

      // 20220807-PE It seems that now the shaders have been compiled before the first window resize - good.
      //             Just leave this here for now.

      if (! g.shaders_have_been_compiled) {
         g.init_shaders();
      } else {
         std::cout << "in new_startup_on_glarea_resize() shaders have already been compiled!" << std::endl;
      }
   }

}

void
new_startup_on_glarea_enter(GtkGLArea *glarea) {

  std::cout << "enter!" << std::endl;
}

// void on_glarea_realize(GtkWidget *widget); // using this give linking problems.

GtkWidget *new_startup_create_glarea_widget() {

   GtkWidget *gl_area = gtk_gl_area_new();

#if (GTK_MAJOR_VERSION > 4) || (GTK_MAJOR_VERSION == 4 && GTK_MINOR_VERSION >= 12)
   // Disable OpenGL ES (if not OpenGL_ES set on the command line)
   if (! graphics_info_t::graphics_is_gl_es)
      gtk_gl_area_set_allowed_apis(GTK_GL_AREA(gl_area), GDK_GL_API_GL);
#endif

   g_signal_connect(gl_area, "realize",   G_CALLBACK(new_startup_realize),   NULL);
   g_signal_connect(gl_area, "unrealize", G_CALLBACK(new_startup_unrealize), NULL);
   g_signal_connect(gl_area, "render",    G_CALLBACK(new_startup_on_glarea_render),  NULL);
   g_signal_connect(gl_area, "resize",    G_CALLBACK(new_startup_on_glarea_resize),  NULL);

   gtk_widget_set_can_focus(gl_area, TRUE);
   gtk_widget_set_focusable(gl_area, TRUE);

   gtk_widget_set_hexpand(gl_area, TRUE);
   gtk_widget_set_vexpand(gl_area, TRUE);

   return gl_area;

}


void on_glarea_scale_changed(GtkGestureZoom* self,
                             gdouble scale,
                             gpointer user_data) {
   graphics_info_t g;
   std::cout << "on_glarea_scale_changed " << scale << std::endl;
   // mouse_zoom() expects args (delta-x, delta-y)
   // we need to convert scale into something like that.
   g.mouse_zoom_by_scale_factor(scale);
}

void on_glarea_drag_begin_primary(GtkGestureDrag *gesture,
                                  double          x,
                                  double          y,
                                  GtkWidget      *area) {

   graphics_info_t g;

   // if (g.using_trackpad)
   //    g.on_glarea_drag_begin_secondary(gesture, x, y, area);
   // else
   //    g.on_glarea_drag_begin_primary(gesture, x, y, area);

   g.on_glarea_drag_begin_primary(gesture, x, y, area);

}

void on_glarea_drag_update_primary(GtkGestureDrag *gesture,
                                   double          delta_x,
                                   double          delta_y,
                                   GtkWidget      *area) {

   graphics_info_t g;

   // if (g.using_trackpad)
   //    // Hack for mac. Needs more thought.
   //    g.on_glarea_drag_update_secondary(gesture, delta_x, delta_y, area);
   // else
   //    g.on_glarea_drag_update_primary(gesture, delta_x, delta_y, area);

   g.on_glarea_drag_update_primary(gesture, delta_x, delta_y, area);

}

void on_glarea_drag_end_primary(GtkGestureDrag *gesture,
                                double          x,
                                double          y,
                                GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_end_primary(gesture, x, y, area);
}


void on_glarea_drag_begin_secondary(GtkGestureDrag *gesture,
                                    double          x,
                                    double          y,
                                    GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_begin_secondary(gesture, x, y, area);
}

void on_glarea_drag_update_secondary(GtkGestureDrag *gesture,
                                     double          delta_x,
                                     double          delta_y,
                                     GtkWidget      *area) {

   graphics_info_t g;
   g.on_glarea_drag_update_secondary(gesture, delta_x, delta_y, area);
}

void on_glarea_drag_end_secondary(GtkGestureDrag *gesture,
                                  double          x,
                                  double          y,
                                  GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_end_secondary(gesture, x, y, area);
}



void on_glarea_drag_begin_middle(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_begin_middle(gesture, x, y, area);
}

void on_glarea_drag_update_middle(GtkGestureDrag *gesture,
                                  double          delta_x,
                                  double          delta_y,
                                  GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_update_middle(gesture, delta_x, delta_y, area);
}

void on_glarea_drag_end_middle(GtkGestureDrag *gesture,
                               double          x,
                               double          y,
                               GtkWidget      *area) {
   graphics_info_t g;
   g.on_glarea_drag_end_middle(gesture, x, y, area);
}


gboolean
on_glarea_key_controller_key_pressed(GtkEventControllerKey *controller,
                                     guint                  keyval,
                                     guint                  keycode,
                                     guint                  modifiers,
                                     GtkButton             *button) {

   graphics_info_t g;
   // allow other controllers to act (say TAB has been pressed)
   gboolean handled = g.on_glarea_key_controller_key_pressed(controller, keyval, keycode, modifiers);
   return gboolean(handled);
}


void
on_glarea_key_controller_key_released(GtkEventControllerKey *controller,
                                      guint                  keyval,
                                      guint                  keycode,
                                      guint                  modifiers,
                                      GtkButton             *button) {

   graphics_info_t g;
   g.on_glarea_key_controller_key_released(controller, keyval, keycode, modifiers);
}


void
on_glarea_click(GtkGestureClick* click_gesture,
                gint n_press,
                gdouble x,
                gdouble y,
                gpointer user_data) {

   graphics_info_t g;
   g.on_glarea_click(click_gesture, n_press, x, y, user_data);

   // Not in Gtk4.
   // GtkWidget *w;
   // GtkWindow *window = gtk_widget_get_window(w);
   // GtkWidget *focused_widget = gtk_window_get_focus(window);

}

void
on_glarea_scrolled(GtkEventControllerScroll *controller,
                   double                    dx,
                   double                    dy,
                   gpointer                  user_data) {

   graphics_info_t g;
   g.on_glarea_scrolled(controller, dx, dy, user_data);

}

void
on_glarea_swipe(GtkGestureSwipe *controller,
                double           vel_x,
                double           vel_y,
                gpointer         user_data) {

   graphics_info_t g;
   // swipes happend a lot - a click-and-drag seems to be a swipe

   // std::cout << "------------------ swipe " << vel_x << " " << vel_y << std::endl;

   // 20240414-PE no, because click-drag-release is a swipe.
   // g.using_trackpad = true;

}

void
on_glarea_motion(GtkEventControllerMotion *controller,
                 gdouble x,
                 gdouble y,
                 gpointer user_data) {

   graphics_info_t g;
   g.on_glarea_motion(controller, x, y, user_data);
}

void
on_glarea_motion_enter(GtkEventControllerMotion *controller,
                       gdouble                   x,
                       gdouble                   y,
                       GdkCrossingMode           mode,
                       gpointer                  user_data) {
   // std::cout << "Motion enter" << std::endl;
}

void
on_glarea_motion_leave(GtkEventControllerMotion *controller,
                       GdkCrossingMode           mode,
                       gpointer                  user_data) {

   // std::cout << "Motion leave" << std::endl;
}

void setup_gestures_for_opengl_widget_in_main_window(GtkWidget *glarea) {

   // std::cout << "========== start setting up GTK4 style event controlllers" << std::endl;

   GtkGesture *zoom_controller           = gtk_gesture_zoom_new();
   g_signal_connect(zoom_controller, "scale-changed", G_CALLBACK(on_glarea_scale_changed), glarea);
   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(zoom_controller));
   GtkEventController *key_controller = gtk_event_controller_key_new();

   g_signal_connect(key_controller, "key-pressed",  G_CALLBACK(on_glarea_key_controller_key_pressed),  glarea);
   g_signal_connect(key_controller, "key-released", G_CALLBACK(on_glarea_key_controller_key_released), glarea);
   gtk_widget_add_controller(GTK_WIDGET(glarea), key_controller);

   GtkGesture *drag_controller_secondary = gtk_gesture_drag_new();
   GtkGesture *drag_controller_primary   = gtk_gesture_drag_new();
   GtkGesture *drag_controller_middle    = gtk_gesture_drag_new();
   GtkGesture *click_controller          = gtk_gesture_click_new();
   GtkGesture *swipe_controller          = gtk_gesture_swipe_new();

   GtkEventControllerScrollFlags scroll_flags = GTK_EVENT_CONTROLLER_SCROLL_VERTICAL;
   GtkEventController *scroll_controller = gtk_event_controller_scroll_new(scroll_flags);

   // #ifdef __APPLE__
   //    mouse_view_rotate_button_mask = GDK_BUTTON1_MASK; // GDK_BUTTON_PRIMARY
   //    mouse_pick_button_mask        = GDK_BUTTON1_MASK; // GDK_BUTTON_PRIMARY
   // #endif

   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_primary), GDK_BUTTON_PRIMARY);

   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(drag_controller_primary));
   g_signal_connect(drag_controller_primary, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_primary),  glarea);
   g_signal_connect(drag_controller_primary, "drag-update", G_CALLBACK(on_glarea_drag_update_primary), glarea);
   g_signal_connect(drag_controller_primary, "drag-end",    G_CALLBACK(on_glarea_drag_end_primary),    glarea);

   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_secondary), GDK_BUTTON_SECONDARY);

   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(drag_controller_secondary));
   g_signal_connect(drag_controller_secondary, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_secondary),  glarea);
   g_signal_connect(drag_controller_secondary, "drag-update", G_CALLBACK(on_glarea_drag_update_secondary), glarea);
   g_signal_connect(drag_controller_secondary, "drag-end",    G_CALLBACK(on_glarea_drag_end_secondary),    glarea);

   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_middle), GDK_BUTTON_MIDDLE);

   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(drag_controller_middle));
   g_signal_connect(drag_controller_middle, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_middle),  glarea);
   g_signal_connect(drag_controller_middle, "drag-update", G_CALLBACK(on_glarea_drag_update_middle), glarea);
   g_signal_connect(drag_controller_middle, "drag-end",    G_CALLBACK(on_glarea_drag_end_middle),    glarea);

   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(click_controller));
   g_signal_connect(click_controller, "pressed",  G_CALLBACK(on_glarea_click),  glarea);

   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(scroll_controller));
   g_signal_connect(scroll_controller, "scroll",  G_CALLBACK(on_glarea_scrolled),  glarea);

   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(swipe_controller));
   g_signal_connect(swipe_controller, "swipe",  G_CALLBACK(on_glarea_swipe),  glarea);

   GtkEventController *motion_controller = gtk_event_controller_motion_new();
   gtk_event_controller_set_propagation_phase(motion_controller, GTK_PHASE_CAPTURE);
   gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(motion_controller));
   g_signal_connect(motion_controller, "motion", G_CALLBACK(on_glarea_motion),       glarea);
   g_signal_connect(motion_controller, "enter",  G_CALLBACK(on_glarea_motion_enter), glarea);
   g_signal_connect(motion_controller, "leave",  G_CALLBACK(on_glarea_motion_leave), glarea);

}

void
install_icons_into_theme(GtkWidget *w) {


   // This is how to lookup an icon
   //
   // const char *fallbacks[] = {"one", "two", "three"};
   // GtkIconTheme *icon_theme = gtk_icon_theme_get_for_display(gtk_widget_get_display (w));
   // GtkTextDirection td = GTK_TEXT_DIR_LTR;
   // GtkIconLookupFlags icon_lookup_flags = GTK_ICON_LOOKUP_FORCE_REGULAR;
   //
   // GtkIconPaintable *icon = gtk_icon_theme_lookup_icon(icon_theme,
   //                                                     "my-icon-name", // icon name
   //                                                     fallbacks,
   //                                                     48, // icon size
   //                                                     1,  // scale
   //                                                     td,
   //                                                     icon_lookup_flags);
   // GdkPaintable *paintable = GDK_PAINTABLE (icon);
   // // Use the paintable
   // g_object_unref(icon);


#if 0 // just testing
   const char *theme_name = gtk_icon_theme_get_theme_name(icon_theme);
   if (theme_name)
      std::cout << "===== theme-name: " << theme_name << " === " << std::endl;
   else
      std::cout << "===== null theme-name === " << std::endl;
#endif

   GtkIconTheme *icon_theme = gtk_icon_theme_get_for_display(gtk_widget_get_display(w));
   std::string pkg_data_dir = coot::package_data_dir();
   std::string pixmap_dir = coot::util::append_dir_dir(pkg_data_dir, "pixmaps");
   std::string pixmap_dark_dir = coot::util::append_dir_dir(pixmap_dir, "dark");
   gtk_icon_theme_add_search_path(icon_theme, pixmap_dark_dir.c_str());
   gtk_icon_theme_add_search_path(icon_theme, pixmap_dir.c_str());

   // This is only necessary when coot is installed in a non-standard location
   // i.e. other than /usr or ~/.local (or /usr/local perhaps?)
   // which makes it convenient for testing without system-wide coot installation
   std::string prefix_dir = coot::prefix_dir();
   std::string icons_dir = coot::util::append_dir_dir(prefix_dir, "share/icons");
   gtk_icon_theme_add_search_path(icon_theme, icons_dir.c_str());

}


void
on_go_to_residue_keyboarding_mode_entry_key_controller_key_released(GtkEventControllerKey *controller,
                                                                    guint                  keyval,
                                                                    guint                  keycode,
                                                                    guint                  modifiers,
                                                                    GtkEntry              *entry) {

   std::cout << "in on_go_to_residue_keyboarding_mode_entry_key_controller_key_released() "
             << keycode << std::endl;
   GtkWidget *window = widget_from_builder("keyboard_go_to_residue_window");

   if (keycode == 36) {
      std::string s = gtk_editable_get_text(GTK_EDITABLE(entry));
      graphics_info_t g;
      g.apply_go_to_residue_keyboading_string(s);
      gtk_editable_set_text(GTK_EDITABLE(entry), "");
      gtk_widget_set_visible(GTK_WIDGET(window), FALSE);
   }

   if (keycode == 53) {
      gtk_widget_set_visible(GTK_WIDGET(window), FALSE);
      gtk_editable_set_text(GTK_EDITABLE(entry), "");
   }
}

void setup_go_to_residue_keyboarding_mode_entry_signals() {
   GtkWidget *entry = widget_from_builder("keyboard_go_to_residue_entry");
   if (entry) {
      GtkEventController *key_controller = gtk_event_controller_key_new();
      g_signal_connect(key_controller, "key-released", G_CALLBACK(on_go_to_residue_keyboarding_mode_entry_key_controller_key_released), entry);
      gtk_widget_add_controller(GTK_WIDGET(entry), key_controller);
   }
}

// maybe this function can be in a better home.
void
handle_start_scripts() {

   // 20240609-PE note to self scm_c_primitive_load() fails with a crash because we
   // have not done the scm_boot_guile() call (g_application_run() is called where
   // scm_boot_guile() should be called. I don't know what to do).

   xdg_t xdg;
   std::vector<std::filesystem::path> scripts;
#ifdef USE_GUILE
   scripts = xdg.get_scheme_config_scripts();
   for (const auto &script : scripts) {
      std::cout << "Load scheme config script " << script.c_str() << " (ignored)" << std::endl;
      // scm_c_primitive_load(script.c_str());
   }
#endif
   scripts = xdg.get_python_config_scripts();
   for (const auto &script : scripts) {
      // std::cout << "Load python config script " << script.c_str() << std::endl;
      logger.log(log_t::INFO, logging::function_name_t(__FUNCTION__),
		 "Load python script", script);
      run_python_script(script.string().c_str());
   }
#ifdef USE_GUILE
   if (graphics_info_t::run_startup_scripts_flag) {
      if (graphics_info_t::run_state_file_status) {
         std::pair<bool, std::filesystem::path> script = xdg.get_scheme_state_script();
         if (script.first) {
            // scm_c_primitive_load(script.second.c_str());
         }
      }
   }
#endif
   if (graphics_info_t::run_startup_scripts_flag) {
      if (graphics_info_t::run_state_file_status) {
         std::pair<bool, std::filesystem::path> script = xdg.get_python_state_script();
         if (script.first) {
            run_python_script(script.second.string().c_str());
         }
      }
   }
}


GtkWidget *
create_local_picture(const std::string &local_filename) {

   GtkWidget *picture = 0;

   std::string pdd = coot::package_data_dir();
   // std::cout << "pdd " << pdd << std::endl;
   std::string icon_dir = coot::util::append_dir_file(pdd, "images");
   std::vector<std::string> pixmap_directories_gtk4 = {};
   pixmap_directories_gtk4.push_back(icon_dir);

   if (local_filename.empty())
      return 0;

   if (coot::file_exists(std::string("./") + local_filename)) {
      picture = gtk_picture_new_for_filename(local_filename.c_str());
   }
   else {
        for (unsigned int i=0; i<pixmap_directories_gtk4.size(); i++) {
           const std::string &d = pixmap_directories_gtk4[i];
           std::string fn = coot::util::append_dir_file(d, local_filename);
           if (coot::file_exists(fn)) {
              picture = gtk_picture_new_for_filename(fn.c_str());
              break;
           }
        }
   }
   if (picture) {
#if GTK_MAJOR_VERSION == 4 && GTK_MINOR_VERSION >= 8
      gtk_picture_set_content_fit(GTK_PICTURE(picture), GTK_CONTENT_FIT_FILL);
#else
      g_warning("gtk_picture_set_content_fit() not available in your version of GTK.");
#endif
   }
   return picture;
}

GtkWidget*
new_startup_create_splash_screen_window() {

   GtkWidget *splash_screen_window = gtk_window_new();
   gtk_window_set_title(GTK_WINDOW(splash_screen_window), "Coot-Splash");
   gtk_window_set_decorated(GTK_WINDOW(splash_screen_window), FALSE);
   GtkWidget *picture = create_local_picture("coot-1.1.20.png");

   gtk_widget_set_hexpand(GTK_WIDGET(picture),TRUE);
   gtk_widget_set_vexpand(GTK_WIDGET(picture),TRUE);

   gtk_widget_set_halign(GTK_WIDGET(picture),GTK_ALIGN_FILL);
   gtk_widget_set_halign(GTK_WIDGET(picture),GTK_ALIGN_FILL);

   gtk_widget_set_size_request(picture, 660, 371);
   // std::cout << "@@@@@@@@@@@@@@ create_pixmap_gtk4_version() returned image " << image << std::endl;
   gtk_widget_set_visible(picture, TRUE);

   gtk_window_set_child(GTK_WINDOW(splash_screen_window), picture);
   return splash_screen_window;
}

// needs to be packed in a gpointer (use dynamic allocation)
struct application_activate_data {
   int argc;
   char** argv;
   GtkWidget* splash_screen;
   GtkApplication* application;
   GtkWidget* app_window;
   command_line_data cld;

   application_activate_data(int _argc, char** _argv, command_line_data&& cld) {
      argc = _argc;
      argv = _argv;
      splash_screen = nullptr;
      application = nullptr;
      app_window = nullptr;
      this->cld = std::move(cld);
   }
};

gboolean
on_app_window_key_controller_key_pressed(GtkEventControllerKey *controller,
                                     guint                  keyval,
                                     guint                  keycode,
                                     guint                  modifiers,
                                     GtkButton             *button) {

   std::cout << "app window keypressed! " << keyval << " " << keycode << std::endl;
   gboolean handled = FALSE;
   return handled;
}

void
add_key_bindings_for_application_window(GtkWidget *app_window) {

   GtkEventController *key_controller = gtk_event_controller_key_new();
   g_signal_connect(key_controller, "key-pressed",  G_CALLBACK(on_app_window_key_controller_key_pressed), app_window);
   gtk_widget_add_controller(app_window, key_controller);
}

// drag and drop code needs to be reworked. Add this here for now.
int handle_drag_and_drop_string(const std::string &file_name);

void
new_startup_application_activate(GtkApplication *application,
                                 gpointer user_data) {

   application_activate_data* activate_data = (application_activate_data*) user_data;

   activate_data->application = application;

   std::string window_name = "Coot-" + std::string(VERSION);

#ifdef WINDOWS_MINGW
   window_name = "WinCoot-" + std::string(VERSION);
#endif

   GtkWidget *app_window = gtk_application_window_new(application);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_title(GTK_WINDOW(app_window), window_name.c_str());

   add_key_bindings_for_application_window(app_window);

   graphics_info_t::set_main_window(app_window);

   activate_data->app_window = app_window;

   g_idle_add(+[](gpointer user_data) -> gboolean {

      application_activate_data* activate_data = static_cast<application_activate_data*>(user_data);

      GtkWindow* splash_screen = GTK_WINDOW(activate_data->splash_screen);
      GtkWidget* app_window = activate_data->app_window;
      GtkApplication* application = GTK_APPLICATION(activate_data->application);
      int argc = activate_data->argc;
      char** argv = activate_data->argv;

      auto python_init = [argc, argv] () {
         setup_python_basic(argc, argv);
         setup_python_coot_module();

         // this needs the gtkbuilder to have read the ui file
         // because it needs to look up  the coot_main_window
         // and main_toolbar and main_hbox and main_statusbar.
         // setup_python_with_coot_modules(argc, argv);
         // So it is done in new_startup_application_activate().

      };

      graphics_info_t graphics_info;
      graphics_info.init(); // added 20241231

      // use this to look up things - and it is used to attach the lidia
      // application window
      graphics_info.application = application;

      // 20230526-PE this now happens in init_coot_as_python_module()
      // Let's not do it (calling geom.init_standard()) twice.
      // graphics_info.init();

      // but let's do it once at least!

      // this is done in the python startup now.
      // std::cout << "#################### new_startup_application_activate()  calling graphics_info.init() "
      //           << std::endl;
      // graphics_info.init();

      GtkBuilder *builder = gtk_builder_new();
      if (GTK_IS_BUILDER(builder)) {
      } else {
         std::cout << "ERROR:: in new_startup_application_activate() builder was NOT a builder"
                  << std::endl;
         coot_no_state_real_exit(0);
      }

      install_icons_into_theme(GTK_WIDGET(app_window));
      gtk_window_set_icon_name(GTK_WINDOW(app_window), "coot");

      // the main application builder

      std::string dir = coot::package_data_dir();
      // change "glade" to "ui" one day.
      // 20240218-PE today is that day!
      std::string dir_ui = coot::util::append_dir_dir(dir, "ui");
      std::string ui_file_name = "coot-gtk4.ui";
      std::string ui_file_full = coot::util::append_dir_file(dir_ui, ui_file_name);
      if (coot::file_exists(ui_file_name))
         ui_file_full = ui_file_name;

      GError* error = NULL;
      gboolean status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);
      if (status == FALSE) {
         std::cout << "ERROR:: Failure to read or parse " << ui_file_full << std::endl;
         std::cout << error->message << std::endl;
         coot_no_state_real_exit(0);
      }

      // the preferences builder:
      std::string preferences_ui_file_name = "preferences-gtk4.ui";
      std::string preferences_ui_file_name_full = coot::util::append_dir_file(dir_ui, preferences_ui_file_name);
      if (coot::file_exists(preferences_ui_file_name))
         preferences_ui_file_name_full = preferences_ui_file_name;
      GtkBuilder *preferences_builder = gtk_builder_new();
      // std::cout << "::::::::::::::::::::::: reading " << preferences_ui_file_name_full << std::endl;
      status = gtk_builder_add_from_file(preferences_builder, preferences_ui_file_name_full.c_str(), &error);
      // std::cout << "::::::::::::::::::::::: done reading " << preferences_ui_file_name_full << std::endl;
      if (status == FALSE) {
         std::cout << "ERROR:: Failure to read or parse " << preferences_ui_file_name_full << std::endl;
         std::cout << error->message << std::endl;
         coot_no_state_real_exit(0);
      }
      graphics_info_t::set_preferences_gtkbuilder(preferences_builder);

      // set the version in the about dialog
      GtkWidget *about_dialog = GTK_WIDGET(gtk_builder_get_object(builder, "about_dialog"));
      std::string version_str = std::string(VERSION);
      if (version_str.find("-pre") != std::string::npos) {
         version_str += "\n";
         version_str += git_commit();
         std::string s = COOT_BUILD_INFO_STRING;
         if (! s.empty()) {
            version_str += "\n";
            version_str += s;
         }
      }
      // override the value in the coot-gtk4.ui file.
      gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(about_dialog), version_str.c_str());


      python_init();

      // 20231114-PE we can't handle the command line data until the graphics have started.
      // so this should be in the realize() function for the graphics widget.
      // handle_command_line_data(activate_data->cld);

      if (activate_data->cld.do_graphics)
         graphics_info.use_graphics_interface_flag = true;

      // create the preference defaults
      make_preferences_internal();

      guint id = gtk_application_window_get_id(GTK_APPLICATION_WINDOW(app_window));
      // std::cout << "debug:: new_startup_application_activate(): Window id: " << id << std::endl;

      graphics_info_t::set_gtkbuilder(builder);

      // GMenu *menu = create_menu_by_hand(application);
      GMenu *menubar = G_MENU(graphics_info_t::get_gobject_from_builder("menubar"));
      gtk_application_set_menubar(application, G_MENU_MODEL(menubar));
      gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(app_window), TRUE);

      // GtkWidget *graphics_hbox = widget_from_builder("crows_graphics_hbox", builder);
      // GtkWidget *main_window   = widget_from_builder("crows_main_window",   builder);
      GtkWidget *graphics_hbox = widget_from_builder("main_window_graphics_hbox");
      GtkWidget *graphics_vbox = widget_from_builder("main_window_vbox");
      gtk_window_set_child(GTK_WINDOW(app_window), graphics_vbox);

      gtk_window_present(GTK_WINDOW(app_window));

      g_signal_connect(app_window, "destroy", G_CALLBACK(+[](GtkWidget *w, gpointer user_data){
         if(coot::is_layla_initialized()) {
            g_info("De-initializing Layla so that GtkApplication can exit...");
            coot::deinitialize_layla();
         }
      }), nullptr);

      // gtk_widget_set_visible(window, TRUE);

      GtkWidget *gl_area = new_startup_create_glarea_widget();
      graphics_info_t::glareas.push_back(gl_area);
      gtk_widget_set_visible(gl_area, TRUE);
      gtk_box_prepend(GTK_BOX(graphics_hbox), gl_area);
      gtk_window_set_application(GTK_WINDOW(app_window), application);
#ifdef __APPLE__
      gtk_widget_set_size_request(gl_area, 600, 600); // Hmm
      gtk_window_set_default_size(GTK_WINDOW(app_window), 700, 700);
      gtk_window_set_default_widget(GTK_WINDOW(app_window), gl_area);
#else
      gtk_window_set_focus_visible(GTK_WINDOW(app_window), TRUE);

      // 20230729-PE
      // gtk_widget_set_size_request() does't seem to work on the gl_area.
      // So expand the gl_area by setting thw window size just so. This makes the
      // gl_area 900x900 on my desktop. Maybe there is a better way.
      // The console show that new_startup_on_glarea_resize() is called several times:
      // DEBUG:: --- new_startup_on_glarea_resize() 900 900
      // DEBUG:: --- new_startup_on_glarea_resize() 900 710
      // DEBUG:: --- new_startup_on_glarea_resize() 900 900
      // Curious.
      gtk_window_set_default_size(GTK_WINDOW(app_window), 1076, 1023);
      gtk_window_set_default_widget(GTK_WINDOW(app_window), gl_area);
      gtk_widget_set_visible(app_window, TRUE);
      gtk_window_set_focus_visible(GTK_WINDOW(app_window), TRUE);
#endif

      // ---------------------  -----------------------

      // drag and drop: well, just drop for the moment:
      //
      // Set up drop target.
      // GType types[2] = { GDK_TYPE_RGBA, G_TYPE_STRING };
      GType types[7] = { GDK_TYPE_RGBA, G_TYPE_STRING, G_TYPE_PARAM,
                         G_TYPE_OBJECT, G_TYPE_VARIANT, G_TYPE_FILE, G_TYPE_STRV};

      GtkDropTarget *drop_target = gtk_drop_target_new(G_TYPE_STRING, GDK_ACTION_COPY);
      gtk_drop_target_set_gtypes (drop_target, types, G_N_ELEMENTS (types));
      gtk_widget_add_controller(GTK_WIDGET(gl_area), GTK_EVENT_CONTROLLER(drop_target));

      auto on_drop_performed = +[] (GtkDropTarget *drop_target, const GValue *value, double x, double y) {

         // return a gboolean
         gboolean status = FALSE;

         g_print("DEBUG:: Drop performed!\n");
         GType type = G_VALUE_TYPE(value);
         std::cout << "DEBUG:: value is of type " << type << std::endl;

         if (G_VALUE_HOLDS(value, G_TYPE_FILE)) {
            std::cout << "!!!!!!!!!!!!! value holds a file!" << std::endl;
            GFile *file = (GFile *)g_value_get_object(value);
            if (file) {
               std::cout << "DEBUG:: got file: " << file << std::endl;
               const gchar *filename = g_file_get_path(file);
               std::cout << "DEBUG:: got filename: " << filename << std::endl;
               handle_drag_and_drop_string(filename);
               status = TRUE;
            } else {
               std::cout << "got null file " << std::endl;
            }
         }

         if (type == G_TYPE_OBJECT) {
            std::cout << "G_TYPE_OBJECT! " << std::endl;
         } else {
            std::cout << "not type G_TYPE_OBJECT! " << std::endl;
         }

         if (type == G_TYPE_STRV) {
            std::cout << "G_TYPE_STRV! " << std::endl;
         } else {
            std::cout << "not type G_TYPE_STRV! " << std::endl;
         }

         if (type == G_TYPE_FILE) {
            std::cout << "G_TYPE_FILE! " << std::endl;
         } else {
            std::cout << "not type G_TYPE_FILE! " << std::endl;
         }

         if (type == G_TYPE_PARAM) {
            std::cout << "G_TYPE_PARAM! " << std::endl;
         } else {
            std::cout << "not type G_TYPE_PARAM! " << std::endl;
         }

         if (type == G_TYPE_VARIANT) {
            std::cout << "G_TYPE_VARIANT! " << std::endl;
         } else {
            std::cout << "not type G_TYPE_VARIANT! " << std::endl;
         }

         if (type == G_TYPE_POINTER) {
            std::cout << "G_TYPE_POINTER! " << std::endl;
         } else {
            std::cout << "not type G_TYPE_POINTER! " << std::endl;
         }

         if (type == G_TYPE_STRING) {
            std::cout << "G_TYPE_STRING! " << std::endl;
            const char *text = g_value_get_string(value);
            if (text) {
               unsigned long ll = strlen(text);
               std::cout << "DEBUG:: text has length " << ll << std::endl;
               if (ll > 0) {
                  handle_drag_and_drop_string(text);
                  status = TRUE;
               }
            } else {
               std::cout << "DEBUG:: on_drop_performed(): text: was null" << std::endl;
            }
         } else {
            std::cout << "not type G_TYPE_STRING! " << std::endl;
         }
         std::cout << "DEBUG:: returning from on_drop_performed()." << std::endl;
         return status;
      };
      g_signal_connect(drop_target, "drop", G_CALLBACK(on_drop_performed), NULL);

      // ------------------ no screenshot for macOS  -----------------------

#ifdef __APPLE__
      // GtkWidget *menu_item = widget_from_builder("screenshot-menu-item");
      // if (menu_item)
      // gtk_label_set_text(GTK_LABEL(menu_item), "Screenshot Not Available");
#endif

      // ---------------------  -----------------------

      gtk_widget_grab_focus(gl_area); // at the start, fixes focus problem
      setup_gestures_for_opengl_widget_in_main_window(gl_area);

      create_actions(application);

      setup_python_with_coot_modules(argc, argv);

      setup_gui_components();
      setup_go_to_residue_keyboarding_mode_entry_signals();

      handle_start_scripts(); // what used to be in ~/.coot/*.py.
                              // 20260124-PE Not to self: these are not
                              // command-line scripts.

      // now we are ready to show graphical objects made from reading files:
      handle_command_line_data(activate_data->cld);


      // 20251019-PE is this the first time Coot-1 has been started?
      {
         bool show_first_startup_dialog = true;
         xdg_t xdg;
         std::filesystem::path state_home = xdg.get_state_home();
         if (std::filesystem::exists(state_home)) {
            std::filesystem::path state_py = state_home / "0-coot.state.py";
            if (std::filesystem::exists(state_py)) {
               show_first_startup_dialog = false;
            }
         }
         if (show_first_startup_dialog) {
            GtkWidget *dialog = widget_from_builder("first-startup-dialog");
            GtkWidget *main_window_widget = graphics_info_t::get_main_window();
            if (main_window_widget) {
               GtkWindow *main_window = GTK_WINDOW(main_window_widget);
               gtk_window_set_transient_for(GTK_WINDOW(dialog), main_window);
            }
            gtk_widget_set_visible(dialog, TRUE);
         }
      }

      // load_tutorial_model_and_data();
      delete activate_data;

      auto destroy_splash_screen_callback = +[] (gpointer data) {
         if (data) {
            GtkWindow* splash_screen = GTK_WINDOW(data);
            gtk_window_destroy(splash_screen);
         }
         return G_SOURCE_REMOVE;
      };
      g_idle_add(destroy_splash_screen_callback, splash_screen);

      return G_SOURCE_REMOVE;
   }, activate_data);


}

// move these to the top.
void setup_symm_lib();
void check_reference_structures_dir();

void load_css() {

   std::string fn = "coot.css";
   if (coot::file_exists(fn)) {
      GdkDisplay *display = gdk_display_get_default();
      GtkCssProvider *provider = gtk_css_provider_new();
      GFile *gf = g_file_new_for_path(fn.c_str());
      gtk_css_provider_load_from_file(provider, gf);
      gtk_style_context_add_provider_for_display(display,
                                                 GTK_STYLE_PROVIDER(provider),
                                                 GTK_STYLE_PROVIDER_PRIORITY_FALLBACK);
      g_object_unref(provider);
   }

}

void window_removed(GtkApplication* self,GtkWindow* window, gpointer user_data) {

   // this is not needed because closing the main window using the window manager
   // causes g_application_run() in the function below to return. Hence we
   // just fall out at the end of main().
   //
   // Or that's what *should* happen.

   // std::cout << "quit here" << std::endl;
   // g_application_quit(self);

}

int do_no_graphics_mode(command_line_data& cld, int argc, char** argv) {

   handle_command_line_data(cld);
   setup_python_basic(argc, argv);
   setup_python_coot_module();
   setup_python_with_coot_modules(argc, argv);
   run_command_line_scripts();
   start_command_line_python_maybe(true, argc, argv);

   return 0;
}

int
do_self_tests() {

   std::cout << "INFO:: Running internal self tests" << std::endl;

   // return true on success
   clipper::Test_core test_core;       bool result_core    = test_core();
   clipper::Test_contrib test_contrib; bool result_contrib = test_contrib();
   std::cout<<" INFO:: Test Clipper core   : "<<(result_core   ?"OK":"FAIL")<<std::endl;
   std::cout<<" INFO:: Test Clipper contrib: "<<(result_contrib?"OK":"FAIL")<<std::endl;

   // 20240309-PE I need tests
   //   1: internal tests (that can use tutorial-modern and rnasa)
   //   2: internal tests that use the monomer library
   //   3: internal tests that use greg data
   // executed by --self-test-internal, --self-test-monomer-library --self-test-greg-data.
   // Currently I only have type 1:

   // return 1 on success
   int gis = test_internal();
   int shell_exit_code = 1;
   if (result_core)
      if (result_contrib)
         if (gis == 1)
            shell_exit_code = 0;
   return shell_exit_code;

}


int new_startup(int argc, char **argv) {

#ifdef USE_LIBCURL
   curl_global_init(CURL_GLOBAL_NOTHING); // nothing extra (e.g. ssl or WIN32)
#endif

   coot::set_realpath_for_coot_executable(argv[0]);

   mmdb::InitMatType();

   // setup_symm_lib();
   // check_reference_structures_dir();

   command_line_data cld = parse_command_line(argc, argv);

   if (cld.run_internal_tests_and_exit)
      return do_self_tests();

   if(!cld.do_graphics) {
      return do_no_graphics_mode(cld, argc, argv);
   }

   if (cld.use_opengl_es)
      graphics_info_t::graphics_is_gl_es = true;

   gtk_init();

   load_css();

   // Tell us the GTK version
   std::string gtk_version_string =
      std::to_string(GTK_MAJOR_VERSION) + "." +
      std::to_string(GTK_MINOR_VERSION) + "." +
      std::to_string(GTK_MICRO_VERSION);
   logger.log(log_t::INFO, "Built with GTK", gtk_version_string);

   GtkWidget *splash_screen = nullptr;
   if (cld.use_splash_screen) {
      splash_screen = new_startup_create_splash_screen_window();
      gtk_widget_set_visible(splash_screen, TRUE);
   }

   g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);
   // Here's how you access that:
   // gboolean dark_mode_flag = FALSE;
   // g_object_get(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", &dark_mode_flag, NULL);

   // dark mode vs light-mode switch
   // 20251215-PE was widget_from_preferences_builder("light-mode-dark-mode-switch");
   // I think libadwaita is needed for mode switch
   GtkWidget *mode_switch = nullptr;
   if (mode_switch) {

      auto mode_switch_callback = +[] (GtkSwitch *sw, gboolean state, gpointer user_data) {

        GtkSettings *settings = gtk_settings_get_default();

        // 'state' is TRUE if the user just clicked "On"
        // 'state' is FALSE if the user just clicked "Off"
        if (state) {
            //TO Dark Mode
            g_object_set(settings, "gtk-application-prefer-dark-theme", TRUE, NULL);
        } else {
            // TO Light Mode
            g_object_set(settings, "gtk-application-prefer-dark-theme", FALSE, NULL);
        }

        // Return FALSE to allow the switch to complete the animation/toggle
        return gboolean(FALSE);
        };

         gpointer *user_data = nullptr;
         g_signal_connect(G_OBJECT(mode_switch), "state-set", G_CALLBACK(mode_switch_callback), user_data);
   }

   GError *error = NULL;
#if GLIB_MAJOR_VERSION == 2 && GLIB_MINOR_VERSION >= 74 || GLIB_MAJOR_VERSION > 2
   GtkApplication *app = gtk_application_new ("org.emsley.coot",
      (GApplicationFlags) (G_APPLICATION_DEFAULT_FLAGS | G_APPLICATION_NON_UNIQUE));
#else
   GtkApplication *app = gtk_application_new ("org.emsley.coot",
      (GApplicationFlags) (G_APPLICATION_NON_UNIQUE));
#endif
   g_application_register(G_APPLICATION(app), NULL, &error);

   application_activate_data *activate_data = new application_activate_data(argc,argv,std::move(cld));
   activate_data->splash_screen = splash_screen;
   // this destroys active_data
   g_signal_connect(app, "activate", G_CALLBACK(new_startup_application_activate), activate_data);

   // how about this - needed for Bernie/Windows?
   // void window_removed ( GtkApplication* self, GtkWindow* window, gpointer user_data )
   g_signal_connect(app, "window-removed", G_CALLBACK(window_removed), nullptr);

   // delete activate_data; Nope. This is used in new_startup_application_activate.
   // Delete it there if you want to delete it.

   // read in inchikeys - is this the right place for this?
   // 2026-02-06-PE No. It should be in graphics_info_t::init(), as is ptm_database init.
   graphics_info_t::read_inchikeys();

   int status = g_application_run(G_APPLICATION(app), 1, argv);
   std::cout << "--- g_application_run() returns with status " << status << std::endl;
   g_object_unref(app);
   return status;
}
