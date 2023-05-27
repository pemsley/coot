
#include <iostream>
#include <string>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "graphics-info.h"
#include "create-menu-item-actions.hh"
#include "setup-gui-components.hh"
#include "coot-setup-python.hh"

void print_opengl_info();

void init_framebuffers(GtkWidget *glarea) {

   // put this into graphics-info I suppose.

   // std::cout << "DEBUG:: use_framebuffers: " << graphics_info_t::use_framebuffers << std::endl;

   std::cout << "----- start init_framebuffers() ----" << std::endl;

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

   std::cout << "----- done init_framebuffers() ----" << std::endl;
}


#include "text-rendering-utils.hh"

void
new_startup_realize(GtkWidget *gl_area) {

   // std::cout << "new_startup_realize() ------------------- start ------------------"
   //              << std::endl;

   gtk_gl_area_make_current(GTK_GL_AREA (gl_area));

   if (gtk_gl_area_get_error(GTK_GL_AREA (gl_area)) != NULL)
      return;

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
   g.init_framebuffers(w, h);
   g.init_joey_ssao_stuff(w, h);

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
   g.setup_draw_for_bad_nbc_atom_pair_markers();
   g.setup_draw_for_anchored_atom_markers_init();
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

   g.setup_key_bindings();

   GLenum err = glGetError();
   if (err)
      std::cout << "ERROR:: new_startup_realize() --start-- err is " << err << std::endl;

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

   if (false)
      std::cout << "DEBUG:: --- new_startup_on_glarea_resize() " <<  width << " " << height << std::endl;

   graphics_info_t g;
   // for the GL widget, not the window.
   g.graphics_x_size = width;
   g.graphics_y_size = height;
   g.reset_frame_buffers(width, height); // currently makes the widget blank (not drawn)
   g.resize_framebuffers_textures_renderbuffers(width, height); // 20220131-PE added from crows merge
   g.reset_hud_buttons_size_and_position();

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

// void on_glarea_realize(GtkWidget *widget); // using this give linking problems.

GtkWidget *new_startup_create_glarea_widget() {

   GtkWidget *gl_area = gtk_gl_area_new();
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
   double s = pow(scale, 1.02);
   std::cout << "on_glarea_scale_changed " << scale << " " << s << std::endl;
   g.mouse_zoom(s, 0.0);
}

void on_glarea_drag_begin_primary(GtkGestureDrag *gesture,
                                  double          x,
                                  double          y,
                                  GtkWidget      *area) {
   graphics_info_t g;

#ifdef __APPLE__
   g.on_glarea_drag_begin_secondary(gesture, x, y, area);
#else
   g.on_glarea_drag_begin_primary(gesture, x, y, area);
#endif
}

void on_glarea_drag_update_primary(GtkGestureDrag *gesture,
                                   double          delta_x,
                                   double          delta_y,
                                   GtkWidget      *area) {

   graphics_info_t g;

#ifdef __APPLE__
   // Hack for mac. Needs more thought.
   g.on_glarea_drag_update_secondary(gesture, delta_x, delta_y, area);
#else
   g.on_glarea_drag_update_primary(gesture, delta_x, delta_y, area);
#endif

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
   std::string icon_dir = coot::util::append_dir_dir(pkg_data_dir, "icons/hicolor/16x16/actions");
   gtk_icon_theme_add_search_path(icon_theme, icon_dir.c_str());
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
      gtk_widget_hide(GTK_WIDGET(window));
   }

   if (keycode == 53) {
      gtk_widget_hide(GTK_WIDGET(window));
      gtk_editable_set_text(GTK_EDITABLE(entry), "");
   }
}

// in screen-utils.cc
void setup_application_icon(GtkWindow *window);

void setup_go_to_residue_keyboarding_mode_entry_signals() {
   GtkWidget *entry = widget_from_builder("keyboard_go_to_residue_entry");
   if (entry) {
      GtkEventController *key_controller = gtk_event_controller_key_new();
      g_signal_connect(key_controller, "key-released", G_CALLBACK(on_go_to_residue_keyboarding_mode_entry_key_controller_key_released), entry);
      gtk_widget_add_controller(GTK_WIDGET(entry), key_controller);

   }
}


GtkWidget *
create_local_picture(const std::string &local_filename) {

   GtkWidget *picture = 0;

   std::string pdd = coot::package_data_dir();
   std::cout << "pdd " << pdd << std::endl;
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
#endif
   }
   return picture;
}

GtkWidget*
new_startup_create_splash_screen_window() {

   GtkWidget *splash_screen_window = gtk_window_new();
   gtk_window_set_title(GTK_WINDOW(splash_screen_window), "Coot-Splash");
   gtk_window_set_decorated(GTK_WINDOW(splash_screen_window), FALSE);
   GtkWidget *picture = create_local_picture("coot-1.png");

   gtk_widget_set_hexpand(GTK_WIDGET(picture),TRUE);
   gtk_widget_set_vexpand(GTK_WIDGET(picture),TRUE);

   gtk_widget_set_halign(GTK_WIDGET(picture),GTK_ALIGN_FILL);
   gtk_widget_set_halign(GTK_WIDGET(picture),GTK_ALIGN_FILL);

   gtk_widget_set_size_request(picture, 660, 371);
   // std::cout << "@@@@@@@@@@@@@@ create_pixmap_gtk4_version() returned image " << image << std::endl;
   gtk_widget_show(picture);

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

   application_activate_data(int _argc, char** _argv) {
      argc = _argc;
      argv = _argv;
      splash_screen = nullptr;
      application = nullptr;
      app_window = nullptr;
   }
};

#include "command-line.hh"

void
new_startup_application_activate(GtkApplication *application,
                                 gpointer user_data) {

   application_activate_data* activate_data = (application_activate_data*) user_data;

   activate_data->application = application;

   std::string window_name = "GTK4 Coot-" + std::string(VERSION);
   GtkWidget *app_window = gtk_application_window_new(application);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_title(GTK_WINDOW(app_window), window_name.c_str());
   setup_application_icon(GTK_WINDOW(app_window)); // 20220807-PE not sure what this does in gtk4 or if it works.
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

      graphics_info.application = application;

      // 20230526-PE this now happens i init_coot_as_python_module()
      // Let's not do it (including geom.init_standard()) twice.
      // graphics_info.init();

      GtkBuilder *builder = gtk_builder_new();
      if (GTK_IS_BUILDER(builder)) {
      } else {
         std::cout << "ERROR:: in new_startup_application_activate() builder was NOT a builder"
                  << std::endl;
         exit(0);
      }

      std::string dir = coot::package_data_dir();
      std::string dir_glade = coot::util::append_dir_dir(dir, "glade");
      std::string glade_file_name = "coot-gtk4.ui";
      std::string glade_file_full = coot::util::append_dir_file(dir_glade, glade_file_name);
      if (coot::file_exists(glade_file_name))
         glade_file_full = glade_file_name;

      GError* error = NULL;
      gboolean status = gtk_builder_add_from_file(builder, glade_file_full.c_str(), &error);
      if (status == FALSE) {
         std::cout << "ERROR:: Failure to read or parse " << glade_file_full << std::endl;
         std::cout << error->message << std::endl;
         exit(0);
      }

      python_init();

      // set this by parsing the command line arguments
      graphics_info.use_graphics_interface_flag = true;


      GtkWidget *sb = GTK_WIDGET(gtk_builder_get_object(builder, "main_window_statusbar"));
      graphics_info_t::statusbar = sb;
      // std::cout << "debug:: statusbar: " << sb << std::endl;

      install_icons_into_theme(GTK_WIDGET(sb));


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
      // gtk_widget_show(window);

      GtkWidget *gl_area = new_startup_create_glarea_widget();
      graphics_info_t::glareas.push_back(gl_area);
      gtk_widget_show(gl_area);
      gtk_box_prepend(GTK_BOX(graphics_hbox), gl_area);
      gtk_window_set_application(GTK_WINDOW(app_window), application);
      gtk_window_set_default_size(GTK_WINDOW(app_window), 1000, 900);
      gtk_window_set_default_widget(GTK_WINDOW(app_window), gl_area);
      gtk_widget_set_size_request(gl_area, 900, 900); // Hmm
      gtk_widget_show(app_window);
      gtk_window_set_focus_visible(GTK_WINDOW(app_window), TRUE);

      gtk_widget_grab_focus(gl_area); // at the start, fixes focus problem
      setup_gestures_for_opengl_widget_in_main_window(gl_area);

      create_actions(application);

      setup_python_with_coot_modules(argc, argv);

      setup_gui_components();
      setup_go_to_residue_keyboarding_mode_entry_signals();

      // if there is no command line arguments, the the function that sets this data is not run
      // so cld is null
      command_line_data *cld = static_cast<command_line_data *>(g_object_get_data(G_OBJECT(application),
                                                                                  "command-line-data"));
      if (cld) {
         handle_command_line_data(*cld);
         run_command_line_scripts();
      }

      // load_tutorial_model_and_data();

      g_idle_add(+[](gpointer data)-> gboolean {
         GtkWindow* splash_screen = GTK_WINDOW(data);
         gtk_window_destroy(splash_screen);
         return G_SOURCE_REMOVE;
      }, splash_screen);

      return G_SOURCE_REMOVE;
   }, activate_data);

   // delete activate_data; // 20230515-PE restore this when other command line stuff is working OK

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

void
application_open_callback(GtkApplication *app,
                          GFile          **files,
                          gint            n_files,
                          gchar          *hint,
                          gpointer        user_data) {

   
   command_line_data cld;

   for (gint i=0; i<n_files; i++) {
      GFile *file = files[i];
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_NAME,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      if (file_info) {
         const char *file_name = g_file_info_get_name(file_info);
         if (file_name) {
            std::cout << "Handle " << file_name << std::endl;
            cld.add(std::string(file_name));
         } else {
            std::cout << "file_name was null " << std::endl;
         }
      } else {
         std::cout << "application_open_callback() error " << i << " " << error->message << std::endl;
      }
   }

   // make a pointer to that stuff
   command_line_data *cld_p = new command_line_data(cld);
   g_object_set_data(G_OBJECT(app), "command-line-data", cld_p);

   // 20230515-PE Is this really what I want to do?
   // This seems like a bit of a hack.
   // Perhaps put the contentx of new_startup_application_activate() is a new function
   // That both this function and new_startup_application_activate() call.
   //
   new_startup_application_activate(app, user_data);

}

typedef struct {
  gboolean switch_option;
} AppOptions;
static AppOptions app_options;

void
command_line_stuff(GApplication *app, AppOptions *options) {

   const GOptionEntry cmd_params[] =
  {
    {
      .long_name = "my_switch_option",
      .short_name = 'm',
      .flags = G_OPTION_FLAG_NONE,     // see `GOptionFlags`
      .arg = G_OPTION_ARG_NONE,        // type of option (see `GOptionArg`)
      .arg_data = &(options->switch_option),// store data here
      .description = "<my description>",
      .arg_description = NULL,
    },
    {NULL}
  };

  g_application_add_main_option_entries(G_APPLICATION (app), cmd_params);
}


void application_command_line_callback(GtkApplication *app, GVariant *parameters, gpointer user_data) {

#if 0
   GVariantIter iter;
   GVariant *argument;
   gchar *key;
   gsize length;
   g_variant_iter_init(&iter, arguments);
   while (g_variant_iter_next(&iter, "{sv}", &key, &argument)) {
      std::string ss = g_variant_get_string(argument, &length);
      std::cout << "command line argument: " << key << " " << ss << std::endl;
      g_variant_unref(argument);
   };
   g_variant_unref(arguments);
#endif

   GVariant *argument;
   GVariantIter iter;
   const char *arg;

   return;

   /* Convert the command line arguments to a GVariant */
   // variant = g_variant_new_strv((const gchar * const *)argv, argc);

   // variant = arguments;

   std::cout << "Here A " << parameters << std::endl;
   /* Create an iterator for the GVariant */
   // iter = g_variant_iter_new(parameters);

   g_variant_iter_init(&iter, parameters);

   /* Loop over the arguments and print them */
   std::cout << "Here B " << &iter << std::endl;
   while (g_variant_iter_next(&iter, "{sv}", &arg, &argument)) {
      std::cout << "Here C " << &iter << std::endl;
      g_print("Argument: %s\n", arg);
   }

   /* Free the iterator and GVariant */
   // g_variant_iter_free(iter);

   // g_variant_unref(variant);
  
}


int new_startup(int argc, char **argv) {

#ifdef USE_LIBCURL
   curl_global_init(CURL_GLOBAL_NOTHING); // nothing extra (e.g. ssl or WIN32)
#endif

   mmdb::InitMatType();

   // setup_symm_lib();
   // check_reference_structures_dir();

   gtk_init();

   load_css();

   GtkWidget *splash_screen = new_startup_create_splash_screen_window();
   gtk_widget_show(splash_screen);

   g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);

   GError *error = NULL;
   // GtkApplication *app = gtk_application_new ("org.emsley.coot", G_APPLICATION_HANDLES_COMMAND_LINE);
   GtkApplication *app = gtk_application_new ("org.emsley.coot", G_APPLICATION_HANDLES_OPEN);
   g_application_register(G_APPLICATION(app), NULL, &error);
   // g_application_set_flags(G_APPLICATION(app), G_APPLICATION_HANDLES_COMMAND_LINE);


   // command_line_stuff(G_APPLICATION(app), &app_options);

   // g_signal_connect(app, "command-line", G_CALLBACK(application_command_line_callback), nullptr);

   application_activate_data *activate_data = new application_activate_data(argc,argv);
   activate_data->splash_screen = splash_screen;
   g_signal_connect(app, "open",     G_CALLBACK(application_open_callback), activate_data); // passed on
   // this destroys active_data
   g_signal_connect(app, "activate", G_CALLBACK(new_startup_application_activate), activate_data);

   // delete activate_data; Nope. This is used in new_startup_application_activate.
   // Delete it there if you want to delete it.

   int status = g_application_run(G_APPLICATION(app), argc, argv);
   std::cout << "--- g_application_run() returns with status " << status << std::endl;
   g_object_unref (app);
   return status;
}
