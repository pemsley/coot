
#include <iostream>
#include <string>
#include <gtk/gtk.h>
#include <epoxy/gl.h>
extern "C" { void load_tutorial_model_and_data(); }

GtkWidget *widget_from_builder(const std::string &w_name, GtkBuilder *builder) {
   GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(GTK_BUILDER(builder), w_name.c_str()));
   return  w;
}

void print_opengl_info();

#include "graphics-info.h"

void init_framebuffers(GtkWidget *glarea) {

   // put this into graphics-info I suppose.

   std::cout << "DEBUG:: use_framebuffers: " << graphics_info_t::use_framebuffers << std::endl;

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
   int w = allocation.width;
   int h = allocation.height;

      if (graphics_info_t::use_framebuffers) {
         unsigned int index_offset = 0;
         GLenum err;
         graphics_info_t::screen_framebuffer.init(w, h, index_offset, "screen/occlusion");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post screen_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 1;
         graphics_info_t::blur_y_framebuffer.init(w, h, index_offset, "blur-y");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_y_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 2;
         graphics_info_t::blur_x_framebuffer.init(w, h, index_offset, "blur-x");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_x_framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 3;
         graphics_info_t::combine_textures_using_depth_framebuffer.init(w, h, index_offset, "new-blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_combine framebuffer init() err is "
                                                << err << std::endl;
         index_offset = 4;
         graphics_info_t::blur_framebuffer.init(w, h, index_offset, "blur");
         err = glGetError(); if (err) std::cout << "start on_glarea_realize() post blur_framebuffer init() err is "
                                                << err << std::endl;
      }

}



void
new_startup_realize(GtkWidget *gl_area) {

   std::cout << "new_startup_realize() ------------------- start ------------------"
             << std::endl;

   GdkGLContext *context;
   gtk_gl_area_make_current(GTK_GL_AREA (gl_area));

   if (gtk_gl_area_get_error(GTK_GL_AREA (gl_area)) != NULL)
      return;

   context = gtk_gl_area_get_context(GTK_GL_AREA(gl_area));

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(gl_area), TRUE);

   print_opengl_info();

   graphics_info_t g;
   init_framebuffers(gl_area); // Hmm - I don't know what this does compared to below.
   g.init_framebuffers();
   g.init_buffers();
   g.init_joey_ssao_stuff();
   g.init_shaders();
   g.setup_lights();

   g.setup_key_bindings();

   GLenum err = glGetError();
   if (err)
      std::cout << "ERROR:: new_startup_realize() --start-- err is " << err << std::endl;
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


void
new_startup_on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   std::cout << "DEBUG: new_startup_on_glarea_resize() " <<  width << " " << height << std::endl;
   graphics_info_t g;
   // for the GL widget, not the window.
   g.graphics_x_size = width;
   g.graphics_y_size = height;
   // g.reset_frame_buffers(width, height); // currently makes the widget blank (not drawn)
   g.init_shaders();
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

void
on_open_clicked(GSimpleAction *action,
                GVariant *parameter,
                gpointer data) {

   std::cout << "open clicked" << std::endl;

}

void
on_close_clicked(GSimpleAction *action,
                 GVariant *parameter,
                 gpointer data) {

   std::cout << "close clicked" << std::endl;

}


GMenu *create_menu_by_hand(const GtkApplication *application) {
   const GActionEntry entries[] = {
      { "open",  on_open_clicked,  NULL, NULL, NULL, { 0, 0, 0 } },
      { "close", on_close_clicked, NULL, NULL, NULL, { 0, 0, 0 } }
   };

   g_action_map_add_action_entries(G_ACTION_MAP(application), entries, G_N_ELEMENTS(entries), NULL);

   GMenu *menu      = g_menu_new();
   GMenu *file_menu = g_menu_new();

   GMenuItem *item;
   item = g_menu_item_new("Open", "app.open");
   g_menu_append_item(file_menu, item);

   item = g_menu_item_new("Close", "app.close");
   g_menu_append_item(file_menu, item);

   g_menu_append_submenu(menu, "File", G_MENU_MODEL(file_menu));

   return menu;
}


void on_glarea_drag_begin_primary(GtkGestureDrag *gesture,
                          double          x,
                          double          y,
                          GtkWidget      *area) {

   // display_info_t di;
   // di.mouse_x = x;
   // di.mouse_y = y;
   // di.drag_begin_x = x;
   // di.drag_begin_y = y;

   graphics_info_t g;
   g.on_glarea_drag_begin_primary(gesture, x, y, area);

}

void on_glarea_drag_update_primary(GtkGestureDrag *gesture,
                           double          delta_x,
                           double          delta_y,
                           GtkWidget      *area) {

   graphics_info_t g;
   // g.on_glarea_drag_update_primary(gesture, delta_x, delta_y, area);

   g.on_glarea_drag_update_secondary(gesture, delta_x, delta_y, area);

}

void on_glarea_drag_end_primary(GtkGestureDrag *gesture,
                                double          x,
                                double          y,
                                GtkWidget      *area) {

   // std::cout << "drag end" << std::endl;
   // do nothing at the moment.
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

}

void
on_glarea_scrolled(GtkEventControllerScroll *controller,
                   double                    dx,
                   double                    dy,
                   gpointer                  user_data) {

   graphics_info_t g;
   g.on_glarea_scrolled(controller, dx, dy, user_data);

}


void setup_gestures(GtkWidget *glarea) {

      std::cout << "================= setting up GTK4 style event controlllers ====================" << std::endl;

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

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER (drag_controller_primary));
      g_signal_connect(drag_controller_primary, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_primary),  glarea);
      g_signal_connect(drag_controller_primary, "drag-update", G_CALLBACK(on_glarea_drag_update_primary), glarea);
      g_signal_connect(drag_controller_primary, "drag-end",    G_CALLBACK(on_glarea_drag_end_primary),    glarea);

      gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_secondary), GDK_BUTTON_SECONDARY);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER (drag_controller_secondary));
      g_signal_connect(drag_controller_secondary, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_secondary),  glarea);
      g_signal_connect(drag_controller_secondary, "drag-update", G_CALLBACK(on_glarea_drag_update_secondary), glarea);
      g_signal_connect(drag_controller_secondary, "drag-end",    G_CALLBACK(on_glarea_drag_end_secondary),    glarea);

      gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_middle), GDK_BUTTON_MIDDLE);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER (drag_controller_middle));
      g_signal_connect(drag_controller_middle, "drag-begin",  G_CALLBACK(on_glarea_drag_begin_middle),  glarea);
      g_signal_connect(drag_controller_middle, "drag-update", G_CALLBACK(on_glarea_drag_update_middle), glarea);
      g_signal_connect(drag_controller_middle, "drag-end",    G_CALLBACK(on_glarea_drag_end_middle),    glarea);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(click_controller));
      g_signal_connect(click_controller, "pressed",  G_CALLBACK(on_glarea_click),  glarea);

      gtk_widget_add_controller(GTK_WIDGET(glarea), GTK_EVENT_CONTROLLER(scroll_controller));
      g_signal_connect(scroll_controller, "scroll",  G_CALLBACK(on_glarea_scrolled),  glarea);

}

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"

void on_coords_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file   = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);

#if 0
      GSList *files_list = gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(dialog));
      while (files_list) {

         const char *fnc = static_cast<const char *>(files_list->data);
         if (fnc) {
            std::string fn(fnc);
            handle_read_draw_molecule_with_recentre(fn, 0);
         }
         files_list = g_slist_next(files_list);
      }
#endif

      GtkWidget *recentre_combobox = widget_from_builder("coords_filechooserdialog_recentre_combobox");
      int active_item_index = gtk_combo_box_get_active(GTK_COMBO_BOX(recentre_combobox));
      bool move_molecule_here_flag = false;
      bool recentre_on_read_pdb_flag = false;
      if (active_item_index == 0)
         recentre_on_read_pdb_flag = true;
      if (active_item_index == 1)
         recentre_on_read_pdb_flag = false;
      if (active_item_index == 2)
         move_molecule_here_flag = true;

      // open_file (file);
      if (file_name) {
         std::cout << "info: " << file_name << " " << move_molecule_here_flag << " " << recentre_on_read_pdb_flag
                   << std::endl;

         if (move_molecule_here_flag) {
            handle_read_draw_molecule_and_move_molecule_here(file_name);
         } else {
            if (recentre_on_read_pdb_flag)
               handle_read_draw_molecule_with_recentre(file_name, 1);
            else
               handle_read_draw_molecule_with_recentre(file_name, 0); // no recentre
         }
      }
   }
   gtk_window_close(GTK_WINDOW(dialog));
}


void open_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   // Ancient GTK3
   // open_coords_dialog();

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File",
                                                   parent_window,
                                                   action,
                                                   _("_Cancel"),
                                                   GTK_RESPONSE_CANCEL,
                                                   _("_Open"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);

   g_signal_connect(dialog, "response", G_CALLBACK(on_coords_filechooser_dialog_response_gtk4), NULL);
   gtk_widget_show(dialog);

}

void open_dataset_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *dataset_chooser = widget_from_builder("coords_filechooser_dialog");
   GtkWidget *main_window = graphics_info_t::get_main_window();
   gtk_window_set_transient_for(GTK_WINDOW(dataset_chooser), GTK_WINDOW(main_window));

   set_directory_for_filechooser(dataset_chooser);
   set_file_selection_dialog_size(dataset_chooser);
   add_filechooser_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);

}

void auto_open_mtz_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {


   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   int is_auto_read_fileselection = 1;
   set_directory_for_filechooser(dataset_chooser);
   add_filename_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   g_object_set_data(G_OBJECT(dataset_chooser), "imol", GINT_TO_POINTER(-1)); // 20220627-PE do I need this?
   g_object_set_data(G_OBJECT(dataset_chooser), "is_auto", GINT_TO_POINTER(is_auto_read_fileselection));
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);

}

void load_tutorial_model_and_data_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                         G_GNUC_UNUSED GVariant *parameter,
                                         G_GNUC_UNUSED gpointer user_data) {
   load_tutorial_model_and_data();
}

void exit_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   coot_checked_exit(0);
}

#include "curlew.h" // 20220628-PE why does this exist? why is curlew() not in the .hh file?

void curlew_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                   G_GNUC_UNUSED GVariant *parameter,
                   G_GNUC_UNUSED gpointer user_data) {
   curlew();
}

void
create_actions(GtkApplication *application) {

   GSimpleAction *simple_action;

   GtkWindow *application_window = gtk_application_get_active_window(application);

   simple_action = g_simple_action_new("open_coordinates_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(open_coordinates_action), application_window);

   simple_action = g_simple_action_new("auto_open_mtz_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(auto_open_mtz_action), application_window);

   simple_action = g_simple_action_new("open_dataset_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(open_dataset_action), application_window);

   simple_action = g_simple_action_new("load_tutorial_model_and_data_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(load_tutorial_model_and_data_action), NULL);

   simple_action = g_simple_action_new("curlew_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(curlew_action), NULL);

   simple_action = g_simple_action_new("exit_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(exit_action), NULL);
}

void
new_startup_application_activate(GtkApplication *application,
                                 gpointer        user_data) {

   GtkBuilder *builder = gtk_builder_new();
   if (GTK_IS_BUILDER(builder)) {
   } else {
      std::cout << "in new_startup_application_activate() builder was NOT a builder" << std::endl;
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

   GtkWidget *app_window = gtk_application_window_new(application);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_title(GTK_WINDOW(app_window), "Coot App Main Window");
   graphics_info_t::set_main_window(app_window);

   guint id = gtk_application_window_get_id(GTK_APPLICATION_WINDOW(app_window));
   std::cout << "debug:: new_startup_application_activate(): Window id: " << id << std::endl;

   graphics_info_t g;
   g.set_gtkbuilder(builder);

   //GMenu *menu = create_menu_by_hand(application);
   GMenu *menubar = G_MENU(g.get_gobject_from_builder("menubar"));
   gtk_application_set_menubar(application, G_MENU_MODEL(menubar));
   gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(app_window), TRUE);

   // GtkWidget *graphics_hbox = widget_from_builder("crows_graphics_hbox", builder);
   // GtkWidget *main_window   = widget_from_builder("crows_main_window",   builder);
   GtkWidget *graphics_hbox = widget_from_builder("main_window_hbox", builder);
   GtkWidget *graphics_vbox = widget_from_builder("main_window_vbox", builder);
   // GObject *menubar  = g.get_gobject_from_builder("main_window_menubar");

   gtk_window_set_child(GTK_WINDOW(app_window), graphics_vbox);

   gtk_window_present(GTK_WINDOW(app_window));
   // gtk_widget_show(window);

   std::cout << "debug:: new_startup_application_activate(): setting the menubar: " << menubar << std::endl;
   gtk_application_set_menubar(application, G_MENU_MODEL(menubar));
   gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(app_window), TRUE);

   GtkWidget *gl_area = new_startup_create_glarea_widget();
   graphics_info_t::glareas.push_back(gl_area);
   gtk_widget_show(gl_area);
   // gtk_box_prepend(GTK_BOX(graphics_hbox), gl_area); // crows
   gtk_box_prepend(GTK_BOX(graphics_hbox), gl_area);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_default_size(GTK_WINDOW(app_window), 300, 300);
   gtk_widget_set_size_request(gl_area, 700, 400); // bigger than the window size - for testing.
   gtk_widget_show(app_window);

   setup_gestures(gl_area);

   create_actions(application);

   // load_tutorial_model_and_data();

}

// move these to the top.
void setup_symm_lib();
void check_reference_structures_dir();

int new_startup(int argc, char **argv) {

   graphics_info_t graphics_info;
   setup_symm_lib();
   check_reference_structures_dir();
   graphics_info.init();
   gtk_init();

   // set this by parsing the command line arguments
   graphics_info.use_graphics_interface_flag = true;

   g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);

   GError *error = NULL;
   GtkApplication *app = gtk_application_new ("org.coot.crows", G_APPLICATION_FLAGS_NONE);
   g_signal_connect(app, "activate", G_CALLBACK(new_startup_application_activate), NULL);
   g_application_register(G_APPLICATION(app), NULL, &error);
   int status = g_application_run (G_APPLICATION (app), argc, argv);
   std::cout << "--- g_application_run() returns with status " << status << std::endl;
   g_object_unref (app);
   return status;
}
