
#include "graphics-info.h"
#include "dynamic-menus.hh"
#include "startup-utils.hh"
#include "draw-2.hh" // for gl area widget
#include "init-from-gtkbuilder.hh"

// return success status
bool init_from_gtkbuilder(GtkWidget *app_window) { // i.e. main window

   // get the right file first...

   bool status = true;

   std::string dir = coot::package_data_dir();
   std::string dir_glade = coot::util::append_dir_dir(dir, "glade");
   std::string ui_file_name = "coot-gtk4.ui";
   // std::string ui_file_name = "test-fragment.ui";
   std::string ui_file_full = coot::util::append_dir_file(dir_glade, ui_file_name);

   // local directory override
   if (coot::file_exists(ui_file_name))
      ui_file_full = ui_file_name;

   const char *env = getenv("COOT_GLADE");
   if (env)
      ui_file_full = std::string(env);

   GtkBuilder *builder = gtk_builder_new();

   GError* error = NULL;
   gboolean add_from_file_status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);
   if (add_from_file_status == FALSE) {
      std::cout << "DEBUG:: init_from_gtkbuilder(): glade file: " << ui_file_full
                << " add_from_file_status: " << add_from_file_status << std::endl;
      std::cout << "ERROR:: Failure to read or parse " << ui_file_full << std::endl;
      std::cout << "ERROR:: " << error->message << std::endl;
      exit(1);
   }

   GtkWidget *graphics_hbox = GTK_WIDGET(gtk_builder_get_object(builder, "main_window_graphics_hbox"));

   // 20220310-PE and the preferences builder too now

   GtkBuilder *preferences_builder = get_builder_for_preferences_dialog();
   graphics_info_t::set_preferences_gtkbuilder(preferences_builder);

   if (graphics_hbox) {

      graphics_info_t::set_gtkbuilder(builder); // store for future widget queries

      // GtkWidget *main_window = GTK_WIDGET(gtk_builder_get_object(builder, "main_window"));
      GtkWidget *main_window = app_window;
      GtkWidget *sb          = GTK_WIDGET(gtk_builder_get_object(builder, "main_window_statusbar"));
      GtkWidget *main_window_deletable_label = GTK_WIDGET(gtk_builder_get_object(builder, "main_window_deletable_label"));
      graphics_info_t::statusbar = sb;
      graphics_info_t::add_status_bar_text("Locked and loaded.");

      if (main_window_deletable_label) // it might not be looked up correctly when testing
         gtk_widget_hide(main_window_deletable_label); // 20220531-PE GTK4: can't delete it.

      if (true) {
         std::cout << "DEBUG:: main_window:   " << main_window << std::endl;
         std::cout << "DEBUG:: graphics_hbox: " << graphics_hbox << std::endl;
         std::cout << "DEBUG:: statusbar:     " << sb << std::endl;
      }

      // 20220612-PE is this needed now?
      if (main_window)
         graphics_info_t::set_main_window(main_window);

      std::string main_title = make_main_window_title();
      gtk_window_set_title(GTK_WINDOW(main_window), main_title.c_str());

      create_dynamic_menus(main_window);

      GtkWidget *glarea = create_gtkglarea_widget();
      if (glarea) {
         graphics_info_t::glareas.push_back(glarea);

         gtk_box_append(GTK_BOX(graphics_hbox), glarea);
         GError *err = gtk_gl_area_get_error(GTK_GL_AREA(glarea));
         if (err)
            std::cout << "ERROR:: GL error in init_from_gtkbuilder()" << err << std::endl;

         // gtk_widget_show(main_window);
         // gtk_widget_realize(glarea); // Don't do this explicitly

      } else {
         std::cout << "WARNING:: init_from_gtkbuilder(): glarea null" << std::endl;
         status = false;
      }
   } else {
      std::cout << "WARNING:: init_from_gtkbuilder(): graphics_hbox was null" << std::endl;
      status = false;
   }
   return status;
}
