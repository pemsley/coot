
#define GIO_COMPILATION
#include "graphics-info.h"
#include "coot-application.hh"
#include "init-from-gtkbuilder.hh"
void setup_application_icon(GtkWindow *window);
#include "startup-utils.hh" // nothing yet

void
application_activate(GtkApplication *app,
                     gpointer        user_data) {


#if (GTK_MAJOR_VERSION == 4)

   auto do_window_resizing_widgets = [] () {

      GtkWidget *box = widget_from_builder("main_window_resize_window_button_box");
#ifdef __APPLE__
      // gtk_widget_show(box);
      // GtkWidget *window = widget_from_builder("main_window");
      // 20220407-PE this causes a crash
      // gtk_window_set_has_resize_grip(GTK_WINDOW(main_window), FALSE);
      // this expands the window fully in height - I don't want that.
      // gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
#else
      if (box)
         gtk_widget_hide(box);
#endif
   };

   GtkWidget *app_window = gtk_application_window_new(app);
   graphics_info_t::main_window = app_window;

   bool success = init_from_gtkbuilder(app_window);
   if (success) {

      GtkWidget *main_window_vbox = widget_from_builder("main_window_vbox");
      std::cout << "-------------------- found main_window_vbox " << main_window_vbox << std::endl;
      if (main_window_vbox) {
         std::cout << "-------------------- calling gtk_window_set_child " << app_window << std::endl;
         gtk_window_set_child(GTK_WINDOW(app_window), main_window_vbox);
         setup_application_icon(GTK_WINDOW(app_window)); // put this in init_from_gtkbuilder()
         std::cout << "-------------------- calling do_window_resizing_widgets() " << main_window_vbox << std::endl;
         do_window_resizing_widgets();
         gtk_widget_show(main_window_vbox);
         gtk_widget_show(app_window);
      }
   }

#endif
}


int start_using_application(int argc, char **argv) {

   int status = 0;
#if (GTK_MAJOR_VERSION >= 4)

   gtk_init();

   if (graphics_info_t::use_graphics_interface_flag) {
      GError *error = NULL;
      GtkApplication *app = gtk_application_new ("org.emsley.coot", G_APPLICATION_FLAGS_NONE);

      std::cout << "DEBUG:: startup is g_application " << G_IS_APPLICATION(app) << std::endl;
      std::cout << "DEBUG:: startup is g_application argc " << argc << std::endl;
      std::cout << "DEBUG:: startup is g_application argv " << argv << std::endl;

      g_signal_connect(app, "activate", G_CALLBACK(application_activate), NULL);
      gboolean register_status = g_application_register(G_APPLICATION(app), NULL, &error);
      std::cout << "g_application_register() return status " << register_status << std::endl;
      if (error)
         std::cout << "ERROR:: post-register error message " << error->message << std::endl;
      std::cout << ":::::::::::::::::::::::::::::::::: calling g_application_run()" << std::endl;
      status = g_application_run(G_APPLICATION (app), argc, argv);
      std::cout << "::::::::::::::::::::::::::::::::::::: done g_application_run()" << std::endl;
      std::cout << "---------------- g_application_run() returns " << status << std::endl;
      if (error)
         std::cout << "ERROR:: post run error message " << error->message << std::endl;
      g_object_unref(app);
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "------------ start_using_application() returns --------" << std::endl;
      std::cout << "-------------------------------------------------------" << std::endl;
      std::cout << "-------------------------------------------------------" << std::endl;
   }

#endif
   return status;
   

}
