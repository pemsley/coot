#ifndef INIT_FROM_GTKBUILDER_HH
#define INIT_FROM_GTKBUILDER_HH

#include <gtk/gtk.h>

// we get the main_app_window from gtk_application_window_new() not from the builder file.
bool init_from_gtkbuilder(GtkWidget *main_app_window);

#endif // INIT_FROM_GTKBUILDER_HH
