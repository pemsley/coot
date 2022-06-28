
#ifndef WIDGET_FROM_BUILDER_HH
#define WIDGET_FROM_BUILDER_HH

#include <string>
#include <gtk/gtk.h>

GtkWidget *widget_from_builder(const std::string &w_name, GtkBuilder *builder);

GtkWidget *widget_from_builder(const std::string &w_name);

GtkBuilder *get_builder_for_preferences_dialog();

GtkWidget *widget_from_preferences_builder(const std::string &w_name);

#endif // WIDGET_FROM_BUILDER_HH
