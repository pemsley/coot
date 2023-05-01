//
#include "widget-from-builder.hh"
#include "graphics-info.h"

GtkWidget *widget_from_builder(const std::string &w_name, GtkBuilder *builder) {
   GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(GTK_BUILDER(builder), w_name.c_str()));
   return  w;
}

GtkWidget *widget_from_builder(const std::string &w_name) {

   GtkWidget *w = graphics_info_t::get_widget_from_builder(w_name);
   return w;
}

GtkWidget *widget_from_preferences_builder(const std::string &w_name) {

   GtkWidget *w = graphics_info_t::get_widget_from_preferences_builder(w_name);
   return w;
}
