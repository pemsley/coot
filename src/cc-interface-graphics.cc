
#include "cc-interface.hh"

#include "widget-from-builder.hh"

#include "graphics-info.h"

// 20211019-PE these have moved into graphics_info_t now because I want to add a
// key-binding "Esc" to do an unfullscreen()

void fullscreen() {
   graphics_info_t::fullscreen();

   GtkWidget *vbox       = widget_from_builder("main_window_vbox");
   GtkWidget *overlay    = widget_from_builder("main_window_graphics_overlay");
   GtkWidget *status_bar = widget_from_builder("main_window_statusbar");
   GtkWidget *tool_bar   = widget_from_builder("main_window_toolbar");
   GtkWidget *menu_bar   = widget_from_builder("main_window_menubar");

   gtk_container_remove(GTK_CONTAINER(vbox), status_bar);
   gtk_overlay_add_overlay(GTK_OVERLAY(overlay), status_bar);
   // gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), status_bar, TRUE);
   gtk_widget_set_halign(status_bar, GTK_ALIGN_START);
   gtk_widget_set_valign(status_bar, GTK_ALIGN_END);

   gtk_container_remove(GTK_CONTAINER(vbox), tool_bar);
   gtk_overlay_add_overlay(GTK_OVERLAY(overlay), tool_bar);
   // gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), tool_bar, TRUE);
   gtk_widget_set_halign(tool_bar, GTK_ALIGN_START);
   gtk_widget_set_valign(tool_bar, GTK_ALIGN_START);
   
   // gtk_container_remove(GTK_CONTAINER(vbox), menu_bar);
   // gtk_overlay_add_overlay(GTK_OVERLAY(overlay), menu_bar);
   // gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), menu_bar, TRUE);
   // gtk_widget_set_halign(menu_bar, GTK_ALIGN_START);
   // gtk_widget_set_valign(menu_bar, GTK_ALIGN_START);
   

}

void unfullscreen() {
   graphics_info_t::unfullscreen();

   // unoverlay and reparent here.

   GtkWidget *vbox       = widget_from_builder("main_window_vbox");
   GtkWidget *overlay    = widget_from_builder("main_window_graphics_overlay");
   GtkWidget *status_bar = widget_from_builder("main_window_statusbar");
   GtkWidget *tool_bar   = widget_from_builder("main_window_toolbar");
   GtkWidget *menu_bar   = widget_from_builder("main_window_menubar");
#if 0
   gtk_overlay_remove_overlay(GTK_OVERLAY(overlay), tool_bar);
   gtk_container_add(GTK_CONTAINER(vbox), tool_bar);

   gtk_overlay_remove_overlay(GTK_OVERLAY(overlay), status_bar);
   gtk_container_add(GTK_CONTAINER(vbox), status_bar);
#endif
}
