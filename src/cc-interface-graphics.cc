
#include "cc-interface.hh"

#include "widget-from-builder.hh"

void fullscreen() {

   GtkWidget *window = widget_from_builder("main_window");

   if (GTK_IS_WINDOW(window)) {
      gtk_window_fullscreen(GTK_WINDOW(window));
   }
                                                           

}

void unfullscreen() {

   GtkWidget *window = widget_from_builder("main_window");

   if (GTK_IS_WINDOW(window)) {
      gtk_window_unfullscreen(GTK_WINDOW(window));
   }

}
