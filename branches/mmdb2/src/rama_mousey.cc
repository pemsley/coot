/* src/rama_mousey.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006 by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

#include <iostream>

#include "rama_mousey.hh" // has gtk/gtk.h
#include "interface.h"    // needs gtk/gtk.h

#include "rama_plot.hh"

// The motion callback was attached at the canvas, so widget is a
// canvas here.
// 
gint rama_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   double x,y;
   int x_as_int, y_as_int;
   GdkModifierType state;

   if (event->is_hint) {
      // std::cout << "was hint";  usually true
      gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
      x = x_as_int;
      y = y_as_int;
   } else {
      x = event->x;
      y = event->y;
      state = (GdkModifierType) event->state;
   }

   coot::rama_plot *plot =
      static_cast<coot::rama_plot *> (gtk_object_get_user_data(GTK_OBJECT(widget)));

   if (plot) { 
      // plot->map_mouse_pos(x,y); 
      plot->mouse_motion_notify(event,x,y);
   } else { 
      std::cout << "ERROR:: in getting user data from canvas\n";
   } 

   // we need to convert x, y to phi, psi space.
   // 
   // Now, x and y are in widget space, so thats 0->319 in x and y
   // right now.
   //
   // That needs to get mapped by the scaling factor which is in
   // setup_canvas().
   //
   // That brings us to widget space, 319 -> 255.
   // 
   return 0; 
}


gint rama_button_press (GtkWidget *widget, GdkEventButton *event) {
   coot::rama_plot *plot =
      static_cast<coot::rama_plot *> (gtk_object_get_user_data(GTK_OBJECT(widget)));

   plot->button_press(widget, event);

   return 0; 
   
}

gint rama_key_release_event(GtkWidget *widget, GdkEventKey *event) {

   coot::rama_plot *plot =
      static_cast<coot::rama_plot *> (gtk_object_get_user_data(GTK_OBJECT(widget)));
   gint i = plot->key_release_event(widget, event);
   return i; 
}

gint rama_key_press_event(GtkWidget *widget, GdkEventKey *event) {

   return 0; 

}

void rama_show_preferences() {
 
   //    GtkWidget *widget = create_propertybox1();
   GtkWidget *widget = create_dynarama_properties_window(); 
   gtk_widget_show(widget); 

}

void rama_zoom_out(GtkWidget *widget) {

   coot::rama_plot *plot =
      static_cast<coot::rama_plot *> (gtk_object_get_user_data(GTK_OBJECT(widget)));

   plot->zoom_out(); 

}

void rama_zoom_in(GtkWidget *widget) {

   coot::rama_plot *plot =
      (coot::rama_plot *) gtk_object_get_user_data(GTK_OBJECT(widget));

   plot->zoom_in(); 

}

#endif // HAVE_GTK_CANVAS or HAVE_GNOME_CANVAS
