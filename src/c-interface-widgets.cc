/* src/c-interface-build.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by the University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <string>
#include <vector>

#include <gtk/gtk.h>
#include "c-interface.h"

#include "graphics-info.h"

void do_rot_trans_adjustments(GtkWidget *dialog) { 
   graphics_info_t g;
   g.do_rot_trans_adjustments(dialog);
}

short int delete_item_widget_is_being_shown() {
   short int r = 0; 
   if (graphics_info_t::delete_item_widget != NULL) {
      r = 1;
   }
   return r;
}

short int delete_item_widget_keep_active_on() {
   short int r = 0;
   if (delete_item_widget_is_being_shown()) { 
      GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
					     "delete_item_keep_active_checkbutton");
      if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	 r = 1;
      }
   }
}

void store_delete_item_widget_position() {

   gint upositionx, upositiony;
   gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
			       &upositionx, &upositiony);
   graphics_info_t::delete_item_widget_x_position = upositionx;
   graphics_info_t::delete_item_widget_y_position = upositiony;
   gtk_widget_destroy(graphics_info_t::delete_item_widget);
   clear_delete_item_widget();
}

void clear_delete_item_widget() {

   graphics_info_t::delete_item_widget = NULL;
}


void store_delete_item_widget(GtkWidget *widget) {
   graphics_info_t::delete_item_widget = widget;
}
