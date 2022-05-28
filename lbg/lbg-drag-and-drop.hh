/* lbg/lbg-drag-and-drop.cc
 * 
 * Copyright 2010 by the University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifdef HAVE_GOOCANVAS

enum { TARGET_STRING };

#include <string>
#include <gtk/gtk.h>

extern "C" G_MODULE_EXPORT gboolean
on_lbg_drag_drop (GtkWidget *widget,
	      GdkDragContext *context,
	      gint x, gint y,
	      guint time,
	      gpointer user_data);


extern "C" G_MODULE_EXPORT void
on_lbg_drag_data_received (GtkWidget *widget, 
		       GdkDragContext *context, 
		       gint x, gint y,
		       GtkSelectionData *selection_data, 
		       guint target_type, 
		       guint time,
		       gpointer data);
#endif 
