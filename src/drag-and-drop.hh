/* src/drag-and-drop.hh
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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

#ifndef DRAG_AND_DROP_HH
#define DRAG_AND_DROP_HH

#include <string>
#include <gtk/gtk.h>

enum { TARGET_STRING, TEXT_URL, TEXT_URI };

gboolean
on_gl_canvas_drag_drop(GtkWidget *widget,
		       GdkDragContext *context,
		       gint x, gint y,
		       guint time,
		       gpointer user_data);


void
on_drag_data_received (GtkWidget *widget, 
		       GdkDragContext *context, 
		       gint x, gint y,
		       GtkSelectionData *selection_data, 
		       guint target_type, 
		       guint time,
		       gpointer data);

int handle_drag_and_drop_single_item(const std::string &file_name);


#endif // DRAG_AND_DROP_HH

