/* src/coot-tw.hh
 * 
 * Copyright 2004, 2005 by Paul Emsley, The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

#ifdef HAVE_GTK_CANVAS
#include <gtk/gtk.h>
#include <gdk_imlib.h>
#include <gtk-canvas.h>

#else 
   #ifdef HAVE_GNOME_CANVAS
      #include <gtk/gtk.h>
      #include <libgnomecanvas/libgnomecanvas.h>
      #include <gdk-pixbuf/gdk-pixbuf.h>
      typedef GnomeCanvas     GtkCanvas;
      typedef GnomeCanvasItem GtkCanvasItem;
   #endif
#endif


#include <vector>
#include <string>


GtkWidget *
tw_lookup_widget (GtkWidget       *widget,
		  const gchar     *widget_name);
   
void
on_tw_close_button_clicked     (GtkButton *button,
				gpointer         user_data);

namespace coot {

   class tw {
      std::vector<GtkCanvasItem *> canvas_item_vec;
      GtkCanvas *create_canvas(); // and store it.
      GtkWidget *create_tw_dialog();
      void add_elements_to_canvas();
      int  add_image_to_canvas(const std::string &filename, int position_index);
      static void free_imlib_image (GtkObject *object, gpointer data);
   public:
      tw() { canvas = 0; /* NULL */ }
      GtkCanvas *canvas;
      GtkCanvas *Canvas() const {return canvas;}
      GtkWidget *create_filled_toplevel();

      // item clicked callback
      static gint on_tw_canvas_item_clicked(GtkCanvasItem *item, GdkEvent *event,
					    gpointer data);

      static gint tw_motion_notify (GtkWidget *widget, GdkEventMotion *event); 

   };
}

#endif // HAVE_GTK_CANVAS or HAVE_GNOME_CANVAS

