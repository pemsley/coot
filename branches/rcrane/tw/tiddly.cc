
#ifdef HAVE_GTK_CANVAS
#include <gtk/gtk.h>
#include <gdk_imlib.h>
#include <gtk-canvas.h>

#else

#ifdef HAVE_GNOME_CANVAS
#include <gtk/gtk.h>
#include <libgnomecanvas/libgnomecanvas.h>

#endif
#endif



#include "coot-tw.hh"

int main(int argc, char **argv) {

#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   
   gtk_init(&argc, &argv);

   coot::tw tiddly;
   GtkWidget *w = tiddly.create_filled_toplevel();

   gtk_widget_show(w);
   gtk_main();
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   return 0; 
}

