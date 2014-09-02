
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
