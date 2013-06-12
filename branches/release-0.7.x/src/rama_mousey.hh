

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

#include <gtk/gtk.h>
 
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS
#endif

BEGIN_C_DECLS

// Callback Functions from declared in setup_canvas
// (i.e. can be c++ functions)
// 
gint rama_button_press (GtkWidget *widget, GdkEventButton *event);
gint rama_motion_notify(GtkWidget *widget, GdkEventMotion *event);
gint rama_key_press_event(GtkWidget *widget, GdkEventKey *event); // not used.
gint rama_key_release_event(GtkWidget *widget, GdkEventKey *event); 


// Extern C Functions called from callbacks.c
//
// (i.e.) must be extern c functions
void rama_show_preferences();
void rama_zoom_in(GtkWidget *widget);
void rama_zoom_out(GtkWidget *widget); 

END_C_DECLS

#endif // HAVE_GTK_CANVAS or HAVE_GNOME_CANVAS
