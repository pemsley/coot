#include <gtk/gtk.h>

extern "C" G_MODULE_EXPORT
void
layla_on_close(GtkButton* button, gpointer user_data) {
    // todo: this should probably do some checks before just closing
    gtk_window_close(GTK_WINDOW(user_data));
}