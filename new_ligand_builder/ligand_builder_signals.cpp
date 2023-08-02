#include <gtk/gtk.h>
#include "ligand_builder_state.hpp"

extern "C" G_MODULE_EXPORT
void
layla_on_close(GtkButton* button, gpointer user_data) {
    // todo: this should probably do some checks before just closing
    gtk_window_close(GTK_WINDOW(user_data));
}

extern "C" G_MODULE_EXPORT
void
layla_on_apply(GtkButton* button, gpointer user_data) {
    g_warning("TODO: Implement 'Apply'");
}
