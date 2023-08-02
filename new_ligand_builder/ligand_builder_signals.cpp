#include <gtk/gtk.h>
#include "ligand_builder_state.hpp"

using namespace coot::ligand_editor;

extern "C" G_MODULE_EXPORT
void
layla_on_close(GtkButton* button, gpointer user_data) {
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");

    state->file_exit();
}

extern "C" G_MODULE_EXPORT
void
layla_on_apply(GtkButton* button, gpointer user_data) {
    g_warning("TODO: Implement 'Apply'");
}
