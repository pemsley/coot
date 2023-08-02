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

extern "C" G_MODULE_EXPORT
void
layla_on_invalid_molecule_toggled(GtkCheckButton* check_button, gpointer user_data) {
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    coot_ligand_editor_set_allow_invalid_molecules(canvas, gtk_check_button_get_active(check_button));
}

extern "C" G_MODULE_EXPORT
void
layla_on_show_alerts_toggled(GtkCheckButton* check_button, gpointer user_data) {
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    g_warning("TODO: Implement 'Show Alerts'");
}

extern "C" G_MODULE_EXPORT
void
layla_on_scale_spinbutton_value_changed(GtkSpinButton* self,gpointer user_data){
    LigandBuilderState* state = (LigandBuilderState*)g_object_get_data(G_OBJECT(user_data), "ligand_builder_instance");
    CootLigandEditorCanvas* canvas = state->get_canvas();
    double new_scale = gtk_spin_button_get_value(self);
    // This should prevent infinite cascade of signals being emited
    if (coot_ligand_editor_get_scale(canvas) != new_scale) {
        coot_ligand_editor_set_scale(canvas, new_scale);
    }
}