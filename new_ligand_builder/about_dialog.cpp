#include "about_dialog.hpp"


GtkWidget* coot::ligand_editor::build_about_dialog() noexcept {
    GtkWidget* dialog = gtk_about_dialog_new();
    return dialog;
}