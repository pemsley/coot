#include "about_dialog.hpp"


GtkWidget* coot::ligand_editor::build_about_dialog() noexcept {
    GtkAboutDialog* dialog = GTK_ABOUT_DIALOG(gtk_about_dialog_new());
    const char* authors[3];
    authors[0] = "Jakub Smulski";
    authors[1] = "Paul Emsley";
    authors[2] = nullptr;
    gtk_about_dialog_set_authors(dialog, authors);
    gtk_about_dialog_set_program_name(dialog, "Coot Layla");
    gtk_about_dialog_set_version(dialog, "v1.0");
    gtk_about_dialog_set_comments(dialog, "Ligand builder, successor to Lidia");
    gtk_about_dialog_set_license_type(dialog, GTK_LICENSE_GPL_3_0);
    gtk_about_dialog_set_copyright(dialog, "Who owns the copyrights now?");
    return GTK_WIDGET(dialog);
}