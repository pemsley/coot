#include "setup-gui-components.hh"

#include <gtk/gtk.h>
#include "graphics-info.h"

void setup_fetch_pdb_using_code_action() {
   int n = COOT_ACCESSION_CODE_WINDOW_OCA;
   GtkWidget *frame = widget_from_builder("accession_code_frame");
   GtkWidget *label = widget_from_builder("accession_code_label");
   //gtk_label_set_text(GTK_LABEL(label), "PDB Accession Code: ");
   g_object_set_data(G_OBJECT(frame), "mode", GINT_TO_POINTER(n));

}

void setup_gui_components() {
   g_info("Initializing UI components...");
   setup_fetch_pdb_using_code_action();
   g_info("Done initializing UI components.");
}
