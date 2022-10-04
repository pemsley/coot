#include "setup-gui-components.hh"

#include <gtk/gtk.h>
#include "graphics-info.h"

gboolean
generic_hide_on_escape_controller_cb(
      GtkEventControllerKey  *controller,
      guint                  keyval,
      guint                  keycode,
      GdkModifierType        modifiers,
      GtkWidget              *to_be_hidden) {
   gboolean handled = TRUE;
   switch (keyval) {
      case GDK_KEY_Escape: {
         gtk_widget_hide(to_be_hidden);
         break;
      }
      default: {
         g_debug("generic_hide_on_escape_controller_cb: unhandled key: %s",gdk_keyval_name(keyval));
         handled = FALSE;
         break;
      }
   }
   return gboolean(handled);
}

void setup_generic_hide_on_escape_controller(GtkWidget* target_widget, GtkWidget* to_be_hidden) {
   GtkEventController *key_controller = gtk_event_controller_key_new();
   g_debug("Setting up hide-on-Escape controller on %p to hide %p",target_widget,to_be_hidden);
   g_signal_connect(key_controller, "key-pressed",G_CALLBACK(generic_hide_on_escape_controller_cb), to_be_hidden);
   gtk_widget_add_controller(target_widget, key_controller);
}

void setup_fetch_pdb_using_code() {
   int n = COOT_ACCESSION_CODE_WINDOW_OCA;
   GtkWidget *frame = widget_from_builder("accession_code_frame");
   //GtkWidget *label = widget_from_builder("accession_code_label");
   //gtk_label_set_text(GTK_LABEL(label), "PDB Accession Code: ");
   g_object_set_data(G_OBJECT(frame), "mode", GINT_TO_POINTER(n));
   setup_generic_hide_on_escape_controller(widget_from_builder("accession_code_entry"),frame);
}

void setup_get_monomer() {
   GtkWidget* frame = widget_from_builder("get_monomer_frame");
   GtkWidget* entry = widget_from_builder("get_monomer_entry");
   setup_generic_hide_on_escape_controller(entry,frame);
}

void setup_gui_components() {
   g_info("Initializing UI components...");
   setup_get_monomer();
   setup_fetch_pdb_using_code();
   g_info("Done initializing UI components.");
}
