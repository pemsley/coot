#include "setup-gui-components.hh"

#include <gtk/gtk.h>
#include "graphics-info.h"
#include "c-interface-gtk-widgets.h"

// this function is both defined and implemented here.
// No other files should ever need it.
inline GMenuModel* menu_model_from_builder(const std::string& m_name) {
   GMenuModel *m = G_MENU_MODEL(graphics_info_t::get_gobject_from_builder(m_name));
   return m;
}

void setup_menubuttons() {
   GtkWidget* add_module_menubutton = widget_from_builder("add_module_menubutton");
   GMenuModel *modules_menu = menu_model_from_builder("modules-menu");
   gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(add_module_menubutton), modules_menu);
   
   // toolbar button - connect the refine menu to the GtkMenuButton
   GtkWidget *refine_menubutton = widget_from_builder("refine_menubutton");
   GMenuModel *refine_menu = menu_model_from_builder("refine-menu");
   gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(refine_menubutton), refine_menu);

   GtkWidget *fixed_atoms_menubutton = widget_from_builder("fixed_atoms_menubutton");
   GMenuModel *fixed_atoms_menu = menu_model_from_builder("fixed-atoms-menu");
   gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(fixed_atoms_menubutton), fixed_atoms_menu);

   GtkWidget *delete_menubutton = widget_from_builder("delete_menubutton");
   GMenuModel *delete_item_menu = menu_model_from_builder("delete-item-menu");
   gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(delete_menubutton), delete_item_menu);

   // move this function to where it can be called when we click on the "Mutate"
   // button (both of them, I suppose).
   auto add_typed_menu_to_mutate_menubutton = [] (const std::string &residue_type) {
      if (residue_type == "PROTEIN") {
         GtkWidget *mutate_menubutton = widget_from_builder("simple_mutate_menubutton");
         GMenuModel *mutate_menu = menu_model_from_builder("mutate-protein-menu");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);

         mutate_menubutton = widget_from_builder("mutate_and_autofit_menubutton");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);
      }
      if (residue_type == "NUCLEIC-ACID") {
         GtkWidget *mutate_menubutton = widget_from_builder("simple_mutate_menubutton");
         GMenuModel *mutate_menu = menu_model_from_builder("mutate-nucleic-acid-menu");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);

         mutate_menubutton = widget_from_builder("mutate_and_autofit_menubutton");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);
      }
   };


   add_typed_menu_to_mutate_menubutton("PROTEIN");
}

gboolean generic_hide_on_escape_controller_cb(
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

void setup_accession_code_frame() {
   GtkWidget *frame = widget_from_builder("accession_code_frame");
   GtkWidget* entry = widget_from_builder("accession_code_entry");
   g_signal_connect(entry,"activate",G_CALLBACK(+[](GtkEntry* entry, gpointer user_data){
      GtkWidget* frame = GTK_WIDGET(user_data);
      handle_get_accession_code(frame, GTK_WIDGET(entry));
   }),frame);
   setup_generic_hide_on_escape_controller(entry,frame);
}

void setup_get_monomer() {
   GtkWidget* frame = widget_from_builder("get_monomer_frame");
   GtkWidget* entry = widget_from_builder("get_monomer_entry");
   g_signal_connect(entry,"activate",G_CALLBACK(+[](GtkEntry* entry, gpointer user_data){
      handle_get_monomer_code(GTK_WIDGET(entry));
   }),NULL);
   setup_generic_hide_on_escape_controller(entry,frame);
}

void setup_gui_components() {
   g_info("Initializing UI components...");
   setup_menubuttons();
   setup_get_monomer();
   setup_accession_code_frame();
   g_info("Done initializing UI components.");
}
