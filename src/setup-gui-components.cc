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

void attach_css_style_class_to_overlays() {

   GtkCssProvider *provider = gtk_css_provider_new();
   gtk_css_provider_load_from_data (provider, ".mainWindowOverlayChild { background: rgba(0,0,0,0.7); }", -1);

   auto set_transparency_on_widget = [provider](GtkWidget* widget){
      GtkStyleContext *context = gtk_widget_get_style_context(widget);
      gtk_style_context_add_provider (context, GTK_STYLE_PROVIDER (provider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
      gtk_style_context_add_class (context, "mainWindowOverlayChild");
      g_debug("'mainWindowOverlayChild' CSS class set for: %p %s",widget,G_OBJECT_CLASS_NAME(widget));
   };

   GtkWidget* overlay = widget_from_builder("main_window_graphics_overlay");
   GtkWidget* to_skip = widget_from_builder("main_window_graphics_hbox");
   for(GtkWidget* child = gtk_widget_get_first_child(overlay); 
       child != nullptr; 
       child = gtk_widget_get_next_sibling(child)) {
      if(child != to_skip) 
         set_transparency_on_widget(child);
   }
}

gboolean
on_python_scripting_entry_key_pressed(GtkEventControllerKey *controller,
                                                      guint                  keyval,
                                                      guint                  keycode,
                                                      GdkModifierType        modifiers,
                                                      GtkEntry              *entry) {
   gboolean handled = TRUE;
   
   switch(keyval) {
      case GDK_KEY_Up: {
         const char *entry_txt = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (entry_txt) {
            std::string t = graphics_info_t::command_history.get_previous_command();
            gtk_editable_set_text(GTK_EDITABLE(entry), t.c_str());
            g_debug("Setting command entry text to '%s'",t.c_str());
         }
         break;
      }
      case GDK_KEY_Down: {
         std::string t = graphics_info_t::command_history.get_next_command();
         gtk_editable_set_text(GTK_EDITABLE(entry), t.c_str());
         g_debug("Setting command entry text to '%s'",t.c_str());
         break;
      }
      case GDK_KEY_Escape: {
          g_idle_add(+[](gpointer data)-> gboolean {
            GtkRevealer* revealer = GTK_REVEALER(widget_from_builder("python_scripting_revealer"));
            gtk_revealer_set_reveal_child(revealer,FALSE);
            return G_SOURCE_REMOVE;
         },NULL);
         break;
      }
      default: {
         handled = FALSE;
         g_debug("Python scripting entry: Unhandled key: %s",gdk_keyval_name(keyval));
      }
   }
   return gboolean(handled);
}


void on_python_scripting_entry_activated(GtkEntry* entry, gpointer user_data) {
   const char *entry_txt = gtk_editable_get_text(GTK_EDITABLE(entry));
   g_info("Running python command: '%s'",entry_txt);
   PyRun_SimpleString(entry_txt);

   // add a copy of the text to history
   graphics_info_t::command_history.add_to_history(std::string(entry_txt));
   // clear the entry
   gtk_editable_set_text(GTK_EDITABLE(entry), "");
}

void setup_python_scripting_entry() {
   GtkWidget *entry = widget_from_builder("python_scripting_entry");
   if(entry == NULL) {
      g_error("'python_scripting_entry' from builder is NULL");
      return;
   }
   GtkEventController *key_controller_entry = gtk_event_controller_key_new();

   // for 'Up' and 'Down' keys, i.e. history lookup
   // and for 'Esc' key to hide the revealer
   g_signal_connect(key_controller_entry, "key-pressed",
                    G_CALLBACK(on_python_scripting_entry_key_pressed), entry);

   // for executing Python commands
   g_signal_connect(entry, "activate",G_CALLBACK(on_python_scripting_entry_activated), entry);

   gtk_widget_add_controller(entry, key_controller_entry);
}

void set_vertical_toolbar_internal_alignment() {
   GtkWidget *toolbar = widget_from_builder("main_window_vbox_inner");
   for(GtkWidget* child = gtk_widget_get_first_child(toolbar); 
       child != nullptr; 
       child = gtk_widget_get_next_sibling(child)) {
         // No need to do this for plain buttons.
         if(!(GTK_IS_MENU_BUTTON(child)||GTK_IS_TOGGLE_BUTTON(child))) {
            g_debug("set_vertical_toolbar_internal_alignment: Skippping toolbar item %p of type %s.",child,G_OBJECT_TYPE_NAME(child));
            continue;
         }
         GtkWidget* target;
         if(GTK_IS_MENU_BUTTON(child)) {
            target = gtk_menu_button_get_child(GTK_MENU_BUTTON(child));
         } else {
            target = gtk_button_get_child(GTK_BUTTON(child));
         }
         if(!target) {
            g_debug("set_vertical_toolbar_internal_alignment: Skippping toolbar item %p of type %s because its' \"child\" property is not set.",
            child,G_OBJECT_TYPE_NAME(child));
            continue;
         }
         // This is a hack. The parent that we get isn't our button but an internal GtkBox 
         // which is designated for storing GtkButton's 'child' widget.
         //
         // Unfortunately, currently there seems to be no other way to set this.
         // Gtk4 removed the necessary APIs.
         GtkWidget* parent_widget = gtk_widget_get_parent(target);
         if(!GTK_IS_BOX(parent_widget)) {
            if(GTK_IS_BOX(target)) {
               g_warning("set_vertical_toolbar_internal_alignment: Toolbar item %p of type %s: "
               "The parent widget that wraps %s::child is not a GtkBox but a %s. "
               "%s::child however is a GtkBox. Attempt will be made to align it. It might not work.",
               child,G_OBJECT_TYPE_NAME(child),G_OBJECT_TYPE_NAME(child),G_OBJECT_TYPE_NAME(parent_widget),G_OBJECT_TYPE_NAME(child));
               parent_widget = target;
            } else {
               g_warning("set_vertical_toolbar_internal_alignment: Skippping toolbar item %p of type %s: "
               "The parent widget that wraps %s::child is not a GtkBox but a %s",
               child,G_OBJECT_TYPE_NAME(child),G_OBJECT_TYPE_NAME(child),G_OBJECT_TYPE_NAME(parent_widget));
               continue;
            }
         }
         g_info("set_vertical_toolbar_internal_alignment: Aligning toolbar item %p of type %s.",child,G_OBJECT_TYPE_NAME(child));
         gtk_widget_set_halign(parent_widget, GTK_ALIGN_START);
   }
}

void setup_gui_components() {
   g_info("Initializing UI components...");
   setup_menubuttons();
   setup_get_monomer();
   setup_accession_code_frame();
   setup_python_scripting_entry();
   attach_css_style_class_to_overlays();
   set_vertical_toolbar_internal_alignment();
   g_info("Done initializing UI components.");
}
