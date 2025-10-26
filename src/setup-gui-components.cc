/*
 * src/setup-gui-components.cc
 *
 * Copyright 2022 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#include <gtk/gtk.h>
#include "graphics-info.h"
#include "c-interface-gtk-widgets.h"
#include "setup-gui-components.hh"
#include "gtk/gtkshortcut.h"
#include "utils/coot-utils.hh"
#include "widget-from-builder.hh"

// this function is both defined and implemented here.
// No other files should ever need it.
inline GMenuModel* menu_model_from_builder(const std::string& m_name) {
   GMenuModel *m = G_MENU_MODEL(graphics_info_t::get_gobject_from_builder(m_name));
   return m;
}

// this function is both defined and implemented here.
// No other files should ever need it.
inline GMenu* menu_from_builder(const std::string& m_name) {
   GMenu *m = G_MENU(graphics_info_t::get_gobject_from_builder(m_name));
   return m;
}

// button (both of them, I suppose).
void
add_typed_menu_to_mutate_menubutton(const std::string &action_type, const std::string &residue_type) {

   // should I (do I need to) remove the menu model that is already attached to the menu button?

   if (action_type == "AUTOFIT") {
      GtkWidget *mutate_menubutton = widget_from_builder("mutate_and_autofit_menubutton");
      if (residue_type == "PROTEIN") {
         GMenuModel *mutate_menu = menu_model_from_builder("mutate-protein-menu");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);
      }
   }

   if (action_type == "SIMPLE") {
      GtkWidget *mutate_menubutton = widget_from_builder("simple_mutate_menubutton");
      if (residue_type == "PROTEIN") {
         GMenuModel *mutate_menu = menu_model_from_builder("mutate-protein-menu");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);
      }
      if (residue_type == "NUCLEIC-ACID") {
         GMenuModel *mutate_menu = menu_model_from_builder("mutate-nucleic-acid-menu");
         gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(mutate_menubutton), mutate_menu);
      }
   }
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

   GtkWidget *rotate_translate_button = widget_from_builder("rotate_translate_menubutton");
   GMenuModel *rotate_translate_menu = menu_model_from_builder("rotate-translate-menu");
   gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(rotate_translate_button), rotate_translate_menu);

   GtkWidget *rigid_body_button = widget_from_builder("vertical_toolbar_rigid_body_menubutton");
   GMenuModel *rigid_body_menu = menu_model_from_builder("rigid-body-fit-menu");
   gtk_menu_button_set_menu_model(GTK_MENU_BUTTON(rigid_body_button), rigid_body_menu);

   add_typed_menu_to_mutate_menubutton("AUTOFIT", "PROTEIN");
}

void setup_mutate_residue_range_dialog() {

   auto callback_func = +[] (GtkTextBuffer* buf, gpointer user_data) {
        std::cout << "on_mutate_molecule_sequence_text:buffer:changed --- start --- " << std::endl;
        GtkWidget *res_no_1_widget = widget_from_builder("mutate_molecule_resno_1_entry");
        GtkWidget *res_no_2_widget = widget_from_builder("mutate_molecule_resno_2_entry");
        GtkWidget *text_widget     = widget_from_builder("mutate_molecule_sequence_text");
        GtkWidget *label_widget    = widget_from_builder("mutate_residue_range_counts_label");
        mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
   };

   GtkWidget* mutate_molecule_sequence_text = widget_from_builder("mutate_molecule_sequence_text");
   GtkTextBuffer* buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(mutate_molecule_sequence_text));
   g_signal_connect(buffer, "changed", G_CALLBACK(callback_func), nullptr);
}

gboolean generic_hide_on_escape_controller_cb(GtkEventControllerKey  *controller,
                                              guint                  keyval,
                                              guint                  keycode,
                                              GdkModifierType        modifiers,
                                              GtkWidget              *to_be_hidden) {
   gboolean handled = TRUE;
   switch (keyval) {
      case GDK_KEY_Escape: {
         gtk_widget_set_visible(to_be_hidden, FALSE);
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

void setup_validation_graph_dialog() {

   GtkWidget *model_combobox = widget_from_builder("validation_graph_model_combobox");
   gtk_combo_box_set_model(GTK_COMBO_BOX(model_combobox),GTK_TREE_MODEL(graphics_info_t::validation_graph_model_list));
   gtk_combo_box_set_id_column(GTK_COMBO_BOX(model_combobox),0);

}

void setup_ramachandran_plot_chooser_dialog() {

   GtkWidget *model_combobox = widget_from_builder("ramachandran_plot_molecule_chooser_model_combobox");
   gtk_combo_box_set_model(GTK_COMBO_BOX(model_combobox), GTK_TREE_MODEL(graphics_info_t::ramachandran_plot_model_list));
   gtk_combo_box_set_id_column(GTK_COMBO_BOX(model_combobox),0);

}

void setup_get_monomer() {

   GtkWidget* frame = widget_from_builder("get_monomer_frame");
   GtkWidget* entry = widget_from_builder("get_monomer_entry");
   g_signal_connect(entry,"activate",G_CALLBACK(+[](GtkEntry* entry, gpointer user_data){
      handle_get_monomer_code(GTK_WIDGET(entry));
      GtkWidget *frame = widget_from_builder("get_monomer_frame");
      gtk_widget_set_visible(frame, FALSE);
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
   for (GtkWidget* child = gtk_widget_get_first_child(overlay);
       child != nullptr; 
       child = gtk_widget_get_next_sibling(child)) {
       if (child != to_skip)
          set_transparency_on_widget(child);
   }
}


void
add_python_scripting_entry_completion(GtkWidget *entry) {

   // call this *after* python has been setup!

   graphics_info_t g; // for history

   GtkEntryCompletion *completion = gtk_entry_completion_new();
   gtk_entry_completion_set_popup_completion(completion, TRUE);
   gtk_entry_completion_set_text_column(completion, 0);
   gtk_entry_completion_set_minimum_key_length(completion, 2);
   gtk_entry_set_completion(GTK_ENTRY(entry), completion);

   std::vector<std::string> completions;
   std::vector<std::string> module_coot_completions;
   std::vector<std::string> module_coot_utils_completions;

   PyErr_Clear();

   Py_ssize_t pos = 0;
   PyObject *key;
   PyObject *value;

   // note to future self... if you get a crash here then
   // it's because you've messed up the Python startup.
   // Perhaps by calling a function that no longer exists.

   // Get the module object for the `coot` module.
   PyObject *module = PyImport_ImportModule("coot");

   // this protects against crash on looking up the dictionary
   // when module is null
   if (module) {
      // Get the dictionary object for the `coot` module.
      PyObject *dict = PyModule_GetDict(module);
      // Iterate over the keys and values in the dictionary.
      while (PyDict_Next(dict, &pos, &key, &value)) {
	 // Do something interesting with the key and value.
	 // printf("Key: %s, Value: %s\n", PyUnicode_AsUTF8AndSize(key, NULL), PyUnicode_AsUTF8AndSize(value, NULL));
	 std::string key_c = std::string("coot.") +  (PyUnicode_AsUTF8AndSize(key, NULL));
	 module_coot_completions.push_back(key_c);
      }
      Py_DECREF(module);
      // Get the module object for the `sys` module.
      module = PyImport_ImportModule("coot_utils");

      if (module) {
	 if (false)
	    std::cout << "INFO:: add_python_scripting_entry_completion() coot_utils imported successfully" << std::endl;
      } else {
	 std::cout << "ERROR:: add_python_scripting_entry_completion() coot_utils import failure" << std::endl;
	 if (PyErr_Occurred())
	    PyErr_PrintEx(0);
	 return;
      }

      // Get the dictionary object for the `sys` module.
      dict = PyModule_GetDict(module);
      // Iterate over the keys and values in the dictionary.
      while (PyDict_Next(dict, &pos, &key, &value)) {
	 // Do something interesting with the key and value.
	 // printf("Key: %s, Value: %s\n", PyUnicode_AsUTF8AndSize(key, NULL), PyUnicode_AsUTF8AndSize(value, NULL));
	 std::string key_c = std::string("coot_utils.") +  (PyUnicode_AsUTF8AndSize(key, NULL));
	 module_coot_utils_completions.push_back(key_c);
      }
      Py_DECREF(module);

      // command history
      std::vector<std::string> chv = g.command_history.commands;

      chv = g.command_history.unique_commands(); // there *were* unique already

      if (false) chv.clear(); // 20230516-PE while testing.

      // add together the completions
      completions.push_back("import coot");
      completions.push_back("import coot_utils");
      completions.insert(completions.end(), chv.begin(),                           chv.end());
      completions.insert(completions.end(), module_coot_completions.begin(),       module_coot_completions.end());
      completions.insert(completions.end(), module_coot_utils_completions.begin(), module_coot_utils_completions.end());

      // maybe only once!
      GtkListStore *store = gtk_list_store_new(1, G_TYPE_STRING);
      GtkTreeIter iter;

      for (unsigned int i=0; i<completions.size(); i++) {
	 gtk_list_store_append( store, &iter );
	 std::string c = completions[i];
	 // std::cout << "adding to gtk-completion: " << c << std::endl;
	 gtk_list_store_set( store, &iter, 0, c.c_str(), -1 );
      }

      gtk_entry_completion_set_model(completion, GTK_TREE_MODEL(store));
   } else {
      std::cout << "ERROR:: in add_python_scripting_entry_completion() no module named 'coot'" << std::endl;
   }

}


gboolean
on_python_scripting_entry_key_pressed(GtkEventControllerKey *controller,
                                      guint                  keyval,
                                      guint                  keycode,
                                      GdkModifierType        modifiers,
                                      GtkEntry              *entry) {

   // This function is called on Ctrl and Shift, and Arrowkey Up and Down key presses

   gboolean handled = TRUE;
   bool control_is_pressed = (modifiers & GDK_CONTROL_MASK);

   std::cout << "on_python_scripting_entry_key_pressed() keyval: " << keyval << " keycode: " << keycode << std::endl;

   switch(keyval) {
      case GDK_KEY_Up: {
         handled = TRUE;
         if (control_is_pressed) {
            graphics_info_t g;
            std::string t = g.command_history.get_previous_command();
            gtk_editable_set_text(GTK_EDITABLE(entry), t.c_str());
         }
         break;
      }
      case GDK_KEY_Down: {
         handled = TRUE;
         if (control_is_pressed) {
            graphics_info_t g;
            std::string t = g.command_history.get_next_command();
            gtk_editable_set_text(GTK_EDITABLE(entry), t.c_str());
         }
         break;
      }
      case GDK_KEY_Escape: {
         auto func = +[] (gpointer data) {
            GtkRevealer* revealer = GTK_REVEALER(widget_from_builder("python_scripting_revealer"));
            gtk_revealer_set_reveal_child(revealer,FALSE);
            return gboolean(G_SOURCE_REMOVE);
         };
         g_idle_add(func, NULL);
         break;
      }
      default: {
         handled = FALSE;
         g_debug("Python scripting entry: Unhandled key: %s",gdk_keyval_name(keyval));
      }
   }

   return gboolean(handled);
}

void
on_python_scripting_entry_key_released(GtkEventControllerKey *controller,
                                       guint                  keyval,
                                       guint                  keycode,
                                       guint                  modifiers,
                                       GtkButton             *button) {

   graphics_info_t g;
   std::cout << "on_python_scripting_entry_key_released() keyval: " << keyval << " keycode: " << keycode << std::endl;

}


// 20230516-PE trying to add back the python completion and history that was in
// gtk3 coot into gtk4 coot.
//
// 20230516-PE I am, for the moment, not adding the header coot-setup-python.hh here because
// it doesn't include gtk stuff (for now).
void add_python_scripting_entry_completion(GtkWidget *entry);

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

   g_signal_connect(key_controller_entry, "key-pressed",  G_CALLBACK(on_python_scripting_entry_key_pressed),  entry);
   // g_signal_connect(key_controller_entry, "key-released", G_CALLBACK(on_python_scripting_entry_key_released), entry);

   gtk_widget_add_controller(entry, key_controller_entry);

   // for executing Python commands
   g_signal_connect(entry, "activate", G_CALLBACK(on_python_scripting_entry_activated), entry);

   // PE adds history and completions
   add_python_scripting_entry_completion(entry);
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
         GtkWidget* target = nullptr;
         if(GTK_IS_MENU_BUTTON(child)) {
#if GTK_MAJOR_VERSION == 4 && GTK_MINOR_VERSION >= 6
            target = gtk_menu_button_get_child(GTK_MENU_BUTTON(child));
#endif
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
#if GTK_MINOR_VERSION < 14 // 20241001-PE so that I don't see this when using fatal warnings
               g_warning("set_vertical_toolbar_internal_alignment: Toolbar item %p of type %s: "
               "The parent widget that wraps %s::child is not a GtkBox but a %s. "
               "%s::child however is a GtkBox. Attempt will be made to align it. It might not work.",
               child,G_OBJECT_TYPE_NAME(child),G_OBJECT_TYPE_NAME(child),G_OBJECT_TYPE_NAME(parent_widget),G_OBJECT_TYPE_NAME(child));
#endif
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

void setup_curlew_banner() {
    GtkWidget* curlew_banner = widget_from_builder("curlew_banner");
    std::string dir = coot::package_data_dir();
    std::string pixmaps_dir = coot::util::append_dir_dir(dir, "pixmaps");
    std::string banner_filepath = coot::util::append_dir_file(pixmaps_dir, "curlew-long.png");
    gtk_picture_set_filename(GTK_PICTURE(curlew_banner), banner_filepath.c_str());
}


// put this in a header
void tomo_scale_adjustment_changed(GtkAdjustment *adj, gpointer user_data);

void setup_tomo_widgets() {
   GtkWidget *scale = widget_from_builder("tomo_scale");
   GtkAdjustment *adjustment_current = gtk_range_get_adjustment(GTK_RANGE(scale));
   g_signal_connect(G_OBJECT(adjustment_current), "value_changed", G_CALLBACK(tomo_scale_adjustment_changed), NULL);
}


void setup_preferences() {

	// This stuff is no longer in preferences - perhaps they should be but they are not today,
	// so #ifdef out the code.

   // 20230627-PE put this in setup-gui-components - it should only happen once.
   // 20240916-PE done!
   {
#if 0
      // fill the bond combobox
      GtkComboBoxText *combobox = GTK_COMBO_BOX_TEXT(widget_from_preferences_builder("preferences_bond_width_combobox"));
      if (combobox) {
         for (int j = 1; j < 21; j++) {
            std::string s = graphics_info_t::int_to_string(j);
            gtk_combo_box_text_append_text(combobox, s.c_str());
         }
      } else {
         std::cout << "ERROR:: failed to find preferences_bond_width_combobox " << std::endl;
      }
#endif

#if 0
      // fill the font combobox
      combobox = GTK_COMBO_BOX_TEXT(widget_from_preferences_builder("preferences_font_size_combobox"));
      // 20230926-PE there was a crash here - maybe combobox was not looked up correctly.
      // Needs investigation, but add protection for now
      if (combobox) {
         std::vector<std::string> fonts;
         // fonts.push_back("Times Roman 10");
         // fonts.push_back("Times Roman 24");
         fonts.push_back("Fixed 8/13");
         fonts.push_back("Fixed 9/15");
         for (unsigned int j = 0; j < fonts.size(); j++)
            gtk_combo_box_text_append_text(combobox, fonts[j].c_str());
      } else {
         std::cout << "ERROR:: failed to find preferences_font_size_combobox" << std::endl;
      }
#endif
   }
}

void setup_aniso_hscale() {

   GtkWidget *scale = widget_from_builder("aniso_probability_hscale");
   GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(scale));
   if (adjustment) {
      gtk_adjustment_set_lower(adjustment, 0.0);
      gtk_adjustment_set_upper(adjustment, 1.0);
      gtk_adjustment_set_step_increment(adjustment, 1.0);
      gtk_adjustment_set_page_increment(adjustment, 5.0);
      gtk_adjustment_set_page_size(adjustment, 0.0);
      gtk_adjustment_set_value(adjustment, 0.5f);
      gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
      gtk_scale_set_digits(GTK_SCALE(scale), 2);
      gtk_scale_add_mark(GTK_SCALE(scale), 0.0, GTK_POS_BOTTOM, "0.0");
      gtk_scale_add_mark(GTK_SCALE(scale), 1.0, GTK_POS_BOTTOM, "1.0");
   }
}

void setup_gui_components() {

   g_info("Initializing UI components...");
   setup_menubuttons();
   setup_validation_graph_dialog();
   setup_mutate_residue_range_dialog();
   setup_ramachandran_plot_chooser_dialog();
   setup_get_monomer();
   setup_accession_code_frame();
   setup_python_scripting_entry();
   setup_curlew_banner();
   setup_tomo_widgets();
   setup_preferences();
   attach_css_style_class_to_overlays();
   set_vertical_toolbar_internal_alignment();
   setup_aniso_hscale();
   g_info("Done initializing UI components.");
}
