/*
 * src/command-palette.cc
 *
 * Copyright 2026 by Medical Research Council
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

#ifdef USE_PYTHON
#include <Python.h>
#endif

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

#include <gtk/gtk.h>

#include "graphics-info.h"
#include "widget-from-builder.hh"
#include "command-palette.hh"

#include "utils/logging.hh"
extern logging logger;

// The idea of searching the menus with a GtkSearchEntry is from
// Oliver Clarke's coot1_trimmings.

namespace {

   class command_palette_item_t {
   public:
      std::string label;
      std::string path;       // e.g. "Calculate > Modelling" or "Key Binding"
      std::string search_text; // lower-case path and label
      std::string doc_text;    // lower-case name and description of the "api-function" of the menu item
      std::string action;     // without the "app." prefix. Empty for a key binding.
      GVariant *target;       // can be null
      keyboard_key_t key;     // used if action is empty
      command_palette_item_t(const std::string &label_in, const std::string &path_in) :
         label(label_in), path(path_in), target(nullptr), key(0) {
         gchar *lc = g_utf8_strdown((path + " " + label).c_str(), -1);
         search_text = lc;
         g_free(lc);
      }
   };

   // refilled every time the palette is shown
   std::vector<command_palette_item_t> command_palette_items;

   // GMenu labels use "_" for the mnemonic and "__" for a real underscore
   std::string strip_mnemonic(const std::string &s) {
      std::string r;
      for (unsigned int i=0; i<s.length(); i++) {
         if (s[i] == '_') {
            if (i+1 < s.length() && s[i+1] == '_') {
               r += '_';
               i++;
            }
         } else {
            r += s[i];
         }
      }
      return r;
   }

   // Menu items can have an "api-function" attribute: the name of the scripting function that does what the
   // menu item does. It is a function in the coot python module unless it is given as "module.function".
   // Return the lower-case name and description (the docstring, without the signature and the parameters)
   // of that function, so that we can search on those too.
   std::string api_function_search_text(const std::string &api_function) {

      static std::map<std::string, std::string> cache;
      std::map<std::string, std::string>::const_iterator it = cache.find(api_function);
      if (it != cache.end()) return it->second;

      std::string text = api_function;
#ifdef USE_PYTHON
      std::string module_name = "coot";
      std::string function_name = api_function;
      std::string::size_type dot_pos = api_function.rfind('.');
      if (dot_pos != std::string::npos) {
         module_name   = api_function.substr(0, dot_pos);
         function_name = api_function.substr(dot_pos+1);
      }
      // Don't import the module here (that can have side-effects), use it only if it is already loaded.
      PyObject *module_name_py = PyUnicode_FromString(module_name.c_str());
      PyObject *module = PyImport_GetModule(module_name_py);
      if (module) {
         PyObject *function = PyObject_GetAttrString(module, function_name.c_str());
         if (function) {
            PyObject *doc = PyObject_GetAttrString(function, "__doc__");
            const char *doc_c = (doc && PyUnicode_Check(doc)) ? PyUnicode_AsUTF8(doc) : nullptr;
            if (doc_c) {
               std::istringstream doc_stream(doc_c);
               std::string line;
               while (std::getline(doc_stream, line)) {
                  std::string::size_type first_char_pos = line.find_first_not_of(" \t");
                  if (first_char_pos == std::string::npos) continue;
                  line = line.substr(first_char_pos);
                  if (line == "Parameters") break;                              // the rest is the arguments
                  if (line.compare(0, function_name.length() + 1, function_name + "(") == 0) continue; // signature
                  text += " ";
                  text += line;
               }
            }
            Py_XDECREF(doc);
            Py_DECREF(function);
         } else {
            logger.log(log_t::WARNING, "command palette: api-function not found:", api_function);
         }
         Py_DECREF(module);
      }
      Py_XDECREF(module_name_py);
      PyErr_Clear();
#endif
      gchar *lc = g_utf8_strdown(text.c_str(), -1);
      text = lc;
      g_free(lc);
      cache[api_function] = text;
      return text;
   }

   void add_menu_items(GMenuModel *menu, const std::string &path, std::vector<command_palette_item_t> *items_p) {

      GActionGroup *action_group = G_ACTION_GROUP(graphics_info_t::application);
      int n = g_menu_model_get_n_items(menu);
      for (int i=0; i<n; i++) {
         std::string label;
         gchar *l = nullptr;
         if (g_menu_model_get_item_attribute(menu, i, G_MENU_ATTRIBUTE_LABEL, "s", &l)) {
            label = strip_mnemonic(l);
            g_free(l);
         }
         GMenuModel *sub_menu = g_menu_model_get_item_link(menu, i, G_MENU_LINK_SUBMENU);
         if (! sub_menu)
            sub_menu = g_menu_model_get_item_link(menu, i, G_MENU_LINK_SECTION);
         if (sub_menu) {
            std::string sub_path = path;
            if (! label.empty()) {
               if (! sub_path.empty()) sub_path += " > ";
               sub_path += label;
            }
            add_menu_items(sub_menu, sub_path, items_p);
            g_object_unref(sub_menu);
         } else {
            gchar *a = nullptr;
            if (g_menu_model_get_item_attribute(menu, i, G_MENU_ATTRIBUTE_ACTION, "s", &a)) {
               std::string action(a);
               g_free(a);
               if (! label.empty() && action.compare(0, 4, "app.") == 0) {
                  action = action.substr(4);
                  if (g_action_group_has_action(action_group, action.c_str())) {
                     command_palette_item_t item(label, path);
                     item.action = action;
                     item.target = g_menu_model_get_item_attribute_value(menu, i, G_MENU_ATTRIBUTE_TARGET, nullptr);
                     gchar *af = nullptr;
                     if (g_menu_model_get_item_attribute(menu, i, "api-function", "s", &af)) {
                        item.doc_text = api_function_search_text(af);
                        g_free(af);
                     }
                     items_p->push_back(item);
                  }
               }
            }
         }
      }
   }

   void fill_command_palette_items() {

      for (auto &item : command_palette_items)
         if (item.target)
            g_variant_unref(item.target);
      command_palette_items.clear();

      // the main menubar - which includes the menu items added by scripting
      GMenuModel *menubar = gtk_application_get_menubar(graphics_info_t::application);
      if (menubar)
         add_menu_items(menubar, "", &command_palette_items);

      // the toolbar and module menus
      std::vector<std::pair<std::string, std::string> > other_menus = {
         {"refine-menu",           "Refine"},
         {"fixed-atoms-menu",      "Fixed Atoms"},
         {"delete-item-menu",      "Delete"},
         {"rotate-translate-menu", "Rotate/Translate"},
         {"rigid-body-fit-menu",   "Rigid-Body Fit"},
         {"ccp4-menu",             "CCP4"},
         {"cryo-em-menu",          "Cryo-EM"},
         {"ligand-menu",           "Ligand"},
         {"restraints-menu",       "Restraints"}};
      for (const auto &m : other_menus) {
         GObject *o = graphics_info_t::get_gobject_from_builder(m.first);
         if (o)
            add_menu_items(G_MENU_MODEL(o), m.second, &command_palette_items);
      }

      for (const auto &kb : graphics_info_t::key_bindings_map) {
         if (! kb.second.description.empty()) {
            command_palette_item_t item(kb.second.description, "Key Binding");
            item.key = kb.first;
            command_palette_items.push_back(item);
         }
      }
   }

   bool item_is_enabled(const command_palette_item_t &item) {
      if (item.action.empty()) return true;
      return g_action_group_get_action_enabled(G_ACTION_GROUP(graphics_info_t::application), item.action.c_str());
   }

   std::string key_binding_key_label(const keyboard_key_t &key) {
      std::string s;
      if (key.ctrl_is_pressed) s = "Ctrl ";
      const gchar *n = gdk_keyval_name(key.gdk_key);
      if (n) s += n;
      return s;
   }

   void hide_command_palette() {
      gtk_widget_set_visible(widget_from_builder("command_palette_frame"), FALSE);
      graphics_info_t::graphics_grab_focus();
   }

   void run_row(GtkListBoxRow *row) {

      if (! row) return;
      unsigned int idx = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(row), "item_index"));
      if (idx >= command_palette_items.size()) return;
      const command_palette_item_t &item = command_palette_items[idx];
      if (! item_is_enabled(item)) return;
      hide_command_palette();
      if (item.action.empty()) {
         std::map<keyboard_key_t, key_bindings_t>::const_iterator it = graphics_info_t::key_bindings_map.find(item.key);
         if (it != graphics_info_t::key_bindings_map.end())
            it->second.run();
         graphics_info_t::graphics_draw();
      } else {
         g_action_group_activate_action(G_ACTION_GROUP(graphics_info_t::application), item.action.c_str(), item.target);
      }
   }

   GtkWidget *make_row(unsigned int item_index) {

      const command_palette_item_t &item = command_palette_items[item_index];
      GtkWidget *row = gtk_list_box_row_new();
      g_object_set_data(G_OBJECT(row), "item_index", GUINT_TO_POINTER(item_index));
      GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 12);
      GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 1);
      GtkWidget *label = gtk_label_new(item.label.c_str());
      GtkWidget *path_label = gtk_label_new(item.path.c_str());
      gtk_label_set_xalign(GTK_LABEL(label), 0.0);
      gtk_label_set_xalign(GTK_LABEL(path_label), 0.0);
      gtk_label_set_ellipsize(GTK_LABEL(path_label), PANGO_ELLIPSIZE_START);
      gtk_widget_add_css_class(path_label, "dim-label");
      gtk_widget_add_css_class(path_label, "caption");
      gtk_box_append(GTK_BOX(vbox), label);
      gtk_box_append(GTK_BOX(vbox), path_label);
      gtk_widget_set_hexpand(vbox, TRUE);
      gtk_box_append(GTK_BOX(hbox), vbox);
      if (item.action.empty()) {
         GtkWidget *key_label = gtk_label_new(key_binding_key_label(item.key).c_str());
         gtk_widget_add_css_class(key_label, "dim-label");
         gtk_box_append(GTK_BOX(hbox), key_label);
      }
      gtk_widget_set_margin_top(hbox, 3);
      gtk_widget_set_margin_bottom(hbox, 3);
      gtk_widget_set_margin_start(hbox, 6);
      gtk_widget_set_margin_end(hbox, 6);
      // grey out, rather than make the row itself insensitive, so that the selection can still move over it
      gtk_widget_set_sensitive(hbox, item_is_enabled(item));
      gtk_list_box_row_set_child(GTK_LIST_BOX_ROW(row), hbox);
      return row;
   }

   void on_search_changed(GtkSearchEntry *entry, gpointer user_data) {

      GtkWidget *list_box = widget_from_builder("command_palette_listbox");
      GtkWidget *scrolled_window = widget_from_builder("command_palette_scrolled_window");
      GtkWidget *child = gtk_widget_get_first_child(list_box);
      while (child) {
         GtkWidget *next = gtk_widget_get_next_sibling(child);
         gtk_list_box_remove(GTK_LIST_BOX(list_box), child);
         child = next;
      }

      gchar *lc = g_utf8_strdown(gtk_editable_get_text(GTK_EDITABLE(entry)), -1);
      std::string query(lc);
      g_free(lc);
      std::vector<std::string> terms = coot::util::split_string_no_blanks(query, " ");
      std::vector<unsigned int> matches;     // all the terms are in the path and label
      std::vector<unsigned int> doc_matches; // the api-function name or description was needed for the match
      if (! terms.empty()) {
         for (unsigned int i=0; i<command_palette_items.size(); i++) {
            const command_palette_item_t &item = command_palette_items[i];
            bool all_found = true;
            bool needed_doc = false;
            for (const auto &term : terms) {
               if (item.search_text.find(term) == std::string::npos) {
                  if (item.doc_text.find(term) == std::string::npos) {
                     all_found = false;
                     break;
                  }
                  needed_doc = true;
               }
            }
            if (all_found) {
               if (needed_doc)
                  doc_matches.push_back(i);
               else
                  matches.push_back(i);
            }
         }
      }

      // items for which the label starts with the first search term go first
      if (! terms.empty()) {
         const std::string &t = terms[0];
         auto label_starts_with_term = [&t] (unsigned int idx) {
            const command_palette_item_t &item = command_palette_items[idx];
            // the label is the end of the search text
            std::string::size_type label_start = item.search_text.length() - item.label.length();
            return item.search_text.compare(label_start, t.length(), t) == 0;
         };
         std::stable_partition(matches.begin(), matches.end(), label_starts_with_term);
      }
      // the matches on the label and path go before those that needed the function documentation
      matches.insert(matches.end(), doc_matches.begin(), doc_matches.end());

      for (unsigned int idx : matches)
         gtk_list_box_append(GTK_LIST_BOX(list_box), make_row(idx));

      gtk_widget_set_visible(scrolled_window, ! matches.empty());
      GtkListBoxRow *first_row = gtk_list_box_get_row_at_index(GTK_LIST_BOX(list_box), 0);
      if (first_row)
         gtk_list_box_select_row(GTK_LIST_BOX(list_box), first_row);
      GtkAdjustment *adj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrolled_window));
      gtk_adjustment_set_value(adj, 0.0);
   }

   void move_selection(int delta) {

      GtkWidget *list_box = widget_from_builder("command_palette_listbox");
      GtkWidget *scrolled_window = widget_from_builder("command_palette_scrolled_window");
      GtkListBoxRow *row = gtk_list_box_get_selected_row(GTK_LIST_BOX(list_box));
      if (! row) return;
      GtkListBoxRow *new_row = gtk_list_box_get_row_at_index(GTK_LIST_BOX(list_box),
                                                             gtk_list_box_row_get_index(row) + delta);
      if (! new_row) return;
      gtk_list_box_select_row(GTK_LIST_BOX(list_box), new_row);
      // keep the selected row in view (the focus stays in the search entry)
      graphene_rect_t bounds;
      if (gtk_widget_compute_bounds(GTK_WIDGET(new_row), list_box, &bounds)) {
         GtkAdjustment *adj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(scrolled_window));
         gtk_adjustment_clamp_page(adj, bounds.origin.y, bounds.origin.y + bounds.size.height);
      }
   }

   gboolean on_search_entry_key_pressed(GtkEventControllerKey *controller, guint keyval, guint keycode,
                                        GdkModifierType modifiers, gpointer user_data) {
      gboolean handled = FALSE;
      if (keyval == GDK_KEY_Down) {
         move_selection(1);
         handled = TRUE;
      }
      if (keyval == GDK_KEY_Up) {
         move_selection(-1);
         handled = TRUE;
      }
      if (keyval == GDK_KEY_Escape) {
         hide_command_palette();
         handled = TRUE;
      }
      return handled;
   }

   void setup_command_palette() {

      GtkWidget *frame    = widget_from_builder("command_palette_frame");
      GtkWidget *entry    = widget_from_builder("command_palette_search_entry");
      GtkWidget *list_box = widget_from_builder("command_palette_listbox");

#if GTK_CHECK_VERSION(4,8,0)
      gtk_search_entry_set_search_delay(GTK_SEARCH_ENTRY(entry), 100); // ms
#endif

      g_signal_connect(entry, "search-changed", G_CALLBACK(on_search_changed), nullptr);

      auto on_entry_activate = +[] (GtkSearchEntry *entry, gpointer user_data) {
         GtkWidget *list_box = widget_from_builder("command_palette_listbox");
         run_row(gtk_list_box_get_selected_row(GTK_LIST_BOX(list_box)));
      };
      g_signal_connect(entry, "activate", G_CALLBACK(on_entry_activate), nullptr);

      auto on_row_activated = +[] (GtkListBox *list_box, GtkListBoxRow *row, gpointer user_data) {
         run_row(row);
      };
      g_signal_connect(list_box, "row-activated", G_CALLBACK(on_row_activated), nullptr);

      // capture phase, so that we see the arrow keys and Escape before the text widget does
      GtkEventController *key_controller = gtk_event_controller_key_new();
      gtk_event_controller_set_propagation_phase(key_controller, GTK_PHASE_CAPTURE);
      g_signal_connect(key_controller, "key-pressed", G_CALLBACK(on_search_entry_key_pressed), nullptr);
      gtk_widget_add_controller(entry, key_controller);

      // hide the palette when the focus goes elsewhere (e.g. a click in the graphics). The controller
      // is on the frame, not the entry, so that clicking on a row of the list doesn't count.
      // Don't grab the focus for the graphics here, the focus is wherever the user put it.
      auto on_focus_leave = +[] (GtkEventControllerFocus *controller, gpointer user_data) {
         gtk_widget_set_visible(GTK_WIDGET(user_data), FALSE);
      };
      GtkEventController *focus_controller = gtk_event_controller_focus_new();
      g_signal_connect(focus_controller, "leave", G_CALLBACK(on_focus_leave), frame);
      gtk_widget_add_controller(frame, focus_controller);
   }
}

void show_command_palette() {

   static bool is_setup = false;
   if (! is_setup) {
      setup_command_palette();
      is_setup = true;
   }

   GtkWidget *frame = widget_from_builder("command_palette_frame");
   GtkWidget *entry = widget_from_builder("command_palette_search_entry");
   fill_command_palette_items();
   gtk_editable_set_text(GTK_EDITABLE(entry), "");
   on_search_changed(GTK_SEARCH_ENTRY(entry), nullptr); // clear the results of the previous search
   gtk_widget_set_visible(frame, TRUE);
   gtk_widget_grab_focus(entry);
}
