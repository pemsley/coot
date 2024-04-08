/*
 * validation-graphs/gtk4-test-sequence-view.cc
 *
 * Copyright 2023 by Medical Research Council
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

#include <string>
#include "coot-utils/atom-selection-container.hh"
#include "sequence-view-widget.hh"

void fill(GtkWidget *window, mmdb::Manager *mol) {

   GtkWidget *scrolled_window = gtk_scrolled_window_new();
   GtkWidget *overlay = gtk_overlay_new();
   gtk_overlay_set_child(GTK_OVERLAY(overlay), GTK_WIDGET(scrolled_window));
   gtk_widget_set_hexpand(scrolled_window, TRUE);
   gtk_widget_set_vexpand(scrolled_window, TRUE);
   gtk_widget_set_size_request(scrolled_window, 200, 130);
   gtk_window_set_child(GTK_WINDOW(window), overlay);
   CootSequenceView *sv = coot_sequence_view_new();
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(sv));
   gtk_widget_set_margin_start(GTK_WIDGET(sv), 10);
   gtk_widget_set_margin_end(GTK_WIDGET(sv), 10);
   gtk_widget_set_margin_top(GTK_WIDGET(sv), 10);
   gtk_widget_set_margin_bottom(GTK_WIDGET(sv), 10);
   int imol = 0; // for now
   coot_sequence_view_set_structure(sv, imol, mol);

   auto click_function_callback = +[] (CootSequenceView* self, int imol, const coot::residue_spec_t &spec, gpointer user_data) {
      std::cout << " clicked: " << imol << " " << spec << std::endl;
   };
   g_signal_connect(sv, "residue-clicked", G_CALLBACK(click_function_callback), nullptr);

   GtkWidget *button = gtk_button_new_with_label("Close");
   
   gtk_overlay_add_overlay(GTK_OVERLAY(overlay), button);
   gtk_widget_set_halign(GTK_WIDGET(button),GTK_ALIGN_START);
   gtk_widget_set_valign(GTK_WIDGET(button),GTK_ALIGN_END);

}

int main(int argc, char **argv) {

   int status = 0;
   if (argc > 1) {
      std::string pdb_file_name = argv[1];
      auto atom_sel = get_atom_selection(pdb_file_name, false, true, false);
      if (atom_sel.read_success) {
         gtk_init();
         g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);
         GtkApplication* app = gtk_application_new("org.pemsley.test-sequence-view", G_APPLICATION_DEFAULT_FLAGS);
         GError *error = NULL;
         g_application_register(G_APPLICATION(app), NULL, &error);
         if (error) {
            std::cout << "ERROR" << error->message << std::endl;
         } else {
            // let's go
            auto activate_callback = +[] (GtkApplication *app, gpointer user_data) {
               mmdb::Manager *mol = static_cast<mmdb::Manager *>(user_data);
               GtkWidget *win = gtk_application_window_new(app);
               gtk_application_add_window(app, GTK_WINDOW(win));
               gtk_window_set_application(GTK_WINDOW(win), app);
               gtk_window_set_title(GTK_WINDOW(win), "Coot: Sequence View");
               int width  = 500;
               int height = 50;
               gtk_window_set_default_size(GTK_WINDOW(win), width, height);
               fill(win, mol);
               gtk_widget_set_visible(win, TRUE);
            };
            g_signal_connect(app, "activate", G_CALLBACK(activate_callback), atom_sel.mol);

            // Ideally the arguments should be argc, argv and
            // call above g_application_new with G_APPLICATION_HANDLES_OPEN
            // and connect a signal to "open":
            // g_signal_connect(app, "open", G_CALLBACK(open_callback), activate_data).
            return g_application_run(G_APPLICATION(app), 0, 0);
         }
      
      } else {
         status = 1; // bad input
      }
   } else {
      status = 1; // bad input
   }

   return status;
}
