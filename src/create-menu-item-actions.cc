
#include <iostream>
#include <gtk/gtk.h>

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"
#include "widget-from-builder.hh"

extern "C" { void load_tutorial_model_and_data(); }

void on_coords_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file   = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);

#if 0
      GSList *files_list = gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(dialog));
      while (files_list) {

         const char *fnc = static_cast<const char *>(files_list->data);
         if (fnc) {
            std::string fn(fnc);
            handle_read_draw_molecule_with_recentre(fn, 0);
         }
         files_list = g_slist_next(files_list);
      }
#endif

      GtkWidget *recentre_combobox = widget_from_builder("coords_filechooserdialog_recentre_combobox");
      int active_item_index = gtk_combo_box_get_active(GTK_COMBO_BOX(recentre_combobox));
      bool move_molecule_here_flag = false;
      bool recentre_on_read_pdb_flag = false;
      if (active_item_index == 0)
         recentre_on_read_pdb_flag = true;
      if (active_item_index == 1)
         recentre_on_read_pdb_flag = false;
      if (active_item_index == 2)
         move_molecule_here_flag = true;

      // open_file (file);
      if (file_name) {
         std::cout << "info: " << file_name << " " << move_molecule_here_flag << " " << recentre_on_read_pdb_flag
                   << std::endl;

         if (move_molecule_here_flag) {
            handle_read_draw_molecule_and_move_molecule_here(file_name);
         } else {
            if (recentre_on_read_pdb_flag)
               handle_read_draw_molecule_with_recentre(file_name, 1);
            else
               handle_read_draw_molecule_with_recentre(file_name, 0); // no recentre
         }
      }
   }
   gtk_window_close(GTK_WINDOW(dialog));
}


void open_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   // Ancient GTK3
   // open_coords_dialog();

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File",
                                                   parent_window,
                                                   action,
                                                   _("_Cancel"),
                                                   GTK_RESPONSE_CANCEL,
                                                   _("_Open"),
                                                   GTK_RESPONSE_ACCEPT,
                                                   NULL);

   g_signal_connect(dialog, "response", G_CALLBACK(on_coords_filechooser_dialog_response_gtk4), NULL);
   gtk_widget_show(dialog);

}

void open_dataset_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *dataset_chooser = widget_from_builder("coords_filechooser_dialog");
   GtkWidget *main_window = graphics_info_t::get_main_window();
   gtk_window_set_transient_for(GTK_WINDOW(dataset_chooser), GTK_WINDOW(main_window));

   set_directory_for_filechooser(dataset_chooser);
   set_file_selection_dialog_size(dataset_chooser);
   add_filechooser_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);

}

void auto_open_mtz_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {


   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   int is_auto_read_fileselection = 1;
   set_directory_for_filechooser(dataset_chooser);
   add_filename_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   g_object_set_data(G_OBJECT(dataset_chooser), "imol", GINT_TO_POINTER(-1)); // 20220627-PE do I need this?
   g_object_set_data(G_OBJECT(dataset_chooser), "is_auto", GINT_TO_POINTER(is_auto_read_fileselection));
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);

}

void load_tutorial_model_and_data_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                         G_GNUC_UNUSED GVariant *parameter,
                                         G_GNUC_UNUSED gpointer user_data) {
   load_tutorial_model_and_data();
}

void exit_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   coot_checked_exit(0);
}

#include "curlew.h" // 20220628-PE why does this exist? why is curlew() not in the .hh file?

void curlew_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                   G_GNUC_UNUSED GVariant *parameter,
                   G_GNUC_UNUSED gpointer user_data) {
   curlew();
}

void
create_actions(GtkApplication *application) {

   GSimpleAction *simple_action;

   GtkWindow *application_window = gtk_application_get_active_window(application);

   simple_action = g_simple_action_new("open_coordinates_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(open_coordinates_action), application_window);

   simple_action = g_simple_action_new("auto_open_mtz_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(auto_open_mtz_action), application_window);

   simple_action = g_simple_action_new("open_dataset_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(open_dataset_action), application_window);

   simple_action = g_simple_action_new("load_tutorial_model_and_data_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(load_tutorial_model_and_data_action), NULL);

   simple_action = g_simple_action_new("curlew_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(curlew_action), NULL);

   simple_action = g_simple_action_new("exit_action", NULL);
   g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
   g_signal_connect(simple_action, "activate", G_CALLBACK(exit_action), NULL);
}

