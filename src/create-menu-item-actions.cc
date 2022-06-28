
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

   // void gtk_file_chooser_add_choice (GtkFileChooser* chooser,
   //                                   const char* id,
   //                                   const char* label,
   //                                   const char** options,
   //                                   const char** option_labels)


   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "recentre", "Centre on New Molecule", NULL, NULL);
   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "recentre", "No Recentre", NULL, NULL);
   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "recentre", "Move Molecule Here", NULL, NULL);

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

void open_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *dataset_chooser = widget_from_builder("map_name_filechooser_dialog");
   set_directory_for_filechooser(dataset_chooser);
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

void get_monomer_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                   G_GNUC_UNUSED GVariant *parameter,
                   G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *dialog = wrapped_create_libcheck_monomer_dialog();
  gtk_widget_show(dialog);
}

void import_cif_dictionary_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                  G_GNUC_UNUSED GVariant *parameter,
                                  G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *chooser = widget_from_builder("cif_dictionary_filechooser_dialog");
   set_directory_for_filechooser(chooser);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, chooser);
   gtk_widget_show(chooser);
}

void toggle_display_frames_per_second_action (G_GNUC_UNUSED GSimpleAction *simple_action,
                                              G_GNUC_UNUSED GVariant *parameter,
                                             G_GNUC_UNUSED gpointer user_data) {

int state = get_fps_flag();
if (state == 1)
   set_show_fps(0);
else
   set_show_fps(1);
}

void
search_monomer_library_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = widget_from_builder("monomer_search_dialog");
   gtk_widget_show(w);
}

void
fetch_pdb_using_code_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {

   int n = COOT_ACCESSION_CODE_WINDOW_OCA;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}

void
fetch_pdb_and_map_using_eds_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {

   int n = COOT_ACCESSION_CODE_WINDOW_EDS;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}

void
fetch_pdb_and_map_using_pdb_redo_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {

   int n = COOT_ACCESSION_CODE_WINDOW_PDB_REDO;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}

void
save_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {
      GCallback callback_func = G_CALLBACK(save_molecule_coords_combobox_changed);
   int imol = first_coords_imol();
   int imol_unsaved = first_unsaved_coords_imol();
   if (imol_unsaved != -1)
      imol = imol_unsaved;
   std::cout << "DEBUG:: in on_save_coordinates1_activate() with imol_unsaved "
             << imol_unsaved << std::endl;
   set_save_molecule_number(imol); /* set *save* molecule number */

   // this is the molecule chooser, not the file chooser
   //
   GtkWidget *widget = widget_from_builder("save_coords_dialog");
   GtkWidget *combobox = widget_from_builder("save_coordinates_combobox");

   if (combobox) {
      fill_combobox_with_coordinates_options(combobox, callback_func, imol);
      set_transient_and_position(COOT_UNDEFINED_WINDOW, widget);
      gtk_widget_show(widget);
      gtk_window_present(GTK_WINDOW(widget));
   } else {
      std::cout << "ERROR:: in on_save_coordinates1_activate() bad combobox!\n";
   }
}

void
save_symmetry_coordinates_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {

   setup_save_symmetry_coords();

}

void
save_state_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *file_chooser = coot_save_state_chooser();
   gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(file_chooser), save_state_file_name_raw());
   add_filename_filter_button(file_chooser, COOT_SCRIPTS_FILE_SELECTION);
   set_file_selection_dialog_size(file_chooser);
   gtk_widget_show(file_chooser);
}

void
recover_session_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   recover_session();
}

void
file_export_map_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   short int is_fragment = false;
   export_map_gui(is_fragment);

}

void
file_export_map_fragment_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                G_GNUC_UNUSED GVariant *parameter,
                                G_GNUC_UNUSED gpointer user_data) {
   short int is_fragment = true;
   export_map_gui(is_fragment);

}

void
close_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *widget = wrapped_create_new_close_molecules_dialog(); // uses builder
   gtk_widget_show(widget);
}

void
change_chain_ids_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_change_chain_id_dialog(); // uses builder
   gtk_widget_show(w);
}

void
copy_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
   do_edit_copy_molecule();
}

void
copy_molecule_fragment_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {
   do_edit_copy_fragment();
}

void
merge_molecules_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   gtk_widget_show(w);
}

void
move_molecule_here_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = widget_from_builder("move_molecule_here_dialog");
   fill_move_molecule_here_dialog(w);
   gtk_widget_show(w);
}

void
mutate_molecule_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_mutate_sequence_dialog();
   gtk_widget_show(w);
}

void
edit_replace_fragment_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {
     do_edit_replace_fragment();
}

void
edit_replace_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {
     do_edit_replace_residue();
}

void
renumber_residues_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_renumber_residue_range_dialog();
   gtk_widget_show(w);
}

void
residue_info_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {
   do_residue_info_dialog();
}

void
edit_restraints_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w =  wrapped_create_residue_editor_select_monomer_type_dialog();
   gtk_widget_show(w);
}

#include "c-interface-preferences.h"
void
show_preferences_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   show_preferences();
   update_preference_gui();
}

void fill_and_show_shader_preferences(); // should this be in a header?

void
show_shader_preferences_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {
   fill_and_show_shader_preferences();
}

void
create_actions(GtkApplication *application) {

   auto add_action = [application] (const std::string &action_name,
                                    void (*action_function) (GSimpleAction *simple_action,
                                                             GVariant *parameter,
                                                             gpointer user_data)) {
      GSimpleAction *simple_action = g_simple_action_new(action_name.c_str(), NULL);
      g_action_map_add_action(G_ACTION_MAP(application), G_ACTION(simple_action));
      g_signal_connect(simple_action, "activate", G_CALLBACK(action_function), NULL);
   };

   add_action("open_coordinates_action", open_coordinates_action);
   add_action(   "auto_open_mtz_action",    auto_open_mtz_action);
   add_action(    "open_dataset_action",     open_dataset_action);
   add_action(        "open_map_action",         open_map_action);

   add_action("get_monomer_action", get_monomer_action);
   add_action(     "curlew_action",      curlew_action);
   add_action(       "exit_action",        exit_action);

   add_action(            "import_cif_dictionary_action",            import_cif_dictionary_action);
   add_action( "toggle_display_frames_per_second_action", toggle_display_frames_per_second_action);
   add_action(           "search_monomer_library_action",           search_monomer_library_action);
   add_action(             "fetch_pdb_using_code_action",             fetch_pdb_using_code_action);
   add_action(      "fetch_pdb_and_map_using_eds_action",      fetch_pdb_and_map_using_eds_action);
   add_action( "fetch_pdb_and_map_using_pdb_redo_action", fetch_pdb_and_map_using_pdb_redo_action);
   add_action(                 "save_coordinates_action",                 save_coordinates_action);
   add_action(        "save_symmetry_coordinates_action",        save_symmetry_coordinates_action);
   add_action(                       "save_state_action",                       save_state_action);
   add_action(                  "recover_session_action",                  recover_session_action);
   add_action(                  "file_export_map_action",                  file_export_map_action);
   add_action(         "file_export_map_fragment_action",                  file_export_map_action);
   add_action(                   "close_molecule_action",                   close_molecule_action);

   add_action(       "change_chain_ids_action",        change_chain_ids_action);
   add_action(          "copy_molecule_action",           copy_molecule_action);
   add_action( "copy_molecule_fragment_action",  copy_molecule_fragment_action);
   add_action(        "merge_molecules_action",         merge_molecules_action);
   add_action(     "move_molecule_here_action",      move_molecule_here_action);
   add_action(        "mutate_molecule_action",         mutate_molecule_action);
   add_action(  "edit_replace_fragment_action",   edit_replace_fragment_action);
   add_action(   "edit_replace_residue_action",    edit_replace_residue_action);
   add_action(      "renumber_residues_action",       renumber_residues_action);
   add_action(           "residue_info_action",            residue_info_action);
   add_action(        "edit_restraints_action",         edit_restraints_action);
   add_action(       "show_preferences_action",        show_preferences_action);
   add_action("show_shader_preferences_action", show_shader_preferences_action);

   add_action("load_tutorial_model_and_data_action", load_tutorial_model_and_data_action);
}

