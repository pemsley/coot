
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

void on_dataset_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                 int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);
      bool auto_read_flag = false, ismtz = false, ismtzauto = false, iscnsauto = false;
      auto_read_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(chooser), "auto_read_flag"));
      ismtz = is_mtz_file_p(file_name);

      if (ismtz)
         ismtzauto = mtz_file_has_phases_p(file_name);
      else
         iscnsauto = cns_file_has_phases_p(file_name);

      if ( ismtzauto || iscnsauto ) {
         if (auto_read_flag)
            wrapped_auto_read_make_and_draw_maps(file_name);
         else
            manage_column_selector(file_name); /* calls create_column_label_window(), fills and displays.*/
      } else {
         manage_column_selector(file_name); // try read a cif (strangely)
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


   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "recentre",      "Centre on New Molecule", NULL, NULL);
   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "no-recentre",   "No Recentre",            NULL, NULL);
   gtk_file_chooser_add_choice(GTK_FILE_CHOOSER(dialog), "move-mol-here", "Move Molecule Here",     NULL, NULL);

   g_signal_connect(dialog, "response", G_CALLBACK(on_coords_filechooser_dialog_response_gtk4), NULL);
   gtk_widget_show(dialog);

}

void open_dataset_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {
#if 0
   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   GtkWidget *main_window = graphics_info_t::get_main_window();
   gtk_window_set_transient_for(GTK_WINDOW(dataset_chooser), GTK_WINDOW(main_window));

   set_directory_for_filechooser(dataset_chooser);
   set_file_selection_dialog_size(dataset_chooser);
   add_filechooser_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);
#endif

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", parent_window, action,
                                                   _("_Cancel"), GTK_RESPONSE_CANCEL,
                                                   _("_Open"), GTK_RESPONSE_ACCEPT,
                                                   NULL);
   g_signal_connect(dialog, "response", G_CALLBACK(on_dataset_filechooser_dialog_response_gtk4), NULL);
   g_object_set_data(G_OBJECT(dialog), "auto_read_flag", GINT_TO_POINTER(FALSE));
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.mtz");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_show(dialog);
}

void auto_open_mtz_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {

#if 0
   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   int is_auto_read_fileselection = 1;
   set_directory_for_filechooser(dataset_chooser);
   add_filename_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   g_object_set_data(G_OBJECT(dataset_chooser), "imol", GINT_TO_POINTER(-1)); // 20220627-PE do I need this?
   g_object_set_data(G_OBJECT(dataset_chooser), "is_auto", GINT_TO_POINTER(is_auto_read_fileselection));
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);
#endif

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", parent_window, action,
                                                   _("_Cancel"), GTK_RESPONSE_CANCEL,
                                                   _("_Open"), GTK_RESPONSE_ACCEPT,
                                                   NULL);
   g_signal_connect(dialog, "response", G_CALLBACK(on_dataset_filechooser_dialog_response_gtk4), NULL);
   g_object_set_data(G_OBJECT(dialog), "auto_read_flag", GINT_TO_POINTER(TRUE));
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.mtz");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_show(dialog);
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


void on_cif_dictionary_filechooser_dialog_response_gtk4(GtkDialog *dialog,
                                                 int        response) {

   if (response == GTK_RESPONSE_ACCEPT) {
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);
      int imol_enc = -999997;
      short int new_molecule_checkbutton_state = 0;
      handle_cif_dictionary_for_molecule(file_name, imol_enc, new_molecule_checkbutton_state);
   }
   gtk_widget_hide(GTK_WIDGET(dialog));
}

void import_cif_dictionary_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                  G_GNUC_UNUSED GVariant *parameter,
                                  G_GNUC_UNUSED gpointer user_data) {

#if 0
   GtkWidget *chooser = widget_from_builder("cif_dictionary_filechooser_dialog");
   set_directory_for_filechooser(chooser);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, chooser);
   gtk_widget_show(chooser);
#endif

   GtkWindow *parent_window = GTK_WINDOW(user_data);
   GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_OPEN;
   GtkWidget *dialog = gtk_file_chooser_dialog_new("Open File", parent_window, action,
                                                   _("_Cancel"), GTK_RESPONSE_CANCEL,
                                                   _("_Open"), GTK_RESPONSE_ACCEPT,
                                                   NULL);
   g_signal_connect(dialog, "response", G_CALLBACK(on_cif_dictionary_filechooser_dialog_response_gtk4), NULL);
   g_object_set_data(G_OBJECT(dialog), "auto_read_flag", GINT_TO_POINTER(FALSE));
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.cif");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(dialog), filterselect);
   gtk_widget_show(dialog);

}

void toggle_display_frames_per_second_action(G_GNUC_UNUSED GSimpleAction *simple_action,
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
align_and_mutate_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_align_and_mutate_dialog();
   gtk_widget_show(w);
}


void
ligand_builder_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {

#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "FIXME:: start_ligand_builder_gui_internal() " << std::endl;
#else
   start_ligand_builder_gui_internal(menuitem, user_data);
#endif
}


void
run_script_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *widget = coot_run_script_chooser();
   add_filename_filter_button(widget, COOT_SCRIPTS_FILE_SELECTION);
   gtk_widget_show(widget);
}




void
lsq_superpose_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_least_squares_dialog(); // uses builder
   gtk_widget_show(w);
}

// ---------------- where is LSQ Plane? --------------------------


void
other_modelling_tools_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_other_model_tools_dialog();
   gtk_widget_show(w);
}



void
ssm_superposition_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_superpose_dialog(); // uses builder

   /* we get returned w = 0 when there is no MMDBSSM. (We are doing it
      this way because we don't have to introduce HAVE_MMDBSSM into the
      *c* compiler arguments (this is simpler)).  */
  if (w)
     gtk_widget_show(w);
}


void
calculate_updating_maps_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {
   show_calculate_updating_maps_gui();
}



void
background_colour_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {
   // black or white options
   std::cout << "black or white options here" << std::endl;
}


void
bond_colours_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                    G_GNUC_UNUSED GVariant *parameter,
                    G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = widget_from_builder("coords_colour_control_dialog");
   graphics_info_t g;
   g.fill_bond_colours_dialog_internal(w);
   gtk_widget_show(w);

}


void
bond_parameters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_bond_parameters_dialog(); // uses builder
   gtk_widget_show(w);
}


void
draw_cell_and_symmetry_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *show_symm_window = wrapped_create_show_symmetry_window();
   gtk_widget_show(show_symm_window);
}


void
display_only_active_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   display_only_active();
}

#include "cc-interface.hh" // for fullscreen()

void
fullscreen_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {
   fullscreen();
}


#include "generic-display-objects-c.h"
void
generic_objects_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {
   generic_objects_gui_wrapper();
}



void
go_to_atom_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *widget = wrapped_create_goto_atom_window(); // uses gtkbuilder

				/* now we need to fill the entry boxes
				   with default vaules and the option
				   menu according to molecules that
				   have coordinates. */

   gtk_widget_show(widget);
}


void
label_atoms_in_residue_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   label_atoms_in_residue();
}


void
label_CA_atoms_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {

   auto label_all_CAs = [] (int imol) {
      graphics_info_t::molecules[imol].add_labels_for_all_CAs();
   };

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      label_all_CAs(imol);
      graphics_draw();
   }
}


void
label_neighbours_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                        G_GNUC_UNUSED GVariant *parameter,
                        G_GNUC_UNUSED gpointer user_data) {
   label_neighbours();
}


void show_map_parameters_dialog(); // in glade-callbacks

void
map_parameters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                      G_GNUC_UNUSED GVariant *parameter,
                      G_GNUC_UNUSED gpointer user_data) {
   show_map_parameters_dialog();
}


void
ghost_control_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {
   GtkWidget *w = wrapped_create_ncs_control_dialog(); // uses builder
   gtk_widget_show(w);
}


void
spin_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   toggle_idle_spin_function();
}


void
rock_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                 G_GNUC_UNUSED GVariant *parameter,
                 G_GNUC_UNUSED gpointer user_data) {
   toggle_idle_rock_function();
}


void
ribbons_colour_by_chain_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                               G_GNUC_UNUSED GVariant *parameter,
                               G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "Chain";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}


void
ribbons_colour_rainbow_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "colorRampChainsScheme";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}


void
ribbons_colour_by_secondary_structure_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                             G_GNUC_UNUSED GVariant *parameter,
                                             G_GNUC_UNUSED gpointer user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "colorBySecondaryScheme";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}

void
screenshot_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {
    GtkWidget *file_chooser = coot_screendump_chooser();

    set_transient_and_position(COOT_UNDEFINED_WINDOW, file_chooser);
    g_object_set_data(G_OBJECT(file_chooser), "image_type", GINT_TO_POINTER(COOT_SCREENDUMP_SIMPLE));
    gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(file_chooser), "coot-screendump.tga");
    gtk_widget_show(file_chooser);
    check_for_dark_blue_density(); /* give a dialog if density it too dark (blue) */
}



void
sequence_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                     G_GNUC_UNUSED GVariant *parameter,
                     G_GNUC_UNUSED gpointer user_data) {

   // 20220628-PE this is different now?

   // add_on_sequence_view_choices()
   std::cout << "add on sequence_view options here " << std::endl;
}

void
undo_last_navigation_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

   undo_last_move();
}

void
undo_symmetry_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {
   undo_symmetry_view();
}

void
clear_atom_labels_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {
   remove_all_atom_labels();
}


void
distances_and_angles_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

  GtkWidget *widget = wrapped_create_geometry_dialog();
  set_transient_and_position(COOT_DISTANCES_ANGLES_WINDOW, widget);
  store_geometry_dialog(widget); /* needed to deactivate the distance
				    togglebutton after 2nd atoms
				    clicked in graphics */
  gtk_widget_show(widget);

}

void
pointer_distances_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                         G_GNUC_UNUSED GVariant *parameter,
                         G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = widget_from_builder("pointer_distances_dialog");
   fill_pointer_distances_widget(w);
   gtk_widget_show(w);
}

void
environment_distances_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *widget = widget_from_builder("environment_distance_dialog");
   fill_environment_widget(widget);
   gtk_widget_show(widget);

}


void
check_delete_waters_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                           G_GNUC_UNUSED GVariant *parameter,
                           G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_check_waters_dialog();
   int imol_map = imol_refinement_map();
   gtk_widget_show(w);
   if (imol_map < 0)
      show_select_map_dialog();

}


void
density_fit_analysis_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for density fit" << std::endl;
}


void
difference_map_peaks_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for difference map peaks fit" << std::endl;
}

void
geometry_analysis_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for geometry analysis" << std::endl;
}

void
gln_and_asn_b_factor_outlier_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                         G_GNUC_UNUSED GVariant *parameter,
                                         G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for GLN ASN" << std::endl;
}

void
chiral_volumes_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for chiral volumes" << std::endl;
}

void
ncs_differences_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                            G_GNUC_UNUSED GVariant *parameter,
                            G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for ncs differences" << std::endl;
}


void
peptide_flips_from_diff_map_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for peptide flips from difference map" << std::endl;
}

void
peptide_omega_analysis_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                   G_GNUC_UNUSED GVariant *parameter,
                                   G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for omega analysis " << std::endl;
}

void
pukka_puckers_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                          G_GNUC_UNUSED GVariant *parameter,
                          G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for Pukka Puckers" << std::endl;
}

void
ramachandran_plot_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                              G_GNUC_UNUSED GVariant *parameter,
                              G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for ramachandran plot" << std::endl;
}

void
rotamer_analysis_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                             G_GNUC_UNUSED GVariant *parameter,
                             G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for rotamer analysis" << std::endl;
}

void
temp_factor_analysis_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                 G_GNUC_UNUSED GVariant *parameter,
                                 G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for temperature factor analysis" << std::endl;
}

void
temp_factor_variance_analysis_item_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                          G_GNUC_UNUSED GVariant *parameter,
                                          G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus for temperature factor variance analysis" << std::endl;
}

void
unmodelled_blobs_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                                        G_GNUC_UNUSED GVariant *parameter,
                                        G_GNUC_UNUSED gpointer user_data) {

   std::cout << "dynamic menus unmodelled blobs" << std::endl;
}


void
remarks_browser_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                       G_GNUC_UNUSED GVariant *parameter,
                       G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *w = wrapped_create_remarks_browser_molecule_chooser_dialog();
   gtk_widget_show(w);
}


void
about_coot_action(G_GNUC_UNUSED GSimpleAction *simple_action,
                  G_GNUC_UNUSED GVariant *parameter,
                  G_GNUC_UNUSED gpointer user_data) {

   std::cout << "About Coot" << std::endl;

   GtkWidget *dialog = widget_from_builder("about_dialog");
   if (dialog) {
      gtk_widget_show(dialog);
   }
}

void
orthographic_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {
   set_use_perspective_projection(0);
}


void
perspective_view_action(G_GNUC_UNUSED GSimpleAction *simple_action,
G_GNUC_UNUSED GVariant *parameter,
G_GNUC_UNUSED gpointer user_data) {
   set_use_perspective_projection(1);
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

   // File

   add_action(     "open_coordinates_action",      open_coordinates_action);
   add_action(        "auto_open_mtz_action",         auto_open_mtz_action);
   add_action(         "open_dataset_action",          open_dataset_action);
   add_action(             "open_map_action",              open_map_action);
   add_action("import_cif_dictionary_action", import_cif_dictionary_action);

   add_action("get_monomer_action", get_monomer_action);
   add_action(     "curlew_action",      curlew_action);
   add_action(       "exit_action",        exit_action);

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

   // Edit

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

   // Calculate

   add_action(       "align_and_mutate_action",        align_and_mutate_action);
   add_action(         "ligand_builder_action",          ligand_builder_action);
   add_action(          "lsq_superpose_action",           lsq_superpose_action);
   add_action(             "run_script_action",              run_script_action);
   add_action(      "ssm_superposition_action",       ssm_superposition_action);
   add_action(  "other_modelling_tools_action",   other_modelling_tools_action);
   add_action("calculate_updating_maps_action", calculate_updating_maps_action);

   add_action("load_tutorial_model_and_data_action", load_tutorial_model_and_data_action);

   // Draw

   add_action(   "display_only_active_action",    display_only_active_action);
   add_action(     "background_colour_action",      background_colour_action);
   add_action(       "bond_parameters_action",        bond_parameters_action);
   add_action(          "bond_colours_action",           bond_colours_action);
   add_action(            "fullscreen_action",             fullscreen_action);
   add_action(            "go_to_atom_action",             go_to_atom_action);
   add_action(        "label_CA_atoms_action",         label_CA_atoms_action);
   add_action(        "map_parameters_action",         map_parameters_action);
   add_action(       "generic_objects_action",        generic_objects_action);
   add_action(      "label_neighbours_action",       label_neighbours_action);
   add_action("label_atoms_in_residue_action", label_atoms_in_residue_action);
   add_action("draw_cell_and_symmetry_action", draw_cell_and_symmetry_action);

   add_action(   "ghost_control_action",      ghost_control_action);
   add_action(        "spin_view_action",         spin_view_action);
   add_action(        "rock_view_action",         rock_view_action);
   add_action( "perspective_view_action",  perspective_view_action);
   add_action("orthographic_view_action", orthographic_view_action);

   add_action("screenshot_action",                           screenshot_action);
   add_action("sequence_view_action",                     sequence_view_action);
   add_action("undo_symmetry_view_action",           undo_symmetry_view_action);
   add_action("undo_last_navigation_action",       undo_last_navigation_action);
   add_action("ribbons_colour_by_chain_action", ribbons_colour_by_chain_action);
   add_action("ribbons_colour_rainbow_action",   ribbons_colour_rainbow_action);
   add_action("ribbons_colour_by_secondary_structure_action", ribbons_colour_by_secondary_structure_action);

   add_action( "toggle_display_frames_per_second_action", toggle_display_frames_per_second_action);

   // Measures

   add_action(    "clear_atom_labels_action",     clear_atom_labels_action);
   add_action(    "pointer_distances_action",     pointer_distances_action);
   add_action( "distances_and_angles_action",  distances_and_angles_action);
   add_action("environment_distances_action", environment_distances_action);

   // Validate

   add_action(                  "unmodelled_blobs_action",                   unmodelled_blobs_action);
   add_action(               "check_delete_waters_action",                check_delete_waters_action);
   add_action(         "density_fit_analysis_item_action",          density_fit_analysis_item_action);
   add_action(         "difference_map_peaks_item_action",          difference_map_peaks_item_action);
   add_action(            "geometry_analysis_item_action",             geometry_analysis_item_action);
   add_action(               "chiral_volumes_item_action",                chiral_volumes_item_action);
   add_action(              "ncs_differences_item_action",               ncs_differences_item_action);
   add_action(                "pukka_puckers_item_action",                 pukka_puckers_item_action);
   add_action(            "ramachandran_plot_item_action",             ramachandran_plot_item_action);
   add_action(             "rotamer_analysis_item_action",              rotamer_analysis_item_action);
   add_action(         "temp_factor_analysis_item_action",          temp_factor_analysis_item_action);
   add_action(       "peptide_omega_analysis_item_action",        peptide_omega_analysis_item_action);
   add_action(  "peptide_flips_from_diff_map_item_action",   peptide_flips_from_diff_map_item_action);
   add_action( "gln_and_asn_b_factor_outlier_item_action",  gln_and_asn_b_factor_outlier_item_action);
   add_action("temp_factor_variance_analysis_item_action", temp_factor_variance_analysis_item_action);

   // About

   add_action("remarks_browser_action", remarks_browser_action);
   add_action("about_coot_action", about_coot_action);
}

