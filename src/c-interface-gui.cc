/* src/c-interface-gui.cc
 *
 * Copyright 2003, 2004, 2007 The University of York
 * Copyright 2008 The University of Oxford
 * Author: Paul Emsley
 * Copyright 2007 The University of York
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 * Author: Bernhard Lohkamp
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#ifndef PYTHONH
#include <Python.h>
#endif
#endif

#include "compat/coot-sysdep.h"


#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include <string.h> // strlen, strncpy
#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define snprintf _snprintf
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

#if !defined(_MSC_VER)
#include <glob.h> // for globbing.
#endif

#include "guile-fixups.h"

#include "c-interface-scm.hh"
#include "coot-fileselections.h"
#include "coot-references.h"
#include "coot-preferences.h"
#include "rotate-translate-modes.hh"
#include "curl-utils.hh"

#include "graphics-info.h"
#include "interface.h"
#include "c-interface.h"
#include "c-interface-generic-objects.h"
#include "c-interface-gtk-widgets.h"
#include "c-interface-preferences.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "file-system-utils.hh"
#include "cmtz-interface.hh"
#include "cmtz-interface-gui.hh"
#include "coords/mmdb.hh"  // for centre of molecule
#include "clipper/core/clipper_instance.h"

#include "c-interface-gui.hh"
#include "utils/win-compat.hh"

#include "c-interface-widgets.hh"

#include "widget-from-builder.hh"
#include "utils/coot-utils.hh"
#include "read-molecule.hh"

#include "utils/logging.hh"
extern logging logger;


void set_show_paths_in_display_manager(int i) {
   std::string cmd = "set-show-paths-in-display-manager";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
   graphics_info_t::show_paths_in_display_manager_flag = i;
}

int show_paths_in_display_manager_state() {
   add_to_history_simple("show-paths-in-display-manager-state");
   return graphics_info_t::show_paths_in_display_manager_flag;
}

/*! \brief set the GUI dark mode state
 */
void set_use_dark_mode(short int state) {

   if (state)
      g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", TRUE, NULL);
   else
      g_object_set(gtk_settings_get_default(), "gtk-application-prefer-dark-theme", FALSE, NULL);

}



/*! \brief display the open coordinates dialog */
void open_coords_dialog() {

   if (graphics_info_t::use_graphics_interface_flag) {

      /* This split was here because the buttons don't work. They act on the
	 file list, using the file list as a CList.  And CList is deprecated
	 in GTk+2.  So the button-press callback code needs to be adjusted. */
      GtkWidget *coords_filechooser = coot_file_chooser();

      // add_ccp4i_project_optionmenu(coords_filechooser1, COOT_COORDS_FILE_SELECTION);

      // GtkWidget *file_filter_button = add_filename_filter_button(coords_filechooser, COOT_COORDS_FILE_SELECTION);

      add_filechooser_filter_button(coords_filechooser, COOT_COORDS_FILE_SELECTION);

      // sort_button = add_sort_button_filechooser(coords_filechooser1); // FIXME
      add_recentre_on_read_pdb_combobox(coords_filechooser);
      set_directory_for_coot_file_chooser(coords_filechooser);
      set_file_selection_dialog_size(coords_filechooser);
      set_transient_and_position(COOT_UNDEFINED_WINDOW, coords_filechooser);
      gtk_widget_set_visible (coords_filechooser, TRUE);

      /* in gtk2 we have to push the buttons after we show the selection */
      // push_the_buttons_on_filechooser(file_filter_button, sort_button, coords_filechooser1);
      // FIXME

   }
}


void
open_cif_dictionary_file_selector_dialog() {

   if (graphics_info_t::use_graphics_interface_flag) {

      GtkWidget *filechooser = coot_cif_dictionary_chooser(); // a chooser or a fileselection
      add_sort_button_fileselection(filechooser);
      set_directory_for_coot_file_chooser(filechooser);
      set_file_selection_dialog_size(filechooser);

      // add_ccp4i_project_optionmenu(fileselection, COOT_CIF_DICTIONARY_FILE_SELECTION);
      // add_filename_filter_button(fileselection, COOT_CIF_DICTIONARY_FILE_SELECTION);
      // add_sort_button_fileselection(fileselection);
      // set_directory_for_fileselection(fileselection);
      // set_file_selection_dialog_size(fileselection);

      if (true) {

	 GtkWidget *aa_hbutton_box = gtk_dialog_get_content_area(GTK_DIALOG(filechooser));
	 // if (GTK_IS_HBUTTON_BOX(aa_hbutton_box)) {
         if (aa_hbutton_box) { // check the button type
	    add_cif_dictionary_selector_molecule_selector(filechooser, aa_hbutton_box);
	    add_cif_dictionary_selector_create_molecule_checkbutton(filechooser, aa_hbutton_box);
	 }
      }
      gtk_widget_set_visible(filechooser, TRUE);
   }
}

void
add_cif_dictionary_selector_molecule_selector(GtkWidget *fileselection, // maybe it's a chooser
					      GtkWidget *aa_hbox) {

   std::cout << "GTK-FIXME --- delete this function add_cif_dictionary_selector_molecule_selector"
	     << std::endl;

}


void cif_dictionary_molecule_menu_item_select(GtkWidget *item, GtkPositionType pos) {

   // pos is the value stored in with GINT_TO_POINTER() in the signal connect.
   //
   // std::cout << "select menu item " << item << " pos " << pos << std::endl;
}


void
add_cif_dictionary_selector_create_molecule_checkbutton(GtkWidget *fileselection,
							GtkWidget *aa_hbox) {

   // if we came from a chooser, aa_hbox is an hbutton_box
   // if we came from a selector, aa_hbox is an hbox.

   GtkWidget *frame = gtk_frame_new("Make a Molecule");
   GtkWidget *checkbutton = gtk_check_button_new_with_label(" Generate a Molecule");
   g_object_set_data_full(G_OBJECT(fileselection),
			  "cif_dictionary_file_selector_create_molecule_checkbutton",
			  checkbutton, // remove gtk_widget_ref
			  (GDestroyNotify) NULL);

   graphics_info_t g;

   if (g.cif_dictionary_file_selector_create_molecule_flag)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);

   // do we need to connect this signal?
   GCallback callback_func =
      G_CALLBACK(on_cif_dictionary_file_selector_create_molecule_checkbutton_toggled);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   gtk_box_append(GTK_BOX(aa_hbox), frame);
#else
   gtk_box_pack_start(GTK_BOX(aa_hbox), frame, FALSE, TRUE, 0);
#endif
   gtk_frame_set_child(GTK_FRAME(frame), checkbutton);
   gtk_widget_set_visible(checkbutton, TRUE);
   gtk_widget_set_visible(frame, TRUE);
 }



void
on_cif_dictionary_file_selector_create_molecule_checkbutton_toggled (GtkButton       *button,
								     gpointer         user_data) {

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
      std::cout << "Make a molecule after dictionary" << std::endl;
   } else {
      std::cout << "on_cif_dictionary_file_selector_create_molecule_checkbutton_toggled() "
		<< "Do nothing" << std::endl;
   }
}



GtkWidget *wrapped_nothing_bad_dialog(const std::string &label) {

   graphics_info_t g;
   return g.wrapped_nothing_bad_dialog(label);
}


#include "widget-headers.hh"
#include "widget-from-builder.hh"

GtkWidget *wrapped_create_remarks_browser_molecule_chooser_dialog() {

   // GtkWidget *w = create_remarks_browser_molecule_chooser_dialog();
   GtkWidget *w = widget_from_builder("remarks_browser_molecule_chooser_dialog");
   fill_remarks_browswer_chooser(w);
   return w;
}

void fill_remarks_browswer_chooser(GtkWidget *w) {

   auto get_model_molecule_vector = [] () {
                                       graphics_info_t g;
                                       std::vector<int> vec;
                                       int n_mol = g.n_molecules();
                                       for (int i=0; i<n_mol; i++)
                                          if (g.is_valid_model_molecule(i))
                                             vec.push_back(i);
                                       return vec;
                                    };

   GtkWidget *combobox = widget_from_builder("remarks_browser_molecule_chooser_combobox_text");
   if (combobox) {
      graphics_info_t g;
      gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox));
      // GCallback callback_func = G_CALLBACK(remarks_browswer_molecule_item_select);
      GCallback callback_func = G_CALLBACK(remarks_browswer_molecule_combobox_changed);
      int imol_active = first_coords_imol();
      g.imol_remarks_browswer = imol_active;
      auto mv = get_model_molecule_vector();
      g.fill_combobox_with_molecule_options(combobox, callback_func, imol_active, mv);
   } else {
      std::cout << "fill_remarks_browswer_chooser() failed to get combobox" << std::endl;
   }
}

void remarks_browswer_molecule_combobox_changed(GtkWidget *combobox, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::imol_remarks_browswer = imol;
}



void remarks_browswer_molecule_item_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t::imol_remarks_browswer = pos;
}

void show_remarks_browswer() {

   if (graphics_info_t::use_graphics_interface_flag) {
      remarks_dialog(graphics_info_t::imol_remarks_browswer);
   }
}

// To be used to (typically) get the menu item text label from chain
// option menus (rather than the ugly/broken casting of
// GtkPositionType data.  A wrapper to a static graphics_info_t
// function.
std::string menu_item_label(GtkWidget *menu_item) {


//    char *data = NULL;
//    data = (char *)pos;
//    // this can fail when more than one sequence mutate is used at the same time:
//    if (data)
//       graphics_info_t::superpose_imol1_chain = data;

   return graphics_info_t::menu_item_label(menu_item);

}



/*! \brief show the Undo Molecule chooser - i.e. choose the molecule
  to which the "Undo" button applies. */
void show_set_undo_molecule_chooser() {

   GtkWidget *w = wrapped_create_undo_molecule_chooser_dialog();
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);

}

GtkWidget *wrapped_create_undo_molecule_chooser_dialog() {

   GtkWidget *dialog         = widget_from_builder("undo_molecule_chooser_dialog");
   GtkWidget *model_combobox = widget_from_builder("undo_molecule_chooser_comboboxtext");
   graphics_info_t g;

   // g.fill_combobox_with_undo_options(combobox);

   auto get_model_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_model_molecule(i))
                                           vec.push_back(i);
                                     return vec;
                                  };

   auto combobox_changed_func = +[] (GtkWidget *combobox, gpointer user_data) {
      graphics_info_t g;
      int imol_coords = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
      g.set_undo_molecule_number(imol_coords);
   };

   int imol_active = g.Undo_molecule(coot::UNDO);
   GCallback func = G_CALLBACK(combobox_changed_func);
   auto model_list = get_model_molecule_vector();
   g.fill_combobox_with_molecule_options(model_combobox, func, imol_active, model_list);
   return dialog;
}


/* We try as .phs and .cif files.

   Called after we select the mtz filename, sets up the column label
   widget and displays it. */
void
manage_column_selector(const char *filename) {

   if (graphics_info_t::use_graphics_interface_flag) {
      // try to read as phs, cif etc, if not, return a selection
      // widget

      GtkWidget *w = coot::column_selector_using_cmtz(filename);

      if (w) {
         gtk_widget_set_visible(w, TRUE);
         gtk_window_present(GTK_WINDOW(w));
      }
   }

   std::string cmd = "manage-column-selector";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(filename));
   add_to_history_typed(cmd, args);

}

void
manage_refmac_column_selection(GtkWidget *run_refmac_dialog) {

   // called by an mtz file chooser response

   if (graphics_info_t::use_graphics_interface_flag) {
     coot::setup_refmac_parameters_from_file(run_refmac_dialog);
   }
}

void
store_refmac_mtz_file_label(GtkWidget *w) {

  graphics_info_t::refmac_dialog_mtz_file_label = w;
}

GtkWidget *get_refmac_mtz_file_label() {

  return graphics_info_t::refmac_dialog_mtz_file_label;
}

// we want to have an interface to save refmac parameters in map objects,
// so that we can save the original mtz file and labels in a map file
void save_refmac_params_to_map(int imol_map,
                               const char *mtz_filename,
                               const char *fobs_col,
                               const char *sigfobs_col,
                               const char *r_free_col,
                               int r_free_flag_sensible) {

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t::molecules[imol_map].store_refmac_params(std::string(mtz_filename),
                                                               std::string(fobs_col),
                                                               std::string(sigfobs_col),
                                                               std::string(r_free_col),
                                                               r_free_flag_sensible);
   } else {
      std::cout << "WARNGING:: invalid map molecule!" <<std::endl;
   }

}

void save_refmac_phase_params_to_map(int imol_map,
                                     const char *phi,
                                     const char *fom,
                                     const char *hla,
                                     const char *hlb,
                                     const char *hlc,
                                     const char *hld) {

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t::molecules[imol_map].store_refmac_phase_params(std::string(phi),
                                                                     std::string(fom),
                                                                     std::string(hla),
                                                                     std::string(hlb),
                                                                     std::string(hlc),
                                                                     std::string(hld));
   } else {
      std::cout << "WARNGING:: invalid map molecule!" <<std::endl;
   }

}

// get string for column 0 (which are strings)
std::string
get_active_label_in_combobox(GtkComboBox *combobox) {

   graphics_info_t g;
   return g.get_active_label_in_combobox(combobox);
}

void handle_column_label_make_fourier_v2(GtkWidget *column_label_window) {

   // does any of this refmac code work now?

   std::cout << ":::::::::::::::::::::::: handle_column_label_make_fourier_v2() " << std::endl;

   // GtkWidget *weights_checkbutton     = lookup_widget(GTK_WIDGET(column_label_window), "use_weights_checkbutton");
   // GtkWidget *is_diff_map_checkbutton = lookup_widget(GTK_WIDGET(column_label_window), "difference_map_checkbutton");
   // GtkWidget *reso_limit_checkbutton  = lookup_widget(GTK_WIDGET(column_label_window), "column_labels_use_resolution_limits_checkbutton");
   GtkWidget *weights_checkbutton     = widget_from_builder("use_weights_checkbutton");
   GtkWidget *is_diff_map_checkbutton = widget_from_builder("difference_map_checkbutton");
   GtkWidget *reso_limit_checkbutton  = widget_from_builder("column_labels_use_resolution_limits_checkbutton");

   bool use_weights_flag       = false;
   bool is_difference_map_flag = false;
   bool limit_reso_flag = false;

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(weights_checkbutton)))     use_weights_flag       = true;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(is_diff_map_checkbutton))) is_difference_map_flag = true;

   // GtkWidget *amplitudes_combobox = lookup_widget(column_label_window, "column_selector_amplitudes_combobox");
   // GtkWidget *phases_combobox     = lookup_widget(column_label_window, "column_selector_phases_combobox");
   // GtkWidget *weights_combobox    = lookup_widget(column_label_window, "column_selector_weights_combobox");
   GtkWidget *amplitudes_combobox = widget_from_builder("column_selector_amplitudes_combobox");
   GtkWidget *phases_combobox     = widget_from_builder("column_selector_phases_combobox");
   GtkWidget *weights_combobox    = widget_from_builder("column_selector_weights_combobox");

   std::string f_label   = get_active_label_in_combobox(GTK_COMBO_BOX(amplitudes_combobox));
   std::string phi_label = get_active_label_in_combobox(GTK_COMBO_BOX(phases_combobox));
   std::string w_label;
   if (use_weights_flag)
      w_label = get_active_label_in_combobox(GTK_COMBO_BOX(weights_combobox));

   std::string fobs_col, sigfobs_col, r_free_col;
   bool have_refmac_params = false;
   bool sensible_r_free_col = false;
   bool is_anomalous_flag = false;
   float low_res_limit  = -1.0;
   float high_res_limit = -1.0;

   /* --------- Resolution limits --------- */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(reso_limit_checkbutton))) {

      //GtkEntry  *low_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(column_label_window), "column_labels_reso_low_entry"));
      // GtkEntry *high_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(column_label_window), "column_labels_reso_high_entry"));

      GtkEntry  *low_entry = GTK_ENTRY(widget_from_builder("column_labels_reso_low_entry"));
      GtkEntry *high_entry = GTK_ENTRY(widget_from_builder("column_labels_reso_high_entry"));

      std::string l = gtk_editable_get_text(GTK_EDITABLE(low_entry));
      std::string h = gtk_editable_get_text(GTK_EDITABLE(high_entry));
      bool low_OK = true;
      bool high_OK = true;

      // It's OK not to set a low resolution limit, but not OK to not to set a high resolution limit
      // (if that is the case, act as if no resolution limit was enabled)

      if (! l.empty()) {
         try {
            float ll = coot::util::string_to_float(l);
            low_res_limit = ll;
         }
         catch (const std::runtime_error &rte) {
            std::cout << "WARNING:: " << rte.what() << std::endl;
            low_OK = false;
         }
      } else {
         low_res_limit = 9999.9;
      }

      if (! h.empty()) {
         try {
            float hh = coot::util::string_to_float(h);
            high_res_limit = hh;
         }
         catch (const std::runtime_error &rte) {
            std::cout << "WARNING:: " << rte.what() << std::endl;
         }
      } else {
         high_OK = false;
      }

      if (low_OK && high_OK)
         limit_reso_flag = true;

      if (false) {// debugging, force non-sane limits
         low_res_limit  = 999.9;
         high_res_limit = 999.9;
         limit_reso_flag = true;
      }
   }


   /* --------- Refmac label stuff --------- */

   // GtkWidget *refmac_columns_checkbutton = lookup_widget(GTK_WIDGET(column_label_window), "refmac_column_labels_checkbutton");
   GtkWidget *refmac_columns_checkbutton = widget_from_builder("refmac_column_labels_checkbutton");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(refmac_columns_checkbutton))) {
      have_refmac_params = 1;
      // GtkWidget *fobs_combobox    = lookup_widget(column_label_window, "column_label_selector_refmac_fobs_combobox");
      // GtkWidget *sigfobs_combobox = lookup_widget(column_label_window, "column_label_selector_refmac_sigfobs_combobox");
      // GtkWidget *r_free_combobox  = lookup_widget(column_label_window, "column_label_selector_refmac_rfree_combobox");
      GtkWidget *fobs_combobox    = widget_from_builder("column_label_selector_refmac_fobs_combobox");
      GtkWidget *sigfobs_combobox = widget_from_builder("column_label_selector_refmac_sigfobs_combobox");
      GtkWidget *r_free_combobox  = widget_from_builder("column_label_selector_refmac_rfree_combobox");

      fobs_col    = get_active_label_in_combobox(GTK_COMBO_BOX(   fobs_combobox));
      sigfobs_col = get_active_label_in_combobox(GTK_COMBO_BOX(sigfobs_combobox));
      r_free_col  = get_active_label_in_combobox(GTK_COMBO_BOX( r_free_combobox));
      if (! r_free_col.empty()) sensible_r_free_col = true;
   }

   /* --------- Save the column data in the attached pointer --------- */

   std::string mtz_filename;
   gpointer d = g_object_get_data(G_OBJECT(column_label_window), "f_phi_columns"); // set in column_selector_using_cmtz()
   coot::mtz_column_types_info_t *saved_f_phi_columns = 0;
   if (d) {
      saved_f_phi_columns = static_cast<coot::mtz_column_types_info_t *> (d);
      mtz_filename = saved_f_phi_columns->mtz_filename;
   }

   /* --------- make and draw --------- */

   if (false)
      std::cout << "debug:: calling make_and_draw_map_with_reso_with_refmac_params() with the refmac params "
                << have_refmac_params << " " << fobs_col << " " << sigfobs_col << " "
                << r_free_col << " " << sensible_r_free_col << std::endl;

   make_and_draw_map_with_reso_with_refmac_params(mtz_filename.c_str(),
                                                  f_label.c_str(),
                                                  phi_label.c_str(),
                                                  w_label.c_str(),
                                                  use_weights_flag, is_difference_map_flag,
                                                  have_refmac_params,
                                                  fobs_col.c_str(),
                                                  sigfobs_col.c_str(),
                                                  r_free_col.c_str(),
                                                  sensible_r_free_col,
                                                  is_anomalous_flag,
                                                  limit_reso_flag,
                                                  low_res_limit, high_res_limit);

   gtk_widget_set_visible(column_label_window, FALSE);

}

void handle_column_label_make_fourier(GtkWidget *column_label_window) {

   if (false)
      std::cout << "---- handle_column_label_make_fourier() with column_label_window "
                << column_label_window << std::endl;

   GtkWidget *refmac_checkbutton;
   int icol;
   int use_weights = 0;
   int is_diff_map;
   short int sensible_r_free_col = 0;
   short int have_refmac_params = 0; /* default not */
   short int use_resolution_limits_flag = 0;
   float low_reso_lim = -1.0;  /* unset */
   float high_reso_lim = -1.0; /* unset */
   short int is_anomalous_flag = 0;

   GtkWidget *fobs_option_menu;
   GtkWidget *sigfobs_option_menu;
   GtkWidget *r_free_option_menu;

   GtkWidget *fobs_menu;
   GtkWidget *sigfobs_menu;
   GtkWidget *r_free_menu;

   GtkCheckButton *check_weights;
   GtkCheckButton *is_diff_map_checkbutton;
   GtkCheckButton *resolution_limit_check_button;
   GtkEntry *low_entry;
   GtkEntry *high_entry;


   /* Was the "Use Weights checkbutton clicked?  */

   // check_weights = GTK_CHECK_BUTTON(lookup_widget(GTK_WIDGET(column_label_window), "use_weights_checkbutton"));
   check_weights = GTK_CHECK_BUTTON(widget_from_builder("use_weights_checkbutton"));


   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_weights))) {
      use_weights = 1;
   } else {
      use_weights = 0;
   }

   /* Similarly, we ask, was the "Is difference map" checkbutton active?  */

   is_diff_map_checkbutton = GTK_CHECK_BUTTON(widget_from_builder("difference_map_checkbutton"));

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton))) {
      is_diff_map = 1;
   } else {
      is_diff_map = 0;
   }
   void *t = g_object_get_data(G_OBJECT(column_label_window), "f_phi_columns");
   coot::mtz_column_types_info_t *saved_f_phi_columns = static_cast<coot::mtz_column_types_info_t *> (t);

   if (! saved_f_phi_columns)
      return;

   const char *object_mtz_filename = saved_f_phi_columns->mtz_filename.c_str();

   /* Get the values that the user has selected in the option menu
      buttons. */

   {

      // GtkWidget *amplitudes_combobox = lookup_widget(column_label_window, "column_selector_amplitudes_combobox");
      // GtkWidget *phases_combobox     = lookup_widget(column_label_window, "column_selector_phases_combobox");
      // GtkWidget *weights_combobox    = lookup_widget(column_label_window, "column_selector_weights_combobox");
      GtkWidget *amplitudes_combobox = widget_from_builder("column_selector_amplitudes_combobox");
      GtkWidget *phases_combobox     = widget_from_builder("column_selector_phases_combobox");
      GtkWidget *weights_combobox    = widget_from_builder("column_selector_weights_combobox");

      std::string phi_label;
      std::string f_label;
      std::string w_label;
      std::string fobs_col;
      std::string sigfobs_col;
      std::string r_free_col;

      f_label = get_active_label_in_combobox(GTK_COMBO_BOX(amplitudes_combobox));
      phi_label = get_active_label_in_combobox(GTK_COMBO_BOX(phases_combobox));

      if (use_weights) {
	 w_label = get_active_label_in_combobox(GTK_COMBO_BOX(weights_combobox));
	 std::cout << " Making map from " << f_label << " " << phi_label << " and "
		   << w_label << std::endl;
      } else {
	 std::cout << " Making map from " << f_label << " and " << phi_label << std::endl;
      }

      /* is the resolution limit check button in use? */
      // resolution_limit_check_button = GTK_CHECK_BUTTON(lookup_widget(GTK_WIDGET(column_label_window), "column_labels_use_resolution_limits_checkbutton"));
      resolution_limit_check_button = GTK_CHECK_BUTTON(widget_from_builder("column_labels_use_resolution_limits_checkbutton"));
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(resolution_limit_check_button))) {

	 /* yes, it is.. */

	 // low_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(column_label_window), "column_labels_reso_low_entry"));
	 // high_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(column_label_window), "column_labels_reso_high_entry"));
	 low_entry  = GTK_ENTRY(widget_from_builder("column_labels_reso_low_entry"));
	 high_entry = GTK_ENTRY(widget_from_builder("column_labels_reso_high_entry"));

	 low_reso_lim  = get_positive_float_from_entry(low_entry);
	 high_reso_lim = get_positive_float_from_entry(high_entry);
	 std::cout << "Resolution limits: low: " << low_reso_lim << " and high: "
		   << high_reso_lim << std::endl;
	 if (high_reso_lim > 0.0001)
	    use_resolution_limits_flag = 1;
	 // if low_reso_lim is not set, presume that it was 999.9;
	 if (low_reso_lim < 0.0)
	    low_reso_lim = 999.9;
      }

      /* Refmac label stuff */

      // refmac_checkbutton = lookup_widget(GTK_WIDGET(column_label_window), "refmac_column_labels_checkbutton");
      refmac_checkbutton = widget_from_builder("refmac_column_labels_checkbutton");

      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(refmac_checkbutton))) {

	 have_refmac_params = 1;

	 // GtkWidget *fobs_combobox    = lookup_widget(column_label_window, "column_label_selector_refmac_fobs_combobox");
	 // GtkWidget *sigfobs_combobox = lookup_widget(column_label_window, "column_label_selector_refmac_sigfobs_combobox");
	 // GtkWidget *r_free_combobox  = lookup_widget(column_label_window, "column_label_selector_refmac_rfree_combobox");
	 GtkWidget *fobs_combobox    = widget_from_builder("column_label_selector_refmac_fobs_combobox");
	 GtkWidget *sigfobs_combobox = widget_from_builder("column_label_selector_refmac_sigfobs_combobox");
	 GtkWidget *r_free_combobox  = widget_from_builder("column_label_selector_refmac_rfree_combobox");

	 fobs_col    = get_active_label_in_combobox(GTK_COMBO_BOX(fobs_combobox));
	 sigfobs_col = get_active_label_in_combobox(GTK_COMBO_BOX(fobs_combobox));
	 r_free_col  = get_active_label_in_combobox(GTK_COMBO_BOX(fobs_combobox));
      }

      std::cout << "---------------------- Here" << std::endl;

      /* And proceed with the actual map-making.
	 If use_weights is 1, then weights should be used.*/
      make_and_draw_map_with_reso_with_refmac_params(object_mtz_filename,
						     f_label.c_str(),
						     phi_label.c_str(),
						     w_label.c_str(),
						     use_weights, is_diff_map,
						     have_refmac_params,
						     fobs_col.c_str(),
						     sigfobs_col.c_str(),
						     r_free_col.c_str(),
						     sensible_r_free_col,
						     is_anomalous_flag,
						     use_resolution_limits_flag,
						     low_reso_lim, high_reso_lim);
   }
   /* We can destroy the column_label_window top level widget now. */
   gtk_widget_set_visible(column_label_window, FALSE);

}


void fill_combobox_with_expert_options(GtkWidget *amplitudes_combobox) {

   // data for this widget is set in... coot::column_selector_using_cmtz(const std::string &filename)

   // GtkWidget *column_label_window = lookup_widget(amplitudes_combobox, "column_label_window");
   GtkWidget *column_label_window = widget_from_builder("column_label_window");
   coot::mtz_column_types_info_t *saved_f_phi_columns
      = static_cast<coot::mtz_column_types_info_t *> (g_object_get_data(G_OBJECT(column_label_window),
									"f_phi_columns"));
   if (saved_f_phi_columns) {
      const coot::mtz_column_types_info_t &col_labs = *saved_f_phi_columns;
      int f_prefered_idx = col_labs.get_prefered_f_col_idx();
      std::vector<coot::mtz_type_label> labels = col_labs.f_cols;
      std::vector<coot::mtz_type_label> d_labels = col_labs.d_cols;
      labels.insert(labels.end(), d_labels.begin(), d_labels.end());
      my_combo_box_text_add_items(GTK_COMBO_BOX(amplitudes_combobox), labels, f_prefered_idx);
   } else {
      std::cout << "failed to lookup" << std::endl;
   }

}


void fill_about_window(GtkWidget *widget) {

   GtkWidget *text_widget = widget_from_builder("about_window_text");

   std::string body_text("\n\n   Brought to you by:\n\n   Paul Emsley & Kevin Cowtan\n\n   Using the dictionaries of:\n    Alexei Vagin\n");

#ifdef USE_DUNBRACK_ROTAMERS
   body_text += "    Roland Dunbrack & co-workers\n\n";
#else
   body_text += "    Jane and David Richardson\n";
   body_text += "    & co-workers\n\n";
#endif

#ifdef WINDOWS_MINGW
   body_text += "    Ported to Windows for you by:\n";
   body_text += "    Bernhard Lohkamp\n\n";
#endif // MINGW

   body_text += "  Using the libraries of:\n   Eugene Krissinel\n   Kevin Cowtan\n   Stuart McNicholas\n   Ralf W. Grosse-Kunstleve\n   Janne Lof\n   Raghavendra Chandrashekara\n   Paul Bourke & Cory Gene Bloyd\n   Matteo Frigo & Steven G. Johnson\n   & many others.\n\n  Windows 2000 Binaries\n   Bernhard Lohkamp\n\n  Macintosh Binaries\n   William Scott\n\n";

   std::string widget_text("\n   Coot version ");
   widget_text += VERSION;
   widget_text += body_text;

  gtk_text_view_set_editable (GTK_TEXT_VIEW (text_widget), FALSE);
  gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (text_widget), GTK_WRAP_WORD);
  gtk_text_buffer_set_text (gtk_text_view_get_buffer (GTK_TEXT_VIEW (text_widget)),
			    widget_text.c_str(), -1);

}

void add_coot_references_button(GtkWidget *widget) {

   if (! widget) return;

   // hbox = GTK_DIALOG(widget)->action_area;
   GtkWidget *hbox = gtk_dialog_get_header_bar(GTK_DIALOG(widget));
   GtkWidget *button = gtk_button_new_with_label("References");
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE-FIXME box packing and reorder and secondary                                                 
   gtk_box_append(GTK_BOX(hbox), button);
#else
   gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, TRUE, 0);
   gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox), button, TRUE);
   gtk_box_reorder_child(GTK_BOX(hbox), button, 2);
#endif

   g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(wrapped_create_coot_references_dialog), NULL);
   gtk_widget_set_visible(button, TRUE);

}

GtkWidget *wrapped_create_coot_references_dialog() {

  GtkWidget *coot_reference_button;
  // GtkWidget *references_dialog = create_coot_references_dialog();
  GtkWidget *references_dialog = widget_from_builder("coot_references_dialog");
  // coot_reference_button = lookup_widget(references_dialog, "coot_references_coot_toolbutton");
  coot_reference_button = widget_from_builder("coot_references_coot_toolbutton");
  g_signal_emit_by_name(G_OBJECT(coot_reference_button), "clicked");
  gtk_widget_set_visible(references_dialog, TRUE);
  return references_dialog;

}


void fill_references_notebook(GtkButton *toolbutton, int reference_id) {

  GtkWidget *notebook;
  GtkWidget *ref_text_view;
  GtkWidget *bib_text_view;
  GtkTextBuffer *ref_buffer;
  GtkTextBuffer *bib_buffer;
  GtkTextIter end_iter;
  std::string ref_description;
  std::string ref_text;
  std::string bib_text;
  std::string title;
  std::string author;
  std::string journal;
  std::string year;
  std::string volume;
  std::string number;
  std::string pages;
  std::string bib_type;
  std::string bib_id;
  std::string bib_title;
  std::string bib_author;
  std::string bib_journal;

  notebook      = widget_from_builder("coot_references_notebook");
  ref_text_view = widget_from_builder("coot_references_textview");
  bib_text_view = widget_from_builder("coot_bibtext_textview");

  ref_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(ref_text_view));
  bib_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(bib_text_view));
  //  gtk_widget_unref(GTK_WIDGET(ref_buffer));
  //  gtk_widget_unref(GTK_WIDGET(bib_buffer));

  ref_buffer = gtk_text_buffer_new(NULL);
  bib_buffer = gtk_text_buffer_new(NULL);

  gtk_text_view_set_buffer(GTK_TEXT_VIEW(ref_text_view), ref_buffer);
  gtk_text_view_set_buffer(GTK_TEXT_VIEW(bib_text_view), bib_buffer);
  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(ref_text_view), GTK_WRAP_WORD);

  // simple text fill for now
  ref_text = "dummy";
  bib_text = "dummy";
  ref_description = "dummy";
  bib_title = "dummy";
  bib_author = "dummy";
  bib_journal = "dummy";

  if (reference_id == COOT_REFERENCE_COOT) {

     if (0) {
	ref_description = "If you have found this software to be useful, you are requested to cite:\n\n";

	title       = "Coot: model-building tools for molecular graphics";
	author      = "Emsley P, Cowtan K";
	journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
	year        = "2004";
	volume      = "60";
	pages       = "2126-2132";
	number      = "Part 12 Sp. Iss. 1 DEC";

	bib_type = "Article";
	bib_id   = "emsley04:coot";
	bib_author = "Paul Emsley and Kevin Cowtan";
	bib_journal = "Acta Crystallographica Section D - Biological Crystallography";
     }

     ref_description = "If you have found this software to be useful, you are requested to cite:\n\n";
     title       = "Features and Development of Coot";
     author      = "Emsley P, Lohkamp B, Scott W, Cowtan K";
     journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
     year        = "2010";
     volume      = "66";
     pages       = "486-501";
     number      = "4";

     bib_type = "Article";
     bib_id   = "emsley10:coot";
     bib_author = "Paul Emsley, Bernhard Lohkamp, William Scott and Kevin Cowtan";
     bib_journal = "Acta Crystallographica Section D - Biological Crystallography";

  }

  if (reference_id == COOT_REFERENCE_WINCOOT) {

    ref_description = "Please cite as for Coot now. You find additional information on WinCoot in the following. Feel free to cite as well as above:\n\n";

    title      = "Coot News";
    author     = "Lohkamp B, Emsley P, Cowtan K";
    journal    = "CCP4 Newsletter";
    year       = "2005";
    volume     = "42";
    pages      = "";
    number     = "Contribution 7";

    bib_type   = "Article";
    bib_id     = "lohkamp05:wincoot";
    bib_author = "Bernhard Lohkamp and Paul Emsley and Kevin Cowtan";

  }
  if (reference_id == COOT_REFERENCE_REFMAC) {

    ref_description = "The reference for the REFMAC5 Dictionary is:\n\n";

    title      = "REFMAC5 dictionary: organization of prior chemical knowledge and guidelines for its use";
    author     = "Vagin AA, Steiner RA, Lebedev AA, Potterton L, McNicholas S, Long F, Murshudov GN";
    journal    = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year       = "2004";
    volume     = "60";
    number     = "12 Part 1";
    pages      = "2184-2195";

    bib_type   = "Article";
    bib_id     = "Vagin:ba5073";
    bib_title  = "{\\it REFMAC}5 dictionary: organization of prior chemical knowledge and guidelines for its use";
    bib_author = "Vagin, Alexei A. and Steiner, Roberto A. and Lebedev, Andrey A. and Potterton, Liz and McNicholas, Stuart and Long, Fei and Murshudov, Garib N.";
    bib_journal = "Acta Crystallographica Section D";

  }
  if (reference_id == COOT_REFERENCE_SSM) {

    ref_description = "If using \"SSM Superposition\", please cite:\n\n";

    title       = "Secondary-structure matching (SSM), a new tool for fast protein structure alignment in three dimensions";
    author      = "Krissinel E, Henrick K";
    journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year        = "2004";
    volume      = "60";
    number      = "12 Part 1";
    pages       = "2256-2268";

    bib_type  = "Article";
    bib_id    = "Krissinel:ba5056";
    bib_author = "Krissinel, E. and Henrick, K.";
    bib_journal = "Acta Crystallographica Section D";
  }

  if (reference_id == COOT_REFERENCE_MMDB) {

    ref_description = "The reference for the Macromolecular Database (MMDB) is:\n\n";

    title       = "The new CCP4 Coordinate Library as a toolkit for the design of coordinate-related applications in protein crystallography";
    author      = "Krissinel EB, Winn MD, Ballard CC, Ashton AW, Patel P, Potterton EA, McNicholas SJ, Cowtan KD, Emsley P.";
    journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year        = "2004";
    volume      = "60";
    pages       = "2250-2255";

    bib_type    = "Article";
    bib_id      = "Krissinel:ba5055";
    bib_author  = "Krissinel, E. B. and Winn, M. D. and Ballard, C. C. and Ashton, A. W. and Patel, P. and Potterton, E. A. and McNicholas, S. J. and Cowtan, K. D. and Emsley, P.";
    bib_journal = "Acta Crystallographica Section D";
  }
  if (reference_id == COOT_REFERENCE_CLIPPER) {

    ref_description = "The reference for clipper is:\n\n";

    title       = "The Clipper C++ libraries for X-ray crystallography";
    author      = "Cowtan K";
    journal     = "IUCr Computing Commission Newsletter";
    year        = "2003";
    volume      = "2";
    pages       = "4-9";

    bib_type   = "Article";
    bib_id     = "cowtan03:clipper";
    bib_author = "Kevin Cowtan";
  }

  if (reference_id == COOT_REFERENCE_BUCCANEER) {

    ref_description = "If using the Coot sequencing tool (\"Cootaneer\") or Fast Secondary Structure Search, please cite:\n\n";

    title       = "Fitting molecular fragments into electron density";
    author      = "Cowtan K";
    journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year        = "2008";
    volume      = "64";
    number      = "1";
    pages       = "83-89";

    bib_type    = "Article";
    bib_id      = "Cowtan:ba5104";
    bib_author  = "Cowtan, Kevin";
    bib_journal = "Acta Crystallographica Section D";
  }

  if (reference_id == COOT_REFERENCE_MOLPROBITY) {

    ref_description = "If using Molprobity tools (probe, reduce, etc.), please cite:\n\n";

    title       = "MolProbity: all-atom contacts and structure validation for proteins and nucleic acids";
    author      = "Davis IW, Leaver-Fay A, Chen VB, Block JN, Kapral GJ, Wang X, Murray LW, Arendall WB 3rd, Snoeyink J, Richardson JS, Richardson DC.";
    journal     = "Nucleic Acids Research";
    year        = "2007";
    volume      = "35";
    pages       = "375-383";

    bib_type    = "Article";
    bib_id      = "Davis:molprobity07";
    bib_author  = "Ian W. Davis and Andrew Leaver-Fay and Vincent B. Chen and Jeremy N. Block and Gary J. Kapral and Xueyi Wang and Laura W. Murray and W. Bryan Arendall III and Jack Snoeyink and Jane S. Richardson and David C. Richardson";
  }

  if (reference_id == COOT_REFERENCE_CALPHA) {

    ref_description = "The reference for CALPHA (idea as used in Coot when converting batons to main-chain atoms) is:\n\n";

    title       = "Polyalanine Reconstruction from C[alpha] Positions Using the Program CALPHA Can Aid Initial Phasing of Data by Molecular Replacement Procedures";
    author      = "Esnouf RM";
    journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year        = "1997";
    volume      = "53";
    number      = "6";
    pages       = "665--672";

    bib_type    = "Article";
    bib_id      = "Esnouf:ad0021";
    bib_title   = "Polyalanine Reconstruction from C{$\\alpha$} Positions Using the Program {\\it CALPHA} Can Aid Initial Phasing of Data by Molecular Replacement Procedures";
    bib_author  = "Esnouf, R. M.";
    bib_journal = "Acta Crystallographica Section D";
  }

  if (reference_id == COOT_REFERENCE_XLIGAND) {

    ref_description = "The reference for X-LIGAND (algorithm similar to is used in ligand fitting in Coot):\n\n";

    title       = "X-LIGAND: an application for the automated addition of flexible ligands into electron density";
    author = "Oldfield, TJ";
    journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year        = "2001";
    volume      = "57";
    pages       = "696-705";

    bib_type   = "Article";
    bib_id     = "Oldfield:be0006";
    bib_title  = "{\\it X-LIGAND}: an application for the automated addition of flexible ligands into electron density";
    bib_author = "Oldfield, T. J.";
    bib_journal = "Acta Crystallographica Section D";
  }

  if (reference_id == COOT_REFERENCE_EDS) {

    ref_description = "The reference for The Electron Density Server:\n\n";

    title       = "The Uppsala Electron-Density Server";
    author      = "Kleywegt GJ, Harris MR, Zou JY, Taylor TC, Wahlby A, Jones TA";
    journal     = "ACTA CRYSTALLOGRAPHICA SECTION D-BIOLOGICAL CRYSTALLOGRAPHY";
    year        = "2004";
    volume      = "60";
    pages       = "2240-2249";

    bib_type   = "Article";
    bib_id     = "kleywegt:eds";
    bib_title  = "The Uppsala Electron-Density Server";
    bib_author = "G. J.,Kleywegt and M. R. Harris and J. Y. Zou and T. C. Taylor and A. Wahlby A and T. A. Jones";
    bib_journal = "Acta Crystallographica Section D";
  }

  if (reference_id == COOT_REFERENCE_OTHERS) {
    ref_text = "others ref";
    bib_text = "otheres bib";
  }

  // assemble ref_text
  ref_text  = author + "\n";
  ref_text += title + "\n";
  ref_text += journal + " " + volume + ", " + pages + ", " + year + ".\n";

  // assemble bib_text
  bib_text  = "@" + bib_type + "{" + bib_id + ",\n";
  bib_text += "  author  =  {" + bib_author + "},\n";
  if (bib_title == "dummy") {
    bib_text += "  title   =  {" + title + "},\n";
  } else {
    bib_text += "  title   =  {" + bib_title + "},\n";
  }
  if (bib_journal == "dummy") {
    bib_text += "  journal =  {" + journal + "},\n";
  } else {
    bib_text += "  journal =  {" + bib_journal + "},\n";
  }
  bib_text += "  year    =  " + year + ",\n";
  bib_text += "  volume  =  " + volume + ",\n";
  if (pages == "") {
    bib_text += "\n";
  } else {
    bib_text += "  pages   =  " + pages + "\n";
  }
  bib_text += "}";

  gtk_text_buffer_get_end_iter(ref_buffer, &end_iter);
  if (ref_description != "dummy") {
    gtk_text_buffer_insert(ref_buffer, &end_iter, ref_description.c_str(), -1);
  }
  gtk_text_buffer_insert(ref_buffer, &end_iter, ref_text.c_str(), -1);
  gtk_text_buffer_set_text(bib_buffer, bib_text.c_str(), -1);

  gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), 0);

}

void set_graphics_window_size(int x_size, int y_size) {
   set_graphics_window_size_internal(x_size, y_size, 0);
}

void set_graphics_window_size_internal(int x_size, int y_size, int as_widget_flag) {

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      g.graphics_x_size = x_size;
      g.graphics_y_size = y_size;
      GtkWidget *win = g.get_main_window();
      if (win) {
	 GtkWindow *window = GTK_WINDOW(win);

#if (GTK_MAJOR_VERSION >= 4)
	 if (!as_widget_flag) {
	    // this is the default, gtk_window should be resized this way
            gtk_window_set_default_size(GTK_WINDOW(window), x_size, y_size);
	 } else {
	    // resize using this function at your own risk. It may do and/or result in unexpected things
            gtk_widget_set_size_request(GTK_WIDGET(window), x_size, y_size);
	 }
#else
         gtk_window_resize(window, x_size, y_size);
	 while (gtk_events_pending())
	    gtk_main_iteration();
	 while (gdk_events_pending()) // not a thing in GTK4
	    gtk_main_iteration();
#endif

      }
      graphics_draw();
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("set-graphics-window-size");
   command_strings.push_back(graphics_info_t::int_to_string(x_size));
   command_strings.push_back(graphics_info_t::int_to_string(y_size));
   add_to_history(command_strings);
}


void add_recentre_on_read_pdb_combobox(GtkWidget *filechooser) {

   GtkWidget *combobox = widget_from_builder("coords_filechooserdialog_recentre_combobox");

   // gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), "Recentre on Molecule");
   // gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), "Don't Recentre");
   // gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), "Recentre Molecule Here");

   if (graphics_info_t::recentre_on_read_pdb)
      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 0);
   if (!graphics_info_t::recentre_on_read_pdb)
      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 1);
}


void store_graphics_window_position(int x_pos, int y_pos) {

   graphics_info_t g;

   if ((x_pos == g.graphics_x_position) && (y_pos == g.graphics_y_position)) {
      // do nothing
   } else {
      g.graphics_x_position = x_pos;
      g.graphics_y_position = y_pos;

      std::string cmd = "store-graphics-window-position";
      std::vector<coot::command_arg_t> args;
      args.push_back(x_pos);
      args.push_back(y_pos);
      // add_to_history_typed(cmd, args); enough!
   }
}

void set_graphics_window_position(int x_pos, int y_pos) {

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      GtkWidget *main = g.get_main_window();
      if (main) {

#if (GTK_MAJOR_VERSION >= 4)

         std::cout << "20220602-PE FIXME in set_graphics_window_position() " << std::endl;
#else
         gtk_window_move(GTK_WINDOW(main), x_pos, y_pos); // window manager can ignore this
	 while (gtk_events_pending())
	    gtk_main_iteration();
#endif

      }
      graphics_draw();
   }
   std::string cmd = "set-graphics-window-position";
   std::vector<coot::command_arg_t> args;
   args.push_back(x_pos);
   args.push_back(y_pos);
   add_to_history_typed(cmd, args);
}



/*! \brief store the graphics window position and size to zenops-graphics-window-size-and-postion.scm in
 *         the preferences directory. */
void graphics_window_size_and_position_to_preferences() {

   // Note to self: is there a "get preferences dir" function?
   std::string h = coot::get_home_dir();
   if (!h.empty()) {
      // 20220507-PE pref_dir is now .coot
      std::string pref_dir = coot::util::append_dir_dir(h, ".coot");
      if (! coot::is_directory_p(pref_dir)) {
         // make it
	 // pref_dir = coot::get_directory(pref_dir); // oops not in this branch.
	 struct stat s;
	 int fstat = stat(pref_dir.c_str(), &s);
	 if (fstat == -1 ) { // file not exist
	    int status = coot::util::create_directory(pref_dir);
            if (status != 0) {
               std::cout << "status " << status << std::endl;
               std::string m("WARNING:: failed to create directory ");
               m += pref_dir;
               info_dialog(m.c_str()); // 20220507-PE make this argument a string one rainy day
            }
	 }
      }
      if (coot::is_directory_p(pref_dir)) {
         graphics_info_t g;
         int x  = g.graphics_x_position;
         int y  = g.graphics_y_position;
         int xs = g.graphics_x_size;
         int ys = g.graphics_y_size;

         GtkWidget *main_window = g.get_main_window();
         if (main_window) {
#if (GTK_MAJOR_VERSION >= 4)
            std::cout << "in graphics_window_size_and_position_to_preferences() find the window position and size" << std::endl;
#else
            gtk_window_get_position(GTK_WINDOW(main_window), &x, &y);
            gtk_window_get_size(GTK_WINDOW(main_window), &xs, &ys);
#endif

            std::string file_name = coot::util::append_dir_file(pref_dir, "xenops-graphics.scm");
            std::ofstream f(file_name.c_str());
            if (f) {
               f << "(set-graphics-window-position " << x  << " " << y  << ")\n";
               f << "(set-graphics-window-size     " << xs << " " << ys << ")\n";
            }
            f.close();
            file_name = coot::util::append_dir_file(pref_dir, "xenops-graphics.py");
            std::ofstream fp(file_name.c_str());
            if (fp) {
               fp << "import coot\n";
               fp << "coot.set_graphics_window_position(" << x  << ", " << y << ")\n";
               fp << "coot.set_graphics_window_size(" << xs << ", " << ys << ")\n";
            }
            fp.close();
         }
      } else {
         std::cout << "WARNING:: $HOME/.coot is not a directory - settings not saved" << std::endl;
         info_dialog("WARNING:: $HOME/.coot is not a directory - settings not saved");
      }
   }

}


/* a general purpose version of the above, where we pass a widget flag */
void
store_window_size(int window_type, GtkWidget *widget) {

   if (window_type == COOT_FILESELECTION_DIALOG) { // 20220319-PE bleugh
      GtkAllocation allocation;
      gtk_widget_get_allocation(widget, &allocation);
      graphics_info_t::file_chooser_dialog_x_size = allocation.width;
      graphics_info_t::file_chooser_dialog_y_size = allocation.height;
   }
}

void set_file_selection_dialog_size(GtkWidget *dialog) {

   if (graphics_info_t::file_chooser_dialog_x_size > 0) {
      graphics_info_t g;
      std::cout << "DEBUG:: set size request for dialog "
                << g.file_chooser_dialog_x_size << " " << g.file_chooser_dialog_y_size << std::endl;
      gtk_widget_set_size_request(dialog, g.file_chooser_dialog_x_size, g.file_chooser_dialog_y_size);
   }
}


/* return negative if fail */
float get_positive_float_from_entry(GtkEntry *w) {

   float f = -1.0;
   if (graphics_info_t::use_graphics_interface_flag) {
      const gchar *text = gtk_editable_get_text(GTK_EDITABLE(w));
      if (strlen(text) > 0) {
	 float tmp = atof(text);
	 if (tmp > 0) {
	    if (tmp < 9e10) {
	       f = tmp;
	    }
	 }
      }
   }
   return f;
}

// return TRUE if we don't want the window destroyed.
int
coot_checked_exit(int retval) {

   graphics_info_t g;

   // 20200822-PE save the (new) python history
   g.command_history.write_history();

   int i_unsaved = g.check_for_unsaved_changes();
   std::string cmd = "coot-checked-exit";
   std::vector<coot::command_arg_t> args;
   args.push_back(retval);
   add_to_history_typed(cmd, args);
   if (i_unsaved == 0) { // no unsaved.
      coot_real_exit(retval);
   }
   return TRUE; // path where there were unsaved changes, we don't
		// want to exit.
}

// return TRUE if we don't want the window destroyed.
gboolean
coot_checked_exit_gtk2(int retval) {

   graphics_info_t g;

   // 20200822-PE save the (new) python history
   g.command_history.write_history();

   int i_unsaved = g.check_for_unsaved_changes();
   std::string cmd = "coot-checked-exit";
   std::vector<coot::command_arg_t> args;
   args.push_back(retval);
   add_to_history_typed(cmd, args);
   if (i_unsaved == 0) { // no unsaved.
#ifdef USE_GUILE
#  ifdef USE_GUILE_GTK
      run_clear_backups(retval);
#  else
#    ifdef USE_PYGTK // MacOSX fink-build path
        run_clear_backups_py(retval);
#    endif
#  endif
#else
#  ifdef USE_PYGTK
      run_clear_backups_py(retval);
#  else
      coot_real_exit(retval);
#  endif // USE_PYGTK
#endif // USE_GUILE
   }
   return TRUE; // path where there were unsaved changes, we don't
		// want to exit.
}


#ifdef USE_GUILE
void run_clear_backups(int retval) {

   // just exit if we don't have guile-gtk

#ifndef USE_GUILE_GTK
   coot_real_exit(retval);
#else
   SCM r = safe_scheme_command("(clear-backups-maybe)");

   if (scm_is_undefined(r)) {
      // not false and not not false, function didn't run then...
      std::cout << "WARNING:: (clear-backups-maybe) returns "
		<< scm_to_locale_string(display_scm(r))
		<< std::endl;
      coot_real_exit(retval);
   } else {
//       std::cout << "DEBUG:: (clear-babckups-maybe) returned: "
// 		<< scm_to_locale_string(display_scm(r))
// 		<< std::endl;
   }

   // if r was #f then exit.
   //
   if (! SCM_NFALSEP(r)) { // backup gui was not needed/shown
      coot_real_exit(retval);
   }
#endif
}
#endif


#ifdef USE_PYTHON
void run_clear_backups_py(int retval) {

   PyObject *r = safe_python_command_with_return("clear_backups_maybe()");

   if (r == NULL || r == Py_None) {
      // not false and not not false, function didn't run then...
      std::cout << "WARNING:: clear_backups_maybe() returns "
		<< PyUnicode_AsUTF8String(PyObject_Str(r))
		<< std::endl;
      coot_real_exit(retval);
   } else {
//       std::cout << "DEBUG:: (clear-babckups-maybe) returned: "
// 		<< scm_to_locale_string(display_scm(r))
// 		<< std::endl;
   }

   // if r was #f then exit.
   //
   if (r == Py_False) { // backup gui was not needed/shown
      coot_real_exit(retval);
   }
}
#endif // PYTHON

void coot_clear_backup_or_real_exit(int retval) {

#ifdef USE_GUILE
   run_clear_backups(retval);
#else
#ifdef USE_PYGTK
   run_clear_backups_py(retval);
#else
   coot_real_exit(retval);
#endif // USE_PYGTK
#endif // USE_GUILE
}

void
coot_real_exit(int retval) {
   coot_save_state_and_exit(retval, 1);
}

void
coot_no_state_real_exit(int retval) {

   graphics_info_t g;
   g.command_history.write_history();

   // this is called (only) from on_window1_delete_event()
   coot_save_state_and_exit(retval, 0);
}

void
coot_save_state_and_exit(int retval, int save_state_flag) {

   graphics_info_t::static_thread_pool.stop(true);

   // wait for refinement to finish (c.f conditionally_wait_for_refinement_to_finish())
   while (graphics_info_t::restraints_lock) {
      std::this_thread::sleep_for(std::chrono::milliseconds(30));
   }

   if (save_state_flag) {
      save_state(); // we get error message in save_state()
   }

   // save the history
   if (save_state_flag) {
      graphics_info_t g;
      if (! g.disable_state_script_writing)
         g.save_history();
   }

#ifdef USE_MYSQL_DATABASE
   db_finish_up();
#endif // USE_MYSQL_DATABASE

   // #ifdef USE_PYTHON
   // Py_Finalize();
   // #endif

   for (int imol=0; imol<graphics_n_molecules(); imol++)
      graphics_info_t::molecules[imol].close_yourself();

   // why is this windows-only?

#ifdef WINDOWS_MINGW
   clipper::ClipperInstantiator::instance().destroy();
   graphics_info_t::static_thread_pool.stop(true);
#endif

   exit(retval);
}

// This is called by when the save coordinates menu item is pressed,
// for the save coordinates fileselection (other fileselections do not
// have the name on the action area vbox.
//
void
add_file_dialog_action_area_vbox(GtkWidget *fileselection) {

   std::cout << "GTK-FIXME no fileselection" << std::endl;
}


// where data type:
// 0 coords
// 1 mtz etc
// 2 maps
// 3 cif dictionary
// 4 scripting files
//
GtkWidget *add_filename_filter_button(GtkWidget *fileselection,
				      short int data_type) {

   GtkWidget *button = 0;
   add_filechooser_filter_button(fileselection, data_type);
   return button;
}

void add_save_coordinates_include_hydrogens_and_aniso_checkbutton(GtkWidget *fileselection) {

   bool no_chooser_filter = 1;

   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      	no_chooser_filter = 0;
   }

}


// Paul requested a new function for filechooser filter
// only available in gtk2
// here we go


// where data type:
// 0 coords
// 1 mtz etc
// 2 maps
// 3 cif dictionary
// 4 scripting files
//
void add_filechooser_filter_button(GtkWidget *filechooser, short int data_type) {

  int d = data_type;

  std::vector<std::string> globs;

  GtkFileFilter *filterall    = gtk_file_filter_new();
  GtkFileFilter *filterselect = gtk_file_filter_new();

  gtk_file_filter_set_name(filterall, "All Files");
  gtk_file_filter_add_pattern(filterall, "*");

  // actually only Choosers - I should change the define.

  if (d == COOT_COORDS_FILE_SELECTION || d == COOT_SAVE_COORDS_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "CoordinatesFiles");

    globs = *graphics_info_t::coordinates_glob_extensions;
  };

  if (d == COOT_DATASET_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "Data Files");

    globs = *graphics_info_t::data_glob_extensions;
  };

  if (d == COOT_MAP_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "Map Files");

    globs = *graphics_info_t::map_glob_extensions;
  };

  if (d == COOT_CIF_DICTIONARY_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "Dictionary Files");

    globs = *graphics_info_t::dictionary_glob_extensions;
  };

  if (d == COOT_SCRIPTS_FILE_SELECTION) {
    // BL says:: we dont have a script extensions (yet)
    // so we make one just here (no adding of extensions etc as yet)

    std::vector<std::string> script_glob_extension;
#ifdef USE_PYTHON
    script_glob_extension.push_back("*.py");
#endif // USE_PYTHON
#ifdef USE_GUILE
    script_glob_extension.push_back("*.scm");
#endif // USE_GUILE

    gtk_file_filter_set_name(filterselect, "scripting-files");
    g_object_set_data(G_OBJECT(filechooser), "filter", filterselect);

    globs = script_glob_extension;

  };

  std::string s;
  for (unsigned int i=0; i<globs.size(); i++) {
    s = "*";
    s += globs[i];
    gtk_file_filter_add_pattern (filterselect, s.c_str());
  };

  gtk_file_chooser_add_filter(GTK_FILE_CHOOSER (filechooser),
                              GTK_FILE_FILTER (filterall));
  gtk_file_chooser_add_filter(GTK_FILE_CHOOSER (filechooser),
                              GTK_FILE_FILTER (filterselect));

  if (filter_fileselection_filenames_state() == 1) {
    // filter automatically
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER (filechooser),
				GTK_FILE_FILTER (filterselect));
  } else {
    // show all
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER (filechooser),
				GTK_FILE_FILTER (filterall));
  }

}

// and lets have a function for extra custom filters//

void add_filechooser_extra_filter_button(GtkWidget *fileselection,
				      const gchar *filtername,
                                      const gchar *globname) {

  GtkFileFilter *filter = gtk_file_filter_new ();

  gtk_file_filter_set_name (filter, filtername);
  gtk_file_filter_add_pattern (filter, globname);
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER (fileselection),
			       GTK_FILE_FILTER (filter));

}


void
on_read_map_difference_map_toggle_button_toggled (GtkButton       *button,
						  gpointer         user_data)
{
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
      std::cout << "is a difference map...!\n";
   }
}


#include <gdk/gdkkeysyms.h> // for keyboarding.

#ifdef FIX_THE_KEY_PRESS_EVENTS
gboolean
on_filename_filter_key_press_event (GtkWidget       *widget,
				    GdkEventKey     *event,
				    gpointer         user_data) {

   //    if (event->keyval == GDK_Return || event->keyval == GDK_Tab)
   //    { // Tab is not good.  It takes you to the next widget too,
   //    which is not want we want.  It's confusing, so let's just use
   //    return.

   if (event->keyval == GDK_KEY_Return) {
      handle_filename_filter_gtk2(widget);
   }
   return FALSE;
}
#endif



std::string pre_directory_file_selection(GtkWidget *sort_button) {


   std::cout << "GTK-FIXME pre_directory_file_selection() " << std::endl;

   std::string pre_directory("");

   /*
   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(gtk_object_get_user_data(GTK_OBJECT(sort_button)));

   // The menu item is a container than contains a label.
   // How do we get to the label given the container?
   // Strangely enough we use the history_pulldown.

   GList *dlist = gtk_container_children(GTK_CONTAINER(history_pulldown));
   GList *free_list = dlist;

   while (dlist) {
      gchar *t = GTK_LABEL(dlist->data)->label;
      if (t != NULL) {
	 pre_directory = t;
      } else {
	 std::cout << "WARNING:: null label t " << std::endl;
      }
      dlist = dlist->next;
   }
   g_list_free(free_list);
   */

   return pre_directory;
}




void push_the_buttons_on_fileselection(GtkWidget *filter_button,
				       GtkWidget *sort_button,
				       GtkWidget *fileselection) {

   std::cout << "GTK-FIXME no fileselection C push the buttons" << std::endl;
}

void
filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v) {

   std::cout << "GTK-FIXME no fileselection filelist_into_fileselection_clist " << std::endl;
}

/*  Eleanor likes to sort her files by date when selecting a file
*/
GtkWidget *add_sort_button_fileselection(GtkWidget *fileselection) {

   std::cout << "GTK-FIXME no fileselection add_sort_button_fileselection" << std::endl;
   return 0;
}

void add_is_difference_map_checkbutton(GtkWidget *fileselection) {

  bool add_map_button = 1;

  if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      add_map_button = 0;
  }

}


void
on_recentre_on_read_pdb_toggle_button_toggled (GtkButton       *button,
					       gpointer         user_data)
{
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
      std::cout << "INFO:: activated recentering on new coordinates.\n";
   } else {
      std::cout << "INFO:: de-activated recentering on new coordinates.\n";
   }
}


/*  ------------------------------------------------------------------------ */
/*              scripting gtk interface                                      */
/*  ------------------------------------------------------------------------ */

#ifdef FIX_THE_KEY_PRESS_EVENTS
#endif



void reveal_python_scripting_entry() {
   GtkRevealer* revealer = GTK_REVEALER(widget_from_builder("python_scripting_revealer"));
   gtk_revealer_set_reveal_child(revealer,TRUE);
}

void toggle_reveal_python_scripting_entry() {
   GtkRevealer* revealer = GTK_REVEALER(widget_from_builder("python_scripting_revealer"));
   gtk_revealer_set_reveal_child(revealer,!gtk_revealer_get_reveal_child(revealer));
}

// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_guile_window_entry(GtkWidget *entry) {

#ifdef USE_GUILE
#endif //  USE_GUILE

}

// Similar to fill_option_menu_with_coordinates_options, but I moved
// it to graphics_info_t because it is also used when there is an
// ambiguity in the map for refinement (graphics_info_t::refine)
//
int fill_combobox_with_map_options(GtkWidget *combobox, GCallback signalfunc) {

   graphics_info_t g;
   int imol_active = -1;
   int ii = imol_refinement_map();
   if (is_valid_map_molecule(ii)) {
      imol_active = ii;
   } else {
      for (int i=0; i<g.n_molecules(); i++) {
	 if (is_valid_map_molecule(i)) {
	    imol_active = i;
	    break;
	 }
      }
   }
   g.fill_combobox_with_map_options(combobox, signalfunc, imol_active);
   return imol_active;

}


// This is for maps which come from mtz (i.e. have SFs)
int fill_combobox_with_map_mtz_options(GtkWidget *combobox, GCallback signalfunc) {
   graphics_info_t g;
   return g.fill_combobox_with_map_mtz_options(combobox, signalfunc, 0);
}


void set_on_off_single_map_skeleton_radio_buttons(GtkWidget *skeleton_frame,
						  int imol) {
   graphics_info_t g;
   g.set_on_off_single_map_skeleton_radio_buttons(skeleton_frame, imol);

}

void set_contour_sigma_button_and_entry(GtkWidget *window, int imol) {
   graphics_info_t g;
   g.set_contour_sigma_button_and_entry(window, imol);
}

// and the reverse function (button -> value)
void set_contour_by_sigma_step_maybe(GtkWidget *window, int imol) {

   // GtkWidget *button = lookup_widget(window, "single_map_sigma_checkbutton");
   // GtkWidget *entry  = lookup_widget(window, "single_map_sigma_step_entry");
   GtkWidget *button = widget_from_builder("single_map_sigma_checkbutton");
   GtkWidget *entry  = widget_from_builder("single_map_sigma_step_entry");

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
      const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
      if (text) {
	 float v = atof(text);
// 	 graphics_info_t::molecules[imol].contour_by_sigma_flag = 1;
// 	 graphics_info_t::molecules[imol].contour_sigma_step = v;
	 graphics_info_t::molecules[imol].set_contour_by_sigma_step(v, 1);
      }
   } else {
      // 0.0 is ignored.
      graphics_info_t::molecules[imol].set_contour_by_sigma_step(0.0, 0);
   }
}

/*  ------------------------------------------------------------------------ */
/*              widget utilities                                             */
/*  ------------------------------------------------------------------------ */
void set_transient_and_position(int widget_type, GtkWidget *window) {

   GtkWidget *main_window_widget = graphics_info_t::get_main_window();
   if (main_window_widget) {
      GtkWindow *main_window = GTK_WINDOW(main_window_widget);
      gtk_window_set_transient_for(GTK_WINDOW(window), main_window);
      if (widget_type == COOT_DELETE_WINDOW) {

         bool done_set_pos = false;
         if (graphics_info_t::delete_item_widget_x_position > -100) {
            if (graphics_info_t::delete_item_widget_y_position > -100) {

#if (GTK_MAJOR_VERSION >= 4)
               std::cout << "in set_transient_and_position() FIXME gtk_window_move()" << std::endl;
#else
               gtk_window_move(GTK_WINDOW(window),
                                 graphics_info_t::delete_item_widget_x_position,
                                 graphics_info_t::delete_item_widget_y_position);
#endif

                     // 	       gtk_widget_set_uposition(window,
                     // 					graphics_info_t::delete_item_widget_x_position,
                     // 					graphics_info_t::delete_item_widget_y_position);
               done_set_pos = true;
            }
         }
         if (! done_set_pos) {
            int x_pos = graphics_info_t::graphics_x_position - 100;
            int y_pos = graphics_info_t::graphics_y_position + 100;
            if (x_pos < 5) x_pos = 5;
               std::cout << "GTK-FIXME no gtk_widget_set_uposition D" << std::endl;
            // gtk_widget_set_uposition(window, x_pos, y_pos);
         }
      }
   }
}

void set_transient_for_main_window(GtkWidget *dialog) {

   graphics_info_t::set_transient_for_main_window(dialog);
}


GtkWidget *coot_file_chooser() {

   GtkWidget *w = widget_from_builder("coords_filechooser_dialog");
   // gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(w), TRUE);

   return w;
}

GtkWidget *coot_dataset_chooser() {

   // GtkWidget *w = create_dataset_filechooserdialog1();
   GtkWidget *w = widget_from_builder("dataset_filechooser_dialog");
   return w;
}

GtkWidget *coot_map_name_chooser() {

   GtkWidget *w = widget_from_builder("map_name_filechooser_dialog");
   return w;
}

GtkWidget *coot_save_coords_chooser() {

   GtkWidget *w = widget_from_builder("save_coordinates_filechooser_dialog");

   // gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (w), TRUE);

   return w;
}

GtkWidget *coot_cif_dictionary_chooser() {

   GtkWidget *w = widget_from_builder("cif_dictionary_filechooser_dialog");
   return w;
}

GtkWidget *coot_run_script_chooser() {

   GtkWidget *w = widget_from_builder("run_script_filechooser_dialog");
   return w;
}

GtkWidget *coot_save_state_chooser() {

   GtkWidget *w = widget_from_builder("save_state_filechooserdialog");
   // gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (w), TRUE);
   return w;
}

GtkWidget *coot_save_symmetry_chooser() {

   // 20220625-PE
   // dialog widgets need to be able to default so that gtk_widget_grab_default(), which
   // is called by gtk_file_chooser_set_do_overwrite_confirmation(), can work
   // without a critcal (console) error.
   //
   // Does that need to be set in the glade file? I did so and it seemed to work
   // (whereas setting it here did not).

   GtkWidget *w = widget_from_builder("save_symmetry_coords_filechooser_dialog");
   // gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER (w), TRUE);
   return w;
}

GtkWidget *coot_screendump_chooser() {

   GtkWidget *w = widget_from_builder("screendump_filechooser_dialog");
   // gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER (w), TRUE);
   return w;
}


void set_directory_for_coot_file_chooser(GtkWidget *coords_chooser) {

   set_directory_for_filechooser(coords_chooser);
}

const char *coot_file_chooser_file_name(GtkWidget *widget) {

   const char *f = 0;
   return f;
}

// this is not a gui function - it is a network function
//
void network_get_accession_code_entity(const std::string &text, int mode) {

   auto join = [] (const std::string &d, const std::string &f) {
      return d + std::string("/") + f;
   };

   auto join_2d = [] (const std::string &d1, const std::string &d2, const std::string &f) {
      return d1 + std::string("/") + d2 + std::string("/") + f;
   };

   // 20240630-PE need to check that the file already exists before downloading it
   xdg_t xdg;
   std::string download_dir = xdg.get_download_dir();
   make_directory_maybe(download_dir);
   std::string dld = coot::get_directory(download_dir);
   if (! dld.empty()) {
      download_dir = dld;
      std::string down_id = coot::util::downcase(text);
      std::string pdbe_server = "https://www.ebi.ac.uk";
      std::string pdbe_pdb_file_dir = "pdbe/entry-files/download";
      std::string pdb_url_dir = pdbe_server + "/" + pdbe_pdb_file_dir;

      std::string pdb_file_name = std::string("pdb") + down_id + std::string(".ent");
      std::string cif_file_name =                      down_id + std::string(".cif");
      std::string pdb_filepath = coot::util::append_dir_file(download_dir, pdb_file_name);
      std::string cif_filepath = coot::util::append_dir_file(download_dir, cif_file_name);

      std::string pdb_url = join(pdb_url_dir, pdb_file_name);
      std::string cif_url = join(pdb_url_dir, cif_file_name);

      if (mode == 1) { // mtz mode
         std::string mtz_file_name = down_id + std::string("_map.mtz");
         std::string mtz_filepath = coot::util::append_dir_file(download_dir, mtz_file_name);
         if (coot::file_exists(mtz_filepath)) {
            auto_read_make_and_draw_maps(mtz_filepath.c_str()); // get rid of .c_str() one day
         } else {
            std::string mtz_url = join_2d(pdbe_server, pdbe_pdb_file_dir, mtz_file_name);
            int status = coot_get_url(mtz_url, mtz_filepath);
            if (status == 0) {
               auto_read_make_and_draw_maps(mtz_filepath.c_str());
            }
         }
      } else {
         // blocking!
         int status = coot_get_url(pdb_url, pdb_filepath);
         // coot_get_url() returns the return value of curl_easy_perform()
         // CURLE_OK is 0
         if (coot::file_exists(pdb_filepath)) {
            read_pdb(pdb_filepath);
         } else {
            if (status == 0) {
               read_pdb(pdb_filepath);
            } else {
               if (coot::file_exists(cif_filepath)) {
                  read_pdb(cif_filepath);
               } else {
                  status = coot_get_url(cif_url, cif_filepath);
                  if (status == 0) {
                     read_pdb(cif_filepath);
                  }
               }
            }
         }
      }
   }
}



/*  ----------------------------------------------------------------------- */
/*                  get by accession code:                                  */
/*  ----------------------------------------------------------------------- */

void handle_get_accession_code(GtkWidget *frame, GtkWidget *entry) {

   auto join = [] (const std::string &d, const std::string &f) {
      return d + std::string("/") + f;
   };

   auto fetch_pdb_redo = [join] (const std::string &code) {

      // 20240630-PE need to check that the file already exists before downloading it
      xdg_t xdg;
      std::string download_dir = join(xdg.get_cache_home().string(), "coot-download");
      std::string dld = coot::get_directory(download_dir);
      if (! dld.empty()) {
         download_dir = dld;
         std::string down_id = coot::util::downcase(code);
         std::string server = "https://pdb-redo.eu";
         std::string server_dir = std::string("db") + "/" + code;
         std::string pdb_file_name = code + "_final.pdb";
         std::string mtz_file_name = code + "_final.mtz";
         // make a "join()" function
         std::string pdb_url = server + "/" + server_dir + "/" + pdb_file_name;
         std::string mtz_url = server + "/" + server_dir + "/" + mtz_file_name;
         std::string pdb_filepath = coot::util::append_dir_dir(download_dir, pdb_file_name);
         std::string mtz_filepath = coot::util::append_dir_dir(download_dir, mtz_file_name);
         int status = coot_get_url(pdb_url, pdb_filepath);
         if (status == 0) {
            read_pdb(pdb_filepath);
            status = coot_get_url(mtz_url, mtz_filepath);
            if (status == 0) {
               // why is auto_read_mtz() not a thing? Use a std::string arg
               auto_read_make_and_draw_maps(mtz_filepath.c_str());
            }
         }
      } else {
         std::cout << "WARNING:: failed to make directory " << download_dir
                   << std::endl;
      }

   };

   auto network_get = [fetch_pdb_redo] (const std::string &text, int n) {

      if (n == COOT_ACCESSION_CODE_WINDOW_OCA) {
         network_get_accession_code_entity(text, 0); // coords
      }
      if (n == COOT_ACCESSION_CODE_WINDOW_EDS) {
         network_get_accession_code_entity(text, 0); // coords
         network_get_accession_code_entity(text, 1); // mtz
      }
      if (n == COOT_ACCESSION_CODE_WINDOW_OCA_WITH_SF) {
         std::cout << "WARNING:: OCA+SF no longer supported" << std::endl;
      }
      if (n == COOT_ACCESSION_CODE_WINDOW_PDB_REDO) {
         fetch_pdb_redo(text);
      }
      if (n == COOT_UNIPROT_ID) {
         fetch_alphafold_model_for_uniprot_id(text);
      }
   };

   const gchar *text_c = gtk_editable_get_text(GTK_EDITABLE(entry));

   if (! text_c) {
      std::cout << "WARNING:: handle_get_accession_code no text " << std::endl;
   } else {
      std::string text_s = std::string(text_c);
      std::string text = coot::util::remove_trailing_whitespace(text_s);
      std::cout << "PDB Accession Code: " << text << std::endl;
      std::cout << "frame: " << frame << std::endl;
      int n = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(frame), "mode"));
      std::cout << "DEBUG:: extracted accession code handle mode n " << n << std::endl;
      if (n == COOT_EMDB_CODE) {
         fetch_emdb_map(text);
      } else {
         if (n == COOT_COD_CODE) {
#ifdef USE_LIBCURL
            fetch_cod_entry(text);
#endif
         } else {
            network_get(text, n);
         }
      }
   }

   // and hide the accession code window
   gtk_widget_set_visible(frame, FALSE);
}



// Set the internal state of the torsion restraints
// (hang-over from old interface)
void do_torsions_toggle(GtkWidget *checkbutton) {

   graphics_info_t g;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton)))
      g.do_torsion_restraints = 1;
   else
      g.do_torsion_restraints = 0;
}

void set_refine_params_comboboxes() {

   graphics_info_t g;
   GtkWidget *cb1 = widget_from_builder("refine_params_geman_mcclure_alpha_combobox");
   GtkWidget *cb2 = widget_from_builder("refine_params_rama_restraints_combobox");
   GtkWidget *cb3 = widget_from_builder("refine_params_lennard_jones_epsilon_combobox");
   GtkWidget *cb4 = widget_from_builder("refine_params_torsion_weight_combobox");
   GtkWidget *cb5 = widget_from_builder("refine_params_overall_weight_combobox");
   GtkWidget *tb  = widget_from_builder("refine_params_more_control_togglebutton");

   if (cb1) gtk_combo_box_set_active(GTK_COMBO_BOX(cb1), g.refine_params_dialog_geman_mcclure_alpha_combobox_position);
   if (cb2) gtk_combo_box_set_active(GTK_COMBO_BOX(cb2), g.refine_params_dialog_rama_restraints_weight_combobox_position);
   if (cb3) gtk_combo_box_set_active(GTK_COMBO_BOX(cb3), g.refine_params_dialog_lennard_jones_epsilon_combobox_position);
   if (cb4) gtk_combo_box_set_active(GTK_COMBO_BOX(cb4), g.refine_params_dialog_torsions_weight_combox_position);


   if (tb) {
      if (g.refine_params_dialog_extra_control_frame_is_visible) {
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tb), TRUE);
         // GtkWidget *frame = lookup_widget(button, "refine_params_more_control_frame");
         // gtk_widget_set_visible(frame, TRUE);
      }
   }

}

// void fill_chiral_volume_molecule_option_menu(GtkWidget *w) {
// }

// 20211012-PE temporary arrangement
void
new_fill_combobox_with_coordinates_options(GtkWidget *combobox_molecule, GCallback callback_func, int imol_active);


void fill_chiral_volume_molecule_combobox(GtkWidget *dialog) {

   // GtkWidget *combobox = lookup_widget(dialog, "check_chiral_volumes_molecule_combobox");

   GtkWidget *combobox = widget_from_builder("check_chiral_volumes_molecule_combobox");

   graphics_info_t g;
   // int imol = graphics_info_t::check_chiral_volume_molecule;
   GCallback callback_func = G_CALLBACK(g.check_chiral_volume_molecule_combobox_changed);

   GtkWidget *vbox = widget_from_builder("check_chiral_volumes_dialog_vbox");

   auto my_delete_box_items = [] (GtkWidget *widget, void *data) {

#if (GTK_MAJOR_VERSION >= 4)
      std::cout << "FIXME my_delete_box_items() in fill_chiral_volume_molecule_combobox() " << std::endl;
#else
      gtk_container_remove(GTK_CONTAINER(data), widget);
#endif
   };

#if (GTK_MAJOR_VERSION >= 4)
      std::cout << "FIXME the-foreach in fill_chiral_volume_molecule_combobox() " << std::endl;
#else
      gtk_container_foreach(GTK_CONTAINER(vbox), my_delete_box_items, vbox);
#endif

   // 20211011-PE the code says not to use this function - I don't know why
   // g.fill_combobox_with_coordinates_options(combobox, callback_func, imol);
   //
   // Use fill_combobox_with_molecule_options() instead.

   std::vector<int> molecule_indices;
   std::vector<int> maps_vec;
   for (int i=0; i<g.n_molecules(); i++)
      if (is_valid_model_molecule(i))
         molecule_indices.push_back(i);

   if (! molecule_indices.empty()) {
      int imol_first = molecule_indices[0];
      // g.fill_combobox_with_molecule_options(combobox, callback_func, imol_first, molecule_indices);
      GtkWidget *combobox_new = gtk_combo_box_new();
      gtk_widget_set_visible(combobox_new, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      gtk_box_append(GTK_BOX(vbox), combobox_new);
#else
      gtk_box_pack_start(GTK_BOX(vbox), combobox_new, FALSE, FALSE, 4);
#endif
      g.new_fill_combobox_with_coordinates_options(combobox_new, callback_func, imol_first);
   }

}

void
pepflips_by_difference_map_dialog() {

   graphics_info_t g;
   GtkWidget *dialog = widget_from_builder("pepflips_by_difference_map_dialog");
   GtkWidget *model_combobox = widget_from_builder("pepflips_by_difference_map_dialog_model_comboboxtext");
   GtkWidget *map_combobox   = widget_from_builder("pepflips_by_difference_map_dialog_map_comboboxtext");

   std::cout << "debug model_combobox " << model_combobox << std::endl;
   std::cout << "debug   map_combobox " <<   map_combobox << std::endl;

   // clear combox boxes from that vbox:
   //
   gtk_cell_layout_clear(GTK_CELL_LAYOUT(model_combobox));
   gtk_cell_layout_clear(GTK_CELL_LAYOUT(map_combobox));

   GtkWidget *entry = widget_from_builder("pepflips_by_difference_map_dialog_entry");
   gtk_editable_set_text(GTK_EDITABLE(entry), "3.6");
   // create new comboboxes

   int imol_active = 0; // doesn't matter
   int imol_map = imol_refinement_map();
   GCallback callback = G_CALLBACK(NULL); // combox box is only read on Apply button press
   g.new_fill_combobox_with_coordinates_options(model_combobox, callback, imol_active);
   g.fill_combobox_with_difference_map_options(map_combobox, callback, imol_map);


   set_transient_for_main_window(dialog);
   gtk_widget_set_visible(dialog, TRUE);

   g_object_set_data(G_OBJECT(dialog), "model_combobox", model_combobox);
   g_object_set_data(G_OBJECT(dialog),   "map_combobox",   map_combobox);
}

#include "coot-utils/pepflip-using-difference-map.hh"

void
on_pepflip_residue_button_clicked(GtkButton *button, gpointer user_data) {

   coot::residue_spec_t *spec = reinterpret_cast<coot::residue_spec_t *>(user_data);
   graphics_info_t g;
   int imol = spec->int_user_data;
   g.go_to_residue(imol, *spec);

}


void pepflips_by_difference_map_results_dialog(int imol_coords, int imol_difference_map, float n_sigma) {

   typedef std::tuple<std::string, GCallback, gpointer> button_tuple;

   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_difference_map)) {
         graphics_info_t g;
         if (g.molecules[imol_difference_map].is_difference_map_p()) {
            const clipper::Xmap<float> &diff_xmap = g.molecules[imol_difference_map].xmap;
            mmdb::Manager *mol = g.molecules[imol_coords].atom_sel.mol;
            coot::pepflip_using_difference_map pf(mol, diff_xmap);
            std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);

            if (! flips.empty()) {
               std::vector<button_tuple> buttons;
               for (unsigned int i=0; i<flips.size(); i++) {
                  mmdb::Residue *residue_p = flips[i].get_residue(mol);
                  if (residue_p) {
                     GCallback cb = G_CALLBACK(on_pepflip_residue_button_clicked);
                     std::string res_name = residue_p->GetResName();
                     std::string button_label = flips[i].label(res_name);
                     coot::residue_spec_t *spec = new coot::residue_spec_t(flips[i]);
                     spec->int_user_data = imol_coords;
                     button_tuple bt = std::tuple<std::string, GCallback, gpointer>(button_label, cb, spec);
                     buttons.push_back(bt);
                  }
               }
               GtkWidget *dialog = g.dialog_box_of_buttons_internal("Pepflips", buttons, " Close ");
               gtk_widget_set_visible(dialog, TRUE);
            } else {
               info_dialog("No pepflips found");
            }
         }
      }
   }
}


/*  ------------------------------------------------------------------------ */
//            enviromnent and other distances
/*  ------------------------------------------------------------------------ */
//


void toggle_environment_show_distances(GtkCheckButton *button) {

   graphics_info_t g;

   GtkWidget *hbox                    = widget_from_builder("environment_distance_distances_frame");
   GtkWidget *distance_type_frame     = widget_from_builder("environment_distances_type_selection");
   GtkWidget *label_atom_check_button = widget_from_builder("environment_distance_label_atom_checkbutton");

   if (gtk_check_button_get_active(button)) {

      g.environment_show_distances = 1;
      gtk_widget_set_sensitive(hbox, TRUE);
      gtk_widget_set_sensitive(label_atom_check_button, TRUE);
      gtk_widget_set_sensitive(distance_type_frame, TRUE);

      std::pair<int, int> r = g.get_closest_atom();
      if (r.first >= 0) {
	 g.mol_no_for_environment_distances = r.second;
	 g.update_environment_distances_maybe(r.first, r.second);
	 graphics_draw();
      }

   } else {
      // std::cout << "toggled evironment distances off" << std::endl;
      g.environment_show_distances = 0;
      gtk_widget_set_sensitive(hbox, FALSE);
      gtk_widget_set_sensitive(distance_type_frame, FALSE);
      graphics_draw();
      // gtk_widget_set_sensitive(label_atom_check_button, FALSE); // keep it always active
   }
}

/* a graphics_info_t function wrapper: */
void residue_info_release_memory(GtkWidget *widget) {

   graphics_info_t g;
   g.residue_info_release_memory(widget);

}

void unset_residue_info_widget() {
   graphics_info_t g;
   g.residue_info_dialog = NULL;
}


/*  ----------------------------------------------------------------------- */
/*                  pointer distances                                      */
/*  ----------------------------------------------------------------------- */
void fill_pointer_distances_widget(GtkWidget *widget) {

   // GtkWidget *min_entry   = lookup_widget(widget, "pointer_distances_min_dist_entry");
   // GtkWidget *max_entry   = lookup_widget(widget, "pointer_distances_max_dist_entry");
   // GtkWidget *checkbutton = lookup_widget(widget, "pointer_distances_checkbutton");
   // GtkWidget *frame       = lookup_widget(widget, "pointer_distances_frame");
   GtkWidget *min_entry   = widget_from_builder("pointer_distances_min_dist_entry");
   GtkWidget *max_entry   = widget_from_builder("pointer_distances_max_dist_entry");
   GtkWidget *checkbutton = widget_from_builder("pointer_distances_checkbutton");
   GtkWidget *frame       = widget_from_builder("pointer_distances_frame");
   GtkWidget *grid        = widget_from_builder("show_pointer_distances_grid");

   float min_dist = graphics_info_t::pointer_min_dist;
   float max_dist = graphics_info_t::pointer_max_dist;

   std::cout << "here A with min_entry "   << min_entry << std::endl;
   std::cout << "here A with max_entry "   << max_entry << std::endl;
   std::cout << "here A with checkbutton " << checkbutton << std::endl;
   std::cout << "here A with frame "       << frame << std::endl;
   std::cout << "here A with grid "        << grid << std::endl;

   gtk_editable_set_text(GTK_EDITABLE(min_entry), graphics_info_t::float_to_string(min_dist).c_str());
   gtk_editable_set_text(GTK_EDITABLE(max_entry), graphics_info_t::float_to_string(max_dist).c_str());

   if (graphics_info_t::show_pointer_distances_flag) {
      gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);
      gtk_widget_set_sensitive(grid, TRUE);
   } else {
      gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), FALSE);
      gtk_widget_set_sensitive(grid, FALSE);
   }

}

void execute_pointer_distances_settings(GtkWidget *widget) {

   GtkWidget *min_entry   = widget_from_builder("pointer_distances_min_dist_entry");
   GtkWidget *max_entry   = widget_from_builder("pointer_distances_max_dist_entry");

   float min_dist = 0.0;
   float max_dist = 0.0;

   try {

      const gchar *tt = gtk_editable_get_text(GTK_EDITABLE(min_entry));
      float t = coot::util::string_to_float(std::string(tt));
      if ((t >= 0.0) && (t < 999.9))
         min_dist = t;

      tt = gtk_editable_get_text(GTK_EDITABLE(max_entry));
      t = coot::util::string_to_float(std::string(tt));
      if ((t >= 0.0) && (t < 999.9))
         max_dist = t;

      graphics_info_t::pointer_max_dist = max_dist;
      graphics_info_t::pointer_min_dist = min_dist;

      graphics_info_t g;
      g.make_pointer_distance_objects();
      g.graphics_draw();
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

}

void toggle_pointer_distances_show_distances(GtkCheckButton *checkbutton) {

   GtkWidget *grid = widget_from_builder("show_pointer_distances_grid");

   if (gtk_check_button_get_active(checkbutton)) {
      set_show_pointer_distances(1);
      gtk_widget_set_sensitive(grid, TRUE);
   } else {
      set_show_pointer_distances(0);
      gtk_widget_set_sensitive(grid, FALSE);
   }
}


/*  ------------------------------------------------------------------------ */
//            popup-menu for model_toolbar
/*  ------------------------------------------------------------------------ */
//

#ifdef FIX_THE_KEY_PRESS_EVENTS
void
toolbar_popup_menu (GtkToolbar *toolbar,
		    GdkEventButton *event_button,
		    gpointer user_data)
{
   // deprecated GTK-FIXME


   // I don't know how to fix this

   /*
   GtkHandleBox *hdlbox = GTK_HANDLE_BOX(GTK_WIDGET(toolbar)->parent);
   GtkWidget *menu = gtk_menu_new ();
   GtkWidget *item;

   static struct {
      char const *text;
      coot::model_toolbar::toolbar_position_type pos;
   } const pos_items[] = {
      { N_("Display to the right"),  coot::model_toolbar::RIGHT },
      { N_("Display to the left"),   coot::model_toolbar::LEFT },
      { N_("Display on the top"),    coot::model_toolbar::TOP },
      { N_("Display on the bottom"), coot::model_toolbar::BOTTOM },
   };

   if (hdlbox->child_detached) {
      item = gtk_menu_item_new_with_label (_("Reattach to main window"));
      gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);
      g_signal_connect (G_OBJECT (item), "activate",
			G_CALLBACK (reattach_modelling_toolbar),
			NULL);
   } else {
      size_t ui;
      GSList *group = NULL;

      for (ui = 0; ui < G_N_ELEMENTS (pos_items); ui++) {
	 char const *text = _(pos_items[ui].text);
	 coot::model_toolbar::toolbar_position_type pos = pos_items[ui].pos;

	 item = gtk_radio_menu_item_new_with_label(group, text);
	 group = gtk_radio_menu_item_get_group(GTK_RADIO_MENU_ITEM (item));

	 if (graphics_info_t::model_toolbar_position_state == pos) {
	    GTK_CHECK_MENU_ITEM(item)->active = 1;
	 } else {
	    GTK_CHECK_MENU_ITEM(item)->active = 0;
	 }

	 gtk_menu_shell_append(GTK_MENU_SHELL (menu), item);
	 g_object_set_data(G_OBJECT (item), "position", GINT_TO_POINTER (pos));
	 g_signal_connect(G_OBJECT (item), "toggled",
			  G_CALLBACK (set_model_toolbar_docked_position_callback),
			  item);
      }

      //    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(save_w), TRUE);

   }

   item = gtk_menu_item_new ();
   gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);
   gtk_widget_set_sensitive (item, FALSE);

   item = gtk_menu_item_new_with_label (_("Hide"));
   gtk_menu_shell_append (GTK_MENU_SHELL (menu), item);
   g_signal_connect (G_OBJECT (item), "activate",
		     G_CALLBACK (hide_modelling_toolbar),
		     NULL);

   gtk_widget_set_visible_all (menu, TRUE);
   gtk_menu_popup (GTK_MENU(menu), NULL, NULL, NULL, NULL, 0,
		   (event_button != NULL) ? event_button->time
		   : gtk_get_current_event_time());
   */
}
#endif


int suck_model_fit_dialog_bl() {

   /*

   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *main_window_hbox = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),  // commented
						  "main_window_hbox");
      GtkWidget *dialog = graphics_info_t::model_fit_refine_dialog;

      if (main_window_hbox) {
	 if (dialog) {
	    GtkWidget *handlebox2 = gtk_handle_box_new();
	    GtkWidget *hbox = lookup_widget(dialog, "model_fit_refine_dialog_vbox");  // commented
	    gtk_container_add (GTK_CONTAINER(handlebox2), hbox);
	    gtk_widget_reparent(hbox, handlebox2);
	    gtk_box_pack_start(GTK_BOX(main_window_hbox), handlebox2, FALSE, TRUE, 0);
	    gtk_widget_set_visible(handlebox2, TRUE);
	    gtk_widget_destroy(dialog);
	 } else {
	    std::cout << "no dialog\n";
	 }
      } else {
	 std::cout << "no hbox\n";
      }
      // I mean that it *should be* sucked.
      graphics_info_t::model_fit_refine_dialog_was_sucked = 1;
   }
   add_to_history_simple("suck-model-fit-dialog-bl");

   */
   return 0;
}

int suck_model_fit_dialog() {

   /*
   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *main_window_hbox = lookup_widget(GTK_WIDGET(graphics_info_t::glarea), // commented
						  "main_window_hbox");
      GtkWidget *main_window_side_frame = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),  // commented
							"main_window_model_fit_dialog_frame");
      GtkWidget *dialog = graphics_info_t::model_fit_refine_dialog;

      if (main_window_hbox) {
	 if (dialog) {
	    GtkWidget *hbox = lookup_widget(dialog, "model_fit_refine_dialog_vbox");  // commented
	    gtk_widget_reparent(hbox, main_window_side_frame);
	    gtk_widget_set_visible(main_window_side_frame, TRUE);
	    gtk_widget_destroy(dialog);
	 } else {
	    std::cout << "no dialog\n";
	 }
      } else {
	 std::cout << "no hbox\n";
      }
      // I mean that it *should be* sucked.
      graphics_info_t::model_fit_refine_dialog_was_sucked = 1;
   }
   add_to_history_simple("suck-model-fit-dialog");
   */
   return 0;
}


/* Return the dialog if it exists, else null */
// Actually, this doesn't do what it says on the tin.  It closes the
// dialog_hbox and returns the dialog.
//
GtkWidget *close_model_fit_dialog(GtkWidget *dialog_hbox) {

   GtkWidget *w = NULL;
   if (graphics_info_t::model_fit_refine_dialog_was_sucked) {
      GtkWidget *main_window_side_frame = widget_from_builder("main_window_model_fit_dialog_frame");
      gtk_widget_set_visible(dialog_hbox, FALSE);
      gtk_widget_set_visible(main_window_side_frame, FALSE);
   } else {
      w = widget_from_builder("model_refine_dialog");
   }
   // graphics_info_t::model_fit_refine_dialog_was_sucked = 0;

   return w;
}


/*  ------------------------------------------------------------------------ */
//            main_toolbar things
/*  ------------------------------------------------------------------------ */
//

/*! \brief hide the horizontal main toolbar in the GTK2 version */
void hide_main_toolbar() {
   if (graphics_info_t::use_graphics_interface_flag) {
      // GtkWidget *w = lookup_widget(graphics_info_t::get_main_window(), "main_toolbar");
      GtkWidget *w = widget_from_builder("main_toolbar");
      if (!w) {
	 std::cout << "hide_main_toolbar(): failed to lookup main toolbar" << std::endl;
      } else {
	 graphics_info_t::main_toolbar_show_hide_state = 0;
	 gtk_widget_set_visible(w, FALSE);
      }
   }
}

/*! \brief show the horizontal maub toolbar in the GTK2 version
  (the toolbar is shown by default) */
void show_main_toolbar() {

   // main_toolbar no longer exists - do I still want this function?
   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *w = widget_from_builder("main_toolbar");

      if (!w) {
	 std::cout << "show_main_toolbar(): failed to lookup main toolbar" << std::endl;
      } else {
	 graphics_info_t::main_toolbar_show_hide_state = 1;
	 gtk_widget_set_visible(w, TRUE);
      }
   }
}


/*  ------------------------------------------------------------------------ */
// other modelling tools
/*  ------------------------------------------------------------------------ */

GtkWidget *wrapped_create_other_model_tools_dialog() {

   GtkWidget *widget = graphics_info_t::other_modelling_tools_dialog;
   if (!widget) {
      // GtkWidget *w = create_other_model_tools_dialog();
      GtkWidget *w = widget_from_builder("other_model_tools_dialog");
      graphics_info_t::other_modelling_tools_dialog = w;
      graphics_info_t::set_other_modelling_tools_button_names(w);
      widget = w;
   }
   return widget;
}

void unset_other_modelling_tools_dialog() {

   graphics_info_t::other_modelling_tools_dialog = NULL;
}


void unset_model_fit_refine_dialog() {
   graphics_info_t::model_fit_refine_dialog = NULL;
}

void unset_refine_params_dialog() {
   graphics_info_t::refine_params_dialog = NULL;
}



// We should allow labels that are simply "FWT" and "PHWT" without
// dataset and xtal info.
int valid_labels(const std::string &mtz_file_name, const std::string &f_col,
		 const std::string &phi_col,
		 const std::string &weight_col,
		 bool use_weights) {

   int valid = 0;

   short int have_f = 0;
   short int have_phi = 0;
   short int have_weight = 1; // later turn on test if we have weights.

   std::string f_col_str(f_col);
   std::string phi_col_str(phi_col);
   std::string weight_col_str("");

   if (use_weights)
      weight_col_str = weight_col;

   // These now return have 0 members on failure
   //
//    char **f_cols      = get_f_cols(mtz_file_name, &n_f);
//    char **phi_cols    = get_phi_cols(mtz_file_name, &n_phi);
//    char **weight_cols = get_weight_cols(mtz_file_name, &n_weight);
//    char **d_cols      = get_d_cols(mtz_file_name, &n_d); // anom

   coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);

   // Check first the MTZ column labels that don't have a slash
   for (unsigned int i=0; i<r.f_cols.size(); i++) {
      std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(r.f_cols[i].column_label);
      if (p.second.length() > 0)
	 if (p.second == f_col_str) {
	    have_f = 1;
	    break;
	 }
   }
   for (unsigned int i=0; i<r.phi_cols.size(); i++) {
      std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(r.phi_cols[i].column_label);
      if (p.second.length() > 0)
	 if (p.second == phi_col_str) {
	    have_phi = 1;
	    break;
	 }
   }
   if (use_weights) {
      for (unsigned int i=0; i<r.weight_cols.size(); i++) {
	 std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(r.weight_cols[i].column_label);
	 if (p.second.length() > 0)
	    if (p.second == weight_col_str) {
	       have_weight = 1;
	       break;
	    }
      }
   }


   // Now check the MTZ column labels that *do* have a slash
   if (r.f_cols.size() > 0) {
      for (unsigned int i=0; i< r.f_cols.size(); i++) {
	 if (f_col_str == r.f_cols[i].column_label) {
	    have_f = 1;
	    break;
	 }
      }
   } else {
      std::cout << "ERROR: no f_cols! " << std::endl;
   }

   // We can be trying to make an anomalous fourier.
   if (! have_f) {
      if (r.d_cols.size() > 0) {
	 for (unsigned int i=0; i< r.d_cols.size(); i++) {
	    std::cout << "comparing " << f_col_str << " " << r.d_cols[i].column_label << std::endl;
	    if (f_col_str == r.d_cols[i].column_label) {
	       have_f = 1;
	       break;
	    }
	    std::pair<std::string, std::string> p =
	       coot::util::split_string_on_last_slash(r.d_cols[i].column_label);
	    if (p.second.length() > 0) {
	       if (f_col_str == p.second) {
		  have_f = 1;
		  break;
	       }
	    }
	 }
      }
   }

   if (r.phi_cols.size() > 0) {
      for (unsigned int i=0; i< r.phi_cols.size(); i++) {
	 if (phi_col_str == r.phi_cols[i].column_label) {
	    have_phi = 1;
	    break;
	 }
      }
   } else {
      std::cout << "ERROR: no phi_cols! " << std::endl;
   }

   if (use_weights) {
      have_weight = 0;
      weight_col_str = std::string(weight_col);
      if (r.weight_cols.size() > 0) {
	 for (unsigned int i=0; i< r.weight_cols.size(); i++) {
	    if (weight_col_str == r.weight_cols[i].column_label) {
	       have_weight = 1;
	       break;
	    }
	 }
      } else {
	 std::cout << "ERROR: bad (null) weight_cols! " << std::endl;
      }
   }

   if (have_f && have_phi && have_weight)
      valid = 1;

   if (false)  // debug
      std::cout << "DEBUG:: done checking for valid column labels... returning "
		<< valid << " have-f: " << have_f << " have_phi: " << have_phi << " "
		<< "have_weight: " << have_weight << std::endl;
   return valid;
}

/* We need to know if an mtz file has phases.  If it doesn't then we */
/*  go down a (new 20060920) different path. */
int mtz_file_has_phases_p(const char *mtz_file_name) {

   coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);
//    std::cout << "DEBUG:: mtz_file_has_phases_p: " << mtz_file_name << " has "
// 	     << r.phi_cols.size() << " phasing columns" << std::endl;
   if (r.phi_cols.size() > 0)
      return 1;
   else
      return 0;
}

int is_mtz_file_p(const char *mtz_file_name) {

   // Let's say that we have to find some F columns in an file for it
   // to be counted as an mtz file.
   //
   if (coot::file_exists(mtz_file_name)) {

      coot::mtz_column_types_info_t r = coot::get_mtz_columns(mtz_file_name);
      if (r.f_cols.size() > 0)
	 return 1;
      else
	 return 0;
   } else {
      return 0;
   }
}


int cns_file_has_phases_p(const char *cns_file_name) {

   int r = 0;
   if (coot::file_exists(cns_file_name)) {
      FILE* file = fopen( cns_file_name, "r" );
      if (! file) {
	 std::cout << "WARNING:: oops! failed to open " << cns_file_name << std::endl;
      } else {
	 char buf[4096];
	 for ( int i = 0; i < 4096; i++ ) buf[i] = std::toupper(fgetc(file));
	 fclose( file );
	 buf[4095] = 0;
	 if ( strstr( buf, "ALPHA" ) != NULL && strstr( buf, "BETA"  ) != NULL &&
	      strstr( buf, "GAMMA" ) != NULL && strstr( buf, "SYMOP" ) != NULL &&
	      strstr( buf, " F1="  ) != NULL && strstr( buf, " F2="  ) != NULL )
	    r = 1;
	 else
	    r = 0;
      }
   } else {
      r = 0;
   }
   return r;
}




/*  ----------------------------------------------------------------------- */
/*                        go to atom widget                                 */
/*  ----------------------------------------------------------------------- */

void save_go_to_atom_widget(GtkWidget *widget) { /* store in a static */
   graphics_info_t::go_to_atom_window = widget;
}

void unset_go_to_atom_widget() {
   graphics_info_t::go_to_atom_window = NULL;
}


// not really a button select, its a menu item select
void
save_molecule_coords_combobox_changed(GtkWidget *combobox, gpointer data) {

   // graphics_info_t g;

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));

   std::cout << "INFO:: save_molecule_coords_button_select(): Save coords molecule save_imol now: "
	     << imol << std::endl;

   graphics_info_t::save_imol = imol;
}




/* a c callable wrapper to the graphics_info_t function */
void fill_option_menu_with_coordinates_options(GtkWidget *option_menu,
					       GCallback signal_func,
					       int imol_active_position) {

   graphics_info_t g;
//    g.fill_option_menu_with_coordinates_options(option_menu,
// 					       signal_func,
// 					       imol_active_position);

   std::cout << "100% full of wrongability: fill_option_menu_with_coordinates_options"
	     << std::endl;
}

void fill_combobox_with_coordinates_options(GtkWidget *combobox,
					    GCallback signal_func,
					    int imol_active_position) {
   graphics_info_t g;
   g.fill_combobox_with_coordinates_options(combobox, signal_func, imol_active_position);

}



// void store_refmac_params(const char *fobs_col, const char *sigfobs_col,
// 			 const char *r_free_col, int sensible_f_free_col) {

//    graphics_info_t g;
//    int imol = g.n_molecules;


// }

/*  ----------------------------------------------------------------------- */
/*              new close molecule                                          */
/*  ----------------------------------------------------------------------- */
void
old_new_close_molecules(GtkWidget *window) {

//    // GtkWidget *vbox = lookup_widget(window, "new_delete_molecules_vbox");
//    GtkWidget *vbox = widget_from_builder("new_delete_molecules_vbox");
//    short int closed_something_flag = 0;
//    std::vector<int> closed_molecules;

//    if (GTK_IS_BOX(vbox)) {

// #if (GTK_MAJOR_VERSION >= 4)

//       // 20220602-PE
//       std::cout << "in new_close_molecules() FIXME deleting marked molecules" << std::endl;
// #else
//       GList *dlist = gtk_container_get_children(GTK_CONTAINER(vbox));
//       GList *free_list = dlist;

//       while (dlist) {
//          GtkWidget *list_item = GTK_WIDGET(dlist->data);
//          if (GTK_IS_TOGGLE_BUTTON(list_item)) {
//             if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(list_item))) {
//                int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(list_item), "imol"));
//                closed_molecules.push_back(imol);
//             }
//          } else {
//             std::cout << "not a toggle button" << std::endl;
//          }
//          dlist = dlist->next;
//       }
//       g_list_free(free_list);
// #endif

//    }

//    if (! closed_molecules.empty()) {
//       for (const auto &imol : closed_molecules) {
//          graphics_info_t::molecules[imol].close_yourself();
//       }
//    }

//    // update go to atom molecule now that we may have deleted the
//    // currently set one.
//    if (! closed_molecules.empty()) {
//       graphics_info_t g;
//       for (unsigned int i=0; i<closed_molecules.size(); i++) {
// 	 if (closed_molecules[i] == g.go_to_atom_molecule()) {
// 	    // set it to the bottom model molecule:
// 	    for (int imol=graphics_info_t::n_molecules()-1; imol>=0; imol--) {
// 	       if (is_valid_model_molecule(imol)) {
// 	          g.set_go_to_atom_molecule(imol);
// 	          break;
//                }
//             }
//          }
//       }
//       closed_something_flag = true;
//    }

//    // ------ here ----- if there is a sequence view of a closed molecule being displayed then
//    // hide it. - that should be done in close_yourself().


//    if (closed_something_flag) {
//       if (graphics_info_t::go_to_atom_window) {
// 	 graphics_info_t g;
// 	 // GtkWidget *combobox = lookup_widget(graphics_info_t::go_to_atom_window, "go_to_atom_molecule_combobox");
// 	 GtkWidget *combobox = widget_from_builder("go_to_atom_molecule_combobox");
// 	 int gimol = g.go_to_atom_molecule();

// 	 GCallback callback_func = G_CALLBACK(graphics_info_t::go_to_atom_mol_combobox_changed);
// 	 g.fill_combobox_with_coordinates_options(combobox, callback_func, gimol);

//       }
//       graphics_draw();
//    }
}


void
close_molecules_gtk4(GtkWidget *dialog) {

   GtkWidget *vbox = widget_from_builder("new_delete_molecules_vbox");
   if (vbox) {

      std::vector<int> closed_molecules;
      GtkWidget *item_widget = gtk_widget_get_first_child(vbox);
      while (item_widget) {

         if (gtk_check_button_get_active(GTK_CHECK_BUTTON(item_widget))) {
            int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget), "imol"));
            closed_molecules.push_back(imol);
         }
         item_widget = gtk_widget_get_next_sibling(item_widget);
      };

      if (! closed_molecules.empty()) {
         for (const auto &imol : closed_molecules) {
            graphics_info_t::molecules[imol].close_yourself();
         }
      }
   }
   graphics_info_t g;
   g.update_scroll_wheel_map_on_molecule_close();
   graphics_draw();
}

GtkWidget *wrapped_create_new_close_molecules_dialog() {

   GtkWidget *dialog = widget_from_builder("new_close_molecules_dialog");
   GtkWidget *vbox   = widget_from_builder("new_delete_molecules_vbox"); // nice and consistent...

   graphics_info_t::clear_out_container(vbox);

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model() ||
	  graphics_info_t::molecules[imol].has_xmap() ||
	  graphics_info_t::molecules[imol].has_nxmap()) { // NXMAP-FIXME, check
	 std::string button_name("delete_molecule_checkbutton_");
	 std::string mol_name("   ");
	 mol_name += graphics_info_t::int_to_string(imol);
	 mol_name += "  ";
	 mol_name += graphics_info_t::molecules[imol].name_for_display_manager();
	 button_name += graphics_info_t::int_to_string(imol);
         GtkWidget *checkbutton = gtk_check_button_new_with_label(mol_name.c_str());
         g_object_set_data(G_OBJECT(checkbutton), "imol", GINT_TO_POINTER(imol));
	 gtk_widget_set_visible(checkbutton, TRUE);
         gtk_box_append(GTK_BOX(vbox), checkbutton);
      }
   }
   return dialog;
}


            // try deleting this at some stage

/*  ----------------------------------------------------------------------- */
/*              old close molecule                                          */
/*  ----------------------------------------------------------------------- */
/* get the molecule to delete from the optionmenu */
void
close_molecule_by_widget(GtkWidget *optionmenu) {

   std::cout << "GTK-FIXME no gtk_option_menu_get_menu" << std::endl;

   /*
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
   GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));

   if (active_item) {
      int *imol_p = (int *)g_object_get_data(G_OBJECT(active_item));

      if (imol_p) {

	 std::cout << " Closing molecule number " << *imol_p << std::endl;
	 graphics_info_t g;
	 if (is_valid_model_molecule(*imol_p) || is_valid_map_molecule(*imol_p)) {
	    // g.molecules[*imol_p].close_yourself();
	    close_molecule(*imol_p);
	    fill_close_option_menu_with_all_molecule_options(optionmenu);
	    graphics_draw();

	 } else {
	    std::cout << "ERROR: Closing invalid molecule number" << *imol_p
		      << std::endl;
	 }
      } else {
	 std::cout << "ERROR: trapped an error in close_molecule!\n";
	 std::cout << "       This should never happen!\n";
      }
   } else {
      std::cout << "WARNING:: menu has no active item - nothing to close.\n";
      }
   */
}


void close_molecule(int imol) {

   graphics_info_t g;
   int old_go_to_atom_molecule = g.go_to_atom_molecule();
   bool was_map = false;
   if (is_valid_map_molecule(imol))
      was_map = true;

   // we put this here so that this deletes display manager
   // frames/molecules that are not valid models. For example:
   // 1) open display manager
   // 2) add an MG,
   // 3) delete an MG
   // 4) Try to use the "Delete model" button in the display manager
   //    -> Fail (i.e. nothing happens) without this line
   //
   // g.delete_molecule_from_from_display_manager(imol, was_map);

   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      g.delete_pointers_to_map_in_other_molecules(imol);
      g.molecules[imol].close_yourself();
      // and close the graphics ligand view if it was a residue of this molecule
      g.close_graphics_ligand_view_for_mol(imol);
   }
   if (! was_map) {
      int go_to_atom_imol_new = g.update_go_to_atom_molecule_on_go_to_atom_molecule_deleted();
      if (graphics_info_t::go_to_atom_window) {
         // std::cout << ".....re fill go to atom window here" << std::endl;
         if (imol == old_go_to_atom_molecule) {
	    g.update_go_to_atom_window_on_other_molecule_chosen(go_to_atom_imol_new);
	    g.update_go_to_atom_window_on_changed_mol(go_to_atom_imol_new);
         }
      }
   }

   g.clear_up_moving_atoms_maybe(imol);
   g.update_scroll_wheel_map_on_molecule_close();

   graphics_draw();
   std::string cmd = "close-molecule";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);

}



void fill_close_option_menu_with_all_molecule_options(GtkWidget *optionmenu) {

   std::cout << "GTK-FIXME no gtk_option_menu_get_menu" << std::endl;

   /*
   graphics_info_t g;
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));

   // Strangely enough, you can't set active for menuitems of menus
   // that have already been attached to optionmenus.  So we get this
   // optionmenu, delete it, create a new menu, and add menu items to
   // that (as we did before), set active and *then* attach the menu
   // to the optionmenu.

   gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   GtkWidget *menuitem;

   for (int imol=0; imol<g.n_molecules(); imol++) {

      if (g.molecules[imol].atom_sel.n_selected_atoms > 0 ||
	  g.molecules[imol].has_xmap() || g.molecules[imol].has_nxmap()) {
	 char s[200];
	 snprintf(s,199,"%d",imol);
	 std::string ss(s);
	 ss += " ";
	 ss += g.molecules[imol].name_;
	 menuitem = gtk_menu_item_new_with_label (ss.c_str());
	 g_signal_connect (G_OBJECT (menuitem), "activate",
			   GTK_SIGNAL_FUNC(close_molecule_item_select),
			   GINT_TO_POINTER(imol));
	 int *ip = new int;
	 *ip = imol;
	 gtk_object_set_user_data(GTK_OBJECT(menuitem), ip );
	 gtk_menu_shell_append(GTK_MENU_SHELL(menu), menuitem);
	 gtk_widget_set_visible(menuitem, TRUE);
      }
   }

   if (g.n_molecules() > 0) {
      gtk_menu_set_active(GTK_MENU(menu), 0);
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu),
			    menu);
   */

}



void
close_molecule_item_select(GtkWidget *item, GtkPositionType pos) {

   std::cout << "DEBUG:: activating closing position/imol " << pos << std::endl;

}

//
char* get_text_for_density_size_widget() {

   // convert the density in the graphics_info_t

   graphics_info_t g;
   char *text;

   // we are interfacing with a c function, so we will need to malloc
   // the space for the returned char*.
   //
   // The c function should delete it.
   //
   // Dontcha just *love" this sort coding! (yeuch)
   //
   text = (char *) malloc(100);
   snprintf(text,100,"%-5.1f", g.box_radius_xray);

   return text;

}

char *
get_text_for_density_size_em_widget() {
   graphics_info_t g;
   char *text = (char *) malloc(100);
   snprintf(text,100,"%-5.1f", g.box_radius_em);
   return text;
}

GtkWidget *wrapped_create_show_symmetry_window() {

   // GtkWidget *show_symm_window = create_show_symmetry_window();
   GtkWidget *show_symm_window = widget_from_builder("show_symmetry_window");

   /* Symmetry Search Radius Entry */
   GtkWidget *entry;
   char *text;
   int imol = -1;

   for (int ii=0; ii<graphics_n_molecules(); ii++) {
      if (is_valid_model_molecule(ii)) {
	 imol = ii;
	 break;
      }
      if (is_valid_map_molecule(ii)) {
	 imol = ii;
	 break;
      }
   }

   /* The Show Symmetry RadioButtons */

   // GtkCheckButton *button = nullptr;
   // if (get_show_symmetry() == 1) {
   //    button = GTK_CHECK_BUTTON(widget_from_builder("show_symmetry_yes_radiobutton"));
   // } else {
   //    button = GTK_CHECK_BUTTON(widget_from_builder("show_symmetry_no_radiobutton"));
   // }

   // gtk_check_button_set_active(GTK_CHECK_BUTTON(button), TRUE);

   // 20230516-PE now we use a switch

   GtkWidget *switch_button = widget_from_builder("show_symmetry_switch");
   if (get_show_symmetry() == 1) {
      gtk_switch_set_active(GTK_SWITCH(switch_button), TRUE);
   } else {
      gtk_switch_set_active(GTK_SWITCH(switch_button), FALSE);
   }

#if 0 // 20230513-PE not now

   /* Colour Merge */

   GtkScale *hscale = GTK_SCALE(widget_from_builder("hscale_symmetry_colour"));
   GtkAdjustment *adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0.5, 0.0, 3.0, 0.02, 0.05, 2.0));
   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   g_signal_connect(G_OBJECT (adjustment), "value_changed",
		    G_CALLBACK(symmetry_colour_adjustment_changed),
		    NULL);

#endif

   /*  The Symmetry Search Radius Entry */

    entry = widget_from_builder("symmetry_radius_entry");
    text = get_text_for_symmetry_size_widget(); /* const gchar *text */
    gtk_editable_set_text(GTK_EDITABLE(entry), text);
    free (text);

    /* The Unit Cell Radiobuttons - pre GtkSwitch*/

    // If (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
    //    GtkWidget *check_button = nullptr;
    //    if (get_show_unit_cell(imol) == 1) {
    //       check_button = widget_from_builder("unit_cell_yes_radiobutton");
    //    } else {
    //       check_button = widget_from_builder("unit_cell_no_radiobutton");
    //    }
    //    gtk_check_button_set_active(GTK_CHECK_BUTTON(check_button), TRUE);
    // }


    if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
       GtkWidget *switch_button = widget_from_builder("show_unit_cell_switch");
       if (get_show_unit_cell(imol) == 1)
          gtk_switch_set_active(GTK_SWITCH(switch_button), TRUE);
       else
          gtk_switch_set_active(GTK_SWITCH(switch_button), FALSE);
    }


    //  The Expanded Atoms Label checkbutton

    // checkbutton = lookup_widget(show_symm_window, "show_symmetry_expanded_labels_checkbutton");
    GtkWidget *checkbutton = widget_from_builder("show_symmetry_expanded_labels_checkbutton");
    if (graphics_info_t::symmetry_atom_labels_expanded_flag)
       gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);

#if 0
    // GtkWidget *colour_button = lookup_widget(show_symm_window, "symmetry_colorbutton");
    GtkWidget *colour_button = widget_from_builder("symmetry_colorbutton"); // a GtkButton
    if (colour_button) {
       GdkRGBA bg_colour;
       bg_colour.red   = (guint)(graphics_info_t::symmetry_colour[0] * 65535);
       bg_colour.green = (guint)(graphics_info_t::symmetry_colour[1] * 65535);
       bg_colour.blue  = (guint)(graphics_info_t::symmetry_colour[2] * 65535);
       // gtk_color_button_set_color(GTK_COLOR_BUTTON(colour_button), &bg_colour);
       // gtk_color_button_set_rgba(GTK_COLOR_BUTTON(colour_button), &bg_colour);
       gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(colour_button), &bg_colour);
    } else {
       std::cout << "failed to lookup colourbutton" << std::endl;
    }
#endif

    GtkWidget *box_for_colour_button = widget_from_builder("hbox_for_colour_button");
    if (box_for_colour_button) {
       GtkWidget *child_item = gtk_widget_get_first_child(box_for_colour_button);
       if (child_item) {
          // the colour button has already been added
       } else {

#if 0 // 20230513-PE color dialog is not in GTK 4.4.0 (it is in 4.10)
          GtkWidget *col_dialog = gtk_color_dialog_new();
          GtkWidget *colour_button_dialog = gtk_color_dialog_button_new(col_dialog);
          gtk_box_append(GTK_BOX(box_for_colour_button), colour_button_dialog);
#endif

          auto on_color_set_func = +[] (GtkColorButton *self, gpointer user_data) {
             GdkRGBA rgba;
             gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(self), &rgba);
             std::cout << "Selected color: " << gdk_rgba_to_string(&rgba) << std::endl;
             graphics_info_t::rgba_to_symmetry_colour(rgba);
             graphics_info_t::update_symmetry();
             graphics_info_t::graphics_draw();
          };

          GdkRGBA rgba = graphics_info_t::symmetry_colour_to_rgba();
          if (false)
             std::cout << " colours " << rgba.red << " " << rgba.green << " " << rgba.blue << std::endl;
          GtkWidget *colour_button = gtk_color_button_new_with_rgba(&rgba);
          gtk_box_append(GTK_BOX(box_for_colour_button), colour_button);
          g_signal_connect(G_OBJECT(colour_button), "color-set", G_CALLBACK(on_color_set_func), nullptr);
       }
    }

//     // The symmetry colour molecule checkbutton
//     checkbutton = lookup_widget(show_symm_window,
// 				"show_symmetry_molecule_rotate_colour_map_checkbutton");

//     if (graphics_info_t::symmetry_rotate_colour_map_flag) {
//        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
//     }


//     // The symmetry colour symop checkbutton
//     checkbutton = lookup_widget(show_symm_window,
// 				"show_symmetry_colour_by_symop_checkbutton");

//     if (graphics_info_t::symmetry_colour_by_symop_flag) {
//        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
//     }

//     // Whole chain checkbutton
//     checkbutton = lookup_widget(show_symm_window,
// 				"show_symmetry_whole_molecule_checkbutton");

//     if (graphics_info_t::symmetry_whole_chain_flag) {
//        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
//     }

    return show_symm_window;

}


void symmetry_colour_adjustment_changed (GtkAdjustment *adj,
					 GtkWidget *window) {

   float f = gtk_adjustment_get_value(adj);

   // does a graphics_draw() for us...
   set_symmetry_colour_merge(f); /* this adjusts graphics_info_t::
				    symm_colour_merge_weight, which is a double
				    array.  But we only ever use the 0th
				    position of it in combine_colour() */
}


GtkWidget *symmetry_molecule_controller_dialog() {

   graphics_info_t g;
   return g.wrapped_create_symmetry_controller_dialog();
}


/* used by destroy callback, needed because there should only be one of these. */
void set_symmetry_controller_dialog_widget(GtkWidget *w) {

   graphics_info_t g;
   g.symmetry_controller_dialog = w;
}


void
handle_map_colour_change(int imol, GdkRGBA map_col) {

   graphics_info_t::molecules[imol].handle_map_colour_change(map_col,
                                                             graphics_info_t::swap_difference_map_colours,
                                                             graphics_info_t::GL_CONTEXT_MAIN,
                                                             graphics_info_t::get_rotation_centre_co(),
                                                             graphics_info_t::box_radius_xray);

   if (false) { // for the moment
      if (graphics_info_t::display_mode_use_secondary_p()) {
         graphics_info_t g;
         g.make_gl_context_current(graphics_info_t::GL_CONTEXT_SECONDARY);
         g.molecules[imol].handle_map_colour_change(map_col,
                                                    g.swap_difference_map_colours,
                                                    graphics_info_t::GL_CONTEXT_SECONDARY,
                                                    graphics_info_t::get_rotation_centre_co(),
                                                    graphics_info_t::box_radius_xray);
         g.make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
      }
   }

   //cout << "using map colours:   "
   //	<< map_col[0] << " "  << map_col[1] << " "  << map_col[2] << endl;
   graphics_draw();
}

void set_last_map_colour(double f1, double f2, double f3) {

   graphics_info_t g;
   g.set_last_map_colour(f1, f2, f3);
}

void set_map_colour(int imol, float red, float green, float blue) {

   if (is_valid_map_molecule(imol)) {
      GdkRGBA colour;

      // why were these multipliers here - what were they doing?
      // 0.0 to 1.0 is the range.

      // colour.red   = red   * 65535.0;
      // colour.green = green * 65535.0;
      // colour.blue   = blue * 65535.0;

      colour.red   = red;
      colour.green = green;
      colour.blue  = blue;

      short int swap_col = graphics_info_t::swap_difference_map_colours;
      graphics_info_t::molecules[imol].handle_map_colour_change(colour, swap_col,
                                                                graphics_info_t::GL_CONTEXT_MAIN,
                                                                graphics_info_t::get_rotation_centre_co(),
                                                                graphics_info_t::box_radius_xray);
      if (graphics_info_t::display_mode_use_secondary_p()) {
         graphics_info_t g;
         g.make_gl_context_current(graphics_info_t::GL_CONTEXT_SECONDARY);
         graphics_info_t::molecules[imol].handle_map_colour_change(colour, swap_col,
                                                                   graphics_info_t::GL_CONTEXT_SECONDARY,
                                                                   graphics_info_t::get_rotation_centre_co(),
                                                                   graphics_info_t::box_radius_xray);
         g.make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
      }

      graphics_draw();
   }
}

/*! \brief set the colour of the imolth map using a (7-character) hex colour */
void set_map_hexcolour(int imol, const char *hex_colour) {

   coot::colour_holder ch(hex_colour);
   set_map_colour(imol, ch.red, ch.green, ch.blue);

}



// void add_on_map_colour_choices(GtkWidget *menu) {

//    // GtkWidget *sub_menu = lookup_widget(menu, sub_menu_name.c_str());
//    // GtkWidget *sub_menu = widget_from_builder(sub_menu_name); // No, because it's dynamically added
//    //                                                           // in create_initial_map_color_submenu()

//    GtkWidget *sub_menu = gtk_menu_item_get_submenu(GTK_MENU_ITEM(menu));
//    if (!sub_menu) {
//       std::cout << "ERROR:: in add_on_map_colour_choices() sub menu map_colour1_menu not found in add_on_map_colour_choices()\n";
//    } else {
//       gtk_container_foreach(GTK_CONTAINER(sub_menu),
//                             my_delete_menu_items,
//                             (gpointer) sub_menu);
//       GCallback callback = G_CALLBACK(map_colour_mol_selector_activate);
//       for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
//          if (graphics_info_t::molecules[imol].has_xmap() ||
//              graphics_info_t::molecules[imol].has_nxmap()) { // NXMAP-FIXME
//             std::string name = graphics_info_t::molecules[imol].dotted_chopped_name();
//             add_map_colour_mol_menu_item(imol, name, sub_menu, callback);
//          }
//       }
//    }
// }

// 20220602-PE don't call this function
#if 0
void my_delete_menu_items(GtkWidget *widget, void *user_data) {
   gtk_container_remove(GTK_CONTAINER(data), widget);
}
#endif



// where is this called from?

// void map_colour_mol_selector_activate(GtkMenuItem     *menuitem,
//                                       gpointer         user_data) {

//    int imol = GPOINTER_TO_INT(user_data);
//    show_map_colour_selector(imol);

// }


// ---------------------------------------------------------
// Scroll wheel, similar
// ---------------------------------------------------------
//

int scroll_wheel_map() {
   return graphics_info_t::scroll_wheel_map;
}

void set_scroll_wheel_map(int imap) {

   if (is_valid_map_molecule(imap)) {
      graphics_info_t g;
      g.set_scrollable_map(imap);
   }
}

void clear_out_container(GtkWidget *vbox); // in c-interface.cc

GtkWidget *wrapped_create_bond_parameters_dialog() {

   // move this into graphics_info_t I think

   graphics_info_t g;

   // GtkWidget *widget = create_bond_parameters_dialog();
   GtkWidget *dialog = widget_from_builder("bond_parameters_dialog");
   GtkWidget *combobox = widget_from_builder("bond_parameters_molecule_comboboxtext");

   // fill the colour map rotation entry

   // check the Carbons only check button

   // fill the molecule bond width option menu

   // check the draw hydrogens check button

   // Consider the case were we set the bond_parameters_molecule and
   // then close that molecule Usually, we want imol to be set to
   // bond_parameters_molecule but in the case of a closed molecule we
   // want bond_parameters_molecule to be set to first_coords_imol()
   // i.e. imol.

   int imol = first_coords_imol(); // can be -1;
   if (g.bond_parameters_molecule >= 0)
      if (g.molecules[g.bond_parameters_molecule].has_model())
	 imol = g.bond_parameters_molecule;
      else
	 g.bond_parameters_molecule = imol;
   else
      // g.bond_parameters_molecule not set yet.
      g.bond_parameters_molecule = imol;

   auto get_model_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_model_molecule(i))
                                           vec.push_back(i);
                                     return vec;
                                  };
   int imol_active = g.bond_parameters_molecule;
   auto model_list = get_model_molecule_vector();
   GCallback callback_func = G_CALLBACK(nullptr);
   g.fill_combobox_with_molecule_options(combobox, callback_func, imol_active, model_list);

   g.fill_bond_parameters_internals(combobox, imol);

   return dialog;
}

void apply_bond_parameters(GtkWidget *w) {

   graphics_info_t g;

   int imol = -1; // 20240713-PE was g.bond_parameters_molecule;

   GtkWidget *bond_parameters_molecule_comboboxtext = widget_from_builder("bond_parameters_molecule_comboboxtext");
   if (bond_parameters_molecule_comboboxtext) {
      imol = g.combobox_get_imol(GTK_COMBO_BOX(bond_parameters_molecule_comboboxtext));
   }

   if (imol >= 0) {
      if (imol < g.n_molecules()) {
	 if (graphics_info_t::molecules[imol].has_model()) {

	    // bond thickness
	    //
	    // note that the bond_thickness_intermediate_value should
	    // be (is?) set to -1 when the widget is created.  The
	    // value is changed when we activate the bond width menu
	    // item.
	    //
	    if (g.bond_thickness_intermediate_value > 0) {
	       set_bond_thickness(imol, g.bond_thickness_intermediate_value);
	    }

	    // intermediate atom bond thickness: We don't want fat
	    // bonds obscuring the intermediate atoms, so we set them
	    // to be 2 pixels fatter than these fat bonds if necesary.
	    //
	    if (g.bond_thickness_intermediate_atoms < (g.bond_thickness_intermediate_value + 2)) {
	       g.bond_thickness_intermediate_atoms =
		  g.bond_thickness_intermediate_value + 2;
	    }

	    // draw hydrogens?

	    GtkWidget *check_button = widget_from_builder("draw_hydrogens_yes_radiobutton");
	    if (gtk_check_button_get_active(GTK_CHECK_BUTTON(check_button))) {
	       set_draw_hydrogens(imol, 1);
	    } else {
	       set_draw_hydrogens(imol, 0);
	    }
	    g.update_environment_distances_by_rotation_centre_maybe(imol);
	 }
      }
   }
   graphics_draw();
}

// void skeletonize_map_by_optionmenu(GtkWidget *optionmenu) {

//    graphics_info_t g;
//    g.skeletonize_map_by_optionmenu(optionmenu);

// }

void skeletonize_map_by_combobox(GtkWidget *combobox) {

   graphics_info_t g;
   g.skeletonize_map_by_combobox(combobox);
}


void
skeletonize_map_single_map_maybe(GtkWidget *window, int imol) {

   // GtkWidget *on_radio_button = lookup_widget(window, "single_map_skeleton_on_radiobutton");
   GtkWidget *on_radio_button = widget_from_builder("single_map_skeleton_on_radiobutton");

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(on_radio_button))) {

      graphics_info_t::skeletonize_map(imol, 0);
      if (graphics_info_t::map_for_skeletonize < 0) {
	 // it was unset, so set it...
	 graphics_info_t::map_for_skeletonize = imol;
      }
   } else {
      graphics_info_t::unskeletonize_map(imol);
   }
}

void set_file_for_save_filechooser(GtkWidget *fileselection) {

   graphics_info_t g;
   g.set_file_for_save_filechooser(fileselection);
}



GtkWidget *wrapped_create_skeleton_dialog() {

   // from the menu item callback, we don't need to display the ca_mode label.
   bool display_ca_model_label = 0;
   graphics_info_t g;
   return g.wrapped_create_skeleton_dialog(display_ca_model_label);
}


void save_coordinates_using_widget(GtkWidget *dialog) {

   // the widget that we get passed is the filechooser dialog
   // the data was set in on_save_coords_dialog_save_button_clicked.

   {

      gpointer data = g_object_get_data(G_OBJECT(dialog), "imol");
      int imol = GPOINTER_TO_INT(data);

      bool save_hydrogens = 1;
      bool save_aniso_records = 1;

      // get the filename?

      // GtkWidget *chk_but = lookup_widget(GTK_WIDGET(widget), "checkbutton_hydrogens");
      GtkWidget *chk_but = widget_from_builder("checkbutton_hydrogens");
      if (! gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(chk_but)))
	 save_hydrogens = 0;
      chk_but = widget_from_builder("checkbutton_aniso");
      if (! gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(chk_but)))
	 save_aniso_records = 0;

      // const gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(widget));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *filename = g_file_info_get_name(file_info);

      std::cout << "INFO:: save coordinates for molecule " << imol << " to file " << filename << std::endl;

      graphics_info_t g;
      if (is_valid_model_molecule(imol)) {
         int save_conect_records = g.write_conect_records_flag;
	 int ierr = g.molecules[imol].save_coordinates(filename, save_hydrogens, save_aniso_records, save_conect_records);
	 if (! ierr) {
	    std::string s = "Saved coordinates file ";
	    s += filename;
	    s += ".";
	    g.add_status_bar_text(s);
	 }
      }
   }
}

void save_symmetry_coords_from_filechooser(GtkWidget *filechooser) {

   coot::Symm_Atom_Pick_Info_t *symm_info = (coot::Symm_Atom_Pick_Info_t *) g_object_get_data(G_OBJECT(filechooser), "symm_info");

   // const gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));

   std::cout << "debug:: symm_info: " << symm_info << std::endl;
   std::cout << "debug:: symm_info->imol: " << symm_info->imol << std::endl;

   GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(filechooser));
   if (file) {
      const char *file_name = g_file_get_path(file);
      if (file_name) {
         if (symm_info) {
            // std::cout << "Preshift to origin:  " << symm_info->pre_shift_to_origin << std::endl;
            save_symmetry_coords(symm_info->imol,
                                 file_name,
                                 symm_info->symm_trans.isym(),
                                 symm_info->symm_trans.x(),
                                 symm_info->symm_trans.y(),
                                 symm_info->symm_trans.z(),
                                 symm_info->pre_shift_to_origin.us,
                                 symm_info->pre_shift_to_origin.vs,
                                 symm_info->pre_shift_to_origin.ws);
         }
      }
   } else {
      logger.log(log_t::WARNING, logging::function_name_t("save_symmetry_coords_from_filechooser"),
                 "No file");
   }
}


GtkWidget *wrapped_create_goto_atom_window() {

   GtkWidget *widget = graphics_info_t::go_to_atom_window;
   if (widget) {

#if (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME widget raise
#else
      if (!gtk_widget_get_mapped(widget))
         gtk_widget_set_visible(widget, TRUE);
      else
         gdk_window_raise(GDK_WINDOW(gtk_widget_get_window(widget))); // can I just cast it like this?
      // std::cout << "GTK-FIXME no raise wrapped_create_goto_atom_window()" << std::endl;
#endif

   } else {
      // widget = create_goto_atom_window();
      widget = widget_from_builder("goto_atom_window");
      graphics_info_t::go_to_atom_window = widget;
      if (graphics_info_t::go_to_atom_window_x_position > -1) {
         if (false) // 20220315-PE too much noise - I should fix this later
            std::cout << "GTK-FIXME no gtk_widget_set_uposition F2" << std::endl;

         // GTK3 - you can't set widget/window positions.
         // 	 gtk_widget_set_uposition(widget,
         // 				  graphics_info_t::go_to_atom_window_x_position,
         // 				  graphics_info_t::go_to_atom_window_y_position);

      }
      graphics_info_t g;
      g.fill_go_to_atom_window_gtk3(widget);
   }
   return widget;
}


void post_go_to_atom_window() {

   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *widget = wrapped_create_goto_atom_window();
      gtk_widget_set_visible(widget, TRUE);
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("post-go-to-atom-window");
   add_to_history(command_strings);
}

void fill_go_to_atom_window(GtkWidget *widget) {

   graphics_info_t g;
   g.fill_go_to_atom_window_gtk3(widget);

}


/* used by keypress (return) callbacks */
//
// read the widget values and apply them to the graphics
//
int apply_go_to_atom_values(GtkWidget * window) {

  const gchar *chain_str;
  const gchar *res_str;
  const gchar *atom_name_str;

  GtkEntry *entry;
				/* FIXME */

				/* Dealloc the strings here */

  entry = GTK_ENTRY(widget_from_builder("go_to_atom_chain_entry"));
  chain_str = gtk_editable_get_text(GTK_EDITABLE(entry));

  entry = GTK_ENTRY(widget_from_builder("go_to_atom_residue_entry"));
  res_str = gtk_editable_get_text(GTK_EDITABLE(entry));

  entry = GTK_ENTRY(widget_from_builder("go_to_atom_atom_name_entry"));
  atom_name_str = gtk_editable_get_text(GTK_EDITABLE(entry));

  set_go_to_atom_chain_residue_atom_name_strings(chain_str, res_str, atom_name_str);

  return 0;
}


// Delete me if compiling with GTK1 works.  It seems that nothing
// calls this function.
//
// // called by a function in callback.c
// //
// // void on_go_to_atom_residue_list_selection_changed (GtkList         *gtklist,
// // 						   gpointer         user_data) {

// void on_go_to_atom_residue_list_selection_changed (GtkList         *gtklist,
// 						   gpointer         user_data) {

// //    graphics_info_t g;
// //    g.on_go_to_atom_residue_list_selection_changed(gtklist, user_data);

//    // old residue list function stub.
// }



void clear_atom_list(GtkWidget *atom_gtklist) {

}

// void on_go_to_atom_residue_list_select_child (GtkList         *list,
// 					      GtkWidget       *widget,
// 					      gpointer         user_data) {
//    std::cout << "child selected.\n";
// }

// void on_go_to_atom_residue_list_unselect_child (GtkList         *list,
// 						GtkWidget       *widget,
// 						gpointer         user_data) {
//    std::cout << "child unselected.\n";
// }


void apply_go_to_atom_from_widget(GtkWidget *widget) {

   graphics_info_t g;
   g.apply_go_to_atom_from_widget(widget);
}

// This is a different action to wrapped_create_goto_atom_window,
// because we don't want to raise an already existing dialog, we need
// to destroy it.
//
GtkWidget *wrapped_create_residue_info_dialog() {

   GtkWidget *widget = graphics_info_t::residue_info_dialog;
   if (widget) {
      // raise/uniconify (or whatever) what we have:
      //
// not this widget...
//       if (!GTK_WIDGET_MAPPED(widget))
// 	 gtk_widget_set_visible(widget, TRUE);
//       else
// 	 gdk_window_raise(widget->window);

      widget = widget_from_builder("residue_info_dialog");
      graphics_info_t::residue_info_dialog = widget;

   } else {

      // create (then store) a new one.

      widget = widget_from_builder("residue_info_dialog");
      graphics_info_t::residue_info_dialog = widget;
   }
   return widget;

}

int residue_info_dialog_is_displayed() {

   int r = 0;
   if (graphics_info_t::residue_info_dialog)
      r = 1;
   return r;
}

// #include "c-interface-gui.hh"

void ideal_nucleic_acid_by_widget(GtkWidget *builder_dialog) {

   std::string type = "RNA";
   std::string form = "A";
   short int single_stranded_flag = 0;
   // GtkWidget *entry = lookup_widget(builder_dialog, "nucleotide_sequence");
   GtkWidget *entry = widget_from_builder("nucleotide_sequence");

   // GtkWidget *type_combobox   = lookup_widget(builder_dialog, "nucleotide_builder_type_combobox");
   // GtkWidget *form_combobox   = lookup_widget(builder_dialog, "nucleotide_builder_form_combobox");
   // GtkWidget *strand_combobox = lookup_widget(builder_dialog, "nucleotide_builder_strand_combobox");
   GtkWidget *type_combobox   = widget_from_builder("nucleotide_builder_type_combobox");
   GtkWidget *form_combobox   = widget_from_builder("nucleotide_builder_form_combobox");
   GtkWidget *strand_combobox = widget_from_builder("nucleotide_builder_strand_combobox");

   type = get_active_label_in_combobox(GTK_COMBO_BOX(type_combobox));
   form = get_active_label_in_combobox(GTK_COMBO_BOX(form_combobox));
   std::string strand = get_active_label_in_combobox(GTK_COMBO_BOX(strand_combobox));
   if (strand == "Single")
      single_stranded_flag = 1;
   const char *txt = gtk_editable_get_text(GTK_EDITABLE(entry));
   if (txt) {
      ideal_nucleic_acid(type.c_str(), form.c_str(), single_stranded_flag, txt);
   }
}



GtkWidget *wrapped_create_display_control_window() {

   graphics_info_t g;
   return g.wrapped_create_display_control_window();
}


void
align_labels_checkbutton_toggled(GtkToggleButton *togglebutton) {

   float align = 0.0;
   if (gtk_toggle_button_get_active(togglebutton))
      align = 1.0;

   graphics_info_t g;

   int n_mols = graphics_info_t::n_molecules();
   for (int i=0; i<n_mols; i++) {
      if (is_valid_model_molecule(i)) {
         std::string name_stub = "display_mol_entry_";
         std::string name = name_stub + coot::util::int_to_string(i);
         // GtkWidget *entry = lookup_widget(g.display_control_window(), name.c_str());
         GtkWidget *entry = 0; // 20220309-PE fixme one day. Not today
         if (entry) {
            // 20180304
            // This only changes the alignment for entries that have smaller text than
            // the widget.  I want all widgets to adjust their text.
            // Maybe use PangoLayout. read gtkentry.c to see if you can work out what
            // happens when the user does a click-drag (left or right) on the text.
            // So I will make the checkbutton invisible for now.
            //
            // gtk_entry_set_alignment(GTK_ENTRY(entry), align);

            // gtk_misc_set_alignment(GTK_MISC(entry), align, 0.5); // no. an entry is not a misc.
         }
      }
   }
}



// BL things for file_chooser
void set_file_chooser_selector(int istate) {
   graphics_info_t::gtk2_file_chooser_selector_flag = istate;
}

int file_chooser_selector_state(){
   return graphics_info_t::gtk2_file_chooser_selector_flag;
}

void set_file_chooser_overwrite(int istate) {
   graphics_info_t::gtk2_chooser_overwrite_flag = istate;
}

int file_chooser_overwrite_state(){
   return graphics_info_t::gtk2_chooser_overwrite_flag;
}


void export_map_gui(short int export_map_fragment) {

   // GtkWidget *w = create_export_map_dialog();

   // this is the widget that chooses the map molecule, not the file chooser.
   //
   GtkWidget *w = widget_from_builder("export_map_frame");

   GtkWidget *hbox = widget_from_builder("export_map_fragment_hbox");
   gtk_widget_set_visible(hbox, export_map_fragment == 0 ? FALSE : TRUE);

   // GtkWidget *combobox = lookup_widget(w, "export_map_map_combobox");
   GtkWidget *combobox = widget_from_builder("export_map_map_combobox");

   graphics_info_t g;

   // g.fill_option_menu_with_map_options(option_menu, NULL);

   // we don't want to do anything when the menu is
   // pressed. We do want to know what the active
   // item was.

   g_object_set_data(G_OBJECT(w), "is_map_fragment", GINT_TO_POINTER(export_map_fragment));
   int imol_active = imol_refinement_map();
   g.fill_combobox_with_map_options(combobox, NULL, imol_active);
   gtk_widget_set_visible(w, TRUE);

}

void on_export_map_dialog_ok_button_clicked_cc(GtkButton *button) {

   // GtkWidget *w = lookup_widget(GTK_WIDGET(button), "export_map_dialog");
   GtkWidget *w = widget_from_builder("export_map_dialog");

   // GtkWidget *combobox = lookup_widget(GTK_WIDGET(button), "export_map_map_combobox");
   GtkWidget *combobox = widget_from_builder("export_map_map_combobox");

   // GtkWidget *text_entry  = lookup_widget(GTK_WIDGET(button), "export_map_radius_entry");
   GtkWidget *text_entry = widget_from_builder("export_map_radius_entry");
   int is_map_fragment = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "is_map_fragment"));
   const gchar *entry_text = gtk_editable_get_text(GTK_EDITABLE(text_entry));

   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(combobox));

   if (true) {

      // GtkWidget *file_chooser_dialog = create_export_map_filechooserdialog();
      GtkWidget *file_chooser_dialog = widget_from_builder("export_map_file_chooser_dialog");
      unsigned int l = std::string(entry_text).length();
      char *c = new char [l + 1];
      strncpy(c, entry_text, l+1);
      g_object_set_data(G_OBJECT(file_chooser_dialog), "is_map_fragment",  GINT_TO_POINTER(is_map_fragment));
      g_object_set_data(G_OBJECT(file_chooser_dialog), "export_map_radius_entry_text",  c);
      g_object_set_data(G_OBJECT(file_chooser_dialog), "map_molecule_number",  GINT_TO_POINTER(imol_map));
      set_transient_and_position(COOT_UNDEFINED_WINDOW, file_chooser_dialog);
      gtk_widget_set_visible(file_chooser_dialog, TRUE);
   }

   gtk_widget_set_visible(w, FALSE);
}





// we add a universal function to set the file names
// in file chooser or selector

// void set_filename_for_filechooserselection(GtkWidget *fileselection,
// const gchar *filename) {
//
// bool chooser = 0;
   // if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      // chooser = 1;
      // gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(fileselection),
					// filename);
   // }
//
// }


void store_geometry_dialog(GtkWidget *w) {

   graphics_info_t g;
   g.geometry_dialog = w;
}


#if 0 // 20211202-PE we no longer want to do this
void store_fixed_atom_dialog(GtkWidget *w) {
   graphics_info_t::fixed_atom_dialog = w;
}
#endif


/* not for user consumption, this finds (from itself) the residue type
   and calls the graphics_info_t function. */
void fill_chi_angles_vbox(GtkWidget *vbox) {

   graphics_info_t g;
   gchar *strval = (gchar *) g_object_get_data(G_OBJECT(vbox), "strval");
   g.fill_chi_angles_vbox(vbox, strval, graphics_info_t::EDIT_CHI);
}

GtkWidget *wrapped_create_add_additional_representation_gui() {

   // 20220309-PE what should this function do these days?
   // (Just passing through fixing lookup_widgets....)
   // i.e. think about this another time.

   std::cout << "::::::::::::::: wrapped_create_add_additional_representation_gui() " << std::endl;

   GtkWidget *w = 0;
   if (graphics_info_t::use_graphics_interface_flag) {

      // w = create_add_reps_dialog();

      graphics_info_t g;
      w = widget_from_builder("add_reps_dialog");
      
      // GtkWidget *combobox = lookup_widget(w, "add_reps_molecule_combobox");
      GtkWidget *combobox = widget_from_builder("add_reps_molecule_combobox");

      int imol_active_position = g.get_active_atom().first;
      fill_combobox_with_coordinates_options(combobox, nullptr, imol_active_position);

#ifdef HAVE_GTK_COMBO_BOX_GET_ACTIVE_TEXT
      // set the active item to be the 8
      gtk_combo_box_set_active(GTK_COMBO_BOX(add_rep_bond_width_combobox), 7);
#endif
   }
   return w;
}

void add_reps_molecule_combobox_changed(GtkWidget *combobox, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::add_reps_molecule_combobox_molecule = imol;

}

void add_reps_molecule_option_menu_item_select(GtkWidget *item, GtkPositionType pos) {
   graphics_info_t::add_reps_molecule_option_menu_item_select_molecule = pos;
}

void
add_additional_representation_by_dialog(GtkDialog *dialog) {

   GtkWidget *combobox = widget_from_builder("add_reps_molecule_combobox");

   GtkWidget *chain_id_entry    = widget_from_builder("add_rep_chain_id_entry");
   GtkWidget *resno_start_entry = widget_from_builder("add_rep_resno_start_entry");
   GtkWidget *resno_end_entry   = widget_from_builder("add_rep_resno_end_entry");
   GtkWidget *ins_code_entry    = widget_from_builder("add_rep_ins_code_entry");
   GtkWidget *string_selection_entry = widget_from_builder("add_rep_selection_string_entry");

   GtkWidget *position_radiobutton = widget_from_builder("add_rep_radiobutton_position");
   GtkWidget *resno_radiobutton    = widget_from_builder("add_rep_radiobutton_res_number");
   GtkWidget *selection_string_radiobutton = widget_from_builder("add_rep_radiobutton_selection_string");

   GtkWidget *add_reps_fat_bonds_radiobutton = widget_from_builder("add_rep_rep_fat_bonds_radiobutton");
   GtkWidget *add_rep_bond_width_combobox = widget_from_builder("add_rep_bond_width_combobox");
   GtkWidget *add_reps_ball_and_stick_radiobutton = widget_from_builder("add_rep_rep_ball_and_stick_radiobutton");

   int imol_active = -1;

   float bond_width = 8;
   int bonds_box_type = coot::NORMAL_BONDS;
   short int representation_type = coot::SIMPLE_LINES;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(add_reps_ball_and_stick_radiobutton))) {
      representation_type = coot::BALL_AND_STICK;
   }
   bool draw_H_flag = 1;

   graphics_info_t g;
   std::string bond_width_txt = g.get_active_label_in_comboboxtext(GTK_COMBO_BOX_TEXT(combobox));

   if (representation_type == coot::BALL_AND_STICK)
      bond_width = 0.15; // not 8

   GtkWidget *dcw = g.display_control_window();

   // int imol = graphics_info_t::add_reps_molecule_option_menu_item_select_molecule;
   int imol = graphics_info_t::add_reps_molecule_combobox_molecule;

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(position_radiobutton))) {
      std::pair<bool, std::pair<int, coot::atom_spec_t> > aas = active_atom_spec();
      if (aas.first) {
	 int imol_active = aas.second.first;
	 coot::atom_selection_info_t asi(aas.second.second.chain_id,
					 aas.second.second.res_no,
					 aas.second.second.res_no,
					 aas.second.second.ins_code);
         GtkWidget *glarea_0 = 0;
         GtkWidget *glarea_1 = 0;
         if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
         if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
	 gl_context_info_t glci(glarea_0, glarea_1);
	 g.molecules[imol_active].add_additional_representation(representation_type,
								bonds_box_type,
								bond_width,
								draw_H_flag,
								asi, dcw, glci, g.Geom_p());
      }
   }
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(resno_radiobutton))) {
      // std::cout << "By chainid/resno" << std::endl;
      std::string chain_id = gtk_editable_get_text(GTK_EDITABLE(chain_id_entry));
      std::string resno_1s = gtk_editable_get_text(GTK_EDITABLE(resno_start_entry));
      std::string resno_2s = gtk_editable_get_text(GTK_EDITABLE(resno_end_entry));
      std::string ins_code = gtk_editable_get_text(GTK_EDITABLE(ins_code_entry));
      if (is_valid_model_molecule(imol)) {
	 int resno_1 = atoi(resno_1s.c_str());
	 int resno_2 = atoi(resno_2s.c_str());
	 coot::atom_selection_info_t asi(chain_id, resno_1, resno_2, ins_code);
         GtkWidget *glarea_0 = 0;
         GtkWidget *glarea_1 = 0;
         if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
         if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
	 gl_context_info_t glci(glarea_0, glarea_1);
	 graphics_info_t::molecules[imol].add_additional_representation(representation_type,
									bonds_box_type,
									bond_width,
									draw_H_flag,
									asi, dcw, glci, g.Geom_p());
      }
   }
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(selection_string_radiobutton))) {
      // std::cout << "By selection string" << std::endl;
      std::string s = gtk_editable_get_text(GTK_EDITABLE(string_selection_entry));
      coot::atom_selection_info_t asi(s);
      GtkWidget *glarea_0 = 0;
      GtkWidget *glarea_1 = 0;
      if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
      if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
      gl_context_info_t glci(glarea_0, glarea_1);
      graphics_info_t::molecules[imol].add_additional_representation(representation_type,
								     bonds_box_type,
								     bond_width,
								     draw_H_flag,
								     asi, dcw, glci, g.Geom_p());
   }
   graphics_draw();
   

}


// old
void add_additional_representation_by_widget(GtkWidget *dialog) {

   // decode the widgets and add the representation here.

   // GtkWidget *option_menu = lookup_widget(w, "add_rep_molecule_optionmenu");
   GtkWidget *combobox = widget_from_builder("add_reps_molecule_combobox");

   GtkWidget *chain_id_entry    = widget_from_builder("add_rep_chain_id_entry");
   GtkWidget *resno_start_entry = widget_from_builder("add_rep_resno_start_entry");
   GtkWidget *resno_end_entry   = widget_from_builder("add_rep_resno_end_entry");
   GtkWidget *ins_code_entry    = widget_from_builder("add_rep_ins_code_entry");
   GtkWidget *string_selection_entry = widget_from_builder("add_rep_selection_string_entry");

   GtkWidget *position_radiobutton = widget_from_builder("add_rep_radiobutton_position");
   GtkWidget *resno_radiobutton    = widget_from_builder("add_rep_radiobutton_res_number");
   GtkWidget *selection_string_radiobutton = widget_from_builder("add_rep_radiobutton_selection_string");

   GtkWidget *add_reps_fat_bonds_radiobutton = widget_from_builder("add_rep_rep_fat_bonds_radiobutton");
   GtkWidget *add_rep_bond_width_combobox = widget_from_builder("add_rep_bond_width_combobox");
   GtkWidget *add_reps_ball_and_stick_radiobutton = widget_from_builder("add_rep_rep_ball_and_stick_radiobutton");

   int imol_active = 0;

   float bond_width = 8;
   int bonds_box_type = coot::NORMAL_BONDS;
   short int representation_type = coot::SIMPLE_LINES;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(add_reps_ball_and_stick_radiobutton))) {
      representation_type = coot::BALL_AND_STICK;
   }
   bool draw_H_flag = 1;
   gchar* bond_width_text = 0;

#ifdef HAVE_GTK_COMBO_BOX_GET_ACTIVE_TEXT
   bond_width_text = gtk_combo_box_get_active_text(GTK_COMBO_BOX(add_rep_bond_width_combobox));
#else
   // 2.4 does not have gtk_combo_box_get_active_text() and 1.x does
   // not have gtk_combo_box.
#endif

   if (bond_width_text) {
      bond_width = atof(bond_width_text);
   } else {
      // Not currently an error.  It happens in gtk1 version.
      // std::cout << "ERROR:: null bond_width_text, using default of 8" << std::endl;
   }

   if (representation_type == coot::BALL_AND_STICK)
      bond_width = 0.15; // not 8

   graphics_info_t g;
   GtkWidget *dcw = g.display_control_window();

   // int imol = graphics_info_t::add_reps_molecule_option_menu_item_select_molecule;
   int imol = graphics_info_t::add_reps_molecule_combobox_molecule;

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(position_radiobutton))) {
      std::pair<bool, std::pair<int, coot::atom_spec_t> > aas = active_atom_spec();
      if (aas.first) {
	 int imol_active = aas.second.first;
	 coot::atom_selection_info_t asi(aas.second.second.chain_id,
					 aas.second.second.res_no,
					 aas.second.second.res_no,
					 aas.second.second.ins_code);
         GtkWidget *glarea_0 = 0;
         GtkWidget *glarea_1 = 0;
         if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
         if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
	 gl_context_info_t glci(glarea_0, glarea_1);
	 g.molecules[imol_active].add_additional_representation(representation_type,
								bonds_box_type,
								bond_width,
								draw_H_flag,
								asi, dcw, glci, g.Geom_p());
      }
   }
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(resno_radiobutton))) {
      // std::cout << "By chainid/resno" << std::endl;
      std::string chain_id = gtk_editable_get_text(GTK_EDITABLE(chain_id_entry));
      std::string resno_1s = gtk_editable_get_text(GTK_EDITABLE(resno_start_entry));
      std::string resno_2s = gtk_editable_get_text(GTK_EDITABLE(resno_end_entry));
      std::string ins_code = gtk_editable_get_text(GTK_EDITABLE(ins_code_entry));
      if (is_valid_model_molecule(imol)) {
	 int resno_1 = atoi(resno_1s.c_str());
	 int resno_2 = atoi(resno_2s.c_str());
	 coot::atom_selection_info_t asi(chain_id, resno_1, resno_2, ins_code);
         GtkWidget *glarea_0 = 0;
         GtkWidget *glarea_1 = 0;
         if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
         if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
	 gl_context_info_t glci(glarea_0, glarea_1);
	 graphics_info_t::molecules[imol].add_additional_representation(representation_type,
									bonds_box_type,
									bond_width,
									draw_H_flag,
									asi, dcw, glci, g.Geom_p());
      }
   }
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(selection_string_radiobutton))) {
      // std::cout << "By selection string" << std::endl;
      std::string s = gtk_editable_get_text(GTK_EDITABLE(string_selection_entry));
      coot::atom_selection_info_t asi(s);
      GtkWidget *glarea_0 = 0;
      GtkWidget *glarea_1 = 0;
      if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
      if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
      gl_context_info_t glci(glarea_0, glarea_1);
      graphics_info_t::molecules[imol].add_additional_representation(representation_type,
								     bonds_box_type,
								     bond_width,
								     draw_H_flag,
								     asi, dcw, glci, g.Geom_p());
   }
   graphics_draw();
}

#include "widget-from-builder.hh"

GtkWidget *wrapped_create_residue_editor_select_monomer_type_dialog() {

   // GtkWidget *w = create_residue_editor_select_monomer_type_dialog();
   GtkWidget *w = widget_from_builder("residue_editor_select_monomer_type_dialog");
   GtkWidget *combo_box = widget_from_builder("residue_editor_select_monomer_type_combobox");

   if (combo_box)
      gtk_combo_box_text_remove_all(GTK_COMBO_BOX_TEXT(combo_box));

   graphics_info_t g;
   std::vector<std::string> v = g.Geom_p()->monomer_types();

   // fill the combobox

   for (unsigned int i=0; i<v.size(); i++) {
      std::string s = v[i];
      gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT(combo_box), s.c_str());
      gtk_combo_box_set_active(GTK_COMBO_BOX(combo_box), i);
   }
   return w;
}


void show_restraints_editor_by_index(int menu_item_index) {

   graphics_info_t g;
   std::vector<std::string> v = g.Geom_p()->monomer_types();
   for (unsigned int i=0; i<v.size(); i++) {
      int i_int = i;
      if (i_int==menu_item_index)
	 show_restraints_editor(v[i_int].c_str());
   }
}


void clear_restraints_editor_by_dialog(GtkWidget *dialog) { /* close button pressed */
   graphics_info_t g;
   g.clear_restraints_editor_by_dialog(dialog);
}

void show_restraints_editor(std::string monomer_type) {

   int imol = 0; // maybe this should be passed? Pretty esoteric though.

   if (graphics_info_t::use_graphics_interface_flag) {

      if (false) {
	 std::cout << "ERROR:: null monomer_type - no restraints editor" << std::endl;
      } else {
	 graphics_info_t g;
	 coot::protein_geometry *pg = g.Geom_p();

	 std::pair<bool, coot::dictionary_residue_restraints_t> p =
	    pg->get_monomer_restraints(monomer_type, imol);

	 if (p.first) {
	    coot::dictionary_residue_restraints_t restraints = p.second;
	    coot::restraints_editor r;
	    r.fill_dialog(restraints);
            GtkWidget *dialog = r.get_dialog();
            std::cout << "DEBUG:: show_restraints_editor(): here with dialog " << dialog << std::endl;
	    set_transient_and_position(COOT_EDIT_RESTRAINTS_DIALOG, dialog);
	    g.restraints_editors.push_back(r);
	 }
      }
   }
}



// ===================================================================
//                   sequence view
// ===================================================================

void set_sequence_view_is_docked(short int state) {
   graphics_info_t::sequence_view_is_docked_flag = state;
}


void nsv(int imol) {

#ifdef HAVE_GOOCANVAS

   if (is_valid_model_molecule(imol)) {

      GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);

      if (w) {

	 // it already exists... just raise it and map it.

	 GtkWidget *canvas = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);

	 // so what is the window (which we shall call widget)?
	 // GtkWidget *widget = lookup_widget(canvas, "nsv_dialog");
	 GtkWidget *widget = widget_from_builder("nsv_dialog");

	 if (widget) {
	    // if (!GTK_WIDGET_MAPPED(widget)) { // gone
            if (true) {
	       gtk_widget_set_visible(widget, TRUE);
	    } else {
	       std::cout << "GTK-FIXME no raise" << std::endl;
	       // gdk_window_raise(widget->window);
	    }
	 } else {

            std::cout << "WARNING:: in nsv() lookup of nsv_dialog from canvas failed " << std::endl;

	    // widget = lookup_widget(canvas, "sequence_view_dialog");
	    widget = widget_from_builder("sequence_view_dialog");

	    if (widget) {
	       //if (!GTK_WIDGET_MAPPED(widget)) { // gone
               if (true) {
		  gtk_widget_set_visible(widget, TRUE);
	       } else {
#if (GTK_MAJOR_VERSION < 4)
                  GdkWindow *ww = gtk_widget_get_window(widget);
		  gdk_window_raise(ww);
#endif
	       }
	    }
	 }

      } else {

	 graphics_info_t g;
         GtkWidget *main_window_vbox = 0;
         if (g.sequence_view_is_docked_flag) {
            main_window_vbox = widget_from_builder("main_window_vbox");
         }

         std::cout << "::::::::::::::::::: debug:: sequence_view_is_docked_flag " << g.sequence_view_is_docked_flag
                   << " main_window_vbox " << main_window_vbox << std::endl;
	 std::string name = g.molecules[imol].name_for_display_manager();
	 exptl::nsv *seq_view =
	    new exptl::nsv(g.molecules[imol].atom_sel.mol, name, imol,
                           main_window_vbox,
			   g.use_graphics_interface_flag,
			   g.nsv_canvas_pixel_limit);
	 // I think that there is a false positive for scan-build here.
	 // The memory for sequence view is deleted before the new
	 // pointer is assigned.
	 g.set_sequence_view_is_displayed(seq_view->Canvas(), imol);
      }
   }
#endif // GOOCANVAS
}

void set_nsv_canvas_pixel_limit(int cpl) {
   graphics_info_t::nsv_canvas_pixel_limit = cpl;
}



void sequence_view_old_style(int imol) {

#if 0  // fix this another time - sequence_view needs to be converted to goocanvas
   graphics_info_t g;
   if (g.molecules[imol].has_model()) {
      graphics_info_t g;

      GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
      if (w) {

	 // it already exists... just raise it and map it.

	 // GtkWidget *canvas = g.sequence_view_is_displayed[imol];
	 GtkWidget *canvas = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
	 // so what is the window (which we shall call widget)?
	 // GtkWidget *widget = lookup_widget(canvas, "sequence_view_dialog");
	 GtkWidget *widget = widget_from_builder("sequence_view_dialog");

	 if (widget) {
            std::cout << "sequence_view_dialog() raise the widget here " << std::endl;
	    // if (!GTK_WIDGET_MAPPED(widget)) {
	    //    gtk_widget_set_visible(widget, TRUE);
	    // } else {
	    //    // gdk_window_raise(widget->window);
	    //    std::cout << "sequence_view_dialog() raise the widget here " << std::endl;
	    // }
	 }

      } else {

	 // create a new one

	 // Let's have a name that has the leading / and .pdb stripped
	 //
	 std::string short_name;

	 std::string::size_type islash = g.molecules[imol].name_.find_last_of("/");
	 std::string tstring;
	 if (islash == std::string::npos) {
	    // no slash found
	    tstring = g.molecules[imol].name_;
	 } else {
	    tstring = g.molecules[imol].name_.substr(islash + 1);
	 }

	 std::string::size_type ipdb = tstring.rfind(".pdb");

	 if (ipdb == std::string::npos) {
	    std::cout << "INFO .pdb not found in filename" << std::endl;
	    short_name = tstring;
	 } else {
	    short_name = tstring.substr(0, ipdb);
	 }
	 coot::sequence_view *seq_view =
	    new coot::sequence_view(g.molecules[imol].atom_sel.mol,
				    short_name, imol);
	 // I think that there is a false positive for scan-build here.
	 // The memory for sequence view is deleted before the new
	 // pointer is assigned.
	 g.set_sequence_view_is_displayed(seq_view->Canvas(), imol);
      }
   }
#endif
}


void set_sequence_view_is_displayed(GtkWidget *widget, int imol) {
   graphics_info_t g;
   g.set_sequence_view_is_displayed(widget, imol);
}

/*  ----------------------------------------------------------------------- */
/* Multirefine interface (because in guile-gtk there is no way to
		    insert toolbuttons into the toolbar) so this
		    rather kludgy interface.  It should go when we
		    move to guile-gnome, I think. */
/*  ----------------------------------------------------------------------- */
/* BL comment:: at the moment this has a python implementation too, it should
   however be using the toolbuttons on the main toolbar, as we can! Furthermore
   we have the general problem that in this way we cannot control a python
   started function if guile is available too. We would need to know if it was
   started from python or guile... */
/* well I have put the python buttons into the scripting layer using the main
   toolbar. It's still here for completeness. */
void toolbar_multi_refine_stop() {

#if (GTK_MAJOR_VERSION > 1)

#ifdef USE_GUILE

   // the idle function looks at this value
   std::string s = "(set! *continue-multi-refine* #f)";
   safe_scheme_command(s.c_str());

   set_visible_toolbar_multi_refine_continue_button(1);
   set_visible_toolbar_multi_refine_cancel_button(1);
   toolbar_multi_refine_button_set_sensitive("continue", 1);
   toolbar_multi_refine_button_set_sensitive("cancel",   1);
   toolbar_multi_refine_button_set_sensitive("stop",     0); // it's already stopped.

#else

#ifdef USE_PYTHON
   // the idle function looks at this value
   std::string s = "global continue_multi_refine; continue_multi_refine = False";
   safe_python_command(s.c_str());

   set_visible_toolbar_multi_refine_continue_button(1);
   set_visible_toolbar_multi_refine_cancel_button(1);
   toolbar_multi_refine_button_set_sensitive("continue", 1);
   toolbar_multi_refine_button_set_sensitive("cancel",   1);
   toolbar_multi_refine_button_set_sensitive("stop",     0); // it's already stopped.
#endif // USE_PYTHON
#endif // USE_GUILE
#endif // GTK_MAJOR_VERSION

}


void toolbar_multi_refine_continue() {

#if (GTK_MAJOR_VERSION > 1)

#ifdef USE_GUILE

   toolbar_multi_refine_button_set_sensitive("stop",     1);
   toolbar_multi_refine_button_set_sensitive("cancel",   0);
   toolbar_multi_refine_button_set_sensitive("continue", 0);
   std::string s = "(set! *continue-multi-refine* #t)";
   safe_scheme_command(s.c_str());
   s = "(gtk-idle-add *multi-refine-idle-proc*)";
   safe_scheme_command(s.c_str());

#else

#ifdef USE_PYTHON

   toolbar_multi_refine_button_set_sensitive("stop",     1);
   toolbar_multi_refine_button_set_sensitive("cancel",   0);
   toolbar_multi_refine_button_set_sensitive("continue", 0);
   std::string s = "global continue_multi_refine; continue_multi_refine = True";
   safe_python_command(s.c_str());
   s = "global multi_refine_idle_proc; gobject.idle_add(multi_refine_idle_proc)";
   safe_python_command(s.c_str());

#endif // USE_PYTHON
#endif // USE_GUILE
#endif // GTK_MAJOR_VERSION

}

void toolbar_multi_refine_cancel() {

#if (GTK_MAJOR_VERSION > 1)
#ifdef USE_GUILE

   // the idle function looks at this value
   std::string s = "(set! *continue-multi-refine* #f)";
   safe_scheme_command(s.c_str());
   toolbar_multi_refine_button_set_sensitive("stop", 1); // for next time
   set_visible_toolbar_multi_refine_continue_button(0);
   set_visible_toolbar_multi_refine_stop_button(0);
   set_visible_toolbar_multi_refine_cancel_button(0);

#else

#ifdef USE_PYTHON

   // the idle function looks at this value
   std::string s = "global continue_multi_refine; continue_multi_refine = False";
   safe_python_command(s.c_str());
   toolbar_multi_refine_button_set_sensitive("stop", 1); // for next time
   set_visible_toolbar_multi_refine_continue_button(0);
   set_visible_toolbar_multi_refine_stop_button(0);
   set_visible_toolbar_multi_refine_cancel_button(0);

#endif // USE_PYTHON
#endif // USE_GUILE
#endif // GTK_MAJOR_VERSION

}



void set_visible_toolbar_multi_refine_stop_button(short int state) {

   graphics_info_t g;
   if (graphics_info_t::use_graphics_interface_flag) {
      // GtkWidget *w = lookup_widget(g.get_main_window(), "toolbar_multi_refine_stop_button");
      GtkWidget *w = widget_from_builder("toolbar_multi_refine_stop_button");
      if (w) {
	 if (state) {
	    gtk_widget_set_visible(w, TRUE);
	 } else {
	    gtk_widget_set_visible(w, FALSE);
	 }
      }
   }
}

void set_visible_toolbar_multi_refine_continue_button(short int state) {

   graphics_info_t g;
   if (graphics_info_t::use_graphics_interface_flag) {
      // GtkWidget *w = lookup_widget(g.get_main_window(), "toolbar_multi_refine_continue_button");
      GtkWidget *w = widget_from_builder("toolbar_multi_refine_continue_button");
      if (w) {
	 if (state) {
	    gtk_widget_set_visible(w, TRUE);
	 } else {
	    gtk_widget_set_visible(w, FALSE);
	 }
      }
      toolbar_multi_refine_button_set_sensitive("cancel", 0);
   }
}

void set_visible_toolbar_multi_refine_cancel_button(short int state) {

   graphics_info_t g;
   if (graphics_info_t::use_graphics_interface_flag) {
      // GtkWidget *w = lookup_widget(g.get_main_window(), "toolbar_multi_refine_cancel_button");
      GtkWidget *w = widget_from_builder("toolbar_multi_refine_cancel_button");
      if (w) {
	 if (state) {
	    gtk_widget_set_visible(w, TRUE);
	 } else {
	    gtk_widget_set_visible(w, FALSE);
	 }
      }
   }
}

/* button_type is one of "stop", "continue", "cancel"
   state is 1 for on, 0 for off. */
void toolbar_multi_refine_button_set_sensitive(const char *button_type, short int state) {
   GtkWidget *w = NULL;

   std::string bt(button_type);

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      if (bt == "cancel")
	 w = widget_from_builder("toolbar_multi_refine_cancel_button");
      if (bt == "continue")
	 w = widget_from_builder("toolbar_multi_refine_continue_button");
      if (bt == "stop")
	 w = widget_from_builder("toolbar_multi_refine_stop_button");

      if (w) {
	 if (state) {
	    gtk_widget_set_sensitive(w, TRUE);
	 } else {
	    gtk_widget_set_sensitive(w, FALSE);
	 }
      }
   }
}




// ---------------------------------------------
//        Map Sharpening dialog
// ---------------------------------------------

GtkWidget *wrapped_create_map_sharpening_dialog() {

   std::cout << ":::::::::::::::::::::: wrapped_create_map_sharpening_dialog()" << std::endl;

   float sharpening_limit = graphics_info_t::map_sharpening_scale_limit;

   // GtkWidget *w = create_map_sharpening_dialog();

   GtkWidget *w = widget_from_builder("map_sharpening_dialog");

   graphics_info_t g;
   GCallback signal_func = G_CALLBACK(map_sharpening_map_select_combobox_changed);
   GtkWidget *combobox = widget_from_builder("map_sharpening_molecule_combobox");

   int imol_prefered = imol_refinement_map();
   int imol = g.fill_combobox_with_map_mtz_options(combobox, signal_func, imol_prefered); // map options now

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::imol_map_sharpening = imol;

      GtkWidget *h_scale = widget_from_builder("map_sharpening_hscale");

      GtkAdjustment *adj = gtk_adjustment_new(0.0, -sharpening_limit, 2*sharpening_limit,
					      0.05, 0.2, (sharpening_limit+0.1));
      gtk_range_set_adjustment(GTK_RANGE(h_scale), GTK_ADJUSTMENT(adj));
      g_object_set_data_full(G_OBJECT (w), "map_sharpening_adjustment",
			     g_object_ref (adj),
			     (GDestroyNotify) g_object_unref);

      g_signal_connect(G_OBJECT(adj), "value_changed", G_CALLBACK(map_sharpening_value_changed), NULL);

      // set to sharpening value
      gtk_adjustment_set_value(GTK_ADJUSTMENT(adj), g.molecules[imol].sharpen_b_factor());

      int ticks = 3;  // number of ticks on the (one) side (not including centre tick)
      for (int i=0; i<=2*ticks; i++) {
	 float p = float (i-ticks) * (1.0/float(ticks)) * sharpening_limit;
	 std::string pos_string = coot::util::float_to_string_using_dec_pl(p,0);
	 gtk_scale_add_mark(GTK_SCALE(h_scale), p, GTK_POS_BOTTOM, pos_string.c_str());
      }
      gtk_scale_add_mark(GTK_SCALE(h_scale), -sharpening_limit, GTK_POS_BOTTOM, "\n  Sharpen");
      gtk_scale_add_mark(GTK_SCALE(h_scale),  sharpening_limit, GTK_POS_BOTTOM, "\nBlur");

   }

   return w;
}

void
calc_and_set_optimal_b_factor ( GtkWidget *w ) {

   float sharpening_limit = graphics_info_t::map_sharpening_scale_limit;
   int imol = graphics_info_t::imol_map_sharpening;
   float Bopt = optimal_B_kurtosis(imol);
   if (fabs(Bopt-graphics_info_t::map_sharpening_scale_limit) <= 0.1) {
      std::string txt;
      txt = "INFO:: Optimisation did NOT converge.\n The value may be bogus.";
      info_dialog_and_text(txt.c_str());
   }
   // GtkWidget *h_scale = lookup_widget(w, "map_sharpening_hscale");
   GtkWidget *h_scale = widget_from_builder("map_sharpening_hscale");
   GtkAdjustment *adj = gtk_range_get_adjustment(GTK_RANGE(h_scale));
   gtk_adjustment_set_value(adj, Bopt);
}

void map_sharpening_map_select_combobox_changed(GtkWidget *widget, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(widget));
   graphics_info_t::imol_map_sharpening = imol;

}


void
map_sharpening_map_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t::imol_map_sharpening = pos;

   GtkWidget *adj = widget_from_builder("map_sharpening_adjustment");
   gtk_adjustment_set_value(GTK_ADJUSTMENT(adj), graphics_info_t::molecules[pos].sharpen_b_factor());

}

void map_sharpening_value_changed (GtkAdjustment *adj,
				   GtkWidget *window) {

   int imol = graphics_info_t::imol_map_sharpening;
   float value = gtk_adjustment_get_value(adj);
   // std::cout << "sharpen " << imol << " by " << value << std::endl;
   if (is_valid_map_molecule(imol)) {
      sharpen(imol, value);
   }
}

void set_baton_build_params_from_widget(GtkWidget *params_dialog) {

   GtkWidget *ent_res_no   = widget_from_builder("baton_build_params_residue_number_entry");
   GtkWidget *ent_chain_id = widget_from_builder("baton_build_params_chain_id_entry");
   GtkWidget *check_button = widget_from_builder("baton_build_params_backwards_checkbutton");

   const char *resno_txt = gtk_editable_get_text(GTK_EDITABLE(ent_res_no));
   const char *chain_id  = gtk_editable_get_text(GTK_EDITABLE(ent_chain_id));

   const char *direction = "forwards";
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_button)))
      direction = "backwards";

   try {
      int rn = coot::util::string_to_int(resno_txt);
      set_baton_build_params(rn, chain_id, direction);
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << " aborting." << std::endl;
   }

}

// keyboarding mode
void show_go_to_residue_keyboarding_mode_window() {

   // GtkWidget *w = create_keyboard_goto_residue_window();
   GtkWidget *w = widget_from_builder("keyboard_go_to_residue_window");
   graphics_info_t g;
   // g.graphics_x_position, graphics_x_size
   int x_pos = g.graphics_x_position + 5;
   int y_pos = g.graphics_y_position + g.graphics_y_size + 65;

   set_transient_and_position(COOT_UNDEFINED_WINDOW, w);
   gtk_widget_set_visible(w, TRUE);

}


void handle_go_to_residue_keyboarding_mode(const char *text) {
   graphics_info_t::apply_go_to_residue_keyboading_string(text);
}


void
on_instanced_mesh_generic_objects_dialog_object_check_button_toggled(GtkCheckButton *button,
                                                                     gpointer user_data) {

   int combo_ints = GPOINTER_TO_INT(user_data);
   int imol = combo_ints/1000;
   int obj_no = combo_ints - 1000 * imol;
   bool state = false;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(button)))
      state = 1;

   std::cout << "debug:: on_instanced_mesh_generic_objects_dialog_object_check_button_toggled() imol " << imol
             << " obj_no " << obj_no << std::endl;
   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      int n_meshes = m.instanced_meshes.size();
      if (obj_no >=0 && obj_no < n_meshes) {
         m.instanced_meshes[obj_no].set_draw_status(state);
         graphics_draw();
      }
   }
}


void
generic_objects_dialog_grid_add_object_for_molecule_internal(int imol,
                                                             int mesh_index,
                                                             int grid_row_offset, // because there are already items in the grid
                                                             const Instanced_Markup_Mesh &imm,
                                                             GtkWidget *dialog,
                                                             GtkWidget *grid) {

   if (! imm.is_closed()) {
      GtkWidget *checkbutton = gtk_check_button_new_with_mnemonic (("Display"));
      std::string label_str = imm.get_name();
      // label_str = "NMO: "  + label_str;
      GtkWidget *label = gtk_label_new(label_str.c_str());

      int i_row = grid_row_offset;

      std::string stub = "generic_object_" + std::to_string(i_row);
      std::string toggle_button_name = stub + "_toggle_button"; // is it a toggle button?
      std::string label_name = stub + "_label";

      // set the names of these widgets so that they can be
      // looked up and toggled/hidden dynamically.

      g_object_set_data(G_OBJECT(dialog), toggle_button_name.c_str(), checkbutton);
      g_object_set_data(G_OBJECT(dialog), label_name.c_str(), label);

      std::cout << "debug:: imm with name " << label_str << " at row " << i_row << std::endl;

      // grid child left top width height
      gtk_grid_attach (GTK_GRID (grid), label,       0, i_row, 1, 1);
      gtk_grid_attach (GTK_GRID (grid), checkbutton, 1, i_row, 1, 1);

      if (imm.get_draw_status())
	 gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);

      g_signal_connect(G_OBJECT(checkbutton), "toggled",
		       G_CALLBACK(on_instanced_mesh_generic_objects_dialog_object_check_button_toggled),
		       GINT_TO_POINTER(imol * 1000 + mesh_index));

      gtk_widget_set_visible (label, TRUE);
      gtk_widget_set_visible (checkbutton, TRUE);

   }
}




GtkWidget *wrapped_create_generic_objects_dialog() {

   auto clear_the_grid = [] (GtkWidget *grid) {
      GtkWidget *child = gtk_widget_get_first_child(GTK_WIDGET(grid));
      unsigned int child_count = 0;
      while (child) {
         child_count += 1;
         GtkWidget *next = gtk_widget_get_next_sibling(child);
         // gtk_grid_remove_row(GTK_GRID(grid), 0);
         gtk_grid_remove(GTK_GRID(grid), child);
         child = next; // for next round
      };
      // std::cout << "There were " << child_count << " children in the grid" << std::endl;
   };

   graphics_info_t g;

   GtkWidget *dialog = widget_from_builder("generic_objects_dialog");
   g.generic_objects_dialog = dialog;

   GtkWidget *grid = widget_from_builder("generic_objects_dialog_grid");
   if (grid) {

      clear_the_grid(grid);
      unsigned int io_count = 0;
      unsigned int n_objs = g.generic_display_objects.size();
      for (unsigned int io=0; io<n_objs; io++) {
	 const meshed_generic_display_object &gdo = g.generic_display_objects.at(io);
         if (! gdo.mesh.is_closed()) {
            g.generic_objects_dialog_grid_add_object_internal(gdo, dialog, grid, io);
            io_count++;
         }
      }

      for (unsigned int imol=0; imol<g.molecules.size(); imol++) {
         const molecule_class_info_t &m = g.molecules[imol];
         for (unsigned int j=0; j<m.instanced_meshes.size(); j++) {
            const Instanced_Markup_Mesh &imm = m.instanced_meshes[j];
            if (! imm.is_closed()) {
               generic_objects_dialog_grid_add_object_for_molecule_internal(imol, j, io_count, imm, dialog, grid);
               io_count++;
            }
         }
      }
   } else {
      std::cout << "failed to get grid " << std::endl;
   }

   return dialog;
}


/* return a new object number (so that we can set it to be displayed). */
int add_generic_display_object(const meshed_generic_display_object &gdo) {

   graphics_info_t g;
   int n_objs = g.generic_display_objects.size();
   g.generic_display_objects.push_back(gdo);
   if (g.generic_objects_dialog) {
      // GtkWidget *table = lookup_widget(g.generic_objects_dialog, "generic_objects_dialog_grid");
      GtkWidget *table = widget_from_builder("generic_objects_dialog_grid");
      if (table) {
         // auto resize now
	 // gtk_table_resize(GTK_TABLE(table), n_objs+1, 2);
	 g.generic_objects_dialog_grid_add_object_internal(gdo,
                                                           g.generic_objects_dialog,
                                                           table,
                                                           n_objs);
      }
   }
   return n_objs;
}


// -------------- You can't install extensions if you don't have boost ----------

#ifdef HAVE_BOOST
#include <boost/crc.hpp>
#endif

std::pair<bool, std::string>
checksums_match(const std::string &file_name, const std::string &checksum) {

   bool state = false;
   std::string message;

#ifdef HAVE_BOOST
   std::ifstream f(file_name.c_str());
   if (f) {
      std::string dl_str((std::istreambuf_iterator<char>(f)),
			 std::istreambuf_iterator<char>());

      // boost::crc_basic<16> crc_ccitt1( 0x1021, 0xFFFF, 0, false, false );
      boost::crc_basic<16> crc_ccitt1(0xffff, 0x0, 0, false, false );
      crc_ccitt1.process_bytes(dl_str.c_str(), dl_str.size());
      // std::cout << "DEBUG:: checksum compare " << crc_ccitt1.checksum() << " " << checksum << std::endl;
      std::string s = coot::util::int_to_string(crc_ccitt1.checksum());
      if (s == checksum)
	 state = true;
      else
         message = s + " vs " + checksum;
   }
#endif // HAVE_BOOST
   return std::pair<bool, std::string> (state, message);
}

#include "cc-interface.hh"

// 20200302 Be'er Sheva new style extension installation
// Put these in c-interface-curlew?
void
curlew_install_extension_file(const std::string &file_name, const std::string &checksum,
                              GtkWidget *install_button, GtkWidget *uninstall_button) {

   if (!file_name.empty()) {

#ifndef WINDOWS_MINGW
      // std::string url_prefix = "https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/";
      std::string url_prefix = "https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/curlew-extensions/Coot-1/";
#else
       std::string url_prefix = "https://bernhardcl.github.io/coot/";
#endif
      url_prefix += "scripts";
      url_prefix += "/";
      url_prefix += file_name;

      std::string download_dir = "coot-download";
      download_dir = coot::get_directory(download_dir.c_str());
      std::string dl_fn = download_dir + "/";
      dl_fn += file_name;

      if (false) // debug
         std::cout << "get this " << url_prefix << " as this " << dl_fn << std::endl;

      int r = coot_get_url(url_prefix.c_str(), dl_fn.c_str());

      if (r) {
         std::cout << "WARNING:: bad URL retrieve " << file_name << std::endl;
      } else {
         // Happy path
         if (coot::file_exists(dl_fn)) {
            std::pair<bool, std::string> checksum_result = checksums_match(dl_fn, checksum);
            if (checksum_result.first) {
               // I want a function that returns preferences_dir
               std::string home_directory = coot::get_home_dir();
               if (!home_directory.empty()) {
                  std::string preferences_dir = coot::util::append_dir_dir(home_directory, ".coot");
                  std::string preferences_file_name = coot::util::append_dir_file(preferences_dir, file_name);
                  std::cout << "debug:: attempting to copy " << dl_fn << " as " << preferences_file_name << std::endl;
                  int status = coot::copy_file(dl_fn, preferences_file_name); // it returns a bool actually
                  if (status != 0) {
                     std::cout << "WARNING:: rename status " << status << " failed to install " << file_name << std::endl;
                     std::cout << "WARNING:: rename error: " << strerror(errno) << std::endl;
                     std::cout << "WARNING:: fall-back: run the script from download-dir: " << dl_fn << std::endl;
                     run_script(dl_fn.c_str());
                  } else {
                     std::cout << "debug:: renaming successful" << std::endl;
                     std::cout << "debug:: run_script() called on " << preferences_file_name << std::endl;
                     run_script(preferences_file_name.c_str());
                     //  std::cout << "hiding install_button " << install_button << std::endl;
                     gtk_widget_set_visible(install_button, FALSE);
                     //  std::cout << "show uninstall_button  " << uninstall_button << std::endl;
                     gtk_widget_set_visible(uninstall_button, TRUE);
                     
                  }
               } else {
                  std::cout << "No HOME env var" << std::endl;
               }
            } else {
               std::cout << "WARNING:: Failure in checksum match " << dl_fn << " " << checksum_result.second << std::endl;
            }
         } else {
            std::cout << "WARNING:: download target file " << dl_fn << " does not exist" << std::endl;
         }
      }
   }
}

bool
curlew_uninstall_extension_file(const std::string &file_name) {

   bool r_status = false;

   // I want a function that returns preferences_dir
   std::string home = coot::get_home_dir();
   if (!home.empty()) {
      std::string home_directory(home);
      std::string preferences_dir = coot::util::append_dir_dir(home_directory, ".coot");
      std::string preferences_file_name = coot::util::append_dir_file(preferences_dir, file_name);
      std::string renamed_file_name = preferences_file_name + "_uninstalled";
      if (coot::file_exists(preferences_file_name)) {
#ifndef WINDOWS_MINGW
         int status = rename(preferences_file_name.c_str(), renamed_file_name.c_str());
#else
          int status = coot::rename_win(preferences_file_name.c_str(), renamed_file_name.c_str());
#endif
         if (status != 0) {
            std::cout << "WARNING:: rename status " << status << " failed to uninstall " << file_name << std::endl;
         } else {
            // OK
            r_status = true;
         }
      }
   }
   return r_status;
}


// If I put this in c-interface-curlew.cc, then it doesn't resolve when linking.
// I don't understand why. It's called from callbacks.c
//
/* curlew install button clicked callback action */
void curlew_dialog_install_extensions(GtkWidget *curlew_dialog, int n_extensions) {

   // Look up the checkbuttow widgets, find the attached filenames and checksums
   // download them, check the checksums and, if they match install them to preferences
   // and run them.

   if (curlew_dialog) {
      for (int i=0; i<n_extensions; i++) {
	 std::string cb_name = "curlew_selected_check_button_";
	 cb_name += coot::util::int_to_string(i);
	 std::string uninstall_button_name = "curlew_uninstall_button_";
	 uninstall_button_name += coot::util::int_to_string(i);
	 std::string hbox_name = "curlew_extension_hbox_";
	 hbox_name += coot::util::int_to_string(i);

	 // GtkWidget *check_button     = lookup_widget(curlew_dialog, cb_name.c_str());
	 // GtkWidget *uninstall_button = lookup_widget(curlew_dialog, uninstall_button_name.c_str());
	 // GtkWidget *hbox             = lookup_widget(curlew_dialog, hbox_name.c_str());

         // I am not sure about this! Curlew doesn't work yet anyway! (scripts not converted or checked)
	 GtkWidget *check_button     = widget_from_builder(cb_name.c_str());
	 GtkWidget *uninstall_button = widget_from_builder(uninstall_button_name.c_str());
	 GtkWidget *hbox             = widget_from_builder(hbox_name.c_str());

	 if (check_button) {
	    int status = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_button));
	    if (status) { // selected for download/install
	       if (true)
		  std::cout << "Got check_button " << check_button << " for i " << cb_name << " " << status
		            <<std::endl;

               gchar *file_name_cstr = static_cast<gchar *> (g_object_get_data(G_OBJECT(check_button), "file-name"));
               gchar *checksum_cstr  = static_cast<gchar *> (g_object_get_data(G_OBJECT(check_button), "checksum"));

               if (file_name_cstr) {

		  std::string file_name(file_name_cstr);

		  if (!file_name.empty()) {

#ifndef WINDOWS_MINGW
              std::string url_prefix = "https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/";
#else
              std::string url_prefix = "https://bernhardcl.github.io/coot/";
#endif
		     url_prefix += "extensions";
		     url_prefix += "/";
		     url_prefix += file_name;

		     std::string download_dir = "coot-download"; // FIXME
		     std::string dl_fn = download_dir + "/";
		     dl_fn += file_name;

		     if (false) // debug
			std::cout << "get this " << url_prefix << " as this " << dl_fn << std::endl;

		     int r = coot_get_url(url_prefix.c_str(), dl_fn.c_str());

		     if (r) {
			std::cout << "WARNING:: bad URL retrieve " << file_name << std::endl;
		     } else {
			if (coot::file_exists(dl_fn)) {
			   std::string checksum;
			   if (checksum_cstr) checksum = checksum_cstr;
                           std::pair<bool, std::string> checksum_result = checksums_match(dl_fn, checksum);
                           if (checksum_result.first) {
			      // I want a function that returns preferences_dir
                              std::string home_directory = coot::get_home_dir();
                              if (!home_directory.empty()) {
				 std::string preferences_dir = coot::util::append_dir_dir(home_directory, ".coot");
				 std::string preferences_file_name = coot::util::append_dir_file(preferences_dir, file_name);
                                 std::cout << "debug:: attempting to rename " << dl_fn << " as " << preferences_file_name << std::endl;
#ifndef WINDOWS_MINGW
				 int status = rename(dl_fn.c_str(), preferences_file_name.c_str());
#else
                 int status = coot::rename_win(dl_fn.c_str(), preferences_file_name.c_str());
#endif
				 if (status != 0) {
				    std::cout << "WARNING:: rename status " << status << " failed to install " << file_name << std::endl;
				 } else {
                                    std::cout << "debug:: AA  renaming successful" << std::endl;
				    std::cout << "debug:: AA run_script() on " << preferences_file_name << std::endl;
				    run_script(preferences_file_name.c_str());

                                    std::cout << "hiding check_button " << check_button << std::endl;
                                    gtk_widget_set_visible(check_button, FALSE);
                                    std::cout << "show uninstall_button  " << uninstall_button << std::endl;
                                    gtk_widget_set_visible(uninstall_button, TRUE);
				    // make the hbox insensitive
				    if (hbox) {
				       gtk_widget_set_sensitive(hbox, FALSE);
				    }
				 }
			      } else {
                                 std::cout << "No HOME env var" << std::endl;
                              }
			   } else {
			      std::cout << "WARNING:: Failure in checksum match " << dl_fn << " " << checksum_result.second << std::endl;
			   }
			} else {
			   std::cout << "WARNING:: file does not exist " << file_name << std::endl;
			}
		     }
		  } else {
		     std::cout << "WARNING:: file_name data was empty" << std::endl;
		  }
	       } else {
		  std::cout << "WARNING:: No file name data" << std::endl;
	       }
	    }
	 }
      }
   }
}

