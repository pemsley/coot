/* src/c-interface-gui.cc
 * 
 * Copyright 2003, 2004, 2007 The University of York
 * Copyright 2008 The University of Oxford
 * Author: Paul Emsley
 * Copyright 2007 The University of York
 * Author: Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

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
#include "rotate-translate-modes.hh"
#include "nsv.hh"

#include "graphics-info.h"
#include "interface.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "cmtz-interface.hh"
#include "mmdb.h"  // for centre of molecule
#include "clipper/core/clipper_instance.h"

#if (GTK_MAJOR_VERSION > 2 || (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION > 5))
#define HAVE_GTK_COMBO_BOX_GET_ACTIVE_TEXT
#endif

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


/*! \brief display the open coordinates dialog */
void open_coords_dialog() {


   if (graphics_info_t::use_graphics_interface_flag) { 
   
#if (GTK_MAJOR_VERSION == 1)

      GtkWidget *coords_fileselection1 = coot_file_chooser();
      GtkWidget *file_filter_button;
      GtkWidget *sort_button;
      add_ccp4i_project_optionmenu(coords_fileselection1, COOT_COORDS_FILE_SELECTION);

      file_filter_button = add_filename_filter_button(coords_fileselection1, 
						      COOT_COORDS_FILE_SELECTION);
      sort_button = add_sort_button_fileselection(coords_fileselection1);
      add_recentre_on_read_pdb_checkbutton(coords_fileselection1);
      set_directory_for_fileselection(coords_fileselection1);
      push_the_buttons_on_fileselection(file_filter_button, sort_button, 
					coords_fileselection1);
      set_file_selection_dialog_size(coords_fileselection1);
      gtk_widget_show (coords_fileselection1);

#else

      /* This split was here because the buttons don't work. They act on the
	 file list, using the file list as a CList.  And CList is deprecated
	 in GTk+2.  So the button-press callback code needs to be adjusted. */
      GtkWidget *file_filter_button;
      GtkWidget *sort_button;
      GtkWidget *coords_fileselection1 = coot_file_chooser(); // a chooser or selector, depends.
      add_ccp4i_project_optionmenu(coords_fileselection1, COOT_COORDS_FILE_SELECTION);
      file_filter_button = add_filename_filter_button(coords_fileselection1, 
						      COOT_COORDS_FILE_SELECTION);
      sort_button = add_sort_button_fileselection(coords_fileselection1); 
      add_recentre_on_read_pdb_checkbutton(coords_fileselection1);
      set_directory_for_coot_file_chooser(coords_fileselection1);
      set_file_selection_dialog_size(coords_fileselection1);
      gtk_widget_show (coords_fileselection1);
      /* in gtk2 we have to push the buttons after we show the selection */
      push_the_buttons_on_fileselection(file_filter_button, sort_button, 
					coords_fileselection1);
#endif
   }
}
  


void
set_graphics_rotamer_dialog(GtkWidget *w) {
   graphics_info_t::rotamer_dialog = w;
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
   gtk_widget_show(w);

}

GtkWidget *wrapped_create_undo_molecule_chooser_dialog() {

   GtkWidget *w = create_undo_molecule_chooser_dialog();
   GtkWidget *option_menu = lookup_widget(w, "undo_molecule_chooser_option_menu");
   graphics_info_t g;
   
   g.fill_option_menu_with_undo_options(option_menu);
   return w;
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
      
      if (w) 
	 gtk_widget_show(w);
   }
   std::string cmd = "manage-column-selector";
   std::vector<coot::command_arg_t> args;
   args.push_back(filename);
   add_to_history_typed(cmd, args);
}

void
manage_refmac_column_selection(GtkWidget *w) {

   if (graphics_info_t::use_graphics_interface_flag) {
     coot::setup_refmac_parameters_from_file(w);
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


void handle_column_label_make_fourier(GtkWidget *column_label_window) {

  GtkWidget *refmac_checkbutton;
   int icol; 
   int use_weights = 0;
   int is_diff_map;
   short int sensible_r_free_col = 0;
   short int have_refmac_params = 0; /* default not */
   short int use_resolution_limits_flag = 0;
   float low_reso_lim = -1.0;	/* unset */
   float high_reso_lim = -1.0;	/* unset */
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

  std::string phi_label;
  std::string f_label;
  std::string w_label;
  std::string fobs_col;
  std::string sigfobs_col;
  std::string r_free_col;

  /* Was the "Use Weights checkbutton clicked?  */

   check_weights = GTK_CHECK_BUTTON(lookup_widget(GTK_WIDGET(column_label_window), 
					    "use_weights_checkbutton"));

   
   if (GTK_TOGGLE_BUTTON(check_weights)->active == 1) { 
      use_weights = 1;
   } else { 
      use_weights = 0;
   }

 /* Similarly, we ask, was the "Is difference map" checkbutton active?  */

   is_diff_map_checkbutton = GTK_CHECK_BUTTON(lookup_widget(GTK_WIDGET(column_label_window),
							    "difference_map_checkbutton"));

   if (GTK_TOGGLE_BUTTON(is_diff_map_checkbutton)->active) { 
     is_diff_map = 1;
   } else{ 
     is_diff_map = 0;
   }
     
   coot::mtz_column_types_info_t *saved_f_phi_columns
      = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(column_label_window));

   const char *object_mtz_filename = saved_f_phi_columns->mtz_filename.c_str();
   /* Get the values that the user has selected in the option menu
      buttons. */ 

   icol = saved_f_phi_columns->selected_phi_col; 
   if (icol == -1) { 
     printf("WARNING!!!! No phi col was set!!!!!!! \n");
   } else { 
      phi_label = saved_f_phi_columns->phi_cols[icol].column_label; 

     icol = saved_f_phi_columns->selected_f_col;
     if (icol < int(saved_f_phi_columns->f_cols.size()))
	f_label = saved_f_phi_columns->f_cols[icol].column_label;
     else {
	f_label = saved_f_phi_columns->d_cols[icol-saved_f_phi_columns->f_cols.size()].column_label;
	is_anomalous_flag = 1;
     }
   
     if (use_weights) { 
	std::cout << " Making map from " << f_label << " " << phi_label << " and "
		  << w_label << std::endl;
	icol = saved_f_phi_columns->selected_weight_col;
	w_label = saved_f_phi_columns->weight_cols[icol].column_label;
     } else { 
	std::cout << " Making map from " << f_label << " and " << phi_label << std::endl;
     } 

     /* is the resolution limit check button in use? */
     resolution_limit_check_button = 
       GTK_CHECK_BUTTON(lookup_widget(GTK_WIDGET(column_label_window), 
				       "column_labels_use_resolution_limits_checkbutton"));
     if (GTK_TOGGLE_BUTTON(resolution_limit_check_button)->active) { 

       /* yes, it is.. */

       low_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(column_label_window),
			     "column_labels_reso_low_entry"));
       high_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(column_label_window),
			      "column_labels_reso_high_entry"));

       low_reso_lim  = get_positive_float_from_entry(low_entry);
       high_reso_lim = get_positive_float_from_entry(high_entry);
       if (low_reso_lim > 0.001) 
	 if (high_reso_lim > 0.0001)
	   use_resolution_limits_flag = 1;
     }

     /* Refmac label stuff */

     refmac_checkbutton = lookup_widget(GTK_WIDGET(column_label_window),
					"refmac_column_labels_checkbutton");

     if (GTK_TOGGLE_BUTTON(refmac_checkbutton)->active) { 

       have_refmac_params = 1; 

       /* find the refmac option menu */
       fobs_option_menu    = lookup_widget(column_label_window, "refmac_fobs_optionmenu");
       sigfobs_option_menu = lookup_widget(column_label_window, "refmac_sigfobs_optionmenu");
       r_free_option_menu  = lookup_widget(column_label_window, "refmac_rfree_optionmenu");
  
       /* find the refmac menus */
       fobs_menu    = gtk_option_menu_get_menu(GTK_OPTION_MENU(fobs_option_menu));
       sigfobs_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(sigfobs_option_menu));
       r_free_menu  = gtk_option_menu_get_menu(GTK_OPTION_MENU(r_free_option_menu));

       /* now assign the columns */
       icol = saved_f_phi_columns->selected_refmac_fobs_col;
       fobs_col = saved_f_phi_columns->f_cols[icol].column_label;

       icol = saved_f_phi_columns->selected_refmac_sigfobs_col;
       sigfobs_col = saved_f_phi_columns->sigf_cols[icol].column_label;

       icol = saved_f_phi_columns->selected_refmac_r_free_col; /* magic -1 if not set */
       if (icol >= 0) { 
	 // 
	 sensible_r_free_col = 1;
	 r_free_col = saved_f_phi_columns->r_free_cols[icol].column_label;
       } else { 
	 sensible_r_free_col = 0;
	 r_free_col = "";
       }
     }


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
   gtk_widget_destroy(column_label_window);

}



void fill_f_optionmenu_with_expert_options(GtkWidget *f_optionmenu) {

   coot::fill_f_optionmenu(f_optionmenu, 1);

}

void
wrapped_create_run_refmac_dialog() {
  GtkWidget *window = create_run_refmac_dialog();
  GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(refmac_molecule_button_select);
  GtkWidget *diff_map_button = lookup_widget(window, "run_refmac_diff_map_checkbutton");
  GtkWidget *optionmenu;
  int imol_coords = first_coords_imol();
  int have_file = 0;
  /* only GTK 2!? */
#if (GTK_MAJOR_VERSION > 1)
  GtkWidget *labels = lookup_widget(window, "run_refmac_column_labels_frame");
  GtkWidget *ncs_button = lookup_widget(window, "run_refmac_ncs_checkbutton");
  GtkWidget *mtz_file_radiobutton = lookup_widget(window, "run_refmac_mtz_file_radiobutton");
  optionmenu = lookup_widget(window, "run_refmac_method_optionmenu");
  fill_option_menu_with_refmac_methods_options(optionmenu);

  optionmenu = lookup_widget(window, "run_refmac_phase_input_optionmenu");
  fill_option_menu_with_refmac_phase_input_options(optionmenu);
  if (GTK_TOGGLE_BUTTON(mtz_file_radiobutton)->active) have_file = 1;
#endif

  set_refmac_molecule(imol_coords);

  optionmenu = lookup_widget(window, "run_refmac_coords_optionmenu");
  fill_option_menu_with_coordinates_options(optionmenu, callback_func, imol_coords);

  optionmenu = lookup_widget(window, "run_refmac_map_optionmenu");
  /*  fill_option_menu_with_refmac_options(optionmenu); */
  fill_option_menu_with_refmac_labels_options(optionmenu);

  /* to set the labels set the active item; only if not twin and
     if we really want the labels from map mtz*/
  if (refmac_use_twin_state() == 0 && have_file == 0) {
    GtkWidget *active_menu_item = gtk_menu_get_active(GTK_MENU(gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu))));
    if (active_menu_item) {
      gtk_menu_item_activate(GTK_MENU_ITEM(active_menu_item));
    }
  }

  /* fill optionmenu for no label refmac and show if refmac version is new enough */
  /* only GTK 2!? */
#if (GTK_MAJOR_VERSION > 1)
  if (refmac_runs_with_nolabels()) {
    GtkWidget *checkbutton = lookup_widget(window, "run_refmac_nolabels_checkbutton");
    gtk_widget_show(checkbutton);
    if (get_refmac_phase_input()) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), FALSE);
      gtk_widget_show(labels);
    } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
      gtk_widget_hide(labels);
    }
    GtkWidget *extra_options = lookup_widget(window, "run_refmac_extra_refinement_options_frame");
    GtkWidget *twin_check_button = lookup_widget(window, "run_refmac_twin_checkbutton");
    GtkWidget *sad_extras = lookup_widget(window, "run_refmac_sad_extra_hbox");
    gtk_widget_show(extra_options);
    gtk_widget_hide(twin_check_button);

    if (refmac_runs_with_nolabels() >= 2) {
      /* add the tls, twin and sad buttons */
       gtk_widget_show(twin_check_button);
      /* update the check buttons */
      GtkWidget *mtz_file_label = lookup_widget(window, "run_refmac_mtz_file_label");
      store_refmac_mtz_file_label(mtz_file_label);
      /* set the filename if there */
      const gchar *mtz_filename = get_saved_refmac_file_filename();
      if (mtz_filename) {
	gtk_label_set_text(GTK_LABEL(mtz_file_label), mtz_filename);
	fill_option_menu_with_refmac_file_labels_options(optionmenu);
      }
      if (refmac_use_twin_state()) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(twin_check_button), TRUE);
	gtk_widget_hide(sad_extras);
      } else {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(twin_check_button), FALSE);
      }
      if (get_refmac_phase_input() == 3) {
	gtk_widget_set_sensitive(twin_check_button, FALSE);
	gtk_widget_show(sad_extras);
	/* fill the entry with 1st existing atom */
	fill_refmac_sad_atom_entry(window);
      } else {
	gtk_widget_set_sensitive(twin_check_button, TRUE);
	gtk_widget_hide(sad_extras);
      }

    }
  } else {
    gtk_widget_show(labels);
  }

  optionmenu = lookup_widget(window, "run_refmac_ncycle_optionmenu");
  fill_option_menu_with_refmac_ncycle_options(optionmenu);
  //#endif

  /* set the ncs button depending on state */
  if (refmac_use_ncs_state()) {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ncs_button), TRUE);
  } else {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ncs_button), FALSE);
  }
#endif

  optionmenu = lookup_widget(window, "run_refmac_ccp4i_optionmenu");
  clear_refmac_ccp4i_project();
  add_ccp4i_projects_to_optionmenu(optionmenu, 
				   COOT_COORDS_FILE_SELECTION,
				   GTK_SIGNAL_FUNC(run_refmac_ccp4i_option_menu_signal_func));

  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(diff_map_button), TRUE);

  gtk_widget_show(window);
}


// Not needed? because we look at the active menu item at OK button-press time?
//
// Well, that was indeeed the way that we used to do it, now (WDW)
// that we rationalize the coordinates molecules option menu filling
// we have to set the the active molecule in this callback.
// 
void
refmac_molecule_button_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t::refmac_molecule = pos;

}

void fill_about_window(GtkWidget *widget) {

   GtkWidget *text_widget;
   text_widget = lookup_widget(GTK_WIDGET(widget), "about_window_text");

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

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
   gtk_text_insert (GTK_TEXT (text_widget), NULL, NULL, NULL,
		    _(widget_text.c_str()), -1);
#else   
  gtk_text_view_set_editable (GTK_TEXT_VIEW (text_widget), FALSE);
  gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (text_widget), GTK_WRAP_WORD);
  gtk_text_buffer_set_text (gtk_text_view_get_buffer (GTK_TEXT_VIEW (text_widget)),
			    _(widget_text.c_str()), -1);
#endif    
}

void add_coot_references_button(GtkWidget *widget) {

  GtkWidget *hbox;
  GtkWidget *button;
  hbox = GTK_DIALOG(widget)->action_area;
  button = gtk_button_new_with_label("References");
  gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, TRUE, 0);
#if (GTK_MAJOR_VERSION > 1)
  gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox), button, TRUE);
  gtk_box_reorder_child(GTK_BOX(hbox), button, 2);
  g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(wrapped_create_coot_references_dialog), NULL);
#endif // GTK_MAJOR_VERSION
  gtk_widget_show(button);
  
}

GtkWidget *wrapped_create_coot_references_dialog() {
  
#if (((GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 5)) || GTK_MAJOR_VERSION > 2)
  GtkWidget *references_dialog;
  GtkWidget *coot_reference_button;
  references_dialog = create_coot_references_dialog();
  coot_reference_button = lookup_widget(references_dialog, "coot_references_coot_toolbutton");
  g_signal_emit_by_name(G_OBJECT(coot_reference_button), "clicked");
  gtk_widget_show(references_dialog);
  return references_dialog;
#else
  GtkWidget *w = 0;
  return w;
#endif // GTK_MAJOR_VERSION

}


#ifdef COOT_USE_GTK2_INTERFACE
void fill_references_notebook(GtkToolButton *toolbutton, int reference_id) {

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

  notebook      = lookup_widget(GTK_WIDGET(toolbutton), "coot_references_notebook");
  ref_text_view = lookup_widget(GTK_WIDGET(toolbutton), "coot_references_textview");
  bib_text_view = lookup_widget(GTK_WIDGET(toolbutton), "coot_bibtext_textview");

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
     volume      = "";
     pages       = "in press";
     number      = "";
     
     bib_type = "Article";
     bib_id   = "emsley10:coot";
     bib_author = "Paul Emsley, Bernhard Lohkamp, William Scott and Kevin Cowtan";
     bib_journal = "Acta Crystallographica Section D - Biological Crystallography";
     
  }

  if (reference_id == COOT_REFERENCE_WINCOOT) {

    ref_description = "If you used WinCoot please cite Coot and the following:\n\n";

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
#endif // GTK_MAJOR_VERSION

void set_graphics_window_size(int x_size, int y_size) {

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      g.graphics_x_size = x_size;
      g.graphics_y_size = y_size;
      if (g.glarea) {
	 GtkWidget *win = lookup_widget(g.glarea, "window1");
	 GtkWindow *window = GTK_WINDOW(win);

#if (GTK_MAJOR_VERSION > 1)
         gtk_window_resize(window, x_size, y_size);
#else
	 // does this do a configure_event?  If so, then we don't need
	 // to do the graphics_draw() below.
         gtk_window_set_default_size(window, x_size, y_size);
#endif
	 while (gtk_events_pending())
	    gtk_main_iteration();
	 while (gdk_events_pending())
	    gtk_main_iteration();
// 	 std::cout << "DEBUG:: set " << window << " to size "
// 		   << x_size << " " << y_size << std::endl;
	 graphics_draw();
      }
      graphics_draw();
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("set-graphics-window-size");
   command_strings.push_back(graphics_info_t::int_to_string(x_size));
   command_strings.push_back(graphics_info_t::int_to_string(y_size));
   add_to_history(command_strings);
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
      GtkWidget *main = lookup_widget(g.glarea, "window1");
      if (main) { 
	 gtk_widget_set_uposition(main, x_pos, y_pos);
	 while (gtk_events_pending())
	    gtk_main_iteration();
      }
      graphics_draw();
   }
   std::string cmd = "set-graphics-window-position";
   std::vector<coot::command_arg_t> args;
   args.push_back(x_pos);
   args.push_back(y_pos);
   add_to_history_typed(cmd, args);
}

/* a general purpose version of the above, where we pass a widget flag */
void
store_window_position(int window_type, GtkWidget *widget) {

   // Note that we can't call this function on xxx_dialog_destroy
   // because we only have an gtkobject there, and if we cast there:
   // GTK_WIDGET(object) to get the 'widget' of this function, then
   // widget->window is NULL, and gdk_window_get_root_origin fails.

   gint upositionx, upositiony;

// BL says:: in windows root is not properly defined as in X11, so ok to use
// simple gdk_window_get_position function, I hope!
#ifdef WINDOWS_MINGW
   gdk_window_get_position (widget->window, &upositionx, &upositiony);
#else
   gdk_window_get_root_origin (widget->window, &upositionx, &upositiony);
#endif // MINGW

//    std::cout << "in store_window_position, widget is " << widget
//      	     << " widget->window is " << widget->window << std::endl;

   if (window_type == COOT_MODEL_REFINE_DIALOG) {
      graphics_info_t::model_fit_refine_x_position = upositionx;
      graphics_info_t::model_fit_refine_y_position = upositiony;
   } 

   if (window_type == COOT_DELETE_WINDOW) {
      // notice that for delete item, this does not get called from
      // the destroy callback, because we can't get to widget (the
      // dialog) that has a window
      graphics_info_t::delete_item_widget_x_position = upositionx;
      graphics_info_t::delete_item_widget_y_position = upositiony;
   }

   if (window_type == COOT_GO_TO_ATOM_WINDOW) {
      graphics_info_t::go_to_atom_window_x_position = upositionx;
      graphics_info_t::go_to_atom_window_y_position = upositiony;
   }
   
   if (window_type == COOT_ACCEPT_REJECT_WINDOW) {
      graphics_info_t::accept_reject_dialog_x_position = upositionx;
      graphics_info_t::accept_reject_dialog_y_position = upositiony;
   }
   
   if (window_type == COOT_ROTATE_TRANSLATE_DIALOG) {
      graphics_info_t::rotate_translate_x_position = upositionx;
      graphics_info_t::rotate_translate_y_position = upositiony;
   }

   if (window_type == COOT_DISPLAY_CONTROL_WINDOW) {
      graphics_info_t::display_manager_x_position = upositionx;
      graphics_info_t::display_manager_y_position = upositiony;
      graphics_info_t::display_manager_x_size = widget->allocation.width;
      graphics_info_t::display_manager_y_size = widget->allocation.height;
   }

   if (window_type == COOT_DISPLAY_CONTROL_MAPS_VBOX) {
      graphics_info_t::display_manager_maps_vbox_x_size =
	 widget->allocation.width;
      graphics_info_t::display_manager_maps_vbox_y_size =
	 widget->allocation.height;
   }
   if (window_type == COOT_DISPLAY_CONTROL_MOLECULES_VBOX) {
      graphics_info_t::display_manager_molecules_vbox_x_size =
	 widget->allocation.width;
      graphics_info_t::display_manager_molecules_vbox_y_size =
	 widget->allocation.height;
   }
   if (window_type == COOT_DISPLAY_CONTROL_PANE) {
      // This is a klude because this version of gtk doesn't seem to
      // have a nice accessor such as gtk_pane_get_position()
      GtkPaned *paned = GTK_PANED(widget);
      graphics_info_t::display_manager_paned_position = paned->child1_size;
   }

   if (window_type == COOT_EDIT_CHI_DIALOG) {
      graphics_info_t::edit_chi_angles_dialog_x_position = upositionx;
      graphics_info_t::edit_chi_angles_dialog_y_position = upositiony;
   }

   if (window_type == COOT_ROTAMER_SELECTION_DIALOG) {
      graphics_info_t::rotamer_selection_dialog_x_position = upositionx;
      graphics_info_t::rotamer_selection_dialog_y_position = upositiony;
   }

   if (window_type == COOT_RAMACHANDRAN_PLOT_WINDOW) {
      graphics_info_t::ramachandran_plot_x_position = upositionx;
      graphics_info_t::ramachandran_plot_y_position = upositiony;
   }

   if (window_type == COOT_DISTANCES_ANGLES_WINDOW) {
      graphics_info_t::distances_and_angles_dialog_x_position = upositionx;
      graphics_info_t::distances_and_angles_dialog_y_position = upositiony;
   }
}

/* a general purpose version of the above, where we pass a widget flag */
void
store_window_size(int window_type, GtkWidget *widget) {

   if (window_type == COOT_FILESELECTION_DIALOG) {
      graphics_info_t::file_selection_dialog_x_size = widget->allocation.width;
      graphics_info_t::file_selection_dialog_y_size = widget->allocation.height;
   }
}

void set_file_selection_dialog_size(GtkWidget *dialog) {

   if (graphics_info_t::file_selection_dialog_x_size > 0) {
#if (GTK_MAJOR_VERSION > 2 || (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION > 2))
      if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
//          gtk_widget_set_size_request(dialog,
//  			   graphics_info_t::file_selection_dialog_x_size,
// 			   graphics_info_t::file_selection_dialog_y_size);
         gtk_window_set_default_size(GTK_WINDOW(dialog),
				     graphics_info_t::file_selection_dialog_x_size,
				     graphics_info_t::file_selection_dialog_y_size);
      } else {
         gtk_window_resize(GTK_WINDOW(dialog),
 			   graphics_info_t::file_selection_dialog_x_size,
			   graphics_info_t::file_selection_dialog_y_size);
      }
#else
      gtk_widget_set_usize(dialog,
			   graphics_info_t::file_selection_dialog_x_size,
			   graphics_info_t::file_selection_dialog_y_size);
#endif // GTK_MAJOR
   }
}


/* return negative if fail */
float get_positive_float_from_entry(GtkEntry *w) {

   float f = -1.0;
   if (graphics_info_t::use_graphics_interface_flag) { 
      const gchar *text = gtk_entry_get_text(w);
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
gboolean
coot_checked_exit(int retval) {

   //    cout << "exitting with status " << retval << endl;


   graphics_info_t g;
   int i_unsaved = g.check_for_unsaved_changes();
   std::string cmd = "coot-checked-exit";
   std::vector<coot::command_arg_t> args;
   args.push_back(retval);
   add_to_history_typed(cmd, args);
   if (i_unsaved == 0) { // no unsaved.
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
   return TRUE; // path where there were unsaved changes, we don't
		// want to exit.
}

#ifdef USE_GUILE
void run_clear_backups(int retval) {

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
} 
#endif


#ifdef USE_PYTHON
void run_clear_backups_py(int retval) {

   PyObject *r = safe_python_command_with_return("clear_backups_maybe()");

   if (r == NULL || r == Py_None) { 
      // not false and not not false, function didn't run then...
      std::cout << "WARNING:: clear_backups_maybe() returns "
		<< PyString_AsString(PyObject_Str(r))
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

   save_state(); // we get error message in save_state()

   // save the history
   graphics_info_t g;
   if (! g.disable_state_script_writing)
      g.save_history();

   // std::cout << "mapview_real_exit" << std::endl;

#ifdef USE_MYSQL_DATABASE

   db_finish_up();

#endif // USE_MYSQL_DATABASE   

#ifdef USE_PYTHON
   // Py_Finalize();
#endif

   clipper::ClipperInstantiator::instance().destroy();

   gtk_exit(retval); 

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

   bool no_chooser_filter = 1;
   GtkWidget *button = 0;

#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      	no_chooser_filter = 0;
	add_filechooser_filter_button(fileselection, data_type);
   }
#endif

   if (no_chooser_filter) {
    GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
    GtkWidget *frame = gtk_frame_new("File-name filter:");
    int d = data_type;
    button = gtk_toggle_button_new_with_label("Filter");
 //       std::cout << "in add_filename_filter_button data_type is " << data_type << std::endl;

    gtk_widget_ref(button);
    gtk_widget_show(button);
    gtk_container_add(GTK_CONTAINER(aa),frame);
    gtk_container_add(GTK_CONTAINER(frame), button);
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
    gtk_signal_connect (GTK_OBJECT (button), "toggled",
			GTK_SIGNAL_FUNC (on_filename_filter_toggle_button_toggled_gtk1),
			GINT_TO_POINTER(d));
#else
    // callback in c-interface-gtk2.cc
    gtk_signal_connect (GTK_OBJECT (button), "toggled",
			GTK_SIGNAL_FUNC (on_filename_filter_toggle_button_toggled),
			GINT_TO_POINTER(d));
#endif   
    gtk_widget_show(frame);
   }
   
   return button;
}


// Paul requested a new function for filechooser filter
// only available in gtk2
// here we go

#if (GTK_MAJOR_VERSION > 1)
// where data type:
// 0 coords
// 1 mtz etc
// 2 maps
// 3 cif dictionary
// 4 scripting files
// 
void add_filechooser_filter_button(GtkWidget *fileselection, 
				      short int data_type) { 
   
  int d = data_type;

  int i = 0;
  std::vector<std::string> globs;

  GtkFileFilter *filterall = gtk_file_filter_new ();
  GtkFileFilter *filterselect = gtk_file_filter_new ();

  gtk_file_filter_set_name (filterall, "all-files");
  gtk_file_filter_add_pattern (filterall, "*");
      
  if (d == COOT_COORDS_FILE_SELECTION || d == COOT_SAVE_COORDS_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "coordinate-files");

    globs = *graphics_info_t::coordinates_glob_extensions;
  };

  if (d == COOT_DATASET_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "data-files");

    globs = *graphics_info_t::data_glob_extensions;
  };

  if (d == COOT_MAP_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "map-files");

    globs = *graphics_info_t::map_glob_extensions;
  };

  if (d == COOT_CIF_DICTIONARY_FILE_SELECTION) {

    gtk_file_filter_set_name (filterselect, "dictionary-files");

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

    gtk_file_filter_set_name (filterselect, "scripting-files");

    globs = script_glob_extension;

  };

  std::string s;
  for (unsigned int i=0; i<globs.size(); i++) {
    s = "*";
    s += globs[i];
    gtk_file_filter_add_pattern (filterselect, s.c_str());
  };

  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER (fileselection),
			       GTK_FILE_FILTER (filterall));
  gtk_file_chooser_add_filter (GTK_FILE_CHOOSER (fileselection),
			       GTK_FILE_FILTER (filterselect));

  if (filter_fileselection_filenames_state() == 1) {
    // filter automatically
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER (fileselection),
				GTK_FILE_FILTER (filterselect));
  } else {
    // show all
    gtk_file_chooser_set_filter(GTK_FILE_CHOOSER (fileselection),
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


#endif // GTK_MAJOR_VERSION

void
on_read_map_difference_map_toggle_button_toggled (GtkButton       *button,
						  gpointer         user_data)
{
   if (GTK_TOGGLE_BUTTON(button)->active) { 
      std::cout << "is a difference map...!\n";
   } 
}

#if (GTK_MAJOR_VERSION > 1) || defined (GTK_ENABLE_BROKEN)
#else 
void
on_filename_filter_toggle_button_toggled_gtk1(GtkButton       *button,
					      gpointer         user_data)
{
   int int_user_data = GPOINTER_TO_INT(user_data);
   int data_type = int_user_data & 31; // lower 5 bits
   int file_selection_type = data_type; 

   // We need to add text to the string of the dictectory we are in
   // (pre_directory), so first we need to find pre_directory (as per
   // fileselection_sort_button_clicked()
   // 
   GtkWidget *sort_button = lookup_widget(GTK_WIDGET(button),
					  "fileselection_sort_button");
   if (sort_button) { 
      // std::cout << "Hooray! we found the sort button!\n";
      // usually, we do.
   } else { 
      std::cout << "Boo we failed to find the sort button!\n";
   } 
   std::string pre_directory = pre_directory_file_selection(sort_button);
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button,
							    file_selection_type);
   
   if (fileselection) { 
      if (GTK_TOGGLE_BUTTON(button)->active) { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Unfilter");
	 
	 // so now we have pre_directory
	 // 
	 // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
	 std::vector<std::string> v = filtered_by_glob(pre_directory, data_type);
	 // std::cout << "DEBUG:: filtering by glob using data type: " << data_type
	 // << " returns" << std::endl;
	 // for (unsigned int iv=0; iv< v.size(); iv++)
	 // std::cout << iv << " " << v[iv] << std::endl;

	 filelist_into_fileselection_clist(fileselection, v);

      } else { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Filter");
	 gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
					 (pre_directory + "/").c_str());
      }
   } else {
      std::cout << "ERROR:: no fileselection found from sort button\n";
      std::cout << "ERROR:: lookup of file selection from " << sort_button
		<< " failed for file selection type: " << file_selection_type
		<< std::endl;
   }
}
#endif // GTK_MAJOR_VERSION > 1 or BROKEN

#include <gdk/gdkkeysyms.h> // for keyboarding.

gboolean
on_filename_filter_key_press_event (GtkWidget       *widget,
				    GdkEventKey     *event,
				    gpointer         user_data)
{
   //    if (event->keyval == GDK_Return || event->keyval == GDK_Tab)
   //    { // Tab is not good.  It takes you to the next widget too,
   //    which is not want we want.  It's confusing, so let's just use
   //    return.
   if (event->keyval == GDK_Return) { 
      // std::cout << "Return pressed!\n";
#ifdef COOT_USE_GTK2_INTERFACE
      handle_filename_filter_gtk2(widget);
#else       
      handle_filename_filter_gtk1(widget);
#endif       
   } 
   return FALSE;
}


// This is for the "bespoke" filename filtering, reading the glob from
// an entry added to the action area.
void
handle_filename_filter_gtk1(GtkWidget *entry_widget) {
   
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   // nothing
#else 

   std::cout << "running handle_filename_filter\n";

   const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry_widget));

   // We need to add text to the string of the dictectory we are in
   // (pre_directory), so first we need to find pre_directory (as per
   // fileselection_sort_button_clicked()
   // 
   GtkWidget *sort_button = lookup_widget(entry_widget, "fileselection_sort_button");
   if (sort_button) { 
      // std::cout << "Hooray! we found the sort button!\n";
      // usually, we do.
   } else { 
      std::cout << "Boo we failed to find the sort button!\n";
   } 
   std::string pre_directory = pre_directory_file_selection(sort_button);

   int file_selection_type = COOT_SCRIPTS_FILE_SELECTION; 

   // so now we have pre_directory
   // 
   // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button,
							    file_selection_type);
   if (fileselection) { 
      GtkCList  *file_list = GTK_CLIST(GTK_FILE_SELECTION(fileselection)->file_list);
      std::string file_name_glob;
      file_name_glob = pre_directory;
      file_name_glob += "/";
      file_name_glob += text;

      glob_t myglob;
      int flags = 0;
      glob(file_name_glob.c_str(), flags, 0, &myglob);
      size_t count;


      char **p;
      std::vector<std::string> v;
      for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) { 
	 v.push_back(std::string(*p));
      }
      globfree(&myglob);

      if (v.size() > 0) { 
	 gtk_clist_clear(file_list);
	 std::string::size_type islash;
	 std::string t;
	 for (unsigned int i=0; i<v.size(); i++) {
	    islash = v[i].find_last_of("/");
	    if (islash == string::npos) { 
	       // no slash found:
	       t = v[i];
	    } else {
	       t = v[i].substr(islash + 1);
	    }
	    char *text = new char[t.length()+1];
	    strncpy(text, t.c_str(), t.length()+1);
	    gtk_clist_append(file_list, &text);
	 }
      }
   } else { 
      std::cout << "ERROR:: couldn't find fileselection\n";
   } 

#endif // WINDOWS_MINGW
}

std::string pre_directory_file_selection(GtkWidget *sort_button) { 

   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(gtk_object_get_user_data(GTK_OBJECT(sort_button)));

   // The menu item is a container than contains a label.
   // How do we get to the label given the container?
   // Strangely enough we use the history_pulldown.

   GList *dlist = gtk_container_children(GTK_CONTAINER(history_pulldown));
   GList *free_list = dlist;
   std::string pre_directory("");
   
   while (dlist) {
      // GtkWidget *list_item;
      // list_item = (GtkWidget *) (dlist->data);
      gchar *t = GTK_LABEL(dlist->data)->label;
      if (t != NULL) {
	 pre_directory = t; 
      } else {
	 std::cout << "WARNING:: null label t " << std::endl;
      } 
      dlist = dlist->next;
   }
   g_list_free(free_list); 

   return pre_directory;
}




void push_the_buttons_on_fileselection(GtkWidget *filter_button, 
				       GtkWidget *sort_button,
				       GtkWidget *fileselection) {

  bool no_chooser = 1;
#if (GTK_MAJOR_VERSION > 1)
  if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE)
    {
      no_chooser = 0;
    }
#endif
  if (filter_fileselection_filenames_state() && no_chooser) {
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
    gtk_signal_emit_by_name(GTK_OBJECT(filter_button), "clicked");    
#else
    g_signal_emit_by_name(G_OBJECT(filter_button), "clicked");
#endif
    std::cout << "INFO:: Filtering file names \n";
  }
  if (graphics_info_t::sticky_sort_by_date) {
    GtkWidget *file_list;
    if (no_chooser) {
      std::cout << "INFO:: Sorting files by date\n";
      file_list = GTK_FILE_SELECTION(fileselection)->file_list;
    }
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
     fileselection_sort_button_clicked_gtk1(sort_button, (GtkCList *) file_list);
#else
     if (no_chooser) {
       g_signal_emit_by_name(G_OBJECT(sort_button), "clicked", file_list);
     }
#endif     
  }
}

void 
filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v) {

   GtkCList  *file_list = GTK_CLIST(GTK_FILE_SELECTION(fileselection)->file_list);
   gtk_clist_clear(file_list);
   std::string::size_type islash;
   std::string t;
   for (unsigned int i=0; i<v.size(); i++) {
      islash = v[i].find_last_of("/");
      if (islash == string::npos) { 
	 // no slash found:
	 t = v[i];
      } else {
	 t = v[i].substr(islash + 1);
      }
      char *text = new char[t.length()+1];
      strncpy(text, t.c_str(), t.length()+1);
      gtk_clist_append(file_list, &text);
   }
}

/*  Eleanor likes to sort her files by date when selecting a file
*/
GtkWidget *add_sort_button_fileselection(GtkWidget *fileselection) {

   GtkWidget *button = 0;
   bool doit = 1;

#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::gtk2_file_chooser_selector_flag == 1) {
      doit = 0;
   }
#endif
   
   if (doit) {
      GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
      GtkWidget *frame = gtk_frame_new("File Order");
      button = gtk_button_new_with_label("  Sort by Date  ");
      gtk_widget_ref(button);
      gtk_object_set_data_full(GTK_OBJECT(fileselection),
			       "fileselection_sort_button",
			       button,
			       (GtkDestroyNotify) gtk_widget_unref);
      GtkWidget *file_list = GTK_FILE_SELECTION(fileselection)->file_list;
      
      GtkOptionMenu *history_pulldown =
	 GTK_OPTION_MENU(GTK_FILE_SELECTION(fileselection)->history_pulldown);
      
      gtk_object_set_user_data(GTK_OBJECT(button), (char *) history_pulldown); 
      
      
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
      gtk_signal_connect (GTK_OBJECT(button), "clicked",
			  (GtkSignalFunc) fileselection_sort_button_clicked_gtk1,
			  file_list);
#else
      g_signal_connect (G_OBJECT(button), "clicked",
			G_CALLBACK(fileselection_sort_button_clicked),
			file_list);
#endif
      
      gtk_container_add(GTK_CONTAINER(aa),frame);
      gtk_container_add(GTK_CONTAINER(frame), button);
      gtk_widget_show(frame);
      gtk_widget_show(button);

   } else {
	// we have the chooser and dont need a sort button
   }
   return button;
}

void add_is_difference_map_checkbutton(GtkWidget *fileselection) { 

  bool add_map_button = 1;
#if (GTK_MAJOR_VERSION > 1)
  if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      add_map_button = 0;
  }
#endif // GTK_MAJOR_VERSION
  if (add_map_button) {
   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *button = gtk_check_button_new_with_label("Is Difference Map");
   GtkWidget *frame = gtk_frame_new("Difference Map?");

   gtk_widget_ref(button);
   gtk_object_set_data_full(GTK_OBJECT(fileselection),
			    "map_fileselection_is_difference_map_button",
			    button,
			    (GtkDestroyNotify) gtk_widget_unref);
   gtk_widget_show(button);
   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), button);
   gtk_signal_connect (GTK_OBJECT (button), "toggled",
		       GTK_SIGNAL_FUNC (on_read_map_difference_map_toggle_button_toggled),
		       NULL);
   gtk_widget_show(frame);
  }

}

void add_recentre_on_read_pdb_checkbutton(GtkWidget *fileselection) { 

   bool doit = 1;
#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      doit = 0;
      GtkWidget *combobox =
	 lookup_widget(GTK_WIDGET(fileselection),
		       "coords_filechooserdialog1_recentre_combobox");
      gtk_combo_box_append_text(GTK_COMBO_BOX(combobox), "Recentre on Molecule");
      gtk_combo_box_append_text(GTK_COMBO_BOX(combobox), "Dont't Recentre");
      gtk_combo_box_append_text(GTK_COMBO_BOX(combobox), "Recentre Molecule Here");
      if (graphics_info_t::recentre_on_read_pdb)
     	  gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 0);
      if (!graphics_info_t::recentre_on_read_pdb)
     	  gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 1);
   }
#endif
   
   if (doit) {

      if (0) { 
      
	 GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
	 GtkWidget *button = gtk_check_button_new_with_label("Recentre");
	 GtkWidget *frame = gtk_frame_new("Recentre?");
	 GtkTooltips *tooltips = gtk_tooltips_new();

	 gtk_widget_ref(button);
	 gtk_object_set_data_full(GTK_OBJECT(fileselection),
				  "coords_fileselection1_recentre_checkbutton",
				  button,
				  (GtkDestroyNotify) gtk_widget_unref);
	 gtk_tooltips_set_tip(tooltips, button,
			      _("Deactivate this checkbutton if you don't want to change the view centre when these new coordinates are read"), NULL);

	 // shall we activate the button?
	 if (graphics_info_t::recentre_on_read_pdb)
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	 gtk_widget_show(button);
	 gtk_container_add(GTK_CONTAINER(aa),frame);
	 gtk_container_add(GTK_CONTAINER(frame), button);
	 gtk_signal_connect (GTK_OBJECT (button), "toggled",
			     GTK_SIGNAL_FUNC (on_recentre_on_read_pdb_toggle_button_toggled),
			     NULL);
	 gtk_widget_show(frame);
      } else {

	 GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
	 GtkWidget *om = gtk_option_menu_new();
	 GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(om));
	 if (menu)
	    gtk_widget_destroy(menu);
	 menu = GTK_WIDGET(gtk_menu_new());
	 GtkWidget *frame = gtk_frame_new("Recentre?");
	 GtkTooltips *tooltips = gtk_tooltips_new();

	 gtk_widget_ref(om);
	 gtk_object_set_data_full(GTK_OBJECT(fileselection),
				  "coords_fileselection1_recentre_optionmenu",
				  om,
				  (GtkDestroyNotify) gtk_widget_unref);
// 	 gtk_tooltips_set_tip(tooltips, om,
// 			      _("Deactivate this checkbutton if you don't want to change the view centre when these new coordinates are read"), NULL);

	 GtkWidget *menuitem;

	 menuitem = gtk_menu_item_new_with_label("Recentre on Molecule");
	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 gtk_widget_show(menuitem);
	 // shall we activate the button?
	 if (graphics_info_t::recentre_on_read_pdb)
	    gtk_menu_item_activate(GTK_MENU_ITEM(menuitem));
	 //
	 menuitem = gtk_menu_item_new_with_label("Don't Recentre");
	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 if (!graphics_info_t::recentre_on_read_pdb)
	    gtk_menu_item_activate(GTK_MENU_ITEM(menuitem));
	 gtk_widget_show(menuitem);
	 // 
	 menuitem = gtk_menu_item_new_with_label("Recentre Molecule Here");
	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 gtk_widget_show(menuitem);
	 


	 gtk_widget_show(om);
	 gtk_container_add(GTK_CONTAINER(aa),frame);
	 gtk_container_add(GTK_CONTAINER(frame), om);
	 gtk_option_menu_set_menu(GTK_OPTION_MENU(om), menu);
	 gtk_widget_show(frame);
      }
   }
}


void
on_recentre_on_read_pdb_toggle_button_toggled (GtkButton       *button,
					       gpointer         user_data)
{
   if (GTK_TOGGLE_BUTTON(button)->active) { 
      std::cout << "INFO:: activated recentering on new coordinates.\n";
   } else {
      std::cout << "INFO:: de-activated recentering on new coordinates.\n";
   }
}


/*  ------------------------------------------------------------------------ */
/*              scripting gtk interface                                      */
/*  ------------------------------------------------------------------------ */

// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_python_window_entry(GtkWidget *entry) { 

#ifdef USE_PYTHON

   // add python entry in entry callback code here...

#if (GTK_MAJOR_VERSION > 1) 
    g_signal_connect(G_OBJECT(entry), "activate",
		     G_CALLBACK(python_window_enter_callback),
		     (gpointer) entry);
#else
    gtk_signal_connect(GTK_OBJECT(entry), "activate",
		       GTK_SIGNAL_FUNC(python_window_enter_callback),
		       entry);
# endif // GTK_MAJOR_VERSION

#endif // PYTHON

}

// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_guile_window_entry(GtkWidget *entry) { 

#ifdef USE_GUILE
   gtk_signal_connect(GTK_OBJECT(entry), "activate",
		      GTK_SIGNAL_FUNC(guile_window_enter_callback),
		      entry);

#endif //  USE_GUILE

}

#ifdef USE_PYTHON
void python_window_enter_callback( GtkWidget *widget,
				   GtkWidget *entry )
{
#if (GTK_MAJOR_VERSION > 1) 
  const gchar *entry_text;
#else
  char *entry_text;
#endif
  entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
  printf("Entry contents: %s\n", entry_text);

  // Sigh. PyRun_SimpleString needs a (char *), not a (const gchar *):
  size_t new_length = strlen(entry_text)+1;
  char *new_text;
  new_text = new char[new_length];
  strncpy(new_text, entry_text, new_length);
  printf("Running string: %s\n", new_text);
  
  PyRun_SimpleString(new_text);
  
  // clear the entry
  gtk_entry_set_text(GTK_ENTRY(entry),"");

  delete [] new_text;

}
#endif


#ifdef USE_GUILE
void guile_window_enter_callback( GtkWidget *widget,
				  GtkWidget *entry )
{
  const gchar *entry_text;
  entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
  printf("Entry contents: %s\n", entry_text);

  // scm_c_eval_string(entry_text); 

  // extern SCM scm_catch (SCM tag, SCM thunk, SCM handler);
  //

  //SCM handler = scm_c_eval_string("'(lambda (key . args))"); 

   SCM handler = scm_c_eval_string ("(lambda (key . args) "
     "(display (list \"Error in proc:\" key \" args: \" args)) (newline))"); 
			   // "(newline))"); 

  // scm_catch(SCM_BOOL_T, scm_c_eval_string(entry_text), handler); 

  std::string thunk("(lambda() "); 
  thunk += entry_text; 
  thunk += " )";


  SCM scm_thunk = scm_c_eval_string(thunk.c_str()); 
  scm_catch(SCM_BOOL_T, scm_thunk, handler);


  // clear the entry
  // strcpy(entry_text, "");  surely not needed.
  gtk_entry_set_text(GTK_ENTRY(entry),"");
}
#endif //  USE_GUILE


// Similar to fill_option_menu_with_coordinates_options, but I moved
// it to graphics_info_t because it is also used when there is an
// ambiguity in the map for refinement (graphics_info_t::refine)
// 
int fill_option_menu_with_map_options(GtkWidget *option_menu, GtkSignalFunc signalfunc) {

   graphics_info_t g;

   return g.fill_option_menu_with_map_options(option_menu, signalfunc);
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

   GtkWidget *button = lookup_widget(window, "single_map_sigma_checkbutton");
   GtkWidget *entry  = lookup_widget(window, "single_map_sigma_step_entry");

   if (GTK_TOGGLE_BUTTON(button)->active) { 
      const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
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

   if (graphics_info_t::glarea) { 
      GtkWindow *main_window =
	 GTK_WINDOW(lookup_widget(graphics_info_t::glarea, "window1"));
      gtk_window_set_transient_for(GTK_WINDOW(window), main_window);
      if (widget_type == COOT_DELETE_WINDOW) {
	 if (graphics_info_t::delete_item_widget_x_position > -100) {
	    if (graphics_info_t::delete_item_widget_y_position > -100) {
	       gtk_widget_set_uposition(window,
					graphics_info_t::delete_item_widget_x_position,
					graphics_info_t::delete_item_widget_y_position);
	    }
	 }
      }
      if (widget_type == COOT_DISTANCES_ANGLES_WINDOW) {
	 if (graphics_info_t::distances_and_angles_dialog_x_position > -100) {
	    if (graphics_info_t::distances_and_angles_dialog_y_position > -100) {
	       gtk_widget_set_uposition(window,
					graphics_info_t::distances_and_angles_dialog_x_position,
					graphics_info_t::distances_and_angles_dialog_y_position);
	    }
	 }
      }
   }
}

GtkWidget *coot_file_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_coords_fileselection1 ();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_coords_fileselection1();
      gtk_file_selection_set_select_multiple(GTK_FILE_SELECTION(w), TRUE);
   } else {
      w = create_coords_filechooserdialog1(); 
      gtk_file_chooser_set_select_multiple(GTK_FILE_CHOOSER(w), TRUE);
   }
#endif
   return w;
}

GtkWidget *coot_dataset_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_dataset_fileselection1 ();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_dataset_fileselection1();
   } else {
      w = create_dataset_filechooserdialog1(); 
   }
#endif
   return w;
}

GtkWidget *coot_map_name_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_map_name_fileselection1();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_map_name_fileselection1();
   } else {
      w = create_map_name_filechooserdialog1(); 
   }
#endif
   return w;
}

GtkWidget *coot_save_coords_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_save_coords_fileselection1 ();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_save_coords_fileselection1();
   } else {
      w = create_save_coords_filechooserdialog1();
#if (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION < 10)
      // we don't have confirmation overwrite
#else      
      gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (w), TRUE);
#endif      
   }
#endif
   return w;
}

GtkWidget *coot_cif_dictionary_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_cif_dictionary_fileselection ();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_cif_dictionary_fileselection();
   } else {
      w = create_cif_dictionary_filechooserdialog1(); 
   }
#endif
   return w;
}

GtkWidget *coot_run_script_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_run_script_fileselection();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_run_script_fileselection();
   } else {
      w = create_run_script_filechooserdialog1(); 
   }
#endif
   return w;
}

GtkWidget *coot_save_state_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_save_state_fileselection();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_save_state_fileselection();
   } else {
      w = create_save_state_filechooserdialog1(); 

#if (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION < 10)
      // we don't have confirmation overwrite
#else      
      gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (w), TRUE);
#endif      
      
   }
#endif
   return w;
}

GtkWidget *coot_save_symmetry_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_save_symmetry_coords_fileselection();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_save_symmetry_coords_fileselection();
   } else {
      w = create_save_symmetry_coords_filechooserdialog1();

#if (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION < 10)
      // we don't have confirmation overwrite
#else      
      gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (w), TRUE);
#endif      
   }
#endif
   return w;
}

GtkWidget *coot_screendump_chooser() {

   GtkWidget *w;

#if (GTK_MAJOR_VERSION == 1)
   w = create_screendump_fileselection();
#else
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::OLD_STYLE) {
      w = create_screendump_fileselection();
   } else {
      w = create_screendump_filechooserdialog1(); 

#if (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION < 10)
      // we don't have confirmation overwrite
#else      
      gtk_file_chooser_set_do_overwrite_confirmation (GTK_FILE_CHOOSER (w), TRUE);
#endif      

   }
#endif
   return w;

}


void set_directory_for_coot_file_chooser(GtkWidget *coords_fileselection1) {

#if (GTK_MAJOR_VERSION == 1)
      set_directory_for_fileselection(coords_fileselection1);
#else
      if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
	set_directory_for_filechooser(coords_fileselection1);
      } else {
        set_directory_for_fileselection(coords_fileselection1);
      }
#endif
}

const char *coot_file_chooser_file_name(GtkWidget *widget) {

   const char *f = 0;
#if (GTK_MAJOR_VERSION == 1)
#else

#endif
   return f;
}



/*  ----------------------------------------------------------------------- */
/*                  get by accession code:                                  */
/*  ----------------------------------------------------------------------- */


/* Accession code, and dispatch guile command to download and display
   the model.  Hmmm.  */
void handle_get_accession_code(GtkWidget *widget) { 

   const gchar *text = gtk_entry_get_text(GTK_ENTRY(widget));
   cout << "PDB Accession Code: " << text << endl;
   int *n = (int *) gtk_object_get_user_data(GTK_OBJECT(lookup_widget(GTK_WIDGET(widget),
								      "accession_code_window")));
#ifdef USE_GUILE
   string scheme_command;

   if (*n == COOT_ACCESSION_CODE_WINDOW_EDS) {
      // 20050725 EDS code:

      scheme_command = "(get-eds-pdb-and-mtz ";
      scheme_command += single_quote(text);
      scheme_command += ")";
   } else { 
      if (*n == 1) { 
	 scheme_command = "(get-ebi-pdb \"";
      } else { 
	 // *n == 2 see callbacks.c on_get_pdb_and_sf_using_code1_activate
	 scheme_command = "(get-ebi-pdb-and-sfs \"";
      }
      scheme_command += text;
      scheme_command += "\")";
   }

   safe_scheme_command(scheme_command); 

#else 
   
#ifdef USE_PYTHON
      string python_command;

   if (*n == COOT_ACCESSION_CODE_WINDOW_EDS) {
      // 20050725 EDS code:

      python_command = "get_eds_pdb_and_mtz(";
      python_command += single_quote(text);
      python_command += ")";
   } else {
      if (*n == 1) {
         python_command = "get_ebi_pdb(";
      } else {
         // *n == 2 see callbacks.c on_get_pdb_and_sf_using_code1_activate
         python_command = "get_ebi_pdb_and_sfs(";
      }
      python_command += single_quote(text);
      python_command += ")";
   }

   safe_python_command(python_command);
#else 
   std::cout << "WARING:: Executable not compiled with guile or python." << std::endl;
   std::cout << "         This won't work." << std::endl; 

#endif // USE_PYTHON

#endif // USE_GUILE

   // and kill the accession code window
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(widget),
				    "accession_code_window")); 
} 




// Set the internal state of the torsion restraints and de-sensitize
// the peptide restraints if the torsion toggle button goes off.
// 
void do_torsions_toggle(GtkWidget *togglebutton) {

   graphics_info_t g;
   GtkWidget *peptide_checkbutton =
      lookup_widget(togglebutton,
		    "refine_params_use_peptide_torsions_checkbutton");

   if (GTK_TOGGLE_BUTTON(togglebutton)->active)
      g.do_torsion_restraints = 1;
   else
      g.do_torsion_restraints = 0;

   if (g.do_torsion_restraints == 0) {
      gtk_widget_set_sensitive(peptide_checkbutton, FALSE);
   } else {
      gtk_widget_set_sensitive(peptide_checkbutton, TRUE);
   }
}

GtkWidget *wrapped_create_refine_params_dialog() {

   GtkWidget *w = create_refine_params_dialog();
   set_refine_params_toggle_buttons(w);
   return w;
}

void set_refine_params_toggle_buttons(GtkWidget *button) {

   // initiallly buttons are inactive and sensitive

   graphics_info_t g;
   GtkWidget *checkbutton =
      lookup_widget(button, "refine_params_use_torsions_checkbutton");
   GtkWidget *planar_peptide_restraints_checkbutton =
      lookup_widget(button, "refine_params_use_planar_peptides_checkbutton");
   GtkWidget *phi_psi_peptide_checkbutton =
      lookup_widget(button, "refine_params_use_peptide_torsions_checkbutton");
   GtkWidget *link_torsion_type_vbox =
      lookup_widget(button, "peptide_torsions_restraints_vbox");
   // refine_params_use_ramachandran_goodness_torsions_checkbutton
   GtkWidget *rama_button =
      lookup_widget(button, "refine_params_use_ramachandran_goodness_torsions_checkbutton");
   
   if (g.do_torsion_restraints) {
      g.do_torsion_restraints = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(phi_psi_peptide_checkbutton), FALSE);
   }


   GtkWidget *omega = lookup_widget(button,
		       "refine_params_use_peptide_omegas_checkbutton");
   if (g.do_peptide_omega_torsion_restraints) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(omega), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(omega), FALSE);
   }

   // planar peptide restraints:
   int pps = planar_peptide_restraints_state();
   if (pps) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(planar_peptide_restraints_checkbutton), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(planar_peptide_restraints_checkbutton), FALSE);
   }      

   // refine_params_use_helix_peptide_torsions_radiobutton
   // refine_params_use_beta_strand_peptide_torsions_radiobutton
   // refine_params_use_ramachandran_goodness_torsions_radiobutton

   GtkWidget *sec_str_rest_no_rest_radiobutton =
      lookup_widget(button, "sec_str_rest_no_rest_radiobutton");
   GtkWidget *sec_str_rest_helix_rest_radiobutton = 
      lookup_widget(button, "sec_str_rest_helix_rest_radiobutton");
   GtkWidget *sec_str_rest_strand_rest_radiobutton = 
      lookup_widget(button, "sec_str_rest_strand_rest_radiobutton");

#ifdef HAVE_GSL   
   if (graphics_info_t::pseudo_bonds_type == coot::NO_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_no_rest_radiobutton), TRUE);
   if (graphics_info_t::pseudo_bonds_type == coot::HELIX_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_helix_rest_radiobutton), TRUE);
   if (graphics_info_t::pseudo_bonds_type == coot::STRAND_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_strand_rest_radiobutton), TRUE);
   if (graphics_info_t::do_rama_restraints)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rama_button), TRUE);
#endif // HAVE_GSL

   GtkWidget *refinement_weight_entry = lookup_widget(button, "refine_params_weight_matrix_entry");
   if (refinement_weight_entry) {
      gtk_entry_set_text(GTK_ENTRY(refinement_weight_entry),
			 coot::util::float_to_string(g.geometry_vs_map_weight).c_str());
   } 
} 

void fill_chiral_volume_molecule_option_menu(GtkWidget *w) { 

   GtkWidget *optionmenu = lookup_widget(w, "check_chiral_volumes_molecule_optionmenu");

   // now set chiral_volume_molecule_option_menu_item_select_molecule to the top of the list
   for (int i=0; i<graphics_info_t::n_molecules(); i++) { 
      if (graphics_info_t::molecules[i].has_model()) {
	 graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule = i;
	 break;
      } 
   }
   int imol = graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule;
   GtkSignalFunc callback_func =
      GTK_SIGNAL_FUNC(chiral_volume_molecule_option_menu_item_select);

   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(optionmenu, callback_func, imol);

}


/*  ------------------------------------------------------------------------ */
//            enviromnent and other distances
/*  ------------------------------------------------------------------------ */
//


void toggle_environment_show_distances(GtkToggleButton *button) {

   graphics_info_t g;

   //    if (g.environment_show_distances == 0) {
   
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(button),
				   "environment_distance_distances_frame");
   GtkWidget *label_atom_check_button =
      lookup_widget(GTK_WIDGET(button),
		    "environment_distance_label_atom_checkbutton");
   
   if (button->active) { 
      // std::cout << "toggled evironment distances on" << std::endl;
      g.environment_show_distances = 1;
      gtk_widget_set_sensitive(hbox, TRUE);
      gtk_widget_set_sensitive(label_atom_check_button, TRUE);

      //
      std::pair<int, int> r =  g.get_closest_atom();
//       std::cout << "DEBUG:: got close info: " 
// 		<< r.first << " " << r.second << std::endl;
      if (r.first >= 0) {
	 g.mol_no_for_environment_distances = r.second;
	 g.update_environment_distances_maybe(r.first, r.second);
	 graphics_draw();
      }
      
   } else {
      // std::cout << "toggled evironment distances off" << std::endl;
      g.environment_show_distances = 0;
      gtk_widget_set_sensitive(hbox, FALSE);
      gtk_widget_set_sensitive(label_atom_check_button, FALSE);
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

   GtkWidget *min_entry   = lookup_widget(widget, "pointer_distances_min_dist_entry");
   GtkWidget *max_entry   = lookup_widget(widget, "pointer_distances_max_dist_entry");
   GtkWidget *checkbutton = lookup_widget(widget, "pointer_distances_checkbutton");
   GtkWidget *frame       = lookup_widget(widget, "pointer_distances_frame");

   float min_dist = graphics_info_t::pointer_min_dist;
   float max_dist = graphics_info_t::pointer_max_dist;
   
   gtk_entry_set_text(GTK_ENTRY(min_entry),
		      graphics_info_t::float_to_string(min_dist).c_str());
   gtk_entry_set_text(GTK_ENTRY(max_entry),
		      graphics_info_t::float_to_string(max_dist).c_str());

   if (graphics_info_t::show_pointer_distances_flag) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
      gtk_widget_set_sensitive(frame, TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), FALSE);
      gtk_widget_set_sensitive(frame, FALSE);
   }

}

void execute_pointer_distances_settings(GtkWidget *widget) {

   GtkWidget *min_entry   = lookup_widget(widget, "pointer_distances_min_dist_entry");
   GtkWidget *max_entry   = lookup_widget(widget, "pointer_distances_max_dist_entry");
   // GtkWidget *checkbutton = lookup_widget(widget, "pointer_distances_checkbutton");

   float min_dist = 0.0;
   float max_dist = 0.0;

   float t;

   const gchar *tt = gtk_entry_get_text(GTK_ENTRY(min_entry));
   t = atof(tt);

   if ((t >= 0.0) && (t < 999.9))
      min_dist = t;

   tt = gtk_entry_get_text(GTK_ENTRY(max_entry));
   t = atof(tt);

   if ((t >= 0.0) && (t < 999.9))
      max_dist = t;

   graphics_info_t::pointer_max_dist = max_dist;
   graphics_info_t::pointer_min_dist = min_dist;

   
}


void toggle_pointer_distances_show_distances(GtkToggleButton *togglebutton) {

   GtkWidget *frame = lookup_widget(GTK_WIDGET(togglebutton),
				    "pointer_distances_frame");
   if (togglebutton->active) {
      set_show_pointer_distances(1);
      gtk_widget_set_sensitive(frame, TRUE);
   } else {
      set_show_pointer_distances(0);
      gtk_widget_set_sensitive(frame, FALSE);
   }
   
}




/*  ------------------------------------------------------------------------ */
//            model_toolbar things
/*  ------------------------------------------------------------------------ */
//

/*! \brief hide the vertical modelling toolbar in the GTK2 version */
void hide_modelling_toolbar() {
   if (graphics_info_t::use_graphics_interface_flag) { 
      GtkWidget *w = 0;
      GtkWidget *handle_box = lookup_widget(graphics_info_t::glarea,
					"model_fit_refine_toolbar_handlebox");
#if (GTK_MAJOR_VERSION > 1)
      if (graphics_info_t::model_toolbar_position_state == coot::model_toolbar::TOP ||
	  graphics_info_t::model_toolbar_position_state == coot::model_toolbar::BOTTOM) {
	w = handle_box;
      } else {
	// get the frame of the left/right toolbar
	w = gtk_widget_get_parent(handle_box);
      }
#endif // GKT_MAJOR_VERSION
      if (!w) {
	 std::cout << "failed to lookup toolbar" << std::endl;
      } else {
	 graphics_info_t::model_toolbar_show_hide_state = 0;
	 gtk_widget_hide(w);
      }
   }
}

/*! \brief show the vertical modelling toolbar in the GTK2 version
  (the toolbar is shown by default) */
void show_modelling_toolbar() {
   if (graphics_info_t::use_graphics_interface_flag) { 
      GtkWidget *w = 0;
      GtkWidget *handle_box = lookup_widget(graphics_info_t::glarea,
					    "model_fit_refine_toolbar_handlebox");

#if (GTK_MAJOR_VERSION > 1)
      if (graphics_info_t::model_toolbar_position_state == coot::model_toolbar::TOP ||
	  graphics_info_t::model_toolbar_position_state == coot::model_toolbar::BOTTOM) {
	w = handle_box;
      } else {
	w = gtk_widget_get_parent(handle_box);
      }
#endif // GKT_MAJOR_VERSION
      if (!w) {
	 std::cout << "failed to lookup toolbar" << std::endl;
      } else {
	 graphics_info_t::model_toolbar_show_hide_state = 1;
	 gtk_widget_show(w);
      }
   }
}


void
show_model_toolbar_all_icons() {

  GtkWidget *hsep           = lookup_widget(graphics_info_t::glarea,
					    "model_toolbar_hsep_toolitem2");
  GtkWidget *vsep           = lookup_widget(graphics_info_t::glarea,
					    "model_toolbar_vsep_toolitem2");
  GtkWidget *toolbar_radiobutton = lookup_widget(graphics_info_t::glarea,
						 "model_toolbar_all_icons");
    
  for (int i=0; i<(*graphics_info_t::model_toolbar_icons).size(); i++) {
    show_model_toolbar_icon(i);
  }
#if (GTK_MAJOR_VERSION >1)
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(toolbar_radiobutton), TRUE);
  if (graphics_info_t::preferences_widget) {
    // have preferences open and shall update the tree model for icon
    GtkWidget *icons_treeview = lookup_widget(graphics_info_t::preferences_widget, "preferences_model_toolbar_icon_tree");
    GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(icons_treeview));
    graphics_info_t::update_model_toolbar_icons(model);
  }
#endif // GTK2

  gtk_widget_show(hsep);
  gtk_widget_show(vsep);
}

void
show_model_toolbar_main_icons() {

  GtkWidget *hsep           = lookup_widget(graphics_info_t::glarea,
					    "model_toolbar_hsep_toolitem2");
  GtkWidget *vsep           = lookup_widget(graphics_info_t::glarea,
					    "model_toolbar_vsep_toolitem2");
  GtkWidget *toolbar_radiobutton = lookup_widget(graphics_info_t::glarea,
						 "model_toolbar_main_icons");
  
  for (int i=0; i<(*graphics_info_t::model_toolbar_icons).size(); i++) {
    if ((*graphics_info_t::model_toolbar_icons)[i].default_show_flag == 1) {
      show_model_toolbar_icon(i);
    } else {
      hide_model_toolbar_icon(i);
    }
  }
#if (GTK_MAJOR_VERSION >1)
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(toolbar_radiobutton), TRUE);
  if (graphics_info_t::preferences_widget) {
    // have preferences open and shall update the tree model for icon
    GtkWidget *icons_treeview = lookup_widget(graphics_info_t::preferences_widget, "preferences_model_toolbar_icon_tree");
    GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(icons_treeview));
    graphics_info_t::update_model_toolbar_icons(model);
  }
#endif // GTK2
  gtk_widget_hide(hsep);
  gtk_widget_hide(vsep);
}

void
update_model_toolbar_icons_menu() {

  GtkWidget *user_defined_button = lookup_widget(graphics_info_t::glarea,
						 "model_toolbar_user_defined1");
  GtkWidget *main_icons_button   = lookup_widget(graphics_info_t::glarea,
						 "model_toolbar_main_icons");
  GtkWidget *all_icons_button    = lookup_widget(graphics_info_t::glarea,
						 "model_toolbar_all_icons");

  int activate = 1;   // 0 is user defined, 1 all icons, 2 main/default icons

  for (int i=0; i<(*graphics_info_t::model_toolbar_icons).size(); i++) {
    if ((*graphics_info_t::model_toolbar_icons)[i].show_hide_flag == 0) {
      if ((*graphics_info_t::model_toolbar_icons)[i].show_hide_flag == (*graphics_info_t::model_toolbar_icons)[i].default_show_flag) {
	activate = 2;
      } else {
	activate = 0;
	break;
      }
    }
  }

  if (activate) {
    gtk_widget_hide(user_defined_button);
    if (activate == 1) {
      gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(all_icons_button), TRUE);
    } else {
      gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(main_icons_button), TRUE);
    }
  } else {
    gtk_widget_show(user_defined_button);
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(user_defined_button), TRUE);
  }

}

/*  ------------------------------------------------------------------------ */
//            popup-menu for model_toolbar
/*  ------------------------------------------------------------------------ */
//

void
toolbar_popup_menu (GtkToolbar *toolbar, 
		    GdkEventButton *event_button,
		    gpointer user_data)
{
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

#if (GTK_MAJOR_VERSION > 1)
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

  gtk_widget_show_all (menu);
  gtk_menu_popup (GTK_MENU(menu), NULL, NULL, NULL, NULL, 0, 
		  (event_button != NULL) ? event_button->time 
		  : gtk_get_current_event_time());
#endif // GTK_MAJOR_VERSION
}

void
set_model_toolbar_docked_position_callback(GtkWidget *w, gpointer user_data) {

#if (GTK_MAJOR_VERSION > 1)
  int pos = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (w), "position"));
  set_model_toolbar_docked_position(pos);
#endif //GTK_MAJOR_VERSION

}

void
reattach_modelling_toolbar() {

#if (GTK_MAJOR_VERSION > 1)
  GtkWidget *handlebox = lookup_widget(graphics_info_t::glarea,
				       "model_fit_refine_toolbar_handlebox");
  GdkEvent *event = gdk_event_new (GDK_DELETE);
  event->any.type = GDK_DELETE;
  event->any.window = GTK_HANDLE_BOX(handlebox)->float_window;
  event->any.send_event = TRUE;
  gtk_main_do_event (event);
  gdk_event_free (event);
  g_object_ref (GTK_HANDLE_BOX(handlebox)->float_window);
  g_object_ref (GTK_HANDLE_BOX(handlebox));
#endif // GTK_MAJOR_VERSION
}

/*  ------------------------------------------------------------------------ */
//            reparenting
/*  ------------------------------------------------------------------------ */
//

void
set_model_toolbar_docked_position(int state) {
  
  if (graphics_info_t::use_graphics_interface_flag) {
    GtkWidget *left_frame  = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
					   "main_window_model_fit_dialog_frame_left");
    GtkWidget *right_frame = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
					   "main_window_model_fit_dialog_frame");
    GtkWidget *handle = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
				       "model_fit_refine_toolbar_handlebox");
    GtkWidget *vbox    = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
				       "vbox1");
    GtkWidget *toolbar = lookup_widget(handle, "model_toolbar");
    GtkWidget *vsep    = lookup_widget(handle, "model_toolbar_vsep_toolitem");
    GtkWidget *hsep    = lookup_widget(handle, "model_toolbar_hsep_toolitem");
    GtkWidget *style   = lookup_widget(handle, "model_toolbar_style_toolitem");

    // reattach first, in case it wasn't and then change the mode
    if (GTK_HANDLE_BOX(handle)->child_detached) {
      reattach_modelling_toolbar();
    }

    switch (state) {

    case coot::model_toolbar::RIGHT:
      // dock to right frame
      gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar), 
				  GTK_ORIENTATION_VERTICAL);
      gtk_handle_box_set_handle_position(GTK_HANDLE_BOX(handle),
					 GTK_POS_TOP);
      gtk_widget_reparent(handle, right_frame);
      if (graphics_info_t::model_toolbar_show_hide_state) {
	gtk_widget_show(right_frame);
      }
      gtk_widget_hide(left_frame);
      graphics_info_t::model_toolbar_position_state = 0;
      gtk_widget_show(hsep);
      gtk_widget_show(style);
      gtk_widget_hide(vsep);
      break;

    case coot::model_toolbar::LEFT:
      // dock to left frame
      gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar), 
				  GTK_ORIENTATION_VERTICAL);
      gtk_handle_box_set_handle_position(GTK_HANDLE_BOX(handle),
					 GTK_POS_TOP);
      gtk_widget_reparent(handle, left_frame);
      if (graphics_info_t::model_toolbar_show_hide_state) {
	gtk_widget_show(left_frame);
      }
      gtk_widget_hide(right_frame);
      graphics_info_t::model_toolbar_position_state = 1;
      gtk_widget_show(hsep);
      gtk_widget_show(style);
      gtk_widget_hide(vsep);
      break;

    case coot::model_toolbar::TOP:
      // dock to the top
      gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar),
				  GTK_ORIENTATION_HORIZONTAL);
      gtk_handle_box_set_handle_position(GTK_HANDLE_BOX(handle),
					 GTK_POS_LEFT);
      gtk_widget_reparent(handle, vbox);
      gtk_box_set_child_packing(GTK_BOX(vbox), handle,
				FALSE, FALSE, 0, GTK_PACK_START);
      gtk_box_reorder_child(GTK_BOX(vbox), handle, 1);
      
      gtk_widget_hide(left_frame);
      gtk_widget_hide(right_frame);
      graphics_info_t::model_toolbar_position_state = 2;
      gtk_widget_hide(hsep);
      gtk_widget_hide(style);
      gtk_widget_show(vsep);
      break;

    case coot::model_toolbar::BOTTOM:
      // dock to the bottom
      gtk_toolbar_set_orientation(GTK_TOOLBAR(toolbar),
				  GTK_ORIENTATION_HORIZONTAL);
      gtk_handle_box_set_handle_position(GTK_HANDLE_BOX(handle),
					 GTK_POS_LEFT);
      gtk_widget_reparent(handle, vbox);
      gtk_box_set_child_packing(GTK_BOX(vbox), handle,
				FALSE, FALSE, 0, GTK_PACK_START);
      gtk_box_reorder_child(GTK_BOX(vbox), handle, 4);
      gtk_widget_hide(left_frame);
      gtk_widget_hide(right_frame);
      graphics_info_t::model_toolbar_position_state = 3;
      gtk_widget_hide(hsep);
      gtk_widget_hide(style);
      gtk_widget_show(vsep);
      break;

    default: 
      std::cout <<"INFO:: invalid position "<< state <<std::endl;
      break;

    }
  }
}

int suck_model_fit_dialog_bl() {

   if (graphics_info_t::use_graphics_interface_flag) { 
      GtkWidget *main_window_hbox = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
						  "main_window_hbox");
      GtkWidget *dialog = graphics_info_t::model_fit_refine_dialog;

      if (main_window_hbox) {
	 if (dialog) {
	    GtkWidget *handlebox2 = gtk_handle_box_new();
	    GtkWidget *hbox = lookup_widget(dialog, "model_fit_refine_dialog_vbox");
	    gtk_container_add (GTK_CONTAINER(handlebox2), hbox);
	    gtk_widget_reparent(hbox, handlebox2);
	    gtk_box_pack_start(GTK_BOX(main_window_hbox), handlebox2, FALSE, TRUE, 0);
	    gtk_widget_show(handlebox2);
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
   return 0;
}

int suck_model_fit_dialog() {

   if (graphics_info_t::use_graphics_interface_flag) { 
      GtkWidget *main_window_hbox = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
						  "main_window_hbox");
      GtkWidget *main_window_side_frame = lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
							"main_window_model_fit_dialog_frame");
      GtkWidget *dialog = graphics_info_t::model_fit_refine_dialog;

      if (main_window_hbox) {
	 if (dialog) {
	    GtkWidget *hbox = lookup_widget(dialog, "model_fit_refine_dialog_vbox");
	    gtk_widget_reparent(hbox, main_window_side_frame);
	    gtk_widget_show(main_window_side_frame);
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
   return 0;
}


/* Return the dialog if it exists, else null */
// Actually, this doesn't do what it says on the tin.  It closes the
// dialog_hbox and returns the dialog.
// 
GtkWidget *close_model_fit_dialog(GtkWidget *dialog_hbox) {

   GtkWidget *w = NULL;
   if (graphics_info_t::model_fit_refine_dialog_was_sucked) {
      GtkWidget *main_window_side_frame =
	 lookup_widget(GTK_WIDGET(graphics_info_t::glarea),
		       "main_window_model_fit_dialog_frame");
      gtk_widget_destroy(dialog_hbox);
      gtk_widget_hide(main_window_side_frame);
   } else {
      w = lookup_widget(GTK_WIDGET(dialog_hbox), "model_refine_dialog");
   }
   // graphics_info_t::model_fit_refine_dialog_was_sucked = 0;

   return w;
}


GtkWidget *wrapped_create_model_fit_refine_dialog() {

   GtkWidget *widget = graphics_info_t::model_fit_refine_dialog;
   if (widget) {
      
      // raise/uniconify (or whatever) what we have:
      // 
      if (!GTK_WIDGET_MAPPED(widget))
	 gtk_widget_show(widget);
      else
	 gdk_window_raise(widget->window);
   } else {
      
      // create (then store) a new one.
      
      widget = create_model_refine_dialog();
      graphics_info_t::model_fit_refine_dialog = widget;
      if (graphics_info_t::model_fit_refine_dialog_was_sucked) {
	 suck_model_fit_dialog();
      } else {
	 if (graphics_info_t::model_fit_refine_dialog_stays_on_top_flag == 1) {
	    gtk_window_set_transient_for(GTK_WINDOW(widget),
					 GTK_WINDOW(lookup_widget(graphics_info_t::glarea,
								  "window1")));
	    
	    if (graphics_info_t::model_fit_refine_x_position > -1) { 
	       gtk_widget_set_uposition(widget,
					graphics_info_t::model_fit_refine_x_position,
					graphics_info_t::model_fit_refine_y_position);
	    }
	 }
      }
   }
   
   GtkWidget *button;
   
   button = lookup_widget(widget, "model_refine_dialog_pointer_atom_button");
   if (button) {
      if (graphics_info_t::model_fit_refine_place_atom_at_pointer_string  != "")
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),
			    graphics_info_t::model_fit_refine_place_atom_at_pointer_string.c_str());
   }
   
   update_model_fit_refine_dialog_menu(widget);
   update_model_fit_refine_dialog_buttons(widget);

   graphics_info_t::set_model_fit_refine_button_names(widget);

   // Refmac button
   GtkWidget *refmac_button = lookup_widget(widget, "model_refine_dialog_refmac_button");
   if (refmac_button) {
      if (graphics_info_t::external_refinement_program_button_label != "*-*") {
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(refmac_button)->child),
			    graphics_info_t::external_refinement_program_button_label.c_str());
      }
   } 

   return widget; 
}

// function to update selected menu items (currently rot/trans only
void
update_model_fit_refine_dialog_menu(GtkWidget *dialog) {

#if (GTK_MAJOR_VERSION > 1)
   GtkWidget *menu_item;
   // update the menu for the rot/trans menu
   if (graphics_info_t::rot_trans_object_type == ROT_TRANS_TYPE_CHAIN) {
     menu_item = lookup_widget(dialog, "model_refine_dialog_rot_trans_by_chain");
     gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
   } else {
     if (graphics_info_t::rot_trans_object_type == ROT_TRANS_TYPE_MOLECULE) {
       menu_item = lookup_widget(dialog, "model_refine_dialog_rot_trans_by_molecule");
       gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
     } else {
       menu_item = lookup_widget(dialog, "model_refine_dialog_rot_trans_by_residue_range");
       gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
     }
   }
#endif
}

// function to update button labels (currently rot/trans only)
void
update_model_fit_refine_dialog_buttons(GtkWidget *dialog) {

  GtkWidget *button;
  GList *children;

  button = lookup_widget(dialog,
			 "model_refine_dialog_rot_trans_togglebutton");
  if (button) {
    if (graphics_info_t::model_fit_refine_rotate_translate_zone_string != "") {
	//	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),
#if (GTK_MAJOR_VERSION > 1)
      children = gtk_container_get_children(GTK_CONTAINER(GTK_BIN(button)->child));
#else
      children = gtk_container_children(GTK_CONTAINER(GTK_BIN(button)->child));
#endif
      children = children->next;
      gtk_label_set_text(GTK_LABEL(children->data),
			 graphics_info_t::model_fit_refine_rotate_translate_zone_string.c_str());
    }
      
  }
  graphics_info_t::set_model_fit_refine_button_names(dialog);

}


GtkWidget *wrapped_create_other_model_tools_dialog() {

   GtkWidget *widget = graphics_info_t::other_modelling_tools_dialog;
   if (!widget) {
      GtkWidget *w = create_other_model_tools_dialog();
      graphics_info_t::other_modelling_tools_dialog = w;
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
int valid_labels(const char *mtz_file_name, const char *f_col, const char *phi_col, 
		 const char *weight_col, int use_weights) {
   
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

   coot::mtz_column_types_info_t r = coot::get_f_phi_columns(mtz_file_name);

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
	    if (f_col_str == r.d_cols[i].column_label) { 
	       have_f = 1;
	       break;
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

   
   // deallocate (free) f_cols, phi_cols, weight_cols here
   //    for (int i=0; i<n_f; i++)
      //  free(f_cols[i]);
      // free(f_cols);
      //    for (int i=0; i<n_phi; i++)
      // free(phi_cols[i]);
      // free(phi_cols);
      //    for (int i=0; i<n_weight; i++)
      // free(weight_cols[i]);
      // free(weight_cols);

   if (have_f && have_phi && have_weight) 
      valid = 1;

//     std::cout << "INFO:: done checking for valid column labels... "
//  	     << valid << " " << have_f << " " << have_phi << " "
//  	     << have_weight << std::endl;
   return valid;
}

/* We need to know if an mtz file has phases.  If it doesn't then we */
/*  go down a (new 20060920) different path. */
int mtz_file_has_phases_p(const char *mtz_file_name) {

   coot::mtz_column_types_info_t r = coot::get_f_phi_columns(mtz_file_name);
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
      
      coot::mtz_column_types_info_t r = coot::get_f_phi_columns(mtz_file_name);
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
   std::cout << "cns_file_has_phases_p called " << std::endl;
   if (coot::file_exists(cns_file_name)) { 
      FILE* file = fopen( cns_file_name, "r" );
      char buf[4096];
      for ( int i = 0; i < 4096; i++ ) buf[i] = toupper(fgetc(file));
      fclose( file );
      buf[4095] = 0;
      if ( strstr( buf, "ALPHA" ) != NULL && strstr( buf, "BETA"  ) != NULL &&
	   strstr( buf, "GAMMA" ) != NULL && strstr( buf, "SYMOP" ) != NULL &&
	   strstr( buf, " F1="  ) != NULL && strstr( buf, " F2="  ) != NULL )
	 r = 1;
      else
	 r = 0; 
   } else {
      r = 0;
   }
   std::cout << "cns_file_has_phases_p returns " << r << std::endl;
   return r;
}




/*  ----------------------------------------------------------------------- */
/*                        go to atom widget                                 */
/*  ----------------------------------------------------------------------- */

// return -1 on error
//
int go_to_atom_molecule_optionmenu_active_molecule(GtkWidget *widget) { 

   graphics_info_t g;
   return g.go_to_atom_molecule_optionmenu_active_molecule(widget);
}


void save_go_to_atom_widget(GtkWidget *widget) { /* store in a static */
   graphics_info_t::go_to_atom_window = widget;
}

void unset_go_to_atom_widget() {
   graphics_info_t::go_to_atom_window = NULL;
} 


// not really a button select, its a menu item select
void
save_molecule_coords_button_select(GtkWidget *item, GtkPositionType pos) { 

   // graphics_info_t g;
   // std::cout << "INFO:: Save coords molecule now: " << pos << std::endl;
   graphics_info_t::save_imol = pos;
}




/* a c callable wrapper to the graphics_info_t function */
void fill_option_menu_with_coordinates_options(GtkWidget *option_menu, 
					       GtkSignalFunc signal_func,
					       int imol_active_position) {

   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(option_menu,
					       signal_func,
					       imol_active_position);
}

void fill_option_menu_with_coordinates_options_unsaved_first(GtkWidget *option_menu, 
							     GtkSignalFunc signal_func,
							     int imol_active_position) {
   // a mess.  This function could be deleted.
   fill_option_menu_with_coordinates_options(option_menu, signal_func, imol_active_position);
} 



// void store_refmac_params(const char *fobs_col, const char *sigfobs_col,
// 			 const char *r_free_col, int sensible_f_free_col) { 

//    graphics_info_t g;
//    int imol = g.n_molecules;
   

// } 

void free_memory_run_refmac(GtkWidget *window) {

   GtkWidget *option_menu = lookup_widget(window,
					  "run_refmac_coords_optionmenu");
   void *imol_ptr;
   GtkWidget *menu;
   GtkWidget *active_item;

   if (option_menu) {
      menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
      active_item = gtk_menu_get_active(GTK_MENU(menu));
      if (active_item) { 
	 imol_ptr = gtk_object_get_user_data(GTK_OBJECT(active_item));
      } else { 
	 std::cout << "no active item in coords option_menu\n";
      } 

      // free each of the menu items in menu
      
      // run over items in menu somehow:
   } else { 
      std::cout << "ERROR:: can't find coords option_menu in free_memory_run_refmac\n";
   } 


   option_menu = lookup_widget(window, "run_refmac_map_optionmenu");
   if (refmac_use_twin_state() == 0) {
     if (option_menu) { 
       menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
       active_item = gtk_menu_get_active(GTK_MENU(menu));
       if (active_item) { 
	 imol_ptr = gtk_object_get_user_data(GTK_OBJECT(active_item));
       } else { 
	 std::cout << "no active item in maps option_menu\n";
       }

       // free each of the menu items in menu
      
       // run over items in menu somehow:
     } else { 
       std::cout << "ERROR:: can't find map option_menu in free_memory_run_refmac\n";
     }
   }
//    std::cout << "debugging bad window got to end of free_memory_run_refmac"
// 	     << std::endl;
} 

/*  ----------------------------------------------------------------------- */
/*              new close molecule                                          */
/*  ----------------------------------------------------------------------- */
void 
new_close_molecules(GtkWidget *window) {
   
   GtkWidget *vbox = lookup_widget(window, "new_delete_molecules_vbox");
   short int closed_something_flag = 0;
   std::vector<int> closed_model_molecules;

   GtkWidget *checkbutton;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) { 
      if (graphics_info_t::molecules[imol].has_model() || 
	  graphics_info_t::molecules[imol].has_map()) { 
	 std::string button_name("delete_molecule_checkbutton_");
	 button_name += graphics_info_t::int_to_string(imol);
	 checkbutton = lookup_widget(vbox, button_name.c_str());
	 if (checkbutton) { 
	    if (GTK_TOGGLE_BUTTON(checkbutton)->active) { 
	       if (is_valid_model_molecule(imol))
		  closed_model_molecules.push_back(imol);
	       
	       // Close the sequence view of that molecule if it was
	       // displayed (fixes a bug):
	       //
	       // graphics_info_t::sequence_view_is_displayed[imol] returns the canvas
// 	       std::cout << "DEBUG:: graphics_info_t::sequence_view_is_displayed["
// 			 << imol << "] is "
// 			 << graphics_info_t::sequence_view_is_displayed[imol] << std::endl;
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
	       GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
	       if (w) {
		  // gtk_widget_destroy(graphics_info_t::sequence_view_is_displayed[imol]);
		  //
		  // use undisplay().  But how do we get to the object?
		  graphics_info_t g;
		  coot::sequence_view *seq_view = g.get_sequence_view(imol);
		  if (seq_view) { 
		     seq_view->undisplay(imol);
		  } else {
		     std::cout << "ERROR:! missing seq_view for molecule number : "
			       << imol << std::endl;
		  }
		  GtkWidget *window = lookup_widget(w, "sequence_view_dialog");
		  if (window) { 
		     gtk_widget_destroy(window);
		  } else { 
		     window = lookup_widget(w, "nsv_dialog");
		     if (window) 
			gtk_widget_destroy(window);
		  }
	       } 
#endif
	       //graphics_info_t::molecules[imol].close_yourself();
	       close_molecule(imol);
	       closed_something_flag = 1;
	    }
	 }
      }
   }


   // update go to atom molecule now that we may have deleted the
   // currently set one.
   if (closed_something_flag) {
      graphics_info_t g;
      for (unsigned int i=0; i<closed_model_molecules.size(); i++) {
	 if (closed_model_molecules[i] == g.go_to_atom_molecule()) {
	    // set it to the bottom model molecule:
	    for (int imol=graphics_info_t::n_molecules()-1; imol>=0; imol--) {
	       if (is_valid_model_molecule(imol)) {
		  g.set_go_to_atom_molecule(imol);
		  break;
	       }
	    }
	 }
      }
   }
      

   if (closed_something_flag) { 
      if (graphics_info_t::go_to_atom_window) { 
	 graphics_info_t g;
	 GtkWidget *optionmenu = lookup_widget(graphics_info_t::go_to_atom_window, 
					       "go_to_atom_molecule_optionmenu");
	 // g.fill_go_to_atom_option_menu(optionmenu);
	 int gimol = g.go_to_atom_molecule();
	 
	 GtkSignalFunc callback_func =
	    GTK_SIGNAL_FUNC(graphics_info_t::go_to_atom_mol_menu_item_select);
	 g.fill_option_menu_with_coordinates_options(optionmenu, callback_func, gimol);
      }
      graphics_draw();
   }
}

GtkWidget *wrapped_create_new_close_molecules_dialog() { 

   GtkWidget *w = create_new_close_molecules_dialog();

   GtkWidget *vbox = lookup_widget(w, "new_delete_molecules_vbox");
   GtkWidget *checkbutton;
   GtkWidget *frame;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) { 
      if (graphics_info_t::molecules[imol].has_model() || 
	  graphics_info_t::molecules[imol].has_map()) { 
	 std::string button_name("delete_molecule_checkbutton_");
	 std::string mol_name("   ");
	 mol_name += graphics_info_t::int_to_string(imol);
	 mol_name += "  ";
	 mol_name += graphics_info_t::molecules[imol].name_for_display_manager();
	 button_name += graphics_info_t::int_to_string(imol);
	 checkbutton = gtk_check_button_new_with_label(mol_name.c_str());
	 gtk_widget_ref (checkbutton);
	 gtk_object_set_data_full (GTK_OBJECT (w), 
				   button_name.c_str(), checkbutton,
				   (GtkDestroyNotify) gtk_widget_unref);
	 gtk_widget_show (checkbutton);
	 frame = gtk_frame_new(NULL);
	 gtk_container_add (GTK_CONTAINER(frame), checkbutton);

	 gtk_box_pack_start (GTK_BOX (vbox), frame, FALSE, FALSE, 0);
	 gtk_container_set_border_width (GTK_CONTAINER (frame), 2);
	 gtk_widget_show (frame);
      }
   }

   return w;
} 

/*  ----------------------------------------------------------------------- */
/*              old close molecule                                          */
/*  ----------------------------------------------------------------------- */
/* get the molecule to delete from the optionmenu */
void
close_molecule_by_widget(GtkWidget *optionmenu) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
   GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));

   if (active_item) { 
      int *imol_p = (int *) gtk_object_get_user_data(GTK_OBJECT(active_item));

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
}


void close_molecule(int imol) {


   graphics_info_t g;
   int old_go_to_atom_molecule = g.go_to_atom_molecule();
   if (is_valid_model_molecule(imol) ||
       is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].close_yourself();
   }
   int go_to_atom_imol_new = g.update_go_to_atom_molecule_on_go_to_atom_molecule_deleted();
   if (graphics_info_t::go_to_atom_window) {
      // std::cout << ".....re fill go to atom window here" << std::endl;
      if (imol == old_go_to_atom_molecule) {
	 g.update_go_to_atom_window_on_other_molecule_chosen(go_to_atom_imol_new);
	 g.update_go_to_atom_window_on_changed_mol(go_to_atom_imol_new);
      } 
   }
   graphics_draw();
   std::string cmd = "close-molecule";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);

}

void fill_close_option_menu_with_all_molecule_options(GtkWidget *optionmenu) {

   
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
	  g.molecules[imol].xmap_is_filled[0]) {
	 char s[200];
	 snprintf(s,199,"%d",imol);
	 std::string ss(s);
	 ss += " ";
	 ss += g.molecules[imol].name_;
	 // char *txt = (char *)malloc(ss.length()+1);
	 // snprintf(txt, 3, "%d", imol); 
	 menuitem = gtk_menu_item_new_with_label (ss.c_str());
	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     GTK_SIGNAL_FUNC(close_molecule_item_select),
			     GINT_TO_POINTER(imol));
	 int *ip = new int;
	 *ip = imol;
	 gtk_object_set_user_data(GTK_OBJECT(menuitem), ip );
	 gtk_menu_append(GTK_MENU(menu), menuitem); 
	 gtk_widget_show(menuitem); 
      }
   }

   if (g.n_molecules() > 0) {
      gtk_menu_set_active(GTK_MENU(menu), 0);
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu),
			    menu);
   
}



void
close_molecule_item_select(GtkWidget *item, GtkPositionType pos) {

   std::cout << "activating closing position/imol " << pos << std::endl;
   
}

void add_ccp4i_project_shortcut(GtkWidget *fileselection) {

#if (GTK_MAJOR_VERSION > 1)
   // Paul likes to have a current dir shortcut, here we go then:
   gchar *current_dir = g_get_current_dir();
   gtk_file_chooser_add_shortcut_folder(GTK_FILE_CHOOSER(fileselection),
					current_dir,
					NULL);
   g_free(current_dir);
//    std::cout << "DEBUG:: adding a short cut..." << std::endl;
//    std::cout << "DEBUG:: widget is filechooser: " << GTK_IS_FILE_CHOOSER(fileselection) << std::endl;
   // BL says: we simply add a short cut to ccp4 project folder
   // based on ccp4_defs_file_name()
   // add all projects to shortcut (the easiest option for now)
   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   
   std::vector<std::pair<std::string, std::string> > project_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);
   
   for (unsigned int i=0; i<project_pairs.size(); i++) {
      const char *folder = project_pairs[i].second.c_str();
      int len = strlen(folder);
      if (len > 0) {
	 gtk_file_chooser_add_shortcut_folder(GTK_FILE_CHOOSER(fileselection),
					      project_pairs[i].second.c_str(),
					      NULL);
      }
   }
#endif   // GTK_MAJOR_VERSION
}

void add_ccp4i_project_optionmenu(GtkWidget *fileselection, int file_selection_type) {

   bool add_shortcut = 0;

#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      add_shortcut = 1;
   }
#endif

   if (add_shortcut) {
//       std::cout << "in add_ccp4i_project_optionmenu widget is fileselection "
// 		<< GTK_IS_FILE_CHOOSER(fileselection) << std::endl;
      add_ccp4i_project_shortcut(fileselection);

   } else {

      GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;

      GtkWidget *optionmenu = gtk_option_menu_new();
      gtk_widget_ref(optionmenu);
      gtk_widget_show(optionmenu);
      gtk_object_set_data(GTK_OBJECT(fileselection), "ccp4i_project_optionmenu", optionmenu);
      gtk_object_set_user_data(GTK_OBJECT(optionmenu), GINT_TO_POINTER(file_selection_type));
      GtkSignalFunc project_signal_func =
	 GTK_SIGNAL_FUNC(option_menu_refmac_ccp4i_project_signal_func);
      add_ccp4i_projects_to_optionmenu(optionmenu, file_selection_type, project_signal_func);

      // Let's put the optionmenu in a frame with a label
      GtkWidget *frame = gtk_frame_new("CCP4i Project Directory");
      gtk_container_add(GTK_CONTAINER(aa), frame);
      gtk_widget_show(frame);
      gtk_container_add(GTK_CONTAINER(frame),optionmenu);
   }
}

void add_ccp4i_projects_to_optionmenu(GtkWidget *optionmenu,
				      int file_selection_type,
				      GtkSignalFunc func) {

   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();

   std::vector<std::pair<std::string, std::string> > project_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);

   // project_pairs.push_back(std::pair<std::string, std::string> ("thingy", "X"));

   GtkWidget *menu = gtk_menu_new();
   
   for (unsigned int i=0; i<project_pairs.size(); i++) {
      GtkWidget *menu_item = gtk_menu_item_new_with_label(project_pairs[i].first.c_str());
      gtk_widget_show(menu_item);
      gtk_menu_append(GTK_MENU(menu), menu_item);
      int i_info = file_selection_type;
      i_info <<= 10;
      i_info += i;
      gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
			 GTK_SIGNAL_FUNC(func),
			 GINT_TO_POINTER(i_info));
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu), menu);

   // set the active menu item here...
   
   gtk_widget_show(menu);
}

// This happens when someone actives the option menu of the CCP4i
// project in file selectors.
void
option_menu_refmac_ccp4i_project_signal_func(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t g;
   int file_selection_type =  pos >> 10;
   int idx = pos - (file_selection_type << 10);
   g.ccp4_projects_index_last = idx;
   
   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > pr_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);
   if (idx < int(pr_pairs.size())) {
      g.set_directory_for_fileselection_string(pr_pairs[idx].second);

      // we need to call set_directory_for_fileselection with an
      // argument that is the fileselection widget.  The question is:
      // what is that widget?
      //
      // Here, we don't know, so we try to look up each of the
      // fileselection widgets
      //
      GtkWidget *optionmenu = lookup_widget(item, "ccp4i_project_optionmenu");

      if (! optionmenu) {
	 std::cout << "WARNING:: failed to find ccp4i optionmenu in "
		   << "option_menu_refmac_ccp4i_project_signal_func" << std::endl;
      } else { 

	 GtkWidget *fileselection;
	 fileselection = lookup_file_selection_widgets(item, file_selection_type);
	 
	 if (fileselection) {
	    g.set_directory_for_fileselection(fileselection);
	 } else {
	    std::cout << "WARNING:: failed to find filesection in "
		      << "option_menu_refmac_ccp4i_project_signal_func" << std::endl;
	 }
      }
   } else {
      std::cout << "ERROR:: error in indexing in option_menu_refmac_ccp4i_project_signal_func"
		<< std::endl;
   }
}

void
run_refmac_ccp4i_option_menu_signal_func(GtkWidget *item, GtkPositionType pos) {

   // graphics_info_t g;
   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > pr_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);

   if (pos < int(pr_pairs.size())) {
      std::cout << "CCP4i Project Directory for refmac: " 
		<< pr_pairs[pos].second << std::endl;
      graphics_info_t::refmac_ccp4i_project_dir = pr_pairs[pos].second;

   } else {
      std::cout << "ERROR:: error in indexing in run_refmac_ccp4i_option_menu_signal_func"
		<< std::endl;
   }
}

void clear_refmac_ccp4i_project() { 

   graphics_info_t::refmac_ccp4i_project_dir = std::string("");

} 

void add_filename_filter(GtkWidget *fileselection) { 

  bool add_filter = 1;
#if (GTK_MAJOR_VERSION > 1)
  if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      add_filter = 0;
      // maybe we use it to set the scripting filter for now!
      // not general enough but may do for now!
      add_filechooser_filter_button(fileselection, COOT_SCRIPTS_FILE_SELECTION);
  }
#endif
  if (add_filter) {
   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;

   GtkWidget *frame = gtk_frame_new("File-name filter:");
   GtkWidget *entry = gtk_entry_new();
   gtk_widget_set_usize (entry, 60, -2);

   gtk_widget_ref(entry);
   gtk_widget_show(entry);
   gtk_widget_ref(frame);

   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), entry);
   // I want the entry to be not-expandable:  How do I do that?
   // gtk_box_pack_start (GTK_BOX (hbox), entry, FALSE, TRUE, 0);

   gtk_signal_connect (GTK_OBJECT (entry), "key_press_event",
		       GTK_SIGNAL_FUNC (on_filename_filter_key_press_event),
		       NULL);
   gtk_widget_show(frame);
  }

}


// 
GtkWidget *lookup_file_selection_widgets(GtkWidget *item, int file_selection_type) {

   GtkWidget *w = 0;
   if (file_selection_type == COOT_COORDS_FILE_SELECTION)
      w = lookup_widget(GTK_WIDGET(item), "coords_fileselection1");

   if (file_selection_type == COOT_DATASET_FILE_SELECTION) 
      w = lookup_widget(GTK_WIDGET(item), "dataset_fileselection1");

   if (file_selection_type == COOT_MAP_FILE_SELECTION) 
      w = lookup_widget(GTK_WIDGET(item), "map_name_fileselection1");

   // phs_coordinates_fileselection doesn't have a filter button yet.
   if (file_selection_type == COOT_PHS_COORDS_FILE_SELECTION) 
      w = lookup_widget(GTK_WIDGET(item), "phs_coordinates_fileselection");

   if (file_selection_type == COOT_SAVE_COORDS_FILE_SELECTION) 
      w = lookup_widget(GTK_WIDGET(item), "save_coords_fileselection1");

   if (file_selection_type == COOT_CIF_DICTIONARY_FILE_SELECTION) 
      w = lookup_widget(GTK_WIDGET(item), "cif_dictionary_fileselection");

   if (file_selection_type == COOT_SCRIPTS_FILE_SELECTION) 
      w = lookup_widget(GTK_WIDGET(item), "run_script_fileselection");

   return w;
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
   snprintf(text,100,"%-5.1f", g.box_radius);

   return text;

}
GtkWidget *wrapped_create_show_symmetry_window() {

   GtkWidget *show_symm_window = create_show_symmetry_window();
   GtkWidget *checkbutton; 
   GtkButton       *button; 
   

   /* Colour Merge */
   GtkAdjustment *adjustment;
   GtkScale *hscale;

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

   if (get_show_symmetry() == 1) { 
				    
      button = GTK_BUTTON(lookup_widget(show_symm_window,
					"show_symmetry_yes_radiobutton"));
   } else { 
      button = GTK_BUTTON(lookup_widget(show_symm_window,
					"show_symmetry_no_radiobutton"));
   }
      
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

/* Show Symmetry as Calphas checkbutton */

//    checkbutton = lookup_widget(GTK_WIDGET(button),
// 			       "show_symmetry_as_calphas_checkbutton");
//    if (get_symmetry_as_calphas_state()) { 
//      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
//    }
   
/* The Colour Merge hscale */

   hscale = GTK_SCALE(lookup_widget(show_symm_window,
				    "hscale_symmetry_colour"));
   
   adjustment = GTK_ADJUSTMENT 
      (gtk_adjustment_new(0.5, 0.0, 3.0, 0.02, 0.05, 2.0));

   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   gtk_signal_connect (GTK_OBJECT (adjustment), "value_changed",
		       GTK_SIGNAL_FUNC (symmetry_colour_adjustment_changed), 
		       NULL);

/*  The Symmetry Search Radius Entry */

    entry = lookup_widget(show_symm_window, "symmetry_radius_entry");
    
    text = get_text_for_symmetry_size_widget(); /* const gchar *text */
    gtk_entry_set_text(GTK_ENTRY(entry), text);

    free (text); 

/* The Unit Cell Radiobuttons */

    if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) { 
       if (get_show_unit_cell(imol) == 1) { 
	  button = GTK_BUTTON(lookup_widget(show_symm_window,
					    "unit_cell_yes_radiobutton"));
       } else { 
	  button = GTK_BUTTON(lookup_widget(show_symm_window,
					    "unit_cell_no_radiobutton"));
       }
       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    } 

       
    //  The Expanded Atoms Label checkbutton

    checkbutton = lookup_widget(show_symm_window,
				"show_symmetry_expanded_labels_checkbutton");

    if (graphics_info_t::symmetry_atom_labels_expanded_flag) {
       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
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


   // does a graphics_draw() for us...
   set_symmetry_colour_merge(adj->value); /* this adjusts
				             graphics_info_t::
					     symm_colour_merge_weight,
					     which is a double
					     array.  But we only
					     ever use the 0th
					     position of it in
					     combine_colour() */
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
handle_map_colour_change(int imol, gdouble* map_col) {
   
   graphics_info_t::molecules[imol].handle_map_colour_change(map_col,
							    graphics_info_t::swap_difference_map_colours,
							    graphics_info_t::GL_CONTEXT_MAIN);

   if (graphics_info_t::display_mode_use_secondary_p()) {
      graphics_info_t g;
      g.make_gl_context_current(graphics_info_t::GL_CONTEXT_SECONDARY);
      g.molecules[imol].handle_map_colour_change(map_col,
						 g.swap_difference_map_colours,
						 graphics_info_t::GL_CONTEXT_SECONDARY);
      g.make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
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
      double *colours = new double[4];
      colours[0] = red;
      colours[1] = green;
      colours[2] = blue;
      short int swap_col = graphics_info_t::swap_difference_map_colours;
      graphics_info_t::molecules[imol].handle_map_colour_change(colours, swap_col,
								graphics_info_t::GL_CONTEXT_MAIN);
      if (graphics_info_t::display_mode_use_secondary_p()) {
	 graphics_info_t g;
	 g.make_gl_context_current(graphics_info_t::GL_CONTEXT_SECONDARY);
	 graphics_info_t::molecules[imol].handle_map_colour_change(colours, swap_col,
								   graphics_info_t::GL_CONTEXT_SECONDARY);
	 g.make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
      } 
      delete [] colours;
   }
}




void add_on_map_colour_choices(GtkWidget *menu) {

   GtkSignalFunc callback = GTK_SIGNAL_FUNC(map_colour_mol_selector_activate);
   std::string sub_menu_name = "map_colour1_menu";
   GtkWidget *sub_menu = lookup_widget(menu, sub_menu_name.c_str());
   if (!sub_menu) {
      std::cout << "ERROR: sub menu not found in add_on_map_colour_choices\n";
   } else {
      gtk_container_foreach(GTK_CONTAINER(sub_menu),
			    my_delete_menu_items,
			    (gpointer) sub_menu);
      for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
	 if (graphics_info_t::molecules[imol].has_map()) {
	    std::string name;
	    name = graphics_info_t::molecules[imol].dotted_chopped_name();
	    add_map_colour_mol_menu_item(imol, name, sub_menu, callback);
	 }
      }
   }
}

void
add_map_colour_mol_menu_item(int imol, const std::string &name,
			     GtkWidget *menu, GtkSignalFunc callback) {

   int *imol_data = new int;
   *imol_data = imol;
   GtkWidget *menu_item = gtk_menu_item_new_with_label(name.c_str());
   gtk_container_add(GTK_CONTAINER(menu), menu_item);
   gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
		      callback, (gpointer) imol_data);
   gtk_widget_show(menu_item);

}

void my_delete_menu_items(GtkWidget *widget, void *data) {
   gtk_container_remove(GTK_CONTAINER(data), widget);

}


void map_colour_mol_selector_activate (GtkMenuItem     *menuitem,
				       gpointer         user_data) {

   GtkWidget *col_sel_window;
   GtkWidget  *colorseldlg;
   GtkColorSelection *colorsel;
   gdouble *colour;

   struct map_colour_data_type *map_colour_data; 
   map_colour_data = (struct map_colour_data_type *) user_data; 

   col_sel_window = create_map_colour_selection_window(map_colour_data); 
   colorseldlg = GTK_WIDGET(lookup_widget(col_sel_window, "map_colour_selection"));
   colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(colorseldlg)->colorsel);   

   colour = get_map_colour(map_colour_data->imol);
   gtk_color_selection_set_color(colorsel, colour);
   gtk_widget_show(col_sel_window); 
   free(colour); 
}

// ---------------------------------------------------------
// Scroll wheel, similar
// ---------------------------------------------------------
//

int scroll_wheel_map() {
   return graphics_info_t::scroll_wheel_map;
}

void set_scroll_wheel_map(int imap) {

   if (is_valid_map_molecule(imap)) {
      graphics_info_t::scroll_wheel_map = imap;
   }
}


void add_on_map_scroll_whell_choices(GtkWidget *menu) {

   GtkSignalFunc callback =
      GTK_SIGNAL_FUNC(map_scroll_wheel_mol_selector_activate);
   //    std::string sub_menu_name = "attach_scroll_wheel_to_which_map_1";
   // std::string sub_menu_name = "mapscroll_wheelmap1";
   std::string sub_menu_name = "map_scroll_wheel_menu";
   GtkWidget *sub_menu = lookup_widget(menu, sub_menu_name.c_str());
   if (!sub_menu) {
      std::cout << "ERROR: sub menu not found in add_on_map_scroll_whell_choices\n";
   } else {
      gtk_container_foreach(GTK_CONTAINER(sub_menu),
			    my_delete_menu_items,
			    (gpointer) sub_menu);
      for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
	 if (graphics_info_t::molecules[imol].has_map()) {
	    std::string name;
	    name = graphics_info_t::molecules[imol].dotted_chopped_name();
	    add_map_scroll_wheel_mol_menu_item(imol, name, sub_menu, callback);
	 }
      }
   }
}


void
map_scroll_wheel_mol_selector_activate (GtkMenuItem     *menuitem,
					gpointer         user_data)
{
   int *imol = (int *) user_data;
   std::cout << "INFO:: scroll map change to map number "
	     << *imol << std::endl;

   graphics_info_t g;
   g.scroll_wheel_map = *imol; // needs to be set before we can active
			       // the button.
   g.activate_scroll_radio_button_in_display_manager(*imol);
   set_scrollable_map(*imol);  // 0 is ignored.
}

void
add_map_scroll_wheel_mol_menu_item(int imol, const std::string &name,
				    GtkWidget *menu, GtkSignalFunc callback) {

   int *imol_data = new int;
   *imol_data = imol;
   GtkWidget *menu_item = gtk_menu_item_new_with_label(name.c_str());
   gtk_container_add(GTK_CONTAINER(menu), menu_item);
   gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
		      callback, (gpointer) imol_data);
   gtk_widget_show(menu_item);

}

GtkWidget *wrapped_create_bond_parameters_dialog() {

   GtkWidget *widget = create_bond_parameters_dialog();

   GtkWidget *optionmenu =
      lookup_widget(widget, "bond_parameters_molecule_optionmenu");

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(graphics_info_t::bond_parameters_molecule_menu_item_select);


   // fill the colour map rotation entry

   // check the Carbons only check button

   // fill the molecule bond width option menu

   // check the draw hydrogens check button

   // Consider the case were we set the bond_parameters_molecule and
   // then close that molecule Usually, we want imol to be set to
   // bond_parameters_molecule but in the case of a closed molecule we
   // want bond_parameters_molecule to be set to first_coords_imol()
   // i.e. imol.

   graphics_info_t g;
   int imol = first_coords_imol(); // can be -1;
   if (g.bond_parameters_molecule >= 0)
      if (g.molecules[g.bond_parameters_molecule].has_model())
	 imol = g.bond_parameters_molecule;
      else
	 g.bond_parameters_molecule = imol;
   else
      // g.bond_parameters_molecule not set yet.
      g.bond_parameters_molecule = imol;

   std::cout << "DEBUG:: in wrapped_create_bond_parameters_dialog imol is "
 	     << imol << " g.bond_parameters_molecule "
	     << g.bond_parameters_molecule << std::endl;

   g.fill_option_menu_with_coordinates_options(optionmenu, callback_func, imol);
   graphics_info_t::fill_bond_parameters_internals(widget, imol);

   return widget;
}

void apply_bond_parameters(GtkWidget *w) {

   graphics_info_t g;
   int imol = g.bond_parameters_molecule;

   
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

	    GtkWidget *toggle_button = lookup_widget(w, "draw_hydrogens_yes_radiobutton");
	    if (GTK_TOGGLE_BUTTON(toggle_button)->active) {
	       set_draw_hydrogens(imol, 1);
	    } else {
	       set_draw_hydrogens(imol, 0);
	    }
	 }
      }
   }
   graphics_draw();
}

void skeletonize_map_by_optionmenu(GtkWidget *optionmenu) { 

   GtkWidget *window = lookup_widget(GTK_WIDGET(optionmenu), "skeleton_dialog");

   GtkWidget *on_radio_button;
   GtkWidget *prune_check_button;

   on_radio_button = lookup_widget(window, "skeleton_on_radiobutton");

   short int do_it = 0; 
   short int prune_it = 0;
   if (! is_valid_map_molecule(graphics_info_t::map_for_skeletonize)) {
      std::cout << "ERROR:: Trapped a bad map for skeletoning!" << std::endl;
   } else {
      if (GTK_TOGGLE_BUTTON(on_radio_button)->active) { 
	 do_it = 1;
      }
      prune_check_button = lookup_widget(window,"skeleton_prune_and_colour_checkbutton");
      if (GTK_TOGGLE_BUTTON(prune_check_button)->active) { 
	 prune_it = 1;
      }

      if (do_it)
	 graphics_info_t::skeletonize_map(prune_it, graphics_info_t::map_for_skeletonize);
      else {
	 std::cout << "INFO:: unskeletonizing map number "
		   << graphics_info_t::map_for_skeletonize << std::endl;
	 graphics_info_t::unskeletonize_map(graphics_info_t::map_for_skeletonize);
      }
   }
} 

void
skeletonize_map_single_map_maybe(GtkWidget *window, int imol) { 
   GtkWidget *on_radio_button = 
      lookup_widget(window, "single_map_skeleton_on_radiobutton");

   if (GTK_TOGGLE_BUTTON(on_radio_button)->active) { 
      graphics_info_t::skeletonize_map(0, imol);
      if (graphics_info_t::map_for_skeletonize < 0) {
	 // it was unset, so set it...
	 graphics_info_t::map_for_skeletonize = imol;
      }
   } else { 
      graphics_info_t::unskeletonize_map(imol);
   } 
}

void set_file_for_save_fileselection(GtkWidget *fileselection) { 

   graphics_info_t g;
   bool no_chooser = 1;
#if (GTK_MAJOR_VERSION > 1)
   if (g.gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      no_chooser = 0;
      g.set_file_for_save_filechooser(fileselection);
   }
#endif
   if (no_chooser) {
      g.set_file_for_save_fileselection(fileselection);
   }
}



GtkWidget *wrapped_create_skeleton_dialog() { 

   graphics_info_t g;
   return g.wrapped_create_skeleton_dialog();
}


void save_coordinates_using_widget(GtkWidget *widget) {

   // the widget that we get passed is the fileselection widget

   char *stuff =  (char *) gtk_object_get_user_data(GTK_OBJECT(widget));

   if (! stuff) {

      std::cout << "Ooops no data associated with that widget - "
		<< " no molecules with coordinates?" << std::endl;

   } else { 

      int imol = *((int *) stuff);

      // How do we get the filename?

      const gchar *filename;
#if (GTK_MAJOR_VERSION > 1)
      if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE)  {
      	filename = gtk_file_chooser_get_filename
	 (GTK_FILE_CHOOSER(widget));
      } else {
	filename = gtk_file_selection_get_filename
	 (GTK_FILE_SELECTION(widget));
     }
#else
      filename = gtk_file_selection_get_filename
	 (GTK_FILE_SELECTION(widget));
#endif

      std::cout << "save coordinates for molecule "
		<< imol << " to file " << filename << std::endl;

      graphics_info_t g;
      int ierr = save_coordinates(imol, filename);
      if (! ierr) { 
	 std::string s = "Saved coordinates file ";
	 s += filename;
	 s += ".";
	 g.statusbar_text(s);
      }
   }
}

void save_symmetry_coords_from_fileselection(GtkWidget *fileselection) {

   coot::Symm_Atom_Pick_Info_t *symm_info =
      (coot::Symm_Atom_Pick_Info_t *) gtk_object_get_user_data(GTK_OBJECT(fileselection));

   const gchar *filename;
#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
	 filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));
   } else {
	 filename = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fileselection));
   }
#else
   filename = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fileselection));
#endif // GTK_MAJOR_VERSION

   if (symm_info) {
      // std::cout << "Preshift to origin:  " << symm_info->pre_shift_to_origin << std::endl;
      save_symmetry_coords(symm_info->imol,
			   filename,
			   symm_info->symm_trans.isym(),
			   symm_info->symm_trans.x(), 
			   symm_info->symm_trans.y(), 
			   symm_info->symm_trans.z(), 
			   symm_info->pre_shift_to_origin.us,			   
			   symm_info->pre_shift_to_origin.vs,
			   symm_info->pre_shift_to_origin.ws);
   } else {
      std::cout << "ERROR:: failed to get user data from save symmetry coords fileselection"
		<< std::endl;
      std::cout << "ERROR:: saving of symmetry coordinates failed" << std::endl;
   }
}


GtkWidget *wrapped_create_goto_atom_window() {

   GtkWidget *widget = graphics_info_t::go_to_atom_window; 
   if (widget) {
      if (!GTK_WIDGET_MAPPED(widget))
	 gtk_widget_show(widget);
      else
	 gdk_window_raise(widget->window);
   } else {
      widget = create_goto_atom_window();
      graphics_info_t::go_to_atom_window = widget;
      if (graphics_info_t::go_to_atom_window_x_position > -1) { 
	 gtk_widget_set_uposition(widget,
				  graphics_info_t::go_to_atom_window_x_position,
				  graphics_info_t::go_to_atom_window_y_position);
      }
      fill_go_to_atom_window(widget);
   } 
   return widget;
} 


void post_go_to_atom_window() {

   std::cout << "DEBUG:: post_go_to_atom_window()" << std::endl;
   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *widget = wrapped_create_goto_atom_window();
      gtk_widget_show(widget);
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("post-go-to-atom-window");
   add_to_history(command_strings);
}

void fill_go_to_atom_window(GtkWidget *widget) {

     GtkWidget *option_menu; 
     GtkWidget *chain_entry;
     GtkWidget *residue_entry; 
     GtkWidget *atom_name_entry; 
     GtkWidget *residue_gtklist;
     gchar *text; 
     GtkWidget *scrolled_window;

/* First lets do the molecule optionmenu. */

     graphics_info_t g;
     GtkSignalFunc callback_func =
      GTK_SIGNAL_FUNC(graphics_info_t::go_to_atom_mol_menu_item_select);
     option_menu = lookup_widget(GTK_WIDGET(widget), 
				 "go_to_atom_molecule_optionmenu");
     
     /* These are in a special order: The residue is done first
	because it is set to a magic number (-9999 (or so)) initially.
	In that case, we do magic in
	get_text_for_go_to_atom_residue_entry(), i.e. look up a real
	atom of a molecule and set also the go to chain and the go to
	atom name  */

     /* The residue entry */

     residue_entry = lookup_widget(GTK_WIDGET(widget),
				   "go_to_atom_residue_entry"); 


     // text = get_text_for_go_to_atom_residue_entry();  // old

     // on startup, tinkers with
     // go to atom params, yuck,
     // I think.
     std::string rt = coot::util::int_to_string(g.go_to_atom_residue());

     gtk_entry_set_text(GTK_ENTRY(residue_entry), rt.c_str()); 

     /* Now that the go to atom molecule has been set, we can use it
	to fill the molecule option menu */
     int gimol = g.go_to_atom_molecule();
     g.fill_option_menu_with_coordinates_options(option_menu,
						 callback_func,
						 gimol);

     /* The chain entry */

     chain_entry = lookup_widget(GTK_WIDGET(widget),
				 "go_to_atom_chain_entry");
     
     // text = get_text_for_go_to_atom_chain_entry(); 
     gtk_entry_set_text(GTK_ENTRY(chain_entry), g.go_to_atom_chain()); 

     
     /* The Atom Name entry */

     atom_name_entry = lookup_widget(GTK_WIDGET(widget),
				     "go_to_atom_atom_name_entry"); 
     // text = get_text_for_go_to_atom_atom_name_entry(); 
     gtk_entry_set_text(GTK_ENTRY(atom_name_entry), g.go_to_atom_atom_name()); 

     /* The Residue List */

     /* The residue list cant be added to a scrolled window in glade,
	so we create only a scrolled window in glade
	(go_to_atom_residue_scrolledwindow) and add the list to it
	like is done in examples/list/list.c */

     scrolled_window = lookup_widget(GTK_WIDGET(widget),
				     "go_to_atom_residue_scrolledwindow");
     residue_gtklist=gtk_list_new();

#if (GTK_MAJOR_VERSION == 1)

     GtkWidget *residue_tree = gtk_tree_new();
     gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
					   residue_tree);
     gtk_tree_set_selection_mode (GTK_TREE(residue_tree),
				  GTK_SELECTION_SINGLE);
     gtk_widget_show(residue_tree);

     // now set the adjustment of the viewport/scrolledwindow to the
     // gtklist for residue
     // 
     // This bit of magic took 2 full days to find and is necessary
     // for well-formed ajustments on the residue list, which means
     // that gtk_list_scroll_vertical GTK_SCROLL_JUMP works! (See
     // make_synthetic_select_on_residue_list). The $1000 feature.
     // Jan, I hope you are satisfied - it was a struggle.

//      gtk_container_set_focus_vadjustment(GTK_CONTAINER(residue_tree), 
// 					 gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW (scrolled_window)));
     
     gtk_widget_ref(residue_tree);
     gtk_object_set_data_full(GTK_OBJECT(widget), "go_to_atom_residue_tree",
			      residue_tree, 
			      (GtkDestroyNotify) gtk_widget_unref);

     gtk_signal_connect(GTK_OBJECT(residue_tree),
  			"selection_changed",
  			GTK_SIGNAL_FUNC(on_go_to_atom_residue_tree_selection_changed_gtk1),
  			NULL);
     
     /* The atom list */
     scrolled_window = lookup_widget(GTK_WIDGET(widget),
				     "go_to_atom_atom_scrolledwindow");
     GtkWidget *atom_gtklist=gtk_list_new();
     gtk_scrolled_window_add_with_viewport( GTK_SCROLLED_WINDOW(scrolled_window),
					    atom_gtklist);
     /* attach the name to the widget (by hand (as interface.c does
	it) so that we can look it up in the callback of residue selection changed */
     gtk_widget_ref(atom_gtklist);
     gtk_object_set_data_full(GTK_OBJECT(widget), "go_to_atom_atom_list", 
			      atom_gtklist, 
			      (GtkDestroyNotify) gtk_widget_unref);

     gtk_widget_show(atom_gtklist);


      gtk_signal_connect(GTK_OBJECT(atom_gtklist),
 			"selection_changed",
 			GTK_SIGNAL_FUNC(on_go_to_atom_atom_list_selection_changed_gtk1),
 			NULL);

     /* fill those atom and residue lists (which uses
	graphics_info_t::go_to_atom_residue()) */
      g.fill_go_to_atom_residue_list_gtk1(residue_tree);

#else
     // -----------------------------------------------------------------
     //                GTK2 path
     // -----------------------------------------------------------------
     GtkWidget *atom_list_scrolled_window =
	lookup_widget(GTK_WIDGET(widget), "go_to_atom_atom_scrolledwindow");
     g.fill_go_to_atom_window_gtk2(widget, // the go to atom window
				   scrolled_window,
				   atom_list_scrolled_window);

#endif      

     /* store the widget */
     save_go_to_atom_widget(widget);

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

  entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_chain_entry"));
  chain_str = gtk_entry_get_text(entry);

  entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_residue_entry"));
  res_str = gtk_entry_get_text(entry);

  entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_atom_name_entry"));
  atom_name_str = gtk_entry_get_text(entry);

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


#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

void on_go_to_atom_residue_tree_selection_changed_gtk1 (GtkList         *gtktree,
							gpointer         user_data) {
   
   graphics_info_t g;
   g.on_go_to_atom_residue_tree_selection_changed_gtk1(gtktree, user_data);
}
#endif


void clear_atom_list(GtkWidget *atom_gtklist) {

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
   gtk_list_clear_items(GTK_LIST(atom_gtklist), 0, -1);
#else
   // FILL ME
#endif 
}

void on_go_to_atom_residue_list_select_child (GtkList         *list,
					      GtkWidget       *widget,
					      gpointer         user_data) {
   std::cout << "child selected.\n";
}

void on_go_to_atom_residue_list_unselect_child (GtkList         *list,
						GtkWidget       *widget,
						gpointer         user_data) {
   std::cout << "child unselected.\n"; 
}


#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

void on_go_to_atom_atom_list_selection_changed_gtk1(GtkList         *atom_gtklist,
					       gpointer         user_data) {

   graphics_info_t g;
   GList   *dlist;
   coot::model_view_atom_button_info_t *item_data;
   
   /* Fetch the doubly linked list of selected items
    * of the List, remember to treat this as read-only!
    */
   dlist=GTK_LIST(atom_gtklist)->selection;
   /* If there are no selected items there is nothing more
    * to do than just telling the user so
    */
   if (!dlist) {
      // g_print("Selection cleared\n");
      return;
   }  
   /* Ok, we got a selection and so we print it
    */
   // std::cout << "The selection is a ";
    
   /* Get the list item from the doubly linked list
    * and then query the data associated with list_item_data_key.
    * We then just print it */
   while (dlist) {
      GtkObject       *list_item;
      // gchar           *item_data_string;

      if (! dlist->data) 
	 std::cout << "ERROR: (coot) on_go_to_atom_atom_list_selection_changed: "
		   << "no dlist->data!\n";
      
      list_item=GTK_OBJECT(dlist->data);
//       item_data_string=gtk_object_get_data(list_item,
// 					   list_item_data_key);

      if (list_item) { 
	 item_data = (coot::model_view_atom_button_info_t *) gtk_object_get_user_data(list_item);
	 //       std::cout << "DEBUG:: buton_label " << item_data->button_label << std::endl;

	 if (item_data) { 
	    CAtom *at = item_data->atom;

	    // Try looking up at in g.go_to_atom_molecule() here?
	    if (at) {
	       g.set_go_to_atom_chain_residue_atom_name(at->GetChainID(),
							at->GetSeqNum(),
							at->name,
							at->altLoc);
	       g.update_widget_go_to_atom_values(g.go_to_atom_window, at);
	    }
	 } else { 
	    std::cout << "ERROR:: Oops! Can't get item data in "
		      << "on_go_to_atom_atom_list_selection_changed\n";
	 }
      } else { 
	    std::cout << "ERROR:: Oops! NULL list_item in "
		      << "on_go_to_atom_atom_list_selection_changed\n";
      }
      dlist=dlist->next;
   }
   // g_print("\n");

}

#endif // GTK version


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
// 	 gtk_widget_show(widget);
//       else
// 	 gdk_window_raise(widget->window);
      
      gtk_widget_destroy(widget);
      widget = create_residue_info_dialog();
      graphics_info_t::residue_info_dialog = widget;
   } else {

      // create (then store) a new one.
      
      widget = create_residue_info_dialog();
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

GtkWidget *wrapped_nucleotide_builder_dialog() {

   GtkWidget *w = create_nucleotide_builder_dialog(); 
   return w;
} 

void ideal_nucleic_acid_by_widget(GtkWidget *builder_dialog) {

   std::string type = "RNA";
   std::string form = "A";
   short int single_stranded_flag = 0;
   GtkWidget *entry = lookup_widget(builder_dialog, "nucleotide_sequence");
   GtkWidget *type_optionmenu = lookup_widget(builder_dialog,
					      "nucleotide_builder_type_optionmenu");
   GtkWidget *form_optionmenu = lookup_widget(builder_dialog,
					      "nucleotide_builder_form_optionmenu");
   GtkWidget *strand_optionmenu = lookup_widget(builder_dialog,
						"nucleotide_builder_strand_optionmenu");


   GtkWidget *menu;
   GtkWidget *active_item;
   int active_index;

   menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(type_optionmenu));
   active_item = gtk_menu_get_active(GTK_MENU(menu));
   active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);
   std::cout << "DEBUG:: active_index for type: " << active_index << std::endl;
   if (active_index == 1)
      type = "DNA";

   menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(form_optionmenu));
   active_item = gtk_menu_get_active(GTK_MENU(menu));
   active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);
   std::cout << "DEBUG:: active_index for form: " << active_index << std::endl;
   if (active_index == 1)
      form = "B";

   menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(strand_optionmenu));
   active_item = gtk_menu_get_active(GTK_MENU(menu));
   active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);
   std::cout << "DEBUG:: active_index for strand: " << active_index << std::endl;
   if (active_index == 1)
      single_stranded_flag = 1;

   
   const char *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      ideal_nucleic_acid(type.c_str(), form.c_str(), single_stranded_flag, txt);
   }
}



GtkWidget *wrapped_create_display_control_window() { 

   graphics_info_t g;
   GtkWidget *widget = g.display_control_window();
   
   if (widget) { 
      // raise/uniconify (or whatever) what we have:
      // 
      if (!GTK_WIDGET_MAPPED(widget))
	 gtk_widget_show(widget);
      else
	 gdk_window_raise(widget->window);
   } else {

      // create (then store) a new one.
      
      widget = create_display_control_window_glade();
      save_display_control_widget_in_graphics(widget); 
      add_map_and_mol_display_control_widgets(); /* uses just saved widget */

      // set the size and shape of the window and paned components:
      if (graphics_info_t::display_manager_x_size != -1) {
	 gtk_window_set_default_size(GTK_WINDOW(widget),
			      graphics_info_t::display_manager_x_size,
				     graphics_info_t::display_manager_y_size);
	 if (graphics_info_t::display_manager_paned_position != -1) {
	    GtkPaned *paned = GTK_PANED(lookup_widget(widget,
						      "display_control_vpaned"));
	    gtk_paned_set_position(paned,
				   graphics_info_t::display_manager_paned_position);
	 }
      }
      if (graphics_info_t::display_manager_x_position != -1) {
	 gtk_widget_set_uposition(widget,
				  graphics_info_t::display_manager_x_position,
				  graphics_info_t::display_manager_y_position);
      }
   }
   return widget;
}

// BL things for file_chooser
#ifdef COOT_USE_GTK2_INTERFACE
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

#endif

// we add a universal function to set the file names
// in file chooser or selector

void set_filename_for_filechooserselection(GtkWidget *fileselection,
					   const gchar *filename) {

   bool chooser = 0;
#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
      chooser = 1;
      gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(fileselection),
					filename);     
   }
#endif

   if (chooser) {
	// dont do anything, as already done
   } else {
	   gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),  
					   filename);
   }
}

// functions to dock the accept/reject dialog
void set_accept_reject_dialog_docked(int istate){
   if (graphics_info_t::use_graphics_interface_flag) {
	  graphics_info_t::accept_reject_dialog_docked_flag = istate;
	  // reset the widget upon change of mode
	  set_accept_reject_dialog(0);
   }
}

int accept_reject_dialog_docked_state(){
  return graphics_info_t::accept_reject_dialog_docked_flag;
}

// functions to show/hide/sensitise docked accept/reject dialog
void set_accept_reject_dialog_docked_show(int state){
   if (graphics_info_t::use_graphics_interface_flag) {
	  graphics_info_t::accept_reject_dialog_docked_show_flag = state;
	  if (state == 0) {
		 GtkWidget *dialog = lookup_widget(GTK_WIDGET(graphics_info_t::glarea), "accept_reject_dialog_frame_docked");
		 // hide the widget and make sensitive again
		 gtk_widget_set_sensitive(dialog, TRUE);
		 gtk_widget_hide(dialog);
		 // reset the widget
		 set_accept_reject_dialog(0);
	  }
   }
}

int accept_reject_dialog_docked_show_state() {
  return graphics_info_t::accept_reject_dialog_docked_show_flag;
}

// functions for the refinement toolbar style
void set_model_toolbar_style(int istate) {
   graphics_info_t::model_toolbar_style_state = istate;
   if (graphics_info_t::use_graphics_interface_flag) { 
      GtkWidget *menuitem;
      if (istate <= 1) {
	 menuitem = lookup_widget(main_window(), "model_toolbar_icons1");
	 gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menuitem), TRUE);
      } else if (istate == 2) {
	 menuitem = lookup_widget(main_window(), "model_toolbar_icons_and_text1");
	 gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menuitem), TRUE);
      } else {
	 menuitem = lookup_widget(main_window(), "model_toolbar_text1");
	 gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menuitem), TRUE);
      }
   }
}

int model_toolbar_style_state() {
  return graphics_info_t::model_toolbar_style_state;
}

GtkWidget *wrapped_create_geometry_dialog() {
   graphics_info_t g;
   GtkWidget *w = NULL;
   if (g.geometry_dialog) { 
      w = g.geometry_dialog;
      // I'm not sure this magic does anything - it's a transient (and
      // I don't have a minimize handle for it).
      if (!GTK_WIDGET_MAPPED(w))
	 gtk_widget_show(w);
      else
	 gdk_window_raise(w->window);
   } else {
      w = create_geometry_dialog();
   }
   return w;
} 

void store_geometry_dialog(GtkWidget *w) { 

   graphics_info_t g;
   g.geometry_dialog = w;
}



void store_fixed_atom_dialog(GtkWidget *w) {

   graphics_info_t::fixed_atom_dialog = w;

} 

GtkWidget *wrapped_create_fixed_atom_dialog() {

   GtkWidget *w = create_fixed_atom_dialog();
   graphics_info_t::fixed_atom_dialog = w;
   return w;
}

int do_align_mutate_sequence(GtkWidget *w) {

   //
   int handled_state = 0;  // initially unhandled (return value).

   bool renumber_residues_flag = 0; // make this derived from the GUI one day
   int imol = graphics_info_t::align_and_mutate_imol;
   std::string chain_id = graphics_info_t::align_and_mutate_chain_from_optionmenu;
   GtkWidget *autofit_checkbutton = lookup_widget(w, "align_and_mutate_autofit_checkbutton");

   short int do_auto_fit = 0;
   if (GTK_TOGGLE_BUTTON(autofit_checkbutton)->active)
      do_auto_fit = 1;

   graphics_info_t g;
   int imol_refinement_map = g.Imol_Refinement_Map();

   short int early_stop = 0;
   if (do_auto_fit == 1)
      if (imol_refinement_map == -1)
	 early_stop = 1;

   if (early_stop) {
      std::string s = "WARNING:: autofit requested, but \n   refinement map not set!";
      std::cout << s << "\n";
      GtkWidget *warn = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(warn);
   } else { 

      handled_state = 1;
      if (imol >= 0) {
	 GtkWidget *text = lookup_widget(w, "align_and_mutate_sequence_text");
	 char *txt = NULL;
      
#if (GTK_MAJOR_VERSION == 1) 
	 gint start_pos = 0;
	 gint end_pos = -1;
	 txt = gtk_editable_get_chars(GTK_EDITABLE(text), start_pos, end_pos);
#else
	 // std::cout << "Gtk2 text view code... " << std::endl;
	 // text is a GtkTextView in GTK2
	 GtkTextView *tv = GTK_TEXT_VIEW(text);
	 GtkTextBuffer* tb = gtk_text_view_get_buffer(tv);
	 GtkTextIter startiter;
	 GtkTextIter enditer;
	 gtk_text_buffer_get_iter_at_offset(tb, &startiter, 0);
	 gtk_text_buffer_get_iter_at_offset(tb, &enditer, -1);
	 txt = gtk_text_buffer_get_text(tb, &startiter, &enditer, 0);
#endif 	 
      
	 if (txt) {
	    std::string sequence(txt);

	    if (is_valid_model_molecule(imol)) {
	       graphics_info_t g;
	       g.mutate_chain(imol, chain_id, sequence, do_auto_fit, renumber_residues_flag);
	       graphics_draw();
	    }
	 }
      } else {
	 std::cout << "WARNING:: inapproproate molecule number " << imol << std::endl;
      }
   }
   return handled_state;
}


void align_and_mutate_molecule_menu_item_activate(GtkWidget *item, 
						  GtkPositionType pos) {

   GtkWidget *chain_optionmenu = lookup_widget(item, "align_and_mutate_chain_optionmenu");
   GtkSignalFunc chain_callback = GTK_SIGNAL_FUNC(align_and_mutate_chain_option_menu_item_activate);
   graphics_info_t::align_and_mutate_imol = pos;
   int imol = pos;
   std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_optionmenu, imol,
								   chain_callback);
}

void align_and_mutate_chain_option_menu_item_activate (GtkWidget *item,
						       GtkPositionType pos) {

   graphics_info_t::align_and_mutate_chain_from_optionmenu = menu_item_label(item);
   std::cout << "align_and_mutate_chain_from_optionmenu is now "
	     << graphics_info_t::align_and_mutate_chain_from_optionmenu
	     << std::endl;
}


/* not for user consumption, this finds (from itself) the residue type
   and calls the graphics_info_t function. */
void fill_chi_angles_vbox(GtkWidget *vbox) {

   graphics_info_t g;
   gchar *strval = (gchar *) gtk_object_get_user_data(GTK_OBJECT(vbox));
   g.fill_chi_angles_vbox(vbox, strval);
}



GtkWidget *wrapped_create_add_additional_representation_gui() {

   GtkWidget *w = 0;
   if (graphics_info_t::use_graphics_interface_flag) {
      w = create_add_reps_dialog();
      // update/generate the option menu menu as usual.
      GtkWidget *option_menu = lookup_widget(w, "add_rep_molecule_optionmenu");
      GtkWidget *chain_id_entry = lookup_widget(w, "add_rep_chain_id_entry");
      GtkWidget *resno_start_entry = lookup_widget(w, "add_rep_resno_start_entry");
      GtkWidget *resno_end_entry = lookup_widget(w, "add_rep_resno_end_entry");
      GtkWidget *ins_code_entry = lookup_widget(w, "add_rep_ins_code_entry");
      GtkWidget *string_selection_entry = lookup_widget(w, "add_rep_selection_string_entry");
      
      GtkWidget *position_radiobutton = lookup_widget(w, "add_rep_radiobutton_position");
      GtkWidget *resno_radiobutton = lookup_widget(w, "add_rep_radiobutton_res_number");
      GtkWidget *selection_string_radiobutton = lookup_widget(w, "add_rep_radiobutton_selection_string");
      
      GtkWidget *add_reps_fat_bonds_radiobutton = lookup_widget(w, "add_rep_rep_fat_bonds_radiobutton");
      GtkWidget *add_reps_ball_and_stick_radiobutton = lookup_widget(w, "add_rep_rep_ball_and_stick_radiobutton");
      GtkWidget *add_rep_bond_width_combobox = lookup_widget(w, "add_rep_bond_width_combobox");

      GtkSignalFunc signal_func = 
	 GTK_SIGNAL_FUNC(add_reps_molecule_option_menu_item_select);
      int imol_active_position = graphics_info_t::add_reps_molecule_option_menu_item_select_molecule;
      fill_option_menu_with_coordinates_options(option_menu,  signal_func, imol_active_position);

#ifdef HAVE_GTK_COMBO_BOX_GET_ACTIVE_TEXT
      // set the active item to be the 8
      gtk_combo_box_set_active(GTK_COMBO_BOX(add_rep_bond_width_combobox), 7);
#endif      
   } 
   return w;
}


void add_reps_molecule_option_menu_item_select(GtkWidget *item, GtkPositionType pos) { 
   graphics_info_t::add_reps_molecule_option_menu_item_select_molecule = pos;
}


void add_additional_representation_by_widget(GtkWidget *w) {
   GtkWidget *option_menu = lookup_widget(w, "add_rep_molecule_optionmenu");
   
   GtkWidget *chain_id_entry = lookup_widget(w, "add_rep_chain_id_entry");
   GtkWidget *resno_start_entry = lookup_widget(w, "add_rep_resno_start_entry");
   GtkWidget *resno_end_entry = lookup_widget(w, "add_rep_resno_end_entry");
   GtkWidget *ins_code_entry = lookup_widget(w, "add_rep_ins_code_entry");
   GtkWidget *string_selection_entry = lookup_widget(w, "add_rep_selection_string_entry");

   GtkWidget *position_radiobutton = lookup_widget(w, "add_rep_radiobutton_position");
   GtkWidget *resno_radiobutton = lookup_widget(w, "add_rep_radiobutton_res_number");
   GtkWidget *selection_string_radiobutton = lookup_widget(w, "add_rep_radiobutton_selection_string");

   GtkWidget *add_reps_fat_bonds_radiobutton = lookup_widget(w, "add_rep_rep_fat_bonds_radiobutton");
   GtkWidget *add_reps_ball_and_stick_radiobutton = lookup_widget(w, "add_rep_rep_ball_and_stick_radiobutton");
   GtkWidget *add_rep_bond_width_combobox = lookup_widget(w, "add_rep_bond_width_combobox");

   float bond_width = 8;
   int bonds_box_type = coot::NORMAL_BONDS;
   short int representation_type = coot::SIMPLE_LINES;
   if (GTK_TOGGLE_BUTTON(add_reps_ball_and_stick_radiobutton)->active) {
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
      
   graphics_info_t g;
   GtkWidget *dcw = g.display_control_window();
   int imol = graphics_info_t::add_reps_molecule_option_menu_item_select_molecule;
   if (GTK_TOGGLE_BUTTON(position_radiobutton)->active) {
      // std::cout << "By position" << std::endl;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > aas = active_atom_spec();
      if (aas.first) {
	 int imol_active = aas.second.first;
	 coot::atom_selection_info_t asi(aas.second.second.chain,
					 aas.second.second.resno,
					 aas.second.second.resno,
					 aas.second.second.insertion_code);
	 gl_context_info_t glci(graphics_info_t::glarea, graphics_info_t::glarea_2);
	 g.molecules[imol_active].add_additional_representation(representation_type,
								bonds_box_type,
								bond_width,
								draw_H_flag,
								asi, dcw, glci);
      } 
   } 
   if (GTK_TOGGLE_BUTTON(resno_radiobutton)->active) {
      // std::cout << "By chainid/resno" << std::endl;
      std::string chain_id = gtk_entry_get_text(GTK_ENTRY(chain_id_entry));
      std::string resno_1s = gtk_entry_get_text(GTK_ENTRY(resno_start_entry));
      std::string resno_2s = gtk_entry_get_text(GTK_ENTRY(resno_end_entry));
      std::string ins_code = gtk_entry_get_text(GTK_ENTRY(ins_code_entry));
      if (is_valid_model_molecule(imol)) {
	 int resno_1 = atoi(resno_1s.c_str());
	 int resno_2 = atoi(resno_2s.c_str());
	 coot::atom_selection_info_t asi(chain_id, resno_1, resno_2, ins_code);
	 gl_context_info_t glci(graphics_info_t::glarea, graphics_info_t::glarea_2);
	 graphics_info_t::molecules[imol].add_additional_representation(representation_type,
									bonds_box_type,
									bond_width,
									draw_H_flag,
									asi, dcw, glci);
      } 
   } 
   if (GTK_TOGGLE_BUTTON(selection_string_radiobutton)->active) {
      // std::cout << "By selection string" << std::endl;
      std::string s = gtk_entry_get_text(GTK_ENTRY(string_selection_entry));
      coot::atom_selection_info_t asi(s);
      gl_context_info_t glci(graphics_info_t::glarea, graphics_info_t::glarea_2);
      graphics_info_t::molecules[imol].add_additional_representation(representation_type,
								     bonds_box_type,
								     bond_width,
								     draw_H_flag,
								     asi, dcw, glci);
   }
   graphics_draw();
} 


#if (GTK_MAJOR_VERSION > 1) 
GtkWidget *wrapped_create_residue_editor_select_monomer_type_dialog() {
   GtkWidget *w = create_residue_editor_select_monomer_type_dialog();
   GtkWidget *combo_box = lookup_widget(w, "residue_editor_select_monomer_type_combobox");
   graphics_info_t g;
   std::vector<std::string> v = g.Geom_p()->monomer_types();
   // remove the 2 items that are already there from the glade interface (I suppose).
   gtk_combo_box_remove_text(GTK_COMBO_BOX(combo_box), 0);
   gtk_combo_box_remove_text(GTK_COMBO_BOX(combo_box), 0);
   for (unsigned int i=0; i<v.size(); i++) {
      // std::string s = coot::util::int_to_string(i);
      // s += " ";
      // s += v[i];
      std::string s = v[i];
      gtk_combo_box_append_text (GTK_COMBO_BOX (combo_box), s.c_str());
      gtk_combo_box_set_active(GTK_COMBO_BOX(combo_box), i);
   }
   return w;
}
#endif


#if (GTK_MAJOR_VERSION > 1) 
void show_restraints_editor_by_index(int menu_item_index) {
   
   graphics_info_t g;
   std::vector<std::string> v = g.Geom_p()->monomer_types();
   for (unsigned int i=0; i<v.size(); i++) {
      if (i==menu_item_index)
	 show_restraints_editor(v[i].c_str());
   }
} 
#endif 


#if (GTK_MAJOR_VERSION > 1) 
void show_restraints_editor(const char *monomer_type) {

   if (graphics_info_t::use_graphics_interface_flag) {

      if (! monomer_type) {
	 std::cout << "ERROR:: null monomer_type - no restraints editor" << std::endl;
      } else { 
	 graphics_info_t g;
	 coot::protein_geometry *pg = g.Geom_p();

	 std::pair<short int, coot::dictionary_residue_restraints_t> p =
	    pg->get_monomer_restraints(monomer_type);
   
	 if (p.first) { 
	    coot::dictionary_residue_restraints_t restraints = p.second;
      
	    coot::restraints_editor r;
	    r.fill_dialog(restraints);
	    g.restraints_editors.push_back(r);
	 }
      }
   }
}

#endif


// ===================================================================
//                   sequence view
// ===================================================================


void nsv(int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   if (is_valid_model_molecule(imol)) {
      GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
      if (w) {

	 // it already exists... just raise it and map it.

	 // GtkWidget *canvas = g.sequence_view_is_displayed[imol];
	 GtkWidget *canvas = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
	 // so what is the window (which we shall call widget)?

	 
	 GtkWidget *widget = lookup_widget(canvas, "nsv_dialog");

	 if (widget) { 
	    if (!GTK_WIDGET_MAPPED(widget)) {
	       gtk_widget_show(widget);
	    } else {
	       gdk_window_raise(widget->window);
	    }
	 } else { 

	    widget = lookup_widget(canvas, "sequence_view_dialog");
	    if (widget) { 
	       if (!GTK_WIDGET_MAPPED(widget)) {
		  gtk_widget_show(widget);
	       } else {
		  gdk_window_raise(widget->window);
	       }
	    }
	 }

      } else {
	 graphics_info_t g;
	 std::string name = g.molecules[imol].name_for_display_manager();
	 exptl::nsv *seq_view =
	    new exptl::nsv(g.molecules[imol].atom_sel.mol, name, imol,
			   g.use_graphics_interface_flag);
	 // 
	 g.set_sequence_view_is_displayed(seq_view->Canvas(), imol);
      }
   }
#endif // defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
}


void sequence_view_old_style(int imol) {
   
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   graphics_info_t g;
   if (g.molecules[imol].has_model()) {
      graphics_info_t g;

      GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
      if (w) {

	 // it already exists... just raise it and map it.

	 // GtkWidget *canvas = g.sequence_view_is_displayed[imol];
	 GtkWidget *canvas = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
	 // so what is the window (which we shall call widget)?
	 GtkWidget *widget = lookup_widget(canvas, "sequence_view_dialog");

	 if (widget) { 
	    if (!GTK_WIDGET_MAPPED(widget)) {
	       gtk_widget_show(widget);
	    } else {
	       gdk_window_raise(widget->window);
	    }
	 }

      } else {

	 // create a new one

	 // Let's have a name that has the leading / and .pdb stripped
	 //
	 std::string short_name;

	 std::string::size_type islash = g.molecules[imol].name_.find_last_of("/");
	 std::string tstring;
	 if (islash == string::npos) { 
	    // no slash found
	    tstring = g.molecules[imol].name_;
	 } else { 
	    tstring = g.molecules[imol].name_.substr(islash + 1);
	 }
      
	 std::string::size_type ipdb = tstring.rfind(".pdb");

	 if (ipdb == string::npos) {
	    std::cout << "INFO .pdb not found in filename" << std::endl;
	    short_name = tstring;
	 } else {
	    short_name = tstring.substr(0, ipdb);
	 }
	 coot::sequence_view *seq_view =
	    new coot::sequence_view(g.molecules[imol].atom_sel.mol,
				    short_name, imol);
	 g.set_sequence_view_is_displayed(seq_view->Canvas(), imol);
      }
   }
#endif // HAVE_GTK_CANVAS   
}

void
add_on_sequence_view_choices() {

   graphics_info_t g;
   GtkWidget *menu = lookup_widget(GTK_WIDGET(g.glarea), "seq_view_menu");

   if (menu) { 
      gtk_container_foreach(GTK_CONTAINER(menu),
			    my_delete_ramachandran_mol_option,
			    (gpointer) menu);
      for(int i=0; i<g.n_molecules(); i++) {
	 if (g.molecules[i].has_model()) {
	    std::string name;
	    name = graphics_info_t::molecules[i].dotted_chopped_name();
	    update_sequence_view_menu_manual(i, name.c_str());
	 }
      }
   }
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
      GtkWidget *w = lookup_widget(g.glarea, "toolbar_multi_refine_stop_button");
      if (w) {
	 if (state) { 
	    gtk_widget_show(w);
	 } else { 
	    gtk_widget_hide(w);
	 } 
      }
   }
}

void set_visible_toolbar_multi_refine_continue_button(short int state) {

   graphics_info_t g;
   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *w = lookup_widget(g.glarea, "toolbar_multi_refine_continue_button");
      if (w) {
	 if (state) {
	    gtk_widget_show(w);
	 } else { 
	    gtk_widget_hide(w);
	 } 
      } 
      toolbar_multi_refine_button_set_sensitive("cancel", 0);
   }
} 

void set_visible_toolbar_multi_refine_cancel_button(short int state) {

   graphics_info_t g;
   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *w = lookup_widget(g.glarea, "toolbar_multi_refine_cancel_button");
      if (w) {
	 if (state) {
	    gtk_widget_show(w);
	 } else { 
	    gtk_widget_hide(w);
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
	 w = lookup_widget(g.glarea, "toolbar_multi_refine_cancel_button");
      if (bt == "continue")
	 w = lookup_widget(g.glarea, "toolbar_multi_refine_continue_button");
      if (bt == "stop")
	 w = lookup_widget(g.glarea, "toolbar_multi_refine_stop_button");
      
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

GtkWidget *wrapped_create_map_shapening_dialog() {

   float sharpening_limit = graphics_info_t::map_sharpening_scale_limit;
   GtkWidget *w = create_map_sharpening_dialog();
   GtkSignalFunc signal_func = GTK_SIGNAL_FUNC(map_sharpening_map_select);

   GtkWidget *option_menu = lookup_widget(w, "map_sharpening_optionmenu");

   int imol = fill_option_menu_with_map_options(option_menu, signal_func);
   graphics_info_t::imol_map_sharpening = imol;

   std::cout << "DEBUG:: imol from fill_option_menu_with_map_options() "
	     << imol << std::endl;

   GtkWidget *h_scale = lookup_widget(w, "map_sharpening_hscale");
   //GtkObject *adj = gtk_adjustment_new(0.0, -sharpening_limit, 2*sharpening_limit,
   // 0.05, 2, 30.1);
   GtkObject *adj = gtk_adjustment_new(0.0, -sharpening_limit, 2*sharpening_limit,
				       0.05, 0.2, (sharpening_limit+0.1));
   gtk_range_set_adjustment(GTK_RANGE(h_scale), GTK_ADJUSTMENT(adj));

   gtk_signal_connect(GTK_OBJECT(adj), "value_changed",
		      GTK_SIGNAL_FUNC(map_sharpening_value_changed), NULL);

#if (GTK_MAJOR_VERSION > 1)
#if (GTK_MAJOR_VERSION > 2) || (GTK_MINOR_VERSION > 14)
   int ticks = 3;  // number of ticks on the (one) side (not including centre tick)
   for (int i=0; i<=2*ticks; i++) {
      float p = float (i-ticks) * (1.0/float(ticks)) * sharpening_limit;
      std::string pos_string = coot::util::float_to_string_using_dec_pl(p,1);
      gtk_scale_add_mark(GTK_SCALE(h_scale),
			 p,
			 GTK_POS_BOTTOM, pos_string.c_str());
   }
   gtk_scale_add_mark(GTK_SCALE(h_scale), -sharpening_limit, GTK_POS_BOTTOM, "\nSharpen");
   gtk_scale_add_mark(GTK_SCALE(h_scale),  sharpening_limit, GTK_POS_BOTTOM, "\nBlur");
#endif   
#endif

   // Don't display the cancel button.
   GtkWidget *c = lookup_widget(w, "map_sharpening_cancel_button");
   gtk_widget_hide(c);

   return w;
}

void
map_sharpening_map_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t::imol_map_sharpening = pos;
}

void map_sharpening_value_changed (GtkAdjustment *adj, 
				   GtkWidget *window) {

   int imol = graphics_info_t::imol_map_sharpening;
   if (is_valid_map_molecule(imol)) {
      // std::cout << "sharpen " << imol << " by " << adj->value << std::endl;
      sharpen(imol, adj->value);
   } 
}

void set_baton_build_params_from_widget(GtkWidget *params_dialog) {

   GtkWidget *ent_res_no   = lookup_widget(params_dialog, "baton_build_params_residue_number_entry");
   GtkWidget *ent_chain_id = lookup_widget(params_dialog, "baton_build_params_chain_id_entry");
   GtkWidget *check_button = lookup_widget(params_dialog, "baton_build_params_backwards_checkbutton");

   const char *resno_txt = gtk_entry_get_text(GTK_ENTRY(ent_res_no));
   const char *chain_id  = gtk_entry_get_text(GTK_ENTRY(ent_chain_id));

   const char *direction = "forwards";
   if (GTK_TOGGLE_BUTTON(check_button)->active)
      direction = "backwards";

   try {
      int rn = coot::util::string_to_int(resno_txt);
      set_baton_build_params(rn, chain_id, direction);
   }
   catch (std::runtime_error rte) {
      std::cout << rte.what() << " aborting." << std::endl;
   }

} 

