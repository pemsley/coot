/* src/c-interface-gui.cc
 * 
 * Copyright 2007 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
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
#define sleep Sleep
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER


#include "graphics-info.h"
#include "interface.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "cmtz-interface.hh"
#include "mmdb.h"  // for centre of molecule

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


void handle_column_label_make_fourier(GtkWidget *column_label_window) {

  GtkWidget *refmac_checkbutton;
   int icol; 
   int use_weights = 0;
   int is_diff_map;
   short int sensible_r_free_col;
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
     
   GtkWidget *menuitem;
   GtkWidget *active_menu;
   int imol;
   int *n;
   
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
     if (icol < saved_f_phi_columns->f_cols.size()) 
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


/* Return 1 if we moved to a molecule centre, else go to origin and
   return 0. */
/* centre on last-read (and displayed) molecule with zoom 100. */
// 
// However, if we are already *at* that molecule centre, Reset View
// moves to the centre of the next displayed molecule (with wrapping).
// 
int reset_view() {

   int istat = 0;

   // Question:  Are we are a molecule centre?
   // If we are, move to the next available centre.
   // If not, then move to the last read molecule.
   //
   graphics_info_t g;
   std::vector<coot::Cartesian> molecule_centres(graphics_info_t::n_molecules,
						 coot::Cartesian(0,0,0));
   int centred_on_molecule_number = -1;
   coot::Cartesian current_centre = g.RotationCentre();
   coot::Cartesian new_centre(0,0,0); // gets set.
   int last_molecule = -1;
   
   for (int imol=(graphics_info_t::n_molecules -1); imol>=0; imol--) {
      if (graphics_info_t::molecules[imol].is_displayed_p()) {
	 if (last_molecule == -1)
	    last_molecule = imol;
	 coot::Cartesian mc = centre_of_molecule(graphics_info_t::molecules[imol].atom_sel);
	 molecule_centres[imol] = mc;
	 if ((mc - current_centre).length() < 0.1) {
	    if (centred_on_molecule_number == -1) {
	       centred_on_molecule_number = imol;
	    }
	 }
      }
   }

   // If we were not centred on a molecule then centre on the last
   // available molecule.
   //
   if (centred_on_molecule_number == -1) {
      if (last_molecule != -1) { 
	 new_centre = molecule_centres[last_molecule];
	 std::string s = "Centring no molecule number ";
	 s += g.int_to_string(last_molecule);
	 s += " ";
	 s += graphics_info_t::molecules[last_molecule].name_for_display_manager();
	 g.statusbar_text(s);
	 istat = 1;
      } else {
	 std::string s = "No displayed molecules";
	 g.statusbar_text(s);
	 new_centre = current_centre;
      }
   } else {

      // OK, we were centred on a molecule... which is the next one we
      // want to centre on?
      // Let's make a list of the available molecules:
      std::vector<int> available_molecules;
      for(int imol=0; imol<graphics_info_t::n_molecules; imol++) {
	 if (graphics_info_t::molecules[imol].is_displayed_p()) {
	    available_molecules.push_back(imol);
	 }
      }
      
      if (available_molecules.size() == 1) {
	 // no other molecule to centre on.
	 new_centre = current_centre;
      } else {

	 // we want the molecule after the molecule that we are
	 // currently centred on, if not that (which may be because we
	 // are at the last molecule in the list) then the first
	 // molecule in the list.
	 // 
	 int first_in_list = available_molecules[0];
	 int next = -1;
	 for (unsigned int iav=0; iav<available_molecules.size(); iav++) {
	    if (available_molecules[iav] > centred_on_molecule_number) {
	       next = available_molecules[iav];
	       break;
	    }
	 }
	 if (next == -1)
	    next = first_in_list;
	 new_centre = molecule_centres[next];
	 std::string s = "Centring on molecule number ";
	 s += g.int_to_string(next);
	 s += " ";
	 s += graphics_info_t::molecules[next].name_for_display_manager();
	 g.statusbar_text(s);
	 istat = 1;
      }
   }

   g.setRotationCentreAndZoom(new_centre, 100.0);
   g.zoom = 100.0;
      
   for(int ii=0; ii<graphics_info_t::n_molecules; ii++) {
      graphics_info_t::molecules[ii].update_map();
      graphics_info_t::molecules[ii].update_symmetry();
   }
   graphics_draw();
   add_to_history_simple("reset-view");
   return istat;
   
}

void fill_about_window(GtkWidget *widget) {

   GtkWidget *text_widget;
   text_widget = lookup_widget(GTK_WIDGET(widget), "about_window_text");

   std::string body_text("\n\n   Brought to you by:\n\n   Paul Emsley & Kevin Cowtan\n\n   Using the dictionaries of:\n    Alexei Vagin\n    Roland Dunbrack & co-workers\n\n  Using the libraries of:\n   Eugene Krissinel\n   Kevin Cowtan\n   Stuart McNicholas\n   Ralf W. Grosse-Kunstleve\n   Janne Lof\n   Raghavendra Chandrashekara\n   Paul Bourke & Cory Gene Bloyd\n   Matteo Frigo & Steven G. Johnson\n   & many others.\n\n  Windows 2000 Binaries\n   Bernhard Lohkamp\n\n  Macintosh Binaries\n   William Scott\n\n"); 

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

void set_graphics_window_size(int x_size, int y_size) {

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      g.graphics_x_size = x_size;
      g.graphics_y_size = y_size;
      if (g.glarea) { 
	 GtkWindow *window = GTK_WINDOW(lookup_widget(g.glarea, "window1"));
	 if (window) { 
	    gtk_window_set_default_size (window, x_size, y_size);
	    while (gtk_events_pending())
	       gtk_main_iteration();
	 }
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("set-graphics-window-size");
   command_strings.push_back(graphics_info_t::int_to_string(x_size));
   command_strings.push_back(graphics_info_t::int_to_string(y_size));
   add_to_history(command_strings);
}

void store_graphics_window_position(int x_pos, int y_pos) { 

   graphics_info_t g;
   g.graphics_x_position = x_pos;
   g.graphics_y_position = y_pos;

   std::string cmd = "store-graphics-window-position";
   std::vector<coot::command_arg_t> args;
   args.push_back(x_pos);
   args.push_back(y_pos);
   add_to_history_typed(cmd, args);
} 

void set_graphics_window_position(int x_pos, int y_pos) {

   if (graphics_info_t::use_graphics_interface_flag) { 
      graphics_info_t g;
      GtkWidget *main = lookup_widget(g.glarea, "window1");
      if (main) { 
	 gtk_widget_set_uposition(main, x_pos, y_pos);
	 while (gtk_events_pending())
	    gtk_main_iteration();
	 while (gdk_events_pending())
	    gtk_main_iteration();
      }
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
   gdk_window_get_root_origin (widget->window, &upositionx, &upositiony);

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
      gtk_widget_set_usize(dialog,
			   graphics_info_t::file_selection_dialog_x_size,
			   graphics_info_t::file_selection_dialog_y_size);
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
   int i = g.check_for_unsaved_changes();
   std::string cmd = "coot-checked-exit";
   std::vector<coot::command_arg_t> args;
   args.push_back(retval);
   add_to_history_typed(cmd, args);
   if (i == 0)
      coot_real_exit(retval);
   return TRUE; // path where there were unsaved changes, we don't
		// want to exit.

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
   gtk_exit(retval); 

}

void add_recentre_on_read_pdb_checkbutton(GtkWidget *fileselection) { 

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

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
   w = create_coords_fileselection1 ();
#else
   w = create_coords_fileselection1();
   // w = 0;
#endif
   return w;
}

void set_directory_for_coot_file_chooser(GtkWidget *coords_fileselection1) {

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
      set_directory_for_fileselection(coords_fileselection1);
#else

#endif
}

const char *coot_file_chooser_file_name(GtkWidget *widget) {

   const char *f = 0;
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
#else

#endif
   return f;
} 



/*  ------------------------------------------------------------------------ */
//            reparenting
/*  ------------------------------------------------------------------------ */
// 
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
   
   button = lookup_widget(widget,
			  "model_refine_dialog_rot_trans_togglebutton");
   if (button) {
      if (graphics_info_t::model_fit_refine_rotate_translate_zone_string != "")
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),
			    graphics_info_t::model_fit_refine_rotate_translate_zone_string.c_str());
      
   }
   graphics_info_t::set_model_fit_refine_button_names(widget);

   return widget; 
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
   coot::mtz_column_types_info_t r = coot::get_f_phi_columns(mtz_file_name);
   if (r.f_cols.size() > 0)
      return 1;
   else
      return 0;
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
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) { 
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
	       if (graphics_info_t::sequence_view_is_displayed[imol]) {
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
		  GtkWidget *canvas = graphics_info_t::sequence_view_is_displayed[imol];
		  GtkWidget *window = lookup_widget(canvas, "sequence_view_dialog");
		  gtk_widget_destroy(window);
	       } 
#endif
	       graphics_info_t::molecules[imol].close_yourself();
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
	    for (int imol=graphics_info_t::n_molecules-1; imol>=0; imol--) {
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
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) { 
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
      int *imol = (int *) gtk_object_get_user_data(GTK_OBJECT(active_item));

      if (imol) { 

	 std::cout << " Closing molecule number " << *imol << std::endl;
	 graphics_info_t g;
	 if ((*imol >= 0) && (*imol < g.n_molecules)) { 
	    g.molecules[*imol].close_yourself();
	    fill_close_option_menu_with_all_molecule_options(optionmenu);
	    graphics_draw();
	 
	 } else { 
	    std::cout << "ERROR: Closing invalid molecule number" << *imol 
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

   if (is_valid_model_molecule(imol) ||
       is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].close_yourself();
   }
   if (graphics_info_t::go_to_atom_window) {
      std::cout << ".....re fill go to atom window here" << std::endl;
      graphics_info_t g;
      if (imol == g.go_to_atom_molecule()) {
	 g.update_go_to_atom_window_on_new_mol();
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

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      
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

   if (graphics_n_molecules() > 0) {
      gtk_menu_set_active(GTK_MENU(menu), 0);
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu),
			    menu);
   
}



void
close_molecule_item_select(GtkWidget *item, GtkPositionType pos) {

   std::cout << "activating closing position/imol " << pos << std::endl;
   
}

void add_ccp4i_project_optionmenu(GtkWidget *fileselection) {

   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;

   GtkWidget *optionmenu = gtk_option_menu_new();
   gtk_widget_ref(optionmenu);
   gtk_widget_show(optionmenu);
   GtkSignalFunc project_signal_func =
      GTK_SIGNAL_FUNC(option_menu_refmac_ccp4i_project_signal_func);
   add_ccp4i_projects_to_optionmenu(optionmenu, project_signal_func);

   // lets put the optionmenu in a frame with a label
   GtkWidget *frame = gtk_frame_new("CCP4i Project Directory");
   gtk_container_add(GTK_CONTAINER(aa), frame);
   gtk_widget_show(frame);
   gtk_container_add(GTK_CONTAINER(frame),optionmenu);
}

void add_ccp4i_projects_to_optionmenu(GtkWidget *optionmenu,
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
      gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
			 GTK_SIGNAL_FUNC(func),
			 GINT_TO_POINTER(i));
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu), menu);

   // set the active menu item here...
   
   gtk_widget_show(menu);
}

void
option_menu_refmac_ccp4i_project_signal_func(GtkWidget *item, GtkPositionType pos) {
   graphics_info_t g;
   g.ccp4_projects_index_last = pos;
   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > pr_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);
   if (pos < int(pr_pairs.size())) {
      g.set_directory_for_fileselection_string(pr_pairs[pos].second);

      // we need to call set_directory_for_fileselection with an
      // argument that is the fileselection widget.  The question is:
      // what is that widget?
      //
      // Here, we don't know, so we try to look up each of the
      // fileselection widgets
      //
      GtkWidget *fileselection;
      fileselection = lookup_file_selection_widgets(item);

      if (fileselection) {
	 g.set_directory_for_fileselection(fileselection);
      } else {
	 std::cout << "WARNING:: failed to find filesection in "
		   << "option_menu_refmac_ccp4i_project_signal_func" << std::endl;
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


GtkWidget *lookup_file_selection_widgets(GtkWidget *item) {

   GtkWidget *w;
   w = lookup_widget(GTK_WIDGET(item), "coords_fileselection1");
   if (! w) {
      w = lookup_widget(GTK_WIDGET(item), "dataset_fileselection1");
      if (! w) {
	 w = lookup_widget(GTK_WIDGET(item), "map_name_fileselection1");
	 if (! w) {
	    w = lookup_widget(GTK_WIDGET(item), "phs_coordinates_fileselection");
	    if (! w) {
	       w = lookup_widget(GTK_WIDGET(item), "save_coords_fileselection1");
	       if (! w) {
		  w = lookup_widget(GTK_WIDGET(item), "cif_dictionary_fileselection");
		  if (! w) {
		     w = lookup_widget(GTK_WIDGET(item), "run_script_fileselection");
		  }
	       }
	    }
	 }
      }
   }
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
   int imol = 0; 

   
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
 
    if (get_show_unit_cell(imol) == 1) { 
       button = GTK_BUTTON(lookup_widget(show_symm_window,
					 "unit_cell_yes_radiobutton"));
    } else { 
       button = GTK_BUTTON(lookup_widget(show_symm_window,
					 "unit_cell_no_radiobutton"));
    }

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
       
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

/*    printf("Clipping adjustment: %f\n", adj->value); */

/* We will have a mol_no at some stage, set it to 1 for now. */

   set_symmetry_colour_merge(0, adj->value); /* this adjusts
						graphics_info_t::
						symm_colour_merge_weight,
						which is a double
						array.  But we only
						ever use the 0th
						position of it in
						combine_colour() */

}

void
handle_map_colour_change(int mol, gdouble* map_col) {
   
   graphics_info_t::molecules[mol].handle_map_colour_change(map_col,
							    graphics_info_t::swap_difference_map_colours);

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
      graphics_info_t::molecules[imol].handle_map_colour_change(colours, graphics_info_t::swap_difference_map_colours);
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
      for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
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
      for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
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

//    std::cout << "DEBUG:: in wrapped_create_bond_parameters_dialog imol is "
// 	     << imol << std::endl;

   g.fill_option_menu_with_coordinates_options(optionmenu, callback_func, imol);
   graphics_info_t::fill_bond_parameters_internals(widget, imol);

   return widget;
}

void apply_bond_parameters(GtkWidget *w) {

   graphics_info_t g;
   int imol = g.bond_parameters_molecule;

   
   if (imol >= 0) {
      if (imol < g.n_molecules) {
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

	    // Draw NCS ghosts?

	    GtkWidget *ncs_toggle_button =
	       lookup_widget(w, "draw_ncs_ghosts_yes_radiobutton");
	    if (GTK_TOGGLE_BUTTON(ncs_toggle_button)->active) {
	       std::cout << "set_draw_ncs_ghosts " << imol << " " << "1" << std::endl;
	       set_draw_ncs_ghosts(imol, 1);
	    } else {
	       std::cout << "set_draw_ncs_ghosts " << imol << " " << "0" << std::endl;
	       set_draw_ncs_ghosts(imol, 0);
	    }
	    
	    // bye bye colour map rotation entry

// 	    GtkWidget *entry = lookup_widget(w, "bond_parameters_colour_map_rotation_entry");
// 	    std::pair<short int, float> p = float_from_entry(entry);
// 	    // 	    std::cout << "DEBUG:: colour map rotation: float from entry pair: "
// 	    //	     << p.first << " " << p.second
// 	    // << std::endl;
// 	    if (p.first) {
// 	       set_colour_map_rotation_on_read_pdb(p.second);
// 	    } else {
// 	       set_colour_map_rotation_on_read_pdb(32.0);
// 	    }

	    // bye bye colour map rotation button

	    // colour map rotate carbon only?

// 	    GtkWidget *checkbutton =
// 	       lookup_widget(w, "bond_parameters_rotate_colour_map_c_only_checkbutton");
// 	    if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
// 	       set_colour_map_rotation_on_read_pdb_c_only_flag(1);
// 	    } else {
// 	       set_colour_map_rotation_on_read_pdb_c_only_flag(0);
// 	    } 
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
   g.set_file_for_save_fileselection(fileselection);
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
   
      const gchar *filename = gtk_file_selection_get_filename
	 (GTK_FILE_SELECTION(widget));

      std::cout << "save coordinates for molecule "
		<< imol << " to file " << filename << std::endl;

      graphics_info_t g;
      int ierr = g.molecules[imol].save_coordinates(filename);
      if (! ierr) { 
	 std::string s = "Saved coordinates file ";
	 s += filename;
	 s += ".";
	 g.statusbar_text(s);
      }
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
      if (graphics_info_t::model_fit_refine_x_position > -1) { 
	 gtk_widget_set_uposition(widget,
				  graphics_info_t::go_to_atom_window_x_position,
				  graphics_info_t::go_to_atom_window_y_position);
      }
      fill_go_to_atom_window(widget);
   } 
   return widget;
} 

gchar *get_text_for_go_to_atom_chain_entry() { 

   //
   graphics_info_t g; 

   gchar *txt = (gchar *)malloc(100);
   strcpy(txt, g.go_to_atom_chain()); 

   return txt; 
   
}

gchar *get_text_for_go_to_atom_residue_entry() { 

   // 
   graphics_info_t g;

   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%d", g.go_to_atom_residue()); 
   
   return txt; 

}
 
gchar *get_text_for_go_to_atom_atom_name_entry() { 

   // 
   graphics_info_t g;

   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_atom_name()); 
   
   return txt; 

}

void post_go_to_atom_window() { 
   GtkWidget *widget = wrapped_create_goto_atom_window();
   gtk_widget_show(widget);
   std::vector<std::string> command_strings;
   command_strings.push_back("post-go-to-atom-window");
   add_to_history(command_strings);
}

void fill_go_to_atom_window(GtkWidget *widget) {

     GtkWidget *option_menu; 
     GtkWidget *chain_entry;
     GtkWidget *residue_entry; 
     GtkWidget *atom_name_entry; 
     GtkWidget *atom_gtklist;
     GtkWidget *residue_gtklist;
     GtkWidget *residue_tree;
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


     text = get_text_for_go_to_atom_residue_entry(); // on startup, tinkers with
                                                     // go to atom params, yuck,
                                                     // I think.

     gtk_entry_set_text(GTK_ENTRY(residue_entry), text); 

     /* Now that the go to atom molecule has been set, we can use it
	to fill the molecule option menu */
     int gimol = g.go_to_atom_molecule();
     g.fill_option_menu_with_coordinates_options(option_menu,
						 callback_func,
						 gimol);

     /* The chain entry */

     chain_entry = lookup_widget(GTK_WIDGET(widget),
				 "go_to_atom_chain_entry");

     text = get_text_for_go_to_atom_chain_entry(); 
     gtk_entry_set_text(GTK_ENTRY(chain_entry), text); 

     
     /* The Atom Name entry */

     atom_name_entry = lookup_widget(GTK_WIDGET(widget),
				     "go_to_atom_atom_name_entry"); 
     text = get_text_for_go_to_atom_atom_name_entry(); 
     gtk_entry_set_text(GTK_ENTRY(atom_name_entry), text); 

     /* The Residue List */

     /* The residue list cant be added to a scrolled window in glade,
	so we create only a scrolled window in glade
	(go_to_atom_residue_scrolledwindow) and add the list to it
	like is done in examples/list/list.c */

     scrolled_window = lookup_widget(GTK_WIDGET(widget),
				     "go_to_atom_residue_scrolledwindow");
     residue_gtklist=gtk_list_new();

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

     residue_tree=gtk_tree_new();
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
     atom_gtklist=gtk_list_new();
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

   std::cout << "applying go to atom values!" << std::endl;

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

  std::cout << "set_go_to_atom_chain_residue_atom_name_strings!" << std::endl;
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
