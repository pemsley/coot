/* src/c-interface-build-gui.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007,2008, 2009 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009 The University of Oxford
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

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <stdlib.h>
#include <iostream>

#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#include <windows.h>
#endif
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.

#include <vector>
#include <string>

#include <mmdb/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h"

#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "graphics-info.h"

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-fasta.hh"

#include "skeleton/BuildCas.h"
#include "ligand/helix-placement.hh"
#include "ligand/fast-ss-search.hh"

#include "trackball.h" // adding exportable rotate interface

#include "utils/coot-utils.hh"  // for is_member_p
#include "coot-utils/coot-map-heavy.hh"  // for fffear

#include "guile-fixups.h"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: but in order that we do not get error: "_POSIX_C_SOURCE"
// redefined problems, we should include python at the
// beginning. Double grr!
//
// #ifdef USE_PYTHON
// #include "Python.h"
// #endif // USE_PYTHON


#include "c-interface.h"
#include "cc-interface.hh"

#include "ligand/ligand.hh" // for rigid body fit by atom selection.

#include "cmtz-interface.hh" // for valid columns mtz_column_types_info_t
#include "c-interface-mmdb.hh"
#include "c-interface-scm.hh"
#include "c-interface-python.hh"

#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else 
#include "ligand/richardson-rotamer.hh"
#endif 

void do_regularize_kill_delete_dialog() {
   graphics_info_t g;
   if (g.delete_item_widget) { 
      gtk_widget_destroy(g.delete_item_widget);
      g.delete_item_widget = NULL;
      // hopefully superfluous:
      g.delete_item_atom = 0;
      g.delete_item_residue = 0;
      g.delete_item_residue_hydrogens = 0;
   }
}
   

/* moving gtk functionn out of build functions, delete_atom() updates
   the go to atom atom list on deleting an atom  */
void update_go_to_atom_residue_list(int imol) {

   graphics_info_t g;
      if (g.go_to_atom_window) {
	 int go_to_atom_imol = g.go_to_atom_molecule();
	 if (go_to_atom_imol == imol) { 

	    // The go to atom molecule matched this molecule, so we
	    // need to regenerate the residue and atom lists.
	    GtkWidget *gtktree = lookup_widget(g.go_to_atom_window,
					       "go_to_atom_residue_tree");
	    GtkWidget *gtk_atom_list = lookup_widget(g.go_to_atom_window,
						     "go_to_atom_atom_list");
#if (GTK_MAJOR_VERSION == 1)
	    g.fill_go_to_atom_residue_list_gtk1(gtktree);
#else 	    
	    g.fill_go_to_atom_residue_tree_and_atom_list_gtk2(imol, gtktree, gtk_atom_list);
#endif	    
	 } 
      }
}

/* utility function, moving widget work out of c-interface-build.cc */
void delete_object_handle_delete_dialog(short int do_delete_dialog) {
   if (graphics_info_t::delete_item_widget != NULL) {
      if (do_delete_dialog) { // via ctrl

	 // another check is needed, is the check button active?
	 // 
	 // If not we can go ahead and delete the dialog
	 //
	 GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
						"delete_item_keep_active_checkbutton");
	 if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	    // don't kill the widget
	    pick_cursor_maybe(); // it was set to normal_cursor() in
                                 // graphics-info-define's delete_item().
	 } else {
	 
	    gint upositionx, upositiony;
	    gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
					&upositionx, &upositiony);
	    graphics_info_t::delete_item_widget_x_position = upositionx;
	    graphics_info_t::delete_item_widget_y_position = upositiony;
	    gtk_widget_destroy(graphics_info_t::delete_item_widget);
	    graphics_info_t::delete_item_widget = NULL;
	    graphics_draw();
	 }
      }
   }
} 



/* We need to set the pending delete flag and that can't be done in
   callback, so this wrapper does it */
GtkWidget *wrapped_create_delete_item_dialog() {

   GtkWidget *widget = create_delete_item_dialog();
   GtkWidget *atom_toggle_button;

   if (delete_item_mode_is_atom_p()) { 
      atom_toggle_button = lookup_widget(GTK_WIDGET(widget),
					 "delete_item_atom_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(atom_toggle_button), TRUE);
      std::cout << "Click on the atom that you wish to delete\n";
   } else {
      if (delete_item_mode_is_water_p()) {
	 GtkWidget *water_toggle_button = lookup_widget(widget,
							"delete_item_water_radiobutton");
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(water_toggle_button), TRUE);
      } else { 
	 set_delete_residue_mode(); // The widget default radio button
	 std::cout << "Click on an atom in the residue that you wish to delete\n";
      }
   }
   graphics_info_t::pick_pending_flag = 1;
   pick_cursor_maybe();
   set_transient_and_position(COOT_DELETE_WINDOW, widget);
   store_delete_item_widget(widget);
   return widget; 
}

// -----------------------------------------------------
//  move molecule here widget
// -----------------------------------------------------
GtkWidget *wrapped_create_move_molecule_here_dialog() {

   GtkWidget *w = create_move_molecule_here_dialog();
   fill_move_molecule_here_dialog(w);
   return w;
}

// called (also) by the callback of toggling the
// move_molecule_here_big_molecules_checkbutton.
// 
void
fill_move_molecule_here_dialog(GtkWidget *w) {

   GtkWidget *option_menu  = lookup_widget(w, "move_molecule_here_optionmenu");
   GtkWidget *check_button = lookup_widget(w, "move_molecule_here_big_molecules_checkbutton");

   bool fill_with_small_molecule_only_flag = 1;
   int imol_active = first_small_coords_imol();
   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(graphics_info_t::move_molecule_here_item_select);

   if (check_button) {
      if (GTK_TOGGLE_BUTTON(check_button)->active) { 
	 fill_with_small_molecule_only_flag = 0;
	 imol_active = first_coords_imol();
      }
   }

   graphics_info_t g;
   graphics_info_t::move_molecule_here_molecule_number = imol_active;
   g.fill_option_menu_with_coordinates_options_possibly_small(option_menu, callback_func, imol_active,
							      fill_with_small_molecule_only_flag);
} 


void move_molecule_here_by_widget(GtkWidget *w) {

   int imol = graphics_info_t::move_molecule_here_molecule_number;
   move_molecule_to_screen_centre_internal(imol);
   std::vector<std::string> command_strings;
   command_strings.push_back("move-molecule-here");
   command_strings.push_back(clipper::String(imol));
   add_to_history(command_strings);
}

// This is a copy - more or less - of
// fill_option_menu_with_coordinates_options, except we also add at
// the top "New Molecule" if a molecule by the name of "Pointer Atoms"
// is not found.
// 
// Note that we can't use fill_option_menu_with_coordinates_options
// and add to it because gtk_menu_set_active fails/is ignored.
// 
// fill_pointer_atom_molecule_option_menu
// 
void fill_place_atom_molecule_option_menu(GtkWidget *optionmenu) { 

   GtkSignalFunc callback_func = 
      GTK_SIGNAL_FUNC(graphics_info_t::pointer_atom_molecule_menu_item_activate);
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
   if (menu) 
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   GtkWidget *menuitem;
   int pointer_atoms_mol = -1;

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) { 
	 if (graphics_info_t::molecules[imol].name_ == "Pointer Atoms") { 
	    pointer_atoms_mol = imol;
	 } 
      }
   }

   int menu_index = 0;

   if (pointer_atoms_mol == -1) { 
      // There were no pointer atoms so let's create "New Molecule" at
      // the top of the list.
      // 
      GtkWidget *menu_item = gtk_menu_item_new_with_label("New Molecule");
      int imol_new = -10;
      gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
			 callback_func,
			 GINT_TO_POINTER(imol_new));
      gtk_menu_append(GTK_MENU(menu), menu_item); 
      gtk_widget_show(menu_item); 
      gtk_menu_set_active(GTK_MENU(menu), 0);
      menu_index = 0;
   }

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 std::string ss = graphics_info_t::int_to_string(imol);
	 ss += " " ;
	 int ilen = graphics_info_t::molecules[imol].name_.length();
	 int left_size = ilen-graphics_info_t::go_to_atom_menu_label_n_chars_max;
	 if (left_size <= 0) {
	    // no chop
	    left_size = 0;
	 } else {
	    // chop
	    ss += "...";
	 }
	 ss += graphics_info_t::molecules[imol].name_.substr(left_size, ilen);
	 menuitem =  gtk_menu_item_new_with_label (ss.c_str());
	 menu_index++;
	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     callback_func,
			     GINT_TO_POINTER(imol));
	 gtk_menu_append(GTK_MENU(menu), menuitem); 
	 gtk_widget_show(menuitem);

	 // set any previously saved active position:
	 if (graphics_info_t::user_pointer_atom_molecule == imol) {
	    std::cout << "setting active menu item to "
		      << menu_index << std::endl;
	    gtk_menu_set_active(GTK_MENU(menu), menu_index);
	 }
      }
   }
   
   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu),
			    menu);
}

/* Now the refinement weight can be set from an entry in the refine_params_dialog. */
void set_refinemenent_weight_from_entry(GtkWidget *entry) {

   const char *text = gtk_entry_get_text(GTK_ENTRY(entry));
   try {
      float f = coot::util::string_to_float(text);
      graphics_info_t::geometry_vs_map_weight = f;
   }
   catch (std::runtime_error rte) {
      std::cout << "in set_refinemenent_weight_from_entry " << rte.what() << std::endl;
   } 
} 

void place_atom_at_pointer_by_window() { 

   // put up a widget which has a OK callback button which does a 
   // g.place_typed_atom_at_pointer();
   GtkWidget *window = create_pointer_atom_type_dialog();
   
   GtkWidget *optionmenu = lookup_widget(window, "pointer_atom_molecule_optionmenu");
   //       GtkSignalFunc callback_func =
   // 	 GTK_SIGNAL_FUNC(graphics_info_t::pointer_atom_molecule_menu_item_activate);
   
   fill_place_atom_molecule_option_menu(optionmenu);
   gtk_widget_show(window);

}


// User data has been placed in the window - we use it to get the
// molecule number.
void baton_mode_calculate_skeleton(GtkWidget *window) {

   int imol = -1;

   int *i;

   std::cout << "getting intermediate data in baton_mode_calculate_skeleton "
	     << std::endl;
   i = (int *) gtk_object_get_user_data(GTK_OBJECT(window));

   std::cout << "got intermediate int: " << i << " " << *i << std::endl;

   imol = *i;

   std::cout << "calculating map for molecule " << imol << std::endl;
   if (imol < graphics_info_t::n_molecules() && imol >= 0) { 
      skeletonize_map(0, imol);
   }
}


GtkWidget *wrapped_create_renumber_residue_range_dialog() {

   GtkWidget *w = create_renumber_residue_range_dialog();
   int imol = first_coords_imol();
   graphics_info_t::renumber_residue_range_molecule = imol;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.fill_renumber_residue_range_dialog(w);  // fills the coordinates option menu
      g.fill_renumber_residue_range_internal(w, imol); // fills the chain option menu
   }
   return w;
}

void renumber_residues_from_widget(GtkWidget *window) {

   int imol = graphics_info_t::renumber_residue_range_molecule;

   GtkWidget *e1 = lookup_widget(window, "renumber_residue_range_resno_1_entry");
   GtkWidget *e2 = lookup_widget(window, "renumber_residue_range_resno_2_entry");
   GtkWidget *offent = lookup_widget(window, "renumber_residue_range_offset_entry");
   

   std::pair<short int, int> r1  = int_from_entry(e1);
   std::pair<short int, int> r2  = int_from_entry(e2);
   std::pair<short int, int> off = int_from_entry(offent);

   if (r1.first && r2.first && off.first) {
      int start_res = r1.second;
      int last_res =  r2.second;
      int offset   = off.second;

      if (imol >= 0) {
	 if (imol < graphics_info_t::n_molecules()) {
	    if (graphics_info_t::molecules[imol].has_model()) {
	       std::string chain = graphics_info_t::renumber_residue_range_chain;
	       
	       renumber_residue_range(imol, chain.c_str(),
				      start_res, last_res, offset);
	    }
	 }
      }
   } else {
      std::cout << "Sorry. Couldn't read residue or offset from entry widget\n";
   } 
}

/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */

void execute_refmac(GtkWidget *window) {  /* lookup stuff here. */

   // The passed window, is the refmac dialog, where one selects the
   // coords molecule and the map molecule.

   GtkWidget *option_menu = lookup_widget(window,
					  "run_refmac_coords_optionmenu");

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));

   GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));

   int imol_coords = graphics_info_t::refmac_molecule; 
   if (imol_coords < 0) {
      std::cout << "No coordinates selected for refmac\n";
   } else { 

//       std::cout << " Running refmac coords molecule number "
// 		<< imol_coords << std::endl;

      option_menu = lookup_widget(window, "run_refmac_map_optionmenu");
      GtkWidget *mtz_file_radiobutton;
#if (GTK_MAJOR_VERSION > 1)
      mtz_file_radiobutton = lookup_widget(window, "run_refmac_mtz_file_radiobutton");
#else
      mtz_file_radiobutton = NULL;
#endif
      menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
      active_item = gtk_menu_get_active(GTK_MENU(menu));
      int imol_map_refmac = -1;
      int have_mtz_file = 0;
      if (mtz_file_radiobutton && GTK_TOGGLE_BUTTON(mtz_file_radiobutton)->active) {
	have_mtz_file = 1;
      }

      // active_item is set if there was at least one map with refmac params:
      // if none, it is null.
      
      if (active_item == 0 && have_mtz_file == 0) {
	 add_status_bar_text("No map has associated Refmac Parameters - no REFMAC!");
      } else {
	int imol_window = -1;
	std::string mtz_in_filename = "";
	if (!have_mtz_file) {
	  // we get imol from a map mtz file
	   imol_window = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(active_item)));
	} else {
	   // check the filename of the button
	   GtkWidget *button_mtz_label = lookup_widget(window, "run_refmac_mtz_file_label");
#if (GTK_MAJOR_VERSION > 1)
	   const gchar *mtz_filename = gtk_label_get_text(GTK_LABEL(button_mtz_label));
	   mtz_in_filename = mtz_filename;
#else
	   gchar **mtz_filename = 0;
	   gtk_label_get(GTK_LABEL(button_mtz_label), mtz_filename);
	   mtz_in_filename = (char *)mtz_filename;
#endif // GTK
	   if (mtz_in_filename == "(None)") {
	      have_mtz_file = 0;
	      std::cout << "WARNING:: no mtz file given" <<std::endl;
	   } else {
	      if (! coot::file_exists(mtz_in_filename)) {
		 have_mtz_file = 0;
		 std::cout << "WARNING:: mtz file " << mtz_in_filename << " does not exist" <<std::endl;
	      }
	   }
	}

	if (imol_window < 0 && have_mtz_file == 0) {
	  if (have_mtz_file == 0) {
	    std::cout << "No mtz file selected for refmac\n";
	  } else {
	    std::cout << "No map data selected for refmac\n";
	  }
	} else { 
	  imol_map_refmac = imol_window;
	  if (!is_valid_map_molecule(imol_map_refmac) && have_mtz_file == 0) {
	    std::string s = "Invalid molecule number: ";
	    s += graphics_info_t::int_to_string(imol_map_refmac);
	    std::cout << s << std::endl;
	    graphics_info_t g;
	    g.add_status_bar_text(s);
	  } else {
	    // normal path
	    //	       if (graphics_info_t::molecules[imol_map_refmac].Have_sensible_refmac_params()) { 
	    // just check for refmac mtz file now (either from map or direct
	    if (have_mtz_file == 1 ||
		graphics_info_t::molecules[imol_map_refmac].Refmac_mtz_filename().size() > 0) {
	      if (!have_mtz_file) {
		std::cout << " Running refmac refmac params molecule number "
			  << imol_map_refmac << std::endl;
	      } else {
		std::cout << " Running refmac from mtz file " << mtz_in_filename << "(not map)" << std::endl;
	      }

	      graphics_info_t g;

	      std::string refmac_dir("coot-refmac");
	      short int have_ccp4i_project = 0;
	      if (graphics_info_t::refmac_ccp4i_project_dir != "") { 
		refmac_dir = graphics_info_t::refmac_ccp4i_project_dir;
		have_ccp4i_project = 1;
	      }
	      int istat = make_directory_maybe(refmac_dir.c_str());
	      if (istat != 0) { // fails
		std::cout << "WARNING failed to make directory for refmac -"
			  << " run refmac fails\n" << std::endl;
	      } else {
		// now lookup the active state of the difference map and
		// the phase combine buttons:
		//
		int diff_map_flag;
		int phase_combine_flag;
		GtkWidget *checkbutton;
		// phase_combine_flag is set in refmac_phase_input now.
		// 0: no phase
		// 1: combine (with phase and FOM, as before)
		// 2: combine with HL
		// 3: SAD
		phase_combine_flag = get_refmac_phase_input();
	       
		checkbutton =  lookup_widget(window,"run_refmac_diff_map_checkbutton");
		if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
		  diff_map_flag = 1;
		} else {
		  diff_map_flag = 0;
		}

		// g.molecules[imol_coords].increment_refmac_count();
      
		std::string pdb_in_filename  = refmac_dir;
		std::string pdb_out_filename = refmac_dir;
		std::string mtz_out_filename = refmac_dir;
		std::cout << "DEBUG:: pdb_in_filename is now 1 " <<  pdb_in_filename << std::endl;
		if (! have_ccp4i_project) { 
		  pdb_in_filename  += "/";
		  pdb_out_filename += "/";
		  mtz_out_filename += "/";
		}
		std::cout << "DEBUG:: pdb_in_filename is now 2 " <<  pdb_in_filename << std::endl;
		pdb_in_filename += g.molecules[imol_coords].Refmac_in_name();
		std::cout << "DEBUG:: pdb_in_filename is now 3 " <<  pdb_in_filename << std::endl;

		// cleverness happens in Refmac_out_name:
		pdb_out_filename += g.molecules[imol_coords].Refmac_out_name();
		mtz_out_filename += g.molecules[imol_coords].Refmac_mtz_out_name();

		if (! have_mtz_file) {
		  mtz_in_filename = g.molecules[imol_map_refmac].Refmac_mtz_filename();
		}
		std::string refmac_count_string =
		  g.int_to_string(g.molecules[imol_coords].Refmac_count());

		std::cout << "DEBUG:: mtz_out_filename: " <<
		  mtz_out_filename << std::endl;
		std::cout << "DEBUG:: pdb_out_filename: " <<
		  pdb_out_filename << std::endl;

		// now get the column labels

		// before running refmac we may want to set refmac parameters from the GUI
		// this should overwrite whatever has been set as refmac parameters before
		// we do it before checking for phases, so that these can be included later
		coot::mtz_column_types_info_t *saved_f_phi_columns
		  = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
		     
		std::string phib_string = "";
		std::string fom_string  = "";
		std::string fobs_col;
		std::string sigfobs_col;

		int icol;
		int sensible_r_free_col = 0;
		std::string fiobs_col;
		std::string sigfiobs_col;
		std::string r_free_col;

		// different in GTK1/2
#if (GTK_MAJOR_VERSION > 1)

		// and we need to get the column labels if we have an input mtz (not map)
		if (have_mtz_file) {

		  if (refmac_use_twin_state() == 0) {
		    // for now we only use Is in twin not in 'normal' refinement
		      icol = saved_f_phi_columns->selected_refmac_fobs_col;
		      fobs_col = saved_f_phi_columns->f_cols[icol].column_label;
		  } else {
		    // for twin we check both Is and Fs
		    // first check the I
		    icol = saved_f_phi_columns->selected_refmac_iobs_col;
		    if (icol > -1) {
		      fiobs_col = saved_f_phi_columns->i_cols[icol].column_label;
		      set_refmac_use_intensities(1);
		    } else {
		      // we must have F/sigF then
		      icol = saved_f_phi_columns->selected_refmac_fobs_col;
		      fiobs_col = saved_f_phi_columns->f_cols[icol].column_label;
		      set_refmac_use_intensities(0);
		    }
		    icol = saved_f_phi_columns->selected_refmac_sigfobs_col;
		    sigfiobs_col = saved_f_phi_columns->sigf_cols[icol].column_label;
		  }
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

		std::string phi_label = "";
		std::string fom_label = "";
		std::string hla_label = "";
		std::string hlb_label = "";
		std::string hlc_label = "";
		std::string hld_label = "";

		if (saved_f_phi_columns->selected_refmac_fobs_col > -1 &&
		    saved_f_phi_columns->selected_refmac_sigfobs_col > -1) {
		  // assign the labels if we have them (from the GUI)
		  //std::cout << "BL DEBUG:: selected fobs col " << saved_f_phi_columns->selected_refmac_fobs_col<<std::endl;
		  //std::cout << "BL DEBUG:: selected sigf col " << saved_f_phi_columns->selected_refmac_sigfobs_col<<std::endl;
		  //std::cout << "BL DEBUG:: selected rfree col " << saved_f_phi_columns->selected_refmac_r_free_col<<std::endl;
		  //std::cout << "BL DEBUG:: selected phi col " << saved_f_phi_columns->selected_refmac_phi_col<<std::endl;
		  //std::cout << "BL DEBUG:: selected fom col " << saved_f_phi_columns->selected_refmac_fom_col<<std::endl;
		  //std::cout << "BL DEBUG:: selected hla col " << saved_f_phi_columns->selected_refmac_hla_col<<std::endl;

		  if (phase_combine_flag == 3) {
		    // SAD F/sigF columns
		    // we make F=/F- and sigF+/f- a list to pass to scripting refmac
		    std::string fp_col;
		    std::string fm_col;
		    std::string sigfp_col;
		    std::string sigfm_col;
		    // F+
		    icol = saved_f_phi_columns->selected_refmac_fp_col;
		    fp_col = saved_f_phi_columns->fpm_cols[icol].column_label;
		    // F-
		    icol = saved_f_phi_columns->selected_refmac_fm_col;
		    fm_col = saved_f_phi_columns->fpm_cols[icol].column_label;
		    // sigF+
		    icol = saved_f_phi_columns->selected_refmac_sigfp_col;
		    sigfp_col = saved_f_phi_columns->sigfpm_cols[icol].column_label;
		    // sigF-
		    icol = saved_f_phi_columns->selected_refmac_sigfm_col;
		    sigfm_col = saved_f_phi_columns->sigfpm_cols[icol].column_label;
		    // make lists
#ifdef USE_GUILE
		    fobs_col  = "(cons ";
		    fobs_col += single_quote(fp_col);
		    fobs_col += " ";
		    fobs_col += single_quote(fm_col);
		    fobs_col += ")";
		    sigfobs_col  = "(cons ";
		    sigfobs_col += single_quote(sigfp_col);
		    sigfobs_col += " ";
		    sigfobs_col += single_quote(sigfm_col);
		    sigfobs_col += ")";
#else
#ifdef USE_PYTHON
		    fobs_col = "[";
		    fobs_col += single_quote(fp_col);
		    fobs_col += ", ";
		    fobs_col += single_quote(fm_col);
		    fobs_col += "]";
		    sigfobs_col = "[";
		    sigfobs_col += single_quote(sigfp_col);
		    sigfobs_col += ", ";
		    sigfobs_col += single_quote(sigfm_col);
		    sigfobs_col += "]";
#endif // USE_GUILE
#endif // USE_PYTHON
		    // now get the information about anomalous atom
		    GtkWidget *atom_entry    = lookup_widget(window, "run_refmac_sad_atom_entry");
		    GtkWidget *fp_entry      = lookup_widget(window, "run_refmac_sad_fp_entry");
		    GtkWidget *fpp_entry     = lookup_widget(window, "run_refmac_sad_fpp_entry");
		    GtkWidget *lambda_entry  = lookup_widget(window, "run_refmac_sad_lambda_entry");
		    const gchar *atom_str  = gtk_entry_get_text(GTK_ENTRY(atom_entry));
		    std::string fp_str     = gtk_entry_get_text(GTK_ENTRY(fp_entry));
		    std::string fpp_str    = gtk_entry_get_text(GTK_ENTRY(fpp_entry));
		    std::string lambda_str = gtk_entry_get_text(GTK_ENTRY(lambda_entry));
		    float fp, fpp, lambda;
		    if (fp_str != "") {
		      fp = atof(fp_str.c_str());
		    } else {
		      fp = -9999; // magic unset
		    }
		    if (fpp_str != "") {
		      fpp = atof(fpp_str.c_str());
		    } else {
		      fpp = -9999; // magic unset
		    }
		    if (lambda_str != "") {
		      lambda = atof(lambda_str.c_str());
		    } else {
		      lambda = -9999; // magic unset
		    }
		    add_refmac_sad_atom(atom_str, fp, fpp, lambda);

		  } else {
		    icol = saved_f_phi_columns->selected_refmac_fobs_col;
		    fobs_col = saved_f_phi_columns->f_cols[icol].column_label;

		    icol = saved_f_phi_columns->selected_refmac_sigfobs_col;
		    sigfobs_col = saved_f_phi_columns->sigf_cols[icol].column_label;
		  }

		  icol = saved_f_phi_columns->selected_refmac_r_free_col; /* magic -1 if not set */
		  if (icol >= 0) { 
		    // 
		    sensible_r_free_col = 1;
		    r_free_col = saved_f_phi_columns->r_free_cols[icol].column_label;
		  } else { 
		    sensible_r_free_col = 0;
		    r_free_col = "";
		  }

		  // We save the phase and FOM as 'fourier_*_labels' too, so that they are saved!?
		  if (phase_combine_flag == 1) {
		    icol = saved_f_phi_columns->selected_refmac_phi_col;
		    if (icol == -1) { 
		      printf("INFO:: no phase available (phi/fom)! \n");
		    } else { 
		      phi_label = saved_f_phi_columns->phi_cols[icol].column_label; 
		      icol = saved_f_phi_columns->selected_refmac_fom_col;
		      fom_label = saved_f_phi_columns->weight_cols[icol].column_label;
		      if (! have_mtz_file) {
			graphics_info_t::molecules[imol_map_refmac].store_refmac_phase_params(std::string(phi_label),
											      std::string(fom_label),
											      std::string(hla_label),
											      std::string(hlb_label),
											      std::string(hlc_label),
											      std::string(hld_label));
		      }
		    }
		  }

		  // check the HLs
		  if (phase_combine_flag == 2) {
		    icol = saved_f_phi_columns->selected_refmac_hla_col;
		    if (icol == -1) {
		      printf("INFO:: no phase available (HLs)! \n");
		    } else { 
		      hla_label = saved_f_phi_columns->hl_cols[icol].column_label;
		      icol = saved_f_phi_columns->selected_refmac_hlb_col;
		      hlb_label = saved_f_phi_columns->hl_cols[icol].column_label;
		      icol = saved_f_phi_columns->selected_refmac_hlc_col;
		      hlc_label = saved_f_phi_columns->hl_cols[icol].column_label;
		      icol = saved_f_phi_columns->selected_refmac_hld_col;
		      hld_label = saved_f_phi_columns->hl_cols[icol].column_label;
		      g_print("BL DEBUG:: have HLs \n");
		      if (! have_mtz_file) {
			graphics_info_t::molecules[imol_map_refmac].store_refmac_phase_params(std::string(phi_label),
											      std::string(fom_label),
											      std::string(hla_label),
											      std::string(hlb_label),
											      std::string(hlc_label),
											      std::string(hld_label));
		      }
		    }
		  }
			    
		  if (have_mtz_file){
		    graphics_info_t g;
		    g.store_refmac_params(std::string(mtz_in_filename),
					  std::string(fobs_col), 
					  std::string(sigfobs_col), 
					  std::string(r_free_col),
					  sensible_r_free_col);
		    set_refmac_used_mtz_file(1);
		  } else {
		    graphics_info_t::molecules[imol_map_refmac].store_refmac_params(std::string(mtz_in_filename),
										    std::string(fobs_col), 
										    std::string(sigfobs_col), 
										    std::string(r_free_col),
										    sensible_r_free_col);
		    set_refmac_used_mtz_file(0);
		  }
		}

		//if (g.molecules[imol_map_refmac].Fourier_weight_label() != "") {
		//  phib_string = g.molecules[imol_map_refmac].Fourier_phi_label();
		//  fom_string  = g.molecules[imol_map_refmac].Fourier_weight_label();
		//} else {
		if (phase_combine_flag == 1) {
		  if (! have_mtz_file) {
		    if (g.molecules[imol_map_refmac].Refmac_phi_col() != "") {
		      phib_string = g.molecules[imol_map_refmac].Refmac_phi_col();
		      fom_string  = g.molecules[imol_map_refmac].Refmac_fom_col();
		    } else {
		      std::cout << "WARNING:: Can't do phase combination if we don't use FOMs ";
		      std::cout << "to make the map" << std::endl;
		      std::cout << "WARNING:: Turning off phase combination." << std::endl;
		      phase_combine_flag = 0;
		    }
		  } else {
		    if (phi_label != "" && fom_label != "") {
		      phib_string = phi_label;
		      fom_string  = fom_label;
		    } else {
		      std::cout << "WARNING:: Can't do phase combination if we don't use FOMs ";
		      std::cout << "to make the map" << std::endl;
		      std::cout << "WARNING:: Turning off phase combination." << std::endl;
		      phase_combine_flag = 0;		      
		    }
		  }
		}
		// 	    std::cout << "DEBUG:: fom_string " << fom_string << " "
		// 		      << g.molecules[imol_map_refmac].Fourier_weight_label()
		// 		      << std::endl;

		// now check for HLs
		if (phase_combine_flag == 2) {
		  std::string hla_string = "";
		  std::string hlb_string;
		  std::string hlc_string;
		  std::string hld_string;
		  if (! have_mtz_file && g.molecules[imol_map_refmac].Refmac_hla_col() != "") {
		    hla_string = g.molecules[imol_map_refmac].Refmac_hla_col();
		    hlb_string = g.molecules[imol_map_refmac].Refmac_hlb_col();
		    hlc_string = g.molecules[imol_map_refmac].Refmac_hlc_col();
		    hld_string = g.molecules[imol_map_refmac].Refmac_hld_col();
		  } else {
		    if (have_mtz_file && hla_label != ""){
		      hla_string = hla_label;
		      hlb_string = hlb_label;
		      hlc_string = hlc_label;
		      hld_string = hld_label;
		    } else {
		      std::cout << "WARNING:: no valid HL columns found" <<std::endl;
		    }
		  }
		  if (hla_string != "") {
		    // now save the HLs in a list string (phib) for refmac, depending on scripting
		    std::vector<std::string> hl_list;
		    hl_list.push_back(hla_string);
		    hl_list.push_back(hlb_string);
		    hl_list.push_back(hlc_string);
		    hl_list.push_back(hld_string);
#ifdef USE_GUILE
		    //phib_string = "(list ";
		    // here we pass it as a string, scheme will later make a list
		    phib_string  = hla_string;
		    phib_string += " ";
		    phib_string += hlb_string;
		    phib_string += " ";
		    phib_string += hlc_string;
		    phib_string += " ";
		    phib_string += hld_string;
		    //phib_string += ")";
#else
#ifdef USE_PYTHON
		    phib_string = "[";
		    phib_string += single_quote(hla_string);
		    phib_string += ", ";
		    phib_string += single_quote(hlb_string);
		    phib_string += ", ";
		    phib_string += single_quote(hlc_string);
		    phib_string += ", ";
		    phib_string += single_quote(hld_string);
		    phib_string += "]";
#endif // USE_GUILE
#endif // USE_PYTHON
		    fom_string = "";
		  } else {
		    std::cout << "WARNING:: Can't do phase combination if we don't have HLs " << std::endl;
		    std::cout << "WARNING:: Turning off phase combination." << std::endl;
		    phase_combine_flag = 0;
		  }
		}
		// for TWIN we reset the flags as we dont have phase combination for twin yet
		if (refmac_use_twin_state() == 1) {
		  phase_combine_flag = 0;
		  phib_string = "";
		  fom_string = "";
		}

#else
		// now gtk1
		if (graphics_info_t::molecules[imol_map_refmac].Have_sensible_refmac_params()) {

		  std::cout << " Running refmac refmac params molecule number "
			    << imol_map_refmac << std::endl;
		  
		  fobs_col = g.molecules[imol_map_refmac].Refmac_fobs_col();
		  sigfobs_col = g.molecules[imol_map_refmac].Refmac_sigfobs_col();
		  r_free_col = g.molecules[imol_map_refmac].Refmac_r_free_col();
		  sensible_r_free_col = g.molecules[imol_map_refmac].Refmac_r_free_sensible();

		}
#endif //GTK1/2
		     
		std::string cif_lib_filename = ""; // default, none
		if (graphics_info_t::cif_dictionary_filename_vec->size() > 0) {
		  cif_lib_filename = (*graphics_info_t::cif_dictionary_filename_vec)[0];
		}


		// 	    std::cout << "DEBUG:: attempting to write pdb input file "
		// 		      << pdb_in_filename << std::endl;
		int ierr = g.molecules[imol_coords].write_pdb_file(pdb_in_filename);
		if (!ierr) { 
		  std::cout << "refmac ccp4i project dir " 
			    << graphics_info_t::refmac_ccp4i_project_dir 
			    << std::endl;
		  int run_refmac_with_no_labels = 0;
#if (GTK_MAJOR_VERSION > 1)
		  if (refmac_runs_with_nolabels()) {
		    GtkWidget *nolabels_checkbutton = lookup_widget(window,
								    "run_refmac_nolabels_checkbutton");
		    if (GTK_TOGGLE_BUTTON(nolabels_checkbutton)->active) {
		      run_refmac_with_no_labels = 1;
		      fobs_col    = "";
		      sigfobs_col = "";
		      r_free_col  = "";
		      sensible_r_free_col = 0;
		    }
		  }

#endif // GTK2

		  // And finally run refmac
		  if (run_refmac_with_no_labels == 1 || fobs_col != "") {
		     short int make_molecules_flag = 1; // not a sub-thread, (so do things
		                                        // the normal/old way).
		     execute_refmac_real(pdb_in_filename, pdb_out_filename,
					 mtz_in_filename, mtz_out_filename,
					 cif_lib_filename,
					 fobs_col, sigfobs_col, r_free_col, sensible_r_free_col,
					 make_molecules_flag,
					 refmac_count_string,
					 g.swap_pre_post_refmac_map_colours_flag,
					 imol_map_refmac,
					 diff_map_flag,
					 phase_combine_flag, phib_string, fom_string,
					 graphics_info_t::refmac_ccp4i_project_dir);
		  } else {

		    std::cout << "WARNING:: we cannot run Refmac without without valid labels" <<std::endl;
		  }
		} else {
		  std::cout << "WARNING:: fatal error in writing pdb input file"
			    << pdb_in_filename << " for refmac.  Can't run refmac"
			    << std::endl;
		}
	      }
	    }
	  }
	}
	
      }
   }
}


void set_refmac_phase_input(int phase_flag) {

  graphics_info_t g;
#if (GTK_MAJOR_VERSION > 1)
  g.set_refmac_phase_input(phase_flag);
#else
  std::cout<<"BL INFO:: cannot use phase recombination in GTK1 (yet), so wont change flag!\n"<<std::endl;;
  g.set_refmac_phase_input(0);
#endif // GTK2
}

void set_refmac_use_sad(int state) {

  graphics_info_t g;
#if (GTK_MAJOR_VERSION > 1)
  g.set_refmac_use_sad(state);
#else
  std::cout<<"BL INFO:: cannot use SAD refinement in GTK1 (yet), so wont change flag!\n"<<std::endl;;
  g.set_refmac_use_sad(0);
#endif // GTK2
}


void fill_option_menu_with_refmac_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_options(optionmenu);

} 

void fill_option_menu_with_refmac_methods_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_methods_options(optionmenu);

} 

void fill_option_menu_with_refmac_phase_input_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_phase_input_options(optionmenu);

} 

void fill_option_menu_with_refmac_labels_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_labels_options(optionmenu);

}

void fill_option_menu_with_refmac_file_labels_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_file_labels_options(optionmenu);

} 

void fill_option_menu_with_refmac_ncycle_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_ncycle_options(optionmenu);

}

void update_refmac_column_labels_frame(GtkWidget *optionmenu, 
				       GtkWidget *fobs_menu, GtkWidget *fiobs_menu, GtkWidget *fpm_menu,
				       GtkWidget *r_free_menu,
				       GtkWidget *phases_menu, GtkWidget *fom_menu, GtkWidget *hl_menu) {
  graphics_info_t g;
  g.update_refmac_column_labels_frame(optionmenu,
				      fobs_menu, fiobs_menu, fpm_menu,
				      r_free_menu,
				      phases_menu, fom_menu, hl_menu);

}

void
fill_refmac_sad_atom_entry(GtkWidget *w) {

  GtkWidget *atom_entry    = lookup_widget(w, "run_refmac_sad_atom_entry");
  GtkWidget *fp_entry      = lookup_widget(w, "run_refmac_sad_fp_entry");
  GtkWidget *fpp_entry     = lookup_widget(w, "run_refmac_sad_fpp_entry");
  GtkWidget *lambda_entry  = lookup_widget(w, "run_refmac_sad_lambda_entry");
  if (graphics_info_t::refmac_sad_atoms.size() > 0) {
    std::string atom_name  = graphics_info_t::refmac_sad_atoms[0].atom_name;
    float fp     = graphics_info_t::refmac_sad_atoms[0].fp;
    float fpp    = graphics_info_t::refmac_sad_atoms[0].fpp;
    float lambda = graphics_info_t::refmac_sad_atoms[0].lambda;
    std::string fp_str = "";
    std::string fpp_str = "";
    std::string lambda_str = "";
    if (fabs(fp + 9999) >= 0.1) {
      fp_str = graphics_info_t::float_to_string(fp);
    } 
    if (fabs(fpp + 9999) >= 0.1) {
      fpp_str = graphics_info_t::float_to_string(fpp);
    }
    if (fabs(lambda + 9999) >= 0.1) {
      lambda_str = graphics_info_t::float_to_string(lambda);
    }
    gtk_entry_set_text(GTK_ENTRY(atom_entry), atom_name.c_str());
    gtk_entry_set_text(GTK_ENTRY(fp_entry), fp_str.c_str());
    gtk_entry_set_text(GTK_ENTRY(fpp_entry), fpp_str.c_str());
    gtk_entry_set_text(GTK_ENTRY(lambda_entry), lambda_str.c_str());
  }
}


void apply_add_OXT_from_widget(GtkWidget *w) {

   int imol = graphics_info_t::add_OXT_molecule;
   int resno = -9999;
   std::string chain_id = graphics_info_t::add_OXT_chain;

   GtkWidget *terminal_checkbutton = lookup_widget(w, "add_OXT_c_terminus_radiobutton");
   GtkWidget *residue_number_entry = lookup_widget(w, "add_OXT_residue_entry");

   if (GTK_TOGGLE_BUTTON(terminal_checkbutton)->active) {
      std::cout << "DEBUG:: auto determine C terminus..." << std::endl;
      // we need to determine the last residue in this chain:
      if (is_valid_model_molecule(imol)) { 
	 if (graphics_info_t::molecules[imol].has_model()) {
	    std::pair<short int, int> p =
	       graphics_info_t::molecules[imol].last_residue_in_chain(chain_id);
	    if (p.first) {
	       resno = p.second;
	    } 
	 }
      }
   } else {
      // we get the resno from the widget
      std::pair<short int, int> p = int_from_entry(residue_number_entry);
      if (p.first) {
	 resno = p.second;
      }
   }

   if (resno > -9999) { 
      if (is_valid_model_molecule(imol)) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    std::cout << "DEBUG:: adding OXT to " << imol << " "
		      << chain_id << " " << resno << std::endl;
	    
	    add_OXT_to_residue(imol, resno, "", chain_id.c_str());
	 }
      }
   } else {
      std::cout << "WARNING:: Could not determine last residue - not adding OXT\n";
   } 
}


GtkWidget *wrapped_create_add_OXT_dialog() {

   GtkWidget *w = create_add_OXT_dialog();

   GtkWidget *option_menu = lookup_widget(w, "add_OXT_molecule_optionmenu");

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(graphics_info_t::add_OXT_molecule_item_select);

   graphics_info_t g;
   int imol = first_coords_imol();
   graphics_info_t::add_OXT_molecule = imol;
   g.fill_option_menu_with_coordinates_options(option_menu, callback_func, imol);

   g.fill_add_OXT_dialog_internal(w, imol); // function needs object (not static)

   return w;
}

void setup_alt_conf_with_dialog(GtkWidget *dialog) {

   GtkWidget *widget_ca = lookup_widget(dialog, 
					"add_alt_conf_ca_radiobutton");
   GtkWidget *widget_whole = lookup_widget(dialog, 
					   "add_alt_conf_whole_single_residue_radiobutton");
   GtkWidget *widget_range = lookup_widget(dialog, 
					   "add_alt_conf_residue_range_radiobutton");

   if (graphics_info_t::alt_conf_split_type_number() == 0)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget_ca), TRUE);
   if (graphics_info_t::alt_conf_split_type_number() == 1)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget_whole), TRUE);
   if (graphics_info_t::alt_conf_split_type_number() == 2)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget_range), TRUE);

   if (graphics_info_t::alt_conf_split_type_number() < 2) { 
      std::cout << "Click on the residue you want to split" << std::endl;
   } else { 
      std::cout << "Click on a residue range you want to split" << std::endl;
   }

   graphics_info_t::in_add_alt_conf_define = 1;
   graphics_info_t::pick_cursor_maybe();
   graphics_info_t::pick_pending_flag = 1;
   graphics_info_t::add_alt_conf_dialog = dialog;
} 


void altconf() { 

   GtkWidget *widget = create_add_alt_conf_dialog();
   setup_alt_conf_with_dialog(widget);
   gtk_widget_show(widget);
}

/*  ------------------------------------------------------------------------ */
/*                         recover session:                                  */
/*  ------------------------------------------------------------------------ */
/* section Recover Session Function */
/* After a crash (shock horror!) we provide this convenient interface
   to restore the session.  It runs through all the molecules with
   models and looks at the coot backup directory looking for related
   backup files that are more recent that the read file. */
void recover_session() { 

   int i_rec = 0;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) { 
      if (graphics_info_t::molecules[imol].has_model()) { 
	 coot::backup_file_info info = 
	    graphics_info_t::molecules[imol].recent_backup_file_info();
	 if (info.status) { 

	    coot::backup_file_info *info_copy = new coot::backup_file_info;
	    *info_copy = info;
	    info_copy->imol = imol;
	    
	    GtkWidget *widget = create_recover_coordinates_dialog();
	    gtk_object_set_user_data(GTK_OBJECT(widget), info_copy);
	    
	    GtkWidget *label1, *label2;
	    label1 = lookup_widget(widget, "recover_coordinates_read_coords_label");
	    label2 = lookup_widget(widget, "recover_coordinates_backup_coordinates_label");

	    gtk_label_set_text(GTK_LABEL(label1), info.name.c_str());
	    gtk_label_set_text(GTK_LABEL(label2), info.backup_file_name.c_str());

	    gtk_widget_show(widget);
	    i_rec++;
	 }
      }
   }
   if (i_rec == 0) {
      GtkWidget *w = create_nothing_to_recover_dialog();
      gtk_widget_show(w);
   }
}

// widget needed for lookup of user data:
// 
void execute_recover_session(GtkWidget *widget) { 

   coot::backup_file_info *info = (coot::backup_file_info *) gtk_object_get_user_data(GTK_OBJECT(widget));

   if (info) { 
      
      graphics_info_t g;
      if (info->imol >= 0 && info->imol < g.n_molecules()) {
	 std::string cwd = coot::util::current_working_dir();
	 g.molecules[info->imol].execute_restore_from_recent_backup(info->backup_file_name, cwd);
	 graphics_draw();
      }
   } else { 
      std::cout << "ERROR:: couldn't find user data in execute_recover_session\n";
   } 
} 


/*  ----------------------------------------------------------------------- */
/*                         Merge Molecules                                  */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_merge_molecules_dialog() {

   GtkWidget *w = create_merge_molecules_dialog();
   // fill the dialog here
   GtkWidget *molecule_option_menu = lookup_widget(w, "merge_molecules_optionmenu");
   GtkWidget *molecules_vbox       = lookup_widget(w, "merge_molecules_vbox");

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(merge_molecules_menu_item_activate);
   GtkSignalFunc checkbox_callback_func = GTK_SIGNAL_FUNC(on_merge_molecules_check_button_toggled);


   fill_vbox_with_coordinates_options(molecules_vbox, checkbox_callback_func);

   int imol_master = graphics_info_t::merge_molecules_master_molecule;
   if (imol_master == -1) { 
      for (int i=0; i<graphics_info_t::n_molecules(); i++) {
	 if (graphics_info_t::molecules[i].has_model()) {
	    graphics_info_t::merge_molecules_master_molecule = i;
	    imol_master = i;
	    break;
	 }
      }
   }

   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(molecule_option_menu,
					       callback_func, imol_master);
   return w;
}

void merge_molecules_menu_item_activate(GtkWidget *item, 
					GtkPositionType pos) {

   graphics_info_t::merge_molecules_master_molecule = pos;
}

void fill_vbox_with_coordinates_options(GtkWidget *dialog,
					GtkSignalFunc checkbox_callback_func) {

   GtkWidget *checkbutton;
   std::string button_label;
   GtkWidget *molecules_vbox = lookup_widget(dialog, "merge_molecules_vbox");

   // Unset any preconcieved notion of merging molecules:
   // 
   graphics_info_t::merge_molecules_merging_molecules->resize(0);

   for (int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 button_label = graphics_info_t::int_to_string(i);
	 button_label += " ";
	 button_label += graphics_info_t::molecules[i].name_for_display_manager();
	 std::string button_name = "merge_molecules_checkbutton_";
	 button_name += graphics_info_t::int_to_string(i);

	 checkbutton = gtk_check_button_new_with_label(button_label.c_str());
  	 gtk_widget_ref (checkbutton);
  	 gtk_object_set_data_full (GTK_OBJECT (dialog),
  				   button_name.c_str(), checkbutton,
  				   (GtkDestroyNotify) gtk_widget_unref);
	 // The callback (if active) adds this molecule to the merging molecules list.
	 // If not active, it tries to remove it from the list.
	 //
	 // Why am I doing it like this instead of just looking at the
	 // state of the checkbutton when the OK button is pressed?
	 // Because (for the life of me) I can't seem to correctly
	 // lookup the checkbuttons from the button (or dialog for
	 // that matter).
	 // 
	 //  We look at the state when the
	 // "Merge" button is pressed - we don't need a callback to do
	 // that.
	 // 
  	 gtk_signal_connect (GTK_OBJECT (checkbutton), "toggled",
  			     GTK_SIGNAL_FUNC (checkbox_callback_func),
  			     GINT_TO_POINTER(i));
	 gtk_widget_show (checkbutton);
	 gtk_box_pack_start (GTK_BOX (molecules_vbox), checkbutton, FALSE, FALSE, 0);
	 gtk_container_set_border_width (GTK_CONTAINER (checkbutton), 2);
      }
   }
}

// The callback (if active) adds this molecule to the merging molecules list.
// If not active, it tries to remove it from the list.
// 
void on_merge_molecules_check_button_toggled (GtkToggleButton *togglebutton,
					      gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   if (togglebutton->active) {
      std::cout << "INFO:: adding molecule " << imol << " to merging list\n";
      graphics_info_t::merge_molecules_merging_molecules->push_back(imol);
   } else {
      std::cout << "INFO:: removing molecule " << imol << " from merging list\n";
      if (coot::is_member_p(*graphics_info_t::merge_molecules_merging_molecules, imol)) {
	 // passing a pointer
	 coot::remove_member(graphics_info_t::merge_molecules_merging_molecules, imol);
      }
   }
}


// Display the gui
void do_merge_molecules_gui() {

   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   gtk_widget_show(w);
} 

// The action on Merge button press:
// 
void do_merge_molecules(GtkWidget *dialog) {

   std::vector<int> add_molecules = *graphics_info_t::merge_molecules_merging_molecules;
   if (add_molecules.size() > 0) { 
      std::pair<int, std::vector<std::string> > stat =
	 merge_molecules_by_vector(add_molecules, graphics_info_t::merge_molecules_master_molecule);
      if (stat.first)
	 graphics_draw();
   }
}

/*  ----------------------------------------------------------------------- */
/*                         Mutate Sequence GUI                              */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_mutate_sequence_dialog() {

   GtkWidget *w = create_mutate_sequence_dialog();

   set_transient_and_position(COOT_MUTATE_RESIDUE_RANGE_WINDOW, w);

   GtkWidget *molecule_option_menu = lookup_widget(w, "mutate_molecule_optionmenu");
   GtkWidget *chain_option_menu    = lookup_widget(w, "mutate_molecule_chain_optionmenu");
//    GtkWidget *entry1 = lookup_widget(w, "mutate_molecule_resno_1_entry");
//    GtkWidget *entry2 = lookup_widget(w, "mutate_molecule_resno_2_entry");
//    GtkWidget *textwindow = lookup_widget(w, "mutate_molecule_sequence_text");
   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(mutate_sequence_molecule_menu_item_activate);


   // Get the default molecule and fill chain optionmenu with the molecules chains:
   int imol = -1; 
   for (int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }
   if (imol >= 0) {
      graphics_info_t::mutate_sequence_imol = imol;
      GtkSignalFunc callback =
	 GTK_SIGNAL_FUNC(mutate_sequence_chain_option_menu_item_activate);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu, imol,
								      callback);
      graphics_info_t::mutate_sequence_chain_from_optionmenu = set_chain;
   } else {
      graphics_info_t::mutate_sequence_imol = -1; // flag for can't mutate
   }
   graphics_info_t g;
   // std::cout << "DEBUG:: filling option menu with default molecule " << imol << std::endl;
   g.fill_option_menu_with_coordinates_options(molecule_option_menu, callback_func, imol);
   return w;
}

void mutate_sequence_molecule_menu_item_activate(GtkWidget *item, 
						 GtkPositionType pos) {

   // change the chain id option menu here...
   std::cout << "DEBUG:: mutate_sequence_molecule_menu_item_activate got pos:"
	     << pos << std::endl;

   graphics_info_t::mutate_sequence_imol = pos;

   GtkWidget *chain_option_menu =
      lookup_widget(item, "mutate_molecule_chain_optionmenu");

   GtkSignalFunc callback_func =
      GTK_SIGNAL_FUNC(mutate_sequence_chain_option_menu_item_activate);
   
   std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu,
								   pos, callback_func);

   graphics_info_t::mutate_sequence_chain_from_optionmenu = set_chain;
}


void mutate_sequence_chain_option_menu_item_activate (GtkWidget *item,
						      GtkPositionType pos) { 

   graphics_info_t::mutate_sequence_chain_from_optionmenu = menu_item_label(item);
}


// Don't use this function.  Use the one in graphics_info_t which you
// pass the callback function and get back a chain id.
// nvoid fill_chain_option_menu(GtkWidget *chain_option_menu, int imol) {

//   GtkSignalFunc callback_func =
//      GTK_SIGNAL_FUNC(mutate_sequence_chain_option_menu_item_activate);

   // fill_chain_option_menu_with_callback(chain_option_menu, imol, callback_func);
// }

// the generic form of the above
// void fill_chain_option_menu_with_callback(GtkWidget *chain_option_menu, int imol,
//  					  GtkSignalFunc callback_func) {

   // junk this function and use the one that returns a string.
// }




// The "Mutate" button action:
// 
void do_mutate_sequence(GtkWidget *dialog) {

#ifdef USE_PYTHON
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = coot::STATE_PYTHON;
#endif
#else // python not used
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = 0;
#endif
#endif   
   
   
   // decode the dialog here

   GtkWidget *entry1 = lookup_widget(dialog, "mutate_molecule_resno_1_entry");
   GtkWidget *entry2 = lookup_widget(dialog, "mutate_molecule_resno_2_entry");

   int t;
   int res1 = -9999, res2 = -99999;
   graphics_info_t g;
   
   const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry1));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res1 = t;
   entry_text = gtk_entry_get_text(GTK_ENTRY(entry2));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res2 = t;

// BL says: we should set a flag tha we swapped the direction and swap back
// before we call fit-gap to actually build backwards!!
   int swap_flag = 0;
   if (res2 < res1) {
      t = res1;
      res1 = res2;
      res2 = t;
      swap_flag = 1;
   }


   // set the imol and chain_id:
   // 
   int imol = graphics_info_t::mutate_sequence_imol;
   std::string chain_id = graphics_info_t::mutate_sequence_chain_from_optionmenu;

   // Auto fit?
   GtkWidget *checkbutton = lookup_widget(dialog, "mutate_sequence_do_autofit_checkbutton"); 
   short int autofit_flag = 0;

   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      autofit_flag = 1;
      

   if (imol>= 0) { // redundant
      if (is_valid_model_molecule(imol)) { 

	 // get the sequence:
	 GtkWidget *text = lookup_widget(dialog, "mutate_molecule_sequence_text");
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
	    sequence = coot::util::plain_text_to_sequence(sequence);
	    std::cout << "we got the sequence: " << sequence << std::endl;

	    if (int(sequence.length()) == (res2 - res1 + 1)) {
	       std::vector<std::string> cmd_strings;
	       if (autofit_flag)
		  cmd_strings.push_back("mutate-and-autofit-residue-range");
	       else 
		  cmd_strings.push_back("mutate-residue-range");
	       cmd_strings.push_back(graphics_info_t::int_to_string(imol));
	       cmd_strings.push_back(single_quote(chain_id));
	       cmd_strings.push_back(graphics_info_t::int_to_string(res1));
	       cmd_strings.push_back(graphics_info_t::int_to_string(res2));
	       cmd_strings.push_back(single_quote(sequence));
	       std::string cmd = g.state_command(cmd_strings, state_lang);
// BL says: I believe we should distinguish between python and guile here!?
#ifdef USE_GUILE
	       if (state_lang == coot::STATE_SCM) {
		  safe_scheme_command(cmd);
	       }
#else
#ifdef USE_PYTHON
              if (state_lang == coot::STATE_PYTHON) {
                 safe_python_command(cmd);
              }
#endif // PYTHON
#endif // GUILE
	       update_go_to_atom_window_on_changed_mol(imol);
	    } else {
	       std::cout << "WARNING:: can't mutate.  Sequence of length: "
			 << sequence.length() << " but residue range size: "
			 << res2 - res1 + 1 << "\n";
	    } 
	 } else {
	    std::cout << "WARNING:: can't mutate.  No sequence\n";
	 } 
      } else {
	 std::cout << "WARNING:: Bad molecule number: " << imol << std::endl;
	 std::cout << "          Can't mutate." << std::endl;
      }
   } else {
      std::cout << "WARNING:: unassigned molecule number: " << imol << std::endl;
      std::cout << "          Can't mutate." << std::endl;
   }
}

GtkWidget *wrapped_fit_loop_dialog() {

   GtkWidget *w = wrapped_create_mutate_sequence_dialog();

   GtkWidget *label              = lookup_widget(w, "function_for_molecule_label");
   GtkWidget *method_frame       = lookup_widget(w, "loop_fit_method_frame");
   GtkWidget *mutate_ok_button   = lookup_widget(w, "mutate_sequence_ok_button");
   GtkWidget *fit_loop_ok_button = lookup_widget(w, "fit_loop_ok_button");
   GtkWidget *checkbutton        = lookup_widget(w, "mutate_sequence_do_autofit_checkbutton");
#if (GTK_MAJOR_VERSION > 1)
   GtkWidget *rama_checkbutton   = lookup_widget(w, "mutate_sequence_use_ramachandran_restraints_checkbutton");
#endif
   
   gtk_label_set_text(GTK_LABEL(label), "\nFit loop in Molecule:\n");
   gtk_widget_hide(mutate_ok_button);
   gtk_widget_hide(checkbutton);
   gtk_widget_show(fit_loop_ok_button);
#if (GTK_MAJOR_VERSION > 1)
   gtk_widget_show(rama_checkbutton);
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rama_checkbutton), TRUE);
#endif

   gtk_widget_show(method_frame);

   return w;
}

// And the function called by the Fit Loop (OK) button.
// 
void fit_loop_from_widget(GtkWidget *dialog) {

#ifdef USE_PYTHON
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = coot::STATE_PYTHON;
#endif
#else // python not used
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = 0;
#endif
#endif

   // decode the dialog here

   GtkWidget *entry1 = lookup_widget(dialog, "mutate_molecule_resno_1_entry");
   GtkWidget *entry2 = lookup_widget(dialog, "mutate_molecule_resno_2_entry");

   int t;
   int res1 = -9999, res2 = -99999;
   graphics_info_t g;
   
   const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry1));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res1 = t;
   entry_text = gtk_entry_get_text(GTK_ENTRY(entry2));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res2 = t;

// BL says: we should set a flag that we swapped the direction and swap back
// before we call fit-gap to actually build backwards!!
   int swap_flag = 0;
   if (res2 < res1) {
      t = res1;
      res1 = res2;
      res2 = t;
      swap_flag = 1;
   }


   // set the imol and chain_id:
   // 
   int imol = graphics_info_t::mutate_sequence_imol;
   std::string chain_id = graphics_info_t::mutate_sequence_chain_from_optionmenu;

   // Auto fit?
   GtkWidget *checkbutton = lookup_widget(dialog, "mutate_sequence_do_autofit_checkbutton"); 
   short int autofit_flag = 0;

   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      autofit_flag = 1;

   // use Ramachandran restraints?
   int use_rama_restraints = 0;
#if (GTK_MAJOR_VERSION > 1)   
   GtkWidget *rama_checkbutton   = lookup_widget(dialog, "mutate_sequence_use_ramachandran_restraints_checkbutton");
   if (GTK_TOGGLE_BUTTON(rama_checkbutton)->active) 
      use_rama_restraints = 1;
#endif   

   if (imol>= 0) { // redundant
      if (is_valid_model_molecule(imol)) {

	 // get the sequence:
	 GtkWidget *text = lookup_widget(dialog, "mutate_molecule_sequence_text");
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
	    sequence = coot::util::plain_text_to_sequence(sequence);
	    int text_widget_sequence_length = sequence.length();
	    std::cout << "INFO:: mutating to the sequence :" << sequence
		      << ":" << std::endl;

	    if (int(sequence.length()) == (res2 - res1 + 1) ||
	        sequence == "") {
	    } else {
	       // so set sequence to poly-ala and give us a message:
	       sequence = "";
	       for (int i=0; i<(res2 - res1 + 1); i++)
		  sequence += "A";

	       std::cout << "WARNING:: Sequence of length: "
			 << text_widget_sequence_length << " but residue range size: "
			 << res2 - res1 + 1 << ".  Using Poly-Ala\n";
	       std::string s("WARNING:: Mis-matched sequence length\nUsing Poly Ala");
	       GtkWidget *w = wrapped_nothing_bad_dialog(s);
	       gtk_widget_show(w);
	    }
            if (swap_flag == 1) {
               t = res1;
               res1 = res2;
               res2 = t;
            }

	    std::vector<std::string> cmd_strings;
	    cmd_strings.push_back("fit-gap");
	    cmd_strings.push_back(graphics_info_t::int_to_string(imol));
	    cmd_strings.push_back(single_quote(chain_id));
	    cmd_strings.push_back(graphics_info_t::int_to_string(res1));
	    cmd_strings.push_back(graphics_info_t::int_to_string(res2));
	    cmd_strings.push_back(single_quote(sequence));
	    cmd_strings.push_back(graphics_info_t::int_to_string(use_rama_restraints));
	    std::string cmd = g.state_command(cmd_strings, state_lang);

#ifdef USE_GUILE
	    if (state_lang == coot::STATE_SCM) {
	       safe_scheme_command(cmd);
	    }
#else
#ifdef USE_PYTHON
            if (state_lang == coot::STATE_PYTHON) {
               safe_python_command(cmd);
            }
#endif // PYTHON
#endif // GUILE
	 }
      }
   }
}


/*  ----------------------------------------------------------------------- */
/*                         Align and Mutate GUI                             */
/*  ----------------------------------------------------------------------- */
GtkWidget *wrapped_create_align_and_mutate_dialog() {

   graphics_info_t g;
   GtkWidget *w = create_align_and_mutate_dialog();
   // fill w
   GtkWidget *mol_optionmenu   = lookup_widget(w, "align_and_mutate_molecule_optionmenu");
   GtkWidget *chain_optionmenu = lookup_widget(w, "align_and_mutate_chain_optionmenu");
   GtkSignalFunc callback = GTK_SIGNAL_FUNC(align_and_mutate_molecule_menu_item_activate);
   GtkSignalFunc chain_callback = GTK_SIGNAL_FUNC(align_and_mutate_chain_option_menu_item_activate);

   int imol = graphics_info_t::align_and_mutate_imol;
   if (imol == -1 || (! g.molecules[imol].has_model())) { 
      for (int i=0; i<g.n_molecules(); i++) {
	 if (g.molecules[i].has_model()) {
	    imol = i;
	    break;
	 }
      }
   }

   if (imol >= 0) {
      g.fill_option_menu_with_coordinates_options(mol_optionmenu, callback, imol);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_optionmenu, imol,
								      chain_callback);
      graphics_info_t::align_and_mutate_chain_from_optionmenu = set_chain;
   }
   
   return w;
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




/*  ----------------------------------------------------------------------- */
/*                  Change chain ID                                         */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_change_chain_id_dialog() {

   GtkWidget *w = create_change_chain_id_dialog();
   GtkWidget *mol_option_menu =  lookup_widget(w, "change_chain_id_molecule_optionmenu");
   GtkWidget *chain_option_menu =  lookup_widget(w, "change_chain_id_chain_optionmenu");
   GtkWidget *residue_range_no_radiobutton =  lookup_widget(w, "change_chain_residue_range_no_radiobutton");

   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(residue_range_no_radiobutton), TRUE);

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(change_chain_ids_mol_option_menu_item_activate);

   int imol = first_coords_imol();
   if (imol >= 0) {
      graphics_info_t::change_chain_id_molecule = imol;
      GtkSignalFunc chain_callback_func =
	 GTK_SIGNAL_FUNC(change_chain_ids_chain_menu_item_activate);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu,
								       imol,
								      chain_callback_func);
      graphics_info_t::change_chain_id_from_chain = set_chain;
   }
   graphics_info_t g; 
   g.fill_option_menu_with_coordinates_options(mol_option_menu, callback_func, imol);
   return w;
}

void
change_chain_ids_mol_option_menu_item_activate(GtkWidget *item,
					       GtkPositionType pos) {
   graphics_info_t::change_chain_id_molecule = pos;
   int imol = pos;
   GtkWidget *chain_option_menu =  lookup_widget(item, "change_chain_id_chain_optionmenu");
   GtkSignalFunc chain_callback_func =
      GTK_SIGNAL_FUNC(change_chain_ids_chain_menu_item_activate);
   std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu,
								   imol,
								   chain_callback_func);
   graphics_info_t::change_chain_id_from_chain = set_chain;
}

void
change_chain_ids_chain_menu_item_activate(GtkWidget *item,
					  GtkPositionType pos) {
   graphics_info_t::change_chain_id_from_chain = menu_item_label(item);
}


void
change_chain_id_by_widget(GtkWidget *w) {

   GtkWidget *residue_range_yes_radiobutton =  lookup_widget(w, "change_chain_residue_range_yes_radiobutton");

   GtkWidget *residue_range_from_entry  =  lookup_widget(w, "change_chain_residues_from_entry");
   GtkWidget *residue_range_to_entry  =  lookup_widget(w, "change_chains_residues_to_entry");
   GtkWidget *change_chains_new_chain_entry  =  lookup_widget(w, "change_chains_new_chain_id");

   int imol = graphics_info_t::change_chain_id_molecule;
   short int use_res_range_flag = 0;
   int from_resno = -9999;
   int to_resno = -9999;

   if (GTK_TOGGLE_BUTTON(residue_range_yes_radiobutton)->active) { 
      use_res_range_flag = 1;
      std::pair<short int, int> p1 = int_from_entry(residue_range_from_entry);
      std::pair<short int, int> p2 = int_from_entry(residue_range_to_entry);
      if (p1.first)
	 from_resno = p1.second;
      if (p2.first)
	 to_resno = p2.second;
   }

   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(change_chains_new_chain_entry));

   if (txt) { 
   
      if (is_valid_model_molecule(imol)) {
	 std::string to_chain_id(txt);
	 std::string from_chain_id = graphics_info_t::change_chain_id_from_chain;
	 std::pair<int, std::string> r = 
	    graphics_info_t::molecules[imol].change_chain_id(from_chain_id,
							     to_chain_id,
							     use_res_range_flag,
							     from_resno,
							     to_resno);
	 if (r.first == 1) { // it went OK
	    update_go_to_atom_window_on_changed_mol(imol);
	    graphics_draw();
	 } else {
	    GtkWidget *ws = wrapped_nothing_bad_dialog(r.second);
	    gtk_widget_show(ws);
	 }
      }
   } else {
      std::cout << "ERROR: Couldn't get txt in change_chain_id_by_widget\n";
   }
}


void fill_option_menu_with_refine_options(GtkWidget *option_menu) { 

   graphics_info_t g;

   g.fill_option_menu_with_map_options(option_menu, 
				       GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select));
}

void
set_rigid_body_fit_acceptable_fit_fraction(float f) {
   if (f >= 0.0 && f<= 1.0) { 
      graphics_info_t::rigid_body_fit_acceptable_fit_fraction = f;
   } else {
      std::cout << "ignoring set_rigid_body_fit_acceptable_fit_fraction"
		<< " of " << f << std::endl;
   } 
} 


void
my_delete_ramachandran_mol_option(GtkWidget *widget, void *data) {
   gtk_container_remove(GTK_CONTAINER(data), widget);
}


void
show_fix_nomenclature_errors_gui(int imol,
				 const std::vector<std::pair<std::string, coot::residue_spec_t> > &nomenclature_errors) {
   if (graphics_info_t::use_graphics_interface_flag) {
      if (is_valid_model_molecule(imol)) {

	 GtkWidget *w = create_fix_nomenclature_errors_dialog();

	 GtkWidget *label = lookup_widget(w, "fix_nomenclature_errors_label");

	 std::string s = "\n  Molecule number ";
	 s += coot::util::int_to_string(imol);
	 s += " has ";
	 s += coot::util::int_to_string(nomenclature_errors.size());
	 s += " nomenclature error";
	 if (nomenclature_errors.size() > 1)
	    s += "s";
	 s += ".\n";
	 if (nomenclature_errors.size() > 1)
	    s += "  Correct them?\n";
	 else 
	    s += "  Correct it?\n";

	 gtk_object_set_user_data(GTK_OBJECT(w), GINT_TO_POINTER(imol));
	 
	 gtk_label_set_text(GTK_LABEL(label), s.c_str());

	 GtkWidget *box = lookup_widget(w, "nomenclature_errors_vbox");

	 if (box) {
	    // fill box
	    for (unsigned int i=0; i<nomenclature_errors.size(); i++) {
	       s = nomenclature_errors[i].first; // the residue type
	       s += " ";
	       s += nomenclature_errors[i].second.format();
	       GtkWidget *label = gtk_label_new(s.c_str());
	       gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(label), FALSE, FALSE, 2);
	       gtk_widget_show(GTK_WIDGET(label));
	    }
	 }
	 gtk_widget_show(w);

      }
   }
}



/*  ----------------------------------------------------------------------- */
/*                  get molecule by libcheck/refmac code                    */
/*  ----------------------------------------------------------------------- */

/* Libcheck monomer code */
void 
handle_get_libcheck_monomer_code(GtkWidget *widget) { 

   const gchar *text = gtk_entry_get_text(GTK_ENTRY(widget));
   std::cout << "Refmac monomer Code: " << text << std::endl;
   get_monomer(text);

   // and kill the libcheck code window
   GtkWidget *window = lookup_widget(GTK_WIDGET(widget), "libcheck_monomer_dialog");
   if (window)
      gtk_widget_destroy(window);
   else 
      std::cout << "failed to lookup window in handle_get_libcheck_monomer_code" 
		<< std::endl;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                               skeleton                                   */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

GtkWidget *
create_skeleton_colour_selection_window() { 

   GtkWidget  *colorseldialog;
   GtkWidget *colorsel;

   colorseldialog = 
      gtk_color_selection_dialog_new("Skeleton Colour Selection"); 

/* How do we get to the buttons? */

   colorsel = GTK_COLOR_SELECTION_DIALOG(colorseldialog)->colorsel;

  /* Capture "color_changed" events in col_sel_window */

  gtk_signal_connect (GTK_OBJECT (colorsel), "color_changed",
                      (GtkSignalFunc)on_skeleton_color_changed, 
		      (gpointer)colorsel);
  
  gtk_signal_connect(GTK_OBJECT(GTK_COLOR_SELECTION_DIALOG(colorseldialog)->
				ok_button), "clicked",
		     GTK_SIGNAL_FUNC(on_skeleton_col_sel_cancel_button_clicked),
		     colorseldialog);

  gtk_signal_connect(GTK_OBJECT(GTK_COLOR_SELECTION_DIALOG(colorseldialog)->
				cancel_button), "clicked",
		     GTK_SIGNAL_FUNC(on_skeleton_col_sel_cancel_button_clicked), 
		     colorseldialog);

  gtk_color_selection_set_color(GTK_COLOR_SELECTION(colorsel),
				graphics_info_t::skeleton_colour);

  return GTK_WIDGET(colorseldialog);

}

/*! \brief show the strand placement gui.

  Choose the python version in there, if needed.  Call scripting
  function, display it in place, don't return a widget. */
void   place_strand_here_dialog() {

   if (graphics_info_t::use_graphics_interface_flag) {
      if (graphics_info_t::prefer_python) {
#ifdef USE_PYTHON
	 std::cout << "safe python commaond place_strand_here_gui()"
		   << std::endl;
	 safe_python_command("place_strand_here_gui()");
#endif // PYTHON
      } else {
#ifdef USE_GUILE
	 safe_scheme_command("(place-strand-here-gui)");
#endif 	 
      } 
   } 
}


/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                               fast secondary structure search            */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

GtkWidget *
wrapped_create_fast_ss_search_dialog() {

#if (GTK_MAJOR_VERSION > 1)
  GtkWidget *dialog;
  GtkWidget *helix_temp_combobox;
  GtkWidget *strand_temp_combobox;
  GtkWidget *helix_noaa_combobox;
  GtkWidget *strand_noaa_combobox;
  GtkWidget *radius_combobox;

  dialog = create_fast_ss_search_dialog();

  helix_temp_combobox = lookup_widget(dialog, "fast_sss_dialog_helix_template_combobox");
  helix_noaa_combobox = lookup_widget(dialog, "fast_sss_dialog_helix_no_aa_combobox");
  strand_temp_combobox = lookup_widget(dialog, "fast_sss_dialog_strand_template_combobox");
  strand_noaa_combobox = lookup_widget(dialog, "fast_sss_dialog_strand_no_aa_combobox");
  radius_combobox = lookup_widget(dialog, "fast_sss_dialog_radius_combobox");

  // fill the comboboxes (done automatically, set the active ones)
  gtk_combo_box_set_active(GTK_COMBO_BOX(helix_temp_combobox), 0);
  gtk_combo_box_set_active(GTK_COMBO_BOX(helix_noaa_combobox), 1);
  gtk_combo_box_set_active(GTK_COMBO_BOX(strand_temp_combobox), 1);
  gtk_combo_box_set_active(GTK_COMBO_BOX(strand_noaa_combobox), 0);
  gtk_combo_box_set_active(GTK_COMBO_BOX(radius_combobox),1);

  return dialog;
#else
  GtkWidget *w = 0;
  return w;
#endif // GTK_MAJOR_VERSION
}
