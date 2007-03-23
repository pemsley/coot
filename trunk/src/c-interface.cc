/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by Paul Emsley, The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

// $Id: c-interface.cc 1458 2007-01-26 20:20:18Z emsley $
// $LastChangedDate: 2007-01-26 20:20:18 +0000 (Fri, 26 Jan 2007) $
// $Rev: 1458 $
 

#include <stdlib.h>
#include <iostream>

#if !defined(WINDOWS_MINGW) && !defined(_MSC_VER)
#include <glob.h> // for globbing.
#endif

#ifdef USE_GUILE
#include <guile/gh.h>
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON


#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
#define GTK_ENABLE_BROKEN
#if defined _MSC_VER
#include <windows.h>
#endif
#endif

#define HAVE_CIF  // will become unnessary at some stage.

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

#include "clipper/ccp4/ccp4_map_io.h"
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.

#include <vector>
#include <string>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"
#include "coot-utils.hh"
#include "coot-map-utils.hh"
#include "coot-database.hh"

// #include "xmap-interface.h"
#include "graphics-info.h"

#include "atom-utils.h" // asc_to_graphics

#include "BuildCas.h"

#include "trackball.h" // adding exportable rotate interface


#include "c-interface.h"
#include "cc-interface.hh"
#include "ppmutil.h"

#include "positioned-widgets.h"

// moving column_label selection to c-interface from mtz bits.
#include "cmtz-interface.hh"
// #include "mtz-bits.h" stuff from here moved to cmtz-interface

char *coot_version() { 
   char *r = new char[40];
   strncpy(r, VERSION, 38);
   return r;
}

// Return 0 if not a valid name ( -> #f in scheme)
// e.g. /a/b/c.pdb "d/e/f.mtz FWT PHWT"
// 
const char *molecule_name(int imol) {

   const char *r = NULL;
   if (is_valid_map_molecule(imol)) {
      r = graphics_info_t::molecules[imol].name_.c_str();
      return r;
   }
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].name_.c_str();
   }
   std::string cmd = "molecule-name";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   
   return r;
} 

//  Display characteristics:
//
// 
void set_esoteric_depth_cue(int istate) {

   graphics_info_t::esoteric_depth_cue_flag = istate;
   std::string cmd = "set-esoteric-depth-cue";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

int  esoteric_depth_cue_state() {
   add_to_history_simple("esoteric-depth-cue-state");
   return graphics_info_t::esoteric_depth_cue_flag;
}

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


/*! \brief shall we start up the Gtk and the graphics window? 

   if passed the command line argument --no-graphics, coot will not
   start up gtk itself.

   An interface function for Ralf.
*/
short int use_graphics_interface_state() {

   add_to_history_simple("use-graphics-interface-state");
   return graphics_info_t::use_graphics_interface_flag; 

}

/*! \brief start Gtk (and graphics) 

   This function is useful if it was not started already (which can be
   achieved by using the command line argument --no-graphics).

   An interface for Ralf */
void start_graphics_interface() {
   add_to_history_simple("start-graphics-interface");
   gtk_main(); 
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

// Is this a repeat of something?  I don't know (doesn't look like it)
// 
void run_generic_script(const std::vector<std::string> &cmd_strings) {
   
   graphics_info_t g;

#ifdef USE_GUILE   
   std::string s = g.state_command(cmd_strings, coot::STATE_SCM);
   safe_scheme_command(s);
#endif

#ifdef USE_PYTHON
#ifndef USE_GUILE
   std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
   safe_python_command(s);
#endif    
#endif    

   std::string cmd = "run-generic-script";
   std::vector<coot::command_arg_t> args;
   for(unsigned int i=0; i<cmd_strings.size(); i++) 
      args.push_back(clipper::String(cmd_strings[i]));
   add_to_history_typed(cmd, args);
}

   

/*  ----------------------------------------------------------------------- */
/*                  Display lists                                           */
/*  ----------------------------------------------------------------------- */
void set_display_lists_for_maps(int istat) {

   graphics_info_t::display_lists_for_maps_flag = istat;
   for (int i=0; i<graphics_info_t::n_molecules; i++)
      if (graphics_info_t::molecules[i].has_map())
	 graphics_info_t::molecules[i].update_map();
   std::string cmd = "set-display-lists-for-maps";
   std::vector<coot::command_arg_t> args;
   args.push_back(istat);
   add_to_history_typed(cmd, args);
   graphics_draw();
}


//
// Return the molecule number of the molecule that we just filled.
// Return -1 if there was a failure.
// 
int handle_read_draw_molecule(const char *filename) {

   return handle_read_draw_molecule_with_recentre(filename,
						  graphics_info_t::recentre_on_read_pdb);
}


int handle_read_draw_molecule_with_recentre(const char *filename,
					   int recentre_on_read_pdb_flag) {

   //
   // cout << "handle_read_draw_molecule: handling " << filename << endl;
   
   graphics_info_t g;
   int imol = g.n_molecules;
   if (! filename)
      return -1;
   
   std::string f(filename);

   g.expand_molecule_space_maybe();

   // returns e.g. ".ins"

   int  istat = -1;
   std::string extention = coot::util::file_name_extension(filename);
   if (coot::util::extension_is_for_shelx_coords(extention)) {
      return read_shelx_ins_file(filename);

   } else { 
      // recentre and not a backup-restore
      // -1 is for failure strangely.
      istat = g.molecules[imol].handle_read_draw_molecule(f, recentre_on_read_pdb_flag, 0);
   }
   if (istat > -1) {
      std::cout << "Molecule " << g.n_molecules << " read successfully\n";
      g.n_molecules++;

      // if the go to atom widget exists, update its optionmenu to
      // reflect the existance of this new molecule.

      if (g.go_to_atom_window) {
	 g.set_go_to_atom_molecule(imol);
	 g.update_go_to_atom_window_on_new_mol();
	 // g.update_go_to_atom_window_on_changed_mol(imol);
      } else {
	 // The Go To Atom window is not displayed.
	 g.set_go_to_atom_molecule(imol);
      }
      
      // now force a draw of the molecule
      //
      graphics_draw();

      std::vector<std::string> command_strings;
      command_strings.push_back("handle-read-draw-molecule-with-recentre");
      command_strings.push_back(single_quote(coot::util::intelligent_debackslash(filename)));
      command_strings.push_back(graphics_info_t::int_to_string(recentre_on_read_pdb_flag));
      add_to_history(command_strings);
      std::string s("Successfully read coordinates file ");
      s += filename;
      s += ".  Molecule number ";
      s += coot::util::int_to_string(imol);
      s += " created.";
      g.statusbar_text(s);
      return imol;
   } else {
      std::string s("Failed to read coordinates file ");
      s += filename;
      g.statusbar_text(s);
      return -1;
   }

   std::string cmd = "handle-read-draw-molecule-with-recentre";
   std::vector<coot::command_arg_t> args;
   args.push_back(filename);
   args.push_back(recentre_on_read_pdb_flag);
   add_to_history_typed(cmd, args);
}


int read_pdb(const char *filename) {
   return handle_read_draw_molecule(filename); 
} 

void set_draw_zero_occ_markers(int status) { 
   
   graphics_info_t g;
   g.draw_zero_occ_spots_flag = status;
   std::string cmd = "set-draw-zero-occ-markers";
   std::vector<coot::command_arg_t> args;
   args.push_back(status);
   add_to_history_typed(cmd, args);
   graphics_draw();
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
#endif    
}

int first_coords_imol() {

   int imol = -1;
   for (int i=0; i<graphics_n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }
   add_to_history_simple("first-coords-imol");
   return imol;
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


void hardware_stereo_mode() {

   if (graphics_info_t::use_graphics_interface_flag) { 
      if (graphics_info_t::display_mode != coot::HARDWARE_STEREO_MODE) {
	 int previous_mode = graphics_info_t::display_mode;
	 graphics_info_t::display_mode = coot::HARDWARE_STEREO_MODE;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in hardware_stereo_mode!\n";
	 } else { 
	    short int try_hardware_stereo_flag = 1;
	    GtkWidget *glarea = gl_extras(vbox, try_hardware_stereo_flag);
	    if (glarea) { 
	       std::cout << "INFO:: switch to hardware_stereo_mode succeeded\n";
	       if (graphics_info_t::i_fn_token) { 
		  toggle_idle_function(); // turn it off;
	       }
	       gtk_widget_destroy(graphics_info_t::glarea);
	       graphics_info_t::glarea = glarea;
	       gtk_widget_show(glarea);
	       graphics_draw();
	    } else {
	       std::cout << "WARNING:: switch to hardware_stereo_mode failed\n";
	       graphics_info_t::display_mode = previous_mode;
	    }
	 }
      } else {
	 std::cout << "Already in hardware stereo mode" << std::endl;
      }
   }
   add_to_history_simple("hardware-stereo-mode");

}

void mono_mode() {

   if (graphics_info_t::use_graphics_interface_flag) { 
      if (graphics_info_t::display_mode != coot::MONO_MODE) { 
	 int previous_mode = graphics_info_t::display_mode;
	 graphics_info_t::display_mode = coot::MONO_MODE;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in mono mode!\n";
	 } else {
	    short int try_hardware_stereo_flag = 0;
	    GtkWidget *glarea = gl_extras(vbox, try_hardware_stereo_flag);
	    if (glarea) { 
	       std::cout << "INFO:: switch to mono_mode succeeded\n";
	       if (graphics_info_t::i_fn_token) { 
		  toggle_idle_function(); // turn it off;
	       }
	       gtk_widget_destroy(graphics_info_t::glarea);
	       if (graphics_info_t::glarea_2) { 
		  gtk_widget_destroy(graphics_info_t::glarea_2);
		  graphics_info_t::glarea_2 = NULL;
	       }
	       graphics_info_t::glarea = glarea;
	       gtk_widget_show(glarea);
	       graphics_draw();
	    } else {
	       graphics_info_t::display_mode = previous_mode;
	       std::cout << "WARNING:: switch to mono mode failed\n";
	    }
	 } 
      } else {
	 // std::cout << "Already in mono mode" << std::endl; // we know.
      }
   }
   add_to_history_simple("mono-mode");
}

/*! \brief turn on side bye side stereo mode */
void side_by_side_stereo_mode(short int use_wall_eye_flag) {

   if (graphics_info_t::use_graphics_interface_flag) {
      if (graphics_info_t::display_mode != coot::SIDE_BY_SIDE_STEREO) {
	 if (use_wall_eye_flag == 1)
	    graphics_info_t::in_wall_eyed_side_by_side_stereo_mode = 1;
	 int previous_mode = graphics_info_t::display_mode;
	 short int stereo_mode = coot::SIDE_BY_SIDE_STEREO;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 GtkWidget *glarea = gl_extras(vbox, stereo_mode);
	 if (glarea) {
	    if (graphics_info_t::i_fn_token) { 
	       toggle_idle_function(); // turn it off;
	    }
	    gtk_widget_destroy(graphics_info_t::glarea);
	    graphics_info_t::glarea = glarea; // glarea_2 is stored by gl_extras()
	    gtk_widget_show(glarea);
	    gtk_widget_show(graphics_info_t::glarea_2);
	    graphics_draw();
	 } else {
	    std::cout << "WARNING:: switch to side by side mode failed!\n";
	 } 
      }
   }
   // add_to_history_simple("side-by-side-stereo-mode");
} 

/* DTI stereo mode - undocumented, secret interface for testing, currently */
// when it works, call it dti_side_by_side_stereo_mode()
void set_dti_stereo_mode(short int state) {

   if (graphics_info_t::use_graphics_interface_flag) {
      if (graphics_info_t::display_mode != coot::DTI_SIDE_BY_SIDE_STEREO) {
	 int previous_mode = graphics_info_t::display_mode;
	 short int stereo_mode = coot::DTI_SIDE_BY_SIDE_STEREO;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 GtkWidget *glarea = gl_extras(vbox, stereo_mode);
	 if (glarea) {
	    if (graphics_info_t::i_fn_token) { 
	       toggle_idle_function(); // turn it off;
	    }
	    gtk_widget_destroy(graphics_info_t::glarea);
	    graphics_info_t::glarea = glarea; // glarea_2 is stored by gl_extras()
	    gtk_widget_show(glarea);
	    gtk_widget_show(graphics_info_t::glarea_2);
	    graphics_draw();
	 } else {
	    std::cout << "WARNING:: switch to side by side mode failed!\n";
	 } 
      }
   }
   // add_to_history_simple("dti-side-by-side-stereo-mode");
} 


int stereo_mode_state() {
   add_to_history_simple("stereo-mode-state");
   return graphics_info_t::display_mode;
}

void set_hardware_stereo_angle_factor(float f) {
   graphics_info_t::hardware_stereo_angle_factor = f;
   std::string cmd = "set-hardware-stereo-angel-factor";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

float hardware_stereo_angle_factor_state() {
   add_to_history_simple("hardware-stereo-angle-factor-state");
   return graphics_info_t::hardware_stereo_angle_factor;
}


void graphics_draw() {
   graphics_info_t::graphics_draw();
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
   Py_Finalize();
#endif
   gtk_exit(retval); 

}

void 
set_run_state_file_status(short int istat) {
   std::string cmd = "set-run-state-file-status";
   std::vector<coot::command_arg_t> args;
   args.push_back(istat);
   add_to_history_typed(cmd, args);
   
   graphics_info_t::run_state_file_status = istat;
}


void set_sticky_sort_by_date() {

   add_to_history_simple("set-sticky-sort-by-date");
   graphics_info_t g;
   g.sticky_sort_by_date = 1;

}

void add_filename_filter(GtkWidget *fileselection) { 

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

// where data type:
// 0 coords
// 1 mtz etc
// 2 maps
// 
GtkWidget *add_filename_filter_button(GtkWidget *fileselection, 
				      short int data_type) { 
   
   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *button = gtk_toggle_button_new_with_label("Filter");
   GtkWidget *frame = gtk_frame_new("File-name filter:");
   int d = data_type;
   
   gtk_widget_ref(button);
   gtk_widget_show(button);
   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), button);
   gtk_signal_connect (GTK_OBJECT (button), "toggled",
		       GTK_SIGNAL_FUNC (on_filename_filter_toggle_button_toggled),
		       GINT_TO_POINTER(d));
   gtk_widget_show(frame);
   return button;
}

void add_is_difference_map_checkbutton(GtkWidget *fileselection) { 

   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *button = gtk_check_button_new_with_label("Is Difference Map");
   GtkWidget *frame = gtk_frame_new("Difference Map?");
   int imol = 0; // FIXME

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
		       GINT_TO_POINTER(imol));
   gtk_widget_show(frame);

}

void set_filter_fileselection_filenames(int istate) {
   std::string cmd = "set-filter-fileselection-filenames";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);

   graphics_info_t::filter_fileselection_filenames_flag = istate;

}

int filter_fileselection_filenames_state() {
   add_to_history_simple("filter_fileselection_filenames_state");
   return graphics_info_t::filter_fileselection_filenames_flag;
}


void swap_map_colours(int imol1, int imol2) {

   if (is_valid_map_molecule(imol1)) {
      if (is_valid_map_molecule(imol2)) {
	 graphics_info_t g;
	 std::vector<float> map_1_colours = g.molecules[imol1].map_colours();
	 std::vector<float> map_2_colours = g.molecules[imol2].map_colours();
	 double *colours1 = new double[4];
	 colours1[0] = map_1_colours[0];
	 colours1[1] = map_1_colours[1];
	 colours1[2] = map_1_colours[2];
	 double *colours2 = new double[4];
	 colours2[0] = map_2_colours[0];
	 colours2[1] = map_2_colours[1];
	 colours2[2] = map_2_colours[2];
	 g.molecules[imol1].handle_map_colour_change(colours2, g.swap_difference_map_colours);
	 g.molecules[imol2].handle_map_colour_change(colours1, g.swap_difference_map_colours);
	 delete [] colours1;
	 delete [] colours2;
      }
   }
   std::string cmd = "swap-map-colours";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol1);
   args.push_back(imol2);
   add_to_history_typed(cmd, args);
}

void set_keep_map_colour_after_refmac(int istate) {
   std::string cmd = "set-keep-map-colour-after-refmac";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
   graphics_info_t::swap_pre_post_refmac_map_colours_flag = istate;
}

int keep_map_colour_after_refmac_state() {
   add_to_history_simple("keep_map_colour_after_refmac_state");
   return graphics_info_t::swap_pre_post_refmac_map_colours_flag;
} 

void
on_read_map_difference_map_toggle_button_toggled (GtkButton       *button,
						  gpointer         user_data)
{
   if (GTK_TOGGLE_BUTTON(button)->active) { 
      std::cout << "is a difference map...!\n";
   } 
}

void
on_filename_filter_toggle_button_toggled (GtkButton       *button,
					  gpointer         user_data)
{
   
   int data_type = GPOINTER_TO_INT(user_data);

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
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button);
   
   std::vector<std::string> v;
   
   if (fileselection) { 
      if (GTK_TOGGLE_BUTTON(button)->active) { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Unfilter");
	 
	 // so now we have pre_directory
	 // 
	 // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
	 std::vector<std::string> v = filtered_by_glob(pre_directory, data_type);
	 // we want to stat the directory and add all the files in it:

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
      } else { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Filter");
	 gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
					 (pre_directory + "/").c_str());
      }
   } else {
      std::cout << "no fileselection found from sort button\n";
   }
}

std::vector<std::string> filtered_by_glob(const std::string &pre_directory, 
					  int data_type) { 

   std::vector<std::string> v; // returned object
   std::vector<std::string> globs;

#if !defined(WINDOWS_MINGW) && !defined(_MSC_VER)

   // std::map<std::string, int, std::less<std::string> >  files;

   if (data_type == 0) 
      globs = *graphics_info_t::coordinates_glob_extensions;
   if (data_type == 1) 
      globs = *graphics_info_t::data_glob_extensions;
   if (data_type == 2) 
      globs = *graphics_info_t::map_glob_extensions;
   if (data_type == 3) 
      globs = *graphics_info_t::dictionary_glob_extensions;

   for (unsigned int i=0; i<globs.size(); i++) { 

      std::string file_name_glob = pre_directory;
      file_name_glob += "/";

      file_name_glob += "*";
      file_name_glob += globs[i];
      glob_t myglob;
      int flags = 0;
      glob(file_name_glob.c_str(), flags, 0, &myglob);
      size_t count;

      char **p;
      for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) {
	 std::string f(*p);
	 if (! string_member(f, v))
	    v.push_back(f);
      }
      globfree(&myglob);
   }

   // now we need to sort v;
   std::sort(v.begin(), v.end(), compare_strings);

#endif // WINDOWS_MINGW

   return v;
}

bool
compare_strings(const std::string &a, const std::string &b) { 
   return a < b ? 1 : 0;
} 


// Return 1 if search appears in list, 0 if not)
// 
short int 
string_member(const std::string &search, const std::vector<std::string> &list) { 
   
   short int v = 0;
   for (unsigned int i=0; i<list.size(); i++) { 
      if (list[i] == search) { 
	 v = 1;
	 break;
      }
   }
   return v;
} 

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
      handle_filename_filter(widget);
   } 
   return FALSE;
}

void
handle_filename_filter(GtkWidget *entry_widget) { 

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


   // so now we have pre_directory
   // 
   // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button);
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

  if (filter_fileselection_filenames_state()) { 
    gtk_signal_emit_by_name(GTK_OBJECT(filter_button), "clicked");
  }
  if (graphics_info_t::sticky_sort_by_date) {
     GtkWidget *file_list = GTK_FILE_SELECTION(fileselection)->file_list;
     std::cout << "INFO:: Sorting files by date\n";
     fileselection_sort_button_clicked(sort_button, (GtkCList *) file_list);
  }
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

// --------------------------------------------------------------------
// Ctrl for rotate or pick:
// --------------------------------------------------------------------
// 
// Coot mailing list discussion: users want Ctrl for Rotation or Ctrl
// for picking, so that accidental picking when rotation is meant is
// avoided.
// 
void set_control_key_for_rotate(int state) {
   graphics_info_t::control_key_for_rotate_flag = state;
}

int control_key_for_rotate_state() {
   return graphics_info_t::control_key_for_rotate_flag;
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



/*  ------------------------------------------------------------------------ */
/*                         Model/Fit/Refine Functions:                       */
/*  ------------------------------------------------------------------------ */
void post_model_fit_refine_dialog() {

   GtkWidget *widget = wrapped_create_model_fit_refine_dialog();
   if (graphics_info_t::use_graphics_interface_flag) { 
      gtk_widget_show(widget);
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("post-model-fit-refine-dialog");
   add_to_history(command_strings);
}

void post_other_modelling_tools_dialog() {

   GtkWidget *widget = wrapped_create_model_fit_refine_dialog();
   if (graphics_info_t::use_graphics_interface_flag) { 
      gtk_widget_show(widget);
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("post-other-modelling-tools-dialog");
   add_to_history(command_strings);

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



void show_select_map_dialog() {
   graphics_info_t g;
   g.show_select_map_dialog();
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


// 
int make_and_draw_map(const char* mtz_file_name,
		      const char *f_col, const char *phi_col,
		      const char *weight_col, int use_weights,
		      int is_diff_map) {
    
   graphics_info_t g;
   int imol = -1; // failure initially.
   struct stat buf;

   std::string f_col_str(f_col);
   std::string phi_col_str(phi_col);
   std::string weight_col_str("");
   if (use_weights)
      weight_col_str = std::string(weight_col);

   int status = stat(mtz_file_name, &buf);
   // 
   if (status != 0) {
      std::cout << "WARNING:: Can't find file " << mtz_file_name << std::endl;
      if (S_ISDIR(buf.st_mode)) {
	 std::cout << mtz_file_name << " is a directory! " << std::endl;
      }
   } else {

      if (valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)) { 
      
	 std::cout << "making map from mtz filename " << mtz_file_name << std::endl;
	 imol = g.n_molecules;
	 g.molecules[imol].map_fill_from_mtz(std::string(mtz_file_name),
					     f_col_str,
					     phi_col_str,
					     weight_col_str,
					     use_weights, is_diff_map);
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 graphics_draw();
	 g.activate_scroll_radio_button_in_display_manager(imol);
	 std::vector<std::string> command_strings;
	 command_strings.push_back("make-and-draw-map");
	 command_strings.push_back(single_quote(mtz_file_name));
	 command_strings.push_back(single_quote(f_col));
	 command_strings.push_back(single_quote(phi_col));
	 command_strings.push_back(single_quote(weight_col));
	 command_strings.push_back(graphics_info_t::int_to_string(use_weights));
	 command_strings.push_back(graphics_info_t::int_to_string(is_diff_map));
	 add_to_history(command_strings);
	 
      } else {
	 std::cout << "WARNING:: label(s) not found in mtz file " 
		   << mtz_file_name << " " << f_col_str << " " 
		   <<  phi_col_str << " ";
	 if (use_weights)
	    std::cout << weight_col_str << std::endl;
	 else 
	    std::cout << std::endl;
      }      
   }
   return imol; // possibly -1
}

int  make_and_draw_map_with_refmac_params(const char *mtz_file_name, 
					  const char *a, const char *b,
					  const char *weight,
					  int use_weights, int is_diff_map,
					  short int have_refmac_params,
					  const char *fobs_col,
					  const char *sigfobs_col,
					  const char *r_free_col,
					  short int sensible_f_free_col) { 

   graphics_info_t g;
   int imol = g.n_molecules;

   // this is order dependent.  the restore-state comand that is
   // constructed in make_and_draw_map checks to see if we
   // have_refmac_params, so we need to set them before making the map.
   // 
   if (have_refmac_params)
      g.molecules[imol].store_refmac_params(std::string(mtz_file_name),
					    std::string(fobs_col), 
					    std::string(sigfobs_col), 
					    std::string(r_free_col),
					    sensible_f_free_col);
   
   return make_and_draw_map(mtz_file_name, a, b, weight, use_weights, is_diff_map);
}

// return imol, possibly -1;
int make_and_draw_map_with_reso_with_refmac_params(const char *mtz_file_name, 
						   const char *f_col,
						   const char *phi_col,
						   const char *weight_col,
						   int use_weights, int is_diff_map,
						   short int have_refmac_params,
						   const char *fobs_col,
						   const char *sigfobs_col,
						   const char *r_free_col,
						   short int sensible_r_free_col,
						   short int is_anomalous_flag,
						   short int use_reso_limits,
						   float low_reso_limit,
						   float high_reso_limit) {

   graphics_info_t g;
   int imol = -1;

   // this is order dependent.  the restore-state comand that is
   // constructed in make_and_draw_map checks to see if we
   // have_refmac_params, so we need to set them before making the map.
   // 
   if (have_refmac_params)
      g.molecules[g.n_molecules].store_refmac_params(std::string(mtz_file_name),
					    std::string(fobs_col), 
					    std::string(sigfobs_col), 
					    std::string(r_free_col),
					    sensible_r_free_col);
   
   struct stat buf;
   int status = stat(mtz_file_name, &buf);

   if (status != 0) {
      std::cout << "Error finding MTZ file " << mtz_file_name << std::endl;
      if (S_ISDIR(buf.st_mode)) {
	 std::cout << mtz_file_name << " is a directory! " << std::endl;
      }
   } else {
      std::string map_type;
      if (is_diff_map)
	 map_type = "difference";
      else
	 map_type = "conventional";
      
      std::cout << "making " << map_type << " map from MTZ filename "
		<< mtz_file_name << " using " << f_col << " "
		<< phi_col << std::endl;

      if (valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)) {
	 std::string weight_col_str("");
	 if (use_weights)
	    weight_col_str = std::string(weight_col);
	 imol = g.n_molecules;
	 g.molecules[imol].map_fill_from_mtz_with_reso_limits(std::string(mtz_file_name),
							      std::string(f_col),
							      std::string(phi_col),
							      weight_col_str,
							      use_weights,
							      is_anomalous_flag,
							      is_diff_map,
							      use_reso_limits,
							      low_reso_limit,
							      high_reso_limit);
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 g.activate_scroll_radio_button_in_display_manager(imol);
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: label(s) not found in MTZ file " 
		   << mtz_file_name << " " << f_col << " " 
		   <<  phi_col << " ";
	 if (use_weights)
	    std::cout << weight_col << std::endl;
	 else 
	    std::cout << std::endl;
      }
   }
   if (imol != -1) {

      // We reset some strings if we weren't given refmac params -
      // otherwise we quote garbage or unallocated memory.
      std::string weight_col_str;
      std::string fobs_col_str;
      std::string sigfobs_col_str;
      std::string r_free_col_str;
      if (weight_col)
	 weight_col_str = single_quote(weight_col);
      else
	 weight_col_str = single_quote("Weight:None-specified");
      
      if (! have_refmac_params) {
	 fobs_col_str    = single_quote("Fobs:None-specified");
	 sigfobs_col_str = single_quote("SigF:None-specified");
	 r_free_col_str  = single_quote("RFree:None-specified");
	 sensible_r_free_col = 0;
      } else {
	 fobs_col_str    = single_quote(fobs_col);
	 sigfobs_col_str = single_quote(sigfobs_col_str);
	 r_free_col_str  = single_quote(r_free_col);
      }
      
      std::vector<std::string> command_strings;
      command_strings.push_back("make-and-draw-map-with-reso-with-refmac-params");
      command_strings.push_back(single_quote(mtz_file_name));
      command_strings.push_back(single_quote(f_col));
      command_strings.push_back(single_quote(phi_col));
      command_strings.push_back(weight_col_str);
      command_strings.push_back(graphics_info_t::int_to_string(use_weights));
      command_strings.push_back(graphics_info_t::int_to_string(is_diff_map));
      command_strings.push_back(graphics_info_t::int_to_string(have_refmac_params));
      command_strings.push_back(fobs_col_str);
      command_strings.push_back(sigfobs_col_str);
      command_strings.push_back(r_free_col_str);
      command_strings.push_back(graphics_info_t::int_to_string(sensible_r_free_col));
      command_strings.push_back(graphics_info_t::int_to_string(is_anomalous_flag));
      command_strings.push_back(graphics_info_t::int_to_string(use_reso_limits));
      command_strings.push_back(graphics_info_t::float_to_string(low_reso_limit));
      command_strings.push_back(graphics_info_t::float_to_string(high_reso_limit));
      add_to_history(command_strings);
   }
   return imol;
}


int auto_read_make_and_draw_maps(const char *mtz_file_name) {

   const char *f_col   = graphics_info_t::auto_read_MTZ_FWT_col.c_str();
   const char *phi_col = graphics_info_t::auto_read_MTZ_PHWT_col.c_str();
   const char *weight_col = "";
   
   const char *f_col_diff   = graphics_info_t::auto_read_MTZ_DELFWT_col.c_str();
   const char *phi_col_diff = graphics_info_t::auto_read_MTZ_PHDELWT_col.c_str();
   int imol1, imol2;
   
   imol1 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							  f_col,
							  phi_col,
							  weight_col,
							  0, 0, //    use_weights,  is_diff_map,
							  0,  //   short int have_refmac_params,
							  "", //   const char *fobs_col,
							  "", //   const char *sigfobs_col,
							  "", //   const char *r_free_col,
							  0,  //   short int sensible_f_free_col,
							  0,  //   short int is_anomalous_flag,
							  0,  //   short int use_reso_limits,
							  0,  //   float low_reso_limit,
							  0); //   float high_reso_limit
   
   imol2 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							  f_col_diff,
							  phi_col_diff,
							  weight_col,
							  0, 1, //    use_weights,  is_diff_map,
							  0,  //   short int have_refmac_params,
							  "", //   const char *fobs_col,
							  "", //   const char *sigfobs_col,
							  "", //   const char *r_free_col,
							  0,  //   short int sensible_f_free_col,
							  0,  //   short int is_anomalous_flag,
							  0,  //   short int use_reso_limits,
							  0,  //   float low_reso_limit,
							  0); //   float high_reso_limit
   
   std::string s;
   if (imol1 < 0) { 
      s += "Failed to find columns ";
      s += graphics_info_t::auto_read_MTZ_FWT_col;
      s += " and ";
      s += graphics_info_t::auto_read_MTZ_PHWT_col;
      s += "in that mtz file\n";
   }
   if (imol2 < 0) { 
      s += "Failed to find columns ";
      s += f_col_diff;
      s += " and ";
      s += phi_col_diff;
      s += " in that mtz file\n";
   }
   if (imol1<0 || imol2 <0) {

      imol1 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							     "2FOFCWT",
							     "PH2FOFCWT",
							     weight_col,
							     0, 0, //    use_weights,  is_diff_map,
							     0,  //   short int have_refmac_params,
							     "", //   const char *fobs_col,
							     "", //   const char *sigfobs_col,
							     "", //   const char *r_free_col,
							     0,  //   short int sensible_f_free_col,
							     0,  //   short int is_anomalous_flag,
							     0,  //   short int use_reso_limits,
							     0,  //   float low_reso_limit,
							     0); //   float high_reso_limit
      
      imol2 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							     "FOFCWT",
							     "PHFOFCWT",
							     weight_col,
							     0, 1, //    use_weights,  is_diff_map,
							     0,  //   short int have_refmac_params,
							     "", //   const char *fobs_col,
							     "", //   const char *sigfobs_col,
							     "", //   const char *r_free_col,
							     0,  //   short int sensible_f_free_col,
							     0,  //   short int is_anomalous_flag,
							     0,  //   short int use_reso_limits,
							     0,  //   float low_reso_limit,
							     0); //   float high_reso_limit
   

      if (imol1<0 || imol2 <0) {
	 GtkWidget *w = wrapped_nothing_bad_dialog(s);
	 gtk_widget_show(w);
      }
   }
   
   // we like the 2fo-fc map to be the scrolled one after autoreading,
   if (imol1 >= 0) { 
      graphics_info_t g;
      g.scroll_wheel_map = imol1;
      g.activate_scroll_radio_button_in_display_manager(imol1);
   }
   
   return imol2; // return to callbacks.c (currently ignored)
} 

int auto_read_do_difference_map_too_state() {

   int i = graphics_info_t::auto_read_do_difference_map_too_flag; 

   return i;

} 
void set_auto_read_do_difference_map_too(int i) {

   graphics_info_t::auto_read_do_difference_map_too_flag = i;
} 


void set_auto_read_column_labels(const char *fwt, const char *phwt, 
				 int is_for_diff_map_flag) {

   if (is_for_diff_map_flag) {
      graphics_info_t::auto_read_MTZ_DELFWT_col = fwt;
      graphics_info_t::auto_read_MTZ_PHDELWT_col = phwt;
   } else {
      graphics_info_t::auto_read_MTZ_FWT_col = fwt;
      graphics_info_t::auto_read_MTZ_PHWT_col = phwt;
   }
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


void toggle_idle_function() { 

   graphics_info_t g; 

   if (g.i_fn_token == 0) { 
      g.i_fn_token = gtk_idle_add((GtkFunction)animate_idle, g.glarea);
   } else {
      gtk_idle_remove(g.i_fn_token);
      g.i_fn_token = 0; 
   }
   add_to_history_simple("toggle_idle_function");
}

// in degrees
void set_idle_function_rotate_angle(float f) {

   graphics_info_t g;
   std::string cmd = "set-idle-function-rotate-angle";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
   g.idle_function_rotate_angle = f; 
}

float idle_function_rotate_angle() {

   std::string cmd = "idle-function-rotate-angle";
   add_to_history_simple(cmd);
   return graphics_info_t::idle_function_rotate_angle;
} 



void do_tw() {

}

// Another name for wrapped_nothing_bad_dialog, but this function also
// displays the widget so nothing is returned.
void info_dialog(const char *txt) {

   if (graphics_info_t::use_graphics_interface_flag) { 
      if (txt) {
	 std::string s(txt);
	 GtkWidget *w = wrapped_nothing_bad_dialog(s);
	 gtk_widget_show(w);
      }
   }
   std::string cmd = "info-dialog";
   std::vector<coot::command_arg_t> args;
   args.push_back(txt);
   add_to_history_typed(cmd, args);
}


GtkWidget *main_menubar() {

   GtkWidget *w = lookup_widget(graphics_info_t::statusbar, "menubar1");
   return w; 
}

GtkWidget *main_statusbar() {
   return graphics_info_t::statusbar;
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


/*  ------------------------------------------------------------------------ */
/*                         file selection                                    */
/*  ------------------------------------------------------------------------ */

void
set_directory_for_fileselection(GtkWidget *fileselection1) { 
   graphics_info_t g;
   g.set_directory_for_fileselection(fileselection1);
}

void
save_directory_from_fileselection(const GtkWidget *fileselection) {
   graphics_info_t g;
   g.save_directory_from_fileselection(fileselection);
}


/*  Eleanor likes to sort her files by date when selecting a file
*/
GtkWidget *add_sort_button_fileselection(GtkWidget *fileselection) {

   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *frame = gtk_frame_new("File Order");
   GtkWidget *button = gtk_button_new_with_label("  Sort by Date  ");
   gtk_widget_ref(button);
   gtk_object_set_data_full(GTK_OBJECT(fileselection),
			    "fileselection_sort_button",
			    button,
			    (GtkDestroyNotify) gtk_widget_unref);
   GtkWidget *file_list = GTK_FILE_SELECTION(fileselection)->file_list;
 
   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(GTK_FILE_SELECTION(fileselection)->history_pulldown);

   gtk_object_set_user_data(GTK_OBJECT(button), (char *) history_pulldown); 

   gtk_signal_connect (GTK_OBJECT(button), "clicked",
		       (GtkSignalFunc) fileselection_sort_button_clicked,
		       file_list);

   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), button);
   gtk_widget_show(frame);
   gtk_widget_show(button);
   return button;
}

short int compare_mtimes(str_mtime a, str_mtime b) {
   return a.mtime > b.mtime;
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




void fileselection_sort_button_clicked( GtkWidget *sort_button,
					GtkCList  *file_list) {

   std::vector<str_mtime> v;
   char *text;
   struct stat buf;
   time_t mtime;

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
	 std::cout << "null label t " << std::endl;
      } 
      dlist = dlist->next;
   }
   g_list_free(free_list); 

   int status;
   std::string file_name;
   for(int i=0; i<file_list->rows; i++) {
      gtk_clist_get_text(file_list,i,0,&text);
      file_name = pre_directory;
      if (file_name != "") {
	 file_name += "/";
	 file_name += text;
      } else {
	 // Oh dear.  This directory only.  Else failure
	 // this shouldn't happen hopefully.
	 file_name = text;
      } 
      status = stat(file_name.c_str(),&buf);
      if (status == 0) { 
	 mtime = buf.st_mtime;
	 v.push_back(str_mtime(text,mtime));
      } else {
	 std::cout << "error stating " << file_name << std::endl;
      }
   }
//    v.push_back(str_mtime("Cabbage",324232));
//    v.push_back(str_mtime("zipped",345234523452));
   std::sort(v.begin(), v.end(), compare_mtimes);

   gtk_clist_clear(file_list);

   for(unsigned int i=0; i<v.size(); i++) {
      text = (char *) v[i].file.c_str(); 
      gtk_clist_append(file_list, &text);
   }
}

void quanta_buttons() {
   graphics_info_t g;
   g.quanta_buttons();
}

void quanta_like_zoom() { 
   graphics_info_t::quanta_like_zoom_flag = 1;
} 

void set_scroll_by_wheel_mouse(int istate) {
   graphics_info_t::do_scroll_by_wheel_mouse_flag = istate;
}

int scroll_by_wheel_mouse_state() {
   return graphics_info_t::do_scroll_by_wheel_mouse_flag;
}

/*! \brief set the default inital contour for 2FoFc-style map

in sigma */
void set_default_initial_contour_level_for_map(float n_sigma) {

   graphics_info_t::default_sigma_level_for_map = n_sigma;

} 

/*! \brief set the default inital contour for FoFc-style map

in sigma */
void set_default_initial_contour_level_for_difference_map(float n_sigma) {

   graphics_info_t::default_sigma_level_for_fofc_map = n_sigma;

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

void set_map_line_width(int w) {
   graphics_info_t::map_line_width = w;
   // update the maps because they may be being draw as graphical
   // objects.
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++)
      graphics_info_t::molecules[imol].update_map();
   graphics_draw();
}

int map_line_width_state() {
   return graphics_info_t::map_line_width;
} 


//
char* get_text_for_symmetry_size_widget() {

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
   snprintf(text,100,"%-5.1f", g.symmetry_search_radius);

   return text;

}

void set_density_size_from_widget(const char *text) {

   float tmp;
   graphics_info_t g;

   tmp = atof(text);

   if ((tmp > 0.0) && (tmp < 9999.9)) {
      g.box_radius = tmp;
   } else {

      cout << "Cannot interpret " << text << ".  Assuming 10A" << endl;
      g.box_radius = 10.0;
   }
   //
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_map();
   }
   graphics_draw();
}

void set_density_size(float f) {

   graphics_info_t g;
   g.box_radius = f;
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_map();
   }
   graphics_draw();
}

/*! \brief set the extent of the box/radius of electron density contours */
void set_map_radius(float f) {

   set_density_size(f);

} 



void set_display_intro_string(const char *str) {

   if (graphics_info_t::use_graphics_interface_flag) { 
      if (str) { 
	 std::string s(str);
	 graphics_info_t g;
	 g.display_density_level_screen_string = s;
	 g.statusbar_text(s);
      }
   }
} 

void set_swap_difference_map_colours(int i) { 
   graphics_info_t::swap_difference_map_colours = i;
}

/* return success status 0 = failure (imol does not have a map) */
int set_map_is_difference_map(int imol) { 

   int istatus = 0;
   if (imol< graphics_n_molecules()) { 
      if (graphics_info_t::molecules[imol].has_map()) { 
	 graphics_info_t::molecules[imol].set_map_is_difference_map();
	 istatus = 1;
	 graphics_draw();
      } else { 
	 std::cout << "WARNING:: molecule " << imol << " does not have a map." <<  std::endl;
      } 

   } else { 
      std::cout << "WARNING:: No such molecule as " << imol << std::endl;
   } 
   return istatus;

} 

/* return the index of the new molecule or -1 on failure */
int another_level() {

   int istat = -1;
   int imap = -1;
   for (int i=0; i<graphics_info_t::n_molecules; i++) {
      if (is_valid_map_molecule(i)) {
	 if (! graphics_info_t::molecules[i].is_difference_map_p()) { 
	    imap = i;
	 }
      }
   }
   if (imap > -1) {
      istat = another_level_from_map_molecule_number(imap);
   }

   return istat;
}

int another_level_from_map_molecule_number(int imap) {
   int istat = -1;
   if (is_valid_map_molecule(imap)) {
      // create another map with the same parameters as imap and then
      // push up the contour level a sigma.
//       std::cout << "DEBUG:: calling make_and_draw_map_with_reso_with_refmac_params"
// 		<< std::endl;
      istat = make_and_draw_map_with_reso_with_refmac_params(
	  graphics_info_t::molecules[imap].save_mtz_file_name.c_str(),
	  graphics_info_t::molecules[imap].save_f_col.c_str(),
	  graphics_info_t::molecules[imap].save_phi_col.c_str(),
          graphics_info_t::molecules[imap].save_weight_col.c_str(),
          graphics_info_t::molecules[imap].save_use_weights,
          graphics_info_t::molecules[imap].save_is_diff_map_flag,
	  0, "None", "None", "None", 0, // refmac params
          graphics_info_t::molecules[imap].save_is_anomalous_map_flag,
          graphics_info_t::molecules[imap].save_use_reso_limits,
          graphics_info_t::molecules[imap].save_low_reso_limit,
          graphics_info_t::molecules[imap].save_high_reso_limit);

      if (istat != -1) { 

	 float map_sigma = graphics_info_t::molecules[istat].map_sigma();
	 float current_contour_level = graphics_info_t::molecules[istat].contour_level[0];
	 graphics_info_t::molecules[istat].set_contour_level(current_contour_level +
							     map_sigma*1.0);
	 graphics_info_t::molecules[istat].update_map();
	 graphics_draw();
      }
   }
   return istat;
}



void set_map_radius_slider_max(float f) {
   graphics_info_t::map_radius_slider_max = f;
} 


// return 1 on "yes, it has a cell".
// 
int has_unit_cell_state(int imol) { 

   int istate = 0;
   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model() ||
	     graphics_info_t::molecules[imol].has_map()) {
	    istate = graphics_info_t::molecules[imol].have_unit_cell;
	 }
      }
   }
   return istate;
}

void set_symmetry_size_from_widget(const char *text) {
   
   float tmp;
   graphics_info_t g;

   tmp = atof(text);

   if ((tmp > 0.0) && (tmp < 9999.9)) {
      g.symmetry_search_radius = tmp;
   } else {

      cout << "Cannot interpret " << text << ".  Assuming 10A" << endl;
      g.symmetry_search_radius = 10.0;
   }
   //
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_symmetry();
   }
   graphics_draw();
   
}

void set_symmetry_size(float f) {
   graphics_info_t g;
   g.symmetry_search_radius = f;
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_symmetry();
   }
   graphics_draw();
}

/* When the coordinates for one (or some) symmetry operator are missing
   (which happens sometimes, but rarely), try changing setting this to 2
   (default is 1).  It slows symmetry searching, which is why it is not
   set to 2 by default.  */
void set_symmetry_shift_search_size(int shift) {
   
   graphics_info_t::symmetry_shift_search_size = shift;
}


void set_symmetry_molecule_rotate_colour_map(int imol, int state) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      g.molecules[imol].symmetry_rotate_colour_map_flag = state;
   }
   graphics_draw();
} 

int symmetry_molecule_rotate_colour_map_state(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].symmetry_rotate_colour_map_flag;
   }
   return r;
} 

void set_symmetry_colour_by_symop(int imol, int state) {

   if (graphics_info_t::use_graphics_interface_flag) { 
      graphics_info_t g;
      if (is_valid_model_molecule(imol)) { 
	 g.molecules[imol].symmetry_colour_by_symop_flag = state;
	 graphics_draw();
      }
   }
}

void set_symmetry_whole_chain(int imol, int state) {
   
   if (graphics_info_t::use_graphics_interface_flag) { 
      graphics_info_t g;
      if (is_valid_model_molecule(imol)) { 
	 g.molecules[imol].symmetry_whole_chain_flag = state;
	 if (g.glarea)
	    g.update_things_on_move_and_redraw();
      }
   }
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



void set_fps_flag(int thing) {

   graphics_info_t g;
   g.SetShowFPS(thing); 

} 

// For people without PCs with fast graphics cards :)
// 
void set_active_map_drag_flag(int t) {

   graphics_info_t g;

   g.SetActiveMapDrag(t);
}

int get_fps_flag() {

   graphics_info_t g;

   return g.GetFPSFlag();
} 


//
short int get_active_map_drag_flag() {

   graphics_info_t g;

   return g.GetActiveMapDrag();
}

void set_draw_hydrogens(int imol, int istate) {

   graphics_info_t g;
   
   if ((imol < g.n_molecules) && (imol >= 0)) {
      g.molecules[imol].set_draw_hydrogens_state(istate);
      graphics_draw();
   } else { 
      std::cout << "WARNING:: No such molecule number " << imol << "\n";
   } 
} 

void set_show_origin_marker(int istate) {
   graphics_info_t::show_origin_marker_flag = istate;
   graphics_draw();
} 

int  show_origin_marker_state() {
   return graphics_info_t::show_origin_marker_flag;
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





void set_last_map_contour_level(float level) {

   graphics_info_t g;
   g.set_last_map_contour_level(level);

}

void set_last_map_contour_level_by_sigma(float level) {

   graphics_info_t g;
   g.set_last_map_contour_level_by_sigma(level);
}

void set_last_map_sigma_step(float f) { 

   graphics_info_t g;
   g.set_last_map_sigma_step(f);

} 

 


void
handle_symmetry_colour_change(int mol, gdouble* col) {

   //
   graphics_info_t::symm_colour[0][0] = col[0];
   graphics_info_t::symm_colour[0][1] = col[1];
   graphics_info_t::symm_colour[0][2] = col[2];

   graphics_draw();

}

gdouble*
get_map_colour(int imol) {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   if (imol < graphics_info_t::n_molecules) { 
      if (graphics_info_t::molecules[imol].has_map()) { 
	 colour[0] = graphics_info_t::molecules[imol].map_colour[0][0]; 
	 colour[1] = graphics_info_t::molecules[imol].map_colour[0][1]; 
	 colour[2] = graphics_info_t::molecules[imol].map_colour[0][2];
      }
   }

   return colour;

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

//! \brief return the colour of the imolth map (e.g.: (list 0.4 0.6
//0.8). If invalid imol return #f.
// 
#ifdef USE_GUILE
SCM map_colour_components(int imol) {

   SCM r = SCM_BOOL(0);
   if (is_valid_map_molecule(imol)) {
      double rc = graphics_info_t::molecules[imol].map_colour[0][0]; 
      double gc = graphics_info_t::molecules[imol].map_colour[0][1]; 
      double bc = graphics_info_t::molecules[imol].map_colour[0][2]; 
      r = SCM_CAR(scm_listofnull);
      // put red at the front of the resulting list
      r = scm_cons(scm_double2num(bc), r);
      r = scm_cons(scm_double2num(gc), r);
      r = scm_cons(scm_double2num(rc), r);
   }
   return r; 
}
#endif


/* Functions for Cancel button on map colour selection  */
void save_previous_map_colour(int imol) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].save_previous_map_colour();
   } 

} 

void restore_previous_map_colour(int imol) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].restore_previous_map_colour();
   }
   graphics_draw();
}


// -------------------------------------------------------------------

gdouble*
get_symmetry_bonds_colour(int idummy) {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   colour[0] = graphics_info_t::symm_colour[0][0];
   colour[1] = graphics_info_t::symm_colour[0][1];
   colour[2] = graphics_info_t::symm_colour[0][2];
   return colour;
}





// In future the gui will usefully set a mol number and we
// will use that.
//
void set_show_symmetry_master(short int state) {

   // 
   graphics_info_t g;

   // show symmetry state is no longer part of the molecule(s).
   
      
//       g.molecules[ii].show_symmetry = state;

//       if ( state == 1 ) {
// 	 g.molecules[mol_no].update_symmetry();
//       }
//    }

   g.show_symmetry = state; 
   for (int ii=0; ii<g.n_molecules; ii++)
      if (is_valid_model_molecule(ii))
	 graphics_info_t::molecules[ii].update_symmetry();
   graphics_draw();

   if (state) { 
      // Now count the number of model molecules that have symmetry
      // available.  If there are none, then pop up a warning.

      int n_has_symm = 0;
      int n_model_molecules = 0;
      for (int ii=0; ii<g.n_molecules; ii++)
	 if (is_valid_model_molecule(ii)) {
	    n_model_molecules++;
	    mat44 my_matt;
	    int err = graphics_info_t::molecules[ii].atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
	    if (err == SYMOP_Ok) {
	       n_has_symm++;
	       break;
	    }
	 }
      if ((n_has_symm == 0) && (n_model_molecules > 0)) {
	 std::string s = "WARNING:: there are no model molecules\n"; 
	 s += " that can display symmetry.  \n\nCRYST1 problem?";
	 if (graphics_info_t::use_graphics_interface_flag) { 
	    GtkWidget *w = g.wrapped_nothing_bad_dialog(s);
	    gtk_widget_show(w);
	 }
      }
   }

}

void set_show_symmetry_molecule(int imol, short int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_symmetry(state);
      if (state)
	 graphics_info_t::molecules[imol].update_symmetry();
      graphics_draw();
   }

}


void symmetry_as_calphas(int mol_no, short int state) {

   graphics_info_t g;
   if (is_valid_model_molecule(mol_no)) { 
      g.molecules[mol_no].symmetry_as_calphas = state;
      g.molecules[mol_no].update_symmetry();
   }
   graphics_draw();

}

short int get_symmetry_as_calphas_state(int imol) {

   graphics_info_t g;
   int r = -1;
   if (is_valid_model_molecule(imol))
      r = g.molecules[imol].symmetry_as_calphas;
       
   return r;
} 


//  There is no dependence on the mol_no.  The intereface where we
//  pass (int mol_no) is kept, because molecule symmetry dependence
//  may re-surface in future.
// 
short int get_show_symmetry() {

   return graphics_info_t::show_symmetry; // master

}

GtkWidget *symmetry_molecule_controller_dialog() {

   graphics_info_t g;
   return g.wrapped_create_symmetry_controller_dialog();
} 

   
void
set_clipping_front(float v) {

   graphics_info_t::clipping_front = v;

   graphics_draw();

}


void
set_clipping_back(float v) {

   graphics_info_t::clipping_back = v;

   graphics_draw();
}


/*  ----------------------------------------------------------------------- */
/*                         Colour                                           */
/*  ----------------------------------------------------------------------- */

void
set_symmetry_colour_merge(int mol_no, float v) {

   graphics_info_t::symm_colour_merge_weight[mol_no] = v;
   graphics_draw();
}

void set_colour_map_rotation_on_read_pdb(float f) {
   graphics_info_t::rotate_colour_map_on_read_pdb = f; 
}

void set_colour_map_rotation_on_read_pdb_flag(short int i) {
   graphics_info_t::rotate_colour_map_on_read_pdb_flag = i; 
}

void set_colour_map_rotation_on_read_pdb_c_only_flag(short int i) {

   graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag = i;
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      if (is_valid_model_molecule(imol)) {
	 if (graphics_info_t::molecules[imol].Bonds_box_type() == coot::COLOUR_BY_CHAIN_BONDS) {
	    graphics_info_t::molecules[imol].make_bonds_type_checked();
	 }
      }
   }
   graphics_draw();
}

void set_symmetry_atom_labels_expanded(int state) {
   graphics_info_t::symmetry_atom_labels_expanded_flag = state;
   graphics_draw();
}


/* widget work */
GtkWidget *wrapped_create_coords_colour_control_dialog() {

   GtkWidget *w = create_coords_colour_control_dialog();

   graphics_info_t g;
   g.fill_bond_colours_dialog_internal(w);
   return w;
}


float get_molecule_bonds_colour_map_rotation(int imol) {

   float r = -1.0;
   if (is_valid_model_molecule(imol))
      r = graphics_info_t::molecules[imol].bonds_colour_map_rotation;
   return r;
}

void  set_molecule_bonds_colour_map_rotation(int imol, float f) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].bonds_colour_map_rotation = f;

}



// ----------------- Rotation Centre ----------------------

void set_rotation_centre(float x, float y, float z) {
   graphics_info_t g;
   g.setRotationCentre(coot::Cartesian(x,y,z));
   if (g.glarea)
      g.update_things_on_move_and_redraw();
}

// The redraw happens somewhere else...
void set_rotation_centre_internal(float x, float y, float z) {
   graphics_info_t g;
   g.setRotationCentre(coot::Cartesian(x,y,z));
}

float rotation_centre_position(int axis) {  /* only return one value: x=0, y=1, z=2 */
   graphics_info_t g;
   coot::Cartesian p = g.RotationCentre();
   // std::cout << "DEBUG:: rotation centre : " << p << " axis: " << axis << "\n";
   float r = 0.0;
   if (axis == 0)
      r = p.x();
   if (axis == 1)
      r = p.y();
   if (axis == 2)
      r = p.z();
   return r;
}

void set_colour_by_chain(int imol) { 
   
   if (is_valid_model_molecule(imol)) {
      short int f = graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag;
      graphics_info_t::molecules[imol].make_colour_by_chain_bonds(f);
      graphics_draw();
   }
}
void set_colour_by_molecule(int imol) { 

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].make_colour_by_molecule_bonds();
      graphics_draw();
   }
}


/*  Section Map colour*/
/* default for maps is 31 degrees. */
void set_colour_map_rotation_for_map(float f) {

   graphics_info_t::rotate_colour_map_for_map = f;
}



/*  ----------------------------------------------------------------------- */
/*                         Unit Cell                                        */
/*  ----------------------------------------------------------------------- */
short int
get_show_unit_cell(int imol) {

   return graphics_info_t::molecules[imol].show_unit_cell_flag;

}

void
set_show_unit_cell(int imol, short int state) {


   //    for (int imol=0; imol<graphics_n_molecules(); imol++) {
   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].show_unit_cell_flag = state;
   }
   //    }
   graphics_draw();
}

void set_show_unit_cells_all(short int istate) {

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) { 
	 graphics_info_t::molecules[imol].show_unit_cell_flag = istate;
      }
      if (is_valid_map_molecule(imol)) { 
	 graphics_info_t::molecules[imol].show_unit_cell_flag = istate;
      }
   }
   graphics_draw();

} 



// -----------------------------------------------------------------------
//                       Anisotropic Atoms
// -----------------------------------------------------------------------

void 
set_limit_aniso(short int state) {
   //
   graphics_info_t g; 
   
   g.show_aniso_atoms_radius_flag = state;
} 

void
set_aniso_limit_size_from_widget(const char *text) {
   
   float tmp;
   graphics_info_t g; 

   tmp = atof(text);

   if ((tmp >= 0.0) && (tmp < 99999.9)) {

      g.show_aniso_atoms_radius = tmp;
   } else {
      cout << "Cannot interpret " << text << ".  Assuming 10A" << endl;
      g.show_aniso_atoms_radius = 10.0;
   }
}

float
get_limit_aniso() {

   graphics_info_t g;
   
   return g.show_aniso_atoms_radius;

}

// Do if for all molecule, if we do it for one.
// The Anisotropic Atoms widget has no ability to select molecules.
// 
// This then is not a property of a molecule, but is a property of the
// graphics. 
// 
short int
get_show_aniso() {

   return graphics_info_t::show_aniso_atoms_flag;
}

void
set_show_aniso(int state) {

   graphics_info_t::show_aniso_atoms_flag = state;
   graphics_draw();
}

char *get_text_for_aniso_limit_radius_entry() {
   char *text;
   graphics_info_t g;
   
   text = (char *) malloc(100);
   snprintf(text, 99, "%-5.1f", g.show_aniso_atoms_radius);

   return text;
}

short int
get_show_limit_aniso() {
   
   return graphics_info_t::show_aniso_atoms_radius_flag;
}

void
set_aniso_probability(float f) {

   graphics_info_t::show_aniso_atoms_probability = f;
   graphics_draw();
}

// return e.g. 47.9 (%). 
float
get_aniso_probability() {

   return graphics_info_t::show_aniso_atoms_probability;
}


/*  ---------------------------------------------------------------------- */
/*                         Display Functions                               */
/*  ---------------------------------------------------------------------- */

void set_default_bond_thickness(int t) {

   graphics_info_t g;
   g.default_bond_width = t;

} 

void set_bond_thickness(int imol, float t) {

   graphics_info_t g;
   g.set_bond_thickness(imol, t);

}

void set_bond_thickness_intermediate_atoms(float t) { 

   graphics_info_t g;
   g.set_bond_thickness_intermediate_atoms(t);

} 


void set_unbonded_atom_star_size(float f) {
   graphics_info_t g;
   g.unbonded_atom_star_size = f;
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

int make_ball_and_stick(int imol, const char *atom_selection_str,
			float bond_thickness, float sphere_size,
			int do_spheres_flag) {

   int i = imol;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].make_ball_and_stick(std::string(atom_selection_str),
							   bond_thickness,
							   sphere_size, do_spheres_flag);
      graphics_draw();
   }
   return i;
}


int clear_ball_and_stick(int imol) {

   if (is_valid_model_molecule(imol)) {
      GLuint dummy_tag;
      graphics_info_t::molecules[imol].clear_display_list_object(dummy_tag);
      graphics_draw();
   }
   return 0;
} 

/*  ----------------------------------------------------------------------- */
/*                  dots display                                            */
/*  ----------------------------------------------------------------------- */
int dots(int imol,
	  const char *atom_selection_str,
	  float dot_density, float sphere_size_scale) {

   int idots = -1;
   if (is_valid_model_molecule(imol)) {
      if (atom_selection_str) { 
	 idots = graphics_info_t::molecules[imol].make_dots(std::string(atom_selection_str),
							    dot_density,
							    sphere_size_scale);
      }
   }
   graphics_draw();
   return idots;
}

void clear_dots(int imol, int dots_handle) {

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules)) { 
      graphics_info_t::molecules[imol].clear_dots(dots_handle);
      graphics_draw();
   }
}

/* return the number of dots sets for molecule number imol */
int n_dots_sets(int imol) {

   int r = -1;

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules)) { 
      r = graphics_info_t::molecules[imol].n_dots_sets();
   } else {
      std::cout << "WARNING:: Bad molecule number: " << imol << std::endl;
   } 
   return r;
} 


std::pair<short int, float> float_from_entry(GtkWidget *entry) {

   std::pair<short int, float> p(0,0);
   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      float f = atof(txt);
      p.second = f;
      p.first = 1;
   }
   return p;
}

std::pair<short int, int> int_from_entry(GtkWidget *entry) {

   std::pair<short int, int> p(0,0);
   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      int i = atoi(txt);
      p.second = i;
      p.first = 1;
   }
   return p;
}




// -----------------------------------------------------------------------
//                       Smooth Scrolling
// -----------------------------------------------------------------------

void set_smooth_scroll_flag(int v) {

   graphics_info_t::smooth_scroll = v;
}

int  get_smooth_scroll() {
   
   return graphics_info_t::smooth_scroll;
}

// useful interface for gui (entry)
void set_smooth_scroll_steps_str(const char *text) {

   int v;
   v = atoi(text);
   if (v > 0 && v < 10000000) {
      set_smooth_scroll_steps(v);
   } else {
      cout << "Cannot interpret " << text << ".  Assuming 10 steps" << endl;
      set_smooth_scroll_steps(10);
   }
}

// useful interface for scripting
void set_smooth_scroll_steps(int v) {
      graphics_info_t::smooth_scroll_steps = v;
} 

   
char *get_text_for_smooth_scroll_steps() {

   char *text;

   text = (char *) malloc(100);
   snprintf(text, 99, "%-5d", graphics_info_t::smooth_scroll_steps);

   return text;
}

// useful interface for gui (entry)
void  set_smooth_scroll_limit_str(const char *text) {

   float v;

   v = atof(text);

   if (v >0 && v < 1000) { 
      graphics_info_t::smooth_scroll_limit = v;
   } else {
      cout << text << " out of range: using 10A" << endl;
      graphics_info_t::smooth_scroll_limit = 10;
   }
}

// useful for scripting
void set_smooth_scroll_limit(float lim) {
   graphics_info_t::smooth_scroll_limit = lim; 
} 

char *get_text_for_smooth_scroll_limit() {

   char *text;
   
   text = (char *) malloc(100);
   snprintf(text, 99, "%-5.1f", graphics_info_t::smooth_scroll_limit);

   return text;
}

void set_stop_scroll_diff_map(int i) { 
   graphics_info_t::stop_scroll_diff_map_flag = i;
} 

void set_stop_scroll_iso_map(int i) { 
   graphics_info_t::stop_scroll_iso_map_flag = i;
}


void set_stop_scroll_diff_map_level(float f) { 
   graphics_info_t::stop_scroll_diff_map_level = f;
} 

void set_stop_scroll_iso_map_level(float f) { 
   graphics_info_t::stop_scroll_iso_map_level = f;
}



// -----------------------------------------------------------

void set_font_size(int i) { 

   graphics_info_t g;

   g.set_font_size(i);

}

int get_font_size() {

   return graphics_info_t::atom_label_font_size;

}


/*  ---------------------------------------------------------------------- */
/*                         Rotation Centre Cube Size                       */
/*  ---------------------------------------------------------------------- */

void set_rotation_centre_size_from_widget(const gchar *text) { 
   
   float val;
   graphics_info_t g;

   val = atof(text); 
   if ((val > 1000) || (val < 0)) { 
      cout << "Invalid cube size: " << text << ". Assuming 1.0A" << endl; 
      val = 1.0; 
   } 
   g.rotation_centre_cube_size = val; 
   graphics_draw();
}

void set_rotation_centre_size(float f) {
   graphics_info_t g;
   g.rotation_centre_cube_size = f;
   graphics_draw();
} 

gchar *get_text_for_rotation_centre_cube_size() { 
   
   char *text; 
   graphics_info_t g; 

   text = (char *)  malloc (100);
   snprintf(text, 90, "%-6.3f", g.rotation_centre_cube_size); 
   return text; 
}

short int
recentre_on_read_pdb() {
   return graphics_info_t::recentre_on_read_pdb; 
}

void
set_recentre_on_read_pdb(short int i) {
   graphics_info_t::recentre_on_read_pdb = i;
}

/*  ---------------------------------------------------------------------- */
/*                         orthogonal axes                                 */
/*  ---------------------------------------------------------------------- */
void set_draw_axes(int i) {
   graphics_info_t::draw_axes_flag = i;
} 


// -------------------------------------------------------------------------
//                        (density) iso level increment entry
// -------------------------------------------------------------------------
//

// imol is ignored.
//
char* get_text_for_iso_level_increment_entry(int imol) {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-6.4f", g.iso_level_increment);

   return text;

}

void set_iso_level_increment_from_text(const char *text, int imol) {

   float val;

   graphics_info_t g;

   val = atof(text);

   if ((val > 10000) || (val < -10000)) {
      cout << "Cannot interpret: " << text
	   << ".  Assuming 0.05 for increment" << endl;
      val  = 0.05;

   }

   cout << "setting iso_level_increment to " << val << endl; 
   g.iso_level_increment = val;

   graphics_draw();
}

void set_iso_level_increment(float val) { 
   graphics_info_t g;
   g.iso_level_increment = val;
} 

// imol is ignored.
//
char* get_text_for_diff_map_iso_level_increment_entry(int imol) {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-6.4f", g.diff_map_iso_level_increment);
   return text;

}

void set_diff_map_iso_level_increment_from_text(const char *text, int imol) {

   float val;
   graphics_info_t g;

   val = atof(text);

   if ((val > 10000) || (val < -10000)) {
      cout << "Cannot interpret: " << text
	   << ".  Assuming 0.005 for increment" << endl;
      val  = 0.005;
   } 
   g.diff_map_iso_level_increment = val;
   graphics_draw();
}

void set_diff_map_iso_level_increment(float val) { 
   graphics_info_t::diff_map_iso_level_increment = val;
} 

void set_map_sampling_rate_text(const char *text) {

   float val;
   val = atof(text);

   if ((val > 100) || (val < 1)) {
      cout << "Nonsense value: " << text
	   << ".  Assuming 1.5 for increment" << endl;
      val  = 1.5;
   }
   set_map_sampling_rate(val);

}

void set_map_sampling_rate(float r) {

   graphics_info_t g;
   g.map_sampling_rate = r;

}

char* get_text_for_map_sampling_rate_text() {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-5.4f", g.map_sampling_rate);
   return text;


}

float get_map_sampling_rate() {
   graphics_info_t g;
   return g.map_sampling_rate;
}


/* applies to the current map */
void change_contour_level(short int is_increment) { // else is decrement. 

   graphics_info_t g; 
   int s = g.scroll_wheel_map;

   if (g.molecules[s].max_xmaps > 0) {
      if (g.molecules[s].is_difference_map_p()) {
	 g.molecules[s].contour_level[0] +=
	    g.diff_map_iso_level_increment;
      } else {
	 // normal case
	 if (is_increment) { 
	    g.molecules[s].contour_level[0] += g.iso_level_increment;
	 } else {
	    g.molecules[s].contour_level[0] -= g.iso_level_increment;
	 }
      }
      g.molecules[s].update_map();
      graphics_draw();
      std::cout << "contour level of molecule [" << s << "]:  "
		<< g.molecules[s].contour_level[0] << std::endl;
   }
} 

GtkWidget *main_window() {
   return graphics_info_t::glarea; 
}; 

int graphics_n_molecules() {
   return graphics_info_t::n_molecules;
}

int next_map_for_molecule(int imol) { /* return a map number */
   return graphics_info_t::molecules[imol].next_free_map();
}

// imol is used imap is ignored.
// You can fix this anachronism one day if you like.  FIXME.
// 
void set_scrollable_map(int imol) {

   graphics_info_t g;
   if (is_valid_map_molecule(imol)) {
      int imap = 0; // ignored
      g.set_Scrollable_Map(imol, imap); // in graphics-info.h
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid molecule"
		<< " in set_scrollable_map\n";
   }
     
}
 
/*  ----------------------------------------------------------------------- */
/*                  utility function                                        */
/*  ----------------------------------------------------------------------- */
// return -1 if atom not found.
int atom_index(int imol, const char *chain_id, int iresno, const char *atom_id) {

   int index = -1;
   graphics_info_t g;
   std::string altconf("");
   std::string inscode("");

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules) { 
	 // return g.molecules[imol].atom_index(chain_id, iresno, atom_id);
	 index = g.molecules[imol].full_atom_spec_to_atom_index(std::string(chain_id),
								iresno,
								inscode,
								std::string(atom_id),
								altconf);
      }
   }

   return index;
}

// Refine zone needs to be passed atom indexes (which it then converts
// to residue numbers - sigh).  So we need a function to get an atom
// index from a given residue to use with refine_zone()
// 
int atom_index_first_atom_in_residue(int imol, const char *chain_id, 
				     int iresno, const char *ins_code) {

   int index = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g; 
      index = g.molecules[imol].atom_index_first_atom_in_residue(std::string(chain_id),
								 iresno,
								 std::string(ins_code));
   }
   return index;
} 

float median_temperature_factor(int imol) {

   float low_cut = 2.0;
   float high_cut = 100.0;
   short int low_cut_flag = 0;
   short int high_cut_flag = 0;

   float median = -1.0;
   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 median = coot::util::median_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
							graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
							low_cut, high_cut,
							low_cut_flag,
							high_cut_flag);
      } else {
	 std::cout << "WARNING:: molecule " << imol << " has no model\n";
      }
   } else {
      std::cout << "WARNING:: no such molecule as " << imol << "\n";
   }
   return median;
}

float average_temperature_factor(int imol) { 

   float low_cut = 2.0;
   float high_cut = 100.0;
   short int low_cut_flag = 0;
   short int high_cut_flag = 0;

   float av = -1.0;
   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 av = coot::util::average_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
						     graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
						     low_cut, high_cut,
						     low_cut_flag,
						     high_cut_flag);
      } else {
	 std::cout << "WARNING:: molecule " << imol << " has no model\n";
      }
   } else {
      std::cout << "WARNING:: no such molecule as " << imol << "\n";
   }
   return av;
}

char *centre_of_mass_string(int imol) {

   char *s = 0; // guile/SWIG sees this as #f
   if (is_valid_model_molecule(imol)) {
      realtype x, y, z;
      GetMassCenter(graphics_info_t::molecules[imol].atom_sel.atom_selection,
		    graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
		    x, y, z);
      std::string sc = "(";
      sc += graphics_info_t::float_to_string(x);
      sc += " ";
      sc += graphics_info_t::float_to_string(y);
      sc += " ";
      sc += graphics_info_t::float_to_string(z);
      sc += ")";
      s = new char [sc.length() + 1];
      strcpy(s, sc.c_str());
      return s;
   }
   return s;
}



void clear_pending_picks() {
   graphics_info_t g;
   g.clear_pending_picks();
}

#include "gl-matrix.h"
void print_view_matrix() { 		/* print the view matrix */

   graphics_info_t g;
   GL_matrix m;
   m.from_quaternion(g.quat);
   std::cout << "View Matrix:" << std::endl;
   m.print_matrix();
}

float get_view_matrix_element(int row, int col) {

   graphics_info_t g;
   GL_matrix m;
   m.from_quaternion(g.quat);
   return m.matrix_element(row, col);
}


float get_view_quaternion_internal(int element) {

   if ((element >= 0) &&
       (element < 4)) {
      return graphics_info_t::quat[element];
   } else {
      std::cout << "Bad element for quaternion: " << element
		<< " returning dummy -9999" << std::endl;
      return -9999;
   }
}

void set_view_quaternion(float i, float j, float k, float l) {

   double mag2 = i*i + j*j + k*k + l*l;
   double mag=sqrt(mag2);

   if (fabs(mag) > 0.5) {
      graphics_info_t::quat[0] = i/mag;
      graphics_info_t::quat[1] = j/mag;
      graphics_info_t::quat[2] = k/mag;
      graphics_info_t::quat[3] = l/mag;
      graphics_draw();
   } else {
      std::cout << "Bad view quaternion" << std::endl;
   } 
}



// ------------------------------------------------------
//                   Skeleton
// ------------------------------------------------------


void
handle_skeleton_colour_change(int mol, gdouble* map_col) {

   graphics_info_t::skeleton_colour[0] = map_col[0];
   graphics_info_t::skeleton_colour[1] = map_col[1];
   graphics_info_t::skeleton_colour[2] = map_col[2];

   graphics_draw();

}

gdouble*
get_skeleton_colour() {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   colour[0] = graphics_info_t::skeleton_colour[0]; 
   colour[1] = graphics_info_t::skeleton_colour[1]; 
   colour[2] = graphics_info_t::skeleton_colour[2]; 

   return colour;
}

void set_skeleton_colour(int imol, float r, float g, float b) {

   graphics_info_t::skeleton_colour[0] = r;
   graphics_info_t::skeleton_colour[1] = g;
   graphics_info_t::skeleton_colour[2] = b;

   graphics_draw();
}

void
skel_greer_on() {

   int i_skel_set = 0;
   graphics_info_t g; 
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    g.molecules[imol].greer_skeleton_draw_on = 1;
	    // g.molecules[imol].update_skeleton(); // withdrawn
	    i_skel_set = 1;
	    break;
	 }
      }
      if (i_skel_set) break;
   }
   graphics_draw();
}

void
skel_greer_off() {

   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      for (int imap=0; imap<graphics_info_t::molecules[imol].max_xmaps; imap++) {

	 if (graphics_info_t::molecules[imol].xmap_is_filled[imap] &&
	     graphics_info_t::molecules[imol].xmap_is_diff_map[imap] != 1) {

	    graphics_info_t::molecules[imol].greer_skeleton_draw_on = 0;
	 }
      }
   }
}

#include "graphical_skel.h"
#include "xmap-utils.h"

// For some as yet unknown reason, this code is executed when we
// select skeleton off.
// (After it has run, skel_foadi_off is executed).
// 
void
skel_foadi_on() {

   int i_found_skeletonizable_map = 0;
   graphics_info_t g;

   // we use break so that we get only one skeleton. 
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    i_found_skeletonizable_map = 1;

	    // so that we don't do this when the skeleton is on already:
	    //
	    if (g.molecules[imol].fc_skeleton_draw_on == 0) {
	       g.molecules[imol].fc_skeleton_draw_on = 1;

	       mean_and_variance<float> mv = 
		  map_density_distribution(g.molecules[imol].xmap_list[imap],0); 

	       cout << "Mean and sigma of map: " << mv.mean 
		    << " and " << sqrt(mv.variance) << endl; 

	       float map_cutoff = mv.mean + 1.5*sqrt(mv.variance); 
	       g.skeleton_level = map_cutoff; 
	    
	       // derived from sktest:
	       // 
	       g.molecules[imol].xskel_cowtan.init(g.molecules[imol].xmap_list[imap].spacegroup(), 
						   g.molecules[imol].xmap_list[imap].cell(),
						   g.molecules[imol].xmap_list[imap].grid_sampling());

	       cout << "INFO:: making skeleton cowtan..." << endl; 
	       GraphicalSkel cowtan(g.molecules[imol].xmap_list[0],
				    g.molecules[imol].xskel_cowtan); //fill xskel_cowtan

	       g.molecules[imol].xskel_is_filled = 1; // TRUE

	       // various experiments....

	       // cowtan.tip_filter(xmap_list[0], &xskl); // tinker with xskel_cowtan

	       //cowtan.prune(g.molecules[imol].xmap_list[imap],
	       //	 &g.molecules[imol].xskel_cowtan);

	       //
	       cout << "INFO:: pruning cowtan..." << endl; 
	       //
	       cowtan.Pprune(g.molecules[imol].xmap_list[imap],
			     &g.molecules[imol].xskel_cowtan,
			     map_cutoff);


	       // --------------------------------------------------------------
	       // --------------------------------------------------------------
	       // 	    cout << "cuckoo code:" << endl; 

	       // 	    BuildCas bc(g.molecules[imol].xmap_list[imap], map_cutoff); 

	       // 	    // mark segments by connectivity
	       // 	    // 
	       // 	    int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
	       // 						       g.molecules[imol].xmap_list[imap],
	       // 						       map_cutoff); 

	       // 	    cout << "There were " << nsegments << "segment in cuckoo map" << endl; 
	       // 	    // now transfer the segment map to xskel_cowtan
	       // 	    // 
	       // 	    bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	    


	       // --------------------------------------------------------------
	       // --------------------------------------------------------------


	       // now display the skeleton
	    
	       g.molecules[imol].update_clipper_skeleton();

	       // now create a new molecule and put the bonds
	       // into it
	       // 
	       // g.n_molecules++;
	       // 	    std::cout << "g.n_molecules is now " << g.n_molecules << endl;

	       // 
	       break;   // only do one (the first one). 
	    }
	 }
      }
      if (i_found_skeletonizable_map) break;
   }
   graphics_draw();
}

void
skel_foadi_off() {

   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      for (int imap=0; imap<graphics_info_t::molecules[imol].max_xmaps; imap++) {
	 
	 if (graphics_info_t::molecules[imol].xmap_is_filled[imap] &&
	     graphics_info_t::molecules[imol].xmap_is_diff_map[imap] != 1) {
	    
	    graphics_info_t::molecules[imol].fc_skeleton_draw_on = 0;
	 }
      }
   }
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

int 
skeletonize_map(int prune_flag, int imol) {

   graphics_info_t::skeletonize_map(prune_flag, imol);
   return 0;
}

int unskeletonize_map(int imol) {
   graphics_info_t::unskeletonize_map(imol);
   return imol;
} 


// Old unused code
// 
void 
autobuild_ca_on() { 

   graphics_info_t g; 

   g.autobuild_flag = 1;

   // we use break so that we get only one skeleton. 
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    if (g.molecules[imol].xskel_is_filled == 0) { 
	       
	       cout << "----------------------------------------" << endl; 
	       cout << "must Calculate the Map Skeleton first..." << endl; 
	       cout << "----------------------------------------" << endl; 
	       
	    } else {

	       float map_cutoff  = g.skeleton_level;
	       //
	       // save a pointer to this map btw.
	       // 
	       BuildCas bc(g.molecules[imol].xmap_list[imap], map_cutoff); 


	       // mark segments by connectivity
	       // 
	       int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
							  g.molecules[imol].xmap_list[imap],
							  map_cutoff); 

	       cout << "INFO:: There were " << nsegments << " different segments" << endl; 
	       
	       // --------------------------------------------------------------
	       // --------------------------------------------------------------
	       cout << "cuckoo code:" << endl; 

	       bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	       g.molecules[imol].update_clipper_skeleton();

	       // --------------------------------------------------------------
	       // --------------------------------------------------------------


	       // bc.depth_search_skeleton_testing_2(); 

	       // add branch points
	       //
	       cout << "INFO:: finding branch points..." << endl; 
	       //
	       vector<coot::Cartesian> branch_pts =
		  bc.find_branch_points(g.molecules[imol].xmap_list[imap],
					g.molecules[imol].xskel_cowtan,
					map_cutoff);
	       // 
	       cout << "INFO:: converting branch points to asc..." << endl; 
	       // 
	       atom_selection_container_t branch_pts_as_asc =
		  bc.convert_to_atoms(g.molecules[imol].xmap_list[imap],
				      branch_pts, "branch points");
 
	       // slight tangle here (internally to BuildCas,
	       // branch_points_symm_expanded is a vector<coot::Cartesian>, but
	       // we want an asc, so need to convert, so we need cell and
	       // symm (so we pass a const reference to the map). 
	       // 
	       atom_selection_container_t s_e_branch_pts_as_asc = 
		  bc.symmetry_expanded_branch_points(g.molecules[imol].xmap_list[imap]); 

	       cout << "INFO:: c-interface branch points conversion to atoms done!" << endl; 

	       cout << "INFO:: c-interface converting skeleton points to atoms..." << endl; 

// 	       asc_and_grids all_skels_pts_in_asu =
// 		  bc.all_skel_pts_in_asu(g.molecules[imol].xmap_list[imap],
// 					 g.molecules[imol].xskel_cowtan,
// 					 map_cutoff); // was 0.2

	       asc_and_grids all_skels_pts_in_asu =
		  bc.toplevel_skel_pts_in_asu(); // use internal segment_map

	       cout << "INFO:: c-interface expanding skeleton points by symmetry..." << endl; 

	       atom_selection_container_t big_ball =
		  bc.build_big_ball(g.molecules[imol].xmap_list[imap],
				    all_skels_pts_in_asu.asc, 
				    all_skels_pts_in_asu.grid_points); 

	       GraphicalSkel cowtan; 

	       int n_tips = cowtan.N_tips(g.molecules[imol].xmap_list[imap],
					  g.molecules[imol].xskel_cowtan,
					  map_cutoff);

 	       bc.interconnectedness(n_tips);

	       // now make that atom_selection_container for branch points:
	       // 

	       // Turn this back on when we have filled in molecule and map control.
	       // As it stood this code makes the molecule "active" (i.e. clickable)
	       // but the atoms were not displayed. 
	       // 
//  	       g.molecules[imol_new].initialize_coordinate_things_on_read_molecule(label);
//  	       g.molecules[imol_new].atom_sel = s_e_branch_pts_as_asc; 
//  	       g.molecules[imol_new].makebonds(0.1); // we don't want to join branch points
// 	       g.n_molecules++;

 	       std::string label = "branch points (symm expanded)";
	       asc_to_graphics(s_e_branch_pts_as_asc, label, ATOM_BONDS, 0.1); 

	       // debug_atom_selection_container(g.molecules[imol].atom_sel); 



	       // now make the atom_selection_container for the big ball
	       // viewable:
	       // 

// 	       g.molecules[imol_new+1].initialize_coordinate_things_on_read_molecule("big ball");
// 	       g.molecules[imol_new+1].atom_sel = big_ball; 
// 	       g.molecules[imol_new+1].make_ca_bonds(3.72, 3.85); //uses atom_sel
// 	       g.n_molecules++;
	       
	       // CA_BONDS !?
	       // 
	       asc_to_graphics(big_ball, "big ball", CA_BONDS, 3.72, 3.85); 


	       
	       // bc.ca_grow(g.molecules[imol].xmap_list[imap]); 
	       bc.ca_grow_recursive(); 

	       // show the cas
// 	       g.molecules[imol_new+2].initialize_coordinate_things_on_read_molecule("Auto-built C-alphas"); 
// 	       g.molecules[imol_new+2].atom_sel = bc.grown_Cas(); 
// 	       g.molecules[imol_new+2].make_ca_bonds(2.6, 4.3);  // Yikes! :-)
// 	       g.n_molecules++;

	       //
	       asc_to_graphics(bc.grown_Cas(), "Auto-built C-alphas", CA_BONDS, 2.6,4.3); 

	       bc.grown_Cas().mol->WritePDBASCII("autobuilt.pdb"); 

	    }
	 }
      }
   }
   graphics_draw();
}


void
do_skeleton_prune() { 

   graphics_info_t g;
   float map_cutoff  = g.skeleton_level;

   // we use break so that we get only one skeleton. 
   short int done_skel = 0;
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {
	    
	    if (g.molecules[imol].xskel_is_filled == 1) { 
	       
	       BuildCas bc(g.molecules[imol].xmap_list[imap], map_cutoff); 


	       // mark segments by connectivity
	       // 
	       int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
							  g.molecules[imol].xmap_list[imap],
							  map_cutoff); 

	       cout << "INFO:: There were " << nsegments << " different segments" << endl; 
	       
	       bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	       g.molecules[imol].update_clipper_skeleton();
	       done_skel = 1;
	       break;
	    }
	 }
      }
      if (done_skel) break;
   }
} 

void test_fragment() {

#ifdef HAVE_GSL
   graphics_info_t g;
   g.rotamer_graphs(0);

#endif // HAVE_GSL
}

// we redefine TRUE here somewhere...
// #include <gdk/gdkglconfig.h>
// #include <gtk/gtkgl.h>
// #include <gdk/x11/gdkglx.h>
// #include <gdk/x11/gdkglglxext.h>

int test_function(int i, int j) {

//    graphics_info_t g;
//    g.wrapped_create_symmetry_controller_dialog();
//    return 0;

   if (1) {
      graphics_info_t g;
      g.Geom_p()->hydrogens_connect_file("THH", "thh_connect.txt");
   }

   if (0) {

      // GTK2 GTkGLExt code
//       GdkGLContext *glcontext = gtk_widget_get_gl_context(graphics_info_t::glarea);
//       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(graphics_info_t::glarea);
//       GdkGLConfig *glconfig = gtk_widget_get_gl_config(graphics_info_t::glarea);
//       Display *dpy = gdk_x11_gl_config_get_xdisplay (glconfig);
      // Bool glXMakeCurrent(Display * dpy,
      //                     GLXDrawable  Drawable,
      //                     GLXContext  Context)
      // gdk_gl_glXMakeContextCurrent(dpy, gldrawable, glcontext);

      // bwah!
      // glXMakeCurrent(dpy, gldrawable, glcontext);

      // another way?
//       GtkWidget *w = graphics_info_t::glarea;
//       GdkGLContext *glcontext = gtk_widget_get_gl_context (w);
//       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (w);
//       int i = gdk_gl_drawable_gl_begin (gldrawable, glcontext);
//       std::cout << "DEBUG gdk_gl_drawable_gl_begin returns state: "
// 		<< i << std::endl;
//       return i;
   } 

   if (0) {
      int imol = i;
      if (is_valid_model_molecule(imol)) { 
	 const coot::residue_spec_t clicked_residue("A", 1);
	 short int is_n_term_addition = 1;
	 CAtom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[10];
	 CChain *chain_p = at->GetChain();
	 std::pair<bool, std::string> p = 
	    graphics_info_t::molecules[imol].residue_type_next_residue_by_alignment(clicked_residue, chain_p, is_n_term_addition);
	 if (p.first == 1) { 
	    std::cout << "next residue: " << p.second << std::endl;
	 } else {
	    std::cout << "no next residue found." << std::endl;
	 }
      }
   } 


   if (0) { 
      GtkWidget *w = wrapped_create_least_squares_dialog();
      gtk_widget_show(w);
   }
      

   if (0) { 
      std::vector<std::string> s;
      s.push_back("");
      s.push_back("123");
      s.push_back("123/456");
      s.push_back("123/456/");
      
      for (unsigned int i=0; i<s.size(); i++) { 
	 std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(s[i]);
	 std::cout << "For string :" << s[i] << ": split is :"
		   << p.first << ": :" << p.second << ":" << std::endl; 
      }

      std::string t = "/my/thing/int.mtz data/crystal/FWT data/crystal/PHWT";
      std::vector<std::string> v = coot::util::split_string(t, " ");

      std::cout << "splitting :" << t << ": on " << " " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
	 std::cout << "split " << i << " :" << v[i] << ":\n";
      }
   }
   return 0;
}

int write_connectivity(const char *monomer_name, const char *filename) {

   graphics_info_t g;
   return g.Geom_p()->hydrogens_connect_file(monomer_name, filename);
} 


void screendump_image(const char *filename) {

   GLint viewport[4];
   glGetIntegerv(GL_VIEWPORT, viewport);
   glPixelTransferi(GL_MAP_COLOR, GL_FALSE);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   glPixelStorei(GL_PACK_ALIGNMENT, 1);
   
   unsigned char* pixels = new unsigned char[viewport[2]*viewport[3]*IMAGEINFO_RGBA_SIZE];
   glReadPixels(0, 0, viewport[2], viewport[3], GL_RGBA, GL_UNSIGNED_BYTE, pixels);
   image_info iinfo(viewport[2], viewport[3], pixels, IMAGEINFO_RGBA);

   // file should be a ppm file for now, png when we add it to configure
   iinfo.invert();
   int istatus = 0;
   try { 
      istatus = iinfo.write(filename);
   }
   catch (...) {
      std::string s("Can't write that image format at the moment.\n");
      s += "ppm is suggested instead.";
      wrapped_nothing_bad_dialog(s);
   }
   delete [] pixels; // does iinfo copy the data or the pointer? possible crash.
   std::cout << "screendump_image status " << istatus << std::endl;
   if (istatus) {
      std::string s = "Screendump image ";
      s += filename;
      s += " written";
      graphics_info_t g;
      g.statusbar_text(s);
   }
}

void make_image_raster3d(const char *filename) {

   std::string r3d_name = filename;
   r3d_name += ".r3d";
   raster3d(r3d_name.c_str());
#ifdef USE_GUILE

   std::string cmd("(raytrace 'raster3d ");
   cmd += single_quote(r3d_name);
   cmd += " ";
   cmd += single_quote(filename);
   cmd += "'dummy 'dummy)";
   safe_scheme_command(cmd);
   
#endif   
}

void make_image_povray(const char *filename) {
   std::string pov_name = filename;
   pov_name += ".pov";
   povray(pov_name.c_str());
#ifdef USE_GUILE

   int x_size = graphics_info_t::glarea->allocation.width;
   int y_size = graphics_info_t::glarea->allocation.height;
   std::string cmd("(raytrace 'povray ");
   cmd += single_quote(pov_name);
   cmd += " ";
   cmd += single_quote(filename);
   cmd += " ";
   cmd += graphics_info_t::int_to_string(x_size);
   cmd += " ";
   cmd += graphics_info_t::int_to_string(y_size);
   cmd += ")";
   safe_scheme_command(cmd);

#endif   
}



void
autobuild_ca_off() { 

   graphics_info_t g; 
   g.autobuild_flag = 0; 

}

void
handle_read_ccp4_map(const char* filename, int is_diff_map_flag) {

   if (filename) { 
      std::string str(filename); 
      graphics_info_t g;

      int istate = g.molecules[g.n_molecules].read_ccp4_map(str, is_diff_map_flag); 

      if (istate > -1) { // not a failure
	 // std::cout << "successfully read map into molecule " << g.n_molecules << std::endl;
	 int imol = graphics_n_molecules();
	 std::string name = g.molecules[imol].dotted_chopped_name();
	 g.scroll_wheel_map = imol;  // change the current scrollable map.
	 g.n_molecules++; 
      } else { 
	 std::cout << "Read map " << str << " failed" << std::endl;
	 std::string s = "Read map ";
	 s += str;
	 s += " failed.";
	 g.statusbar_text(s);
      } 
      graphics_draw();
   } else {
      // error
      std::cout << "ERROR:: filename null in handle_read_ccp4_map\n";
   } 
}


/*  ------------------------------------------------------------------------ */
/*                         clipping */
/*  ------------------------------------------------------------------------ */

void do_clipping1_activate(){

   GtkWidget *clipping_window;
   GtkScale *hscale;
   GtkAdjustment *adjustment;
 
   /* connect this to displaying the new clipping window */

   clipping_window = create_clipping_window();

   hscale = GTK_SCALE(lookup_widget(clipping_window, "hscale1")); 
/*    gtk_scale_set_draw_value(hscale, TRUE);  already does */

   adjustment = GTK_ADJUSTMENT 
      (gtk_adjustment_new(0.0, -10.0, 20.0, 0.05, 4.0, 10.1)); 

   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   gtk_signal_connect (GTK_OBJECT (adjustment), "value_changed",
		       GTK_SIGNAL_FUNC (clipping_adjustment_changed), NULL);
   
   gtk_widget_show(clipping_window); 
}

void clipping_adjustment_changed (GtkAdjustment *adj, GtkWidget *window) { 

   /*    printf("Clipping adjustment: %f\n", adj->value); */

   set_clipping_front(adj->value);
   set_clipping_back (adj->value);
}
 


/*  ----------------------------------------------------------------------- */
/*                        virtual trackball                                 */
/*  ----------------------------------------------------------------------- */

void
vt_surface(int v){ 

   graphics_info_t g;
   g.set_vt_surface(v);
   std::vector<std::string> command_strings;
//    command_strings.push_back("vt-surface");
//    command_strings.push_back(graphics_info_t::int_to_string(v));
//    add_to_history(command_strings);
}

int vt_surface_status() { 

   graphics_info_t g;
   return g.vt_surface_status();
}

/*  ----------------------------------------------------------------------- */
/*                        save coordintes                                   */
/*  ----------------------------------------------------------------------- */

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

// return status 1 is good, 0 is fail.
int save_coordinates(int imol, const char *filename) { 

   int ierr = 0;
   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    ierr = graphics_info_t::molecules[imol].save_coordinates(filename);
	 }
      } 
   } 
   std::vector<std::string> command_strings;
   command_strings.push_back("save-coordinates");
   command_strings.push_back(coot::util::int_to_string(imol));
   command_strings.push_back(filename);
   add_to_history(command_strings);
   return ierr;
}

void set_save_coordinates_in_original_directory(int i) {

   graphics_info_t::save_coordinates_in_original_dir_flag = i;

}


/* access to graphics_info_t::save_imol for use in callback.c */
int save_molecule_number_from_option_menu() {
   return graphics_info_t::save_imol;
}

/* access from callback.c, not to be used in scripting, I suggest. */
void set_save_molecule_number(int imol) {
   graphics_info_t::save_imol = imol;
}



/*  ----------------------------------------------------------------------- */
/*                        .phs file reading                                 */
/*  ----------------------------------------------------------------------- */

void
read_phs_and_coords_and_make_map(const gchar *pdb_filename){

   // This function is the .phs equivalent of c.f. make_and_draw_map,
   // map_fill_from_mtz.  We have previously stored the phs_filename
   // in the static graphics_info_t.
   // 
   graphics_info_t g; 

   int imol = g.n_molecules;

   // don't forget that this is a map.
   //
   int istat = g.molecules[imol].make_map_from_phs(std::string(pdb_filename),
						   g.get_phs_filename());

   if (istat != -1) { 
      g.n_molecules++;
      graphics_draw();
   } else {
      // give us a warning message then
      std::string w = "Sadly, the cell or space group is not comprehensible in\n";
      w += "the pdb file: ";
      w += pdb_filename;
      w += "\n";
      w += "Can't make map from phs file.";
      GtkWidget *widget = wrapped_nothing_bad_dialog(w);
      gtk_widget_show(widget);
   } 
}

/*! \brief read a phs file, the cell and symm information is from
  previously read (most recently read) coordinates file

 For use with phs data filename provided on the command line */
int 
read_phs_and_make_map_using_cell_symm_from_previous_mol(const char *phs_filename) {

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int done_flag = 0;
   int r = -1;

   int imol_ref = -1;

   for (int i=graphics_info_t::n_molecules-1; i>=0; i--) {
      if (is_valid_model_molecule(i)) {
	 imol_ref = i;
	 break;
      }
   }

   if (imol_ref > -1) 
      r = read_phs_and_make_map_using_cell_symm_from_mol(phs_filename, imol_ref);

   return r;
}


/*! \brief read a phs file and use the cell and symm in molecule
  number imol and use the resolution limits reso_lim_low and
  reso_lim_high  */
int
read_phs_and_make_map_with_reso_limits(int imol_ref, const char* phs_filename,
				       float reso_lim_low, float reso_lim_high) {
   // This function is the .phs equivalent of c.f. make_and_draw_map,
   // map_fill_from_mtz.  We have previously stored the phs_filename
   // in the static graphics_info_t.
   // 
   graphics_info_t g;
   int imol = g.n_molecules;

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int done_flag = 0;
   int istat = -1; // returned value


   if (g.molecules[imol_ref].have_unit_cell) {
      // convert from a set of coordinates

      spacegroup.init(clipper::Spgr_descr(g.molecules[imol_ref].atom_sel.mol->GetSpaceGroup()));

      clipper::Cell_descr cell_d(g.molecules[imol_ref].atom_sel.mol->get_cell().a,
				 g.molecules[imol_ref].atom_sel.mol->get_cell().b,
				 g.molecules[imol_ref].atom_sel.mol->get_cell().c,
				 clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().alpha),
				 clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().beta),
				 clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().gamma));

      cell.init(cell_d);
      done_flag = 1;
   } else {
      // no conversion needed, just get from map
      if (g.molecules[imol_ref].has_map()) {
	 cell.init(g.molecules[imol_ref].xmap_list[0].cell().descr());
	 spacegroup.init(g.molecules[imol_ref].xmap_list[0].spacegroup().descr()); 
	 done_flag = 1;
      } 
   }

   if (done_flag) { 

      // don't forget that this is a map.
      //
      std::string phs_file(phs_filename);
      istat = g.molecules[imol].make_map_from_phs_using_reso(phs_file,
								 spacegroup,
								 cell, 
								 reso_lim_low, reso_lim_high);

      if (istat != -1) {
	 imol = istat;
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 std::string w = "Sadly, something bad happened reading phs file using\n";
	 w += "the molecule number ";
	 w += coot::util::int_to_string(imol_ref); 
	 w += "\n";
	 w += "Can't make map from phs file.";
	 GtkWidget *widget = wrapped_nothing_bad_dialog(w);
	 gtk_widget_show(widget);
      }
   } else {
      // give us a warning message then
      std::string w = "Sadly, the cell or space group is not comprehensible in\n";
      w += "the molecule number ";
      w += coot::util::int_to_string(imol_ref); 
      w += "\n";
      w += "Can't make map from phs file.";
      GtkWidget *widget = wrapped_nothing_bad_dialog(w);
      gtk_widget_show(widget);
   }

   return istat;
}



int
read_phs_and_make_map_using_cell_symm_from_mol(const char *phs_filename_str, int imol_ref) { 

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int got_cell_symm_flag = 0;
   int imol = -1;// set bad molecule initally
   
   graphics_info_t g; 
//       std::cout << "DEBUG:: read_phs_and_make_map_using_cell_symm_from_mol "
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().a << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().b << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().c << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().alpha << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().beta << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().gamma << "  " 
// 		<< std::endl;

   if (is_valid_model_molecule(imol_ref) || is_valid_map_molecule(imol_ref)) {
      if (g.molecules[imol_ref].have_unit_cell) { 

	 std::string s(g.molecules[imol_ref].atom_sel.mol->GetSpaceGroup());
	 if (s == "R 3") {
	    std::cout << "debug --------- symbol transformation ------\n";
	    s = "P 3*";
	 }
	 std::cout << "------------- :" << s << ":\n";
	 s = "P 3*";
	 clipper::Spgr_descr sgd(s);
	 spacegroup.init(sgd);

	 clipper::Cell_descr cell_d(g.molecules[imol_ref].atom_sel.mol->get_cell().a,
				    g.molecules[imol_ref].atom_sel.mol->get_cell().b,
				    g.molecules[imol_ref].atom_sel.mol->get_cell().c,
				    clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().alpha),
				    clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().beta),
				    clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().gamma));

	 cell.init(cell_d);
	 got_cell_symm_flag = 1;
      } else {
	 // get the cell/symm from a map:
	 if (g.molecules[imol_ref].has_map()) {
	    cell.init(g.molecules[imol_ref].xmap_list[0].cell().descr());
	    spacegroup.init(g.molecules[imol_ref].xmap_list[0].spacegroup().descr()); 
	    got_cell_symm_flag = 1;
	 } else { 

	 } 
      }

      if (got_cell_symm_flag) {
	 std::string phs_filename(phs_filename_str); 

	 imol = g.n_molecules; 
	 g.molecules[g.n_molecules].make_map_from_phs(spacegroup, cell, phs_filename);
	 g.n_molecules++; 
	 graphics_draw();
      }
   }

   return imol;
}


int
read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(int imol_ref) { 

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int done_flag = 0;
   int imol = -1; // bad molecule
   
   graphics_info_t g; 
//       std::cout << "DEBUG:: read_phs_and_make_map_using_cell_symm_from_mol "
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().a << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().b << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().c << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().alpha << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().beta << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().gamma << "  " 
// 		<< std::endl;

   if (is_valid_model_molecule(imol_ref) || is_valid_map_molecule(imol_ref)) { 

      if (g.molecules[imol_ref].have_unit_cell) {
	 // Coordinates block:
	 MyCMMDBManager *mol = g.molecules[imol_ref].atom_sel.mol;;
	 clipper::Spacegroup t_spgr = coot::util::get_spacegroup_from_symops(mol);
	 if (t_spgr.is_null()) {
	    std::cout << "WARNING:: Cant get spacegroup from coordinates!\n";
	 } else {
	    spacegroup = t_spgr; 
	    // If you don't want to use get_cell, use mol->GetCell();
	    clipper::Cell_descr cell_d(mol->get_cell().a,
				       mol->get_cell().b,
				       mol->get_cell().c,
				       clipper::Util::d2rad(mol->get_cell().alpha),
				       clipper::Util::d2rad(mol->get_cell().beta),
				       clipper::Util::d2rad(mol->get_cell().gamma));
   
	    cell.init(cell_d);
	    done_flag = 1;
	 }

      } else {
	 // Map block
	 if (g.molecules[imol_ref].has_map()) {
	    cell.init(g.molecules[imol_ref].xmap_list[0].cell().descr());
	    spacegroup.init(g.molecules[imol_ref].xmap_list[0].spacegroup().descr()); 
	    done_flag = 1;
	 }
      }

      if (done_flag) {
	 std::string phs_filename(graphics_get_phs_filename()); 

	 imol = g.n_molecules; 
	 g.molecules[g.n_molecules].make_map_from_phs(spacegroup, cell, phs_filename);
	 g.n_molecules++; 
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: Failed to get cell/symm - skipping.\n";
      } 
   }
   return imol;
}

int
read_phs_and_make_map_using_cell_symm(const char *phs_file_name,
				      const char *hm_spacegroup, float a, float b, float c,
				      float alpha, float beta, float gamma) { /*! in degrees */

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   graphics_info_t g;

   spacegroup.init(clipper::Spgr_descr(std::string(hm_spacegroup)));
   cell.init (clipper::Cell_descr(a, b, c, 
				  clipper::Util::d2rad(alpha), 
				  clipper::Util::d2rad(beta), 
				  clipper::Util::d2rad(gamma)));
   
   
   std::string phs_filename(phs_file_name); 

   int imol = g.n_molecules;
   g.molecules[g.n_molecules].make_map_from_phs(spacegroup, cell, phs_filename);
   g.n_molecules++; 
   graphics_draw();

   return imol;
}


void
graphics_store_phs_filename(const gchar *phs_filename) {

   graphics_info_t g;
   g.set_phs_filename(std::string(phs_filename));
}


const char *
graphics_get_phs_filename() {

   graphics_info_t g;
   return g.get_phs_filename().c_str(); 
}

short int possible_cell_symm_for_phs_file() {

   if (graphics_info_t::n_molecules == 0) { 
      return 0; 
   } else {
      return 1; 
   } 
}

// return a string to each of the cell parameters in molecule imol.
// 
gchar *get_text_for_phs_cell_chooser(int imol, char *field) { 

   // we first look in atomseletion

   graphics_info_t g;
   gchar *retval = NULL;
   retval = (gchar *) malloc(12); 
   int ihave_cell = 0; 
   float cell[6];
   const char *spgrp = NULL; 

   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    if (g.molecules[imol].have_unit_cell) { 

	       ihave_cell = 1; 

	       cell[0] = g.molecules[imol].atom_sel.mol->get_cell().a;
	       cell[1] = g.molecules[imol].atom_sel.mol->get_cell().b;
	       cell[2] = g.molecules[imol].atom_sel.mol->get_cell().c;
	       cell[3] = g.molecules[imol].atom_sel.mol->get_cell().alpha;
	       cell[4] = g.molecules[imol].atom_sel.mol->get_cell().beta;
	       cell[5] = g.molecules[imol].atom_sel.mol->get_cell().gamma;
	       spgrp   = g.molecules[imol].atom_sel.mol->GetSpaceGroup(); 

	    } else { 

	       if (g.molecules[imol].max_xmaps > 0) { 

		  ihave_cell = 1; 

		  clipper::Spacegroup spacegroup = g.molecules[imol].xmap_list[0].spacegroup(); 
		  clipper::Cell       ccell      = g.molecules[imol].xmap_list[0].cell();

		  cell[0] = g.molecules[imol].xmap_list[0].cell().a(); 
		  cell[1] = g.molecules[imol].xmap_list[0].cell().b(); 
		  cell[2] = g.molecules[imol].xmap_list[0].cell().c(); 
		  cell[3] = g.molecules[imol].xmap_list[0].cell().alpha() * RADTODEG; 
		  cell[4] = g.molecules[imol].xmap_list[0].cell().beta()  * RADTODEG; 
		  cell[5] = g.molecules[imol].xmap_list[0].cell().gamma() * RADTODEG; 

		  spgrp = spacegroup.descr().symbol_hm().c_str(); 
	       } 
	    } 
   

	    if (spgrp) { 
	       if ( ! (strcmp(field, "symm") ) ) { 
		  snprintf(retval, 11, "%-s", spgrp);
	       }
	       if ( ! (strcmp(field, "a") ) ) { 
		  snprintf(retval, 11, "%7.3f", cell[0]);   
	       } 
	       if ( ! (strcmp(field, "b") ) ) { 
		  snprintf(retval, 11, "%7.2f",  cell[1]);  
	       } 
	       if ( ! (strcmp(field, "c") ) ) { 
		  snprintf(retval, 11, "%7.2f",  cell[2]);  
	       } 
	       if ( ! (strcmp(field, "alpha") ) ) { 
		  snprintf(retval, 11, "%6.2f",   cell[3]); 
	       } 
	       if ( ! (strcmp(field, "beta") ) ) { 
		  snprintf(retval, 11, "%6.2f",  cell[4]);  
	       } 
	       if ( ! (strcmp(field, "gamma") ) ) { 
		  snprintf(retval, 11, "%6.2f",   cell[5]); 
	       }


	       if (! ihave_cell) { 
		  strcpy(retval, "  -  "); 
	       }
	    }
	 }
      }
   }
   return retval; 
}


/*  ----------------------------------------------------------------------- */
/*                        undo last move                                    */
/*  ----------------------------------------------------------------------- */
void undo_last_move() {

   graphics_info_t g;
   g.undo_last_move(); // does a redraw
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

// // A misnamed function.  The atom list is not filled, it is cleared.
// // 
void fill_go_to_atom_residue_and_atom_lists(GtkWidget *residue_gtklist,
 					    GtkWidget *atom_gtklist) {

    graphics_info_t g;
    g.fill_go_to_atom_residue_list(residue_gtklist);
}


// not really a button select, its a menu item select
void
save_molecule_coords_button_select(GtkWidget *item, GtkPositionType pos) { 

   // graphics_info_t g;
   std::cout << "INFO:: Save coords molecule now: " << pos << std::endl;
   graphics_info_t::save_imol = pos;
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

int go_to_atom_molecule_number() {
   graphics_info_t g;
   return g.go_to_atom_molecule();
}

char *go_to_atom_chain_id() {
   graphics_info_t g; 
   gchar *txt = (gchar *)malloc(100);
   strcpy(txt, g.go_to_atom_chain()); 
   return txt; 
}

char *go_to_atom_atom_name() {
   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_atom_name()); 
   return txt; 
}

int go_to_atom_residue_number() {
   graphics_info_t g;
   return g.go_to_atom_residue();
}

char *go_to_atom_ins_code() {
   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_ins_code()); 
   return txt; 
}

char *go_to_atom_alt_conf() {
   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_alt_conf()); 
   return txt; 
}


// Note that t3 is an atom name with (possibly) an altLoc tag (after the comma).
// 
int set_go_to_atom_chain_residue_atom_name(const char *t1, int iresno, const char *t3) {

   graphics_info_t g; 

   // so we need to split t3 if it has a comma
   // 
   std::string t3s(t3);
   std::string::size_type icomma = t3s.find_last_of(",");
   if (icomma == string::npos) {

      // there was no comma, conventional usage:
      g.set_go_to_atom_chain_residue_atom_name(t1, iresno, t3); 

   } else { 

      std::string atname = t3s.substr(0,icomma);
      std::string altloc = t3s.substr(icomma+1, t3s.length());
      g.set_go_to_atom_chain_residue_atom_name(t1, iresno,
					       atname.c_str(),
					       altloc.c_str());

   }

   int success = g.try_centre_from_new_go_to_atom(); 
   if (success) { 
      update_things_on_move_and_redraw();
      CAtom *at = 0; // passed but not used, it seems.
      GtkWidget *window = graphics_info_t::go_to_atom_window;
      if (window)
	 g.update_widget_go_to_atom_values(window, at);
   }
   return success; 
}

// A C++ function interface:
// 
int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec) {

   graphics_info_t g;

   g.set_go_to_atom_chain_residue_atom_name(atom_spec.chain.c_str(), 
					    atom_spec.resno,
					    atom_spec.insertion_code.c_str(), 
					    atom_spec.atom_name.c_str(),
					    atom_spec.alt_conf.c_str());

   int success = g.try_centre_from_new_go_to_atom(); 
   if (success)
      update_things_on_move_and_redraw(); 

   return success; 
}


int set_go_to_atom_chain_residue_atom_name_strings(const char *t1, const char *t2, const char *t3)
{
   int it2 = atoi(t2); 
   return set_go_to_atom_chain_residue_atom_name(t1, it2, t3); 
}

// FIXME to use altconf.
// 
// int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec) { 
   
//    return set_go_to_atom_chain_residue_atom_name(atom_spec.chain, 
// 						 atom_spec.resno,
// 						 atom_spec.atom_name);

// }


int 
goto_next_atom_maybe_new(GtkWidget *goto_atom_window) { 

//    int it2 = atoi(t2); 
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, +1);

   graphics_info_t g;
   return g.intelligent_next_atom_centring(goto_atom_window);
   
} 


int 
goto_previous_atom_maybe_new(GtkWidget *goto_atom_window) { 

//    int it2 = atoi(t2); 
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, +1);

   graphics_info_t g;
   return g.intelligent_previous_atom_centring(goto_atom_window);
   
} 
int 
goto_prev_atom_maybe(const gchar *t1, const gchar *t2, const gchar *t3,
		     GtkEntry *res_entry) { 

   int it2 = atoi(t2); 
   return goto_near_atom_maybe(t1, it2, t3, res_entry, -1); 
   
} 

int 
goto_near_atom_maybe(const char *t1, int ires, const char *t3,
		     GtkEntry *res_entry, int idiff) { 

   graphics_info_t g; 

   int ires_l = ires + idiff ; // for next residue, or previous.
   
   g.set_go_to_atom_chain_residue_atom_name(t1, ires_l, t3); 
   
   int success = g.try_centre_from_new_go_to_atom(); 

   if (success) {
      char *txt = (char *)malloc(6);
      snprintf(txt, 5, "%d", ires_l); 
      gtk_entry_set_text(GTK_ENTRY(res_entry), txt);
      update_things_on_move_and_redraw(); 
   }
   return success; 
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
     int gimol = g.go_to_atom_molecule();

     g.fill_option_menu_with_coordinates_options(option_menu,
						 callback_func,
						 gimol);

     /* These are in a special order: The residue is done first
	because it is set to a magic number (-9999 (or so)) initially.
	In that case, we do magic in
	get_text_for_go_to_atom_residue_entry(), i.e. look up a real
	atom of a molecule and set also the go to chain and the go to
	atom name  */

     /* The residue entry */

     residue_entry = lookup_widget(GTK_WIDGET(widget),
				   "go_to_atom_residue_entry"); 
     text = get_text_for_go_to_atom_residue_entry(); 
     gtk_entry_set_text(GTK_ENTRY(residue_entry), text); 

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

#else

     residue_tree = 0; // FIXME residue_tree = gtk_tree_view_new(); perhaps

#endif      


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
 			GTK_SIGNAL_FUNC(on_go_to_atom_residue_tree_selection_changed),
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
 			GTK_SIGNAL_FUNC(on_go_to_atom_atom_list_selection_changed),
 			NULL);

     /* fill those atom and residue lists (which uses
	graphics_info_t::go_to_atom_residue()) */
     fill_go_to_atom_residue_and_atom_lists(residue_tree,
					    atom_gtklist);

     /* store the widget */
     save_go_to_atom_widget(widget);

}

void apply_go_to_atom_from_widget(GtkWidget *widget) { 

   graphics_info_t g;
   g.apply_go_to_atom_from_widget(widget);
} 



/* For dynarama callback sake. The widget/class knows which coot
   molecule that it was generated from, so in order to go to the
   molecule from dynarama, we first need to the the molecule - because
   set_go_to_atom_chain_residue_atom_name() does not mention the
   molecule (see "Next/Previous Residue" for reasons for that).  This
   function simply calls the graphics_info_t function of the same
   name. */
void set_go_to_atom_molecule(int imol) {

   graphics_info_t g;
   g.set_go_to_atom_molecule(imol); 
   std::vector<std::string> command_strings;
   command_strings.push_back("set-go-to-atom-molecule");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);

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



// called by a function in callback.c
//
// void on_go_to_atom_residue_list_selection_changed (GtkList         *gtklist,
// 						   gpointer         user_data) {

void on_go_to_atom_residue_list_selection_changed (GtkList         *gtklist,
						   gpointer         user_data) {
   
//    graphics_info_t g;
//    g.on_go_to_atom_residue_list_selection_changed(gtklist, user_data);

   // old residue list function stub.
}

void on_go_to_atom_residue_tree_selection_changed (GtkList         *gtktree,
						   gpointer         user_data) {
   
   graphics_info_t g;
   g.on_go_to_atom_residue_tree_selection_changed(gtktree, user_data);
}

void clear_atom_list(GtkWidget *atom_gtklist) {

   gtk_list_clear_items(GTK_LIST(atom_gtklist), 0, -1);

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



void on_go_to_atom_atom_list_selection_changed(GtkList         *atom_gtklist,
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

void on_go_to_atom_atom_list_select_child (GtkList         *list,
					   GtkWidget       *widget,
					   gpointer         user_data) {
   std::cout << "child selected.\n";
}

void on_go_to_atom_atom_list_unselect_child (GtkList         *list,
					     GtkWidget       *widget,
					     gpointer         user_data) {
   std::cout << "child unselected.\n"; 
}


/*  ----------------------------------------------------------------------- */
/*                  bond representation                                     */
/*  ----------------------------------------------------------------------- */

void graphics_to_ca_representation(int imol) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol))
      g.molecules[imol].ca_representation();
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_ca_representation" << std::endl;
   graphics_draw();
   
   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
} 

void graphics_to_ca_plus_ligands_representation   (int imol) { 
   graphics_info_t g;
   g.molecules[imol].ca_plus_ligands_representation();
   graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-plus-ligands-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}


void graphics_to_bonds_representation(int imol) {
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      g.molecules[imol].bond_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-ca-plus-ligands-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_bonds_representation" << std::endl;
   graphics_draw();

}

void graphics_to_bonds_no_waters_representation(int imol) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)){ 
      g.molecules[imol].bonds_no_waters_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-no-waters-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_bonds_no_waters_representation"
		<< std::endl;
   graphics_draw();

} 

void graphics_to_sec_struct_bonds_representation(int imol) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      g.molecules[imol].bonds_sec_struct_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-sec-struct-bonds-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_sec_struct_bonds_representation"
		<< std::endl;
   graphics_draw();
} 

void graphics_to_ca_plus_ligands_sec_struct_representation(int imol) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      g.molecules[imol].ca_plus_ligands_sec_struct_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-ca-plus-ligands-sec-struct-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_ca_plus_ligands_sec_struct_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_rainbow_representation(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].ca_plus_ligands_rainbow_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-ca-plus-ligands-rainbow-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_ca_plus_ligands_rainbow_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_b_factor_representation(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].b_factor_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("b-factor-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_b_factor_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_occupancy_represenation(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].occupancy_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("occupancy-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_occupancy_representation"
		<< std::endl;
   graphics_draw();
}



int
graphics_molecule_bond_type(int imol) { 

   graphics_info_t g;
   // std::cout << "graphics_molecule_bond_type for mol: " << imol << std::endl;
   if (is_valid_model_molecule(imol)) { 
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-molecule-bond-type");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
      return g.molecules[imol].Bonds_box_type();
   }
   return -1;
}



// -------------------------------------------------------------------------
//                        skeletonization level
// -------------------------------------------------------------------------
//

gchar *get_text_for_skeletonization_level_entry() { 

   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10); 
   snprintf(txt, 9, "%f", g.skeleton_level); 

   return txt; 
} 

void set_skeletonization_level_from_widget(const char *txt) { 

   float tmp; 
   graphics_info_t g;

   tmp = atof(txt); 

   if (tmp > 0.0 &&  tmp < 9999.9) { 
      g.skeleton_level = tmp; 
   } else { 
      
      cout << "Cannot interpret " << txt << " using 0.2 instead" << endl; 
      g.skeleton_level = 0.2; 
   } 

         
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {
	 
	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    // 
	    g.molecules[imol].update_clipper_skeleton();

	 } 
      }
   }
   graphics_draw();
}


gchar *get_text_for_skeleton_box_size_entry() { 

   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10); 

   snprintf(txt, 9, "%f", g.skeleton_box_radius); 
   return txt;
} 

void set_skeleton_box_size_from_widget(const char *txt) { 
   float tmp; 
   graphics_info_t g;

   tmp = atof(txt); 

   if (tmp > 0.0 &&  tmp < 9999.9) { 
      g.skeleton_box_radius = tmp; 
   } else { 
      
      cout << "Cannot interpret " << txt << " using 0.2 instead" << endl; 
      g.skeleton_box_radius = 0.2; 
   }

   set_skeleton_box_size(g.skeleton_box_radius);
}

void set_skeleton_box_size(float f) {

   graphics_info_t g;
   g.skeleton_box_radius = f;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-skeleton-box-size");
   command_strings.push_back(graphics_info_t::float_to_string(f));
   add_to_history(command_strings);
      
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {
	 
	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    // 
	    g.molecules[imol].update_clipper_skeleton();

	 } 
      }
   }
   graphics_draw();
} 

/*  ----------------------------------------------------------------------- */
/*                  map and molecule control                                */
/*  ----------------------------------------------------------------------- */

void save_display_control_widget_in_graphics(GtkWidget *widget) { 

   graphics_info_t g; 
   g.save_display_control_widget_in_graphics(widget);
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

void 
post_display_control_window() { 

   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
   std::vector<std::string> command_strings;
   command_strings.push_back("post-display-control-window");
   add_to_history(command_strings);
}

GSList **gslist_for_scroll_in_display_manager_p() {

   return &graphics_info_t::gslist_for_scroll_in_display_manager;
}



 
void add_map_display_control_widgets() { 

   graphics_info_t g; 

   for (int ii=0; ii<g.n_molecules; ii++) {
      if (g.molecules[ii].max_xmaps > 0){ 
	 
	 g.molecules[ii].update_map_in_display_control_widget(); 
      } 
   }
}


void add_mol_display_control_widgets() { 

   graphics_info_t g; 
   
   for (int ii=0; ii<g.n_molecules; ii++) {
      if (! (g.molecules[ii].atom_sel.atom_selection == NULL)) { 
	 
	 g.molecules[ii].new_mol_in_display_control_widget(); 
      } 
   } 
} 


void add_map_and_mol_display_control_widgets() { 

   add_mol_display_control_widgets(); 
   add_map_display_control_widgets(); 
}


void reset_graphics_display_control_window() { 

   graphics_info_t g; 

   g.save_display_control_widget_in_graphics(NULL); 

} 

// widget (toggle button) call-backs
// 
int toggle_display_map(int imol, int imap) { 
   
   graphics_info_t g;
   int i = g.molecules[imol].toggle_display_map(imap); 
   graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("toggle-display-map");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   command_strings.push_back(graphics_info_t::int_to_string(imap));
   add_to_history(command_strings);
   return i;
} 


int toggle_display_mol(int imol) { 

   graphics_info_t g; 
   
   int i = g.molecules[imol].toggle_display_mol();
   graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("toggle-display-mol");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
   return i; 
}

// Return the new pickable? state.
// 
int toggle_active_mol(int imol) { 

   graphics_info_t g; 

   int i = g.molecules[imol].toggle_active_mol(); 

   return i; 

} 


int mol_is_displayed(int imol) { 

   graphics_info_t g; 

   return g.molecules[imol].drawit; 

} 

int mol_is_active(int imol) { 

   graphics_info_t g; 

   return g.molecules[imol].atom_selection_is_pickable(); 

} 

int map_is_displayed(int imol) { 

   graphics_info_t g;
   return g.molecules[imol].drawit_for_map; 

} 

char *show_spacegroup(int imol) { 

   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_map() || 
	  graphics_info_t::molecules[imol].has_model()) { 
	 std::string spg =  graphics_info_t::molecules[imol].show_spacegroup();
	 std::cout << "INFO:: spacegroup: " << spg << std::endl;
	 char *s = new char(spg.length()+1);
	 strncpy(s, spg.c_str(), spg.length()+1);
	 return s;
      } else { 
	 std::string spg("No spacegroup for this molecule");
	 std::cout << "INFO:: spacegroup: " << spg << std::endl;
	 char *s = new char(spg.length()+1);
	 strncpy(s, spg.c_str(), spg.length()+1);
	 return s;
      }
   }

   // If it was a bad molecule, return pointer to null.
   char *s = new char[1];
   s[0] = 0;
   return s;
}


/*  ----------------------------------------------------------------------- */
/*                  zoom                                                    */
/*  ----------------------------------------------------------------------- */

void 
scale_zoom_internal(float f) {

   graphics_info_t g;
   g.zoom *= fabs(f);
}

void scale_zoom(float f) {

   graphics_info_t g;
   scale_zoom_internal(f);
   graphics_draw();

}

float zoom_factor() {
   graphics_info_t g;
   return g.zoom;
}

void set_smooth_scroll_do_zoom(int i) {
   graphics_info_t g;
   g.smooth_scroll_do_zoom = i; 
}


int  smooth_scroll_do_zoom() {

   return graphics_info_t::smooth_scroll_do_zoom;
} 


float smooth_scroll_zoom_limit() {

   return graphics_info_t::smooth_scroll_zoom_limit; 
}


void set_smooth_scroll_zoom_limit(float f) {

   graphics_info_t::smooth_scroll_zoom_limit = f; 
}

void set_zoom_adjustment(GtkWidget *w) { 
   graphics_info_t::set_zoom_adjustment(w);
} 




// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_python_window_entry(GtkWidget *entry) { 

#ifdef USE_PYTHON

   // add python entery in entry callback code here...

      gtk_signal_connect(GTK_OBJECT(entry), "activate",
		      GTK_SIGNAL_FUNC(python_window_enter_callback),
		      entry);

#endif

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
  const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
  printf("Entry contents: %s\n", entry_text);
  {
     char *py_text;
     py_text = new char [strlen(entry_text)+1];
     strncpy(py_text, entry_text, 500); 
     //
     // Py_Initialize();
     PyRun_SimpleString(py_text);

     // clear the entry
     gtk_entry_set_text(GTK_ENTRY(entry),"");
  }

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


void clear_moving_atoms_object() {  /* redraw done here. */

   graphics_info_t g;
   g.clear_moving_atoms_object();
   
}

void clear_up_moving_atoms() {

   graphics_info_t g;
   std::cout << "c-interface clear_up_moving_atoms..." << std::endl;
   g.clear_up_moving_atoms();
   g.clear_moving_atoms_object();

} 

// Similar to fill_option_menu_with_coordinates_options, but I moved
// it to graphics_info_t because it is also used when there is an
// ambiguity in the map for refinement (graphics_info_t::refine)
// 
void fill_option_menu_with_map_options(GtkWidget *option_menu, GtkSignalFunc signalfunc) {

   graphics_info_t g;

   g.fill_option_menu_with_map_options(option_menu, signalfunc);
}

void fill_option_menu_with_skeleton_options(GtkWidget *option_menu) {  /* a wrapper */

   graphics_info_t g;
   GtkSignalFunc signalfunc = GTK_SIGNAL_FUNC(graphics_info_t::skeleton_map_select);
   g.fill_option_menu_with_map_options(option_menu, signalfunc,
				       graphics_info_t::map_for_skeletonize);

}


void set_initial_map_for_skeletonize() { 
   
   graphics_info_t::set_initial_map_for_skeletonize();

}

void set_max_skeleton_search_depth(int v) { 
   graphics_info_t g;
   g.set_max_skeleton_search_depth(v);
} 

/* Set the radio buttons in the frame to the be on or off for the map
   that is displayed in the optionmenu (those menu items "activate"
   callbacks (graphics_info::skeleton_map_select change
   g.map_for_skeletonize).  */
void set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame) { 

   graphics_info_t g;
   g.set_on_off_skeleton_radio_buttons(skeleton_frame);
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

void set_contour_by_sigma_step_by_mol(float f, short int state, int imol) { 
   
   if (imol < graphics_info_t::n_molecules) { 
      if (imol >= 0) { 
	 if (graphics_info_t::molecules[imol].has_map()) {
	    graphics_info_t::molecules[imol].set_contour_by_sigma_step(f, state);
	 }
      }
   }
}

void export_map(int imol, const char *filename) {

   if (is_valid_map_molecule(imol)) {

      clipper::CCP4MAPfile mapout;
      mapout.open_write(std::string(filename));
      mapout.export_xmap(graphics_info_t::molecules[imol].xmap_list[0]);
      mapout.close_write(); 
      
   } else {
      graphics_info_t g;
      g.statusbar_text("Invalid map molecule number");
   }

}

int transform_map_raw(int imol, 
		      double r00, double r01, double r02, 
		      double r10, double r11, double r12, 
		      double r20, double r21, double r22, 
		      double t0, double t1, double t2,
		      double pt1, double pt2, double pt3, double box_size) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol)) {
      clipper::Mat33<double> m(r00, r01, r02, r10, r11, r12, r20, r21, r22);
      clipper::Coord_orth c(t0, t1, t2);
      clipper::RTop_orth rtop(m,c);
      clipper::RTop_orth rtop_inv = rtop.inverse();
      clipper::Coord_orth pt(pt1, pt2, pt3);
      clipper::Xmap<float> new_map =
	 coot::util::transform_map(graphics_info_t::molecules[imol].xmap_list[0],
				   rtop, pt, box_size);

      const coot::ghost_molecule_display_t ghost_info;
      // int is_diff_map_flag = graphics_info_t::molecules[imol].is_difference_map_p();
      // int swap_colours_flag = graphics_info_t::swap_difference_map_colours;
      mean_and_variance<float> mv = map_density_distribution(new_map, 0);
      std::string name = "Transformed map";
      imol_new = graphics_info_t::n_molecules;
      graphics_info_t::molecules[imol_new].new_map(new_map, name);
      
      graphics_info_t::n_molecules++;
      graphics_draw();

   } else {
      std::cout << "molecule " << imol << " is not a valid map" << std::endl;
   }
   return imol_new;
}



void do_torsions_toggle(GtkWidget *button) {

   graphics_info_t g;
   GtkWidget *peptide_checkbutton =
      lookup_widget(button,
		    "refine_params_use_peptide_torsions_checkbutton");


   if (g.do_torsion_restraints) {
      g.do_torsion_restraints = 0;
      gtk_widget_set_sensitive(peptide_checkbutton, FALSE);
   } else {
      g.do_torsion_restraints = 1;
      gtk_widget_set_sensitive(peptide_checkbutton, TRUE);
   }
}

void do_peptide_torsions_toggle() {
   graphics_info_t g;
   if (g.do_peptide_torsion_restraints) {
      g.do_peptide_torsion_restraints = 0;
   } else {
      g.do_peptide_torsion_restraints = 1;
   }
}

void set_refine_params_toggle_buttons(GtkWidget *button) {

   // initiallly buttons are inactive and sensitive

   graphics_info_t g;
   GtkWidget *checkbutton =
      lookup_widget(button, "refine_params_use_torsions_checkbutton");
   GtkWidget *phi_psi_peptide_checkbutton =
      lookup_widget(button, "refine_params_use_peptide_torsions_checkbutton");
   GtkWidget *link_torsion_type_vbox =
      lookup_widget(button, "peptide_torsions_restraints_vbox");
   

   if (g.do_torsion_restraints) {
      g.do_torsion_restraints = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(phi_psi_peptide_checkbutton), FALSE);
   }

   if (g.do_peptide_torsion_restraints) {
      g.do_peptide_torsion_restraints = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(phi_psi_peptide_checkbutton), TRUE);
      gtk_widget_set_sensitive(link_torsion_type_vbox, TRUE);
   } else {
      gtk_widget_set_sensitive(link_torsion_type_vbox, FALSE);
   } 

   GtkWidget *omega = lookup_widget(button,
		       "refine_params_use_peptide_omegas_checkbutton");
   if (g.do_peptide_omega_torsion_restraints) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(omega), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(omega), FALSE);
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

   if (graphics_info_t::pseudo_bonds_type == coot::NO_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_no_rest_radiobutton), TRUE);
   if (graphics_info_t::pseudo_bonds_type == coot::HELIX_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_helix_rest_radiobutton), TRUE);
   if (graphics_info_t::pseudo_bonds_type == coot::STRAND_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_strand_rest_radiobutton), TRUE);

} 

// either alpha helix, beta strand or ramachandran goodness
// (see ideal/simple_restraint.hh link torsions)
void set_refine_params_phi_psi_restraints_type(int restraints_type) { 

  // do_peptide_omega_torsion_restraints are for phi/psi restraints 
  // (omega not included)
  graphics_info_t::do_peptide_torsion_restraints = restraints_type;

} 



void set_fix_chiral_volumes_before_refinement(int istate) {
   graphics_info_t::fix_chiral_volume_before_refinement_flag = istate;
} 

void check_chiral_volumes(int imol) { 
   graphics_info_t g;
   if (imol < graphics_info_t::n_molecules) { 
      if (graphics_info_t::molecules[imol].has_model()) { 
	 g.check_chiral_volumes(imol);
      } else {
	 std::cout << "WARNING:: molecule " << imol 
		   <<  " does not have coordinates\n";
      }
   } else {
      std::cout << "WARNING:: no such molecule " << imol << std::endl;
   } 
}

void check_chiral_volumes_from_widget(GtkWidget *window) { 
   
   check_chiral_volumes(graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule);
}


void fill_chiral_volume_molecule_option_menu(GtkWidget *w) { 

   GtkWidget *optionmenu = lookup_widget(w, "check_chiral_volumes_molecule_optionmenu");

   // now set chiral_volume_molecule_option_menu_item_select_molecule to the top of the list
   for (int i=0; i<graphics_info_t::n_molecules; i++) { 
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

void chiral_volume_molecule_option_menu_item_select(GtkWidget *item, GtkPositionType pos) { 

   graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule = pos;

}


void set_dragged_refinement_steps_per_frame(int v) {

   graphics_info_t g;
   g.dragged_refinement_steps_per_frame = v;
}

int dragged_refinement_steps_per_frame() {
   return graphics_info_t::dragged_refinement_steps_per_frame; 
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
   
   cout << "not compiled with guile.  This won't work" << endl; 

#endif // USE_GUILE

   // and kill the accession code window
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(widget),
				    "accession_code_window")); 
   
} 



/*  ----------------------------------------------------------------------- */
/*                  get molecule by libcheck/refmac code                    */
/*  ----------------------------------------------------------------------- */

/* Libcheck monomer code */
void 
handle_get_libcheck_monomer_code(GtkWidget *widget) { 

   const gchar *text = gtk_entry_get_text(GTK_ENTRY(widget));
   std::cout << "Refmac monomer Code: " << text << std::endl;
   int imol = get_monomer(text);

   // and kill the libcheck code window
   GtkWidget *window = lookup_widget(GTK_WIDGET(widget), "libcheck_monomer_dialog");
   if (window)
      gtk_widget_destroy(window);
   else 
      std::cout << "failed to lookup window in handle_get_libcheck_monomer_code" 
		<< std::endl;
}

// Return the new molecule number, or else a negitive error code.
// 
int get_monomer(const char *three_letter_code) {

   int imol = -1;

#ifdef USE_GUILE
   string scheme_command;

   scheme_command = "(monomer-molecule-from-3-let-code \"";

   scheme_command += three_letter_code;
   scheme_command += "\"";

   // now add in the bespoke cif library if it was given
   std::string cif_lib_filename = "";
   if (graphics_info_t::cif_dictionary_filename_vec->size() > 0)
      cif_lib_filename = (*graphics_info_t::cif_dictionary_filename_vec)[0];

   scheme_command += " ";
   std::string quoted_cif_lib_filename = single_quote(cif_lib_filename);
   scheme_command += quoted_cif_lib_filename;

   if (graphics_info_t::libcheck_ccp4i_project_dir != "") { 
      scheme_command += " ";
      scheme_command += single_quote(graphics_info_t::libcheck_ccp4i_project_dir);
   }

   scheme_command += ")";

   SCM v = safe_scheme_command(scheme_command);

   int was_int_flag = gh_scm2bool(scm_integer_p(v));

   if (was_int_flag)
      imol = gh_scm2int(v);

#else 
   
   std::cout << "not compiled with guile.  This won't work \n"
	     << "Need function to be coded in python..." << std::endl; 

#endif // USE_GUILE

   return imol;
} 


#ifdef USE_GUILE
SCM safe_scheme_command(const std::string &scheme_command) { 

   // FIXME!
   SCM handler = scm_c_eval_string ("(lambda (key . args) (display (list \"(safe_scheme_command) Error in proc: key: \" key \" args: \" args)) (newline))"); 

   // I am undecided if I want this or not:
   std::cout << "safe running: " << scheme_command << std::endl; 
   std::string thunk("(lambda() "); 
   thunk += scheme_command; 
   thunk += " )";

   SCM scm_thunk = scm_c_eval_string(thunk.c_str()); 
   SCM v = scm_catch(SCM_BOOL_T, scm_thunk, handler);

//   int is_int_p = scm_integer_p(v);
//   if (is_int_p) { 
//      std::cout << "returned value was int: " <<  std::endl;
//   } 

   add_to_history_simple(thunk);
  // std::cout << "INFO:: finished scheme command " << std::endl;
   return v;
}
#else  // not guile
// dummy function
void safe_scheme_command(const std::string &scheme_command) { /* do nothing */
   // here only for compilation purposes.
}
#endif // USE_GUILE

void safe_python_command(const std::string &python_cmd) {

#ifdef USE_PYTHON
   PyRun_SimpleString((char *)python_cmd.c_str());
#endif   
} 


void post_scripting_window() {

#ifdef USE_GUILE
  // window = create_guile_window(); 
//   entry = lookup_widget(window, "guile_window_entry"); 
//   setup_guile_window_entry(entry); // USE_PYTHON and USE_GUILE used here

  if (graphics_info_t::guile_gui_loaded_flag == TRUE) { 

     scm_c_eval_string("(coot-gui)");

  } else { 
     // we don't get a proper status from guile_gui_loaded_flag so
     // lets check again here whether MAPVIEW_GUI_DIR was defined.
     char *t; 
     t = getenv(COOT_SCHEME_DIR); // was #defined
     if (t) { 
	std::cout << COOT_SCHEME_DIR << " was defined to be " << t << std::endl
		  << "   but loading of scripting window scheme code failed." 
		  << std::endl; 
     } else { 
	std::cout << COOT_SCHEME_DIR << " was not defined - cannot open ";
	std::cout << "scripting window" << std::endl; 
     } 
  } 
#endif


#ifdef USE_PYTHON

  GtkWidget *window; 
  GtkWidget *entry; 
  window = create_python_window();

  entry = lookup_widget(window, "python_window_entry");
  setup_python_window_entry(entry); // USE_PYTHON and USE_GUILE used here
  gtk_widget_show(window);

  // clear the entry here
#endif

}

/* called from c-inner-main */
void run_command_line_scripts() {

//     std::cout << "There are " << graphics_info_t::command_line_scripts->size() 
//  	     << " command line scripts to run\n";

    for (unsigned int i=0; i<graphics_info_t::command_line_scripts->size(); i++)
       run_script((*graphics_info_t::command_line_scripts)[i].c_str());
}


void
set_guile_gui_loaded_flag() { 
   
   graphics_info_t g; 
   g.guile_gui_loaded_flag = TRUE; 
} 

void set_found_coot_gui() { 
   
   cout << "Coot Scripting GUI code found and loaded." << endl; 
   graphics_info_t g; 
   g.guile_gui_loaded_flag = TRUE; 
}

// return an atom index
int atom_spec_to_atom_index(int imol, char *chain, int resno, char *atom_name) { 
   graphics_info_t g; 
   if (imol < graphics_n_molecules()) 
      return g.molecules[imol].atom_spec_to_atom_index(chain, resno, atom_name);
   else
      return -1;
}

int full_atom_spec_to_atom_index(int imol, const char *chain, int resno,
				 const char *inscode, const char *atom_name,
				 const char *altloc) {

   if (imol < graphics_n_molecules()) 
      return graphics_info_t::molecules[imol].full_atom_spec_to_atom_index(std::string(chain), resno, std::string(inscode), std::string(atom_name), std::string(altloc));
   else
      return -1;
}



// ??? FIXME for the future, this doesn't work when we have both guile
// and python.  We need to choose which script interpretter to use
// based on filename extension.
void
run_script(const char *filename) { 

   struct stat buf;
   int status = stat(filename, &buf);
   std::string fn(filename);
   if (status == 0) {

      short int is_python = 0;

      std::string::size_type ipy = fn.rfind(".py");
      if (ipy != std::string::npos) {
	 if (fn.substr(ipy) == ".py")
	    is_python = 1;
      }
	 
      if (is_python) { 
	 run_python_script(filename);
      } else { 
	 run_guile_script(filename);
      }
   } else { 
      std::cout  << "WARNING:: Can't run script: " << filename 
		 << " no such file." << std::endl;
   }
}

// If we have both GUILE and PYTHON, use the state file as if it were GUILE
// 
void
run_state_file() { 
   std::string filename;
#ifdef USE_GUILE
   filename = "0-coot.state.scm";
   struct stat buf;
   int status = stat(filename.c_str(), &buf);
   if (status == 0) { 
      run_guile_script(filename.c_str());
   }
#else 
#ifdef USE_PYTHON
   filename = "0-coot.state.py";
   struct stat buf;
   int status = stat(filename.c_str(), &buf);
   if (status == 0) { 
      run_python_script(filename.c_str());
   }
#endif
#endif
}


void
run_state_file_maybe() { 

   std::string filename("0-coot.state.scm");
#ifdef USE_PYTHON
#ifndef USE_GUILE
   filename = "0-coot.state.py";
#endif
#endif
   graphics_info_t g;

   /*  0: never run it */
   /*  1: ask to run it */
   /*  2: always run it */
   if (g.run_state_file_status == 1 || g.run_state_file_status == 2) { 
      
      // can we stat a status file?
      // 
      struct stat buf;
      int status = stat(filename.c_str(), &buf);
      if (status == 0) { 
	 if (g.run_state_file_status == 2) {
	    run_script(filename.c_str());
	 } else {
	    if (graphics_info_t::use_graphics_interface_flag) { 
	       GtkWidget *dialog = wrapped_create_run_state_file_dialog();
	       gtk_widget_show(dialog);
	    } 
	 }
      }
   }
}

GtkWidget *wrapped_create_run_state_file_dialog() {

   std::string filename("0-coot.state.scm");
   short int il = 1;
   GtkWidget *w = create_run_state_file_dialog();

   GtkWidget *vbox_mols = lookup_widget(w, "mols_vbox");

   graphics_info_t g;
   std::vector<std::string> v = g.save_state_data_and_models(filename, il);
   for (unsigned int i=0; i<v.size(); i++) { 
      //       std::cout << "Got molecule: " << v[i] << std::endl;
      std::string s = "    ";
      s += v[i];
      GtkWidget *label = gtk_label_new(s.c_str());
      gtk_misc_set_alignment (GTK_MISC (label), 0.0, 0.5);
      gtk_box_pack_start(GTK_BOX(vbox_mols), label, FALSE, FALSE, 2);
      gtk_widget_show(label);
   } 
   return w;
}


void
run_guile_script(const char *filename) { 

#ifdef USE_GUILE
   std::string thunk("(lambda() "); 
   thunk += "(load ";
   thunk += "\"";
   thunk += filename;
   thunk += "\"))";

   SCM handler = scm_c_eval_string ("(lambda (key . args) "
     "(display (list \"Error in proc:\" key \" args: \" args)) (newline))"); 

   SCM scm_thunk = scm_c_eval_string(thunk.c_str());
   scm_catch(SCM_BOOL_T, scm_thunk, handler);
#endif // USE_GUILE   

} 

void
run_python_script(const char *filename_in) { 

#ifdef USE_PYTHON

   std::string s = coot::util::intelligent_debackslash(filename_in);
   std::string simple = "execfile(";
   simple += single_quote(s);
   simple += ")";
   std::cout << "Running python script " << s  << std::endl;
   // not a const argument?  Dear oh dear....
   PyRun_SimpleString((char *)simple.c_str());

#endif // USE_PYTHON
}



/*  ----------------------------------------------------------------------- */
/*                  ramachandran plot                                       */
/*  ----------------------------------------------------------------------- */



#include "rama_plot.hh"
// More dynarama stuff:
// call with value non-zero for on, 0 for off/not.
//
// This should not be used for scripting.
// 
void set_dynarama_is_displayed(GtkWidget *dyna_toplev, int imol) {

   graphics_info_t g;
   g.set_dynarama_is_displayed(dyna_toplev, imol);
}

GtkWidget *dynarama_is_displayed_state(int imol) {

   GtkWidget *w = NULL;
   if (is_valid_model_molecule(imol)) {
      w = graphics_info_t::dynarama_is_displayed[imol];
   }
   return w;
}


// window is the dynarama window.
// 
// Return -1 on error, return -9999 with a phi/psi edit window.
// 
int get_mol_from_dynarama(GtkWidget *window) {

   int imol = -1;
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS) 
   // graphics_info_t g; 
   if (window) {

      GtkWidget *canvas = lookup_widget(GTK_WIDGET(window), "canvas");

      if (canvas) { 
	 coot::rama_plot * plot =
	    (coot::rama_plot *)
	    gtk_object_get_user_data(GTK_OBJECT(canvas));
	 if (plot) 
	    imol = plot->molecule_number();

	 // This is no longer an error now that I have also added the
	 // destroy_window callback (which also does a
	 // set_dynarama_is_displayed)
	 // So just return -1 if the plot has already been deleted.
//c  	 else
//c  	    std::cout << "ERROR could not get user data from Ramachandran window" 
//c  		      << std::endl;
      } else {
	 std::cout << "ERROR:: failed to find canvas in window" << std::endl;
      } 
   }
#endif // HAVE_GTK_CANVAS   
   return imol; 
}
   


GtkWidget *dynarama_widget(int imol) {

   GtkWidget *w = NULL;
   if (imol < graphics_info_t::n_molecules) {
      w = graphics_info_t::dynarama_is_displayed[imol];
   }
   return w;
}


void add_on_rama_choices(){  // the the menu

   // std::cout << "adding rama molecule options:" << std::endl;

   // first delete all the current menu items.
   //
   graphics_info_t g;
   GtkWidget* menu = lookup_widget(GTK_WIDGET(g.glarea), "rama_plot_menu");

   if (menu) {
      gtk_container_foreach(GTK_CONTAINER(menu),
			    my_delete_ramachandran_mol_option,
			    (gpointer) menu);
   
      std::string name;
      for (int i=0; i<g.n_molecules; i++) {
	 if (g.molecules[i].has_model() > 0) {
	    name = graphics_info_t::molecules[i].dotted_chopped_name();
	    update_ramachandran_plot_menu_manual(i, name.c_str());
	 }
      }
   }
}

void
my_delete_ramachandran_mol_option(GtkWidget *widget, void *data) {
   gtk_container_remove(GTK_CONTAINER(data), widget);
}



void
set_moving_atoms(double phi, double psi) { 

   graphics_info_t g;
   g.set_edit_phi_psi_to(phi, psi);
}

void
accept_phi_psi_moving_atoms() { 

   graphics_info_t g;
   g.accept_moving_atoms();
   clear_moving_atoms_object();

}

void
setup_edit_phi_psi(short int state) {

   graphics_info_t g;
   g.in_edit_phi_psi_define = state;
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;

      std::cout << "click on an atom in the residue for phi/psi editting"
		<< std::endl;
   } else {
      g.normal_cursor();
   } 
}

void destroy_edit_backbone_rama_plot() { 

   graphics_info_t g;
   g.destroy_edit_backbone_rama_plot();

} 


/*  ----------------------------------------------------------------------- */
/*           sequence_view                                                  */
/*  ----------------------------------------------------------------------- */

// A pure sequence function, not sequence view, so that people can cut
// n paste the sequence of a pdb file from the console.  There will be
// a scripting level function called print-sequence that gets called
// for every chain in the mol
// 
void print_sequence_chain(int imol, const char *chain_id) {

   std::string seq;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (std::string(chain_p->GetChainID()) == chain_id) { 
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    int residue_count_block = 0;
	    int residue_count_line = 0;
	    if (nres > 0 ) {
	       residue_count_block = chain_p->GetResidue(0)->GetSeqNum();
	       residue_count_line  = residue_count_block;
	       if (residue_count_block > 0)
		  while (residue_count_block > 10)
		     residue_count_block -= 10;
	       if (residue_count_line > 0)
		  while (residue_count_line > 50)
		     residue_count_line -= 50;
	    }
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       seq += coot::util::three_letter_to_one_letter(residue_p->GetResName());
	       if (residue_count_block == 10) {
		  seq += " ";
		  residue_count_block = 0;
	       }
	       if (residue_count_line == 50) {
		  seq += "\n";
		  residue_count_line = 0;
	       }
	       residue_count_block++;
	       residue_count_line++;
	    }
	 }
      }
      std::cout << "> " << graphics_info_t::molecules[imol].name_sans_extension(0)
		<< " chain " << chain_id << std::endl;
      std::cout << seq << std::endl;
   }
}

void do_sequence_view(int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   graphics_info_t g;
   if (g.molecules[imol].has_model()) {
      graphics_info_t g;

      if (g.sequence_view_is_displayed[imol] != 0) {

	 // it already exists... just raise it and map it.

	 GtkWidget *canvas = g.sequence_view_is_displayed[imol];
	 // so what is the window (which we shall call widget)?
	 GtkWidget *widget = lookup_widget(canvas, "sequence_view_dialog");

	 if (!GTK_WIDGET_MAPPED(widget)) {
	    gtk_widget_show(widget);
	 } else {
	    gdk_window_raise(widget->window);
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
      for(int i=0; i<g.n_molecules; i++) {
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
/*           rotate moving atoms peptide                                    */
/*  ----------------------------------------------------------------------- */

void change_peptide_carbonyl_by(double angle) { /* in degrees. */
   graphics_info_t g;
   g.change_peptide_carbonyl_by(angle);
} 

void change_peptide_peptide_by(double angle) {   /* in degress */
   graphics_info_t g;
   g.change_peptide_peptide_by(angle);
}

void execute_setup_backbone_torsion_edit(int imol, int atom_index) {
   graphics_info_t g;
   g.execute_setup_backbone_torsion_edit(imol, atom_index);
}

void setup_backbone_torsion_edit(short int state) { 

   graphics_info_t g;
   graphics_info_t::in_backbone_torsion_define = state;
   if (state) { 
      std::cout << "click on an atom in the peptide to change" << std::endl; 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
   } else { 
      g.normal_cursor();
   }
}

void set_refine_with_torsion_restraints(int istate) {

   graphics_info_t::do_torsion_restraints = istate;

} 

   

void set_backbone_torsion_peptide_button_start_pos(int ix, int iy) { 
   
   graphics_info_t g;
   g.set_backbone_torsion_peptide_button_start_pos(ix, iy);
} 

void change_peptide_peptide_by_current_button_pos(int ix, int iy) { 

   graphics_info_t g;
   g.change_peptide_peptide_by_current_button_pos(ix, iy);
}

void set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy) { 

   graphics_info_t g;
   g.set_backbone_torsion_carbonyl_button_start_pos(ix, iy);

} 

void change_peptide_carbonyl_by_current_button_pos(int ix, int iy) { 

   graphics_info_t g;
   g.change_peptide_carbonyl_by_current_button_pos(ix, iy);

} 

/*  ----------------------------------------------------------------------- */
/*                  cif stuff                                               */
/*  ----------------------------------------------------------------------- */

// and make (and display) a sigma_a map.
// 
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
// 
int read_cif_data(const char *filename, int imol_coordinates) {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      cout << "Reading cif file: " << filename << endl; 
      graphics_info_t g; 
      int imol = g.n_molecules;
      int istat = g.molecules[imol].make_map_from_cif(std::string(filename), imol_coordinates);

      // std::cout << "DEBUG:: in read_cif_data, istat is " << istat << std::endl;
      if (istat != -1) { 
	 g.n_molecules++;
	 graphics_draw();
      }
      return imol;
   }
   return -1; // which is status in an error
}

// and make (and display) a 2fofc map.
// 
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
// 
int read_cif_data_2fofc_map(const char *filename, int imol_coordinates) {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      cout << "Reading cif file: " << filename << endl; 

      graphics_info_t g; 

      int imol = g.n_molecules;

      int istat = g.molecules[imol].make_map_from_cif_2fofc(std::string(filename), imol_coordinates);

      if (istat != -1) { 
	 g.n_molecules++;
	 graphics_draw();
	 return imol;
      }
      return -1; // an error
   }
}


// and make (and display) a fofc map.
// 
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
// 
int read_cif_data_fofc_map(const char *filename, int imol_coordinates) {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      cout << "Reading cif file: " << filename << endl; 

      graphics_info_t g; 

      int imol = g.n_molecules;

      int istat = g.molecules[imol].make_map_from_cif_fofc(std::string(filename), imol_coordinates);

      if (istat != -1) { 
	 g.n_molecules++;
	 graphics_draw();
	 return imol;
      }
      return -1; // an error
   }
}



// This cif file, we presume, has phases.
// So we don't need a molecule to calculate them from.
// 
int auto_read_cif_data_with_phases(const char *filename) {

   int returned_mol_index = read_cif_data_with_phases_sigmaa(filename);
   read_cif_data_with_phases_diff_sigmaa(filename);
   return returned_mol_index;
}

int read_cif_data_with_phases_sigmaa(const char *filename) {

   graphics_info_t g; 
   int imol = g.n_molecules;
   
   int returned_mol_index = -1;
   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      std::cout << "Reading cif file: " << filename << std::endl; 

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.
      std::string fn = filename;
      int istat = g.molecules[imol].make_map_from_cif(fn);
      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 imol = -1;
      }
   }
   return imol;
}

int read_cif_data_with_phases_diff_sigmaa(const char *filename) {

   graphics_info_t g; 
   int imol = g.n_molecules;
   
   int returned_mol_index = -1;
   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      std::cout << "Reading cif file: " << filename << std::endl; 

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.
      std::string fn = filename;
      int istat = g.molecules[imol].make_map_from_cif_diff_sigmaa(fn);
      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 imol = -1;
      }
   }
   return imol;
} 


int read_cif_data_with_phases_fo_fc(const char *filename) {

   return read_cif_data_with_phases_nfo_fc(filename, molecule_map_type::TYPE_FO_FC);

} 

int read_cif_data_with_phases_2fo_fc(const char *filename) {

   return read_cif_data_with_phases_nfo_fc(filename, molecule_map_type::TYPE_2FO_FC);
}

int read_cif_data_with_phases_fo_alpha_calc(const char *filename) {
   return read_cif_data_with_phases_nfo_fc(filename, molecule_map_type::TYPE_FO_ALPHA_CALC);
}

int read_cif_data_with_phases_nfo_fc(const char *filename,
				     int map_type) {
   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

      graphics_info_t g; 

      int imol = g.n_molecules;
      std::string f(filename);
      short int swap_col = graphics_info_t::swap_difference_map_colours;

      int istat = g.molecules[imol].make_map_from_cif_nfofc(f, map_type, swap_col);

      if (istat != -1) {

	 g.scroll_wheel_map = g.n_molecules; // change the current scrollable map.
	 g.n_molecules++;
	 graphics_draw();
	 return imol;
      }
      return -1; // error
   }
}

int handle_shelx_fcf_file_internal(const char *filename) {

   graphics_info_t g;
   std::vector<std::string> cmd;
   cmd.push_back("handle-shelx-fcf-file");
   cmd.push_back(single_quote(filename));

#ifdef USE_GUILE   
   std::string s = g.state_command(cmd, coot::STATE_SCM);
   safe_scheme_command(s);
#endif

#ifdef USE_PYTHON
#ifndef USE_GUILE
   std::string s = g.state_command(cmd, coot::STATE_PYTHON);
   safe_python_command(s);
#endif    
#endif    
   return 1; // needed so that try_read_cif_file succeeds, so that we
	     // don't try to read this file as CNS data.
}

/*  ----------------------------------------------------------------------- */
/*                  CNS data stuff                                          */
/*  ----------------------------------------------------------------------- */
int handle_cns_data_file(const char *filename) {

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      return -1; // which is status in an error
   } else {
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << std::endl;
      } else {
	 std::cout << "FIXME:: Fill this CNS data stub." << std::endl;
      } 
   }
   return 0;
}


void set_residue_density_fit_scale_factor(float f) {

   graphics_info_t::residue_density_fit_scale_factor = f;
}

float residue_density_fit_scale_factor() {
   return graphics_info_t::residue_density_fit_scale_factor; 
}


// dictionary
void handle_cif_dictionary(const char *filename) {

   graphics_info_t g;
   g.add_cif_dictionary(filename, 1); // show dialog if no bonds

}

void read_cif_dictionary(const char *filename) { 
   
   handle_cif_dictionary(filename);

} 

/* Use the environment variable COOT_REFMAC_LIB_DIR to find cif files
   in subdirectories and import them all. */
void import_all_refmac_cifs() {

   graphics_info_t g;
   g.import_all_refmac_cifs();

} 



/*  ----------------------------------------------------------------------- */
/*                  atom labelling                                          */
/*  ----------------------------------------------------------------------- */
/* The guts happens in molecule_class_info_t, here is just the
   exported interface */
int add_atom_label(int imol, char *chain_id, int iresno, char *atom_id) {

   int i = 0;
   if (is_valid_model_molecule(imol)) { 
      graphics_info_t g;
      i = g.molecules[imol].add_atom_label(chain_id, iresno, atom_id);
   }
   return i;
} 

int remove_atom_label(int imol, char *chain_id, int iresno, char *atom_id) {
   graphics_info_t g;
   return g.molecules[imol].remove_atom_label(chain_id, iresno, atom_id);
} 

void remove_all_atom_labels() { 
   graphics_info_t g;
   g.remove_all_atom_labels();
} 

void set_label_on_recentre_flag(int i) { 
   graphics_info_t::label_atom_on_recentre_flag = i;
} 

int centre_atom_label_status() { 
   return graphics_info_t::label_atom_on_recentre_flag;
} 

void set_brief_atom_labels(int istat) {
   graphics_info_t::brief_atom_labels_flag = istat;
   graphics_draw();
}

int brief_atom_labels_state() {
   return graphics_info_t::brief_atom_labels_flag;
}

/*  ----------------------------------------------------------------------- */
/*                  scene rotation (by script)                              */
/*  ----------------------------------------------------------------------- */
/* stepsize in degrees */
void rotate_y_scene(int nsteps, float stepsize) { 

  float spin_quat[4];
   graphics_info_t g;

   // spin it 1 degree
   float tbs =  g.get_trackball_size(); 
   for(int i=0; i<nsteps; i++) { 
     trackball(spin_quat, 0, 0, 0.0174*stepsize, 0.000, tbs);
     add_quats(spin_quat, g.quat, g.quat);
     graphics_draw();
   }
} 

/* stepsize in degrees */
void rotate_x_scene(int nsteps, float stepsize) { 

  float spin_quat[4];
   graphics_info_t g;

   // spin it 1 degree
   float tbs =  g.get_trackball_size(); 
   for(int i=0; i<nsteps; i++) { 
     trackball(spin_quat, 0, 0, 0.0, 0.0174*stepsize, tbs);
     add_quats(spin_quat, g.quat, g.quat);
     graphics_draw();
   }
} 

void rotate_z_scene(int nsteps, float stepsize) { 

   // c.f globjects.cc:do_screen_z_rotate()
   // 

   float spin_quat[4];
   graphics_info_t g;
   for(int i=0; i<nsteps; i++) { 
      trackball(spin_quat, 
		1.0, 1.0,
		1.0, 1.0 + 0.0174*stepsize,
		0.4);
      add_quats(spin_quat, g.quat, g.quat);
      graphics_draw();
   }
}

/*! \brief Bells and whistles rotation 

    spin, zoom and translate.

    where axis is either x,y or z,
    stepsize is in degrees, 
    zoom_by and x_rel etc are how much zoom, x,y,z should 
            have changed by after nstep steps.
*/
void spin_zoom_trans(int axis, int nsteps, float stepsize, float zoom_by, 
		     float x_rel, float y_rel, float z_rel) {

   float spin_quat[4];
   graphics_info_t g;
   float tbs =  g.get_trackball_size();
   float x_frag = 0.0;
   float y_frag = 0.0;
   float z_frag = 0.0;
   if (nsteps != 0) { 
      x_frag = x_rel/float(nsteps);
      y_frag = y_rel/float(nsteps);
      z_frag = z_rel/float(nsteps);
   }
   float zoom_init = g.zoom;
   float zoom_final = g.zoom * zoom_by;
   float zoom_frag = 1.0;
   if (nsteps !=0) {
      zoom_frag = (zoom_final - zoom_init)/float(nsteps);
   }
   
   std::cout << "zoom_frag is " << zoom_frag << std::endl;
   for(int i=0; i<nsteps; i++) {
      if (axis == 1) { 
	 trackball(spin_quat, 0, 0, 0.0, 0.0174*stepsize, tbs);
	 add_quats(spin_quat, g.quat, g.quat);
      }
      if (axis == 2) { 
	 trackball(spin_quat, 0, 0, 0.0174*stepsize, 0.000, tbs);
	 add_quats(spin_quat, g.quat, g.quat);
      }
      if (axis == 3) { 
	 trackball(spin_quat, 
		   1.0, 1.0,
		   1.0, 1.0 + 0.0174*stepsize,
		   0.4);
	 add_quats(spin_quat, g.quat, g.quat);
      }
      g.zoom = zoom_init + float(i+1)*zoom_frag;
      coot::Cartesian c(g.X() + x_frag, g.Y() + y_frag, g.Z() + z_frag);
      g.setRotationCentre(c);
      graphics_draw();
   }
} 




/*  ----------------------------------------------------------------------- */
/*                  graphics background colour                              */
/*  ----------------------------------------------------------------------- */
/* stepsize in degrees */
void set_background_colour(double red, double green, double blue) {

   graphics_info_t g;

   glClearColor(red,green,blue,1.0);
   g.background_colour[0] = red; 
   g.background_colour[1] = green; 
   g.background_colour[2] = blue; 
   glFogfv(GL_FOG_COLOR, g.background_colour);
   graphics_draw();

}

int  background_is_black_p() {

   int v = 0;
   graphics_info_t g;
   if (g.background_colour[0] < 0.1)
      if (g.background_colour[1] < 0.1)
	 if (g.background_colour[2] < 0.1)
	    v = 1;

   return v;
}


/*  ----------------------------------------------------------------------- */
/*                  pepflip                                                 */
/*  ----------------------------------------------------------------------- */
// use the values that are in graphics_info
void do_pepflip(short int state) {

   graphics_info_t g;

   g.set_in_pepflip_define(state);
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
      std::cout << "click on a atom in the peptide you wish to flip: "
		<< std::endl;
   } else {
      g.normal_cursor();
   } 
      
} 

void pepflip(int ires, const char *chain_id, int imol) { /* the residue with CO,
							   for scripting interface. */

   if (imol < graphics_n_molecules()) { 
      graphics_info_t g;
      g.molecules[imol].pepflip_residue(ires, std::string(""), std::string(chain_id));
      graphics_draw();
   } 
} 


// ------------------------------------------------------------------
//                                Utility
// ------------------------------------------------------------------
// 
// File system Utility function: maybe there is a better place for it...
// Return like mkdir: mkdir returns zero on success, or -1 if an  error  occurred
//
// if it already exists as a dir, return 0 of course.
// 
int
make_directory_maybe(const char *dir) {
   return coot::util::create_directory(std::string(dir));
} 


void add_coordinates_glob_extension(const char *ext) { 
   
   graphics_info_t g;
   g.add_coordinates_glob_extension(std::string(ext));
} 

void add_data_glob_extension(const char *ext) { 
   graphics_info_t g;
   g.add_data_glob_extension(std::string(ext));
} 

void add_dictionary_glob_extension(const char *ext) { 
   graphics_info_t g;
   g.add_dictionary_glob_extension(std::string(ext));
} 

void add_map_glob_extension(const char *ext) { 
   graphics_info_t g;
   g.add_map_glob_extension(std::string(ext));
} 

int do_anti_aliasing_state() {
   return graphics_info_t::do_anti_aliasing_flag;
} 


void set_do_anti_aliasing(int state) {

   graphics_info_t g;
   g.set_do_anti_aliasing(state);
} 


void set_do_GL_lighting(int state) {
   graphics_info_t::do_lighting_flag = state;
   setup_lighting(state);
   graphics_draw();
}


int do_GL_lighting_state() {
   return graphics_info_t::do_lighting_flag;
}


// Glaxo people: Chuang Chang (surely not how you spell it), Moira X.
// 

/*  ----------------------------------------------------------------------- */
/*                  crosshairs                                              */
/*  ----------------------------------------------------------------------- */
void set_draw_crosshairs(short int i) { 

   graphics_info_t g;
   g.draw_crosshairs_flag = i;
   if (i > 0 ) { 
      g.crosshairs_text(); 
      graphics_draw();
   }
} 

short int draw_crosshairs_state() {
   return graphics_info_t::draw_crosshairs_flag; 
}


/*  ----------------------------------------------------------------------- */
/*                  citation notice                                         */
/*  ----------------------------------------------------------------------- */
void citation_notice_off() { 

   graphics_info_t::show_citation_notice = 0;

} 

/*  ----------------------------------------------------------------------- */
/*                  cursor function                                         */
/*  ----------------------------------------------------------------------- */
void normal_cursor() {

   graphics_info_t g;
   g.normal_cursor();
   graphics_draw();
}

void fleur_cursor() {
   graphics_info_t g;
   g.fleur_cursor();
   graphics_draw();

}

void pick_cursor_maybe() {

   graphics_info_t g;
   g.pick_cursor_maybe();
   graphics_draw();
}

void rotate_cursor() {
   normal_cursor();
}

void set_pick_cursor_index(int i) {
   graphics_info_t::pick_cursor_index = GdkCursorType(i);
}


/*  ------------------------------------------------------------------------ */
/*                       povray/raster3d interface                           */
/*  ------------------------------------------------------------------------ */
void raster3d(const char *filename) {

   graphics_info_t g;
   g.raster3d(std::string(filename));
}

void povray(const char *filename) {

   graphics_info_t g;
   g.povray(std::string(filename));
}

void set_raster3d_bond_thickness(float f) { 

   graphics_info_t::raster3d_bond_thickness = f;

} 

void set_raster3d_density_thickness(float f) {

   graphics_info_t::raster3d_density_thickness = f;

} 

void
raster_screen_shot() {  // run raster3d or povray and guile
                		         // script to render and display image

   // do some checking for povray/render here:

   std::string cmd("(render-image)");  // this is a render function 
   
   // cmd = "(povray-image)";

   safe_scheme_command(cmd);
}


void set_renderer_show_atoms(int istate) {

   graphics_info_t::renderer_show_atoms_flag = istate;
} 


/*  ----------------------------------------------------------------------- */
/*                  browser url                                          */
/*  ----------------------------------------------------------------------- */
void browser_url(const char *url) {

   if (url) { 
      std::string u(url);
      std::vector<std::string> commands;
      commands.push_back("system");
      std::string s = graphics_info_t::browser_open_command;
      if (s == "firefox" || s == "mozilla" || s == "netscape") { 
	 s += " -remote 'openURL(\\\"";
	 s += u;
	 s += "\\\",new-window)'";
	 commands.push_back(single_quote(s));
      } else {
	 if (s == "open") {
	    s += " ";
	    s += url;
	 } else {
	    s += " ";
	    s += url;
	 }
	 commands.push_back(single_quote(s));
      }

      std::string c = languagize_command(commands);
#ifdef USE_GUILE
      safe_scheme_command(c);
#else
#ifdef USE_PYTHON
      safe_python_command(c);
#endif
#endif
   }
}


void set_browser_interface(const char *browser) {

   if (browser) {
      graphics_info_t::browser_open_command = browser;
   } 
} 

void handle_online_coot_search_request(const char *entry_text) {

   if (entry_text) {
      clipper::String text(entry_text);
      std::vector<clipper::String> bits = text.split(" ");
      if (bits.size() > 0) { 
	 std::string s = "http://www.google.co.uk/search?q=";
	 s += bits[0];
	 for (int i=1; i<bits.size(); i++) {
	    s += "+";
	    s += bits[i];
	 }
	 s += "+coot+site%3Awww.ysbl.york.ac.uk";
	 browser_url(s.c_str());
      }
   } 
} 


/*  ----------------------------------------------------------------------- */
/*                  remote control                                          */
/*  ----------------------------------------------------------------------- */
/* section Remote Control */

// called by c_inner_main() if we have guile
void make_socket_listener_maybe() {


   std::vector<std::string> cmd;

   if (graphics_info_t::try_port_listener) { 
      cmd.push_back("open-coot-listener-socket");
      cmd.push_back(graphics_info_t::int_to_string(graphics_info_t::remote_control_port_number));
      cmd.push_back(single_quote(graphics_info_t::remote_control_hostname));

      // This is not a static function, perhaps it should be:
      graphics_info_t g;
      std::string scm_command = g.state_command(cmd, coot::STATE_SCM);

      safe_scheme_command(scm_command);

      if (graphics_info_t::coot_socket_listener_idle_function_token == -1)
	 if (graphics_info_t::listener_socket_have_good_socket_state) 
	    graphics_info_t::coot_socket_listener_idle_function_token =
	       gtk_idle_add((GtkFunction) coot_socket_listener_idle_func,
			    graphics_info_t::glarea);
   }
}

void set_coot_listener_socket_state_internal(int sock_state) {
   graphics_info_t::listener_socket_have_good_socket_state = sock_state;
}



int coot_socket_listener_idle_func(GtkWidget *w) { 

#ifdef USE_GUILE
   std::cout << "DEBUG:: running socket idle function" << std::endl;
   if (graphics_info_t::listener_socket_have_good_socket_state) { 
      std::cout << "DEBUG:: running guile function" << std::endl;
      safe_scheme_command("(coot-listener-idle-function-proc)");
   }
#endif
   return 1;
}

/*  ----------------------------------------------------------------------- */
/*                  Surfaces                                                */
/*  ----------------------------------------------------------------------- */
void do_surface(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].make_surface(state);
      graphics_draw();
   }
}


/*  ----------------------------------------------------------------------- */
/*                  SHELX stuff                                             */
/*  ----------------------------------------------------------------------- */

/* section SHELXL Functions */
// return 
int read_shelx_ins_file(const char *filename) {

   int istat = -1;
   graphics_info_t g;
   if (filename) { 
      int imol = graphics_info_t::n_molecules;
      g.expand_molecule_space_maybe();

      istat = g.molecules[imol].read_shelx_ins_file(std::string(filename));
      if (istat != 1) {
	 std::cout << "ERROR:: " << istat << " on read_shelx_ins_file "
		   << filename << std::endl;
      } else {
	 std::cout << "Molecule " << g.n_molecules << " read successfully\n";
	 istat = g.n_molecules; // for return status 
	 g.n_molecules++;
	 if (g.go_to_atom_window) {
	    g.set_go_to_atom_molecule(imol);
	    g.update_go_to_atom_window_on_new_mol();
	 }
	 graphics_draw();
	 std::vector<std::string> command_strings;
	 command_strings.push_back("read-shelx-ins-file");
	 command_strings.push_back(single_quote(filename));
	 add_to_history(command_strings);
      }
   } else {
      std::cout << "ERROR:: null filename in read_shelx_ins_file" << std::endl;
   }
   return istat;
   
}

int write_shelx_ins_file(int imol, const char *filename) {

   int istat = 0;
   if (filename) { 
      if (is_valid_model_molecule(imol)) {
	 std::pair<int, std::string> stat = graphics_info_t::molecules[imol].write_shelx_ins_file(std::string(filename));
	 istat = stat.first;
	 graphics_info_t g;
	 g.statusbar_text(stat.second);
      } else {
	 std::cout << "WARNING:: invalid molecule (" << imol
		   << ") for write_shelx_ins_file" << std::endl;
      }
   }
   return istat;
}


/*  ----------------------------------------------------------------------- */
/*                  SMILES                                                  */
/*  ----------------------------------------------------------------------- */
void do_smiles_gui() {

#ifdef USE_GUILE

   safe_scheme_command("(smiles-gui)");


#endif // USE_GUILE

} 

