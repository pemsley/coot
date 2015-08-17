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

#include "compat/coot-sysdep.h"

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

#include <mmdb2/mmdb_manager.h>
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
#include "cc-interface-scripting.hh"

#include "ligand/ligand.hh" // for rigid body fit by atom selection.

#include "coot-fileselections.h"
#include "cmtz-interface.hh" // for valid columns mtz_column_types_info_t
#include "c-interface-mmdb.hh"
#include "c-interface-scm.hh"
#include "c-interface-python.hh"

/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */

void execute_refmac(GtkWidget *window) { 

   // The passed window, is the refmac dialog, where one selects the
   // coords molecule and the map molecule.

   graphics_info_t g;
   GtkWidget *option_menu = lookup_widget(window, "run_refmac_coords_optionmenu");
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));

   int imol_coords = graphics_info_t::refmac_molecule;

   std::cout << "------------------ in execute_refmac: imol_coords is " << imol_coords << std::endl;
   
   if (! is_valid_model_molecule(imol_coords)) {

      std::cout << "INFO:: No coordinates molecule selected for running refmac\n";

   } else {

      option_menu = lookup_widget(window, "run_refmac_map_optionmenu");
      GtkWidget *mtz_file_radiobutton = lookup_widget(window, "run_refmac_mtz_file_radiobutton");
      menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
      active_item = gtk_menu_get_active(GTK_MENU(menu));
      int imol_map_refmac = -1;
      bool have_mtz_file = false;
      if (mtz_file_radiobutton && GTK_TOGGLE_BUTTON(mtz_file_radiobutton)->active) {
	 have_mtz_file = true;
      }

      // active_item is set if there was at least one map with refmac params:
      // if none, it is null.
      
      if ((active_item == 0) && (have_mtz_file == false)) {
	 add_status_bar_text("No map has associated Refmac Parameters - no REFMAC!");
      } else {
	 std::cout << "DEBUG:: Happy-path 1 active_item: " << active_item
		   << " and have_mtz_file " << have_mtz_file << std::endl;
	 int imol_window = -1;
	 std::string mtz_in_filename = "";
	 if (! have_mtz_file) {
	    // we get imol from a map mtz file
	    // get imol_window from active item of run_refmac_map_optionmenu.
	    //
	    imol_window = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(active_item)));
	 } else {
	    // check the filename of the button
	    GtkWidget *button_mtz_label = lookup_widget(window, "run_refmac_mtz_file_label");
	    const gchar *mtz_filename = gtk_label_get_text(GTK_LABEL(button_mtz_label));
	    std::cout << "DEBUG:: set mtz_in_filename to " << mtz_in_filename << std::endl;
	    mtz_in_filename = mtz_filename;
	    if (mtz_in_filename == "(None)") {
	       have_mtz_file = false;
	       std::cout << "WARNING:: [A] no mtz file given" <<std::endl;
	    } else {
	       if (! coot::file_exists(mtz_in_filename)) {
		  have_mtz_file = false;
		  std::cout << "WARNING:: mtz file " << mtz_in_filename << " does not exist" <<std::endl;
	       }
	    }
	 }

	 std::cout << "DEBUG:: imol_window: " << imol_window << std::endl;

	 if ((imol_window < 0) && (have_mtz_file == false)) { // do you mean || here?
	    if (have_mtz_file == false) {
	       std::cout << "INFO:: No mtz file selected for refmac\n";
	    } else {
	       std::cout << "INFO:: No map data selected for refmac\n";
	    }
	 } else {
	    std::cout << "Happy-path-2 with imol_window " << imol_window << std::endl;
	    imol_map_refmac = imol_window;
	    if (! is_valid_map_molecule(imol_map_refmac) && (! have_mtz_file)) {
	       std::string s = "Invalid molecule number: ";
	       s += graphics_info_t::int_to_string(imol_map_refmac);
	       std::cout << s << std::endl;
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
		     if (! have_ccp4i_project) { 
			pdb_in_filename  += "/";
			pdb_out_filename += "/";
			mtz_out_filename += "/";
		     }
		     pdb_in_filename += g.molecules[imol_coords].Refmac_in_name();

		     // cleverness happens in Refmac_out_name:
		     pdb_out_filename += g.molecules[imol_coords].Refmac_out_name();
		     mtz_out_filename += g.molecules[imol_coords].Refmac_mtz_out_name();

		     if (! have_mtz_file) {
			mtz_in_filename = g.molecules[imol_map_refmac].Refmac_mtz_filename();
		     }
		     std::string refmac_count_string =
			g.int_to_string(g.molecules[imol_coords].Refmac_count());

		     std::cout << "DEBUG:: mtz_out_filename: " << mtz_out_filename << std::endl;
		     std::cout << "DEBUG:: pdb_out_filename: " << pdb_out_filename << std::endl;

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

		     // and we need to get the column labels if we have an input mtz (not map)
		     if (have_mtz_file) {

			if (refmac_use_twin_state() == 0) {
			   // for now we only use Is in twin not in 'normal' refinement
			   icol = saved_f_phi_columns->selected_refmac_fobs_col;
			   if ((icol >=0 ) && (icol < int(saved_f_phi_columns->f_cols.size()))) { 
			      fobs_col = saved_f_phi_columns->f_cols[icol].column_label;  // Minmin crash
			   } else {
			      std::cout << "ERROR:: trapped inappropriate access of f_cols in execute_refmac() "
					<< icol << " " << saved_f_phi_columns->f_cols.size() << std::endl;
			   } 
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
   g.set_refmac_phase_input(phase_flag);

}

void set_refmac_use_sad(int state) {

   graphics_info_t g;
   g.set_refmac_use_sad(state);

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

void
wrapped_create_run_refmac_dialog() {
   
   GtkWidget *window = create_run_refmac_dialog();
   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(refmac_molecule_button_select);
   GtkWidget *diff_map_button = lookup_widget(window, "run_refmac_diff_map_checkbutton");
   GtkWidget *optionmenu;
   int imol_coords = first_coords_imol();
   int have_file = 0;

   GtkWidget *labels = lookup_widget(window, "run_refmac_column_labels_frame");
   GtkWidget *ncs_button = lookup_widget(window, "run_refmac_ncs_checkbutton");
   GtkWidget *mtz_file_radiobutton = lookup_widget(window, "run_refmac_mtz_file_radiobutton");
   optionmenu = lookup_widget(window, "run_refmac_method_optionmenu");
   fill_option_menu_with_refmac_methods_options(optionmenu);

   optionmenu = lookup_widget(window, "run_refmac_phase_input_optionmenu");
   fill_option_menu_with_refmac_phase_input_options(optionmenu);
   if (GTK_TOGGLE_BUTTON(mtz_file_radiobutton)->active) have_file = 1;

   set_refmac_molecule(imol_coords);

   optionmenu = lookup_widget(window, "run_refmac_coords_optionmenu");
   fill_option_menu_with_coordinates_options(optionmenu, callback_func, imol_coords);

   optionmenu = lookup_widget(window, "run_refmac_map_optionmenu");
   
   /*  fill_option_menu_with_refmac_options(optionmenu); */
   fill_option_menu_with_refmac_labels_options(optionmenu); // change the name of this function -
                                                            // they are molecules with refmac mtz files.

   /* to set the labels set the active item; only if not twin and
      if we really want the labels from map mtz*/
   if (refmac_use_twin_state() == 0 && have_file == 0) {
      GtkWidget *active_menu_item =
	 gtk_menu_get_active(GTK_MENU(gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu))));
      if (active_menu_item) {
	 gtk_menu_item_activate(GTK_MENU_ITEM(active_menu_item));
      }
   }

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

   /* set the ncs button depending on state */
   if (refmac_use_ncs_state()) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ncs_button), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ncs_button), FALSE);
   }

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
