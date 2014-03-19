/* src/c-interface-ligands.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2008, 2009 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <stdlib.h>
#include <iostream>
#include <stdexcept>

#if defined _MSC_VER
#include <windows.h>
#endif

#include "lbg/lbg.hh" // it matters where this is - rdkit issues...

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

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/peak-search.hh"

#include "ligand/wligand.hh"

#include "guile-fixups.h"

/* in here we check if libcheck is available (if scripting is available) */
GtkWidget *wrapped_create_libcheck_monomer_dialog() {

   GtkWidget *w = create_libcheck_monomer_dialog();

#ifdef USE_GUILE

   std::string c = "(command-in-path? libcheck-exe)";
   SCM v = safe_scheme_command(c.c_str());

   if (!scm_is_true(v)) {
      GtkWidget *l = lookup_widget(w, "no_libcheck_frame");
      if (l) {
	 gtk_widget_show(l);
      }
   } 
   
#endif // USE_GUILE

#ifdef USE_PYTHON

   // something similar here.
   
#endif // USE_PYTHON   
   
   return w;
} 

/*  ----------------------------------------------------------------------- */
/*                  ligand fitting stuff                                    */
/*  ----------------------------------------------------------------------- */

// We need to look up the vboxes and add items, like we did in code in
// gtk-manual.c
// 
// This is called by callback.c's find ligand button clicked callback.
// 
// Return the state.  If there is not ligands, protein or map, then
// return 0.  This is checked by the calling function (in callbacks.c)
// and the dialog is not displayed if 0 is returned.
// 
int fill_ligands_dialog(GtkWidget *find_ligand_dialog) {

   int ifound_map, ifound_coords, ifound_ligand; 
   short int diff_maps_only_flag = 0;
   ifound_map = fill_ligands_dialog_map_bits(find_ligand_dialog, diff_maps_only_flag);
   if (ifound_map == 0) {
      std::cout << "WARNING:: you must have a map to search for ligands!"
		<< std::endl;
      std::string s("WARNING:: you must have a map to\n search for ligands!");
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   } 
   ifound_coords = fill_ligands_dialog_protein_bits(find_ligand_dialog);
   if (ifound_coords == 0) {
      std::cout << "Error: you must have a protein to mask the map!"
		<< std::endl;
      std::string s("WARNING:: you must have a protein\n to mask the map");
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   } 
   ifound_ligand = fill_ligands_dialog_ligands_bits(find_ligand_dialog);
   if (ifound_ligand == 0) {
      std::cout << "Error: you must have at least one  ligand to search for!"
		<< std::endl;
      std::string s("WARNING:: you must have at least one\n          ligand to search for!\n");
      s += "         Ligands have less than ";
      s += graphics_info_t::int_to_string(graphics_info_t::find_ligand_ligand_atom_limit);
      s += " atoms\n";
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   }

   // The mask waters toggle buttons:

   GtkWidget *togglebutton;
   togglebutton = lookup_widget(find_ligand_dialog,
				"find_ligand_mask_waters_yes_radiobutton");
   if (graphics_info_t::find_ligand_mask_waters_flag)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(togglebutton), TRUE);
   else 
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(togglebutton), FALSE);


   // The Search/Here toggle buttons:
   // 
   GtkWidget *search_here_toggle_button;
   search_here_toggle_button = lookup_widget(find_ligand_dialog,
					     "find_ligands_search_here_radiobutton");
   if (search_here_toggle_button) 
      if (graphics_info_t::find_ligand_here_cluster_flag)
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(search_here_toggle_button), TRUE);
   
   // 040211: fill new sigma level entry
   fill_ligands_sigma_level_entry(find_ligand_dialog);

   // expert options
   fill_ligands_expert_options(find_ligand_dialog);
   // shall we see the expert option frame?
   if (graphics_info_t::ligand_expert_flag == 0) {
      GtkWidget *frame = lookup_widget(find_ligand_dialog, "ligand_expert_frame");
      gtk_widget_hide(frame);
   }
      
   return ifound_ligand * ifound_map * ifound_coords;

   // 050924 New "Expert Options" entries:
   

}

// 110204: fill new sigma level entry
// 
void fill_ligands_sigma_level_entry(GtkWidget *dialog) {
   
   GtkWidget *entry = lookup_widget(dialog, "find_ligand_sigma_level_entry");
   gtk_entry_set_text(GTK_ENTRY(entry), graphics_info_t::float_to_string(graphics_info_t::ligand_cluster_sigma_level).c_str());
   
}

void fill_ligands_expert_options(GtkWidget *find_ligand_dialog) {

   GtkWidget *entry = lookup_widget(find_ligand_dialog, "ligand_n_samples_entry");
   gtk_entry_set_text(GTK_ENTRY(entry),
		      graphics_info_t::int_to_string(graphics_info_t::ligand_wiggly_ligand_n_samples).c_str());

   entry = lookup_widget(find_ligand_dialog, "ligand_n_top_ligands_entry");
   gtk_entry_set_text(GTK_ENTRY(entry),
		      graphics_info_t::int_to_string(graphics_info_t::find_ligand_n_top_ligands).c_str());

}


int fill_ligands_dialog_map_bits(GtkWidget *find_ligand_dialog, 
				 short int diff_maps_only_flag) {

   return fill_ligands_dialog_map_bits_by_dialog_name(find_ligand_dialog,
						      "find_ligand_map",
						      diff_maps_only_flag);

}

// dialog_name is typically find_ligand_map
// 
int fill_ligands_dialog_map_bits_by_dialog_name(GtkWidget *find_ligand_dialog,
						const char *dialog_name, 
						short int diff_maps_only_flag) {

   int ifound = 0;
   graphics_info_t g;
   // Add map elements:
   GSList *find_ligand_map_group = NULL;
   //
   std::string vbox_name = dialog_name;
   vbox_name += "_vbox";
   
   GtkWidget *find_ligand_map_vbox =
      lookup_widget(find_ligand_dialog,vbox_name.c_str());
   if (find_ligand_map_vbox == NULL) {
      std::cout << "disaster! find_ligand map vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules(); imol++) {
	 if (g.molecules[imol].has_xmap()) {

	    if ((!diff_maps_only_flag) ||
		g.molecules[imol].is_difference_map_p()) {

	       ifound++; // there was a map
	       std::string map_str(dialog_name);
	       map_str += "_radiobutton_";
	       map_str += g.int_to_string(imol);
	       std::string map_button_label = g.int_to_string(imol);
	       map_button_label += " ";
	       map_button_label += g.molecules[imol].name_;

	       GtkWidget *find_ligand_map_radiobutton_imol =
		  gtk_radio_button_new_with_label (find_ligand_map_group,
						   map_button_label.c_str());
	       find_ligand_map_group =
		  gtk_radio_button_group (GTK_RADIO_BUTTON (find_ligand_map_radiobutton_imol));
	       gtk_widget_ref (find_ligand_map_radiobutton_imol);
	       gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
					 map_str.c_str(),
					 find_ligand_map_radiobutton_imol,
					 (GtkDestroyNotify) gtk_widget_unref);
	       gtk_widget_show (find_ligand_map_radiobutton_imol);
	       gtk_box_pack_start (GTK_BOX (find_ligand_map_vbox),
				   find_ligand_map_radiobutton_imol, FALSE, FALSE, 0);
	    }
	 }
      }
   }
   return ifound; 
}

int fill_ligands_dialog_protein_bits(GtkWidget *find_ligand_dialog) {
   
   return fill_ligands_dialog_protein_bits_by_dialog_name(find_ligand_dialog,
							  "find_ligand_protein");
   
}


int fill_ligands_dialog_protein_bits_by_dialog_name(GtkWidget *find_ligand_dialog,
						    const char *dialog_name) {
   
   
   int ifound = 0;
   graphics_info_t g;
   // Add protein elements:
   GSList *find_ligand_protein_group = NULL;
   //
   std::string vbox_name(dialog_name);
   vbox_name += "_vbox";
   
   GtkWidget *find_ligand_protein_vbox =
      lookup_widget(find_ligand_dialog,vbox_name.c_str());
   if (find_ligand_protein_vbox == NULL) {
      std::cout << "disaster! find_ligand protein vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules(); imol++) {
	 if (g.molecules[imol].atom_sel.n_selected_atoms) {
	    //
	    ifound = 1; // there was a protein
	    std::string protein_str(dialog_name);
	    protein_str += "_radiobutton_";
	    protein_str += g.int_to_string(imol);
	    std::string protein_button_label = g.int_to_string(imol);
	    protein_button_label += " ";
	    protein_button_label += g.molecules[imol].name_;
	    GtkWidget *find_ligand_protein_radiobutton_imol =
	       gtk_radio_button_new_with_label (find_ligand_protein_group,
						protein_button_label.c_str());
	    find_ligand_protein_group =
	       gtk_radio_button_group (GTK_RADIO_BUTTON
				       (find_ligand_protein_radiobutton_imol));
	    gtk_widget_ref (find_ligand_protein_radiobutton_imol);
	    gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
				      protein_str.c_str(),
				      find_ligand_protein_radiobutton_imol,
				      (GtkDestroyNotify) gtk_widget_unref);
	    gtk_widget_show (find_ligand_protein_radiobutton_imol);
	    gtk_box_pack_start (GTK_BOX (find_ligand_protein_vbox),
				find_ligand_protein_radiobutton_imol, FALSE, FALSE, 0);
	 }
      }
   }
   return ifound; 
}

int fill_vbox_with_coords_options_by_dialog_name(GtkWidget *find_ligand_dialog,
						 const char *dialog_name,
						 short int have_ncs_flag) {
   
   
   int ifound = 0;
   graphics_info_t g;
   // Add protein elements:
   GSList *find_ligand_protein_group = NULL;
   //
   std::string vbox_name(dialog_name);
   vbox_name += "_vbox";
   
   GtkWidget *find_ligand_protein_vbox =
      lookup_widget(find_ligand_dialog,vbox_name.c_str());
   if (find_ligand_protein_vbox == NULL) {
      std::cout << "disaster! fill_vbox_with_coords_options_by_dialog_name coords"
		<< " vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules(); imol++) {
	 if (g.molecules[imol].has_model()) {

	    if (!have_ncs_flag || g.molecules[imol].has_ncs_p()) { 
	       //
	       ifound = 1; // there was a protein
	       std::string protein_str(dialog_name);
	       protein_str += "_radiobutton_";
	       protein_str += g.int_to_string(imol);
	       std::string protein_button_label = g.int_to_string(imol);
	       protein_button_label += " ";
	       protein_button_label += g.molecules[imol].name_;
	       GtkWidget *find_ligand_protein_radiobutton_imol =
		  gtk_radio_button_new_with_label (find_ligand_protein_group,
						   protein_button_label.c_str());
	       find_ligand_protein_group =
		  gtk_radio_button_group (GTK_RADIO_BUTTON
					  (find_ligand_protein_radiobutton_imol));
	       gtk_widget_ref (find_ligand_protein_radiobutton_imol);
	       gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
					 protein_str.c_str(),
					 find_ligand_protein_radiobutton_imol,
					 (GtkDestroyNotify) gtk_widget_unref);
	       gtk_widget_show (find_ligand_protein_radiobutton_imol);
	       gtk_box_pack_start (GTK_BOX (find_ligand_protein_vbox),
				   find_ligand_protein_radiobutton_imol, FALSE, FALSE, 0);
	    }
	 }
      }
   }
   return ifound; 
}


int fill_ligands_dialog_ligands_bits(GtkWidget *find_ligand_dialog) {

   int ifound = 0;
   graphics_info_t g;
   // Add ligand elements
   //
   GtkWidget *hbox;

   GtkWidget *find_ligand_ligands_vbox =
      lookup_widget(find_ligand_dialog, "find_ligand_ligands_vbox");
   if (find_ligand_ligands_vbox == NULL) {
      std::cout << "disaster! find_ligand protein vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules(); imol++) {
	 if ((g.molecules[imol].atom_sel.n_selected_atoms < graphics_info_t::find_ligand_ligand_atom_limit) &&
	     g.molecules[imol].has_model()) {
	    ifound = 1; // there was a ligand

	    // create an hbox:
	    hbox = gtk_hbox_new (FALSE, 0);
	    
	    std::string ligands_str("find_ligand_ligand_checkbutton_");
	    ligands_str += g.int_to_string(imol);
	    std::string ligands_button_label = g.int_to_string(imol);
	    ligands_button_label += " ";
	    ligands_button_label += g.molecules[imol].name_;

	    // ligand on/off check button
	    // 
	    GtkWidget *find_ligand_ligands_checkbutton_imol =
	       gtk_check_button_new_with_label (ligands_button_label.c_str());
	    gtk_widget_ref (find_ligand_ligands_checkbutton_imol);
	    gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
				      ligands_str.c_str(),
				      find_ligand_ligands_checkbutton_imol,
				      (GtkDestroyNotify) gtk_widget_unref);

	    gtk_widget_show (find_ligand_ligands_checkbutton_imol);
	    gtk_box_pack_start (GTK_BOX (hbox),
				find_ligand_ligands_checkbutton_imol, FALSE, FALSE, 0);

	    // wligand on/off check button
	    //
	    std::string wligands_str("find_ligand_wligand_checkbutton_");
	    wligands_str += g.int_to_string(imol);
	    GtkWidget *find_ligand_wligands_checkbutton_imol =
	       gtk_check_button_new_with_label ("Flexible?");
	    gtk_widget_ref (find_ligand_wligands_checkbutton_imol);
	    gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
				      wligands_str.c_str(),
				      find_ligand_wligands_checkbutton_imol,
				      (GtkDestroyNotify) gtk_widget_unref);

	    gtk_widget_show (find_ligand_wligands_checkbutton_imol);
	    gtk_box_pack_start (GTK_BOX (hbox),
				find_ligand_wligands_checkbutton_imol, FALSE, FALSE, 0);

	    // pack the hbox into the vbox
	    gtk_box_pack_start (GTK_BOX (find_ligand_ligands_vbox),
				hbox, FALSE, FALSE, 0);
	    gtk_widget_show(hbox);
	 }
      }
   }
   return ifound; 
}


/* get which map to search, protein mask and ligands from button and
   then do it*/
// Recall that the map and the protein are radio buttons (i.e. we want
// exactly one of each).  Hoewever, the ligands are check buttons, we
// can have as many of those on that we like.
//
// We should add a check to the install_ligand that checks the number
// of atoms in the ligand and gives a warning... hmmm... that invovles
// another popup-sigh...  Oh well, plan for it, even if you don't
// implement it now - i.e. invoke a post-check function that creates
// the coot::ligand object.
//
void execute_get_mols_ligand_search(GtkWidget *button) {

   graphics_info_t g;
   GtkWidget *ligand_button;
   std::vector<int> chief_ligand_many_atoms; // caches imols with lots
					     // of atoms.
   std::vector<std::pair<int, bool> > wiggly_ligand_info; 

   // extract the sigma level and stick it in
   // graphics_info_t::ligand_cluster_sigma_level
   // 
   set_ligand_cluster_sigma_level_from_widget(button);
   set_ligand_expert_options_from_widget(button); // ligand_n_top_sites
						  // and wiggly ligand
						  // n (conformer)
						  // samples
   
   // flags to say if we have found things:
   short int found_active_button_for_map = 0;
   short int found_active_button_for_protein = 0;
   short int found_active_button_for_ligands = 0;
   // This is where we store the mols:
   int find_ligand_map_mol = -1; // gets set?
   int find_ligand_protein_mol = -1; // gets set.
   // std::vector<int> find_ligand_ligand_mols;  old.
   
   // Find the first active map radiobutton
   found_active_button_for_map = 0;
   for (int imol=0; imol<g.n_molecules(); imol++) {
      if (g.molecules[imol].has_xmap()) { 
	 std::string map_str = "find_ligand_map_radiobutton_";
	 map_str += g.int_to_string(imol);
	 ligand_button = lookup_widget(button, map_str.c_str());
	 if (ligand_button) { 
	    if (GTK_TOGGLE_BUTTON(ligand_button)->active) {
	       find_ligand_map_mol = imol;
	       found_active_button_for_map = 1;
	       break;
	    }
	 } else {
	    std::cout << map_str << " widget not found in "
		      << "execute_get_mols_ligand_search" << std::endl;
	 }
      }
   }

   // Find the first active protein radiobutton
   found_active_button_for_protein = 0;
   for (int imol=0; imol<g.n_molecules(); imol++) {
      if (g.molecules[imol].atom_sel.n_selected_atoms > 0) { 
	 std::string protein_str = "find_ligand_protein_radiobutton_";
	 protein_str += g.int_to_string(imol);
	 ligand_button = lookup_widget(button, protein_str.c_str());
	 if (ligand_button) { 
	    if (GTK_TOGGLE_BUTTON(ligand_button)->active) {
	       find_ligand_protein_mol = imol;
	       found_active_button_for_protein = 1;
	       break;
	    }
	 } else {
	    std::cout << protein_str << " widget not found in "
		      << "execute_get_mols_ligand_search" << std::endl;
	 }
      }
   }

   // Now, do we mask waters for the protein mask?
   
   GtkWidget *togglebutton;
   togglebutton = lookup_widget(button,
				"find_ligand_mask_waters_yes_radiobutton");
   if (GTK_TOGGLE_BUTTON(togglebutton)->active)
      graphics_info_t::find_ligand_mask_waters_flag = 1;
   else 
      graphics_info_t::find_ligand_mask_waters_flag = 0;


   // The Search/Here toggle buttons:
   
   GtkWidget *search_here_toggle_button;
   search_here_toggle_button = lookup_widget(button,
					     "find_ligands_search_here_radiobutton");
   if (search_here_toggle_button) {
      if (GTK_TOGGLE_BUTTON(search_here_toggle_button)->active) {
	 std::cout << " Activating SEARCH HERE in ligand fitting" << std::endl;
	 graphics_info_t::find_ligand_here_cluster_flag = 1;
      } else {
	 graphics_info_t::find_ligand_here_cluster_flag = 0;
      }
   }
   

   // For each imol in molecules, if we have selected coordinates,
   // construct a string begining "find_ligand_ligands_imol_" then the
   // molecule number
   
   found_active_button_for_ligands = 0;
   for (int imol=0; imol<g.n_molecules(); imol++) {
      if (g.molecules[imol].has_model() &&
	  g.molecules[imol].atom_sel.n_selected_atoms < graphics_info_t::find_ligand_ligand_atom_limit) {
	 std::string ligand_str = "find_ligand_ligand_checkbutton_";
	 ligand_str += g.int_to_string(imol);
	 ligand_button = lookup_widget(button, ligand_str.c_str());

	 std::string wiggly_str = "find_ligand_wligand_checkbutton_";
	 wiggly_str += g.int_to_string(imol);
	 GtkWidget *wiggly_button = lookup_widget(button, wiggly_str.c_str());
		  
	 if (ligand_button && wiggly_button) { 
	    if (GTK_TOGGLE_BUTTON(ligand_button)->active) {
	       
	       bool wiggly_state = 0;
	       if (GTK_TOGGLE_BUTTON(wiggly_button)->active)
		  wiggly_state = 1;
	       wiggly_ligand_info.push_back(std::pair<int, bool> (imol, wiggly_state));
// 	       std::cout << "DEBUG:: wiggly info: " << imol <<  " " << wiggly_state
// 			 << " pushed back" << std::endl;
	       found_active_button_for_ligands = 1;
	       int n_atoms = g.molecules[imol].atom_sel.n_selected_atoms;
	       if (n_atoms > 100) {
		  std::cout << "WARNING:: molecule " << imol
			    << " has unexpectedly many atoms ("
			    << n_atoms << ")" << std::endl;
		  chief_ligand_many_atoms.push_back(imol);
	       } 
// 	    } else {
// 	       std::cout << "DEBUG:: button " << ligand_str
// 			 << " was not active" << std::endl;
	    }
	 } else {
	    std::cout << ligand_str << " widget not found in "
		      << "execute_get_mols_ligand_search" << std::endl;
	 }
      }
   }

   if ( found_active_button_for_map &&
	found_active_button_for_protein &&
	found_active_button_for_ligands) {

      // So store the clicked info in a graphics_info_t static:
      // because we need this info when we click the OK button of the
      // create_find_ligand_many_atoms_dialog() widget.  We don't want
      // to mess with set_user_data for many data.
      // 
      graphics_info_t g;
      g.set_find_ligands_mols(find_ligand_map_mol,
			      find_ligand_protein_mol,
			      wiggly_ligand_info);

      if (chief_ligand_many_atoms.size() == 0 ) {
#if defined USE_GUILE && !defined WINDOWS_MINGW
	 execute_ligand_search();
#else
#ifdef USE_PYTHON
	 execute_ligand_search_py();
#endif // USE_PYTHON
#endif // USE_GUILE
      } else {
	 // we need to delete this widget when OK and cancel of the
	 // many atoms widget is pressed, so let's attach it as user
	 // data.
	 GtkWidget *widget = lookup_widget(button, "find_ligand_dialog");
	 do_find_ligand_many_atoms_in_ligands(widget);
      }
   } else {
	 std::cout << "Something wrong in the selection of map/molecules"
		   << std::endl;
   }

}

// q_ligands is questionable ligands (i.e. very large)
//   
void do_find_ligand_many_atoms_in_ligands(GtkWidget *find_ligand_dialog) {

   GtkWidget *widget = create_find_ligand_many_atoms_dialog();
   gtk_object_set_user_data(GTK_OBJECT(widget), (char *) find_ligand_dialog);
   gtk_widget_show(widget); 
}

void set_ligand_expert_options_from_widget(GtkWidget *button) {

   GtkWidget *entry = lookup_widget(button, "ligand_n_samples_entry");
   const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
   if (text) {
      int isample = atoi(text);
      if ((isample > 0) && (isample < 1000000))
	 graphics_info_t::ligand_wiggly_ligand_n_samples = isample;
   }
   entry = lookup_widget(button, "ligand_n_top_ligands_entry");
   text = gtk_entry_get_text(GTK_ENTRY(entry));
   if (text) {
      int itop = atoi(text);
      if ((itop > 0) && (itop < 1000000))
	 graphics_info_t::find_ligand_n_top_ligands = itop;
   }
} 


/*  extract the sigma level and stick it in */
/*  graphics_info_t::ligand_cluster_sigma_level */
void set_ligand_cluster_sigma_level_from_widget(GtkWidget *button) { 
   
   GtkWidget *entry = lookup_widget(button, "find_ligand_sigma_level_entry");
   short int setit = 0;

   if (entry) {
      const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
      if (text) { 
	 float f = atof(text);
	 if (f > 0.0 && f < 1000.0) { 
	    graphics_info_t::ligand_cluster_sigma_level = f;
	    setit = 1;
	 } 
      }
   }
   if (setit == 0) 
      std::cout << "INFO:: ignoring bogus attempt to set "
		<< "the ligand search sigma level" 
		<< std::endl;
} 



// The name has beend changed because the function we want at the
// scripting layer (start_ligand_builder_gui()) should not need
// arguments.
// 
void
start_ligand_builder_gui_internal(GtkMenuItem     *menuitem,
				  gpointer         user_data) {

      start_ligand_builder_gui();
}

void
start_ligand_builder_gui() { 

   if (graphics_info_t::use_graphics_interface_flag) { 

#ifdef HAVE_GOOCANVAS
      lig_build::molfile_molecule_t mm;
      CMMDBManager *mol = NULL;
      std::string molecule_file_name = "coot-lidia.mol"; // non-null file name passed to lbg, used
      // in save function
      std::string view_name;
      std::pair<bool, coot::residue_spec_t> dummy_pair(0, coot::residue_spec_t());
      bool use_graphics_interface_flag = 1;
      bool stand_alone_flag = 0;
      int imol_dummy = -1;

      int (*get_url_func_pointer) (const char *s1, const char *s2) = NULL;
#ifdef USE_LIBCURL
      get_url_func_pointer= coot_get_url;
#endif    

      lbg(mm, dummy_pair, mol, view_name, molecule_file_name, imol_dummy,
	  graphics_info_t::Geom_p(),
	  use_graphics_interface_flag, stand_alone_flag,
	  get_url_func_pointer,
	  prodrg_import_function,
	  sbase_import_function,
	  get_drug_mdl_via_wikipedia_and_drugbank
	  );
#else
      std::cout << "No goocanvas" << std::endl;
#endif // HAVE_GOOCANVAS

   }
}
