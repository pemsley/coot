/* src/c-interface-validate.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009, 2010 The University of Oxford
 * Author: Paul Emsley
 * Copyright 2006, 2007 by Bernhard Lohkamp
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
 
#if defined _MSC_VER
#include <windows.h>
#endif

#include <stdlib.h>
#include <iostream>

 
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


#include "graphics-info.h"

#ifdef USE_GUILE
#include <guile/gh.h>

#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    

#endif // SCM version
#endif // USE_GUILE

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include "c-interface.h"
#include "cc-interface.hh"
#include "ligand.hh"

#include "peak-search.hh"
#include "user-mods.hh"

/*  ----------------------------------------------------------------------- */
/*                  check waters interface                                  */
/*  ----------------------------------------------------------------------- */

void set_check_waters_b_factor_limit(float f) {
   graphics_info_t::check_waters_b_factor_limit = f;
}


void set_check_waters_map_sigma_limit(float f) {
   graphics_info_t::check_waters_map_sigma_limit = f;
}


void set_check_waters_min_dist_limit(float f) {
   graphics_info_t::check_waters_min_dist_limit = f;
} 


void set_check_waters_max_dist_limit(float f) {
   graphics_info_t::check_waters_max_dist_limit = f;
}

GtkWidget *wrapped_create_check_waters_dialog() { 

   // Need to fill:
   // 
   // Molecule number optionmenu (check_waters_molecule_optionmenu)
   // 
   // check_waters_b_factor_entry
   // check_waters_map_sigma_entry
   // check_waters_min_dist_entry
   // check_waters_max_dist_entry

   // There is a imol_refinement map check done in the callbacks.c
   // function that calls this function

   GtkWidget *dialog = create_check_waters_dialog();

   // Opps - this (logical OR) should be on by default:
   GtkWidget *check_waters_OR_radiobutton  = lookup_widget(dialog, "check_waters_OR_radiobutton");

   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_waters_OR_radiobutton), TRUE);

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(check_waters_molecule_menu_item_activate);

   GtkWidget *optionmenu = lookup_widget(dialog, "check_waters_molecule_optionmenu");
//    std::cout << "optionmenu: " << optionmenu << std::endl;
//    std::cout << "optionmenu is widget: " << GTK_IS_WIDGET(optionmenu) << std::endl;
//    std::cout << "optionmenu is option menu: " << GTK_IS_OPTION_MENU(optionmenu) << std::endl;

   // now fill that dialog's optionmenu with coordinate options.
   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) { 
	 graphics_info_t::check_waters_molecule = imol;
	 break;
      }
   }
   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(optionmenu, callback_func,
					       graphics_info_t::check_waters_molecule);

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
   
   GtkWidget *entry;
   // char text[100];
   std::string text_str;

   // b-factor
   entry = lookup_widget(dialog, "check_waters_b_factor_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_b_factor_limit);
   gtk_entry_set_text(GTK_ENTRY(entry), text_str.c_str());
      

   // map sigma
   entry = lookup_widget(dialog, "check_waters_map_sigma_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_map_sigma_limit);
   gtk_entry_set_text(GTK_ENTRY(entry), text_str.c_str());

   // min_dist
   entry = lookup_widget(dialog, "check_waters_min_dist_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_min_dist_limit);
   gtk_entry_set_text(GTK_ENTRY(entry), text_str.c_str());

   // max_dist
   entry = lookup_widget(dialog, "check_waters_max_dist_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_max_dist_limit);
   gtk_entry_set_text(GTK_ENTRY(entry), text_str.c_str());

   // 20100131 We have put the variance map check into this dialog
   // too, to better organise the menus.
   //
   GtkWidget *diff_map_option_menu =
      lookup_widget(dialog, "check_water_by_difference_map_optionmenu");

   // It's not in the gtk1 version
   // 
   if (diff_map_option_menu) {
      // set imol_active to the first difference map.
      int imol_active = -1; // unset (no difference maps found yet)
      for (int i=0; i<graphics_n_molecules(); i++) {
	 if (is_valid_map_molecule(i)) {
	    if (map_is_difference_map(i)) {
	       imol_active = i;
	       break;
	    }
	 }
      }

      if (imol_active != -1) {
	 graphics_info_t::check_waters_by_difference_map_map_number = imol_active;
	 GtkSignalFunc signal_func =
	    GTK_SIGNAL_FUNC(check_water_by_difference_maps_option_menu_item_select);
	 g.fill_option_menu_with_difference_map_options(diff_map_option_menu,
							signal_func, imol_active);
      }
   }

   return dialog;

}

void
check_water_by_difference_maps_option_menu_item_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t::check_waters_by_difference_map_map_number = pos;
}

// The OK button was pressed on the dialog, so read the dialog and do
// the check.
//
// called by a callbacks.c function.
// 
void do_check_waters_by_widget(GtkWidget *dialog) {


   // GtkWidget *optionmenu = lookup_widget(dialog, "check_waters_molecule_optionmenu");
   GtkWidget *action_optionmenu = lookup_widget(dialog, "check_waters_action_optionmenu");
   // GtkWidget *checklogic_AND_radiobutton = lookup_widget(dialog, "check_waters_AND_radiobutton");
   GtkWidget *checklogic_OR_radiobutton  = lookup_widget(dialog, "check_waters_OR_radiobutton");

   GtkWidget *entry1, *entry2, *entry3, *entry4;
   entry1 = lookup_widget(dialog, "check_waters_b_factor_entry");
   entry2 = lookup_widget(dialog, "check_waters_map_sigma_entry");
   entry3 = lookup_widget(dialog, "check_waters_min_dist_entry");
   entry4 = lookup_widget(dialog, "check_waters_max_dist_entry");

   //
   float b_factor_lim  = get_positive_float_from_entry(GTK_ENTRY(entry1));
   float map_sigma_lim = get_positive_float_from_entry(GTK_ENTRY(entry2));
   float min_dist      = get_positive_float_from_entry(GTK_ENTRY(entry3));
   float max_dist      = get_positive_float_from_entry(GTK_ENTRY(entry4));

   graphics_info_t::check_waters_b_factor_limit = b_factor_lim;
   graphics_info_t::check_waters_map_sigma_limit = map_sigma_lim;
   graphics_info_t::check_waters_min_dist_limit = min_dist;
   graphics_info_t::check_waters_max_dist_limit = max_dist;

   bool use_b_factor_limit_test = 1;
   bool use_map_sigma_limit_test = 1;
   bool use_min_dist_test = 1;
   bool use_max_dist_test = 1;
   bool use_difference_map_test = 1;

   GtkWidget *hbox1 = lookup_widget(dialog, "check_waters_b_factor_hbox");
   GtkWidget *hbox2 = lookup_widget(dialog, "check_waters_sigma_level_hbox");
   GtkWidget *hbox3 = lookup_widget(dialog, "check_waters_min_dist_hbox");
   GtkWidget *hbox4 = lookup_widget(dialog, "check_waters_max_dist_hbox");

   GtkToggleButton *checkbutton1 =
      GTK_TOGGLE_BUTTON(lookup_widget(dialog, "check_waters_b_factor_entry_active_checkbutton")); 
   GtkToggleButton *checkbutton2 =
      GTK_TOGGLE_BUTTON(lookup_widget(dialog, "check_waters_map_sigma_entry_active_checkbutton"));
   GtkToggleButton *checkbutton3 =
      GTK_TOGGLE_BUTTON(lookup_widget(dialog, "check_waters_min_dist_entry_active_checkbutton"));
   GtkToggleButton *checkbutton4 =
      GTK_TOGGLE_BUTTON(lookup_widget(dialog, "check_waters_max_dist_entry_active_checkbutton"));
   GtkToggleButton *checkbutton5 =
      GTK_TOGGLE_BUTTON(lookup_widget(dialog, "check_waters_by_difference_map_active_checkbutton"));

   if (! checkbutton1->active)
      use_b_factor_limit_test = 0; 
   if (! checkbutton2->active)
      use_map_sigma_limit_test = 0; 
   if (! checkbutton3->active)
      use_min_dist_test = 0; 
   if (! checkbutton4->active)
      use_max_dist_test = 0;
   if (checkbutton5) 
      if (! checkbutton5->active)
	 use_difference_map_test = 0;
      
   GtkWidget *zero_occ_checkbutton = lookup_widget(dialog, "check_waters_zero_occ_checkbutton");
   GtkWidget *partial_occ_close_contact_checkbutton =
      lookup_widget(dialog, "check_waters_low_occ_dist_checkbutton");

   short int zero_occ_flag = 0;
   short int part_occ_dist_flag = 0;
   if (GTK_TOGGLE_BUTTON(zero_occ_checkbutton)->active)
      zero_occ_flag = 1;
   if (GTK_TOGGLE_BUTTON(partial_occ_close_contact_checkbutton)->active)
      part_occ_dist_flag = 1;
   
   //
   short int logical_operator_and_or_flag = 0; // logical AND
   if (GTK_TOGGLE_BUTTON(checklogic_OR_radiobutton)->active) {
      logical_operator_and_or_flag = 1;
   }

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(action_optionmenu));
   GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));
   int active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);

   // This will give us another dialog
   // 
   if (use_difference_map_test) {
      int imol_diff_map = graphics_info_t::check_waters_by_difference_map_map_number;
      check_waters_by_difference_map(graphics_info_t::check_waters_molecule, imol_diff_map, 1);
   } 

   if (use_b_factor_limit_test == 0)
      b_factor_lim = -100.0;
   if (use_map_sigma_limit_test == 0)
      map_sigma_lim = -100.0;
   if (use_min_dist_test == 0)
      min_dist = -100.0;
   if (use_max_dist_test == 0)
      max_dist = -100.0;  // sets a flag in find_water_baddies_OR
   if (active_index == 0) {
      GtkWidget *w = wrapped_checked_waters_baddies_dialog(graphics_info_t::check_waters_molecule,
							   b_factor_lim,
							   map_sigma_lim,
							   min_dist,
							   max_dist,
							   part_occ_dist_flag,
							   zero_occ_flag,
							   logical_operator_and_or_flag);
      gtk_widget_show(w);

   } else {

      // delete those baddies:
      delete_checked_waters_baddies(graphics_info_t::check_waters_molecule,
				    b_factor_lim,
				    map_sigma_lim,
				    min_dist,
				    max_dist,
				    part_occ_dist_flag,
				    zero_occ_flag,
				    logical_operator_and_or_flag); // calls graphics_draw()
   }
   
}

void store_checked_waters_baddies_dialog(GtkWidget *w) {
   graphics_info_t::checked_waters_baddies_dialog = w;
}



void check_waters_molecule_menu_item_activate(GtkWidget *item, 
					      GtkPositionType pos) {

   graphics_info_t::check_waters_molecule = pos;
}


// On check OK, we fire up this widget which is a vbox of baddy radio buttons.
// 
GtkWidget *wrapped_checked_waters_baddies_dialog(int imol, float b_factor_lim, float map_sigma_lim, float min_dist, float max_dist, short int part_occ_contact_flag, short int zero_occ_flag, short int logical_operator_and_or_flag) {

   GtkWidget *w = NULL;
   if (graphics_info_t::use_graphics_interface_flag) { 
      w = create_checked_waters_baddies_dialog();

      graphics_info_t g;
      int imol_for_map = g.Imol_Refinement_Map();

      if (is_valid_model_molecule(imol)) {
	 if (!is_valid_map_molecule(imol_for_map)) {
	    std::cout << "WARNING:: Not a valid map for density testing "
		      << imol_for_map << std::endl;
	 } else {

	    std::vector<coot::atom_spec_t> baddies = graphics_info_t::molecules[imol].find_water_baddies(b_factor_lim, graphics_info_t::molecules[imol_for_map].xmap_list[0], graphics_info_t::molecules[imol_for_map].map_sigma(), map_sigma_lim, min_dist, max_dist, part_occ_contact_flag, zero_occ_flag, logical_operator_and_or_flag);

	    // User data is used to keyboard up and down baddie water
	    // list (in graphics_info_t::checked_waters_next_baddie).
	    gtk_object_set_user_data(GTK_OBJECT(w), GINT_TO_POINTER(baddies.size()));
	 
	    GtkWidget *button;
	    GtkWidget *vbox = lookup_widget(w, "checked_waters_baddies_vbox");
	    GSList *gr_group = NULL;
	    
	    if (baddies.size() > 0 ) { 
	       for (int i=0; i<int(baddies.size()); i++) {
		  
		  // 	       std::cout << "Suspicious water: "
		  // 			 << baddies[i].atom_name
		  // 			 << baddies[i].alt_conf << " "
		  // 			 << baddies[i].resno << " "
		  // 			 << baddies[i].insertion_code << " "
		  // 			 << baddies[i].chain << "\n";
		  
		  std::string button_label(" ");
		  button_label += baddies[i].chain;
		  button_label += " " ;
		  button_label += graphics_info_t::int_to_string(baddies[i].resno);
		  button_label += " " ;
		  button_label += baddies[i].atom_name;
		  button_label += " " ;
		  button_label += baddies[i].alt_conf;
		  button_label += " [Occ: " ;
		  button_label += graphics_info_t::float_to_string(baddies[i].float_user_data);
		  button_label += "]   " ;
		  button_label += baddies[i].string_user_data;
		  button_label += " " ;
		  
		  button = gtk_radio_button_new_with_label(gr_group, button_label.c_str());
		  gr_group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
		  coot::atom_spec_t *atom_spec = new coot::atom_spec_t(baddies[i]);
		  atom_spec->int_user_data = imol;

		  std::string button_name = "checked_waters_baddie_button_";
		  button_name += coot::util::int_to_string(i);

 		  gtk_object_set_data_full(GTK_OBJECT(w),
 					   button_name.c_str(), button,
					   NULL);
		  
		  gtk_signal_connect(GTK_OBJECT(button), "clicked",
				     GTK_SIGNAL_FUNC (graphics_info_t::on_generic_atom_spec_button_clicked),
				     atom_spec);
		  
		  GtkWidget *frame = gtk_frame_new(NULL);
		  gtk_container_add(GTK_CONTAINER(frame), button);
		  
		  gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, 0);
		  gtk_container_set_border_width(GTK_CONTAINER(frame), 2);
		  gtk_widget_show(button);
		  gtk_widget_show(frame);
	       }
	    } else {
	 
	       std::string s = "There were no suspicious waters \nmatching those criteria in\n";
	       s += " Molecule ";
	       s += graphics_info_t::molecules[imol].dotted_chopped_name();
	       w = wrapped_nothing_bad_dialog(s);
	    }
	 }
      }
   }
   store_checked_waters_baddies_dialog(w);
   return w;
}


void delete_checked_waters_baddies(int imol, float b_factor_lim, float map_sigma_lim,
				   float min_dist,
				   float max_dist,
				   short int part_occ_contact_flag,
				   short int zero_occ_flag,
				   short int logical_operator_and_or_flag) {

   graphics_info_t g;
   int imol_for_map = g.Imol_Refinement_Map();

   if (is_valid_model_molecule(imol)) {
      if (!is_valid_map_molecule(imol_for_map)) {
	 std::cout << "WARNING:: Not a valid map for density testing " << imol_for_map << std::endl;
	 show_select_map_dialog();
      } else {
	 std::vector<coot::atom_spec_t> baddies =
	    graphics_info_t::molecules[imol].find_water_baddies(b_factor_lim,
								graphics_info_t::molecules[imol_for_map].xmap_list[0],
								graphics_info_t::molecules[imol_for_map].map_sigma(),
								map_sigma_lim,
								min_dist,
								max_dist,
								part_occ_contact_flag,
								zero_occ_flag,
								logical_operator_and_or_flag);
	 
	 int ideleted = graphics_info_t::molecules[imol].delete_atoms(baddies);
	 std::string s = "Deleted ";
	 s += graphics_info_t::int_to_string(ideleted);
	 s += " waters";
	 if (graphics_info_t::use_graphics_interface_flag) { 
	    GtkWidget *w = wrapped_nothing_bad_dialog(s);
	    gtk_widget_show(w);
	    graphics_draw();
	 }
      }
   }

}

void check_chiral_volumes(int imol) { 
   graphics_info_t g;
   if (imol < graphics_info_t::n_molecules()) { 
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



void set_fix_chiral_volumes_before_refinement(int istate) {
   graphics_info_t::fix_chiral_volume_before_refinement_flag = istate;
} 


void check_chiral_volumes_from_widget(GtkWidget *window) { 
   
   check_chiral_volumes(graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule);
}



// --------------------------------------------------------------------------
//                     difference map variance check 
// --------------------------------------------------------------------------
// 
void check_waters_by_difference_map(int imol_waters, int imol_diff_map,
				    int interactive_flag) {

   graphics_info_t g;
   g.check_waters_by_difference_map(imol_waters, imol_diff_map, interactive_flag);

}

void deviant_geometry(int imol) {

   if (is_valid_model_molecule(imol)) { 
      float strictness = 3.0;
      graphics_info_t::molecules[imol].find_deviant_geometry(strictness);
   }
}


short int is_valid_model_molecule(int imol) { 
   short int v = graphics_info_t::is_valid_model_molecule(imol);
   //   std::cout << "DEBUG:: in is_valid_model_molecule for " << imol << " returning "
   // << v << std::endl;
   return v;
} 


short int is_valid_map_molecule(int imol) { 

   return graphics_info_t::is_valid_map_molecule(imol);
}

#include "geometry-graphs.hh"

// Fee the local_block_info_p s
// 
void free_geometry_graph(GtkWidget *dialog) { 

#ifdef HAVE_GSL
   
   if (dialog) {
      GtkWidget *w = lookup_widget(dialog, "geometry_graph_canvas");
      if (w) { 
	 GtkObject *obj = GTK_OBJECT(w);
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
	 coot::geometry_graphs *graphs = (coot::geometry_graphs *) gtk_object_get_user_data(obj);
	 if (!graphs) {
	    std::cout << "ERROR:: NULL graphs in free_geometry_graph\n";
	 } else {
	    graphs->close_yourself();
	 }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
      }
   }
#endif // HAVE_GSL   
} 

void unset_geometry_graph(GtkWidget *dialog) {  /* set the graphics info
						 static to NULL, so
						 that we on longer try
						 to update the
						 widget*/
#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   graphics_info_t g;

   int imol;
   if (dialog) {
      GtkWidget *w = lookup_widget(dialog, "geometry_graph_canvas");
      if (w) { 
	 GtkObject *obj = GTK_OBJECT(w);
	 coot::geometry_graphs *graphs = (coot::geometry_graphs *) gtk_object_get_user_data(obj);
	 if (!graphs) {

	    std::cout << "ERROR:: NULL graphs in unset_geometry_graph\n";

	 } else { 
	 
	    imol = graphs->Imol();
	 
	    if (is_valid_model_molecule(imol)) {
	       
	       coot::set_validation_graph(imol, graphs->Graph_Type(), NULL);

// Old style array of molecules code
// 	       if (graphs->Graph_Type() == coot::GEOMETRY_GRAPH_GEOMETRY) {
// 		  g.geometry_graph[imol] = NULL;
// 	       }
// 	       if (graphs->Graph_Type() == coot::GEOMETRY_GRAPH_B_FACTOR) {
// 		  g.b_factor_variance_graph[imol] = NULL;
// 	       }
// 	       if (graphs->Graph_Type() == coot::GEOMETRY_GRAPH_OMEGA_DISTORTION) {
// 		  g.omega_distortion_graph[imol] = NULL;
// 	       }
// 	       if (graphs->Graph_Type() == coot::GEOMETRY_GRAPH_ROTAMER) {
// 		  g.rotamer_graph[imol] = NULL;
// 	       }
// 	       if (graphs->Graph_Type() == coot::GEOMETRY_GRAPH_DENSITY_FIT) {
// 		  g.residue_density_fit_graph[imol] = NULL;
// 	       }
// 	       if (graphs->Graph_Type() == coot::GEOMETRY_GRAPH_NCS_DIFFS) {
// 		  g.ncs_diffs_graph[imol] = NULL;
// 	       }

	    }
	 }
      } else {
	 std::cout << "Failed to find w in unset_geometry_graph\n";
      } 
   }
   // std::cout << "Done unset_geometry_graph\n";
#endif 
#endif // HAVE_GSL   
}


#ifdef USE_GUILE 
/*! \brief Activate rotamer graph analysis for molecule number imol.  

Return rotamer info - function used in testing.  */
SCM rotamer_graphs(int imol) {

   SCM r = SCM_BOOL_F;

   graphics_info_t g;
   coot::rotamer_graphs_info_t results = g.rotamer_graphs(imol);
   if (results.info.size() > 0) {
      r = SCM_EOL;
      for (int ir=int(results.info.size()-1); ir>=0; ir--) {
	 SCM ele = SCM_EOL;
	 SCM name = scm_makfrom0str(results.info[ir].rotamer_name.c_str());
	 ele = scm_cons(name, ele);
	 SCM pr = scm_double2num(results.info[ir].probability);
	 ele = scm_cons(pr, ele);
	 SCM inscode = scm_makfrom0str(results.info[ir].inscode.c_str());
	 ele = scm_cons(inscode, ele);
	 SCM resno = scm_int2num(results.info[ir].resno);
	 ele = scm_cons(resno, ele);
	 SCM chainid = scm_makfrom0str(results.info[ir].chain_id.c_str());
	 ele = scm_cons(chainid, ele);

	 r = scm_cons(ele, r);
      }
   } 

   return r;
} 
#endif

#ifdef USE_PYTHON
/*! \brief Activate rotamer graph analysis for molecule number imol.  

Return rotamer info - function used in testing.  */
PyObject *rotamer_graphs_py(int imol) {

   PyObject *r;
   r = Py_False;

   graphics_info_t g;
   coot::rotamer_graphs_info_t results = g.rotamer_graphs(imol);
   if (results.info.size() > 0) {
      r = PyList_New(results.info.size());
      for (int ir=int(results.info.size()-1); ir>=0; ir--) {
	 PyObject *ele = PyList_New(5);
	 PyObject *name = PyString_FromString(results.info[ir].rotamer_name.c_str());
	 PyList_SetItem(ele, 4, name);;
	 PyObject *pr = PyFloat_FromDouble(results.info[ir].probability);
	 PyList_SetItem(ele, 3, pr);
	 PyObject *inscode = PyString_FromString(results.info[ir].inscode.c_str());
	 PyList_SetItem(ele, 2, inscode);
	 PyObject *resno = PyInt_FromLong(results.info[ir].resno);
	 PyList_SetItem(ele, 1, resno);
	 PyObject *chainid = PyString_FromString(results.info[ir].chain_id.c_str());
	 PyList_SetItem(ele, 0, chainid);

         PyList_SetItem(r, ir, ele);
      }
   } 

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }

   return r;
}
#endif // USE_PYTHON


// -----------------------------------------------------
// The geometry graphs have a home on the range:
// -----------------------------------------------------

void add_on_validation_graph_mol_options(GtkWidget *menu, const char *type_in) {

#ifdef HAVE_GSL
   graphics_info_t g;
   std::string validation_type(type_in);
   std::string sub_menu_name;
   GtkSignalFunc callback = 0; // depends on type
   short int found_validation_type = 0;

   if (validation_type == "b factor") {
      callback = GTK_SIGNAL_FUNC(validation_graph_b_factor_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "temp_factor_variance_submenu";
   }
   if (validation_type == "geometry") {
      callback = GTK_SIGNAL_FUNC(validation_graph_geometry_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "geometry_submenu";
   }
   if (validation_type == "omega") {
      callback = GTK_SIGNAL_FUNC(validation_graph_omega_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "omega_submenu";
   }
   if (validation_type == "rotamer") {
      callback = GTK_SIGNAL_FUNC(validation_graph_rotamer_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "rotamer_submenu";
   }
   if (validation_type == "density-fit") {
      callback = GTK_SIGNAL_FUNC(validation_graph_density_fit_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "density_fit_submenu";
   }
   if (validation_type == "probe") {
      callback = GTK_SIGNAL_FUNC(probe_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "probe_submenu";
   }
   if (validation_type == "gln_and_asn_b_factor_outliers") {
      callback = GTK_SIGNAL_FUNC(gln_and_asn_b_factor_outlier_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "gln_and_asn_b_factor_outliers_submenu";
   }
   if (validation_type == "ncs-diffs") {
      callback = GTK_SIGNAL_FUNC(validation_graph_ncs_diffs_mol_selector_activate);
      found_validation_type = 1;
      sub_menu_name = "ncs_diffs_submenu";
   }

   GtkWidget *sub_menu = lookup_widget(menu, sub_menu_name.c_str());

   if (sub_menu) { 
   
      gtk_container_foreach(GTK_CONTAINER(sub_menu),
			    my_delete_validaton_graph_mol_option,
			    (gpointer) sub_menu);

      for(int i=0; i<g.n_molecules(); i++) {
	 if (g.molecules[i].has_model()) {
	    std::string name;
	    name = graphics_info_t::molecules[i].dotted_chopped_name();
	    add_validation_mol_menu_item(i, name, sub_menu, callback);
	 }
      }
   } else {
      std::cout << "ERROR:: sub menu not found: " << sub_menu_name << std::endl;
   } 
#endif // HAVE_GSL
}

void
add_validation_mol_menu_item(int imol,
			     const std::string &name,
			     GtkWidget *menu,
			     GtkSignalFunc callback) {

   GtkWidget *menu_item = gtk_menu_item_new_with_label(name.c_str());
   gtk_container_add(GTK_CONTAINER(menu), menu_item);
   gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
		      callback, GINT_TO_POINTER(imol));
   gtk_widget_show(menu_item);
}


void
my_delete_validaton_graph_mol_option(GtkWidget *widget, void *data) {
   gtk_container_remove(GTK_CONTAINER(data), widget);
}

void validation_graph_b_factor_mol_selector_activate (GtkMenuItem     *menuitem,
						      gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
      graphics_info_t g;
      g.b_factor_graphs(imol);
#else    
      printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */

}

void validation_graph_geometry_mol_selector_activate (GtkMenuItem     *menuitem,
						      gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
      graphics_info_t g;
      g.geometric_distortion(imol);
#else    
      printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */

}

void validation_graph_omega_mol_selector_activate (GtkMenuItem     *menuitem,
						   gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
      graphics_info_t g;
      g.omega_graphs(imol);
#else    
      printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */

}

void validation_graph_rotamer_mol_selector_activate (GtkMenuItem     *menuitem,
						     gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   graphics_info_t g;
   g.rotamer_graphs(imol);
#else    
   printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */

}

void validation_graph_density_fit_mol_selector_activate (GtkMenuItem     *menuitem,
							 gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   graphics_info_t g;
   g.density_fit_graphs(imol);
#else    
   printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */
}

void probe_mol_selector_activate (GtkMenuItem     *menuitem,
				  gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
      graphics_info_t g;

      std::vector<std::string> cmd_strings;
      cmd_strings.push_back("probe");
      cmd_strings.push_back(coot::util::int_to_string(imol));
      run_generic_script(cmd_strings);
#else    
      printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */
}

// is the probe executable available?
// 1 for yes, 0 for no.
// 
int probe_available_p() {
   int r=0;

#if defined(USE_GUILE) && !defined(WINDOWS_MINGW)

   std::string command("(command-in-path? *probe-command*)");

   SCM scm_thunk = safe_scheme_command(command); 

   int was_boolean_flag = gh_scm2bool(scm_boolean_p(scm_thunk));

   if (was_boolean_flag)
      if (gh_scm2bool(scm_thunk) == 1)
	 r = 1;

#else
// BL says:: here comes some (experimental) code to do the same thing in python
// it's a bit longer and uses compiled python code but dont see another option
// currently (and there might be none..)
#ifdef USE_PYTHON

    PyObject *result;
    result = safe_python_command_with_return("command_in_path_qm(probe_command)");

    int was_boolean_flag = PyInt_AsLong(result);
    if (was_boolean_flag) {
              r = 1;
//	      std::cout << "BL DEBUG:: flag here " << was_boolean_flag << std::endl;
    }

//    std::cout << "BL DEBUG:: r here " << r << std::endl;

#endif // USE_PYTHON
#endif // USE_GUILE   

   return r;
} 

#ifdef USE_PYTHON
// is the probe executable available?
// 1 for yes, 0 for no.
// 
int probe_available_p_py() {
    int r=0;

    PyObject *result;
    result = safe_python_command_with_return("command_in_path_qm(probe_command)");

    int was_boolean_flag = PyInt_AsLong(result);
    if (was_boolean_flag) {
              r = 1;
    }
   return r;
} 
#endif // USE_PYTHON

void gln_and_asn_b_factor_outlier_mol_selector_activate (GtkMenuItem     *menuitem,
							 gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   gln_asn_b_factor_outliers(imol);
}


void
create_initial_validation_graph_submenu_generic(GtkWidget *widget,
						const std::string &menu_name,
						const std::string &sub_menu_name) {

   GtkWidget *b_factor_menu_item = lookup_widget(widget, menu_name.c_str());
   GtkWidget *b_factor_sub_menu = gtk_menu_new();
   gtk_widget_ref(b_factor_sub_menu);
   gtk_object_set_data_full(GTK_OBJECT(widget),
			    sub_menu_name.c_str(),
			    b_factor_sub_menu,
			    (GtkDestroyNotify) gtk_widget_unref);

   gtk_menu_item_set_submenu(GTK_MENU_ITEM(b_factor_menu_item),
			     b_factor_sub_menu);
   
}

// ---------------------------------------------------------------------
//                difference map
// ---------------------------------------------------------------------
//
// where level is in sigma
// 
void
difference_map_peaks(int imol, int imol_coords,
		     float n_sigma,
		     int do_positive_level_flag,
		     int do_negative_levels_flag) {

   // Notice that we make wrapped_create_check_waters_dialog be part
   // of graphics_info_t, because it uses clipper data in the
   // interface - and c-interface.h does not know about clipper.


   // I don't think we want ligand/cluster search.  We just want peak
   // searching.
   // 
   if (is_valid_map_molecule(imol)) {
      if (graphics_info_t::molecules[imol].is_difference_map_p()) {

	 // c.f. trace-high-res.cc
	 coot::peak_search ps(graphics_info_t::molecules[imol].xmap_list[0]);
	 std::vector<std::pair<clipper::Coord_orth, float> > centres;

	 if (is_valid_model_molecule(imol_coords)) {
	    centres =
	       ps.get_peaks(graphics_info_t::molecules[imol].xmap_list[0],
			    graphics_info_t::molecules[imol_coords].atom_sel.mol,
			    n_sigma, do_positive_level_flag, do_negative_levels_flag);
	 } else { 
	    centres =
	       ps.get_peaks(graphics_info_t::molecules[imol].xmap_list[0],
			    n_sigma, do_positive_level_flag, do_negative_levels_flag);
	 }
	 
	 if (centres.size() == 0) {
	    std::string info_string("No difference map peaks\nat ");
	    info_string += graphics_info_t::float_to_string(n_sigma);
	    info_string += " sigma";
	    GtkWidget *w = wrapped_nothing_bad_dialog(info_string);
	    gtk_widget_show(w);
	 } else {
	    float map_sigma = graphics_info_t::molecules[imol].map_sigma();
	    GtkWidget *w = graphics_info_t::wrapped_create_diff_map_peaks_dialog(centres, map_sigma);
	    gtk_widget_show(w);

	    std::cout << "\n   Found these peak positions:\n";
	    for (unsigned int i=0; i<centres.size(); i++) {
	       std::cout << "   " << i << " " << centres[i].second << " "
			 << centres[i].first.format() << std::endl;
	    }
	    std::cout << "\n   Found " << centres.size() << " peak positions:\n";
	 }
      }
   } else {
      std::cout << "Molecule number " << imol
		<< " is not a valid map molecule" << std::endl;
   }
}

void set_difference_map_peaks_widget(GtkWidget *w) {
   graphics_info_t::difference_map_peaks_dialog = w;
} 


void
clear_diff_map_peaks() {
   
   graphics_info_t g;
   g.clear_diff_map_peaks();
}


GtkWidget *wrapped_create_generate_diff_map_peaks_dialog() {

   // c.f. wrapped_create_check_waters_diff_map_dialog()

   GtkWidget *dialog = create_generate_diff_map_peaks_dialog();

   int ifound;
   short int diff_maps_only_flag = 1;

   ifound = fill_ligands_dialog_map_bits_by_dialog_name(dialog, "generate_diff_map_peaks_map", 
							diff_maps_only_flag);
   if (ifound == 0) {
      std::cout << "Error: you must have a difference map to analyse!" << std::endl;
      GtkWidget *none_frame = lookup_widget(dialog, "no_difference_maps_frame");
      gtk_widget_show(none_frame);
   }

   // the name of the vbox which is looked up is "generate_diff_map_peaks_model_vbox".
   ifound = fill_ligands_dialog_protein_bits_by_dialog_name(dialog, "generate_diff_map_peaks_model");
   if (ifound == 0) {
      std::cout << "Difference map checker is better having specified coordinates...\n";
   }

   // the sigma entry:
   GtkWidget *entry = lookup_widget(dialog, "generate_diff_map_peaks_sigma_level_entry");
   gtk_entry_set_text(GTK_ENTRY(entry),
		      graphics_info_t::float_to_string(graphics_info_t::difference_map_peaks_sigma_level).c_str());

   
return dialog;

}

void difference_map_peaks_by_widget(GtkWidget *dialog) {

   // c.f. check_waters_by_difference_map_by_widget(GtkWidget *dialog)

   short int found_active_button_for_map = 0;
   // short int found_active_button_for_coords = 0; // but not terminal
   int imol_diff_map = -1;
   int imol_coords = -1;

   // Check the difference map:

   // the strings correspond to the above function with _radiobuton_ tagged on.
   
   GtkWidget *map_button;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_map()) {
	 if (graphics_info_t::molecules[imol].is_difference_map_p()) {
	    std::string map_str = "generate_diff_map_peaks_map_radiobutton_";
	    map_str += graphics_info_t::int_to_string(imol);
	    map_button = lookup_widget(dialog, map_str.c_str());
	    if (map_button) {
	       if (GTK_TOGGLE_BUTTON(map_button)->active) {
		  imol_diff_map = imol;
		  found_active_button_for_map = 1;
	       }
	    } else {
	       std::cout << "ooops failed to find button" << map_str << std::endl;
	    }
	 }
      }
   }

   // Check the coords:
   
   GtkWidget *coords_button;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 std::string coords_str = "generate_diff_map_peaks_model_radiobutton_";
	 coords_str += graphics_info_t::int_to_string(imol);
	 coords_button = lookup_widget(dialog, coords_str.c_str());
	 if (coords_button) {
	    if (GTK_TOGGLE_BUTTON(coords_button)->active) {
	       imol_coords = imol;
	       found_active_button_for_map = 1;
	    }
	 } else {
	    std::cout << "ooops failed to find button" << coords_str << std::endl;
	 }
      }
   }

   // Check the level:

   GtkWidget *sigma_entry =
      lookup_widget(dialog, "generate_diff_map_peaks_sigma_level_entry");

   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(sigma_entry));
   float v = atof(txt);

   short int good_sigma = 0;
   if (v > -1000 && v < 1000) {
      good_sigma = 1;
   } else {
      std::cout << "WARNING:: Invalid sigma level: " << v
		<< " can't do peak search." << std::endl;
   }

   // Check the negative level checkbutton:

   short int do_negative_level = 0;
   short int do_positive_level = 0;
   GtkWidget *checkbutton_negative =
      lookup_widget(dialog, "generate_diff_map_peaks_negative_level_checkbutton");
   GtkWidget *checkbutton_positive =
      lookup_widget(dialog, "generate_diff_map_peaks_positive_level_checkbutton");

   if (GTK_TOGGLE_BUTTON(checkbutton_negative)->active)
      do_negative_level = 1;

   if (GTK_TOGGLE_BUTTON(checkbutton_positive)->active)
      do_positive_level = 1;

   if (found_active_button_for_map) {
      if (good_sigma)
	 // if imol_coords is -1 it is ignored in difference_map_peaks
	 difference_map_peaks(imol_diff_map, imol_coords, v,
			      do_positive_level, do_negative_level);
   } else {
      std::cout << "WARNING:: failed to find a difference map "
		<< "Can't do peak search" << std::endl;
      GtkWidget *w = wrapped_nothing_bad_dialog("WARNING:: failed to find difference map\nCan't do peak search");
      gtk_widget_show(w);
   }
}


// ----------------------------------------------------------------------------------
//                            ramachandran plot
// ----------------------------------------------------------------------------------


// fill the widget:
GtkWidget *wrapped_ramachandran_plot_differences_dialog() {

   GtkWidget *w = 0; // Not NULL, compiler (maybe).

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   w = create_ramachandran_plot_differences_dialog();

   // We don't have to worry about chains because they are not active on startup.

   GtkWidget *optionmenu1 =
      lookup_widget(w, "ramachandran_plot_differences_first_mol_optionmenu");
   GtkWidget *optionmenu2 =
      lookup_widget(w, "ramachandran_plot_differences_second_mol_optionmenu");

   GtkSignalFunc signal_func1 =
      GTK_SIGNAL_FUNC(ramachandran_plot_differences_mol_option_menu_activate_first);
   GtkSignalFunc signal_func2 =
      GTK_SIGNAL_FUNC(ramachandran_plot_differences_mol_option_menu_activate_second);

   int imol = -1;
   for (int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }

   if (imol >= 0) {
      graphics_info_t g;
      g.fill_option_menu_with_coordinates_options(optionmenu1, signal_func1, imol);
      g.fill_option_menu_with_coordinates_options(optionmenu2, signal_func2, imol);
      graphics_info_t::ramachandran_plot_differences_imol1 = imol;
      graphics_info_t::ramachandran_plot_differences_imol2 = imol;
   }

#endif
   return w;
}


// the OK button callback.
// (read the widget)
int do_ramachandran_plot_differences_by_widget(GtkWidget *w) {

   std::string first_chain, second_chain;
   int istat = 0; // incomprehensible

   int imol1 = graphics_info_t::ramachandran_plot_differences_imol1;
   int imol2 = graphics_info_t::ramachandran_plot_differences_imol2;

   first_chain  = graphics_info_t::ramachandran_plot_differences_imol1_chain;
   second_chain = graphics_info_t::ramachandran_plot_differences_imol2_chain;

   GtkWidget *checkbutton1 = lookup_widget(GTK_WIDGET(w),
					  "ramachandran_plot_differences_first_chain_checkbutton");
   GtkWidget *checkbutton2 = lookup_widget(GTK_WIDGET(w),
					  "ramachandran_plot_differences_second_chain_checkbutton");

   if (GTK_TOGGLE_BUTTON(checkbutton1)->active && GTK_TOGGLE_BUTTON(checkbutton2)->active) {
      istat = 1;
      ramachandran_plot_differences_by_chain(imol1, imol2,
					     first_chain.c_str(),
					     second_chain.c_str());
   } else {
      if (!GTK_TOGGLE_BUTTON(checkbutton1)->active && !GTK_TOGGLE_BUTTON(checkbutton2)->active) {
	 istat = 1;
	 ramachandran_plot_differences(imol1, imol2);
      } else {
	 std::cout << "INFO:: incomprehensible molecule/chain selection" << std::endl;
	 std::string s = "Can't make sense of chain selection.  Try again?";
	 GtkWidget *w = wrapped_nothing_bad_dialog(s);
	 gtk_widget_show(w);
      }
   }
   return istat;
}

/*! \brief set the contour levels for theremachandran plot, default
  values are 0.02 (prefered) 0.002 (allowed) */
void set_ramachandran_plot_contour_levels(float level_prefered, float level_allowed) {

   graphics_info_t::rama_level_prefered = level_prefered;
   graphics_info_t::rama_level_allowed  = level_allowed;
} 

void set_ramachandran_plot_background_block_size(float blocksize) {
   graphics_info_t::rama_plot_background_block_size = blocksize;
}




// OK, the molecule was changed in the option menu, so if the checkbutton is on,
// then change the elements of the chain option menu
// 
void ramachandran_plot_differences_mol_option_menu_activate_first(GtkWidget *item, GtkPositionType pos) {
   graphics_info_t::ramachandran_plot_differences_imol1 = pos;
   GtkWidget *chain_optionmenu = lookup_widget(GTK_WIDGET(item),
						"ramachandran_plot_differences_first_chain_optionmenu");
   GtkWidget *checkbutton = lookup_widget(GTK_WIDGET(item),
					  "ramachandran_plot_differences_first_chain_checkbutton");
   if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
      fill_ramachandran_plot_differences_option_menu_with_chain_options(chain_optionmenu, 1);
   }
}

void ramachandran_plot_differences_mol_option_menu_activate_second(GtkWidget *item, GtkPositionType pos) {
   graphics_info_t::ramachandran_plot_differences_imol2 = pos; 
   GtkWidget *chain_optionmenu = lookup_widget(GTK_WIDGET(item),
						"ramachandran_plot_differences_second_chain_optionmenu");
   GtkWidget *checkbutton = lookup_widget(GTK_WIDGET(item),
					  "ramachandran_plot_differences_second_chain_checkbutton");
   if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
      fill_ramachandran_plot_differences_option_menu_with_chain_options(chain_optionmenu, 0);
   }
}


void ramachandran_plot_differences_chain_option_menu_activate_first(GtkWidget *item, GtkPositionType pos){

   graphics_info_t::ramachandran_plot_differences_imol1_chain = menu_item_label(item);
}

void ramachandran_plot_differences_chain_option_menu_activate_second(GtkWidget *item, GtkPositionType pos){

      graphics_info_t::ramachandran_plot_differences_imol2_chain = menu_item_label(item);
}

void do_ramachandran_plot(int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules()) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    coot::rama_plot *rama = new coot::rama_plot;
	    if (graphics_info_t::ramachandran_plot_x_position > 0)
	       rama->set_position(graphics_info_t::ramachandran_plot_x_position,
				  graphics_info_t::ramachandran_plot_y_position);
	    short int is_kleywegt_plot_flag = 0;
	    rama->set_n_diffs(graphics_info_t::rama_n_diffs);
	    rama->init(imol,
		       graphics_info_t::molecules[imol].dotted_chopped_name(), 
		       graphics_info_t::rama_level_prefered,
		       graphics_info_t::rama_level_allowed,
		       graphics_info_t::rama_plot_background_block_size,
		       is_kleywegt_plot_flag);
	    rama->draw_it(graphics_info_t::molecules[imol].atom_sel.mol); 
	 }
      }
   }
#endif // HAVE_GTK_CANVAS
}

void set_kleywegt_plot_n_diffs(int ndiffs) {
   graphics_info_t::rama_n_diffs = ndiffs;
}


void
ramachandran_plot_differences(int imol1, int imol2) { 

   
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *rama_widget = dynarama_is_displayed_state(imol1); // the underlying variable is
                                                                // set by the init_internal() function
                                                                // of rama_plot.
   if (rama_widget) {
      std::string s = "Sorry. Can't have Kleywegt Plot reference molecule the same\n";
      s += "an existing displayed Ramachandran Plot";
      info_dialog(s.c_str());
   } else { 
      if (imol1 >= 0) {
	 if (imol1 < graphics_info_t::n_molecules()) {  
	    if (graphics_info_t::molecules[imol1].has_model()) { 
	       if (imol2 >= 0) {
		  if (imol2 < graphics_info_t::n_molecules()) { 
		     if (graphics_info_t::molecules[imol2].has_model()) { 
			coot::rama_plot *rama = new coot::rama_plot;
			short int is_kleywegt_plot_flag = 1;
			rama->set_n_diffs(graphics_info_t::rama_n_diffs);
			rama->init(imol1,
				   graphics_info_t::molecules[imol1].dotted_chopped_name(),  
				   graphics_info_t::rama_level_prefered,
				   graphics_info_t::rama_level_allowed,
				   graphics_info_t::rama_plot_background_block_size,
				   is_kleywegt_plot_flag);
			rama->draw_it(imol1, imol2,
				      graphics_info_t::molecules[imol1].atom_sel.mol,
				      graphics_info_t::molecules[imol2].atom_sel.mol);
		     }
		  }
	       }
	    } 
	 }
      }
   }
#endif // HAVE_GTK_CANVAS
}

void ramachandran_plot_differences_by_chain(int imol1, int imol2, 
					    const char *a_chain, const char *b_chain) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   GtkWidget *rama_widget = dynarama_is_displayed_state(imol1); // the underlying variable is
                                                                // set by the init_internal() function
                                                                // of rama_plot.
   if (rama_widget) {
      std::string s = "Sorry. Can't have Kleywegt Plot reference molecule the same\n";
      s += "an existing displayed Ramachandran Plot";
      info_dialog(s.c_str());
   } else { 
      if (is_valid_model_molecule(imol1)) {
	 if (is_valid_model_molecule(imol2)) {
	    short int is_kleywegt_plot_flag = 0;
	    coot::rama_plot *rama = new coot::rama_plot; 
	    rama->set_n_diffs(graphics_info_t::rama_n_diffs);
	    rama->init(imol1,
		       graphics_info_t::molecules[imol1].dotted_chopped_name(),
		       graphics_info_t::rama_level_prefered,
		       graphics_info_t::rama_level_allowed,
		       graphics_info_t::rama_plot_background_block_size,
		       is_kleywegt_plot_flag);
	    if (graphics_info_t::molecules[imol1].is_from_shelx_ins() ||
		graphics_info_t::molecules[imol2].is_from_shelx_ins())
	       rama->allow_seqnum_offset();
// 	    std::cout << "rama differences on mols: " << imol1 << " " << a_chain
// 		      << " to " << imol2 << " " << b_chain << std::endl;
	    rama->draw_it(imol1, imol2,
			  graphics_info_t::molecules[imol1].atom_sel.mol,
			  graphics_info_t::molecules[imol2].atom_sel.mol,
			  std::string(a_chain), std::string(b_chain));
	    rama->set_kleywegt_plot_uses_chain_ids();
	 } else {
	    std::cout << "WARNING no molecule number " << imol2 << " in molecule (or closed)\n";
	 }
      } else {
	 std::cout << "WARNING no molecule number " << imol1 << " in molecule (or closed)\n";
      }
   }
#endif // HAVE_GTK_CANVAS
}


// This is called when the "Use chain" option button is pressed.
// 
void fill_ramachandran_plot_differences_option_menu_with_chain_options(GtkWidget *chain_optionmenu, 
								       int is_first_mol_flag) {
   GtkWidget *mol_optionmenu = NULL;

   if (is_first_mol_flag) {
      mol_optionmenu =
	 lookup_widget(chain_optionmenu,
		       "ramachandran_plot_differences_first_mol_optionmenu");
   } else {
      mol_optionmenu =
	 lookup_widget(chain_optionmenu,
		       "ramachandran_plot_differences_second_mol_optionmenu");
   }

   GtkSignalFunc callback_func;
   int imol;

   if (is_first_mol_flag) { 
      imol = graphics_info_t::ramachandran_plot_differences_imol1;
      callback_func =
	 GTK_SIGNAL_FUNC(ramachandran_plot_differences_chain_option_menu_activate_first);
   } else {
      imol = graphics_info_t::ramachandran_plot_differences_imol2;
      callback_func =
	 GTK_SIGNAL_FUNC(ramachandran_plot_differences_chain_option_menu_activate_second);
   }

   if (imol >=0 && imol< graphics_info_t::n_molecules()) {
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_optionmenu, imol, callback_func);
      if (is_first_mol_flag) {
	 graphics_info_t::ramachandran_plot_differences_imol1_chain = set_chain;
      } else {
	 graphics_info_t::ramachandran_plot_differences_imol2_chain = set_chain;
      }
   } else {
      std::cout << "ERROR:: in imol in fill_rama plot diffs: " << imol << std::endl;
   } 

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

#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   if (is_valid_model_molecule(imol)) {
      // w = graphics_info_t::dynarama_is_displayed[imol];
      w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
   }
#endif   
   return w;
}


GtkWidget *dynarama_widget(int imol) {

   GtkWidget *w = NULL;
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   if (imol < graphics_info_t::n_molecules()) {
      // w = graphics_info_t::dynarama_is_displayed[imol];
      w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
   }
#endif    
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


// ----------------------------------------------------------------------------------
//                            Alignment
// ----------------------------------------------------------------------------------
/*! \brief do a internal alignment of all the assigned sequences,
  return a list of mismatches.

Return scheme false on failure to align (e.g. not assigned sequence)
and the empty list on no alignment mismatches.

If there are any differences of any kind, return a list of 3 list:
mutations deletions insertions.

*/

#ifdef USE_GUILE
SCM alignment_mismatches_scm(int imol) {

   SCM r = SCM_BOOL_F;
   std::vector<std::pair<coot::residue_spec_t,std::string> > mutations;
   std::vector<std::pair<coot::residue_spec_t,std::string> > insertions;
   std::vector<std::pair<coot::residue_spec_t,std::string> > deletions;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::vector<coot::chain_mutation_info_container_t> > ar = 
      	 graphics_info_t::molecules[imol].residue_mismatches();
      if (ar.first)
	 r = SCM_EOL;
      for (unsigned int ir=0; ir<ar.second.size(); ir++) {
	 for (unsigned int im=0; im<ar.second[ir].mutations.size(); im++) {
	    mutations.push_back(ar.second[ir].mutations[im]);
	 }
	 for (unsigned int is=0; is<ar.second[ir].single_insertions.size(); is++) {
	    insertions.push_back(ar.second[ir].single_insertions[is]);
	 }
	 for (unsigned int id=0; id<ar.second[ir].deletions.size(); id++) {
	    coot::residue_spec_t del = ar.second[ir].deletions[id];
	    std::pair<coot::residue_spec_t, std::string> d(del, "Delete");
	    deletions.push_back(d);
	 }
      }
   }

   if ((mutations.size() > 0) || (insertions.size() > 0) || (deletions.size() > 0)) {
      SCM insertions_scm = SCM_EOL;
      SCM deletions_scm = SCM_EOL;
      SCM mutations_scm = SCM_EOL;
      for (unsigned int i=0; i<mutations.size(); i++) {
	 SCM rs_scm = scm_residue(mutations[i].first);
	 SCM str = scm_makfrom0str(mutations[i].second.c_str());
	 SCM c = SCM_EOL;
	 c = scm_cons(str, c);
	 c = scm_cons(str, rs_scm);
	 mutations_scm = scm_cons(c, mutations_scm);
      }
      for (unsigned int i=0; i<insertions.size(); i++) {
	 SCM rs_scm = scm_residue(insertions[i].first);
	 SCM str = scm_makfrom0str(insertions[i].second.c_str());
	 SCM c = SCM_EOL;
	 c = scm_cons(str, c);
	 c = scm_cons(str, rs_scm);
	 insertions_scm = scm_cons(c, insertions_scm);
      }
      for (unsigned int i=0; i<deletions.size(); i++) {
	 SCM rs_scm = scm_residue(deletions[i].first);
	 SCM str = scm_makfrom0str(deletions[i].second.c_str());
	 SCM c = SCM_EOL;
	 c = scm_cons(str, c);
	 c = scm_cons(str, rs_scm);
	 deletions_scm = scm_cons(c, deletions_scm);
      }
      r = SCM_EOL;
      // These are reversed so that the residue numbers come out in
      // numerical order (not backwards) and the returned list is
      // (list mutations deletions insertions).
      r = scm_cons(scm_reverse(insertions_scm), r);
      r = scm_cons(scm_reverse(deletions_scm),  r);
      r = scm_cons(scm_reverse(mutations_scm),  r);
   } 
   return r;
}
#endif // USE_GUILE


/*! \brief do a internal alignment of all the assigned sequences,
  return a list of mismatches that need to be made to model number
  imol to match the input sequence.

Return a list of mutations deletions insetions.
Return  False on failure to align (e.g. not assigned sequence)
and the empty list on no alignment mismatches.*/
#ifdef USE_PYTHON
PyObject *alignment_mismatches_py(int imol) {

   PyObject *r = Py_False;

   std::vector<std::pair<coot::residue_spec_t,std::string> > mutations;
   std::vector<std::pair<coot::residue_spec_t,std::string> > insertions;
   std::vector<std::pair<coot::residue_spec_t,std::string> > deletions;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::vector<coot::chain_mutation_info_container_t> > ar = 
      	 graphics_info_t::molecules[imol].residue_mismatches();
      if (ar.first)
	r = PyList_New(0);
      for (unsigned int ir=0; ir<ar.second.size(); ir++) {
	 for (unsigned int im=0; im<ar.second[ir].mutations.size(); im++) {
	    mutations.push_back(ar.second[ir].mutations[im]);
	 }
	 for (unsigned int is=0; is<ar.second[ir].single_insertions.size(); is++) {
	    insertions.push_back(ar.second[ir].single_insertions[is]);
	 }
	 for (unsigned int id=0; id<ar.second[ir].deletions.size(); id++) {
	    coot::residue_spec_t del = ar.second[ir].deletions[id];
	    std::pair<coot::residue_spec_t, std::string> d(del, "Delete");
	    deletions.push_back(d);
	 }
      }
   }

   if ((mutations.size() > 0) || (insertions.size() > 0) || (deletions.size() > 0)) {
     PyObject *insertions_py = PyList_New(0);
     PyObject *deletions_py = PyList_New(0);
     PyObject * mutations_py = PyList_New(0);
     for (unsigned int i=0; i<mutations.size(); i++) {
	 PyObject *rs_py = py_residue(mutations[i].first);
	 PyObject *str = PyString_FromString(mutations[i].second.c_str());
	 PyList_Insert(rs_py, 0, str);
	 PyList_Append(mutations_py, rs_py);
	 Py_XDECREF(str);
	 Py_XDECREF(rs_py);
      }
      for (unsigned int i=0; i<insertions.size(); i++) {
	 PyObject *rs_py = py_residue(insertions[i].first);
	 PyObject *str = PyString_FromString(insertions[i].second.c_str());
	 PyList_Insert(rs_py, 0, str);
	 PyList_Append(insertions_py, rs_py);
	 Py_XDECREF(str);
	 Py_XDECREF(rs_py);
      }
      for (unsigned int i=0; i<deletions.size(); i++) {
	 PyObject *rs_py = py_residue(deletions[i].first);
	 PyObject *str = PyString_FromString(deletions[i].second.c_str());
	 PyList_Insert(rs_py, 0, str);
	 PyList_Append(deletions_py, rs_py);
	 Py_XDECREF(str);
	 Py_XDECREF(rs_py);
      }
      r = PyList_New(3);
      // These are reversed so that the residue numbers come out in
      // numerical order (not backwards) and the returned list is
      // [mutations, deletions, insertions].
      PyList_SetItem(r, 0, mutations_py);
      PyList_SetItem(r, 1, deletions_py);
      PyList_SetItem(r, 2, insertions_py);
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON


// ----------------------------------------------------------------------------------
//                            LSQ
// ----------------------------------------------------------------------------------

void setup_lsq_deviation(int state) {

   graphics_info_t::in_lsq_plane_deviation = state;
} 

void setup_lsq_plane_define(int state) {

   std::cout << "DEBUG:: in_lsq_plane_define " << state << std::endl;
   graphics_info_t::in_lsq_plane_define = state;

} 

GtkWidget *wrapped_create_lsq_plane_dialog() {

   graphics_info_t g;
   return g.wrapped_create_lsq_plane_dialog();

}

void unset_lsq_plane_dialog() /* callback from destroy of widget */
{

   graphics_info_t::lsq_plane_dialog = 0;
   graphics_info_t::lsq_plane_atom_positions->resize(0);
} 

void remove_last_lsq_plane_atom() {

   graphics_info_t g;
   g.remove_last_lsq_plane_atom();
} 


void set_interactive_probe_dots_molprobity_radius(float r) { 

  graphics_info_t::probe_dots_on_chis_molprobity_radius = r;

} 

/*! \brief return the radius over which we can run interactive probe.
*/
float interactive_probe_dots_molprobity_radius() { 
  return graphics_info_t::probe_dots_on_chis_molprobity_radius; 
} 


void gln_asn_b_factor_outliers(int imol) {

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::use_graphics_interface_flag) { 
	 std::vector<std::pair<coot::atom_spec_t, std::string> > v = 
	    coot::util::gln_asn_b_factor_outliers(graphics_info_t::molecules[imol].atom_sel.mol);
	 
	 std::cout << "Found " << v.size() << " GLN/ASN B-factor outliers" << std::endl;
	 if (v.size() > 0) {
	    // c-interface-preferences unformatted dots reader had
	    // interesting code to make an interesting-places gui:
	    for (int i=0; i<v.size(); i++) {
	       std::cout << v[i].second << std::endl;
	    }
#if defined USE_GUILE && !defined WINDOWS_MINGW
	    graphics_info_t g;
	    std::vector<coot::util::atom_spec_and_button_info_t> outlier_atoms;
	    for (int i=0; i<v.size(); i++) {
	       std::string callback_func = "(lambda() (do-180-degree-side-chain-flip ";
	       callback_func += coot::util::int_to_string(imol);
	       callback_func += " ";
	       callback_func += single_quote(v[i].first.chain);
	       callback_func += " ";
	       callback_func += coot::util::int_to_string(v[i].first.resno);
	       callback_func += " ";
	       callback_func += single_quote(v[i].first.insertion_code);
	       callback_func += " ";
	       callback_func += single_quote(v[i].first.alt_conf);
	       callback_func += "))";
	       v[i].first.int_user_data = imol; // kludge in the imol, used for callback.
	       coot::util::atom_spec_and_button_info_t asi(v[i].first, v[i].second, callback_func);
	       outlier_atoms.push_back(asi);
	    }
	    std::string error_type = "Z score: ";
	    std::vector<std::string> cmd_strings;
	    cmd_strings.push_back("interesting-things-with-fix-maybe");
	    cmd_strings.push_back(single_quote("GLN and ASN B-factor Outliers"));
	    std::string ls = coot::util::interesting_things_list_with_fix(outlier_atoms, error_type);
	    cmd_strings.push_back(ls);
	    std::string s = g.state_command(cmd_strings, coot::STATE_SCM);
	    std::cout << "scheme command: " << s << std::endl;
	    safe_scheme_command(s);
#else
#ifdef USE_PYGTK
            graphics_info_t g;
            std::vector<coot::util::atom_spec_and_button_info_t> outlier_atoms;
            for (int i=0; i<v.size(); i++) {
               std::string callback_func = "[do_180_degree_side_chain_flip,";
               callback_func += coot::util::int_to_string(imol);
               callback_func += ",";
               callback_func += single_quote(v[i].first.chain);
               callback_func += ",";
               callback_func += coot::util::int_to_string(v[i].first.resno);
               callback_func += ",";
               callback_func += single_quote(v[i].first.insertion_code);
               callback_func += ",";
               callback_func += single_quote(v[i].first.alt_conf);
               callback_func += "]";
               v[i].first.int_user_data = imol; // kludge in the imol, used for callback.
               coot::util::atom_spec_and_button_info_t asi(v[i].first, v[i].second, callback_func);
               outlier_atoms.push_back(asi);
            }
            std::string error_type = "Z score: ";
            std::vector<std::string> cmd_strings;
            cmd_strings.push_back("interesting_things_with_fix_maybe");
            cmd_strings.push_back(single_quote("GLN and ASN B-factor Outliers"));
            std::string ls = coot::util::interesting_things_list_with_fix_py(outlier_atoms, error_type);
            cmd_strings.push_back(ls);
            std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
            std::cout << "python command: " << s << std::endl;
            safe_python_command(s);

#endif // PYGTK
#endif // USE_GUILE
	 } else {
	    std::string label = "Coot detected no GLN or ASN B-factor Outliers";
	    GtkWidget *w = wrapped_nothing_bad_dialog(label);
	    gtk_widget_show(w);
	 }
      }
   } 
} 

#ifdef USE_PYTHON
void gln_asn_b_factor_outliers_py(int imol) {

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::use_graphics_interface_flag) { 
	 std::vector<std::pair<coot::atom_spec_t, std::string> > v = 
	    coot::util::gln_asn_b_factor_outliers(graphics_info_t::molecules[imol].atom_sel.mol);
	 
	 std::cout << "Found " << v.size() << " GLN/ASN B-factor outliers" << std::endl;
	 if (v.size() > 0) {
	    // c-interface-preferences unformatted dots reader had
	    // interesting code to make an interesting-places gui:
	    for (int i=0; i<v.size(); i++) {
	       std::cout << v[i].second << std::endl;
	    }
#ifdef USE_PYGTK
            graphics_info_t g;
            std::vector<coot::util::atom_spec_and_button_info_t> outlier_atoms;
            for (int i=0; i<v.size(); i++) {
               std::string callback_func = "[do_180_degree_side_chain_flip,";
               callback_func += coot::util::int_to_string(imol);
               callback_func += ",";
               callback_func += single_quote(v[i].first.chain);
               callback_func += ",";
               callback_func += coot::util::int_to_string(v[i].first.resno);
               callback_func += ",";
               callback_func += single_quote(v[i].first.insertion_code);
               callback_func += ",";
               callback_func += single_quote(v[i].first.alt_conf);
               callback_func += "]";
               v[i].first.int_user_data = imol; // kludge in the imol, used for callback.
               coot::util::atom_spec_and_button_info_t asi(v[i].first, v[i].second, callback_func);
               outlier_atoms.push_back(asi);
            }
            std::string error_type = "Z score: ";
            std::vector<std::string> cmd_strings;
            cmd_strings.push_back("interesting_things_with_fix_maybe");
            cmd_strings.push_back(single_quote("GLN and ASN B-factor Outliers"));
            std::string ls = coot::util::interesting_things_list_with_fix_py(outlier_atoms, error_type);
            cmd_strings.push_back(ls);
            std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
            std::cout << "python command: " << s << std::endl;
            safe_python_command(s);
#endif // PYGTK

	 } else {
	    std::string label = "Coot detected no GLN or ASN B-factor Outliers";
	    GtkWidget *w = wrapped_nothing_bad_dialog(label);
	    gtk_widget_show(w);
	 }
      }
   } 
} 
#endif // PYTHON

/*! \brief For experienced Cooters who don't like Coot nannying about
  chiral volumes during refinement. */
void set_show_chiral_volume_errors_dialog(short int istate) {

   graphics_info_t::show_chiral_volume_errors_dialog_flag = istate;

}


/* does this live here really? */

#ifdef USE_GUILE
SCM get_torsion_scm(int imol, SCM atom_spec_1, SCM atom_spec_2, SCM atom_spec_3, SCM atom_spec_4) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as1 = atom_spec_from_scm_expression(atom_spec_1);
      coot::atom_spec_t as2 = atom_spec_from_scm_expression(atom_spec_2);
      coot::atom_spec_t as3 = atom_spec_from_scm_expression(atom_spec_3);
      coot::atom_spec_t as4 = atom_spec_from_scm_expression(atom_spec_4);

      graphics_info_t g;
      bool suc = g.set_angle_tors(imol, as1, as2, as3, as4); // just set a return value (to
							     // be used later, don't move
							     // the coordinates).
      if (suc) { 
	 double tors = g.get_geometry_torsion();
	 r = scm_double2num(tors);
      } else {
	 std::cout << "   WARNING:: (some) atoms not found in molecule #"
		   << imol << " "
		   << as1 << " " 
		   << as2 << " "
		   << as3 << " "
		   << as4 << std::endl;
      }
   }
   return r;
} 
#endif /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *get_torsion_py(int imol, PyObject *atom_spec_1, PyObject *atom_spec_2, PyObject *atom_spec_3, PyObject *atom_spec_4) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as1 = atom_spec_from_python_expression(atom_spec_1);
      coot::atom_spec_t as2 = atom_spec_from_python_expression(atom_spec_2);
      coot::atom_spec_t as3 = atom_spec_from_python_expression(atom_spec_3);
      coot::atom_spec_t as4 = atom_spec_from_python_expression(atom_spec_4);

      graphics_info_t g;
      bool suc = g.set_angle_tors(imol, as1, as2, as3, as4);
      if (suc) { 
	 double tors = g.get_geometry_torsion();
	 r = PyFloat_FromDouble(tors);
      } else {
	 std::cout << "   WARNING:: (some) atoms not found in molecule #"
		   << imol << " "
		   << as1 << " " 
		   << as2 << " "
		   << as3 << " "
		   << as4 << std::endl;
	    } 
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif /* USE_PYTHON */

// This is model manipulation, it does not belong here.
//
// tors is in degrees.
#ifdef USE_GUILE
SCM set_torsion_scm(int imol, const char *chain_id, int res_no, const char *insertion_code,
		    const char *alt_conf,
		    const char *atom_name_1,
		    const char *atom_name_2,
		    const char *atom_name_3,
		    const char *atom_name_4,
		    double tors) { 

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      double new_tors = g.molecules[imol].set_torsion(chain_id, res_no,
						      insertion_code, alt_conf,
						      atom_name_1,
						      atom_name_2,
						      atom_name_3,
						      atom_name_4, tors,
						      *g.Geom_p());
      r = scm_double2num(new_tors);
   }
   return r;
}
#endif // USE_GUILE

// This is model manipulation, it does not belong here (as above).
//
// tors is in degrees.
#ifdef USE_PYTHON
PyObject *set_torsion_py(int imol, const char *chain_id, int res_no, const char *insertion_code,
                         const char *alt_conf,
                         const char *atom_name_1,
                         const char *atom_name_2,
                         const char *atom_name_3,
                         const char *atom_name_4,
                         double tors) { 

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      double new_tors = g.molecules[imol].set_torsion(chain_id, res_no,
						      insertion_code, alt_conf,
						      atom_name_1,
						      atom_name_2,
						      atom_name_3,
						      atom_name_4, tors,
						      *g.Geom_p());
      r = PyFloat_FromDouble(new_tors);
   }
   if (PyBool_Check(r))
     Py_XINCREF(r);
   return r;
}
#endif // USE_PYTHON


#ifdef USE_GUILE
/*! \brief return the parsed user mod fields from the PDB file
  file_name (output by reduce most likely) */
SCM user_mods_scm(const char *file_name) {

   coot::flips_container f(file_name);
   return f.user_mods();
}

#endif // USE_GUILE
#ifdef USE_PYTHON
PyObject *user_mods_py(const char *file_name) {

   coot::flips_container f(file_name);

   return f.user_mods_py();
}
#endif // USE_PYTHON

