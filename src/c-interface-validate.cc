
/* src/c-interface-validate.cc
 *
 * Copyright 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009, 2010 The University of Oxford
 * Author: Paul Emsley
 * Copyright 2006, 2007 by Bernhard Lohkamp
 * Copyright 2015 by Medical Research Council
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"

#if defined _MSC_VER
#include <windows.h>
#endif

#include <stdlib.h>
#include <iostream>

#include "globjects.h" //includes gtk/gtk.h

#include <vector>
#include <string>

#include "clipper/ccp4/ccp4_mtz_io.h" // pathology plots

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"

#include "graphics-info.h"

#ifdef USE_GUILE
#include <libguile.h>
#endif // USE_GUILE

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: Python.h needs to come before to stop"_POSIX_C_SOURCE" redefined problems
//
// #ifdef USE_PYTHON
// #include "Python.h"
// #endif // USE_PYTHON

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "ligand/ligand.hh"

#include "coot-utils/peak-search.hh"
#include "user-mods.hh"

#include "data-pair-remover.hh"

#include "c-interface-gui.hh"
#include "widget-headers.hh"

#include "widget-from-builder.hh"

#include "rama_plot_with_canvas.hh"

#include "utils/logging.hh"
extern logging logger;


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

   // GtkWidget *dialog = create_check_waters_dialog();
   GtkWidget *dialog = widget_from_builder("check_waters_dialog");

   set_transient_and_position(COOT_CHECK_WATERS_DIALOG, dialog);

   // Opps - this (logical OR) should be on by default:
   GtkWidget *check_waters_OR_radiobutton  = widget_from_builder("check_waters_OR_radiobutton");

   gtk_check_button_set_active(GTK_CHECK_BUTTON(check_waters_OR_radiobutton), TRUE);

   GtkWidget *check_waters_action_combobox = widget_from_builder("check_waters_action_comboboxtext");

   if (check_waters_action_combobox) {
      gtk_combo_box_set_active(GTK_COMBO_BOX(check_waters_action_combobox), 0); // "Check"
   }

   GCallback callback_func = G_CALLBACK(nullptr);

   GtkWidget *combobox_molecule = widget_from_builder("check_waters_molecule_comboboxtext");

   gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox_molecule));

   // now fill that dialog's optionmenu with coordinate options.
   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 graphics_info_t::check_waters_molecule = imol;
	 break;
      }
   }

   graphics_info_t g;
   if (combobox_molecule)
      g.fill_combobox_with_coordinates_options(combobox_molecule, callback_func, g.check_waters_molecule);

   GtkWidget *entry;
   // char text[100];
   std::string text_str;

   // b-factor
   entry = widget_from_builder("check_waters_b_factor_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_b_factor_limit);
   gtk_editable_set_text(GTK_EDITABLE(entry), text_str.c_str());


   // map sigma
   entry = widget_from_builder("check_waters_map_sigma_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_map_sigma_limit);
   gtk_editable_set_text(GTK_EDITABLE(entry), text_str.c_str());

   // min_dist
   entry = widget_from_builder("check_waters_min_dist_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_min_dist_limit);
   gtk_editable_set_text(GTK_EDITABLE(entry), text_str.c_str());

   // max_dist
   entry = widget_from_builder("check_waters_max_dist_entry");
   text_str = graphics_info_t::float_to_string(graphics_info_t::check_waters_max_dist_limit);
   gtk_editable_set_text(GTK_EDITABLE(entry), text_str.c_str());

   // 20100131 We have put the variance map check into this dialog
   // too, to better organise the menus.
   //
   //GtkWidget *diff_map_option_menu =
   // widget_from_builder("check_water_by_difference_map_optionmenu");

   GtkWidget *diff_map_combobox = widget_from_builder("check_waters_by_difference_map_combobox");

   if (diff_map_combobox) {
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
	 GCallback signal_func = G_CALLBACK(check_water_by_difference_maps_combobox_changed);
	 g.fill_combobox_with_difference_map_options(diff_map_combobox, signal_func, imol_active);
      }
   }

   return dialog;

}

void
check_water_by_difference_maps_combobox_changed(GtkWidget *combobox, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::check_waters_by_difference_map_map_number = imol;

}


// The OK button was pressed on the dialog, so read the dialog and do
// the check.
//
// called by a callbacks.c function.
//
void do_check_waters_by_widget(GtkWidget *dialog) {

   GtkWidget *checklogic_OR_radiobutton  = widget_from_builder("check_waters_OR_radiobutton");

   GtkWidget *entry1, *entry2, *entry3, *entry4;
   entry1 = widget_from_builder("check_waters_b_factor_entry");
   entry2 = widget_from_builder("check_waters_map_sigma_entry");
   entry3 = widget_from_builder("check_waters_min_dist_entry");
   entry4 = widget_from_builder("check_waters_max_dist_entry");

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

   GtkWidget *hbox1 = widget_from_builder("check_waters_b_factor_hbox");
   GtkWidget *hbox2 = widget_from_builder("check_waters_sigma_level_hbox");
   GtkWidget *hbox3 = widget_from_builder("check_waters_min_dist_hbox");
   GtkWidget *hbox4 = widget_from_builder("check_waters_max_dist_hbox");

   GtkCheckButton *checkbutton1 = GTK_CHECK_BUTTON(widget_from_builder("check_waters_b_factor_entry_active_checkbutton"));
   GtkCheckButton *checkbutton2 = GTK_CHECK_BUTTON(widget_from_builder("check_waters_map_sigma_entry_active_checkbutton"));
   GtkCheckButton *checkbutton3 = GTK_CHECK_BUTTON(widget_from_builder("check_waters_min_dist_entry_active_checkbutton"));
   GtkCheckButton *checkbutton4 = GTK_CHECK_BUTTON(widget_from_builder("check_waters_max_dist_entry_active_checkbutton"));
   GtkCheckButton *checkbutton5 = GTK_CHECK_BUTTON(widget_from_builder("check_waters_by_difference_map_active_checkbutton"));

   if (! gtk_check_button_get_active(checkbutton1)) use_b_factor_limit_test  = 0;
   if (! gtk_check_button_get_active(checkbutton2)) use_map_sigma_limit_test = 0;
   if (! gtk_check_button_get_active(checkbutton3)) use_min_dist_test        = 0;
   if (! gtk_check_button_get_active(checkbutton4)) use_max_dist_test        = 0;

   if (checkbutton5)
      if (! gtk_check_button_get_active(checkbutton5))
	 use_difference_map_test = 0;

   GtkWidget *zero_occ_checkbutton = widget_from_builder("check_waters_zero_occ_checkbutton");
   GtkWidget *partial_occ_close_contact_checkbutton =
      widget_from_builder("check_waters_low_occ_dist_checkbutton");

   short int zero_occ_flag = 0;
   short int part_occ_dist_flag = 0;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(zero_occ_checkbutton)))
      zero_occ_flag = 1;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(partial_occ_close_contact_checkbutton)))
      part_occ_dist_flag = 1;

   //
   short int logical_operator_and_or_flag = 0; // logical AND
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checklogic_OR_radiobutton))) {
      logical_operator_and_or_flag = 1;
   }

   // Check or Delete?
   GtkWidget *action_combobox = widget_from_builder("check_waters_action_comboboxtext");
   std::string action_string;
   gchar *at = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(action_combobox));
   if (at) action_string = at;
   if (!at) std::cout << "ERROR: null from action combobox " << action_combobox << std::endl;

   // This will give us another dialog
   //
   if (use_difference_map_test) {
      int imol_diff_map = graphics_info_t::check_waters_by_difference_map_map_number;
      check_waters_by_difference_map(graphics_info_t::check_waters_molecule, imol_diff_map, 1);
   }

   GtkWidget *model_combobox = widget_from_builder("check_waters_molecule_comboboxtext");
   int imol_check_waters = my_combobox_get_imol(GTK_COMBO_BOX(model_combobox));

   if (use_b_factor_limit_test == 0)
      b_factor_lim = -100.0;
   if (use_map_sigma_limit_test == 0)
      map_sigma_lim = -100.0;
   if (use_min_dist_test == 0)
      min_dist = -100.0;
   if (use_max_dist_test == 0)
      max_dist = -100.0;  // sets a flag in find_water_baddies_OR
   if (action_string == "Check") {
      GtkWidget *w = wrapped_checked_waters_baddies_dialog(imol_check_waters,
							   b_factor_lim,
							   map_sigma_lim,
							   min_dist,
							   max_dist,
							   part_occ_dist_flag,
							   zero_occ_flag,
							   logical_operator_and_or_flag);
      set_transient_for_main_window(w);
      gtk_widget_set_visible(w, TRUE);
   }


   if (action_string == "Delete") {

      // delete those baddies:
      delete_checked_waters_baddies(imol_check_waters,
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



void check_waters_molecule_combobox_changed(GtkWidget *combobox,
					  gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::check_waters_molecule = imol;
}


std::vector<coot::atom_spec_t>
check_waters_baddies(int imol, float b_factor_lim, float map_sigma_lim, float min_dist, float max_dist, short int part_occ_contact_flag, short int zero_occ_flag, short int logical_operator_and_or_flag) {

   std::vector<coot::atom_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_for_map = g.Imol_Refinement_Map();
      v = g.molecules[imol].find_water_baddies(b_factor_lim, g.molecules[imol_for_map].xmap, g.molecules[imol_for_map].map_sigma(), map_sigma_lim, min_dist, max_dist, part_occ_contact_flag, zero_occ_flag, logical_operator_and_or_flag);
      if (graphics_info_t::use_graphics_interface_flag) {
	 GtkWidget *w = wrapped_checked_waters_baddies_dialog(imol,
							      b_factor_lim,
							      map_sigma_lim,
							      min_dist,
							      max_dist,
							      part_occ_contact_flag,
							      zero_occ_flag,
							      logical_operator_and_or_flag);
	 gtk_widget_set_visible(w, TRUE);
      }
   }
   return v;
}


// On check OK, we fire up this widget which is a vbox of baddy radio buttons.
//
GtkWidget *wrapped_checked_waters_baddies_dialog(int imol, float b_factor_lim, float map_sigma_lim,
                                                 float min_dist, float max_dist,
                                                 short int part_occ_contact_flag, short int zero_occ_flag,
                                                 short int logical_operator_and_or_flag) {

   GtkWidget *w = NULL;
   if (graphics_info_t::use_graphics_interface_flag) {
      // w = create_checked_waters_baddies_dialog();
      w = widget_from_builder("checked_waters_baddies_dialog");

      set_transient_and_position(COOT_CHECKED_WATERS_BADDIES_DIALOG, w);

      graphics_info_t g;
      int imol_for_map = g.Imol_Refinement_Map();

      GtkWidget *button_group = nullptr;
      if (is_valid_model_molecule(imol)) {
	 if (!is_valid_map_molecule(imol_for_map)) {
	    std::cout << "WARNING:: Not a valid map for density testing "
		      << imol_for_map << std::endl;
	 } else {

	    std::vector<coot::atom_spec_t> baddies =
	       g.molecules[imol].find_water_baddies(b_factor_lim, graphics_info_t::molecules[imol_for_map].xmap, graphics_info_t::molecules[imol_for_map].map_sigma(), map_sigma_lim, min_dist, max_dist, part_occ_contact_flag, zero_occ_flag, logical_operator_and_or_flag);

	    // User data is used to keyboard up and down baddie water
	    // list (in graphics_info_t::checked_waters_next_baddie).
	    g_object_set_data(G_OBJECT(w), "baddies_size", GINT_TO_POINTER(baddies.size()));

	    GtkWidget *vbox = widget_from_builder("checked_waters_baddies_vbox");

            g.clear_out_container(vbox);

	    if (baddies.size() > 0 ) {
	       for (int i=0; i<int(baddies.size()); i++) {

		  // 	       std::cout << "Suspicious water: "
		  // 			 << baddies[i].atom_name
		  // 			 << baddies[i].alt_conf << " "
		  // 			 << baddies[i].resno << " "
		  // 			 << baddies[i].insertion_code << " "
		  // 			 << baddies[i].chain << "\n";

		  std::string button_label(" ");
		  button_label += baddies[i].chain_id;
		  button_label += " " ;
		  button_label += graphics_info_t::int_to_string(baddies[i].res_no);
		  button_label += " " ;
		  button_label += baddies[i].atom_name;
		  button_label += " " ;
		  button_label += baddies[i].alt_conf;
		  button_label += " [Occ: " ;
		  button_label += graphics_info_t::float_to_string(baddies[i].float_user_data);
		  button_label += "]   " ;
		  button_label += baddies[i].string_user_data;
		  button_label += " " ;

		  GtkWidget *toggle_button = gtk_toggle_button_new_with_label(button_label.c_str());
                  gtk_widget_set_margin_top(toggle_button, 2);
                  gtk_widget_set_margin_bottom(toggle_button, 2);
                  gtk_widget_set_margin_start(toggle_button, 6);
                  gtk_widget_set_margin_end(toggle_button, 6);

                  // set the group here.
                  if (button_group)
                     gtk_toggle_button_set_group(GTK_TOGGLE_BUTTON(toggle_button), GTK_TOGGLE_BUTTON(button_group));
                  else
                     button_group = toggle_button;

		  coot::atom_spec_t *atom_spec = new coot::atom_spec_t(baddies[i]);
		  atom_spec->int_user_data = imol;

		  std::string button_name = "checked_waters_baddie_button_";
		  button_name += coot::util::int_to_string(i);

		  g_signal_connect(G_OBJECT(toggle_button), "toggled",
				   G_CALLBACK(graphics_info_t::on_generic_atom_spec_toggle_button_toggled),
				   atom_spec);
		  gtk_box_append(GTK_BOX(vbox), toggle_button);
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
	 show_select_map_frame();
      } else {
	 std::vector<coot::atom_spec_t> baddies =
	    graphics_info_t::molecules[imol].find_water_baddies(b_factor_lim,
								graphics_info_t::molecules[imol_for_map].xmap,
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
	    gtk_widget_set_visible(w, TRUE);
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

#ifdef USE_GUILE
SCM chiral_volume_errors_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      r = SCM_EOL;
      graphics_info_t g;
      std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > v = g.molecules[imol].bad_chiral_volumes();
      for (std::size_t i=0; i<v.second.size(); i++) {
	 SCM atom_spec_scm = atom_spec_to_scm(v.second[i]);
	 r = scm_cons(atom_spec_scm, r);
      }
      r = scm_reverse(r);
   }
   return r;

}
#endif /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *chiral_volume_errors_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > v = g.molecules[imol].bad_chiral_volumes();
      r = PyList_New(v.second.size());
      for (std::size_t i=0; i<v.second.size(); i++) {
	 PyObject *atom_spec_py = atom_spec_to_py(v.second[i]);
	 PyList_SetItem(r, i, atom_spec_py);
      }
   }
   if (PyBool_Check(r))
     Py_INCREF(r);
   return r;
}
#endif	/* USE_PYTHON */




void set_fix_chiral_volumes_before_refinement(int istate) {
   graphics_info_t::fix_chiral_volume_before_refinement_flag = istate;
}


void check_chiral_volumes_from_widget(GtkWidget *window) {

   check_chiral_volumes(graphics_info_t::check_chiral_volume_molecule);
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

   std::cout << "GTK-FIXME ------------ free_geometry_graph() " << std::endl;
#if 0

   if (dialog) {
      GtkWidget *w = widget_from_builder("geometry_graph_canvas");
      if (w) {
	 GtkObject *obj = GTK_OBJECT(w);

         // Note to self: FIXME - where is this set?
	 // coot::geometry_graphs *graphs = (coot::geometry_graphs *) gtk_object_get_user_data(obj);
         GObject *o = g_object_get_data(obj, "graphs");
         if (o) {
            coot::geometry_graphs *graphs = static_cast<coot::geometry_graphs *> (o);
            if (!graphs) {
               std::cout << "ERROR:: NULL graphs in free_geometry_graph\n";
            } else {
               graphs->close_yourself();
            }
         } else {
            std::cout << "ERROR:: failed convert old gtk lookup to new gobject lookup " << std::endl;
         }

      }
   }

#endif
}

void unset_geometry_graph(GtkWidget *dialog) {  /* set the graphics info
						 static to NULL, so
						 that we on longer try
						 to update the
						 widget*/

   // kill this off
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
	 SCM name = scm_from_locale_string(results.info[ir].rotamer_name.c_str());
	 ele = scm_cons(name, ele);
	 SCM pr = scm_from_double(results.info[ir].probability);
	 ele = scm_cons(pr, ele);
	 SCM inscode = scm_from_locale_string(results.info[ir].inscode.c_str());
	 ele = scm_cons(inscode, ele);
	 SCM resno = scm_from_int(results.info[ir].resno);
	 ele = scm_cons(resno, ele);
	 SCM chainid = scm_from_locale_string(results.info[ir].chain_id.c_str());
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
   if (! results.info.empty()) {
      r = PyList_New(results.info.size());
      for (int ir=int(results.info.size()-1); ir>=0; ir--) {
	 PyObject *ele = PyList_New(5);
	 PyObject *name = myPyString_FromString(results.info[ir].rotamer_name.c_str());
	 PyList_SetItem(ele, 4, name);;
	 PyObject *pr = PyFloat_FromDouble(results.info[ir].probability);
	 PyList_SetItem(ele, 3, pr);
	 PyObject *inscode = myPyString_FromString(results.info[ir].inscode.c_str());
	 PyList_SetItem(ele, 2, inscode);
	 PyObject *resno = PyLong_FromLong(results.info[ir].resno);
	 PyList_SetItem(ele, 1, resno);
	 PyObject *chainid = myPyString_FromString(results.info[ir].chain_id.c_str());
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


void
my_delete_validaton_graph_mol_option(GtkWidget *widget, void *data) {
   // gtk_container_remove(GTK_CONTAINER(data), widget);
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
		               float max_closeness,
		               int do_positive_level_flag,
		               int do_negative_levels_flag,
                               int around_model_only_flag) {

   // Notice that we make wrapped_create_check_waters_dialog be part
   // of graphics_info_t, because it uses clipper data in the
   // interface - and c-interface.h does not know about clipper.

   // I don't think we want ligand/cluster search.  We just want peak
   // searching.
   //
   if (! is_valid_map_molecule(imol)) {
      // std::cout << "Molecule number " << imol << " is not a valid map molecule" << std::endl;
      logger.log(log_t::WARNING, "Molecule number", imol, "is not a valid molecule");
      return;
   }

   if (!graphics_info_t::molecules[imol].is_difference_map_p()) {
      return;
   }

   // c.f. trace-high-res.cc
   coot::peak_search ps(graphics_info_t::molecules[imol].xmap);
   ps.set_max_closeness(max_closeness);
   std::vector<std::pair<clipper::Coord_orth, float> > centres;

   if (is_valid_model_molecule(imol_coords)) {
      centres =
         ps.get_peaks(graphics_info_t::molecules[imol].xmap,
            graphics_info_t::molecules[imol_coords].atom_sel.mol,
            n_sigma, do_positive_level_flag, do_negative_levels_flag,
                           around_model_only_flag);
   } else {
      centres =
         ps.get_peaks(graphics_info_t::molecules[imol].xmap,
            n_sigma, do_positive_level_flag, do_negative_levels_flag);
   }

   if (centres.size() == 0) {
      if (graphics_info_t::use_graphics_interface_flag) {
         std::string info_string("No difference map peaks\nat ");
         info_string += graphics_info_t::float_to_string(n_sigma);
         info_string += " sigma";
         GtkWidget *w = wrapped_nothing_bad_dialog(info_string);
         gtk_widget_set_visible(w, TRUE);
      }
   } else {
      float map_sigma = graphics_info_t::molecules[imol].map_sigma();
      if (graphics_info_t::use_graphics_interface_flag) {
         graphics_info_t::show_diff_map_peaks_vbox(imol, imol_coords, centres,
                                                               n_sigma,
                                                               do_positive_level_flag,
                                                               do_negative_levels_flag,
                                                               around_model_only_flag);
         GtkWidget *peaks_vbox = widget_from_builder("diff_map_peaks_vbox");
         gtk_widget_set_visible(peaks_vbox, TRUE);
      }

      // std::cout << "\nFound these peak positions:\n";
      logger.log(log_t::INFO, logging::function_name_t("difference_map_peaks"),
		 "Found these peak positions");
      for (unsigned int i=0; i<centres.size(); i++) {
         // std::cout << "   " << i << " dv: "
	 // 	   << centres[i].second << " n-rmsd: "
	 // 	   << centres[i].second/map_sigma << " "
	 // 	   << centres[i].first.format() << std::endl;
	 auto dv  = centres[i].second;
	 auto pos = centres[i].first.format();
	 auto nrmsd = dv/map_sigma;
	 logger.log(log_t::INFO, {"      ", i, "dv", dv, "n-rmsd:", nrmsd, "at", pos});
      }
      // std::cout << "\n   Found " << centres.size() << " peak positions:\n";
      logger.log(log_t::INFO, "Found", std::to_string(centres.size()), "peak positions");
   }
}

void set_difference_map_peaks_widget(GtkWidget *w) {
   graphics_info_t::difference_map_peaks_dialog = w;
}

// In the GUI for difference map peaks, there is not a means to set the
// max_closeness, so here is a means to set it and query it. */
void set_difference_map_peaks_max_closeness(float m) {
   graphics_info_t::difference_map_peaks_max_closeness = m;
}

float difference_map_peaks_max_closeness() {
   return graphics_info_t::difference_map_peaks_max_closeness;
}



void
clear_diff_map_peaks() {

   graphics_info_t g;
   g.clear_diff_map_peaks();
}


GtkWidget *wrapped_create_generate_diff_map_peaks_dialog() {

   // c.f. wrapped_create_check_waters_diff_map_dialog()
   // GtkWidget *dialog = widget_from_builder("generate_diff_map_peaks_dialog");
   // short int diff_maps_only_flag = 1;
   // int ifound = fill_ligands_dialog_map_bits_by_dialog_name(dialog, "generate_diff_map_peaks_map",
   //                                                         diff_maps_only_flag);
   // if (ifound == 0) {
   //    std::cout << "Error: you must have a difference map to analyse!" << std::endl;
   //    GtkWidget *none_frame = widget_from_builder("no_difference_maps_frame");
   //    gtk_widget_set_visible(none_frame, TRUE);
   // }
   // the name of the vbox which is looked up is "generate_diff_map_peaks_model_vbox".
   // ifound = fill_ligands_dialog_protein_bits_by_dialog_name(dialog, "generate_diff_map_peaks_model");
   // if (ifound == 0) {
   //    std::cout << "Difference map checker is better having specified coordinates...\n";
   // }


   graphics_info_t g;
   GtkWidget *dialog = widget_from_builder("generate_diff_map_peaks_dialog");

   GtkWidget *model_combobox = widget_from_builder("generate_diff_map_peaks_molecule_combobox");
   GtkWidget   *map_combobox = widget_from_builder("generate_diff_map_peaks_map_combobox");
   GtkWidget *frame_1 = widget_from_builder("no_difference_maps_frame1");
   GtkWidget *frame_2 = widget_from_builder("generate_diff_maps_peaks_no_models_frame");

   auto get_model_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_model_molecule(i))
                                           vec.push_back(i);
                                     return vec;
                                  };

   auto get_map_molecule_vector = [] () {
                                     graphics_info_t g;
                                     std::vector<int> vec;
                                     int n_mol = g.n_molecules();
                                     for (int i=0; i<n_mol; i++)
                                        if (g.is_valid_map_molecule(i))
                                           if (g.molecules[i].xmap_is_diff_map)
                                             vec.push_back(i);
                                     return vec;
                                  };
   int imol_active = -1;
   int imol_active_map = -1;
   GCallback func = G_CALLBACK(nullptr); // we don't care until this dialog is read
   auto model_list = get_model_molecule_vector();
   auto map_list   =   get_map_molecule_vector();
   g.fill_combobox_with_molecule_options(model_combobox, func, imol_active,  model_list);
   g.fill_combobox_with_molecule_options(  map_combobox, func, imol_active_map, map_list);

   // std::cout << "::::::::   map_list size " <<   map_list.size() << std::endl;
   // std::cout << ":::::::: model_list size " << model_list.size() << std::endl;

   if (model_list.empty()) {
      gtk_widget_set_visible(model_combobox, FALSE);
      gtk_widget_set_visible(frame_2, TRUE);
   } else {
      gtk_widget_set_visible(model_combobox, TRUE);
      gtk_widget_set_visible(frame_2, FALSE);
   }

   if (map_list.empty()) {
      gtk_widget_set_visible(map_combobox, FALSE);
      gtk_widget_set_visible(frame_1, TRUE);
   } else {
      gtk_widget_set_visible(map_combobox, TRUE);
      gtk_widget_set_visible(frame_1, FALSE);
   }

   // the sigma entry:
   GtkWidget *entry = widget_from_builder("generate_diff_map_peaks_sigma_level_entry");
   std::string s = g.float_to_string(g.difference_map_peaks_sigma_level);
   gtk_editable_set_text(GTK_EDITABLE(entry), s.c_str());

   return dialog;
}

void difference_map_peaks_from_dialog() {

   // Check the level:

   GtkWidget *sigma_entry = widget_from_builder("generate_diff_map_peaks_sigma_level_entry");
   const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(sigma_entry));
   float v = coot::util::string_to_float(txt);
   bool good_sigma = false;
   if (v > -1000 && v < 1000) {
      good_sigma = true;
   } else {
      std::cout << "WARNING:: Invalid sigma level: " << v << " can't do peak search." << std::endl;
   }


   bool do_negative_level = false;
   bool do_positive_level = false;
   GtkWidget *checkbutton_negative = widget_from_builder("generate_diff_map_peaks_negative_level_checkbutton");
   GtkWidget *checkbutton_positive = widget_from_builder("generate_diff_map_peaks_positive_level_checkbutton");
   //
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton_negative)))
      do_negative_level = true;

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton_positive)))
      do_positive_level = true;

   GtkWidget *model_combobox = widget_from_builder("generate_diff_map_peaks_molecule_combobox");
   GtkWidget   *map_combobox = widget_from_builder("generate_diff_map_peaks_map_combobox");

   int imol_coords   = my_combobox_get_imol(GTK_COMBO_BOX(model_combobox));
   int imol_diff_map = my_combobox_get_imol(GTK_COMBO_BOX(map_combobox));
   bool around_model_only = true;

   if (good_sigma)
      difference_map_peaks(imol_diff_map, imol_coords, v,
                           graphics_info_t::difference_map_peaks_max_closeness,
                           do_positive_level, do_negative_level, around_model_only);
}


#ifdef USE_PYTHON
PyObject *map_peaks_around_molecule_py(int imol_map, float n_sigma, int do_negative_also_flag, int imol_coords) {

   PyObject *r = Py_False;
   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_model_molecule(imol_coords)) {
	 const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
	 coot::peak_search ps(xmap);
	 ps.set_max_closeness(graphics_info_t::difference_map_peaks_max_closeness);
	 std::cout << "using max_closeness " << graphics_info_t::difference_map_peaks_max_closeness
		   << std::endl;
	 int do_positive_level_flag = 1;
	 std::cout << "getting centres with negative-flag " << do_negative_also_flag
		   << std::endl;
         int around_model_only_flag = false;
	 std::vector<std::pair<clipper::Coord_orth, float> > centres =
	    ps.get_peaks(graphics_info_t::molecules[imol_map].xmap,
			 graphics_info_t::molecules[imol_coords].atom_sel.mol,
			 n_sigma, do_positive_level_flag, do_negative_also_flag,
                         around_model_only_flag);
	 r = PyList_New(centres.size());
	 for (unsigned int i=0; i<centres.size(); i++) {
	    PyObject *coords = PyList_New(3);
	    PyObject *pair_py = PyList_New(2);
	    PyList_SetItem(coords, 0, PyFloat_FromDouble(centres[i].first.x()));
	    PyList_SetItem(coords, 1, PyFloat_FromDouble(centres[i].first.y()));
	    PyList_SetItem(coords, 2, PyFloat_FromDouble(centres[i].first.z()));
	    PyList_SetItem(pair_py, 0, PyFloat_FromDouble(centres[i].second));
	    PyList_SetItem(pair_py, 1, coords);
	    PyList_SetItem(r, i, pair_py);
	 }
      }
   }
   if (PyBool_Check(r))
     Py_INCREF(r);
   return r;
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
PyObject *map_peaks_py(int imol_map, float n_sigma) {

   PyObject *r = Py_False;

   if (is_valid_map_molecule(imol_map)) {
      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
      int do_positive_levels_flag = 1;
      int also_negative_levels_flag = 0;
      coot::peak_search ps(xmap);
      ps.set_max_closeness(0.0f);
      std::vector<std::pair<clipper::Coord_orth, float> > peaks =
	 ps.get_peaks(xmap, n_sigma, do_positive_levels_flag, also_negative_levels_flag);
      r = PyList_New(peaks.size());
      for (unsigned int i=0; i<peaks.size(); i++) {
	 PyObject *coords = PyList_New(3);
	 PyList_SetItem(coords, 0, PyFloat_FromDouble(peaks[i].first.x()));
	 PyList_SetItem(coords, 1, PyFloat_FromDouble(peaks[i].first.y()));
	 PyList_SetItem(coords, 2, PyFloat_FromDouble(peaks[i].first.z()));
	 PyList_SetItem(r, i, coords);
      }
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }

   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *map_peaks_near_point_py(int imol_map, float n_sigma, float x, float y, float z,
				  float radius) {

   PyObject *r = Py_False;

   if (is_valid_map_molecule(imol_map)) {

      mmdb::Atom *at = new mmdb::Atom;
      at->SetCoordinates(x,y,z, 1.0, 10.0);
      at->SetAtomName(" CA ");
      at->SetElementName(" C");

      graphics_info_t g;
      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_atom(at);
      mol->SetSpaceGroup(g.molecules[imol_map].xmap.spacegroup().symbol_hm().c_str());
      coot::util::set_mol_cell(mol, g.molecules[imol_map].xmap.cell());

      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
      int do_positive_levels_flag = 1;
      int also_negative_levels_flag = 0;
      coot::peak_search ps(xmap);
      int around_model_only_flag = 0;
      std::vector<std::pair<clipper::Coord_orth, float> > peaks =
	 ps.get_peaks(xmap, mol, n_sigma, do_positive_levels_flag, also_negative_levels_flag, around_model_only_flag);
      clipper::Coord_orth ref_pt(x,y,z);
      std::vector<std::pair<clipper::Coord_orth, float> > close_peaks;
      for (unsigned int i=0; i<peaks.size(); i++) {
	 if (clipper::Coord_orth::length(ref_pt, peaks[i].first) < radius) {
	    close_peaks.push_back(peaks[i]);
	 }
      }
      r = PyList_New(close_peaks.size());
      for (unsigned int i=0; i<close_peaks.size(); i++) {
	 PyObject *coords = PyList_New(4);
	 PyList_SetItem(coords, 0, PyFloat_FromDouble(close_peaks[i].first.x()));
	 PyList_SetItem(coords, 1, PyFloat_FromDouble(close_peaks[i].first.y()));
	 PyList_SetItem(coords, 2, PyFloat_FromDouble(close_peaks[i].first.z()));
	 PyList_SetItem(coords, 3, PyFloat_FromDouble(close_peaks[i].second));
	 PyList_SetItem(r, i, coords);
      }
      delete mol;
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }

   return r;
}
#endif


#ifdef USE_GUILE
SCM map_peaks_scm(int imol_map, float n_sigma) {

   SCM r = SCM_BOOL_F;

   if (is_valid_map_molecule(imol_map)) {
      r = SCM_EOL;
      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
      int do_positive_levels_flag = 1;
      int also_negative_levels_flag = 0;
      coot::peak_search ps(xmap);
      ps.set_max_closeness(0.0f);
      std::vector<std::pair<clipper::Coord_orth, float> > peaks =
	 ps.get_peaks(xmap, n_sigma, do_positive_levels_flag, also_negative_levels_flag);
      for (unsigned int i=0; i<peaks.size(); i++) {
	 SCM pt = SCM_EOL;
	 SCM x_scm = scm_from_double(peaks[i].first.x());
	 SCM y_scm = scm_from_double(peaks[i].first.y());
	 SCM z_scm = scm_from_double(peaks[i].first.z());
	 pt = scm_cons(x_scm, pt);
	 pt = scm_cons(y_scm, pt);
	 pt = scm_cons(z_scm, pt);
	 pt = scm_reverse(pt);
	 r = scm_cons(pt, r);
      }
   }
   return r;
}
#endif

#ifdef USE_GUILE
SCM map_peaks_near_point_scm(int imol_map, float n_sigma, float x, float y, float z,
			     float radius) {

   SCM r = SCM_BOOL_F;
   if (is_valid_map_molecule(imol_map)) {

      graphics_info_t g;
      mmdb::Atom *at = new mmdb::Atom;
      at->SetCoordinates(x,y,z,1.0,10.0);
      at->SetAtomName(" CA ");
      at->SetElementName(" C");

      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_atom(at);
      mol->SetSpaceGroup(g.molecules[imol_map].xmap.spacegroup().symbol_hm().c_str());
      coot::util::set_mol_cell(mol, g.molecules[imol_map].xmap.cell());

      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
      int do_positive_levels_flag = 1;
      int also_negative_levels_flag = 0;
      coot::peak_search ps(xmap);
      int around_model_only_flag = 0;
      std::vector<std::pair<clipper::Coord_orth, float> > peaks =
	 ps.get_peaks(xmap, mol, n_sigma, do_positive_levels_flag, also_negative_levels_flag, around_model_only_flag);
      clipper::Coord_orth ref_pt(x,y,z);
      r = SCM_EOL;
      std::vector<std::pair<clipper::Coord_orth, float> > close_peaks;
      for (unsigned int i=0; i<peaks.size(); i++) {
	 if (clipper::Coord_orth::length(ref_pt, peaks[i].first) < radius) {
	    close_peaks.push_back(peaks[i]);
	 }
      }

      if (1) { // debug
	 for (unsigned int i=0; i<close_peaks.size(); i++) {
	    std::cout << "close peak " << i << " " << close_peaks[i].first.format() << "   "
		      << close_peaks[i].second << std::endl;
	 }
      }

      for (unsigned int i=0; i<close_peaks.size(); i++) {
	 SCM pt = SCM_EOL;
	 SCM d     = scm_from_double(close_peaks[i].second);
	 SCM x_scm = scm_from_double(close_peaks[i].first.x());
	 SCM y_scm = scm_from_double(close_peaks[i].first.y());
	 SCM z_scm = scm_from_double(close_peaks[i].first.z());
	 pt = scm_cons(d, pt);
	 pt = scm_cons(x_scm, pt);
	 pt = scm_cons(y_scm, pt);
	 pt = scm_cons(z_scm, pt);
	 pt = scm_reverse(pt);
	 r = scm_cons(pt, r);
      }
      r = scm_reverse(r);
      delete mol;
   }
   return r;
}
#endif
// (map-peaks-near-point-scm 1 4 44.7 11 13.6 6)


// ----------------------------------------------------------------------------------
//                            ramachandran plot
// ----------------------------------------------------------------------------------


// fill the widget:
GtkWidget *wrapped_ramachandran_plot_differences_dialog() {

   GtkWidget *w = 0; // Not NULL, compiler (maybe).

   std::cout << "ERROR:: using unconverted wrapped_ramachandran_plot_differences_dialog()" << std::endl;

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   // w = create_ramachandran_plot_differences_dialog();
   w = widget_from_builder("ramachandran_plot_differences_dialog");

   // We don't have to worry about chains because they are not active on startup.

   GtkWidget *combobox1 = widget_from_builder("ramachandran_plot_differences_first_mol_combobox");
   GtkWidget *combobox2 = widget_from_builder("ramachandran_plot_differences_second_mol_combobox");

   if (! combobox1) std::cout << "null combobox1" << std::endl;
   if (! combobox2) std::cout << "null combobox2" << std::endl;

   GCallback signal_func1 = G_CALLBACK(ramachandran_plot_differences_mol_combobox_first_changed);
   GCallback signal_func2 = G_CALLBACK(ramachandran_plot_differences_mol_combobox_second_changed);

   int imol = -1;
   for (int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }

   if (imol >= 0) {
      graphics_info_t g;
      g.fill_combobox_with_coordinates_options(combobox1, signal_func1, imol);
      g.fill_combobox_with_coordinates_options(combobox2, signal_func2, imol);
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

   GtkWidget *checkbutton1 = widget_from_builder("ramachandran_plot_differences_first_chain_checkbutton");
   GtkWidget *checkbutton2 = widget_from_builder("ramachandran_plot_differences_second_chain_checkbutton");

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton1)) &&
       gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton2))) {
      istat = 1;
      ramachandran_plot_differences_by_chain(imol1, imol2,
					     first_chain.c_str(),
					     second_chain.c_str());
   } else {
      if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton1)) &&
	  !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton2))) {
	 istat = 1;
	 ramachandran_plot_differences(imol1, imol2);
      } else {
	 std::cout << "INFO:: incomprehensible molecule/chain selection" << std::endl;
	 std::string s = "Can't make sense of chain selection.  Try again?";
	 GtkWidget *nbd = wrapped_nothing_bad_dialog(s);
	 gtk_widget_set_visible(nbd, TRUE);
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



void ramachandran_plot_differences_mol_combobox_first_changed(GtkWidget *combobox, gpointer pos) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::ramachandran_plot_differences_imol1 = imol;

   GtkWidget *chain_combobox = widget_from_builder("ramachandran_plot_differences_first_chain_combobox");
   GtkWidget *checkbutton    = widget_from_builder("ramachandran_plot_differences_first_chain_checkbutton");
   if (!chain_combobox) {
      std::cout << "first bad combobox" << std::endl;
   } else {
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) {
	 fill_ramachandran_plot_differences_combobox_with_chain_options(chain_combobox, true);
      }
   }
}

void ramachandran_plot_differences_mol_combobox_second_changed(GtkWidget *combobox, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::ramachandran_plot_differences_imol2 = imol;

   GtkWidget *chain_combobox = widget_from_builder("ramachandran_plot_differences_second_chain_combobox");
   GtkWidget *checkbutton    = widget_from_builder("ramachandran_plot_differences_second_chain_checkbutton");
   if (!chain_combobox) {
      std::cout << "first bad combobox" << std::endl;
   } else {
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) {
	 fill_ramachandran_plot_differences_combobox_with_chain_options(chain_combobox, false);
      }
   }

}

void ramachandran_plot_differences_chain_combobox_first_changed(GtkWidget *combobox, gpointer data) {

   std::string label = get_active_label_in_combobox(GTK_COMBO_BOX(combobox));
   graphics_info_t::ramachandran_plot_differences_imol1_chain = label;
   std::cout << "changed the first chain combobox " << label << std::endl;
}

void ramachandran_plot_differences_chain_combobox_second_changed(GtkWidget *combobox, gpointer data) {

   std::string label = get_active_label_in_combobox(GTK_COMBO_BOX(combobox));
   graphics_info_t::ramachandran_plot_differences_imol2_chain = label;
   std::cout << "changed the second chain combobox " << label << std::endl;

}


void do_ramachandran_plot(int imol) {

#ifdef HAVE_GOOCANVAS

   if (is_valid_model_molecule(imol)) {

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
                 is_kleywegt_plot_flag,
                 graphics_info_t::rama_psi_axis_mode);
      if (rama->dynawin) {
         rama->draw_it(graphics_info_t::molecules[imol].atom_sel.mol);
      } else {
         std::cout<<"WARNING:: could not initialise ramachandran\n"<<std::endl;
      }
   }
#endif // HAVE_GOOCANVAS
}

void set_kleywegt_plot_n_diffs(int ndiffs) {
   graphics_info_t::rama_n_diffs = ndiffs;
}


void
ramachandran_plot_differences(int imol1, int imol2) {

   std::cout << "ERROR:: called unconverted ramachandran_plot_differences()" << std::endl;


#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *rama_widget = dynarama_is_displayed_state(imol1); // the underlying variable is
                                                                // set by the init_internal() function
                                                                // of rama_plot.
   if (rama_widget) {
      GtkWidget *w = coot::get_validation_graph(imol1, coot::RAMACHANDRAN_PLOT);
      if (w && (imol1 == imol2)) {
         coot::rama_plot *plot =
               (coot::rama_plot *) gtk_object_get_user_data(GTK_OBJECT(w));
         gtk_widget_set_visible(w, TRUE);
         plot->make_kleywegt_plot(1);
      } else {
         // dont have the widget, make a new one?! info for now
         // FIXME
         std::string s = "Sorry. Cant find an existing plot\n";
         s += "not making a new one, yet.";
         info_dialog(s.c_str());
      }
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
               is_kleywegt_plot_flag,
               graphics_info_t::rama_psi_axis_mode);
         if (rama->dynawin) {
            rama->draw_it(imol1, imol2,
                     graphics_info_t::molecules[imol1].atom_sel.mol,
                     graphics_info_t::molecules[imol2].atom_sel.mol);
           } else {
              std::cout<<"WARNING:: could not initialise ramachandran\n"<<std::endl;
         }
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
      GtkWidget *w = coot::get_validation_graph(imol1, coot::RAMACHANDRAN_PLOT);
      if (w) {
         coot::rama_plot *plot =
               (coot::rama_plot *) gtk_object_get_user_data(GTK_OBJECT(w));
         gtk_widget_set_visible(w, TRUE);
         plot->make_kleywegt_plot(1);
      } else {
         // dont have the widget, make a new one?! info for now
         // FIXME
         std::string s = "Sorry. Cant find an existing plot\n";
         s += "not making a new one, yet.";
         info_dialog(s.c_str());
      }
   } else {
      if (is_valid_model_molecule(imol1)) {
	 if (is_valid_model_molecule(imol2)) {
       short int is_kleywegt_plot_flag = 1;
	    coot::rama_plot *rama = new coot::rama_plot;
	    rama->set_n_diffs(graphics_info_t::rama_n_diffs);
	    rama->init(imol1,
		       graphics_info_t::molecules[imol1].dotted_chopped_name(),
		       graphics_info_t::rama_level_prefered,
		       graphics_info_t::rama_level_allowed,
		       graphics_info_t::rama_plot_background_block_size,
             is_kleywegt_plot_flag,
             graphics_info_t::rama_psi_axis_mode);
	    if (graphics_info_t::molecules[imol1].is_from_shelx_ins() ||
		graphics_info_t::molecules[imol2].is_from_shelx_ins())
	       rama->allow_seqnum_offset();
// 	    std::cout << "rama differences on mols: " << imol1 << " " << a_chain
// 		      << " to " << imol2 << " " << b_chain << std::endl;
       if (rama->dynawin) {
	    rama->draw_it(imol1, imol2,
			  graphics_info_t::molecules[imol1].atom_sel.mol,
			  graphics_info_t::molecules[imol2].atom_sel.mol,
			  std::string(a_chain), std::string(b_chain));
	    rama->set_kleywegt_plot_uses_chain_ids();
       } else {
          std::cout<<"WARNING:: could not initialise ramachandran\n"<<std::endl;
       }
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
void fill_ramachandran_plot_differences_combobox_with_chain_options(GtkWidget *chain_combobox,
								    int is_first_mol_flag) {

   GtkWidget *mol_combobox = 0;

   if (is_first_mol_flag) {
      mol_combobox = widget_from_builder("ramachandran_plot_differences_first_mol_combobox");
   } else {
      mol_combobox = widget_from_builder("ramachandran_plot_differences_second_mol_combobox");
   }

   GCallback callback_func;
   int imol;

   if (is_first_mol_flag) {
      imol = graphics_info_t::ramachandran_plot_differences_imol1;
      callback_func = G_CALLBACK(ramachandran_plot_differences_chain_combobox_first_changed);
   } else {
      imol = graphics_info_t::ramachandran_plot_differences_imol2;
      callback_func = G_CALLBACK(ramachandran_plot_differences_chain_combobox_second_changed);
   }

   if (imol >=0 && imol< graphics_info_t::n_molecules()) {
      std::string set_chain =
	 graphics_info_t::fill_combobox_with_chain_options(chain_combobox, imol, callback_func);
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



#include "rama_plot_with_canvas.hh"
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

#ifdef HAVE_GOOCANVAS
   if (is_valid_model_molecule(imol)) {
      // w = graphics_info_t::dynarama_is_displayed[imol];
      w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
   }
#endif

   return w;
}


GtkWidget *dynarama_widget(int imol) {

   GtkWidget *w = NULL;
#ifdef DO_GEOMETRY_GRAPHS
   if (imol < graphics_info_t::n_molecules()) {
      // w = graphics_info_t::dynarama_is_displayed[imol];
#ifdef HAVE_GOOCANVAS
      w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
#endif
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

   std::cout << "ERROR:: using unconverted get_mol_from_dynarama() FIXME " << std::endl;

#ifdef HAVE_GOOCANVAS
   // graphics_info_t g;
   if (window) {

      GtkWidget *canvas = widget_from_builder("canvas");

      if (canvas) {
	 coot::rama_plot * plot = static_cast<coot::rama_plot *>(gtk_object_get_user_data(GTK_OBJECT(canvas)));
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
#endif // HAVE_GOOCANVAS
   return imol;
}

// resize
#if 0
// 20220602-PE this is no longer the way to handle resize events.
void
resize_rama_canvas(GtkWidget *widget, GdkEventConfigure *event) {

#ifdef HAVE_GOOCANVAS
    coot::rama_plot rp;
    rp.resize_rama_canvas_internal(widget, event);
#endif

}
#endif


void toggle_dynarama_outliers(GtkWidget *window, int state) {

#ifdef HAVE_GOOCANVAS
   int imol = get_mol_from_dynarama(window);
   // GtkWidget *canvas = lookup_widget(GTK_WIDGET(window), "canvas");
   GtkWidget *canvas = widget_from_builder("canvas");
   if (canvas) {
      // 20211201-PE FIXME - where was (should have been) this set?
      coot::rama_plot *plot = static_cast<coot::rama_plot *>(g_object_get_data(G_OBJECT(canvas), "rama_plot"));
      if (plot) {
	 if (is_valid_model_molecule(imol)) {
	    mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	    if (mol) {
	       plot->show_outliers_only(mol, state);
	    }
	 }
      }
   }
#endif
}

void set_ramachandran_psi_axis_mode(int mode) {
}

int
ramachandran_psi_axis_mode() {
   return 0;
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
   SCM list_of_alignments_as_text = SCM_EOL;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::vector<coot::chain_mutation_info_container_t> > ar =
      	 graphics_info_t::molecules[imol].residue_mismatches(graphics_info_t::alignment_wgap,
							     graphics_info_t::alignment_wspace);
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

      if ((mutations.size() > 0) || (insertions.size() > 0) || (deletions.size() > 0)) {
	 SCM insertions_scm = SCM_EOL;
	 SCM deletions_scm = SCM_EOL;
	 SCM mutations_scm = SCM_EOL;
	 for (unsigned int i=0; i<mutations.size(); i++) {
	    SCM rs_scm = residue_spec_to_scm(mutations[i].first);
	    SCM str = scm_from_locale_string(mutations[i].second.c_str());
	    SCM c = SCM_EOL;
	    c = scm_cons(str, c);
	    c = scm_cons(str, rs_scm);
	    mutations_scm = scm_cons(c, mutations_scm);
	 }
	 for (unsigned int i=0; i<insertions.size(); i++) {
	    SCM rs_scm = residue_spec_to_scm(insertions[i].first);
	    SCM str = scm_from_locale_string(insertions[i].second.c_str());
	    SCM c = SCM_EOL;
	    c = scm_cons(str, c);
	    c = scm_cons(str, rs_scm);
	    insertions_scm = scm_cons(c, insertions_scm);
	 }
	 for (unsigned int i=0; i<deletions.size(); i++) {
	    SCM rs_scm = residue_spec_to_scm(deletions[i].first);
	    SCM str = scm_from_locale_string(deletions[i].second.c_str());
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

	 for (std::size_t i=0; i<ar.second.size(); i++) {

	    // and the dialog text

	    const coot::chain_mutation_info_container_t &mic = ar.second[i];

	    std::cout << ":::::::::::::: Here with alignment_string: " << mic.alignment_string
		      << std::endl;

	    SCM alignment_as_text_scm = scm_from_locale_string(mic.alignment_string.c_str());

	    list_of_alignments_as_text = scm_cons(alignment_as_text_scm, list_of_alignments_as_text);
	 }


	 // Put list_of_alignments_as_text at the end of r
	 r = scm_reverse(scm_cons(list_of_alignments_as_text, scm_reverse(r)));
      }
   }
   return r;
}
#endif // USE_GUILE


/*! \brief do a internal alignment of all the assigned sequences,
  return a list of mismatches that need to be made to model number
  imol to match the input sequence.

Return a list of mutations deletions insetions.
Return  False on failure to align (e.g. not assigned sequence)
and the empty list on no alignment mismatches.

Returns a list of alignment as text as 4th element
*/
#ifdef USE_PYTHON
PyObject *alignment_mismatches_py(int imol) {

   PyObject *r = Py_False;

   std::vector<std::pair<coot::residue_spec_t,std::string> > mutations;
   std::vector<std::pair<coot::residue_spec_t,std::string> > insertions;
   std::vector<std::pair<coot::residue_spec_t,std::string> > deletions;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::vector<coot::chain_mutation_info_container_t> > ar =
      	 graphics_info_t::molecules[imol].residue_mismatches(graphics_info_t::alignment_wgap,
							     graphics_info_t::alignment_wspace);
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
	 PyObject *rs_py = residue_spec_to_py(mutations[i].first);
	 PyObject *str = myPyString_FromString(mutations[i].second.c_str());
	 PyList_Insert(rs_py, 0, str);
	 PyList_Append(mutations_py, rs_py);
	 Py_XDECREF(str);
	 Py_XDECREF(rs_py);
      }
      for (unsigned int i=0; i<insertions.size(); i++) {
	 PyObject *rs_py = residue_spec_to_py(insertions[i].first);
	 PyObject *str = myPyString_FromString(insertions[i].second.c_str());
	 PyList_Insert(rs_py, 0, str);
	 PyList_Append(insertions_py, rs_py);
	 Py_XDECREF(str);
	 Py_XDECREF(rs_py);
      }
      for (unsigned int i=0; i<deletions.size(); i++) {
	 PyObject *rs_py = residue_spec_to_py(deletions[i].first);
	 PyObject *str = myPyString_FromString(deletions[i].second.c_str());
	 PyList_Insert(rs_py, 0, str);
	 PyList_Append(deletions_py, rs_py);
	 Py_XDECREF(str);
	 Py_XDECREF(rs_py);
      }

      // 20220228-PE merge resolution - I think I want this block but it uses old style python API
      // fix it later.
#if 0
      r = PyList_New(4);
      // These are reversed so that the residue numbers come out in
      // numerical order (not backwards) and the returned list is
      // [mutations, deletions, insertions].
      PyList_SetItem(r, 0, mutations_py);
      PyList_SetItem(r, 1, deletions_py);
      PyList_SetItem(r, 2, insertions_py);
      
      PyObject *list_of_alignments_as_text = PyList_New(0);
      for (std::size_t i=0; i<ar.second.size(); i++) {
         
         // and the dialog text
         
         const coot::chain_mutation_info_container_t &mic = ar.second[i];
         
         PyObject *alignment_as_text_py = PyString_FromString(mic.alignment_string.c_str());
         PyList_Append(list_of_alignments_as_text, alignment_as_text_py);
      }
      
      // Put list_of_alignments_as_text at the end of r
      PyList_SetItem(r, 3, list_of_alignments_as_text);
#endif

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
	    for (unsigned int i=0; i<v.size(); i++) {
	       std::cout << v[i].second << std::endl;
	    }
#if 0 // 20211003-PE Not today, guile gui
	    graphics_info_t g;
	    std::vector<coot::util::atom_spec_and_button_info_t> outlier_atoms;
	    for (unsigned int i=0; i<v.size(); i++) {
	       std::string callback_func = "(lambda() (do-180-degree-side-chain-flip ";
	       callback_func += coot::util::int_to_string(imol);
	       callback_func += " ";
	       callback_func += single_quote(v[i].first.chain_id);
	       callback_func += " ";
	       callback_func += coot::util::int_to_string(v[i].first.res_no);
	       callback_func += " ";
	       callback_func += single_quote(v[i].first.ins_code);
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
#endif

            graphics_info_t g;
            std::vector<coot::util::atom_spec_and_button_info_t> outlier_atoms;
            for (unsigned int i=0; i<v.size(); i++) {
               std::string callback_func = "[coot.do_180_degree_side_chain_flip,";
               callback_func += coot::util::int_to_string(imol);
               callback_func += ",";
               callback_func += single_quote(v[i].first.chain_id);
               callback_func += ",";
               callback_func += coot::util::int_to_string(v[i].first.res_no);
               callback_func += ",";
               callback_func += single_quote(v[i].first.ins_code);
               callback_func += ",";
               callback_func += single_quote(v[i].first.alt_conf);
               callback_func += "]";
               v[i].first.int_user_data = imol; // kludge in the imol, used for callback.
               coot::util::atom_spec_and_button_info_t asi(v[i].first, v[i].second, callback_func);
               outlier_atoms.push_back(asi);
            }
            std::string error_type = "Z score: ";
            std::vector<std::string> cmd_strings;
            // "import coot_gui" should happen when coot starts with gui.
            cmd_strings.push_back("import coot_gui ; coot_gui.interesting_things_with_fix_maybe");
            cmd_strings.push_back(single_quote("GLN and ASN B-factor Outliers"));
            std::string ls = coot::util::interesting_things_list_with_fix_py(outlier_atoms, error_type);
            cmd_strings.push_back(ls);
            std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
            std::cout << "python command: " << s << std::endl;
            safe_python_command(s);

	 } else {
	    std::string label = "Coot detected no GLN or ASN B-factor Outliers";
	    GtkWidget *w = wrapped_nothing_bad_dialog(label);
	    gtk_widget_set_visible(w, TRUE);
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
	    for (unsigned int i=0; i<v.size(); i++) {
	       std::cout << v[i].second << std::endl;
	    }
#ifdef USE_PYGTK
            graphics_info_t g;
            std::vector<coot::util::atom_spec_and_button_info_t> outlier_atoms;
            for (unsigned int i=0; i<v.size(); i++) {
               std::string callback_func = "[do_180_degree_side_chain_flip,";
               callback_func += coot::util::int_to_string(imol);
               callback_func += ",";
               callback_func += single_quote(v[i].first.chain_id);
               callback_func += ",";
               callback_func += coot::util::int_to_string(v[i].first.res_no);
               callback_func += ",";
               callback_func += single_quote(v[i].first.ins_code);
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
	    gtk_widget_set_visible(w, TRUE);
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
	 r = scm_from_double(tors);
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
      r = scm_from_double(new_tors);
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
   if (PyBool_Check(r)) {
      Py_XINCREF(r);
   }
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


#ifdef USE_GUILE
SCM all_molecule_rotamer_score(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::rotamer_score_t rs = g.all_molecule_rotamer_score(imol);
      SCM a_scm = scm_from_double(rs.score);
      SCM b_scm = scm_from_int(rs.n_rotamer_residues());
      r = scm_list_2(a_scm, b_scm);
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *all_molecule_rotamer_score_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::rotamer_score_t rs = g.all_molecule_rotamer_score(imol);
      PyObject *a_py = PyFloat_FromDouble(rs.score);
      PyObject *b_py = PyLong_FromLong(rs.n_rotamer_residues());
      r = PyList_New(2);
      PyList_SetItem(r, 0, a_py);
      PyList_SetItem(r, 1, b_py);
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}
#endif

#ifdef USE_GUILE
SCM all_molecule_ramachandran_score(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      coot::rama_score_t rs = graphics_info_t::molecules[imol].get_all_molecule_rama_score();
      SCM a_scm = scm_from_double(rs.score);
      SCM b_scm = scm_from_int(rs.n_residues());
      SCM c_scm = scm_from_double(rs.score_non_sec_str);
      SCM d_scm = scm_from_int(rs.n_residues_non_sec_str());
      SCM e_scm = scm_from_int(rs.n_zeros);
      SCM by_residue_scm = SCM_EOL;
      for (std::size_t ii=0; ii<rs.scores.size(); ii++) {
	 SCM residue_spec_scm = residue_spec_to_scm(rs.scores[ii].res_spec);
	 SCM d_scm = scm_from_double(rs.scores[ii].score);
	 SCM phi_scm = scm_from_double(rs.scores[ii].phi_psi.phi());
	 SCM psi_scm = scm_from_double(rs.scores[ii].phi_psi.psi());
	 SCM phi_psi_scm = scm_list_2(phi_scm, psi_scm);
	 if (false)
	    std::cout << "here with residue pointers "
		      << rs.scores[ii].residue_prev << " "
		      << rs.scores[ii].residue_this << " "
		      << rs.scores[ii].residue_next << " "
		      << std::endl;
	 if (rs.scores[ii].residue_prev &&
	     rs.scores[ii].residue_this &&
	     rs.scores[ii].residue_next) {
	    SCM res_names_scm = scm_list_3(scm_from_locale_string(rs.scores[ii].residue_prev->GetResName()),
					  scm_from_locale_string(rs.scores[ii].residue_this->GetResName()),
					  scm_from_locale_string(rs.scores[ii].residue_next->GetResName()));
	    SCM residue_results_scm = scm_list_4(phi_psi_scm, residue_spec_scm, d_scm, res_names_scm);
	    by_residue_scm = scm_cons(residue_results_scm, by_residue_scm);
	 } else {
	    SCM residue_results_scm = scm_list_3(phi_psi_scm, residue_spec_scm, d_scm);
	    by_residue_scm = scm_cons(residue_results_scm, by_residue_scm);
	 }
      }
      // r = SCM_LIST6(a_scm, b_scm, c_scm, d_scm, e_scm, scm_reverse(by_residue_scm));
      // new guile
      r = scm_cons2(a_scm, b_scm, scm_list_4(c_scm, d_scm, e_scm, scm_reverse(by_residue_scm)));
   }

   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *all_molecule_ramachandran_score_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      coot::rama_score_t rs = graphics_info_t::molecules[imol].get_all_molecule_rama_score();
      PyObject *a_py = PyFloat_FromDouble(rs.score);
      PyObject *b_py = PyLong_FromLong(rs.n_residues());
      PyObject *c_py = PyFloat_FromDouble(rs.score_non_sec_str);
      PyObject *d_py = PyLong_FromLong(rs.n_residues_non_sec_str());
      PyObject *e_py = PyLong_FromLong(rs.n_zeros);
      PyObject *info_by_residue_py = PyList_New(rs.scores.size());
      for (std::size_t ii=0; ii<rs.scores.size(); ii++) {
	 PyObject *info_for_residue_py = PyList_New(4);
	 PyObject *residue_spec_py = residue_spec_to_py(rs.scores[ii].res_spec);
	 if (rs.scores[ii].residue_prev &&
	     rs.scores[ii].residue_this &&
	     rs.scores[ii].residue_next) {
	    PyObject *phi_py = PyFloat_FromDouble(rs.scores[ii].phi_psi.phi());
	    PyObject *psi_py = PyFloat_FromDouble(rs.scores[ii].phi_psi.psi());
            PyObject *residue_score_py = PyFloat_FromDouble(rs.scores[ii].score);
            PyObject *phi_psi_py = PyList_New(2);
	    PyObject *res_names_py = PyList_New(3);
	    PyList_SetItem(phi_psi_py, 0, phi_py);
	    PyList_SetItem(phi_psi_py, 1, psi_py);
	    PyList_SetItem(res_names_py, 0, myPyString_FromString(rs.scores[ii].residue_prev->GetResName()));
	    PyList_SetItem(res_names_py, 1, myPyString_FromString(rs.scores[ii].residue_this->GetResName()));
	    PyList_SetItem(res_names_py, 2, myPyString_FromString(rs.scores[ii].residue_next->GetResName()));
	    PyList_SetItem(info_for_residue_py, 0, phi_psi_py);
	    PyList_SetItem(info_for_residue_py, 1, residue_spec_py);
            PyList_SetItem(info_for_residue_py, 2, residue_score_py);
	    PyList_SetItem(info_for_residue_py, 3, res_names_py);
	    PyList_SetItem(info_by_residue_py, ii, info_for_residue_py);
	 } else {
	    PyList_SetItem(info_by_residue_py, ii, PyLong_FromLong(-1)); // shouldn't happen
	 }
      }
      r = PyList_New(6);
      PyList_SetItem(r, 0, a_py);
      PyList_SetItem(r, 1, b_py);
      PyList_SetItem(r, 2, c_py);
      PyList_SetItem(r, 3, d_py);
      PyList_SetItem(r, 4, e_py);
      PyList_SetItem(r, 5, info_by_residue_py);
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }

   return r;
}

PyObject *all_molecule_ramachandran_region_py(int imol) {

   PyObject *r = Py_False;
   PyObject *pair;
   if (is_valid_model_molecule(imol)) {
      coot::rama_score_t rs = graphics_info_t::molecules[imol].get_all_molecule_rama_score();
      std::vector<std::pair<coot::residue_spec_t, int> > rama_region = rs.region;
      int region_size = rama_region.size();
      if (region_size > 0) {
        r = PyList_New(region_size);
        for (unsigned int i=0; i<rama_region.size(); i++) {
          pair = PyTuple_New(2);
          PyTuple_SetItem(pair, 0, residue_spec_to_py(rama_region[i].first));
          PyTuple_SetItem(pair, 1, PyLong_FromLong(rama_region[i].second));
          PyList_SetItem(r, i, pair);
        }
      } else {
        std::cout << "INFO:: empty ramachandran region list" << std::endl;
      }
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}

#include "coot-utils/atom-overlaps.hh"

//! \brief return 1 if this residue clashes with the symmetry-related
//!  atoms of the same molecule.
//!
//! 0 means that it did not clash,
//! -1 means that the residue or molecule could not be found or that there
//!    was no cell and symmetry.
int clashes_with_symmetry(int imol, const char *chain_id, int res_no, const char *ins_code,
			  float clash_dist) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::residue_spec_t spec(chain_id, res_no, ins_code);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(spec);
      if (residue_p) {
	 if (mol) {
	    std::vector<mmdb::Residue *> dummy; // neighbours
	    coot::atom_overlaps_container_t ao(residue_p, dummy, mol, g.Geom_p());
	    std::vector<coot::atom_overlap_t> v = ao.symmetry_contacts(clash_dist);
	    if (v.empty())
	       r = 0;
	    else
	       r = 1;
	 }
      }

   }
   return r;
}

#include "analysis/b-factor-histogram.hh"

#ifdef HAVE_GOOCANVAS
#include "goograph/goograph.hh"
#endif

//! B-factor distribution histogram
void b_factor_distribution_graph(int imol) {

#ifdef HAVE_GOOCANVAS
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      coot::b_factor_histogram b(mol);
      b.model();
      std::vector<std::pair<double, double> > data  = b.get_data();
      std::vector<std::pair<double, double> > model = b.get_model();

      coot::goograph* g = new coot::goograph();
      int trace_data  = g->trace_new();
      int trace_model = g->trace_new();

      g->set_trace_type(trace_data, coot::graph_trace_info_t::PLOT_TYPE_BAR);
      g->set_trace_colour(trace_data, "#88bb88");
      g->set_data(trace_data, data);

      // Don't draw the model until it is correctly calculated
      // g->set_trace_type(trace_model, coot::graph_trace_info_t::PLOT_TYPE_LINE);
      // g->set_trace_colour(trace_model, "#333333");
      // g->set_data(trace_model, model);

      g->set_plot_title("B-factor histogram");
      // g->set_extents(coot::goograph::X_AXIS, 0,  140);
      // g->set_extents(coot::goograph::Y_AXIS, 0, 2000);

      g->draw_graph();
      g->show_dialog();

   }
#endif // HAVE_GOOCANVASHAVE_GOOCANVAS
}



#endif // USE_PYTHON

// SWIG testing.
//
// Foo foo_test() {
//    Foo f;
//    return f;
// }

// std::vector<Foo>
// foo_vec() {
//    std::vector<Foo> v;
//    v.push_back(Foo());
//    v.push_back(Foo());
//    v.push_back(Foo());
//    return v;
// }

#include "python-classes.hh"

/*  ----------------------------------------------------------------------- */
/*                  Pathology Plots                                         */
/*  ----------------------------------------------------------------------- */
#ifdef USE_PYTHON
PyObject *pathology_data(const std::string &mtz_file_name,
			 const std::string &fp_col,
			 const std::string &sigfp_col) {

   PyObject *r = Py_False;

   std::vector<std::pair<double, double> >   fp_vs_reso_data;
   std::vector<std::pair<double, double> > fosf_vs_reso_data;
   std::vector<std::pair<double, double> >      sf_vs_f_data;
   std::vector<std::pair<double, double> >    fosf_vs_f_data;
   // how about sigF vs resolution also?
   double invresolsq_max = 0;

   try {
      clipper::CCP4MTZfile mtz;
      std::cout << "INFO:: reading mtz file " << mtz_file_name << std::endl;
      mtz.open_read(mtz_file_name);
      clipper::HKL_data< clipper::datatypes::F_sigF<float> > fsigf;
      std::string dataname = "/*/*/[" + fp_col + " " + sigfp_col + "]";
      mtz.import_hkl_data(fsigf, dataname);
      mtz.close_read();

      int n_reflns = fsigf.num_obs();
      clipper::HKL_info::HKL_reference_index hri;
      for (hri=fsigf.first(); !hri.last(); hri.next()) {
         if (! clipper::Util::isnan(fsigf[hri].f())) {
            double invresolsq = hri.invresolsq();
            std::pair<double, double> p(invresolsq, fsigf[hri].f());
            fp_vs_reso_data.push_back(p);
            const double &f    = fsigf[hri].f();
            const double &sigf = fsigf[hri].sigf();
            if (! clipper::Util::isnan(sigf)) {
               if (sigf != 0) {
                  std::pair<double, double> p1(invresolsq, f/sigf);
                  fosf_vs_reso_data.push_back(p1);
                  std::pair<double, double> p2(f, sigf);
                  sf_vs_f_data.push_back(p2);
                  std::pair<double, double> p3(f, f/sigf);
                  fosf_vs_f_data.push_back(p3);
                  if (invresolsq > invresolsq_max)
                  invresolsq_max = invresolsq;
               }
            }
         }
      }
   }
   catch (const clipper::Message_fatal &e) {
      std::cout << "error: " << e.text() << std::endl;
   }

   std::cout << "INFO:: pathology_plots() found "
   << fp_vs_reso_data.size() << " data" << std::endl;


   // this is just a bit of fun - looking for large FP outliers.
   //
   if (false) {
      int n_data = fp_vs_reso_data.size();
      for (std::size_t i=0; i<fp_vs_reso_data.size(); i++) {
	 int i_int = i;
	 int low_lim = i-20;
	 int high_lim = i+20;
	 if (low_lim < 0) low_lim = 0;
	 if (high_lim >= n_data) high_lim = n_data;
	 double sum = 0;
	 double n = 0;
	 for (int j=low_lim; j<high_lim; j++) {
	    if (j != i_int) {
	       sum += fp_vs_reso_data[j].second;
	       n += 1;
	    }
	 }
	 double local_average =  sum/n;
	 double this_f    = fp_vs_reso_data[i].second;
	 double this_reso = fp_vs_reso_data[i].first;
	 if (this_f > 3 * local_average) {
	    std::cout << "   " << this_reso << " " << this_f << " / " << local_average
		      << " = " << this_f/local_average << std::endl;
	 }
      }
   }

   if (  fp_vs_reso_data.size() > 0 &&
       fosf_vs_reso_data.size() > 0 &&
  	    sf_vs_f_data.size() > 0 &&
          fosf_vs_f_data.size() > 0) {

      if (false) {
	 PyTypeObject *type = NULL; // should be something
	 PyObject *args = NULL;
	 PyObject *kwds = NULL;
	 // PyObject *test_object = PathologyData_new(type, args, kwds);
      }

      unsigned int data_size = fp_vs_reso_data.size();
      if (data_size > 20000) {
	 double r = double(20000)/double(data_size);
	 fp_vs_reso_data.erase(std::remove_if(fp_vs_reso_data.begin(),
					      fp_vs_reso_data.end(),
					      data_pair_remover(r)),
			       fp_vs_reso_data.end());
	 fosf_vs_reso_data.erase(std::remove_if(fosf_vs_reso_data.begin(),
					      fosf_vs_reso_data.end(),
					      data_pair_remover(r)),
			       fosf_vs_reso_data.end());
	 sf_vs_f_data.erase(std::remove_if(sf_vs_f_data.begin(),
					      sf_vs_f_data.end(),
					      data_pair_remover(r)),
			       sf_vs_f_data.end());
	 fosf_vs_f_data.erase(std::remove_if(fosf_vs_f_data.begin(),
					      fosf_vs_f_data.end(),
					      data_pair_remover(r)),
			       fosf_vs_f_data.end());

	 if (0) {
	    std::cout << "  now data size " << fp_vs_reso_data.size() << "" << std::endl;
	    std::cout << "  now data size " << fosf_vs_reso_data.size() << "" << std::endl;
	    std::cout << "  now data size " << sf_vs_f_data.size() << "" << std::endl;
	    std::cout << "  now data size " << fosf_vs_f_data.size() << "" << std::endl;
	 }

      }

      PyObject *r0 = PyList_New(fp_vs_reso_data.size());
      for (unsigned int i=0; i<fp_vs_reso_data.size(); i++) {
	 PyObject *o = PyTuple_New(2);
	 PyTuple_SetItem(o, 0, PyFloat_FromDouble(fp_vs_reso_data[i].first));
	 PyTuple_SetItem(o, 1, PyFloat_FromDouble(fp_vs_reso_data[i].second));
	 PyList_SetItem(r0, i, o);
      }
      PyObject *r1 = PyList_New(fosf_vs_reso_data.size());
      for (unsigned int i=0; i<fosf_vs_reso_data.size(); i++) {
	 PyObject *o = PyTuple_New(2);
	 PyTuple_SetItem(o, 0, PyFloat_FromDouble(fosf_vs_reso_data[i].first));
	 PyTuple_SetItem(o, 1, PyFloat_FromDouble(fosf_vs_reso_data[i].second));
	 PyList_SetItem(r1, i, o);
      }
      PyObject *r2 = PyList_New(sf_vs_f_data.size());
      for (unsigned int i=0; i<sf_vs_f_data.size(); i++) {
	 PyObject *o = PyTuple_New(2);
	 PyTuple_SetItem(o, 0, PyFloat_FromDouble(sf_vs_f_data[i].first));
	 PyTuple_SetItem(o, 1, PyFloat_FromDouble(sf_vs_f_data[i].second));
	 PyList_SetItem(r2, i, o);
      }
      PyObject *r3 = PyList_New(fosf_vs_f_data.size());
      for (unsigned int i=0; i<fosf_vs_f_data.size(); i++) {
	 PyObject *o = PyTuple_New(2);
	 PyTuple_SetItem(o, 0, PyFloat_FromDouble(fosf_vs_f_data[i].first));
	 PyTuple_SetItem(o, 1, PyFloat_FromDouble(fosf_vs_f_data[i].second));
	 PyList_SetItem(r3, i, o);
      }
      r = PyList_New(5);
      PyList_SetItem(r, 0, PyFloat_FromDouble(invresolsq_max));
      PyList_SetItem(r, 1, r0);
      PyList_SetItem(r, 2, r1);
      PyList_SetItem(r, 3, r2);
      PyList_SetItem(r, 4, r3);


   }

   if (PyBool_Check(r))
     Py_INCREF(r);

   return r;
}
#endif // USE_PYTHON

#include "coot-utils/cablam-markup.hh"
#include "c-interface-generic-objects.h"

std::vector<std::pair<coot::residue_spec_t, double> >
add_cablam_markup(int imol, const std::string &cablam_log_file_name) {

   std::vector<std::pair<coot::residue_spec_t, double> > residues_vec;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<coot::cablam_markup_t> v = coot::make_cablam_markups(mol, cablam_log_file_name);

      std::cout << "INFO:: Made " << v.size() << " cablam markups " << std::endl;
      std::vector<coot::cablam_markup_t>::const_iterator it;
      int idx_cablam = generic_object_index("xxCaBLAM");
      if (idx_cablam == -1)
         idx_cablam = new_generic_object_number("xxCaBLAM");
      else
         generic_object_clear(idx_cablam);
      set_display_generic_object(idx_cablam, 0); // don't display while we are adding
      if (v.size() > 0) {
         for (it=v.begin(); it!=v.end(); ++it) {
            std::pair<coot::residue_spec_t, double> p(coot::residue_spec_t(it->residue), it->score);
            residues_vec.push_back(p);
         }
      }
      for (it=v.begin(); it!=v.end(); ++it) {
         const coot::cablam_markup_t &cm(*it);
         to_generic_object_add_point(idx_cablam, "pink", 14, cm.O_prev_pos.x(), cm.O_prev_pos.y(), cm.O_prev_pos.z());
         to_generic_object_add_point(idx_cablam, "pink", 14, cm.O_this_pos.x(), cm.O_this_pos.y(), cm.O_this_pos.z());
         to_generic_object_add_point(idx_cablam, "pink", 14, cm.O_next_pos.x(), cm.O_next_pos.y(), cm.O_next_pos.z());

         if (false) {
            std::cout << "line 1: " << cm.O_this_pos.format() << " to " << cm.CA_proj_point_this.format() << std::endl;
            std::cout << "line 2: " << cm.O_prev_pos.format() << " to " << cm.CA_proj_point_prev.format() << std::endl;
            std::cout << "line 3: " << cm.O_next_pos.format() << " to " << cm.CA_proj_point_next.format() << std::endl;
         }

         to_generic_object_add_line(idx_cablam, "hotpink", 2,
                                    cm.O_this_pos.x(), cm.O_this_pos.y(), cm.O_this_pos.z(),
                                    cm.CA_proj_point_this.x(), cm.CA_proj_point_this.y(), cm.CA_proj_point_this.z());

         to_generic_object_add_line(idx_cablam, "hotpink", 2,
                                    cm.O_prev_pos.x(), cm.O_prev_pos.y(), cm.O_prev_pos.z(),
                                    cm.CA_proj_point_prev.x(), cm.CA_proj_point_prev.y(), cm.CA_proj_point_prev.z());

         to_generic_object_add_line(idx_cablam, "hotpink", 2,
                                    cm.O_next_pos.x(), cm.O_next_pos.y(), cm.O_next_pos.z(),
                                    cm.CA_proj_point_next.x(), cm.CA_proj_point_next.y(), cm.CA_proj_point_next.z());

         to_generic_object_add_line(idx_cablam, "hotpink", 2,
                                    cm.CA_proj_point_this.x(), cm.CA_proj_point_this.y(), cm.CA_proj_point_this.z(),
                                    cm.CA_proj_point_prev.x(), cm.CA_proj_point_prev.y(), cm.CA_proj_point_prev.z());

         to_generic_object_add_line(idx_cablam, "hotpink", 2,
                                    cm.CA_proj_point_this.x(), cm.CA_proj_point_this.y(), cm.CA_proj_point_this.z(),
                                    cm.CA_proj_point_next.x(), cm.CA_proj_point_next.y(), cm.CA_proj_point_next.z());

      }
      set_display_generic_object(idx_cablam, 1); // now we can see it
      graphics_draw();
   }
   return residues_vec;
}


#ifdef USE_GUILE
SCM add_cablam_markup_scm(int imol, const std::string &cablam_log_file_name) {

   std::vector<std::pair<coot::residue_spec_t, double> > v = add_cablam_markup(imol, cablam_log_file_name);
   SCM r = SCM_EOL;
   std::vector<std::pair<coot::residue_spec_t, double> >::const_iterator it;
   for (it=v.begin(); it!=v.end(); it++) {
      SCM item_scm = scm_list_2(residue_spec_to_scm(it->first), scm_from_double(it->second));
      r = scm_cons(item_scm, r);
   }
   r = scm_reverse(r);
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *add_cablam_markup_py(int imol, const std::string &cablam_log_file_name) {

   std::vector<std::pair<coot::residue_spec_t, double> > v = add_cablam_markup(imol, cablam_log_file_name);
   PyObject *r = PyList_New(v.size());
   for (unsigned int i=0; i<v.size(); i++) {
      const double score = v[i].second;
      const coot::residue_spec_t &spec(v[i].first);
      PyObject *item_py = PyList_New(2);
      PyList_SetItem(item_py, 0, residue_spec_to_py(spec));
      PyList_SetItem(item_py, 1, PyFloat_FromDouble(score));
      PyList_SetItem(r, i, item_py);
   }
   return r;

}
#endif // USE_PYTHON


float atom_overlap_score(int imol) {

   float s = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      graphics_info_t g;
      bool ignore_waters_flag = false;
      coot::atom_overlaps_container_t ao(mol, g.Geom_p(), ignore_waters_flag);
      ao.make_all_atom_overlaps();
      s = ao.score();
   }


   return s;
}
