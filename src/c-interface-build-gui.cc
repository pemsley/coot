/* src/c-interface-build-gui.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007,2008, 2009 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009 The University of Oxford
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
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
#include <windows.h>
#endif

#include "globjects.h" //includes gtk/gtk.h

#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include "graphics-info.h"

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-fasta.hh"

#include "skeleton/BuildCas.h"
#include "ligand/helix-placement.hh"
#include "ligand/fast-ss-search.hh"

#include "utils/coot-utils.hh"  // for is_member_p
#include "coot-utils/coot-map-heavy.hh"  // for fffear

#include "widget-headers.hh"
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
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"

#include "widget-from-builder.hh"

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

#if 0 // 20220601-PE don't use this function - kept for linking failure diagnostics
void do_regularize_kill_delete_dialog() {
}
#endif


/* moving gtk functionn out of build functions, delete_atom() updates
   the go to atom atom list on deleting an atom  */
void update_go_to_atom_residue_list(int imol) {

   // called from delete_atom and nowhere else.

#if 0 // lookup_widget cleansing

   std::cout << "::::::::::::::::::::::::::::::::: update_go_to_atom_residue_list() " << std::endl;

   graphics_info_t g;
   if (g.go_to_atom_window) {
      int go_to_atom_imol = g.go_to_atom_molecule();
      if (go_to_atom_imol == imol) {

         // The go to atom molecule matched this molecule, so we
         // need to regenerate the residue and atom lists.
         // GtkWidget *gtktree = lookup_widget(g.go_to_atom_window, "go_to_atom_residue_tree");
         // GtkWidget *gtk_atom_list = lookup_widget(g.go_to_atom_window, "go_to_atom_atom_list");
         g.fill_go_to_atom_residue_tree_and_atom_list_gtk2(imol, gtktree, gtk_atom_list);
      }
   }
#endif
}


void
post_delete_item_dialog() {

   // GtkWidget *w = wrapped_create_delete_item_dialog();
   // gtk_widget_set_visible(w, TRUE);

}




/* We need to set the pending delete flag and that can't be done in
   callback, so this wrapper does it */
GtkWidget *wrapped_create_delete_item_dialog() {

   // GtkWidget *widget = create_delete_item_dialog();
   // GtkWidget *atom_toggle_button;

   // if (delete_item_mode_is_atom_p()) {
   //    atom_toggle_button = lookup_widget(GTK_WIDGET(widget), "delete_item_atom_radiobutton");
   //    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(atom_toggle_button), TRUE);
   //    std::cout << "Click on the atom that you wish to delete\n";
   // } else {
   //    if (delete_item_mode_is_water_p()) {
   //       GtkWidget *water_toggle_button = lookup_widget(widget,
   //      						"delete_item_water_radiobutton");
   //       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(water_toggle_button), TRUE);
   //    } else {
   //       if (delete_item_mode_is_sidechain_p()) {
   //       GtkWidget *sidechain_toggle_button = lookup_widget(widget,
   //      						"delete_item_sidechain_radiobutton");
   //       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sidechain_toggle_button), TRUE);

   //       set_delete_sidechain_mode();
   //       std::cout << "Click on an atom in the residue that you wish to delete\n";
   //       } else {
   //          if (delete_item_mode_is_chain_p()) {
   //             GtkWidget *chain_toggle_button = lookup_widget(widget,
   //      							  "delete_item_chain_radiobutton");
   //             set_delete_chain_mode();
   //             gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chain_toggle_button), TRUE);
   //             std::cout << "Click on an atom in the chain that you wish to delete\n";
   //          } else {

   //             if (delete_item_mode_is_sidechain_range_p()) {

   //      	  GtkWidget *sidechain_range_toggle_button = lookup_widget(widget,
   //      								   "delete_item_sidechain_range_radiobutton");
   //      	  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sidechain_range_toggle_button), TRUE);

   //      	  set_delete_sidechain_range_mode();
   //             } else {

   //      	  // if (delete_item_mode_is_residue_p()) {
   //      	  // if nothing else, let's choose delete residue mode
   //      	  if (true) {
   //      	     GtkWidget *chain_toggle_button = lookup_widget(widget,
   //      							    "delete_item_residue_radiobutton");
   //      	     set_delete_residue_mode();
   //      	     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chain_toggle_button), TRUE);
   //      	  }
   //             }
   //          }
   //       }
   //    }
   // }
   // graphics_info_t::pick_pending_flag = 1;
   // pick_cursor_maybe();
   // set_transient_and_position(COOT_DELETE_WINDOW, widget);
   // store_delete_item_widget(widget);
   // return widget;

   return 0;
}

// -----------------------------------------------------
//  move molecule here widget
// -----------------------------------------------------
// GtkWidget *wrapped_create_move_molecule_here_dialog() {

   // GtkWidget *w = create_move_molecule_here_dialog();
   //    fill_move_molecule_here_dialog(w);
   // return w;

//   return 0;
// }

// called (also) by the callback of toggling the
// move_molecule_here_big_molecules_checkbutton.
//
void
fill_move_molecule_here_frame(GtkWidget *w) {

   // // GtkWidget *option_menu  = lookup_widget(w, "move_molecule_here_optionmenu");

   graphics_info_t g;

   GtkWidget *combobox = widget_from_builder("move_molecule_here_combobox");

   // GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(graphics_info_t::move_molecule_here_item_select);

   GCallback callback_func = G_CALLBACK(graphics_info_t::move_molecule_here_combobox_changed);

   int imol_active = first_coords_imol();

   g.move_molecule_here_molecule_number = imol_active;

   gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox));
   g.fill_combobox_with_coordinates_options(combobox, callback_func, imol_active);

}


void move_molecule_here_by_widget() { // no widget needed.

   int imol = graphics_info_t::move_molecule_here_molecule_number;
   move_molecule_to_screen_centre_internal(imol);
   std::vector<std::string> command_strings;
   command_strings.push_back("move-molecule-here");
   command_strings.push_back(clipper::String(imol));
   add_to_history(command_strings);
}


void fill_place_atom_molecule_combobox(GtkWidget *combobox) {

   if (! combobox) {
      std::cout << "ERROR fill_place_atom_molecule_combobox() null combobox" << std::endl;
      return;
   }

   graphics_info_t g;
   GCallback callback_func = G_CALLBACK(g.pointer_atom_molecule_combobox_changed);
   int imol_active = g.user_pointer_atom_molecule;
   if (! is_valid_model_molecule(imol_active))
      imol_active = first_coords_imol();
   g.fill_combobox_with_coordinates_options(combobox, callback_func, imol_active);

}


/* Now the refinement weight can be set from an entry in the refine_params_dialog. */
void set_refinement_weight_from_entry(GtkWidget *entry) {

   const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      float f = coot::util::string_to_float(text);
      graphics_info_t::geometry_vs_map_weight = f;
   }
   catch (const std::runtime_error &rte) {
      std::cout << "in set_refinemenent_weight_from_entry " << rte.what() << std::endl;
   }
}

void add_estimated_map_weight_to_entry(GtkWidget *entry) {

   int imol_map = imol_refinement_map();
   if (is_valid_map_molecule(imol_map)) {
      float v = estimate_map_weight(imol_map);
      graphics_info_t::geometry_vs_map_weight = v;
      std::string t = coot::util::float_to_string(v);
      gtk_editable_set_text(GTK_EDITABLE(entry), t.c_str());
   }

}

float estimate_map_weight(int imol_map) {

   graphics_info_t g;
   float w = g.get_estimated_map_weight(imol_map);
   return w;
}


void place_atom_at_pointer_by_window() {

   // put up a widget which has a OK callback button which does a
   // g.place_typed_atom_at_pointer();
   // GtkWidget *window = create_pointer_atom_type_dialog();

   GtkWidget *window = widget_from_builder("pointer_atom_type_dialog");

   //   GtkSignalFunc callback_func =
   // 	GTK_SIGNAL_FUNC(graphics_info_t::pointer_atom_molecule_menu_item_activate);

   // GtkWidget *optionmenu = lookup_widget(window, "pointer_atom_molecule_optionmenu");
   //   fill_place_atom_molecule_option_menu(optionmenu);

   GtkWidget *combobox = widget_from_builder("pointer_atom_molecule_combobox");
   fill_place_atom_molecule_combobox(combobox);
   gtk_widget_set_visible(window, TRUE);

}


// User data has been placed in the window - we use it to get the
// molecule number.
void baton_mode_calculate_skeleton(GtkWidget *window) {

   std::cout << "getting intermediate data in baton_mode_calculate_skeleton "
	     << std::endl;

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));
   std::cout << "calculating map for molecule " << imol << std::endl;

   if (imol < graphics_info_t::n_molecules() && imol >= 0) {
      skeletonize_map(imol, 0);
   }
}


void change_the_contents_of_the_chain_id_combobox(GtkWidget *w, gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(w));
   GtkWidget *chain_id_combobox = widget_from_builder("renumber_residue_range_chain_id_combobox");
   // std::cout << "::::::::::: change the contents of chain_id_combobox " << chain_id_combobox
   // << " using imol " << imol << std::endl;
   GCallback null_func(NULL); // we don't do anything when the Chain ID changes. The chain-id combox is only interesting on *read*.
   g.fill_combobox_with_chain_options(chain_id_combobox, imol, null_func);
};

GtkWidget *wrapped_create_renumber_residue_range_dialog() {

   // GtkWidget *w = create_renumber_residue_range_dialog();
   GtkWidget *w = widget_from_builder("renumber_residue_range_dialog");
   GtkWidget      *mol_combobox = widget_from_builder("renumber_residue_range_molecule_combobox");
   GtkWidget *chain_id_combobox = widget_from_builder("renumber_residue_range_chain_id_combobox");

   int imol = first_coords_imol();

   graphics_info_t g;

   graphics_info_t::renumber_residue_range_molecule = imol;
   if (is_valid_model_molecule(imol)) {
      GCallback func = G_CALLBACK(change_the_contents_of_the_chain_id_combobox);
      g.new_fill_combobox_with_coordinates_options(mol_combobox, func, imol);

      GCallback null_func(NULL); // we don't do anything when the Chain ID changes. The chain-id combox is only interesting on *read*.
      // NULL is tested for in fill_combobox_with_chain_options().
      g.fill_combobox_with_chain_options(chain_id_combobox, imol, null_func);
   
      // by default, now the N-term button is off for the first choice
      // (and C-term is on for the second)
      GtkWidget *entry_1 = widget_from_builder("renumber_residue_range_resno_1_entry");
      GtkWidget *entry_2 = widget_from_builder("renumber_residue_range_resno_2_entry");

      // gtk_widget_set_sensitive(entry_2, FALSE); pre-gtk3

      // but anyway, let's put the residue number of the active residue there, just in case
      // the user wanted to start from there.
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
	 int res_no = pp.second.second.res_no;
	 gtk_editable_set_text(GTK_EDITABLE(entry_1), coot::util::int_to_string(res_no).c_str());
      }
   }
   return w;
}

bool renumber_residues_from_widget(GtkWidget *window) {

   bool status = true; // OK

   int imol = graphics_info_t::renumber_residue_range_molecule;

   if (is_valid_model_molecule(imol)) {

      GtkWidget *e1     = widget_from_builder("renumber_residue_range_resno_1_entry");
      GtkWidget *e2     = widget_from_builder("renumber_residue_range_resno_2_entry");
      GtkWidget *offent = widget_from_builder("renumber_residue_range_offset_entry");
      GtkWidget *rb1    = widget_from_builder("renumber_residue_range_radiobutton_1"); // N-term button
      GtkWidget *rb4    = widget_from_builder("renumber_residue_range_radiobutton_4"); // C-term button
      GtkWidget *chain_id_combobox = widget_from_builder("renumber_residue_range_chain_id_combobox");

      std::pair<short int, int> r1  = int_from_entry(e1);
      std::pair<short int, int> r2  = int_from_entry(e2);
      std::pair<short int, int> off = int_from_entry(offent);

      std::string chain_id = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(chain_id_combobox));
      mmdb::Chain *chain_p = graphics_info_t::molecules[imol].get_chain(chain_id);

      if (chain_p) {
	 if (gtk_check_button_get_active(GTK_CHECK_BUTTON(rb1))) {
	    // use N-terminus of chain
	    std::pair<bool, int> nt_resno = coot::util::min_resno_in_chain(chain_p);
	    if (nt_resno.first) {
	       r1.first = 1;
	       r1.second = nt_resno.second;
	    }
	 }

	 if (gtk_check_button_get_active(GTK_CHECK_BUTTON(rb4))) {
	    // use C-terminus of chain
	    std::pair<bool, int> ct_resno = coot::util::max_resno_in_chain(chain_p);
	    if (ct_resno.first) {
	       r2.first = 1;
	       r2.second = ct_resno.second;
	    }
	 }

	 if (r1.first && r2.first && off.first) {
	    int start_res = r1.second;
	    int last_res =  r2.second;
	    int offset   = off.second;

	    if (imol >= 0) {
	       if (imol < graphics_info_t::n_molecules()) {
		  if (graphics_info_t::molecules[imol].has_model()) {

		     // renumber_residue_range returns 0 upon fail
		     // including overlap, so test for this here?!
		     int status;
		     status = renumber_residue_range(imol, chain_id.c_str(), start_res, last_res, offset);
		     if (!status) {
			// error of sorts
			std::string s = "WARNING:: could not renumber residue range.\n";
			s += "Maybe your selection overlaps with existing residues.\n";
			s += "Please revise your selection.";
			info_dialog(s.c_str());
                        status = false;
		     }

		  }
	       }
	    }
	 } else {
	    std::cout << "WARNING:: Sorry. Couldn't read residue or offset from entry widget\n";
            // status = false; // hmmm.. maybe
	 }
      } else {
	 std::cout << "ERROR:: missing chain" << chain_id << std::endl;
      }
   }
   return status;
}



void apply_add_OXT_from_widget(GtkWidget *ok_button) {

   // GtkWidget *combobox = lookup_widget(ok_button, "add_OXT_molecule_combobox");
   GtkWidget *combobox = widget_from_builder("add_OXT_molecule_combobox");

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));

   std::cout << "DEBUG:: apply_add_OXT_from_widget() combobox " << combobox << " imol " << imol << std::endl;
   int resno = -9999;
   std::string chain_id = graphics_info_t::add_OXT_chain;

   GtkWidget *terminal_checkbutton = widget_from_builder("add_OXT_c_terminus_radiobutton");
   GtkWidget *residue_number_entry = widget_from_builder("add_OXT_residue_entry");

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(terminal_checkbutton))) {
      std::cout << "DEBUG:: auto determine C terminus for imol " << imol << std::endl;
      // we need to determine the last residue in this chain:
      if (is_valid_model_molecule(imol)) {
	 std::cout << "in apply_add_OXT_from_widget() here with chain_id :" << chain_id <<  ":" << std::endl;
	 graphics_info_t g;
	 std::pair<bool, int> p = g.molecules[imol].last_protein_residue_in_chain(chain_id);
	 if (p.first) {
	    resno = p.second;
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
	    if (false)
	       std::cout << "DEBUG:: adding OXT to " << imol << " "
			 << chain_id << " " << resno << std::endl;
	    add_OXT_to_residue(imol, chain_id.c_str(), resno, "");
	 }
      }
   } else {
      std::cout << "WARNING:: Could not determine last residue - not adding OXT "
		<< imol << " " << resno << "\n";
   }
}


// uses builder
GtkWidget *wrapped_create_add_OXT_dialog() {

   graphics_info_t g;

   // GtkWidget *w = create_add_OXT_dialog();
   GtkWidget *w = widget_from_builder("add_OXT_dialog");

   GtkWidget *combobox = widget_from_builder("add_OXT_molecule_combobox");
   GCallback callback_func = G_CALLBACK(g.add_OXT_molecule_combobox_changed);

   int imol = first_coords_imol();
   g.add_OXT_molecule = imol;

   if (combobox) {
      g.fill_combobox_with_coordinates_options(combobox, callback_func, imol);
      g.fill_add_OXT_dialog_internal(w, imol); // function needs object (not static)
   } else {
      std::cout << "ERROR:: in wrapped_create_add_OXT_dialog() failed to find combobox!" << std::endl;
   }

   return w;
}

void setup_alt_conf_with_dialog(GtkWidget *dialog) {

   GtkWidget *widget_ca    = widget_from_builder("add_alt_conf_ca_radiobutton");
   GtkWidget *widget_whole = widget_from_builder("add_alt_conf_whole_single_residue_radiobutton");
   GtkWidget *widget_range = widget_from_builder("add_alt_conf_residue_range_radiobutton");

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

   // GtkWidget *widget = create_add_alt_conf_dialog();
   GtkWidget *widget = widget_from_builder("add_alt_conf_dialog");
   setup_alt_conf_with_dialog(widget);
   gtk_widget_set_visible(widget, TRUE);
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
	 coot::backup_file_info_t info = graphics_info_t::molecules[imol].recent_backup_file_info();
	 if (info.valid_status) {

	    coot::backup_file_info_t *info_copy = new coot::backup_file_info_t;
	    *info_copy = info;
	    info_copy->imol = imol;

	    // GtkWidget *widget = create_recover_coordinates_dialog();
	    GtkWidget *widget = widget_from_builder("recover_coordinates_dialog");
	    g_object_set_data(G_OBJECT(widget), "backup_file_info", info_copy);

	    GtkWidget *label1, *label2;
	    label1 = widget_from_builder("recover_coordinates_read_coords_label");
	    label2 = widget_from_builder("recover_coordinates_backup_coordinates_label");

	    gtk_label_set_text(GTK_LABEL(label1), info.name.c_str());
	    gtk_label_set_text(GTK_LABEL(label2), info.backup_file_name.c_str());

	    gtk_widget_set_visible(widget, TRUE);
	    i_rec++;
	 }
      }
   }
   if (i_rec == 0) {
      // GtkWidget *w = create_nothing_to_recover_dialog();
      GtkWidget *w = widget_from_builder("nothing_to_recover_dialog");
      gtk_widget_set_visible(w, TRUE);
   }
}

// widget needed for lookup of user data:
//
void execute_recover_session(GtkWidget *widget) {

   coot::backup_file_info_t *info =
      (coot::backup_file_info_t *) g_object_get_data(G_OBJECT(widget), "backup_file_info");

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

   GtkWidget *w = widget_from_builder("merge_molecules_dialog");

   // fill the dialog here
   // GtkWidget *molecule_option_menu = lookup_widget(w, "merge_molecules_optionmenu");
   // GtkWidget *combobox = lookup_widget(w, "merge_molecules_combobox");
   // GtkWidget *molecules_vbox       = lookup_widget(w, "merge_molecules_vbox");

   GtkWidget *combobox       = widget_from_builder("merge_molecules_combobox");
   GtkWidget *molecules_vbox = widget_from_builder("merge_molecules_vbox");

   // GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(merge_molecules_menu_item_activate);
   GCallback callback_func = G_CALLBACK(merge_molecules_master_molecule_combobox_changed);

   GCallback checkbox_callback_func = G_CALLBACK(nullptr);

   graphics_info_t g;
   g.clear_out_container(molecules_vbox);

   // the molecules vbox
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

   // 20220326-PE  this should be a member function of graphics_info_t and indeed should be
   // used in fill_combobox_with_model_molecule_options which wraps
   // fill_combobox_with_molecule_options(). That function doesn't exist yet, just to be clear.
   auto get_model_molecule_vector = [] () {
                                       graphics_info_t g;
                                       std::vector<int> vec;
                                       int n_mol = g.n_molecules();
                                       for (int i=0; i<n_mol; i++)
                                          if (g.is_valid_model_molecule(i))
                                             vec.push_back(i);
                                       return vec;
                                    };

   std::vector<int> molecules_index_vec = get_model_molecule_vector();
   g.fill_combobox_with_molecule_options(combobox, callback_func, imol_master, molecules_index_vec);
   return w;
}

// void merge_molecules_menu_item_activate(GtkWidget *item,
// 					GtkPositionType pos) {
//    graphics_info_t::merge_molecules_master_molecule = pos;
// }

// #include "c-interface-gui.hh"

void merge_molecules_master_molecule_combobox_changed(GtkWidget *combobox,
						      gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::merge_molecules_master_molecule = imol;
}

void fill_vbox_with_coordinates_options(GtkWidget *dialog,
					                         GCallback checkbox_callback_func) {

   GtkWidget *molecules_vbox = widget_from_builder("merge_molecules_vbox");

   // Unset any preconceived notion of merging molecules:
   //
   // graphics_info_t::merge_molecules_merging_molecules->clear();

   gtk_box_set_spacing(GTK_BOX(molecules_vbox), 4);

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
         std::string button_label;
	      button_label = graphics_info_t::int_to_string(imol);
	      button_label += " ";
	      button_label += graphics_info_t::molecules[imol].name_for_display_manager();
	      std::string button_name = "merge_molecules_checkbutton_";
	      button_name += graphics_info_t::int_to_string(imol);

	      GtkWidget *checkbutton = gtk_check_button_new_with_label(button_label.c_str());
         g_object_set_data(G_OBJECT(checkbutton), "imol", GINT_TO_POINTER(imol));
         gtk_widget_set_visible(checkbutton, TRUE);
         gtk_box_append(GTK_BOX(molecules_vbox), checkbutton);
      }
   }
}

// The callback (if active) adds this molecule to the merging molecules list.
// If not active, it tries to remove it from the list.
//
void on_merge_molecules_check_button_toggled(GtkCheckButton *checkbutton,
					                              gpointer        user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   if (gtk_check_button_get_active(checkbutton)) {
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

   GtkWidget *w = wrapped_create_merge_molecules_dialog(); // uses builder
   gtk_widget_set_visible(w, TRUE);
}

// The action on Merge button press:
//
void do_merge_molecules(GtkWidget *dialog) {

   auto get_molecules_for_merging = [] () {
      std::vector<int> v;
      GtkWidget *box = widget_from_builder("merge_molecules_vbox");
      GtkWidget *item_widget = gtk_widget_get_first_child(box);
      while (item_widget) {
         if (gtk_check_button_get_active(GTK_CHECK_BUTTON(item_widget))) {
            int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget), "imol"));
            v.push_back(imol);
         }
         item_widget = gtk_widget_get_next_sibling(item_widget);
      };
      return v;
   };

   // std::vector<int> add_molecules = *graphics_info_t::merge_molecules_merging_molecules;
   std::vector<int> add_molecules = get_molecules_for_merging();

   if (!add_molecules.empty()) {

      if (true)
         std::cout << "calling merge_molecules_by_vector into "
                  << graphics_info_t::merge_molecules_master_molecule
                  << " n-molecules " << add_molecules.size()
                  << " starting with " << add_molecules[0]
                  << std::endl;

      std::pair<int, std::vector<merge_molecule_results_info_t> > stat =
         merge_molecules_by_vector(add_molecules, graphics_info_t::merge_molecules_master_molecule);
         if (stat.first)
            graphics_draw();
   }
}

/*  ----------------------------------------------------------------------- */
/*                         Mutate Sequence GUI                              */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_mutate_sequence_dialog() {

   // also used by wrapped_fit_loop_rama_search_dialog();

   printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- start --------------\n");

   graphics_info_t g;

   // GtkWidget *w = create_mutate_sequence_dialog();
   GtkWidget *dialog = widget_from_builder("mutate_sequence_dialog");
   printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- dialog: %p\n", dialog);
   set_transient_and_position(COOT_MUTATE_RESIDUE_RANGE_WINDOW, dialog);

   // GtkWidget *molecule_option_menu = lookup_widget(w, "mutate_molecule_optionmenu");
   // GtkWidget *chain_option_menu    = lookup_widget(w, "mutate_molecule_chain_optionmenu");
   //   GtkWidget *chain_option_menu    = lookup_widget(w, "mutate_molecule_chain_optionmenu");
   //    GtkWidget *entry1 = lookup_widget(w, "mutate_molecule_resno_1_entry");
   //    GtkWidget *entry2 = lookup_widget(w, "mutate_molecule_resno_2_entry");
   //    GtkWidget *textwindow = lookup_widget(w, "mutate_molecule_sequence_text");

   GtkWidget *combobox_molecule = widget_from_builder("mutate_sequence_molecule_combobox");
   GtkWidget *combobox_chain    = widget_from_builder("mutate_sequence_chain_combobox_text");
   // GCallback callback_func      = G_CALLBACK(mutate_sequence_molecule_menu_item_activate);
   GCallback callback_func      = G_CALLBACK(mutate_sequence_molecule_combobox_changed);

   GtkWidget *mutate_ok_button   = widget_from_builder("mutate_sequence_ok_button");
   GtkWidget *fit_loop_ok_button = widget_from_builder("fit_loop_ok_button");
   gtk_widget_set_visible(  mutate_ok_button, TRUE);
   gtk_widget_set_visible(fit_loop_ok_button, FALSE);

   printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- combobox_molecule: %p\n", combobox_molecule);
   printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- combobox_chain   : %p\n", combobox_chain);
   // Get the default molecule and fill chain combobox with the molecules chains:
   int imol = -1;
   for (int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }
   if (imol >= 0) {
      graphics_info_t::mutate_sequence_imol = imol;
      // GCallback callback = G_CALLBACK(mutate_sequence_chain_option_menu_item_activate);

      printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- calling fill_combobox_with_coordinates_options()\n");
      g.fill_combobox_with_coordinates_options(combobox_molecule, callback_func, imol);
      printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): --------------------------    done fill_combobox_with_coordinates_options()\n");

      GCallback callback = G_CALLBACK(mutate_sequence_chain_combobox_changed);
      printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- calling fill_combobox_with_chain_options()\n");
      std::string set_chain = graphics_info_t::fill_combobox_with_chain_options(combobox_chain, imol, callback);
      graphics_info_t::mutate_sequence_chain_from_combobox = set_chain;

   } else {
      graphics_info_t::mutate_sequence_imol = -1; // flag for can't mutate
   }
   return dialog;
}


void mutate_sequence_molecule_combobox_changed(GtkWidget *combobox, gpointer data) { // don't use this

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));

   graphics_info_t::mutate_sequence_imol = imol;
   GCallback chain_callback_func = G_CALLBACK(mutate_sequence_chain_combobox_changed);
   GtkWidget *chain_combobox = widget_from_builder("mutate_sequence_chain_combobox_text");
   graphics_info_t g;
   std::string set_chain = g.fill_combobox_with_chain_options(chain_combobox, imol, chain_callback_func);
   // graphics_info_t::mutate_sequence_chain_from_optionmenu = set_chain;
   graphics_info_t::mutate_sequence_chain_from_combobox = set_chain;

   printf("DEBUG:: wrapped_fit_loop_rama_search_dialog(): -------------------------- end --------------\n");
}

void mutate_sequence_molecule_menu_item_activate(GtkWidget *item,
						 GtkPositionType pos) {

   // change the chain id option menu here...
   std::cout << "DEBUG:: mutate_sequence_molecule_menu_item_activate got pos:"
	     << pos << std::endl;

   graphics_info_t::mutate_sequence_imol = pos;

   // GtkWidget *chain_option_menu = lookup_widget(item, "mutate_molecule_chain_optionmenu");

   GtkWidget *chain_combobox = widget_from_builder("mutate_molecule_chain_combobox");

   GCallback callback_func = G_CALLBACK(mutate_sequence_chain_option_menu_item_activate);

   std::string set_chain = graphics_info_t::fill_combobox_with_chain_options(chain_combobox, pos, callback_func);

   // graphics_info_t::mutate_sequence_chain_from_optionmenu = set_chain;
}

void mutate_sequence_chain_combobox_changed(GtkWidget *combobox, gpointer data) {

   char *atc = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   if (atc)
      graphics_info_t::mutate_sequence_chain_from_combobox = atc;
}

void mutate_sequence_chain_option_menu_item_activate (GtkWidget *item,
						      GtkPositionType pos) {

   // is this a callback of a combobox? I have my doubts

   // graphics_info_t::mutate_sequence_chain_from_optionmenu = menu_item_label(item);
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

   GtkWidget *entry1 = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *entry2 = widget_from_builder("mutate_molecule_resno_2_entry");

   int t;
   int res1 = -9999, res2 = -99999;
   graphics_info_t g;

   const gchar *entry_text = gtk_editable_get_text(GTK_EDITABLE(entry1));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res1 = t;
   entry_text = gtk_editable_get_text(GTK_EDITABLE(entry2));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res2 = t;

// BL says: we should set a flag that we swapped the direction and swap back
// before we call fit-gap to actually build backwards
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

   std::string chain_id = graphics_info_t::mutate_sequence_chain_from_combobox;


   // Auto fit?
   GtkWidget *checkbutton = widget_from_builder("mutate_sequence_do_autofit_checkbutton");
   short int autofit_flag = 0;

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton)))
      autofit_flag = 1;

   if (imol>= 0) { // redundant
      if (is_valid_model_molecule(imol)) {

	 // get the sequence:
	 GtkWidget *text = widget_from_builder("mutate_molecule_sequence_text");
	 char *txt = NULL;

	 GtkTextView *tv = GTK_TEXT_VIEW(text);
	 GtkTextBuffer* tb = gtk_text_view_get_buffer(tv);
	 GtkTextIter startiter;
	 GtkTextIter enditer;
	 gtk_text_buffer_get_iter_at_offset(tb, &startiter, 0);
	 gtk_text_buffer_get_iter_at_offset(tb, &enditer, -1);
	 txt = gtk_text_buffer_get_text(tb, &startiter, &enditer, 0);

	 std::string mutate_scripting_function = "mutate-and-autofit-residue-range";
	 if (! autofit_flag)
	    mutate_scripting_function = "mutate-residue-range";

	 if (txt) {
	    std::string sequence(txt);
	    sequence = coot::util::plain_text_to_sequence(sequence);
	    std::cout << "we got the sequence: " << sequence << std::endl;

	    if (int(sequence.length()) == (res2 - res1 + 1)) {
               // let's not use scripting.
	       std::vector<std::string> cmd_strings;
	       if (autofit_flag) {
                  mutate_and_autofit_residue_range(imol, chain_id.c_str(), res1, res2, sequence.c_str());
               } else {
                  mutate_residue_range(imol, chain_id.c_str(), res1, res2, sequence.c_str());
               }
               update_go_to_atom_window_on_changed_mol(imol);
               g.update_validation(imol);

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

GtkWidget *wrapped_fit_loop_rama_search_dialog() {

   GtkWidget *w = wrapped_create_mutate_sequence_dialog();

   GtkWidget *label              = widget_from_builder("function_for_molecule_label");
   GtkWidget *method_frame       = widget_from_builder("loop_fit_method_frame");
   GtkWidget *mutate_ok_button   = widget_from_builder("mutate_sequence_ok_button");
   GtkWidget *fit_loop_ok_button = widget_from_builder("fit_loop_ok_button");
   GtkWidget *checkbutton        = widget_from_builder("mutate_sequence_do_autofit_checkbutton");

   GtkWidget *rama_checkbutton   = widget_from_builder("mutate_sequence_use_ramachandran_restraints_checkbutton");

   gtk_label_set_text(GTK_LABEL(label), "\nFit loop in Molecule:\n");
   gtk_widget_set_visible(checkbutton, FALSE);
   gtk_widget_set_visible(mutate_ok_button,   FALSE);
   gtk_widget_set_visible(fit_loop_ok_button, TRUE);
   gtk_widget_set_visible(rama_checkbutton, TRUE);
   gtk_check_button_set_active(GTK_CHECK_BUTTON(rama_checkbutton), TRUE);

   gtk_widget_set_visible(method_frame, TRUE);

   return w;
}

void wrapped_fit_loop_db_loop_dialog() {

   std::vector<std::string> v;
   v.push_back("click-protein-db-loop-gui");

   if (graphics_info_t::prefer_python) {
#ifdef USE_PYTHON
      safe_python_command("import coot_gui");
      std::cout << "debug:: wrapped_fit_loop_db_loop_dialog() safe_python_command coot_gui.click_protein_db_loop_gui()"
                << std::endl;
      std::string c = graphics_info_t::pythonize_command_strings(v);
      safe_python_command("coot_gui.click_protein_db_loop_gui()");
#endif
   } else {
#ifdef USE_GUILE
      std::string c = graphics_info_t::schemize_command_strings(v);
      safe_scheme_command(c);
#endif
   }
}


/*  ----------------------------------------------------------------------- */
/*                         Align and Mutate GUI                             */
/*  ----------------------------------------------------------------------- */
GtkWidget *wrapped_create_align_and_mutate_dialog() {

   graphics_info_t g;
   // GtkWidget *w = create_align_and_mutate_dialog();
   GtkWidget *w = widget_from_builder("align_and_mutate_dialog");

   GtkWidget *mol_combobox   = widget_from_builder("align_and_mutate_molecule_combobox");
   GtkWidget *chain_combobox = widget_from_builder("align_and_mutate_chain_combobox");

   GCallback molecule_callback = G_CALLBACK(align_and_mutate_molecule_combobox_changed);
   GCallback    chain_callback = G_CALLBACK(align_and_mutate_chain_combobox_changed);

   int imol = graphics_info_t::align_and_mutate_imol;
   bool try_again = false;

   if (imol == -1) try_again = true;
   if (! is_valid_model_molecule(imol)) try_again = true;

   if (try_again) {
      for (int i=0; i<g.n_molecules(); i++) {
	 if (g.molecules[i].has_model()) {
	    imol = i;
	    break;
	 }
      }
   }

   if (imol >= 0) {
      g.fill_combobox_with_coordinates_options(mol_combobox, molecule_callback, imol);
      std::string set_chain = g.fill_combobox_with_chain_options(chain_combobox, imol,
								 chain_callback);
      g.align_and_mutate_chain_from_combobox = set_chain;
   }

   return w;
}


GtkWidget *wrapped_create_fixed_atom_dialog() {

   // GtkWidget *w = create_fixed_atom_dialog();
   // graphics_info_t::fixed_atom_dialog = w;

   GtkWidget *w = widget_from_builder("fixed_atom_dialog");
   return w;
}

#include "c-interface-gui.hh"

int do_align_mutate_sequence(GtkWidget *w) {

   //
   int handled_state = 0;  // initially unhandled (return value).

   bool renumber_residues_flag = 0; // make this derived from the GUI one day
   int imol = graphics_info_t::align_and_mutate_imol;

   GtkWidget *molecule_combobox = widget_from_builder("align_and_mutate_molecule_combobox");
   GtkWidget    *chain_combobox = widget_from_builder("align_and_mutate_chain_combobox");
   std::string chain_id = get_active_label_in_combobox(GTK_COMBO_BOX(chain_combobox));
   imol = my_combobox_get_imol(GTK_COMBO_BOX(molecule_combobox));

   GtkWidget *autofit_checkbutton = widget_from_builder("align_and_mutate_autofit_checkbutton");

   // std::cout << "--- in do_align_mutate_sequence(): combobox " << molecule_combobox
   //           << " " << GTK_IS_COMBO_BOX(molecule_combobox) << " chain-id:" << chain_id << ":"
   //           << std::endl;

   bool do_auto_fit = false;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(autofit_checkbutton)))
      do_auto_fit = true;

   graphics_info_t g;
   int imol_refinement_map = g.Imol_Refinement_Map();

   bool early_stop = false;
   if (do_auto_fit)
      if (imol_refinement_map == -1)
	 early_stop = true;

   if (early_stop) {
      std::string s = "WARNING:: autofit requested, but \n   refinement map not set!";
      std::cout << s << "\n";
      GtkWidget *warn = wrapped_nothing_bad_dialog(s);
      gtk_widget_set_visible(warn, TRUE);

   } else {

      handled_state = 1;
      if (imol >= 0) {
	 GtkWidget *text = widget_from_builder("align_and_mutate_sequence_text");
	 char *txt = NULL;

	 // text is a GtkTextView in GTK2
	 GtkTextView *tv = GTK_TEXT_VIEW(text);
	 GtkTextBuffer* tb = gtk_text_view_get_buffer(tv);
	 GtkTextIter startiter;
	 GtkTextIter enditer;
	 gtk_text_buffer_get_iter_at_offset(tb, &startiter, 0);
	 gtk_text_buffer_get_iter_at_offset(tb, &enditer, -1);
	 txt = gtk_text_buffer_get_text(tb, &startiter, &enditer, 0);

	 if (txt) {
	    std::string sequence(txt);

	    if (is_valid_model_molecule(imol)) {

	       std::cout << "debug:: calling mutate_chain " << imol << " chain-id: " << chain_id << " "
			 << sequence << " " << do_auto_fit << std::endl;
	       g.mutate_chain(imol, chain_id, sequence, do_auto_fit, renumber_residues_flag);
	       g.update_validation(imol);
	       graphics_draw();

	    }
	 }
      } else {
	 std::cout << "WARNING:: inapproproate molecule number " << imol << std::endl;
      }
   }
   return handled_state;
}


// void align_and_mutate_molecule_menu_item_activate(GtkWidget *item,
// 						  GtkPositionType pos) {

//    // GtkWidget *chain_optionmenu = lookup_widget(item, "align_and_mutate_chain_optionmenu");
//    GtkWidget *chain_combobox = lookup_widget(item, "align_and_mutate_chain_combobox");
//    GCallback chain_callback = GCallback(align_and_mutate_chain_option_menu_item_activate);
//    graphics_info_t::align_and_mutate_imol = pos;
//    int imol = pos;
//    std::string set_chain = graphics_info_t::fill_combobox_with_chain_options(chain_combobox, imol,
// 									     chain_callback);
// }

// // needs combobox version GTK-FIXME
// void align_and_mutate_chain_option_menu_item_activate (GtkWidget *item,
// 						       GtkPositionType pos) {

//    graphics_info_t::align_and_mutate_chain_from_combobox = menu_item_label(item);
//    std::cout << "align_and_mutate_chain_from_combobox is now "
// 	     << graphics_info_t::align_and_mutate_chain_from_combobox
// 	     << std::endl;
// }

void align_and_mutate_molecule_combobox_changed(GtkWidget *combobox,
						gpointer data) {

}

void align_and_mutate_chain_combobox_changed(GtkWidget *combobox,
					     gpointer data) {

}


/*  ----------------------------------------------------------------------- */
/*                  Change chain ID                                         */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_change_chain_id_dialog() {

   graphics_info_t g;

   GtkWidget *w = widget_from_builder("change_chain_id_dialog");

   // GtkWidget *mol_option_menu =  lookup_widget(w, "change_chain_id_molecule_optionmenu");
   // GtkWidget *chain_option_menu =  lookup_widget(w, "change_chain_id_chain_optionmenu");

   GtkWidget *mol_combobox   =  widget_from_builder("change_chain_id_molecule_combobox");
   GtkWidget *chain_combobox =  widget_from_builder("change_chain_id_chain_combobox");
   GtkWidget *residue_range_no_radiobutton = widget_from_builder("change_chain_residue_range_no_radiobutton");

   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(residue_range_no_radiobutton), TRUE);

   // GCallback callback_func = G_CALLBACK(change_chain_ids_mol_option_menu_item_activate);
   GCallback molecule_callback_func = G_CALLBACK(change_chain_ids_molecule_combobox_changed);

   int imol = first_coords_imol();
   if (imol >= 0) {
      g.change_chain_id_molecule = imol;
      GCallback chain_callback_func = NULL; // G_CALLBACK(change_chain_ids_chain_menu_item_activate);
      std::string set_chain = g.fill_combobox_with_chain_options(chain_combobox,
								 imol,
								 chain_callback_func);
      g.change_chain_id_from_chain = set_chain;
   }

   std::cout << "----------- fill_combobox_with_coordinates_options " << mol_combobox << std::endl;
   g.fill_combobox_with_coordinates_options(mol_combobox, molecule_callback_func, imol);
   return w;
}


// GTK3 dump
// void
// change_chain_ids_mol_option_menu_item_activate(GtkWidget *item,
// 					       GtkPositionType pos) {
//    graphics_info_t::change_chain_id_molecule = pos;
//    int imol = pos;
//    GtkWidget *chain_option_menu =  lookup_widget(item, "change_chain_id_chain_optionmenu");
//    GCallback chain_callback_func = G_CALLBACK(change_chain_ids_chain_menu_item_activate);
//    std::string set_chain = graphics_info_t::fill_combobox_with_chain_options(chain_option_menu,
// 									     imol,
// 									     chain_callback_func);
//    graphics_info_t::change_chain_id_from_chain = set_chain;
// }

void
change_chain_ids_molecule_combobox_changed(GtkWidget *combobox, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t::change_chain_id_molecule = imol;
   // GtkWidget *chain_combobox = lookup_widget(combobox, "change_chain_id_chain_combobox");
   GtkWidget *chain_combobox = widget_from_builder("change_chain_id_chain_combobox");
   if (chain_combobox) {
      graphics_info_t g;
      g.fill_combobox_with_chain_options(chain_combobox, imol, NULL);
   }
}


// // needs a combobox version
// void
// change_chain_ids_chain_menu_item_activate(GtkWidget *item,
// 					  GtkPositionType pos) {
//    graphics_info_t::change_chain_id_from_chain = menu_item_label(item);
// }


void
change_chain_id_by_widget(GtkWidget *w) {

   GtkWidget *residue_range_yes_radiobutton = widget_from_builder("change_chain_residue_range_yes_radiobutton");
   GtkWidget *residue_range_from_entry      = widget_from_builder("change_chain_residues_from_entry");
   GtkWidget *residue_range_to_entry        = widget_from_builder("change_chains_residues_to_entry");
   GtkWidget *change_chains_new_chain_entry = widget_from_builder("change_chains_new_chain_id");
   GtkWidget *change_chain_id_from_chain_combobox = widget_from_builder("change_chain_id_chain_combobox");

   int imol = graphics_info_t::change_chain_id_molecule;
   bool use_res_range_flag = false;
   int from_resno = -9999;
   int to_resno   = -9999;

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(residue_range_yes_radiobutton))) {
      use_res_range_flag = true;
      std::pair<short int, int> p1 = int_from_entry(residue_range_from_entry);
      std::pair<short int, int> p2 = int_from_entry(residue_range_to_entry);
      if (p1.first)
	 from_resno = p1.second;
      if (p2.first)
	 to_resno = p2.second;
   }

   const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(change_chains_new_chain_entry));

   if (txt) {

      if (is_valid_model_molecule(imol)) {
	 std::string to_chain_id(txt);

         // 20190810-PE use the widget to find the value now, rather than looking up a stored value
	 std::string from_chain_id = get_active_label_in_combobox(GTK_COMBO_BOX(change_chain_id_from_chain_combobox));

	 if (false)
	    std::cout << "in change_chain_id_molecule() with " << imol << " "
		      << from_chain_id << " " << to_chain_id<< std::endl;

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
	    gtk_widget_set_visible(ws, TRUE);
	 }
	 graphics_info_t g;
	 g.update_validation(imol);
      }
   } else {
      std::cout << "ERROR: Couldn't get txt in change_chain_id_by_widget\n";
   }
}

// replace this - we will find this problem again at link time, I guess.
// void fill_option_menu_with_refine_options(GtkWidget *option_menu) {
//    graphics_info_t g;
//    g.fill_option_menu_with_map_options(option_menu,
// 				       G_CALLBACK(graphics_info_t::refinement_map_select));
// }

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
   // gtk_container_remove(GTK_CONTAINER(data), widget);
   std::cout << "FIXME in my_delete_ramachandran_mol_option() " << std::endl;
}


void
show_fix_nomenclature_errors_gui(int imol,
				 const std::vector<std::pair<std::string, coot::residue_spec_t> > &nomenclature_errors) {

   return; // add this hack because this function crashes, I don't know why -fix later.

   if (graphics_info_t::use_graphics_interface_flag) {
      if (is_valid_model_molecule(imol)) {

	 // GtkWidget *w = create_fix_nomenclature_errors_dialog();
	 GtkWidget *w = widget_from_builder("fix_nomenclature_errors_dialog");
	 GtkWidget *label = widget_from_builder("fix_nomenclature_errors_label");

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

	 g_object_set_data(G_OBJECT(w), "imol", GINT_TO_POINTER(imol));

	 gtk_label_set_text(GTK_LABEL(label), s.c_str());

	 // GtkWidget *box = lookup_widget(w, "nomenclature_errors_vbox");
	 GtkWidget *box = widget_from_builder("nomenclature_errors_vbox");

	 if (box) {
	    // fill box
            std::cout << "Fix box with " << nomenclature_errors.size() << " nomenclature_errors " << std::endl;

	    for (unsigned int i=0; i<nomenclature_errors.size(); i++) {
	       s = nomenclature_errors[i].first; // the residue type
	       s += " ";
	       s += nomenclature_errors[i].second.format();
	       GtkWidget *l = gtk_label_new(s.c_str());
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
               // 20220528-PE-FIXME box packing
#else
	       gtk_box_pack_start(GTK_BOX(box), GTK_WIDGET(l), FALSE, FALSE, 2);
#endif
	       gtk_widget_set_visible(GTK_WIDGET(l), TRUE);
	    }
	 }
	 gtk_widget_set_visible(w, TRUE);

      }
   }
}


#include "get-monomer.hh"

/*  ----------------------------------------------------------------------- */
/*                  get monomer                                             */
/*  ----------------------------------------------------------------------- */

/* Get monomer code */
void
handle_get_monomer_code(GtkWidget *entry_widget) {

   GtkWidget *failed_to_get_monomer_frame = widget_from_builder("get_monomer_no_entry_frame");
   const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry_widget));

   if (! failed_to_get_monomer_frame) return;

   std::string text_s(text);
   text_s = coot::util::Upper(text_s);

   // this is set below
   int no_entry_frame_shown = gtk_widget_is_visible(failed_to_get_monomer_frame);

   if (no_entry_frame_shown == 0) { // normal case

      int imol = get_monomer(text_s);

      if (is_valid_model_molecule(imol)) {
      } else {
         gtk_widget_set_visible(failed_to_get_monomer_frame, TRUE);
      }

   } else {

      std::cout << "INFO:: handle_get_monomer_code(): Get-by-network method " << text << std::endl;

      int imol = get_monomer_molecule_by_network_and_dict_gen(text_s);
      if (! is_valid_model_molecule(imol)) {
         info_dialog("WARNING:: Failed to import molecule");
      }
   }
}


/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                               skeleton                                   */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

GtkWidget *
create_skeleton_colour_selection_window() {

   std::cout << "--------------- fix up color selection " << std::endl;
#if 0
   GtkWidget  *colorseldialog;
   GtkWidget *colorsel;

   colorseldialog =
      gtk_color_selection_dialog_new("Skeleton Colour Selection");

/* How do we get to the buttons? */

   colorsel = GTK_COLOR_SELECTION_DIALOG(colorseldialog)->colorsel;

  /* Capture "color_changed" events in col_sel_window */

  gtk_signal_connect (GTK_OBJECT (colorsel), "color_changed",
                      G_CALLBACK(on_skeleton_color_changed),
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
#endif
  return 0;
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

  GtkWidget *dialog;
  GtkWidget *helix_temp_combobox;
  GtkWidget *strand_temp_combobox;
  GtkWidget *helix_noaa_combobox;
  GtkWidget *strand_noaa_combobox;
  GtkWidget *radius_combobox;

  // dialog = create_fast_ss_search_dialog();
  dialog = widget_from_builder("fast_ss_search_dialog");

  helix_temp_combobox = widget_from_builder("fast_sss_dialog_helix_template_combobox");
  helix_noaa_combobox = widget_from_builder("fast_sss_dialog_helix_no_aa_combobox");
  strand_temp_combobox = widget_from_builder("fast_sss_dialog_strand_template_combobox");
  strand_noaa_combobox = widget_from_builder("fast_sss_dialog_strand_no_aa_combobox");
  radius_combobox = widget_from_builder("fast_sss_dialog_radius_combobox");

  // fill the comboboxes (done automatically, set the active ones)
  gtk_combo_box_set_active(GTK_COMBO_BOX(helix_temp_combobox), 0);
  gtk_combo_box_set_active(GTK_COMBO_BOX(helix_noaa_combobox), 1);
  gtk_combo_box_set_active(GTK_COMBO_BOX(strand_temp_combobox), 1);
  gtk_combo_box_set_active(GTK_COMBO_BOX(strand_noaa_combobox), 0);
  gtk_combo_box_set_active(GTK_COMBO_BOX(radius_combobox),1);

  return dialog;
}


/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
// Edit Functions that have been promoted from Extensions -> Modelling
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
void do_edit_copy_molecule() {

   auto get_molecule_vector = [] () {

      graphics_info_t g;
      std::vector<int> vec;
      int n_mol = g.n_molecules();
      for (int i=0; i<n_mol; i++) {
	 if (g.is_valid_model_molecule(i))
	    vec.push_back(i);
	 if (g.is_valid_map_molecule(i))
	    vec.push_back(i);
      }
      return vec;
   };

   GtkWidget *frame    = widget_from_builder("copy-molecule-frame");
   GtkWidget *combobox = widget_from_builder("copy_molecule_comboboxtext");
   gtk_widget_set_visible(frame, TRUE);
   auto molecule_list = get_molecule_vector();
   int imol_active = -1;
   GCallback func = G_CALLBACK(nullptr); // we don't care until this dialog is read
   graphics_info_t g;
   g.fill_combobox_with_molecule_options(combobox, func, imol_active, molecule_list);

}

void  do_edit_copy_fragment() {

   graphics_info_t g;
   GtkWidget *frame = widget_from_builder("copy_fragment_frame");
   GtkWidget *vbox   = widget_from_builder("copy_fragment_vbox");
   int imol_active = g.get_active_atom().first;

   auto my_delete_box_items = [] (GtkWidget *widget, void *data) {
   };


   GtkWidget *combobox_molecule = widget_from_builder("copy_fragment_combobox"); // its a GtkComboBoxText
   GCallback callback_func = G_CALLBACK(NULL); // combobox is only used when it's read on OK response

   // new_fill_combobox_with_coordinates_options() doesn't set the active item - I don't understand why.
   g.new_fill_combobox_with_coordinates_options(combobox_molecule, callback_func, imol_active);
   g_object_set_data(G_OBJECT(frame), "combobox", combobox_molecule); // for reading. 20220828-PE still needed?
   set_transient_for_main_window(frame);
   gtk_widget_set_visible(frame, TRUE);

   // the dialog response callback for this is on_copy_fragment_dialog_response()

}


// 20240930-PE remove the usage of scripting from this function
// Make it an overlay
void do_edit_replace_fragment() {

   // Molecule Working     [molecule chooser] # needs updating
   // Molecule Reference   [molecule chooser] # contains the fragment to be copied
   // Atom Selection       [________________]
   //                        Cancel   Replace


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

#ifdef USE_GUILE
   std::string cmd =
      "(molecule-chooser-gui \"Define the molecule that needs updating\" (lambda (imol-base) (generic-chooser-and-entry \"Molecule that contains the new fragment:\" \"Atom Selection\" \"//\" (lambda (imol-fragment atom-selection-str) (replace-fragment imol-base imol-fragment atom-selection-str)))))";
   if (state_lang == coot::STATE_SCM) {
      safe_scheme_command(cmd);
   }
#else
#ifdef USE_PYTHON
   if (state_lang == coot::STATE_PYTHON) {

      std::string cmd =
         "import coot_gui\ncoot_gui.molecule_chooser_gui(\"Define the molecule that needs updating\", lambda imol_base: coot_gui.generic_chooser_and_entry(\"Molecule that contains the new fragment:\", \"Atom Selection\", \"//\", lambda imol_fragment, atom_selection_str: coot.replace_fragment(imol_base, imol_fragment, atom_selection_str)))";
      safe_python_command(cmd);
   }
#endif // PYTHON
#endif // GUILE
}
