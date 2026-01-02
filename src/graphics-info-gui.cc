/* src/graphics-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
 * Copyright 2014, 2016 by Medical Research Council
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
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <algorithm>

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
// #include <GL/glut.h>  // for some reason...  // Eh?

#include <iostream>
#ifdef _MSC_VER
#include "compat/dirent.h"
#else
#include <dirent.h>   // for refmac dictionary files
#endif

#include <string.h> // strncpy

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER && !defined WINDOWS_MINGW
#include <unistd.h>
#endif

#include <clipper/core/map_utils.h> // Map_stats
#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"

#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Bond_lines.hh"

#include "skeleton/graphical_skel.h"

#include "graphics-info.h"
#include "interface.h"

#include "molecule-class-info.h"
#include "skeleton/BuildCas.h"

#include "coot-utils/gl-matrix.h" // for baton rotation

#include "analysis/bfkurt.hh"

#include "globjects.h"
#include "ligand/ligand.hh"

#include "ligand/dunbrack.hh"

#include "utils/coot-utils.hh"

//temp
#include "cmtz-interface.hh"

#include "manipulation-modes.hh"
#include "guile-fixups.h"
#include "widget-from-builder.hh"
#include "c-interface-gtk-widgets.h"

#include "utils/logging.hh"
extern logging logger;

void do_accept_reject_dialog(std::string fit_type, const coot::refinement_results_t &rr) {

   do_accept_reject_hud_buttons(fit_type, rr); // not that we use the args (yet?)

}

void do_accept_reject_hud_buttons(std::string fit_type, const coot::refinement_results_t &rr) {

   graphics_info_t g;
   g.show_atom_pull_toolbar_buttons();
   g.show_accept_reject_hud_buttons();

}


void
graphics_info_t::show_refinement_and_regularization_parameters_frame() {

   GtkWidget *frame = widget_from_builder("refinement_and_regularization_parameters_frame");
   auto visible = gtk_widget_get_visible(frame);
   gtk_widget_set_visible(frame,!visible);

   // GtkWidget *vbox_container = widget_from_builder("refinement_and_regularization_vbox_container");
   // GtkWidget *vbox_outer     = widget_from_builder("refinement_and_regularization_vbox_outer"); // inside vbox_container

   // we don't want to gtk_overlay_add_overlay() for vbox_outer if it has already been added.
   // The proxy for testing if that is the case is testing if it is realized:
   //
   // if (gtk_widget_get_realized(vbox_outer)) {
   //    // std::cout << "vbox_outer already realized" << std::endl;
   // } else {

   //    // 20211027-PE this is how the old dialog was filled.
   //    // set_refine_params_toggle_buttons(dialog);


   set_refine_params_comboboxes();

   GtkWidget *overall_weight_combobox = widget_from_builder("refine_params_overall_weight_combobox");
   gtk_combo_box_text_remove_all(GTK_COMBO_BOX_TEXT(overall_weight_combobox));
   std::vector<float> mv = {0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 20.0};
   graphics_info_t g;
   float w = g.geometry_vs_map_weight;
   for (auto m : mv) {
      std::string t = coot::util::float_to_string_using_dec_pl(w * m, 2);
      gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(overall_weight_combobox), t.c_str());
   }
   gtk_combo_box_set_active(GTK_COMBO_BOX(overall_weight_combobox), 4);

   GtkWidget *use_torsions_checkbutton = widget_from_builder("refine_params_use_torsions_checkbutton");
   GtkWidget *use_planepep_checkbutton = widget_from_builder("refine_params_use_planar_peptides_checkbutton");
   GtkWidget *use_transpep_checkbutton = widget_from_builder("refine_params_use_trans_peptide_restraints_checkbutton");
   GtkWidget *use_rama_restr_checkbutton = widget_from_builder("refine_params_use_ramachandran_goodness_torsions_checkbutton");

   if (false) {
      std::cout << "debug:: do_torsions " << g.do_torsion_restraints << std::endl;
      std::cout << "debug:: do_trans_peptide_restraints " << g.do_trans_peptide_restraints << std::endl;
      std::cout << "debug:: planar peptides " << Geom_p()->planar_peptide_restraint_state() << std::endl;
   }

   if (g.do_torsion_restraints)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_torsions_checkbutton), TRUE);
   else
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_torsions_checkbutton), FALSE);

   if (g.do_trans_peptide_restraints)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_transpep_checkbutton), TRUE);
   else
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_transpep_checkbutton), FALSE);

   if (Geom_p()->planar_peptide_restraint_state())
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_planepep_checkbutton), TRUE);
   else
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_planepep_checkbutton), FALSE);

   if (g.do_rama_restraints)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_rama_restr_checkbutton), TRUE);
   else
      gtk_check_button_set_active(GTK_CHECK_BUTTON(use_rama_restr_checkbutton), FALSE);

   // } //vbox outer not realized
}


// static
// 20230218-PE webassembly merge: make this void
void
graphics_info_t::info_dialog(const std::string &s, bool use_markup) {
   wrapped_nothing_bad_dialog(s, use_markup);
}

void
graphics_info_t::info_dialog_and_text(const std::string &s, bool use_markup) {

  if (graphics_info_t::use_graphics_interface_flag) {
     info_dialog(s, use_markup);
   }
   std::cout<< s <<std::endl;
}

// This uses an info_dialog but uses gtk_label_set_markup() to display
// in a monospace font.
void
graphics_info_t::info_dialog_alignment(coot::chain_mutation_info_container_t mutation_info) const {

   if (! use_graphics_interface_flag) return;

   std::string s = mutation_info.alignment_string;

   info_dialog(s); // get trashed by markup text
   // GtkWidget *label = lookup_widget(dialog, "nothing_bad_label");
   GtkWidget *label = widget_from_builder("nothing_bad_label");
   gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_LEFT);

   // guessing that we need > 6, could be more than 6.
   gtk_label_set_markup(GTK_LABEL(label), s.c_str());
}

void
graphics_info_t::info_dialog_refinement_non_matching_atoms(std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > nma) {

   std::string s = "WARNING:: Failed to match (to the dictionary) the following model atom names:\n";
   for (unsigned int i=0; i<nma.size(); i++) {
      s += "   ";
      s += nma[i].first->GetChainID();
      s += " ";
      s += int_to_string(nma[i].first->GetSeqNum());
      s += " ";
      s += nma[i].first->GetResName();
      s += "  ";
      for (unsigned int j=0; j<nma[i].second.size(); j++) {
	 s += "\"";
	 s += nma[i].second[j];
	 s += "\"   ";
      }
      s += "\n";
   }

   if(! nma.empty()) {
      s += "\n";
      s += "   That would cause exploding atoms, so the refinement didn't start\n";
   }

   info_dialog(s); // get trashed by markup text
   // GtkWidget *label = lookup_widget(dialog, "nothing_bad_label");
   GtkWidget *label = widget_from_builder("nothing_bad_label");
   gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_LEFT);
}



// To be used to (typically) get the menu item text label from chain
// option menus (rather than the ugly/broken casting of
// GtkPositionType data.
// static
std::string
graphics_info_t::menu_item_label(GtkWidget *menu_item) {

   // do we need this function? try commenting it out :-)

   gchar *text = 0;
#if 0
   if (GTK_BIN (menu_item)->child) {
      GtkWidget *child = GTK_BIN (menu_item)->child;
      if (GTK_IS_LABEL (child)) {
	 gtk_label_get (GTK_LABEL (child), &text);
	 // g_print ("menu item text: %s\n", text);
      }
   }
#endif
   if (text) {
      std::string s(text);
      return s;
   } else {
      return std::string("");
   }
}




// static
void
graphics_info_t::set_transient_and_position(int widget_type, GtkWidget *window) {

   GtkWindow *main_window = GTK_WINDOW(get_main_window());
   gtk_window_set_transient_for(GTK_WINDOW(window), main_window);

   if (widget_type == COOT_EDIT_CHI_DIALOG) {

      if (graphics_info_t::edit_chi_angles_dialog_x_position > -100) {
	 if (graphics_info_t::edit_chi_angles_dialog_y_position > -100) {
	    std::cout << "GTK3 FIXME set_transient_and_position() no gtk_widget_set_uposition" << std::endl;
	    // gtk_widget_set_uposition(window,
	    // graphics_info_t::edit_chi_angles_dialog_x_position,
	    //graphics_info_t::edit_chi_angles_dialog_y_position);
	 }
      }
   }


   // perhaps we should do something like this for all transient dialogs, hmmm...
   //
   if (widget_type == COOT_ROTAMER_SELECTION_DIALOG) {
      bool done_set_pos = false;
      if (graphics_info_t::rotamer_selection_dialog_x_position > -100) {
	 if (graphics_info_t::rotamer_selection_dialog_y_position > -100) {
	    std::cout << "GTK3 FIXME set_transient_and_position() no gtk_widget_set_uposition" << std::endl;
	    //gtk_widget_set_uposition(window,
	    // graphics_info_t::rotamer_selection_dialog_x_position,
	    // graphics_info_t::rotamer_selection_dialog_y_position);
	    done_set_pos = true;
	 }
      }
      if (! done_set_pos) {
	 int x_pos = graphics_info_t::graphics_x_position - 100;
	 int y_pos = graphics_info_t::graphics_y_position + 100;
	 if (x_pos < 5) x_pos = 5;

	    std::cout << "GTK3 FIXME set_transient_and_position() no gtk_widget_set_uposition" << std::endl;
	 // gtk_widget_set_uposition(window, x_pos, y_pos);
      }
   }
}


// static
void
graphics_info_t::add_status_bar_text(const std::string &text) {

   if (use_graphics_interface_flag) {

      GtkWidget *statusbar = widget_from_builder("main_window_statusbar");
      if (statusbar) {
	 std::string sbt = text;

         // std::cout << "pushing statusbar text: " << sbt << std::endl;
         gtk_statusbar_push(GTK_STATUSBAR(statusbar), statusbar_context_id, sbt.c_str());

      } else {
         std::cout << "no statusbar" << std::endl;
      }
   }
}



void
graphics_info_t::set_directory_for_fileselection_string(std::string filename) {

   directory_for_fileselection = filename;
}


void
graphics_info_t::save_directory_from_filechooser(const GtkWidget *filechooser) {

   if (filechooser) {
      if (GTK_IS_FILE_CHOOSER(filechooser)) {
	 // gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(filechooser));
         GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(filechooser));
         GError *error = NULL;
         GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                                  G_FILE_QUERY_INFO_NONE, NULL, &error);
         const char *filename = g_file_info_get_name(file_info);
         
	 if (filename) {
	    directory_for_filechooser = coot::util::file_name_directory(filename);
	    // g_free(filename);
	 }
      }
   }
}

void
graphics_info_t::save_directory_for_saving_from_filechooser(const GtkWidget *filechooser) {

   if (GTK_IS_FILE_CHOOSER(filechooser)) {
      // gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(filechooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *filename = g_file_info_get_name(file_info);

      if (filename) {
         directory_for_saving_for_filechooser = coot::util::file_name_directory(filename);
         // g_free(filename);
      }
   }
}

void
graphics_info_t::set_directory_for_filechooser_string(std::string filename) {

   directory_for_filechooser = filename;
}

void
graphics_info_t::set_directory_for_filechooser(GtkWidget *filechooser) const {

   if (directory_for_filechooser != "") {
      std::cout << "INFO:: set directory_for_filechooser " << directory_for_filechooser << std::endl;

      // 20220602-PE FIXME
      std::cout << "INFO:: in set_directory_for_filechooser() FIXME" << std::endl;
      GFile *f = g_file_new_for_path(directory_for_filechooser.c_str());
      GError *err = NULL;
      gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filechooser), f, &err);

   } else {
      std::string cwd = coot::util::current_working_dir();
      GFile *f = g_file_new_for_path(cwd.c_str());
      GError *err = NULL;
      gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filechooser), f, &err);
   }
}

void
graphics_info_t::set_file_for_save_filechooser(GtkWidget *filechooser) const {

   // just like set_directory_for_filechooser actually, but we give
   // it the full filename, not just the directory.

   int imol = save_imol;
   if (imol >= 0 && imol < n_molecules()) {
      std::string stripped_name = molecules[imol].stripped_save_name_suggestion();
      std::string full_name = stripped_name;

      if (directory_for_saving_for_filechooser != "") {
	full_name = directory_for_saving_for_filechooser + stripped_name;
      } else {

	// if we have a directory in the filechooser path we take this
	if (directory_for_saving_for_filechooser != "") {

	} else {
	  // otherwise we make one
	  gchar *current_dir = g_get_current_dir();
	  full_name = g_build_filename(current_dir, stripped_name.c_str(), NULL);
	  directory_for_saving_for_filechooser = current_dir;
	  g_free(current_dir);
	}
      }

      if (true)
	 std::cout << "DEBUG:: Setting filechooser with file: " << full_name
		   << std::endl;

      // https://docs.gtk.org/gtk3/method.FileChooser.set_current_name.html
      //
      // use gtk_file_chooser_set_current_name() if the file is "made up" e.g. "Untitled.doc"
      // and use gtk_file_chooser_set_file() if the file exists.

      if (g_file_test(full_name.c_str(), G_FILE_TEST_EXISTS)) {

         GFile *file = g_file_new_for_path(full_name.c_str());
         GError *error = 0;
         gtk_file_chooser_set_file(GTK_FILE_CHOOSER(filechooser), file, &error); // deprecated
         g_object_unref(file);

      } else {

         GFile *file_dir = g_file_new_for_path(directory_for_saving_for_filechooser.c_str());
         GError *error = 0;

         // this is the right directory, but the dialog gives us a warning messge. Don't understand
         // std::cout << "not using file dir for " << directory_for_saving_for_filechooser << std::endl;
         // gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(filechooser), file_dir, &error);

         gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(filechooser), stripped_name.c_str());
         g_object_unref(file_dir);
      }
   }
}

void
graphics_info_t::set_scrollable_map(int imol) {

   scroll_wheel_map = imol;
   activate_scroll_radio_button_in_display_manager(imol);

}
   

// i find this somewhat asthetically pleasing (maybe because the
// display control widgets are uniquely named [which was a bit of a
// struggle in C, I seem to recall]).
//
// scroll_wheel_map must be set correctly before coming here.
//
// static
void
graphics_info_t::activate_scroll_radio_button_in_display_manager(int imol) {

//    graphics_info_t g;
//    if (g.display_control_window()) {
//       std::string wname = "map_scroll_button_";
//       wname += graphics_info_t::int_to_string(imol);
//       GtkWidget *w = lookup_widget(g.display_control_window(), wname.c_str());
//       if (w) {
// 	 if (g.scroll_wheel_map == imol) {
// 	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
// 	 }
//       }
//    }

   auto get_scrollable_map_checkbutton = [] () {

      GtkWidget *w = nullptr;
      GtkWidget *display_map_vbox = widget_from_builder("display_map_vbox");
      if (display_map_vbox) {
         GtkWidget *item_widget = gtk_widget_get_first_child(display_map_vbox);
         while (item_widget) {
            // the child is a hbox
            GtkWidget *item_widget_inner = gtk_widget_get_first_child(item_widget);
            while (item_widget_inner) {
               // the checkbuttons are the Display checkbutton or the Scroll checkbutton
               if (GTK_IS_CHECK_BUTTON(item_widget_inner)) {
                  const char *l = gtk_check_button_get_label(GTK_CHECK_BUTTON(item_widget_inner));
                  if (l) {
                     std::string ls(l);
                     if (ls == "Scroll") {
                        int imol_button = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget_inner), "imol"));
                        if (imol_button == scroll_wheel_map) {
                           w = item_widget_inner;
                        }
                     }
                  }
               }
               item_widget_inner = gtk_widget_get_next_sibling(item_widget_inner);
            };
            item_widget = gtk_widget_get_next_sibling(item_widget);
         }
      }
      return w;
   };


   graphics_info_t g;
   if (use_graphics_interface_flag) {
      GtkWidget *w = get_scrollable_map_checkbutton();
      if (w) {
         // it's a radio button
         gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
      }
   }
}

float
graphics_info_t::get_estimated_map_weight(int imol_map) {

   float v = -1.0; // invalid value (negative)

   if (is_valid_map_molecule(imol_map)) {
      float mean = graphics_info_t::molecules[imol_map].map_mean();
      float sd   = graphics_info_t::molecules[imol_map].map_sigma();

      v = 50*0.3/sd;
      if (graphics_info_t::molecules[imol_map].is_EM_map())
	 v *= 0.35;
   }
   return v;
}

void
on_select_fitting_map_dialog_estimate_weight_button_clicked(GtkButton *button, gpointer user_data) {

   GtkWidget *entry = GTK_WIDGET(user_data);
   graphics_info_t g;
   float e = g.get_estimated_map_weight(g.Imol_Refinement_Map());
   std::string t = coot::util::float_to_string_using_dec_pl(e, 2);
   g.geometry_vs_map_weight = e;

   gtk_editable_set_text(GTK_EDITABLE(entry), t.c_str());
}

// static
void
graphics_info_t::select_refinement_map_combobox_changed(GtkWidget *combobox, gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
   g.set_refinement_map(imol);

}

void
graphics_info_t::show_select_map_frame() {

   show_select_map_frame_gtkbuilder();

}

void
graphics_info_t::show_select_map_frame_gtkbuilder() {
   std::cout << "here in show_select_map_dialog_gtkbuilder() " << std::endl;

   if (use_graphics_interface_flag) {

      GtkWidget *frame = get_widget_from_builder("select_map_for_fitting_frame");

      int imol_map = Imol_Refinement_Map();

      // If no map has been set before, set the map to the top of the
      // list (if there are maps in the list)
      //
      if (imol_map == -1) {
         for (int imol=0; imol<n_molecules(); imol++) {
            if (molecules[imol].has_xmap()) {
               imol_refinement_map = imol;
               break;
            }
         }
      }

      GtkWidget *combobox = get_widget_from_builder("select_map_for_fitting_combobox");
      GCallback callback = G_CALLBACK(select_refinement_map_combobox_changed);
      fill_combobox_with_map_options(combobox, callback, imol_refinement_map);

      // now fill the weight entry:
      GtkWidget *weight_entry = get_widget_from_builder("select_fitting_map_dialog_weight_entry");
      std::string ws = coot::util::float_to_string_using_dec_pl(geometry_vs_map_weight, 3);
      if (weight_entry)
         gtk_editable_set_text(GTK_EDITABLE(weight_entry), ws.c_str());
      else
         std::cout << "ERROR:: show_select_map_dialog_gtkbuilder() failed to get weight entry" << std::endl;

      // add a callback for the weight estimate button
      GtkWidget *estimate_buton = get_widget_from_builder("select_fitting_map_dialog_estimate_button");
      if (estimate_buton) {
         g_signal_connect(G_OBJECT(estimate_buton), "clicked",
                          G_CALLBACK(on_select_fitting_map_dialog_estimate_weight_button_clicked),
                          weight_entry);
      } else {
         std::cout << "ERROR:: show_select_map_dialog_gtkbuilder() failed to get estimate button" << std::endl;
      }

      // set_transient_for_main_window(dialog);
      gtk_widget_set_visible(frame, TRUE);


   }
}

void
on_dialog_box_of_buttons_response(GtkDialog *dialog,
                                  gint       response_id,
                                  gpointer   user_data) {
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

GtkWidget *
graphics_info_t::dialog_box_of_buttons_internal(const std::string &window_title,
                                                const std::vector<std::tuple<std::string, GCallback, gpointer> > &buttons,
                                                const std::string &close_button_label) {


   std::string full_title = std::string("Coot: ") + window_title;
   GtkWidget *dialog = gtk_dialog_new();
   GtkWidget *scrolled_window = gtk_scrolled_window_new();
   gtk_window_set_default_size(GTK_WINDOW(dialog), 180, 300);
   gtk_window_set_title(GTK_WINDOW(dialog), full_title.c_str());
   GtkWidget *vbox_outer = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
   GtkWidget *vbox       = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
   GtkWidget *close_button = gtk_button_new_with_label(close_button_label.c_str());
   for (unsigned int i=0; i<buttons.size(); i++) {
      GtkWidget *button  = gtk_button_new_with_label(std::get<0>(buttons[i]).c_str());
      GCallback callback = std::get<1>(buttons[i]);
      gtk_box_append(GTK_BOX(vbox), button);
      gtk_widget_set_margin_start (button, 4);
      gtk_widget_set_margin_end   (button, 4);
      gtk_widget_set_margin_top   (button, 4);
      gtk_widget_set_margin_bottom(button, 4);
      g_signal_connect(G_OBJECT(button), "clicked", callback, std::get<2>(buttons[i]));
      gtk_widget_set_visible(button, TRUE);
   }
   gtk_box_append(GTK_BOX(vbox_outer), scrolled_window);
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), vbox);
   gtk_widget_set_vexpand(vbox, TRUE);               // not sure that all of these are needed.
   gtk_widget_set_hexpand(scrolled_window, TRUE);    //
   gtk_widget_set_vexpand(scrolled_window, TRUE);    //

   gtk_widget_set_margin_start (vbox_outer, 4);
   gtk_widget_set_margin_end   (vbox_outer, 4);
   gtk_widget_set_margin_top   (vbox_outer, 4);
   gtk_widget_set_margin_bottom(vbox_outer, 4);

   g_signal_connect(G_OBJECT(dialog), "response", G_CALLBACK(on_dialog_box_of_buttons_response), dialog);
   gtk_dialog_add_button(GTK_DIALOG(dialog), close_button_label.c_str(), GTK_RESPONSE_CLOSE);

   gtk_widget_set_visible(scrolled_window, TRUE);
   gtk_widget_set_visible(vbox, TRUE);
   gtk_widget_set_visible(vbox_outer, TRUE);
   gtk_widget_set_visible(close_button, TRUE);
   set_transient_for_main_window(dialog);
   return dialog;
}



GtkWidget *
graphics_info_t::wrapped_create_skeleton_dialog(bool show_ca_mode_needs_skel_label) {

   // GtkWidget *w = create_skeleton_dialog();
   // GtkWidget *combobox    = lookup_widget(w, "skeleton_map_combobox");
   // GtkWidget *frame = lookup_widget(w, "skeleton_dialog_on_off_frame");
   // GtkWidget *label = lookup_widget(w, "ca_baton_mode_needs_skel_label");
   // GtkWidget *ok_button = lookup_widget(w, "skeleton_ok_button");
   GtkWidget *w         = widget_from_builder("skeleton_dialog");
   GtkWidget *combobox  = widget_from_builder("skeleton_map_combobox");
   GtkWidget *frame     = widget_from_builder("skeleton_dialog_on_off_frame");
   GtkWidget *label     = widget_from_builder("ca_baton_mode_needs_skel_label");
   GtkWidget *ok_button = widget_from_builder("skeleton_ok_button");

   // add user data to the OK button, we use it to determine if we go
   // on to display the baton dialog
   int show_baton_dialog = 0;

   if (show_ca_mode_needs_skel_label) {
      show_baton_dialog = 1;
   }

   g_signal_connect(G_OBJECT(ok_button),
		    "clicked",
		    G_CALLBACK(on_skeleton_ok_button_dynamic_clicked),
		    GINT_TO_POINTER(show_baton_dialog));

   if (show_ca_mode_needs_skel_label) {
      gtk_widget_set_visible(label, TRUE);
   }
   set_initial_map_for_skeletonize();
   fill_combobox_with_skeleton_options(combobox);
   set_on_off_skeleton_radio_buttons(frame);
   return w;
}


// static
void
graphics_info_t::on_skeleton_ok_button_dynamic_clicked (GtkButton       *button,
							gpointer         user_data) {

   // GtkWidget *window = lookup_widget(GTK_WIDGET(button), "skeleton_dialog");
   GtkWidget *window = widget_from_builder("skeleton_dialog");

   // GtkWidget *optionmenu = lookup_widget(window, "skeleton_map_optionmenu");
   GtkWidget *combobox = widget_from_builder("skeleton_map_combobox");

   int do_baton_mode = GPOINTER_TO_INT(user_data);

   std::cout << "do_baton_mode:: " << do_baton_mode << std::endl;

   graphics_info_t g;
   // g.skeletonize_map_by_optionmenu(optionmenu);
   g.skeletonize_map_by_combobox(combobox);
   gtk_widget_set_visible(window, FALSE);
   if (do_baton_mode) {
      int state = g.try_set_draw_baton(1);
      if (state) {
	 //  GtkWidget *w = create_baton_dialog();
	 GtkWidget *w = widget_from_builder("baton_dialog");
	 gtk_widget_set_visible(w, TRUE);
      }
   }
}

// void
// graphics_info_t::skeletonize_map_by_optionmenu(GtkWidget *optionmenu) {

//    GtkWidget *window = lookup_widget(GTK_WIDGET(optionmenu), "skeleton_dialog");

//    GtkWidget *on_radio_button = NULL;
//    GtkWidget *prune_check_button = NULL;

//    on_radio_button = lookup_widget(window, "skeleton_on_radiobutton");

//    short int do_it = 0;
//    short int prune_it = 0;
//    if (! is_valid_map_molecule(graphics_info_t::map_for_skeletonize)) {
//       std::cout << "ERROR:: Trapped a bad map for skeletoning!" << std::endl;
//    } else {
//       if (GTK_TOGGLE_BUTTON(on_radio_button)->active) {
// 	 do_it = 1;
//       }
//       prune_check_button = lookup_widget(window,"skeleton_prune_and_colour_checkbutton");
//       if (GTK_TOGGLE_BUTTON(prune_check_button)->active) {
// 	 prune_it = 1;
//       }

//       if (do_it)
// 	 graphics_info_t::skeletonize_map(graphics_info_t::map_for_skeletonize, prune_it);
//       else {
// 	 std::cout << "INFO:: unskeletonizing map number "
// 		   << graphics_info_t::map_for_skeletonize << std::endl;
// 	 graphics_info_t::unskeletonize_map(graphics_info_t::map_for_skeletonize);
//       }
//    }
// }

void
graphics_info_t::skeletonize_map_by_combobox(GtkWidget *combobox) {

   // GtkWidget *window = lookup_widget(GTK_WIDGET(combobox), "skeleton_dialog");
   GtkWidget *window = widget_from_builder("skeleton_dialog");

   GtkWidget *on_radio_button = NULL;
   GtkWidget *prune_check_button = NULL;

   on_radio_button = widget_from_builder("skeleton_on_radiobutton");

   short int do_it = 0;
   short int prune_it = 0;
   if (! is_valid_map_molecule(map_for_skeletonize)) {
      std::cout << "ERROR:: Trapped a bad map for skeletoning!" << std::endl;
   } else {
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(on_radio_button))) {
	 do_it = 1;
      }
      // prune_check_button = lookup_widget(window,"skeleton_prune_and_colour_checkbutton");
      prune_check_button = widget_from_builder("skeleton_prune_and_colour_checkbutton");
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(prune_check_button))) {
	 prune_it = 1;
      }

      if (do_it)
	 graphics_info_t::skeletonize_map(graphics_info_t::map_for_skeletonize, prune_it);
      else {
	 std::cout << "INFO:: unskeletonizing map number "
		   << graphics_info_t::map_for_skeletonize << std::endl;
	 graphics_info_t::unskeletonize_map(graphics_info_t::map_for_skeletonize);
      }
   }
}




int
graphics_info_t::try_set_draw_baton(short int i) {
   graphics_info_t g;
   if (i) {
      bool have_skeled_map_state = g.start_baton_here();
      if (have_skeled_map_state)
	 g.draw_baton_flag = 1;
   } else {
      g.draw_baton_flag = 0;
   }
   graphics_draw();
   return g.draw_baton_flag;
}


void
graphics_info_t::fill_combobox_with_skeleton_options(GtkWidget *combobox) {

   graphics_info_t g;
   GCallback signalfunc = G_CALLBACK(skeleton_map_combobox_changed);
   fill_combobox_with_map_options(combobox, signalfunc, imol_refinement_map);
}

void
graphics_info_t::skeleton_map_combobox_changed(GtkWidget *combobox, gpointer data) {

   graphics_info_t g;
   g.map_for_skeletonize = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
}



#if 0
void
graphics_info_t::fill_option_menu_with_coordinates_options(GtkWidget *option_menu,
							   GtkSignalFunc signal_func,
							   int imol_active_position) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
//    std::cout << "option_menu: " << option_menu << std::endl;
//    std::cout << "menu: " << menu << std::endl;
//    std::cout << "menu is widget: " << GTK_IS_WIDGET(menu) << std::endl;
//    std::cout << "menu is menu: " << GTK_IS_MENU(menu) << std::endl;

   // menu is not GTK_MENU on Gtk2 Ubuntu kalypso 64 bit
   if (GTK_IS_MENU(menu))
      gtk_widget_destroy(menu);

   menu = gtk_menu_new();

   GtkWidget *menuitem;
   int item_count = 0;

   for (int imol=0; imol<n_molecules(); imol++) {

      if (molecules[imol].has_model() > 0) {

	 std::string ss = int_to_string(imol);
	 ss += " " ;
	 int ilen = molecules[imol].name_.length();
	 int left_size = ilen-go_to_atom_menu_label_n_chars_max;
	 if (left_size <= 0) {
	    // no chop
	    left_size = 0;
	 } else {
	    // chop
	    ss += "...";
	 }

	 ss += molecules[imol].name_.substr(left_size, ilen);

	 menuitem = gtk_menu_item_new_with_label (ss.c_str());

	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     signal_func,
			     GINT_TO_POINTER(imol));

	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 gtk_widget_set_visible(menuitem, TRUE);
 	 if (imol == imol_active_position) {
 	    gtk_menu_set_active(GTK_MENU(menu), item_count);
	 }
	 item_count++;
      }
   }

//    gtk_menu_set_active(GTK_MENU(menu), imol_active_position);

   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu),
			    menu);

}
#endif

#if 0
void
graphics_info_t::fill_option_menu_with_coordinates_options_possibly_small(GtkWidget *option_menu,
									  GtkSignalFunc callback_func,
									  int imol_active,
									  bool fill_with_small_molecule_only_flag) {

   int n_atoms_means_big_molecule = N_ATOMS_MEANS_BIG_MOLECULE;  // 400
   short int set_last_active_flag = 0;
   std::vector<int> fill_with_these_molecules;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model()) {
	 int n_atoms = molecules[imol].atom_sel.n_selected_atoms;
	 if ( (!fill_with_small_molecule_only_flag) || (n_atoms < n_atoms_means_big_molecule) )
	    fill_with_these_molecules.push_back(imol);
      }
   }
   fill_option_menu_with_coordinates_options_internal_3(option_menu, callback_func,
							fill_with_these_molecules,
							set_last_active_flag, imol_active);
}
#endif



void
graphics_info_t::set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame) {

   // GtkWidget *on_button = lookup_widget(skeleton_frame, "skeleton_on_radiobutton");
   // GtkWidget *off_button = lookup_widget(skeleton_frame, "skeleton_off_radiobutton");

   GtkWidget *on_button = widget_from_builder("skeleton_on_radiobutton");
   GtkWidget *off_button = widget_from_builder("skeleton_off_radiobutton");

   int imol = map_for_skeletonize;
   if (imol >= 0) {
      if (molecules[imol].xskel_is_filled) {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button),  TRUE);
      } else {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), TRUE);
      }
   }
}

void
graphics_info_t::set_on_off_single_map_skeleton_radio_buttons(GtkWidget *skeleton_frame,
							      int imol) {
   // GtkWidget *on_button = lookup_widget(skeleton_frame, "single_map_skeleton_on_radiobutton");
   // GtkWidget *off_button = lookup_widget(skeleton_frame, "single_map_skeleton_off_radiobutton");

   GtkWidget *on_button = widget_from_builder("single_map_skeleton_on_radiobutton");
   GtkWidget *off_button = widget_from_builder("single_map_skeleton_off_radiobutton");

   if (imol >= 0) {
      if (molecules[imol].xskel_is_filled) {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button),  TRUE);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), FALSE);
      } else {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button),  FALSE);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), TRUE);
      }
   }
}

void
graphics_info_t::set_contour_sigma_button_and_entry(GtkWidget *window, int imol) {

   GtkWidget *entry = widget_from_builder("single_map_sigma_step_entry");
   GtkWidget *checkbutton = widget_from_builder("single_map_sigma_checkbutton");

   if (imol < n_molecules()) {
      if (molecules[imol].has_xmap()) {
	 float v = molecules[imol].contour_sigma_step;
	 gtk_editable_set_text(GTK_EDITABLE(entry), float_to_string(v).c_str());
	 if (molecules[imol].contour_by_sigma_flag) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
	 } else {
	    gtk_widget_set_sensitive(entry, FALSE);
	 }


	 GtkWidget *level_entry = widget_from_builder("single_map_properties_contour_level_entry");
	 float lev = molecules[imol].contour_level;
	 gtk_editable_set_text(GTK_EDITABLE(level_entry), float_to_string(lev).c_str());
      }
   }

}

// coot::rama_plot is an unknown type if we don't have canvas
#ifdef HAVE_GOOCANVAS // 20230501-PE old (Goocanvas) version of the Rama plot
void
graphics_info_t::handle_rama_plot_update(coot::rama_plot *plot) {

   if (plot) {
      // if it's a normal plot: update it
      if (plot->is_kleywegt_plot()) {
         // are the molecule numbers from which the kleywegt plot
         // was generated still valid?
         std::pair<int, int> p = plot->molecule_numbers();
         if (graphics_info_t::molecules[p.first].has_model() &&
             graphics_info_t::molecules[p.second].has_model()) {
            std::pair<std::string, std::string> chain_ids = plot->chain_ids();
            std::cout << "updating kleywegt plot with chain ids :" << chain_ids.first
                      << ": :" << chain_ids.second << ":" << std::endl;
            if (plot->kleywegt_plot_uses_chain_ids_p())
               plot->draw_it(p.first, p.second,
                             graphics_info_t::molecules[p.first].atom_sel.mol,
                     graphics_info_t::molecules[p.second].atom_sel.mol,
                     chain_ids.first, chain_ids.second);
            else
               plot->draw_it(p.first, p.second,
                             graphics_info_t::molecules[p.first].atom_sel.mol,
                     graphics_info_t::molecules[p.second].atom_sel.mol);
         } else {
            plot->hide_yourself();
         }
      } else {
         // check if selection is there
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(plot->selection_checkbutton))) {
            plot->apply_selection_from_widget();
         } else {
            plot->draw_it(molecules[imol_moving_atoms].atom_sel.mol);
         }
      }
   } else {
      std::cout << "ERROR:: (trapped) in handle_rama_plot_update() attempt to draw to null plot\n";
   }
}
#endif

// static
void
graphics_info_t::set_transient_for_main_window(GtkWidget *dialog) {

   GtkWidget *main_window_widget = graphics_info_t::get_main_window();
   if (main_window_widget) {
      GtkWindow *main_window = GTK_WINDOW(main_window_widget);
      gtk_window_set_transient_for(GTK_WINDOW(dialog), main_window);
   }
}


// --------------------------------------------------------------------------------
//                 residue info widget
// --------------------------------------------------------------------------------

#include "widget-from-builder.hh"
#include "ligand/primitive-chi-angles.hh"


// 23 Oct 2003: Why is this so difficult?  Because we want to attach
// atom info (what springs to mind is a pointer to the atom) for each
// entry, so that when the text in the entry is changed, we know to
// modify the atom.
//
// The problem with that is that behind our backs, that atom could
// disappear (e.g close molecule or delete residue, mutate or
// whatever), we are left with a valid looking (i.e. non-NULL)
// pointer, but the memory to which is points is invalid -> crash when
// we try to reference it.
//
// How shall we get round this?  refcounting?
//
// Instead, let's make a trivial class that contains the information
// we need to do a SelectAtoms to find the pointer to the atom, that
// class shall be called select_atom_info, it shall contain the
// molecule number, the chain id, the residue number, the insertion
// code, the atom name, the atom altconf.
//


void
graphics_info_t::output_residue_info_as_text(int atom_index, int imol) {

   // It would be cool to flash the residue here.
   // (heh - it is).
   // 20230111-PE gone now.
   //
   graphics_info_t g;
   mmdb::Atom *picked_atom = g.molecules[imol].atom_sel.atom_selection[atom_index];

   if (picked_atom) {

      // g.flash_selection(imol,
      //                   picked_atom->residue->seqNum,
      //                   picked_atom->GetInsCode(),
      //                   picked_atom->residue->seqNum,
      //                   picked_atom->GetInsCode(),
      //                   picked_atom->altLoc,
      //                   picked_atom->residue->GetChainID());

      mmdb::PAtom *atoms = NULL;
      int n_atoms = 0;
      mmdb::Residue *residue_p = picked_atom->residue;
      if (residue_p) {
         residue_p->GetAtomTable(atoms,n_atoms);
         if (atoms) {
            for (int i=0; i<n_atoms; i++) {
               mmdb::Atom *at = atoms[i];
               if (at) {
                  std::string segid = atoms[i]->segID;
                  std::cout << "(" << imol << ") \""
                            << at->name << "\"/"
                            << at->GetModelNum()
                            << "/\""
                            << at->GetChainID()  << "\"/"
                            << at->GetSeqNum()   << "/\""
                            << at->GetResName()
                            << "\", \""
                            << segid
                            << "\" occ: "
                            << at->occupancy
                            << " with B-factor: "
                            << at->tempFactor
                            << " element: \""
                            << at->element
                            << "\""
                            << " at " << "("
                            << at->x << ","
                            << at->y << ","
                            << at->z << ")" << std::endl;
               }
            }
         }
      }

      // chi angles:
      coot::primitive_chi_angles chi_angles(picked_atom->residue);
      try {
         std::vector<coot::alt_confed_chi_angles> chis = chi_angles.get_chi_angles();
         if (chis.size() > 0) {
            unsigned int i_chi_set = 0;
            std::cout << "   Chi Angles:" << std::endl;
            for (unsigned int ich=0; ich<chis[i_chi_set].chi_angles.size(); ich++) {
               std::cout << "     chi "<< chis[i_chi_set].chi_angles[ich].first << ": "
                         << chis[i_chi_set].chi_angles[ich].second
                         << " degrees" << std::endl;
            }
         } else {
            std::cout << "No Chi Angles for this residue" << std::endl;
         }
      }
      catch (const std::runtime_error &mess) {
         std::cout << mess.what() << std::endl;
      }

   }
}

//static
void
graphics_info_t::output_residue_info_dialog(int imol, const coot::residue_spec_t &rs) {

   // This is a kludge - really the main output_residue_info_dialog function should take a Residue *.
   // For now I will just find the atom index of the first atom in rs;

   graphics_info_t g;
   mmdb::Residue *residue_p = g.get_residue(imol, rs); // get_residue() is non-static.
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *rat = residue_atoms[iat];
         if (! rat->isTer()) {
            // now find at in the atom selection

            const auto &atom_sel = molecules[imol].atom_sel;
            for (int i=0; i<atom_sel.n_selected_atoms; i++) {
               if (atom_sel.atom_selection[i] == rat) {
                  output_residue_info_dialog(imol, i);
                  break;
               }
            }
            break;
         }
      }
   }

}



// The reader is graphics_info_t::apply_residue_info_changes();
//
void
graphics_info_t::output_residue_info_dialog(int imol, int atom_index) {

   auto add_chi_angles = [] (mmdb::Residue *residue, GtkWidget *dialog) {
                            try {
                               coot::primitive_chi_angles chi_angles(residue);
                               std::vector<coot::alt_confed_chi_angles> chis = chi_angles.get_chi_angles();
                               // GtkWidget *chi_angles_frame = lookup_widget(dialog, "chi_angles_frame");
                               GtkWidget *chi_angles_frame = widget_from_builder("residue_info_chi_angles_frame");
                               gtk_widget_set_visible(chi_angles_frame, TRUE);
                               if (chis.size() > 0) {
                                  unsigned int i_chi_set = 0;
                                  // undisplay all chi angles
                                  for (unsigned int i=1; i<=5; i++) {
                                     std::string label_name = "residue_info_chi_" + std::to_string(i) + "_label";
                                     GtkWidget *label = widget_from_builder(label_name.c_str());
                                     if (label) {
                                        gtk_widget_set_visible(label, FALSE);
                                     }
                                  }
                                  // now display
                                  for (unsigned int ich=0; ich<chis[i_chi_set].chi_angles.size(); ich++) {

                                     int ic = chis[i_chi_set].chi_angles[ich].first;
                                     std::string label_name = "residue_info_chi_";
                                     label_name += coot::util::int_to_string(ic);
                                     label_name += "_label";
                                     // GtkWidget *label = lookup_widget(dialog, label_name.c_str());
                                     GtkWidget *label = widget_from_builder(label_name.c_str());
                                     if (label) {
                                        std::string text = "Chi ";
                                        text += coot::util::int_to_string(ic);
                                        text += ":  ";
                                        if (chis[i_chi_set].alt_conf != "") {
                                           text += " alt conf: ";
                                           text += chis[i_chi_set].alt_conf;
                                           text += " ";
                                        }
                                        text += coot::util::float_to_string(chis[i_chi_set].chi_angles[ich].second);
                                        text += " degrees";
                                        gtk_label_set_text(GTK_LABEL(label), text.c_str());
                                        gtk_widget_set_visible(label, TRUE);
                                     } else {
                                        std::cout << "WARNING:: chi label not found " << label_name << std::endl;
                                     }
                                  }
                               }
                            }
                            catch (const std::runtime_error &mess) {
                               std::cout << mess.what() << std::endl;
                            }
                         };

   if (graphics_info_t::residue_info_edits.size() > 0) {

      std::string s =  "You have pending (un-Applied) residue edits.\n";
      s += "Deal with them first.";
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_set_visible(w, TRUE);

   } else {

      if (imol <graphics_info_t::n_molecules()) {
         if (graphics_info_t::molecules[imol].has_model()) {
            if (atom_index < graphics_info_t::molecules[imol].atom_sel.n_selected_atoms) {

               graphics_info_t g;
               output_residue_info_as_text(atom_index, imol);
               mmdb::Atom *selected_atom = g.molecules[imol].atom_sel.atom_selection[atom_index];
               std::string residue_name = selected_atom->GetResName();
               mmdb::PPAtom atoms;
               int n_atoms;
               selected_atom->residue->GetAtomTable(atoms,n_atoms);
               // GtkWidget *dialog = wrapped_create_residue_info_dialog(); // just the (unfilled) dialog
               GtkWidget *dialog = widget_from_builder("residue_info_dialog");

               mmdb::Residue *residue = selected_atom->residue;
               coot::residue_spec_t *res_spec_p = new coot::residue_spec_t(residue->GetChainID(), residue->GetSeqNum(), residue->GetInsCode());

               // fill the master atom
               // GtkWidget *master_occ_entry      = lookup_widget(widget, "residue_info_master_atom_occ_entry");
               // GtkWidget *master_b_factor_entry = lookup_widget(widget, "residue_info_master_atom_b_factor_entry");

               GtkWidget *master_occ_entry      = widget_from_builder("residue_info_master_atom_occ_entry");
               GtkWidget *master_b_factor_entry = widget_from_builder("residue_info_master_atom_b_factor_entry");

               // Do I need to clear the signals from the previous time that this widget was shown?
               // Or how do I add these just once?

               gtk_editable_set_text(GTK_EDITABLE(master_occ_entry), "1.0");

               std::string b_entry_text = graphics_info_t::float_to_string(graphics_info_t::default_new_atoms_b_factor);
               gtk_editable_set_text(GTK_EDITABLE(master_b_factor_entry), b_entry_text.c_str());

               g_object_set_data(G_OBJECT(dialog), "res_spec_p",  res_spec_p);
               g.fill_output_residue_info_widget(dialog, imol, residue_name, atoms, n_atoms);
               set_transient_for_main_window(dialog);
               gtk_widget_set_visible(dialog, TRUE);
               g.reset_residue_info_edits();

               add_chi_angles(residue, dialog);

            }
         }
      }
   }

}



void
graphics_info_t::fill_output_residue_info_widget(GtkWidget *dialog, int imol,
						 const std::string residue_name,
						 mmdb::PPAtom atoms, int n_atoms) {

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

   // first do the label of the dialog
   // GtkWidget *label_widget = lookup_widget(widget, "residue_info_residue_label");
   // GtkWidget *residue_name_widget = lookup_widget(widget, "residue_info_residue_name_label");

   GtkWidget *label_widget        = widget_from_builder("residue_info_residue_label");
   GtkWidget *residue_name_widget = widget_from_builder("residue_info_residue_name_label");

   // GtkWidget *table = lookup_widget(widget, "residue_info_atom_table");
   GtkWidget *grid = widget_from_builder("residue_info_atom_grid");
   g_object_set_data(G_OBJECT(grid), "imol", GINT_TO_POINTER(imol));
   clear_the_grid(grid);

   std::cout << "::::::::::::::::: fill_output_residue_info_widget() grid " << grid << std::endl;

   // set the column labels of the grid
   gint top_attach = 0;
   GtkWidget *atom_info_label = gtk_label_new(" Atom Info ");
   GtkWidget *occupancy_label = gtk_label_new(" Occupancy ");
   GtkWidget  *b_factor_label = gtk_label_new(" Temperature Factor ");
   GtkWidget  *alt_conf_label = gtk_label_new(" Alt Conf ");
   gtk_grid_attach(GTK_GRID(grid), atom_info_label, 0, top_attach, 1, 1);
   gtk_grid_attach(GTK_GRID(grid), occupancy_label, 1, top_attach, 1, 1);
   gtk_grid_attach(GTK_GRID(grid),  b_factor_label, 2, top_attach, 1, 1);
   gtk_grid_attach(GTK_GRID(grid),  alt_conf_label, 4, top_attach, 1, 1);
   gtk_widget_set_margin_bottom(atom_info_label, 8);
   gtk_widget_set_margin_bottom(occupancy_label, 8);
   gtk_widget_set_margin_bottom( b_factor_label, 8);

   // name
   graphics_info_t g;
   std::string res_name_string = residue_name + std::string(" ");
   std::pair<bool, std::string> p = g.Geom_p()->get_monomer_name(residue_name, imol);
   if (p.first) {
      res_name_string += p.second;
      gtk_label_set_text(GTK_LABEL(residue_name_widget), res_name_string.c_str());
   }

   gtk_widget_set_size_request(dialog, -1, 600); // so that we can see a few more rows

   residue_info_n_atoms = n_atoms;
   for (int i=0; i<n_atoms; i++)
      fill_output_residue_info_widget_atom(dialog, grid, imol, atoms[i], i);

   if (n_atoms > 0) {

      std::string chain_id = atoms[0]->GetChainID();
      int res_no = atoms[0]->GetSeqNum();
      std::string ins_code = atoms[0]->residue->GetInsCode();

      std::string label = "Molecule: ";
      label += int_to_string(imol);
      label += " ";
      label += molecules[imol].name_;
      label += "\n";
      label += chain_id;
      label += " ";
      label += std::to_string(res_no);
      label += " ";
      label += ins_code;

      gtk_label_set_text(GTK_LABEL(label_widget), label.c_str());

   }
   
}

void
graphics_info_t::fill_output_residue_info_widget_atom(GtkWidget *dialog, GtkWidget *grid, int imol, mmdb::PAtom atom,
						      int iatom) {

   // We pass the dialog so that the residue occ and b-factor entries can be
   // added via g_object_set_data().

   // GtkWidget *residue_info_dialog_local = lookup_widget(table, "residue_info_dialog");
   GtkWidget *residue_info_dialog_local = widget_from_builder("residue_info_dialog");

   gint left_attach = 0;
   gint top_attach = iatom+1; // because we have to top row begin the column labels

   // The text label of the atom name:
   left_attach = 0;
   std::string label_str = "  ";
   label_str += atom->GetChainID();
   label_str += "/";
   label_str += graphics_info_t::int_to_string(atom->GetSeqNum());
   if (std::string(atom->GetInsCode()) != "") {
      label_str += atom->GetInsCode();
   }
   label_str += " ";
   label_str += atom->GetResName();
   label_str += "/";
   label_str += atom->name;
   if (false) { // do alt locs in the atom label (20090914, not now that we have altloc entries)
      if (std::string(atom->altLoc) != std::string("")) {
	 label_str += ",";
	 label_str += atom->altLoc;
      }
   }
   if (std::string(atom->segID) != std::string("")) {
      label_str += " ";
      label_str += atom->segID;
   }
   label_str += "  ";

   GtkWidget *residue_info_atom_info_label = gtk_label_new (label_str.c_str());
   // gtk_table_attach(GTK_TABLE(table), residue_info_atom_info_label,
   //      	    left_attach, right_attach, top_attach, bottom_attach,
   //      	    xopt, yopt, xpad, ypad);

   gtk_grid_attach(GTK_GRID(grid), residue_info_atom_info_label, left_attach, top_attach, 1, 1);
   // gtk_widget_ref (residue_info_atom_info_label);
   g_object_set_data_full(G_OBJECT (residue_info_dialog_local), "residue_info_atom_info_label",
                          residue_info_atom_info_label, NULL);
   gtk_widget_set_visible (residue_info_atom_info_label, TRUE);

   // The Occupancy entry:
   left_attach = 1;
   coot::select_atom_info *ai = new coot::select_atom_info;
   *ai = coot::select_atom_info(iatom, imol,
				std::string(atom->GetChainID()),
				atom->GetSeqNum(),
				std::string(atom->GetInsCode()),
				std::string(atom->name),
				std::string(atom->altLoc));

   std::string widget_name = "residue_info_occ_entry_";
   widget_name += int_to_string(iatom);
   //
   GtkWidget *residue_info_occ_entry = gtk_entry_new();
   // gtk_widget_ref (residue_info_occ_entry);
   
   // g_object_set_data_full(G_OBJECT (residue_info_dialog_local),
   // widget_name.c_str(), residue_info_occ_entry, NULL);

   g_object_set_data(G_OBJECT(dialog), widget_name.c_str(), residue_info_occ_entry);

   // gtk_widget_set_size_request(residue_info_occ_entry, 20, -1);

   gtk_editable_set_width_chars(GTK_EDITABLE(residue_info_occ_entry), 6);
   gtk_widget_set_visible(residue_info_occ_entry, TRUE);
   g_object_set_data(G_OBJECT(residue_info_occ_entry), "select_atom_info", ai);
   gtk_editable_set_text(GTK_EDITABLE(residue_info_occ_entry),
                         graphics_info_t::float_to_string(atom->occupancy).c_str());
   // gtk_table_attach(GTK_TABLE(table), residue_info_occ_entry,
   //      	    left_attach, right_attach, top_attach, bottom_attach,
   //      	    xopt, yopt, xpad, ypad);

   gtk_grid_attach(GTK_GRID(grid), residue_info_occ_entry, left_attach, top_attach, 1, 1);

      // Note that we have to use key_release_event because if we use
   // key_press_event, when we try to get the value from the widget
   // (gtk_entry_get_text) then that does not see the key/number that
   // was just pressed.

   // B-factor entry:
   left_attach = 2;

   widget_name = "residue_info_b_factor_entry_";
   widget_name += int_to_string(iatom);

   GtkWidget *residue_info_b_factor_entry = gtk_entry_new ();
   // gtk_widget_ref (residue_info_b_factor_entry);
   // g_object_set_data_full(G_OBJECT (residue_info_dialog_local),
   // widget_name.c_str(), residue_info_b_factor_entry, NULL);

   g_object_set_data(G_OBJECT(dialog), widget_name.c_str(), residue_info_b_factor_entry);

   // gtk_widget_set_size_request(residue_info_b_factor_entry, 20, -1);

   gtk_editable_set_width_chars(GTK_EDITABLE(residue_info_b_factor_entry), 6);

   gtk_widget_set_visible (residue_info_b_factor_entry, TRUE);
   gtk_editable_set_text(GTK_EDITABLE(residue_info_b_factor_entry),
                         graphics_info_t::float_to_string(atom->tempFactor).c_str());
   g_object_set_data(G_OBJECT(residue_info_b_factor_entry), "select_atom_info", ai);

   gtk_grid_attach(GTK_GRID(grid), residue_info_b_factor_entry, left_attach, top_attach, 1, 1);

   // Alt Conf label:
   GtkWidget *alt_conf_label = gtk_label_new(" ");
   gtk_widget_set_visible(alt_conf_label, TRUE);
   left_attach = 3;
   // gtk_table_attach(GTK_TABLE(table), alt_conf_label,
   //      	    left_attach, right_attach, top_attach, bottom_attach,
   //      	    xopt, yopt, xpad, ypad);
   gtk_grid_attach(GTK_GRID(grid), alt_conf_label, left_attach, top_attach, 1, 1);


   // The Alt Conf entry:
   left_attach = 4;
   ai = new coot::select_atom_info;
   *ai = coot::select_atom_info(iatom, imol,
				std::string(atom->GetChainID()),
				atom->GetSeqNum(),
				std::string(atom->GetInsCode()),
				std::string(atom->name),
				std::string(atom->altLoc));
   widget_name = "residue_info_altloc_entry_";
   widget_name += int_to_string(iatom);
   //
   GtkWidget *residue_info_altloc_entry = gtk_entry_new ();
   // gtk_widget_ref (residue_info_altloc_entry);
   g_object_set_data_full(G_OBJECT (residue_info_dialog_local),
			  widget_name.c_str(), residue_info_altloc_entry,
			  NULL);

   // gtk_widget_set_size_request(residue_info_altloc_entry, 20, -1);
   gtk_editable_set_width_chars(GTK_EDITABLE(residue_info_altloc_entry), 6);

   gtk_widget_set_visible (residue_info_altloc_entry, TRUE);
   g_object_set_data(G_OBJECT(residue_info_altloc_entry), "select_atom_info", ai);
   gtk_editable_set_text(GTK_EDITABLE(residue_info_altloc_entry), atom->altLoc);
   // gtk_table_attach(GTK_TABLE(table), residue_info_altloc_entry,
   //      	    left_attach, right_attach, top_attach, bottom_attach,
   //      	    xopt, yopt, xpad, ypad);
   gtk_grid_attach(GTK_GRID(grid), residue_info_altloc_entry, left_attach, top_attach, 1, 1);

}


//static
void
graphics_info_t::residue_info_edit_b_factor_apply_to_other_entries_maybe(GtkWidget *master_b_factor_entry) {

   GtkWidget *checkbutton = widget_from_builder("residue_info_b_factor_apply_all_checkbutton");
   GtkWidget *grid        = widget_from_builder("residue_info_atom_grid");

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) {
      const char *txt = gtk_editable_get_text(GTK_EDITABLE(master_b_factor_entry));
      // the first line is labels
      for (int iat=1; iat<10000; iat++) {
         GtkWidget *w = gtk_grid_get_child_at(GTK_GRID(grid), 2, iat);
         if (!w) break;
         gtk_editable_set_text(GTK_EDITABLE(w), txt);
      }
   }
}


//static
void
graphics_info_t::residue_info_edit_occ_apply_to_other_entries_maybe(GtkWidget *master_occ_entry) {

   GtkWidget *occ_checkbutton = widget_from_builder("residue_info_occ_apply_all_checkbutton");
   GtkWidget *alt_checkbutton = widget_from_builder("residue_info_occ_apply_to_altconf_checkbutton");
   GtkWidget *alt_entry       = widget_from_builder("residue_info_occ_apply_to_alt_conf_entry");
   GtkWidget *grid            = widget_from_builder("residue_info_atom_grid");

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(occ_checkbutton))) {
      const char *txt = gtk_editable_get_text(GTK_EDITABLE(master_occ_entry));
      // the first line is labels
      for (int iat=1; iat<10000; iat++) {
         GtkWidget *w = gtk_grid_get_child_at(GTK_GRID(grid), 1, iat);
         if (!w) {
            // std::cout << "null editable at iat " << iat << std::endl;
            break;
         } else {
            gtk_editable_set_text(GTK_EDITABLE(w), txt);
         }
      }
   }
}


// static
void
graphics_info_t::residue_info_add_b_factor_edit(coot::select_atom_info sai, float val) {

   graphics_info_t g;
   short int made_substitution_flag = 0;
   for (unsigned int i=0; i<g.residue_info_edits.size(); i++) {
      if (sai.udd == g.residue_info_edits[i].udd) {
	 g.residue_info_edits[i].add_b_factor_edit(val);
	 made_substitution_flag = 1;
	 break;
      }
   }
   if (! made_substitution_flag) {
      sai.add_b_factor_edit(val);
      g.residue_info_edits.push_back(sai);
   }
}

// static
void
graphics_info_t::residue_info_add_occ_edit(coot::select_atom_info sai, float val) {

   graphics_info_t g;
   short int made_substitution_flag = 0;
   for (unsigned int i=0; i<g.residue_info_edits.size(); i++) {
      if (sai.udd == g.residue_info_edits[i].udd) {
	 g.residue_info_edits[i].add_occ_edit(val);
	 made_substitution_flag = 1;
	 break;
      }
   }
   if (! made_substitution_flag) {
      sai.add_occ_edit(val);
      g.residue_info_edits.push_back(sai);
   }
}

void
graphics_info_t::apply_residue_info_changes() {

   std::cout << "New apply_residue_info_changes(): n_atoms: " << residue_info_n_atoms << std::endl;

   GtkWidget *grid = widget_from_builder("residue_info_atom_grid");

   reset_residue_info_edits(); // clear them

   if (residue_info_n_atoms > 0) {
      // the first line in the grid is labels
      for (int iat=1; iat <= residue_info_n_atoms; iat++) {
         GtkWidget *occ_entry      = gtk_grid_get_child_at(GTK_GRID(grid), 1, iat);
         GtkWidget *b_factor_entry = gtk_grid_get_child_at(GTK_GRID(grid), 2, iat);
         GtkWidget *alt_conf_entry = gtk_grid_get_child_at(GTK_GRID(grid), 4, iat); //    in above: left_attach = 4;
         if (GTK_IS_EDITABLE(occ_entry)) {
            if (GTK_IS_EDITABLE(b_factor_entry)) {
               const gchar *t_occ  = gtk_editable_get_text(GTK_EDITABLE(occ_entry));
               const gchar *t_bfac = gtk_editable_get_text(GTK_EDITABLE(b_factor_entry));
               const gchar *t_altconfc = gtk_editable_get_text(GTK_EDITABLE(alt_conf_entry));
               try {
                  float f_occ  = coot::util::string_to_float(std::string(t_occ));
                  float f_bfac = coot::util::string_to_float(std::string(t_bfac));
                  std::string t_alt_conf;
                  if (t_altconfc)
                     t_alt_conf = t_altconfc;
                  coot::select_atom_info *ai_p = static_cast<coot::select_atom_info *>(g_object_get_data(G_OBJECT(occ_entry), "select_atom_info"));
                  if (ai_p) {
                     coot::select_atom_info &ai(*ai_p);
                     ai.add_b_factor_edit(f_bfac);
                     ai.add_occ_edit(f_occ);
                     std::cout << "add_altloc_edit() " << t_alt_conf << std::endl;
                     ai.add_altloc_edit(t_alt_conf);
                     residue_info_edits.push_back(ai);
                  }
               }
               catch (const std::runtime_error &e) {
                  std::cout << "WARNING::" << e.what() << std::endl;
               }
            }
         }
      }
   }

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(grid), "imol"));
   if (is_valid_model_molecule(imol))
      molecules[imol].apply_atom_edits(residue_info_edits);

   reset_residue_info_edits();
}


// static
// (should be called by the destroy event and the close button)
void
graphics_info_t::residue_info_release_memory(GtkWidget *dialog) {

   return; // 20220401-PE let's not delete anything for now

   GtkWidget *entry;
   for (int i=0; i<residue_info_n_atoms; i++) {
      std::string widget_name = "residue_info_b_factor_entry_";
      widget_name += int_to_string(i);
      // entry = lookup_widget(dialog, widget_name.c_str());
      entry = nullptr;
      std::cout << "FIXME:: in residue_info_release_memory() look up entry correctly" << std::endl;
      if (entry) {
	 coot::select_atom_info *sai_p =
	    (coot::select_atom_info *) g_object_get_data(G_OBJECT(entry), "select_atom_info");
	 if (sai_p) {
	    // delete sai_p; // memory bug.  We cant do this
	    // std::cout << "hmmm.. not deleting sai_p" << std::endl;
	 } else {
	    std::cout << "ERROR:: no user data in b-factor entry widget\n";
	 }
      } else {
         std::cout << "ERROR:: in residue_info_release_memory() failed to find entry " << std::endl;
      }
      // same for occ entry:
      widget_name = "residue_info_occ_entry_";
      widget_name += int_to_string(i);
      // entry = lookup_widget(dialog, widget_name.c_str());
      entry = 0;
      std::cout << "FIXME in residue_info_release_memory() look up entry correctly" << std::endl;
      if (entry) {
	 coot::select_atom_info *sai_p =
	    (coot::select_atom_info *) g_object_get_data(G_OBJECT(entry), "select_atom_info");
	 if (sai_p) {
	    // std::cout << "not deleting sai_p" << std::endl;
	    // delete sai_p; // memory bug.  We cant do this
	 } else {
	    std::cout << "ERROR:: no user data in occ entry widget\n";
	 }
      } else {
         std::cout << "ERROR:: in residue_info_release_memory() failed to find entry 2 " << std::endl;
      }
   }
}

// static
void
graphics_info_t::pointer_atom_molecule_combobox_changed(GtkWidget *combobox, gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
   std::cout << "debug:: changed to imol " << imol << std::endl;
   g.user_pointer_atom_molecule = imol;

}


void
graphics_info_t::new_fill_combobox_with_coordinates_options(GtkWidget *combobox_molecule,
                                                            GCallback callback_func,
                                                            int imol_active) {

   // This presumes that the combobox_molecule is fresh and packed into a widget already

   auto get_molecule_indices = [] () {
                                  std::vector<int> molecule_indices;
                                  for (int i=0; i<graphics_info_t::n_molecules(); i++) {
                                     if (graphics_info_t::molecules[i].has_model()) {
                                        molecule_indices.push_back(i);
                                     }
                                  }
                                  return molecule_indices;
                               };

   std::vector<int> molecule_indices = get_molecule_indices();

   GtkTreeModel *model_1 = gtk_combo_box_get_model(GTK_COMBO_BOX(combobox_molecule));
   std::cout << "debug:: new_fill_combobox_with_coordinates_options() model_1 " << model_1 << std::endl;
   GtkListStore *list_store = GTK_LIST_STORE(model_1);
   std::cout << "debug:: new_fill_combobox_with_coordinates_options() list_store " << list_store << std::endl;
   // gtk_list_store_clear(list_store);

   GtkListStore *store = gtk_list_store_new(2, G_TYPE_INT, G_TYPE_STRING);

   std::cout << "debug:: new_fill_combobox_with_coordinates_options() list_store " << store << std::endl;

   gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox_molecule));

   GtkTreeIter iter;
   for (unsigned int ii=0; ii<molecule_indices.size(); ii++) {
      const auto &imol = molecule_indices[ii];
      const molecule_class_info_t &m = graphics_info_t::molecules[imol];
      std::string ss = std::to_string(imol) + " " + m.name_for_display_manager();
      gtk_list_store_append(store, &iter);
      gtk_list_store_set(store, &iter, 0, imol, 1, ss.c_str(), -1);
   }
   
   GtkTreeModel *model = GTK_TREE_MODEL(store);
   GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
   gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combobox_molecule), renderer, true);
   gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combobox_molecule), renderer, "text", 1, NULL);
   gtk_combo_box_set_model(GTK_COMBO_BOX(combobox_molecule), model);

   for (unsigned int ii=0; ii<molecule_indices.size(); ii++) {
      const auto &imol = molecule_indices[ii];
      const molecule_class_info_t &m = graphics_info_t::molecules[imol];
      if (imol == imol_active) {
         // 20220415-PE this doesn't work (for renumber residues - annoying)

         std::cout << "!!!!!!!!!!! setting active on a gtk combobox " << imol_active << std::endl;
         gtk_combo_box_set_active(GTK_COMBO_BOX(combobox_molecule), imol_active);
         std::cout << "!!!!!!!!!!! combobox get_active() returns " <<  gtk_combo_box_get_active(GTK_COMBO_BOX(combobox_molecule)) << std::endl;
         if (GTK_IS_COMBO_BOX(combobox_molecule))
            std::cout << "!!!!!!!!!!! " << "combobox is a combobox" << std::endl;
         if (GTK_IS_COMBO_BOX_TEXT(combobox_molecule))
            std::cout << "!!!!!!!!!!! " << "combobox is a comboboxtext" << std::endl;
      }
   }


   if (callback_func)
      g_signal_connect(combobox_molecule, "changed", callback_func, NULL);

}

void
graphics_info_t::fill_combobox_with_molecule_options(GtkWidget *combobox,
                                                     GCallback signal_func,
                                                     int imol_active_position,
                                                     const std::vector<int> &molecules_index_vec) {

   // This might be OK if combobox had been created by us, not glade code, and has been added to an hbox already.

   gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox));

   GtkListStore *store = gtk_list_store_new(2, G_TYPE_INT, G_TYPE_STRING);
   GtkTreeIter iter;
   int active_idx = 0;
   unsigned int n_mol = molecules_index_vec.size();
   for (unsigned int imap=0; imap<n_mol; imap++) {

      //std::cout << "::::::::: debug in fill_combobox_with_molecule_options() imap " << imap << std::endl;
      int imol = molecules_index_vec[imap];
      std::string ss; // = int_to_string(imol); done in renderer now.
      ss = int_to_string(imol); // 20220321-PE well, the map number is missing in the export map molecule chooser
                                // combobox - add it back and sort out duplication later (maybe model molecules?)
      ss += " " ;
      int ilen = molecules[imol].name_.length();
      int left_size = ilen-go_to_atom_menu_label_n_chars_max;
      if (left_size <= 0)
	 left_size = 0;
      else
	 ss += "...";
      ss += molecules[imol].name_.substr(left_size, ilen);

      gtk_list_store_append(store, &iter);
      gtk_list_store_set(store, &iter, 0, imol, 1, ss.c_str(), -1);

      if (imol == imol_active_position)
	 active_idx = imap;
   }

   if (signal_func)
      g_signal_connect(combobox, "changed", signal_func, NULL);
   GtkTreeModel *model = GTK_TREE_MODEL(store);
   GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
   gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combobox), renderer, TRUE);
   gtk_cell_layout_set_attributes (GTK_CELL_LAYOUT(combobox), renderer, "text", 1, NULL);
   gtk_combo_box_set_model(GTK_COMBO_BOX(combobox), model);

   // maybe this can go into the above loop?
   if (molecules_index_vec.size() > 0)
      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), active_idx);

}

void
graphics_info_t::fill_combobox_with_coordinates_options(GtkWidget *combobox,
							GCallback callback_func,
							int imol_active) {

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- don't use this function -----\n");
   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- use fill_combobox_with_molecule_options --\n");

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- start -----\n");

   std::vector<int> fill_with_these_molecules;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model()) {
         fill_with_these_molecules.push_back(imol);
      }
   }

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- Here A -----\n");

   if (false)
      std::cout << "debug:: --- in fill_combobox_with_coordinates_options() n_molecules: "
                << fill_with_these_molecules.size() << " and imol_active " << imol_active
                << std::endl;

   GtkListStore *store = gtk_list_store_new(2, G_TYPE_INT, G_TYPE_STRING);
   GtkTreeIter iter;
   int active_idx = 0;
   int n_mol = fill_with_these_molecules.size();

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- Here B -----\n");

   for (int idx=0; idx<n_mol; idx++) {
      int imol = fill_with_these_molecules[idx];
      std::string ss; // = int_to_string(imol); done in renderer now.
      ss += " " ;
      int ilen = molecules[imol].name_.length();
      int left_size = ilen-go_to_atom_menu_label_n_chars_max;
      if (left_size <= 0)
         left_size = 0;
      else
         ss += "...";
      ss += molecules[imol].name_.substr(left_size, ilen);

      //  std::cout << "debug:: --- in fill_combobox_with_coordinates_options() "
      //            << imol << " " << ss << std::endl;
      gtk_list_store_append(store, &iter);
      gtk_list_store_set(store, &iter, 0, imol, 1, ss.c_str(), -1);

      if (imol == imol_active)
         active_idx = idx;

   }

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- Here C -----\n");

   if (callback_func)
      g_signal_connect(combobox, "changed", callback_func, NULL);
   GtkTreeModel *model = GTK_TREE_MODEL(store);
   GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
   gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combobox), renderer, TRUE);
   gtk_cell_layout_set_attributes (GTK_CELL_LAYOUT(combobox), renderer, "text", 1, NULL);
   gtk_combo_box_set_model(GTK_COMBO_BOX(combobox), model);

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- Here D with combobox %p\n", combobox);

   // maybe this can go into the above loop?
   if (! fill_with_these_molecules.empty())
      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), active_idx);

   printf("DEBUG:: fill_combobox_with_coordinates_options(): -------------------------- end -----\n");
}



// The callback_func pass here is connected to the combobox, not the menu items.
//
void
graphics_info_t::fill_combobox_with_coordinates_options_with_set_last(GtkWidget *combobox,
								      GCallback callback_func,
								      bool set_last_active_flag) {

   int imol_active = -1;
   std::vector<int> fill_with_these_molecules;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model()) {
         fill_with_these_molecules.push_back(imol);
      }
   }
   if (fill_with_these_molecules.size() > 0) {
      imol_active = fill_with_these_molecules[0];

      if (set_last_active_flag)
       	 imol_active = fill_with_these_molecules.back();
   }

   fill_combobox_with_coordinates_options(combobox, callback_func, imol_active);

}

#if 0
// 20100629 this was used in fill_renumber_residue_range_dialog() - the modern way, I guess.
void
graphics_info_t::fill_option_menu_with_coordinates_options_internal_2(GtkWidget *option_menu,
								      GtkSignalFunc callback_func,
								      short int set_last_active_flag,
								      int imol_active) {

   std::vector<int> fill_with_these_molecules;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model()) {
	 fill_with_these_molecules.push_back(imol);
      }
   }
   fill_option_menu_with_coordinates_options_internal_3(option_menu, callback_func, fill_with_these_molecules,
							set_last_active_flag, imol_active);
}
#endif

#if 0
void
graphics_info_t::fill_option_menu_with_coordinates_options_internal_3(GtkWidget *option_menu,
								      GtkSignalFunc callback_func,
								      std::vector<int> fill_with_these_molecules,
								      short int set_last_active_flag,
								      int imol_active) {

   // like the column labels from an mtz file, similarly fill this
   // option_menu with items that correspond to molecules that have
   // coordinates.
   //

   // Get the menu of the optionmenu (which was set in interface.c:
   // gtk_option_menu_set_menu (GTK_OPTION_MENU (go_to_atom_molecule_optionmenu),
   //                           go_to_atom_molecule_optionmenu_menu);
   //
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));


   // for the strangeness of destroying and re-adding the menu to the
   // option menu, set the comments in the
   // fill_close_option_menu_with_all_molecule_options function

   // menu is not GTK_MENU on Gtk2 Ubuntu kalypso 64 bit
   if (menu)
      gtk_widget_destroy(menu);

   /* Create a menu for the optionmenu button.  The various molecule
    numbers will be added to this menu as menuitems*/
   menu = gtk_menu_new();

   // GtkWidget *optionmenu_menu = gtk_menu_new();
   GtkWidget *menuitem;
   // int last_imol = 0;
   int last_menu_item_index = 0;

   if (false) {  // debug
      std::cout << "fill_option_menu_with_coordinates_options_internal_3 with these: "
		<< std::endl;
      for (unsigned int idx=0; idx<fill_with_these_molecules.size(); idx++) {
	 std::cout << fill_with_these_molecules[idx] << " ";
      }
      std::cout << std::endl;
   }

   int menu_index = 0; // for setting of imol_active as active mol in go to atom
   for (unsigned int idx=0; idx<fill_with_these_molecules.size(); idx++) {

      int jmol = fill_with_these_molecules[idx];

      if (molecules[jmol].has_model()) {

	 std::string ss = int_to_string(jmol);
	 ss += " " ;
	 int ilen = molecules[jmol].name_.length();
	 int left_size = ilen-go_to_atom_menu_label_n_chars_max;
	 if (left_size <= 0) {
	    // no chop
	    left_size = 0;
	 } else {
	    // chop
	    ss += "...";
	 }
	 ss += molecules[jmol].name_.substr(left_size, ilen);
	 menuitem = gtk_menu_item_new_with_label (ss.c_str());

	 // std::cout << "user pointer to int on callback set to " << jmol << std::endl;

	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     // GTK_SIGNAL_FUNC(go_to_atom_mol_button_select),
			     callback_func,
			     GINT_TO_POINTER(jmol));

	 // Note that we probably don't need to do the following
	 // because we already pass a GINT_TO_POINTER(jmol) in the
	 // signal connect.
	 //
	 // But on reflection.. perhaps we do because we do a
	 // menu_get_active in save_go_to_atom_mol_menu_active_position
	 //
	 // we set user data on the menu item, so that when this goto
	 // Atom widget is cancelled, we can whatever was the molecule
	 // number corresponding to the active position of the menu
	 //
	 // Should be freed in on_go_to_atom_cancel_button_clicked
	 // (callbacks.c)
	 //

	 gtk_object_set_user_data(GTK_OBJECT(menuitem), GINT_TO_POINTER(jmol));
	 gtk_menu_append(GTK_MENU(menu), menuitem);

	 // std::cout << " comparing " << jmol << " and imol_active " << imol_active << std::endl;
	 if (jmol == imol_active)
	    gtk_menu_set_active(GTK_MENU(menu), menu_index);

	 // we do need this bit of course:
	 gtk_widget_set_visible(menuitem, TRUE);
	 last_menu_item_index++;
	 menu_index++;
      }
   }

   // set any previously saved active position:
   // but overridden by flag:
   if (set_last_active_flag) {
      // explanation of -1 offset: when there are 2 menu items,
      // last_menu_item_index is 2, but the item index of the last
      // item is 2 - 1.
      //
      gtk_menu_set_active(GTK_MENU(menu), (last_menu_item_index-1));
   } else {
      // the old way (ie. not ..._with_active_mol() mechanism)
      if (imol_active == -1) {
	 gtk_menu_set_active(GTK_MENU(menu), (last_menu_item_index-1));
      }
   }

   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu),
			    menu);
}
#endif

#if 0
void
graphics_info_t::fill_option_menu_with_coordinates_options_internal_with_active_mol(GtkWidget *option_menu,
										    GtkSignalFunc callback_func,
										    int imol_active) {

   short int set_last_active_flag = 0;
   fill_option_menu_with_coordinates_options_internal_2(option_menu, callback_func,
							set_last_active_flag, imol_active);
}
#endif

#if 0
// not const
void
graphics_info_t::fill_option_menu_with_undo_options(GtkWidget *option_menu) {


   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   GtkWidget *menuitem;
   int active_mol_no = -1;
   int undo_local = undo_molecule;  // defaults to -1; If it is unset (-1)
			       // then pick the first undoable
			       // molecule as the undo_molecule.
   int pos_count = 0;

   for (int i=0; i<n_molecules(); i++) {
      // if (molecules[i].has_model()) {
      if (molecules[i].atom_sel.mol) {
	 if (molecules[i].Have_modifications_p()) {
	    if (undo_local == -1)
	       undo_local = i;
	    char s[200];
	    snprintf(s, 199, "%d", i);
	    std::string ss(s);
	    ss += " ";
	    ss += molecules[i].name_;
	    menuitem = gtk_menu_item_new_with_label(ss.c_str());
	    gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			       GTK_SIGNAL_FUNC(graphics_info_t::undo_molecule_select),
			       GINT_TO_POINTER(i));
	    gtk_menu_append(GTK_MENU(menu), menuitem);
	    gtk_widget_set_visible(menuitem, TRUE);
	    if (i == undo_molecule)
	       gtk_menu_set_active(GTK_MENU(menu), pos_count);
	    pos_count++;
	 }
      }
   }

   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

   // undo_local should have been set (this function should not be called
   // if there are not more than 1 molecule with modifications).
   if (undo_local > -1) {
      set_undo_molecule_number(undo_local);
   }
}
#endif

void
graphics_info_t::fill_combobox_with_undo_options(GtkWidget *combobox_molecule) {

   // make the first undo molecule (a molecule with changes) be the active one.

   // 20220708-PE this is how to clear a combobox

   if (combobox_molecule) {
      gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox_molecule));

      int imol_active = -1;
      for (int i=0; i<n_molecules(); i++) {
         if (molecules[i].has_model()) {
            if (molecules[i].atom_sel.mol) {
               if (molecules[i].Have_modifications_p()) {
                  imol_active = i;
                  break;
               }
            }
         }
      }

      GCallback callback = G_CALLBACK(undo_molecule_combobox_changed);
      fill_combobox_with_coordinates_options(combobox_molecule, callback, imol_active);
   } else {
      std::cout << "ERROR:: in fill_combobox_with_undo_options() combobox_molecule is null" << std::endl;
   }
}



// static
void
graphics_info_t::refinement_map_combobox_changed(GtkWidget *c, gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(c));
   // now what?
}



// a static function
void
graphics_info_t::undo_molecule_combobox_changed(GtkWidget *combobox, gpointer data) {
   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
   g.set_undo_molecule_number(imol);
   std::cout << "INFO:: undo molecule number set to " << imol << std::endl;
}

void
graphics_info_t::set_baton_build_params(int istart_resno,
					const char *chain_id,
					const char *backwards) {

   baton_build_params_active = 1; // don't ignore baton_build_params
				  // in placing atom.
   baton_build_start_resno = istart_resno;
   std::string dir(backwards);
   if (dir == "backwards") {
      baton_build_direction_flag = -1;
   } else {
      if (dir == "forwards") {
	 baton_build_direction_flag = 1;
      } else {
	 baton_build_direction_flag = 0; // unknown.
      }
   }
   baton_build_chain_id = std::string(chain_id);
   // std::cout << "DEBUG:: baton_build_chain_id set to " << baton_build_chain_id << std::endl;
}

void
graphics_info_t::model_fit_refine_unactive_togglebutton(const std::string &button_name) const {

   std::cout << "-------------------- debug in model_fit_refine_unactive_togglebutton() " << button_name << std::endl;

   if (model_fit_refine_dialog) {
      // GtkWidget *toggle_button = lookup_widget(model_fit_refine_dialog, button_name.c_str());
      GtkWidget *toggle_button = widget_from_builder(button_name.c_str());
      if (toggle_button)
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
      else
	 std::cout << "ERROR:: in model_fit_refine_unactive_togglebutton() failed to find button: "
                   << button_name << std::endl;

   } else {
      // std::cout << "DEBUG:: model_fit_refine_dialog not found" << std::endl;
   }


   std::string toolbar_button_name = "not-found";
   if (button_name == "model_refine_dialog_refine_togglebutton")
      toolbar_button_name = "model_toolbar_refine_togglebutton";
   if (button_name == "model_refine_dialog_regularize_zone_togglebutton")
      toolbar_button_name = "model_toolbar_regularize_togglebutton";
   if (button_name == "model_refine_dialog_rigid_body_togglebutton")
      toolbar_button_name = "model_toolbar_rigid_body_fit_togglebutton";
   if (button_name == "model_refine_dialog_rot_trans_togglebutton")
      toolbar_button_name = "model_toolbar_rot_trans_toolbutton"; // 20220304-PE what should this be?
   if (button_name == "model_refine_dialog_auto_fit_rotamer_togglebutton")
      toolbar_button_name = "model_toolbar_auto_fit_rotamer_togglebutton";
   if (button_name == "model_refine_dialog_rotamer_togglebutton")
      toolbar_button_name = "model_toolbar_rotamers_togglebutton";
   if (button_name == "model_refine_dialog_edit_chi_angles_togglebutton")
      toolbar_button_name = "model_toolbar_edit_chi_angles_togglebutton";
   if (button_name == "model_refine_dialog_torsion_general_togglebutton")
      toolbar_button_name = "model_toolbar_torsion_general_toggletoolbutton";
   if (button_name == "model_refine_dialog_pepflip_togglebutton")
      toolbar_button_name = "model_toolbar_flip_peptide_togglebutton";
   if (button_name == "model_refine_dialog_do_180_degree_sidechain_flip_togglebutton")
      toolbar_button_name = "model_toolbar_sidechain_180_togglebutton";
   if (button_name == "model_refine_dialog_edit_backbone_torsions_togglebutton")
      toolbar_button_name = "model_toolbar_edit_backbone_torsions_toggletoolbutton";
   if (button_name == "model_refine_dialog_mutate_auto_fit_togglebutton")
      toolbar_button_name = "model_toolbar_mutate_and_autofit_togglebutton";
   if (button_name == "model_refine_dialog_mutate_togglebutton")
      toolbar_button_name = "model_toolbar_simple_mutate_togglebutton";
   if (button_name == "model_refine_dialog_fit_terminal_residue_togglebutton")
      toolbar_button_name = "model_toolbar_add_terminal_residue_togglebutton";

   std::cout << "-------------------- debug in model_fit_refine_unactive_togglebutton() toolbar_button_name "
             << toolbar_button_name << std::endl;

   // now, button_name may have been
   // model_refine_dialog_edit_phi_psi_togglebutton or
   // model_refine_dialog_edit_backbone_torsions_togglebutton, we
   // don't have toolbar equivalents of those.
   //
   if (toolbar_button_name != "not-found") {
      // GtkWidget *toggle_button = lookup_widget(graphics_info_t::get_main_window(), toolbar_button_name.c_str());
      GtkWidget *toggle_button = widget_from_builder(toolbar_button_name.c_str());
//       std::cout << "DEBUG:: toggle_button for gtk2 toolbar: " << button_name << "->"
// 		<< toolbar_button_name << " " << toggle_button << std::endl;

      if (toggle_button) {
	// somehow we cannot use ->active on the toggle_tool_buttons?!
	gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggle_button));
	if (active)
	  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
      }
   }

}


void
graphics_info_t::other_modelling_tools_unactive_togglebutton(const std::string &button_name) const {

   if (other_modelling_tools_dialog) {
      // GtkWidget *toggle_button = lookup_widget(other_modelling_tools_dialog, button_name.c_str());
      GtkWidget *toggle_button = widget_from_builder(button_name.c_str());
      if (toggle_button)
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
      else
	 std::cout << "ERROR:: failed to find button: " << button_name
		   << std::endl;
   }
}

int
graphics_info_t::wrapped_create_edit_chi_angles_dialog(const std::string &res_type,
						       edit_chi_edit_type mode) {

   // mode is either EDIT_CHI or RESIDUE_PARTIAL_ALT_LOCS.

   // GtkWidget *dialog = create_edit_chi_angles_dialog();
   GtkWidget *dialog = widget_from_builder("edit_chi_angles_dialog");
   if (mode == RESIDUE_PARTIAL_ALT_LOCS) {
      gtk_window_set_title(GTK_WINDOW(dialog), "Add Alternative Conformer Split by Torsion");
   }

   set_transient_and_position(COOT_EDIT_CHI_DIALOG, dialog);

   // Fill the vbox with buttons with atom labels about which there
   // are rotatable torsions:
   //
   GtkWidget *vbox = widget_from_builder("edit_chi_angles_vbox");

   clear_out_container(vbox);

   std::cout << "debug:: calling fill_chi_angles_vbox() with mode " << mode << std::endl;
   int n_boxes = fill_chi_angles_vbox(vbox, res_type, mode);

   // this needs to be deleted when dialog is destroyed, I think,
   // (currently it isn't).
   int NCHAR = 100;
   char *s = new char[NCHAR];
   for (int i=0; i<NCHAR; i++) s[i] = 0;
   strncpy(s, res_type.c_str(), NCHAR-1);
   g_object_set_data(G_OBJECT(vbox), "res_type", s);
   // we get this in c-interface-gui.cc's fill_chi_angles_vbox();

   gtk_widget_set_visible(dialog, TRUE);

   // and now the hydrogen torsion checkbutton:
   GtkWidget *checkbutton = widget_from_builder("edit_chi_angles_add_hydrogen_torsions_checkbutton");

   if (find_hydrogen_torsions_flag)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);

   // "Reverse Fragment" button
   edit_chi_angles_reverse_fragment = 0; // reset the static

   return n_boxes;
}

void
graphics_info_t::clear_out_container(GtkWidget *box) {

   if (box) {
      GtkWidget *item_widget = gtk_widget_get_first_child(box);
      while (item_widget) {
         GtkWidget *w = item_widget;
         item_widget = gtk_widget_get_next_sibling(item_widget);
         gtk_box_remove(GTK_BOX(box), w);
      };
   }
}


int
graphics_info_t::fill_chi_angles_vbox(GtkWidget *vbox, std::string monomer_type,
				      edit_chi_edit_type mode) {

   int imol = 0; // FIXME - extract this from the vbox
   //
   // g_object_get_data(G_OBJECT(vbox), ...);
   // g_object_set_data(..) elsewhere is needed of course

   int n_non_const_torsions = -1; // unset

   clear_out_container(vbox);

   std::pair<short int, coot::dictionary_residue_restraints_t> p =
      Geom_p()->get_monomer_restraints(monomer_type, imol);

   if (p.first) {

      std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
	 p.second.get_non_const_torsions(find_hydrogen_torsions_flag);
      n_non_const_torsions = torsion_restraints.size();

      // We introduce here ichi (which gets incremented if the current
      // torsion is not const), we do that so that we have consistent
      // indexing in the torsions vector with chi_angles's change_by()
      // (see comments above execute_edit_chi_angles()).

      int ichi = 0;
      for (unsigned int i=0; i<torsion_restraints.size(); i++) {
	 if (!torsion_restraints[i].is_const()) {
	    std::string label = "  ";
	    label += torsion_restraints[i].id();
	    label += "  ";
	    label += torsion_restraints[i].atom_id_2_4c();
	    label += " <--> ";
	    label += torsion_restraints[i].atom_id_3_4c();
	    label += "  ref: ";
	    label += coot::util::float_to_string(torsion_restraints[i].angle());
	    label += "  per: ";
	    label += coot::util::int_to_string(torsion_restraints[i].periodicity());
	    GtkWidget *button = gtk_button_new_with_label(label.c_str());
            gtk_widget_set_margin_top(button, 2);
            gtk_widget_set_margin_bottom(button, 2);
            gtk_widget_set_margin_start(button, 6);
            gtk_widget_set_margin_end(button, 6);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
            // 20220528-PE FIXME events
#else
	    gtk_widget_set_events(GTK_WIDGET(button),
				  GDK_EXPOSURE_MASK      |
				  GDK_BUTTON_PRESS_MASK  |
				  GDK_BUTTON_RELEASE_MASK|
				  GDK_SCROLL_MASK        |
				  GDK_POINTER_MOTION_MASK|
				  GDK_POINTER_MOTION_HINT_MASK);
#endif
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(on_change_current_chi_button_clicked),
			     GINT_TO_POINTER(ichi));
	    g_signal_connect(G_OBJECT(button), "enter",
			     G_CALLBACK(on_change_current_chi_button_entered),
			     GINT_TO_POINTER(ichi));

#ifdef FIX_THE_KEY_PRESS_EVENTS
	    g_signal_connect(G_OBJECT(button), "motion_notify_event",
			     G_CALLBACK(on_change_current_chi_motion_notify),
			     GINT_TO_POINTER(ichi));
#endif

	    int int_mode = static_cast<int> (mode);
	    g_object_set_data(G_OBJECT(button), "i_bond", GINT_TO_POINTER(ichi));
	    g_object_set_data(G_OBJECT(button), "chi_edit_mode", GINT_TO_POINTER(int_mode));

	    gtk_widget_set_name(button, "edit_chi_angles_button");

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
	    gtk_box_append(GTK_BOX(vbox), button);
#else
	    gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
	    gtk_container_set_border_width(GTK_CONTAINER(button), 2);
#endif
	    gtk_widget_set_visible(button, TRUE);
	    ichi++;
	 }
      }
   }
   return n_non_const_torsions;
}

// static
void
graphics_info_t::on_change_current_chi_button_clicked(GtkButton *button,
						      gpointer user_data) {

   graphics_info_t g;
   int i = GPOINTER_TO_INT(user_data);
   g.edit_chi_current_chi = i + 1;
   g.in_edit_chi_mode_flag = 1;

   int i_mode = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "chi_edit_mode"));
   edit_chi_edit_type mode = static_cast<edit_chi_edit_type> (i_mode);

   int i_bond = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "i_bond"));
   std::cout << "DEBUG:: g.on_change_current_chi_button_clicked() "
	     << g.edit_chi_current_chi << " mode " << mode
	     << " i_bond " << i_bond << std::endl;

   if (mode == RESIDUE_PARTIAL_ALT_LOCS) {

      // wag the dog: when true, split the large fragment, not the small one.
      bool wag_the_dog = false; // test for Ctrl pressed to set true

      // to get access to event, I need to use a callback for button-release
      // not button clicked
      //
      // if (event->state & GDK_CONTROL_MASK)
      //    wag_the_dog = true;

      g.residue_partial_alt_locs_split_residue(i_bond, wag_the_dog);
      graphics_draw();
   }
}

// static
void
graphics_info_t::on_change_current_chi_button_entered(GtkButton *button,
						      gpointer user_data) {

   // graphics_info_t g;
   // int ibond_user = GPOINTER_TO_INT(user_data);
   // g.setup_flash_bond_internal(ibond_user);

}


// Create a moving atoms molecule, consisting of the Ca(n), Ca(n+1) of
// the peptide, N(n) C(n+1), O(n+1).  Note the alt conf should be the
// same (we can have altconfed mainchain).
//
void
graphics_info_t::execute_setup_backbone_torsion_edit(int imol, int atom_index) {

   // if not the same altconf for all atoms, give up (for now).
   //
   // Do an atom selection of this residue and the previous one
   // (unless this was a N, then do this one and the next)
   //
   // Run through the atom selection finding the atoms.  put the atoms
   // into new residues (appropriate) and construct a mol and and asc
   // and make that the moving atoms asc.  Be sure to put the residue
   // with the peptide N first so that mmdb finds it first when it
   // does the atom selection that's part of make_asc(mol) - which
   // means that get_first_atom_with_atom_name() will work like we
   // want it to.

   if (imol < n_molecules()) {
      if (molecules[imol].has_model()) {
	 if (atom_index < molecules[imol].atom_sel.n_selected_atoms) {
	    mmdb::Atom *this_atom_p = molecules[imol].atom_sel.atom_selection[atom_index];
	    int offset = 0; // usually this and next residue
	    std::string this_atname(this_atom_p->name);
	    if (this_atname == " N  ") {
	       offset = -1;
	    }
	    int this_res = this_atom_p->GetSeqNum();
	    char *chain_id = this_atom_p->GetChainID();
	    int SelectionHandle = molecules[imol].atom_sel.mol->NewSelection();
	    char *ins_code = this_atom_p->GetInsCode();
	    char *altconf  = this_atom_p->altLoc;
	    std::string a_tmp(altconf);
	    if (a_tmp != "") {
	       a_tmp += ",";
	    }
	    char *search_altconf = (char *) a_tmp.c_str();
	    molecules[imol].atom_sel.mol->SelectAtoms (SelectionHandle, 0, chain_id,
					      this_res + offset , // starting resno
					      ins_code, // any insertion code
					      this_res + 1 + offset, // ending resno
					      ins_code, // ending insertion code
					      "*", // any residue name
					      "*", // atom name
					      "*", // elements
					      search_altconf  // alt loc.
					      );
	    int nSelAtoms;
	    mmdb::PPAtom SelectAtoms;
	    molecules[imol].atom_sel.mol->GetSelIndex(SelectionHandle, SelectAtoms, nSelAtoms);
	    if (nSelAtoms < 4) { // the very min
	       std::cout <<  "WARNING:: not enough atoms in atom selection in "
			 <<  "execute_setup_backbone_torsion_edit" << std::endl;
	    } else {


	       // We construct a moving atom asc atom by atom...
	       mmdb::Atom *next_ca = NULL, *next_n = NULL;
	       mmdb::Atom *this_c = NULL, *this_o = NULL, *this_ca = NULL;

	       // and the extra atoms that we need (we extract the
	       // coordinates) to construct a pair of ramachandran
	       // points in the rama_plot:

	       mmdb::Atom *prev_c = NULL;
	       mmdb::Atom *this_n = NULL;
	       mmdb::Atom *next_c = NULL;
	       mmdb::Atom *next_plus_1_n = NULL;

	       // to get prev_c and next_plus_1_n we do atom selections:

	       // You can add (this) hydrogen to that list if you want.
	       mmdb::Atom *at;
	       //
	       for (int iat=0; iat<nSelAtoms; iat++) {
		  at = SelectAtoms[iat];
		  if (at->GetSeqNum() == (this_res + offset) ) {
		     std::string n(at->name);
		     if (n == " CA ") {
			this_ca = at;
		     }
		     if (n == " O  ") {
			this_o = at;
		     }
		     if (n == " C  ") {
			this_c = at;
		     }
		     if (n == " N  ") {
			this_n = at;
		     }
		  }
		  if (at->GetSeqNum() == (this_res + 1 + offset) ) {
		     std::string n(at->name);
		     if (n == " CA ") {
			next_ca = at;
		     }
		     if (n == " N  ") {
			next_n = at;
		     }
		     if (n == " C  ") {
			next_c = at;
		     }
		  }
	       }

	       int SelHnd_prev_c = molecules[imol].atom_sel.mol->NewSelection();
	       // to get prev_c and next_plus_1_n we do atom selections:
	       molecules[imol].atom_sel.mol->SelectAtoms (SelHnd_prev_c, 0, chain_id,
							  this_res - 1 + offset , // starting resno
							  ins_code, // any insertion code
							  this_res - 1 + offset, // ending resno
							  ins_code, // ending insertion code
							  "*",    // any residue name
							  " C  ", // atom name
							  "*",    // elements
							  search_altconf  // alt loc.
							  );

	       int nSelAtoms_prev_c;
	       mmdb::PPAtom SelectAtoms_prev_c;

	       molecules[imol].atom_sel.mol->GetSelIndex(SelHnd_prev_c,
							 SelectAtoms_prev_c,
							 nSelAtoms_prev_c);
	       if (nSelAtoms_prev_c > 0) {
		  prev_c = SelectAtoms_prev_c[0];
	       } else {
		  std::cout << "Oops:: didn't find prev_c\n";
	       }
	       molecules[imol].atom_sel.mol->DeleteSelection(SelHnd_prev_c);

	       // And similarly for next_plus_1_n:
	       int SelHnd_next_plus_1_n = molecules[imol].atom_sel.mol->NewSelection();
	       molecules[imol].atom_sel.mol->SelectAtoms (SelHnd_next_plus_1_n, 0, chain_id,
							  this_res + 2 + offset , // starting resno
							  ins_code, // any insertion code
							  this_res + 2 + offset, // ending resno
							  ins_code, // ending insertion code
							  "*",    // any residue name
							  " N  ", // atom name
							  "*",    // elements
							  search_altconf  // alt loc.
							  );
	       int nSelAtoms_next_plus_1_n;
	       mmdb::PPAtom SelectAtoms_next_plus_1_n;

	       molecules[imol].atom_sel.mol->GetSelIndex(SelHnd_next_plus_1_n,
							 SelectAtoms_next_plus_1_n,
							 nSelAtoms_next_plus_1_n);
	       if (nSelAtoms_next_plus_1_n > 0) {
		  next_plus_1_n = SelectAtoms_next_plus_1_n[0];
	       } else {
		  std::cout << "Oops:: didn't find next + 1 N\n";
	       }
	       molecules[imol].atom_sel.mol->DeleteSelection(SelHnd_next_plus_1_n);


	       if (next_ca && next_n && this_ca && this_o && this_c) {

		  // new addition 25Feb2004
		  rama_plot_for_2_phi_psis(imol, atom_index);

		  mmdb::Manager *mol = new mmdb::Manager;
		  mmdb::Model *model = new mmdb::Model;
		  mmdb::Chain *chain = new mmdb::Chain;
		  mmdb::Residue *res1 = new mmdb::Residue;
		  mmdb::Residue *res2 = new mmdb::Residue;
		  res1->SetResName(this_ca->GetResName());
		  res2->SetResName(next_ca->GetResName());
		  res1->seqNum = this_ca->GetSeqNum();
		  res2->seqNum = next_ca->GetSeqNum();
		  chain->SetChainID(this_ca->GetChainID());

		  at = new mmdb::Atom;
		  at->Copy(this_ca);
		  res1->AddAtom(at);

		  at = new mmdb::Atom;
		  at->Copy(this_c);
		  res1->AddAtom(at);

		  at = new mmdb::Atom;
		  at->Copy(this_o);
		  res1->AddAtom(at);

		  at = new mmdb::Atom;
		  at->Copy(next_ca);
		  res2->AddAtom(at);

		  at = new mmdb::Atom;
		  at->Copy(next_n);
		  res2->AddAtom(at);

		  chain->AddResidue(res1);
		  chain->AddResidue(res2);
		  model->AddChain(chain);
		  mol->AddModel(model);
		  mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
		  mol->FinishStructEdit();
		  imol_moving_atoms = imol;
		  moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
		  atom_selection_container_t asc = make_asc(mol);
		  regularize_object_bonds_box.clear_up();
		  make_moving_atoms_graphics_object(imol, asc);

		  // save the fixed end points:
		  backbone_torsion_end_ca_1 =
		     clipper::Coord_orth(this_ca->x, this_ca->y, this_ca->z);
		  backbone_torsion_end_ca_2 =
		     clipper::Coord_orth(next_ca->x, next_ca->y, next_ca->z);


		  // add to rama_points:
		  rama_points.clear();
		  rama_points.add("this_ca", backbone_torsion_end_ca_1);
		  rama_points.add("next_ca", backbone_torsion_end_ca_2);
		  // the moving atoms:
		  rama_points.add("this_c", clipper::Coord_orth(this_c->x,
								this_c->y,
								this_c->z));
		  rama_points.add("this_o", clipper::Coord_orth(this_o->x,
								this_o->y,
								this_o->z));
		  rama_points.add("next_n", clipper::Coord_orth(next_n->x,
								next_n->y,
								next_n->z));

		  if (this_n)
		     rama_points.add("this_n",  clipper::Coord_orth(this_n->x,
								    this_n->y,
								    this_n->z));

		  if (next_c)
		     rama_points.add("next_c",  clipper::Coord_orth(next_c->x,
								    next_c->y,
								    next_c->z));

		  // next_plus_1_n, prev_c;
		  if (next_plus_1_n)
		     rama_points.add("next+1_n", clipper::Coord_orth(next_plus_1_n->x,
								     next_plus_1_n->y,
								     next_plus_1_n->z));
		  if (prev_c)
		     rama_points.add("prev_c", clipper::Coord_orth(prev_c->x,
								   prev_c->y,
								   prev_c->z));

// 		  std::cout << "DEBUG:: backbone_torsion_end_ca_1: "
// 			    << backbone_torsion_end_ca_1.format() << std::endl;
// 		  std::cout << "DEBUG:: backbone_torsion_end_ca_2: "
// 			    << backbone_torsion_end_ca_2.format() << std::endl;

		  graphics_draw();
		  // GtkWidget *widget = create_edit_backbone_torsions_dialog();
		  GtkWidget *widget = widget_from_builder("edit_backbone_torsions_dialog");
		  set_edit_backbone_adjustments(widget);
		  gtk_widget_set_visible(widget, TRUE);
        // update the graph to show both
        GtkWidget *adj;
        adj = widget_from_builder("edit_backbone_torsions_rotate_carbonyl_adjustment");
        g_signal_emit_by_name(G_OBJECT(adj), "value_changed");

	       } else {
		  std::cout << "WARNING:: not all atoms found in "
			    << "execute_setup_backbone_torsion_edit" << std::endl;

	       }
	    }
	    //
	    molecules[imol].atom_sel.mol->DeleteSelection(SelectionHandle);
	 }
      }
   }
}

void
graphics_info_t::set_edit_backbone_adjustments(GtkWidget *widget) {

   GtkWidget *hscale_peptide = widget_from_builder("edit_backbone_torsions_rotate_peptide_hscale");

   GtkWidget *hscale_carbonyl = widget_from_builder("edit_backbone_torsions_rotate_carbonyl_hscale");

//    gfloat value,
//    gfloat lower,
//    gfloat upper,
//    gfloat step_increment,
//    gfloat page_increment,
//    gfloat page_size

   GtkAdjustment *adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -180.0, 360.0, 0.1, 1.0, 180));
   gtk_range_set_adjustment(GTK_RANGE(hscale_peptide), adjustment);
   g_signal_connect(G_OBJECT(adjustment), "value_changed",
		    G_CALLBACK(graphics_info_t::edit_backbone_peptide_changed_func), NULL);


   // and the carbonyl:
   adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -180.0, 360.0, 0.1, 1.0, 180));
   gtk_range_set_adjustment(GTK_RANGE(hscale_carbonyl), adjustment);
   g_signal_connect(G_OBJECT(adjustment), "value_changed",
		    G_CALLBACK(graphics_info_t::edit_backbone_carbonyl_changed_func), NULL);

   // add a name to be able to lookup
   g_object_set_data(G_OBJECT(widget), "edit_backbone_torsions_rotate_carbonyl_adjustment",
                     adjustment);

}

// static
void
graphics_info_t::edit_backbone_peptide_changed_func(GtkAdjustment *adj, GtkWidget *window) {

#ifdef HAVE_GOOCANVAS

   graphics_info_t g;
   // std::cout << "change backbone peptide by: " << adj->value << std::endl;

   std::pair<short int, clipper::Coord_orth> this_c = rama_points.get("this_c");
   std::pair<short int, clipper::Coord_orth> this_o = rama_points.get("this_o");
   std::pair<short int, clipper::Coord_orth> next_n = rama_points.get("next_n");

   if (this_c.first && this_o.first && next_n.first) {
      mmdb::Atom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);
      mmdb::Atom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
      mmdb::Atom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);

      double rad_angle = clipper::Util::d2rad(gtk_adjustment_get_value(adj));
      clipper::Coord_orth new_c =
	 coot::util::rotate_around_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
				this_c.second,
				backbone_torsion_end_ca_1, rad_angle);
      clipper::Coord_orth new_o =
	 coot::util::rotate_around_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
					  this_o.second,
					  backbone_torsion_end_ca_1, rad_angle);
      clipper::Coord_orth new_n =
	 coot::util::rotate_around_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
					  next_n.second,
					  backbone_torsion_end_ca_1, rad_angle);
      n_atom_p->x = new_n.x();
      n_atom_p->y = new_n.y();
      n_atom_p->z = new_n.z();

      c_atom_p->x = new_c.x();
      c_atom_p->y = new_c.y();
      c_atom_p->z = new_c.z();

      o_atom_p->x = new_o.x();
      o_atom_p->y = new_o.y();
      o_atom_p->z = new_o.z();

      std::pair<std::pair<double, double>, std::pair<double, double> > pp =
	 g.phi_psi_pairs_from_moving_atoms();

//    std::cout << pp.first.first  << " " << pp.first.second << "      "
// 	     << pp.second.first << " " << pp.second.second << std::endl;

      if (edit_phi_psi_plot) {
	 std::vector <coot::util::phi_psi_t> vp;
	 std::string label = int_to_string(c_atom_p->GetSeqNum());
	 if (pp.first.first > -200) {
	    label += c_atom_p->GetChainID();
	    coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first),
					  clipper::Util::rad2d(pp.first.second),
					  "resname", label, 1, "inscode", "chainid");
	    vp.push_back(phipsi1);
	 }
	 if (pp.second.first > -200) {
	    label = int_to_string(n_atom_p->GetSeqNum());
	    label += n_atom_p->GetChainID();
	    coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first),
					  clipper::Util::rad2d(pp.second.second),
					  "resname", label, 1, "inscode", "chainid");

	    vp.push_back(phipsi2);
	 }
	 if (vp.size() > 0)
	    edit_phi_psi_plot->draw_it(vp);
      }

      regularize_object_bonds_box.clear_up();
      int imol = 0; // should be fine for backbone edits
      g.make_moving_atoms_graphics_object(imol, *moving_atoms_asc);
      graphics_draw();

   } else {
      std::cout << "ERROR:: can't find rama points in edit_backbone_peptide_changed_func"
		<< std::endl;
   }

#endif
}

// static
void
graphics_info_t::edit_backbone_carbonyl_changed_func(GtkAdjustment *adj, GtkWidget *window) {

#ifdef HAVE_GOOCANVAS
   graphics_info_t g;
   // std::cout << "change backbone peptide by: " << adj->value << std::endl;

   std::pair<short int, clipper::Coord_orth> this_c = rama_points.get("this_c");
   std::pair<short int, clipper::Coord_orth> this_o = rama_points.get("this_o");

   if (this_c.first && this_o.first) {
      mmdb::Atom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
      mmdb::Atom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);
      mmdb::Atom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);

      clipper::Coord_orth carbonyl_n_pos(n_atom_p->x, n_atom_p->y, n_atom_p->z);

      double rad_angle = clipper::Util::d2rad(gtk_adjustment_get_value(adj));
      clipper::Coord_orth new_c =
	 coot::util::rotate_around_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
					  this_c.second,
					  backbone_torsion_end_ca_1, rad_angle);
      clipper::Coord_orth new_o =
	 coot::util::rotate_around_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
					  this_o.second,
					  backbone_torsion_end_ca_1, rad_angle);

      c_atom_p->x = new_c.x();
      c_atom_p->y = new_c.y();
      c_atom_p->z = new_c.z();

      o_atom_p->x = new_o.x();
      o_atom_p->y = new_o.y();
      o_atom_p->z = new_o.z();

      std::pair<std::pair<double, double>, std::pair<double, double> > pp =
	 g.phi_psi_pairs_from_moving_atoms();

//    std::cout << pp.first.first  << " " << pp.first.second << "      "
// 	     << pp.second.first << " " << pp.second.second << std::endl;


      if (edit_phi_psi_plot) {
	 std::vector <coot::util::phi_psi_t> vp;
	 std::string label = int_to_string(c_atom_p->GetSeqNum());
	 if (pp.first.first > -200) {
       label += " ";
       label += c_atom_p->GetChainID();
       label += " ";
       label += c_atom_p->GetResName();
	    coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first),
					  clipper::Util::rad2d(pp.first.second),
					  "resname", label, 1, "inscode", "chainid");
	    vp.push_back(phipsi1);
	 }
	 if (pp.second.first > -200) {
	    label = int_to_string(n_atom_p->GetSeqNum());
       label += " ";
	    label += n_atom_p->GetChainID();
       label += " ";
       label += n_atom_p->GetResName();
	    coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first),
					  clipper::Util::rad2d(pp.second.second),
					  "resname", label, 1, "inscode", "chainid");

	    vp.push_back(phipsi2);
	 }
	 if (vp.size() > 0)
	    edit_phi_psi_plot->draw_it(vp);
      }


      regularize_object_bonds_box.clear_up();
      int imol = 0; // should be fine for backbone edits
      g.make_moving_atoms_graphics_object(imol, *moving_atoms_asc);
      graphics_draw();

   } else {
      std::cout << "ERROR:: can't find rama points in edit_backbone_peptide_changed_func"
		<< std::endl;
   }

#endif
}


// #include "mmdb-extras.h"

// Tinker with the moving atoms
void
graphics_info_t::change_peptide_carbonyl_by(double angle) {

#ifdef HAVE_GOOCANVAS

//    std::cout << "move carbonyl by " << angle << std::endl;
   mmdb::Atom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);
   mmdb::Atom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
   mmdb::Atom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);

   clipper::Coord_orth carbonyl_n_pos(n_atom_p->x, n_atom_p->y, n_atom_p->z);
   clipper::Coord_orth carbonyl_c_pos(c_atom_p->x, c_atom_p->y, c_atom_p->z);
   clipper::Coord_orth carbonyl_o_pos(o_atom_p->x, o_atom_p->y, o_atom_p->z);

   double rad_angle = clipper::Util::d2rad(angle);

   clipper::Coord_orth new_c =
      coot::util::rotate_around_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
			  carbonyl_c_pos,
			  carbonyl_n_pos, rad_angle);
   clipper::Coord_orth new_o =
      coot::util::rotate_around_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
			  carbonyl_o_pos,
			  carbonyl_n_pos, rad_angle);

   c_atom_p->x = new_c.x();
   c_atom_p->y = new_c.y();
   c_atom_p->z = new_c.z();

   o_atom_p->x = new_o.x();
   o_atom_p->y = new_o.y();
   o_atom_p->z = new_o.z();

   std::pair<std::pair<double, double>, std::pair<double, double> > pp =
      phi_psi_pairs_from_moving_atoms();

//    std::cout << pp.first.first  << " " << pp.first.second << "      "
// 	     << pp.second.first << " " << pp.second.second << std::endl;



   if (edit_phi_psi_plot) {
      std::vector <coot::util::phi_psi_t> vp;
      std::string label = int_to_string(c_atom_p->GetSeqNum());
      label += c_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first),
				    clipper::Util::rad2d(pp.first.second),
				    "resname", label, 1, "inscode", "chainid");
      label = int_to_string(n_atom_p->GetSeqNum());
      label += n_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first),
				    clipper::Util::rad2d(pp.second.second),
				    "resname", label, 1, "inscode", "chainid");

      vp.push_back(phipsi1);
      vp.push_back(phipsi2);
      edit_phi_psi_plot->draw_it(vp);
   }


   regularize_object_bonds_box.clear_up();
   int imol = 0; // should be fine for backbone edits
   make_moving_atoms_graphics_object(imol, *moving_atoms_asc);
   graphics_draw();

#endif // HAVE_GOOCANVAS
}


// Return a pair of phi psi's for the residue of the carbonyl and the next
// residue.
//
// If for some reason the pair is incalculable, put phi (first) to -2000.
//
std::pair<std::pair<double, double>, std::pair<double, double> >
graphics_info_t::phi_psi_pairs_from_moving_atoms() {

   std::pair<std::pair<double, double>, std::pair<double, double> > p;

   // There are 2 atoms in the moving atoms that we need for this calculation,
   // this_C and next_N.  The rest we will pull out of a pre-stored vector of
   // pairs (made in execute_setup_backbone_torsion_edit): rama_points.

   mmdb::Atom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
   mmdb::Atom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);

   clipper::Coord_orth this_c (c_atom_p->x, c_atom_p->y, c_atom_p->z);
   clipper::Coord_orth next_n (n_atom_p->x, n_atom_p->y, n_atom_p->z);

   std::pair<short int, clipper::Coord_orth> prev_c      = rama_points.get("prev_c");
   std::pair<short int, clipper::Coord_orth> this_ca     = rama_points.get("this_ca");
   std::pair<short int, clipper::Coord_orth> this_n      = rama_points.get("this_n");
   std::pair<short int, clipper::Coord_orth> next_ca     = rama_points.get("next_ca");
   std::pair<short int, clipper::Coord_orth> next_c      = rama_points.get("next_c");
   std::pair<short int, clipper::Coord_orth> next_plus_n = rama_points.get("next+1_n");

   if (prev_c.first && this_ca.first && this_n.first) {

      // we can calculate the first ramachadran phi/psi pair
      double phi = clipper::Coord_orth::torsion(prev_c.second, this_n.second,  this_ca.second, this_c);
      double psi = clipper::Coord_orth::torsion(this_n.second, this_ca.second, this_c,         next_n);

      p.first.first  = phi;
      p.first.second = psi;

   } else {

      // can't get the first ramachandran point
      p.first.first = -2000;

   }

   if (next_ca.first && next_c.first && next_plus_n.first) {

      // we can calculate the first ramachadran phi/psi pair

      double phi = clipper::Coord_orth::torsion(this_c, next_n, next_ca.second, next_c.second);
      double psi = clipper::Coord_orth::torsion(next_n, next_ca.second, next_c.second, next_plus_n.second);

      p.second.first  = phi;
      p.second.second = psi;

   } else {

      // can't get the second ramachandran point
      p.second.first = -2000;
   }


   return p;
}
// Tinker with the moving atoms
void
graphics_info_t::change_peptide_peptide_by(double angle) {

#ifdef HAVE_GOOCANVAS

   //    std::cout << "move peptide by " << angle << std::endl;

   mmdb::Atom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);
   mmdb::Atom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
   mmdb::Atom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);

   clipper::Coord_orth carbonyl_n_pos(n_atom_p->x, n_atom_p->y, n_atom_p->z);
   clipper::Coord_orth carbonyl_c_pos(c_atom_p->x, c_atom_p->y, c_atom_p->z);
   clipper::Coord_orth carbonyl_o_pos(o_atom_p->x, o_atom_p->y, o_atom_p->z);

   double rad_angle = clipper::Util::d2rad(angle);

   clipper::Coord_orth new_c =
      coot::util::rotate_around_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			  carbonyl_c_pos,
			  backbone_torsion_end_ca_1, rad_angle);
   clipper::Coord_orth new_o =
      coot::util::rotate_around_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			  carbonyl_o_pos,
			  backbone_torsion_end_ca_1, rad_angle);

   clipper::Coord_orth new_n =
      coot::util::rotate_around_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			  carbonyl_n_pos,
			  backbone_torsion_end_ca_1, rad_angle);

   n_atom_p->x = new_n.x();
   n_atom_p->y = new_n.y();
   n_atom_p->z = new_n.z();

   c_atom_p->x = new_c.x();
   c_atom_p->y = new_c.y();
   c_atom_p->z = new_c.z();

   o_atom_p->x = new_o.x();
   o_atom_p->y = new_o.y();
   o_atom_p->z = new_o.z();

   std::pair<std::pair<double, double>, std::pair<double, double> > pp =
      phi_psi_pairs_from_moving_atoms();

//    std::cout << pp.first.first  << " " << pp.first.second << "      "
// 	     << pp.second.first << " " << pp.second.second << std::endl;

#ifdef DO_RAMA_PLOT

   if (edit_phi_psi_plot) {
      std::vector <coot::util::phi_psi_t> vp;
      std::string label = int_to_string(c_atom_p->GetSeqNum());
      label += c_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first),
				    clipper::Util::rad2d(pp.first.second),
				    "resname", label, 1, "inscode", "chainid");
      label = int_to_string(n_atom_p->GetSeqNum());
      label += n_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first),
				    clipper::Util::rad2d(pp.second.second),
				    "resname", label, 1, "inscode", "chainid");

      vp.push_back(phipsi1);
      vp.push_back(phipsi2);
      edit_phi_psi_plot->draw_it(vp);
   }
#endif


   regularize_object_bonds_box.clear_up();
   int imol = 0; // should be fine for backbone edits
   make_moving_atoms_graphics_object(imol, *moving_atoms_asc);
   graphics_draw();

#endif // HAVE_GOOCANVAS
}


void
graphics_info_t::set_backbone_torsion_peptide_button_start_pos(int ix, int iy) {
   backbone_torsion_peptide_button_start_pos_x = ix;
   backbone_torsion_peptide_button_start_pos_y = iy;
}

void
graphics_info_t::set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy) {
   backbone_torsion_carbonyl_button_start_pos_x = ix;
   backbone_torsion_carbonyl_button_start_pos_y = iy;
}


void
graphics_info_t::change_peptide_peptide_by_current_button_pos(int ix, int iy) {

   double diff = 0.05* (ix - backbone_torsion_peptide_button_start_pos_x);
   change_peptide_peptide_by(diff);
}

void
graphics_info_t::change_peptide_carbonyl_by_current_button_pos(int ix, int iy) {

   double diff = 0.05 * (ix - backbone_torsion_carbonyl_button_start_pos_x);
   change_peptide_carbonyl_by(diff);
}

//      ----------------- sequence view ----------------

// 20211201-PE do we need this function!?

#if 0 // was defined (HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
//
// return NULL on failure to find the class for this molecule
coot::sequence_view *
graphics_info_t::get_sequence_view(int imol) {

   GtkWidget *w = NULL;
   coot::sequence_view *r = NULL;

   if (imol < n_molecules()) {
      if (molecules[imol].has_model()) {
         // w = sequence_view_is_displayed[imol];
         w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
         r = (coot::sequence_view *) gtk_object_get_user_data(GTK_OBJECT(w));
         // std::cout << "DEBUG:: user data from " << w << " is " << r << std::endl;
      }
   }
   return r;
}
#endif // HAVE_GTK_CANVAS

void
graphics_info_t::set_sequence_view_is_displayed(GtkWidget *widget, int imol) {

#ifdef HAVE_GOOCANVAS

   if (imol < n_molecules()) {

      // first delete the old sequence view if it exists
      if (false) {  // we don't need to do this because destroying the scrolled window will do it
                    // for us by proxy
         GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
         if (w) {
            gpointer o = g_object_get_data(G_OBJECT(w), "nsv"); // w is the canvas
            std::cout << "got o " << o << std::endl;
            coot::sequence_view *sv = reinterpret_cast<coot::sequence_view *>(o);
            std::cout << "in set_sequence_view_is_displayed() extracted sv " << sv << std::endl;
            delete sv;
         }
      }

//       coot::sequence_view *sv = (coot::sequence_view *)
// 	 gtk_object_get_user_data(GTK_OBJECT(sequence_view_is_displayed[imol]));
//       std::cout << "DEBUG:: seting sequence_view_is_displayed[" << imol
// 		<< "] " << widget << std::endl;
//       sequence_view_is_displayed[imol] = widget; // ols style


      coot::set_validation_graph(imol, coot::SEQUENCE_VIEW, widget);
   }

#endif // HAVE_GOOCANVAS

}

GtkWidget *
graphics_info_t::get_sequence_view_is_displayed(int imol) const {

   GtkWidget *w = 0;
#ifdef HAVE_GOOCANVAS
   if (is_valid_model_molecule(imol))
       w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
#endif
   return w;
}


void
graphics_info_t::unset_geometry_dialog_distance_togglebutton() {

   if (geometry_dialog) {
      GtkWidget *toggle_button = widget_from_builder("geometry_distance_togglebutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
   }
}

void
graphics_info_t::unset_geometry_dialog_dynamic_distance_togglebutton() {

   GtkWidget *toggle_button = widget_from_builder("geometry_dynamic_distance_togglebutton");
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
}


void
graphics_info_t::unset_geometry_dialog_angle_togglebutton() {

   GtkWidget *toggle_button = widget_from_builder("geometry_angle_togglebutton");
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
}


void
graphics_info_t::unset_geometry_dialog_torsion_togglebutton() {

   if (geometry_dialog) {
      GtkWidget *toggle_button = widget_from_builder("geometry_torsion_togglebutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
   }
}


void
graphics_info_t::set_zoom_adjustment(GtkWidget *widget) {

   GtkWidget *zoom_hscale = widget_from_builder("zoom_hscale");

   //  GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -18.0, 36.0, 0.01, 1.0, 18));
   //  GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -2.0, 4.0, 0.01, 1.0, 2));

   GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(zoom, zoom*0.125, zoom*8, 0.01, 0.5, zoom));

   gtk_range_set_adjustment(GTK_RANGE(zoom_hscale), adj);
   g_signal_connect(G_OBJECT(adj), "value_changed",
		    G_CALLBACK(graphics_info_t::zoom_adj_changed),
		    NULL);
}


// static
void
graphics_info_t::zoom_adj_changed(GtkAdjustment *adj, GtkWidget *window) {
   graphics_info_t g;

//    double scaled = pow(M_E, -adj->value);
//    std::cout <<  adj->value << " " << scaled << std::endl;

   // g.zoom *= scaled;
   g.zoom = gtk_adjustment_get_value(adj);
   graphics_draw();

}


// ----------------------------------------------------------------------------
//                check waters
// ----------------------------------------------------------------------------
//
void
graphics_info_t::check_waters_by_difference_map(int imol_waters, int imol_diff_map,
						int interactive_flag) {

   if (is_valid_model_molecule(imol_waters)) {
      if (is_valid_map_molecule(imol_diff_map)) {
	 if (molecules[imol_diff_map].is_difference_map_p()) {
	    std::vector <coot::atom_spec_t> v = molecules[imol_waters].check_waters_by_difference_map(molecules[imol_diff_map].xmap, check_waters_by_difference_map_sigma_level);
	    if (interactive_flag) {
	       GtkWidget *w = wrapped_create_checked_waters_by_variance_dialog(v, imol_waters);
	       gtk_widget_set_visible(w, TRUE);
	    }
	 } else {
	    std::cout << "molecule " <<  imol_diff_map
		      << " is not a difference map\n";
	 }
      } else {
	 std::cout << "molecule " <<  imol_diff_map << "has no map\n";
      }
   } else {
      std::cout << "molecule " <<  imol_waters << "has no model\n";
   }
}

// results widget:
//
//

GtkWidget *
graphics_info_t::wrapped_create_checked_waters_by_variance_dialog(const std::vector <coot::atom_spec_t> &v, int imol) {

   GtkWidget *w;

   if (v.size() > 0) {
      // w = create_interesting_waters_by_difference_map_check_dialog();
      w = widget_from_builder("interesting_waters_by_difference_map_check_dialog");
      GtkWidget *vbox = widget_from_builder("interesting_waters_by_difference_map_check_vbox");
      GtkWidget *button = 0;
      coot::atom_spec_t *atom_spec;
      // GSList *gr_group = NULL;
      GtkWidget *toggle_button_group = nullptr;

      for (unsigned int i=0; i<v.size(); i++) {

	 std::cout << "INFO:: Suspicious water: "
		   << v[i].atom_name
		   << v[i].alt_conf << " "
		   << v[i].res_no << " "
		   << v[i].ins_code << " "
		   << v[i].chain_id << "\n";

	 std::string button_label(" ");
	 button_label += v[i].chain_id;
	 button_label += " " ;
	 button_label += int_to_string(v[i].res_no);
	 button_label += " " ;
	 button_label += v[i].atom_name;
	 button_label += " " ;
	 button_label += v[i].alt_conf;
	 button_label += " " ;

         GtkWidget *toggle_button = gtk_toggle_button_new_with_label(button_label.c_str());

         if (toggle_button_group)
            gtk_toggle_button_set_group(GTK_TOGGLE_BUTTON(toggle_button), GTK_TOGGLE_BUTTON(toggle_button_group));
         else
            toggle_button_group = toggle_button;

	 atom_spec = new coot::atom_spec_t(v[i]);
	 atom_spec->int_user_data = imol;

	 g_signal_connect(G_OBJECT(button), "toggled",
			  G_CALLBACK(on_generic_atom_spec_toggle_button_toggled),
			  atom_spec);

	 gtk_box_append(GTK_BOX(vbox), toggle_button);
         gtk_widget_set_margin_top(toggle_button, 2);
         gtk_widget_set_margin_bottom(toggle_button, 2);
         gtk_widget_set_margin_start(toggle_button, 6);
         gtk_widget_set_margin_end(toggle_button, 6);
      }
   } else {
      std::cout << "INFO:: There are no unusual waters\n";
      std::string s = "There were no strange/anomalous waters\n";
      s += "(in relation to the difference map).";
      w = wrapped_nothing_bad_dialog(s);
   }
   return w;
}


// static
void
   graphics_info_t::on_generic_atom_spec_toggle_button_toggled(GtkToggleButton *toggle_button, gpointer user_data) {

   if (gtk_toggle_button_get_active(toggle_button)) {
      graphics_info_t g;
      coot::atom_spec_t *atom_spec_p = static_cast<coot::atom_spec_t *>(user_data);
      const auto &atom_spec(*atom_spec_p);
      
      g.set_go_to_atom_molecule(atom_spec.int_user_data);
      g.set_go_to_atom_chain_residue_atom_name(atom_spec.chain_id.c_str(),
					       atom_spec.res_no,
					       atom_spec.atom_name.c_str(),
					       atom_spec.alt_conf.c_str());
      g.try_centre_from_new_go_to_atom();
      g.update_things_on_move_and_redraw();
   }

}


GtkWidget *
graphics_info_t::wrapped_create_chiral_restraints_problem_dialog(const std::vector<std::string> &sv) const {

   // 20211011-PE FIXME This needs testing

   // GtkWidget *w = create_chiral_restraints_problem_dialog();
   GtkWidget *w = widget_from_builder("chiral_restraints_problem_dialog");

   GtkWidget *label = widget_from_builder("chiral_volume_restraints_problem_label");
   std::string s = "\n   Problem finding restraints for the following residues:   \n\n";
   for (unsigned int i=0; i<sv.size(); i++) {
      s += sv[i];
      s += "  ";
      if (10*((i+1)/10) == (i+1))
         s += "\n";
   }
   s += "\n";
   gtk_label_set_text(GTK_LABEL(label), s.c_str());
   return w;
}


GtkWidget *
graphics_info_t::wrapped_check_chiral_volumes_dialog(const std::vector <coot::atom_spec_t> &v,
						     int imol) {

   auto make_button_label = [] (const coot::atom_spec_t &as) {
                               // c.f. how we add rotamers: (fill_rotamer_selection_buttons)
                               std::string button_label(" ");
                               button_label += as.chain_id;
                               button_label += " " ;
                               button_label += int_to_string(as.res_no);
                               button_label += " " ;
                               button_label += as.atom_name;
                               button_label += " " ;
                               button_label += as.alt_conf;
                               button_label += " " ;
                               return button_label;
                            };

#if (GTK_MAJOR_VERSION >= 4)
#else
   auto my_delete_items = [] (GtkWidget *widget, void *data) {
                             gtk_container_remove(GTK_CONTAINER(data), widget);
                          };
   auto clear_vbox = [my_delete_items] (GtkWidget *box) {
                        gtk_container_foreach(GTK_CONTAINER(box), my_delete_items, box);
                   };

#endif

   GtkWidget *dialog = NULL;

   std::cout  << "There were " << v.size() << " bad chiral volumes: " << std::endl;

   // 20220602-PE FIXME delete widgets
   std::cout << "in wrapped_check_chiral_volumes_dialog() FIXME delete widgets" << std::endl;

   if (v.size() > 0) {
      dialog = widget_from_builder("bad_chiral_volumes_dialog");
      GtkWidget *bad_chiral_volume_atom_vbox = widget_from_builder("chiral_volume_baddies_vbox");
      // clear_vbox(bad_chiral_volume_atom_vbox);
      for (unsigned int i=0; i<v.size(); i++) {
         if (false)
            std::cout << "  " << v[i].chain_id << " " << v[i].res_no << " " << v[i].atom_name << " "
                      << v[i].alt_conf << " " << "\n";
         std::string button_label(make_button_label(v[i]));
         GtkWidget *button = gtk_button_new_with_label(button_label.c_str());
         coot::atom_spec_t *atom_spec = new coot::atom_spec_t(v[i]);
	 atom_spec->int_user_data = imol;
	 g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(on_inverted_chiral_volume_button_clicked), atom_spec);
	 gtk_box_append(GTK_BOX(bad_chiral_volume_atom_vbox), button);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
	 gtk_container_set_border_width(GTK_CONTAINER(button), 6);
#endif
	 gtk_widget_set_visible(button, TRUE);
      }
      gtk_widget_set_visible(dialog, TRUE);

   } else {
      std::cout << "Congratulations: there are no bad chiral volumes in this molecule.\n";
      dialog = widget_from_builder("no_bad_chiral_volumes_dialog");
   }
   return dialog;
}

// static
void
graphics_info_t::on_inverted_chiral_volume_button_clicked (GtkButton       *button,
							   gpointer         user_data) {


   // This function may get called for an "unclick" too.  It may need
   // an active? test like the generic atom spec button callback.  But
   // perhaps not.

   graphics_info_t g;
   coot::atom_spec_t *atom_spec = (coot::atom_spec_t *) user_data;
//    std::cout << "atom_spec: "
// 	     << atom_spec->chain << " "
// 	     << atom_spec->resno << " "
// 	     << atom_spec->atom_name << " "
// 	     << atom_spec->alt_conf << " "
// 	     << std::endl;

   g.set_go_to_atom_molecule(atom_spec->int_user_data);
   g.set_go_to_atom_chain_residue_atom_name(atom_spec->chain_id.c_str(),
					    atom_spec->res_no,
					    atom_spec->atom_name.c_str(),
					    atom_spec->alt_conf.c_str());
   g.try_centre_from_new_go_to_atom();
   g.update_things_on_move_and_redraw();

}

// static
void
graphics_info_t::check_chiral_volume_molecule_combobox_changed(GtkWidget *w, gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(w));
   check_chiral_volume_molecule = imol;
}


void
graphics_info_t::fill_bond_parameters_internals(GtkWidget *combobox_for_molecule, int imol_active) {

   graphics_info_t g;

   // GtkWidget *bond_width_option_menu = lookup_widget(w, "bond_parameters_bond_width_optionmenu");
   GtkWidget *bond_width_combobox = widget_from_builder("bond_parameters_bond_width_combobox_text");
   gtk_combo_box_text_remove_all(GTK_COMBO_BOX_TEXT(bond_width_combobox));

   GtkWidget *draw_hydrogens_yes_radiobutton  = widget_from_builder("draw_hydrogens_yes_radiobutton");
   GtkWidget *draw_hydrogens_no_radiobutton   = widget_from_builder("draw_hydrogens_no_radiobutton");
   GtkWidget *draw_ncs_ghosts_yes_radiobutton = widget_from_builder("draw_ncs_ghosts_yes_radiobutton");
   GtkWidget *draw_ncs_ghosts_no_radiobutton  = widget_from_builder("draw_ncs_ghosts_no_radiobutton");

   bond_thickness_intermediate_value = -1;

   int idx_active = -1; // changed if the current bond width is the same as a comboxbox item
   int current_bond_width = 3;
   if (is_valid_model_molecule(imol_active))
      current_bond_width = molecules[imol_active].get_bond_thickness();

   for (int i=1; i<21; i++) {
      std::string s(int_to_string(i));
      gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(bond_width_combobox), nullptr, s.c_str());
      if (i == current_bond_width)
         idx_active = i-1; // not very elegant
   }
   if (idx_active >= 0)
      gtk_combo_box_set_active(GTK_COMBO_BOX(bond_width_combobox), idx_active);

   GCallback combobox_changed_func = G_CALLBACK(bond_parameters_bond_width_combobox_changed);
   g_signal_connect(bond_width_combobox, "changed", combobox_changed_func, NULL);

   // std::cout << "debug:: g_object set data on bond_width_combobox " << bond_width_combobox
   //              << " to  combobox_for_molecule " << combobox_for_molecule << std::endl;
   g_object_set_data(G_OBJECT(bond_width_combobox), "bond_parameters_molecule_combobox", combobox_for_molecule);


   // Draw Hydrogens?
   if (imol_active >= 0 ) {
      if (imol_active < n_molecules()) {
	 if (molecules[imol_active].has_model()) {
	    if (molecules[imol_active].draw_hydrogens()) {
	       gtk_check_button_set_active(GTK_CHECK_BUTTON(draw_hydrogens_yes_radiobutton), TRUE);
	    } else {
	       gtk_check_button_set_active(GTK_CHECK_BUTTON(draw_hydrogens_no_radiobutton), TRUE);
	    }
	 }
      }
   }


   // Draw NCS ghosts?
   if (imol_active >= 0 ) {
      if (imol_active < n_molecules()) {
	 if (molecules[imol_active].has_model()) {
	    if (molecules[imol_active].draw_ncs_ghosts_p()) {
	       if (molecules[imol_active].ncs_ghosts_have_rtops_p()) {
		  gtk_check_button_set_active(GTK_CHECK_BUTTON(draw_ncs_ghosts_yes_radiobutton), TRUE);
	       } else {
		  gtk_check_button_set_active(GTK_CHECK_BUTTON(draw_ncs_ghosts_no_radiobutton), TRUE);
	       }
	    } else {
	       gtk_check_button_set_active(GTK_CHECK_BUTTON(draw_ncs_ghosts_no_radiobutton), TRUE);
	    }
	 }
      }
   }

   // Make the frame be insensitive if there is no NCS.
   GtkWidget *frame = widget_from_builder("ncs_frame");
   short int make_insensitive = 1;
   if (imol_active >= 0 ) {
      if (imol_active < n_molecules()) {
	 if (molecules[imol_active].has_model()) {
	    if (molecules[imol_active].has_ncs_p()) {
	       make_insensitive = 0;
	    } else {
	       std::cout << "INFO:: in fill_bond_parameters_internals no NCS for  "
			 << imol_active << "\n";
	    }
	 } else {
	    std::cout << "ERROR:: bad imol in fill_bond_parameters_internals no model "
		      << imol_active << "\n";
	 }
      } else {
	 std::cout << "ERROR:: bad imol in fill_bond_parameters_internals i " << imol_active << "\n";
      }
   } else {
      std::cout << "ERROR:: bad imol in fill_bond_parameters_internals " << imol_active << "\n";
   }
   if (make_insensitive)
      gtk_widget_set_sensitive(frame, FALSE);
   else
      gtk_widget_set_sensitive(frame, TRUE);

}

// static
void
graphics_info_t::bond_parameters_bond_width_combobox_changed(GtkWidget *bond_width_combobox, gpointer data) {

   int id = gtk_combo_box_get_active(GTK_COMBO_BOX(bond_width_combobox));

   // we can't treat this as a comboboxtext if we have added numbers to it (not text)
   //
   // char *txt = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));

   // because of the way items were added. I don't know how to connect the active item
   // to an int in the general case (of a combobox). Frustration. Question for EB.
   //
   if (id >= 0) {
      graphics_info_t g;
      int bw = 1 + id;
      GtkWidget *molecule_combobox = GTK_WIDGET(g_object_get_data(G_OBJECT(bond_width_combobox), "bond_parameters_molecule_combobox"));
      std::cout << "debug:: g_object get data on bond_width_combobox " << bond_width_combobox
                << " for molecule_combobox " << molecule_combobox << std::endl;

      if (GTK_IS_COMBO_BOX(molecule_combobox)) {
         // std::cout << "debug:: " << molecule_combobox << " IS a combobox" << std::endl;
         int imol = g.combobox_get_imol(GTK_COMBO_BOX(molecule_combobox));
         // std::cout << "debug:: imol  from " << molecule_combobox << " is " << imol << std::endl;
         g.set_bond_thickness(imol, bw);
      } else {
         std::cout << "debug:: " << molecule_combobox << " is NOT a combobox" << std::endl;
      }
   }

}



void
graphics_info_t::fill_bond_colours_dialog_internal(GtkWidget *w) {

   // First the (global) step adjustment:
   GtkScale *hscale = GTK_SCALE(widget_from_builder("bond_parameters_colour_rotation_hscale"));
   std::cout << "in fill_bond_colours_dialog_internal() hscale " << hscale << std::endl;
   GtkAdjustment *adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(rotate_colour_map_on_read_pdb, 0.0, 370.0, 1.0, 20.0, 10.1));
   std::cout << "in fill_bond_colours_dialog_internal() adjustment " << adjustment << std::endl;
   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   g_signal_connect(G_OBJECT(adjustment), "value_changed",
                    G_CALLBACK(bond_parameters_colour_rotation_adjustment_changed), NULL);


   // Dead now.
   //
   // Now the "C only" checkbutton:
   // GtkWidget *checkbutton = widget_from_builder("bond_parameters_rotate_colour_map_c_only_checkbutton");
   // n   if (rotate_colour_map_on_read_pdb_c_only_flag)
   // gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);

   // 20220315-PE  but it doesn't work and who wants this mode anyway? Pink nitrogens?
   //
   // GtkWidget *frame_c_only = widget_from_builder("bond_parameters_rotate_colour_map_c_only_frame");
   // gtk_widget_set_visible(frame_c_only, FALSE);


   // Now the tricky bit, fill the scrolled vbox of molecule colour rotation step sliders:
   //
   GtkWidget *coords_colours_vbox = widget_from_builder("coords_colours_vbox");
   GtkWidget *label269;
   GtkWidget *label270;

   clear_out_container(coords_colours_vbox);

   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model()) {

	 std::string m = "Molecule ";
	 m += coot::util::int_to_string(imol);
	 m += " ";
	 m += molecules[imol].name_for_display_manager();
         GtkWidget *frame_molecule_N = gtk_frame_new(m.c_str());
         gtk_widget_set_hexpand(frame_molecule_N, TRUE);
         gtk_widget_set_margin_top(frame_molecule_N, 2);
         gtk_widget_set_margin_bottom(frame_molecule_N, 2);
         gtk_widget_set_margin_start(frame_molecule_N, 6);
         gtk_widget_set_margin_end(frame_molecule_N, 6);

	 // gtk_widget_ref (frame_molecule_N);
	 // g_object_set_data_full(G_OBJECT (coords_colour_control_dialog), "frame_molecule_N", frame_molecule_N, NULL);

	 gtk_box_append(GTK_BOX(coords_colours_vbox), frame_molecule_N);
	 gtk_widget_set_size_request(frame_molecule_N, 171, -1);

         GtkWidget *hbox136 = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
         gtk_widget_set_hexpand(hbox136, TRUE);

	 // gtk_widget_ref (hbox136);
	 // g_object_set_data_full (G_OBJECT (coords_colour_control_dialog), "hbox136", hbox136, NULL);

	 gtk_widget_set_visible (hbox136, TRUE);
	 gtk_frame_set_child(GTK_FRAME(frame_molecule_N), hbox136);

	 // g_object_set_data_full(G_OBJECT (coords_colour_control_dialog), "label269", label269, NULL);

	 label269 = gtk_label_new("    ");
	 gtk_box_append(GTK_BOX(hbox136), label269);

	 GtkAdjustment *adjustment_mol = GTK_ADJUSTMENT(gtk_adjustment_new(molecules[imol].bonds_colour_map_rotation,
                                                                           0.0, 370.0, 1.0, 20.0, 10.1));
         GtkWidget *coords_colour_hscale_mol_N = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, adjustment_mol);
         gtk_widget_set_hexpand(coords_colour_hscale_mol_N, TRUE);

	 gtk_range_set_adjustment(GTK_RANGE(coords_colour_hscale_mol_N), adjustment_mol);
	 g_signal_connect(G_OBJECT(adjustment_mol), "value_changed",
			  G_CALLBACK(bonds_colour_rotation_adjustment_changed), NULL);
	 g_object_set_data(G_OBJECT(adjustment_mol), "imol", GINT_TO_POINTER(imol));

	 gtk_widget_set_visible (coords_colour_hscale_mol_N, TRUE);

	 gtk_box_append(GTK_BOX (hbox136), coords_colour_hscale_mol_N);

	 label270 = gtk_label_new("  degrees  ");

	 gtk_widget_set_visible (label270, TRUE);
	 gtk_box_append(GTK_BOX(hbox136), label270);

         // needed?
         gtk_label_set_xalign(GTK_LABEL(label270), 0.5);
         gtk_label_set_yalign(GTK_LABEL(label270), 0.56);

	 gtk_widget_set_visible(frame_molecule_N, TRUE);
      }
   }

}

// static
void graphics_info_t::bond_parameters_colour_rotation_adjustment_changed(GtkAdjustment *adj,
									 GtkWidget *window) {

   graphics_info_t g;
   g.rotate_colour_map_on_read_pdb = gtk_adjustment_get_value(adj);
   graphics_draw(); // unnecessary.

}


// static
void graphics_info_t::bonds_colour_rotation_adjustment_changed(GtkAdjustment *adj,
							       GtkWidget *window) {

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(adj), "imol"));

   if (molecules[imol].has_model()) {
      float f =  gtk_adjustment_get_value(adj);
      molecules[imol].update_bonds_colour_using_map_rotation(f);
   }
   graphics_draw();

}

// // static
// void
// graphics_info_t::bond_width_item_select(GtkWidget *item, GtkPositionType pos) {

//    graphics_info_t g;
//    g.bond_thickness_intermediate_value = pos;
//    if (g.bond_thickness_intermediate_value > 0) {
//       int imol = g.bond_parameters_molecule;
//       if (is_valid_model_molecule(imol))
// 	 g.set_bond_thickness(imol, g.bond_thickness_intermediate_value);
//    }
// }


void
graphics_info_t::fill_add_OXT_dialog_internal(GtkWidget *dialog, int imol) {

   GtkWidget *chain_combobox = widget_from_builder("add_OXT_chain_combobox");
   GCallback signal_func = G_CALLBACK(add_OXT_chain_combobox_changed);

   std::string a = fill_combobox_with_chain_options(chain_combobox, imol, signal_func);
   if (a != "no-chain") {
      graphics_info_t::add_OXT_chain = a;
   }
}

// static
void
graphics_info_t::add_OXT_molecule_combobox_changed(GtkWidget *widget, gpointer data) {

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(widget));
   add_OXT_molecule = imol;

   std::cout << "DEBUG:: in add_OXT_molecule_combobox_changed() " << widget << " imol " << imol << std::endl;

   GtkWidget *chain_combobox = widget_from_builder("add_OXT_chain_combobox");
   GCallback func = G_CALLBACK(add_OXT_chain_combobox_changed);
   fill_combobox_with_chain_options(chain_combobox, imol, func);

}


// static
void
graphics_info_t::add_OXT_chain_combobox_changed(GtkWidget *combobox,
						gpointer data) {

   graphics_info_t g;
   add_OXT_chain = g.get_active_label_in_comboboxtext(GTK_COMBO_BOX_TEXT(combobox));

}


std::string
graphics_info_t::get_active_label_in_comboboxtext(GtkComboBoxText *combobox) {

   std::string s;
   gchar *at = gtk_combo_box_text_get_active_text(combobox);
   if (at)
      s = at;
   return s;
}

// get string for column 0 (which are strings)
std::string
graphics_info_t::get_active_label_in_combobox(GtkComboBox *combobox) const {

   std::string f_label;
   GtkTreeModel *model = gtk_combo_box_get_model(GTK_COMBO_BOX(combobox));
   GtkTreeIter iter;
   gboolean state = gtk_combo_box_get_active_iter(GTK_COMBO_BOX(combobox), &iter);
   if (state) {
      GValue f_label_as_value = { 0, };
      // g_value_init (&f_label_as_value, G_TYPE_STRING); init is done in the get below
      gtk_tree_model_get_value(model, &iter, 0, &f_label_as_value);
      const char *f_label_cstr = g_value_get_string(&f_label_as_value);
      f_label = f_label_cstr;
   } else {
      std::cout << "WARNING:: in get_active_label_in_combobox(): Bad state for get_active_iter"
		<< std::endl;
   }
   return f_label;
}



// static
std::string
graphics_info_t::fill_combobox_with_chain_options(GtkWidget *combobox_text,
						  int imol,
						  GCallback f) {

   return fill_combobox_with_chain_options(combobox_text, imol, f, "unset-chain");
}

// static
std::string
graphics_info_t::fill_combobox_with_chain_options(GtkWidget *combobox_text,
						  int imol,
						  GCallback f,
						  const std::string &acid) {

   // Note to self: Glade-2 generates comboboxes with gtk_combo_box_new_text()
   // This does not exist in gtk3.
   // gtk_combo_box_text_new() exists in gtk3, that creates a GtkComboBoxText, not a GtkComboBox
   //
   // All the GtkComboBoxes that I created with glade can be used as GtkComboBoxTexts.
   // Which will mean post-processing gtk2-interface.c to change
   // GtkComboBox to GtkComboBoxText
   //
   // So let treat the combobox that is passed to this function as a GtkComboBoxText.

   std::string r("no-chain");

   GtkComboBoxText *cb_as_text = GTK_COMBO_BOX_TEXT(combobox_text);

   gtk_combo_box_text_remove_all(cb_as_text);

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      std::vector<std::string> chains = coot::util::chains_in_molecule(mol);
      for (std::size_t i=0; i<chains.size(); i++) {
	 const std::string &ch_id = chains[i];
	 gtk_combo_box_text_append_text(cb_as_text, ch_id.c_str());
	 if ((i==0) || (ch_id == acid)) {
	    gtk_combo_box_set_active(GTK_COMBO_BOX(combobox_text), i);
	    r = ch_id;
	 }
      }
   }

   if (f)
      g_signal_connect(combobox_text, "changed", f, NULL);

   return r;
}



#if 0
// return the string at the top of the list:
//
//static
std::string graphics_info_t::fill_option_menu_with_chain_options(GtkWidget *option_menu,
								 int imol,
								 GCallback signal_func,
								 const std::string &active_chain_id) {

   std::pair<bool, std::string> top_string(false, "no-chain-set-in-fill_option_menu_with_chain_options");

   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> chains = coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
      GtkWidget *menu_item;
      GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
      if (menu)
	 gtk_widget_destroy(menu);
      menu = gtk_menu_new();
      gtk_widget_set_visible(menu, TRUE);
      for (unsigned int i=0; i<chains.size(); i++) {
	 menu_item = gtk_menu_item_new_with_label(chains[i].c_str());
	 g_signal_connect(G_OBJECT(menu_item), "activate", signal_func, NULL);
	 gtk_menu_append(GTK_MENU(menu), menu_item);
	 // g_object_set_data(G_OBJECT(menu_item), "chain-id", chains[i].c_str());
	 gtk_widget_set_visible(menu_item, TRUE);
	 if ((i == 0) || (chains[i] == active_chain_id)) {

	    if (i==0) {
	       top_string.first = true;
	       top_string.second = chains[i];
	    }
	    gtk_menu_set_active(GTK_MENU(menu), i);
	 }
      }
      /* Link the new menu to the optionmenu widget */
      gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
   }
   return top_string.second;
}
#endif


// static
void
graphics_info_t::fill_renumber_residue_range_dialog(GtkWidget *window) {

   // surely deletable?
#if 0
   graphics_info_t g;

   // GtkWidget *molecule_option_menu = lookup_widget(window, "renumber_residue_range_molecule_optionmenu");

   GtkWidget *molecule_combobox = widget_from_builder("renumber_residue_range_combobox");

   // GtkWidget *chain_option_menu = lookup_widget(window, "renumber_residue_range_chain_optionmenu");

   // renumber_residue_range_resno_1_entry
   // renumber_residue_range_resno_2_entry
   // renumber_residue_range_offset_entry

   // fill molecules option menu
   //    GtkSignalFunc callback_func =
   //       GTK_SIGNAL_FUNC(graphics_info_t::renumber_residue_range_molecule_menu_item_select);
   GCallback callback_func = G_CALLBACK(renumber_residue_range_molecule_combobox_changed);

   // g.fill_option_menu_with_coordinates_options(molecule_option_menu, callback_func);

   short int set_last_active_flag = 0;
   int imol_active = renumber_residue_range_molecule;

   // g.fill_option_menu_with_coordinates_options_internal_2(molecule_option_menu, callback_func,
   // set_last_active_flag, imol);

   g.fill_combobox_with_coordinates_options(molecule_combobox, callback_func, imol_active);

#endif
}

void
graphics_info_t::fill_renumber_residue_range_internal(GtkWidget *w, int imol) {

   GtkWidget *chain_combobox = widget_from_builder("renumber_residue_range_chain_combobox");
   GCallback callback_func = G_CALLBACK(renumber_residue_range_chain_combobox_changed);
   std::string a = fill_combobox_with_chain_options(chain_combobox, imol, callback_func);
   if (a != "no-chain") {
      graphics_info_t::renumber_residue_range_chain = a;
   }
}

void
graphics_info_t::renumber_residue_range_molecule_combobox_changed(GtkWidget *combobox,
								  gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
   renumber_residue_range_molecule = imol;
   GtkWidget *window = widget_from_builder("renumber_residue_range_dialog");
   g.fill_renumber_residue_range_internal(window, imol);

}

// static
void
graphics_info_t::renumber_residue_range_chain_combobox_changed(GtkWidget *combobox, gpointer data) {

   graphics_info_t g;
   std::string c = g.get_active_label_in_comboboxtext(GTK_COMBO_BOX_TEXT(combobox));
   g.renumber_residue_range_chain = c;
}

#include "coot-utils/peak-search.hh"

// static
// force_fill default false
void
graphics_info_t::fill_difference_map_peaks_button_box() {

   // does nothing if the diff map peaks dialog is not realized.

   auto make_map_parts_of_label = [] (unsigned int i_peak, const std::vector<std::pair<clipper::Coord_orth, float> > &centres,
                                      float map_sigma) {

      std::string label = "Peak ";
      float f = centres[i_peak].second/map_sigma;
      label += int_to_string(i_peak+1);
      label += ": ";
      label += float_to_string(centres[i_peak].second);
      label += " (";
      label += float_to_string(f);
      label += " rmsd) ";
#if 0
      label += "at (";
      label += coot::util::float_to_string_using_dec_pl(centres[i_peak].first.x(), 2);
      label += ", ";
      label += coot::util::float_to_string_using_dec_pl(centres[i_peak].first.y(), 2);
      label += ", ";
      label += coot::util::float_to_string_using_dec_pl(centres[i_peak].first.z(), 2);
      label += ")";
#endif
      // 20251215-PE instead of that, find the nearest atom to centres[i_peak]
      return label;
   };

   auto make_model_parts_of_label = [] (unsigned int i_peak,
                                        const std::vector<std::pair<clipper::Coord_orth, float> > &centres,
                                        int imol_coords) {

      std::string l;
      const clipper::Coord_orth &c = centres[i_peak].first;
      if (is_valid_model_molecule(imol_coords)) {
         coot::Cartesian p(c.x(), c.y(), c.z());
         coot::at_dist_info_t at_info = molecules[imol_coords].closest_atom(p, false);
         if (at_info.dist < 5.0) {
            if (at_info.atom) {
               l = " ";
               l += at_info.atom->GetChainID();
               l += " ";
               l += std::to_string(at_info.atom->GetSeqNum());
               l += " (";
               l += at_info.atom->GetResName();
               l += ") ";
               l += at_info.atom->GetAtomName();
               l += " d: ";
               l += coot::util::float_to_string_using_dec_pl(at_info.dist, 2);
            }
         }
      }
      return l;
   };

   auto fill_difference_map_button_box_inner = [make_map_parts_of_label, make_model_parts_of_label]
      (GtkWidget *button_vbox,
       const std::vector<std::pair<clipper::Coord_orth, float> > &centres,
       float map_sigma,
       int imol_coords) {

      clear_out_container(button_vbox);
      // a cutn'paste jobby from fill_rotamer_selection_buttons().
      GtkWidget *group = nullptr; // initially
      for (unsigned int i=0; i<centres.size(); i++) {
         std::string label = make_map_parts_of_label(i, centres, map_sigma);
         std::string m_lab = make_model_parts_of_label(i, centres, imol_coords);
         if (! m_lab.empty())
            label += m_lab;
         GtkWidget *radio_button = gtk_toggle_button_new_with_label(label.c_str());
         std::string button_name = "difference_map_peaks_button_";
         button_name += int_to_string(i);
         if (group)
            gtk_toggle_button_set_group(GTK_TOGGLE_BUTTON(radio_button), GTK_TOGGLE_BUTTON(group));
         else
            group = radio_button;

         gtk_widget_set_margin_bottom(radio_button, 4);
         gtk_widget_set_margin_top(radio_button, 4);
         gtk_widget_set_margin_start(radio_button, 6);
         gtk_widget_set_margin_end(radio_button, 6);

         coot::diff_map_peak_helper_data *hd = new coot::diff_map_peak_helper_data;
         hd->ipeak = i;
         hd->pos = centres[i].first;

         g_signal_connect(G_OBJECT (radio_button), "toggled",
                          G_CALLBACK(on_diff_map_peak_button_selection_toggled), hd);
         gtk_box_append(GTK_BOX(button_vbox), radio_button);

      }
   };

   auto make_diff_map_peaks = [] (GtkWidget *peaks_vbox) {

      bool do_positive_level_flag = true;
      bool do_negative_level_flag = true;
      bool around_model_only_flag = true;
      int imol_map    = -1;
      int imol_coords = -1;
      float n_sigma   = 5;
      do_positive_level_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(peaks_vbox), "do_positive_level_flag"));
      do_negative_level_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(peaks_vbox), "do_negative_level_flag"));
      around_model_only_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(peaks_vbox), "around_model_only_flag"));
      around_model_only_flag = false; // 20231228-PE this gets ignored currently
      char *n_sigma_cs = static_cast<char *>  (g_object_get_data(G_OBJECT(peaks_vbox), "n_sigma_str"));
      imol_map    = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(peaks_vbox), "imol_map"));
      imol_coords = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(peaks_vbox), "imol_model"));
      if (n_sigma_cs) {
         std::string n_sigma_str(n_sigma_cs);
         n_sigma = coot::util::string_to_float(n_sigma_str);
      }
      coot::peak_search ps(molecules[imol_map].xmap);
      ps.set_max_closeness(difference_map_peaks_max_closeness);
      std::vector<std::pair<clipper::Coord_orth, float> > centres;
      if (is_valid_model_molecule(imol_coords)) {
         if (is_valid_map_molecule(imol_map)) {
            centres = ps.get_peaks(molecules[imol_map].xmap,
                                    molecules[imol_coords].atom_sel.mol,
                                    n_sigma, do_positive_level_flag, do_negative_level_flag,
                                    around_model_only_flag);
         }
      }
      return centres;
   };


   GtkWidget *pane_to_show  = widget_from_builder("main_window_ramchandran_and_validation_pane");
   gtk_widget_set_visible(pane_to_show,  TRUE);

   GtkWidget *pane = widget_from_builder("main_window_graphics_rama_vs_graphics_pane");
   int pos = gtk_paned_get_position(GTK_PANED(pane));
   if (pos < 300)
      gtk_paned_set_position(GTK_PANED(pane), 380);

   GtkWidget *vbox78 = widget_from_builder("dialog-vbox78");
   gtk_widget_set_visible(vbox78,   TRUE);

   GtkWidget *outer_vbox = widget_from_builder("diff_map_peaks_outer_vbox");
   gtk_widget_set_visible(outer_vbox,   TRUE);

   GtkWidget *buttons_vbox = widget_from_builder("diff_map_peaks_vbox");

   GtkWidget *vboxes_vbox = widget_from_builder("validation_boxes_vbox");
   gtk_widget_set_visible(vboxes_vbox,   TRUE);

   std::vector<std::pair<clipper::Coord_orth, float> > centres = make_diff_map_peaks(buttons_vbox);
   // std::cout << "make_diff_map_peaks() made " << centres.size() << " centres" << std::endl;
   logger.log(log_t::INFO, logging::function_name_t("make_diff_map_peaks()"),
	      {"made", std::to_string(centres.size()), "centres"});
   float map_sigma = 0.5;
   int imol_map = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(buttons_vbox), "imol_map"));
   if (is_valid_map_molecule(imol_map))
      map_sigma = molecules[imol_map].map_sigma();
   int imol_coords = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(buttons_vbox), "imol_model"));
   fill_difference_map_button_box_inner(buttons_vbox, centres, map_sigma, imol_coords);
}


void
graphics_info_t::show_diff_map_peaks_vbox(int imol_map, int imol_coords,
                                                      const std::vector<std::pair<clipper::Coord_orth, float> > &centres_in,
                                                      float n_sigma,
                                                      bool do_positive_level_flag,
                                                      bool do_negative_level_flag,
                                                      bool around_model_only_flag) {

   std::vector<std::pair<clipper::Coord_orth, float> > centres = centres_in;

   // GtkWidget *w = create_diff_map_peaks_dialog();
   GtkWidget *peaks_vbox = widget_from_builder("diff_map_peaks_vbox");

   if (false) {
      GtkAllocation allocation;
      gtk_widget_get_allocation(peaks_vbox, &allocation);
      std::cout << "DEBUG:: peaks vbox allocation "
                << allocation.width << " " << allocation.height << std::endl;
   }

   char *n_sigma_cs = new char[20];
   std::string ns = std::to_string(n_sigma);
   unsigned int l = ns.length();
   strncpy(n_sigma_cs, ns.c_str(), l+1);

   g_object_set_data(G_OBJECT(peaks_vbox), "imol_map",               GINT_TO_POINTER(imol_map));
   g_object_set_data(G_OBJECT(peaks_vbox), "imol_model",             GINT_TO_POINTER(imol_coords));
   g_object_set_data(G_OBJECT(peaks_vbox), "do_positive_level_flag", GINT_TO_POINTER(do_positive_level_flag));
   g_object_set_data(G_OBJECT(peaks_vbox), "do_negative_level_flag", GINT_TO_POINTER(do_negative_level_flag));
   g_object_set_data(G_OBJECT(peaks_vbox), "around_model_only_flag", GINT_TO_POINTER(around_model_only_flag));
   g_object_set_data(G_OBJECT(peaks_vbox), "n_sigma_str", n_sigma_cs);

   // We don't need to do this if was already done (last time that the window was opened). How do I check
   // if there is already a signal attached to this button?  Why not just do it in glade?
   // g_signal_connect(G_OBJECT(update_button), "clicked", G_CALLBACK(diff_map_peaks_dialog_update_button_clicked_func), nullptr);

   // for . and , synthetic clicking.
   g_object_set_data(G_OBJECT(peaks_vbox), "centres_size", GINT_TO_POINTER(centres.size()));

   fill_difference_map_peaks_button_box(); // with position buttons

   // not used in the callback now that the button contains a pointer
   // to this info:
   diff_map_peaks->clear();
   for (unsigned int i=0; i<centres.size(); i++)
      diff_map_peaks->push_back(centres[i].first);

   max_diff_map_peaks = centres.size();

   if (! centres.empty()) {
      graphics_info_t g;
      coot::Cartesian c(centres[0].first.x(), centres[0].first.y(), centres[0].first.z());
      g.setRotationCentre(c, true);
      g.update_things_on_move();
      graphics_draw();
   }
}


// static
void
graphics_info_t::on_diff_map_peak_button_selection_toggled(GtkToggleButton  *button,
							   gpointer         user_data) {

   coot::diff_map_peak_helper_data *hd = static_cast<coot::diff_map_peak_helper_data *>(user_data);

   coot::Cartesian c(hd->pos.x(), hd->pos.y(), hd->pos.z());
   int button_state = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button));
   if (false)
      std::cout << "debug:: Here in on_diff_map_peak_button_selection_toggled() " << c << " "
                <<  button_state << std::endl;

   graphics_info_t g;

   if (button_state) {

      bool force_jump = true; // No slide, so that the following updates are done
                              // when we reach the new centre.
      bool have_jumped = g.setRotationCentre(c, force_jump);
      bool do_updates_now = have_jumped;
      if (do_updates_now) {
         g.update_things_on_move();
      }
      graphics_draw();
      std::string s = "Difference map peak number ";
      s += int_to_string(hd->ipeak);
      g.add_status_bar_text(s);
   }
   g.graphics_grab_focus();
}


// static
std::pair<short int, float>
graphics_info_t::float_from_entry(GtkWidget *entry) {

   std::pair<short int, float> p(0,0);
   const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(entry));
   if (txt) {
      float f = atof(txt);
      p.second = f;
      p.first = 1;
   }
   return p;
}

std::pair<short int, int>
graphics_info_t::int_from_entry(GtkWidget *entry) {

   std::pair<short int, int> p(0,0);
   const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(entry));
   if (txt) {
      int i = atoi(txt);
      p.second = i;
      p.first = 1;
   }
   return p;
}


// symmetry control dialog:
GtkWidget *graphics_info_t::wrapped_create_symmetry_controller_dialog() const {

   GtkWidget *w = symmetry_controller_dialog;
   if (! w) {
      //  w = create_symmetry_controller_dialog();
      w = widget_from_builder("symmetry_controller_dialog");
      symmetry_controller_dialog = w;
      for (int imol=0; imol<n_molecules(); imol++) {
	 if (molecules[imol].has_model())
	    molecules[imol].fill_symmetry_control_frame(w);
      }
   }

   //gtk_window_deiconify(GTK_WINDOW(w));
   return w;
}



GtkWidget *
graphics_info_t::wrapped_create_lsq_plane_dialog() {

   //GtkWidget *w = create_lsq_plane_dialog();
   GtkWidget *w = widget_from_builder("lsq_plane_dialog");
   pick_cursor_maybe();
   lsq_plane_dialog = w;
   GtkWindow *main_window = GTK_WINDOW(get_main_window());
   gtk_window_set_transient_for(GTK_WINDOW(w), main_window);

   return w;
}


//
GtkWidget *
wrapped_create_multi_residue_torsion_dialog(const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &pairs) {

   // GtkWidget *w = create_multi_residue_torsion_dialog();
   GtkWidget *w = widget_from_builder("multi_residue_torsion_dialog");
   GtkWidget *vbox = widget_from_builder("multi_residue_torsion_vbox");
   graphics_info_t::multi_residue_torsion_reverse_fragment_mode = 0; // reset every time

   for (unsigned int i=0; i<pairs.size(); i++) {
      std::string s;
      s += pairs[i].first->name;
      s += " ";
      s += coot::util::int_to_string(pairs[i].first->GetSeqNum());
      s += "  ->  ";
      s += pairs[i].second->name;
      s += " ";
      s += coot::util::int_to_string(pairs[i].second->GetSeqNum());
      GtkWidget *button = gtk_button_new_with_label(s.c_str());
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      gtk_box_append(GTK_BOX(vbox), button);
#else
      gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 2);
      gtk_container_set_border_width(GTK_CONTAINER(button), 2);
#endif
      g_signal_connect(G_OBJECT(button), "clicked",
		       G_CALLBACK(graphics_info_t::on_multi_residue_torsion_button_clicked),
		       GINT_TO_POINTER(i));
      gtk_widget_set_visible(button, TRUE);
      coot::atom_spec_t atom_spec_1(pairs[i].first);
      coot::atom_spec_t atom_spec_2(pairs[i].second);
      std::pair<coot::atom_spec_t, coot::atom_spec_t> *atom_spec_pair =
	 new std::pair<coot::atom_spec_t, coot::atom_spec_t>(atom_spec_1, atom_spec_2);
      g_object_set_data_full(G_OBJECT(button), "spec_pair", atom_spec_pair, g_free);
   }

   return w;
}

// static
void
graphics_info_t::on_multi_residue_torsion_button_clicked(GtkButton *button,
							 gpointer user_data) {

   graphics_info_t g;
   int i = GPOINTER_TO_INT(user_data);
   GtkWidget *check_button = widget_from_builder("multi_residue_torsion_reverse_checkbutton");
   std::pair<coot::atom_spec_t, coot::atom_spec_t> *atom_spec_pair =
      static_cast<std::pair<coot::atom_spec_t, coot::atom_spec_t> *> (g_object_get_data (G_OBJECT (button), "spec_pair"));
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_button)))
      g.multi_residue_torsion_reverse_fragment_mode = 1;
   else
      g.multi_residue_torsion_reverse_fragment_mode = 0;

   if (atom_spec_pair) {
      if (g.moving_atoms_asc->n_selected_atoms) {
	 if (moving_atoms_asc->mol) {
	    int index_1 = -1; // unset
	    int index_2 = -1; // unset
	    for (int ii=0; ii<g.moving_atoms_asc->n_selected_atoms; ii++) {
	       coot::atom_spec_t moving_spec_1(moving_atoms_asc->atom_selection[ii]);
	       if (moving_spec_1 == atom_spec_pair->first)
		  index_1 = ii;
	       if (moving_spec_1 == atom_spec_pair->second)
		  index_2 = ii;
	       if (index_1 != -1)
		  if (index_2 != -1)
		     break;
	    }

	    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_button)))
	       g.multi_residue_torsion_reverse_fragment_mode = 1;

	    if (index_1 == -1) {
	       std::cout << "ERROR:: index_1 not found " << std::endl;
	    } else {
	       if (index_2 == -1) {
		  std::cout << "ERROR:: index_2 not found " << std::endl;
	       } else {
		  std::pair<int, int> p(index_1, index_2);
		  g.multi_residue_torsion_rotating_atom_index_pair = p;
	       }
	    }
	 }
      }
   }
}



std::pair<double, double>
graphics_info_t::get_pointer_position_frac() const {

   double x = GetMouseBeginX();
   double y = GetMouseBeginY();

   GtkAllocation allocation;
   GtkWidget *glarea = glareas[0];
   gtk_widget_get_allocation(glarea, &allocation);

   double x_max = allocation.width;
   double y_max = allocation.height;

   double xf = x/x_max;
   double yf = y/y_max;
   return std::pair<double, double> (xf, yf);
}



int
graphics_info_t::add_molecular_representation(int imol,
                                              const std::string &atom_selection,
                                              const std::string &colour_scheme,
                                              const std::string &style,
                                              int secondary_structure_usage_flag) {

   std::cout << "g.add_molecular_representation(): atom_sel: \"" << atom_selection << "\" "
             << "colour-scheme: \"" << colour_scheme << "\" "
             << "style \"" << style << "\"" << std::endl;

   GtkWidget *w = widget_from_builder("molecular_representations_dialog");
   gtk_widget_set_visible(w, TRUE);
   set_transient_for_main_window(w);

   attach_buffers();

   int status = molecules[imol].add_molecular_representation(atom_selection, colour_scheme, style, secondary_structure_usage_flag);

   update_molecular_representation_widgets();
   graphics_draw();
   return status;
}

int
graphics_info_t::add_ribbon_representation_with_user_defined_colours(int imol, const std::string &name) {

   GtkWidget *w = widget_from_builder("molecular_representation_meshes_frame");
   gtk_widget_set_visible(w, TRUE);

   attach_buffers();

   int status = -1;
   molecules[imol].add_ribbon_representation_with_user_defined_residue_colours(user_defined_colours, name);

   update_molecular_representation_widgets();
   graphics_draw();
   return status;
}

void
graphics_info_t::remove_molecular_representation(int imol, int idx) {

   GtkWidget *w = widget_from_builder("main_window_meshes_frame");
   unsigned int n_mesh = 0;
   for (unsigned int i=0; i<molecules.size(); i++)
      n_mesh += molecules[i].meshes.size();

   if (n_mesh == 0)
      gtk_widget_set_visible(w, FALSE);

   molecules[imol].remove_molecular_representation(idx);
}



// -------- Meshes control (i.e. the Meshes of molecule_class_info)
void
graphics_info_t::set_show_molecular_representation(int imol, unsigned int mesh_idx, bool on_off) {

   GtkWidget *w = widget_from_builder("main_window_meshes_frame");
   gtk_widget_set_visible(w, TRUE);

   if (is_valid_model_molecule(imol)) {
      auto &meshes = molecules[imol].meshes;
      if (mesh_idx < meshes.size()) {
         auto mesh = meshes[mesh_idx];
         mesh.set_draw_mesh_state(on_off);
      }
   } 
}

// static
void
graphics_info_t::undisplay_all_molecule_meshes(int imol) {

   if (is_valid_model_molecule(imol)) {
      auto &m = molecules[imol];
      for (unsigned int mesh_idx=0; mesh_idx<m.meshes.size(); mesh_idx++) {
	 auto &mesh = m.meshes[mesh_idx];
	 mesh.set_draw_mesh_state(false);
      }
   }
}

// static
void
graphics_info_t::molecular_representation_meshes_checkbutton_toggled(GtkCheckButton *button, gpointer *user_data) {

   const char *n = static_cast<const char *>(g_object_get_data(G_OBJECT(button), "name"));
   if (n) {
      std::string sn(n);
      int imol     = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "imol"));
      int mesh_idx = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "mesh_idx"));
      if (is_valid_model_molecule(imol)) {
         auto &m = molecules[imol];
         if (mesh_idx < static_cast<int>(m.meshes.size())) {
            auto &mesh = m.meshes[mesh_idx];
            if (gtk_check_button_get_active(button)) {
               mesh.set_draw_mesh_state(true);
            } else {
               mesh.set_draw_mesh_state(false);
            }
         }
      } else {
         std::cout << "ERROR:: main_window_meshes_togglebutton_toggled() not a valid molecule" << std::endl;
      }
   }
   graphics_draw();
}

void
graphics_info_t::update_molecular_representation_widgets() {

   // find the display toggle button for mesh idx_mesh for molecule imol
   auto find_button = [] (GtkWidget *box, unsigned int imol, unsigned int idx_mesh) {
                         GtkWidget *w = 0;
                         std::string test_name = "molecular_representation_meshes_dialog_mesh_button_" +
                            std::to_string(imol) + "_" + std::to_string(idx_mesh);

                         GtkWidget *item_widget = gtk_widget_get_first_child(box);
                         while (item_widget) {
                            char *n = static_cast<char *>(g_object_get_data(G_OBJECT(item_widget), "name"));
                            if (n) {
                               std::string sn(n);
                               if (sn == test_name) {
                                  w = item_widget;
                                  break;
                               }
                            }
                            item_widget = gtk_widget_get_next_sibling(item_widget);
                         };
                         return w;
                      };

   GtkWidget *frame = widget_from_builder("molecular_representations_frame");
   GtkWidget *vbox  = widget_from_builder("molecular_representations_vbox");

   std::cout << "in update_molecular_representation_widgets(): vbox " << vbox << std::endl;

   unsigned int n_mesh = 0;
   for (unsigned int i=0; i<molecules.size(); i++)
      n_mesh += molecules[i].meshes.size();

   if (frame) {
      if (n_mesh == 0) {
         gtk_widget_set_visible(frame, FALSE);
      } else {
         gtk_widget_set_visible(frame, TRUE);
      }
   }

   // std::cout << "DEBUG:: update_molecular_representation_widgets() n_mesh " << n_mesh << std::endl;

   if (!vbox) return;

   for (unsigned int imol=0; imol<molecules.size(); imol++) {
      for (unsigned int j=0; j<molecules[imol].meshes.size(); j++) {
         const auto &mesh = molecules[imol].meshes[j];
         std::string m_name = mesh.name;
         std::string name = std::to_string(imol) + ": " + m_name;
         GtkWidget *w = find_button(vbox, imol, j);
         if (! w) {
            std::string label = name;
            w = gtk_check_button_new_with_label(label.c_str());
            gtk_widget_set_can_focus(w, FALSE);
            // gtk_widget_set_can_default(w, FALSE);
            std::string widget_name = "molecular_representation_meshes_dialog_mesh_button_" +
               std::to_string(imol) + "_" + std::to_string(j);
            unsigned int l = widget_name.length();
            char *widget_name_cstr = new char[l + 1];             // bleugh!
            strncpy(widget_name_cstr, widget_name.c_str(), l+1);  // bleugh!
            g_object_set_data(G_OBJECT(w), "name", widget_name_cstr);
            g_object_set_data(G_OBJECT(w), "imol",     GINT_TO_POINTER(imol));
            g_object_set_data(G_OBJECT(w), "mesh_idx", GINT_TO_POINTER(j));
            gtk_check_button_set_active(GTK_CHECK_BUTTON(w), TRUE);
            g_signal_connect(G_OBJECT(w), "toggled", G_CALLBACK(molecular_representation_meshes_checkbutton_toggled), nullptr);
            gtk_box_append(GTK_BOX(vbox), w);
            gtk_widget_set_visible(w, TRUE);
            gtk_widget_set_margin_start (w, 4);
            gtk_widget_set_margin_end   (w, 4);
            gtk_widget_set_margin_top   (w, 2);
            gtk_widget_set_margin_bottom(w, 2);
         }
      }
   }
}

void
graphics_info_t::fill_generic_validation_box_of_buttons(const std::string &dialog_label,
                                                        const std::vector<labelled_button_info_t> &v) {

    auto cb = +[] (GtkButton *button, gpointer user_data) {
       clipper::Coord_orth *co = reinterpret_cast<clipper::Coord_orth *>(user_data);
       set_rotation_centre(*co);
    };

    if (! v.empty()) {
       GtkWidget *box = widget_from_builder("generic_validation_box_of_buttons_box");
       if (box) {
          clear_out_container(box);
          for (unsigned int i = 0; i < v.size(); i++) {
             GtkWidget *box_for_item = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
             GtkWidget *button = gtk_button_new_with_label(v[i].label.c_str());
             gtk_widget_set_hexpand(button, TRUE);

             // I can't do this:
             // std::shared_ptr<clipper::Coord_orth> sco = std::make_shared<clipper::Coord_orth>(v[i].position);
             // void *user_data = reinterpret_cast<void *>(sco);
             // I should use a GObject?

             clipper::Coord_orth *co = new clipper::Coord_orth(v[i].position); // never deleted
             void *user_data = reinterpret_cast<void *>(co);
             g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(cb), user_data);
             gtk_box_append(GTK_BOX(box_for_item), button);
             gtk_box_append(GTK_BOX(box), box_for_item);
          }
       }
       GtkWidget *dialog = widget_from_builder("generic_validation_box_of_buttons_dialog");
       std::string title = std::string("Coot: ") + dialog_label;
       gtk_window_set_title(GTK_WINDOW(dialog), title.c_str());
       set_transient_for_main_window(dialog);
       gtk_window_present(GTK_WINDOW(dialog));
    }
}

void
graphics_info_t::fill_atoms_with_zero_occupancy_box_of_buttons(const std::vector<labelled_button_info_t> &lbv) {


    auto cb = +[] (GtkButton *button, gpointer user_data) {
       clipper::Coord_orth *co = reinterpret_cast<clipper::Coord_orth *>(user_data);
       set_rotation_centre(*co);
    };

    if (! lbv.empty()) {

       GtkWidget *frame = widget_from_builder("main_window_validation_graph_frame");
       gtk_widget_set_visible(frame, TRUE);

       GtkWidget *validation_graph_vbx = widget_from_builder("main_window_validation_graph_vbox");
       gtk_widget_set_visible(validation_graph_vbx, TRUE);

       GtkWidget *pane_to_show  = widget_from_builder("main_window_ramchandran_and_validation_pane");
       gtk_widget_set_visible(pane_to_show,  TRUE);

       GtkWidget *pane = widget_from_builder("main_window_graphics_rama_vs_graphics_pane");
       int pos = gtk_paned_get_position(GTK_PANED(pane));
       if (pos < 300)
	  gtk_paned_set_position(GTK_PANED(pane), 380);

       GtkWidget *validation_boxes_vbox = widget_from_builder("validation_boxes_vbox");
       gtk_widget_set_visible(validation_boxes_vbox, TRUE);

       GtkWidget *outer_vbox = widget_from_builder("dialog-vbox78");
       gtk_widget_set_visible(outer_vbox,   TRUE);

       GtkWidget *awzo_outer_vbox = widget_from_builder("atoms_with_zero_occupancy_outer_vbox");
       gtk_widget_set_visible(awzo_outer_vbox, TRUE);

       GtkWidget *box = widget_from_builder("atoms_with_zero_occupancy_vbox");
       if (box) {
	  clear_out_container(box);
	  for (unsigned int i = 0; i < lbv.size(); i++) {
	     GtkWidget *box_for_item = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 4);
	     GtkWidget *button = gtk_button_new_with_label(lbv[i].label.c_str());
	     gtk_widget_set_hexpand(button, TRUE);

	     clipper::Coord_orth *co = new clipper::Coord_orth(lbv[i].position); // never deleted
	     void *user_data = reinterpret_cast<void *>(co);
	     g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(cb), user_data);
	     gtk_box_append(GTK_BOX(box_for_item), button);
	     gtk_box_append(GTK_BOX(box), box_for_item);
	  }
       }
    }
}




void
graphics_info_t::add_shortcuts_to_window(GtkWidget *shortcuts_window) {

#if 0
   {
      auto shortcut_activated = +[] (GtkWidget *widget,
                                     GVariant  *unused,
                                     gpointer   row) {
         g_print ("activated %s\n", gtk_label_get_label(GTK_LABEL(row)));
         return gboolean(TRUE);
      };

      auto create_ctrl_g = +[] () {
         return gtk_keyval_trigger_new(GDK_KEY_g, GDK_CONTROL_MASK);
      };

      auto create_x = +[] () {
         return gtk_keyval_trigger_new(GDK_KEY_x, GdkModifierType(0));
      };

      struct {
         const char *description;
         GtkShortcutTrigger * (* create_trigger_func) (void);
      } shortcuts[] = {
         { "Press Ctrl-G", create_ctrl_g },
         { "Press X", create_x },
      };

      for (int i=0; i<2; i++) {
         GtkEventController *controller = gtk_shortcut_controller_new();
         GtkWidget *row = gtk_label_new (shortcuts[i].description);

         GtkShortcut *shortcut = gtk_shortcut_new(shortcuts[i].create_trigger_func(),
                                                  gtk_callback_action_new(shortcut_activated, row, NULL));

         gtk_shortcut_controller_add_shortcut(GTK_SHORTCUT_CONTROLLER(controller), shortcut);
      }

      // 4.14!
      // GtkShortcutsSection* section = ...
      // gtk_shortcuts_window_add_section(GTK_SHORTCUTS_WINDOW(shortcuts_window), section);
   }
#endif

}


// gui stuff = this is the function that is called when the slider changes
// g_object_get_data() is used to extract the imol
//
void
graphics_info_t::set_tomo_section_view_section(int imol, int section_index) {

   auto _ = [] (int err) {
      std::string s = std::to_string(err);
      if (err == GL_INVALID_ENUM)      s = "GL_INVALID_ENUM";
      if (err == GL_INVALID_OPERATION) s = "GL_INVALID_OPERATION";
      if (err == GL_INVALID_VALUE)     s = "GL_INVALID_VALUE";
      return s;
   };

   bool use_z_translation = true;
   bool do_X_and_Y_sections = true;

   if (is_valid_map_molecule(imol)) {

      // auto tp_start = std::chrono::high_resolution_clock::now();
      const auto &xmap = molecules[imol].xmap;
      int axis = 2; // for now
      clipper::Cell c_cell = xmap.cell();
      coot::Cell cell(c_cell.a(), c_cell.b(), c_cell.c(), c_cell.alpha(), c_cell.beta(), c_cell.gamma());
      clipper::Grid_sampling gs = xmap.grid_sampling();
      tomo_view_info = tomo_view_info_t(imol, cell, section_index, axis);

      float mean =    molecules[imol].map_mean();
      float std_dev = molecules[imol].map_sigma();
      float data_value_for_top    = mean + 3.5f * std_dev; // was  2.5
      float data_value_for_bottom = mean - 2.0f * std_dev; // was -1.5

      // data_value_for_bottom -= 3.0 * std_dev;

      if (false) {
         std::cout << "-------- current texture-meshes: " << std::endl;
         for (const auto &tm : texture_meshes)
            std::cout << "    " << tm.get_name() << std::endl;
      }

      // maybe I should replace the texture rather than delete all and create a new one.
      // texture_meshes.clear();

      auto eraser = [] (const TextureMesh &tm) {
         return (tm.get_name().find("-section") != std::string::npos);
      };

      texture_meshes.erase(std::remove_if(texture_meshes.begin(), texture_meshes.end(), eraser),
                           texture_meshes.end());

      // auto tp_1 = std::chrono::high_resolution_clock::now();
      // auto tp_2 = std::chrono::high_resolution_clock::now();


      if (true) {
         mini_texture_t m(xmap, section_index, axis, data_value_for_bottom, data_value_for_top); // 128 for 11729
         float x_len = m.x_size;
         float y_len = m.y_size;
         float z_pos = 0.0f;
         if (use_z_translation) z_pos = m.z_position;

         // std::cout << "Z-section x_len " << x_len << " y_len " << y_len << std::endl;

         attach_buffers();
         GLenum err = glGetError();
         if (err) std::cout << "GL ERROR:: tomo_section() A " << _(err) << "\n";
         Texture t(m, "mini-texture Z-section");
         err = glGetError();
         if (err) std::cout << "GL ERROR:: tomo_section() B " << _(err) << "\n";
         TextureInfoType ti(t, "mini-texture Z-section");
         err = glGetError();
         if (err) std::cout << "GL ERROR:: tomo_section() C " << _(err) << "\n";
         ti.unit = 0; // what is this?
         err = glGetError();
         if (err) std::cout << "GL ERROR:: tomo_section() D " << _(err) << "\n";
         TextureMesh tm("Tomo texture-mesh Z-section");
         tm.add_texture(ti);
         err = glGetError();
         if (err) std::cout << "GL ERROR:: tomo_section() E " << _(err) << "\n";
         tm.setup_tomo_quad(x_len, y_len, 0.0f, 0.0f, z_pos, false);
         // auto tp_2 = std::chrono::high_resolution_clock::now();
         err = glGetError();
         if (err) std::cout << "GL ERROR:: tomo_section() F " << _(err) << "\n";
         texture_meshes.push_back(tm);
         m.clear();
      }

      // auto tp_now = std::chrono::high_resolution_clock::now();
      // auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - tp_start);
      // auto td    = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1);
      // std::cout << "graphics_info_t::tomo_section() all  " << delta.count() << " ms" << std::endl;
      // std::cout << "graphics_info_t::tomo_section() test " <<    td.count() << " ms" << std::endl;

      if (do_X_and_Y_sections) {
         attach_buffers();
         // section index is now a mess.
         int section_index_X = gs.nu()/2;
         int section_index_Y = gs.nv()/2;
         section_index_X = section_index;
         section_index_Y = section_index;
         mini_texture_t m_x(xmap, section_index_X, 0, data_value_for_bottom, data_value_for_top);
         mini_texture_t m_y(xmap, section_index_Y, 1, data_value_for_bottom, data_value_for_top);

         Texture t_x(m_x, "mini-texture X-section");
         Texture t_y(m_y, "mini-texture Y-section");

         TextureMesh tm_x("Tomo texture-mesh X-section");
         TextureMesh tm_y("Tomo texture-mesh Y-section");
         TextureInfoType ti_x(t_x, "Tomo texture X section");
         TextureInfoType ti_y(t_y, "Tomo texture Y section");
         ti_x.unit = 0;
         ti_y.unit = 0;
         tm_x.add_texture(ti_x);
         tm_y.add_texture(ti_y);

         float offset_x_X_section = - c_cell.c() - c_cell.a() * 0.05f;
         float offset_y_X_section = 0.0f;
         float offset_x_Y_section = 0.0f;
         float offset_y_Y_section = c_cell.b() + c_cell.b() * 0.1f;

         tm_x.setup_tomo_quad(m_x.x_size, m_x.y_size, offset_x_X_section, offset_y_X_section, m_x.z_position, true);
         tm_y.setup_tomo_quad(m_y.x_size, m_y.y_size, offset_x_Y_section, offset_y_Y_section, m_y.z_position, false);

         texture_meshes.push_back(tm_x);
         texture_meshes.push_back(tm_y);
         m_x.clear();
         m_y.clear();

         GLenum err = glGetError();
         if (err) std::cout << "GL ERROR:: set_tomo_section_view_section() G " << _(err) << "\n";
      }

      graphics_draw();
   };
}
