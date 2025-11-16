/* src/glade-callbacks.cc
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2008 The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
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


#include <iostream>
#include <gtk/gtk.h>

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"
#include "positioned-widgets.h"
#include "interface.h"
#include "coot-references.h"

// put preferences functions into their own file, not here.
#include "coot-preferences.h"
#include "rotate-translate-modes.hh"
#include "read-phs.h"
#include "gtk-manual.h"
#include "c-interface-refine.h"
#include "widget-from-builder.hh"

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;

#include <vector>
#include "utils/coot-utils.hh"
#include "ideal/add-linked-cho.hh"
#include "graphics-info.h"

#include "cc-interface.hh" // for read_ccp4_map()

#include "utils/logging.hh"
extern logging logger;

// Let's put the new refinement and regularization control tools together here
// (although not strictly main window)


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refine_control_button_clicked(GtkButton       *button,
                                               gpointer         user_data) {
   graphics_info_t::show_refinement_and_regularization_parameters_frame();
}




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_select_map_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                           G_GNUC_UNUSED gpointer         user_data) {

   show_select_map_frame();
}

// <signal name="toggled" handler="on_model_toolbar_range_define_togglebutton_toggled" swapped="no"/>

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_range_define_togglebutton_toggled(GtkToggleButton *togglebutton,
                                                   G_GNUC_UNUSED gpointer user_data) {

   gboolean active = gtk_toggle_button_get_active(togglebutton);
   graphics_info_t g;
   if (active) {
      if (g.in_range_define == 0) g.in_range_define = 1;
   } else {
      if (g.in_range_define == 1) g.in_range_define = 0;
      if (g.in_range_define == 2) g.in_range_define = 0;
   }
   std::cout << "here now with active " << active << " in_range_define " << g.in_range_define << std::endl;

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_fixed_atoms_menubutton_activated(GtkMenuButton *button,
                                                  gpointer       user_data) {
   GtkWidget *w = wrapped_create_fixed_atom_dialog();
   gtk_widget_set_visible(w, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rigid_body_fit_togglebutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data) {

   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active) {
      printf("Rigid Body:\n");
      do_rigid_body_refine(1);
   } else {
      do_rigid_body_refine(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rot_trans_togglebutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active) {
      do_rot_trans_setup(1);
   } else {
      do_rot_trans_setup(0);
   }
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_auto_fit_rotamer_button_clicked(GtkButton *button,
                                                 gpointer   user_data) {

#if 0 // 20220813-PE there is no setup now
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active)
      setup_auto_fit_rotamer(1);
   else
      setup_auto_fit_rotamer(0);
#endif

   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      std::string alt_conf = aa.second->altLoc;
      coot::residue_spec_t res_spec(coot::atom_spec_t(aa.second));
      g.auto_fit_rotamer_ng(imol, res_spec, alt_conf);
   }

}


// 20220812-PE this used to be a toggle-button - it expected an atom pick after toggling the
// button. These days we find the rotamers for the residue at the centre of the screen.
// There is no longer any "setup_rotamers()"
//
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rotamers_button_clicked(GtkButton *button,
                                         gpointer  user_data) {

   // Find rotamers for the residue at the centre of the screen
   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      g.do_rotamers(imol, aa.second);
   }
}

// delete this when done
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_chi_angles_togglebutton_toggled(GtkToggleButton *toggletoolbutton,
                                                                          gpointer         user_data) {

   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active) {
      setup_edit_chi_angles(1);
   } else {
      setup_edit_chi_angles(0);
      set_show_chi_angle_bond(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_chi_angles_button_clicked(GtkButton *button,
                                                gpointer   user_data) {

   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   if (is_valid_model_molecule(imol)) {
      auto &m = g.molecules[imol];
      mmdb:: Atom *at = active_atom.second;
      int idx = m.get_atom_index(at);
      g.execute_edit_chi_angles(idx, imol);
   }
}



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_flip_peptide_button_clicked(GtkButton *button,
                                             gpointer   user_data) {
   graphics_info_t g;
   // we want the actual atom, not the CA of the residue
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = g.active_atom_spec_simple();
   if (aa.first) {
      int imol = aa.second.first;
      if (is_valid_model_molecule(imol)) {
         auto &m = g.molecules[imol];
         coot::atom_spec_t atom_spec(aa.second.second);
         m.pepflip(atom_spec);
         g.graphics_draw();
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_side_chain_180_button_clicked(GtkButton *button,
                                               gpointer   user_data) {

   // look at check_if_in_180_degree_flip_define, it checks for
   // intermediate atoms - we should do that here too.

   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   std::cout << "imol = " << imol << " atom = " << active_atom.second << std::endl;
   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *at = active_atom.second;
      std::string alt_conf(at->altLoc);
      coot::atom_spec_t at_spec(at);
      coot::residue_spec_t spec(at_spec);
      auto &m = g.molecules[imol];
      // change this signature to use a residue spec and an alt_conf.
      int istatus = m.do_180_degree_side_chain_flip(spec.chain_id, spec.res_no,
						    spec.ins_code, alt_conf, g.Geom_p());
      g.graphics_draw();
   }

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_mutate_and_autofit_menubutton_activate(GtkMenuButton *menubutton,
                                                        gpointer       user_data) {

   // this function seems not to be called on menu button click. Hmmm.
   std::cout << "on_model_toolbar_mutate_and_autofit_menubutton_active "
             << " select the right menu for mutate_and_autofit_menubutton here" << std::endl;
}


extern "C" G_MODULE_EXPORT
void
on_simple_mutate_menubutton_activate(GtkMenuButton *menubutton,
                                     gpointer       user_data) {

   // "activate" happens when the popover menu is shown. That is not what we want.

   std::cout << "on_simple_mutate_menubutton_activate "
             << " select the right menu for mutate_and_autofit_menubutton here" << std::endl;
}

#include "setup-gui-components.hh"

// try again with just a button

extern "C" G_MODULE_EXPORT
void
on_simple_mutate_button_clicked(GtkButton *button,
                                gpointer   user_data) {

   GtkWidget *menu_button = widget_from_builder("simple_mutate_menubutton");
   if (menu_button) {
      graphics_info_t g;
      auto active_atom = g.get_active_atom();
      int imol = active_atom.first;
      if (is_valid_model_molecule(imol)) {
         gtk_widget_set_visible(menu_button, TRUE);
         mmdb::Atom *atom = active_atom.second;
         mmdb::Residue *residue_p = atom->residue;
         if (residue_p) {
            if (coot::util::is_nucleotide_by_dict(residue_p, *g.Geom_p())) {
               add_typed_menu_to_mutate_menubutton("SIMPLE", "NUCLEIC-ACID");
            } else {
               add_typed_menu_to_mutate_menubutton("SIMPLE", "PROTEIN");
            }
            g_signal_emit_by_name(menu_button, "activate", NULL);
         }
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_terminal_residue_button_clicked(GtkButton *button,
                                                     gpointer   user_data)
{
   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = active_atom.second->residue;
      coot::atom_spec_t atom_spec(active_atom.second);
      std::string chain_id = atom_spec.chain_id;
      mmdb::Atom *atom = active_atom.second;
      std::string terminus_type = g.molecules[imol].get_term_type(atom);

      if (coot::util::is_nucleotide_by_dict_dynamic_add(residue_p, g.Geom_p())) {
         std::cout << "add nucleotide" << std::endl;
         g.execute_simple_nucleotide_addition(imol, terminus_type, residue_p, chain_id);
      } else {
         coot::residue_spec_t residue_spec(atom_spec);
         std::string new_type = "ALA";
         g.execute_add_terminal_residue(imol, terminus_type, residue_p, chain_id, new_type, true);
      }
   }

}




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_alt_conf_button_clicked(GtkButton *button,
                                             gpointer   user_data) {

   // altconf(); 20240304-PE this was this old "dialog before atom pick" method

   graphics_info_t g;
   auto active_atom_pair = g.get_active_atom();
   int imol = active_atom_pair.first;
   if (is_valid_model_molecule(imol)) {
      // split_residue_range(imol, add_alt_conf_atom_index, naii.atom_index);  // from graphics-info-defines
      // For today, we will make this just a single residue.
      const auto &active_atom = active_atom_pair.second;
      coot::atom_spec_t aas(active_atom);
      g.split_residue(imol, aas.chain_id, aas.res_no, aas.ins_code, aas.alt_conf);
   }
}




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_undo_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {
   apply_undo();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_redo_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {
  apply_redo();
}


extern "C" G_MODULE_EXPORT
void
on_save_coordinates_filechooser_dialog_response(GtkDialog       *dialog,
                                                gint             response_id,
                                                gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fnc = g_file_info_get_name(file_info);
      if (fnc) {

         // imol set in
         // on_save_coords_dialog_save_button_clicked(GtkButton       *button,
         //                                           gpointer         user_data)
         int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "imol"));
         save_coordinates(imol, fnc);
      }
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_dataset_filechooser_dialog_response(GtkDialog       *dialog,
                                       gint             response_id,
                                       gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      // const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fnc = g_file_info_get_name(file_info);
      if (fnc) {
         int is_auto = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "is_auto"));
         std::string fn(fnc);
         if (is_auto == 1) {
            auto_read_make_and_draw_maps_from_mtz(fnc);
         } else {
            manage_column_selector(fnc); // move the function declaration into a c++ header one day
         }
      }
      save_directory_from_filechooser(GTK_WIDGET(dialog));
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_map_name_filechooser_dialog_response(GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      // const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fnc = g_file_info_get_name(file_info);
      if (fnc) {
         std::string fn(fnc);
         // std::cout << "Now do something with " << fn << std::endl;
         bool is_diff_map = false;
         GtkWidget *checkbutton = widget_from_builder("map_filechooser_is_difference_map_button");
         // 20220809-PE GTK4 (post merge) - this is a checkbutton FIXME
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
            is_diff_map = true;
         std::cout << "----------------- read_ccp4_map() " << std::endl;
         read_ccp4_map(fn, is_diff_map);
         std::cout << "----------------- done read_ccp4_map() " << std::endl;
      }
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   
}


#include "rsr-functions.hh"


extern "C" G_MODULE_EXPORT
void
on_add_an_atom_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_water_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                    G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Water");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_lithium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Li");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_sodium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Na");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_potasium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                       G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("K");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_magnesium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                        G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Mg");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_calcium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Ca");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_strontium_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                        G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Sr");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_nickel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Ni");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_copper_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Cu");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_zinc_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Zn");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_chlorine_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                       G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Cl");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_bromine_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("Br");
   gtk_widget_set_visible(box, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_add_an_atom_iodine_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("add_an_atom_box");
   add_an_atom("I");
   gtk_widget_set_visible(box, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_get_monomer_ok_button_clicked(GtkButton       *button,
                                 gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("get_monomer_entry");
   if (entry) {
      handle_get_monomer_code(entry);
   }
   GtkWidget *frame = widget_from_builder("get_monomer_frame");
   gtk_widget_set_visible(frame, FALSE);
   graphics_info_t::graphics_grab_focus();
}


extern "C" G_MODULE_EXPORT
void on_generic_overlay_frame_cancel_button_clicked(GtkButton       *button,
                                                    gpointer         user_data) {
   GtkWidget* frame_widget = GTK_WIDGET(user_data);
   if(frame_widget) {
      gtk_widget_set_visible(frame_widget, FALSE);
   } else {
      g_error("ERROR:: in on_generic_overlay_frame_cancel_button_clicked() 'user_data' is NULL. Cannot hide overlay frame.");
   }
}

// 20240518-PE this is the only function in this file that uses python. Hmm.
// So rewrite mutate_by_overlap() into C++ (non-trivial)
//
#include "cc-interface-scripting.hh"

extern "C" G_MODULE_EXPORT
void on_replace_residue_ok_button_clicked(GtkButton *button,
                                          gpointer user_data) {

   GtkWidget *frame = widget_from_builder("replace_residue_frame");
   GtkWidget *entry = widget_from_builder("replace_residue_entry");
   std::string new_residue_type = gtk_editable_get_text(GTK_EDITABLE(entry));

   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = aa.second->residue;
      if (residue_p) {
	 std::string chain_id = residue_p->GetChainID();
	 int res_no = residue_p->GetSeqNum();
	 g.molecules[imol].mutate_by_overlap(chain_id, res_no, new_residue_type);
	 g.graphics_draw();
      }
   }
   gtk_widget_set_visible(frame, FALSE);

}

extern "C" G_MODULE_EXPORT
void on_replace_residue_cancel_button_clicked(GtkButton *button,
                                              gpointer user_data) {

   GtkWidget *frame = widget_from_builder("replace_residue_frame");
   gtk_widget_set_visible(frame, FALSE);

}

#include "c-interface-ligands.hh"

extern "C" G_MODULE_EXPORT
void
on_smiles_to_simple_3d_ok_button_clicked(GtkButton       *button,
                                         gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("smiles_to_simple_3d_entry");
   if (entry) {
      std::string smiles_string = gtk_editable_get_text(GTK_EDITABLE(entry));
      smiles_to_simple_3d(smiles_string);
   }
   GtkWidget *frame = widget_from_builder("smiles_to_simple_3d_frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_diff_map_peaks_close_button_clicked(GtkButton       *button,
                                       gpointer         user_data) {

   GtkWidget *vbox = widget_from_builder("diff_map_peaks_outer_vbox");
   clear_diff_map_peaks();
   gtk_widget_set_visible(vbox, FALSE);
   graphics_info_t::hide_vertical_validation_frame_if_appropriate();
}

extern "C" G_MODULE_EXPORT
void
on_diff_map_peaks_update_button_clicked(GtkButton *button,
                                       gpointer         user_data) {
   graphics_info_t g;
   g.fill_difference_map_peaks_button_box();

}

#include "dynamic-validation.hh"

extern "C" G_MODULE_EXPORT
void
on_dynamic_validation_update_button_clicked(GtkButton *button,
                                            gpointer  user_data) {

   update_dynamic_validation();
}

extern "C" G_MODULE_EXPORT
void
on_dynamic_validation_close_button_clicked(GtkButton *button,
                                           gpointer  user_data) {

   GtkWidget *validation_vbox = widget_from_builder("dynamic_validation_vbox");
   gtk_widget_set_visible(validation_vbox, FALSE);
   graphics_info_t::hide_vertical_validation_frame_if_appropriate();
}

extern "C" G_MODULE_EXPORT
void
on_dynamic_validation_include_missing_sidechains_checkbutton_toggled(GtkCheckButton *button,
                                                                     gpointer  user_data) {

   std::cout << "on_dynamic_validation_include_missing_sidechains_button_toggled() " << std::endl;

}

extern "C" G_MODULE_EXPORT
void
on_atoms_with_zero_occupancy_close_button_clicked(GtkButton *button, gpointer data) {

   GtkWidget *outer_vbox = widget_from_builder("atoms_with_zero_occupancy_outer_vbox");
   gtk_widget_set_visible(outer_vbox, FALSE);
   graphics_info_t g;
   g.graphics_grab_focus();
}


extern "C" G_MODULE_EXPORT
void
on_go_to_ligand_button_clicked(GtkButton *button,
                               gpointer   user_data) {
   go_to_ligand();
   graphics_info_t g;
   g.graphics_grab_focus();
}


extern "C" G_MODULE_EXPORT
void
on_graphics_grab_focus_button_clicked (GtkButton       *button,
                                       gpointer         user_data) {
   graphics_info_t g;
   g.graphics_grab_focus();

}

#include "cc-interface-graphics.hh"

extern "C" G_MODULE_EXPORT
void
on_coot_points_button_clicked(GtkButton       *button,
                              gpointer         user_data) {
   show_coot_points_frame();
   graphics_info_t g;
   g.graphics_grab_focus();
}



extern "C" G_MODULE_EXPORT
void
on_gaussian_surface_cancel_button_clicked(GtkButton       *button,
                                          gpointer         user_data) {
   GtkWidget *frame = widget_from_builder("gaussian_surface_frame");
   gtk_widget_set_visible(frame, FALSE);
   graphics_info_t g;
   g.graphics_grab_focus();
}

extern "C" G_MODULE_EXPORT
void
on_gaussian_surface_ok_button_clicked(GtkButton       *button,
                                      gpointer         user_data) {

   std::cout << "read the gui - make the surface" << std::endl;
   GtkWidget *frame = widget_from_builder("gaussian_surface_frame");
   GtkWidget *mol_chooser_combobox = widget_from_builder("gaussian_surface_molecule_chooser_combobox");
   GtkWidget *e_sigma          = widget_from_builder("gaussian_surface_sigma_entry");
   GtkWidget *e_radius         = widget_from_builder("gaussian_surface_radius_entry");
   GtkWidget *e_contour_level  = widget_from_builder("gaussian_surface_contour_level_entry");
   GtkWidget *e_b_factor       = widget_from_builder("gaussian_surface_b_factor_entry");
   GtkWidget *e_chain_col_mode = widget_from_builder("gaussian_surface_chain_colour_entry");

   int imol = my_combobox_get_imol(GTK_COMBO_BOX(mol_chooser_combobox));

   try {
      float sigma  = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_sigma)));
      float radius = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_radius)));
      float cl     = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_contour_level)));
      float bf     = coot::util::string_to_float(gtk_editable_get_text(GTK_EDITABLE(e_b_factor)));
      int cc_mode  = coot::util::string_to_int(gtk_editable_get_text(GTK_EDITABLE(e_chain_col_mode)));
      set_gaussian_surface_sigma(sigma);
      set_gaussian_surface_box_radius(radius);
      set_gaussian_surface_contour_level(cl);
      set_gaussian_surface_fft_b_factor(bf);
      set_gaussian_surface_chain_colour_mode(cc_mode);
      gaussian_surface(imol);
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

   gtk_widget_set_visible(frame, FALSE);

}

void
fill_comboboxtext_with_atom_of_residue_type(const char *rn, GtkWidget *comboboxtext) {

   std::cout << "fill_comboboxtext_with_atom_of_residue_type() --- start --- " << std::endl;

   if (rn) {
      int imol = 0;
      std::string residue_type(rn);
      auto &geom = *graphics_info_t::Geom_p();
      bool state = geom.have_dictionary_for_residue_type(residue_type, imol, graphics_info_t::cif_dictionary_read_number, true);
      graphics_info_t::cif_dictionary_read_number++;
      std::pair<bool, coot::dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(residue_type, imol);
      if (rp.first) {
         gtk_combo_box_text_remove_all(GTK_COMBO_BOX_TEXT(comboboxtext));
         const auto &restraints = rp.second;
         const auto &atoms = restraints.atom_info;
         for (const auto &atom : atoms) {
            if (atom.type_symbol != "H") {
               std::cout << "  type " << rn << " fill comboboxtext " << comboboxtext << " with " << atom.atom_id << std::endl;
               gtk_combo_box_text_append(GTK_COMBO_BOX_TEXT(comboboxtext), atom.atom_id.c_str(), atom.atom_id.c_str());
            }
         }
      } else {
         std::cout << "ERROR get_monomer_restraints() failed for type " << residue_type << std::endl;
      }
   } else {
      std::cout << "ERROR:: null rn in fill_comboboxtext_with_atom_of_residue_type()" << std::endl;
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_delete_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
						     gpointer         user_data) {

   std::cout << "delete atom first toggled" << std::endl;
   GtkWidget *delete_atom_combobox = widget_from_builder("acedrg_link_first_delete_atom_chooser_combobox");
   GtkWidget *entry                = widget_from_builder("acedrg_link_first_residue_name_entry");
   if (delete_atom_combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(delete_atom_combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, delete_atom_combobox);
            } else {
               gtk_widget_set_sensitive(delete_atom_combobox, FALSE);
            }
         }
      }
   }

}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_change_bond_order_checkbutton_toggled(GtkCheckButton *checkbutton,
							   gpointer         user_data) {

   std::cout << "cbo first toggled" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_first_change_bond_order_combobox");
   GtkWidget *combobox_atom_1 = widget_from_builder("acedrg_link_first_change_bond_order_atom_1_combobox");
   GtkWidget *combobox_atom_2 = widget_from_builder("acedrg_link_first_change_bond_order_atom_2_combobox");
   GtkWidget *entry           = widget_from_builder("acedrg_link_first_residue_name_entry");
   if (checkbutton) {
      if (gtk_check_button_get_active(checkbutton)) {
         gtk_widget_set_sensitive(combobox, TRUE);
         gtk_widget_set_sensitive(combobox_atom_1, TRUE);
         gtk_widget_set_sensitive(combobox_atom_2, TRUE);
         if (entry) {
            const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
            if (t) {
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_1);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_2);
            }
         }
      } else {
         gtk_widget_set_sensitive(combobox, FALSE);
         gtk_widget_set_sensitive(combobox_atom_1, FALSE);
         gtk_widget_set_sensitive(combobox_atom_2, FALSE);
      }
   }

}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_delete_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
						     gpointer         user_data) {

   std::cout << "delete atom second toggled" << std::endl;
   GtkWidget *delete_atom_combobox = widget_from_builder("acedrg_link_second_delete_atom_chooser_combobox");
   GtkWidget *entry                = widget_from_builder("acedrg_link_second_residue_name_entry");
   if (delete_atom_combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(delete_atom_combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, delete_atom_combobox);
            } else {
               gtk_widget_set_sensitive(delete_atom_combobox, FALSE);
            }
         }
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_change_bond_order_checkbutton_toggled(GtkCheckButton *checkbutton,
							   gpointer         user_data) {

   std::cout << "cbo toggled second" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_second_change_bond_order_combobox");
   GtkWidget *combobox_atom_1 = widget_from_builder("acedrg_link_second_change_bond_order_atom_1_combobox");
   GtkWidget *combobox_atom_2 = widget_from_builder("acedrg_link_second_change_bond_order_atom_2_combobox");
   GtkWidget *entry           = widget_from_builder("acedrg_link_second_residue_name_entry");
   if (checkbutton) {
      if (gtk_check_button_get_active(checkbutton)) {
         gtk_widget_set_sensitive(combobox, TRUE);
         gtk_widget_set_sensitive(combobox_atom_1, TRUE);
         gtk_widget_set_sensitive(combobox_atom_2, TRUE);
         if (entry) {
            const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
            if (t) {
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_1);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox_atom_2);
            }
         }
      } else {
         gtk_widget_set_sensitive(combobox, FALSE);
         gtk_widget_set_sensitive(combobox_atom_1, FALSE);
         gtk_widget_set_sensitive(combobox_atom_2, FALSE);
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_change_charge_on_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                               gpointer         user_data) {

   std::cout << "change charge atom first toggled" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_first_change_charge_on_atom_chooser_combobox");
   GtkWidget *entry    = widget_from_builder("acedrg_link_first_residue_name_entry");
   if (combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox);
            } else {
               gtk_widget_set_sensitive(combobox, FALSE);
            }
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_change_charge_on_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                                gpointer         user_data) {

   std::cout << "change charge atom second toggled" << std::endl;
   GtkWidget *combobox = widget_from_builder("acedrg_link_second_change_charge_on_atom_chooser_combobox");
   GtkWidget *entry    = widget_from_builder("acedrg_link_second_residue_name_entry");
   if (combobox) {
      if (entry) {
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            if (gtk_check_button_get_active(checkbutton)) {
               gtk_widget_set_sensitive(combobox, TRUE);
               fill_comboboxtext_with_atom_of_residue_type(t, combobox);
            } else {
               gtk_widget_set_sensitive(combobox, FALSE);
            }
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_ok_button_clicked(GtkButton       *button,
                                 gpointer         user_data) {

   auto simple_link = [] (const std::string &residue_name_first,
                          const std::string &residue_name_second,
                          const std::string &atom_name_first,
                          const std::string &atom_name_second,
                          const std::string &cif_file_name_1,
                          const std::string &cif_file_name_2,
                          const std::string &bond_order) {

      std::string ss = "LINK: ";

      ss += "RES-NAME-1 ";
      ss += residue_name_first;
      ss += " ";
      ss += "ATOM-NAME-1 ";
      ss += atom_name_first;
      ss += " ";
      if (!cif_file_name_1.empty())
         ss += std::string("FILE-1 ") + cif_file_name_1;

      ss += "RES-NAME-2 ";
      ss += residue_name_second;
      ss += " ";
      ss += "ATOM-NAME-2 ";
      ss += atom_name_second;
      ss += " ";
      if (!cif_file_name_2.empty())
         ss += std::string("FILE-2 ") + cif_file_name_2;

      ss += std::string(" BOND-TYPE ");
      ss += coot::util::upcase(bond_order);
      std::cout << ss << std::endl;
      run_acedrg_link_generation(ss);
   };

   // Here's an example:
   //
   // LINK: RES-NAME-1 LYS ATOM-NAME-1 NZ RES-NAME-2 PLP ATOM-NAME-2 C4A BOND-TYPE DOUBLE DELETE ATOM O 1 DELETE ATOM O4A 2 CHANGE BOND C OXT double 1 CHANGE BOND C4 C4A triple 2


   auto link = [] (const std::string &residue_name_first, const std::string &residue_name_second,
                   const std::string &atom_name_first, const std::string &atom_name_second,
                   const std::string &cif_file_name_1, const std::string &cif_file_name_2,
                   const std::string &bond_order,
                   bool delete_atom_first, const char *da_first,
                   bool delete_atom_second, const char *da_second,
                   bool change_bond_order_first, const char *cbo_first,
                   const char *change_bond_order_first_atom_1,  const char *change_bond_order_first_atom_2,
                   bool change_bond_order_second, const char *cbo_second,
                   const char *change_bond_order_second_atom_1, const char *change_bond_order_second_atom_2,
                   bool change_charge_on_first_residue_atom,  const char *change_charge_on_first_atom,
                   bool change_charge_on_second_residue_atom, const char *change_charge_on_second_atom) {

      std::string ss = "LINK: ";
      ss += "RES-NAME-1 ";
      ss += residue_name_first;
      ss += " ";
      ss += "ATOM-NAME-1 ";
      ss += atom_name_first;
      ss += " ";
      if (!cif_file_name_1.empty())
         ss += "FILE-1 " + cif_file_name_1 + std::string(" ");
      if (delete_atom_first)
         if (da_first)
         ss += "DELETE ATOM " + std::string(da_first) + std::string(" 1 ");
      if (change_charge_on_first_residue_atom)
         ss += std::string("CHANGE CHARGE ") + std::string(change_charge_on_first_atom) + std::string(" 1 ");
      if (change_bond_order_first)
         if (cbo_first)
            if (change_bond_order_first_atom_1)
               if (change_bond_order_first_atom_1)
                  ss += std::string("CHANGE BOND ") + std::string(change_bond_order_first_atom_1) + std::string(" ") +
                     std::string(change_bond_order_first_atom_2) + std::string(" ") + std::string(cbo_first) + " 1 ";

      ss += "RES-NAME-2 ";
      ss += residue_name_second;
      ss += " ";
      ss += "ATOM-NAME-2 ";
      ss += " ";
      ss += atom_name_second;
      ss += " ";
      if (!cif_file_name_2.empty())
         ss += std::string("FILE-2 ") + cif_file_name_2 + std::string(" ");
      if (delete_atom_second)
         if (da_second)
         ss += "DELETE ATOM " + std::string(da_second) + std::string(" 2 ");
      if (change_charge_on_second_residue_atom)
         ss += std::string("CHANGE CHARGE ") + std::string(change_charge_on_second_atom) + std::string(" 2 ");
      if (change_bond_order_second)
         if (cbo_second)
            if (change_bond_order_second_atom_1)
               if (change_bond_order_second_atom_1)
                  ss += std::string("CHANGE BOND ") + std::string(change_bond_order_second_atom_1) + std::string(" ") +
                  change_bond_order_second_atom_2 + std::string(" ") + std::string(cbo_second) + " 2 ";

      ss += std::string(" BOND-TYPE ");
      ss += coot::util::upcase(bond_order);
      run_acedrg_link_generation(ss);
      std::cout << ss << std::endl;
   };

   auto get_cif_file_name = [] (const std::string &residue_type) {

      std::string file_name;
      auto &geom = *graphics_info_t::Geom_p();
      int imol = 0;
      bool state = geom.have_dictionary_for_residue_type(residue_type, imol, graphics_info_t::cif_dictionary_read_number, true);
      graphics_info_t::cif_dictionary_read_number++;
      std::pair<bool, coot::dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(residue_type, imol);
      if (rp.first) {
          file_name = rp.second.cif_file_name;
      }
      return file_name;
   };

   GtkWidget *w = widget_from_builder("acedrg_link_interface_frame");
   gtk_widget_set_visible(w, FALSE);

   GtkWidget *bond_order_combobox = widget_from_builder("acedrg_link_bond_order_combobox");

   GtkWidget *entry_first                         = widget_from_builder("acedrg_link_first_residue_name_entry");
   GtkWidget *atom_name_combobox_first            = widget_from_builder("acedrg_link_first_atom_name_chooser_combobox");

   GtkWidget *entry_second                         = widget_from_builder("acedrg_link_second_residue_name_entry");
   GtkWidget *atom_name_combobox_second            = widget_from_builder("acedrg_link_second_atom_name_chooser_combobox");

   // delete atom
   GtkWidget *delete_atom_checkbutton_second       = widget_from_builder("acedrg_link_second_delete_atom_checkbutton");
   GtkWidget *delete_atom_combobox_second          = widget_from_builder("acedrg_link_second_delete_atom_chooser_combobox");
   GtkWidget *delete_atom_checkbutton_first       = widget_from_builder("acedrg_link_first_delete_atom_checkbutton");
   GtkWidget *delete_atom_combobox_first          = widget_from_builder("acedrg_link_first_delete_atom_chooser_combobox");

   // change bond order
   GtkWidget *change_bond_order_checkbutton_first = widget_from_builder("acedrg_link_first_change_bond_order_checkbutton");
   GtkWidget *change_bond_order_combobox_first    = widget_from_builder("acedrg_link_first_change_bond_order_combobox");
   GtkWidget *change_bond_order_checkbutton_second = widget_from_builder("acedrg_link_second_change_bond_order_checkbutton");
   GtkWidget *change_bond_order_combobox_second    = widget_from_builder("acedrg_link_second_change_bond_order_combobox");
   GtkWidget *cbcbof1                              = widget_from_builder("acedrg_link_first_change_bond_order_atom_1_combobox");
   GtkWidget *cbcbof2                              = widget_from_builder("acedrg_link_first_change_bond_order_atom_2_combobox");
   GtkWidget *cbcbos1                              = widget_from_builder("acedrg_link_second_change_bond_order_atom_1_combobox");
   GtkWidget *cbcbos2                              = widget_from_builder("acedrg_link_second_change_bond_order_atom_2_combobox");

   // change charge
   GtkWidget *cc_combobox_first  = widget_from_builder("acedrg_link_first_change_charge_on_atom_chooser_combobox");
   GtkWidget *cc_combobox_second = widget_from_builder("acedrg_link_second_change_charge_on_atom_chooser_combobox");
   GtkWidget *change_charge_first_checkbutton  = widget_from_builder("on_acedrg_link_first_change_charge_on_atom_checkbutton");
   GtkWidget *change_charge_second_checkbutton = widget_from_builder("on_acedrg_link_second_change_charge_on_atom_checkbutton");

   // need to add a pair of comboboxes for change bond order atom names for both first and second.

   if (!bond_order_combobox) return;

   if (entry_first && atom_name_combobox_first && delete_atom_combobox_first && delete_atom_checkbutton_first) {
      std::cout << "Here A " << change_bond_order_combobox_first << " " << change_bond_order_checkbutton_first << std::endl;
      if (change_bond_order_combobox_first && change_bond_order_checkbutton_first) {
         std::cout << "Here B " << std::endl;
         if (entry_second && atom_name_combobox_second && delete_atom_combobox_second && delete_atom_checkbutton_second) {
            std::cout << "Here C " << std::endl;
         if (change_bond_order_combobox_second && change_bond_order_checkbutton_second) {
               std::cout << "Here D " << atom_name_combobox_first << " " << atom_name_combobox_second << std::endl;
               std::cout << "Here D first  " << GTK_IS_COMBO_BOX_TEXT(atom_name_combobox_first)  << std::endl;
               std::cout << "Here D second " << GTK_IS_COMBO_BOX_TEXT(atom_name_combobox_second) << std::endl;
               char *bond_order = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(bond_order_combobox));
               if (bond_order) {
                  std::cout << "Here D 0" << std::endl;
                  char *atom_name_first  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(atom_name_combobox_first));
                  std::cout << "Here D 1" << std::endl;
                  char *atom_name_second = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(atom_name_combobox_second));
                  std::cout << "Here D 2" << std::endl;
                  if (atom_name_first  == 0) std::cout << "Here D with atom_name_first null"  << std::endl;
                  if (atom_name_second == 0) std::cout << "Here D with atom_name_second null" << std::endl;
                  std::cout << "Here Done D" << std::endl;
                  if (atom_name_first && atom_name_second) {
                     std::cout << "Here E atom_name_first " << atom_name_first << std::endl;
                     std::cout << "Here E atom_name_second " << atom_name_second << std::endl;
                     bool change_bond_order_first  = false;
                     bool change_bond_order_second = false;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_bond_order_checkbutton_first)))  change_bond_order_first  = true;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_bond_order_checkbutton_second))) change_bond_order_second = true;
                     bool delete_atom_first  = false;
                     bool delete_atom_second = false;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(delete_atom_checkbutton_first)))  delete_atom_first  = true;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(delete_atom_checkbutton_second))) delete_atom_second = true;
                     const char *residue_name_first  = gtk_editable_get_text(GTK_EDITABLE(entry_first));
                     const char *residue_name_second = gtk_editable_get_text(GTK_EDITABLE(entry_second));
                     bool change_charge_first  = false;
                     bool change_charge_second = false;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_charge_first_checkbutton)))  change_charge_first  = true;
                     if (gtk_check_button_get_active(GTK_CHECK_BUTTON(change_charge_second_checkbutton))) change_charge_second = true;
                     if (residue_name_first) {
                        if (residue_name_second) {
                           std::string cif_file_name_1 = get_cif_file_name(residue_name_first);
                           std::string cif_file_name_2 = get_cif_file_name(residue_name_second);
                           char *da_first   = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(delete_atom_combobox_first));
                           char *da_second  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(delete_atom_combobox_second));
                           char *cbo_first  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(change_bond_order_combobox_first));
                           char *cbo_second = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(change_bond_order_combobox_second));
                           std::cout << "Here F " << change_bond_order_first << " " << change_bond_order_second << " "
                                     << delete_atom_first << " " << delete_atom_second << std::endl;
                           if (atom_name_first && atom_name_second) {
                              if (!change_bond_order_first && !change_bond_order_second && !delete_atom_first && !delete_atom_second) {
                                 simple_link(residue_name_first, residue_name_second, atom_name_first, atom_name_second,
                                             cif_file_name_1, cif_file_name_2, std::string(bond_order));
                              } else {
                                 if (cbcbof1 && cbcbof2 && cbcbos1 && cbcbos2) {
                                    const char *change_bond_order_first_atom_1  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbof1));
                                    const char *change_bond_order_first_atom_2  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbof2));
                                    const char *change_bond_order_second_atom_1 = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbos1));
                                    const char *change_bond_order_second_atom_2 = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cbcbos2));
                                    const char *cc_first_atom                   = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cc_combobox_first));
                                    const char *cc_second_atom                  = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(cc_combobox_second));
                                    link(residue_name_first, residue_name_second, atom_name_first, atom_name_second,
                                         cif_file_name_1, cif_file_name_2, std::string(bond_order),
                                         delete_atom_first, da_first,
                                         delete_atom_second, da_second,
                                         change_bond_order_first, cbo_first,   change_bond_order_first_atom_1,  change_bond_order_first_atom_2,
                                         change_bond_order_second, cbo_second, change_bond_order_second_atom_1, change_bond_order_second_atom_2,
                                         change_charge_first, cc_first_atom,
                                         change_charge_second, cc_second_atom);
                                 } else {
                                    std::cout << "combobox cbo lookup failure" << std::endl;
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     std::cout << "WARNING:: BAD input: atom_name_first && atom_name_second failed" << std::endl;
                  }
               }
            }
         }
      }
   }
   graphics_info_t g;
   g.graphics_grab_focus();
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_first_residue_activate(GtkEntry *entry, gpointer user_data) {

   GtkWidget *atom_name_combobox_first = widget_from_builder("acedrg_link_first_atom_name_chooser_combobox");
   std::cout << "Fill this: " << atom_name_combobox_first << std::endl;
   const char *rn = gtk_editable_get_text(GTK_EDITABLE(entry));
   gtk_widget_set_sensitive(atom_name_combobox_first, TRUE);
   if (rn) {
      fill_comboboxtext_with_atom_of_residue_type(rn, atom_name_combobox_first);
   } else {
      std::cout << "No residue name first " << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_acedrg_link_second_residue_activate(GtkEntry *entry, gpointer user_data) {

   GtkWidget *atom_name_combobox_second = widget_from_builder("acedrg_link_second_atom_name_chooser_combobox");
   std::cout << "Fill this: " << atom_name_combobox_second << std::endl;
   const char *rn = gtk_editable_get_text(GTK_EDITABLE(entry));
   gtk_widget_set_sensitive(atom_name_combobox_second, TRUE);
   if (rn) {
      fill_comboboxtext_with_atom_of_residue_type(rn, atom_name_combobox_second);
   } else {
      std::cout << "No residue name second" << std::endl;
   }
}

extern "C" G_MODULE_EXPORT
void
on_acedrg_link_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
				     G_GNUC_UNUSED gpointer         user_data) {

   std::cout << "Cancel" << std::endl;
   GtkWidget *w = widget_from_builder("acedrg_link_interface_frame");
   gtk_widget_set_visible(w, FALSE);
   graphics_info_t g;
   g.graphics_grab_focus();
}

extern "C" G_MODULE_EXPORT
void
on_flip_hand_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("flip_map_hand_frame");
   gtk_widget_set_visible(frame, FALSE);
   graphics_info_t g;
   g.graphics_grab_focus();

}

extern "C" G_MODULE_EXPORT
void
on_flip_hand_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                               G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame                = widget_from_builder("flip_map_hand_frame");
   GtkWidget *mol_chooser_combobox = widget_from_builder("flip_map_hand_comboboxtext");
   gtk_widget_set_visible(frame, FALSE);
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(mol_chooser_combobox));
   flip_hand(imol);
   graphics_info_t g;
   g.graphics_grab_focus();

}

extern "C" G_MODULE_EXPORT
void
on_make_masked_maps_by_chain_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("make_masked_maps_by_chain_frame");
   GtkWidget *mol_chooser_combobox = widget_from_builder("make_masked_maps_by_chain_model_comboboxtext");
   GtkWidget *map_chooser_combobox = widget_from_builder("make_masked_maps_by_chain_map_comboboxtext");
   int imol     = my_combobox_get_imol(GTK_COMBO_BOX(mol_chooser_combobox));
   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(map_chooser_combobox));
   make_masked_maps_split_by_chain(imol, imol_map);
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_make_masked_maps_by_chain_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("make_masked_maps_by_chain_frame");
   gtk_widget_set_visible(frame, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_sharpen_blur_map_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame          = widget_from_builder("sharpen_blur_map_frame");
   GtkWidget *combobox       = widget_from_builder("sharpen_blur_map_comboboxtext");
   GtkWidget *checkbutton    = widget_from_builder("sharpen_blur_map_resample_checkbutton");
   GtkWidget *resample_entry = widget_from_builder("sharpen_blur_map_resample_entry");
   GtkWidget *b_factor_entry = widget_from_builder("sharpen_blur_map_entry");
   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   float resample_factor = 1.0;
   float b_factor = 0.0;
   if (b_factor_entry) {
      const char *t = gtk_editable_get_text(GTK_EDITABLE(b_factor_entry));
      try {
         b_factor = coot::util::string_to_float(std::string(t));
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
   }
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) {
      const char *t = gtk_editable_get_text(GTK_EDITABLE(resample_entry));
      try {
         resample_factor = coot::util::string_to_float(std::string(t));
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      // 20250115-PE make this non-blocking if you can (non-trivial)
      // you will need to split the calculation from the update of the gui and graphics.
      sharpen_blur_map_with_resampling(imol_map, b_factor, resample_factor);
   } else {
      // 20250115-PE make this non-blocking if you can
      sharpen_blur_map(imol_map, b_factor);
   }
   if (frame)
      gtk_widget_set_visible(frame, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_sharpen_blur_map_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                          G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("sharpen_blur_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_sharpen_blur_map_resample_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                 G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *resample_entry = widget_from_builder("sharpen_blur_map_resample_entry");
   GtkWidget *resample_label = widget_from_builder("sharpen_blur_map_resample_label");
   if (resample_entry) {
      if (gtk_check_button_get_active(checkbutton)) {
         gtk_widget_set_sensitive(resample_entry, TRUE);
         gtk_widget_set_sensitive(resample_label, TRUE);
      } else {
         gtk_widget_set_sensitive(resample_entry, FALSE);
         gtk_widget_set_sensitive(resample_label, FALSE);
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_ccp4i2_save_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                              G_GNUC_UNUSED gpointer         user_data) {

   auto get_first_model_molecule = [] () {
      graphics_info_t g;
      int imol = -1;
      int n_mol = g.molecules.size();
      for (int ii=0; ii<n_mol; ii++) {
         if (is_valid_model_molecule(ii))
            return ii;
      }
      return imol;
   };

   int imol = get_first_model_molecule();
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::filesystem::path save_dir("coot-ccp4i2");
      if (! std::filesystem::exists(save_dir))
         std::filesystem::create_directory(save_dir);
      std::string s = g.molecules[imol].stripped_save_name_suggestion();
      std::filesystem::path fn = save_dir / s;
      g.molecules[imol].save_coordinates(fn);
   }
}

extern "C" G_MODULE_EXPORT
void
on_copy_map_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                  G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_map_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                              G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_map_frame");
   GtkWidget *map_chooser_combobox = widget_from_builder("copy_map_comboboxtext");
   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(map_chooser_combobox));
   copy_molecule(imol_map);
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_molecule_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
				       G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy-molecule-frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_molecule_copy_button_clicked(G_GNUC_UNUSED GtkButton       *button,
				     G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame    = widget_from_builder("copy-molecule-frame");
   GtkWidget *combobox = widget_from_builder("copy_molecule_comboboxtext");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   copy_molecule(imol);
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_fragment_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame    = widget_from_builder("copy_fragment_frame");
   GtkWidget *combobox = widget_from_builder("copy_fragment_combobox");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   GtkWidget *entry = widget_from_builder("copy_fragment_atom_selection_entry");
   std::string text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));

   int imol_new = new_molecule_by_atom_selection(imol, text.c_str());
   GtkWidget *checkbutton = widget_from_builder("copy_fragment_move_molecule_here_checkbutton");
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton)))
      move_molecule_to_screen_centre_internal(imol_new);

   if (is_valid_model_molecule(imol_new))
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_fragment_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                       G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_fragment_frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_ncs_chain_copy_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                      G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_ncs_chain_frame");
   GtkWidget *combobox = widget_from_builder("copy_ncs_chain_molecule_combobox");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   GtkWidget *entry = widget_from_builder("copy_ncs_chain_entry");
   std::string text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
   copy_from_ncs_master_to_others(imol, text.c_str());
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_ncs_chain_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                        G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_ncs_chain_frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_ncs_residue_range_copy_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                                G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_ncs_residue_range_frame");
   GtkWidget *combobox = widget_from_builder("copy_ncs_residue_range_molecule_combobox");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   GtkWidget *entry_chain_id    = widget_from_builder("copy_ncs_residue_range_chain_entry");
   GtkWidget *entry_resno_start = widget_from_builder("copy_ncs_residue_range_start_residue_number_entry");
   GtkWidget *entry_resno_end   = widget_from_builder("copy_ncs_residue_range_end_residue_number_entry");
   std::string chain_id_text    = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry_chain_id)));
   std::string resno_start_text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry_resno_start)));
   std::string resno_end_text   = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry_resno_end)));
   try {
      int resno_start = coot::util::string_to_int(resno_start_text);
      int resno_end   = coot::util::string_to_int(resno_end_text);
      copy_residue_range_from_ncs_master_to_others(imol, chain_id_text.c_str(), resno_start, resno_end);
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
      logger.log(log_t::WARNING, logging::function_name_t("on_copy_ncs_residue_range_copy_button_clicked"),
                 "bad resno range", resno_start_text, resno_end_text);
   }
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_copy_ncs_residue_range_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                                G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("copy_ncs_residue_range_frame");
   gtk_widget_set_visible(frame, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_make_smooth_map_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                         G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("make_smooth_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_make_smooth_map_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                         G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *combobox_map    = widget_from_builder("make_smooth_map_comboboxtext");
   GtkWidget *combobox_factor = widget_from_builder("make_smooth_map_factor_comboboxtext");
   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(combobox_map));
   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox_factor));

   try {
      float factor = coot::util::string_to_float(std::string(t));
      smooth_map(imol_map, factor);
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

   GtkWidget *frame = widget_from_builder("make_smooth_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_make_difference_map_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                             G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("make_difference_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_make_difference_map_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                         G_GNUC_UNUSED gpointer         user_data) {


   GtkWidget *map_combobox_1 = widget_from_builder("make_difference_map_map_1_comboboxtext");
   GtkWidget *map_combobox_2 = widget_from_builder("make_difference_map_map_2_comboboxtext");
   int imol_map_1 = my_combobox_get_imol(GTK_COMBO_BOX(map_combobox_1));
   int imol_map_2 = my_combobox_get_imol(GTK_COMBO_BOX(map_combobox_2));
   const char *t = gtk_editable_get_text(GTK_EDITABLE(widget_from_builder("make_difference_map_scale_entry")));
   try{
      float scale = coot::util::string_to_float(std::string(t));
      difference_map(imol_map_1, imol_map_2, scale);
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

   GtkWidget *frame = widget_from_builder("make_difference_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_set_map_is_difference_map_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("set_map_is_difference_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_set_map_is_difference_map_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                               G_GNUC_UNUSED gpointer         user_data) {


   GtkWidget *map_combobox = widget_from_builder("set_map_is_difference_map_map_comboboxtext");
   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(map_combobox));
   // I could have another widget here - a switch to turn to back to not a difference map
   // in the future.
   set_map_is_difference_map(imol_map, TRUE);

   GtkWidget *frame = widget_from_builder("set_map_is_difference_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_make_an_average_map_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                             G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("make_an_average_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_make_an_average_map_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                         G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *box = widget_from_builder("make_an_average_map_box");
   GtkWidget *item_widget = gtk_widget_get_first_child(box);
   std::vector<std::pair<int, float> > imol_map_and_scale_vec;
   while (item_widget) {
      GtkWidget *map_combobox = gtk_widget_get_first_child(item_widget);
      GtkWidget *label = gtk_widget_get_next_sibling(map_combobox);
      GtkWidget *entry =gtk_widget_get_next_sibling(label);
      int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(map_combobox));
      const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
      if (t) {
         try {
            float scale = coot::util::string_to_float(std::string(t));
            imol_map_and_scale_vec.push_back(std::make_pair(imol_map, scale));
         }
         catch (const std::runtime_error &e) {
            std::cout << "WARNING::" << e.what() << std::endl;
         }
      } else {
         std::cout << "null t in on_make_an_average_map_ik_button_clicked()" << std::endl;
      }
      item_widget = gtk_widget_get_next_sibling(item_widget);
   };


   GtkWidget *frame = widget_from_builder("make_an_average_map_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_glyco_wta_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                   G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *w = widget_from_builder("glyco-wta-frame");
   if (w)
      gtk_widget_set_visible(w, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_glyco_wta_fit_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame    = widget_from_builder("glyco-wta-frame");
   GtkWidget *combobox = widget_from_builder("glyco_wta_glycosylation_name_comboboxtext");
   graphics_info_t g;
   std::string t = g.get_active_label_in_comboboxtext(GTK_COMBO_BOX_TEXT(combobox));
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      int imol_map = g.Imol_Refinement_Map();
      if (is_valid_map_molecule(imol_map)) {
         std::string tt;
         if (t == "NAG-NAG-BMA")            tt = "NAG-NAG-BMA";
         if (t == "High Mannose")           tt = "high-mannose";
         if (t == "Hybrid")                 tt = "hybrid";
         if (t == "Mammalian Bianntennary") tt = "mammalian-biantennary";
         if (t == "Plant Bianntennary")     tt = "plant-biantennary";
         clipper::Xmap<float> xmap = g.molecules[imol_map].xmap;
         coot::residue_spec_t res_spec(coot::atom_spec_t(aa.second));
         // coot::cho::add_named_glyco_tree(tt, &g.molecules[imol].atom_sel, imol, xmap, g.Geom_p(), as.chain_id, as.res_no);
         g.molecules[imol].add_named_glyco_tree(tt, g.Geom_p(), res_spec, xmap); // needs bonds update
         g.graphics_draw();
      } else {
         std::cout << "not a valid map " << imol_map << std::endl;
      }
   } else {
      std::cout << "not a valid model " << imol << std::endl;
   }
   if (frame)
      gtk_widget_set_visible(frame, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_tmblmf_cancel_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("transform_map_by_lsq_model_fit_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_tmblmf_ok_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                            G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *combobox_1 = widget_from_builder("tmblmf_comboboxtext_1");
   GtkWidget *combobox_2 = widget_from_builder("tmblmf_comboboxtext_2");
   int imol_model_1 = my_combobox_get_imol(GTK_COMBO_BOX(combobox_1));
   int imol_model_2 = my_combobox_get_imol(GTK_COMBO_BOX(combobox_2));

   GtkWidget *chain_id_1_entry    = widget_from_builder("tmblmf_chain_id_1_entry");
   GtkWidget *chain_id_2_entry    = widget_from_builder("tmblmf_chain_id_2_entry");
   GtkWidget *resno_start_1_entry = widget_from_builder("tmblmf_resno_start_1_entry");
   GtkWidget *resno_start_2_entry = widget_from_builder("tmblmf_resno_start_1_entry");
   GtkWidget *resno_end_1_entry   = widget_from_builder("tmblmf_resno_end_1_entry");
   GtkWidget *resno_end_2_entry   = widget_from_builder("tmblmf_resno_end_1_entry");

   GtkWidget *radius_entry = widget_from_builder("tmblmf_radius_entry");

   const char *chain_id_1_text    = gtk_editable_get_text(GTK_EDITABLE(chain_id_1_entry));
   const char *chain_id_2_text    = gtk_editable_get_text(GTK_EDITABLE(chain_id_2_entry));
   const char *resno_start_1_text = gtk_editable_get_text(GTK_EDITABLE(resno_start_1_entry));
   const char *resno_start_2_text = gtk_editable_get_text(GTK_EDITABLE(resno_start_2_entry));
   const char *resno_end_1_text   = gtk_editable_get_text(GTK_EDITABLE(resno_end_1_entry));
   const char *resno_end_2_text   = gtk_editable_get_text(GTK_EDITABLE(resno_end_2_entry));

   const char *radius_text = gtk_editable_get_text(GTK_EDITABLE(radius_entry));
   try {
      float radius = coot::util::string_to_float(std::string(radius_text));
      int resno_start_1 = coot::util::string_to_int(std::string(resno_start_1_text));
      int resno_start_2 = coot::util::string_to_int(std::string(resno_start_2_text));
      int resno_end_1   = coot::util::string_to_int(std::string(resno_end_1_text));
      int resno_end_2   = coot::util::string_to_int(std::string(resno_end_2_text));

      // get the matrix and call transform_map_raw()

      graphics_info_t g;
      int match_type = 1; // main - does this also work for RNA?
      coot::lsq_range_match_info_t match_info(resno_start_1, resno_end_1, chain_id_1_text,
                                              resno_start_2, resno_end_2, chain_id_2_text,
                                              match_type);

      g.lsq_matchers->push_back(match_info);

      // 20250128-PE note to self: use a transform_map() function where
      // I can pass the cell and rtop and spacegroup as objects,
      // c.f. coot::util::transform_map()

      std::pair<int, clipper::RTop_orth> status_and_rtop =
         g.apply_lsq(imol_model_1, imol_model_2, *g.lsq_matchers);
      int status = status_and_rtop.first;
      if (status == 1) {
         int imol_refinement_map = g.Imol_Refinement_Map();
         if (is_valid_map_molecule(imol_refinement_map)) {
            g.lsq_matchers->clear();
            const clipper::Xmap<float> &xmap = g.molecules[imol_refinement_map].xmap;
            clipper::Spacegroup spacegroup = xmap.spacegroup();
            clipper::Cell cell = xmap.cell();
            const clipper::RTop_orth &rtop = status_and_rtop.second;
            clipper::Mat33<double>  mat = rtop.rot();
            clipper::Vec3<double> trans = rtop.trn();
            clipper::Coord_orth rotation_centre = g.get_rotation_centre_co();
            float alpha = clipper::Util::rad2d(cell.alpha());  // need degrees
            float beta  = clipper::Util::rad2d(cell.beta());
            float gamma = clipper::Util::rad2d(cell.gamma());
            transform_map_raw(imol_refinement_map,
                              mat(0,0), mat(0,1), mat(0,2),
                              mat(1,0), mat(1,1), mat(1,2),
                              mat(2,0), mat(2,1), mat(2,2),
                              trans[0], trans[1], trans[2],
                              rotation_centre.x(), rotation_centre.y(), rotation_centre.z(),
                              radius, spacegroup.symbol_hm().c_str(),
                              cell.a(), cell.b(), cell.c(),
                              alpha, beta, gamma);
         }
      } else {
         std::cout << "WARNING:: in on_tmblmf_ok_button_clicked() bad matrix status" << std::endl;
      }
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

   GtkWidget *frame = widget_from_builder("transform_map_by_lsq_model_fit_frame");
   if (frame)
      gtk_widget_set_visible(frame, FALSE);
}

// this does not belong in main window. Move it.

extern "C" G_MODULE_EXPORT
void
on_hole_dialog_response(GtkDialog       *dialog,
                        gint             response_id,
                        gpointer         user_data) {

   std::cout << "-------------------------------- on_hole_dialog_response() " << response_id << std::endl;

   if (response_id == GTK_RESPONSE_OK) {
      std::cout << "here.... " << std::endl;
      graphics_info_t g;
      clipper::Coord_orth start = g.hole_start;
      clipper::Coord_orth   end = g.hole_end;
      GtkWidget *commbobox = widget_from_builder("hole_model_comboboxtext");
      int imol = my_combobox_get_imol(GTK_COMBO_BOX(commbobox));
      std::cout << "here.... with imol " << imol << std::endl;
      if (is_valid_model_molecule(imol)) {
         float colour_map_multiplier = 1.0;
         float colour_map_offset = 0.0;
         int n_runs = 3000;
         bool show_probe_radius_graph_flag = true;
         GtkWidget *entry = widget_from_builder("hole_surface_dots_file_name_entry");
         const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (t) {
            std::string export_surface_dots_file_name = std::string(t);
            hole(imol, start.x(), start.y(), start.z(), end.x(), end.y(), end.z(),
                 colour_map_multiplier, colour_map_offset,
                 n_runs, show_probe_radius_graph_flag,
                 export_surface_dots_file_name);
         }
      }
   }

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_hole_set_start_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                 G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("hole_start_entry");
   const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      graphics_info_t g;
      g.set_hole_start();
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }
}

extern "C" G_MODULE_EXPORT
void
on_hole_set_end_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                               G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("hole_end_entry");
   const char *t = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      graphics_info_t g;
      g.set_hole_end();
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }
}

void fill_logging_text_view() {

   GtkWidget *textview = widget_from_builder("logging_textview");
   GtkTextBuffer *buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textview));
   GtkTextIter end_iter;

   std::vector<logging::log_item> lv = logger.get_log_history_from(graphics_info_t::logging_line_index);

   if (! lv.empty()) {

      // Create a text tag for the color
      //
      GtkTextTagTable *tag_table = gtk_text_buffer_get_tag_table(buffer);

      GtkTextTag *datetime_color_tag;

      datetime_color_tag = gtk_text_tag_table_lookup(tag_table, "datetime-color");
      if (datetime_color_tag == NULL) {
         datetime_color_tag = gtk_text_tag_new("datetime-color");
         std::string datetime_colour = "#8888aa";
         g_object_set(datetime_color_tag, "foreground", datetime_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), datetime_color_tag);
      }
      //
      GtkTextTag *info_color_tag = gtk_text_tag_table_lookup(tag_table, "info-type-color");
      if (info_color_tag == NULL) {
         info_color_tag = gtk_text_tag_new("info-type-color");
         std::string info_type_colour = "#00eeee";
         g_object_set(info_color_tag, "foreground", info_type_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), info_color_tag);
      }
      //
      GtkTextTag *debug_color_tag = gtk_text_tag_table_lookup(tag_table, "debug-type-color");
      if (debug_color_tag == NULL) {
         debug_color_tag = gtk_text_tag_new("debug-type-color");
         std::string debug_type_colour = "#999999";
         g_object_set(debug_color_tag, "foreground", debug_type_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), debug_color_tag);
      }
      //
      GtkTextTag *warning_color_tag = gtk_text_tag_table_lookup(tag_table, "warning-type-color");
      if (warning_color_tag == NULL) {
         GtkTextTag *warning_color_tag = gtk_text_tag_new("warning-type-color");
         std::string warning_type_colour = "#ff9900";
         g_object_set(warning_color_tag, "foreground", warning_type_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), warning_color_tag);
      }
      //
      GtkTextTag *error_color_tag = gtk_text_tag_table_lookup(tag_table, "error-type-color");
      if (error_color_tag == NULL) {
         GtkTextTag *error_color_tag = gtk_text_tag_new("error-type-color");
         std::string error_type_colour = "#ee2222";
         g_object_set(error_color_tag, "foreground", error_type_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), error_color_tag);
      }
      //
      GtkTextTag *gl_error_color_tag = gtk_text_tag_table_lookup(tag_table, "gl-error-type-color");
      if (gl_error_color_tag == NULL) {
         GtkTextTag *gl_error_color_tag = gtk_text_tag_new("gl-error-type-color");
         std::string gl_error_type_colour = "#cc8822";
         g_object_set(gl_error_color_tag, "foreground", gl_error_type_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), gl_error_color_tag);
      }
      //
      GtkTextTag *default_color_tag = gtk_text_tag_table_lookup(tag_table, "default-type-color");
      if (default_color_tag == NULL) {
         GtkTextTag *default_color_tag = gtk_text_tag_new("default-type-color");
         std::string default_type_colour = "#bbbbbb";
         g_object_set(default_color_tag, "foreground", default_type_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), default_color_tag);
      }
      //
      GtkTextTag *function_name_color_tag = gtk_text_tag_table_lookup(tag_table, "function-name-color");
      if (function_name_color_tag == NULL) {
         GtkTextTag *function_name_color_tag = gtk_text_tag_new("function-name-color");
         std::string function_name_colour = "#cc99cc";
         g_object_set(function_name_color_tag, "foreground", function_name_colour.c_str(), NULL);
         gtk_text_tag_table_add(gtk_text_buffer_get_tag_table(buffer), function_name_color_tag);
      }

      for (const auto &item : lv) {

	 // Insert the colored text

	 gtk_text_buffer_get_end_iter(buffer, &end_iter);

	 std::string log_type_text = item.type_as_string();
	 std::string color_tag_name = "default-type-color";
	 if (item.type == log_t::WARNING)  color_tag_name = "warning-type-color";
	 if (item.type == log_t::DEBUG)    color_tag_name = "debug-type-color";
	 if (item.type == log_t::INFO)     color_tag_name = "info-type-color";
	 if (item.type == log_t::ERROR)    color_tag_name = "error-type-color";
	 if (item.type == log_t::GL_ERROR) color_tag_name = "gl-error-type-color";

	 gtk_text_buffer_insert_with_tags_by_name(buffer, &end_iter, log_type_text.c_str(),
						  log_type_text.size(), color_tag_name.c_str(),
						  NULL);

	 // insert the date

	 std::string d = item.get_date_string();
	 gtk_text_buffer_get_end_iter(buffer, &end_iter);
	 // gtk_text_buffer_insert(buffer, &end_iter, d.c_str(), d.size());
	 gtk_text_buffer_insert_with_tags_by_name(buffer, &end_iter, d.c_str(), d.size(),
						  "datetime-color", NULL);

	 // insert the function
	 if (! item.function_name.empty()) {
	    std::string fn = std::string(" ") + item.function_name.get_name();
	    gtk_text_buffer_get_end_iter(buffer, &end_iter);
	    gtk_text_buffer_insert_with_tags_by_name(buffer, &end_iter, fn.c_str(), fn.size(),
						     "function-name-color", NULL);
	 }

	 // insert the message

	 std::string t = std::string(" ") + item.message;
	 gtk_text_buffer_get_end_iter(buffer, &end_iter);
	 gtk_text_buffer_insert(buffer, &end_iter, t.c_str(), t.size());

	 // insert a new-line

	 gtk_text_buffer_get_end_iter(buffer, &end_iter);
	 gtk_text_buffer_insert(buffer, &end_iter, "\n", 1); // Add newline character

      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_show_logging_menubutton_checkbutton_toggled(GtkCheckButton *checkbutton,
					       G_GNUC_UNUSED gpointer user_data) {

   GtkWidget *vbox = widget_from_builder("logging_vbox");
   if (gtk_check_button_get_active(checkbutton)) {
      gtk_widget_set_visible(vbox, TRUE);
      fill_logging_text_view();
   } else {
      gtk_widget_set_visible(vbox, FALSE);
   }

}
