/* src/glade-callbacks.cc
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2008 The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
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


#include "Python.h"

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
#include "c-interface-preferences.h"
#include "rotate-translate-modes.hh"
#include "restraints-editor-c.h"
#include "generic-display-objects-c.h"
#include "c-interface-refmac.h"
#include "gtk-widget-conversion-utils.h"
#include "curlew.h"
#include "read-phs.h"
#include "gtk-manual.h"
#include "c-interface-refine.h"
#include "widget-from-builder.hh"

// from support.h
// GtkWidget* lookup_widget (GtkWidget *widget, const gchar *widget_name);
#include "support.h"

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;

#include <vector>
#include "utils/coot-utils.hh"
#include "graphics-info.h"

#include "cc-interface.hh" // for read_ccp4_map()


// Let's put the new refinement and regularization control tools together here
// (although not strictly main window)


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refine_control_button_clicked(GtkButton       *button,
                                                                   gpointer         user_data) {
   graphics_info_t::show_refinement_and_regularization_parameters_frame();
}



// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_torsion_weight_combobox_changed(GtkComboBox     *combobox,
//                                                                      gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_rama_restraints_combobox_changed(GtkComboBox     *combobox,
//                                                                       gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_strand_rest_radiobutton_toggled(GtkToggleButton *togglebutton,
//                                                                     gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_helix_rest_radiobutton_toggled(GtkToggleButton *togglebutton,
//                                                                    gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_no_rest_radiobutton_toggled(GtkToggleButton *togglebutton,
//                                                                 gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_use_torsions_checkbutton_toggled(GtkToggleButton *togglebutton,
//                                                                       gpointer         user_data) {
// }



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_select_map_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                           G_GNUC_UNUSED gpointer         user_data) {

   show_select_map_dialog();
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
on_model_toolbar_refine_togglebutton_toggled (GtkToggleButton *toggletoolbutton,
                                              gpointer         user_data) {

   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active)
      do_refine(1);
   else
      do_refine(0);		/* unclick button */
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_regularize_togglebutton_toggled(GtkToggleButton *toggletoolbutton,
                                                                     gpointer         user_data) {

   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active)
      do_regularize(1);
   else
      do_regularize(0);		/* unclick button */
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_fixed_atoms_menubutton_activated(GtkMenuButton *button,
                                                  gpointer       user_data) {
   GtkWidget *w = wrapped_create_fixed_atom_dialog();
   gtk_widget_show(w);
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
on_model_toolbar_rotamers_button_clicked(GtkButton *toggletoolbutton,
                                         gpointer         user_data) {

   // Find rotamers for the residue at the centre of the screen
   graphics_info_t g;
   std::pair<int, mmdb::Atom *> aa = g.get_active_atom();
   int imol = aa.first;
   if (is_valid_model_molecule(imol)) {
      g.do_rotamers(imol, aa.second);
   }

}


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
on_model_toolbar_torsion_general_toggletoolbutton_toggled(GtkToggleButton *toggletoolbutton,
                                                                              gpointer         user_data) {

  gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
  if (active) {
    setup_torsion_general(1);
  } else {
    setup_torsion_general(0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_flip_peptide_button_clicked(GtkButton *button,
                                             gpointer   user_data) {
   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   if (is_valid_model_molecule(imol)) {
      auto &m = g.molecules[imol];
      coot::atom_spec_t atom_spec(active_atom.second);
      m.pepflip(atom_spec);
      g.graphics_draw();
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_side_chain_180_button_clicked(GtkButton *button,
                                               gpointer   user_data) {

#if 0 // do it directly these days
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
  if (active)
    setup_180_degree_flip(1);
  else
    setup_180_degree_flip(0);
#endif

   // look at check_if_in_180_degree_flip_define, it checks for
   // intermediate atoms - we should do that here too.

   std::cout << "Here in on_model_toolbar_add_terminal_residue_button_clicked()" << std::endl;
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
      int istatus = m.do_180_degree_side_chain_flip(spec.chain_id, spec.res_no, spec.ins_code, alt_conf, g.Geom_p());
      g.graphics_draw();
   }

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_backbone_torsions_toggletoolbutton_toggled(GtkToggleButton *toggletoolbutton,
                                                                                     gpointer         user_data) {

  gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
  if (active) {
    setup_backbone_torsion_edit(1);
  } else {
    setup_backbone_torsion_edit(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_mutate_and_autofit_menubutton_activate
                                        (GtkMenuButton *menubutton,
                                        gpointer         user_data)
{
   // this function seems not to be called on menu button click. Hmmm.
   std::cout << "on_model_toolbar_mutate_and_autofit_menubutton_active "
            << " select the right menu for mutate_and_autofit_menubutton here" << std::endl;
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_simple_mutate_togglebutton_toggled
                                        (GtkToggleButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(toggletoolbutton));
   if (active)
      setup_mutate(1);
   else
      setup_mutate(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_terminal_residue_button_clicked(GtkButton *button,
                                                     gpointer   user_data)
{
   graphics_info_t g;
   auto active_atom = g.get_active_atom();
   int imol = active_atom.first;
   std::cout << "on_model_toolbar_add_terminal_residue_button_clicked 1" << std::endl;
   if (is_valid_model_molecule(imol)) {
      std::cout << "on_model_toolbar_add_terminal_residue_button_clicked 2" << std::endl;
      mmdb::Residue *residue_p = active_atom.second->residue;
      coot::atom_spec_t atom_spec(active_atom.second);
      std::string chain_id = atom_spec.chain_id;
      coot::residue_spec_t residue_spec(atom_spec);
      std::string new_type = "ALA";
      mmdb::Atom *atom = active_atom.second;
      std::string terminus_type = g.molecules[imol].get_term_type(atom);
      g.execute_add_terminal_residue(imol, terminus_type, residue_p,
                                     chain_id, new_type, true);
   }

}



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_alt_conf_toolbutton_clicked
                                        (GtkToggleButton   *toolbutton,
                                         gpointer         user_data) {
  altconf();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_atom_button_clicked(GtkButton       *button,
                                                             gpointer         user_data) {
   place_atom_at_pointer();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_clear_pending_picks_button_clicked(GtkButton       *button,
                                                                        gpointer         user_data) {
   clear_pending_picks();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_delete_button_clicked(GtkButton       *button,
                                                           gpointer         user_data) {
  GtkWidget *widget = wrapped_create_delete_item_dialog();
  gtk_widget_show(widget);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_undo_button_clicked   (GtkButton       *button,
                                                            gpointer         user_data) {
   apply_undo();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_redo_button_clicked (GtkButton       *button,
                                                         gpointer         user_data)
{
  apply_redo();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refmac_button_clicked (GtkToggleButton   *toolbutton,
                                                            gpointer         user_data)
{
  /* wrapped_create_run_refmac_dialog(); */
  wrapped_create_simple_refmac_dialog();

}


extern "C" G_MODULE_EXPORT
void
on_refine_params_torsion_weight_combobox_changed(GtkComboBox     *combobox,
                                                                     gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_torsion_weight_from_text(active_item_idx, t);
}

extern "C" G_MODULE_EXPORT
void
on_refine_params_rama_restraints_combobox_changed (GtkComboBox     *combobox,
                                                                       gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_ramachandran_restraints_weight_from_text(active_item_idx, t);
}



extern "C" G_MODULE_EXPORT
void
on_ssm_superposition1_activate         (GMenuItem     *menuitem,
                                                            gpointer         user_data) {
   GtkWidget *w = wrapped_create_superpose_dialog(); // uses builder

   /* we get returned w = 0 when there is no MMDBSSM. (We are doing it
      this way because we don't have to introduce HAVE_MMDBSSM into the
      *c* compiler arguments (this is simpler)).  */
  if (w)
     gtk_widget_show(w);
}



// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_strand_rest_radiobutton_toggled(GtkToggleButton *togglebutton,
//n                                                                    gpointer         user_data) {
//  if (gtk_toggle_button_get_active(togglebutton))
//    set_secondary_structure_restraints_type(2);
//}

// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_helix_rest_radiobutton_toggled(GtkToggleButton *togglebutton,
//                                                                    gpointer         user_data) {
//   if (gtk_toggle_button_get_active(togglebutton))
//     set_secondary_structure_restraints_type(1);
// }

// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_no_rest_radiobutton_toggled(GtkToggleButton *togglebutton,
//                                                                 gpointer         user_data) {
//   if (gtk_toggle_button_get_active(togglebutton))
//     set_secondary_structure_restraints_type(0);
// }

// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_use_torsions_checkbutton_toggled(GtkToggleButton *togglebutton,
//                                                                       gpointer         user_data) {
//    do_torsions_toggle(GTK_WIDGET(togglebutton));
// }



extern "C" G_MODULE_EXPORT
void
on_open_coordinates1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  open_coords_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_open_dataset1_activate (GMenuItem     *menuitem,
                                               gpointer         user_data) {

   GtkWidget *dataset_chooser = widget_from_builder("dataset_filechooser_dialog");
   set_directory_for_filechooser(dataset_chooser);
   set_file_selection_dialog_size(dataset_chooser);
   add_filechooser_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);
   gtk_widget_show(dataset_chooser);
   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);

}

extern "C" G_MODULE_EXPORT
void
on_auto_open_mtz_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data) {

   int is_auto_read_fileselection = 1;
   int is;
   GtkWidget *file_filter_button;
   GtkWidget *dataset_chooser = coot_dataset_chooser();

   set_directory_for_filechooser(dataset_chooser);

   // add_ccp4i_project_optionmenu(dataset_fileselection1, COOT_DATASET_FILE_SELECTION);

   file_filter_button = add_filename_filter_button(dataset_chooser, COOT_DATASET_FILE_SELECTION);

   // sort_button = add_sort_button_fileselection(dataset_chooser);
   /*   set_directory_for_fileselection(dataset_fileselection1); */

   /* stuff in user data saying if this is autoread or not... */
   is = is_auto_read_fileselection;
   g_object_set_data(G_OBJECT(dataset_chooser), "imol", GINT_TO_POINTER(is));

   // set_file_selection_dialog_size(dataset_chooser);

   g_object_set_data(G_OBJECT(dataset_chooser), "is_auto", GINT_TO_POINTER(is_auto_read_fileselection));

   set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_chooser);
   gtk_widget_show(dataset_chooser);

   // what does this do?
   // push_the_buttons_on_fileselection(file_filter_button, sort_button, dataset_chooser);
}


extern "C" G_MODULE_EXPORT
void
on_save_coordinates1_activate          (GMenuItem     *menuitem,
                                                            gpointer         user_data)
{
   GCallback callback_func = G_CALLBACK(save_molecule_coords_combobox_changed);
   int imol = first_coords_imol();
   int imol_unsaved = first_unsaved_coords_imol();
   if (imol_unsaved != -1)
      imol = imol_unsaved;
   std::cout << "DEBUG:: in on_save_coordinates1_activate() with imol_unsaved "
             << imol_unsaved << std::endl;
   set_save_molecule_number(imol); /* set *save* molecule number */

   // this is the molecule chooser, not the file chooser
   //
   GtkWidget *widget = widget_from_builder("save_coords_dialog");
   GtkWidget *combobox = widget_from_builder("save_coordinates_combobox");

   if (combobox) {
      fill_combobox_with_coordinates_options(combobox, callback_func, imol);
      set_transient_and_position(COOT_UNDEFINED_WINDOW, widget);
      gtk_widget_show(widget);
      gtk_window_present(GTK_WINDOW(widget));
   } else {
      std::cout << "ERROR:: in on_save_coordinates1_activate() bad combobox!\n";
   }
}


extern "C" G_MODULE_EXPORT
void
on_save_coordinates_filechooser_dialog_response(GtkDialog       *dialog,
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

         // imol set in
         // on_save_coords_dialog_save_button_clicked(GtkButton       *button,
         //                                                               gpointer         user_data)
         int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "imol"));
         save_coordinates(imol, fnc);
      }
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
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
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
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
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   gtk_widget_hide(GTK_WIDGET(dialog));
   
}


extern "C" G_MODULE_EXPORT
void
on_map_name_filechooser_dialog_file_activated(GtkFileChooser* dialog,
                                                                  gpointer user_data) {

#if 0 // 20220809-PE well, today it seems to cause a double read of the map! Strange

   // 20220319-PE shouldn't need to connect to this says the documentation - hmmm.....
   // const char *fn = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
   GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
   GError *error = NULL;
   GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                            G_FILE_QUERY_INFO_NONE, NULL, &error);
   const char *fnc = g_file_info_get_name(file_info);
   bool is_diff_map = false;
   GtkWidget *checkbutton = widget_from_builder("map_filechooser_is_difference_map_button");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
      is_diff_map = true;
   read_ccp4_map(fnc, is_diff_map);

   gtk_widget_hide(GTK_WIDGET(dialog));
#endif

}


// This is not a main-window callback! Move it - and others like it.
//
extern "C" G_MODULE_EXPORT
void
on_coords_filechooser_dialog_response(GtkDialog       *dialog,
                                                          gint             response_id,
                                                          gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {

      GtkWidget *recentre_combobox = widget_from_builder("coords_filechooserdialog_recentre_combobox");
      int active_item_index = gtk_combo_box_get_active(GTK_COMBO_BOX(recentre_combobox));

      // 20220601-PE I should read multiple GFiles here, I suppose

#if 0
      GSList *files_list = gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(dialog));
      while (files_list) {

         const char *fnc = static_cast<const char *>(files_list->data);
         if (fnc) {
            std::string fn(fnc);
            handle_read_draw_molecule_with_recentre(fn, 0);
         }
         files_list = g_slist_next(files_list);
      }
#endif

      // const char *fn = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *fn = g_file_info_get_name(file_info);
      save_directory_from_filechooser(GTK_WIDGET(dialog));

      bool recentre_on_read_pdb_flag = false;
      bool move_molecule_here_flag = false;
      if (active_item_index == 0)
         recentre_on_read_pdb_flag = true;
      if (active_item_index == 1)
         recentre_on_read_pdb_flag = false;
      if (active_item_index == 2)
         move_molecule_here_flag = true;

      if (move_molecule_here_flag) {
         handle_read_draw_molecule_and_move_molecule_here(fn);
      } else {
         if (recentre_on_read_pdb_flag)
            handle_read_draw_molecule_with_recentre(fn, 1);
         else
            handle_read_draw_molecule_with_recentre(fn, 0); // no recentre
      }

      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}

extern "C" G_MODULE_EXPORT
void
on_coords_filechooser_dialog_file_activated(GtkFileChooser* dialog,
                                                                gpointer user_data) {

   // 20220427-PE I am not sure that this is needed now. Double click on a filename
   // seems to work for datasets - and they don't have a file_activated callback.

   // 20220319-PE shouldn't need to connect to this says the documentation - hmmm.....
   // const char *fn = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
   //  handle_read_draw_molecule_with_recentre(fn, 1);
   save_directory_from_filechooser(GTK_WIDGET(dialog));
   gtk_widget_hide(GTK_WIDGET(dialog));

}


#include "rsr-functions.hh"

extern "C" G_MODULE_EXPORT
void
on_menubar_regularize_residue_activate(GMenuItem *menuitem,
                                                           gpointer     user_data) {

   regularize_residue();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_regularize_tandem_3_activate(GMenuItem *menuitem,
                                                            gpointer     user_data) {
   regularize_tandem_3();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_regularize_sphere_activate(GMenuItem *menuitem,
                                                          gpointer     user_data) {
   regularize_sphere();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_single_residue_activate(GMenuItem *menuitem,
                                                           gpointer     user_data) {

   rsr_refine_residue();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_tandem_5_activate(GMenuItem *menuitem,
                                                     gpointer     user_data) {

   rsr_refine_tandem_5();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_tandem_3_activate(GMenuItem *menuitem,
                                                     gpointer     user_data) {
   rsr_refine_tandem_3();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_sphere_plus_activate(GMenuItem *menuitem,
                                                        gpointer     user_data) {

   rsr_sphere_refine_plus();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_chain_activate(GMenuItem *menuitem,
                                                  gpointer     user_data) {

   rsr_refine_chain();
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_residue_range_activate(GMenuItem *menuitem,
                                      gpointer     user_data) {

   do_refine(1);
}

extern "C" G_MODULE_EXPORT
void
on_menubar_rsr_sphere_activate(GMenuItem *menuitem,
                                                   gpointer     user_data) {
   rsr_sphere_refine();
}


extern "C" G_MODULE_EXPORT
void
on_delete_item_atom_item_activate(GMenuItem *menuitem,
                                                      gpointer     user_data) {
   set_delete_atom_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_water_item_activate(GMenuItem *menuitem,
                                                      gpointer     user_data) {

   set_delete_water_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_sidechain_item_activate(GMenuItem *menuitem,
                                                           gpointer     user_data) {
   set_delete_sidechain_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_sidechains_in_residue_range_item_activate(GMenuItem *menuitem,
                                                                             gpointer     user_data) {
   set_delete_sidechain_range_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_residue_hydrogen_atoms_item_activate(GMenuItem *menuitem,
                                                                        gpointer     user_data) {

   set_delete_residue_hydrogens_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_residue_item_activate(GMenuItem *menuitem,
                                                         gpointer     user_data) {

   set_delete_residue_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_residue_range_item_activate(GMenuItem *menuitem,
                                           gpointer     user_data) {
   set_delete_residue_zone_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_chain_item_activate(GMenuItem *menuitem,
                                   gpointer     user_data) {

   set_delete_chain_mode();
}

extern "C" G_MODULE_EXPORT
void
on_calculate_updating_maps1_activate(GMenuItem *menuitem,
                                     gpointer     user_data) {

   show_calculate_updating_maps_pythonic_gui();

}



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_icons_menubar_icons_item_activate(GMenuItem *menuitem,
                                                                       gpointer     user_data) {

   GtkWidget *tb;

   tb = widget_from_builder("main_window_model_toolbar_second_top");
   // gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);
   tb = widget_from_builder("main_window_model_toolbar_lower");
   // gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);
   tb = widget_from_builder("main_window_model_toolbar_bottom");
   // gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);

   GtkWidget *mi = widget_from_builder("rotate_translate_item_menu_item_top");
   // gtk_menu_item_set_label(GTK_MENU_ITEM(mi), "");
   mi = widget_from_builder("menubar_for_rsr_item_menuitem");
   // gtk_menu_item_set_label(GTK_MENU_ITEM(mi), "");
   mi = widget_from_builder("menubar_for_delete_items_menu_item_top");
   // gtk_menu_item_set_label(GTK_MENU_ITEM(mi), "");
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_icons_menubar_icons_and_text_item_activate(GMenuItem *menuitem,
                                                                                gpointer     user_data) {

   GtkWidget *tb;

   tb = widget_from_builder("main_window_model_toolbar_second_top");
   // gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_BOTH_HORIZ);
   tb = widget_from_builder("main_window_model_toolbar_lower");
   // gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_BOTH_HORIZ);
   tb = widget_from_builder("main_window_model_toolbar_bottom");
   // gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_BOTH_HORIZ);

   GtkWidget *mi = widget_from_builder("rotate_translate_item_menu_item_top");
   // gtk_menu_item_set_label(GTK_MENU_ITEM(mi), "Rotate/Translate");
   mi = widget_from_builder("menubar_for_rsr_item_menuitem");
   // gtk_menu_item_set_label(GTK_MENU_ITEM(mi), "Real Space Refinement");
   mi = widget_from_builder("menubar_for_delete_items_menu_item_top");
   // gtk_menu_item_set_label(GTK_MENU_ITEM(mi), "   Delete");

}


extern "C" G_MODULE_EXPORT
void
on_ribbons_colour_by_chain_menu_item_activate(GMenuItem *menuitem,
                                                                  gpointer     user_data) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "Chain";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}

extern "C" G_MODULE_EXPORT
void
on_ribbons_colour_rainbow_menu_item_activate(GMenuItem *menuitem,
                                                                 gpointer     user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "colorRampChainsScheme";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}

extern "C" G_MODULE_EXPORT
void
on_ribbons_colour_by_secondary_structure_menu_item_activate(GMenuItem *menuitem,
                                                                 gpointer     user_data) {
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string colour_scheme = "colorBySecondaryScheme";
      std::string atom_selection = "//";
      std::string style = "Ribbon";
      graphics_info_t g;
      int status = g.add_molecular_representation(imol, atom_selection, colour_scheme, style);
   }
}

extern "C" G_MODULE_EXPORT
void
on_draw_perspective_perspective_menu_item_activate(GMenuItem *menuitem,
                                                                       gpointer     user_data) {
   set_use_perspective_projection(1);
}

extern "C" G_MODULE_EXPORT
void
on_draw_perspective_orthographic_menu_item_activate(GMenuItem *menuitem,
                                                                        gpointer     user_data) {
   set_use_perspective_projection(0);
}

extern "C" G_MODULE_EXPORT
void
on_about1_activate(GMenuItem     *menuitem,
                                       gpointer         user_data) {

   // GtkWidget *about_window = create_aboutdialog();
   GtkWidget *about_window = widget_from_builder("aboutdialog");
   add_coot_references_button(about_window);
   gtk_widget_show(about_window);
}


extern "C" G_MODULE_EXPORT
void
on_main_window_resize_window_up_arrow_clicked(GtkButton       *button,
                                                                  gpointer         user_data) {

   GtkWidget *window = graphics_info_t::get_main_window();
   GtkAllocation allocation;
   gtk_widget_get_allocation(window, &allocation);
   int w = allocation.width;
   int h = allocation.height;
   int h_new = h - 30;
   gtk_window_set_resizable(GTK_WINDOW(window), TRUE);
#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "Resizing: A number of GdkWindow APIs are no longer available. "
   "This includes gdk_window_reparent(), gdk_window_set_geometry_hints(), "
   "gdk_window_raise(), gdk_window_restack(), gdk_window_move(), "
   "gdk_window_resize(). If you need to manually control the position "
   "or stacking of your X11 windows, you you will have to use Xlib apis." << std::endl;
#else
   gtk_window_resize(GTK_WINDOW(window), w, h_new);
#endif
   // gtk_window_set_resizable(GTK_WINDOW(window), FALSE);

   // Note to self gtk_window_set_resizable() also expands the window fully in Y. Bleugh.
}

extern "C" G_MODULE_EXPORT
void
on_main_window_resize_window_down_arrow_clicked(GtkButton       *button,
                                                                  gpointer         user_data) {
   GtkWidget *window = graphics_info_t::get_main_window();
   GtkAllocation allocation;
   gtk_widget_get_allocation(window, &allocation);
   int w = allocation.width;
   int h = allocation.height;
   int h_new = h + 30;
   gtk_window_set_resizable(GTK_WINDOW(window), TRUE);
#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "Resizing: A number of GdkWindow APIs are no longer available. "
   "This includes gdk_window_reparent(), gdk_window_set_geometry_hints(), "
   "gdk_window_raise(), gdk_window_restack(), gdk_window_move(), "
   "gdk_window_resize(). If you need to manually control the position "
   "or stacking of your X11 windows, you you will have to use Xlib apis." << std::endl;
#else
   gtk_window_resize(GTK_WINDOW(window), w, h_new);
#endif
   // gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_main_window_resize_window_left_arrow_clicked(GtkButton       *button,
                                                                  gpointer         user_data) {
   GtkWidget *window = graphics_info_t::get_main_window();
   GtkAllocation allocation;
   gtk_widget_get_allocation(window, &allocation);
   int w = allocation.width;
   int h = allocation.height;
   int w_new = w - 30;
   gtk_window_set_resizable(GTK_WINDOW(window), TRUE);
#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "Resizing: A number of GdkWindow APIs are no longer available. "
   "This includes gdk_window_reparent(), gdk_window_set_geometry_hints(), "
   "gdk_window_raise(), gdk_window_restack(), gdk_window_move(), "
   "gdk_window_resize(). If you need to manually control the position "
   "or stacking of your X11 windows, you you will have to use Xlib apis." << std::endl;
#else
   gtk_window_resize(GTK_WINDOW(window), w_new, h);
#endif
   // gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_main_window_resize_window_right_arrow_clicked(GtkButton       *button,
                                                 gpointer         user_data) {

   GtkWidget *window = graphics_info_t::get_main_window();
   GtkAllocation allocation;
   gtk_widget_get_allocation(window, &allocation);
   int w = allocation.width;
   int h = allocation.height;
   int w_new = w + 30;
   gtk_window_set_resizable(GTK_WINDOW(window), TRUE);
#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "Resizing: A number of GdkWindow APIs are no longer available. "
   "This includes gdk_window_reparent(), gdk_window_set_geometry_hints(), "
   "gdk_window_raise(), gdk_window_restack(), gdk_window_move(), "
   "gdk_window_resize(). If you need to manually control the position "
   "or stacking of your X11 windows, you you will have to use Xlib apis." << std::endl;
#else
   gtk_window_resize(GTK_WINDOW(window), w_new, h);
#endif
   // gtk_window_set_resizable(GTK_WINDOW(window), FALSE);
}
