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


// Let's put the new refinement and regularization control tools together here
// (although not strictly main window)


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refine_control_button_clicked_gtkbuilder_callback(GtkButton       *button,
                                                                   gpointer         user_data)
{
   graphics_info_t g;
   g.show_refinement_and_regularization_parameters_dialog();
}




extern "C" G_MODULE_EXPORT
void
on_refinement_and_regularization_parameters_dialog_close_gtkbuilder_callback(GtkDialog *dialog,
                                                                              gpointer   user_data) {
   gtk_widget_hide(GTK_WIDGET(dialog));
}


extern "C" G_MODULE_EXPORT
void
on_refinement_and_regularization_parameters_dialog_response_gtkbuilder_callback(GtkDialog       *dialog,
                                                                                gint             response_id,
                                                                                gpointer         user_data) {

   std::cout << "in on_refinement_and_regularization_parameters_dialog_response_gtkbuilder_callback() got response_id " << response_id << std::endl;
   if (response_id == GTK_RESPONSE_CLOSE) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}


// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_torsion_weight_combobox_changed_gtkbuilder_callback(GtkComboBox     *combobox,
//                                                                      gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_rama_restraints_combobox_changed_gtkbuilder_callback(GtkComboBox     *combobox,
//                                                                       gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_strand_rest_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                     gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_helix_rest_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                    gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_no_rest_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                 gpointer         user_data) {
// }


// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_use_torsions_checkbutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                       gpointer         user_data) {
// }



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_select_map_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   show_select_map_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refine_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_refine(1);
  else
    do_refine(0);		/* unclick button */

}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_regularize_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active)
      do_regularize(1);
   else
      do_regularize(0);		/* unclick button */
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_fixed_atoms_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_fixed_atom_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rigid_body_fit_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active) {
      printf("Rigid Body:\n");
      do_rigid_body_refine(1);
   } else {
      do_rigid_body_refine(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rot_trans_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active) {
      do_rot_trans_setup(1);
   } else {
      do_rot_trans_setup(0);
   }
}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_auto_fit_rotamer_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
     setup_auto_fit_rotamer(1);
  else
    setup_auto_fit_rotamer(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rotamers_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active)
      setup_rotamers(1);
   else
      setup_rotamers(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_chi_angles_togglebutton_toggled_gtkbuilder_callback
                                              (GtkToggleToolButton *toggletoolbutton,
                                               gpointer         user_data)
{
   gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active) {
      setup_edit_chi_angles(1);
   } else {
      setup_edit_chi_angles(0);
      set_show_chi_angle_bond(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_torsion_general_toggletoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{

  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active) {
    setup_torsion_general(1);
  } else {
    setup_torsion_general(0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_flip_peptide_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_pepflip(1);
  else
     do_pepflip(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_sidechain_180_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    setup_180_degree_flip(1);
  else
    setup_180_degree_flip(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_edit_backbone_torsions_toggletoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{

  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active) {
    setup_backbone_torsion_edit(1);
  } else {
    setup_backbone_torsion_edit(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_mutate_and_autofit_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    setup_mutate_auto_fit(1);
  else
     setup_mutate_auto_fit(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_simple_mutate_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active)
      setup_mutate(1);
   else
      setup_mutate(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_terminal_residue_togglebutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_add_terminal_residue(1);
  else
    do_add_terminal_residue(0);
}



extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_alt_conf_toolbutton_clicked_gtkbuilder_callback
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  altconf();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_add_atom_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   place_atom_at_pointer();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_clear_pending_picks_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   clear_pending_picks();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_delete_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = wrapped_create_delete_item_dialog();
  gtk_widget_show(widget);
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_undo_button_clicked_gtkbuilder_callback   (GtkButton       *button,
                                                            gpointer         user_data)
{
   apply_undo();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_redo_button_clicked_gtkbuilder_callback (GtkButton       *button,
                                                         gpointer         user_data)
{
  apply_redo();
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refmac_button_clicked_gtkbuilder_callback (GtkToolButton   *toolbutton,
                                                            gpointer         user_data)
{
  /* wrapped_create_run_refmac_dialog(); */
  wrapped_create_simple_refmac_dialog();

}


extern "C" G_MODULE_EXPORT
void
on_refine_params_torsion_weight_combobox_changed_gtkbuilder_callback(GtkComboBox     *combobox,
                                                                     gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_torsion_weight_from_text(active_item_idx, t);
}

extern "C" G_MODULE_EXPORT
void
on_refine_params_rama_restraints_combobox_changed_gtkbuilder_callback (GtkComboBox     *combobox,
                                                                       gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_ramachandran_restraints_weight_from_text(active_item_idx, t);
}

// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_strand_rest_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//n                                                                    gpointer         user_data) {
//  if (gtk_toggle_button_get_active(togglebutton))
//    set_secondary_structure_restraints_type(2);
//}

// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_helix_rest_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                    gpointer         user_data) {
//   if (gtk_toggle_button_get_active(togglebutton))
//     set_secondary_structure_restraints_type(1);
// }

// extern "C" G_MODULE_EXPORT
// void
// on_sec_str_rest_no_rest_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                 gpointer         user_data) {
//   if (gtk_toggle_button_get_active(togglebutton))
//     set_secondary_structure_restraints_type(0);
// }

// extern "C" G_MODULE_EXPORT
// void
// on_refine_params_use_torsions_checkbutton_toggled_gtkbuilder_callback(GtkToggleButton *togglebutton,
//                                                                       gpointer         user_data) {
//    do_torsions_toggle(GTK_WIDGET(togglebutton));
// }



extern "C" G_MODULE_EXPORT
void
on_open_coordinates1_activate_gtkbuilder_callback          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  open_coords_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_open_dataset1_activate_gtkbuilder_callback (GtkMenuItem     *menuitem,
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
on_auto_open_mtz_activate_gtkbuilder_callback              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  int is_auto_read_fileselection = 1;
  int is;
  GtkWidget *file_filter_button;
  GtkWidget *sort_button;
  GtkWidget *dataset_fileselection1 = coot_dataset_chooser ();
  add_ccp4i_project_optionmenu(dataset_fileselection1, COOT_DATASET_FILE_SELECTION);
  file_filter_button = add_filename_filter_button(dataset_fileselection1,
						  COOT_DATASET_FILE_SELECTION);
  sort_button = add_sort_button_fileselection(dataset_fileselection1);
  /*   set_directory_for_fileselection(dataset_fileselection1); */

  /* stuff in user data saying if this is autoread or not... */
  is = is_auto_read_fileselection;
  g_object_set_data(G_OBJECT(dataset_fileselection1), "imol", GINT_TO_POINTER(is));
  set_file_selection_dialog_size(dataset_fileselection1);

  set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_fileselection1);
  gtk_widget_show (dataset_fileselection1);
  push_the_buttons_on_fileselection(file_filter_button, sort_button, dataset_fileselection1);
}

extern "C" G_MODULE_EXPORT
void
on_save_coordinates_filechooser_dialog_response_gtkbuilder_callback(GtkDialog       *dialog,
                                                                    gint             response_id,
                                                                    gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {
      const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      if (fnc) {
         std::string fn(fnc);
         std::cout << "Now do something with " << fn << std::endl;
      }
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}

extern "C" G_MODULE_EXPORT
void
on_dataset_filechooser_dialog_response_gtkbuilder_callback(GtkDialog       *dialog,
                                                           gint             response_id,
                                                           gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      if (fnc) {
         std::string fn(fnc);
         manage_column_selector(fnc); // move the function declaration into a c++ header one day
      }
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}

extern "C" G_MODULE_EXPORT
void
on_map_name_filechooser_dialog_response_gtkbuilder_callback(GtkDialog       *dialog,
                                                            gint             response_id,
                                                            gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {
      const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      if (fnc) {
         std::string fn(fnc);
         std::cout << "Now do something with " << fn << std::endl;
      }
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}

extern "C" G_MODULE_EXPORT
void
on_coords_filechooser_dialog_response_gtkbuilder_callback(GtkDialog       *dialog,
                                                          gint             response_id,
                                                          gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {
      const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      if (fnc) {
         std::string fn(fnc);
         std::cout << "Now do something with " << fn << std::endl;
      }
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}

