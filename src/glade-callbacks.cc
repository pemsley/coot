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


#include "Python.h"

#include <iostream>
#include <gtk/gtk.h>

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"
#include "positioned-widgets.h"
#include "interface.h"
#include "coot-references.h"
#include "c-interface-gui.hh"

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
#include "cc-interface.hh"
#include "cfc-widgets-c-interface.h"

#include "graphics-info.h"
#include "validation-graphs/validation-graphs.hh"

#include "widget-from-builder.hh"

#include "fit-loop-gui.hh"
#include "c-interface-refine.hh"
#include "gtkglarea-rama-plot.hh"
#include "cc-interface-scripting.hh"

void get_monomer_dictionary_in_subthread(const std::string &comp_id, bool state);


// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;

extern "C" G_MODULE_EXPORT
gboolean on_about_dialog_close_request(GtkAboutDialog *dialog, gpointer user_data) {
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   return TRUE; // Prevent the default close behavior (destruction)
}

extern "C" G_MODULE_EXPORT
gboolean
on_select_fitting_map_dialog_close_request(GtkAboutDialog *dialog, gpointer user_data) {
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   return TRUE; // Prevent the default close behavior (destruction)
}

extern "C" G_MODULE_EXPORT
void
on_column_label_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *column_label_window = widget_from_builder("column_label_window");
  handle_column_label_make_fourier_v2(column_label_window);
}


extern "C" G_MODULE_EXPORT
void
on_column_label_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   printf("column_label_window CANCEL button clicked\n");
   gtk_widget_set_visible(widget_from_builder("column_label_window"), FALSE);
}



/* Density Radius Window Widgets  */

extern "C" G_MODULE_EXPORT
void
on_density_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data) {

   GtkEntry      *entry;
   const char *text;

   GtkEntry *entry_xray = GTK_ENTRY(widget_from_builder("map_parameters_x_ray_radius_entry"));
   GtkEntry *entry_em   = GTK_ENTRY(widget_from_builder("map_parameters_em_radius_entry"));
   const char *text_xray = gtk_editable_get_text(GTK_EDITABLE(entry_xray));
   const char *text_em   = gtk_editable_get_text(GTK_EDITABLE(entry_em));
   int imol = -1;

   set_density_size_from_widget(text_xray);
   set_density_size_em_from_widget(text_em);

   /* Now the increment of the iso level entry */

   entry = GTK_ENTRY(widget_from_builder("iso_level_increment_entry"));

   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   set_iso_level_increment_from_text(text, imol);  /* imol is ignored. */
   /* we change the iso level for all maps */

   /* As above, except a difference map */

   entry = GTK_ENTRY(widget_from_builder("diff_map_iso_level_increment_entry"));

   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   set_diff_map_iso_level_increment_from_text(text, imol); /* imol is ignored. */
                                       /* we change the iso level for all maps */

   entry = GTK_ENTRY(widget_from_builder("map_sampling_rate_entry"));

   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   set_map_sampling_rate_text(text);

   // what a mish-mash of naming schemes. (it was one of the very first dialogs)
   GtkWidget *dialog = widget_from_builder("global_map_properties_window");
   gtk_widget_set_visible(dialog, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_density_cancel_clicked              (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_set_visible(widget_from_builder("global_map_properties_window"), FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_fps_window_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data) {

   /* Lets look up the buttons and see if they are active. If the yes
      button is active then set the flag to 1, if no is active, set it
      to 0.  */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("radiobutton1"))))
      set_fps_flag(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("radiobutton2"))))
      set_fps_flag(0);

   gtk_widget_set_visible(widget_from_builder("fps_window"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_active_map_ok_button_clicked(GtkButton       *button,
                                gpointer         user_data) {

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("active_map_radiobutton_yes"))))
      set_active_map_drag_flag(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("active_map_radiobutton_no"))))
      set_active_map_drag_flag(0);

   gtk_widget_set_visible(widget_from_builder( "active_map_window"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data) {

   // 20230513-PE now the widget is active - we don't do anything other than hide
   // (the actual button label is "Close")

   gtk_widget_set_visible(widget_from_builder( "show_symmetry_window"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_apply_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  /* Just a copy of the above funtion, without the destroy. */

   GtkEntry      *entry;
   GtkWidget     *checkbutton;
   const char *text;

/* Show Symmetry Radiobuttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("show_symmetry_yes_radiobutton"))))
      set_show_symmetry_master(1);


   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("show_symmetry_no_radiobutton"))))
      set_show_symmetry_master(0);

/* Symmetry Radius Entry */

   entry = (GTK_ENTRY(widget_from_builder(
				    "symmetry_radius_entry")));

   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   /* printf("ok button symmetry: radius text: %s\n",text); */

   set_symmetry_size_from_widget(text);


/* Show UnitCell Radiobuttons */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("unit_cell_yes_radiobutton"))))
      set_show_unit_cells_all(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("unit_cell_no_radiobutton"))))
      set_show_unit_cells_all(0);

/* The Symmetry Colour Checkbutton */

/*    checkbutton = widget_from_builder( */
/* 			       "show_symmetry_molecule_rotate_colour_map_checkbutton"); */

/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_molecule_rotate_colour_map(1); */
/*    else */
/*      set_symmetry_molecule_rotate_colour_map(0); */


/* The Symmetry Colour by Symop Checkbutton */

/*    checkbutton = widget_from_builder( */
/* 			       "show_symmetry_colour_by_symop_checkbutton"); */
/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_colour_by_symop(1); */
/*    else */
/*      set_symmetry_colour_by_symop(0); */


/* The Symmetry Whole Chain Checkbutton */

/*    checkbutton = widget_from_builder( */
/* 			       "show_symmetry_whole_molecule_checkbutton"); */
/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_whole_chain(1); */
/*    else */
/*      set_symmetry_whole_chain(0); */

/* The Expanded Atom Label Checkbutton */

   checkbutton = widget_from_builder("show_symmetry_expanded_labels_checkbutton");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
     set_symmetry_atom_labels_expanded(1);
   else
     set_symmetry_atom_labels_expanded(0);
}



extern "C" G_MODULE_EXPORT
void
on_symmetry_colour_patch_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
   /* GtkWidget *colorseldlg;
   gdouble *colour;
   GtkColorSelection *colorsel;
*/
   // GTK-FIXME fix the colours

/*    colorseldlg = create_symmetry_colour_selection_window(); */

/*    colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(colorseldlg)->colorsel); */
/*    colour = get_symmetry_bonds_colour(0); */
/*    gtk_color_selection_set_color(colorsel, colour); */

/*    gtk_widget_set_visible(colorseldlg, TRUE); */
}



// old code - delete on a rainy day
extern "C" G_MODULE_EXPORT
void
on_show_aniso_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkEntry      *entry;
   const char *text;

/* Limit Display Atoms? */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder(
								    "limit_display_radius_yes_radiobutton"
								    ))))
      set_limit_aniso(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder(
								    "limit_display_radius_no_radiobutton"
								    ))))
      set_limit_aniso(0);


/* Limit Display Radius Entry */

   entry = GTK_ENTRY(widget_from_builder(
				   "limit_display_radius_entry"));

   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   set_aniso_limit_size_from_widget(text);

 /* Show Aniso Radiobuttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder(
								    "show_aniso_yes_radiobutton"))))

      set_show_aniso(1);	/* model, state */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("show_aniso_no_radiobutton"))))
      set_show_aniso(0);
}


extern "C" G_MODULE_EXPORT
void
on_show_aniso_close_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
/* Goodbye Widget */
   gtk_widget_set_visible(widget_from_builder("aniso_window"), FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_smooth_scrolling_window_ok_button_clicked (GtkButton       *button,
                                              gpointer         user_data) {
   GtkEntry *entry;
   const char *text;

/* Show Smooth Scrolling Radio Buttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("smooth_scroll_yes_radiobutton"))))

      set_smooth_scroll_flag(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("smooth_scroll_no_radiobutton"))))

      set_smooth_scroll_flag(0);

/* Smooth Scroll Distance Limit */

   entry = GTK_ENTRY(widget_from_builder("smooth_scroll_limit_entry"));
   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   set_smooth_scroll_limit_str(text);


 /*  Smooth Scroll Steps */

   entry = GTK_ENTRY(widget_from_builder( "smooth_scroll_steps_entry"));
   text = gtk_editable_get_text(GTK_EDITABLE(entry));

   set_smooth_scroll_steps_str(text);

   gtk_widget_set_visible(widget_from_builder("smooth_scroll_window"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_phs_info_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   phs_pdb_cell_symm(); /* which runs
			   create_phs_coordinates_fileselection() */

   gtk_widget_set_visible(widget_from_builder( "phs_info_box"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_phs_info_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_set_visible(widget_from_builder("phs_info_box"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_map_and_mol_control1_activate       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window();
   set_transient_for_main_window(widget);
   gtk_widget_set_visible(widget, TRUE);
}



extern "C" G_MODULE_EXPORT
void
on_go_to_atom_apply_button_clicked (GtkButton       *button,
				    gpointer         user_data)
{

  GtkWidget *widget;

  widget = widget_from_builder("goto_atom_window");

  apply_go_to_atom_from_widget(widget);


}


extern "C" G_MODULE_EXPORT
void
on_go_to_atom_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget;

  widget     = widget_from_builder("goto_atom_window");

  unset_go_to_atom_widget();
  gtk_widget_set_visible(widget, FALSE);
                                /* There is something that had been
				   added to the Go To Atom window that
				   is not a widget.  It has been
				   cleared perhaps but not destroyed.
				   I don't think that it is the dialog
				   itself.  The problem does not
				   happen in the GTK1 path.  */
}



extern "C" G_MODULE_EXPORT
void
on_go_to_atom_button_clicked (GtkButton       *button,
                              gpointer         user_data) {

   GtkWidget *dialog = wrapped_create_goto_atom_window(); // uses gtkbuilder
   gtk_widget_set_visible(dialog, TRUE);
}



extern "C" G_MODULE_EXPORT
void
on_go_to_atom_next_residue_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{

/*   GtkEntry *entry;  */
/*   GtkEntry *residue_entry;  */
/*   gchar *chain_str; */
/*   gchar *res_str;  */
/*   gchar *atom_name_str;  */

/*   entry = GTK_ENTRY(widget_from_builder( */
/* 				  "go_to_atom_chain_entry")); */
/*   chain_str = gtk_entry_get_text(entry); */

/*   residue_entry = GTK_ENTRY(widget_from_builder( */
/* 					  "go_to_atom_residue_entry")); */
/*   res_str = gtk_entry_get_text(residue_entry); */

/*   entry = GTK_ENTRY(widget_from_builder( */
/* 				  "go_to_atom_atom_name_entry")); */
/*   atom_name_str = gtk_entry_get_text(entry); */

/*   goto_next_atom_maybe(chain_str, res_str, atom_name_str, residue_entry); */

  GtkWidget *window = widget_from_builder("goto_atom_window");
  goto_next_atom_maybe_new();
}


extern "C" G_MODULE_EXPORT
void
on_go_to_atom_previous_residue_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{

/*   GtkEntry *entry;  */
/*   GtkEntry *residue_entry;  */
/*   gchar *chain_str; */
/*   gchar *res_str;  */
/*   gchar *atom_name_str;  */

/*   entry = GTK_ENTRY(widget_from_builder( */
/* 				  "go_to_atom_chain_entry")); */
/*   chain_str = gtk_entry_get_text(entry); */

/*   residue_entry = GTK_ENTRY(widget_from_builder( */
/* 					  "go_to_atom_residue_entry")); */
/*   res_str = gtk_entry_get_text(residue_entry); */

/*   entry = GTK_ENTRY(widget_from_builder( */
/* 				  "go_to_atom_atom_name_entry")); */
/*   atom_name_str = gtk_entry_get_text(entry); */
/*   goto_prev_atom_maybe(chain_str, res_str, atom_name_str, residue_entry);  */

  GtkWidget *window = widget_from_builder("goto_atom_window");
  if (! window)
     printf("ERROR:: in on_go_to_atom_previous_residue_button_clicked NULL window\n");
  goto_previous_atom_maybe_new();
}

extern "C" G_MODULE_EXPORT
void
on_go_to_atom_show_waters_checkbutton_toggled(GtkCheckButton *check_button,
                                              gpointer user_data) {

   graphics_info_t g;
   g.fill_go_to_atom_window_residue_and_atom_lists_gtk4();
}

extern "C" G_MODULE_EXPORT
void
on_go_to_atom_show_ligands_only_checkbutton_toggled(GtkCheckButton *check_button,
                                                    gpointer user_data) {

   graphics_info_t g;
   g.fill_go_to_atom_window_residue_and_atom_lists_gtk4();
}


/* void */
/* on_skeleton_box_size_ok_button_clicked (GtkButton       *button, */
/*                                         gpointer         user_data) */
/* { */
/*   GtkEntry *entry;  */
/* } */


extern "C" G_MODULE_EXPORT
void
on_skel_box_radius_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   const char    *txt;
   GtkEntry *entry = GTK_ENTRY(widget_from_builder("skeleton_box_size_entry"));
   txt = gtk_editable_get_text(GTK_EDITABLE(entry));
   set_skeleton_box_size_from_widget(txt);
   gtk_widget_set_visible(widget_from_builder("skeletonization_box_radius_window"), FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_skel_box_radius_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{

   gtk_widget_set_visible(widget_from_builder("skeletonization_box_radius_window"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_skeletonization_level_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{

  const char *txt;

  GtkEntry *entry = GTK_ENTRY(widget_from_builder("skeleton_level_entry"));

  txt = gtk_editable_get_text(GTK_EDITABLE(entry));

  set_skeletonization_level_from_widget(txt); /* does a skeleton
						 update and redraw */

  gtk_widget_set_visible(widget_from_builder("skeletonization_level_window"), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_skeletonization_level_apply_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{


  const char *txt;

  GtkEntry *entry = GTK_ENTRY(widget_from_builder("skeleton_level_entry"));

  txt = gtk_editable_get_text(GTK_EDITABLE(entry));

  set_skeletonization_level_from_widget(txt); /* does a skeleton
						 update and redraw */


}


extern "C" G_MODULE_EXPORT
void
on_skeletonization_level_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data)
{

  gtk_widget_set_visible(widget_from_builder("skeletonization_level_window"), FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_display_control_ok_button_clicked   (GtkButton       *button,
                                                            gpointer         user_data) {

   reset_graphics_display_control_window(); /* Needed! (also resets the scroll group) */
   GtkWidget *w = widget_from_builder("display_control_window_glade");
   if (w)
      gtk_widget_set_visible(w, FALSE);
   graphics_info_t::graphics_grab_focus();
}

extern "C" G_MODULE_EXPORT
void
on_rotation_centre_size_ok_button_clicked(GtkButton       *button,
                                          gpointer         user_data) {

   /* pass back the value from the entry */
  GtkEntry *entry = GTK_ENTRY(widget_from_builder("rotation_centre_cube_size_entry"));
  const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  set_rotation_centre_size_from_widget(text);
  gtk_widget_set_visible(widget_from_builder("rotation_centre_cube_size_window"), FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_phs_cell_choice_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   handle_phs_cell_choice_ok_button_clicked(GTK_WIDGET(button));
}


extern "C" G_MODULE_EXPORT
void
on_phs_cell_choice_cancel_button_clicked (GtkButton       *button,
                                                              gpointer         user_data)
{
   printf("Cancel trying to read a phs file.\n");
   GtkWidget *window = widget_from_builder("phs_cell_choice_window");
   gtk_widget_set_visible(window, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_fetch_pdb_and_map_using_pdb_redo1_activate(GMenuItem     *menuitem,
                                              gpointer         user_data) {

   int n = COOT_ACCESSION_CODE_WINDOW_PDB_REDO;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_set_visible(window, TRUE);

}

extern "C" G_MODULE_EXPORT
void
on_fetch_corresponding_alphafold_model_and_superpose1_activate
                                        (GMenuItem     *menuitem,
					 gpointer       user_data) {

  fetch_and_superpose_alphafold_models_using_active_molecule();
}




extern "C" G_MODULE_EXPORT
void
on_accession_code_get_it_button_clicked(GtkButton *button, gpointer user_data) {

   GtkWidget *entry = widget_from_builder("accession_code_entry");
   GtkWidget *frame = widget_from_builder("accession_code_frame");
   // 20230731-PE don't do this here, you dozy pillock.
   // int n = COOT_ACCESSION_CODE_WINDOW_OCA;
   // g_object_set_data(G_OBJECT(frame), "mode", GINT_TO_POINTER(n));
   handle_get_accession_code(frame, entry);
}

extern "C" G_MODULE_EXPORT
void
on_emdb_map_code_get_it_button_clicked(GtkButton *button, gpointer user_data) {

   GtkWidget *entry = widget_from_builder("emdb_map_code_entry");
   GtkWidget *frame = widget_from_builder("emdb_map_code_frame");
   int n = COOT_EMDB_CODE;
   g_object_set_data(G_OBJECT(frame), "mode", GINT_TO_POINTER(n));
   handle_get_accession_code(frame, entry);
}

extern "C" G_MODULE_EXPORT
void
on_dynarama_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				     "dynarama_window");
   int imol = get_mol_from_dynarama(window); /*  return -9999 on edit rama window */

/* No need to free dynarama plot data here, because it is done on
   window destroy callback */
/*     if (imol >= 0) */
/*       set_dynarama_is_displayed(0, imol); */ /*  which frees/deletes the */
     /*  			             memory of the user data. */


   if (imol == -9999)
     accept_phi_psi_moving_atoms();

   gtk_widget_set_visible(window, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_dynarama_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				     "dynarama_window");
   int imol = get_mol_from_dynarama(window); // return -9999 on edit rama window

/*    printf("get_mol_from_dynarama return %d\n", imol); */

/* No need to free dynarama plot data here, because it is done on
   window destroy callback */
/*     if (imol >= 0) */
/*        set_dynarama_is_displayed(0, imol);  which frees/deletes the
                                               memory of the user data. */

   if (imol == -9999)
     clear_moving_atoms_object();

   gtk_widget_set_visible(window, FALSE);

   /* maybe we should also destroy the edit backbone torsion dialog, if it exists. */

}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                                            gpointer         user_data) {

   // execute_get_mols_ligand_search() no longer returns the number of ligands
   execute_get_mols_ligand_search(GTK_WIDGET(button)); /* which then runs execute_ligand_search */
   GtkWidget *window = widget_from_builder("find_ligand_dialog");
   // free_ligand_search_user_data(GTK_WIDGET(button)); // not if not destroyed? Needs checking.
   gtk_widget_set_visible(window, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder("find_ligand_dialog");
   free_ligand_search_user_data(GTK_WIDGET(button));
   gtk_widget_set_visible(window, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_many_atoms_continue_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   GtkWidget *window = widget_from_builder("find_ligand_many_atoms_dialog");
   GtkWidget *find_ligand_dialog = widget_from_builder("find_ligand_dialog");

   execute_ligand_search();
   gtk_widget_set_visible(window, FALSE);
   gtk_widget_set_visible(find_ligand_dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_many_atoms_cancel_button_clicked(GtkButton       *button,
                                                                    gpointer         user_data) {

   GtkWidget *window = widget_from_builder("find_ligand_many_atoms_dialog");
   gtk_widget_set_visible(window, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_dismiss_button_clicked (GtkButton       *button,
					       gpointer         user_data) {

   GtkWidget *hbox   = widget_from_builder("model_fit_refine_dialog_vbox");
   GtkWidget *dialog = widget_from_builder("model_refine_dialog");

   clear_pending_picks();
   normal_cursor();
   dialog = close_model_fit_dialog(hbox);
   /* was it a top-level?  If so, kill off the top-level. */
   if (dialog)
      gtk_widget_set_visible(dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void on_save_coords_filechooser_dialog_response(GtkDialog *dialog,
                                                int        response) {

   if (response == GTK_RESPONSE_YES) { // maybe not the right one, but it is the one set in
                                       // on_save_coords_dialog_save_button_clicked()
      GtkFileChooser *chooser = GTK_FILE_CHOOSER (dialog);
      GFile *file   = gtk_file_chooser_get_file(chooser);
      char *file_name = g_file_get_path(file);
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "imol"));
      graphics_info_t g;
      if (g.is_valid_model_molecule(imol)) {
         bool save_hydrogens = true;
         bool save_conect_records = true;
         bool save_aniso_records = true;
         g.molecules[imol].save_coordinates(file_name, save_hydrogens, save_aniso_records, save_conect_records);
      }
   }

   // destroy this?
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_save_coords_cancel_clicked(G_GNUC_UNUSED GtkButton       *button,
                              G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *mol_selector_frame = widget_from_builder("save_coords_frame");
   gtk_widget_set_visible(mol_selector_frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_save_coords_save_button_clicked(G_GNUC_UNUSED GtkButton       *button,
                                   G_GNUC_UNUSED gpointer         user_data) {

   // we need to select the molecule to save - this is someone clicking on the
   // "Save Molecule" button in the save molecule chooser - not in a file selector

   GtkWidget *combobox = widget_from_builder("save_coordinates_combobox");
   GtkWidget *mol_selector_frame = widget_from_builder("save_coords_frame");
   if (! combobox) {
      std::cout << "ERROR:: on_save_coords_dialog_save_button_clicked: bad combobox\n";
   } else {

      int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
      GtkWindow *parent_window = GTK_WINDOW(graphics_info_t::get_main_window());
      GtkFileChooserAction action = GTK_FILE_CHOOSER_ACTION_SAVE;

      // gtk_file_chooser_dialog_new() is deprecrated!
      //
      GtkWidget *file_chooser_dialog = gtk_file_chooser_dialog_new("Save Coordinates",
                                                                   parent_window,
                                                                   action,
                                                                   ("_Cancel"),
                                                                   GTK_RESPONSE_CANCEL,
                                                                   ("_Save"),
                                                                   GTK_RESPONSE_YES,
                                                                   NULL);

      // 20230910-PE imol is used in on_save_coords_filechooser_dialog_response().
      // Maybe I could/should use user-data instead of get/set data on the widget...
      //
      g_object_set_data(G_OBJECT(file_chooser_dialog), "imol", GINT_TO_POINTER(imol));
      g_signal_connect(file_chooser_dialog, "response",
                       G_CALLBACK(on_save_coords_filechooser_dialog_response), NULL);
      gtk_widget_set_visible(file_chooser_dialog, TRUE);
      set_file_for_save_filechooser(file_chooser_dialog);
      add_filename_filter_button(file_chooser_dialog, COOT_SAVE_COORDS_FILE_SELECTION);

   }
   gtk_widget_set_visible(mol_selector_frame, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_refine_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    do_refine(1);
  else
    do_refine(0);		/* unclick button */

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_regularize_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    do_regularize(1);
  else
    do_regularize(0);		/* unclick button */

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_fixed_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_fixed_atom_dialog();
  gtk_widget_set_visible(w, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_rigid_body_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
    printf("Rigid Body:\n");
    do_rigid_body_refine(1);
  } else {
     do_rigid_body_refine(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_refine_params_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {

   // 20211026-PE I don't think that this is used. Call back is in glade-callbacks-main-window.cc
   // GtkWidget *widget = wrapped_create_refine_params_dialog();
   // GtkWidget *widget = widget_from_builder("refinement_and_regularization_parameters_dialog");
   // gtk_widget_set_visible(widget, TRUE);
}

extern "C" G_MODULE_EXPORT
void
on_refinement_and_regularization_vbox_close_button_clicked(GtkButton       *button,
                                                           gpointer         user_data) {

   GtkWidget *frame = widget_from_builder("refinement_and_regularization_parameters_frame");
   gtk_widget_set_visible(frame,FALSE);

   // pressing the button means that the focus goes elsewhere (not sure where). So bring it back to the graphics
   // widget;

   GtkWidget *glarea = graphics_info_t::glareas[0];
   gtk_widget_grab_focus(glarea);
}

extern "C" G_MODULE_EXPORT
void
on_refine_params_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder(
				    "refine_params_dialog");
  GtkWidget *entry = widget_from_builder(
				   "refine_params_weight_matrix_entry");
  if (entry) {
    set_refinement_weight_from_entry(entry);
  }

  gtk_widget_set_visible(widget, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_goto_atom_window_destroy            (GtkWidget       *object,
                                        gpointer         user_data)
{
/*   printf("on_goto_atom_window_destroy is executed. unsetting the static in graphics_info\n"); */

/* we can get here from a WM dialog kill (as well as the conventional
   "Cancel" button callback. Fixes July 12 bug? */
  unset_go_to_atom_widget();
}



extern "C" G_MODULE_EXPORT
void
on_find_waters_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   execute_find_waters();
   GtkWidget *widget = widget_from_builder("find_waters_dialog");
   gtk_widget_set_visible(widget, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_find_waters_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget;

   widget = widget_from_builder(
			  "find_waters_dialog");

   gtk_widget_set_visible(widget, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_fast_sss_dialog_citation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  GtkWidget *toolbutton = widget_from_builder("coot_references_buccaneer_toolbutton");
  fill_references_notebook(GTK_BUTTON(toolbutton), COOT_REFERENCE_BUCCANEER);

}


/* void */
/* on_enviroment_distance_checkbutton_toggled (GtkToggleButton *togglebutton, */
/* 					    gpointer         user_data) */
/* { */
/*    toggle_environment_show_distances(); */
/* } */


/* void */
/* on_enviroment_distance_dialog_ok_button_clicked (GtkButton       *button, */
/* 						 gpointer         user_data) */
/* { */
/*    GtkWidget *widget; */
/*    widget = widget_from_builder("environment_distance_dialog"); */
/*    execute_environment_settings(GTK_WIDGET(button)); */
/*    gtk_widget_destroy(widget); */

/* } */


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_planar_peptides_checkbutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data)
{
   if (gtk_check_button_get_active(checkbutton)) {
      add_planar_peptide_restraints();
   } else {
      remove_planar_peptide_restraints();
   }
}



extern "C" G_MODULE_EXPORT
void
on_refine_params_use_trans_peptide_restraints_checkbutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data) {

   if (gtk_check_button_get_active(checkbutton)) {
      std::cout << "debug:: in on_refine_params_use_trans_peptide_restraints_checkbutton_toggled() "
                << " active" << std::endl;
      set_use_trans_peptide_restraints(1);
   } else {
      std::cout << "debug:: in on_refine_params_use_trans_peptide_restraints_checkbutton_toggled() "
                << " inactive" << std::endl;
      set_use_trans_peptide_restraints(0);
   }
}




extern "C" G_MODULE_EXPORT
void
on_map_dynamic_map_sampling_checkbutton_toggled (GtkToggleButton *togglebutton,
						 gpointer         user_data)
{
   toggle_dynamic_map_sampling();
}


extern "C" G_MODULE_EXPORT
void
on_map_dynamic_map_size_display_checkbutton_toggled (GtkToggleButton *togglebutton,
						     gpointer         user_data)
{
   toggle_dynamic_map_display_size();
}


extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *fileselection;

  fileselection = widget_from_builder("cif_dictionary_fileselection");
  gtk_widget_set_visible(fileselection, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_fit_terminal_residue_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    do_add_terminal_residue(1);
  else
    do_add_terminal_residue(0);
}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_ALA_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;

   gtk_widget_set_visible(window, FALSE);
   do_mutation("ALA", istate);
}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_ARG_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("ARG", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_ASN_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("ASN", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_ASP_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("ASP", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_CYS_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("CYS", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_GLN_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("GLN", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_GLU_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("GLU", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_GLY_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("GLY", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_HIS_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("HIS", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_ILE_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("ILE", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_LEU_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("LEU", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_LYS_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("LYS", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_MET_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("MET", istate);

}

extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_MSE_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("MSE", istate);

}



extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_PHE_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("PHE", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_PRO_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("PRO", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_SER_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("SER", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_THR_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("THR", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_TRP_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("TRP", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_TYR_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("TYR", istate);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_VAL_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *window = widget_from_builder(
				    "residue_type_chooser_window");
   GtkWidget *stub_button =
     widget_from_builder("residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_set_visible(window, FALSE);
   do_mutation("VAL", istate);

}




extern "C" G_MODULE_EXPORT
void
on_rotate_translate_obj_ok_button_clicked (GtkButton       *button,
                                           gpointer         user_data) {
  GtkWidget *widget = widget_from_builder("rotate_translate_obj_dialog");

  // graphics_unsetup_rotate_translate_buttons(widget);
  rot_trans_reset_previous();
  accept_regularizement();
  gtk_widget_set_visible(widget, FALSE);
  clear_up_moving_atoms(); // redraw done here

}


extern "C" G_MODULE_EXPORT
void
on_rotate_translate_obj_cancel_button_clicked (GtkButton       *button,
                                               gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder("rotate_translate_obj_dialog");

 /*   clear_moving_atoms_object(); // redraw done here */
  std::cout << "clear_up_moving_atoms()" << std::endl;
  clear_up_moving_atoms();
  std::cout << "gtk_widget_set_visible()" << std::endl;
  gtk_widget_set_visible(widget, FALSE);
  std::cout << "normal_cursor()" << std::endl;
  normal_cursor();
}





extern "C" G_MODULE_EXPORT
void
on_delete_item_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
					    gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("delete_item_residue_radiobutton"))))
    set_delete_residue_mode();

}


extern "C" G_MODULE_EXPORT
void
on_delete_item_atom_radiobutton_toggled (GtkToggleButton *togglebutton,
					 gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("delete_item_atom_radiobutton"))))
    set_delete_atom_mode();
}



extern "C" G_MODULE_EXPORT
void
on_delete_item_sidechain_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_sidechain_mode();
}

extern "C" G_MODULE_EXPORT
void
on_delete_item_sidechain_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data) {

  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_sidechain_range_mode();
}



extern "C" G_MODULE_EXPORT
void
on_delete_item_chain_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_chain_mode();
}




extern "C" G_MODULE_EXPORT
void
on_delete_item_residue_hydrogens_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_residue_hydrogens_mode();

}


extern "C" G_MODULE_EXPORT
void
on_delete_item_water_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_delete_water_mode();
}


/* Old unused function, see new_close_molecules_dialog() */
extern "C" G_MODULE_EXPORT
void
on_close_molecule_close_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "close_molecule_dialog");
  GtkWidget *optionmenu = widget_from_builder("close_molecule_optionmenu");
  close_molecule_by_widget(optionmenu);
/*   gtk_widget_set_visible(window, FALSE); */
}


extern "C" G_MODULE_EXPORT
void
on_close_molecule_cancel_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "close_molecule_dialog");
  gtk_widget_set_visible(window, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_delete_item_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data) {

   // 20220602-PE what does this do these days?
   GtkWidget *widget = widget_from_builder("delete_item_dialog");
   clear_pending_delete_item();
   clear_pending_picks(); 	/* hmmm.. not sure 20050610 */
   normal_cursor();
   // store_delete_item_widget(NULL);
   gtk_widget_set_visible(widget, FALSE);
}


/* this is the Apply button now */
extern "C" G_MODULE_EXPORT
void
on_residue_info_apply_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {

   apply_residue_info_changes();
   // GtkWidget *widget = widget_from_builder("residue_info_dialog");
   // gtk_widget_set_visible(widget, FALSE);
   add_status_bar_text("Occupancies and B-factors have been updated");
   GtkWidget *label = widget_from_builder("occupancy_b_factors_updated_label");
   if (label) {
      gtk_widget_set_visible(label, TRUE);

      auto label_callback = +[] (gpointer user_data) {
         GtkWidget *w = GTK_WIDGET(user_data);
         gtk_widget_set_visible(w, FALSE);
         return 0;
      };
      g_timeout_add(1000, G_SOURCE_FUNC(label_callback), label);

   }
}


extern "C" G_MODULE_EXPORT
void
on_residue_info_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *widget = widget_from_builder("residue_info_dialog");

   residue_info_release_memory(widget);  // Hmmm! that seems dangerous
   gtk_widget_set_visible(widget, FALSE);
   unset_residue_info_widget(); // 20230515-PE seems an ancient thing to do. Needed?
}

extern "C" G_MODULE_EXPORT
void
on_residue_info_master_atom_occ_entry_changed (GtkEntry     *entry,
                                               gpointer      user_data) {

   const char *txt = gtk_editable_get_text(GTK_EDITABLE(entry));

   if (txt) {
      std::string s(txt);
      std::cout << "master atom occ changed to " << s << std::endl;
      graphics_info_t g;
      g.residue_info_edit_occ_apply_to_other_entries_maybe(GTK_WIDGET(entry));
   }
}


extern "C" G_MODULE_EXPORT
void
on_residue_info_master_atom_occ_entry_activate(GtkWidget *entry, gpointer user_data) {

}

extern "C" G_MODULE_EXPORT
void
on_residue_info_master_atom_b_factor_entry_changed(GtkEntry     *entry,
                                                   gpointer      user_data) {

   const char *txt = gtk_editable_get_text(GTK_EDITABLE(entry));

   if (txt) {
      std::string s(txt);
      std::cout << "master atom B-factor changed to " << s << std::endl;
      graphics_info_t g;
      g.residue_info_edit_b_factor_apply_to_other_entries_maybe(GTK_WIDGET(entry));
   }
}

extern "C" G_MODULE_EXPORT
void
on_residue_info_master_atom_b_factor_entry_activate(GtkWidget *entry, gpointer user_data) {

}

extern "C" G_MODULE_EXPORT
void
on_keyboard_mutate_entry_changed(GtkEntry     *entry,
                                 gpointer      user_data) {
}

extern "C" G_MODULE_EXPORT
void
on_keyboard_mutate_entry_activate(GtkWidget *entry, gpointer user_data) {

   std::string s(gtk_editable_get_text(GTK_EDITABLE(entry)));
   mutate_active_residue_to_single_letter_code(s);
   GtkWidget *frame = widget_from_builder("keyboard_mutate_frame");
   gtk_widget_set_visible(frame, FALSE);
   graphics_info_t::graphics_grab_focus();
}

extern "C" G_MODULE_EXPORT
void
on_hints_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("hints_dialog");
   gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_find_ligands_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  do_find_ligands_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_rotamer_selection_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *dialog = widget_from_builder(
				    "rotamer_selection_dialog");
   accept_regularizement();
   clear_moving_atoms_object();
   gtk_widget_set_visible(dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_rotamer_selection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("rotamer_selection_dialog");
   int type;
   int imol;

   clear_up_moving_atoms();
   type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "type"));

   /* if this is a add-alt conf dialog, undo (the addition of the alt
      conf) on Cancel click */
   if (type == 1) {
     imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "imol"));
     if (is_valid_model_molecule(imol)) {
       set_undo_molecule(imol);
       apply_undo();
     }
   }
   gtk_widget_set_visible(dialog, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_rotamer_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
      setup_rotamers(1);
   else
      setup_rotamers(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_mutate_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
      setup_mutate(1);
   else
      setup_mutate(0);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_pointer_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   place_atom_at_pointer();
}


extern "C" G_MODULE_EXPORT
void
on_unsaved_changes_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *dialog = widget_from_builder("unsaved_changes_dialog");
   gtk_widget_set_visible(dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_unsaved_changes_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("unsaved_changes_dialog");
   gtk_widget_set_visible(dialog, FALSE);
   coot_clear_backup_or_real_exit(0);
}


extern "C" G_MODULE_EXPORT
void
on_environment_distance_label_atom_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                       gpointer        user_data) {

  if (gtk_check_button_get_active(checkbutton)) {
     set_environment_distances_label_atom(1);
     graphics_info_t g;
     g.try_label_unlabel_active_atom();
  } else {
     set_environment_distances_label_atom(0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_baton_accept_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  accept_baton_position();
}


extern "C" G_MODULE_EXPORT
void
on_baton_try_again_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  baton_tip_try_another();
}


extern "C" G_MODULE_EXPORT
void
on_baton_lengthen_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  lengthen_baton();
}

extern "C" G_MODULE_EXPORT
void
on_baton_shorten_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
  shorten_baton();
}

/* dismiss button */
extern "C" G_MODULE_EXPORT
void
on_baton_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder(
				    "baton_dialog");
  try_set_draw_baton(0);
  set_baton_mode(0);		/* if you can't see it, there's no
				   point in trying to move it */
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_baton_button_clicked (GtkButton       *button,
                                                                 gpointer         user_data)
{
  /* return whether the baton was drawn or not, We don't want a baton
     dialog if the baton is not drawn (e.g. when there is no
     skeletonized map). */
  int state = try_set_draw_baton(1);
  if (state) {
     // GtkWidget *widget = create_baton_dialog();
     GtkWidget *widget = widget_from_builder("baton_dialog()");
    gtk_widget_set_visible(widget, TRUE);
  }
}




extern "C" G_MODULE_EXPORT
void
on_use_weights_checkbutton_toggled     (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = widget_from_builder("column_label_window_weights_hbox");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      gtk_widget_set_sensitive(GTK_WIDGET(hbox), TRUE);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(hbox), FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_environment_distance_checkbutton_toggled (GtkCheckButton *checkbutton,
                                             gpointer         user_data) {
   toggle_environment_show_distances(checkbutton);
}


extern "C" G_MODULE_EXPORT
void
on_environment_distance_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("environment_distance_dialog");
   execute_environment_settings(GTK_WIDGET(button));
   gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_read_pdb_recentre_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    set_recentre_on_read_pdb(1);
   } else {
    set_recentre_on_read_pdb(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_read_pdb_recentre_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
     set_recentre_on_read_pdb(0);
  } else {
     set_recentre_on_read_pdb(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_read_pdb_recentre_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *widget = widget_from_builder("read_pdb_recentre_dialog");
   // gtk_widget_set_visible(widget, FALSE);

   GtkWidget *w = widget_from_builder("read_pdb_recentre_dialog");
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_pointer_atom_type_cancel_button_clicked(GtkButton       *button,
                                                               gpointer         user_data) {

   // GtkWidget *widget = widget_from_builder("pointer_atom_type_dialog");
   // gtk_widget_set_visible(widget, FALSE);

   GtkWidget *w = widget_from_builder("pointer_atom_type_dialog");
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_pointer_atom_type_ok_button_clicked (GtkButton       *button,
                                                            gpointer         user_data) {

   // GtkWidget *dialog = widget_from_builder("pointer_atom_type_dialog");

  GtkWidget *dialog = widget_from_builder("pointer_atom_type_dialog");

  GtkWidget *entry = widget_from_builder("pointer_atom_type_other_entry");
  const char *entry_text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));

  if (strlen(entry_text) > 0) {
     place_typed_atom_at_pointer(entry_text);
  } else {
     /* Adding something here?
        Remember to change also molecule_class_info_t::add_typed_pointer_atom(). */

     GtkToggleButton *tbut;
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_water"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Water");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_ca"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Ca");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_mg"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Mg");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_na"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Na");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_cl"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Cl");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_br"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Br");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_so4"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("SO4");
     tbut = GTK_TOGGLE_BUTTON(widget_from_builder("pointer_atom_type_radiobutton_po4"));
     if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("PO4");
  }
  /* Recall that the molecule is set by the callback from menu item "activate" */

  // gtk_widget_set_visible(dialog, FALSE);
  gtk_widget_set_visible(dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_refmac_column_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *frame = widget_from_builder("column_label_refmac_frame");

   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(frame, TRUE);
   } else {
      gtk_widget_set_sensitive(frame, FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_run_button_clicked       (GtkButton       *button,
                                                            gpointer         user_data) {
   // use the simple refmac interface.
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_cancel_button_clicked    (GtkButton       *button,
                                                            gpointer         user_data) {

   // not used
}



extern "C" G_MODULE_EXPORT
void
on_run_refmac_map_mtz_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // not used
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_mtz_filechooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   // not usedl
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_mtz_filechooserdialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
   // not used
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_mtz_filechooser_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   // not used
}



extern "C" G_MODULE_EXPORT
void
on_run_refmac_file_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   // not used
}



extern "C" G_MODULE_EXPORT
void
on_run_refmac_sad_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   // not used
}


/* Actually, we only are interested in the state of this when the
   "Run" button is pressed */
extern "C" G_MODULE_EXPORT
void
on_run_refmac_diff_map_checkbutton_toggled (GtkToggleButton *togglebutton,
					    gpointer         user_data)
{

}


/* Actually, we only are interested in the state of this when the
   "Run" button is pressed */
/* not any more! We should show the phase column options depending on what exists */
extern "C" G_MODULE_EXPORT
void
on_run_refmac_phase_combine_checkbutton_toggled (GtkToggleButton *togglebutton,
						 gpointer         user_data)
{
   // not used
}


extern "C" G_MODULE_EXPORT
void
on_baton_undo_button_clicked(GtkButton       *button,
                             gpointer         user_data) {

   baton_build_delete_last_residue();
}



extern "C" G_MODULE_EXPORT
void
on_undo_molecule_chooser_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data) {

   GtkWidget *widget   = widget_from_builder("undo_molecule_chooser_dialog");
   GtkWidget *combobox = widget_from_builder("undo_molecule_chooser_comboboxtext");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol))
      set_undo_molecule(imol);
   gtk_widget_set_visible(widget, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_skeleton_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *window = widget_from_builder("skeleton_dialog");
   // GtkWidget *combobox   = widget_from_builder("skeleton_map_combobox");
  GtkWidget *window   = widget_from_builder("skeleton_dialog");
  GtkWidget *combobox = widget_from_builder("skeleton_map_combobox");
  int do_baton_mode = GPOINTER_TO_INT(user_data);
  /*
  GtkWidget *optionmenu = widget_from_builder("skeleton_map_optionmenu");
  skeletonize_map_by_optionmenu(optionmenu);
  */
  skeletonize_map_by_combobox(combobox);
  gtk_widget_set_visible(window, FALSE);
  if (do_baton_mode)
    try_set_draw_baton(1);

}


extern "C" G_MODULE_EXPORT
void
on_skeleton_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *window = widget_from_builder("skeleton_dialog");
  GtkWidget *window = widget_from_builder("skeleton_dialog");
  gtk_widget_set_visible(window, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_residue_info_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   residue_info_apply_all_checkbutton_toggled();
}



extern "C" G_MODULE_EXPORT
void
on_add_alt_conf_ca_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      set_add_alt_conf_split_type_number(0);
  }
}

extern "C" G_MODULE_EXPORT
void
on_add_alt_conf_whole_single_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
							  gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      set_add_alt_conf_split_type_number(1);
  }
}

extern "C" G_MODULE_EXPORT
void
on_add_alt_conf_residue_range_radiobutton_toggled (GtkToggleButton *togglebutton,
						   gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      set_add_alt_conf_split_type_number(2);
  }
}

extern "C" G_MODULE_EXPORT
void
on_add_alt_conf_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog = widget_from_builder("add_alt_conf_dialog");
  unset_add_alt_conf_dialog();
  unset_add_alt_conf_define();
  gtk_widget_set_visible(dialog, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_validation_dialog_next_button_clicked (GtkButton       *button,
					  gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_validation_dialog_prev_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_validation_dialog_cancel_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("");

  gtk_widget_set_visible(widget, FALSE);

}




extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_reverse_fragment_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_edit_chi_angles_reverse_fragment_state(1);
  else
    set_edit_chi_angles_reverse_fragment_state(0);
}




extern "C" G_MODULE_EXPORT
void
on_run_state_file_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *dialog = widget_from_builder("run_state_file_dialog");    // old glade style
  GtkWidget *dialog = widget_from_builder("run_state_file_dialog");
  gtk_widget_set_visible(dialog, FALSE);
  run_state_file();
}


extern "C" G_MODULE_EXPORT
void
on_run_state_file_cancel_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
   // GtkWidget *dialog = widget_from_builder("run_state_file_dialog"); // old glade style
  GtkWidget *dialog = widget_from_builder("run_state_file_dialog");
  gtk_widget_set_visible(dialog, FALSE);
  gtk_widget_set_visible(dialog, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsions_dialog_destroy(GtkWidget       *object,
                                        gpointer         user_data) {

   clear_moving_atoms_object();
  /* FIXME: also clear out the edib backbone ramaplot, if it exists. */
  /*   destroy_edit_backbone_rama_plot(); */
}


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsion_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "edit_backbone_torsions_dialog");
  accept_regularizement();	/* does a clear too. */
  destroy_edit_backbone_rama_plot();
  gtk_widget_set_visible(window, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsion_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder("edit_backbone_torsions_dialog");
   /*   clear_moving_atoms_object(); done as part of window destroy
        callback */
   destroy_edit_backbone_rama_plot(); // 20211006-PE this function name should be changed
   //gtk_widget_set_visible(window, FALSE);
   gtk_widget_set_visible(window, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_workflow_radiobutton_coords_only_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_workflow_radiobutton_coords_and_dataset_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_workflow_radiobutton_coords_and_map_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_workflow_radiobutton_map_from_dataset_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_workstate_radiobutton_clean_slate_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_workflow_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{

}

extern "C" G_MODULE_EXPORT
void
on_workflow_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_select_map_for_fitting_cancel_button_clicked(GtkButton       *button,  // OK button
                                                gpointer         user_data) {

   GtkWidget *frame = widget_from_builder( "select_map_for_fitting_frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_select_map_for_fitting_button_clicked(GtkButton       *button,  // OK button
                                         gpointer         user_data) {

   GtkWidget *frame       = widget_from_builder( "select_map_for_fitting_frame");
   GtkWidget *weight_entry = widget_from_builder("select_fitting_map_dialog_weight_entry");

   if (weight_entry) {
      std::string t = gtk_editable_get_text(GTK_EDITABLE(weight_entry));
      float f = coot::util::string_to_float(t);
      graphics_info_t g;
      g.geometry_vs_map_weight = f;
   }
   gtk_widget_set_visible(frame, FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_map_select_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   show_select_map_frame();
}


extern "C" G_MODULE_EXPORT
void
on_select_fitting_map_dialog_estimate_button_clicked (GMenuItem     *menuitem,
                                                                          gpointer         user_data)
{
   std::cout << "on_select_fitting_map_dialog_estimate_button_clicked()" << std::endl;
}



extern "C" G_MODULE_EXPORT
void
on_add_alt_conf_dialog_destroy         (GtkWidget       *object,
                                        gpointer         user_data)
{

  unset_add_alt_conf_define();
  unset_add_alt_conf_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_help_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  /* GtkWidget *widget = create_run_refmac_help_dialog();
     we use the 'new' help now, but preserve the old dialog
     as backup for now */
   // GtkWidget *widget = create_run_refmac_nolabels_help_dialog();
   GtkWidget *widget = widget_from_builder("run_refmac_nolabels_help_dialog");
   gtk_widget_set_visible(widget, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("run_refmac_help_dialog");
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_nolabels_help_button_clicked
					(GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *widget = create_run_refmac_nolabels_help_dialog();
   GtkWidget *widget = widget_from_builder("run_refmac_nolabels_help_dialog");
   gtk_widget_set_visible(widget, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_nolabels_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("run_refmac_nolabels_help_dialog");
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_no_restraints_info_dialog_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = widget_from_builder("no_restraints_info_dialog");
  gtk_widget_set_visible(window, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_no_cif_dictionary_bonds_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = widget_from_builder("no_cif_dictionary_bonds_dialog");
  gtk_widget_set_visible(window, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_pointer_atom_type_other_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *entry = widget_from_builder("pointer_atom_type_other_entry");
  gtk_widget_set_sensitive(entry, TRUE);

}


extern "C" G_MODULE_EXPORT
void
on_ligand_big_blob_dismiss_button_clicked(GtkButton       *button,
                                          gpointer         user_data) {

  GtkWidget *window = widget_from_builder( "ligand_big_blob_dialog");
  gtk_widget_set_visible(window, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  GtkWidget *widget = widget_from_builder("edit_chi_angles_dialog");
  accept_regularizement();
  unset_moving_atom_move_chis();
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_cancel_button_clicked(GtkButton       *button,
                                        gpointer         user_data) {

  GtkWidget *widget = widget_from_builder("edit_chi_angles_dialog");
  clear_up_moving_atoms();	/* and remove the graphics object */
  unset_moving_atom_move_chis();
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_ok_button_clicked(GtkButton       *button,
                                  gpointer         user_data) {

   GtkWidget *widget = widget_from_builder("check_waters_dialog");
   do_check_waters_by_widget(widget);
   gtk_widget_set_visible(widget, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_check_waters_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("check_waters_dialog");
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_normal_rotation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  set_graphics_edit_current_chi(0);

}


extern "C" G_MODULE_EXPORT
void
on_geometry_distance_togglebutton_toggled(GtkToggleButton *togglebutton,
                                          gpointer         user_data) {
  if (gtk_toggle_button_get_active(togglebutton))
    do_distance_define();

}




extern "C" G_MODULE_EXPORT
void
on_geometry_clear_last_distance_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   std::cout << "debug:: in on_geometry_clear_last_distance_button_clicked()" << std::endl;
   clear_last_measure_distance();

}


extern "C" G_MODULE_EXPORT
void
on_geometry_clear_all_distances_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  clear_measure_distances();
}


extern "C" G_MODULE_EXPORT
void
on_geometry_clear_atom_labels_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  remove_all_atom_labels();
}


extern "C" G_MODULE_EXPORT
void
on_geometry_dialog_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

#if 0 // 20211006-PE save for reference
  GtkWidget *dialog = widget_from_builder("geometry_dialog");
  /* should we clear geometry on close dialog?  Currently, I think not. */
  /* it is the COOT_DISTANCES_ANGLES_WINDOW, hmm. */
  store_geometry_dialog(NULL);
  gtk_widget_set_visible(dialog, FALSE);
#endif

  GtkWidget *frame = widget_from_builder("geometry_frame");
  store_geometry_dialog(NULL);
  gtk_widget_set_visible(frame, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_geometry_angle_togglebutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

   if (gtk_toggle_button_get_active(togglebutton))
      do_angle_define();

}

extern "C" G_MODULE_EXPORT
void
on_geometry_torsion_togglebutton_toggled(GtkToggleButton *togglebutton,
                                         gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    do_torsion_define();
}

extern "C" G_MODULE_EXPORT
void
on_new_ligands_info_dialog_ok_button_clicked (GtkButton       *button,
                                              gpointer         user_data) {
   // GtkWidget *w = widget_from_builder("new_ligands_info_dialog");
   // gtk_widget_set_visible(w, FALSE);
   GtkWidget *w = widget_from_builder("new_ligands_info_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_no_new_ligands_info_dialog_ok_button_clicked (GtkButton       *button,
                                                 gpointer         user_data) {
   // GtkWidget *w = widget_from_builder("no_new_ligands_info_dialog");
   // gtk_widget_set_visible(w, FALSE);
   GtkWidget *w = widget_from_builder("no_new_ligands_info_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_zoom_dialog_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data) {

  GtkWidget *w = widget_from_builder("zoom_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data) {
   /* needs to set widget */
  unset_moving_atom_move_chis();
  set_show_chi_angle_bond(0);
}


extern "C" G_MODULE_EXPORT
void
on_check_waters_low_occ_dist_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_zero_occ_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

}



extern "C" G_MODULE_EXPORT
void
on_recover_coordinates_ok_button_clicked(GtkButton       *button,
                                         gpointer         user_data) {

   GtkWidget *widget = widget_from_builder("recover_coordinates_dialog");

  execute_recover_session(widget); /* widget needed for lookup of user data */
  gtk_widget_set_visible(widget, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_recover_coordinates_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {
   GtkWidget *widget = widget_from_builder("recover_coordinates_dialog");
   gtk_widget_set_visible(widget, FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_centre_atom_label_ok_button_clicked (GtkButton       *button,
                                                            gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("centre_atom_label_dialog");
   GtkWidget *on  = widget_from_builder("centre_atom_label_radiobutton_on");

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(on))) {
      set_label_on_recentre_flag(1);
   } else {
      set_label_on_recentre_flag(0);
   }
   gtk_widget_set_visible(dialog, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_centre_atom_label_radiobutton_on_toggled (GtkToggleButton       *button,
                                                                 gpointer         user_data)
{
   GtkWidget *on  = widget_from_builder("centre_atom_label_radiobutton_on");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(on))) {
      set_label_on_recentre_flag(1);
   } else {
      set_label_on_recentre_flag(0);
   }
}



extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_help_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
   //   GtkWidget *w = create_chi_angle_help_dialog();
   GtkWidget *w = widget_from_builder("chi_angle_help_dialog");
   gtk_widget_set_visible(w, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_help_chi_angles_dismiss_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("chi_angle_help_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_rotate_translate_obj_dialog_destroy (GtkWidget       *object,
                                                            gpointer         user_data)
{
   /* need to save the position coordinates of dialog */
   rot_trans_reset_previous();
}


extern "C" G_MODULE_EXPORT
void
on_no_symmetry_warning_ok_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("no_symmetry_warning_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_nothing_to_recover_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("nothing_to_recover_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_superpose_dialog_superpose_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_dialog");
  execute_superpose(w);
  gtk_widget_set_visible(w, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_superpose_dialog_cancel_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_superpose_nonsense_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_nonsense_dialog");
  gtk_widget_set_visible(w, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_superpose_nonsense_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_dialog");
  gtk_widget_set_visible(w, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_add_terminal_residue_finds_none_ok_button_clicked (GtkButton       *button,
						      gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("add_terminal_residue_finds_none_dialog");
  gtk_widget_set_visible(w, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_single_map_properties_contour_level_entry_activate(GtkWidget *entry, gpointer user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   std::string s = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      float f = coot::util::string_to_float(s);
      set_contour_by_sigma_step_by_mol(imol, f, 1);
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR in on_single_map_properties_contour_level_entry_activate() failed to comprehend "
                << s << std::endl;
   }

}


extern "C" G_MODULE_EXPORT
void
on_single_map_properties_absolute_radiobutton_toggled (GtkCheckButton *checkbutton,
                                                       gpointer         user_data) {

   // std::cout << "Absolute button toggled - doing nothing" << std::endl;
}

extern "C" G_MODULE_EXPORT
void
handle_map_properties_fresnel_change(int imol, GtkWidget *checkbutton) {

   std::cout << "Here 0 in handle_map_properties_fresnel_change() " << checkbutton << std::endl;
   if (! graphics_info_t::is_valid_map_molecule(imol)) return;
   std::cout << "Here A in handle_map_properties_fresnel_change() " << checkbutton << std::endl;

   molecule_class_info_t &m = graphics_info_t::molecules[imol];
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) {
      std::cout << "Here B in handle_map_properties_fresnel_change() " << checkbutton << std::endl;
      GtkWidget *bias_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton),  "bias_entry"));
      GtkWidget *scale_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton), "scale_entry"));
      GtkWidget *power_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton), "power_entry"));
      std::string  bias_entry_text  = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(bias_entry)));
      std::string scale_entry_text  = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(scale_entry)));
      std::string power_entry_text  = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(power_entry)));
      try {
         float bias  = coot::util::string_to_float(bias_entry_text);
         float scale = coot::util::string_to_float(scale_entry_text);
         float power = coot::util::string_to_float(power_entry_text);
         std::cout << "Here C in handle_map_properties_fresnel_change() " << bias << " " << scale << " " << power << std::endl;
         m.fresnel_settings.update_settings(true, bias, scale, power);
      }
      catch (const std::runtime_error &rte) {
         std::cout << "WARNING:: Failed to parse " << rte.what() << " " << bias_entry_text<< std::endl;
      }
   } else {
      m.fresnel_settings.state = false; // something touched me deep inside
   }
   graphics_draw();
}


extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_fresnel_state_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                           gpointer         user_data) {

   std::cout << "Toggled! " << std::endl;
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "imol"));
   handle_map_properties_fresnel_change(imol, GTK_WIDGET(checkbutton));
}

extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_fresnel_bias_entry_activate(GtkEntry* self, gpointer user_data) {

   std::cout << "bias entry key press activate" << std::endl;
   GtkWidget *checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "fresnel_checkbutton"));
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(self), "imol"));
   handle_map_properties_fresnel_change(imol, checkbutton);
}

extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_fresnel_scale_entry_activate(GtkEntry* self, gpointer user_data) {

   std::cout << "scale entry key press activate" << std::endl;
   GtkWidget *checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "fresnel_checkbutton"));
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(self), "imol"));
   handle_map_properties_fresnel_change(imol, checkbutton);
}

extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_fresnel_power_entry_activate(GtkEntry* self, gpointer user_data) {

   std::cout << "popwer entry key press activate" << std::endl;
   GtkWidget *checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "fresnel_checkbutton"));
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(self), "imol"));
   handle_map_properties_fresnel_change(imol, checkbutton);
}


extern "C" G_MODULE_EXPORT
void
on_no_bad_chiral_volumes_dialog_ok_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("no_bad_chiral_volumes_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_check_chiral_volumes_ok_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("check_chiral_volumes_dialog");
   check_chiral_volumes_from_widget(w);
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_check_chiral_volumes_cancel_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("check_chiral_volumes_dialog");
   gtk_widget_set_visible(dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_bad_chiral_volumes_dialog_response(GtkDialog       *dialog,
                                                          gint             response_id,
                                                          gpointer         user_data) {

   if (response_id == GTK_RESPONSE_CLOSE)
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_rigid_body_refinement_failed_dialog_ok_button_clicked (GtkButton *button,
                                                                              gpointer user_data)
{
  GtkWidget *w = widget_from_builder(
			       "rigid_body_refinement_failed_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_baton_mode_calculate_skeleton_ok_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  GtkWidget *w = widget_from_builder(
			       "baton_mode_make_skeleton_dialog");
  baton_mode_calculate_skeleton(w); /* get the imol from here */
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_baton_mode_calculate_skeleton_cancel_button_clicked (GtkButton       *button,
							gpointer         user_data)
{
  GtkWidget *w = widget_from_builder(
			       "baton_mode_make_skeleton_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_column_label_expert_mode_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *frame = widget_from_builder(
				   "column_labels_resolution_limit_frame");

/*   GtkWidget *f_optionmenu = widget_from_builder("optionmenu1"); */

  GtkWidget *combobox = widget_from_builder("column_selector_amplitudes_combobox");

  /* we also need to redo the F column label chooser to include anomalous option  */


   if (gtk_widget_get_visible(frame)) {
    gtk_widget_set_visible(frame, FALSE);
   } else {
     gtk_widget_set_visible(frame, TRUE);
/*      fill_f_optionmenu_with_expert_options(f_optionmenu); */
     fill_combobox_with_expert_options(combobox);
   }

}


extern "C" G_MODULE_EXPORT
void
on_column_labels_use_resolution_limits_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *frame = widget_from_builder("resolution_limits_hbox");
  if (gtk_toggle_button_get_active(togglebutton))
     gtk_widget_set_sensitive(frame, TRUE);
  else
     gtk_widget_set_sensitive(frame, FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_merge_molecules_ok_button_clicked(GtkButton       *button,
                                     gpointer         user_data) {

   GtkWidget *w = widget_from_builder("merge_molecules_dialog");
   do_merge_molecules(w);
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
gboolean
on_merge_molecules_dialog_close_request(GtkAboutDialog *dialog, gpointer user_data) {
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   return TRUE; // Prevent the default close behavior (destruction)
}


extern "C" G_MODULE_EXPORT
void
on_merge_molecules_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("merge_molecules_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("mutate_sequence_dialog");
   do_mutate_sequence(w);
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("mutate_sequence_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_draw_hydrogens_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  // Applying the bond parameters applies the bond width too, which is
  // a bit of extra overhead.
   GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   apply_bond_parameters(w);
}


extern "C" G_MODULE_EXPORT
void
on_draw_hydrogens_no_radiobutton_toggled(GtkToggleButton *togglebutton,
                                         gpointer         user_data) {

  /*

    20150803-PE:
    We don't need to do this - because this is a radiobutton with 2
    options, if we turn on the no, then the yes is turned off so we
    run the above on_draw_hydrogens_no_radiobutton_toggled().

    If this code is activated, then we call apply_bond_parameters()
    and thus set_draw_hydrogens() twice.

  GtkWidget *w = widget_from_builder("bond_parameters_dialog");
  apply_bond_parameters(w);
  */
}

extern "C" G_MODULE_EXPORT
void
on_renumber_residues_molecule_combobox_changed(GtkComboBox     *combobox,
                                                                   gpointer         user_data) {

   // This chaged signal is attached to the combobox items in new_fill_combobox_with_coordinates_options().
   // We could do it here, I suppose, attached the the combobox, but renumber_residue_range is from old code
   // and it seems to work.
   
}

extern "C" G_MODULE_EXPORT
void
on_renumber_residue_range_radiobutton_3_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  GtkWidget *entry_2 = widget_from_builder("renumber_residue_range_resno_2_entry");
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    gtk_widget_set_sensitive(GTK_WIDGET(entry_2), TRUE);
  } else {
    gtk_widget_set_sensitive(GTK_WIDGET(entry_2), FALSE);
  }
}


extern "C" G_MODULE_EXPORT
void
on_renumber_residue_range_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("renumber_residue_range_dialog");
  bool status = renumber_residues_from_widget(w);
  if (status)
     gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_renumber_residue_range_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("renumber_residue_range_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_add_OXT_c_terminus_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_add_OXT_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_add_OXT_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("add_OXT_dialog");
  apply_add_OXT_from_widget(w);
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_add_OXT_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *w = widget_from_builder("add_OXT_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_add_OXT_button_clicked(GtkButton       *button,
                                              gpointer         user_data) {

   GtkWidget *w = wrapped_create_add_OXT_dialog(); // uses builder
   set_transient_for_main_window(w);
   gtk_widget_set_visible(w, TRUE);
}



extern "C" G_MODULE_EXPORT
void
on_bond_parameters_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data) {

   // GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   // apply_bond_parameters(w);
   // gtk_widget_set_visible(w, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_bond_parameters_apply_button_clicked(GtkButton       *button,
                                        gpointer         user_data) {

   // GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   // apply_bond_parameters(w);
}




extern "C" G_MODULE_EXPORT
void
on_bond_parameters_close_button_clicked (GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   gtk_widget_set_visible(w, FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_ligand_no_blobs_OK_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   //   GtkWidget *dialog = widget_from_builder("ligand_no_blobs_dialog");
   // gtk_widget_set_visible(dialog, FALSE);

   GtkWidget *dialog = widget_from_builder("ligand_no_blobs_dialog");
   gtk_widget_set_visible(dialog, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_new_delete_molecules_ok_button_clicked(GtkButton       *button,
                                          gpointer         user_data) {

   GtkWidget *w = widget_from_builder("new_close_molecules_dialog");
   close_molecules_gtk4(w);
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_new_delete_molecules_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("new_close_molecules_dialog");
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_find_blobs_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("unmodelled_blobs_dialog");
  execute_find_blobs_from_widget(w);
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_find_blobs_cancel_button_clicked    (GtkButton       *button,
                                                            gpointer         user_data) {
  GtkWidget *w = widget_from_builder("unmodelled_blobs_dialog");
  gtk_widget_set_visible(w, FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_chiral_restraints_problem_ok_button_clicked (GtkButton       *button,
						gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("chiral_restraints_problem_dialog");
  gtk_widget_set_visible(w, FALSE);
}

/* We'll keep this for now because it is used by
   create_check_waters_diff_map_dialog() and that is still in the
   interface definition (glade file).
 */
extern "C" G_MODULE_EXPORT
void
on_check_waters_diff_map_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{

/*   GtkWidget *w = widget_from_builder("check_waters_diff_map_dialog"); */
/*   check_waters_by_difference_map_by_widget(w); */
/*   gtk_widget_set_visible(w, FALSE); */


}


extern "C" G_MODULE_EXPORT
void
on_check_waters_diff_map_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("check_waters_diff_map_dialog");
  gtk_widget_set_visible(w, FALSE);

}




extern "C" G_MODULE_EXPORT
void
on_interesting_waters_by_difference_map_check_ok_button_clicked (GtkButton       *button,
								 gpointer         user_data)
{

  GtkWidget *w = widget_from_builder(
			       "interesting_waters_by_difference_map_check_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_nothing_bad_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nothing_bad_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_skeletonize_map_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_skeletonize_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_antialias_dialog_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

  if (gtk_toggle_button_get_active(togglebutton))
    set_do_anti_aliasing(1);
  else
    set_do_anti_aliasing(0);

}


extern "C" G_MODULE_EXPORT
void
on_antialias_dialog_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_do_anti_aliasing(0);
  else
    set_do_anti_aliasing(1);

}


extern "C" G_MODULE_EXPORT
void
on_antialiasing_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("antialiasing_dialog");
  gtk_widget_set_visible(w, FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coords_fileselection_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("save_symmetry_coords_fileselection");
  save_symmetry_coords_from_filechooser(w);
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coords_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("save_symmetry_coords_fileselection");
  gtk_widget_set_visible(w, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_pepflips_by_difference_map_dialog_close (GtkDialog *dialog,
                                            gpointer   user_data) {
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

// use a header for this - which one? c-interface-gui.hh?
void pepflips_by_difference_map_results_dialog(int imol_coords, int imol_map, float n_sigma);

extern "C" G_MODULE_EXPORT
void
on_pepflips_by_difference_map_dialog_response(GtkDialog       *dialog,
                                              gint             response_id,
                                              gpointer         user_data) {
   if (response_id == GTK_RESPONSE_OK) {
      std::cout << "............... repsonse OK " << std::endl;
      GtkWidget *model_combobox = GTK_WIDGET(g_object_get_data(G_OBJECT(dialog), "model_combobox"));
      GtkWidget   *map_combobox = GTK_WIDGET(g_object_get_data(G_OBJECT(dialog),   "map_combobox"));
      GtkWidget *entry = widget_from_builder("pepflips_by_difference_map_dialog_entry");
      std::string s = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      try {
         int imol_coords = my_combobox_get_imol(GTK_COMBO_BOX(model_combobox));
         int imol_map    = my_combobox_get_imol(GTK_COMBO_BOX(  map_combobox));
         float n_sigma = coot::util::string_to_float(s);
         pepflips_by_difference_map_results_dialog(imol_coords, imol_map, n_sigma);
      }
      catch (const std::runtime_error &rte) {
         // log this
         std::cout << "Failed to convert " << s << " to a number" << std::endl;
      }
   } else {
      std::cout << "response was not OK" << std::endl;
   }
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("stereo_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_mono_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    mono_mode();
}


extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_hardware_stereo_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *nothing_bad_dialog;
  GtkWidget *label_widget;
  GtkWidget *mono_togglebutton;

  if (gtk_toggle_button_get_active(togglebutton)) {
    hardware_stereo_mode();

    if (stereo_mode_state() != 1) { /* coot::HARDWARE_STEREO_MODE */
      mono_togglebutton = widget_from_builder("stereo_dialog_mono_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mono_togglebutton), TRUE);
      /* do this in C rather than mess about calling c++ function: */
      // nothing_bad_dialog = create_nothing_bad_dialog();
      nothing_bad_dialog = widget_from_builder("nothing_bad_dialog");
      label_widget = widget_from_builder("nothing_bad_label");
      gtk_label_set_text(GTK_LABEL(label_widget), "This computer appears not to be able\nto do hardware stereo");
      gtk_widget_set_visible(nothing_bad_dialog, TRUE);
    }
  }
}

extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_side_by_side_stereo_crosseyed_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    side_by_side_stereo_mode(0); /* passed used_wall_eye flag */
}


extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_side_by_side_stereo_walleyed_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    side_by_side_stereo_mode(1); /* passed used_wall_eye flag */
  }
}

extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_zalman_stereo_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    zalman_stereo_mode();
  }

}


extern "C" G_MODULE_EXPORT
gboolean
on_shader_settings_dialog_close_request(GtkWidget       *widget,
                                  gpointer         user_data) {

   // 2026-02-21-PE this now does happen
   gtk_widget_set_visible(widget, FALSE);
   return TRUE;
}

extern "C" G_MODULE_EXPORT
void
on_shader_settings_dialog_close(GtkDialog *dialog,
                                                    gpointer   user_data) {

   // The ::close signal is a keybinding signal which gets emitted when the user uses a keybinding to close the dialog.
   // The default binding for this signal is the Escape key.

   std::cout << "-------------- on_shader_settings_dialog close  " << std::endl;
}


/* end preferences */

extern "C" G_MODULE_EXPORT
void
on_generate_diff_map_peaks_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("generate_diff_map_peaks_dialog");
   difference_map_peaks_from_dialog(); // make the results (and show them)
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_generate_diff_map_peaks_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("generate_diff_map_peaks_dialog");
   gtk_widget_set_visible(w, FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_superpose_reference_chain_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                 gpointer user_data) {

   GtkWidget *combobox = widget_from_builder("superpose_dialog_reference_chain_combobox");
   if (gtk_check_button_get_active(checkbutton)) {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);

      GtkWidget *combobox_1 = widget_from_builder("superpose_dialog_reference_mol_combobox");
      GtkWidget *combobox_2 = widget_from_builder("superpose_dialog_moving_mol_combobox");

      int imol1 = my_combobox_get_imol(GTK_COMBO_BOX(combobox_1));
      int imol2 = my_combobox_get_imol(GTK_COMBO_BOX(combobox_2));
      int imol_active = imol1;

      fill_superpose_combobox_with_chain_options(imol_active, 1);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_superpose_moving_chain_checkbutton_toggled(GtkCheckButton *checkbutton,
                                              gpointer        user_data) {

   GtkWidget *combobox = widget_from_builder("superpose_dialog_moving_chain_combobox");
   if (gtk_check_button_get_active(checkbutton)) {
      GtkWidget *combobox_2 = widget_from_builder("superpose_dialog_moving_mol_combobox");
      int imol2 = my_combobox_get_imol(GTK_COMBO_BOX(combobox_2));
      int imol_active = imol2;

      fill_superpose_combobox_with_chain_options(imol_active, 0);
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
   }
}

// is this used now?
extern "C" G_MODULE_EXPORT
void
on_draw_ncs_ghosts_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {
/* Function no longer used. Handled in another dialog

   Kept in glade (not visible) for historical reasons

   GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   if (gtk_toggle_button_get_active(togglebutton)) {
      printf("yes radiobutton toggled on.\n");
      make_ncs_ghosts_maybe(w);
   }
*/
}


// is this used now?
extern "C" G_MODULE_EXPORT
void
on_draw_ncs_ghosts_no_radiobutton_toggled(GtkToggleButton *togglebutton,
                                          gpointer         user_data) {

}

extern "C" G_MODULE_EXPORT
void
aniso_probability_hscale_value_changed(GtkScale* range,
                                       gpointer user_data) {

   GtkWidget *bond_parameters_molecule_comboboxtext  = widget_from_builder("bond_parameters_molecule_comboboxtext");
   GtkWidget *draw_anisotropic_atoms_yes_radiobutton = widget_from_builder("draw_anisotropic_atoms_yes_radiobutton");
   if (bond_parameters_molecule_comboboxtext) {
      if (draw_anisotropic_atoms_yes_radiobutton) {
         GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
         float fvalue = gtk_adjustment_get_value(adjustment);
         graphics_info_t g;
         g.show_aniso_atoms_probability = fvalue;
         int imol = g.combobox_get_imol(GTK_COMBO_BOX(bond_parameters_molecule_comboboxtext));
         if (gtk_check_button_get_active(GTK_CHECK_BUTTON(draw_anisotropic_atoms_yes_radiobutton))) {
            // std::cout << "\ncalling set_show_atoms_as_aniso() with prob " << g.show_aniso_atoms_probability << std::endl;
            graphics_info_t::molecules[imol].make_bonds_type_checked("aniso_probability_hscale_value_changed");
            g.graphics_draw();
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_draw_anisotropic_atoms_yes_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                  gpointer         user_data) {

   GtkWidget *bond_parameters_molecule_comboboxtext =
      widget_from_builder("bond_parameters_molecule_comboboxtext");
   if (bond_parameters_molecule_comboboxtext) {
      graphics_info_t g;
      int imol = g.combobox_get_imol(GTK_COMBO_BOX(bond_parameters_molecule_comboboxtext));
      if (gtk_check_button_get_active(checkbutton)) {
         set_show_aniso_atoms(imol, 1);
      } else {
         set_show_aniso_atoms(imol, 0);
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_draw_anisotropic_atoms_no_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                 gpointer         user_data) {

   // handled above

}

extern "C" G_MODULE_EXPORT
void
show_anisotropic_atoms_as_ortep_switch_state_set(GtkSwitch *switch_widget,
                                                 gboolean   state,
                                                 gpointer   user_data) {

   GtkWidget *bond_parameters_molecule_comboboxtext =
      widget_from_builder("bond_parameters_molecule_comboboxtext");

   if (bond_parameters_molecule_comboboxtext) {
      graphics_info_t g;
      int imol = g.combobox_get_imol(GTK_COMBO_BOX(bond_parameters_molecule_comboboxtext));
      set_show_aniso_atoms_as_ortep(imol, state);
      if (state) {
         // for the aniso button on:
         GtkWidget *checkbutton = widget_from_builder("draw_anisotropic_atoms_yes_radiobutton");
         if (checkbutton) {
            gtk_check_button_set_active(GTK_CHECK_BUTTON(checkbutton), TRUE);
         }
      }
   }
}



extern "C" G_MODULE_EXPORT
void
on_ncs_maps_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ncs_maps_dialog");
   make_dynamically_transformed_ncs_maps_by_widget(w);
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_ncs_maps_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ncs_maps_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_pointer_distances_checkbutton_toggled(GtkCheckButton *checkbutton,
                                         gpointer        user_data) {

   if (gtk_check_button_get_active(checkbutton)) {
      printf("pointer distances toggle button toggled on\n");
   } else {
      printf("pointer distances toggle button toggled off\n");
   }
   toggle_pointer_distances_show_distances(checkbutton);
}


extern "C" G_MODULE_EXPORT
void
on_pointer_distances_min_dist_entry_activate(GtkEntry        *entry,
                                             gpointer         user_data) {

   const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      float f = coot::util::string_to_float(std::string(text));
      graphics_info_t g;
      g.pointer_min_dist = f;
      g.make_pointer_distance_objects();
      g.graphics_draw();
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }

}



extern "C" G_MODULE_EXPORT
void
on_pointer_distances_max_dist_entry_activate(GtkEntry        *entry,
                                             gpointer         user_data) {

   const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   try {
      float f = coot::util::string_to_float(std::string(text));
      graphics_info_t g;
      g.pointer_max_dist = f;
      g.make_pointer_distance_objects();
      g.graphics_draw();
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_pointer_distances_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("pointer_distances_dialog");
   execute_pointer_distances_settings(dialog); // upadte and graphics_draw()
   gtk_widget_set_visible(dialog, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_align_and_mutate_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *dialog = widget_from_builder("align_and_mutate_dialog");
   GtkWidget *dialog = widget_from_builder("align_and_mutate_dialog");
   int handled_state = do_align_mutate_sequence(dialog);
   if (handled_state == 1)
      gtk_widget_set_visible(dialog, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_align_and_mutate_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("align_and_mutate_dialog");
   gtk_widget_set_visible(dialog, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_differences_ok_button_clicked(GtkButton       *button,
                                                   gpointer         user_data) {

   GtkWidget *w = widget_from_builder("ramachandran_plot_differences_dialog");
   int istat = do_ramachandran_plot_differences_by_widget(w);
   if (istat) 			/* the plot was drawn (i.e. no chain selection funnies) */
      gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_differences_cancel_button_clicked (GtkButton       *button,
                                                        gpointer         user_data) {

   GtkWidget *w = widget_from_builder("ramachandran_plot_differences_dialog");
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_differences_first_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *combobox = widget_from_builder(
				       "ramachandran_plot_differences_first_chain_combobox");
   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
      fill_ramachandran_plot_differences_combobox_with_chain_options(combobox, 1);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
   }

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_differences_second_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *combobox = widget_from_builder(
				       "ramachandran_plot_differences_second_chain_combobox");

   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
      fill_ramachandran_plot_differences_combobox_with_chain_options(combobox, 0);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
   }

}


extern "C" G_MODULE_EXPORT
void
on_checked_waters_baddies_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("checked_waters_baddies_dialog");
   gtk_widget_set_visible(w, FALSE);

}



extern "C" G_MODULE_EXPORT
void
on_delete_item_keep_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
}




extern "C" G_MODULE_EXPORT
void
on_main_window_statusbar_text_pushed   (GtkStatusbar    *statusbar,
                                        guint            context_id,
                                        gchar           *text,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_fit_loop_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("mutate_sequence_dialog");
   fit_loop_using_dialog();
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_A_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
   do_base_mutation("A");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_C_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
   do_base_mutation("C");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_G_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
   do_base_mutation("G");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_T_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
   do_base_mutation("T");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_U_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
   do_base_mutation("U");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_nucleic_acid_base_chooser_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
   clear_pending_picks();
}



extern "C" G_MODULE_EXPORT
void
on_change_chains_rechain_button_clicked (GtkButton       *button,
                                        gpointer         user_data) {

   // GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   change_chain_id_by_widget(w);
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_change_chain_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   // gtk_widget_set_visible(w, FALSE);

   GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_delete_item_residue_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_delete_residue_zone_mode();

}


extern "C" G_MODULE_EXPORT
void
on_on_line_documentation_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("doc_urls_dialog");
   gtk_widget_set_visible(w, FALSE);

}




extern "C" G_MODULE_EXPORT
void
on_change_chain_residue_range_no_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                     gpointer        user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_change_chain_residue_range_yes_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                      gpointer        user_data) {

   GtkWidget *hbox = widget_from_builder("change_chain_id_residue_range_hbox");

   if (gtk_check_button_get_active(checkbutton)) {
      gtk_widget_set_sensitive(hbox, TRUE);
      std::cout << "residue range hbox sensitive TRUE" << std::endl;
   } else {
      gtk_widget_set_sensitive(hbox, FALSE);
      std::cout << "residue range hbox sensitive FALSE" << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_do_autofit_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                  gpointer        user_data) {

   int imol_map = -1;

   if (gtk_check_button_get_active(checkbutton)) {
      imol_map = imol_refinement_map();
      if (imol_map == -1) {
	 gtk_check_button_set_active(checkbutton, FALSE);
	 show_select_map_frame();
	 info_dialog("A map has not yet been assigned for Refinement/Fitting");
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_use_ramachandran_restraints_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                                   gpointer        user_data) {
   /* not doing anything because the button state read at execution time  */
}

extern "C" G_MODULE_EXPORT
void
on_check_waters_b_factor_entry_active_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                          gpointer         user_data) {

   // GtkWidget *hbox = widget_from_builder("check_waters_b_factor_hbox");
   GtkWidget *entry = widget_from_builder("check_waters_b_factor_entry");
   if (gtk_check_button_get_active(checkbutton))
      gtk_widget_set_sensitive(entry, TRUE);
   else
      gtk_widget_set_sensitive(entry, FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_check_waters_min_dist_entry_active_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                          gpointer        user_data) {

   GtkWidget *hbox = widget_from_builder("check_waters_min_dist_hbox");
   GtkWidget *entry = widget_from_builder("check_waters_min_dist_entry");
   if (gtk_check_button_get_active(checkbutton))
      gtk_widget_set_sensitive(entry, TRUE);
   else
      gtk_widget_set_sensitive(entry, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_max_dist_entry_active_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                          gpointer         user_data) {

   // GtkWidget *hbox = widget_from_builder("check_waters_max_dist_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_max_dist_hbox");
   GtkWidget *entry = widget_from_builder("check_waters_max_dist_entry");
   if (gtk_check_button_get_active(checkbutton))
      gtk_widget_set_sensitive(entry, TRUE);
   else
      gtk_widget_set_sensitive(entry, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_map_sigma_entry_active_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                           gpointer         user_data)
{
   // GtkWidget *hbox = widget_from_builder("check_waters_sigma_level_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_sigma_level_hbox");
   if (gtk_check_button_get_active(checkbutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_check_waters_by_difference_map_active_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                             gpointer        user_data) {

   // GtkWidget *hbox = widget_from_builder("check_waters_by_difference_map_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_by_difference_map_hbox");
   if (gtk_check_button_get_active(checkbutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_residue_info_occ_apply_all_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                  gpointer        user_data) {

   GtkWidget *entry = widget_from_builder("residue_info_master_atom_occ_entry");
   GtkWidget *alt_conf_checkbutton = widget_from_builder("residue_info_occ_apply_to_altconf_checkbutton");

   if (gtk_check_button_get_active(checkbutton)) {
      gtk_widget_set_sensitive(entry, TRUE);
   } else {
      if (! gtk_check_button_get_active(GTK_CHECK_BUTTON(alt_conf_checkbutton)))
         gtk_widget_set_sensitive(entry, FALSE);
   }
}



extern "C" G_MODULE_EXPORT
void
on_residue_info_b_factor_apply_all_checkbutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("residue_info_master_atom_b_factor_entry");
   if (gtk_check_button_get_active(checkbutton))
      gtk_widget_set_sensitive(entry, TRUE);
   else
      gtk_widget_set_sensitive(entry, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_other_modelling_tools_close_button_clicked(GtkButton       *button,
                                              gpointer         user_data) {

   GtkWidget *w = widget_from_builder("other_model_tools_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_cis_trans_conversion_toggle_button_toggled(GtkToggleButton *togglebutton,
                                              gpointer         user_data) {

  if (gtk_toggle_button_get_active(togglebutton))
      do_cis_trans_conversion_setup(1);
   else
      do_cis_trans_conversion_setup(0);

}


extern "C" G_MODULE_EXPORT
void
on_other_model_tools_dialog_destroy    (GtkWidget       *object,
                                        gpointer         user_data) {

   std::cout << "---------------- this should not happen on_other_model_tools_dialog_destroy " << std::endl;
   do_cis_trans_conversion_setup(0);
   unset_other_modelling_tools_dialog();
}

extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_do_180_degree_sidechain_flip_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    setup_180_degree_flip(1);
  else
    setup_180_degree_flip(0);

}


extern "C" G_MODULE_EXPORT
void
on_screendump_image_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *filechooser = widget_from_builder(
					    "screendump_filechooser"); /* now consistent with above */

   if (filechooser)
     gtk_widget_set_visible(filechooser, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_reverse_fragment_direction_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
      setup_reverse_direction(1);
   else
      setup_reverse_direction(0);
}


extern "C" G_MODULE_EXPORT
void
on_place_helix_here_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  place_helix_here();
}


extern "C" G_MODULE_EXPORT
void
on_diff_map_peaks_dialog_destroy(GtkWidget       *object,
                                        gpointer         user_data) {
  set_difference_map_peaks_widget(0); /* a null pointer */
}


extern "C" G_MODULE_EXPORT
void
on_symmetry_controller_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("symmetry_controller_dialog");
  gtk_widget_set_visible(w, FALSE);

}





extern "C" G_MODULE_EXPORT
void
on_show_symmetry_molecule_control_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = symmetry_molecule_controller_dialog();
  gtk_widget_set_visible(w, TRUE);
}



extern "C" G_MODULE_EXPORT
void
on_ncs_control_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ncs_control_dialog");
   gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_lsq_plane_add_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
      setup_lsq_plane_define(1);
      setup_lsq_deviation(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_lsq_plane_deviant_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
      setup_lsq_deviation(1);
      setup_lsq_plane_define(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_lsq_plane_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("lsq_plane_dialog");
  gtk_widget_set_visible(w, FALSE);
  normal_cursor();
}


extern "C" G_MODULE_EXPORT
void
on_lsq_plane_dialog_destroy            (GtkWidget       *object,
                                        gpointer         user_data)
{
  unset_lsq_plane_dialog();	/* which clears the plane points too */
  normal_cursor();
}


extern "C" G_MODULE_EXPORT
void
on_lsq_plane_delete_last_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  remove_last_lsq_plane_atom();
}


extern "C" G_MODULE_EXPORT
void
on_coord_colour_control_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  GtkWidget *w = widget_from_builder("coords_colour_control_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_coot_doc_url_monolithic_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  browser_url("https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/web/docs/coot.html");

}


extern "C" G_MODULE_EXPORT
void
on_coot_doc_url_sectioned_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  browser_url("https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/web/docs/coot.html");
}


extern "C" G_MODULE_EXPORT
void
on_coot_online_doc_search_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkEntry *entry = GTK_ENTRY(widget_from_builder("coot_online_doc_search_entry"));
  const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  handle_online_coot_search_request(text);
}


extern "C" G_MODULE_EXPORT
void
on_bond_parameters_rotate_colour_map_c_only_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  short int i=0;
  if (gtk_toggle_button_get_active(togglebutton))
    i=1;
  set_colour_map_rotation_on_read_pdb_c_only_flag(i);
}



extern "C" G_MODULE_EXPORT
void
on_generic_objects_dialog_response(GtkDialog       *dialog,
                                   gint             response_id,
                                   gpointer         user_data) {

    if (response_id == GTK_RESPONSE_CLOSE) {
       gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
    }
}

extern "C" G_MODULE_EXPORT
void
on_generic_objects_dialog_close_button_clicked(GtkButton       *button,
                                               gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("generic_objects_dialog");
   gtk_widget_set_visible(dialog, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_map_radius_x_ray_apply_button_clicked(GtkButton       *button,
                                                               gpointer         user_data) {

   GtkEntry *entry_xray = GTK_ENTRY(widget_from_builder("map_parameters_x_ray_radius_entry"));
   GtkEntry *entry_em   = GTK_ENTRY(widget_from_builder("map_parameters_em_radius_entry"));
   const char *text = gtk_editable_get_text(GTK_EDITABLE(entry_xray));
   set_density_size_from_widget(text);
}

extern "C" G_MODULE_EXPORT
void
on_map_radius_em_apply_button_clicked(GtkButton       *button,
                                                          gpointer         user_data) {

   GtkEntry *entry_xray = GTK_ENTRY(widget_from_builder("map_parameters_x_ray_radius_entry"));
   GtkEntry *entry_em   = GTK_ENTRY(widget_from_builder("map_parameters_em_radius_entry"));
   const char *text = gtk_editable_get_text(GTK_EDITABLE(entry_em));
   set_density_size_em_from_widget(text);
}




extern "C" G_MODULE_EXPORT
void
on_refine_params_use_helix_peptide_torsions_radiobutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data)
{
  /* not visible */
  printf("helix togglebutton toggled - ignored\n");
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_beta_strand_peptide_torsions_radiobutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data)
{
  /* not visible */
  printf("beta strand togglebutton toggled - ignored\n");

}


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_ramachandran_goodness_torsions_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                                        gpointer         user_data) {

  int state = 0;
  if (gtk_check_button_get_active(checkbutton)) {
    state = 1;
  }
  set_refine_ramachandran_angles(state);
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_peptide_omegas_checkbutton_toggled
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data)
{
  if (gtk_check_button_get_active(checkbutton)) {
    add_omega_torsion_restriants();
  } else {
    remove_omega_torsion_restriants();
  }
}

extern "C" G_MODULE_EXPORT
void
on_refine_params_use_torsions_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                  gpointer         user_data) {

   do_torsions_toggle(GTK_WIDGET(checkbutton)); // 20231007-PE very old function
}

extern "C" G_MODULE_EXPORT
void
on_other_tools_RNA_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   // non-run function
  // GtkWidget *w = wrapped_nucleotide_builder_dialog();
  // gtk_widget_set_visible(w, TRUE);
}



extern "C" G_MODULE_EXPORT
void
on_ideal_rna_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("nucleotide_builder_dialog");
  ideal_nucleic_acid_by_widget(w);
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_ideal_rna_cancel_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("nucleotide_builder_dialog");
  gtk_widget_set_visible(w, FALSE);
}


// extern "C" G_MODULE_EXPORT
// void
// on_unit_cell_yes_radiobutton_toggled   (GtkCheckButton *checkbutton,
//                                         gpointer        user_data)
// {
//    if (gtk_check_button_get_active(checkbutton))
//       set_show_unit_cells_all(1);
//   else
//       set_show_unit_cells_all(0);
// }


// extern "C" G_MODULE_EXPORT
// void
// on_unit_cell_no_radiobutton_toggled(GtkCheckButton *checkbutton,
//                                     gpointer        user_data) {

//    if (gtk_check_button_get_active(checkbutton))
//       set_show_unit_cells_all(0);
//    else
//       set_show_unit_cells_all(1);
// }

extern "C" G_MODULE_EXPORT
void
show_unit_cell_switch_state_set(GtkSwitch *switch_widget,
                                gboolean   state,
                                gpointer   user_data) {

   if (state)
      set_show_unit_cells_all(1);
  else
     set_show_unit_cells_all(0);
}



extern "C" G_MODULE_EXPORT
void
on_move_molecule_here_ok_button_clicked(GtkButton       *button,
                                        gpointer         user_data) {

  move_molecule_here_by_widget();
  GtkWidget *w = widget_from_builder("move_molecule_here_frame");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_move_molecule_here_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("move_molecule_here_frame");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_monomer_library_search_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("monomer_search_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_monomer_library_search_button_clicked(GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("monomer_search_entry");
   GtkWidget *viewport = widget_from_builder("monomer_search_results_viewport");

   if (entry) {
      entry_char_type *text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      if (text) {
         handle_make_monomer_search(text, viewport);
      }
   }
}


extern "C" G_MODULE_EXPORT
void
on_least_squares_match_type_all_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_least_squares_match_type_main_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_least_squares_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("least_squares_dialog");
  apply_lsq_matches_by_widget(w);

}



extern "C" G_MODULE_EXPORT
void
on_least_squares_cancel_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("least_squares_dialog");
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_residue_type_chooser_stub_checkbutton_toggled (GtkToggleButton *togglebutton,
                                                  gpointer         user_data) {
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton)))
     set_residue_type_chooser_stub_state(1);
  else
     set_residue_type_chooser_stub_state(0);
}


extern "C" G_MODULE_EXPORT
void
on_set_undo_molecule_button_clicked    (GtkButton       *button,
                                        gpointer         user_data) {
  show_set_undo_molecule_chooser();
}

extern "C" G_MODULE_EXPORT
void
on_sec_str_rest_no_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_secondary_structure_restraints_type(0);
}


extern "C" G_MODULE_EXPORT
void
on_sec_str_rest_helix_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_secondary_structure_restraints_type(1);
}


extern "C" G_MODULE_EXPORT
void
on_sec_str_rest_strand_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_secondary_structure_restraints_type(2);
}


extern "C" G_MODULE_EXPORT
void
on_update_go_to_atom_from_current_position_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   update_go_to_atom_from_current_position();
}


extern "C" G_MODULE_EXPORT
void
on_symmetry_colorbutton_color_set      (GtkColorButton  *colorbutton,
                                        gpointer         user_data) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

   // need to update color usage
#else
  GdkColor colour;
  gdouble color[4]; // use first 3
  double r = 1.0 / 65535.0;
  gtk_color_button_get_color(colorbutton, &colour);
  color[0] = colour.red   * r;
  color[1] = colour.green * r;
  color[2] = colour.blue  * r;
  handle_symmetry_colour_change(1,color);
#endif

}

extern "C" G_MODULE_EXPORT
void
on_display_control_all_maps_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

   if (gtk_toggle_button_get_active(togglebutton))
      set_all_maps_displayed(1);
   else
      set_all_maps_displayed(0);

}


extern "C" G_MODULE_EXPORT
void
on_display_control_all_models_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

   if (gtk_toggle_button_get_active(togglebutton))
      set_all_models_displayed_and_active(1);
   else
      set_all_models_displayed_and_active(0);

}


extern "C" G_MODULE_EXPORT
void
on_single_map_properties_contour_level_apply_button_clicked(GtkButton       *apply_button,
                                                            gpointer         user_data) {

   //single_map_properties_apply_contour_level_to_map(w); /* check now
   //							  made here
   //							  for valid
   //							  map
   //							  molecule. */

   int imol    = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(apply_button), "imol"));
   GtkWidget *entry = GTK_WIDGET(g_object_get_data(G_OBJECT(apply_button), "contour_level_entry"));
   GtkWidget *checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(apply_button), "single_map_properties_absolute_radiobutton"));

   std::cout << "imol: " << imol << std::endl;
   std::cout << "entry " << entry << std::endl;
   std::cout << "checkbutton " << checkbutton << std::endl;

   if (is_valid_map_molecule(imol)) {
      std::string t = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      try {
         float f = coot::util::string_to_float(t);
         if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) {
            set_contour_level_absolute(imol, f);
         } else {
            set_contour_level_in_sigma(imol, f);
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "Failed to interpret " << t << std::endl;
      }
   }
   
}

extern "C" G_MODULE_EXPORT
void
on_display_map_style_as_lines_radiobutton_toggled(GtkCheckButton *checkbutton,
                                                  gpointer        user_data) {

   
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "imol"));
   std::cout << "on_display_map_style_as_lines_radiobutton_toggled() imol " << imol << std::endl;
   if (gtk_check_button_get_active(checkbutton)) {
      set_draw_map_standard_lines(imol, 1);
      set_draw_solid_density_surface(imol, 0);
   } else {
      set_draw_map_standard_lines(imol, 0);
      set_draw_solid_density_surface(imol, 1);
   }
}


extern "C" G_MODULE_EXPORT
void
on_checked_waters_baddies_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data) {
  store_checked_waters_baddies_dialog(NULL);
}




extern "C" G_MODULE_EXPORT
void
on_phs_coordinates_filechooserdialog1_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{

   GtkWidget *phs_filechooser = widget_from_builder("phs_coordinates_filechooserdialog1");
   if (response_id == GTK_RESPONSE_OK) {
     // filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(phs_fileselection));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(phs_filechooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);
     save_directory_from_filechooser(phs_filechooser);
     read_phs_and_coords_and_make_map(file_name);
   }
   gtk_widget_set_visible(phs_filechooser, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_phs_coordinates_filechooserdialog1_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{

  GtkWidget *phs_fileselection1 = widget_from_builder("phs_coordinates_filechooserdialog1");
  gtk_widget_set_visible(phs_fileselection1, FALSE);
}


#if (GTK_MAJOR_VERSION >= 4)
// GtkFileChooserConfirmation is no longer a thing.
#else
GtkFileChooserConfirmation
on_save_coords_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser,
					gpointer user_data) {

   if (file_chooser_overwrite_state() == 1) {
      return GTK_FILE_CHOOSER_CONFIRMATION_CONFIRM;
   } else {
      return GTK_FILE_CHOOSER_CONFIRMATION_ACCEPT_FILENAME;
   }
}
#endif


extern "C" G_MODULE_EXPORT
void
on_save_coords_filechooserdialog1_response
					(GtkDialog * dialog,
					gint response_id,
					gpointer user_data)
{
  if (response_id == GTK_RESPONSE_OK) {
    GtkWidget *fileselection = widget_from_builder("save_coords_filechooserdialog1");
    save_directory_for_saving_from_filechooser(fileselection);
    const char *stuff = static_cast<const char *>(g_object_get_data(G_OBJECT(fileselection), "stuff"));
    save_coordinates_using_widget(fileselection);
    // this is wrong // FIXME-LATER-PE
    // free(stuff);
    gtk_widget_set_visible(fileselection, FALSE);
  } else {
    GtkWidget *fileselection = widget_from_builder("save_coords_filechooserdialog1");
    gtk_widget_set_visible(fileselection, FALSE);
  }
}


extern "C" G_MODULE_EXPORT
void
on_save_coords_filechooserdialog1_destroy(GtkWidget * object,
                                          gpointer user_data) {

  GtkWidget *fileselection = widget_from_builder("save_coords_filechooserdialog1");
  gtk_widget_set_visible(fileselection, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_filechooserdialog1_destroy(GtkWidget * object,
                                             gpointer user_data) {

  GtkWidget *fileselection = widget_from_builder("cif_dictionary_filechooserdialog1");
  gtk_widget_set_visible(fileselection, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coords_filechooserdialog1_destroy
					(GtkWidget * object,
					gpointer user_data) {

  GtkWidget *coords_fileselection1 = widget_from_builder("save_symmetry_coords_filechooserdialog1");
  gtk_widget_set_visible(coords_fileselection1, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_screendump_filechooser_dialog_response (GtkDialog * dialog,
                                           gint response_id,
                                           gpointer user_data) {

   //GtkWidget *file_chooser = widget_from_builder("screendump_filechooser_dialog");
   if (response_id == GTK_RESPONSE_OK) {

      int image_type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "image_type"));
      // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
//      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
//                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_get_path(file);

      if (image_type == COOT_SCREENDUMP_SIMPLE) {
         screendump_tga(file_name);
      }
      if (image_type == COOT_SCREENDUMP_POVRAY) { // I doubt that this will ever work again.
         make_image_povray(file_name);
      }
      if (image_type == COOT_SCREENDUMP_RASTER3D) {
         make_image_raster3d(file_name);
      }
   }
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_torsion_general_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    setup_torsion_general(1);
  else
    setup_torsion_general(0);

}

extern "C" G_MODULE_EXPORT
void
on_accept_reject_reverse_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  toggle_torsion_general_reverse();
}

extern "C" G_MODULE_EXPORT
void
on_accept_reject_refinement_atom_pull_autoclear_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

   int state = get_auto_clear_atom_pull_restraint_state();
   if (state) {
     set_auto_clear_atom_pull_restraint(0);
   } else {
     set_auto_clear_atom_pull_restraint(1);
   }

}

extern "C" G_MODULE_EXPORT
void
on_accept_reject_atom_pull_clear_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {
   clear_all_atom_pull_restraints();
}



extern "C" G_MODULE_EXPORT
void
on_geometry_dynamic_distance_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    setup_dynamic_distances(1);
  else
    setup_dynamic_distances(0);
}



extern "C" G_MODULE_EXPORT
void
on_fixed_atom_dialog_close(GtkDialog       *dialog,
                                               gpointer         user_data) {
}

extern "C" G_MODULE_EXPORT
void
on_fixed_atom_dialog_response(GtkDialog       *dialog,
                                                  gint             response_id,
                                                  gpointer         user_data) {
   if (response_id == GTK_RESPONSE_CLOSE) {
      GtkWidget *dialog = widget_from_builder("fixed_atom_dialog");
      gtk_widget_set_visible(dialog, FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_fixed_atom_dialog_clear_all_fixed_atoms_button_clicked(GtkButton *togglebutton,
                                                                       gpointer         user_data) {

   clear_fixed_atoms_all(); // or maybe clear_all_fixed_atom(imol);
}


extern "C" G_MODULE_EXPORT
void
on_fixed_atom_dialog_fix_atom_togglebutton_toggled(GtkButton *togglebutton,
                                                                       gpointer         user_data) {

   setup_fixed_atom_pick(1, 0); // set a pending pick (not an unpick)
}

extern "C" G_MODULE_EXPORT
void
on_fixed_atom_dialog_unfix_atom_togglebutton_toggled(GtkButton *togglebutton,
                                                                       gpointer         user_data) {
   setup_fixed_atom_pick(1, 1);
}


extern "C" G_MODULE_EXPORT
void
on_add_reps_dialog_response(GtkDialog       *dialog,
                            gint             response_id,
                            gpointer         user_data)  {

   std::cout << "response " << response_id << std::endl;

   if (response_id == GTK_RESPONSE_OK) {

      // now read the dialog widgets to see what to represent...

      std::cout << "... read the dialog " << std::endl;

      add_additional_representation_by_dialog(dialog);

   }

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

}

// 
extern "C" G_MODULE_EXPORT
void
on_add_reps_dialog_close (GtkDialog *dialog,
                                              gpointer   user_data) {

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_display_additional_representations_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}

extern "C" G_MODULE_EXPORT
void
on_residue_editor_select_monomer_type_ok_button_clicked (G_GNUC_UNUSED GtkButton       *button,
                                                         G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("residue_editor_select_monomer_type_dialog");
   GtkWidget *combo_box = widget_from_builder("residue_editor_select_monomer_type_combobox");
   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo_box));
   show_restraints_editor(t);
   gtk_widget_set_visible(dialog, FALSE);

}


//   GtkWidget *dialog = widget_from_builder("residue_editor_select_monomer_type_dialog");
//   GtkWidget *combo_box = widget_from_builder("residue_editor_select_monomer_type_combobox");
//   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo_box));
//   show_restraints_editor(t);
//   gtk_widget_set_visible(dialog, FALSE);
// }


extern "C" G_MODULE_EXPORT
void
on_residue_editor_select_monomer_type_cancel_button_clicked (G_GNUC_UNUSED GtkButton       *button,
                                                             G_GNUC_UNUSED gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("residue_editor_select_monomer_type_dialog");
   gtk_widget_set_visible(dialog, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_fix_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_fixed_atom_dialog();
  gtk_widget_set_visible(w, TRUE);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_fast_sss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog;
  dialog = wrapped_create_fast_ss_search_dialog();
  gtk_widget_set_visible(dialog, TRUE);

}


extern "C" G_MODULE_EXPORT
void
on_build_na_dialog_cancelbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("build_na_dialog");
   gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_add_hydrogen_torsions_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                             gpointer        user_data) {

   GtkWidget *vbox = widget_from_builder("edit_chi_angles_vbox");

   if (gtk_check_button_get_active(checkbutton)) {
      set_find_hydrogen_torsions(1);
   } else {
      set_find_hydrogen_torsions(0);
   }
   fill_chi_angles_vbox(vbox);
}

extern "C" G_MODULE_EXPORT
void
on_mask_map_by_atom_selection_cancel_button_clicked(GtkButton *button,
                                                    gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("mask_map_by_atom_selection_dialog");
   gtk_widget_set_visible(dialog, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_mask_map_by_atom_selection_ok_button_clicked(GtkButton *button,
                                                gpointer         user_data) {

   graphics_info_t g;
   std::cout << "OK!" << std::endl;
   GtkWidget *dialog         = widget_from_builder("mask_map_by_atom_selection_dialog");
   GtkWidget *checkbutton    = widget_from_builder("mask_map_by_atom_selection_invert_checkbutton");
   GtkWidget *entry_1        = widget_from_builder("mask_map_by_atom_selection_atom_selection_entry");
   GtkWidget *entry_2        = widget_from_builder("mask_map_by_atom_selection_radius_entry");
   GtkWidget *model_combobox = widget_from_builder("mask_map_by_atom_selection_model_combobox");
   GtkWidget *map_combobox   = widget_from_builder("mask_map_by_atom_selection_map_combobox");

   std::string sel_string    = gtk_editable_get_text(GTK_EDITABLE(entry_1));
   std::string radius_string = gtk_editable_get_text(GTK_EDITABLE(entry_2));
   bool invert_flag = false;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) invert_flag = true;

   int imol_model   = g.combobox_get_imol(GTK_COMBO_BOX(model_combobox));
   int imol_for_map = g.combobox_get_imol(GTK_COMBO_BOX(map_combobox));

   // 20230427-PE maybe the radius should be in the function call?
   if (! radius_string.empty()) {
      try {
         float radius = coot::util::string_to_float(radius_string);
         g.map_mask_atom_radius = radius;
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

   mask_map_by_atom_selection(imol_for_map, imol_model, sel_string.c_str(), invert_flag);
   gtk_widget_set_visible(dialog, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_find_ligands_search_here_radiobutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  // 2023-05-07-PE this function has gone? fix this one day.
  // set_ligand_dialog_number_of_sites_sensitivity(GTK_WIDGET(button));

}


extern "C" G_MODULE_EXPORT
void
on_toolbar_multi_refine_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  toolbar_multi_refine_continue();
}


extern "C" G_MODULE_EXPORT
void
on_toolbar_multi_refine_stop_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  toolbar_multi_refine_stop();
}

extern "C" G_MODULE_EXPORT
void
on_toolbar_multi_refine_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

  toolbar_multi_refine_cancel();
}


extern "C" G_MODULE_EXPORT
void
on_map_sharpening_dialog_response() {

   // there is only one response
   GtkWidget *dialog = widget_from_builder("map_sharpening_dialog");
   gtk_widget_set_visible(dialog, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_map_sharpening_cancel_button_clicked(GtkButton       *button,
                                        gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_map_partition_by_chain_dialog_response(GtkDialog       *dialog,
                                          gint             response_id,
                                          gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      std::cout << "read the dialog - do the partitioning" << std::endl;
      GtkWidget *combobox_1 = widget_from_builder("map_partition_by_chain_map_combobox");
      GtkWidget *combobox_2 = widget_from_builder("map_partition_by_chain_model_combobox");
      int imol_model = my_combobox_get_imol(GTK_COMBO_BOX(combobox_2));
      int imol_map   = my_combobox_get_imol(GTK_COMBO_BOX(combobox_1));
      map_partition_by_chain_threaded(imol_map, imol_model);
   }

   // if (response_id == GTK_RESPONSE_CANCEL)
   // std::cout << "just close" << std::endl;
      
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_baton_build_params_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("baton_build_params_dialog");
  set_baton_build_params_from_widget(w);
  gtk_widget_set_visible(w, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_baton_build_params_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("baton_build_params_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_baton_build_set_params_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

   // GtkWidget *w = create_baton_build_params_dialog();
   GtkWidget *w = widget_from_builder("baton_build_params_dialog");
   gtk_widget_set_visible(w, TRUE);

}


extern "C" G_MODULE_EXPORT
void
on_move_molecule_here_big_molecules_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *frame = widget_from_builder("move_molecule_here_frame");
  fill_move_molecule_here_frame(frame);
}

extern "C" G_MODULE_EXPORT
void
on_python_scripting_button(GtkToggleButton *togglebutton, gpointer user_data) {
   toggle_reveal_python_scripting_entry();
}


extern "C" G_MODULE_EXPORT
void
on_environment_distances_h_bonds_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton)))
    set_show_environment_distances_h_bonds(1);
  else
    set_show_environment_distances_h_bonds(0);
}


extern "C" G_MODULE_EXPORT
void
on_environment_distances_bumps_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton)))
    set_show_environment_distances_bumps(1);
  else
    set_show_environment_distances_bumps(0);

}

extern "C" G_MODULE_EXPORT
void
on_displayed_map_style_as_lines_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                         gpointer         user_data) {

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton), "imol"));
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      set_draw_map_standard_lines(imol, 1);
      set_draw_solid_density_surface(imol, 0);
   }
}

extern "C" G_MODULE_EXPORT
void
on_map_opacity_hscale_value_changed(GtkRange        *range,
                                    gpointer         user_data) {


  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(range), "imol"));

  GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
  float fvalue = 0.01 * gtk_adjustment_get_value(adjustment);
  if (fvalue > 0.99)
    fvalue = 1.0;

  set_solid_density_surface_opacity(imol, fvalue);
}



extern "C" G_MODULE_EXPORT
void
on_remarks_browser_molecule_chooser_dialog_response(GtkDialog       *dialog,
                                                    gint             response_id,
                                                    gpointer         user_data) {

   std::cout << "here A in on_remarks_browser_molecule_chooser_dialog_response() with response_id " << response_id << std::endl;
   if (response_id == GTK_RESPONSE_OK) {
      std::cout << "here B in on_remarks_browser_molecule_chooser_dialog_response() with response_id " << response_id << std::endl;
      show_remarks_browswer(); // makes an on-the-fly dialog!
   }
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_remarks_browser_molecule_chooser_ok_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

  GtkWidget *w = widget_from_builder("remarks_browser_molecule_chooser_dialog");
  gtk_widget_set_visible(w, FALSE);
  show_remarks_browswer(); // there we look up which molecule to show.

}


extern "C" G_MODULE_EXPORT
void
on_fix_nomenclature_errors_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("fix_nomenclature_errors_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "imol"));
  fix_nomenclature_errors(imol);
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_fix_nomenclature_errors_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("fix_nomenclature_errors_dialog");
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("multi_residue_torsion_dialog");
  clear_up_moving_atoms();
  clear_pending_picks(); /* emcompasses in_multi_residue_torsion_define (but not mode) */
  clear_multi_residue_torsion_mode();
  normal_cursor();
  gtk_widget_set_visible(w, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_OK_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("multi_residue_torsion_dialog");
  gtk_widget_set_visible(w, FALSE);
  accept_regularizement();
  clear_multi_residue_torsion_mode();
}

extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_reverse_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_multi_residue_torsion_reverse_mode(1);
  else
    set_multi_residue_torsion_reverse_mode(0);

}



extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_pick_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("multi_residue_torsion_pick_dialog");
  gtk_widget_set_visible(w, FALSE);
  clear_pending_picks(); /* emcompasses in_multi_residue_torsion_mode */
  normal_cursor();
  show_multi_residue_torsion_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_pick_cancel_button_activate
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
}


extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_pick_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("multi_residue_torsion_pick_dialog");
  gtk_widget_set_visible(w, FALSE);
  clear_pending_picks(); /* emcompasses in_multi_residue_torsion_define (but not mode) */
  clear_multi_residue_torsion_mode();
  normal_cursor();
}



/* This is the call-back for the button on the Other Modelling Tools dialog. */
extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_start_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  setup_multi_residue_torsion(); // shows a dialog

}


/* wrong callback possibly */
extern "C" G_MODULE_EXPORT
void
on_keyboard_go_to_residue_entry_changed(GtkEditable     *editable,
                                        gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_generic_objects_display_all_togglebutton_toggled(GtkToggleButton *togglebutton,
                                                    gpointer         user_data) {

  int state = 0;
  if (gtk_toggle_button_get_active(togglebutton))
     state = 1;
  set_display_all_generic_objects(state);

}


extern "C" G_MODULE_EXPORT
void
on_generic_objects_close_all_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {
  close_all_generic_objects();
}


extern "C" G_MODULE_EXPORT
void
on_export_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

  GtkWidget *w = widget_from_builder("export_map_dialog");
  gtk_widget_set_visible(w, FALSE);
}

/* void */
/* on_export_map_filechooserdialog_cancel_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data) */
/* { */

/*   GtkWidget *w = widget_from_builder("export_map_filechooserdialog"); */
/*   gtk_widget_set_visible(w, FALSE); */


/* } */

/* void */
/* on_export_map_filechooserdialog_save_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data) */
/* { */

/*   GtkWidget *w = widget_from_builder("export_map_filechooserdialog"); */
/*   int imol_map = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "map_molecule_number")); */
/*   short int is_map_fragment = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "is_map_fragment")); */
/*   char *txt = (char *) g_object_get_data(G_OBJECT(w), "export_map_radius_entry_text"); */
/*   const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w)); */

/*   printf("got text: %s\n", txt); */
/*   printf("is_map_fragment: %d\n", is_map_fragment); */
/*   printf("imol_map: %d\n", imol_map); */

/*   if (is_map_fragment) {  */
/*     export_map_fragment_with_text_radius(imol_map, txt, filename); */
/*   } else {  */
/*     export_map(imol_map, filename); */
/*   }  */
/*   gtk_widget_set_visible(w, FALSE); */

/* } */


extern "C" G_MODULE_EXPORT
void
on_export_map_frame_cancel_button_clicked(GtkButton* self, gpointer user_data) {
   GtkWidget *frame = widget_from_builder("export_map_frame");
   gtk_widget_set_visible(frame, FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_export_map_frame_ok_button_clicked(GtkButton* self, gpointer user_data) {
   GtkWidget *file_chooser_dialog = widget_from_builder("export_map_file_chooser_dialog");
   GtkWidget *combobox            = widget_from_builder("export_map_map_combobox");
   GtkWidget *radius_entry        = widget_from_builder("export_map_radius_entry");
   GtkWidget *frame               = widget_from_builder("export_map_frame");
   int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   int is_map_fragment = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(frame), "is_map_fragment"));
   // std::cout << "debug:: in on_export_map_dialog_response() imol_map is " << imol_map << std::endl;
   const char *entry_text = gtk_editable_get_text(GTK_EDITABLE(radius_entry));
   // std::cout << "debug:: in on_export_map_dialog_response() got entry_text \"" << entry_text << "\"" << std::endl;
   GString* text_copy   = g_string_new(entry_text);
   gtk_widget_set_visible(GTK_WIDGET(frame), FALSE);
   gtk_widget_set_visible(file_chooser_dialog, TRUE);
   g_object_set_data(G_OBJECT(file_chooser_dialog), "map_molecule_number", GINT_TO_POINTER(imol_map));
   g_object_set_data(G_OBJECT(file_chooser_dialog), "is_map_fragment",     GINT_TO_POINTER(is_map_fragment));
   // std::cout << "debug:: in on_export_map_dialog_response() storing entry text " << text_copy << std::endl;
   g_object_set_data(G_OBJECT(file_chooser_dialog), "export_map_radius_entry_text", text_copy);
}

extern "C" G_MODULE_EXPORT
void
on_export_map_file_chooser_dialog_response (GtkDialog       *dialog,
                                                                gint             response_id,
                                                                gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) { 
      int imol_map        = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "map_molecule_number"));
      int is_map_fragment = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "is_map_fragment"));

      // std::cout << "extracted imol_map " << imol_map << " from file chooser dialog " << std::endl;
      // std::cout << "extracted is_map_fragment " << is_map_fragment << " from file chooser dialog " << std::endl;

      // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      const char *file_name = g_file_get_path(file);

      if (GTK_IS_FILE_CHOOSER(dialog)) {
         if (is_map_fragment > 0) {
            GString *txt_radius_str = static_cast<GString *>(g_object_get_data(G_OBJECT(dialog), "export_map_radius_entry_text"));
            const char *entry_text = g_string_free(txt_radius_str, FALSE); // leaking entry_text - ho hum.
            if (entry_text == 0) {
               std::cout << "ERROR:: entry_text is null " << std::endl;
            }
            export_map_fragment_with_text_radius(imol_map, entry_text, file_name);
            
         } else {
            export_map(imol_map, file_name);
         }
      }
   }
   
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_find_ligand_real_space_refine_solutions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  printf("toggled\n");
  if (gtk_toggle_button_get_active(togglebutton))
    set_find_ligand_do_real_space_refinement(1);
  else
    set_find_ligand_do_real_space_refinement(0);

}



extern "C" G_MODULE_EXPORT
void
on_replace_fragment_dialog_response(GtkDialog *dialog,
				    gint response_id,
				    gpointer user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      graphics_info_t g;
      GtkWidget *entry = widget_from_builder("replace_fragment_atom_selection_entry");
      GtkWidget *combobox_from = widget_from_builder("replace_fragment_from_molecule_combobox");
      GtkWidget *combobox_to   = widget_from_builder("replace_fragment_to_molecule_combobox");
      std::string text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      int imol_from = g.combobox_get_imol(GTK_COMBO_BOX(combobox_from));
      int imol_to   = g.combobox_get_imol(GTK_COMBO_BOX(combobox_to));
      replace_fragment(imol_to, imol_from, text.c_str()); // move this to C++ api one day
   }
   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_general_coot_molecule_chooser_dialog_response (GtkDialog       *dialog,
                                                  gpointer         user_data) {

   std::cout << "on_general_coot_molecule_chooser_dialog_response()" << std::endl;
}

extern "C" G_MODULE_EXPORT
void
on_general_coot_molecule_chooser_dialog_close (GtkDialog       *dialog,
                                               gpointer         user_data) {

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_general_coot_molecule_chooser_with_entry_and_checkbutton_dialog_close (GtkDialog       *dialog,
                                                  gpointer         user_data) {

   std::cout << "on_general_coot_molecule_chooser_with_entry_and_checkbutton_dialog_close"
             << std::endl;
}

extern "C" G_MODULE_EXPORT
void
on_general_coot_molecule_chooser_with_entry_and_checkbutton_dialog_response (GtkDialog       *dialog,
                                                                             gpointer         user_data) {

   std::cout << "on_general_coot_molecule_chooser_with_entry_and_checkbutton_dialog_response"
             << std::endl;
}



extern "C" G_MODULE_EXPORT
void
on_weight_maxtrix_estimate_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *entry = widget_from_builder("refine_params_weight_matrix_entry");
  /*  and set geometry_vs_map_weight */
  add_estimated_map_weight_to_entry(entry);

}


extern "C" G_MODULE_EXPORT
void
on_mutate_molecule_resno_1_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
   GtkWidget *res_no_1_widget = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *res_no_2_widget = widget_from_builder("mutate_molecule_resno_2_entry");
   GtkWidget *text_widget     = widget_from_builder("mutate_molecule_sequence_text");
   GtkWidget *label_widget    = widget_from_builder("mutate_residue_range_counts_label");
   mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
}


extern "C" G_MODULE_EXPORT
void
on_mutate_molecule_resno_2_entry_changed(GtkEditable     *editable,
                                         gpointer         user_data)
{
   GtkWidget *res_no_1_widget = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *res_no_2_widget = widget_from_builder("mutate_molecule_resno_2_entry");
   GtkWidget *text_widget     = widget_from_builder("mutate_molecule_sequence_text");
   GtkWidget *label_widget    = widget_from_builder("mutate_residue_range_counts_label");
   mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
}


extern "C" G_MODULE_EXPORT
void
on_display_control_last_model_only_button_clicked (GtkButton       *button,
                                                   gpointer         user_data) {
  set_only_last_model_molecule_displayed();

}


extern "C" G_MODULE_EXPORT
void
on_curlew_close_button_clicked(GtkButton *button, gpointer user_data) {

   GtkWidget *dialog = widget_from_builder("curlew_dialog");
   gtk_widget_set_visible(dialog, FALSE);

}



extern "C" G_MODULE_EXPORT
void
 on_draw_central_atom_label_activate (GMenuItem     *menuitem,
                                      gpointer         user_data) {
   std::cout << "draw central atom " << std::endl;
}



extern "C" G_MODULE_EXPORT
void
on_calculate_load_tutorial_model_and_data1_activate
                                        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  load_tutorial_model_and_data();
}

extern "C" G_MODULE_EXPORT
void
on_refine_params_geman_mcclure_alpha_combobox_changed(GtkComboBox     *combobox,
                                                      gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   try {
      float v = coot::util::string_to_float(t);
      set_refinement_geman_mcclure_alpha(v); // in cc-interface.hh
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_lennard_jones_epsilon_combobox_changed(GtkComboBox     *combobox,
                                                        gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox); // save it for set active item next time
   set_refinement_lennard_jones_epsilon_from_text(active_item_idx, t);
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_rama_restraints_weight_combobox_changed(GtkComboBox     *combobox,
                                                         gpointer         user_data) {

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_ramachandran_restraints_weight_from_text(active_item_idx, t);

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
on_refine_params_overall_weight_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   if (t)
      std::cout << "in on_refine_params_overall_weight_combobox_changed() " << t << std::endl;
   else
      std::cout << "in on_refine_params_overall_weight_combobox_changed() t was null "  << std::endl;
   set_refinement_overall_weight_from_text(t);
}




extern "C" G_MODULE_EXPORT
void
on_refine_params_more_control_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   if (togglebutton) {
      GtkWidget *frame = widget_from_builder("refine_params_more_control_frame");
   }
}


extern "C" G_MODULE_EXPORT
void
on_accept_reject_flip_this_peptide_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  pepflip_intermediate_atoms();
}


extern "C" G_MODULE_EXPORT
void
on_accept_reject_flip_next_peptide_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  pepflip_intermediate_atoms_other_peptide();
}


extern "C" G_MODULE_EXPORT
void
on_accept_reject_crankshaft_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  crankshaft_peptide_rotation_optimization_intermediate_atoms();
}


extern "C" G_MODULE_EXPORT
void
on_accept_reject_backrub_rotamer_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  backrub_rotamer_intermediate_atoms();
}

extern "C" G_MODULE_EXPORT
void
on_symmetry_always_on_checkbutton_toggled (GtkCheckButton *checkbutton,
					   gpointer        user_data) {

   GtkWidget *symmetry_on_radio_button = NULL;
   if (gtk_check_button_get_active(checkbutton)) {
      add_symmetry_on_to_preferences_and_apply();
      symmetry_on_radio_button = widget_from_builder("show_symmetry_yes_radiobutton");
      if (! gtk_check_button_get_active(GTK_CHECK_BUTTON(symmetry_on_radio_button)))
	 gtk_check_button_set_active(GTK_CHECK_BUTTON(symmetry_on_radio_button), TRUE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_yes_radiobutton_toggled(GtkCheckButton *checkbutton,
                                          gpointer       user_data) {

   // we don't need a callback for the "no" button - just to the opposite

   std::cout << "on_show_symmetry_yes_radiobutton_toggled() "
             << gtk_check_button_get_active(checkbutton) << std::endl;

  if (gtk_check_button_get_active(checkbutton)) {
      set_show_symmetry_master(1);
  } else {
     set_show_symmetry_master(0);
  }
}

extern "C" G_MODULE_EXPORT
void
on_symmetry_radius_entry_activate(GtkEntry* self,
                                  gpointer user_data) {

   const char *text = gtk_editable_get_text(GTK_EDITABLE(self));
   if (text) {
      std::string t(text);
      try {
         float f = coot::util::string_to_float(t);
         set_symmetry_size(f);
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
   }
}

extern "C" G_MODULE_EXPORT
void
show_symmetry_switch_state_set(GtkSwitch *switch_widget,
                               gboolean   state,
                               gpointer   user_data) {

   if (state)
      set_show_symmetry_master(1);
  else
      set_show_symmetry_master(0);
}

// #ifdef FIX_THE_KEY_PRESS_EVENTS
// extern "C" G_MODULE_EXPORT
// gboolean
// on_symmetry_radius_entry_key_release_event(GtkWidget       *widget,
//                                            GdkEventKey     *event,
//                                            gpointer         user_data) {

//    const char *text;
//    if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
//       text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(widget)));
//       set_symmetry_size_from_widget(text);
//    }
//    return TRUE;
// }
// #endif


extern "C" G_MODULE_EXPORT
void
on_hscale_symmetry_colour_value_changed(GtkRange        *range,
                                        gpointer         user_data) {

   gdouble f = gtk_range_get_value(range);
   set_symmetry_colour_merge(f);
}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_expanded_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton)))
     set_symmetry_atom_labels_expanded(1);
   else
     set_symmetry_atom_labels_expanded(0);
}



extern "C" G_MODULE_EXPORT
void
on_simple_refmac_dialog_response       (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   if (response_id == GTK_RESPONSE_CLOSE) {
      /* do I need to do this? */
      /* gtk_widget_set_visible(dialog, FALSE); */
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_OK) {
      simple_refmac_run_refmac(GTK_WIDGET(dialog));
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_simple_refmac_mtz_file_button_clicked (GtkButton       *button,
                                          gpointer         user_data) {

   GtkWidget *w = widget_from_builder("simple_refmac_filechooser_dialog");
   GtkWidget *simple_refmac_dialog = widget_from_builder("simple_refmac_dialog");
   /* automtically file filter only mtz files */
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.mtz");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(w), filterselect);
   /* we need to find the file combo box in the simple refmac dialog */
   g_object_set_data(G_OBJECT(w), "simple_refmac_dialog", simple_refmac_dialog);
   gtk_widget_set_visible(w, TRUE);
}

extern "C" G_MODULE_EXPORT
void
on_simple_refmac_filechooser_dialog_response (GtkDialog       *dialog,
                                              gint             response_id,
                                              gpointer         user_data) {

   // const gchar *file_name = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

   GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
   GError *error = NULL;
   GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                            G_FILE_QUERY_INFO_NONE, NULL, &error);
   const char *file_name = g_file_info_get_name(file_info);

   if (response_id == GTK_RESPONSE_CLOSE) {
      std::cout << "on_simple_refmac_filechooserdialog_response() Close\n";
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      std::cout << "on_simple_refmac_filechooserdialog_response() Cancel\n";
   }

   if (response_id == GTK_RESPONSE_OK) {
      // GtkWidget *simple_refmac_dialog = GTK_WIDGET(g_object_get_data(G_OBJECT(dialog), "simple_refmac_dialog"));
      GtkWidget *simple_refmac_dialog = widget_from_builder("simple_refmac_dialog");
      if (simple_refmac_dialog) {
         GtkWidget *file_combobox = widget_from_builder("simple_refmac_mtz_file_combobox");
         if (file_combobox) {
            gtk_combo_box_text_remove_all(GTK_COMBO_BOX_TEXT(file_combobox));
            gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(file_combobox), file_name);
            gtk_combo_box_set_active(GTK_COMBO_BOX(file_combobox), 0);
         }
      }
   }

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);

}


void
handle_map_properties_specularity_change(int imol, GtkWidget *checkbutton) {

   molecule_class_info_t &m = graphics_info_t::molecules[imol];

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) {
      // std::cout << "Turn on specularity " << std::endl;
      GtkWidget *strength_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton),  "strength_entry"));
      GtkWidget *shininess_entry = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton), "shininess_entry"));
      if (! strength_entry)  return;
      if (! shininess_entry) return;
      std::string strength_entry_text  = gtk_editable_get_text(GTK_EDITABLE(strength_entry));
      std::string shininess_entry_text = gtk_editable_get_text(GTK_EDITABLE(shininess_entry));
      float f1 = coot::util::string_to_float(strength_entry_text);
      float f2 = coot::util::string_to_float(shininess_entry_text);
      m.material_for_maps.specular_strength = f1;
      m.material_for_maps.shininess         = f2;
      m.material_for_maps.turn_specularity_on(true);
      std::cout << "in handle_map_properties_specularity_change() imol: " << imol << " do: " <<  m.material_for_maps.do_specularity
                << " strength " << m.material_for_maps.specular_strength << " shiny " << m.material_for_maps.shininess << std::endl;
   } else {
      // std::cout << "Turn off specularity " << std::endl;
      m.material_for_maps.turn_specularity_on(false);
   }
   graphics_draw();
}

extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_specularity_strength_entry_activate(GtkEntry* self, gpointer user_data) {

   std::cout << "strength entry key press activate" << std::endl;
   GtkWidget *checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "specularity_checkbutton"));
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(self), "imol"));
   handle_map_properties_specularity_change(imol, checkbutton);
}

extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_specularity_shininess_entry_activate(GtkEntry* self, gpointer user_data) {

   std::cout << "shininess entry key press activate" << std::endl;
   GtkWidget *checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "specularity_checkbutton"));
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(self), "imol"));
   handle_map_properties_specularity_change(imol, checkbutton);
}




extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_filechooser_dialog_response(GtkDialog       *dialog,
                                              gint             response_id,
                                              gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      //  const char *fnc = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);

      if (file_name) {
         read_cif_dictionary(file_name);
      }
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_validation_graph_model_combobox_changed(GtkComboBox* self, gpointer user_data) {

   GtkTreeIter iter;
   if (gtk_combo_box_get_active_iter(self, &iter)) {
      int new_active_model;
      gtk_tree_model_get(gtk_combo_box_get_model(self),&iter,1,&new_active_model,-1);
      graphics_info_t::update_active_validation_graph_model(new_active_model);
   } else {
      // this is noisy - and it seems possible that the active iter has not yet been set.
      if (false) {
         std::string mess = "on_validation_graph_model_combobox_changed(): ";
         mess += "Could not get active iter in validation graph model ComboBox";
         g_warning(mess.c_str());
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_validation_graph_chain_id_combobox_changed(GtkComboBoxText* self, gpointer user_data) {
   auto chain_id = std::string(gtk_combo_box_text_get_active_text(self));
   graphics_info_t::change_validation_graph_chain(chain_id);
}

void
on_validation_graph_checkbutton_toggled(GtkCheckButton* self, coot::validation_graph_type graph_type) {

   if (gtk_check_button_get_active(self)) {
      graphics_info_t g;
      // read imol from the widget, but now now, let's use active_validation_graph_model_idx
      int imol = g.active_validation_graph_model_idx;
      graphics_info_t::create_validation_graph(imol, graph_type);
   } else {
      graphics_info_t::destroy_validation_graph(graph_type);
   }
}

extern "C" G_MODULE_EXPORT
void
on_density_fit_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::density_fit);
}

extern "C" G_MODULE_EXPORT
void
on_temp_factor_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::temp_factor);
}

extern "C" G_MODULE_EXPORT
void
on_rota_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::rota);
}

extern "C" G_MODULE_EXPORT
void
on_rama_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::rama);
}

extern "C" G_MODULE_EXPORT
void
on_omega_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::omega);
}

extern "C" G_MODULE_EXPORT
void
on_geometry_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::geometry);
}

extern "C" G_MODULE_EXPORT
void
on_ncs_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::ncs);
}

extern "C" G_MODULE_EXPORT
void
on_density_correlation_graph_toggled(GtkCheckButton* self, gpointer user_data) {
   on_validation_graph_checkbutton_toggled(self,coot::validation_graph_type::density_correlation);
}

extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_molecule_chooser_ok_button_clicked(GtkButton       *button,
                                                        gpointer         user_data) {

   GtkWidget *dialog          = widget_from_builder("ramachandran_plot_molecule_chooser_dialog");
   GtkWidget *combobox        = widget_from_builder("ramachandran_plot_molecule_chooser_model_combobox");
   GtkWidget *selection_entry = widget_from_builder("ramachandran_plot_molecule_chooser_residue_selection_entry");
   GtkWidget *scrolled        = widget_from_builder("ramachandran_plots_scrolled_window");
   GtkWidget *pane            = widget_from_builder("main_window_ramchandran_and_validation_pane");

   std::string residue_selection_string = gtk_editable_get_text(GTK_EDITABLE(selection_entry));

   GtkTreeIter iter;
   if (gtk_combo_box_get_active_iter(GTK_COMBO_BOX(combobox), &iter)) {
      int imol_active;
      gtk_tree_model_get(gtk_combo_box_get_model(GTK_COMBO_BOX(combobox)),&iter,1,&imol_active,-1);
      // imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox)); // but use Jakub-style comboboxes
      show_opengl_ramachandran_plot(imol_active, residue_selection_string);
      gtk_widget_set_visible(dialog, FALSE);
      gtk_widget_set_visible(scrolled, TRUE);

      // Make sure that the pane is big enough
      int pos = gtk_paned_get_position(GTK_PANED(pane));
      if (pos < 200)
         gtk_paned_set_position(GTK_PANED(pane), 480);
   } else {
      std::cout << "ERROR:: on_ramachandran_plot_molecule_chooser_ok_button_clicked() get active iter failed"
                << std::endl;
   }
   graphics_info_t::graphics_grab_focus();
}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_molecule_chooser_cancel_button_clicked (GtkButton       *button,
                                                             gpointer         user_data) {
   GtkWidget *w = widget_from_builder("ramachandran_plot_molecule_chooser_dialog");
   gtk_widget_set_visible(w, FALSE);
   graphics_info_t::graphics_grab_focus();

}



extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_specularity_state_checkbutton_toggled(GtkCheckButton *checkbutton,
                                                               gpointer         user_data) {

   // was it set?
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "imol"));
   handle_map_properties_specularity_change(imol, GTK_WIDGET(checkbutton));
   graphics_info_t::graphics_grab_focus();

}

// ----------------------------------- updating maps -----

extern "C" G_MODULE_EXPORT
void
on_updating_maps_cancel_button_clicked(GtkButton       *button,
                                       gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("updating_maps_dialog");
   gtk_widget_set_visible(dialog, FALSE);
   graphics_info_t::graphics_grab_focus();

}


extern "C" G_MODULE_EXPORT
void
on_updating_maps_ok_button_clicked(GtkButton       *button,
                                   gpointer         user_data) {

   GtkWidget *model_combobox    = widget_from_builder("updating_maps_model_combobox");
   GtkWidget *map_combobox      = widget_from_builder("updating_maps_map_combobox");
   GtkWidget *diff_map_combobox = widget_from_builder("updating_maps_diff_map_combobox");
   GtkWidget *check_button      = widget_from_builder("updating_maps_auto_update_checkbutton");

   int imol          = my_combobox_get_imol(GTK_COMBO_BOX(model_combobox));
   int imol_map      = my_combobox_get_imol(GTK_COMBO_BOX(map_combobox));
   int imol_diff_map = my_combobox_get_imol(GTK_COMBO_BOX(diff_map_combobox));

   bool auto_update_flag = false;
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(check_button))) auto_update_flag = true;

   if (auto_update_flag) {
      set_auto_updating_sfcalc_genmap(imol, imol_map, imol_diff_map);
   } else {
      calculate_maps_and_stats_py(imol, imol_map, imol_map, imol_diff_map);
   }

   GtkWidget *dialog = widget_from_builder("updating_maps_dialog");
   gtk_widget_set_visible(dialog, FALSE);

   GtkWidget *points_button = widget_from_builder("coot-points-button");
   gtk_widget_set_visible(points_button, TRUE);
   graphics_info_t::graphics_grab_focus();

}


extern "C" G_MODULE_EXPORT
void
on_ligand_check_dialog_close_button_clicked(GtkButton       *button,
                                            gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("ligand_check_dialog");
   gtk_widget_set_visible(dialog, FALSE);
   graphics_info_t::graphics_grab_focus();

}

extern "C" G_MODULE_EXPORT
void
on_generic_validation_box_of_buttons_close_button_clicked(GtkButton       *button,
                                                          gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("generic_validation_box_of_buttons_dialog");
   GtkWidget *box = widget_from_builder("generic_validation_box_of_buttons_box");
   if (box) {
      graphics_info_t g;
      g.clear_out_container(box);
   }
   gtk_widget_set_visible(dialog, FALSE);
   graphics_info_t::graphics_grab_focus();
}

extern "C" G_MODULE_EXPORT
void
on_download_monomers_cancel_button_clicked(GtkButton       *button,
					   gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("download_monomers_dialog");
   gtk_widget_set_visible(dialog, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_download_monomers_ok_button_clicked(GtkButton       *button,
				       gpointer         user_data) {

   GtkWidget *vbox = widget_from_builder("download_monomers_dialog_vbox_inner");
   if (vbox) {
      GtkWidget *item_widget = gtk_widget_get_first_child(vbox);
      while (item_widget) {
	 gchar *comp_id = static_cast<gchar *>(g_object_get_data(G_OBJECT(item_widget), "comp_id"));
         std::cout << "debug:: on_download_monomers_ok_button_clicked comp_id is " << comp_id << std::endl;
	 if (comp_id) {
	    GtkWidget *dialog = widget_from_builder("download_monomers_dialog");
	    int run_get_monomer_post_fetch_flag =
	       GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "run_get_monomer_post_fetch_flag"));
	    get_monomer_dictionary_in_subthread(comp_id, run_get_monomer_post_fetch_flag);
	 }
	 item_widget = gtk_widget_get_next_sibling(item_widget);
      };

   }

   GtkWidget *dialog = widget_from_builder("download_monomers_dialog");
   gtk_widget_set_visible(dialog, FALSE);
}


extern "C" G_MODULE_EXPORT
void
on_add_other_solvent_molecules_new_residue_type_button_clicked(GtkButton  *button,
                                                               gpointer    user_data) {

   GtkWidget *entry = widget_from_builder("add_other_solvent_molecules_new_residue_type_entry");
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(entry), "imol"));
   std::cout << "DEBUG:: transfering imol " << imol << std::endl;
   g_object_set_data(G_OBJECT(entry), "imol", GINT_TO_POINTER(imol));
   gtk_widget_set_visible(entry, TRUE);
}

#include "get-monomer.hh"

extern "C" G_MODULE_EXPORT
void on_add_other_solvent_molecules_new_residue_type_entry_activate(GtkEntry *entry,
                                                                    gpointer  user_data) {

   // this callback is similar to the callback for the build-in residue type
   // buttons in the dialog
   const char *s = gtk_editable_get_text(GTK_EDITABLE(entry));
   graphics_info_t g;
   std::string type(s);
   int imol_ligand = get_monomer(type);
   fit_to_map_by_random_jiggle(imol_ligand, "A", 1, "", 100, 2.0);
   if (g.is_valid_model_molecule(imol_ligand)) {
      coot::residue_spec_t rspec("A", 1, "");
      mmdb::Residue *residue_p = g.molecules[imol_ligand].get_residue(rspec);
      if (residue_p) {
         mmdb::Manager *mol = g.molecules[imol_ligand].atom_sel.mol;
         std::vector<mmdb::Residue *> v;
         std::string alt_conf;
         v.push_back(residue_p);
         short int save_state = g.refinement_immediate_replacement_flag;
         g.refinement_immediate_replacement_flag = 1;
         coot::refinement_results_t rr = g.refine_residues_vec(imol_ligand, v, alt_conf, mol);
         g.refinement_immediate_replacement_flag = save_state;
         c_accept_moving_atoms();
         delete_hydrogen_atoms(imol_ligand);
         int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(entry), "imol"));
         if (g.is_valid_model_molecule(imol)) {
            std::vector<atom_selection_container_t> add_molecules_at_sels;
            add_molecules_at_sels.push_back(g.molecules[imol_ligand].atom_sel);
            g.molecules[imol].merge_molecules(add_molecules_at_sels);
            close_molecule(imol_ligand);
         }
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_add_other_solvent_molecules_close_button_clicked(GtkButton       *button,
						    gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("add_other_solvent_molecules_dialog");
   if (dialog) {
      gtk_widget_set_visible(dialog, FALSE);
   }
}

extern "C" G_MODULE_EXPORT
void
on_first_startup_cancel_button_clicked(GtkButton       *button,
                                       gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("first-startup-dialog");
   gtk_widget_set_visible(dialog, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_first_startup_use_left_button_clicked(GtkButton       *button,
                                       gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("first-startup-dialog");
   gtk_widget_set_visible(dialog, FALSE);

   preferences_internal_change_value_int(PREFERENCES_VIEW_ROTATION_MOUSE_BUTTON, 1);
   set_use_primary_mouse_button_for_view_rotation(1);
}

extern "C" G_MODULE_EXPORT
void
on_first_startup_use_right_button_clicked(GtkButton       *button,
                                      gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("first-startup-dialog");
   gtk_widget_set_visible(dialog, FALSE);

   preferences_internal_change_value_int(PREFERENCES_VIEW_ROTATION_MOUSE_BUTTON, 0);
   set_use_primary_mouse_button_for_view_rotation(0);
}

extern "C" G_MODULE_EXPORT
void
on_material_lighting_ambient_colorbutton_color_set(GtkColorButton *colorbutton,
                                                   gpointer        user_data) {

   GdkRGBA rgba;
   gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(colorbutton), &rgba);
   GtkWidget *combobox = widget_from_builder("material_lighting_molecule_comboboxtext");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      glm::vec4 ambient(rgba.red, rgba.green, rgba.blue, 1.0f);
      g.molecules[imol].material_for_models.ambient = ambient;
      g.molecules[imol].model_molecule_meshes.set_material_ambient(ambient);
      g.graphics_draw();
   }
}

extern "C" G_MODULE_EXPORT
void
on_material_lighting_diffuse_colorbutton_color_set(GtkColorButton *colorbutton,
                                                   gpointer        user_data) {

   GdkRGBA rgba;
   gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(colorbutton), &rgba);
   GtkWidget *combobox = widget_from_builder("material_lighting_molecule_comboboxtext");
   int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      glm::vec4 diffuse(rgba.red, rgba.green, rgba.blue, 1.0f);
      g.molecules[imol].material_for_models.ambient = diffuse;
      g.molecules[imol].model_molecule_meshes.set_material_diffuse(diffuse);
      g.graphics_draw();
   }
}

extern "C" G_MODULE_EXPORT
void
on_button_clicked(GtkButton       *button,
                  gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("ligand_check_dialog");
   gtk_widget_set_visible(dialog, FALSE);

}

