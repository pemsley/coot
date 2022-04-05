/* src/gtk-manual.hh
 *
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Copyright 2008, 2009 by The University of Oxford
 * Author: Paul Emsley
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

#ifndef GTK_MANUAL_H
#define GTK_MANUAL_H

#include <gtk/gtk.h>
#ifndef HAVE_SUPPORT_H
#define HAVE_SUPPORT_H
#include "support.h"
#endif /* HAVE_SUPPORT_H */

#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif
#endif

BEGIN_C_DECLS


#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


void
on_map_color_changed(GtkWidget *w,   gpointer tmd);
/* 		 GtkColorSelection *cs); */

void
on_map_col_sel_ok_button_clicked(GtkButton       *button,
                                 gpointer         user_data);

void
on_map_col_sel_cancel_button_clicked(GtkButton       *button,
                                     gpointer         *tmd);


/* Map Colour */

struct map_colour_data_type {
   int imol;
#if (GTK_MAJOR_VERSION >= 4) || GTK_DISABLE_DEPRECATED
   GtkColorChooser* color_chooser;
#else
   GtkColorSelection *color_selection;
#endif

};

GtkWidget* create_map_colour_selection_window(struct map_colour_data_type *mcdt);


/* Symmetry Colour */
#if GTK_MAJOR_VERSION >=4 || GTK_DISABLE_DEPRECATED
void on_symmetry_color_changed(GtkWidget *w, GtkColorChooser *cs);
#else
void on_symmetry_color_changed(GtkWidget *w, GtkColorSelection *cs);
#endif

void
on_symm_col_sel_ok_button_clicked(GtkButton       *button,
                                  gpointer         user_data);

void
on_symm_col_sel_cancel_button_clicked(GtkButton       *button,
                                      gpointer         user_data);

GtkWidget *
create_symmetry_colour_selection_window();


void
create_initial_map_color_submenu(GtkWidget *widget);

void
update_map_colour_menu_manual(int imol, const char *label);

void
my_map_colour_activate (GtkMenuItem     *menuitem,
                        gpointer         user_data);

/* similar stuff for the scroll wheel */
void
create_initial_map_scroll_wheel_submenu(GtkWidget *widget);


void
update_map_scroll_wheel_menu_manual(int imol, const char *name);

void
my_map_scroll_wheel_activate(GtkMenuItem     *menuitem,
			     gpointer         user_data);

/* And similar for ramachandran plot: */
void create_initial_ramachandran_mol_submenu(GtkWidget *widget);
void update_ramachandran_plot_menu_manual(int imol, const char *name);
void rama_plot_mol_selector_activate (GtkMenuItem     *menuitem,
				      gpointer         user_data);

/* And similar for sequence view: */
void create_initial_sequence_view_mol_submenu(GtkWidget *widget);
void update_sequence_view_menu_manual(int imol, const char *name);
void sequence_view_mol_selector_activate (GtkMenuItem     *menuitem,
					  gpointer         user_data);


/* skeleton colour */


GtkWidget *create_skeleton_colour_selection_window();

#if GTK_MAJOR_VERSION >=4 || GTK_DISABLE_DEPRECATED
void on_skeleton_color_changed(GtkWidget *w, GtkColorChooser *colorsel);
#else
void on_skeleton_color_changed(GtkWidget *w, GtkColorSelection *colorsel);
#endif

void on_skeleton_col_sel_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data);

void on_skeleton_col_sel_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data);

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*              map and molecule display control                            */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */



GtkWidget *get_radio_button_in_scroll_group(int imol_this);

void update_name_in_display_control_molecule_combo_box(GtkWidget *display_control_window_glade, const gchar *name, int n); 
void simple_display_control_mol_menu_item(GtkWidget *model_menu, int imol_no, int map_coords_mol_flag);

void render_as_bonds_button_select(int imol);
void render_as_ca_bonds_button_select(int imol);
void render_as_ca_plus_ligands_bonds_button_select(int imol);
void render_as_bonds_colored_by_chain_button_select(int imol);
void render_as_bonds_goodsell_colored_by_chain_button_select(int imol);
void render_as_bonds_colored_by_molecule_button_select(int imol);
void render_as_bonds_no_waters(int imol);
void render_as_ca_plus_ligands_sec_str_bonds_button_select(int imol);
void render_as_sec_struct_bonds_button_select(int imol);
void render_as_rainbow_representation_button_select(int imol);
void render_as_b_factor_representation_button_select(int imol);
void render_as_b_factor_cas_representation_button_select(int imol);
void render_as_occupancy_representation_button_select(int imol);


/* void */
/* on_display_control_map_displayed_button_clicked   (GtkButton       *button, */
/* 						   gpointer         user_data);  */

/* map */
void
on_display_control_map_displayed_button_toggled   (GtkToggleButton       *button,
						   gpointer         user_data);
void
on_display_control_map_properties_button_clicked   (GtkButton       *button,
						   gpointer         user_data);

void display_control_add_delete_molecule_button(int imol, GtkWidget *hbox32,
						short int is_map_molecule);

void
on_display_control_delete_molecule_button_clicked   (GtkButton       *button,
						   gpointer         user_data);

/* molecule */
void
on_display_control_mol_displayed_button_toggled   (GtkToggleButton       *button,
						   gpointer         user_data);
void
on_display_control_mol_active_button_toggled   (GtkToggleButton       *button,
						gpointer         user_data);
void
on_display_control_mol_properties_button_toggled   (GtkButton       *button,
						    gpointer         user_data);
void
on_display_control_map_scroll_radio_button_toggled (GtkToggleButton    *button,
						    gpointer         user_data);

void
on_display_control_map_scroll_radio_button_group_changed (GtkRadioButton *button,
							  gpointer         user_data);


void fill_map_colour_patch(GtkWidget *patch_frame, int imol);


GtkWidget *selections_and_colours_combobox();


/* ------------------------------------------------------------------------------ */
/*              cell chooser for phs                                              */
/* ------------------------------------------------------------------------------ */


GSList *display_cell_chooser_box(GtkWidget *phs_cell_choice_window,
			      GSList *phs_cell_group, int n);

void display_none_cell_chooser_box(GtkWidget *phs_cell_choice_window,
				   GSList *phs_cell_group);


void
on_add_rep_all_on_check_button_toggled   (GtkToggleButton       *button,
					  gpointer         user_data);



/* ------------------------------------------------------------------------------ */
/*  manual splash screen                                                          */
/* ------------------------------------------------------------------------------ */

GtkWidget* create_splash_screen_window_for_file(const char *file_name);


/* ------------------------------------------------------------------------------ */
/*  get the imol from a combobox.                                                 */
/* ------------------------------------------------------------------------------ */

// maybe this needs its own header?
int my_combobox_get_imol(GtkComboBox *combobox);



END_C_DECLS

#endif /* GTK_MANUAL_H */
