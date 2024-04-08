/*
 * src/single-map-properties-dialog.hh
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

//
#ifndef SINGLE_MAP_PROPERTIES_HH
#define SINGLE_MAP_PROPERTIES_HH

#include <utility>
#include <gtk/gtk.h>

std::pair<GtkWidget *, GtkBuilder *> create_single_map_properties_dialog_gtk3();

void fill_single_map_properties_dialog_gtk3(std::pair<GtkWidget *, GtkBuilder *> w_and_b, int imol);
GtkWidget *wrapped_create_single_map_properties_dialog_gtk3(int imol);

// and the callbacks

void on_displayed_map_style_as_lines_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data);
void on_displayed_map_style_as_transparent_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data);
void on_displayed_map_style_as_cut_glass_radiobutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data);
void on_map_opacity_hscale_value_changed_gtkbuilder_callback(GtkRange *range);
void on_single_map_properties_contour_level_apply_button_clicked_gtkbuilder_callback(GtkButton *button, gpointer user_data);
void on_single_map_sigma_checkbutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data);
void on_map_properties_dialog_specularity_state_checkbutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data);
void on_map_properties_dialog_fresnel_state_checkbutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data);
void on_map_properties_dialog_specularity_strength_entry_activate_gtkbuilder_callback(GtkEntry *entry, gpointer user_data);
void on_map_properties_dialog_specularity_shininess_entry_activate_gtkbuilder_callback(GtkEntry *entry, gpointer user_data);
void on_map_properties_dialog_fresnel_bias_entry_activate_gtkbuilder_callback(GtkEntry *entry, gpointer user_data);
void on_map_properties_dialog_fresnel_scale_entry_activate_gtkbuilder_callback(GtkEntry *entry, gpointer user_data);
void on_map_properties_dialog_fresnel_power_entry_activate_gtkbuilder_callback(GtkEntry *entry, gpointer user_data);
extern "C" G_MODULE_EXPORT void on_single_map_properties_colour_button_clicked_gtkbuilder_callback(GtkButton *button, gpointer user_data);
extern "C" G_MODULE_EXPORT void on_single_map_properties_ok_button_clicked_gtkbuilder_callback(GtkButton *button, gpointer user_data);

#endif // SINGLE_MAP_PROPERTIES_HH
