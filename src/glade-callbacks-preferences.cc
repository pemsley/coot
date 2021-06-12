/* src/glade-callbacks-preferences.cc
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

// from support.h
// GtkWidget* lookup_widget (GtkWidget *widget, const gchar *widget_name);
#include "support.h"

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;


extern "C" G_MODULE_EXPORT
void
on_preferences_general_radiotoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_GENERAL_PREFERENCES);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bond_radiotoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_BOND_PREFERENCES);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_radiotoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_MAP_PREFERENCES);
}



extern "C" G_MODULE_EXPORT
void
on_preferences_geometry_radiotoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_GEOMETRY_PREFERENCES);
}




extern "C" G_MODULE_EXPORT
void
on_preferences_colour_radiotoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_COLOUR_PREFERENCES);
}




extern "C" G_MODULE_EXPORT
void
on_preferences_other_radiotoolbutton_toggled_gtkbuilder_callback
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_OTHER_PREFERENCES);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_ok_button_clicked_gtkbuilder_callback       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "preferences");
  save_preferences();
  gtk_widget_destroy(w);
  clear_preferences();
}

extern "C" G_MODULE_EXPORT
void
on_preferences_reset_button_clicked_gtkbuilder_callback    (GtkButton       *button,
                                        gpointer         user_data)
{
  reset_preferences();

}

extern "C" G_MODULE_EXPORT
void
on_preferences_destroy_gtkbuilder_callback                 (GtkWidget       *object,
                                        gpointer         user_data)
{
  clear_preferences();
}

extern "C" G_MODULE_EXPORT
void
on_preferences_geometry_cis_peptide_bad_yes_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MARK_CIS_BAD, 1);
    set_mark_cis_peptides_as_bad(1);
  }
}

extern "C" G_MODULE_EXPORT
void
on_preferences_geometry_cis_peptide_bad_no_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MARK_CIS_BAD, 0);
    set_mark_cis_peptides_as_bad(0);
  }
}

extern "C" G_MODULE_EXPORT
void
on_preferences_bond_colours_hscale_value_changed_gtkbuilder_callback
                                        (GtkRange        *range,
                                        gpointer         user_data)
{
  GtkAdjustment *adjustment;
  float fvalue;
  adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
  fvalue = gtk_adjustment_get_value(adjustment);
  preferences_internal_change_value_float(PREFERENCES_BOND_COLOURS_MAP_ROTATION, fvalue);
  set_colour_map_rotation_on_read_pdb(fvalue);
}

extern "C" G_MODULE_EXPORT
void
on_preferences_bond_colours_checkbutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_BOND_COLOUR_ROTATION_C_ONLY, 1);
    set_colour_map_rotation_on_read_pdb_c_only_flag(1);
  } else {
    preferences_internal_change_value_int(PREFERENCES_BOND_COLOUR_ROTATION_C_ONLY, 0);
    set_colour_map_rotation_on_read_pdb_c_only_flag(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_black_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, 0, 0, 0);
    set_background_colour(0, 0, 0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_white_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, 1, 1, 1);
    set_background_colour(1, 1, 1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_own_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  GtkWidget *w;
  GdkColor bg_colour;
  float fval1;
  float fval2;
  float fval3;
  w = lookup_widget(GTK_WIDGET(togglebutton), "preferences_bg_colour_colorbutton");
  gtk_color_button_get_color(GTK_COLOR_BUTTON(w), &bg_colour);
  fval1 = (float)bg_colour.red / 65535;
  fval2 = (float)bg_colour.green / 65535;
  fval3 = (float)bg_colour.blue / 65535;

  preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, fval1, fval2, fval3);
  set_background_colour(fval1, fval2, fval3);

}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_colorbutton_color_set_gtkbuilder_callback
                                        (GtkColorButton  *colorbutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  w = lookup_widget(GTK_WIDGET(colorbutton), "preferences_bg_colour_own_radiobutton");
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
    GdkColor bg_colour;
    float fval1;
    float fval2;
    float fval3;
    gtk_color_button_get_color(colorbutton, &bg_colour);
    fval1 = (float)bg_colour.red / 65535;
    fval2 = (float)bg_colour.green / 65535;
    fval3 = (float)bg_colour.blue / 65535;

    preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, fval1, fval2, fval3);
    set_background_colour(fval1, fval2, fval3);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_bg_colour_colorbutton_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w;
  w = lookup_widget(GTK_WIDGET(button), "preferences_bg_colour_own_radiobutton");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);

}



extern "C" G_MODULE_EXPORT
void
on_preferences_map_radius_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
					 gpointer         user_data)
{
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval > 0) && (fval <1000)) {
    preferences_internal_change_value_float(PREFERENCES_MAP_RADIUS, fval);
    set_map_radius(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_radius_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
					 gpointer         user_data)
{
  GtkEntry *entry;
  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable), "preferences_map_radius_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval > 0) && (fval <200)) {
    preferences_internal_change_value_float(PREFERENCES_MAP_RADIUS, fval);
    set_map_radius(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_increment_size_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
					 gpointer         user_data)
{
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if (fval > 0) {
    preferences_internal_change_value_float(PREFERENCES_MAP_ISOLEVEL_INCREMENT, fval);
    set_iso_level_increment(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_increment_size_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry *entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable),
					    "preferences_map_increment_size_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if (fval > 0) {
    preferences_internal_change_value_float(PREFERENCES_MAP_ISOLEVEL_INCREMENT, fval);
    set_iso_level_increment(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_diff_increment_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if (fval > 0) {
    preferences_internal_change_value_float(PREFERENCES_DIFF_MAP_ISOLEVEL_INCREMENT, fval);
    set_diff_map_iso_level_increment(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_diff_increment_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry *entry;
  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable), "preferences_map_diff_increment_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if (fval > 0) {
    preferences_internal_change_value_float(PREFERENCES_DIFF_MAP_ISOLEVEL_INCREMENT, fval);
    set_diff_map_iso_level_increment(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_sampling_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval < 100) && (fval > 1)) {
    preferences_internal_change_value_float(PREFERENCES_MAP_SAMPLING_RATE, fval);
    set_map_sampling_rate(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_sampling_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry *entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable), "preferences_map_sampling_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval < 100) && (fval > 1)) {
    preferences_internal_change_value_float(PREFERENCES_MAP_SAMPLING_RATE, fval);
    set_map_sampling_rate(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_dynamic_sampling_checkbutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_DYNAMIC_MAP_SAMPLING, 1);
    set_dynamic_map_sampling_on();
  } else {
    preferences_internal_change_value_int(PREFERENCES_DYNAMIC_MAP_SAMPLING, 0);
    set_dynamic_map_sampling_off();
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_dynamic_size_checkbutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_DYNAMIC_MAP_SIZE_DISPLAY, 1);
    set_dynamic_map_size_display_on();
  } else {
    preferences_internal_change_value_int(PREFERENCES_DYNAMIC_MAP_SIZE_DISPLAY, 0);
    set_dynamic_map_size_display_off();
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_diff_map_colours_coot_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SWAP_DIFF_MAP_COLOURS, 0);
    set_swap_difference_map_colours(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_diff_map_colours_o_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SWAP_DIFF_MAP_COLOURS, 1);
    set_swap_difference_map_colours(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_colours_hscale_value_changed_gtkbuilder_callback
                                        (GtkRange        *range,
                                        gpointer         user_data)
{
  GtkAdjustment *adjustment;
  float fvalue;
  adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
  fvalue = gtk_adjustment_get_value(adjustment);
  preferences_internal_change_value_float(PREFERENCES_MAP_COLOURS_MAP_ROTATION, fvalue);
  set_colour_map_rotation_for_map(fvalue);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SMOOTH_SCROLL, 1);
    set_smooth_scroll_flag(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SMOOTH_SCROLL, 0);
    set_smooth_scroll_flag(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_steps_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  const gchar *text = gtk_entry_get_text(entry);
  int ival = 0;
  ival = atoi(text);
  if ((ival < 10000000) && (ival > 0)) {
    preferences_internal_change_value_int(PREFERENCES_SMOOTH_SCROLL_STEPS, ival);
    set_smooth_scroll_steps(ival);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_steps_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry *entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable),
					    "preferences_smooth_scroll_steps_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  int ival = 0;
  ival = atoi(text);
  if ((ival < 10000000) && (ival > 0)) {
    preferences_internal_change_value_int(PREFERENCES_SMOOTH_SCROLL_STEPS, ival);
    set_smooth_scroll_steps(ival);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_limit_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval < 1000) && (fval > 0)) {
    preferences_internal_change_value_float(PREFERENCES_SMOOTH_SCROLL_LIMIT, fval);
    set_smooth_scroll_limit(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_smooth_scroll_limit_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  GtkEntry *entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable),
					    "preferences_smooth_scroll_limit_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval < 1000) && (fval > 0)) {
    preferences_internal_change_value_float(PREFERENCES_SMOOTH_SCROLL_LIMIT, fval);
    set_smooth_scroll_limit(fval);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_drag_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAP_DRAG, 1);
    set_active_map_drag_flag(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_map_drag_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAP_DRAG, 0);
    set_active_map_drag_flag(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_antialias_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ANTIALIAS, 1);
    set_do_anti_aliasing(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_antialias_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ANTIALIAS, 0);
    set_do_anti_aliasing(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_hid_spherical_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_VT_SURFACE, 2);
    vt_surface(2);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_hid_flat_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_VT_SURFACE, 1);
    vt_surface(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_filechooser_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_CHOOSER, 0);
    set_file_chooser_selector(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_filechooser_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_CHOOSER, 1);
    set_file_chooser_selector(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_file_overwrite_yes_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_OVERWRITE, 1);
    set_file_chooser_overwrite(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_file_overwrite_no_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_OVERWRITE, 0);
    set_file_chooser_overwrite(0);
  }

}

extern "C" G_MODULE_EXPORT
void
on_preferences_file_filter_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_FILTER, 1);
    set_filter_fileselection_filenames(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_file_filter_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_FILTER, 0);
    set_filter_fileselection_filenames(0);
  }

}

extern "C" G_MODULE_EXPORT
void
on_preferences_file_sort_by_date_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_SORT_DATE, 1);
    set_sticky_sort_by_date();
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_file_sort_by_date_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_SORT_DATE, 0);
    unset_sticky_sort_by_date();
  }

}

extern "C" G_MODULE_EXPORT
void
on_preferences_dialog_accept_docked_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *hbox;
  GtkWidget *show_checkbutton;
  GtkWidget *hide_checkbutton;
  hbox             = lookup_widget(GTK_WIDGET(togglebutton), "preferences_dialog_accept_docked_hbox");
  show_checkbutton = lookup_widget(GTK_WIDGET(togglebutton), "preferences_dialog_accept_docked_show_radiobutton");
  hide_checkbutton = lookup_widget(GTK_WIDGET(togglebutton), "preferences_dialog_accept_docked_hide_radiobutton");
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ACCEPT_DIALOG_DOCKED, 1);
    set_accept_reject_dialog_docked(1);
    /* shall update the hbox */
    if (accept_reject_dialog_docked_show_state() == 1) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(show_checkbutton), TRUE);
    } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hide_checkbutton), TRUE);
    }
    gtk_widget_show(hbox);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_dialog_accept_detouched_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *hbox;
  hbox             = lookup_widget(GTK_WIDGET(togglebutton), "preferences_dialog_accept_docked_hbox");
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ACCEPT_DIALOG_DOCKED, 0);
    set_accept_reject_dialog_docked(0);
    /* shall update the hbox */
    if (accept_reject_dialog_docked_show_state() == 1) {
      set_accept_reject_dialog_docked_show(0);
    }
    gtk_widget_hide(hbox);
  }

}

extern "C" G_MODULE_EXPORT
void
on_preferences_dialog_accept_docked_show_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ACCEPT_DIALOG_DOCKED_SHOW, 1);
    set_accept_reject_dialog_docked_show(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_dialog_accept_docked_hide_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ACCEPT_DIALOG_DOCKED_SHOW, 0);
    set_accept_reject_dialog_docked_show(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_dialog_accept_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_IMMEDIATE_REPLACEMENT, 0);
    set_refinement_immediate_replacement(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_dialog_accept_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_IMMEDIATE_REPLACEMENT, 1);
    set_refinement_immediate_replacement(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_recentre_pdb_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_RECENTRE_PDB, 1);
    set_recentre_on_read_pdb(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_recentre_pdb_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_RECENTRE_PDB, 0);
    set_recentre_on_read_pdb(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_console_info_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_CONSOLE_COMMANDS, 1);
    set_console_display_commands_state(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_console_info_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_CONSOLE_COMMANDS, 0);
    set_console_display_commands_state(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_tips_on_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_TIPS, 1);
    set_tip_of_the_day_flag(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_tips_off_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_TIPS, 0);
    set_tip_of_the_day_flag(0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_refinement_speed_molasses_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, 4);
    set_dragged_refinement_steps_per_frame(4);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_refinement_speed_crock_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, 120);
    set_dragged_refinement_steps_per_frame(120);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_refinement_speed_default_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, 80);
    set_dragged_refinement_steps_per_frame(80);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_refinement_speed_own_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  w = lookup_widget(GTK_WIDGET(togglebutton), "preferences_refinement_speed_entry");
  if (gtk_toggle_button_get_active(togglebutton)) {
    const gchar* entry_text = gtk_entry_get_text(GTK_ENTRY(w));
    int val;
    val = atoi(entry_text);
    if ((val > 10000) || (val < 1)) {
      printf("Cannot interpret: %s Assuming default 80 \n", entry_text);
      val  = 80;
      gtk_entry_set_text(GTK_ENTRY(w), "80");
    }
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, val);
    set_dragged_refinement_steps_per_frame(val);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_refinement_speed_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  GtkWidget *w;
  const gchar* entry_text;
  int val;
  w = lookup_widget(GTK_WIDGET(entry), "preferences_refinement_speed_own_radiobutton");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
  entry_text = gtk_entry_get_text(GTK_ENTRY(w));
  val = atoi(entry_text);
  if ((val > 10000) || (val < 1)) {
    printf("Cannot interpret: %s Assuming default 80 \n", entry_text);
    val  = 80;
    gtk_entry_set_text(entry, "80");
  }
  preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, val);
  set_dragged_refinement_steps_per_frame(val);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_refinement_speed_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
					 gpointer         user_data)
{
  GtkWidget *w;
  GtkWidget *togglebutton;
  const gchar* entry_text;
  int val;

  w = lookup_widget(GTK_WIDGET(editable), "preferences_refinement_speed_entry");
  togglebutton = lookup_widget(GTK_WIDGET(editable), "preferences_refinement_speed_own_radiobutton");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(togglebutton), TRUE);

  entry_text = gtk_entry_get_text(GTK_ENTRY(w));
  val = atoi(entry_text);
  if ((val > 10000) || (val < 1)) {
    printf("Cannot interpret: %s Assuming default 80 \n", entry_text);
    val  = 80;
    gtk_entry_set_text(GTK_ENTRY(w), "80");
  }
  preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, val);
  set_dragged_refinement_steps_per_frame(val);

}


extern "C" G_MODULE_EXPORT
void
on_preferences_spin_speed_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  float fval;
  const gchar* entry_text = gtk_entry_get_text(entry);
  fval = atof(entry_text);
  if ((fval > 360) || (fval < 0)) {
    printf("Cannot interpret: %s Assuming default 1.0 \n", entry_text);
    fval  = 1.0;
    gtk_entry_set_text(entry, "1.0");
  }
  preferences_internal_change_value_float(PREFERENCES_SPIN_SPEED, fval);
  set_idle_function_rotate_angle(fval);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_spin_speed_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  float fval;
  GtkWidget *w = lookup_widget(GTK_WIDGET(editable), "preferences_spin_speed_entry");
  const gchar* entry_text = gtk_entry_get_text(GTK_ENTRY(w));
  fval = atof(entry_text);
  if ((fval > 360) || (fval < 0)) {
    printf("Cannot interpret: %s Assuming default 1.0 \n", entry_text);
    fval  = 1.0;
    gtk_entry_set_text(GTK_ENTRY(w), "1.0");
  }
  preferences_internal_change_value_float(PREFERENCES_SPIN_SPEED, fval);
  set_idle_function_rotate_angle(fval);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_font_size_small_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, 1);
    set_font_size(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_font_size_medium_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, 2);
    set_font_size(2);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_font_size_large_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, 3);
    set_font_size(3);
  }

}

extern "C" G_MODULE_EXPORT
void
on_preferences_font_size_others_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
    GtkWidget *w;
    w = lookup_widget(GTK_WIDGET(togglebutton), "preferences_font_size_combobox");
    gint ival = gtk_combo_box_get_active(GTK_COMBO_BOX(w));
    ival += 4;
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, ival);
    set_font_size(ival);
  }
}

extern "C" G_MODULE_EXPORT
void
on_preferences_font_size_combobox_changed_gtkbuilder_callback
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(combobox),
			       "preferences_font_size_others_radiobutton");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
  gint ival = gtk_combo_box_get_active(combobox);
  ival += 4;
  preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, ival);
  set_font_size(ival);

}

extern "C" G_MODULE_EXPORT
void
on_preferences_font_colour_default_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
      preferences_internal_change_value_float3(PREFERENCES_FONT_COLOUR, 1.0, 0.8, 0.8);
      set_font_colour(1.0, 0.8, 0.8);
   }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_font_colour_own_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  GdkColor font_colour;
  float fval1;
  float fval2;
  float fval3;
  int previous_state;

  if (gtk_toggle_button_get_active(togglebutton)) {

    previous_state = preferences_internal_font_own_colour_flag();

    if (previous_state != -1) { 	/* not unset */
      w = lookup_widget(GTK_WIDGET(togglebutton), "preferences_font_colorbutton");
      gtk_color_button_get_color(GTK_COLOR_BUTTON(w), &font_colour);
      fval1 = (float) font_colour.red   / (float) 65535;
      fval2 = (float) font_colour.green / (float) 65535;
      fval3 = (float) font_colour.blue  / (float) 65535;

      preferences_internal_change_value_float3(PREFERENCES_FONT_COLOUR, fval1, fval2, fval3);
      printf("     set_font_colour() - path B\n");
      set_font_colour(fval1, fval2, fval3);
      preferences_internal_change_value_int(PREFERENCES_FONT_OWN_COLOUR_FLAG, 1);
    }
  }
}



extern "C" G_MODULE_EXPORT
void
on_preferences_font_colorbutton_color_set_gtkbuilder_callback
                                        (GtkColorButton  *colorbutton,
                                        gpointer         user_data)
{

  GdkColor font_colour;
  GtkWidget *w;
  float fval1;
  float fval2;
  float fval3;
  gtk_color_button_get_color(colorbutton, &font_colour);
  fval1 = (float) font_colour.red   / (float) 65535;
  fval2 = (float) font_colour.green / (float) 65535;
  fval3 = (float) font_colour.blue  / (float) 65535;

  preferences_internal_change_value_float3(PREFERENCES_FONT_COLOUR, fval1, fval2, fval3);
  set_font_colour(fval1, fval2, fval3);
  /* should set own colour button active (if colours away from default) */
  if (fval1  >= 0.999 &&
      fval2 >= 0.799 && fval2 <= 0.801 &&
      fval3 >= 0.799 && fval3 <= 0.801) {
     // set default button active
     w = lookup_widget(GTK_WIDGET(colorbutton),
                       "preferences_font_colour_default_radiobutton");
     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
  } else {
    /* set own font colour button active */
    w = lookup_widget(GTK_WIDGET(colorbutton), "preferences_font_colour_own_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
  }
}

extern "C" G_MODULE_EXPORT
void
on_preferences_font_colorbutton_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  /* actually not doing anything */
}


extern "C" G_MODULE_EXPORT
void
on_preferences_pink_pointer_entry_activate_gtkbuilder_callback
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
  float fval;
  const gchar* entry_text = gtk_entry_get_text(entry);
  fval = atof(entry_text);
  if ((fval > 1000) || (fval < 0)) {
    printf("Invalid cube size: %s Assuming default 0.1 A \n", entry_text);
    fval  = 0.1;
    gtk_entry_set_text(entry, "0.1");
  }
  preferences_internal_change_value_float(PREFERENCES_PINK_POINTER, fval);
  set_rotation_centre_size(fval);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_pink_pointer_entry_changed_gtkbuilder_callback
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{
  float fval;
  GtkWidget *w = lookup_widget(GTK_WIDGET(editable), "preferences_pink_pointer_entry");
  const gchar* entry_text = gtk_entry_get_text(GTK_ENTRY(w));
  fval = atof(entry_text);
  if ((fval > 1000) || (fval < 0)) {
    printf("Invalid cube size: %s Assuming default 0.1 A \n", entry_text);
    fval  = 0.1;
    gtk_entry_set_text(GTK_ENTRY(w), "0.1");
  }
  preferences_internal_change_value_float(PREFERENCES_PINK_POINTER, fval);
  set_rotation_centre_size(fval);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_bond_width_combobox_changed_gtkbuilder_callback
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
  gint val;
  val = gtk_combo_box_get_active(combobox);
  val += 1;  /* offset */
  preferences_internal_change_value_int(PREFERENCES_BONDS_THICKNESS, val);
  set_default_bond_thickness(val);
}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_show_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    show_modelling_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_SHOW, 1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_right_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(0);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_left_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(1);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_top_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(2);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 2);
  }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_bottom_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(3);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 3);
  }
}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_hide_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    hide_modelling_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_SHOW, 0);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_main_icons_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_all_icons_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_style_icons_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_STYLE, 1);
    set_model_toolbar_style(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_style_both_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_STYLE, 2);
    set_model_toolbar_style(2);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_style_text_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_STYLE, 3);
    set_model_toolbar_style(3);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_show_icon_all_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  show_model_toolbar_all_icons();

}


extern "C" G_MODULE_EXPORT
void
on_preferences_model_toolbar_show_icon_selection_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  show_model_toolbar_main_icons();

}

extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_show_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    show_main_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_SHOW, 1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_hide_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    hide_main_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_SHOW, 0);
  }

}


/* not used currently */
extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_top_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_bottom_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_right_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_left_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_style_icons_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_STYLE, 1);
    set_main_toolbar_style(1);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_style_both_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_STYLE, 2);
    set_main_toolbar_style(2);
  }

}


extern "C" G_MODULE_EXPORT
void
on_preferences_main_toolbar_style_text_radiobutton_toggled_gtkbuilder_callback
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_STYLE, 3);
    set_main_toolbar_style(3);
  }

}

