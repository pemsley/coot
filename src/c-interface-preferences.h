/* src/atom-selection.hh
 * 
 * Copyright 2011 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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

#ifndef C_INTERFACE_PREFERENCES_H
#define C_INTERFACE_PREFERENCES_H

#include <gtk/gtk.h>

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

/*  ----------------------------------------------------------------------- */
/*                  Preferences Notebook                                    */
/*  ----------------------------------------------------------------------- */
/* section Preferences */
void preferences();
void show_preferences();
void clear_preferences();
void set_mark_cis_peptides_as_bad(int istate); /* in geometry graph */
int show_mark_cis_peptides_as_bad_state();
#if (GTK_MAJOR_VERSION >= 4)
void show_hide_preferences_tabs(GtkToggleButton *toggletoolbutton, int preference_type);
#else
void show_hide_preferences_tabs(GtkToggleToolButton *toggletoolbutton, int preference_type);
#endif
void update_preference_gui();
void make_preferences_internal();
void make_preferences_internal_default();
void reset_preferences();
void save_preferences();
void preferences_internal_change_value_int(int preference_type, int ivalue);
void preferences_internal_change_value_int2(int preference_type, int ivalue1, int ivalue2);
void preferences_internal_change_value_float(int preference_type, float fvalue);
void preferences_internal_change_value_float3(int preference_type, 
					float fvalue1, float fvalue2, float fvalue3);
void show_model_toolbar_icon(int pos);
void hide_model_toolbar_icon(int pos);
void fill_preferences_model_toolbar_icons(GtkWidget *preferences,
				     	  GtkWidget *scrolled_window);
void update_model_toolbar_icons_menu();

void show_main_toolbar_icon(int pos);
void hide_main_toolbar_icon(int pos);
void fill_preferences_main_toolbar_icons(GtkWidget *preferences,
				     	  GtkWidget *scrolled_window);
void update_main_toolbar_icons_menu();
void update_toolbar_icons_menu(int toolbar_index);

int preferences_internal_font_own_colour_flag();

END_C_DECLS

#endif /* C_INTERFACE_PREFERENCES_H */

