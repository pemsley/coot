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

// from support.h
// GtkWidget* lookup_widget (GtkWidget *widget, const gchar *widget_name);
#include "support.h"

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_refine_control_button_clicked_gtkbuilder_callback
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = wrapped_create_refine_params_dialog();
  gtk_widget_show(widget);
}


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

