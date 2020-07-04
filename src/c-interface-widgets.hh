/* src/c-interface-widgets.hh
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2015 by Medical Research Council
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

/* svn $Id: c-interface.h 1458 2007-01-26 20:20:18Z emsley $ */

#ifndef C_INTERFACE_WIDGETS_HH
#define C_INTERFACE_WIDGETS_HH

/*
  The following extern stuff here because we want to return the
  filename from the file entry box.  That code (e.g.) 
  on_ok_button_coordinates_clicked (callback.c), is written and
  compiled in c.
 
  But, we need that function to set the filename in mol_info, which 
  is a c++ class.
 
p  So we need to have this function external for c++ linking.
 
*/

#define COOT_SCHEME_DIR "COOT_SCHEME_DIR"
#define COOT_PYTHON_DIR "COOT_PYTHON_DIR"

/* I think that this should this be a .hh file */

void add_ligand_builder_menu_item_maybe();
void start_ligand_builder_gui_XXX(GtkMenuItem     *menuitem,
				  gpointer         user_data);


/* ------------------------------------------------------------------------- */
/*                    Cif dictionary                                         */
/* ------------------------------------------------------------------------- */
// Add a molecule chooser for a new cif file
//
void add_cif_dictionary_selector_molecule_selector(GtkWidget *fileselection,
						   GtkWidget *aa_hbox);

void add_cif_dictionary_selector_create_molecule_checkbutton(GtkWidget *fileselection,
							     GtkWidget *aa_hbutton_box);
void on_cif_dictionary_file_selector_create_molecule_checkbutton_toggled(GtkButton *button,
									 gpointer user_data);


/* ------------------------------------------------------------------------- */
/*                    Generic Objects                                        */
/* ------------------------------------------------------------------------- */

GtkWidget *wrapped_create_generic_objects_dialog();
/* and this uses callbacks: */
void on_generic_objects_dialog_object_toggle_button_toggled(GtkButton       *button,
							    gpointer         user_data);
/* and... */
void
generic_objects_dialog_table_add_object_internal(const meshed_generic_display_object &gdo,
						 GtkWidget *dialog,
						 GtkWidget *table,
						 int io);

/* return a new object number (so that we can set it to be displayed). */
int add_generic_display_object(const coot::old_generic_display_object_t &gdo);


/*  ----------------------------------------------------------------------- */
/*                  GUIL Utility Functions                                  */
/*  ----------------------------------------------------------------------- */

// gui nuts and bolts
void on_simple_text_dialog_close_button_pressed( GtkWidget *button,
						 GtkWidget *dialog);

#endif // C_INTERFACE_WIDGETS_HH
