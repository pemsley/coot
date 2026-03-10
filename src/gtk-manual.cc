/* src/gtk-manual.cc
 *
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Copyright 2008, 2009 by The University of Oxford
 * Author: Paul Emsley
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


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#include <string.h>
#include <functional>

// NOTE:: The order of these 6 include files seems fragile
#include <gtk/gtk.h>
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"

#include "gtk-manual.h"
#include "gtk-manual.hh"

#include "utils/coot-utils.hh"
#include "widget-from-builder.hh"

#include "graphics-info.h" // 20220409-PE now we check use_graphics_interface_flag

#include "utils/logging.hh"
extern logging logger;

/* This is the signal handler for a color change event created when
   the colorseldialog has had its colour changed. */

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                                 map                                      */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void
on_map_color_changed(GtkWidget *w, gpointer tmd) {

   struct map_colour_data_type* t = static_cast<struct map_colour_data_type*> (tmd);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE FIXME color
   GdkRGBA colour;
#else
   GdkColor color;
   gtk_color_selection_get_current_color(t->color_selection, &color);
   GdkRGBA map_color;
   map_color.red   = color.red    /65535.0;;
   map_color.green = color.green  /65535.0;;
   map_color.blue  = color.blue   /65535.0;;
   handle_map_colour_change(t->imol, map_color);
#endif

}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                               symmetry                                   */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

GtkWidget *
create_symmetry_colour_selection_window() {


/*    GtkColorSelectionDialog *colorseldialog;  */
   GtkWidget  *colorseldialog;

   GtkButton *ok_button;
   GtkButton *cancel_button;
   GtkButton *help_button;
   GtkWidget *colorsel;

   // colorseldialog = gtk_color_selection_dialog_new("Symmetry Colour Selection");
   GtkWindow *parent_window = NULL; // FIXME
   colorseldialog = gtk_color_chooser_dialog_new("Symmetry Colour Selection", parent_window);


/* How do we get to the buttons? */

   // colorsel = GTK_COLOR_SELECTION_DIALOG(colorseldialog)->colorsel;

  /* Capture "color_changed" events in col_sel_window */

   /*
  g_signal_connect (G_OBJECT (colorsel), "color_changed",
		    G_CALLBACK(on_symmetry_color_changed),
		    (gpointer)colorsel);

  g_signal_connect(GOBJECT(GTK_COLOR_SELECTION_DIALOG(colorseldialog)->
			   ok_button), "clicked",
		   G_CALLBACK(on_symm_col_sel_cancel_button_clicked),
		   colorseldialog);

  g_signal_connect(G_OBJECT(GTK_COLOR_SELECTION_DIALOG(colorseldialog)->
			    cancel_button), "clicked",
		     G_CALLBACK(on_symm_col_sel_cancel_button_clicked),
		     colorseldialog);
   */

  /* give this widget a name so that we can look it up? */
  g_object_set_data (G_OBJECT (colorseldialog), "symmetry_bonds_colour_selection",
		     colorseldialog);

  return GTK_WIDGET(colorseldialog);

}


#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
void
on_symmetry_color_changed(GtkWidget *w, GtkColorSelection *colorsel) {
   gdouble color[4];

   // gtk_color_selection_get_color(colorsel,color);

   /* we pass back the model number */
   handle_symmetry_colour_change(1,color);
}
#endif

/*  The colour selection dialog has had its OK button pressed */
void
on_symm_col_sel_ok_button_clicked (GtkButton       *button,
				   gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_set_visible(w, FALSE);
}

/*  The colour selection dialog has had its cancel button pressed */
void
on_symm_col_sel_cancel_button_clicked (GtkButton       *button,
				      gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_set_visible(w, FALSE);
				/* we should put the colour back to
				   how it used to be then. */
}

void
create_initial_ramachandran_mol_submenu(GtkWidget *widget) {

#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "in create_initial_ramachandran_mol_submenu() FIXME" << std::endl;
#else
   GtkWidget *window1 = widget;

   GtkWidget *rama_menu = widget_from_builder("ramachandran_plot1");
   // GtkWidget *rama_draw_submenu = gtk_menu_new();
   GMenu *rama_draw_submenu = g_menu_new();
   g_object_set_data_full (G_OBJECT (window1), "rama_plot_submenu", rama_draw_submenu, NULL); // 20211002-PE what does this do now?

   g_menu_item_set_submenu (G_MENU_ITEM (rama_menu), rama_draw_submenu);
   g_object_set_data(G_OBJECT(rama_menu), "rama_plot_submenu", rama_draw_submenu);

#endif
}



void
rama_plot_mol_selector_activate (GMenuItem     *menuitem,
				 gpointer       user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   /* We should come here and be given imol.  New molecules should insert
      themselves into the Ramachandran Plot menu(item). */

   GtkWidget *rama_widget = dynarama_is_displayed_state(imol);
   if (rama_widget == NULL) {
      do_ramachandran_plot(imol);
   } else {
      gtk_widget_set_visible(rama_widget, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME window raise
#else
      gdk_window_raise(GDK_WINDOW(gtk_widget_get_window(rama_widget)));
#endif
   }
}

/* And similar for sequence view: */
void create_initial_sequence_view_mol_submenu(GtkWidget *widget) {

#if (GTK_MAJOR_VERSION >= 4)

   std::cout << "in create_initial_sequence_view_mol_submenu() FIXME" << std::endl;

#else

   GtkWidget *seq_view_draw = widget_from_builder("sequence_view1");
   GtkWidget *seq_view_menu = gtk_menu_new();

   g_object_set_data_full(G_OBJECT(widget), "sequence_view_menu", seq_view_menu, NULL); // 20211002-PE what does this do (likewsie to rama)
   gtk_menu_item_set_submenu (GTK_MENU_ITEM(seq_view_draw), seq_view_menu);
   g_object_set_data(G_OBJECT(seq_view_draw), "sequence_view_submenu", seq_view_menu);
#endif
}

void update_sequence_view_menu_manual(int imol, const char *name) {

   std::cout << "error:: update_sequence_view_menu_manual(): Don't use this " << std::endl;

#if 0 // 20220602-PE and with GTK4 it doesn't compile either.

   // GtkWidget *window1 = lookup_widget(main_window(), "window1");
   GtkWidget *seq_view_menu = widget_from_builder("seq_view_menu");
   GtkWidget *menu_item;

   menu_item = gtk_menu_item_new_with_label (name);

   // 20220309-PE why do I need to do this these days?
   // gtk_widget_ref (menu_item);
   //  g_object_set_data_full (G_OBJECT(window1), "seq_view_menu_item",
   // 			   menu_item,
   // NULL);

   gtk_widget_set_visible(menu_item, TRUE);
   gtk_container_add(GTK_CONTAINER(seq_view_menu), menu_item);
   g_signal_connect (G_OBJECT(menu_item), "activate",
		     G_CALLBACK(sequence_view_mol_selector_activate),
		     GINT_TO_POINTER(imol));
#endif
}

void sequence_view_mol_selector_activate (GMenuItem     *menuitem,
					  gpointer         user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  std::cout << "debug:: sequence_view_mol_selector_activate() calling do_sequence_view() " << imol  << std::endl;
   do_sequence_view(imol);

}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                               skeleton                                   */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

#include "c-interface-widgets.hh"

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
void
on_skeleton_color_changed(GtkWidget *w, GtkColorSelection *colorsel) {
   gdouble color[4];
   for (int i=0; i<4; i++) color[i] = 0.0;

   std::cout << "fix the colour" << std::endl;

   // gtk_color_selection_get_color(colorsel,color);

   /* we pass back the model number */
   handle_skeleton_colour_change(1,color);
}
#endif


/*  The colour selection dialog has had its OK button pressed */
void
on_skeleton_col_sel_ok_button_clicked (GtkButton       *button,
				       gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   // gtk_widget_destroy(w);
   gtk_widget_set_visible(w, FALSE);
}

/*  The colour selection dialog has had its cancel button pressed */
void
on_skeleton_col_sel_cancel_button_clicked (GtkButton       *button,
				      gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   // gtk_widget_destroy(w);
   gtk_widget_set_visible(w, FALSE);
				/* we should put the colour back to
				   how it used to be then. */
}


void
on_display_manager_selections_and_colours_combobox_changed(GtkComboBox     *combo_box,
                                                           gpointer         user_data) {


   // 20220305-PE note to self: opening the Display Manager calls this function for every model

   if (false)
      std::cout << "DEBUG:: :::::: on_display_manager_selections_and_colours_combobox_changed() "
                << combo_box << " " << GPOINTER_TO_INT(user_data) << std::endl;

   int imol = GPOINTER_TO_INT(user_data);
   gchar *txt = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo_box));
   // std::cout << "DEBUG:: text: \"" << txt << "\" user data (imol) " << imol << std::endl;
   logger.log(log_t::DEBUG, "text:", txt, "user data (imol)", imol);

   if (txt) {
      std::string at(txt);
      // std::cout << "at: " << at << std::endl;
      if (at == ("Bonds (Colour by Atom)")) {
         // std::cout << "on_display_manager_selections_and_colours_combobox_changed() display as bonds " <<std::endl;
         graphics_to_bonds_representation(imol);
      }
      if (at == ("C-alphas/Backbone")) {
         // std::cout << "display as CA " <<std::endl;
         graphics_to_ca_representation(imol);
      }
      if (at == ("Bonds (Colour by Chain)")) {
         render_as_bonds_colored_by_chain_button_select(imol);
      }
      if (at == ("Bonds (Colour by Molecule)")) {
         render_as_bonds_colored_by_molecule_button_select(imol);
      }
      // 20211019-PE currently this can't happen
      if (at == ("Bonds (Goodsell Colour by Chain)")) {
         render_as_bonds_goodsell_colored_by_chain_button_select(imol);
      }
      if (at == ("Colour by Sec. Str. Bonds")) {
         render_as_sec_struct_bonds_button_select(imol);
      }
      if (at == ("CAs + Ligands")) {
         render_as_ca_plus_ligands_bonds_button_select(imol);
      }
      if (at ==  ("CAs+Ligs SecStr Col")) {
         render_as_ca_plus_ligands_sec_str_bonds_button_select(imol);
      }
      if (at == ("Jones' Rainbow")) {
         render_as_rainbow_representation_button_select(imol);
      }
      if (at ==  ("Colour by Atom - No Waters")) {
         render_as_bonds_no_waters(imol);
      }
      if (at ==  ("Colour by B-factor - Backbone")) {
         render_as_b_factor_cas_representation_button_select(imol);
      }
      if (at ==  ("Colour by B-factor - All")) {
         render_as_b_factor_representation_button_select(imol);
      }
      if (at ==  ("Colour by Occupancy")) {
         render_as_occupancy_representation_button_select(imol);
      }
   }
}


// make and show (will need to be packed)
//
GtkWidget *selections_and_colours_combobox(int imol) {

   GtkWidget *combobox = gtk_combo_box_text_new();

   gtk_widget_set_margin_start (combobox, 2);
   gtk_widget_set_margin_end   (combobox, 2);
   gtk_widget_set_margin_top   (combobox, 1);
   gtk_widget_set_margin_bottom(combobox, 1);
   
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Bonds (Colour by Atom)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Bonds (Colour by Chain)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Bonds (Colour by Molecule)"));
   // gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Bonds (Goodsell Colour by Chain)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("C-alphas/Backbone"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("CAs + Ligands"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("CAs+Ligs SecStr Col"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Jones' Rainbow"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Colour by Sec. Str. Bonds"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Colour by Atom - No Waters"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Colour by B-factor - Backbone"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Colour by B-factor - All"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), ("Colour by Occupancy"));

   gpointer user_data = GINT_TO_POINTER(imol);
   g_signal_connect(G_OBJECT(combobox), "changed",
                    G_CALLBACK(on_display_manager_selections_and_colours_combobox_changed), user_data);

   int index = 0; // needs to be dynamic (at some stage)

   int bbt = get_graphics_molecule_bond_type(imol);

   // This is ugly.
   // bonds colour by atom: 1
   // bonds colour by molecule: 8
   // bonds colour by chain: 3
   // bonds colour by goodsell: 21
   // bonds colour by Sec Str 6
   // CA: 2
   // CA ligands: 4
   // CA S S: 7
   // Rainbow: 9
   // No waters: 5
   // B-factor backbone: 1... wrong
   // occ: 11
   //

   // enum { UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2,
   //        COLOUR_BY_CHAIN_BONDS=3,
   //        CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
   //        BONDS_NO_HYDROGENS=15,
   //        CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
   //        CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
   //        CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS=17,
   //        COLOUR_BY_MOLECULE_BONDS=8,
   //        COLOUR_BY_RAINBOW_BONDS=9,
   //        COLOUR_BY_B_FACTOR_BONDS=10,
   //        COLOUR_BY_OCCUPANCY_BONDS=11,
   //        COLOUR_BY_USER_DEFINED_COLOURS_BONDS=12 };

   if (bbt ==  3) index =  1;
   if (bbt ==  8) index =  2;
   if (bbt ==  2) index =  3;
   if (bbt ==  4) index =  4;
   if (bbt ==  7) index =  5;
   if (bbt ==  9) index =  6;
   if (bbt ==  6) index =  7;
   if (bbt ==  5) index =  8;
   if (bbt == 14) index =  9;
   if (bbt == 10) index = 10;
   if (bbt == 11) index = 11;

   // std::cout << "Here with bbt " << bbt << " index " << index << std::endl;

   gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), index);
   gtk_widget_set_visible(combobox, TRUE);
   return combobox;

}


/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*              map and molecule display control                            */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */



/* Coordinates  */
/* n is the molecule number (not necessarily the nth element in the
   molecule display VBox) */
void display_control_molecule_combo_box(const std::string &name, int imol,
                                        G_GNUC_UNUSED bool show_add_reps_frame_flag) {

   GtkWidget *display_control_molecule_vbox = widget_from_builder("display_molecule_vbox");

   // Wrap each molecule in a vbox: the hbox for the molecule controls, and below it a mesh toggles section
   GtkWidget *mol_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
   g_object_set_data(G_OBJECT(mol_vbox), "imol", GINT_TO_POINTER(imol));
   gtk_box_append(GTK_BOX(display_control_molecule_vbox), mol_vbox);
   gtk_widget_set_visible(mol_vbox, TRUE);

   GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2); // 2 pixels between widgets inside the box
   gtk_widget_set_margin_start(hbox, 2);
   gtk_widget_set_margin_end(hbox, 8);
   g_object_set_data(G_OBJECT(hbox), "imol", GINT_TO_POINTER(imol)); // so that we can delete this box on delete molecule
   gtk_box_append(GTK_BOX(mol_vbox), hbox);
   gtk_widget_set_visible(hbox, TRUE);

   // Create a vbox for mesh toggle buttons (initially hidden, shown when meshes are added)
   GtkWidget *mesh_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
   gtk_widget_set_margin_start(mesh_vbox, 30);
   g_object_set_data(G_OBJECT(mol_vbox), "mesh_vbox", mesh_vbox);
   gtk_box_append(GTK_BOX(mol_vbox), mesh_vbox);
   gtk_widget_set_visible(mesh_vbox, FALSE);

   // We need to add thesee items:
   // 1: molecule number label
   // 2: entry with molecule name
   // 3: "Display" checkbutton
   // 4: "Active" checkbutton
   // 5: "DrawingMode" combobox
   // 6: "Delete model" button -  done in display_control_add_delete_molecule_button()

   // 1: molecule number label
   std::string imol_str = std::to_string(imol);
   GtkWidget *molecule_number_label = gtk_label_new(imol_str.c_str());
   gtk_widget_set_size_request(molecule_number_label, 20, -1); // not using margins
   gtk_widget_set_visible(molecule_number_label, TRUE);
   gtk_box_append(GTK_BOX(hbox), molecule_number_label);

   // 2: entry with molecule name
   GtkWidget *entry = gtk_entry_new();
   gtk_widget_set_hexpand(entry, TRUE);
   gtk_editable_set_text(GTK_EDITABLE(entry), name.c_str());
   gtk_widget_set_visible(entry, TRUE);
   gtk_box_append(GTK_BOX(hbox), entry);
   std::string widget_name = "display_model_entry_" + std::to_string(imol);
   gtk_widget_set_name(entry, widget_name.c_str());

   // 3: "Display" checkbutton
   GtkWidget *display_checkbutton = gtk_check_button_new_with_label("Display");
   gtk_widget_set_visible(display_checkbutton, TRUE);
   g_object_set_data(G_OBJECT(display_checkbutton), "imol", GINT_TO_POINTER(imol));
   gtk_box_append(GTK_BOX(hbox), display_checkbutton);
   gtk_check_button_set_active(GTK_CHECK_BUTTON(display_checkbutton), mol_is_displayed(imol));

   // 4: "Active" checkbutton
   GtkWidget *active_checkbutton = gtk_check_button_new_with_label("Active");
   gtk_widget_set_visible(active_checkbutton, TRUE);
   g_object_set_data(G_OBJECT(active_checkbutton), "imol", GINT_TO_POINTER(imol));
   gtk_box_append(GTK_BOX(hbox), active_checkbutton);
   gtk_check_button_set_active(GTK_CHECK_BUTTON(active_checkbutton), mol_is_active(imol));
   // when Display is untoggled we need to untoggle this active_checkbutton too
   g_object_set_data(G_OBJECT(display_checkbutton), "active_check_button", active_checkbutton);

   // 5: Drawing mode
   GtkWidget *sel_and_col_combobox = selections_and_colours_combobox(imol);
   gtk_box_append(GTK_BOX(hbox), sel_and_col_combobox);

   // when Display is untoggled via the API, we need to get to the buttons to change
   // the state in the gui (give the hbox) (see set_display_control_button_state()).
   g_object_set_data(G_OBJECT(hbox), "display_check_button", display_checkbutton);
   g_object_set_data(G_OBJECT(hbox),  "active_check_button",  active_checkbutton);

   // 6: "Delete" map button
   display_control_add_delete_molecule_button(imol, hbox, display_control_molecule_vbox, false);

   // connect signals to display and active
   g_signal_connect(G_OBJECT(display_checkbutton), "toggled",
                    G_CALLBACK(on_display_control_mol_displayed_button_toggled),
                    GINT_TO_POINTER(imol));
   g_signal_connect(G_OBJECT(active_checkbutton), "toggled",
                    G_CALLBACK(on_display_control_mol_active_button_toggled),
                    GINT_TO_POINTER(imol));

}

void display_control_add_delete_molecule_button(int imol,
                                                GtkWidget *hbox32,
                                                GtkWidget *vbox_for_molecules,
						bool is_map_molecule) {

   if (! hbox32) {
      std::cout << "ERROR:: in display_control_add_delete_molecule_button() null hbox32" << std::endl;
      return;
   }

   std::string button_string = "Delete Model";
   if (is_map_molecule)
      button_string = "Delete Map";
   GtkWidget *delete_button = gtk_button_new_with_label((button_string.c_str()));
   gtk_widget_set_visible(delete_button, TRUE);

   // used in callback on_display_control_mol_active_button_toggled()
   // GtkWidget *vbox_for_molecules     = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "vbox_for_molecules"));
   // GtkWidget *vbox_for_this_molecule = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "vbox_for_this_molecule"));

   g_object_set_data(G_OBJECT(delete_button), "hbox_for_this_molecule", hbox32); // bad name
   g_object_set_data(G_OBJECT(delete_button), "vbox_for_molecules",     vbox_for_molecules);

   gtk_box_append(GTK_BOX (hbox32), delete_button);
   gtk_widget_set_margin_start (delete_button, 2);
   gtk_widget_set_margin_end   (delete_button, 2);
   gtk_widget_set_margin_top   (delete_button, 1);
   gtk_widget_set_margin_bottom(delete_button, 1);
   g_signal_connect(G_OBJECT(delete_button), "clicked",
		      G_CALLBACK(on_display_control_delete_molecule_button_clicked),
		      GINT_TO_POINTER(imol));

}

// Create a "All On" check button that acts on all add reps of this molecule.
// Underneath that, create a frame and in it put a vbox.
//
void add_add_reps_frame_and_vbox(GtkWidget *display_control_window_glade,
				 GtkWidget *hbox_for_single_molecule, int imol_no,
				 bool show_add_reps_frame_flag) {

   GtkWidget *frame = gtk_frame_new("Additional Representations");
   GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
   // show the frame if there were additional representations
   if (show_add_reps_frame_flag)
      gtk_widget_set_visible(frame, TRUE);
   // gtk_widget_ref (v);
   // gtk_widget_ref (frame);

   std::string arl = "   Show Additional Representations  ";
   GtkWidget *all_on_check_button = gtk_check_button_new_with_label(arl.c_str());
   if (show_add_reps_frame_flag)
      gtk_widget_set_visible(all_on_check_button, TRUE);
   gtk_box_append(GTK_BOX(hbox_for_single_molecule), all_on_check_button);
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(all_on_check_button), TRUE);
   std::string widget_name = "add_rep_all_on_check_button_";
   widget_name += coot::util::int_to_string(imol_no);
   g_object_set_data_full (G_OBJECT (display_control_window_glade),
                           widget_name.c_str(),
                           all_on_check_button, NULL);
   g_signal_connect(G_OBJECT (all_on_check_button),  "toggled",
                    G_CALLBACK (on_add_rep_all_on_check_button_toggled),
                    GINT_TO_POINTER(imol_no));

   // set the name so that it can be looked up.
   widget_name = "add_rep_display_control_frame_vbox_";
   widget_name += coot::util::int_to_string(imol_no);

   // g_object_set_data_full(G_OBJECT (display_control_window_glade), widget_name.c_str(), v, NULL);

   widget_name = "add_rep_display_control_frame_";
   widget_name += coot::util::int_to_string(imol_no);
   g_object_set_data_full (G_OBJECT (display_control_window_glade), widget_name.c_str(), frame, NULL);

   gtk_box_append(GTK_BOX(hbox_for_single_molecule), frame);
   gtk_frame_set_child(GTK_FRAME(frame), vbox);

}


void
on_add_rep_all_on_check_button_toggled(GtkToggleButton    *button,
                                       gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   if (gtk_toggle_button_get_active(button)) {
      // std::cout << "Turn on all add reps of " << imol << std::endl;
      set_show_all_additional_representations(imol, 1);
   } else {
      // std::cout << "Turn off all add reps of " << imol << std::endl;
      set_show_all_additional_representations(imol, 0);
   }
}

static void
on_display_control_mesh_toggle_toggled(GtkCheckButton *button, gpointer user_data) {

   int imol     = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "imol"));
   int mesh_idx = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "mesh_idx"));
   bool state = gtk_check_button_get_active(button);
   if (is_valid_model_molecule(imol)) {
      auto &m = graphics_info_t::molecules[imol];
      if (mesh_idx >= 0 && mesh_idx < static_cast<int>(m.meshes.size())) {
         m.meshes[mesh_idx].set_draw_mesh_state(state);
      }
   }
   graphics_draw();
}

static void
on_display_control_generic_object_toggle_toggled(GtkCheckButton *button, gpointer user_data) {

   int obj_idx = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "generic_object_idx"));
   bool state = gtk_check_button_get_active(button);
   auto &gdo = graphics_info_t::generic_display_objects;
   if (obj_idx >= 0 && obj_idx < static_cast<int>(gdo.size())) {
      gdo[obj_idx].mesh.set_draw_mesh_state(state);
   }
   graphics_draw();
}

// Find the mol_vbox for a given molecule in the display control window
static GtkWidget *find_display_control_mol_vbox(int imol) {

   GtkWidget *display_control_molecule_vbox = widget_from_builder("display_molecule_vbox");
   if (!display_control_molecule_vbox) return nullptr;
   GtkWidget *child = gtk_widget_get_first_child(display_control_molecule_vbox);
   while (child) {
      int imol_child = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(child), "imol"));
      if (imol_child == imol)
         return child;
      child = gtk_widget_get_next_sibling(child);
   }
   return nullptr;
}

void update_display_control_mesh_toggles(int imol) {

   if (! graphics_info_t::use_graphics_interface_flag) return;

   GtkWidget *mol_vbox = find_display_control_mol_vbox(imol);
   if (!mol_vbox) return;

   GtkWidget *mesh_vbox = GTK_WIDGET(g_object_get_data(G_OBJECT(mol_vbox), "mesh_vbox"));
   if (!mesh_vbox) return;

   if (!is_valid_model_molecule(imol)) return;
   auto &mol = graphics_info_t::molecules[imol];
   auto &gdo = graphics_info_t::generic_display_objects;

   // Count existing toggle buttons in mesh_vbox
   int n_existing = 0;
   for (GtkWidget *w = gtk_widget_get_first_child(mesh_vbox); w; w = gtk_widget_get_next_sibling(w))
      n_existing++;

   // Count how many toggles we need: molecule meshes + generic display objects for this imol
   int n_mol_meshes = mol.meshes.size();
   int n_gdo_for_imol = 0;
   for (unsigned int i = 0; i < gdo.size(); i++)
      if (gdo[i].imol == imol && !gdo[i].mesh.name.empty())
         n_gdo_for_imol++;

   int n_total = n_mol_meshes + n_gdo_for_imol;
   if (n_total == 0) {
      gtk_widget_set_visible(mesh_vbox, FALSE);
      return;
   }

   // Only add new toggles beyond what already exists.
   // Molecule meshes come first, then generic display objects.
   int toggle_idx = 0;

   // Add toggle buttons for molecule meshes
   for (unsigned int j = 0; j < mol.meshes.size(); j++) {
      if (toggle_idx >= n_existing) {
         const auto &mesh = mol.meshes[j];
         std::string label = mesh.name;
         GtkWidget *cb = gtk_check_button_new_with_label(label.c_str());
         gtk_check_button_set_active(GTK_CHECK_BUTTON(cb), mesh.get_draw_this_mesh());
         g_object_set_data(G_OBJECT(cb), "imol",     GINT_TO_POINTER(imol));
         g_object_set_data(G_OBJECT(cb), "mesh_idx", GINT_TO_POINTER(j));
         g_signal_connect(G_OBJECT(cb), "toggled",
                          G_CALLBACK(on_display_control_mesh_toggle_toggled), nullptr);
         gtk_box_append(GTK_BOX(mesh_vbox), cb);
         gtk_widget_set_visible(cb, TRUE);
         gtk_widget_set_margin_start(cb, 4);
      }
      toggle_idx++;
   }

   // Add toggle buttons for generic display objects belonging to this molecule
   for (unsigned int i = 0; i < gdo.size(); i++) {
      if (gdo[i].imol == imol && !gdo[i].mesh.name.empty()) {
         if (toggle_idx >= n_existing) {
            std::string label = gdo[i].mesh.name;
            GtkWidget *cb = gtk_check_button_new_with_label(label.c_str());
            gtk_check_button_set_active(GTK_CHECK_BUTTON(cb), gdo[i].mesh.get_draw_this_mesh());
            g_object_set_data(G_OBJECT(cb), "imol",               GINT_TO_POINTER(imol));
            g_object_set_data(G_OBJECT(cb), "generic_object_idx", GINT_TO_POINTER(i));
            g_signal_connect(G_OBJECT(cb), "toggled",
                             G_CALLBACK(on_display_control_generic_object_toggle_toggled), nullptr);
            gtk_box_append(GTK_BOX(mesh_vbox), cb);
            gtk_widget_set_visible(cb, TRUE);
            gtk_widget_set_margin_start(cb, 4);
         }
         toggle_idx++;
      }
   }
   gtk_widget_set_visible(mesh_vbox, TRUE);
}

GtkWidget *molecule_index_to_display_manager_entry(int imol) {

   // see also set_display_control_button_state()

   // Search for a widget by name up to 3 levels deep in the hierarchy
   std::function<GtkWidget *(const std::string &, GtkWidget *, int)> find_widget_by_name;
   find_widget_by_name = [&find_widget_by_name] (const std::string &widget_name, GtkWidget *parent, int depth) -> GtkWidget * {
      if (depth <= 0) return nullptr;
      GtkWidget *child = gtk_widget_get_first_child(parent);
      while (child) {
         const char *wn = gtk_widget_get_name(child);
         if (wn) {
            std::string wns(wn);
            if (widget_name == wns) return child;
         }
         GtkWidget *found = find_widget_by_name(widget_name, child, depth - 1);
         if (found) return found;
         child = gtk_widget_get_next_sibling(child);
      }
      return nullptr;
   };
   auto get_inner_entry = [&find_widget_by_name] (const std::string &widget_name, GtkWidget *vbox) {
      return find_widget_by_name(widget_name, vbox, 3);
   };

   GtkWidget *entry = nullptr;
   std::string idx_part = std::to_string(imol);

   if (is_valid_map_molecule(imol)) {
      std::string widget_name = "display_map_entry_" + idx_part;
      GtkWidget *vbox = widget_from_builder("display_map_vbox");
      entry = get_inner_entry(widget_name, vbox);
   }

   if (is_valid_model_molecule(imol)) {
      std::string widget_name = "display_model_entry_" + idx_part;
      GtkWidget *vbox = widget_from_builder("display_molecule_vbox");
      entry = get_inner_entry(widget_name, vbox);
   }

   return entry;
}

// we are using const char * here because update_name_in_display_control_molecule_combo_box is
// declared in gtk-manual.h. Sort that out another day.
void
update_name_in_display_control_molecule_combo_box(int imol, const char *display_name) {

   GtkWidget *entry = molecule_index_to_display_manager_entry(imol);
   if (entry)
      gtk_editable_set_text(GTK_EDITABLE(entry), display_name);
   else
      std::cout << "DEBUG:: update_name_in_display_control_molecule_combo_box() entry lookup failed" << std::endl;
}

void
render_as_bonds_button_select(int imol) {

   graphics_to_bonds_representation(imol);
}

void
render_as_bonds_colored_by_chain_button_select(int imol) {

  set_colour_by_chain(imol);
}

void
render_as_bonds_goodsell_colored_by_chain_button_select(int imol) {

  set_colour_by_chain_goodsell_mode(imol);
}


void
render_as_bonds_colored_by_molecule_button_select(int imol) {

  set_colour_by_molecule(imol);
}

void
render_as_bonds_no_waters(int imol) {

  graphics_to_bonds_no_waters_representation(imol);
}


void
render_as_ca_bonds_button_select(int imol) {

   graphics_to_ca_representation(imol);
}

void
render_as_ca_plus_ligands_bonds_button_select(int imol) {
   graphics_to_ca_plus_ligands_representation(imol);
}

void
render_as_ca_plus_ligands_sec_str_bonds_button_select(int imol) {
   graphics_to_ca_plus_ligands_sec_struct_representation(imol);
}

void
render_as_sec_struct_bonds_button_select(int imol) {
   graphics_to_sec_struct_bonds_representation(imol);
}

void render_as_rainbow_representation_button_select(int imol) {
   graphics_to_rainbow_representation(imol);
}

void render_as_b_factor_representation_button_select(int imol) {
   graphics_to_b_factor_representation(imol);

}

void render_as_b_factor_cas_representation_button_select(int imol) {
   graphics_to_b_factor_cas_representation(imol);

}

void render_as_occupancy_representation_button_select(int imol) {
   graphics_to_occupancy_representation(imol);
}



/* n is the nth element in the molecule display VBox */
void
display_control_map_combo_box(const std::string &name, int imol) {

   if (! graphics_info_t::use_graphics_interface_flag) return;

   GtkWidget *display_map_vbox = widget_from_builder("display_map_vbox");
   GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2); // 2 pixels between widgets
   gtk_widget_set_margin_start(hbox, 2);
   gtk_widget_set_margin_end(hbox, 8);
   g_object_set_data(G_OBJECT(hbox), "imol", GINT_TO_POINTER(imol)); // so that we can delete this box on delete molecule
   gtk_box_append(GTK_BOX(display_map_vbox), hbox);
   gtk_widget_set_visible(hbox, TRUE);

   // Do I need to clear_out_container(display_map_vbox) here?

   // We need to add thesee items:
   // 1: molecule number label
   // 2: entry with molecule name
   // 3: "Display" checkbutton
   // 4: "Scroll" radiobutton
   // 5: "Properties" button
   // 6: "Delete Map" button -  done in display_control_add_delete_molecule_button()

   // 1: molecule number label
   std::string imol_str = std::to_string(imol);
   GtkWidget *molecule_number_label = gtk_label_new(imol_str.c_str());
   gtk_widget_set_visible(molecule_number_label, TRUE);
   gtk_widget_set_size_request(molecule_number_label, 20, -1); // not using margins (testing)
   gtk_box_append(GTK_BOX(hbox), molecule_number_label);

   // 2: entry with molecule name
   GtkWidget *entry = gtk_entry_new();
   gtk_widget_set_hexpand(entry, TRUE);
   gtk_editable_set_text(GTK_EDITABLE(entry), name.c_str());
   gtk_widget_set_visible(entry, TRUE);
   gtk_box_append(GTK_BOX(hbox), entry);
   std::string widget_name = "display_map_entry_" + std::to_string(imol);
   gtk_widget_set_name(entry, widget_name.c_str());

   // 3: "Display" checkbutton
   GtkWidget *display_checkbutton = gtk_check_button_new_with_label("Display");
   gtk_widget_set_margin_start (display_checkbutton, 2);
   gtk_widget_set_margin_end   (display_checkbutton, 2);
   gtk_widget_set_margin_top   (display_checkbutton, 1);
   gtk_widget_set_margin_bottom(display_checkbutton, 1);
   gtk_widget_set_visible(display_checkbutton, TRUE);
   gtk_box_append(GTK_BOX(hbox), display_checkbutton);
   gtk_check_button_set_active(GTK_CHECK_BUTTON(display_checkbutton), map_is_displayed(imol));
   g_object_set_data(G_OBJECT(hbox), "display_check_button", display_checkbutton); // for set_display_control_button_state()

   // 4: "Scroll" checkbutton
   GtkWidget *scroll_button = gtk_check_button_new_with_label("Scroll");
   GtkWidget *scroll_group  = get_radio_button_in_scroll_group(imol);
   if (scroll_group)
      gtk_check_button_set_group(GTK_CHECK_BUTTON(scroll_button), GTK_CHECK_BUTTON(scroll_group));
   g_object_set_data(G_OBJECT(scroll_button), "imol", GINT_TO_POINTER(imol));

   // std::cout << ":::::::::::::::: scroll wheel map " << graphics_info_t::scroll_wheel_map << std::endl;
   // maybe scroll_wheel_map was not set yet? So set it now, for this map
   if (graphics_info_t::scroll_wheel_map == -1)
      graphics_info_t::scroll_wheel_map = imol;
   if (imol == graphics_info_t::scroll_wheel_map)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(scroll_button), TRUE);
   gtk_box_append(GTK_BOX(hbox), scroll_button);

   // 5: "Properties" button
   GtkWidget *properties_button = gtk_button_new_with_label("Properties");
   gtk_widget_set_margin_start (properties_button, 2);
   gtk_widget_set_margin_end   (properties_button, 2);
   gtk_widget_set_margin_top   (properties_button, 1);
   gtk_widget_set_margin_bottom(properties_button, 1);
   gtk_widget_set_visible(properties_button, TRUE);
   gtk_box_append(GTK_BOX(hbox), properties_button);

   // 6: "Delete" map button
   display_control_add_delete_molecule_button(imol, hbox, display_map_vbox, true);

   // connect signals to display, scroll and properties
   g_signal_connect(G_OBJECT(display_checkbutton), "toggled",
                    G_CALLBACK(on_display_control_map_displayed_button_toggled),
                    GINT_TO_POINTER(imol));
   g_signal_connect(G_OBJECT(scroll_button), "toggled",
                    G_CALLBACK(on_display_control_map_scroll_radio_button_toggled),
                    GINT_TO_POINTER(imol));
   g_signal_connect(G_OBJECT(properties_button), "clicked",
		     G_CALLBACK (on_display_control_map_properties_button_clicked),
		     GINT_TO_POINTER(imol));

}


void
on_display_control_map_displayed_button_toggled(GtkCheckButton       *button,
                                                gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   if (gtk_check_button_get_active(button)) {
      set_map_displayed(imol, 1);
   } else {
      set_map_displayed(imol, 0);
   }
}

/* Added 20050316 (Bangalore) */
void
on_display_control_map_scroll_radio_button_toggled (GtkCheckButton *button,
						    gpointer         user_data) {
   int imol = GPOINTER_TO_INT(user_data);
   if (gtk_check_button_get_active(button)) {
      set_scrollable_map(imol);
   }
}

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME radio buttons
#else
void
on_display_control_map_scroll_radio_button_group_changed (GtkRadioButton *button,
							  gpointer         user_data) {
   // do nothing these days. Scroll group is handled gtk-internally via radio buttons.
}
#endif

GtkWidget *get_radio_button_in_scroll_group(int imol_this) {

   GtkWidget *w = nullptr;
   GtkWidget *display_map_vbox = widget_from_builder("display_map_vbox");

   // The children of the display_map_vbox are hboxes
   GtkWidget *hbox_child = gtk_widget_get_first_child(display_map_vbox);
   if (hbox_child) {
      GtkWidget *item_widget = gtk_widget_get_first_child(hbox_child);
      if (item_widget) {
         unsigned int inner_child_count = 1;
         while (item_widget && !w) {
            item_widget = gtk_widget_get_next_sibling(item_widget);
            inner_child_count++;
            if (inner_child_count == 4) {
               if (GTK_IS_CHECK_BUTTON(item_widget)) {
                  w = item_widget;
               }
            }
         }
      }
   }

   return w;
}

// return NULL when there are no radio buttons found that have a group
// (i.e. there are no maps displayed in the Display Manager.
//
GtkWidget *get_radio_button_in_scroll_group_old(int imol_this) {

   // 20090517: Previously, we'd been looking at all
   // map_scroll_buttons (i.e. for every map), but in the case of 2
   // maps (say 1 and 3) when we are rendering the first combo_box the
   // combo_box for map 3 doesn't exist - so we get a widget not found
   // warning (ugly).
   //
   // This is only called from one place, (in the combo_box renderer)
   // so let's only look for map molecule with a molecule number less
   // than imol_this.  (if not found, return null of course).

   //
   //
   GtkWidget *w = NULL;
   if (true) {
      // for (int i=0; i<graphics_n_molecules(); i++) {
      //
      // 20120124 but surely when we are creating a new dialog and we
      // are here from the first map, we don't want to check all the
      // next maps // for previous scroll buttons!?  - Those widgets
      // haven't been created yet.
      for (int i=0; i<imol_this; i++) { // 20120124
	 if (is_valid_map_molecule(i)) {
	    if (i != imol_this) {
	       std::string test_name = "map_scroll_button_";
	       test_name += coot::util::int_to_string(i);
	       // w = lookup_widget(display_manager_dialog, test_name.c_str());
               std::cout << "get_radio_button_in_scroll_group(): do a proper lookup of w here "  << std::endl;
	       w = 0;
	       if (w)
		  break;
	    }
	 }
      }
   }
   return w;
}


void
on_display_control_delete_molecule_button_clicked(GtkButton       *button,
                                                  gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   if (true)
      std::cout << "DEBUG:: calling close_molecule() for " << imol << " from "
		<< "on_display_control_delete_molecule_button_clicked"
		<< std::endl;

   // these are set in ...
   GtkWidget *vbox_for_molecules     = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "vbox_for_molecules"));
   GtkWidget *hbox_for_this_molecule = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "hbox_for_this_molecule"));

   // std::cout << "here are the widgets! " << vbox_for_molecules << " " << hbox_for_this_molecule << std::endl;
   if (vbox_for_molecules) {
      gtk_box_remove(GTK_BOX(vbox_for_molecules), GTK_WIDGET(hbox_for_this_molecule));
   }

   close_molecule(imol);
}

#include "single-map-properties-dialog.hh"

void set_transient_for_main_window(GtkWidget *dialog);

void
on_display_control_map_properties_button_clicked(GtkButton       *button,
                                                 gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   // GtkWidget *window = create_single_mmdb:: Atom *atp_properties_dialog();
   // GtkWidget *patch_frame = lookup_widget(window, "single_map_colour_button_frame");

  // GtkWidget *single_map_properties_colour_button = lookup_widget(window, "single_map_properties_colour_button");
  // GtkWidget *label = lookup_widget(window, "label114");

  // FIXME? block commented.  Did it do anything?
//   // chunk from the glade FAQ:
//   GdkColor red = { 0, 65535, 0, 0 };
//   GtkRcStyle *rc_style = gtk_rc_style_new ();
//   rc_style->bg[GTK_STATE_NORMAL] = red;
//   rc_style->color_flags[GTK_STATE_NORMAL] |= GTK_RC_BG;
//   // gtk_widget_modify_style (single_map_properties_colour_button, rc_style);
//   // can't set patch frame fg usefully.
//   gtk_widget_modify_style (label, rc_style);
//   gtk_rc_style_unref (rc_style);


   GtkWidget *dialog = wrapped_create_single_map_properties_dialog_gtk3(imol);
   if (dialog) {
      set_transient_for_main_window(dialog);
      gtk_widget_set_visible(dialog, TRUE);
   }

}

void
fill_map_colour_patch(GtkWidget *patch_frame, int imol){

  GdkRGBA color;
  int width, height;
  GtkWidget *widget;
  GdkRGBA mol_colour = get_map_colour(imol);
  gushort red, green, blue;
  GtkWidget *widget_thing;

  red   =  gushort (mol_colour.red   * 254.0);
  green =  gushort (mol_colour.green * 254.0);
  blue  =  gushort (mol_colour.blue  * 254.0);

  widget = patch_frame;

  widget = gtk_drawing_area_new();
  //  widget_thing = lookup_widget(GTK_WIDGET(patch_frame), "single_map_colour_hbox");
  widget_thing = widget_from_builder("single_map_colour_hbox");
  widget_thing = gtk_window_new();


  printf("adding widget to patch_frame\n");
  gtk_window_set_child(GTK_WINDOW(widget_thing), widget);

  // printf("gdk_gc_new\n");
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
  GdkWindow *window = 0;
  window = gtk_widget_get_window(widget);
  if (! window)
     return;
#endif


  // gc = gdk_gc_new(window);

  printf("get window size\n");
  /* find proper dimensions for rectangle */
  // gtk_window_get_size(window, &width, &height);

  /* the color we want to use */
  color = mol_colour;

  /* red, green, and blue are passed values, indicating the RGB triple
   * of the color we want to draw. Note that the values of the RGB components
   * within the GdkColor are taken from 0 to 65535, not 0 to 255.
   */
  color.red = red * (65535/255);
  color.green = green * (65535/255);
  color.blue = blue * (65535/255);


  /* the pixel value indicates the index in the colormap of the color.
   * it is simply a combination of the RGB values we set earlier
   */
  // color.pixel = (gulong)(red*65536 + green*256 + blue);

  /* However, the pixel valule is only truly valid on 24-bit (TrueColor)
   * displays. Therefore, this call is required so that GDK and X can
   * give us the closest color available in the colormap
   */
  printf("colour alloc\n");

  // gdk_color_alloc(gtk_widget_get_colormap(widget), color);

  /* set the foreground to our color */
  printf("set background\n");

  // gdk_gc_set_background(gc, color);

  /* draw the rectangle */
  printf("draw rectangle:\n");

  // gdk_draw_rectangle(window, gc, 1, 0, 0, width, height);
}



void
on_display_control_mol_displayed_button_toggled(GtkCheckButton *check_button,
                                                gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   // 20220807-PE what does this do!?
   GtkWidget *active_check_button = GTK_WIDGET(g_object_get_data(G_OBJECT(check_button), "active_check_button"));

   if (gtk_check_button_get_active(check_button)) {
      gtk_check_button_set_active(GTK_CHECK_BUTTON(active_check_button), TRUE);
      set_mol_displayed(imol, 1);
      // set_mol_active(imol, 1);
   } else {
      set_mol_displayed(imol, 0);
      gtk_check_button_set_active(GTK_CHECK_BUTTON(active_check_button), FALSE);
      // set_mol_active(imol, 0);
   }

}


void
on_display_control_mol_active_button_toggled(GtkCheckButton *check_button,
                                             gpointer        user_data) {

   // why not pass imol as user data?
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(check_button), "imol"));
   if (gtk_check_button_get_active(check_button)) {
      set_mol_active(imol, 1);
   } else {
      set_mol_active(imol, 0);
   }
}


/* ------------------------------------------------------------------------------ */
/*              cell chooser for phs                                              */
/* ------------------------------------------------------------------------------ */

GSList *display_cell_chooser_box(GtkWidget *phs_cell_choice_window,
				 GSList *phs_cell_group,
				 int n) {

  GtkWidget *vbox39;
  GtkWidget *phs_cell_chooser_vbox;
  GtkWidget *hbox33;
  GtkWidget *phs_cell_radiobutton_1 = nullptr;
  GtkWidget *label53;
  GtkWidget *phs_cell_symm_entry_1;
  GtkWidget *label54;
  GtkWidget *phs_cell_a_entry_1;
  GtkWidget *label55;
  GtkWidget *phs_cell_b_entry_1;
  GtkWidget *label56;
  GtkWidget *phs_cell_c_entry_1;
  GtkWidget *label57;
  GtkWidget *phs_cell_alpha_entry_1;
  GtkWidget *label58;
  GtkWidget *phs_cell_beta_entry_1;
  GtkWidget *label59;
  GtkWidget *phs_cell_gamma_entry_1;
  GtkWidget *hbox34;
  GtkWidget *phs_cell_none_radiobutton = nullptr;
  GtkWidget *hbox32;
  GtkWidget *phs_cell_choice_ok_button;
  GtkWidget *phs_cell_choice_cancel_button;

/* messing about with string variables */
  gchar *widget_name;
  gchar *tmp_name;

  widget_name = (gchar *) malloc(100);

  // phs_cell_chooser_vbox = lookup_widget(phs_cell_choice_window, "phs_cell_chooser_vbox");
  phs_cell_chooser_vbox = widget_from_builder("phs_cell_chooser_vbox");

  hbox33 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox33);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "hbox33", hbox33,
			  NULL);
  gtk_widget_set_visible (hbox33, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX(phs_cell_chooser_vbox), hbox33);
#else
  gtk_box_pack_start (GTK_BOX (phs_cell_chooser_vbox), hbox33, TRUE, TRUE, 4);
  gtk_container_set_border_width (GTK_CONTAINER (hbox33), 6);
#endif

  strcpy(widget_name, "phs_cell_radiobutton_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


#if (GTK_MAJOR_VERSION >= 4)
  // 20220602-PE FIXME radio buttons
#else
  phs_cell_radiobutton_1 = gtk_radio_button_new_with_label (phs_cell_group, "");
  phs_cell_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (phs_cell_radiobutton_1));
#endif
  // gtk_widget_ref (phs_cell_radiobutton_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_radiobutton_1,
			  NULL);
  gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (phs_cell_radiobutton_1), FALSE);

  gtk_widget_set_visible (phs_cell_radiobutton_1, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX(hbox33), phs_cell_radiobutton_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_radiobutton_1, FALSE, FALSE, 4);
#endif

  label53 = gtk_label_new (("Symm"));
  // gtk_widget_ref (label53);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label53", label53, NULL);
  gtk_widget_set_visible (label53, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX(hbox33), label53);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label53, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_symm_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_symm_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_symm_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_symm_entry_1,
			  NULL);
  gtk_widget_set_visible (phs_cell_symm_entry_1, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), phs_cell_symm_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_symm_entry_1, TRUE, TRUE, 0);
#endif

  gtk_widget_set_size_request(GTK_WIDGET(phs_cell_symm_entry_1), 80, -2);

  label54 = gtk_label_new (("a"));
  // gtk_widget_ref (label54);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label54", label54, NULL);
  gtk_widget_set_visible (label54, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), label54);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label54, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_a_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_a_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_a_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_a_entry_1, NULL);

  gtk_widget_set_visible (phs_cell_a_entry_1, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), phs_cell_a_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_a_entry_1, TRUE, TRUE, 0);
#endif
  gtk_widget_set_size_request (phs_cell_a_entry_1, 65, -2);

  label55 = gtk_label_new (("b"));
  // gtk_widget_ref (label55);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label55", label55,
			  NULL);
  gtk_widget_set_visible (label55, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), label55);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label55, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_b_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_b_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_b_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_b_entry_1, NULL);
  gtk_widget_set_visible (phs_cell_b_entry_1, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX(hbox33), phs_cell_b_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_b_entry_1, TRUE, TRUE, 0);
#endif
  gtk_widget_set_size_request (phs_cell_b_entry_1, 65, -2);

  label56 = gtk_label_new (("c"));
  // gtk_widget_ref (label56);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label56", label56, NULL);
  gtk_widget_set_visible (label56, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), label56);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label56, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_c_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_c_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_c_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_c_entry_1, NULL);
  gtk_widget_set_visible (phs_cell_c_entry_1, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append (GTK_BOX (hbox33), phs_cell_c_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_c_entry_1, TRUE, TRUE, 0);
#endif
  gtk_widget_set_size_request (phs_cell_c_entry_1, 65, -2);

  label57 = gtk_label_new (("alpha"));
  // gtk_widget_ref (label57);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label57", label57,
			  NULL);
  gtk_widget_set_visible (label57, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), label57);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label57, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_alpha_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_alpha_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_alpha_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_alpha_entry_1,
			  NULL);
  gtk_widget_set_visible (phs_cell_alpha_entry_1, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), phs_cell_alpha_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_alpha_entry_1, TRUE, TRUE, 0);
#endif
  gtk_widget_set_size_request (phs_cell_alpha_entry_1, 60, -2);

  label58 = gtk_label_new (("beta"));
  // gtk_widget_ref (label58);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label58", label58,
			  NULL);
  gtk_widget_set_visible (label58, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), label58);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label58, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_beta_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_beta_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_beta_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_beta_entry_1, NULL);
  gtk_widget_set_visible (phs_cell_beta_entry_1, TRUE);
#if  (GTK_MAJOR_VERSION == 4)
  gtk_box_append (GTK_BOX (hbox33), phs_cell_beta_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_beta_entry_1, TRUE, TRUE, 0);
#endif
  gtk_widget_set_size_request (phs_cell_beta_entry_1, 60, -2);

  label59 = gtk_label_new (("gamma"));
  // gtk_widget_ref (label59);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label59", label59,
			  NULL);
  gtk_widget_set_visible (label59, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), label59);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), label59, FALSE, FALSE, 2);
#endif


  strcpy(widget_name, "phs_cell_gamma_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_gamma_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_gamma_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_gamma_entry_1, NULL);
  gtk_widget_set_visible (phs_cell_gamma_entry_1, TRUE);
#if (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox33), phs_cell_gamma_entry_1);
#else
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_gamma_entry_1, TRUE, TRUE, 0);
#endif
  gtk_widget_set_size_request (phs_cell_gamma_entry_1, 60, -2);

  return phs_cell_group;

}


/* ------------------------------------------------------------------------------ */
/*     None of the above cell chooser for phs                                     */
/* ------------------------------------------------------------------------------ */

void display_none_cell_chooser_box(GtkWidget *phs_cell_choice_window,
				   GSList *phs_cell_group) {

  GtkWidget *hbox34;
  GtkWidget *phs_cell_chooser_vbox;
  GtkWidget *phs_cell_none_radiobutton = nullptr;

  phs_cell_chooser_vbox = widget_from_builder("phs_cell_chooser_vbox");

  hbox34 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox34);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "hbox34", hbox34, NULL);
  gtk_widget_set_visible (hbox34, TRUE);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX(phs_cell_chooser_vbox), hbox34);
#else
  gtk_box_pack_start (GTK_BOX (phs_cell_chooser_vbox), hbox34, TRUE, TRUE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (hbox34), 6);
#endif

#if (GTK_MAJOR_VERSION >= 4)
  // 20220602-PE FIXME radio buttons
#else
  phs_cell_none_radiobutton = gtk_radio_button_new_with_label (phs_cell_group, _("None of the Above"));
  phs_cell_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (phs_cell_none_radiobutton));
#endif
  // gtk_widget_ref (phs_cell_none_radiobutton);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "phs_cell_none_radiobutton",
			  phs_cell_none_radiobutton,
			  NULL);

  gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (phs_cell_none_radiobutton), TRUE);

  gtk_widget_set_visible (phs_cell_none_radiobutton, TRUE);

#if  (GTK_MAJOR_VERSION == 4)
  gtk_box_append(GTK_BOX (hbox34), phs_cell_none_radiobutton);
#else
  gtk_box_pack_start (GTK_BOX (hbox34), phs_cell_none_radiobutton, FALSE, FALSE, 4);
#endif

}



/* ------------------------------------------------------------------------------ */
/* Additional Representation Handling                                             */
/* ------------------------------------------------------------------------------ */

// Return the vbox for the Add Reps, (which will then be used in
// display_control_add_reps().
//
// It's a bit ugly to send back a widget to the molecule_info_class_t
// function.  If you want to remove that in future (make this function
// a void), then you'll have to set a name on the vbox that can be
// looked up by display_control_add_reps().
//
//
// Function no longer needed, frame is made by
// display_control_molecule_combo_box().
//
GtkWidget *
display_control_add_reps_container(GtkWidget *display_control_window_glade,
				   int imol_no) {

   GtkWidget *w = NULL;

   if (display_control_window_glade) {
      std::string name = "add_rep_display_control_frame_vbox_";
      name += coot::util::int_to_string(imol_no);
      // GtkWidget *t = lookup_widget(display_control_window_glade, name.c_str());
      GtkWidget *t = 0;
      std::cout << "display_control_add_reps_container(): Do a proper lookup of t here" << std::endl;
      if (t)
	 w = t;
      else
	 std::cout << "ERROR:: in display_control_add_reps_frame failed to lookup "
		   << name << " widget" << std::endl;
   }
   return w;
}

GtkWidget *
display_control_add_reps_frame(GtkWidget *display_control_window_glade,
			       int imol_no) {

   GtkWidget *w = NULL;

   if (display_control_window_glade) {
      std::string name = "add_rep_display_control_frame_";
      name += coot::util::int_to_string(imol_no);
      // GtkWidget *t = lookup_widget(display_control_window_glade, name.c_str());
      GtkWidget *t = 0;
      std::cout << "display_control_add_reps_frame(): Do a proper lookup of t here" << std::endl;
      if (t)
	 w = t;
      else
	 std::cout << "ERROR:: in display_control_add_reps_frame failed to lookup "
		   << name << " widget" << std::endl;
   }
   return w;
}

GtkWidget *
display_control_add_reps_all_on_check_button(GtkWidget *display_control_window_glade,
			       int imol_no) {
   GtkWidget *w = NULL;
   if (display_control_window_glade) {
      std::string name = "add_rep_all_on_check_button_";
      name += coot::util::int_to_string(imol_no);
      // GtkWidget *t = lookup_widget(display_control_window_glade, name.c_str());
      GtkWidget *t = 0;
      std::cout << "display_control_add_reps_all_on_check_button(): Do a proper lookup of t here" << std::endl;
      if (t)
	 w = t;
      else
	 std::cout << "ERROR:: in display_control_add_reps_all_on_check_button failed to lookup "
		   << name << " widget" << std::endl;
   }
   return w;
}

void display_control_add_reps(GtkWidget *display_control_window_glade,
			      int imol_no, int add_rep_no,
			      bool show_it,
			      int bonds_box_type, const std::string &name) {

   if (display_control_window_glade) {
      GtkWidget *add_rep_vbox  = display_control_add_reps_container(display_control_window_glade, imol_no);
      GtkWidget *add_rep_frame = display_control_add_reps_frame(    display_control_window_glade, imol_no);
      GtkWidget *add_rep_all_on_checkbutton =
	 display_control_add_reps_all_on_check_button(display_control_window_glade, imol_no);

//       std::cout <<  "DEBUG::  add_rep_vbox is " << add_rep_vbox  << std::endl;
//       std::cout <<  "DEBUG:: add_rep_frame is " << add_rep_frame << std::endl;

      GtkWidget *hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
      // gtk_widget_ref (hbox);
#if (GTK_MAJOR_VERSION == 4)
      gtk_box_append(GTK_BOX(add_rep_vbox), hbox);
#else
      gtk_box_pack_start(GTK_BOX(add_rep_vbox), hbox, FALSE, FALSE, 0);
#endif
      std::string label = name;
      GtkWidget *toggle_button_show_it = gtk_check_button_new_with_label(label.c_str());
      // gtk_widget_ref (toggle_button_show_it);
      if (show_it) {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button_show_it), TRUE);
	 gtk_widget_set_visible(add_rep_all_on_checkbutton, TRUE);
      } else {
	 gtk_widget_set_visible(add_rep_all_on_checkbutton, FALSE);
      }
      int cc = encode_ints(imol_no, add_rep_no);
      g_signal_connect(G_OBJECT(toggle_button_show_it), "toggled",
			 G_CALLBACK(add_rep_toggle_button_toggled),
			 GINT_TO_POINTER(cc));
#if (GTK_MAJOR_VERSION == 4)
      gtk_box_append(GTK_BOX(hbox), toggle_button_show_it);
#else
      gtk_box_pack_start(GTK_BOX(hbox), toggle_button_show_it, FALSE, FALSE, 0);
#endif
      gtk_widget_set_visible(toggle_button_show_it, TRUE);
      gtk_widget_set_visible(hbox, TRUE);
      gtk_widget_set_visible(add_rep_vbox, TRUE);
      gtk_widget_set_visible(add_rep_frame, TRUE);
   }
}

void add_rep_toggle_button_toggled(GtkToggleButton       *button,
				   gpointer         user_data) {

   std::pair<int, int> p = decode_ints(GPOINTER_TO_INT(user_data));
   int on_off_flag = 0;
   if (gtk_toggle_button_get_active(button))
      on_off_flag = 1;
   set_show_additional_representation(p.first, p.second, on_off_flag);
}

#if 0
GtkWidget*
create_splash_screen_window_for_file(const char *file_name) {

   // I have another version of this in main.cc? Hmm...

   GtkWidget *splash_screen_window = gtk_window_new();
   gtk_widget_set_name (splash_screen_window, "splash_screen_window");
   std::string window_name = "Coot " + std::string(VERSION);
   gtk_window_set_title (GTK_WINDOW (splash_screen_window), _(window_name.c_str()));
#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "in create_splash_screen_window_for_file() set position " << std::endl;
#else
   gtk_window_set_position (GTK_WINDOW (splash_screen_window), GTK_WIN_POS_CENTER);
#endif

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   std::cout << "in create_splash_screen_window_for_file() what has this become now?" << std::endl;
   // gtk_window_set_type_hint (GTK_WINDOW (splash_screen_window), GDK_SURFACE_TYPE_HINT_SPLASHSCREEN);
#else
   gtk_window_set_type_hint (GTK_WINDOW (splash_screen_window), GDK_WINDOW_TYPE_HINT_SPLASHSCREEN);
#endif

   GtkWidget *image6807 = create_pixmap(splash_screen_window, file_name); // create_pixmap() I can keep/copy over
   gtk_widget_set_visible (image6807, TRUE);

   gtk_window_set_child(GTK_WINDOW(splash_screen_window), image6807);

   return splash_screen_window;
}
#endif


#ifdef COOT_USE_GTK2_INTERFACE
/* Hack in a function, because it's missing somehow from Bernhard's
   commits */

/* I dont think this should be here at all!
   I think it should be in gtk2-interface.c !
   But somehow it disappeared from there, so I will put it back!!!
   Where did all my things go........ */

#endif	/* COOT_USE_GTK2_INTERFACE */
