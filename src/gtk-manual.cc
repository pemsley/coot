/* src/gtk-manual.cc
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


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#include <string.h>

// NOTE:: The order of these 6 include files seems fragile
#include <gtk/gtk.h>
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"

#include "gtk-manual.h"
#include "gtk-manual.hh"

#include "interface.h"  /* for create_single_map_properties_dialog() */
#include "utils/coot-utils.hh"

#include "widget-from-builder.hh"


/* This is the signal handler for a color change event created when
   the colorseldialog has had its colour changed. */

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                                 map                                      */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

void
on_map_color_changed(GtkWidget *w,
		     gpointer tmd) {

   struct map_colour_data_type* t = static_cast<struct map_colour_data_type*> (tmd);
   GdkColor color;
   gtk_color_selection_get_current_color(t->color_selection, &color);
   GdkRGBA map_color;
   map_color.red   = color.red    /65535.0;;
   map_color.green = color.green  /65535.0;;
   map_color.blue  = color.blue   /65535.0;;
   handle_map_colour_change(t->imol, map_color);

}


/* 		 GtkColorSelection *cs); */
/*  The colour selection dialog has had its OK button pressed */
void
on_map_col_sel_ok_button_clicked        (GtkButton       *button,
					 gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_destroy(w);
}

/*  The colour selection dialog has had its cancel button pressed */
void
on_map_col_sel_cancel_button_clicked (GtkButton       *button,
				      gpointer         *tmd)
{

  struct map_colour_data_type* t = (struct map_colour_data_type*) tmd;
  int imol = t->imol;
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "map_color_selection_dialog");

  /* we should put the colour back to
     how it used to be then. */
  restore_previous_map_colour(t->imol);
  gtk_widget_destroy(dialog);
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


void
on_symmetry_color_changed(GtkWidget *w,
			  GtkColorSelection *colorsel) {
   gdouble color[4];

   // gtk_color_selection_get_color(colorsel,color);

   /* we pass back the model number */
   handle_symmetry_colour_change(1,color);
}

/*  The colour selection dialog has had its OK button pressed */
void
on_symm_col_sel_ok_button_clicked (GtkButton       *button,
				   gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_destroy(w);
}

/*  The colour selection dialog has had its cancel button pressed */
void
on_symm_col_sel_cancel_button_clicked (GtkButton       *button,
				      gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_destroy(w);
				/* we should put the colour back to
				   how it used to be then. */
}

/* ----------------------------------------------------------------- */
/* dynamic map color menu */
/* ----------------------------------------------------------------- */

void
create_initial_map_color_submenu(GtkWidget *window1) {

   GtkWidget *map_colour1 = widget_from_builder("map_colour1");
   GtkWidget *map_colour1_menu = gtk_menu_new();

   // gtk_widget_ref (map_colour1_menu);
   g_object_set_data_full (G_OBJECT (window1), "map_colour1_menu", map_colour1_menu, NULL);
   gtk_menu_item_set_submenu (GTK_MENU_ITEM (map_colour1), map_colour1_menu);
}



void
create_initial_ramachandran_mol_submenu(GtkWidget *widget) {


 /* We need to get to window1 */

   GtkWidget *window1 = widget;

   // rama_draw = GTK_WIDGET(lookup_widget(window1, "ramachandran_plot1"));
   GtkWidget *rama_draw = widget_from_builder("ramachandran_plot1");

   GtkWidget *rama_draw_menu = gtk_menu_new();

   g_object_set_data_full (G_OBJECT (window1), "rama_plot_menu", rama_draw_menu, NULL); // 20211002-PE what does this do now?

   gtk_menu_item_set_submenu (GTK_MENU_ITEM (rama_draw), rama_draw_menu);
   g_object_set_data(G_OBJECT(rama_draw), "rama_plot_menu", rama_draw_menu);
}


void
update_ramachandran_plot_menu_manual(int imol, const char *name) {

   GtkWidget *menu_item;
   GtkWidget *window1;
   GtkWidget *rama_plot_menu;
   char *text;

   window1 = GTK_WIDGET(lookup_widget(main_window(), "window1"));

   menu_item = gtk_menu_item_new_with_label (name);

   rama_plot_menu = GTK_WIDGET(lookup_widget(window1, "rama_plot_menu"));

   // gtk_widget_ref (menu_item);
   g_object_set_data_full (G_OBJECT (window1), "rama_plot_menu_item",
			   menu_item, NULL);

  gtk_widget_show (menu_item);
  gtk_container_add (GTK_CONTAINER (rama_plot_menu), menu_item);

  g_signal_connect (G_OBJECT (menu_item), "activate",
		    G_CALLBACK (rama_plot_mol_selector_activate),
		    GINT_TO_POINTER(imol));
}


void
rama_plot_mol_selector_activate (GtkMenuItem     *menuitem,
				 gpointer         user_data)
{
  int imol = GPOINTER_TO_INT(user_data);
  GtkWidget *rama_widget = 0;
/*   printf("selector activate: do rama plot for molecule: %d\n", *imol); */

/* We should come here and be given imol.  New molecules should insert
   themselves into the Ramachandran Plot menu(item). */


  rama_widget = dynarama_is_displayed_state(imol);
  if (rama_widget == NULL) {
    do_ramachandran_plot(imol);
  } else {
     // if (!GTK_WIDGET_MAPPED(rama_widget))
     if (true) {
        gtk_widget_show(rama_widget);
     } else {
#if (GTK_MAJOR_VERSION < 4)
        gdk_window_raise(GDK_WINDOW(gtk_widget_get_window(rama_widget)));
#endif
     }
  }

}

/* And similar for sequence view: */
void create_initial_sequence_view_mol_submenu(GtkWidget *widget) {

   GtkWidget *seq_view_draw = widget_from_builder("sequence_view1");
   GtkWidget *seq_view_menu = gtk_menu_new();

   g_object_set_data_full(G_OBJECT(widget), "seq_view_menu", seq_view_menu, NULL); // 20211002-PE what does this do (likewsie to rama)
   gtk_menu_item_set_submenu (GTK_MENU_ITEM(seq_view_draw), seq_view_menu);
   g_object_set_data(G_OBJECT(seq_view_draw), "seq_view_menu", seq_view_menu);
}

void update_sequence_view_menu_manual(int imol, const char *name) {

   char *text;
   GtkWidget *window1 = lookup_widget(main_window(), "window1");
   GtkWidget *seq_view_menu = lookup_widget(window1, "seq_view_menu");
   GtkWidget *menu_item;

   menu_item = gtk_menu_item_new_with_label (name);
   // gtk_widget_ref (menu_item);
   g_object_set_data_full (G_OBJECT(window1), "seq_view_menu_item",
			   menu_item,
			   NULL);
   gtk_widget_show(menu_item);
   gtk_container_add(GTK_CONTAINER(seq_view_menu), menu_item);
   g_signal_connect (G_OBJECT(menu_item), "activate",
		     G_CALLBACK(sequence_view_mol_selector_activate),
		     GINT_TO_POINTER(imol));
}

void sequence_view_mol_selector_activate (GtkMenuItem     *menuitem,
					  gpointer         user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  std::cout << "calling do_sequence_view() " << imol  << std::endl;
   do_sequence_view(imol);

}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*                               skeleton                                   */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

#include "c-interface-widgets.hh"

void
on_skeleton_color_changed(GtkWidget *w,
			  GtkColorSelection *colorsel) {
   gdouble color[4];
   for (int i=0; i<4; i++) color[i] = 0.0;

   std::cout << "fix the colour" << std::endl;

   // gtk_color_selection_get_color(colorsel,color);

   /* we pass back the model number */
   handle_skeleton_colour_change(1,color);
}



/*  The colour selection dialog has had its OK button pressed */
void
on_skeleton_col_sel_ok_button_clicked (GtkButton       *button,
				       gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_destroy(w);
}

/*  The colour selection dialog has had its cancel button pressed */
void
on_skeleton_col_sel_cancel_button_clicked (GtkButton       *button,
				      gpointer         user_data)
{
   GtkWidget *w = GTK_WIDGET(user_data);
   gtk_widget_destroy(w);
				/* we should put the colour back to
				   how it used to be then. */
}


void
on_display_manager_selections_and_colours_combobox_changed(GtkComboBox     *combo_box,
                                                           gpointer         user_data) {

   std::cout << "DEBUG:: :::::: on_display_manager_selections_and_colours_combobox_changed() "
             << std::endl;

   int imol = GPOINTER_TO_INT(user_data);
   gchar *txt = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo_box));
   std::cout << "DEBUG:: text: \"" << txt << "\" user data " << imol << std::endl;

   if (txt) {
      std::string at(txt);
      if (at == _("Bonds (Colour by Atom)")) {
         std::cout << "display as bonds " <<std::endl;
         graphics_to_bonds_representation(imol);
      }
      if (at == _("C-alphas/Backbone")) {
         std::cout << "display as CA " <<std::endl;
         graphics_to_ca_representation(imol);
      }
      if (at == _("Bonds (Colour by Molecule)")) {
         render_as_bonds_colored_by_molecule_button_select(imol);
      }
      if (at == _("Bonds (Colour by Chain)")) {
         render_as_bonds_colored_by_chain_button_select(imol);
      }
      if (at == _("Bonds (Goodsell Colour by Chain)")) {
         render_as_bonds_goodsell_colored_by_chain_button_select(imol);
      }
      if (at == _("Bonds (Colour by Sec. Str.)")) {
         render_as_sec_struct_bonds_button_select(imol);
      }
      if (at == _("CAs + Ligands")) {
         render_as_ca_plus_ligands_bonds_button_select(imol);
      }
      if (at ==  _("CAs+Ligs SecStr Col")) {
         render_as_ca_plus_ligands_sec_str_bonds_button_select(imol);
      }
      if (at == _("Jones' Rainbow")) {
         render_as_rainbow_representation_button_select(imol);
      }
      if (at ==  _("Colour by Atom - No Waters")) {
         render_as_bonds_no_waters(imol);
      }
      if (at ==  _("Colour by B-factor - CAs")) {
         render_as_b_factor_cas_representation_button_select(imol);
      }
      if (at ==  _("Colour by B-factor - All")) {
         render_as_b_factor_representation_button_select(imol);
      }
      if (at ==  _("Colour by Occupancy")) {
         render_as_occupancy_representation_button_select(imol);
      }
   }
}


// make and show (will need to be packed)
//
GtkWidget *selections_and_colours_combobox(int imol) {

   GtkWidget *combobox = gtk_combo_box_text_new();

   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Bonds (Colour by Atom)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Bonds (Colour by Molecule)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Bonds (Colour by Chain)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Bonds (Goodsell Colour by Chain)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Bonds (Colour by Sec. Str.)"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("C-alphas/Backbone"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("CAs + Ligands"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("CAs+Ligs SecStr Col"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Jones' Rainbow"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Colour by Atom - No Waters"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Colour by B-factor - Backbone"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Colour by B-factor - All"));
   gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), _("Colour by Occupancy"));

   gpointer user_data = GINT_TO_POINTER(imol);
   g_signal_connect(G_OBJECT(combobox), "changed",
                    G_CALLBACK(on_display_manager_selections_and_colours_combobox_changed), user_data);

   int index = 0; // needs to be dynamic (at some stage)

   int bbt = graphics_molecule_bond_type(imol);

   // std::cout << "Here with bbt " << bbt << std::endl;
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
   if (bbt ==  8) index =  1;
   if (bbt ==  3) index =  2;
   if (bbt == 21) index =  3;
   if (bbt ==  6) index =  4;
   if (bbt ==  2) index =  5;
   if (bbt ==  4) index =  6;
   if (bbt ==  7) index =  7;
   if (bbt ==  9) index =  8;
   if (bbt == 11) index = 12;

   gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), index);
   gtk_widget_show(combobox);
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
void display_control_molecule_combo_box(GtkWidget *display_control_window_glade,
					const gchar *name,
					int imol, bool show_add_reps_frame_flag) {

   std::cout << "start display_control_molecule_combo_box() " << std::endl;
   GtkWidget *display_molecule_vbox;
   /*   GtkWidget *display_control_window_glade; passed parameter */

   GtkWidget *display_mol_frame_1;
   GtkWidget *hbox31;
   GtkWidget *hbox32;
   GtkWidget *entry2;
   GtkWidget *displayed_button_1;
   GtkWidget *active_button_1;
   GtkWidget *render_optionmenu_1;
   GtkWidget *render_optionmenu_1_menu;
   GtkWidget *glade_menuitem;
   // GtkWidget *menu;
   int bond_type;
   GtkWidget *active_item;
   GtkWidget *mol_number_label;

   /* messing about with string variables for unique lookup values/name of the widgets */
   std::string widget_name;

  display_molecule_vbox = lookup_widget(display_control_window_glade,
					"display_molecule_vbox");

  display_mol_frame_1 = gtk_frame_new (NULL);
  // gtk_widget_ref (display_mol_frame_1);

  widget_name = "display_mol_frame_";
  int nn = 3;
  if (imol > 9) nn = 2;
  if (imol > 99) nn = 1;
  std::string four_char_imol(nn, '0');
  // widget_name += four_char_imol + coot::util::int_to_string(imol);
  widget_name += coot::util::int_to_string(imol);

  std::cout << "debug:: frame widget_name " << widget_name << std::endl;

  g_object_set_data_full (G_OBJECT (display_control_window_glade),
			  widget_name.c_str(),
			  display_mol_frame_1,
                          NULL);
  gtk_widget_show (display_mol_frame_1);
  gtk_box_pack_start (GTK_BOX (display_molecule_vbox), display_mol_frame_1,
		      FALSE, FALSE, 0);

  hbox31 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox31);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), "hbox31", hbox31, NULL);
  gtk_widget_show (hbox31);

  GtkWidget *vbox_single_molecule_all_attribs = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  gtk_widget_show(vbox_single_molecule_all_attribs);
  gtk_container_add (GTK_CONTAINER (display_mol_frame_1), vbox_single_molecule_all_attribs);
  gtk_box_pack_start(GTK_BOX(vbox_single_molecule_all_attribs), hbox31, FALSE, FALSE, 2);


/* -- molecule number label */

  widget_name = "display_mol_number_";
  widget_name += four_char_imol + coot::util::int_to_string(imol);

  mol_number_label = gtk_label_new (coot::util::int_to_string(imol).c_str());;
  g_object_set_data_full (G_OBJECT (display_control_window_glade),
			  widget_name.c_str(), mol_number_label, NULL);
  gtk_widget_show (mol_number_label);
  gtk_box_pack_start (GTK_BOX (hbox31), mol_number_label, FALSE, FALSE, 3);

/* -- done molecule number label */

  widget_name = "display_mol_entry_";
  widget_name += four_char_imol + coot::util::int_to_string(imol);

  entry2 = gtk_entry_new ();
  // gtk_widget_ref (entry2);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), widget_name.c_str(), entry2,
			  NULL);
  if (name) {
    gtk_entry_set_text(GTK_ENTRY(entry2), name);
    /* these 2 seem not to do what I want :-( */
    // gtk_entry_set_position(GTK_ENTRY(entry2), strlen(name)-1);
    // gtk_entry_append_text(GTK_ENTRY(entry2), "");
  }
  gtk_editable_set_editable(GTK_EDITABLE (entry2), FALSE);


  gtk_widget_show (entry2);
  gtk_box_pack_start (GTK_BOX (hbox31), entry2, TRUE, TRUE, 3);

  hbox32 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox32);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), "hbox32", hbox32,
                            NULL);
  gtk_widget_show (hbox32);
  gtk_box_pack_start (GTK_BOX (hbox31), hbox32, TRUE, TRUE, 0);

  widget_name = "display_mol_button_";
  widget_name += four_char_imol + coot::util::int_to_string(imol);

  displayed_button_1 = gtk_check_button_new_with_label (_("Display"));
  // gtk_widget_ref (displayed_button_1);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(displayed_button_1),
			       mol_is_displayed(imol));

  g_object_set_data_full (G_OBJECT (display_control_window_glade),
                          widget_name.c_str(), displayed_button_1, NULL);

  gtk_widget_show (displayed_button_1);
  gtk_box_pack_start (GTK_BOX (hbox32), displayed_button_1, FALSE, FALSE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (displayed_button_1), 2);

  widget_name = "active_mol_button_";
  widget_name += four_char_imol + coot::util::int_to_string(imol);

  active_button_1 = gtk_check_button_new_with_label (_("Active"));
  // gtk_widget_ref (active_button_1);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(active_button_1), mol_is_active(imol));
  g_object_set_data_full (G_OBJECT (display_control_window_glade),
                          widget_name.c_str(), active_button_1, NULL);
  gtk_widget_show (active_button_1);
  gtk_box_pack_start (GTK_BOX (hbox32), active_button_1, FALSE, FALSE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (active_button_1), 2);

  GtkWidget *sel_and_col_combobox = selections_and_colours_combobox(imol);
  gtk_box_pack_start(GTK_BOX(hbox32), sel_and_col_combobox, FALSE, FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(sel_and_col_combobox), 4);


/* Set User Data, the molecule which this button(s) is attached to
   (casting (int *) to (char *)).
*/
  g_object_set_data(G_OBJECT(displayed_button_1), "imol", GINT_TO_POINTER(imol));
  g_object_set_data(G_OBJECT(   active_button_1), "imol", GINT_TO_POINTER(imol));

/* Add signals for the Active and Display toggle buttons */

  g_signal_connect(G_OBJECT (displayed_button_1), "toggled",
		     G_CALLBACK (on_display_control_mol_displayed_button_toggled),
		     NULL);

  g_signal_connect(G_OBJECT (active_button_1), "toggled",
		     G_CALLBACK (on_display_control_mol_active_button_toggled),
		     NULL);



/* Set User Data, the molecule which this button(s) is attached to
   (casting (int *) to (char *)).
*/
  g_object_set_data(G_OBJECT(displayed_button_1), "imol", GINT_TO_POINTER(imol));
  g_object_set_data(G_OBJECT(   active_button_1), "imol", GINT_TO_POINTER(imol));



  /* Which menu item should be displayed in the selector?  If we
     turned to a C-alpha representation and the closed the display
     control window.  when it comes back, we want it to be C-alpha
     too, not the default (1) - "Bonds"  */

  bond_type = graphics_molecule_bond_type(imol);

/* And finally connect the menu to the optionmenu */
  // gtk_option_menu_set_menu (GTK_OPTION_MENU (render_optionmenu_1),
  // render_optionmenu_1_menu);

  // A delete molecule button
  display_control_add_delete_molecule_button(imol, hbox32, false);

  // Now add the additional representations frame and vbox
  add_add_reps_frame_and_vbox(display_control_window_glade,
			      vbox_single_molecule_all_attribs, imol, show_add_reps_frame_flag);

}

void display_control_add_delete_molecule_button(int imol, GtkWidget *hbox32,
						short int is_map_molecule) {

   std::string delete_button_name = "delete_molecule_";
   delete_button_name += coot::util::int_to_string(imol);
   std::string button_string = "Delete Model";
   if (is_map_molecule)
      button_string = "Delete Map";
   GtkWidget *delete_button = gtk_button_new_with_label(_(button_string.c_str()));
   gtk_widget_show(delete_button);
   gtk_box_pack_start (GTK_BOX (hbox32), delete_button, FALSE, FALSE, 0);
   gtk_container_set_border_width (GTK_CONTAINER (delete_button), 2);
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
   GtkWidget *v = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
   // show the frame if there were additional representations
   if (show_add_reps_frame_flag)
      gtk_widget_show(frame);
   // gtk_widget_ref (v);
   // gtk_widget_ref (frame);

   std::string arl = "   Show Additional Representations  ";
   GtkWidget *all_on_check_button = gtk_check_button_new_with_label(arl.c_str());
   if (show_add_reps_frame_flag)
      gtk_widget_show(all_on_check_button);
   gtk_box_pack_start(GTK_BOX(hbox_for_single_molecule), all_on_check_button, FALSE, FALSE, 2);
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
   g_object_set_data_full (G_OBJECT (display_control_window_glade),
			     widget_name.c_str(),
			   v, NULL);

   widget_name = "add_rep_display_control_frame_";
   widget_name += coot::util::int_to_string(imol_no);
   g_object_set_data_full (G_OBJECT (display_control_window_glade),
			     widget_name.c_str(),
			   frame, NULL);

   gtk_container_add(GTK_CONTAINER(hbox_for_single_molecule), frame);
   gtk_container_add(GTK_CONTAINER(frame), v);

}


void
on_add_rep_all_on_check_button_toggled   (GtkToggleButton       *button,
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


void
update_name_in_display_control_molecule_combo_box(GtkWidget *display_control_window_glade,
						  const gchar *display_name,
						  int n) {
  int i;
  char entry_name[1024];
  GtkWidget *entry;
  gchar *tmp_name;
  int imol = n;

  for (i=0; i<1024; i++)
    entry_name[i]= 0;
  if (is_valid_map_molecule(imol)) {
    memcpy(entry_name, "display_map_entry_", 18);
      } else {
    memcpy(entry_name, "display_mol_entry_", 18);
      }
  tmp_name = entry_name + strlen(entry_name);
  snprintf(tmp_name, 4, "%-d", imol);

/*   printf("debug:: molecule number (dereferenced): %d\n", *n); */
/*   printf("debug:: pointer string: %s\n", tmp_name); */
/*   printf("debug:: searching for entry name %s\n", entry_name); */
  entry = lookup_widget(display_control_window_glade, entry_name);

  if (entry)
    gtk_entry_set_text(GTK_ENTRY(entry), display_name);
  else
    printf("oops no entry found with name %s\n", entry_name);
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
GtkWidget *display_control_map_combo_box(GtkWidget *display_control_window_glade,
					 const gchar *name,
					 int n) {

  GtkWidget *my_combo_box = 0;

  GtkWidget *display_map_vbox;
/*   GtkWidget *display_control_window_glade; passed parameter */

  GtkWidget *display_map_frame_1;
  GtkWidget *hbox31;
  GtkWidget *hbox32;
  GtkWidget *entry2;
  GtkWidget *displayed_button_1;
  GtkWidget *mol_label;

/* messing about with string variables */
  gchar *widget_name;
  gchar *tmp_name;

  GtkWidget *scroll_radio_button_1;

  /* we need to find references to objects that contain this widget */

  display_map_vbox = lookup_widget(display_control_window_glade,
				   "display_map_vbox");

  display_map_frame_1 = gtk_frame_new (NULL);
  // gtk_widget_ref (display_map_frame_1);

  widget_name = (gchar *) malloc(100); /* should be enough */

  strcpy(widget_name, "display_map_frame_");
  tmp_name = widget_name + strlen(widget_name);

  snprintf(tmp_name, 4, "%-d", n);

/*   printf("display_map_frame_{thing} name constructed as: :%s:\n", widget_name);  */


  g_object_set_data_full (G_OBJECT (display_control_window_glade), widget_name,
			  display_map_frame_1, NULL);
  gtk_widget_show (display_map_frame_1);
  /* setting to true means that the buttons etc in the box can expand
     vertically to fill the box  */
  gtk_box_pack_start (GTK_BOX (display_map_vbox), display_map_frame_1, FALSE, FALSE, 0);

  hbox31 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox31);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), "hbox31", hbox31, NULL);

  gtk_widget_show (hbox31);
  gtk_container_add (GTK_CONTAINER (display_map_frame_1), hbox31);

/* -- molecule number label */

  strcpy(widget_name, "display_map_number_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);

  mol_label = gtk_label_new (_(tmp_name));
  // gtk_widget_ref (mol_label);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), widget_name, mol_label,
			  NULL);

  gtk_widget_show (mol_label);
  gtk_box_pack_start (GTK_BOX (hbox31), mol_label, FALSE, FALSE, 3);

/* -- molecule number label */

/* -- */

  strcpy(widget_name, "display_map_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);

  entry2 = gtk_entry_new ();
  // gtk_widget_ref (entry2);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), widget_name, entry2,
			  NULL);
  if (name) {
    gtk_entry_set_text(GTK_ENTRY(entry2), name);
  }

  std::cout << "set entry not editable" << std::endl;
  // gtk_entry_set_editable(GTK_ENTRY (entry2), FALSE);

  gtk_widget_show (entry2);
  gtk_box_pack_start (GTK_BOX (hbox31), entry2, TRUE, TRUE, 0);

  hbox32 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox32);
  g_object_set_data_full (G_OBJECT (display_control_window_glade), "hbox32", hbox32,
			  NULL);
  gtk_widget_show (hbox32);
  gtk_box_pack_start (GTK_BOX (hbox31), hbox32, TRUE, TRUE, 0);

/* -- */

  strcpy(widget_name, "displayed_button_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);

  displayed_button_1 = gtk_check_button_new_with_label (_("Display"));
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(displayed_button_1), map_is_displayed(n));
  // gtk_widget_ref (displayed_button_1);
  g_object_set_data_full (G_OBJECT (display_control_window_glade),
			    widget_name,
			    displayed_button_1,
			  NULL);
  gtk_widget_show (displayed_button_1);
  gtk_box_pack_start (GTK_BOX (hbox32), displayed_button_1, FALSE, FALSE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (displayed_button_1), 2);

  g_signal_connect(G_OBJECT (displayed_button_1),  "toggled",
		     G_CALLBACK (on_display_control_map_displayed_button_toggled),
		     GINT_TO_POINTER(n));

/*   // associate with the button a pointer to the variable which */
/*   // contains the passed variable int n */
/*   //  */
/*   // we cast as a (char *) to make the compiler happy. */
/*   //  */
    g_object_set_data (G_OBJECT (displayed_button_1), "imol", GINT_TO_POINTER(n));

/* -- */
  /* 20050316 Today I add scroll check-button, as Charlie asked for ages ago. */

  strcpy(widget_name, "map_scroll_button_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);

  GtkWidget *previous_radio_button =
     get_radio_button_in_scroll_group(display_control_window_glade, n);

  if (previous_radio_button)
     scroll_radio_button_1 =
	gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(previous_radio_button),
						    _("Scroll"));
  else
     scroll_radio_button_1 = gtk_radio_button_new_with_label(NULL, _("Scroll"));

  g_object_set_data_full(G_OBJECT(display_control_window_glade),
			   widget_name,
 			   scroll_radio_button_1,
 			   NULL);


  gtk_widget_show(scroll_radio_button_1);
  gtk_box_pack_start(GTK_BOX(hbox32), scroll_radio_button_1, FALSE,FALSE, 2);
  g_signal_connect(G_OBJECT(scroll_radio_button_1), "toggled",
		     G_CALLBACK (on_display_control_map_scroll_radio_button_toggled),
		     GINT_TO_POINTER(n));
  g_signal_connect(G_OBJECT(scroll_radio_button_1), "group_changed",
		     G_CALLBACK (on_display_control_map_scroll_radio_button_group_changed),
		     GINT_TO_POINTER(n));

  if (scroll_wheel_map() == n) {
     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(scroll_radio_button_1), TRUE);
  } else {
     if (previous_radio_button) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(scroll_radio_button_1), FALSE);
     }
  }


/* -- */

  strcpy(widget_name, "properties_button_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);

  displayed_button_1 = gtk_button_new_with_label (_("Properties"));
  // gtk_widget_ref (displayed_button_1);
  g_object_set_data_full (G_OBJECT (display_control_window_glade),
			  widget_name,
			  displayed_button_1,
			  NULL);
  gtk_widget_show (displayed_button_1);
  gtk_box_pack_start (GTK_BOX (hbox32), displayed_button_1, FALSE, FALSE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (displayed_button_1), 2);

  g_signal_connect(G_OBJECT (displayed_button_1),  "clicked",
		     G_CALLBACK (on_display_control_map_properties_button_clicked),
		     GINT_TO_POINTER(n));

/*   // associate with the button a pointer to the variable which */
/*   // contains the passed variable int n */
/*   //  */
  /*   // we cast as a (char *) to make the compiler happy. */
/*   //  */
    g_object_set_data (G_OBJECT (displayed_button_1), "imol", GINT_TO_POINTER(n));

  // A delete molecule button
  display_control_add_delete_molecule_button(n, hbox32, true);

  free(widget_name);
  return my_combo_box;
}


void
on_display_control_map_displayed_button_toggled   (GtkToggleButton       *button,
						   gpointer         user_data)
{

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "imol"));

   if (gtk_toggle_button_get_active(button))
      set_map_displayed(imol, 1);
   else
      set_map_displayed(imol, 0);
}

/* Added 20050316 (Bangalore) */
void
on_display_control_map_scroll_radio_button_toggled (GtkToggleButton *button,
						    gpointer         user_data) {
  int imol = GPOINTER_TO_INT(user_data);
  const char *state = "inactive";
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
    state = "active";
    set_scrollable_map(imol);
  }
}

void
on_display_control_map_scroll_radio_button_group_changed (GtkRadioButton *button,
							  gpointer         user_data) {
   // do nothing these days. Scroll group is handled gtk-internally via radio buttons.
}

// return NULL when there are no radio buttons found that have a group
// (i.e. there are no maps displayed in the Display Manager.
//
GtkWidget *get_radio_button_in_scroll_group(GtkWidget *display_manager_dialog,
					    int imol_this) {

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
   if (display_manager_dialog) {
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
	       w = lookup_widget(display_manager_dialog, test_name.c_str());
	       if (w)
		  break;
	    }
	 }
      }
   }
   return w;
}


void
on_display_control_delete_molecule_button_clicked   (GtkButton       *button,
						   gpointer         user_data)
{

   int imol = GPOINTER_TO_INT(user_data);
   if (0)
      std::cout << "DEBUG:: calling close_molecule() for " << imol << " from "
		<< "on_display_control_delete_molecule_button_clicked"
		<< std::endl;
   close_molecule(imol);
}

#include "single-map-properties-dialog.hh"

void
on_display_control_map_properties_button_clicked   (GtkButton       *button,
						   gpointer         user_data)
{

/* Remove (comment out) archaic use of casting int * for user data. */
  int imol = GPOINTER_TO_INT(user_data);
  // GtkWidget *window = create_single_map_properties_dialog();
//   GtkWidget *patch_frame = lookup_widget(window,
// 					 "single_map_colour_button_frame");


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

  gtk_widget_show(dialog);
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
  widget_thing = lookup_widget(GTK_WIDGET(patch_frame), "single_map_colour_hbox");
  widget_thing = gtk_window_new(GTK_WINDOW_TOPLEVEL);


  printf("adding widget to patch_frame\n");
  gtk_container_add(GTK_CONTAINER(widget_thing), widget);

  // printf("gdk_gc_new\n");
  GdkWindow *window = 0;
#if (GTK_MAJOR_VERSION < 4)
  window = gtk_widget_get_window(widget);
#endif

  if (! window)
     return;

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
on_display_control_mol_displayed_button_toggled(GtkToggleButton *button,
                                                gpointer         user_data)
{
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "imol"));
   GtkWidget *active_toggle_button;

   // char *widget_name = (char *) malloc(100);
   // char *tmp_name;

/*   printf("DEBUG::  display toggle of molecule: %d\n", imol); */

   // strcpy(widget_name, "active_button_");
   // tmp_name = widget_name + strlen(widget_name);
   // snprintf(tmp_name, 4, "%-d", imol);
   
   int nn = 3;
   if (imol > 9) nn = 2;
   if (imol > 99) nn = 1;
   std::string four_char_imol(nn, '0');
   four_char_imol += coot::util::int_to_string(imol);
   std::string widget_name = "active_mol_button_" + four_char_imol;

/*   printf("mol display button clicked %d, active: %d\n", *imol, button->active); */

  if (imol >= 0 && imol < graphics_n_molecules()) {
     if (gtk_toggle_button_get_active(button))
        set_mol_displayed(imol, 1);
     else
        set_mol_displayed(imol, 0);

     /*     printf("looking up widget name %s\n", widget_name); */
     active_toggle_button = lookup_widget(GTK_WIDGET(button), widget_name.c_str());
     if (active_toggle_button) {
      /*  printf("INFO:: Got active_toggle_button from name: %s\n", widget_name); */

      if (mol_is_displayed(imol)) {
	// activate the button
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(active_toggle_button), TRUE);
      } else {
	/* deactivate the button */
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(active_toggle_button), FALSE);
      }


     } else {
        std::cout << "ERROR:: Failed to find active_toggle_button from name:"
                  << widget_name << std::endl;
     }
  } else {
     printf("ERROR:: (ignoring) display toggle of bogus molecule: %d\n", imol);
  }

}


void
on_display_control_mol_active_button_toggled   (GtkToggleButton  *toggle_button,
						gpointer         user_data)
{
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(toggle_button), "imol"));
   int iactive;
   if (gtk_toggle_button_get_active(toggle_button)) {
      set_mol_active(imol, 1);
/*     iactive = toggle_active_mol(*imol);  */
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
  GtkWidget *phs_cell_radiobutton_1;
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
  GtkWidget *phs_cell_none_radiobutton;
  GtkWidget *hbox32;
  GtkWidget *phs_cell_choice_ok_button;
  GtkWidget *phs_cell_choice_cancel_button;

/* messing about with string variables */
  gchar *widget_name;
  gchar *tmp_name;

  widget_name = (gchar *) malloc(100);

  phs_cell_chooser_vbox = lookup_widget(phs_cell_choice_window, "phs_cell_chooser_vbox");

  hbox33 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox33);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "hbox33", hbox33,
			  NULL);
  gtk_widget_show (hbox33);
  gtk_box_pack_start (GTK_BOX (phs_cell_chooser_vbox), hbox33, TRUE, TRUE, 4);
  gtk_container_set_border_width (GTK_CONTAINER (hbox33), 6);

  strcpy(widget_name, "phs_cell_radiobutton_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_radiobutton_1 = gtk_radio_button_new_with_label (phs_cell_group, "");
  phs_cell_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (phs_cell_radiobutton_1));
  // gtk_widget_ref (phs_cell_radiobutton_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_radiobutton_1,
			  NULL);
  gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (phs_cell_radiobutton_1), FALSE);

  gtk_widget_show (phs_cell_radiobutton_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_radiobutton_1, FALSE, FALSE, 4);

  label53 = gtk_label_new (_("Symm"));
  // gtk_widget_ref (label53);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label53", label53,
			  NULL);
  gtk_widget_show (label53);
  gtk_box_pack_start (GTK_BOX (hbox33), label53, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_symm_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_symm_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_symm_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_symm_entry_1,
			  NULL);
  gtk_widget_show (phs_cell_symm_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_symm_entry_1, TRUE, TRUE, 0);

  gtk_widget_set_size_request(GTK_WIDGET(phs_cell_symm_entry_1), 80, -2);

  label54 = gtk_label_new (_("a"));
  // gtk_widget_ref (label54);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label54", label54,
			  NULL);
  gtk_widget_show (label54);
  gtk_box_pack_start (GTK_BOX (hbox33), label54, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_a_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_a_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_a_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_a_entry_1, NULL);

  gtk_widget_show (phs_cell_a_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_a_entry_1, TRUE, TRUE, 0);
  gtk_widget_set_size_request (phs_cell_a_entry_1, 65, -2);

  label55 = gtk_label_new (_("b"));
  // gtk_widget_ref (label55);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label55", label55,
			  NULL);
  gtk_widget_show (label55);
  gtk_box_pack_start (GTK_BOX (hbox33), label55, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_b_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_b_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_b_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_b_entry_1, NULL);
  gtk_widget_show (phs_cell_b_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_b_entry_1, TRUE, TRUE, 0);
  gtk_widget_set_size_request (phs_cell_b_entry_1, 65, -2);

  label56 = gtk_label_new (_("c"));
  // gtk_widget_ref (label56);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label56", label56,
			  NULL);
  gtk_widget_show (label56);
  gtk_box_pack_start (GTK_BOX (hbox33), label56, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_c_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_c_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_c_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_c_entry_1, NULL);
  gtk_widget_show (phs_cell_c_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_c_entry_1, TRUE, TRUE, 0);
  gtk_widget_set_size_request (phs_cell_c_entry_1, 65, -2);

  label57 = gtk_label_new (_("alpha"));
  // gtk_widget_ref (label57);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label57", label57,
			  NULL);
  gtk_widget_show (label57);
  gtk_box_pack_start (GTK_BOX (hbox33), label57, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_alpha_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_alpha_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_alpha_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_alpha_entry_1,
			  NULL);
  gtk_widget_show (phs_cell_alpha_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_alpha_entry_1, TRUE, TRUE, 0);
  gtk_widget_set_size_request (phs_cell_alpha_entry_1, 60, -2);

  label58 = gtk_label_new (_("beta"));
  // gtk_widget_ref (label58);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label58", label58,
			  NULL);
  gtk_widget_show (label58);
  gtk_box_pack_start (GTK_BOX (hbox33), label58, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_beta_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_beta_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_beta_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_beta_entry_1, NULL);
  gtk_widget_show (phs_cell_beta_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_beta_entry_1, TRUE, TRUE, 0);
  gtk_widget_set_size_request (phs_cell_beta_entry_1, 60, -2);

  label59 = gtk_label_new (_("gamma"));
  // gtk_widget_ref (label59);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "label59", label59,
			  NULL);
  gtk_widget_show (label59);
  gtk_box_pack_start (GTK_BOX (hbox33), label59, FALSE, FALSE, 2);


  strcpy(widget_name, "phs_cell_gamma_entry_");
  tmp_name = widget_name + strlen(widget_name);
  snprintf(tmp_name, 4, "%-d", n);


  phs_cell_gamma_entry_1 = gtk_entry_new ();
  // gtk_widget_ref (phs_cell_gamma_entry_1);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), widget_name,
			  phs_cell_gamma_entry_1, NULL);
  gtk_widget_show (phs_cell_gamma_entry_1);
  gtk_box_pack_start (GTK_BOX (hbox33), phs_cell_gamma_entry_1, TRUE, TRUE, 0);
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
  GtkWidget *phs_cell_none_radiobutton;


  phs_cell_chooser_vbox = lookup_widget(phs_cell_choice_window, "phs_cell_chooser_vbox");

  hbox34 = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
  // gtk_widget_ref (hbox34);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "hbox34", hbox34, NULL);
  gtk_widget_show (hbox34);
  gtk_box_pack_start (GTK_BOX (phs_cell_chooser_vbox), hbox34, TRUE, TRUE, 0);
  gtk_container_set_border_width (GTK_CONTAINER (hbox34), 6);

  phs_cell_none_radiobutton = gtk_radio_button_new_with_label (phs_cell_group, _("None of the Above"));
  phs_cell_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (phs_cell_none_radiobutton));
  // gtk_widget_ref (phs_cell_none_radiobutton);
  g_object_set_data_full (G_OBJECT (phs_cell_choice_window), "phs_cell_none_radiobutton",
			  phs_cell_none_radiobutton,
			  NULL);

  gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (phs_cell_none_radiobutton), TRUE);

  gtk_widget_show (phs_cell_none_radiobutton);
  gtk_box_pack_start (GTK_BOX (hbox34), phs_cell_none_radiobutton, FALSE, FALSE, 4);

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
      GtkWidget *t = lookup_widget(display_control_window_glade, name.c_str());
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
      GtkWidget *t = lookup_widget(display_control_window_glade, name.c_str());
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
      GtkWidget *t = lookup_widget(display_control_window_glade, name.c_str());
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
      gtk_box_pack_start(GTK_BOX(add_rep_vbox), hbox, FALSE, FALSE, 0);
      std::string label = name;
      GtkWidget *toggle_button_show_it = gtk_check_button_new_with_label(label.c_str());
      // gtk_widget_ref (toggle_button_show_it);
      if (show_it) {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button_show_it), TRUE);
	 gtk_widget_show(add_rep_all_on_checkbutton);
      } else {
	 gtk_widget_hide(add_rep_all_on_checkbutton);
      }
      int cc = encode_ints(imol_no, add_rep_no);
      g_signal_connect(G_OBJECT(toggle_button_show_it), "toggled",
			 G_CALLBACK(add_rep_toggle_button_toggled),
			 GINT_TO_POINTER(cc));
      gtk_box_pack_start(GTK_BOX(hbox), toggle_button_show_it, FALSE, FALSE, 0);

      gtk_widget_show(toggle_button_show_it);
      gtk_widget_show(hbox);
      gtk_widget_show(add_rep_vbox);
      gtk_widget_show(add_rep_frame);
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

GtkWidget*
create_splash_screen_window_for_file(const char *file_name)
{
  GtkWidget *splash_screen_window;

  splash_screen_window = gtk_window_new (GTK_WINDOW_POPUP);
  gtk_widget_set_name (splash_screen_window, "splash_screen_window");
  gtk_window_set_title (GTK_WINDOW (splash_screen_window), _("Coot"));
  gtk_window_set_position (GTK_WINDOW (splash_screen_window), GTK_WIN_POS_CENTER);

  gtk_window_set_type_hint (GTK_WINDOW (splash_screen_window),
			    GDK_WINDOW_TYPE_HINT_SPLASHSCREEN);

  GtkWidget *image6807 = create_pixmap(splash_screen_window, file_name);
  gtk_widget_show (image6807);

  gtk_container_add (GTK_CONTAINER (splash_screen_window), image6807);

  return splash_screen_window;
}


#ifdef COOT_USE_GTK2_INTERFACE
/* Hack in a function, because it's missing somehow from Bernhard's
   commits */

/* I dont think this should be here at all!
   I think it should be in gtk2-interface.c !
   But somehow it disappeared from there, so I will put it back!!!
   Where did all my things go........ */

#endif	/* COOT_USE_GTK2_INTERFACE */
