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

#include "widget-from-builder.hh"

void add_on_validation_graph_mol_options(GtkWidget *menu, const char *type_in);


// from support.h
// GtkWidget* lookup_widget (GtkWidget *widget, const gchar *widget_name);
#include "support.h"

// this from callbacks.h (which I don't want to include here)
typedef const char entry_char_type;


#if 0 // testing
extern "C" G_MODULE_EXPORT
void
on_calculate_load_tutorial_model_and_data1_activate(GtkMenuItem *menuitem, gpointer user_data) {
   std::cout << "Load tutorial model and data" << std::endl;
   load_tutorial_model_and_data();
}
#endif


#if 0 // testing
// testing function for 2021 main window
extern "C" G_MODULE_EXPORT
void on___glade_unnamed_94_activate(GtkMenuItem *menuitem, gpointer user_data) {
   std::cout << "Load tutorial model and data" << std::endl;
   load_tutorial_model_and_data();
}
#endif


extern "C" G_MODULE_EXPORT
void
on_window1_destroy (GtkWidget       *object,
                                        gpointer         user_data) {

   // 20220327-PE the window manager killing the Coot window doesn't call this

   std::cout << "on_window1_destroy()" << std::endl;

   // this no longer exists
   // gtk_main_quit();
}

/* When the user uses the window manager to close coot, this gets called. */
/* When the window manager "close window" events happens it send the
   application a delete_event event.  If we return FALSE from this
   attached function, then a "destroy" signal will be emitted.
   Returning TRUE means that we don't want the window to be destroyed.
   See helloworld.c in the gtk_tut.txt */
extern "C" G_MODULE_EXPORT
gboolean
on_window1_delete_event                (GtkWidget       *widget,
                                                            GdkEvent        *event,
                                                            gpointer         user_data)
{
   printf("---------------------------- on_window1_delete_event() ------------------\n");
  /* coot_checked_exit() calls coot_real_exit() and that calls exit(),
     so we don't (normally?) return from this function. */
  coot_no_state_real_exit(0);
  return 0;
}


extern "C" G_MODULE_EXPORT
void on_file1_activate (GMenuItem     *menuitem,
                                            gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_exit1_activate (GMenuItem     *menuitem,
                                       gpointer         user_data)
{
/*   printf("calling gtk_main_quit()\n"); */
/*   gtk_main_quit(); */

/*    printf("---------------------------- on_exit1_activate() ------------------\n"); */
   coot_checked_exit(0); /* without error */
}


extern "C" G_MODULE_EXPORT
void
on_clipping1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_clipping1_activate();
}


extern "C" G_MODULE_EXPORT
void
on_cancel_coords_button1_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *coords_fileselection1 = widget_from_builder( "coords_fileselection1");
   gtk_widget_hide(coords_fileselection1);
}


extern "C" G_MODULE_EXPORT
void
on_cancel_dataset_button1_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dataset_fileselection1 = widget_from_builder("dataset_fileselection1");
   gtk_widget_hide(dataset_fileselection1);
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
   gtk_widget_hide(widget_from_builder("column_label_window"));
}


extern "C" G_MODULE_EXPORT
void
on_clipping_button_clicked             (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_hide(widget_from_builder( "clipping_window"));
}




/* Density Radius Window Widgets  */

extern "C" G_MODULE_EXPORT
void
on_density_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkEntry      *entry;

/*    entry_char_type *text; */
   const char *text;

   GtkEntry *entry_xray = GTK_ENTRY(widget_from_builder("map_parameters_xray_radius_entry"));
   GtkEntry *entry_em   = GTK_ENTRY(widget_from_builder("map_paameters_em_radius_entry"));
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
   gtk_widget_hide(dialog);

}

void
show_map_parameters_dialog() {

   char *text;
   int imol = 0;		/* FIXME */

   // this widget is looked up in
   // on_density_ok_button_clicked()

   GtkWidget *density_window = widget_from_builder("global_map_properties_window");

   // 20220315-PE archaic but OK for now
   GtkEntry *entry_xray = GTK_ENTRY(widget_from_builder("map_parameters_x_ray_radius_entry"));
   GtkEntry *entry_em   = GTK_ENTRY(widget_from_builder("map_parameters_em_radius_entry"));
   text = get_text_for_density_size_widget(); /* const gchar *text */
   gtk_editable_set_text(GTK_EDITABLE(entry_xray), text);
   text = get_text_for_density_size_em_widget(); /* const gchar *text */
   gtk_editable_set_text(GTK_EDITABLE(entry_em), text);
   free (text);
   text = 0;

   /* Now the iso level increment entry  */

   GtkWidget *entry;
   entry = widget_from_builder("iso_level_increment_entry");
   text = get_text_for_iso_level_increment_entry(imol);

   gtk_editable_set_text(GTK_EDITABLE(GTK_ENTRY(entry)), text);

   /* Now the iso level for the differenece map increment entry  */

   entry = widget_from_builder("diff_map_iso_level_increment_entry");
   text = get_text_for_diff_map_iso_level_increment_entry(imol);

   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   /* Now the map rate multiplier: */
   entry = widget_from_builder("map_sampling_rate_entry");
   text = get_text_for_map_sampling_rate_text();

   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   GtkWidget *checkbutton = widget_from_builder("map_dynamic_map_sampling_checkbutton");
   set_map_dynamic_map_sampling_checkbutton(checkbutton);
   checkbutton = widget_from_builder("map_dynamic_map_size_display_checkbutton");
   set_map_dynamic_map_display_size_checkbutton(checkbutton);

 /* Show the widget */
   gtk_widget_show(density_window);
}


/* In the menubar, Edit Density size has been selected. */
// this is now map_parameters "Map Parameters"
extern "C" G_MODULE_EXPORT
void
on_density_size1_activate (GMenuItem     *menuitem,
                           gpointer         user_data) {

   show_map_parameters_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_density_cancel_clicked              (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_hide(widget_from_builder("global_map_properties_window"));
}



extern "C" G_MODULE_EXPORT
gboolean
on_hscale1_button_press_event(GtkWidget  *widget,
                              GtkButton  *event,
                              gpointer    user_data) {

   std::cout << "debug:: on_hscale1_button_press_event - what even is this!?" << std::endl;
   return FALSE;
}


extern "C" G_MODULE_EXPORT
gboolean
on_hscale1_button_release_event(GtkWidget *widget,
                                GtkButton *event,
                                gpointer   user_data) {

   std::cout << "debug:: on_hscale1_button_release_event - what even is this!?" << std::endl;
   return FALSE;
}


extern "C" G_MODULE_EXPORT
void
on_fps1_activate(GMenuItem     *menuitem,
                 gpointer         user_data) {

   // GtkWidget *fps_window = create_fps_window();
   GtkWidget *fps_window = widget_from_builder("fps_window");
   GtkButton       *button;

   if ( get_fps_flag() == 1) {
      button = GTK_BUTTON(widget_from_builder(
					"radiobutton1"));
   } else {
      button = GTK_BUTTON(widget_from_builder(
					"radiobutton2"));
   }

   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
   gtk_widget_show(fps_window);

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

   gtk_widget_hide(widget_from_builder("fps_window"));

}


extern "C" G_MODULE_EXPORT
void
on_active_map_ok_button_clicked(GtkButton       *button,
                                gpointer         user_data) {

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("active_map_radiobutton_yes"))))
      set_active_map_drag_flag(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("active_map_radiobutton_no"))))
      set_active_map_drag_flag(0);

   gtk_widget_hide(widget_from_builder( "active_map_window"));

}


extern "C" G_MODULE_EXPORT
void
on_dragged_map1_activate (GMenuItem     *menuitem,
                                              gpointer         user_data)
{
   // GtkWidget *active_map_window = create_active_map_window();
   GtkWidget *active_map_window = widget_from_builder("active_map_window");
   GtkButton       *button;

   if ( get_active_map_drag_flag() == 1 ) {
      button = GTK_BUTTON(widget_from_builder("active_map_radiobutton_yes"));
   } else {
      button = GTK_BUTTON(widget_from_builder("active_map_radiobutton_no"));
   }

   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
   gtk_widget_show(active_map_window);
}

// Draw -> Map Colour has been removed.
//
// extern "C" G_MODULE_EXPORT
// void
// on_map_colour1_activate (GMenuItem     *menuitem,
//                                              gpointer         user_data)
// {
//    // GtkWidget *menu = widget_from_builder("rotamer_analysis1");

//    //std::cout << "::::::::::::::::::::::::: on_map_colour1_activate() " << std::endl;

//    GtkWidget *menu = widget_from_builder("map_colour1");
//    if (menu) {
//       add_on_map_colour_choices(menu);
//    } else {
//       printf("ERROR:: failed to get map_colour1 menu in on_map_colour1_activate\n");
//    }
// }

extern "C" G_MODULE_EXPORT
void
on_show_symmetry1_activate             (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *show_symm_window;

   show_symm_window = wrapped_create_show_symmetry_window();

/* Now show the popup eventually */

   gtk_widget_show(show_symm_window);

}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_ok_button_clicked     (GtkButton       *button,
                                                            gpointer         user_data) {
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

/* Goodbye Mr Widget */
   gtk_widget_hide(widget_from_builder( "show_symmetry_window"));

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

/*    gtk_widget_show(colorseldlg);  */
}



extern "C" G_MODULE_EXPORT
void
on_about_ok_button_clicked             (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_hide(widget_from_builder("about_window"));
}


extern "C" G_MODULE_EXPORT
void
on_anisotropic_atoms1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *aniso_window;
   GtkWidget *entry;
   char *text;
   GtkButton       *button;
   GtkScale *hscale;
   GtkAdjustment *adjustment;
   float hscale_initial;

   // aniso_window = create_aniso_window();
   aniso_window = widget_from_builder("aniso_window");

/* Show Aniso Radiobuttons */
   if (get_show_aniso() == 1) {
      button = GTK_BUTTON(widget_from_builder(
					"show_aniso_yes_radiobutton"));
   } else {
      button = GTK_BUTTON(widget_from_builder(
					"show_aniso_no_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


/* Set Limit Radiobuttons */
   if (get_show_limit_aniso() == 1) {
      button = GTK_BUTTON(widget_from_builder(
					"limit_display_radius_yes_radiobutton"));
   } else {
      button = GTK_BUTTON(widget_from_builder(
					"limit_display_radius_no_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


/* The Aniso Display Limit Entry */

   entry = widget_from_builder("limit_display_radius_entry");

   text = get_text_for_aniso_limit_radius_entry();
   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   free(text);

/* The Probability Hscale */

   hscale = GTK_SCALE(widget_from_builder(
				    "aniso_probability_hscale"));

   hscale_initial = get_aniso_probability();

   adjustment = GTK_ADJUSTMENT
         (gtk_adjustment_new(hscale_initial, 0.0, 110.0, 0.01, 4.0, 10.));

   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);

#if 0 // FIXME-LATER-PE header issues
   g_signal_connect (G_OBJECT(adjustment), "value_changed",
		       G_CALLBACK(aniso_probability_adjustment_changed),
		       NULL);
#endif


/* And finally show the widget */
   gtk_widget_show(aniso_window);

}


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
   gtk_widget_hide(widget_from_builder("aniso_window"));
}

extern "C" G_MODULE_EXPORT
void aniso_probability_adjustment_changed(GtkAdjustment *adj, GtkWidget *window) {

  set_aniso_probability(gtk_adjustment_get_value(adj));
}


extern "C" G_MODULE_EXPORT
void on_smooth_scrolling_window_ok_button_clicked (GtkButton       *button,
                                                                       gpointer         user_data)
{
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

   gtk_widget_hide(widget_from_builder("smooth_scroll_window"));

}


extern "C" G_MODULE_EXPORT
void
on_recentring1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *smooth_scroll_window;
   GtkWidget *entry;
   char *text;
   GtkButton *button;

   // smooth_scroll_window = create_smooth_scroll_window();
   smooth_scroll_window = widget_from_builder("smooth_scroll_window"); // what is this?

/* Show Smooth Scrolling Radio Buttons */
   if (get_smooth_scroll() == 1) {

      button = GTK_BUTTON(widget_from_builder(
					"smooth_scroll_yes_radiobutton"));
   } else {
      button = GTK_BUTTON(widget_from_builder(
					"smooth_scroll_no_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


/* Smooth Scroll Distance Limit */

   entry = widget_from_builder("smooth_scroll_limit_entry");

   text = get_text_for_smooth_scroll_limit();
   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   free(text);

/*  Smooth Scroll Steps */


   entry = widget_from_builder("smooth_scroll_steps_entry");

   text = get_text_for_smooth_scroll_steps();

   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   free(text);


/* And finally show the widget */
   gtk_widget_show(smooth_scroll_window);

}


extern "C" G_MODULE_EXPORT
void
on_font_size1_activate                 (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *font_size_window;

   GtkButton *button = 0;

   // font_size_window = create_font_size_window();
   font_size_window = widget_from_builder("font_size_window"); // what is this?

/* The Font Size RadioButtons */

   if (get_font_size() == 1) {
            button = GTK_BUTTON(widget_from_builder(
					"font_size_small_radiobutton"));
   }
   if (get_font_size() == 2) {
            button = GTK_BUTTON(widget_from_builder(
					"font_size_medium_radiobutton"));
   }
   if (get_font_size() == 3) {
            button = GTK_BUTTON(widget_from_builder(
					"font_size_large_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

/* show widget */
   gtk_widget_show(font_size_window);
}

extern "C" G_MODULE_EXPORT
void
on_font_size_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{

/* The Font Size RadioButtons */
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("font_size_small_radiobutton"))))

      set_font_size(1);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("font_size_medium_radiobutton"))))

      set_font_size(2);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget_from_builder("font_size_large_radiobutton"))))

      set_font_size(3);

/* goodbye widget */
   gtk_widget_hide(widget_from_builder("font_size_window"));

}


/* not used? */
extern "C" G_MODULE_EXPORT
void
on_greer_skeleton1_activate            (GMenuItem     *menuitem,
                                        gpointer         user_data)
{


}

/* not used? */
extern "C" G_MODULE_EXPORT
void
on_cowtan_foadi_skeleton1_activate     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_these_are_placeholders1_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_mapcolourmap1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{


}


extern "C" G_MODULE_EXPORT
void
on_map1_activate                       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   printf("ping the colourchange widget for map1\n");
}


extern "C" G_MODULE_EXPORT
void
on_map2_activate                       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

extern "C" G_MODULE_EXPORT
void
on_attach_scroll_wheel_to_which_map_1_activate (GMenuItem     *menuitem,
						gpointer         user_data) {

#if 0 // 20220602-PE I don't think that this widget exists any more

   /* it doesn't matter what we pass back, a lookup is done to find the
      submenu */
   GtkWidget *menu = widget_from_builder("attach_scroll_wheel_to_which_map_1");
   if (menu) {
      add_on_map_scroll_wheel_choices(menu);
   } else {
      printf("ERROR:: failed to get menu in on_attach_scroll_wheel_to_which_map_1_activate\n");
   }
#endif

}


extern "C" G_MODULE_EXPORT
void
on_greer_on_activate                   (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   printf("Making skel_greer....\n");
   skel_greer_on();

}


extern "C" G_MODULE_EXPORT
void
on_greer_off_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   printf("skel_greer OFF\n");
   skel_greer_off();
}

/* Notice, with this sort of button, for some reason, when we select
   "on" we run on_foadi_on_activate() then on_foadi_off_activate().

   When we select on, we run on_foadi_off_activate()
   on_foadi_on_activate().

*/
extern "C" G_MODULE_EXPORT
void
on_foadi_on_activate                   (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_foadi_off_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
/*    printf("skel_foadi OFF\n");  */

/*   skel_foadi_off(); dead */

}


extern "C" G_MODULE_EXPORT
void
on_open_map1_activate (GMenuItem     *menuitem,
                                           gpointer         user_data) {

   GtkWidget *map_name_chooser = widget_from_builder("map_name_filechooser_dialog");
   GtkWidget *filter_button = add_filename_filter_button(map_name_chooser, COOT_MAP_FILE_SELECTION);

   gtk_widget_show (map_name_chooser);
   add_is_difference_map_checkbutton(map_name_chooser);
   set_directory_for_filechooser(map_name_chooser);

}



extern "C" G_MODULE_EXPORT
void
on_cancel_button_map_name_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{

   gtk_widget_hide(widget_from_builder("map_name_fileselection1"));

}


extern "C" G_MODULE_EXPORT
void
on_vt_flat1_activate                   (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   vt_surface(VT_FLAT);

}


extern "C" G_MODULE_EXPORT
void
on_vt_spherical_surface1_activate      (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   vt_surface(VT_SPHERICAL);


}


extern "C" G_MODULE_EXPORT
void
on_skeleton_colour1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *col_sel_window = create_skeleton_colour_selection_window();
   GtkWidget *col_sel_window = widget_from_builder("skeleton_colour_selection_window");
   gtk_widget_show(col_sel_window);
}



extern "C" G_MODULE_EXPORT
void
on_phs_info_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   phs_pdb_cell_symm(); /* which runs
			   create_phs_coordinates_fileselection() */

   gtk_widget_hide(widget_from_builder( "phs_info_box"));

}


extern "C" G_MODULE_EXPORT
void
on_phs_info_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_hide(widget_from_builder("phs_info_box"));

}


extern "C" G_MODULE_EXPORT
void
on_cancel_phs_coord_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{

  gtk_widget_hide(widget_from_builder( "phs_coordinates_fileselection"));


}


extern "C" G_MODULE_EXPORT
void
on_map_and_mol_control1_activate       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
}

extern "C" G_MODULE_EXPORT
void
on_draw_display_only_active_activate (GMenuItem     *menuitem,
                                                          gpointer         user_data)
{

  /* display only the active mol */
  display_only_active();

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
  gtk_widget_hide(widget);	/* There is something that had been
				   added to the Go To Atom window that
				   is not a widget.  It has been
				   cleared perhaps but not destroyed.
				   I don't think that it is the dialog
				   itself.  The problem does not
				   happen in the GTK1 path.  */
}


extern "C" G_MODULE_EXPORT
void
on_go_to_atom1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data) {
     GtkWidget *widget;

     /* fills the window too */
     widget = wrapped_create_goto_atom_window(); // uses gtkbuilder

				/* now we need to fill the entry boxes
				   with default vaules and the option
				   menu according to molecules that
				   have coordinates. */


     /* and finally show the widget: */
     gtk_widget_show(widget);
}

extern "C" G_MODULE_EXPORT
void on_go_to_atom_button_clicked (GtkButton       *button,
                                   gpointer         user_data) {

   GtkWidget *dialog = wrapped_create_goto_atom_window(); // uses gtkbuilder
   gtk_widget_show(dialog);
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
  goto_next_atom_maybe_new(window);
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
  goto_previous_atom_maybe_new(window);
}


// -----------------Toolbar buttons - Go to Atom and Go to ligand


extern "C" G_MODULE_EXPORT
void
on_auto_clear_atom_pull_restraints_togglebutton_toggled(GtkToggleButton *toggle_button,
                                                        gpointer user_data) {

   std::cout << "on_auto_clear_atom_pull_restraints_togglebutton_toggled()" << std::endl;

}


extern "C" G_MODULE_EXPORT
void
on_clear_atom_pull_restraints_toolbutton_clicked(GtkButton   *toolbutton,
                                                 gpointer         user_data) {

   std::cout << "on_clear_atom_pull_restraints_toolbutton_clicked()" << std::endl;
}

// -----------------------------------------

extern "C" G_MODULE_EXPORT
void
on_skeleton_box_radius1_activate       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *entry;
  gchar *text;

  // GtkWidget *widget = create_skeletonization_box_radius_window();
  GtkWidget *widget = widget_from_builder("skeletonization_box_radius_window");

  text = get_text_for_skeleton_box_size_entry();
  entry = widget_from_builder("skeleton_box_size_entry");
  gtk_editable_set_text(GTK_EDITABLE(entry), text);
  g_free(text);

	/* show the widget */
  gtk_widget_show(widget);

}


extern "C" G_MODULE_EXPORT
void
on_skeletonization_level1_activate     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *entry;
   gchar *text;

   // GtkWidget *widget = create_skeletonization_level_window();
   GtkWidget *widget = widget_from_builder("skeletonization_level_window");

   /* Fill the entry: */
   text = get_text_for_skeletonization_level_entry();
   entry = widget_from_builder("skeleton_level_entry");
   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   g_free(text);
   /* show the widget */
   gtk_widget_show(widget);
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
   gtk_widget_hide(widget_from_builder("skeletonization_box_radius_window"));
}


extern "C" G_MODULE_EXPORT
void
on_skel_box_radius_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{

   gtk_widget_hide(widget_from_builder("skeletonization_box_radius_window"));

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

  gtk_widget_hide(widget_from_builder("skeletonization_level_window"));

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

  gtk_widget_hide(widget_from_builder("skeletonization_level_window"));
}


extern "C" G_MODULE_EXPORT
void
on_autobuild_ca_on_activate            (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

/*   autobuild_ca_on();  - removed to junk now */

}


extern "C" G_MODULE_EXPORT
void
on_autobuild_ca_off_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
    autobuild_ca_off();
}


extern "C" G_MODULE_EXPORT
void
on_display_control_ok_button_clicked   (GtkButton       *button,
                                                            gpointer         user_data) {

   // GtkWidget *w = widget_from_builder("display_control_window_glade");
   GtkWidget *w = widget_from_builder("display_control_window_glade");
   GtkWidget *maps_vbox = 0;
   GtkWidget *molecules_vbox = 0;
   GtkWidget *pane = 0;

   reset_graphics_display_control_window(); /* Needed! (also resets the scroll group) */
   if (w) {
      store_window_position(COOT_DISPLAY_CONTROL_WINDOW, w);
      maps_vbox      = widget_from_builder("display_map_vbox");
      molecules_vbox = widget_from_builder("display_molecule_vbox");
      pane           = widget_from_builder("display_control_vpaned");

      /* store the size, actually */
      if (maps_vbox)
         store_window_position(COOT_DISPLAY_CONTROL_MAPS_VBOX, maps_vbox);
      if (molecules_vbox)
         store_window_position(COOT_DISPLAY_CONTROL_MOLECULES_VBOX, molecules_vbox);
      if (pane)
         store_window_position(COOT_DISPLAY_CONTROL_PANE, pane);

      // gtk_widget_destroy(w); // 20220309-PE not these days, buddy-boy
      gtk_widget_hide(w);
   } else {
      printf("Error:: in on_display_control_ok_button_clicked() failed to lookup display_control_window_glade\n");
   }
}


extern "C" G_MODULE_EXPORT
void
on_display_control_window_glade_destroy (GtkWidget       *object,
                                                             gpointer         user_data) {

   // Do nothing the the poor old widget.

   //    std::cout << "-------------------------------------------- on_display_control_window_glade_destroy()"
   // << std::endl;

   // we don't want to use the stored pointer in graphics_info_t any more. We just use
   // widget_from_builder() to get the display control dialog.

   // reset_graphics_display_control_window(); /* (also resets the scroll group) */
}


extern "C" G_MODULE_EXPORT
gboolean
on_display_control_window_glade_delete_event(GtkWidget       *widget,
                                                                 GdkEvent        *event,
                                                                 gpointer         user_data) {

   if (false)
      std::cout << "------------------ on_display_control_window_glade_delete_event()"
                << std::endl;

   gtk_widget_hide(widget);
   return TRUE;
}


extern "C" G_MODULE_EXPORT
void
on_rotation_centre_size_ok_button_clicked(GtkButton       *button,
                                                              gpointer         user_data) {

				/* pass back the value from the entry */
  GtkEntry *entry = GTK_ENTRY(widget_from_builder("rotation_centre_cube_size_entry"));
  const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  set_rotation_centre_size_from_widget(text);
  gtk_widget_hide(widget_from_builder("rotation_centre_cube_size_window"));

}


extern "C" G_MODULE_EXPORT
void
on_pink_pointer_size1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   // Will this ever get used?

   GtkWidget *window;
   GtkWidget *entry;
   gchar *text;

   // window = create_rotation_centre_cube_size_window();
   window = widget_from_builder("rotation_centre_cube_size_window");
   entry = widget_from_builder("rotation_centre_cube_size_entry");
   text = get_text_for_rotation_centre_cube_size();
   gtk_editable_set_text(GTK_EDITABLE(entry), text);
   gtk_widget_show(window);
   free(text);

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
   gtk_widget_hide(window);
}

extern "C" G_MODULE_EXPORT
void
on_test_thing1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  test_fragment();
}


/* void */
/* on_regularize1_activate                (GMenuItem     *menuitem, */
/*                                         gpointer         user_data) */
/* { */
/*   do_regularize(); */
/* } */


extern "C" G_MODULE_EXPORT
void
on_scripting_window_activate(GMenuItem     *menuitem,
                                                 gpointer         user_data) {
}


extern "C" G_MODULE_EXPORT
void
on_scheme_window_close_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   // 20220309-PE when will this get used...?
   //  GtkWidget *window = widget_from_builder("scheme_window");
  GtkWidget *window = widget_from_builder("scheme_window");
  gtk_widget_hide(window);

}


extern "C" G_MODULE_EXPORT
void
on_fetch_pdb_using_code1_activate (GMenuItem     *menuitem,
                                                       gpointer         user_data)
{

   int n = COOT_ACCESSION_CODE_WINDOW_OCA;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}


extern "C" G_MODULE_EXPORT
void
on_fetch_pdb_and_sf_using_code1_activate (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   int n = COOT_ACCESSION_CODE_WINDOW_OCA_WITH_SF;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}


extern "C" G_MODULE_EXPORT
void
on_fetch_pdb_and_map_using_eds1_activate (GMenuItem     *menuitem,
                                                              gpointer         user_data)
{
   int n = COOT_ACCESSION_CODE_WINDOW_EDS;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}


extern "C" G_MODULE_EXPORT
void
on_fetch_pdb_and_map_using_pdb_redo1_activate
                                        (GMenuItem     *menuitem,
					 gpointer         user_data) {
   int n = COOT_ACCESSION_CODE_WINDOW_PDB_REDO;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);

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
on_fetch_alphafold_model_using_uniprot_id1_activate(GMenuItem     *menuitem,
                                                    gpointer         user_data) {
   int n = COOT_UNIPROT_ID;
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);

}



#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_accession_code_entry_key_press_event (GtkWidget       *widget,
                                         GdkEventKey     *event,
                                         gpointer         user_data) {

 /* go somewhere if keypress was a carriage return  */

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
     GtkWidget *entry = widget;
     GtkWidget *dialog = widget_from_builder("accession_code_window");
     handle_get_accession_code(dialog, entry);
  }
  return FALSE;
}
#endif

extern "C" G_MODULE_EXPORT
void
on_accession_code_get_it_button_clicked(GtkButton *button, gpointer user_data) {

   GtkWidget *entry = widget_from_builder("accession_code_entry");
   GtkWidget *frame = widget_from_builder("accession_code_frame");
   handle_get_accession_code(frame, entry);
}



extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  add_on_rama_choices();
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

   gtk_widget_hide(window);
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

   gtk_widget_hide(window);

   /* maybe we should also destroy the edit backbone torsion dialog, if it exists. */

}

extern "C" G_MODULE_EXPORT
void
on_dynarama_window_destroy             (GtkWidget       *object,
                                        gpointer         user_data)
{
  /* Maybe object is window? */
   GtkWidget *window = widget_from_builder(
				     "dynarama_window");
/* we can't use store_window_position here because object->window is
   null in the store_window_position() function, heyho. */
/* Hopefully the configure-events for this dialog will save the
   position. */
/*    store_window_position(COOT_RAMACHANDRAN_PLOT_WINDOW, GTK_WIDGET(object)); */
   int imol = get_mol_from_dynarama(window); // return -9999 on edit rama window
   if (imol >= 0) {
      set_dynarama_is_displayed(0, imol); // which frees/deletes the
			      	          // memory of the user data.
   }
}



extern "C" G_MODULE_EXPORT
void
on_background_black1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  set_background_colour(0.0,0.0,0.0);

}


extern "C" G_MODULE_EXPORT
void
on_background_white1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  set_background_colour(1.0,1.0,1.0);
}


extern "C" G_MODULE_EXPORT
void
on_find_ligands1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
				/* Not used any more */
   // see on_model_refine_dialog_find_ligands_button_clicked()
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                                            gpointer         user_data) {

   int n_ligands = execute_get_mols_ligand_search(GTK_WIDGET(button));
			                    	/* which then runs execute_ligand_search */
   if (n_ligands > 0) {
      GtkWidget *window = widget_from_builder("find_ligand_dialog");
      // free_ligand_search_user_data(GTK_WIDGET(button)); // not if not destroyed? Needs checking.
      gtk_widget_hide(window);
   } else {
      info_dialog("WARNING:: No ligands were selected");
   }
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = widget_from_builder("find_ligand_dialog");
   free_ligand_search_user_data(GTK_WIDGET(button));
   gtk_widget_hide(window);
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_many_atoms_continue_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   GtkWidget *window = widget_from_builder("find_ligand_many_atoms_dialog");
   GtkWidget *find_ligand_dialog = widget_from_builder("find_ligand_dialog");

   execute_ligand_search();
   gtk_widget_hide(window);
   gtk_widget_hide(find_ligand_dialog);
}


extern "C" G_MODULE_EXPORT
void
on_find_ligand_many_atoms_cancel_button_clicked(GtkButton       *button,
                                                                    gpointer         user_data) {

   GtkWidget *window = widget_from_builder("find_ligand_many_atoms_dialog");
   gtk_widget_hide(window);
}


extern "C" G_MODULE_EXPORT
void
on_prune_and_colour1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

  do_skeleton_prune();

}

extern "C" G_MODULE_EXPORT
void
on_residue_info1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_residue_info_dialog();
}

extern "C" G_MODULE_EXPORT
gboolean
on_residue_info_dialog_delete_event(GtkWidget       *widget,
                                                        GdkEvent        *event,
                                                        gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}


/* the model fit refine menu item was activated */
extern "C" G_MODULE_EXPORT
void
on_model_refine_activate               (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *window_or_frame;
   window_or_frame = wrapped_create_model_fit_refine_dialog();
   gtk_widget_show(window_or_frame);
}




extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_dismiss_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{

   GtkWidget *hbox   = widget_from_builder("model_fit_refine_dialog_vbox");
   GtkWidget *dialog = widget_from_builder("model_refine_dialog");

   clear_pending_picks();
   normal_cursor();
   dialog = close_model_fit_dialog(hbox);
   store_window_position(COOT_MODEL_REFINE_DIALOG, dialog);
   /* was it a top-level?  If so, kill off the top-level. */
   if (dialog)
      gtk_widget_hide(dialog);
}



extern "C" G_MODULE_EXPORT
void
on_save_coords_dialog_save_button_clicked (GtkButton       *button,
                                                               gpointer         user_data) {

   // we need to select the molecule to save - this is someone clicking on the
   // "Save Molecule" button in the save molecule chooser - not in a file selector

   GtkWidget *combobox = widget_from_builder("save_coordinates_combobox");
   GtkWidget *dialog = widget_from_builder("save_coords_dialog");
   if (! combobox) {
      std::cout << "ERROR:: on_save_coords_dialog_save_button_clicked: bad combobox\n";
   } else {
      int imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
      GtkWidget *chooser = coot_save_coords_chooser(); // uses builder
      g_object_set_data(G_OBJECT(chooser), "imol", GINT_TO_POINTER(imol));
      set_file_for_save_filechooser(chooser);
      gtk_widget_show(chooser);
      set_transient_and_position(COOT_UNDEFINED_WINDOW, chooser);
   }
   gtk_widget_hide(dialog);

}

/*                   do we need this function?                  */
extern "C" G_MODULE_EXPORT
void
on_save_coord_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
}


extern "C" G_MODULE_EXPORT
void
on_save_coords_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{

   // Is this function used now?
   
   GtkWidget *widget = widget_from_builder("save_coords_fileselection1");
   // this looks wrong
   gpointer o = g_object_get_data(G_OBJECT(widget), "stuff");
   free(o);
   gtk_widget_hide(widget);

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
  gtk_widget_show(w);
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
                                        gpointer         user_data)
{

   // 20211026-PE I don't think that this is used. Call back is in glade-callbacks-main-window.cc
   // GtkWidget *widget = wrapped_create_refine_params_dialog();
   // GtkWidget *widget = widget_from_builder("refinement_and_regularization_parameters_dialog");
   // gtk_widget_show(widget);
}

#include "graphics-info.h"

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

  gtk_widget_hide(widget);
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
on_distance1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   do_distance_define();
}


extern "C" G_MODULE_EXPORT
void
on_angle1_activate                     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   do_angle_define();

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_pepflip_togglebutton_toggled (GtkButton       *button,
					       gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    do_pepflip(1);
  else
     do_pepflip(0);

}

/* accept_reject_refinement_accept_button */

extern "C" G_MODULE_EXPORT
void
on_accept_reject_refinement_accept_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "accept_reject_refinement_dialog");

  /* Pressing Return while focus is on the Accept/Reject dialog brings us here. */

  stop_refinement_internal();
  c_accept_moving_atoms();
  save_accept_reject_dialog_window_position(window);
  set_accept_reject_dialog(NULL);
  gtk_widget_hide(window);
}


extern "C" G_MODULE_EXPORT
void
on_accept_reject_refinement_reject_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "accept_reject_refinement_dialog");
  save_accept_reject_dialog_window_position(window);
  /*   clear_up_moving_atoms(); done in destroy of the window */

  stop_refinement_internal();
  clear_up_moving_atoms();
  gtk_widget_hide(window);
}

extern "C" G_MODULE_EXPORT
void
on_accept_reject_refinement_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{


  /* Pressing Escape while focus is on the Accept/Reject dialog brings us here. */

  /* 20070801 To Fix a crash reported by "Gajiwala, Ketan", we need to
     reset the value for graphics_info_t::accept_reject_dialog (it's
     gone now).  And I suppose that we should clean up (and undisplay)
     the intermediate atoms too.
 */

  set_accept_reject_dialog(NULL);
  stop_refinement_internal();
  /* I want to merely clear the stick restraint, not refine again after I did that */

  /* clear_atom_pull_restraint_on_accept_reject_destroy(); */

  clear_up_moving_atoms();
}


/* accept_reject_refinement_docked_accept_button */

extern "C" G_MODULE_EXPORT
void
on_accept_reject_refinement_docked_accept_button_clicked (GtkButton       *button,
							  gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "accept_reject_dialog_frame_docked");
  accept_regularizement();
  if (accept_reject_dialog_docked_show_state() == 1) {
    gtk_widget_set_sensitive(window, FALSE);
    set_accept_reject_dialog(NULL);
    clear_up_moving_atoms();
    GtkWidget *p = main_window();
    gtk_widget_grab_focus(p);
  } else {
    gtk_widget_hide(window);
  }
}


extern "C" G_MODULE_EXPORT
void
on_accept_reject_refinement_docked_reject_button_clicked (GtkButton       *button,
							  gpointer         user_data)
{
   GtkWidget *window = widget_from_builder(
                                     "accept_reject_dialog_frame_docked");
   /* we only hide the widget, I guess we may have to clean up as well */
   /* done below in hide callback */
   if (accept_reject_dialog_docked_show_state() == 1) {
      gtk_widget_set_sensitive(window, FALSE);
      set_accept_reject_dialog(NULL);
      clear_up_moving_atoms();
      GtkWidget *p = main_window();
      gtk_widget_grab_focus(p);
   } else {
      gtk_widget_hide(window);
   }
}

extern "C" G_MODULE_EXPORT
void
on_accept_reject_dialog_frame_docked_hide
                                        (GtkWidget       *widget,
                                        gpointer         user_data)
{
  set_accept_reject_dialog(NULL);
  clear_up_moving_atoms();
  GtkWidget *p = main_window();
  gtk_widget_grab_focus(p);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_find_waters_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   // GtkWidget *widget = create_find_waters_dialog();
   GtkWidget *widget = widget_from_builder("find_waters_dialog");

   fill_find_waters_dialog(widget);
   gtk_widget_show(widget);

}


extern "C" G_MODULE_EXPORT
void
on_find_waters_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("find_waters_dialog");
   execute_find_waters();
   gtk_widget_hide(widget);
}


extern "C" G_MODULE_EXPORT
void
on_find_waters_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget;

   widget = widget_from_builder(
			  "find_waters_dialog");

   gtk_widget_hide(widget);
}


extern "C" G_MODULE_EXPORT
void
on_fast_sss_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  widget = widget_from_builder("fast_ss_search_dialog");
  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_fast_sss_dialog_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *dialog;
  GtkWidget *helix_checkbutton;
  GtkWidget *helix_temp_combobox;
  GtkWidget *helix_noaa_combobox;
  GtkWidget *strand_checkbutton;
  GtkWidget *strand_temp_combobox;
  GtkWidget *strand_noaa_combobox;
  GtkWidget *radius_checkbutton;
  GtkWidget *radius_combobox;
  int use_helix  = 0;
  int helix_target;
  int helix_length;
  int use_strand = 0;
  int strand_target;
  int strand_length;
  float radius = 0.;


  dialog = widget_from_builder("fast_ss_search_dialog");

  helix_checkbutton   = widget_from_builder("fast_sss_dialog_helix_checkbutton");
  helix_temp_combobox = widget_from_builder("fast_sss_dialog_helix_template_combobox");
  helix_noaa_combobox = widget_from_builder("fast_sss_dialog_helix_no_aa_combobox");
  strand_checkbutton   = widget_from_builder("fast_sss_dialog_strand_checkbutton");
  strand_temp_combobox = widget_from_builder("fast_sss_dialog_strand_template_combobox");
  strand_noaa_combobox = widget_from_builder("fast_sss_dialog_strand_no_aa_combobox");
  radius_checkbutton   = widget_from_builder("fast_sss_dialog_local_checkbutton");
  radius_combobox = widget_from_builder("fast_sss_dialog_radius_combobox");

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(helix_checkbutton))) {
    use_helix = 1;
  }
  helix_length = gtk_combo_box_get_active(GTK_COMBO_BOX(helix_noaa_combobox)) *2 + 5;
  helix_target = gtk_combo_box_get_active(GTK_COMBO_BOX(helix_temp_combobox));

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(strand_checkbutton))) {
    use_strand = 1;
  }
  strand_length = gtk_combo_box_get_active(GTK_COMBO_BOX(strand_noaa_combobox)) *2 + 5;
  strand_target = gtk_combo_box_get_active(GTK_COMBO_BOX(strand_temp_combobox));

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(radius_checkbutton))) {
    radius = (gtk_combo_box_get_active(GTK_COMBO_BOX(radius_combobox)) + 1. ) * 10.;
  }



  //g_print("INFO:: run fast secondary structure search with params:\n");
  //g_print("INFO:: helix: %i %i %i\n", use_helix, helix_length, helix_target);
  //g_print("INFO:: strand: %i %i %i\n", use_strand, strand_length, strand_target);
  //g_print("INFO:: radius: %f \n", radius);

  find_secondary_structure_local(use_helix, helix_length, helix_target,
				 use_strand, strand_length, strand_target,
				 radius);

  gtk_widget_hide(dialog);

}


extern "C" G_MODULE_EXPORT
void
on_fast_sss_dialog_citation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *dialog;
  GtkWidget *toolbutton;
  dialog = wrapped_create_coot_references_dialog();
  toolbutton = widget_from_builder("coot_references_buccaneer_toolbutton");
  fill_references_notebook(GTK_BUTTON(toolbutton), COOT_REFERENCE_BUCCANEER);

}


extern "C" G_MODULE_EXPORT
void
on_environment_distances1_activate     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   //  GtkWidget *widget = create_environment_distance_dialog();
   GtkWidget *widget = widget_from_builder("environment_distance_dialog");
   fill_environment_widget(widget);
   gtk_widget_show(widget);

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
on_refine_params_use_torsions_checkbutton_toggled (GtkToggleButton *togglebutton,
						   gpointer         user_data)
{
#if 0 // 20211026-PE  old
   GtkWidget *omega_checkbutton = widget_from_builder("refine_params_use_peptide_omegas_checkbutton");
   GtkWidget *phi_psi_restraints_vbox = widget_from_builder("peptide_torsions_restraints_vbox");

   do_torsions_toggle(GTK_WIDGET(togglebutton));

   /* We don't see these widgets, currently */
   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(omega_checkbutton,       TRUE);
   } else {
      gtk_widget_set_sensitive(omega_checkbutton,       FALSE);
      gtk_widget_set_sensitive(phi_psi_restraints_vbox, FALSE);
      /* no rama restraints */
   }
#endif

   do_torsions_toggle(GTK_WIDGET(togglebutton));
}

extern "C" G_MODULE_EXPORT
void
on_refine_params_use_planar_peptides_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   if (gtk_toggle_button_get_active(togglebutton)) {
      add_planar_peptide_restraints();
   } else {
      remove_planar_peptide_restraints();
   }
}



extern "C" G_MODULE_EXPORT
void
on_refine_params_use_trans_peptide_restraints_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   if (gtk_toggle_button_get_active(togglebutton)) {
      set_use_trans_peptide_restraints(1);
   } else {
      set_use_trans_peptide_restraints(0);
   }
}



extern "C" G_MODULE_EXPORT
void
on_refine_params_use_initial_pos_checkbutton_toggled (GtkToggleButton *togglebutton,
						      gpointer         user_data)
{

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

/* I mean phi/psi torsions */
extern "C" G_MODULE_EXPORT
void
on_refine_params_use_peptide_torsions_checkbutton_toggled (GtkToggleButton *togglebutton,
							   gpointer         user_data)
{

   GtkWidget *frame = widget_from_builder("peptide_torsions_restraints_vbox");
   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(frame, TRUE);
   } else {
      gtk_widget_set_sensitive(frame, FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_import_cif_dictionary1_activate     (GMenuItem     *menuitem,
                                                            gpointer         user_data)
{
  open_cif_dictionary_file_selector_dialog();
}



extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *fileselection;

  fileselection = widget_from_builder("cif_dictionary_fileselection");
  gtk_widget_hide(fileselection);

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

   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
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
   gtk_widget_hide(window);
   do_mutation("VAL", istate);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_rot_trans_togglebutton_toggled (GtkButton       *button,
						 gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
    do_rot_trans_setup(1);
  } else {
    do_rot_trans_setup(0);
  }

}




extern "C" G_MODULE_EXPORT
void on_rotate_translate_obj_ok_button_clicked (GtkButton       *button,
                                                                    gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder("rotate_translate_obj_dialog");

  // graphics_unsetup_rotate_translate_buttons(widget);
  rot_trans_reset_previous();
  accept_regularizement();
  store_window_position(COOT_ROTATE_TRANSLATE_DIALOG, widget);
  gtk_widget_hide(widget);
  clear_up_moving_atoms(); // redraw done here

}


extern "C" G_MODULE_EXPORT
void on_rotate_translate_obj_cancel_button_clicked (GtkButton       *button,
                                                                        gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder("rotate_translate_obj_dialog");

 /*   clear_moving_atoms_object(); // redraw done here */
  clear_up_moving_atoms();
  store_window_position(COOT_ROTATE_TRANSLATE_DIALOG, widget);
  gtk_widget_hide(widget);
  normal_cursor();
}


extern "C" G_MODULE_EXPORT
void
on_run_script1_activate (GMenuItem     *menuitem,
                                             gpointer         user_data) {

   GtkWidget *widget = coot_run_script_chooser();
   add_filename_filter_button(widget, COOT_SCRIPTS_FILE_SELECTION);
   gtk_widget_show(widget);
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


/* Old unused function, see new_close_molecules_dialog(), new_close_molecules() */
extern "C" G_MODULE_EXPORT
void
on_close_molecule1_activate            (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

/* Old Stylee */
/*   GtkWidget *widget = create_close_molecule_dialog(); */
/*   GtkWidget *optionmenu = widget_from_builder("close_molecule_optionmenu"); */
/*    fill_close_option_menu_with_all_molecule_options(optionmenu); */
/*   gtk_widget_show(widget); */

   GtkWidget *widget = wrapped_create_new_close_molecules_dialog(); // uses builder
   gtk_widget_show(widget);
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
/*   gtk_widget_hide(window); */
}


extern "C" G_MODULE_EXPORT
void
on_close_molecule_cancel_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *window = widget_from_builder(
				    "close_molecule_dialog");
  gtk_widget_hide(window);

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
   store_window_position(COOT_DELETE_WINDOW, widget);
   // store_delete_item_widget(NULL);
   gtk_widget_hide(widget);
}


extern "C" G_MODULE_EXPORT
void
on_hints1_activate                     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *widget = create_hints_dialog();
   GtkWidget *widget = widget_from_builder("hints_dialog");
   gtk_widget_show(widget);
}


/* this is the Apply button now */
extern "C" G_MODULE_EXPORT
void
on_residue_info_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *widget = widget_from_builder("residue_info_dialog");
   GtkWidget *widget = widget_from_builder("residue_info_dialog");
   apply_residue_info_changes(widget);
   /*    gtk_widget_hide(widget); not now that it's the Apply button*/

}


extern "C" G_MODULE_EXPORT
void
on_residue_info_cancel_button_clicked  (GtkButton       *button,
                                                            gpointer         user_data)
{
   // GtkWidget *widget = widget_from_builder("residue_info_dialog");
   GtkWidget *widget = widget_from_builder("residue_info_dialog");

   residue_info_release_memory(widget);  // Hmmm! that seems dangerous

   // gtk_widget_hide(widget); not now we use builder
   gtk_widget_hide(widget);
   unset_residue_info_widget();
}

extern "C" G_MODULE_EXPORT
void
on_residue_info_dialog_destroy         (GtkWidget       *object,
                                        gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("residue_info_dialog");
   residue_info_release_memory(widget);
   clear_residue_info_edit_list();
   unset_residue_info_widget();

}



extern "C" G_MODULE_EXPORT
void
on_hints_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("hints_dialog");
   gtk_widget_hide(widget);

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
   store_window_position(COOT_ROTAMER_SELECTION_DIALOG, dialog);
   gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_rotamer_selection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder(
				    "rotamer_selection_dialog");
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

   store_window_position(COOT_ROTAMER_SELECTION_DIALOG, dialog);
   gtk_widget_hide(dialog);
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


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_go_to_atom_chain_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
    apply_go_to_atom_values(widget_from_builder("goto_atom_window"));
  }
  return FALSE;
}
#endif


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_go_to_atom_residue_entry_key_press_event (GtkWidget       *widget,
					     GdkEventKey     *event,
					     gpointer         user_data)
{


  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
    apply_go_to_atom_values(widget_from_builder("goto_atom_window"));
  }
  return FALSE;
}
#endif


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_go_to_atom_atom_name_entry_key_press_event (GtkWidget       *widget,
					       GdkEventKey     *event,
					       gpointer         user_data)
{

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
    apply_go_to_atom_values(widget_from_builder("goto_atom_window"));
  }
  return FALSE;
}
#endif

extern "C" G_MODULE_EXPORT
gboolean
on_go_to_atom_show_HOH_checkbutton_toggled(GtkToggleButton *togglebutton,
                                           gpointer         user_data) {

   if (gtk_toggle_button_get_active(togglebutton))
      std::cout << "Now show the waters" << std::endl;
   else
      std::cout << "Now hide the waters" << std::endl;

   return FALSE;
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
   gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_unsaved_changes_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("unsaved_changes_dialog");
   gtk_widget_hide(dialog);
   coot_clear_backup_or_real_exit(0);
}


extern "C" G_MODULE_EXPORT
void
on_environment_distance_label_atom_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
    set_environment_distances_label_atom(1);
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
on_baton_tip_previous_button_clicked      (GtkButton       *button,
					   gpointer         user_data) {

   baton_tip_previous();
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
  gtk_widget_hide(widget);

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
    gtk_widget_show(widget);
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
on_environment_distance_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   toggle_environment_show_distances(togglebutton);
}


extern "C" G_MODULE_EXPORT
void
on_environment_distance_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget;
   widget = widget_from_builder("environment_distance_dialog");
   execute_environment_settings(GTK_WIDGET(button));
   gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_as_calphas_checkbutton_toggled (GtkToggleButton *togglebutton,
						 gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    symmetry_as_calphas(0, 1);
  } else {
    symmetry_as_calphas(0, 0);
  }
}

extern "C" G_MODULE_EXPORT
gboolean
on_show_symmetry_window_delete_event(GtkWidget       *widget,
                                                         GdkEvent        *event,
                                                         gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}



extern "C" G_MODULE_EXPORT
void
on_coordinates_recentring1_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   // GtkWidget *widget = create_read_pdb_recentre_dialog();
   GtkWidget *widget = widget_from_builder("read_pdb_recentre_dialog");
   /* lookup the toggle widgets here and
      set the acording to
      recentre_on_read_pdb()  */
   GtkWidget *yes_radio_button = widget_from_builder(
                                                     "read_pdb_recentre_yes_radiobutton");
   GtkWidget *no_radio_button = widget_from_builder(
                                                    "read_pdb_recentre_no_radiobutton");

   if (recentre_on_read_pdb()) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(yes_radio_button), TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(no_radio_button), FALSE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(yes_radio_button), FALSE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(no_radio_button), TRUE);
   }


   gtk_widget_show(widget);
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
   // gtk_widget_hide(widget);

   GtkWidget *w = widget_from_builder("read_pdb_recentre_dialog");
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_pointer_atom_type_cancel_button_clicked(GtkButton       *button,
                                                               gpointer         user_data) {

   // GtkWidget *widget = widget_from_builder("pointer_atom_type_dialog");
   // gtk_widget_hide(widget);

   GtkWidget *w = widget_from_builder("pointer_atom_type_dialog");
   gtk_widget_hide(w);

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

  // gtk_widget_hide(dialog);
  gtk_widget_hide(dialog);
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
on_model_refine_dialog_refmac_button_clicked (GtkButton       *button,
					      gpointer         user_data)
{

   wrapped_create_simple_refmac_dialog(); // uses builder

}



extern "C" G_MODULE_EXPORT
void
on_run_refmac_tls_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // not used
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_twin_checkbutton_toggled (GtkToggleButton *togglebutton,
                                                            gpointer         user_data) {

   // not used
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_sad_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
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
on_run_refmac_mtz_file_radiobutton_toggled
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
on_run_refmac_file_help_button_clicked (GtkButton       *button,
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
on_run_refmac_sad_help_button_clicked  (GtkButton       *button,
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


extern "C" G_MODULE_EXPORT
void
on_run_refmac_nolabels_checkbutton_toggled (GtkToggleButton *togglebutton,
                                            gpointer         user_data)
{
   // not used
}


/* we want to update the phases/hl boxes if the mtz changes */
#if 0				/* don't know what to do with optionmenu changed */
extern "C" G_MODULE_EXPORT
void
on_run_refmac_map_optionmenu_changed   (GtkOptionMenu   *optionmenu,
                                        gpointer         user_data)
{
  //update_refmac_column_labels_frame(optionmenu);

}
#endif

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
extern "C" G_MODULE_EXPORT
void
on_run_refmac_ncs_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // not used
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
on_baton_undo_button_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
   baton_build_delete_last_residue();
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_undo_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{
   apply_undo();
}


extern "C" G_MODULE_EXPORT
void
on_undo_molecule_chooser_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("undo_molecule_chooser_dialog");
   gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_clear_pending_button_clicked (GtkButton *button,
						     gpointer user_data)
{
   clear_pending_picks();
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
  gtk_widget_hide(window);
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
  gtk_widget_hide(window);
}


extern "C" G_MODULE_EXPORT
void
on_skeleton1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   /* we need to fill the optionmenu.  Make it a graphics_info_t this time. */
  GtkWidget *w = wrapped_create_skeleton_dialog();
  gtk_widget_show(w);
}


#if (GTK_MAJOR_VERSION >= 4)
// 20220601-PE  delete this (ancient) function, I think.
#else
extern "C" G_MODULE_EXPORT
void
on_virtual_trackball_menu_pops         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *flat_check_menu_item;
  GtkWidget *spherical_check_menu_item;

  flat_check_menu_item      = widget_from_builder("flat1");
  spherical_check_menu_item = widget_from_builder("spherical_surface1");

/*   printf("got menu items 0x%x 0x%x\n",  */
/* 	 flat_check_menu_item, spherical_check_menu_item); */

  if (vt_surface_status() == 1) { /* 1 is VT_FLAT */
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(flat_check_menu_item), TRUE);
  }

  if (vt_surface_status() == 2) { /* 2 is VT_SPHERICAL */
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(spherical_check_menu_item), TRUE);
  }

/*   if (GTK_CHECK_MENU_ITEM(flat_check_menu_item)->active) {  */
/*     printf("flat is active\n"); */
/*   }  */
/*   if (GTK_CHECK_MENU_ITEM(spherical_check_menu_item)->active) {  */
/*     printf("spherical is active\n"); */
/*   }  */
}
#endif

// not use now, I think.
extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_auto_fit_rotamer_togglebutton_toggled (GtkButton *button,
							gpointer user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
     setup_auto_fit_rotamer(1);
  else
     setup_auto_fit_rotamer(0);
}




extern "C" G_MODULE_EXPORT
void
on_residue_info2_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_residue_info_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_import_all_dictionary_cifs1_activate (GMenuItem     *menuitem,
					 gpointer         user_data)
{
  import_all_refmac_cifs();
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_destroy         (GtkWidget       *object,
                                        gpointer         user_data)
{
   unset_model_fit_refine_dialog();
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
on_residue_parameters1_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_residue_info_dialog();
}



extern "C" G_MODULE_EXPORT
void
on_crosshairs1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *widget = widget_from_builder("crosshairs_dialog");
   GtkWidget *button = widget_from_builder("crosshairs_on_radiobutton");

   if (draw_crosshairs_state()) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
   }
   gtk_widget_show(widget);
}


extern "C" G_MODULE_EXPORT
void
on_crosshairs_on_radiobutton_toggled   (GtkToggleButton *togglebutton,
                                                            gpointer         user_data) {

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      set_draw_crosshairs(1);
   }
}


extern "C" G_MODULE_EXPORT
void
on_crosshairs_off_radiobutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    set_draw_crosshairs(0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_display_crosshairs_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *dialog = widget_from_builder("crosshairs_dialog");

  gtk_widget_hide(dialog);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_mutate_auto_fit_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
    setup_mutate_auto_fit(1);
  else
     setup_mutate_auto_fit(0);
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
  gtk_widget_hide(dialog);

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

  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_edit_phi_psi_togglebutton_toggled (GtkToggleButton *togglebutton,
							  gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    setup_edit_phi_psi(1);
  } else {
    setup_edit_phi_psi(0);
  }
}




extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_redo_button_clicked (GtkButton       *button,
					    gpointer         user_data) {
  apply_redo();
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_edit_chi_angles_togglebutton_toggled (GtkToggleButton *togglebutton,
							     gpointer         user_data) {
  if (gtk_toggle_button_get_active(togglebutton)) {
    setup_edit_chi_angles(1);
  } else {
    setup_edit_chi_angles(0);
    set_show_chi_angle_bond(0);
  }
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





#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_window1_configure_event             (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  /*
  printf("on_window1_configure_event start\n");
  gint upositionx, upositiony;
  GtkWindow *window = gtk_widget_get_window(widget);
  printf("on_window1_configure_event with window: 0x%d\n", window);
  gtk_window_get_position(window, &upositionx, &upositiony); // not a window strangely.
  store_graphics_window_position(upositionx, upositiony);
  printf("on_window1_configure_event done\n");
  */
  return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_run_state_file_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *dialog = widget_from_builder("run_state_file_dialog");    // old glade style
  GtkWidget *dialog = widget_from_builder("run_state_file_dialog");
  gtk_widget_hide(dialog);
  run_state_file();
}


extern "C" G_MODULE_EXPORT
void
on_run_state_file_cancel_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
   // GtkWidget *dialog = widget_from_builder("run_state_file_dialog"); // old glade style
  GtkWidget *dialog = widget_from_builder("run_state_file_dialog");
  gtk_widget_hide(dialog);
  gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
gboolean
on_run_state_file_dialog_delete_event(GtkWidget       *widget,
                                                          GdkEvent        *event,
                                                          gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsions_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
  clear_moving_atoms_object();
  /* FIXME: also clear out the edib backbone ramaplot, if it exists. */
  /*   destroy_edit_backbone_rama_plot(); */
}


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsion_rotate_peptide_button_pressed
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
   int ix, iy;
   GdkModifierType state;
   GdkWindow *window = gtk_widget_get_window(GTK_WIDGET(button));
   gdk_window_get_pointer(window, &ix, &iy, &state);
   set_backbone_torsion_peptide_button_start_pos(ix, iy);
#endif
}


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsion_rotate_peptide_button_released
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  /* ignored */
}


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_edit_backbone_torsion_rotate_peptide_button_motion_notify_event
                                        (GtkWidget       *widget,
                                        GdkEventMotion  *event,
                                        gpointer         user_data)
{
   int ix, iy;
   GdkModifierType state;
   GdkWindow *window = gtk_widget_get_window(widget);
   gdk_window_get_pointer(window, &ix, &iy, &state);
   change_peptide_peptide_by_current_button_pos(ix, iy);
   return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_pressed
                                        (GtkButton       *button,
                                         gpointer         user_data) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
#else
   int ix, iy;
   GdkModifierType state;
   GdkWindow *window = gtk_widget_get_window(GTK_WIDGET(button));
   gdk_window_get_pointer(window, &ix, &iy, &state);
   set_backbone_torsion_carbonyl_button_start_pos(ix, iy);
#endif
}


extern "C" G_MODULE_EXPORT
void
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_released
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  /* ignored */
}


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_motion_notify_event
                                        (GtkWidget       *widget,
                                        GdkEventMotion  *event,
                                        gpointer         user_data)
{
  int ix, iy;
  GdkModifierType state;
  GdkWindow *window = gtk_widget_get_window(widget);
  gdk_window_get_pointer(window, &ix, &iy, &state);
/*   printf("button moved to %d %d \n", ix, iy); */

  change_peptide_carbonyl_by_current_button_pos(ix, iy);
  return FALSE;
}
#endif


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
  gtk_widget_hide(window);
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
   //gtk_widget_hide(window);
   gtk_widget_hide(window);

}


extern "C" G_MODULE_EXPORT
void
on_sequence_view1_activate             (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   add_on_sequence_view_choices();
}


extern "C" G_MODULE_EXPORT
void
on_clear_simple_distances2_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   clear_measure_distances();
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
on_select_map_for_fitting_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  // this doesn't do anything because gtk_dialog_run() is used

//  GtkWidget *widget = widget_from_builder(
//				    "select_fitting_map_dialog");

 // gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
gboolean
on_select_fitting_map_dialog_delete_event(GtkWidget       *widget,
                                                              GdkEvent        *event,
                                                              gpointer         user_data) {
   gtk_widget_hide(widget);
   return TRUE;
}



extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_map_select_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   show_select_map_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_clear_atom_labels1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  remove_all_atom_labels();
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
on_model_refine_dialog_add_alt_conf_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  altconf();
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
   gtk_widget_show(widget);
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("run_refmac_help_dialog");
  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_nolabels_help_button_clicked
					(GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *widget = create_run_refmac_nolabels_help_dialog();
   GtkWidget *widget = widget_from_builder("run_refmac_nolabels_help_dialog");
   gtk_widget_show(widget);
}


extern "C" G_MODULE_EXPORT
void
on_run_refmac_nolabels_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("run_refmac_nolabels_help_dialog");
  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_no_restraints_info_dialog_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = widget_from_builder("no_restraints_info_dialog");
  gtk_widget_hide(window);

}


extern "C" G_MODULE_EXPORT
void
on_no_cif_dictionary_bonds_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = widget_from_builder("no_cif_dictionary_bonds_dialog");
  gtk_widget_hide(window);

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
on_ligand_big_blob_dismiss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = widget_from_builder( "ligand_big_blob_dialog");
  gtk_widget_hide(window);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("edit_chi_angles_dialog");
  accept_regularizement();
  unset_moving_atom_move_chis();
  store_window_position(COOT_EDIT_CHI_DIALOG, widget);
  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("edit_chi_angles_dialog");
  clear_up_moving_atoms();	/* and remove the graphics object */
  unset_moving_atom_move_chis();
  store_window_position(COOT_EDIT_CHI_DIALOG, widget);
  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder("check_waters_dialog");
  do_check_waters_by_widget(widget);
  gtk_widget_hide(widget);
}


extern "C" G_MODULE_EXPORT
void
on_check_waters_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = widget_from_builder("check_waters_dialog");
  gtk_widget_hide(widget);

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
on_geometry_distance_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
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
  store_window_position(COOT_DISTANCES_ANGLES_WINDOW, dialog);
  store_geometry_dialog(NULL);
  gtk_widget_hide(dialog);
#endif

  GtkWidget *frame = widget_from_builder("geometry_frame");
  store_geometry_dialog(NULL);
  gtk_widget_hide(frame);
}


extern "C" G_MODULE_EXPORT
void
on_geometry_angle_togglebutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    do_angle_define();

}


extern "C" G_MODULE_EXPORT
void
on_new_ligands_info_dialog_ok_button_clicked (GtkButton       *button,
                                                                  gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("new_ligands_info_dialog");
   // gtk_widget_hide(w);
   GtkWidget *w = widget_from_builder("new_ligands_info_dialog");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_no_new_ligands_info_dialog_ok_button_clicked (GtkButton       *button,
                                                                     gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("no_new_ligands_info_dialog");
   // gtk_widget_hide(w);
   GtkWidget *w = widget_from_builder("no_new_ligands_info_dialog");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_zoom1_activate                      (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   // GtkWidget *w = create_zoom_dialog();
   GtkWidget *w = widget_from_builder("zoom_dialog");
   set_zoom_adjustment(w);
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_zoom_dialog_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("zoom_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data)
{
   /* needs to set widget */
   /*  store_window_position(COOT_EDIT_CHI_DIALOG, widget); */
  unset_moving_atom_move_chis();
  set_show_chi_angle_bond(0);
}


extern "C" G_MODULE_EXPORT
void
on_check_waters_low_occ_dist_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_zero_occ_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_get_monomer1_activate               (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *dialog = wrapped_create_get_monomer_dialog();
  gtk_widget_show(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_get_monomer_ok_button_clicked(GtkButton       *button,
                                 gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("get_monomer_entry");
   if (entry) {
      handle_get_monomer_code(entry);
   }
   GtkWidget *frame = widget_from_builder("get_monomer_frame");
   gtk_widget_hide(frame);
}

extern "C" G_MODULE_EXPORT
void on_generic_overlay_frame_cancel_button_clicked(GtkButton       *button,
                                          gpointer         user_data) {
   GtkWidget* frame_widget = GTK_WIDGET(user_data);
   if(frame_widget) {
      gtk_widget_hide(frame_widget);
   } else {
      g_error("'user_data' is NULL. Cannot hide overlay frame.");
   }
   
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_libcheck_monomer_entry_key_press_event (GtkWidget       *widget,
					   GdkEventKey     *event,
					   gpointer         user_data)
{

   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      handle_get_libcheck_monomer_code(widget);
  }

  return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_recover_coordinates_ok_button_clicked(GtkButton       *button,
                                         gpointer         user_data) {

   GtkWidget *widget = widget_from_builder("recover_coordinates_dialog");

  execute_recover_session(widget); /* widget needed for lookup of user data */
  gtk_widget_hide(widget);

}


extern "C" G_MODULE_EXPORT
void
on_recover_coordinates_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = widget_from_builder("recover_coordinates_dialog");
  gtk_widget_hide(widget);
}


extern "C" G_MODULE_EXPORT
void
on_recover_session1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  recover_session();
}


extern "C" G_MODULE_EXPORT
void
on_centre_atom_label1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *widget = create_centre_atom_label_dialog();
   GtkWidget *widget = widget_from_builder("centre_atom_label_dialog");
   GtkWidget *on  = widget_from_builder("centre_atom_label_radiobutton_on");
   GtkWidget *off = widget_from_builder("centre_atom_label_radiobutton_off");
   int v = centre_atom_label_status();
   if (v) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off), TRUE);
   }
   gtk_widget_show(widget);
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
   gtk_widget_hide(dialog);
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
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_help_chi_angles_dismiss_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("chi_angle_help_dialog");
  gtk_widget_hide(w);
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
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_nothing_to_recover_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("nothing_to_recover_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_superpose_dialog_superpose_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_dialog");
  execute_superpose(w);
  gtk_widget_hide(w);

}

extern "C" G_MODULE_EXPORT
void
on_superpose_dialog_cancel_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_superpose_nonsense_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_nonsense_dialog");
  gtk_widget_hide(w);

}

extern "C" G_MODULE_EXPORT
void
on_superpose_nonsense_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("superpose_dialog");
  gtk_widget_hide(w);

}

extern "C" G_MODULE_EXPORT
void
on_add_terminal_residue_finds_none_ok_button_clicked (GtkButton       *button,
						      gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("add_terminal_residue_finds_none_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_single_map_sigma_checkbutton_toggled (GtkToggleButton *togglebutton,
                                                             gpointer         user_data)
{
   GtkWidget *window = widget_from_builder("single_map_properties_dialog");
   GtkWidget *entry  = widget_from_builder("single_map_sigma_step_entry");
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));
   const char *text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
   float v;

   if (gtk_toggle_button_get_active(togglebutton)) {
      if (text) {
         v = atof(text);
         set_contour_by_sigma_step_by_mol(v, 1, imol);
         gtk_widget_set_sensitive(entry, TRUE);
      }
   } else {
      /* 0.0 is ignored. */
      set_contour_by_sigma_step_by_mol(0.0, 0, imol);
      gtk_widget_set_sensitive(entry, FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_single_map_properties_absolute_radiobutton_toggled (GtkToggleButton *togglebutton,
                                                                           gpointer         user_data)
{
   std::cout << "Absolute button toggled" << std::endl;
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton), "imol"));
   if (gtk_toggle_button_get_active(togglebutton)) {
      GtkWidget *entry = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton), "contour_level_entry"));
      const char *text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      float f = coot::util::string_to_float(text);
      set_contour_by_sigma_step_by_mol(f, 1, imol);
   }

}

void handle_map_properties_fresnel_change(int imol, GtkWidget *togglebutton) {

   molecule_class_info_t &m = graphics_info_t::molecules[imol];
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      std::cout << "Here B in handle_map_properties_fresnel_change() " << togglebutton << std::endl;
      GtkWidget *bias_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton),  "bias_entry"));
      GtkWidget *scale_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton), "scale_entry"));
      GtkWidget *power_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton), "power_entry"));
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

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_map_properties_dialog_fresnel_bias_entry_key_press_event (GtkWidget       *widget,
                                                                                 GdkEventKey     *event,
                                                                                 gpointer         user_data)
{

   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "imol"));
      GtkWidget *togglebutton = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "fresnel_checkbutton"));
      handle_map_properties_fresnel_change(imol, togglebutton);
   }
   return FALSE;
}
#endif

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_map_properties_dialog_fresnel_scale_entry_key_press_event (GtkWidget       *widget,
                                                              GdkEventKey     *event,
                                                              gpointer         user_data)
{

   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "imol"));
      GtkWidget *togglebutton = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "fresnel_checkbutton"));
      handle_map_properties_fresnel_change(imol, togglebutton);
   }
   return FALSE;
}
#endif


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_map_properties_dialog_fresnel_power_entry_key_press_event (GtkWidget       *widget,
                                                                                  GdkEventKey     *event,
                                                                                  gpointer         user_data)
{
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "imol"));
      GtkWidget *togglebutton = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "fresnel_checkbutton"));
      handle_map_properties_fresnel_change(imol, togglebutton);
   }
   return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_fresnel_state_checkbutton_toggled (GtkToggleButton *togglebutton,
                                                            gpointer         user_data)
{

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton), "imol"));
   handle_map_properties_fresnel_change(imol, GTK_WIDGET(togglebutton));
}



extern "C" G_MODULE_EXPORT
void
on_no_bad_chiral_volumes_dialog_ok_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("no_bad_chiral_volumes_dialog");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_check_chiral_volumes_ok_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("check_chiral_volumes_dialog");
   check_chiral_volumes_from_widget(w);
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_check_chiral_volumes_cancel_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("check_chiral_volumes_dialog");
   gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_find_bad_chiral_atoms1_activate     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   // GtkWidget *w = create_check_chiral_volumes_dialog();
   GtkWidget *w = widget_from_builder("check_chiral_volumes_dialog");
   fill_chiral_volume_molecule_combobox(w);
   gtk_widget_show(w);
}



extern "C" G_MODULE_EXPORT
void
on_chiral_volume_baddies_dialog_cancel_button_clicked (GtkButton       *button,
						       gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("bad_chiral_volumes_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_bad_chiral_volumes_dialog_response(GtkDialog       *dialog,
                                                          gint             response_id,
                                                          gpointer         user_data) {

   if (response_id == GTK_RESPONSE_CLOSE)
      gtk_widget_hide(GTK_WIDGET(dialog));

}

extern "C" G_MODULE_EXPORT
void
on_rigid_body_refinement_failed_dialog_ok_button_clicked (GtkButton *button,
                                                                              gpointer user_data)
{
  GtkWidget *w = widget_from_builder(
			       "rigid_body_refinement_failed_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_baton_mode_calculate_skeleton_ok_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  GtkWidget *w = widget_from_builder(
			       "baton_mode_make_skeleton_dialog");
  baton_mode_calculate_skeleton(w); /* get the imol from here */
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_baton_mode_calculate_skeleton_cancel_button_clicked (GtkButton       *button,
							gpointer         user_data)
{
  GtkWidget *w = widget_from_builder(
			       "baton_mode_make_skeleton_dialog");
  gtk_widget_hide(w);
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
    gtk_widget_hide(frame);
   } else {
     gtk_widget_show(frame);
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
  GtkWidget *frame = widget_from_builder(
				   "resolution_limits_hbox");
  if (gtk_toggle_button_get_active(togglebutton))
     gtk_widget_set_sensitive(frame, TRUE);
  else
     gtk_widget_set_sensitive(frame, FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_merge_molecules_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("merge_molecules_dialog");
  do_merge_molecules(w);
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_merge_molecules_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("merge_molecules_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("mutate_sequence_dialog");
   do_mutate_sequence(w);
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("mutate_sequence_dialog");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_merge_molecules1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_mutate_molecule1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_mutate_sequence_dialog();
   gtk_widget_show(w);
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
on_draw_hydrogens_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

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

/* radio button is the N-terminal button */
extern "C" G_MODULE_EXPORT
void
on_renumber_residue_range_radiobutton_1_toggled(GtkToggleButton *togglebutton,
                                                                    gpointer         user_data) {

   GtkWidget *entry_1 = widget_from_builder("renumber_residue_range_resno_1_entry");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      gtk_widget_set_sensitive(GTK_WIDGET(entry_1), FALSE);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(entry_1), TRUE);
   }
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
     gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_renumber_residue_range_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("renumber_residue_range_dialog");
  gtk_widget_hide(w);

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
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_add_OXT_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data) {

   GtkWidget *w = widget_from_builder("add_OXT_dialog");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_add_OXT_button_clicked(GtkButton       *button,
                                              gpointer         user_data) {

   GtkWidget *w = wrapped_create_add_OXT_dialog(); // uses builder
   set_transient_for_main_window(w);
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_renumber_residues1_activate(GMenuItem    *menuitem,
                               gpointer      user_data) {

   GtkWidget *w = wrapped_create_renumber_residue_range_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_bond_parameters1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_bond_parameters_dialog();
   gtk_widget_show(w);
}

extern "C" G_MODULE_EXPORT
void
on_bond_parameters_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data) {

   // GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   // apply_bond_parameters(w);
   // gtk_widget_hide(w);
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
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_background_colour1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data) {

   // GtkWidget *wb = widget_from_builder("background_black1");
   // GtkWidget *ww = widget_from_builder("background_white1");

   GtkWidget *wb = widget_from_builder("background_black1");
   GtkWidget *ww = widget_from_builder("background_white1");

#if (GTK_MAJOR_VERSION >= 4)
#else
   if (background_is_black_p())
      /*     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE); */
      /* a GtkRadioMenuItem not a toggle button */
      gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(wb), TRUE);
   else
      gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(ww), TRUE);
#endif
}



extern "C" G_MODULE_EXPORT
void
on_ligand_no_blobs_OK_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   //   GtkWidget *dialog = widget_from_builder("ligand_no_blobs_dialog");
   // gtk_widget_hide(dialog);

   GtkWidget *dialog = widget_from_builder("ligand_no_blobs_dialog");
   gtk_widget_hide(dialog);
}



extern "C" G_MODULE_EXPORT
void
on_new_delete_molecules_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   // GtkWidget *w = widget_from_builder("new_close_molecules_dialog");
   GtkWidget *w = widget_from_builder("new_close_molecules_dialog");
   new_close_molecules(w);

  /*   gtk_widget_hide(w); */
}


extern "C" G_MODULE_EXPORT
void
on_new_delete_molecules_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("new_close_molecules_dialog");
   GtkWidget *w = widget_from_builder("new_close_molecules_dialog");
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot2_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  add_on_rama_choices();

}


extern "C" G_MODULE_EXPORT
void
on_incorrect_chiral_volumes1_activate  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *w = create_check_chiral_volumes_dialog();
   GtkWidget *w = widget_from_builder("check_chiral_volumes_dialog");
   fill_chiral_volume_molecule_combobox(w);
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_unmodelled_blobs1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_unmodelled_blobs_dialog();
  gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_find_blobs_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("unmodelled_blobs_dialog");
  execute_find_blobs_from_widget(w);
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_find_blobs_cancel_button_clicked    (GtkButton       *button,
                                                            gpointer         user_data) {
  GtkWidget *w = widget_from_builder("unmodelled_blobs_dialog");
  gtk_widget_hide(w);

}



extern "C" G_MODULE_EXPORT
void
on_chiral_restraints_problem_ok_button_clicked (GtkButton       *button,
						gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("chiral_restraints_problem_dialog");
  gtk_widget_hide(w);
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
/*   gtk_widget_hide(w); */


}


extern "C" G_MODULE_EXPORT
void
on_check_waters_diff_map_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("check_waters_diff_map_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_by_difference_map_variance1_activate (GMenuItem     *menuitem,
						      gpointer         user_data)
{
  /* GtkWidget *w = wrapped_create_check_waters_diff_map_dialog(); */
  /* gtk_widget_show(w); */
}


extern "C" G_MODULE_EXPORT
void
on_interesting_waters_by_difference_map_check_ok_button_clicked (GtkButton       *button,
								 gpointer         user_data)
{

  GtkWidget *w = widget_from_builder(
			       "interesting_waters_by_difference_map_check_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_nothing_bad_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nothing_bad_dialog");
   gtk_widget_hide(w);
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
on_antialiasing1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *w = create_antialiasing_dialog();
   GtkWidget *w = widget_from_builder("antialiasing_dialog");
   GtkWidget *checkbutton;
   if (do_anti_aliasing_state()) {
      checkbutton = widget_from_builder("antialias_dialog_yes_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
   gtk_widget_show(w);
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
  gtk_widget_hide(w);

}

// ##########################################################################################

//                 is this a good place for geometry_graphs callbacks?

extern "C" G_MODULE_EXPORT
void
on_geometry_graphs_close_button_clicked(GtkButton       *button,
                                                            gpointer         user_data) {

   // this usage of user_data was set using glade.
   GtkWidget *dialog = GTK_WIDGET(user_data);
   if (dialog)
      gtk_widget_hide(dialog);
   else
      std::cout << "ERROR getting dialog in on_geometry_graphs_ok_button_clicked\n";

}


extern "C" G_MODULE_EXPORT
void
on_geometry_graphs_dialog_destroy      (GtkWidget       *object,
                                                            gpointer         user_data) {

  /* This is not causing the GTK_IS_WIDGET (widget) unref failure */

   GtkWidget *w = widget_from_builder("geometry_graphs_dialog");
   if (! w) {
     printf("ERROR:: in on_geometry_graphs_dialog_destroy() failure getting dialog \n");
   } else {
      unset_geometry_graph(w);
      free_geometry_graph(w);
   }

}

// ##########################################################################################


extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coords_fileselection_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("save_symmetry_coords_fileselection");
  save_symmetry_coords_from_filechooser(w);
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coords_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("save_symmetry_coords_fileselection");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coordinates1_activate (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  setup_save_symmetry_coords();
}


extern "C" G_MODULE_EXPORT
void
on_experimental1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_geometry_analysis1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *menu = widget_from_builder("geometry_analysis1");
   if (menu) {
      const char *type = "geometry";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_geometry_analysis1_activate\n");
   }

}


extern "C" G_MODULE_EXPORT
void
on_peptide_omega_analysis1_activate    (GMenuItem     *menuitem,
                                                            gpointer         user_data) {

   GtkWidget *menu = widget_from_builder("peptide_omega_analysis1");
   if (menu) {
      const char *type = "omega";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      std::cout << "ERROR:: failed to get menu in on_peptide_omega_analysis1_activate" << std::endl;
   }

}

extern "C" G_MODULE_EXPORT
void
on_pukka_puckers_1_activate(GMenuItem     *menuitem,
                                                gpointer         user_data) {

   GtkWidget *menu = widget_from_builder("pukka_puckers_1");
   if (menu) {
      const char *type = "puckers";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      std::cout << "ERROR:: failed to get menu in on_pukka_puckers_activate" << std::endl;
   }
}


extern "C" G_MODULE_EXPORT
void
on_peptide_flips_from_difference_map1_activate_gtkbuilder_glade(GMenuItem     *menuitem,
                                                                gpointer         user_data)
{
   pepflips_by_difference_map_dialog(); // in c-interface-gui. Sets data for the
                                        // model_combobox and the map_combobox
}

extern "C" G_MODULE_EXPORT
void
on_pepflips_by_difference_map_dialog_close (GtkDialog *dialog,
                                                                gpointer   user_data) {
   gtk_widget_hide(GTK_WIDGET(dialog));
}

// use a header for this - which one? c-interface-gui.hh?
void pepflips_by_difference_map_results_dialog(int imol_coords, int imol_map, float n_sigma);

extern "C" G_MODULE_EXPORT
void
on_pepflips_by_difference_map_dialog_response(GtkDialog       *dialog,
                                                                  gint             response_id,
                                                                  gpointer         user_data) {
   if (response_id == GTK_RESPONSE_APPLY) {
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
   }
   gtk_widget_hide(GTK_WIDGET(dialog));
}

extern "C" G_MODULE_EXPORT
void
on_ncs_differences1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *menu = widget_from_builder("ncs_differences1");
   if (menu) {
      const char *type = "ncs-diffs";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_ncs_differences1_activate\n");
   }
}

////B B FACTOR
extern "C" G_MODULE_EXPORT
void
on_temp_fact_analysis1_activate
                                        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *menu = widget_from_builder("temp_fact_analysis1");
  if (menu) {
     const char *type = "calc b factor";
     add_on_validation_graph_mol_options(menu, type);
  } else {
     std::cout << "ERROR:: failed to get menu in on_temp_fact_analysis1_activate\n";
  }

}
////E B FACTOR

extern "C" G_MODULE_EXPORT
void
on_temp_fact_variance_analysis1_activate
                                        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *menu = widget_from_builder("temp_fact_variance_analysis1");
   if (menu) {
      const char *type = "b factor";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_temp_fact_variance_analysis1_activate\n");
   }

}


extern "C" G_MODULE_EXPORT
void
on_rotamer_analysis1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *menu = widget_from_builder("rotamer_analysis1");
   if (menu) {
      const char *type = "rotamer";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_rotamer_analysis1_activate\n");
   }
}


extern "C" G_MODULE_EXPORT
void
on_density_fit_analysis1_activate      (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *menu = widget_from_builder("density_fit_analysis1");
   if (menu) {
      const char *type = "density-fit";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_density_fit1_activate\n");
   }

}


extern "C" G_MODULE_EXPORT
void
on_stereo1_activate(GMenuItem     *menuitem,
                                       gpointer         user_data) {

   //  GtkWidget *w = create_stereo_dialog();
   GtkWidget *w = widget_from_builder("stereo_dialog");
   GtkWidget *checkbutton;

   if (stereo_mode_state() == 1) { /* coot::HARDWARE_STEREO_MODE */
      checkbutton = widget_from_builder("stereo_dialog_hardware_stereo_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
   if (stereo_mode_state() == 2) { /* coot::SIDE_BY_SIDE_STEREO */
      checkbutton = widget_from_builder("stereo_dialog_side_by_side_stereo_crosseyed_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
   if (stereo_mode_state() == 4) { /* coot::SIDE_BY_SIDE_STEREO_WALL_EYE */
      checkbutton = widget_from_builder("stereo_dialog_side_by_side_stereo_walleyed_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
   if (stereo_mode_state() == 3) { /* coot::DTI_SIDE_BY_SIDE_STEREO */
      checkbutton = widget_from_builder("stereo_dialog_dti_side_by_side_stereo_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
   if (stereo_mode_state() == 5) { /* coot::ZALMAN_STEREO */
      checkbutton = widget_from_builder("stereo_dialog_zalman_stereo_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
   gtk_widget_show(w);
}



extern "C" G_MODULE_EXPORT
void
on_stereo_dialog_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("stereo_dialog");
  gtk_widget_hide(w);
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
      gtk_widget_show(nothing_bad_dialog);
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

/* Preference section */
extern "C" G_MODULE_EXPORT
void
on_preferences1_activate               (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  preferences();
}

void fill_and_show_shader_preferences() {

   GtkWidget *w  = widget_from_builder("shader_settings_dialog");
   GtkWidget *r1 = widget_from_builder("shader_settings_ssao_strength_scale");
   GtkWidget *r2 = widget_from_builder("shader_settings_ssao_radius_scale");
   GtkWidget *r3 = widget_from_builder("shader_settings_ssao_n_kernel_samples_scale");
   GtkWidget *r4 = widget_from_builder("shader_settings_shadow_strength_scale");
   GtkWidget *r5 = widget_from_builder("shader_settings_depth_blur_focus_depth_scale");
   GtkWidget *r6 = widget_from_builder("shader_settings_depth_blur_strength_scale");
   GtkWidget *r7 = widget_from_builder("shader_settings_ssao_bias_scale");
   GtkWidget *r8 = widget_from_builder("shader_settings_brightness_scale");
   GtkWidget *r9 = widget_from_builder("shader_settings_gamma_scale");

   GtkWidget *sssb_0 = widget_from_builder("shader_settings_ssao_smoothing_blur_size_0_radiobutton");
   GtkWidget *sssb_1 = widget_from_builder("shader_settings_ssao_smoothing_blur_size_1_radiobutton");
   GtkWidget *sssb_2 = widget_from_builder("shader_settings_ssao_smoothing_blur_size_2_radiobutton");

   GtkWidget *sss_1 = widget_from_builder("shader_settings_shadow_softness_1_radiobutton");
   GtkWidget *sss_2 = widget_from_builder("shader_settings_shadow_softness_2_radiobutton");
   GtkWidget *sss_3 = widget_from_builder("shader_settings_shadow_softness_3_radiobutton");

   GtkWidget *strm_1 = widget_from_builder("shader_settings_shadow_texture_resolution_multiplier_1_radiobutton");
   GtkWidget *strm_2 = widget_from_builder("shader_settings_shadow_texture_resolution_multiplier_2_radiobutton");
   GtkWidget *strm_3 = widget_from_builder("shader_settings_shadow_texture_resolution_multiplier_3_radiobutton");
   GtkWidget *strm_4 = widget_from_builder("shader_settings_shadow_texture_resolution_multiplier_4_radiobutton");
   GtkWidget *strm_5 = widget_from_builder("shader_settings_shadow_texture_resolution_multiplier_5_radiobutton");
   GtkWidget *strm_6 = widget_from_builder("shader_settings_shadow_texture_resolution_multiplier_6_radiobutton");

   GtkWidget *do_blur_checkbutton         = widget_from_builder("shader_settings_depth_blur_outline_depth_blur_radiobutton");
   GtkWidget *do_outline_checkbutton      = widget_from_builder("shader_settings_depth_blur_outline_outline_radiobutton");
   GtkWidget *do_blur_outline_checkbutton = widget_from_builder("shader_settings_depth_blur_outline_off_radiobutton");

   GtkWidget    *basic_mode_togglebutton = widget_from_builder("shader_settings_basic_mode_togglebutton");
   GtkWidget    *fancy_mode_togglebutton = widget_from_builder("shader_settings_fancy_mode_togglebutton");
   GtkWidget *standard_mode_togglebutton = widget_from_builder("shader_settings_standard_mode_togglebutton");

   GtkWidget *do_depth_fog_checkbutton = widget_from_builder("shader_settings_do_depth_fog_checkbutton");

   graphics_info_t g;

   std::cout << "fill_and_show_shader_preferences()    fancy_mode_togglebutton " << fancy_mode_togglebutton << std::endl;
   std::cout << "fill_and_show_shader_preferences() standard_mode_togglebutton " << standard_mode_togglebutton << std::endl;

   // oh dear... labels and variables inconsistent
   if (g.displayed_image_type == g.SHOW_AO_SCENE)    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(fancy_mode_togglebutton), TRUE);
   if (g.displayed_image_type == g.SHOW_BASIC_SCENE) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(standard_mode_togglebutton), TRUE);

   if (g.ssao_blur_size == 0) gtk_check_button_set_active(GTK_CHECK_BUTTON(sssb_0), TRUE);
   if (g.ssao_blur_size == 1) gtk_check_button_set_active(GTK_CHECK_BUTTON(sssb_1), TRUE);
   if (g.ssao_blur_size == 2) gtk_check_button_set_active(GTK_CHECK_BUTTON(sssb_2), TRUE);

   if (g.shadow_softness == 1) gtk_check_button_set_active(GTK_CHECK_BUTTON(sss_1), TRUE);
   if (g.shadow_softness == 2) gtk_check_button_set_active(GTK_CHECK_BUTTON(sss_2), TRUE);
   if (g.shadow_softness == 3) gtk_check_button_set_active(GTK_CHECK_BUTTON(sss_3), TRUE);

   if (g.shadow_texture_multiplier == 1) gtk_check_button_set_active(GTK_CHECK_BUTTON(strm_1), TRUE);
   if (g.shadow_texture_multiplier == 2) gtk_check_button_set_active(GTK_CHECK_BUTTON(strm_2), TRUE);
   if (g.shadow_texture_multiplier == 3) gtk_check_button_set_active(GTK_CHECK_BUTTON(strm_3), TRUE);
   if (g.shadow_texture_multiplier == 4) gtk_check_button_set_active(GTK_CHECK_BUTTON(strm_4), TRUE);
   if (g.shadow_texture_multiplier == 5) gtk_check_button_set_active(GTK_CHECK_BUTTON(strm_5), TRUE);
   if (g.shadow_texture_multiplier == 6) gtk_check_button_set_active(GTK_CHECK_BUTTON(strm_6), TRUE);

   if (! g.shader_do_outline_flag && !g.shader_do_depth_of_field_blur_flag)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(do_blur_outline_checkbutton), TRUE);

   if (g.shader_do_depth_of_field_blur_flag) // not shader_do_depth_blur_flag (what's that used for? - delete it)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(do_blur_checkbutton), TRUE);

   if (g.shader_do_outline_flag)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(do_outline_checkbutton), TRUE);

   if (graphics_info_t::shader_do_depth_fog_flag)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(do_depth_fog_checkbutton), TRUE);
   else
      gtk_check_button_set_active(GTK_CHECK_BUTTON(do_depth_fog_checkbutton), FALSE);

   // make this insensitve if mode is not fancy
   GtkWidget *fancy_vbox1 = widget_from_builder("shader_settings_fancy_vbox1");
   GtkWidget *fancy_vbox2 = widget_from_builder("shader_settings_fancy_vbox2");
   std::cout << "fill_and_show_shader_preferences() fancy_vbox1 " << fancy_vbox1 << std::endl;
   bool is_fancy_mode = true;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(basic_mode_togglebutton)))    is_fancy_mode = false;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(standard_mode_togglebutton))) is_fancy_mode = false;
   if (! is_fancy_mode) {
      gtk_widget_set_sensitive(fancy_vbox1, FALSE);
      gtk_widget_set_sensitive(fancy_vbox2, FALSE);
   }

   double v1 = graphics_info_t::ssao_strength;
   double v2 = graphics_info_t::SSAO_radius;
   double v3 = graphics_info_t::n_ssao_kernel_samples;
   double v4 = graphics_info_t::shadow_strength;
   double v5 = graphics_info_t::focus_blur_z_depth;
   double v6 = graphics_info_t::focus_blur_strength;
   double v7 = graphics_info_t::SSAO_bias;
   double v8 = graphics_info_t::effects_brightness;
   double v9 = graphics_info_t::effects_gamma;

   gtk_range_set_range(GTK_RANGE(r1), 0.0, 2.0);
   gtk_range_set_value(GTK_RANGE(r1), v1);
   gtk_range_set_range(GTK_RANGE(r2), 0.0, 100.0);
   gtk_range_set_value(GTK_RANGE(r2), v2);
   gtk_range_set_range(GTK_RANGE(r3), 0.0, 256.0);
   gtk_range_set_value(GTK_RANGE(r3), v3);
   gtk_range_set_range(GTK_RANGE(r4), 0.0, 1.0);
   gtk_range_set_value(GTK_RANGE(r4), v4);
   gtk_range_set_range(GTK_RANGE(r5), 0.0, 1.0);
   gtk_range_set_value(GTK_RANGE(r5), v5);
   gtk_range_set_range(GTK_RANGE(r6), 0.0, 6.0);
   gtk_range_set_value(GTK_RANGE(r6), v6);
   gtk_range_set_range(GTK_RANGE(r7), 0.0, 0.4);
   gtk_range_set_value(GTK_RANGE(r7), v7);
   gtk_range_set_range(GTK_RANGE(r8), 0.0, 3.0);
   gtk_range_set_value(GTK_RANGE(r8), v8);
   gtk_range_set_range(GTK_RANGE(r9), 0.0, 2.0);
   gtk_range_set_value(GTK_RANGE(r9), v9);

   gtk_widget_show(w);

}

extern "C" G_MODULE_EXPORT
void
on_shader_preferences_activate (GMenuItem *menuitem,
                                gpointer   user_data) {
   fill_and_show_shader_preferences();
}


extern "C" G_MODULE_EXPORT
gboolean
on_shader_settings_dialog_delete_event(GtkWidget       *widget,
                                       GdkEvent        *event,
                                       gpointer         user_data) {

   // this happens when the window manager closes the window (or tries to)
   // std::cout << "-------------- on_shader_settings_dialog_delete_event() " << std::endl;
   gtk_widget_hide(widget);
   return TRUE;
}

extern "C" G_MODULE_EXPORT
gboolean
on_shader_settings_dialog_destroy_event(GtkWidget       *widget,
                                                           GdkEvent        *event,
                                                           gpointer         user_data) {

   // this doesn't happen
   std::cout << "-------------- on_shader_settings_dialog_destroy_event() " << std::endl;
   return TRUE;
}

extern "C" G_MODULE_EXPORT
gboolean
on_shader_settings_dialog_destroy(GtkWidget       *widget,
                                                           gpointer         user_data) {

   // this doesn't happen
   std::cout << "-------------- on_shader_settings_dialog_ destroy " << std::endl;
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
on_diff_map_peaks_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("diff_map_peaks_dialog");
   clear_diff_map_peaks();
   gtk_widget_hide(dialog);
}

extern "C" G_MODULE_EXPORT
void
on_diff_map_peaks_dialog_update_button_clicked(GtkButton       *button,
                                                                   gpointer         user_data) {
   graphics_info_t g;
   g.fill_difference_map_peaks_button_box(true); // force fill.
   
}

extern "C" G_MODULE_EXPORT
void
on_generate_diff_map_peaks_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("generate_diff_map_peaks_dialog");
   difference_map_peaks_by_widget(w); // make the results (and show the results dialog)
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_generate_diff_map_peaks_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("generate_diff_map_peaks_dialog");
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_difference_map_peaks1_activate (GMenuItem     *menuitem,
                                                       gpointer         user_data)
{
  GtkWidget *w = wrapped_create_generate_diff_map_peaks_dialog();
  gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_superpose_reference_chain_checkbutton_toggled(GtkToggleButton *togglebutton,
                                                                     gpointer user_data) {

  GtkWidget *combobox = widget_from_builder("superpose_dialog_reference_chain_combobox");
  if (gtk_toggle_button_get_active(togglebutton)) {
    gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
    printf("calling fill_superpose_combobox_with_chain_options()\n");
    fill_superpose_combobox_with_chain_options(combobox, 1);
    printf("done fill_superpose_combobox_with_chain_options()\n");
  } else {
    gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
  }

  printf("done on_superpose_reference_chain_checkbutton_toggled()\n");

}


extern "C" G_MODULE_EXPORT
void
on_superpose_moving_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *combobox = widget_from_builder("superpose_dialog_moving_chain_combobox");

  if (gtk_toggle_button_get_active(togglebutton)) {
    fill_superpose_combobox_with_chain_options(combobox, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
  } else {
    gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
  }
}

extern "C" G_MODULE_EXPORT
void
on_draw_fullscreen_activate(GMenuItem *menuitem,
                                                gpointer     user_data) {

   fullscreen();
   set_show_modelling_toolbar(0);
}


extern "C" G_MODULE_EXPORT
void
on_draw_ncs_ghosts_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
/* Function no longer used.  Kept in glade (not visible) for historical reasons

   GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   if (gtk_toggle_button_get_active(togglebutton)) {
      printf("yes radiobutton toggled on.\n");
      make_ncs_ghosts_maybe(w);
   }
*/
}


extern "C" G_MODULE_EXPORT
void
on_draw_ncs_ghosts_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_ncs_maps_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ncs_maps_dialog");
   make_dynamically_transformed_ncs_maps_by_widget(w);
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_ncs_maps_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ncs_maps_dialog");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_ncs_maps1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_ncs_maps_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_rotamer_selection_dialog_destroy    (GtkWidget       *object,
                                        gpointer         user_data)
{
   /* set the dialog
      store_window_position(COOT_ROTAMER_SELECTION_DIALOG, dialog); */
   set_graphics_rotamer_dialog(NULL);
}


extern "C" G_MODULE_EXPORT
void
on_pointer_distances_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
      printf("pointer distances toggle button toggled on\n");
   } else {
      printf("pointer distances toggle button toggled off\n");
   }
   toggle_pointer_distances_show_distances(togglebutton);
}


extern "C" G_MODULE_EXPORT
void
on_pointer_distances_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("pointer_distances_dialog");
   execute_pointer_distances_settings(dialog);
   gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_pointer_distances1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *w = create_pointer_distances_dialog();
   GtkWidget *w = widget_from_builder("pointer_distances_dialog");
   fill_pointer_distances_widget(w);
   gtk_widget_show(w);
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
      gtk_widget_hide(dialog);

}


extern "C" G_MODULE_EXPORT
void
on_align_and_mutate_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = widget_from_builder("align_and_mutate_dialog");
   gtk_widget_hide(dialog);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_differences_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ramachandran_plot_differences_dialog");
   int istat = do_ramachandran_plot_differences_by_widget(w);
   if (istat) 			/* the plot was drawn (i.e. no chain selection funnies) */
      gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_plot_differences_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ramachandran_plot_differences_dialog");
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_ramachandran_differences_plot1_activate
                                        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_ramachandran_plot_differences_dialog();
   gtk_widget_show(w);

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
on_check_waters1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   /* we do this (show the map selection widget) here because we want
	  it to appear after the water check dialog so that it appears
	  on top rather than underneath. */
   GtkWidget *w = wrapped_create_check_waters_dialog();
   int imol_map = imol_refinement_map();
   gtk_widget_show(w);
   if (imol_map < 0)
      show_select_map_dialog();

}



extern "C" G_MODULE_EXPORT
void
on_checked_waters_baddies_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("checked_waters_baddies_dialog");
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_align_and_mutate1_activate            (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_align_and_mutate_dialog();
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_about_other_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_delete_item_keep_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
}


#if 0
extern "C" G_MODULE_EXPORT
void
on_delete_item_dialog_destroy          (GtkWidget       *object,
                                       gpointer         user_data) {

   // this function should not be linked even

   clear_pending_delete_item();
   clear_pending_picks();
   normal_cursor();
   store_delete_item_widget(NULL);
}
#endif


extern "C" G_MODULE_EXPORT
void
on_save_state1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
/*    save_state(); old inteface - before DK. */

   GtkWidget *file_chooser = coot_save_state_chooser();

   gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(file_chooser), save_state_file_name_raw());
   add_filename_filter_button(file_chooser, COOT_SCRIPTS_FILE_SELECTION);
   /* add_ccp4i_project_optionmenu(fileselection, COOT_SCRIPTS_FILE_SELECTION); */
   set_file_selection_dialog_size(file_chooser);
   gtk_widget_show(file_chooser);
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
on_main_window_statusbar_text_popped   (GtkStatusbar    *statusbar,
                                        guint            context_id,
                                        gchar           *text,
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
on_fit_loop1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

/*    GtkWidget *w = wrapped_fit_loop_dialog(); */
/*    gtk_widget_show(w); */
}

#include "fit-loop-gui.hh"

extern "C" G_MODULE_EXPORT
void
on_fit_loop_by_rama_search1_activate (GMenuItem     *menuitem,
                                                          gpointer         user_data)
{

   GtkWidget *w = create_fit_loop_rama_search_dialog_gtkbuilder_version(); // fixed
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_fit_loop_by_database_search1_activate (GMenuItem     *menuitem,
					  gpointer         user_data)
{
  wrapped_fit_loop_db_loop_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_fit_loop_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("mutate_sequence_dialog");
   fit_loop_using_dialog();
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_reset_view1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   reset_view();
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_A_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_hide(w);
   do_base_mutation("A");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_C_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_hide(w);
   do_base_mutation("C");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_G_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_hide(w);
   do_base_mutation("G");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_T_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_hide(w);
   do_base_mutation("T");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_U_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_hide(w);
   do_base_mutation("U");
}


extern "C" G_MODULE_EXPORT
void
on_base_chooser_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("nucleic_acid_base_chooser_dialog");
   gtk_widget_hide(w);
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
on_change_chain_ids2_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_change_chain_id_dialog(); // uses builder
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_change_chains_rechain_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   // GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   change_chain_id_by_widget(w);
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_change_chain_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   // gtk_widget_hide(w);

   GtkWidget *w = widget_from_builder("change_chain_id_dialog");
   gtk_widget_hide(w);

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
on_on_line_docs_url1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *w = create_doc_urls_dialog();
   GtkWidget *w = widget_from_builder("doc_urls_dialog");
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_on_line_documentation_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("doc_urls_dialog");
   gtk_widget_hide(w);

}




extern "C" G_MODULE_EXPORT
void
on_save_state_cancel_button1_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder(
				"save_state_fileselection");
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_change_chain_residue_range_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{


}


extern "C" G_MODULE_EXPORT
void
on_change_chain_residue_range_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

   // GtkWidget *hbox = widget_from_builder("change_chain_id_residue_range_hbox");
   GtkWidget *hbox = widget_from_builder("change_chain_id_residue_range_hbox");

   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);
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
	 show_select_map_dialog();
	 info_dialog("A map has not yet been assigned for Refinement/Fitting");
      }
   }
}

extern "C" G_MODULE_EXPORT
void
on_mutate_sequence_use_ramachandran_restraints_checkbutton_toggled(GtkToggleButton *togglebutton,
                                                                   gpointer         user_data) {
   /* not doing anything because the button state read at execution time  */
}

extern "C" G_MODULE_EXPORT
void
on_check_waters_b_factor_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // GtkWidget *hbox = widget_from_builder("check_waters_b_factor_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_b_factor_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);
}




extern "C" G_MODULE_EXPORT
void
on_check_waters_min_dist_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = widget_from_builder(
				   "check_waters_min_dist_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_max_dist_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // GtkWidget *hbox = widget_from_builder("check_waters_max_dist_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_max_dist_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}


extern "C" G_MODULE_EXPORT
void
on_check_waters_map_sigma_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // GtkWidget *hbox = widget_from_builder("check_waters_sigma_level_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_sigma_level_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_check_waters_by_difference_map_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // GtkWidget *hbox = widget_from_builder("check_waters_by_difference_map_hbox");
   GtkWidget *hbox = widget_from_builder("check_waters_by_difference_map_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}




extern "C" G_MODULE_EXPORT
void
on_residue_info_occ_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   // GtkWidget *entry = widget_from_builder("residue_info_master_atom_occ_entry");
   // GtkWidget *alt_conf_checkbutton = widget_from_builder("residue_info_occ_apply_to_altconf_checkbutton");

   GtkWidget *entry = widget_from_builder("residue_info_master_atom_occ_entry");
   GtkWidget *alt_conf_checkbutton = widget_from_builder("residue_info_occ_apply_to_altconf_checkbutton");

   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(entry, TRUE);
   } else {
      if (! gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(alt_conf_checkbutton)))
         gtk_widget_set_sensitive(entry, FALSE);
   }
}


extern "C" G_MODULE_EXPORT
void
on_residue_info_occ_apply_to_altconf_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   // GtkWidget *occ_entry = widget_from_builder("residue_info_master_atom_occ_entry");
   // GtkWidget *alt_conf_entry = widget_from_builder("residue_info_occ_apply_to_alt_conf_entry");

   GtkWidget *occ_entry      = widget_from_builder("residue_info_master_atom_occ_entry");
   GtkWidget *alt_conf_entry = widget_from_builder("residue_info_occ_apply_to_alt_conf_entry");

   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(occ_entry,      TRUE);
      gtk_widget_set_sensitive(alt_conf_entry, TRUE);
   } else {
      gtk_widget_set_sensitive(occ_entry,      FALSE);
      gtk_widget_set_sensitive(alt_conf_entry, FALSE);
   }

}




extern "C" G_MODULE_EXPORT
void
on_residue_info_b_factor_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   // GtkWidget *entry = widget_from_builder("residue_info_master_atom_b_factor_entry");

   GtkWidget *entry = widget_from_builder("residue_info_master_atom_b_factor_entry");

  if (gtk_toggle_button_get_active(togglebutton))
    gtk_widget_set_sensitive(entry, TRUE);
  else
    gtk_widget_set_sensitive(entry, FALSE);
}



extern "C" G_MODULE_EXPORT
void
on_other_modelling_tools_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("other_model_tools_dialog");
   GtkWidget *w = widget_from_builder("other_model_tools_dialog");

   // gtk_widget_hide(w);
   gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_other_modelling_tools1_activate     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_other_model_tools_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_cis_trans_conversion_toggle_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
      do_cis_trans_conversion_setup(1);
   else
      do_cis_trans_conversion_setup(0);

}


extern "C" G_MODULE_EXPORT
void
on_other_model_tools_dialog_destroy    (GtkWidget       *object,
                                        gpointer         user_data)
{
   do_cis_trans_conversion_setup(0);
   unset_other_modelling_tools_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_undo_last_navigation1_activate      (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   undo_last_move();
}


extern "C" G_MODULE_EXPORT
void
on_undo_symmetry_view1_activate (GMenuItem     *menuitem,
                                                     gpointer         user_data)
{
   undo_symmetry_view();
}


extern "C" G_MODULE_EXPORT
void
on_get_pdb_and_map_using_eds1_activate (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   int n = COOT_ACCESSION_CODE_WINDOW_EDS;
   // GtkWidget *window = create_accession_code_window();
   GtkWidget *window = widget_from_builder("accession_code_window");
   g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
   gtk_widget_show(window);
}


extern "C" G_MODULE_EXPORT
void
on_ligand_big_blob_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data)
{
  free_blob_dialog_memory(GTK_WIDGET(object));
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
on_simple1_activate(GMenuItem     *menuitem,
                                        gpointer         user_data) {

   GtkWidget *file_chooser = coot_screendump_chooser();
   set_transient_and_position(COOT_UNDEFINED_WINDOW, file_chooser);
   /*    set_directory_for_fileselection(fileselection); */
   /*    set_filename_for_filechooserselection(fileselection, "coot.png"); */
   /*    set_file_selection_dialog_size(fileselection); */

   g_object_set_data(G_OBJECT(file_chooser), "image_type", GINT_TO_POINTER(COOT_SCREENDUMP_SIMPLE));

   gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(file_chooser), "coot-screendump.tga");

   gtk_widget_show(file_chooser);

   check_for_dark_blue_density(); /* give a dialog if density it too dark (blue) */

}


extern "C" G_MODULE_EXPORT
void
on_povray1_activate                    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *fileselection = coot_screendump_chooser();
   set_transient_and_position(COOT_UNDEFINED_WINDOW, fileselection);
   g_object_set_data(G_OBJECT(fileselection), "mode", GINT_TO_POINTER(COOT_SCREENDUMP_POVRAY));
   /*    set_directory_for_fileselection(fileselection); */
   gtk_widget_show(fileselection);
}


extern "C" G_MODULE_EXPORT
void
on_raster3d1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *file_chooser = coot_screendump_chooser();
   set_transient_and_position(COOT_UNDEFINED_WINDOW, file_chooser);
   g_object_set_data(G_OBJECT(file_chooser), "image_type", GINT_TO_POINTER(COOT_SCREENDUMP_RASTER3D));
   gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(file_chooser), "coot-screendump.png.r3d");
   gtk_widget_show(file_chooser);
}


extern "C" G_MODULE_EXPORT
void
on_screendump_image_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *filechooser = widget_from_builder("screendump_filechooser");
  /*

    To restore this then I need to fix up the usage of gtk_object_get_user_data()
    i.e. is the user data set properly. This is a merge issue, fix later.

   int image_type = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(filechooser)));
   const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(filechooser));

   if (image_type == COOT_SCREENDUMP_SIMPLE) {
      screendump_image(filename);
   }
   if (image_type == COOT_SCREENDUMP_POVRAY) {
      make_image_povray(filename);
   }
   if (image_type == COOT_SCREENDUMP_RASTER3D) {
      make_image_raster3d(filename);
   }
  */
   if (filechooser)
     gtk_widget_hide(filechooser);
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
     gtk_widget_hide(filechooser);

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
on_other_tools_place_strand_here_button_clicked     (GtkButton       *button,
                                        gpointer         user_data) {
  place_strand_here_dialog(); 	/* choose the python version in there, if needed. */
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
  gtk_widget_hide(w);

}





extern "C" G_MODULE_EXPORT
void
on_show_symmetry_molecule_control_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = symmetry_molecule_controller_dialog();
  gtk_widget_show(w);
}



extern "C" G_MODULE_EXPORT
void
on_ncs_control_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("ncs_control_dialog");
   gtk_widget_hide(w);
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
  gtk_widget_hide(w);
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
on_plane_distances1_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_lsq_plane_dialog();
   setup_lsq_deviation(0);
   setup_lsq_plane_define(1);
   gtk_widget_show(w);
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
on_ncs_ghost_control1_activate         (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_ncs_control_dialog(); // uses builder
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
gboolean
on_ncs_control_dialog_delete_event(GtkWidget       *widget,
                                                       GdkEvent        *event,
                                                       gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}


extern "C" G_MODULE_EXPORT
gboolean
on_coords_colour_control_dialog_delete_event(GtkWidget       *widget,
                                                                 GdkEvent        *event,
                                                                 gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}


extern "C" G_MODULE_EXPORT
void
on_coords_colour_control_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
   // do nothing!
}


extern "C" G_MODULE_EXPORT
void
on_coord_colour_control_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  GtkWidget *w = widget_from_builder("coords_colour_control_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_bond_colours1_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   // note on_coord_colour_control_ok_button_clicked()
   // (above)
   GtkWidget *w = widget_from_builder("coords_colour_control_dialog");
   graphics_info_t g;
   g.fill_bond_colours_dialog_internal(w);
   gtk_widget_show(w);

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


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_coot_online_doc_search_entry_key_press_event (GtkWidget       *widget,
                                                 GdkEventKey     *event,
                                                 gpointer         user_data)
{

  const char *text;

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
     text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(widget)));
     handle_online_coot_search_request(text);
  }

  return FALSE;
}
#endif


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
on_smiles1_activate                    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_smiles_gui();
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
on_generic_display_objects1_activate   (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

  generic_objects_gui_wrapper();

}


extern "C" G_MODULE_EXPORT
void on_generic_objects_dialog_response(GtkDialog       *dialog,
                                                            gint             response_id,
                                                            gpointer         user_data) {

    if (response_id == GTK_RESPONSE_CLOSE) {
       gtk_widget_hide(GTK_WIDGET(dialog));
    }
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
// what a terrible function name!
extern "C" G_MODULE_EXPORT
gboolean
on_entry1_key_press_event(GtkWidget       *widget,
                          GdkEventKey     *event,
                          gpointer         user_data) {

   // 1: I don't like the name of this callback.
   // 2: There is no EM version of this.

   GtkEntry *entry = (GTK_ENTRY(widget_from_builder("map_parameters_x_ray_radius_entry")));
   const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      set_density_size_from_widget(text);
   }
   return FALSE;
}
#endif

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
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  /* not visible */
  printf("helix togglebutton toggled - ignored\n");
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_beta_strand_peptide_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  /* not visible */
  printf("beta strand togglebutton toggled - ignored\n");

}


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_ramachandran_goodness_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int state = 0;
  if (gtk_toggle_button_get_active(togglebutton)) {
    state = 1;
  }
  set_refine_ramachandran_angles(state);
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_use_peptide_omegas_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    add_omega_torsion_restriants();
  } else {
    remove_omega_torsion_restriants();
  }
}


extern "C" G_MODULE_EXPORT
void
on_probe_clashes1_activate             (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   const char *type = "probe";
   GtkWidget *menu = widget_from_builder("probe_clashes1");
   if (menu) {
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_probe_clashes1_activate\n");
   }

}


extern "C" G_MODULE_EXPORT
void
on_validate1_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
#if 0 // 20211002-PE I no longer wish to do this
   GtkWidget *menu_item = 0;

   if (probe_available_p() == 0) { /* no */
      menu_item = widget_from_builder("probe_clashes1");
      if (!menu_item) {
         printf("Failed to get probe_clashes1 menu item :-(\n");
      } else {
         /* desensitize it */
         gtk_widget_set_sensitive(menu_item, FALSE);
      }
   }
#endif
}


extern "C" G_MODULE_EXPORT
void
on_spin_view_on_off1_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   std::cout << "found on_spin_view_on_off1_activate() " << std::endl;
   toggle_idle_spin_function();
}


extern "C" G_MODULE_EXPORT
void
on_other_tools_RNA_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_nucleotide_builder_dialog();
  gtk_widget_show(w);
}

extern "C" G_MODULE_EXPORT
void
on_other_tools_base_pair_toggle_button_toggled      (GtkToggleButton       *button,
                                        gpointer         user_data) {
  if (gtk_toggle_button_get_active(button))
    setup_base_pairing(1);
  else
    setup_base_pairing(0);
}



extern "C" G_MODULE_EXPORT
void
on_ideal_rna_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("nucleotide_builder_dialog");
  ideal_nucleic_acid_by_widget(w);
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_ideal_rna_cancel_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("nucleotide_builder_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_unit_cell_yes_radiobutton_toggled   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   if (gtk_toggle_button_get_active(togglebutton))
      set_show_unit_cells_all(1);
  else
      set_show_unit_cells_all(0);
}


extern "C" G_MODULE_EXPORT
void
on_unit_cell_no_radiobutton_toggled    (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_show_unit_cells_all(0);
  else
      set_show_unit_cells_all(1);
}




extern "C" G_MODULE_EXPORT
void
on_move_molecule_here1_activate        (GMenuItem     *menuitem,
                                                            gpointer         user_data) {

   // GtkWidget *w = wrapped_create_move_molecule_here_dialog();
   GtkWidget *w = widget_from_builder("move_molecule_here_dialog");
   fill_move_molecule_here_dialog(w);
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_move_molecule_here_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  move_molecule_here_by_widget();
  GtkWidget *w = widget_from_builder("move_molecule_here_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_move_molecule_here_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("move_molecule_here_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_monomer_library_search_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("monomer_search_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_search_monomer_library1_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // GtkWidget *w = create_monomer_search_dialog();
   GtkWidget *w = widget_from_builder("monomer_search_dialog");
   gtk_widget_show(w);
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
on_monomer_search_entry_changed (GtkEditable     *editable,
                                 gpointer         user_data) {

}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_monomer_search_entry_key_press_event (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data) {

   GtkWidget *entry = widget_from_builder("monomer_search_entry");
   entry_char_type *text;
   GtkWidget *viewport = widget_from_builder("monomer_search_results_viewport");

   if (entry) {
      if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
         text = gtk_editable_get_text(GTK_EDITABLE(entry));
         if (text) {
            handle_make_monomer_search(text, viewport);
         }
      }
   }
   return FALSE;
}
#endif


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
on_least_squares_close_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("least_squares_dialog");
  update_lsq_dialog_store_values(w);
  gtk_widget_hide(w);
}




extern "C" G_MODULE_EXPORT
void
on_least_squares_cancel_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("least_squares_dialog");
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_lsq_superpose1_activate             (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_least_squares_dialog(); // uses builder
   gtk_widget_show(w);

}


#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

   // no GtkStateType
#else

extern "C" G_MODULE_EXPORT
void
on_goto_atom_window_state_changed      (GtkWidget       *widget,
                                                            GtkStateType     state,
                                                            gpointer         user_data)
{
  /* how does this window get a state changed event?  Not at all, as
     far as I can see.  Nothing happens here. */
  printf("Go To Atom window state changed!\n");
}
#endif


#if (GTK_MAJOR_VERSION >= 4)
// 20220601-PE  I don't know what GdkEventConfigure has become (if anything)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_goto_atom_window_configure_event    (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{

  store_window_position(COOT_GO_TO_ATOM_WINDOW, widget);
  return FALSE;
}
#endif


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_model_refine_dialog_configure_event (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{

  store_window_position(COOT_MODEL_REFINE_DIALOG, widget);
  return FALSE;
}
#endif


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_display_control_window_glade_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_DISPLAY_CONTROL_WINDOW, widget);
  return FALSE;
}
#endif

// void
//nuser_function (GMenuItem *menuitem,
//               gpointer     user_data)

extern "C" G_MODULE_EXPORT
void
on_rotate_translate_item_residue_range_item_activate(GMenuItem *menuitem,
                                                     gpointer user_data) {

   std::cout << "on_rotate_translate_item_residue_range_item_activate() " << std::endl;
   set_rot_trans_object_type(ROT_TRANS_TYPE_ZONE);
   do_rot_trans_setup(1);
}

extern "C" G_MODULE_EXPORT
void
on_rotate_translate_item_residue_item_activate(GMenuItem *menuitem,
                                               gpointer user_data) {

   std::cout << "on_rotate_translate_item_residue_item_activate" << std::endl;
   set_rot_trans_object_type(ROT_TRANS_TYPE_RESIDUE);
   do_rot_trans_setup(1);
}


extern "C" G_MODULE_EXPORT
void
on_rotate_translate_item_molecule_item_activate(GMenuItem *menuitem,
                                                gpointer user_data) {

   std::cout << "on_rotate_translate_item_molecule_item_activate" << std::endl;
   set_rot_trans_object_type(ROT_TRANS_TYPE_MOLECULE);
   do_rot_trans_setup(1);
}

extern "C" G_MODULE_EXPORT
void
on_rotate_translate_item_chain_item_activate(GMenuItem *menuitem,
                                             gpointer user_data) {

   std::cout << "on_rotate_translate_item_chain_item_activate()" << std::endl;
   set_rot_trans_object_type(ROT_TRANS_TYPE_CHAIN);
   do_rot_trans_setup(1);
}


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_delete_item_dialog_configure_event  (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_DELETE_WINDOW, widget);
  return FALSE;
}
#endif



#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_rotate_translate_obj_dialog_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_ROTATE_TRANSLATE_DIALOG, widget);
  return FALSE;
}
#endif



#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_accept_reject_refinement_dialog_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_ACCEPT_REJECT_WINDOW, widget);
  return FALSE;
}
#endif


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_dynarama_window_configure_event     (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  resize_rama_canvas(widget, event);
   store_window_position(COOT_RAMACHANDRAN_PLOT_WINDOW, widget);
   return FALSE;
}
#endif


#if (GTK_MAJOR_VERSION >= 4)
// delete-event is no longer a thing, I think
#else
extern "C" G_MODULE_EXPORT
gboolean
on_screendump_filechooser_dialog_delete_event(GtkWidget       *widget,
                                              GdkEvent        *event,
                                              gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}
#endif



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
                                        gpointer         user_data)
{
  show_set_undo_molecule_chooser();
}


extern "C" G_MODULE_EXPORT
void
on_gln_and_asn_b_factor_outliers1_activate
                                        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *menu = widget_from_builder("gln_and_asn_b_factor_outliers1");
   if (menu) {
      const char *type = "gln_and_asn_b_factor_outliers";
      add_on_validation_graph_mol_options(menu, type);
   } else {
      printf("failed to get menu in on_gln_and_asn_b_factor_outliers1_activate\n");
   }
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
on_refine_params_dialog_destroy        (GtkWidget       *object,
                                        gpointer         user_data)
{
  unset_refine_params_dialog();
}


extern "C" G_MODULE_EXPORT
void
on_update_go_to_atom_from_current_position_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   update_go_to_atom_from_current_position();
}


/* gtk2 extras */
extern "C" G_MODULE_EXPORT
void
on_rz_simple_activate                  (GMenuItem     *menuitem,
                                        gpointer         user_data) {
}

extern "C" G_MODULE_EXPORT
void
on_rz_start_multizone1_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data) {
}

extern "C" G_MODULE_EXPORT
void
on_rz_end_multizone_activate           (GMenuItem     *menuitem,
                                        gpointer         user_data){

}

extern "C" G_MODULE_EXPORT
void
on_display_manager_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window(); // uses widget_from_builder()
   gtk_widget_show(widget);
}


extern "C" G_MODULE_EXPORT
void
on_reset_view_button_clicked           (GtkButton       *button,
                                        gpointer         user_data) {
   reset_view();
}

extern "C" G_MODULE_EXPORT
void
on_display_manager_toolbutton_clicked  (GtkButton   *toolbutton,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
}

extern "C" G_MODULE_EXPORT
void
on_reset_view_toolbutton_clicked       (GtkButton   *toolbutton,
                                        gpointer         user_data)
{
   reset_view();
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
on_single_map_properties_contour_level_apply_button_clicked (GtkButton       *button,
                                                                                 gpointer         user_data)
{

   //single_map_properties_apply_contour_level_to_map(w); /* check now
   //							  made here
   //							  for valid
   //							  map
   //							  molecule. */

   int imol    = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "imol"));
   GtkWidget *entry = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "contour_level_entry"));
   GtkWidget *togglebutton = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "single_map_properties_absolute_radiobutton"));

   std::cout << "imol: " << imol << std::endl;
   std::cout << "entry " << entry << std::endl;
   std::cout << "togglebutton " << togglebutton << std::endl;

   if (is_valid_map_molecule(imol)) {
      std::string t = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      try {
         float f = coot::util::string_to_float(t);
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
            set_contour_level_in_sigma(imol, f);
         } else {
            set_contour_level_absolute(imol, f);
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "Failed to interpret " << t << std::endl;
      }
      
   }
   
}

extern "C" G_MODULE_EXPORT
void
on_display_map_style_as_lines_radiobutton_toggled (GtkToggleButton *togglebutton,
                                                                       gpointer         user_data) {

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton), "imol"));
   std::cout << "on_display_map_style_as_lines_radiobutton_toggled() imol " << imol << std::endl;
   if (gtk_toggle_button_get_active(togglebutton)) {
      set_draw_map_standard_lines(imol, 1);
      set_draw_solid_density_surface(imol, 0);
   } else {
      set_draw_map_standard_lines(imol, 0);
      set_draw_solid_density_surface(imol, 1);
   }
}

extern "C" G_MODULE_EXPORT
void
on_display_map_style_surface_radiobutton_toggled(GtkToggleButton *togglebutton,
                                                                     gpointer         user_data) {

   // we don't need to do anything because it's all handled by the other callback
}



extern "C" G_MODULE_EXPORT
void
on_checked_waters_baddies_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
  store_checked_waters_baddies_dialog(NULL);
}


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_style_changed (GtkToolbar      *toolbar,
                                GtkToolbarStyle  style,
                                gpointer         user_data) {
  /* this does not do anything and doesnt need to */

   // 20220601-PE so delete whatever calls this function (if anything now)
}
#endif


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_model_toolbar_button_press_event    (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data)
{

   if (event->type == GDK_BUTTON_PRESS && event->button == 3) {
      GdkEventButton *eb = static_cast<GdkEventButton *>(user_data);
      toolbar_popup_menu(GTK_TOOLBAR(widget), eb, event);
      return TRUE;
   }

   return FALSE;
}
#endif




extern "C" G_MODULE_EXPORT
void
on_model_toolbar_icons_and_text1_activate (GMenuItem     *menuitem,
                                           gpointer         user_data) {
  /*
  GtkWidget *toolbar = widget_from_builder("model_toolbar");
  if (GTK_CHECK_MENU_ITEM(menuitem)->active){
      gtk_toolbar_set_style (GTK_TOOLBAR (toolbar), GTK_TOOLBAR_BOTH_HORIZ);
      GtkWidget *button;
      button = widget_from_builder("model_toolbar_refine_control_button");
      gtk_button_set_label(GTK_BUTTON(button), "Refine/Regularize Control...");
      button = widget_from_builder("model_toolbar_select_map_button");
      gtk_button_set_label(GTK_BUTTON(button), "Select Map...");
  }
  */

}

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_icons1_activate       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  /*
  GtkWidget *toolbar = widget_from_builder("model_toolbar");
  if (GTK_CHECK_MENU_ITEM(menuitem)->active){
    gtk_toolbar_set_style (GTK_TOOLBAR (toolbar), GTK_TOOLBAR_ICONS);
    GtkWidget *button;
    button = widget_from_builder("model_toolbar_refine_control_button");
    gtk_button_set_label(GTK_BUTTON(button), "R/RC");
    button = widget_from_builder("model_toolbar_select_map_button");
    gtk_button_set_label(GTK_BUTTON(button), "Map");
  }
  */
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_text1_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  /*
  GtkWidget *toolbar = widget_from_builder("model_toolbar");
  if (GTK_CHECK_MENU_ITEM(menuitem)->active){
    gtk_toolbar_set_style (GTK_TOOLBAR (toolbar), GTK_TOOLBAR_TEXT);
    GtkWidget *button;
    button = widget_from_builder("model_toolbar_refine_control_button");
    gtk_button_set_label(GTK_BUTTON(button), "Refine/Regularize Control...");
    button = widget_from_builder("model_toolbar_select_map_button");
    gtk_button_set_label(GTK_BUTTON(button), "Select Map...");
  }
  */
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_main_icons_activate   (GMenuItem     *menuitem,
                                        gpointer         user_data) {

#if (GTK_MAJOR_VERSION >= 4)
#else
   if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
      show_model_toolbar_main_icons();
   }
#endif

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_all_icons_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
#if (GTK_MAJOR_VERSION >= 4)
#else
   if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
      show_model_toolbar_all_icons();
   }
#endif
}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_user_defined1_activate(GMenuItem     *menuitem,
                                        gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_model_toolbar_setting1_activate(GMenuItem     *menuitem,
                                   gpointer         user_data)
{
  update_model_toolbar_icons_menu();
}


#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_menutoolbutton1_show_menu (GtkMenuToolButton *menutoolbutton,
                                            gpointer         user_data) {
  /* I dont think anything needs to happen here */
   // 20220601-PE delete this function
}
#endif

extern "C" G_MODULE_EXPORT
void
on_model_toolbar_display_manager_togglebutton_toggled(GtkToggleButton *toggletoolbutton,
                                                      gpointer         user_data) {

   if (gtk_toggle_button_get_active(toggletoolbutton)) {
      g_print("BL DEBUG:: display menu toggled");
   }
}


extern "C" G_MODULE_EXPORT
void
on_toolbar_display_manager_maps_all_activate(GMenuItem     *menuitem,
                                             gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_toolbar_display_manager_molecules_all_activate(GMenuItem     *menuitem,
                                                  gpointer         user_data)
{

}

extern "C" G_MODULE_EXPORT
void
on_calculate_scripting_python1_activate (GMenuItem     *menuitem,
                                         gpointer         user_data) {

  post_python_scripting_window();

}



extern "C" G_MODULE_EXPORT
void
on_calculate_scripting_scheme1_activate (GMenuItem     *menuitem,
                                         gpointer         user_data) {

   post_scheme_scripting_window();
}

extern "C" G_MODULE_EXPORT
void
on_aboutdialog_close                   (GtkDialog       *dialog,
                                        gpointer         user_data)
{
/* not this callback... */
}

extern "C" G_MODULE_EXPORT
void
on_aboutdialog_response (GtkDialog       *dialog,
                         gint             response_id,
                         gpointer         user_data) {
  gtk_widget_hide(GTK_WIDGET(dialog));
}

extern "C" G_MODULE_EXPORT
void
on_check_waters_check1_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data) {
}

extern "C" G_MODULE_EXPORT
void
on_check_waters_delete1_activate       (GMenuItem     *menuitem,
                                        gpointer         user_data){
}


#if (GTK_MAJOR_VERSION >= 4)
#else
/* start of chooser insert */
extern "C" G_MODULE_EXPORT
void
on_coords_filechooserdialog1_response  (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   // I don't think that this is used now - see
   // on_coords_filechooser_dialog_response()

   if (response_id == GTK_RESPONSE_OK) {
      const gchar *filename;
      GtkWidget *coords_fileselection1;
      GtkWidget *combobox;
      int recentre_on_read_pdb_flag = 0;
      short int move_molecule_here_flag = 0;
      int active_index;
      GSList *sel_files;
      /*   GFile  *gfile; for xxx_get_files(), which we don't use (too modern) */

      coords_fileselection1 = widget_from_builder( "coords_filechooserdialog1");

      combobox = widget_from_builder("coords_filechooserdialog1_recentre_combobox");

      if (combobox) {
         active_index = gtk_combo_box_get_active(GTK_COMBO_BOX(combobox));
         /* active items: 0: recentre view on new molecule
            1: don't recentre new molecule
            2: recentre new molecule on view
         */
         if (active_index == 0)
            recentre_on_read_pdb_flag = 1;
         if (active_index == 1)
            recentre_on_read_pdb_flag = 0;
         if (active_index == 2)
            move_molecule_here_flag = 1;
      }

      save_directory_from_filechooser(coords_fileselection1);

      /* no longer used now that we use ..._get_files (multiple files
         (GFiles) returned in a GSList). */
      /*   filename = gtk_file_chooser_get_filename  */
      /*      (GTK_FILE_CHOOSER(coords_fileselection1)); */

      sel_files = gtk_file_chooser_get_files(GTK_FILE_CHOOSER(coords_fileselection1));
      while (sel_files) {
         filename = (char *) sel_files->data;

         // printf("DEBUG:: filename: %s\n", filename);

         /*     From here, we go into c++ (that's why the c++ function
                handle_read_draw needs to be declared external) and read the
                molecule and display it. */

         if (move_molecule_here_flag) {
            handle_read_draw_molecule_and_move_molecule_here(filename);
         } else {
            if (recentre_on_read_pdb_flag)
               handle_read_draw_molecule_with_recentre(filename, 1);
            else
               handle_read_draw_molecule_with_recentre(filename, 0); // no recentre
         }
         sel_files = sel_files->next;
      }

      gtk_widget_hide(coords_fileselection1);

   } else {
      GtkWidget *coords_fileselection1 = widget_from_builder("coords_filechooserdialog1");
      gtk_widget_hide(coords_fileselection1);
   }
}
#endif

extern "C" G_MODULE_EXPORT
void
on_coords_filechooserdialog1_destroy(GtkWidget       *object,
                                     gpointer         user_data) {

   // hides, doesn't destroy

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *coords_fileselection1 = widget_from_builder("coords_filechooserdialog1");
  gtk_widget_hide(coords_fileselection1);
}


extern "C" G_MODULE_EXPORT
void
on_coords_filechooserdialog1_recentre_checkbutton_toggled(GtkToggleButton *togglebutton,
                                                          gpointer         user_data) {

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      set_recentre_on_read_pdb(1);
   } else {
      set_recentre_on_read_pdb(0);
   }
}


extern "C" G_MODULE_EXPORT
void
on_dataset_filechooserdialog1_response (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      gchar *copied_filename;
      int auto_read_flag = 0, ismtz = 0, ismtzauto = 0, iscnsauto = 0;

      GtkWidget *dataset_filechooser = widget_from_builder("dataset_filechooserdialog1");
      save_directory_from_filechooser(dataset_filechooser);
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dataset_filechooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);

      /*    printf("dataset filename: %s\n", filename); */

      copied_filename = (char *) malloc(strlen(file_name) + 1);
      strcpy(copied_filename, file_name);

      auto_read_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dataset_filechooser), "auto_read_flag"));
      ismtz = is_mtz_file_p(file_name);
      if (ismtz) ismtzauto = mtz_file_has_phases_p(file_name);
      else       iscnsauto = cns_file_has_phases_p(file_name);

      if ( ismtzauto || iscnsauto ) {

	 if (auto_read_flag)
	    wrapped_auto_read_make_and_draw_maps(file_name);
	 else
	    /* this does a create_column_label_window, fills and displays it. */
	    manage_column_selector(copied_filename);

      } else {

	 /* no phases path */
	 if (auto_read_flag) printf ("INFO:: This file is not a map coefficient file. Coot can auto-read\nINFO::  - MTZ files from refmac, phenix.refine, phaser, parrot, dm.\nINFO::  - CNS files (new 2009 format only) with cell, symops, F1, F2.\n");
	 if ( ismtz )
	    calc_phases_generic(file_name);
	 else
	    /* try to read as a phs, cif etc... */
	    manage_column_selector(copied_filename);
      }
      gtk_widget_hide(dataset_filechooser);
      free(copied_filename);
   } else {
      GtkWidget *dataset_filechooser = widget_from_builder("dataset_filechooserdialog1");
      gtk_widget_hide(dataset_filechooser);
   }

}

extern "C" G_MODULE_EXPORT
gboolean
on_dataset_filechooser_dialog_delete_event(GtkWidget       *widget,
                                           GdkEvent        *event,
                                           gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}


extern "C" G_MODULE_EXPORT
void
on_dataset_filechooserdialog1_destroy (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *dataset_fileselection1 = widget_from_builder(
                                                "dataset_filechooserdialog1");

  gtk_widget_hide(dataset_fileselection1);
}

// I don't think that this is used now - delete it
extern "C" G_MODULE_EXPORT
void
on_map_name_filechooserdialog1_response(GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   // ------------------------------------- this function is not used. ----------------------------
   // 20220809-PE Delete it.

   std::cout << "------------------------ this function is not used ---------------------------"
             << std::endl;

   if (response_id == GTK_RESPONSE_OK) {

      const char* filename;
      char *sfile;
      GtkWidget* map_name_fileselection1;
      short int is_diff_map_flag = 0;

      GtkWidget *map_name_filechooser = widget_from_builder("map_name_filechooserdialog1");
      save_directory_from_filechooser(map_name_filechooser);

      /* I don't think that we need to malloc this. */

      GtkWidget *checkbutton = widget_from_builder("map_filechooser_is_difference_map_button");
      if (checkbutton)
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
	    is_diff_map_flag = 1;

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(map_name_filechooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);
      // 20220601-PE old/simple
      // filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(map_name_fileselection1));

      std::cout << "INFO :: CCP4 map filename: " << file_name << std::endl;
      sfile = (char *) malloc (1001);
      strncpy(sfile, file_name, 1000);

      gtk_widget_hide(map_name_filechooser); /* the file browser,
						      when destroyed,
						      scribbles over
						      filename. */
      handle_read_ccp4_map_internal(sfile, is_diff_map_flag);
      free(sfile);

   } else {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

}


extern "C" G_MODULE_EXPORT
void
on_map_filechooser_is_difference_map_button_toggled (GtkToggleButton *togglebutton,
                                                     gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_map_name_filechooserdialog1_destroy (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *map_name_fileselection1 = widget_from_builder(
                                                "map_name_filechooserdialog1");

  gtk_widget_hide(map_name_fileselection1);
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
   gtk_widget_hide(phs_filechooser);

}


extern "C" G_MODULE_EXPORT
void
on_phs_coordinates_filechooserdialog1_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *phs_fileselection1 = widget_from_builder("phs_coordinates_filechooserdialog1");
  gtk_widget_hide(phs_fileselection1);
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
    gtk_widget_hide(fileselection);
  } else {
    GtkWidget *fileselection = widget_from_builder("save_coords_filechooserdialog1");
    gtk_widget_hide(fileselection);
  }
}


extern "C" G_MODULE_EXPORT
void
on_save_coords_filechooserdialog1_destroy
					(GtkWidget * object,
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = widget_from_builder(
                                                "save_coords_filechooserdialog1");

  gtk_widget_hide(fileselection);
}


extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_filechooserdialog1_response(GtkDialog * dialog,
                                              gint response_id,
                                              gpointer user_data) {

   /*   int new_compid_idx; */
   int imol_enc = -999997;	/* unset */
   GtkWidget *dictionary_molecule_selector_option_menu;
   GtkWidget *menu;
   GtkWidget *active_menu_item;
   GtkWidget *checkbutton = widget_from_builder("cif_dictionary_file_selector_create_molecule_checkbutton");
   short int new_molecule_checkbutton_state = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
      new_molecule_checkbutton_state = 1;

   if (response_id == GTK_RESPONSE_OK) {

      GtkWidget *cif_dictionary_filechooser = widget_from_builder("cif_dictionary_filechooserdialog1");
      save_directory_from_filechooser(cif_dictionary_filechooser);
      // filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(cif_dictionary_filechooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);
     
      // dictionary_molecule_selector_option_menu =
      // widget_from_builder("cif_dictionary_file_selector_molecule_select_option_menu");

      /* GTK3 FIXME
         if (dictionary_molecule_selector_option_menu) {
         menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(dictionary_molecule_selector_option_menu));
         if (menu) {
         active_menu_item = gtk_menu_get_active(GTK_MENU(menu));
         if (active_menu_item) {
         imol_enc = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(active_menu_item),
         "select_molecule_number"));
         }
         }
         } else {
         printf("-------- missing dictionary_molecule_selector_option_menu ---\n");
         }
      */

      handle_cif_dictionary_for_molecule(file_name, imol_enc, new_molecule_checkbutton_state);

      gtk_widget_hide(cif_dictionary_filechooser);
   } else {
      GtkWidget *filechooser = widget_from_builder("cif_dictionary_filechooserdialog1");
      gtk_widget_hide(filechooser);
   }

}


extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_filechooserdialog1_destroy
					(GtkWidget * object,
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = widget_from_builder("cif_dictionary_filechooserdialog1");
  gtk_widget_hide(fileselection);
}


extern "C" G_MODULE_EXPORT
void
on_run_script_filechooser_dialog_response(GtkDialog * dialog,
                                                              gint response_id,
                                                              gpointer user_data) {

   GtkWidget *file_chooser = widget_from_builder("run_script_filechooser_dialog");
   if (response_id == GTK_RESPONSE_OK) {
      // const char *script_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(file_chooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *script_filename = g_file_info_get_name(file_info);
      run_script(script_filename);
      gtk_widget_hide(file_chooser);
   }
   gtk_widget_hide(file_chooser);

}

extern "C" G_MODULE_EXPORT
void
on_run_script_filechooser_dialog_file_activated(GtkFileChooser* dialog,
                                                gpointer user_data) {

   GtkWidget *file_chooser = widget_from_builder("run_script_filechooser_dialog");
   // const char *script_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
   GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(file_chooser));
   GError *error = NULL;
   GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                            G_FILE_QUERY_INFO_NONE, NULL, &error);
   const char *script_filename = g_file_info_get_name(file_info);
   run_script(script_filename);
   gtk_widget_hide(file_chooser);
}



extern "C" G_MODULE_EXPORT
void
on_run_script_filechooserdialog1_destroy
					(GtkWidget * object,
					gpointer user_data) {

   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
   GtkWidget *file_chooser = widget_from_builder("run_script_filechooser_dialog");
   gtk_widget_hide(file_chooser);
}

#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_run_script_filechooser_dialog_delete_event(GtkWidget       *widget,
                                              GdkEvent        *event,
                                              gpointer         user_data) {

   return gboolean(TRUE);
}
#endif



#if (GTK_MAJOR_VERSION >= 4)
// GtkFileChooserConfirmation is not a thing
#else
GtkFileChooserConfirmation
on_save_symmetry_coords_filechooserdialog1_confirm_overwrite (GtkFileChooser * filechooser,
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
on_save_symmetry_coords_filechooserdialog1_response(GtkDialog * dialog,
                                                    gint response_id,
                                                    gpointer user_data) {
   if (response_id == GTK_RESPONSE_OK) {
      GtkWidget *w = widget_from_builder("save_symmetry_coords_filechooserdialog1");
      save_symmetry_coords_from_filechooser(w);
      gtk_widget_hide(w);
   } else {
      GtkWidget *coords_fileselection1 = widget_from_builder("save_symmetry_coords_filechooserdialog1");
      gtk_widget_hide(coords_fileselection1);
   }
}


extern "C" G_MODULE_EXPORT
void
on_save_symmetry_coords_filechooserdialog1_destroy
					(GtkWidget * object,
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *coords_fileselection1 = widget_from_builder(
                                                "save_symmetry_coords_filechooserdialog1");

  gtk_widget_hide(coords_fileselection1);
}


#if (GTK_MAJOR_VERSION >= 4)
#else
GtkFileChooserConfirmation
on_save_state_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser,
					gpointer user_data)
{

  if (file_chooser_overwrite_state() == 1) {

    return GTK_FILE_CHOOSER_CONFIRMATION_CONFIRM;

  } else {

    return GTK_FILE_CHOOSER_CONFIRMATION_ACCEPT_FILENAME;

  }
}
#endif

extern "C" G_MODULE_EXPORT
void
on_save_state_filechooserdialog1_response (GtkDialog * dialog,
					gint response_id,
					gpointer user_data)
{
   if (response_id == GTK_RESPONSE_OK) {
      GtkWidget *w = widget_from_builder("save_state_filechooserdialog1");

      // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(w));

      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);
   
      save_state_file(file_name);    /* write the file */
      set_save_state_file_name(file_name); /* save as a static in graphics_info_t */
      gtk_widget_hide(w);
   } else {
      GtkWidget *coords_fileselection1 = widget_from_builder("save_state_filechooserdialog1");
      gtk_widget_hide(coords_fileselection1);
   }
}


extern "C" G_MODULE_EXPORT
void
on_save_state_filechooserdialog1_destroy (GtkWidget * object,
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *coords_fileselection1 = widget_from_builder(
                                                "save_state_filechooserdialog1");

  gtk_widget_hide(coords_fileselection1);
}


#if (GTK_MAJOR_VERSION >= 4)
#else
GtkFileChooserConfirmation
on_screendump_filechooser_dialog_confirm_overwrite
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
on_screendump_filechooser_dialog_response (GtkDialog * dialog,
                                           gint response_id,
                                           gpointer user_data) {

   GtkWidget *file_chooser = widget_from_builder("screendump_filechooser_dialog");
   if (response_id == GTK_RESPONSE_OK) {

      int image_type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(file_chooser), "image_type"));
      // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));
      GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(file_chooser));
      GError *error = NULL;
      GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                               G_FILE_QUERY_INFO_NONE, NULL, &error);
      const char *file_name = g_file_info_get_name(file_info);

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
   gtk_widget_hide(file_chooser);
}




extern "C" G_MODULE_EXPORT
void
on_screendump_filechooser_dialog_file_activated(GtkFileChooser* dialog,
                                                gpointer user_data) {

   GtkWidget *file_chooser = widget_from_builder("screendump_filechooser_dialog");
   int image_type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(file_chooser), "image_type"));
   // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(file_chooser));

   GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(file_chooser));
   GError *error = NULL;
   GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                            G_FILE_QUERY_INFO_NONE, NULL, &error);
   const char *file_name = g_file_info_get_name(file_info);
   
   if (image_type == COOT_SCREENDUMP_SIMPLE) {
      screendump_image(file_name);
   }
   if (image_type == COOT_SCREENDUMP_POVRAY) {
      make_image_povray(file_name);
   }
   if (image_type == COOT_SCREENDUMP_RASTER3D) {
      make_image_raster3d(file_name);
   }
   gtk_widget_hide(file_chooser);
}


// not used now?
extern "C" G_MODULE_EXPORT
void
on_screendump_filechooserdialog1_destroy (GtkWidget * object,
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = widget_from_builder("screendump_filechooserdialog1");

  gtk_widget_hide(fileselection);
}
/* end of chooser insert */


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
on_fix_atom_togglebutton_toggled       (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    setup_fixed_atom_pick(1, 0);
  } else {
    setup_fixed_atom_pick(0, 0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_unfix_atom_togglebutton_toggled     (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    setup_fixed_atom_pick(1, 1);
  } else {
    setup_fixed_atom_pick(0, 1);
  }
}


extern "C" G_MODULE_EXPORT
void
on_clear_fixed_atoms_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  clear_fixed_atoms_all();
}


extern "C" G_MODULE_EXPORT
void
on_fixed_atom_close_button_clicked     (GtkButton       *button,
                                                            gpointer         user_data) {
  GtkWidget *dialog = widget_from_builder("fixed_atom_dialog");
  gtk_widget_hide(dialog);
}


#if 0 // 20211202-PE there is no destroy on this now. It has been remade (indeed,
      // it was not even in the glade file!)
extern "C" G_MODULE_EXPORT
void
on_fixed_atom_dialog_destroy           (GtkWidget       *object,
                                        gpointer         user_data)
{
  store_fixed_atom_dialog(0);
}
#endif

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
      gtk_widget_hide(dialog);
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


// is this better done with a dialog response now?
extern "C" G_MODULE_EXPORT
void
on_add_rep_add_rep_button_clicked (GtkButton       *button,
                                                       gpointer         user_data)
{
   // GtkWidget *w = widget_from_builder("add_reps_dialog");
   // add_additional_representation_by_widget(w);
   // gtk_widget_hide(w);
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

   gtk_widget_hide(GTK_WIDGET(dialog));

}

// 
extern "C" G_MODULE_EXPORT
void
on_add_reps_dialog_close (GtkDialog *dialog,
                                              gpointer   user_data) {

   gtk_widget_hide(GTK_WIDGET(dialog));

}

extern "C" G_MODULE_EXPORT
void
on_add_rep_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("add_reps_dialog");
  gtk_widget_hide(w);
}

extern "C" G_MODULE_EXPORT
void
on_additional_representation1_activate (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   // 20210929-PE is this connected now?
   std::cout << "debug:: on_additional_representation1_activate() called" << std::endl;
   GtkWidget *w = wrapped_create_add_additional_representation_gui();
   gtk_widget_show(w);
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
on_all1_activate                       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_all2_activate                       (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

// extern "C" G_MODULE_EXPORT
// void
// on_residue_editor_select_monomer_type_ok_button_clicked (GtkButton       *button,
//                                                          gpointer         user_data) {


//   GtkWidget *dialog = widget_from_builder("residue_editor_select_monomer_type_dialog");
//   GtkWidget *combo_box = widget_from_builder("residue_editor_select_monomer_type_combobox");
//   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo_box));
//   show_restraints_editor(t);
//   gtk_widget_hide(dialog);
// }


extern "C" G_MODULE_EXPORT
void
on_residue_editor_select_monomer_type_cancel_button_clicked (GtkButton       *button,
							gpointer         user_data) {

  GtkWidget *dialog = widget_from_builder("residue_editor_select_monomer_type_dialog");
  gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_restraints1_activate                (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_residue_editor_select_monomer_type_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_restraint_editor_add_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("restraints_editor_dialog");
  restraints_editor_add_restraint_by_widget(w);
}


extern "C" G_MODULE_EXPORT
void
on_restraints_editor_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("restraints_editor_dialog");
  if (w) {
    clear_restraints_editor_by_dialog(w);
    gtk_widget_hide(w);
  }

}

extern "C" G_MODULE_EXPORT
void
on_restraints_editor_save_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("restraints_editor_dialog");
  restraints_editor_save_restraint_by_widget(w);
}

extern "C" G_MODULE_EXPORT
void
on_restraints_editor_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("restraints_editor_dialog");
  apply_restraint_by_widget(w);
}

extern "C" G_MODULE_EXPORT
void
on_restraint_editor_delete_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = widget_from_builder("restraints_editor_dialog");
  restraints_editor_delete_restraint_by_widget(w);
}



#if (GTK_MAJOR_VERSION >= 4)
#else
GtkFileChooserConfirmation
on_save_restraint_chooserdialog_confirm_overwrite
                                        (GtkFileChooser  *filechooser,
					 gpointer         user_data) {

   if (file_chooser_overwrite_state() == 1) {
      return GTK_FILE_CHOOSER_CONFIRMATION_CONFIRM;
   } else {
      return GTK_FILE_CHOOSER_CONFIRMATION_ACCEPT_FILENAME;
   }
}
#endif

extern "C" G_MODULE_EXPORT
void
on_save_restraint_chooserdialog_response(GtkDialog       *dialog,
					 gint             response_id,
					 gpointer         user_data) {
/* Maybe there are responses other than OK and cancel, so don't factor
   out the destroy() */
  GtkWidget *w = widget_from_builder("save_restraint_chooserdialog");
  if (response_id == GTK_RESPONSE_OK) {
    save_monomer_restraints_by_widget(dialog);
    gtk_widget_hide(w);
  }
  if (response_id == GTK_RESPONSE_CANCEL) {
    gtk_widget_hide(w);
  }
}

/* This is not the way. */
extern "C" G_MODULE_EXPORT
void
on_save_restraint_chooserdialog_close  (GtkDialog       *dialog,
                                                            gpointer         user_data) {
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_fix_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_fixed_atom_dialog();
  gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_model_refine_dialog_fast_sss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog;
  dialog = wrapped_create_fast_ss_search_dialog();
  gtk_widget_show(dialog);

}


extern "C" G_MODULE_EXPORT
void
on_other_tools_build_na_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
   // GtkWidget *w = create_build_na_dialog();
   GtkWidget *w = widget_from_builder("build_na_dialog");
   gtk_widget_show(w);

}


extern "C" G_MODULE_EXPORT
void
on_build_na_dialog_cancelbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("build_na_dialog");
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_build_na_dialog_okbutton_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w;
   GtkEntry *entry;
   const char *text;
   float r;
   w = widget_from_builder("build_na_dialog");
   entry = (GTK_ENTRY(widget_from_builder("build_na_dialog_radius_entry")));
   text = gtk_editable_get_text(GTK_EDITABLE(entry));
   r = atof(text);
   find_nucleic_acids_local(r);
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_build_na_dialog_radius_entry_activate(GtkEntry        *entry,
                                         gpointer         user_data) {

   /* BL note:: this is almost a repetetion of
      on_build_na_dialog_okbutton_clicked
      we could/should have a function for this?!
    */
   GtkWidget *w;
   const char *text;
   float r;
   w = widget_from_builder("build_na_dialog");
   text = gtk_editable_get_text(GTK_EDITABLE(entry));
   r = atof(text);
   find_nucleic_acids_local(r);
   gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_coot_references_coot_toolbutton_clicked(GtkButton   *toolbutton,
                                           gpointer         user_data) {

   fill_references_notebook(toolbutton, COOT_REFERENCE_COOT);

}

extern "C" G_MODULE_EXPORT
void
on_coot_references_wincoot_toolbutton_clicked(GtkButton   *toolbutton,
                                              gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_WINCOOT);
}


extern "C" G_MODULE_EXPORT
void
on_coot_references_refmac_toolbutton_clicked(GtkButton   *toolbutton,
                                             gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_REFMAC);

}


extern "C" G_MODULE_EXPORT
void
on_coot_references_ssm_toolbutton_clicked(GtkButton   *toolbutton,
                                          gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_SSM);

}


extern "C" G_MODULE_EXPORT
void
on_coot_references_mmdb_toolbutton_clicked(GtkButton   *toolbutton,
                                           gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_MMDB);
}

extern "C" G_MODULE_EXPORT
void
on_coot_references_clipper_toolbutton_clicked(GtkButton   *toolbutton,
                                              gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_CLIPPER);

}


extern "C" G_MODULE_EXPORT
void
on_coot_references_buccaneer_toolbutton_clicked(GtkButton   *toolbutton,
                                                gpointer         user_data) {

   fill_references_notebook(toolbutton, COOT_REFERENCE_BUCCANEER);
}


extern "C" G_MODULE_EXPORT
void
on_coot_references_molprobity_toolbutton_clicked(GtkButton   *toolbutton,
                                                 gpointer         user_data) {

   fill_references_notebook(toolbutton, COOT_REFERENCE_MOLPROBITY);
}


extern "C" G_MODULE_EXPORT
void
on_coot_references_calpha_toolbutton_clicked(GtkButton   *toolbutton,
                                             gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_CALPHA);
}


extern "C" G_MODULE_EXPORT
void
on_coot_references_xligand_toolbutton_clicked(GtkButton   *toolbutton,
                                              gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_XLIGAND);
}

extern "C" G_MODULE_EXPORT
void
on_coot_references_eds_toolbutton_clicked(GtkButton   *toolbutton,
                                          gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_EDS);
}


extern "C" G_MODULE_EXPORT
void
on_coot_references_others_toolbutton_clicked(GtkButton   *toolbutton,
                                             gpointer         user_data) {
  fill_references_notebook(toolbutton, COOT_REFERENCE_OTHERS);
}


extern "C" G_MODULE_EXPORT
void
on_coot_references_closebutton_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog;
  dialog = widget_from_builder("coot_references_dialog");
  gtk_widget_hide(dialog);

}




extern "C" G_MODULE_EXPORT
void
on_refine_params_use_ramachandran_goodness_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton)))
    set_refine_ramachandran_angles(1);
  else
    set_refine_ramachandran_angles(0);

}


extern "C" G_MODULE_EXPORT
void
on_edit_chi_angles_add_hydrogen_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *vbox = widget_from_builder("edit_chi_angles_vbox");

   if (gtk_toggle_button_get_active(togglebutton)) {
      set_find_hydrogen_torsions(1);
   } else {
      set_find_hydrogen_torsions(0);
   }
   fill_chi_angles_vbox(vbox);
}

#if (GTK_MAJOR_VERSION >= 4)
// is this new or ancient - haha! I don't know
#else
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rot_trans_toolbutton_show_menu
                                        (GtkMenuToolButton *toolbutton,
					 gpointer         user_data) {

   GtkWidget *menu_item = NULL;
   GtkWidget *menu = gtk_menu_tool_button_get_menu(toolbutton);
   GList *children = gtk_container_get_children(GTK_CONTAINER(menu));
   if (get_rot_trans_object_type() == ROT_TRANS_TYPE_ZONE) {
      menu_item = GTK_WIDGET(children->data);
   }
   if (get_rot_trans_object_type() == ROT_TRANS_TYPE_CHAIN) {
      children = children->next;
      menu_item = GTK_WIDGET(children->data);
   }
   if (get_rot_trans_object_type() == ROT_TRANS_TYPE_MOLECULE) {
      children = children->next;
      children = children->next;
      menu_item = GTK_WIDGET(children->data);
   }
   gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
}
#endif

#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
void
on_model_toolbar_rot_trans_toolbutton_clicked (GtkMenuToolButton *toolbutton,
                                               gpointer         user_data) {

  do_rot_trans_setup(1);
}
#endif


extern "C" G_MODULE_EXPORT
void
on_find_ligands_search_all_radiobutton_toggled(GtkButton       *button,
                                               gpointer         user_data) {
}


extern "C" G_MODULE_EXPORT
void
on_find_ligands_search_here_radiobutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  set_ligand_dialog_number_of_sites_sensitivity(GTK_WIDGET(button));

}

extern "C" G_MODULE_EXPORT
void
on_symmetry_controller_dialog_destroy  (GtkWidget       *object,
                                        gpointer         user_data)
{
  set_symmetry_controller_dialog_widget(0);
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
on_map_sharpening1_activate            (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_map_sharpening_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_map_sharpening_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = widget_from_builder("map_sharpening_dialog");
   gtk_widget_hide(w);

}

extern "C" G_MODULE_EXPORT
void
on_map_sharpening_optimize_button_clicked ( GtkButton       *button,
                                                                gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("map_sharpening_dialog");
   calc_and_set_optimal_b_factor(w);
}

extern "C" G_MODULE_EXPORT
void
on_map_sharpening_reset_button_clicked(GtkButton       *button,
                                                           gpointer         user_data)
{
    // reset to zero!?
    GtkWidget *h_scale = widget_from_builder("map_sharpening_hscale");
    GtkAdjustment *adj = gtk_range_get_adjustment(GTK_RANGE(h_scale));
    gtk_adjustment_set_value(adj, 0.);

}

extern "C" G_MODULE_EXPORT
void
on_map_sharpening_dialog_response() {

   // there is only one response
   GtkWidget *dialog = widget_from_builder("map_sharpening_dialog");
   gtk_widget_hide(dialog);

}


extern "C" G_MODULE_EXPORT
void
on_map_sharpening_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_baton_build_params_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("baton_build_params_dialog");
  set_baton_build_params_from_widget(w);
  gtk_widget_hide(w);

}


extern "C" G_MODULE_EXPORT
void
on_baton_build_params_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("baton_build_params_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_baton_build_set_params_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

   // GtkWidget *w = create_baton_build_params_dialog();
   GtkWidget *w = widget_from_builder("baton_build_params_dialog");
   gtk_widget_show(w);

}

extern "C" G_MODULE_EXPORT
void
on_coords_toolbutton_clicked           (GtkButton   *toolbutton,
                                        gpointer         user_data)
{
  open_coords_dialog();
}

extern "C" G_MODULE_EXPORT
void
on_go_to_atom_toolbutton_clicked(GtkButton   *toolbutton,
                                                     gpointer         user_data) {
   GtkWidget *widget = wrapped_create_goto_atom_window(); // uses builder
   gtk_widget_show(widget);
}

extern "C" G_MODULE_EXPORT
void
on_go_to_ligand_button_clicked(GtkButton *button,
                               gpointer   user_data) {
  go_to_ligand();
}

extern "C" G_MODULE_EXPORT
void
on_move_molecule_here_big_molecules_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *dialog = widget_from_builder("move_molecule_here_dialog");
  fill_move_molecule_here_dialog(dialog);
}

extern "C" G_MODULE_EXPORT
void
on_python_scripting_button(GtkToggleButton *togglebutton, gpointer user_data) {
   reveal_python_scripting_entry();
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

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_environment_distance_max_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data) {

   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
       execute_environment_settings(widget);
  }
  return FALSE;
}
#endif

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_environment_distance_min_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data) {

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
       execute_environment_settings(widget);
  }
  return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_pisa_interfces_close_button_clicked (GtkButton       *button,
                                                            gpointer         user_data) {

  GtkWidget *w = widget_from_builder("pisa_interfaces_dialog");
  gtk_widget_hide(w);

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
on_map_opacity_hscale_value_changed (GtkRange        *range,
                                                         gpointer         user_data) {

  GtkAdjustment *adjustment;
  float fvalue;
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(range), "imol"));

  adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
  fvalue = 0.01 * gtk_adjustment_get_value(adjustment);
  if (fvalue > 0.99)
    fvalue = 1.0;

  set_solid_density_surface_opacity(imol, fvalue);

}


extern "C" G_MODULE_EXPORT
void
on_refine_params_weight_matrix_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{

  GtkWidget *entry = widget_from_builder("refine_params_weight_matrix_entry");
  struct entry_info_t ei = coot_entry_to_val(GTK_ENTRY(entry));
  if (ei.float_is_set)
    set_matrix(ei.val_as_float);

}


// extern "C" G_MODULE_EXPORT
// void
// on_remarks_browser1_activate (GMenuItem     *menuitem,
//                                                   gpointer         user_data) {

//   GtkWidget *w = wrapped_create_remarks_browser_molecule_chooser_dialog();
//   std::cout << "in on_remarks_browser1_activate got w " << w << std::endl;
//   gtk_widget_show(w);

// }

// extern "C" G_MODULE_EXPORT
// void
// on_remarks_browser_molecule_chooser_cancel_button_clicked
//                                         (GtkButton       *button,
// 					 gpointer         user_data) {
//    // GtkWidget *w = widget_from_builder("remarks_browser_molecule_chooser_dialog");
//    // gtk_widget_hide(w);

//    GtkWidget *dialog = widget_from_builder("remarks_browser_molecule_chooser_dialog");
//    gtk_widget_hide(dialog);

// }


extern "C" G_MODULE_EXPORT
void
on_remarks_browser1_activate (GMenuItem     *menuitem,
                              gpointer         user_data) {

   GtkWidget *w = wrapped_create_remarks_browser_molecule_chooser_dialog(); // uses builder
   gtk_widget_show(w);

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
   gtk_widget_hide(GTK_WIDGET(dialog));
}

extern "C" G_MODULE_EXPORT
void
on_remarks_browser_molecule_chooser_ok_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

  GtkWidget *w = widget_from_builder("remarks_browser_molecule_chooser_dialog");
  gtk_widget_hide(w);
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
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_fix_nomenclature_errors_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("fix_nomenclature_errors_dialog");
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_ligand_builder1_activate(GMenuItem *menuitem,
                            gpointer   user_data) {
#if (GTK_MAJOR_VERSION >= 4)
   std::cout << "FIXME:: start_ligand_builder_gui_internal() " << std::endl;
#else
   start_ligand_builder_gui_internal(menuitem, user_data);
#endif
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
  gtk_widget_hide(w);
}


extern "C" G_MODULE_EXPORT
void
on_multi_residue_torsion_OK_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = widget_from_builder("multi_residue_torsion_dialog");
  gtk_widget_hide(w);
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
  gtk_widget_hide(w);
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
  gtk_widget_hide(w);
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

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_keyboard_go_to_residue_entry_key_release_event(GtkWidget       *widget,
                                                  GdkEventKey     *event,
                                                  gpointer         user_data) {
   return FALSE; // use the key press event.
}
#endif

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_keyboard_go_to_residue_entry_key_press_event (GtkWidget       *widget,
                                                 GdkEventKey     *event,
                                                 gpointer         user_data) {

   GtkWidget *w = widget_from_builder("keyboard_go_to_residue_window");
   const gchar *text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(widget)));
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      handle_go_to_residue_keyboarding_mode(text);
      gtk_widget_hide(w);
      return TRUE;
   }
   if (event->keyval == GDK_KEY_Escape) {
      gtk_widget_hide(w);
      return TRUE;
   }
   return FALSE;
}
#endif

#if (GTK_MAJOR_VERSION >= 4)
#else
extern "C" G_MODULE_EXPORT
gboolean
on_keyboard_go_to_residue_window_delete_event(GtkWidget       *widget,
                                              GdkEvent        *event,
                                              gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_mogul_geometry_dialog_close_button_clicked (GtkButton       *button,
                                               gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("mogul_geometry_results_table_dialog");
   /* And the histogram?  How do I look that up? */
   gtk_widget_hide(dialog);
}


extern "C" G_MODULE_EXPORT
void
on_ligand_check_okbutton_clicked(GtkButton       *button,
                                 gpointer         user_data) {

  GtkWidget *w = widget_from_builder("ligand_check_dialog");
  gtk_widget_hide(w);

}

extern "C" G_MODULE_EXPORT
void
on_generic_objects_dialog_closebutton_clicked(GtkButton       *button,
                                              gpointer         user_data) {

   std::cout << "------------------------ no longer a thing: "
             << "on_generic_objects_dialog_closebutton_clicked" << std::endl;

#if 0 // 20211007-PE
   GtkWidget *w = widget_from_builder("generic_objects_dialog");
   gtk_widget_hide(w);
   clear_generic_objects_dialog_pointer();
   graphics_draw();
#endif

}

/* I don't know how this function gets activated, it's not the close button of the dialog */
extern "C" G_MODULE_EXPORT
void
on_generic_objects_dialog_close (GtkDialog       *dialog,
                                                     gpointer         user_data) {

  clear_generic_objects_dialog_pointer(); /* needed here? */
  graphics_draw();

}

extern "C" G_MODULE_EXPORT
void
on_generic_objects_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data) {

  clear_generic_objects_dialog_pointer();

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
on_file_export_map1_activate                (GMenuItem     *menuitem,
                                                                 gpointer         user_data) {

   short int is_fragment = false;
   export_map_gui(is_fragment);

}

extern "C" G_MODULE_EXPORT
void
on_file_export_map_fragment1_activate       (GMenuItem     *menuitem,
                                                            gpointer         user_data) {

   short int is_fragment = true;
   export_map_gui(is_fragment);
}

extern "C" G_MODULE_EXPORT
void
on_export_map_dialog_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data) {

  on_export_map_dialog_ok_button_clicked_cc(button);
}

extern "C" G_MODULE_EXPORT
void
on_export_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

  GtkWidget *w = widget_from_builder("export_map_dialog");
  gtk_widget_hide(w);
}

/* void */
/* on_export_map_filechooserdialog_cancel_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data) */
/* { */

/*   GtkWidget *w = widget_from_builder("export_map_filechooserdialog"); */
/*   gtk_widget_hide(w); */


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
/*   gtk_widget_hide(w); */

/* } */


extern "C" G_MODULE_EXPORT
void
on_export_map_dialog_response (GtkDialog       *dialog,
                                                   gint             response_id,
                                                   gpointer         user_data) {


   if (response_id == GTK_RESPONSE_OK) {
      GtkWidget *file_chooser_dialog = widget_from_builder("export_map_file_chooser_dialog");
      GtkWidget *combobox            = widget_from_builder("export_map_map_combobox");
      GtkWidget *radius_entry        = widget_from_builder("export_map_radius_entry");
      int imol_map = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
      int is_map_fragment = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "is_map_fragment"));
      // std::cout << "debug:: in on_export_map_dialog_response() imol_map is " << imol_map << std::endl;
      const char *entry_text = gtk_editable_get_text(GTK_EDITABLE(radius_entry));
      // std::cout << "debug:: in on_export_map_dialog_response() got entry_text \"" << entry_text << "\"" << std::endl;
      GString* text_copy   = g_string_new(entry_text);
      GString* text_copy_2 = g_string_new(entry_text);
      // gtk_widget_hide(GTK_WIDGET(dialog));
      gtk_widget_show(file_chooser_dialog);
      g_object_set_data(G_OBJECT(file_chooser_dialog), "map_molecule_number", GINT_TO_POINTER(imol_map));
      g_object_set_data(G_OBJECT(file_chooser_dialog), "is_map_fragment",     GINT_TO_POINTER(is_map_fragment));
      // std::cout << "debug:: in on_export_map_dialog_response() storing entry text " << text_copy << std::endl;
      g_object_set_data(G_OBJECT(file_chooser_dialog), "export_map_radius_entry_text", text_copy);
      // decoded:
      // int imol_map = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "map_molecule_number"));

   }
   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
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

      if (is_map_fragment > 0) {
         if (GTK_IS_FILE_CHOOSER(dialog)) {
            // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

            GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
            GError *error = NULL;
            GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                                     G_FILE_QUERY_INFO_NONE, NULL, &error);
            const char *file_name = g_file_info_get_name(file_info);

            GString *txt_radius_str = static_cast<GString *>(g_object_get_data(G_OBJECT(dialog), "export_map_radius_entry_text"));
            const char *entry_text = g_string_free(txt_radius_str, FALSE); // leaking entry_text - ho hum.
            if (entry_text == 0) {
               std::cout << "ERROR:: entry_text is null " << std::endl;
            }
            export_map_fragment_with_text_radius(imol_map, entry_text, file_name);
         }
      } else {
         if (GTK_IS_FILE_CHOOSER(dialog)) {
            // const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
            GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
            GError *error = NULL;
            GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                                     G_FILE_QUERY_INFO_NONE, NULL, &error);
            const char *file_name = g_file_info_get_name(file_info);
            export_map(imol_map, file_name);
         }
      }
   }
   
   gtk_widget_hide(GTK_WIDGET(dialog));
}

#include "cfc-widgets-c-interface.h"

extern "C" G_MODULE_EXPORT
void
on_cfc_dialog_response                 (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   printf("handle response id %d\n", response_id);

   if (response_id == GTK_RESPONSE_CLOSE) {
/*       close_cfc_dialog(GTK_WIDGET(dialog)); */
      gtk_widget_hide(GTK_WIDGET(dialog)); /* FIXME - not a widget for gtk_widget_unref no destroy */

   }

}


extern "C" G_MODULE_EXPORT
void
on_dynarama_outliers_only_togglebutton_toggled (GtkToggleButton *togglebutton,
						gpointer         user_data)
{
   GtkWidget *window = widget_from_builder("dynarama_window");
   toggle_dynarama_outliers(window, gtk_toggle_button_get_active(togglebutton)); /* get the imol from window */
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
on_edit_copy_molecule_activate        (GMenuItem     *menuitem,
                                                           gpointer         user_data)
{
  do_edit_copy_molecule();
}


extern "C" G_MODULE_EXPORT
void
on_edit_copy_fragment_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_edit_copy_fragment();
}

extern "C" G_MODULE_EXPORT
void
on_copy_fragment_dialog_response(GtkDialog *dialog,
                                                     gint response_id,
                                                     gpointer user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      graphics_info_t g;
      GtkWidget *entry = widget_from_builder("copy_fragment_atom_selection_entry");
      std::string text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(entry)));
      GtkWidget *combobox = GTK_WIDGET(g_object_get_data(G_OBJECT(dialog), "combobox"));
      int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
      int imol_new = new_molecule_by_atom_selection(imol, text.c_str());
      GtkWidget *checkbutton = widget_from_builder("copy_fragment_move_molecule_here_checkbutton");
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
         move_molecule_to_screen_centre_internal(imol_new);
      if (is_valid_model_molecule(imol_new))
         gtk_widget_hide(GTK_WIDGET(dialog));
   }
   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

}


extern "C" G_MODULE_EXPORT
void
on_edit_replace_residue_activate      (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_edit_replace_residue();
}


extern "C" G_MODULE_EXPORT
void
on_edit_replace_fragment_activate     (GMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_edit_replace_fragment();
}


extern "C" G_MODULE_EXPORT
void
on_edit_renumber_residues_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_renumber_residue_range_dialog();
  gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_edit_change_chain_ids1_activate (GMenuItem     *menuitem,
                                                        gpointer         user_data) {

   GtkWidget *w = wrapped_create_change_chain_id_dialog(); // uses builder
   gtk_widget_show(w);
}

extern "C" G_MODULE_EXPORT
void
on_residue_editor_select_monomer_type_dialog_close (GtkDialog       *dialog,
                                                    gpointer         user_data) {
   // need to add a "Response ID" to the button in glade
   gtk_widget_hide(GTK_WIDGET(dialog));
}

extern "C" G_MODULE_EXPORT
void
on_residue_editor_select_monomer_type_dialog_response (GtkDialog       *dialog,
                                                                           gint response_id,
                                                                           gpointer         user_data) {

   // need to add a "Response ID" to the button in glade

   std::cout << "on_residue_editor_select_monomer_type_dialog_response()" << std::endl;
  if (response_id == GTK_RESPONSE_OK) {
     GtkWidget *combo_box = widget_from_builder("residue_editor_select_monomer_type_combobox");
     std::string rt = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combo_box));
     show_restraints_editor(rt);
  }

  if (response_id == GTK_RESPONSE_CANCEL) {
  }
  gtk_widget_hide(GTK_WIDGET(dialog));
  
  
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

   gtk_widget_hide(GTK_WIDGET(dialog));
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
on_edit_restraints_activate(GMenuItem     *menuitem,
                                                gpointer         user_data)
{
   std::cout << "on_edit_restraints_activate() " << std::endl;
   GtkWidget *w =  wrapped_create_residue_editor_select_monomer_type_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_edit_merge_molecules1_activate(GMenuItem *menuitem,
                                  gpointer   user_data) {

   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   set_transient_for_main_window(w);
   gtk_widget_show(w);
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
on_mutate_molecule_sequence_text_insert_at_cursor
                                        (GtkTextView     *textview,
                                        gchar           *string,
                                        gpointer         user_data)
{
   GtkWidget *res_no_1_widget = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *res_no_2_widget = widget_from_builder("mutate_molecule_resno_2_entry");
   GtkWidget *text_widget     = widget_from_builder("mutate_molecule_sequence_text");
   GtkWidget *label_widget    = widget_from_builder("mutate_residue_range_counts_label");
   mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_mutate_molecule_sequence_text_key_release_event(GtkWidget       *widget,
                                                   GdkEventKey     *event,
                                                   gpointer         user_data) {

   GtkWidget *res_no_1_widget = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *res_no_2_widget = widget_from_builder("mutate_molecule_resno_2_entry");
   GtkWidget *text_widget     = widget_from_builder("mutate_molecule_sequence_text");
   GtkWidget *label_widget    = widget_from_builder("mutate_residue_range_counts_label");
   mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
   return FALSE;
}
#endif


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_mutate_molecule_sequence_text_button_release_event
                                        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data) {

  GtkWidget *res_no_1_widget = widget_from_builder("mutate_molecule_resno_1_entry");
  GtkWidget *res_no_2_widget = widget_from_builder("mutate_molecule_resno_2_entry");
  GtkWidget *text_widget     = widget_from_builder("mutate_molecule_sequence_text");
  GtkWidget *label_widget    = widget_from_builder("mutate_residue_range_counts_label");
  mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
  return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_display_control_last_model_only_button_clicked (GtkButton       *button,
                                                   gpointer         user_data) {
  set_only_last_model_molecule_displayed();

}


extern "C" G_MODULE_EXPORT
void
on_display_control_align_labels_checkbutton_toggled (GtkToggleButton *togglebutton,
                                                     gpointer         user_data) {

  align_labels_checkbutton_toggled(togglebutton);

}



extern "C" G_MODULE_EXPORT
void
on_curlew_install_button_clicked(GtkButton *button,
                                 gpointer   user_data) {

  GtkWidget *dialog = widget_from_builder("curlew_dialog");
  if (dialog) {
     int n_items = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "n_extensions"));
    curlew_dialog_install_extensions(dialog, n_items); /* some of which were selected */
  }
}



extern "C" G_MODULE_EXPORT
void
on_curlew_dialog_close(GtkDialog       *dialog,
                       gpointer         user_data) {
  gtk_widget_hide(GTK_WIDGET(dialog)); /* or maybe hide */
}


extern "C" G_MODULE_EXPORT
void
on_curlew_dialog_response (GtkDialog       *dialog,
                           gint             response_id,
                           gpointer         user_data) {

  /*
  printf("in on_curlew_dialog_response with response_id %d\n", response_id);
  printf("   cf response_id %d\n", GTK_RESPONSE_CLOSE);
  printf("   cf response_id %d\n", GTK_RESPONSE_OK);
  printf("   cf response_id %d\n", GTK_RESPONSE_CANCEL);
  */

  if (response_id == GTK_RESPONSE_CLOSE)
    gtk_widget_hide(GTK_WIDGET(dialog));

}


extern "C" G_MODULE_EXPORT
void
on_modelling_activate(GMenuItem     *menuitem,
                      gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_edit_settings_activate              (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

extern "C" G_MODULE_EXPORT
void
on_calculate_all_molecule_activate (GMenuItem     *menuitem,
                                    gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_calculate_dock_sequence_activate    (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_calculate_map_tools_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_calculate_modules_activate          (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_calculate_ncs_tools_activate        (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


extern "C" G_MODULE_EXPORT
void
on_calculate_pisa_activate             (GMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

#include "cc-interface-scripting.hh"


extern "C" G_MODULE_EXPORT
void
on_calculate_pisa_assemblies_activate(GMenuItem     *menuitem,
                                      gpointer         user_data) {

   safe_python_command("import parse_pisa_xml");
   safe_python_command("parse_pisa_xml.pisa_molecule_chooser_gui('assemblies')");
}

extern "C" G_MODULE_EXPORT
void
on_calculate_pisa_interfaces_activate (GMenuItem     *menuitem,
                                       gpointer         user_data) {

   safe_python_command("import parse_pisa_xml");
   safe_python_command("parse_pisa_xml.pisa_molecule_chooser_gui('interfaces')");
}


extern "C" G_MODULE_EXPORT
void
on_draw_representation_tools_activate(GMenuItem     *menuitem,
                                      gpointer         user_data) {

}


extern "C" G_MODULE_EXPORT
void
on_draw_rock_view_activate(GMenuItem     *menuitem,
                           gpointer         user_data) {
   toggle_idle_rock_function();
}



extern "C" G_MODULE_EXPORT
void
 on_draw_central_atom_label_activate (GMenuItem     *menuitem,
                                      gpointer         user_data) {
   std::cout << "draw central atom " << std::endl;
}

extern "C" G_MODULE_EXPORT
void
on_draw_label_neighbours_activate  (GMenuItem     *menuitem,
                                    gpointer         user_data)
{
   label_neighbours();
}

extern "C" G_MODULE_EXPORT
void
on_draw_label_atoms_in_residue_activate  (GMenuItem     *menuitem,
                                          gpointer         user_data) {
   label_atoms_in_residue();
}

extern "C" G_MODULE_EXPORT
void
on_draw_label_CA_atoms_activate(GMenuItem     *menuitem,
                                gpointer         user_data) {

   // make this function built-in now - this is the first time that activer_atom() is used in this file.

   auto label_all_CAs = [] (int imol) {
                           graphics_info_t::molecules[imol].add_labels_for_all_CAs();
                   };

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      label_all_CAs(imol);
      graphics_draw();
   }
   
}

extern "C" G_MODULE_EXPORT
void
on_draw_additional_representations_activate  (GMenuItem     *menuitem,
                                                 gpointer         user_data)
{
   std::cout << "on draw additional representations menu item activate" << std::endl;
  GtkWidget *w = wrapped_create_add_additional_representation_gui();
  gtk_widget_show(w);
   
}

extern "C" G_MODULE_EXPORT
void
on_calculate_align_and_mutate_activate (GMenuItem     *menuitem,
                                                           gpointer         user_data) {
   GtkWidget *w = wrapped_create_align_and_mutate_dialog();
   gtk_widget_show(w);
}


extern "C" G_MODULE_EXPORT
void
on_calculate_find_ligands_item_activate() {
   do_find_ligands_dialog();
}

extern "C" G_MODULE_EXPORT
void
on_calculate_find_waters_item_activate() {

   std::cout << "################################ on_calculate_find_waters_item_activate() " << std::endl;
   GtkWidget *dialog = widget_from_builder("find_waters_dialog");
   std::cout << "find_waters_dialog " << dialog << std::endl;
   fill_find_waters_dialog(dialog);
   gtk_widget_show(dialog);
      
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
on_refine_params_geman_mcclure_alpha_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));

   printf("GTK3 FIXME on_refine_params_geman_mcclure_alpha_combobox_changed\n");
   // int active_item_idx = gtk_combo_box_text_get_active(combobox);
   // set_refinement_geman_mcclure_alpha_from_text(active_item_idx, t);
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_lennard_jones_epsilon_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox); // save it for set active item next time
   set_refinement_lennard_jones_epsilon_from_text(active_item_idx, t);
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_rama_restraints_weight_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
   const char *t = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_ramachandran_restraints_weight_from_text(active_item_idx, t);
}


extern "C" G_MODULE_EXPORT
void
on_refine_params_torsions_weight_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
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
on_symmetry_always_on_checkbutton_toggled (GtkToggleButton *togglebutton,
					   gpointer         user_data) {

   GtkWidget *symmetry_on_radio_button = NULL;
   if (gtk_toggle_button_get_active(togglebutton)) {
      add_symmetry_on_to_preferences_and_apply();
      symmetry_on_radio_button = widget_from_builder("show_symmetry_yes_radiobutton");
      if (! gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(symmetry_on_radio_button)))
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(symmetry_on_radio_button), TRUE);
   }

}


extern "C" G_MODULE_EXPORT
void
on_curlew1_activate              (GMenuItem     *menuitem,
                                                      gpointer         user_data) {
  curlew();
}

extern "C" G_MODULE_EXPORT
void
on_show_symmetry_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data) {

  if (gtk_toggle_button_get_active(togglebutton)) {
      set_show_symmetry_master(0);
  }
}


extern "C" G_MODULE_EXPORT
void
on_show_symmetry_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
      set_show_symmetry_master(1);
  }
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_symmetry_radius_entry_key_release_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  const char *text;
  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
     text = gtk_editable_get_text(GTK_EDITABLE(GTK_ENTRY(widget)));
     set_symmetry_size_from_widget(text);
  }
  return TRUE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_hscale_symmetry_colour_value_changed
                                        (GtkRange        *range,
                                        gpointer         user_data)
{
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
on_map_radius_em_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkEntry *entry = GTK_ENTRY(widget_from_builder("map_radius_em_entry"));
  const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  printf("set_density_size_em_from_widget() %s\n", text);
  set_density_size_em_from_widget(text);
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_map_radius_em_entry_key_press_event (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  GtkEntry *entry = (GTK_ENTRY(widget_from_builder("map_radius_em_entry")));
  const char *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter)
     set_density_size_em_from_widget(text);
  return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_simple_refmac_dialog_response       (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data) {

   if (response_id == GTK_RESPONSE_CLOSE) {
      /* do I need to do this? */
      /* gtk_widget_hide(dialog); */
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_OK) {
      simple_refmac_run_refmac(GTK_WIDGET(dialog));
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

}


extern "C" G_MODULE_EXPORT
void
on_simple_refmac_dialog_close (GtkDialog       *dialog,
                                                   gpointer         user_data)
{
   /* Do I need to do anything here? */
}

extern "C" G_MODULE_EXPORT
gboolean
on_simple_refmac_dialog_delete_event(GtkWidget       *widget,
                                                         GdkEvent        *event,
                                                         gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}



extern "C" G_MODULE_EXPORT
void
on_simple_refmac_mtz_file_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = widget_from_builder("simple_refmac_filechooser_dialog");
   GtkWidget *simple_refmac_dialog = widget_from_builder("simple_refmac_dialog");
   /* automtically file filter only mtz files */
   GtkFileFilter *filterselect = gtk_file_filter_new();
   gtk_file_filter_add_pattern(filterselect, "*.mtz");
   gtk_file_chooser_set_filter(GTK_FILE_CHOOSER(w), filterselect);
   /* we need to find the file combo box in the simple refmac dialog */
   g_object_set_data(G_OBJECT(w), "simple_refmac_dialog", simple_refmac_dialog);
   gtk_widget_show(w);
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

   gtk_widget_hide(GTK_WIDGET(dialog));

}

extern "C" G_MODULE_EXPORT
void
on_label_neighbours1_activate(GMenuItem     *menuitem,
                                                  gpointer         user_data)
{
  label_neighbours();
}

extern "C" G_MODULE_EXPORT
void
on_label_atoms_in_residue1_activate    (GMenuItem     *menuitem,
                                                            gpointer         user_data) {

  label_atoms_in_residue();
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_residue_type_chooser_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{
   const char *entry_text = gtk_editable_get_text(GTK_EDITABLE(widget));
   GtkWidget *stub_button = widget_from_builder("residue_type_chooser_stub_checkbutton");
   GtkWidget *window = widget_from_builder("residue_type_chooser_window");
   short int istate = 0;

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
      istate = 1;
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      handle_residue_type_chooser_entry_chose_type(entry_text, istate);
      gtk_widget_hide(window);
   }
   return FALSE;
}
#endif


void handle_map_properties_specularity_change(int imol, GtkWidget *togglebutton) {

   molecule_class_info_t &m = graphics_info_t::molecules[imol];

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
      std::cout << "Turn on specularity " << std::endl;
      GtkWidget *strength_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton),  "strength_entry"));
      GtkWidget *shininess_entry = GTK_WIDGET(g_object_get_data(G_OBJECT(togglebutton), "shininess_entry"));
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
      std::cout << "Turn off specularity " << std::endl;
      m.material_for_maps.turn_specularity_on(false);
   }
   graphics_draw();
}


#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_map_properties_dialog_specularity_strength_entry_key_press_event (GtkWidget       *widget,
                                                                     GdkEventKey     *event,
                                                                     gpointer         user_data)
{

   std::cout << "strength entry key press callback" << std::endl;
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      // the g_object_set_data() for these is done in fill_single_map_properties_dialog_gtk3()
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "imol"));
      GtkWidget *togglebutton = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "specularity_checkbutton"));
      std::cout << "call handle_map_properties_specularity_change() " << std::endl;
      handle_map_properties_specularity_change(imol, togglebutton);
   }
   return FALSE; // otherwise the text can't edited!
}
#endif

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_map_properties_dialog_specularity_shininess_entry_key_press_event (GtkWidget       *widget,
                                                                      GdkEventKey     *event,
                                                                      gpointer         user_data) {
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "imol"));
      GtkWidget *togglebutton = GTK_WIDGET(g_object_get_data(G_OBJECT(widget), "specularity_checkbutton"));
      handle_map_properties_specularity_change(imol, togglebutton);
   }
   return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
void
on_map_properties_dialog_specularity_state_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(togglebutton), "imol"));
   handle_map_properties_specularity_change(imol, GTK_WIDGET(togglebutton));
}

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_single_map_properties_step_size_entry_key_press_event (GtkWidget       *widget,
                                                                              GdkEventKey     *event,
                                                                              gpointer         user_data)
{
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(widget), "imol"));
      std::string t = gtk_editable_get_text(GTK_EDITABLE(widget));
      try {
         float f = coot::util::string_to_float(t);
      }
      catch (const std::runtime_error &rte) {
         std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }
   return FALSE;
}
#endif


extern "C" G_MODULE_EXPORT
gboolean
on_cif_dictionary_filechooser_dialog_delete_event(GtkWidget       *widget,
                                                  GdkEvent        *event,
                                                  gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}



extern "C" G_MODULE_EXPORT
void
on_cif_dictionary_filechooser_dialog_file_activated(GtkFileChooser* dialog,
                                                    gpointer user_data) {

   // 20220319-PE shouldn't need to connect to this says the documentation - hmmm.....
   // const char *fn = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

   GFile *file = gtk_file_chooser_get_file(GTK_FILE_CHOOSER(dialog));
   GError *error = NULL;
   GFileInfo *file_info = g_file_query_info(file, G_FILE_ATTRIBUTE_STANDARD_CONTENT_TYPE,
                                            G_FILE_QUERY_INFO_NONE, NULL, &error);
   const char *file_name = g_file_info_get_name(file_info);
   read_cif_dictionary(file_name);
   gtk_widget_hide(GTK_WIDGET(dialog));
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
      gtk_widget_hide(GTK_WIDGET(dialog));
   }

   if (response_id == GTK_RESPONSE_CANCEL) {
      gtk_widget_hide(GTK_WIDGET(dialog));
   }
}


extern "C" G_MODULE_EXPORT
gboolean
on_keyboard_mutate_dialog_delete_event(GtkWidget       *widget,
                                                           GdkEvent        *event,
                                                           gpointer         user_data) {

   gtk_widget_hide(widget);
   return TRUE;
}

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_keyboard_mutate_entry_key_press_event (GtkWidget       *widget,
                                          GdkEventKey     *event,
                                          gpointer         user_data) {

   GtkWidget *dialog = widget_from_builder("keyboard_mutate_dialog");
   GtkWidget *entry  = widget_from_builder("keyboard_mutate_entry");
   if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
      const char *txt = gtk_editable_get_text(GTK_EDITABLE(entry));
      if (txt) {
         std::string ss(txt);
         gtk_editable_set_text(GTK_EDITABLE(entry), ""); // for next time!
         mutate_active_residue_to_single_letter_code(ss);
      }
      gtk_widget_hide(dialog);
   }
   return FALSE;
}
#endif

#ifdef FIX_THE_KEY_PRESS_EVENTS
extern "C" G_MODULE_EXPORT
gboolean
on_keyboard_mutate_entry_key_release_event(GtkWidget       *widget,
                                           GdkEventKey     *event,
                                           gpointer         user_data) {
   return FALSE;
}
#endif
