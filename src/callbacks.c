/* src/callbacks.c
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


#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 

#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif


#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h> // for keyboarding.

#include "coot-fileselections.h"
#include "coot-preferences.h"
#include "coot-references.h"
#include "rotate-translate-modes.hh"

#include "callbacks.h"
#include "interface.h"
#include "gtk-manual.h"
#include "restraints-editor-c.h"

#include "gtk-widget-conversion-utils.h"


#include "read-phs.h"

/* map colour selection stuff. */


#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "c-interface-refmac.h"
#include "c-interface-refine.h"	/* do things with intermediate atoms */
#include "c-interface-preferences.h"
#include "positioned-widgets.h"

#include "generic-display-objects-c.h"
// #include "curlew.h"

/* This is our data identification string to store
 * data in list items
 */
const gchar *list_item_data_key="list_item_data";

/* This is our data identification string to store
 * data in list items
 */
const gchar *list_item_data_key_for_atoms="list_item_data_for_atoms";




void
on_window1_destroy                     (GtkWidget       *object,
                                        gpointer         user_data)
{
  gtk_main_quit();
}

/* When the user uses the window manager to close coot, this gets called. */
/* When the window manager "close window" events happens it send the
   application a delete_event event.  If we return FALSE from this
   attached function, then a "destroy" signal will be emitted.
   Returning TRUE means that we don't want the window to be destroyed.
   See helloworld.c in the gtk_tut.txt */
gboolean
on_window1_delete_event                (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
/*   printf("---------------------------- on_window1_delete_event() ------------------\n"); */
  /* coot_checked_exit() calls coot_real_exit() and that calls exit(),
     so we don't (normally?) return from this function. */
  coot_no_state_real_exit(0);
  return 0;
}


void
on_open_coordinates1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{ 
  open_coords_dialog();
}


void
on_open_dataset1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  int is_auto_read_fileselection = 0;
  int is;
  GtkWidget *file_filter_button;
  GtkWidget *sort_button;
  GtkWidget *dataset_fileselection1 = coot_dataset_chooser();
  add_ccp4i_project_optionmenu(dataset_fileselection1, COOT_DATASET_FILE_SELECTION);
  file_filter_button = add_filename_filter_button(dataset_fileselection1, 
						  COOT_DATASET_FILE_SELECTION);
  sort_button = add_sort_button_fileselection(dataset_fileselection1); 
  /* add_filename_filter_button needs the sort button to be in place (for now) */
  set_directory_for_fileselection(dataset_fileselection1);

  /* stuff in user data saying if this is autoread or not... */
  is = is_auto_read_fileselection;
  g_object_set_data(G_OBJECT(dataset_fileselection1), "imol", GINT_TO_POINTER(is));
  set_file_selection_dialog_size(dataset_fileselection1);

  set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_fileselection1);
  gtk_widget_show (dataset_fileselection1);
  /* in gtk2 we have to push the buttons after we show the selection */
  push_the_buttons_on_fileselection(file_filter_button, sort_button, 
				    dataset_fileselection1);

}

void
on_auto_open_mtz_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  int is_auto_read_fileselection = 1;
  int is;
  GtkWidget *file_filter_button;
  GtkWidget *sort_button;
  GtkWidget *dataset_fileselection1 = coot_dataset_chooser ();
  add_ccp4i_project_optionmenu(dataset_fileselection1, COOT_DATASET_FILE_SELECTION);
  file_filter_button = add_filename_filter_button(dataset_fileselection1, 
						  COOT_DATASET_FILE_SELECTION);
  sort_button = add_sort_button_fileselection(dataset_fileselection1); 
  set_directory_for_fileselection(dataset_fileselection1);

  /* stuff in user data saying if this is autoread or not... */
  is = is_auto_read_fileselection;
  g_object_set_data(G_OBJECT(dataset_fileselection1), "imol", GINT_TO_POINTER(is));
  set_file_selection_dialog_size(dataset_fileselection1);

  set_transient_and_position(COOT_UNDEFINED_WINDOW, dataset_fileselection1);
  gtk_widget_show (dataset_fileselection1);
  push_the_buttons_on_fileselection(file_filter_button, sort_button, 
				    dataset_fileselection1);
}


void on_file1_activate (GtkMenuItem     *menuitem,
			gpointer         user_data)
{

}


void
on_exit1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
/*   printf("calling gtk_main_quit()\n"); */
/*   gtk_main_quit(); */

/*    printf("---------------------------- on_exit1_activate() ------------------\n"); */
   coot_checked_exit(0); /* without error */
}


void
on_clipping1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_clipping1_activate(); 
}




void
on_cancel_coords_button1_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(button),
						    "coords_fileselection1");

   gtk_widget_destroy(coords_fileselection1);

}


void
on_cancel_dataset_button1_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dataset_fileselection1 = lookup_widget(GTK_WIDGET(button),
						   "dataset_fileselection1");
   gtk_widget_destroy(dataset_fileselection1);
}


void
on_column_label_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *column_label_window = lookup_widget(GTK_WIDGET(button), "column_label_window");

  handle_column_label_make_fourier(column_label_window);

}


void
on_column_label_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   printf("column_label_window CANCEL button clicked\n"); 
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				     "column_label_window"));
}


void
on_about1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *about_window = create_aboutdialog();
   add_coot_references_button(about_window);
   gtk_widget_show(about_window);   
}


void
on_clipping_button_clicked             (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				    "clipping_window"));
}




/* Density Radius Window Widgets  */

void
on_density_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkEntry      *entry; 

/*    entry_char_type *text; */
   const char *text;
   int imol = 0;		/*  FIXME */

   entry = (GTK_ENTRY(lookup_widget(GTK_WIDGET(button), "entry1")));
   text = gtk_entry_get_text(entry);

 /*   printf("We found entry text: %s\n", text); */

   set_density_size_from_widget(text);

   /* Now the increment of the iso level entry */

   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				   "iso_level_increment_entry"));

   text = gtk_entry_get_text(entry);

   set_iso_level_increment_from_text(text, imol);  /* imol is ignored. */
                             /* we change the iso level for all maps */

   /* As above, except a difference map */

   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				   "diff_map_iso_level_increment_entry"));

   text = gtk_entry_get_text(entry);

   set_diff_map_iso_level_increment_from_text(text, imol); /* imol is ignored. */
                                       /* we change the iso level for all maps */


   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				   "map_sampling_rate_entry"));

   text = gtk_entry_get_text(entry);

   set_map_sampling_rate_text(text);
                               


 /* Goodbye Mr Widget */
   
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				     "global_map_properties_window"));
}
/* In the menubar, Edit Density size has been selected. */
void
on_density_size1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *density_window;
   GtkWidget *entry; 
   GtkWidget *checkbutton;
   char *text;
   int imol = 0;		/* FIXME */
   
   density_window = create_global_map_properties_window();
   entry = lookup_widget(density_window, "entry1");

   text = get_text_for_density_size_widget(); /* const gchar *text */
   
  /* Question to self: is this GTK_ENTRY necessary 
     if entry is an GtkEntry? */
   gtk_entry_set_text(GTK_ENTRY(entry), text);

   free (text); 
  
 /* Now the iso level increment entry  */

   entry = lookup_widget(density_window, "iso_level_increment_entry");
   text = get_text_for_iso_level_increment_entry(imol);

   gtk_entry_set_text(GTK_ENTRY(entry), text);

 /* Now the iso level for the differenece map increment entry  */

   entry = lookup_widget(density_window, "diff_map_iso_level_increment_entry");
   text = get_text_for_diff_map_iso_level_increment_entry(imol);

   gtk_entry_set_text(GTK_ENTRY(entry), text);

   /* Now the map rate multiplier: */
   entry = lookup_widget(density_window, "map_sampling_rate_entry");
   text = get_text_for_map_sampling_rate_text();

   gtk_entry_set_text(GTK_ENTRY(entry), text);

   checkbutton = lookup_widget(density_window, "map_dynamic_map_sampling_checkbutton");
   set_map_dynamic_map_sampling_checkbutton(checkbutton);
   checkbutton = lookup_widget(density_window, "map_dynamic_map_size_display_checkbutton");
   set_map_dynamic_map_display_size_checkbutton(checkbutton);

 /* Show the widget */
   gtk_widget_show(density_window);
   
}


void
on_density_cancel_clicked              (GtkButton       *button,
                                        gpointer         user_data)
{

   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				     "global_map_properties_window"));

}



gboolean
on_hscale1_button_press_event          (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data)
{
   
   printf("button press\n");
   return FALSE;
}


gboolean
on_hscale1_button_release_event        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data)
{

   printf("button release\n");
   return FALSE;
}


void
on_fps1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *fps_window = create_fps_window(); 
   GtkButton       *button; 

   if ( get_fps_flag() == 1) { 
      button = GTK_BUTTON(lookup_widget(fps_window,
					"radiobutton1"));
   } else { 
      button = GTK_BUTTON(lookup_widget(fps_window,
					"radiobutton2"));
   } 

   
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

   gtk_widget_show(fps_window);

}


void
on_fps_window_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
 /* Lets look up the buttons and see if they are active. If the yes
    button is active then set the flag to 1, if no is active, set it
    to 0.  */
   
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "radiobutton1"))))
      set_fps_flag(1);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "radiobutton2"))))
      set_fps_flag(0);

   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), "fps_window"));
   
}


void
on_active_map_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "active_map_radiobutton_yes"))))
      set_active_map_drag_flag(1);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "active_map_radiobutton_no"))))
      set_active_map_drag_flag(0);

   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				    "active_map_window"));

}


void
on_dragged_map1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *active_map_window = create_active_map_window(); 
   GtkButton       *button; 


   if ( get_active_map_drag_flag() == 1 ) { 
      button = GTK_BUTTON(lookup_widget(active_map_window, 
					"active_map_radiobutton_yes")); 
   } else { 
      button = GTK_BUTTON(lookup_widget(active_map_window, 
					"active_map_radiobutton_no")); 
   }

   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


   gtk_widget_show(active_map_window); 


}


void
on_map_colour1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "rotamer_analysis1");
  if (menu) { 
    add_on_map_colour_choices(menu);
  } else { 
    printf("ERROR:: failed to get menu in on_map_colour1_activate\n");
  }
}    

void
on_show_symmetry1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *show_symm_window;

   show_symm_window = wrapped_create_show_symmetry_window();

/* Now show the popup eventually */

   gtk_widget_show(show_symm_window);

}


void
on_show_symmetry_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkEntry      *entry; 
   GtkWidget     *checkbutton;
   const char *text;

/* Show Symmetry Radiobuttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "show_symmetry_yes_radiobutton"))))
      set_show_symmetry_master(1);
				  

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "show_symmetry_no_radiobutton"))))
      set_show_symmetry_master(0);

/* Symmetry Radius Entry */
   
   entry = (GTK_ENTRY(lookup_widget(GTK_WIDGET(button), 
				    "symmetry_radius_entry")));

   text = gtk_entry_get_text(entry);

   /* printf("ok button symmetry: radius text: %s\n",text); */

   set_symmetry_size_from_widget(text);


/* Show UnitCell Radiobuttons */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "unit_cell_yes_radiobutton"))))
      set_show_unit_cells_all(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "unit_cell_no_radiobutton"))))
      set_show_unit_cells_all(0);

/* The Symmetry Colour Checkbutton */

/*    checkbutton = lookup_widget(GTK_WIDGET(button),  */
/* 			       "show_symmetry_molecule_rotate_colour_map_checkbutton"); */

/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_molecule_rotate_colour_map(1); */
/*    else */
/*      set_symmetry_molecule_rotate_colour_map(0); */

   
/* The Symmetry Colour by Symop Checkbutton */

/*    checkbutton = lookup_widget(GTK_WIDGET(button),  */
/* 			       "show_symmetry_colour_by_symop_checkbutton"); */
/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_colour_by_symop(1); */
/*    else */
/*      set_symmetry_colour_by_symop(0); */


/* The Symmetry Whole Chain Checkbutton */

/*    checkbutton = lookup_widget(GTK_WIDGET(button),  */
/* 			       "show_symmetry_whole_molecule_checkbutton"); */
/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_whole_chain(1); */
/*    else */
/*      set_symmetry_whole_chain(0); */

/* The Expanded Atom Label Checkbutton */

   checkbutton = lookup_widget(GTK_WIDGET(button), 
			       "show_symmetry_expanded_labels_checkbutton");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) 
     set_symmetry_atom_labels_expanded(1);
   else
     set_symmetry_atom_labels_expanded(0);

/* Goodbye Mr Widget */
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				    "show_symmetry_window"));

}


void
on_show_symmetry_apply_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  /* Just a copy of the above funtion, without the destroy. */

   GtkEntry      *entry; 
   GtkWidget     *checkbutton;
   const char *text;

/* Show Symmetry Radiobuttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "show_symmetry_yes_radiobutton"))))
      set_show_symmetry_master(1);
				  

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "show_symmetry_no_radiobutton"))))
      set_show_symmetry_master(0);

/* Symmetry Radius Entry */
   
   entry = (GTK_ENTRY(lookup_widget(GTK_WIDGET(button), 
				    "symmetry_radius_entry")));

   text = gtk_entry_get_text(entry);

   /* printf("ok button symmetry: radius text: %s\n",text); */

   set_symmetry_size_from_widget(text);


/* Show UnitCell Radiobuttons */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "unit_cell_yes_radiobutton"))))
      set_show_unit_cells_all(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "unit_cell_no_radiobutton"))))
      set_show_unit_cells_all(0);

/* The Symmetry Colour Checkbutton */

/*    checkbutton = lookup_widget(GTK_WIDGET(button),  */
/* 			       "show_symmetry_molecule_rotate_colour_map_checkbutton"); */

/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_molecule_rotate_colour_map(1); */
/*    else */
/*      set_symmetry_molecule_rotate_colour_map(0); */

   
/* The Symmetry Colour by Symop Checkbutton */

/*    checkbutton = lookup_widget(GTK_WIDGET(button),  */
/* 			       "show_symmetry_colour_by_symop_checkbutton"); */
/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_colour_by_symop(1); */
/*    else */
/*      set_symmetry_colour_by_symop(0); */


/* The Symmetry Whole Chain Checkbutton */

/*    checkbutton = lookup_widget(GTK_WIDGET(button),  */
/* 			       "show_symmetry_whole_molecule_checkbutton"); */
/*    if (GTK_TOGGLE_BUTTON(checkbutton)->active)  */
/*      set_symmetry_whole_chain(1); */
/*    else */
/*      set_symmetry_whole_chain(0); */

/* The Expanded Atom Label Checkbutton */

   checkbutton = lookup_widget(GTK_WIDGET(button), 
			       "show_symmetry_expanded_labels_checkbutton");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) 
     set_symmetry_atom_labels_expanded(1);
   else
     set_symmetry_atom_labels_expanded(0);
}



void
on_symmetry_colour_patch_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
   GtkWidget *colorseldlg;
   gdouble *colour;
   GtkColorSelection *colorsel;

   // GTK-FIXME fix the colours

/*    colorseldlg = create_symmetry_colour_selection_window(); */

/*    colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(colorseldlg)->colorsel); */
/*    colour = get_symmetry_bonds_colour(0); */
/*    gtk_color_selection_set_color(colorsel, colour); */

/*    gtk_widget_show(colorseldlg);  */
}



void
on_about_ok_button_clicked             (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button),
				    "about_window"));
}


void
on_anisotropic_atoms1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *aniso_window;
   GtkWidget *entry; 
   char *text;
   GtkButton       *button; 
   GtkScale *hscale;
   GtkAdjustment *adjustment;
   float hscale_initial;
 
   aniso_window = create_aniso_window();   

/* Show Aniso Radiobuttons */
   if (get_show_aniso() == 1) { 
      button = GTK_BUTTON(lookup_widget(aniso_window, 
					"show_aniso_yes_radiobutton"));
   } else { 
      button = GTK_BUTTON(lookup_widget(aniso_window,
					"show_aniso_no_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


/* Set Limit Radiobuttons */
   if (get_show_limit_aniso() == 1) { 
      button = GTK_BUTTON(lookup_widget(aniso_window,			
					"limit_display_radius_yes_radiobutton"));
   } else { 
      button = GTK_BUTTON(lookup_widget(aniso_window,			
					"limit_display_radius_no_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


/* The Aniso Display Limit Entry */

   entry = lookup_widget(aniso_window, "limit_display_radius_entry");

   text = get_text_for_aniso_limit_radius_entry();
   gtk_entry_set_text(GTK_ENTRY(entry), text);

   free(text);

/* The Probability Hscale */

   hscale = GTK_SCALE(lookup_widget(GTK_WIDGET(button),
				    "aniso_probability_hscale"));
   
   hscale_initial = get_aniso_probability();

   adjustment = GTK_ADJUSTMENT 
         (gtk_adjustment_new(hscale_initial, 0.0, 110.0, 0.01, 4.0, 10.)); 

   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   g_signal_connect (G_OBJECT(adjustment), "value_changed",
		       G_CALLBACK(aniso_probability_adjustment_changed), 
		       NULL);


/* And finally show the widget */
   gtk_widget_show(aniso_window);

}


void
on_show_aniso_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkEntry      *entry; 
   const char *text;

/* Limit Display Atoms? */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button),
								    "limit_display_radius_yes_radiobutton"
								    ))))
      set_limit_aniso(1);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button),
								    "limit_display_radius_no_radiobutton"
								    ))))
      set_limit_aniso(0);


/* Limit Display Radius Entry */

   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				   "limit_display_radius_entry"));

   text = gtk_entry_get_text(entry); 

   set_aniso_limit_size_from_widget(text);

 /* Show Aniso Radiobuttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button),
								    "show_aniso_yes_radiobutton"))))

      set_show_aniso(1);	/* model, state */

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "show_aniso_no_radiobutton"))))
      set_show_aniso(0); 
}


void
on_show_aniso_close_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
/* Goodbye Widget */
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button),"aniso_window"));
}

void aniso_probability_adjustment_changed(GtkAdjustment *adj, GtkWidget *window) { 

  set_aniso_probability(gtk_adjustment_get_value(adj));
}


void on_smooth_scrolling_window_ok_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{
   GtkEntry *entry;
   const char *text;

/* Show Smooth Scrolling Radio Buttons */
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "smooth_scroll_yes_radiobutton"))))
       
      set_smooth_scroll_flag(1);
   
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "smooth_scroll_no_radiobutton"))))
      
      set_smooth_scroll_flag(0);

/* Smooth Scroll Distance Limit */

   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				   "smooth_scroll_limit_entry"));
   text = gtk_entry_get_text(entry);

   set_smooth_scroll_limit_str(text);


 /*  Smooth Scroll Steps */

   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				   "smooth_scroll_steps_entry"));

   text = gtk_entry_get_text(entry);

   set_smooth_scroll_steps_str(text);


   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button),"smooth_scroll_window"));

}


void
on_recentring1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *smooth_scroll_window;
   GtkWidget *entry;
   char *text;
   GtkButton *button;

   smooth_scroll_window = create_smooth_scroll_window();
   
/* Show Smooth Scrolling Radio Buttons */
   if (get_smooth_scroll() == 1) { 

      button = GTK_BUTTON(lookup_widget(smooth_scroll_window, 
					"smooth_scroll_yes_radiobutton"));
   } else { 
      button = GTK_BUTTON(lookup_widget(smooth_scroll_window,
					"smooth_scroll_no_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);


/* Smooth Scroll Distance Limit */

   entry = lookup_widget(smooth_scroll_window, "smooth_scroll_limit_entry");

   text = get_text_for_smooth_scroll_limit();
   gtk_entry_set_text(GTK_ENTRY(entry), text);

   free(text);

/*  Smooth Scroll Steps */


   entry = lookup_widget(smooth_scroll_window, "smooth_scroll_steps_entry");

   text = get_text_for_smooth_scroll_steps();

   gtk_entry_set_text(GTK_ENTRY(entry), text);

   free(text);


/* And finally show the widget */
   gtk_widget_show(smooth_scroll_window);

}


void
on_font_size1_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *font_size_window;

   GtkButton *button = 0;

   font_size_window = create_font_size_window();

/* The Font Size RadioButtons */

   if (get_font_size() == 1) { 
            button = GTK_BUTTON(lookup_widget(font_size_window, 
					"font_size_small_radiobutton"));
   }
   if (get_font_size() == 2) { 
            button = GTK_BUTTON(lookup_widget(font_size_window, 
					"font_size_medium_radiobutton"));
   }
   if (get_font_size() == 3) { 
            button = GTK_BUTTON(lookup_widget(font_size_window, 
					"font_size_large_radiobutton"));
   }
   gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);

/* show widget */
   gtk_widget_show(font_size_window);
}

void
on_font_size_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{   

/* The Font Size RadioButtons */
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "font_size_small_radiobutton"))))
       
      set_font_size(1);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "font_size_medium_radiobutton"))))
       
      set_font_size(2);

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "font_size_large_radiobutton"))))
       
      set_font_size(3);

/* goodbye widget */
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button),"font_size_window"));

}


/* not used? */
void
on_greer_skeleton1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{


}

/* not used? */
void
on_cowtan_foadi_skeleton1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_these_are_placeholders1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_mapcolourmap1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{


}


void
on_map1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   printf("ping the colourchange widget for map1\n"); 
}


void
on_map2_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

void
on_attach_scroll_wheel_to_which_map_1_activate (GtkMenuItem     *menuitem,
						gpointer         user_data)
{
  /* it doesn't matter what we pass back, a lookup is done to find the
     submenu */
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "attach_scroll_wheel_to_which_map_1");
  if (menu) { 
    add_on_map_scroll_wheel_choices(menu);
  } else { 
    printf("ERROR:: failed to get menu in on_attach_scroll_wheel_to_which_map_1_activate\n");
  }

}


void
on_greer_on_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   printf("Making skel_greer....\n"); 
   skel_greer_on(); 

}


void
on_greer_off_activate                  (GtkMenuItem     *menuitem,
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
void
on_foadi_on_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_foadi_off_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
/*    printf("skel_foadi OFF\n");  */

/*   skel_foadi_off(); dead */ 

}


void
on_open_map1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *map_name_fileselection1; 
  GtkWidget *filter_button;
  GtkWidget *sort_button;
   
   map_name_fileselection1 = coot_map_name_chooser(); 

   gtk_widget_show (map_name_fileselection1);
   add_is_difference_map_checkbutton(map_name_fileselection1);
   add_ccp4i_project_optionmenu(map_name_fileselection1, COOT_MAP_FILE_SELECTION);
   filter_button = add_filename_filter_button(map_name_fileselection1,
					      COOT_MAP_FILE_SELECTION);
   sort_button = add_sort_button_fileselection(map_name_fileselection1); 
   set_directory_for_fileselection(map_name_fileselection1);

   push_the_buttons_on_fileselection(filter_button, sort_button, 
				     map_name_fileselection1);

   set_file_selection_dialog_size(map_name_fileselection1);
   gtk_widget_show (map_name_fileselection1);

   push_the_buttons_on_fileselection(filter_button, sort_button, 
				     map_name_fileselection1);

}




void
on_cancel_button_map_name_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{

   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				     "map_name_fileselection1"));
   
}


void
on_vt_flat1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   vt_surface(VT_FLAT); 

}


void
on_vt_spherical_surface1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   vt_surface(VT_SPHERICAL); 


}


void
on_skeleton_colour1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *col_sel_window;
  
   col_sel_window = create_skeleton_colour_selection_window();
   gtk_widget_show(col_sel_window); 


}



void
on_phs_info_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   phs_pdb_cell_symm(); /* which runs
			   create_phs_coordinates_fileselection() */

   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				    "phs_info_box"));

}


void
on_phs_info_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				    "phs_info_box"));

}


void
on_cancel_phs_coord_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{

  gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				     "phs_coordinates_fileselection"));


}


void
on_map_and_mol_control1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
}

void
on_display_only_active1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  /* display only the active mol and the refinement map */
  display_only_active();

}





void
on_go_to_atom_apply_button_clicked (GtkButton       *button,
				    gpointer         user_data)
{

  GtkWidget *widget; 

  widget = lookup_widget(GTK_WIDGET(button), "goto_atom_window");
  
  apply_go_to_atom_from_widget(widget);


}


void
on_go_to_atom_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget; 
  
  widget     = lookup_widget(GTK_WIDGET(button), "goto_atom_window");

  unset_go_to_atom_widget();
  gtk_widget_destroy(widget);	/* There is something that had been
				   added to the Go To Atom window that
				   is not a widget.  It has been
				   cleared perhaps but not destroyed.
				   I don't think that it is the dialog
				   itself.  The problem does not
				   happen in the GTK1 path.  */
}


void
on_go_to_atom1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
     GtkWidget *widget;

     /* fills the window too */
     widget = wrapped_create_goto_atom_window(); 

				/* now we need to fill the entry boxes
				   with default vaules and the option
				   menu according to molecules that
				   have coordinates. */


     /* and finally show the widget: */
     gtk_widget_show(widget); 
}


void
on_go_to_atom_next_residue_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{

/*   GtkEntry *entry;  */
/*   GtkEntry *residue_entry;  */
/*   gchar *chain_str; */
/*   gchar *res_str;  */
/*   gchar *atom_name_str;  */
 
/*   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),  */
/* 				  "go_to_atom_chain_entry")); */
/*   chain_str = gtk_entry_get_text(entry); */

/*   residue_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),  */
/* 					  "go_to_atom_residue_entry")); */
/*   res_str = gtk_entry_get_text(residue_entry); */

/*   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),  */
/* 				  "go_to_atom_atom_name_entry")); */
/*   atom_name_str = gtk_entry_get_text(entry); */

/*   goto_next_atom_maybe(chain_str, res_str, atom_name_str, residue_entry); */

  GtkWidget *window = lookup_widget(GTK_WIDGET(button), "goto_atom_window");
  goto_next_atom_maybe_new(window);
}


void
on_go_to_atom_previous_residue_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{

/*   GtkEntry *entry;  */
/*   GtkEntry *residue_entry;  */
/*   gchar *chain_str; */
/*   gchar *res_str;  */
/*   gchar *atom_name_str;  */
 
/*   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),  */
/* 				  "go_to_atom_chain_entry")); */
/*   chain_str = gtk_entry_get_text(entry); */

/*   residue_entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),  */
/* 					  "go_to_atom_residue_entry")); */
/*   res_str = gtk_entry_get_text(residue_entry); */

/*   entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),  */
/* 				  "go_to_atom_atom_name_entry")); */
/*   atom_name_str = gtk_entry_get_text(entry); */
/*   goto_prev_atom_maybe(chain_str, res_str, atom_name_str, residue_entry);  */

  GtkWidget *window = lookup_widget(GTK_WIDGET(button), "goto_atom_window");
  if (! window) 
     printf("ERROR:: in on_go_to_atom_previous_residue_button_clicked NULL window\n");
  goto_previous_atom_maybe_new(window);
}


void
on_skeleton_box_radius1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  GtkWidget *entry; 
  gchar *text; 

  widget = create_skeletonization_box_radius_window(); 

  text = get_text_for_skeleton_box_size_entry(); 
  entry = lookup_widget(widget, "skeleton_box_size_entry"); 
  gtk_entry_set_text(GTK_ENTRY(entry), text); 
  g_free(text); 
  
	/* show the widget */
  gtk_widget_show(widget); 

}


void
on_skeletonization_level1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  GtkWidget *entry; 
  gchar *text; 
  
  widget = create_skeletonization_level_window(); 

				/* Fill the entry: */
  text = get_text_for_skeletonization_level_entry(); 
  entry = lookup_widget(widget, "skeleton_level_entry"); 
  gtk_entry_set_text(GTK_ENTRY(entry), text); 

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


void
on_skel_box_radius_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkEntry *entry;
  const char    *txt; 

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				  "skeleton_box_size_entry")); 

  txt = gtk_entry_get_text(entry); 

  set_skeleton_box_size_from_widget(txt); 

  gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				   "skeletonization_box_radius_window")); 
  
}


void
on_skel_box_radius_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{

  gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				   "skeletonization_box_radius_window")); 
  
}


void
on_skeletonization_level_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{

  GtkEntry *entry; 
  const char *txt; 

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				  "skeleton_level_entry")); 

  txt = gtk_entry_get_text(entry); 

  set_skeletonization_level_from_widget(txt); /* does a skeleton
						 update and redraw */

  gtk_widget_destroy(lookup_widget(GTK_WIDGET(button),
				   "skeletonization_level_window")); 

}


void
on_skeletonization_level_apply_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{

  
  GtkEntry *entry; 
  const char *txt; 

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button),
				  "skeleton_level_entry")); 

  txt = gtk_entry_get_text(entry); 

  set_skeletonization_level_from_widget(txt); /* does a skeleton
						 update and redraw */


}


void
on_skeletonization_level_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data)
{

  gtk_widget_destroy(lookup_widget(GTK_WIDGET(button),
				   "skeletonization_level_window")); 
}


void
on_autobuild_ca_on_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

/*   autobuild_ca_on();  - removed to junk now */

}


void
on_autobuild_ca_off_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
    autobuild_ca_off(); 
}


void
on_display_control_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"display_control_window_glade");
   GtkWidget *maps_vbox = 0;
   GtkWidget *molecules_vbox = 0;
   GtkWidget *pane = 0;

  reset_graphics_display_control_window(); /* Needed! (also resets the scroll group) */
  if (w) { 
     store_window_position(COOT_DISPLAY_CONTROL_WINDOW, w);
     maps_vbox      = lookup_widget(w, "display_map_vbox"); 
     molecules_vbox = lookup_widget(w, "display_molecule_vbox"); 
     pane           = lookup_widget(w, "display_control_vpaned");
     /* store the size, actually */
     if (maps_vbox)
	store_window_position(COOT_DISPLAY_CONTROL_MAPS_VBOX, maps_vbox);
     if (molecules_vbox)
	store_window_position(COOT_DISPLAY_CONTROL_MOLECULES_VBOX,
			      molecules_vbox);
     if (pane) 
	store_window_position(COOT_DISPLAY_CONTROL_PANE, pane);

     gtk_widget_destroy(w);
  }
}


void
on_display_control_window_glade_destroy (GtkWidget       *object,
					 gpointer         user_data)
{
  reset_graphics_display_control_window(); /* (also resets the scroll group) */
}


void
on_rotation_centre_size_ok_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{

				/* pass back the value from the entry */
  GtkEntry *entry; 
  const char *text; 

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button), 
				  "rotation_centre_cube_size_entry")); 

  text = gtk_entry_get_text(entry); 

  set_rotation_centre_size_from_widget(text); 


  gtk_widget_destroy(lookup_widget(GTK_WIDGET(button), 
				   "rotation_centre_cube_size_window")); 

}


void
on_pink_pointer_size1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *window; 
  GtkWidget *entry; 
  gchar *text; 

  window = create_rotation_centre_cube_size_window(); 

  entry = lookup_widget(GTK_WIDGET(window), "rotation_centre_cube_size_entry");

  text = get_text_for_rotation_centre_cube_size(); 

  gtk_entry_set_text(GTK_ENTRY(entry), text); 
  

  gtk_widget_show(window);
  free(text); 

}


void
on_phs_cell_choice_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window; 
   GtkWidget *info_window; 
   int i; 

   /* messing about with string variables */
   gchar *widget_name; 
   gchar *tmp_name; 

   widget_name = (gchar *) malloc(25); /* freed */

   window = lookup_widget(GTK_WIDGET(button), "phs_cell_choice_window");

   for (i=0; i< graphics_n_molecules(); i++) {

      if (has_unit_cell_state(i)) { 

	 strcpy(widget_name, "phs_cell_radiobutton_"); 
	 tmp_name = widget_name + strlen(widget_name); 
	 snprintf(tmp_name, 3, "%-d", i); 

	 if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), widget_name)))) { 
	    printf("proceeding with phs reading using cell from molecule %d.\n", i);
	
	    read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(i); 
	    break; 
	 }
      }
   }
   free(widget_name);

   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "phs_cell_none_radiobutton")))) { 
      printf("special value for none for phs_cell radiobuton active\n");
      info_window = create_phs_info_box();
      gtk_widget_show(info_window); 
   } 
   gtk_widget_destroy(window); 
}


void on_phs_cell_choice_cancel_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{
  GtkWidget *window; 

   printf("Cancel trying to read a phs file.\n");

  window = lookup_widget(GTK_WIDGET(button), "phs_cell_choice_window"); 
  gtk_widget_destroy(window); 
}

void
on_test_thing1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  test_fragment(); 
}


/* void */
/* on_regularize1_activate                (GtkMenuItem     *menuitem, */
/*                                         gpointer         user_data) */
/* { */
/*   do_regularize(); */
/* } */


void
on_scripting_window_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  /* What does this function do!? */

#if (GTK_MAJOR_VERSION == 1) 
  post_scheme_scripting_window();
#endif

}


void
on_scheme_window_close_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				    "scheme_window");
  gtk_widget_destroy(window);

}


void
on_get_pdb_using_code1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *window; 
  int n = 1;
  window = create_accession_code_window(); 
  g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
  gtk_widget_show(window); 
}


void
on_get_pdb_and_sf_using_code1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *window; 
  int n = 2;
  window = create_accession_code_window(); 
  g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
  gtk_widget_show(window); 
}

void
on_fetch_pdb_and_map_using_pdbredo1_activate
                                        (GtkMenuItem     *menuitem,
					 gpointer         user_data) {

  GtkWidget *window; 
  int n = 3;
  window = create_accession_code_window(); 
  g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
  gtk_widget_show(window); 

}


gboolean on_accession_code_entry_key_press_event (GtkWidget       *widget,
						  GdkEventKey     *event,
						  gpointer         user_data)
{

 /* go somewhere if keypress was a carriage return  */

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) { 
    handle_get_accession_code(widget); 
  } 

  return FALSE;
}

void
on_accession_code_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data) {
  
  GtkWidget *entry = lookup_widget(GTK_WIDGET(button), "accession_code_entry");
  handle_get_accession_code(entry); 
}



void
on_ramachandran_plot1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  add_on_rama_choices();
}


void
on_dynarama_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				     "dynarama_window"); 
   int imol = get_mol_from_dynarama(window); /*  return -9999 on edit rama window */

/* No need to free dynarama plot data here, because it is done on
   window destroy callback */
/*     if (imol >= 0) */
/*       set_dynarama_is_displayed(0, imol); */ /*  which frees/deletes the */
     /*  			             memory of the user data. */


   if (imol == -9999) 
     accept_phi_psi_moving_atoms();

   gtk_widget_destroy(window); 
}

void
on_dynarama_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
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

   gtk_widget_destroy(window); 

   /* maybe we should also destroy the edit backbone torsion dialog, if it exists. */

}

void
on_dynarama_window_destroy             (GtkWidget       *object,
                                        gpointer         user_data)
{
  /* Maybe object is window? */
   GtkWidget *window = lookup_widget(GTK_WIDGET(object),
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



void
on_background_black1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  set_background_colour(0.0,0.0,0.0);
  
}


void
on_background_white1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  set_background_colour(1.0,1.0,1.0); 
}


void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
				/* Not used any more */
}


void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window;

   execute_get_mols_ligand_search(GTK_WIDGET(button));
				/* which then runs
				   execute_ligand_search */
   window = lookup_widget(GTK_WIDGET(button), "find_ligand_dialog");
   free_ligand_search_user_data(GTK_WIDGET(button));
   gtk_widget_destroy(window);
}


void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button), "find_ligand_dialog");
   free_ligand_search_user_data(GTK_WIDGET(button));
   gtk_widget_destroy(window);
}


void
on_find_ligand_many_atoms_continue_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				     "find_ligand_many_atoms_dialog");
   GtkWidget *find_ligand_dialog = (GtkWidget *) g_object_get_data(G_OBJECT(window), "dialog");

/* Needed at all, the if? */
#if defined USE_GUILE && defined USE_PYTHON
   execute_ligand_search();
#endif
   gtk_widget_destroy(window);
   gtk_widget_destroy(find_ligand_dialog); 
}


void
on_find_ligand_many_atoms_cancel_button_clicked (GtkButton       *button,
						 gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				     "find_ligand_many_atoms_dialog");
   gtk_widget_destroy(window);
}


void
on_prune_and_colour1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  
  do_skeleton_prune();

}

void
on_residue_info1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_residue_info_dialog();
}

/* the model fit refine menu item was activated */
void
on_model_refine_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *window_or_frame;
   window_or_frame = wrapped_create_model_fit_refine_dialog();
   gtk_widget_show(window_or_frame);
}




void
on_model_refine_dialog_dismiss_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{

   GtkWidget *hbox   = lookup_widget(GTK_WIDGET(button), "model_fit_refine_dialog_vbox");
   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "model_refine_dialog");

   clear_pending_picks();
   normal_cursor();
   dialog = close_model_fit_dialog(hbox);
   store_window_position(COOT_MODEL_REFINE_DIALOG, dialog);
   /* was it a top-level?  If so, kill off the top-level. */
   if (dialog) 
      gtk_widget_destroy(dialog);
}



void
on_save_coordinates1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  GtkWidget *combobox; 
  GCallback callback_func = G_CALLBACK(save_molecule_coords_button_select);
  int imol = first_coords_imol();
  int imol_unsaved = first_unsaved_coords_imol();
  if (imol_unsaved != -1) 
    imol = imol_unsaved;
  set_save_molecule_number(imol); /* set *save* molecule number */

  widget = create_save_coords_dialog(); 

  combobox = lookup_widget(GTK_WIDGET(widget), "save_coordinates_combobox");

  if (combobox) {

/*     fill_option_menu_with_coordinates_options_unsaved_first(option_menu, callback_func, imol); */
    fill_combobox_with_coordinates_options(combobox, callback_func, imol);
    set_transient_and_position(COOT_UNDEFINED_WINDOW, widget);
    gtk_widget_show(widget);
    gtk_window_present(GTK_WINDOW(widget));
  } else {
    printf("bad combobox!\n");
  }
}


void
on_save_coords_dialog_save_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *combobox = lookup_widget(GTK_WIDGET(button), "save_coordinates_combobox");
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "save_coords_dialog");
  GtkWidget *chooser;
  int imol;
  if (! combobox) {
    printf("on_save_coords_dialog_save_button_clicked: bad combobox\n");
  } else {
    imol = my_combobox_get_imol(GTK_COMBO_BOX(combobox));
    chooser = coot_save_coords_chooser();
    g_object_set_data(G_OBJECT(chooser), "imol", GINT_TO_POINTER(imol));
    gtk_widget_show(chooser);
  }
  gtk_widget_destroy(dialog);

}


void
on_save_coord_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  char *stuff;
  GtkWidget *fileselection = lookup_widget(GTK_WIDGET(button), "save_coords_fileselection1");

  widget = lookup_widget(GTK_WIDGET(button), "save_coords_fileselection1");
  save_directory_for_saving_from_fileselection(fileselection);
  stuff = g_object_get_data(G_OBJECT(widget), "stuff"); /* probably needs fixing GTK-FIXME */
  save_coordinates_using_widget(widget);
  free(stuff);
  gtk_widget_destroy(widget);
}


void
on_save_coords_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  char *stuff;

  widget = lookup_widget(GTK_WIDGET(button),
			 "save_coords_fileselection1");
  stuff = g_object_get_data(G_OBJECT(widget), "stuff");
  free(stuff);
  gtk_widget_destroy(widget);

}


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


void
on_model_refine_dialog_fixed_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_fixed_atom_dialog();
  gtk_widget_show(w);
}


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


void
on_model_refine_dialog_refine_params_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = wrapped_create_refine_params_dialog();
  gtk_widget_show(widget);
}


void
on_refine_params_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button),
				    "refine_params_dialog");
  GtkWidget *entry = lookup_widget(GTK_WIDGET(button), 
				   "refine_params_weight_matrix_entry");
  if (entry) { 
    set_refinement_weight_from_entry(entry);
  }

  gtk_widget_destroy(widget);
}


void
on_goto_atom_window_destroy            (GtkWidget       *object,
                                        gpointer         user_data)
{
/*   printf("on_goto_atom_window_destroy is executed. unsetting the static in graphics_info\n"); */

/* we can get here from a WM dialog kill (as well as the conventional
   "Cancel" button callback. Fixes July 12 bug? */
  unset_go_to_atom_widget(); 
}


void
on_distance1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   do_distance_define();
}


void
on_angle1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   do_angle_define();

}


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

void
on_accept_reject_refinement_accept_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "accept_reject_refinement_dialog");

  /* Pressing Return while focus is on the Accept/Reject dialog brings us here. */

  stop_refinement_internal();
  accept_moving_atoms();
  save_accept_reject_dialog_window_position(window);
  set_accept_reject_dialog(NULL);
  gtk_widget_destroy(window);
}


void
on_accept_reject_refinement_reject_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "accept_reject_refinement_dialog");
  save_accept_reject_dialog_window_position(window);
  /*   clear_up_moving_atoms(); done in destroy of the window */

  stop_refinement_internal();
  gtk_widget_destroy(window);
}

void
on_accept_reject_refinement_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{

  /* Pressing Esc while focus is on the Accept/Reject dialog brings us here. */

  /* 20070801 To Fix a crash reported by "Gajiwala, Ketan", we need to
     reset the value for graphics_info_t::accept_reject_dialog (it's
     gone now).  And I suppose that we should clean up (and undisplay)
     the intermediate atoms too.
 */

  set_accept_reject_dialog(NULL);
  stop_refinement_internal();
  /* I want to merely clear the stick restraint, not refine again after I did that */
  clear_atom_pull_restraint_on_accept_reject_destroy();
  clear_up_moving_atoms();
}


/* accept_reject_refinement_docked_accept_button */

void
on_accept_reject_refinement_docked_accept_button_clicked (GtkButton       *button,
							  gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
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


void
on_accept_reject_refinement_docked_reject_button_clicked (GtkButton       *button,
							  gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
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


void
on_model_refine_dialog_find_waters_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

   GtkWidget *widget;

   widget = create_find_waters_dialog();
   
   fill_find_waters_dialog(widget);
   gtk_widget_show(widget);

}


void
on_find_waters_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget;
   
   widget = lookup_widget(GTK_WIDGET(button), 
			  "find_waters_dialog");

   execute_find_waters(GTK_WIDGET(button));
   gtk_widget_destroy(widget);

}


void
on_find_waters_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget;
   
   widget = lookup_widget(GTK_WIDGET(button), 
			  "find_waters_dialog");

   gtk_widget_destroy(widget);
}


void
on_fast_sss_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget;
  widget = lookup_widget(GTK_WIDGET(button), "fast_ss_search_dialog");
  gtk_widget_destroy(widget);

}


void
on_fast_sss_dialog_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
#if (GTK_MAJOR_VERSION > 1)
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

  dialog = lookup_widget(GTK_WIDGET(button), "fast_ss_search_dialog");

  helix_checkbutton   = lookup_widget(dialog, "fast_sss_dialog_helix_checkbutton");
  helix_temp_combobox = lookup_widget(dialog, "fast_sss_dialog_helix_template_combobox");
  helix_noaa_combobox = lookup_widget(dialog, "fast_sss_dialog_helix_no_aa_combobox");
  strand_checkbutton   = lookup_widget(dialog, "fast_sss_dialog_strand_checkbutton");
  strand_temp_combobox = lookup_widget(dialog, "fast_sss_dialog_strand_template_combobox");
  strand_noaa_combobox = lookup_widget(dialog, "fast_sss_dialog_strand_no_aa_combobox");
  radius_checkbutton   = lookup_widget(dialog, "fast_sss_dialog_local_checkbutton");
  radius_combobox = lookup_widget(dialog, "fast_sss_dialog_radius_combobox");

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

  gtk_widget_destroy(dialog);
#endif /* GTK_MAJOR_VERSION */
}


void
on_fast_sss_dialog_citation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
#if (GTK_MAJOR_VERSION > 1)
  GtkWidget *dialog;
  GtkWidget *toolbutton;
  dialog = wrapped_create_coot_references_dialog();
  toolbutton = lookup_widget(dialog, "coot_references_buccaneer_toolbutton");
  fill_references_notebook(GTK_TOOL_BUTTON(toolbutton), COOT_REFERENCE_BUCCANEER);
#else
  g_print("INFO:: sorry not in GTK+ 1.2\n");
#endif /* GTK_MAJOR_VERSION */

}


void
on_environment_distances1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *widget;
   widget = create_environment_distance_dialog();
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
/*    widget = lookup_widget(GTK_WIDGET(button), "environment_distance_dialog"); */
/*    execute_environment_settings(GTK_WIDGET(button)); */
/*    gtk_widget_destroy(widget); */
   
/* } */


void
on_refine_params_use_torsions_checkbutton_toggled (GtkToggleButton *togglebutton, 
						   gpointer         user_data)
{
   GtkWidget *omega_checkbutton = 
      lookup_widget(GTK_WIDGET(togglebutton), 
		    "refine_params_use_peptide_omegas_checkbutton");
   GtkWidget *phi_psi_restraints_vbox = 
      lookup_widget(GTK_WIDGET(togglebutton), 
		    "peptide_torsions_restraints_vbox");

   do_torsions_toggle(GTK_WIDGET(togglebutton));

   /* We don't see these widgets, currently */
   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(omega_checkbutton,       TRUE);
   } else {
      gtk_widget_set_sensitive(omega_checkbutton,       FALSE);
      gtk_widget_set_sensitive(phi_psi_restraints_vbox, FALSE);
      /* no rama restraints */
   }
}

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

				

void
on_refine_params_use_initial_pos_checkbutton_toggled (GtkToggleButton *togglebutton,
						      gpointer         user_data)
{

}


void
on_map_dynamic_map_sampling_checkbutton_toggled (GtkToggleButton *togglebutton,
						 gpointer         user_data)
{
   toggle_dynamic_map_sampling();
}


void
on_map_dynamic_map_size_display_checkbutton_toggled (GtkToggleButton *togglebutton,
						     gpointer         user_data)
{
   toggle_dynamic_map_display_size();
}

/* I mean phi/psi torsions */
void
on_refine_params_use_peptide_torsions_checkbutton_toggled (GtkToggleButton *togglebutton,
							   gpointer         user_data)
{

   GtkWidget *frame = lookup_widget(GTK_WIDGET(togglebutton), 
				    "peptide_torsions_restraints_vbox");
   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(frame, TRUE);
   } else { 
      gtk_widget_set_sensitive(frame, FALSE);
   }
}


void on_import_cif_dictionary1_activate     (GtkMenuItem     *menuitem,
					     gpointer         user_data)
{
  open_cif_dictionary_file_selector_dialog();
}



void
on_cif_dictionary_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *fileselection;

  fileselection = lookup_widget(GTK_WIDGET(button), "cif_dictionary_fileselection");
  gtk_widget_destroy(fileselection);

}


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


void
on_residue_type_chooser_ALA_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   
   gtk_widget_destroy(window);
   do_mutation("ALA", istate);
}


void
on_residue_type_chooser_ARG_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("ARG", istate);

}


void
on_residue_type_chooser_ASN_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("ASN", istate);

}


void
on_residue_type_chooser_ASP_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("ASP", istate);

}


void
on_residue_type_chooser_CYS_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("CYS", istate);

}


void
on_residue_type_chooser_GLN_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("GLN", istate);

}


void
on_residue_type_chooser_GLU_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("GLU", istate);

}


void
on_residue_type_chooser_GLY_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("GLY", istate);

}


void
on_residue_type_chooser_HIS_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("HIS", istate);

}


void
on_residue_type_chooser_ILE_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("ILE", istate);

}


void
on_residue_type_chooser_LEU_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("LEU", istate);

}


void
on_residue_type_chooser_LYS_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("LYS", istate);

}


void
on_residue_type_chooser_MET_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("MET", istate);

}

void
on_residue_type_chooser_MSE_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("MSE", istate);

}



void
on_residue_type_chooser_PHE_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("PHE", istate);

}


void
on_residue_type_chooser_PRO_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("PRO", istate);

}


void
on_residue_type_chooser_SER_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("SER", istate);

}


void
on_residue_type_chooser_THR_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("THR", istate);

}


void
on_residue_type_chooser_TRP_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("TRP", istate);

}


void
on_residue_type_chooser_TYR_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("TYR", istate);

}


void
on_residue_type_chooser_VAL_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "residue_type_chooser_window");
   GtkWidget *stub_button = 
     lookup_widget(window, "residue_type_chooser_stub_checkbutton");
   short int istate = 0;
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(stub_button)))
     istate = 1;
   gtk_widget_destroy(window);
   do_mutation("VAL", istate);

}


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


void
on_model_refine_dialog_rot_trans_by_residue_range_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
    set_model_fit_refine_rotate_translate_zone_label("Rotate/Translate Zone");
    set_rot_trans_object_type(12);
  }

}


void
on_model_refine_dialog_rot_trans_by_chain_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
    set_model_fit_refine_rotate_translate_zone_label("Rotate/Translate Chain");
    set_rot_trans_object_type(13);
  }

}


void
on_model_refine_dialog_rot_trans_by_molecule_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
    set_model_fit_refine_rotate_translate_zone_label("Rotate/Translate Molecule");
    set_rot_trans_object_type(14);
  }

}


void on_rotate_translate_obj_ok_button_clicked (GtkButton       *button,
						gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), 
				    "rotate_translate_obj_dialog");

  // graphics_unsetup_rotate_translate_buttons(widget);
  rot_trans_reset_previous();
  accept_regularizement();
  store_window_position(COOT_ROTATE_TRANSLATE_DIALOG, widget);
  gtk_widget_destroy(widget);
  clear_up_moving_atoms(); // redraw done here

}


void on_rotate_translate_obj_cancel_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), 
				    "rotate_translate_obj_dialog");

 /*   clear_moving_atoms_object(); // redraw done here */
  clear_up_moving_atoms();
  store_window_position(COOT_ROTATE_TRANSLATE_DIALOG, widget);
  gtk_widget_destroy(widget);
  normal_cursor();
}


void
on_run_script1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget = coot_run_script_chooser();
  add_sort_button_fileselection(widget);
  add_filename_filter_button(widget, COOT_SCRIPTS_FILE_SELECTION);
  add_ccp4i_project_optionmenu(widget, 
                               COOT_SCRIPTS_FILE_SELECTION);
  gtk_widget_show(widget);
}



void
on_model_refine_dialog_db_main_togglebutton_toggled (GtkButton       *button,
					       gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
     do_db_main(1);
  else 
     do_db_main(0);
}


void
on_model_refine_dialog_delete_button_clicked (GtkButton       *button,
					      gpointer         user_data)
{
  GtkWidget *widget = wrapped_create_delete_item_dialog();
  gtk_widget_show(widget);
}


void
on_delete_item_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
					    gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(togglebutton), "delete_item_residue_radiobutton"))))
    set_delete_residue_mode();

}


void
on_delete_item_atom_radiobutton_toggled (GtkToggleButton *togglebutton, 
					 gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(togglebutton), "delete_item_atom_radiobutton"))))
    set_delete_atom_mode();
}



void
on_delete_item_sidechain_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_sidechain_mode();
}

void
on_delete_item_sidechain_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data) {

  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_sidechain_range_mode();
}



void
on_delete_item_chain_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_chain_mode();
}




void
on_delete_item_residue_hydrogens_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_delete_residue_hydrogens_mode();

}


void
on_delete_item_water_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_delete_water_mode();
}


/* Old unused function, see new_close_molecules_dialog(), new_close_molecules() */
void
on_close_molecule1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

/* Old Stylee */
/*   GtkWidget *widget = create_close_molecule_dialog(); */
/*   GtkWidget *optionmenu = lookup_widget(widget, "close_molecule_optionmenu"); */
/*    fill_close_option_menu_with_all_molecule_options(optionmenu); */
/*   gtk_widget_show(widget); */

  GtkWidget *widget = wrapped_create_new_close_molecules_dialog();
  gtk_widget_show(widget);
}


/* Old unused function, see new_close_molecules_dialog() */
void
on_close_molecule_close_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "close_molecule_dialog");
  GtkWidget *optionmenu = lookup_widget(window, "close_molecule_optionmenu");
  close_molecule_by_widget(optionmenu);
/*   gtk_widget_destroy(window); */
}


void
on_close_molecule_cancel_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "close_molecule_dialog");
  gtk_widget_destroy(window);

}


void
on_delete_item_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget;

   widget = lookup_widget(GTK_WIDGET(button), "delete_item_dialog"); 
   clear_pending_delete_item();
   clear_pending_picks(); 	/* hmmm.. not sure 20050610 */
   normal_cursor();
   store_window_position(COOT_DELETE_WINDOW, widget);
   store_delete_item_widget(NULL);
   gtk_widget_destroy(widget); 
}


void
on_hints1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *widget = create_hints_dialog();
   gtk_widget_show(widget);
}


/* this is the Apply button now */
void
on_residue_info_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "residue_info_dialog");
   apply_residue_info_changes(widget);
/*    gtk_widget_destroy(widget); not now that it's the Apply button*/

}


void
on_residue_info_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "residue_info_dialog");
   residue_info_release_memory(widget);
   gtk_widget_destroy(widget);
   unset_residue_info_widget();

}

void
on_residue_info_dialog_destroy         (GtkWidget       *object,
                                        gpointer         user_data)
{
   GtkWidget *widget = lookup_widget(GTK_WIDGET(object), "residue_info_dialog");
   residue_info_release_memory(widget);
   clear_residue_info_edit_list();
   unset_residue_info_widget();

}



void
on_hints_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "hints_dialog");
   gtk_widget_destroy(widget);

}


void 
on_model_refine_dialog_find_ligands_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  do_find_ligands_dialog();
}


void
on_rotamer_selection_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button),
				    "rotamer_selection_dialog");
   accept_regularizement();
   clear_moving_atoms_object();
   store_window_position(COOT_ROTAMER_SELECTION_DIALOG, dialog);
   gtk_widget_destroy(dialog);
}


void
on_rotamer_selection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button),
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
   gtk_widget_destroy(dialog);
}



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


gboolean
on_go_to_atom_chain_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{
  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) { 
    apply_go_to_atom_values(lookup_widget(widget, "goto_atom_window"));
  } 
  return FALSE;
}


gboolean
on_go_to_atom_residue_entry_key_press_event (GtkWidget       *widget,
					     GdkEventKey     *event,
					     gpointer         user_data)
{


  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) { 
    apply_go_to_atom_values(lookup_widget(widget, "goto_atom_window"));
  } 
  return FALSE;
}


gboolean
on_go_to_atom_atom_name_entry_key_press_event (GtkWidget       *widget,
					       GdkEventKey     *event,
					       gpointer         user_data)
{

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) { 
    apply_go_to_atom_values(lookup_widget(widget, "goto_atom_window"));
  } 

  return FALSE;
}


void
on_model_refine_dialog_pointer_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   place_atom_at_pointer();
}


void
on_unsaved_changes_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button),
				     "unsaved_changes_dialog");
   gtk_widget_destroy(dialog);
}


void
on_unsaved_changes_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button),
				     "unsaved_changes_dialog");
   gtk_widget_destroy(dialog);
   coot_clear_backup_or_real_exit(0);
}


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


void
on_baton_accept_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  accept_baton_position();
}


void
on_baton_try_again_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  baton_tip_try_another();
}

void
on_baton_tip_previous_button_clicked      (GtkButton       *button,
					   gpointer         user_data) { 

   baton_tip_previous();
} 



void
on_baton_lengthen_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  lengthen_baton();
}

void
on_baton_shorten_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
  shorten_baton();
}

/* dismiss button */
void
on_baton_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button),
				    "baton_dialog");
  try_set_draw_baton(0);
  set_baton_mode(0);		/* if you can't see it, there's no
				   point in trying to move it */
  gtk_widget_destroy(widget);

}


void on_model_refine_dialog_baton_button_clicked (GtkButton       *button,
						  gpointer         user_data)
{
  GtkWidget *widget;
  /* return whether the baton was drawn or not, We don't want a baton
     dialog if the baton is not drawn (e.g. when there is no
     skeletonized map). */
  int state = try_set_draw_baton(1);
  if (state) { 
    widget = create_baton_dialog();
    gtk_widget_show(widget);
  }
}




void
on_use_weights_checkbutton_toggled     (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton),
				  "column_label_window_weights_hbox");
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
    gtk_widget_set_sensitive(GTK_WIDGET(hbox), TRUE);
  } else { 
    gtk_widget_set_sensitive(GTK_WIDGET(hbox), FALSE);
  }

}


void
on_environment_distance_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   toggle_environment_show_distances(togglebutton);
}


void
on_environment_distance_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget;
   widget = lookup_widget(GTK_WIDGET(button), "environment_distance_dialog");
   execute_environment_settings(GTK_WIDGET(button));
   gtk_widget_destroy(widget);

}


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



void
on_coordinates_recentring1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget = create_read_pdb_recentre_dialog();
				/* lookup the toggle widgets here and
				   set the acording to
				   recentre_on_read_pdb()  */
  GtkWidget *yes_radio_button = lookup_widget(widget, 
					      "read_pdb_recentre_yes_radiobutton");
  GtkWidget *no_radio_button = lookup_widget(widget, 
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


void
on_read_pdb_recentre_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), 
				   "read_pdb_recentre_dialog");
  gtk_widget_destroy(widget);

}


void
on_pointer_atom_type_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button),
				    "pointer_atom_type_dialog");
  gtk_widget_destroy(widget);

}


void
on_pointer_atom_type_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button),
				    "pointer_atom_type_dialog");
  GtkToggleButton *tbut;

  GtkWidget *entry = lookup_widget(GTK_WIDGET(button), "pointer_atom_type_other_entry");
  const char *entry_text = gtk_entry_get_text(GTK_ENTRY(entry));


  if (strlen(entry_text) > 0) { 
    place_typed_atom_at_pointer(entry_text);
  } else { 
    /* Adding something here? 
       Remember to change also molecule_class_info_t::add_typed_pointer_atom(). */
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_water"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Water");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_ca"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Ca");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_mg"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Mg");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_na"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Na");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_cl"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Cl");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_br"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("Br");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_so4"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("SO4");
    tbut = GTK_TOGGLE_BUTTON(lookup_widget(GTK_WIDGET(button), "pointer_atom_type_radiobutton_po4"));
    if (gtk_toggle_button_get_active(tbut)) place_typed_atom_at_pointer("PO4");
  }
  /* Recall that the molecule is set by the callback from menu item "activate" */

  gtk_widget_destroy(dialog);
}


void
on_refmac_column_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *frame = lookup_widget(GTK_WIDGET(togglebutton),
				   "column_label_refmac_frame");
  if (gtk_toggle_button_get_active(togglebutton)) {
    gtk_widget_set_sensitive(frame, TRUE);
  } else { 
    gtk_widget_set_sensitive(frame, FALSE);
  }
}


void
on_run_refmac_run_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				    "run_refmac_dialog");
/*   printf("debugging bad window ------------- starting on button click callback.... \n"); */

  execute_refmac(window);
/*   free_memory_run_refmac(window); */
/*   printf("debugging bad window about to refmac window destroy\n"); */
  /* This lookup is for testing/fixing? double running wierdness from
     the scm_catch */
  if (GTK_IS_BUTTON(button)) { 
     window = lookup_widget(GTK_WIDGET(button), "run_refmac_dialog");
     if (window)
	gtk_widget_destroy(window);
  }
/*   printf("debugging bad window done on_run refmac callback...\n"); */
}


void
on_run_refmac_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				    "run_refmac_dialog");
  free_memory_run_refmac(window);
  gtk_widget_destroy(window);

}


void
on_model_refine_dialog_refmac_button_clicked (GtkButton       *button,
					      gpointer         user_data)
{
  
  wrapped_create_run_refmac_dialog();

}


void
on_single_map_properties_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "single_map_properties_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));

  if (is_valid_map_molecule(imol)) { 
    set_contour_by_sigma_step_maybe(window, imol);
    skeletonize_map_single_map_maybe(window, imol);
  } 
  gtk_widget_destroy(window);

}

/* Not sure that this exists any more... */
void
on_single_map_properties_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "single_map_properties_dialog");
  gtk_widget_destroy(window);


}



void
on_single_map_properties_colour_button_clicked (GtkButton       *button,
						gpointer         user_data)
{
  
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "single_map_properties_dialog");
  struct map_colour_data_type *map_colour_data; 
  GtkWidget *col_sel_window;   
  GtkWidget  *colorseldlg;
  GtkColorSelection *colorsel;
  gdouble *colour;
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));

  if (is_valid_map_molecule(imol)) {

#if 0				/* color stuff GTK-FIXME */
      if (1) { 
	/*     printf("pop up the colour chooser for map: %d \n", imol);  */
	map_colour_data = (struct map_colour_data_type *) 
	  malloc(sizeof(struct map_colour_data_type));
	map_colour_data->imol = imol; 
	map_colour_data->imap = 0; /* historical cruft */
	col_sel_window = create_map_colour_selection_window(map_colour_data);  
	/* I think this lookup is unnecessary because the
	   create_map_colour_selection_window returns the dialog
	   (c.f. symmetry_colour_selection)*/
	colorseldlg = GTK_WIDGET(lookup_widget(col_sel_window, "map_colour_selection")); 
	colorsel = GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(colorseldlg)->colorsel);
	colour = get_map_colour(imol);
	gtk_color_selection_set_color(colorsel, colour);
	gtk_widget_show(col_sel_window);

      } else { 
	printf("no imol for the colour chooser for map:\n"); 
      }
#endif
  }
}

#if 0				/* don't know what to do with optionmenu changed */
void
on_run_refmac_phase_input_optionmenu_changed
                                        (GtkOptionMenu   *optionmenu,
                                        gpointer         user_data)
{
  GtkWidget *phases_hbox;
  GtkWidget *hl_hbox;
  GtkWidget *no_labels_checkbutton;
  GtkWidget *twin_checkbutton;
  GtkWidget *sad_extras;
  GtkWidget *fobs_hbox;
  GtkWidget *fpm_hbox;
  int phase_combine_flag;

  phases_hbox = lookup_widget(GTK_WIDGET(optionmenu), "refmac_dialog_phases_hbox");
  hl_hbox     = lookup_widget(GTK_WIDGET(optionmenu), "refmac_dialog_hl_hbox");
  no_labels_checkbutton = lookup_widget(GTK_WIDGET(optionmenu), "run_refmac_nolabels_checkbutton");
  twin_checkbutton = lookup_widget(GTK_WIDGET(optionmenu), "run_refmac_twin_checkbutton");
  sad_extras = lookup_widget(GTK_WIDGET(optionmenu), "run_refmac_sad_extra_hbox");
  fobs_hbox  = lookup_widget(GTK_WIDGET(optionmenu), "refmac_dialog_fobs_hbox");
  fpm_hbox   = lookup_widget(GTK_WIDGET(optionmenu), "refmac_dialog_fpm_hbox");

  phase_combine_flag = get_refmac_phase_input();

  if (phase_combine_flag == 1) {
    gtk_widget_show(phases_hbox);
    gtk_widget_hide(hl_hbox);
  } else {
    if (phase_combine_flag == 2) {
      gtk_widget_hide(phases_hbox);
      gtk_widget_show(hl_hbox);
    } else {
      gtk_widget_hide(phases_hbox);
      gtk_widget_hide(hl_hbox);
    }
  }
  if (refmac_runs_with_nolabels() == 2) {
    if (phase_combine_flag) {
      /* current version of refmac (5.5) doesnt allow phase input for twin or SAD refinement
	 so we de-sensitise and uncheck the buttons */
      gtk_widget_set_sensitive(twin_checkbutton, FALSE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(twin_checkbutton), FALSE);
      /* current version doesnt allow phase input without giving labels */
      gtk_widget_set_sensitive(no_labels_checkbutton, FALSE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(no_labels_checkbutton), FALSE);
    } else {
      gtk_widget_set_sensitive(twin_checkbutton, TRUE);
      gtk_widget_set_sensitive(no_labels_checkbutton, TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(no_labels_checkbutton), TRUE);
    }
    if (phase_combine_flag == 3) {
      /* SAD */
      gtk_widget_show(sad_extras);
      /* change label box from fobs to f+/- as SAD needs this */
      gtk_widget_hide(fobs_hbox);
      gtk_widget_show(fpm_hbox);
    } else {
      gtk_widget_hide(sad_extras);
      /* change label box back from f+/- to fobs for 'normal' refinement */
      gtk_widget_hide(fpm_hbox);
      gtk_widget_show(fobs_hbox);
    }

  }
}
#endif

void
on_run_refmac_tls_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_refmac_use_tls(1);
  } else {
    set_refmac_use_tls(0);
  }
}


void
on_run_refmac_twin_checkbutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *map_optionmenu  = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_map_optionmenu");
  GtkWidget *sad_extras      = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_sad_extra_hbox");
  GtkWidget *mtz_button      = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_map_mtz_radiobutton");
  GtkWidget *mtz_frame       = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_map_mtz_hbox");
  GtkWidget *mtz_file_radiobutton = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_mtz_file_radiobutton");
  GtkWidget *fobs_hbox       = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_fobs_hbox");
  GtkWidget *fiobs_hbox      = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_fiobs_hbox");
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_refmac_use_twin(1);
    /* update the column labels */
    fill_option_menu_with_refmac_labels_options(map_optionmenu);
    /* hide SAD extras, no need to switch off use_sad as done in use_twin */
    gtk_widget_hide(sad_extras);
    /* make the the mtz frame and label insensitive */
    gtk_widget_set_sensitive(mtz_button, FALSE);
    gtk_widget_set_sensitive(mtz_frame, FALSE);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mtz_file_radiobutton), TRUE);
    /* change label box from fobs to fiobs as twin takes both */
    gtk_widget_hide(fobs_hbox);
    gtk_widget_show(fiobs_hbox);
  } else {
    set_refmac_use_twin(0);
    /* update the column labels */
    fill_option_menu_with_refmac_labels_options(map_optionmenu);
    /* hide the twin mtz frame and label but show the 'normal' ones */
    gtk_widget_set_sensitive(mtz_button, TRUE);
    gtk_widget_set_sensitive(mtz_frame, TRUE);
    /* change label box from fiobs to fobs for 'normal' refinement */
    gtk_widget_hide(fiobs_hbox);
    gtk_widget_show(fobs_hbox);
  }

}


void
on_run_refmac_sad_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *twin_checkbutton = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_twin_checkbutton");
  GtkWidget *sad_extras = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_sad_extra_hbox");
  GtkWidget *fobs_hbox  = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_fobs_hbox");
  GtkWidget *fpm_hbox = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_fpm_hbox");
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_refmac_use_sad(1);
    /* de-sensitise the TWIN button, no need to switch off use_twin as done in use_sad */
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(twin_checkbutton), FALSE);
    gtk_widget_set_sensitive(twin_checkbutton, FALSE);
    gtk_widget_show(sad_extras);
    /* fill the entry with 1st existing atom */
    fill_refmac_sad_atom_entry(GTK_WIDGET(togglebutton));
    /* change label box from fobs to f+/- as SAD needs this */
    gtk_widget_hide(fobs_hbox);
    gtk_widget_show(fpm_hbox);
  } else {
    set_refmac_use_sad(0);
    gtk_widget_set_sensitive(twin_checkbutton, TRUE);
    gtk_widget_hide(sad_extras);
    /* change label box back from f+/- to fobs for 'normal' refinement */
    gtk_widget_hide(fpm_hbox);
    gtk_widget_show(fobs_hbox);
  }

}


void
on_run_refmac_map_mtz_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *map_optionmenu  = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_map_optionmenu");
  GtkWidget *active_menu_item;
  if (gtk_toggle_button_get_active(togglebutton)) {
    printf("GTK3 FIXME\n");
  }

}


void
on_run_refmac_mtz_file_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *map_optionmenu  = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_map_optionmenu");
  if (gtk_toggle_button_get_active(togglebutton)) {
    fill_option_menu_with_refmac_file_labels_options(map_optionmenu);
  }
}


void
on_run_refmac_mtz_filechooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{
  GtkWidget *mtz_fileselection;
  GtkWidget *mtz_label;
  GtkWidget *refmac_dialog;
  const gchar *filename;

  mtz_fileselection = lookup_widget(GTK_WIDGET(dialog),
				    "run_refmac_mtz_filechooserdialog");
  if (response_id == GTK_RESPONSE_OK) {
    save_directory_from_filechooser(mtz_fileselection);

    filename = gtk_file_chooser_get_filename 
      (GTK_FILE_CHOOSER(mtz_fileselection));
 
    /* save the label name */
    mtz_label = get_refmac_mtz_file_label();
    gtk_label_set_text(GTK_LABEL(mtz_label), filename);
    refmac_dialog = lookup_widget(mtz_label, "run_refmac_dialog");
    manage_refmac_column_selection(refmac_dialog);
  } 
  gtk_widget_destroy(mtz_fileselection);
}


void
on_run_refmac_mtz_filechooserdialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{

  GtkWidget *mtz_fileselection = lookup_widget(GTK_WIDGET(object),
					       "run_refmac_mtz_filechooserdialog");

  gtk_widget_destroy(mtz_fileselection);

}


void
on_run_refmac_mtz_filechooser_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *mtz_file_chooser;
  mtz_file_chooser = create_run_refmac_mtz_filechooserdialog();
  /* add file filter to filechooserbutton */
  add_filechooser_filter_button(mtz_file_chooser, COOT_DATASET_FILE_SELECTION);
  add_ccp4i_project_optionmenu(mtz_file_chooser, COOT_DATASET_FILE_SELECTION);
  gtk_widget_show(mtz_file_chooser);

}


void
on_run_refmac_file_help_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = create_run_refmac_file_help_dialog();
  gtk_widget_show(widget);


}

void
on_run_refmac_file_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "run_refmac_file_help_dialog");
  gtk_widget_destroy(widget);

}


void
on_run_refmac_sad_help_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = create_run_refmac_sad_help_dialog();
  gtk_widget_show(widget);

}


void
on_run_refmac_sad_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "run_refmac_sad_help_dialog");
  gtk_widget_destroy(widget);

}


void
on_run_refmac_nolabels_checkbutton_toggled (GtkToggleButton *togglebutton,
                                            gpointer         user_data)
{
  GtkWidget *labels = lookup_widget(GTK_WIDGET(togglebutton), "run_refmac_column_labels_frame");
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    gtk_widget_hide(labels);
  } else {
    gtk_widget_show(labels);
  }
}


/* we want to update the phases/hl boxes if the mtz changes */
#if 0				/* don't know what to do with optionmenu changed */
void
on_run_refmac_map_optionmenu_changed   (GtkOptionMenu   *optionmenu,
                                        gpointer         user_data)
{
  //update_refmac_column_labels_frame(optionmenu);

}
#endif

/* Actually, we only are interested in the state of this when the
   "Run" button is pressed */
void
on_run_refmac_diff_map_checkbutton_toggled (GtkToggleButton *togglebutton,
					    gpointer         user_data)
{
  
}


/* Actually, we only are interested in the state of this when the
   "Run" button is pressed */
void
on_run_refmac_ncs_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_refmac_use_ncs(1);
  } else {
    set_refmac_use_ncs(0);
  }

}

/* Actually, we only are interested in the state of this when the
   "Run" button is pressed */
/* not any more! We should show the phase column options depending on what exists */
void
on_run_refmac_phase_combine_checkbutton_toggled (GtkToggleButton *togglebutton,
						 gpointer         user_data)
{

#if 0				/* complete rewrite */
  GtkWidget *phases_hbox;
  GtkWidget *hl_hbox;
  GtkWidget *phases_optionmenu;
  GtkWidget *hl_optionmenu;
  GtkWidget *active_item;
  GtkWidget *phases_menu;
  GtkWidget *hl_menu;
  GtkWidget *phases_button;
  GtkWidget *hl_button;
  phases_hbox = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_phases_hbox");
  hl_hbox     = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_hl_hbox");
  phases_optionmenu = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_phases_optionmenu");
  hl_optionmenu     = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_hl_optionmenu");
  phases_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(phases_optionmenu));
  hl_menu     = gtk_option_menu_get_menu(GTK_OPTION_MENU(hl_optionmenu));
  phases_button = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_phases_radiobutton");
  hl_button     = lookup_widget(GTK_WIDGET(togglebutton), "refmac_dialog_hl_radiobutton");
  if (gtk_toggle_button_get_active(togglebutton)) {
    gtk_widget_show(phases_hbox);
    gtk_widget_show(hl_hbox);
    active_item = gtk_menu_get_active(GTK_MENU(phases_menu));
    if (active_item) {
      gtk_widget_set_sensitive(phases_hbox, TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(phases_button), TRUE);
    } else {
      gtk_widget_set_sensitive(phases_hbox, FALSE);
    }
    active_item = gtk_menu_get_active(GTK_MENU(hl_menu));
    if (active_item) {
      gtk_widget_set_sensitive(hl_hbox, TRUE);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(hl_button), TRUE);
    } else {
      gtk_widget_set_sensitive(hl_hbox, FALSE);
    }
    
  } else {
    gtk_widget_hide(phases_hbox);
    gtk_widget_hide(hl_hbox);
  }

#endif
}


void
on_baton_undo_button_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
   baton_build_delete_last_residue();
}


void
on_model_refine_dialog_undo_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{
   apply_undo();
}


void
on_undo_molecule_chooser_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{
   GtkWidget *widget = lookup_widget(GTK_WIDGET(button),
				     "undo_molecule_chooser_dialog");
   gtk_widget_destroy(widget);

}


void
on_model_refine_dialog_clear_pending_button_clicked (GtkButton *button,
						     gpointer user_data)
{
   clear_pending_picks();
}


void
on_skeleton_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "skeleton_dialog");
  GtkWidget *optionmenu = lookup_widget(window, "skeleton_map_optionmenu");
  int do_baton_mode = GPOINTER_TO_INT(user_data);

  skeletonize_map_by_optionmenu(optionmenu);
  gtk_widget_destroy(window);
  if (do_baton_mode)
    try_set_draw_baton(1);

}


void
on_skeleton_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				    "skeleton_dialog");
  gtk_widget_destroy(window);
}


void
on_skeleton1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   /* we need to fill the optionmenu.  Make it a graphics_info_t this time. */
  GtkWidget *w = wrapped_create_skeleton_dialog();
  gtk_widget_show(w);
}


void
on_virtual_trackball_menu_pops         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *flat_check_menu_item;
  GtkWidget *spherical_check_menu_item;

  flat_check_menu_item = lookup_widget(GTK_WIDGET(menuitem), 
				       "flat1");

  spherical_check_menu_item = lookup_widget(GTK_WIDGET(menuitem), 
				       "spherical_surface1");

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


void
on_model_refine_dialog_auto_fit_rotamer_togglebutton_toggled (GtkButton *button,
							gpointer user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button)))
     setup_auto_fit_rotamer(1);
  else 
     setup_auto_fit_rotamer(0);
}




void
on_residue_info2_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_residue_info_dialog();
}


void
on_import_all_dictionary_cifs1_activate (GtkMenuItem     *menuitem,
					 gpointer         user_data)
{
  import_all_refmac_cifs();
}


void
on_model_refine_dialog_destroy         (GtkWidget       *object,
                                        gpointer         user_data)
{
   unset_model_fit_refine_dialog();
}


void
on_residue_info_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   residue_info_apply_all_checkbutton_toggled();
}


void
on_residue_parameters1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_residue_info_dialog();
}



void
on_crosshairs1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget = create_crosshairs_dialog();
  GtkWidget *button = lookup_widget(widget, "crosshairs_on_radiobutton");

  if (draw_crosshairs_state()) { 
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
  } 
  gtk_widget_show(widget); 
}


void
on_crosshairs_on_radiobutton_toggled   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
    set_draw_crosshairs(1);
  } 
}


void
on_crosshairs_off_radiobutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
    set_draw_crosshairs(0);
  } 
}


void
on_display_crosshairs_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "crosshairs_dialog");

  gtk_widget_destroy(dialog);

}


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



void
on_add_alt_conf_ca_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
      set_add_alt_conf_split_type_number(0);
  }
}

void
on_add_alt_conf_whole_single_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
							  gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
      set_add_alt_conf_split_type_number(1);
  }
}

void
on_add_alt_conf_residue_range_radiobutton_toggled (GtkToggleButton *togglebutton,
						   gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
      set_add_alt_conf_split_type_number(2);
  }
}

void
on_add_alt_conf_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "add_alt_conf_dialog");
  unset_add_alt_conf_dialog();
  unset_add_alt_conf_define();
  gtk_widget_destroy(dialog);

}


void
on_validation_dialog_next_button_clicked (GtkButton       *button,
					  gpointer         user_data) {

}


void
on_validation_dialog_prev_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{

}


void
on_validation_dialog_cancel_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{

  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "");

  gtk_widget_destroy(widget); 
  
}


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




void
on_model_refine_dialog_redo_button_clicked (GtkButton       *button,
					    gpointer         user_data) {
  apply_redo();
}


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





gboolean
on_window1_configure_event             (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  gint upositionx, upositiony;
  GtkWindow *window = gtk_widget_get_window(widget);
  gtk_window_get_position(window, &upositionx, &upositiony);	    
  store_graphics_window_position(upositionx, upositiony);
  return FALSE;
}


void
on_run_state_file_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "run_state_file_dialog");
  gtk_widget_destroy(dialog);
  run_state_file();
}


void
on_run_state_file_cancel_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "run_state_file_dialog");
  gtk_widget_destroy(dialog);
}


void
on_edit_backbone_torsions_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
  clear_moving_atoms_object();
  /* FIXME: also clear out the edib backbone ramaplot, if it exists. */
/*   destroy_edit_backbone_rama_plot(); */
}


void
on_edit_backbone_torsion_rotate_peptide_button_pressed
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  int ix, iy;
  GdkModifierType state;
  GtkWindow *window = gtk_widget_get_window(button);
  gdk_window_get_pointer(window, &ix, &iy, &state);
  set_backbone_torsion_peptide_button_start_pos(ix, iy);
}


void
on_edit_backbone_torsion_rotate_peptide_button_released
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  /* ignored */
}


gboolean
on_edit_backbone_torsion_rotate_peptide_button_motion_notify_event
                                        (GtkWidget       *widget,
                                        GdkEventMotion  *event,
                                        gpointer         user_data)
{
  int ix, iy;
  GdkModifierType state;
  GtkWindow *window = gtk_widget_get_window(widget);
  gdk_window_get_pointer(window, &ix, &iy, &state);
  change_peptide_peptide_by_current_button_pos(ix, iy);
  return FALSE;
}


void
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_pressed
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  int ix, iy;
  GdkModifierType state;
  GtkWindow *window = gtk_widget_get_window(button);
  gdk_window_get_pointer(window, &ix, &iy, &state);
  set_backbone_torsion_carbonyl_button_start_pos(ix, iy);

}


void
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_released
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  /* ignored */
}


gboolean
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_motion_notify_event
                                        (GtkWidget       *widget,
                                        GdkEventMotion  *event,
                                        gpointer         user_data)
{

  int ix, iy;
  GdkModifierType state;
  GtkWindow *window = gtk_widget_get_window(widget);
  gdk_window_get_pointer(window, &ix, &iy, &state);
/*   printf("button moved to %d %d \n", ix, iy); */

  change_peptide_carbonyl_by_current_button_pos(ix, iy);
  return FALSE;
}


void
on_model_refine_dialog_edit_backbone_torsions_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) { 
    setup_backbone_torsion_edit(1);
  } else {
    setup_backbone_torsion_edit(0);
  }
}

void
on_edit_backbone_torsion_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				    "edit_backbone_torsions_dialog");
  accept_regularizement();	/* does a clear too. */
  destroy_edit_backbone_rama_plot();
  gtk_widget_destroy(window);
}


void
on_edit_backbone_torsion_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(button), 
				    "edit_backbone_torsions_dialog");
/*   clear_moving_atoms_object(); done as part of window destroy
     callback */
  destroy_edit_backbone_rama_plot();
  gtk_widget_destroy(window);
}


void
on_sequence_view1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   add_on_sequence_view_choices();
}


void
on_clear_simple_distances2_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   clear_simple_distances();
}


void
on_workflow_radiobutton_coords_only_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_workflow_radiobutton_coords_and_dataset_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_workflow_radiobutton_coords_and_map_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_workflow_radiobutton_map_from_dataset_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_workstate_radiobutton_clean_slate_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_workflow_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{

}

void
on_workflow_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{

}


void
on_select_map_for_fitting_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
//  GtkWidget *widget = lookup_widget(GTK_WIDGET(button),
//				    "select_fitting_map_dialog");

 // gtk_widget_destroy(widget);

}


void
on_model_refine_dialog_map_select_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   show_select_map_dialog();
}


void
on_clear_atom_labels1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  remove_all_atom_labels();
}



void
on_model_refine_dialog_add_alt_conf_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  altconf();
}

void
on_add_alt_conf_dialog_destroy         (GtkWidget       *object,
                                        gpointer         user_data)
{

  unset_add_alt_conf_define();
  unset_add_alt_conf_dialog();
}


void
on_run_refmac_help_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  /* GtkWidget *widget = create_run_refmac_help_dialog();
     we use the 'new' help now, but preserve the old dialog
     as backup for now */
  GtkWidget *widget = create_run_refmac_nolabels_help_dialog();
  gtk_widget_show(widget);
}


void
on_run_refmac_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "run_refmac_help_dialog");
  gtk_widget_destroy(widget);

}


void
on_run_refmac_nolabels_help_button_clicked      
					(GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = create_run_refmac_nolabels_help_dialog();
  gtk_widget_show(widget);
}


void
on_run_refmac_nolabels_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "run_refmac_nolabels_help_dialog");
  gtk_widget_destroy(widget);

}


void
on_no_restraints_info_dialog_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = lookup_widget(GTK_WIDGET(button), "no_restraints_info_dialog");
  gtk_widget_destroy(window);

}


void
on_no_cif_dictionary_bonds_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = lookup_widget(GTK_WIDGET(button), "no_cif_dictionary_bonds_dialog");
  gtk_widget_destroy(window);

}


void
on_pointer_atom_type_other_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *entry = lookup_widget(GTK_WIDGET(button), "pointer_atom_type_other_entry");
  gtk_widget_set_sensitive(entry, TRUE);

}


void
on_ligand_big_blob_dismiss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				   "ligand_big_blob_dialog");
  gtk_widget_destroy(window);

}


void
on_edit_chi_angles_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "edit_chi_angles_dialog");
  accept_regularizement();
  unset_moving_atom_move_chis();
  store_window_position(COOT_EDIT_CHI_DIALOG, widget);
  gtk_widget_destroy(widget);

}


void
on_edit_chi_angles_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "edit_chi_angles_dialog");
  clear_up_moving_atoms();	/* and remove the graphics object */
  unset_moving_atom_move_chis();
  store_window_position(COOT_EDIT_CHI_DIALOG, widget);
  gtk_widget_destroy(widget);

}


void
on_check_waters_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "check_waters_dialog");
  do_check_waters_by_widget(widget);
  gtk_widget_destroy(widget);
}


void
on_check_waters_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "check_waters_dialog");
  gtk_widget_destroy(widget);

}


void
on_edit_chi_angles_normal_rotation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  set_graphics_edit_current_chi(0);

}


void
on_geometry_distance_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    do_distance_define();

}




void
on_geometry_clear_last_distance_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  clear_last_simple_distance();

}


void
on_geometry_clear_all_distances_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  clear_simple_distances();
}


void
on_geometry_clear_atom_labels_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  remove_all_atom_labels();
}


void
on_geometry_dialog_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "geometry_dialog");
  /* should we clear geometry on close dialog?  Currently, I think not. */
  /* it is the COOT_DISTANCES_ANGLES_WINDOW, hmm. */
  store_window_position(COOT_DISTANCES_ANGLES_WINDOW, dialog);
  store_geometry_dialog(NULL);
  gtk_widget_destroy(dialog);

}


void
on_geometry_angle_togglebutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    do_angle_define();

}


void
on_distances_and_angles1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *widget = wrapped_create_geometry_dialog();
  
  set_transient_and_position(COOT_DISTANCES_ANGLES_WINDOW, widget);

  store_geometry_dialog(widget); /* needed to deactivate the distance
				    togglebutton after 2nd atoms
				    clicked in graphics */
  set_transient_and_position(COOT_UNDEFINED_WINDOW, widget);
  gtk_widget_show(widget);

}


void
on_geometry_dialog_destroy             (GtkWidget       *object,
                                        gpointer         user_data)
{

  /* OK, so the user moved the dialog somewhere, we want to store that
     position before we go. */

/* Nope!  we can't do the cast of object->widget, (then widget->window
   is NULL) and store function fails. */
/*   store_window_position(COOT_DISTANCES_ANGLES_WINDOW, GTK_WIDGET(object)); */
  
  /* However, we do want to unset the geometry_dialog pointer */
  store_geometry_dialog(NULL);

}


void
on_new_ligands_info_dialog_ok_button_clicked (GtkButton       *button,
					      gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "new_ligands_info_dialog");
  gtk_widget_destroy(w);
}


void on_no_new_ligands_info_dialog_ok_button_clicked (GtkButton       *button,
						      gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "no_new_ligands_info_dialog");
  gtk_widget_destroy(w);
}


void
on_zoom1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = create_zoom_dialog();
  set_zoom_adjustment(w);
  gtk_widget_show(w);

}


void
on_zoom_dialog_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "zoom_dialog");
  gtk_widget_destroy(w);

}


void
on_edit_chi_angles_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data)
{
   /* needs to set widget */
   /*  store_window_position(COOT_EDIT_CHI_DIALOG, widget); */
  unset_moving_atom_move_chis();
  set_show_chi_angle_bond(0);
}


void
on_check_waters_low_occ_dist_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_check_waters_zero_occ_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_get_monomer1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *dialog = wrapped_create_libcheck_monomer_dialog();
  
  gtk_widget_show(dialog);
}


void
on_libcheck_monomer_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "libcheck_monomer_dialog");
  GtkWidget *entry  = lookup_widget(widget, "libcheck_monomer_entry");
  handle_get_libcheck_monomer_code(entry);
  
/*   gtk_widget_destroy(widget);  done in handle_get_libcheck_monomer_code */
}


void on_libcheck_monomer_cancel_button_clicked(GtkButton       *button,
					       gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "libcheck_monomer_dialog");
  gtk_widget_destroy(widget);
}


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


void
on_recover_coordinates_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "recover_coordinates_dialog");

  execute_recover_session(widget); /* widget needed for lookup of user data */
  gtk_widget_destroy(widget);

}


void
on_recover_coordinates_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = lookup_widget(GTK_WIDGET(button), "recover_coordinates_dialog");
  gtk_widget_destroy(widget);
}


void
on_recover_session1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  recover_session();
}


void
on_centre_atom_label1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *widget = create_centre_atom_label_dialog();
  GtkWidget *on  = lookup_widget(widget, "centre_atom_label_radiobutton_on");
  GtkWidget *off = lookup_widget(widget, "centre_atom_label_radiobutton_off");
  int v = centre_atom_label_status();
  if (v) {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on), TRUE);
  } else { 
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off), TRUE);
  }
  gtk_widget_show(widget);
}


void
on_centre_atom_label_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "centre_atom_label_dialog");
  GtkWidget *on  = lookup_widget(dialog, "centre_atom_label_radiobutton_on");

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(on))) {
    set_label_on_recentre_flag(1);
  } else { 
    set_label_on_recentre_flag(0);
  }
  gtk_widget_destroy(dialog);
}


void
on_edit_chi_angles_help_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = create_chi_angle_help_dialog();
  gtk_widget_show(w);
}


void
on_help_chi_angles_dismiss_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "chi_angle_help_dialog");
  gtk_widget_destroy(w);
}


void
on_rotate_translate_obj_dialog_destroy (GtkWidget       *object,
                                        gpointer         user_data)
{
   /* need to save the position coordinates of dialog */
   rot_trans_reset_previous();
}


void
on_no_symmetry_warning_ok_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "no_symmetry_warning_dialog");
  gtk_widget_destroy(w);
}


void
on_nothing_to_recover_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "nothing_to_recover_dialog");
  gtk_widget_destroy(w);

}


void
on_superpose_dialog_superpose_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "superpose_dialog");
  execute_superpose(w);
  gtk_widget_destroy(w);

}

void
on_superpose_dialog_cancel_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "superpose_dialog");
  gtk_widget_destroy(w);
}


void
on_superpose_nonsense_ok_button_clicked (GtkButton       *button,
					 gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "superpose_nonsense_dialog");
  gtk_widget_destroy(w);

}

void
on_superpose_nonsense_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "superpose_dialog");
  gtk_widget_destroy(w);

}



void
on_ssm_superposition1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_superpose_dialog();

  /* we get returned w = 0 when there is no MMDBSSM. (We are doing it
     this way because we don't have to introduce HAVE_MMDBSSM into the
     *c* compiler arguments (this is simpler)).  */
  if (w) 
    gtk_widget_show(w);
  else 
    printf("Function not available - need to recompile with libmmdbssm\n");
}

void
on_add_terminal_residue_finds_none_ok_button_clicked (GtkButton       *button,
						      gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "add_terminal_residue_finds_none_dialog");
  gtk_widget_destroy(w);
}


void
on_single_map_sigma_checkbutton_toggled (GtkToggleButton *togglebutton,
					 gpointer         user_data)
{
  GtkWidget *window = lookup_widget(GTK_WIDGET(togglebutton),
				    "single_map_properties_dialog");
  GtkWidget *entry  = lookup_widget(window, "single_map_sigma_step_entry");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));
  const char *text = gtk_entry_get_text(GTK_ENTRY(entry));
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


void
on_no_bad_chiral_volumes_dialog_ok_button_clicked (GtkButton       *button,
						   gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "no_bad_chiral_volumes_dialog");
  gtk_widget_destroy(w);
}


void
on_check_chiral_volumes_ok_button_clicked (GtkButton       *button,
					   gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "check_chiral_volumes_dialog");
  check_chiral_volumes_from_widget(w);
  gtk_widget_destroy(w);

}


void
on_check_chiral_volumes_cancel_button_clicked (GtkButton       *button,
					       gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "check_chiral_volumes_dialog");
  gtk_widget_destroy(w);
}


void
on_find_bad_chiral_atoms1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = create_check_chiral_volumes_dialog();
  fill_chiral_volume_molecule_combobox(w);
  gtk_widget_show(w);
}



void
on_chiral_volume_baddies_dialog_cancel_button_clicked (GtkButton       *button,
						       gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "bad_chiral_volumes_dialog");
  gtk_widget_destroy(w);
}


void
on_rigid_body_refinement_failed_dialog_ok_button_clicked (GtkButton *button,
							  gpointer user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
			       "rigid_body_refinement_failed_dialog");
  gtk_widget_destroy(w);
}


void
on_baton_mode_calculate_skeleton_ok_button_clicked (GtkButton       *button,
						    gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
			       "baton_mode_make_skeleton_dialog");
  baton_mode_calculate_skeleton(w); /* get the imol from here */
  gtk_widget_destroy(w);
}


void
on_baton_mode_calculate_skeleton_cancel_button_clicked (GtkButton       *button,
							gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
			       "baton_mode_make_skeleton_dialog");
  gtk_widget_destroy(w);
}


void
on_column_label_expert_mode_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *frame = lookup_widget(GTK_WIDGET(button),
				   "column_labels_resolution_limit_frame");

/*   GtkWidget *f_optionmenu = lookup_widget(GTK_WIDGET(button), "optionmenu1"); */

  GtkWidget *combobox = lookup_widget(GTK_WIDGET(button), "column_selector_amplitudes_combobox");

  /* we also need to redo the F column label chooser to include anomalous option  */


   if (gtk_widget_get_visible(frame)) {
    gtk_widget_hide(frame);
   } else {
     gtk_widget_show(frame);
/*      fill_f_optionmenu_with_expert_options(f_optionmenu); */
     fill_combobox_with_expert_options(combobox);
   }

}


void
on_column_labels_use_resolution_limits_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *frame = lookup_widget(GTK_WIDGET(togglebutton),
				   "resolution_limits_hbox");
  if (gtk_toggle_button_get_active(togglebutton))
     gtk_widget_set_sensitive(frame, TRUE);
  else
     gtk_widget_set_sensitive(frame, FALSE);
}




void
on_merge_molecules_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "merge_molecules_dialog");
  do_merge_molecules(w);
  gtk_widget_destroy(w);

}


void
on_merge_molecules_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "merge_molecules_dialog");
  gtk_widget_destroy(w);

}


void
on_mutate_sequence_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "mutate_sequence_dialog");
  do_mutate_sequence(w);
  gtk_widget_destroy(w);

}


void
on_mutate_sequence_cancel_button_clicked (GtkButton       *button,
					  gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "mutate_sequence_dialog");
  gtk_widget_destroy(w);

}


void
on_merge_molecules1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   gtk_widget_show(w);
}


void
on_mutate_molecule1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_mutate_sequence_dialog();
   gtk_widget_show(w);
}



void
on_draw_hydrogens_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  // Applying the bond parameters applies the bond width too, which is
  // a bit of extra overhead.
  GtkWidget *w = lookup_widget(GTK_WIDGET(togglebutton), "bond_parameters_dialog");
  apply_bond_parameters(w);
}


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

  GtkWidget *w = lookup_widget(GTK_WIDGET(togglebutton), "bond_parameters_dialog");
  apply_bond_parameters(w);
  */
}

/* radio button is the N-terminal button */
void
on_renumber_residue_range_radiobutton_1_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  GtkWidget *entry_1 = lookup_widget(GTK_WIDGET(togglebutton), "renumber_residue_range_resno_1_entry");
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    gtk_widget_set_sensitive(GTK_WIDGET(entry_1), FALSE);
  } else {
    gtk_widget_set_sensitive(GTK_WIDGET(entry_1), TRUE);
  }
}

void
on_renumber_residue_range_radiobutton_3_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  GtkWidget *entry_2 = lookup_widget(GTK_WIDGET(togglebutton), "renumber_residue_range_resno_2_entry");
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    gtk_widget_set_sensitive(GTK_WIDGET(entry_2), TRUE);
  } else {
    gtk_widget_set_sensitive(GTK_WIDGET(entry_2), FALSE);
  }
}


void
on_renumber_residue_range_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
			       "renumber_residue_range_dialog");
  renumber_residues_from_widget(w);
  gtk_widget_destroy(w);

}


void
on_renumber_residue_range_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
			       "renumber_residue_range_dialog");
  gtk_widget_destroy(w);

}


void
on_add_OXT_c_terminus_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_add_OXT_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_add_OXT_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "add_OXT_dialog");
  apply_add_OXT_from_widget(w);
  gtk_widget_destroy(w);
}


void
on_add_OXT_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "add_OXT_dialog");
  gtk_widget_destroy(w);
}


void
on_model_refine_dialog_add_OXT_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_add_OXT_dialog();
  gtk_widget_show(w);
}


void
on_renumber_residues1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_renumber_residue_range_dialog();
  gtk_widget_show(w);
}


void
on_bond_parameters1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_bond_parameters_dialog();
  gtk_widget_show(w);
}

void
on_bond_parameters_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "bond_parameters_dialog");
  apply_bond_parameters(w);
  gtk_widget_destroy(w);
}

void
on_bond_parameters_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "bond_parameters_dialog");
  apply_bond_parameters(w);
}




void
on_bond_parameters_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "bond_parameters_dialog");
  gtk_widget_destroy(w);

}


void
on_background_colour1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *wb = lookup_widget(GTK_WIDGET(menuitem), "background_black1");
  GtkWidget *ww = lookup_widget(GTK_WIDGET(menuitem), "background_white1");

  if (background_is_black_p())
/*     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE); */
    				/* a GtkRadioMenuItem not a toggle button */
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(wb), TRUE);
  else 
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(ww), TRUE);

}



void
on_ligand_no_blobs_OK_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "ligand_no_blobs_dialog");
  gtk_widget_destroy(dialog);
}



void
on_new_delete_molecules_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "new_close_molecules_dialog");
  new_close_molecules(w);
/*   gtk_widget_destroy(w); */
}


void
on_new_delete_molecules_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "new_close_molecules_dialog");
  gtk_widget_destroy(w);

}


void
on_ramachandran_plot2_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  add_on_rama_choices("");

}


void
on_incorrect_chiral_volumes1_activate  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = create_check_chiral_volumes_dialog();
  fill_chiral_volume_molecule_combobox(w);
  gtk_widget_show(w);

}


void
on_unmodelled_blobs1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_unmodelled_blobs_dialog();
  gtk_widget_show(w);

}


void
on_find_blobs_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "unmodelled_blobs_dialog");
  execute_find_blobs_from_widget(w);
  gtk_widget_destroy(w);

}


void
on_find_blobs_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "unmodelled_blobs_dialog");
  gtk_widget_destroy(w);

}



void
on_chiral_restraints_problem_ok_button_clicked (GtkButton       *button,
						gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "chiral_restraints_problem_dialog");
  gtk_widget_destroy(w);
}

/* We'll keep this for now because it is used by
   create_check_waters_diff_map_dialog() and that is still in the
   interface definition (glade file).
 */
void
on_check_waters_diff_map_ok_button_clicked (GtkButton       *button,
					    gpointer         user_data)
{

/*   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "check_waters_diff_map_dialog"); */
/*   check_waters_by_difference_map_by_widget(w); */
/*   gtk_widget_destroy(w); */
  

}


void
on_check_waters_diff_map_cancel_button_clicked (GtkButton       *button,
						gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "check_waters_diff_map_dialog");
  gtk_widget_destroy(w);

}


/* void */
/* on_check_waters_by_difference_map_variance1_activate (GtkMenuItem     *menuitem, */
/* 						      gpointer         user_data) */
/* { */
/*   GtkWidget *w = wrapped_create_check_waters_diff_map_dialog(); */
/*   gtk_widget_show(w); */
/* } */


void
on_interesting_waters_by_difference_map_check_ok_button_clicked (GtkButton       *button,
								 gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
			       "interesting_waters_by_difference_map_check_dialog");
  gtk_widget_destroy(w);

}


void
on_nothing_bad_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "nothing_bad_dialog");
  gtk_widget_destroy(w);
}


void
on_skeletonize_map_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}


void
on_skeletonize_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}


void
on_antialiasing1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = create_antialiasing_dialog();
  GtkWidget *checkbutton;
  if (do_anti_aliasing_state()) { 
    checkbutton = lookup_widget(w, "antialias_dialog_yes_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
  } 
  gtk_widget_show(w);
}


void
on_antialias_dialog_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_do_anti_aliasing(1);
  else 
    set_do_anti_aliasing(0);

}


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


void
on_antialiasing_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "antialiasing_dialog");
  gtk_widget_destroy(w);

}


void
on_geometry_graphs_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "geometry_graphs_dialog");
  if (w)
    gtk_widget_destroy(w);
  else 
    printf("ERROR getting dialog in on_geometry_graphs_ok_button_clicked\n");

}


void
on_geometry_graphs_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data)
{

  /* This is not causing the GTK_IS_WIDGET (widget) unref failure */

   GtkWidget *w = lookup_widget(GTK_WIDGET(object), "geometry_graphs_dialog"); 
   if (! w) { 
     printf("ERROR getting dialog in on_geometry_graphs_dialog_destroy\n"); 
   } else {  
/*      printf("DEBUG:: unsetting and freeing graph:\n");  */
     unset_geometry_graph(w); 
     free_geometry_graph(w); 
   }

}


void
on_save_symmetry_coords_fileselection_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "save_symmetry_coords_fileselection");
  save_symmetry_coords_from_fileselection(w);
  gtk_widget_destroy(w);
}


void
on_save_symmetry_coords_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "save_symmetry_coords_fileselection");
  gtk_widget_destroy(w);

}


void
on_save_symmetry_coordinates1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  setup_save_symmetry_coords();
}


void
on_experimental1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_geometry_analysis1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "geometry";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "geometry_analysis1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_peptide_omega_analysis1_activate\n");
  }

}


void
on_peptide_omega_analysis1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "omega";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "peptide_omega_analysis1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_peptide_omega_analysis1_activate\n");
  }

}

void
on_ncs_differences1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "ncs-diffs";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "ncs_differences1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_peptide_omega_analysis1_activate\n");
  }
}

////B B FACTOR
void
on_temp_fact_analysis1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "calc b factor";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "temp_fact_analysis1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_temp_fact_analysis1_activate\n");
  }

}
////E B FACTOR

void
on_temp_fact_variance_analysis1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "b factor";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "temp_fact_variance_analysis1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_temp_fact_variance_analysis1_activate\n");
  }

}


void
on_rotamer_analysis1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "rotamer";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "rotamer_analysis1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_rotamer_analysis1_activate\n");
  }
}


void
on_density_fit_analysis1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "density-fit";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "density_fit_analysis1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_density_fit1_activate\n");
  }

}


void
on_stereo1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = create_stereo_dialog();
  GtkWidget *checkbutton;

  if (stereo_mode_state() == 1) { /* coot::HARDWARE_STEREO_MODE */
    checkbutton = lookup_widget(w, "stereo_dialog_hardware_stereo_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
  }
  if (stereo_mode_state() == 2) { /* coot::SIDE_BY_SIDE_STEREO */
    checkbutton = lookup_widget(w, "stereo_dialog_side_by_side_stereo_crosseyed_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
  }
  if (stereo_mode_state() == 4) { /* coot::SIDE_BY_SIDE_STEREO_WALL_EYE */
    checkbutton = lookup_widget(w, "stereo_dialog_side_by_side_stereo_walleyed_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
  }
  if (stereo_mode_state() == 3) { /* coot::DTI_SIDE_BY_SIDE_STEREO */
    checkbutton = lookup_widget(w, "stereo_dialog_dti_side_by_side_stereo_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
  }
  if (stereo_mode_state() == 5) { /* coot::ZALMAN_STEREO */
    checkbutton = lookup_widget(w, "stereo_dialog_zalman_stereo_radiobutton");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
  }
  gtk_widget_show(w);
}



void
on_stereo_dialog_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "stereo_dialog");
  gtk_widget_destroy(w);
}


void
on_stereo_dialog_mono_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    mono_mode();
}


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
      mono_togglebutton = lookup_widget(GTK_WIDGET(togglebutton), "stereo_dialog_mono_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(mono_togglebutton), TRUE);
      /* do this in C rather than mess about calling c++ function: */
      nothing_bad_dialog = create_nothing_bad_dialog();
      label_widget = lookup_widget(nothing_bad_dialog, "nothing_bad_label");
      gtk_label_set_text(GTK_LABEL(label_widget), "This computer appears not to be able\nto do hardware stereo");
      gtk_widget_show(nothing_bad_dialog);
    }
  }
}

void
on_stereo_dialog_side_by_side_stereo_crosseyed_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    side_by_side_stereo_mode(0); /* passed used_wall_eye flag */
}


void
on_stereo_dialog_side_by_side_stereo_walleyed_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    side_by_side_stereo_mode(1); /* passed used_wall_eye flag */
  }
}

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
void
on_preferences1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  preferences();
}


void
on_preferences_general_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_GENERAL_PREFERENCES);
}


void
on_preferences_bond_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_BOND_PREFERENCES);
}


void
on_preferences_map_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_MAP_PREFERENCES);
}



void
on_preferences_geometry_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_GEOMETRY_PREFERENCES);
}




void
on_preferences_colour_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_COLOUR_PREFERENCES);
}




void
on_preferences_other_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  show_hide_preferences_tabs(toggletoolbutton, COOT_OTHER_PREFERENCES);
}


void
on_preferences_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "preferences");
  save_preferences();
  gtk_widget_destroy(w);
  clear_preferences();
}

void
on_preferences_reset_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  reset_preferences();

}

void
on_preferences_destroy                 (GtkWidget       *object,
                                        gpointer         user_data)
{
  clear_preferences();
}

void
on_preferences_geometry_cis_peptide_bad_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MARK_CIS_BAD, 1);
    set_mark_cis_peptides_as_bad(1);
  }
}

void
on_preferences_geometry_cis_peptide_bad_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MARK_CIS_BAD, 0);
    set_mark_cis_peptides_as_bad(0);
  }
}

void
on_preferences_bond_colours_hscale_value_changed
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

void
on_preferences_bond_colours_checkbutton_toggled
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


void
on_preferences_bg_colour_black_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, 0, 0, 0);
    set_background_colour(0, 0, 0);
  }
}


void
on_preferences_bg_colour_white_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_float3(PREFERENCES_BG_COLOUR, 1, 1, 1);
    set_background_colour(1, 1, 1);
  }

}


void
on_preferences_bg_colour_own_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
#if (GTK_MAJOR_VERSION > 1)
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
#endif

}


void
on_preferences_bg_colour_colorbutton_color_set
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


void
on_preferences_bg_colour_colorbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w;
  w = lookup_widget(GTK_WIDGET(button), "preferences_bg_colour_own_radiobutton");
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);

}



void
on_preferences_map_radius_entry_activate
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


void
on_preferences_map_radius_entry_changed
                                        (GtkEditable     *editable,
					 gpointer         user_data)
{
  GtkEntry *entry;
  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(editable), "preferences_map_radius_entry"));
  const gchar *text = gtk_entry_get_text(entry);
  float fval = 0;
  fval = atof(text);
  if ((fval > 0) && (fval <1000)) {
    preferences_internal_change_value_float(PREFERENCES_MAP_RADIUS, fval);
    set_map_radius(fval);
  }
  
}


void
on_preferences_map_increment_size_entry_activate
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


void
on_preferences_map_increment_size_entry_changed
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


void
on_preferences_map_diff_increment_entry_activate
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


void
on_preferences_map_diff_increment_entry_changed
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


void
on_preferences_map_sampling_entry_activate
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


void
on_preferences_map_sampling_entry_changed
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


void
on_preferences_map_dynamic_sampling_checkbutton_toggled
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


void
on_preferences_map_dynamic_size_checkbutton_toggled
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


void
on_preferences_diff_map_colours_coot_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SWAP_DIFF_MAP_COLOURS, 0);
    set_swap_difference_map_colours(0);
  }

}


void
on_preferences_diff_map_colours_o_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SWAP_DIFF_MAP_COLOURS, 1);
    set_swap_difference_map_colours(1);
  }

}


void
on_preferences_map_colours_hscale_value_changed
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


void
on_preferences_smooth_scroll_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SMOOTH_SCROLL, 1);
    set_smooth_scroll_flag(1);
  }

}


void
on_preferences_smooth_scroll_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_SMOOTH_SCROLL, 0);
    set_smooth_scroll_flag(0);
  }

}


void
on_preferences_smooth_scroll_steps_entry_activate
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


void
on_preferences_smooth_scroll_steps_entry_changed
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


void
on_preferences_smooth_scroll_limit_entry_activate
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


void
on_preferences_smooth_scroll_limit_entry_changed
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


void
on_preferences_map_drag_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAP_DRAG, 1);
    set_active_map_drag_flag(1);
  }

}


void
on_preferences_map_drag_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAP_DRAG, 0);
    set_active_map_drag_flag(0);
  }

}


void
on_preferences_antialias_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ANTIALIAS, 1);
    set_do_anti_aliasing(1);
  }

}


void
on_preferences_antialias_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ANTIALIAS, 0);
    set_do_anti_aliasing(0);
  }

}


void
on_preferences_hid_spherical_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_VT_SURFACE, 2);
    vt_surface(2);
  }

}


void
on_preferences_hid_flat_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_VT_SURFACE, 1);
    vt_surface(1);
  }

}


void
on_preferences_filechooser_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_CHOOSER, 0);
    set_file_chooser_selector(0);
  }

}


void
on_preferences_filechooser_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_CHOOSER, 1);
    set_file_chooser_selector(1);
  }

}


void
on_preferences_file_overwrite_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_OVERWRITE, 1);
    set_file_chooser_overwrite(1);
  }

}


void
on_preferences_file_overwrite_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_OVERWRITE, 0);
    set_file_chooser_overwrite(0);
  }

}

void
on_preferences_file_filter_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_FILTER, 1);
    set_filter_fileselection_filenames(1);
  }

}


void
on_preferences_file_filter_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_FILTER, 0);
    set_filter_fileselection_filenames(0);
  }

}

void
on_preferences_file_sort_by_date_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_SORT_DATE, 1);
    set_sticky_sort_by_date();
  }

}


void
on_preferences_file_sort_by_date_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FILE_SORT_DATE, 0);
    unset_sticky_sort_by_date();
  }

}

void
on_preferences_dialog_accept_docked_radiobutton_toggled
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


void
on_preferences_dialog_accept_detouched_radiobutton_toggled
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

void
on_preferences_dialog_accept_docked_show_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ACCEPT_DIALOG_DOCKED_SHOW, 1);
    set_accept_reject_dialog_docked_show(1);
  }

}


void
on_preferences_dialog_accept_docked_hide_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_ACCEPT_DIALOG_DOCKED_SHOW, 0);
    set_accept_reject_dialog_docked_show(0);
  }

}


void
on_preferences_dialog_accept_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_IMMEDIATE_REPLACEMENT, 0);
    set_refinement_immediate_replacement(0);
  }

}


void
on_preferences_dialog_accept_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_IMMEDIATE_REPLACEMENT, 1);
    set_refinement_immediate_replacement(1);
  }

}


void
on_preferences_recentre_pdb_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_RECENTRE_PDB, 1);
    set_recentre_on_read_pdb(1);
  }

}


void
on_preferences_recentre_pdb_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_RECENTRE_PDB, 0);
    set_recentre_on_read_pdb(0);
  }

}


void
on_preferences_console_info_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_CONSOLE_COMMANDS, 1);
    set_console_display_commands_state(1);
  }

}


void
on_preferences_console_info_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_CONSOLE_COMMANDS, 0);
    set_console_display_commands_state(0);
  }

}


void
on_preferences_tips_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_TIPS, 1);
    set_tip_of_the_day_flag(1);
  }

}


void
on_preferences_tips_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_TIPS, 0);
    set_tip_of_the_day_flag(0);
  }

}


void
on_preferences_refinement_speed_molasses_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, 4);
    set_dragged_refinement_steps_per_frame(4);
  }

}


void
on_preferences_refinement_speed_crock_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, 120);
    set_dragged_refinement_steps_per_frame(120);
  }

}


void
on_preferences_refinement_speed_default_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_REFINEMENT_SPEED, 80);
    set_dragged_refinement_steps_per_frame(80);
  }

}


void
on_preferences_refinement_speed_own_radiobutton_toggled
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


void
on_preferences_refinement_speed_entry_activate
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


void
on_preferences_refinement_speed_entry_changed
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


void
on_preferences_spin_speed_entry_activate
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


void
on_preferences_spin_speed_entry_changed
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


void
on_preferences_font_size_small_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, 1);
    set_font_size(1);
  }

}


void
on_preferences_font_size_medium_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, 2);
    set_font_size(2);
  }

}


void
on_preferences_font_size_large_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_FONT_SIZE, 3);
    set_font_size(3);
  }

}

void
on_preferences_font_size_others_radiobutton_toggled
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

void
on_preferences_font_size_combobox_changed
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

void
on_preferences_font_colour_default_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
      preferences_internal_change_value_float3(PREFERENCES_FONT_COLOUR, 1.0, 0.8, 0.8);
      set_font_colour(1.0, 0.8, 0.8);
   }
}


void
on_preferences_font_colour_own_radiobutton_toggled
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



void
on_preferences_font_colorbutton_color_set
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

void
on_preferences_font_colorbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  /* actually not doing anything */
}


void
on_preferences_pink_pointer_entry_activate
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


void
on_preferences_pink_pointer_entry_changed
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


void
on_preferences_bond_width_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
  gint val;
  val = gtk_combo_box_get_active(combobox);
  val += 1;  /* offset */
  preferences_internal_change_value_int(PREFERENCES_BONDS_THICKNESS, val);
  set_default_bond_thickness(val);
}


void
on_preferences_model_toolbar_show_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    show_modelling_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_SHOW, 1);
  }
 
}


void
on_preferences_model_toolbar_right_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(0);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 0);
  }

}


void
on_preferences_model_toolbar_left_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(1);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 1);
  }

}


void
on_preferences_model_toolbar_top_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(2);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 2);
  }
}


void
on_preferences_model_toolbar_bottom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_model_toolbar_docked_position(3);
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_POSITION, 3);
  }
}


void
on_preferences_model_toolbar_hide_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    hide_modelling_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_SHOW, 0);
  }

}


void
on_preferences_model_toolbar_main_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_preferences_model_toolbar_all_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_preferences_model_toolbar_style_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_STYLE, 1);
    set_model_toolbar_style(1);
  }

}


void
on_preferences_model_toolbar_style_both_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_STYLE, 2);
    set_model_toolbar_style(2);
  }

}


void
on_preferences_model_toolbar_style_text_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MODEL_TOOLBAR_STYLE, 3);
    set_model_toolbar_style(3);
  }

}


void
on_preferences_model_toolbar_show_icon_all_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  show_model_toolbar_all_icons();
  
}


void
on_preferences_model_toolbar_show_icon_selection_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  show_model_toolbar_main_icons();

}

void
on_preferences_main_toolbar_show_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    show_main_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_SHOW, 1);
  }

}


void
on_preferences_main_toolbar_hide_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    hide_main_toolbar();
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_SHOW, 0);
  }

}


/* not used currently */
void
on_preferences_main_toolbar_top_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_preferences_main_toolbar_bottom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_preferences_main_toolbar_right_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_preferences_main_toolbar_left_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_preferences_main_toolbar_style_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_STYLE, 1);
    set_main_toolbar_style(1);
  }

}


void
on_preferences_main_toolbar_style_both_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_STYLE, 2);
    set_main_toolbar_style(2);
  }

}


void
on_preferences_main_toolbar_style_text_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) {
    preferences_internal_change_value_int(PREFERENCES_MAIN_TOOLBAR_STYLE, 3);
    set_main_toolbar_style(3);
  }

}


/* end preferences */

void
on_diff_map_peaks_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "diff_map_peaks_dialog");
  clear_diff_map_peaks();
  gtk_widget_destroy(w);
}


void
on_generate_diff_map_peaks_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "generate_diff_map_peaks_dialog");
  difference_map_peaks_by_widget(w);
  gtk_widget_destroy(w);

}


void
on_generate_diff_map_peaks_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "generate_diff_map_peaks_dialog");
  gtk_widget_destroy(w);

}


void
on_difference_map_peaks1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_generate_diff_map_peaks_dialog();
  gtk_widget_show(w);
}


void
on_superpose_reference_chain_checkbutton_toggled(GtkToggleButton *togglebutton,
						 gpointer         user_data) {

  GtkWidget *combobox = lookup_widget(GTK_WIDGET(togglebutton), 
					"superpose_dialog_reference_chain_combobox");
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


void
on_superpose_moving_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *combobox = lookup_widget(GTK_WIDGET(togglebutton), 
				      "superpose_dialog_moving_chain_combobox");

  if (gtk_toggle_button_get_active(togglebutton)) {
    fill_superpose_combobox_with_chain_options(combobox, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
  } else {
    gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
  }
}


void
on_draw_ncs_ghosts_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
/* Function no longer used.  Kept in glade (not visible) for historical reasons

   GtkWidget *w = lookup_widget(GTK_WIDGET(togglebutton), "bond_parameters_dialog");
   if (gtk_toggle_button_get_active(togglebutton)) { 
      printf("yes radiobutton toggled on.\n");
      make_ncs_ghosts_maybe(w);
   }
*/
}


void
on_draw_ncs_ghosts_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_ncs_maps_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "ncs_maps_dialog");
   make_dynamically_transformed_ncs_maps_by_widget(w);
   gtk_widget_destroy(w);
}


void
on_ncs_maps_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "ncs_maps_dialog");
   gtk_widget_destroy(w);
}


void
on_ncs_maps1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_ncs_maps_dialog();
   gtk_widget_show(w);
}


void
on_rotamer_selection_dialog_destroy    (GtkWidget       *object,
                                        gpointer         user_data)
{
   /* set the dialog
      store_window_position(COOT_ROTAMER_SELECTION_DIALOG, dialog); */
   set_graphics_rotamer_dialog(NULL);
}


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


void
on_pointer_distances_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "pointer_distances_dialog");
   execute_pointer_distances_settings(dialog);
   gtk_widget_destroy(dialog);
}


void
on_pointer_distances1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = create_pointer_distances_dialog();
   fill_pointer_distances_widget(w);
   gtk_widget_show(w);
}


void
on_align_and_mutate_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "align_and_mutate_dialog");
   int handled_state = do_align_mutate_sequence(dialog);
   if (handled_state == 1) 
      gtk_widget_destroy(dialog);

}


void
on_align_and_mutate_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "align_and_mutate_dialog");
   gtk_widget_destroy(dialog);

}


void
on_ramachandran_plot_differences_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "ramachandran_plot_differences_dialog");
   int istat = do_ramachandran_plot_differences_by_widget(w);
   if (istat) 			/* the plot was drawn (i.e. no chain selection funnies) */
      gtk_widget_destroy(w);

}


void
on_ramachandran_plot_differences_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "ramachandran_plot_differences_dialog");
   gtk_widget_destroy(w);

}


void
on_ramachandran_differences_plot1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_ramachandran_plot_differences_dialog();
   gtk_widget_show(w);

}


void
on_ramachandran_plot_differences_first_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *combobox = lookup_widget(GTK_WIDGET(togglebutton), 
				       "ramachandran_plot_differences_first_chain_combobox");
   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
      fill_ramachandran_plot_differences_combobox_with_chain_options(combobox, 1);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
   }

}


void
on_ramachandran_plot_differences_second_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *combobox = lookup_widget(GTK_WIDGET(togglebutton), 
				       "ramachandran_plot_differences_second_chain_combobox");

   if (gtk_toggle_button_get_active(togglebutton)) {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), TRUE);
      fill_ramachandran_plot_differences_combobox_with_chain_options(combobox, 0);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(combobox), FALSE);
   }

}


void
on_check_waters1_activate              (GtkMenuItem     *menuitem,
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



void
on_checked_waters_baddies_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "checked_waters_baddies_dialog");
   gtk_widget_destroy(w);

}


void
on_align_and_mutate1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_align_and_mutate_dialog();
   gtk_widget_show(w);

}


void
on_about_other_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{

}


void
on_delete_item_keep_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
}


void
on_delete_item_dialog_destroy          (GtkWidget       *object,
                                        gpointer         user_data)
{

   clear_pending_delete_item();
   clear_pending_picks(); 
   normal_cursor();
   store_delete_item_widget(NULL);
}


void
on_save_state1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
/*    save_state(); old inteface - before DK. */

   GtkWidget *fileselection = coot_save_state_chooser();
   /*   gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),  
	save_state_file_name_raw()); */
   set_filename_for_filechooserselection(fileselection,
					 save_state_file_name_raw());

   add_filename_filter_button(fileselection, COOT_SCRIPTS_FILE_SELECTION);
   add_ccp4i_project_optionmenu(fileselection, 
                                COOT_SCRIPTS_FILE_SELECTION);
   set_file_selection_dialog_size(fileselection);
   gtk_widget_show(fileselection);
}


void
on_geometry_torsion_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    do_torsion_define();
}


void
on_main_window_statusbar_text_popped   (GtkStatusbar    *statusbar,
                                        guint            context_id,
                                        gchar           *text,
                                        gpointer         user_data)
{

}


void
on_main_window_statusbar_text_pushed   (GtkStatusbar    *statusbar,
                                        guint            context_id,
                                        gchar           *text,
                                        gpointer         user_data)
{

}


void
on_fit_loop1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

/*    GtkWidget *w = wrapped_fit_loop_dialog(); */
/*    gtk_widget_show(w); */
}


void
on_fit_loop_by_rama_search1_activate (GtkMenuItem     *menuitem,
				      gpointer         user_data)
{

   GtkWidget *w = wrapped_fit_loop_rama_search_dialog();
   gtk_widget_show(w);
}


void
on_fit_loop_by_database_search1_activate (GtkMenuItem     *menuitem,
					  gpointer         user_data)
{
  wrapped_fit_loop_db_loop_dialog();
}


void
on_fit_loop_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "mutate_sequence_dialog");
   fit_loop_from_widget(w);
   gtk_widget_destroy(w);

}


void
on_reset_view1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   reset_view();
}


void
on_base_chooser_A_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"nucleic_acid_base_chooser_dialog");
   gtk_widget_destroy(w);
   do_base_mutation("A");
}


void
on_base_chooser_C_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"nucleic_acid_base_chooser_dialog");
   gtk_widget_destroy(w);
   do_base_mutation("C");
}


void
on_base_chooser_G_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"nucleic_acid_base_chooser_dialog");
   gtk_widget_destroy(w);
   do_base_mutation("G");
}


void
on_base_chooser_T_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"nucleic_acid_base_chooser_dialog");
   gtk_widget_destroy(w);
   do_base_mutation("T");
}


void
on_base_chooser_U_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"nucleic_acid_base_chooser_dialog");
   gtk_widget_destroy(w);
   do_base_mutation("U");
}


void
on_base_chooser_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"nucleic_acid_base_chooser_dialog");
   gtk_widget_destroy(w);
}


void
on_nucleic_acid_base_chooser_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
   clear_pending_picks();
}


void
on_change_chain_ids2_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_change_chain_id_dialog();
   gtk_widget_show(w);

}


void
on_change_chains_rechain_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "change_chain_id_dialog"); 
   change_chain_id_by_widget(w);
   gtk_widget_destroy(w);
}


void
on_change_chain_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "change_chain_id_dialog"); 
   gtk_widget_destroy(w);

}


void
on_delete_item_residue_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_delete_residue_zone_mode();

}


void
on_on_line_docs_url1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = create_doc_urls_dialog();
   gtk_widget_show(w);

}


void
on_on_line_documentation_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "doc_urls_dialog");
   gtk_widget_destroy(w);

}




void
on_save_state_cancel_button1_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button),
				"save_state_fileselection");
   gtk_widget_destroy(w);
}


void
on_change_chain_residue_range_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{


}


void
on_change_chain_residue_range_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton), 
				   "change_chain_id_residue_range_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else 
      gtk_widget_set_sensitive(hbox, FALSE);
}


void
on_mutate_sequence_do_autofit_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   int imol_map = -1;

   if (gtk_toggle_button_get_active(togglebutton)) { 
      imol_map = imol_refinement_map();
      if (imol_map == -1) { 
	 gtk_toggle_button_set_active(togglebutton, FALSE);
	 show_select_map_dialog();
	 info_dialog("A map has not yet been assigned for Refinement/Fitting");
      }
   }
}

void
on_mutate_sequence_use_ramachandran_restraints_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   /* not doing anything?! */
}

void
on_check_waters_b_factor_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton), 
				   "check_waters_b_factor_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);
}




void
on_check_waters_min_dist_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton), 
				   "check_waters_min_dist_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}


void
on_check_waters_max_dist_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton), 
				   "check_waters_max_dist_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}


void
on_check_waters_map_sigma_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton), 
				   "check_waters_sigma_level_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);
}



void
on_check_waters_by_difference_map_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(togglebutton), 
				   "check_waters_by_difference_map_hbox");
   if (gtk_toggle_button_get_active(togglebutton))
      gtk_widget_set_sensitive(hbox, TRUE);
   else
      gtk_widget_set_sensitive(hbox, FALSE);

}




void
on_residue_info_occ_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *entry = lookup_widget(GTK_WIDGET(togglebutton), 
				   "residue_info_master_atom_occ_entry");
  GtkWidget *alt_conf_checkbutton = lookup_widget(GTK_WIDGET(togglebutton),
						  "residue_info_occ_apply_to_altconf_checkbutton");
  
  if (gtk_toggle_button_get_active(togglebutton)) { 
    gtk_widget_set_sensitive(entry, TRUE);
  } else { 
    if (! gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(alt_conf_checkbutton)))
      gtk_widget_set_sensitive(entry, FALSE);
  } 
}


void
on_residue_info_occ_apply_to_altconf_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  GtkWidget *occ_entry = lookup_widget(GTK_WIDGET(togglebutton), 
				       "residue_info_master_atom_occ_entry");
  GtkWidget *alt_conf_entry = lookup_widget(GTK_WIDGET(togglebutton), 
					    "residue_info_occ_apply_to_alt_conf_entry");
  if (gtk_toggle_button_get_active(togglebutton)) { 
    gtk_widget_set_sensitive(occ_entry,      TRUE);
    gtk_widget_set_sensitive(alt_conf_entry, TRUE);
  } else {
    gtk_widget_set_sensitive(occ_entry,      FALSE);
    gtk_widget_set_sensitive(alt_conf_entry, FALSE);
  }

}




void
on_residue_info_b_factor_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  GtkWidget *entry = lookup_widget(GTK_WIDGET(togglebutton), 
				   "residue_info_master_atom_b_factor_entry");
  
  if (gtk_toggle_button_get_active(togglebutton)) 
    gtk_widget_set_sensitive(entry, TRUE);
  else 
    gtk_widget_set_sensitive(entry, FALSE);
}



void
on_other_modelling_tools_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), 
				"other_model_tools_dialog");
   gtk_widget_destroy(w);
}


void
on_other_modelling_tools1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_other_model_tools_dialog();
   gtk_widget_show(w);
}


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


void
on_other_model_tools_dialog_destroy    (GtkWidget       *object,
                                        gpointer         user_data)
{
   do_cis_trans_conversion_setup(0);
   unset_other_modelling_tools_dialog();
}


void
on_undo_last_navigation1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   undo_last_move();
}


void
on_get_pdb_and_map_using_eds1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *window; 
  int n = COOT_ACCESSION_CODE_WINDOW_EDS;
  window = create_accession_code_window(); 
  g_object_set_data(G_OBJECT(window), "mode", GINT_TO_POINTER(n));
  gtk_widget_show(window); 
}


void
on_ligand_big_blob_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data)
{
  free_blob_dialog_memory(GTK_WIDGET(object));
}


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


void
on_simple1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *fileselection = coot_screendump_chooser();
   g_object_set_data(G_OBJECT(fileselection), "mode", GINT_TO_POINTER(COOT_SCREENDUMP_SIMPLE));
   set_directory_for_fileselection(fileselection);
   set_filename_for_filechooserselection(fileselection, "coot.png");
   set_file_selection_dialog_size(fileselection);
   add_ccp4i_project_optionmenu(fileselection, COOT_IMAGE_FILE_SELECTION);
   gtk_widget_show(fileselection);

   check_for_dark_blue_density(); /* give a dialog if density it too dark (blue) */

}


void
on_povray1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *fileselection = coot_screendump_chooser();
   g_object_set_data(G_OBJECT(fileselection), "mode", GINT_TO_POINTER(COOT_SCREENDUMP_POVRAY));
   set_directory_for_fileselection(fileselection);
   set_filename_for_filechooserselection(fileselection,
				   "coot-povray");
   set_file_selection_dialog_size(fileselection);
   add_ccp4i_project_optionmenu(fileselection, 
                                COOT_IMAGE_FILE_SELECTION);
   gtk_widget_show(fileselection);
}


void
on_raster3d1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *fileselection = coot_screendump_chooser();
   g_object_set_data(G_OBJECT(fileselection), "mode", GINT_TO_POINTER(COOT_SCREENDUMP_RASTER3D));
   set_directory_for_fileselection(fileselection);
   set_filename_for_filechooserselection(fileselection, 
				   "coot.png");
   set_file_selection_dialog_size(fileselection);
   add_ccp4i_project_optionmenu(fileselection, 
                                COOT_IMAGE_FILE_SELECTION);
   gtk_widget_show(fileselection);
}



void
on_screendump_image_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *fileselection = lookup_widget(GTK_WIDGET(button), 
					    "screendump_fileselection");
   gtk_widget_destroy(fileselection);

}


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


void
on_place_helix_here_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  place_helix_here();
}


void
on_other_tools_place_strand_here_button_clicked     (GtkButton       *button,
                                        gpointer         user_data) { 
  place_strand_here_dialog(); 	/* choose the python version in there, if needed. */
} 

void
on_diff_map_peaks_dialog_destroy       (GtkWidget       *object,
                                        gpointer         user_data)
{
  set_difference_map_peaks_widget(0); /* a null pointer */
}


void
on_symmetry_controller_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "symmetry_controller_dialog");
  gtk_widget_destroy(w);

}


void
on_molecule_0_checkbutton_toggled      (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) 
    set_show_symmetry_molecule(imol, 1);
  else 
    set_show_symmetry_molecule(imol, 0);

}


void
on_display_sphere_radiobutton_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) { 
    set_symmetry_whole_chain(imol, 0);
    symmetry_as_calphas(imol, 0); /* does an update_symmetry() */
  }
}


void
on_display_all_radiobutton_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) { 
    symmetry_as_calphas(imol, 0);
    set_symmetry_whole_chain(imol, 1);
/*   } else { */
/*     symmetry_as_calphas(imol, 1); */
/*     printf("DEBUG:: all for molecule %d CA state 1\n", imol); */
   }
}


void
on_display_CA_radiobutton_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  
  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) { 
     symmetry_as_calphas(imol, 1);
/*  } else { */
/*     printf("DEBUG:: CA for molecule %d CA state 0\n", imol); */
/*     symmetry_as_calphas(imol, 0); */
  } /* the off toggle of this button is deal with by the active state
       of other radio buttons. */
}


void
on_colour_symm_std_molecule_0_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_symmetry_colour_by_symop(imol, 0);
    set_symmetry_molecule_rotate_colour_map(imol, 0);
  }
}


void
on_colour_symm_by_symop_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_symmetry_molecule_rotate_colour_map(imol, 1); /* yes, I mean this */
    set_symmetry_colour_by_symop(imol, 1);
  }
}


void
on_colour_symm_by_molecule_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_toggle_button_get_active(togglebutton)) {
    set_symmetry_colour_by_symop(imol, 0);
    set_symmetry_molecule_rotate_colour_map(imol, 1);
  }
}


void
on_show_symmetry_molecule_control_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = symmetry_molecule_controller_dialog();
  gtk_widget_show(w);
}



void
on_ncs_controller_molecule_n_display_ncs_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  int state = 0;
  if (gtk_toggle_button_get_active(togglebutton)) { 
    state = 1;
    make_ncs_ghosts_maybe(imol);
  }
  /*    printf("NCS_controller Display NCS ghosts for imol %d %d\n", imol, state); */
  set_draw_ncs_ghosts(imol, state);
}


void
on_ncs_controller_molecule_n_display_chain_ich_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   int imol_chain = GPOINTER_TO_INT(user_data);
   int imol = imol_chain/1000;
   int ich = imol_chain - imol*1000;
   int state = 0;
   if (gtk_toggle_button_get_active(togglebutton)) { 
     state = 1;
   }
   printf("\nNCS_controller display chain toggled for imol %d chain %d state %d\n",  
	  imol, ich, state);
   ncs_control_display_chain(imol, ich, state);
}


void
on_ncs_controller_ncs_master_chain_ich_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(togglebutton), "ncs_control_dialog"); 
   int imol_chain = GPOINTER_TO_INT(user_data);
   int imol = imol_chain/1000;
   int ich = imol_chain - imol*1000;
/*    printf("==== DEBUG:: chain raiobutton toggled: imol %d ich %d active-state: %d \n",  */
/* 	  imol, gtk_toggle_button_get_active(ich, togglebutton)); */
   if (gtk_toggle_button_get_active(togglebutton)) { 
/*      printf("NCS_controller_ncs_master_chain_ich_radiobutton_toggled on for imol %d %d %d\n",  */
/* 	    imol_chain, imol, ich); */

/*      ncs_control_change_ncs_master_to_chain(imol, ich); (done in the following function) */

     ncs_control_change_ncs_master_to_chain_update_widget(w, imol, ich);
   } 
}


void
on_ncs_control_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "ncs_control_dialog"); 
   gtk_widget_destroy(w);
}


void
on_lsq_plane_add_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) { 
      setup_lsq_deviation(0);
      setup_lsq_plane_define(1);
   } else { 
      setup_lsq_deviation(0);
      setup_lsq_plane_define(1);
   } 
}


void
on_lsq_plane_deviant_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton)) { 
      setup_lsq_deviation(1);
      setup_lsq_plane_define(0);
   } else { 
      setup_lsq_deviation(1);
      setup_lsq_plane_define(0);
   } 
}


void
on_lsq_plane_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "lsq_plane_dialog");
  gtk_widget_destroy(w);
  normal_cursor();
}


void
on_lsq_plane_dialog_destroy            (GtkWidget       *object,
                                        gpointer         user_data)
{
  unset_lsq_plane_dialog();	/* which clears the plane points too */
  normal_cursor();
}


void
on_plane_distances1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_lsq_plane_dialog();
   setup_lsq_deviation(0);
   setup_lsq_plane_define(1);
   gtk_widget_show(w);
}


void
on_lsq_plane_delete_last_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  remove_last_lsq_plane_atom();
}


void
on_ncs_ghost_control1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_ncs_control_dialog();
   gtk_widget_show(w);
}


void
on_coords_colour_control_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{

}


void
on_coord_colour_control_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "coords_colour_control_dialog");
  gtk_widget_destroy(w);

}


void
on_bond_colours1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_coords_colour_control_dialog();
  gtk_widget_show(w);

}



void
on_coot_doc_url_monolithic_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  browser_url("https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/web/docs/coot.html");

}


void
on_coot_doc_url_sectioned_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  browser_url("https://www2.mrc-lmb.cam.ac.uk/Personal/pemsley/coot/web/docs/coot.html");
}


gboolean
on_coot_online_doc_search_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  const char *text;

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
    text = gtk_entry_get_text(GTK_ENTRY(widget));
    handle_online_coot_search_request(text);
  }
  
  return FALSE;
}


void
on_coot_online_doc_search_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkEntry *entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(button), "coot_online_doc_search_entry"));
  const char *text = gtk_entry_get_text(entry);
  handle_online_coot_search_request(text);
}


void
on_smiles1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_smiles_gui();
}


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


void
on_generic_display_objects1_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  generic_objects_gui_wrapper();

}


gboolean
on_entry1_key_press_event              (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{
  GtkEntry *entry = (GTK_ENTRY(lookup_widget(widget, "entry1")));
  const char *text = gtk_entry_get_text(entry);
  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) { 
    set_density_size_from_widget(text);
  }
  return FALSE;
}


void
on_map_radius_apply_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkEntry *entry = (GTK_ENTRY(lookup_widget(GTK_WIDGET(button), "entry1")));
   const char *text = gtk_entry_get_text(entry);
   set_density_size_from_widget(text);
}




void
on_refine_params_use_helix_peptide_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  /* not visible */
  printf("helix togglebutton toggled - ignored\n");
}


void
on_refine_params_use_beta_strand_peptide_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{				
  /* not visible */
  printf("beta strand togglebutton toggled - ignored\n");

}


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


void
on_probe_clashes1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "probe";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "probe_clashes1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_probe_clashes1_activate\n");
  }

}


void
on_validate1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *menu_item = 0;

  if (probe_available_p() == 0) { /* no */
    menu_item = lookup_widget(GTK_WIDGET(menuitem), "probe_clashes1");
    if (!menu_item) { 
      printf("Failed to get probe_clashes1 menu item :-(\n");
    } else { 
      /* desensitize it */
      gtk_widget_set_sensitive(menu_item, FALSE);
    }
  } 
}


void
on_spin_view_on_off1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  toggle_idle_spin_function();
}


void
on_rock_view_on_off1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 
  toggle_idle_rock_function();
} 


void
on_other_tools_RNA_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_nucleotide_builder_dialog(); 
  gtk_widget_show(w);
}

void
on_other_tools_base_pair_toggle_button_toggled      (GtkToggleButton       *button,
                                        gpointer         user_data) { 
  if (gtk_toggle_button_get_active(button))
    setup_base_pairing(1);
  else
    setup_base_pairing(0);
} 



void
on_ideal_rna_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "nucleotide_builder_dialog");
  ideal_nucleic_acid_by_widget(w);
  gtk_widget_destroy(w);
}


void
on_ideal_rna_cancel_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "nucleotide_builder_dialog");
  gtk_widget_destroy(w);
}


void
on_unit_cell_yes_radiobutton_toggled   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_show_unit_cells_all(1);
  else 
      set_show_unit_cells_all(0);
}


void
on_unit_cell_no_radiobutton_toggled    (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
      set_show_unit_cells_all(0);
  else 
      set_show_unit_cells_all(1);
}




void
on_move_molecule_here1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_move_molecule_here_dialog();
  gtk_widget_show(w);

}


void
on_move_molecule_here_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "move_molecule_here_dialog"); 
  move_molecule_here_by_widget(w);
  gtk_widget_destroy(w);
}


void
on_move_molecule_here_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "move_molecule_here_dialog"); 
  gtk_widget_destroy(w);

}


void
on_monomer_library_search_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "monomer_search_dialog");
  gtk_widget_destroy(w);
}


void
on_search_monomer_library1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = create_monomer_search_dialog();
  gtk_widget_show(w);
}


void
on_monomer_library_search_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *entry = lookup_widget(GTK_WIDGET(button), 
				   "monomer_search_entry");
  entry_char_type *text;
  GtkWidget *viewport = lookup_widget(GTK_WIDGET(button), 
				   "monomer_search_results_viewport");

  if (entry) { 
    text = gtk_entry_get_text(GTK_ENTRY(entry));
    if (text) {
      handle_make_monomer_search(text, viewport);
    }
  } 
}


void
on_monomer_search_entry_changed        (GtkEditable     *editable,
                                        gpointer         user_data)
{

}


gboolean
on_monomer_search_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  GtkWidget *entry = lookup_widget(GTK_WIDGET(widget), 
				   "monomer_search_entry");
  entry_char_type *text;
  GtkWidget *viewport = lookup_widget(GTK_WIDGET(widget), 
				   "monomer_search_results_viewport");

  if (entry) { 
    if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) { 
      text = gtk_entry_get_text(GTK_ENTRY(entry));
      if (text) {
	handle_make_monomer_search(text, viewport);
      }
    }
  } 
  return FALSE;
}


void
on_least_squares_match_type_all_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_least_squares_match_type_main_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_least_squares_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "least_squares_dialog");
  apply_lsq_matches_by_widget(w);

}


void
on_least_squares_close_button_clicked  (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "least_squares_dialog");
  update_lsq_dialog_store_values(w);
  gtk_widget_destroy(w);
}




void
on_least_squares_cancel_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "least_squares_dialog");
  gtk_widget_destroy(w);

}


void
on_lsq_superpose1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_least_squares_dialog();
  gtk_widget_show(w);

}


void
on_goto_atom_window_state_changed      (GtkWidget       *widget,
                                        GtkStateType     state,
                                        gpointer         user_data)
{
  /* how does this window get a state changed event?  Not at all, as
     far as I can see.  Nothing happens here. */
  printf("Go To Atom window state changed!\n");
}


gboolean
on_goto_atom_window_configure_event    (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{

  store_window_position(COOT_GO_TO_ATOM_WINDOW, widget);
  return FALSE;
}


gboolean
on_model_refine_dialog_configure_event (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{

  store_window_position(COOT_MODEL_REFINE_DIALOG, widget);
  return FALSE;
}


gboolean
on_display_control_window_glade_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_DISPLAY_CONTROL_WINDOW, widget);
  return FALSE;
}


gboolean
on_delete_item_dialog_configure_event  (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_DELETE_WINDOW, widget);
  return FALSE;
}


gboolean
on_rotate_translate_obj_dialog_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_ROTATE_TRANSLATE_DIALOG, widget);
  return FALSE;
}


gboolean
on_accept_reject_refinement_dialog_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  store_window_position(COOT_ACCEPT_REJECT_WINDOW, widget);
  return FALSE;
}


gboolean
on_dynarama_window_configure_event     (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data)
{
  resize_rama_canvas(widget, event);
   store_window_position(COOT_RAMACHANDRAN_PLOT_WINDOW, widget);
   return FALSE;
}


void
on_coords_fileselection1_destroy       (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_dataset_fileselection1_destroy      (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_map_name_fileselection1_destroy     (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_phs_coordinates_fileselection_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_save_coords_fileselection1_destroy  (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_save_symmetry_coords_fileselection_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_save_state_fileselection_destroy    (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}


void
on_screendump_fileselection_destroy    (GtkWidget       *object,
                                        gpointer         user_data)
{
   store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
}

void
on_residue_type_chooser_stub_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton)))
    set_residue_type_chooser_stub_state(1);
  else
    set_residue_type_chooser_stub_state(0);
}


void
on_set_undo_molecule_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  show_set_undo_molecule_chooser();
}


void
on_gln_and_asn_b_factor_outliers1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  char *type = "gln_and_asn_b_factor_outliers";
  GtkWidget *menu = lookup_widget(GTK_WIDGET(menuitem), "gln_and_asn_b_factor_outliers1");
  if (menu) { 
    add_on_validation_graph_mol_options(menu, type);
  } else { 
    printf("failed to get menu in on_gln_and_asn_b_factor_outliers1_activate\n");
  }

}




void
on_sec_str_rest_no_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_secondary_structure_restraints_type(0);
}


void
on_sec_str_rest_helix_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_secondary_structure_restraints_type(1);
}


void
on_sec_str_rest_strand_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_secondary_structure_restraints_type(2);
}



void
on_refine_params_dialog_destroy        (GtkWidget       *object,
                                        gpointer         user_data)
{
  unset_refine_params_dialog();
}


void
on_update_go_to_atom_from_current_position_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   update_go_to_atom_from_current_position(); 
}


/* gtk2 extras */
void
on_rz_simple_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 
}

void
on_rz_start_multizone1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data) {
}

void
on_rz_end_multizone_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data){

}

void
on_display_manager_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
}


void
on_reset_view_button_clicked           (GtkButton       *button,
                                        gpointer         user_data) { 
   reset_view();
} 

void
on_display_manager_toolbutton_clicked  (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
}

void
on_reset_view_toolbutton_clicked       (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
   reset_view();
}

void
on_symmetry_colorbutton_color_set      (GtkColorButton  *colorbutton,
                                        gpointer         user_data) {

  GdkColor colour;
  gdouble color[4]; // use first 3
  double r = 1.0 / 65535.0;
  gtk_color_button_get_color(colorbutton, &colour);
  color[0] = colour.red   * r;
  color[1] = colour.green * r;
  color[2] = colour.blue  * r;
  handle_symmetry_colour_change(1,color);

}

void
on_display_control_all_maps_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(togglebutton))
    set_all_maps_displayed(1);
  else
    set_all_maps_displayed(0);

}


void
on_display_control_all_models_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_button_get_active(togglebutton))
    set_all_models_displayed_and_active(1);
  else 
    set_all_models_displayed_and_active(0);

}


void
on_single_map_properties_contour_level_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "single_map_properties_dialog");
  single_map_properties_apply_contour_level_to_map(w); /* check now
							  made here
							  for valid
							  map
							  molecule. */
}


void
on_checked_waters_baddies_dialog_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{
  store_checked_waters_baddies_dialog(NULL);
}


void
on_model_toolbar_style_changed         (GtkToolbar      *toolbar,
                                        GtkToolbarStyle  style,
                                        gpointer         user_data)
{

  /* this does not do anything and doesnt need to */


}


gboolean
on_model_toolbar_button_press_event    (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data)
{

  if (event->type == GDK_BUTTON_PRESS && event->button == 3) {
    toolbar_popup_menu(GTK_TOOLBAR(widget), user_data, event);
    return TRUE;
  }

  return FALSE;
}


void
on_model_toolbar_refine_control_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = wrapped_create_refine_params_dialog();
  gtk_widget_show(widget);
}


void
on_model_toolbar_select_map_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   show_select_map_dialog();
}


void
on_model_toolbar_refine_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_refine(1);
  else 
    do_refine(0);		/* unclick button */
    
}

void
on_model_toolbar_regularize_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_regularize(1);
  else 
    do_regularize(0);		/* unclick button */
}

void
on_model_toolbar_fixed_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_fixed_atom_dialog();
  gtk_widget_show(w);
}



void
on_model_toolbar_rigid_body_fit_togglebutton_toggled
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


void
on_model_toolbar_rot_trans_togglebutton_toggled
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

void
on_model_toolbar_auto_fit_rotamer_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
     setup_auto_fit_rotamer(1);
  else 
    setup_auto_fit_rotamer(0);
}


void
on_model_toolbar_rotamers_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active)
      setup_rotamers(1);
   else 
      setup_rotamers(0);
}


#if (GTK_MAJOR_VERSION > 1) 
void
on_model_toolbar_edit_chi_angles_togglebutton_toggled
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
#endif	/* GTK_MAJOR_VERSION */


void
on_model_toolbar_torsion_general_toggletoolbutton_toggled
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


void
on_model_toolbar_flip_peptide_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_pepflip(1);
  else 
     do_pepflip(0);
}


void
on_model_toolbar_sidechain_180_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    setup_180_degree_flip(1);
  else 
    setup_180_degree_flip(0);
}


void
on_model_toolbar_edit_backbone_torsions_toggletoolbutton_toggled
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


void
on_model_toolbar_mutate_and_autofit_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    setup_mutate_auto_fit(1);
  else 
     setup_mutate_auto_fit(0);
}


void
on_model_toolbar_simple_mutate_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
   gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
   if (active)
      setup_mutate(1);
   else 
      setup_mutate(0);
}



void
on_model_toolbar_find_water_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *widget = create_find_waters_dialog();
   fill_find_waters_dialog(widget);
   gtk_widget_show(widget);
}


void
on_model_toolbar_add_terminal_residue_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggletoolbutton));
  if (active)
    do_add_terminal_residue(1);
  else 
    do_add_terminal_residue(0);
}



void
on_model_toolbar_add_alt_conf_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  altconf();

}


void
on_model_toolbar_add_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   place_atom_at_pointer();
}


void
on_model_toolbar_clear_pending_picks_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   clear_pending_picks();
}


void
on_model_toolbar_delete_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *widget = wrapped_create_delete_item_dialog();
  gtk_widget_show(widget);
}


void
on_model_toolbar_undo_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
   apply_undo();
}


void
on_model_toolbar_redo_button_clicked   (GtkButton       *button,
                                        gpointer         user_data)
{
  apply_redo();
}


void
on_model_toolbar_refmac_button_clicked (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  wrapped_create_run_refmac_dialog();

}


void
on_model_toolbar_icons_and_text1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  /*
  GtkWidget *toolbar = lookup_widget(GTK_WIDGET(menuitem), "model_toolbar");
  if (GTK_CHECK_MENU_ITEM(menuitem)->active){
      gtk_toolbar_set_style (GTK_TOOLBAR (toolbar), GTK_TOOLBAR_BOTH_HORIZ);
      GtkWidget *button;
      button = lookup_widget(GTK_WIDGET(toolbar), "model_toolbar_refine_control_button");
      gtk_button_set_label(GTK_BUTTON(button), "Refine/Regularize Control...");
      button = lookup_widget(GTK_WIDGET(toolbar), "model_toolbar_select_map_button");
      gtk_button_set_label(GTK_BUTTON(button), "Select Map...");
  }
  */

}

void
on_model_toolbar_icons1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  /*
  GtkWidget *toolbar = lookup_widget(GTK_WIDGET(menuitem), "model_toolbar");
  if (GTK_CHECK_MENU_ITEM(menuitem)->active){
    gtk_toolbar_set_style (GTK_TOOLBAR (toolbar), GTK_TOOLBAR_ICONS); 
    GtkWidget *button;
    button = lookup_widget(GTK_WIDGET(toolbar), "model_toolbar_refine_control_button");
    gtk_button_set_label(GTK_BUTTON(button), "R/RC");
    button = lookup_widget(GTK_WIDGET(toolbar), "model_toolbar_select_map_button");
    gtk_button_set_label(GTK_BUTTON(button), "Map");
  }
  */
}


void
on_model_toolbar_text1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  /*
  GtkWidget *toolbar = lookup_widget(GTK_WIDGET(menuitem), "model_toolbar");
  if (GTK_CHECK_MENU_ITEM(menuitem)->active){
    gtk_toolbar_set_style (GTK_TOOLBAR (toolbar), GTK_TOOLBAR_TEXT); 
    GtkWidget *button;
    button = lookup_widget(GTK_WIDGET(toolbar), "model_toolbar_refine_control_button");
    gtk_button_set_label(GTK_BUTTON(button), "Refine/Regularize Control...");
    button = lookup_widget(GTK_WIDGET(toolbar), "model_toolbar_select_map_button");
    gtk_button_set_label(GTK_BUTTON(button), "Select Map...");
  }
  */
}


void
on_model_toolbar_main_icons_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
    show_model_toolbar_main_icons();
  }
}


void
on_model_toolbar_all_icons_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(menuitem))) {
    show_model_toolbar_all_icons();
  }
}


void
on_model_toolbar_user_defined1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data) {

}


void
on_model_toolbar_setting1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  update_model_toolbar_icons_menu();
}


void
on_model_toolbar_menutoolbutton1_show_menu
                                        (GtkMenuToolButton *menutoolbutton,
                                        gpointer         user_data)
{
  /* I dont think anything needs to happen here */

}

void
on_model_toolbar_display_manager_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data)
{
  if (gtk_toggle_tool_button_get_active(toggletoolbutton)) {
	g_print("BL DEBUG:: display menu toggled");
  }
  
}


void
on_toolbar_display_manager_maps_all_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_toolbar_display_manager_molecules_all_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

void
on_scripting_python1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 

  post_python_scripting_window();

} 


void
on_scripting_scheme1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 
  post_scheme_scripting_window();

} 

void
on_aboutdialog_close                   (GtkDialog       *dialog,
                                        gpointer         user_data)
{
/* not this callback... */
}
void
on_aboutdialog_response                (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{
  gtk_widget_destroy(GTK_WIDGET(dialog));
}
void
on_check_waters_check1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data) {
}

void
on_check_waters_delete1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data){ 
}


/* start of chooser insert */
void
on_coords_filechooserdialog1_response  (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data){

 if (response_id == GTK_RESPONSE_OK) {
  const gchar *filename; 
  GtkWidget *coords_fileselection1;
  GtkWidget *combobox;
  int recentre_on_read_pdb_flag = 0;
  short int move_molecule_here_flag = 0;
  int active_index;
  GSList *sel_files;  
/*   GFile  *gfile; for xxx_get_files(), which we don't use (too modern) */

  coords_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
                                        "coords_filechooserdialog1");
  
  combobox = lookup_widget(GTK_WIDGET(dialog), 
                              "coords_filechooserdialog1_recentre_combobox");
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


  sel_files = 
    gtk_file_chooser_get_filenames(GTK_FILE_CHOOSER(coords_fileselection1));

  while (sel_files) { 

    filename = (char *) sel_files->data;

    printf("DEBUG:: filename: %s\n", filename);
   
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
  
  gtk_widget_destroy(coords_fileselection1);

 } else {
  GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
                                                "coords_filechooserdialog1");

  gtk_widget_destroy(coords_fileselection1);

 }
}

void
on_coords_filechooserdialog1_destroy  (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(object),
                                                "coords_filechooserdialog1");

  gtk_widget_destroy(coords_fileselection1);
}


void
on_coords_filechooserdialog1_recentre_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
    set_recentre_on_read_pdb(1);
   } else { 
    set_recentre_on_read_pdb(0);
   }
}


void
on_dataset_filechooserdialog1_response (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{

   if (response_id == GTK_RESPONSE_OK) {
      const gchar *filename; 
      gchar *copied_filename;
      int auto_read_flag = 0, ismtz = 0, ismtzauto = 0, iscnsauto = 0;

      GtkWidget *dataset_fileselection1;

      dataset_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
					     "dataset_filechooserdialog1");

      save_directory_from_filechooser(dataset_fileselection1);
      filename = gtk_file_chooser_get_filename 
	 (GTK_FILE_CHOOSER(dataset_fileselection1));
   
      /*    printf("dataset filename: %s\n", filename); */

      copied_filename = (char *) malloc(strlen(filename) + 1);
      strcpy(copied_filename, filename);

      auto_read_flag = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dataset_fileselection1), "auto_read_flag"));
      ismtz = is_mtz_file_p(filename);
      if (ismtz) ismtzauto = mtz_file_has_phases_p(filename);
      else       iscnsauto = cns_file_has_phases_p(filename);

      if ( ismtzauto || iscnsauto ) {

	 if (auto_read_flag) 
	    wrapped_auto_read_make_and_draw_maps(filename);
	 else 
	    /* this does a create_column_label_window, fills and displays it. */
	    manage_column_selector(copied_filename);
      
      } else { 

	 /* no phases path */
	 if (auto_read_flag) printf ("INFO:: This file is not a map coefficient file. Coot can auto-read\nINFO::  - MTZ files from refmac, phenix.refine, phaser, parrot, dm.\nINFO::  - CNS files (new 2009 format only) with cell, symops, F1, F2.\n");
	 if ( ismtz )
	    calc_phases_generic(filename);
	 else 
	    /* try to read as a phs, cif etc... */
	    manage_column_selector(copied_filename);
      }
      gtk_widget_destroy(dataset_fileselection1);
      free(copied_filename);
   } else {
      GtkWidget *dataset_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
							"dataset_filechooserdialog1");

      gtk_widget_destroy(dataset_fileselection1);
   }

}


void
on_dataset_filechooserdialog1_destroy (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *dataset_fileselection1 = lookup_widget(GTK_WIDGET(object),
                                                "dataset_filechooserdialog1");

  gtk_widget_destroy(dataset_fileselection1);
}


void
on_map_name_filechooserdialog1_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{

   if (response_id == GTK_RESPONSE_OK) { 
      const char* filename; 
      char *sfile; 
      GtkWidget* map_name_fileselection1; 
      GtkWidget *checkbutton;
      short int is_diff_map_flag = 0;
   
      map_name_fileselection1 = GTK_WIDGET(lookup_widget(GTK_WIDGET(dialog), 
							 "map_name_filechooserdialog1"));
      save_directory_from_filechooser(map_name_fileselection1);

      /* I don't think that we need to malloc this. */

      checkbutton = lookup_widget(GTK_WIDGET(dialog), 
				  "map_filechooser_is_difference_map_button");

      if (checkbutton)
	if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) 
	    is_diff_map_flag = 1;

      filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(map_name_fileselection1)); 
   

      printf("CCP4 map filename: %s\n", filename); 
      sfile = (char *) malloc (1001); 
      strncpy(sfile, filename, 1000); 

      gtk_widget_destroy(map_name_fileselection1); /* the file browser,
						      when destroyed,
						      scribbles over
						      filename. */
      handle_read_ccp4_map(sfile, is_diff_map_flag);
      free(sfile);

   } else {
      GtkWidget *map_name_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
							 "map_name_filechooserdialog1");

      gtk_widget_destroy(map_name_fileselection1);
   }

}


void
on_map_filechooser_is_difference_map_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_map_name_filechooserdialog1_destroy (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *map_name_fileselection1 = lookup_widget(GTK_WIDGET(object),
                                                "map_name_filechooserdialog1");

  gtk_widget_destroy(map_name_fileselection1);
}


void
on_phs_coordinates_filechooserdialog1_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{

   GtkWidget *phs_fileselection; 
   phs_fileselection = lookup_widget(GTK_WIDGET(dialog), 
                                     "phs_coordinates_filechooserdialog1");
   if (response_id == GTK_RESPONSE_OK) {
     const char *filename;    
     
     filename = gtk_file_chooser_get_filename
       (GTK_FILE_CHOOSER(phs_fileselection));
     
     save_directory_from_filechooser(phs_fileselection);
     read_phs_and_coords_and_make_map(filename); 
   } 
   gtk_widget_destroy(phs_fileselection);

}


void
on_phs_coordinates_filechooserdialog1_destroy
                                        (GtkWidget       *object,
                                        gpointer         user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *phs_fileselection1 = lookup_widget(GTK_WIDGET(object),
                                                "phs_coordinates_filechooserdialog1");

  gtk_widget_destroy(phs_fileselection1);
}


GtkFileChooserConfirmation
on_save_coords_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser, 
					gpointer user_data)
{

  if (file_chooser_overwrite_state() == 1) {
    return GTK_FILE_CHOOSER_CONFIRMATION_CONFIRM;
  } else {
    return GTK_FILE_CHOOSER_CONFIRMATION_ACCEPT_FILENAME;
  }

}


void
on_save_coords_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data)
{
  if (response_id == GTK_RESPONSE_OK) {
    char *stuff;
    GtkWidget *fileselection = lookup_widget(GTK_WIDGET(dialog), "save_coords_filechooserdialog1");

    save_directory_for_saving_from_filechooser(fileselection);
    stuff = g_object_get_data(G_OBJECT(fileselection), "stuff");
    save_coordinates_using_widget(fileselection);
    free(stuff);
    gtk_widget_destroy(fileselection); 
  } else {
    GtkWidget *fileselection = lookup_widget(GTK_WIDGET(dialog),
                                                "save_coords_filechooserdialog1");

    gtk_widget_destroy(fileselection);
  }
}


void
on_save_coords_filechooserdialog1_destroy
					(GtkWidget * object, 
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = lookup_widget(GTK_WIDGET(object),
                                                "save_coords_filechooserdialog1");

  gtk_widget_destroy(fileselection);
}


void
on_cif_dictionary_filechooserdialog1_response(GtkDialog * dialog,
					      gint response_id,
					      gpointer user_data) {

  int imol_enc = -999997;	/* unset */
  const char *filename;
  GtkWidget *fileselection;
  GtkWidget *dictionary_molecule_selector_option_menu;
  GtkWidget *menu;
  GtkWidget *active_menu_item;
  GtkWidget *checkbutton = lookup_widget(GTK_WIDGET(dialog),
					 "cif_dictionary_file_selector_create_molecule_checkbutton");
  short int new_molecule_checkbutton_state = 0;
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton)))
     new_molecule_checkbutton_state = 1;

  if (response_id == GTK_RESPONSE_OK) {

     fileselection = lookup_widget(GTK_WIDGET(dialog), "cif_dictionary_filechooserdialog1");
     save_directory_from_filechooser(fileselection);
     filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));

     dictionary_molecule_selector_option_menu =
	lookup_widget(GTK_WIDGET(dialog),
		      "cif_dictionary_file_selector_molecule_select_option_menu");

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

     handle_cif_dictionary_for_molecule(filename, imol_enc,
					new_molecule_checkbutton_state);

     gtk_widget_destroy(fileselection);
} else {
   fileselection = lookup_widget(GTK_WIDGET(dialog),
				 "cif_dictionary_filechooserdialog1");

   gtk_widget_destroy(fileselection);
 }

}


void
on_cif_dictionary_filechooserdialog1_destroy
					(GtkWidget * object, 
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = lookup_widget(GTK_WIDGET(object),
                                                "cif_dictionary_filechooserdialog1");

  gtk_widget_destroy(fileselection);
}


void
on_run_script_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data)
{

  if (response_id == GTK_RESPONSE_OK) {
    GtkWidget *fileselection = lookup_widget(GTK_WIDGET(dialog),
					   "run_script_filechooserdialog1");

    const char *script_filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));
    run_script(script_filename);
    gtk_widget_destroy(fileselection);

  } else {
    GtkWidget *fileselection = lookup_widget(GTK_WIDGET(dialog),
                                                "run_script_filechooserdialog1");

    gtk_widget_destroy(fileselection);
  }

}


void
on_run_script_filechooserdialog1_destroy
					(GtkWidget * object, 
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = lookup_widget(GTK_WIDGET(object),
                                                "run_script_filechooserdialog1");

  gtk_widget_destroy(fileselection);
}


#if (GTK_MAJOR_VERSION > 1) && (GTK_MINOR_VERSION > 9)
GtkFileChooserConfirmation
on_save_symmetry_coords_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser, 
					gpointer user_data)
{

  if (file_chooser_overwrite_state() == 1) {

    return GTK_FILE_CHOOSER_CONFIRMATION_CONFIRM;

  } else {

    return GTK_FILE_CHOOSER_CONFIRMATION_ACCEPT_FILENAME;

  }

}
#endif /* GTK_MAJOR_VERSION */


void
on_save_symmetry_coords_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data)
{
  if (response_id == GTK_RESPONSE_OK) {
    GtkWidget *w = lookup_widget(GTK_WIDGET(dialog), "save_symmetry_coords_filechooserdialog1");
    save_symmetry_coords_from_fileselection(w);
    gtk_widget_destroy(w);
  } else {
    GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
                                                "save_symmetry_coords_filechooserdialog1");

    gtk_widget_destroy(coords_fileselection1);

  }
}


void
on_save_symmetry_coords_filechooserdialog1_destroy
					(GtkWidget * object, 
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(object),
                                                "save_symmetry_coords_filechooserdialog1");

  gtk_widget_destroy(coords_fileselection1);
}


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


void
on_save_state_filechooserdialog1_response (GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data)
{
  if (response_id == GTK_RESPONSE_OK) {
   GtkWidget *w = lookup_widget(GTK_WIDGET(dialog),
				"save_state_filechooserdialog1");

   const char *filename = gtk_file_chooser_get_filename 
     (GTK_FILE_CHOOSER(w));

   save_state_file(filename);    /* write the file */
   set_save_state_file_name(filename); /* save as a static in graphics_info_t */
   gtk_widget_destroy(w);
  } else {
    GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(dialog),
						   "save_state_filechooserdialog1");

    gtk_widget_destroy(coords_fileselection1);
  }
}


void
on_save_state_filechooserdialog1_destroy (GtkWidget * object, 
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *coords_fileselection1 = lookup_widget(GTK_WIDGET(object),
                                                "save_state_filechooserdialog1");

  gtk_widget_destroy(coords_fileselection1);
}


#if (GTK_MAJOR_VERSION > 1) && (GTK_MINOR_VERSION > 9)
GtkFileChooserConfirmation
on_screendump_filechooserdialog1_confirm_overwrite 
					(GtkFileChooser * filechooser, 
					gpointer user_data)
{

  if (file_chooser_overwrite_state() == 1) {

    return GTK_FILE_CHOOSER_CONFIRMATION_CONFIRM;

  } else {

    return GTK_FILE_CHOOSER_CONFIRMATION_ACCEPT_FILENAME;

  }

}
#endif /* GTK_MAJOR_VERSION */


void
on_screendump_filechooserdialog1_response (GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data)
{
#if (GTK_MAJOR_VERSION > 1)
  if (response_id == GTK_RESPONSE_OK) {

   GtkWidget *fileselection = lookup_widget(GTK_WIDGET(dialog), 
					    "screendump_filechooserdialog1");
   int image_type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(fileselection), "image_type"));
   const char *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));

   if (image_type == COOT_SCREENDUMP_SIMPLE) { 
      screendump_image(filename);
   } 
   if (image_type == COOT_SCREENDUMP_POVRAY) { 
      make_image_povray(filename);
   } 
   if (image_type == COOT_SCREENDUMP_RASTER3D) { 
      make_image_raster3d(filename);
   } 
   gtk_widget_destroy(fileselection);

  } else {
    GtkWidget *fileselection = lookup_widget(GTK_WIDGET(dialog),
                                                "screendump_filechooserdialog1");

    gtk_widget_destroy(fileselection);
  }
#endif /* GTK_MAJOR_VERSION  */
}


void
on_screendump_filechooserdialog1_destroy (GtkWidget * object, 
					gpointer user_data)
{

  store_window_size(COOT_FILESELECTION_DIALOG, GTK_WIDGET(object));
  GtkWidget *fileselection = lookup_widget(GTK_WIDGET(object),
                                                "screendump_filechooserdialog1");

  gtk_widget_destroy(fileselection);
}
/* end of chooser insert */


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

void
on_accept_reject_reverse_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  toggle_torsion_general_reverse();
}

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

void
on_accept_reject_atom_pull_clear_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data) { 
   clear_all_atom_pull_restraints();
}



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


void
on_clear_fixed_atoms_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
  clear_fixed_atoms_all();
}


void
on_fixed_atom_close_button_clicked     (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button),
				    "fixed_atom_dialog");
  gtk_widget_destroy(dialog);
}


void
on_fixed_atom_dialog_destroy           (GtkWidget       *object,
                                        gpointer         user_data)
{
  store_fixed_atom_dialog(0);
}

void
on_add_rep_add_rep_button_clicked      (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "add_reps_dialog");
  add_additional_representation_by_widget(w);
  gtk_widget_destroy(w);
}


void
on_add_rep_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "add_reps_dialog");
  gtk_widget_destroy(w);
}

void
on_additional_representation1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_add_additional_representation_gui();
  gtk_widget_show(w);
}

void
on_display_additional_representations_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}

void
on_all1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_all2_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

void
on_residue_editor_select_monomer_type_ok_button_clicked (GtkButton       *button,
						    gpointer         user_data) { 
  

  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "residue_editor_select_monomer_type_dialog");
  GtkWidget *combo_box = lookup_widget(GTK_WIDGET(button), "residue_editor_select_monomer_type_combobox");
  gchar *t = 0;

  t = gtk_combo_box_get_active_text(GTK_COMBO_BOX(combo_box));
  show_restraints_editor(t);
  gtk_widget_destroy(dialog);
}


void
on_residue_editor_select_monomer_type_cancel_button_clicked (GtkButton       *button,
							gpointer         user_data) {

  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "residue_editor_select_monomer_type_dialog");
  gtk_widget_destroy(dialog);
}


void
on_restraints1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_residue_editor_select_monomer_type_dialog();
  gtk_widget_show(w);
}


void
on_restraint_editor_add_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  restraints_editor_add_restraint_by_widget(w);
}


void
on_restraints_editor_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  if (w) { 
    clear_restraints_editor_by_dialog(w);
    gtk_widget_destroy(w);
  }

}

void
on_restraints_editor_save_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  restraints_editor_save_restraint_by_widget(w);
}

void
on_restraints_editor_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  apply_restraint_by_widget(w);
}

void
on_restraint_editor_delete_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  restraints_editor_delete_restraint_by_widget(w);
}



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

void
on_save_restraint_chooserdialog_response(GtkDialog       *dialog,
					 gint             response_id,
					 gpointer         user_data) { 
/* Maybe there are responses other than OK and cancel, so don't factor
   out the destroy() */
  GtkWidget *w = lookup_widget(GTK_WIDGET(dialog), "save_restraint_chooserdialog");
  if (response_id == GTK_RESPONSE_OK) {
    save_monomer_restraints_by_widget(dialog);
    gtk_widget_destroy(w);
  }
  if (response_id == GTK_RESPONSE_CANCEL) {
    gtk_widget_destroy(w);
  }
}

/* This is not the way. */
void
on_save_restraint_chooserdialog_close  (GtkDialog       *dialog,
                                        gpointer         user_data) { 
}


void
on_model_refine_dialog_fix_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = wrapped_create_fixed_atom_dialog();
  gtk_widget_show(w);
}


void
on_model_refine_dialog_fast_sss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog;
  dialog = wrapped_create_fast_ss_search_dialog();
  gtk_widget_show(dialog);

}


void
on_other_tools_build_na_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w;
   w = create_build_na_dialog();
   gtk_widget_show(w); 

}


void
on_build_na_dialog_cancelbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w;
   w = lookup_widget(GTK_WIDGET(button), "build_na_dialog");
   gtk_widget_destroy(w);

}


void
on_build_na_dialog_okbutton_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{
   GtkWidget *w;
   GtkEntry *entry;
   const char *text;
   float r;
   w = lookup_widget(GTK_WIDGET(button), "build_na_dialog");
   entry = (GTK_ENTRY(lookup_widget(GTK_WIDGET(button), 
                                    "build_na_dialog_radius_entry")));
   text = gtk_entry_get_text(entry);
   r = atof(text);
   find_nucleic_acids_local(r);
   gtk_widget_destroy(w);

}


void
on_build_na_dialog_radius_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data)
{
   /* BL note:: this is almost a repetetion of 
      on_build_na_dialog_okbutton_clicked
      we could/should have a function for this?!
    */
   GtkWidget *w;
   const char *text;
   float r;
   w = lookup_widget(GTK_WIDGET(entry), "build_na_dialog");
   text = gtk_entry_get_text(entry);
   r = atof(text);
   find_nucleic_acids_local(r);
   gtk_widget_destroy(w);

}


void
on_coot_references_coot_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_COOT);

}

void
on_coot_references_wincoot_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_WINCOOT);

}


void
on_coot_references_refmac_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_REFMAC);

}


void
on_coot_references_ssm_toolbutton_clicked  
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_SSM);

}


void
on_coot_references_mmdb_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_MMDB);

}

void
on_coot_references_clipper_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_CLIPPER);

}


void
on_coot_references_buccaneer_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_BUCCANEER);

}


void
on_coot_references_molprobity_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_MOLPROBITY);

}


void
on_coot_references_calpha_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_CALPHA);

}


void
on_coot_references_xligand_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_XLIGAND);

}

void
on_coot_references_eds_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_EDS);

}


void
on_coot_references_others_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  fill_references_notebook(toolbutton, COOT_REFERENCE_OTHERS);

}


void
on_coot_references_closebutton_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *dialog;
  dialog = lookup_widget(GTK_WIDGET(button), "coot_references_dialog");
  gtk_widget_destroy(dialog);

}




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


void
on_edit_chi_angles_add_hydrogen_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
   GtkWidget *vbox = lookup_widget(GTK_WIDGET(togglebutton),"edit_chi_angles_vbox");

   if (gtk_toggle_button_get_active(togglebutton)) { 
      set_find_hydrogen_torsions(1);
   } else { 
      set_find_hydrogen_torsions(0);
   } 
   fill_chi_angles_vbox(vbox);
}

#ifdef GTK_TYPE_MENU_TOOL_BUTTON
void
on_model_toolbar_rot_trans_toolbutton_show_menu
                                        (GtkMenuToolButton *toolbutton,
					 gpointer         user_data) { 

  GtkWidget *menu;
  GList *children;
  GtkWidget *menu_item = NULL;
  menu = gtk_menu_tool_button_get_menu(toolbutton);
  children = gtk_container_get_children(GTK_CONTAINER(menu));
  if (get_rot_trans_object_type() == ROT_TRANS_TYPE_ZONE) {
    menu_item = children->data;
  }
  if (get_rot_trans_object_type() == ROT_TRANS_TYPE_CHAIN) {
    children = children->next;
    menu_item = children->data;
  }
  if (get_rot_trans_object_type() == ROT_TRANS_TYPE_MOLECULE) {
    children = children->next;
    children = children->next;
    menu_item = children->data;
  }
  gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(menu_item), TRUE);
  
}
#endif

#ifdef GTK_TYPE_MENU_TOOL_BUTTON
void
on_model_toolbar_rot_trans_toolbutton_clicked
                                        (GtkMenuToolButton *toolbutton,
					 gpointer         user_data) { 

  printf("clicked!\n");
  do_rot_trans_setup(1);
} 
#endif


void
on_find_ligands_search_all_radiobutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data) { 
}


void
on_find_ligands_search_here_radiobutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data) {

  set_ligand_dialog_number_of_sites_sensitivity(GTK_WIDGET(button));

}

void
on_symmetry_controller_dialog_destroy  (GtkWidget       *object,
                                        gpointer         user_data)
{
  set_symmetry_controller_dialog_widget(0);
}


void
on_toolbar_multi_refine_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  toolbar_multi_refine_continue();
}


void
on_toolbar_multi_refine_stop_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  toolbar_multi_refine_stop();
}

void
on_toolbar_multi_refine_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) { 

  toolbar_multi_refine_cancel();
} 


  

void
on_map_sharpening1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_map_sharpening_dialog();
   gtk_widget_show(w);
}


void
on_map_sharpening_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data)
{

   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "map_sharpening_dialog");
   gtk_widget_destroy(w);

}

void
on_map_sharpening_optimize_button_clicked    ( GtkButton       *button,
                                        	gpointer         user_data)
{
	GtkWidget *w = lookup_widget(GTK_WIDGET(button), "map_sharpening_dialog");
	calc_and_set_optimal_b_factor(w);
}

void
on_map_sharpening_reset_button_clicked (GtkButton       *button,
                                        gpointer         user_data)
{
    // reset to zero!?
    GtkWidget *h_scale = lookup_widget(GTK_WIDGET(button), "map_sharpening_hscale");
    GtkAdjustment *adj = gtk_range_get_adjustment(GTK_RANGE(h_scale));
    gtk_adjustment_set_value(adj, 0.);

}


void
on_map_sharpening_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

}


void
on_baton_build_params_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "baton_build_params_dialog");
  set_baton_build_params_from_widget(w);
  gtk_widget_destroy(w);

}


void
on_baton_build_params_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "baton_build_params_dialog");
  gtk_widget_destroy(w);
}


/* BL says:: comment out whole function for GTK1 as we dont have the dialog */
void
on_baton_build_set_params_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) { 

  GtkWidget *w = create_baton_build_params_dialog();
  gtk_widget_show(w);

} 


void
on_coords_toolbutton_clicked           (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  open_coords_dialog();
}

void
on_go_to_atom_toolbutton_clicked       (GtkToolButton   *toolbutton,
                                        gpointer         user_data) { 
  
  GtkWidget *widget = wrapped_create_goto_atom_window(); 
  gtk_widget_show(widget); 
}

void
on_go_to_ligand_toolbutton_clicked     (GtkToolButton   *toolbutton,
                                        gpointer         user_data) { 
  go_to_ligand();
} 


void
on_move_molecule_here_big_molecules_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *dialog = lookup_widget(GTK_WIDGET(togglebutton), "move_molecule_here_dialog");
  fill_move_molecule_here_dialog(dialog);
}



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

gboolean
on_environment_distance_max_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
       execute_environment_settings(widget);
  }
  return FALSE;
}

gboolean
on_environment_distance_min_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
       execute_environment_settings(widget);
  }
  return FALSE;
}


void
on_pisa_interfces_close_button_clicked (GtkButton       *button,
                                        gpointer         user_data) { 

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "pisa_interfaces_dialog");
  gtk_widget_destroy(w);
} 

void
on_displayed_map_style_as_lines_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data) {
  

  GtkWidget *window = lookup_widget(GTK_WIDGET(togglebutton),
				    "single_map_properties_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
    set_draw_map_standard_lines(imol, 1);
    set_draw_solid_density_surface(imol, 0);
  }

} 

/* we should call this "third-map-mode" or something */
void
on_displayed_map_style_as_cut_glass_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data) { 

  GtkWidget *window = lookup_widget(GTK_WIDGET(togglebutton),
				    "single_map_properties_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
    set_draw_map_standard_lines(imol, 0);
    set_draw_solid_density_surface(imol, 1);
    set_flat_shading_for_solid_density_surface(1);
  }
} 


void
on_displayed_map_style_as_transparent_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data) { 
  GtkWidget *window = lookup_widget(GTK_WIDGET(togglebutton),
				    "single_map_properties_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) { 
    set_draw_map_standard_lines(imol, 0);
    set_draw_solid_density_surface(imol, 1);
    set_flat_shading_for_solid_density_surface(0);
  }
  
} 

void
on_map_opacity_hscale_value_changed    (GtkRange        *range,
                                        gpointer         user_data) { 

  GtkAdjustment *adjustment;
  float fvalue;
  GtkWidget *window = lookup_widget(GTK_WIDGET(range),
				    "single_map_properties_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));

  adjustment = gtk_range_get_adjustment(GTK_RANGE(range));
  fvalue = 0.01 * gtk_adjustment_get_value(adjustment);
  set_solid_density_surface_opacity(imol, fvalue);

}


void
on_refine_params_weight_matrix_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{

  GtkWidget *entry = lookup_widget(GTK_WIDGET(editable), 
				   "refine_params_weight_matrix_entry");
  struct entry_info_t ei = coot_entry_to_val(GTK_ENTRY(entry));
  if (ei.float_is_set)
    set_matrix(ei.val_as_float);
  
}


void
on_remarks_browser1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 

  GtkWidget *w = wrapped_create_remarks_browser_molecule_chooser_dialog();
  gtk_widget_show(w);

} 

void
on_remarks_browser_molecule_chooser_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) { 
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "remarks_browser_molecule_chooser_dialog");
  gtk_widget_destroy(w);

} 

void
on_remarks_browser_molecule_chooser_ok_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) { 

  show_remarks_browswer(); // there we look up which molecule to show.
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "remarks_browser_molecule_chooser_dialog");
  gtk_widget_destroy(w);

}


void
on_fix_nomenclature_errors_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "fix_nomenclature_errors_dialog");
  int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "imol"));
  fix_nomenclature_errors(imol);
  gtk_widget_destroy(w);
}


void
on_fix_nomenclature_errors_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "fix_nomenclature_errors_dialog");
  gtk_widget_destroy(w);
}


void
on_ligand_builder1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 
   start_ligand_builder_gui_internal(menuitem, user_data);
} 


void
on_multi_residue_torsion_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "multi_residue_torsion_dialog");
  clear_up_moving_atoms();
  clear_pending_picks(); /* emcompasses in_multi_residue_torsion_define (but not mode) */
  clear_multi_residue_torsion_mode();
  normal_cursor();
  gtk_widget_destroy(w);
}


void
on_multi_residue_torsion_OK_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "multi_residue_torsion_dialog");
  gtk_widget_destroy(w);
  accept_regularizement();
  clear_multi_residue_torsion_mode();
}

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



void
on_multi_residue_torsion_pick_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "multi_residue_torsion_pick_dialog");
  gtk_widget_destroy(w);
  clear_pending_picks(); /* emcompasses in_multi_residue_torsion_mode */
  normal_cursor();
  show_multi_residue_torsion_dialog();
}


void
on_multi_residue_torsion_pick_cancel_button_activate
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
}


void
on_multi_residue_torsion_pick_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "multi_residue_torsion_pick_dialog");
  gtk_widget_destroy(w);
  clear_pending_picks(); /* emcompasses in_multi_residue_torsion_define (but not mode) */
  clear_multi_residue_torsion_mode();
  normal_cursor();
}



/* This is the call-back for the button on the Other Modelling Tools dialog. */
void
on_multi_residue_torsion_start_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  
  setup_multi_residue_torsion(); // shows a dialog

}


/* wrong callback possibly */
void
on_keyboard_go_to_residue_entry_changed   (GtkEditable     *editable,
                                        gpointer         user_data)
{

}




gboolean
on_keyboard_go_to_residue_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(widget, "keyboard_goto_residue_window");
  const gchar *text = gtk_entry_get_text(GTK_ENTRY(widget));
  if (event->keyval == GDK_KEY_Return || event->keyval == GDK_KEY_KP_Enter) {
    handle_go_to_residue_keyboarding_mode(text);
    gtk_widget_destroy(w);
    return TRUE;
  }
  if (event->keyval == GDK_KEY_Escape) {
    gtk_widget_destroy(w);
    return TRUE;
  }
  return FALSE;
}

void
on_mogul_geometry_dialog_close_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) { 

   GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "mogul_geometry_results_table_dialog");

   /* And the histogram?  How do I look that up? */
   gtk_widget_destroy(dialog);
}


void 
on_ligand_check_okbutton_clicked(GtkButton       *button,
				 gpointer         user_data) { 

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "ligand_check_dialog");
  gtk_widget_destroy(w);

} 

void
on_generic_objects_dialog_closebutton_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) {

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "generic_objects_dialog");
  gtk_widget_destroy(w);
  clear_generic_objects_dialog_pointer();
  graphics_draw();

} 

/* I don't know how this function gets activated, it's not the close button of the dialog */
void
on_generic_objects_dialog_close        (GtkDialog       *dialog,
                                        gpointer         user_data) { 

  clear_generic_objects_dialog_pointer(); /* needed here? */
  graphics_draw();

}

void
on_generic_objects_dialog_destroy      (GtkWidget       *object,
                                        gpointer         user_data) { 

  clear_generic_objects_dialog_pointer();

} 


void
on_generic_objects_display_all_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  int state = 0;
  if (gtk_toggle_button_get_active(togglebutton))
    state = 1;
  set_display_all_generic_objects(state);

}


void
on_generic_objects_close_all_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  close_all_generic_objects();
}



void
on_export_map1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 

  short int is_fragment = 0;
  export_map_gui(is_fragment);

} 

void
on_export_map_fragment1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data) { 

  short int is_fragment = 1;
  export_map_gui(is_fragment);
}

void
on_export_map_dialog_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data) {

  on_export_map_dialog_ok_button_clicked_cc(button);
} 

void
on_export_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data) { 

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "export_map_dialog");
  gtk_widget_destroy(w);
} 

/* void */
/* on_export_map_filechooserdialog_cancel_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data) */
/* { */

/*   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "export_map_filechooserdialog"); */
/*   gtk_widget_destroy(w); */
  

/* } */

/* void */
/* on_export_map_filechooserdialog_save_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data) */
/* { */

/*   GtkWidget *w = lookup_widget(GTK_WIDGET(button), "export_map_filechooserdialog"); */
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
/*   gtk_widget_destroy(w); */

/* } */

void
on_export_map_filechooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
					 gpointer         user_data) { 

  int imol_map;
  int is_map_fragment;
  char *txt;
  const char *filename;

  if (response_id == GTK_RESPONSE_OK) {
    imol_map = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "map_molecule_number")); 
    is_map_fragment = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "is_map_fragment")); 
    txt = g_object_get_data(G_OBJECT(dialog), "export_map_radius_entry_text");
    filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog)); 

    if (is_map_fragment) { 
      export_map_fragment_with_text_radius(imol_map, txt, filename);
    } else { 
      export_map(imol_map, filename);
    } 
  }
  gtk_widget_destroy(GTK_WIDGET(dialog));
} 

#include "cfc-widgets-c-interface.h"

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


void
on_dynarama_outliers_only_togglebutton_toggled (GtkToggleButton *togglebutton,
						gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(togglebutton), "dynarama_window");
   toggle_dynarama_outliers(window, gtk_toggle_button_get_active(togglebutton)); /* get the imol from window */
}



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

void
on_edit_copy_molecule1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  do_edit_copy_molecule();
}


void
on_edit_copy_fragment1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_edit_copy_fragment();
}


void
on_edit_replace_residue1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_edit_replace_residue();
}


void
on_edit_replace_fragment1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  do_edit_replace_fragment();
}


void
on_edit_renumber_residues1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

  GtkWidget *w = wrapped_create_renumber_residue_range_dialog();
  gtk_widget_show(w);
}


void
on_edit_change_chain_ids1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

   GtkWidget *w = wrapped_create_change_chain_id_dialog();
   gtk_widget_show(w);
}


void
on_edit_merge_molecules1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   gtk_widget_show(w);
}


void
on_weight_maxtrix_estimate_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *entry = lookup_widget(GTK_WIDGET(button), "refine_params_weight_matrix_entry");
  estimate_map_weight(entry);

}


void
on_mutate_molecule_resno_1_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{

  GtkWidget *res_no_1_widget = lookup_widget(GTK_WIDGET(editable), "mutate_molecule_resno_1_entry");
  GtkWidget *res_no_2_widget = lookup_widget(GTK_WIDGET(editable), "mutate_molecule_resno_2_entry");
  GtkWidget *text_widget     = lookup_widget(GTK_WIDGET(editable), "mutate_molecule_sequence_text");
  GtkWidget *label_widget    = lookup_widget(GTK_WIDGET(editable), "mutate_residue_range_counts_label");
  mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);

}


void
on_mutate_molecule_resno_2_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data)
{

  GtkWidget *res_no_1_widget = lookup_widget(GTK_WIDGET(editable), "mutate_molecule_resno_1_entry");
  GtkWidget *res_no_2_widget = lookup_widget(GTK_WIDGET(editable), "mutate_molecule_resno_2_entry");
  GtkWidget *text_widget     = lookup_widget(GTK_WIDGET(editable), "mutate_molecule_sequence_text");
  GtkWidget *label_widget    = lookup_widget(GTK_WIDGET(editable), "mutate_residue_range_counts_label");
  mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
}


void
on_mutate_molecule_sequence_text_insert_at_cursor
                                        (GtkTextView     *textview,
                                        gchar           *string,
                                        gpointer         user_data)
{

  GtkWidget *res_no_1_widget = lookup_widget(GTK_WIDGET(textview), "mutate_molecule_resno_1_entry");
  GtkWidget *res_no_2_widget = lookup_widget(GTK_WIDGET(textview), "mutate_molecule_resno_2_entry");
  GtkWidget *text_widget     = lookup_widget(GTK_WIDGET(textview), "mutate_molecule_sequence_text");
  GtkWidget *label_widget    = lookup_widget(GTK_WIDGET(textview), "mutate_residue_range_counts_label");
  mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
}


gboolean
on_mutate_molecule_sequence_text_key_release_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data)
{

  GtkWidget *res_no_1_widget = lookup_widget(GTK_WIDGET(widget), "mutate_molecule_resno_1_entry");
  GtkWidget *res_no_2_widget = lookup_widget(GTK_WIDGET(widget), "mutate_molecule_resno_2_entry");
  GtkWidget *text_widget     = lookup_widget(GTK_WIDGET(widget), "mutate_molecule_sequence_text");
  GtkWidget *label_widget    = lookup_widget(GTK_WIDGET(widget), "mutate_residue_range_counts_label");
  mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
  return FALSE;
}


gboolean
on_mutate_molecule_sequence_text_button_release_event
                                        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data)
{

  GtkWidget *res_no_1_widget = lookup_widget(GTK_WIDGET(widget), "mutate_molecule_resno_1_entry");
  GtkWidget *res_no_2_widget = lookup_widget(GTK_WIDGET(widget), "mutate_molecule_resno_2_entry");
  GtkWidget *text_widget     = lookup_widget(GTK_WIDGET(widget), "mutate_molecule_sequence_text");
  GtkWidget *label_widget    = lookup_widget(GTK_WIDGET(widget), "mutate_residue_range_counts_label");
  mutate_molecule_dialog_check_counts(res_no_1_widget, res_no_2_widget, text_widget, label_widget);
  return FALSE;
}


void
on_display_control_last_model_only_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  set_only_last_model_molecule_displayed();

}


void
on_display_control_align_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

  align_labels_checkbutton_toggled(togglebutton);

}



void
on_curlew_install_button_clicked(GtkButton *button,
				 gpointer   user_data) {

  GtkWidget *dialog = lookup_widget(GTK_WIDGET(button), "curlew_dialog");
  int n_items = 0;
  if (dialog) {
    n_items = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "n_extensions"));
    curlew_dialog_install_extensions(dialog, n_items); /* some of which were selected */
  }
}



void
on_curlew_dialog_close                 (GtkDialog       *dialog,
                                        gpointer         user_data)
{
  gtk_widget_destroy(GTK_WIDGET(dialog)); /* or maybe hide */
}


void
on_curlew_dialog_response              (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{

  /* 
  printf("in on_curlew_dialog_response with response_id %d\n", response_id);
  printf("   cf response_id %d\n", GTK_RESPONSE_CLOSE);
  printf("   cf response_id %d\n", GTK_RESPONSE_OK);
  printf("   cf response_id %d\n", GTK_RESPONSE_CANCEL);
  */

  if (response_id == GTK_RESPONSE_CLOSE)
    gtk_widget_destroy(GTK_WIDGET(dialog));

}


void
on_modelling_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data) {

}


void
on_edit_settings_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}

void
on_calculate_all_molecule_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_dock_sequence_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_map_tools_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_modules_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_ncs_tools_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_pisa_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_draw_representation_tools_activate  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_views_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_calculate_load_tutorial_model_and_data1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  load_tutorial_model_and_data();
}


void
on_refine_params_geman_mcclure_alpha_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

   const char *t = gtk_combo_box_get_active_text(GTK_COMBO_BOX(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_geman_mcclure_alpha_from_text(active_item_idx, t);
}


void
on_refine_params_lennard_jones_epsilon_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
   const char *t = gtk_combo_box_get_active_text(GTK_COMBO_BOX(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_lennard_jones_epsilon_from_text(active_item_idx, t);
}


void
on_refine_params_rama_restraints_weight_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
   const char *t = gtk_combo_box_get_active_text(GTK_COMBO_BOX(combobox));
   int active_item_idx = gtk_combo_box_get_active(combobox);
   set_refinement_ramachandran_restraints_weight_from_text(active_item_idx, t);
}


void
on_refine_params_more_control_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

   if (togglebutton) {
      GtkWidget *frame = lookup_widget(GTK_WIDGET(togglebutton), "refine_params_more_control_frame");
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(togglebutton))) {
         gtk_widget_show(frame);
         set_refine_params_dialog_more_control_frame_is_active(1);
      } else {
         gtk_widget_hide(frame);
         set_refine_params_dialog_more_control_frame_is_active(0);
      }
   }
}


void
on_accept_reject_flip_this_peptide_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  pepflip_intermediate_atoms();
}


void
on_accept_reject_flip_next_peptide_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  pepflip_intermediate_atoms_other_peptide();
}


void
on_accept_reject_crankshaft_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  crankshaft_peptide_rotation_optimization_intermediate_atoms();
}


void
on_accept_reject_backrub_rotamer_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  backrub_rotamer_intermediate_atoms();
}

void
on_symmetry_always_on_checkbutton_toggled (GtkToggleButton *togglebutton,
					   gpointer         user_data) {

   GtkWidget *symmetry_on_radio_button = NULL;
   if (gtk_toggle_button_get_active(togglebutton)) {
      add_symmetry_on_to_preferences_and_apply();
      symmetry_on_radio_button = lookup_widget(GTK_WIDGET(togglebutton), "show_symmetry_yes_radiobutton");
      if (! gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(symmetry_on_radio_button)))
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(symmetry_on_radio_button), TRUE);
   }

}


void
on_curlew1_activate              (GtkMenuItem     *menuitem,
                                  gpointer         user_data) {

  curlew();

}

