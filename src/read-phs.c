/* src/read-phs.c
 * 
 * Copyright 2002, 2006 by The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifdef _MSC_VER
#define snprintf _snprintf
#endif

#include <string.h>
#include <stdio.h>


#include <gtk/gtk.h>


#include "c-interface.h"
#include "gtk-manual.h" // for display_cell_chooser_box
#include "read-phs.h"
#include "interface.h"

#include "c-interface-gtk-widgets.h" // for main_window()

#include "coot-fileselections.h" // for file chooser filter button

int
try_read_phs_file(const char *filename) { 

 /* if the extention of filename was ".phs", then try to read it as if
    it were a XtalView .phs file. */

   char *f_pt;
   GtkWidget *widget; 
   /* int is_pha_extension;   pha is a phs file too (shelx) */

   f_pt = strrchr(filename, '.'); 

   if ( f_pt != NULL ) { 

      if ( (! (strncmp(f_pt, ".phs", 4)))
	   ||
	   (! (strncmp(f_pt, ".pha", 4)))) {

	 /* we have a match */
	 printf("%s is a .phs file\n", filename); 

	 if (! possible_cell_symm_for_phs_file()) {  /* question */
	   widget = create_phs_info_box();	/* this returns a GtkWidget */
	   gtk_widget_show(widget); 
	 } else { 
	   do_phs_cell_choice_window(); 
	 }
	 graphics_store_phs_filename(filename); 
	 
	 return 1;			/* we try to deal with a .phs file */
	 
      } else {
	printf("%s is not a .phs file\n", filename); 
      }
   }
   
   return 0;			/* we do not try to deal with a .phs file */
}

void do_phs_cell_choice_window() { 

  GtkWidget *window; 
  GtkEntry  *entry; 
  int i; 
  int ii;
  gchar *txt; 
  GSList *phs_cell_group = NULL;
  GtkWidget *toggle_button;
  GtkWidget *mw_widget;
  GtkWindow *mw;

/* messing about with string variables */
  gchar entry_name[100]; 
  gchar widget_name[100];  		/* for the radiobutton */
  gchar *tmp_name; 
  int mol_count = 0; 
  int mol_with_cell = -1; /* unset */

  for (ii=0; ii<100; ii++) { 
     entry_name[ii] = 0;
    widget_name[ii] = 0;
  }
    
  
/*   entry_name = (gchar *) malloc(100); */
/*   widget_name = (gchar *) malloc(25);  */

  window = create_phs_cell_choice_window();

  printf("..... debug:: we got window 0x%p\n", window);
  fflush(stdout);
  printf("..... debug:: we got main_window 0x%p\n", main_window());
  fflush(stdout);
  if (! window) {
    printf("ERROR:: failed to get valid window from create_phs_cell_choice_window()\n");
    return;
  }
  if (! main_window()) {
    printf("ERROR:: failed to get main_window() from create_phs_cell_choice_window()\n");
    return;
  }
  
  mw_widget = lookup_widget(main_window(), "window1");

  if (! mw_widget) {
    printf("ERROR:: failed to get valid mw_widget in create_phs_cell_choice_window()\n");
    return;
  }

  
  printf("..... debug:: we got mw_widget 0x%p\n", mw_widget);
  fflush(stdout);

  mw = GTK_WINDOW(mw_widget);  

  for (i=0; i<graphics_n_molecules(); i++) {

    // only draw the add data if we have cell and symmetry info:
    if (has_unit_cell_state(i)) {

      phs_cell_group = display_cell_chooser_box(window, phs_cell_group, i); /* manipulate phs_cell_group */

      strcpy(entry_name, "phs_cell_symm_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i, "symm"); 
      gtk_entry_set_text(entry, txt); 
    
      strcpy(entry_name, "phs_cell_a_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i,"a"); 
      gtk_entry_set_text(entry, txt); 
    
      strcpy(entry_name, "phs_cell_b_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i,"b"); 
      gtk_entry_set_text(entry, txt); 
    
      strcpy(entry_name, "phs_cell_c_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i,"c"); 
      gtk_entry_set_text(entry, txt); 
    
      strcpy(entry_name, "phs_cell_alpha_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i,"alpha"); 
      gtk_entry_set_text(entry, txt); 
    
      strcpy(entry_name, "phs_cell_beta_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i,"beta"); 
      gtk_entry_set_text(entry, txt); 
    
      strcpy(entry_name, "phs_cell_gamma_entry_"); 
      tmp_name = entry_name + strlen(entry_name); 
      snprintf(tmp_name, 3, "%-d", i); 

      entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(window), entry_name));
      txt = get_text_for_phs_cell_chooser(i,"gamma"); 
      gtk_entry_set_text(entry, txt); 

      mol_count++;
      mol_with_cell = i;
    }
  } 
  display_none_cell_chooser_box(window, phs_cell_group); 

  /* For George Sheldrick: if there is only one molecule, then make
     that radio button be active (rather than None of the Above).  */
  if (mol_count == 1) { 
     strcpy(widget_name, "phs_cell_radiobutton_"); 
     tmp_name = widget_name + strlen(widget_name); 
     snprintf(tmp_name, 3, "%-d", mol_with_cell); 
     toggle_button = lookup_widget(GTK_WIDGET(window), widget_name);
     if (toggle_button) { 
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), TRUE);
     }
  }
  gtk_widget_show(window); 
  // force to the top? (even on Mac?) 
  // maybe the above show may not be needed any more

   gtk_window_set_transient_for(GTK_WINDOW(window), mw);
   gtk_window_present(GTK_WINDOW(window));
} 


/* Called back from the OK Button being pressed in the phs_info_box
   (which is the one that says that we need cell and symmetry from a
   pdb file).

   We need to go on to copy what make_and_draw_map() does, which
   includes map_fill_from_mtz(). 

   So actually, what we want is a wrapper function that calls
   molecule_class_info_t::map_fill_from_phs()

*/
int phs_pdb_cell_symm() {

   GtkWidget *widget; 

   if (1) {
     GtkWidget *file_filter_button;
     GtkWidget *sort_button;
     widget = create_phs_coordinates_filechooserdialog1();
     // we use th coords filter button, that should be ok
     // add_ccp4i_project_optionmenu(widget, COOT_COORDS_FILE_SELECTION);
     file_filter_button = add_filename_filter_button(widget, 
                                                     COOT_COORDS_FILE_SELECTION);
     sort_button = add_sort_button_fileselection(widget);
     push_the_buttons_on_fileselection(file_filter_button, sort_button, 
                                       widget);
   }

   gtk_widget_show(widget); 

   return 0; 

}


/* debugging function (that doesn't use C++) (the problem was that I
   was destroying the widget *then* useing the string (which of course
   had been deallocated) - nothing to do with c/c++ interface). */
void 
test_read_coords(const gchar *filename) { 

   printf("in the c test test_read_coords, filename is %s", filename);
}
