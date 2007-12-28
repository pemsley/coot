/* src/cmtz-interface.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include <iostream>

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

// We are not using NLS yet.
// #ifndef WINDOWS_MINGW
// #define ENABLE_NLS
// #endif
// #ifdef DATADIR
// #endif // DATADIR

#include <gtk/gtk.h>
#include <GL/glut.h> // for glutInit()


#include "interface.h"
#ifndef HAVE_SUPPORT_H
#define HAVE_SUPPORT_H
#include "support.h"
#endif /* HAVE_SUPPORT_H */

#include "read-phs.h"
#include "read-cif.h"

#include "clipper/ccp4/ccp4_mtz_io.h"
#include "cmtz-interface.hh"
#include "coot-utils.hh"


coot::mtz_column_types_info_t
coot::get_f_phi_columns(const std::string &filename) {

   // std::cout << "DEBUG:: getting f_phi_columns..." << std::endl;

   coot::mtz_column_types_info_t a;
   a.read_success = 0;
   a.selected_f_col      = 0;
   a.selected_phi_col    = -1; /* unset */
   a.selected_weight_col = 0;
   a.selected_refmac_fobs_col = 0;
   a.selected_refmac_sigfobs_col = 0;
   a.selected_refmac_r_free_col = 0;

   clipper::CCP4MTZfile f;
   short int is_mtz_file = 1; 
   // new try catch here
   
   try { 
      f.open_read(filename);
   }

   catch (...) {
      std::cout << "INFO:: not an mtz file: " << filename << std::endl;
      is_mtz_file = 0;
   } 

   if (is_mtz_file) { 
      std::vector<clipper::String> v = f.column_labels();
      // std::cout << "INFO:: found " << v.size() << " column labels in " << filename << "\n";
      if (v.size() > 1) { 
	 a.read_success = 1;
	 a.mtz_filename = filename; 
	 for (unsigned int i=0; i<v.size(); i++) {
	    // std::cout << i << " " << v[i] << "\n";
	    std::string label;
	    std::string type;
	    std::string::size_type ispace = v[i].find_last_of(" ");
	    if (ispace == std::string::npos) {
	       std::cout <<  "WARNING:: uninterprettable label \"" << v[i] << "\" of "
			 << filename << "\n";
	    } else {
	       label = v[i].substr(0, ispace);
	       type  = v[i].substr(ispace+1);
	       // std::cout << "Got label :" << label << ": and type :" << type << ":\n";
	       if (type == "F")
		  a.f_cols.push_back(coot::mtz_type_label(label, 'F'));
	       if (type == "G")
		  a.f_cols.push_back(coot::mtz_type_label(label, 'F'));
	       if (type == "L")
		  a.f_cols.push_back(coot::mtz_type_label(label, 'F'));
	       if (type == "Q")
		  a.sigf_cols.push_back(coot::mtz_type_label(label, 'Q'));
	       if (type == "P")
		  a.phi_cols.push_back(coot::mtz_type_label(label, 'P'));
	       if (type == "D")
		  a.d_cols.push_back(coot::mtz_type_label(label, 'D'));
	       if (type == "W")
		  a.weight_cols.push_back(coot::mtz_type_label(label, 'W'));
	       if (type == "I")
		  a.r_free_cols.push_back(coot::mtz_type_label(label, 'I'));
	    }
	 }
      }
   }
   return a;
}


/* used when the column label widget is being created   */
void
coot::setup_refmac_parameters(GtkWidget *window, 
			      const coot::mtz_column_types_info_t &col_labs) {
  
  int i;
  GtkWidget *fobs_option_menu    = lookup_widget(window, "refmac_fobs_optionmenu");
  GtkWidget *sigfobs_option_menu = lookup_widget(window, "refmac_sigfobs_optionmenu");
  GtkWidget *r_free_option_menu  = lookup_widget(window, "refmac_rfree_optionmenu");
  
  GtkWidget *fobs_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(fobs_option_menu));
  GtkWidget *sigfobs_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(sigfobs_option_menu));
  GtkWidget *r_free_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(r_free_option_menu));

  GtkWidget *menuitem;

  if (fobs_menu)
     gtk_widget_destroy(fobs_menu);
  if (sigfobs_menu)
     gtk_widget_destroy(sigfobs_menu);
  if (r_free_menu)
     gtk_widget_destroy(r_free_menu);

  fobs_menu = gtk_menu_new();
  sigfobs_menu = gtk_menu_new();
  r_free_menu = gtk_menu_new();


  /* Fobs */
  for (i=0; i<col_labs.f_cols.size(); i++) { 
     menuitem = make_menu_item((gchar *) col_labs.f_cols[i].column_label.c_str(),
			       GTK_SIGNAL_FUNC(refmac_f_button_select),
			       GINT_TO_POINTER(i));
     gtk_menu_append(GTK_MENU(fobs_menu), menuitem);
     gtk_widget_show(menuitem);
  } 
  /* Sig Fobs */
  for (i=0; i<col_labs.sigf_cols.size(); i++) {
     menuitem = make_menu_item( (gchar *) col_labs.sigf_cols[i].column_label.c_str(),
				GTK_SIGNAL_FUNC(refmac_sigf_button_select),
				GINT_TO_POINTER(i));
     gtk_menu_append(GTK_MENU(sigfobs_menu), menuitem);
     gtk_widget_show(menuitem);
  }


  /* R free */

  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));

  if (col_labs.r_free_cols.size() == 0) 
    save_f_phi_columns->selected_refmac_r_free_col = -1; /* magic -1 */
  
				/* see on_column_label_ok_button_clicked in callbacks.c */
  for (i=0; i<col_labs.r_free_cols.size(); i++) {
     menuitem = make_menu_item((gchar *) col_labs.r_free_cols[i].column_label.c_str(),
			       GTK_SIGNAL_FUNC(refmac_r_free_button_select),
			       GINT_TO_POINTER(i));
     gtk_menu_append(GTK_MENU(r_free_menu), menuitem);
     gtk_widget_show(menuitem);
  }

  /* Link the menus to the optionmenus */
  gtk_option_menu_set_menu(GTK_OPTION_MENU(fobs_option_menu), fobs_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(sigfobs_option_menu), sigfobs_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(r_free_option_menu), r_free_menu);

  gtk_widget_show(fobs_menu);
  gtk_widget_show(sigfobs_menu);
  gtk_widget_show(r_free_menu);
  
}


std::vector<std::string> coot::get_f_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;

   return v;
}

std::vector<std::string> coot::get_phi_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;

   return v;
}

std::vector<std::string> coot::get_weight_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;

   return v;
}

std::vector<std::string> coot::get_d_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;

   return v;
}


void 
f_button_select(GtkWidget *item, GtkPositionType pos) { 
   
   GtkWidget *window;
   GtkWidget *checkbutton;
   std::string lab;
   short int make_diff_map_flag = 0;
   
   /* If this was an anomalous label, we want a difference map... */
   
   window = lookup_widget(GTK_WIDGET(item), "column_label_window");
   coot::mtz_column_types_info_t *save_f_phi_columns
      = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
   
   save_f_phi_columns->selected_f_col = pos;
//    std::cout << "debug:: pos " << pos << " f_cols size "
// 	     << save_f_phi_columns->f_cols.size() << std::endl;
   if (pos >= save_f_phi_columns->f_cols.size()) {
      //       printf("%d was an anomalous label\n", pos);
      make_diff_map_flag = 1;
      lab = save_f_phi_columns->d_cols[pos-save_f_phi_columns->f_cols.size()].column_label;
   } else {
      lab = save_f_phi_columns->f_cols[pos].column_label;
   }
   
   /* also add code that checks to see if the column label begins with
      "DEL" and if so, it changes the difference-map? checkbutton on
      this window to be active (default is inactive).  */

   std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(lab);
   
   if (p.second.length() > 2) {
      // std::cout << "DEBUG DEL test :" << p.second.substr(0,3) << ":\n";
      if ( p.second.substr(0,3) == "DEL") { 
	 make_diff_map_flag = 1;
      }
   }
   
   if (p.second.length() > 3) {
      if (p.second.substr(0,4) == "FOFC") {
	 make_diff_map_flag = 1; 
      }
   }
   
   if (make_diff_map_flag) {
      checkbutton = lookup_widget(window, "difference_map_checkbutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }
}


GtkWidget *make_menu_item( gchar         *name,
                           GtkSignalFunc  callback,
                           gpointer       data )
{
    GtkWidget *item;
  
    item = gtk_menu_item_new_with_label (name);
    gtk_signal_connect (GTK_OBJECT (item), "activate",
                        callback, data);
    gtk_widget_show (item);

    return(item);
}


void
phase_button_select(GtkWidget *item, GtkPositionType pos) { 
   
   /*     printf("setting phase position %d\n", pos);   */
   GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));

  save_f_phi_columns->selected_phi_col = pos;
}

void
weight_button_select(GtkWidget *item, GtkPositionType pos) { 
   
/*    printf("setting weight  position %d\n", pos);  */
   GtkWidget *window = lookup_widget(item, "column_label_window");
   coot::mtz_column_types_info_t *save_f_phi_columns
      = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
   save_f_phi_columns->selected_weight_col = pos;
}


void 
refmac_f_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting refmac f obs position %d\n", pos); 
   GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_fobs_col = pos;
}

void 
refmac_sigf_button_select(GtkWidget *item, GtkPositionType pos) { 
   
   printf("setting refmac sigf position %d\n", pos);  
   GtkWidget *window = lookup_widget(item, "column_label_window");
   coot::mtz_column_types_info_t *save_f_phi_columns
      = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
   save_f_phi_columns->selected_refmac_sigfobs_col = pos;
}

void 
refmac_r_free_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting r free position %d\n", pos); 
   GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_r_free_col = pos;
}


void
coot::fill_f_optionmenu(GtkWidget *optionmenu_f, short int is_expert_mode_flag) { 

   /* Tinker with optionmenu1, the selection of Fs */

   /* create a menu for the optionmenu button.  The various column
    labels will be added to this menu as menuitems*/
  GtkWidget *optionmenu1_menu;
  GtkWidget *menuitem;
  int i;
  int n_lab = 0;

  optionmenu1_menu = gtk_menu_new();
   
  GtkWidget *window = lookup_widget(optionmenu_f, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));

//    std::cout << "DEBUG:: (get) user data save_f_phi_columns pointer: " << save_f_phi_columns
// 	     << std::endl;

   if (! save_f_phi_columns) {
     std::cout << "ERROR:: null save_f_phi_columns in fill_f_optionmenu\n";
     return; 

  } else { 

      for (i=0; i< save_f_phi_columns->f_cols.size(); i++) {
// 	 std::cout << "Making f menu item for  " << i << "/"
// 		   << save_f_phi_columns->f_cols.size() << " " 
// 		   << save_f_phi_columns->f_cols[i].column_label << "\n";
	 menuitem = make_menu_item((gchar *) save_f_phi_columns->f_cols[i].column_label.c_str(),
				   GTK_SIGNAL_FUNC(f_button_select),
				   GINT_TO_POINTER(i));
	 gtk_menu_append(GTK_MENU(optionmenu1_menu), menuitem);
	 gtk_widget_show(menuitem);
	 // If label was xxxxFWT then make it be the active menu item
	 int l = save_f_phi_columns->f_cols[i].column_label.length();
	 if (l > 3) {
	    std::string last_bit = save_f_phi_columns->f_cols[i].column_label.substr(l-4, 4);
	    // std::cout << "DEBUG:: last_bit :" << last_bit << ":\n";
	    if (last_bit == "/FWT") { 
	       /* was FWT */
	       gtk_menu_set_active(GTK_MENU(optionmenu1_menu), i);
	       save_f_phi_columns->selected_f_col = i;
	    }
	 }
	 n_lab++;
      }

      if (is_expert_mode_flag) { 
	 for (i=0; i< save_f_phi_columns->d_cols.size(); i++) { 
	    menuitem = make_menu_item((gchar *) save_f_phi_columns->d_cols[i].column_label.c_str(),
				      GTK_SIGNAL_FUNC(f_button_select),
				      GINT_TO_POINTER(i+n_lab));
	    gtk_menu_append(GTK_MENU(optionmenu1_menu), menuitem);
	    gtk_widget_show(menuitem);
	 } 
      }
      /* Link the menu  to the optionmenu widget */
      gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu_f), optionmenu1_menu);
   }
}

////////////////////////////////////////////////////////////////////////////
//              column_selector_using_cmtz()
////////////////////////////////////////////////////////////////////////////
// 
GtkWidget *
coot::column_selector_using_cmtz(const std::string &filename) { 
   
   int i;
   GtkWidget *column_label_window;
   GtkWidget *optionmenu2_menu, *optionmenu3_menu;
   GtkWidget *menuitem;
   
   GtkWidget *optionmenu_f, *optionmenu_phi, *optionmenu_weight;
   GtkCheckButton *check_weights;
   int is_phs = 0;
   int is_cif = 0;
   int is_cns_data = 0;
   int is_expert_mode_flag = 0;

   // We use this indirection because f_phi_columns gets attached to
   // the column chooser dialog
   coot::mtz_column_types_info_t *f_phi_columns = new coot::mtz_column_types_info_t;
   *f_phi_columns = coot::get_f_phi_columns(filename);
   f_phi_columns->mtz_filename = filename;

//    std::cout << "DEBUG:: in column_selector_using_cmtz got read success of "
// 	     << f_phi_columns->read_success << std::endl;
      
   if (f_phi_columns->read_success == 0 ) { /*  not a valid mtz file */
      std::cout << "INFO:: data file " << filename << " is not a valid mtz file\n";
      is_phs = try_read_phs_file(filename.c_str()); /* Try reading the data file 
						   as an XtalView .phs file */ 
      if (is_phs == 0) { 
	 is_cif = try_read_cif_file(filename.c_str()); 

	 // This no longer makes sense given that
	 // try_read_cns_data_file takes an imol argument also now.
// 	 if (is_cif == 0) 
// 	    is_cns_data = try_read_cns_data_file(filename.c_str());
      }

      return 0;
   } 

   /* Else filename was OK */


/*    printf("The F columns: \n"); */
/*    for (i=0; i< a.n_f_cols; i++) {  */
/*       printf("%s\n", a.f_cols[i].column_label); */
/*    } */

/*    printf("The P columns: \n"); */
/*    for (i=0; i< a.n_phi_cols; i++) {  */
/*       printf("%s\n", a.phi_cols[i].column_label); */
/*    } */


/* Recall that save_f_phi_columns is now attached to this widget */

   /* Stuff a pointer to mtz info into the dialog: */
   column_label_window = create_column_label_window();
   gtk_object_set_user_data(GTK_OBJECT(column_label_window), f_phi_columns);

   /* The default column labels are at the top of the list.  The
      selcted_{f,phi}_cols get changed by the menubutton function
      (callbacks I suppose, but they are not listed there, they are in
      this file). */

   optionmenu_f = lookup_widget(column_label_window, "optionmenu1"); 
   coot::fill_f_optionmenu(optionmenu_f, is_expert_mode_flag); /* flag: 0 */


   /* And now do the same for the P (phase) columns */

   optionmenu2_menu = gtk_menu_new(); 
   
   /* set the default to the top for phases (if we have phase columns) */
   if (f_phi_columns->phi_cols.size() > 0) 
     f_phi_columns->selected_phi_col    = 0; /* set to top */

   for (i=0; i< f_phi_columns->phi_cols.size(); i++) {
      menuitem = make_menu_item((gchar *)f_phi_columns->phi_cols[i].column_label.c_str(),
				GTK_SIGNAL_FUNC(phase_button_select),
				GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(optionmenu2_menu), menuitem);
      gtk_widget_show(menuitem);
      // If label was xxxxFWT then make it be the active menu item
      int l = f_phi_columns->phi_cols[i].column_label.length();
      if (l > 5) { 
	 if  (f_phi_columns->phi_cols[i].column_label.substr(l-5, l) == "/PHWT") {
				/* was PHWT */
	    gtk_menu_set_active(GTK_MENU(optionmenu2_menu), i);
	    f_phi_columns->selected_phi_col = i;
	 }
      }
   }

   
   optionmenu_phi = lookup_widget(column_label_window, "optionmenu2"); 

   /* Link the menu  to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu_phi), optionmenu2_menu);


   /* And now do the same for the W (weight) columns */


   optionmenu3_menu = gtk_menu_new(); 
   
   for (i=0; i< f_phi_columns->weight_cols.size(); i++) {
      menuitem = make_menu_item( (gchar *) f_phi_columns->weight_cols[i].column_label.c_str(),
				 GTK_SIGNAL_FUNC(weight_button_select),
				 GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(optionmenu3_menu), menuitem);
      gtk_widget_show(menuitem);
   }
   
   optionmenu_weight = lookup_widget(column_label_window, 
				    "optionmenu3"); 
 
   /* Link the menu  to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu_weight), 
			    optionmenu3_menu);


   /* By default, we want the use weights checkbutton to be off */

   check_weights = GTK_CHECK_BUTTON(lookup_widget(column_label_window, 
						  "use_weights_checkbutton"));
   
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_weights), FALSE);


   /* New addition: the refmac buttons  */
   coot::setup_refmac_parameters(column_label_window, *f_phi_columns);

   return column_label_window;
}



