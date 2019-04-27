/* src/cmtz-interface.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005 by The University of York
 * Copyright 2013 by Medical Research Council
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

#include "compat/coot-sysdep.h"


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

#include <clipper/ccp4/ccp4_mtz_io.h>

#include "interface.h"
#ifndef HAVE_SUPPORT_H
#define HAVE_SUPPORT_H
#include "support.h"
#endif /* HAVE_SUPPORT_H */

#include "read-phs.h"
#include "read-cif.h"

#include "cmtz-interface.hh"
#include "cmtz-interface-gui.hh"
#include "utils/coot-utils.hh"

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "c-interface-refmac.h"
#include "graphics-info.h"


coot::mtz_column_types_info_t
coot::get_mtz_columns(const std::string &filename) {

   coot::mtz_column_types_info_t a;
   a.read_success = 0;
   a.selected_f_col      = 0;
   a.selected_phi_col    = -1; /* unset */
   a.selected_weight_col = 0;
   a.selected_refmac_fobs_col = 0;
   a.selected_refmac_sigfobs_col = 0;
   a.selected_refmac_r_free_col = 0;
   a.selected_refmac_phi_col = -1; /* unset */
   a.selected_refmac_fom_col = 0;
   a.selected_refmac_hla_col = -1; /* unset */
   a.selected_refmac_hlb_col = 0;
   a.selected_refmac_hlc_col = 0;
   a.selected_refmac_hld_col = 0;
   // new ones for SAD and twin refinement
   // I, I+/-, F+/- and appropriate sigs
   a.selected_refmac_fp_col = 0;
   a.selected_refmac_sigfp_col = 0;
   a.selected_refmac_fm_col = 0;
   a.selected_refmac_sigfm_col = 0;
   a.selected_refmac_iobs_col = -1; /* unset */
   a.selected_refmac_sigiobs_col = 0;
   a.selected_refmac_ip_col = 0;
   a.selected_refmac_sigip_col = 0;
   a.selected_refmac_im_col = 0;
   a.selected_refmac_sigim_col = 0;

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
	       //std::cout << "Got label :" << label << ": and type :" << type << ":\n";
	       if (type == "F")
		 a.f_cols.push_back(coot::mtz_type_label(label, 'F', i));
	       if (type == "G")
		 a.fpm_cols.push_back(coot::mtz_type_label(label, 'G', i));
	       if (type == "L")
		 a.sigfpm_cols.push_back(coot::mtz_type_label(label, 'L', i));
	       if (type == "Q")
		 a.sigf_cols.push_back(coot::mtz_type_label(label, 'Q', i));
	       if (type == "P")
		 a.phi_cols.push_back(coot::mtz_type_label(label, 'P', i));
	       if (type == "D")
		 a.d_cols.push_back(coot::mtz_type_label(label, 'D', i));
	       if (type == "W")
		 a.weight_cols.push_back(coot::mtz_type_label(label, 'W', i));
	       if (type == "I")
		 a.r_free_cols.push_back(coot::mtz_type_label(label, 'I', i));
	       if (type == "A")
		 a.hl_cols.push_back(coot::mtz_type_label(label, 'A', i));
	       if (type == "J")
		 a.i_cols.push_back(coot::mtz_type_label(label, 'J', i));
	       // for completeness; not used yet
	       if (type == "K")
		 a.ipm_cols.push_back(coot::mtz_type_label(label, 'K', i));
	       if (type == "M")
		 a.sigipm_cols.push_back(coot::mtz_type_label(label, 'M', i));
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
  
  unsigned int i;
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


/* used when the column label widget is being created   */
void
coot::setup_refmac_parameters_from_file(GtkWidget *window) {
  
  GtkWidget *option_menu;
  std::string filename;
  GtkWidget *file_mtz_radiobutton = lookup_widget(window, "run_refmac_mtz_file_radiobutton");
  if (GTK_TOGGLE_BUTTON(file_mtz_radiobutton)->active) {
    // twin/mtz filename given
    GtkWidget *file_mtz_label = lookup_widget(window, "run_refmac_mtz_file_label");
    // in twin we use the label as a dummy widget
    option_menu = file_mtz_label;

    const gchar *mtz_filename = gtk_label_get_text(GTK_LABEL(file_mtz_label));
    filename = mtz_filename;
    if (coot::file_exists(filename)) {
      graphics_info_t::saved_refmac_file_filename = filename.c_str();
    } else {
      filename = "";
    }
  } else {
    // not twin/mtz filename but map
    // first get the filename from the currently selected mtz/map
    option_menu = lookup_widget(window, "run_refmac_map_optionmenu");
    GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
    GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));
    if (active_item == 0) {
      add_status_bar_text("No map has associated Refmac Parameters - no REFMAC!");
    } else { 
      int imol_window = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(active_item)));
      if (imol_window < 0) {
	std::cout << "No map data selected for refmac\n";
      } else { 
      
	int imol_map_refmac = imol_window;
	if (!is_valid_map_molecule(imol_map_refmac)) {
	  std::string s = "Invalid molecule number: ";
	  s += graphics_info_t::int_to_string(imol_map_refmac);
	  std::cout << s << std::endl;
	  graphics_info_t g;
	  g.add_status_bar_text(s);
	} else {
	  // just check for refmac mtz file now
	  if (graphics_info_t::molecules[imol_map_refmac].Refmac_mtz_filename().size() > 0) {
	    filename = graphics_info_t::molecules[imol_map_refmac].Refmac_mtz_filename();
	  } else {
	    std::cout << "No valid mtz file" <<std::endl;
	  }
	}
      }
    }
  }

  if (filename.size() > 0) {
   // get the column labels
   coot::mtz_column_types_info_t *f_phi_columns = new coot::mtz_column_types_info_t;
   *f_phi_columns = coot::get_mtz_columns(filename);
   f_phi_columns->mtz_filename = filename;
   const coot::mtz_column_types_info_t &col_labs = *f_phi_columns;

   /* Stuff a pointer to mtz info into the dialog: */
   GtkWidget *refmac_dialog = lookup_widget(window, "run_refmac_dialog");
   gtk_object_set_user_data(GTK_OBJECT(refmac_dialog), f_phi_columns);

   unsigned int i;

   GtkWidget *fobs_option_menu    = lookup_widget(window, "refmac_dialog_fobs_optionmenu");
   GtkWidget *fpm_option_menu     = lookup_widget(window, "refmac_dialog_fpm_optionmenu");
   GtkWidget *fiobs_option_menu   = lookup_widget(window, "refmac_dialog_fiobs_optionmenu");
   GtkWidget *ipm_option_menu     = lookup_widget(window, "refmac_dialog_ipm_optionmenu");
   GtkWidget *r_free_option_menu  = lookup_widget(window, "refmac_dialog_rfree_optionmenu");
   GtkWidget *phases_option_menu  = lookup_widget(window, "refmac_dialog_phases_optionmenu");
   GtkWidget *fom_option_menu     = lookup_widget(window, "refmac_dialog_fom_optionmenu");
   GtkWidget *hl_option_menu      = lookup_widget(window, "refmac_dialog_hl_optionmenu");

   GtkWidget *fobs_menu   = gtk_option_menu_get_menu(GTK_OPTION_MENU(fobs_option_menu));
   GtkWidget *fpm_menu    = gtk_option_menu_get_menu(GTK_OPTION_MENU(fpm_option_menu));
   /* we dont have a separate I optionsmenu, but append all Fs and Is into one combined one */
   GtkWidget *fiobs_menu  = gtk_option_menu_get_menu(GTK_OPTION_MENU(fiobs_option_menu));
   GtkWidget *ipm_menu    = gtk_option_menu_get_menu(GTK_OPTION_MENU(ipm_option_menu));
   GtkWidget *r_free_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(r_free_option_menu));
   GtkWidget *phases_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(phases_option_menu));
   GtkWidget *fom_menu    = gtk_option_menu_get_menu(GTK_OPTION_MENU(fom_option_menu));
   GtkWidget *hl_menu     = gtk_option_menu_get_menu(GTK_OPTION_MENU(hl_option_menu));

   GtkWidget *menuitem;
   GtkWidget *menuitem2;

   if (fobs_menu)
      gtk_widget_destroy(fobs_menu);
   if (fpm_menu)
      gtk_widget_destroy(fpm_menu);
   if (fiobs_menu)
      gtk_widget_destroy(fiobs_menu);
   if (ipm_menu)
      gtk_widget_destroy(ipm_menu);
   if (r_free_menu)
      gtk_widget_destroy(r_free_menu);
   if (phases_menu)
      gtk_widget_destroy(phases_menu);
   if (fom_menu)
      gtk_widget_destroy(fom_menu);
   if (hl_menu)
      gtk_widget_destroy(hl_menu);

   fobs_menu = gtk_menu_new();
   fpm_menu = gtk_menu_new();
   fiobs_menu = gtk_menu_new();
   ipm_menu = gtk_menu_new();
   r_free_menu = gtk_menu_new();
   phases_menu = gtk_menu_new();
   fom_menu = gtk_menu_new();
   hl_menu = gtk_menu_new();


  /* Fobs and Sig Fobs*/
  /* we always assume following pairs of Fobs and sigFobs */
  int fobs_pos;
  int sigfobs_pos;
  int islash;
  std::string label;
  std::string tmp;
  coot::mtz_column_types_info_t a;
  for (unsigned ii=0; ii<col_labs.f_cols.size(); ii++) {
     fobs_pos = col_labs.f_cols[ii].column_position;
     /* Sig Fobs */
     for (unsigned int j=0; j<col_labs.sigf_cols.size(); j++) {
       sigfobs_pos = col_labs.sigf_cols[j].column_position;
       if ((fobs_pos + 1) == sigfobs_pos) {
	 // matching f sigf pair
	 tmp = col_labs.sigf_cols[j].column_label;
	 islash = tmp.find_last_of("/");
	 label = col_labs.f_cols[ii].column_label;
	 label += "  ";
	 label += tmp.substr(islash + 1);
	 menuitem  = make_menu_item((gchar *) label.c_str(),
				   GTK_SIGNAL_FUNC(refmac_dialog_f_button_select),
				   GINT_TO_POINTER(ii));
	 menuitem2 = make_menu_item((gchar *) label.c_str(),
				   GTK_SIGNAL_FUNC(refmac_dialog_f_button_select),
				   GINT_TO_POINTER(ii));
	 gtk_menu_append(GTK_MENU(fobs_menu), menuitem);
	 // append to the combined one too
	 gtk_menu_append(GTK_MENU(fiobs_menu), menuitem2);
	 gtk_widget_show(menuitem);
       }
     }
  }


  /* F+, sigF+ and F-, sigF-*/
  /* assume in pair of pairs (F+,sF+,F-,sF-) */
  /* NB: we take F+,F-,sF+,sF- into account too */
  int fp_pos;
  int fm_pos;
  int sigfp_pos;
  int sigfm_pos;
  //  coot::mtz_column_types_info_t a;
  // why do I use mmdb functions here?
  int good_no_of_fpm    = mmdb::mod(col_labs.fpm_cols.size(), 2);
  int good_no_of_sigfpm = mmdb::mod(col_labs.sigfpm_cols.size(), 2);
  if (good_no_of_fpm && good_no_of_sigfpm) {
    std::cout << "WARNING:: inconsistent number of F+/F- and or sigF+/sigF-, i.e. not multiple of 2.\n" << std::endl;
    std::cout << "Detection of F+/F- and associated sigmas may be screwed!" << std::endl;
  }
  for (unsigned icol=0; icol<col_labs.fpm_cols.size(); icol++) {
     fp_pos = col_labs.fpm_cols[icol].column_position;
     if (icol+1 < col_labs.fpm_cols.size()) {
       fm_pos = col_labs.fpm_cols[icol+1].column_position;
       /* Sig F+/- */
       for (unsigned int j=0; j<col_labs.sigfpm_cols.size(); j++) {
	 sigfp_pos = col_labs.sigfpm_cols[j].column_position;
	 if (j+1 < col_labs.sigfpm_cols.size()) {
	   sigfm_pos = col_labs.sigfpm_cols[icol+1].column_position;
	   if ((fm_pos == fp_pos+2 && sigfp_pos == fp_pos+1 && sigfm_pos == fm_pos+1) ||
	       (fm_pos == fp_pos+1 && sigfp_pos == fp_pos+2 && sigfm_pos == fm_pos+2)) {
	     // matching f/f- sigf+/- pairs
	     label = col_labs.fpm_cols[icol].column_label;
	     label += "  ";
	     tmp = col_labs.fpm_cols[icol+1].column_label;
	     islash = tmp.find_last_of("/");
	     label += tmp.substr(islash + 1);
	     label += "  ";
	     tmp = col_labs.sigfpm_cols[j].column_label;
	     islash = tmp.find_last_of("/");
	     label += tmp.substr(islash + 1);
	     label += "  ";
	     tmp = col_labs.sigfpm_cols[j+1].column_label;
	     islash = tmp.find_last_of("/");
	     label += tmp.substr(islash + 1);
	     menuitem = make_menu_item((gchar *) label.c_str(),
				       GTK_SIGNAL_FUNC(refmac_dialog_fpm_button_select),
				       GINT_TO_POINTER(icol));
	     gtk_menu_append(GTK_MENU(fpm_menu), menuitem);
	     gtk_widget_show(menuitem);
	   }
	 }
       }
     }
  }


  /* Iobs and Sig Iobs*/
  /* we always asume following pairs of Iobs and sigIobs */
  /* same as for F/sigF BUT we dont have an extra column for SigI
     so we get it from the sigF column */
  int iobs_pos = 0;
  int sigiobs_pos;
  for (i=0; i<col_labs.i_cols.size(); i++) {
     fobs_pos = col_labs.i_cols[i].column_position;
     /* Sig Iobs */
     for (unsigned int j=0; j<col_labs.sigf_cols.size(); j++) {
       sigiobs_pos = col_labs.sigf_cols[j].column_position;
       if ((iobs_pos + 1) == sigiobs_pos) {
	 // matching i sigi pair
	 tmp = col_labs.sigf_cols[j].column_label;
	 islash = tmp.find_last_of("/");
	 label = col_labs.i_cols[i].column_label;
	 label += "  ";
	 label += tmp.substr(islash + 1);
	 menuitem = make_menu_item((gchar *) label.c_str(),
				   GTK_SIGNAL_FUNC(refmac_dialog_i_button_select),
				   GINT_TO_POINTER(i));
	 // only append to combined menu
	 gtk_menu_append(GTK_MENU(fiobs_menu), menuitem);
	 gtk_widget_show(menuitem);
       }
     }
  }


  /* same as for F+/F-; Not used currently */
  /* I+, sigI+ and I-, sigI-*/
  /* assume in pair of pairs (F+,sF+,F-,sF-) */
  /* NB: we take F+,F-,sF+,sF- into account too */
  int ip_pos;
  int im_pos;
  int sigip_pos;
  int sigim_pos;
  int good_no_of_ipm    = mmdb::mod(col_labs.ipm_cols.size(),    2);
  int good_no_of_sigipm = mmdb::mod(col_labs.sigipm_cols.size(), 2);
  if (good_no_of_ipm && good_no_of_sigipm) {
     std::cout << "WARNING:: inconsistent number of I+/I- and or sigI+/sigI-, i.e. not multiple of 2.\n"
	       << std::endl;
     std::cout << "Detection of I+/I- and associated sigmas may be screwed!" << std::endl;
  }
  for (i=0; i<col_labs.ipm_cols.size(); i++) {
     ip_pos = col_labs.ipm_cols[i].column_position;
     if (i+1 < col_labs.ipm_cols.size()) {
       im_pos = col_labs.ipm_cols[i+1].column_position;
       /* Sig I+/- */
       for (unsigned int j=0; j<col_labs.sigipm_cols.size(); j++) {
	 sigip_pos = col_labs.sigipm_cols[j].column_position;
	 if (j+1 < col_labs.sigipm_cols.size()) {
	   sigim_pos = col_labs.sigipm_cols[i+1].column_position;
	   if ((im_pos == ip_pos+2 && sigip_pos == ip_pos+1 && sigim_pos == im_pos+1) ||
	       (im_pos == ip_pos+1 && sigip_pos == ip_pos+2 && sigim_pos == im_pos+2)) {
	     // matching i+/i- sigi+/- pairs
	     label = col_labs.ipm_cols[i].column_label;
	     label += "  ";
	     tmp = col_labs.ipm_cols[i+1].column_label;
	     islash = tmp.find_last_of("/");
	     label += tmp.substr(islash + 1);
	     label += "  ";
	     tmp = col_labs.sigipm_cols[j].column_label;
	     islash = tmp.find_last_of("/");
	     label += tmp.substr(islash + 1);
	     label += "  ";
	     tmp = col_labs.sigipm_cols[j+1].column_label;
	     islash = tmp.find_last_of("/");
	     label += tmp.substr(islash + 1);
	     menuitem = make_menu_item((gchar *) label.c_str(),
				       GTK_SIGNAL_FUNC(refmac_dialog_ipm_button_select),
				       GINT_TO_POINTER(i));
	     gtk_menu_append(GTK_MENU(ipm_menu), menuitem);
	     gtk_widget_show(menuitem);
	   }
	 }
       }
     }
  }


   /* R free */

   if (col_labs.r_free_cols.size() == 0) {
     f_phi_columns->selected_refmac_r_free_col = -1; /* magic -1 */
   }

   /* see on_column_label_ok_button_clicked in callbacks.c */
   for (i=0; i<col_labs.r_free_cols.size(); i++) {
      menuitem = make_menu_item((gchar *) col_labs.r_free_cols[i].column_label.c_str(),
				GTK_SIGNAL_FUNC(refmac_dialog_r_free_button_select),
				GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(r_free_menu), menuitem);
      gtk_widget_show(menuitem);
   }

   /* Phases and FOM */
  for (i=0; i<col_labs.phi_cols.size(); i++) {
    menuitem = make_menu_item((gchar *) col_labs.phi_cols[i].column_label.c_str(),
			      GTK_SIGNAL_FUNC(refmac_dialog_phases_button_select),
			      GINT_TO_POINTER(i));
    gtk_menu_append(GTK_MENU(phases_menu), menuitem);
    gtk_widget_show(menuitem);
  }

  for (i=0; i<col_labs.weight_cols.size(); i++) {
    menuitem = make_menu_item((gchar *) col_labs.weight_cols[i].column_label.c_str(),
			      GTK_SIGNAL_FUNC(refmac_dialog_fom_button_select),
			      GINT_TO_POINTER(i));
    gtk_menu_append(GTK_MENU(fom_menu), menuitem);
    gtk_widget_show(menuitem);
  }

   /* HL coefficients */
  int hl_pos;
  int hl_next_pos;
  int good_no_of_hls = mmdb::mod(col_labs.hl_cols.size(), 4);
  if (good_no_of_hls) {
    std::cout << "WARNING:: inconsistent number of HL coefficients, i.e. not multiple of 4.\nDetection of HL sets may be screwed!" << std::endl;
  }
  for (i=0; i<col_labs.hl_cols.size(); i++) {
    int hl_err = 0;
    // find the first one and add the next 3 and increment i
    hl_pos = col_labs.hl_cols[i].column_position;
    label = col_labs.hl_cols[i].column_label;
    if (i+1 < col_labs.hl_cols.size()) {
       for (unsigned int j=i+1; j<i+4; j++) {
	// check if next position is HL too
	hl_next_pos = col_labs.hl_cols[j].column_position;
	if (hl_next_pos == (hl_pos + 1)) {
	  tmp = col_labs.hl_cols[j].column_label;
	  islash = tmp.find_last_of("/");
	  label += "  ";
	  label += tmp.substr(islash + 1);
	  hl_pos += 1;
	} else {
	  std::cout<<"WARNING:: inconsistent HL coefficient set! Only " << j-i << " coefficients found!" << std::endl;
	  hl_err = 1;
	  break;
	}
      }
    } else {
      hl_err = 1;
    }
    if (! hl_err) {
      menuitem = make_menu_item((gchar *) label.c_str(),
				GTK_SIGNAL_FUNC(refmac_dialog_hl_button_select),
				GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(hl_menu), menuitem);
      gtk_widget_show(menuitem);
      i += 3;
    }
  }

  update_refmac_column_labels_frame(option_menu,
				    fobs_menu, fiobs_menu, fpm_menu,
				    r_free_menu, 
				    phases_menu, fom_menu, hl_menu);
  
  /* Link the menus to the optionmenus */
  gtk_option_menu_set_menu(GTK_OPTION_MENU(fobs_option_menu), fobs_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(fiobs_option_menu), fiobs_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(fpm_option_menu), fpm_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(r_free_option_menu), r_free_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(phases_option_menu), phases_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(fom_option_menu), fom_menu);
  gtk_option_menu_set_menu(GTK_OPTION_MENU(hl_option_menu), hl_menu);

  //graphics_info_t::update_refmac_column_labels_frame(option_menu);

  //std::cout << "BL DEBUG:: selected fobs col " << save_f_phi_columns->selected_refmac_fobs_col<<std::endl;
  //std::cout << "BL DEBUG:: selected sigf col " << save_f_phi_columns->selected_refmac_sigfobs_col<<std::endl;
  //std::cout << "BL DEBUG:: selected phi col " << save_f_phi_columns->selected_refmac_phi_col<<std::endl;
  //std::cout << "BL DEBUG:: selected fom col " << save_f_phi_columns->selected_refmac_fom_col<<std::endl;
  //std::cout << "BL DEBUG:: selected hla col " << save_f_phi_columns->selected_refmac_hla_col<<std::endl;

  }
}


std::vector<std::string>
coot::get_f_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;
   mtz_column_types_info_t ti = get_mtz_columns(mtz_file_name);
   for (unsigned int i=0; i<ti.f_cols.size(); i++) { 
      v.push_back(ti.f_cols[i].column_label);
   }
   return v;
}


std::vector<std::string>
coot::get_sigf_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;
   mtz_column_types_info_t ti = get_mtz_columns(mtz_file_name);
   for (unsigned int i=0; i<ti.sigf_cols.size(); i++) {
      v.push_back(ti.sigf_cols[i].column_label);
   }
   return v;
}

std::vector<std::string>
coot::get_r_free_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;
   mtz_column_types_info_t ti = get_mtz_columns(mtz_file_name);
   for (unsigned int i=0; i<ti.r_free_cols.size(); i++) { 
      v.push_back(ti.r_free_cols[i].column_label);
   }
   return v;
} 



std::vector<std::string>
coot::get_phi_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;
   mtz_column_types_info_t ti = get_mtz_columns(mtz_file_name);
   for (unsigned int i=0; i<ti.phi_cols.size(); i++) { 
      v.push_back(ti.phi_cols[i].column_label);
   }
   return v;
}

std::vector<std::string>
coot::get_weight_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;
   mtz_column_types_info_t ti = get_mtz_columns(mtz_file_name);
   for (unsigned int i=0; i<ti.weight_cols.size(); i++) { 
      v.push_back(ti.weight_cols[i].column_label);
   }
   return v;
}

std::vector<std::string>
coot::get_d_cols(const std::string &mtz_file_name) {

   std::vector<std::string> v;
   mtz_column_types_info_t ti = get_mtz_columns(mtz_file_name);
   for (unsigned int i=0; i<ti.d_cols.size(); i++) { 
      v.push_back(ti.d_cols[i].column_label);
   }
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
   if (pos >= int(save_f_phi_columns->f_cols.size())) {
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

// copy of above for newer interface
void 
refmac_dialog_f_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting refmac f obs position %d\n", pos); 
  GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_fobs_col = pos;
  // get the corresponding Sig Fobs
  int sigfobs_pos;
  int fobs_pos = save_f_phi_columns->f_cols[pos].column_position;
  for (unsigned int i=0; i<save_f_phi_columns->sigf_cols.size(); i++) {
    sigfobs_pos = save_f_phi_columns->sigf_cols[i].column_position;
    if (sigfobs_pos == (fobs_pos + 1)) {
      save_f_phi_columns->selected_refmac_sigfobs_col = i;
    }
  }
}

// obsolete?
//void 
//refmac_dialog_sigf_button_select(GtkWidget *item, GtkPositionType pos) { 
//   
//   printf("setting refmac sigf position %d\n", pos);  
//   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
//   coot::mtz_column_types_info_t *save_f_phi_columns
//      = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
//   save_f_phi_columns->selected_refmac_sigfobs_col = pos;
//}

void 
refmac_dialog_fpm_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting refmac f+/- obs position %d\n", pos); 
  GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_fp_col = pos;
  save_f_phi_columns->selected_refmac_fm_col = pos+1;
  // get the corresponding F- and SigF+/F-
  int sigfp_pos;
  int sigfm_pos;
  int fp_pos = save_f_phi_columns->fpm_cols[pos].column_position;
  int fm_pos = save_f_phi_columns->fpm_cols[pos+1].column_position;
  for (unsigned int i=0; i<save_f_phi_columns->sigfpm_cols.size()-1; i++) {
    sigfp_pos = save_f_phi_columns->sigfpm_cols[i].column_position;
    sigfm_pos = save_f_phi_columns->sigfpm_cols[i+1].column_position;
    if ((fm_pos == fp_pos+2 && sigfp_pos == fp_pos+1 && sigfm_pos == fm_pos+1) ||
	(fm_pos == fp_pos+1 && sigfp_pos == fp_pos+2 && sigfm_pos == fm_pos+2)) {
      save_f_phi_columns->selected_refmac_sigfp_col = i;
      save_f_phi_columns->selected_refmac_sigfm_col = i+1;
    }
  }
}

// as for F but I
void 
refmac_dialog_i_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting refmac i obs position %d\n", pos); 
  GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_iobs_col = pos;
  // get the corresponding Sig Fobs
  int sigiobs_pos;
  int iobs_pos = save_f_phi_columns->i_cols[pos].column_position;
  for (unsigned int i=0; i<save_f_phi_columns->sigf_cols.size(); i++) {
    sigiobs_pos = save_f_phi_columns->sigf_cols[i].column_position;
    if (sigiobs_pos == (iobs_pos + 1)) {
      save_f_phi_columns->selected_refmac_sigiobs_col = i;
    }
  }
}

// as for F+/F- but with I
void 
refmac_dialog_ipm_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting refmac i+/- obs position %d\n", pos); 
  GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_ip_col = pos;
  save_f_phi_columns->selected_refmac_im_col = pos+1;
  // get the corresponding I- and SigI+/I-
  int sigip_pos;
  int sigim_pos;
  int ip_pos = save_f_phi_columns->ipm_cols[pos].column_position;
  int im_pos = save_f_phi_columns->ipm_cols[pos+1].column_position;
  for (unsigned int i=0; i<save_f_phi_columns->sigipm_cols.size()-1; i++) {
    sigip_pos = save_f_phi_columns->sigf_cols[i].column_position;
    sigim_pos = save_f_phi_columns->sigf_cols[i+1].column_position;
    if ((im_pos == ip_pos+2 && sigip_pos == ip_pos+1 && sigim_pos == im_pos+1) ||
	(im_pos == ip_pos+1 && sigip_pos == ip_pos+2 && sigim_pos == im_pos+2)) {
      save_f_phi_columns->selected_refmac_sigip_col = i;
      save_f_phi_columns->selected_refmac_sigim_col = i+1;
    }
  }
}

void 
refmac_dialog_r_free_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting r free position %d\n", pos); 
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_r_free_col = pos;
}

void 
refmac_dialog_phases_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting phases position %d\n", pos); 
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_phi_col = pos;
}

void 
refmac_dialog_fom_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting fom position %d\n", pos); 
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_fom_col = pos;
}

void 
refmac_dialog_hl_button_select(GtkWidget *item, GtkPositionType pos) { 
   
  printf("setting hl position %d\n", pos); 
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(window));
  save_f_phi_columns->selected_refmac_hla_col = pos;
  // get the 3 following HLs
  save_f_phi_columns->selected_refmac_hlb_col = pos + 1;
  save_f_phi_columns->selected_refmac_hlc_col = pos + 2;
  save_f_phi_columns->selected_refmac_hld_col = pos + 3;
}

// end new 

void
coot::fill_f_optionmenu(GtkWidget *optionmenu_f, short int is_expert_mode_flag) { 

   /* Tinker with optionmenu1, the selection of Fs */

   /* create a menu for the optionmenu button.  The various column
    labels will be added to this menu as menuitems*/
  GtkWidget *optionmenu1_menu;
  GtkWidget *menuitem;
  unsigned int i;
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
   
   unsigned int i;
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
   *f_phi_columns = coot::get_mtz_columns(filename);
   f_phi_columns->mtz_filename = filename;

//    std::cout << "DEBUG:: in column_selector_using_cmtz got read success of "
// 	     << f_phi_columns->read_success << std::endl;
      
   if (f_phi_columns->read_success == 0 ) { /*  not a valid mtz file */
      std::cout << "INFO:: data file " << filename << " is not a valid mtz file\n";
      is_phs = try_read_phs_file(filename.c_str()); /* Try reading the data file 
						   as an XtalView .phs file */ 
      if (is_phs == 0) { 
	 int imol_new = try_read_cif_file(filename.c_str());

	 if (coot::util::file_name_extension(filename) == ".fcf") { 
	    if (is_valid_map_molecule(imol_new)) {
	       graphics_info_t g;
	       g.scroll_wheel_map = imol_new;  // change the current scrollable map.
	       g.activate_scroll_radio_button_in_display_manager(imol_new);
	    }
	 }
	 

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
   set_transient_and_position(COOT_MTZ_COLUMN_SELECTOR_DIALOG, column_label_window);
   gtk_object_set_user_data(GTK_OBJECT(column_label_window), f_phi_columns);


   // ----------------------- comboboxes! ----------------------

   coot::column_selector_using_cmtz_setup_comboboxes(column_label_window,
						     f_phi_columns);


   // ------------------- back to option menus ------------------

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

// where should this go?
namespace coot {
   void on_column_label_combobox_changed(GtkComboBox *combobox, gpointer user_data);
}

// Thsee are for labels that don't have callbacks
// The active values of these combo boxes are read when some other widget
// is activated.
//
// Make a copy of this for std::vector<std::string>.
void
my_combo_box_text_add_labels(GtkComboBox *combobox,
			     const std::vector<coot::mtz_type_label> &labels,
			     int active_label_index) {

   GtkListStore *store = gtk_list_store_new(1, G_TYPE_STRING);
   GtkTreeIter iter;
   for (unsigned int ii=0; ii<labels.size(); ii++) {
      const std::string &col_lab = labels[ii].column_label;
      gtk_list_store_append(store, &iter);
      gtk_list_store_set(store, &iter, 0, col_lab.c_str(), -1);
   }

   // does't do anything
   g_signal_connect(combobox, "changed", G_CALLBACK(coot::on_column_label_combobox_changed), NULL);
   
   GtkTreeModel *model = GTK_TREE_MODEL(store);
   GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
   gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combobox), renderer, TRUE);

   gtk_combo_box_set_model(GTK_COMBO_BOX(combobox), model);

   if (active_label_index >= 0)
      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), active_label_index);

}

void
coot::column_selector_using_cmtz_setup_comboboxes(GtkWidget *column_label_window,
						  coot::mtz_column_types_info_t *f_phi_columns) {

   GtkWidget *amplitudes_combobox_w = lookup_widget(column_label_window,
						  "column_selector_amplitudes_combobox");
   GtkWidget *phases_combobox_w = lookup_widget(column_label_window,
					      "column_selector_phases_combobox");
   GtkWidget *weights_combobox_w = lookup_widget(column_label_window,
						 "column_selector_weights_combobox");
   GtkComboBox *amplitudes_combobox = GTK_COMBO_BOX(amplitudes_combobox_w);
   GtkComboBox *phases_combobox     = GTK_COMBO_BOX(phases_combobox_w);
   GtkComboBox *weights_combobox    = GTK_COMBO_BOX(weights_combobox_w);

   const mtz_column_types_info_t &col_labs = *f_phi_columns;

   std::vector<std::string> labels;
   my_combo_box_text_add_labels(amplitudes_combobox, col_labs.f_cols,      0);
   my_combo_box_text_add_labels(phases_combobox,     col_labs.phi_cols,    0);
   my_combo_box_text_add_labels(weights_combobox,    col_labs.weight_cols, 0);

}

void
coot::on_column_label_combobox_changed(GtkComboBox *combobox,
				       gpointer user_data) {

   // std::cout << "changed" << std::endl; // happens on set-active
}
