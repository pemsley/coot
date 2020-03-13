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
// #include <GL/glut.h> // for glutInit()

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

#if 0
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

#endif // GTK2 option-menus

  // comboboxes
  GtkWidget *fobs_combobox    = lookup_widget(window, "column_label_selector_refmac_fobs_combobox");
  GtkWidget *sigfobs_combobox = lookup_widget(window, "column_label_selector_refmac_sigfobs_combobox");
  GtkWidget *rfree_combobox   = lookup_widget(window, "column_label_selector_refmac_rfree_combobox");

  my_combo_box_text_add_items(GTK_COMBO_BOX(fobs_combobox), col_labs.f_cols, 0);
  my_combo_box_text_add_items(GTK_COMBO_BOX(sigfobs_combobox), col_labs.sigf_cols, 0);
  my_combo_box_text_add_items(GTK_COMBO_BOX(rfree_combobox), col_labs.r_free_cols, 0);


}


/* used when the column label widget is being created   */
void
coot::setup_refmac_parameters_from_file(GtkWidget *window) {

   // this function got trashed in a merge. Was it a gtk2 thing?
   // Replace if needed (obviously).

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
   coot::mtz_column_types_info_t *save_f_phi_columns =
      (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "save_f_phi_columns");

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


GtkWidget *make_menu_item(gchar        *name,
                          GCallback callback,
                          gpointer      data ) {
    GtkWidget *item = gtk_menu_item_new_with_label(name);
    g_signal_connect (G_OBJECT (item), "activate", callback, data);
    gtk_widget_show (item);
    return(item);
}


void
phase_button_select(GtkWidget *item, GtkPositionType pos) {

   /*     printf("setting phase position %d\n", pos);   */
   GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");

  save_f_phi_columns->selected_phi_col = pos;
}

void
weight_button_select(GtkWidget *item, GtkPositionType pos) {

/*    printf("setting weight  position %d\n", pos);  */
   GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
   save_f_phi_columns->selected_weight_col = pos;
}


void
refmac_f_button_select(GtkWidget *item, GtkPositionType pos) {

  printf("setting refmac f obs position %d\n", pos);
   GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
  save_f_phi_columns->selected_refmac_fobs_col = pos;
}

void
refmac_sigf_button_select(GtkWidget *item, GtkPositionType pos) {

   printf("setting refmac sigf position %d\n", pos);
   GtkWidget *window = lookup_widget(item, "column_label_window");
   coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
   save_f_phi_columns->selected_refmac_sigfobs_col = pos;
}

void
refmac_r_free_button_select(GtkWidget *item, GtkPositionType pos) {

  printf("setting r free position %d\n", pos);
  GtkWidget *window = lookup_widget(item, "column_label_window");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
  save_f_phi_columns->selected_refmac_r_free_col = pos;
}

// copy of above for newer interface
void
refmac_dialog_f_button_select(GtkWidget *item, GtkPositionType pos) {

  printf("setting refmac f obs position %d\n", pos);
  GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
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
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
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
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
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
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
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
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
  save_f_phi_columns->selected_refmac_r_free_col = pos;
}

void
refmac_dialog_phases_button_select(GtkWidget *item, GtkPositionType pos) {

  printf("setting phases position %d\n", pos);
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
  save_f_phi_columns->selected_refmac_phi_col = pos;
}

void
refmac_dialog_fom_button_select(GtkWidget *item, GtkPositionType pos) {

  printf("setting fom position %d\n", pos);
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
  save_f_phi_columns->selected_refmac_fom_col = pos;
}

void
refmac_dialog_hl_button_select(GtkWidget *item, GtkPositionType pos) {

  printf("setting hl position %d\n", pos);
   GtkWidget *window = lookup_widget(item, "run_refmac_dialog");
  coot::mtz_column_types_info_t *save_f_phi_columns
     = (coot::mtz_column_types_info_t *) g_object_get_data(G_OBJECT(window), "f_phi_columns");
  save_f_phi_columns->selected_refmac_hla_col = pos;
  // get the 3 following HLs
  save_f_phi_columns->selected_refmac_hlb_col = pos + 1;
  save_f_phi_columns->selected_refmac_hlc_col = pos + 2;
  save_f_phi_columns->selected_refmac_hld_col = pos + 3;
}

// end new


////////////////////////////////////////////////////////////////////////////
//              column_selector_using_cmtz()
////////////////////////////////////////////////////////////////////////////
//
GtkWidget *
coot::column_selector_using_cmtz(const std::string &filename) {

   std::cout << "debug:: in coot::column_selector_using_cmtz " << filename << std::endl;

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

   if (false)
      std::cout << "----------- f_phi_columns->mtz_filename "
		<< f_phi_columns->mtz_filename << " was attached to f_phi_columns "
		<< f_phi_columns << std::endl;

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

/* Recall that save_f_phi_columns is now attached to this widget */

   /* Stuff a pointer to mtz info into the dialog: */
   column_label_window = create_column_label_window();
   set_transient_and_position(COOT_MTZ_COLUMN_SELECTOR_DIALOG, column_label_window);

   // modern: g_object_set_data(G_OBJECT(column_label_window), "f_phi_columns", f_phi_columns);
   // old but consistend with gtk_object_get_user_data
   // gtk_object_set_user_data(GTK_OBJECT(column_label_window), f_phi_columns);
   g_object_set_data(G_OBJECT(column_label_window), "f_phi_columns", f_phi_columns);



   // ----------------------- comboboxes! ----------------------

   coot::column_selector_using_cmtz_setup_comboboxes(column_label_window,
						     f_phi_columns);


   /* By default, we want the use weights checkbutton to be off */

   check_weights = GTK_CHECK_BUTTON(lookup_widget(column_label_window,
						  "use_weights_checkbutton"));
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_weights), FALSE);

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
// Pass active_label_index == -1 if you don't have a prefered index.
//
// Make a copy of this for std::vector<std::string>.
void
my_combo_box_text_add_items(GtkComboBox *combobox,
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

int
my_combobox_get_imol(GtkComboBox *combobox) {

   int imol = -1;

   if (combobox) {
      graphics_info_t g;
      imol = g.combobox_get_imol(combobox);
   }
   return imol;
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

   // If the f_col   is called xxx/FWT  prefer that
   // If the phi_col is called xxx/PHWT prefer that

   int f_prefered_idx = col_labs.get_prefered_f_col_idx();
   int p_prefered_idx = col_labs.get_prefered_phi_col_idx();

   std::vector<std::string> labels;
   my_combo_box_text_add_items(amplitudes_combobox, col_labs.f_cols,      f_prefered_idx);
   my_combo_box_text_add_items(phases_combobox,     col_labs.phi_cols,    p_prefered_idx);
   my_combo_box_text_add_items(weights_combobox,    col_labs.weight_cols, 0);

}

void
coot::on_column_label_combobox_changed(GtkComboBox *combobox, gpointer user_data) {

   // std::cout << "on_column_label_combobox_changed() changed" << std::endl;

   // Change the "is-difference-map" checkbutton when we choose a amplitude label that begins with
   // "DEL"

   GtkCheckButton *is_diff_map_checkbutton = GTK_CHECK_BUTTON(lookup_widget(GTK_WIDGET(combobox), "difference_map_checkbutton"));

   GtkWidget *amplidues_combobox = lookup_widget(GTK_WIDGET(combobox), "column_selector_amplitudes_combobox");
   if (GTK_WIDGET(combobox) != amplidues_combobox) return;

   graphics_info_t g;
   std::string al = g.get_active_label_in_combobox(GTK_COMBO_BOX(combobox));
   std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(al);
   if (p.second.length() > 3) {
      std::string ss = p.second.substr(0,3);
      if (ss == "DEL") {
         if (!gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton))) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton), TRUE);
         }
      } else {
         if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton))) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton), FALSE);
         }
      }
   } else {
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton))) {
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(is_diff_map_checkbutton), FALSE);
      }
   }

}
