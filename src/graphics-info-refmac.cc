/* src/graphics-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
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

// See also the below function, which should be used in future.
//
// c.f. the other function:
// graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu,
// 						   GtkSignalFunc signal_func,
//						   int imol_active_position).
//

// There is very little here worth saving.

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
// #include <GL/glut.h>  // for some reason...  // Eh?

#include <iostream>
#include <dirent.h>   // for refmac dictionary files

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER && !defined WINDOWS_MINGW
#include <unistd.h>
#else
//#include "coot-sysdep.h"
#endif

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "clipper/core/map_utils.h" // Map_stats
#include "skeleton/graphical_skel.h"

#include "interface.h"

#include "molecule-class-info.h"
#include "skeleton/BuildCas.h"

#include "gl-matrix.h" // for baton rotation
#include "trackball.h" // for baton rotation

#include "analysis/bfkurt.hh"

#include "globjects.h"
#include "ligand/ligand.hh"
#include "graphics-info.h"

#include "ligand/dunbrack.hh"

#include "utils/coot-utils.hh"

#include "cmtz-interface.hh"
#include "cmtz-interface-gui.hh"

#include "manipulation-modes.hh"

#include "guile-fixups.h"

#if 0
int
graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu,
						   GtkSignalFunc signal_func) {

   return fill_option_menu_with_map_options_generic(option_menu, signal_func);
}
#endif

#if 0
int
graphics_info_t::fill_option_menu_with_map_mtz_options(GtkWidget *option_menu,
						       GtkSignalFunc signal_func) {

   int imol_active = imol_refinement_map;
   return fill_option_menu_with_map_options_generic(option_menu, signal_func, imol_active);
}
#endif

int
graphics_info_t::fill_combobox_with_map_mtz_options(GtkWidget *combobox, GCallback signal_func,
						    int imol_active) {

   int imol = fill_combobox_with_map_options(combobox, signal_func, imol_active);
   return imol;
}


#if 0

int
graphics_info_t::fill_option_menu_with_map_options_generic(GtkWidget *option_menu,
                                                           GtkSignalFunc signal_func,
                                                           int mtz_only) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   GtkWidget *menuitem;
   int active_map_mol_no = -1;

   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();
   int menu_index = 0;

   for (int i=0; i<n_molecules(); i++) {
      if (is_valid_map_molecule(i)) {
         if (not mtz_only ||
             (molecules[i].Refmac_mtz_filename().size() > 0 && mtz_only)) {
                std::string label = clipper::String(i);
                label += " ";
                label += molecules[i].name_for_display_manager();


                menuitem = gtk_menu_item_new_with_label(label.c_str());
                if (signal_func)
                   gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
                                      GTK_SIGNAL_FUNC(signal_func),
                                      GINT_TO_POINTER(i));
                g_object_set_data(G_OBJECT(menuitem), "map_molecule_number", GINT_TO_POINTER(i));
                gtk_menu_append(GTK_MENU(menu), menuitem);
                gtk_widget_show(menuitem);
                if (active_map_mol_no == -1) {
                   active_map_mol_no = i;
                   gtk_menu_set_active(GTK_MENU(menu), menu_index);
                }
             }
      }
      menu_index++;
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
   return active_map_mol_no;
}
#endif


#if 0
// c.f. the other function:
// graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu,
// 						   GtkSignalFunc signal_func)
//
void
graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu,
						   GtkSignalFunc signal_func,
						   int imol_active_position) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   GtkWidget *menuitem;

   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();
   int menu_index = 0;

   for (int i=0; i<n_molecules(); i++) {
      if (is_valid_map_molecule(i)) {
	 char s[200];
	 snprintf(s,199,"%d", i);
	 std::string ss(s);
	 ss += " ";
	 ss += molecules[i].name_;
	 menuitem = gtk_menu_item_new_with_label(ss.c_str());
	 gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			    GTK_SIGNAL_FUNC(signal_func),
			    GINT_TO_POINTER(i));
	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 gtk_widget_show(menuitem);
	 if (i == imol_active_position)
	    gtk_menu_set_active(GTK_MENU(menu), menu_index);
	 menu_index++; // setup for next round
      }
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
}
#endif


#if 0
void
graphics_info_t::fill_option_menu_with_difference_map_options(GtkWidget *option_menu,
							      GtkSignalFunc signal_func,
							      int imol_active_position) {

   std::vector<int> maps_vec;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].is_difference_map_p())
	 maps_vec.push_back(i);
   }

   fill_option_menu_with_map_options_internal(option_menu, signal_func, maps_vec,
					      imol_active_position);

}
#endif


#if 0
void
graphics_info_t::fill_option_menu_with_map_options_internal(GtkWidget *option_menu,
							    GtkSignalFunc signal_func,
							    std::vector<int> map_molecule_numbers,
							    int imol_active_position) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   GtkWidget *menuitem;

   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();
   int menu_index = 0;

   for (unsigned int imap=0; imap<map_molecule_numbers.size(); imap++) {
      int imap_int = imap;
      int i = map_molecule_numbers[imap];
     if (is_valid_map_molecule(i)) {
	 char s[200];
	 snprintf(s,199,"%d", i);
	 std::string ss(s);
	 ss += " ";
	 ss += molecules[i].name_;
	 menuitem = gtk_menu_item_new_with_label(ss.c_str());
	 gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			    GTK_SIGNAL_FUNC(signal_func),
			    GINT_TO_POINTER(i));
	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 gtk_widget_show(menuitem);
	 if (imap_int == imol_active_position)
	    gtk_menu_set_active(GTK_MENU(menu), menu_index);
	 menu_index++; // setup for next round
      }
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
}
#endif


int
graphics_info_t::fill_combobox_with_map_options(GtkWidget *combobox,
						GCallback signal_func,
						int imol_active_position) {

   // delete this function on merge - hmm what did I mean by that?
   std::vector<int> maps_vec;
   for (int i=0; i<n_molecules(); i++)
      if (is_valid_map_molecule(i))
	 maps_vec.push_back(i);

   fill_combobox_with_molecule_options(combobox, signal_func, imol_active_position,
				       maps_vec);

   return -1; // Hmm. Needs checking
}


void
graphics_info_t::fill_combobox_with_difference_map_options(GtkWidget *combobox,
							   GCallback signal_func,
							   int imol_active_position) {

   std::vector<int> maps_vec;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].is_difference_map_p())
	 maps_vec.push_back(i);
   }

   fill_combobox_with_molecule_options(combobox, signal_func, imol_active_position, maps_vec);

}


#if 0
// These are of course *maps*.
void
graphics_info_t::fill_option_menu_with_refmac_options(GtkWidget *option_menu) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();


   GtkWidget *menuitem;

   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].Have_sensible_refmac_params()) {
	 char s[200];
	 snprintf(s, 199, "%d", i);
	 std::string ss(s);
	 ss += " ";
	 ss += molecules[i].name_;
	 menuitem = gtk_menu_item_new_with_label(ss.c_str());

	 // We do a menu_get_active in
	 // save_go_to_atom_mol_menu_active_position.  Hmmm... Does
	 // that function exist?  I don't see it!
	 //
	 // we set user data on the menu item, so that when this goto
	 // Atom widget is cancelled, we can whatever was the molecule
	 // number corresponding to the active position of the menu
	 //
	 // Should be freed in on_go_to_atom_cancel_button_clicked
	 // (callbacks.c)
	 //
	 gtk_object_set_user_data(GTK_OBJECT(menuitem), GINT_TO_POINTER(i));

	 gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			    GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select),
			    GINT_TO_POINTER(i));
	 gtk_menu_append(GTK_MENU(menu), menuitem);
	 gtk_widget_show(menuitem);
      }
   }
//    gtk_menu_set_active(GTK_MENU(menu), 0);
   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu),
			    menu);

}
#endif

void
graphics_info_t::fill_combobox_with_refmac_methods_options(GtkWidget *combobox) {

   std::cout << "in fill_combobox_with_refmac_methods_options " << combobox << std::endl;

   std::vector<std::string> v;
   v.push_back("restrained refinement ");
   v.push_back("rigid body refinement ");
   v.push_back("TLS & restrained refinement ");
   for (unsigned int i=0; i<v.size(); i++)
      gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), v[i].c_str());

   gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 0);
   g_signal_connect(combobox, "changed",
		    G_CALLBACK(refmac_refinement_method_combobox_changed), NULL);
}

void
graphics_info_t::fill_combobox_with_refmac_phase_input_options(GtkWidget *combobox) {

   std::vector<std::string> v;
   v.push_back("no prior phase information");
   v.push_back("phase and FOM");
   v.push_back("Hendrickson-Lattman coefficients");
   v.push_back("SAD data directly");
   for (unsigned int i=0; i<v.size(); i++)
      gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), v[i].c_str());

   gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 0);
   g_signal_connect(combobox, "changed",
		    G_CALLBACK(refmac_refinement_phase_info_combobox_changed), NULL);
}

// put this in the class
void
fill_combobox_with_refmac_mtz_file_options_changed(GtkWidget *combobox, gpointer data) {

   std::cout << "mtz file option changed " << data << std::endl;

}

// These are mtz files actually.  Change the name of this function
void
graphics_info_t::fill_combobox_with_refmac_mtz_file_options(GtkWidget *combobox) {

   // maybe we need gtk_combo_box_text_remove_all here.

   std::vector<std::pair<int, std::string> > mtz_files;
   for (int i=0; i<n_molecules(); i++) {
      if (is_valid_map_molecule(i)) {
	 std::string mtz_filename = molecules[i].Refmac_mtz_filename();
	 if (! mtz_filename.empty()) {
	    bool already_in_list = false;
	    for (unsigned int k=0; k<mtz_files.size(); k++) {
	       if (mtz_filename == mtz_files[k].second) {
		  already_in_list = true;
		  break;
	       }
	    }
	    if (! already_in_list)
	       mtz_files.push_back(std::pair<int, std::string> (i, mtz_filename));
	 }
      }
   }
   // now, using mtz_files, fill the menu and connect signals
   //
   if (mtz_files.size() > 0) {
      for (unsigned int j=0; j<mtz_files.size(); j++) {
	 int i = mtz_files[j].first;
	 const std::string &mtz_filename = mtz_files[j].second;
	 gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox), mtz_filename.c_str());
      }

      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 0);

      GCallback signal_func = G_CALLBACK(fill_combobox_with_refmac_mtz_file_options_changed);
      g_signal_connect(combobox, "changed", signal_func, NULL);

   }

}

// to fill the labels directly from a from an mtz file (used in TWIN refinement)
void
graphics_info_t::fill_combobox_with_refmac_file_labels_options(GtkWidget *combobox) {

   GtkWidget *mtz_file_label = lookup_widget(combobox, "run_refmac_mtz_file_label");

   for (int i=0; i<n_molecules(); i++) {
      std::string fn = molecules[i].Refmac_file_mtz_filename();
      if (! fn.empty()) {
	 gtk_label_set_text(GTK_LABEL(mtz_file_label), fn.c_str());
      }
   }

}

void
graphics_info_t::fill_combobox_with_refmac_ncycles_options(GtkWidget *combobox) {

   GtkComboBoxText *cb = GTK_COMBO_BOX_TEXT(combobox);
   gtk_combo_box_text_append_text(cb, "1");
   gtk_combo_box_text_append_text(cb, "2");
   gtk_combo_box_text_append_text(cb, "5");
   gtk_combo_box_text_append_text(cb, "8");
   gtk_combo_box_text_append_text(cb, "10");
   gtk_combo_box_text_append_text(cb, "16");
   gtk_combo_box_text_append_text(cb, "24");

   gtk_combo_box_set_active(GTK_COMBO_BOX(combobox), 3);

}

#if 0
void
graphics_info_t::fill_option_menu_with_refmac_methods_options(GtkWidget *option_menu) {

  GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
  GtkWidget *menuitem;

  if (menu)
    gtk_widget_destroy(menu);
  menu = gtk_menu_new();

  std::vector<std::string> v;
  v.push_back("restrained refinement ");
  v.push_back("rigid body refinement ");
  v.push_back("TLS & restrained refinement ");

  for (unsigned int i=0; i<v.size(); i++) {
    menuitem = gtk_menu_item_new_with_label((char *) v[i].c_str());
    gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
    		       GTK_SIGNAL_FUNC(graphics_info_t::refmac_change_refinement_method),
    		       GINT_TO_POINTER(i));
    gtk_menu_append(GTK_MENU(menu), menuitem);
    gtk_widget_show(menuitem);
  }

  // active the right setting
  gtk_menu_set_active(GTK_MENU(menu), refmac_refinement_method);

  gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

}
#endif

#if 0
void
graphics_info_t::fill_option_menu_with_refmac_phase_input_options(GtkWidget *option_menu) {

  GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
  GtkWidget *menuitem;

  if (menu)
    gtk_widget_destroy(menu);
  menu = gtk_menu_new();

  std::vector<std::string> v;
  v.push_back("no prior phase information");
  v.push_back("phase and FOM");
  v.push_back("Hendrickson-Lattman coefficients");
  v.push_back("SAD data directly");

  for (unsigned int i=0; i<v.size(); i++) {
    menuitem = gtk_menu_item_new_with_label((char *) v[i].c_str());
    gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
    		       GTK_SIGNAL_FUNC(graphics_info_t::refmac_change_phase_input),
    		       GINT_TO_POINTER(i));
    gtk_menu_append(GTK_MENU(menu), menuitem);
    gtk_widget_show(menuitem);
  }

  // active the right setting
  gtk_menu_set_active(GTK_MENU(menu), refmac_phase_input);

  gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

}
#endif

#if 0
// These are mtz files actually.  Change the name of this function
void
graphics_info_t::fill_option_menu_with_refmac_labels_options(GtkWidget *option_menu) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   GtkWidget *menuitem;
   std::vector<std::pair<int, std::string> > mtz_files;
   for (int i=0; i<n_molecules(); i++) {
      if (is_valid_map_molecule(i)) {

	 // first make a list with all mtz files and at the same time filter out duplicates
	 //
	 std::string mtz_filename = molecules[i].Refmac_mtz_filename();
	 if (! mtz_filename.empty()) {
	    bool already_in_list = false;
	    for (unsigned int k=0; k<mtz_files.size(); k++) {
	       if (mtz_filename == mtz_files[k].second) {
		  already_in_list = true;
		  break;
	       }
	    }
	    if (! already_in_list) {
	       mtz_files.push_back(std::pair<int, std::string> (i, mtz_filename));
	    }
	 }
      }
   }

   // now, using mtz_files, fill the menu and connect signals
   //
   for (unsigned int j=0; j<mtz_files.size(); j++) {
      int i = mtz_files[j].first;
      std::string mtz_filename = mtz_files[j].second;

      menuitem = gtk_menu_item_new_with_label(mtz_filename.c_str());

      // We do a menu_get_active in
      // save_go_to_atom_mol_menu_active_position.  Hmmm... Does
      // that function exist?  I don't see it!
      //
      // we set user data on the menu item, so that when this goto
      // Atom widget is cancelled, we can whatever was the molecule
      // number corresponding to the active position of the menu
      //
      // Should be freed in on_go_to_atom_cancel_button_clicked
      // (callbacks.c)
      //
      gtk_object_set_user_data(GTK_OBJECT(menuitem), GINT_TO_POINTER(i));

      gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			 GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select_add_columns),
			 GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(menu), menuitem);
      gtk_widget_show(menuitem);
   }

   // set the first one active if there is at least one
   if (mtz_files.size() > 0) {
      gtk_menu_set_active(GTK_MENU(menu), 0);
   }

   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

}
#endif


#if 0
// to fill the labels directly from a from an mtz file (used in TWIN refinement)
void
graphics_info_t::fill_option_menu_with_refmac_file_labels_options(GtkWidget *option_menu) {

  GtkWidget *mtz_file_label = lookup_widget(option_menu, "run_refmac_mtz_file_label");

  std::string file_mtz_filename;
  std::string label_filename = gtk_label_get_text(GTK_LABEL(mtz_file_label));
  if (coot::file_exists(label_filename)) {
    coot::setup_refmac_parameters_from_file(option_menu); // check the widget?!
  } else {
    // check if we have a saved filename
    const gchar *saved_filename = saved_refmac_file_filename;
    if (!saved_filename) {
      // pre-select a filename if we have an old twin_mtz_file
      for (int i=0; i<n_molecules(); i++) {
	// first make a list with all mtz files and at the same time filter out dublicates
	if (molecules[i].Refmac_file_mtz_filename().size() > 0) {
	  file_mtz_filename = molecules[i].Refmac_file_mtz_filename();
	  gtk_label_set_text(GTK_LABEL(mtz_file_label), file_mtz_filename.c_str());
	}
      }
    } else {
      file_mtz_filename = saved_filename;
    }
    if (coot::file_exists(file_mtz_filename)) {
      coot::setup_refmac_parameters_from_file(option_menu);
      gtk_label_set_text(GTK_LABEL(mtz_file_label), file_mtz_filename.c_str());
    } else {
      // we dont have any mtz files given
      // delete the contents of the menu(s)
      GtkWidget *fiobs_optionmenu  = lookup_widget(option_menu, "refmac_dialog_fiobs_optionmenu");
      GtkWidget *fiobs_menu        = gtk_option_menu_get_menu(GTK_OPTION_MENU(fiobs_optionmenu));
      GtkWidget *r_free_optionmenu = lookup_widget(option_menu, "refmac_dialog_rfree_optionmenu");
      GtkWidget *r_free_menu       = gtk_option_menu_get_menu(GTK_OPTION_MENU(r_free_optionmenu));
      GtkWidget *fobs_optionmenu   = lookup_widget(option_menu, "refmac_dialog_fobs_optionmenu");
      GtkWidget *fobs_menu         = gtk_option_menu_get_menu(GTK_OPTION_MENU(fobs_optionmenu));
      GtkWidget *fpm_optionmenu    = lookup_widget(option_menu, "refmac_dialog_fpm_optionmenu");
      GtkWidget *fpm_menu          = gtk_option_menu_get_menu(GTK_OPTION_MENU(fpm_optionmenu));
      GtkWidget *ipm_optionmenu    = lookup_widget(option_menu, "refmac_dialog_ipm_optionmenu");
      GtkWidget *ipm_menu          = gtk_option_menu_get_menu(GTK_OPTION_MENU(ipm_optionmenu));
      GtkWidget *phases_optionmenu = lookup_widget(option_menu, "refmac_dialog_phases_optionmenu");
      GtkWidget *phases_menu       = gtk_option_menu_get_menu(GTK_OPTION_MENU(phases_optionmenu));
      GtkWidget *fom_optionmenu    = lookup_widget(option_menu, "refmac_dialog_fom_optionmenu");
      GtkWidget *fom_menu          = gtk_option_menu_get_menu(GTK_OPTION_MENU(fom_optionmenu));
      GtkWidget *hl_optionmenu     = lookup_widget(option_menu, "refmac_dialog_hl_optionmenu");
      GtkWidget *hl_menu           = gtk_option_menu_get_menu(GTK_OPTION_MENU(hl_optionmenu));

      GtkWidget *menu;
      if (fiobs_menu) {
	gtk_widget_destroy(fiobs_menu);
      }
      fiobs_menu = gtk_menu_new();

      if (r_free_menu) {
	gtk_widget_destroy(r_free_menu);
      }
      r_free_menu = gtk_menu_new();

      if (fobs_menu) {
	gtk_widget_destroy(fobs_menu);
      }
      fobs_menu = gtk_menu_new();

      if (fpm_menu) {
	gtk_widget_destroy(fpm_menu);
      }
      fpm_menu = gtk_menu_new();

      if (ipm_menu) {
	gtk_widget_destroy(ipm_menu);
      }
      ipm_menu = gtk_menu_new();

      if (phases_menu) {
	gtk_widget_destroy(phases_menu);
      }
      phases_menu = gtk_menu_new();

      if (fom_menu) {
	gtk_widget_destroy(fom_menu);
      }
      fom_menu = gtk_menu_new();

      if (hl_menu) {
	gtk_widget_destroy(hl_menu);
      }
      hl_menu = gtk_menu_new();

    }
  }
}
#endif

#if 0
void
graphics_info_t::fill_option_menu_with_refmac_ncycle_options(GtkWidget *option_menu) {

  GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
  GtkWidget *menuitem;

  if (menu)
    gtk_widget_destroy(menu);
  menu = gtk_menu_new();

  std::vector<int> v = *preset_number_refmac_cycles;

  for (unsigned int i=0; i<v.size(); i++) {
    menuitem = gtk_menu_item_new_with_label((char *) int_to_string(v[i]).c_str());
    gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
    		       GTK_SIGNAL_FUNC(graphics_info_t::refmac_change_ncycles),
    		       GINT_TO_POINTER(i));
    gtk_menu_append(GTK_MENU(menu), menuitem);
    gtk_widget_show(menuitem);
  }

  // activate the correct setting:
  int found = 0;
  int target_cycle = refmac_ncycles;
  for (unsigned int i=0; i<v.size(); i++) {
    if (v[i] == target_cycle) {
      gtk_menu_set_active(GTK_MENU(menu), i);
      found = 1;
      break;
    }
  }
  if (! found) {
    std::cout <<"INFO:: could not find given no of cycles " << target_cycle <<
      " in preset list. Set to default 5!" <<std::endl;
    gtk_menu_set_active(GTK_MENU(menu), 4);
  }

  gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
}
#endif

void
graphics_info_t::add_refmac_ncycle_no(int &cycle) {

  preset_number_refmac_cycles->push_back(cycle);
}

// a static function
void
graphics_info_t::refinement_map_select(GtkWidget *item, GtkPositionType pos) {
   graphics_info_t g;
   g.set_refinement_map(pos);
}

void
graphics_info_t::refinement_map_select_add_columns(GtkWidget *item, GtkPositionType pos) {

   coot::setup_refmac_parameters_from_file(item);
   graphics_info_t g;
   g.set_refinement_map(pos);
}

// static
void
graphics_info_t::select_refinement_map_combobox_changed(GtkWidget *combobox, gpointer data) {

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox));
   g.set_refinement_map(imol);

}


void
graphics_info_t::set_refmac_phase_input(int phase_flag) {

  graphics_info_t g;

  switch (phase_flag) {

  case coot::refmac::NO_PHASES:
    g.refmac_phase_input = coot::refmac::NO_PHASES;
    break;

  case coot::refmac::PHASE_FOM:
    g.refmac_phase_input = coot::refmac::PHASE_FOM;
    break;

  case coot::refmac::HL:
    g.refmac_phase_input = coot::refmac::HL;
    break;

  case coot::refmac::SAD:
    g.refmac_phase_input = coot::refmac::SAD;
    break;

  default:
    g.refmac_phase_input = coot::refmac::NO_PHASES;
    break;
  }
}

void
graphics_info_t::refmac_change_phase_input(GtkWidget *item, GtkPositionType pos) {

  set_refmac_phase_input(pos);

}

// static
void
graphics_info_t::refmac_refinement_phase_info_combobox_changed(GtkWidget *combobox, gpointer data) {

   std::vector<std::string> v;
   v.push_back("no prior phase information");
   v.push_back("phase and FOM");
   v.push_back("Hendrickson-Lattman coefficients");
   v.push_back("SAD data directly");

   gchar *at = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   if (at) {
      std::string active_text(at);
      for (std::size_t i=0; i<v.size(); i++) {
	 if (active_text == v[i]) {
	    set_refmac_phase_input(i);
	 }
      }
   }
}


void
graphics_info_t::set_refmac_refinement_method(int method) {

  graphics_info_t g;

  switch (method) {

  case coot::refmac::RESTRAINED:
    g.refmac_refinement_method = coot::refmac::RESTRAINED;
    break;

  case coot::refmac::RIGID_BODY:
    g.refmac_refinement_method = coot::refmac::RIGID_BODY;
    break;

  case coot::refmac::RESTRAINED_TLS:
    g.refmac_refinement_method = coot::refmac::RESTRAINED_TLS;
    break;

  default:
    g.refmac_refinement_method = coot::refmac::RESTRAINED;
    break;
  }
}

void
graphics_info_t::refmac_change_refinement_method(GtkWidget *item, GtkPositionType pos) {

  set_refmac_refinement_method(pos);

}


// static
void
graphics_info_t::refmac_refinement_method_combobox_changed(GtkWidget *combobox, gpointer data) {

   gchar *at = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));
   if (at) {
      std::string active_text(at);
      std::vector<std::string> v;
      v.push_back("restrained refinement ");
      v.push_back("rigid body refinement ");
      v.push_back("TLS & restrained refinement ");
      for (std::size_t i=0; i<v.size(); i++) {
	 if (active_text == v[i]) {
	    set_refmac_refinement_method(i);
	 }
      }
   }
}


void
graphics_info_t::set_refmac_n_cycles(int no_cycles) {

  if (no_cycles < 0 || no_cycles > 100) {
    std::cout<< "INFO:: number of cycles out of 'normal' range (0-100). Reset to 5." << std::endl;
    no_cycles = 5;
  }
  refmac_ncycles = no_cycles;

}

void
graphics_info_t::refmac_change_ncycles(GtkWidget *item, GtkPositionType pos) {

  std::string no_str;
  int ncycles;
  std::vector<int> v;
  v = *preset_number_refmac_cycles;
  ncycles = v[pos];
  if (ncycles < 0 || ncycles > 100) {
    std::cout<< "INFO:: number of cycles out of 'normal' range (0-100). Reset to 5." << std::endl;
    ncycles = 5;
  }
  set_refmac_n_cycles(ncycles);

}


void
graphics_info_t::set_refmac_use_tls(int state) {

  graphics_info_t g;

  switch (state) {

  case coot::refmac::TLS_ON:
    g.refmac_use_tls_flag = coot::refmac::TLS_ON;
    break;

  case coot::refmac::TLS_OFF:
    g.refmac_use_tls_flag = coot::refmac::TLS_OFF;
    break;

  default:
    g.refmac_use_tls_flag = coot::refmac::TLS_ON;
    break;
  }
}


void
graphics_info_t::set_refmac_use_twin(int state) {

  graphics_info_t g;

  switch (state) {

  case coot::refmac::TWIN_ON:
    g.refmac_use_twin_flag = coot::refmac::TWIN_ON;
    // we switch off SAD then (for now)
    g.refmac_use_sad_flag = coot::refmac::SAD_OFF;
    break;

  case coot::refmac::TWIN_OFF:
    g.refmac_use_twin_flag = coot::refmac::TWIN_OFF;
    break;

  default:
    g.refmac_use_twin_flag = coot::refmac::TWIN_OFF;
    break;
  }
}


void
graphics_info_t::set_refmac_use_sad(int state) {

  graphics_info_t g;

  switch (state) {

  case coot::refmac::SAD_ON:
    g.refmac_use_sad_flag = coot::refmac::SAD_ON;
    // we switch off TWIN then (for now)
    g.refmac_use_twin_flag = coot::refmac::TWIN_OFF;
    break;

  case coot::refmac::SAD_OFF:
    g.refmac_use_sad_flag = coot::refmac::SAD_OFF;
    break;

  default:
    g.refmac_use_sad_flag = coot::refmac::SAD_OFF;
    break;
  }
}


void
graphics_info_t::set_refmac_use_ncs(int state) {

  graphics_info_t g;

  switch (state) {

  case coot::refmac::NCS_ON:
    g.refmac_use_ncs_flag = coot::refmac::NCS_ON;
    break;

  case coot::refmac::NCS_OFF:
    g.refmac_use_ncs_flag = coot::refmac::NCS_OFF;
    break;

  default:
    g.refmac_use_ncs_flag = coot::refmac::NCS_ON;
    break;
  }
}

void
graphics_info_t::set_refmac_use_intensities(int state) {

  graphics_info_t g;

  switch (state) {

  case coot::refmac::AMPLITUDES:
    g.refmac_use_intensities_flag = coot::refmac::AMPLITUDES;
    break;

  case coot::refmac::INTENSITIES:
    g.refmac_use_intensities_flag = coot::refmac::INTENSITIES;
    break;

  default:
    g.refmac_use_intensities_flag = coot::refmac::AMPLITUDES;
    break;
  }
}

void
graphics_info_t::set_refmac_used_mtz_file(int state) {

  graphics_info_t g;

  switch (state) {

  case coot::refmac::MTZ:
    g.refmac_used_mtz_file_flag = coot::refmac::MTZ;
    break;

  case coot::refmac::MAP:
    g.refmac_used_mtz_file_flag = coot::refmac::MAP;
    break;

  default:
    g.refmac_used_mtz_file_flag = coot::refmac::MTZ;
    break;
  }
}



void
graphics_info_t::add_refmac_sad_atom(const char *atom_name, float fp, float fpp, float lambda) {

  coot::refmac::sad_atom_info_t refmac_sad_atom_info(atom_name, fp, fpp, lambda);
  int replaced = 0;
  for (unsigned int i=0; i<graphics_info_t::refmac_sad_atoms.size(); i++) {
    // check for existing name
    if (graphics_info_t::refmac_sad_atoms[i].atom_name == atom_name) {
      graphics_info_t::refmac_sad_atoms[i] = refmac_sad_atom_info;
      replaced = 1;
      break;
    }
  }
  if (not replaced) {
    graphics_info_t::refmac_sad_atoms.push_back(refmac_sad_atom_info);
  }

}

void
graphics_info_t::store_refmac_params(const std::string &mtz_filename,
				     const std::string &fobs_col,
				     const std::string &sigfobs_col,
				     const std::string &r_free_col,
				     int r_free_flag) {

  have_sensible_refmac_params = 1; // true
  refmac_mtz_file_filename = mtz_filename;
  refmac_fobs_col = fobs_col;
  refmac_sigfobs_col = sigfobs_col;
  refmac_r_free_col = r_free_col;
  refmac_r_free_flag_sensible = r_free_flag;

  std::cout << "INFO:: Stored refmac parameters (for file): "
	    << refmac_fobs_col << " "
	    << refmac_sigfobs_col;
  if (r_free_flag)
    std::cout << " " << refmac_r_free_col << " is sensible." << std::endl;
  else
    std::cout << " the r-free-flag is not sensible" << std::endl;
}

#if 0
void
graphics_info_t::update_refmac_column_labels_frame(GtkWidget *map_optionmenu,
						   GtkWidget *fobs_menu, GtkWidget *fiobs_menu, GtkWidget *fpm_menu,
						   GtkWidget *r_free_menu,
						   GtkWidget *phases_menu, GtkWidget *fom_menu, GtkWidget *hl_menu) {

  GtkWidget *optionmenu;
  GtkWidget *menu;
  GtkWidget *dialog = lookup_widget(map_optionmenu, "run_refmac_dialog");
  GtkWidget *mtz_file_radiobutton = lookup_widget(map_optionmenu, "run_refmac_mtz_file_radiobutton");
  int imol_map_refmac = -1;

  coot::mtz_column_types_info_t *saved_f_phi_columns
    = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(dialog));

  if (not refmac_use_twin_flag && not GTK_TOGGLE_BUTTON(mtz_file_radiobutton)->active) {
    menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(map_optionmenu));
    GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));
    imol_map_refmac = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(active_item)));
  } else {
    // if given mtz file: get the parameters from coot::refmac::saved_refmac_parameters
    std::string file_mtz_filename;
    GtkWidget *twin_mtz_label = lookup_widget(map_optionmenu, "run_refmac_mtz_file_label");
#if (GTK_MAJOR_VERSION > 1)
    const gchar *mtz_filename = gtk_label_get_text(GTK_LABEL(twin_mtz_label));
    file_mtz_filename = mtz_filename;
#else
    gchar **mtz_filename = 0;
    gtk_label_get(GTK_LABEL(twin_mtz_label), mtz_filename);
    file_mtz_filename = (char *)mtz_filename;
#endif // GTK
    std::string tmp_mtz;
    for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].Refmac_file_mtz_filename().size() > 0) {
	 std::string tmp_mtz_inner = molecules[i].Refmac_file_mtz_filename();
	 if (tmp_mtz_inner == file_mtz_filename) {
	    g_print("INFO:: update the labels based on map %i\n", i);
	    imol_map_refmac = i;
	 }
      }
    }
  }

  // get existing refmac parameters and refmac phase parameters and set active
  // otherwise default to first elements.
  if (imol_map_refmac > -1 && molecules[imol_map_refmac].Have_sensible_refmac_params()) {
    std::string fobs_string    = molecules[imol_map_refmac].Refmac_fobs_col();
    std::string r_free_string  = molecules[imol_map_refmac].Refmac_r_free_col();
    // now find 'em
    int fobs_sigfobs_pair = 0;
    int f_col_pos;
    bool break_flag = false;
    for (unsigned int i=0; i<saved_f_phi_columns->f_cols.size(); i++) {
      f_col_pos = saved_f_phi_columns->f_cols[i].column_position;
      for (unsigned int j=0; j<saved_f_phi_columns->sigf_cols.size(); j++) {
	if (saved_f_phi_columns->sigf_cols[j].column_position == f_col_pos + 1) {
	  if (saved_f_phi_columns->f_cols[i].column_label == fobs_string) {
	    saved_f_phi_columns->selected_refmac_fobs_col = i;
	    saved_f_phi_columns->selected_refmac_sigfobs_col = j;
	    gtk_menu_set_active(GTK_MENU(fobs_menu), fobs_sigfobs_pair);
	    break_flag = true;
	    break;
	  }
	  fobs_sigfobs_pair += 1;
	}
      }
      if (break_flag) {
	break;
      }
    }

    for (unsigned int i=0; i<saved_f_phi_columns->r_free_cols.size(); i++) {
      if (saved_f_phi_columns->r_free_cols[i].column_label == r_free_string) {
	gtk_menu_set_active(GTK_MENU(r_free_menu), i);
	saved_f_phi_columns->selected_refmac_r_free_col = i;
	break;
      }
    }
  } else {
    // set the F/Is for twin?!
    if (saved_f_phi_columns->f_cols.size() > 0 || saved_f_phi_columns->i_cols.size() > 0) {
      gtk_menu_set_active(GTK_MENU(fiobs_menu), 0);
      // default set first I col (not F col) if exists
      if (saved_f_phi_columns->i_cols.size() > 0) {
	saved_f_phi_columns->selected_refmac_iobs_col = 0;
      } else {
	saved_f_phi_columns->selected_refmac_fobs_col = 0;
      }
    }
    // F+/F- ignoring saved position?! FIXME
    if (saved_f_phi_columns->fpm_cols.size() > 0 && saved_f_phi_columns->sigfpm_cols.size() > 0) {
      gtk_menu_set_active(GTK_MENU(fpm_menu), 0);
      saved_f_phi_columns->selected_refmac_fp_col = 0;
      saved_f_phi_columns->selected_refmac_fm_col = 1;
      saved_f_phi_columns->selected_refmac_sigfp_col = 0;
      saved_f_phi_columns->selected_refmac_sigfm_col = 1;
    }
    // default setting F, SigF and freeR (doublication?)
    if (saved_f_phi_columns->f_cols.size() > 0) {
      gtk_menu_set_active(GTK_MENU(fobs_menu), 0);
      saved_f_phi_columns->selected_refmac_fobs_col = 0;
    }
    if (saved_f_phi_columns->r_free_cols.size() > 0) {
      gtk_menu_set_active(GTK_MENU(r_free_menu), 0);
      saved_f_phi_columns->selected_refmac_r_free_col = 0;
    }
  }

  // update the phase info, no matter if used or not (as it may when changing the phase input)
  // phase & fom
  if (saved_f_phi_columns->phi_cols.size() > 0) {
    gtk_menu_set_active(GTK_MENU(phases_menu), 0);
    saved_f_phi_columns->selected_refmac_phi_col = 0;
  }
  if (saved_f_phi_columns->weight_cols.size() > 0) {
    gtk_menu_set_active(GTK_MENU(fom_menu), 0);
    saved_f_phi_columns->selected_refmac_fom_col = 0;
  }
  // update the phase info, no matter if used or not (as it may when changing the phase input)
  // HLs
  if (saved_f_phi_columns->hl_cols.size() > 3) {
    gtk_menu_set_active(GTK_MENU(hl_menu), 0);
    saved_f_phi_columns->selected_refmac_hla_col = 0;
    saved_f_phi_columns->selected_refmac_hlb_col = 1;
    saved_f_phi_columns->selected_refmac_hlc_col = 2;
    saved_f_phi_columns->selected_refmac_hld_col = 3;
  }

  // if we have saved parameters use these
  if (imol_map_refmac > -1) {
    // find phases
    std::string phib_string = molecules[imol_map_refmac].Refmac_phi_col();
    if (molecules[imol_map_refmac].Have_refmac_phase_params() && phib_string != "") {
      std::string fom_string  = molecules[imol_map_refmac].Refmac_fom_col();
      // now find 'em
      // phase and FOM
      for (unsigned int i=0; i<saved_f_phi_columns->phi_cols.size(); i++) {
	if (saved_f_phi_columns->phi_cols[i].column_label == phib_string) {
	  gtk_menu_set_active(GTK_MENU(phases_menu), i);
	  saved_f_phi_columns->selected_refmac_phi_col = i;
	  break;
	}
      }
      optionmenu = lookup_widget(map_optionmenu, "refmac_dialog_fom_optionmenu");
      menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
      for (unsigned int i=0; i<saved_f_phi_columns->weight_cols.size(); i++) {
	if (saved_f_phi_columns->weight_cols[i].column_label == fom_string) {
	  gtk_menu_set_active(GTK_MENU(fom_menu), i);
	  saved_f_phi_columns->selected_refmac_fom_col = i;
	  break;
	}
      }
    }

    // find HLs
    std::string hla_string  = molecules[imol_map_refmac].Refmac_hla_col();
    int hl_set_pos = 0;
    int hla_pos;
    if (molecules[imol_map_refmac].Have_refmac_phase_params() && hla_string != "") {
      for (unsigned int i=0; i<saved_f_phi_columns->hl_cols.size(); i++) {
	hla_pos = saved_f_phi_columns->hl_cols[i].column_position;
	//check if we have a consecutive set of 4 HLs
	if (saved_f_phi_columns->hl_cols[i+1].column_position == hla_pos + 1 &&
	    saved_f_phi_columns->hl_cols[i+2].column_position == hla_pos + 2 &&
	    saved_f_phi_columns->hl_cols[i+3].column_position == hla_pos + 3) {

	  if (saved_f_phi_columns->hl_cols[i].column_label == hla_string) {
	    saved_f_phi_columns->selected_refmac_hla_col = i;
	    saved_f_phi_columns->selected_refmac_hlb_col = i + 1;
	    saved_f_phi_columns->selected_refmac_hlc_col = i + 2;
	    saved_f_phi_columns->selected_refmac_hld_col = i + 3;
	    gtk_menu_set_active(GTK_MENU(hl_menu), hl_set_pos);
	    break;
	  }
	  hl_set_pos += 1;
	  i += 3;
	}
      }
    }
  }
}
#endif
