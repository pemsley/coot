/* src/graphics-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008 by the University of Oxford
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


#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
#include <GL/glut.h>  // for some reason...  // Eh?

#include <iostream>
#include <dirent.h>   // for refmac dictionary files

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER && !defined WINDOWS_MINGW
#include <unistd.h>
#else
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define	snprintf _snprintf
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#define AddAtomA AddAtom
#endif

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"

#include "clipper/core/map_utils.h" // Map_stats
#include "graphical_skel.h"


#include "interface.h"

#include "molecule-class-info.h"
#include "BuildCas.h"

#include "gl-matrix.h" // for baton rotation
#include "trackball.h" // for baton rotation

#include "bfkurt.hh"

#include "globjects.h"
#include "ligand.hh"
#include "graphics-info.h"

#include "dunbrack.hh"

#include "coot-utils.hh"

//temp
#include "cmtz-interface.hh"



// A few non-class members - should be somewhere else, I guess.
// 
void initialize_graphics_molecules() { 
  graphics_info_t g;
  g.initialize_molecules(); 
}

// return the new molecule number
// static
int graphics_info_t::create_molecule() { 
   int imol = molecules.size();
//    std::cout << "========================== creating molecule number " 
// 	     << imol << " ===========" << std::endl;
   molecules.push_back(molecule_class_info_t(imol));
   return imol;
}

void
graphics_info_t::post_recentre_update_and_redraw() { 

   //
   int t0 = glutGet(GLUT_ELAPSED_TIME);
   for (int ii=0; ii<n_molecules(); ii++) {
      molecules[ii].update_clipper_skeleton();
      molecules[ii].update_map();  // uses statics in graphics_info_t
      // and redraw the screen using the new map
   }
	 
   int t1 = glutGet(GLUT_ELAPSED_TIME);
   std::cout << "Elapsed time for map contouring: " << t1-t0 << "ms" << std::endl;

   for (int ii=0; ii<n_molecules(); ii++) {
      molecules[ii].update_symmetry();
   }
   make_pointer_distance_objects();
   graphics_draw();
}






// e.g. fit type is "Rigid Body Fit" or "Regularization" etc.
//
// if fit_type is "Torsion General" show the Reverse button.
// 
void do_accept_reject_dialog(std::string fit_type, const coot::refinement_results_t &rr) {

   GtkWidget *window = wrapped_create_accept_reject_refinement_dialog();
   GtkWindow *main_window = GTK_WINDOW(lookup_widget(graphics_info_t::glarea, 
						     "window1"));
   GtkWidget *label;
   if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG_DOCKED){
     label = lookup_widget(GTK_WIDGET(window),
				      "accept_dialog_accept_docked_label_string");
   } else {
     label = lookup_widget(GTK_WIDGET(window),
				      "accept_dialog_accept_label_string");
     gtk_window_set_transient_for(GTK_WINDOW(window), main_window);

     // now set the position, if it was set:
     if ((graphics_info_t::accept_reject_dialog_x_position > -100) && 
	 (graphics_info_t::accept_reject_dialog_y_position > -100)) {
       gtk_widget_set_uposition(window,
				graphics_info_t::accept_reject_dialog_x_position,
				graphics_info_t::accept_reject_dialog_y_position);
     }
   }

   update_accept_reject_dialog_with_results(window, coot::CHI_SQUAREDS, rr);
   if (rr.lights.size() > 0){
     if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG_DOCKED){
       add_accept_reject_lights(GTK_WIDGET(main_window), rr);
     } else {
       add_accept_reject_lights(window, rr);
     }
   }
   
   std::string txt = "";
   txt += "Accept ";
   txt += fit_type;
   txt += "?";

   gtk_label_set_text(GTK_LABEL(label), txt.c_str());

   // Was this a torsion general, in which we need to active the reverse button?
   if (fit_type == "Torsion General") {
      GtkWidget *reverse_button = lookup_widget(window, "accept_reject_reverse_button");
      gtk_widget_show(reverse_button);	
   }
   
   if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG_DOCKED){
     // we need to show some individual widget to make sure we get the right amount
     // of light boxes
     GtkWidget *button_box = lookup_widget(GTK_WIDGET(main_window), "hbuttonbox1");
     gtk_widget_show_all(button_box);
     gtk_widget_show(label);
     gtk_widget_show(window);
   } else {
     gtk_widget_show(window);
   }
}

void add_accept_reject_lights(GtkWidget *window, const coot::refinement_results_t &ref_results) {

#if (GTK_MAJOR_VERSION > 1) 
   GtkWidget *frame;
   if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG_DOCKED) {
     frame = lookup_widget(window, "accept_reject_lights_frame_docked");
   } else {
     frame = lookup_widget(window, "accept_reject_lights_frame");
   }
   gtk_widget_show(frame);

   std::vector<std::pair<std::string, std::string> > boxes;
   boxes.push_back(std::pair<std::string, std::string>("Bonds",                    "bonds_"));
   boxes.push_back(std::pair<std::string, std::string>("Angles",                  "angles_"));
   boxes.push_back(std::pair<std::string, std::string>("Torsions",              "torsions_"));
   boxes.push_back(std::pair<std::string, std::string>("Planes",                  "planes_"));
   boxes.push_back(std::pair<std::string, std::string>("Chirals",                "chirals_"));
   boxes.push_back(std::pair<std::string, std::string>("Non-bonded", "non_bonded_contacts_"));
   boxes.push_back(std::pair<std::string, std::string>("Rama",                      "rama_"));
 
   for (unsigned int i_rest_type=0; i_rest_type<ref_results.lights.size(); i_rest_type++) {
      // std::cout << "Lights for " << ref_results.lights[i_rest_type].second << std::endl;
      for (unsigned int ibox=0; ibox<boxes.size(); ibox++) {
	 if (ref_results.lights[i_rest_type].name == boxes[ibox].first) {
	    std::string stub = boxes[ibox].second.c_str();
	    std::string event_box_name;
	    if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG_DOCKED) {
	      event_box_name = stub + "eventbox_docked";
		} else {
	      event_box_name = stub + "eventbox";
	    }
	    GtkWidget *w = lookup_widget(frame, event_box_name.c_str());
	    if (w) { 
	       GtkWidget *p = w->parent;
	       if (boxes[ibox].first != "Rama") { 
		  GdkColor color = colour_by_distortion(ref_results.lights[i_rest_type].value);
		  set_colour_accept_reject_event_box(w, &color);
	       } else {
		  GdkColor color = colour_by_rama_plot_distortion(ref_results.lights[i_rest_type].value);
		  set_colour_accept_reject_event_box(w, &color);
	       } 
	       gtk_widget_show(p);
	    } else {
	       std::cout << "ERROR:: lookup of event_box_name: " << event_box_name
			 << " failed" << std::endl;
	    } 

	    // we do not add labels for the docked box
	    if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG){
	    
	      std::string label_name = stub + "label";
	      GtkWidget *label = lookup_widget(frame, label_name.c_str());
	      gtk_label_set_text(GTK_LABEL(label), ref_results.lights[i_rest_type].label.c_str());
	      gtk_widget_show(label);
	    } else {
		GtkTooltips *tooltips;
		tooltips = GTK_TOOLTIPS(lookup_widget(window, "tooltips"));
		std::string tips_info = ref_results.lights[i_rest_type].label;
		gtk_tooltips_set_tip(tooltips, w, tips_info.c_str(), NULL);
	    }
	 }
      }
   }
#endif // GTK_MAJOR_VERSION   
}

// Actually, it seems that this does not do anything for GTK == 1. So
// the function that calls it is not compiled (for Gtk1).
// 
void set_colour_accept_reject_event_box(GtkWidget *eventbox, GdkColor *col) {

#if (GTK_MAJOR_VERSION == 1)    
  GtkRcStyle *rc_style = gtk_rc_style_new ();
  rc_style->fg[GTK_STATE_NORMAL] = *col;
  // rc_style->color_flags[GTK_STATE_NORMAL] |= GTK_RC_FG; // compiler failure, try...
  GtkRcFlags new_colour_flags = GtkRcFlags(GTK_RC_FG | rc_style->color_flags[GTK_STATE_NORMAL]);
  rc_style->color_flags[GTK_STATE_NORMAL] = new_colour_flags;
  gtk_widget_modify_style (eventbox, rc_style);
#else    
   gtk_widget_modify_bg(eventbox, GTK_STATE_NORMAL, col);
#endif
   
}

GdkColor colour_by_distortion(float dist) {

   GdkColor col;

   col.pixel = 1;
   col.blue  = 0;

   if (dist < 0.0) { 
      // black for negative numbers
      col.red   = 0;
      col.green = 0;
   } else {
      if (dist < 2.0) { 
	 col.red   = 0;
	 col.green = 55535;
      } else {
	 if (dist < 5.0) {
	    col.red   = 55000;
	    col.green = 55000;
	    // col.blue  = 22000;
      } else {
	    if (dist < 8.0) {
	       col.red   = 64000;
	       col.green = 32000;
	    } else {
	       col.red   = 65535;
	       col.green = 0;
	    }
	 }
      }
   }
   return col;
}

GdkColor colour_by_rama_plot_distortion(float plot_value) {

   GdkColor col;

   col.pixel = 1;
   col.blue  = 0;

   if (plot_value < -15.0) { 
      col.red   = 0;
      col.green = 55535;
   } else {
      if (plot_value < -13.0) {
	 col.red   = 55000;
	 col.green = 55000;
	 // col.blue  = 22000;
      } else {
	 if (plot_value < -10.0) {
	    col.red   = 64000;
	    col.green = 32000;
	 } else {
	    col.red   = 65535;
	    col.green = 0;
	 }
      }
   }
   return col;
} 




void
update_accept_reject_dialog_with_results(GtkWidget *accept_reject_dialog,
					 coot::accept_reject_text_type text_type,
					 const coot::refinement_results_t &rr) {

    std::string extra_text = rr.info;
    if (extra_text != "") {
      if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG){
      
	// now look up the label in window and change it.
	GtkWidget *extra_label = lookup_widget(GTK_WIDGET(accept_reject_dialog),
					       "extra_text_label");
      
	if (text_type == coot::CHIRAL_CENTRES)
	  extra_label = lookup_widget(GTK_WIDGET(accept_reject_dialog), "chiral_centre_text_label");
      
	gtk_label_set_text(GTK_LABEL(extra_label), extra_text.c_str());
      } else {

	// we have a docked accept/reject dialog
	GtkWidget *window = lookup_widget(accept_reject_dialog, "window1");
	GtkTooltips *tooltips;
	GtkTooltipsData *td;
	tooltips = GTK_TOOLTIPS(lookup_widget(GTK_WIDGET(window), "tooltips"));

	int cis_pep_warn = extra_text.find("CIS");
	int chirals_warn = extra_text.find("chiral");

	GtkWidget *cis_eventbox = lookup_widget(GTK_WIDGET(accept_reject_dialog),
						       "cis_peptides_eventbox_docked");
	GtkWidget *p = cis_eventbox->parent;
	GtkWidget *chirals_eventbox = lookup_widget(GTK_WIDGET(accept_reject_dialog),"chirals_eventbox_docked");

	td = gtk_tooltips_data_get(chirals_eventbox);
	std::string old_tip = td->tip_text;

	std::string tips_info_chirals;
        std::string tips_info_cis;
	
	if (cis_pep_warn > -1) {
	  // we have extra cis peptides

	  if (chirals_warn > -1) {
	    // remove the chirals warn text
	    string::size_type start_warn = extra_text.find("WARN");
	    string::size_type end_warn = extra_text.size();
	    tips_info_cis = extra_text.substr(start_warn, end_warn); // list the extra cis peptides here?
	  } else {
	    tips_info_cis = extra_text; // list the extra cis peptides?
	  }

	  gtk_tooltips_set_tip(tooltips, cis_eventbox, tips_info_cis.c_str(), NULL);
	  GdkColor red = colour_by_distortion(1000.0);	// red
	  set_colour_accept_reject_event_box(cis_eventbox, &red);
	  gtk_widget_show(p);

	}
	if (chirals_warn > -1) {
	  // we may have some extra chiral text

	  if (cis_pep_warn > -1) {
	    // remove the cis warn text
	    string::size_type start_warn = extra_text.find("WARN");
	    tips_info_chirals = extra_text.substr(0,start_warn) + old_tip;
	  } else {
	    tips_info_chirals = extra_text + old_tip;
	  }

	  gtk_tooltips_set_tip(tooltips, chirals_eventbox, tips_info_chirals.c_str(), NULL);

	}
	if (extra_text == " "){
	// reset the tips
	  gtk_widget_hide(p);
	  gtk_tooltips_set_tip(tooltips, chirals_eventbox, old_tip.c_str(), NULL);
	}
      }		   
    }
    if (rr.lights.size() > 0)
      add_accept_reject_lights(accept_reject_dialog, rr);
}


GtkWidget *
wrapped_create_accept_reject_refinement_dialog() {

  GtkWidget *w;
  if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG_DOCKED){
    w = lookup_widget(GTK_WIDGET(graphics_info_t::glarea), "accept_reject_dialog_frame_docked");
  } else {
    w = create_accept_reject_refinement_dialog();
  }
  graphics_info_t::accept_reject_dialog = w;
  return w;
 }

// static
void
graphics_info_t::info_dialog(const std::string &s) {
   
   if (graphics_info_t::use_graphics_interface_flag) { 
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   }
}


// static 
void
graphics_info_t::store_window_position(int window_type, GtkWidget *widget) {

   // I am not sure that this function is worth anything now.  It will
   // only be useful for dialog destroys that happen from within
   // graphics-info class... i.e. not many.  Most destroys happen from
   // callbacks.c and thus will use the c-interface.cc's version of
   // this function.
   //
   // If you find yourself in future wanting to add stuff here,
   // instead just move the body of the c-interface.cc's function here
   // and make the c-interface.cc just a one line which calls this
   // function.
   
   gint upositionx, upositiony;
   gdk_window_get_root_origin (widget->window, &upositionx, &upositiony);

   if (window_type == COOT_EDIT_CHI_DIALOG) {
      graphics_info_t::edit_chi_angles_dialog_x_position = upositionx;
      graphics_info_t::edit_chi_angles_dialog_y_position = upositiony;
   }

   if (window_type == COOT_ROTAMER_SELECTION_DIALOG) {
      graphics_info_t::rotamer_selection_dialog_x_position = upositionx;
      graphics_info_t::rotamer_selection_dialog_y_position = upositiony;
   }
}

// static 
void
graphics_info_t::set_transient_and_position(int widget_type, GtkWidget *window) {

   GtkWindow *main_window =
      GTK_WINDOW(lookup_widget(graphics_info_t::glarea, "window1"));
   gtk_window_set_transient_for(GTK_WINDOW(window), main_window);
   
   if (widget_type == COOT_EDIT_CHI_DIALOG) {
//       std::cout << "set_transient_and_position 2" 
// 		<< graphics_info_t::edit_chi_angles_dialog_x_position << "\n";
	 
      if (graphics_info_t::edit_chi_angles_dialog_x_position > -100) {
	 if (graphics_info_t::edit_chi_angles_dialog_y_position > -100) {
	    gtk_widget_set_uposition(window,
				     graphics_info_t::edit_chi_angles_dialog_x_position,
				     graphics_info_t::edit_chi_angles_dialog_y_position);
	 }
      }
   }

   if (widget_type == COOT_ROTAMER_SELECTION_DIALOG) {
      if (graphics_info_t::rotamer_selection_dialog_x_position > -100) {
	 if (graphics_info_t::rotamer_selection_dialog_y_position > -100) {
	    gtk_widget_set_uposition(window,
				     graphics_info_t::rotamer_selection_dialog_x_position,
				     graphics_info_t::rotamer_selection_dialog_y_position);
	 }
      }
   }


}


void
graphics_info_t::statusbar_text(const std::string &text) const {

   if (use_graphics_interface_flag) { 
      if (statusbar) {
	 std::string sbt = text;
	 // If it is "too long" chop it down.
	 unsigned int max_width = 130;
	 GtkWidget *main_window = lookup_widget(glarea, "window1");
	 // some conversion between the window width and the max text length
	 max_width = main_window->allocation.width/4 -38;
	 if (sbt.length() > max_width) { // some number
	    // -------------------------
	    //        |                |  
	    //     200-130            200
	    int l = sbt.length();
	    std::string short_text = text.substr(l-max_width, max_width);
	    // std::cout << "short_text length: " << short_text.length() << std::endl;
	    sbt = "..." + short_text;
	 }
	 gtk_statusbar_push(GTK_STATUSBAR(statusbar),
			    statusbar_context_id,
			    sbt.c_str());
      }
   }
}



double graphics_info_t::GetMouseBeginX() const { return mouse_begin_x; };

double graphics_info_t::GetMouseBeginY() const { return mouse_begin_y; };

void graphics_info_t::SetMouseBegin(double x, double y) {

   mouse_begin_x = x;
   mouse_begin_y = y;
}

// static 
GtkWidget *graphics_info_t::wrapped_nothing_bad_dialog(const std::string &label) { 

   GtkWidget *w = NULL;
   if (use_graphics_interface_flag) { 
      w = create_nothing_bad_dialog();
      GtkWidget *label_widget = lookup_widget(w, "nothing_bad_label");
      gtk_label_set_text(GTK_LABEL(label_widget), label.c_str());
   }
      return w;
}

void
graphics_info_t::set_do_anti_aliasing(int state) { 

   if (do_anti_aliasing_flag != state) {
      if (glarea) { 
	 if (make_current_gl_context(glarea)) { 
	       if (state) { 
	       // should we also add a (quality) hint here?
	       glEnable(GL_LINE_SMOOTH);
	       glEnable(GL_BLEND);
	       glBlendFunc(GL_SRC_ALPHA,GL_ZERO);
	       // glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // Thanks Stuart McN.
	    } else {
	       glDisable(GL_LINE_SMOOTH);
	       glDisable(GL_BLEND);
	       glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // Thanks Stuart McN.
	    }
	    graphics_draw();
	 }
      }
   }
   graphics_info_t::do_anti_aliasing_flag = state;
}

void
graphics_info_t::add_cif_dictionary(std::string cif_dictionary_filename,
				    short int show_no_bonds_dialog_maybe_flag) {

   int nbonds = geom_p->init_refmac_mon_lib(cif_dictionary_filename,
					    cif_dictionary_read_number);
   cif_dictionary_read_number++; 
   if (nbonds > 0) { 
      cif_dictionary_filename_vec->push_back(cif_dictionary_filename);
      if (show_no_bonds_dialog_maybe_flag) {
	 display_density_level_this_image = 1;
	 std::string s;
	 s = "Read ";
	 s += int_to_string(nbonds);
	 s += " bond restraints from ";
	 s += cif_dictionary_filename;
	 display_density_level_screen_string = s;
	 statusbar_text(s);
	 graphics_draw();
      }
      std::cout << display_density_level_screen_string << std::endl;
   } else {
      std::cout << "init_refmac_mon_lib "  << cif_dictionary_filename
		<< " had no bond restraints\n";
      if (use_graphics_interface_flag) { 
	 if (show_no_bonds_dialog_maybe_flag) {
	    GtkWidget *widget = create_no_cif_dictionary_bonds_dialog();
	    gtk_widget_show(widget);
	 }
      }
   }
}

void
graphics_info_t::import_all_refmac_cifs() {

   char *env = getenv("COOT_REFMAC_LIB_DIR");
   if (! env) {
      std::cout << "Can't import dictionary because COOT_REFMAC_LIB_DIR is not defined\n";
   } else {

      std::string coot_refmac_lib_dir(env);

      // is coot_refmac_lib_dir a directory?

      struct stat buf;
      int status = stat(coot_refmac_lib_dir.c_str(), &buf);
      if (status != 0) {
	 std::cout << "Error finding directory " << coot_refmac_lib_dir << std::endl;
      } else {
	 if (S_ISDIR(buf.st_mode)) {
	    std::cout << coot_refmac_lib_dir << " is a directory (good). " << std::endl;

	    std::string data_dir = add_dir_file(coot_refmac_lib_dir, "data");
	    std::string monomer_dir = add_dir_file(data_dir, "monomers");

	    // good

	    DIR *lib_dir = opendir(monomer_dir.c_str());
	    if (lib_dir == NULL) {
	       std::cout << "An ERROR occured on opening the directory "
			 << monomer_dir << std::endl;
	    } else {

	       struct dirent *dir_ent;

	       // loop until the end of the filelist (readdir returns NULL)
	       // 
	       while (1) {
		  dir_ent = readdir(lib_dir);
		  if (dir_ent == NULL) {
		     break;
		  } else {
		     
		     std::string sub_dir_part(std::string(dir_ent->d_name));

		     if ( ! (sub_dir_part == ".") ) { 
			std::string subdirname = add_dir_file(monomer_dir, sub_dir_part);
			
			// we need to test that sub_dir_part is a directory:
			// (if not, silently skip over it)
			//
			status = stat(subdirname.c_str() , &buf);
			if (S_ISDIR(buf.st_mode)) { 

			   DIR *sub_dir = opendir(subdirname.c_str());
			
			   if (sub_dir == NULL) {
			      std::cout << "An ERROR occured on opening the subdirectory "
					<< subdirname << std::endl;
			   } else {
			   
			      struct dirent *sub_dir_ent;
			   
			      while (1) {
				 sub_dir_ent = readdir(sub_dir);
				 if (sub_dir_ent == NULL) {
				    break;
				 } else {
				    std::string cif_filename =
				       add_dir_file(subdirname, std::string(sub_dir_ent->d_name));
				    status = stat(cif_filename.c_str(), &buf);
				    if (status == 0) {
				       if (S_ISREG(buf.st_mode)) { 
					  add_cif_dictionary(cif_filename, 0);
				       }
				    }
				 }
			      }
			   }
			   closedir(sub_dir);
			}
		     } // not "."
		  }
	       }
	       closedir(lib_dir);
	    }
	 } else {
	    std::cout << "Failure to import - " << coot_refmac_lib_dir
		      << " is not a directory\n";
	 } 
      }
   } 
}

// a static
std::string
graphics_info_t::add_dir_file(const std::string &dirname, const std::string &filename) {

   std::string r = dirname;
   r += "/";
   r += filename;
   return r;
}


void
graphics_info_t::save_directory_from_fileselection(const GtkWidget *fileselection) {

   const gchar *filename = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fileselection));
   directory_for_fileselection = coot::util::file_name_directory(filename);
   // std::cout << "saved directory name: " << directory_for_fileselection << std::endl;
}

void
graphics_info_t::save_directory_for_saving_from_fileselection(const GtkWidget *fileselection) {

   const gchar *filename = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fileselection));
   directory_for_saving_for_fileselection = coot::util::file_name_directory(filename);
   // std::cout << "saved directory name: " << directory_for_fileselection << std::endl;
}


void
graphics_info_t::set_directory_for_fileselection_string(std::string filename) {

   directory_for_fileselection = filename;
}

void
graphics_info_t::set_directory_for_fileselection(GtkWidget *fileselection) const {

#if (GTK_MAJOR_VERSION > 1)
  if (graphics_info_t::gtk2_file_chooser_selector_flag == coot::CHOOSER_STYLE) {
    set_directory_for_filechooser(fileselection);
  } else {
    if (directory_for_fileselection != "") {
//       std::cout << "set directory_for_fileselection "
// 		<< directory_for_fileselection << std::endl;
      gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
				      directory_for_fileselection.c_str());
    } else {
      // std::cout << "not setting directory_for_fileselection" << std::endl;
    }
  }
#else
   if (directory_for_fileselection != "") {
//       std::cout << "set directory_for_fileselection "
// 		<< directory_for_fileselection << std::endl;
      gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
				      directory_for_fileselection.c_str());
   } else {
      // std::cout << "not setting directory_for_fileselection" << std::endl;
   } 
#endif // GTK_MAJOR_VERSION
}

void 
graphics_info_t::set_file_for_save_fileselection(GtkWidget *fileselection) const { 

   // just like set_directory_for_fileselection actually, but we give
   // it the full filename, not just the directory.

   int imol = save_imol;
   if (imol >= 0 && imol < graphics_info_t::n_molecules()) { 
      std::string stripped_name = 
	 graphics_info_t::molecules[imol].stripped_save_name_suggestion();
      std::string full_name = stripped_name;
      //       if (graphics_info_t::save_coordinates_in_original_dir_flag != 0) 
      if (graphics_info_t::directory_for_saving_for_fileselection != "")
	 full_name = directory_for_saving_for_fileselection + stripped_name;

      std::cout << "INFO:: Setting fileselection with file: " << full_name
		<< std::endl;
      gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
				      full_name.c_str());
   }
}

#if (GTK_MAJOR_VERSION > 1)
void
graphics_info_t::save_directory_from_filechooser(const GtkWidget *fileselection) {

   gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));
   directory_for_filechooser = coot::util::file_name_directory(filename);
   // std::cout << "saved directory name: " << directory_for_fileselection << std::endl;
   std::cout << "saved directory name: " << directory_for_filechooser << std::endl;
   g_free(filename);
}

void
graphics_info_t::save_directory_for_saving_from_filechooser(const GtkWidget *fileselection) {

   gchar *filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileselection));
   directory_for_saving_for_filechooser = coot::util::file_name_directory(filename);
   // std::cout << "saved directory name: " << directory_for_fileselection << std::endl;
   std::cout << "saved directory name and filename: " << directory_for_filechooser <<" "<<filename << std::endl;
   g_free(filename);
}

void
graphics_info_t::set_directory_for_filechooser_string(std::string filename) {

   directory_for_filechooser = filename;
}

void
graphics_info_t::set_directory_for_filechooser(GtkWidget *fileselection) const {

   if (directory_for_filechooser != "") {
       std::cout << "set directory_for_filechooser "
                 << directory_for_filechooser << std::endl;
      gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(fileselection),
                                      directory_for_filechooser.c_str());
   } else {
      // std::cout << "not setting directory_for_fileselection" << std::endl;
   }

}

void
graphics_info_t::set_file_for_save_filechooser(GtkWidget *fileselection) const {
   // just like set_directory_for_filechooser actually, but we give
   // it the full filename, not just the directory.

   int imol = save_imol;
   if (imol >= 0 && imol < graphics_info_t::n_molecules()) {
      std::string stripped_name =
         graphics_info_t::molecules[imol].stripped_save_name_suggestion();
      std::string full_name = stripped_name;

      if (graphics_info_t::directory_for_saving_for_filechooser != "") {
	full_name = directory_for_saving_for_filechooser + stripped_name;
      } else {
	// if we have a directory in the fileselection path we take this
	if (graphics_info_t::directory_for_saving_for_fileselection != "") {
	  directory_for_saving_for_filechooser = graphics_info_t::directory_for_saving_for_fileselection;
	  
	} else {
	  // otherwise we make one
	  gchar *current_dir = g_get_current_dir();
	  full_name = g_build_filename(current_dir, stripped_name.c_str(), NULL);
	  directory_for_saving_for_filechooser = current_dir;
	  g_free(current_dir);
	}
      }

      std::cout << "INFO:: Setting fileselection with file: " << full_name
                << std::endl;
      if (g_file_test(full_name.c_str(), G_FILE_TEST_EXISTS)) {
         gtk_file_chooser_set_filename(GTK_FILE_CHOOSER(fileselection),
                                      full_name.c_str());
	 // we shouldnt need to call set_current_name and the filename
	 // should be automatically set in the entry field, but this seems
	 // to be buggy at least on gtk+-2.0 v. 1.20.13
	 gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(fileselection),
					   stripped_name.c_str());
      } else {
         gtk_file_chooser_set_current_folder(GTK_FILE_CHOOSER(fileselection),
                                      directory_for_saving_for_filechooser.c_str());
         gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(fileselection),
                                      stripped_name.c_str());
      }
   }
}
#endif // GTK_MAJOR_VERSION


// I find this somewhat asthetically pleasing (maybe because the
// display control widgets are uniquely named [which was a bit of a
// struggle in C, I seem to recall]).
// 
// static
void
graphics_info_t::activate_scroll_radio_button_in_display_manager(int imol) {

   graphics_info_t g;
   if (g.display_control_window()) {
      std::string wname = "map_scroll_button_";
      wname += graphics_info_t::int_to_string(imol);
      GtkWidget *w = lookup_widget(g.display_control_window(), wname.c_str());
      if (w) {
	 if (g.scroll_wheel_map == imol) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	 }
      }
   }
}



void
graphics_info_t::setRotationCentre(int index, int imol) {

   PCAtom atom = molecules[imol].atom_sel.atom_selection[index];
   
   float x = atom->x; 
   float y = atom->y; 
   float z = atom->z;

   old_rotation_centre_x = rotation_centre_x; 
   old_rotation_centre_y = rotation_centre_y; 
   old_rotation_centre_z = rotation_centre_z;
   short int do_zoom_flag = 0;

   if (smooth_scroll == 1)
      smooth_scroll_maybe(x,y,z, do_zoom_flag, 100.0); 

   rotation_centre_x = x; 
   rotation_centre_y = y; 
   rotation_centre_z = z;

//    std::cout << "DEBUG:: dynarama_is_displayed[" << imol << "] is "
// 	     << dynarama_is_displayed[imol] << std::endl;

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
   if (w) {
      coot::rama_plot *plot = (coot::rama_plot *)
	 gtk_object_get_user_data(GTK_OBJECT(w));
      int ires = atom->GetSeqNum();
      std::string chainid(atom->GetChainID());
//       std::cout << "DEBUG:: about to plot big square " << plot
//   		<< " with ires=" << ires << " chain: " << chainid << std::endl;
      std::string ins_code(atom->GetInsCode());
      plot->big_square(ires, ins_code, chainid);
   } 
#endif // HAVE_GTK_CANVAS      

   if (environment_show_distances) { 
      update_environment_graphics_object(index, imol);
      // label new centre
      if (environment_distance_label_atom) {
	 molecules[imol].unlabel_last_atom();
	 molecules[imol].add_to_labelled_atom_list(index);
      }
      if (show_symmetry)
	 update_symmetry_environment_graphics_object(index, imol);
   } else { 

      if (label_atom_on_recentre_flag) { 
	 molecules[imol].unlabel_last_atom();
	 molecules[imol].add_to_labelled_atom_list(index);
      }
   }
}

// We need to know the atom index and imol of the last centred atom
// before this gets activated. The atom index and imol must be right
// because they are used without protection i) to get to
// the right molecule_class_info_t and then the right atom.
// 
void
graphics_info_t::update_environment_distances_maybe(int index, int imol) {

   if (environment_show_distances) { 
      if (go_to_atom_molecule() < n_molecules()) {
 	 if (is_valid_model_molecule(imol)) { 
 	    update_environment_graphics_object(index, imol);
 	    if (show_symmetry)
  	       update_symmetry_environment_graphics_object(index, imol);
 	 }
      }
   }
}

// Return imol = -1 if no (close) atoms found.
// 
// index, imol
std::pair<int, int> 
graphics_info_t::get_closest_atom() const {
   // int index, int imol
   
   std::pair <float, int> dist_info;
   float dist_min = 999999999;
   coot::Cartesian rc = RotationCentre();
   int imol_close = -1;
   int index_close = -1;

   for (int imol=0; imol<n_molecules(); imol++) { 

      if (molecules[imol].has_model()) { 
	 dist_info = molecules[imol].nearest_atom(rc);
	 if (dist_info.first < dist_min) { 
	    imol_close = imol;
	    index_close = dist_info.second;
	    dist_min = dist_info.first;
	 } 
      }
   }
   return std::pair<int, int>(index_close, imol_close);
} 


void
graphics_info_t::setRotationCentre(const symm_atom_info_t &symm_atom_info) {

   std::cout << "setRotationCentre by symmetry atom" << std::endl;
   
   // Invalid read according to valgrind
   PCAtom atom = symm_atom_info.trans_sel[symm_atom_info.atom_index];

   if (atom) { 
      float x = atom->x; // invalid read according to valgrind.
      float y = atom->y; // ditto
      float z = atom->z; // ditto
      
      rotation_centre_x = x; 
      rotation_centre_y = y; 
      rotation_centre_z = z;
   } else { 
      std::cout << "ERROR:: NULL atom in setRotationCentre(symm_atom_info_t)\n";
   } 

}

void
graphics_info_t::setRotationCentre(const coot::clip_hybrid_atom &hybrid_atom) { 

   std::cout << "setRotationCentre by symmetry hybrid atom " 
	     << hybrid_atom.atom << " at " 
	     << hybrid_atom.pos << std::endl;
   
   rotation_centre_x = hybrid_atom.pos.x();
   rotation_centre_y = hybrid_atom.pos.y();
   rotation_centre_z = hybrid_atom.pos.z();

} 


void
graphics_info_t::smooth_scroll_maybe(float x, float y, float z,
				     short int do_zoom_and_move_flag,
				     float target_zoom) {

   float xd = x - rotation_centre_x;
   float yd = y - rotation_centre_y;
   float zd = z - rotation_centre_z;

   float pre_zoom = zoom;
   
   float frac = 1/float (smooth_scroll_steps);
   float stepping_x = frac*xd;
   float stepping_y = frac*yd;
   float stepping_z = frac*zd;

   int n_extra_steps = 10;
   float zoom_in = 1.0 + 1/float(2.0*n_extra_steps + smooth_scroll_steps);
   float zoom_out;
   if (zoom_in != 0.0)      // silly user/divide by zero protection
      zoom_out = 1/zoom_in;
   else 
      zoom_out = 1;

   // This bit of code doesn't get executed (practically ever).
   if (smooth_scroll_do_zoom && (pre_zoom < 70.0)) { 
      if ( (xd*xd + yd*yd + zd*zd) > smooth_scroll_limit*smooth_scroll_limit ) {
	 for (int ii=0; ii<n_extra_steps+smooth_scroll_steps; ii++) {
	    graphics_info_t::zoom *= zoom_in; // typically 1.1
	    graphics_draw();
	 }
      }
   }

   if ( (xd*xd + yd*yd + zd*zd) < smooth_scroll_limit*smooth_scroll_limit ) {
      smooth_scroll_on = 1; // flag to stop wirecube being drawn.

      if (0) { 

	 for (int ii=0; ii<smooth_scroll_steps; ii++) {
	    rotation_centre_x += stepping_x;
	    rotation_centre_y += stepping_y;
	    rotation_centre_z += stepping_z;
	    if (do_zoom_and_move_flag)
	       graphics_info_t::zoom = pre_zoom +
		  float(ii+1)*frac*(target_zoom - pre_zoom);
	    graphics_draw();
	 }

      } else {

	 // -6x^2 +6x parametric function
	 if (smooth_scroll_steps > 0) { 
	    float rotation_centre_x_start = rotation_centre_x;
	    float rotation_centre_y_start = rotation_centre_y;
	    float rotation_centre_z_start = rotation_centre_z;
	    for (int ii=0; ii<smooth_scroll_steps; ii++) {
	       float range_frac = float(ii)/float(smooth_scroll_steps);
	       float f_x = (-2*range_frac*range_frac*range_frac + 3*range_frac*range_frac);
	       
	       rotation_centre_x = rotation_centre_x_start + f_x*xd;
	       rotation_centre_y = rotation_centre_y_start + f_x*yd;
	       rotation_centre_z = rotation_centre_z_start + f_x*zd;

	       if (do_zoom_and_move_flag)
		  graphics_info_t::zoom = pre_zoom +
		     float(ii+1)*frac*(target_zoom - pre_zoom);
	       graphics_draw();

	    }
	 }
      }
   }


   // Also not executed generally.
   if (smooth_scroll_do_zoom && pre_zoom < 70.0 ) { 
      if ( (xd*xd + yd*yd + zd*zd) > smooth_scroll_limit*smooth_scroll_limit ) {
	 for (int ii=0; ii<(n_extra_steps + smooth_scroll_steps); ii++) {
	    graphics_info_t::zoom *= zoom_out; // typically 0.9
	    graphics_draw();
	 }
      }
   }

   smooth_scroll_on = 0;
   zoom = pre_zoom;
}

std::vector<int>
graphics_info_t::displayed_map_imols() const {

   std::vector<int> is;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].drawit_for_map == 1) {
	 is.push_back(i);
      }
   }
   return is;
} 

   
void
graphics_info_t::setRotationCentre(coot::Cartesian centre) {

   old_rotation_centre_x = rotation_centre_x; 
   old_rotation_centre_y = rotation_centre_y; 
   old_rotation_centre_z = rotation_centre_z; 

   if (graphics_info_t::smooth_scroll == 1)
      smooth_scroll_maybe(centre.x(), centre.y(), centre.z(),
			  0, 100.0); // don't zoom and dummy value

   rotation_centre_x = centre.get_x();
   rotation_centre_y = centre.get_y();
   rotation_centre_z = centre.get_z();
}

void
graphics_info_t::setRotationCentreAndZoom(coot::Cartesian centre,
					  float target_zoom) {

   old_rotation_centre_x = rotation_centre_x; 
   old_rotation_centre_y = rotation_centre_y; 
   old_rotation_centre_z = rotation_centre_z; 

   if (graphics_info_t::smooth_scroll == 1)
      smooth_scroll_maybe(centre.x(), centre.y(), centre.z(),
			  1, 100.0);

   rotation_centre_x = centre.get_x();
   rotation_centre_y = centre.get_y();
   rotation_centre_z = centre.get_z();
}




void
graphics_info_t::ShowFPS(){

   long t = 0;

   t = glutGet(GLUT_ELAPSED_TIME);
   if (t - graphics_info_t::T0 >= 5000) {
      GLfloat seconds = (t-T0)/1000.0;
      GLfloat fps = GLfloat (Frames)/seconds;

      std::string s = "INFO:: ";
      s += int_to_string(Frames);
      s += " frames in ";
      s += float_to_string(seconds);
      s += " seconds = ";
      s += float_to_string(fps);
      s += " frames/sec";

      graphics_info_t g;
      g.statusbar_text(s);
      std::cout << s << std::endl;
      graphics_info_t::T0=t;
      graphics_info_t::Frames=0;
   }
}

// We need to reset the Frames so that the first time we get a FPS
// response we are not including all those frames that were made
// without the timer being on.
// 
void
graphics_info_t::SetShowFPS(int t) {

   show_fps_flag = t;
   Frames = 0;
} 
   
//
void
graphics_info_t::SetActiveMapDrag(int t) {

   active_map_drag_flag = t;

} 


void
graphics_info_t::set_font_size(int size) {

//    cout << "graphics_info_t::setting atom_label_font_size to "
// 	<< size << endl;
   
   atom_label_font_size = size;
   
   if (size == 2) {
      atom_label_font = GLUT_BITMAP_HELVETICA_12;
   } else {
      if (size < 2) { 
      atom_label_font = GLUT_BITMAP_HELVETICA_10;
      } else {
	if (size == 3) {
	  atom_label_font = GLUT_BITMAP_HELVETICA_18;
	} else {
	  // no all other fonts
	  if (size == 4) {
	    atom_label_font = GLUT_BITMAP_TIMES_ROMAN_10;
	  } else if (size == 5) {
	    atom_label_font = GLUT_BITMAP_TIMES_ROMAN_24;
	  } else if (size == 6) {
	    atom_label_font = GLUT_BITMAP_8_BY_13;
	  } else if (size == 7) {
	    atom_label_font = GLUT_BITMAP_9_BY_15;
	  } else {
	    // somethign above 7 -> reset default
	    atom_label_font = GLUT_BITMAP_HELVETICA_12;
	  }
	}
      }
   }

   // make the labels (if there are any) change now

   graphics_draw();
}


// For the bin?
void
graphics_info_t::update_map_colour_menu()
{
   for (int ii=0; ii<n_molecules(); ii++)
      molecules[ii].update_map_colour_menu_maybe(ii);
}


// virtual trackball
void
graphics_info_t::set_vt_surface(int v){ 

   if (v == 1) { // VT_FLAT
      //
      trackball_size = 8.8;

   } else {  // VT_SPHERICAL
      trackball_size = 0.8; // as it was in the original code
   }
}

int
graphics_info_t::vt_surface_status() const { 
   int status = 2;
   if (trackball_size > 5)
      status = 1;
   return status;
} 

// phs reading
//
std::string
graphics_info_t::get_phs_filename() const {

   return phs_filename; 

}

void
graphics_info_t::set_phs_filename(std::string filename) {

   phs_filename = filename; 

} 



void
graphics_info_t::show_select_map_dialog() {

   if (use_graphics_interface_flag) {

      GtkWidget *widget = create_select_fitting_map_dialog();
      GtkWidget *optionmenu = lookup_widget(GTK_WIDGET(widget),
					    "select_map_for_fitting_optionmenu");
      
      int imol_map = Imol_Refinement_Map();
      
      // If no map has been set before, set the map to the top of the
      // list (if there are maps in the list)
      // 
      if (imol_map == -1) { 
	 for (int imol=0; imol<n_molecules(); imol++) { 
	    if (molecules[imol].has_map()) {
	       imol_refinement_map = imol;
	       break;
	    }
	 }
      }
      // note that this uses one of 2 similarly named function:
      fill_option_menu_with_map_options(optionmenu,
					GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select),
					imol_refinement_map);
      
      
      // Old notes:
      // now activate the first menu item, i.e. creating this menu is as
      // if the first item in the menu was activated:
      //    for (int i=0; i<n_molecules; i++) {
      //       if (molecules[i].xmap_is_filled[0] == 1) {
      // 	 set_refinement_map(i);
      // 	 break;
      //       }
      //    }
      //
      // 10/6/2004: Well, that is in fact totally crap.  What the *user*
      // expects is that when they open this dialog that the default map
      // (the one they have previously chosen) will be in the active
      // position. i.e. the mechanism described above was totaly wrong.
      //
      // What we really want to do is call
      // fill_option_menu_with_map_options and pass it the active item
      // as an argument.
      //
      gtk_widget_show(widget);
   } else {
      std::cout << "No graphics!  Can't make Map Selection dialog.\n";
   } 
   
} 

// See also the below function, which should be used in future.
// 
// c.f. the other function:
// graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu, 
// 						   GtkSignalFunc signal_func,
//						   int imol_active_position).
// 
void
graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu, 
						   GtkSignalFunc signal_func) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   GtkWidget *menuitem;

   //   std::cout << "DEBUG:: menu: " << menu << std::endl;
   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();
   
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].xmap_is_filled[0] == 1) {
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
      }
   }
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);
}

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
      if (molecules[i].xmap_is_filled[0] == 1) {
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

  for (int i=0; i<v.size(); i++) {
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
  v.push_back("Hendrickson-Lattman coefficients ");
  
  for (int i=0; i<v.size(); i++) {
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

// These are mtz files actually.
void
graphics_info_t::fill_option_menu_with_refmac_labels_options(GtkWidget *option_menu) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   if (menu)
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   
   GtkWidget *menuitem;
   
   std::vector<std::pair<int, std::string> > mtz_files;
   for (int i=0; i<n_molecules(); i++) {
      // first make a list with all mtz files and at the same time filter out dublicates
      if (molecules[i].Refmac_mtz_filename().size() > 0) {
	 std::string mtz_filename = molecules[i].Refmac_mtz_filename();
	 char s[200];
	 snprintf(s, 199, "%d", i);
	 std::string ss(s);
	 ss += " ";
	 ss += molecules[i].name_;
	 ss += "";
	 ss += mtz_filename;
	 int exists = 0;
	 for (unsigned int k=0; k<mtz_files.size(); k++) {
	   if (mtz_filename == mtz_files[k].second) {
	     exists = 1;
	     break;
	   }
	 }
	 if (!exists) {
	   mtz_files.push_back(std::pair<int, std::string> (i, mtz_filename));
	 }
      }
   }

   // now fill the menu and connect signals
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

#if (GTK_MAJOR_VERSION > 1)     
     gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select_add_columns),
			GINT_TO_POINTER(i));
#else
     gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select),
			GINT_TO_POINTER(i));
#endif
     gtk_menu_append(GTK_MENU(menu), menuitem);
     gtk_widget_show(menuitem);
   }

   // set the first one active if there is at least one
   if (mtz_files.size() > 0) {
     gtk_menu_set_active(GTK_MENU(menu), 0);
   }
   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu),
			    menu);

}

// to fill the labels directly from a from an mtz file (used in TWIN refinement)
void
graphics_info_t::fill_option_menu_with_refmac_twin_labels_options(GtkWidget *option_menu) {

#if (GTK_MAJOR_VERSION > 1)
  GtkWidget *mtz_file_label = lookup_widget(option_menu, "run_refmac_mtz_file_label");

  std::string twin_mtz_filename;
  std::string label_filename = gtk_label_get_text(GTK_LABEL(mtz_file_label));
  if (coot::file_exists(label_filename)) {
    coot::setup_refmac_parameters_from_file(option_menu); // check the widget?!
  } else {
    // check if we have a saved filename
    const gchar *saved_filename = saved_refmac_twin_filename;
    if (!saved_filename) {
      // pre-select a filename if we have an old twin_mtz_file
      for (int i=0; i<n_molecules(); i++) {
	// first make a list with all mtz files and at the same time filter out dublicates
	if (molecules[i].Refmac_twin_mtz_filename().size() > 0) {
	  twin_mtz_filename = molecules[i].Refmac_twin_mtz_filename();
	  gtk_label_set_text(GTK_LABEL(mtz_file_label), twin_mtz_filename.c_str());
	}
      }
    } else {
      twin_mtz_filename = saved_filename;
    }
    if (coot::file_exists(twin_mtz_filename)) {
      coot::setup_refmac_parameters_from_file(option_menu);
      gtk_label_set_text(GTK_LABEL(mtz_file_label), twin_mtz_filename.c_str());
    } else {
      // we dont have any mtz files given
      // delete the contents of the menu(s), currently only fiobs and r_free
      GtkWidget *fiobs_optionmenu  = lookup_widget(option_menu, "refmac_dialog_fiobs_optionmenu");
      GtkWidget *fiobs_menu  = gtk_option_menu_get_menu(GTK_OPTION_MENU(fiobs_optionmenu));
      GtkWidget *r_free_optionmenu = lookup_widget(option_menu, "refmac_dialog_rfree_optionmenu");
      GtkWidget *r_free_menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(r_free_optionmenu));
      GtkWidget *menu;
      if (fiobs_menu) {
	gtk_widget_destroy(fiobs_menu);
      }
      fiobs_menu = gtk_menu_new();
      if (r_free_menu) {
	gtk_widget_destroy(r_free_menu);
      }
      r_free_menu = gtk_menu_new();
    }
  }

#endif // GTK2
}

void
graphics_info_t::fill_option_menu_with_refmac_ncycle_options(GtkWidget *option_menu) {

  GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
  GtkWidget *menuitem;

  if (menu)
    gtk_widget_destroy(menu);
  menu = gtk_menu_new();
   
  std::vector<int> v = *preset_number_refmac_cycles;

  for (int i=0; i<v.size(); i++) {
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
  for (int i=0; i<v.size(); i++) {
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

  default:
    g.refmac_phase_input = coot::refmac::NO_PHASES;
    break;
  }
}

void
graphics_info_t::refmac_change_phase_input(GtkWidget *item, GtkPositionType pos) {

  set_refmac_phase_input(pos);

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

  default:
    g.refmac_refinement_method = coot::refmac::RESTRAINED;
    break;
  }
}

void
graphics_info_t::refmac_change_refinement_method(GtkWidget *item, GtkPositionType pos) {

  set_refmac_refinement_method(pos);

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
graphics_info_t::add_refmac_sad_atom(const char *atom_name, float fp, float fpp, float lambda) {

  coot::refmac::sad_atom_info_t refmac_sad_atom_info(atom_name, fp, fpp, lambda);
  int replaced = 0;
  for (int i=0; i<graphics_info_t::refmac_sad_atoms.size(); i++) {
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
graphics_info_t::update_refmac_column_labels_frame(GtkWidget *map_optionmenu,
						   GtkWidget *fobs_menu, GtkWidget *fiobs_menu, GtkWidget *fpm_menu,
						   GtkWidget *r_free_menu,
						   GtkWidget *phases_menu, GtkWidget *fom_menu, GtkWidget *hl_menu) {

  std::cout <<"BL DEBUG:: update refmac columns"<< std::endl;
  GtkWidget *optionmenu;
  GtkWidget *menu;
  GtkWidget *dialog = lookup_widget(map_optionmenu, "run_refmac_dialog");
  int imol_map_refmac = -1;

  coot::mtz_column_types_info_t *saved_f_phi_columns
    = (coot::mtz_column_types_info_t *) gtk_object_get_user_data(GTK_OBJECT(dialog));
  
  if (not refmac_use_twin_flag) {
    menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(map_optionmenu));
    GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));
    imol_map_refmac = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(active_item)));
  } else {
    // in TWIN: find (last) map which may have refmac twin parameters corresponding 
    // to given filename.
    std::string twin_mtz_filename;
    GtkWidget *twin_mtz_label = lookup_widget(map_optionmenu, "run_refmac_mtz_file_label");
#if (GTK_MAJOR_VERSION > 1)
    const gchar *mtz_filename = gtk_label_get_text(GTK_LABEL(twin_mtz_label));
    twin_mtz_filename = mtz_filename;
#else
    gchar **mtz_filename = 0;
    gtk_label_get(GTK_LABEL(twin_mtz_label), mtz_filename);
    twin_mtz_filename = (char *)mtz_filename;
#endif // GTK
    std::cout <<"BL DEBUG:: have filename from label in update "<< twin_mtz_filename<<std::endl;
    std::string tmp_mtz;
    for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].Refmac_twin_mtz_filename().size() > 0) {
	std::string tmp_mtz = molecules[i].Refmac_twin_mtz_filename();
	if (tmp_mtz == twin_mtz_filename) {
	  std::cout <<"BL DEBUG:: found existing Refmac ywin file with imol "<< i <<std::endl;
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
    for (int i=0; i<saved_f_phi_columns->f_cols.size(); i++) {
      f_col_pos = saved_f_phi_columns->f_cols[i].column_position;
      for (int j=0; j<saved_f_phi_columns->sigf_cols.size(); j++) {
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
      
    for (int i=0; i<saved_f_phi_columns->r_free_cols.size(); i++) {
      if (saved_f_phi_columns->r_free_cols[i].column_label == fobs_string) {
	gtk_menu_set_active(GTK_MENU(r_free_menu), i);
	saved_f_phi_columns->selected_refmac_r_free_col = i;
	break;
      }
    }
  } else {
    std::cout <<"BL DEBUG:: set default (first) fs" <<std::endl;
    std::cout <<"BL DEBUG:: size of f_cols" <<saved_f_phi_columns->f_cols.size() <<std::endl;
    std::cout <<"BL DEBUG:: size of i_cols" <<saved_f_phi_columns->i_cols.size() <<std::endl;
    if (refmac_use_twin_flag) {
      if (saved_f_phi_columns->f_cols.size() > 0 || saved_f_phi_columns->i_cols.size() > 0) {
	gtk_menu_set_active(GTK_MENU(fiobs_menu), 0);
      }
    } else {
      // default setting F, SigF and freeR (doublication?)
      if (saved_f_phi_columns->f_cols.size() > 0) {
	gtk_menu_set_active(GTK_MENU(fobs_menu), 0);
	// no need to set?!
	//saved_f_phi_columns->selected_refmac_fobs_col = 0;
      }
    }
    if (saved_f_phi_columns->r_free_cols.size() > 0) {
      gtk_menu_set_active(GTK_MENU(r_free_menu), 0);
      // no need to set?!
      //saved_f_phi_columns->selected_refmac_r_free_col = 0;
    } 
  }

  if (imol_map_refmac > -1) {
    // find phases
    std::string phib_string = molecules[imol_map_refmac].Refmac_phi_col();
    if (molecules[imol_map_refmac].Have_refmac_phase_params() && phib_string != "") {
      std::string fom_string  = molecules[imol_map_refmac].Refmac_fom_col();
      // now find 'em
      // phase and FOM
      for (int i=0; i<saved_f_phi_columns->phi_cols.size(); i++) {
	if (saved_f_phi_columns->phi_cols[i].column_label == phib_string) {
	  gtk_menu_set_active(GTK_MENU(phases_menu), i);
	  saved_f_phi_columns->selected_refmac_phi_col = i;
	  break;
	}
      }
      optionmenu = lookup_widget(map_optionmenu, "refmac_dialog_fom_optionmenu");
      menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
      for (int i=0; i<saved_f_phi_columns->weight_cols.size(); i++) {
	if (saved_f_phi_columns->weight_cols[i].column_label == fom_string) {
	  gtk_menu_set_active(GTK_MENU(fom_menu), i);
	  saved_f_phi_columns->selected_refmac_fom_col = i;
	  break;
	}
      }
    } else {
      // update the phase info, no matter if used or not (as it may when changing the phase input)
      // phase & fom
      if (saved_f_phi_columns->phi_cols.size() > 0) {
	gtk_menu_set_active(GTK_MENU(phases_menu), 0);
	saved_f_phi_columns->selected_refmac_phi_col = 0;
      }
      if (saved_f_phi_columns->weight_cols.size() > 0) {
	gtk_menu_set_active(GTK_MENU(fom_menu), 0);
      }
    }

    // find HLs
    std::string hla_string  = molecules[imol_map_refmac].Refmac_hla_col();
    int hl_set_pos = 0;
    int hla_pos;
    if (molecules[imol_map_refmac].Have_refmac_phase_params() && hla_string != "") {
      for (int i=0; i<saved_f_phi_columns->hl_cols.size(); i++) {
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
    } else {
      // update the phase info, no matter if used or not (as it may when changing the phase input)
      // HLs
      if (saved_f_phi_columns->hl_cols.size() > 3) {
	gtk_menu_set_active(GTK_MENU(hl_menu), 0);
	saved_f_phi_columns->selected_refmac_hla_col = 0;
	saved_f_phi_columns->selected_refmac_hlb_col = 1;
	saved_f_phi_columns->selected_refmac_hlc_col = 2;
	saved_f_phi_columns->selected_refmac_hld_col = 3;
      }
    }
  }

}


// static
void
graphics_info_t::skeleton_map_select(GtkWidget *item, GtkPositionType pos) { 

   graphics_info_t g;
   g.map_for_skeletonize = pos;

   // So now we have changed the map to skeletonize (or that has been
   // skeletonized).  We need to look up the radio buttons and set
   // them according to whether this map has a skeleton or not.
   //
   // Recall that for radio buttons, you have to set them on, setting
   // them off doesn't work.
   //
   GtkWidget *on_button  = lookup_widget(item, "skeleton_on_radiobutton");
   GtkWidget *off_button = lookup_widget(item, "skeleton_off_radiobutton");
   if (graphics_info_t::molecules[g.map_for_skeletonize].has_map()) {
      if (graphics_info_t::molecules[g.map_for_skeletonize].fc_skeleton_draw_on) {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button), TRUE);
      } else {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), TRUE);
      }
   }

}


void
graphics_info_t::fill_option_menu_with_coordinates_options(GtkWidget *option_menu, 
							   GtkSignalFunc signal_func,
							   int imol_active_position) {

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
//    std::cout << "option_menu: " << option_menu << std::endl;
//    std::cout << "menu: " << menu << std::endl;
//    std::cout << "menu is widget: " << GTK_IS_WIDGET(menu) << std::endl;
//    std::cout << "menu is menu: " << GTK_IS_MENU(menu) << std::endl;

   // menu is not GTK_MENU on Gtk2 Ubuntu kalypso 64 bit
   if (GTK_IS_MENU(menu))
      gtk_widget_destroy(menu);
   
   menu = gtk_menu_new();

   GtkWidget *menuitem;
   int item_count = 0; 

   for (int imol=0; imol<n_molecules(); imol++) {

      if (molecules[imol].has_model() > 0) { 

	 std::string ss = int_to_string(imol);
	 ss += " " ;
	 int ilen = molecules[imol].name_.length();
	 int left_size = ilen-go_to_atom_menu_label_n_chars_max;
	 if (left_size <= 0) {
	    // no chop
	    left_size = 0;
	 } else {
	    // chop
	    ss += "...";
	 } 

	 ss += molecules[imol].name_.substr(left_size, ilen);

	 menuitem = gtk_menu_item_new_with_label (ss.c_str());

	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     signal_func,
			     GINT_TO_POINTER(imol)); 

	 gtk_menu_append(GTK_MENU(menu), menuitem); 
	 gtk_widget_show(menuitem);
 	 if (imol == imol_active_position) {
 	    gtk_menu_set_active(GTK_MENU(menu), item_count);
	 }
	 item_count++;
      }
   }
   
//    gtk_menu_set_active(GTK_MENU(menu), imol_active_position); 

   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu),
			    menu);

}


void
graphics_info_t::set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame) { 
   GtkWidget *on_button = lookup_widget(skeleton_frame, 
					"skeleton_on_radiobutton");
   GtkWidget *off_button = lookup_widget(skeleton_frame, 
					 "skeleton_off_radiobutton");

   int imol = map_for_skeletonize;
   if (imol >= 0) { 
      if (molecules[imol].xskel_is_filled) { 
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button),  TRUE);
      } else {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), TRUE);
      }
   } 
} 

void 
graphics_info_t::set_on_off_single_map_skeleton_radio_buttons(GtkWidget *skeleton_frame, 
							      int imol) { 
   GtkWidget *on_button = lookup_widget(skeleton_frame, 
					"single_map_skeleton_on_radiobutton");
   GtkWidget *off_button = lookup_widget(skeleton_frame, 
					 "single_map_skeleton_off_radiobutton");

   if (imol >= 0) { 
      if (molecules[imol].xskel_is_filled) { 
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button),  TRUE);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), FALSE);
      } else {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(on_button),  FALSE);
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(off_button), TRUE);
      }
   } 
}

void
graphics_info_t::set_contour_sigma_button_and_entry(GtkWidget *window, int imol) { 

   GtkWidget *entry = lookup_widget(window, "single_map_sigma_step_entry");
   GtkWidget *checkbutton = lookup_widget(window, "single_map_sigma_checkbutton");

   if (imol < n_molecules()) { 
      if (molecules[imol].has_map()) { 
	 float v = molecules[imol].contour_sigma_step;
	 gtk_entry_set_text(GTK_ENTRY(entry), float_to_string(v).c_str());
	 if (molecules[imol].contour_by_sigma_flag) { 
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
	 } else { 
	    gtk_widget_set_sensitive(entry, FALSE);
	 }

	 
	 GtkWidget *level_entry =
	    lookup_widget(window, "single_map_properties_contour_level_entry");
	 float lev = molecules[imol].contour_level[0];
	 gtk_entry_set_text(GTK_ENTRY(level_entry), float_to_string(lev).c_str());
      }
   }

}


// static
void
graphics_info_t::skeletonize_map(short int prune_it, int imol) { 

   graphics_info_t g;

	
   if (imol < 0) { // i.e. is -1

      std::cout << "INFO:: map for skeletonization was not selected" << std::endl;
      
   } else {

      // so that we don't do this when the skeleton is on already:
      //
      if (g.molecules[imol].fc_skeleton_draw_on == 0) {
	 g.molecules[imol].fc_skeleton_draw_on = 1;

	 //       mean_and_variance<float> mv = 
	 // 	 map_density_distribution(g.molecules[imol].xmap_list[0],0); 

	 clipper::Map_stats stats(g.molecules[imol].xmap_list[0]);

	 std::cout << "Mean and sigma of map: " << stats.mean() << " and " 
		   << stats.std_dev() << std::endl; 
      
	 float map_cutoff = stats.mean() + 1.5*stats.std_dev(); 
	 g.skeleton_level = map_cutoff; 
	    
	 // derived from sktest:
	 // 
	 g.molecules[imol].xskel_cowtan.init(g.molecules[imol].xmap_list[0].spacegroup(), 
					     g.molecules[imol].xmap_list[0].cell(),
					     g.molecules[imol].xmap_list[0].grid_sampling());
      
	 std::cout << "INFO:: making skeleton cowtan..." << std::endl; 
	 GraphicalSkel cowtan(g.molecules[imol].xmap_list[0],
			      g.molecules[imol].xskel_cowtan); //fill xskel_cowtan

	 g.molecules[imol].xskel_is_filled = 1; // TRUE
      
	 // various experiments....
      
	 // cowtan.tip_filter(xmap_list[0], &xskl); // tinker with xskel_cowtan
      
	 //cowtan.prune(g.molecules[imol].xmap_list[imap],
	 //	 &g.molecules[imol].xskel_cowtan);
      
	 //
	 cowtan.Pprune(g.molecules[imol].xmap_list[0],
		       &g.molecules[imol].xskel_cowtan,
		       map_cutoff);

	 if (prune_it) { 
	    BuildCas bc(g.molecules[imol].xmap_list[0], map_cutoff); 

	    // mark segments by connectivity
	    // 
	    int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
						       g.molecules[imol].xmap_list[0],
						       map_cutoff); 
	 
	    cout << "INFO:: There were " << nsegments << " different segments" << endl; 
	 
	    bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	    g.molecules[imol].set_colour_skeleton_by_segment(); // use random colours

	 } else { 
	    g.molecules[imol].set_colour_skeleton_by_level(); // use conventional
							      // colouring, (just
							      // sets a flag)
	 } 
      
      
	 // now display the skeleton
      
	 g.molecules[imol].update_clipper_skeleton();
	 graphics_draw();

      } else { 
	 std::cout << "This map has a skeleton already" << std::endl;
      }
   }
}

// static
void
graphics_info_t::set_initial_map_for_skeletonize() { 

   // Initially map_for_skeletonize is -1;
   
   if (graphics_info_t::map_for_skeletonize == -1) { 
      for (int imol=0; imol<n_molecules();imol++) { 
	 if (graphics_info_t::molecules[imol].has_map()) { 
	    graphics_info_t::map_for_skeletonize = imol;
	    break;
	 } 
      }
   }
}

// static
void
graphics_info_t::unskeletonize_map(int imol) {

   graphics_info_t g;

   if (imol >= 0) { 
      g.molecules[imol].unskeletonize_map();
      graphics_draw();
   } else { 
      std::cout << "Map skeleton not selected from optionmenu." << std::endl;
      std::cout << "Please try again and this time, select "
		<< "a map from the optionmenu" << std::endl;
   } 
} 




// Do we need to delete the old regularize_object_bonds_box?
// Yes, well, clear_up() it.
// 
void
graphics_info_t::clear_moving_atoms_object() {

   in_edit_chi_mode_flag = 0;
   in_edit_torsion_general_flag = 0;
   in_moving_atoms_drag_atom_mode_flag = 0; // no more dragging atoms
   have_fixed_points_sheared_drag_flag = 0;
   // and set the rotation translation atom index to unknown again:
   // rot_trans_atom_index_rotation_origin_atom = -1;
   rot_trans_rotation_origin_atom = NULL;

   std::cout << "clearing intermediate object..." << std::endl;
   graphical_bonds_container empty_box;
   regularize_object_bonds_box.clear_up();
   regularize_object_bonds_box = empty_box;

   dynamic_distances.clear();
   
   graphics_draw();
}

void
graphics_info_t::set_refinement_map(int i) {
   imol_refinement_map = i;
}



// Recall that moving_atoms_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
// We want to put the moving parts back into the object that we were
// regularizing (imol_moving_atoms).
// 
void
graphics_info_t::accept_moving_atoms() {

//    std::cout << ":::: INFO:: imol moving atoms is "
// 	     << imol_moving_atoms << std::endl;
   
   if (moving_atoms_asc_type == coot::NEW_COORDS_ADD) {
      molecules[imol_moving_atoms].add_coords(*moving_atoms_asc);
   } else { 
      if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF) { 
	 molecules[imol_moving_atoms].replace_coords(*moving_atoms_asc, 1);
	 update_geometry_graphs(*moving_atoms_asc, imol_moving_atoms);
      } else {
	 if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE) {
	    molecules[imol_moving_atoms].replace_coords(*moving_atoms_asc, 0);
	    update_geometry_graphs(*moving_atoms_asc, imol_moving_atoms);
	 } else { 
	    if (moving_atoms_asc_type == coot::NEW_COORDS_INSERT) {
	       molecules[imol_moving_atoms].insert_coords(*moving_atoms_asc);
	    } else { 
	       if  (moving_atoms_asc_type == coot::NEW_COORDS_INSERT_CHANGE_ALTCONF) {
		  molecules[imol_moving_atoms].insert_coords_change_altconf(*moving_atoms_asc);
	       } else {
		  std::cout << "------------ ERROR! -------------------" << std::endl;
		  std::cout << "       moving_atoms_asc_type not known: ";
		  std::cout << moving_atoms_asc_type << std::endl;
		  std::cout << "------------ ERROR! -------------------" << std::endl;
	       }
	    }
	 }
      }
   }
   
   if (do_probe_dots_post_refine_flag) {
      setup_for_probe_dots_on_chis_molprobity(imol_moving_atoms);
   }
   clear_up_moving_atoms();

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *w = coot::get_validation_graph(imol_moving_atoms, coot::RAMACHANDRAN_PLOT);
   if (w) {
      coot::rama_plot *plot = (coot::rama_plot *)
	 gtk_object_get_user_data(GTK_OBJECT(w));
      // std::cout << "updating rama plot for " << imol_moving_atoms << std::endl;
      handle_rama_plot_update(plot);
   }
#endif // HAVE_GTK_CANVAS || HAVE_GNOME_CANVAS

   update_environment_distances_by_rotation_centre_maybe(imol_moving_atoms);

   normal_cursor(); // we may have had fleur cursor.
   // and set the rotation translation atom index to unknown again:
   // rot_trans_atom_index_rotation_origin_atom = -1;
   rot_trans_rotation_origin_atom = NULL;
   in_moving_atoms_drag_atom_mode_flag = 0; // no more dragging atoms
   have_fixed_points_sheared_drag_flag = 0;
   in_edit_chi_mode_view_rotate_mode = 0;

   if (do_probe_dots_post_refine_flag) {
      do_interactive_probe();
   }
}


void
graphics_info_t::update_environment_distances_by_rotation_centre_maybe(int imol_moving_atoms) {
   
   // Oh this is grimly "long hand".
   graphics_info_t g;
   if (g.environment_show_distances) {
      coot::at_dist_info_t at_d_i = g.molecules[imol_moving_atoms].closest_atom(RotationCentre());
      if (at_d_i.atom) {
	 int atom_index;
	 if (at_d_i.atom->GetUDData(g.molecules[imol_moving_atoms].atom_sel.UDDAtomIndexHandle,
				    atom_index) == UDDATA_Ok) {
	    g.update_environment_distances_maybe(atom_index, imol_moving_atoms);
	 }
      }
   }
}


#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
// coot::rama_plot is an unknown type if we don't have canvas
void
graphics_info_t::handle_rama_plot_update(coot::rama_plot *plot) {

   if (plot) {
      // if it's a normal plot: update it
      if (plot->is_kleywegt_plot()) {
	 // are the molecule numbers from which the kleywegt plot
	 // was generated still valid?
	 std::pair<int, int> p = plot->molecule_numbers();
	 if (graphics_info_t::molecules[p.first].has_model() && 
	     graphics_info_t::molecules[p.second].has_model()) { 
	    std::pair<std::string, std::string> chain_ids = plot->chain_ids();
	    std::cout << "updating kleywegt plot with chain ids :" << chain_ids.first
		      << ": :" << chain_ids.second << ":" << std::endl;
	    if (plot->kleywegt_plot_uses_chain_ids_p())
	       plot->draw_it(p.first, p.second,
			     graphics_info_t::molecules[p.first].atom_sel.mol,
			     graphics_info_t::molecules[p.second].atom_sel.mol,
			     chain_ids.first, chain_ids.second);
	    else 
	       plot->draw_it(p.first, p.second,
			     graphics_info_t::molecules[p.first].atom_sel.mol,
			     graphics_info_t::molecules[p.second].atom_sel.mol);
	 } else {
	    // close down the plot
	    plot->destroy_yourself();
	 }
      } else {
	 plot->draw_it(molecules[imol_moving_atoms].atom_sel.mol);
      } 
   } else {
      std::cout << "ERROR:: (trapped) in accept_moving_atoms attempt to draw to null plot\n";
   } 
}
#endif // HAVE_GNOME_CANVAS or HAVE_GTK_CANVAS


void 
graphics_info_t::clear_up_moving_atoms() { 

   std::cout << "INFO:: graphics_info_t::clear_up_moving_atoms..." << std::endl;
   moving_atoms_asc_type = 0; // unset
   in_moving_atoms_drag_atom_mode_flag = 0; // no more dragging atoms
   have_fixed_points_sheared_drag_flag = 0;
   // and take out any drag refine idle function:
//    std::cout << "DEBUG:: removing token " << drag_refine_idle_function_token
// 	     << " (clear_up_moving_atoms) " << std::endl;
   gtk_idle_remove(drag_refine_idle_function_token); 
   drag_refine_idle_function_token = -1; // magic "not in use" value
   
   if (moving_atoms_asc->atom_selection != NULL) {
      if (moving_atoms_asc->n_selected_atoms > 0) { 
	 moving_atoms_asc->mol->DeleteSelection(moving_atoms_asc->SelectionHandle);
	 moving_atoms_asc->atom_selection = NULL;
      } else {
	 std::cout << "WARNING:: attempting to delete non-NULL ";
	 std::cout << "moving_atoms_asc.atom_selection" << std::endl;
	 std::cout << "but moving_atoms_asc.n_selected_atoms == 0" << std::endl;
	 std::cout << "ignoring " << std::endl;
      }
   } else {
      std::cout << "WARNING:: attempting to delete NULL moving_atoms_asc.atom_selection"
		<< std::endl;
      std::cout << "Ignoring. " << std::endl;
   }
   if (moving_atoms_asc->mol != NULL) {
      if (moving_atoms_asc->n_selected_atoms > 0) { 
	 moving_atoms_asc->mol = NULL; 
      } else {
	 std::cout << "attempting to delete non-NULL moving_atoms_asc.mol" << std::endl;
	 std::cout << "but moving_atoms_asc.n_selected_atoms == 0" << std::endl;
	 std::cout << "ignoring " << std::endl;
      } 
   } else {
      std::cout << "attempting to delete NULL moving_atoms_asc.mol" << std::endl;
      std::cout << "ignoring " << std::endl;
   }

   dynamic_distances.clear();

   // and now the signal that moving_atoms_asc has been cleared:
   //
   moving_atoms_asc->n_selected_atoms = 0;

#ifdef HAVE_GSL
   last_restraints = coot::restraints_container_t(); // last_restraints.size() = 0;
#endif // HAVE_GSL   
}

// static
gint
graphics_info_t::drag_refine_idle_function(GtkWidget *widget) {

#ifdef HAVE_GSL

   int retprog = graphics_info_t::drag_refine_refine_intermediate_atoms();

   if (retprog != GSL_CONTINUE) {
      if (retprog == GSL_ENOPROG)
	 std::cout << " NOPROGRESS" << std::endl;
      if (retprog == GSL_SUCCESS)
	 std::cout << " SUCCESS" << std::endl;
      long t1 = glutGet(GLUT_ELAPSED_TIME);
      std::cout << " TIME:: (dragged refinement): " << float(t1-T0)/1000.0 << std::endl;

      graphics_info_t g;
      g.check_and_warn_bad_chirals_and_cis_peptides();

      gtk_idle_remove(graphics_info_t::drag_refine_idle_function_token);
      graphics_info_t::drag_refine_idle_function_token = -1; // magic "not in use" value
   } 
   // 
   graphics_draw();

#else
   int retprog = -1;
#endif // HAVE_GSL   
   return retprog;
}

// static 
gint
graphics_info_t::drag_refine_refine_intermediate_atoms() {

   int retprog = -1;
#ifdef HAVE_GSL

   graphics_info_t g;

   // While the results of the refinement are a conventional result
   // (unrefined), let's continue.  However, there are return values
   // that we will stop refining and remove the idle function is on a
   // GSL_ENOPROG(RESS) and GSL_SUCCESS.... actually, we will remove
   // it on anything other than a GSL_CONTINUE
   //

   // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
   coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
   // coot::restraint_usage_Flags flags = coot::BONDS_AND_PLANES;

   if (do_torsion_restraints)
      // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED;
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;

   if (do_rama_restraints)
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
   
   if (do_torsion_restraints && do_rama_restraints)
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	    

   // print_initial_chi_squareds_flag is 1 the first time then we turn it off.
   graphics_info_t::saved_dragged_refinement_results = 
      g.last_restraints.minimize(flags, dragged_refinement_steps_per_frame,
				 print_initial_chi_squareds_flag);
   retprog = graphics_info_t::saved_dragged_refinement_results.progress;
   print_initial_chi_squareds_flag = 0;
   int do_disulphide_flag = 0;
   int draw_hydrogens_flag = 0;
   if (molecules[imol_moving_atoms].draw_hydrogens())
      draw_hydrogens_flag = 1;
   Bond_lines_container bonds(*(g.moving_atoms_asc), do_disulphide_flag, draw_hydrogens_flag);
   g.regularize_object_bonds_box.clear_up();
   g.regularize_object_bonds_box = bonds.make_graphical_bonds();


   // Update the Accept/Reject Dialog if it exists (and it should do,
   // if we are doing dragged refinement).
   if (accept_reject_dialog) {
      if (saved_dragged_refinement_results.lights.size() > 0) { 
	 update_accept_reject_dialog_with_results(accept_reject_dialog,
						  coot::CHI_SQUAREDS,
						  saved_dragged_refinement_results);
      }
   }
   
#endif // HAVE_GSL

   return retprog;
}


// static
void
graphics_info_t::add_drag_refine_idle_function() {

   // add a idle function if there isn't one in operation already.
   graphics_info_t g;
   if (g.drag_refine_idle_function_token == -1) {
      g.drag_refine_idle_function_token =
	 gtk_idle_add((GtkFunction)graphics_info_t::drag_refine_idle_function, g.glarea);
      T0 = glutGet(GLUT_ELAPSED_TIME);
      print_initial_chi_squareds_flag = 1;
   }
} 



void
graphics_info_t::set_dynarama_is_displayed(GtkWidget *dyna_toplev, int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   // first delete the old plot for this molecule (if it exists)
   // 
   if (imol < graphics_info_t::n_molecules() && imol >= 0) {

      // Clear out the old one if it was there.
      GtkWidget *w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
      if (w) {
	 coot::rama_plot *plot =
	    (coot::rama_plot *) gtk_object_get_user_data(GTK_OBJECT(w));
	 delete plot;
      }
      coot::set_validation_graph(imol, coot::RAMACHANDRAN_PLOT, dyna_toplev);
   }
#endif // HAVE_GTK_CANVAS   
}



// As these are moving atoms, they cannot point to atoms of a real
// molecule in the molecules array.  [Otherwise as we move moving
// atoms, the static molecule's atoms move (but note that their bonds
// are not updated)].
// 
void
graphics_info_t::make_moving_atoms_graphics_object(const atom_selection_container_t &asc) {


   if (! moving_atoms_asc) {
      moving_atoms_asc = new atom_selection_container_t;
   }
   *moving_atoms_asc = asc;

   // these not needed now, 
//    moving_atoms_asc->mol = asc.mol;
//    moving_atoms_asc->n_selected_atoms = asc.n_selected_atoms;
//    moving_atoms_asc->atom_selection = asc.atom_selection;
//    moving_atoms_asc->read_success = asc.read_success;
//    moving_atoms_asc->SelectionHandle = asc.SelectionHandle;
//    moving_atoms_asc->UDDAtomIndexHandle = asc.UDDAtomIndexHandle;
//    moving_atoms_asc->UDDOldAtomIndexHandle = asc.UDDOldAtomIndexHandle;

   int do_disulphide_flag = 0;

//    std::cout << "DEBUG:: There are " << moving_atoms_asc->n_selected_atoms
// 	     << " atoms in the graphics object" << std::endl;
//    for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
//       std::cout << moving_atoms_asc->atom_selection[i] << std::endl;
//    }
//
//    Needed that to debug doubly clear_up on a bonds box.  Now points
//    reset in clear_up().
   
//    std::cout << "DEBUG:: bonds box type of molecule " << imol_moving_atoms
// 	     << " is " << molecules[imol_moving_atoms].Bonds_box_type()
// 	     << std::endl;
   if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS) {
      
      Bond_lines_container bonds;
      bonds.do_Ca_bonds(*moving_atoms_asc, 1.0, 4.7);
      std::cout << "done CA bonds" << std::endl;
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
      
   } else {
      int draw_hydrogens_flag = 0;
      if (molecules[imol_moving_atoms].draw_hydrogens())
	 draw_hydrogens_flag = 1;
      Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag, draw_hydrogens_flag);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
   }
}


// Display the graphical object of the regularization.
// static
void 
graphics_info_t::moving_atoms_graphics_object() { 
   
   // very much most of the time, this will be zero
   // 
   if (graphics_info_t::regularize_object_bonds_box.num_colours > 0) { 

      // now we want to draw out our bonds in white, 
      glColor3f (0.9, 0.9, 0.9);
      coot::CartesianPair pair;
      Lines_list ll;
      
      glLineWidth(graphics_info_t::bond_thickness_intermediate_atoms);
      for (int i=0; i< graphics_info_t::regularize_object_bonds_box.num_colours; i++) {

	 switch(i) {
	 case blue:
	    glColor3f (0.60, 0.6, 0.99);
	    break;
	 case red:
	    glColor3f (0.95, 0.65, 0.65);
	    break;
	 default:
	    glColor3f (0.8, 0.8, 0.8);
	 }
	 ll = graphics_info_t::regularize_object_bonds_box.bonds_[i];

	 glBegin(GL_LINES); 
	 for (int j=0; j< graphics_info_t::regularize_object_bonds_box.bonds_[i].num_lines; j++) {
	   
	    pair = ll.pair_list[j];
	    
	    glVertex3f(pair.getStart().get_x(),
		       pair.getStart().get_y(),
		       pair.getStart().get_z());
	    glVertex3f(pair.getFinish().get_x(),
		       pair.getFinish().get_y(),
		       pair.getFinish().get_z());
	 }
	 glEnd();
      }
   }
}

// This does (draws) symmetry too.
// 
// static
void
graphics_info_t::environment_graphics_object() {

   graphics_info_t g;
   g.environment_graphics_object_internal(environment_object_bonds_box);
   if (g.show_symmetry)
      g.environment_graphics_object_internal(symmetry_environment_object_bonds_box);
}

// This is the GL rendering of the environment bonds box
// 
void
graphics_info_t::environment_graphics_object_internal(const graphical_bonds_container &env_bonds_box) const {

   if (environment_show_distances == 1) {
      
      if (env_bonds_box.num_colours > 0) {
	 
	 coot::CartesianPair pair;
	 Lines_list ll;
	 coot::Cartesian text_pos;
	 float dist;
	 
	 glEnable(GL_LINE_STIPPLE);
	 glLineStipple (1, 0x00FF);
	 glLineWidth(2.0);
	 for (int i=0; i< env_bonds_box.num_colours; i++) {
	    ll = env_bonds_box.bonds_[i];
	    float it = float(i);
	    if (it > 1.0) 
	       it = 1.0;

	    // now we want to draw out our bonds in various colour,
	    // according to if they have a carbon or not.
	    // 
	    glColor3f (0.8, 0.8-0.4*it, 0.4+0.5*it);
	    
	    for (int j=0; j< env_bonds_box.bonds_[i].num_lines; j++) {
	   
	       pair = ll.pair_list[j];
	    
	       glBegin(GL_LINES); 
	       glVertex3f(pair.getStart().get_x(),
			  pair.getStart().get_y(),
			  pair.getStart().get_z());
	       glVertex3f(pair.getFinish().get_x(),
			  pair.getFinish().get_y(),
			  pair.getFinish().get_z());
	       glEnd();
	       text_pos = pair.getFinish().mid_point(pair.getStart()) +
		  coot::Cartesian(0.0, 0.1, 0.1);
	       glRasterPos3f(text_pos.x(), text_pos.y(), text_pos.z());
	       dist = (pair.getStart() - pair.getFinish()).amplitude();
	       printString(float_to_string(dist));
	    }
	 }
	 glDisable(GL_LINE_STIPPLE);
      }
   }
}


void
graphics_info_t::update_environment_graphics_object(int atom_index, int imol) {

   environment_object_bonds_box =
      molecules[imol].make_environment_bonds_box(atom_index);
}

void
graphics_info_t::update_symmetry_environment_graphics_object(int atom_index, int imol) {

   symmetry_environment_object_bonds_box =
      molecules[imol].make_symmetry_environment_bonds_box(atom_index);
} 

// Just a stub.  Not used currently.
// 
std::string graphics_info_t::make_mmdb_atom_string_from_go_to_atom() { 

   std::string a;
   return a; 

}




// ------------------------------------------------------------------
//                   density level
// ------------------------------------------------------------------
//
void
graphics_info_t::set_density_level_string(int imol, float dlevel) {

   graphics_info_t g;
   float map_sigma = g.molecules[imol].map_sigma();
   
   display_density_level_screen_string = "map " + int_to_string(imol);
   display_density_level_screen_string += " level = ";
   display_density_level_screen_string += float_to_string(dlevel);
   display_density_level_screen_string += "e/A^3 (";
   display_density_level_screen_string += float_to_string(dlevel/map_sigma);
   display_density_level_screen_string += "sigma)";
   
}

// ------------------------------------------------------------------
//                   geometry
// ------------------------------------------------------------------
//
void
graphics_info_t::display_geometry_distance() {

   CAtom *atom1 = molecules[geometry_atom_index_1_mol_no].atom_sel.atom_selection[geometry_atom_index_1];
   CAtom *atom2 = molecules[geometry_atom_index_2_mol_no].atom_sel.atom_selection[geometry_atom_index_2];
      
   clipper::Coord_orth p1(atom1->x, atom1->y, atom1->z);
   clipper::Coord_orth p2(atom2->x, atom2->y, atom2->z);

   double dist = clipper::Coord_orth::length(p1, p2);

   std::cout << "        distance atom 1: "
	     << "(" << geometry_atom_index_1_mol_no << ") " 
	     << atom1->name << "/"
	     << atom1->GetChainID()  << "/"
	     << atom1->GetSeqNum()   << "/"
	     << atom1->GetResName() << std::endl;
   std::cout << "        distance atom 2: "
	     << "(" << geometry_atom_index_2_mol_no << ") " 
	     << atom2->name << "/"
	     << atom2->GetChainID()  << "/"
	     << atom2->GetSeqNum()   << "/"
	     << atom2->GetResName() << std::endl;

   std::pair<clipper::Coord_orth, clipper::Coord_orth> p(p1, p2);
   distance_object_vec->push_back(p);
   graphics_draw();

   std::cout << "        distance: " << dist << " Angstroems" << std::endl;
   std::string s = "Distance: ";
   s += float_to_string(dist);
   s += " A";
   statusbar_text(s);
}

void 
graphics_info_t::display_geometry_distance_symm(const coot::Cartesian &p1, const coot::Cartesian &p2) { 

   std::pair<clipper::Coord_orth, clipper::Coord_orth> p(clipper::Coord_orth(p1.x(), p1.y(), p1.z()), 
							 clipper::Coord_orth(p2.x(), p2.y(), p2.z()));
   distance_object_vec->push_back(p);
   graphics_draw();

} 

void
graphics_info_t::display_geometry_angle() const {

   /// old way, where we dont do symmetry atoms
//    CAtom *atom1 = molecules[geometry_atom_index_1_mol_no].atom_sel.atom_selection[geometry_atom_index_1];
//    CAtom *atom2 = molecules[geometry_atom_index_2_mol_no].atom_sel.atom_selection[geometry_atom_index_2];
//    CAtom *atom3 = molecules[geometry_atom_index_3_mol_no].atom_sel.atom_selection[geometry_atom_index_3];
      
   clipper::Coord_orth p1(angle_tor_pos_1.x(), angle_tor_pos_1.y(), angle_tor_pos_1.z());
   clipper::Coord_orth p2(angle_tor_pos_2.x(), angle_tor_pos_2.y(), angle_tor_pos_2.z());
   clipper::Coord_orth p3(angle_tor_pos_3.x(), angle_tor_pos_3.y(), angle_tor_pos_3.z());


   clipper::Coord_orth v1 = p2 - p1;
   clipper::Coord_orth v2 = p2 - p3;

//    std::cout << "positions for angles"
// 	     << "   " << p1.format() << std::endl
// 	     << "   " << p2.format() << std::endl
// 	     << "   " << p3.format() << std::endl;

//    std::cout << "display_geometry_angle: " << std::endl
//  	     << "      " << v1.format() << std::endl
//  	     << "      " << v2.format() << std::endl;

   double dp = clipper::Coord_orth::dot(v1,v2);
   double len_v1 = sqrt(v1.lengthsq());
   double len_v2 = sqrt(v2.lengthsq());
   len_v1 = len_v1 < 0.0001 ? 0.0001 : len_v1;
   len_v2 = len_v2 < 0.0001 ? 0.0001 : len_v2;
   double cos_theta = dp/(len_v1 * len_v2);
   double theta = acos(cos_theta);

   // no symmetry version, we only have pos now
//    std::cout << "       angle atom 1: "
// 	     << "(" << geometry_atom_index_2_mol_no << ") " 
// 	     << atom1->name << "/"
// 	     << atom1->GetChainID()  << "/"
// 	     << atom1->GetSeqNum()   << "/"
// 	     << atom1->GetResName() << std::endl;
//    std::cout << "       angle atom 2: "
// 	     << "(" << geometry_atom_index_2_mol_no << ") " 
// 	     << atom2->name << "/"
// 	     << atom2->GetChainID()  << "/"
// 	     << atom2->GetSeqNum()   << "/"
// 	     << atom2->GetResName() << std::endl;
//    std::cout << "       angle atom 3: "
// 	     << "(" << geometry_atom_index_2_mol_no << ") " 
// 	     << atom3->name << "/"
// 	     << atom3->GetChainID()  << "/"
// 	     << atom3->GetSeqNum()   << "/"
// 	     << atom3->GetResName() << std::endl;

   std::cout << "       angle: " << theta*57.29578 << " degrees "
	     << std::endl;

   display_density_level_this_image = 1;
   display_density_level_screen_string = "  Angle:  ";
   display_density_level_screen_string += float_to_string(theta*57.29578);
   display_density_level_screen_string += " degrees";
   statusbar_text(display_density_level_screen_string);
   // redraw is in calling function for angles.
}


void
graphics_info_t::display_geometry_torsion() const {

   clipper::Coord_orth p1(angle_tor_pos_1.x(), angle_tor_pos_1.y(), angle_tor_pos_1.z());
   clipper::Coord_orth p2(angle_tor_pos_2.x(), angle_tor_pos_2.y(), angle_tor_pos_2.z());
   clipper::Coord_orth p3(angle_tor_pos_3.x(), angle_tor_pos_3.y(), angle_tor_pos_3.z());
   clipper::Coord_orth p4(angle_tor_pos_4.x(), angle_tor_pos_4.y(), angle_tor_pos_4.z());

   clipper::Coord_orth v1  = p2 - p1;
   clipper::Coord_orth v2  = p2 - p3;
   clipper::Coord_orth v3  = p3 - p4;
   clipper::Coord_orth v3r = p4 - p3;

//    std::cout << ":p1: " << p1.format() << std::endl;
//    std::cout << ":p2: " << p2.format() << std::endl;
//    std::cout << ":p3: " << p3.format() << std::endl;
//    std::cout << ":p4: " << p4.format() << std::endl;

//     std::cout << "display_geometry_angle: " << std::endl
//  	     << "      " << v1.format() << std::endl
//  	     << "      " << v2.format() << std::endl;

   // double dp = clipper::Coord_orth::dot(v1,v2);
   double len_v1 = sqrt(v1.lengthsq());
   double len_v2 = sqrt(v2.lengthsq());
   double len_v3 = sqrt(v3.lengthsq());
   len_v1 = len_v1 < 0.0001 ? 0.0001 : len_v1;
   len_v2 = len_v2 < 0.0001 ? 0.0001 : len_v2;
   len_v3 = len_v3 < 0.0001 ? 0.0001 : len_v3;
   // double cos_theta = dp/(len_v1 * len_v2);
   // double theta = acos(cos_theta);

   // we could do this when we kept the atom indcies for the torsion.
   // Now we just save the positions, so we can't give this nice
   // geometry list, ho hum...
//    std::cout << "       angle atom 1: "
// 	     << "(" << geometry_atom_index_1_mol_no << ") " 
// 	     << atom1->name << "/"
// 	     << atom1->GetChainID()  << "/"
// 	     << atom1->GetSeqNum()   << "/"
// 	     << atom1->GetResName() << std::endl;
//    std::cout << "       angle atom 2: "
// 	     << "(" << geometry_atom_index_2_mol_no << ") " 
// 	     << atom2->name << "/"
// 	     << atom2->GetChainID()  << "/"
// 	     << atom2->GetSeqNum()   << "/"
// 	     << atom2->GetResName() << std::endl;
//    std::cout << "       angle atom 3: "
// 	     << "(" << geometry_atom_index_3_mol_no << ") " 
// 	     << atom3->name << "/"
// 	     << atom3->GetChainID()  << "/"
// 	     << atom3->GetSeqNum()   << "/"
// 	     << atom3->GetResName() << std::endl;
//    std::cout << "       angle atom 4: "
// 	     << "(" << geometry_atom_index_4_mol_no << ") " 
// 	     << atom4->name << "/"
// 	     << atom4->GetChainID()  << "/"
// 	     << atom4->GetSeqNum()   << "/"
// 	     << atom4->GetResName() << std::endl;

   double tors = clipper::Coord_orth::torsion(p1, p2, p3, p4);
   double torsion = clipper::Util::rad2d(tors);
   std::cout << "       torsion: " << torsion << " degrees "
	     << std::endl;
   
   display_density_level_this_image = 1;
   display_density_level_screen_string = "  Torsion:  ";
   display_density_level_screen_string += float_to_string(torsion);
   display_density_level_screen_string += " degrees";
   statusbar_text(display_density_level_screen_string);
   graphics_draw();

}




void
graphics_info_t::pepflip() {

   molecules[imol_pepflip].pepflip(atom_index_pepflip);
   normal_cursor();
   model_fit_refine_unactive_togglebutton("model_refine_dialog_pepflip_togglebutton");
}
 
// ----------------------------------------------------------
//                     some utilities:
// ----------------------------------------------------------

std::string
graphics_info_t::int_to_string(int i) {
   char s[100];
   for (int ii=0; ii<100; ii++) s[ii]=0;
   snprintf(s, 99, "%d", i);
   return std::string(s);
}

std::string
graphics_info_t::float_to_string(float f) {
   char s[100];
   // initial s, stop valgrind complaining
   for (int i=0; i<100; i++) s[i]=0;
   snprintf(s,99,"%5.2f",f);
   return std::string(s);
}

std::string
graphics_info_t::float_to_string_using_dec_pl(float f, unsigned short int n_dec_pl) {
   char s[100];
   for (int i=0; i<100; i++) s[i]=0;
   snprintf(s,99,"%7.4f",f); // haha, FIXME. (use n_dec_pl, not 4)
   return std::string(s);
}



int
graphics_info_t::Imol_Refinement_Map() const { 

   int nmaps = 0;
   int only_map = -1;  // gets set

   if (imol_refinement_map != -1) { // has been set by user

      if (imol_refinement_map < n_molecules())
	 if (imol_refinement_map >= 0)
	    if (molecules[imol_refinement_map].has_map())
	       return imol_refinement_map;
   }
      
   for (int imol=0; imol<n_molecules(); imol++) { 
      if (molecules[imol].has_map()) { 
	 nmaps++;
	 only_map = imol;
      } 
   }
   if (nmaps == 1) {
      // let's set it then (trial code)
      imol_refinement_map = only_map;
      return only_map;
   } else { 
      return -1;
   }
}

int
graphics_info_t::set_imol_refinement_map(int imol) {

   int r = -1;
   if (molecules[imol].has_map()) { 
      imol_refinement_map = imol;
      r = imol;
   }
   return r;
} 



void
graphics_info_t::add_vector_to_RotationCentre(const coot::Cartesian &vec) { 
   
   rotation_centre_x += vec.x();
   rotation_centre_y += vec.y();
   rotation_centre_z += vec.z();

   if (GetActiveMapDrag() == 1) {
      for (int ii=0; ii<n_molecules(); ii++) { 
	 if (molecules[ii].has_map()) { 
	    molecules[ii].update_map(); // to take account
	                                // of new rotation centre.
	 }
      }
   }
   for (int ii=0; ii<n_molecules(); ii++) { 
      molecules[ii].update_symmetry();
   }
   graphics_draw();
}



// return NULL on not found:
CAtom *
graphics_info_t::find_atom_in_moving_atoms(const coot::atom_spec_t &at) const {

   CAtom *cat = NULL;
   if (moving_atoms_asc->mol != NULL) { 
      int SelHnd = coot::get_selection_handle(moving_atoms_asc->mol, at);
      int nSelAtoms; 
      PPCAtom local_SelAtom = NULL; 
      moving_atoms_asc->mol->GetSelIndex(SelHnd, local_SelAtom, nSelAtoms);
      if (nSelAtoms > 0)
	 cat = local_SelAtom[0];
//       std::cout << "DEBUG:: in find_atom_in_moving_atoms: here are the "
// 		<< nSelAtoms << " qualifying atoms..." << std::endl;
//       for(int i=0; i<nSelAtoms; i++)
// 	 std::cout << "      " << i << "  " << local_SelAtom[i] << std::endl;
      moving_atoms_asc->mol->DeleteSelection(SelHnd);
   }
   return cat;
}

// Rotate position round direction
// 
clipper::Coord_orth
graphics_info_t::rotate_round_vector(const clipper::Coord_orth &direction,
				     const clipper::Coord_orth &position,
				     const clipper::Coord_orth &origin_shift,
				     double angle) const {

   // moved this function down to utils
   return coot::util::rotate_round_vector(direction, position, origin_shift, angle);
}



int
graphics_info_t::load_db_main() { 
   
   return 1;

}


// --------------------------------------------------------------------------------
//                 residue info widget
// --------------------------------------------------------------------------------


void
graphics_info_t::fill_output_residue_info_widget(GtkWidget *widget, int imol, PPCAtom atoms, int n_atoms) {

   // first do the label of the dialog
   GtkWidget *label_widget = lookup_widget(widget, "residue_info_residue_label");

   std::string label = "Molecule: ";
   label += int_to_string(imol);
   label += " ";
   label += molecules[imol].name_;

   gtk_label_set_text(GTK_LABEL(label_widget), label.c_str());
   GtkWidget *table = lookup_widget(widget, "residue_info_atom_table");
   
   residue_info_n_atoms = n_atoms; 
   for (int i=0; i<n_atoms; i++)
      graphics_info_t::fill_output_residue_info_widget_atom(table, imol, atoms[i], i);
}

void
graphics_info_t::fill_output_residue_info_widget_atom(GtkWidget *table, int imol, PCAtom atom,
						      int iatom) {

   GtkWidget *residue_info_dialog = lookup_widget(table, "residue_info_dialog");
   gint left_attach = 0;
   gint right_attach = 1;
   gint top_attach = iatom;
   gint bottom_attach = top_attach + 1;
   GtkAttachOptions xopt = GTK_FILL, yopt=GTK_FILL;
   guint xpad=0, ypad=0;

   // The text label of the atom name:
   left_attach = 0;
   right_attach = left_attach + 1;
   std::string label_str = "  ";
   label_str += atom->GetChainID();
   label_str += "/";
   label_str += graphics_info_t::int_to_string(atom->GetSeqNum());
   if (std::string(atom->GetInsCode()) != "") {
      label_str += atom->GetInsCode();
   }
   label_str += " ";
   label_str += atom->GetResName();
   label_str += "/";
   label_str += atom->name;
   if (std::string(atom->altLoc) != std::string("")) {
      label_str += ",";
      label_str += atom->altLoc;
   }
   if (std::string(atom->segID) != std::string("")) {
      label_str += " ";
      label_str += atom->segID;
   }
   label_str += "  ";

   GtkWidget *residue_info_atom_info_label = gtk_label_new (label_str.c_str());
   gtk_table_attach(GTK_TABLE(table), residue_info_atom_info_label,
		    left_attach, right_attach, top_attach, bottom_attach,
		    xopt, yopt, xpad, ypad);
   gtk_widget_ref (residue_info_atom_info_label);
   gtk_object_set_data_full (GTK_OBJECT (residue_info_dialog),
			     "residue_info_atom_info_label", residue_info_atom_info_label,
			     (GtkDestroyNotify) gtk_widget_unref);
   gtk_widget_show (residue_info_atom_info_label);


   // The Occupancy entry:
   left_attach = 1;
   right_attach = left_attach + 1;
   coot::select_atom_info *ai = new coot::select_atom_info;
   *ai = coot::select_atom_info(iatom, imol, 
				std::string(atom->GetChainID()), 
				atom->GetSeqNum(),
				std::string(atom->GetInsCode()),
				std::string(atom->name),
				std::string(atom->altLoc));

   std::string widget_name = "residue_info_occ_entry_";
   widget_name += int_to_string(iatom);
   // 
   GtkWidget *residue_info_occ_entry = gtk_entry_new ();
   gtk_widget_ref (residue_info_occ_entry);
   gtk_object_set_data_full (GTK_OBJECT (residue_info_dialog),
			     widget_name.c_str(), residue_info_occ_entry,
			     (GtkDestroyNotify) gtk_widget_unref);
   gtk_widget_show (residue_info_occ_entry);
   gtk_object_set_user_data(GTK_OBJECT(residue_info_occ_entry), ai);
   gtk_entry_set_text(GTK_ENTRY(residue_info_occ_entry),
		      graphics_info_t::float_to_string(atom->occupancy).c_str());
   gtk_table_attach(GTK_TABLE(table), residue_info_occ_entry,
		    left_attach, right_attach, top_attach, bottom_attach,
		    xopt, yopt, xpad, ypad);


      // Note that we have to use key_release_event because if we use
   // key_press_event, when we try to get the value from the widget
   // (gtk_entry_get_text) then that does not see the key/number that
   // was just pressed.

   // B-factor entry:
   left_attach = 2;
   right_attach = left_attach + 1;

   widget_name = "residue_info_b_factor_entry_";
   widget_name += int_to_string(iatom);

   GtkWidget *residue_info_b_factor_entry = gtk_entry_new ();
   gtk_widget_ref (residue_info_b_factor_entry);
   gtk_object_set_data_full (GTK_OBJECT (residue_info_dialog),
			     widget_name.c_str(), residue_info_b_factor_entry,
			     (GtkDestroyNotify) gtk_widget_unref);
   gtk_widget_show (residue_info_b_factor_entry);
   gtk_entry_set_text(GTK_ENTRY(residue_info_b_factor_entry),
		      graphics_info_t::float_to_string(atom->tempFactor).c_str());
   // gtk_widget_set_usize(residue_info_b_factor_entry, 40, -2);
   gtk_object_set_user_data(GTK_OBJECT(residue_info_b_factor_entry), ai);
   gtk_widget_set_events(residue_info_b_factor_entry, 
			 GDK_KEY_PRESS_MASK     |
			 GDK_KEY_RELEASE_MASK);
   gtk_table_attach(GTK_TABLE(table), residue_info_b_factor_entry,
		    left_attach, right_attach, top_attach, bottom_attach,
		    xopt, yopt, xpad, ypad);

}

// static 
gboolean
graphics_info_t::on_residue_info_master_atom_occ_changed (GtkWidget       *widget,
							  GdkEventKey     *event,
							  gpointer         user_data) {
   const gchar *s = gtk_entry_get_text(GTK_ENTRY(widget));
   if (s) { 
      // consider strtof:
      // 
      // double f = atof(s);
      graphics_info_t::residue_info_pending_edit_occ = 1;
      
      graphics_info_t g;
      g.residue_info_edit_occ_apply_to_other_entries_maybe(widget);
   }
   return TRUE; 
}

// static
gboolean
graphics_info_t::on_residue_info_master_atom_b_factor_changed (GtkWidget       *widget,
							       GdkEventKey     *event,
							       gpointer         user_data) {

   // Let's get the entry value:
   //
   const gchar *s = gtk_entry_get_text(GTK_ENTRY(widget));
   if (s) { 
      // consider strtof:
      // 
      // float f = atof(s);
      graphics_info_t::residue_info_pending_edit_b_factor = 1;
      
      graphics_info_t g;
      g.residue_info_edit_b_factor_apply_to_other_entries_maybe(widget);
   }
   return TRUE;
}




//static 
void
graphics_info_t::residue_info_edit_b_factor_apply_to_other_entries_maybe(GtkWidget *start_entry) {

   // first find the checkbox:
   GtkWidget *dialog = lookup_widget(start_entry, "residue_info_dialog"); 
   GtkWidget *checkbutton = lookup_widget(dialog, "residue_info_b_factor_apply_all_checkbutton");
   std::string widget_name;
   GtkWidget *entry;

   if (! checkbutton) { 
      std::cout << "ERROR:: could not find checkbutton" << std::endl; 
   } else { 
      if (GTK_TOGGLE_BUTTON(checkbutton)->active) { 
	 
	 // propogate the change to the other b-factor widgets
	 std::string entry_text(gtk_entry_get_text(GTK_ENTRY(start_entry)));
	 for (int i=0; i<graphics_info_t::residue_info_n_atoms; i++) {
	    widget_name = "residue_info_b_factor_entry_";
	    widget_name += int_to_string(i);
	    entry = lookup_widget(dialog, widget_name.c_str());
	    if (entry) { 
	       gtk_entry_set_text(GTK_ENTRY(entry), entry_text.c_str());
	    } else { 
	       std::cout << "ERROR: no entry\n";
	    } 
	 }
      } 
   }
}


//static 
void
graphics_info_t::residue_info_edit_occ_apply_to_other_entries_maybe(GtkWidget *start_entry) {

   // first find the checkbox:
   GtkWidget *dialog = lookup_widget(start_entry, "residue_info_dialog"); 
   GtkWidget *checkbutton = lookup_widget(dialog, "residue_info_occ_apply_all_checkbutton");
   std::string widget_name;
   GtkWidget *entry;

   if (! checkbutton) { 
      std::cout << "ERROR:: could not find checkbutton" << std::endl; 
   } else { 
      if (GTK_TOGGLE_BUTTON(checkbutton)->active) { 
	 
	 // propogate the change to the other b-factor widgets
	 std::string entry_text(gtk_entry_get_text(GTK_ENTRY(start_entry)));
	 for (int i=0; i<graphics_info_t::residue_info_n_atoms; i++) {
	    widget_name = "residue_info_occ_entry_";
	    widget_name += int_to_string(i);
	    entry = lookup_widget(dialog, widget_name.c_str());
	    if (entry) { 
	       gtk_entry_set_text(GTK_ENTRY(entry), entry_text.c_str());
	    } else { 
	       std::cout << "ERROR: no entry\n";
	    } 
	 }
      } 
   }
}


// static
void
graphics_info_t::residue_info_add_b_factor_edit(coot::select_atom_info sai, 
						float val) { 

   graphics_info_t g;
   short int made_substitution_flag = 0;
   for (unsigned int i=0; i<g.residue_info_edits->size(); i++) {
      if (sai.udd == (*g.residue_info_edits)[i].udd) {
	 (*g.residue_info_edits)[i].add_b_factor_edit(val);
	 made_substitution_flag = 1;
	 break;
      }
   }
   if (! made_substitution_flag) { 
      sai.add_b_factor_edit(val);
      g.residue_info_edits->push_back(sai);
   }
}

// static
void
graphics_info_t::residue_info_add_occ_edit(coot::select_atom_info sai, 
					   float val) { 

   graphics_info_t g;
   short int made_substitution_flag = 0;
   for (unsigned int i=0; i<g.residue_info_edits->size(); i++) {
      if (sai.udd == (*g.residue_info_edits)[i].udd) {
	 (*g.residue_info_edits)[i].add_occ_edit(val);
	 made_substitution_flag = 1;
	 break;
      }
   }
   if (! made_substitution_flag) { 
      sai.add_occ_edit(val);
      g.residue_info_edits->push_back(sai);
   }
}

// This is the callback when the OK button of the residue info was pressed.
//
// The new way with a table:
void
graphics_info_t::apply_residue_info_changes(GtkWidget *dialog) {


   int imol = -1;
   // This is where we accumulate the residue edits:
   std::vector<coot::select_atom_info> local_atom_edits;

   GtkWidget *table = lookup_widget(dialog, "residue_info_atom_table");

#if (GTK_MAJOR_VERSION > 1)
   GList *container_list = gtk_container_get_children(GTK_CONTAINER(table));
#else    
   GList *container_list = gtk_container_children(GTK_CONTAINER(table));
#endif
   
   // The children are a list, gone in "backward", just like we'd been
   // consing onto a list as we added widgets to the table.
   // 
   int len = g_list_length(container_list);
   for(int i=0; i < len; i+=3) {
      if (i+1<len) { 
	 GtkWidget *widget_b = (GtkWidget*) g_list_nth_data(container_list, i);
	 GtkWidget *widget_o = (GtkWidget*) g_list_nth_data(container_list, i+1);
	 std::string b_text = gtk_entry_get_text(GTK_ENTRY(widget_b));
	 std::string o_text = gtk_entry_get_text(GTK_ENTRY(widget_o));
// 	 std::cout << "b_text :" <<b_text << std::endl;
// 	 std::cout << "o_text :" <<o_text << std::endl;

	 // Handle OCCUPANCY edits
	 // 
	 coot::select_atom_info *ai =
	    (coot::select_atom_info *) gtk_object_get_user_data(GTK_OBJECT(widget_o));
	 if (ai) {
	    imol = ai->molecule_number;  // hehe
	    CAtom *at = ai->get_atom(graphics_info_t::molecules[imol].atom_sel.mol);
	    // std::cout << "got atom at " << at << std::endl;
	    std::pair<short int, float>  occ_entry =
	       graphics_info_t::float_from_entry(GTK_WIDGET(widget_o));
	    if (occ_entry.first) {
	       if (at) { 
// 		  std::cout << "    occ comparison " << occ_entry.second << " "
// 		  << at->occupancy << std::endl;
		  if (abs(occ_entry.second - at->occupancy) > 0.009) {
		     coot::select_atom_info local_at = *ai;
		     local_at.add_occ_edit(occ_entry.second);
		     local_atom_edits.push_back(local_at);
		  } 
	       }
	    }
	 } else {
	    std::cout << "no user data found for widget_o" << std::endl;
	 }

	 // HANDLE B-FACTOR edits
	 ai = (coot::select_atom_info *) gtk_object_get_user_data(GTK_OBJECT(widget_b));
	 if (ai) {
	    imol = ai->molecule_number;  // hehe
	    CAtom *at = ai->get_atom(graphics_info_t::molecules[imol].atom_sel.mol);
	    std::pair<short int, float>  temp_entry =
	       graphics_info_t::float_from_entry(GTK_WIDGET(widget_b));
	    if (temp_entry.first) {
	       if (at) { 
		  // std::cout << "    temp comparison " << temp_entry.second
		  // << " " << at->tempFactor << std::endl;
		  if (abs(temp_entry.second - at->tempFactor) > 0.009) {
		     coot::select_atom_info local_at = *ai;
		     local_at.add_b_factor_edit(temp_entry.second);
		     local_atom_edits.push_back(local_at);
		  }
	       }
	    }
	 }
	 
      } else {
	 std::cout << "Programmer error in decoding table." << std::endl;
	 std::cout << "  Residue Edits not applied!" << std::endl;
      }
   }

   if (local_atom_edits.size() >0)
      if (imol >= 0)
	 graphics_info_t::molecules[imol].apply_atom_edits(local_atom_edits);

   residue_info_edits->clear();
   // delete res_spec; // can't do this: the user may press the button twice
}

// static
// (should be called by the destroy event and the close button)
void
graphics_info_t::residue_info_release_memory(GtkWidget *dialog) { 

   GtkWidget *entry;
   for (int i=0; i<residue_info_n_atoms; i++) { 
      std::string widget_name = "residue_info_b_factor_entry_";
      widget_name += int_to_string(i);
      entry = lookup_widget(dialog, widget_name.c_str());
      if (entry) { 
	 coot::select_atom_info *sai_p = (coot::select_atom_info *) gtk_object_get_user_data(GTK_OBJECT(entry));
	 if (sai_p) { 
	    // delete sai_p; // memory bug.  We cant do this
	    // std::cout << "hmmm.. not deleting sai_p" << std::endl;
	 } else { 
	    std::cout << "ERROR:: no user data in b-factor entry widget\n"; 
	 } 
      }
      // same for occ entry:
      widget_name = "residue_info_occ_entry_";
      widget_name += int_to_string(i);
      entry = lookup_widget(dialog, widget_name.c_str());
      if (entry) { 
	 coot::select_atom_info *sai_p = (coot::select_atom_info *) gtk_object_get_user_data(GTK_OBJECT(entry));
	 if (sai_p) { 
	    // std::cout << "not deleting sai_p" << std::endl;
	    // delete sai_p; // memory bug.  We cant do this
	 } else { 
	    std::cout << "ERROR:: no user data in occ entry widget\n"; 
	 }
      }
   } 
} 



// 
int
graphics_info_t::pointer_atom_molecule() const {

   int i = -1; // must be changed by this function.

   if (user_pointer_atom_molecule >= 0) {
      if (user_pointer_atom_molecule < n_molecules()) {
	 if (molecules[user_pointer_atom_molecule].open_molecule_p()) {
	    i = user_pointer_atom_molecule;
	 }
      }
   }
      
   if (i == -1) { 
      // user did not explicictly set the molecule, the usual case:

      for (int imol=0; imol<n_molecules(); imol++) {
	 if (molecules[imol].open_molecule_p()) { // not closed
	    if (molecules[imol].name_ == "Pointer Atoms") {
	       return imol;
	    }
	 }
      }
	     
      // If we get here, it was not found, let's create one:
	 
      std::cout << "Creating a molecule for Pointer Atoms" << std::endl;
      MyCMMDBManager *MMDBManager = new MyCMMDBManager();

      // do we attach a model, chain etc here?
      // Yes.
      CModel *model_p = new CModel;
      CChain *chain_p = new CChain;
   
      model_p->AddChain(chain_p);
      MMDBManager->AddModel(model_p);
   
      atom_selection_container_t asc = make_asc(MMDBManager);
      int imol = create_molecule();
      molecules[imol].install_model(imol, asc, "Pointer Atoms", 1);
      return imol;
   }
   return i;
}

void
graphics_info_t::update_things_on_move_and_redraw() {

   update_things_on_move();
   graphics_draw();
}


void
graphics_info_t::update_things_on_move() {
   
   for (int ii=0; ii<n_molecules(); ii++) { 
      molecules[ii].update_map();
      molecules[ii].update_clipper_skeleton(); 
      molecules[ii].update_symmetry();
   }
} 

void
graphics_info_t::start_baton_here() {

   baton_root = RotationCentre();
   
   int imol_for_skel = imol_for_skeleton(); // if unset, sets and
					    // returns if only one
					    // map, else return -1
   if (imol_for_skel < 0) {

      std::cout << "WARNING: no skeleton found " << std::endl;

      GtkWidget *w = create_baton_mode_make_skeleton_dialog();
      int *imol_copy = new int;
      *imol_copy = imol_for_skel;
      gtk_object_set_user_data(GTK_OBJECT(w), (char *)imol_copy);
      gtk_widget_show(w);

   } else {

      molecules[imol_for_skel].fill_skeleton_treenodemap(); // filled only if not filled.
      clipper::Coord_grid dummy_cg;
      short int use_dummy_cg = 0;
      baton_next_directions(imol_for_skel, NULL, baton_root, dummy_cg, use_dummy_cg);
      baton_next_ca_options_index = 0; 
      baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
   }
}

void
graphics_info_t::baton_next_directions(int imol_for_skel, const CAtom *latest_atom,
				       const coot::Cartesian &baton_root,
				       const clipper::Coord_grid &cg_start,
				       short int use_cg_start) {


   std::cout << "DEBUG in baton_next_directions imol_for_skel is "
	     << imol_for_skel << std::endl;
         
   std::vector<clipper::Coord_orth> previous_ca_positions;
   // store the position of the just accepted atom as a previous atom
   // 
   // previous_ca_positions.push_back(to_coord_orth(baton_root));
   int imol_baton_atoms = baton_build_atoms_molecule();


   // std::cout << "DEBUG INFO:::::: latest_atom is " << latest_atom << std::endl;

   if (latest_atom == NULL) {
      previous_ca_positions.push_back(to_coord_orth(baton_root));
   } else {
      previous_ca_positions = molecules[imol_baton_atoms].previous_baton_atom(latest_atom, 
									      baton_build_direction_flag);
   }
//    std::cout << "DEBUG: in graphics: cg_start: " << cg_start.format() << "  "
// 	     << use_cg_start << std::endl;
   *baton_next_ca_options = molecules[imol_for_skel].next_ca_by_skel(previous_ca_positions,
								     cg_start,
								     use_cg_start,
								     3.8,
								     skeleton_level, 
								     max_skeleton_search_depth);

   // Print out the baton_next_ca_options
   // 
   std::cout << "-- baton_next_ca_options" << std::endl;
   for(unsigned int i=0; i<baton_next_ca_options->size(); i++) {
      std::cout << "   " << (*baton_next_ca_options)[i].score  << "  "
		<< (*baton_next_ca_options)[i].position.format() << std::endl;
   }
   std::cout << "--" << std::endl;

   // Graphics the baton_next_ca_options
   // 
   std::string molname("Baton Atom Guide Points");
   if (baton_tmp_atoms_to_new_molecule)
      create_molecule_and_display(*baton_next_ca_options, molname);
   else 
      update_molecule_to(*baton_next_ca_options, molname);

}



void
graphics_info_t::baton_object() {

//    std::cout << "baton from " << baton_root << " to " << baton_tip
// 	     << " draw_baton_flag: " << draw_baton_flag << std::endl;
   
   if (graphics_info_t::draw_baton_flag) {
      glColor3f (0.8, 0.8, 0.9);
      glBegin(GL_LINES);
      glVertex3f(graphics_info_t::baton_root.x(),
		 graphics_info_t::baton_root.y(),
		 graphics_info_t::baton_root.z());
      glVertex3f( graphics_info_t::baton_tip.x(),
		  graphics_info_t::baton_tip.y(),
		  graphics_info_t::baton_tip.z());
      glEnd();
   }
   
} 

// aka imol_for_baton_atoms()
int
graphics_info_t::baton_build_atoms_molecule() const {

   int imol = -1;
   
   for (int im=0; im<n_molecules(); im++) {
      if (molecules[im].name_ == "Baton Atoms") {
	 return im;
      }
   }

   
   std::cout << "Creating a molecule for Baton Atoms" << std::endl;
   // not found, let's create one:
   MyCMMDBManager *MMDBManager = new MyCMMDBManager();
   

   // do we attach a model, chain etc here?
   // Yes.
   CModel *model_p = new CModel;
   CChain *chain_p = new CChain;
   chain_p->SetChainID(baton_build_chain_id.c_str());
   
   model_p->AddChain(chain_p);
   MMDBManager->AddModel(model_p);

   // now lets add an mmdbcryst to the mmdbmanger, which is generated
   // from the skeleton map.
   //
   int  imol_for_skel = imol_for_skeleton();
   std::cout << "DEBUG:: imol_for_skeleton in baton_build_atoms_molecule is " << imol_for_skel << std::endl;
   if (imol_for_skel >= 0) {
      // CMMDBCryst *cryst = new CMMDBCryst;
      MMDBManager->SetCell(molecules[imol_for_skel].xskel_cowtan.cell().descr().a(),
			   molecules[imol_for_skel].xskel_cowtan.cell().descr().b(),
			   molecules[imol_for_skel].xskel_cowtan.cell().descr().c(),
			   clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().alpha()),
			   clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().beta()),
			   clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().gamma()), 1);

      std::string spacegroup = molecules[imol_for_skel].xskel_cowtan.spacegroup().symbol_hm();
      std::cout << "DEBUG got spacegroup " << spacegroup << std::endl;

      std::cout << "setting spacegroup of Baton Atoms to be: " << spacegroup << std::endl;
      std::cout << "setting cell of Baton Atoms to be: "
		<< molecules[imol_for_skel].xskel_cowtan.cell().format() << std::endl;
      
      int istat_spgr = MMDBManager->SetSpaceGroup((char *)spacegroup.c_str()); // bleugh!
      std::cout << "DEBUG:: status from SetSpaceGroup: " << istat_spgr << std::endl;
      if (istat_spgr != 0) {
	 std::cout << "Problem:: mmdb does not understand space group: " << spacegroup << std::endl;
      }  
      
   } else {
      std::cout << "WARNING: skeleton not found - no symmetry for Baton Atoms " << std::endl;
   }
   
   atom_selection_container_t asc = make_asc(MMDBManager);
   asc.SelectionHandle = -1;
   imol = create_molecule();
   molecules[imol].install_model(imol, asc, "Baton Atoms", 1);
   std::cout << "INFO:: Baton Build molecule has SelectionHandle: " 
	     << molecules[imol].atom_sel.SelectionHandle << " " 
	     <<  asc.SelectionHandle << std::endl;
   return imol;
} 

void
graphics_info_t::accept_baton_position() {

   int imol_for_skel = imol_for_skeleton();

   // First add the atom to the baton build molecule
   //
   //
   CAtom *baton_atom = NULL; // for baton_next_directions() usage?
   int imol = baton_build_atoms_molecule();
   if (imol >= 0) {
      baton_atom = molecules[imol].add_baton_atom(baton_tip, 
						  baton_build_start_resno, 
						  baton_build_params_active,
						  baton_build_direction_flag);
      baton_build_params_active = 0; // This flag was set after
				     // set_baton_build_params.  We
				     // clear it now so that we don't
				     // any more force the start resno -
				     // molecule_class_info_t::add_baton_atom
				     // can work it out.
   }
   std::cout << "setting screen rotation centre to " << baton_tip << std::endl; 
   setRotationCentre(baton_tip);
   for(int ii=0; ii<n_molecules(); ii++) {
      // but not skeleton, lets do skeleton only on a middle-mouse recentre
      molecules[ii].update_map();
      molecules[ii].update_symmetry();
   }

   // debug
   // 
   // std::cout << "DEBUG:: in accept_baton_position baton_atom is " 
   // << baton_atom << std::endl;

   // But first show us the options for the next point as dummy atoms
   // 

   // int imol_for_skel = imol_for_skeleton();
   // internally set the baton_next_ca_options
   //
   if (imol_for_skel < 0) {
      std::cout << "Ooops:: must have a skeleton first" << std::endl;
   } else {
      short int use_cg = 1;
      clipper::Coord_grid cg = (*baton_next_ca_options)[baton_next_ca_options_index].near_grid_pos;
      baton_next_directions(imol_for_skel, baton_atom, baton_tip, cg, use_cg); // old tip
   }
						   
   // Now set the baton tip to the next point:

   
   // set these for redraw
   baton_root = baton_tip;
   baton_next_ca_options_index = 0; // for next
   baton_length = 3.8; // reset to to optimal

   baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
      
   graphics_draw();

}



coot::Cartesian
graphics_info_t::baton_tip_by_ca_option(int index) const {

   coot::Cartesian tip_pos(0.0, 0.0, 0.0);
   unsigned int uindex = index;

   if (!baton_next_ca_options) {
      std::cout << "ERROR: baton_next_ca_options is NULL\n";
   } else { 
      if (uindex >= baton_next_ca_options->size()) {
	 if ((uindex == 0) && (baton_next_ca_options->size() == 0)) {
	    std::cout << "INFO:: no baton next positions from here\n";
	 } else { 
	    std::cout << "ERROR: bad baton_next_ca_options index: "
		      << index << " size " << baton_next_ca_options->size()
		      << std::endl;
	 }
      } else {
	 // now we want a vector baton_length in the direction starting
	 // at baton_root to baton_next_ca_options[index]
	 //
	 coot::Cartesian target_point = to_cartesian((*baton_next_ca_options)[index].position);
	 std::cout << "Ca option " << index << " score: "
		   << (*baton_next_ca_options)[index].score << std::endl;
	 
	 coot::Cartesian target_dir = target_point - baton_root;

	 target_dir.unit_vector_yourself();

	 target_dir *= baton_length;

	 tip_pos = target_dir + baton_root;
      }
   }
   return tip_pos;
} 


// Gawd, Kevin's having an effect on me.. (my programming at least...)
// "const double &" indeed?  Tsk! Whatever next?
// 
void
graphics_info_t::rotate_baton(const double &x, const double &y) {

   mouse_current_x = x;
   mouse_current_y = y;
   double diff;

   diff  = mouse_current_x - GetMouseBeginX();
   diff += mouse_current_y - GetMouseBeginY();
   
   coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
   coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
   coot::Cartesian screen_z = (front - centre);

   clipper::Coord_orth new_pos = rotate_round_vector(to_coord_orth(screen_z),
						     to_coord_orth(baton_tip),
						     to_coord_orth(baton_root),
						     0.01*diff);

   baton_tip = to_cartesian(new_pos);
   graphics_draw();
} 

void
graphics_info_t::toggle_baton_mode() {

   if (baton_mode == 0) { 
      baton_mode = 1;
      std::cout << "INFO::baton rotation mode on." << std::endl;
   } else {
      baton_mode = 0;
      std::cout << "INFO::baton rotation mode off." << std::endl;
   }
}

void
graphics_info_t::baton_try_another() {

   baton_next_ca_options_index++;
   
   // make baton_next_ca_options_index an unsigned int
   if (baton_next_ca_options_index >= int(baton_next_ca_options->size())) {
      std::cout << "info: cycling back to start of ca options" << std::endl;
      baton_next_ca_options_index = 0;
   }
   baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
   graphics_draw();
}

void
graphics_info_t::shorten_baton() {

   double short_factor = 0.952;
   baton_length *= short_factor;
//    baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
   coot::Cartesian baton_vec = baton_tip - baton_root;
   baton_vec *= short_factor;
   baton_tip = baton_root + baton_vec;
   graphics_draw();
}


void
graphics_info_t::lengthen_baton() {

   double lengthen_factor = 1.05;
   baton_length *= lengthen_factor;
   coot::Cartesian baton_vec = baton_tip - baton_root;
   baton_vec *= lengthen_factor;
   baton_tip = baton_root + baton_vec;
   graphics_draw();
}

void
graphics_info_t::baton_build_delete_last_residue() {

   int imol = baton_build_atoms_molecule();
   if (imol > -1) {
      std::pair<short int, CAtom *> new_centre
	 = molecules[imol].baton_build_delete_last_residue();
      if (new_centre.first) {
	 coot::Cartesian new_centre_cart(new_centre.second->x,
				   new_centre.second->y,
				   new_centre.second->z);
	 setRotationCentre(new_centre_cart);
	 baton_root = new_centre_cart;
	 int imol_for_skel = imol_for_skeleton();
	 if (imol_for_skel >= 0) {
	    short int use_cg = 1;
	    // recall std::vector<coot::scored_skel_coord> *baton_next_ca_options;
	    // 
	    std::pair <short int,clipper::Coord_grid> cg =
	       molecules[imol_for_skel].search_for_skeleton_near(new_centre_cart);

	    if (cg.first)
	       baton_next_directions(imol_for_skel, new_centre.second, new_centre_cart, cg.second, use_cg);
	 } 
	 baton_next_ca_options_index = 0;
	 baton_length = 3.8;
	 baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
	 graphics_draw();
      }
   } 
} 


// should be similar to Imol for refinement?
int
graphics_info_t::imol_for_skeleton() const {

//    for (int imol=0; imol<n_molecules; imol++) {
//       if (molecules[imol].xskel_is_filled) {
// 	 return imol;
//       }
//    }

//    int nmaps = 0;
//    int only_map;

//     if (map_for_skeletonize == -1) { 
//        for (int imol=0; imol<n_molecules; imol++) {
//  	 if (molecules[imol].has_map()) {
//  	    nmaps++;
//  	 only_map = imol;
//  	 }
//        }
//        if (nmaps == 1) {
//  	 map_for_skeletonize = only_map;
//        }
//     }
   
   return map_for_skeletonize;
} 

;
void
graphics_info_t::create_molecule_and_display(std::vector<coot::scored_skel_coord> &pos_position,
					     const std::string &molname) {

   int imol = create_empty_molecule(molname);
   std::vector<coot::Cartesian> cv;
   // now add atoms:
   for (unsigned int i=0; i<pos_position.size(); i++) {
      coot::Cartesian c(pos_position[i].position.x(),
		  pos_position[i].position.y(),
		  pos_position[i].position.z());
      cv.push_back(c);
   }
   molecules[imol].add_multiple_dummies(cv);

}

// as above, except we update molecule with name molname to the
// pos_positions (and delete everything else).
// 
void
graphics_info_t::update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position,
					  const std::string &molname) {

   int imol = lookup_molecule_name(molname);

   if (imol >= 0) {
      if (pos_position.size() > 0) 
	 graphics_info_t::molecules[imol].update_molecule_to(pos_position);
      else
	 std::cout << "WARNING:: No atoms guide points in update_molecule_to."
		   << "  Not updating guide points molecule" << std::endl;
   } else {
      create_molecule_and_display(pos_position, molname);
   } 
}

// return -1 on error
int
graphics_info_t::lookup_molecule_name(const std::string &molname) const {
   
   for (int imol=0; imol<n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].name_ == molname) {
	 return imol;
      }
   }
   return -1;
}


  
int
graphics_info_t::create_empty_molecule(const std::string &molname) {

   std::cout << "Creating a molecule for " << molname << std::endl;

   MyCMMDBManager *MMDBManager = new MyCMMDBManager();

   CModel *model_p = new CModel;
   CChain *chain_p = new CChain;
   
   model_p->AddChain(chain_p);
   MMDBManager->AddModel(model_p);
   
   atom_selection_container_t asc = make_asc(MMDBManager);
   int imol = create_molecule();
   molecules[imol].install_model(imol, asc, molname, 1);
   asc.read_error_message = "No error";
   asc.read_success = 1;
   return imol;
} 


// ------------------------------------------------------------------
//                        undo functions
// ------------------------------------------------------------------

// This is the callback when "Undo" is pressed:
//
// There is a problem now that we are including redo code: the
// question is what to do when we reach the end of the modifications
// in the current undo_molecule.  This is hard... let's go with the
// current functionality for the moment, which is to unset it when we
// get to the end of the modifications list.
//
// But we need to modify the response of Undo_molecule() depending on
// whether it was asked by apply_undo or apply_redo, we do that by
// passing an enumerated type.
// 
// 
void
graphics_info_t::apply_undo() {

   int umol = Undo_molecule(coot::UNDO);
   // std::cout << "DEBUG:: undo molecule : " << umol << std::endl;
   if (umol == -2) {
      if (use_graphics_interface_flag) { 
	 GtkWidget *dialog = create_undo_molecule_chooser_dialog();
	 GtkWidget *option_menu = lookup_widget(dialog,
						"undo_molecule_chooser_option_menu");
	 fill_option_menu_with_undo_options(option_menu);
	 gtk_widget_show(dialog);
      }
   } else {
      if (umol == -1) {
	 std::cout << "There are no molecules with modifications "
		   << "that can be undone" << std::endl;
      } else {

	 if (molecules[umol].Have_modifications_p()) { 
	    if (molecules[umol].is_displayed_p()) { 
	       molecules[umol].apply_undo();
	       if (use_graphics_interface_flag) { 
		  graphics_draw();
		  
		  // need to update the atom and residue list in Go To Atom widget
		  // (maybe)
		  update_go_to_atom_window_on_changed_mol(umol);
		  
		  // update the ramachandran, if there was one
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
		  GtkWidget *w = coot::get_validation_graph(umol, coot::RAMACHANDRAN_PLOT);
		  if (w) {
		     coot::rama_plot *plot = (coot::rama_plot *)
			gtk_object_get_user_data(GTK_OBJECT(w));
		     handle_rama_plot_update(plot);
		  }
#endif // HAVE_GTK_CANVAS   
	       }
	    } else {
	       if (use_graphics_interface_flag) { 
		  std::string s = "WARNING:: Coot will not undo modifications on a \n";
		  s += "molecule that is not displayed";
		  GtkWidget *w = wrapped_nothing_bad_dialog(s);
		  gtk_widget_show(w);
	       }
	    }
	 } else {
	    if (use_graphics_interface_flag) { 
	       undo_molecule = -1; // reset it
	       std::cout << "WARNING:: !!!  Changing the molecule to which "
			 << "\"Undo\"s are done." << std::endl;
	       std::string s = "WARNING:: Changing to Undo molecule";
	       statusbar_text(s);
	    }
	    apply_undo();       // find another molecule to undo
	 }
      }
   }

   // and now tinker with the Redo button to make it active
   //
   activate_redo_button();  // has protection for --no-graphics
}

void
graphics_info_t::apply_redo() { 
   
   int umol = Undo_molecule(coot::REDO);
   if (umol == -2) { // ambiguity
      GtkWidget *dialog = create_undo_molecule_chooser_dialog();
      GtkWidget *option_menu = lookup_widget(dialog,
					     "undo_molecule_chooser_option_menu");
      fill_option_menu_with_undo_options(option_menu);
      gtk_widget_show(dialog);
   } else {
      if (umol == -1) { // unset
	 std::cout << "There are no molecules with modifications "
		   << "that can be re-done" << std::endl;
      } else {

	 if (molecules[umol].Have_redoable_modifications_p()) {
	    // std::cout << "DEBUG:: applying redo" << std::endl;
	    molecules[umol].apply_redo();
	    graphics_draw();
	    
	    // need to update the atom and residue list in Go To Atom widget
	    // (maybe)
	    update_go_to_atom_window_on_changed_mol(umol);
	 } else {
	    // std::cout << "DEBUG:: not applying redo" << std::endl;
	 }
      }
   }
}

// This is a noddy - better is needed.
//
// No account is taken to deactivate the button again when we have run
// out of redos.
//
// When the widget gets created, the redo button is always insensitve
// - it should depend on if there are redoable molecules.
// 
void
graphics_info_t::activate_redo_button() {

   GtkWidget *dialog = model_fit_refine_dialog;

   
   if (dialog) { 
      // which it should be!
      GtkWidget *button = lookup_widget(dialog, "model_refine_dialog_redo_button");
      gtk_widget_set_sensitive(button, TRUE);
   }
} 


// static 
void
graphics_info_t::pointer_atom_molecule_menu_item_activate(GtkWidget *item, 
							  GtkPositionType pos) {

   graphics_info_t g;
   //    std::cout << "DEBUG:: pointer_atom_molecule_menu_item_activate sets user_pointer_atom_molecule to " << pos << std::endl;
   g.user_pointer_atom_molecule = pos;

}


   
// We are passed an GtkOptionMenu *option_menu
//
void
graphics_info_t::fill_option_menu_with_coordinates_options(GtkWidget *option_menu,
							   GtkSignalFunc callback_func) { 

   fill_option_menu_with_coordinates_options_internal(option_menu, callback_func, 0);

}

// See Changelog 2004-05-05
// 
// We are passed an GtkOptionMenu *option_menu
//
void
graphics_info_t::fill_option_menu_with_coordinates_options_internal(GtkWidget *option_menu,
								    GtkSignalFunc callback_func, 
								    short int set_last_active_flag) {

   int imol_active = -1; // To allow the function to work as it used to.
   fill_option_menu_with_coordinates_options_internal_2(option_menu, callback_func, set_last_active_flag, imol_active);

}

void
graphics_info_t::fill_option_menu_with_coordinates_options_internal_2(GtkWidget *option_menu,
								      GtkSignalFunc callback_func, 
								      short int set_last_active_flag,
								      int imol_active) { 

   // like the column labels from an mtz file, similarly fill this
   // option_menu with items that correspond to molecules that have
   // coordinates.
   //

   // Get the menu of the optionmenu (which was set in interface.c:
   // gtk_option_menu_set_menu (GTK_OPTION_MENU (go_to_atom_molecule_optionmenu),
   //                           go_to_atom_molecule_optionmenu_menu);
   //
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   

   // for the strangeness of destroying and re-adding the menu to the
   // option menu, set the comments in the
   // fill_close_option_menu_with_all_molecule_options function

   // menu is not GTK_MENU on Gtk2 Ubuntu kalypso 64 bit
   if (menu) 
      gtk_widget_destroy(menu);

   /* Create a menu for the optionmenu button.  The various molecule
    numbers will be added to this menu as menuitems*/
   menu = gtk_menu_new();

   // GtkWidget *optionmenu_menu = gtk_menu_new();
   GtkWidget *menuitem;
   // int last_imol = 0;
   int last_menu_item_index = 0;

   int menu_index = 0; // for setting of imol_active as active mol in go to atom
   for (int imol=0; imol<n_molecules(); imol++) {

//       std::cout << "in fill_option_menu_with_coordinates_options, "
// 		<< "g.molecules[" << imol << "].atom_sel.n_selected_atoms is "
// 		<< g.molecules[imol].atom_sel.n_selected_atoms << std::endl;
      
      if (molecules[imol].has_model()) { 

	 std::string ss = int_to_string(imol);
	 ss += " " ;
	 int ilen = molecules[imol].name_.length();
	 int left_size = ilen-go_to_atom_menu_label_n_chars_max;
	 if (left_size <= 0) {
	    // no chop
	    left_size = 0;
	 } else {
	    // chop
	    ss += "...";
	 } 
	 ss += molecules[imol].name_.substr(left_size, ilen);
	 menuitem = gtk_menu_item_new_with_label (ss.c_str());

	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     // GTK_SIGNAL_FUNC(go_to_atom_mol_button_select),
			     callback_func,
			     GINT_TO_POINTER(imol)); 

	 // Note that we probably don't need to do the following
	 // because we already pass a GINT_TO_POINTER(imol) in the
	 // signal connect.
	 //
	 // But on reflection.. perhaps we do because we do a
	 // menu_get_active in save_go_to_atom_mol_menu_active_position
	 // 
	 // we set user data on the menu item, so that when this goto
	 // Atom widget is cancelled, we can whatever was the molecule
	 // number corresponding to the active position of the menu
	 //
	 // Should be freed in on_go_to_atom_cancel_button_clicked
	 // (callbacks.c)
	 // 
	  
	 gtk_object_set_user_data(GTK_OBJECT(menuitem), GINT_TO_POINTER(imol));
	 gtk_menu_append(GTK_MENU(menu), menuitem); 

	 if (imol == imol_active)
	    gtk_menu_set_active(GTK_MENU(menu), menu_index);

	 // we do need this bit of course:
	 gtk_widget_show(menuitem); 
	 last_menu_item_index++;
	 menu_index++; 
      }
   }
   
   // set any previously saved active position:
   // but overridden by flag:
   if (set_last_active_flag) {
      // explanation of -1 offset: when there are 2 menu items,
      // last_menu_item_index is 2, but the item index of the last
      // item is 2 - 1.
      //
      gtk_menu_set_active(GTK_MENU(menu), (last_menu_item_index-1)); 
   } else {
      // the old way (ie. not ..._with_active_mol() mechanism)
      if (imol_active == -1) { 
	 if (go_to_atom_mol_menu_active_position >= 0) {
	    gtk_menu_set_active(GTK_MENU(menu), go_to_atom_mol_menu_active_position); 
	 }
      }
   }

   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu),
			    menu);
}


void
graphics_info_t::fill_option_menu_with_coordinates_options_internal_with_active_mol(GtkWidget *option_menu,
										    GtkSignalFunc callback_func, 
										    int imol_active) {

   short int set_last_active_flag = 0;
   fill_option_menu_with_coordinates_options_internal_2(option_menu, callback_func,
							set_last_active_flag, imol_active);
} 


// Return -2 on ambiguity, -1 on unset and a molecule number >=0 for
// no ambiguity (or undo_molecule has been set already).
//
int
graphics_info_t::Undo_molecule(coot::undo_type undo_type) const {

   int r = -1;
   if (undo_molecule > -1)
      r = undo_molecule;
   else {
      int n_mol = 0; 
      for (int imol=0; imol<n_molecules(); imol++) {
	 // Argh.  I want to store an function (name) as a variable.
	 // How do I do that?

	 if (undo_type == coot::UNDO) { 
	    if (molecules[imol].Have_modifications_p()) {
	       n_mol++;
	       r = imol;
	    }
	 }

	 if (undo_type == coot::REDO) {
	    if (molecules[imol].Have_redoable_modifications_p()) {
	       n_mol++;
	       r = imol;
	    }
	 }
      }
      if (n_mol > 1) {
	 r = -2;
      }
   }
   return r;
}

// not const
void
graphics_info_t::fill_option_menu_with_undo_options(GtkWidget *option_menu) {

   
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
   if (menu) 
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   GtkWidget *menuitem;
   int first = -1;
   
   for (int i=0; i<n_molecules(); i++) {
      // if (molecules[i].has_model()) { 
      if (molecules[i].atom_sel.mol) { 
	 if (molecules[i].Have_modifications_p()) {
	    if (first == -1)
	       first = i;
	    char s[200];
	    snprintf(s, 199, "%d", i);
	    std::string ss(s);
	    ss += " ";
	    ss += molecules[i].name_;
	    menuitem = gtk_menu_item_new_with_label(ss.c_str());
	    gtk_signal_connect(GTK_OBJECT(menuitem), "activate",
			       GTK_SIGNAL_FUNC(graphics_info_t::undo_molecule_select),
			       GINT_TO_POINTER(i));
	    gtk_menu_append(GTK_MENU(menu), menuitem);
	    gtk_widget_show(menuitem);
	 }
      }
   }

   gtk_option_menu_set_menu(GTK_OPTION_MENU(option_menu), menu);

   // first should have been set (this function should not be called
   // if there are not more than 1 molecule with modifications).
   if (first > -1) {
      set_undo_molecule_number(first);
   }
} 

// a static function
void
graphics_info_t::undo_molecule_select(GtkWidget *item, GtkPositionType pos) {
   graphics_info_t g;
   g.set_undo_molecule_number(pos);
}

void
graphics_info_t::set_baton_build_params(int istart_resno,
					const char *chain_id, 
					const char *backwards) { 

   baton_build_params_active = 1; // don't ignore baton_build_params
				  // in placing atom.
   baton_build_start_resno = istart_resno;
   std::string dir(backwards);
   if (dir == "backwards") { 
      baton_build_direction_flag = -1;
   } else { 
      if (dir == "forwards") { 
	 baton_build_direction_flag = 1;
      } else { 
	 baton_build_direction_flag = 0; // unknown.
      }
   }
   baton_build_chain_id = std::string(chain_id);

} 

void
graphics_info_t::set_bond_thickness(int imol, float t) {

   if (imol < n_molecules()) {
      if (imol >= 0) {
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    molecules[imol].set_bond_thickness(t);
	    graphics_draw();
	 }
      }
   } else {
      std::cout << "Ignoring attempt to set bond with for molecule "
		<< imol << std::endl;
   }
} 

void
graphics_info_t::set_bond_thickness_intermediate_atoms(float t) {

   bond_thickness_intermediate_atoms = t;

} 

void
graphics_info_t::crosshairs_text() const { 

   if (draw_crosshairs_flag > 0) { 
      std::cout << "Crosshair ticks: 1.54A (C-C bond), 2.7A (H-bond), 3.8A (Ca-Ca)\n";
   }
} 

// Thank you for this idea Stuart Makay
// 
void
graphics_info_t::pick_cursor_maybe() {

   if (control_key_for_rotate_flag) {
      pick_cursor_real();
   }
}

void
graphics_info_t::pick_cursor_real() {

   if (use_graphics_interface_flag) { 
      //    GdkCursorType c = GDK_CROSSHAIR;
      GdkCursorType c = pick_cursor_index;
      GdkCursor *cursor;
      cursor = gdk_cursor_new (c);
      gdk_window_set_cursor (glarea->window, cursor);
      gdk_cursor_destroy (cursor);
   }
}

// static
void
graphics_info_t::normal_cursor() { 

   if (use_graphics_interface_flag) { 
      if (control_key_for_rotate_flag) { 
	 GdkCursorType c = GDK_LEFT_PTR;
	 GdkCursor *cursor;
	 cursor = gdk_cursor_new (c);
	 gdk_window_set_cursor (glarea->window, cursor);
	 gdk_cursor_destroy (cursor);
      }
   }
}

// static
void
graphics_info_t::watch_cursor() { 

   if (use_graphics_interface_flag) { 
      GdkCursorType c = GDK_WATCH;
      GdkCursor *cursor;
      cursor = gdk_cursor_new (c);
      gdk_window_set_cursor (glarea->window, cursor);
      gdk_cursor_destroy (cursor);
      while (gtk_events_pending()) {
	 gtk_main_iteration();
      }
   }
}

// static
void
graphics_info_t::fleur_cursor() { 

   if (use_graphics_interface_flag) { 
      GdkCursorType c = GDK_FLEUR;
      GdkCursor *cursor;
      cursor = gdk_cursor_new (c);
      gdk_window_set_cursor (glarea->window, cursor);
      gdk_cursor_destroy (cursor);
   }
}



// static
short int 
graphics_info_t::alt_conf_split_type_number() { 
   
   return graphics_info_t::alt_conf_split_type;

} 

void
graphics_info_t::model_fit_refine_unactive_togglebutton(const std::string &button_name) const { 

   if (model_fit_refine_dialog) {
      GtkWidget *toggle_button = lookup_widget(model_fit_refine_dialog, button_name.c_str());
      if (toggle_button)
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
      else 
	 std::cout << "ERROR:: failed to find button: " << button_name << std::endl;
      
   } else {
      // std::cout << "DEBUG:: model_fit_refine_dialog not found" << std::endl;
   }

#if (GTK_MAJOR_VERSION > 1)

   std::string toolbar_button_name = "not-found";
   if (button_name == "model_refine_dialog_refine_togglebutton")
      toolbar_button_name = "model_toolbar_refine_togglebutton";
   if (button_name == "model_refine_dialog_regularize_zone_togglebutton")
      toolbar_button_name = "model_toolbar_regularize_togglebutton";
   if (button_name == "model_refine_dialog_rigid_body_togglebutton")
      toolbar_button_name = "model_toolbar_rigid_body_fit_togglebutton";
   if (button_name == "model_refine_dialog_rot_trans_togglebutton")
      toolbar_button_name = "model_toolbar_rot_trans_togglebutton";
   if (button_name == "model_refine_dialog_auto_fit_rotamer_togglebutton")
      toolbar_button_name = "model_toolbar_auto_fit_rotamer_togglebutton";
   if (button_name == "model_refine_dialog_rotamer_togglebutton")
      toolbar_button_name = "model_toolbar_rotamers_togglebutton";
   if (button_name == "model_refine_dialog_edit_chi_angles_togglebutton")
      toolbar_button_name = "model_toolbar_edit_chi_angles_togglebutton";
   if (button_name == "model_refine_dialog_pepflip_togglebutton")
      toolbar_button_name = "model_toolbar_flip_peptide_togglebutton";
   if (button_name == "model_refine_dialog_do_180_degree_sidechain_flip_togglebutton")
      toolbar_button_name = "model_toolbar_sidechain_180_togglebutton";
   if (button_name == "model_refine_dialog_mutate_auto_fit_togglebutton")
      toolbar_button_name = "model_toolbar_mutate_and_autofit_togglebutton";
   if (button_name == "model_refine_dialog_mutate_togglebutton")
      toolbar_button_name = "model_toolbar_simple_mutate_togglebutton";
   if (button_name == "model_refine_dialog_fit_terminal_residue_togglebutton")
      toolbar_button_name = "model_toolbar_add_terminal_residue_togglebutton";

   // now, button_name may have been
   // model_refine_dialog_edit_phi_psi_togglebutton or
   // model_refine_dialog_edit_backbone_torsions_togglebutton, we
   // don't have toolbar equivalents of those.
   // 
   if (toolbar_button_name != "not-found") { 
      GtkWidget *toggle_button = lookup_widget(graphics_info_t::glarea,
					       toolbar_button_name.c_str());
//       std::cout << "DEBUG:: toggle_button for gtk2 toolbar: " << button_name << "->"
// 		<< toolbar_button_name << " " << toggle_button << std::endl;

      if (toggle_button) {
	// somehow we cannot use ->active on the toggle_tool_buttons?!
	gboolean active = gtk_toggle_tool_button_get_active(GTK_TOGGLE_TOOL_BUTTON(toggle_button));
	if (active)
	  gtk_toggle_tool_button_set_active(GTK_TOGGLE_TOOL_BUTTON(toggle_button), FALSE);
      }
   }
#endif    
   
} 


void
graphics_info_t::other_modelling_tools_unactive_togglebutton(const std::string &button_name) const { 

   if (other_modelling_tools_dialog) {
      GtkWidget *toggle_button = lookup_widget(other_modelling_tools_dialog,
					       button_name.c_str());
      if (toggle_button)
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
      else 
	 std::cout << "ERROR:: failed to find button: " << button_name
		   << std::endl;
   }
}


// The start of the edit one residue with the pullable point in the 
// Ramachandran.  Currently disabled, but will be a nice feature in future.
// 
void
graphics_info_t::execute_edit_phi_psi(int atom_index, int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS) 
   std::pair<double, double> phi_psi = molecules[imol].get_phi_psi(atom_index);

   if (phi_psi.first > -200.0) { 
      coot::rama_plot *plot = new coot::rama_plot;

      plot->init("phi/psi-edit");  // magic string
      coot::util::phi_psi_t phipsi(phi_psi.first, phi_psi.second, "resname", 
				   "moving residue", 1, "inscode", "chainid");
      plot->draw_it(phipsi);
      
      moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
      imol_moving_atoms = imol;
      short int whole_res_flag = 1;
      atom_selection_container_t residue_asc = 
	 graphics_info_t::molecules[imol].edit_residue_pull_residue(atom_index,
								    whole_res_flag);
      make_moving_atoms_graphics_object(residue_asc);
      
      graphics_draw();

   } else { 
      std::cout << "Can't find ramachandran angles for this residue" << std::endl;
   } 
#endif // HAVE_GTK_CANVAS   
}

void
graphics_info_t::rama_plot_for_single_phi_psi(int imol, int atom_index) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   std::pair<double, double> phi_psi = molecules[imol].get_phi_psi(atom_index);

   if (phi_psi.first > -200.0) { 
      // coot::rama_plot *plot = new coot::rama_plot;
      edit_phi_psi_plot = new coot::rama_plot;

      edit_phi_psi_plot->init("backbone-edit");  // magic string

      // construct label:
      std::string label;
      CAtom *at = molecules[imol].atom_sel.atom_selection[atom_index];
      int resno            = at->GetSeqNum();
      std::string chain_id = at->GetChainID();
      label = int_to_string(resno);
      label += chain_id;

      coot::util::phi_psi_t phipsi(phi_psi.first, phi_psi.second, "resname", 
				   label, 1, "inscode", "chainid");
      edit_phi_psi_plot->draw_it(phipsi);
      
   }
#endif // HAVE_GTK_CANVAS
}

void
graphics_info_t::rama_plot_for_2_phi_psis(int imol, int atom_index) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   std::pair<double, double> phi_psi = molecules[imol].get_phi_psi(atom_index);

   if (phi_psi.first > -200.0) { 

      edit_phi_psi_plot = new coot::rama_plot;

      edit_phi_psi_plot->init("backbone-edit");  // magic string

      // construct label:
      std::string label;
      CAtom *at = molecules[imol].atom_sel.atom_selection[atom_index];
      int resno            = at->GetSeqNum();
      std::string chain_id = at->GetChainID();
      label = int_to_string(resno);
      label += chain_id;

      coot::util::phi_psi_t phipsi(phi_psi.first, phi_psi.second, "resname", 
				   label, 1, "inscode", "chainid");
      edit_phi_psi_plot->draw_it(phipsi);
      
   }
#endif // HAVE_GTK_CANVAS
}

// activated from the edit torsion angles cancel button (and OK button, I
// suppose).
// 
void 
graphics_info_t::destroy_edit_backbone_rama_plot() {  // only one of these.

   printf("start graphics_info_t::destroy_edit_backbone_rama_plot()\n");

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   if (edit_phi_psi_plot) { 
      // we need to get to the widget "dynarama_window" and destroy it.
      edit_phi_psi_plot->destroy_yourself();
      edit_phi_psi_plot = 0; // Richard Baxter bug
   } else { 
      std::cout << "WARNING:: edit_phi_psi_plot is NULL\n";
   } 
#endif // HAVE_GTK_CANVAS
   printf("done in graphics_info_t::destroy_edit_backbone_rama_plot()\n");
}


// This is called as part of the callback of moving a point in the ramachandran
// edit window
// 
// We (Kevin and I) currently don't like the way the atoms move, so the button
// is invisible now.  However, it will be reinstated in the future, when we
// change residue_edit_phi_psi to use torsion restraints.  residue_edit_phi_psi
// should return a double which is the distortion at the end of the refinement.
// The box should be coloured according to distortion.  We should add a help
// button underneath the canvas (perhaps not visible unless made so, when we are
// in a edit-one-residue situation) [25Feb2004].
// 
void
graphics_info_t::set_edit_phi_psi_to(double phi, double psi) { 

   // tinker with the coordinates of the moving_atoms_asc

   short int istat = molecules[imol_moving_atoms].residue_edit_phi_psi(*moving_atoms_asc, edit_phi_psi_atom_index, phi, psi);

   if (istat) {
      // You would do this much less heay-weightedly if you have your
      // bonding already sorted out here.
      // 
      int do_disulphide_flag = 0;
      Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
      graphics_draw();
   }
}


// There are 2 places that have to be united in how they get their torsions.
//
// Here is 1) - (or more precisely in
// wrapped_create_edit_chi_angles_dialog) i.e. where we make create
// the buttons for the torsions.
// 
// These torsion labels/indices must match:
// 2) update_residue_by_chi_change()/ chi_angles::change_by()
// 
// Particularly we should exclude hydrogen rotations consistently.
// Note that setup_flash_bond_internal() needs to be consistent too.
// 
void
graphics_info_t::execute_edit_chi_angles(int atom_index, int imol) {

   // check that we have chis for this residue:
   int n_chis = molecules[imol].N_chis(atom_index);

   if (n_chis) {

      std::string res_type(molecules[imol].atom_sel.atom_selection[atom_index]->residue->GetResName());
      chi_angles_clicked_atom_spec =
	 coot::atom_spec_t(molecules[imol].atom_sel.atom_selection[atom_index]);
      chi_angles_clicked_atom_spec.int_user_data = 1; // not magic "don't use" value

      // Make Phil Evans happy (well, slightly happier.. :-)
      if (res_type == "MSE")
	 chi_angles_clicked_atom_spec.atom_name = " C  ";
      if (res_type == "ARG")
	 chi_angles_clicked_atom_spec.atom_name = " C  ";
      if (res_type == "PHE")
	 chi_angles_clicked_atom_spec.atom_name = " C  ";
      if (res_type == "TYR")
	 chi_angles_clicked_atom_spec.atom_name = " C  ";

      // belt and braces:
      if ( (res_type == "GLY") || (res_type == "ALA") ) {
	 std::cout << "This residue does not have chi angles." << std::endl;
      } else { 

	 // copy the residue, just like we do in execute_edit_phi_psi:
	 // 
	 moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
	 imol_moving_atoms = imol;
	 short int whole_res_flag = 0; // We only want to pull the
				       // atoms of *this* alternative
				       // conformation (and this
				       // includes atoms with altconf
				       // "").
	 atom_selection_container_t residue_asc = 
	    graphics_info_t::molecules[imol].edit_residue_pull_residue(atom_index,
								       whole_res_flag);

	 regularize_object_bonds_box.clear_up();
   
	 int ires = wrapped_create_edit_chi_angles_dialog(res_type);
	 if (ires > 0) { 
	    std::cout << "Use the 1,2,3,4 keys to select rotamers, 0 for "
		      << "normal rotation mode" << std::endl;
	    make_moving_atoms_graphics_object(residue_asc);

	    if (do_probe_dots_on_rotamers_and_chis_flag) {
	       setup_for_probe_dots_on_chis_molprobity(imol);
	    }
	 } else {
	    std::cout << "WARNING:: couldn't find torsions in the dictionary "
		      << "for this residue: " << res_type << std::endl;
	 }
	 graphics_draw();
      }
   } else {
      std::cout << "This residue does not have chi angles." << std::endl;
      std::cout << "Missing dictionary, perhaps? " << std::endl;
   }
}

void
graphics_info_t::setup_for_probe_dots_on_chis_molprobity(int imol) {

   if (moving_atoms_asc->n_selected_atoms) {

      // make a directory where the probe dots will go:
      int dir_status = coot::util::create_directory("coot-molprobity");
      
      int n_atoms = moving_atoms_asc->n_selected_atoms;
      molecules[imol].atom_sel.mol->WritePDBASCII("molprobity-tmp-reference-file.pdb");
      
      // so find the centre and radius of the set of moving atoms:
      coot::Cartesian acc(0,0,0);
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 coot::Cartesian pt(moving_atoms_asc->atom_selection[i]->x,
			    moving_atoms_asc->atom_selection[i]->y,
			    moving_atoms_asc->atom_selection[i]->z);
	 acc += pt;
      }
      coot::Cartesian av(acc.x()/float(n_atoms),
			 acc.y()/float(n_atoms),
			 acc.z()/float(n_atoms));
      probe_dots_on_chis_molprobity_centre = av;
      float max_d = 0;
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 coot::Cartesian pt(moving_atoms_asc->atom_selection[i]->x,
			    moving_atoms_asc->atom_selection[i]->y,
			    moving_atoms_asc->atom_selection[i]->z);
	 float this_d = (pt - av).amplitude();
	 if (this_d > max_d)
	    max_d = this_d;
      }

      // so we have the maximum radius of the atoms from the centre
      // point, but we should enlargen the because we are swinging chi
      // and we need to make contact to the reference protein atoms:
      probe_dots_on_chis_molprobity_radius = (max_d + 2) * 1.7;
      if (dir_status == 0) { // success
	 do_probe_dots_on_rotamers_and_chis();
      }
   }
}

int
graphics_info_t::wrapped_create_edit_chi_angles_dialog(const std::string &res_type) {

   GtkWidget *dialog = create_edit_chi_angles_dialog();
   set_transient_and_position(COOT_EDIT_CHI_DIALOG, dialog);
   
   // Fill the vbox with buttons with atom labels about which there
   // are rotatable torsions:
   //
   GtkWidget *vbox = lookup_widget(dialog,"edit_chi_angles_vbox");

   std::vector <coot::dict_torsion_restraint_t> v =
      get_monomer_torsions_from_geometry(res_type);

   // We introduce here ichi (which gets incremented if the current
   // torsion is not const), we do that so that we have consistent
   // indexing in the torsions vector with chi_angles's change_by()
   // (see comments above execute_edit_chi_angles()).

   int ichi = 0;
   for (unsigned int i=0; i<v.size(); i++) {
      if (!v[i].is_const()) {
	 std::string label = "  ";
	 label += v[i].id(); 
	 label += "  ";
	 label += v[i].atom_id_2_4c();
	 label += " <--> ";
	 label += v[i].atom_id_3_4c();
	 label += "  ";
	 GtkWidget *button = gtk_button_new_with_label(label.c_str());
	 gtk_signal_connect(GTK_OBJECT(button), "clicked",
			    GTK_SIGNAL_FUNC(on_change_current_chi_button_clicked),
			    GINT_TO_POINTER(ichi));
	 gtk_signal_connect(GTK_OBJECT(button), "enter",
			    GTK_SIGNAL_FUNC(on_change_current_chi_button_entered),
			    GINT_TO_POINTER(ichi));
	 gtk_widget_set_name(button, "edit_chi_angles_button");
	 gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
	 gtk_container_set_border_width(GTK_CONTAINER(button), 2);
	 gtk_widget_show(button);
	 ichi++;
      }
   }
   gtk_widget_show(dialog);

   return v.size();
}

// static
void
graphics_info_t::on_change_current_chi_button_clicked(GtkButton *button,
						      gpointer user_data) {

   graphics_info_t g;
   int i = GPOINTER_TO_INT(user_data);
   g.edit_chi_current_chi = i + 1;
   g.in_edit_chi_mode_flag = 1;

} 

// static
void
graphics_info_t::on_change_current_chi_button_entered(GtkButton *button,
						      gpointer user_data) {

   graphics_info_t g;
   int ibond_user = GPOINTER_TO_INT(user_data);

   g.setup_flash_bond_internal(ibond_user+1);

}


void
graphics_info_t::setup_flash_bond_internal(int ibond_user) {

   int bond = ibond_user + 1;
   // turn it off first, only enable it if we find a pair:
   draw_chi_angle_flash_bond_flag = 0; // member data item

   std::cout << "flash bond ibond_user: " << ibond_user << std::endl;

   // std::cout << "highlight ligand bond " << bond << std::endl;

   // get the residue type and from that the atom name pairs:
   // 
   if (! moving_atoms_asc) {
      std::cout << "ERROR: moving_atoms_asc is NULL" << std::endl;
   } else { 
      if (moving_atoms_asc->n_selected_atoms == 0) {
	 std::cout << "ERROR: no atoms in moving_atoms_asc" << std::endl;
      } else { 
	 CModel *model_p = moving_atoms_asc->mol->GetModel(1);
	 if (model_p) {
	    CChain *chain_p = model_p->GetChain(0);
	    if (chain_p) {
	       CResidue *residue_p = chain_p->GetResidue(0);
	       if (residue_p) {

		  std::string residue_type(residue_p->GetResName());
		  std::vector<std::vector<int> > contact_indices =
		     coot::util::get_contact_indices_from_restraints(residue_p, geom_p, 0);
		  coot::chi_angles c(residue_p, 0);

		  std::pair<std::string, std::string> atom_names;

		  // c.f. get_torsion_bonds_atom_pairs in chi-angles.cc 
		  std::vector <coot::dict_torsion_restraint_t> monomer_torsions = 
		     geom_p->get_monomer_torsions_from_geometry(residue_type);

		  std::pair<short int, coot::dictionary_residue_restraints_t> r =
		     geom_p->get_monomer_restraints(residue_type);

// 		  std::cout << "DEBUG:: monomer_torsions.size(): "
// 			    << monomer_torsions.size() << std::endl;

		  int hydrogen_or_const_torsion_count = 0;
// debugging		  
// 		  for(unsigned int i=0; i<monomer_torsions.size(); i++) {
// 		     std::string this_one_str;
// 		     std::string atom1 = monomer_torsions[i].atom_id_1();
// 		     std::string atom2 = monomer_torsions[i].atom_id_4();
// 		     if ( (!r.second.is_hydrogen(atom1) && !r.second.is_hydrogen(atom2))
// 			  || find_hydrogen_torsions) {
// 		     } else {
// 			hydrogen_torsion_count++;
// 		     }
// 		     if ((bond - 2 + hydrogen_torsion_count) == i)
// 			this_one_str = " <==";
// 		     std::cout << "   DEBUG:: torsion i:" << i << " htc:"
// 			       << hydrogen_torsion_count << " match:"
// 			       << (bond - 2 + hydrogen_torsion_count) << " "
// 			       << monomer_torsions[i].atom_id_1() << " "
// 			       << monomer_torsions[i].atom_id_2() << " "
// 			       << monomer_torsions[i].atom_id_3() << " "
// 			       << monomer_torsions[i].atom_id_4() << " "
// 			       << monomer_torsions[i].periodicity() << " "
// 			       << monomer_torsions[i].angle() << " "
// 			       << monomer_torsions[i].esd() << this_one_str
// 			       << std::endl;
		     
// 		  }

		  std::cout << "    got " << monomer_torsions.size() << " monomer torsions "
			    << std::endl;

		  hydrogen_or_const_torsion_count = 0;
		  if (monomer_torsions.size() > 0) { 
		     for(unsigned int i=0; i<monomer_torsions.size(); i++) {

			if (!monomer_torsions[i].is_const()) { 
			   std::string atom1 = monomer_torsions[i].atom_id_1();
			   std::string atom2 = monomer_torsions[i].atom_id_4();
			   int hydrogen_torsion_val = 0;
			   if (r.second.is_hydrogen(atom1)) hydrogen_torsion_val = 1;
			   if (r.second.is_hydrogen(atom2)) hydrogen_torsion_val = 1;
			   if (!hydrogen_torsion_val || find_hydrogen_torsions) {
			      std::cout << "   comparing (" << bond << ") "
					<< bond - 2 + hydrogen_or_const_torsion_count
					<< " vs " << i << std::endl;
			      if ((bond - 2 + hydrogen_or_const_torsion_count) == int(i)) {
				 std::string atom2 = monomer_torsions[i].atom_id_2();
				 std::string atom3 = monomer_torsions[i].atom_id_3();
				 atom_names = std::pair<std::string, std::string> (atom2, atom3);
				 std::cout << "   match atom names :" << atom_names.first
					   << ": :" << atom_names.second << ":" << std::endl;
				 break; 
			      }
			   }
			   if (hydrogen_torsion_val)
			      hydrogen_or_const_torsion_count++;
			} else {
			   hydrogen_or_const_torsion_count++;
			}
		     }
		  }

		  std::cout << "         :" << atom_names.first << ": :"
			    << atom_names.second << ":" << std::endl;

		  if ((atom_names.first != "") &&
		      (atom_names.second != "") && 
		      (atom_names.first != "empty") &&
		      (atom_names.second != "empty")) {
		     
		     PPCAtom residue_atoms;
		     int nResidueAtoms;
		     residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
		     
		     if (nResidueAtoms > 0) { // of course it is!
			for (int iat1=0; iat1<nResidueAtoms; iat1++) {
			   std::string ra1=residue_atoms[iat1]->name;
			   if (ra1 == atom_names.first) {
			      for (int iat2=0; iat2<nResidueAtoms; iat2++) {
				 std::string ra2=residue_atoms[iat2]->name;
				 if (ra2 == atom_names.second) {
				    
				    draw_chi_angle_flash_bond_flag = 1;
				    clipper::Coord_orth p1(residue_atoms[iat1]->x,
							   residue_atoms[iat1]->y,
							   residue_atoms[iat1]->z);
				    clipper::Coord_orth p2(residue_atoms[iat2]->x,
							   residue_atoms[iat2]->y,
							   residue_atoms[iat2]->z);

				    std::pair<clipper::Coord_orth, clipper::Coord_orth> cp(p1, p2);

				    graphics_info_t g;
				    g.add_flash_bond(cp);
				    graphics_draw();
				 }
			      }
			   }
			}
		     } 
		  } 
	       }
	    }
	 }
      }
   }
}

void
graphics_info_t::add_flash_bond(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &bond) {

   draw_chi_angle_flash_bond_flag = 1;
   flash_bond = bond;

} 

// static
void graphics_info_t::draw_chi_angles_flash_bond() {

   if (draw_chi_angle_flash_bond_flag) {
      glLineWidth(10);
      glColor3f(0.3,1.0,0.3);
      glBegin(GL_LINES);
      glVertex3f(graphics_info_t::flash_bond.first.x(),
		 graphics_info_t::flash_bond.first.y(),
		 graphics_info_t::flash_bond.first.z());
      glVertex3f(graphics_info_t::flash_bond.second.x(),
		 graphics_info_t::flash_bond.second.y(),
		 graphics_info_t::flash_bond.second.z());
      glEnd();
   }
} 


// ----------------------------------------------------------------------------
//                          map colour stuff
// ----------------------------------------------------------------------------

void
graphics_info_t::set_last_map_colour(double f1, double f2, double f3) const { 

   int imap = -1; 
   for (int i=0; i<n_molecules(); i++) { 
      if (molecules[i].has_map()) { 
	 imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of colour" << std::endl;
   } else {
      double *colours = new double[4];
      colours[0] = f1;
      colours[1] = f2;
      colours[2] = f3;
      molecules[imap].handle_map_colour_change(colours, swap_difference_map_colours,
					       GL_CONTEXT_MAIN);
      if (display_mode_use_secondary_p()) {
	 make_gl_context_current(GL_CONTEXT_SECONDARY);
	 molecules[imap].handle_map_colour_change(colours, swap_difference_map_colours,
					       GL_CONTEXT_SECONDARY);
	 make_gl_context_current(GL_CONTEXT_MAIN);
      } 
      delete [] colours;
   }
}

void
graphics_info_t::set_last_map_contour_level(float level) {

   int imap = -1; 
   for (int i=0; i<n_molecules(); i++) { 
      if (molecules[i].has_map()) { 
	 imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of contour" << std::endl;
   } else {
      molecules[imap].set_contour_level(level);
   }
}

void
graphics_info_t::set_last_map_contour_level_by_sigma(float f) {

   int imap = -1; 
   for (int i=0; i<n_molecules(); i++) { 
      if (molecules[i].has_map()) { 
	 imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of contour" << std::endl;
   } else {
      molecules[imap].set_contour_level_by_sigma(f);
   }
} 


// And turn it on.
void 
graphics_info_t::set_last_map_sigma_step(float f) { 

   int imap = -1;
   for (int i=0; i<n_molecules(); i++) { 
      if (molecules[i].has_map()) { 
	 imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of contour step" 
		<< std::endl;
   } else {
//       molecules[imap].contour_by_sigma_flag = 1;
//       molecules[imap].contour_sigma_step = f;
      molecules[imap].set_contour_by_sigma_step(f, 1);
   }
}


// Create a moving atoms molecule, consisting of the Ca(n), Ca(n+1) of
// the peptide, N(n) C(n+1), O(n+1).  Note the alt conf should be the
// same (we can have altconfed mainchain).
// 
void
graphics_info_t::execute_setup_backbone_torsion_edit(int imol, int atom_index) { 
   
   // if not the same altconf for all atoms, give up (for now).
   // 
   // Do an atom selection of this residue and the previous one
   // (unless this was a N, then do this one and the next)
   // 
   // Run through the atom selection finding the atoms.  put the atoms
   // into new residues (appropriate) and construct a mol and and asc
   // and make that the moving atoms asc.  Be sure to put the residue
   // with the peptide N first so that mmdb finds it first when it
   // does the atom selection that's part of make_asc(mol) - which
   // means that get_first_atom_with_atom_name() will work like we
   // want it to.

   if (imol < n_molecules()) { 
      if (molecules[imol].has_model()) { 
	 if (atom_index < molecules[imol].atom_sel.n_selected_atoms) { 
	    CAtom *this_atom_p = molecules[imol].atom_sel.atom_selection[atom_index];
	    int offset = 0; // usually this and next residue
	    std::string this_atname(this_atom_p->name);
	    if (this_atname == " N  ") { 
	       offset = -1;
	    }
	    int this_res = this_atom_p->GetSeqNum();
	    char *chain_id = this_atom_p->GetChainID();
	    int SelectionHandle = molecules[imol].atom_sel.mol->NewSelection();
	    char *ins_code = this_atom_p->GetInsCode();
	    char *altconf  = this_atom_p->altLoc;
	    std::string a_tmp(altconf);
	    if (a_tmp != "") { 
	       a_tmp += ",";
	    }
	    char *search_altconf = (char *) a_tmp.c_str();
	    molecules[imol].atom_sel.mol->SelectAtoms (SelectionHandle, 0, chain_id,
					      this_res + offset , // starting resno
					      ins_code, // any insertion code
					      this_res + 1 + offset, // ending resno
					      ins_code, // ending insertion code
					      "*", // any residue name
					      "*", // atom name
					      "*", // elements
					      search_altconf  // alt loc.
					      );
	    int nSelAtoms;
	    PPCAtom SelectAtoms;
	    molecules[imol].atom_sel.mol->GetSelIndex(SelectionHandle, SelectAtoms, nSelAtoms);
	    if (nSelAtoms < 4) { // the very min
	       std::cout <<  "WARNING:: not enough atoms in atom selection in "
			 <<  "execute_setup_backbone_torsion_edit" << std::endl;
	    } else {


	       // We construct a moving atom asc atom by atom...
	       CAtom *next_ca = NULL, *next_n = NULL;
	       CAtom *this_c = NULL, *this_o = NULL, *this_ca = NULL;
	       
	       // and the extra atoms that we need (we extract the
	       // coordinates) to construct a pair of ramachandran
	       // points in the rama_plot:

	       CAtom *prev_c = NULL;
	       CAtom *this_n = NULL;
	       CAtom *next_c = NULL;
	       CAtom *next_plus_1_n = NULL;

	       // to get prev_c and next_plus_1_n we do atom selections:

	       // You can add (this) hydrogen to that list if you want.
	       CAtom *at;
	       // 
	       for (int iat=0; iat<nSelAtoms; iat++) { 
		  at = SelectAtoms[iat];
		  if (at->GetSeqNum() == (this_res + offset) ) { 
		     std::string n(at->name);
		     if (n == " CA ") { 
			this_ca = at;
		     }
		     if (n == " O  ") { 
			this_o = at;
		     }
		     if (n == " C  ") { 
			this_c = at;
		     }
		     if (n == " N  ") { 
			this_n = at;
		     }
		  }
		  if (at->GetSeqNum() == (this_res + 1 + offset) ) { 
		     std::string n(at->name);
		     if (n == " CA ") { 
			next_ca = at;
		     }
		     if (n == " N  ") { 
			next_n = at;
		     }
		     if (n == " C  ") { 
			next_c = at;
		     }
		  }
	       }

	       int SelHnd_prev_c = molecules[imol].atom_sel.mol->NewSelection();
	       // to get prev_c and next_plus_1_n we do atom selections:
	       molecules[imol].atom_sel.mol->SelectAtoms (SelHnd_prev_c, 0, chain_id,
							  this_res - 1 + offset , // starting resno
							  ins_code, // any insertion code
							  this_res - 1 + offset, // ending resno
							  ins_code, // ending insertion code
							  "*",    // any residue name
							  " C  ", // atom name
							  "*",    // elements
							  search_altconf  // alt loc.
							  );

	       int nSelAtoms_prev_c;
	       PPCAtom SelectAtoms_prev_c;

	       molecules[imol].atom_sel.mol->GetSelIndex(SelHnd_prev_c, 
							 SelectAtoms_prev_c, 
							 nSelAtoms_prev_c);
	       if (nSelAtoms_prev_c > 0) { 
		  prev_c = SelectAtoms_prev_c[0];
	       } else { 
		  std::cout << "Oops:: didn't find prev_c\n";
	       } 
	       molecules[imol].atom_sel.mol->DeleteSelection(SelHnd_prev_c);
	       
	       // And similarly for next_plus_1_n:
	       int SelHnd_next_plus_1_n = molecules[imol].atom_sel.mol->NewSelection();
	       molecules[imol].atom_sel.mol->SelectAtoms (SelHnd_next_plus_1_n, 0, chain_id,
							  this_res + 2 + offset , // starting resno
							  ins_code, // any insertion code
							  this_res + 2 + offset, // ending resno
							  ins_code, // ending insertion code
							  "*",    // any residue name
							  " N  ", // atom name
							  "*",    // elements
							  search_altconf  // alt loc.
							  );
	       int nSelAtoms_next_plus_1_n;
	       PPCAtom SelectAtoms_next_plus_1_n;

	       molecules[imol].atom_sel.mol->GetSelIndex(SelHnd_next_plus_1_n, 
							 SelectAtoms_next_plus_1_n, 
							 nSelAtoms_next_plus_1_n);
	       if (nSelAtoms_next_plus_1_n > 0) { 
		  next_plus_1_n = SelectAtoms_next_plus_1_n[0];
	       } else {
		  std::cout << "Oops:: didn't find next + 1 N\n";
	       } 
	       molecules[imol].atom_sel.mol->DeleteSelection(SelHnd_next_plus_1_n);

	       
	       if (next_ca && next_n && this_ca && this_o && this_c) { 

		  // new addition 25Feb2004
		  rama_plot_for_2_phi_psis(imol, atom_index);

		  CMMDBManager *mol = new CMMDBManager;
		  CModel *model = new CModel;
		  CChain *chain = new CChain;
		  CResidue *res1 = new CResidue;
		  CResidue *res2 = new CResidue;
		  res1->SetResName(this_ca->GetResName());
		  res2->SetResName(next_ca->GetResName());
		  res1->seqNum = this_ca->GetSeqNum();
		  res2->seqNum = next_ca->GetSeqNum();
		  chain->SetChainID(this_ca->GetChainID());

		  at = new CAtom; 
		  at->Copy(this_ca);
		  res1->AddAtom(at);

		  at = new CAtom; 
		  at->Copy(this_c);
		  res1->AddAtom(at);

		  at = new CAtom; 
		  at->Copy(this_o);
		  res1->AddAtom(at);

		  at = new CAtom; 
		  at->Copy(next_ca);
		  res2->AddAtom(at);

		  at = new CAtom; 
		  at->Copy(next_n);
		  res2->AddAtom(at);

		  chain->AddResidue(res1);
		  chain->AddResidue(res2);
		  model->AddChain(chain);
		  mol->AddModel(model);
		  mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
		  mol->FinishStructEdit();
		  imol_moving_atoms = imol;
		  moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
		  atom_selection_container_t asc = make_asc(mol);
		  regularize_object_bonds_box.clear_up();
		  make_moving_atoms_graphics_object(asc);

		  // save the fixed end points:
		  backbone_torsion_end_ca_1 = 
		     clipper::Coord_orth(this_ca->x, this_ca->y, this_ca->z);
		  backbone_torsion_end_ca_2 = 
		     clipper::Coord_orth(next_ca->x, next_ca->y, next_ca->z);
		  

		  // add to rama_points:
		  rama_points.clear();
		  rama_points.add("this_ca", backbone_torsion_end_ca_1);
		  rama_points.add("next_ca", backbone_torsion_end_ca_2);
		  // the moving atoms:
		  rama_points.add("this_c", clipper::Coord_orth(this_c->x,
								this_c->y,
								this_c->z));
		  rama_points.add("this_o", clipper::Coord_orth(this_o->x,
								this_o->y,
								this_o->z));
		  rama_points.add("next_n", clipper::Coord_orth(next_n->x,
								next_n->y,
								next_n->z));

		  if (this_n) 
		     rama_points.add("this_n",  clipper::Coord_orth(this_n->x,
								    this_n->y,
								    this_n->z));

		  if (next_c) 
		     rama_points.add("next_c",  clipper::Coord_orth(next_c->x,
								    next_c->y,
								    next_c->z));

		  // next_plus_1_n, prev_c;
		  if (next_plus_1_n)
		     rama_points.add("next+1_n", clipper::Coord_orth(next_plus_1_n->x,
								     next_plus_1_n->y,
								     next_plus_1_n->z));
		  if (prev_c)
		     rama_points.add("prev_c", clipper::Coord_orth(prev_c->x,
								   prev_c->y,
								   prev_c->z));

// 		  std::cout << "DEBUG:: backbone_torsion_end_ca_1: " 
// 			    << backbone_torsion_end_ca_1.format() << std::endl;
// 		  std::cout << "DEBUG:: backbone_torsion_end_ca_2: " 
// 			    << backbone_torsion_end_ca_2.format() << std::endl;

		  graphics_draw();
		  GtkWidget *widget = create_edit_backbone_torsions_dialog();
		  set_edit_backbone_adjustments(widget);
		  gtk_widget_show(widget);
		  
	       } else { 
		  std::cout << "WARNING:: not all atoms found in " 
			    << "execute_setup_backbone_torsion_edit" << std::endl;
		     
	       } 
	    }
	    // 
	    molecules[imol].atom_sel.mol->DeleteSelection(SelectionHandle);
	 } 
      }
   }
}

void
graphics_info_t::set_edit_backbone_adjustments(GtkWidget *widget) { 

   GtkWidget *hscale_peptide = lookup_widget(widget, 
					     "edit_backbone_torsions_rotate_peptide_hscale");

   GtkWidget *hscale_carbonyl = lookup_widget(widget, 
					     "edit_backbone_torsions_rotate_carbonyl_hscale");

//    gfloat value,
//    gfloat lower,
//    gfloat upper,
//    gfloat step_increment,
//    gfloat page_increment,
//    gfloat page_size

   GtkAdjustment *adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -180.0, 360.0, 0.1, 1.0, 180));
   gtk_range_set_adjustment(GTK_RANGE(hscale_peptide), adjustment);
   gtk_signal_connect(GTK_OBJECT(adjustment), "value_changed",
		      GTK_SIGNAL_FUNC(graphics_info_t::edit_backbone_peptide_changed_func), NULL);

   
   // and the carbonyl:
   adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -180.0, 360.0, 0.1, 1.0, 180));
   gtk_range_set_adjustment(GTK_RANGE(hscale_carbonyl), adjustment);
   gtk_signal_connect(GTK_OBJECT(adjustment), "value_changed",
		      GTK_SIGNAL_FUNC(graphics_info_t::edit_backbone_carbonyl_changed_func), NULL);

}

// static 
void
graphics_info_t::edit_backbone_peptide_changed_func(GtkAdjustment *adj, GtkWidget *window) { 
   
   graphics_info_t g;
   // std::cout << "change backbone peptide by: " << adj->value << std::endl;
   
   std::pair<short int, clipper::Coord_orth> this_c = rama_points.get("this_c");
   std::pair<short int, clipper::Coord_orth> this_o = rama_points.get("this_o");
   std::pair<short int, clipper::Coord_orth> next_n = rama_points.get("next_n");

   if (this_c.first && this_o.first && next_n.first) { 
      CAtom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);
      CAtom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
      CAtom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);
      
      double rad_angle = clipper::Util::d2rad(adj->value);
      clipper::Coord_orth new_c = 
	 g.rotate_round_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			       this_c.second,
			       backbone_torsion_end_ca_1, rad_angle);
      clipper::Coord_orth new_o = 
	 g.rotate_round_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			       this_o.second,
			       backbone_torsion_end_ca_1, rad_angle);
      clipper::Coord_orth new_n = 
	 g.rotate_round_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			       next_n.second,
			       backbone_torsion_end_ca_1, rad_angle);
      n_atom_p->x = new_n.x();
      n_atom_p->y = new_n.y();
      n_atom_p->z = new_n.z();
  
      c_atom_p->x = new_c.x();
      c_atom_p->y = new_c.y();
      c_atom_p->z = new_c.z();
      
      o_atom_p->x = new_o.x();
      o_atom_p->y = new_o.y();
      o_atom_p->z = new_o.z();

      std::pair<std::pair<double, double>, std::pair<double, double> > pp = 
	 g.phi_psi_pairs_from_moving_atoms();
      
//    std::cout << pp.first.first  << " " << pp.first.second << "      " 
// 	     << pp.second.first << " " << pp.second.second << std::endl;

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

      if (edit_phi_psi_plot) { 
	 std::vector <coot::util::phi_psi_t> vp;
	 std::string label = int_to_string(c_atom_p->GetSeqNum());
	 if (pp.first.first > -200) { 
	    label += c_atom_p->GetChainID();
	    coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first), 
					  clipper::Util::rad2d(pp.first.second),
					  "resname", label, 1, "inscode", "chainid");
	    vp.push_back(phipsi1);
	 }
	 if (pp.second.first > -200) { 
	    label = int_to_string(n_atom_p->GetSeqNum());
	    label += n_atom_p->GetChainID();
	    coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first), 
					  clipper::Util::rad2d(pp.second.second),
					  "resname", label, 1, "inscode", "chainid");

	    vp.push_back(phipsi2);
	 }
	 if (vp.size() > 0) 
	    edit_phi_psi_plot->draw_it(vp);
      } 
#endif // HAVE_GTK_CANVAS
      regularize_object_bonds_box.clear_up();
      g.make_moving_atoms_graphics_object(*moving_atoms_asc);
      graphics_draw();

   } else { 
      std::cout << "ERROR:: can't find rama points in edit_backbone_peptide_changed_func" 
		<< std::endl;
   } 
} 

// static 
void
graphics_info_t::edit_backbone_carbonyl_changed_func(GtkAdjustment *adj, GtkWidget *window) { 

   graphics_info_t g;
   // std::cout << "change backbone peptide by: " << adj->value << std::endl;
   
   std::pair<short int, clipper::Coord_orth> this_c = rama_points.get("this_c");
   std::pair<short int, clipper::Coord_orth> this_o = rama_points.get("this_o");

   if (this_c.first && this_o.first) { 
      CAtom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
      CAtom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);
      CAtom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);

      clipper::Coord_orth carbonyl_n_pos(n_atom_p->x, n_atom_p->y, n_atom_p->z);

      double rad_angle = clipper::Util::d2rad(adj->value);
      clipper::Coord_orth new_c = 
	 g.rotate_round_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
			       this_c.second,
			       backbone_torsion_end_ca_1, rad_angle);
      clipper::Coord_orth new_o = 
	 g.rotate_round_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
			       this_o.second,
			       backbone_torsion_end_ca_1, rad_angle);

      c_atom_p->x = new_c.x();
      c_atom_p->y = new_c.y();
      c_atom_p->z = new_c.z();
      
      o_atom_p->x = new_o.x();
      o_atom_p->y = new_o.y();
      o_atom_p->z = new_o.z();

      std::pair<std::pair<double, double>, std::pair<double, double> > pp = 
	 g.phi_psi_pairs_from_moving_atoms();
      
//    std::cout << pp.first.first  << " " << pp.first.second << "      " 
// 	     << pp.second.first << " " << pp.second.second << std::endl;

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

      if (edit_phi_psi_plot) { 
	 std::vector <coot::util::phi_psi_t> vp;
	 std::string label = int_to_string(c_atom_p->GetSeqNum());
	 if (pp.first.first > -200) { 
	    label += c_atom_p->GetChainID();
	    coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first), 
					  clipper::Util::rad2d(pp.first.second),
					  "resname", label, 1, "inscode", "chainid");
	    vp.push_back(phipsi1);
	 }
	 if (pp.second.first > -200) { 
	    label = int_to_string(n_atom_p->GetSeqNum());
	    label += n_atom_p->GetChainID();
	    coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first), 
					  clipper::Util::rad2d(pp.second.second),
					  "resname", label, 1, "inscode", "chainid");

	    vp.push_back(phipsi2);
	 }
	 if (vp.size() > 0) 
	    edit_phi_psi_plot->draw_it(vp);
      } 

#endif // HAVE_GTK_CANVAS
      regularize_object_bonds_box.clear_up();
      g.make_moving_atoms_graphics_object(*moving_atoms_asc);
      graphics_draw();

   } else { 
      std::cout << "ERROR:: can't find rama points in edit_backbone_peptide_changed_func" 
		<< std::endl;
   } 
}


// #include "mmdb-extras.h"

// Tinker with the moving atoms
void
graphics_info_t::change_peptide_carbonyl_by(double angle) { 

//    std::cout << "move carbonyl by " << angle << std::endl;
   CAtom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);
   CAtom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
   CAtom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);

   clipper::Coord_orth carbonyl_n_pos(n_atom_p->x, n_atom_p->y, n_atom_p->z);
   clipper::Coord_orth carbonyl_c_pos(c_atom_p->x, c_atom_p->y, c_atom_p->z);
   clipper::Coord_orth carbonyl_o_pos(o_atom_p->x, o_atom_p->y, o_atom_p->z);

   double rad_angle = clipper::Util::d2rad(angle);
   
   clipper::Coord_orth new_c = 
      rotate_round_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
			  carbonyl_c_pos,
			  carbonyl_n_pos, rad_angle);
   clipper::Coord_orth new_o = 
      rotate_round_vector(carbonyl_n_pos - backbone_torsion_end_ca_1,
			  carbonyl_o_pos,
			  carbonyl_n_pos, rad_angle);

   c_atom_p->x = new_c.x();
   c_atom_p->y = new_c.y();
   c_atom_p->z = new_c.z();

   o_atom_p->x = new_o.x();
   o_atom_p->y = new_o.y();
   o_atom_p->z = new_o.z();

   std::pair<std::pair<double, double>, std::pair<double, double> > pp = 
      phi_psi_pairs_from_moving_atoms();
   
//    std::cout << pp.first.first  << " " << pp.first.second << "      " 
// 	     << pp.second.first << " " << pp.second.second << std::endl;
   

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   if (edit_phi_psi_plot) { 
      std::vector <coot::util::phi_psi_t> vp;
      std::string label = int_to_string(c_atom_p->GetSeqNum());
      label += c_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first), 
				    clipper::Util::rad2d(pp.first.second),
				    "resname", label, 1, "inscode", "chainid");
      label = int_to_string(n_atom_p->GetSeqNum());
      label += n_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first), 
				    clipper::Util::rad2d(pp.second.second),
				    "resname", label, 1, "inscode", "chainid");

      vp.push_back(phipsi1);
      vp.push_back(phipsi2);
      edit_phi_psi_plot->draw_it(vp);
   } 

#endif // HAVE_GTK_CANVAS

   regularize_object_bonds_box.clear_up();
   make_moving_atoms_graphics_object(*moving_atoms_asc);
   graphics_draw();
} 


// Return a pair of phi psi's for the residue of the carbonyl and the next
// residue.  
// 
// If for some reason the pair is incalculable, put phi (first) to -2000.
//
std::pair<std::pair<double, double>, std::pair<double, double> >
graphics_info_t::phi_psi_pairs_from_moving_atoms() {

   std::pair<std::pair<double, double>, std::pair<double, double> > p;

   // There are 2 atoms in the moving atoms that we need for this calculation,
   // this_C and next_N.  The rest we will pull out of a pre-stored vector of
   // pairs (made in execute_setup_backbone_torsion_edit): rama_points.

   CAtom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
   CAtom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);

   clipper::Coord_orth this_c (c_atom_p->x, c_atom_p->y, c_atom_p->z);
   clipper::Coord_orth next_n (n_atom_p->x, n_atom_p->y, n_atom_p->z);

   std::pair<short int, clipper::Coord_orth> prev_c      = rama_points.get("prev_c");
   std::pair<short int, clipper::Coord_orth> this_ca     = rama_points.get("this_ca");
   std::pair<short int, clipper::Coord_orth> this_n      = rama_points.get("this_n");
   std::pair<short int, clipper::Coord_orth> next_ca     = rama_points.get("next_ca");
   std::pair<short int, clipper::Coord_orth> next_c      = rama_points.get("next_c");
   std::pair<short int, clipper::Coord_orth> next_plus_n = rama_points.get("next+1_n");
   
   if (prev_c.first && this_ca.first && this_n.first) { 

      // we can calculate the first ramachadran phi/psi pair
      double phi = clipper::Coord_orth::torsion(prev_c.second, this_n.second,  this_ca.second, this_c);
      double psi = clipper::Coord_orth::torsion(this_n.second, this_ca.second, this_c,         next_n);

      p.first.first  = phi;
      p.first.second = psi;

   } else { 

      // can't get the first ramachandran point
      p.first.first = -2000;

   }

   if (next_ca.first && next_c.first && next_plus_n.first) { 

      // we can calculate the first ramachadran phi/psi pair
      
      double phi = clipper::Coord_orth::torsion(this_c, next_n, next_ca.second, next_c.second);
      double psi = clipper::Coord_orth::torsion(next_n, next_ca.second, next_c.second, next_plus_n.second);

      p.second.first  = phi;
      p.second.second = psi;

   } else { 

      // can't get the second ramachandran point
      p.second.first = -2000;
   } 


   return p;
} 

// Tinker with the moving atoms
void
graphics_info_t::change_peptide_peptide_by(double angle) { 

//    std::cout << "move peptide by " << angle << std::endl;

   CAtom *n_atom_p = coot::get_first_atom_with_atom_name(" N  ", *moving_atoms_asc);
   CAtom *c_atom_p = coot::get_first_atom_with_atom_name(" C  ", *moving_atoms_asc);
   CAtom *o_atom_p = coot::get_first_atom_with_atom_name(" O  ", *moving_atoms_asc);

   clipper::Coord_orth carbonyl_n_pos(n_atom_p->x, n_atom_p->y, n_atom_p->z);
   clipper::Coord_orth carbonyl_c_pos(c_atom_p->x, c_atom_p->y, c_atom_p->z);
   clipper::Coord_orth carbonyl_o_pos(o_atom_p->x, o_atom_p->y, o_atom_p->z);

   double rad_angle = clipper::Util::d2rad(angle);

   clipper::Coord_orth new_c = 
      rotate_round_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			  carbonyl_c_pos,
			  backbone_torsion_end_ca_1, rad_angle);
   clipper::Coord_orth new_o = 
      rotate_round_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			  carbonyl_o_pos,
			  backbone_torsion_end_ca_1, rad_angle);

   clipper::Coord_orth new_n = 
      rotate_round_vector(backbone_torsion_end_ca_2 - backbone_torsion_end_ca_1,
			  carbonyl_n_pos,
			  backbone_torsion_end_ca_1, rad_angle);

   n_atom_p->x = new_n.x();
   n_atom_p->y = new_n.y();
   n_atom_p->z = new_n.z();

   c_atom_p->x = new_c.x();
   c_atom_p->y = new_c.y();
   c_atom_p->z = new_c.z();

   o_atom_p->x = new_o.x();
   o_atom_p->y = new_o.y();
   o_atom_p->z = new_o.z();

   std::pair<std::pair<double, double>, std::pair<double, double> > pp = 
      phi_psi_pairs_from_moving_atoms();
   
//    std::cout << pp.first.first  << " " << pp.first.second << "      " 
// 	     << pp.second.first << " " << pp.second.second << std::endl;

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   if (edit_phi_psi_plot) { 
      std::vector <coot::util::phi_psi_t> vp;
      std::string label = int_to_string(c_atom_p->GetSeqNum());
      label += c_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi1(clipper::Util::rad2d(pp.first.first), 
				    clipper::Util::rad2d(pp.first.second),
				    "resname", label, 1, "inscode", "chainid");
      label = int_to_string(n_atom_p->GetSeqNum());
      label += n_atom_p->GetChainID();
      coot::util::phi_psi_t phipsi2(clipper::Util::rad2d(pp.second.first), 
				    clipper::Util::rad2d(pp.second.second),
				    "resname", label, 1, "inscode", "chainid");

      vp.push_back(phipsi1);
      vp.push_back(phipsi2);
      edit_phi_psi_plot->draw_it(vp);
   } 

#endif // HAVE_GTK_CANVAS

   regularize_object_bonds_box.clear_up();
   make_moving_atoms_graphics_object(*moving_atoms_asc);
   graphics_draw();
} 


void 
graphics_info_t::set_backbone_torsion_peptide_button_start_pos(int ix, int iy) { 
   backbone_torsion_peptide_button_start_pos_x = ix;
   backbone_torsion_peptide_button_start_pos_y = iy;
} 

void 
graphics_info_t::set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy) { 
   backbone_torsion_carbonyl_button_start_pos_x = ix;
   backbone_torsion_carbonyl_button_start_pos_y = iy;
} 


void
graphics_info_t::change_peptide_peptide_by_current_button_pos(int ix, int iy) { 

   double diff = 0.05* (ix - backbone_torsion_peptide_button_start_pos_x);
   change_peptide_peptide_by(diff);
} 

void
graphics_info_t::change_peptide_carbonyl_by_current_button_pos(int ix, int iy) { 

   double diff = 0.05 * (ix - backbone_torsion_carbonyl_button_start_pos_x);
   change_peptide_carbonyl_by(diff);
} 

//      ----------------- sequence view ----------------
#if defined (HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
//
// return NULL on failure to find the class for this molecule
coot::sequence_view *
graphics_info_t::get_sequence_view(int imol) {

   GtkWidget *w = NULL;
   coot::sequence_view *r = NULL;

   if (imol < n_molecules()) {
      if (molecules[imol].has_model()) {
	 // w = sequence_view_is_displayed[imol];
	 w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
	 r = (coot::sequence_view *) gtk_object_get_user_data(GTK_OBJECT(w));
	 // std::cout << "DEBUG:: user data from " << w << " is " << r << std::endl;
      }
   } 
   return r;
}
#endif // HAVE_GTK_CANVAS

void
graphics_info_t::set_sequence_view_is_displayed(GtkWidget *widget, int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

   // first delete the old sequence view if it exists
   GtkWidget *w = coot::get_validation_graph(imol, coot::SEQUENCE_VIEW);
   if (w) {
      coot::sequence_view *sv = (coot::sequence_view *)
	 gtk_object_get_user_data(GTK_OBJECT(w));
      delete sv;
   }

   if (imol < n_molecules()) {
//       coot::sequence_view *sv = (coot::sequence_view *)
// 	 gtk_object_get_user_data(GTK_OBJECT(sequence_view_is_displayed[imol]));
//       std::cout << "DEBUG:: seting sequence_view_is_displayed[" << imol
// 		<< "] " << widget << std::endl;
      // sequence_view_is_displayed[imol] = widget; // ols style
      coot::set_validation_graph(imol, coot::SEQUENCE_VIEW, widget);
   }
#endif // HAVE_GTK_CANVAS   
} 

// ---------------------- geometry objects -----------------------------
void
graphics_info_t::geometry_objects() {

   int ndist = distance_object_vec->size();
   double dist;
   clipper::Coord_orth text_pos;

   if (ndist > 0) {
      glEnable(GL_LINE_STIPPLE);
      glLineStipple (1, 0x00FF);
      glLineWidth(2.0);
      glColor3f(0.5, 0.8, 0.6);

      for (int i=0; i<ndist; i++) {
	 glBegin(GL_LINES);
	 glVertex3d( (*distance_object_vec)[i].first.x(),
		     (*distance_object_vec)[i].first.y(),
		     (*distance_object_vec)[i].first.z());
	 glVertex3d( (*distance_object_vec)[i].second.x(),
		     (*distance_object_vec)[i].second.y(),
		     (*distance_object_vec)[i].second.z());
	 glEnd();
	 text_pos = (*distance_object_vec)[i].first + 
	    0.5 * ( (*distance_object_vec)[i].second - (*distance_object_vec)[i].first + 
		   clipper::Coord_orth(0.0, 0.1, 0.1));
	 glRasterPos3d(text_pos.x(), text_pos.y(), text_pos.z());
	 dist = clipper::Coord_orth::length( (*distance_object_vec)[i].first, (*distance_object_vec)[i].second);
	 printString(float_to_string(dist));

      }
      glDisable(GL_LINE_STIPPLE);
   }

   if (dynamic_distances.size() > 0) {
      draw_dynamic_distances();
   }
}

// static
void
graphics_info_t::draw_dynamic_distances() {

   if (dynamic_distances.size() > 0) {
      glLineWidth(2.0);
      glColor3f(0.5, 0.8, 0.6);

      glEnable(GL_LINE_STIPPLE);
      for (unsigned int i=0; i<dynamic_distances.size(); i++) {
	 dynamic_distances[i].draw_dynamic_distance();
      }
      glDisable(GL_LINE_STIPPLE);
   }
}

void
coot::intermediate_atom_distance_t::draw_dynamic_distance() const {
   
   //glEnable(GL_LINE_STIPPLE);
   glBegin(GL_LINES);
   glVertex3d(dynamic_atom->x,
	      dynamic_atom->y,
	      dynamic_atom->z);
   
   glVertex3d(static_position.x(),
	      static_position.y(),
	      static_position.z());
   glEnd();
   // glDisable(GL_LINE_STIPPLE);

   coot::Cartesian at_pt(dynamic_atom->x,
			 dynamic_atom->y,
			 dynamic_atom->z);

   // The length of the intermediate distance:
   coot::Cartesian text_pos = at_pt + static_position;
   text_pos *= 0.5;
   text_pos += coot::Cartesian(0.0, 0.1, 0.1);
   coot::Cartesian vec_diff = at_pt - static_position;
   float dist = vec_diff.length();
   glRasterPos3d(text_pos.x(), text_pos.y(), text_pos.z());
   std::string t = coot::util::float_to_string(dist);
   printString(t);
}

void
graphics_info_t::pointer_distances_objects() {

   // and pointer distances:
   
   if (pointer_distances_object_vec->size() > 0) {
      double dist;
      clipper::Coord_orth text_pos;
      glEnable(GL_LINE_STIPPLE);
      glLineStipple (1, 0x00FF);
      glLineWidth(2.0);
      glColor3f(0.5, 0.7, 0.8);
      int ndist = pointer_distances_object_vec->size();
      for (int i=0; i<ndist; i++) {
	 glBegin(GL_LINES);
	 glVertex3d( (*pointer_distances_object_vec)[i].first.x(),
		     (*pointer_distances_object_vec)[i].first.y(),
		     (*pointer_distances_object_vec)[i].first.z());
	 glVertex3d( (*pointer_distances_object_vec)[i].second.x(),
		     (*pointer_distances_object_vec)[i].second.y(),
		     (*pointer_distances_object_vec)[i].second.z());
	 glEnd();
	 text_pos = (*pointer_distances_object_vec)[i].first + 
	    0.5 * ( (*pointer_distances_object_vec)[i].second - (*pointer_distances_object_vec)[i].first + 
		   clipper::Coord_orth(0.0, 0.1, 0.1));
	 glRasterPos3d(text_pos.x(), text_pos.y(), text_pos.z());
	 dist = clipper::Coord_orth::length( (*pointer_distances_object_vec)[i].first, (*pointer_distances_object_vec)[i].second);
	 printString(float_to_string(dist));
      }
      glDisable(GL_LINE_STIPPLE);
   }
} 

void
graphics_info_t::make_pointer_distance_objects() {

   clipper::Coord_orth cen(rotation_centre_x,
			   rotation_centre_y,
			   rotation_centre_z);

   std::vector<clipper::Coord_orth> distances;
   std::vector<clipper::Coord_orth> mol_distances;

   if (show_pointer_distances_flag) { 
      for (int imol=0; imol<n_molecules(); imol++) {
	 if (molecules[imol].has_model()) {
	    if (molecules[imol].is_displayed_p()) { 
	       mol_distances = molecules[imol].distances_to_point(cen,
								  pointer_min_dist,
								  pointer_max_dist);
	       if (mol_distances.size() > 0) {
		  // append
		  for (unsigned int id=0; id<mol_distances.size(); id++) { 
		     distances.push_back(mol_distances[id]);
		  }
	       }
	    }
	 }
      }
   
      pointer_distances_object_vec->resize(0);
      for (unsigned int id=0; id<distances.size(); id++) {
	 pointer_distances_object_vec->push_back(std::pair<clipper::Coord_orth, clipper::Coord_orth> (distances[id], cen));
      }
   }
}

void
graphics_info_t::clear_pointer_distances() {

   pointer_distances_object_vec->resize(0);
   graphics_draw();

}

void
graphics_info_t::clear_simple_distances() {

   distance_object_vec->resize(0);
   graphics_draw();
}

void
graphics_info_t::clear_last_simple_distance() {

   int n = distance_object_vec->size();
   if (n > 0) {
      distance_object_vec->resize(n-1);
      graphics_draw();
   }
}


// ---------------------- generic objects -----------------------------

void
coot::generic_display_object_t::add_line(const coot::colour_holder &colour_in,
					 const std::string &colour_name,
					 const int &width_in, 
					 const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords_in) {

   int lines_set_index = -1; // magic unset value
   
   for (unsigned int ils=0; ils<lines_set.size(); ils++) {
      if (lines_set[ils].colour_name == colour_name) {
	 if (lines_set[ils].width == width_in) {
	    lines_set_index = ils;
	    break;
	 }
      }
   }

   if (lines_set_index == -1) {
      generic_display_line_set_t t(colour_in, colour_name, width_in);
      lines_set.push_back(t);
      lines_set_index = lines_set.size() -1;
   }

   coot::generic_display_line_t line(coords_in);
   lines_set[lines_set_index].add_line(line);
   
}

void coot::generic_display_object_t::add_point(const coot::colour_holder &colour_in,
					       const std::string &colour_name,
					       const int &size_in, 
					       const clipper::Coord_orth &coords_in) {
   int points_set_index = -1; // magic unset number
   for (unsigned int ips=0; ips<points_set.size(); ips++) {
      if (points_set[ips].colour_name == colour_name) {
	 if (points_set[ips].size == size_in) {
	    points_set_index = ips;
	    break;
	 } 
      }
   }
   if (points_set_index == -1) {
      coot::generic_display_point_set_t point_set(colour_in, colour_name, size_in);
      points_set.push_back(point_set);
      points_set_index = points_set.size() -1;
   }
   points_set[points_set_index].add_point(coords_in);
}

// static
void
graphics_info_t::draw_generic_objects() {

   // std::cout << "debug:: drawing " << generic_objects_p->size() << " generic objects" << std::endl;
   for (unsigned int i=0; i<generic_objects_p->size(); i++) {

      //       std::cout << "debug:: drawing generic object - outer " << i << std::endl;
      if ((*generic_objects_p)[i].is_displayed_flag) {

	 // std::cout << "debug:: drawing generic  " << (*generic_objects_p)[i].lines_set.size()
	 // << " lines " << std::endl;

	 // Lines
	 for (unsigned int ils=0; ils< (*generic_objects_p)[i].lines_set.size(); ils++) {
	    glLineWidth((*generic_objects_p)[i].lines_set[ils].width);
	    glColor3f((*generic_objects_p)[i].lines_set[ils].colour.red,
		      (*generic_objects_p)[i].lines_set[ils].colour.green,
		      (*generic_objects_p)[i].lines_set[ils].colour.blue);
	    glBegin(GL_LINES);
	    unsigned int s = (*generic_objects_p)[i].lines_set[ils].lines.size();
	    for (unsigned int iline=0; iline<s; iline++) {
	       glVertex3f((*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first.x(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first.y(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first.z());
	       glVertex3f((*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second.x(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second.y(),
			  (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second.z());
	    }
	    glEnd();
	 }

	 // Points
	 for (unsigned int ips=0; ips<(*generic_objects_p)[i].points_set.size(); ips++) {
	    glPointSize((*generic_objects_p)[i].points_set[ips].size);
	    glColor3f((*generic_objects_p)[i].points_set[ips].colour.red,
		      (*generic_objects_p)[i].points_set[ips].colour.green,
		      (*generic_objects_p)[i].points_set[ips].colour.blue);
	    glBegin(GL_POINTS);
	    unsigned int npoints = (*generic_objects_p)[i].points_set[ips].points.size();
	    for (unsigned int ipoint=0; ipoint<npoints; ipoint++) { 
	       glVertex3f((*generic_objects_p)[i].points_set[ips].points[ipoint].x(),
			  (*generic_objects_p)[i].points_set[ips].points[ipoint].y(),
			  (*generic_objects_p)[i].points_set[ips].points[ipoint].z());
	    }
	    glEnd();
	 }

         // Display lists
	 for (unsigned int idl=0; idl<(*generic_objects_p)[i].GL_display_list_handles.size(); idl++) {
             glCallList((*generic_objects_p)[i].GL_display_list_handles[idl]);
         }
      }
   }
}

void
graphics_info_t::draw_generic_text() {

   // should be const, I think.

   if (generic_texts_p->size() > 0 ) { 
      GLfloat pink[3] =  { 1.0, 0.8, 0.8 };
      glColor3fv(pink);
      glPushAttrib (GL_LIST_BIT);
      void *font = graphics_info_t::atom_label_font;
      font = GLUT_BITMAP_TIMES_ROMAN_24;
      for (unsigned int i=0; i<generic_texts_p->size(); i++) {
	 
	 glRasterPos3f((*generic_texts_p)[i].x, (*generic_texts_p)[i].y, (*generic_texts_p)[i].z);
	 for (unsigned int is = 0; is < (*generic_texts_p)[i].s.length(); is++)
	    glutBitmapCharacter (font, (*generic_texts_p)[i].s[is]);
      }
      glPopAttrib ();
   }
} 

// static
coot::colour_holder
coot::generic_display_object_t::colour_values_from_colour_name(const std::string &c) {

   coot::colour_holder colour; 
   colour.red = 0.4; 
   colour.green = 0.4; 
   colour.blue = 0.4;

   if (c == "blue") {
      colour.red = 0.1; 
      colour.green = 0.1; 
      colour.blue = 0.8;
   } else {
      if (c == "sky") {
	 colour.red = 0.4; 
	 colour.green = 0.4; 
	 colour.blue = 0.6;
      } else {
	 if (c == "green") {
	    colour.red   = 0.1; 
	    colour.green = 0.8; 
	    colour.blue  = 0.1;
	 } else {
	    if (c == "greentint") {
	       colour.red = 0.4; 
	       colour.green = 0.65; 
	       colour.blue = 0.4;
	    } else { 
	       if (c == "sea") {
		  colour.red = 0.1; 
		  colour.green = 0.6; 
		  colour.blue = 0.6;
	       } else {
		  if (c == "yellow") {
		     colour.red = 0.8; 
		     colour.green = 0.8; 
		     colour.blue = 0.0;
		  } else {
		     if (c == "yellowtint") {
			colour.red = 0.65; 
			colour.green = 0.65; 
			colour.blue = 0.4;
		     } else {
			if (c == "orange") {
			   colour.red = 0.9; 
			   colour.green = 0.6; 
			   colour.blue = 0.1;
			} else {
			   if (c == "red") {
			      colour.red = 0.9; 
			      colour.green = 0.1; 
			      colour.blue = 0.1;
			   } else {
			      if (c == "hotpink") {
				 colour.red = 0.9; 
				 colour.green = 0.2; 
				 colour.blue = 0.6;
			      } else {
				 if (c == "cyan") {
				    colour.red = 0.1; 
				    colour.green = 0.7; 
				    colour.blue = 0.7;
				 } else {
				    if (c == "aquamarine") {
				       colour.red = 0.1; 
				       colour.green = 0.8; 
				       colour.blue = 0.6;
				    } else {
				       if (c == "forestgreen") {
					  colour.red = 0.6; 
					  colour.green = 0.8; 
					  colour.blue = 0.1;
				       } else {
					  if (c == "yellowgreen") {
					     colour.red   = 0.6; 
					     colour.green = 0.8; 
					     colour.blue  = 0.2;
					  } else {
					     if (c == "goldenrod") {
						colour.red   = 0.85; 
						colour.green = 0.65; 
						colour.blue  = 0.12;
					     } else {
						if (c == "orangered") {
						   colour.red   = 0.9; 
						   colour.green = 0.27; 
						   colour.blue  = 0.0;
						}
					     }
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 } 
      } 
   }

//    std::cout << "debug:: in colour_values_from_colour_name from colour " << c
// 	     << " we assign colour values "
// 	     << colour[0] << " " 
// 	     << colour[1] << " " 
// 	     << colour[2] << "\n";
   return colour;
} 


// ----------------- done generic objects -----------------------------

void 
graphics_info_t::unset_geometry_dialog_distance_togglebutton() { 

   if (geometry_dialog) { 
      GtkWidget *toggle_button = lookup_widget(geometry_dialog, 
					       "geometry_distance_togglebutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
   }
}

void
graphics_info_t::unset_geometry_dialog_dynamic_distance_togglebutton() {

   if (geometry_dialog) { 
      GtkWidget *toggle_button = lookup_widget(geometry_dialog, 
					       "geometry_dynamic_distance_togglebutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
   }
}


void
graphics_info_t::unset_geometry_dialog_angle_togglebutton() { 

   if (geometry_dialog) { 
      GtkWidget *toggle_button = lookup_widget(geometry_dialog, 
					       "geometry_angle_togglebutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
   }
}


void
graphics_info_t::unset_geometry_dialog_torsion_togglebutton() { 

   if (geometry_dialog) { 
      GtkWidget *toggle_button = lookup_widget(geometry_dialog, 
					       "geometry_torsion_togglebutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), FALSE);
   }
}


void
graphics_info_t::remove_all_atom_labels() { 

   for (int i=0; i<n_molecules(); i++) { 
      if (molecules[i].has_model()) { 
	 molecules[i].remove_atom_labels();
      } 
   } 
   graphics_draw();

} 


// static
std::string
graphics_info_t::ccp4_defs_file_name() {

#if defined WIN32
#ifdef WINDOWS_MINGW
// BL says:: in my windows it's found in $USERPROFILE
// would guess that's true for other win32 too....
    char *home = getenv("USERPROFILE");
#else
    char *home = getenv("HOMEPATH");
#endif // WINDOWS_MINGW
#else
    char *home = getenv("HOME");
#endif // WIN32

#if defined(WINDOWS_MINGW)|| defined(_MSC_VER)
    std::string path = "/CCP4/windows/directories.def";
#else    
    std::string path = "/.CCP4/unix/directories.def";
#endif     
   std::string filename = home + path;

   return filename;
}


void
graphics_info_t::set_zoom_adjustment(GtkWidget *widget) { 

   GtkWidget *zoom_hscale = lookup_widget(widget, "zoom_hscale");

   //  GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -18.0, 36.0, 0.01, 1.0, 18));
   //  GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -2.0, 4.0, 0.01, 1.0, 2));

   GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(zoom, zoom*0.125, zoom*8, 0.01, 0.5, zoom));

   gtk_range_set_adjustment(GTK_RANGE(zoom_hscale), adj);
   gtk_signal_connect(GTK_OBJECT(adj), "value_changed",
		      GTK_SIGNAL_FUNC(graphics_info_t::zoom_adj_changed),
		      NULL);
}


// static
void
graphics_info_t::zoom_adj_changed(GtkAdjustment *adj, GtkWidget *window) { 
   graphics_info_t g;

//    double scaled = pow(M_E, -adj->value);
//    std::cout <<  adj->value << " " << scaled << std::endl;

   // g.zoom *= scaled;
   g.zoom = adj->value;
   graphics_draw();

} 

void
graphics_info_t::add_coordinates_glob_extension(const std::string &extension) { 

   coordinates_glob_extensions->push_back(extension);
} 
  

void
graphics_info_t::add_data_glob_extension(const std::string &extension) { 
   data_glob_extensions->push_back(extension);

} 

void
graphics_info_t::add_map_glob_extension(const std::string &extension) { 
   map_glob_extensions->push_back(extension);

} 

void
graphics_info_t::add_dictionary_glob_extension(const std::string &extension) { 
   dictionary_glob_extensions->push_back(extension);

} 

void
graphics_info_t::remove_coordinates_glob_extension(const std::string &extension) { 

  std::vector<std::string>::iterator it;
  for (it = coordinates_glob_extensions->begin(); it<coordinates_glob_extensions->end(); it++) {
    if (*it == extension) {
      coordinates_glob_extensions->erase(it);
      // could put in break here!?
      // avoid since it could happen that you have multiples of same entry
    }
  }
} 
  

void
graphics_info_t::remove_data_glob_extension(const std::string &extension) { 

  std::vector<std::string>::iterator it;
  for (it = data_glob_extensions->begin(); it<data_glob_extensions->end(); it++) {
    if (*it == extension) {
      data_glob_extensions->erase(it);
    }
  }
} 

void
graphics_info_t::remove_map_glob_extension(const std::string &extension) { 

  std::vector<std::string>::iterator it;
  for (it = map_glob_extensions->begin(); it<map_glob_extensions->end(); it++) {
    if (*it == extension) {
      map_glob_extensions->erase(it);
    }
  }
} 

void
graphics_info_t::remove_dictionary_glob_extension(const std::string &extension) { 

  std::vector<std::string>::iterator it;
  for (it = dictionary_glob_extensions->begin(); it<dictionary_glob_extensions->end(); it++) {
    if (*it == extension) {
      dictionary_glob_extensions->erase(it);
    }
  }
}

// ----------------------------------------------------------------------------
//                check waters 
// ----------------------------------------------------------------------------
// 
void
graphics_info_t::check_waters_by_difference_map(int imol_waters, int imol_diff_map, 
						int interactive_flag) {
   
   if (imol_waters < n_molecules()) {
      if (imol_diff_map < n_molecules()) {
	 if (molecules[imol_waters].has_model()) {
	    if (molecules[imol_diff_map].has_map()) {
	       if (molecules[imol_diff_map].is_difference_map_p()) {
		  std::vector <coot::atom_spec_t> v = molecules[imol_waters].check_waters_by_difference_map(molecules[imol_diff_map].xmap_list[0], check_waters_by_difference_map_sigma_level);
		  if (interactive_flag) { 
		     GtkWidget *w = wrapped_create_checked_waters_by_variance_dialog(v, imol_waters);
		     gtk_widget_show(w);
		  }
	       } else {
		  std::cout << "molecule " <<  imol_diff_map
			    << " is not a difference map\n";
	       }
	    } else {
	       std::cout << "molecule " <<  imol_diff_map << "has no map\n";
	    }
	 } else {
	    std::cout << "molecule " <<  imol_waters << "has no model\n";
	 } 
      } else {
	 std::cout << "no molecule for difference map\n";
      }
   } else {
      std::cout << "no molecule for difference (water) coordinates\n";
   } 
}

// results widget:
// 
// 

GtkWidget *
graphics_info_t::wrapped_create_checked_waters_by_variance_dialog(const std::vector <coot::atom_spec_t> &v, int imol) {

   GtkWidget *w;

   if (v.size() > 0) {
      w = create_interesting_waters_by_difference_map_check_dialog();
      GtkWidget *vbox = lookup_widget(w, "interesting_waters_by_difference_map_check_vbox");
      GtkWidget *button;
      coot::atom_spec_t *atom_spec;
      GSList *gr_group = NULL;
      
      for (unsigned int i=0; i<v.size(); i++) {

	 std::cout << "Suspicious water: "
		   << v[i].atom_name
		   << v[i].alt_conf << " "
		   << v[i].resno << " "
		   << v[i].insertion_code << " "
		   << v[i].chain << "\n";

	 std::string button_label(" ");
	 button_label += v[i].chain;
	 button_label += " " ;
	 button_label += int_to_string(v[i].resno);
	 button_label += " " ;
	 button_label += v[i].atom_name;
	 button_label += " " ;
	 button_label += v[i].alt_conf;
	 button_label += " " ;

	 button = gtk_radio_button_new_with_label(gr_group, button_label.c_str());
	 gr_group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
	 atom_spec = new coot::atom_spec_t(v[i]);
	 atom_spec->int_user_data = imol;
      
	 gtk_signal_connect(GTK_OBJECT(button), "clicked",
			    GTK_SIGNAL_FUNC (on_generic_atom_spec_button_clicked),
			    atom_spec);

	 GtkWidget *frame = gtk_frame_new(NULL);
	 gtk_container_add(GTK_CONTAINER(frame), button);
	 gtk_box_pack_start(GTK_BOX(vbox), frame, FALSE, FALSE, 0);
	 gtk_container_set_border_width(GTK_CONTAINER(frame), 2);
	 gtk_widget_show(button);
	 gtk_widget_show(frame);
      } 
   } else {
      std::cout << "There are no unusual waters\n";
      w = wrapped_nothing_bad_dialog("There were no strange/anomalous waters");
   }
   return w;
}


// static 
void
graphics_info_t::on_generic_atom_spec_button_clicked (GtkButton *button,
						      gpointer user_data) {

   if (GTK_TOGGLE_BUTTON(button)->active) { 

      graphics_info_t g;
      coot::atom_spec_t *atom_spec = (coot::atom_spec_t *) user_data;
//       std::cout << "atom_spec: " 
// 		<< atom_spec->chain << " " << atom_spec->resno << " " << atom_spec->atom_name
// 		<< std::endl;
   
      g.set_go_to_atom_molecule(atom_spec->int_user_data);
      g.set_go_to_atom_chain_residue_atom_name(atom_spec->chain.c_str(),
					       atom_spec->resno,
					       atom_spec->atom_name.c_str(),
					       atom_spec->alt_conf.c_str());
      g.try_centre_from_new_go_to_atom();
      g.update_things_on_move_and_redraw();
   }
}

void
graphics_info_t::check_chiral_volumes(int imol) {

   if (imol < n_molecules()) {
      if (molecules[imol].has_model()) {
	 // return a pair: first is the residues for which no
	 // restraints were found second is a vector of atom specs
	 // that violate chiral volume constraint.
	 std::pair<std::vector<std::string>, std::vector <coot::atom_spec_t> > v =
	    molecules[imol].bad_chiral_volumes();
	 GtkWidget *w = wrapped_check_chiral_volumes_dialog(v.second, imol);
	 if (w) 
	    gtk_widget_show(w);
	 if (v.first.size() != 0) { // bad, there was at least one residue not found in dic.
	    GtkWidget *w = wrapped_create_chiral_restraints_problem_dialog(v.first);
	    gtk_widget_show(w);
	 }
      }
   }
}

GtkWidget *
graphics_info_t::wrapped_create_chiral_restraints_problem_dialog(const std::vector<std::string> &sv) const {

   GtkWidget *w = create_chiral_restraints_problem_dialog();

   GtkWidget *label = lookup_widget(w, "chiral_volume_restraints_problem_label");
   std::string s = "\n   Problem finding restraints for the following residues:   \n\n";
   for (unsigned int i=0; i<sv.size(); i++) {
      s += sv[i];
      s += "  ";
      if (10*((i+1)/10) == (i+1))
	 s += "\n";
   }
   s += "\n";
   gtk_label_set_text(GTK_LABEL(label), s.c_str());
   return w;
}


GtkWidget *
graphics_info_t::wrapped_check_chiral_volumes_dialog(const std::vector <coot::atom_spec_t> &v,
						     int imol) {

   GtkWidget *w = NULL;
   
   std::cout  << "There were " << v.size() << " bad chiral volumes: " << std::endl;

   if (v.size() > 0) {
      GtkWidget *button;
      w = create_bad_chiral_volumes_dialog ();
      GtkWidget *bad_chiral_volume_atom_vbox = 
	 lookup_widget(w, "chiral_volume_baddies_vbox");
      coot::atom_spec_t *atom_spec;
      for (unsigned int i=0; i<v.size(); i++) { 
	 std::cout << "  "
		   << v[i].chain << " " 
		   << v[i].resno << " " 
		   << v[i].atom_name << " " 
		   << v[i].alt_conf << " " 
		   << "\n";

	 // c.f. how we add rotamers: (fill_rotamer_selection_buttons)
	 std::string button_label(" ");
	 button_label += v[i].chain;
	 button_label += " " ;
	 button_label += int_to_string(v[i].resno);
	 button_label += " " ;
	 button_label += v[i].atom_name;
	 button_label += " " ;
	 button_label += v[i].alt_conf;
	 button_label += " " ;

	 button = gtk_button_new_with_label(button_label.c_str());
	 atom_spec = new coot::atom_spec_t(v[i]);
	 atom_spec->int_user_data = imol;
      
	 gtk_signal_connect(GTK_OBJECT(button), "clicked",
			    GTK_SIGNAL_FUNC (on_bad_chiral_volume_button_clicked),
			    atom_spec);

	 gtk_box_pack_start(GTK_BOX(bad_chiral_volume_atom_vbox),
			    button, FALSE, FALSE, 0);
	 gtk_container_set_border_width(GTK_CONTAINER(button), 2);
	 gtk_widget_show(button);
      }

   } else { 
      std::cout << "Congratulations: there are no bad chiral volumes in this molecule.\n";
      w = create_no_bad_chiral_volumes_dialog();
   } 
   return w;
}

// static 
void
graphics_info_t::on_bad_chiral_volume_button_clicked (GtkButton       *button,
						      gpointer         user_data) {


   // This function may get called for an "unclick" too.  It may need
   // an active? test like the generic atom spec button callback.  But
   // perhaps not.
   
   graphics_info_t g;
   coot::atom_spec_t *atom_spec = (coot::atom_spec_t *) user_data;
//    std::cout << "atom_spec: " 
// 	     << atom_spec->chain << " "
// 	     << atom_spec->resno << " "
// 	     << atom_spec->atom_name << " "
// 	     << atom_spec->alt_conf << " "
// 	     << std::endl;
   
   g.set_go_to_atom_molecule(atom_spec->int_user_data);
   g.set_go_to_atom_chain_residue_atom_name(atom_spec->chain.c_str(),
					    atom_spec->resno,
					    atom_spec->atom_name.c_str(),
					    atom_spec->alt_conf.c_str());
   g.try_centre_from_new_go_to_atom();
   g.update_things_on_move_and_redraw();


}


void
graphics_info_t::set_moving_atoms(atom_selection_container_t asc,
				  int imol, int new_coords_type) {

   imol_moving_atoms = imol;
   make_moving_atoms_graphics_object(asc);
   moving_atoms_asc_type = new_coords_type;
} 

//   static
void graphics_info_t::bond_parameters_molecule_menu_item_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t g;
   g.bond_parameters_molecule = pos;
   GtkWidget *w = lookup_widget(GTK_WIDGET(item), "bond_parameters_dialog");
   fill_bond_parameters_internals(w, pos);

} 

// static
//
// imol can be -1 (for no molecules available)
// 
void graphics_info_t::fill_bond_parameters_internals(GtkWidget *w,
						    int imol) {

   // We don't do these 2 any more because they have been moved to
   // bond colour dialog:
   // fill the colour map rotation entry
   // check the Carbons only check button

   // fill the molecule bond width option menu

   // check the draw hydrogens check button

   // and now also set the adjustment on the hscale

   graphics_info_t g;

   GtkWidget *bond_width_option_menu = lookup_widget(w, "bond_parameters_bond_width_optionmenu");
   GtkWidget *draw_hydrogens_yes_radiobutton  = lookup_widget(w, "draw_hydrogens_yes_radiobutton");
   GtkWidget *draw_hydrogens_no_radiobutton   = lookup_widget(w, "draw_hydrogens_no_radiobutton");
   GtkWidget *draw_ncs_ghosts_yes_radiobutton = lookup_widget(w, "draw_ncs_ghosts_yes_radiobutton");
   GtkWidget *draw_ncs_ghosts_no_radiobutton  = lookup_widget(w, "draw_ncs_ghosts_no_radiobutton");

   // bye bye entry
//    gtk_entry_set_text(GTK_ENTRY(entry),
// 		      float_to_string(rotate_colour_map_on_read_pdb).c_str());

   g.bond_thickness_intermediate_value = -1;
   
   // Fill the bond width option menu.
   // Put a redraw on the menu item activate callback.
   // We do the thing with the new menu for the option_menu
   // 
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(bond_width_option_menu));
   GtkSignalFunc signal_func = GTK_SIGNAL_FUNC(graphics_info_t::bond_width_item_select);
   if (menu) 
      gtk_widget_destroy(menu);
   menu = gtk_menu_new();
   GtkWidget *menu_item;
   int current_bond_width = 3;
   if (imol >= 0 ) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    current_bond_width = molecules[imol].bond_thickness();
	 }
      }
   }

   for (int i=1; i<21; i++) {
      std::string s = int_to_string(i);
      menu_item = gtk_menu_item_new_with_label(s.c_str());
      gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
			 GTK_SIGNAL_FUNC(signal_func),
			 GINT_TO_POINTER(i));
      gtk_menu_append(GTK_MENU(menu), menu_item);
      gtk_widget_show(menu_item);
      if (i == current_bond_width) {
	 gtk_menu_set_active(GTK_MENU(menu), i);
      }
   }
   gtk_menu_set_active(GTK_MENU(menu), current_bond_width-1); // 0 offset
   gtk_option_menu_set_menu(GTK_OPTION_MENU(bond_width_option_menu), menu);

   // Draw Hydrogens?
   if (imol >= 0 ) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    if (molecules[imol].draw_hydrogens()) {
	       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_hydrogens_yes_radiobutton), TRUE);
	    } else {
	       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_hydrogens_no_radiobutton), TRUE);
	    }
	 }
      }
   }

   // Draw NCS ghosts?
   if (imol >= 0 ) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    if (molecules[imol].draw_ncs_ghosts_p()) {
	       if (molecules[imol].ncs_ghosts_have_rtops_p()) { 
		  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_ncs_ghosts_yes_radiobutton), TRUE);
	       } else {
		  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_ncs_ghosts_no_radiobutton), TRUE);
	       }
	    } else {
	       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(draw_ncs_ghosts_no_radiobutton), TRUE);
	    }
	 }
      }
   }
   // Make the frame be insensitive if there is no NCS.
   GtkWidget *frame = lookup_widget(w, "ncs_frame");
   short int make_insensitive = 1;
   if (imol >= 0 ) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    if (molecules[imol].has_ncs_p()) {
	       make_insensitive = 0;
	    } else {
	       std::cout << "INFO:: in fill_bond_parameters_internals no NCS for  "
			 << imol << "\n"; 
	    }
	 } else {
	    std::cout << "ERROR:: bad imol in fill_bond_parameters_internals no model "
		      << imol << "\n"; 
	 }
      } else {
	 std::cout << "ERROR:: bad imol in fill_bond_parameters_internals i " << imol << "\n"; 
      }
   } else {
      std::cout << "ERROR:: bad imol in fill_bond_parameters_internals " << imol << "\n"; 
   }
   if (make_insensitive)
      gtk_widget_set_sensitive(frame, FALSE);
   else 
      gtk_widget_set_sensitive(frame, TRUE);
   
}


void
graphics_info_t::fill_bond_colours_dialog_internal(GtkWidget *w) {

   // First the (global) step adjustment:
   GtkScale *hscale = GTK_SCALE(lookup_widget(w, "bond_parameters_colour_rotation_hscale"));
   GtkAdjustment *adjustment = GTK_ADJUSTMENT
      (gtk_adjustment_new(rotate_colour_map_on_read_pdb, 0.0, 370.0, 1.0, 20.0, 10.1));
   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   gtk_signal_connect(GTK_OBJECT(adjustment), "value_changed",
		      GTK_SIGNAL_FUNC(bond_parameters_colour_rotation_adjustment_changed), NULL);


   // Now the "C only" checkbutton:
   GtkWidget *checkbutton = lookup_widget(w, "bond_parameters_rotate_colour_map_c_only_checkbutton");
   if (rotate_colour_map_on_read_pdb_c_only_flag) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   }


   // Now the tricky bit, fill the scrolled vbox of molecule colour rotation step sliders:

   GtkWidget *frame_molecule_N;
   GtkWidget *coords_colour_control_dialog = w;
   GtkWidget *coords_colours_vbox = lookup_widget(w, "coords_colours_vbox");
   GtkWidget *hbox136;
   GtkWidget *label269;
   GtkWidget *label270;
   GtkWidget *coords_colour_hscale_mol_N;
   
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model()) { 

	 std::string m = "Molecule ";
	 m += coot::util::int_to_string(imol);
	 m += " ";
	 m += molecules[imol].name_for_display_manager();
	 frame_molecule_N = gtk_frame_new (m.c_str());
	 gtk_widget_ref (frame_molecule_N);
	 gtk_object_set_data_full (GTK_OBJECT (coords_colour_control_dialog),
				   "frame_molecule_N", frame_molecule_N,
				   (GtkDestroyNotify) gtk_widget_unref);
	 gtk_box_pack_start (GTK_BOX (coords_colours_vbox), frame_molecule_N, TRUE, TRUE, 0);
	 gtk_widget_set_usize (frame_molecule_N, 171, -2);
	 gtk_container_set_border_width (GTK_CONTAINER (frame_molecule_N), 6);

	 hbox136 = gtk_hbox_new (FALSE, 0);
	 gtk_widget_ref (hbox136);
	 gtk_object_set_data_full (GTK_OBJECT (coords_colour_control_dialog), "hbox136", hbox136,
				   (GtkDestroyNotify) gtk_widget_unref);
	 gtk_widget_show (hbox136);
	 gtk_container_add (GTK_CONTAINER (frame_molecule_N), hbox136);

	 label269 = gtk_label_new (_("    "));
	 gtk_widget_ref (label269);
	 gtk_object_set_data_full (GTK_OBJECT (coords_colour_control_dialog), "label269", label269,
				   (GtkDestroyNotify) gtk_widget_unref);
	 gtk_widget_show (label269);
	 gtk_box_pack_start (GTK_BOX (hbox136), label269, FALSE, FALSE, 0);

	 GtkAdjustment *adjustment_mol = GTK_ADJUSTMENT
	    (gtk_adjustment_new(molecules[imol].bonds_colour_map_rotation,
				0.0, 370.0, 1.0, 20.0, 10.1));
	 coords_colour_hscale_mol_N = gtk_hscale_new (adjustment_mol);
	 gtk_range_set_adjustment(GTK_RANGE(coords_colour_hscale_mol_N), adjustment_mol);
	 gtk_signal_connect(GTK_OBJECT(adjustment_mol), "value_changed",
			    GTK_SIGNAL_FUNC(bonds_colour_rotation_adjustment_changed), NULL);
	 gtk_object_set_user_data(GTK_OBJECT(adjustment_mol), GINT_TO_POINTER(imol));
	 
	 gtk_widget_ref (coords_colour_hscale_mol_N);
	 gtk_object_set_data_full (GTK_OBJECT (coords_colour_control_dialog),
				   "coords_colour_hscale_mol_N",
				   coords_colour_hscale_mol_N,
				   (GtkDestroyNotify) gtk_widget_unref);
	 gtk_widget_show (coords_colour_hscale_mol_N);
	 gtk_box_pack_start (GTK_BOX (hbox136), coords_colour_hscale_mol_N, TRUE, TRUE, 0);

	 label270 = gtk_label_new (_("  degrees  "));
	 gtk_widget_ref (label270);
	 gtk_object_set_data_full (GTK_OBJECT (coords_colour_control_dialog), "label270", label270,
				   (GtkDestroyNotify) gtk_widget_unref);
	 gtk_widget_show (label270);
	 gtk_box_pack_start (GTK_BOX (hbox136), label270, FALSE, FALSE, 0);
	 gtk_misc_set_alignment (GTK_MISC (label270), 0.5, 0.56);

	 gtk_widget_show(frame_molecule_N);
      }
   }
   
}

// static 
void graphics_info_t::bond_parameters_colour_rotation_adjustment_changed(GtkAdjustment *adj,
									 GtkWidget *window) {

   graphics_info_t g;
   g.rotate_colour_map_on_read_pdb = adj->value;
   graphics_draw(); // unnecessary.
   
}


// static 
void graphics_info_t::bonds_colour_rotation_adjustment_changed(GtkAdjustment *adj,
							       GtkWidget *window) {
   int imol = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(adj)));

   if (molecules[imol].has_model()) {
      molecules[imol].bonds_colour_map_rotation = adj->value;
   }
   graphics_draw();
   
}

// static
void
graphics_info_t::bond_width_item_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t g;
   g.bond_thickness_intermediate_value = pos;
   if (g.bond_thickness_intermediate_value > 0) {
      int imol = g.bond_parameters_molecule;
      if (is_valid_model_molecule(imol)) 
	 g.set_bond_thickness(imol, g.bond_thickness_intermediate_value);
   }
}


void
graphics_info_t::fill_add_OXT_dialog_internal(GtkWidget *widget, int imol) {

   GtkWidget *chain_optionmenu = lookup_widget(widget, "add_OXT_chain_optionmenu");
//    GtkWidget *at_c_terminus_radiobutton =
//       lookup_widget(widget, "add_OXT_c_terminus_radiobutton");
   GtkSignalFunc signal_func = GTK_SIGNAL_FUNC(graphics_info_t::add_OXT_chain_menu_item_activate);

   std::string a = fill_chain_option_menu(chain_optionmenu, imol, signal_func);
   if (a != "no-chain") {
      graphics_info_t::add_OXT_chain = a;
   }
}

// static (menu item ativate callback)
void
graphics_info_t::add_OXT_chain_menu_item_activate (GtkWidget *item,
						   GtkPositionType pos) {

   char *data = NULL;
   data = (char *)pos;
   if (data)
      add_OXT_chain = std::string(data);
}

// a copy of the c-interface function which does not pass the signal
// function.  We also return the string at the top of the list:
// (return "no-chain" if it was not assigned (nothing in the list)).
//
// static
// 
std::string 
graphics_info_t::fill_chain_option_menu(GtkWidget *chain_option_menu, int imol,
					GtkSignalFunc signal_func) {

   std::string r("no-chain");

   if (imol<graphics_info_t::n_molecules()) {
      if (imol >= 0) { 
	 if (graphics_info_t::molecules[imol].has_model()) {
	    std::vector<std::string> chains = coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
	    GtkWidget *menu_item;
	    GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(chain_option_menu));
	    if (menu)
	       gtk_widget_destroy(menu);
	    menu = gtk_menu_new();
	    short int first_chain_set_flag = 0;
	    std::string first_chain;
	    for (unsigned int i=0; i<chains.size(); i++) {
	       if (first_chain_set_flag == 0) {
		  first_chain_set_flag = 1;
		  first_chain = chains[i];
	       }
	       menu_item = gtk_menu_item_new_with_label(chains[i].c_str());
	       char *v = new char[chains[i].length() + 1];
	       strcpy(v, chains[i].c_str());
	       gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
				  signal_func, v);
	       gtk_menu_append(GTK_MENU(menu), menu_item);
	       gtk_widget_show(menu_item);
	    }
	    /* Link the new menu to the optionmenu widget */
	    gtk_option_menu_set_menu(GTK_OPTION_MENU(chain_option_menu), menu);

	    if (first_chain_set_flag) {
	       r = first_chain;
	    }
	 }
      }
   }
   return r;
} 

// static
void
graphics_info_t::add_OXT_molecule_item_select(GtkWidget *item,
					      GtkPositionType pos) {

   graphics_info_t g;
   g.add_OXT_molecule = pos;
   GtkWidget *w = lookup_widget(item, "add_OXT_dialog");
   g.fill_add_OXT_dialog_internal(w,pos);

}


// static
void
graphics_info_t::fill_renumber_residue_range_dialog(GtkWidget *window) {

   graphics_info_t g;

   GtkWidget *molecule_option_menu =
      lookup_widget(window, "renumber_residue_range_molecule_optionmenu");
//    GtkWidget *chain_option_menu =
//       lookup_widget(window, "renumber_residue_range_chain_optionmenu");
   
   // renumber_residue_range_resno_1_entry
   // renumber_residue_range_resno_2_entry
   // renumber_residue_range_offset_entry

   // fill molecules option menu
   GtkSignalFunc callback_func = 
      GTK_SIGNAL_FUNC(graphics_info_t::renumber_residue_range_molecule_menu_item_select);
   g.fill_option_menu_with_coordinates_options(molecule_option_menu, callback_func);

}

void
graphics_info_t::fill_renumber_residue_range_internal(GtkWidget *w, int imol) {

   GtkWidget *chain_option_menu =
      lookup_widget(w, "renumber_residue_range_chain_optionmenu");
   GtkSignalFunc callback_func = 
      GTK_SIGNAL_FUNC(graphics_info_t::renumber_residue_range_chain_menu_item_select);
   std::string a = fill_chain_option_menu(chain_option_menu, imol, callback_func);
   if (a != "no-chain") {
      graphics_info_t::renumber_residue_range_chain = a;
   } 
}


void
graphics_info_t::renumber_residue_range_molecule_menu_item_select(GtkWidget *item,
								  GtkPositionType pos) {
   graphics_info_t::renumber_residue_range_molecule = pos;
   GtkWidget *window = lookup_widget(GTK_WIDGET(item),
				     "renumber_residue_range_dialog");
   graphics_info_t g;
   g.fill_renumber_residue_range_internal(window, pos);
}

void
graphics_info_t::renumber_residue_range_chain_menu_item_select(GtkWidget *item,
					        	       GtkPositionType pos) {

   char *data = NULL;
   data = (char *)pos;
   if (data)
      graphics_info_t::renumber_residue_range_chain = std::string(data);
}


// static
GtkWidget *
graphics_info_t::wrapped_create_diff_map_peaks_dialog(const std::vector<std::pair<clipper::Coord_orth, float> > &centres, float map_sigma) {

   GtkWidget *w = create_diff_map_peaks_dialog();
   difference_map_peaks_dialog = w; // save it for use with , and .
                                    // (globjects key press callback)
   GtkWidget *radio_button;
   GSList *diff_map_group = NULL;
   GtkWidget *button_vbox = lookup_widget(w, "diff_map_peaks_vbox");
   GtkWidget *frame;

   // for . and , synthetic clicking.
   gtk_object_set_user_data(GTK_OBJECT(w), GINT_TO_POINTER(centres.size()));

   // a cutn'paste jobby from fill_rotamer_selection_buttons().
   for (unsigned int i=0; i<centres.size(); i++) {
      std::string label = int_to_string(i+1);
      label += " ";
      label += float_to_string(centres[i].second);
      label += " (";
      label += float_to_string(centres[i].second/map_sigma);
      label += " sigma)";
      label += " at ";
      label += centres[i].first.format();
      radio_button = gtk_radio_button_new_with_label(diff_map_group,
						     label.c_str());
      std::string button_name = "difference_map_peaks_button_";
      button_name += int_to_string(i);
      
      diff_map_group = gtk_radio_button_group (GTK_RADIO_BUTTON (radio_button));
      gtk_widget_ref (radio_button);
      gtk_object_set_data_full (GTK_OBJECT (w),
				button_name.c_str(), radio_button,
				(GtkDestroyNotify) gtk_widget_unref);
      // int *iuser_data = new int;
      coot::diff_map_peak_helper_data *hd = new coot::diff_map_peak_helper_data;
      hd->ipeak = i;
      hd->pos = centres[i].first;
	 
      // *iuser_data = i;
      gtk_signal_connect (GTK_OBJECT (radio_button), "toggled",
			  GTK_SIGNAL_FUNC (on_diff_map_peak_button_selection_toggled),
			  hd);
             
       gtk_widget_show (radio_button);
       frame = gtk_frame_new(NULL);
       gtk_container_add(GTK_CONTAINER(frame), radio_button);
       gtk_box_pack_start (GTK_BOX (button_vbox),
			   frame, FALSE, FALSE, 0);
       gtk_container_set_border_width (GTK_CONTAINER (frame), 2);
       gtk_widget_show(frame);
   }

   // not used in the callback now that the button contains a pointer
   // to this info:
   diff_map_peaks->resize(0);
   for (unsigned int i=0; i<centres.size(); i++) 
      diff_map_peaks->push_back(centres[i].first);
   max_diff_map_peaks = centres.size();

   if (centres.size() > 0) {
      graphics_info_t g;
      coot::Cartesian c(centres[0].first.x(), centres[0].first.y(), centres[0].first.z());
      g.setRotationCentre(c);
      for(int ii=0; ii<n_molecules(); ii++) {
	 molecules[ii].update_map();
	 molecules[ii].update_symmetry();
      }
      graphics_draw();
   }
   return w;
}


// static 
void
graphics_info_t::on_diff_map_peak_button_selection_toggled (GtkButton       *button,
							    gpointer         user_data) {

   coot::diff_map_peak_helper_data *hd = (coot::diff_map_peak_helper_data *) user_data;
   // int i = hd->ipeak;
   
   graphics_info_t g;
   // std::cout << "button number " << i << " pressed\n";
   if (GTK_TOGGLE_BUTTON(button)->active) { 
      // std::cout << "button number " << i << " was active\n";
      coot::Cartesian c(hd->pos.x(), hd->pos.y(), hd->pos.z());
      g.setRotationCentre(c);
      for(int ii=0; ii<n_molecules(); ii++) {
	 molecules[ii].update_map();
	 molecules[ii].update_symmetry();
      }
      g.make_pointer_distance_objects();
      graphics_draw();
      std::string s = "Difference map peak number ";
      s += int_to_string(hd->ipeak);
      g.statusbar_text(s);
   } 
}

void
graphics_info_t::clear_diff_map_peaks() {

   diff_map_peaks->resize(0);
   max_diff_map_peaks = 0;

   // Also need to clear the user data on the buttons of the widget -
   // so I should pass the widget pointer then, shouldn't I?

}



// -------- keyboard rotamer control: ---------
// static
void
graphics_info_t::rotamer_dialog_neighbour_rotamer(int istep) {

   graphics_info_t g;
   if (g.rotamer_dialog) {
      // void *t  = (void *) (gtk_object_get_user_data(GTK_OBJECT(g.rotamer_dialog)));
      // std::cout << "user data: " << t << std::endl;
      int n_rotamers = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(g.rotamer_dialog)));
      // std::cout << "We find " << n_rotamers << " rotamers in the widget\n";
      GtkWidget *button;
      short int ifound_active_button = 0;
      int active_button_number = 0;
      int new_active_button_number;
      for (int i=0; i<n_rotamers; i++) {
	 std::string button_name = "rotamer_selection_button_rot_";
	 button_name += int_to_string(i);
	 button = lookup_widget(g.rotamer_dialog, button_name.c_str());
	 if (button) { 
	    if (GTK_TOGGLE_BUTTON(button)->active) {
	       ifound_active_button = 1;
	       active_button_number = i;
	       break;
	    }
	 } else {
	    std::cout << "ERROR:: rotamer button not found " << button_name << std::endl;
	 }
      }
      if (ifound_active_button) {
	 if (istep == 1) {
	    new_active_button_number = active_button_number + 1;
	    if (new_active_button_number == n_rotamers) {
	       new_active_button_number = 0;
	    }
	 } else {
	    new_active_button_number = active_button_number - 1;
	    if (new_active_button_number < 0) {
	       new_active_button_number = n_rotamers -1;
	    }
	 }
	 std::string button_name = "rotamer_selection_button_rot_";
	 button_name += int_to_string(new_active_button_number);
	 GtkWidget *new_button = lookup_widget(g.rotamer_dialog, button_name.c_str());
	 gtk_signal_emit_by_name(GTK_OBJECT(new_button), "clicked");
      
      } else {
	 std::cout << "ERROR:: not active rotamer button found " << std::endl;
      }
   }
}

void
graphics_info_t::rotamer_dialog_next_rotamer() {

   graphics_info_t::rotamer_dialog_neighbour_rotamer(+1);
}



// static
void
graphics_info_t::rotamer_dialog_previous_rotamer() {

   graphics_info_t::rotamer_dialog_neighbour_rotamer(-1);
}


// -------- keyboard difference map peak control: ------------
// static 
void graphics_info_t::difference_map_peaks_next_peak() {
   graphics_info_t::difference_map_peaks_neighbour_peak(1);
}

// static
void graphics_info_t::difference_map_peaks_previous_peak() {
   graphics_info_t::difference_map_peaks_neighbour_peak(-1);
}

// static
void graphics_info_t::difference_map_peaks_neighbour_peak(int istep) { // could be private

   graphics_info_t g;
   if (g.difference_map_peaks_dialog) {
      int n_peaks = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(g.difference_map_peaks_dialog)));
      GtkWidget *button;
      short int ifound_active_button = 0;
      int active_button_number = -99;     // set later
      int new_active_button_number = -99; // set later
      for (int i=0; i<n_peaks; i++) {
	 std::string button_name = "difference_map_peaks_button_";
	 button_name +=  int_to_string(i);
	 button = lookup_widget(g.difference_map_peaks_dialog, button_name.c_str());
	 if (button) {
	    if (GTK_TOGGLE_BUTTON(button)->active) {
	       ifound_active_button = 1;
	       active_button_number = i;
	    }
	 } else {
	    std::cout << "DEBUG:: Failed to find button " << button_name << "\n";
	 } 
      }
      if (ifound_active_button) {
	 if (istep == 1) {
	    new_active_button_number = active_button_number +1;
	    if (new_active_button_number == n_peaks)
	       new_active_button_number = 0;
	 } else {
	    new_active_button_number = active_button_number - 1;
	    if (new_active_button_number < 0)
	       new_active_button_number = n_peaks -1;
	 }
      }
      std::string button_name = "difference_map_peaks_button_";
      button_name += int_to_string(new_active_button_number);
      GtkWidget *new_button = lookup_widget(g.difference_map_peaks_dialog,
					    button_name.c_str());
      gtk_signal_emit_by_name(GTK_OBJECT(new_button), "clicked");
      
   } else {
      std::cout << "ERROR:: difference_map_peaks_neighbour_peak called in error\n";
   }
}

// static
void
graphics_info_t::checked_waters_next_baddie(int dir) {

   graphics_info_t g;
   GtkWidget *dialog = g.checked_waters_baddies_dialog;
   if (dialog) {
      int n_baddies = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(dialog)));
      GtkWidget *button;
      bool ifound_active_button = 0;
      int active_button_number = -99; // set later
      int new_active_button_number = -99; // set later
      
      for (int i=0; i<n_baddies; i++) {
	 std::string button_name = "checked_waters_baddie_button_";
	 button_name += int_to_string(i);
	 button = lookup_widget(dialog, button_name.c_str());
	 if (button) {
	    if (GTK_TOGGLE_BUTTON(button)->active) {
	       ifound_active_button = 1;
	       active_button_number = i;
	    }
	 } else {
	    std::cout << "failed to find button " << button_name
		      << std::endl;
	 }
      }
      if (ifound_active_button) {
	 if (dir == 1) {
	    new_active_button_number = active_button_number + 1;
	    if (new_active_button_number == n_baddies) {
	       new_active_button_number = 0;
	    }
	 } else {
	    new_active_button_number = active_button_number - 1;
	    if (new_active_button_number < 0)
	       new_active_button_number = n_baddies - 1;
	 }
	 std::string active_button_name = "checked_waters_baddie_button_";
	 active_button_name += int_to_string(new_active_button_number);
	 GtkWidget *new_active_button =
	    lookup_widget(dialog, active_button_name.c_str());
	 gtk_signal_emit_by_name(GTK_OBJECT(new_active_button), "clicked");
      } else {
	 std::cout << "active button not found" << std::endl;
      }
   }
}



#ifdef USE_GUILE
// static
SCM
graphics_info_t::safe_scheme_command(const std::string &scheme_command) { 

   // FIXME!
   SCM handler = scm_c_eval_string ("(lambda (key . args) (display (list \"(safe_scheme_command) Error in proc: key: \" key \" args: \" args)) (newline))"); 

   // I am undecided if I want this or not:
   // std::cout << "safe running :" << scheme_command << ":" << std::endl; 
   std::string thunk("(lambda() "); 
   thunk += scheme_command; 
   thunk += " )";
   SCM scm_thunk = SCM_BOOL_F;

   // try/catch does not make flow control come back here when bad thunk.
   //    std::cout << "..... making thunk... " << std::endl;
   scm_thunk = scm_c_eval_string(thunk.c_str());
   //    std::cout << "..... making thunk done. " << std::endl;

   //    std::cout << "..... scm_catch... " << std::endl;
   SCM v = scm_catch(SCM_BOOL_T, scm_thunk, handler);
   //    std::cout << "..... scm_catch done " << std::endl;

   SCM dest = SCM_BOOL_F;
   SCM mess = scm_makfrom0str("scm_catch returns: ~s\n");
   SCM sf = scm_simple_format(dest, mess, scm_list_1(v));
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
   std::string bad_str = scm_to_locale_string(sf);
#else   
   std::string bad_str = SCM_STRING_CHARS(sf);
#endif    
   // std::cout << bad_str << std::endl;

//   int is_int_p = scm_integer_p(v);
//   if (is_int_p) { 
//      std::cout << "returned value was int: " <<  std::endl;
//   } 

  // std::cout << "INFO:: finished scheme command " << std::endl;
   return v;
}
#endif // USE_GUILE

#ifdef USE_GUILE
// static
SCM
graphics_info_t::process_socket_string_waiting() {

   SCM r = SCM_BOOL_F;   // was unitiailized
   if (graphics_info_t::have_socket_string_waiting_flag) {

      graphics_info_t::have_socket_string_waiting_flag = 0; // draw() looks here
      std::string ss = graphics_info_t::socket_string_waiting;

      // really the right way?  Perhaps we should just stick to scheme
      // internals?
      std::vector<std::string> v;
      v.push_back("eval-socket-string");
      v.push_back(coot::util::single_quote(ss));

      graphics_info_t g;
      std::string s = g.state_command(v, coot::STATE_SCM);
      r = safe_scheme_command(s);
      
   }
   return r;
} 
#endif

// static 
gboolean
graphics_info_t::process_socket_string_waiting_bool(gpointer user_data) {
 
#ifdef USE_GUILE
   if (graphics_info_t::have_socket_string_waiting_flag) {
      graphics_info_t::have_socket_string_waiting_flag = 0; // draw() looks here
      std::string ss = graphics_info_t::socket_string_waiting;

      // try internal evaluation:
      if (1) {
	 SCM ss_scm = scm_makfrom0str(ss.c_str());

	 std::cout << "DEBUG: evaluting :" << ss << ":" << std::endl;
	 // SCM r = safe_scheme_command(ss);
	 // SCM r = scm_eval_string(ss_scm);
	 scm_eval_string(ss_scm);
	 // should store r.
	 std::cout << "DEBUG: done evaluating" << std::endl;
      }

      // really the right way?  Perhaps we should just stick to scheme
      // internals?
      if (0) { 
	 std::vector<std::string> v;
	 v.push_back("eval-socket-string");
	 v.push_back(coot::util::single_quote(ss));

	 graphics_info_t g;
	 std::string s = g.state_command(v, coot::STATE_SCM);
	 safe_scheme_command(s);
      }
   }
   std::cout << " =============== unsetting mutex lock =========" << std::endl;
   graphics_info_t::socket_string_waiting_mutex_lock = 0; // we're done.  release lock.
   return FALSE; // don't call this function again, idly.

#else // USE_GUILE

#ifdef USE_PYTHON

   // Bernhard to fill this part.

#endif // USE_PYTHON   
   
   return FALSE; 
#endif
}



// static 
std::pair<short int, float>
graphics_info_t::float_from_entry(GtkWidget *entry) {

   std::pair<short int, float> p(0,0);
   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      float f = atof(txt);
      p.second = f;
      p.first = 1;
   }
   return p;
}

std::pair<short int, int>
graphics_info_t::int_from_entry(GtkWidget *entry) {

   std::pair<short int, int> p(0,0);
   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      int i = atoi(txt);
      p.second = i;
      p.first = 1;
   }
   return p;
}

// static
std::string
graphics_info_t::backslash_filename(const std::string &s) { // needed for windows?

   std::string r = s;

   for (unsigned int i=0; i<s.size(); i++) {
      if (s[i] == '/')
	 r[i] = '\\';
   }
   return r;
}



// symmetry control dialog:
GtkWidget *graphics_info_t::wrapped_create_symmetry_controller_dialog() const {

   GtkWidget *w = create_symmetry_controller_dialog();
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_model())
 	 molecules[imol].fill_symmetry_control_frame(w);
   }
   return w;
} 

GtkWidget *
graphics_info_t::wrapped_create_lsq_plane_dialog() {

   GtkWidget *w = create_lsq_plane_dialog();
   pick_cursor_maybe();
   lsq_plane_dialog = w;
   GtkWindow *main_window = GTK_WINDOW(lookup_widget(glarea, "window1"));
   gtk_window_set_transient_for(GTK_WINDOW(w), main_window);
   
   return w;
} 


void
graphics_info_t::render_lsq_plane_atoms() {  // put a blob at atoms in lsq_plane_atom_positions

   if (lsq_plane_atom_positions->size() > 0) { 
      glColor3f(0.6, 0.6, 0.9);
      glPointSize(8.0);
      glBegin(GL_POINTS); 
      for (unsigned int i=0; i<lsq_plane_atom_positions->size(); i++) {
	 glVertex3f((*lsq_plane_atom_positions)[i].x(),
		    (*lsq_plane_atom_positions)[i].y(),
		    (*lsq_plane_atom_positions)[i].z());
      }
      glEnd();
   }
}

int
graphics_info_t::measure_lsq_plane_deviant_atom(int imol, int atom_index) {

   int r = 0;
   if (molecules[imol].has_model()) {
      CAtom *at = molecules[imol].atom_sel.atom_selection[atom_index];
      clipper::Coord_orth p(at->x, at->y, at->z);

      if (lsq_plane_atom_positions->size() > 2) {

	 graphics_draw();
	 std::pair<float,float> d_pair =
	    coot::lsq_plane_deviation(*lsq_plane_atom_positions, p);
	 float d = d_pair.first;
	 
	 std::string s("Atom ");
	 s += at->name;
	 std::string a(at->altLoc);
	 if (a != "") {
	    s += ",";
	    s += a;
	 }
	 s += " ";
	 s += int_to_string(at->GetSeqNum());
	 s += at->GetChainID();
	 s += " is ";
	 s += float_to_string(d);
	 s += "A from the least squares plane";
	 statusbar_text(s);
      } else {
	 std::string s("Not enough atoms to find plane");
	 std::cout << s << "\n";
	 statusbar_text(s);
      }
   }
   return r;
} 


int
graphics_info_t::add_lsq_plane_atom(int imol, int atom_index) {

   if (molecules[imol].has_model()) {
      CAtom *at = molecules[imol].atom_sel.atom_selection[atom_index];
      clipper::Coord_orth p(at->x, at->y, at->z);
      std::string s("Added plane atom ");
      s += at->name;
      s += " ";
      s += int_to_string(at->GetSeqNum());
      s += at->GetChainID();
      std::cout << s << std::endl;
      statusbar_text(s);
      lsq_plane_atom_positions->push_back(p);
      graphics_draw();
   }
   return 0; 
}

int
graphics_info_t::remove_last_lsq_plane_atom() {

   if (lsq_plane_atom_positions->size() > 1) {
      lsq_plane_atom_positions->resize(lsq_plane_atom_positions->size()-1);
      graphics_draw();
   }
   return 0;
}

coot::view_info_t
coot::view_info_t::interpolate(const coot::view_info_t &view1,
			       const coot::view_info_t &view2,
			       int n_steps) {
   coot::view_info_t view;
   graphics_info_t g;

//    std::cout << "start quat interpolation: zooms: " << view1.zoom << " " << view2.zoom
// 	     << " and centres: "
// 	     << view1.rotation_centre << " to " << view2.rotation_centre << std::endl;

//     std::cout << "quaternion interpolation using " << n_steps << " steps"
// 	      << std::endl;

   float total_zoom_by = view2.zoom/view1.zoom;
   float frac_zoom_by =  1;
   int smooth_scroll_state = graphics_info_t::smooth_scroll;
   graphics_info_t::smooth_scroll = 0;
   if (n_steps > 0)
      frac_zoom_by = total_zoom_by/float(n_steps); 

   // non-slerping
   // 
   if (0) { 
      for (int i=0; i<=n_steps; i++) {
	 float frac = float(i)/float(n_steps);
	 coot::Cartesian rct =
	    view1.rotation_centre + (view2.rotation_centre - view1.rotation_centre).by_scalar(frac);
	 for (int iq=0; iq<4; iq++)
	    g.quat[iq] = view1.quat[iq] + frac*(view2.quat[iq]-view1.quat[iq]);
	 g.zoom = view1.zoom + frac*(view2.zoom-view1.zoom);
	 g.setRotationCentre(rct);
	 graphics_info_t::graphics_draw();
      }
   }

   if (1) {
      float omega = acos(coot::view_info_t::dot_product(view1, view2)/(view1.quat_length()*view2.quat_length()));
      if (omega != 0.0) { 
	 // slerping
	 //
	 if (n_steps < 1)
	    n_steps = 1;
	 float t_step = float(1.0)/float(n_steps);
	 for (float t=0; t<=1; t+=t_step) {
	    float one_over_sin_omega = 1/sin(omega);
	    float frac1 = sin((1-t)*omega) * one_over_sin_omega;
	    float frac2 = sin(t*omega) * one_over_sin_omega;
	    for (int iq=0; iq<4; iq++)
	       g.quat[iq] = frac1*view1.quat[iq] + frac2*view2.quat[iq];
	    coot::Cartesian rct =
	       view1.rotation_centre + (view2.rotation_centre - view1.rotation_centre).by_scalar(t);
	    g.setRotationCentre(rct);
	    g.zoom = view1.zoom + t*(view2.zoom-view1.zoom);
	    graphics_info_t::graphics_draw();
	 }
      } else {
	 // non slerping
	 for (int i=0; i<=n_steps; i++) {
	    float frac = float(i)/float(n_steps);
	    coot::Cartesian rct =
	       view1.rotation_centre + (view2.rotation_centre - view1.rotation_centre).by_scalar(frac);
	    for (int iq=0; iq<4; iq++)
	       g.quat[iq] = view1.quat[iq] + frac*(view2.quat[iq]-view1.quat[iq]);
	    g.setRotationCentre(rct);
	    g.zoom = view1.zoom + frac*(view2.zoom-view1.zoom);
	    graphics_info_t::graphics_draw();
	 }
      }
   }

   graphics_info_t::smooth_scroll = smooth_scroll_state;
   return view;
}

float
coot::view_info_t::quat_length() const {

   float d = 0.0;
   for (int i=0; i<4; i++) {
      d += quat[i]*quat[i];
   }
   return sqrt(d);
}


// static
float 
coot::view_info_t::dot_product(const coot::view_info_t &view1,
			       const coot::view_info_t &view2) {

   float d = 0.0;
   for (int i=0; i<4; i++) {
      d += view1.quat[i]*view2.quat[i];
   }
   return d;
}

std::ofstream&
coot::operator<<(std::ofstream &f, coot::view_info_t &view) {

   // position quaternion zoom view-name
   //

#ifdef USE_GUILE
   if (! view.is_simple_spin_view_flag) { 
      f << "(add-view ";
      f << "(list ";
      f << "   ";
      f << view.rotation_centre.x();
      f << " ";
      f << view.rotation_centre.y();
      f << " ";
      f << view.rotation_centre.z();
      f << ")\n";

      f << "   (list ";
      f << view.quat[0]; 
      f << " ";
      f << view.quat[1]; 
      f << " ";
      f << view.quat[2]; 
      f << " ";
      f << view.quat[3];
      f << ")\n";
      
      f << "   ";
      f << view.zoom; 
      f << "\n";

      f << "   ";
      f << coot::util::single_quote(view.view_name);
   
      f << ")\n";
   } else {
      f << "(add-spin-view ";
      f << coot::util::single_quote(view.view_name);
      f << " ";
      f << view.n_spin_steps;
      f << " ";
      f << view.degrees_per_step * view.n_spin_steps;
      f << ")\n";
   }
#else
#ifdef USE_PYTHON
   if (! view.is_simple_spin_view_flag) { 
      f << "add_view(";
      f << "[";
      f << view.rotation_centre.x();
      f << ", ";
      f << view.rotation_centre.y();
      f << ", ";
      f << view.rotation_centre.z();
      f << "],\n";

      f << "   [";
      f << view.quat[0]; 
      f << ", ";
      f << view.quat[1]; 
      f << ", ";
      f << view.quat[2]; 
      f << ", ";
      f << view.quat[3];
      f << "],\n";
      
      f << "   ";
      f << view.zoom; 
      f << ",\n";

      f << "   ";
      f << coot::util::single_quote(view.view_name);
   
      f << ")\n";
   } else {
      f << "add_spin_view(";
      f << coot::util::single_quote(view.view_name);
      f << ", ";
      f << view.n_spin_steps;
      f << ", ";
      f << view.degrees_per_step * view.n_spin_steps;
      f << ")\n";
   }
#endif // USE_PYTHON
#endif // USE_GUILE
   return f;

}

   
bool
coot::view_info_t::matches_view (const coot::view_info_t &view) const { 
   float frac = 0.01;
   bool matches = 0;
   if (zoom < view.zoom*(1+frac)) { 
      if (zoom > view.zoom*(1-frac)) { 
	 if (rotation_centre.x() < view.rotation_centre.x()*(1+frac)) { 
	    if (rotation_centre.x() > view.rotation_centre.x()*(1-frac)) { 
	       if (rotation_centre.y() < view.rotation_centre.y()*(1+frac)) { 
		  if (rotation_centre.y() > view.rotation_centre.y()*(1-frac)) { 
		     if (rotation_centre.z() < view.rotation_centre.z()*(1+frac)) { 
			if (rotation_centre.z() > view.rotation_centre.z()*(1-frac)) { 
			   if (quat[0] < view.quat[0]*(1+frac)) { 
			      if (quat[0] > view.quat[0]*(1-frac) ){ 
				 if (quat[1] < view.quat[1]*(1+frac)) { 
				    if (quat[1] > view.quat[1]*(1-frac) ){ 
				       if (quat[2] < view.quat[2]*(1+frac)) { 
					  if (quat[2] > view.quat[2]*(1-frac) ){ 
					     if (quat[3] < view.quat[3]*(1+frac)) { 
						if (quat[3] > view.quat[3]*(1-frac) ){
						   matches = 1;
						}
					     }
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return matches;
}
