/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
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

// $Id: c-interface.cc 1458 2007-01-26 20:20:18Z emsley $
// $LastChangedDate: 2007-01-26 20:20:18 +0000 (Fri, 26 Jan 2007) $
// $Rev: 1458 $
 

#include <stdlib.h>
#include <iostream>
#include <fstream>

#if !defined(WINDOWS_MINGW) && !defined(_MSC_VER)
#include <glob.h> // for globbing.
#endif

#ifdef USE_GUILE
#include <guile/gh.h>
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON


#if defined(_MSC_VER)
#define usleep(x) Sleep(x/1000)
#include <windows.h>
#endif

// Here we used to define GTK_ENABLE_BROKEN if defined(WINDOWS_MINGW)
// Now we don't want to enable broken stuff.  That is not the way.

#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#define sleep Sleep
#define snprintf _snprintf
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

#include "clipper/ccp4/ccp4_map_io.h"
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.

#include <vector>
#include <string>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"
#include "coot-utils.hh"
#include "coot-map-utils.hh"
#include "coot-database.hh"

// #include "xmap-interface.h"
#include "graphics-info.h"

#include "atom-utils.h" // asc_to_graphics

#include "BuildCas.h"

#include "trackball.h" // adding exportable rotate interface


#include "c-interface.h"
#include "cc-interface.hh"

#include "positioned-widgets.h"

// moving column_label selection to c-interface from mtz bits.
#include "cmtz-interface.hh"
// #include "mtz-bits.h" stuff from here moved to cmtz-interface

char *coot_version() {

   std::string version_string = VERSION;
#ifdef USE_GUILE
   version_string += " [with guile ";
   version_string += coot::util::int_to_string(SCM_MAJOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(SCM_MINOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(SCM_MICRO_VERSION);
   version_string += " embedded]";
#endif    
#ifdef USE_PYTHON
   version_string += " [with python ";
   version_string += coot::util::int_to_string(PY_MAJOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(PY_MINOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(PY_MICRO_VERSION);
   version_string += " embedded]";
#endif    
   int len = version_string.length();
   
   char *r = new char[len+1];
   strncpy(r, version_string.c_str(), len+1);
   return r;
}

// Return 0 if not a valid name ( -> #f in scheme)
// e.g. /a/b/c.pdb "d/e/f.mtz FWT PHWT"
// 
const char *molecule_name(int imol) {

   const char *r = NULL;
   if (is_valid_map_molecule(imol)) {
      r = graphics_info_t::molecules[imol].name_.c_str();
      return r;
   }
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].name_.c_str();
   }
   std::string cmd = "molecule-name";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   
   return r;
} 

//  Display characteristics:
//
// 
void set_esoteric_depth_cue(int istate) {

   graphics_info_t::esoteric_depth_cue_flag = istate;
   std::string cmd = "set-esoteric-depth-cue";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

int  esoteric_depth_cue_state() {
   add_to_history_simple("esoteric-depth-cue-state");
   return graphics_info_t::esoteric_depth_cue_flag;
}


/*! \brief shall we start up the Gtk and the graphics window? 

   if passed the command line argument --no-graphics, coot will not
   start up gtk itself.

   An interface function for Ralf.
*/
short int use_graphics_interface_state() {

   add_to_history_simple("use-graphics-interface-state");
   return graphics_info_t::use_graphics_interface_flag; 

}

/*! \brief start Gtk (and graphics) 

   This function is useful if it was not started already (which can be
   achieved by using the command line argument --no-graphics).

   An interface for Ralf */
void start_graphics_interface() {
   add_to_history_simple("start-graphics-interface");
   gtk_main(); 
} 




// Is this a repeat of something?  I don't know (doesn't look like it)
// 
void run_generic_script(const std::vector<std::string> &cmd_strings) {
   
   graphics_info_t g;

#ifdef USE_GUILE   
   std::string s = g.state_command(cmd_strings, coot::STATE_SCM);
   safe_scheme_command(s);
#endif

#ifdef USE_PYTHON
#ifndef USE_GUILE
   std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
   safe_python_command(s);
#endif    
#endif    

   std::string cmd = "run-generic-script";
   std::vector<coot::command_arg_t> args;
   for(unsigned int i=0; i<cmd_strings.size(); i++) 
      args.push_back(clipper::String(cmd_strings[i]));
   add_to_history_typed(cmd, args);
}

   

/*  ----------------------------------------------------------------------- */
/*                  Display lists                                           */
/*  ----------------------------------------------------------------------- */
void set_display_lists_for_maps(int istat) {

   graphics_info_t::display_lists_for_maps_flag = istat;
   for (int i=0; i<graphics_info_t::n_molecules; i++)
      if (graphics_info_t::molecules[i].has_map())
	 graphics_info_t::molecules[i].update_map();
   std::string cmd = "set-display-lists-for-maps";
   std::vector<coot::command_arg_t> args;
   args.push_back(istat);
   add_to_history_typed(cmd, args);
   graphics_draw();
}


//
// Return the molecule number of the molecule that we just filled.
// Return -1 if there was a failure.
// 
int handle_read_draw_molecule(const char *filename) {

   int r = graphics_info_t::recentre_on_read_pdb;
   return handle_read_draw_molecule_with_recentre(filename, r);
}


int handle_read_draw_molecule_with_recentre(const char *filename,
					   int recentre_on_read_pdb_flag) {

   int r = -1;
   //
   // cout << "handle_read_draw_molecule: handling " << filename << endl;
   
   graphics_info_t g;
   int imol = g.n_molecules;
   if (! filename)
      return -1;
   
   std::string f(filename);

   g.expand_molecule_space_maybe();

   // returns e.g. ".ins"

   int  istat = -1;
   std::string extention = coot::util::file_name_extension(filename);
   if (coot::util::extension_is_for_shelx_coords(extention)) {

      //       if (0) { // Don't do this path.  handle_read_draw_molecule() call
	    // get_atom_selection() which does this check and goes
	    // down the shelx path already.
      // But!  That's not the whole story read_shelx_ins_file calls
      // the read_shelx_ins_file() method of molecule_class_info_t,
      // which does different things to handle_read_draw_molecule,
      // specifically, it doesn't call get_atom_selection() and so the
      // Hydrogen names are not fixed.  Currently that is the case.
      // We'd like the fixing of hydrogen names in this
      // read_shelx_ins_file path too.
      
      r = read_shelx_ins_file(filename);

   } else { 
      // recentre and not a backup-restore
      // -1 is for failure strangely.
      istat = g.molecules[imol].handle_read_draw_molecule(f, recentre_on_read_pdb_flag, 0);

      if (istat > -1) {
	 std::cout << "Molecule " << g.n_molecules << " read successfully\n";
	 g.n_molecules++;

	 // if the go to atom widget exists, update its optionmenu to
	 // reflect the existance of this new molecule.

	 if (g.go_to_atom_window) {
	    g.set_go_to_atom_molecule(imol);
	    g.update_go_to_atom_window_on_new_mol();
	    // g.update_go_to_atom_window_on_changed_mol(imol);
	 } else {
	    // The Go To Atom window is not displayed.
	    g.set_go_to_atom_molecule(imol);
	 }
      
	 // now force a draw of the molecule
	 //
	 graphics_draw();

	 std::string s("Successfully read coordinates file ");
	 s += filename;
	 s += ".  Molecule number ";
	 s += coot::util::int_to_string(imol);
	 s += " created.";
	 g.statusbar_text(s);
	 r =  imol;
      } else {
	 std::string s("Failed to read coordinates file ");
	 s += filename;
	 g.statusbar_text(s);
	 r =  -1;
      }
   }

   std::string cmd = "handle-read-draw-molecule-with-recentre";
   std::vector<coot::command_arg_t> args;
   args.push_back(filename);
   args.push_back(recentre_on_read_pdb_flag);
   add_to_history_typed(cmd, args);
   return r; 
}


int read_pdb(const char *filename) {
   // history is done in the handle function
   return handle_read_draw_molecule(filename); 
} 

void set_draw_zero_occ_markers(int status) { 
   
   graphics_info_t g;
   g.draw_zero_occ_spots_flag = status;
   std::string cmd = "set-draw-zero-occ-markers";
   std::vector<coot::command_arg_t> args;
   args.push_back(status);
   add_to_history_typed(cmd, args);
   graphics_draw();
}



int first_coords_imol() {

   int imol = -1;
   for (int i=0; i<graphics_n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }
   add_to_history_simple("first-coords-imol");
   return imol;
}



void hardware_stereo_mode() {

   if (graphics_info_t::use_graphics_interface_flag) { 
      if (graphics_info_t::display_mode != coot::HARDWARE_STEREO_MODE) {
	 int previous_mode = graphics_info_t::display_mode;
	 graphics_info_t::display_mode = coot::HARDWARE_STEREO_MODE;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in hardware_stereo_mode!\n";
	 } else {

	    if ( (previous_mode == coot::SIDE_BY_SIDE_STEREO) ||
		 (previous_mode == coot::DTI_SIDE_BY_SIDE_STEREO)  || 
		 (previous_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ) {

	       if (graphics_info_t::glarea_2) {
		  gtk_widget_destroy(graphics_info_t::glarea_2);
		  graphics_info_t::glarea_2 = NULL;
	       }
	    }
	    
	    short int try_hardware_stereo_flag = 1;
	    GtkWidget *glarea = gl_extras(vbox, try_hardware_stereo_flag);
	    if (glarea) { 
	       std::cout << "INFO:: switch to hardware_stereo_mode succeeded\n";
	       if (graphics_info_t::i_fn_token) { 
		  toggle_idle_function(); // turn it off;
	       }
	       gtk_widget_destroy(graphics_info_t::glarea);
	       graphics_info_t::glarea = glarea;
	       gtk_widget_show(glarea);
	       graphics_draw();
	    } else {
	       std::cout << "WARNING:: switch to hardware_stereo_mode failed\n";
	       graphics_info_t::display_mode = previous_mode;
	    }
	 }
      } else {
	 std::cout << "Already in hardware stereo mode" << std::endl;
      }
   }
   add_to_history_simple("hardware-stereo-mode");

}

void mono_mode() {

   if (graphics_info_t::use_graphics_interface_flag) {
      
      if (graphics_info_t::display_mode != coot::MONO_MODE) { 
	 int previous_mode = graphics_info_t::display_mode;
	 graphics_info_t::display_mode = coot::MONO_MODE;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in mono mode!\n";
	 } else {
	    short int try_hardware_stereo_flag = 0;
	    GtkWidget *glarea = gl_extras(vbox, try_hardware_stereo_flag);
	    if (glarea) { 
	       std::cout << "INFO:: switch to mono_mode succeeded\n";
	       if (graphics_info_t::i_fn_token) { 
		  toggle_idle_function(); // turn it off;
	       }
	       gtk_widget_destroy(graphics_info_t::glarea);
	       if (graphics_info_t::glarea_2) { 
		  gtk_widget_destroy(graphics_info_t::glarea_2);
		  graphics_info_t::glarea_2 = NULL;
	       }
	       graphics_info_t::glarea = glarea;
	       gtk_widget_show(glarea);
	       graphics_draw();
	    } else {
	       graphics_info_t::display_mode = previous_mode;
	       std::cout << "WARNING:: switch to mono mode failed\n";
	    }
	 }
      } else {
	 // std::cout << "Already in mono mode" << std::endl; // we know.
      }
   }
   add_to_history_simple("mono-mode");
}

/*! \brief turn on side bye side stereo mode */
void side_by_side_stereo_mode(short int use_wall_eye_flag) {

   if (graphics_info_t::use_graphics_interface_flag) {

      
      // If it wasn't in side by side stereo mode, then we need to
      // generated 2 new glares by calling gl_extras().
      // 
      if (!((graphics_info_t::display_mode == coot::SIDE_BY_SIDE_STEREO) ||
	    (graphics_info_t::display_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ||
	    (graphics_info_t::display_mode == coot::DTI_SIDE_BY_SIDE_STEREO))) {
	 
	 if (use_wall_eye_flag == 1) {
	    graphics_info_t::in_wall_eyed_side_by_side_stereo_mode = 1;
	    graphics_info_t::display_mode = coot::SIDE_BY_SIDE_STEREO_WALL_EYE;
	 } else {
	    graphics_info_t::in_wall_eyed_side_by_side_stereo_mode = 0;
	    graphics_info_t::display_mode = coot::SIDE_BY_SIDE_STEREO;
	 }
	 int previous_mode = graphics_info_t::display_mode;
	 short int stereo_mode = coot::SIDE_BY_SIDE_STEREO;
	 if (use_wall_eye_flag)
	    stereo_mode = coot::SIDE_BY_SIDE_STEREO_WALL_EYE;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 GtkWidget *glarea = gl_extras(vbox, stereo_mode);
	 if (glarea) {
	    if (graphics_info_t::i_fn_token) { 
	       toggle_idle_function(); // turn it off;
	    }
	    gtk_widget_destroy(graphics_info_t::glarea);
	    graphics_info_t::glarea = glarea; // glarea_2 is stored by gl_extras()
	    gtk_widget_show(glarea);
	    gtk_widget_show(graphics_info_t::glarea_2);
	    graphics_draw();
	 } else {
	    std::cout << "WARNING:: switch to side by side mode failed!\n";
	 } 
      } else {

	 if (use_wall_eye_flag == 1) {
	    graphics_info_t::in_wall_eyed_side_by_side_stereo_mode = 1;
	    graphics_info_t::display_mode = coot::SIDE_BY_SIDE_STEREO_WALL_EYE;
	 } else {
	    graphics_info_t::in_wall_eyed_side_by_side_stereo_mode = 0;
	    graphics_info_t::display_mode = coot::SIDE_BY_SIDE_STEREO;
	 }
	 // were were already in some sort of side by side stereo mode:
	 graphics_draw();
      }
   }
   std::vector<coot::command_arg_t> args;
   args.push_back(use_wall_eye_flag);
   add_to_history_typed("side-by-side-stereo-mode", args);
}

/* DTI stereo mode - undocumented, secret interface for testing, currently */
// when it works, call it dti_side_by_side_stereo_mode()
void set_dti_stereo_mode(short int state) {

   if (graphics_info_t::use_graphics_interface_flag) {
      if (graphics_info_t::display_mode != coot::DTI_SIDE_BY_SIDE_STEREO) {
	 int previous_mode = graphics_info_t::display_mode;
	 short int stereo_mode = coot::DTI_SIDE_BY_SIDE_STEREO;
	 GtkWidget *vbox = lookup_widget(graphics_info_t::glarea, "vbox1");
	 GtkWidget *glarea = gl_extras(vbox, stereo_mode);
	 if (glarea) {
	    if (graphics_info_t::i_fn_token) { 
	       toggle_idle_function(); // turn it off;
	    }
	    gtk_widget_destroy(graphics_info_t::glarea);
	    graphics_info_t::glarea = glarea; // glarea_2 is stored by gl_extras()
	    gtk_widget_show(glarea);
	    gtk_widget_show(graphics_info_t::glarea_2);
	    graphics_draw();
	 } else {
	    std::cout << "WARNING:: switch to side by side mode failed!\n";
	 } 
      }
   }
   // add_to_history_simple("dti-side-by-side-stereo-mode");
} 


int stereo_mode_state() {
   add_to_history_simple("stereo-mode-state");
   return graphics_info_t::display_mode;
}

void set_hardware_stereo_angle_factor(float f) {
   graphics_info_t::hardware_stereo_angle_factor = f;
   std::string cmd = "set-hardware-stereo-angel-factor";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

float hardware_stereo_angle_factor_state() {
   add_to_history_simple("hardware-stereo-angle-factor-state");
   return graphics_info_t::hardware_stereo_angle_factor;
}


void graphics_draw() {
   graphics_info_t::graphics_draw();
}



void 
set_run_state_file_status(short int istat) {
   std::string cmd = "set-run-state-file-status";
   std::vector<coot::command_arg_t> args;
   args.push_back(istat);
   add_to_history_typed(cmd, args);
   
   graphics_info_t::run_state_file_status = istat;
}


void set_sticky_sort_by_date() {

   add_to_history_simple("set-sticky-sort-by-date");
   graphics_info_t g;
   g.sticky_sort_by_date = 1;

}

void add_filename_filter(GtkWidget *fileselection) { 

   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;

   GtkWidget *frame = gtk_frame_new("File-name filter:");
   GtkWidget *entry = gtk_entry_new();
   gtk_widget_set_usize (entry, 60, -2);

   gtk_widget_ref(entry);
   gtk_widget_show(entry);
   gtk_widget_ref(frame);

   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), entry);
   // I want the entry to be not-expandable:  How do I do that?
   // gtk_box_pack_start (GTK_BOX (hbox), entry, FALSE, TRUE, 0);

   gtk_signal_connect (GTK_OBJECT (entry), "key_press_event",
		       GTK_SIGNAL_FUNC (on_filename_filter_key_press_event),
		       NULL);
   gtk_widget_show(frame);

}

// where data type:
// 0 coords
// 1 mtz etc
// 2 maps
// 
GtkWidget *add_filename_filter_button(GtkWidget *fileselection, 
				      short int data_type) { 
   
   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *button = gtk_toggle_button_new_with_label("Filter");
   GtkWidget *frame = gtk_frame_new("File-name filter:");
   int d = data_type;
   
   gtk_widget_ref(button);
   gtk_widget_show(button);
   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), button);
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
   gtk_signal_connect (GTK_OBJECT (button), "toggled",
		       GTK_SIGNAL_FUNC (on_filename_filter_toggle_button_toggled_gtk1),
		       GINT_TO_POINTER(d));
#else   
   gtk_signal_connect (GTK_OBJECT (button), "toggled",
		       GTK_SIGNAL_FUNC (on_filename_filter_toggle_button_toggled),
		       GINT_TO_POINTER(d));
#endif   
   gtk_widget_show(frame);
   return button;
}

void add_is_difference_map_checkbutton(GtkWidget *fileselection) { 

   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *button = gtk_check_button_new_with_label("Is Difference Map");
   GtkWidget *frame = gtk_frame_new("Difference Map?");
   int imol = 0; // FIXME

   gtk_widget_ref(button);
   gtk_object_set_data_full(GTK_OBJECT(fileselection),
			    "map_fileselection_is_difference_map_button",
			    button,
			    (GtkDestroyNotify) gtk_widget_unref);
   gtk_widget_show(button);
   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), button);
   gtk_signal_connect (GTK_OBJECT (button), "toggled",
		       GTK_SIGNAL_FUNC (on_read_map_difference_map_toggle_button_toggled),
		       GINT_TO_POINTER(imol));
   gtk_widget_show(frame);

}

void set_filter_fileselection_filenames(int istate) {
   std::string cmd = "set-filter-fileselection-filenames";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);

   graphics_info_t::filter_fileselection_filenames_flag = istate;

}

int filter_fileselection_filenames_state() {
   add_to_history_simple("filter_fileselection_filenames_state");
   return graphics_info_t::filter_fileselection_filenames_flag;
}

void set_unit_cell_colour(float red, float green, float blue) {

   coot::colour_holder ch(red, green, blue);
   graphics_info_t::cell_colour = ch;
   graphics_draw();

   std::string cmd = "set-unit-cell-colour";
   std::vector<coot::command_arg_t> args;
   args.push_back(red);
   args.push_back(green);
   args.push_back(blue);
   add_to_history_typed(cmd, args);
} 


void swap_map_colours(int imol1, int imol2) {

   if (is_valid_map_molecule(imol1)) {
      if (is_valid_map_molecule(imol2)) {
	 graphics_info_t g;
	 std::vector<float> map_1_colours = g.molecules[imol1].map_colours();
	 std::vector<float> map_2_colours = g.molecules[imol2].map_colours();
	 double *colours1 = new double[4];
	 colours1[0] = map_1_colours[0];
	 colours1[1] = map_1_colours[1];
	 colours1[2] = map_1_colours[2];
	 double *colours2 = new double[4];
	 colours2[0] = map_2_colours[0];
	 colours2[1] = map_2_colours[1];
	 colours2[2] = map_2_colours[2];
	 g.molecules[imol1].handle_map_colour_change(colours2, g.swap_difference_map_colours);
	 g.molecules[imol2].handle_map_colour_change(colours1, g.swap_difference_map_colours);
	 delete [] colours1;
	 delete [] colours2;
      }
   }
   std::string cmd = "swap-map-colours";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol1);
   args.push_back(imol2);
   add_to_history_typed(cmd, args);
}

void set_keep_map_colour_after_refmac(int istate) {
   std::string cmd = "set-keep-map-colour-after-refmac";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
   graphics_info_t::swap_pre_post_refmac_map_colours_flag = istate;
}

int keep_map_colour_after_refmac_state() {
   add_to_history_simple("keep_map_colour_after_refmac_state");
   return graphics_info_t::swap_pre_post_refmac_map_colours_flag;
} 

void
on_read_map_difference_map_toggle_button_toggled (GtkButton       *button,
						  gpointer         user_data)
{
   if (GTK_TOGGLE_BUTTON(button)->active) { 
      std::cout << "is a difference map...!\n";
   } 
}

#if (GTK_MAJOR_VERSION > 1) || defined (GTK_ENABLE_BROKEN)
#else 
void
on_filename_filter_toggle_button_toggled_gtk1(GtkButton       *button,
					     gpointer         user_data)
{
   int data_type = GPOINTER_TO_INT(user_data);

   // We need to add text to the string of the dictectory we are in
   // (pre_directory), so first we need to find pre_directory (as per
   // fileselection_sort_button_clicked()
   // 
   GtkWidget *sort_button = lookup_widget(GTK_WIDGET(button),
					  "fileselection_sort_button");
   if (sort_button) { 
      // std::cout << "Hooray! we found the sort button!\n";
      // usually, we do.
   } else { 
      std::cout << "Boo we failed to find the sort button!\n";
   } 
   std::string pre_directory = pre_directory_file_selection(sort_button);
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button);
   
   if (fileselection) { 
      if (GTK_TOGGLE_BUTTON(button)->active) { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Unfilter");
	 
	 // so now we have pre_directory
	 // 
	 // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
	 std::vector<std::string> v = filtered_by_glob(pre_directory, data_type);
	 // std::cout << "DEBUG:: filtering by glob using data type: " << data_type
	 // << " returns" << std::endl;
	 // for (unsigned int iv=0; iv< v.size(); iv++)
	 // std::cout << iv << " " << v[iv] << std::endl;

	 filelist_into_fileselection_clist(fileselection, v);

      } else { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Filter");
	 gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
					 (pre_directory + "/").c_str());
      }
   } else {
      std::cout << "no fileselection found from sort button\n";
   }
}
#endif // GTK_MAJOR_VERSION > 1 or BROKEN

void 
filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v) {

   GtkCList  *file_list = GTK_CLIST(GTK_FILE_SELECTION(fileselection)->file_list);
   gtk_clist_clear(file_list);
   std::string::size_type islash;
   std::string t;
   for (unsigned int i=0; i<v.size(); i++) {
      islash = v[i].find_last_of("/");
      if (islash == string::npos) { 
	 // no slash found:
	 t = v[i];
      } else {
	 t = v[i].substr(islash + 1);
      }
      char *text = new char[t.length()+1];
      strncpy(text, t.c_str(), t.length()+1);
      gtk_clist_append(file_list, &text);
   }
}

std::vector<std::string> filtered_by_glob(const std::string &pre_directory, 
					  int data_type) { 

   std::vector<std::string> v; // returned object
   std::vector<std::string> globs;

#if !defined(WINDOWS_MINGW) && !defined(_MSC_VER)

   // std::map<std::string, int, std::less<std::string> >  files;

   if (data_type == 0) 
      globs = *graphics_info_t::coordinates_glob_extensions;
   if (data_type == 1) 
      globs = *graphics_info_t::data_glob_extensions;
   if (data_type == 2) 
      globs = *graphics_info_t::map_glob_extensions;
   if (data_type == 3) 
      globs = *graphics_info_t::dictionary_glob_extensions;

   for (unsigned int i=0; i<globs.size(); i++) { 

      std::string file_name_glob = pre_directory;
      file_name_glob += "/";

      file_name_glob += "*";
      file_name_glob += globs[i];
      glob_t myglob;
      int flags = 0;
      glob(file_name_glob.c_str(), flags, 0, &myglob);
      size_t count;

      char **p;
      for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) {
	 std::string f(*p);
	 if (! string_member(f, v))
	    v.push_back(f);
      }
      globfree(&myglob);
   }

   // now we need to sort v;
   std::sort(v.begin(), v.end(), compare_strings);

#endif // WINDOWS_MINGW

   return v;
}

bool
compare_strings(const std::string &a, const std::string &b) { 
   return a < b ? 1 : 0;
} 


// Return 1 if search appears in list, 0 if not)
// 
short int 
string_member(const std::string &search, const std::vector<std::string> &list) { 
   
   short int v = 0;
   for (unsigned int i=0; i<list.size(); i++) { 
      if (list[i] == search) { 
	 v = 1;
	 break;
      }
   }
   return v;
} 

#include <gdk/gdkkeysyms.h> // for keyboarding.

gboolean
on_filename_filter_key_press_event (GtkWidget       *widget,
				    GdkEventKey     *event,
				    gpointer         user_data)
{
   //    if (event->keyval == GDK_Return || event->keyval == GDK_Tab)
   //    { // Tab is not good.  It takes you to the next widget too,
   //    which is not want we want.  It's confusing, so let's just use
   //    return.
   if (event->keyval == GDK_Return) { 
      // std::cout << "Return pressed!\n";
#ifdef COOT_USE_GTK2_INTERFACE
      handle_filename_filter_gtk2(widget);
#else       

      handle_filename_filter_gtk1(widget);
#endif       
   } 
   return FALSE;
}


void
handle_filename_filter_gtk1(GtkWidget *entry_widget) {
   
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   // nothing
#else 

   std::cout << "running handle_filename_filter\n";

   const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry_widget));

   // We need to add text to the string of the dictectory we are in
   // (pre_directory), so first we need to find pre_directory (as per
   // fileselection_sort_button_clicked()
   // 
   GtkWidget *sort_button = lookup_widget(entry_widget, "fileselection_sort_button");
   if (sort_button) { 
      // std::cout << "Hooray! we found the sort button!\n";
      // usually, we do.
   } else { 
      std::cout << "Boo we failed to find the sort button!\n";
   } 
   std::string pre_directory = pre_directory_file_selection(sort_button);


   // so now we have pre_directory
   // 
   // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button);
   if (fileselection) { 
      GtkCList  *file_list = GTK_CLIST(GTK_FILE_SELECTION(fileselection)->file_list);
      std::string file_name_glob;
      file_name_glob = pre_directory;
      file_name_glob += "/";
      file_name_glob += text;

      glob_t myglob;
      int flags = 0;
      glob(file_name_glob.c_str(), flags, 0, &myglob);
      size_t count;


      char **p;
      std::vector<std::string> v;
      for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) { 
	 v.push_back(std::string(*p));
      }
      globfree(&myglob);

      if (v.size() > 0) { 
	 gtk_clist_clear(file_list);
	 std::string::size_type islash;
	 std::string t;
	 for (unsigned int i=0; i<v.size(); i++) {
	    islash = v[i].find_last_of("/");
	    if (islash == string::npos) { 
	       // no slash found:
	       t = v[i];
	    } else {
	       t = v[i].substr(islash + 1);
	    }
	    char *text = new char[t.length()+1];
	    strncpy(text, t.c_str(), t.length()+1);
	    gtk_clist_append(file_list, &text);
	 }
      }
   } else { 
      std::cout << "ERROR:: couldn't find fileselection\n";
   } 

#endif // WINDOWS_MINGW
}

std::string pre_directory_file_selection(GtkWidget *sort_button) { 

   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(gtk_object_get_user_data(GTK_OBJECT(sort_button)));

   // The menu item is a container than contains a label.
   // How do we get to the label given the container?
   // Strangely enough we use the history_pulldown.

   GList *dlist = gtk_container_children(GTK_CONTAINER(history_pulldown));
   GList *free_list = dlist;
   std::string pre_directory("");
   
   while (dlist) {
      // GtkWidget *list_item;
      // list_item = (GtkWidget *) (dlist->data);
      gchar *t = GTK_LABEL(dlist->data)->label;
      if (t != NULL) {
	 pre_directory = t; 
      } else {
	 std::cout << "WARNING:: null label t " << std::endl;
      } 
      dlist = dlist->next;
   }
   g_list_free(free_list); 

   return pre_directory;
}




void push_the_buttons_on_fileselection(GtkWidget *filter_button, 
				       GtkWidget *sort_button,
				       GtkWidget *fileselection) {

  if (filter_fileselection_filenames_state()) { 
    gtk_signal_emit_by_name(GTK_OBJECT(filter_button), "clicked");
  }
  if (graphics_info_t::sticky_sort_by_date) {
     GtkWidget *file_list = GTK_FILE_SELECTION(fileselection)->file_list;
     std::cout << "INFO:: Sorting files by date\n";
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
     fileselection_sort_button_clicked_gtk1(sort_button, (GtkCList *) file_list);
#else     
     fileselection_sort_button_clicked(sort_button, file_list);
#endif     
  }
}


// --------------------------------------------------------------------
// Ctrl for rotate or pick:
// --------------------------------------------------------------------
// 
// Coot mailing list discussion: users want Ctrl for Rotation or Ctrl
// for picking, so that accidental picking when rotation is meant is
// avoided.
// 
void set_control_key_for_rotate(int state) {
   graphics_info_t::control_key_for_rotate_flag = state;
}

int control_key_for_rotate_state() {
   return graphics_info_t::control_key_for_rotate_flag;
}




/*  ------------------------------------------------------------------------ */
/*                         Model/Fit/Refine Functions:                       */
/*  ------------------------------------------------------------------------ */
void post_model_fit_refine_dialog() {

   GtkWidget *widget = wrapped_create_model_fit_refine_dialog();
   if (graphics_info_t::use_graphics_interface_flag) { 
      gtk_widget_show(widget);
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("post-model-fit-refine-dialog");
   add_to_history(command_strings);
}

void post_other_modelling_tools_dialog() {

   GtkWidget *widget = wrapped_create_model_fit_refine_dialog();
   if (graphics_info_t::use_graphics_interface_flag) { 
      gtk_widget_show(widget);
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("post-other-modelling-tools-dialog");
   add_to_history(command_strings);

} 



void show_select_map_dialog() {
   graphics_info_t g;
   g.show_select_map_dialog();
   add_to_history_simple("show-select-map-dialog");
} 

// 
int make_and_draw_map(const char* mtz_file_name,
		      const char *f_col, const char *phi_col,
		      const char *weight_col, int use_weights,
		      int is_diff_map) {
    
   graphics_info_t g;
   int imol = -1; // failure initially.
   struct stat buf;

   std::string f_col_str(f_col);
   std::string phi_col_str(phi_col);
   std::string weight_col_str("");
   if (use_weights)
      weight_col_str = std::string(weight_col);

   int status = stat(mtz_file_name, &buf);
   // 
   if (status != 0) {
      std::cout << "WARNING:: Can't find file " << mtz_file_name << std::endl;
      if (S_ISDIR(buf.st_mode)) {
	 std::cout << mtz_file_name << " is a directory! " << std::endl;
      }
   } else {

      if (valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)) { 
      
	 std::cout << "making map from mtz filename " << mtz_file_name << std::endl;
	 imol = g.n_molecules;
	 g.molecules[imol].map_fill_from_mtz(std::string(mtz_file_name),
					     f_col_str,
					     phi_col_str,
					     weight_col_str,
					     use_weights, is_diff_map);
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 graphics_draw();
	 g.activate_scroll_radio_button_in_display_manager(imol);
	 std::vector<std::string> command_strings;
	 command_strings.push_back("make-and-draw-map");
	 command_strings.push_back(single_quote(mtz_file_name));
	 command_strings.push_back(single_quote(f_col));
	 command_strings.push_back(single_quote(phi_col));
	 command_strings.push_back(single_quote(weight_col));
	 command_strings.push_back(graphics_info_t::int_to_string(use_weights));
	 command_strings.push_back(graphics_info_t::int_to_string(is_diff_map));
	 add_to_history(command_strings);
	 
      } else {
	 std::cout << "WARNING:: label(s) not found in mtz file " 
		   << mtz_file_name << " " << f_col_str << " " 
		   <<  phi_col_str << " ";
	 if (use_weights)
	    std::cout << weight_col_str << std::endl;
	 else 
	    std::cout << std::endl;
      }      
   }
   return imol; // possibly -1
}

int  make_and_draw_map_with_refmac_params(const char *mtz_file_name, 
					  const char *a, const char *b,
					  const char *weight,
					  int use_weights, int is_diff_map,
					  short int have_refmac_params,
					  const char *fobs_col,
					  const char *sigfobs_col,
					  const char *r_free_col,
					  short int sensible_f_free_col) { 

   graphics_info_t g;
   int imol = g.n_molecules;

   // this is order dependent.  the restore-state comand that is
   // constructed in make_and_draw_map checks to see if we
   // have_refmac_params, so we need to set them before making the map.
   // 
   if (have_refmac_params)
      g.molecules[imol].store_refmac_params(std::string(mtz_file_name),
					    std::string(fobs_col), 
					    std::string(sigfobs_col), 
					    std::string(r_free_col),
					    sensible_f_free_col);
   
   return make_and_draw_map(mtz_file_name, a, b, weight, use_weights, is_diff_map);
}

// return imol, possibly -1;
int make_and_draw_map_with_reso_with_refmac_params(const char *mtz_file_name, 
						   const char *f_col,
						   const char *phi_col,
						   const char *weight_col,
						   int use_weights, int is_diff_map,
						   short int have_refmac_params,
						   const char *fobs_col,
						   const char *sigfobs_col,
						   const char *r_free_col,
						   short int sensible_r_free_col,
						   short int is_anomalous_flag,
						   short int use_reso_limits,
						   float low_reso_limit,
						   float high_reso_limit) {

   graphics_info_t g;
   int imol = -1;

   // this is order dependent.  the restore-state comand that is
   // constructed in make_and_draw_map checks to see if we
   // have_refmac_params, so we need to set them before making the map.
   // 
   if (have_refmac_params)
      g.molecules[g.n_molecules].store_refmac_params(std::string(mtz_file_name),
					    std::string(fobs_col), 
					    std::string(sigfobs_col), 
					    std::string(r_free_col),
					    sensible_r_free_col);
   
   struct stat buf;
   int status = stat(mtz_file_name, &buf);

   if (status != 0) {
      std::cout << "Error finding MTZ file " << mtz_file_name << std::endl;
      if (S_ISDIR(buf.st_mode)) {
	 std::cout << mtz_file_name << " is a directory! " << std::endl;
      }
   } else {
      std::string map_type;
      if (is_diff_map)
	 map_type = "difference";
      else
	 map_type = "conventional";
      
      std::cout << "making " << map_type << " map from MTZ filename "
		<< mtz_file_name << " using " << f_col << " "
		<< phi_col << std::endl;

      if (valid_labels(mtz_file_name, f_col, phi_col, weight_col, use_weights)) {
	 std::string weight_col_str("");
	 if (use_weights)
	    weight_col_str = std::string(weight_col);
	 imol = g.n_molecules;
	 g.molecules[imol].map_fill_from_mtz_with_reso_limits(std::string(mtz_file_name),
							      std::string(f_col),
							      std::string(phi_col),
							      weight_col_str,
							      use_weights,
							      is_anomalous_flag,
							      is_diff_map,
							      use_reso_limits,
							      low_reso_limit,
							      high_reso_limit);
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 g.activate_scroll_radio_button_in_display_manager(imol);
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: label(s) not found in MTZ file " 
		   << mtz_file_name << " " << f_col << " " 
		   <<  phi_col << " ";
	 if (use_weights)
	    std::cout << weight_col << std::endl;
	 else 
	    std::cout << std::endl;
      }
   }
   if (imol != -1) {

      // We reset some strings if we weren't given refmac params -
      // otherwise we quote garbage or unallocated memory.
      std::string weight_col_str;
      std::string fobs_col_str;
      std::string sigfobs_col_str;
      std::string r_free_col_str;
      if (weight_col)
	 weight_col_str = single_quote(weight_col);
      else
	 weight_col_str = single_quote("Weight:None-specified");
      
      if (! have_refmac_params) {
	 fobs_col_str    = single_quote("Fobs:None-specified");
	 sigfobs_col_str = single_quote("SigF:None-specified");
	 r_free_col_str  = single_quote("RFree:None-specified");
	 sensible_r_free_col = 0;
      } else {
	 fobs_col_str    = single_quote(fobs_col);
	 sigfobs_col_str = single_quote(sigfobs_col_str);
	 r_free_col_str  = single_quote(r_free_col);
      }
      
      std::vector<std::string> command_strings;
      command_strings.push_back("make-and-draw-map-with-reso-with-refmac-params");
      command_strings.push_back(single_quote(mtz_file_name));
      command_strings.push_back(single_quote(f_col));
      command_strings.push_back(single_quote(phi_col));
      command_strings.push_back(weight_col_str);
      command_strings.push_back(graphics_info_t::int_to_string(use_weights));
      command_strings.push_back(graphics_info_t::int_to_string(is_diff_map));
      command_strings.push_back(graphics_info_t::int_to_string(have_refmac_params));
      command_strings.push_back(fobs_col_str);
      command_strings.push_back(sigfobs_col_str);
      command_strings.push_back(r_free_col_str);
      command_strings.push_back(graphics_info_t::int_to_string(sensible_r_free_col));
      command_strings.push_back(graphics_info_t::int_to_string(is_anomalous_flag));
      command_strings.push_back(graphics_info_t::int_to_string(use_reso_limits));
      command_strings.push_back(graphics_info_t::float_to_string(low_reso_limit));
      command_strings.push_back(graphics_info_t::float_to_string(high_reso_limit));
      add_to_history(command_strings);
   }
   return imol;
}


int auto_read_make_and_draw_maps(const char *mtz_file_name) {

   const char *f_col   = graphics_info_t::auto_read_MTZ_FWT_col.c_str();
   const char *phi_col = graphics_info_t::auto_read_MTZ_PHWT_col.c_str();
   const char *weight_col = "";
   
   const char *f_col_diff   = graphics_info_t::auto_read_MTZ_DELFWT_col.c_str();
   const char *phi_col_diff = graphics_info_t::auto_read_MTZ_PHDELWT_col.c_str();
   int imol1, imol2;
   
   imol1 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							  f_col,
							  phi_col,
							  weight_col,
							  0, 0, //    use_weights,  is_diff_map,
							  0,  //   short int have_refmac_params,
							  "", //   const char *fobs_col,
							  "", //   const char *sigfobs_col,
							  "", //   const char *r_free_col,
							  0,  //   short int sensible_f_free_col,
							  0,  //   short int is_anomalous_flag,
							  0,  //   short int use_reso_limits,
							  0,  //   float low_reso_limit,
							  0); //   float high_reso_limit
   
   imol2 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							  f_col_diff,
							  phi_col_diff,
							  weight_col,
							  0, 1, //    use_weights,  is_diff_map,
							  0,  //   short int have_refmac_params,
							  "", //   const char *fobs_col,
							  "", //   const char *sigfobs_col,
							  "", //   const char *r_free_col,
							  0,  //   short int sensible_f_free_col,
							  0,  //   short int is_anomalous_flag,
							  0,  //   short int use_reso_limits,
							  0,  //   float low_reso_limit,
							  0); //   float high_reso_limit
   
   std::string s;
   if (imol1 < 0) { 
      s += "Failed to find columns ";
      s += graphics_info_t::auto_read_MTZ_FWT_col;
      s += " and ";
      s += graphics_info_t::auto_read_MTZ_PHWT_col;
      s += "in that mtz file\n";
   }
   if (imol2 < 0) { 
      s += "Failed to find columns ";
      s += f_col_diff;
      s += " and ";
      s += phi_col_diff;
      s += " in that mtz file\n";
   }
   if (imol1<0 || imol2 <0) {

      imol1 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							     "2FOFCWT",
							     "PH2FOFCWT",
							     weight_col,
							     0, 0, //    use_weights,  is_diff_map,
							     0,  //   short int have_refmac_params,
							     "", //   const char *fobs_col,
							     "", //   const char *sigfobs_col,
							     "", //   const char *r_free_col,
							     0,  //   short int sensible_f_free_col,
							     0,  //   short int is_anomalous_flag,
							     0,  //   short int use_reso_limits,
							     0,  //   float low_reso_limit,
							     0); //   float high_reso_limit
      
      imol2 = make_and_draw_map_with_reso_with_refmac_params(mtz_file_name, 
							     "FOFCWT",
							     "PHFOFCWT",
							     weight_col,
							     0, 1, //    use_weights,  is_diff_map,
							     0,  //   short int have_refmac_params,
							     "", //   const char *fobs_col,
							     "", //   const char *sigfobs_col,
							     "", //   const char *r_free_col,
							     0,  //   short int sensible_f_free_col,
							     0,  //   short int is_anomalous_flag,
							     0,  //   short int use_reso_limits,
							     0,  //   float low_reso_limit,
							     0); //   float high_reso_limit
   

      if (imol1<0 || imol2 <0) {
	 GtkWidget *w = wrapped_nothing_bad_dialog(s);
	 gtk_widget_show(w);
      }
   }
   
   // we like the 2fo-fc map to be the scrolled one after autoreading,
   if (imol1 >= 0) { 
      graphics_info_t g;
      g.scroll_wheel_map = imol1;
      g.activate_scroll_radio_button_in_display_manager(imol1);
   }
   
   return imol2; // return to callbacks.c (currently ignored)
} 

int auto_read_do_difference_map_too_state() {

   add_to_history_simple("auto-read-do-difference-map-too-state");
   int i = graphics_info_t::auto_read_do_difference_map_too_flag; 
   return i;

} 
void set_auto_read_do_difference_map_too(int i) {

   graphics_info_t::auto_read_do_difference_map_too_flag = i;
   std::string cmd = "set-auto-read-do-dfference-map-too";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
} 


void set_auto_read_column_labels(const char *fwt, const char *phwt, 
				 int is_for_diff_map_flag) {

   if (is_for_diff_map_flag) {
      graphics_info_t::auto_read_MTZ_DELFWT_col = fwt;
      graphics_info_t::auto_read_MTZ_PHDELWT_col = phwt;
   } else {
      graphics_info_t::auto_read_MTZ_FWT_col = fwt;
      graphics_info_t::auto_read_MTZ_PHWT_col = phwt;
   }
   std::string cmd = "set-auto-read-column-labels";
   std::vector<coot::command_arg_t> args;
   args.push_back(fwt);
   args.push_back(phwt);
   args.push_back(is_for_diff_map_flag);
   add_to_history_typed(cmd, args);

}


void toggle_idle_function() { 

   graphics_info_t g; 

   if (g.i_fn_token == 0) { 
      g.i_fn_token = gtk_idle_add((GtkFunction)animate_idle, g.glarea);
   } else {
      gtk_idle_remove(g.i_fn_token);
      g.i_fn_token = 0; 
   }
   add_to_history_simple("toggle_idle_function");
}

// in degrees
void set_idle_function_rotate_angle(float f) {

   graphics_info_t g;
   std::string cmd = "set-idle-function-rotate-angle";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
   g.idle_function_rotate_angle = f; 
}

float idle_function_rotate_angle() {

   std::string cmd = "idle-function-rotate-angle";
   add_to_history_simple(cmd);
   return graphics_info_t::idle_function_rotate_angle;
} 



void do_tw() {

}

// Another name for wrapped_nothing_bad_dialog, but this function also
// displays the widget so nothing is returned.
// 
void info_dialog(const char *txt) {

   graphics_info_t::info_dialog(txt);
   std::string cmd = "info-dialog";
   std::vector<coot::command_arg_t> args;
   args.push_back(txt);
   add_to_history_typed(cmd, args);
}


GtkWidget *main_menubar() {

   GtkWidget *w = lookup_widget(graphics_info_t::statusbar, "menubar1");
   return w; 
}

GtkWidget *main_statusbar() {
   return graphics_info_t::statusbar;
} 




/*  ------------------------------------------------------------------------ */
/*                         file selection                                    */
/*  ------------------------------------------------------------------------ */

void
set_directory_for_fileselection(GtkWidget *fileselection1) { 
   graphics_info_t g;
   g.set_directory_for_fileselection(fileselection1);
}

void
save_directory_from_fileselection(const GtkWidget *fileselection) {
   graphics_info_t g;
   g.save_directory_from_fileselection(fileselection);
}

void save_directory_for_saving_from_fileselection(const GtkWidget *fileselection) {
   graphics_info_t g;
   g.save_directory_for_saving_from_fileselection(fileselection);
} 


/*  Eleanor likes to sort her files by date when selecting a file
*/
GtkWidget *add_sort_button_fileselection(GtkWidget *fileselection) {

   GtkWidget *aa = GTK_FILE_SELECTION(fileselection)->action_area;
   GtkWidget *frame = gtk_frame_new("File Order");
   GtkWidget *button = gtk_button_new_with_label("  Sort by Date  ");
   gtk_widget_ref(button);
   gtk_object_set_data_full(GTK_OBJECT(fileselection),
			    "fileselection_sort_button",
			    button,
			    (GtkDestroyNotify) gtk_widget_unref);
   GtkWidget *file_list = GTK_FILE_SELECTION(fileselection)->file_list;
 
   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(GTK_FILE_SELECTION(fileselection)->history_pulldown);

   gtk_object_set_user_data(GTK_OBJECT(button), (char *) history_pulldown); 


#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
   gtk_signal_connect (GTK_OBJECT(button), "clicked",
		       (GtkSignalFunc) fileselection_sort_button_clicked_gtk1,
		       file_list);
#else
   gtk_signal_connect (GTK_OBJECT(button), "clicked",
		       (GtkSignalFunc) fileselection_sort_button_clicked,
		       file_list);
#endif

   gtk_container_add(GTK_CONTAINER(aa),frame);
   gtk_container_add(GTK_CONTAINER(frame), button);
   gtk_widget_show(frame);
   gtk_widget_show(button);
   return button;
}

bool compare_mtimes(coot::str_mtime a, coot::str_mtime b) {
   return a.mtime > b.mtime;
}



#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
void fileselection_sort_button_clicked_gtk1( GtkWidget *sort_button,
					     GtkCList  *file_list) {

   std::vector<coot::str_mtime> v;
   char *text;
   struct stat buf;
   time_t mtime;

   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(gtk_object_get_user_data(GTK_OBJECT(sort_button)));

   // The menu item is a container than contains a label.
   // How do we get to the label given the container?
   // Strangely enough we use the history_pulldown.

   GList *dlist = gtk_container_children(GTK_CONTAINER(history_pulldown));
   GList *free_list = dlist;
   std::string pre_directory("");
   
   while (dlist) {
      // GtkWidget *list_item;
      // list_item = (GtkWidget *) (dlist->data);
      gchar *t = GTK_LABEL(dlist->data)->label;
      if (t != NULL) {
	 pre_directory = t; 
      } else {
	 std::cout << "null label t " << std::endl;
      } 
      dlist = dlist->next;
   }
   g_list_free(free_list); 

   int status;
   std::string file_name;
   for(int i=0; i<file_list->rows; i++) {
      gtk_clist_get_text(file_list,i,0,&text);
      file_name = pre_directory;
      if (file_name != "") {
	 file_name += "/";
	 file_name += text;
      } else {
	 // Oh dear.  This directory only.  Else failure
	 // this shouldn't happen hopefully.
	 file_name = text;
      } 
      status = stat(file_name.c_str(),&buf);
      if (status == 0) { 
	 mtime = buf.st_mtime;
	 v.push_back(coot::str_mtime(text, mtime));
      } else {
	 std::cout << "error stating " << file_name << std::endl;
      }
   }
//    v.push_back(coot::str_mtime("Cabbage",324232));
//    v.push_back(coot::str_mtime("zipped",345234523452));
   std::sort(v.begin(), v.end(), compare_mtimes);

   gtk_clist_clear(file_list);

   for(unsigned int i=0; i<v.size(); i++) {
      text = (char *) v[i].file.c_str(); 
      gtk_clist_append(file_list, &text);
   }
}

#endif // GTK_1 etc

void quanta_buttons() {
   graphics_info_t g;
   g.quanta_buttons();
   add_to_history_simple("quanta-buttons");
}

void quanta_like_zoom() { 
   graphics_info_t::quanta_like_zoom_flag = 1;
   add_to_history_simple("quanta-like-zoom");
} 

void set_scroll_by_wheel_mouse(int istate) {
   graphics_info_t::do_scroll_by_wheel_mouse_flag = istate;
   std::string cmd = "set-scroll-by-mouse-wheel";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
}

int scroll_by_wheel_mouse_state() {
   add_to_history_simple("scroll-by-wheel-mouse-state");
   return graphics_info_t::do_scroll_by_wheel_mouse_flag;
}

/*! \brief set the default inital contour for 2FoFc-style map

in sigma */
void set_default_initial_contour_level_for_map(float n_sigma) {

   graphics_info_t::default_sigma_level_for_map = n_sigma;
   std::string cmd = "set-default-initial-contour-level-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(n_sigma);
   add_to_history_typed(cmd, args);
} 

/*! \brief set the default inital contour for FoFc-style map

in sigma */
void set_default_initial_contour_level_for_difference_map(float n_sigma) {

   graphics_info_t::default_sigma_level_for_fofc_map = n_sigma;
   std::string cmd = "set-default-initial-contour-level-for-difference-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(n_sigma);
   add_to_history_typed(cmd, args);
} 




void set_map_line_width(int w) {
   graphics_info_t::map_line_width = w;
   // update the maps because they may be being draw as graphical
   // objects.
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++)
      graphics_info_t::molecules[imol].update_map();
   graphics_draw();
   std::string cmd = "set-map-line-width";
   std::vector<coot::command_arg_t> args;
   args.push_back(w);
   add_to_history_typed(cmd, args);
   
}

int map_line_width_state() {
   add_to_history_simple("map-line-width-state");
   return graphics_info_t::map_line_width;
} 


//
char* get_text_for_symmetry_size_widget() {

   // convert the density in the graphics_info_t

   graphics_info_t g;
   char *text;

   // we are interfacing with a c function, so we will need to malloc
   // the space for the returned char*.
   //
   // The c function should delete it.
   // 
   // Dontcha just *love" this sort coding! (yeuch)
   // 
   text = (char *) malloc(100);
   snprintf(text,100,"%-5.1f", g.symmetry_search_radius);

   return text;

}

void set_density_size_from_widget(const char *text) {

   float tmp;
   graphics_info_t g;

   tmp = atof(text);

   if ((tmp > 0.0) && (tmp < 9999.9)) {
      g.box_radius = tmp;
   } else {

      cout << "Cannot interpret " << text << ".  Assuming 10A" << endl;
      g.box_radius = 10.0;
   }
   //
   for (int ii=0; ii<g.n_molecules; ii++) {
      if (is_valid_map_molecule(ii))
	  g.molecules[ii].update_map();
   }
   graphics_draw();
}

void set_density_size(float f) {

   graphics_info_t g;
   g.box_radius = f;
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_map();
   }
   graphics_draw();
   std::string cmd = "set-density-size";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
   
}

/*! \brief set the extent of the box/radius of electron density contours */
void set_map_radius(float f) {
   set_density_size(f);
} 



void set_display_intro_string(const char *str) {

   if (str) { 
      if (graphics_info_t::use_graphics_interface_flag) { 
	 std::string s(str);
	 graphics_info_t g;
	 g.display_density_level_screen_string = s;
	 g.statusbar_text(s);
      }

      std::string cmd = "set-display-intro-string";
      std::vector<coot::command_arg_t> args;
      args.push_back(str);
      add_to_history_typed(cmd, args);
   }
} 

void set_swap_difference_map_colours(int i) { 
   graphics_info_t::swap_difference_map_colours = i;
   std::string cmd = "set-swap-difference-map-colours";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
   
}

/* return success status 0 = failure (imol does not have a map) */
int set_map_is_difference_map(int imol) { 

   int istatus = 0;
   if (imol< graphics_n_molecules()) { 
      if (graphics_info_t::molecules[imol].has_map()) { 
	 graphics_info_t::molecules[imol].set_map_is_difference_map();
	 istatus = 1;
	 graphics_draw();
      } else { 
	 std::cout << "WARNING:: molecule " << imol << " does not have a map." <<  std::endl;
      } 

   } else { 
      std::cout << "WARNING:: No such molecule as " << imol << std::endl;
   }
   std::string cmd = "set-map-is-difference-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   
   return istatus;
} 

/* return the index of the new molecule or -1 on failure */
int another_level() {

   int istat = -1;
   int imap = -1;
   for (int i=0; i<graphics_info_t::n_molecules; i++) {
      if (is_valid_map_molecule(i)) {
	 if (! graphics_info_t::molecules[i].is_difference_map_p()) { 
	    imap = i;
	 }
      }
   }
   if (imap > -1) {
      istat = another_level_from_map_molecule_number(imap);
   }
   // history elsewhere
   return istat;
}

int another_level_from_map_molecule_number(int imap) {
   int istat = -1;
   if (is_valid_map_molecule(imap)) {
      // create another map with the same parameters as imap and then
      // push up the contour level a sigma.
//       std::cout << "DEBUG:: calling make_and_draw_map_with_reso_with_refmac_params"
// 		<< std::endl;
      istat = make_and_draw_map_with_reso_with_refmac_params(
	  graphics_info_t::molecules[imap].save_mtz_file_name.c_str(),
	  graphics_info_t::molecules[imap].save_f_col.c_str(),
	  graphics_info_t::molecules[imap].save_phi_col.c_str(),
          graphics_info_t::molecules[imap].save_weight_col.c_str(),
          graphics_info_t::molecules[imap].save_use_weights,
          graphics_info_t::molecules[imap].save_is_diff_map_flag,
	  0, "None", "None", "None", 0, // refmac params
          graphics_info_t::molecules[imap].save_is_anomalous_map_flag,
          graphics_info_t::molecules[imap].save_use_reso_limits,
          graphics_info_t::molecules[imap].save_low_reso_limit,
          graphics_info_t::molecules[imap].save_high_reso_limit);

      if (istat != -1) { 

	 float map_sigma = graphics_info_t::molecules[istat].map_sigma();
	 float current_contour_level = graphics_info_t::molecules[istat].contour_level[0];
	 graphics_info_t::molecules[istat].set_contour_level(current_contour_level +
							     map_sigma*1.0);
	 graphics_info_t::molecules[istat].update_map();
	 graphics_draw();
      }
   }
   std::string cmd = "another-level-from-map-molecule-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(imap);
   add_to_history_typed(cmd, args);
   return istat;
}



void set_map_radius_slider_max(float f) {
   graphics_info_t::map_radius_slider_max = f;
   std::string cmd = "set-map-radius-slider-max";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
} 


// return 1 on "yes, it has a cell".
// 
int has_unit_cell_state(int imol) { 

   int istate = 0;
   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model() ||
	     graphics_info_t::molecules[imol].has_map()) {
	    istate = graphics_info_t::molecules[imol].have_unit_cell;
	 }
      }
   }
   std::string cmd = "has-unit-cell-state";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return istate;
}
 
void set_symmetry_size_from_widget(const char *text) {
   
   float tmp;
   graphics_info_t g;

   tmp = atof(text);

   if ((tmp > 0.0) && (tmp < 9999.9)) {
      g.symmetry_search_radius = tmp;
   } else {

      cout << "Cannot interpret " << text << ".  Assuming 10A" << endl;
      g.symmetry_search_radius = 10.0;
   }
   //
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_symmetry();
   }
   graphics_draw();
   
}

void set_symmetry_size(float f) {
   graphics_info_t g;
   g.symmetry_search_radius = f;
   for (int ii=0; ii<g.n_molecules; ii++) {
      g.molecules[ii].update_symmetry();
   }
   graphics_draw();
   std::string cmd = "set-symmetry-size";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
}

/* When the coordinates for one (or some) symmetry operator are missing
   (which happens sometimes, but rarely), try changing setting this to 2
   (default is 1).  It slows symmetry searching, which is why it is not
   set to 2 by default.  */
void set_symmetry_shift_search_size(int shift) {
   
   graphics_info_t::symmetry_shift_search_size = shift;
   std::string cmd = "set-symmetry-shift-search-size";
   std::vector<coot::command_arg_t> args;
   args.push_back(shift);
   add_to_history_typed(cmd, args);
}


void set_symmetry_molecule_rotate_colour_map(int imol, int state) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      g.molecules[imol].symmetry_rotate_colour_map_flag = state;
   }
   std::string cmd = "set-symmetry-molecule-rotate-colour-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(state);
   add_to_history_typed(cmd, args);
   graphics_draw();
} 

int symmetry_molecule_rotate_colour_map_state(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].symmetry_rotate_colour_map_flag;
   }
   std::string cmd = "symmetry-molecule-rotate-colour-map-state";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r;
} 

void set_symmetry_colour_by_symop(int imol, int state) {

   if (graphics_info_t::use_graphics_interface_flag) { 
      graphics_info_t g;
      if (is_valid_model_molecule(imol)) { 
	 g.molecules[imol].symmetry_colour_by_symop_flag = state;
	 graphics_draw();
      }
   }
   std::string cmd = "set-symmetry-colour-by-symop";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(state);
   add_to_history_typed(cmd, args);
}

void set_symmetry_whole_chain(int imol, int state) {
   
   if (graphics_info_t::use_graphics_interface_flag) { 
      graphics_info_t g;
      if (is_valid_model_molecule(imol)) { 
	 g.molecules[imol].symmetry_whole_chain_flag = state;
	 if (g.glarea)
	    g.update_things_on_move_and_redraw();
      }
   }
   std::string cmd = "set-symmetry-whole-chain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(state);
   add_to_history_typed(cmd, args);
}




void set_fps_flag(int thing) {

   graphics_info_t g;
   g.SetShowFPS(thing);
   std::string cmd = "set-fps-flag";
   std::vector<coot::command_arg_t> args;
   args.push_back(thing);
   add_to_history_typed(cmd, args);
} 

// For people without PCs with fast graphics cards :)  [like me]
// 
void set_active_map_drag_flag(int t) {

   graphics_info_t g;
   g.SetActiveMapDrag(t);
   std::string cmd = "set-active-map-drag-flag";
   std::vector<coot::command_arg_t> args;
   args.push_back(t);
   add_to_history_typed(cmd, args);
}

int get_fps_flag() {

   graphics_info_t g;
   add_to_history_simple("get-fps-flag");
   return g.GetFPSFlag();
} 


//
short int get_active_map_drag_flag() {

   graphics_info_t g;
   add_to_history_simple("get-active-map-drag-flag");
   return g.GetActiveMapDrag();
}

void set_draw_hydrogens(int imol, int istate) {

   graphics_info_t g;
   
   if ((imol < g.n_molecules) && (imol >= 0)) {
      g.molecules[imol].set_draw_hydrogens_state(istate);
      graphics_draw();
   } else { 
      std::cout << "WARNING:: No such molecule number " << imol << "\n";
   }
   std::string cmd = "set-draw-hydrogens";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(istate);
   add_to_history_typed(cmd, args);
} 

void set_show_origin_marker(int istate) {
   graphics_info_t::show_origin_marker_flag = istate;
   graphics_draw();
   std::string cmd = "set-show-origin-marker";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
} 

int  show_origin_marker_state() {
   add_to_history_simple("show-origin-marker-state");
   return graphics_info_t::show_origin_marker_flag;
} 





void set_last_map_contour_level(float level) {

   graphics_info_t g;
   g.set_last_map_contour_level(level);
   std::string cmd = "set-last-map-contour-level";
   std::vector<coot::command_arg_t> args;
   args.push_back(level);
   add_to_history_typed(cmd, args);
}

void set_last_map_contour_level_by_sigma(float level) {

   graphics_info_t g;
   g.set_last_map_contour_level_by_sigma(level);
   std::string cmd = "set-last-map-contour-level-by-sigma";
   std::vector<coot::command_arg_t> args;
   args.push_back(level);
   add_to_history_typed(cmd, args);
}

void set_last_map_sigma_step(float f) { 

   graphics_info_t g;
   g.set_last_map_sigma_step(f);
   std::string cmd = "set-last-map-sigma-step";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);

} 

 


void
handle_symmetry_colour_change(int mol, gdouble* col) {

   //
   graphics_info_t::symm_colour[0][0] = col[0];
   graphics_info_t::symm_colour[0][1] = col[1];
   graphics_info_t::symm_colour[0][2] = col[2];

   graphics_draw();

}

gdouble*
get_map_colour(int imol) {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   if (imol < graphics_info_t::n_molecules) { 
      if (graphics_info_t::molecules[imol].has_map()) { 
	 colour[0] = graphics_info_t::molecules[imol].map_colour[0][0]; 
	 colour[1] = graphics_info_t::molecules[imol].map_colour[0][1]; 
	 colour[2] = graphics_info_t::molecules[imol].map_colour[0][2];
      }
   }
   std::string cmd = "get-map-colour";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return colour;
}

//! \brief return the colour of the imolth map (e.g.: (list 0.4 0.6
//0.8). If invalid imol return #f.
// 
#ifdef USE_GUILE
SCM map_colour_components(int imol) {

   SCM r = SCM_BOOL(0);
   if (is_valid_map_molecule(imol)) {
      double rc = graphics_info_t::molecules[imol].map_colour[0][0]; 
      double gc = graphics_info_t::molecules[imol].map_colour[0][1]; 
      double bc = graphics_info_t::molecules[imol].map_colour[0][2]; 
      r = SCM_CAR(scm_listofnull);
      // put red at the front of the resulting list
      r = scm_cons(scm_double2num(bc), r);
      r = scm_cons(scm_double2num(gc), r);
      r = scm_cons(scm_double2num(rc), r);
   }
   std::string cmd = "map-colour-components";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r; 
}
#endif


/* Functions for Cancel button on map colour selection  */
void save_previous_map_colour(int imol) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].save_previous_map_colour();
   } 
   std::string cmd = "save-previous-map-colour";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);

} 

void restore_previous_map_colour(int imol) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].restore_previous_map_colour();
   }
   graphics_draw();
   std::string cmd = "restore-previous-map-colour";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}


// -------------------------------------------------------------------

gdouble*
get_symmetry_bonds_colour(int idummy) {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   colour[0] = graphics_info_t::symm_colour[0][0];
   colour[1] = graphics_info_t::symm_colour[0][1];
   colour[2] = graphics_info_t::symm_colour[0][2];
   return colour;
}





// In future the gui will usefully set a mol number and we
// will use that.
//
void set_show_symmetry_master(short int state) {

   // 
   graphics_info_t g;

   // show symmetry state is no longer part of the molecule(s).
   
      
//       g.molecules[ii].show_symmetry = state;

//       if ( state == 1 ) {
// 	 g.molecules[mol_no].update_symmetry();
//       }
//    }

   g.show_symmetry = state; 
   for (int ii=0; ii<g.n_molecules; ii++)
      if (is_valid_model_molecule(ii))
	 graphics_info_t::molecules[ii].update_symmetry();
   graphics_draw();

   if (state) { 
      // Now count the number of model molecules that have symmetry
      // available.  If there are none, then pop up a warning.

      int n_has_symm = 0;
      int n_model_molecules = 0;
      for (int ii=0; ii<g.n_molecules; ii++)
	 if (is_valid_model_molecule(ii)) {
	    n_model_molecules++;
	    mat44 my_matt;
	    int err = graphics_info_t::molecules[ii].atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
	    if (err == SYMOP_Ok) {
	       n_has_symm++;
	       break;
	    }
	 }
      if ((n_has_symm == 0) && (n_model_molecules > 0)) {
	 std::string s = "WARNING:: there are no model molecules\n"; 
	 s += " that can display symmetry.  \n\nCRYST1 problem?";
	 if (graphics_info_t::use_graphics_interface_flag) { 
	    GtkWidget *w = g.wrapped_nothing_bad_dialog(s);
	    gtk_widget_show(w);
	 }
      }
   }
   std::string cmd = "set-show-symmetry-master";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}

void set_show_symmetry_molecule(int imol, short int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_symmetry(state);
      if (state)
	 graphics_info_t::molecules[imol].update_symmetry();
      graphics_draw();
   }
   std::string cmd = "set-show-symmetry-molecule";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


void symmetry_as_calphas(int mol_no, short int state) {

   graphics_info_t g;
   if (is_valid_model_molecule(mol_no)) { 
      g.molecules[mol_no].symmetry_as_calphas = state;
      g.molecules[mol_no].update_symmetry();
   }
   graphics_draw();
   std::string cmd = "symmetry-as-calphas";
   std::vector<coot::command_arg_t> args;
   args.push_back(mol_no);
   args.push_back(state);
   add_to_history_typed(cmd, args);

}

short int get_symmetry_as_calphas_state(int imol) {

   graphics_info_t g;
   int r = -1;
   if (is_valid_model_molecule(imol))
      r = g.molecules[imol].symmetry_as_calphas;
       
   std::string cmd = "get-symmety-as-calphas-state";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r;
} 


//  There is no dependence on the mol_no.  The intereface where we
//  pass (int mol_no) is kept, because molecule symmetry dependence
//  may re-surface in future.
// 
short int get_show_symmetry() {

   add_to_history_simple("get-show-symmetry");
   return graphics_info_t::show_symmetry; // master
}

GtkWidget *symmetry_molecule_controller_dialog() {

   graphics_info_t g;
   return g.wrapped_create_symmetry_controller_dialog();
} 

   
void
set_clipping_front(float v) {

   graphics_info_t::clipping_front = v;
   graphics_draw();
   std::string cmd = "set-clipping-front";
   std::vector<coot::command_arg_t> args;
   args.push_back(v);
   add_to_history_typed(cmd, args);
   
}


void
set_clipping_back(float v) {

   graphics_info_t::clipping_back = v;
   graphics_draw();
   std::string cmd = "set-clipping-back";
   std::vector<coot::command_arg_t> args;
   args.push_back(v);
   add_to_history_typed(cmd, args);
}


/*  ----------------------------------------------------------------------- */
/*                         Colour                                           */
/*  ----------------------------------------------------------------------- */

void
set_symmetry_colour_merge(int mol_no, float v) {

   graphics_info_t::symm_colour_merge_weight[mol_no] = v;
   graphics_draw();
   std::string cmd = "";
   std::vector<coot::command_arg_t> args;
   args.push_back(mol_no);
   args.push_back(v);
   add_to_history_typed(cmd, args);
}

void set_colour_map_rotation_on_read_pdb(float f) {
   graphics_info_t::rotate_colour_map_on_read_pdb = f; 
   std::string cmd = "set-colour-map-rotation-on-read-pdb";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
}

void set_colour_map_rotation_on_read_pdb_flag(short int i) {
   graphics_info_t::rotate_colour_map_on_read_pdb_flag = i; 
   std::string cmd = "set-colour-map-rotation-on-read-pdb-flag";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}

void set_colour_map_rotation_on_read_pdb_c_only_flag(short int i) {

   graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag = i;
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      if (is_valid_model_molecule(imol)) {
	 if (graphics_info_t::molecules[imol].Bonds_box_type() == coot::COLOUR_BY_CHAIN_BONDS) {
	    graphics_info_t::molecules[imol].make_bonds_type_checked();
	 }
      }
   }
   graphics_draw();
   std::string cmd = "set-colour-map-rotation-on-read-pdb-c-only-flag";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}

void set_symmetry_atom_labels_expanded(int state) {
   graphics_info_t::symmetry_atom_labels_expanded_flag = state;
   graphics_draw();
   std::string cmd = "set-symmetry-atom-labels-expanded";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


/* widget work */
GtkWidget *wrapped_create_coords_colour_control_dialog() {

   GtkWidget *w = create_coords_colour_control_dialog();

   graphics_info_t g;
   g.fill_bond_colours_dialog_internal(w);
   return w;
}


float get_molecule_bonds_colour_map_rotation(int imol) {

   float r = -1.0;
   if (is_valid_model_molecule(imol))
      r = graphics_info_t::molecules[imol].bonds_colour_map_rotation;
   std::string cmd = "get-molecule-bonds-colour-map-rotation";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r;
}

void  set_molecule_bonds_colour_map_rotation(int imol, float f) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].bonds_colour_map_rotation = f;
   std::string cmd = "set-molecule-bonds-colour-map-rotation";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(f);
   add_to_history_typed(cmd, args);
}



// ----------------- Rotation Centre ----------------------

void set_rotation_centre(float x, float y, float z) {
   graphics_info_t g;
   g.setRotationCentre(coot::Cartesian(x,y,z));
   if (g.glarea)
      g.update_things_on_move_and_redraw();
   std::string cmd = "set-rotation-centre";
   std::vector<coot::command_arg_t> args;
   args.push_back(x);
   args.push_back(y);
   args.push_back(z);
   add_to_history_typed(cmd, args);
}

// The redraw happens somewhere else...
void set_rotation_centre_internal(float x, float y, float z) {
   graphics_info_t g;
   g.setRotationCentre(coot::Cartesian(x,y,z));
}

float rotation_centre_position(int axis) {  /* only return one value: x=0, y=1, z=2 */
   graphics_info_t g;
   coot::Cartesian p = g.RotationCentre();
   // std::cout << "DEBUG:: rotation centre : " << p << " axis: " << axis << "\n";
   float r = 0.0;
   if (axis == 0)
      r = p.x();
   if (axis == 1)
      r = p.y();
   if (axis == 2)
      r = p.z();
   std::string cmd = "rotation-centre-position";
   std::vector<coot::command_arg_t> args;
   args.push_back(axis);
   add_to_history_typed(cmd, args);
   return r;
}

void set_colour_by_chain(int imol) { 
   
   if (is_valid_model_molecule(imol)) {
      short int f = graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag;
      graphics_info_t::molecules[imol].make_colour_by_chain_bonds(f);
      graphics_draw();
   }
   std::string cmd = "set-colour-by-chain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

void set_colour_by_molecule(int imol) { 

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].make_colour_by_molecule_bonds();
      graphics_draw();
   }
   std::string cmd = "set-colour-by-molecule";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}


/*  Section Map colour*/
/* default for maps is 31 degrees. */
void set_colour_map_rotation_for_map(float f) {

   graphics_info_t::rotate_colour_map_for_map = f;
   std::string cmd = "set-colour-map-rotation-for-maps";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
}



/*  ----------------------------------------------------------------------- */
/*                         Unit Cell                                        */
/*  ----------------------------------------------------------------------- */
short int
get_show_unit_cell(int imol) {

   std::string cmd = "get-show-unit-cell";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return graphics_info_t::molecules[imol].show_unit_cell_flag;
}

void
set_show_unit_cell(int imol, short int state) {


   //    for (int imol=0; imol<graphics_n_molecules(); imol++) {
   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].show_unit_cell_flag = state;
   }
   //    }
   graphics_draw();
   std::string cmd = "set-show-unit-cell";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(state);
   add_to_history_typed(cmd, args);
}

void set_show_unit_cells_all(short int istate) {

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) { 
	 graphics_info_t::molecules[imol].show_unit_cell_flag = istate;
      }
      if (is_valid_map_molecule(imol)) { 
	 graphics_info_t::molecules[imol].show_unit_cell_flag = istate;
      }
   }
   graphics_draw();
   std::string cmd = "set-show-unit-cells-all";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
} 



// -----------------------------------------------------------------------
//                       Anisotropic Atoms
// -----------------------------------------------------------------------

void 
set_limit_aniso(short int state) {
   //
   graphics_info_t g; 
   g.show_aniso_atoms_radius_flag = state;
   std::string cmd = "set-limit-aniso";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
} 

void
set_aniso_limit_size_from_widget(const char *text) {
   
   float tmp;
   graphics_info_t g; 

   tmp = atof(text);

   if ((tmp >= 0.0) && (tmp < 99999.9)) {

      g.show_aniso_atoms_radius = tmp;
   } else {
      cout << "Cannot interpret " << text << ".  Assuming 10A" << endl;
      g.show_aniso_atoms_radius = 10.0;
   }
}

float
get_limit_aniso() {

   graphics_info_t g;
   add_to_history_simple("get-limit-aniso");
   return g.show_aniso_atoms_radius;

}

// Do if for all molecule, if we do it for one.
// The Anisotropic Atoms widget has no ability to select molecules.
// 
// This then is not a property of a molecule, but is a property of the
// graphics. 
// 
short int
get_show_aniso() {

   return graphics_info_t::show_aniso_atoms_flag;
}

void
set_show_aniso(int state) {

   graphics_info_t::show_aniso_atoms_flag = state;
   graphics_draw();
}

char *get_text_for_aniso_limit_radius_entry() {
   char *text;
   graphics_info_t g;
   
   text = (char *) malloc(100);
   snprintf(text, 99, "%-5.1f", g.show_aniso_atoms_radius);

   return text;
}

short int
get_show_limit_aniso() {
   
   return graphics_info_t::show_aniso_atoms_radius_flag;
}

void
set_aniso_probability(float f) {

   graphics_info_t::show_aniso_atoms_probability = f;
   graphics_draw();
}

// return e.g. 47.9 (%). 
float
get_aniso_probability() {

   return graphics_info_t::show_aniso_atoms_probability;
}


/*  ---------------------------------------------------------------------- */
/*                         Display Functions                               */
/*  ---------------------------------------------------------------------- */

void set_default_bond_thickness(int t) {

   graphics_info_t g;
   g.default_bond_width = t;

} 

void set_bond_thickness(int imol, float t) {

   graphics_info_t g;
   g.set_bond_thickness(imol, t);

}

void set_bond_thickness_intermediate_atoms(float t) { 

   graphics_info_t g;
   g.set_bond_thickness_intermediate_atoms(t);

} 


void set_unbonded_atom_star_size(float f) {
   graphics_info_t g;
   g.unbonded_atom_star_size = f;
} 

int make_ball_and_stick(int imol, const char *atom_selection_str,
			float bond_thickness, float sphere_size,
			int do_spheres_flag) {

   int i = imol;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].make_ball_and_stick(std::string(atom_selection_str),
							   bond_thickness,
							   sphere_size, do_spheres_flag);
      graphics_draw();
   }
   return i;
}


int clear_ball_and_stick(int imol) {

   if (is_valid_model_molecule(imol)) {
      GLuint dummy_tag;
      graphics_info_t::molecules[imol].clear_display_list_object(dummy_tag);
      graphics_draw();
   }
   return 0;
} 

/*  ----------------------------------------------------------------------- */
/*                  dots display                                            */
/*  ----------------------------------------------------------------------- */
int dots(int imol,
	  const char *atom_selection_str,
	  float dot_density, float sphere_size_scale) {

   int idots = -1;
   if (is_valid_model_molecule(imol)) {
      if (atom_selection_str) { 
	 idots = graphics_info_t::molecules[imol].make_dots(std::string(atom_selection_str),
							    dot_density,
							    sphere_size_scale);
      }
   }
   graphics_draw();
   return idots;
}

void clear_dots(int imol, int dots_handle) {

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules)) { 
      graphics_info_t::molecules[imol].clear_dots(dots_handle);
      graphics_draw();
   }
}

/* return the number of dots sets for molecule number imol */
int n_dots_sets(int imol) {

   int r = -1;

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules)) { 
      r = graphics_info_t::molecules[imol].n_dots_sets();
   } else {
      std::cout << "WARNING:: Bad molecule number: " << imol << std::endl;
   } 
   return r;
} 


std::pair<short int, float> float_from_entry(GtkWidget *entry) {

   std::pair<short int, float> p(0,0);
   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      float f = atof(txt);
      p.second = f;
      p.first = 1;
   }
   return p;
}

std::pair<short int, int> int_from_entry(GtkWidget *entry) {

   std::pair<short int, int> p(0,0);
   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      int i = atoi(txt);
      p.second = i;
      p.first = 1;
   }
   return p;
}




// -----------------------------------------------------------------------
//                       Smooth Scrolling
// -----------------------------------------------------------------------

void set_smooth_scroll_flag(int v) {

   graphics_info_t::smooth_scroll = v;
}

int  get_smooth_scroll() {
   
   return graphics_info_t::smooth_scroll;
}

// useful interface for gui (entry)
void set_smooth_scroll_steps_str(const char *text) {

   int v;
   v = atoi(text);
   if (v > 0 && v < 10000000) {
      set_smooth_scroll_steps(v);
   } else {
      cout << "Cannot interpret " << text << ".  Assuming 10 steps" << endl;
      set_smooth_scroll_steps(10);
   }
}

// useful interface for scripting
void set_smooth_scroll_steps(int v) {
      graphics_info_t::smooth_scroll_steps = v;
} 

   
char *get_text_for_smooth_scroll_steps() {

   char *text;

   text = (char *) malloc(100);
   snprintf(text, 99, "%-5d", graphics_info_t::smooth_scroll_steps);

   return text;
}

// useful interface for gui (entry)
void  set_smooth_scroll_limit_str(const char *text) {

   float v;

   v = atof(text);

   if (v >0 && v < 1000) { 
      graphics_info_t::smooth_scroll_limit = v;
   } else {
      cout << text << " out of range: using 10A" << endl;
      graphics_info_t::smooth_scroll_limit = 10;
   }
}

// useful for scripting
void set_smooth_scroll_limit(float lim) {
   graphics_info_t::smooth_scroll_limit = lim; 
} 

char *get_text_for_smooth_scroll_limit() {

   char *text;
   
   text = (char *) malloc(100);
   snprintf(text, 99, "%-5.1f", graphics_info_t::smooth_scroll_limit);

   return text;
}

void set_stop_scroll_diff_map(int i) { 
   graphics_info_t::stop_scroll_diff_map_flag = i;
} 

void set_stop_scroll_iso_map(int i) { 
   graphics_info_t::stop_scroll_iso_map_flag = i;
}


void set_stop_scroll_diff_map_level(float f) { 
   graphics_info_t::stop_scroll_diff_map_level = f;
} 

void set_stop_scroll_iso_map_level(float f) { 
   graphics_info_t::stop_scroll_iso_map_level = f;
}



// -----------------------------------------------------------

void set_font_size(int i) { 

   graphics_info_t g;

   g.set_font_size(i);

}

int get_font_size() {
   return graphics_info_t::atom_label_font_size;
}

/*! \brief set the colour of the atom label font - the arguments are
  in the range 0->1 */
void set_font_colour(float red, float green, float blue) {
   graphics_info_t::font_colour = coot::colour_holder(red, green, blue);
   graphics_draw();
}


/*  ---------------------------------------------------------------------- */
/*                         Rotation Centre Cube Size                       */
/*  ---------------------------------------------------------------------- */

void set_rotation_centre_size_from_widget(const gchar *text) { 
   
   float val;
   graphics_info_t g;

   val = atof(text); 
   if ((val > 1000) || (val < 0)) { 
      cout << "Invalid cube size: " << text << ". Assuming 1.0A" << endl; 
      val = 1.0; 
   } 
   g.rotation_centre_cube_size = val; 
   graphics_draw();
}

void set_rotation_centre_size(float f) {
   graphics_info_t g;
   g.rotation_centre_cube_size = f;
   graphics_draw();
} 

gchar *get_text_for_rotation_centre_cube_size() { 
   
   char *text; 
   graphics_info_t g; 

   text = (char *)  malloc (100);
   snprintf(text, 90, "%-6.3f", g.rotation_centre_cube_size); 
   return text; 
}

short int
recentre_on_read_pdb() {
   return graphics_info_t::recentre_on_read_pdb; 
}

void
set_recentre_on_read_pdb(short int i) {
   graphics_info_t::recentre_on_read_pdb = i;
}

/*  ---------------------------------------------------------------------- */
/*                         orthogonal axes                                 */
/*  ---------------------------------------------------------------------- */
void set_draw_axes(int i) {
   graphics_info_t::draw_axes_flag = i;
} 


// -------------------------------------------------------------------------
//                        (density) iso level increment entry
// -------------------------------------------------------------------------
//

// imol is ignored.
//
char* get_text_for_iso_level_increment_entry(int imol) {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-6.4f", g.iso_level_increment);

   return text;

}

void set_iso_level_increment_from_text(const char *text, int imol) {

   float val;

   graphics_info_t g;

   val = atof(text);

   if ((val > 10000) || (val < -10000)) {
      cout << "Cannot interpret: " << text
	   << ".  Assuming 0.05 for increment" << endl;
      val  = 0.05;

   }

   cout << "setting iso_level_increment to " << val << endl; 
   g.iso_level_increment = val;

   graphics_draw();
}

void set_iso_level_increment(float val) { 
   graphics_info_t g;
   g.iso_level_increment = val;
} 

// imol is ignored.
//
char* get_text_for_diff_map_iso_level_increment_entry(int imol) {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-6.4f", g.diff_map_iso_level_increment);
   return text;

}

void set_diff_map_iso_level_increment_from_text(const char *text, int imol) {

   float val;
   graphics_info_t g;

   val = atof(text);

   if ((val > 10000) || (val < -10000)) {
      cout << "Cannot interpret: " << text
	   << ".  Assuming 0.005 for increment" << endl;
      val  = 0.005;
   } 
   g.diff_map_iso_level_increment = val;
   graphics_draw();
}

void set_diff_map_iso_level_increment(float val) { 
   graphics_info_t::diff_map_iso_level_increment = val;
} 

void set_map_sampling_rate_text(const char *text) {

   float val;
   val = atof(text);

   if ((val > 100) || (val < 1)) {
      cout << "Nonsense value: " << text
	   << ".  Assuming 1.5 for increment" << endl;
      val  = 1.5;
   }
   set_map_sampling_rate(val);

}

void set_map_sampling_rate(float r) {

   graphics_info_t g;
   g.map_sampling_rate = r;

}

char* get_text_for_map_sampling_rate_text() {

   char *text;
   graphics_info_t g;

   text = (char *) malloc (100);
   snprintf(text, 90, "%-5.4f", g.map_sampling_rate);
   return text;


}

float get_map_sampling_rate() {
   graphics_info_t g;
   return g.map_sampling_rate;
}


/* applies to the current map */
void change_contour_level(short int is_increment) { // else is decrement. 

   graphics_info_t g; 
   int s = g.scroll_wheel_map;

   if (g.molecules[s].max_xmaps > 0) {
      if (g.molecules[s].is_difference_map_p()) {
	 g.molecules[s].contour_level[0] +=
	    g.diff_map_iso_level_increment;
      } else {
	 // normal case
	 if (is_increment) { 
	    g.molecules[s].contour_level[0] += g.iso_level_increment;
	 } else {
	    g.molecules[s].contour_level[0] -= g.iso_level_increment;
	 }
      }
      g.molecules[s].update_map();
      graphics_draw();
      std::cout << "contour level of molecule [" << s << "]:  "
		<< g.molecules[s].contour_level[0] << std::endl;
   }
} 

GtkWidget *main_window() {
   return graphics_info_t::glarea; 
}; 

int graphics_n_molecules() {
   return graphics_info_t::n_molecules;
}

int next_map_for_molecule(int imol) { /* return a map number */
   return graphics_info_t::molecules[imol].next_free_map();
}

// imol is used imap is ignored.
// You can fix this anachronism one day if you like.  FIXME.
// 
void set_scrollable_map(int imol) {

   graphics_info_t g;
   if (is_valid_map_molecule(imol)) {
      int imap = 0; // ignored
      g.set_Scrollable_Map(imol, imap); // in graphics-info.h
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid molecule"
		<< " in set_scrollable_map\n";
   }
     
}
 
/*  ----------------------------------------------------------------------- */
/*                  utility function                                        */
/*  ----------------------------------------------------------------------- */
// return -1 if atom not found.
int atom_index(int imol, const char *chain_id, int iresno, const char *atom_id) {

   int index = -1;
   graphics_info_t g;
   std::string altconf("");
   std::string inscode("");

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules) { 
	 // return g.molecules[imol].atom_index(chain_id, iresno, atom_id);
	 index = g.molecules[imol].full_atom_spec_to_atom_index(std::string(chain_id),
								iresno,
								inscode,
								std::string(atom_id),
								altconf);
      }
   }

   return index;
}

// Refine zone needs to be passed atom indexes (which it then converts
// to residue numbers - sigh).  So we need a function to get an atom
// index from a given residue to use with refine_zone()
// 
int atom_index_first_atom_in_residue(int imol, const char *chain_id, 
				     int iresno, const char *ins_code) {

   int index = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g; 
      index = g.molecules[imol].atom_index_first_atom_in_residue(std::string(chain_id),
								 iresno,
								 std::string(ins_code));
   }
   return index;
} 

float median_temperature_factor(int imol) {

   float low_cut = 2.0;
   float high_cut = 100.0;
   short int low_cut_flag = 0;
   short int high_cut_flag = 0;

   float median = -1.0;
   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 median = coot::util::median_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
							graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
							low_cut, high_cut,
							low_cut_flag,
							high_cut_flag);
      } else {
	 std::cout << "WARNING:: molecule " << imol << " has no model\n";
      }
   } else {
      std::cout << "WARNING:: no such molecule as " << imol << "\n";
   }
   return median;
}

float average_temperature_factor(int imol) { 

   float low_cut = 2.0;
   float high_cut = 100.0;
   short int low_cut_flag = 0;
   short int high_cut_flag = 0;

   float av = -1.0;
   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 av = coot::util::average_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
						     graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
						     low_cut, high_cut,
						     low_cut_flag,
						     high_cut_flag);
      } else {
	 std::cout << "WARNING:: molecule " << imol << " has no model\n";
      }
   } else {
      std::cout << "WARNING:: no such molecule as " << imol << "\n";
   }
   return av;
}

char *centre_of_mass_string(int imol) {

   char *s = 0; // guile/SWIG sees this as #f
   if (is_valid_model_molecule(imol)) {
      realtype x, y, z;
      GetMassCenter(graphics_info_t::molecules[imol].atom_sel.atom_selection,
		    graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
		    x, y, z);
      std::string sc = "(";
      sc += graphics_info_t::float_to_string(x);
      sc += " ";
      sc += graphics_info_t::float_to_string(y);
      sc += " ";
      sc += graphics_info_t::float_to_string(z);
      sc += ")";
      s = new char [sc.length() + 1];
      strcpy(s, sc.c_str());
      return s;
   }
   return s;
}



void clear_pending_picks() {
   graphics_info_t g;
   g.clear_pending_picks();
}

#include "gl-matrix.h"
void print_view_matrix() { 		/* print the view matrix */

   graphics_info_t g;
   GL_matrix m;
   m.from_quaternion(g.quat);
   std::cout << "View Matrix:" << std::endl;
   m.print_matrix();
}

float get_view_matrix_element(int row, int col) {

   graphics_info_t g;
   GL_matrix m;
   m.from_quaternion(g.quat);
   return m.matrix_element(row, col);
}


float get_view_quaternion_internal(int element) {

   if ((element >= 0) &&
       (element < 4)) {
      return graphics_info_t::quat[element];
   } else {
      std::cout << "Bad element for quaternion: " << element
		<< " returning dummy -9999" << std::endl;
      return -9999;
   }
}

void set_view_quaternion(float i, float j, float k, float l) {

   double mag2 = i*i + j*j + k*k + l*l;
   double mag=sqrt(mag2);

   if (fabs(mag) > 0.5) {
      graphics_info_t::quat[0] = i/mag;
      graphics_info_t::quat[1] = j/mag;
      graphics_info_t::quat[2] = k/mag;
      graphics_info_t::quat[3] = l/mag;
      graphics_draw();
   } else {
      std::cout << "Bad view quaternion" << std::endl;
   } 
}



// ------------------------------------------------------
//                   Skeleton
// ------------------------------------------------------


void
handle_skeleton_colour_change(int mol, gdouble* map_col) {

   graphics_info_t::skeleton_colour[0] = map_col[0];
   graphics_info_t::skeleton_colour[1] = map_col[1];
   graphics_info_t::skeleton_colour[2] = map_col[2];

   graphics_draw();

}

gdouble*
get_skeleton_colour() {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   colour[0] = graphics_info_t::skeleton_colour[0]; 
   colour[1] = graphics_info_t::skeleton_colour[1]; 
   colour[2] = graphics_info_t::skeleton_colour[2]; 

   return colour;
}

void set_skeleton_colour(int imol, float r, float g, float b) {

   graphics_info_t::skeleton_colour[0] = r;
   graphics_info_t::skeleton_colour[1] = g;
   graphics_info_t::skeleton_colour[2] = b;

   graphics_draw();
}

void
skel_greer_on() {

   int i_skel_set = 0;
   graphics_info_t g; 
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    g.molecules[imol].greer_skeleton_draw_on = 1;
	    // g.molecules[imol].update_skeleton(); // withdrawn
	    i_skel_set = 1;
	    break;
	 }
      }
      if (i_skel_set) break;
   }
   graphics_draw();
}

void
skel_greer_off() {

   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      for (int imap=0; imap<graphics_info_t::molecules[imol].max_xmaps; imap++) {

	 if (graphics_info_t::molecules[imol].xmap_is_filled[imap] &&
	     graphics_info_t::molecules[imol].xmap_is_diff_map[imap] != 1) {

	    graphics_info_t::molecules[imol].greer_skeleton_draw_on = 0;
	 }
      }
   }
}

#include "graphical_skel.h"
#include "xmap-utils.h"

// For some as yet unknown reason, this code is executed when we
// select skeleton off.
// (After it has run, skel_foadi_off is executed).
// 
void
skel_foadi_on() {

   int i_found_skeletonizable_map = 0;
   graphics_info_t g;

   // we use break so that we get only one skeleton. 
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    i_found_skeletonizable_map = 1;

	    // so that we don't do this when the skeleton is on already:
	    //
	    if (g.molecules[imol].fc_skeleton_draw_on == 0) {
	       g.molecules[imol].fc_skeleton_draw_on = 1;

	       mean_and_variance<float> mv = 
		  map_density_distribution(g.molecules[imol].xmap_list[imap],0); 

	       cout << "Mean and sigma of map: " << mv.mean 
		    << " and " << sqrt(mv.variance) << endl; 

	       float map_cutoff = mv.mean + 1.5*sqrt(mv.variance); 
	       g.skeleton_level = map_cutoff; 
	    
	       // derived from sktest:
	       // 
	       g.molecules[imol].xskel_cowtan.init(g.molecules[imol].xmap_list[imap].spacegroup(), 
						   g.molecules[imol].xmap_list[imap].cell(),
						   g.molecules[imol].xmap_list[imap].grid_sampling());

	       cout << "INFO:: making skeleton cowtan..." << endl; 
	       GraphicalSkel cowtan(g.molecules[imol].xmap_list[0],
				    g.molecules[imol].xskel_cowtan); //fill xskel_cowtan

	       g.molecules[imol].xskel_is_filled = 1; // TRUE

	       // various experiments....

	       // cowtan.tip_filter(xmap_list[0], &xskl); // tinker with xskel_cowtan

	       //cowtan.prune(g.molecules[imol].xmap_list[imap],
	       //	 &g.molecules[imol].xskel_cowtan);

	       //
	       cout << "INFO:: pruning cowtan..." << endl; 
	       //
	       cowtan.Pprune(g.molecules[imol].xmap_list[imap],
			     &g.molecules[imol].xskel_cowtan,
			     map_cutoff);


	       // --------------------------------------------------------------
	       // --------------------------------------------------------------
	       // 	    cout << "cuckoo code:" << endl; 

	       // 	    BuildCas bc(g.molecules[imol].xmap_list[imap], map_cutoff); 

	       // 	    // mark segments by connectivity
	       // 	    // 
	       // 	    int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
	       // 						       g.molecules[imol].xmap_list[imap],
	       // 						       map_cutoff); 

	       // 	    cout << "There were " << nsegments << "segment in cuckoo map" << endl; 
	       // 	    // now transfer the segment map to xskel_cowtan
	       // 	    // 
	       // 	    bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	    


	       // --------------------------------------------------------------
	       // --------------------------------------------------------------


	       // now display the skeleton
	    
	       g.molecules[imol].update_clipper_skeleton();

	       // now create a new molecule and put the bonds
	       // into it
	       // 
	       // g.n_molecules++;
	       // 	    std::cout << "g.n_molecules is now " << g.n_molecules << endl;

	       // 
	       break;   // only do one (the first one). 
	    }
	 }
      }
      if (i_found_skeletonizable_map) break;
   }
   graphics_draw();
}

void
skel_foadi_off() {

   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      for (int imap=0; imap<graphics_info_t::molecules[imol].max_xmaps; imap++) {
	 
	 if (graphics_info_t::molecules[imol].xmap_is_filled[imap] &&
	     graphics_info_t::molecules[imol].xmap_is_diff_map[imap] != 1) {
	    
	    graphics_info_t::molecules[imol].fc_skeleton_draw_on = 0;
	 }
      }
   }
}


int 
skeletonize_map(int prune_flag, int imol) {

   graphics_info_t::skeletonize_map(prune_flag, imol);
   return 0;
}

int unskeletonize_map(int imol) {
   graphics_info_t::unskeletonize_map(imol);
   return imol;
} 




void
do_skeleton_prune() { 

   graphics_info_t g;
   float map_cutoff  = g.skeleton_level;

   // we use break so that we get only one skeleton. 
   short int done_skel = 0;
   
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {

	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {
	    
	    if (g.molecules[imol].xskel_is_filled == 1) { 
	       
	       BuildCas bc(g.molecules[imol].xmap_list[imap], map_cutoff); 


	       // mark segments by connectivity
	       // 
	       int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan, 
							  g.molecules[imol].xmap_list[imap],
							  map_cutoff); 

	       cout << "INFO:: There were " << nsegments << " different segments" << endl; 
	       
	       bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
	       g.molecules[imol].update_clipper_skeleton();
	       done_skel = 1;
	       break;
	    }
	 }
      }
      if (done_skel) break;
   }
} 

void test_fragment() {

#ifdef HAVE_GSL
   graphics_info_t g;
   g.rotamer_graphs(0);

#endif // HAVE_GSL
}

// we redefine TRUE here somewhere...
// #include <gdk/gdkglconfig.h>
// #include <gtk/gtkgl.h>
// #include <gdk/x11/gdkglx.h>
// #include <gdk/x11/gdkglglxext.h>

int test_function(int i, int j) {

//    graphics_info_t g;
//    g.wrapped_create_symmetry_controller_dialog();
//    return 0;

   if (1) {
      graphics_info_t g;
      g.Geom_p()->hydrogens_connect_file("THH", "thh_connect.txt");
   }

   if (0) {

      // GTK2 GTkGLExt code
//       GdkGLContext *glcontext = gtk_widget_get_gl_context(graphics_info_t::glarea);
//       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(graphics_info_t::glarea);
//       GdkGLConfig *glconfig = gtk_widget_get_gl_config(graphics_info_t::glarea);
//       Display *dpy = gdk_x11_gl_config_get_xdisplay (glconfig);
      // Bool glXMakeCurrent(Display * dpy,
      //                     GLXDrawable  Drawable,
      //                     GLXContext  Context)
      // gdk_gl_glXMakeContextCurrent(dpy, gldrawable, glcontext);

      // bwah!
      // glXMakeCurrent(dpy, gldrawable, glcontext);

      // another way?
//       GtkWidget *w = graphics_info_t::glarea;
//       GdkGLContext *glcontext = gtk_widget_get_gl_context (w);
//       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (w);
//       int i = gdk_gl_drawable_gl_begin (gldrawable, glcontext);
//       std::cout << "DEBUG gdk_gl_drawable_gl_begin returns state: "
// 		<< i << std::endl;
//       return i;
   } 

   if (0) {
      int imol = i;
      if (is_valid_model_molecule(imol)) { 
	 const coot::residue_spec_t clicked_residue("A", 1);
	 short int is_n_term_addition = 1;
	 CAtom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[10];
	 CChain *chain_p = at->GetChain();
	 std::pair<bool, std::string> p = 
	    graphics_info_t::molecules[imol].residue_type_next_residue_by_alignment(clicked_residue, chain_p, is_n_term_addition);
	 if (p.first == 1) { 
	    std::cout << "next residue: " << p.second << std::endl;
	 } else {
	    std::cout << "no next residue found." << std::endl;
	 }
      }
   } 


   if (0) { 
      GtkWidget *w = wrapped_create_least_squares_dialog();
      gtk_widget_show(w);
   }
      

   if (0) { 
      std::vector<std::string> s;
      s.push_back("");
      s.push_back("123");
      s.push_back("123/456");
      s.push_back("123/456/");
      
      for (unsigned int i=0; i<s.size(); i++) { 
	 std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(s[i]);
	 std::cout << "For string :" << s[i] << ": split is :"
		   << p.first << ": :" << p.second << ":" << std::endl; 
      }

      std::string t = "/my/thing/int.mtz data/crystal/FWT data/crystal/PHWT";
      std::vector<std::string> v = coot::util::split_string(t, " ");

      std::cout << "splitting :" << t << ": on " << " " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
	 std::cout << "split " << i << " :" << v[i] << ":\n";
      }
   }
   return 0;
}

int write_connectivity(const char *monomer_name, const char *filename) {

   graphics_info_t g;
   return g.Geom_p()->hydrogens_connect_file(monomer_name, filename);
} 


void screendump_image(const char *filename) {

   int istatus = graphics_info_t::screendump_image(filename);
   std::cout << "screendump_image status " << istatus << std::endl;
   if (istatus) {
      std::string s = "Screendump image ";
      s += filename;
      s += " written";
      graphics_info_t g;
      g.statusbar_text(s);
   }
}

void make_image_raster3d(const char *filename) {

   std::string r3d_name = filename;
   r3d_name += ".r3d";
   raster3d(r3d_name.c_str());
#ifdef USE_GUILE

   std::string cmd("(raytrace 'raster3d ");
   cmd += single_quote(r3d_name);
   cmd += " ";
   cmd += single_quote(filename);
   cmd += "'dummy 'dummy)";
   safe_scheme_command(cmd);
   
#endif   
}

void make_image_povray(const char *filename) {
   std::string pov_name = filename;
   pov_name += ".pov";
   povray(pov_name.c_str());
#ifdef USE_GUILE

   int x_size = graphics_info_t::glarea->allocation.width;
   int y_size = graphics_info_t::glarea->allocation.height;
   std::string cmd("(raytrace 'povray ");
   cmd += single_quote(pov_name);
   cmd += " ";
   cmd += single_quote(filename);
   cmd += " ";
   cmd += graphics_info_t::int_to_string(x_size);
   cmd += " ";
   cmd += graphics_info_t::int_to_string(y_size);
   cmd += ")";
   safe_scheme_command(cmd);

#endif   
}



void
autobuild_ca_off() { 

   graphics_info_t g; 
   g.autobuild_flag = 0; 

}

void
handle_read_ccp4_map(const char* filename, int is_diff_map_flag) {

   if (filename) { 
      std::string str(filename); 
      graphics_info_t g;

      int istate = g.molecules[g.n_molecules].read_ccp4_map(str, is_diff_map_flag,
							    *graphics_info_t::map_glob_extensions);

      if (istate > -1) { // not a failure
	 // std::cout << "successfully read map into molecule " << g.n_molecules << std::endl;
	 int imol = graphics_n_molecules();
	 std::string name = g.molecules[imol].dotted_chopped_name();
	 g.scroll_wheel_map = imol;  // change the current scrollable map.
	 g.n_molecules++; 
      } else { 
	 std::cout << "Read map " << str << " failed" << std::endl;
	 std::string s = "Read map ";
	 s += str;
	 s += " failed.";
	 g.statusbar_text(s);
      } 
      graphics_draw();
   } else {
      // error
      std::cout << "ERROR:: filename null in handle_read_ccp4_map\n";
   } 
}


/*  ------------------------------------------------------------------------ */
/*                         clipping */
/*  ------------------------------------------------------------------------ */

void do_clipping1_activate(){

   GtkWidget *clipping_window;
   GtkScale *hscale;
   GtkAdjustment *adjustment;
 
   /* connect this to displaying the new clipping window */

   clipping_window = create_clipping_window();

   hscale = GTK_SCALE(lookup_widget(clipping_window, "hscale1")); 
/*    gtk_scale_set_draw_value(hscale, TRUE);  already does */

   adjustment = GTK_ADJUSTMENT 
      (gtk_adjustment_new(0.0, -10.0, 20.0, 0.05, 4.0, 10.1)); 

   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   gtk_signal_connect (GTK_OBJECT (adjustment), "value_changed",
		       GTK_SIGNAL_FUNC (clipping_adjustment_changed), NULL);
   
   gtk_widget_show(clipping_window); 
}

void clipping_adjustment_changed (GtkAdjustment *adj, GtkWidget *window) { 

   /*    printf("Clipping adjustment: %f\n", adj->value); */

   set_clipping_front(adj->value);
   set_clipping_back (adj->value);
}
 


/*  ----------------------------------------------------------------------- */
/*                        virtual trackball                                 */
/*  ----------------------------------------------------------------------- */

void
vt_surface(int v){ 

   graphics_info_t g;
   g.set_vt_surface(v);
   std::vector<std::string> command_strings;
//    command_strings.push_back("vt-surface");
//    command_strings.push_back(graphics_info_t::int_to_string(v));
//    add_to_history(command_strings);
}

int vt_surface_status() { 

   graphics_info_t g;
   return g.vt_surface_status();
}

/*  ----------------------------------------------------------------------- */
/*                        save coordintes                                   */
/*  ----------------------------------------------------------------------- */


// return status 1 is good, 0 is fail.
int save_coordinates(int imol, const char *filename) { 

   int ierr = 0;
   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    ierr = graphics_info_t::molecules[imol].save_coordinates(filename);
	 }
      } 
   } 
   std::vector<std::string> command_strings;
   command_strings.push_back("save-coordinates");
   command_strings.push_back(coot::util::int_to_string(imol));
   command_strings.push_back(filename);
   add_to_history(command_strings);
   return ierr;
}

void set_save_coordinates_in_original_directory(int i) {

   // not used now, I think.
   graphics_info_t::save_coordinates_in_original_dir_flag = i;

}


/* access to graphics_info_t::save_imol for use in callback.c */
int save_molecule_number_from_option_menu() {
   return graphics_info_t::save_imol;
}

/* access from callback.c, not to be used in scripting, I suggest. */
void set_save_molecule_number(int imol) {
   graphics_info_t::save_imol = imol;
}



/*  ----------------------------------------------------------------------- */
/*                        .phs file reading                                 */
/*  ----------------------------------------------------------------------- */

void
read_phs_and_coords_and_make_map(const gchar *pdb_filename){

   // This function is the .phs equivalent of c.f. make_and_draw_map,
   // map_fill_from_mtz.  We have previously stored the phs_filename
   // in the static graphics_info_t.
   // 
   graphics_info_t g; 

   int imol = g.n_molecules;

   // don't forget that this is a map.
   //
   int istat = g.molecules[imol].make_map_from_phs(std::string(pdb_filename),
						   g.get_phs_filename());

   if (istat != -1) { 
      g.n_molecules++;
      graphics_draw();
   } else {
      // give us a warning message then
      std::string w = "Sadly, the cell or space group is not comprehensible in\n";
      w += "the pdb file: ";
      w += pdb_filename;
      w += "\n";
      w += "Can't make map from phs file.";
      GtkWidget *widget = wrapped_nothing_bad_dialog(w);
      gtk_widget_show(widget);
   } 
}

/*! \brief read a phs file, the cell and symm information is from
  previously read (most recently read) coordinates file

 For use with phs data filename provided on the command line */
int 
read_phs_and_make_map_using_cell_symm_from_previous_mol(const char *phs_filename) {

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int done_flag = 0;
   int r = -1;

   int imol_ref = -1;

   for (int i=graphics_info_t::n_molecules-1; i>=0; i--) {
      if (is_valid_model_molecule(i)) {
	 imol_ref = i;
	 break;
      }
   }

   if (imol_ref > -1) 
      r = read_phs_and_make_map_using_cell_symm_from_mol(phs_filename, imol_ref);

   return r;
}


/*! \brief read a phs file and use the cell and symm in molecule
  number imol and use the resolution limits reso_lim_low and
  reso_lim_high  */
int
read_phs_and_make_map_with_reso_limits(int imol_ref, const char* phs_filename,
				       float reso_lim_low, float reso_lim_high) {
   // This function is the .phs equivalent of c.f. make_and_draw_map,
   // map_fill_from_mtz.  We have previously stored the phs_filename
   // in the static graphics_info_t.
   // 
   graphics_info_t g;
   int imol = g.n_molecules;

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int done_flag = 0;
   int istat = -1; // returned value


   if (g.molecules[imol_ref].have_unit_cell) {
      // convert from a set of coordinates

      spacegroup.init(clipper::Spgr_descr(g.molecules[imol_ref].atom_sel.mol->GetSpaceGroup()));

      clipper::Cell_descr cell_d(g.molecules[imol_ref].atom_sel.mol->get_cell().a,
				 g.molecules[imol_ref].atom_sel.mol->get_cell().b,
				 g.molecules[imol_ref].atom_sel.mol->get_cell().c,
				 clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().alpha),
				 clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().beta),
				 clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().gamma));

      cell.init(cell_d);
      done_flag = 1;
   } else {
      // no conversion needed, just get from map
      if (g.molecules[imol_ref].has_map()) {
	 cell.init(g.molecules[imol_ref].xmap_list[0].cell().descr());
	 spacegroup.init(g.molecules[imol_ref].xmap_list[0].spacegroup().descr()); 
	 done_flag = 1;
      } 
   }

   if (done_flag) { 

      // don't forget that this is a map.
      //
      std::string phs_file(phs_filename);
      istat = g.molecules[imol].make_map_from_phs_using_reso(phs_file,
								 spacegroup,
								 cell, 
								 reso_lim_low, reso_lim_high);

      if (istat != -1) {
	 imol = istat;
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 std::string w = "Sadly, something bad happened reading phs file using\n";
	 w += "the molecule number ";
	 w += coot::util::int_to_string(imol_ref); 
	 w += "\n";
	 w += "Can't make map from phs file.";
	 GtkWidget *widget = wrapped_nothing_bad_dialog(w);
	 gtk_widget_show(widget);
      }
   } else {
      // give us a warning message then
      std::string w = "Sadly, the cell or space group is not comprehensible in\n";
      w += "the molecule number ";
      w += coot::util::int_to_string(imol_ref); 
      w += "\n";
      w += "Can't make map from phs file.";
      GtkWidget *widget = wrapped_nothing_bad_dialog(w);
      gtk_widget_show(widget);
   }

   return istat;
}



int
read_phs_and_make_map_using_cell_symm_from_mol(const char *phs_filename_str, int imol_ref) { 

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int got_cell_symm_flag = 0;
   int imol = -1;// set bad molecule initally
   
   graphics_info_t g; 
//       std::cout << "DEBUG:: read_phs_and_make_map_using_cell_symm_from_mol "
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().a << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().b << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().c << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().alpha << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().beta << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().gamma << "  " 
// 		<< std::endl;

   if (is_valid_model_molecule(imol_ref) || is_valid_map_molecule(imol_ref)) {
      if (g.molecules[imol_ref].have_unit_cell) { 

	 std::string s(g.molecules[imol_ref].atom_sel.mol->GetSpaceGroup());
	 if (s == "R 3") {
	    std::cout << "debug --------- symbol transformation ------\n";
	    s = "P 3*";
	 }
	 std::cout << "------------- :" << s << ":\n";
	 s = "P 3*";
	 clipper::Spgr_descr sgd(s);
	 spacegroup.init(sgd);

	 clipper::Cell_descr cell_d(g.molecules[imol_ref].atom_sel.mol->get_cell().a,
				    g.molecules[imol_ref].atom_sel.mol->get_cell().b,
				    g.molecules[imol_ref].atom_sel.mol->get_cell().c,
				    clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().alpha),
				    clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().beta),
				    clipper::Util::d2rad(g.molecules[imol_ref].atom_sel.mol->get_cell().gamma));

	 cell.init(cell_d);
	 got_cell_symm_flag = 1;
      } else {
	 // get the cell/symm from a map:
	 if (g.molecules[imol_ref].has_map()) {
	    cell.init(g.molecules[imol_ref].xmap_list[0].cell().descr());
	    spacegroup.init(g.molecules[imol_ref].xmap_list[0].spacegroup().descr()); 
	    got_cell_symm_flag = 1;
	 } else { 

	 } 
      }

      if (got_cell_symm_flag) {
	 std::string phs_filename(phs_filename_str); 

	 imol = g.n_molecules; 
	 g.molecules[g.n_molecules].make_map_from_phs(spacegroup, cell, phs_filename);
	 g.n_molecules++; 
	 graphics_draw();
      }
   }

   return imol;
}


int
read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(int imol_ref) { 

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   short int done_flag = 0;
   int imol = -1; // bad molecule
   
   graphics_info_t g; 
//       std::cout << "DEBUG:: read_phs_and_make_map_using_cell_symm_from_mol "
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().a << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().b << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().c << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().alpha << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().beta << "  " 
// 		<< g.molecules[imol_ref].atom_sel.mol->get_cell().gamma << "  " 
// 		<< std::endl;

   if (is_valid_model_molecule(imol_ref) || is_valid_map_molecule(imol_ref)) { 

      if (g.molecules[imol_ref].have_unit_cell) {
	 // Coordinates block:
	 MyCMMDBManager *mol = g.molecules[imol_ref].atom_sel.mol;;
	 clipper::Spacegroup t_spgr = coot::util::get_spacegroup_from_symops(mol);
	 if (t_spgr.is_null()) {
	    std::cout << "WARNING:: Cant get spacegroup from coordinates!\n";
	 } else {
	    spacegroup = t_spgr; 
	    // If you don't want to use get_cell, use mol->GetCell();
	    clipper::Cell_descr cell_d(mol->get_cell().a,
				       mol->get_cell().b,
				       mol->get_cell().c,
				       clipper::Util::d2rad(mol->get_cell().alpha),
				       clipper::Util::d2rad(mol->get_cell().beta),
				       clipper::Util::d2rad(mol->get_cell().gamma));
   
	    cell.init(cell_d);
	    done_flag = 1;
	 }

      } else {
	 // Map block
	 if (g.molecules[imol_ref].has_map()) {
	    cell.init(g.molecules[imol_ref].xmap_list[0].cell().descr());
	    spacegroup.init(g.molecules[imol_ref].xmap_list[0].spacegroup().descr()); 
	    done_flag = 1;
	 }
      }

      if (done_flag) {
	 std::string phs_filename(graphics_get_phs_filename()); 

	 imol = g.n_molecules; 
	 g.molecules[g.n_molecules].make_map_from_phs(spacegroup, cell, phs_filename);
	 g.n_molecules++; 
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: Failed to get cell/symm - skipping.\n";
      } 
   }
   return imol;
}

int
read_phs_and_make_map_using_cell_symm(const char *phs_file_name,
				      const char *hm_spacegroup, float a, float b, float c,
				      float alpha, float beta, float gamma) { /*! in degrees */

   clipper::Spacegroup spacegroup; 
   clipper::Cell cell;
   graphics_info_t g;

   spacegroup.init(clipper::Spgr_descr(std::string(hm_spacegroup)));
   cell.init (clipper::Cell_descr(a, b, c, 
				  clipper::Util::d2rad(alpha), 
				  clipper::Util::d2rad(beta), 
				  clipper::Util::d2rad(gamma)));
   
   
   std::string phs_filename(phs_file_name); 

   int imol = g.n_molecules;
   g.molecules[g.n_molecules].make_map_from_phs(spacegroup, cell, phs_filename);
   g.n_molecules++; 
   graphics_draw();

   return imol;
}


void
graphics_store_phs_filename(const gchar *phs_filename) {

   graphics_info_t g;
   g.set_phs_filename(std::string(phs_filename));
}


const char *
graphics_get_phs_filename() {

   graphics_info_t g;
   return g.get_phs_filename().c_str(); 
}

short int possible_cell_symm_for_phs_file() {

   if (graphics_info_t::n_molecules == 0) { 
      return 0; 
   } else {
      return 1; 
   } 
}

// return a string to each of the cell parameters in molecule imol.
// 
gchar *get_text_for_phs_cell_chooser(int imol, char *field) { 

   // we first look in atomseletion

   graphics_info_t g;
   gchar *retval = NULL;
   retval = (gchar *) malloc(12); 
   int ihave_cell = 0; 
   float cell[6];
   const char *spgrp = NULL; 

   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    if (g.molecules[imol].have_unit_cell) { 

	       ihave_cell = 1; 

	       cell[0] = g.molecules[imol].atom_sel.mol->get_cell().a;
	       cell[1] = g.molecules[imol].atom_sel.mol->get_cell().b;
	       cell[2] = g.molecules[imol].atom_sel.mol->get_cell().c;
	       cell[3] = g.molecules[imol].atom_sel.mol->get_cell().alpha;
	       cell[4] = g.molecules[imol].atom_sel.mol->get_cell().beta;
	       cell[5] = g.molecules[imol].atom_sel.mol->get_cell().gamma;
	       spgrp   = g.molecules[imol].atom_sel.mol->GetSpaceGroup(); 

	    } else { 

	       if (g.molecules[imol].max_xmaps > 0) { 

		  ihave_cell = 1; 

		  clipper::Spacegroup spacegroup = g.molecules[imol].xmap_list[0].spacegroup(); 
		  clipper::Cell       ccell      = g.molecules[imol].xmap_list[0].cell();

		  cell[0] = g.molecules[imol].xmap_list[0].cell().a(); 
		  cell[1] = g.molecules[imol].xmap_list[0].cell().b(); 
		  cell[2] = g.molecules[imol].xmap_list[0].cell().c(); 
		  cell[3] = g.molecules[imol].xmap_list[0].cell().alpha() * RADTODEG; 
		  cell[4] = g.molecules[imol].xmap_list[0].cell().beta()  * RADTODEG; 
		  cell[5] = g.molecules[imol].xmap_list[0].cell().gamma() * RADTODEG; 

		  spgrp = spacegroup.descr().symbol_hm().c_str(); 
	       } 
	    } 
   

	    if (spgrp) { 
	       if ( ! (strcmp(field, "symm") ) ) { 
		  snprintf(retval, 11, "%-s", spgrp);
	       }
	       if ( ! (strcmp(field, "a") ) ) { 
		  snprintf(retval, 11, "%7.3f", cell[0]);   
	       } 
	       if ( ! (strcmp(field, "b") ) ) { 
		  snprintf(retval, 11, "%7.2f",  cell[1]);  
	       } 
	       if ( ! (strcmp(field, "c") ) ) { 
		  snprintf(retval, 11, "%7.2f",  cell[2]);  
	       } 
	       if ( ! (strcmp(field, "alpha") ) ) { 
		  snprintf(retval, 11, "%6.2f",   cell[3]); 
	       } 
	       if ( ! (strcmp(field, "beta") ) ) { 
		  snprintf(retval, 11, "%6.2f",  cell[4]);  
	       } 
	       if ( ! (strcmp(field, "gamma") ) ) { 
		  snprintf(retval, 11, "%6.2f",   cell[5]); 
	       }


	       if (! ihave_cell) { 
		  strcpy(retval, "  -  "); 
	       }
	    }
	 }
      }
   }
   return retval; 
}


/*  ----------------------------------------------------------------------- */
/*                        undo last move                                    */
/*  ----------------------------------------------------------------------- */
void undo_last_move() {

   graphics_info_t g;
   g.undo_last_move(); // does a redraw
} 


/*  ----------------------------------------------------------------------- */
/*                        go to atom widget                                 */
/*  ----------------------------------------------------------------------- */

// return -1 on error
//
int go_to_atom_molecule_optionmenu_active_molecule(GtkWidget *widget) { 

   graphics_info_t g;
   return g.go_to_atom_molecule_optionmenu_active_molecule(widget);
}


void save_go_to_atom_widget(GtkWidget *widget) { /* store in a static */
   graphics_info_t::go_to_atom_window = widget;
}

void unset_go_to_atom_widget() {
   graphics_info_t::go_to_atom_window = NULL;
} 


// not really a button select, its a menu item select
void
save_molecule_coords_button_select(GtkWidget *item, GtkPositionType pos) { 

   // graphics_info_t g;
   // std::cout << "INFO:: Save coords molecule now: " << pos << std::endl;
   graphics_info_t::save_imol = pos;
}


int go_to_atom_molecule_number() {
   graphics_info_t g;
   return g.go_to_atom_molecule();
}

char *go_to_atom_chain_id() {
   graphics_info_t g; 
   gchar *txt = (gchar *)malloc(100);
   strcpy(txt, g.go_to_atom_chain()); 
   return txt; 
}

char *go_to_atom_atom_name() {
   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_atom_name()); 
   return txt; 
}

int go_to_atom_residue_number() {
   graphics_info_t g;
   return g.go_to_atom_residue();
}

char *go_to_atom_ins_code() {
   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_ins_code()); 
   return txt; 
}

char *go_to_atom_alt_conf() {
   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   snprintf(txt, 9, "%s", g.go_to_atom_alt_conf()); 
   return txt; 
}


// Note that t3 is an atom name with (possibly) an altLoc tag (after the comma).
// 
int set_go_to_atom_chain_residue_atom_name(const char *t1, int iresno, const char *t3) {

   graphics_info_t g; 

   // so we need to split t3 if it has a comma
   // 
   std::string t3s(t3);
   std::string::size_type icomma = t3s.find_last_of(",");
   if (icomma == string::npos) {

      // there was no comma, conventional usage:
      g.set_go_to_atom_chain_residue_atom_name(t1, iresno, t3); 

   } else { 

      std::string atname = t3s.substr(0,icomma);
      std::string altloc = t3s.substr(icomma+1, t3s.length());
      g.set_go_to_atom_chain_residue_atom_name(t1, iresno,
					       atname.c_str(),
					       altloc.c_str());

   }

   int success = g.try_centre_from_new_go_to_atom(); 
   if (success) { 
      update_things_on_move_and_redraw();
      CAtom *at = 0; // passed but not used, it seems.
      GtkWidget *window = graphics_info_t::go_to_atom_window;
      if (window)
	 g.update_widget_go_to_atom_values(window, at);
   }
   return success; 
}




int set_go_to_atom_chain_residue_atom_name_strings(const char *t1, const char *t2, const char *t3)
{
   int it2 = atoi(t2); 
   return set_go_to_atom_chain_residue_atom_name(t1, it2, t3); 
}

// FIXME to use altconf.
// 
// int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec) { 
   
//    return set_go_to_atom_chain_residue_atom_name(atom_spec.chain, 
// 						 atom_spec.resno,
// 						 atom_spec.atom_name);

// }


int 
goto_next_atom_maybe_new(GtkWidget *goto_atom_window) { 

//    int it2 = atoi(t2); 
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, +1);

   graphics_info_t g;
   return g.intelligent_next_atom_centring(goto_atom_window);
   
} 


int 
goto_previous_atom_maybe_new(GtkWidget *goto_atom_window) { 

//    int it2 = atoi(t2); 
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, +1);

   graphics_info_t g;
   return g.intelligent_previous_atom_centring(goto_atom_window);
   
} 
int 
goto_prev_atom_maybe(const gchar *t1, const gchar *t2, const gchar *t3,
		     GtkEntry *res_entry) { 

   int it2 = atoi(t2); 
   return goto_near_atom_maybe(t1, it2, t3, res_entry, -1); 
   
} 

int 
goto_near_atom_maybe(const char *t1, int ires, const char *t3,
		     GtkEntry *res_entry, int idiff) { 

   graphics_info_t g; 

   int ires_l = ires + idiff ; // for next residue, or previous.
   
   g.set_go_to_atom_chain_residue_atom_name(t1, ires_l, t3); 
   
   int success = g.try_centre_from_new_go_to_atom(); 

   if (success) {
      char *txt = (char *)malloc(6);
      snprintf(txt, 5, "%d", ires_l); 
      gtk_entry_set_text(GTK_ENTRY(res_entry), txt);
      update_things_on_move_and_redraw(); 
   }
   return success; 
}


/* For dynarama callback sake. The widget/class knows which coot
   molecule that it was generated from, so in order to go to the
   molecule from dynarama, we first need to the the molecule - because
   set_go_to_atom_chain_residue_atom_name() does not mention the
   molecule (see "Next/Previous Residue" for reasons for that).  This
   function simply calls the graphics_info_t function of the same
   name. */
void set_go_to_atom_molecule(int imol) {

   graphics_info_t g;
   g.set_go_to_atom_molecule(imol); 
   std::vector<std::string> command_strings;
   command_strings.push_back("set-go-to-atom-molecule");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);

}



/*  ----------------------------------------------------------------------- */
/*                  bond representation                                     */
/*  ----------------------------------------------------------------------- */

void graphics_to_ca_representation(int imol) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol))
      g.molecules[imol].ca_representation();
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_ca_representation" << std::endl;
   graphics_draw();
   
   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
} 

void graphics_to_ca_plus_ligands_representation   (int imol) { 
   graphics_info_t g;
   g.molecules[imol].ca_plus_ligands_representation();
   graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-plus-ligands-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}


void graphics_to_bonds_representation(int imol) {
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      g.molecules[imol].bond_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-ca-plus-ligands-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_bonds_representation" << std::endl;
   graphics_draw();

}

void graphics_to_bonds_no_waters_representation(int imol) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)){ 
      g.molecules[imol].bonds_no_waters_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-no-waters-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_bonds_no_waters_representation"
		<< std::endl;
   graphics_draw();

} 

void graphics_to_sec_struct_bonds_representation(int imol) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      g.molecules[imol].bonds_sec_struct_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-sec-struct-bonds-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_sec_struct_bonds_representation"
		<< std::endl;
   graphics_draw();
} 

void graphics_to_ca_plus_ligands_sec_struct_representation(int imol) { 
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      g.molecules[imol].ca_plus_ligands_sec_struct_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-ca-plus-ligands-sec-struct-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_ca_plus_ligands_sec_struct_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_rainbow_representation(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].ca_plus_ligands_rainbow_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-ca-plus-ligands-rainbow-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_ca_plus_ligands_rainbow_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_b_factor_representation(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].b_factor_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("b-factor-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_b_factor_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_occupancy_represenation(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t::molecules[imol].occupancy_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("occupancy-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_occupancy_representation"
		<< std::endl;
   graphics_draw();
}



int
graphics_molecule_bond_type(int imol) { 

   graphics_info_t g;
   // std::cout << "graphics_molecule_bond_type for mol: " << imol << std::endl;
   if (is_valid_model_molecule(imol)) { 
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-molecule-bond-type");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
      return g.molecules[imol].Bonds_box_type();
   }
   return -1;
}

int
set_b_factor_bonds_scale_factor(int imol, float f) {

   int r = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_b_factor_bonds_scale_factor(f);
      r = 1;
   }
   graphics_draw();
   return r;
}



// -------------------------------------------------------------------------
//                        skeletonization level
// -------------------------------------------------------------------------
//

gchar *get_text_for_skeletonization_level_entry() { 

   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10); 
   snprintf(txt, 9, "%f", g.skeleton_level); 

   return txt; 
} 

void set_skeletonization_level_from_widget(const char *txt) { 

   float tmp; 
   graphics_info_t g;

   tmp = atof(txt); 

   if (tmp > 0.0 &&  tmp < 9999.9) { 
      g.skeleton_level = tmp; 
   } else { 
      
      cout << "Cannot interpret " << txt << " using 0.2 instead" << endl; 
      g.skeleton_level = 0.2; 
   } 

         
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {
	 
	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    // 
	    g.molecules[imol].update_clipper_skeleton();

	 } 
      }
   }
   graphics_draw();
}


gchar *get_text_for_skeleton_box_size_entry() { 

   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10); 

   snprintf(txt, 9, "%f", g.skeleton_box_radius); 
   return txt;
} 

void set_skeleton_box_size_from_widget(const char *txt) { 
   float tmp; 
   graphics_info_t g;

   tmp = atof(txt); 

   if (tmp > 0.0 &&  tmp < 9999.9) { 
      g.skeleton_box_radius = tmp; 
   } else { 
      
      cout << "Cannot interpret " << txt << " using 0.2 instead" << endl; 
      g.skeleton_box_radius = 0.2; 
   }

   set_skeleton_box_size(g.skeleton_box_radius);
}

void set_skeleton_box_size(float f) {

   graphics_info_t g;
   g.skeleton_box_radius = f;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-skeleton-box-size");
   command_strings.push_back(graphics_info_t::float_to_string(f));
   add_to_history(command_strings);
      
   for (int imol=0; imol<g.n_molecules; imol++) {
      for (int imap=0; imap<g.molecules[imol].max_xmaps; imap++) {
	 
	 if (g.molecules[imol].xmap_is_filled[imap] &&
	     g.molecules[imol].xmap_is_diff_map[imap] != 1) {

	    // 
	    g.molecules[imol].update_clipper_skeleton();

	 } 
      }
   }
   graphics_draw();
} 

/*  ----------------------------------------------------------------------- */
/*                  map and molecule control                                */
/*  ----------------------------------------------------------------------- */

void save_display_control_widget_in_graphics(GtkWidget *widget) { 

   graphics_info_t g; 
   g.save_display_control_widget_in_graphics(widget);
}

void 
post_display_control_window() { 

   GtkWidget *widget = wrapped_create_display_control_window();
   gtk_widget_show(widget);
   std::vector<std::string> command_strings;
   command_strings.push_back("post-display-control-window");
   add_to_history(command_strings);
}

GSList **gslist_for_scroll_in_display_manager_p() {

   return &graphics_info_t::gslist_for_scroll_in_display_manager;
}



 
void add_map_display_control_widgets() { 

   graphics_info_t g; 

   for (int ii=0; ii<g.n_molecules; ii++) {
      if (g.molecules[ii].max_xmaps > 0){ 
	 
	 g.molecules[ii].update_map_in_display_control_widget(); 
      } 
   }
}


void add_mol_display_control_widgets() { 

   graphics_info_t g; 
   
   for (int ii=0; ii<g.n_molecules; ii++) {
      if (! (g.molecules[ii].atom_sel.atom_selection == NULL)) { 
	 
	 g.molecules[ii].new_mol_in_display_control_widget(); 
      } 
   } 
} 


void add_map_and_mol_display_control_widgets() { 

   add_mol_display_control_widgets(); 
   add_map_display_control_widgets(); 
}


void reset_graphics_display_control_window() { 

   graphics_info_t g; 
   g.save_display_control_widget_in_graphics(NULL); 

} 

/*! \brief make the map displayed/undisplayed, 0 for off, 1 for on */
void set_map_displayed(int imol, int state) {

   graphics_info_t g;
   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].set_map_is_displayed(state);
      if (g.display_control_window())
	 set_display_control_button_state(imol, "Displayed", state);
      graphics_draw();
   }
}

// button_type is "Displayed" or "Active"
void
set_display_control_button_state(int imol, const std::string &button_type, int state) {

   graphics_info_t g;
   if (g.display_control_window()) {
      
      std::string button_name = "";
      
      if (button_type == "Displayed") {
	 button_name = "displayed_button_";
	 button_name += coot::util::int_to_string(imol);
      }
      
      if (button_type == "Active") {
	 button_name = "active_button_";
	 button_name += coot::util::int_to_string(imol);
      }
      GtkWidget *button = lookup_widget(g.display_control_window(), button_name.c_str());
      if (button) {
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), state);
      } else {
	 std::cout << "Opps failed to find " << button_type << " button for mol "
		   << imol << std::endl;
      }
   }
}

/*! \brief make the coordinates molecule displayed/undisplayed, 0 for off, 1 for on */
void set_mol_displayed(int imol, int state) {
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_mol_is_displayed(state);
      if (g.display_control_window())
	 set_display_control_button_state(imol, "Displayed", state);
      graphics_draw();
   } else {
      std::cout << "not valid molecule" << std::endl;
   } 
} 
/*! \brief make the coordinates molecule active/inactve (clickable), 0
  for off, 1 for on */
void set_mol_active(int imol, int state) {
   
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_mol_is_active(state);
      if (g.display_control_window())
	 set_display_control_button_state(imol, "Active", state);
      graphics_draw();
   } else {
      std::cout << "not valid molecule" << std::endl;
   } 
} 



int mol_is_displayed(int imol) { 

   graphics_info_t g; 

   return g.molecules[imol].drawit; 

} 

int mol_is_active(int imol) { 

   graphics_info_t g; 

   return g.molecules[imol].atom_selection_is_pickable(); 

} 

int map_is_displayed(int imol) { 

   graphics_info_t g;
   return g.molecules[imol].drawit_for_map; 

}

/*! \brief if on_or_off is 0 turn off all maps displayed, for other
  values of on_or_off turn on all maps */
void set_all_maps_displayed(int on_or_off) {

   graphics_info_t g;
   int nm = graphics_info_t::n_molecules;
   for (int imol=0; imol<nm; imol++) {
      if (is_valid_map_molecule(imol)) {
	 graphics_info_t::molecules[imol].set_mol_is_active(on_or_off);
	 if (g.display_control_window())
	    set_display_control_button_state(imol, "Displayed", on_or_off);
      }
   }
   graphics_draw();
} 

/*! \brief if on_or_off is 0 turn off all models displayed and active,
  for other values of on_or_off turn on all models. */
void set_all_models_displayed_and_active(int on_or_off) {

   graphics_info_t g;
   int nm = graphics_info_t::n_molecules;
   for (int imol=0; imol<nm; imol++) {
      if (is_valid_model_molecule(imol)) {
	 graphics_info_t::molecules[imol].set_mol_is_active(on_or_off);
	 if (g.display_control_window())
	    set_display_control_button_state(imol, "Active", on_or_off);
	 if (g.display_control_window())
	    set_display_control_button_state(imol, "Displayed", on_or_off);
      }
   }
   graphics_draw();
}


char *show_spacegroup(int imol) { 

   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_map() || 
	  graphics_info_t::molecules[imol].has_model()) { 
	 std::string spg =  graphics_info_t::molecules[imol].show_spacegroup();
	 std::cout << "INFO:: spacegroup: " << spg << std::endl;
	 unsigned int l = spg.length();
	 char *s = new char[l+1];
	 strncpy(s, spg.c_str(), l+1);
	 return s;
      } else { 
	 std::string spg("No spacegroup for this molecule");
	 std::cout << "INFO:: spacegroup: " << spg << std::endl;
	 unsigned int l = spg.length();
	 char *s = new char[spg.length()+1];
	 strncpy(s, spg.c_str(), l+1);
	 return s;
      }
   }

   // If it was a bad molecule, return pointer to null.
   std::cout << "Unknown molecule " << imol << std::endl;
   char *s = new char[1];
   s[0] = 0;
   return s;
}


/*  ----------------------------------------------------------------------- */
/*                  zoom                                                    */
/*  ----------------------------------------------------------------------- */

void 
scale_zoom_internal(float f) {

   graphics_info_t g;
   g.zoom *= fabs(f);
}

void scale_zoom(float f) {

   graphics_info_t g;
   scale_zoom_internal(f);
   graphics_draw();

}

float zoom_factor() {
   graphics_info_t g;
   return g.zoom;
}

void set_smooth_scroll_do_zoom(int i) {
   graphics_info_t g;
   g.smooth_scroll_do_zoom = i; 
}


int  smooth_scroll_do_zoom() {

   return graphics_info_t::smooth_scroll_do_zoom;
} 


float smooth_scroll_zoom_limit() {

   return graphics_info_t::smooth_scroll_zoom_limit; 
}


void set_smooth_scroll_zoom_limit(float f) {

   graphics_info_t::smooth_scroll_zoom_limit = f; 
}

void set_zoom_adjustment(GtkWidget *w) { 
   graphics_info_t::set_zoom_adjustment(w);
} 




// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_python_window_entry(GtkWidget *entry) { 

#ifdef USE_PYTHON

   // add python entery in entry callback code here...

      gtk_signal_connect(GTK_OBJECT(entry), "activate",
		      GTK_SIGNAL_FUNC(python_window_enter_callback),
		      entry);

#endif

}

// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_guile_window_entry(GtkWidget *entry) { 

#ifdef USE_GUILE
   gtk_signal_connect(GTK_OBJECT(entry), "activate",
		      GTK_SIGNAL_FUNC(guile_window_enter_callback),
		      entry);

#endif //  USE_GUILE

}

#ifdef USE_PYTHON
void python_window_enter_callback( GtkWidget *widget,
				   GtkWidget *entry )
{
  const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
  printf("Entry contents: %s\n", entry_text);
  {
     char *py_text;
     py_text = new char [strlen(entry_text)+1];
     strncpy(py_text, entry_text, 500); 
     //
     // Py_Initialize();
     PyRun_SimpleString(py_text);

     // clear the entry
     gtk_entry_set_text(GTK_ENTRY(entry),"");
  }

}
#endif


#ifdef USE_GUILE
void guile_window_enter_callback( GtkWidget *widget,
				  GtkWidget *entry )
{
  const gchar *entry_text;
  entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
  printf("Entry contents: %s\n", entry_text);

  // scm_c_eval_string(entry_text); 

  // extern SCM scm_catch (SCM tag, SCM thunk, SCM handler);
  //

  //SCM handler = scm_c_eval_string("'(lambda (key . args))"); 

   SCM handler = scm_c_eval_string ("(lambda (key . args) "
     "(display (list \"Error in proc:\" key \" args: \" args)) (newline))"); 
			   // "(newline))"); 

  // scm_catch(SCM_BOOL_T, scm_c_eval_string(entry_text), handler); 

  std::string thunk("(lambda() "); 
  thunk += entry_text; 
  thunk += " )";


  SCM scm_thunk = scm_c_eval_string(thunk.c_str()); 
  scm_catch(SCM_BOOL_T, scm_thunk, handler);


  // clear the entry
  // strcpy(entry_text, "");  surely not needed.
  gtk_entry_set_text(GTK_ENTRY(entry),"");
}
#endif //  USE_GUILE


void clear_moving_atoms_object() {  /* redraw done here. */

   graphics_info_t g;
   g.clear_moving_atoms_object();
   
}

void clear_up_moving_atoms() {

   graphics_info_t g;
   std::cout << "c-interface clear_up_moving_atoms..." << std::endl;
   g.clear_up_moving_atoms();
   g.clear_moving_atoms_object();

} 

// Similar to fill_option_menu_with_coordinates_options, but I moved
// it to graphics_info_t because it is also used when there is an
// ambiguity in the map for refinement (graphics_info_t::refine)
// 
void fill_option_menu_with_map_options(GtkWidget *option_menu, GtkSignalFunc signalfunc) {

   graphics_info_t g;

   g.fill_option_menu_with_map_options(option_menu, signalfunc);
}

void fill_option_menu_with_skeleton_options(GtkWidget *option_menu) {  /* a wrapper */

   graphics_info_t g;
   GtkSignalFunc signalfunc = GTK_SIGNAL_FUNC(graphics_info_t::skeleton_map_select);
   g.fill_option_menu_with_map_options(option_menu, signalfunc,
				       graphics_info_t::map_for_skeletonize);

}


void set_initial_map_for_skeletonize() { 
   
   graphics_info_t::set_initial_map_for_skeletonize();

}

void set_max_skeleton_search_depth(int v) { 
   graphics_info_t g;
   g.set_max_skeleton_search_depth(v);
} 

/* Set the radio buttons in the frame to the be on or off for the map
   that is displayed in the optionmenu (those menu items "activate"
   callbacks (graphics_info::skeleton_map_select change
   g.map_for_skeletonize).  */
void set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame) { 

   graphics_info_t g;
   g.set_on_off_skeleton_radio_buttons(skeleton_frame);
}

void set_on_off_single_map_skeleton_radio_buttons(GtkWidget *skeleton_frame, 
						  int imol) { 
   
   graphics_info_t g;
   g.set_on_off_single_map_skeleton_radio_buttons(skeleton_frame, imol);

} 

void set_contour_sigma_button_and_entry(GtkWidget *window, int imol) { 
   graphics_info_t g;
   g.set_contour_sigma_button_and_entry(window, imol);
}

// and the reverse function (button -> value)
void set_contour_by_sigma_step_maybe(GtkWidget *window, int imol) { 

   GtkWidget *button = lookup_widget(window, "single_map_sigma_checkbutton");
   GtkWidget *entry  = lookup_widget(window, "single_map_sigma_step_entry");

   if (GTK_TOGGLE_BUTTON(button)->active) { 
      const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
      if (text) { 
	 float v = atof(text);
// 	 graphics_info_t::molecules[imol].contour_by_sigma_flag = 1;
// 	 graphics_info_t::molecules[imol].contour_sigma_step = v;
	 graphics_info_t::molecules[imol].set_contour_by_sigma_step(v, 1);
      }
   } else { 
      // 0.0 is ignored.
      graphics_info_t::molecules[imol].set_contour_by_sigma_step(0.0, 0);
   }
}

void set_contour_by_sigma_step_by_mol(float f, short int state, int imol) { 
   
   if (imol < graphics_info_t::n_molecules) { 
      if (imol >= 0) { 
	 if (graphics_info_t::molecules[imol].has_map()) {
	    graphics_info_t::molecules[imol].set_contour_by_sigma_step(f, state);
	 }
      }
   }
}

void export_map(int imol, const char *filename) {

   if (is_valid_map_molecule(imol)) {

      clipper::CCP4MAPfile mapout;
      mapout.open_write(std::string(filename));
      mapout.export_xmap(graphics_info_t::molecules[imol].xmap_list[0]);
      mapout.close_write(); 
      
   } else {
      graphics_info_t g;
      g.statusbar_text("Invalid map molecule number");
   }

}

int transform_map_raw(int imol, 
		      double r00, double r01, double r02, 
		      double r10, double r11, double r12, 
		      double r20, double r21, double r22, 
		      double t0, double t1, double t2,
		      double pt1, double pt2, double pt3, double box_size) {

   int imol_new = -1;
   if (is_valid_map_molecule(imol)) {
      clipper::Mat33<double> m(r00, r01, r02, r10, r11, r12, r20, r21, r22);
      clipper::Coord_orth c(t0, t1, t2);
      clipper::RTop_orth rtop(m,c);
      clipper::RTop_orth rtop_inv = rtop.inverse();
      clipper::Coord_orth pt(pt1, pt2, pt3);
      clipper::Xmap<float> new_map =
	 coot::util::transform_map(graphics_info_t::molecules[imol].xmap_list[0],
				   rtop, pt, box_size);

      const coot::ghost_molecule_display_t ghost_info;
      // int is_diff_map_flag = graphics_info_t::molecules[imol].is_difference_map_p();
      // int swap_colours_flag = graphics_info_t::swap_difference_map_colours;
      mean_and_variance<float> mv = map_density_distribution(new_map, 0);
      std::string name = "Transformed map";
      imol_new = graphics_info_t::n_molecules;
      graphics_info_t::molecules[imol_new].new_map(new_map, name);
      
      graphics_info_t::n_molecules++;
      graphics_draw();

   } else {
      std::cout << "molecule " << imol << " is not a valid map" << std::endl;
   }
   return imol_new;
}



void do_torsions_toggle(GtkWidget *button) {

   graphics_info_t g;
   GtkWidget *peptide_checkbutton =
      lookup_widget(button,
		    "refine_params_use_peptide_torsions_checkbutton");


   if (g.do_torsion_restraints) {
      g.do_torsion_restraints = 0;
      gtk_widget_set_sensitive(peptide_checkbutton, FALSE);
   } else {
      g.do_torsion_restraints = 1;
      gtk_widget_set_sensitive(peptide_checkbutton, TRUE);
   }
}

void do_peptide_torsions_toggle() {
   graphics_info_t g;
   if (g.do_peptide_torsion_restraints) {
      g.do_peptide_torsion_restraints = 0;
   } else {
      g.do_peptide_torsion_restraints = 1;
   }
}

void set_refine_params_toggle_buttons(GtkWidget *button) {

   // initiallly buttons are inactive and sensitive

   graphics_info_t g;
   GtkWidget *checkbutton =
      lookup_widget(button, "refine_params_use_torsions_checkbutton");
   GtkWidget *phi_psi_peptide_checkbutton =
      lookup_widget(button, "refine_params_use_peptide_torsions_checkbutton");
   GtkWidget *link_torsion_type_vbox =
      lookup_widget(button, "peptide_torsions_restraints_vbox");
   

   if (g.do_torsion_restraints) {
      g.do_torsion_restraints = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
   } else {
      gtk_widget_set_sensitive(GTK_WIDGET(phi_psi_peptide_checkbutton), FALSE);
   }

   if (g.do_peptide_torsion_restraints) {
      g.do_peptide_torsion_restraints = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(phi_psi_peptide_checkbutton), TRUE);
      gtk_widget_set_sensitive(link_torsion_type_vbox, TRUE);
   } else {
      gtk_widget_set_sensitive(link_torsion_type_vbox, FALSE);
   } 

   GtkWidget *omega = lookup_widget(button,
		       "refine_params_use_peptide_omegas_checkbutton");
   if (g.do_peptide_omega_torsion_restraints) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(omega), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(omega), FALSE);
   }


   // refine_params_use_helix_peptide_torsions_radiobutton
   // refine_params_use_beta_strand_peptide_torsions_radiobutton
   // refine_params_use_ramachandran_goodness_torsions_radiobutton

   GtkWidget *sec_str_rest_no_rest_radiobutton =
      lookup_widget(button, "sec_str_rest_no_rest_radiobutton");
   GtkWidget *sec_str_rest_helix_rest_radiobutton = 
      lookup_widget(button, "sec_str_rest_helix_rest_radiobutton");
   GtkWidget *sec_str_rest_strand_rest_radiobutton = 
      lookup_widget(button, "sec_str_rest_strand_rest_radiobutton");

#ifdef HAVE_GSL   
   if (graphics_info_t::pseudo_bonds_type == coot::NO_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_no_rest_radiobutton), TRUE);
   if (graphics_info_t::pseudo_bonds_type == coot::HELIX_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_helix_rest_radiobutton), TRUE);
   if (graphics_info_t::pseudo_bonds_type == coot::STRAND_PSEUDO_BONDS)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sec_str_rest_strand_rest_radiobutton), TRUE);
#endif // HAVE_GSL   

} 

// either alpha helix, beta strand or ramachandran goodness
// (see ideal/simple_restraint.hh link torsions)
void set_refine_params_phi_psi_restraints_type(int restraints_type) { 

  // do_peptide_omega_torsion_restraints are for phi/psi restraints 
  // (omega not included)
  graphics_info_t::do_peptide_torsion_restraints = restraints_type;

} 

void fill_chiral_volume_molecule_option_menu(GtkWidget *w) { 

   GtkWidget *optionmenu = lookup_widget(w, "check_chiral_volumes_molecule_optionmenu");

   // now set chiral_volume_molecule_option_menu_item_select_molecule to the top of the list
   for (int i=0; i<graphics_info_t::n_molecules; i++) { 
      if (graphics_info_t::molecules[i].has_model()) {
	 graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule = i;
	 break;
      } 
   }
   int imol = graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule;
   GtkSignalFunc callback_func =
      GTK_SIGNAL_FUNC(chiral_volume_molecule_option_menu_item_select);

   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(optionmenu, callback_func, imol);

}

void chiral_volume_molecule_option_menu_item_select(GtkWidget *item, GtkPositionType pos) { 

   graphics_info_t::chiral_volume_molecule_option_menu_item_select_molecule = pos;

}


void set_dragged_refinement_steps_per_frame(int v) {

   graphics_info_t g;
   g.dragged_refinement_steps_per_frame = v;
}

int dragged_refinement_steps_per_frame() {
   return graphics_info_t::dragged_refinement_steps_per_frame; 
} 



/*  ----------------------------------------------------------------------- */
/*                  get by accession code:                                  */
/*  ----------------------------------------------------------------------- */


/* Accession code, and dispatch guile command to download and display
   the model.  Hmmm.  */
void handle_get_accession_code(GtkWidget *widget) { 

   const gchar *text = gtk_entry_get_text(GTK_ENTRY(widget));
   cout << "PDB Accession Code: " << text << endl;
   int *n = (int *) gtk_object_get_user_data(GTK_OBJECT(lookup_widget(GTK_WIDGET(widget),
								      "accession_code_window")));

   

#ifdef USE_GUILE
   string scheme_command;

   if (*n == COOT_ACCESSION_CODE_WINDOW_EDS) {
      // 20050725 EDS code:

      scheme_command = "(get-eds-pdb-and-mtz ";
      scheme_command += single_quote(text);
      scheme_command += ")";
   } else { 
      if (*n == 1) { 
	 scheme_command = "(get-ebi-pdb \"";
      } else { 
	 // *n == 2 see callbacks.c on_get_pdb_and_sf_using_code1_activate
	 scheme_command = "(get-ebi-pdb-and-sfs \"";
      }
      scheme_command += text;
      scheme_command += "\")";
   }

   safe_scheme_command(scheme_command); 

#else 
   
   cout << "not compiled with guile.  This won't work" << endl; 

#endif // USE_GUILE

   // and kill the accession code window
   gtk_widget_destroy(lookup_widget(GTK_WIDGET(widget),
				    "accession_code_window")); 
   
} 





#ifdef USE_GUILE
SCM safe_scheme_command(const std::string &scheme_command) { 
   add_to_history_simple(scheme_command);
   return graphics_info_t::safe_scheme_command(scheme_command);
}
SCM safe_scheme_command_test(const char *cmd) {

   std::string s = cmd;
   return safe_scheme_command(s);
}
#else  // not guile
// dummy function
void safe_scheme_command(const std::string &scheme_command) { /* do nothing */
   // here only for compilation purposes.
}
#endif // USE_GUILE

void safe_python_command(const std::string &python_cmd) {

#ifdef USE_PYTHON
   PyRun_SimpleString((char *)python_cmd.c_str());
#endif   
}

void safe_python_command_by_char_star(const char *python_cmd) {

   std::cout <<  "DEBUG in safe_python_command_by_char_star: running: "
	     << python_cmd << std::endl;
   
#ifdef USE_PYTHON
   PyRun_SimpleString((char *)python_cmd);
#endif   
}




void post_scripting_window() {

#ifdef USE_GUILE
  // window = create_guile_window(); 
//   entry = lookup_widget(window, "guile_window_entry"); 
//   setup_guile_window_entry(entry); // USE_PYTHON and USE_GUILE used here

  if (graphics_info_t::guile_gui_loaded_flag == TRUE) { 

     scm_c_eval_string("(coot-gui)");

  } else { 
     // we don't get a proper status from guile_gui_loaded_flag so
     // lets check again here whether MAPVIEW_GUI_DIR was defined.
     char *t; 
     t = getenv(COOT_SCHEME_DIR); // was #defined
     if (t) { 
	std::cout << COOT_SCHEME_DIR << " was defined to be " << t << std::endl
		  << "   but loading of scripting window scheme code failed." 
		  << std::endl; 
     } else { 
	std::cout << COOT_SCHEME_DIR << " was not defined - cannot open ";
	std::cout << "scripting window" << std::endl; 
     } 
  } 
#endif


#ifdef USE_PYTHON

  GtkWidget *window; 
  GtkWidget *entry; 
  window = create_python_window();

  entry = lookup_widget(window, "python_window_entry");
  setup_python_window_entry(entry); // USE_PYTHON and USE_GUILE used here
  gtk_widget_show(window);

  // clear the entry here
#endif

}

/* called from c-inner-main */
void run_command_line_scripts() {

//     std::cout << "There are " << graphics_info_t::command_line_scripts->size() 
//  	     << " command line scripts to run\n";

    for (unsigned int i=0; i<graphics_info_t::command_line_scripts->size(); i++)
       run_script((*graphics_info_t::command_line_scripts)[i].c_str());
}


void
set_guile_gui_loaded_flag() { 
   
   graphics_info_t g; 
   g.guile_gui_loaded_flag = TRUE; 
} 

void set_found_coot_gui() { 
   
   cout << "Coot Scripting GUI code found and loaded." << endl; 
   graphics_info_t g; 
   g.guile_gui_loaded_flag = TRUE; 
}

// return an atom index
int atom_spec_to_atom_index(int imol, char *chain, int resno, char *atom_name) { 
   graphics_info_t g; 
   if (imol < graphics_n_molecules()) 
      return g.molecules[imol].atom_spec_to_atom_index(chain, resno, atom_name);
   else
      return -1;
}

int full_atom_spec_to_atom_index(int imol, const char *chain, int resno,
				 const char *inscode, const char *atom_name,
				 const char *altloc) {

   if (imol < graphics_n_molecules()) 
      return graphics_info_t::molecules[imol].full_atom_spec_to_atom_index(std::string(chain), resno, std::string(inscode), std::string(atom_name), std::string(altloc));
   else
      return -1;
}



// ??? FIXME for the future, this doesn't work when we have both guile
// and python.  We need to choose which script interpretter to use
// based on filename extension.
void
run_script(const char *filename) { 

   struct stat buf;
   int status = stat(filename, &buf);
   std::string fn(filename);
   if (status == 0) {

      short int is_python = 0;

      std::string::size_type ipy = fn.rfind(".py");
      if (ipy != std::string::npos) {
	 if (fn.substr(ipy) == ".py")
	    is_python = 1;
      }
	 
      if (is_python) { 
	 run_python_script(filename);
      } else { 
	 run_guile_script(filename);
      }
   } else { 
      std::cout  << "WARNING:: Can't run script: " << filename 
		 << " no such file." << std::endl;
   }
}

// If we have both GUILE and PYTHON, use the state file as if it were GUILE
// 
void
run_state_file() { 
   std::string filename;
#ifdef USE_GUILE
   filename = "0-coot.state.scm";
   struct stat buf;
   int status = stat(filename.c_str(), &buf);
   if (status == 0) { 
      run_guile_script(filename.c_str());
   }
#else 
#ifdef USE_PYTHON
   filename = "0-coot.state.py";
   struct stat buf;
   int status = stat(filename.c_str(), &buf);
   if (status == 0) { 
      run_python_script(filename.c_str());
   }
#endif
#endif
}


void
run_state_file_maybe() { 

   std::string filename("0-coot.state.scm");
#ifdef USE_PYTHON
#ifndef USE_GUILE
   filename = "0-coot.state.py";
#endif
#endif
   graphics_info_t g;

   /*  0: never run it */
   /*  1: ask to run it */
   /*  2: always run it */
   if (g.run_state_file_status == 1 || g.run_state_file_status == 2) { 
      
      // can we stat a status file?
      // 
      struct stat buf;
      int status = stat(filename.c_str(), &buf);
      if (status == 0) { 
	 if (g.run_state_file_status == 2) {
	    run_script(filename.c_str());
	 } else {
	    if (graphics_info_t::use_graphics_interface_flag) { 
	       GtkWidget *dialog = wrapped_create_run_state_file_dialog();
	       gtk_widget_show(dialog);
	    } 
	 }
      }
   }
}

GtkWidget *wrapped_create_run_state_file_dialog() {

   std::string filename("0-coot.state.scm");
   short int il = 1;
   GtkWidget *w = create_run_state_file_dialog();

   GtkWidget *vbox_mols = lookup_widget(w, "mols_vbox");

   graphics_info_t g;
   std::vector<std::string> v = g.save_state_data_and_models(filename, il);
   for (unsigned int i=0; i<v.size(); i++) { 
      //       std::cout << "Got molecule: " << v[i] << std::endl;
      std::string s = "    ";
      s += v[i];
      GtkWidget *label = gtk_label_new(s.c_str());
      gtk_misc_set_alignment (GTK_MISC (label), 0.0, 0.5);
      gtk_box_pack_start(GTK_BOX(vbox_mols), label, FALSE, FALSE, 2);
      gtk_widget_show(label);
   } 
   return w;
}


void
run_guile_script(const char *filename) { 

#ifdef USE_GUILE
   std::string thunk("(lambda() "); 
   thunk += "(load ";
   thunk += "\"";
   thunk += filename;
   thunk += "\"))";

   SCM handler = scm_c_eval_string ("(lambda (key . args) "
     "(display (list \"Error in proc:\" key \" args: \" args)) (newline))"); 

   SCM scm_thunk = scm_c_eval_string(thunk.c_str());
   scm_catch(SCM_BOOL_T, scm_thunk, handler);
#endif // USE_GUILE   

} 

void
run_python_script(const char *filename_in) { 

#ifdef USE_PYTHON

   std::string s = coot::util::intelligent_debackslash(filename_in);
   std::string simple = "execfile(";
   simple += single_quote(s);
   simple += ")";
   std::cout << "Running python script " << s  << std::endl;
   // not a const argument?  Dear oh dear....
   PyRun_SimpleString((char *)simple.c_str());

#endif // USE_PYTHON
}



void add_on_rama_choices(){  // the the menu

   // std::cout << "adding rama molecule options:" << std::endl;

   // first delete all the current menu items.
   //
   graphics_info_t g;
   GtkWidget* menu = lookup_widget(GTK_WIDGET(g.glarea), "rama_plot_menu");

   if (menu) {
      gtk_container_foreach(GTK_CONTAINER(menu),
			    my_delete_ramachandran_mol_option,
			    (gpointer) menu);
   
      std::string name;
      for (int i=0; i<g.n_molecules; i++) {
	 if (g.molecules[i].has_model() > 0) {
	    name = graphics_info_t::molecules[i].dotted_chopped_name();
	    update_ramachandran_plot_menu_manual(i, name.c_str());
	 }
      }
   }
}


void destroy_edit_backbone_rama_plot() { 

   graphics_info_t g;
   g.destroy_edit_backbone_rama_plot();

} 


/*  ----------------------------------------------------------------------- */
/*           sequence_view                                                  */
/*  ----------------------------------------------------------------------- */

// A pure sequence function, not sequence view, so that people can cut
// n paste the sequence of a pdb file from the console.  There will be
// a scripting level function called print-sequence that gets called
// for every chain in the mol
// 
void print_sequence_chain(int imol, const char *chain_id) {

   std::string seq;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (std::string(chain_p->GetChainID()) == chain_id) { 
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    int residue_count_block = 0;
	    int residue_count_line = 0;
	    if (nres > 0 ) {
	       residue_count_block = chain_p->GetResidue(0)->GetSeqNum();
	       residue_count_line  = residue_count_block;
	       if (residue_count_block > 0)
		  while (residue_count_block > 10)
		     residue_count_block -= 10;
	       if (residue_count_line > 0)
		  while (residue_count_line > 50)
		     residue_count_line -= 50;
	    }
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       seq += coot::util::three_letter_to_one_letter(residue_p->GetResName());
	       if (residue_count_block == 10) {
		  seq += " ";
		  residue_count_block = 0;
	       }
	       if (residue_count_line == 50) {
		  seq += "\n";
		  residue_count_line = 0;
	       }
	       residue_count_block++;
	       residue_count_line++;
	    }
	 }
      }
      std::cout << "> " << graphics_info_t::molecules[imol].name_sans_extension(0)
		<< " chain " << chain_id << std::endl;
      std::cout << seq << std::endl;
   }
}

void do_sequence_view(int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   graphics_info_t g;
   if (g.molecules[imol].has_model()) {
      graphics_info_t g;

      if (g.sequence_view_is_displayed[imol] != 0) {

	 // it already exists... just raise it and map it.

	 GtkWidget *canvas = g.sequence_view_is_displayed[imol];
	 // so what is the window (which we shall call widget)?
	 GtkWidget *widget = lookup_widget(canvas, "sequence_view_dialog");

	 if (!GTK_WIDGET_MAPPED(widget)) {
	    gtk_widget_show(widget);
	 } else {
	    gdk_window_raise(widget->window);
	 }

      } else {

	 // create a new one

	 // Let's have a name that has the leading / and .pdb stripped
	 //
	 std::string short_name;

	 std::string::size_type islash = g.molecules[imol].name_.find_last_of("/");
	 std::string tstring;
	 if (islash == string::npos) { 
	    // no slash found
	    tstring = g.molecules[imol].name_;
	 } else { 
	    tstring = g.molecules[imol].name_.substr(islash + 1);
	 }
      
	 std::string::size_type ipdb = tstring.rfind(".pdb");

	 if (ipdb == string::npos) {
	    std::cout << "INFO .pdb not found in filename" << std::endl;
	    short_name = tstring;
	 } else {
	    short_name = tstring.substr(0, ipdb);
	 }
	 coot::sequence_view *seq_view =
	    new coot::sequence_view(g.molecules[imol].atom_sel.mol,
				    short_name, imol);
	 g.set_sequence_view_is_displayed(seq_view->Canvas(), imol);
      }
   }
#endif // HAVE_GTK_CANVAS   
}

void
add_on_sequence_view_choices() {

   graphics_info_t g;
   GtkWidget *menu = lookup_widget(GTK_WIDGET(g.glarea), "seq_view_menu");

   if (menu) { 
      gtk_container_foreach(GTK_CONTAINER(menu),
			    my_delete_ramachandran_mol_option,
			    (gpointer) menu);
      for(int i=0; i<g.n_molecules; i++) {
	 if (g.molecules[i].has_model()) {
	    std::string name;
	    name = graphics_info_t::molecules[i].dotted_chopped_name();
	    update_sequence_view_menu_manual(i, name.c_str());
	 }
      }
   }
} 


void set_sequence_view_is_displayed(GtkWidget *widget, int imol) {
   graphics_info_t g;
   g.set_sequence_view_is_displayed(widget, imol);
} 

/*  ----------------------------------------------------------------------- */
/*           rotate moving atoms peptide                                    */
/*  ----------------------------------------------------------------------- */

void change_peptide_carbonyl_by(double angle) { /* in degrees. */
   graphics_info_t g;
   g.change_peptide_carbonyl_by(angle);
} 

void change_peptide_peptide_by(double angle) {   /* in degress */
   graphics_info_t g;
   g.change_peptide_peptide_by(angle);
}

void execute_setup_backbone_torsion_edit(int imol, int atom_index) {
   graphics_info_t g;
   g.execute_setup_backbone_torsion_edit(imol, atom_index);
}

void setup_backbone_torsion_edit(short int state) { 

   graphics_info_t g;
   graphics_info_t::in_backbone_torsion_define = state;
   if (state) { 
      std::cout << "click on an atom in the peptide to change" << std::endl; 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
   } else { 
      g.normal_cursor();
   }
}

void set_refine_with_torsion_restraints(int istate) {

   graphics_info_t::do_torsion_restraints = istate;

} 

   

void set_backbone_torsion_peptide_button_start_pos(int ix, int iy) { 
   
   graphics_info_t g;
   g.set_backbone_torsion_peptide_button_start_pos(ix, iy);
} 

void change_peptide_peptide_by_current_button_pos(int ix, int iy) { 

   graphics_info_t g;
   g.change_peptide_peptide_by_current_button_pos(ix, iy);
}

void set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy) { 

   graphics_info_t g;
   g.set_backbone_torsion_carbonyl_button_start_pos(ix, iy);

} 

void change_peptide_carbonyl_by_current_button_pos(int ix, int iy) { 

   graphics_info_t g;
   g.change_peptide_carbonyl_by_current_button_pos(ix, iy);

} 

/*  ----------------------------------------------------------------------- */
/*                  cif stuff                                               */
/*  ----------------------------------------------------------------------- */

// and make (and display) a sigma_a map.
// 
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
// 
int read_cif_data(const char *filename, int imol_coordinates) {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      cout << "Reading cif file: " << filename << endl; 
      graphics_info_t g; 
      int imol = g.n_molecules;
      int istat = g.molecules[imol].make_map_from_cif(std::string(filename), imol_coordinates);

      // std::cout << "DEBUG:: in read_cif_data, istat is " << istat << std::endl;
      if (istat != -1) { 
	 g.n_molecules++;
	 graphics_draw();
      }
      return imol;
   }
   return -1; // which is status in an error
}

// and make (and display) a 2fofc map.
// 
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
// 
int read_cif_data_2fofc_map(const char *filename, int imol_coordinates) {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      cout << "Reading cif file: " << filename << endl; 

      graphics_info_t g; 

      int imol = g.n_molecules;

      int istat = g.molecules[imol].make_map_from_cif_2fofc(std::string(filename), imol_coordinates);

      if (istat != -1) { 
	 g.n_molecules++;
	 graphics_draw();
	 return imol;
      }
      return -1; // an error
   }
}


// and make (and display) a fofc map.
// 
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
// 
int read_cif_data_fofc_map(const char *filename, int imol_coordinates) {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      cout << "Reading cif file: " << filename << endl; 

      graphics_info_t g; 

      int imol = g.n_molecules;

      int istat = g.molecules[imol].make_map_from_cif_fofc(std::string(filename), imol_coordinates);

      if (istat != -1) { 
	 g.n_molecules++;
	 graphics_draw();
	 return imol;
      }
      return -1; // an error
   }
}



// This cif file, we presume, has phases.
// So we don't need a molecule to calculate them from.
// 
int auto_read_cif_data_with_phases(const char *filename) {

   int returned_mol_index = read_cif_data_with_phases_sigmaa(filename);
   read_cif_data_with_phases_diff_sigmaa(filename);
   return returned_mol_index;
}

int read_cif_data_with_phases_sigmaa(const char *filename) {

   graphics_info_t g; 
   int imol = g.n_molecules;
   
   int returned_mol_index = -1;
   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      std::cout << "Reading cif file: " << filename << std::endl; 

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.
      std::string fn = filename;
      int istat = g.molecules[imol].make_map_from_cif(fn);
      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 imol = -1;
      }
   }
   return imol;
}

int read_cif_data_with_phases_diff_sigmaa(const char *filename) {

   graphics_info_t g; 
   int imol = g.n_molecules;
   
   int returned_mol_index = -1;
   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      std::cout << "Reading cif file: " << filename << std::endl; 

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.
      std::string fn = filename;
      int istat = g.molecules[imol].make_map_from_cif_diff_sigmaa(fn);
      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 imol = -1;
      }
   }
   return imol;
} 


int read_cif_data_with_phases_fo_fc(const char *filename) {

   return read_cif_data_with_phases_nfo_fc(filename, molecule_map_type::TYPE_FO_FC);

} 

int read_cif_data_with_phases_2fo_fc(const char *filename) {

   return read_cif_data_with_phases_nfo_fc(filename, molecule_map_type::TYPE_2FO_FC);
}

int read_cif_data_with_phases_fo_alpha_calc(const char *filename) {
   return read_cif_data_with_phases_nfo_fc(filename, molecule_map_type::TYPE_FO_ALPHA_CALC);
}

int read_cif_data_with_phases_nfo_fc(const char *filename,
				     int map_type) {
   // first, does the file exist?
   struct stat s; 
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   } else {
      
      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

      graphics_info_t g; 

      int imol = g.n_molecules;
      std::string f(filename);
      short int swap_col = graphics_info_t::swap_difference_map_colours;

      int istat = g.molecules[imol].make_map_from_cif_nfofc(f, map_type, swap_col);

      if (istat != -1) {

	 g.scroll_wheel_map = g.n_molecules; // change the current scrollable map.
	 g.n_molecules++;
	 graphics_draw();
	 return imol;
      }
      return -1; // error
   }
}

int handle_shelx_fcf_file_internal(const char *filename) {

   graphics_info_t g;
   std::vector<std::string> cmd;
   cmd.push_back("handle-shelx-fcf-file");
   cmd.push_back(single_quote(filename));

#ifdef USE_GUILE   
   std::string s = g.state_command(cmd, coot::STATE_SCM);
   safe_scheme_command(s);
#endif

#ifdef USE_PYTHON
#ifndef USE_GUILE
   std::string s = g.state_command(cmd, coot::STATE_PYTHON);
   safe_python_command(s);
#endif    
#endif    
   return 1; // needed so that try_read_cif_file succeeds, so that we
	     // don't try to read this file as CNS data.
}


/* Use the environment variable COOT_REFMAC_LIB_DIR to find cif files
   in subdirectories and import them all. */
void import_all_refmac_cifs() {

   graphics_info_t g;
   g.import_all_refmac_cifs();

} 



/*  ----------------------------------------------------------------------- */
/*                  atom labelling                                          */
/*  ----------------------------------------------------------------------- */
/* The guts happens in molecule_class_info_t, here is just the
   exported interface */
int add_atom_label(int imol, char *chain_id, int iresno, char *atom_id) {

   int i = 0;
   if (is_valid_model_molecule(imol)) { 
      graphics_info_t g;
      i = g.molecules[imol].add_atom_label(chain_id, iresno, atom_id);
   }
   return i;
} 

int remove_atom_label(int imol, char *chain_id, int iresno, char *atom_id) {
   graphics_info_t g;
   return g.molecules[imol].remove_atom_label(chain_id, iresno, atom_id);
} 

void remove_all_atom_labels() { 
   graphics_info_t g;
   g.remove_all_atom_labels();
} 

void set_label_on_recentre_flag(int i) { 
   graphics_info_t::label_atom_on_recentre_flag = i;
} 

int centre_atom_label_status() { 
   return graphics_info_t::label_atom_on_recentre_flag;
} 

void set_brief_atom_labels(int istat) {
   graphics_info_t::brief_atom_labels_flag = istat;
   graphics_draw();
}

int brief_atom_labels_state() {
   return graphics_info_t::brief_atom_labels_flag;
}

/*  ----------------------------------------------------------------------- */
/*                  scene rotation (by script)                              */
/*  ----------------------------------------------------------------------- */
/* stepsize in degrees */
void rotate_y_scene(int nsteps, float stepsize) { 

   float spin_quat[4];
   graphics_info_t g;

   // spin it 1 degree
   float tbs =  g.get_trackball_size(); 
   for(int i=0; i<nsteps; i++) { 
     trackball(spin_quat, 0, 0, 0.0174*stepsize, 0.000, tbs);
     add_quats(spin_quat, g.quat, g.quat);
     graphics_draw();
   }
} 

/* stepsize in degrees */
void rotate_x_scene(int nsteps, float stepsize) { 

  float spin_quat[4];
   graphics_info_t g;

   // spin it 1 degree
   float tbs =  g.get_trackball_size(); 
   for(int i=0; i<nsteps; i++) { 
     trackball(spin_quat, 0, 0, 0.0, 0.0174*stepsize, tbs);
     add_quats(spin_quat, g.quat, g.quat);
     graphics_draw();
   }
} 

void rotate_z_scene(int nsteps, float stepsize) { 

   // c.f globjects.cc:do_screen_z_rotate()
   // 

   float spin_quat[4];
   graphics_info_t g;
   for(int i=0; i<nsteps; i++) { 
      trackball(spin_quat, 
		1.0, 1.0,
		1.0, 1.0 + 0.0174*stepsize,
		0.4);
      add_quats(spin_quat, g.quat, g.quat);
      graphics_draw();
   }
}

/*! \brief Bells and whistles rotation 

    spin, zoom and translate.

    where axis is either x,y or z,
    stepsize is in degrees, 
    zoom_by and x_rel etc are how much zoom, x,y,z should 
            have changed by after nstep steps.
*/
void spin_zoom_trans(int axis, int nsteps, float stepsize, float zoom_by, 
		     float x_rel, float y_rel, float z_rel) {

   float spin_quat[4];
   graphics_info_t g;
   float tbs =  g.get_trackball_size();
   float x_frag = 0.0;
   float y_frag = 0.0;
   float z_frag = 0.0;
   if (nsteps != 0) { 
      x_frag = x_rel/float(nsteps);
      y_frag = y_rel/float(nsteps);
      z_frag = z_rel/float(nsteps);
   }
   float zoom_init = g.zoom;
   float zoom_final = g.zoom * zoom_by;
   float zoom_frag = 1.0;
   if (nsteps !=0) {
      zoom_frag = (zoom_final - zoom_init)/float(nsteps);
   }
   int sss = graphics_info_t::smooth_scroll;

   // setRotationCentre looks at this and does a smooth scroll if its
   // on.
   graphics_info_t::smooth_scroll = 0;
   
   // std::cout << "zoom_frag is " << zoom_frag << std::endl;
   for(int i=0; i<nsteps; i++) {
      if (axis == 1) { 
	 trackball(spin_quat, 0, 0, 0.0, 0.0174*stepsize, tbs);
	 add_quats(spin_quat, g.quat, g.quat);
      }
      if (axis == 2) { 
	 trackball(spin_quat, 0, 0, 0.0174*stepsize, 0.000, tbs);
	 add_quats(spin_quat, g.quat, g.quat);
      }
      if (axis == 3) { 
	 trackball(spin_quat, 1.0, 1.0, 1.0, 1.0 + 0.0174*stepsize, 0.4);
	 add_quats(spin_quat, g.quat, g.quat);
      }
      g.zoom = zoom_init + float(i+1)*zoom_frag;
      coot::Cartesian c(g.X() + x_frag, g.Y() + y_frag, g.Z() + z_frag);
      g.setRotationCentre(c);
      graphics_draw();
   }
   graphics_info_t::smooth_scroll = sss;
} 




/*  ----------------------------------------------------------------------- */
/*                  graphics background colour                              */
/*  ----------------------------------------------------------------------- */
/* stepsize in degrees */
void set_background_colour(double red, double green, double blue) {

   graphics_info_t g;

   glClearColor(red,green,blue,1.0);
   g.background_colour[0] = red; 
   g.background_colour[1] = green; 
   g.background_colour[2] = blue; 
   glFogfv(GL_FOG_COLOR, g.background_colour);
   graphics_draw();

}

int  background_is_black_p() {

   int v = 0;
   graphics_info_t g;
   if (g.background_colour[0] < 0.1)
      if (g.background_colour[1] < 0.1)
	 if (g.background_colour[2] < 0.1)
	    v = 1;

   return v;
}



// ------------------------------------------------------------------
//                                Utility
// ------------------------------------------------------------------
// 
// File system Utility function: maybe there is a better place for it...
// Return like mkdir: mkdir returns zero on success, or -1 if an  error  occurred
//
// if it already exists as a dir, return 0 of course.
// 
int
make_directory_maybe(const char *dir) {
   return coot::util::create_directory(std::string(dir));
} 


void add_coordinates_glob_extension(const char *ext) { 
   
   graphics_info_t g;
   g.add_coordinates_glob_extension(std::string(ext));
} 

void add_data_glob_extension(const char *ext) { 
   graphics_info_t g;
   g.add_data_glob_extension(std::string(ext));
} 

void add_dictionary_glob_extension(const char *ext) { 
   graphics_info_t g;
   g.add_dictionary_glob_extension(std::string(ext));
} 

void add_map_glob_extension(const char *ext) { 
   graphics_info_t g;
   g.add_map_glob_extension(std::string(ext));
} 

int do_anti_aliasing_state() {
   return graphics_info_t::do_anti_aliasing_flag;
} 


void set_do_anti_aliasing(int state) {

   graphics_info_t g;
   g.set_do_anti_aliasing(state);
} 


void set_do_GL_lighting(int state) {
   graphics_info_t::do_lighting_flag = state;
   setup_lighting(state);
   graphics_draw();
}


int do_GL_lighting_state() {
   return graphics_info_t::do_lighting_flag;
}


// Glaxo people: Chuang Chang (surely not how you spell it), Moira X.
// 

/*  ----------------------------------------------------------------------- */
/*                  crosshairs                                              */
/*  ----------------------------------------------------------------------- */
void set_draw_crosshairs(short int i) { 

   graphics_info_t g;
   g.draw_crosshairs_flag = i;
   if (i > 0 ) { 
      g.crosshairs_text(); 
      graphics_draw();
   }
} 

short int draw_crosshairs_state() {
   return graphics_info_t::draw_crosshairs_flag; 
}


/*  ----------------------------------------------------------------------- */
/*                  citation notice                                         */
/*  ----------------------------------------------------------------------- */
void citation_notice_off() { 

   graphics_info_t::show_citation_notice = 0;

} 

/*  ----------------------------------------------------------------------- */
/*                  cursor function                                         */
/*  ----------------------------------------------------------------------- */
void normal_cursor() {

   graphics_info_t g;
   g.normal_cursor();
   graphics_draw();
}

void fleur_cursor() {
   graphics_info_t g;
   g.fleur_cursor();
   graphics_draw();

}

void pick_cursor_maybe() {

   graphics_info_t g;
   g.pick_cursor_maybe();
   graphics_draw();
}

void rotate_cursor() {
   normal_cursor();
}

void set_pick_cursor_index(int i) {
   graphics_info_t::pick_cursor_index = GdkCursorType(i);
}


/*  ------------------------------------------------------------------------ */
/*                       povray/raster3d interface                           */
/*  ------------------------------------------------------------------------ */
void raster3d(const char *filename) {

   graphics_info_t g;
   g.raster3d(std::string(filename));
}

void povray(const char *filename) {

   graphics_info_t g;
   g.povray(std::string(filename));
}

void set_raster3d_bond_thickness(float f) { 

   graphics_info_t::raster3d_bond_thickness = f;

} 

void set_raster3d_atom_radius(float f) { 

   graphics_info_t::raster3d_atom_radius = f;

} 

void set_raster3d_density_thickness(float f) {

   graphics_info_t::raster3d_density_thickness = f;

} 

void
raster_screen_shot() {  // run raster3d or povray and guile
                		         // script to render and display image

   // do some checking for povray/render here:

   std::string cmd("(render-image)");  // this is a render function 
   
   // cmd = "(povray-image)";

   safe_scheme_command(cmd);
}


void set_renderer_show_atoms(int istate) {

   graphics_info_t::renderer_show_atoms_flag = istate;
} 


/*  ----------------------------------------------------------------------- */
/*                  browser url                                          */
/*  ----------------------------------------------------------------------- */
void browser_url(const char *url) {

   if (url) { 
      std::string u(url);
      std::vector<std::string> commands;
      commands.push_back("system");
      std::string s = graphics_info_t::browser_open_command;
      if (s == "firefox" || s == "mozilla" || s == "netscape") { 
	 s += " -remote 'openURL(\\\"";
	 s += u;
	 s += "\\\",new-window)'";
	 commands.push_back(single_quote(s));
      } else {
	 if (s == "open") {
	    s += " ";
	    s += url;
	 } else {
	    s += " ";
	    s += url;
	 }
	 commands.push_back(single_quote(s));
      }

      std::string c = languagize_command(commands);
#ifdef USE_GUILE
      safe_scheme_command(c);
#else
#ifdef USE_PYTHON
      safe_python_command(c);
#endif
#endif
   }
}


void set_browser_interface(const char *browser) {

   if (browser) {
      graphics_info_t::browser_open_command = browser;
   } 
} 

void handle_online_coot_search_request(const char *entry_text) {

   if (entry_text) {
      clipper::String text(entry_text);
      std::vector<clipper::String> bits = text.split(" ");
      if (bits.size() > 0) { 
	 std::string s = "http://www.google.co.uk/search?q=";
	 s += bits[0];
	 for (int i=1; i<bits.size(); i++) {
	    s += "+";
	    s += bits[i];
	 }
	 s += "+coot+site%3Awww.ysbl.york.ac.uk";
	 browser_url(s.c_str());
      }
   } 
} 


/*  ----------------------------------------------------------------------- */
/*                  remote control                                          */
/*  ----------------------------------------------------------------------- */
/* section Remote Control */

// called by c_inner_main() if we have guile
void make_socket_listener_maybe() {


   std::vector<std::string> cmd;

   if (graphics_info_t::try_port_listener) { 
      cmd.push_back("open-coot-listener-socket");
      cmd.push_back(graphics_info_t::int_to_string(graphics_info_t::remote_control_port_number));
      cmd.push_back(single_quote(graphics_info_t::remote_control_hostname));

      // This is not a static function, perhaps it should be:
      graphics_info_t g;
      std::string scm_command = g.state_command(cmd, coot::STATE_SCM);

      safe_scheme_command(scm_command);

      if (graphics_info_t::coot_socket_listener_idle_function_token == -1)
	 if (graphics_info_t::listener_socket_have_good_socket_state) 
	    graphics_info_t::coot_socket_listener_idle_function_token =
	       gtk_idle_add((GtkFunction) coot_socket_listener_idle_func,
			    graphics_info_t::glarea);
   }
}

void set_coot_listener_socket_state_internal(int sock_state) {
   graphics_info_t::listener_socket_have_good_socket_state = sock_state;
}



int coot_socket_listener_idle_func(GtkWidget *w) { 

#ifdef USE_GUILE
   std::cout << "DEBUG:: running socket idle function" << std::endl;
   if (graphics_info_t::listener_socket_have_good_socket_state) { 
      std::cout << "DEBUG:: running guile function" << std::endl;
      safe_scheme_command("(coot-listener-idle-function-proc)");
   }
#endif
   return 1;
}

/*  ----------------------------------------------------------------------- */
/*                  Surfaces                                                */
/*  ----------------------------------------------------------------------- */
void do_surface(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].make_surface(state);
      graphics_draw();
   }
}



/*  ----------------------------------------------------------------------- */
/*                  Views                                                   */
/*  ----------------------------------------------------------------------- */
int add_view_here(const char *view_name) {

   std::string name(view_name);
   float quat[4];
   for (int i=0; i<4; i++)
      quat[i] = graphics_info_t::quat[i];
   graphics_info_t g;
   coot::Cartesian rc = g.RotationCentre();
   float zoom = graphics_info_t::zoom;
   coot::view_info_t view(quat, rc, zoom, name);
   graphics_info_t::views->push_back(view);
   return (graphics_info_t::views->size() -1);
}

int add_view_raw(float rcx, float rcy, float rcz, float quat0, float quat1, 
		  float quat2, float quat3, float zoom, const char *view_name) {

   float quat[4];
   quat[0] = quat0;
   quat[1] = quat1;
   quat[2] = quat2;
   quat[3] = quat3;
   coot::Cartesian rc(rcx, rcy, rcz);
   coot::view_info_t v(quat, rc, zoom, view_name);
   graphics_info_t::views->push_back(v);
   return (graphics_info_t::views->size() -1);
}


int remove_named_view(const char *view_name) {

   int r=0;
   bool found = 0;
   std::string vn(view_name);
   std::vector<coot::view_info_t> new_views;

   // don't push back the view if it has the same name and if we
   // haven't found a view with that name before.
   for (int iv=0; iv<graphics_info_t::views->size(); iv++) {
      if ((*graphics_info_t::views)[iv].view_name == vn) {
	 if (found == 1) {
	    new_views.push_back((*graphics_info_t::views)[iv]);
	 }
	 found = 1;
      } else {
	 new_views.push_back((*graphics_info_t::views)[iv]);
      }
   }
   if (found) {
      r = 1;
      *graphics_info_t::views = new_views;
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("remove_named_view");
   command_strings.push_back(single_quote(coot::util::intelligent_debackslash(view_name)));
   add_to_history(command_strings);
   return r;
} 

/*! \brief the given view number */
void remove_view(int view_number) {

   std::vector<coot::view_info_t> new_views;
   for (int iv=0; iv<graphics_info_t::views->size(); iv++) {
      if (iv != view_number)
	 new_views.push_back((*graphics_info_t::views)[iv]);
   }
   *graphics_info_t::views = new_views;

   std::string cmd = "remove-view";
   std::vector<coot::command_arg_t> args;
   args.push_back(view_number);
   add_to_history_typed(cmd, args);
}


void play_views() {

   int nsteps = 2000;
   if (graphics_info_t::views_play_speed > 0.000000001)
      nsteps = int(2000.0/graphics_info_t::views_play_speed);
   float play_speed = 1.0; 
   if (graphics_info_t::views_play_speed > 0.0)
      play_speed = graphics_info_t::views_play_speed;
   
//    std::cout << "DEBUG:: # Views "<< graphics_info_t::views->size() << std::endl;
   for (int iv=0; iv<graphics_info_t::views->size(); iv++) {
      coot::view_info_t view1 = (*graphics_info_t::views)[iv];
      std::cout << "DEBUG:: View "<< iv << " " << view1.view_name << std::endl;
      if (! (view1.is_simple_spin_view_flag ||
	     view1.is_action_view_flag)) {
	 if ((iv+1) < graphics_info_t::views->size()) { 
// 	    std::cout << "DEBUG:: interpolating to  "<< iv+1 << " "
// 		      << view1.view_name << std::endl;
	    coot::view_info_t view2 = (*graphics_info_t::views)[iv+1];
	    if (! (view2.is_simple_spin_view_flag ||
		   view2.is_action_view_flag)) { 
	       coot::view_info_t::interpolate(view1, view2, nsteps);
	       update_things_on_move_and_redraw();
	    }
	 }
      } else {
	 // a simple spin  or an action view here:
// 	    std::cout << "DEBUG:: simple spin "
// 		      << view1.view_name << std::endl;
	 int n_spin_steps = int (float (view1.n_spin_steps) / play_speed);
	 float dps = view1.degrees_per_step*0.5 * play_speed;
	 rotate_y_scene(n_spin_steps, dps);
	 if ((iv+1) < graphics_info_t::views->size()) { 
 	    std::cout << "DEBUG:: interpolating to  "<< iv+1 << " "
 		      << view1.view_name << std::endl;
	    coot::view_info_t view2 = (*graphics_info_t::views)[iv+1];
	    if (!view2.is_simple_spin_view_flag && !view2.is_action_view_flag) {
	       // the quat was not set because this is a simple
	       // rotate, so we must generate it from the current
	       // position
	       coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
				  graphics_info_t::RotationCentre_y(),
				  graphics_info_t::RotationCentre_z());
	       coot::view_info_t current_view(graphics_info_t::quat,
					      rc, graphics_info_t::zoom, "dummy");
	       coot::view_info_t::interpolate(current_view, view2, nsteps);
	       update_things_on_move_and_redraw();
	    }
	 }
      }
   }
   add_to_history_simple("play-views");
}

void remove_this_view() {

   graphics_info_t g;
   coot::Cartesian rc = g.RotationCentre();
   float quat[4];
   for (int i=0; i<4; i++)
      quat[i] = graphics_info_t::quat[i];
   float zoom =  g.zoom;

   int r=0;
   bool found = 0;
   coot::view_info_t v(quat, rc, zoom, "");

   // don't push back the view if it has the same name and if we
   // haven't found a view with that name before.
   std::vector<coot::view_info_t> new_views;
   for (int iv=0; iv<graphics_info_t::views->size(); iv++) {
      if ((*graphics_info_t::views)[iv].matches_view(v)) {
	 if (found == 1) {
	    new_views.push_back((*graphics_info_t::views)[iv]);
	 }
      } else {
	 new_views.push_back((*graphics_info_t::views)[iv]);
      }
   }
   if (found) {
      r = 1;
      *graphics_info_t::views = new_views;
   }
   add_to_history_simple("remove-this-view");
}


int go_to_first_view(int snap_to_view_flag) {

   std::string cmd = "go-to-first-view";
   std::vector<coot::command_arg_t> args;
   args.push_back(snap_to_view_flag);
   add_to_history_typed(cmd, args);
   return go_to_view_number(0, snap_to_view_flag);
}

void clear_all_views() {
   graphics_info_t::views->clear();
}

int go_to_view_number(int view_number, int snap_to_view_flag) {

   int r = 0;
   graphics_info_t g;
   if ((graphics_info_t::views->size() > view_number) && (view_number >= 0)) {
      coot::view_info_t view = (*graphics_info_t::views)[view_number];
      if (view.is_simple_spin_view_flag) {
	 int nsteps = 2000;
	 if (graphics_info_t::views_play_speed > 0.000000001)
	    nsteps = int(2000.0/graphics_info_t::views_play_speed);
	 float play_speed = 1.0; 
	 if (graphics_info_t::views_play_speed > 0.0)
	    play_speed = graphics_info_t::views_play_speed;
	 int n_spin_steps = int (float (view.n_spin_steps) / play_speed);
	 float dps = view.degrees_per_step*0.5 * play_speed;
	 rotate_y_scene(n_spin_steps, dps);
      } else {
	 if (view.is_action_view_flag) {
	    // do nothing
	 } else {
	    if (snap_to_view_flag) {
	       g.setRotationCentre(view.rotation_centre);
	       g.zoom = view.zoom;
	       for (int iq=0; iq<4; iq++)
		  g.quat[iq] = view.quat[iq];
	    } else {
	       coot::view_info_t this_view(g.quat, g.RotationCentre(), g.zoom, "");
	       int nsteps = 2000;
	       if (graphics_info_t::views_play_speed > 0.000000001)
		  nsteps = int(2000.0/graphics_info_t::views_play_speed);
	       coot::view_info_t::interpolate(this_view,
					      (*graphics_info_t::views)[view_number], nsteps);
	    }
	 }
	 update_things_on_move_and_redraw();
      }
   }
   std::string cmd = "go-to-view-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(view_number);
   args.push_back(snap_to_view_flag);
   add_to_history_typed(cmd, args);
   return r;
}


/*! \brief return the number of views */
int n_views() {
   add_to_history_simple("n-views");
   return graphics_info_t::views->size();
}

/*! \brief return the name of the given view, if view_number does not
  specify a view return #f */
#ifdef USE_GUILE
SCM view_name(int view_number) {

   SCM r = SCM_BOOL_F;
   int n_view = graphics_info_t::views->size();
   if (view_number < n_view)
      if (view_number >= 0) {
	 std::string name = (*graphics_info_t::views)[view_number].view_name;
	 r = scm_makfrom0str(name.c_str());
      }
   return r;
}
#endif	/* USE_GUILE */

#ifdef USE_GUILE
SCM view_description(int view_number) { 

   SCM r = SCM_BOOL_F;
   if (view_number >= 0)
      if (view_number < graphics_info_t::views->size()) {
	 std::string d = (*graphics_info_t::views)[view_number].description;
	 if (d != "") {
	    r = scm_makfrom0str(d.c_str());
	 }
      }
   return r;
}
#endif	/* USE_GUILE */


/*! \brief save views to view_file_name */
void save_views(const char *view_file_name) {

   unsigned int n_views = graphics_info_t::views->size();
   if (n_views > 0) {
      std::ofstream f(view_file_name);
      if (! f) {
	 std::cout << "Cannot open view output file"
		   << view_file_name << std::endl;
      } else {
	 f << "; Views\n";
	 for (unsigned int i=0; i<n_views; i++) {
	    f << (*graphics_info_t::views)[i];
	 }
	 std::string s = "Views written to file ";
	 s += view_file_name;
	 add_status_bar_text(s.c_str());
	 // "ching!" sound here.
      }
   } else {					
      std::cout << "no views to save" << std::endl;
   }
}

// return the view number
int add_action_view(const char *view_name, const char *action_function) {
   std::string name(view_name);
   std::string func(action_function);
   coot::view_info_t::view_info_t view(name, func);  // an action view
   graphics_info_t::views->push_back(view);
   return (graphics_info_t::views->size() -1);
}

// return the view number of the new view
int insert_action_view_after_view(int view_number, const char *view_name, const char *action_function) {
   int r = -1; 
   std::string name(view_name);
   std::string func(action_function);
   coot::view_info_t::view_info_t view(name, func);  // an action view
   int n_views = graphics_info_t::views->size(); 
   if (view_number >= n_views) {
      graphics_info_t::views->push_back(view);
      r = (graphics_info_t::views->size() -1);
   } else {
      // insert a view
      std::vector <coot::view_info_t> new_views;
      for (unsigned int iview=0; iview<n_views; iview++) {
	 new_views.push_back((*graphics_info_t::views)[iview]);
	 if (iview == view_number)
	    new_views.push_back(view);
      }
      *graphics_info_t::views = new_views; // bleuch.
      r = view_number + 1;
   }
   return r;
}

void add_view_description(int view_number, const char *descr) {

   if (view_number <= graphics_info_t::views->size())
      if (view_number >= 0)
	 (*graphics_info_t::views)[view_number].add_description(descr);

}

#ifdef USE_GUILE
void go_to_view(SCM view) {
   
   SCM len_view_scm = scm_length(view);
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)      
   int len_view = scm_to_int(len_view_scm);
#else    
   int len_view = gh_scm2int(len_view_scm);
#endif
   
   if (len_view == 4) { 

      graphics_info_t g;
      int nsteps = 2000;
      if (graphics_info_t::views_play_speed > 0.000000001)
	 nsteps = int(2000.0/graphics_info_t::views_play_speed);
   
      // What is the current view:
      // 
      std::string name("Current Position");
      float quat[4];
      for (int i=0; i<4; i++)
	 quat[i] = graphics_info_t::quat[i];
      coot::Cartesian rc = g.RotationCentre();
      float zoom = graphics_info_t::zoom;
      coot::view_info_t view_c(quat, rc, zoom, name);


      // view_target is where we want to go
      float quat_target[4];
      SCM quat_scm = scm_list_ref(view, SCM_MAKINUM(0));
      SCM len_quat_scm = scm_length(quat_scm);
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
      int len_quat = scm_to_int(len_quat_scm);
#else
      int len_quat = gh_scm2int(len_quat_scm);
#endif       
      if (len_quat == 4) { 
	 SCM q0_scm = scm_list_ref(quat_scm, SCM_MAKINUM(0));
	 SCM q1_scm = scm_list_ref(quat_scm, SCM_MAKINUM(1));
	 SCM q2_scm = scm_list_ref(quat_scm, SCM_MAKINUM(2));
	 SCM q3_scm = scm_list_ref(quat_scm, SCM_MAKINUM(3));
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
	 quat_target[0] = scm_to_double(q0_scm);
	 quat_target[1] = scm_to_double(q1_scm);
	 quat_target[2] = scm_to_double(q2_scm);
	 quat_target[3] = scm_to_double(q3_scm);
#else
	 quat_target[0] = gh_scm2double(q0_scm);
	 quat_target[1] = gh_scm2double(q1_scm);
	 quat_target[2] = gh_scm2double(q2_scm);
	 quat_target[3] = gh_scm2double(q3_scm);
#endif 	 

	 SCM rc_target_scm = scm_list_ref(view, SCM_MAKINUM(1));
	 SCM len_rc_target_scm = scm_length(rc_target_scm);
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
	 int len_rc_target = scm_to_int(len_rc_target_scm);
#else	 
	 int len_rc_target = gh_scm2int(len_rc_target_scm);
#endif 	 
	 if (len_rc_target == 3) {

	    SCM centre_x = scm_list_ref(rc_target_scm, SCM_MAKINUM(0));
	    SCM centre_y = scm_list_ref(rc_target_scm, SCM_MAKINUM(1));
	    SCM centre_z = scm_list_ref(rc_target_scm, SCM_MAKINUM(2));

#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
	    double x = scm_to_double(centre_x);
	    double y = scm_to_double(centre_y);
	    double z = scm_to_double(centre_z);
#else
	    double x = gh_scm2double(centre_x);
	    double y = gh_scm2double(centre_y);
	    double z = gh_scm2double(centre_z);
#endif 	    
	    coot::Cartesian rc_target(x,y,z);

	    SCM target_zoom_scm = scm_list_ref(view, SCM_MAKINUM(2));
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
	    double zoom_target = scm_to_double(target_zoom_scm);
#else 	    
	    double zoom_target = gh_scm2double(target_zoom_scm);
#endif 	    

	    SCM name_target_scm = scm_list_ref(view, SCM_MAKINUM(3));
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
	    std::string name_target = scm_to_locale_string(name_target_scm);
#else
	    std::string name_target = SCM_STRING_CHARS(name_target_scm);
#endif
	       
	    coot::view_info_t view_target(quat_target, rc_target, zoom_target, name_target);
	    
	    // do the animation
	    coot::view_info_t::interpolate(view_c, view_target, nsteps);
	 } else {
	    std::cout << "WARNING:: bad centre in view" << std::endl;
	 }
      } else {
	 std::cout << "WARNING:: bad quat in view" << std::endl;
      } 
   }
} 
#endif // USE_GUILE


int add_spin_view(const char *view_name, int n_steps, float degrees_total) {

   graphics_info_t g;
   coot::view_info_t v(view_name, n_steps, degrees_total);
   graphics_info_t::views->push_back(v);
   std::string cmd = "add-spin-view";
   std::vector<coot::command_arg_t> args;
   args.push_back(view_name);
   args.push_back(n_steps);
   args.push_back(degrees_total);
   add_to_history_typed(cmd, args);
   return (graphics_info_t::views->size() -1);
}

void set_views_play_speed(float f) {
   graphics_info_t::views_play_speed = f;
   std::string cmd = "set-views-play-speed";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);

} 

float views_play_speed() {
   return graphics_info_t::views_play_speed;
   add_to_history_simple("views-play-speed");
}


/*  ----------------------------------------------------------------------- */
/*                  remote control                                          */
/*  ----------------------------------------------------------------------- */

void set_socket_string_waiting(const char *s) {

   // wait for lock:
   while (graphics_info_t::socket_string_waiting_mutex_lock != 0) {
      std::cout << "Waiting for lock! "
		<< graphics_info_t::socket_string_waiting_mutex_lock
		<< std::endl;
      usleep(1000000);
   }
   
   std::cout << " =============== setting mutex lock =========" << std::endl;
   // 
   // (This mutex lock *and* waiting flag may be overly complex now
   // that we simply use g_idle_add())
   graphics_info_t::socket_string_waiting_mutex_lock = 1;
   graphics_info_t::socket_string_waiting = s;
   graphics_info_t::have_socket_string_waiting_flag = 1;

   GSourceFunc f = graphics_info_t::process_socket_string_waiting_bool;
   g_idle_add(f, NULL); // if f returns FALSE then f is not called again.

   // old way, generates a Xlib async error sometimes?   
//       gtk_widget_queue_draw_area(graphics_info_t::glarea, 0, 0,
// 				 graphics_info_t::glarea->allocation.width,
// 				 graphics_info_t::glarea->allocation.height);
      
//    std::cout << "INFO:: ---- set_socket_string_waiting set to :"
// 	     << graphics_info_t::socket_string_waiting
// 	     << ":" << std::endl;

// another old way:
    //   gint return_val;
   //   GdkEventExpose event;
   //    gtk_signal_emit_by_name(GTK_OBJECT(graphics_info_t::glarea),
   //                           "configure_event",
   // 			   &event, &return_val);
   
}
