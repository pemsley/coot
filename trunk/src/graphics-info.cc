/* src/graphics-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
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
//#include "coot-sysdep.h"
#endif

#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"

#include "clipper/core/map_utils.h" // Map_stats
#include "graphical_skel.h"

#include "coot-sysdep.h"

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

#include "manipulation-modes.hh"
#include "guile-fixups.h"

#ifdef USE_PYTHON
// #include "Python.h" included above now.
#include "cc-interface.hh"
#endif

// A few non-class members - should be somewhere else, I guess.
// 
void initialize_graphics_molecules() { 
  graphics_info_t g;
  g.initialize_molecules(); 
}

// return a vector of the current valid map molecules
std::vector<int>
graphics_info_t::valid_map_molecules() const {
   
   std::vector<int> v;
   for (unsigned int i=0; i<molecules.size(); i++)
      if (is_valid_map_molecule(i))
	 v.push_back(i);
   return v;
}


// return the new molecule number
// static
int graphics_info_t::create_molecule() { 
   int imol = molecules.size();
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





GdkColor colour_by_distortion(float dist) {

   GdkColor col;

   col.pixel = 1;
   col.blue  = 0;

   if (dist < 0.0) { 
      // black for negative numbers
      col.red   = 0;
      col.green = 0;
   } else {
      if (dist < 1.4 /* was 2.0 before Tickle-fix */) { 
	 col.red   = 0;
	 col.green = 55535;
      } else {
	 if (dist < 2.2 /* was 5.0 */ ) {
	    col.red   = 55000;
	    col.green = 55000;
	    // col.blue  = 22000;
      } else {
	    if (dist < 3.0 /* was 8.0 */ ) {
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
   float scale = 10.0; 

   col.pixel = 1;
   col.blue  = 0;

   if (plot_value < -15.0*scale) { 
      col.red   = 0;
      col.green = 55535;
   } else {
      if (plot_value < -13.0*scale) {
	 col.red   = 55000;
	 col.green = 55000;
	 // col.blue  = 22000;
      } else {
	 if (plot_value < -10.0*scale) {
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

  short int old_flag = graphics_info_t::do_anti_aliasing_flag;
  graphics_info_t::do_anti_aliasing_flag = state;
  if (do_anti_aliasing_flag != old_flag) {
    draw_anti_aliasing();
  }
}

bool
graphics_info_t::background_is_black_p() const {

   bool v = 0;
   if (background_colour[0] < 0.1)
      if (background_colour[1] < 0.1)
	 if (background_colour[2] < 0.1)
	    v = 1;

   return v;
} 

void
graphics_info_t::draw_anti_aliasing() {

   // Bernie code?
   // 
   // PE 20090426 JB reports a crash here (I'm guessing it's here -
   // top of the stack was the callback
   // on_background_white1_activate()).  Can't reproduce crash.
   //
   // I did however change to using background_is_black_p().
  
  // first for stereo
  if (glarea_2) {
    if (make_current_gl_context(glarea_2)) { 
      if (do_anti_aliasing_flag) {
	// should we also add a (quality) hint here?
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	if (background_is_black_p()) {
	  glBlendFunc(GL_SRC_ALPHA,GL_ZERO);
	} else {
	  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // Thanks Stuart McN.
	}
      } else {
	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // Thanks Stuart McN.
      }
    }
  }
  // normal path
  if (glarea) {
    if (make_current_gl_context(glarea)) {
      if (do_anti_aliasing_flag) {
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	if (background_is_black_p()) {
	  glBlendFunc(GL_SRC_ALPHA,GL_ZERO);
	} else {
	  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // Thanks Stuart McN.
	}
      } else {
	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // Thanks Stuart McN.
      }
      graphics_draw();
    }
  }
}

int
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
	 s += " atoms/links in restraints from ";
	 s += cif_dictionary_filename;
	 display_density_level_screen_string = s;
	 add_status_bar_text(s);
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

   for (unsigned int i=0; i<molecules.size(); i++) {
      if (is_valid_model_molecule(i)) {
	 molecules[i].make_bonds_type_checked();
      }
   }
   return nbonds;
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

   if (0) {  // Felix test/play code to orient the residue up the
	     // screen on moving to next residue.
       
      GL_matrix m;
      clipper::Mat33<double> mat_in = m.to_clipper_mat();
      clipper::Mat33<double> mat = coot::util::residue_orientation(atom->residue, mat_in);
      coot::util::quaternion q(mat.inverse());
      quat[0] = q.q0;
      quat[1] = q.q1;
      quat[2] = q.q2;
      quat[3] = q.q3;
   }
   

   update_ramachandran_plot_point_maybe(imol, atom);

   if (environment_show_distances) {
      mol_no_for_environment_distances = imol;
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

// update the green square, where we are.
void
graphics_info_t::update_ramachandran_plot_point_maybe(int imol, CAtom *atom) {

   coot::residue_spec_t r(atom->residue);
   update_ramachandran_plot_point_maybe(imol, r);
}

void
graphics_info_t::update_ramachandran_plot_point_maybe(int imol, const coot::residue_spec_t &res_spec) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
   if (w) {
      coot::rama_plot *plot = static_cast<coot::rama_plot *>
	 (gtk_object_get_user_data(GTK_OBJECT(w)));

      plot->big_square(res_spec.chain, res_spec.resno, res_spec.insertion_code);
   } 
#endif // HAVE_GTK_CANVAS      

}

// called from accept_moving_atoms()
void 
graphics_info_t::update_ramachandran_plot_point_maybe(int imol, atom_selection_container_t moving_atoms) {

   // get the centre residue from moving atoms when 1 or 3 residue are
   // being refined and call update_ramachandran_plot_point_maybe()
   // with the residue spec.
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
   if (aa_spec_pair.first) {
      coot::residue_spec_t r(aa_spec_pair.second.second);
      update_ramachandran_plot_point_maybe(imol, r);
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
// Ignore molecules that are not displayed.
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
	 if (molecules[imol].is_displayed_p()) { 
	    dist_info = molecules[imol].nearest_atom(rc);
	    if (dist_info.first < dist_min) { 
	       imol_close = imol;
	       index_close = dist_info.second;
	       dist_min = dist_info.first;
	    }
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

   smooth_scroll_maybe_sinusoidal_acceleration(x,y,z,do_zoom_and_move_flag, target_zoom);
   // smooth_scroll_maybe_stepped_acceleration(x,y,z,do_zoom_and_move_flag, target_zoom);

}

void
graphics_info_t::smooth_scroll_maybe_sinusoidal_acceleration(float x, float y, float z,
							     short int do_zoom_and_move_flag,
							     float target_zoom) {

   // This is more like how PyMol does it (and is better than stepped
   // acceleration).

   // acceleration between istep 0 and smooth_scroll_steps (n_steps) is:
   //
   // acc = sin(istep/nsteps * 2 pi)
   //
   // v   = -cos(istep/nsteps * pi)

   float xd = x - rotation_centre_x;
   float yd = y - rotation_centre_y;
   float zd = z - rotation_centre_z;
   if ( (xd*xd + yd*yd + zd*zd) < smooth_scroll_limit*smooth_scroll_limit ) {
      

      float pre_zoom = zoom;

      float frac = 1;
      if (smooth_scroll_steps > 0)
	 frac = 1/float (smooth_scroll_steps);
      float stepping_x = frac*xd;
      float stepping_y = frac*yd;
      float stepping_z = frac*zd;
   
      float rc_x_start = rotation_centre_x;
      float rc_y_start = rotation_centre_y;
      float rc_z_start = rotation_centre_z;

      smooth_scroll_on = 1; // flag to stop wirecube being drawn.
      double v_acc = 0; // accumulated distance
      for (int istep=0; istep<smooth_scroll_steps; istep++) {
	 if (do_zoom_and_move_flag)
	    zoom = pre_zoom + float(istep+1)*frac*(target_zoom - pre_zoom);
	 double theta = 2 * M_PI * frac * istep;
	 double v = (1-cos(theta))*frac;
	 v_acc += v;
	 rotation_centre_x = rc_x_start + v_acc * xd;
	 rotation_centre_y = rc_y_start + v_acc * yd;
	 rotation_centre_z = rc_z_start + v_acc * zd;
	 graphics_draw();
      }
   
      smooth_scroll_on = 0;
   }
}

void
graphics_info_t::smooth_scroll_maybe_stepped_acceleration(float x, float y, float z,
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
	    int n_steps = smooth_scroll_steps;
	    if (do_zoom_and_move_flag)
	       n_steps *= 3;
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
      if (molecules[i].has_map()) { 
	 if (molecules[i].is_displayed_p()) {
	    is.push_back(i);
	 }
      }
   }
   return is;
} 

   
void
graphics_info_t::setRotationCentre(coot::Cartesian centre) {

   old_rotation_centre_x = rotation_centre_x; 
   old_rotation_centre_y = rotation_centre_y; 
   old_rotation_centre_z = rotation_centre_z; 

   // std::cout << "in setRotationCentre Cartesian" << graphics_info_t::smooth_scroll << std::endl;
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
			  1, target_zoom);

   rotation_centre_x = centre.get_x();
   rotation_centre_y = centre.get_y();
   rotation_centre_z = centre.get_z();
   zoom = target_zoom;
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
      g.add_status_bar_text(s);
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

   // std::cout << "clearing intermediate object..." << std::endl;
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
   
   if (moving_atoms_asc_type == coot::NEW_COORDS_ADD) { // not used!
      molecules[imol_moving_atoms].add_coords(*moving_atoms_asc);
   } else {
      bool mzo = refinement_move_atoms_with_zero_occupancy_flag;
      if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF) { 
	 molecules[imol_moving_atoms].replace_coords(*moving_atoms_asc, 1, mzo); // doesn't dealloc moving_atoms_asc
	 update_geometry_graphs(*moving_atoms_asc, imol_moving_atoms);
      } else {
	 if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE) {
	    molecules[imol_moving_atoms].replace_coords(*moving_atoms_asc, 0, mzo);
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

   // reset the b-factor?
   if (graphics_info_t::reset_b_factor_moved_atoms_flag) {
     molecules[imol_moving_atoms].set_b_factor_atom_selection(*moving_atoms_asc, graphics_info_t::default_new_atoms_b_factor, 1);
   }
   
   if (do_probe_dots_post_refine_flag) {
      setup_for_probe_dots_on_chis_molprobity(imol_moving_atoms);
   }

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *w = coot::get_validation_graph(imol_moving_atoms, coot::RAMACHANDRAN_PLOT);
   if (w) {
      coot::rama_plot *plot = (coot::rama_plot *)
	 gtk_object_get_user_data(GTK_OBJECT(w));
      // std::cout << "updating rama plot for " << imol_moving_atoms << std::endl;
      handle_rama_plot_update(plot);
      update_ramachandran_plot_point_maybe(imol_moving_atoms, *moving_atoms_asc);
   }
#endif // HAVE_GTK_CANVAS || HAVE_GNOME_CANVAS

   clear_up_moving_atoms();
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


   int mode = MOVINGATOMS;
   run_post_manipulation_hook(imol_moving_atoms, mode);
}

void
graphics_info_t::run_post_manipulation_hook(int imol, int mode) {
   
#if defined USE_GUILE && !defined WINDOWS_MINGW
   run_post_manipulation_hook_scm(imol, mode);
#endif // GUILE
#ifdef USE_PYTHON
   run_post_manipulation_hook_py(imol, mode);
#endif

}

#ifdef USE_GUILE
void
graphics_info_t::run_post_manipulation_hook_scm(int imol,
						int mode) {

   std::string pms = "post-manipulation-hook";
   SCM v = safe_scheme_command(pms);

   if (scm_is_true(scm_procedure_p(v))) {
      std::string ss = "(";
      ss += pms;
      ss += " ";
      ss += int_to_string(imol);
      ss += " ";
      ss += int_to_string(mode);
      ss += ")";
      SCM res = safe_scheme_command(ss);
      SCM dest = SCM_BOOL_F;
      SCM mess =  scm_makfrom0str("result: ~s\n");
      SCM p = scm_simple_format(dest, mess, scm_list_1(res));
      std::cout << scm_to_locale_string(p);
   }
}
#endif


#ifdef USE_PYTHON
void
graphics_info_t::run_post_manipulation_hook_py(int imol, int mode) {

   // BL says:: we can do it all in python API or use the 'lazy' method
   // and check in the python layer (which we will do...)
   PyObject *v;
   int ret;
   std::string pms = "post_manipulation_script";
   std::string check_pms = "callable(";
   check_pms += pms;
   check_pms += ")";
   v = safe_python_command_with_return(check_pms);
   ret = PyInt_AsLong(v);
   if (ret == 1) {
     std::string ss = pms;
     ss += "(";
     ss += int_to_string(imol);
     ss += ", ";
     ss += int_to_string(mode);
     ss += ")";
     PyObject *res = safe_python_command_with_return(ss);
     PyObject *fmt =  PyString_FromString("result: \%s");
     PyObject *tuple = PyTuple_New(1);
     PyTuple_SetItem(tuple, 0, res);
     //PyString_Format(p, tuple);
     PyObject *msg = PyString_Format(fmt, tuple);
     
     std::cout << PyString_AsString(msg)<<std::endl;;
     Py_DECREF(msg);
   }
   Py_XDECREF(v);
}
#endif


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
	    g.mol_no_for_environment_distances = imol_moving_atoms;
	    g.update_environment_distances_maybe(atom_index, imol_moving_atoms);
	 }
      }
   }
}


void 
graphics_info_t::clear_up_moving_atoms() { 

   std::cout << "INFO:: graphics_info_t::clear_up_moving_atoms..." << std::endl;
   moving_atoms_asc_type = coot::NEW_COORDS_UNSET; // unset
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

   if (do_torsion_restraints) {
      if (use_only_extra_torsion_restraints_for_torsions_flag) { 
	 flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      } else {
	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      }
   }

   if (do_rama_restraints)
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
   
   if (do_torsion_restraints && do_rama_restraints) { 
      if (use_only_extra_torsion_restraints_for_torsions_flag) { 
	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
      } else {
	 flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
      }
   }
	    

   // print_initial_chi_squareds_flag is 1 the first time then we turn it off.
   int steps_per_frame = dragged_refinement_steps_per_frame;
   if (! g.last_restraints.include_map_terms())
      steps_per_frame *= 3;
   // std::cout << "steps_per_frame " << steps_per_frame << std::endl;
   graphics_info_t::saved_dragged_refinement_results = 
      g.last_restraints.minimize(flags, steps_per_frame, print_initial_chi_squareds_flag);
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

      if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS) { 
	 Bond_lines_container bonds;
	 bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, 1.0, 4.7);
	 // std::cout << "done CA bonds" << std::endl;
	 regularize_object_bonds_box.clear_up();
	 regularize_object_bonds_box = bonds.make_graphical_bonds();
      } else { 
	 Bond_lines_container bonds;
	 bonds.do_Ca_bonds(*moving_atoms_asc, 1.0, 4.7);
	 // std::cout << "done CA bonds" << std::endl;
	 regularize_object_bonds_box.clear_up();
	 regularize_object_bonds_box = bonds.make_graphical_bonds();
      }
      
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
	 case BLUE_BOND:
	    glColor3f (0.60, 0.6, 0.99);
	    break;
	 case RED_BOND:
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
   if (is_valid_model_molecule(mol_no_for_environment_distances)) {
      if (g.molecules[mol_no_for_environment_distances].is_displayed_p()) { 
	 g.environment_graphics_object_internal(environment_object_bonds_box);
	 if (g.show_symmetry)
	    g.environment_graphics_object_internal(symmetry_environment_object_bonds_box);
      } 
   } 
}

// static
void
graphics_info_t::picked_intermediate_atom_graphics_object() {

   if (flash_intermediate_atom_pick_flag) { 
      glPointSize(12.0);
      glColor3f(0.99, 0.99, 0.2);
      // for some reason, I have to add the point twice and use
      // GL_POINTS, GL_POINT with one vertex does not work.
      glBegin(GL_POINTS);
      glVertex3f(intermediate_flash_point.x(),
		 intermediate_flash_point.y(),
		 intermediate_flash_point.z());
      glVertex3f(intermediate_flash_point.x(),
		 intermediate_flash_point.y(),
		 intermediate_flash_point.z());
      glEnd();
   }
} 

// This is the GL rendering of the environment bonds box
// 
void
graphics_info_t::environment_graphics_object_internal(const graphical_bonds_container &env_bonds_box) const {

   if (! display_environment_graphics_object_as_solid_flag) {
      environment_graphics_object_internal_lines(env_bonds_box); // GL lines
   } else {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glEnable(GL_LIGHT1);
      environment_graphics_object_internal_tubes(env_bonds_box); // GL cylinders and disks
      glDisable(GL_LIGHTING);
   } 
}

// This is the GL rendering of the environment bonds box
// 
void
graphics_info_t::environment_graphics_object_internal_lines(const graphical_bonds_container &env_bonds_box) const {
   if (environment_show_distances == 1) {

      if (env_bonds_box.num_colours > 0) {

	 Lines_list ll;
	 coot::Cartesian text_pos;
	 float dist;

	 float dark_bg_cor = 0.0;
	 if (! background_is_black_p())
	    dark_bg_cor = 0.29;
	 
	 glEnable(GL_LINE_STIPPLE);
	 glLineStipple (1, 0x00FF);
	 glLineWidth(2.0);
	 for (int i=0; i< env_bonds_box.num_colours; i++) {

	    bool display_these_distances_flag = 1;
	    if (i==0)
	       if (!environment_distances_show_bumps)
		  display_these_distances_flag = 0;
	    if (i==1)
	       if (!environment_distances_show_h_bonds)
		  display_these_distances_flag = 0;

	    if (display_these_distances_flag) { 
	       ll = env_bonds_box.bonds_[i]; // lightweight
	       float it = float(i);
	       if (it > 1.0) 
		  it = 1.0;

	       // now we want to draw out our bonds in various colour,
	       // according to if they have a carbon or not.
	       // 
	       glColor3f (0.8-dark_bg_cor, 0.8-0.4*it-dark_bg_cor, 0.4+0.5*it-dark_bg_cor);
	    
	       for (int j=0; j< env_bonds_box.bonds_[i].num_lines; j++) {
	   
		  const coot::CartesianPair &pair = ll.pair_list[j];
	    
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
	 }
	 glDisable(GL_LINE_STIPPLE);
      }
   }
}

// This is the GL rendering of the environment bonds box
// 
void
graphics_info_t::environment_graphics_object_internal_tubes(const graphical_bonds_container &env_bonds_box) const {

   if (environment_show_distances == 1) {

      if (env_bonds_box.num_colours > 0) {
	 
	 coot::Cartesian text_pos;
	 float dist;

	 float dark_bg_cor = 0.0;
	 if (! background_is_black_p())
	    dark_bg_cor = 0.29;
	 
	 glEnable(GL_COLOR_MATERIAL);

	 for (int i=0; i< env_bonds_box.num_colours; i++) {

	    bool display_these_distances_flag = 1;
	    if (i==0)
	       if (!environment_distances_show_bumps)
		  display_these_distances_flag = 0;
	    if (i==1)
	       if (!environment_distances_show_h_bonds)
		  display_these_distances_flag = 0;

	    if (display_these_distances_flag) { 
	       Lines_list ll = env_bonds_box.bonds_[i]; // lightweight
	       float it = float(i);
	       if (it > 1.0) 
		  it = 1.0;

	       // now we want to draw out our bonds in various colour,
	       // according to if they have a carbon or not.
	       // 
	       glColor3f (0.8-dark_bg_cor, 0.8-0.4*it-dark_bg_cor, 0.4+0.5*it-dark_bg_cor);
	    
	       for (int j=0; j< env_bonds_box.bonds_[i].num_lines; j++) {
		  const coot::CartesianPair &pair = ll.pair_list[j];

		  int n_parts = 15;
		  for (unsigned int ipart=0; ipart<n_parts; ipart++) {
		     if (coot::util::even_p(ipart)) { 
			environment_graphics_object_internal_tube(pair, ipart, n_parts);
		     }
		  }

		  // the distance text
		  text_pos = pair.getFinish().mid_point(pair.getStart()) +
		     coot::Cartesian(0.0, 0.2, 0.2);
		  glRasterPos3f(text_pos.x(), text_pos.y(), text_pos.z());
		  dist = (pair.getStart() - pair.getFinish()).amplitude();
		  glDisable(GL_LIGHTING);
		  printString(float_to_string(dist));
		  glEnable(GL_LIGHTING);
	       }
	    }
	 }
      }
   }
}

void
graphics_info_t::environment_graphics_object_internal_tube(const coot::CartesianPair &pair,
							   int ipart, int n_parts) const {
   
   coot::Cartesian bond_vec = pair.getFinish() - pair.getStart();
   coot::Cartesian bond_frag = bond_vec * (1.0/double(n_parts));
   coot::Cartesian base_point = pair.getStart() + (bond_frag * float(ipart));
   double radius = 0.04;

   graphics_object_internal_single_tube(base_point, base_point + bond_frag,
					radius, coot::FLAT_ENDS);
}

void
graphics_info_t::graphics_object_internal_single_tube(const coot::Cartesian &base_point,
						      const coot::Cartesian &end_point,
						      const double &radius,
						      const coot::tube_end_t &end_type) const {
   

   double top =  radius;
   double base = radius;
   int slices  = 12;
   int stacks  = 2;
   
   glPushMatrix();
	       
   coot::Cartesian bond_frag = end_point - base_point;
   double height = bond_frag.amplitude();


   glTranslatef(base_point.x(), base_point.y(), base_point.z());


   // 	    This code from ccp4mg's cprimitive.cc (but modified)
   //  	    ----- 
   double ax;
   double rx = 0; 
   double ry = 0;
   double length = height;
   double vz = bond_frag.z();
	 
   bool rot_x = false;
   if(fabs(vz)>1e-7){
      ax = 180.0/M_PI*acos(vz/length);
      if(vz<0.0) ax = -ax;
      rx = -bond_frag.y()*vz;
      ry = bond_frag.x()*vz;
   }else{
      double vx = bond_frag.x();
      double vy = bond_frag.y();
      ax = 180.0/M_PI*acos(vx/length);
      if(vy<0) ax = -ax;
      rot_x = true;
   }
	 
   if (rot_x) { 
      glRotated(90.0, 0.0, 1.0, 0.0);
      glRotated(ax,  -1.0, 0.0, 0.0);
   } else {
      glRotated(ax, rx, ry, 0.0);
   }
   // 	    --------

   GLUquadric* quad = gluNewQuadric();
   gluCylinder(quad, base, top, height, slices, stacks);

   if (end_type == coot::FLAT_ENDS) { 
      glScalef(1.0, 1.0, -1.0);
      gluDisk(quad, 0, base, slices, 2);
      glScalef(1.0, 1.0, -1.0);
      glTranslated(0,0,height);
      gluDisk(quad, 0, base, slices, 2);
   }

   if (end_type == coot::ROUND_ENDS) {
      GLUquadric* sphere_quad = gluNewQuadric();
      int sphere_slices = 10;
      int sphere_stacks = 10;
      gluSphere(sphere_quad, top, sphere_slices, sphere_stacks);
      glTranslated(0,0,height);
      gluSphere(sphere_quad, top, sphere_slices, sphere_stacks);
      gluDeleteQuadric(sphere_quad);
   } 
   
   gluDeleteQuadric(quad);
   glPopMatrix();
} 


void
graphics_info_t::update_environment_graphics_object(int atom_index, int imol) {

   environment_object_bonds_box =
      molecules[imol].make_environment_bonds_box(atom_index, geom_p);
}

void
graphics_info_t::update_symmetry_environment_graphics_object(int atom_index, int imol) {

   symmetry_environment_object_bonds_box =
      molecules[imol].make_symmetry_environment_bonds_box(atom_index, geom_p);
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
   display_density_level_screen_string += float_to_string_using_dec_pl(dlevel, 3);
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

   // std::pair<clipper::Coord_orth, clipper::Coord_orth> p(p1, p2);
   coot::simple_distance_object_t p(geometry_atom_index_1_mol_no, p1,
				    geometry_atom_index_2_mol_no, p2);
   distance_object_vec->push_back(p);
   graphics_draw();

   std::cout << "        distance: " << dist << " Angstroems" << std::endl;
   std::string s = "Distance: ";
   s += float_to_string(dist);
   s += " A";
   add_status_bar_text(s);
}

void 
graphics_info_t::display_geometry_distance_symm(int imol1, const coot::Cartesian &p1,
						int imol2, const coot::Cartesian &p2) { 

   coot::simple_distance_object_t p(imol1, clipper::Coord_orth(p1.x(), p1.y(), p1.z()),
				    imol2, clipper::Coord_orth(p2.x(), p2.y(), p2.z()));
   distance_object_vec->push_back(p);
   graphics_draw();
   std::cout << "Distance: " << (p1-p2).length() << std::endl;

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
   add_status_bar_text(display_density_level_screen_string);
   // redraw is in calling function for angles.
}


// ---- simple torsion ------
bool
graphics_info_t::set_angle_tors(int imol,
				const coot::atom_spec_t &as1,
				const coot::atom_spec_t &as2,
				const coot::atom_spec_t &as3,
				const coot::atom_spec_t &as4) {

   bool r = 0;
   if (is_valid_model_molecule(imol)) { 
      CAtom *at1 = molecules[imol].get_atom(as1);
      CAtom *at2 = molecules[imol].get_atom(as2);
      CAtom *at3 = molecules[imol].get_atom(as3);
      CAtom *at4 = molecules[imol].get_atom(as4);
      if (! at1)
	 std::cout << "   WARNING:: atom not found in molecule #"
		   << imol << " " << as1 << std::endl;
      if (! at2)
	 std::cout << "   WARNING:: atom not found in molecule #"
		   << imol << " " << as2 << std::endl;
      if (! at3)
	 std::cout << "   WARNING:: atom not found in molecule #"
		   << imol << " " << as3 << std::endl;
      if (! at4)
	 std::cout << "   WARNING:: atom not found in molecule #"
		   << imol << " " << as4 << std::endl;
      
      if (at1 && at2 && at3 && at4) {
	 angle_tor_pos_1 = coot::Cartesian(at1->x, at1->y, at1->z);
	 angle_tor_pos_2 = coot::Cartesian(at2->x, at2->y, at2->z);
	 angle_tor_pos_3 = coot::Cartesian(at3->x, at3->y, at3->z);
	 angle_tor_pos_4 = coot::Cartesian(at4->x, at4->y, at4->z);
	 r = 1;
      }
   }
   return r;
} 

void
graphics_info_t::display_geometry_torsion() const {

   double torsion = get_geometry_torsion(); 

   display_density_level_this_image = 1;
   display_density_level_screen_string = "  Torsion:  ";
   display_density_level_screen_string += float_to_string(torsion);
   display_density_level_screen_string += " degrees";
   add_status_bar_text(display_density_level_screen_string);
   graphics_draw();
   
}

double
graphics_info_t::get_geometry_torsion() const { 

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

   return torsion;

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
       std::cout << "DEBUG:: in find_atom_in_moving_atoms: here are the "
 		<< nSelAtoms << " qualifying atoms..." << std::endl;
       for(int i=0; i<nSelAtoms; i++)
 	 std::cout << "      " << i << "  " << local_SelAtom[i] << std::endl;
      moving_atoms_asc->mol->DeleteSelection(SelHnd);
   } else { 
     std::cout << "WARNING:: OOps: moving_atoms_asc->mol is NULL" << std::endl;
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



// 
int
graphics_info_t::create_pointer_atom_molecule_maybe() const {

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

// return the state whether to really show the baton.
bool
graphics_info_t::start_baton_here() {

   baton_root = RotationCentre();
   
   int imol_for_skel = imol_for_skeleton(); // if unset, sets and
					    // returns if only one
					    // map, else return -1
   if (imol_for_skel < 0) {

      std::cout << "WARNING: no skeleton found " << std::endl;

      std::vector<int> map_molecules = valid_map_molecules();

      if (map_molecules.size() > 0) {
	 GtkWidget *w = wrapped_create_skeleton_dialog(1);
	 gtk_widget_show(w);
	 return 0;

      } else {

	 // 20091218 It is as it was - No map.
	 // 
	 GtkWidget *w = create_baton_mode_make_skeleton_dialog();
	 int *imol_copy = new int;
	 *imol_copy = imol_for_skel;
	 gtk_object_set_user_data(GTK_OBJECT(w), (char *)imol_copy);
	 gtk_widget_show(w);
	 return 0;
      }

   } else {

      molecules[imol_for_skel].fill_skeleton_treenodemap(); // filled only if not filled.
      clipper::Coord_grid dummy_cg;
      short int use_dummy_cg = 0;
      baton_next_directions(imol_for_skel, NULL, baton_root, dummy_cg, use_dummy_cg);
      baton_next_ca_options_index = 0; 
      baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
      return 1;
   }
}

// Fill baton_next_ca_options, used in accept_baton_position() and
// several other functions.  
// 
// Generate the "Atom Guide Points" molecule (based on
// baton_next_ca_options) and display it.
// 
void
graphics_info_t::baton_next_directions(int imol_for_skel, const CAtom *latest_atom,
				       const coot::Cartesian &baton_root,
				       const clipper::Coord_grid &cg_start,
				       short int use_cg_start) {


//    std::cout << "DEBUG in baton_next_directions imol_for_skel is "
// 	     << imol_for_skel << std::endl;
         
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
   if (baton_tmp_atoms_to_new_molecule) {
      create_molecule_and_display(*baton_next_ca_options, molname);
   } else {
      update_molecule_to(*baton_next_ca_options, molname);
   }
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

   
   std::cout << "INFO:: Creating a molecule for Baton Atoms" << std::endl;
   // not found, let's create one:
   CMMDBManager *MMDBManager = new CMMDBManager();
   

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
   if (imol_for_skel >= 0) {
      // CMMDBCryst *cryst = new CMMDBCryst;
      MMDBManager->SetCell(molecules[imol_for_skel].xskel_cowtan.cell().descr().a(),
			   molecules[imol_for_skel].xskel_cowtan.cell().descr().b(),
			   molecules[imol_for_skel].xskel_cowtan.cell().descr().c(),
			   clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().alpha()),
			   clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().beta()),
			   clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().gamma()), 1);

      std::string spacegroup = molecules[imol_for_skel].xskel_cowtan.spacegroup().symbol_hm();

      std::cout << "setting spacegroup of Baton Atoms to be: " << spacegroup << std::endl;
      std::cout << "setting cell of Baton Atoms to be: "
		<< molecules[imol_for_skel].xskel_cowtan.cell().format() << std::endl;
      
      int istat_spgr = MMDBManager->SetSpaceGroup((char *)spacegroup.c_str()); // bleugh!
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
						  baton_build_chain_id,
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
      std::cout << "DEBUG:: accept_baton_position: " << baton_next_ca_options->size() << " "
		<< baton_next_ca_options_index << std::endl;
      if (baton_next_ca_options->size() > 0) { 
	 clipper::Coord_grid cg = (*baton_next_ca_options)[baton_next_ca_options_index].near_grid_pos;
	 baton_next_directions(imol_for_skel, baton_atom, baton_tip, cg, use_cg); // old tip
      } else {
	 clipper::Coord_grid cg;
	 use_cg = 0;
	 baton_next_directions(imol_for_skel, baton_atom, baton_tip, cg, use_cg);
      } 
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
	    tip_pos = non_skeleton_tip_pos();
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

coot::Cartesian
graphics_info_t::non_skeleton_tip_pos() const {

   double l=1.56;
   coot::Cartesian new_tip_pos = baton_root + coot::Cartesian(l,l,l);
   // consider using baton_previous_ca_positions, if/when it get set properly.
   return new_tip_pos;
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

   if (pos_position.size() > 0) { 
      if (is_valid_model_molecule(imol)) {
	 graphics_info_t::molecules[imol].update_molecule_to(pos_position);
      } else { 
	 create_molecule_and_display(pos_position, molname);
      }
   } else {
      std::cout << "WARNING:: No atoms guide points in update_molecule_to."
		<< "  Not updating guide points molecule" << std::endl;
   }
}

// return -1 on no such map.
int
graphics_info_t::lookup_molecule_name(const std::string &molname) const {
   
   for (int imol=0; imol<n_molecules(); imol++) {
      if (is_valid_map_molecule(imol) || (is_valid_model_molecule(imol))) {
	 if (0)
	    std::cout << "comparing map names:\n     :"
		      << graphics_info_t::molecules[imol].name_ << ":\n   "
		      << "  :" << molname << ":" << std::endl;
	 if (graphics_info_t::molecules[imol].name_ == molname) {
	    return imol;
	 }
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
int
graphics_info_t::apply_undo() {

   int state = 0;
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

	 std::string cwd = coot::util::current_working_dir();
	 if (molecules[umol].Have_modifications_p()) { 
	    if (molecules[umol].is_displayed_p()) { 
	       state = molecules[umol].apply_undo(cwd);
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
		  // now update the geometry graphs, so get the asc
		  atom_selection_container_t u_asc = molecules[umol].atom_sel;
		  update_geometry_graphs(u_asc, umol);
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
	    undo_molecule = -1; // reset it
	    if (use_graphics_interface_flag) { 
	       std::cout << "WARNING:: !!!  Changing the molecule to which "
			 << "\"Undo\"s are done." << std::endl;
	       std::string s = "WARNING:: Changing to Undo molecule";
	       add_status_bar_text(s);
	    }
	    apply_undo();       // find another molecule to undo
	 }
      }
   }

   // and now tinker with the Redo button to make it active
   //
   activate_redo_button();  // has protection for --no-graphics
   return state;
}

int 
graphics_info_t::apply_redo() { 

   int state = 0;
   
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
	    std::string cwd = coot::util::current_working_dir();
	    state = molecules[umol].apply_redo(cwd);
	    graphics_draw();
	    
	    // need to update the atom and residue list in Go To Atom widget
	    // (maybe)
	    update_go_to_atom_window_on_changed_mol(umol);
	 } else {
	    // std::cout << "DEBUG:: not applying redo" << std::endl;
	 }
      }
   }
   return state;
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

   // set the static variable for the alt conf of this atom: used when
   // we actually move the atoms.
   chi_angle_alt_conf = molecules[imol].atom_sel.atom_selection[atom_index]->altLoc;

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
	 std::cout << "This residue does not have chi angles (GLY/ALA)." << std::endl;
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
      std::string s = "This residue does not have assigned torsions/chi angles.\n";
      s += "Missing dictionary, perhaps?\n";
      info_dialog(s); // checks use_graphics_interface_flag
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


void
graphics_info_t::setup_flash_bond_internal(int i_torsion_index) {

   // turn it off first, only enable it if we find a pair:
   draw_chi_angle_flash_bond_flag = 0; // member data item

   std::cout << "flash bond i_torsion_index: " << i_torsion_index << std::endl;

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
		  bool add_reverse_contacts = 0;
		  
		  std::pair<std::string, std::string> atom_names;

		  std::pair<short int, coot::dictionary_residue_restraints_t> r =
		     geom_p->get_monomer_restraints(residue_type);

		  if (r.first) { 
		     std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
			r.second.get_non_const_torsions(find_hydrogen_torsions_flag);

		     if (i_torsion_index >= 0 && i_torsion_index < torsion_restraints.size()) {

			atom_names.first  = torsion_restraints[i_torsion_index].atom_id_2_4c();
			atom_names.second = torsion_restraints[i_torsion_index].atom_id_3_4c();
		  
			if ((atom_names.first != "") &&
			    (atom_names.second != "")) {
		     
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




// ---------------------- geometry objects -----------------------------
void
graphics_info_t::geometry_objects() {

   // 20090715 We change the type of distance_object_vec, and attach a
   // molecule from which the distance was made.  Don't display the
   // distance if the molecule corresponding to the start or end point
   // is not displayed.

   int ndist = distance_object_vec->size();
   double dist;
   clipper::Coord_orth text_pos;

   if (ndist > 0) {
      glEnable(GL_LINE_STIPPLE);
      glLineStipple (1, 0x00FF);
      glLineWidth(2.0);
      glColor3f(0.5, 0.8, 0.6);

      for (int i=0; i<ndist; i++) {
	 if (is_valid_model_molecule((*distance_object_vec)[i].imol_start)) { 
	    if (is_valid_model_molecule((*distance_object_vec)[i].imol_end)) {
	       if (molecules[(*distance_object_vec)[i].imol_start].is_displayed_p()) { 
		  if (molecules[(*distance_object_vec)[i].imol_end].is_displayed_p()) { 
		     
		     glBegin(GL_LINES);
		     glVertex3d( (*distance_object_vec)[i].start_pos.x(),
				 (*distance_object_vec)[i].start_pos.y(),
				 (*distance_object_vec)[i].start_pos.z());
		     glVertex3d( (*distance_object_vec)[i].end_pos.x(),
				 (*distance_object_vec)[i].end_pos.y(),
				 (*distance_object_vec)[i].end_pos.z());
		     text_pos = (*distance_object_vec)[i].start_pos + 
			0.5 * ( (*distance_object_vec)[i].end_pos -
				(*distance_object_vec)[i].start_pos + 
				clipper::Coord_orth(0.0, 0.1, 0.1));
		     glEnd();
		     glRasterPos3d(text_pos.x(), text_pos.y(), text_pos.z());
		     dist = clipper::Coord_orth::length( (*distance_object_vec)[i].start_pos,
							 (*distance_object_vec)[i].end_pos);
		     printString(float_to_string(dist));
		  }
	       }
	    }
	 }
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

// update_pointer_distances() you might say
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
	       if (molecules[imol].atom_selection_is_pickable()) { 
		  mol_distances = molecules[imol].distances_to_point(cen,
								     pointer_min_dist,
								     pointer_max_dist);
		  if (mol_distances.size() > 0) {
		     // append
		     for (unsigned int id=0; id<mol_distances.size(); id++)
			distances.push_back(mol_distances[id]);
		  }
	       }
	    }
	 }
      }
   
      pointer_distances_object_vec->clear();
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

   distance_object_vec->clear();
   graphics_draw();
}

void
graphics_info_t::clear_last_simple_distance() {

   int n = distance_object_vec->size();
   if (n > 0) {
      // distance_object_vec->resize(n-1); old style.  Can't do with
      // with new simple_distance_object_t
      std::vector<coot::simple_distance_object_t>::iterator it =
	 distance_object_vec->end();
	 distance_object_vec->erase(it);
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
   graphics_info_t g;
   if (! g.display_generic_objects_as_solid_flag) 
      g.draw_generic_objects_simple();
   else 
      g.draw_generic_objects_solid(); // gluCylinders and gluDisks
}



// static
void
graphics_info_t::draw_generic_objects_simple() {

   // std::cout << "debug:: drawing " << generic_objects_p->size() << " generic objects" << std::endl;
   for (unsigned int i=0; i<generic_objects_p->size(); i++) {

//       std::cout << "debug:: drawing generic object - outer " << i << std::endl;
      if ((*generic_objects_p)[i].is_displayed_flag) {

// 	 std::cout << "debug:: drawing generic  " << (*generic_objects_p)[i].lines_set.size()
// 		   << " lines " << std::endl;

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

// static
void
graphics_info_t::draw_generic_objects_solid() {

   graphics_info_t g;
   double radius = 0.05;
   
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT1);
   glEnable(GL_LIGHT0);
   glEnable(GL_COLOR_MATERIAL);

   for (unsigned int i=0; i<generic_objects_p->size(); i++) {

      if ((*generic_objects_p)[i].is_displayed_flag) {

	 // Lines
	 for (unsigned int ils=0; ils< (*generic_objects_p)[i].lines_set.size(); ils++) {
	    glLineWidth((*generic_objects_p)[i].lines_set[ils].width);
	    glColor3f((*generic_objects_p)[i].lines_set[ils].colour.red,
		      (*generic_objects_p)[i].lines_set[ils].colour.green,
		      (*generic_objects_p)[i].lines_set[ils].colour.blue);
	    unsigned int s = (*generic_objects_p)[i].lines_set[ils].lines.size();
	    for (unsigned int iline=0; iline<s; iline++) {

 	       g.graphics_object_internal_single_tube((*generic_objects_p)[i].lines_set[ils].lines[iline].coords.first,
 						      (*generic_objects_p)[i].lines_set[ils].lines[iline].coords.second,
 						      radius, coot::ROUND_ENDS);
	    }
	 }
	 
	 // Points
	 for (unsigned int ips=0; ips<(*generic_objects_p)[i].points_set.size(); ips++) {
	    glColor3f((*generic_objects_p)[i].points_set[ips].colour.red,
		      (*generic_objects_p)[i].points_set[ips].colour.green,
		      (*generic_objects_p)[i].points_set[ips].colour.blue);
	    
	    unsigned int npoints = (*generic_objects_p)[i].points_set[ips].points.size();
	    for (unsigned int ipoint=0; ipoint<npoints; ipoint++) {
	       int sphere_slices = 5;
	       int sphere_stacks = 5;
	       GLUquadric* sphere_quad = gluNewQuadric();
	       glPushMatrix();
	       glTranslatef((*generic_objects_p)[i].points_set[ips].points[ipoint].x(),
			    (*generic_objects_p)[i].points_set[ips].points[ipoint].y(),
			    (*generic_objects_p)[i].points_set[ips].points[ipoint].z());	 
	       gluSphere(sphere_quad, radius, sphere_slices, sphere_stacks);
	       gluDeleteQuadric(sphere_quad);
	       glPopMatrix();	 
	    }
	 }
      }
   }
   glDisable(GL_LIGHTING);
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

   if (c.length() == 7) {
      if (c[0] == '#') {
	 return coot::colour_holder(c); // hex colour string
      } 
   } 

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
	    colour.red   = 0.05; 
	    colour.green = 0.8; 
	    colour.blue  = 0.05;
	 } else {
	    if (c == "greentint") {
	       colour.red = 0.45; 
	       colour.green = 0.63; 
	       colour.blue = 0.45;
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
   fill_bond_parameters_internals(w, pos); // pos is imol

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
   std::cout << "debug:: safe running :" << scheme_command << ":" << std::endl; 
   std::string thunk("(lambda() "); 
   thunk += scheme_command; 
   thunk += " )";
   SCM scm_thunk = SCM_BOOL_F;

   // try/catch does not make flow control come back here when bad thunk.
   // 
   scm_thunk = scm_c_eval_string(thunk.c_str());
   SCM v = scm_catch(SCM_BOOL_T, scm_thunk, handler);

   SCM dest = SCM_BOOL_F;
   SCM mess = scm_makfrom0str("scm_catch returns: ~s\n");
   SCM sf = scm_simple_format(dest, mess, scm_list_1(v));
   std::string bad_str = scm_to_locale_string(sf);

   return v;
}
#endif // USE_GUILE

void graphics_info_t::run_user_defined_click_func() {
   
#if defined USE_GUILE && ! defined WINDOWS_MINGW

   if (scm_is_true(scm_procedure_p(user_defined_click_scm_func))) { 
      SCM arg_list = SCM_EOL;
      for (unsigned int i=0; i<user_defined_atom_pick_specs.size(); i++) {
	 SCM spec_scm = atom_spec_to_scm(user_defined_atom_pick_specs[i]);
	 SCM spec_with_model_num = scm_cons(SCM_MAKINUM(user_defined_atom_pick_specs[i].model_number),
					    spec_scm);
	 arg_list = scm_cons(spec_with_model_num, arg_list);
      } 
      arg_list = scm_reverse(arg_list);
      
      // what are we running? Print it out.
      SCM dest = SCM_BOOL_F;
      SCM mess = scm_makfrom0str("~s");
      SCM ds = scm_simple_format(dest, mess, scm_list_1(user_defined_click_scm_func));
      SCM da = scm_simple_format(dest, mess, scm_list_1(arg_list));
      std::cout << "INFO applying " << scm_to_locale_string(ds) << " on "
		<< scm_to_locale_string(da) << std::endl;
      
      SCM rest = SCM_EOL;
      SCM v = scm_apply_1(user_defined_click_scm_func, arg_list, rest);
   }

#endif // USE_GUILE

#ifdef USE_PYTHON

   if (user_defined_click_py_func) {

      if (!PyCallable_Check(user_defined_click_py_func)) {
	 std::cout<<"(PYTHON) ERROR:: user_defined_click function must be callable, is "
		  << user_defined_click_py_func->ob_type->tp_name<<std::endl;
      } else {
	 // what are we running? Print it out.
	 std::cout << "INFO applying > " 
		   << PyEval_GetFuncName(user_defined_click_py_func)
		   << " < on ";

	 PyObject *arg_list_py = PyTuple_New(user_defined_atom_pick_specs.size());
	 for (unsigned int i=0; i<user_defined_atom_pick_specs.size(); i++) {
	    PyObject *spec_py = atom_spec_to_py(user_defined_atom_pick_specs[i]);
	    // we need to add the model number too
	    PyObject *model_number_py = PyInt_FromLong(user_defined_atom_pick_specs[i].model_number);
	    PyList_Insert(spec_py, 0, model_number_py);
        
	    // continue output from above
	    PyObject *fmt = PyString_FromString("[%i,%i,'%s',%i,'%s','%s','%s']");
	    PyObject *msg = PyString_Format(fmt, PyList_AsTuple(spec_py));
	    std::cout <<PyString_AsString(msg) << " ";
	    PyTuple_SetItem(arg_list_py, i, spec_py);
	    Py_DECREF(fmt);
	    Py_DECREF(msg);
	 }
	 std::cout <<std::endl; // end ouput
      
	 if (PyTuple_Check(arg_list_py)) {
	    if (!PyCallable_Check(user_defined_click_py_func)) {
	       std::cout << "WARNING:: python user click function should have been callable." << std::endl;
	       std::cout << "WARNING:: Ignoring it." << std::endl;
	       return;
	    }
	    PyObject *result = PyEval_CallObject(user_defined_click_py_func, arg_list_py);
	    Py_DECREF(arg_list_py);
	    if (result) {
	       Py_DECREF(result);
	    }
	 } else {
	    Py_DECREF(arg_list_py);
	    std::cout<<"ERROR:: executing user_defined_click" <<std::endl;
	 }      
      }
   }

#endif // USE_PYTHON
} 


#ifdef USE_GUILE
SCM
graphics_info_t::atom_spec_to_scm(const coot::atom_spec_t &spec) const {

   SCM r = SCM_EOL;
   r = scm_cons(scm_makfrom0str(spec.alt_conf.c_str()), r);
   r = scm_cons(scm_makfrom0str(spec.atom_name.c_str()), r);
   r = scm_cons(scm_makfrom0str(spec.insertion_code.c_str()), r);
   r = scm_cons(SCM_MAKINUM(spec.resno), r);
   r = scm_cons(scm_makfrom0str(spec.chain.c_str()), r);
   r = scm_cons(SCM_MAKINUM(spec.int_user_data), r);

   return r;
} 
#endif

#ifdef USE_PYTHON
// lets have it as a tuple not a list
PyObject *
graphics_info_t::atom_spec_to_py(const coot::atom_spec_t &spec) const {

  //  PyObject *r = PyTuple_New(6);
  PyObject *r = PyList_New(6);
  PyList_SetItem(r, 0, PyInt_FromLong(spec.int_user_data));
  PyList_SetItem(r, 1, PyString_FromString(spec.chain.c_str()));
  PyList_SetItem(r, 2, PyInt_FromLong(spec.resno));
  PyList_SetItem(r, 3, PyString_FromString(spec.insertion_code.c_str()));
  PyList_SetItem(r, 4, PyString_FromString(spec.atom_name.c_str()));
  PyList_SetItem(r, 5, PyString_FromString(spec.alt_conf.c_str()));

  return r;
} 
#endif

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
std::string
graphics_info_t::backslash_filename(const std::string &s) { // needed for windows?

   std::string r = s;

   for (unsigned int i=0; i<s.size(); i++) {
      if (s[i] == '/')
	 r[i] = '\\';
   }
   return r;
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
	 add_status_bar_text(s);
      } else {
	 std::string s("Not enough atoms to find plane");
	 std::cout << s << "\n";
	 add_status_bar_text(s);
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
      add_status_bar_text(s);
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

std::ostream&
coot::operator<<(std::ostream &f, coot::view_info_t &view) {

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
