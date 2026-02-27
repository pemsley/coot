/* src/c-interface.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009 by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// $Id: c-interface.cc 1458 2007-01-26 20:20:18Z emsley $
// $LastChangedDate: 2007-01-26 20:20:18 +0000 (Fri, 26 Jan 2007) $
// $Rev: 1458 $

// Load the head if it hasn't been included.
#include <cerrno>
#include <cstddef>
#include <exception>
#include <stdexcept>
#include <utility>
#include "coords/phenix-geo.hh"
#include "glib.h"
#include "gtk/gtk.h"
#include "gtk/gtkshortcut.h"
#ifdef USE_PYTHON
#ifndef PYTHONH
#define PYTHONH
#include <Python.h>
#include "python-3-interface.hh"
#endif
#endif

#include "compat/coot-sysdep.h"

#include <stdlib.h>
#include <string.h> // strncpy, strncmp
#include <iostream>
#include <fstream>
#include <algorithm>

#if !defined(_MSC_VER)
#include <glob.h> // for globbing.  Needed here?
#endif

#ifdef USE_GUILE
#include <libguile.h>
#include "c-interface-scm.hh"
#include "guile-fixups.h"
#endif // USE_GUILE

#ifdef USE_GUILE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#endif

#ifdef USE_PYTHON
#include "c-interface-python.hh"
#endif // USE_PYTHON

#include "compat/sleep-fixups.h"
#include "coot-utils/gl-matrix.h"

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
#define snprintf _snprintf
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

// #include <GL/glut.h> // needed for glutGet(GLUT_ELAPSED_TIME);

#include "clipper/ccp4/ccp4_map_io.h"

#include "globjects.h" //includes gtk/gtk.h

#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/read-amber-trajectory.hh"
#include "coot-database.hh"
#include "coot-fileselections.h"

// #include "xmap-interface.h"
#include "graphics-info.h"

#include "skeleton/BuildCas.h"

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "c-interface-ligands.hh"

#include "coot-version.hh"

#include "widget-headers.hh"
#include "widget-from-builder.hh"

#include "testing.hh"


#include "positioned-widgets.h"

// moving column_label selection to c-interface from mtz bits.
#include "cmtz-interface.hh"
// #include "mtz-bits.h" stuff from here moved to cmtz-interface

#include "read-molecule.hh" // now with std::string args

#include "widget-from-builder.hh"
#include "gtk-manual.hh"
#include "glarea_tick_function.hh"

#include "validation-graphs/sequence-view-widget.hh"

#include "utils/logging.hh"
extern logging logger;

// This is (already) in git-revision-count.cc
//
int svn_revision() {
    return git_revision_count();
}

std::string coot_version() {

   std::string s = VERSION;
   return s;
}

std::string coot_version_extra_info() {

   std::string version_string; //  was = VERSION;

   version_string += "(revision-count ";
   version_string += coot::util::int_to_string(git_revision_count());
   version_string += ")\n";

#ifdef USE_GUILE
   version_string += "[with guile ";
   version_string += coot::util::int_to_string(SCM_MAJOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(SCM_MINOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(SCM_MICRO_VERSION);
   version_string += " embedded]\n";
#endif
#ifdef USE_PYTHON
   version_string += "[with python ";
   version_string += coot::util::int_to_string(PY_MAJOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(PY_MINOR_VERSION);
   version_string += ".";
   version_string += coot::util::int_to_string(PY_MICRO_VERSION);
   version_string += " embedded]\n";
#endif
   std::string s = "COOT_BUILD_INFO_STRING"; // FIXME
   s = COOT_BUILD_INFO_STRING;
   if (! s.empty()) {
      version_string += "Builder_info: ";
      version_string += s;
      version_string += "\n";
   }

   s = git_commit();
   version_string += "git commit: ";
   version_string += s;
   version_string += "\n";

   if (! s.empty()) {
      std::string bt =  "COOT_SYS_BUILD_TYPE";  // FIXME
      bt = COOT_SYS_BUILD_TYPE;
      version_string += "Binary type: ";
      version_string += bt;
      version_string += "\n";
   }

   return version_string;
}

// return the coot_revision as a char *.  note that svn_revision()
// returns the svn revision as an it.  More useful?
char *coot_revision() {

   std::string revision_string = " (revision ";
   revision_string += coot::util::int_to_string(svn_revision());
   revision_string += ") ";
   int len = revision_string.length();

   char *r = new char[len+1];
   strncpy(r, revision_string.c_str(), len+1);
   return r;
}

#ifdef USE_GUILE
SCM coot_sys_build_type_scm() {

   std::string sb = COOT_SYS_BUILD_TYPE;
   std::cout << sb << std::endl;
   SCM r = scm_from_locale_string(sb.c_str());
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *coot_sys_build_type_py() {

   std::string sb = COOT_SYS_BUILD_TYPE;
   PyObject *r = myPyString_FromString(sb.c_str());
   return r;
}
#endif // USE_PYTHON

/*!  \brief tell coot that you prefer to run python scripts if/when
  there is an option to do so. */
void set_prefer_python() {

#ifdef USE_PYTHON
   graphics_info_t::prefer_python = 1;
#endif // USE_PYTHON
}

/*! \brief the python-prefered mode.

This is available so that the scripting functions know whether on not
to put themselves onto in as menu items.

return 1 for python is prefered, 0 for not. */
int prefer_python() {
   return graphics_info_t::prefer_python;
}




/*  -------------------------------------------------------------------- */
/*                     Testing Interface:                                */
/*  -------------------------------------------------------------------- */

#ifdef USE_GUILE
SCM test_internal_scm() {

   SCM r = SCM_BOOL_T;

#ifdef BUILT_IN_TESTING
   int status = greg_internal_tests();
   if (!status)
      r = SCM_BOOL_F;
   else
      status = greg_tests_using_external_data();
   if (!status)
      r = SCM_BOOL_F;
#endif

   return r;
}
#endif // USE_GUILE

#ifdef USE_GUILE
SCM test_internal_single_scm() {

   SCM r = SCM_BOOL_T;

#ifdef BUILT_IN_TESTING
   int status = test_internal_single();
   if (!status)
      r = SCM_BOOL_F;
#endif
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *test_internal_py() {

   PyObject *r = Py_True;

#ifdef BUILT_IN_TESTING
   int status = test_internal();
   if (!status)
      r = Py_False;
#endif

   Py_INCREF(r);
   return r;
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
PyObject *test_internal_single_py() {

   PyObject *r = Py_True;

#ifdef BUILT_IN_TESTING
   int status = test_internal_single();
   if (!status)
      r = Py_False;
#endif
   Py_INCREF(r);
   return r;
}
#endif // USE_PYTHON


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

#ifdef USE_GUILE
SCM molecule_name_stub_scm(int imol, int include_path_flag) {
   std::string r;
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol))
      r = graphics_info_t::molecules[imol].name_sans_extension(include_path_flag);
   return scm_from_locale_string(r.c_str());
}
#endif

#ifdef USE_PYTHON
PyObject *molecule_name_stub_py(int imol, int include_path_flag) {
   std::string r;
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol))
      r = graphics_info_t::molecules[imol].name_sans_extension(include_path_flag);
   return myPyString_FromString(r.c_str());
}
#endif

void set_molecule_name(int imol, const char *new_name) {

   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].set_name(new_name);
   }
}

void set_show_graphics_ligand_view(int state) {
   graphics_info_t::show_graphics_ligand_view_flag = state;
   graphics_draw();
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

short int python_at_prompt_at_startup_state() {
   return graphics_info_t::python_at_prompt_flag;
}


#include "startup-scripts.hh"

bool run_startup_scripts_state() {
   return graphics_info_t::run_startup_scripts_flag;
}




// Is this a repeat of something?  I don't know (doesn't look like it)
//
void run_generic_script(const std::vector<std::string> &cmd_strings) {

   graphics_info_t g;

#if defined(USE_GUILE) && !defined(WINDOWS_MINGW)
   std::string s = g.state_command(cmd_strings, coot::STATE_SCM);
   safe_scheme_command(s);
#else
#ifdef USE_PYTHON
   std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
   safe_python_command(s);
#endif // USE_PYTHON
#endif // USE_GUILE

   std::string cmd = "run-generic-script";
   std::vector<coot::command_arg_t> args;
   for(unsigned int i=0; i<cmd_strings.size(); i++)
      args.push_back(clipper::String(cmd_strings[i]));
   add_to_history_typed(cmd, args);
}



//  20220528-PE from globjects.cc
// ################################ put this in globjects-new.cc ? #############################

// amount is not in degrees, it is in fractions of a circle, e.g. 10/360.
//
std::vector<float> rotate_rgb(std::vector<float> &rgb, float amount) {

#if 0
   if (true) { // print a rotation to colour table
      int n_cols = 100; // either side
      for (int i = -n_cols; i<n_cols; i++) {
         float rotation_size = 0.01f * static_cast<float>(i);
         std::vector<float> orig_colours = { 0.0f,  0.8f, 0.0f };
         std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
         std::cout << "debug colours::" << rgb_new[0] << " " << rgb_new[1] << " " << rgb_new[2]
                   << " using rotation_size " << rotation_size << std::endl;
      }

      // Result:
      //
      // if we start at solid green then rotation_size for "no rotation" is 0.0
      //                                 rotation_size for full rotation is -0.33 (solid red)
   }
#endif

   std::vector<float> hsv = coot::convert_rgb_to_hsv(rgb);

   // add 20 degrees to hue (or whatever)
   // std::cout << "adding " << amount << " to hue " << std::endl;
   hsv[0] += amount;
   while  (hsv[0] > 1.0) {
      hsv[0] -= 1.0;
   }

   std::vector<float> r = coot::convert_hsv_to_rgb(hsv);
   //     std::cout << "rotate from ("
   // << rgb[0] << " " << rgb[1] << " " << rgb[2] << ")\n"
   //  	     << "         to ("
   // << rgb[0] << " " << rgb[1] << " " << rgb[2] << ")\n";
   return r;
   // return convert_hsv_to_rgb(hsv);

}



// ################################ put this in globjects-new.cc ? #############################

//  20220528-PE from globjects.cc
double
get_idle_function_rock_target_angle() {
#if 0
   graphics_info_t g;

   long t = 0; // glutGet(GLUT_ELAPSED_TIME);
   long delta_t = t - g.time_holder_for_rocking;
   double rock_sf = 0.001 * g.idle_function_rock_freq_scale_factor;

   double theta = delta_t * rock_sf;

   while (theta > 2 * M_PI)
      theta -= 2*M_PI;
   while (theta < -2 * M_PI)
      theta += 2*M_PI;

   double target_angle = g.idle_function_rock_amplitude_scale_factor *
      0.015 * 2 * M_PI * sin(theta);

   return target_angle;
#endif
   return 0.0;
}


//
// Return the molecule number of the molecule that we just filled.
// Return -1 if there was a failure.
//
int handle_read_draw_molecule(const std::string &filename) {

   int r = graphics_info_t::recentre_on_read_pdb;
   return handle_read_draw_molecule_with_recentre(filename, r);
}

void
update_things_on_move_and_redraw() {

   graphics_info_t g;
   g.update_things_on_move_and_redraw();
}

int make_updating_model_molecule(const char *filename) {

   int status = 1;

   int imol = handle_read_draw_molecule_with_recentre(filename, 0);

   if (is_valid_model_molecule(imol)) {
      updating_coordinates_molecule_parameters_t *ucp =
	 new updating_coordinates_molecule_parameters_t(filename);
      graphics_info_t::molecules[imol].continue_watching_coordinates_file = true;
      GSourceFunc f = GSourceFunc(graphics_info_t::molecules[imol].watch_coordinates_file);
      guint updating_mol_timeout_idx = g_timeout_add(500, f, ucp);
      graphics_info_t::molecules[imol].continue_watching_mtz = true;
   } else {
      status = 0;
   }
   return status;
}

void show_calculate_updating_maps_pythonic_gui() {

   std::cout << "debug:: in show_calculate_updating_maps_gui()" << std::endl;

   // 20230430-PE Don't use Python.

#ifdef USE_PYTHON
   std::string cmd = "import coot_gui ; coot_gui.show_updating_maps_chooser()";
   std::cout << "debug:: in show_calculate_updating_maps_gui() calling safe_python_command() cmd " << cmd << std::endl;
   safe_python_command(cmd);
#endif

}


// do we need to return the imol for the model and the map? If so, this
// needs to be in cc-interface.hh
//
void updating_refmac_refinement_files(const char *updating_refmac_refinement_files_json_file_name) {

#ifdef MAKE_UPDATING_REFMAC_REFINEMENT_MOLECULES // using JSON

   if (updating_refmac_refinement_files_json_file_name) {

      // set a timeout function that looks for this file
      GSourceFunc f = GSourceFunc(updating_refmac_refinement_json_timeout_function);
      std::string *fn_p = new std::string(updating_refmac_refinement_files_json_file_name);
      guint updating_idx = g_timeout_add(500, f, fn_p); // possibly illegal.
   }
#else
   std::cout << "ERROR:: updating_refmac_refinement_files() is just a stub - needs CXX11"
	     << std::endl;
#endif // MAKE_UPDATING_REFMAC_REFINEMENT_MOLECULES
}

int updating_refmac_refinement_json_timeout_function(gpointer data) {

   int status = 1; // keep going

#ifdef MAKE_UPDATING_REFMAC_REFINEMENT_MOLECULES

   if (! data) return 0;

   std::string *fn = reinterpret_cast<std::string *>(data);
   std::string json_file_name = *fn;

   std::fstream f(json_file_name);
   if (f) {

      std::string s;
      f.seekg(0, std::ios::end);
      s.reserve(f.tellg());
      f.seekg(0, std::ios::beg);

      s.assign((std::istreambuf_iterator<char>(f)),
	       std::istreambuf_iterator<char>());
      unsigned int n_already_done = 0;

      std::vector<coot::mtz_to_map_info_t> mtz_map_infos;

      try {
	 json j = json::parse(s);
         json ls = j["animation"];
         unsigned int n_cycles = ls.size();

         std::cout << "n_cycles " << n_cycles << std::endl;
         std::string mtz_file_name, f_col, phi_col;
         std::string model_file_name;

         json anim_part_id    = ls["id"];
         json anim_part_map   = ls["map"];
         json anim_part_model = ls["model"];
         json anim_part_x     = ls["xmissing"]; // test
         std::cout << "map: " << anim_part_map << std::endl;
         std::cout << "x:   " << anim_part_x   << std::endl;

         if (! anim_part_id.is_null() && !anim_part_model.is_null() && !anim_part_map.is_null()) {

            std::string id = anim_part_id;
            model_file_name = anim_part_model["filepath"];
            mtz_file_name = anim_part_map["filepath"];

            json column_sets = anim_part_map["columns"];
	    coot::mtz_to_map_info_t mmi_2fofc;
	    coot::mtz_to_map_info_t mmi_fofc;
            std::size_t n_column_sets = column_sets.size();
            for (std::size_t i=0; i<n_column_sets; i++) {
               json col_set = column_sets[i];
               if (col_set["id"] == "2Fo-Fc") {
                  std::cout << "Found 2Fo-Fc" << std::endl;
                  json labels = col_set["labels"];
                  coot::mtz_to_map_info_t mmi;
                  mmi.id = id;
                  mmi.mtz_file_name = mtz_file_name;
                  mmi.f_col         = labels[0];
                  mmi.phi_col       = labels[1];
                  mtz_map_infos.push_back(mmi);
		  mmi_2fofc = mmi;
               }
               if (col_set["id"] == "Fo-Fc") {
                  std::cout << "Found Fo-Fc" << std::endl;
                  json labels = col_set["labels"];
                  coot::mtz_to_map_info_t mmi;
                  mmi.id = id;
                  mmi.mtz_file_name = mtz_file_name;
                  mmi.f_col         = labels[0];
                  mmi.phi_col       = labels[1];
                  mmi.is_difference_map = true;
                  mtz_map_infos.push_back(mmi);
		  mmi_fofc = mmi;
               }
            }
            if (true) { // debug
               std::cout << "Here with id: " << id << std::endl;
               std::cout << "Here with model_file_name: " << model_file_name << std::endl;
               std::cout << "Here with mtz_file_name: " << mtz_file_name << std::endl;
               for (unsigned int i=0; i<mtz_map_infos.size(); i++)
                  std::cout << mtz_map_infos[i].mtz_file_name << " " << mtz_map_infos[i].id << " "
                            << mtz_map_infos[i].f_col << " " << mtz_map_infos[i].phi_col << std::endl;

               // OK, but which map and which model needs to be updated?
               std::cout << "Here are the molecules so far " << std::endl;
               for (int i=0; i<graphics_n_molecules(); i++) {
                  std::cout << "molecule name " << i << " " << graphics_info_t::molecules[i].name_
			    << std::endl;
               }
            }

            int imol_model_updating = -1;
            int imol_map_1_updating = -1;
            int imol_map_2_updating = -1;
            for (int i=0; i<graphics_n_molecules(); i++) {
               std::string name = graphics_info_t::molecules[i].get_name();
               std::string id_model = id + "_model";
               std::string id_2fofc = id + "_2fofc";
               std::string id_fofc  = id + "_fofc";
               if (name == id_model) {
                  if (is_valid_model_molecule(i)) {
                     imol_model_updating = i;
                     break;
                  }
               }
               if (name == id_2fofc) {
                  if (is_valid_map_molecule(i)) {
                     imol_map_1_updating = i;
                     break;
                  }
               }
               if (name == id_fofc) {
                  if (is_valid_map_molecule(i)) {
                     imol_map_2_updating = i;
                     break;
                  }
               }
            }

            if (is_valid_model_molecule(imol_model_updating)) {
               std::string cwd = coot::util::current_working_dir();
	       graphics_info_t::molecules[imol_model_updating].update_molecule(model_file_name, cwd);
            }
            if (is_valid_map_molecule(imol_map_1_updating)) {
	       if (! mmi_2fofc.f_col.empty()) {
		  std::string cwd = coot::util::current_working_dir();
		  graphics_info_t::molecules[imol_map_1_updating].map_fill_from_mtz(mmi_2fofc.mtz_file_name,
										    cwd,
										    mmi_2fofc.f_col,
										    mmi_2fofc.phi_col,
										    mmi_2fofc.w_col,
										    mmi_2fofc.use_weights,
										    mmi_2fofc.is_difference_map,
										    graphics_info_t::map_sampling_rate);
	       }
            }
            if (is_valid_map_molecule(imol_map_2_updating)) {
	       if (! mmi_fofc.f_col.empty()) {
		  std::string cwd = coot::util::current_working_dir();
		  graphics_info_t::molecules[imol_map_1_updating].map_fill_from_mtz(mmi_fofc.mtz_file_name,
										    cwd,
										    mmi_fofc.f_col,
										    mmi_fofc.phi_col,
										    mmi_fofc.w_col,
										    mmi_fofc.use_weights,
										    mmi_fofc.is_difference_map,
										    graphics_info_t::map_sampling_rate);
	       }
            }
         }
      }
      catch(const nlohmann::detail::type_error &e) {
	 std::cout << "ERROR:: " << e.what() << std::endl;
      }
      catch(const nlohmann::detail::parse_error &e) {
	 std::cout << "ERROR:: " << e.what() << std::endl;
      }

      status = 0; // no need to run this function again
   }
#endif // MAKE_UPDATING_REFMAC_REFINEMENT_MOLECULES
   return status;
}

void allow_duplicate_sequence_numbers() {

   graphics_info_t::allow_duplseqnum = true;
}


void set_convert_to_v2_atom_names(short int state) {
   graphics_info_t::convert_to_v2_atom_names_flag = state;
}


int handle_read_draw_molecule_with_recentre(const std::string &filename,
					    int recentre_on_read_pdb_flag) {

   int r = -1;
   //

   graphics_info_t g;

   std::string cmd = "handle-read-draw-molecule-with-recentre";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(filename));
   args.push_back(recentre_on_read_pdb_flag);
   add_to_history_typed(cmd, args);

   std::string f(filename);

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

      r = read_shelx_ins_file(filename, recentre_on_read_pdb_flag);

   } else {
      // recentre and not a backup-restore
      // -1 is for failure strangely.
      int imol = g.create_molecule();
//       std::cout << "DEBUG:: c-interface handle_read_draw_molecule on molecule number: "
// 		<< imol << std::endl;
//       std::cout << " DEBUG:: created placeholder molecule number " << imol << std::endl;
      float bw = graphics_info_t::default_bond_width;
      int bonds_box_type = graphics_info_t::default_bonds_box_type;
      istat = g.molecules[imol].handle_read_draw_molecule(imol, f,
							  coot::util::current_working_dir(),
							  graphics_info_t::Geom_p(),
							  recentre_on_read_pdb_flag, 0,
							  g.allow_duplseqnum,
							  g.convert_to_v2_atom_names_flag,
							  bw, bonds_box_type, true);

      if (istat == 1) {
	 // std::cout << "Molecule " << imol << " read successfully\n";
         logger.log(log_t::INFO, "Molecule ", imol, " read successfully");

	 // we do this somewhat awkward in and out thing with the
	 // molecule, because I don't want to (or am not able to) pass a
	 // modifiable dictionary to the handle_read_draw_molecule()
	 // function.  But maybe I *should* add it as an argument... Not
	 // sure.
	 //
	 // So currently we read in the molecule, read the molecule for
	 // types not in the dictionary and then try to read an sbase
	 // description (should be fast) and if found, rerun the bonding
	 // algorithm.


	 std::vector<std::string> types_with_no_dictionary;
         const coot::protein_geometry *geom_p = g.Geom_p();
         if (g.Geom_p())
	    g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
         else
            std::cout << "ERROR:: handle_read_draw_molecule_with_recentre() Geom_p() returns null" << std::endl;

	 int first_n_types_with_no_dictionary = types_with_no_dictionary.size();

         if (false)
            std::cout << "DEBUG:: there were " << types_with_no_dictionary.size()
                      << " types with no dictionary " << std::endl;

	 for (unsigned int i=0; i<types_with_no_dictionary.size(); i++) {
	    if (false)
	       std::cout << "DEBUG:: calling try_dynamic_add: " << types_with_no_dictionary[i]
			 << " with read number " << g.cif_dictionary_read_number << std::endl;
	    g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], g.cif_dictionary_read_number);
	    g.cif_dictionary_read_number++;
	 }

         if (geom_p)
	    types_with_no_dictionary = g.molecules[imol].no_dictionary_for_residue_type_as_yet(*geom_p);

	 if (types_with_no_dictionary.size()) {
	    if (g.Geom_p()->try_load_ccp4srs_description(types_with_no_dictionary))
	       g.molecules[imol].make_bonds_type_checked();
	 } else {
	    // perhaps we have read dictionaries for everything (but
	    // first check that there had been dictionaries to read.
	    if (first_n_types_with_no_dictionary > 0) {
	       g.molecules[imol].make_bonds_type_checked();
	    }
	 }


	 if (graphics_info_t::nomenclature_errors_mode == coot::PROMPT) {
	    // Now, did that PDB file contain nomenclature errors?
	    std::vector<std::pair<std::string,coot::residue_spec_t> > nomenclature_errors;
            if (geom_p)
               nomenclature_errors = g.molecules[imol].list_nomenclature_errors(geom_p);
	    // gui function checks use_graphics_interface_flag
	    if (nomenclature_errors.size())
	       show_fix_nomenclature_errors_gui(imol, nomenclature_errors);
	 }
	 if (graphics_info_t::nomenclature_errors_mode == coot::AUTO_CORRECT) {
	    fix_nomenclature_errors(imol);
	 }



	 // if the go to atom widget exists, update its optionmenu to
	 // reflect the existance of this new molecule.
	 if (g.go_to_atom_window) {
	    //
	    // 20090620:
	    // Are we sure that we want to set the go_to_atom_molecule
	    // on reading a pdb file?
	    //
	    // The code handling the residue trees presumes that we
	    // don't.  i.e. the residue trees are not updated if the
	    // new molecule menu item is selected but has the pos as
	    // the go_to_atom_molecule() which we (used to) set here.

	    // g.set_go_to_atom_molecule(imol);  // No.

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
	 g.add_status_bar_text(s);
	 r =  imol;
      } else {
	 g.erase_last_molecule();
	 std::string s("Failed to read coordinates file ");
	 s += filename;
	 g.add_status_bar_text(s);
	 r =  -1;
      }
   }
   return r;
}

// this is declared in the wrong header - put it in c-interface.hh?
int move_molecule_to_screen_centre_internal(int imol);

/*! \brief read coordinates from filename and recentre the new
  molecule at the scren rotation centre. */
int handle_read_draw_molecule_and_move_molecule_here(const std::string &filename) {

   int imol = handle_read_draw_molecule_with_recentre(filename, 0);
   move_molecule_to_screen_centre_internal(imol);
   return imol;
}

int read_pdb(const std::string &filename) {
   // history is done in the handle function
   return handle_read_draw_molecule(filename);
}

/*! \brief read coordinates from filename */
int read_coordinates(const std::string &filename) {
   return handle_read_draw_molecule(filename);
}

int read_coordinates_as_string(const std::string &file_contents, const std::string &molecule_name) {

#if !defined _MSC_VER
   pid_t pid = getpid();
#else
   DWORD pid = GetCurrentProcessId();
#endif
   std::string pid_str = std::to_string(pid);
   std::string fn("tmp-");
   fn += pid_str;
   fn += ".pdb";
   std::ofstream f(fn);
   f << file_contents;
   f.close();
   int imol = read_coordinates(fn);
   if (is_valid_model_molecule(imol))
      set_molecule_name(imol, molecule_name.c_str());
   return imol;
}


int read_amber_trajectory(int imol_coords,
                          const std::string &trajectory_file_name,
                          int start_frame,
                          int end_frame,
                          int stride) {

   int imol = -1;
   graphics_info_t g;

   if (!is_valid_model_molecule(imol_coords)) {
      std::cout << "WARNING:: read_amber_trajectory: invalid topology molecule " << imol_coords << std::endl;
      return -1;
   }

   mmdb::Manager *topology_mol = g.molecules[imol_coords].atom_sel.mol;
   mmdb::Manager *traj_mol = coot::read_amber_trajectory(topology_mol, trajectory_file_name,
                                                         start_frame, end_frame, stride);
   if (traj_mol) {
      imol = g.create_molecule();
      atom_selection_container_t asc = make_asc(traj_mol);
      std::string name = trajectory_file_name + "_trajectory";
      g.molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
   }

   return imol;
}


//! set (or unset) GEMMI as the molecule parser. Currently by passing an int.
void set_use_gemmi_as_model_molecule_parser(int state) {

   if (state) {
#ifdef USE_GEMMI
      graphics_info_t g;
      g.set_use_gemmi(state);
#else
      std::cout << "WARNING:: this executable was not compiled with gemmi " << std::endl;
      add_status_bar_text("WARNING:: this executable was not compiled with gemmi");
#endif
   }
}




/*! \brief replace pdb.  Fail if molecule_number is not a valid model molecule.
  Return -1 on failure.  Else return molecule_number  */
int clear_and_update_model_molecule_from_file(int molecule_number,
					      const char *file_name) {
   int imol = -1;
   if (is_valid_model_molecule(molecule_number)) {
      atom_selection_container_t asc = get_atom_selection(file_name, true, graphics_info_t::allow_duplseqnum, true);
      mmdb::Manager *mol = asc.mol;
      graphics_info_t::molecules[molecule_number].replace_molecule(mol);
      imol = molecule_number;
      graphics_draw();
   }
   return imol;
}


#ifdef USE_GUILE
/*! \brief - get the name state of the input model. Return false if there was an erro with the molecule index */
SCM get_input_model_was_cif_state_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      if (g.molecules[imol].get_input_molecule_was_in_mmcif_state())
         r = SCM_BOOL_T;
   }
   return r;
}
#endif

#ifdef USE_PYTHON
/*! \brief - get the name state of the input model */
PyObject *get_input_molecule_was_in_mmcif_state_py(int imol) {

   PyObject *r = PyBool_FromLong(false);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      if (g.molecules[imol].get_input_molecule_was_in_mmcif_state())
         r = PyBool_FromLong(true);
   }
   Py_INCREF(r);
   return r;
}
#endif



void set_draw_zero_occ_markers(int status) {

   graphics_info_t g;
   g.draw_zero_occ_spots_flag = status;
   std::string cmd = "set-draw-zero-occ-markers";
   std::vector<coot::command_arg_t> args;
   args.push_back(status);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

/*! \brief set status of drawing cis-peptide markups

  default status is 1. */
void set_draw_cis_peptide_markups(int status) {

   graphics_info_t g;
   g.draw_cis_peptide_markups = status;
   std::string cmd = "set-draw-cis-peptide-markups";
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

/*! \brief molecule number of first small (<400 atoms) molecule.

return -1 on no such molecule
  */
int first_small_coords_imol() {

   int imol = -1;
   for (int i=0; i<graphics_n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 int n_atoms = graphics_info_t::molecules[i].atom_sel.n_selected_atoms;
	 if (n_atoms < N_ATOMS_MEANS_BIG_MOLECULE) {
	    imol = i;
	    break;
	 }
      }
   }
   add_to_history_simple("first-small-coords-imol");
   return imol;
}


int first_unsaved_coords_imol() {

   int imol = -1;
   for (int i=0; i<graphics_n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 if (graphics_info_t::molecules[i].Have_unsaved_changes_p()) {
	    imol = i;
	    break;
	 }
      }
   }
   add_to_history_simple("first-unsaved-coords-imol");
   return imol;
}



void hardware_stereo_mode() {

   // this should be a graphics-info function. Move me. FIXME

   if (graphics_info_t::use_graphics_interface_flag) {
      if (graphics_info_t::display_mode != coot::HARDWARE_STEREO_MODE) {
	 int previous_mode = graphics_info_t::display_mode;
	 graphics_info_t::display_mode = coot::HARDWARE_STEREO_MODE;
         GtkWidget *gl_area = graphics_info_t::glareas[0];
	 // GtkWidget *vbox = lookup_widget(gl_area, "main_window_vbox");
	 GtkWidget *vbox = widget_from_builder("main_window_vbox");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in hardware_stereo_mode!\n";
	 } else {

	    if ( (previous_mode == coot::SIDE_BY_SIDE_STEREO) ||
		 (previous_mode == coot::DTI_SIDE_BY_SIDE_STEREO)  ||
		 (previous_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ) {

               std::cout << "Do some stereo stuff" << std::endl;
            }
         }
      } else {
         std::cout << "Already in hardware stereo mode" << std::endl;
      }
   }
   add_to_history_simple("hardware-stereo-mode");
}

void zalman_stereo_mode() {

  // FIXME this is not really zalman!!!!!
   if (graphics_info_t::use_graphics_interface_flag) {
      if (graphics_info_t::display_mode != coot::HARDWARE_STEREO_MODE) {
	 int previous_mode = graphics_info_t::display_mode;
	 graphics_info_t::display_mode = coot::ZALMAN_STEREO;

	 GtkWidget *vbox = widget_from_builder("main_window_vbox");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in zalman_stereo_mode!\n";
	 } else {

	    if ( (previous_mode == coot::SIDE_BY_SIDE_STEREO) ||
		 (previous_mode == coot::DTI_SIDE_BY_SIDE_STEREO)  ||
		 (previous_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ) {

	       if (graphics_info_t::glareas.size() == 2) {
		  gtk_widget_set_visible(graphics_info_t::glareas[1], FALSE);
		  graphics_info_t::glareas[1] = NULL;
	       }
	    }

#if 0 // 20220528-PE Graphics no longer works like this
	    short int try_hardware_stereo_flag = 5;
	    GtkWidget *glarea = gl_extras(vbox, try_hardware_stereo_flag);
	    if (glarea) {
	       // std::cout << "INFO:: switch to zalman_stereo_mode succeeded\n";
	       logger.log(log_t::INFO, "switch to zalman_stereo_mode succeeded");
	       if (graphics_info_t::idle_function_spin_rock_token) {
		  toggle_idle_spin_function(); // turn it off;
	       }
	       // gtk_widget_destroy(graphics_info_t::glareas[0]);
	       gtk_widget_set_visible(graphics_info_t::glareas[0], FALSE);
	       graphics_info_t::glareas[0] = glarea;
	       gtk_widget_set_visible(glarea, TRUE);
	       // antialiasing?
	       graphics_info_t g;
	       g.draw_anti_aliasing();
	       graphics_draw();
	    } else {
	       std::cout << "WARNING:: switch to zalman_stereo_mode failed\n";
	       graphics_info_t::display_mode = previous_mode;
	    }
#endif
	 }
      } else {
	 std::cout << "Already in zalman stereo mode" << std::endl;
      }
   }
   add_to_history_simple("zalman-stereo-mode");
}

void mono_mode() {

   if (graphics_info_t::use_graphics_interface_flag) {

      if (graphics_info_t::display_mode != coot::MONO_MODE) {
	 int previous_mode = graphics_info_t::display_mode;
         // GtkWidget *gl_widget = lookup_widget(graphics_info_t::get_main_window(), "window1");
         GtkWidget *gl_widget = graphics_info_t::glareas[0];
         // int x_size = gtk_widget_get_allocated_width(gl_widget);
         // int y_size = gtk_widget_get_allocated_height(gl_widget);
         GtkAllocation allocation;
         gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
         int x_size = allocation.width;
         int y_size = allocation.height;
	 graphics_info_t::display_mode = coot::MONO_MODE;
         graphics_info_t::graphics_draw();
         GtkWidget *gl_area = graphics_info_t::glareas[0];
	 // GtkWidget *vbox = lookup_widget(gl_area, "main_window_vbox");
	 GtkWidget *vbox = widget_from_builder("main_window_vbox");
	 if (!vbox) {
	    std::cout << "ERROR:: failed to get vbox in mono mode!\n";
	 } else {
#if 0 // 20220528-PE Graphics no longer works like this
	    short int try_hardware_stereo_flag = 0;
	    GtkWidget *glarea = gl_extras(vbox, try_hardware_stereo_flag);
	    if (glarea) {
	       // std::cout << "INFO:: switch to mono_mode succeeded\n";
	       logger.log(log_t::INFO, "switch to mono_mode succeeded");
	       if (graphics_info_t::idle_function_spin_rock_token) {
		  toggle_idle_spin_function(); // turn it off;
	       }
	       gtk_widget_destroy(graphics_info_t::get_main_window());
	       if (graphics_info_t::glareas.size() == 2) {
		  gtk_widget_destroy(graphics_info_t::glareas[1]);
		  graphics_info_t::glareas[1] = NULL;
	       }
	       graphics_info_t::glareas[0] = glarea;
	       gtk_widget_set_visible(glarea, TRUE);
               // now we shall resize to half the window size if we had
               // side-by-side stereo before
               if ((previous_mode == coot::SIDE_BY_SIDE_STEREO) ||
                   (previous_mode == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ||
                   (previous_mode == coot::DTI_SIDE_BY_SIDE_STEREO)) {
                  set_graphics_window_size(x_size/2, y_size);
               }
	       // std::cout << "DEBUG:: mono_mode() update maps and draw\n";
	       graphics_info_t g;
	       g.setRotationCentre(coot::Cartesian(g.X(), g.Y(), g.Z()));
	       g.post_recentre_update_and_redraw();
	       g.draw_anti_aliasing();
	       redraw_background();
	    } else {
	       graphics_info_t::display_mode = previous_mode;
	       std::cout << "WARNING:: switch to mono mode failed\n";
	    }
#endif
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
      // generated 2 new glareas by calling gl_extras().
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
         // int previous_mode = graphics_info_t::display_mode;
         short int stereo_mode = coot::SIDE_BY_SIDE_STEREO;
         if (use_wall_eye_flag)
            stereo_mode = coot::SIDE_BY_SIDE_STEREO_WALL_EYE;
         // GtkWidget *vbox = lookup_widget(graphics_info_t::get_main_window(), "vbox1");
#if 0 // 20220528-PE Graphics no longer works like this
         GtkWidget *vbox = widget_from_builder("vbox1"); // 20220309-PE is this what you really want? - needs testing
         GtkWidget *glarea = gl_extras(vbox, stereo_mode);
         if (glarea) {
            if (graphics_info_t::idle_function_spin_rock_token) {
               toggle_idle_spin_function(); // turn it off;
            }
            gtk_widget_destroy(graphics_info_t::glareas[0]);
            graphics_info_t::glareas[0] = glarea; // glarea_2 is stored by gl_extras()
            gtk_widget_set_visible(glarea, TRUE);
            gtk_widget_set_visible(graphics_info_t::glareas[1], TRUE);
            update_maps();
            // antialiasing?
            graphics_info_t g;
            g.draw_anti_aliasing();
            redraw_background();
            graphics_draw();
         } else {
            std::cout << "WARNING:: switch to side by side mode failed!\n";
         }
#endif
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
      if (state) {
         short int stereo_mode;
         if (graphics_info_t::display_mode != coot::DTI_SIDE_BY_SIDE_STEREO) {
            // int previous_mode = graphics_info_t::display_mode;
            stereo_mode = coot::DTI_SIDE_BY_SIDE_STEREO;
         } else {
            stereo_mode = coot::SIDE_BY_SIDE_STEREO;
         }
         // GtkWidget *vbox = lookup_widget(graphics_info_t::get_main_window(), "vbox1");
#if 0 // 20220528-PE Graphics no longer works like this
         GtkWidget *vbox = widget_from_builder("vbox1"); // 20220309-PE is this right?
         GtkWidget *glarea = gl_extras(vbox, stereo_mode);
         if (graphics_info_t::use_graphics_interface_flag) {
            if (graphics_info_t::idle_function_spin_rock_token) {
               toggle_idle_spin_function(); // turn it off;
            }
            gtk_widget_destroy(graphics_info_t::glareas[0]);
            graphics_info_t::glareas[0] = glarea;
            gtk_widget_set_visible(glarea, TRUE);
            if (graphics_info_t::glareas.size() == 2)
               gtk_widget_set_visible(graphics_info_t::glareas[1], TRUE);
            graphics_draw();
         } else {
            std::cout << "WARNING:: switch to side by side mode failed!\n";
         }
#endif
      }
   }
   // add_to_history_simple("dti-side-by-side-stereo-mode");
}


int stereo_mode_state() {
   add_to_history_simple("stereo-mode-state");
   return graphics_info_t::display_mode;
}


/*! \brief set the stereo mode (the relative view of the eyes)

0 is 2010-mode
1 is modern mode
*/
void set_stereo_style(int mode) {

   if (mode == 0)
      graphics_info_t::stereo_style_2010 = true;
   else
      graphics_info_t::stereo_style_2010 = false;

   graphics_draw();
}

void set_hardware_stereo_angle_factor(float f) {
}

void set_stereo_angle(float f) {
   graphics_info_t::stereo_angle = f;
   std::string cmd = "set-stereo-angle";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

float hardware_stereo_angle_factor_state() {
   add_to_history_simple("hardware-stereo-angle-factor-state");
   return 0;
}

void set_model_display_radius(int state, float radius) {

   graphics_info_t::model_display_radius.first  = state;
   graphics_info_t::model_display_radius.second = radius;
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


void unset_sticky_sort_by_date() {

   add_to_history_simple("unset-sticky-sort-by-date");
   graphics_info_t g;
   g.sticky_sort_by_date = 0;

}


void set_filter_fileselection_filenames(int istate) {
   std::string cmd = "set-filter-fileselection-filenames";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);

   graphics_info_t::filter_fileselection_filenames_flag = istate;

}

int filter_fileselection_filenames_state() {
   add_to_history_simple("filter-fileselection-filenames-state");
   return graphics_info_t::filter_fileselection_filenames_flag;
}

/*! \brief is the given file name suitable to be read as coordinates? */
short int file_type_coords(const char *file_name) {

   graphics_info_t g;
   return g.file_type_coords(file_name);
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

/*! \brief return the new molecule number */

void get_coords_for_accession_code(const std::string &text) {

   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(text));
   scripting_function("get-ebi-pdb", args);
}





std::vector<std::string> filtered_by_glob(const std::string &pre_directory,
					  int data_type) {

   std::vector<std::string> v; // returned object
   std::vector<std::string> globs;

#if !defined(_MSC_VER)

   // std::map<std::string, int, std::less<std::string> >  files;

   if (data_type == COOT_COORDS_FILE_SELECTION)
      globs = *graphics_info_t::coordinates_glob_extensions;
   if (data_type == COOT_DATASET_FILE_SELECTION)
      globs = *graphics_info_t::data_glob_extensions;
   if (data_type == COOT_MAP_FILE_SELECTION)
      globs = *graphics_info_t::map_glob_extensions;
   if (data_type == COOT_CIF_DICTIONARY_FILE_SELECTION)
      globs = *graphics_info_t::dictionary_glob_extensions;
   if (data_type == COOT_SAVE_COORDS_FILE_SELECTION)
      globs = *graphics_info_t::coordinates_glob_extensions;
   if (data_type == COOT_PHS_COORDS_FILE_SELECTION)
      globs = *graphics_info_t::coordinates_glob_extensions;

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

#endif // MSC

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



void set_auto_read_column_labels(const char *fwt, const char *phwt,
				 int is_for_diff_map_flag) {

   coot::mtz_column_trials_info_t n(fwt, phwt, is_for_diff_map_flag);
   graphics_info_t::user_defined_auto_mtz_pairs.push_back(n);

   std::string cmd = "set-auto-read-column-labels";
   std::vector<coot::command_arg_t> args;
   args.push_back(fwt);
   args.push_back(phwt);
   args.push_back(is_for_diff_map_flag);
   add_to_history_typed(cmd, args);

}

void toggle_idle_spin_function() {

   graphics_info_t g;

   if (g.do_tick_spin)
      g.do_tick_spin = false;
   else
      g.do_tick_spin = true;

   if (g.do_tick_spin) {
      if (g.glareas[0]) {
         int new_tick_id = gtk_widget_add_tick_callback(g.glareas[0], g.glarea_tick_func, 0, 0);
         g.idle_function_spin_rock_token = new_tick_id;
      }
   }
   graphics_draw();

   add_to_history_simple("toggle-idle-function");
}


void toggle_idle_rock_function() {

   graphics_info_t g;

   if (g.do_tick_rock)
      g.do_tick_rock = false;
   else
      g.do_tick_rock = true;

   if (g.do_tick_rock) {
      g.time_holder_for_rocking = std::chrono::high_resolution_clock::now();
      if (g.glareas[0]) {
         int new_tick_id = gtk_widget_add_tick_callback(g.glareas[0], g.glarea_tick_func, 0, 0);
         g.idle_function_spin_rock_token = new_tick_id;
      }
   }
   graphics_draw();
   add_to_history_simple("toggle-idle-rock-function");
}

/*! \brief Settings for the inevitable discontents who dislike the
   default rocking rates (defaults 1 and 1)  */
void set_rocking_factors(float width, float freq_scale) {

   graphics_info_t g;
   g.idle_function_rock_amplitude_scale_factor = width;
   g.idle_function_rock_freq_scale_factor = freq_scale;

}

/* Turn on nice animated ligand interaction display */
void toggle_flev_idle_ligand_interactions() {

   graphics_info_t g;
   if (g.idle_function_ligand_interactions_token == 0) {
      set_flev_idle_ligand_interactions(1);
   } else {
      set_flev_idle_ligand_interactions(0);
   }
   add_to_history_simple("toggle-idle-ligand-interactions");
}

void set_flev_idle_ligand_interactions(int state) {
   graphics_info_t g;
   if (state == 0) {
      // turn them off if they were on
      if (g.idle_function_ligand_interactions_token) {

	 std::cout << "GTK-FIXME set_flev_idle_ligand_interactions" << std::endl;
	 // g_idle_remove(g.idle_function_ligand_interactions_token);
	 //g.idle_function_ligand_interactions_token = 0;

	 for (unsigned int imol=0; imol<g.molecules.size(); imol++) {
	    if (is_valid_model_molecule(imol)) {
	       g.molecules[imol].draw_animated_ligand_interactions_flag = 0;
	    }
	 }
      }
   } else {
      // turn them on if they were off.
      if (g.idle_function_ligand_interactions_token == 0) {
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
         std::cout << "FIXME toggle_flev_idle_ligand_interactions() timer\n";
#else
	 g.idle_function_ligand_interactions_token = gdk_threads_add_timeout(100, animate_idle_ligand_interactions, NULL); // Hmm.
	 // g.time_holder_for_ligand_interactions = glutGet(GLUT_ELAPSED_TIME);
#endif
      }
   }
   g.graphics_draw();
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
   args.push_back(single_quote(txt));
   add_to_history_typed(cmd, args);
}

// As info_dialog, but print to console too. Useful for error and
// warning messages.
//
void info_dialog_and_text(const char *txt) {

   graphics_info_t g;
   g.info_dialog_and_text(txt);
   std::string cmd = "info-dialog-and-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(txt));
   add_to_history_typed(cmd, args);
}

void info_dialog_with_markup(const char *txt) {

   graphics_info_t g;
   g.info_dialog_and_text(txt, true);
   std::string cmd = "info-dialog-and-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(txt));
   add_to_history_typed(cmd, args);
}


/*! \brief created an ephemeral label in the graphics window
 *
 * the text stays on screen for about 2 sesconds.
 *
 * @param txt the text
*/
void ephemeral_overlay_label(const char *txt) {

   graphics_info_t::ephemeral_overlay_label(std::string(txt));
}



void
set_main_window_title(const char *s) {

   graphics_info_t g;
   if (s) {
      if (g.use_graphics_interface_flag){
	 if (g.get_main_window()) {
	    GtkWidget *win = g.get_main_window();
	    if (win) {
	       std::string ss(s);
	       if (! ss.empty()) {
		  g.main_window_title = ss;
		  GtkWindow *window = GTK_WINDOW(win);
		  gtk_window_set_title(window, s);
	       }
	    }
	 }
      }
   }
}





/*  ------------------------------------------------------------------------ */
/*                         file selection                                    */
/*  ------------------------------------------------------------------------ */

// void
// set_directory_for_fileselection(GtkWidget *fileselection1) {
//    graphics_info_t g;
//    g.set_directory_for_fileselection(fileselection1);
// }

// void
// save_directory_from_fileselection(const GtkWidget *fileselection) {
//    graphics_info_t g;
//    g.save_directory_from_fileselection(fileselection);
// }

// void save_directory_for_saving_from_fileselection(const GtkWidget *fileselection) {
//    graphics_info_t g;
//    g.save_directory_for_saving_from_fileselection(fileselection);
// }

/* and the gtk2 equivalents, we dont use most of them any more but keep
   them for gtk2 move maybe  */


void
set_directory_for_filechooser(GtkWidget *fileselection1) {
   graphics_info_t g;
   g.set_directory_for_filechooser(fileselection1);
}

void
save_directory_from_filechooser(const GtkWidget *filechooser) {
   graphics_info_t g;
   g.save_directory_from_filechooser(filechooser);
}

void save_directory_for_saving_from_filechooser(const GtkWidget *fileselection) {
   graphics_info_t g;
   g.save_directory_for_saving_from_filechooser(fileselection);
}


bool compare_mtimes(coot::str_mtime a, coot::str_mtime b) {
   return a.mtime > b.mtime;
}




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

/*! \brief turn off (0) or on (1) auto recontouring (on screen centre change) (default it on) */
void set_auto_recontour_map(int state) {
   graphics_info_t::auto_recontour_map_flag = state;
}

/*! \brief turn off (0) or on (1) auto recontouring (on screen centre change) (default it on) */
int get_auto_recontour_map() {
   return graphics_info_t::auto_recontour_map_flag;
}



void add_symmetry_on_to_preferences_and_apply() {

   set_show_symmetry_master(1);

   graphics_info_t g;
   g.add_to_preferences("xenops-symmetry.scm", "(set-show-symmetry-master 1)");
   g.add_to_preferences("xenops_symmetry.py",   "coot.set_show_symmetry_master(1)");

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

/*! \brief return the density at the given point for the given map */
float density_at_point(int imol, float x, float y, float z) {

   float r = -999.9;
   if (is_valid_map_molecule(imol)) {
      clipper::Coord_orth p(x,y,z);
      r = graphics_info_t::molecules[imol].density_at_point(p);
   }
   return r;
}


#ifdef USE_GUILE
float density_score_residue_scm(int imol, SCM residue_spec_scm, int imol_map) {

   float v = 0.0;
   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_model_molecule(imol)) {
	 graphics_info_t g;
	 coot::residue_spec_t residue_spec = residue_spec_from_scm(residue_spec_scm);
	 mmdb::Residue *r = g.molecules[imol].get_residue(residue_spec);
	 if (r) {
	    mmdb::PPAtom residue_atoms = 0;
	    int n_residue_atoms;
	    r->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       float d_at = density_at_point(imol_map, at->x, at->y, at->z);
	       v += d_at * at->occupancy;
	    }
	 }
      }
   }
   return v;
}
#endif

#ifdef USE_PYTHON
float density_score_residue_py(int imol, PyObject *residue_spec_py, int imol_map) {

   float v = 0.0;
   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_model_molecule(imol)) {
	 graphics_info_t g;
	 coot::residue_spec_t residue_spec = residue_spec_from_py(residue_spec_py);
	 mmdb::Residue *r = g.molecules[imol].get_residue(residue_spec);
	 if (r) {
	    mmdb::PPAtom residue_atoms = 0;
	    int n_residue_atoms;
	    r->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       float d_at = density_at_point(imol_map, at->x, at->y, at->z);
	       v += d_at * at->occupancy;
	    }
	 }
      }
   }
   return v;
}
#endif

/*! \brief simple density score for given residue (over-ridden by scripting function) */
float density_score_residue(int imol, const char *chain_id, int res_no, const char *ins_code, int imol_map) {

   float v = 0.0;
   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_model_molecule(imol)) {
         graphics_info_t g;
         coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
         mmdb::Residue *r = g.molecules[imol].get_residue(residue_spec);
         if (r) {
            mmdb::PPAtom residue_atoms = 0;
            int n_residue_atoms;
            r->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (!at->isTer()) {
                  float d_at = density_at_point(imol_map, at->x, at->y, at->z);
                  v += d_at * at->occupancy;
                }
            }
         }
      }
   }
   return v;

}


#ifdef USE_GUILE
SCM map_mean_scm(int imol) {

  SCM r = SCM_BOOL_F;
  if (is_valid_map_molecule(imol)) {
    float s = graphics_info_t::molecules[imol].map_mean();
    r = scm_from_double(s);
  }
  return r;
}
#endif

#ifdef USE_GUILE
SCM map_sigma_scm(int imol) {

  SCM r = SCM_BOOL_F;
  if (is_valid_map_molecule(imol)) {
    float s = graphics_info_t::molecules[imol].map_sigma();
    r = scm_from_double(s);
  }
  return r;
}
#endif

#ifdef USE_PYTHON
PyObject *map_sigma_py(int imol) {

  PyObject *r = Py_False;
  if (is_valid_map_molecule(imol)) {
    float s = graphics_info_t::molecules[imol].map_sigma();
    r = PyFloat_FromDouble(s);
  }

  if (PyBool_Check(r)) {
    Py_INCREF(r);
  }
  return r;
}
#endif
#ifdef USE_PYTHON
PyObject *map_mean_py(int imol) {

  PyObject *r = Py_False;
  if (is_valid_map_molecule(imol)) {
    float s = graphics_info_t::molecules[imol].map_mean();
    r = PyFloat_FromDouble(s);
  }

  if (PyBool_Check(r)) {
    Py_INCREF(r);
  }
  return r;
}
#endif


#ifdef USE_GUILE
/*! \brief return either scheme false on non-a-map or list (mean, standard-deviation, skew, kurtosis) */
SCM map_statistics_scm(int imol) {

  SCM r = SCM_BOOL_F;
  if (is_valid_map_molecule(imol)) {
     map_statistics_t ms = graphics_info_t::molecules[imol].map_statistics();
     r = scm_list_4(scm_from_double(ms.mean),
                    scm_from_double(ms.sd),
                    scm_from_double(ms.skew),
                    scm_from_double(ms.kurtosis));
  }
  return r;
}
#endif

#ifdef USE_PYTHON
/*! \brief return either python False on not-a-map or list [mean, standard-deviation, skew, kurtosis] */
PyObject *map_statistics_py(int imol) {

   PyObject *r = Py_False;
  if (is_valid_map_molecule(imol)) {
     map_statistics_t ms = graphics_info_t::molecules[imol].map_statistics();
     r = PyList_New(4);
     PyList_SetItem(r, 0, PyFloat_FromDouble(ms.mean));
     PyList_SetItem(r, 1, PyFloat_FromDouble(ms.sd));
     PyList_SetItem(r, 2, PyFloat_FromDouble(ms.skew));
     PyList_SetItem(r, 3, PyFloat_FromDouble(ms.kurtosis));
  }
  if (PyBool_Check(r))
    Py_INCREF(r);
  return r;
}
#endif




void set_density_size_from_widget(const char *text) {

   if (text) {
      try {
         std::string ss(text);
         float f = coot::util::string_to_float(ss);
         if (f > 0.0) {
            if (f < 1999.9) {
               graphics_info_t g;
               g.box_radius_xray = f;
               for (int ii=0; ii<g.n_molecules(); ii++) {
                  if (is_valid_map_molecule(ii))
                     g.molecules[ii].update_map(true);
               }
            }
         }
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
   }
   graphics_draw();
}

void
set_density_size_em_from_widget(const char *text) {

   if (text) {
      try {
         std::string ss(text);
         float f = coot::util::string_to_float(ss);
         if (f > 0.0) {
            // example tomo is 21000x16000x3000
            if (f < 19999.9) {
               graphics_info_t g;
               g.box_radius_em = f;
               for (int ii=0; ii<g.n_molecules(); ii++) {
                  if (is_valid_map_molecule(ii))
                     g.molecules[ii].update_map(true);
               }
            } else {
               std::cout << "over the limit: " << f << std::endl;
            }
         }
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
   }
   graphics_draw();
}


void set_map_radius_em(float radius) {

   graphics_info_t g;
   g.box_radius_em = radius;
   for (int ii=0; ii<g.n_molecules(); ii++)
      g.molecules[ii].update_map(true);
   graphics_draw();
   std::string cmd = "set-map-radius-em";
   std::vector<coot::command_arg_t> args;
   args.push_back(radius);
   add_to_history_typed(cmd, args);
}

void set_density_size(float f) {

   graphics_info_t g;
   g.box_radius_xray = f;
   for (int ii=0; ii<g.n_molecules(); ii++) {
      g.molecules[ii].update_map(true);
   }
   graphics_draw();
   std::string cmd = "set-density-size";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);

}

void set_density_size_em(float f) {

   graphics_info_t g;
   g.box_radius_em = f;
   for (int ii=0; ii<g.n_molecules(); ii++) {
      g.molecules[ii].update_map(true);
   }
   graphics_draw();
   std::string cmd = "set-density-size-em";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
}

/*! \brief set the extent of the box/radius of electron density contours */
void set_map_radius(float f) {
   set_density_size(f);
}

/*! \brief return the extent of the box/radius of electron density contours */
float get_map_radius() {
  float ret = graphics_info_t::box_radius_xray;
  return ret;
}



void set_display_intro_string(const char *str) {

   if (str) {
      if (graphics_info_t::use_graphics_interface_flag) {
	 std::string s(str);
	 graphics_info_t g;
	 g.display_density_level_screen_string = s;
	 g.add_status_bar_text(s);
      }

      std::string cmd = "set-display-intro-string";
      std::vector<coot::command_arg_t> args;
      args.push_back(single_quote(str));
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



// return 1 on "yes, it has a cell".
//
int has_unit_cell_state(int imol) {

   int istate = 0;
   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules()) {
	 if (graphics_info_t::molecules[imol].has_model() ||
	     graphics_info_t::molecules[imol].has_xmap()) {
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

      std::cout << "Cannot interpret " << text << ".  Assuming 10A" << std::endl;
      g.symmetry_search_radius = 10.0;
   }
   //
   for (int ii=0; ii<g.n_molecules(); ii++) {
      g.molecules[ii].update_symmetry();
   }
   graphics_draw();

}

void set_symmetry_size(float f) {
   graphics_info_t g;
   g.symmetry_search_radius = f;
   for (int ii=0; ii<g.n_molecules(); ii++) {
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
         if (! g.glareas.empty())
            g.update_things_on_move_and_redraw();
      }
   }
   std::string cmd = "set-symmetry-whole-chain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


/*! \brief set show frame-per-second flag */
void set_show_fps(int flag) {

   graphics_info_t g;
   g.SetShowFPS(flag);
   std::string cmd = "set-fps-flag";
   std::vector<coot::command_arg_t> args;
   args.push_back(flag);
   add_to_history_typed(cmd, args);
}


void set_fps_flag(int flag) {

   set_show_fps(flag);
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

   if (is_valid_model_molecule(imol)) {
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

/*! \brief draw little coloured balls on atoms

turn off with state = 0

turn on with state = 1 */
void set_draw_stick_mode_atoms(int imol, short int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_display_stick_mode_atoms(state);
   }
   graphics_draw();
}

void set_draw_missing_residues_loops(short int state) {

   bool prev_state = graphics_info_t::draw_missing_loops_flag;
   bool new_state = state;
   if (new_state != prev_state) {
      graphics_info_t::draw_missing_loops_flag = new_state;
      for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
         if (is_valid_model_molecule(imol)) {
            graphics_info_t::molecules[imol].make_bonds_type_checked();
         }
      }
      graphics_draw();
   }
}



/*! \brief the state of draw hydrogens for molecule number imol.

return -1 on bad imol.  */
int draw_hydrogens_state(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].draw_hydrogens();
   }
   return r;
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




void
handle_symmetry_colour_change(int mol, double* col) {

   //
   graphics_info_t::symmetry_colour[0] = col[0];
   graphics_info_t::symmetry_colour[1] = col[1];
   graphics_info_t::symmetry_colour[2] = col[2];

   graphics_draw();
}

GdkRGBA
get_map_colour(int imol) {

   //
   GdkRGBA colour;

   if (imol < graphics_info_t::n_molecules()) {
      if (graphics_info_t::molecules[imol].has_xmap()) {
         colour = graphics_info_t::molecules[imol].map_colour;
         colour.red   *= 65535;
         colour.green *= 65535;
         colour.blue  *= 65535;
         colour.alpha *= 65535;
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
      double rc = graphics_info_t::molecules[imol].map_colour.red;
      double gc = graphics_info_t::molecules[imol].map_colour.green;
      double bc = graphics_info_t::molecules[imol].map_colour.blue;
      r = SCM_EOL;
      // put red at the front of the resulting list
      r = scm_cons(scm_from_double(bc), r);
      r = scm_cons(scm_from_double(gc), r);
      r = scm_cons(scm_from_double(rc), r);
   }
   std::string cmd = "map-colour-components";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r;
}
#endif

// BL says:: this is for python
#ifdef USE_PYTHON
PyObject *map_colour_components_py(int imol) {

   PyObject *r;
   r = Py_False;
   if (is_valid_map_molecule(imol)) {

      double rc = graphics_info_t::molecules[imol].map_colour.red;
      double gc = graphics_info_t::molecules[imol].map_colour.green;
      double bc = graphics_info_t::molecules[imol].map_colour.blue;
      r = PyList_New(3);
      // put red at the front of the resulting list
      PyList_SetItem(r, 0, PyFloat_FromDouble(rc));
      PyList_SetItem(r, 1, PyFloat_FromDouble(gc));
      PyList_SetItem(r, 2, PyFloat_FromDouble(bc));
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON


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


// Noddy version
void
colour_map_by_other_map(int imol_map, int imol_map_used_for_colouring) {

   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_map_molecule(imol_map_used_for_colouring)) {
         graphics_info_t g;
         clipper::Xmap<float> &xmap_for_colouring = g.molecules[imol_map_used_for_colouring].xmap;

         std::cout << "------------- colour_map_by_other_map() API calling molecules colour_map_using_map()"
                   << std::endl;
         g.molecules[imol_map].colour_map_using_map(xmap_for_colouring);
      }
   }
   graphics_draw();

}

#ifdef USE_PYTHON
//! So, if the map has 4 entries covering the range from  0 to 1, then the table_bin_size would be 0.25
//! and the colour_table list would have 4 entries covering the range 0->0.25, 0.25->0.5, 0.5->0.75, 0.75->1.0
void
colour_map_by_other_map_py(int imol_map, int imol_map_used_for_colouring, float table_bin_start, float table_bin_size,
                           PyObject *colour_table_py) {

   if (is_valid_map_molecule(imol_map)) {
      if (is_valid_map_molecule(imol_map_used_for_colouring)) {
         if (PyList_Check(colour_table_py)) {
            graphics_info_t g;
            clipper::Xmap<float> &xmap_for_colouring = g.molecules[imol_map_used_for_colouring].xmap;
            std::vector<coot::colour_t> colour_list;
            unsigned int n = PyObject_Length(colour_table_py);
            for (unsigned int i=0; i<n; i++) {
               PyObject *i_py = PyList_GetItem(colour_table_py, i);
               if (PyList_Check(i_py)) {
                  unsigned int i_l = PyObject_Length(i_py);
                  if (i_l == 3) {
                     double r = PyFloat_AsDouble(PyList_GetItem(i_py, 0));
                     double g = PyFloat_AsDouble(PyList_GetItem(i_py, 1));
                     double b = PyFloat_AsDouble(PyList_GetItem(i_py, 2));
                     coot::colour_t col(r,g,b);
                     colour_list.push_back(col);
                  }
               } else {
                  std::cout << "Not a list " << std::endl;
                  break;
               }
            }

            std::cout << "debug:: in colour_map_by_other_map_py() colour_list size " << colour_list.size() << std::endl;

            if (colour_list.size() == n) {
               // we read the table OK.
               g.molecules[imol_map].colour_map_using_map(xmap_for_colouring,
                                                          table_bin_start, table_bin_size,
                                                          colour_list);
            }
         } else {
            std::cout << "colour table was not a list " << std::endl;
         }
      }
   }
   graphics_draw();
}
#endif


void
colour_map_by_other_map_turn_off(int imol_map) {

   if (is_valid_map_molecule(imol_map)) {
      graphics_info_t::molecules[imol_map].colour_map_using_other_map_flag = false;
      std::cout << "FIXME:: make the map update" << std::endl;
   }

}


// -------------------------------------------------------------------
#ifdef USE_PYTHON
PyObject *export_molecule_as_x3d(int imol) {

   PyObject *r = PyList_New(3);
   PyList_SetItem(r, 0, PyList_New(0));
   PyList_SetItem(r, 1, PyList_New(0));
   PyList_SetItem(r, 2, PyList_New(0));
   // or model, eventually.
   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      coot::density_contour_triangles_container_t tc = g.molecules[imol].export_molecule_as_x3d();
      if (false)
         std::cout << "debug:: tc.points " << tc.points.size() << " tc.normals " << tc.normals.size()
                   << " tc-point_indices " << tc.point_indices.size() << std::endl;
      if (tc.points.size() > 0) {
         if (tc.point_indices.size() > 0) {
            PyObject *indices_list = PyList_New(3 * tc.point_indices.size());
            PyObject *vertex_list  = PyList_New(3 * tc.points.size());
            PyObject *normals_list = PyList_New(3 * tc.normals.size());

            for (unsigned int i=0; i<tc.points.size(); i++) {
               PyList_SetItem(vertex_list, i*3,   PyFloat_FromDouble(tc.points[i].x()));
               PyList_SetItem(vertex_list, i*3+1, PyFloat_FromDouble(tc.points[i].y()));
               PyList_SetItem(vertex_list, i*3+2, PyFloat_FromDouble(tc.points[i].z()));
            }
            for (unsigned int i=0; i<tc.normals.size(); i++) {
               PyList_SetItem(normals_list, i*3,   PyFloat_FromDouble(tc.normals[i].x()));
               PyList_SetItem(normals_list, i*3+1, PyFloat_FromDouble(tc.normals[i].y()));
               PyList_SetItem(normals_list, i*3+2, PyFloat_FromDouble(tc.normals[i].z()));
            }
            for (unsigned int i=0; i<tc.point_indices.size(); i++) {
               PyList_SetItem(indices_list, i*3,   PyLong_FromLong(tc.point_indices[i].pointID[0]));
               PyList_SetItem(indices_list, i*3+1, PyLong_FromLong(tc.point_indices[i].pointID[1]));
               PyList_SetItem(indices_list, i*3+2, PyLong_FromLong(tc.point_indices[i].pointID[2]));
            }

            PyList_SetItem(r, 0, indices_list);
            PyList_SetItem(r, 1, vertex_list);
            PyList_SetItem(r, 2, normals_list);
         }
      }
   }
   return r;

}
#endif

bool export_molecule_as_obj(int imol, const std::string &fn)  {

   bool status = false;
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      status = graphics_info_t::molecules[imol].export_molecule_as_obj(fn);
   }
   return status;
}

bool export_molecule_as_gltf(int imol, const std::string &file_name) {

   bool status = false;
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      status = graphics_info_t::molecules[imol].export_molecule_as_gltf(file_name);
   }
   return status;
}




// -------------------------------------------------------------------

gdouble*
get_symmetry_bonds_colour(int idummy) {

   //
   gdouble* colour;
   colour = (gdouble *) malloc(4*sizeof(gdouble));

   colour[0] = graphics_info_t::symmetry_colour[0];
   colour[1] = graphics_info_t::symmetry_colour[1];
   colour[2] = graphics_info_t::symmetry_colour[2];
   return colour;
}


// In future the gui will usefully set a mol number and we
// will use that.
//
void set_show_symmetry_master(short int state) {

   // std::cout << "set_show_symmetry_master() " << state << std::endl;

   //
   graphics_info_t g;

   // show symmetry state is no longer part of the molecule(s).


//       g.molecules[ii].show_symmetry = state;

//       if ( state == 1 ) {
// 	 g.molecules[mol_no].update_symmetry();
//       }
//    }

   g.show_symmetry = state;
   for (int ii=0; ii<g.n_molecules(); ii++)
      if (is_valid_model_molecule(ii))
	 graphics_info_t::molecules[ii].update_symmetry();
   graphics_draw();

   if (state) {
      // Now count the number of model molecules that have symmetry
      // available.  If there are none, then pop up a warning.

      int n_has_symm = 0;
      int n_model_molecules = 0;
      for (int ii=0; ii<g.n_molecules(); ii++) {
	 if (is_valid_model_molecule(ii)) {
	    n_model_molecules++;
	    mmdb::mat44 my_matt;
	    int err = graphics_info_t::molecules[ii].atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
	    if (err == mmdb::SYMOP_Ok) {
	       n_has_symm++;
	       break;
	    }
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


void
set_clipping_front(float v) {
   graphics_info_t g;
   g.set_clipping_front(v);
   std::string cmd = "set-clipping-front";
   std::vector<coot::command_arg_t> args;
   args.push_back(v);
   // std::cout << "mid-1 adding to historyin in set_clipping_front" << std::endl;
   add_to_history_typed(cmd, args);
   // std::cout << "done set_clipping_front" << std::endl;
}

void
set_clipping_back(float v) {
   graphics_info_t g;
   g.set_clipping_back(v);

   std::string cmd = "set-clipping-back";
   std::vector<coot::command_arg_t> args;
   args.push_back(v);
   add_to_history_typed(cmd, args);
}

/*! \brief get clipping plane front */
float get_clipping_plane_front() {
   graphics_info_t g;
   return g.get_clipping_plane_front();
}

/*! \brief get clipping plane back */
float get_clipping_plane_back() {
   graphics_info_t g;
   return g.get_clipping_plane_back();
}

/*! increase the *amount* of clipping, that is (independent of projection matrix)*/
void increase_clipping_front() {
   graphics_info_t g;
   g.increase_clipping_front();
}

/*! increase the *amount* of clipping, that is (independent of projection matrix)*/
void increase_clipping_back() {
   graphics_info_t g;
   g.increase_clipping_back();
}

/*! decrease the *amount* of clipping, that is (independent of projection matrix)*/
void decrease_clipping_front() {
   graphics_info_t g;
   g.decrease_clipping_front();
}

/*! decrease the *amount* of clipping, that is (independent of projection matrix)*/
void decrease_clipping_back() {
   graphics_info_t g;
   g.decrease_clipping_back();
}


/*  ----------------------------------------------------------------------- */
/*                         Colour                                           */
/*  ----------------------------------------------------------------------- */

void
set_symmetry_colour_merge(float v) {

   graphics_info_t::symmetry_colour_merge_weight = v;
   graphics_draw();

   std::string cmd = "set-symmetry-colour-merge";
   std::vector<coot::command_arg_t> args;
   args.push_back(v);
   add_to_history_typed(cmd, args);
}

/*! \brief set the symmetry colour base */
void set_symmetry_colour(float r, float g, float b) {

   graphics_info_t::symmetry_colour[0] = r;
   graphics_info_t::symmetry_colour[1] = g;
   graphics_info_t::symmetry_colour[2] = b;

   std::string cmd = "set-symmetry-colour";
   std::vector<coot::command_arg_t> args;
   args.push_back(r);
   args.push_back(g);
   args.push_back(b);
   add_to_history_typed(cmd, args);}


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
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
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

int get_colour_map_rotation_on_read_pdb_c_only_flag() {

  int ret = graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag;
  return ret;
}

/* widget work */
GtkWidget *wrapped_create_coords_colour_control_dialog() {

   // GtkWidget *w = create_coords_colour_control_dialog();
   GtkWidget *w = widget_from_builder("coords_colour_control_dialog");

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
   if (! g.glareas.empty())
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
      std::set<int> s; // dummy
      short int f = graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag;
      bool g = false; // goodsell_mode
      bool force_rebonding = false;
      graphics_info_t::molecules[imol].make_colour_by_chain_bonds(s,f,g, force_rebonding);
      graphics_draw();
   }
   std::string cmd = "set-colour-by-chain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

/*! \brief colour molecule number imol by chain type */
void set_colour_by_ncs_chain(int imol, short int goodsell_mode) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].make_colour_by_ncs_related_chains(goodsell_mode);
      graphics_draw();
   }
   std::string cmd = "set-colour-by-ncs-chain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}


void set_colour_by_chain_goodsell_mode(int imol) {

   if (is_valid_model_molecule(imol)) {
      std::set<int> s; // dummy
      short int f = graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag;
      bool g = true; // goodsell_mode
      bool force_rebonding = false;
      graphics_info_t::molecules[imol].make_colour_by_chain_bonds(s,f,g, force_rebonding);
      graphics_draw();
   }
   std::string cmd = "set-colour-by-chain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

//! \brief set the goodsell chain colour colour wheel step (default 0.22)
void set_goodsell_chain_colour_wheel_step(float s) {

   graphics_info_t::goodsell_chain_colour_wheel_rotation_step = s;
   // now generate new bonds for all molecules drawn in goodsell mode.
   for (int i=0; i<graphics_n_molecules(); i++) {
      if (is_valid_model_molecule(i)) {
         // molecules seem not to know that if they are drawn in goodsell mode. Hmm.
      }
   }
   graphics_draw();

}


void set_colour_by_molecule(int imol) {

   if (is_valid_model_molecule(imol)) {
      bool force_rebonding = false;
      graphics_info_t::molecules[imol].make_colour_by_molecule_bonds(force_rebonding);
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
   std::string cmd = "set-colour-map-rotation-for-map";
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

   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_unit_cell(state);
   }

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
	 graphics_info_t::molecules[imol].set_show_unit_cell(istate);
      }
      if (is_valid_map_molecule(imol)) {
	 graphics_info_t::molecules[imol].set_show_unit_cell(istate);
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
      std::cout << "Cannot interpret " << text << ".  Assuming 10A" << std::endl;
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
   logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__), "don't use this");
}

/*! \brief set show aniso atoms */
void set_show_aniso_atoms(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      bool st = state;
      graphics_info_t::molecules[imol].set_show_atoms_as_aniso(st);
   }
   graphics_draw();
}

/*! \brief set show aniso atoms as ortep */
void set_show_aniso_atoms_as_ortep(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_aniso_atoms_as_ortep(state);
   }
   graphics_draw();
}

/*! \brief set show aniso atoms as ortep */
void set_show_aniso_atoms_as_empty(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_aniso_atoms_as_empty(state);
   }
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

/*! \brief set default for the drawing of atoms in stick mode (default is on (1)) */
void set_draw_stick_mode_atoms_default(short int state) {
   graphics_info_t::draw_stick_mode_atoms_default = state;
}

void set_default_bond_thickness(int t) {

   graphics_info_t g;
   g.default_bond_width = t;

}

/*! \brief set the default representation type (default 1).*/
void set_default_representation_type(int type) {
   graphics_info_t g;
   g.default_bonds_box_type = type;
}


void set_bond_thickness(int imol, float t) {
   graphics_info_t g;
   // std::cout << "debug:: -----------------------------------set_bond_thickness() called with imol "
   // << imol << " thickness " << t << std::endl;
   g.set_bond_thickness(imol, t);
}

/*! \brief allow lines that are further away to be thinner */
void set_use_variable_bond_thickness(short int state) {
   graphics_info_t g;
   g.use_variable_bond_width = state;
}


void set_bond_thickness_intermediate_atoms(float t) {
   graphics_info_t g;
   g.set_bond_thickness_intermediate_atoms(t);
}

void set_bond_colour_rotation_for_molecule(int imol, float value) {
   if (is_valid_model_molecule(imol)) {
      // graphics_info_t::molecules[imol].bonds_colour_map_rotation = value;
      graphics_info_t::molecules[imol].update_bonds_colour_using_map_rotation(value);
      graphics_draw();
   }
}

float get_bond_colour_rotation_for_molecule(int imol) {

   float v = -1;
   if (is_valid_model_molecule(imol)) {
      v = graphics_info_t::molecules[imol].bonds_colour_map_rotation;
   }
   return v;
}


void set_unbonded_atom_star_size(float f) {
   graphics_info_t g;
   g.unbonded_atom_star_size = f;
}

int get_default_bond_thickness() {
   graphics_info_t g;
   int ret = g.default_bond_width;
   return ret;
}

int make_ball_and_stick(int imol,
                        const char *atom_selection_str,
                        float bond_thickness,
                        float sphere_size,
                        int do_spheres_flag) {

   int i = imol;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      GtkWidget *glarea_0 = 0;
      GtkWidget *glarea_1 = 0;
      if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
      if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
      gl_context_info_t glci(glarea_0, glarea_1);
      int dloi =
      graphics_info_t::molecules[imol].make_ball_and_stick(std::string(atom_selection_str),
                                                           bond_thickness,
                                                           sphere_size, do_spheres_flag,
                                                           glci, g.Geom_p());
      // std::cout << "dloi: " << dloi << std::endl;
      graphics_draw();
   }
   return i;
}


int
clear_ball_and_stick(int imol) {

  if (graphics_info_t::use_graphics_interface_flag) {
    if (is_valid_model_molecule(imol)) {
      GLuint dummy_tag = 0;
      graphics_info_t::molecules[imol].clear_display_list_object(dummy_tag);
      graphics_draw();
    }
  }
    return 0;
}

/*! \brief set the model molecule representation stye 0 for ball-and-stick/licorice (default) and 1 for ball */
void set_model_molecule_representation_style(int imol, unsigned int mode) {

   // modes defined in Mesh.hh as an enum

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_model_molecule_representation_style(mode);
   }
   graphics_draw();
}


/*! set show a ribbon/mesh for a given molecule */
void set_show_molecular_representation(int imol, int mesh_index, short int state) {

   if (is_valid_model_molecule(imol)) {
      if (mesh_index >= 0)
         if (mesh_index < static_cast<int>(graphics_info_t::molecules[imol].meshes.size())) {
            // graphics_info_t::molecules[imol].meshes[mesh_index].set_draw_mesh_state(state);
            graphics_info_t g;
            g.set_show_molecular_representation(imol, mesh_index, state);
         }
      graphics_draw();
   }
}


/* clear the given additional representation  */
void
set_show_additional_representation(int imol, int representation_number, int on_off_flag) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_additional_representation(representation_number,
									  on_off_flag);
   }
   graphics_draw();
}

/* \brief display/undisplay all the additional representations for the given molecule  */
void
set_show_all_additional_representations(int imol, int on_off_flag) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_all_additional_representations(on_off_flag);
   }
   graphics_draw();
}

/* \brief undisplay all the additional representations for the given
   molecule, except the given representation number (if it is off, leave it off)  */
void all_additional_representations_off_except(int imol, int representation_number,
					       short int ball_and_sticks_off_too_flag) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].all_additional_representations_off_except(representation_number, ball_and_sticks_off_too_flag);
   }
   graphics_draw();
}




/* delete the given additional representation  */
void delete_additional_representation(int imol, int representation_number) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].clear_additional_representation(representation_number);
   }
   graphics_draw();
}

/* return the index of the additional representation.  Return -1 on error */
int additional_representation_by_string(int imol,  const char *atom_selection_str,
					int representation_type,
					int bonds_box_type,
					float bond_width,
					int draw_hydrogens_flag) {
   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_selection_info_t info(atom_selection_str);
      graphics_info_t g;
      GtkWidget *dcw = g.display_control_window();
      GtkWidget *glarea_0 = 0;
      GtkWidget *glarea_1 = 0;
      if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
      if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
      gl_context_info_t glci(glarea_0, glarea_1);
      r = graphics_info_t::molecules[imol].add_additional_representation(representation_type,
									 bonds_box_type,
									 bond_width,
									 draw_hydrogens_flag,
									 info, dcw, glci,
									 g.Geom_p());
   }
   graphics_draw();
   return r;
}

/* return the index of the additional representation.  Return -1 on error */
int additional_representation_by_attributes(int imol,  const char *chain_id,
					    int resno_start, int resno_end,
					    const char *ins_code,
					    int representation_type,
					    int bonds_box_type,
					    float bond_width,
					    int draw_hydrogens_flag) {

   int r = -1;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      GtkWidget *dcw = g.display_control_window();
      coot::atom_selection_info_t info(chain_id, resno_start, resno_end, ins_code);
      GtkWidget *glarea_0 = 0;
      GtkWidget *glarea_1 = 0;
      if (g.glareas.size() > 0) glarea_0 = g.glareas[0];
      if (g.glareas.size() > 1) glarea_1 = g.glareas[1];
      gl_context_info_t glci(glarea_0, glarea_1);

      r = graphics_info_t::molecules[imol].add_additional_representation(representation_type,
									 bonds_box_type,
									 bond_width,
									 draw_hydrogens_flag,
									 info, dcw, glci, g.Geom_p());
   }
   graphics_draw();
   return r;
}


#ifdef USE_GUILE
SCM additional_representation_info_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      r = SCM_EOL;
      for (unsigned int ir=0; ir<graphics_info_t::molecules[imol].add_reps.size(); ir++) {
	 SCM l = SCM_EOL;
	 std::string s = graphics_info_t::molecules[imol].add_reps[ir].info_string();
	 SCM is_show_flag_scm = SCM_BOOL_F;
	 if (graphics_info_t::molecules[imol].add_reps[ir].show_it)
	    is_show_flag_scm = SCM_BOOL_T;
	 float bw = graphics_info_t::molecules[imol].add_reps[ir].bond_width;
	 SCM bond_width_scm = scm_from_double(bw);
	 SCM atom_spec_scm = SCM_BOOL_F;
	 // lines too long, make a rep
	 coot::additional_representations_t rep =
	    graphics_info_t::molecules[imol].add_reps[ir];
	 int type = rep.atom_sel_info.type;
	 if (type == coot::atom_selection_info_t::BY_STRING)
	    atom_spec_scm
	       = scm_from_locale_string(rep.atom_sel_info.atom_selection_str.c_str());
	 else
	    if (type == coot::atom_selection_info_t::BY_ATTRIBUTES) {
	       atom_spec_scm = SCM_EOL;
	       SCM ins_code_scm = scm_from_locale_string(rep.atom_sel_info.ins_code.c_str());
	       SCM resno_end_scm = scm_from_int(rep.atom_sel_info.resno_end);
	       SCM resno_start_scm   = scm_from_int(rep.atom_sel_info.resno_start);
	       SCM chain_id_scm    = scm_from_locale_string(rep.atom_sel_info.chain_id.c_str());
	       atom_spec_scm = scm_cons(ins_code_scm, atom_spec_scm);
	       atom_spec_scm = scm_cons(resno_end_scm, atom_spec_scm);
	       atom_spec_scm = scm_cons(resno_start_scm, atom_spec_scm);
	       atom_spec_scm = scm_cons(chain_id_scm, atom_spec_scm);
	    }

	 l = scm_cons(bond_width_scm, l);
	 l = scm_cons(is_show_flag_scm, l);
	 l = scm_cons(scm_from_locale_string(s.c_str()), l);
	 l = scm_cons(scm_from_int(ir), l);
	 r = scm_cons(l, r);
      }
   }
   return r;
}
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *additional_representation_info_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      r = PyList_New(0);
      for (unsigned int ir=0; ir<graphics_info_t::molecules[imol].add_reps.size(); ir++) {
	 PyObject *l = PyList_New(4);
	 std::string s = graphics_info_t::molecules[imol].add_reps[ir].info_string();
	 PyObject *is_show_flag_py = Py_False;
	 if (graphics_info_t::molecules[imol].add_reps[ir].show_it)
	    is_show_flag_py = Py_True;
	 float bw = graphics_info_t::molecules[imol].add_reps[ir].bond_width;
	 PyObject *bond_width_py = PyFloat_FromDouble(bw);
	 PyObject *atom_spec_py = Py_False;
	 // lines too long, make a rep
	 coot::additional_representations_t rep =
	    graphics_info_t::molecules[imol].add_reps[ir];
	 int type = rep.atom_sel_info.type;
	 if (type == coot::atom_selection_info_t::BY_STRING)
	    atom_spec_py
	       = myPyString_FromString(rep.atom_sel_info.atom_selection_str.c_str());
	 else
	    if (type == coot::atom_selection_info_t::BY_ATTRIBUTES) {
	       atom_spec_py = PyList_New(4);
	       PyObject *chain_id_py    = myPyString_FromString(rep.atom_sel_info.chain_id.c_str());
	       PyObject *resno_start_py = PyLong_FromLong(rep.atom_sel_info.resno_start);
	       PyObject *resno_end_py   = PyLong_FromLong(rep.atom_sel_info.resno_end);
	       PyObject *ins_code_py    = myPyString_FromString(rep.atom_sel_info.ins_code.c_str());
	       PyList_SetItem(atom_spec_py, 0, chain_id_py);
	       PyList_SetItem(atom_spec_py, 1, resno_start_py);
	       PyList_SetItem(atom_spec_py, 2, resno_end_py);
	       PyList_SetItem(atom_spec_py, 3, ins_code_py);
	    }
	 // we dont use the atom_spec_py!? -> decref
	 Py_XDECREF(atom_spec_py);

     Py_XINCREF(is_show_flag_py);
	 PyList_SetItem(l, 0, PyLong_FromLong(ir));
	 PyList_SetItem(l, 1, myPyString_FromString(s.c_str()));
	 PyList_SetItem(l, 2, is_show_flag_py);
	 PyList_SetItem(l, 3, bond_width_py);
	 PyList_Append(r, l);
	 Py_XDECREF(l);
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}

#endif	/* USE_PYTHON */


/*  ----------------------------------------------------------------------- */
/*                  dots display                                            */
/*  ----------------------------------------------------------------------- */
int dots(int imol,
	 const char *atom_selection_str,
	 const char *dots_name,
	 float dot_density, float sphere_size_scale) {

   int idots = -1;
   if (is_valid_model_molecule(imol)) {
      if (atom_selection_str) {
         // the colour is handled internally to make_dots - there the
         // state of molecule dots colour (see set_dots_colour()
         // below) is checked.
         idots = graphics_info_t::molecules[imol].make_dots(std::string(atom_selection_str),
							    dots_name,
							    dot_density,
							    sphere_size_scale);
      }
   }
   graphics_draw();
   return idots;
}

void set_dots_colour(int imol, float r, float g, float b) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_dots_colour(r,g,b);
      graphics_draw();
   }
}

void unset_dots_colour(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].unset_dots_colour();
      graphics_draw();
   }
}


void clear_dots(int imol, int dots_handle) {

   if (is_valid_model_molecule(imol)) {
      bool cleared_p = graphics_info_t::molecules[imol].clear_dots(dots_handle);
      if (cleared_p)
	 graphics_draw();
   }
}

/*! \brief clear the first dots object for imol with given name */
void clear_dots_by_name(int imol, const char *dots_object_name) {

   if (is_valid_model_molecule(imol)) {
      bool cleared = graphics_info_t::molecules[imol].clear_dots(dots_object_name);
      if (cleared)
	 graphics_draw();

   }
}


/* return the number of dots sets for molecule number imol */
int n_dots_sets(int imol) {

   int r = -1;

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules())) {
      r = graphics_info_t::molecules[imol].n_dots_sets();
   } else {
      // std::cout << "WARNING:: Bad molecule number: " << imol << std::endl;
      logger.log(log_t::WARNING, "Bad molecule number:", imol);
   }
   return r;
}

std::pair<short int, float> float_from_entry(GtkWidget *entry) {

   std::pair<short int, float> p(0,0);
   const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(entry));
   if (txt) {
      float f = atof(txt);
      p.second = f;
      p.first = 1;
   }
   return p;
}

std::pair<short int, int> int_from_entry(GtkWidget *entry) {

   std::pair<short int, int> p(0,0);
   const gchar *txt = gtk_editable_get_text(GTK_EDITABLE(entry));
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
      std::cout << "Cannot interpret " << text << ".  Assuming 10 steps" << std::endl;
      set_smooth_scroll_steps(10);
   }
}

// useful interface for scripting
void set_smooth_scroll_steps(int v) {
      graphics_info_t::smooth_scroll_n_steps = v;
}


char *get_text_for_smooth_scroll_steps() {

   char *text;

   text = (char *) malloc(100);
   snprintf(text, 99, "%-5d", graphics_info_t::smooth_scroll_n_steps);

   return text;
}

// useful interface for gui (entry)
void  set_smooth_scroll_limit_str(const char *text) {

   float v;

   v = atof(text);

   if (v >0 && v < 1000) {
      graphics_info_t::smooth_scroll_limit = v;
   } else {
      std::cout << text << " out of range: using 10A" << std::endl;
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

   if (0)
      std::cout << "--------------- set_font_colour called "
		<< red << " " << green << " " << blue << std::endl;

   graphics_info_t::font_colour = coot::colour_holder(red, green, blue);
   graphics_draw();
}

void set_use_stroke_characters(int state) {

   graphics_info_t::stroke_characters = state;
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
      std::cout << "Invalid cube size: " << text << ". Assuming 1.0A" << std::endl;
      val = 1.0;
   }
   g.user_defined_rotation_centre_crosshairs_size_scale_factor = val;
   graphics_draw();
}

void set_rotation_centre_size(float f) {
   graphics_info_t g;
   // g.rotation_centre_cube_size = f;
   g.user_defined_rotation_centre_crosshairs_size_scale_factor = f;
   graphics_draw();
}

void set_user_defined_rotation_centre_crosshairs_size_scale_factor(float f) {
   graphics_info_t g;
   g.user_defined_rotation_centre_crosshairs_size_scale_factor = f;
   graphics_draw();
}

/*! \brief set rotation centre colour

This is the colour for a dark background - if the background colour is not dark,
then the cross-hair colour becomes the inverse colour */
void set_rotation_centre_cross_hairs_colour(float r, float g, float b, float alpha) {

   glm::vec4 c(r,g,b,alpha);
   graphics_info_t gg;
   gg.set_rotation_centre_cross_hairs_colour(c);
   graphics_draw();
}


gchar *get_text_for_rotation_centre_cube_size() {

   char *text;
   graphics_info_t g;

   text = (char *)  malloc (100);
   snprintf(text, 90, "%-6.3f", g.user_defined_rotation_centre_crosshairs_size_scale_factor);
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


GtkWidget *main_window() {
   return graphics_info_t::get_main_window();
}

int get_number_of_molecules() {
   return graphics_info_t::n_molecules();
}


int graphics_n_molecules() {
   return graphics_info_t::n_molecules();
}

/* return either 1 (yes, there is at least one hydrogen) or 0 (no
   hydrogens, or no such molecule) */
int molecule_has_hydrogens_raw(int imol) {

   int r = 0;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].molecule_has_hydrogens();
   }
   return r;
}



/* a testing/debugging function.  Used in a test to make sure that the
   outside number of a molecule (the vector index) is the same as that
   embedded in the molecule description object.  Return -1 on
   non-valid passed imol. */
int own_molecule_number(int imol) {

  int r = -1;
  if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
     r = graphics_info_t::molecules[imol].MoleculeNumber();
  }
  return r;
}


/*  ----------------------------------------------------------------------- */
/*                  utility function                                        */
/*  ----------------------------------------------------------------------- */
// return -1 if atom not found.
int atom_index(int imol, const char *chain_id, int iresno, const char *atom_id) {

   int index;
   index = atom_index_full(imol, chain_id, iresno, "", atom_id, "");

   return index;
}
// return -1 if atom not found.
int atom_index_full(int imol, const char *chain_id, int iresno, const char *inscode, const char *atom_id, const char *altconf) {

   int index = -1;
   graphics_info_t g;

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules()) {
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

int atom_index_first_atom_in_residue_with_altconf(int imol,
						  const char *chain_id,
						  int iresno,
						  const char *ins_code,
						  const char *alt_conf) {

   int index = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      index = g.molecules[imol].atom_index_first_atom_in_residue(std::string(chain_id),
								 iresno,
								 std::string(ins_code),
								 std::string(alt_conf));
   }
   return index;
}

/*! \brief return the minimum residue number for imol chain chain_id */
int min_resno_in_chain(int imol, const char *chain_id) {

   int res_no_min = 999997;
   if (is_valid_model_molecule(imol)) {
      auto p = graphics_info_t::molecules[imol].min_res_no_in_chain(chain_id);
      if (p.first) {
         res_no_min = p.second;
      }
   }
   return res_no_min;

}

/*! \brief return the maximum residue number for imol chain chain_id */
int max_resno_in_chain(int imol, const char *chain_id) {

   int res_no_max = -99999;
   if (is_valid_model_molecule(imol)) {
      auto p = graphics_info_t::molecules[imol].max_res_no_in_chain(chain_id);
      if (p.first) {
         res_no_max = p.second;
      }
   }
   return res_no_max;
}



float median_temperature_factor(int imol) {

   float low_cut = 2.0;
   float high_cut = 100.0;
   bool low_cut_flag = false;
   bool high_cut_flag = false;

   float median = -1.0;
   if (is_valid_model_molecule(imol)) {
      median = coot::util::median_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
                                                     graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
                                                     low_cut, high_cut,
                                                     low_cut_flag,
                                                     high_cut_flag);
   } else {
      logger.log(log_t::WARNING, "No such molecule number:", imol);
   }
   return median;
}

float average_temperature_factor(int imol) {

   float low_cut = 2.0;
   float high_cut = 100.0;
   short int low_cut_flag = 0;
   short int high_cut_flag = 0;

   float av = -1.0;
   if (imol < graphics_info_t::n_molecules()) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 av = coot::util::average_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
						     graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
						     low_cut, high_cut,
						     low_cut_flag,
						     high_cut_flag);
      } else {
	 // std::cout << "WARNING:: molecule " << imol << " has no model\n";
         logger.log(log_t::WARNING, "Molecule:", imol, "has no model");
      }
   } else {
      logger.log(log_t::WARNING, "No such molecule as:", imol);
   }
   return av;
}

float standard_deviation_temperature_factor(int imol) {

   float low_cut = 2.0;
   float high_cut = 100.0;
   short int low_cut_flag = 0;
   short int high_cut_flag = 0;

   float av = -1.0;
   if (is_valid_model_molecule(imol)) {
      av = coot::util::standard_deviation_temperature_factor(graphics_info_t::molecules[imol].atom_sel.atom_selection,
							     graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
							     low_cut, high_cut,
							     low_cut_flag,
							     high_cut_flag);
   } else {
      // std::cout << "WARNING:: molecule " << imol << " is not a valid model\n";
      logger.log(log_t::WARNING, "Molecule:", imol, "is not a valid model");
   }
   return av;
}

char *centre_of_mass_string(int imol) {

#ifdef USE_GUILE
   char *s = 0; // guile/SWIG sees this as #f
   if (is_valid_model_molecule(imol)) {
      mmdb::realtype x, y, z;
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
#else
// BL says:: we need a python version
#ifdef USE_PYTHON
   char *s = 0; // should do for python too
   if (is_valid_model_molecule(imol)) {
      mmdb::realtype x, y, z;
      GetMassCenter(graphics_info_t::molecules[imol].atom_sel.atom_selection,
		    graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
		    x, y, z);
      std::string sc = "[";
      sc += graphics_info_t::float_to_string(x);
      sc += ",";
      sc += graphics_info_t::float_to_string(y);
      sc += ",";
      sc += graphics_info_t::float_to_string(z);
      sc += "]";
      s = new char [sc.length() + 1];
      strcpy(s, sc.c_str());
      return s;
   }
   return s;
#else
   return 0;
#endif // PYTHON
#endif // GUILE
}

#ifdef USE_PYTHON
char *centre_of_mass_string_py(int imol) {

   char *s = 0; // should do for python too
   if (is_valid_model_molecule(imol)) {
      mmdb::realtype x, y, z;
      GetMassCenter(graphics_info_t::molecules[imol].atom_sel.atom_selection,
		    graphics_info_t::molecules[imol].atom_sel.n_selected_atoms,
		    x, y, z);
      std::string sc = "[";
      sc += graphics_info_t::float_to_string(x);
      sc += ",";
      sc += graphics_info_t::float_to_string(y);
      sc += ",";
      sc += graphics_info_t::float_to_string(z);
      sc += "]";
      s = new char [sc.length() + 1];
      strcpy(s, sc.c_str());
      return s;
   }
   return s;
}
#endif // PYTHON


void clear_pending_picks() {
   graphics_info_t g;
   g.clear_pending_picks();
}

/* produce debugging output from problematic atom picking  */
void set_debug_atom_picking(int istate) {
   graphics_info_t::debug_atom_picking = istate;
}


void print_view_matrix() { 		/* print the view matrix */

   graphics_info_t g;
   GL_matrix m;
   // m.from_quaternion(g.quat);
   std::cout << "FIXME:: use glm::quat " << std::endl;
   std::cout << "View Matrix:" << std::endl;
   m.print_matrix();
}

float get_view_matrix_element(int row, int col) {

   // delete this function

   graphics_info_t g;
   GL_matrix m;
   std::cout << "FIXME:: use glm::quat " << std::endl;
   //  m.from_quaternion(g.quat);
   return m.matrix_element(row, col);
}


float get_view_quaternion_internal(int element) {

   // delete this function
#if 0
   if ((element >= 0) &&
       (element < 4)) {
      return graphics_info_t::quat[element];
   } else {
      std::cout << "Bad element for quaternion: " << element
		<< " returning dummy -9999" << std::endl;
      return -9999;
   }
#endif
   return 0;
}

void set_view_quaternion(float i, float j, float k, float l) {

   double mag2 = i*i + j*j + k*k + l*l;
   double mag=sqrt(mag2);

   if (fabs(mag) > 0.5) {
      graphics_info_t g;
      g.set_view_quaternion(l,i,j,k); // weird.
      graphics_draw();
   } else {
      std::cout << "Bad view quaternion" << std::endl;
   }
}



/* Return 1 if we moved to a molecule centre, else go to origin and
   return 0. */
/* centre on last-read (and displayed) molecule with zoom 100. */
// Also, return 0 if there are no molecules to centre on.
//
// However, if we are already *at* that molecule centre, Reset View
// moves to the centre of the next displayed molecule (with wrapping).
//
int reset_view() {

   int istat = 0;
   graphics_info_t g;
   coot::Cartesian new_centre(0,0,0);
   int new_centred_molecule;

   std::vector<std::pair<int, coot::Cartesian> > candidate_centres;
   for (int ii=0; ii<g.n_molecules(); ii++) {
      if (is_valid_model_molecule(ii)) {
	 if (mol_is_displayed(ii)) {
	    coot::Cartesian c = g.molecules[ii].centre_of_molecule();
	    std::pair<int, coot::Cartesian> p(ii, c);
	    candidate_centres.push_back(p);
	 }
      }
   }

   if (candidate_centres.size() > 0) {

      coot::Cartesian current_centre = g.RotationCentre();

      // Were we centred anywhere already?
      //
      bool was_centred = false;
      int centred_mol = -1;
      float best_fit_for_centred = 9001.1;

      for (unsigned int i=0; i<candidate_centres.size(); i++) {
	 float d = (candidate_centres[i].second - current_centre).length();
	 if (d < best_fit_for_centred) {
	    best_fit_for_centred = d;
	    if (d < 0.1) {
	       was_centred = 1;
	       centred_mol = candidate_centres[i].first;
	    }
	 }
      }

      if (! was_centred) {
	 // centre on the first molecule then
	 new_centre = candidate_centres[0].second;
	 new_centred_molecule = candidate_centres[0].first;
      } else {
	 // we we centred somewhere and we need to centre on the next
	 // molecule.

	 // What is the next molecule number? (Let's wrap to beginning
	 // if we are on the last molecule).  That only makes to do if
	 // there is more than 1 molecule to centre on.

	 new_centre = candidate_centres[0].second;
	 new_centred_molecule = candidate_centres[0].first;
	 if (candidate_centres.size() > 1) {
            // 1 -> 2 or 3 -> 0 say.
	                              // overwrite centre info if there is a
				      // candidate centre with a
				      // molecule number greater than
				      // the centred_mol;
	    for (unsigned int i=0; i<candidate_centres.size(); i++) {
	       if (candidate_centres[i].first > centred_mol) {
		  new_centred_molecule = candidate_centres[i].first;
		  new_centre = candidate_centres[i].second;
		  break;
	       }
	    }
	 }
      }


      // float size = g.molecules[new_centred_molecule].size_of_molecule();
      // if (size < 1.0) size = 1.0;
      // float new_zoom = 7.0*size;
      // g.setRotationCentreAndZoom(new_centre, new_zoom);
      g.setRotationCentre(new_centre); // it's a Cartesian

      for(int ii=0; ii<graphics_info_t::n_molecules(); ii++) {
	 graphics_info_t::molecules[ii].update_map(graphics_info_t::auto_recontour_map_flag);
	 graphics_info_t::molecules[ii].update_symmetry();
      }
      graphics_draw();

   }

   add_to_history_simple("reset-view");
   return istat;
}


/*! \brief set the view rotation scale factor

 Useful/necessary for high resolution displayed, where, without this factor
 the view doesn't rotate enough */
void set_view_rotation_scale_factor(float f) {

   graphics_info_t::view_rotation_per_pixel_scale_factor = f;

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

   for (int imol=0; imol<g.n_molecules(); imol++) {

      if (g.molecules[imol].has_xmap() &&
	  ! g.molecules[imol].xmap_is_diff_map) {
	 g.molecules[imol].greer_skeleton_draw_on = 1;
	 i_skel_set = 1;
      }
      if (i_skel_set) break;
   }
   graphics_draw();
}

void
skel_greer_off() {

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {

      if (graphics_info_t::molecules[imol].has_xmap() &&
	  ! graphics_info_t::molecules[imol].xmap_is_diff_map) {
	    graphics_info_t::molecules[imol].greer_skeleton_draw_on = 0;
      }
   }
}


void test_fragment() {

   graphics_info_t g;
   g.rotamer_graphs(0);
}

int write_connectivity(const char *monomer_name, const char *filename) {

   graphics_info_t g;
   return g.Geom_p()->hydrogens_connect_file(monomer_name, filename);
}


void screendump_image(const char *filename) {

   graphics_draw();  // does this improve the artifacts problem?
   graphics_draw();

   int istatus = graphics_info_t::screendump_image(filename);
   // std::cout << "INFO:: screendump_image status " << istatus << std::endl;
   logger.log(log_t::INFO, "screendump_image status", istatus);
   if (istatus == 1) {
      std::string s = "Screendump image ";
      s += filename;
      s += " written";
      graphics_info_t g;
      g.add_status_bar_text(s);
// BL says: we wanna be nice and convert ppm to bmp for windoze user!?
// but not if we have png!!!
// still use that function with png but only to open file
#ifdef WINDOWS_MINGW
#ifdef USE_PYTHON
      std::string cmd("ppm2bmp(");
      cmd += single_quote(coot::util::intelligent_debackslash(filename));
      cmd += ")";
      safe_python_command(cmd);
#endif // USE_PYTHON
#endif // MINGW
   }

   if (istatus == 0) {
      std::string s = "Failed to write screendump image ";
      s += filename;
      graphics_info_t g;
      g.add_status_bar_text(s);
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

#else
#ifdef USE_PYTHON
   std::string cmd("raytrace('raster3d',");
   cmd += single_quote(coot::util::intelligent_debackslash(r3d_name));
   cmd += ",";
   cmd += single_quote(coot::util::intelligent_debackslash(filename));
   cmd += ",1,1)";
   safe_python_command(cmd);
#endif // USE_PYTHON
#endif // USE_GUILE
}

#ifdef USE_PYTHON
void make_image_raster3d_py(const char *filename) {

   std::string r3d_name = filename;
   r3d_name += ".r3d";
   raster3d(r3d_name.c_str());
   std::string cmd("raytrace('raster3d',");
   cmd += single_quote(coot::util::intelligent_debackslash(r3d_name));
   cmd += ",";
   cmd += single_quote(coot::util::intelligent_debackslash(filename));
   cmd += ",1,1)";
   safe_python_command(cmd);
}
#endif // USE_PYTHON

void make_image_povray(const char *filename) {
   std::string pov_name = filename;
   pov_name += ".pov";
   povray(pov_name.c_str());
#ifdef USE_GUILE

   GtkAllocation allocation;
   if (! graphics_info_t::glareas.empty()) {
      gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
      int x_size = allocation.width;
      int y_size = allocation.height;
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
   }

#else
#ifdef USE_PYTHON
   GtkAllocation allocation;
   if (! graphics_info_t::glareas.empty()) {
      gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
      int x_size = allocation.width;
      int y_size = allocation.height;
      std::string cmd("raytrace('povray',");
      cmd += single_quote(coot::util::intelligent_debackslash(pov_name));
      cmd += ",";
      cmd += single_quote(coot::util::intelligent_debackslash(filename));
      cmd += ",";
      cmd += graphics_info_t::int_to_string(x_size);
      cmd += ",";
      cmd += graphics_info_t::int_to_string(y_size);
      cmd += ")";
      safe_python_command(cmd);
   }
#endif // USE_PYTHON
#endif // USE_GUILE
}

#ifdef USE_PYTHON
void make_image_povray_py(const char *filename) {
   std::string pov_name = filename;
   pov_name += ".pov";
   povray(pov_name.c_str());
   GtkAllocation allocation;
   gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
   int x_size = allocation.width;
   int y_size = allocation.height;
   std::string cmd("raytrace('povray',");
   cmd += single_quote(coot::util::intelligent_debackslash(pov_name));
   cmd += ",";
   cmd += single_quote(coot::util::intelligent_debackslash(filename));
   cmd += ",";
   cmd += graphics_info_t::int_to_string(x_size);
   cmd += ",";
   cmd += graphics_info_t::int_to_string(y_size);
   cmd += ")";
   safe_python_command(cmd);
}
#endif // USE_PYTHON


void
autobuild_ca_off() {

   graphics_info_t g;
   g.autobuild_flag = 0;

}



/*  ------------------------------------------------------------------------ */
/*                         clipping */
/*  ------------------------------------------------------------------------ */

void do_clipping1_activate() {

   std::cout << "############## do_clipping1_activate() " << std::endl;

   GtkScale *hscale;
   GtkAdjustment *adjustment;

   /* connect this to displaying the new clipping window */

   // GtkWidget *clipping_window = create_clipping_window();
   GtkWidget *clipping_window = widget_from_builder("clipping_window");

   GtkWidget *hscale1 = widget_from_builder("hscale1"); // really? needs testing

   hscale = GTK_SCALE(hscale1);
   /*    gtk_scale_set_draw_value(hscale, TRUE);  already does */

   adjustment = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -10.0, 20.0, 0.05, 4.0, 10.1));

   gtk_range_set_adjustment(GTK_RANGE(hscale), adjustment);
   g_signal_connect(G_OBJECT(adjustment), "value_changed",
		    G_CALLBACK(clipping_adjustment_changed), NULL);

   gtk_widget_set_visible(clipping_window, TRUE);

}

void clipping_adjustment_changed (GtkAdjustment *adj, GtkWidget *window) {

   /*    printf("Clipping adjustment: %f\n", adj->value); */

   set_clipping_front(gtk_adjustment_get_value(adj));
   set_clipping_back (gtk_adjustment_get_value(adj));
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
   if (is_valid_model_molecule(imol)) {
      ierr = graphics_info_t::molecules[imol].save_coordinates(filename);
   }

   std::vector<std::string> command_strings;
   command_strings.push_back("save-coordinates");
   command_strings.push_back(coot::util::int_to_string(imol));
   command_strings.push_back(single_quote(filename));
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
read_phs_and_coords_and_make_map(const char *pdb_filename){

   // This function is the .phs equivalent of c.f. make_and_draw_map,
   // map_fill_from_mtz.  We have previously stored the phs_filename
   // in the static graphics_info_t.
   //
   graphics_info_t g;

   int imol = graphics_info_t::create_molecule();

   // don't forget that this is a map.
   //
   int istat = g.molecules[imol].make_map_from_phs(std::string(pdb_filename),
						   g.get_phs_filename());

   if (istat != -1) {
      graphics_draw();
   } else {
      // give us a warning message then
      g.erase_last_molecule();
      std::string w = "Sadly, the cell or space group is not comprehensible in\n";
      w += "the pdb file: ";
      w += pdb_filename;
      w += "\n";
      w += "Can't make map from phs file.";
      graphics_info_t g;
      g.info_dialog(w);
   }
}

/*! \brief read a phs file, the cell and symm information is from
  previously read (most recently read) coordinates file

 For use with phs data filename provided on the command line */
int
read_phs_and_make_map_using_cell_symm_from_previous_mol(const char *phs_filename) {

   clipper::Spacegroup spacegroup;
   clipper::Cell cell;
   int r = -1;

   int imol_ref = -1;

   for (int i=graphics_info_t::n_molecules()-1; i>=0; i--) {
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
   int imol = g.create_molecule();

   clipper::Spacegroup spacegroup;
   clipper::Cell cell;
   short int got_cell_symm_flag = 0;
   int istat = -1; // returned value

   if (is_valid_model_molecule(imol_ref) || is_valid_map_molecule(imol_ref)) {
      if (g.molecules[imol_ref].have_unit_cell) {
	 try {
	    std::pair<clipper::Cell,clipper::Spacegroup> xtal =
	       coot::util::get_cell_symm( g.molecules[imol_ref].atom_sel.mol );
	    cell = xtal.first;
	    spacegroup = xtal.second;
	    got_cell_symm_flag = 1;
	 } catch (const std::runtime_error &except ) {
	    std::cout << "WARNING:: Cant get spacegroup from coordinates!\n";
	    // get the cell/symm from a map:
	    if (g.molecules[imol_ref].has_xmap()) {
	       cell = g.molecules[imol_ref].xmap.cell();
	       spacegroup = g.molecules[imol_ref].xmap.spacegroup();
	       got_cell_symm_flag = 1;
	    }
	 }
      }
   }

   if (got_cell_symm_flag) {

      // don't forget that this is a map.
      //
      std::string phs_file(phs_filename);
      istat = g.molecules[imol].make_map_from_phs_using_reso(phs_file,
							     spacegroup,
							     cell,
							     reso_lim_low, reso_lim_high,
							     graphics_info_t::map_sampling_rate);

      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 imol = istat;
	 graphics_draw();
      } else {
	 g.erase_last_molecule();
	 std::string w = "Sadly, something bad happened reading phs file using\n";
	 w += "the molecule number ";
	 w += coot::util::int_to_string(imol_ref);
	 w += "\n";
	 w += "Can't make map from phs file.";
         g.info_dialog(w);
      }
   } else {
      g.erase_last_molecule();
      // give us a warning message then
      std::string w = "Sadly, the cell or space group is not comprehensible in\n";
      w += "the molecule number ";
      w += coot::util::int_to_string(imol_ref);
      w += "\n";
      w += "Can't make map from phs file.";
      g.info_dialog(w);
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
	 try {
	    std::pair<clipper::Cell,clipper::Spacegroup> xtal =
	       coot::util::get_cell_symm( g.molecules[imol_ref].atom_sel.mol );
	    cell = xtal.first;
	    spacegroup = xtal.second;
	    got_cell_symm_flag = 1;
	 }
	 catch (const std::runtime_error &except) {
	    std::cout << "WARNING:: Cant get spacegroup from coordinates!\n";
	    // get the cell/symm from a map:
	    if (g.molecules[imol_ref].has_xmap()) {
	       cell = g.molecules[imol_ref].xmap.cell();
	       spacegroup = g.molecules[imol_ref].xmap.spacegroup();
	       got_cell_symm_flag = 1;
	    }
	 }
      }

      if (got_cell_symm_flag) {
	 std::string phs_filename(phs_filename_str);

	 imol = g.create_molecule();
	 g.molecules[imol].make_map_from_phs(spacegroup, cell, phs_filename);
	 g.scroll_wheel_map = imol;
	 graphics_draw();
      }
   }

   return imol;
}


int
read_phs_and_make_map_using_cell_symm_from_mol_using_implicit_phs_filename(int imol_ref) {

   clipper::Spacegroup spacegroup;
   clipper::Cell cell;
   short int got_cell_symm_flag = 0;
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
	 try {
	    std::pair<clipper::Cell,clipper::Spacegroup> xtal =
	       coot::util::get_cell_symm( g.molecules[imol_ref].atom_sel.mol );
	    cell = xtal.first;
	    spacegroup = xtal.second;
	    got_cell_symm_flag = 1;
	 } catch (const std::runtime_error &except) {
	    std::cout << "WARNING:: Cant get spacegroup from coordinates!\n";
	    // get the cell/symm from a map:
	    if (g.molecules[imol_ref].has_xmap()) {
	       cell = g.molecules[imol_ref].xmap.cell();
	       spacegroup = g.molecules[imol_ref].xmap.spacegroup();
	       got_cell_symm_flag = 1;
	    }
	 }
      }

      if (got_cell_symm_flag) {
	 std::string phs_filename(g.get_phs_filename());

	 imol = g.create_molecule();
	 g.molecules[imol].make_map_from_phs(spacegroup, cell, phs_filename);
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

   int imol = g.create_molecule();
   g.molecules[imol].make_map_from_phs(spacegroup, cell, phs_filename);
   graphics_draw();
   return imol;
}


void
graphics_store_phs_filename(const char *phs_filename) {

   graphics_info_t g;
   g.set_phs_filename(std::string(phs_filename));
}


short int possible_cell_symm_for_phs_file() {

   if (graphics_info_t::n_molecules() == 0) {
      return 0;
   } else {
      return 1;
   }
}

// return a string to each of the cell parameters in molecule imol.
//
gchar *get_text_for_phs_cell_chooser(int imol, const char *field) {

   // we first look in atomseletion

   graphics_info_t g;
   gchar *retval = NULL;
   retval = (gchar *) malloc(12);
   int ihave_cell = 0;
   mmdb::realtype cell[6];
   const char *spgrp = NULL;

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules()) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    if (g.molecules[imol].have_unit_cell) {

	       ihave_cell = 1;

	       mmdb::realtype vol;
	       int orthcode;
	       g.molecules[imol].atom_sel.mol->GetCell(cell[0], cell[1], cell[2],
						       cell[3], cell[4], cell[5],
						       vol, orthcode);
	       spgrp   = g.molecules[imol].atom_sel.mol->GetSpaceGroup();

	    } else {

	       ihave_cell = 1;

	       clipper::Spacegroup spacegroup = g.molecules[imol].xmap.spacegroup();
	       clipper::Cell       ccell      = g.molecules[imol].xmap.cell();

	       cell[0] = g.molecules[imol].xmap.cell().a();
	       cell[1] = g.molecules[imol].xmap.cell().b();
	       cell[2] = g.molecules[imol].xmap.cell().c();
	       cell[3] = g.molecules[imol].xmap.cell().alpha() * RADTODEG;
	       cell[4] = g.molecules[imol].xmap.cell().beta()  * RADTODEG;
	       cell[5] = g.molecules[imol].xmap.cell().gamma() * RADTODEG;

	       spgrp = spacegroup.descr().symbol_hm().c_str();
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

int go_to_atom_molecule_number() {
   graphics_info_t g;
   return g.go_to_atom_molecule();
}


// 20220723-PE make these strings
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------


// Note that t3 is an atom name with (possibly) an altLoc tag (after the comma).
//
int set_go_to_atom_chain_residue_atom_name(const char *t1, int iresno, const char *t3) {

   graphics_info_t g;
   int success = set_go_to_atom_chain_residue_atom_name_no_redraw(t1, iresno, t3, 1);
   if (success) {
      mmdb::Atom *at = 0; // passed but not used, it seems.
      GtkWidget *window = graphics_info_t::go_to_atom_window;
      if (window)
	 g.update_widget_go_to_atom_values(window, at);
   }
   g.update_environment_distances_by_rotation_centre_maybe(go_to_atom_molecule_number());
   graphics_draw();
   return success;
}

int set_go_to_atom_chain_residue_atom_name_full(const char *chain_id,
						int resno,
						const char *ins_code,
						const char *atom_name,
						const char *alt_conf) {

   graphics_info_t g;
   g.set_go_to_atom_chain_residue_atom_name(chain_id, resno, ins_code, atom_name, alt_conf);
   int success = g.try_centre_from_new_go_to_atom();
   if (success) {
      mmdb::Atom *at = 0; // passed but not used, it seems.
      GtkWidget *window = graphics_info_t::go_to_atom_window;
      if (window)
	 g.update_widget_go_to_atom_values(window, at);
   }
   graphics_draw();
   return success;
}




// Note that t3 is an atom name with (possibly) an altLoc tag (after the comma).
//
int set_go_to_atom_chain_residue_atom_name_no_redraw(const char *t1, int iresno, const char *t3,
						     short int make_the_move_flag) {


   graphics_info_t g;
   // so we need to split t3 if it has a comma
   //
   std::string t3s(t3);
   std::string::size_type icomma = t3s.find_last_of(",");
   if (icomma == std::string::npos) {

      // there was no comma, conventional usage:
      g.set_go_to_atom_chain_residue_atom_name(t1, iresno, t3);

   } else {

      std::string atname = t3s.substr(0,icomma);
      std::string altloc = t3s.substr(icomma+1, t3s.length());
      g.set_go_to_atom_chain_residue_atom_name(t1, iresno,
					       atname.c_str(),
					       altloc.c_str());

   }
   mmdb::Atom *at = 0; // passed but not used, it seems.
   GtkWidget *window = graphics_info_t::go_to_atom_window;
   if (window)
      g.update_widget_go_to_atom_values(window, at);

   int success = 0;
   if (make_the_move_flag) {
      success = g.try_centre_from_new_go_to_atom();
   } else {
      success = 1;
   }
   g.update_things_on_move();
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
goto_next_atom_maybe_new() {

//    int it2 = atoi(t2);
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, +1);

   graphics_info_t g;
   return g.intelligent_next_atom_centring();

}

int
goto_previous_atom_maybe_new() {

//    int it2 = atoi(t2);
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, +1);

   graphics_info_t g;
   return g.intelligent_previous_atom_centring();
}


// DELETE ME.
// int
// goto_prev_atom_maybe(const gchar *t1, const gchar *t2, const gchar *t3,
// 		     GtkEntry *res_entry) {

//    int it2 = atoi(t2);
//    return goto_near_atom_maybe(t1, it2, t3, res_entry, -1);

// }

//
// int
// goto_near_atom_maybe(const char *t1, int ires, const char *t3,
// 		     GtkEntry *res_entry, int idiff) {

//    graphics_info_t g;

//    int ires_l = ires + idiff ; // for next residue, or previous.

//    g.set_go_to_atom_chain_residue_atom_name(t1, ires_l, t3);

//    int success = g.try_centre_from_new_go_to_atom();

//    if (success) {
//       char *txt = (char *)malloc(6);
//       snprintf(txt, 5, "%d", ires_l);
//       gtk_entry_set_text(GTK_ENTRY(res_entry), txt);
//       update_things_on_move_and_redraw();
//    }
//    return success;
// }


/* For dynarama callback sake. The widget/class knows which coot
   molecule that it was generated from, so in order to go to the
   molecule from dynarama, we first need to the the molecule - because
   set_go_to_atom_chain_residue_atom_name() does not mention the
   molecule (see "Next/Previous Residue" for reasons for that).  This
   function simply calls the graphics_info_t function of the same
   name. */
void set_go_to_atom_molecule(int imol) {

   graphics_info_t g;
   int current_go_to_atom_molecule = g.go_to_atom_molecule();
   g.set_go_to_atom_molecule(imol);
   if (current_go_to_atom_molecule != imol)
      update_go_to_atom_window_on_other_molecule_chosen(imol);
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
   if (is_valid_model_molecule(imol)) {
      bool force_rebonding = false;
      std::cout << "calling ca_representation() for imol " << imol << std::endl;
      g.molecules[imol].ca_representation(force_rebonding);
   } else {
      std::cout << "WARNING:: no such valid molecule " << imol
                << " in graphics_to_ca_representation" << std::endl;
   }
   graphics_draw();

   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}

/*! \brief draw molecule number imol coloured by chain */
void graphics_to_colour_by_chain(int imol) {

   if (is_valid_model_molecule(imol)) {
      bool force_rebonding = false;
      graphics_info_t::molecules[imol].make_colour_by_chain_bonds(force_rebonding);
      graphics_draw();
   }
}


void graphics_to_ca_plus_ligands_representation(int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      bool force_rebonding = false;
      g.molecules[imol].ca_plus_ligands_representation(g.Geom_p(), force_rebonding);
      graphics_draw();
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-plus-ligands-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}

void graphics_to_ca_plus_ligands_and_sidechains_representation   (int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].ca_plus_ligands_and_sidechains_representation(g.Geom_p());
      graphics_draw();
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("graphics-to-ca-plus-ligands-and-sidechains-representation");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}


void graphics_to_bonds_representation(int imol) {
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      bool force_rebonding = false;
      g.molecules[imol].bond_representation(g.Geom_p(), force_rebonding);
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-bonds-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_bonds_representation" << std::endl;
   graphics_draw();

}

void graphics_to_colour_by_molecule(int imol) {

   if (is_valid_model_molecule(imol)) {
      bool force_rebonding = false;
      graphics_info_t::molecules[imol].make_colour_by_molecule_bonds(force_rebonding);
      graphics_draw();
   }
   std::string cmd = "set-colour-by-molecule";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

void graphics_to_bonds_no_waters_representation(int imol) {
   graphics_info_t g;
   if (is_valid_model_molecule(imol)){
      g.molecules[imol].bonds_no_waters_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-bonds-no-waters-representation");
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
      g.molecules[imol].ca_plus_ligands_sec_struct_representation(g.Geom_p());
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
      graphics_info_t g;
      g.molecules[imol].ca_plus_ligands_rainbow_representation(g.Geom_p());
      std::vector<std::string> command_strings;
      // BL says:: or maybe we want to keep the following name and change above
      //command_strings.push_back("graphics-to-ca-plus-ligands-rainbow-representation");
      command_strings.push_back("graphics-to-rainbow-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		//BL says:: as above
		//<< " in graphics_to_ca_plus_ligands_rainbow_representation"
		<< " in graphics_to_rainbow_representation"
		<< std::endl;
   graphics_draw();
}

void graphics_to_b_factor_representation(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].b_factor_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-b-factor-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   } else {
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_b_factor_representation"
		<< std::endl;
   }
   graphics_draw();
}

void graphics_to_b_factor_cas_representation(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].b_factor_representation_as_cas();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-b-factor-cas-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   } else {
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_b_factor_representation"
		<< std::endl;
   }
   graphics_draw();
}

void graphics_to_occupancy_representation(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].occupancy_representation();
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-occupancy-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   }
   else
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_occupancy_representation"
		<< std::endl;
   graphics_draw();
}

/*! \brief draw molecule number imol coloured by user-defined atom colours */
void graphics_to_user_defined_atom_colours_representation(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      bool all_atoms_flag = false;
      g.molecules[imol].user_defined_colours_representation(g.Geom_p(), all_atoms_flag, g.draw_missing_loops_flag);
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-user-defined-colours-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   } else {
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_occupancy_representation"
		<< std::endl;
   }
   graphics_draw();
}

/*! \brief draw molecule number imol all atoms coloured by user-defined atom colours */
void graphics_to_user_defined_atom_colours_all_atoms_representation(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      bool all_atoms_flag = true;
      g.molecules[imol].user_defined_colours_representation(g.Geom_p(), all_atoms_flag, g.draw_missing_loops_flag);
      std::vector<std::string> command_strings;
      command_strings.push_back("graphics-to-user-defined-colours-representation");
      command_strings.push_back(graphics_info_t::int_to_string(imol));
      add_to_history(command_strings);
   } else {
      std::cout << "WARNING:: no such valid molecule " << imol
		<< " in graphics_to_occupancy_representation"
		<< std::endl;
   }
   graphics_draw();
}



/*! \brief make the carbon atoms for molecule imol be grey
 */
void set_use_grey_carbons_for_molecule(int imol, short int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_use_bespoke_carbon_atom_colour(state);
      graphics_draw();
   }

}
/*! \brief set the colour for the carbon atoms

can be not grey if you desire, r, g, b in the range 0 to 1.
 */
void set_grey_carbon_colour(int imol, float r, float g, float b) {

   if (is_valid_model_molecule(imol)) {
      coot::colour_t col(r,g,b);
      graphics_info_t::molecules[imol].set_bespoke_carbon_atom_colour(col);
      graphics_draw();
   }
}

/* undocumented feature for development. */
void set_draw_moving_atoms_restraints(int state) {
   graphics_info_t::draw_it_for_moving_atoms_restraints_graphics_object_user_control = state;
   graphics_draw();
}

/* undocumented feature for development. */
short int get_draw_moving_atoms_restraints() {
   return graphics_info_t::draw_it_for_moving_atoms_restraints_graphics_object_user_control;
}



int
get_graphics_molecule_bond_type(int imol) {

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

void change_model_molecule_representation_mode(int up_or_down) {

   graphics_info_t g;
   g.change_model_molecule_representation_mode(up_or_down);
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

void graphics_to_phenix_geo_representation(int imol, int mode, const coot::phenix_geo::phenix_geometry &g) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].update_bonds_using_phenix_geo(g);
      graphics_draw();
   }
}


void graphics_to_phenix_geo_representation(int imol, int mode,
					   const std::string &geo_file_name) {

   coot::phenix_geo::phenix_geometry pg;
   pg.parse(geo_file_name);
   graphics_to_phenix_geo_representation(imol, mode, pg);

}

void phenix_geo_validation_buttons(int imol,
                                   const coot::phenix_geo::phenix_geometry &pg,
                                   double residual_cutoff);

//! \brief validate using phenix geo bonds
//!
//! Typically this would be called shortly after
//!
//! @param imol
void validate_using_phenix_geo_bonds(int imol, const std::string &geo_file_name) {

   if (is_valid_model_molecule(imol)) {
      coot::phenix_geo::phenix_geometry pg;
      pg.parse(geo_file_name);
      float residual_criterion = 4.4;
      phenix_geo_validation_buttons(imol, pg, residual_criterion);
   }
}

/*  Not today
void set_ca_bonds_loop_params(float p1, float p2, float p3) {

    graphics_info_t::ca_bonds_loop_param_1 = p1;
    graphics_info_t::ca_bonds_loop_param_2 = p2;
    graphics_info_t::ca_bonds_loop_param_3 = p3;

}

*/


// -------------------------------------------------------------------------
//                        skeletonization level
// -------------------------------------------------------------------------
//
gchar *get_text_for_skeletonization_level_entry() {

   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);
   g_snprintf(txt, 9, "%f", g.skeleton_level);

   return txt;
}

void set_skeletonization_level_from_widget(const char *txt) {

   float tmp;
   graphics_info_t g;

   tmp = atof(txt);

   if (tmp > 0.0 &&  tmp < 9999.9) {
      g.skeleton_level = tmp;
   } else {
      std::cout << "Cannot interpret " << txt << " using 0.2 instead" << std::endl;
      g.skeleton_level = 0.2;
   }

   for (int imol=0; imol<g.n_molecules(); imol++) {
      if (g.molecules[imol].has_xmap() &&
	  g.molecules[imol].xmap_is_diff_map != 1) {
	 //
	 g.molecules[imol].update_clipper_skeleton();
      }
   }
   graphics_draw();
}

gchar *get_text_for_skeleton_box_size_entry() {

   graphics_info_t g;
   gchar *txt = (gchar *)malloc(10);

   g_snprintf(txt, 9, "%f", g.skeleton_box_radius);
   return txt;
}

void set_skeleton_box_size_from_widget(const char *txt) {
   float tmp;
   graphics_info_t g;

   tmp = atof(txt);

   if (tmp > 0.0 &&  tmp < 9999.9) {
      g.skeleton_box_radius = tmp;
   } else {

      std::cout << "Cannot interpret " << txt << " using 0.2 instead" << std::endl;
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

   for (int imol=0; imol<g.n_molecules(); imol++) {
      if (g.molecules[imol].has_xmap() &&
	  g.molecules[imol].xmap_is_diff_map != 1) {
	 //
	 g.molecules[imol].update_clipper_skeleton();
      }
   }
   graphics_draw();
}

/*  ----------------------------------------------------------------------- */
/*                  Utility Functions                                       */
/*  ----------------------------------------------------------------------- */
// These functions are for storing the molecule number and (some other
// number) as an int and used with GPOINTER_TO_INT and GINT_TO_POINTER.
int encode_ints(int i1, int i2) {
   int i = 1000 * i1 + i2;
   return i;
}

std::pair<int, int> decode_ints(int i) {
   int j = i/1000;
   int k = i - j * 1000;
   return std::pair<int, int>(j,k);
}

void store_keyed_user_name(std::string key, std::string user_name, std::string passwd) {

   std::pair<std::string, std::string> p(user_name, passwd);
   graphics_info_t::user_name_passwd_map[key] = p;

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

   GtkWidget *widget = wrapped_create_display_control_window(); // uses gtkbuilder
   gtk_widget_set_visible(widget, TRUE);
   std::vector<std::string> command_strings;
   command_strings.push_back("post-display-control-window");
   add_to_history(command_strings);

}

// add the calling function here.
void
clear_out_container(GtkWidget *vbox) {

   graphics_info_t g;
   g.clear_out_container(vbox);

}


void add_map_display_control_widgets() {

   // 20220808-PE we probably don't need to do this - now that the display control dialog is not destroyed

   graphics_info_t g;

   GtkWidget *map_vbox = widget_from_builder("display_map_vbox");
   clear_out_container(map_vbox);

   for (int ii=0; ii<g.n_molecules(); ii++)
      if (g.molecules[ii].has_xmap() || g.molecules[ii].has_nxmap())
	 g.molecules[ii].update_map_in_display_control_widget();

}


void add_mol_display_control_widgets() {

   // 20220808-PE we probably don't need to do this - now that the display control dialog is not destroyed

   graphics_info_t g;

   GtkWidget *molecule_vbox = widget_from_builder("display_molecule_vbox");
   clear_out_container(molecule_vbox);

   for (int ii=0; ii<g.n_molecules(); ii++) {
      if (g.molecules[ii].has_model()) {
	 g.molecules[ii].new_coords_mol_in_display_control_widget();
      }
   }
}


void add_map_and_mol_display_control_widgets() {

   // 20220808-PE we probably don't need to do this - now that the display control dialog is not destroyed

   add_mol_display_control_widgets();
   add_map_display_control_widgets();
}


// resets to NULL the scroll group too.
void reset_graphics_display_control_window() {
   graphics_info_t g;
   g.save_display_control_widget_in_graphics(NULL);
}

void close_graphics_display_control_window() {
   graphics_info_t g;
   GtkWidget *w = g.display_control_window();
   if (w) {
      // gtk_widget_destroy(w); // nope
      gtk_widget_set_visible(w, FALSE);
      reset_graphics_display_control_window();
   }
}

/*! \brief make the map displayed/undisplayed, 0 for off, 1 for on */
void set_map_displayed(int imol, int state) {

   graphics_info_t g;
   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].set_map_is_displayed(state);
      set_display_control_button_state(imol, "Displayed", state);
      graphics_draw();
   }
}

void set_draw_map_standard_lines(int imol, short int state) {

   graphics_info_t g;
   if (is_valid_map_molecule(imol)) {
      g.molecules[imol].set_map_is_displayed_as_standard_lines(state);
      graphics_draw();
   }
}


#ifdef USE_GUILE
void display_maps_scm(SCM maps_list_scm) {

   graphics_info_t g;
   int n_mol = graphics_n_molecules();
   std::vector<bool> map_on(n_mol, false);
   SCM s_test = scm_list_p(maps_list_scm);
   if (scm_is_true(s_test)) {
      SCM n_scm = scm_length(maps_list_scm);
      int n = scm_to_int(n_scm);
      for (int i=0; i<n; i++) {
	 SCM item_scm = scm_list_ref(maps_list_scm, scm_from_int(i));
	 if (scm_is_true(scm_integer_p(item_scm))) {
	    int imol = scm_to_int(item_scm);
	    if (is_valid_map_molecule(imol)) {
	       map_on[imol] = true;
	    }
	 }
      }
   }

   for (int imol=0; imol<n_mol; imol++) {
      if (is_valid_map_molecule(imol)) {
	 if (map_on[imol])
	    g.molecules[imol].set_map_is_displayed(1);
	 else
	    g.molecules[imol].set_map_is_displayed(0);
      }
   }
   graphics_draw();
}
#endif


#ifdef USE_PYTHON
void display_maps_py(PyObject *pyo) {

   graphics_info_t g;
   int n_mol = graphics_n_molecules();
   std::vector<bool> map_on(n_mol, false);
   if (PyList_Check(pyo)) {
      int n = PyObject_Length(pyo);
      for (int i=0; i<n; i++) {
	 PyObject *item_py = PyList_GetItem(pyo, i);
	 if (PyLong_Check(item_py)) {
	    int imol = PyLong_AsLong(item_py);
	    if (is_valid_map_molecule(imol)) {
	       map_on[imol] = true;
	    }
	 }
      }
   }
   for (int imol=0; imol<n_mol; imol++) {
      if (is_valid_map_molecule(imol)) {
	 if (map_on[imol])
	    g.molecules[imol].set_map_is_displayed(1);
	 else
	    g.molecules[imol].set_map_is_displayed(0);
      }
   }
   graphics_draw();
}
#endif



// button_type is "Displayed" or "Active"
void
set_display_control_button_state(int imol, const std::string &button_type, int state) {

   // button type is "Active" or "Displayed"
   if (false)
      std::cout << "start: set_display_control_button_state() imol " << imol << " type " << button_type
                << " new_state: " << state << std::endl;

   if (! graphics_info_t::use_graphics_interface_flag) return;

   if (is_valid_model_molecule(imol)) {

      GtkWidget *display_control_vbox = nullptr;
      if (is_valid_model_molecule(imol))
         display_control_vbox = widget_from_builder("display_molecule_vbox");

      if (GTK_IS_BOX(display_control_vbox)) {

         GtkWidget *item_widget = gtk_widget_get_first_child(display_control_vbox);
         while (item_widget) {
            int imol_widget = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget), "imol"));
            if (imol_widget == imol) {
               // item_widget is a mol_vbox, the hbox with controls is its first child
               GtkWidget *hbox = gtk_widget_get_first_child(item_widget);
               if (hbox) {
                  GtkWidget *display_check_button = GTK_WIDGET(g_object_get_data(G_OBJECT(hbox), "display_check_button"));
                  GtkWidget *active_check_button  = GTK_WIDGET(g_object_get_data(G_OBJECT(hbox), "active_check_button"));
                  if (button_type == "Displayed" && display_check_button) {
                     gtk_check_button_set_active(GTK_CHECK_BUTTON(display_check_button), state);
                  }
                  if (button_type == "Active" && active_check_button) {
                     gtk_check_button_set_active(GTK_CHECK_BUTTON(active_check_button), state);
                  }
               }
            }
            item_widget = gtk_widget_get_next_sibling(item_widget);
         };
      }
   }

   if (is_valid_map_molecule(imol)) {
      GtkWidget *display_control_vbox = nullptr;
      if (is_valid_map_molecule(imol))
         display_control_vbox = widget_from_builder("display_map_vbox");
      if (GTK_IS_BOX(display_control_vbox)) {
         GtkWidget *item_widget = gtk_widget_get_first_child(display_control_vbox);
         while (item_widget) {
            GtkWidget *child_0 = gtk_widget_get_first_child(item_widget);
            GtkWidget *child_1 = gtk_widget_get_next_sibling(child_0);
            GtkWidget *child_2 = gtk_widget_get_next_sibling(child_1);
            GtkWidget *display_check_button = child_2;
            // std::cout << "child_2 " << child_2 << std::endl;
            int imol_widget = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget), "imol"));
            if (imol_widget == imol) {
               if (button_type == "Displayed") {
                  gtk_check_button_set_active(GTK_CHECK_BUTTON(display_check_button), state);
               }
            }
            item_widget = gtk_widget_get_next_sibling(item_widget);
         };
      }
   }
}

/*! \brief make the coordinates molecule displayed/undisplayed, 0 for off, 1 for on */
void set_mol_displayed(int imol, int state) {

   // std::cout << "debug:: set_mol_displayed called with imol " << imol << " state " << state << std::endl;

   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      int current_state = graphics_info_t::molecules[imol].get_mol_is_displayed();
      if (current_state != state) {
         graphics_info_t::molecules[imol].set_mol_is_displayed(state);
         if (g.use_graphics_interface_flag)
            set_display_control_button_state(imol, "Displayed", state);
         graphics_draw();
      }
   } else {
      std::cout << "not valid molecule" << std::endl;
   }
}


/*! \brief from all the model molecules, display only imol

This stops flashing/delayed animations with many molecules */
void set_display_only_model_mol(int imol) {

   graphics_info_t g;
   g.undisplay_all_model_molecules_except(imol);
   graphics_draw();
}



/*! \brief make the coordinates molecule active/inactve (clickable), 0
  for off, 1 for on */
void set_mol_active(int imol, int state) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_mol_is_active(state);
      set_display_control_button_state(imol, "Active", state);
      // graphics_draw(); // was this needed?
   } else {
      std::cout << "not valid molecule" << std::endl;
   }
}



int mol_is_displayed(int imol) {

   graphics_info_t g;
   return g.molecules[imol].draw_it;
}

int mol_is_active(int imol) {

   graphics_info_t g;

   if (is_valid_model_molecule(imol)) {
      return g.molecules[imol].atom_selection_is_pickable();
   }
   return 0;
}

int map_is_displayed(int imol) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t g;
      return g.molecules[imol].is_displayed_p();
   }
   return 0; // -1 maybe? No, I think not.
}

/*! \brief if on_or_off is 0 turn off all maps displayed, for other
  values of on_or_off turn on all maps */
void set_all_maps_displayed(int on_or_off) {

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      g.mol_displayed_toggle_do_redraw = false;
      int nm = graphics_info_t::n_molecules();
      for (int imol=0; imol<nm; imol++) {
         if (is_valid_map_molecule(imol)) {
            graphics_info_t::molecules[imol].set_map_is_displayed(on_or_off);
            set_display_control_button_state(imol, "Displayed", on_or_off);
         }
      }
      g.mol_displayed_toggle_do_redraw = true;
      graphics_draw();
   }
}

/*! \brief if on_or_off is 0 turn off all models displayed and active,
  for other values of on_or_off turn on all models. */
void set_all_models_displayed_and_active(int on_or_off) {

   // this could/should be moved into graphics_info_t.

   graphics_info_t g;
   g.mol_displayed_toggle_do_redraw = false;
   int nm = graphics_info_t::n_molecules();
   for (int imol=0; imol<nm; imol++) {
      if (is_valid_model_molecule(imol)) {
	      graphics_info_t::molecules[imol].set_mol_is_active(on_or_off);
	      graphics_info_t::molecules[imol].set_mol_is_displayed(on_or_off);
	      set_display_control_button_state(imol, "Active",    on_or_off);
	      set_display_control_button_state(imol, "Displayed", on_or_off);
      }
   }
   g.mol_displayed_toggle_do_redraw = true;
   graphics_draw();
}

/*\brief display only the active mol and the refinement map */
void display_only_active() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();

   // std::cout << "INFO:: display_only_active()" << aa.first << " " << aa.second.first << " " << aa.second.second << std::endl;
   logger.log(log_t::INFO, "display_only_active()", aa.first, aa.second.first, aa.second.second.format());

   if (aa.first) {
      int imol_active = aa.second.first;
      if (is_valid_model_molecule(imol_active)) {
	 coot::atom_spec_t atom_spec = aa.second.second;
	 int nm = graphics_info_t::n_molecules();
	 for (int imol=0; imol<nm; imol++) {
	    if (is_valid_model_molecule(imol)) {
	       if (imol == imol_active) {
		  graphics_info_t::molecules[imol].set_mol_is_active(true);
		  graphics_info_t::molecules[imol].set_mol_is_displayed(true);
	       } else {
		  graphics_info_t::molecules[imol].set_mol_is_displayed(false);
		  graphics_info_t::molecules[imol].set_mol_is_active(false);
	       }
	    }
	 }
      }
   }
   graphics_draw();
}

void set_only_last_model_molecule_displayed() {

   int n_mols = graphics_info_t::n_molecules();
   int imol_last = -1;
   graphics_info_t g;
   std::vector<int> turn_these_off; // can contain last
   for (int i=0; i<n_mols; i++) {
      if (is_valid_model_molecule(i)) {
	 if (mol_is_displayed(i)) {
	    turn_these_off.push_back(i);
	 }
	 imol_last = i; // updates
      }
   }

   g.mol_displayed_toggle_do_redraw = false;

   for (unsigned int j=0; j<turn_these_off.size(); j++) {
      if (turn_these_off[j] != imol_last) {

         // These do a redraw
         // set_mol_displayed(turn_these_off[j], 0);
         // set_mol_active(turn_these_off[j], 0);

         // std::cout << ".....  turning off button for " << turn_these_off[j] << std::endl;

         g.molecules[turn_these_off[j]].set_mol_is_displayed(0);
         g.molecules[turn_these_off[j]].set_mol_is_active(0);
         set_display_control_button_state(turn_these_off[j], "Displayed", 0);
         set_display_control_button_state(turn_these_off[j], "Active", 0);
      }
   }
   if (is_valid_model_molecule(imol_last)) {
      if (! mol_is_displayed(imol_last)) {

	 // set_mol_displayed(imol_last, 1);
	 // set_mol_active(imol_last, 1);

	 g.molecules[imol_last].set_mol_is_displayed(1);
	 g.molecules[imol_last].set_mol_is_active(1);
         set_display_control_button_state(imol_last, "Displayed", 1);
      }
   }
   g.mol_displayed_toggle_do_redraw = true; // back on again
   graphics_draw();

}



// Bleugh.
char *
show_spacegroup(int imol) {

   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      std::string spg = graphics_info_t::molecules[imol].show_spacegroup();
      // std::cout << "INFO:: spacegroup: " << spg << std::endl;
      logger.log(log_t::INFO, "spacegroup:", spg);
      unsigned int l = spg.length();
      char *s = new char[l+1];
      strncpy(s, spg.c_str(), l+1);
      return s;
   } else {

      // If it was a bad molecule, return pointer to null.
      std::cout << "Unknown molecule " << imol << std::endl;
      char *s = new char[1];
      s[0] = 0;
      return s;
   }
}

#ifdef USE_GUILE
/*! \brief return the spacegroup as a string, return scheme false if unable to do so. */
SCM space_group_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      std::string s =  graphics_info_t::molecules[imol].show_spacegroup();
      r = scm_from_locale_string(s.c_str());
   }
   return r;

}
#endif

#ifdef USE_PYTHON
PyObject *space_group_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      std::string s =  graphics_info_t::molecules[imol].show_spacegroup();
      r = myPyString_FromString(s.c_str());
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif





#ifdef USE_GUILE
/*! \brief return a list of symmetry operators as strings - or scheme false if
  that is not possible. */
SCM symmetry_operators_scm(int imol) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      std::pair<bool, clipper::Spacegroup> sg =
	 graphics_info_t::molecules[imol].space_group();
      if (! sg.second.is_null()) {
	 r = SCM_EOL;
	 std::vector<std::string> sv =
	    graphics_info_t::molecules[imol].get_symop_strings();
	 for (unsigned int i=0; i<sv.size(); i++) {
	    SCM s = scm_from_locale_string(sv[i].c_str());
	    r = scm_cons(s, r);
	 }
	 r = scm_reverse(r);
      } else {
	 std::cout << "WARNING:: in symmetry_operators_scm() null space group " << std::endl;
      }
   }
   return r;
}
#endif


#ifdef USE_PYTHON
/*! \brief return a list of symmetry operators as strings - or python false if
  that is not possible. */
PyObject *symmetry_operators_py(int imol) {

   PyObject *o = Py_False;
   if (is_valid_model_molecule(imol) || is_valid_map_molecule(imol)) {
      std::pair<bool, clipper::Spacegroup> sg =
	 graphics_info_t::molecules[imol].space_group();
      if (! sg.second.is_null()) {
	 std::vector<std::string> sv =
	    graphics_info_t::molecules[imol].get_symop_strings();
	 o = PyList_New(sv.size());
	 for (unsigned int i=0; i<sv.size(); i++) {
	    PyList_SetItem(o, i, myPyString_FromString(sv[i].c_str()));
	 }
      } else {
	 std::cout << "WARNING:: in symmetry_operators_py() null space group " << std::endl;
      }
   }
   if (PyBool_Check(o)) {
     Py_INCREF(o);
   }
   return o;
}
#endif


#ifdef USE_GUILE
SCM
symmetry_operators_to_xHM_scm(SCM symmetry_operators) {
   SCM r = SCM_BOOL_F;
   clipper::Spacegroup sg = scm_symop_strings_to_space_group(symmetry_operators);
   if (! sg.is_null())
      r = scm_from_locale_string(sg.symbol_hm().c_str());
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *
symmetry_operators_to_xHM_py(PyObject *symmetry_operators) {
   PyObject *o = Py_False;
   clipper::Spacegroup sg = py_symop_strings_to_space_group(symmetry_operators);
   if (! sg.is_null())
      o = myPyString_FromString(sg.symbol_hm().c_str());
   if (PyBool_Check(o)) {
     Py_INCREF(o);
   }
   return o;
}
#endif



/*  ----------------------------------------------------------------------- */
/*                  zoom                                                    */
/*  ----------------------------------------------------------------------- */

void
scale_zoom_internal(float f) {

   // now we filter out unusual/erroneous changes in zoom

   graphics_info_t g;
   if (f > 0.0)
      if (f < 1.8)
         if (f >= 0.5)
            g.zoom *= f;

}

void scale_zoom(float f) {

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

void set_zoom(float f) {
   graphics_info_t::zoom = f;
   graphics_draw();
}





void clear_moving_atoms_object() {  /* redraw done here. */

   graphics_info_t g;
   g.clear_moving_atoms_object();

}

void clear_up_moving_atoms() {

   graphics_info_t g;
   // std::cout << "c-interface clear_up_moving_atoms..." << std::endl;
   g.clear_up_moving_atoms();
   g.clear_moving_atoms_object();
}

// this is a better function name
void set_refine_ramachandran_torsion_angles(int state) {
   set_refine_ramachandran_angles(state);
}


// either alpha helix, beta strand or ramachandran goodness
// (see ideal/simple_restraint.hh link torsions)
void set_refine_ramachandran_angles(int state) {
   graphics_info_t::do_rama_restraints = state;

   // Adjust the GUI
   if (graphics_info_t::use_graphics_interface_flag) {
      std::string w_name = "main_toolbar_restraints_rama_label";
      // GtkWidget *w = lookup_widget(graphics_info_t::glareas[0], w_name.c_str());
      GtkWidget *w = widget_from_builder(w_name);
      if (w) {
	 if (state) {
	    if (graphics_info_t::restraints_rama_type == coot::RAMA_TYPE_ZO) {
	       // this is not antialiased on my mac - yuck.
	       // std::string l = "<span background=\"white\" foreground=\"brown\" face=\"PilGi\"
	       // size=\"larger\">Rama</span>";
	       std::string l = "<span background=\"white\" foreground=\"brown\">Rama</span>";
	       gtk_label_set_markup(GTK_LABEL(w), l.c_str());
	    }
	    gtk_widget_set_visible(w, TRUE);
	 } else {
	    gtk_widget_set_visible(w, FALSE);
	 }
      }
   }
}


void set_refine_ramachandran_restraints_type(int type) {
   graphics_info_t::restraints_rama_type = type;
   if (type == 0)
      graphics_info_t::rama_plot_restraints_weight = 1.0; // big numbers make the refinement fail (precision?)
}


void set_refine_ramachandran_restraints_weight(float w) {
   graphics_info_t::rama_plot_restraints_weight = w;
}

float refine_ramachandran_restraints_weight() {
   return graphics_info_t::rama_plot_restraints_weight;
}



int refine_ramachandran_angles_state() {
   return graphics_info_t::do_rama_restraints;
}

/* \brief set the weight for torsion restraints (default 1.0)*/
void set_torsion_restraints_weight(double w) {

   graphics_info_t::torsion_restraints_weight = w;
}



void set_refine_rotamers(int state) {
   graphics_info_t::do_rotamer_restraints = state;
}



#include "c-interface-refine.hh"

void set_refinement_geman_mcclure_alpha_from_text(int combobox_index, const char *t) {

   graphics_info_t g;
   float v = coot::util::string_to_float(t);
   set_refinement_geman_mcclure_alpha(v);
   graphics_info_t::refine_params_dialog_geman_mcclure_alpha_combobox_position = combobox_index;
   // poke the refinement if there are moving atoms
}

void set_refinement_lennard_jones_epsilon_from_text(int idx, const char *t) {

   graphics_info_t g;
   float v = coot::util::string_to_float(t);
   set_refinement_lennard_jones_epsilon(v);
   graphics_info_t::refine_params_dialog_lennard_jones_epsilon_combobox_position = idx;
   std::cout << "############################ refine_params_dialog_lennard_jones_epsilon_combobox_position set "
             << idx << std::endl;
   g.poke_the_refinement();
}

void set_refinement_ramachandran_restraints_weight_from_text(int idx, const char *t) {

   float v = coot::util::string_to_float(t);
   set_refine_ramachandran_restraints_weight(v);
   graphics_info_t::refine_params_dialog_rama_restraints_weight_combobox_position = idx;
   graphics_info_t g;
   g.poke_the_refinement();
}

void set_refinement_overall_weight_from_text(const char *t) {

   if (t) {
      float v = coot::util::string_to_float(t);
      graphics_info_t::geometry_vs_map_weight = v;
      graphics_info_t g;
      g.poke_the_refinement();

   } else {
      std::cout << "WARNING:: in set_refinement_overall_weight_from_text() t null " << std::endl;
   }

}

void set_refinement_torsion_weight_from_text(int idx, const char *t) {

   graphics_info_t g;
   float v = coot::util::string_to_float(t);
   graphics_info_t::refine_params_dialog_torsions_weight_combox_position = idx;
   graphics_info_t::torsion_restraints_weight = v;
   // poke the refinement if there are moving atoms
   g.poke_the_refinement();
}


void set_refine_params_dialog_more_control_frame_is_active(int state) {

   graphics_info_t::refine_params_dialog_extra_control_frame_is_visible = state;
}


void set_dragged_refinement_steps_per_frame(int v) {

   graphics_info_t g;
   g.dragged_refinement_steps_per_frame = v;
}

int dragged_refinement_steps_per_frame() {
   return graphics_info_t::dragged_refinement_steps_per_frame;
}




#ifdef USE_GUILE
SCM safe_scheme_command(const std::string &scheme_command) {
   std::vector<std::string> cs;
   cs.push_back(DIRECT_SCM_STRING);
   cs.push_back(scheme_command);
   add_to_history(cs);
   return graphics_info_t::safe_scheme_command(scheme_command);
}
#endif // USE_GUILE


#ifdef USE_GUILE
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
   std::cout << "debug:: safe_python_command() PyRun_SimpleString() " << python_cmd << std::endl;
   PyRun_SimpleString(python_cmd.c_str());
#endif
}


#ifdef USE_PYTHON
// We need a function to clean up the returned types from safe_python_command_with_return
// especially lists and floats. Who knows why? Maybe a Python bug!
// 'local function' currently
PyObject *py_clean_internal(PyObject *o) {

   PyObject *ret = NULL;
   if (PyList_Check(o)) {
     int l = PyObject_Length(o);
      ret = PyList_New(0);
      for (int item=0; item<l; item++) {
	 PyObject *py_item = PyList_GetItem(o, item);
	 py_item = py_clean_internal(py_item);
	 if (py_item == NULL) {
	   PyErr_Print();
	 }
	 PyList_Append(ret, py_item);
	 // Py_XDECREF(py_item); no!
      }
   } else {
      if (PyBool_Check(o)) {
	 // apparently doesnt seem to need resetting
	 int i = PyLong_AsLong(o);
	 ret = o;
      } else {
	 if (PyLong_Check(o)) {
	    // apparently doesnt seem to need resetting
	    int i=PyLong_AsLong(o);
	    ret = o;
	 } else {
	    if (PyFloat_Check(o)) {
	       // float is a copy
	       double f = PyFloat_AsDouble(o);
	       ret = PyFloat_FromDouble(f);
	    } else {
	       if (PyUnicode_Check(o)) {
		  ret = o;
	       } else {
		  if (PyFunction_Check(o)) {
		     ret = PyObject_Str(o);
		  } else {
		     if (o == Py_None) {
			//std::cout << "BL DEBUG:: have PyNone, not sure what to do with it!?" <<std::endl;
			ret = o;
		     } else {
			std::cout <<"WARNING:: py_clean_internal: incomprehensible argument passed  "
				  << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyObject_Str(o)))
                                  <<std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return ret;
}
#endif  // USE_PYTHON

int pyrun_simple_string(const char *python_command) {

#ifdef USE_PYTHON
  return PyRun_SimpleString(python_command);
#endif
  return -1;
}

#ifdef USE_PYTHON

/**
 * Alternative version that extracts a return value from multi-line code.
 * Wraps the code to capture the last expression's value.
 *
 * @param code Python code to execute
 * @return PyObject* result (new reference), or NULL on error
 */
std::pair<PyObject*, std::string>
execute_python_code_with_result_internal_old(const std::string &code) {

   std::string error_message;
   PyObject* result = NULL;

   // Get __main__ namespace
   PyObject* main_module = PyImport_AddModule("__main__");
   if (!main_module) {
      std::cerr << "ERROR: Failed to get __main__ module" << std::endl;
      return std::make_pair(result, error_message);
   }

   PyObject* global_dict = PyModule_GetDict(main_module);
   if (!global_dict) {
      std::cerr << "ERROR: Failed to get __main__ dictionary" << std::endl;
      return std::make_pair(result, error_message);
   }

   // Try as simple expression first
   result = PyRun_String(code.c_str(), Py_eval_input, global_dict, global_dict);
   if (result) {
      return std::make_pair(result, error_message);
   }

   PyErr_Clear();

   // For multi-line code, we need to handle it differently
   // Execute the code
   PyObject* exec_result = PyRun_String(code.c_str(), Py_file_input, global_dict, global_dict);

   if (!exec_result) {
        std::cerr << "ERROR: Python execution failed:" << std::endl;
        if (PyErr_Occurred()) {
           PyObject *ptype, *pvalue, *ptraceback;
           PyErr_Fetch(&ptype, &pvalue, &ptraceback);
           if (pvalue) {
              PyObject* str_obj = PyObject_Str(pvalue);
              if (str_obj) {
                 const char* str = PyUnicode_AsUTF8(str_obj);
                 if (str) {
                     error_message = str;
                 }
                 Py_DECREF(str_obj);
               }
           }

           Py_XDECREF(ptype);
           Py_XDECREF(pvalue);
           Py_XDECREF(ptraceback);
           PyErr_Print();
        }
        return std::make_pair(result, error_message);
   }

   Py_DECREF(exec_result);  // This is typically Py_None

    // Look for a special variable '__result__' that the code may have set
    // Or just return None to indicate successful execution
    PyObject* stored_result = PyDict_GetItemString(global_dict, "__result__");
    if (stored_result) {
        Py_INCREF(stored_result);  // GetItemString returns borrowed reference
        return std::make_pair(stored_result, error_message);
    }

    // No specific result - return None to indicate success
    // Py_RETURN_NONE;
    return std::make_pair(Py_None, error_message);
}

#include "python-results-container.hh"

execute_python_results_container_t execute_python_code_with_result_internal(const std::string &code) {

   // Try as simple (one-line) expression

   execute_python_results_container_t rc;
   // Get __main__ namespace
   PyObject* main_module = PyImport_AddModule("__main__");
   if (!main_module) {
      std::cerr << "ERROR: Failed to get __main__ module" << std::endl;
      rc.error_message = "Failed to get __main__ module";
      return rc;
   }
   PyObject* global_dict = PyModule_GetDict(main_module);
   if (!global_dict) {
      // std::cerr << "ERROR: Failed to get __main__ dictionary" << std::endl;
      // return std::make_pair(result, error_message);
      std::cerr << "ERROR: Failed to get __main__ dictionary" << std::endl;
      rc.error_message = "Failed to get __main__ dictionary";
      return rc;
   }
   // capture stdout - start
   // Save original stdout and redirect to StringIO
   PyRun_SimpleString("import sys, io");
   PyRun_SimpleString("__original_stdout__ = sys.stdout");
   PyRun_SimpleString("__stdout_capture__ = io.StringIO()");
   PyRun_SimpleString("sys.stdout = __stdout_capture__");
   // capture stdout - end

   PyObject *exec_result = PyRun_String(code.c_str(), Py_eval_input, global_dict, global_dict);
   rc.result = exec_result;
   if (exec_result) {
      // get captured output
      PyObject* stdout_obj = PyDict_GetItemString(global_dict, "__stdout_capture__");
      if (stdout_obj) {
         PyObject* captured = PyObject_CallMethod(stdout_obj, "getvalue", NULL);
         if (captured) {
            const char* output_str = PyUnicode_AsUTF8(captured);
            if (output_str && strlen(output_str) > 0) {
               std::cout << output_str;  // Print to terminal
               rc.stdout = output_str;
            }
            Py_DECREF(captured);
         }
      }
      // Restore stdout
      PyRun_SimpleString("sys.stdout = __original_stdout__");

   } else {
      std::cerr << "ERROR: execute_python_code_with_result_internal(): Python execution failed" << std::endl;
      // Import traceback module and format the exception
      PyObject *ptype, *pvalue, *ptraceback;
      PyErr_Fetch(&ptype, &pvalue, &ptraceback);

      std::cerr << "DEBUG: ptype=" << (ptype ? "NOT NULL" : "NULL") << std::endl;
      std::cerr << "DEBUG: pvalue=" << (pvalue ? "NOT NULL" : "NULL") << std::endl;
      std::cerr << "DEBUG: ptraceback=" << (ptraceback ? "NOT NULL" : "NULL") << std::endl;

      // Always use fallback path - simpler and more robust
      if (pvalue) {
         std::cerr << "DEBUG: Converting pvalue to string" << std::endl;
         PyObject* str_obj = PyObject_Str(pvalue);
         if (str_obj) {
            const char* str = PyUnicode_AsUTF8(str_obj);
            if (str) {
               rc.error_message = std::string(str);
               std::cerr << "DEBUG: rc.error_message set to: " << rc.error_message << std::endl;
            } else {
               rc.error_message = "Python error occurred but could not retrieve error message";
               std::cerr << "DEBUG: PyUnicode_AsUTF8 failed" << std::endl;
            }
            Py_DECREF(str_obj);
         } else {
            rc.error_message = "Python error occurred but PyObject_Str failed";
            std::cerr << "DEBUG: PyObject_Str failed" << std::endl;
         }
      } else {
         rc.error_message = "Python error occurred but no exception value available";
         std::cerr << "DEBUG: No pvalue" << std::endl;
      }

      // ALWAYS clean up the error objects
      Py_XDECREF(ptype);
      Py_XDECREF(pvalue);
      Py_XDECREF(ptraceback);

      // Restore stdout even on error path
      PyRun_SimpleString("sys.stdout = __original_stdout__");
      PyErr_Clear();
   }
   return rc;
}

execute_python_results_container_t execute_python_multiline_code_with_result_internal(const std::string &code) {

   execute_python_results_container_t rc;

   // Get __main__ namespace
   PyObject* main_module = PyImport_AddModule("__main__");
   if (!main_module) {
      std::cerr << "ERROR: Failed to get __main__ module" << std::endl;
      rc.error_message = "Failed to get __main__ module";
      return rc;
   }
   PyObject* global_dict = PyModule_GetDict(main_module);
   if (!global_dict) {
      // std::cerr << "ERROR: Failed to get __main__ dictionary" << std::endl;
      // return std::make_pair(result, error_message);
      std::cerr << "ERROR: Failed to get __main__ dictionary" << std::endl;
      rc.error_message = "Failed to get __main__ dictionary";
      return rc;
   }

   // capture stdout - start
   // Save original stdout and redirect to StringIO
   PyRun_SimpleString("import sys, io");
   PyRun_SimpleString("__original_stdout__ = sys.stdout");
   PyRun_SimpleString("__stdout_capture__ = io.StringIO()");
   PyRun_SimpleString("sys.stdout = __stdout_capture__");
   // capture stdout - end

   PyObject *exec_result = PyRun_String(code.c_str(), Py_file_input, global_dict, global_dict);
   std::cout << "DEBUG:: ------------ exec_result " << exec_result << std::endl;
   if (exec_result) {
      rc.result = exec_result;
      // get captured output
      PyObject* stdout_obj = PyDict_GetItemString(global_dict, "__stdout_capture__");
      if (stdout_obj) {
         PyObject* captured = PyObject_CallMethod(stdout_obj, "getvalue", NULL);
         if (captured) {
            const char* output_str = PyUnicode_AsUTF8(captured);
            if (output_str && strlen(output_str) > 0) {
               std::cout << output_str;  // Print to terminal
               rc.stdout = output_str;
            }
            Py_DECREF(captured);
         }
      }
      // Restore stdout
      PyRun_SimpleString("sys.stdout = __original_stdout__");

   } else {

      std::cerr << "ERROR: execute_python_multiline_code_with_result_internal(): Python execution failed" << std::endl;

      // Import traceback module and format the exception
      PyObject *ptype, *pvalue, *ptraceback;
      // PyErr_Fetch() consumes the error.
      PyErr_Fetch(&ptype, &pvalue, &ptraceback); // 2026-01-11-PE use PyErr_GetRaisedException() in future
      PyObject* traceback_module = PyImport_ImportModule("traceback");
      if (traceback_module && ptraceback) {
         PyObject* format_exception = PyObject_GetAttrString(traceback_module, "format_exception");
         if (format_exception) {
             PyObject* formatted = PyObject_CallFunctionObjArgs(format_exception, ptype, pvalue, ptraceback, NULL);
             if (formatted && PyList_Check(formatted)) {
                 // Join all traceback lines into single string
                 std::string full_traceback;
                 Py_ssize_t size = PyList_Size(formatted);
                 for (Py_ssize_t i = 0; i < size; i++) {
                     PyObject* line = PyList_GetItem(formatted, i);
                     const char* line_str = PyUnicode_AsUTF8(line);
                     if (line_str) {
                         full_traceback += line_str;
                     }
                 }
                 rc.error_message = full_traceback;
             }
             Py_XDECREF(formatted);
             Py_DECREF(format_exception);
         }
         Py_DECREF(traceback_module);

      } else {
         // Fallback to simple error message if traceback unavailable
         if (pvalue) {
            PyObject* str_obj = PyObject_Str(pvalue);
            if (str_obj) {
               const char* str = PyUnicode_AsUTF8(str_obj);
               if (str) {
                   rc.error_message = str;
               }
               Py_DECREF(str_obj);
            }
         }

         Py_XDECREF(ptype);
         Py_XDECREF(pvalue);
         Py_XDECREF(ptraceback);
      }
      return rc;
   }

   // **GET CAPTURED OUTPUT (for multi-line case too)**
   PyObject* stdout_obj = PyDict_GetItemString(global_dict, "__stdout_capture__");
   if (stdout_obj) {
      PyObject* captured = PyObject_CallMethod(stdout_obj, "getvalue", NULL);
      if (captured) {
         const char* output_str = PyUnicode_AsUTF8(captured);
         if (output_str && strlen(output_str) > 0) {
            std::cout << output_str;  // Print to terminal
            // You could also save to a file here if needed
            rc.stdout = output_str;
         }
         Py_DECREF(captured);
      }
   }
   // Restore stdout
   PyRun_SimpleString("sys.stdout = __original_stdout__");

   // Look for a special variable '__result__' that the code may have set
   PyObject* stored_result = PyDict_GetItemString(global_dict, "__result__");
   if (stored_result) {
        Py_INCREF(stored_result);
        // return std::make_pair(stored_result, error_message);
        rc.result = stored_result;
        return rc;
   }

   // return std::make_pair(Py_None, error_message);
   rc.result = Py_None;
   return rc;
}


PyObject* execute_python_code_with_result(const std::string &code) {

   // std::pair<PyObject *, std::string> r = execute_python_code_with_result_internal(code);
   execute_python_results_container_t r = execute_python_code_with_result_internal(code);
   return r.result;
}


// BL says:: let's have a python command with can receive return values
// we need to pass the script file containing the funcn and the funcn itself
// returns a PyObject which can then be used further
// returns NULL for failed run
PyObject *safe_python_command_with_return(const std::string &python_cmd) {

   std::cout << "--------------- start safe_python_command_with_return(): " << python_cmd << std::endl;

   // 20220330-PE I think that this is super ricketty now!
   // Does it only find things in dynamic_atom_overlaps_and_other_outliers module?
   // this function was empty before today, returning NULL.

   // std::cout << "in safe_python_command_with_return() A " << python_cmd << std::endl;

   // command = "import coot; " + python_cmd;
   std::string command = python_cmd;

   PyObject* result = nullptr;
   PyObject *am = PyImport_AddModule("__main__");

   if (am) {
      PyObject* d = PyModule_GetDict(am);

      const char *modulename = "coot";
      PyObject *pName = myPyString_FromString(modulename);
      PyObject *pModule_coot = PyImport_Import(pName);

      std::cout << "running command: " << command << std::endl;
      PyObject* source_code = Py_CompileString(command.c_str(), "adhoc", Py_eval_input);
      if (source_code) {
         PyObject* func = PyFunction_New(source_code, d);
         result = PyObject_CallObject(func, PyTuple_New(0));
         std::cout << "--------------- in safe_python_command_with_return() result at: " << result << std::endl;
         if (result) {
            if(!PyUnicode_Check(result)) {
                std::cout << "--------------- in safe_python_command_with_return() result is probably not a string." << std::endl;
            }
            PyObject* displayed = display_python(result);
            PyObject* as_string = PyUnicode_AsUTF8String(displayed);
            std::cout << "--------------- in safe_python_command_with_return() result: "
                      << PyBytes_AS_STRING(as_string) << std::endl;
            Py_XDECREF(displayed);
            Py_XDECREF(as_string);
         }
         else {
            std::cout << "--------------- in safe_python_command_with_return() result was null" << std::endl;
            if(PyErr_Occurred()) {
               std::cout << "--------------- in safe_python_command_with_return() Printing Python exception:" << std::endl;
               PyErr_Print();
            }
         }

         // debugging
         // PyRun_String("import coot; print(dir(coot))", Py_file_input, d, d);
         Py_XDECREF(func);
         Py_XDECREF(source_code);
      } else {
         std::cout << "DEBUG:: in safe_python_command_with_return, null source_code" << std::endl;
      }
   } else {
      std::cout << "ERROR:: Hopeless failure: module for __main__ is null" << std::endl;
   }
   // 20230605-PE frustratingly this is returning None when I hope/expect it to be True.
   std::cout << "--------------- done safe_python_command_with_return() " << python_cmd << std::endl;
   return result;
}
#endif //PYTHON

#ifdef USE_PYTHON
PyObject *safe_python_command_test(const char *cmd) {

   std::string s = cmd;
   return safe_python_command_with_return(s);
}
#endif //PYTHON

#ifdef USE_PYTHON
/* BL says:: try not to use this!!!! */
/* doesnt work anyway!? */
void safe_python_command_with_unsafe_thread(const char *cmd) {

  Py_BEGIN_ALLOW_THREADS;
  PyRun_SimpleString((char *)cmd);
  Py_END_ALLOW_THREADS;
}

#endif //PYTHON

void safe_python_command_by_char_star(const char *python_cmd) {

#ifdef USE_PYTHON
   PyRun_SimpleString((char *)python_cmd);
#endif
}

// functions to run python commands in guile and vice versa
// ignoring return values for now
#ifdef USE_PYTHON
PyObject *run_scheme_command(const char *cmd) {

  PyObject *ret_py = Py_None;

#ifdef USE_GUILE
  SCM ret_scm = safe_scheme_command(cmd);
  ret_py = scm_to_py(ret_scm);
#endif // USE_GUILE

  if (PyBool_Check(ret_py) || ret_py == Py_None) {
    Py_INCREF(ret_py);
  }
  return ret_py;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM run_python_command(const char *python_cmd) {

   // SCM_UNSPECIFIED corresponds to Py_None
   SCM ret_scm = SCM_UNSPECIFIED; // not SCM_UNDEFINED; :)

#ifdef USE_PYTHON
   PyObject *ret = safe_python_command_with_return(python_cmd);
   if (ret != Py_None) {
     ret_scm = py_to_scm(ret);
   } else {
       Py_XDECREF(ret);
   }
   // 15092012 correct this (I think)
   // need to decref if Py_None I (strongly) believe. (increfed in
   // safe_python_command_with_return.
   // This does not crash. No idea why previous was...
   // Note: not all decrefs have increfs! Some increfs happens with python
   // API. There may be some more fixing required with respect to ref counting.
   // FIXME!

   // 20120904 comment this out
   // Py_XDECREF(ret);
   // as it causes a crash
   // e.g. open Extensions->Settings->Keybindings invokes this function
   // then running a python function (e.g. a key-binding) causes a crash.
   // I believe that Py_XDECREF(ret) is wrong here because we have never
   // called Py_INCREF(ref).  So perhaps we should be using both Py_INCREF()
   // and Py_XDECREF().

#endif // USE_PYTHON

  return ret_scm;
}
#endif // USE_GUILE

#ifdef USE_GUILE
// Return a list describing a residue like that returned by
// residues-matching-criteria (list return-val chain-id resno ins-code)
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
// return a null list on problem
SCM residue_spec_to_scm(const coot::residue_spec_t &res) {
   SCM r = SCM_EOL;

//    std::cout <<  "residue_spec_to_scm on: " << res.chain << " " << res.resno << " "
// 	     << res.insertion_code  << std::endl;
   r = scm_cons(scm_from_locale_string(res.ins_code.c_str()), r);
   r = scm_cons(scm_from_int(res.res_no), r);
   r = scm_cons(scm_from_locale_string(res.chain_id.c_str()), r);
   r = scm_cons(SCM_BOOL_T, r);
   return r;
}
#endif

#ifdef USE_PYTHON
// Return a list describing a residue like that returned by
// residues-matching-criteria [return-val, chain-id, resno, ins-code]
// This is a library function really.  There should be somewhere else to put it.
// It doesn't need expression at the scripting level.
PyObject *residue_spec_to_py(const coot::residue_spec_t &res) {

   // 20251216-PE this no longer prefixes a True
   PyObject *r = PyList_New(3);

   PyList_SetItem(r, 0, myPyString_FromString(res.chain_id.c_str()));
   PyList_SetItem(r, 1, PyLong_FromLong(res.res_no));
   PyList_SetItem(r, 2, myPyString_FromString(res.ins_code.c_str()));

   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
int mark_multiple_atoms_as_fixed_scm(int imol, SCM atom_spec_list, int state) {

   // not tested
   int r = 0;
   SCM list_length_scm = scm_length(atom_spec_list);
   int n = scm_to_int(list_length_scm);

   for (int ispec = 0; ispec<n; ispec++) {
      SCM atom_spec_scm = scm_list_ref(atom_spec_list, scm_from_int(ispec));
      coot::atom_spec_t spec = atom_spec_from_scm_expression(atom_spec_scm);
      graphics_info_t::mark_atom_as_fixed(imol, spec, state);
   }

   if (n > 0) {
      graphics_draw();
   }

   return n; //return a count of how many atoms we successfully marked
}
#endif // USE_GUILE


#ifdef USE_PYTHON
// Guaranteed to return a triple list (will return unset-spec if needed).
//
PyObject *residue_spec_make_triple_py(PyObject *res_spec_py) {

   coot::residue_spec_t res_spec_default;
   PyObject *r = PyList_New(3);

   if (PyList_Check(res_spec_py)) {
      long l = PyObject_Length(res_spec_py);
      int offset = 0;
      if (l == 4) {
	 offset = 1;
      }
      PyObject *chain_id_py = PyList_GetItem(res_spec_py, offset);
      PyObject *res_no_py   = PyList_GetItem(res_spec_py, offset+1);
      PyObject *ins_code_py = PyList_GetItem(res_spec_py, offset+2);
      PyList_SetItem(r, 0, chain_id_py);
      PyList_SetItem(r, 1, res_no_py);
      PyList_SetItem(r, 2, ins_code_py);
   } else {
      PyList_SetItem(r, 0, myPyString_FromString(res_spec_default.chain_id.c_str()));
      PyList_SetItem(r, 1, PyLong_FromLong(res_spec_default.res_no));
      PyList_SetItem(r, 2, myPyString_FromString(res_spec_default.ins_code.c_str()));
   }
   return r;
}
#endif // USE_PYTHON


#ifdef USE_GUILE
// Return a SCM list object of (residue1 residue2 omega)
SCM cis_peptides(int imol) {
   SCM r = SCM_EOL;

   // more info on the real cis peptides derived from atom positions:

   if (is_valid_model_molecule(imol)) {

      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      std::vector<coot::util::cis_peptide_info_t> v =
	 coot::util::cis_peptides_info_from_coords(mol);

      for (unsigned int i=0; i<v.size(); i++) {
	 coot::residue_spec_t r1(v[i].chain_id_1,
				 v[i].resno_1,
				 v[i].ins_code_1);
	 coot::residue_spec_t r2(v[i].chain_id_2,
				 v[i].resno_2,
				 v[i].ins_code_2);
	 SCM scm_r1 = residue_spec_to_scm(r1);
	 SCM scm_r2 = residue_spec_to_scm(r2);
	 SCM scm_residue_info = SCM_EOL;
// 	 std::cout << "DEBUG:: cis pep with omega: "
// 		   << v[i].omega_torsion_angle
// 		   << std::endl;
// 	 SCM scm_omega =
// 	    scm_from_double(clipper::Util::rad2d(v[1].omega_torsion_angle));
	 SCM scm_omega = scm_from_double(v[i].omega_torsion_angle);
	 scm_residue_info = scm_cons(scm_omega, scm_residue_info);
	 scm_residue_info = scm_cons(scm_r2, scm_residue_info);
	 scm_residue_info = scm_cons(scm_r1, scm_residue_info);

	 // add scm_residue_info to r
	 r = scm_cons(scm_residue_info, r);
      }
      r = scm_reverse(r);
   }
   return r;
}
#endif //  USE_GUILE

#ifdef USE_GUILE
// Return a SCM list object of (residue1 residue2 omega)
SCM twisted_trans_peptides(int imol) {
   SCM r = SCM_EOL;

   // more info on the real cis peptides derived from atom positions:

   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      std::vector<coot::util::cis_peptide_quad_info_t> v =
	 coot::cis_peptide_quads_from_coords(mol, 0, g.Geom_p(), false);

      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].type == coot::util::cis_peptide_quad_info_t::TWISTED_TRANS) {
	    try {
	       coot::residue_spec_t r1(v[i].quad.atom_1->GetResidue());
	       coot::residue_spec_t r2(v[i].quad.atom_4->GetResidue());
	       SCM scm_r1 = residue_spec_to_scm(coot::residue_spec_t(r1));
	       SCM scm_r2 = residue_spec_to_scm(coot::residue_spec_t(r2));
	       double omega = v[i].quad.torsion();
	       SCM scm_omega = scm_from_double(omega);
	       SCM scm_residue_info = scm_list_3(scm_r1, scm_r2, scm_omega);

	       // add scm_residue_info to r
	       r = scm_cons(scm_residue_info, r);
	    }
	    catch (const std::runtime_error &e) {
	       std::cout << "WARNING:: " << e.what() << std::endl;
	    }
	 }
      }
      r = scm_reverse(r);
   }
   return r;
}
#endif //  USE_GUILE

#ifdef USE_PYTHON
// Return a python list object of [residue1, residue2, omega]
PyObject *cis_peptides_py(int imol) {
   PyObject *r;
   r = PyList_New(0);

   // more info on the real cis peptides derived from atom positions:

   if (is_valid_model_molecule(imol)) {

      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      std::vector<coot::util::cis_peptide_info_t> v =
	 coot::util::cis_peptides_info_from_coords(mol);

      for (unsigned int i=0; i<v.size(); i++) {
	 coot::residue_spec_t r1(v[i].chain_id_1,
				 v[i].resno_1,
				 v[i].ins_code_1);
	 coot::residue_spec_t r2(v[i].chain_id_2,
				 v[i].resno_2,
				 v[i].ins_code_2);
	 PyObject *py_r1, *py_r2, *py_residue_info;
	 py_r1 = residue_spec_to_py(r1);
	 py_r2 = residue_spec_to_py(r2);
	 py_residue_info = PyList_New(3);
// 	 std::cout << "DEBUG:: cis pep with omega: "
// 		   << v[i].omega_torsion_angle
// 		   << std::endl;
// 	 SCM scm_omega =
// 	    scm_from_double(clipper::Util::rad2d(v[1].omega_torsion_angle));
	 PyObject *py_omega;
	 py_omega = PyFloat_FromDouble(v[i].omega_torsion_angle);
	 PyList_SetItem(py_residue_info, 2, py_omega);
	 PyList_SetItem(py_residue_info, 1, py_r2);
	 PyList_SetItem(py_residue_info, 0, py_r1);

	 // add py_residue_info to r
	 PyList_Append(r, py_residue_info);
	 Py_XDECREF(py_residue_info);

      }
   }

   return r;
}
#endif //  USE_PYTHON

/*! \brief cis-trans convert the active residue of the active atom in the
    inermediate atoms, and continue with the refinement  */
int cis_trans_convert_intermediate_atoms() {

   graphics_info_t g;
   return g.cis_trans_conversion_intermediate_atoms();
}


#ifdef USE_PYTHON
// Return a python list object of [residue1, residue2, omega]
PyObject *twisted_trans_peptides_py(int imol) {

   PyObject *r;
   r = PyList_New(0);

   // more info on the real cis peptides derived from atom positions:

   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      std::vector<coot::util::cis_peptide_quad_info_t> v =
	 coot::cis_peptide_quads_from_coords(mol, 0, g.Geom_p(), false);

      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].type == coot::util::cis_peptide_quad_info_t::TWISTED_TRANS) {
	    try {
	       PyObject *py_r1, *py_r2, *py_residue_info;
	       coot::residue_spec_t r1(v[i].quad.atom_1->GetResidue());
               coot::residue_spec_t r2(v[i].quad.atom_4->GetResidue());
               py_r1 = residue_spec_to_py(coot::residue_spec_t(r1));
	       py_r2 = residue_spec_to_py(coot::residue_spec_t(r2));
	       py_residue_info = PyList_New(3);
	       PyObject *py_omega = PyFloat_FromDouble(v[i].quad.torsion());
	       PyList_SetItem(py_residue_info, 0, py_r1);
	       PyList_SetItem(py_residue_info, 1, py_r2);
	       PyList_SetItem(py_residue_info, 2, py_omega);

	       // add py_residue_info to r
	       PyList_Append(r, py_residue_info);
	       // Py_XDECREF(py_residue_info); // is py_residue_info copied? I doubt that this is needed.

	    }
	    catch (const std::runtime_error &e) {
	       std::cout << "WARNING:: " << e.what() << std::endl;
	    }
	 }
      }
   }

   return r;
}
#endif //  USE_PYTHON

void post_scripting_window() {

#ifdef USE_GUILE
   post_scheme_scripting_window();
#endif

}


/*! \brief pop-up a scripting window for scheming */
void post_scheme_scripting_window() {

#ifdef USE_GUILE

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
		  << "   Nevertheless you will get a simple scripting window."
		  << std::endl;
        // load the fallback window if we have COOT_SCHEME_DIR (only Windows?!)
        // only for gtk2!

        GtkWidget *scheme_entry;
        // GtkWidget *window = create_scheme_window();
        GtkWidget *window = widget_from_builder("scheme_window");

        // scheme_entry = lookup_widget(window, "scheme_window_entry");

        // 20220309-PE when will this get used?
        scheme_entry = widget_from_builder("scheme_window_entry");
        setup_guile_window_entry(scheme_entry); // USE_PYTHON and USE_GUILE used here
        gtk_widget_set_visible(window, TRUE);

     } else {
	std::cout << COOT_SCHEME_DIR << " was not defined - cannot open ";
	std::cout << "scripting window" << std::endl;
     }
  }
#endif

}

#include "network-get.hh"
// this needs a header?
// mode 1 means "mtz" and 0 means coords
void network_get_accession_code_entity(const std::string &text, int mode);

/* called from c-inner-main */
// If you edit this again, move it to graphics_info_t.
//
void
run_command_line_scripts() {

   if (false)
      std::cout << "--------------------- run_command_line_scripts() ----------------"
                << std::endl;

   if (graphics_info_t::command_line_scripts.size()) {
      // std::cout << "INFO:: There are " << graphics_info_t::command_line_scripts.size()
      //          << " command line scripts to run\n";
      logger.log(log_t::INFO, "There are", graphics_info_t::command_line_scripts.size(), "command line scripts to run");
      for (unsigned int i=0; i<graphics_info_t::command_line_scripts.size(); i++)
	 std::cout << "    " << graphics_info_t::command_line_scripts[i].c_str()
		   << std::endl;
   }

   // --------- scripts ------------

   for (unsigned int i=0; i<graphics_info_t::command_line_scripts.size(); i++) {
      std::string fn = graphics_info_t::command_line_scripts[i];
      std::cout << "calling run_script() for file " << fn << std::endl;
      run_script(fn.c_str());
   }

   // --------- commands ------------

   for (unsigned int i=0; i<graphics_info_t::command_line_commands.commands.size(); i++)
      if (graphics_info_t::command_line_commands.is_python)
         safe_python_command(graphics_info_t::command_line_commands.commands[i].c_str());
      else
         safe_scheme_command(graphics_info_t::command_line_commands.commands[i].c_str());

    for (unsigned int i=0; i<graphics_info_t::command_line_commands.commands.size(); i++)
       if (graphics_info_t::command_line_commands.is_python)
	  safe_python_command(graphics_info_t::command_line_commands.commands[i].c_str());
       else
	  safe_scheme_command(graphics_info_t::command_line_commands.commands[i].c_str());

   graphics_info_t g;
   for (unsigned int i=0; i<graphics_info_t::command_line_accession_codes.size(); i++) {
      const std::string &code = g.command_line_accession_codes[i];
      std::cout << "run_command_line_scripts(): get accession code " << code << std::endl;
      network_get_accession_code_entity(code, 0); // mode 0 means "not mtz"
      network_get_accession_code_entity(code, 1); // mtz mode
   }
}

void run_update_self_maybe() { // called when --update-self given at command line

#ifdef USE_LIBCURL
   if (graphics_info_t::update_self) {
#ifdef USE_GUILE
      safe_scheme_command("(update-self)");
#endif // USE_GUILE
#ifdef USE_PYTHON
      safe_python_command("update_self()");
#endif // USE_PYTHON
   }

#endif // USE_LIBCURL
}


void
set_guile_gui_loaded_flag() {

   graphics_info_t g;
   g.guile_gui_loaded_flag = TRUE;
}

void
set_python_gui_loaded_flag() {

   graphics_info_t g;
   g.python_gui_loaded_flag = TRUE;
}

void set_found_coot_gui() {

   graphics_info_t g;
#ifdef USE_GUILE
   std::cout << "Coot Scheme Scripting GUI code found and loaded." << std::endl;
   g.guile_gui_loaded_flag = TRUE;
#endif // USE_GUILE
}

void set_found_coot_python_gui() {

#ifdef USE_PYTHON
   graphics_info_t g;
   std::cout << "Coot Python Scripting GUI code found and loaded." << std::endl;
   g.python_gui_loaded_flag = TRUE;
#endif // USE_PYTHON

}

// return an atom index
int atom_spec_to_atom_index(int imol, const char *chain, int resno, const char *atom_name) {
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

      bool is_python = false;

      std::string::size_type ipy = fn.rfind(".py");
      if (ipy != std::string::npos) {
	 if (fn.substr(ipy) == ".py")
	    is_python = true;
      }

      if (is_python) {
	 run_python_script(filename);
      } else {
	 // std::cout << "DEBUG:: calling run_guile_script() " << filename << std::endl;
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
      graphics_info_t::state_file_was_run_flag = true;
   }
#else
#ifdef USE_PYTHON
   filename = "0-coot.state.py";
   struct stat buf;
   int status = stat(filename.c_str(), &buf);
   if (status == 0) {
      run_python_script(filename.c_str());
      graphics_info_t::state_file_was_run_flag = true;
   }
#endif
#endif
}

#ifdef USE_PYTHON
void
run_state_file_py() {
   std::string filename;
   filename = "0-coot.state.py";
   struct stat buf;
   int status = stat(filename.c_str(), &buf);
   if (status == 0) {
      run_python_script(filename.c_str());
      graphics_info_t::state_file_was_run_flag = true;
   }
}
#endif // USE_PYTHON

#include "pre-load.hh"
#ifdef HAVE_CXX_THREAD
#include <thread>
#endif

void pre_load_rotamer_tables() {

   // I don't have this right.  Deactivate for now.

#ifdef HAVE_CXX_THREAD
   // std::thread t(graphics_info_t::fill_rotamer_probability_tables);
   // t.join();
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
	    graphics_info_t::state_file_was_run_flag = true;
	 } else {
	    if (graphics_info_t::use_graphics_interface_flag) {
	       GtkWidget *dialog = wrapped_create_run_state_file_dialog(); // uses builder
               if (dialog)
                  gtk_widget_set_visible(dialog, TRUE);
               else
                  std::cout << "ERROR:: missing dialog" << std::endl;
	    }
	 }
      }
   }
}

GtkWidget *wrapped_create_run_state_file_dialog() {

#ifdef USE_GUILE
   std::string filename("0-coot.state.scm");
#else
// BL says:: we might want to have it in python too
#ifdef USE_PYTHON
   std::string filename("0-coot.state.py");
#else
   std::string filename("0-coot.state.scm"); // unusual
#endif // python
#endif // USE_GUILE
   short int il = 1;
   graphics_info_t g;
   // GtkWidget *w = create_run_state_file_dialog();
   // GtkWidget *vbox_mols = lookup_widget(w, "mols_vbox");

   GtkWidget *w = NULL;
   GtkWidget *vbox_mols = NULL;

   w = widget_from_builder("run_state_file_dialog");
   vbox_mols = widget_from_builder("mols_vbox");

   // std::cout << "########333333333333333 debug:: w " << w << std::endl;
   // std::cout << "########333333333333333 debug:: vbox_mols " << vbox_mols << std::endl;

   if (w) {
      // std::cout << "wrapped_create_run_state_file_dialog():: got widget w " << w << std::endl;
   } else {
      std::cout << "ERROR:: wrapped_create_run_state_file_dialog():: widget w was null " << std::endl;
   }

   std::vector<std::string> v = g.save_state_data_and_models(filename, il);
   for (unsigned int i=0; i<v.size(); i++) {
      std::string s = "    ";
      s += v[i];
      GtkWidget *label = gtk_label_new(s.c_str());
      // gtk_misc_set_alignment (GTK_MISC (label), 0.0, 0.5);

      gtk_label_set_xalign(GTK_LABEL(label), 0.0);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      gtk_box_append(GTK_BOX(vbox_mols), label);
#else
      gtk_box_pack_start(GTK_BOX(vbox_mols), label, FALSE, FALSE, 2);
#endif

      gtk_widget_set_visible(label, TRUE);
   }
   return w;
}

#ifdef USE_PYTHON
GtkWidget *wrapped_create_run_state_file_dialog_py() {

   std::string filename("0-coot.state.py");
   short int il = 1;
   // GtkWidget *w = create_run_state_file_dialog();
   GtkWidget *w = widget_from_builder("run_state_file_dialog");

   // GtkWidget *vbox_mols = lookup_widget(w, "mols_vbox");
   GtkWidget *vbox_mols = widget_from_builder("mols_vbox");

   graphics_info_t g;
   std::vector<std::string> v = g.save_state_data_and_models(filename, il);
   for (unsigned int i=0; i<v.size(); i++) {
      //       std::cout << "Got molecule: " << v[i] << std::endl;
      std::string s = "    ";
      s += v[i];
      GtkWidget *label = gtk_label_new(s.c_str());
      std::cout << "fix the alignment wrapped_create_run_state_file_dialog_py()" << std::endl;
      // gtk_misc_set_alignment (GTK_MISC (label), 0.0, 0.5);
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      gtk_box_append(GTK_BOX(vbox_mols), label);
#else
      gtk_box_pack_start(GTK_BOX(vbox_mols), label, FALSE, FALSE, 2);
#endif
      gtk_widget_set_visible(label, TRUE);
   }
   return w;
}
#endif // USE_PYTHON


void
run_guile_script(const char *filename) {

   std::cout << "debug:: run_guile_script() A on " << filename << std::endl;
#ifdef USE_GUILE
   std::string thunk("(lambda() ");
   thunk += "(load ";
   thunk += "\"";
   thunk += filename;
   thunk += "\"))";

   std::cout << "debug:: run_guile_script() B on " << filename << std::endl;
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
#if 0 // as it was for Python2
   std::string simple = "execfile(";
   simple += single_quote(s);
   simple += ")";
   std::cout << "Running python script " << s  << std::endl;
   PyRun_SimpleString(simple.c_str());
#endif

   if (coot::file_exists(filename_in)) {
      FILE *fp = fopen(filename_in, "r");
      PyRun_SimpleFile(fp, filename_in);
      fclose(fp);
   } else {
      std::cout << "WARNING:: in run_python_script() file " << filename_in
                << " does not exist" << std::endl;
   }

#endif // USE_PYTHON
}

int
import_python_module(const char *module_name, int use_namespace) {

   int err = 1;

#ifdef USE_PYTHON

   std::string simple;
   if (use_namespace) {
      simple = "import ";
      simple += module_name;
   } else {
      simple = "from ";
      simple += module_name;
      simple += " import *";
   }

   if (false)
      std::cout << "Importing python module " << module_name
                << " using command " << simple << std::endl;

   err = PyRun_SimpleString(simple.c_str());
#endif // USE_PYTHON
   return err;
}


// // 20230501-PE during merge, not sure this is a thing now.
// void add_on_rama_choices() {  // the the menu
//    GtkWidget* menu_item = widget_from_builder("ramachandran_plot1");
//    add_on_validation_graph_mol_options(menu_item, "ramachandran");
// }


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

   print_sequence_chain_general(imol, chain_id, 0, 0, "");
}

void print_sequence_chain_general(int imol, const char *chain_id,
                                   short int pir_format,
                                   short int file_output,
                                   const char *file_name) {

   std::string seq;
   bool with_spaces = false; // block spaced output is easier to read

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         if (std::string(chain_p->GetChainID()) == chain_id) {
            int nres = chain_p->GetNumberOfResidues();
            mmdb::PResidue residue_p;
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
                  if (with_spaces)
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

      std::string full_seq;
      if (pir_format) {
         std::string n = graphics_info_t::molecules[imol].name_sans_extension(0);
         full_seq = ">P1;";
         full_seq += n;
         full_seq += " ";
         full_seq += chain_id;
         full_seq += "\n\n";
         full_seq += seq;
         full_seq += "\n*\n";
      } else {
         std::string n = graphics_info_t::molecules[imol].name_sans_extension(0);
         full_seq = "> ";
         full_seq += n;
         full_seq += " ";
         full_seq += chain_id;
         full_seq += "\n";
         full_seq += seq;
         full_seq += "\n";
      }

      if (file_output) {
         std::ofstream f(file_name);
         if (f) {
            f << full_seq;
            f.close();
         } else {
            std::cout << "WARNING:: failed to open " << file_name << std::endl;
         }
      } else {
         std::cout << full_seq;
      }
   }
}

// the old name for the below function
void do_sequence_view(int imol) {
   sequence_view(imol);
}

// This is a gui function - move it to c-interface-gui.cc (and the above function)
void sequence_view(int imol) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

      GtkWidget *scrolled_window = gtk_scrolled_window_new();
      GtkWidget *frame = gtk_frame_new("");
      gtk_widget_set_hexpand(scrolled_window, TRUE);
      gtk_widget_set_vexpand(scrolled_window, TRUE);
      gtk_widget_set_hexpand(frame, TRUE);
      gtk_widget_set_vexpand(frame, TRUE);

      // The sequence-view is in the frame, the frame is in the scrolled window.
      // The scrolled window is in the overlay.
      // In GTK-overlay speak: the scrolled window is the overlay child
      // and the button is the overlay overlay.

      CootSequenceView *sv = coot_sequence_view_new();
      coot_sequence_view_set_structure(sv, imol, mol);

      gtk_frame_set_child(GTK_FRAME(frame), GTK_WIDGET(sv));
      gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(frame));

      {
         auto click_function = +[] (CootSequenceView* self, int imol, const coot::residue_spec_t &spec, gpointer userdata) {
            graphics_info_t g;
            g.go_to_residue(imol, spec);
            update_go_to_atom_from_current_position();
            g.graphics_grab_focus();
         };
         g_signal_connect(sv, "residue-clicked", G_CALLBACK(click_function), nullptr);
      }

      // now add a close button for that sequence view as an overlay of the sequence view vbox
      //
      GtkWidget *button = gtk_button_new_from_icon_name("window-close");
      GtkStyleContext *sc = gtk_widget_get_style_context(button);
      gtk_style_context_add_class(sc, "circular");
      auto close_button_callback = +[] (GtkButton *button, gpointer data) {
         int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(button), "imol"));
         std::cout << "close this sequence view " << imol << std::endl;
         GtkWidget *sequence_view_box = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "sequence_view_box"));
         // now find the child (the overlay) in the sequence_view_box and remove it

         GtkWidget *item_widget = gtk_widget_get_first_child(sequence_view_box);
         int n_children = 0;
         while (item_widget) {
            n_children++;
            GtkWidget *w = item_widget;
            item_widget = gtk_widget_get_next_sibling(item_widget);
            int imol_overlay = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "imol"));
            if (imol_overlay == imol) {
               gtk_box_remove(GTK_BOX(sequence_view_box), w);
               n_children--;
            }
         };
         // if the sequence view box no longer has children, then we should close up the pane
         if (n_children == 0) {
            GtkWidget *vbox = widget_from_builder("main_window_sequence_view_box");
            gtk_widget_set_visible(vbox, FALSE);
         }
      };
      g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(close_button_callback), NULL);
      g_object_set_data(G_OBJECT(button), "imol", GINT_TO_POINTER(imol));

      GtkWidget *overlay = gtk_overlay_new();
      gtk_overlay_set_child(GTK_OVERLAY(overlay), GTK_WIDGET(scrolled_window));
      gtk_overlay_add_overlay(GTK_OVERLAY(overlay), button);
      GtkWidget *vbox = widget_from_builder("main_window_sequence_view_box");
      gtk_widget_set_visible(vbox, TRUE);
      g_object_set_data(G_OBJECT(button), "sequence_view_box", vbox);
      g_object_set_data(G_OBJECT(overlay), "imol", GINT_TO_POINTER(imol));
      // GTK_ALIGN_END works OK/as intended, except the main graphics widget (or window) is too narrow to see it.
      // Make the window wider and change this to GTK_ALIGN_END.
      // gtk_widget_set_halign(GTK_WIDGET(button), GTK_ALIGN_START);
      gtk_widget_set_halign(GTK_WIDGET(button), GTK_ALIGN_END);
      gtk_widget_set_valign(GTK_WIDGET(button), GTK_ALIGN_START);

      gtk_box_append(GTK_BOX(vbox), overlay);

      // int new_height;
      // gtk_widget_measure(GTK_WIDGET(sv), GTK_ORIENTATION_VERTICAL, 0, &new_height, nullptr, nullptr, nullptr);
      // gtk_widget_set_size_request(vbox, -1, new_height);

      int minimum_size;
      int natural_size;
      gtk_widget_measure(GTK_WIDGET(sv), GTK_ORIENTATION_VERTICAL, 0, &minimum_size, &natural_size, nullptr, nullptr);
      int current_height = gtk_widget_get_height(vbox);
      if (current_height < natural_size) {
         gtk_widget_set_size_request(vbox, -1, natural_size);
      }
   }
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
   if (g.moving_atoms_displayed_p()) {
      g.add_status_bar_text("Edit Backbone is not available while moving atoms are active");
   } else {
      g.in_backbone_torsion_define = state;
      if (state) {
         std::cout << "click on an atom in the peptide to change" << std::endl;
         g.pick_cursor_maybe();
         g.pick_pending_flag = 1;
      } else {
         g.normal_cursor();
      }
   }
}

void setup_dynamic_distances(short int state) {

   graphics_info_t::in_dynamic_distance_define = state;
   if (state) {
      pick_cursor_maybe();
   }
}


void set_refine_with_torsion_restraints(int istate) {

   graphics_info_t::do_torsion_restraints = istate;

}


int refine_with_torsion_restraints_state() {

   return graphics_info_t::do_torsion_restraints;

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

   int r = -1;
      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

   // first, does the file exist?

   if (is_valid_model_molecule(imol_coordinates)) {
      struct stat s;
      int status = stat(filename, &s);
      // stat check the link targets not the link itself, lstat stats the
      // link itself.
      //
      if (status != 0 || !S_ISREG (s.st_mode)) {
	 // std::cout << "INFO:: Error reading " << filename << std::endl;
	 logger.log(log_t::INFO, "Error reading " + std::string(filename));
	 if (S_ISDIR(s.st_mode)) {
	    std::cout << filename << " is a directory." << std::endl;
	 }
	 return -1; // which is status in an error
      } else {
	 // std::cout << "INFO:: Reading cif file: " << filename << std::endl;
	 logger.log(log_t::INFO, "Reading cif file: " + std::string(filename));
	 graphics_info_t g;
	 int imol = g.create_molecule();
	 int istat =
	    g.molecules[imol].make_map_from_cif(imol, std::string(filename), imol_coordinates);

	 if (istat != -1) {
	    graphics_draw();
	 } else {
	    g.erase_last_molecule();
	    imol = -1;
	 }
	 return imol;

      }
   } else {
      std::cout << "WARNING:: " << imol_coordinates << " is not a valid model molecule" << std::endl;
   }
   return -1; // which is status in an error
}

// and make (and display) a 2fofc map.
//
// Pass the file name of the cif file and the molecule number for which
// we will calculate sfs.
//
int read_cif_data_2fofc_map(const char *filename, int imol_coordinates) {

   int imol = -1;
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
	 std::cout << filename << " is a directory." << std::endl;
      }
      return -1; // which is status in an error
   } else {


      if (is_valid_model_molecule(imol_coordinates)) {

	 // std::cout << "INFO:: Reading cif file: " << filename << std::endl;
	 logger.log(log_t::INFO, "Reading cif file: " + std::string(filename));

	 graphics_info_t g;

	 imol = g.create_molecule();
	 int istat = g.molecules[imol].make_map_from_cif_2fofc(imol,
							       std::string(filename),
							       imol_coordinates);
	 if (istat != -1) {
	    graphics_draw();
	 } else {
	    g.erase_last_molecule();
	    imol = -1;
	 }
      } else {
	 std::cout << "WARNING:: molecule " << imol_coordinates
		   << " is not a valid model molecule " << std::endl;
      }
   }
   return imol;
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
	 std::cout << filename << " is a directory." << std::endl;
      }
      return -1; // which is status in an error
   } else {

      std::cout << "Reading cif file: " << filename << std::endl;

      graphics_info_t g;

      int imol = g.create_molecule();
      int istat = g.molecules[imol].make_map_from_cif_fofc(imol,
							   std::string(filename),
							   imol_coordinates);

      if (istat != -1) {
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
   int imol = -1;

   // first, does the file exist?
   struct stat s;
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   //
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << std::endl;
      }
      return -1; // which is status in an error
   } else {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.
      std::string fn = filename;
      imol = g.create_molecule();
      int istat = g.molecules[imol].make_map_from_cif(imol, fn);
      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 g.activate_scroll_radio_button_in_display_manager(imol);
	 graphics_draw();
      } else {
	 g.erase_last_molecule();
	 imol = -1;
      }
   }
   return imol;
}

int read_cif_data_with_phases_diff_sigmaa(const char *filename) {

   graphics_info_t g;
   int imol = -1;

   // first, does the file exist?
   struct stat s;
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   //
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << std::endl;
      }
      return -1; // which is status in an error
   } else {

      std::cout << "Reading cif file (with phases - diff) : " << filename << std::endl;

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.
      std::string fn = filename;
      imol = g.create_molecule();
      int istat = g.molecules[imol].make_map_from_cif_diff_sigmaa(imol, fn);
      if (istat != -1) {
	 g.scroll_wheel_map = imol;
	 g.activate_scroll_radio_button_in_display_manager(imol);
	 graphics_draw();
      } else {
	 g.erase_last_molecule();
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
	 std::cout << filename << " is a directory." << std::endl;
      }
      return -1; // which is status in an error
   } else {

      // This function is the .cif equivalent of
      // c.f. read_phs_and_coords_and_make_map or make_and_draw_map,
      // map_fill_from_mtz.

      graphics_info_t g;

      int imol = g.create_molecule();
      std::string f(filename);
      short int swap_col = graphics_info_t::swap_difference_map_colours;

      int istat = g.molecules[imol].make_map_from_cif_nfofc(imol, f, map_type, swap_col);

      if (istat != -1) {

	 g.scroll_wheel_map = imol; // change the current scrollable map.
	 graphics_draw();
	 return imol;
      } else {
	 g.erase_last_molecule();
      }
      return -1; // error
   }
}

int handle_shelx_fcf_file_internal(const char *filename) {

   graphics_info_t g;
   std::vector<std::string> cmd;
   cmd.push_back("handle-shelx-fcf-file");
   cmd.push_back(single_quote(coot::util::intelligent_debackslash(filename)));

// #ifdef USE_GUILE
//    std::string s = g.state_command(cmd, coot::STATE_SCM);
//    safe_scheme_command(s);
// #endif

// #ifdef USE_PYTHON
// #ifndef USE_GUILE
//    std::string s = g.state_command(cmd, coot::STATE_PYTHON);
//    safe_python_command(s);
// #endif
// #endif

   return read_small_molecule_data_cif(filename);
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
int add_atom_label(int imol, const char *chain_id, int iresno, const char *atom_id) {

   int i = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      i = g.molecules[imol].add_atom_label(chain_id, iresno, atom_id);
   }
   return i;
}

int remove_atom_label(int imol, const char *chain_id, int iresno, const char *atom_id) {
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

void set_seg_ids_in_atom_labels(int istat) {

   graphics_info_t::seg_ids_in_atom_labels_flag = istat;
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

#if 0
   // 20101108 [Gatwick airport] Note: there is code in orient_view()
   // that presumes (with good reason) that this actually rotates the
   // view by nstep*stepsize*2
   //

   float spin_quat[4];
   graphics_info_t g;

   // spin it 1 degree
   float tbs =  g.get_trackball_size();
   for(int i=0; i<nsteps; i++) {
     trackball(spin_quat, 0, 0, 0.0174533*stepsize, 0.000, tbs);
     add_quats(spin_quat, g.quat, g.quat);
     graphics_draw();
   }
#endif
}

/* stepsize in degrees */
void rotate_x_scene(int nsteps, float stepsize) {

#if 0
  float spin_quat[4];
   graphics_info_t g;

   // spin it 1 degree
   float tbs =  g.get_trackball_size();
   for(int i=0; i<nsteps; i++) {
     trackball(spin_quat, 0, 0, 0.0, 0.0174533*stepsize, tbs);
     add_quats(spin_quat, g.quat, g.quat);
     graphics_draw();
   }
#endif
}

void rotate_z_scene(int nsteps, float stepsize) {

   // c.f globjects.cc:do_screen_z_rotate()
   //
#if 0

   float spin_quat[4];
   graphics_info_t g;
   for(int i=0; i<nsteps; i++) {
      trackball(spin_quat,
		1.0, 1.0,
		1.0, 1.0 + 0.0174533*stepsize,
		0.4);
      add_quats(spin_quat, g.quat, g.quat);
      graphics_draw();
   }
#endif
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

#if 0
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
#endif
}

void translate_scene_x(int nsteps) {
   std::cout << "placeholder" << std::endl;
}

void translate_scene_y(int nsteps) {
   std::cout << "placeholder" << std::endl;
}
void translate_scene_z(int nsteps) {
   std::cout << "placeholder" << std::endl;
}




/*  ----------------------------------------------------------------------- */
/*                  graphics background colour                              */
/*  ----------------------------------------------------------------------- */
/* stepsize in degrees */
void set_background_colour(double red, double green, double blue) {

   graphics_info_t g;
   g.background_colour[0] = red;
   g.background_colour[1] = green;
   g.background_colour[2] = blue;

   if (g.use_graphics_interface_flag)
      graphics_draw();
}

void
redraw_background() {
   graphics_draw();
}

int background_is_black_p() {

   graphics_info_t g;
   return g.background_is_black_p();
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

void remove_coordinates_glob_extension(const char *ext) {
   graphics_info_t g;
   g.remove_coordinates_glob_extension(std::string(ext));
}

void remove_data_glob_extension(const char *ext) {
   graphics_info_t g;
   g.remove_data_glob_extension(std::string(ext));
}

void remove_dictionary_glob_extension(const char *ext) {
   graphics_info_t g;
   g.remove_dictionary_glob_extension(std::string(ext));
}

void remove_map_glob_extension(const char *ext) {
   graphics_info_t g;
   g.remove_map_glob_extension(std::string(ext));
}



int do_anti_aliasing_state() {
   return graphics_info_t::do_anti_aliasing_flag;
}


void set_do_anti_aliasing(int state) {

   graphics_info_t g;
   g.set_do_anti_aliasing(state);
}


void set_do_GL_lighting(int state) {

   // no longer meaningful

}


int do_GL_lighting_state() {
   return graphics_info_t::do_lighting_flag;
}


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
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   std::cout << "FIXME in set_pick_cursor_index() " << std::endl;
#else
   graphics_info_t::pick_cursor_index = GdkCursorType(i);
#endif
}


/*  ------------------------------------------------------------------------ */
/*                       povray/raster3d interface                           */
/*  ------------------------------------------------------------------------ */
void raster3d(const char *filename) {

   graphics_info_t g;
   g.raster3d(std::string(filename));
}

void renderman(const char *filename) {
   graphics_info_t g;
   g.renderman(filename);
}

void povray(const char *filename) {

   graphics_info_t g;
   g.povray(std::string(filename));
}

void set_raster3d_bond_thickness(float f) {
   graphics_info_t::raster3d_bond_thickness = f;
}

void set_raster3d_bone_thickness(float f) {
   graphics_info_t::raster3d_bone_thickness = f;
}

void set_raster3d_atom_radius(float f) {

   graphics_info_t::raster3d_atom_radius = f;

}

void set_raster3d_density_thickness(float f) {

   graphics_info_t::raster3d_density_thickness = f;

}

void set_raster3d_water_sphere(int state) {

   graphics_info_t::raster3d_water_sphere_flag = state;

}

/*! \brief set the font size (as a string) for raster3d*/
void set_raster3d_font_size(const char *size_in) {

   graphics_info_t::raster3d_font_size = size_in;
}


void
raster_screen_shot() {  // run raster3d or povray and guile
                		         // script to render and display image

   // do some checking for povray/render here:

// BL says: lets make it for python too:
#ifdef USE_GUILE
   std::string cmd("(render-image)");  // this is a render function

   // cmd = "(povray-image)";

   safe_scheme_command(cmd);

#else
#ifdef USE_PYTHON
   std::string cmd("render_image()");  // this is a render function

   // cmd = "(povray-image)";

   safe_python_command(cmd);
#endif // USE_PYTHON
#endif //USE_GUILE
}

#ifdef USE_PYTHON
void
raster_screen_shot_py() {  // run raster3d or povray and guile
                		         // script to render and display image

   // do some checking for povray/render here:

   std::string cmd("render_image()");  // this is a render function

   // cmd = "(povray-image)";

   safe_python_command(cmd);
}
#endif // USE_PYTHON


void set_renderer_show_atoms(int istate) {

   graphics_info_t::renderer_show_atoms_flag = istate;
}

/*! \brief turn off shadows for raster3d output - give argument 0 to turn off  */
void set_raster3d_shadows_enabled(int state) {
   graphics_info_t::raster3d_enable_shadows = state;
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
	 s += " -remote 'openURL(";
	 s += u;
	 s += ",new-window)'";
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
      c = "open_url(";
      c += single_quote(u);
      c += ")";
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
	 for (unsigned int i=1; i<bits.size(); i++) {
	    s += "+";
	    s += bits[i];
	 }
	 // s += "+coot+site%3Awww.ysbl.york.ac.uk";
	 // s += "+coot+site%3Awww.biop.ox.ac.uk";
	 s += "+coot+site%3Awww2.mrc-lmb.cam.ac.uk";
	 browser_url(s.c_str());
      }
   }
}


/*  ----------------------------------------------------------------------- */
/*                  remote control                                          */
/*  ----------------------------------------------------------------------- */
/* section Remote Control */

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <fcntl.h>


/* tooltips */
void
set_tip_of_the_day_flag (int state) {
   graphics_info_t g;
   if (state == 0) {
	g.do_tip_of_the_day_flag = 0;
   } else {
	g.do_tip_of_the_day_flag = 1;
   }
}

/*  ----------------------------------------------------------------------- */
/*                  Surfaces                                                */
/*  ----------------------------------------------------------------------- */
void do_surface(int imol, int state) {

}

#include "c-interface-generic-objects.h"
#include "api/coot-molecule.hh"

/* per-chain functions can be added later */
void make_generic_surface(int imol, const char *selection_str, int mode) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      float col_scale = g.electrostatic_surface_charge_range;
      std::string selection_string(selection_str);
      // use api!
      coot::molecule_t cm(g.molecules[imol].atom_sel, 0, "");
      //  coot::simple_mesh_t smesh = cm.get_molecular_representation_mesh(selection_string, "ByOwnPotential", "Chains", 0);
      coot::simple_mesh_t smesh;
      std::string type = "";

      if (mode == 1) smesh = cm.get_molecular_representation_mesh(selection_string, "ByOwnPotential", "MolecularSurface", 0);
      if (mode == 2) smesh = cm.get_molecular_representation_mesh(selection_string, "Chains", "MolecularSurface", 0);

      if (mode == 1) type = "Electrostatic";
      if (mode == 2) type = "Molecular";

      g.attach_buffers();
      // do I need this vertex type conversion? Yes
      std::vector<s_generic_vertex> vertices(smesh.vertices.size());
      for (unsigned int i = 0; i < smesh.vertices.size(); i++) {
         vertices[i] = s_generic_vertex(smesh.vertices[i].pos,
                                        smesh.vertices[i].normal,
                                        smesh.vertices[i].color);
      }
      std::string object_name = type + " Surface " + std::to_string(imol) +
         std::string(" ") + std::string(selection_string);
      int obj_mesh = new_generic_object_number(object_name);
      meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
      obj.imol = imol;
      obj.mesh.name = object_name;
      obj.mesh.set_draw_mesh_state(true);
      obj.mesh.import(vertices, smesh.triangles);
      obj.mesh.set_material_specularity(0.7, 128);
      obj.mesh.setup_buffers();

      update_display_control_mesh_toggles(imol);
      graphics_draw();
   }
}
/* per-chain functions can be added later */
void make_molecular_surface(int imol, const char *selection_str) {

   make_generic_surface(imol, selection_str, 2);
}

/* per-chain functions can be added later */
void make_electrostatic_surface(int imol, const char *selection_str) {

   make_generic_surface(imol, selection_str, 1);
}

int molecule_is_drawn_as_surface_int(int imol) {

   // useless function now - remove it.
   return 1;
}


#ifdef USE_GUILE
void do_clipped_surface_scm(int imol, SCM residues_specs) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::residue_spec_t> res_specs_vec = scm_to_residue_specs(residues_specs);
      float col_scale = g.electrostatic_surface_charge_range;
      graphics_info_t::molecules[imol].make_surface(res_specs_vec, *g.Geom_p(), col_scale);
      graphics_draw();
   }
}
#endif //USE_GUILE

#ifdef USE_PYTHON
void do_clipped_surface_py(int imol, PyObject *residues_specs) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::residue_spec_t> res_specs_vec = py_to_residue_specs(residues_specs);
      float col_scale = g.electrostatic_surface_charge_range;
      graphics_info_t::molecules[imol].make_surface(res_specs_vec, *g.Geom_p(), col_scale);
      graphics_draw();
   }
}
#endif //USE_PYTHON

void set_electrostatic_surface_charge_range(float v) {
   graphics_info_t::electrostatic_surface_charge_range = v;
}

float get_electrostatic_surface_charge_range() {
   return graphics_info_t::electrostatic_surface_charge_range;
}

/*  ----------------------------------------------------------------------- */
/*           Sharpen                                                        */
/*  ----------------------------------------------------------------------- */
void sharpen(int imol, float b_factor) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].sharpen(b_factor, false, 0);
      graphics_draw();
   }
}

void sharpen_with_gompertz_scaling(int imol, float b_factor,
				   short int do_gompertz_scaling_flag,
				   float gompertz_factor) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].sharpen(b_factor, do_gompertz_scaling_flag,
					       gompertz_factor);
      graphics_draw();
   }
}




/*  ----------------------------------------------------------------------- */
/*                  Views                                                   */
/*  ----------------------------------------------------------------------- */
int add_view_here(const char *view_name) {

   std::cout << "------------------ debug: in add_view_here() with view name " << view_name << std::endl;

   std::string name(view_name);
   graphics_info_t g;
   coot::Cartesian rc = g.RotationCentre();
   float zoom = graphics_info_t::zoom;
   coot::view_info_t view(graphics_info_t::view_quaternion, rc, zoom, name);

   std::cout << "------------ in add_view_here() made a view with name: " << view.view_name << std::endl;
   std::cout << "------------ in add_view_here() made a view: " << view << std::endl;
   int n_views = graphics_info_t::views.size();
   graphics_info_t::views.push_back(view);
   int this_view_index = n_views;

   std::cout << "------------ in add_view_here() back is " << graphics_info_t::views.back() << std::endl;
   return this_view_index;
}

int add_view_raw(float rcx, float rcy, float rcz, float quat0, float quat1,
		  float quat2, float quat3, float zoom, const char *view_name) {

   glm::quat q(quat0, quat1, quat2, quat3);
   coot::Cartesian rc(rcx, rcy, rcz);
   coot::view_info_t v(q, rc, zoom, view_name);
   graphics_info_t::views.push_back(v);
   return (graphics_info_t::views.size() -1);
}


int remove_named_view(const char *view_name) {

   int r=0;
   bool found = 0;
   std::string vn(view_name);
   std::vector<coot::view_info_t> new_views;

   // This should be part of the view container class
   //
   unsigned int n_views = graphics_info_t::views.size();
   std::vector<coot::view_info_t>::iterator it; // needs to be const_iterator? depending on c++ version?
   for (it=graphics_info_t::views.begin(); it!=graphics_info_t::views.end(); it++) {
      if (it->view_name == vn) {
         graphics_info_t::views.erase(it);
         break;
      }
   }


   std::vector<std::string> command_strings;
   command_strings.push_back("remove_named_view");
   command_strings.push_back(single_quote(coot::util::intelligent_debackslash(view_name)));
   add_to_history(command_strings);
   return r;
}

/*! \brief the given view number */
void remove_view(int view_number) {

   int n_views = graphics_info_t::views.size();
   if (view_number >= 0) {
      if (view_number < n_views) {
         std::vector<coot::view_info_t>::iterator it; // needs to be const_iterator? depending on c++ version?
         int idx = 0;
         for (it=graphics_info_t::views.begin(); it!=graphics_info_t::views.end(); it++, idx++) {
            if (idx == view_number) {
               graphics_info_t::views.erase(it);
               break;
            }
         }
      }
   }

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

   for (unsigned int iv=0; iv<graphics_info_t::views.size(); iv++) {
      coot::view_info_t view1 = graphics_info_t::views[iv];
      // std::cout << "DEBUG:: View "<< iv << " " << view1.view_name << std::endl;
      if (! (view1.is_simple_spin_view_flag ||
	     view1.is_action_view_flag)) {
	 if ((iv+1) < graphics_info_t::views.size()) {
	    coot::view_info_t view2 = graphics_info_t::views[iv+1];
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
	 if ((iv+1) < graphics_info_t::views.size()) {
 	    std::cout << "DEBUG:: interpolating to  "<< iv+1 << " "
 		      << view1.view_name << std::endl;
	    coot::view_info_t view2 = graphics_info_t::views[iv+1];
	    if (!view2.is_simple_spin_view_flag && !view2.is_action_view_flag) {
	       // the quat was not set because this is a simple
	       // rotate, so we must generate it from the current
	       // position
	       coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
				  graphics_info_t::RotationCentre_y(),
				  graphics_info_t::RotationCentre_z());
	       coot::view_info_t current_view(graphics_info_t::view_quaternion,
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
   float zoom =  g.zoom;

   int r=0;
   bool found = false;
   coot::view_info_t v(graphics_info_t::view_quaternion, rc, zoom, "");

   std::vector<coot::view_info_t>::iterator it; // needs to be const_iterator? depending on c++ version?
   for (it=graphics_info_t::views.begin(); it!=graphics_info_t::views.end(); it++) {
      if (it->matches_view(v)) {
         graphics_info_t::views.erase(it);
         break;
      }
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

   std::cout << "---------------- clear_all_views() " << std::endl;
   graphics_info_t::views.clear();
}

int go_to_view_number(int view_number, int snap_to_view_flag) {

   int r = 0;
   graphics_info_t g;
   if ((int(graphics_info_t::views.size()) > view_number) && (view_number >= 0)) {
      coot::view_info_t view = graphics_info_t::views[view_number];
      if (view.is_simple_spin_view_flag) {
	 int nsteps = 2000;
         nsteps = 500;
	 if (graphics_info_t::views_play_speed > 0.000000001)
	    nsteps = int(static_cast<float>(nsteps)/graphics_info_t::views_play_speed);
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
               g.view_quaternion = view.quaternion;
	    } else {
	       coot::view_info_t this_view(g.view_quaternion, g.RotationCentre(), g.zoom, "");
	       int nsteps = 2000;
	       if (graphics_info_t::views_play_speed > 0.000000001)
		  nsteps = int(2000.0/graphics_info_t::views_play_speed);
	       coot::view_info_t::interpolate(this_view,
					      graphics_info_t::views[view_number], nsteps);
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

   if (true) {
      std::cout << "debug in n_views(): with n_views " <<  graphics_info_t::views.size() << std::endl;
      unsigned int nv =  graphics_info_t::views.size();
      for (unsigned int i=0; i<nv; i++) {
         std::cout << "debug:: n_views() " << i << " has name " << graphics_info_t::views.at(i).view_name
                   << " " << graphics_info_t::views.at(i) << std::endl;
      }
   }

   add_to_history_simple("n-views");
   return graphics_info_t::views.size();
}

/*! \brief return the name of the given view, if view_number does not
  specify a view return #f */
#ifdef USE_GUILE
SCM view_name(int view_number) {

   SCM r = SCM_BOOL_F;
   int n_view = graphics_info_t::views.size();
   if (view_number < n_view)
      if (view_number >= 0) {
	 std::string name = graphics_info_t::views[view_number].view_name;
	 r = scm_from_locale_string(name.c_str());
      }
   return r;
}
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *view_name_py(int view_number) {

   PyObject *r;
   r = Py_False;
   int n_view = graphics_info_t::views.size();
   if (view_number < n_view)
      if (view_number >= 0) {
         std::string name = graphics_info_t::views[view_number].view_name;
         r = myPyString_FromString(name.c_str());
      }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // PYTHON

#ifdef USE_GUILE
SCM view_description(int view_number) {

   SCM r = SCM_BOOL_F;
   if (view_number >= 0)
      if (view_number < int(graphics_info_t::views.size())) {
	 std::string d = graphics_info_t::views[view_number].description;
	 if (d != "") {
	    r = scm_from_locale_string(d.c_str());
	 }
      }
   return r;
}
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *view_description_py(int view_number) {

   PyObject *r;
   r = Py_False;
   if (view_number >= 0)
      if (view_number < int(graphics_info_t::views.size())) {
         std::string d = graphics_info_t::views[view_number].description;
         if (d != "") {
            r = myPyString_FromString(d.c_str());
         }
      }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // PYTHON

#include "view.hh"

/*! \brief save views to view_file_name */
void save_views(const char *view_file_name) {

   unsigned int n_views = graphics_info_t::views.size();
   if (n_views > 0) {
      std::ofstream f(view_file_name);
      if (! f) {
	 std::cout << "Cannot open view output file"
		   << view_file_name << std::endl;
      } else {
#ifdef USE_GUILE
	 f << "; Views\n";
#else
#ifdef USE_PYTHON
	 f << "# Views\n";
#endif // PYTHON
#endif // GUILE
	 for (unsigned int i=0; i<n_views; i++) {
	    f << graphics_info_t::views[i];
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
   coot::view_info_t view(name, func);  // an action view
   graphics_info_t::views.push_back(view);
   return (graphics_info_t::views.size() -1);
}

// return the view number of the new view
int insert_action_view_after_view(int view_number, const char *view_name, const char *action_function) {

   // FIX this

   int r = -1;
   std::string name(view_name);
   std::string func(action_function);
   coot::view_info_t view(name, func);  // an action view
   int n_views = graphics_info_t::views.size();
   if (view_number >= n_views) {
      graphics_info_t::views.push_back(view);
      r = (graphics_info_t::views.size() -1);
   } else {
      // insert a view
      std::vector <coot::view_info_t> new_views;
      for (int iview=0; iview<n_views; iview++) {
	 new_views.push_back(graphics_info_t::views[iview]);
	 if (iview == view_number)
	    new_views.push_back(view);
      }
      graphics_info_t::views = new_views; // bleuch.
      r = view_number + 1;
   }
   return r;
}

void add_view_description(int view_number, const char *descr) {

   if (view_number < int(graphics_info_t::views.size()))
      if (view_number >= 0)
	 graphics_info_t::views[view_number].add_description(descr);
}

#ifdef USE_GUILE
void go_to_view(SCM view) {

   std::cout << "FIXME go_to_view()" << std::endl;

}
#endif // USE_GUILE

#ifdef USE_PYTHON
void go_to_view_py(PyObject *view) {

   int len_view;

   len_view = PyObject_Length(view);

   if (len_view == 4) {

      PyObject *quat_python;
      graphics_info_t g;
      int nsteps = 2000;
      if (graphics_info_t::views_play_speed > 0.000000001)
         nsteps = int(2000.0/graphics_info_t::views_play_speed);

      // What is the current view:
      //
      std::string name("Current Position");
      coot::Cartesian rc = g.RotationCentre();
      float zoom = graphics_info_t::zoom;
      coot::view_info_t view_c(graphics_info_t::view_quaternion, rc, zoom, name);

      // view_target is where we want to go
      quat_python = PyList_GetItem(view, 0);
      int len_quat = PyObject_Length(quat_python);
      if (len_quat == 4) {
         PyObject *q0_python = PyList_GetItem(quat_python, 0);
         PyObject *q1_python = PyList_GetItem(quat_python, 1);
         PyObject *q2_python = PyList_GetItem(quat_python, 2);
         PyObject *q3_python = PyList_GetItem(quat_python, 3);
         glm::quat quat_target(PyFloat_AsDouble(q0_python),
                               PyFloat_AsDouble(q1_python),
                               PyFloat_AsDouble(q2_python),
                               PyFloat_AsDouble(q3_python));

         PyObject *rc_target_python = PyList_GetItem(view, 1);
         int len_rc_target = PyObject_Length(rc_target_python);
         if (len_rc_target == 3) {

            PyObject *centre_x = PyList_GetItem(rc_target_python, 0);
            PyObject *centre_y = PyList_GetItem(rc_target_python, 1);
            PyObject *centre_z = PyList_GetItem(rc_target_python, 2);

            double x = PyFloat_AsDouble(centre_x);
            double y = PyFloat_AsDouble(centre_y);
            double z = PyFloat_AsDouble(centre_z);
            coot::Cartesian rc_target(x,y,z);

            PyObject *target_zoom_python = PyList_GetItem(view, 2);
            double zoom_target = PyFloat_AsDouble(target_zoom_python);

            PyObject *name_target_python = PyList_GetItem(view, 3);
            std::string name_target = PyBytes_AS_STRING(PyUnicode_AsUTF8String(name_target_python));

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
#endif // PYTHON


int add_spin_view(const char *view_name, int n_steps, float degrees_total) {

   coot::view_info_t v(view_name, n_steps, degrees_total);
   graphics_info_t::views.push_back(v);

   std::string cmd = "add-spin-view";
   std::vector<coot::command_arg_t> args;
   args.push_back(view_name);
   args.push_back(n_steps);
   args.push_back(degrees_total);
   add_to_history_typed(cmd, args);
   return (graphics_info_t::views.size() -1);
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

   std::cout << " =============== setting mutex lock (scheme version) =========" << std::endl;
   //
   // (This mutex lock *and* waiting flag may be overly complex now
   // that we simply use g_idle_add())
   graphics_info_t::socket_string_waiting_mutex_lock = 1;
   graphics_info_t::socket_string_waiting = s;
   graphics_info_t::have_socket_string_waiting_flag = 1;

   std::cout << "DEBUG:: set_socket_string_waiting() socket_string_waiting set to \""
	     << graphics_info_t::socket_string_waiting << "\"" << std::endl;

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

/*! \brief feed the main thread a python script to evaluate */
void set_socket_python_string_waiting(const char *s) {

   graphics_info_t::socket_python_string_waiting = s;
   graphics_info_t::have_socket_python_string_waiting_flag = 1;

   GSourceFunc f = graphics_info_t::process_socket_python_string_waiting_bool;
   g_idle_add(f, NULL); // if f returns FALSE then f is not called again.
}



// should be in c-interface-maps?

void set_map_sharpening_scale_limit(float f) {
   graphics_info_t::map_sharpening_scale_limit = f;
}


// should be in c-interface-maps?
//
/*! \brief sets the density map of the given molecule to be drawn as a
  (transparent) solid surface. */
void set_draw_solid_density_surface(int imol, short int state) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].set_draw_solid_density_surface(state);
      graphics_draw();
   }
}

void
set_solid_density_surface_opacity(int imol, float opacity) {

   if (is_valid_map_molecule(imol)) {
      graphics_info_t::molecules[imol].density_surface_opacity = opacity;
      graphics_draw();
   }
}

float get_solid_density_surface_opacity(int imol) {

   float opacity = -1;
   if (is_valid_map_molecule(imol)) {
      opacity = graphics_info_t::molecules[imol].density_surface_opacity;
   }
   return opacity;
}


void
set_flat_shading_for_solid_density_surface(short int state) {
   graphics_info_t::do_flat_shading_for_solid_density_surface = state;
   graphics_draw();
}

/*! \brief simple on/off screendoor transparency at the moment, an
  opacity > 0.0 will turn on screendoor transparency (stippling). */
void set_transparent_electrostatic_surface(int imol, float opacity) {

   if (is_valid_model_molecule(imol)) {
      bool flag = 0;
      if (opacity > 0.0)
	 if (opacity < 0.9999)
	    flag = 1;
      graphics_info_t::molecules[imol].transparent_molecular_surface_flag = flag;
      graphics_draw();
   }

}

/*! \brief return 1.0 for non transparent and 0.5 if screendoor
  transparency has been turned on. */
float get_electrostatic_surface_opacity(int imol) {

   float r = -1;
   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].transparent_molecular_surface_flag == 1)
	 r = 0.5;
      else
	 r = 1.0;
   }
   return r;
}


/*! \brief load tutorial model and data  */
void load_tutorial_model_and_data() {

   // implement this

   /*
	   (let* ((prefix-dir (getenv "COOT_PREFIX")))

	     (let* ((pkg-data-dir
		     (if (string? prefix-dir)
			 (append-dir-dir (append-dir-dir prefix-dir "share") "coot")
			 (pkgdatadir)))
		    (data-dir (append-dir-dir pkg-data-dir "data"))
		    (pdb-file-name (append-dir-file data-dir "tutorial-modern.pdb"))
		    (mtz-file-name (append-dir-file data-dir "rnasa-1.8-all_refmac1.mtz")))

	       (read-pdb pdb-file-name)
	       (make-and-draw-map mtz-file-name "FWT" "PHWT" "" 0 0)
	       (make-and-draw-map mtz-file-name "DELFWT" "PHDELWT" "" 0 1)))))

   */

   std::string p = coot::package_data_dir();
   std::string d = coot::util::append_dir_dir(p, "data");

   std::string pdb_fn = coot::util::append_dir_file(d, "tutorial-modern.pdb");
   std::string mtz_fn = coot::util::append_dir_file(d, "rnasa-1.8-all_refmac1.mtz");

   std::cout << "--------- load_tutorial_model_and_data() " << pdb_fn << std::endl;
   std::cout << "--------- load_tutorial_model_and_data() " << mtz_fn << std::endl;

   int imol = handle_read_draw_molecule_with_recentre(pdb_fn.c_str(), true);

   int imol_map = make_and_draw_map_with_refmac_params(mtz_fn.c_str(), "FWT", "PHWT", "", 0, 0, 1, "FGMP18", "SIGFGMP18", "FreeR_flag", 1);
   int imol_diff_map = make_and_draw_map(mtz_fn.c_str(), "DELFWT", "PHDELWT", "", 0, 1);

   if (false) {
      std::cout << "--------- imol: " << imol << std::endl;
      std::cout << "--------- imol_map: " << imol_map << std::endl;
      std::cout << "--------- imol_diff_map: " << imol_diff_map << std::endl;
   }

   // 2025-03-26-PE not all GLibs have g_idle_add_once() at the moment
   // gint idle = g_idle_add_once((GSourceOnceFunc)g.graphics_grab_focus, NULL);
   auto callback = +[] (gpointer data) {
     graphics_info_t g;
     g.graphics_grab_focus();
     return gboolean(FALSE);
   };
   gint idle = g_idle_add(callback, NULL);

}

// is the probe executable available?
// 1 for yes, 0 for no.
//
int probe_available_p() {
   int r = graphics_info_t::probe_available;
   return r;
}

#ifdef USE_PYTHON
// is the probe executable available?
// 1 for yes, 0 for no.
//
int probe_available_p_py() {
   int r = graphics_info_t::probe_available;
   return r;
}
#endif // USE_PYTHON


/*! \brief set the state of showing chiral volume outlier markers - of a model molecule that is,
   not the intermediate atoms (derived from restraints) */
void set_show_chiral_volume_outliers(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].draw_chiral_volume_outlier_markers_flag = state;
      graphics_info_t::molecules[imol].fill_chiral_volume_outlier_marker_positions(state);
      graphics_info_t::update_chiral_volume_outlier_marker_positions();
      graphics_draw();
   }

}

/*! \brief set the state of showing non-bonded contact markers - of a model molecule that is,
   not the intermediate atoms (derived from restraints) */
void set_show_non_bonded_contact_baddies_markers(int imol, int state) {

   if (is_valid_model_molecule(imol)) {
      // for the molecule, not the intermediate atoms.
      graphics_info_t::molecules[imol].set_show_non_bonded_contact_baddies_markers(state);
      graphics_draw();
   }
}


#ifdef USE_PYTHON
/* Get model molecule list */
PyObject *get_model_molecule_list_py() {

   std::vector<int> v;
   graphics_info_t g;
   unsigned int n = g.n_molecules();
   for (unsigned int i=0; i<n; i++) {
      if (is_valid_model_molecule(i))
         v.push_back(i);
   }
   PyObject *l_py = PyList_New(v.size());
   for (unsigned int ii=0; ii<v.size(); ii++) {
      PyList_SetItem(l_py, ii, PyLong_FromLong(v[ii]));
   }
   return l_py;
}
#endif

#ifdef USE_GUILE
/* Get model molecule list */
SCM get_model_molecule_list_scm() {

   std::vector<int> v;
   graphics_info_t g;
   unsigned int n = g.n_molecules();
   for (unsigned int i=0; i<n; i++) {
      if (is_valid_model_molecule(i))
         v.push_back(i);
   }
   SCM l_scm = SCM_EOL;
   // 20240802-PE fill me.
   return l_scm;
}
#endif
