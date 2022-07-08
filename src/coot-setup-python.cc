/* src/coot-setup-python.cc
 * 
 * Copyright 2008 by the University of Oxford
 * Copyright 2014, 2016 by Medical Research Council
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

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <cstring>

#include "compat/coot-sysdep.h"
#include <gtk/gtk.h>
#include "coot-setup-python.hh"
#include "python-classes.hh"
#include "coot-glue.hh"
#include "c-interface.h"
#include "c-interface-preferences.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
// #include "c-inner-main.h" // for does_file_exist() 2015025-PE: not needed for modern
                             // guile interface

#include "graphics-info.h"
#include "command-line.hh"

#include <string>
#include <iostream>

#include <sys/stat.h>
#include <glob.h>


void setup_python(int argc, char **argv) {


#ifdef USE_PYTHON
     //  (on Mac OS, call PyMac_Initialize() instead)
     //http://www.python.org/doc/current/ext/embedding.html

#ifdef USE_PYMAC_INIT 
  PyMac_Initialize();
#else  
  Py_Initialize(); // otherwise it core dumps saying python
  // interpreter not initialized (or something).
  PySys_SetArgv(argc, argv);
#endif     

  init_coot(); // i.e. SWIG_init for python, best we do this before
               // running .coot.py, eh?

  
  std::string home_directory = coot::get_home_dir();;

     short int use_graphics_flag = use_graphics_interface_state();

     // First load coot.py, then load the standard startup files, 
     // then load 0-coot.state.py

     std::string pydirectory = PKGPYTHONDIR; /* prefix/lib/python2.7/site-packages/coot */
     
     char *pydirectory_over = getenv("COOT_PYTHON_DIR");
     if (pydirectory_over)
	pydirectory = pydirectory_over;

     int err = import_python_module("coot", 0);
     if (err == -1) {
	std::cout << "ERROR:: could not import coot.py" << std::endl;
     } else {
	std::cout << "INFO:: coot.py imported" << std::endl;

	if (! use_graphics_flag) {
	   safe_python_command("global use_gui_qm; use_gui_qm = False");
	} else { 
	   // we have gui
	   // BL says:: lets initialize glue too but only if we have pygtk 
	   // (and gtk2)

#ifdef USE_PYGTK
	   initcoot_python();
	   std::cout << "INFO:: coot_python initialized" << std::endl;
#   ifdef USE_GUILE_GTK
	   safe_python_command("global use_gui_qm; use_gui_qm = 2");
#   else
	   safe_python_command("global use_gui_qm; use_gui_qm = 1");
#   endif
#else
	   safe_python_command("global use_gui_qm; use_gui_qm = False");
#endif // PYTGK
	}
	
	std::string coot_load_modules_dot_py = "coot_load_modules.py";
	std::string coot_py_file_name = graphics_info_t::add_dir_file(pydirectory, coot_load_modules_dot_py);
	if (coot::file_exists(coot_py_file_name)) { 
	   run_python_script(coot_py_file_name.c_str());
	} else {
	   std::cout << "WARNING:: No coot modules found! Python scripting crippled. " 
		     << std::endl;
	}
     }
     // try to load extra dir files (if exist) do before preferences (as for
     // scheme version
     try_load_python_extras_dir();

     // we only want to run one state file if using both scripting
     // languages.  Let that be the guile one.
     //
     try_load_dot_coot_py_and_python_scripts(home_directory);
     
#ifndef USE_GUILE

     // we get here if there is no guile but there is python and the ifdef is only here so 
     // that we don't run both state script.

     command_line_data cld = parse_command_line(argc, argv);
     handle_command_line_data(cld);

     // BL says::should this still be here?
     run_state_file_maybe(); // run local 0-coot.state.py?

     // run_update_self_maybe(); nope.  Not yet.
     
#endif // USE_GUILE - not both start-up scripts

#endif // USE_PYTHON  

}

void
setup_python_classes() {
#ifdef USE_PYTHON

      init_pathology_data();

#endif

}

void try_load_dot_coot_py_and_python_scripts(const std::string &home_directory) {

   if (graphics_info_t::run_startup_scripts_flag) {

      short int use_graphics_flag = use_graphics_interface_state();

      // load preferences file .coot_preferences.py
      std::string preferences_dir = graphics_info_t::add_dir_file(home_directory, ".coot-preferences");
      struct stat buff;
      int preferences_dir_status = stat(preferences_dir.c_str(), &buff);

      if (preferences_dir_status != 0) {
         std::cout << "INFO:: preferences directory " << preferences_dir
                   << " does not exist. Won't read preferences." << std::endl;;
      } else {
	 // load all .py files
	 glob_t myglob;
	 int flags = 0;
	 //std::string glob_patt = "/*.py";
	 std::string glob_file = preferences_dir;
	 glob_file += "/*.py";
	 glob(glob_file.c_str(), flags, 0, &myglob);
	 // dont load the coot_toolbuttons.py if no graphics
	 // same for key_bindings (and potentially others)
	 std::set<std::string> exclude_py_files;
         if (! use_graphics_flag) {
            exclude_py_files.insert("coot_toolbuttons.py");
            exclude_py_files.insert("template_key_bindings.py");
         }

         // make this split so that we can run curlew scripts after, and xenops scripts after that.
         //
         std::vector<std::string> basic_scripts;
         std::vector<std::string> xenops_scripts;
         std::vector<std::string> curlew_scripts;
         std::string coot_preferences_py_script;

         for (char **p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) {
            std::string preferences_script(*p);
            if (exclude_py_files.find(preferences_script) == exclude_py_files.end()) {
               bool done = false;
               if (preferences_script.length() > 6) {
                  if (preferences_script.substr(0,6) == "xenops") {
                     done = true;
                     xenops_scripts.push_back(preferences_script);
                  }
               }
               if (preferences_script.length() > 6) {
                  if (preferences_script.substr(0,6) == "curlew") {
                     done = true;
                     curlew_scripts.push_back(preferences_script);
                  }
               }
               if (! done)
                  basic_scripts.push_back(preferences_script);
               if (preferences_script == "coot_preferences.py")
                  coot_preferences_py_script = preferences_script;
            }
         }
         globfree(&myglob);

         for(const auto &script_fn : basic_scripts)
            run_python_script(script_fn.c_str()); // bleugh
         if (! coot_preferences_py_script.empty())
            run_python_script(coot_preferences_py_script.c_str());
         for(const auto &script_fn : curlew_scripts)
            run_python_script(script_fn.c_str());
         for(const auto &script_fn : xenops_scripts)
            run_python_script(script_fn.c_str());
      }

      // update the preferences
      make_preferences_internal();

      // load personal coot file .coot.py
      std::string filename = ".coot.py";
      if (! home_directory.empty()) {
	 std::string coot_py_file_name = graphics_info_t::add_dir_file(home_directory, filename);
	 if (coot::file_exists(coot_py_file_name)) {
	    std::cout << "Loading ~/.coot.py..." << std::endl;
	    run_python_script(coot_py_file_name.c_str());
	 }
      }
   }
}
