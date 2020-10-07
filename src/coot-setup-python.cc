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
#ifdef USE_PYMAC_INIT 

  //  (on Mac OS, call PyMac_Initialize() instead)
  //http://www.python.org/doc/current/ext/embedding.html
  //
  PyMac_Initialize();

#else

   wchar_t** _argv = static_cast<wchar_t **>(PyMem_Malloc(sizeof(wchar_t*)*argc));
   for (int i=0; i<argc; i++) {
      wchar_t* arg = Py_DecodeLocale(argv[i], NULL);
      _argv[i] = arg;
   }
   Py_InitializeEx(0);
   PySys_SetArgv(argc, _argv);

#endif // USE_PYMAC_INIT


   std::string pydirectory = PKGPYTHONDIR; /* prefix/lib/python2.7/site-packages/coot */

   std::cout << "debug:: in setup_python() pydirectory is " << pydirectory << std::endl;

   int err = import_python_module("coot", 0);
   if (err == -1) {
      std::cout << "ERROR:: could not import coot.py" << std::endl;
   } else {
      std::cout << "INFO:: coot.py imported" << std::endl;
      std::string coot_load_modules_dot_py = "coot_load_modules.py";
      std::string coot_py_file_name = coot::util::append_dir_file(pydirectory, coot_load_modules_dot_py);
      if (coot::file_exists(coot_py_file_name)) {
         PyRun_SimpleString("global use_gui_qm; use_gui_qm = False");
         FILE *fp = fopen(coot_py_file_name.c_str(), "r");
         PyRun_SimpleFile(fp, coot_py_file_name.c_str());
         fclose(fp);
         std::cout << "read " << coot_py_file_name << std::endl;
      } else {
         std::cout << "WARNING:: No coot modules found! Python scripting crippled. "
                   << std::endl;
      }
   }


#if 1 // Add a test for PyGObject (the new Python GTK interface)

   initcoot_python_gobject(); // this is not a good name for this function. We need to say
                      // this this is the module that wraps the glue to get
                      // the status-bar, menu-bar etc.
                      // Done.

#endif // GObject
#endif // USE_PYTHON

}

void
setup_python_classes() {
#ifdef USE_PYTHON

  init_pathology_data(); 

#endif

} 

void try_load_dot_coot_py_and_preferences(const std::string &home_directory) {

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
	 size_t count;
	 // dont load the coot_toolbuttons.py if no graphics
	 // same for key_bindings (and potentially others)
	 std::vector<std::string> exclude_py_files;
	 exclude_py_files.push_back("coot_toolbuttons.py");
	 exclude_py_files.push_back("template_key_bindings.py");
	 std::size_t found_substr;

	 char **p;
	 for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) { 
	    char *preferences_file(*p);
	    found_substr = std::string::npos;
	    for (unsigned int i=0; i<exclude_py_files.size(); i++) {
	       found_substr = ((std::string)preferences_file).find(exclude_py_files[i]);
	       if (found_substr != std::string::npos)
		  break;
	    }
	    char *found2 = strstr(preferences_file, "coot_preferences.py");
	    if (((found_substr == std::string::npos) || (use_graphics_flag)) && 
		((!found2) || (prefer_python()))) {
	       if (false) // too verbose
		  std::cout << "INFO:: loading preferences file " << preferences_file
			    << std::endl;
	       run_python_script(preferences_file); // run_python_script writes the file-name to
	                                            // the terminal
	    }
	 }
	 globfree(&myglob);
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
