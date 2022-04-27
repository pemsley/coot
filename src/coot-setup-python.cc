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

#include <algorithm>
#include <cctype>
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
#include <vector>

#include <sys/stat.h>

template <typename T, typename U>
static
bool contains(const T &what, const U &where) {
    return std::find(where.cbegin(), where.cend(), what) != where.cend();
}

static
const std::vector<std::string> EXLCUDED_PY_FILES_NOGRAPHICS = {
    "coot_toolbuttons.py",
    "template_key_bindings.py"
};

class PreferenceFiles {
public:
    std::vector<std::string> basic_scripts;
    std::vector<std::string> xenops_scripts;
    std::vector<std::string> curlew_scripts;
    std::string coot_preferences_py_script;
};

static
PreferenceFiles gather_preferences_files(const std::string &preferences_dir, bool use_graphics) {
    const auto &excluded_py_files = use_graphics ? std::vector<std::string>{} : EXLCUDED_PY_FILES_NOGRAPHICS;
    auto candidates = coot::gather_files_by_patterns(preferences_dir, { "*.py" });

    PreferenceFiles prefs{};
    for (auto &file : candidates) {
        if (contains(file, excluded_py_files))
            continue;

        if (file.length() > 6) {
            auto prefix = file.substr(0, 6);
            std::transform(prefix.begin(), prefix.end(), prefix.begin(), [](unsigned char ch) { return std::tolower(ch); });

            if (prefix == "xenops") {
                prefs.xenops_scripts.push_back(std::move(file));
                continue;
            } else if (prefix == "curlew") {
                prefs.curlew_scripts.push_back(std::move(file));
                continue;
            }
        }

        prefs.basic_scripts.push_back(std::move(file));

        if (file == "coot_preferences.py")
            prefs.coot_preferences_py_script = std::move(file);
    }

    return prefs;
}

void setup_python_basic(int argc, char **argv) {

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

   auto get_pythondir = [] () {
                           std::string p = coot::prefix_dir();
                           std::string dp   = coot::util::append_dir_dir(p,   "lib");
                           std::string python_version = "python";
                           python_version += coot::util::int_to_string(PY_MAJOR_VERSION);
                           python_version += ".";
                           python_version += coot::util::int_to_string(PY_MINOR_VERSION);
                           std::string ddp  = coot::util::append_dir_dir(dp,  python_version);
                           std::string dddp = coot::util::append_dir_dir(ddp, "site-packages");
                           return dddp;
                        };
   auto get_pkgpythondir = [get_pythondir] () {
                              std::string d = get_pythondir();
                              std::string dp   = coot::util::append_dir_dir(d, "coot");
                              return dp;
                           };

   // std::string pkgpydirectory = PKGPYTHONDIR;
   // std::string pydirectory = PYTHONDIR;
   // use ${prefix}/lib/python3.9/site-package for PYTHONDIR
   // use ${pythondir}/coot' for PKGPYTHONDIR (i.e. PYTHONDIR + "/coot")

   std::string pkgpydirectory = get_pkgpythondir();
   std::string    pydirectory = get_pythondir();

   if (false) {
      std::cout << "debug:: in setup_python()    pydirectory is " << pydirectory << std::endl;
      std::cout << "debug:: in setup_python() pkgpydirectory is " << pkgpydirectory << std::endl;
   }

   PyObject *sys_path = PySys_GetObject("path");
   PyList_Append(sys_path, PyUnicode_FromString(pydirectory.c_str()));

   // int err = PyRun_SimpleString("import coot");

#endif // USE_PYTHON

}

void setup_python_coot_module() {

  PyObject *coot = PyImport_ImportModule("coot");
   if (! coot) {
      std::cout << "ERROR:: setup_python_coot_module() Null coot" << std::endl;
   }
}

void setup_python_with_coot_modules(int argc, char **argv) {
#ifdef USE_PYTHON

   auto get_pythondir = [] () {
                           std::string p = coot::prefix_dir();
                           std::string dp   = coot::util::append_dir_dir(p,   "lib");
                           std::string python_version = "python";
                           python_version += coot::util::int_to_string(PY_MAJOR_VERSION);
                           python_version += ".";
                           python_version += coot::util::int_to_string(PY_MINOR_VERSION);
                           std::string ddp  = coot::util::append_dir_dir(dp,  python_version);
                           std::string dddp = coot::util::append_dir_dir(ddp, "site-packages");
                           return dddp;
                        };
   auto get_pkgpythondir = [get_pythondir] () {
                              std::string d = get_pythondir();
                              std::string dp   = coot::util::append_dir_dir(d, "coot");
                              return dp;
                           };

   // std::string pkgpydirectory = PKGPYTHONDIR;
   // std::string pydirectory = PYTHONDIR;
   // use ${prefix}/lib/python3.9/site-package for PYTHONDIR
   // use ${pythondir}/coot' for PKGPYTHONDIR (i.e. PYTHONDIR + "/coot")

   std::string pkgpydirectory = get_pkgpythondir();
   std::string pydirectory = get_pythondir();

   if (false) {
      std::cout << "debug:: in setup_python()    pydirectory is " << pydirectory << std::endl;
      std::cout << "debug:: in setup_python() pkgpydirectory is " << pkgpydirectory << std::endl;
   }

   PyObject *sys_path = PySys_GetObject("path");
   PyList_Append(sys_path, PyUnicode_FromString(pydirectory.c_str()));

   // int err = PyRun_SimpleString("import coot");

   PyObject *sys = PyImport_ImportModule("sys");
   if (! sys) {
      std::cout << "ERROR:: setup_python() Null sys" << std::endl;
   } else {
      // std::cout << "sys imported" << std::endl;
   }
   PyObject *coot = PyImport_ImportModule("coot");

   if (! coot) {
      std::cout << "ERROR:: setup_python() Null coot" << std::endl;
   } else {

      if (true) {
         initcoot_python_gobject(); // this is not a good name for this function. We need to say
                                    // this this is the module that wraps the glue to get
                                    // the status-bar, menu-bar etc. i.e. coot_python_api
         PyObject *io = PyImport_ImportModule("coot_utils"); // this imports coot_gui (which seems wrong)

         std::cout << "@@@@@@@@@@@@@@ coot_utils was imported " << io << std::endl;
         
         
         // date  This has do be done carefully - bit by bit. extension.py has many Python2/Python3
         // idioms.
         // PyImport_ImportModule("extensions");
      }
   }

   PyErr_PrintEx(0);

   // Read the startup-scripts elsewhere. The function is just to
   // read the standard coot python *modules*.

   // std::string home_directory = coot::get_home_dir();
   // try_load_dot_coot_py_and_python_scripts(home_directory);

#endif // USE_PYTHON

}

void
setup_python_classes() {
#ifdef USE_PYTHON

  init_pathology_data();

#endif

}

void try_load_dot_coot_py_and_python_scripts(const std::string &home_directory) {

   std::cout << "--------------- try_load_dot_coot_py_and_python_scripts from "
             << home_directory << " " << graphics_info_t::run_startup_scripts_flag << std::endl;

   if (graphics_info_t::run_startup_scripts_flag) {

      // load preferences file .coot_preferences.py
      // std::string preferences_dir = graphics_info_t::add_dir_file(home_directory, ".coot-preferences");
      //
      // now coot startup scripts will be read from and written to the ~/.coot directory
      // (if it is not a file already)
      //
      std::string startup_scripts_dir = graphics_info_t::add_dir_file(home_directory, ".coot");
      struct stat buff;
      int preferences_dir_status = stat(startup_scripts_dir.c_str(), &buff);

      if (preferences_dir_status != 0) {
         std::cout << "INFO:: preferences directory " << startup_scripts_dir
                   << " does not exist. Won't read preferences." << std::endl;;
      } else {
        const auto prefs = gather_preferences_files(startup_scripts_dir, use_graphics_interface_state());

         for(const auto &script_fn : prefs.basic_scripts)
            run_python_script(script_fn.c_str()); // bleugh
         if (!prefs.coot_preferences_py_script.empty())
            run_python_script(prefs.coot_preferences_py_script.c_str());
         for(const auto &script_fn : prefs.curlew_scripts)
            run_python_script(script_fn.c_str());
         for(const auto &script_fn : prefs.xenops_scripts)
            run_python_script(script_fn.c_str());
      }

      // update the preferences
      make_preferences_internal();

#if 0
      // load personal coot file .coot.py
      std::string filename = ".coot.py";
      if (! home_directory.empty()) {
	 std::string coot_py_file_name = graphics_info_t::add_dir_file(home_directory, filename);
	 if (coot::file_exists(coot_py_file_name)) {
	    std::cout << "Loading ~/.coot.py..." << std::endl;
	    run_python_script(coot_py_file_name.c_str());
	 }
      }
#endif
   }
}
