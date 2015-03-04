#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"
#include <gtk/gtk.h>
#include "coot-setup-python.hh"
#include "coot-glue.hh"
#include "c-interface.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "c-inner-main.h" // for does_file_exist()

#include "graphics-info.h"
#include "command-line.hh"

#include <string>
#include <iostream>

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

  char *home_directory = getenv("HOME");

  // I won't mess with this - not sure what it does.
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
     home_directory = getenv("COOT_HOME");
     std::string pkgdirectory = PKGDATADIR;
     if (!home_directory) {
       home_directory = getenv("HOME");
     }
     if (!home_directory) {
       home_directory = (char *)pkgdirectory.c_str();
     }
#endif

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
	char *coot_load_modules_dot_py_checked = 
	   does_file_exist(pydirectory.c_str(), coot_load_modules_dot_py.c_str());
	if (coot_load_modules_dot_py_checked) { 
	   run_python_script(coot_load_modules_dot_py_checked);
	} else {
	   std::cout << "WARNING:: No coot modules found! Python scripting crippled. " 
		     << std::endl;
	}
     }
     // try to load extra dir files (if exist) do before preferences (as for
     // scheme version
     try_load_python_extras_dir();

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
       
       char **p;
       for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) { 
         char *preferences_file(*p);
         // dont load the coot_toolbuttons.py if no graphics
         char *found = strstr(preferences_file, "coot_toolbuttons.py");
         char *found2 = strstr(preferences_file, "coot_preferences.py");
         if (((!found) || (use_graphics_flag)) && 
             ((!found2) || (prefer_python()))) {
            std::cout << "INFO:: loading preferences file " << preferences_file
                      << std::endl;
            run_python_script(preferences_file);
         }
       }
       globfree(&myglob);
     }
     // update the preferences
     make_preferences_internal();

     // load personal coot file .coot.py
     const char *filename = ".coot.py";
     if (home_directory) { 
        char *check_file = does_file_exist(home_directory, filename);
        if (check_file) {
	   std::cout << "Loading ~/.coot.py..." << std::endl;
	   run_python_script(check_file);
        }
     }

     // we only want to run one state file if using both scripting
     // languages.  Let that be the guile one.
     
#ifndef USE_GUILE

     std::cout << "----- debug:: calling parse_command_line() from setup_python()" << std::endl;
     command_line_data cld = parse_command_line(argc, argv);
     handle_command_line_data(cld);

     // BL says::should this still be here?
     run_state_file_maybe(); // run local 0-coot.state.py?

     // run_update_self_maybe(); nope.  Not yet.
     
#endif // USE_GUILE - not both start-up scripts

#if defined USE_PYGTK && !defined USE_GUILE_GTK

     // No more tips.
     
#endif // USE_PYGTK
     
#endif // USE_PYTHON  

}
