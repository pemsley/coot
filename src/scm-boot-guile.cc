/* src/scm-boot-guile.cc
 * 
 * Copyright 2004 by the University of York
 * Copyright 2016 by Medical Research Council
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

#ifdef USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"  // before guile includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <string>
#include <iostream>
#include <cstddef> // needed for gmp.h I think.

#include <libguile.h>
#include <glib.h>
#include <gtk/gtk.h>

#include <sys/types.h>
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#endif

#include "utils/coot-utils.hh"
#include "scm-boot-guile.hh"
#include "boot-python.hh"
#include "c-interface.h"
#include "c-interface-preferences.h"

#include "startup-scripts.hh"

#include "coot-init-glue.hh"
#include "pre-load.hh"

#include <glob.h>
#include "graphics-info.h"

#include "command-line-extern.hh"

void inner_main(void *closure, int argc, char **argv) {

   short int use_graphics_flag = use_graphics_interface_state();

   SWIG_init();

   std::string handler_string = "(lambda (key . args) ";
   handler_string += "(display (list \"Error in proc:\" key \" args: \" args)) (newline))";
   SCM handler = scm_c_eval_string(handler_string.c_str());

   std::string coot_scm;
   coot_scm = coot::package_data_dir();
   coot_scm = coot::util::append_dir_dir(coot_scm, "scheme");
   coot_scm = coot::util::append_dir_file(coot_scm, "coot.scm");
   std::cout << "debug:: coot_scm " << coot_scm << std::endl;

   std::string thunk_str = "(lambda () (load \"" + coot_scm + "\"))";
   SCM thunk = scm_c_eval_string(thunk_str.c_str());

   scm_catch(SCM_BOOL_T, thunk, handler);

  if (run_startup_scripts_state()) {
     try_load_scheme_extras_dir();
     try_load_dot_coot_and_preferences();
  }

   run_command_line_scripts(); // i.e. -c '(do-something)'
   run_state_file_maybe();
   // pre_load_rotamer_tables(); what's this?

   if (use_graphics_flag)
      gtk_main();
   else
      scm_shell(argc, argv);

}


void my_wrap_scm_boot_guile(int argc, char** argv) { 

/* From libguile/init.h:  */
/* extern void scm_boot_guile (int argc, char **argv, */
/*                             void (*main_func) (void *closure, */
/*                                                int argc, */
/*                                                char **argv), */
/*                             void *closure); */

   scm_boot_guile(argc, argv, inner_main, NULL);

   std::cout << "you should not see this, inner_main should never return\n";

}

void try_load_dot_coot_and_preferences() {

   // python versionn in coot-setup-python.cc

   bool run_startup_scripts_flag = run_startup_scripts_state();

   std::string directory = coot::get_home_dir();

   if (! directory.empty()) {

      if (run_startup_scripts_flag) {

	 std::string preferences_dir = graphics_info_t::add_dir_file(directory, ".coot-preferences");
	 std::string full_path_coot_pref_scm = graphics_info_t::add_dir_file(preferences_dir,
									     "coot-preferences.scm");
	 struct stat buff;
	 int preferences_dir_status = stat(preferences_dir.c_str(), &buff);
     
	 if (preferences_dir_status != 0) { 
	    std::cout << "INFO:: preferences directory " << preferences_dir 
		      << " does not exist. Won't read .scm preferences." << std::endl;;
	 } else {

	    char **p;
	    size_t count;
	    int flags = 0;
	    glob_t myglob;
	    std::string glob_file = preferences_dir + "/*.scm";
	    
	    glob(glob_file.c_str(), flags, 0, &myglob);

	    for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) { 
	       char *preferences_file = (*p);
	       if (preferences_file) {
		  std::string preferences_file_str(preferences_file);
		  // if python prefered, don't load coot-preferences.scm
		  if (preferences_file_str == full_path_coot_pref_scm) {
		     if (prefer_python()) {
			// std::cout << "skip coot-prefences.scm " << std::endl;
		     } else {
			std::cout << "load " << preferences_file_str << std::endl;
			scm_c_primitive_load(preferences_file);
		     }
		  } else {
		     // happy path
		     std::cout << "load " << preferences_file << std::endl;
		     scm_c_primitive_load(preferences_file);
		  } 
	       }
	    }
	 }

         // update preferences
         make_preferences_internal();
	 // Now ~/.coot

	 std::string fn = coot::util::append_dir_file(directory, ".coot");
	 if (coot::file_exists(fn)) {
	    std::cout << "Loading ~/.coot" << std::endl;
	    scm_c_primitive_load(fn.c_str()); 
	 }
      }
   }
}


#endif /* USE_GUILE */
