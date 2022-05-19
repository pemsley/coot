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
  
#ifdef USE_GUILE_GTK
  if (use_graphics_flag)
    coot_init_glue();
#endif	/* USE_GUILE_GTK */

  /* debugging voodoo */
  SCM_DEVAL_P = 1;
  SCM_BACKTRACE_P = 1;
  SCM_RECORD_POSITIONS_P = 1;

  // this is not the right place.  We need to be in
  std::string d1 = coot::util::append_dir_dir(PKGDATADIR, "scheme");
  std::string d2;
  std::string scheme_dir = d1;
  char *e = getenv(COOT_SCHEME_DIR);
  if (e) {
     d2 = e;
     scheme_dir = d2;
  }

  struct stat buf;
  int istat = stat(d2.c_str(), &buf);
  if (istat) {
     // fail
     scheme_dir = d1;
  }

  if (coot::is_directory_p(scheme_dir)) {

     // We should put our scheme files in $prefix/share/guile/{coot} and put a loader in
     // $prefix/share/guile/coot.scm

     std::string full_extra_load_path_cmd =  "(set! %load-path (cons \"";
     full_extra_load_path_cmd += scheme_dir;
     full_extra_load_path_cmd +=  "\" %load-path))";
     scm_c_eval_string(full_extra_load_path_cmd.c_str());

     std::string ee = "(lambda (key . args) ";
     ee += "(display (list \"Error in proc:\" key \" args: \" args)) (newline))";
     SCM handler = scm_c_eval_string(ee.c_str());
     std::string thunk_str = "(lambda ()";
     thunk_str += "(let ((f (%search-load-path \"coot.scm\")))";
     thunk_str += "   (if (eq? f #f)";
     thunk_str += "       (begin";
     thunk_str += "         (display ";
     thunk_str += "            (list \"Error finding coot.scm in %load-path\" %load-path))";
     thunk_str +="             (newline))";
     thunk_str += "      (load f)))) ; if we found it in the path";


     SCM thunk = scm_c_eval_string(thunk_str.c_str());

     /*  We need to get the value out of this, did the thunk fail or not?
    We want to be able to pass the value back to
    graphics_info_t::guile_gui_loaded_flag used in
    do_scripting_window_activate().  We now do that by setting a flag
    in graphics_info_t by calling the scheme function
    (set-found-coot-gui).
 */
     scm_catch(SCM_BOOL_T, thunk, handler);

     std::string flag = "#f";
#ifdef USE_GUILE_GTK
     flag = "#t";
#endif

     if (! use_graphics_flag)
        flag = "#f";

     std::string l = "(lambda () (load-all-scheme " + flag + "))";
     thunk = scm_c_eval_string(l.c_str());
     scm_catch(SCM_BOOL_T, thunk, handler);
  }

  if (run_startup_scripts_state()) {
     try_load_scheme_extras_dir();
     try_load_dot_coot_and_preferences();
  }

  /* now handle the command line data */
   handle_command_line_data_argc_argv(argc, argv);

   run_command_line_scripts(); // i.e. -c '(do-something)'
   run_state_file_maybe();
   pre_load_rotamer_tables();

   if (use_graphics_interface_state()) {
     gtk_main();
  } else {
     short int python_at_prompt_flag = python_at_prompt_at_startup_state();
     if (python_at_prompt_flag)
        start_command_line_python_maybe(argv);
     else
        scm_shell(0, argv);		/* can't pass command line args such
                                           as --pdb --no-graphics etc. (guile
                                           doesn't understand them). */
  }

}


void my_wrap_scm_boot_guile(int argc, char** argv) { 

/* From libguile/init.h:  */
/* extern void scm_boot_guile (int argc, char **argv, */
/*                             void (*main_func) (void *closure, */
/*                                                int argc, */
/*                                                char **argv), */
/*                             void *closure); */
  
  scm_boot_guile(argc, argv, inner_main, NULL);

  printf("you should not see this, c_inner_main should have called exit(0)\n"); 

}

void try_load_dot_coot_and_preferences() {

   auto file_is_directory = [] (const std::string &file_name) {
                               struct stat s;
                               int status = stat(file_name.c_str(), &s);
                               if (status == 0) {            /* the file existed */
                                  std::cout << "file-is_directory here a" << std::endl;
                                  if (S_ISDIR(s.st_mode)) {
                                     std::cout << "file-is_directory here B" << std::endl;
                                     return true;
                                  } else {
                                  }
                               } else {
                                  std::cout << "file-is_directory here C" << std::endl;
                                  std::cout << "oops stating " << file_name << " return non-zero status" << std::endl;
                               }
                               std::cout << "file-is_directory here D" << std::endl;
                               return false;
                            };

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
            if (file_is_directory(fn)) {
               std::cout << "INFO:: Not Loading ~/.coot - it's a directory " << std::endl;
            } else {
               std::cout << "asdfasdfLoading ~/.coot" << std::endl;
               scm_c_primitive_load(fn.c_str());
            }
	 }
      }
   }
}


#endif /* USE_GUILE */
