
#ifdef USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"  // before guile includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <string>
#include <iostream>

#include <libguile.h>
#include <glib.h>
#include  <gtk/gtk.h>


#include "utils/coot-utils.hh"
#include "scm-boot-guile.hh"
#include "c-interface.h"

#include "startup-scripts.hh"

// void coot_init_glue(); // Hmmm... needs fixing?

#include "coot-init-glue.hh"

#include <glob.h>
#include "startup-scripts.hh"
#include "graphics-info.h"

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
     
     std::string e = "(lambda (key . args) ";
     e += "(display (list \"Error in proc:\" key \" args: \" args)) (newline))"; 
     SCM handler = scm_c_eval_string(e.c_str());
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

void start_command_line_python_maybe(char **argv) {

#ifdef USE_PYTHON
   
#if PY_MAJOR_VERSION > 2 
   Py_Main(0, argv);
#else
#if PY_MINOR_VERSION > 2 
   Py_Main(0, argv);
#endif     // PY_MINOR_VERSION
#endif     // PY_MAJOR_VERSION

     //  Skip initialization registration of signal handlers, useful
     //  when Python is embedded. Version 2.4 or later. Thanks Stuart
     //  McNicholas for letting me know about this.
     //
     // Question: Do we need to check that we are not using command
     // line python no graphics before we use this?
     // 
#if PY_MAJOR_VERSION > 2
   Py_InitializeEx(0);
#endif     
#if PY_MAJOR_VERSION == 2
#if PY_MINOR_VERSION > 3
   Py_InitializeEx(0);
#endif     
#endif     
#endif     
}

void try_load_dot_coot_and_preferences() {

   // python versionn in coot-setup-python.cc

   bool run_startup_scripts_flag = run_startup_scripts_state();
   
   char *d1 = getenv("COOT_HOME");
   char *d2 = getenv("HOME");

   std::string directory;
   
   if (d1) {
      directory = d1;
   } else {
      if (d2) {
	 directory = d2;
      }
   }

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
			std::cout << "skip" << std::endl;
			// skip
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
      }
   }
}


#endif /* USE_GUILE */
