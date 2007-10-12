/* src/c-inner-main.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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
#include <guile/gh.h>
#endif

/* On some strange systems still no definition of NULL is found.  Sigh!  */
#ifndef NULL
# if defined __STDC__ && __STDC__
#  define NULL ((void *) 0)
# else
#  define NULL 0
# endif
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h> // for getenv
#include <glib.h>

/* for stat */
#include <sys/types.h>
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#endif

#include  <gtk/gtk.h> /* needed for c-interface.h */
#include "c-interface.h" /* needed for set_guile_gui_loaded_flag() and run_command_line_scripts() */
#include "c-inner-main.h"

#if defined (USE_GUILE) || defined (USE_PYTHON)
extern void SWIG_init();

/* #define SWIG_init    init_coot // this fails to init/register guile functions */

void coot_init_glue();
#endif

/* from support.c gchar->char */
char*
does_file_exist                      (const char     *directory,
				      const char     *filename)
{
  char *full_filename;
  struct stat s;
  gint status;

  full_filename = (gchar*) g_malloc (strlen (directory) + 1
                                     + strlen (filename) + 1);
  strcpy (full_filename, directory);
  /*  strcat (full_filename, G_DIR_SEPARATOR_S); */
  strcat (full_filename, "/");
  strcat (full_filename, filename);

  status = stat (full_filename, &s);
  if (status == 0 && S_ISREG (s.st_mode))
    return full_filename;
  g_free (full_filename);
  return NULL;
}




#ifdef USE_GUILE

void
c_inner_main(void *closure, int argc, char** argv) { 

/*    export commands */
   
  char *filename = ".coot"; 
  char *directory;
  char* check_file;
  char *gui_lib = 0; 
  gchar* full_extra_load_path_cmd; 
  char *pre_str = "(set! %load-path (cons \""; 
  char *post_str = "\" %load-path))"; 
  SCM handler; 
  SCM thunk; 
  char *thunk_str; 
  char *tmp_str;
  struct stat buf;
  int istat;
  short int use_graphics_flag = use_graphics_interface_state();

  // printf("::::::::::::::::: c_inner_main() SWIG_init\n");
  SWIG_init();   


#ifdef USE_GUILE_GTK
  if (use_graphics_flag)
    coot_init_glue();
#endif	/* USE_GUILE_GTK */

  /* debugging voodoo */
  SCM_DEVAL_P = 1;
  SCM_BACKTRACE_P = 1;
  SCM_RECORD_POSITIONS_P = 1;
#if (SCM_MAJOR_VERSION == 1) || (SCM_MINOR_VERSION == 6)
  SCM_RESET_DEBUG_MODE;
#else 
  /* version 1.8 or above */
#endif

/*   First try to getenv COOT_SCHEME_DIR.
     Recall that COOT_SCHEME_DIR is #defined to be "COOT_SCHEME_DIR" */
  gui_lib = getenv(COOT_SCHEME_DIR);

  if (gui_lib) { 
    istat = stat(gui_lib, &buf);
  }

  if (istat || !gui_lib) {	/* failed or env var was not defined*/

/*    so that failed, let's build the hard-coded default directory and
      stat it.  If it exists let's set gui_lib to that value. 
*/
    tmp_str = PKGDATADIR;
    tmp_str = (char *) malloc (strlen(tmp_str) + 9);
    strcpy (tmp_str, PKGDATADIR);
    strcat (tmp_str, "/");	/* something else for Windwoes? */
    strcat (tmp_str, "scheme");
    gui_lib = tmp_str;
    istat = stat(gui_lib, &buf);

    if (istat) { 			/* failed */
      gui_lib = NULL;
      printf("WARNING: Scheme directory: %s not found\n", tmp_str);
      printf("         and environment variable %s not defined\n", 
	     COOT_SCHEME_DIR);
      printf("         or directory did not exist.\n");
    } 
  }

  if (gui_lib) {
    printf("Loading scheme files from %s\n", gui_lib);
    full_extra_load_path_cmd = (gchar *) g_malloc(strlen(pre_str) + 1 +
						  strlen(gui_lib) + 1 +
						  strlen(post_str));
    strcpy (full_extra_load_path_cmd, pre_str);
    strcat (full_extra_load_path_cmd, gui_lib); 
    strcat (full_extra_load_path_cmd, post_str); 
    /* printf("evaluating: %s\n", full_extra_load_path_cmd); */
    scm_c_eval_string(full_extra_load_path_cmd);
    g_free(full_extra_load_path_cmd); 

    handler = scm_c_eval_string("(lambda (key . args) "
	"(display (list \"Error in proc:\" key \" args: \" args)) (newline))"); 

    thunk_str = "(lambda () (let ((f (%search-load-path \"coot.scm\")))"
                               "(if (eq? f #f)"
                                   "(begin"
                                      "(display "
                                         "(list \"Error finding coot.scm in\" "
                                                "%load-path))"
                                      "(newline))"
                                   "(load f)))) ; if we found it in the path";


    thunk = scm_c_eval_string(thunk_str); 

/*  We need to get the value out of this, did the thunk fail or not?
    We want to be able to pass the value back to
    graphics_info_t::guile_gui_loaded_flag used in
    do_scripting_window_activate().  We now do that by setting a flag
    in graphics_info_t by calling the scheme function
    (set-found-coot-gui).
 */
    scm_catch(SCM_BOOL_T, thunk, handler);

    if (use_graphics_flag) { 
#ifdef USE_GUILE_GTK
      thunk_str = "(lambda () (load-all-scheme #t))";
#else 
      thunk_str = "(lambda () (load-all-scheme #f))";
#endif // USE_GUILE_GTK
    } else { 
      thunk_str = "(lambda () (load-all-scheme #f))";
    }
    thunk = scm_c_eval_string(thunk_str); 
    scm_catch(SCM_BOOL_T, thunk, handler);
  }

/* And now read the users own initialization code */
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
  directory = getenv("COOT_HOME"); 
#else
  directory = getenv("HOME"); 
#endif
  if (directory) {  
     check_file = does_file_exist(directory, filename); 

     if (check_file) { 
       printf("Loading ~/.coot..."); 
       scm_c_primitive_load(check_file); 
       printf("done.\n");
     }
  
  run_command_line_scripts();	/* this may turn off run-state-file */

  run_state_file_maybe();

/* for now, make the user start the listener explictly. */
/*   make_socket_listener_maybe(); */

  /* tips gui? */
  if (gui_lib) {
    if (use_graphics_flag) { 
      thunk_str = "(lambda () (tips-gui))";
      thunk = scm_c_eval_string(thunk_str); 
      scm_catch(SCM_BOOL_T, thunk, handler);
    }
  }

/*   scm_shell(argc, argv); */

   }

  /* option for Ralf */
  if (use_graphics_interface_state())
    /* OK! over to you Gtk+! */
    gtk_main(); 
  else 
    scm_shell(0, argv);		/* can't pass command line args such
				   as --pdb --no-graphics etc. (guile
				   doesn't understand them). */
}
#endif


#ifdef USE_GUILE

void c_wrapper_scm_boot_guile(int argc, char** argv) { 

/* From libguile/init.h:  */
/* extern void scm_boot_guile (int argc, char **argv, */
/*                             void (*main_func) (void *closure, */
/*                                                int argc, */
/*                                                char **argv), */
/*                             void *closure); */
  
  scm_boot_guile(argc, argv, c_inner_main, NULL);

  printf("you should not see this, c_inner_main should have called exit(0)\n"); 

} 


#endif /* USE_GUILE */
