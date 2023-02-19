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



void
add_python_scripting_entry_completion(GtkWidget *entry) {

   // std::cout << "======================= add_python_scripting_entry_completion() " << std::endl;

   graphics_info_t g; // for history

   GtkEntryCompletion *completion = gtk_entry_completion_new();
   gtk_entry_completion_set_popup_completion(completion, TRUE);
   gtk_entry_completion_set_text_column(completion, 0);
   gtk_entry_completion_set_minimum_key_length(completion, 2);
   gtk_entry_set_completion(GTK_ENTRY(entry), completion);

   std::vector<std::string> completions;
   std::vector<std::string> module_coot_completions;
   std::vector<std::string> module_coot_utils_completions;

   PyErr_Clear();
   // I'd like to do the dir(coot) using the API, something like:
   // PyObject *object = PyObject_Dir(module_py);
   // but I don't know how to get the module_py
   // PyObject *PyState_FindModule(PyModuleDef *def) - what is def?
   //
   PyObject *object = safe_python_command_with_return("dir(coot)");
   std::cout << "object " << object << std::endl;
   if (object) {
      std::string module_name = "coot";
      module_coot_completions.reserve(2000);

      // auto tp_0 = std::chrono::high_resolution_clock::now();
      if (PyList_Check(object)) {
         Py_ssize_t n = PyList_Size(object);
         for (Py_ssize_t i=0; i<n; i++) {
            PyObject *item = PyList_GetItem(object, i);
            if (PyUnicode_Check(item)) {
               std::string str = module_name + std::string(".") +
                  PyBytes_AS_STRING(PyUnicode_AsUTF8String(item)) + std::string("(");
               module_coot_completions.push_back(str);
            }
         }
      }
      // auto tp_1 = std::chrono::high_resolution_clock::now();
      // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      // std::cout << "Timings: dir coot " << d10 << " microseconds" << std::endl;
   }

   // auto tp_2 = std::chrono::high_resolution_clock::now();
   object = safe_python_command_with_return("dir(coot_utils)");
   if (object) {
      if (PyList_Check(object)) {
         std::string module_name = "coot_utils";
         Py_ssize_t n = PyList_Size(object);
         for (Py_ssize_t i=0; i<n; i++) {
            PyObject *item = PyList_GetItem(object, i);
            if (PyUnicode_Check(item)) {
               std::string str = module_name + std::string(".") +
                  PyBytes_AS_STRING(PyUnicode_AsUTF8String(item)) + std::string("(");
               module_coot_utils_completions.push_back(str);
            }
         }
      }
      // auto tp_3 = std::chrono::high_resolution_clock::now();
      // auto d32 = std::chrono::duration_cast<std::chrono::microseconds>(tp_3 - tp_2).count();
      // std::cout << "Timings: dir coot_utils " << d32 << " microseconds" << std::endl;
   }

   // command history
   std::vector<std::string> chv = g.command_history.commands;
   chv = g.command_history.unique_commands(); // there *were* unique already

   // add together the completions
   completions.push_back("import coot");
   completions.push_back("import coot_utils");
   completions.insert(completions.end(), chv.begin(),                           chv.end());
   completions.insert(completions.end(), module_coot_completions.begin(),       module_coot_completions.end());
   completions.insert(completions.end(), module_coot_utils_completions.begin(), module_coot_utils_completions.end());

   // maybe only once!
   GtkListStore *store = gtk_list_store_new(1, G_TYPE_STRING);
   GtkTreeIter iter;

   for (unsigned int i=0; i<completions.size(); i++) {
      gtk_list_store_append( store, &iter );
      std::string c = completions[i];
      // std::cout << "adding to gtk-completion: " << c << std::endl;
      gtk_list_store_set( store, &iter, 0, c.c_str(), -1 );
   }

   gtk_entry_completion_set_model(completion, GTK_TREE_MODEL(store));

}

// used in both the button and the entry callback.
#include "python-scripting-gui.hh"

// extern "C" G_MODULE_EXPORT
gboolean
on_python_window_entry_key_press_event(GtkWidget   *entry,
                                       GdkEventKey *event,
                                       gpointer     user_data) {

   // std::cout << "Key-press " << event->keyval << " " << event->state << std::endl;

   if (event->keyval == GDK_KEY_Up) {
      // allow history scroll up if Ctrl is also pressed
      if (event->state & GDK_CONTROL_MASK) {
         graphics_info_t g;
         std::string t = g.command_history.get_previous_command();
         // std::cout << "previous-command: \"" << t << "\"" << std::endl;
         gtk_entry_set_text(GTK_ENTRY(entry), t.c_str());
         return TRUE;
      }
   }
   if (event->keyval == GDK_KEY_Down) {
      // Allow history scroll down if Ctrl is also pressed.
      // Without Ctrl, this co-incides with navigating into and up and down in
      // the completion dialog.
      if (event->state & GDK_CONTROL_MASK) {
         graphics_info_t g;
         std::string t = g.command_history.get_next_command();
         gtk_entry_set_text(GTK_ENTRY(entry), t.c_str());
         return TRUE;
      }
   }

   if (event->keyval == GDK_KEY_Return) {
      // GdkModifierType state;
      // gboolean es = gtk_get_current_event_state(&state);
      int state = event->state;
      if (state & GDK_CONTROL_MASK) {
         // std::cout << "Return with control!" << std::endl;
         run_python_scripting_window_entry_text(entry); // clears the entry
         return TRUE;
      } else {
         // std::cout << "Return!" << std::endl;
      }
   }

   // std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(display_python(object)));
   // std::cout << s << std::endl;

   return FALSE;
}

void
run_python_scripting_window_entry_text(GtkWidget *entry) {

   if (entry) {
      const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
      std::string entry_text_as_string(entry_text); // important to make a copy
      PyRun_SimpleString(entry_text);

      // clear the entry
      gtk_entry_set_text(GTK_ENTRY(entry), "");

      graphics_info_t g;
      g.command_history.add_to_history(entry_text_as_string);

   }
}

// We want to evaluate the string when we get a carriage return
// in this entry widget
void
setup_python_window_entry(GtkWidget *entry) {

#ifdef USE_PYTHON

#if 0
   // setup the python entry in entry callback code here...

   g_signal_connect(G_OBJECT(entry), "activate",
                    G_CALLBACK(python_window_enter_callback),
                    (gpointer) entry);

#endif

   g_signal_connect(G_OBJECT(entry), "key-press-event",
                    G_CALLBACK(on_python_window_entry_key_press_event),
                    (gpointer) entry);

   add_python_scripting_entry_completion(entry);

#endif // USE_PYTHON
}

#ifdef USE_PYTHON
void python_window_enter_callback(GtkWidget *widget,
                                  GtkWidget *entry ) {

   std::cout << "============================ python_window_enter_callback()! " << std::endl;

#if 0 //20230128-PE - now we have completions, an "enter" on a completion
      // entry calls this function - we don't want that.

   const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry));
   std::string entry_text_as_string(entry_text); // important to make a copy
   PyRun_SimpleString(entry_text);
   
   // clear the entry
   gtk_entry_set_text(GTK_ENTRY(entry), "");
   
   graphics_info_t g;
   g.command_history.add_to_history(entry_text_as_string);
#endif

}
#endif



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
      initcoot_python_gobject(); // this is not a good name for this function. We need to say
                                 // this this is the module that wraps the glue to get
                                 // the status-bar, menu-bar etc.
      PyImport_ImportModule("coot_utils"); // this imports coot_gui (which seems wrong)
      PyImport_ImportModule("extensions");
   }
   PyErr_PrintEx(0);

   if (graphics_info_t::run_startup_scripts_flag) {
      // setup_python_window_entry() is done only once now. So do it here.
      if (graphics_info_t::use_graphics_interface_flag) {
         GtkWidget *python_entry = widget_from_builder("python_window_entry");
         setup_python_window_entry(python_entry); // USE_PYTHON and USE_GUILE used here
      }
   }
   
   std::string home_directory = coot::get_home_dir();
   try_load_dot_coot_py_and_python_scripts(home_directory);

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
	 // load all .py files
	 glob_t myglob;
	 int flags = 0;
	 //std::string glob_patt = "/*.py";
	 std::string glob_file = startup_scripts_dir;
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
            std::string psf = coot::util::file_name_non_directory(preferences_script);
            if (exclude_py_files.find(psf) == exclude_py_files.end()) {
               bool done = false;
               if (preferences_script.length() > 6) {
                  if (psf.substr(0,6) == "xenops") {
                     done = true;
                     xenops_scripts.push_back(preferences_script);
                  }
               }
               if (psf.length() > 6) {
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
         for(const auto &script_fn : xenops_scripts) {
            run_python_script(script_fn.c_str());
         }
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
