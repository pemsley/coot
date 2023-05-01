/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
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

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#include <algorithm>
#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

#ifdef USE_GUILE
#include <cstdio>  // for std::FILE
#include <libguile.h>
#endif 

#if !defined(WINDOWS_MINGW) && !defined(_MSC_VER)
#include <glob.h> // for globbing.
#endif

#include <iostream>
#include <string>
#include <vector>
#include "utils/coot-utils.hh"

#include "support.h"  // for lookup_widget
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh" // for str_mtime
#include "graphics-info.h"  // for go to atom callabacks

#include "widget-headers.hh"


// Function currently not used
// 
gboolean filename_passed_filter(const std::string &file_name, int filter_type) {

   gboolean r = 0;

   std::vector<std::string> globs;
   if (filter_type == 0) 
      globs = *graphics_info_t::coordinates_glob_extensions;
   if (filter_type == 1) 
      globs = *graphics_info_t::data_glob_extensions;
   if (filter_type == 2) 
      globs = *graphics_info_t::map_glob_extensions;
   if (filter_type == 3) 
      globs = *graphics_info_t::dictionary_glob_extensions;

   std::string extension = coot::util::file_name_extension(file_name);
   for (unsigned int i=0; i<globs.size(); i++) {
      // std::cout << "comparing " << extension << " with " << globs[i] << std::endl;
      if (extension == globs[i]) {
	 r = 1;
	 break;
      }
   }
   return r;
}


// push back mtimes on to the vector of str_mtimes passed as user data.
gboolean
fileselection_sort_button_foreach_func
 (GtkTreeModel *model,
  GtkTreePath  *path,
  GtkTreeIter  *iter,
  gpointer      user_data)
{
   gchar *file_name, *tree_path_str;
   gtk_tree_model_get (model, iter,
		       0, &file_name,
		       -1);
   tree_path_str = gtk_tree_path_to_string(path);
   coot::file_attribs_info_t *file_attribs = (coot::file_attribs_info_t *) (user_data);
   struct stat buf;
   std::string directory_prefix = file_attribs->directory_prefix;
   std::string full_file_path =
      coot::util::append_dir_file(directory_prefix, file_name);
   int status = stat(full_file_path.c_str(), &buf);
   if (status == 0) { 
      time_t mtime = buf.st_mtime;
      file_attribs->file_mtimes.push_back(coot::str_mtime(file_name, mtime));
   } else {
      std::cout << " stat returns " << status << std::endl;
   } 
   //    g_print ("Row %s: %s\n", tree_path_str, file_name);
      
   g_free(tree_path_str);
   g_free(file_name); /* gtk_tree_model_get made copies of       */
   return FALSE;
}

// There was a CR in the file name selection entry in the Run Script
// file selection.
// 
void
handle_filename_filter_gtk2(GtkWidget *entry_widget) {

   std::cout << "delete this function? handle_filename_filter_gtk2" << std::endl;
}


// Function currently not used
gboolean
fileselection_filter_button_foreach_func
 (GtkTreeModel *model,
  GtkTreePath  *path,
  GtkTreeIter  *iter,
  gpointer      user_data)
{
   gchar *file_name;
   gtk_tree_model_get (model, iter,
		       0, &file_name,
		       -1);
   int file_selection_filter_type = GPOINTER_TO_INT(user_data);

   if (!filename_passed_filter(file_name, file_selection_filter_type)) {
      std::cout << file_name << " to be deleted" << std::endl;
      // delete it
      if (gtk_tree_model_get_iter(model, iter, path)) {
	 // continue, normal case
	 std::cout << file_name << " deleted" << std::endl;
	 gtk_list_store_remove(GTK_LIST_STORE(model), iter);
      }
   } else {
      std::cout << "    " << file_name << " keep it" << std::endl;
   }
   return FALSE; 
}



void
on_filename_filter_toggle_button_toggled(GtkButton       *button,
					 gpointer         user_data)
{

}
