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

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   // nothing
#else 

   std::vector<std::string> v;
   const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry_widget));
   GtkWidget *sort_button = lookup_widget(entry_widget, "fileselection_sort_button");
   if (sort_button) { 
      // std::cout << "Hooray! we found the sort button!\n";
      // usually, we do.
   } else { 
      std::cout << "Boo we failed to find the sort button!\n";
   } 
   std::string pre_directory = pre_directory_file_selection(sort_button);

   // so now we have pre_directory
   // 
   // std::cout << "DEBUG:: pre_directory: " << pre_directory << std::endl;
   GtkWidget *fileselection = lookup_widget(entry_widget, "run_script_fileselection");
   if (fileselection) {

      std::string file_name_glob = pre_directory;
      file_name_glob += "/";

      file_name_glob += text;
      glob_t myglob;
      int flags = 0;
      glob(file_name_glob.c_str(), flags, 0, &myglob);
      size_t count;

      char **p;
      for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) {
	 std::string f(*p);
	 // std::cout << "glob testing " << f << std::endl;
	 v.push_back(f);
      }
      globfree(&myglob);

      GtkWidget *fl = GTK_FILE_SELECTION(fileselection)->file_list;
      GtkTreeView *tv = GTK_TREE_VIEW(fl);
      GtkTreeModel *model = gtk_tree_view_get_model(tv);
      GtkTreeIter iter;
      gtk_list_store_clear(GTK_LIST_STORE(model));
      for (unsigned int i=0; i<v.size(); i++) {
	 gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	 gtk_list_store_set(GTK_LIST_STORE(model), &iter,
			    0, coot::util::file_name_non_directory(v[i]).c_str(),
			    -1);
      }


   } else { 
      std::cout << "ERROR:: couldn't find fileselection\n";
   } 
#endif // WINDOWS_MINGW
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


void fileselection_sort_button_clicked( GtkWidget *sort_button,
					GtkWidget  *file_list) {

   GtkOptionMenu *history_pulldown =
      GTK_OPTION_MENU(gtk_object_get_user_data(GTK_OBJECT(sort_button)));

   GList *dlist = gtk_container_children(GTK_CONTAINER(history_pulldown));
   GList *free_list = dlist;
   std::string pre_directory("");
   
   while (dlist) {
      // GtkWidget *list_item;
      // list_item = (GtkWidget *) (dlist->data);
      gchar *t = GTK_LABEL(dlist->data)->label;
      if (t != NULL) {
	 pre_directory = t; 
      } else {
	 std::cout << "null label t " << std::endl;
      } 
      dlist = dlist->next;
   }
   g_list_free(free_list);
   // now pre_directory is set...

   GtkTreeView *tv = GTK_TREE_VIEW(file_list);
   GList *collist = gtk_tree_view_get_columns(tv);
   for (GList *node = collist; node != NULL; node = g_list_next(node)) {
      GtkTreeViewColumn *col = GTK_TREE_VIEW_COLUMN(node->data);
      GtkWidget *label = gtk_tree_view_column_get_widget(col);
      std::string lab = gtk_label_get_text(GTK_LABEL(label));
      if (lab == "Files") {
	 GtkTreeModel *model = gtk_tree_view_get_model(tv);
	 GtkTreeIter iter;
	 // GtkListStore *liststore; // model and liststore are the same thing?
	 coot::file_attribs_info_t file_infos;
	 file_infos.directory_prefix = pre_directory;

	 // fill file_mtimes of the file_infos:
	 gtk_tree_model_foreach(GTK_TREE_MODEL(model),
				fileselection_sort_button_foreach_func,
				&file_infos);

	 // sort the files by date
	 std::sort(file_infos.file_mtimes.begin(), file_infos.file_mtimes.end(),
		   compare_mtimes);

	 // debug
// 	 std::cout << "There are " << file_infos.file_mtimes.size()
// 		   << " file attributes" << std::endl;
// 	 for (int i=0; i<file_infos.file_mtimes.size(); i++)
// 	    std::cout << file_infos.file_mtimes[i].file << std::endl;

	 gtk_list_store_clear(GTK_LIST_STORE(model));

	 for (unsigned int i=0; i<file_infos.file_mtimes.size(); i++) {
	    gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	    gtk_list_store_set(GTK_LIST_STORE(model), &iter,
	    		       0, file_infos.file_mtimes[i].file.c_str(), -1);
	 }
      }
	 
   }
   g_list_free(collist);
   
}


void
on_filename_filter_toggle_button_toggled(GtkButton       *button,
					 gpointer         user_data)
{

   // user_data is used to set the glob filter 
   int int_user_data = GPOINTER_TO_INT(user_data);
   int data_type = int_user_data & 31; // lower 5 bits

   //    std::cout << "DEBUG:: data_type:" << data_type << std::endl;
   
//    std::cout << "more user data:" << (32>>5) << " " << (int_user_data>>5)
// 	     << std::endl;
   int file_selection_type = data_type;

   // We need to add text to the string of the dictectory we are in
   // (pre_directory), so first we need to find pre_directory (as per
   // fileselection_sort_button_clicked()
   // 
   GtkWidget *sort_button = lookup_widget(GTK_WIDGET(button),
					  "fileselection_sort_button");
   if (sort_button) { 
      // std::cout << "Hooray! we found the sort button!\n";
      // usually, we do.
   } else { 
      std::cout << "Boo we failed to find the sort button!\n";
   } 
   std::string pre_directory = pre_directory_file_selection(sort_button);
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button,
							    file_selection_type);
   
   std::vector<std::string> v;
   
   if (fileselection) { 
      if (GTK_TOGGLE_BUTTON(button)->active) { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child), "Unfilter");
	 
	 // so now we have pre_directory
	 // 
	 std::vector<std::string> v = filtered_by_glob(pre_directory, data_type);

	 GtkWidget *fl = GTK_FILE_SELECTION(fileselection)->file_list;
	 GtkTreeView *tv = GTK_TREE_VIEW(fl);

	 GtkTreeModel *model = gtk_tree_view_get_model(tv);
	 GtkTreeIter iter;

	 // OK, one more go.  Delete everything and add the files that pass:
	 
	 gtk_list_store_clear(GTK_LIST_STORE(model));
	 std::vector<std::string> file_vec = filtered_by_glob(pre_directory, data_type);
	 for (unsigned int i=0; i<file_vec.size(); i++) {
	    gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	    gtk_list_store_set(GTK_LIST_STORE(model), &iter,
			       0, coot::util::file_name_non_directory(file_vec[i]).c_str(),
			       -1);
	 }
	 
      } else { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Filter");
	 gtk_file_selection_set_filename(GTK_FILE_SELECTION(fileselection),
					 (pre_directory + "/").c_str());
      }
   } else {
      std::cout << "no fileselection found from sort button\n";
   }
}
