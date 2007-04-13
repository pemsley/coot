/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define sleep Sleep
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER


#include <gtk/gtk.h>
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)


# else 
#include <iostream>
#include <string>
#include <vector>
#include "coot-utils.hh"

#include "support.h"  // for lookup_widget
#include "c-interface.h"
#include "cc-interface.hh" // for str_mtime
#include "graphics-info.h"  // for go to atom callabacks

gboolean filename_passed_filter(const std::string &file_name, int filter_type) {

   if (file_name.substr(0,1) == "c") {
      return 1;
   } else {
      return 0;
   }
}


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
   std::vector<str_mtime> *file_vec_p = (std::vector<str_mtime> *)(user_data);
   struct stat buf;
   int status = stat(file_name, &buf);
   if (status == 0) { 
      time_t mtime = buf.st_mtime;
      file_vec_p->push_back(str_mtime(file_name, mtime));
   }
   //    g_print ("Row %s: %s\n", tree_path_str, file_name);
      
   g_free(tree_path_str);
   g_free(file_name); /* gtk_tree_model_get made copies of       */
   return FALSE;
}

gboolean
fileselection_filter_button_foreach_func
 (GtkTreeModel *model,
  GtkTreePath  *path,
  GtkTreeIter  *iter,
  gpointer      user_data)
{
   gchar *file_name, *tree_path_str;
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


   GtkTreeView *tv = GTK_TREE_VIEW(file_list);
   GList *collist = gtk_tree_view_get_columns(tv);
   for (GList *node = collist; node != NULL; node = g_list_next(node)) {
      GtkTreeViewColumn *col = GTK_TREE_VIEW_COLUMN(node->data);
      GtkWidget *label = gtk_tree_view_column_get_widget(col);
      std::string lab = gtk_label_get_text(GTK_LABEL(label));
      if (lab == "Files") {
	 GtkTreeModel *model = gtk_tree_view_get_model(tv);
	 GtkTreeIter iter;
	 GtkListStore *liststore; // model and liststore are the same thing?
	 std::vector<str_mtime> file_attr_vec;

	 gtk_tree_model_foreach(GTK_TREE_MODEL(model),
				fileselection_sort_button_foreach_func,
				&file_attr_vec);

	 // sort the files by date
	 std::sort(file_attr_vec.begin(), file_attr_vec.end(), compare_mtimes);
// 	 for (int i=0; i<file_attr_vec.size(); i++)
// 	    std::cout << file_attr_vec[i].file << std::endl;
	 gtk_list_store_clear(GTK_LIST_STORE(model));

	 for (int i=0; i<file_attr_vec.size(); i++) {
	    gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	    gtk_list_store_set(GTK_LIST_STORE(model), &iter,
			       0, file_attr_vec[i].file.c_str(), -1);
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
   int data_type = GPOINTER_TO_INT(user_data);

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
   GtkWidget *fileselection = lookup_file_selection_widgets(sort_button);
   
   std::vector<std::string> v;
   
   if (fileselection) { 
      if (GTK_TOGGLE_BUTTON(button)->active) { 
	 gtk_label_set_text(GTK_LABEL(GTK_BIN(button)->child),"Unfilter");
	 
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
	 for (int i=0; i<file_vec.size(); i++) {
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



#endif // GTK2

