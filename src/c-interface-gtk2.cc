

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

#include "c-interface.h"
#include "cc-interface.hh" // for str_mtime

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

void fileselection_sort_button_clicked( GtkWidget *sort_button,
					GtkWidget  *file_list) {


   // Play code.  GtkTreeView is complex.
   
   std::cout << "---- sort button pressed!" << std::endl;

   GtkTreeView *tv = GTK_TREE_VIEW(file_list);
   std::cout << "tree model: " << gtk_tree_view_get_model(tv)
	     << std::endl;
   
   GList *collist = gtk_tree_view_get_columns(tv);
   for (GList *node = collist; node != NULL; node = g_list_next(node)) {
      std::cout << "a column" << std::endl;
      GtkTreeViewColumn *col = GTK_TREE_VIEW_COLUMN(node->data);
      GtkWidget *label = gtk_tree_view_column_get_widget(col);
      std::string lab = gtk_label_get_text(GTK_LABEL(label));
      std::cout << "label: " << lab << std::endl;
      if (lab == "Files") {
	 GtkTreeModel *model = gtk_tree_view_get_model(tv);
	 GtkTreeIter iter;
	 GtkListStore *liststore;
	 std::vector<str_mtime> file_attr_vec;

	 gtk_tree_model_foreach(GTK_TREE_MODEL(model),
				fileselection_sort_button_foreach_func, &file_attr_vec);

	 // sort the files by date
	 std::sort(file_attr_vec.begin(), file_attr_vec.end(), compare_mtimes);
	 for (int i=0; i<file_attr_vec.size(); i++) {
	    std::cout << file_attr_vec[i].file << std::endl;
	 }
      }
	 
   }
   g_list_free(collist);
   
}
#endif // GTK2
