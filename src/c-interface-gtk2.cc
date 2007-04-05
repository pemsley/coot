

#include <gtk/gtk.h>
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)


# else 
#include <iostream>
#include <string>

#include "c-interface.h"

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
	 gboolean istat = gtk_tree_model_get_iter_first(model, &iter); // inits iter
	 gint ich = gtk_tree_model_iter_n_children(model, &iter);
	 std::cout << "n children: " << ich << std::endl;
	 while (istat) {

	    gchar *name;
	    gchar *COL_NAME = "Files";
	    std::cout << "next..!" ;
	    // will crash...	    
// 	    gtk_tree_model_get(model, &iter, COL_NAME, &name, -1);
// 	    g_free(name);
	    istat = gtk_tree_model_iter_next(model, &iter);
	 }
      }
   }
   g_list_free(collist);
   
}
#endif // GTK2
