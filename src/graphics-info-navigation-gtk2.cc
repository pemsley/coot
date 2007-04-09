


#include <gtk/gtk.h>
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)


# else 

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

#include "graphics-info.h"

// a static
//
// Recall that the tree is created in c-interface.cc's fill_go_to_atom_window().
void
graphics_info_t::fill_go_to_atom_residue_tree_gtk2(GtkWidget *gtktree) {

   std::cout << "filling residue tree!" << std::endl;
   
   std::string button_string;
   graphics_info_t g;

   g.go_to_atom_residue(); // sets values of unset (magic -1) go to
			   // atom residue number.

   std::vector<coot::model_view_atom_tree_chain_t> residue_chains = 
      molecules[g.go_to_atom_molecule()].model_view_residue_tree_labels();

   // so, clear the current tree:
   GtkTreeView *tv = GTK_TREE_VIEW(gtktree);
   GtkTreeModel *model = gtk_tree_view_get_model(tv);
   gtk_tree_store_clear(GTK_TREE_STORE(model));

}

// static
void
graphics_info_t::on_go_to_atom_residue_tree_selection_changed (GtkTreeView *gtklist,
							       gpointer user_data) {

   std::cout << "tree selection changed - fill me!" << std::endl;

}


#endif // GTK_MAJOR_VERSION etc
