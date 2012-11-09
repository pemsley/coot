
#ifndef CC_INTERFACE_MOGUL_HH
#define CC_INTERFACE_MOGUL_HH

#include "mogul-interface.hh"

void show_mogul_geometry_dialog(const coot::mogul &m);
GtkWidget *wrapped_create_mogul_geometry_dialog(const coot::mogul &m);

namespace coot {
   GtkCellRenderer *
   mogul_results_add_cell_renderer(GtkTreeView *tree_view,
					 GtkTreeStore *tree_store,
					 const std::string &column_title,
					 int pos,
					 int tree_type);
}

#endif // CC_INTERFACE_MOGUL_HH
