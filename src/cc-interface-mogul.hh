
#ifndef CC_INTERFACE_MOGUL_HH
#define CC_INTERFACE_MOGUL_HH

#include "mogul-interface.hh"

void show_mogul_geometry_dialog(const coot::mogul &m, CResidue *residue);
GtkWidget *wrapped_create_mogul_geometry_dialog(const coot::mogul &m, CResidue *residue);


namespace coot {
   GtkCellRenderer *
   mogul_results_add_cell_renderer(GtkTreeView *tree_view,
					 GtkTreeStore *tree_store,
					 const std::string &column_title,
					 int pos,
				   int tree_type);
   void fill_mogul_bonds_tab(GtkTreeView *mogul_bonds_treeview, const coot::mogul &m, CResidue *r);
   void fill_mogul_angles_tab(GtkTreeView *mogul_angles_treeview, const coot::mogul &m, CResidue *r);
   void fill_mogul_torsions_tab(GtkTreeView *mogul_torsions_treeview, const coot::mogul &m, CResidue *r);
   void on_mogul_bonds_selection_changed(GtkTreeSelection *treeselection,
					 gpointer          user_data);
   void on_mogul_angles_selection_changed(GtkTreeSelection *treeselection,
					  gpointer          user_data);
   void on_mogul_torsions_selection_changed(GtkTreeSelection *treeselection,
					    gpointer          user_data);

}

#endif // CC_INTERFACE_MOGUL_HH
