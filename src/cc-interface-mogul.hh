
#ifndef CC_INTERFACE_MOGUL_HH
#define CC_INTERFACE_MOGUL_HH

#include "analysis/mogul-interface.hh"

#ifdef HAVE_GOOCANVAS

#include "goograph/goograph.hh"

void show_mogul_geometry_dialog(const coot::mogul &m, mmdb::Residue *residue);
void set_null_goograph_pointer(GtkWidget *w);
#endif 

GtkWidget *wrapped_create_mogul_geometry_dialog(const coot::mogul &m, mmdb::Residue *residue);


namespace coot {
   GtkCellRenderer *
   mogul_results_add_cell_renderer(GtkTreeView *tree_view,
					 GtkTreeStore *tree_store,
					 const std::string &column_title,
					 int pos,
				   int tree_type);
   void fill_mogul_bonds_tab(GtkTreeView *mogul_bonds_treeview, GtkWidget *dialog, const coot::mogul &m, mmdb::Residue *r);
   void fill_mogul_angles_tab(GtkTreeView *mogul_angles_treeview, GtkWidget *dialog, const coot::mogul &m, mmdb::Residue *r);
   void fill_mogul_torsions_tab(GtkTreeView *mogul_torsions_treeview, GtkWidget *dialog, const coot::mogul &m, mmdb::Residue *r);
   void on_mogul_bonds_selection_changed(GtkTreeSelection *treeselection,
					 gpointer          user_data);
   void on_mogul_angles_selection_changed(GtkTreeSelection *treeselection,
					  gpointer          user_data);
   void on_mogul_torsions_selection_changed(GtkTreeSelection *treeselection,
					    gpointer          user_data);
   void update_mogul_histogram_dialog(GtkWidget *histogram_dialog,
				      const mogul &m, const std::vector<std::string> &atom_ids,
				      minimol::residue *r, const std::string &altconf);
#ifdef HAVE_GOOCANVAS   
   void mogul_histogram_for_item(goograph *gg, const mogul_item &item,
				 const std::string &x_axis_label, const std::string &title);
   static void mogul_histogram_dialog_close_callback(GtkWidget *button, GtkWidget *dialog);

#endif    

   
}

#endif // CC_INTERFACE_MOGUL_HH
