/* src/main.cc
 * 
 * Copyright 2010 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#include "mmdb_manager.h"
#include <gtk/gtk.h>

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh"

#include "pisa-interface.hh"

void
coot::pisa_interfaces_gui(const std::vector<coot::pisa_interface_t> &gui_info) {

#if (GTK_MAJOR_VERSION > 1)
   if (graphics_info_t::use_graphics_interface_flag) {
      GtkWidget *w = create_pisa_interfaces_dialog();
      gtk_widget_show(w);
      GtkWidget *treeview = lookup_widget(w, "pisa_interfaces_treeview");
      // we show:
      // molecule-number molecule-number chain-id chain-id symop Area energy #H-bonds #salt-bridges #Nss #Ncov
      GtkTreeStore *tree_store_interfaces =
	 gtk_tree_store_new (11, G_TYPE_INT, G_TYPE_INT, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_INT, G_TYPE_INT, G_TYPE_INT, G_TYPE_INT);
      GtkTreeView *tv_interfaces = GTK_TREE_VIEW(treeview);
      GtkTreeIter   toplevel;
      gtk_tree_view_set_model(tv_interfaces, GTK_TREE_MODEL(tree_store_interfaces));
      for (unsigned int i=0; i<gui_info.size(); i++) {
	 gtk_tree_store_append(tree_store_interfaces, &toplevel, NULL);
	 gtk_tree_store_set(tree_store_interfaces, &toplevel,
			    0, gui_info[i].imol_1,
			    1, gui_info[i].imol_2,
			    2, gui_info[i].chain_id_1.c_str(),
			    3, gui_info[i].chain_id_2.c_str(),
			    4, gui_info[i].symop.c_str(),
			    5, gui_info[i].bsa,
			    6, gui_info[i].solv_en,
			    7, gui_info[i].n_h_bonds,
			    8, gui_info[i].n_salt_bridges,
			    9, gui_info[i].n_ss_bonds,
			    10, gui_info[i].n_cov_bonds,
			    -1);
      }
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "Mol No", 0);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "Mol No", 1);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "ChainID", 2);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "ChainID", 3);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "Symop",   4);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "Buried Surface Area", 5);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "Solvation Energy",    6);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "#-H-bonds",           7);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "#Salt-bridges",       8);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "#SS-bonds",           9);
      add_pisa_interfaces_cell_renderer(tv_interfaces, tree_store_interfaces, "#Covalent-bonds",    10);
   }
   
#endif   
}

#if (GTK_MAJOR_VERSION > 1)
void
coot::add_pisa_interfaces_cell_renderer(GtkTreeView *tree_view,
					GtkTreeStore *tree_store,
					const std::string &column_title,
					int pos) {

   GtkCellRenderer *cell_renderer = gtk_cell_renderer_text_new();
   GtkTreeViewColumn *column = gtk_tree_view_column_new();
   gtk_tree_view_column_set_title(column, column_title.c_str());
   gtk_tree_view_append_column(tree_view, column);
   gtk_tree_view_column_pack_start(column, cell_renderer, TRUE);
   gtk_tree_view_column_add_attribute(column, cell_renderer, "text", pos);
   g_object_set(cell_renderer, "editable", FALSE, NULL);
   
   gtk_tree_view_column_set_sort_column_id(column, pos);
   
   char *s = new char[column_title.length() + 1];
   strcpy(s, column_title.c_str());
   g_object_set_data (G_OBJECT (cell_renderer), "column", GINT_TO_POINTER (pos));

}
#endif
