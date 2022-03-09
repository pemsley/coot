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

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#include <string>
#include <mmdb2/mmdb_manager.h>
#include <gtk/gtk.h>

#include "graphics-info.h"
#include "interface.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh"

#include "pisa-interface.hh"
#include "widget-from-builder.hh"

void
coot::pisa_interfaces_gui(const std::vector<coot::pisa_interface_t> &gui_info) {

   if (graphics_info_t::use_graphics_interface_flag) {
      // GtkWidget *w = create_pisa_interfaces_dialog();
      GtkWidget *w = widget_from_builder("pisa_interface_dialog");
      gtk_widget_show(w);
      // GtkWidget *treeview = lookup_widget(w, "pisa_interfaces_treeview");
      GtkWidget *treeview = widget_from_builder("pisa_interfaces_treeview");
      // we show:
      // molecule-number molecule-number chain-id chain-id symop Area energy #H-bonds #salt-bridges #Nss #Ncov
      GtkTreeStore *tree_store_interfaces =
	 gtk_tree_store_new (12, G_TYPE_INT, G_TYPE_INT, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING,
			     G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_FLOAT, G_TYPE_INT, G_TYPE_INT, G_TYPE_INT,
			     G_TYPE_INT);
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
			    5, gui_info[i].interface_area,
			    6, gui_info[i].interface_solv_en,
			    7, gui_info[i].interface_pvalue,
			    8, gui_info[i].n_h_bonds,
			    9, gui_info[i].n_salt_bridges,
			    10, gui_info[i].n_ss_bonds,
			    11, gui_info[i].n_cov_bonds,
			    -1);
      }
      add_pisa_interfaces_cell_renderer(tv_interfaces, "Mol No", 0);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "Mol No", 1);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "ChainID", 2);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "ChainID", 3);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "Symop",   4);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "Interf. Area", 5);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "Solv. Energy", 6);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "p-value",      7);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "#H-bond",      8);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "#Salt-br.",    9);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "#SS",         10);
      add_pisa_interfaces_cell_renderer(tv_interfaces, "#Cov",        11);

      GtkTreeSelection *sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));

      // Here we malloc but don't ever give it back
      // 
      std::vector<coot::pisa_interface_t> *gui_info_copy = new std::vector<coot::pisa_interface_t>;
      *gui_info_copy = gui_info;
      
      g_signal_connect(sel, "changed", (GCallback) coot::on_pisa_interfaces_seletion_changed,
		       gui_info_copy);
   }
   
}

void
coot::add_pisa_interfaces_cell_renderer(GtkTreeView *tree_view,
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
   g_object_set_data (G_OBJECT (cell_renderer), "column", GINT_TO_POINTER (pos));
}



void
coot::on_pisa_interfaces_seletion_changed(GtkTreeSelection *treeselection,
					 gpointer          user_data) {

   std::vector<coot::pisa_interface_t> *gui_info_p
      = static_cast<std::vector<coot::pisa_interface_t> *> (user_data);

   GtkTreeIter  iter;
   GtkTreeModel *model; 

   // This?
   // gtk_tree_selection_iter_is_selected (treeselection, iter);

   // set iter and model.
   // 
   gboolean r = gtk_tree_selection_get_selected(treeselection, &model, &iter);

   if (r) { 
      gchar *chain_id_1, *chain_id_2;
      int imol_1;
      int imol_2;
      
      gtk_tree_model_get (model, &iter,
			  0, &imol_1,
			  1, &imol_2,
			  2, &chain_id_1,
			  3, &chain_id_2,
			  -1);

      for (unsigned int i_interf=0; i_interf<gui_info_p->size(); i_interf++) {
	 // std::cout << i_interf << " " << (*gui_info_p)[i_interf].imol_1 << std::endl;
	 if ((*gui_info_p)[i_interf].imol_1 == imol_1) { 
	    if ((*gui_info_p)[i_interf].imol_2 == imol_2) {
	       if ((*gui_info_p)[i_interf].chain_id_1 == std::string(chain_id_1)) {
		  if ((*gui_info_p)[i_interf].chain_id_2 == std::string(chain_id_2)) {
		     clipper::Coord_orth centre_pt = (*gui_info_p)[i_interf].centre;
		     pisa_interfaces_display_only(imol_1, imol_2, centre_pt);
		  }
	       }
	    }
	 }
      } 
      
      g_free(chain_id_1);
      g_free(chain_id_2);

   }
}
