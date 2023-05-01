/* src/c-interface.cc
 * 
 * Copyright 2007 The University of York
 * Copyright 2014 by Medical Research Council
 * written by Paul Emsley
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include "compat/coot-sysdep.h"

#include <gtk/gtk.h>

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

#include "graphics-info.h"

#include "widget-from-builder.hh"

enum {
   CHAIN_COL,
   RESIDUE_COL,
   VISIBILITY_COL
};

void graphics_info_t::fill_go_to_atom_window_gtk3(GtkWidget *dialog) {

   auto get_model_molecule_vector = [] () {
                                       graphics_info_t g;
                                       std::vector<int> vec;
                                       int n_mol = g.n_molecules();
                                       for (int i=0; i<n_mol; i++)
                                          if (g.is_valid_model_molecule(i))
                                             vec.push_back(i);
                                       return vec;
                                    };

   // make this a wrapper function for a graphics_info_t function
   // After the GTK3 build is working - FIXME

   graphics_info_t g;
   int gimol = g.go_to_atom_molecule();

   gchar *text;

   GCallback callback_func = G_CALLBACK(go_to_atom_mol_combobox_changed);

   // GtkWidget *combobox = lookup_widget(widget, "go_to_atom_molecule_combobox");
   GtkWidget *combobox = widget_from_builder("go_to_atom_molecule_combobox");

   // fill_combobox_with_coordinates_options(combobox, callback_func, gimol);

   // use this instead
   // void
   // graphics_info_t::fill_combobox_with_molecule_options(GtkWidget *combobox,
   //                                                  GCallback signal_func,
   //                                                  int imol_active_position,
   //                                                  const std::vector<int> &molecules_index_vec)

   // maybe make a function to wrap this?
   //
   auto mols = get_model_molecule_vector();
   g.fill_combobox_with_molecule_options(combobox, callback_func, gimol, mols); // not static!

   /* These are in a special order: The residue is done first
      because it is set to a magic number (-9999 (or so)) initially.
      In that case, we do magic in
      get_text_for_go_to_atom_residue_entry(), i.e. look up a real
      atom of a molecule and set also the go to chain and the go to
      atom name  */

   /* The residue entry */

   // residue_entry = lookup_widget(GTK_WIDGET(widget), "go_to_atom_residue_entry");
   GtkWidget *residue_entry = widget_from_builder("go_to_atom_residue_entry");

   // text = get_text_for_go_to_atom_residue_entry();  // old

   // on startup, tinkers with
   // go to atom params, yuck,
   // I think.
   std::string rt = coot::util::int_to_string(g.go_to_atom_residue());

   gtk_editable_set_text(GTK_EDITABLE(residue_entry), rt.c_str());

   /* Now that the go to atom molecule has been set, we can use it
      to fill the molecule option menu */

   // Needs new
   // fill_option_menu_with_coordinates_options(option_menu,
   // callback_func,
   // gimol);

   /* The chain entry */

   // GtkWidget *chain_entry = lookup_widget(GTK_WIDGET(widget), "go_to_atom_chain_entry");
   GtkWidget *chain_entry = widget_from_builder("go_to_atom_chain_entry");

   gtk_editable_set_text(GTK_EDITABLE(chain_entry), g.go_to_atom_chain());


   /* The Atom Name entry */

   // atom_name_entry = lookup_widget(GTK_WIDGET(widget), "go_to_atom_atom_name_entry");
   GtkWidget *atom_name_entry = widget_from_builder("go_to_atom_atom_name_entry");

   gtk_editable_set_text(GTK_EDITABLE(atom_name_entry), g.go_to_atom_atom_name());

   /* The Residue List */

   /* The residue list cant be added to a scrolled window in glade,
      so we create only a scrolled window in glade
      (go_to_atom_residue_scrolledwindow) and add the list to it
      like is done in examples/list/list.c */

   // GtkWidget *residue_scrolled_window = lookup_widget(GTK_WIDGET(widget), "go_to_atom_residue_scrolledwindow");
   // GtkWidget *atom_list_scrolled_window = lookup_widget(GTK_WIDGET(widget), "go_to_atom_atom_scrolledwindow");

   // GtkWidget *residue_scrolled_window = widget_from_builder("go_to_atom_residue_scrolledwindow");
   // GtkWidget *atom_list_scrolled_window = widget_from_builder("go_to_atom_atom_scrolledwindow");

   // g.fill_go_to_atom_window_gtk2(dialog, // the go to atom window
   // residue_scrolled_window,
   // atom_list_scrolled_window);

   // std::cout << "--------------------------- fill_go_to_atom_window_gtk3() calling "
   //              << "fill_go_to_atom_window_residue_and_atom_lists_gtk4() " << std::endl;
   fill_go_to_atom_window_residue_and_atom_lists_gtk4();

   /* store the widget */
   go_to_atom_window = dialog;

}

// static
void
graphics_info_t::fill_go_to_atom_window_residue_and_atom_lists_gtk4() {

   auto add_residue_tree_to_scrolled_window_if_needed = [] (GtkWidget *scrolled_window) {

                                                GtkTreeView *residue_tree = nullptr; // return this

                                                GtkWidget *item_widget = gtk_widget_get_first_child(scrolled_window);
                                                while (item_widget) {
                                                   if (GTK_IS_TREE_VIEW(item_widget)) {
                                                      residue_tree = GTK_TREE_VIEW(item_widget);
                                                   }
                                                   item_widget = gtk_widget_get_next_sibling(item_widget);
                                                }
                                                if (! residue_tree) {
                                                   residue_tree = GTK_TREE_VIEW(gtk_tree_view_new());
                                                   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(residue_tree));
                                                }
                                                return residue_tree;
                                             };


      auto add_atom_list_to_scrolled_window_if_needed = [] (GtkWidget *scrolled_window) {

         GtkTreeView *atom_list = nullptr; // return this
         GtkWidget *item_widget = gtk_widget_get_first_child(scrolled_window);
         while (item_widget) {
            if (GTK_IS_TREE_VIEW(item_widget)) {
               atom_list = GTK_TREE_VIEW(item_widget);
            }
            item_widget = gtk_widget_get_next_sibling(item_widget);
         };
         if (! atom_list) {
            atom_list = GTK_TREE_VIEW(gtk_tree_view_new());
            gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(atom_list));
         }
         return atom_list;
      };


   // --------------------- main line -------------------------------------

   GtkWidget *residue_tree_scrolled_window = widget_from_builder("go_to_atom_residue_scrolledwindow");
   GtkWidget *atom_list_scrolled_window    = widget_from_builder("go_to_atom_atom_scrolledwindow");

   if (!   residue_tree_scrolled_window) return;
   if (! atom_list_scrolled_window) return;

   graphics_info_t g;
   g.go_to_atom_residue();  // 20220812-PE this should be const.

   GtkTreeView *atom_list = add_atom_list_to_scrolled_window_if_needed(atom_list_scrolled_window);

   GtkTreeView *residue_tree = add_residue_tree_to_scrolled_window_if_needed(residue_tree_scrolled_window);
   GtkTreeView *gtktree = residue_tree;

   gtk_widget_set_size_request(GTK_WIDGET(gtktree), 200, -1);

   int imol = go_to_atom_molecule_; // function to get this is not static.
   if (!is_valid_model_molecule(imol)) return;

   // so, clear the current tree:
   GtkTreeView *tv = NULL;
   if (gtktree)
      tv = GTK_TREE_VIEW(gtktree);
   if (! tv)
      tv = GTK_TREE_VIEW(gtk_tree_view_new());

   GtkTreeModel *model = gtk_tree_view_get_model(tv);
   bool need_renderer = true;
   if (model) {
      // std::cout << "clearing old tree store" << std::endl;
      gtk_tree_store_clear(GTK_TREE_STORE(model));
      need_renderer = false;
   }

   bool show_waters_flag = false;
   GtkWidget *show_waters_check_button = widget_from_builder("go_to_atom_show_waters_togglebutton");
   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(show_waters_check_button)))
      show_waters_flag = true;
   std::vector<coot::model_view_atom_tree_chain_t> residue_chains = molecules[imol].model_view_residue_tree_labels(show_waters_flag);

   // 20230402-PE actually, I won't use visibility. I will get a new list of residues
   // GtkTreeStore *tree_store = gtk_tree_store_new (3, G_TYPE_STRING, G_TYPE_POINTER, G_TYPE_BOOLEAN);
   GtkTreeStore *tree_store = gtk_tree_store_new (2, G_TYPE_STRING, G_TYPE_POINTER);
   GtkTreeIter   toplevel, child;

   // what is the connection between tree_store and model?
   gtk_tree_view_set_model(GTK_TREE_VIEW(gtktree), GTK_TREE_MODEL(tree_store));

   for (unsigned int ichain=0; ichain<residue_chains.size(); ichain++) {
      gtk_tree_store_append(GTK_TREE_STORE(tree_store), &toplevel, NULL);

      //       std::cout << "Adding tree item " << residue_chains[ichain].chain_id
      //                 << std::endl;
      gtk_tree_store_set (tree_store, &toplevel,
                          CHAIN_COL, residue_chains[ichain].chain_id.c_str(),
                          RESIDUE_COL, NULL,
                          -1);

      for (unsigned int ires=0; ires<residue_chains[ichain].tree_residue.size(); ires++) {

         // memory leak.  We don't give this memory back.
         coot::residue_spec_t *rsp = new coot::residue_spec_t(residue_chains[ichain].tree_residue[ires].residue_spec);
         gpointer res_spec_ptr = static_cast<gpointer> (rsp);
         gtk_tree_store_append(GTK_TREE_STORE(tree_store), &child, &toplevel);
         if (false)
            std::cout << "tree store set for " << ires << " " << residue_chains[ichain].tree_residue[ires].residue_spec
                      << " " << residue_chains[ichain].tree_residue[ires].button_label << std::endl;
         std::string lab =  residue_chains[ichain].tree_residue[ires].button_label;
         gtk_tree_store_set(tree_store, &child,
                            CHAIN_COL, lab.c_str(),
                            RESIDUE_COL, res_spec_ptr,
                            -1);
      }
   }

   if (need_renderer) {  // need_renderer
      GtkCellRenderer *cell_renderer = gtk_cell_renderer_text_new();
      // GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes ("Chains", cell, "text", 0, NULL);
      GtkTreeViewColumn *column = gtk_tree_view_column_new_with_attributes ("Chains", cell_renderer, "markup", 0, NULL);


      gtk_tree_view_append_column(GTK_TREE_VIEW (tv), GTK_TREE_VIEW_COLUMN (column));

      GtkTreeSelection*   tree_sel = gtk_tree_view_get_selection (tv);
      gtk_tree_selection_set_mode(tree_sel, GTK_SELECTION_SINGLE);
      // double clicks
      g_signal_connect(tv, "row-activated", G_CALLBACK(residue_tree_residue_row_activated), NULL);

      gtk_tree_selection_set_select_function(tree_sel, residue_tree_selection_func, NULL, NULL);

   }

   // Now handle the atom list - clear it.
   //
   if (atom_list) {
      GtkTreeModel *model = gtk_tree_view_get_model(GTK_TREE_VIEW(atom_list));
      if (model) {
         GtkListStore *list_store = GTK_LIST_STORE(model);
         gtk_list_store_clear(list_store);
      }
   }


}

void goto_residue_show_hide_water_items() {

   // can this be simplified? (copied from above)
   auto add_residue_tree_to_scrolled_window_if_needed = [] (GtkWidget *scrolled_window) {

                                                GtkTreeView *residue_tree = nullptr; // return this

                                                GtkWidget *item_widget = gtk_widget_get_first_child(scrolled_window);
                                                while (item_widget) {
                                                   if (GTK_IS_TREE_VIEW(item_widget)) {
                                                      residue_tree = GTK_TREE_VIEW(item_widget);
                                                   }
                                                   item_widget = gtk_widget_get_next_sibling(item_widget);
                                                }
                                                if (! residue_tree) {
                                                   residue_tree = GTK_TREE_VIEW(gtk_tree_view_new());
                                                   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(residue_tree));
                                                }
                                                return residue_tree;
                                             };

#if 0
   // http://scentric.net/tutorial/sec-misc-get-renderer-from-click.html

   GtkTreeViewColumn *col = NULL;
	GList             *columns, *cells;
	cells = gtk_tree_view_column_get_cell_renderers(col);
   GtkCellRenderer *cell_renderer = cells[0];


   GtkTreeIter iter;
   GtkTreePath *path; // needs to be set

   GtkWidget *residue_tree_scrolled_window = widget_from_builder("go_to_atom_residue_scrolledwindow");
   GtkTreeView *residue_tree = add_residue_tree_to_scrolled_window_if_needed(residue_tree_scrolled_window);
   GtkTreeView *gtktree = residue_tree;
   GtkTreeView *tv = NULL;
   if (gtktree)
      tv = GTK_TREE_VIEW(gtktree);
   if (! tv)
      tv = GTK_TREE_VIEW(gtk_tree_view_new());

   GtkTreeModel *model = gtk_tree_view_get_model(tv);

   gtk_tree_model_get_iter(model, &iter, path); //  convert the path to an iter.
   int bol;
   gtk_tree_model_get(model, &iter, VISIBILITY_COL, &bol, -1);
   if (bol == 0) {
      g_object_set(cell_renderer, "visible", FALSE, NULL);
   } else {
      g_object_set(cell_renderer, "visible", TRUE, "active", TRUE, NULL);
   }
#endif
}


// This gets called on window close with the residue_tree as the data.
//
// static
void
graphics_info_t::go_to_atom_residue_tree_destroy(gpointer data) {

   // std::cout << "go_to_atom_residue_tree_destroy called " << data << std::endl;
   // the residue_tree is a GtkTreeView (that's what gtk_tree_view_new() returns).
   GtkTreeView *residue_tree = static_cast<GtkTreeView *>(data);
   // gtk_widget_unref(GTK_WIDGET(residue_tree));
}


// This gets called on window close with the atom_list as the data.
//
// static
void
graphics_info_t::go_to_atom_list_destroy(gpointer data) {

   // std::cout << "go_to_atom_list_destroy called " << data << std::endl;
   GtkTreeView *atom_list = static_cast<GtkTreeView *>(data);
}


void
residue_button_info_free(coot::model_view_residue_button_info_t *bi) {
   delete bi;
}

coot::model_view_residue_button_info_t*
residue_button_info_copy(coot::model_view_residue_button_info_t *bi) {
   coot::model_view_residue_button_info_t *n = new coot::model_view_residue_button_info_t(*bi);
   return n;
}


// static
gboolean
graphics_info_t::residue_tree_selection_func(GtkTreeSelection *selection,
					     GtkTreeModel *model,
					     GtkTreePath *path,
					     gboolean path_currently_selected,
					     gpointer data) {

   GtkTreeIter   iter;
   gboolean can_change_selected_status_flag = TRUE;

    if (gtk_tree_model_get_iter(model, &iter, path)) {
       gchar *name;
       gtk_tree_model_get(model, &iter, CHAIN_COL, &name, -1);
       if (!path_currently_selected) {
	  if (true) {  // if this was a residue, not a chain click
	               // update the go to atom residues from the
   	               // characteristics of this row... bleurgh.. how do I do
	               // that!?  Check the level, is how I did it in gtk1...
	     graphics_info_t g;
	     int go_to_imol = g.go_to_atom_molecule();
	     if (is_valid_model_molecule(go_to_imol)) {
		gpointer residue_data = 0;
		gtk_tree_model_get(model, &iter, RESIDUE_COL, &residue_data, -1);
		if (!residue_data) {
		   // This was a "Outer" chain row click (there was no residue)
		   // std::cout << "Null residue " << residue_data << std::endl;
		   // std::cout << "should expand here, perhaps?" << std::endl;
		} else {
		   // mmdb::Residue *res = (mmdb::Residue *) residue_data;
		   
		   coot::residue_spec_t *rsp = static_cast<coot::residue_spec_t *> (residue_data);
		   mmdb::Residue *res = molecules[go_to_imol].get_residue(*rsp);
		   mmdb::Atom *at = molecules[go_to_imol].intelligent_this_residue_mmdb_atom(res);
		   if (!at) {
		      std::cout << "ERROR:: failed to get atom in intelligent_this_residue_mmdb_atom: "
				<< go_to_imol << " " << res << " (tree selected)" << std::endl;
		   } else { 
		      // this does simple setting, nothing else
		      g.set_go_to_atom_chain_residue_atom_name(at->GetChainID(),
							       at->GetSeqNum(),
							       at->GetInsCode(),
							       at->name,
							       at->altLoc);
		      
		      g.update_widget_go_to_atom_values(g.go_to_atom_window, at);
		      
		      // now we want the atom list to contain the atoms of the
		      // newly selected residue:

		      // Fill me...
		      fill_go_to_atom_atom_list_gtk2(g.go_to_atom_window,
						     g.go_to_atom_molecule(),
						     at->GetChainID(),
						     at->GetSeqNum(),
						     at->GetInsCode());
		   }
		}
	     }
	  }
       }
       g_free(name);
    }
    return can_change_selected_status_flag;
}


// static
void
graphics_info_t::residue_tree_residue_row_activated(GtkTreeView        *treeview,
						    GtkTreePath        *path,
						    GtkTreeViewColumn  *col,
						    gpointer            userdata) {

   // a double-click (I think) so that we do a "Go to atom"

   // This gets called on double-clicking, and not on single clicking
   
   GtkTreeModel *model = gtk_tree_view_get_model(treeview);
   GtkTreeIter   iter;

    if (gtk_tree_model_get_iter(model, &iter, path)) {
       gchar *name;
       gtk_tree_model_get(model, &iter, CHAIN_COL, &name, -1);
       // g_print ("Double-clicked row contains name %s\n", name);
       if (true) {
	  graphics_info_t g;
	  int go_to_imol = g.go_to_atom_molecule();
	  if (is_valid_model_molecule(go_to_imol)) {
	     gpointer residue_data;
	     gtk_tree_model_get(model, &iter, RESIDUE_COL, &residue_data, -1);
	     if (!residue_data) {
		// This was a "Outer" chain row click (there was no residue)
		// std::cout << "Null residue " << residue_data << std::endl;
		// std::cout << "should expand here, perhaps?" << std::endl;
	     } else {

		// mmdb::Residue *res = (mmdb::Residue *) residue_data;

		coot::residue_spec_t *rsp = static_cast<coot::residue_spec_t *> (residue_data);

		mmdb::Residue *res = molecules[go_to_imol].get_residue(*rsp);
		
		mmdb::Atom *at = molecules[go_to_imol].intelligent_this_residue_mmdb_atom(res);
		// this does simple setting, nothing else
		if (!at) {
		   std::cout << "ERROR:: failed to get atom in intelligent_this_residue_mmdb_atom: "
			     << go_to_imol << " " << res << " (row_activated)" << std::endl;
		} else {
                   // this only sets variables
		   g.set_go_to_atom_chain_residue_atom_name(at->GetChainID(),
							    at->GetSeqNum(),
							    at->GetInsCode(),
							    at->name,
							    at->altLoc);
		   
		   g.update_widget_go_to_atom_values(g.go_to_atom_window, at);
		   g.apply_go_to_atom_from_widget(go_to_atom_window);
		}
	     }
	  }
       }
       g_free(name);
    }
}


// 20220813-PE this is in the execution path. Change the function name.
// static
void
graphics_info_t::fill_go_to_atom_atom_list_gtk2(GtkWidget *go_to_atom_window, int imol,
                                                char *chain_id, int resno, char *ins_code) {

   // GtkTreeView *atom_tree = GTK_TREE_VIEW(lookup_widget(go_to_atom_window, "go_to_atom_atom_list"));
   // GtkTreeView *atom_tree = GTK_TREE_VIEW(widget_from_builder("go_to_atom_atom_list"));
   GtkWidget *scrolled_window = widget_from_builder("go_to_atom_atom_scrolledwindow");

   GtkTreeView *atom_tree = nullptr;
   // if there is something alreay in scrolled_window, let's use that as the atom tree.
   // If not, make a new one and add it into scrolled_window.

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME container children

   // set the atom_tree

   GtkWidget *item_widget = gtk_widget_get_first_child(scrolled_window);
   while (item_widget) {
      if (GTK_IS_TREE_VIEW(item_widget)) {
         atom_tree = GTK_TREE_VIEW(item_widget);
      }
      item_widget = gtk_widget_get_next_sibling(item_widget);
   };

#else
   // keep this for a bit while testing
   GList *dlist = gtk_container_get_children(GTK_CONTAINER(scrolled_window));
   GList *free_list = dlist;
   while (dlist) {
      GtkWidget *list_item = GTK_WIDGET(dlist->data);
      if (GTK_IS_TREE_VIEW(list_item)) {
         atom_tree = GTK_TREE_VIEW(list_item);
      }
      dlist = dlist->next;
   }
   g_list_free(free_list);
#endif

   bool need_renderer = true;
   if (!atom_tree) {
      // std::cout << "making new atom_tree..." << std::endl;
      // add a new one
      atom_tree = GTK_TREE_VIEW(gtk_tree_view_new());
      // GtkWidget *scrolled_window = lookup_widget(GTK_WIDGET(go_to_atom_window), "go_to_atom_atom_scrolledwindow");
      // gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
      //                                       GTK_WIDGET(atom_tree));
      gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(scrolled_window), GTK_WIDGET(atom_tree));
   }

   std::string ins_code_str(ins_code);
   std::string chain_id_str(chain_id);
   std::vector<coot::model_view_atom_button_info_t> atoms =
      graphics_info_t::molecules[imol].model_view_atom_button_labels(chain_id_str, resno, ins_code_str);

   GtkListStore *list_store = 0;
   GtkTreeModel *model = gtk_tree_view_get_model(atom_tree);
   if (!model) {
      list_store = gtk_list_store_new (2, G_TYPE_STRING, G_TYPE_POINTER);
   } else {
      // clear (any) old atom tree/list
      list_store = GTK_LIST_STORE(model);
      gtk_list_store_clear(list_store);
      need_renderer = 0;
   }

   GtkTreeIter   toplevel;
   gtk_tree_view_set_model(GTK_TREE_VIEW(atom_tree), GTK_TREE_MODEL(list_store));

   for(unsigned int iatom=0; iatom<atoms.size(); iatom++) {
      // std::cout << " storing atom at: " << atoms[iatom].atom << std::endl;
      gtk_list_store_append(GTK_LIST_STORE(list_store), &toplevel);
      gtk_list_store_set(GTK_LIST_STORE(list_store), &toplevel,
                         CHAIN_COL, atoms[iatom].button_label.c_str(),
                         RESIDUE_COL, (gpointer) atoms[iatom].atom,
                         -1);
   }

   if (need_renderer) {
      GtkCellRenderer *cell = gtk_cell_renderer_text_new();
      GtkTreeViewColumn *column =
         gtk_tree_view_column_new_with_attributes ("Atoms", cell, "text", 0, NULL);
      gtk_tree_view_append_column (GTK_TREE_VIEW (atom_tree), GTK_TREE_VIEW_COLUMN (column));

      GtkTreeSelection*   tree_sel = gtk_tree_view_get_selection (atom_tree);
      gtk_tree_selection_set_mode(tree_sel, GTK_SELECTION_SINGLE);
      // double clicks
      g_signal_connect(atom_tree, "row-activated", G_CALLBACK(atom_tree_atom_row_activated), NULL);
      gtk_tree_selection_set_select_function (tree_sel,
                                              graphics_info_t::atom_tree_selection_func,
                                              NULL, NULL);
   }

}

// static
void
graphics_info_t::atom_tree_atom_row_activated(GtkTreeView        *treeview,
					      GtkTreePath        *path,
					      GtkTreeViewColumn  *col,
					      gpointer            userdata) {

   // This gets called on double-clicking, and not on single clicking

   GtkTreeModel *model = gtk_tree_view_get_model(treeview);
   GtkTreeIter   iter;

   if (gtk_tree_model_get_iter(model, &iter, path)) {
      gchar *name;
      gtk_tree_model_get(model, &iter, CHAIN_COL, &name, -1);
      // g_print("Double clicked row contains name: %s\n", name);
      if (1) {
	 graphics_info_t g;
	 int go_to_imol = g.go_to_atom_molecule();
	 if (go_to_imol< n_molecules()) {
	    gpointer atom_data;
	    gtk_tree_model_get(model, &iter, RESIDUE_COL, &atom_data, -1);
	    if (! atom_data) {
	       std::cout << "ERROR:: no atom data!" << std::endl;
	    } else {
	       mmdb::Atom *at = (mmdb::Atom *) atom_data;
	       // std::cout << " reading from atom at: " << at << std::endl;
	       g.set_go_to_atom_chain_residue_atom_name(at->GetChainID(),
							at->GetSeqNum(),
							at->GetInsCode(),
							at->name,
							at->altLoc);
	       g.update_widget_go_to_atom_values(g.go_to_atom_window, at);
	       g.apply_go_to_atom_from_widget(go_to_atom_window);
	    }
	 }
      }
      g_free(name);
   }
}


// static
gboolean
graphics_info_t::atom_tree_selection_func(GtkTreeSelection *selection,
					  GtkTreeModel *model,
					  GtkTreePath *path,
					  gboolean path_currently_selected,
					  gpointer data) {

   // This gets called on single-clicking

   gboolean can_change_selected_status_flag = TRUE;
   GtkTreeIter   iter;
   
   if (gtk_tree_model_get_iter(model, &iter, path)) {
      gchar *name;
      gtk_tree_model_get(model, &iter, CHAIN_COL, &name, -1);
      // g_print("Double clicked row contains name: %s\n", name);
      if (!path_currently_selected) {
	 graphics_info_t g;
	 int go_to_imol = g.go_to_atom_molecule();
	 if (go_to_imol< n_molecules()) {
	    gpointer atom_data;
	    gtk_tree_model_get(model, &iter, RESIDUE_COL, &atom_data, -1);
	    if (! atom_data) {
	       std::cout << "ERROR:: no atom data!" << std::endl;
	    } else {
	       mmdb::Atom *at = (mmdb::Atom *) atom_data;
	       // std::cout << " reading from atom at: " << at << std::endl;
	       g.set_go_to_atom_chain_residue_atom_name(at->GetChainID(),
							at->GetSeqNum(),
							at->GetInsCode(),
							at->name,
							at->altLoc);

	       g.update_widget_go_to_atom_values(g.go_to_atom_window, at);
	    }
	 }
      }
      g_free(name);
   }
   return can_change_selected_status_flag;
}
