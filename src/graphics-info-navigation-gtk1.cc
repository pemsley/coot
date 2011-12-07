/* src/graphics-info-navigation.cc
 * 
 * Copyright 2006 by The University of York
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
 * Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */




#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#ifdef _MSC_VER
#include <windows.h>
#endif
#include <gtk/gtk.h>

#if (GTK_MAJOR_VERSION == 1)

#include <iostream>

#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

// #include "interface.h"

// #include "molecule-class-info.h"

#include "graphics-info.h"

void
graphics_info_t::make_synthetic_select_on_residue_tree_gtk1(GtkWidget *residue_tree, CAtom *atom_p) const {

   std::cout << "synthetic select on next item in residue tree!\n";
   bool count_subtree_flag;
   int n_list = 0;
   int found_index = 0;

   CResidue *residue_p = atom_p->residue;
   coot::model_view_atom_tree_item_info_t *item_data;
      if (residue_p) {

      GList *children = gtk_container_children(GTK_CONTAINER(residue_tree));
      while (children) { 

	 GtkTree *subtree = GTK_TREE(GTK_TREE_ITEM_SUBTREE(children->data));
	 n_list++;

	 if (GTK_WIDGET_VISIBLE(subtree)) { 
	    count_subtree_flag = 1;
	 } else {
	    count_subtree_flag = 0;
	 }
	 
	 GList *sub_children = gtk_container_children(GTK_CONTAINER(subtree));
	 while (sub_children) {

	    if (count_subtree_flag)
	       n_list++;
	    
	    item_data = (coot::model_view_atom_tree_item_info_t *) gtk_object_get_user_data(GTK_OBJECT(sub_children->data));

	    CResidue *item_residue = item_data->residue;

	    if (item_residue == residue_p) {
// 	       std::cout << "found residue in tree! " << sub_children << " "
// 			 << residue_p << " "
// 			 << residue_p->GetChainID() << residue_p->GetSeqNum()
// 			 << std::endl;
	       // Trash (unselect) the old selection:
	       // 
#if !defined(WINDOWS_MINGW) && !defined(_MSC_VER) // should be GKT2 test
	       GList *tree_selection = GTK_TREE_SELECTION(subtree);
#else
	       GList *tree_selection = gtk_tree_selection_get_selected_rows(GTK_TREE_SELECTION(subtree), NULL);
#endif

	       while (tree_selection) {
		  // gtk_signal_emit_by_name(GTK_OBJECT(tree_selection->data), "deselect");
		  // gtk_tree_unselect_item(tree_selection->data);
		  tree_selection = g_list_remove_link(tree_selection, tree_selection);
	       }

	       // Examine the tree's selection for debugging purposes.
	       //
// 	       GList *debug_tree_selection = residue_tree->selection;
// 	       while (debug_tree_selection) {
		  
// 		  debug_item_data = gtk_object_get_user_data(GTK_OBJECT(debug_tree_selection->data));
// 		  coot::model_view_atom_tree_item_info_t *debug_tree_selection;
// 		  debug_tree_selection = g_list_remove_link(debug_tree_selection, debug_tree_selection);
// 	       }
	       
	       gtk_signal_emit_by_name(GTK_OBJECT(sub_children->data), "select");
	       // selection change signal is only emitted by the root tree:
	       //
	       gtk_signal_emit_by_name(GTK_OBJECT(residue_tree), "selection_changed");
	       found_index = n_list;
	       // If we find the residue and its subtree was not
	       // visible, then we should make the subtree visible.
	       //
	       // 
	    }

	    
	    sub_children = g_list_remove_link(sub_children, sub_children);
	 }
	 children = g_list_remove_link(children, children);
      }
      // position gets clamped 0->1
      float position = float(found_index+5)/float(n_list);
      std::cout << found_index << "/" << n_list << " " << position << std::endl;
      // gtk_tree_scroll_vertical(GTK_TREE(residue_tree), GTK_SCROLL_JUMP, position);
      GtkWidget *sc_win = lookup_widget(residue_tree, "go_to_atom_residue_scrolledwindow");

      GtkAdjustment *vadj = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(sc_win));

      gtk_adjustment_set_value(vadj, 8*float(found_index));

      // gtk_list_scroll_vertical(GTK_SCROLLED_WINDOW(sc_win), GTK_SCROLL_JUMP, position);
   }
}


// a static
//
// Recall that the tree is created in c-interface.cc's fill_go_to_atom_window().
void
graphics_info_t::fill_go_to_atom_residue_list_gtk1(GtkWidget *gtktree) {

   graphics_info_t g;

   g.go_to_atom_residue(); // sets values of unset (magic -1) go to
			   // atom residue number.

   // Claus Flensburg reports a crash around here, when reading in a
   // (dodgy?) pdb file at start up on the command line.

   // So let's test that g.go_to_atom_molecule() is sensible before
   // proceeding: Slightly inelegant logic, so that we can make a
   // small patch instead of a big one.
   bool stop_now = 1;
   if (g.go_to_atom_molecule() < n_molecules())
      if (is_valid_model_molecule(g.go_to_atom_molecule()))
	 stop_now = 0;

   if (stop_now) {
      std::cout << "ERROR:: trapped bad go to atom molecule" << std::endl;
      return;
   }
      
   std::vector<coot::model_view_atom_tree_chain_t> residue_chains = 
      molecules[g.go_to_atom_molecule()].model_view_residue_tree_labels();

   int nres = 0;
   for (unsigned int ichain=0; ichain<residue_chains.size(); ichain++) {
      nres += residue_chains[ichain].tree_residue.size();
   }
   // std::cout << "DEBUG:: tree view has " << nres << " residues" << std::endl;

   // std::cout << "DEBUG:: gtktree has name " << gtk_widget_get_name(gtktree) << std::endl;

   // clear the current residue tree
   //
   // gtk_list_clear_items(GTK_TREE(gtktree), 0, -1);  // FIXME
   //
   // we need to remove all items from a subtree and the all parent
   // tree items.
   // 
   // If the parent tree (or subtrees) is not referenced beforehand,
   // deleting all the (sub)trees items will unparent and destroy it
   // (which is what we want).
   //
   // Possibly not, actually, we are passed the tree here (i.e. a tree
   // is not created here).  So we want to delete all the items and
   // subtrees from it and then add a new set.  So I guess that we
   // should ref the parent tree when we create it.
   // 
   gtk_tree_clear_items(GTK_TREE(gtktree), 0, -1);
   
   // and the current atom list
   clear_atom_list(lookup_widget(GTK_WIDGET(gtktree), "go_to_atom_atom_list"));

   // debugging? 
   GtkSignalFunc residue_tree_item_signal_callback =
      GTK_SIGNAL_FUNC(graphics_info_t::residue_tree_view_itemsignal);
      
   int residue_item_count = 0; // we don't want to add more than 2000
			       // or so residues, otherwise the widget
			       // stops working (for some reason).

   // Here is the plan:
   //
   // outer loop:
   //    create an item
   //    append that item to main tree
   //    show item
   //    create a residue_sub_tree
   //    gtk_tree_item_set_subtree(item, residue_sub_tree)
   //    inner loop:
   //       create a sub item
   //       append sub item to sub tree
   //       show subitem
   // 

   for (unsigned int ichain=0; ichain<residue_chains.size(); ichain++) {
      // the chain label item e.g. "A"
      GtkWidget *chain_tree_item =
	 gtk_tree_item_new_with_label(residue_chains[ichain].chain_id.c_str()); 
      gtk_tree_append(GTK_TREE(gtktree), chain_tree_item);
      gtk_widget_show(chain_tree_item);
      
      GtkWidget *residue_sub_tree = gtk_tree_new();  // this is equivalent to subtree
						     // in the example.
      gtk_tree_set_selection_mode (GTK_TREE(residue_sub_tree),
				   GTK_SELECTION_SINGLE);
      gtk_tree_item_set_subtree(GTK_TREE_ITEM(chain_tree_item), residue_sub_tree);

      gtk_signal_connect (GTK_OBJECT(chain_tree_item), "select",
			  GTK_SIGNAL_FUNC(cb_chain_tree_itemsignal), NULL);
//       gtk_signal_connect (GTK_OBJECT(chain_tree_item), "deselect",
// 			  GTK_SIGNAL_FUNC(cb_chain_tree_itemsignal), NULL);
      
      for(unsigned int ires=0; ires<residue_chains[ichain].tree_residue.size(); ires++) { 
	
	 residue_item_count++;

	 GtkWidget *subitem;
	 if (residue_item_count < 20000) { 
      
	    subitem=gtk_tree_item_new_with_label(residue_chains[ichain].tree_residue[ires].button_label.c_str());
	    gtk_tree_append(GTK_TREE(residue_sub_tree), subitem);
      
	    gtk_widget_show(subitem);
	    	    coot::model_view_atom_tree_item_info_t *item_data =
	     	       new coot::model_view_atom_tree_item_info_t(residue_chains[ichain].tree_residue[ires]);

 	    gtk_object_set_user_data(GTK_OBJECT(subitem), item_data);

	    // For double clicks:  
 	    gtk_signal_connect(GTK_OBJECT(subitem),
 			       "button_press_event",
 			       GTK_SIGNAL_FUNC(go_to_atom_residue_tree_signal_handler_event_gtk1),
 			       NULL);

	    // do we get deselected items?
	    // FIXME: this causes a compilation error: const void *
//  	    gtk_signal_connect (GTK_OBJECT(subitem), "deselect",
//  				residue_tree_item_signal_callback, "deselect");
	 }
      }
   }
}



// static
gint
graphics_info_t::go_to_atom_residue_tree_signal_handler_event_gtk1(GtkWidget *widget, 
								   GdkEventButton *event, 
								   gpointer func_data) {
  if (GTK_IS_TREE_ITEM(widget) &&
       (event->type==GDK_2BUTTON_PRESS ||
        event->type==GDK_3BUTTON_PRESS) ) {
//      printf("I feel %s clicked on button %d\n",
// 	    event->type==GDK_2BUTTON_PRESS ? "double" : "triple", 
// 	    event->button); 
     graphics_info_t g;
     g.apply_go_to_atom_from_widget(go_to_atom_window);
  } 
  return FALSE;
}

/* for all the GtkItem:: and GtkTreeItem:: signals */
// static
int
graphics_info_t::cb_chain_tree_itemsignal( GtkWidget *item,
					   GdkEventButton *event, 
					   gpointer func_data) {

   // Note that we can't user event, it's NULL - this is a select
   // callback, above is a button_press_event callback.

   if (GTK_IS_TREE_ITEM(item)) {
      if (GTK_TREE_ITEM(item)->expanded) {
	 GtkLabel *label = GTK_LABEL(GTK_BIN(item)->child);
	 char *label_name;
	 gtk_label_get(label, &label_name);
// 	 std::cout << "collapsing! " << GTK_TREE (item->parent)->level << " "
// 		   << label_name << std::endl;
	 gtk_tree_item_collapse (GTK_TREE_ITEM(item));
      } else {
	 gtk_tree_item_expand (GTK_TREE_ITEM(item));
      }
   }
   return FALSE;
}

// Mainly called by clicking on a residue in the residue tree.
// Callback set in c-iterface.cc's fill_go_to_atom_window - which
// creates the tree in the first place.
// 
// Also called by a function in callback.c and called by the results of the
// signal_emit in make_synthetic_select_on_residue_list.
//
//  A static.
void
graphics_info_t::on_go_to_atom_residue_tree_selection_changed_gtk1 (GtkList *gtktree,
								    gpointer user_data) {

   // std::cout << "selection changed" << std::endl;
   GList *dlist = GTK_TREE(gtktree)->selection;

   if (!dlist) {
      return;
   }

   GtkObject *list_item;
   coot::model_view_atom_tree_item_info_t *item_data;
   graphics_info_t g; // this is a static, don't forget.
   while (dlist) {
      list_item = GTK_OBJECT(dlist->data);
      int level = GTK_TREE(GTK_WIDGET(list_item)->parent)->level;
      if (level == 1) { 
	 item_data = (coot::model_view_atom_tree_item_info_t *) gtk_object_get_user_data(list_item);
// 	 std::cout << item_data->residue << std::endl;
// 	 std::cout << "getting at from molecule:" << g.go_to_atom_molecule() << std::endl;
// 	 std::cout << "residue has n atoms:" << item_data->residue->GetNumberOfAtoms()
// 		   << std::endl;
	 CAtom *at = graphics_info_t::molecules[g.go_to_atom_molecule()].intelligent_this_residue_mmdb_atom(item_data->residue);
	 // std::cout << "residue tree selection changed: got at:" << at << std::endl;
	 if (at) {
	    // We only want to do this if this was a real event, not a
	    // synthetic one, otherwise make_synthetic_select_on_residue_list
	    // gets called by update_widget_go_to_atom_values(), which calls
	    // this function... endless loop.
	    // 

	    // this does simple setting, nothing else
	    g.set_go_to_atom_chain_residue_atom_name(at->GetChainID(),
						     at->GetSeqNum(),
						     at->GetInsCode(),
						     at->name,
						     at->altLoc);
	 
	    g.update_widget_go_to_atom_values(g.go_to_atom_window, at);
	 
	    // now we want the atom list to contain the atoms of the
	    // newly selected residue:
	 
	    GtkWidget *window = lookup_widget(GTK_WIDGET(gtktree),"goto_atom_window");
	    GtkWidget *atom_gtklist = lookup_widget(window, "go_to_atom_atom_list");
	    fill_go_to_atom_atom_list_gtk1(atom_gtklist,
					   g.go_to_atom_molecule(),
					   at->GetChainID(),
					   at->GetSeqNum());
	 }
      } 
      dlist = dlist->next;
   }
}

// The question we have to ask is: is everything that happens here OK
// when called as a result of a synthetic residue select?
// 
void
graphics_info_t::fill_go_to_atom_atom_list_gtk1(GtkWidget *atom_gtklist, int imol, char *chain_id, int seqno) {

   GtkWidget       *label;
   gchar           *string;
   GtkWidget *list_item;
   graphics_info_t g;

   std::vector<coot::model_view_atom_button_info_t> atoms =
      g.molecules[imol].model_view_atom_button_labels(chain_id, seqno);
   coot::model_view_atom_button_info_t *buttons_copy =
      new coot::model_view_atom_button_info_t[atoms.size()];
   for (unsigned int i=0; i<atoms.size(); i++) {
      buttons_copy[i] = atoms[i];
   }

   // first clear out the atoms already there:
   g.clear_atom_list(atom_gtklist);

   
   for(unsigned int iatom=0; iatom<atoms.size(); iatom++) { 
        
      label=gtk_label_new(buttons_copy[iatom].button_label.c_str());
      list_item=gtk_list_item_new();
      gtk_container_add(GTK_CONTAINER(list_item), label);
      gtk_widget_show(label);
      gtk_container_add(GTK_CONTAINER(atom_gtklist), list_item);
      gtk_widget_show(list_item);
      //       gtk_label_get(GTK_LABEL(label), &string); we don't use string.
//       gtk_object_set_data(GTK_OBJECT(list_item),
//                        list_item_data_key,
//                        string);
      gtk_object_set_user_data(GTK_OBJECT(list_item), &buttons_copy[iatom]);

      // For double clicks:
      gtk_signal_connect(GTK_OBJECT(list_item),
			 "button_press_event",
			 GTK_SIGNAL_FUNC(go_to_atom_atom_list_signal_handler_event_gtk1),
			 NULL);
   }
}


// static
gint
graphics_info_t::go_to_atom_atom_list_signal_handler_event_gtk1(GtkWidget *widget, 
							   GdkEventButton *event, 
							   gpointer func_data) {
  if (GTK_IS_LIST_ITEM(widget) &&
       (event->type==GDK_2BUTTON_PRESS ||
        event->type==GDK_3BUTTON_PRESS) ) {
//      printf("I feel %s clicked on button %d\n",
// 	    event->type==GDK_2BUTTON_PRESS ? "double" : "triple", 
// 	    event->button); 
     graphics_info_t g;
     g.apply_go_to_atom_from_widget(go_to_atom_window);
  } 
  return FALSE;
} 

#endif // #if (GTK_MAJOR_VERSION == 1)
