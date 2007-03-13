/* src/graphics-info-navigation.cc
 * 
 * Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
 * Copyright 2005 by Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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


#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
#define GTK_ENABLE_BROKEN
#if defined _MSC_VER
#include <windows.h>
#endif
#endif


#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>

#include <iostream>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

// #include "interface.h"

// #include "molecule-class-info.h"

#include "graphics-info.h"


// After this function, we call try_centre_from_new_go_to_atom(); 
void graphics_info_t::set_go_to_atom_chain_residue_atom_name(const gchar *t1, 
							     int it2, 
							     const gchar *t3){
   // these are strings, don't forget.
   go_to_atom_chain_     = t1; 
   go_to_atom_residue_   = it2; 
   go_to_atom_atom_name_ = t3; 
   go_to_atom_atom_altLoc_ = "empty";  // Reset it to something
                                       // sensible don't keep the old
                                       // altLoc; see
                                       // find_atom_index_from_goto_info()

   go_to_atom_inscode_ = "";

}

// After this function, we call try_centre_from_new_go_to_atom(); 
void 
graphics_info_t::set_go_to_atom_chain_residue_atom_name(const char *t1, 
							int it2, const char *t3, const char *altLoc) { 
   
   // these are strings, don't forget.
   go_to_atom_chain_     = t1; 
   go_to_atom_residue_   = it2; 
   go_to_atom_atom_name_ = t3; 
   go_to_atom_atom_altLoc_ = altLoc;
   go_to_atom_inscode_ = "";
   
} 


// After this function, we call try_centre_from_new_go_to_atom(); 
void 
graphics_info_t::set_go_to_atom_chain_residue_atom_name(const char *chain_id, 
							int resno, const char *ins_code,
							const char *atom_name, const char *altLoc) {
   // these are strings, don't forget.
   go_to_atom_chain_     = chain_id; 
   go_to_atom_residue_   = resno; 
   go_to_atom_atom_name_ = atom_name; 
   go_to_atom_atom_altLoc_ = altLoc;
   go_to_atom_inscode_ = ins_code;
}


const char *graphics_info_t::go_to_atom_chain() {

   return go_to_atom_chain_.c_str(); 
}

// Return the current go_to_atom residue number, if special case
// (-9999) (not set) the set other values too.
//
int graphics_info_t::go_to_atom_residue() {

   if (go_to_atom_residue_ == -9999) { // magic number
      for (int imol=0; imol<n_molecules; imol++) {
	 if (molecules[imol].atom_sel.n_selected_atoms > 0) {

	    if (molecules[imol].drawit) {

	       go_to_atom_residue_ = molecules[imol].atom_sel.atom_selection[0]->GetSeqNum();
	       go_to_atom_chain_   = std::string(molecules[imol].atom_sel.atom_selection[0]->GetChainID());
	       go_to_atom_molecule_ = imol;

	       // and set the atom name by intelligent atom:
	       //
	       CResidue *res = molecules[imol].atom_sel.atom_selection[0]->residue;
	       CAtom *atom = molecules[imol].intelligent_this_residue_mmdb_atom(res);
	       go_to_atom_atom_name_ = std::string(atom->name);
	       break;
	    }
	 }
      }
   }
   return go_to_atom_residue_;
}

int graphics_info_t::go_to_atom_molecule() {

   return go_to_atom_molecule_; 
}

const char *graphics_info_t::go_to_atom_atom_name() {
   return go_to_atom_atom_name_.c_str(); 
}

const char *
graphics_info_t::go_to_atom_ins_code() {
   return go_to_atom_inscode_.c_str();
}

const char *
graphics_info_t::go_to_atom_alt_conf() {
   return go_to_atom_atom_altLoc_.c_str();
}

void 
graphics_info_t::set_go_to_atom_molecule(int pos) { 

   // 
   go_to_atom_molecule_ = pos;

}

int
graphics_info_t::try_centre_from_new_go_to_atom() { 

   // 
   // std::string atom_string = make_mmdb_atom_string_from_go_to_atom(); 

   // find_atom_index returns pick_info using prestored atom
   // information (go_to_atom_* variables).
   //
   int imol = go_to_atom_molecule();
   pick_info pi = find_atom_index_from_goto_info(imol); 
   if (pi.success) {

      setRotationCentre(pi.atom_index, go_to_atom_molecule()); 

   } else { 
      cout << "Sorry atom " << go_to_atom_atom_name() 
	   << " \"" << go_to_atom_atom_altLoc_ << "\""
	   << "/" << go_to_atom_residue()
	   << "\"" << go_to_atom_inscode_ << "\"" 
	   << "/" << go_to_atom_chain()
	   << " not found in molecule " << go_to_atom_molecule() << endl;
   }
   return pi.success; 
}

// Return atom_name, altconf: e.g. (" CA ", "B") from " CA ,B").
// 
// static
std::pair<std::string, std::string> 
graphics_info_t::split_atom_name(const std::string &atom_name) { 

   std::pair<std::string, std::string> v("","");

   std::string::size_type icomma = atom_name.find_last_of(",");
   if (icomma == string::npos) { 
      // no comma 
      v.first = atom_name;
   } else {
      v.first  = atom_name.substr(0, icomma);
      unsigned int an_length = atom_name.length();
      if (an_length > (icomma + 1)) {
	 v.second = atom_name.substr(icomma + 1, an_length);
      }
   }
   
   return v;
}

// static
//
// return a pair, the first of which is the resno as a string, the
// second of which is the inscode as a string.
// 
std::pair<std::string, std::string> 
graphics_info_t::split_resno_inscode(const std::string &entry_str) {

   std::pair<std::string, std::string> v("","");
   char char_0 = '0';
   char char_9 = '9';

   v.first = entry_str; // without inscode this is correct.

   // This test excludes residues with an altconf of"-".  So unlikely
   // that we can ignore them.

   for (int i=entry_str.length()-1; i>=0; i--) {
      // is it non-numeric?
      if (((entry_str[i] < char_0) ||
	   (entry_str[i] > char_9)) &&
	  (entry_str[i] != '-')) // test for residues with negative resno
	 {
	 // we ignore spaces
	 if (entry_str[i] != ' ') {
	    v.second = entry_str.substr(i, i+1);
	    if (i>0) {
	       v.first  = entry_str.substr(0, i);
	    }
	 }
      } 
   }
//    std::cout << "DEBUG:: :" << v.first << ": :"
// 	     << v.second << ":" << std::endl;
   return v;
}


// Return success status:
// It is fine to call this will null go_to_atom_window.
//
int
graphics_info_t::intelligent_next_atom_centring(GtkWidget *go_to_atom_window) {

   return intelligent_near_atom_centring(go_to_atom_window, std::string("next"));

}

// Return success status:
// It is fine to call this will null go_to_atom_window.
//
int
graphics_info_t::intelligent_previous_atom_centring(GtkWidget *go_to_atom_window) {

   return intelligent_near_atom_centring(go_to_atom_window, std::string("previous"));

}

// direction is either "next" or "previous"
// 
int
graphics_info_t::intelligent_near_atom_centring(GtkWidget *go_to_atom_window,
						const std::string &direction) {

   std::string chain =     go_to_atom_chain_;
   std::string atom_name = go_to_atom_atom_name_;
   std::string ins_code =  go_to_atom_inscode_;
   int resno = go_to_atom_residue();
   int imol = go_to_atom_molecule();

   // check how the spaces in the name work in the
   // find_atom_index_from_goto_info function

   // we have have to do something more clever... but whatever it
   // is, shouldn't it be a member function of
   // molecule_class_info_t?  Yes, it should.


   // OK then, return an atom index, -1 on failure.
   //
   int atom_index = -1;

   if (molecules[imol].atom_sel.mol == 0) {

      std::cout << "ERROR:: bad go to atom molecule (" << imol
		<< ") in intelligent_near_atom_centring" << std::endl;

   } else {

      if (direction == "next") { 
	 atom_index = molecules[imol].intelligent_next_atom(chain, resno, atom_name, ins_code);
      } else { // "previous"
	 atom_index = molecules[imol].intelligent_previous_atom(chain, resno, atom_name, ins_code);
      } 

      if (atom_index != -1) {
	 CAtom *next_atom = molecules[imol].atom_sel.atom_selection[atom_index];

	 go_to_atom_chain_       = next_atom->GetChainID();
	 go_to_atom_atom_name_   = next_atom->name;
	 go_to_atom_residue_     = next_atom->GetSeqNum();
	 go_to_atom_atom_altLoc_ = next_atom->altLoc;

	 // now update the widget with the new values of the above (like
	 // c-interface:goto_near_atom_maybe())

	 if (go_to_atom_window) { 
	    update_widget_go_to_atom_values(go_to_atom_window, next_atom);
	    // 	 GtkWidget *residue_tree = lookup_widget(go_to_atom_window, 
	    // 						 "go_to_atom_residue_tree"); 
	    // make_synthetic_select_on_residue_tree(residue_tree, next_atom);
	 } 
	 try_centre_from_new_go_to_atom();
	 
	 // Update the graphics (glarea widget):
	 // 
	 update_things_on_move_and_redraw(); // (symmetry, environment, map) and draw it
      }
   }
   return 1;
}

// This is a function like:
// graphics_info_t::set_go_to_atom_chain_residue_atom_name(const char
// *t1, int it2, const char *t3, const char *altLoc)
// in that it sets the go_to_atom variables.
// (try_centre_from_new_go_to_atom should be called after this function)
// 
// We will set the atom name to " CA " if there is one, if not, then
// first atom in the residue.
//
// Note the non-optimalness: we find the atom index when we find the
// atom in the residue, but then throw it away again and use
// try_centre_from_new_go_to_atom()
//
void
graphics_info_t::set_go_to_residue_intelligent(const std::string &chain_id, int resno,
					       const std::string &ins_code) { 

      
   // OK, we want a molecule function that returns an atom name
   // (either " CA ", or the first atom in the residue or "no-residue"
   // (error-flag)).

   CAtom *at = molecules[go_to_atom_molecule()].atom_intelligent(chain_id, resno, ins_code);

   if (at) { 
      go_to_atom_chain_ = chain_id;
      go_to_atom_residue_ = resno;
      go_to_atom_atom_name_ = std::string(at->name);
      go_to_atom_atom_altLoc_ = std::string(at->altLoc);
      go_to_atom_inscode_ = ins_code;
   } else { 
      std::cout << "Sorry - can't find residue " << resno << " " << chain_id 
		<< " in molecule " << go_to_atom_molecule() << std::endl;
   } 
} 


// We pass atom so that we know the mmdb atom and residue which we
// have just centred on.  Using that, we can also update the residue
// and atom list (like Jan wants) rather than just the entry widgets.
// 
void
graphics_info_t::update_widget_go_to_atom_values(GtkWidget *window, CAtom *atom)  {

   std::string res_str   = int_to_string(go_to_atom_residue_);
   res_str += go_to_atom_inscode_; 

   GtkEntry *entry;

   if (window) { 
      entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_chain_entry"));
      gtk_entry_set_text(entry, go_to_atom_chain_.c_str());

      entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_residue_entry"));
      gtk_entry_set_text(entry,res_str.c_str());

      entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_atom_name_entry"));
      std::string atom_name_txt = go_to_atom_atom_name_;
      if (! (go_to_atom_atom_altLoc_ == "empty")) {
	 if (go_to_atom_atom_altLoc_ != "") { 
	    atom_name_txt += ",";
	    atom_name_txt += go_to_atom_atom_altLoc_;
	 }
      }
      gtk_entry_set_text(entry, atom_name_txt.c_str());

   } else {
      std::cout << "ERROR: Null window in update_widget_go_to_atom_values\n";
   } 
}



// We need to run through the list items looking for a item that
// has an attached user data that is the residue of that atom:
// 
// When we find it, we generate a synthetic signal that that list_item
// has been selected - which runs the callback
// on_go_to_atom_residue_list_selection_changed().
//
// Ah, but now we have a tree not a list, so we have to look into the
// tree, rather than the top list.  We do something similar elsewhere,
// I think... Oh yes, we check the tree level somewhere...
// on_go_to_atom_residue_tree_selection_changed().
//
// We need to go into the list and check if this item has a tree
// attached (it should do) and check *those* items.
// 
void
graphics_info_t::make_synthetic_select_on_residue_tree(GtkWidget *residue_tree, CAtom *atom_p) const {

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

   make_synthetic_select_on_residue_tree_gtk1(residue_tree, atom_p);

#else

#endif   

}



// Return pick_info using prestored atom information
// (go_to_atom_* variables).
//
pick_info
graphics_info_t::find_atom_index_from_goto_info(int imol) { 

   pick_info pi; 

   pi.success = GL_FALSE; 

   // actually, try ignoring the atom_string argument and consider
   // using mmdb atom selection

   if (imol < 0) { 
      std::cout << "WARNING:: no molecule for imol = " << imol << std::endl;
      
   } else { 

      if (imol >= n_molecules) { 
	 std::cout << "WARNING:: no molecule for imol = " << imol << std::endl;
      } else { 
      
	 if (graphics_info_t::molecules[imol].atom_sel.mol == NULL ) {
	    std::cout << "WARNING: (Programmer error) looking for atoms "
		      << "in a molecule " << imol << " which is null" << std::endl;
	 } else {
	    atom_selection_container_t AtomSel = 
	       graphics_info_t::molecules[imol].atom_sel; 

	    int selHnd = AtomSel.mol->NewSelection(); 

	    // 0 -> any (mmdb) model
	    // 
	    // note that we have to do ugly (char *) casting.  Hopefully 
	    // that will go away in new version of mmdb..? (20 Aug 2002 - PE)
	    // 
	    std::pair<std::string, std::string> p = 
	       graphics_info_t::split_atom_name(go_to_atom_atom_name());

	    char altconf[2];
	    strncpy(altconf, go_to_atom_atom_altLoc_.c_str(), 2);
	    if (go_to_atom_atom_altLoc_ == "empty") { 
	       strcpy(altconf, "");
	    }

// 	    std::cout << "FAI:: searching chain :" << go_to_atom_chain() << std::endl;
// 	    std::cout << "FAI:: searching residue no :" << go_to_atom_residue() << std::endl;
// 	    std::cout << "FAI:: searching inscode :" << go_to_atom_inscode_ << std::endl;
// 	    std::cout << "FAI:: searching atom_name :" << p.first << std::endl;
// 	    std::cout << "FAI:: searching altconf   :" << altconf << ":" << std::endl;

	    AtomSel.mol->SelectAtoms (selHnd, 0, (char *) go_to_atom_chain(),
				      go_to_atom_residue(), // starting resno, an int
				      (char *) go_to_atom_inscode_.c_str(), // any insertion code
				      go_to_atom_residue(), // ending resno
				      (char *) go_to_atom_inscode_.c_str(), // ending insertion code
				      "*", // any residue name
				      (char *) p.first.c_str(), // atom name
				      "*", // elements
				      "*"
				      ); 

	    int nSelAtoms; 
	    PPCAtom local_SelAtom = NULL;

	    // modify nSelAtoms
	    // 
	    AtomSel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms); 

	    if (nSelAtoms > 0) {

	       // So, we have just done a selection of all alt confs in this residue.
	       // We want to return the success with the index of the atom with the
	       // correct alt conf, and failing that, we return with the index of the
	       // first atom in the selection.

	       pi.imol = go_to_atom_molecule(); // not used?

	       // So the question is however, what is the index of the atom in
	       // the AtomSel.atom_selection array?
	       // 
	       // Let's use the UDD (c.f. full_atom_spec_to_atom_index)

	       int ic;
	       std::string target_altconf = go_to_atom_atom_altLoc_;
	       for (int iat=0; iat<nSelAtoms; iat++) {
		  if (target_altconf == local_SelAtom[iat]->altLoc) { 
		     if (local_SelAtom[iat]->GetUDData(AtomSel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
			pi.success = GL_TRUE; 
			pi.atom_index = ic;
// 			std::cout << "DEBUG:: matched altconf gives udd index "
// 				  << pi.atom_index << std::endl;
			break;
		     } else {
	 
			// Fall back to comparing pointers.
			// 
			for (int i=0; i<AtomSel.n_selected_atoms; i++) { 
			   if (AtomSel.atom_selection[i] == local_SelAtom[0]) { 
			      pi.success = GL_TRUE; 
			      pi.atom_index = i;
// 			      std::cout << "DEBUG:: pointer comparison 1 hit "
// 					<< pi.atom_index << std::endl;
			      break; 
			   }
			}
		     }
		  }
	       }

	       // The altconf didn't match, so try an atom that is the same except for the altconf:
	       // 
	       if (pi.success != GL_TRUE) {
		  if (local_SelAtom[0]->GetUDData(AtomSel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
		     pi.success = GL_TRUE; 
		     pi.atom_index = ic;
// 		     std::cout << "DEBUG:: non-matched altconf gives udd index "
// 			       << pi.atom_index << std::endl;
// 		     std::cout << "Atom at that index: "
// 			       << molecules[imol].atom_sel.atom_selection[pi.atom_index] << std::endl;
		  } else {
		     
		     // Fall back to comparing pointers.
		     // 
		     for (int i=0; i<AtomSel.n_selected_atoms; i++) { 
			if (AtomSel.atom_selection[i] == local_SelAtom[0]) { 
			   pi.success = GL_TRUE; 
			   pi.atom_index = i;
// 			   std::cout << "DEBUG:: pointer comparison 2 hit "
// 				     << pi.atom_index << std::endl;
			   break; 
			}
		     }
		  }
	       }
	       if (pi.success == GL_TRUE)
		  go_to_atom_atom_altLoc_ = AtomSel.atom_selection[pi.atom_index]->altLoc;

	       // should we update the go to atom widget's atom name entry now? (it
	       // is not updated elsewhere)
		  
	    } else {
	       // this happens when the search atom is not in the molecule

	       int selHnd_check = AtomSel.mol->NewSelection(); 
	       AtomSel.mol->SelectAtoms (selHnd_check, 0, "*",
					 ANY_RES, // starting resno, an int
					 "*", // any insertion code
					 ANY_RES, // ending resno
					 "*", // ending insertion code
					 "*", // any residue name
					 "*", // atom name
					 "*", // elements
					 "*"  // alt loc.
					 );
	       int nSelAtoms_check;
	       PPCAtom local_SelAtom_check;
	       AtomSel.mol->GetSelIndex(selHnd_check,
					local_SelAtom_check, nSelAtoms_check);

	       std::cout << "There are " << nSelAtoms_check << " atoms "
			 << "in the molecule " << std::endl; 
	       // 
	       std::cout << "find_atom_index_from_goto_info(), "
			 << "no matching atoms of the " << AtomSel.n_selected_atoms
			 <<" in this (non-null) molecule ("
			 << imol << ")" << std::endl;
	       AtomSel.mol->DeleteSelection(selHnd_check);
	    }
	    AtomSel.mol->DeleteSelection(selHnd);
	 }
      }
   }
   return pi; 
}

// imol has changed.
// Now fix up the Go_To_Atom window to match:
// 
void
graphics_info_t::update_go_to_atom_window_on_changed_mol(int imol) {

   // now if the go to atom widget was being displayed, we need to
   // redraw the residue list and atom list (if the molecule of the
   // residue and atom list is the molecule that has just been
   // deleted)

   if (go_to_atom_window) { 

      // The go to atom molecule matched this molecule, so we
      // need to regenerate the residue and atom lists.
      GtkWidget *gtktree = lookup_widget(go_to_atom_window,
					 "go_to_atom_residue_tree");
      if (gtktree == NULL) {
	 std::cout << "ERROR:: gtktree (go_to_atom_residue_tree) is null!\n"; 
      } else {
	 graphics_info_t::fill_go_to_atom_residue_list(gtktree);
      }
   } 
}


// a static
//
// Recall that the tree is created in c-interface.cc's fill_go_to_atom_window().
void
graphics_info_t::fill_go_to_atom_residue_list(GtkWidget *gtktree) {

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

   graphics_info_t::fill_go_to_atom_residue_list_gtk1(gtktree);

#else

#endif   

}


#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

/* for all the GtkItem:: and GtkTreeItem:: signals */
// static
void graphics_info_t::residue_tree_view_itemsignal( GtkWidget *item,
						    gchar     *signame )
{
   gchar *name;
   GtkLabel *label;
   
   /* It's a Bin, so it has one child, which we know to be a
      label, so get that */
   label = GTK_LABEL (GTK_BIN (item)->child);
   /* Get the text of the label */
   gtk_label_get (label, &name);
   /* Get the level of the tree which the item is in */
   g_print ("%s called for item %s->%p, level %d\n", signame, name,
	    item, GTK_TREE (item->parent)->level);
}

#endif // #if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)


// static
gint
graphics_info_t::go_to_atom_residue_list_signal_handler_event(GtkWidget *widget, 
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




void 
graphics_info_t::apply_go_to_atom_from_widget(GtkWidget *widget) {

  const gchar *chain_str;
  const gchar *res_str; 
  const gchar *atom_name_str; 

  GtkEntry *entry; 

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(widget), 
				  "go_to_atom_chain_entry"));
  chain_str = gtk_entry_get_text(entry);

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(widget), 
				  "go_to_atom_residue_entry"));
  res_str = gtk_entry_get_text(entry);

  entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(widget), 
				  "go_to_atom_atom_name_entry"));

  const gchar *txt =  gtk_entry_get_text(entry);
  if (txt) {
     std::pair<std::string, std::string> p = 
	graphics_info_t::split_atom_name(std::string(txt));
//      std::cout << "DEBUG: split: " << std::string(txt) << " into :" 
// 	       << p.first << ":  :" << p.second << ":\n" ;

     // we have to use the version of set_go_to.. that has 4 params,
     // because the 3 parameter version sets the altconf to "empty:.
     // 
     atom_name_str = p.first.c_str();

     std::pair<std::string, std::string> resno_inscode =
	split_resno_inscode(std::string(res_str));
     int resno = atoi(resno_inscode.first.c_str());
     std::string inscode = resno_inscode.second;

//      std::cout << "DEBUG:: in apply_go_to_atom_from_widget: chain_str " << chain_str
// 	       << " resno " << resno << " inscode " << inscode << " atom_name_str "
// 	       << atom_name_str << " split thing " << p.second << "\n";

     set_go_to_atom_chain_residue_atom_name(chain_str,
					    resno,
					    inscode.c_str(),
					    atom_name_str, 
					    p.second.c_str());
     int success = try_centre_from_new_go_to_atom();
     if (success) 
	update_things_on_move_and_redraw();
  }
}


// The question we have to ask is: is everything that happens here OK
// when called as a result of a synthetic residue select?
// 
void
graphics_info_t::fill_go_to_atom_atom_list(GtkWidget *atom_gtklist, int imol, char *chain_id, int seqno) {

   GtkWidget       *label;
   gchar           *string;
   GtkWidget *list_item;
   std::string button_string;
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
			 GTK_SIGNAL_FUNC(go_to_atom_atom_list_signal_handler_event),
			 NULL);
   }
}


// static
gint
graphics_info_t::go_to_atom_atom_list_signal_handler_event(GtkWidget *widget, 
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

// Coordinates have been read into a new molecule (imol).
// 
// Now fix up the Go_To_Atom window to match by changing the option menu
// 
void
graphics_info_t::update_go_to_atom_window_on_new_mol() {

   // std::cout << "DEBUG:: update_go_to_atom_window_on_new_mol called" << std::endl;

   if (go_to_atom_window) {
      GtkWidget *option_menu =
	 lookup_widget(GTK_WIDGET(go_to_atom_window), 
		       "go_to_atom_molecule_optionmenu");
      
      GtkSignalFunc callback_func = 
	 GTK_SIGNAL_FUNC(graphics_info_t::go_to_atom_mol_menu_item_select);
      // set last active (1)
      fill_option_menu_with_coordinates_options_internal(option_menu, callback_func, 1);

      // If there was no molecule already, we need to update the atom
      // lists too.
      int nmol = 0;
      int mol_no= -1;
      for (int imol=0; imol<n_molecules; imol++) {
	 if (molecules[imol].has_model()) {
	    nmol++;
	    mol_no = imol;
	 }
      }
      // if (nmol == 1)
      // if (nmol)
      update_go_to_atom_window_on_changed_mol(mol_no);
   }
}

// return -1 on error
//
int
graphics_info_t::go_to_atom_molecule_optionmenu_active_molecule(GtkWidget *widget) { 

   return go_to_atom_molecule();
}

// a static
void
graphics_info_t::clear_atom_list(GtkWidget *atom_gtklist) {

   gtk_list_clear_items(GTK_LIST(atom_gtklist), 0, -1);

}

// void
// graphics_info_t::fill_go_to_atom_option_menu(GtkWidget *option_menu) {

//    std::cout << "This function should not be called!" << std::endl;

// //    GtkSignalFunc callback_func = 
// //       GTK_SIGNAL_FUNC(graphics_info_t::go_to_atom_mol_menu_item_select);
   
// //    fill_option_menu_with_coordinates_options(option_menu, callback_func);

   
// }


// a static
void
graphics_info_t::go_to_atom_mol_menu_item_select(GtkWidget *item, GtkPositionType pos) { 

   std::cout << "INFO:: (menu item select) Go To Atom molecule now: " << pos << std::endl;
   graphics_info_t g;
   int old_go_to_molecule = g.go_to_atom_molecule();
   g.set_go_to_atom_molecule(pos);
   if (pos != old_go_to_molecule) {
      // old style:
//       GtkWidget *residue_gtklist = lookup_widget(GTK_WIDGET(item),
// 						 "go_to_atom_residue_list");
      GtkWidget *residue_gtktree = lookup_widget(GTK_WIDGET(item),
						 "go_to_atom_residue_tree");
      fill_go_to_atom_residue_list(residue_gtktree);
   }
}

void
graphics_info_t::on_go_to_atom_residue_tree_selection_changed(GtkList *gtktree,
							      gpointer user_data) {

#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)

   graphics_info_t::on_go_to_atom_residue_tree_selection_changed_gtk1(gtktree, user_data);

#endif    

}

void
graphics_info_t::undo_last_move() {  // suggested by Frank von Delft


   coot::Cartesian c(old_rotation_centre_x,
		     old_rotation_centre_y,
		     old_rotation_centre_z);
   
   std::cout << "INFO:: Moving back to old centre: " << c << std::endl;
   setRotationCentre(c);
   for(int ii=0; ii<n_molecules; ii++) {
      molecules[ii].update_map();
      molecules[ii].update_symmetry();
   }
   graphics_draw();
} 


// do it if have intermediate atoms and ctrl is pressed.
// 
// axis: 0 for Z, 1 for X.
// 
short int
graphics_info_t::rotate_intermediate_atoms_maybe(short int axis, double angle) {

   short int handled_flag = 0;

   if (rot_trans_rotation_origin_atom) { 
      if (moving_atoms_asc) {
	 if (moving_atoms_asc->n_selected_atoms > 0) {
	    if (control_is_pressed) {
	       if (axis == 0)
		  rotate_intermediate_atoms_round_screen_z(angle);
	       else 
		  rotate_intermediate_atoms_round_screen_x(angle);
	       handled_flag = 1;
	    }
	 }
      }
   }

   return handled_flag;
}
