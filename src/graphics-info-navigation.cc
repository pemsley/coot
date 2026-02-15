
/* src/graphics-info-navigation.cc
 *
 * Copyright 2004, 2005, 2006 by The University of York
 * Copyright 2013, 2016 by Medical Research Council
 * Author Paul Emsley, Bernhard Lohkamp
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

#include "coot-utils/coot-coord-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#if defined(_MSC_VER)
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

#include <string.h> // strncpy

#include <gtk/gtk.h>

#include <iostream>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"

// #include "interface.h"

// #include "molecule-class-info.h"

#include "graphics-info.h"

#include "cc-interface.hh"  // for statusbar_text

#include "widget-from-builder.hh"

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

// 20230611-PE this is indeed a gui callback
// static
void
graphics_info_t::set_go_to_atom(int imol, const coot::atom_spec_t &spec) {

   graphics_info_t g; // crazy
   g.set_go_to_atom_molecule(imol);;
   go_to_atom_chain_      = spec.chain_id;
   go_to_atom_atom_name_  = spec.atom_name;
   go_to_atom_residue_    = spec.res_no;
   go_to_atom_atom_altLoc_= spec.alt_conf;

   // now user should call try_centre_from_new_go_to_atom()

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
graphics_info_t::set_go_to_atom_chain_residue_atom_name(const std::string &chain_id,
							int resno, const std::string &ins_code,
							const std::string &atom_name, const std::string &altLoc) {
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
// (-9999) (not set) then set other go to atom values too.
//
int graphics_info_t::go_to_atom_residue() {

   if (go_to_atom_residue_ == -9999) { // magic number
      for (int imol=0; imol<n_molecules(); imol++) {
	 if (molecules[imol].atom_sel.n_selected_atoms > 0) {

	    if (molecules[imol].draw_it) {

	       go_to_atom_residue_ = molecules[imol].atom_sel.atom_selection[0]->GetSeqNum();
	       go_to_atom_chain_   = std::string(molecules[imol].atom_sel.atom_selection[0]->GetChainID());
	       if (! is_valid_model_molecule(go_to_atom_molecule()))
		  set_go_to_atom_molecule(imol);

	       // and set the atom name by intelligent atom:
	       //
	       mmdb::Residue *res = molecules[imol].atom_sel.atom_selection[0]->residue;
	       mmdb::Atom *atom = molecules[imol].intelligent_this_residue_mmdb_atom(res);
	       go_to_atom_atom_name_ = std::string(atom->name);
	       break;
	    }
	 }
      }
   }
   return go_to_atom_residue_;
}

mmdb::Atom *
graphics_info_t::get_atom(int imol, const coot::atom_spec_t &spec) const {

   mmdb::Atom *at = nullptr;
   if (is_valid_model_molecule(imol)) {
      at = coot::util::get_atom(spec, molecules[imol].atom_sel.mol);
   }
   return at;
}

// static
mmdb::Residue *
graphics_info_t::get_residue(int imol, const coot::residue_spec_t &spec) {

   mmdb::Residue *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      r = molecules[imol].get_residue(spec);
   }
   // find the residue in the intermediate atoms?
   if (imol == -1) {
      if (moving_atoms_asc) {
         if (moving_atoms_asc->mol) {
            r = coot::util::get_residue(spec, moving_atoms_asc->mol);
         }
      }
   }
   return r;
}

void
graphics_info_t::go_to_residue(int imol, const coot::residue_spec_t &spec) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = get_residue(imol, spec);
      if (residue_p) {
         mmdb::Atom *atom = coot::util::intelligent_this_residue_mmdb_atom(residue_p);
         if (atom) {
            clipper::Coord_orth pt = coot::co(atom);
            set_rotation_centre(pt);
         }
      }
   }
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
      std::cout << "WARNING:: atom with name \"" << go_to_atom_atom_name()
                << "\" alt-loc \"" << go_to_atom_atom_altLoc_ << "\","
                << " res-no: " << go_to_atom_residue()
                << ", ins-code \"" << go_to_atom_inscode_ << "\"," 
                << " chain: \"" << go_to_atom_chain()
                << "\" not found in molecule " << go_to_atom_molecule() << std::endl;
      std::string w = "WARNING:: atom ";
      w += go_to_atom_atom_name();
      w += go_to_atom_atom_altLoc_;
      w += " ";
      w += coot::util::int_to_string(go_to_atom_residue());
      w += " ";
      w += go_to_atom_chain();
      w += " not found in molecule ";
      w += coot::util::int_to_string(go_to_atom_molecule());
      add_status_bar_text(w);
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
   if (icomma == std::string::npos) {
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
graphics_info_t::intelligent_next_atom_centring() {

   return intelligent_near_atom_centring(std::string("next"));

}

// Return success status:
//
int
graphics_info_t::intelligent_previous_atom_centring() {

   return intelligent_near_atom_centring(std::string("previous"));

}

// direction is either "next" or "previous"
//
int
graphics_info_t::intelligent_near_atom_centring(const std::string &direction) {

   auto label_unlabel_centre_atom = [] (mmdb::Atom *at_next, mmdb::Manager *mol, int imol) {
      coot::residue_spec_t residue_next(at_next->GetResidue());
      mmdb::Residue *residue_this = coot::util::get_previous_residue(residue_next, mol);
      if (residue_this) {
         molecules[imol].add_atom_label(residue_next.chain_id.c_str(), residue_next.res_no, " CA ");
         molecules[imol].remove_atom_label(residue_this->GetChainID(), residue_this->GetSeqNum(), " CA ");
      }
   };

   std::string chain =     go_to_atom_chain_;
   std::string atom_name = go_to_atom_atom_name_;
   std::string ins_code =  go_to_atom_inscode_;
   int resno = go_to_atom_residue();
   int imol  = go_to_atom_molecule();

   if (false) {
      std::cout << "intelligent_near_atom_centring() " << direction << std::endl;
      std::cout << "intelligent_near_atom_centring() " << imol << std::endl;
      std::cout << "intelligent_near_atom_centring() :" << chain << ":" << std::endl;
      std::cout << "intelligent_near_atom_centring() " << resno << std::endl;
      std::cout << "intelligent_near_atom_centring() :" << ins_code << ":" << std::endl;
      std::cout << "intelligent_near_atom_centring() :" << atom_name << ":" << std::endl;
   }

   // check how the spaces in the name work in the
   // find_atom_index_from_goto_info function

   // we have have to do something more clever... but whatever it
   // is, shouldn't it be a member function of
   // molecule_class_info_t?  Yes, it should.


   // OK then, return an atom index, -1 on failure.
   //
   int atom_index = -1;

   // Pressing space at the start causes (caused) a crash.  Fixed here. 20080518
   if (! is_valid_model_molecule(imol))
      return atom_index;

   if (molecules[imol].atom_sel.mol == 0) {

      std::cout << "ERROR:: bad go to atom molecule (" << imol
                << ") in intelligent_near_atom_centring" << std::endl;

   } else {

      coot::Cartesian rc = RotationCentre();
      if (direction == "next") {
         atom_index = molecules[imol].intelligent_next_atom(chain, resno, atom_name, ins_code, rc);
      } else { // "previous"
         atom_index = molecules[imol].intelligent_previous_atom(chain, resno, atom_name, ins_code, rc);
      }

      if (atom_index != -1) {
         if (atom_index >= molecules[imol].atom_sel.n_selected_atoms) {
            std::cout << "ERROR:: in intelligent_near_atom_centring() atom index error! " 
                      << atom_index << " " << molecules[imol].atom_sel.n_selected_atoms << std::endl;
            return 0; // failure
         }
         mmdb::Atom *next_atom = molecules[imol].atom_sel.atom_selection[atom_index];
         mmdb::Manager *mol    = molecules[imol].atom_sel.mol;

         if (next_atom) {

            go_to_atom_chain_       = next_atom->GetChainID();
            go_to_atom_atom_name_   = next_atom->name;
            go_to_atom_residue_     = next_atom->GetSeqNum();
            go_to_atom_inscode_     = next_atom->GetInsCode();
            go_to_atom_atom_altLoc_ = next_atom->altLoc;

            // now update the widget with the new values of the above (like
            // c-interface:goto_near_atom_maybe())

            label_unlabel_centre_atom(next_atom, mol, imol);

            GtkWidget *go_to_atom_window = widget_from_builder("goto_atom_window");
            if (go_to_atom_window) {
               update_widget_go_to_atom_values(go_to_atom_window, next_atom);
               //          GtkWidget *residue_tree = lookup_widget(go_to_atom_window,
               //                                                  "go_to_atom_residue_tree");
               // make_synthetic_select_on_residue_tree(residue_tree, next_atom);
            }
            try_centre_from_new_go_to_atom();

            // Update the graphics (glarea widget):
            //
            // update_things_on_move_and_redraw(); // (symmetry, environment, map) and draw it
            // and show something in the statusbar
            std::string ai;
            ai = atom_info_as_text_for_statusbar(atom_index, imol);
            add_status_bar_text(ai);

            std::cout << "if sequence view is displayed update highlighted position here C " << std::endl;
         }
      }
   }
   return 1;
}

void
graphics_info_t::add_picked_atom_info_to_status_bar(int imol, int atom_index) {

   std::string ai;
   ai = atom_info_as_text_for_statusbar(atom_index, imol);
   add_status_bar_text(ai);
}


// This is a function like:
// graphics_info_t::set_go_to_atom_chain_residue_atom_name(const char
// *t1, int it2, const char *t3, const char *altLoc)
// in that it sets the go_to_atom variables.
// (try_centre_from_new_go_to_atom should be called after this function)
//
// We will set the atom name to " CA " if there is one, if not, then
// first atom in the residue.
//';m,.'
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

   mmdb::Atom *at = molecules[go_to_atom_molecule()].atom_intelligent(chain_id, resno, ins_code);

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
graphics_info_t::update_widget_go_to_atom_values(GtkWidget *window, mmdb::Atom *atom)  {

   std::string res_str   = int_to_string(go_to_atom_residue_);
   res_str += go_to_atom_inscode_;

   GtkEntry *entry;

   if (window) {
      // entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_chain_entry"));
      entry = GTK_ENTRY(widget_from_builder("go_to_atom_chain_entry"));
      gtk_editable_set_text(GTK_EDITABLE(entry), go_to_atom_chain_.c_str());

      // entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_residue_entry"));
      entry = GTK_ENTRY(widget_from_builder("go_to_atom_residue_entry"));
      gtk_editable_set_text(GTK_EDITABLE(entry), res_str.c_str());

      // entry = GTK_ENTRY(lookup_widget(window, "go_to_atom_atom_name_entry"));
      entry = GTK_ENTRY(widget_from_builder("go_to_atom_atom_name_entry"));
      std::string atom_name_txt = go_to_atom_atom_name_;
      if (! (go_to_atom_atom_altLoc_ == "empty")) {
	 if (go_to_atom_atom_altLoc_ != "") {
	    atom_name_txt += ",";
	    atom_name_txt += go_to_atom_atom_altLoc_;
	 }
      }
      gtk_editable_set_text(GTK_EDITABLE(entry), atom_name_txt.c_str());

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
graphics_info_t::make_synthetic_select_on_residue_tree(GtkWidget *residue_tree, mmdb::Atom *atom_p) const {

#if (GTK_MAJOR_VERSION == 1)

   make_synthetic_select_on_residue_tree_gtk1(residue_tree, atom_p);

#else

   // actually it would be nice to fix this

#endif

}



// Return pick_info using prestored atom information
// (go_to_atom_* variables).
//
pick_info
graphics_info_t::find_atom_index_from_goto_info(int imol) {

   pick_info pi;

   pi.min_dist = 0; // keep compiler happy
   pi.atom_index = -1; // ditto
   pi.imol = -1; // ditto
   pi.success = GL_FALSE;

   // actually, try ignoring the atom_string argument and consider
   // using mmdb atom selection

   if (imol < 0) {
      std::cout << "WARNING:: no molecule for imol = " << imol << std::endl;

   } else {

      if (imol >= n_molecules()) {
	 std::cout << "WARNING:: no molecule for imol = " << imol << std::endl;
      } else {

	 if (graphics_info_t::molecules[imol].atom_sel.mol == NULL ) {
	    std::cout << "WARNING: (Programmer error) looking for atoms "
		      << "in a molecule " << imol << " which is null" << std::endl;
	 } else {
	    atom_selection_container_t AtomSel =
	       graphics_info_t::molecules[imol].atom_sel;

	    int selHnd = AtomSel.mol->NewSelection(); // d

	    // 0 -> any (mmdb) model
	    //
	    // note that we have to do ugly (char *) casting.  Hopefully
	    // that will go away in new version of mmdb..? (20 Aug 2002 - PE)
	    //
	    std::pair<std::string, std::string> p =
	       graphics_info_t::split_atom_name(go_to_atom_atom_name());

// 	    std::cout << "FAI:: searching chain :" << go_to_atom_chain() << std::endl;
// 	    std::cout << "FAI:: searching residue no :" << go_to_atom_residue() << std::endl;
// 	    std::cout << "FAI:: searching inscode :" << go_to_atom_inscode_ << std::endl;
// 	    std::cout << "FAI:: searching atom_name :" << p.first << std::endl;
// 	    std::cout << "FAI:: searching altconf   :" << altconf << ":" << std::endl;

	    AtomSel.mol->SelectAtoms (selHnd, 0, go_to_atom_chain(),
				      go_to_atom_residue(), // starting resno, an int
				      go_to_atom_inscode_.c_str(), // any insertion code
				      go_to_atom_residue(), // ending resno
				      go_to_atom_inscode_.c_str(), // ending insertion code
				      "*", // any residue name
				      p.first.c_str(), // atom name
				      "*", // elements
				      "*"
				      );

	    int nSelAtoms;
	    mmdb::PPAtom local_SelAtom = NULL;

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
		     if (local_SelAtom[iat]->GetUDData(AtomSel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
			pi.success = GL_TRUE;
			pi.atom_index = ic;
			break;
		     } else {

			// Fall back to comparing pointers.
			//
			for (int i=0; i<AtomSel.n_selected_atoms; i++) {
			   if (AtomSel.atom_selection[i] == local_SelAtom[0]) {
			      pi.success = GL_TRUE;
			      pi.atom_index = i;
			      break;
			   }
			}
		     }
		  }
	       }

	       // The altconf didn't match, so try an atom that is the same except for the altconf:
	       //
	       if (pi.success != GL_TRUE) {
		  if (local_SelAtom[0]->GetUDData(AtomSel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
		     pi.success = GL_TRUE;
		     pi.atom_index = ic;
		  } else {

		     // Fall back to comparing pointers.
		     //
		     for (int i=0; i<AtomSel.n_selected_atoms; i++) {
			if (AtomSel.atom_selection[i] == local_SelAtom[0]) {
			   pi.success = GL_TRUE;
			   pi.atom_index = i;
			   break;
			}
		     }
		  }
	       }

	       if (pi.success == GL_TRUE) {
             go_to_atom_atom_altLoc_ = AtomSel.atom_selection[pi.atom_index]->altLoc;
	       }

	       // should we update the go to atom widget's atom name entry now? (it
	       // is not updated elsewhere)

	    } else {

	       // this happens when the search atom is not in the molecule

	       int selHnd_check = AtomSel.mol->NewSelection();
	       AtomSel.mol->SelectAtoms (selHnd_check, 0, "*",
					 mmdb::ANY_RES, // starting resno, an int
					 "*", // any insertion code
					 mmdb::ANY_RES, // ending resno
					 "*", // ending insertion code
					 "*", // any residue name
					 "*", // atom name
					 "*", // elements
					 "*"  // alt loc.
					 );
	       int nSelAtoms_check;
	       mmdb::PPAtom local_SelAtom_check;
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
      // GtkWidget *residue_tree = lookup_widget(go_to_atom_window, "go_to_atom_residue_tree");
      GtkWidget *residue_tree = widget_from_builder("go_to_atom_residue_tree");
      GtkWidget *atom_list    = widget_from_builder("go_to_atom_atom_list");


      if (false)
         std::cout << ".......................... in update_go_to_atom_window_on_changed_mol() " << imol
                   << " residue_tree " << residue_tree << " atom_list " << atom_list << std::endl;
      if (residue_tree == NULL) {
         std::cout << "ERROR:: residue_tree (go_to_atom_residue_tree) is null!\n";
      } else {
         // graphics_info_t::fill_go_to_atom_residue_tree_and_atom_list_gtk2(imol, residue_tree, atom_list);
         fill_go_to_atom_window_residue_and_atom_lists_gtk4();
      }
   }
}



void
graphics_info_t::apply_go_to_atom_from_widget(GtkWidget *widget) {

   GtkEntry *entry = GTK_ENTRY(widget_from_builder("go_to_atom_chain_entry"));
   const gchar *chain_str = gtk_editable_get_text(GTK_EDITABLE(entry));

   // entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(widget), "go_to_atom_residue_entry"));
   entry = GTK_ENTRY(widget_from_builder("go_to_atom_residue_entry"));
   const gchar *res_str = gtk_editable_get_text(GTK_EDITABLE(entry));

   // entry = GTK_ENTRY(lookup_widget(GTK_WIDGET(widget), "go_to_atom_atom_name_entry"));
   entry = GTK_ENTRY(widget_from_builder("go_to_atom_atom_name_entry"));

   const gchar *txt =  gtk_editable_get_text(GTK_EDITABLE(entry));
   if (txt) {
      std::pair<std::string, std::string> p =
         graphics_info_t::split_atom_name(std::string(txt));
      //      std::cout << "DEBUG: split: " << std::string(txt) << " into :"
      // 	       << p.first << ":  :" << p.second << ":\n" ;

      // we have to use the version of set_go_to.. that has 4 params,
      // because the 3 parameter version sets the altconf to "empty:.
      //
      const gchar *atom_name_str = p.first.c_str();

      std::pair<std::string, std::string> resno_inscode = split_resno_inscode(std::string(res_str));
      int resno = atoi(resno_inscode.first.c_str());
      std::string inscode = resno_inscode.second;

      if (false)
         std::cout << "DEBUG:: in apply_go_to_atom_from_widget(): chain_str " << chain_str
                   << " resno " << resno << " inscode " << inscode << " atom_name_str \""
                   << atom_name_str << "\" split thing \"" << p.second << "\"\n";

      set_go_to_atom_chain_residue_atom_name(chain_str,
                                             resno,
                                             inscode.c_str(),
                                             atom_name_str,
                                             p.second.c_str());
      graphics_grab_focus();

      int success = try_centre_from_new_go_to_atom();
   }
}



// Coordinates have been read into a new molecule (imol).
//
// Now fix up the Go_To_Atom window to match by changing the option menu
//
void
graphics_info_t::update_go_to_atom_window_on_new_mol() {

   if (go_to_atom_window) {

      // GtkWidget *option_menu = lookup_widget(GTK_WIDGET(go_to_atom_window),
      //                                        "go_to_atom_molecule_optionmenu");

      GtkWidget *combobox = widget_from_builder("go_to_atom_molecule_combobox");

      std::cout << "debug:: in update_go_to_atom_window_on_new_mol() got molecule combobox "
                << combobox << std::endl;

      // this may not be the write function for a combobox item
      GCallback callback_func = G_CALLBACK(graphics_info_t::go_to_atom_mol_combobox_changed);

      // set last active (1)
      // fill_option_menu_with_coordinates_options_internal(option_menu, callback_func, 0);

      bool set_last_active_flag = 0;
      gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox));
      fill_combobox_with_coordinates_options_with_set_last(combobox, callback_func, set_last_active_flag);

      // If there was no molecule already, we need to update the atom
      // lists too.

      // Update the residue and atom list to the last displayed
      // molecule.

      // However, if there *was/were* a molecule in existance before
      // the new one was added, we don't want to change the go to atom
      // molecule residue tree (do we?)
      //
      // 20090620:
      //
      //    I wonder why I thought the above.  The behaviour before
      //    todays fix was that on reading a new molecule, the
      //    residue/atom trees were not updated.  How can I have
      //    thought that that was the right thing?
      //
      //    Ah, it is the right thing because the go_to_atom_molecule
      //    is not updated by reading in a new PDB file!  That's why
      //    the atom/residue trees should not be changed (except for
      //    special circumstances).

      int mol_no= -1;
      std::vector<int> imols_existing;
      for (int imol=0; imol<n_molecules(); imol++) {
	 if (molecules[imol].has_model()) {
	    mol_no = imol;
	    imols_existing.push_back(imol);
	 }
      }
      if (mol_no != -1)
         if (imols_existing.size() == 1)
            update_go_to_atom_window_on_changed_mol(mol_no);
   }
}


// Like the above, but don't work out which molecule to update to,
// because we are passed it.
void
graphics_info_t::update_go_to_atom_window_on_other_molecule_chosen(int imol) {

   if (go_to_atom_window) {
      GtkWidget *combobox = widget_from_builder("go_to_atom_molecule_combobox");

      GCallback callback_func = G_CALLBACK(go_to_atom_mol_combobox_changed);
      gtk_cell_layout_clear(GTK_CELL_LAYOUT(combobox));
      fill_combobox_with_coordinates_options(combobox, callback_func, imol);
      update_go_to_atom_window_on_changed_mol(imol);
   }
}

int
graphics_info_t::update_go_to_atom_molecule_on_go_to_atom_molecule_deleted() {

   // Choose the first valid molecule
   //
   int mol_no= -1;
   std::vector<int> imols_existing;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
         mol_no = imol;
         break;
      }
   }
   if (mol_no != -1) {
      go_to_atom_molecule_ = mol_no;
      update_go_to_atom_window_on_changed_mol(mol_no);
   }
   return mol_no;
}


// a static
void
graphics_info_t::clear_atom_list(GtkWidget *atom_gtklist) {

   // gtk_list_clear_items(GTK_LIST(atom_gtklist), 0, -1);

}

// void
// graphics_info_t::fill_go_to_atom_option_menu(GtkWidget *option_menu) {

//    std::cout << "This function should not be called!" << std::endl;

// //    GtkSignalFunc callback_func =
// //       GTK_SIGNAL_FUNC(graphics_info_t::go_to_atom_mol_menu_item_select);

// //    fill_option_menu_with_coordinates_options(option_menu, callback_func);


// }


// // a static
// void
// graphics_info_t::go_to_atom_mol_menu_item_select(GtkWidget *item, GtkPositionType pos) {

//    std::cout << "DEBUG:: (menu item select) Go To Atom molecule now: " << pos << std::endl;
//    graphics_info_t g;
//    int old_go_to_molecule = g.go_to_atom_molecule();
//    g.set_go_to_atom_molecule(pos);

//    if (pos != old_go_to_molecule) {

//       GtkWidget *residue_tree = lookup_widget(GTK_WIDGET(item),
// 						 "go_to_atom_residue_tree");
//       GtkWidget *atom_list = lookup_widget(GTK_WIDGET(item),
// 					   "go_to_atom_atom_list");

//       std::cout << "Debug:: fill_go_to_atom_residue_tree_gtk2 " << pos << std::endl;
//       fill_go_to_atom_residue_tree_and_atom_list_gtk2(pos, residue_tree, atom_list);

//    }
// }

void
graphics_info_t::go_to_atom_mol_combobox_changed(GtkWidget *combobox, gpointer data) {

   GtkTreeIter iter;
   gboolean state = gtk_combo_box_get_active_iter(GTK_COMBO_BOX(combobox), &iter);
   if (state) {
      GtkTreeModel *model = gtk_combo_box_get_model(GTK_COMBO_BOX(combobox));
      GValue label_as_value = { 0, };
      gtk_tree_model_get_value(model, &iter, 0, &label_as_value);
      int imol = g_value_get_int(&label_as_value);
      // std::cout << "imol: " << imol << std::endl;
      graphics_info_t g;
      int old_go_to_molecule = g.go_to_atom_molecule();
      g.set_go_to_atom_molecule(imol);
      if (imol != old_go_to_molecule) {
	 // GtkWidget *residue_tree = lookup_widget(GTK_WIDGET(combobox), "go_to_atom_residue_tree");
	 // GtkWidget *atom_list = lookup_widget(GTK_WIDGET(combobox), "go_to_atom_atom_list");
	 // GtkWidget *residue_tree = widget_from_builder("go_to_atom_residue_tree");
	 // GtkWidget *atom_list = widget_from_builder("go_to_atom_atom_list");
	 // std::cout << "Debug:: fill_go_to_atom_residue_tree_gtk2 " << imol << std::endl;
	 // fill_go_to_atom_residue_tree_and_atom_list_gtk2(imol, residue_tree, atom_list);
         fill_go_to_atom_window_residue_and_atom_lists_gtk4();

      }
   } else {
      std::cout << "bad state" << std::endl;
   }

}

int
graphics_info_t::combobox_get_imol(GtkComboBox *combobox) const {

   int imol = -1;

   if (combobox) {
      GtkTreeIter iter;
      gboolean state = gtk_combo_box_get_active_iter(GTK_COMBO_BOX(combobox), &iter);
      if (state) {
	 GtkTreeModel *model = gtk_combo_box_get_model(GTK_COMBO_BOX(combobox));
	 bool g_value_holds_int_flag = true;
	 // adjust the value of g_value_holds_int_flag here
	 if (g_value_holds_int_flag) {
	    GValue label_as_value = { 0, };
	    gtk_tree_model_get_value(model, &iter, 0, &label_as_value);
	    imol = g_value_get_int(&label_as_value);
	    // std::cout << "DEBUG:: combobox_get_imol() imol: " << imol << std::endl;
	 }
      } else {
	 std::cout << "DEBUG:: combobox_get_imol(): bad state (no active iter in combobox)" << std::endl;
      }
   }
   return imol;
}


void
graphics_info_t::undo_last_move() {  // suggested by Frank von Delft


   coot::Cartesian c = get_old_rotation_centre();

   std::cout << "INFO:: Moving back to old centre: " << c << std::endl;
   setRotationCentre(c);
   for(int ii=0; ii<n_molecules(); ii++) {
      molecules[ii].update_map(graphics_info_t::auto_recontour_map_flag);
      molecules[ii].update_symmetry();
   }
   graphics_draw();
}

// static
std::pair<bool, std::pair<int, coot::atom_spec_t> >
graphics_info_t::active_atom_spec(int imol) {

   return active_atom_spec_internal(imol);
}

std::pair<bool, std::pair<int, coot::atom_spec_t> >
graphics_info_t::active_atom_spec() {
   return active_atom_spec_internal(-1); // any molecule
}

// if imol_only is -1 then we can allow any molecule to be active.  If it
// is >=0, only allow molecule imol_only to have the the active.
//
std::pair<bool, std::pair<int, coot::atom_spec_t> >
graphics_info_t::active_atom_spec_internal(int imol_only) {

   coot::atom_spec_t spec;
   bool was_found_flag = 0;

   graphics_info_t g;
   float dist_best = 999999999.9;
   int imol_closest = -1;
   mmdb::Atom *at_close = 0;

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {

      if ((imol_only == -1) || (imol == imol_only)) {
	 if (is_valid_model_molecule(imol)) {
	    if (graphics_info_t::molecules[imol].is_displayed_p()) {
	       if (graphics_info_t::molecules[imol].atom_selection_is_pickable()) {
		  coot::at_dist_info_t at_info =
		     graphics_info_t::molecules[imol].closest_atom(g.RotationCentre());
		  if (at_info.atom) {
		     if (at_info.dist <= dist_best) {
			dist_best = at_info.dist;
			imol_closest = at_info.imol;
			at_close = at_info.atom;
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (at_close) {
      spec = coot::atom_spec_t(at_close);
      was_found_flag = 1;
   }

   std::pair<int, coot::atom_spec_t> p1(imol_closest, spec);
   return std::pair<bool, std::pair<int, coot::atom_spec_t> > (was_found_flag, p1);
}

std::pair<int, mmdb::Atom *>
graphics_info_t::get_active_atom() const {

   mmdb::Atom *at_close = 0;
   float dist_best = 10.0; // atoms must be closer than this to be "active" 10.0 is too generous?
   int imol_closest = -1;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (true) {
	 if (is_valid_model_molecule(imol)) {
	    if (molecules[imol].is_displayed_p()) {
	       if (molecules[imol].atom_selection_is_pickable()) {
		  coot::at_dist_info_t at_info =
		     molecules[imol].closest_atom(RotationCentre());
		  if (at_info.atom) {
		     if (at_info.dist <= dist_best) {
			dist_best = at_info.dist;
			imol_closest = at_info.imol;
			at_close = at_info.atom;
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (at_close)
      return std::pair<int, mmdb::Atom *>(imol_closest, at_close);
   else
      return std::pair<int, mmdb::Atom *>(-1, 0);
}

std::pair<bool, std::pair<int, coot::atom_spec_t> >
graphics_info_t::active_atom_spec_simple() {

   int imol_closest = -1;
   coot::atom_spec_t spec;
   bool was_found_flag = false;
   float dist_best = 999999999.9;
   mmdb::Atom *at_close = 0;

   graphics_info_t g;
   coot::Cartesian rc = g.RotationCentre(); // not static!
   for (int imol=0; imol<n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
         if (molecules[imol].is_displayed_p()) {
            if (molecules[imol].atom_selection_is_pickable()) {
               bool do_ca_check_flag = false;
               coot::at_dist_info_t at_info =
                  molecules[imol].closest_atom(rc, do_ca_check_flag, "", false);
               if (at_info.atom) {
                  if (at_info.dist <= dist_best) {
                     dist_best = at_info.dist;
                     imol_closest = at_info.imol;
                     at_close = at_info.atom;
                  }
               }
            }
         }
      }
   }
   if (at_close) {
      spec = coot::atom_spec_t(at_close);
      was_found_flag = true;
   }

   std::pair<int, coot::atom_spec_t> p1(imol_closest, spec);
   return std::pair<bool, std::pair<int, coot::atom_spec_t> > (was_found_flag, p1);
}


// static
void
graphics_info_t::apply_go_to_residue_keyboading_string(const std::string &ks) {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aas = active_atom_spec();
   coot::Cartesian rc = g.RotationCentre();

   if (aas.first) {
      if (! coot::sequence::is_sequence_triplet(ks)) {
	 int imol = aas.second.first;
	 mmdb::Atom *new_centre_atom = g.molecules[imol].get_atom(ks, aas.second.second, rc);
	 if (new_centre_atom) {
	    g.apply_go_to_residue_keyboading_string_inner(imol,new_centre_atom);
	 } else {
	    new_centre_atom = g.molecules[imol].get_atom(coot::util::upcase(ks), aas.second.second, rc);
	    g.apply_go_to_residue_keyboading_string_inner(imol, new_centre_atom);
	 }
      } else {
	 // handle sequence triplet
	 int imol = aas.second.first;
	 g.apply_go_to_residue_from_sequence_triplet(imol, ks);
      }
   } else {
      std::cout << "WARNING:: No active atom " << std::endl;
   }
}

void
graphics_info_t::apply_go_to_residue_keyboading_string_inner(int imol, mmdb::Atom *new_centre_atom) {

   if (new_centre_atom) {
      int index;
      coot::Cartesian new_pt(new_centre_atom->x,
                             new_centre_atom->y,
                             new_centre_atom->z);
      setRotationCentre(new_pt);
      //for (int ii=0; ii<graphics_info_t::n_molecules(); ii++) {
      //	 graphics_info_t::molecules[ii].update_map();
      //	 graphics_info_t::molecules[ii].update_clipper_skeleton();
      //     }
      //    graphics_draw();
      update_things_on_move_and_redraw();
      set_go_to_atom_molecule(imol);
      set_go_to_atom_chain_residue_atom_name(new_centre_atom->GetChainID(),
                                             new_centre_atom->GetSeqNum(),
                                             new_centre_atom->GetAtomName());
      update_go_to_atom_window_on_other_molecule_chosen(imol);

      // from setRotationCentre (if used with atom index)
      // to get update on distances, labels etc.
      if (new_centre_atom->GetUDData(molecules[imol].atom_sel.UDDAtomIndexHandle, index) == mmdb::UDDATA_Ok) {
        if (environment_show_distances) {
          mol_no_for_environment_distances = imol;
          update_environment_graphics_object(index, imol);
          // label new centre
          if (environment_distance_label_atom) {
            molecules[imol].unlabel_last_atom();
            molecules[imol].add_to_labelled_atom_list(index);
          }
          if (show_symmetry)
            update_symmetry_environment_graphics_object(index, imol);
        } else {

          if (label_atom_on_recentre_flag) {
            molecules[imol].unlabel_last_atom();
            molecules[imol].add_to_labelled_atom_list(index);
          }
        }
      } else {
        std::cout << "WARNING:: failed to find index. No updating of labels and distanced" << std::endl;
      }
   } else {
     std::cout << "WARNING:: failed to find that residue - no new centre atom " << std::endl;
   }
}

// go to the middle residue of the first occurance of the sequence triplet if you can
// seq_trip is of course something like "ACE"
// return the "found the triplet and moved there" status: 0 for fail.
//
int graphics_info_t::apply_go_to_residue_from_sequence_triplet(int imol, const std::string &seq_trip) {

   int status = 0;

   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *new_centre_atom = graphics_info_t::molecules[imol].get_centre_atom_from_sequence_triplet(seq_trip);
      std::cout << "INFO:: new centre atom: " << new_centre_atom << std::endl;
      if (new_centre_atom)
	 apply_go_to_residue_keyboading_string_inner(imol, new_centre_atom);
   }
   return status;
}



// --- unapply symmetry to current view, (we are looking at
// symmetry and we want to get back to the main molecule,
// preserving the orientation, if possible.
//
// pre_translation is the translation that needed to be applied to
// the molecule so that it was close to the origin (from there we
// do the symmetry expansion, and it seems that st generates the
// symmetry-related molecule that we are looking at now).
int
graphics_info_t::unapply_symmetry_to_view(int imol, const std::vector<std::pair<clipper::RTop_orth, clipper::Coord_orth> > &symm_mat_and_pre_shift) {

   int r=0; // a symm was found and used
   float min_dist = 99999999999.9;
   coot::Cartesian rc = RotationCentre();
   clipper::Coord_orth centre_pt(rc.x(), rc.y(), rc.z());
   clipper::RTop_orth  best_rtop_symm;
   clipper::Coord_orth best_molecule_centre(0,0,0);
   for (unsigned int i=0; i<symm_mat_and_pre_shift.size(); i++) {

      clipper::RTop_orth  rtop_symm = symm_mat_and_pre_shift[i].first;
      clipper::Coord_orth pre_shift = symm_mat_and_pre_shift[i].second;

      clipper::Coord_orth pt_1 = centre_pt.transform(rtop_symm.inverse());
      clipper::Coord_orth pt_2 = pt_1 + pre_shift;

      if (false) {
         std::cout << "box =================== " << i << " ======================= " << std::endl;
         std::cout << "rtop_symm:\n" << rtop_symm.format() << std::endl;
         std::cout << "pre_shift: " << pre_shift.format() << std::endl;
      }

      // Now, is pt_2 close to an atom in the imolth molecule?  If so, what is the distance?
      // that function uses a coot::Cartesian
      //
      coot::Cartesian pt_2c(pt_2.x(), pt_2.y(), pt_2.z());
      std::pair<float, int> na = molecules[imol].nearest_atom(pt_2c);
      if (na.second >= 0) {
         if (na.first < min_dist) {
            min_dist = na.first;
            best_rtop_symm = rtop_symm;
            best_molecule_centre = pt_2;
            r = 1;
         }
      }
   }
   coot::Cartesian nrc(best_molecule_centre.x(), best_molecule_centre.y(), best_molecule_centre.z());

   if (r) {
      clipper::Mat33<double> symm_mat = best_rtop_symm.inverse().rot();
      glm::mat3 glm_symm_mat(symm_mat(0,0), symm_mat(0,1), symm_mat(0,2),
                             symm_mat(1,0), symm_mat(1,1), symm_mat(1,2),
                             symm_mat(2,0), symm_mat(2,1), symm_mat(2,2));
      glm::quat glm_sym_mat_quat = glm::toQuat(glm_symm_mat);
      view_quaternion *= glm_sym_mat_quat;
      setRotationCentre(nrc);
      update_things_on_move_and_redraw();
      graphics_draw();
   }
   return r;
}


void graphics_info_t::register_user_defined_interesting_positions(const std::vector<std::pair<clipper::Coord_orth, std::string> > & udip) {

   user_defined_interesting_positions = udip;
   user_defined_interesting_positions_idx = 0;
}


void
graphics_info_t::sequence_view_highlight_residue_maybe(mmdb::Atom *atom, GtkWidget *svc) {

#ifdef HAVE_GOOCANVAS
   if (svc) {
      if (atom) {
         mmdb::Residue *residue_p = atom->residue;
         if (residue_p) {
            exptl::nsv *nsv = static_cast<exptl::nsv *>(g_object_get_data(G_OBJECT(svc), "nsv"));
            if (nsv) {
               nsv->highlight_residue(residue_p);
            }
         }
      }
   }
#endif // HAVE_GOOCANVAS
}


