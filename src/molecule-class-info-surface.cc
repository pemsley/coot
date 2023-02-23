/* src/molecule-class-info-surface.cc
 * 
 * Copyright 2004, 2006 by The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA.
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#ifdef _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include <vector>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.h"

#include "graphics-info.h"

#include "molecule-class-info.h"


// Atom-based LMB (e.g. electrostatic) surface.
// 
void
molecule_class_info_t::draw_surface() {

#if 0
   if (draw_it == 1) { 
      if (cootsurface) {
	 glEnable(GL_LIGHTING);
	 glEnable(GL_LIGHT0);
	 if (transparent_molecular_surface_flag)
	    draw_transparent_molecular_surface();
	 else 
	    glCallList(theSurface);
	 glDisable(GL_LIGHT0);
	 glDisable(GL_LIGHTING);
      }
   }
#endif
}

void
molecule_class_info_t::draw_transparent_molecular_surface() {
#if 0
   float opacity = 0.4; // pass this 
   if (cootsurface) {
      cootsurface->transparent_draw(opacity);
   }
#endif
} 


void
molecule_class_info_t::make_surface(int on_off_flag,
				    const coot::protein_geometry &geom,
				    float col_scale) {

#if 0 // I'm giving up with surfaces for now in 0.9.x
   if (atom_sel.n_selected_atoms > 0) {

      bool charges_applied_flag = apply_charges(geom);
      
      if (on_off_flag == 0) {
	 cootsurface = NULL;
	 glDeleteLists(theSurface, 1);
      } else {
	 glDeleteLists(theSurface, 1);
	 theSurface = glGenLists(1);
	 glNewList(theSurface, GL_COMPILE);
	 cootsurface = new coot::surface;
	 bool need_charges_assigned = 1;
	 if (charges_applied_flag)
	    need_charges_assigned = 0;

	 if (0) 
	    std::cout << "in molecule_class_info_t::make_surface() charges_applied_flag "
		      << charges_applied_flag
		      << " need_charges_assigned " << need_charges_assigned << std::endl;
	 
	 cootsurface->fill_from(atom_sel.mol, atom_sel.SelectionHandle, col_scale,
				need_charges_assigned);
	 if (cootsurface) 
	    cootsurface->draw(0, 0);
	 glEndList();
      }
   }
#endif
}

// a generic function to convert from a residue_spec_vec to a selection handle.
//
void
molecule_class_info_t::fill_residue_selection(int SelHnd_selection,
					      const std::vector<coot::residue_spec_t> &res_specs_vec,
					      bool allow_waters_flag) {

   std::string no_waters = "!HOH";
   if (allow_waters_flag)
      no_waters = "*";
   
   for (unsigned int ir=0; ir<res_specs_vec.size(); ir++) {
      atom_sel.mol->SelectAtoms(SelHnd_selection, 0,
				res_specs_vec[ir].chain_id.c_str(),
				res_specs_vec[ir].res_no,
				res_specs_vec[ir].ins_code.c_str(), 
				res_specs_vec[ir].res_no,
				res_specs_vec[ir].ins_code.c_str(),
				no_waters.c_str(),
				"*", "*", "*", mmdb::SKEY_OR);
   }

   mmdb::PPAtom atoms = NULL;
   int n_atoms;
   atom_sel.mol->GetSelIndex(SelHnd_selection, atoms, n_atoms);
   std::cout << "debug:: fill_residue_selection selected "
	     << n_atoms << " atoms" << std::endl;
} 


void
molecule_class_info_t::make_surface(const std::vector<coot::residue_spec_t> &res_specs_vec,
				    const coot::protein_geometry &geom,
				    float col_scale) {
#if 0
   // res_specs_vec must not contain waters or bad things might
   // happen.

   int SelHnd_selection = atom_sel.mol->NewSelection();
   fill_residue_selection(SelHnd_selection, res_specs_vec, 0);

   // no waters in protein (leaving the waters in the selection makes
   // a rotten surface)
   // 
   int SelHnd_protein = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(SelHnd_protein, 0,
			     "*",
			     mmdb::ANY_RES, "*",
			     mmdb::ANY_RES, "*",
			     "!HOH",    // RNames
			     "*", "*",  // ANames, Elements
			     "*" );     // Alternate locations.
   
   make_surface(SelHnd_selection, SelHnd_protein, geom, col_scale);
   
   atom_sel.mol->DeleteSelection(SelHnd_selection);
   atom_sel.mol->DeleteSelection(SelHnd_protein);

#endif
}


// use the other version of make surface to turn this surface off.
// 
void
molecule_class_info_t::make_surface(int SelHnd_selection, int SelHnd_all,
				    const coot::protein_geometry &geom,
				    float col_scale) {

#if 0
   glDeleteLists(theSurface, 1);
   theSurface = glGenLists(1);
   glNewList(theSurface, GL_COMPILE);
   cootsurface = new coot::surface;
   bool need_charges_assigned = 1;
   cootsurface->fill_surface(atom_sel.mol, SelHnd_selection, SelHnd_all, col_scale,
			     need_charges_assigned);
   if (cootsurface) 
      cootsurface->draw(0, 0);
   glEndList();
#endif
}


