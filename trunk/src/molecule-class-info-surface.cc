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
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifdef _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include <vector>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "graphics-info.h"

#include "molecule-class-info.h"


// Atom-based LMB (e.g. electrostatic) surface.
// 
void
molecule_class_info_t::draw_surface() {

   if (drawit == 1) { 
      if (cootsurface) {
	 glEnable(GL_LIGHTING);
	 glEnable(GL_LIGHT0);
	 glCallList(theSurface);
	 glDisable(GL_LIGHT0);
	 glDisable(GL_LIGHTING);
      }
   }
}


void
molecule_class_info_t::make_surface(int on_off_flag,
				    const coot::protein_geometry &geom,
				    float col_scale) {

   if (atom_sel.n_selected_atoms > 0) {

      apply_charges(geom);
      
      if (on_off_flag == 0) {
	 cootsurface = NULL;
	 glDeleteLists(theSurface, 1);
      } else {
	 glDeleteLists(theSurface, 1);
	 theSurface = glGenLists(1);
	 glNewList(theSurface, GL_COMPILE);
	 cootsurface = new coot::surface;
	 cootsurface->fill_from(atom_sel.mol, atom_sel.SelectionHandle, col_scale);
	 if (cootsurface) 
	    cootsurface->draw(0, 0);
	 glEndList();
      }
   }
}

// a generic function to convert from a residue_spec_vec to a selection handle.
//
void
molecule_class_info_t::fill_residue_selection(int SelHnd_selection,
					      const std::vector<coot::residue_spec_t> &res_specs_vec) {

   for (unsigned int ir=0; ir<res_specs_vec.size(); ir++) {
      atom_sel.mol->SelectAtoms(SelHnd_selection, 0,
				res_specs_vec[ir].chain.c_str(),
				res_specs_vec[ir].resno,
				res_specs_vec[ir].insertion_code.c_str(), 
				res_specs_vec[ir].resno,
				res_specs_vec[ir].insertion_code.c_str(),
				"*", "*", "*", "*", SKEY_OR);
   }

   PPCAtom atoms = NULL;
   int n_atoms;
   atom_sel.mol->GetSelIndex(SelHnd_selection, atoms, n_atoms);
   std::cout << "debug:: fill_residue_selection selected "
	     << n_atoms << " atoms" << std::endl;
} 


void
molecule_class_info_t::make_surface(const std::vector<coot::residue_spec_t> &res_specs_vec,
				    const coot::protein_geometry &geom) {

   float col_scale = 0.2; // pass this 

   if (0) { // this is how it should work
   int SelHnd_selection = atom_sel.mol->NewSelection();

   // fill_residue_selection(SelHnd_selection, res_specs_vec);

   // make_surface(SelHnd_selection, atom_sel.SelectionHandle, geom);
   make_surface(atom_sel.SelectionHandle, atom_sel.SelectionHandle, geom);
   
   atom_sel.mol->DeleteSelection(SelHnd_selection);
   } else {
      make_surface(1, geom, col_scale); // old function.
   } 
} 


// use the other version of make surface to turn this surface off.
// 
void
molecule_class_info_t::make_surface(int SelHnd_selection, int SelHnd_all,
				    const coot::protein_geometry &geom) {

   float col_scale = 0.2; // pass this
   
   glDeleteLists(theSurface, 1);
   theSurface = glGenLists(1);
   glNewList(theSurface, GL_COMPILE);
   cootsurface = new coot::surface;
   // cootsurface->fill_surface(atom_sel.mol, SelHnd_selection, SelHnd_all);
   cootsurface->fill_surface(atom_sel.mol, SelHnd_all, SelHnd_all, col_scale);
   if (cootsurface) 
      cootsurface->draw(0, 0);
   glEndList();
}


// electron density solid surface
void
molecule_class_info_t::do_solid_surface_for_density(short int on_off_flag) {

}



