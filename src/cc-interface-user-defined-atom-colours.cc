/* src/cc-interface-user-defined-colours.cc
 * 
 * Copyright 2016 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "cc-interface-user-defined-atom-colours.hh"

#include "c-interface.h"   // for is_valid_model_molecule
#include "cc-interface.hh" // residue_spec_from_py

/*  ----------------------------------------------------------------------- */
/*                  user-defined atom colours                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
void set_user_defined_atom_colour_by_residue_py(int imol, PyObject *residue_specs_colour_index_tuple_list_py) {

   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(residue_specs_colour_index_tuple_list_py)) {
	 unsigned int l = PyObject_Length(residue_specs_colour_index_tuple_list_py);
	 if (l > 0) {
	    std::vector<std::pair<coot::residue_spec_t, int> > cis;
	    for (unsigned int i=0; i<l; i++) {
	       PyObject *tuple_py = PyList_GetItem(residue_specs_colour_index_tuple_list_py, i);
	       if (PyTuple_Check(tuple_py)) {
		  unsigned int l2 = PyObject_Length(tuple_py);
		  if (l2 == 2) {
		     PyObject *spec_py = PyTuple_GetItem(tuple_py, 0);
		     PyObject *idx_py  = PyTuple_GetItem(tuple_py, 1);
		     if (PyLong_Check(idx_py)) {
			coot::residue_spec_t spec = residue_spec_from_py(spec_py);
			long ci = PyLong_AsLong(idx_py);
			std::pair<coot::residue_spec_t, int> p(spec, ci);
			cis.push_back(p);
		     }
		  }
	       }
	    }
	    graphics_info_t::molecules[imol].set_user_defined_colour_indices_by_residues(cis);
	 }
      }
   }
}

void set_user_defined_atom_colour_py(int imol, PyObject *atom_specs_colour_index_tuple_list_py) {

   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(atom_specs_colour_index_tuple_list_py)) {
	 unsigned int l = PyObject_Length(atom_specs_colour_index_tuple_list_py);
	 if (l > 0) {
	    std::vector<std::pair<coot::atom_spec_t, int> > cis;
	    for (unsigned int i=0; i<l; i++) {
	       PyObject *tuple_py = PyList_GetItem(atom_specs_colour_index_tuple_list_py, i);
	       if (PyTuple_Check(tuple_py)) {
		  unsigned int l2 = PyObject_Length(tuple_py);
		  if (l2 == 2) {
		     PyObject *spec_py = PyTuple_GetItem(tuple_py, 0);
		     PyObject *idx_py  = PyTuple_GetItem(tuple_py, 1);
		     if (PyLong_Check(idx_py)) {
			coot::atom_spec_t spec = atom_spec_from_python_expression(spec_py);
			long ci = PyLong_AsLong(idx_py);
			std::pair<coot::atom_spec_t, int> p(spec, ci);
			cis.push_back(p);
		     }
		  }
	       }
	    }
	    graphics_info_t::molecules[imol].set_user_defined_colour_indices(cis);
	 }
      }
   }
}

#endif // USE_PYTHON

#ifdef USE_GUILE
void set_user_defined_atom_colour_scm(int imol, SCM atom_specs_colour_index_tuple_list_scm) {

}

void set_user_defined_atom_colour_by_residue_scm(int imol, SCM residue_specs_colour_index_tuple_list_scm) {

}
#endif // USE_GUILE



void clear_user_defined_atom_colours(int imol) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].clear_user_defined_atom_colours();
      

}
