/* src/cc-interface-user-defined-colours.cc
 * 
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include <Python.h>

#include "cc-interface-user-defined-atom-colours.hh"

#include "c-interface.h"   // for is_valid_model_molecule
#include "cc-interface.hh" // residue_spec_from_py

#include "graphics-info.h"

/*  ----------------------------------------------------------------------- */
/*                  user-defined atom colours                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
void set_user_defined_atom_colour_by_selection_py(int imol, PyObject *CID_colour_index_tuple_list_py) {

   // 20220707-PE You are passing a list of tuples, right?

   unsigned int n_new_colours = 0;
   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(CID_colour_index_tuple_list_py)) {
	 unsigned int l = PyObject_Length(CID_colour_index_tuple_list_py);
	 if (l > 0) {
	    std::vector<std::pair<std::string, unsigned int> > cis;
	    for (unsigned int i=0; i<l; i++) {
	       PyObject *tuple_py = PyList_GetItem(CID_colour_index_tuple_list_py, i);
	       if (PyTuple_Check(tuple_py)) {
		  unsigned int l2 = PyObject_Length(tuple_py);
		  if (l2 == 2) {
		     PyObject *selection_py  = PyTuple_GetItem(tuple_py, 0);
		     PyObject *idx_py        = PyTuple_GetItem(tuple_py, 1);
                     if (PyUnicode_Check(selection_py)) {
                        std::string CID_selection = PyBytes_AS_STRING(PyUnicode_AsUTF8String(selection_py));
                        if (PyLong_Check(idx_py)) {
                           long ci = PyLong_AsLong(idx_py);
                           if (ci >= 0) {
                              std::pair<std::string , unsigned int> p(CID_selection, ci);
                              cis.push_back(p);
                           }
                        }
                     }
                  }
               }
            }
            // this sets apply_colour_to_non_carbon_atoms_also to true
            graphics_info_t::molecules[imol].set_user_defined_colour_indices_by_selections(cis);
            n_new_colours = cis.size();
         }
      }
   }
   std::cout << "DEBUG:: set_user_defined_atom_colour_by_selection_py() imol: " << imol
             << " n_new_colours: " << n_new_colours << std::endl;
}

#include "c-interface-python.hh"

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


#ifdef USE_PYTHON
void set_user_defined_colours_py(PyObject *colour_list_in_py) {

   if (PyList_Check(colour_list_in_py)) {
      unsigned int l = PyObject_Length(colour_list_in_py);
      if (l > 0) {
         std::vector<std::pair<unsigned int, coot::colour_holder> > colours;
         for (unsigned int i=0; i<l; i++) {
            PyObject *item_py = PyList_GetItem(colour_list_in_py, i);
            if (PyTuple_Check(item_py)) {
               unsigned int l2 = PyObject_Length(item_py);
               if (l2 == 2) {
                  if (false)
                     std::cout << "l2 = 2 for " << item_py << std::endl;
                  PyObject *colour_index_py = PyTuple_GetItem(item_py, 0);
                  PyObject *colour_list_py  = PyTuple_GetItem(item_py, 1);
                  // std::cout << "debug colour_index_py " << colour_index_py << std::endl;
                  // std::cout << "colour_list_py " << colour_list_py << std::endl;
                  if (colour_index_py) {
                     if (colour_list_py) {
                        if (PyLong_Check(colour_index_py)) {
                           long colour_index = PyLong_AsLong(colour_index_py);
                           if (PyList_Check(colour_list_py)) {
                              unsigned int l3 = PyObject_Length(colour_list_py);
                              if (l3 == 3) {
                                 double r = PyFloat_AsDouble(PyList_GetItem(colour_list_py, 0));
                                 double g = PyFloat_AsDouble(PyList_GetItem(colour_list_py, 1));
                                 double b = PyFloat_AsDouble(PyList_GetItem(colour_list_py, 2));
                                 coot::colour_holder ch(r,g,b);
                                 colours.push_back(std::make_pair(colour_index, ch));
                              }
                           }
                        }
                     }
                  }
               }
            } else {
               std::cout << "DEBUG:: not a tuple " << item_py << std::endl;
            }
         }
         graphics_info_t::set_user_defined_colours(colours);
      }
   }
}

#endif
