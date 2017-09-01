#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

// #include "compat/coot-sysdep.h"

// #include <stdlib.h>
#include <iostream>

// #include <vector>
// #include <string>
// #include <algorithm>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#include "coords/mmdb-crystal.h"

#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "graphics-info.h"

#include "coot-utils/coot-coord-utils.hh"

#include "c-interface.h"
// #include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"

#ifdef USE_PYTHON

//! \brief return a Python object for the bonds
//
PyObject *get_bonds_represenation(int imol) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {

      r = PyTuple_New(2);

      // write a class that also returns atom specs for the atom. Possibly use inheritance
      // or encapsulation.
      graphical_bonds_container bonds_box = graphics_info_t::molecules[imol].get_bonds_represenation();
      if (bonds_box.atom_centres_) {
	 PyObject *all_atom_positions_py = PyTuple_New(bonds_box.n_consolidated_atom_centres);
	 for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
	    PyObject *atom_set_py = PyTuple_New(bonds_box.consolidated_atom_centres[icol].num_points);
	    for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
	       const coot::Cartesian &pt = bonds_box.consolidated_atom_centres[icol].points[i].second;
	       bool is_H_flag     = bonds_box.consolidated_atom_centres[icol].points[i].first;
	       PyObject *atom_triple_py = PyTuple_New(3);
	       PyObject *coords_py = PyTuple_New(3);
	       std::string s = "attrib-filled-later";
	       // hack!
	       mmdb::Atom *at = graphics_info_t::molecules[imol].get_atom_at_pos(pt);
	       if (at) {
		  coot::atom_spec_t spec(at);
		  s = spec.format();
	       }
	       PyTuple_SetItem(coords_py, 0, PyFloat_FromDouble(pt.x()));
	       PyTuple_SetItem(coords_py, 1, PyFloat_FromDouble(pt.y()));
	       PyTuple_SetItem(coords_py, 2, PyFloat_FromDouble(pt.z()));
	       PyTuple_SetItem(atom_triple_py, 0, coords_py);
	       PyTuple_SetItem(atom_triple_py, 1, PyBool_FromLong(is_H_flag));
	       PyTuple_SetItem(atom_triple_py, 2, PyString_FromString(s.c_str()));
	       PyTuple_SetItem(atom_set_py, i, atom_triple_py);
	    }
	    PyTuple_SetItem(all_atom_positions_py, icol, atom_set_py);
	 }
	 PyTuple_SetItem(r, 0, all_atom_positions_py);
      } else {
	 PyObject *empty_py = PyTuple_New(0);
	 PyTuple_SetItem(r, 0, empty_py);
      }
      PyObject *bonds_tuple = PyTuple_New(bonds_box.num_colours);
      for (int i=0; i<bonds_box.num_colours; i++) {
	 graphical_bonds_lines_list &ll = bonds_box.bonds_[i];
	 PyObject *line_set_py = PyTuple_New(bonds_box.bonds_[i].num_lines);
	 for (int j=0; j< bonds_box.bonds_[i].num_lines; j++) {
	    PyObject *p0_py   = PyTuple_New(3);
	    PyObject *p1_py   = PyTuple_New(3);
	    PyObject *pair_py = PyTuple_New(2);
	    PyTuple_SetItem(p0_py, 0, PyFloat_FromDouble(ll.pair_list[j].positions.getStart().get_x()));
	    PyTuple_SetItem(p0_py, 1, PyFloat_FromDouble(ll.pair_list[j].positions.getStart().get_y()));
	    PyTuple_SetItem(p0_py, 2, PyFloat_FromDouble(ll.pair_list[j].positions.getStart().get_z()));
	    PyTuple_SetItem(p1_py, 0, PyFloat_FromDouble(ll.pair_list[j].positions.getFinish().get_x()));
	    PyTuple_SetItem(p1_py, 1, PyFloat_FromDouble(ll.pair_list[j].positions.getFinish().get_y()));
	    PyTuple_SetItem(p1_py, 2, PyFloat_FromDouble(ll.pair_list[j].positions.getFinish().get_z()));
	    PyTuple_SetItem(pair_py, 0, p0_py);
	    PyTuple_SetItem(pair_py, 1, p1_py);
	    PyTuple_SetItem(line_set_py, j, pair_py);
	 }
	 PyTuple_SetItem(bonds_tuple, i, line_set_py);
      }
      PyTuple_SetItem(r, 1, bonds_tuple);
   }
   if (PyBool_Check(r)) { 
      Py_INCREF(r);
   }
   return r;

}

#endif
