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
PyObject *get_bonds_representation(int imol) {

   // Because "residue picking/selection" is crucial for so many coot tools, we want to support
   // a mechanism that allows the client-side display of the "active" residue, we do that by providing
   // a "residue-index" for every atom and bond.

   // currently we use get_at_pos() to convert from coordinates to a atom.
   // This is easier to program, but slow and really a silly thing to do.  It would be better
   // to add the residue_index for the atoms and bonds when the 
   // bonds and atom positions are created.
   // if get_bonds_respresenation() is slow, go back and fix it.
   // 3wj7 takes 0.78s for just the atom loop get_at_pos().

   // Carteian doesn't yet work as key of a map
   // std::map<coot::Cartesian, mmdb::Residue *> residue_map;


   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {

      r = PyTuple_New(2);

      // write a class that also returns atom specs for the atom. Possibly use inheritance
      // or encapsulation.
      graphical_bonds_container bonds_box = graphics_info_t::molecules[imol].get_bonds_representation();
      if (bonds_box.atom_centres_) {
	 PyObject *all_atom_positions_py = PyTuple_New(bonds_box.n_consolidated_atom_centres);
	 for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
	    PyObject *atom_set_py = PyTuple_New(bonds_box.consolidated_atom_centres[icol].num_points);
	    for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
	       const coot::Cartesian &pt = bonds_box.consolidated_atom_centres[icol].points[i].position;
	       bool is_H_flag     = bonds_box.consolidated_atom_centres[icol].points[i].is_hydrogen_atom;
	       long residue_index = bonds_box.consolidated_atom_centres[icol].points[i].residue_index;
	       PyObject *atom_info_quad_py = PyTuple_New(4);
	       PyObject *coords_py = PyTuple_New(3);
	       std::string s = "attrib-filled-later";
	       // Hmm.
	       // Perhaps we actually want the atom spec (as a python object).
	       // In that case, graphical_bonds_atom_info_t should store the atom, not the residue
	       //
	       if (bonds_box.consolidated_atom_centres[icol].points[i].residue_p) {
		  coot::residue_spec_t spec(bonds_box.consolidated_atom_centres[icol].points[i].residue_p);
		  s = spec.format();
	       }

	       // looking up the atoms takes ~1s for a reasonable-sized protein (3jw7)
	       if (false) { 	       // hack!
		  mmdb::Atom *at = graphics_info_t::molecules[imol].get_atom_at_pos(pt);
		  if (at) {
		     // std::cout << "successful atom residue-index lookup " << coot::atom_spec_t(at)
		     // << " at " << pt << "\n";
		     coot::atom_spec_t spec(at);
		     s = spec.format();
		  }
	       }

	       PyObject *residue_index_py = PyInt_FromLong(residue_index);
	       PyTuple_SetItem(coords_py, 0, PyFloat_FromDouble(pt.x()));
	       PyTuple_SetItem(coords_py, 1, PyFloat_FromDouble(pt.y()));
	       PyTuple_SetItem(coords_py, 2, PyFloat_FromDouble(pt.z()));
	       PyTuple_SetItem(atom_info_quad_py, 0, coords_py);
	       PyTuple_SetItem(atom_info_quad_py, 1, PyBool_FromLong(is_H_flag));
	       PyTuple_SetItem(atom_info_quad_py, 2, PyString_FromString(s.c_str()));
	       PyTuple_SetItem(atom_info_quad_py, 3, residue_index_py);
	       PyTuple_SetItem(atom_set_py, i, atom_info_quad_py);
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
	 graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[i];
	 PyObject *line_set_py = PyTuple_New(bonds_box.bonds_[i].num_lines);
	 for (int j=0; j< bonds_box.bonds_[i].num_lines; j++) {
	    const graphics_line_t::cylinder_class_t &cc = ll.pair_list[j].cylinder_class;
            int ri = ll.pair_list[j].residue_index; // set to -1 by constructor, overwite if possible

	    PyObject *p0_py   = PyTuple_New(3);
	    PyObject *p1_py   = PyTuple_New(3);
	    PyObject *positions_and_order_py = PyTuple_New(4);
	    PyObject *order_py = PyInt_FromLong(cc);
            PyObject *residue_index_py = PyInt_FromLong(ri);
	    PyTuple_SetItem(p0_py, 0, PyFloat_FromDouble(ll.pair_list[j].positions.getStart().get_x()));
	    PyTuple_SetItem(p0_py, 1, PyFloat_FromDouble(ll.pair_list[j].positions.getStart().get_y()));
	    PyTuple_SetItem(p0_py, 2, PyFloat_FromDouble(ll.pair_list[j].positions.getStart().get_z()));
	    PyTuple_SetItem(p1_py, 0, PyFloat_FromDouble(ll.pair_list[j].positions.getFinish().get_x()));
	    PyTuple_SetItem(p1_py, 1, PyFloat_FromDouble(ll.pair_list[j].positions.getFinish().get_y()));
	    PyTuple_SetItem(p1_py, 2, PyFloat_FromDouble(ll.pair_list[j].positions.getFinish().get_z()));
	    PyTuple_SetItem(positions_and_order_py, 0, p0_py);
	    PyTuple_SetItem(positions_and_order_py, 1, p1_py);
	    PyTuple_SetItem(positions_and_order_py, 2, order_py);
	    PyTuple_SetItem(positions_and_order_py, 3, residue_index_py);
	    PyTuple_SetItem(line_set_py, j, positions_and_order_py);
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
