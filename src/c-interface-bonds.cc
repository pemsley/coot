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

#include "c-interface-bonds.hh"

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

   // Cartesian doesn't yet work as key of a map
   // std::map<coot::Cartesian, mmdb::Residue *> residue_map;


   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      // write a class that also returns atom specs for the atom. Possibly use inheritance
      // or encapsulation.
      graphical_bonds_container bonds_box = graphics_info_t::molecules[imol].get_bonds_representation();
      r = g.pyobject_from_graphical_bonds_container(imol, bonds_box);
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;

}

//! \brief return a Python object for the bonds
//
PyObject *get_intermediate_atoms_bonds_representation() {

   graphics_info_t g;
   return g.get_intermediate_atoms_bonds_representation();
}




PyObject *get_environment_distances_representation_py(int imol, PyObject *residue_spec_py) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t residue_spec = residue_spec_from_py(residue_spec_py);
      graphics_info_t g;
      graphical_bonds_container gbc = g.molecules[imol].make_environment_bonds_box(residue_spec, g.Geom_p());
      r = g.pyobject_from_graphical_bonds_container(imol, gbc);
   }

   if (PyBool_Check(r)) { 
      Py_INCREF(r);
   }
   return r;

}

//! \brief return a Python object for the radii of the atoms in the dictionary
//
PyObject *get_dictionary_radii() {

   bool use_vdwH_flag = false; // I suppose caller can pass this

   PyObject *r = Py_False;
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

   graphics_info_t g;
   PyObject *residue_dict_py = PyDict_New();
   for (std::size_t i=0; i<g.Geom_p()->size(); i++) {
      const coot::dictionary_residue_restraints_t &rest = g.Geom_p()->get_monomer_restraints(i);
      const std::string residue_name = rest.comp_id();
      const std::vector<coot::dict_atom> &ai = rest.atom_info;
      PyObject *atoms_dict_py = PyDict_New();
      for (std::size_t j=0; j<ai.size(); j++) {
	 const std::string &atom_name = ai[j].atom_id_4c;
	 double r = g.Geom_p()->get_vdw_radius(atom_name, residue_name, imol_enc, use_vdwH_flag);
	 PyDict_SetItem(atoms_dict_py, PyString_FromString(atom_name.c_str()), PyFloat_FromDouble(r));
      }
      PyDict_SetItem(residue_dict_py, PyString_FromString(residue_name.c_str()), atoms_dict_py);
   }
   r = residue_dict_py;
   return r;
}



#endif // Python
