

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <cstddef>

// #include "globjects.h" //includes gtk/gtk.h
#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface-molecular-representation.hh"

#ifdef USE_PYTHON

// Martin's Triangles

// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
      // check that these are strings
      std::string atom_selection = PyString_AsString(atom_selection_py);
      std::string ColorScheme    = PyString_AsString(ColorScheme_py);
      std::string style          = PyString_AsString(style_py);
      status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
   }
   return status;

}

void remove_molecular_represenation(int imol, int rep_no) {

}
#endif // USE_PYTHON
