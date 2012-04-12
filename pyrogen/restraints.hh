
#include "Python.h"
#include <string>
#include <iostream>
#include <vector>

#include "protein-geometry.hh"
#include "mogul-interface.hh"

coot::dictionary_residue_restraints_t
monomer_restraints_py(const char *monomer_type, PyObject *restraints);

namespace coot { 
   mogul make_a_mogul(PyObject* pyo);
} 
