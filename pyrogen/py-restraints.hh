
#include "Python.h"
#include <string>
#include <iostream>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Bond.h>

#include "protein-geometry.hh"
#include "mogul-interface.hh"


// this coot:: for swig
coot::dictionary_residue_restraints_t monomer_restraints_from_python(PyObject *restraints);

//
namespace coot { 
   PyObject *monomer_restraints_to_python(const dictionary_residue_restraints_t &restraints);
}
