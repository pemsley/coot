/* pyrogen/py-restraints.hh
 * 
 * Copyright 2011 by the University of Oxford
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

#include "Python.h"
#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Bond.h>

#include <geometry/protein-geometry.hh>
#include <analysis/mogul-interface.hh>


// this coot:: for swig
coot::dictionary_residue_restraints_t monomer_restraints_from_python(PyObject *restraints);

//
namespace coot {
   PyObject *monomer_restraints_to_python(const dictionary_residue_restraints_t &restraints);
}

