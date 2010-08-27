/* src/c-interface-scm.hh
 * 
 * Copyright 2008 by The University of York
 * Author: Bernhard Lohkamp
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

#ifndef C_INTERFACE_PYTHON_HH
#define C_INTERFACE_PYTHON_HH

#ifdef USE_PYTHON

// Load the head if it hasn't been included.
#ifndef PYTHONH
#include <Python.h>
#endif

#include "coot-coord-utils.hh"
// This is a common denominator really.  It does not depend on mmdb,
// but it can't be declared in c-interface.h because then we'd have to
// include c-interface.h which would cause (resolvable, I think, not
// checked) problems.
// 
// return a python string, decode to c++ using scm_to_locale_string();
PyObject *display_python(PyObject *o);

std::pair<bool, coot::atom_spec_t> make_atom_spec_py(PyObject *spec);
std::pair<bool, coot::residue_spec_t> make_residue_spec_py(PyObject *spec);

#endif  // USE_PYTHON

#endif // C_INTERFACE_PYTHON_HH
