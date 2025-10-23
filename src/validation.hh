/*
 * src/validation.hh
 *
 * Copyright 2018 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

// Many more things should go or be moved to here

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#ifdef USE_GUILE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <libguile.h>
#pragma GCC diagnostic pop
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *c_beta_deviations_py(int imol);
#endif


#ifdef USE_GUILE
SCM c_beta_deviations_scm(int imol);
#endif

// atom_spec_list is a 5-member atom spec
void add_unhappy_atom_marker_py(int imol, PyObject *atom_spec_list);
void remove_unhappy_atom_marker_py(int imol, PyObject *atom_spec_list);


