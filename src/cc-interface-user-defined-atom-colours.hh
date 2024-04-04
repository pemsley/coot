/*
 * src/cc-interface-user-defined-atom-colours.hh
 *
 * Copyright 2016 by Medical Research Council
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

#ifdef USE_PYTHON
#include "Python.h"
#endif

#ifdef USE_GUILE
#include <libguile.h>
#endif

// #include<graphics-info.h>

/*  ----------------------------------------------------------------------- */
/*                  user-defined atom colours                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
//! \brief for the given molecule imol, pass a list of tuples of items with the first parameter begin the residue/selection CID
//!        and the second begin the colour index (an int).
void set_user_defined_atom_colour_by_selection_py(int imol, PyObject *residue_specs_colour_index_tuple_list_py);
//! \brief for the given molecule imol, pass a list of tuples of items with the first parameter being the atom spec
//!        and the second being the colour index (an int).
void set_user_defined_atom_colour_py(int imol, PyObject *atom_specs_colour_index_tuple_list_py);
#endif // USE_PYTHON

#ifdef USE_GUILE
void set_user_defined_atom_colour_scm(int imol, SCM atom_specs_colour_index_tuple_list_scm);
void set_user_defined_atom_colour_by_residue_scm(int imol, SCM residue_specs_colour_index_tuple_list_scm);
#endif // USE_GUILE

//! \brief by default user-defined colour indices will use a colour wheel
//  but here we can set them ourselves by passing something like: [1, [0.2, 0.3, 0.4]], [2, [0.6, 0.4, 0.2]]]
//  i.e. a colour for each colour index. If the colour index goes beyond the limit, or there is a
//  missing value, then colour will be [0.5, 0.5, 0.5]
//
#ifdef USE_PYTHON
void set_user_defined_colours_py(PyObject *colour_list_py);
#endif

void clear_user_defined_atom_colours(int imol);

