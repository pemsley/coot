/*
 * src/cc-interface-molecular-representation.hh
 *
 * Copyright 2017 by Medical Research Council
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


#ifdef USE_MOLECULES_TO_TRIANGLES
#ifdef USE_PYTHON
#include <Python.h>

#include <string>

// Martin's Triangles

//! \brief add a molecular representation, e.g. Ribbon diagram or surface
//!
//! e.g. 0, "//C", "RampChainsScheme", "Ribbon"
//!
//! @param imol the molecule index
//! @param atom_selection_py the atom-selection in mmdb-style CID
//! @param ColorScheme_py colour-scheme-name, e.g. "colorRampChainsScheme" or "Chain" or "colorBySecondaryScheme"
//! @param style_py the representation style, e.g. "Ribbon" or "MolecularSurface"
//! @return -1 on failure, the representation index on success
int add_molecular_representation_py(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py);
#endif // USE_PYTHON

#ifdef USE_GUILE
int add_molecular_representation_scm(int imol, SCM atom_selection_scm, SCM ColorScheme_scm, SCM style_scm);
#endif // USE_GUILE

//! \brief add ribbon representation with user-defined colours
//!
//! not dependent on scm or python
//!
//! @param imol the molecule index
//! @param name the name to appear in the display control widget
//! @return will be -1 - the return value currently has no meaning.
//!
int add_ribbon_representation_with_user_defined_colours(int imol, const std::string &name);

//! \brief remove molecular representation
//!
//! @param imol the molecule index
//! @param rep_no the representation-index
void remove_molecular_representation(int imol, int rep_no);

#endif // USE_MOLECULES_TO_TRIANGLES
