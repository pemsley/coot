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
//! \brief get C-beta deviation information
//!
//! based on geometry defined in the CCP4 monomer library
//!
//! @param imol the molecule index
//! @return a list of per-residue C-beta deviation information, one item
//!         for each residue with a C-beta atom. Each item is a list of 2 items:
//!         [residue-spec, alt-conf-dict], where residue-spec is
//!         [chain-id, res-no, ins-code] and alt-conf-dict maps the
//!         alt-conf (e.g. "" or "A") to the deviation distance in Angstroms.
//!         Return False if imol is not a valid model molecule.
PyObject *c_beta_deviations_py(int imol);
#endif


#ifdef USE_GUILE
SCM c_beta_deviations_scm(int imol);
#endif

#ifdef USE_PYTHON
//! \brief add an "unhappy atom" marker on the given atom
//!
//! The marker is a camera-facing sad-face image drawn at the position
//! of the specified atom, used to flag problem atoms. If imol is not a
//! valid model molecule, or the atom spec is malformed or does not
//! match an atom, this function does nothing.
//!
//! @param imol the molecule index
//! @param atom_spec_list a 5-member atom spec:
//!        [chain-id, res-no, ins-code, atom-name, alt-conf],
//!        e.g. ["A", 145, "", "CE", ""]
void add_unhappy_atom_marker_py(int imol, PyObject *atom_spec_list);

//! \brief remove the "unhappy atom" marker from the given atom
//!
//! The marker at the current position of the specified atom is removed.
//! Note that markers store positions, not atom specs - so if the atom
//! has moved since the marker was added (e.g. by refinement) the marker
//! will not be found (use remove_all_unhappy_atom_markers() in that case).
//!
//! @param imol the molecule index
//! @param atom_spec_list a 5-member atom spec:
//!        [chain-id, res-no, ins-code, atom-name, alt-conf]
void remove_unhappy_atom_marker_py(int imol, PyObject *atom_spec_list);
#endif

//! \brief remove the "unhappy atom" markers from all molecules
void remove_all_unhappy_atom_markers();


