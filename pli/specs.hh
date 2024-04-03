/*
 * pli/specs.hh
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

#ifndef PY_SPECS_HH
#define PY_SPECS_HH

#include <Python.h>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class py_atom_spec_t : public atom_spec_t {
   public:
      py_atom_spec_t(PyObject *obj);
      py_atom_spec_t(const atom_spec_t &spec_in) : atom_spec_t(spec_in) {}
      PyObject *pyobject() const;
   };

   class py_residue_spec_t : public residue_spec_t {
   public:
      py_residue_spec_t(PyObject *obj);
      py_residue_spec_t(const residue_spec_t &spec_in) : residue_spec_t(spec_in) {}
      PyObject *pyobject() const;
   };

}

#endif // PY_SPECS_HH
