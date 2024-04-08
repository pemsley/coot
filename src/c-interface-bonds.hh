/*
 * src/c-interface-bonds.hh
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

#ifndef C_INTERFACE_BONDS_HH
#define C_INTERFACE_BONDS_HH

// change the name of this - it's not just bonds

// c-interface-unexported or some such

#include <clipper/core/coords.h>

// not to be exported to the API
#ifdef USE_PYTHON

PyObject *go_to_ligand_py();
clipper::Coord_orth go_to_ligand_inner();

#endif // USE_PYTHON

#endif // C_INTERFACE_BONDS_HH
