/*
 * src/cc-interface-merge.cc
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

// --------------------------------------------
#ifdef USE_PYTHON
#include <Python.h> // add first for _XOPEN_SOURCE order issues
#endif

#include <cstddef> // define std::ptrdiff_t

#ifdef USE_GUILE
#include <libguile.h>
#endif

#include "graphics-info.h"

// --------------------------------------------
#include "c-interface.h" // for is_valid_model_molecule()
#include "cc-interface.hh"

/*  ------------------------------------------------------------------------ */
/*                         merge fragments                                   */
/*  ------------------------------------------------------------------------ */
//! \name Merge Fragments
//! \{
//! \brief merge fragments
// 
//! each fragment is presumed to be in its own chain.
//
int merge_fragments(int imol) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = graphics_info_t::molecules[imol].merge_fragments();
      graphics_draw();
      graphics_info_t g;
      g.update_validation(imol);
   }

   return status;

}
//! \}
