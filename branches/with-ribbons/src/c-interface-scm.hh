/* src/c-interface-scm.hh
 * 
 * Copyright 2007 by The University of Oxford
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

#ifndef C_INTERFACE_SCM_HH
#define C_INTERFACE_SCM_HH

#ifdef USE_GUILE

// Load the head if it hasn't been included.
#ifndef LIBGUILEH
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#endif

#include "coot-coord-utils.hh"
// This is a common denominator really.  It does not depend on mmdb,
// but it can't be declared in c-interface.h because then we'd have to
// include c-interface.h which would cause (resolvable, I think, not
// checked) problems.
// 
// return a scm string, decode to c++ using scm_to_locale_string();
SCM display_scm(SCM o);

bool scm_is_undefined(SCM o); // ?

#define DIRECT_SCM_STRING ";; # DIRECT SCHEME"

std::pair<bool, coot::atom_spec_t> make_atom_spec(SCM spec);
std::pair<bool, coot::residue_spec_t> make_residue_spec(SCM spec);

#endif  // USE_GUILE

#endif // C_INTERFACE_SCM_HH
