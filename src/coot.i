/*
 * src/coot.i
 *
 * Copyright 2007 by University of York
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

%module coot

%{
#include <cstdio>
#include <cstddef>
#include "globjects.h"  //includes gtk/gtk.h
#include "coot-utils/coot-coord-utils.hh"
#include "c-interface.h"
#include "c-interface-preferences.h"
#include "c-interface-refine.hh"
#include "c-interface-database.hh"
#include "c-interface-python.hh"
#include "c-interface-ligands-swig.hh"
#include "c-interface-mogul.hh"
#include "c-interface-sequence.hh"
#include "c-interface-refine.h"
#include "read-molecule.hh"
#include "file-system-utils.hh"
#include "cc-interface.hh"
#include "cc-interface-image.hh"
#include "cc-interface-scripting.hh"
#include "c-interface-refmac.h"
#include "cc-interface-ncs.hh"
#include "c-interface-generic-objects.h"
#include "generic-display-objects-c.h"
#include "manipulation-modes.hh"
#include "network-get.hh"
#include "rotamer-search-modes.hh"
#include "lbg-interface.hh"
#include "sdf-interface.hh"
#include "probe-clash-score.hh"
#include "cmtz-interface.hh"
#include "coot-version.hh"
#include "get-monomer.hh"
#include "globularize.hh"
#include "curlew.hh"
#include "validation.hh"
#include "pythonic-cfc.hh"
#include "rsr-functions.hh"
#include "cc-interface-user-defined-atom-colours.hh"
#include "cc-interface-alignment.hh"
#include "init_coot_as_python_module.hh"
#include "curl-utils.hh"
#include "cfc-2025.hh"
#include "M2T-interface.hh"
#ifdef USE_MOLECULES_TO_TRIANGLES
#include "cc-interface-molecular-representation.hh"
#endif // USE_MOLECULES_TO_TRIANGLES
%}


%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"

namespace std {
%template(vector_string) vector<std::string>;
%template(pairbf) pair<bool, float>;
%template(IntVector) vector<int>;
%template(vector_atom_spec)      vector<coot::atom_spec_t>;
%template(vector_mtz_type_label) vector<coot::mtz_type_label>;
}

%init %{
  init_coot_as_python_module();
%}

// %feature("autodoc", "1"); // add doc string for Intellisense (hopefully)
%feature("autodoc", "1"); // add type hints hopefully

// If the following line is enabled, then the comments extracted
// using -doxygen on the command line are lost.
// %feature("docstring");
// 2025-12-14-PE because, I guess, %feature("docstring") needs to be
// specified for each function!

// 2025-12-15-PE modern method
%include "coot_cc_interface_docs.i"
%include "coot_c_interface_docs.i"
%include "coot_rsr_functions_docs.i"

#include "globjects.h"  //includes gtk/gtk.h
#include "coot-utils/coot-coord-utils.hh"
/* actually we should ignore all GtkWidgets or similar.... */
%ignore main_window();
%ignore main_menubar();
%ignore main_statusbar();
%ignore main_toolbar();
/* conflicts with redefinition */
%ignore list_nomenclature_errors(int);

%include "c-interface.h"
%include "c-interface-preferences.h"
%include "c-interface-refine.hh"
%include "cc-interface-ncs.hh"
%include "cc-interface-molecular-representation.hh"
%include "cc-interface.hh"
%include "read-molecule.hh"
%include "file-system-utils.hh"
%include "cc-interface-image.hh"
%include "cc-interface-scripting.hh"
%include "c-interface-database.hh"
%include "c-interface-python.hh"
%include "c-interface-ligands-swig.hh"
%include "c-interface-mogul.hh"
%include "c-interface-sequence.hh"
%include "c-interface-refmac.h"
%include "c-interface-refine.h"
%include "c-interface-generic-objects.h"
%include "cmtz-interface.hh"
%include "manipulation-modes.hh"
%include "network-get.hh"
%include "rotamer-search-modes.hh"
%include "lbg-interface.hh"
%include "sdf-interface.hh"
%include "geometry/residue-and-atom-specs.hh" // for atom_spec_t
%include "probe-clash-score.hh"
%include "coot-version.hh"
%include "get-monomer.hh"
%include "globularize.hh"
%include "curlew.hh"
%include "validation.hh"
%include "generic-display-objects-c.h"
%include "pythonic-cfc.hh"
%include "rsr-functions.hh"
%include "cc-interface-user-defined-atom-colours.hh"
%include "../coords/phenix-geo.hh"
%include "curl-utils.hh"
%include "M2T-interface.hh"
%include "cc-interface-alignment.hh"
%include "cfc-2025.hh"

// # ignore some functions with c++ objects
%ignore rigid_body_fit_with_residue_range();
