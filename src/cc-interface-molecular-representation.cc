/*
 * src/cc-interface-molecular-representation.cc
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
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <cstddef>

// #include "globjects.h" //includes gtk/gtk.h
#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "cc-interface-molecular-representation.hh"

// Martin's Triangles

// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_py(int imol, PyObject *atom_selection_py, PyObject *ColorScheme_py, PyObject *style_py) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
      // check that these are strings
      std::string atom_selection = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_selection_py));
      std::string ColorScheme    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ColorScheme_py));
      std::string style          = PyBytes_AS_STRING(PyUnicode_AsUTF8String(style_py));
      int secondary_structure_usage_flag = CALC_SECONDARY_STRUCTURE;
      // status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
      graphics_info_t g;
      status = g.add_molecular_representation(imol, atom_selection, ColorScheme, style, secondary_structure_usage_flag);
      graphics_draw();
   }
   return status;
}

#ifdef USE_GUILE
// e.g. 0, "//C", "RampChainsScheme", "Ribbon"
int add_molecular_representation_scm(int imol, SCM atom_selection_scm, SCM ColorScheme_scm, SCM style_scm) {

   int status = -1;
   if (is_valid_model_molecule(imol)) {
#ifdef USE_MOLECULES_TO_TRIANGLES
      // check that these are strings
      std::string atom_selection = scm_to_locale_string(atom_selection_scm);
      std::string ColorScheme    = scm_to_locale_string(ColorScheme_scm);
      std::string style          = scm_to_locale_string(style_scm);
      graphics_info_t g;
      int secondary_structure_usage_flag = CALC_SECONDARY_STRUCTURE;
      status = g.add_molecular_representation(imol, atom_selection, ColorScheme, style, secondary_structure_usage_flag);
      graphics_draw();
#endif
   }
   return status;
}
#endif // USE_GUILE

int
add_ribbon_representation_with_user_defined_colours(int imol, const std::string &name) {

   // std::string name = "AlphaFold " + std::to_string(imol)

   int status = -1; // nothing useful

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.add_ribbon_representation_with_user_defined_colours(imol, name);
      graphics_draw();
   }
   return status;
}

void remove_molecular_representation(int imol, int rep_no) {

   if (is_valid_model_molecule(imol)) {
      // graphics_info_t::molecules[imol].remove_molecular_representation(rep_no);
      graphics_info_t g;
      g.remove_molecular_representation(imol, rep_no);
      graphics_draw();
   }
}

extern "C" void add_molecular_representation_test() {
#ifdef USE_MOLECULES_TO_TRIANGLES
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
   if (active_atom.first) {
      int imol = active_atom.second.first;
      std::cout << "Ribbons on molecule " << imol << std::endl;
      if (is_valid_model_molecule(imol)) {
         std::string atom_selection = "//A";
         std::string ColorScheme = "colorRampChainsScheme";
         std::string style = "Ribbon";
         int secondary_structure_usage_flag = CALC_SECONDARY_STRUCTURE;
         // status = graphics_info_t::molecules[imol].add_molecular_representation(atom_selection, ColorScheme, style);
         graphics_info_t g;
         g.add_molecular_representation(imol, atom_selection, ColorScheme, style, secondary_structure_usage_flag);
         graphics_info_t::graphics_draw();
      }
   }
#endif
}

#else

#endif // USE_MOLECULES_TO_TRIANGLES

