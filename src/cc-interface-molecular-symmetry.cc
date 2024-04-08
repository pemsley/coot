/*
 * src/cc-interface-molecular-symmetry.cc
 *
 * Copyright 2021 by Medical Research Council
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


// don't forget to git add this file

#ifdef USE_PYTHON  // fixes Python include file order problems
#include "Python.h"
#endif

#include "cc-interface.hh"
#include "c-interface.h"

#include "graphics-info.h"

void add_molecular_symmetry(int imol,
                            double r_00, double r_01, double r_02,
                            double r_10, double r_11, double r_12,
                            double r_20, double r_21, double r_22,
                            double about_origin_x,
                            double about_origin_y,
                            double about_origin_z) {

   if (is_valid_model_molecule(imol)) {
      clipper::Coord_orth molecule_origin(about_origin_x, about_origin_y, about_origin_z);
      clipper::Mat33<double> mol_symm_matrix(r_00, r_01, r_02,
                                             r_10, r_11, r_12,
                                             r_20, r_21, r_22);
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      m.add_molecular_symmetry(mol_symm_matrix, molecule_origin);
      graphics_draw();
   }

}


int add_molecular_symmetry_from_mtrix_from_self_file(int imol) {

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].add_molecular_symmetry_from_mtrix_from_self_file();
   }
   graphics_draw();
   return istat;
}


int add_molecular_symmetry_from_mtrix_from_file(int imol, const std::string &file_name) {

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].add_molecular_symmetry_from_mtrix_from_file(file_name);
   }
   graphics_draw();
   return istat;
}
