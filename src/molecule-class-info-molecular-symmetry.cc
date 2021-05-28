
/* src/molecule-class-info.cc
 *
 * Copyright 2021 by Medical Research Council
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

#include "molecule-class-info.h"

void
molecule_class_info_t::add_molecular_symmetry(const clipper::Mat33<double> &mol_symm,
                                              const clipper::Coord_orth &mol_origin) {

   // std::cout << "#### pushing back mol_origin " << mol_origin.format() << std::endl;
   molecular_symmetry_matrices.push_back(std::make_pair(mol_symm, mol_origin));

}


void
molecule_class_info_t::add_molecular_symmetry_from_mtrix_from_self_file() {

   if (has_model()) {
      std::string file_name = name_;
      if (coot::file_exists(file_name)) {
         add_molecular_symmetry_from_mtrix_from_file(file_name);
      }
   }
}

void
molecule_class_info_t::add_molecular_symmetry_from_mtrix_from_file(const std::string &file_name) {

   std::vector<clipper::RTop_orth> mv = coot::mtrix_info(file_name);
   for (unsigned int i=0; i<mv.size(); i++) {
      const clipper::RTop_orth &rt = mv[i];
      const clipper::Mat33<double> &mat = rt.rot();
      clipper::Coord_orth trn_2(rt.trn());
      clipper::Coord_orth trn(0.5 * trn_2);
      add_molecular_symmetry(mat, trn);
   }
}

