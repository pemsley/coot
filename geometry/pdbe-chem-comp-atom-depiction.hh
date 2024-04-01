/*
 * geometry/pdbe-chem-comp-atom-depiction.hh
 *
 * Copyright 2023 by Medical Research Council
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
#ifndef PDBE_CHEM_COMP_ATOM_DEPICTION_HH
#define PDBE_CHEM_COMP_ATOM_DEPICTION_HH

#include <vector>
#include <string>

namespace coot {

   class depiction_atom_t {
   public:
      depiction_atom_t(const std::string &atom_name, const std::string &element, double xx, double yy, int idx) :
         atom_name(atom_name), element(element), x(xx), y(yy), pdbe_ordinal(idx) {}
      std::string atom_name;
      std::string element;
      double x, y;
      int pdbe_ordinal;
   };

   class chem_comp_atom_depiction_t {
      std::string comp_id;
   public:
      chem_comp_atom_depiction_t() : comp_id("unset") {}
      chem_comp_atom_depiction_t(const std::string &comp_id_in, const std::vector<depiction_atom_t> &atoms_in) :
         comp_id(comp_id_in), atoms(atoms_in) {}
      std::vector<depiction_atom_t> atoms;
      bool empty() const { return atoms.empty(); }
   };

}

#endif // PDBE_CHEM_COMP_ATOM_DEPICTION_HH
