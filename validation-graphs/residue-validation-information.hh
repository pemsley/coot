/*
 * validation-graphs/residue-validation-information.hh
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
#ifndef RESIDUE_VALIDATION_INFORMATION_HH
#define RESIDUE_VALIDATION_INFORMATION_HH

#include <string>
#include "../geometry/residue-and-atom-specs.hh"


namespace coot {
   /// Represents graph data for a single residue in validation graphs
   class residue_validation_information_t {
   public:
      residue_validation_information_t(const coot::residue_spec_t &rs,
                                       const coot::atom_spec_t &atom_spec_in,
                                       double d, const std::string &l) :
      residue_spec(rs), atom_spec(atom_spec_in), function_value(d), label(l) {}
      residue_spec_t residue_spec;
      atom_spec_t atom_spec;
      double function_value;
      std::string label;
   };
}

#endif // RESIDUE_VALIDATION_INFORMATION_HH
