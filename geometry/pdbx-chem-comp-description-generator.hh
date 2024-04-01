/*
 * geometry/pdbx-chem-comp-description-generator.hh
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
#ifndef PDBX_CHEM_COMP_DESCRIPTOR_HH
#define PDBX_CHEM_COMP_DESCRIPTOR_HH

#include <string>

namespace coot {

   class pdbx_chem_comp_description_generator_t {
   public:
      pdbx_chem_comp_description_generator_t() {}
      pdbx_chem_comp_description_generator_t(const std::string &pn, const std::string &pv, const std::string &d) :
         program_name(pn), program_version(pv), descriptor(d) {}
      std::string program_name;
      std::string program_version;
      std::string descriptor;
   };

}

#endif // PDBX_CHEM_COMP_DESCRIPTOR_HH
