/* lidia-core/cod-types.hh
 * 
 * Copyright 2016 by Medical Research Council
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#ifndef COD_TYPES_HH
#define COD_TYPES_HH

#include <string>
#include "use-rdkit.hh"

#include "bond-table-record-t.hh"
#include "bond-record-container-t.hh"

namespace cod {

   // we need info for neighbours (count) and aromatic status
   class ring_info_t {
      std::vector<int> atom_indices; // mirrors the atomRings() data (type) in RDKit
      bool aromaticity;
   public:
      ring_info_t(const std::vector<int> &ai) { atom_indices = ai; aromaticity = false; }
      ring_info_t(const std::vector<int> &ai, bool aromaticity_in)
      { atom_indices = ai; aromaticity = aromaticity_in; }
      void add_atom(unsigned int i) { atom_indices.push_back(i); }
      void set_aromaticity(bool arom) { aromaticity = arom; }
      unsigned int size() const { return atom_indices.size(); }
      bool get_aromaticity() const { return aromaticity; }
   };

   bond_record_container_t read_acedrg_table(const std::string &file_name);
   
} 

#endif // COD_TYPES_HH

#endif // MAKE_ENHANCED_LIGAND_TOOLS
