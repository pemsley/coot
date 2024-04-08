/*
 * src/coot-nomenclature.hh
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

#include <mmdb2/mmdb_manager.h>
#ifndef HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR
#include "geometry/protein-geometry.hh"

namespace coot {

   class nomenclature {
      mmdb::Manager *mol_;
      // return swap status.
      int test_and_fix_PHE_TYR_nomenclature_errors(mmdb::Residue *residue_p,
						   bool apply_swap_when_found);

      // 20110725 new function - don't use rotamers, simply use torsions -90 < tor < 90.
      // 
      int test_and_fix_ASP_GLU_nomenclature_errors(mmdb::Residue *residue_p,
						   bool apply_swap_when_found); 

      std::vector<mmdb::Residue *>  fix_and_swap_maybe(protein_geometry *Geom_p,
						  bool apply_swaps); // adjust mol as needed
                                                                     // if flagged.

   public:
      nomenclature(mmdb::Manager *mol) {
	 mol_ = mol;
      }
      // Here we rename atoms to fix nomeclature errors. Note ILEs are not fixed
      // by renaming atoms.
      // 
      std::vector<mmdb::Residue *>  fix(protein_geometry *Geom_p); // adjust mol as needed.

      std::vector<mmdb::Residue *> list(protein_geometry *Geom_p); // Don't touch the molecule,
                                                              // just analyse it.
   }; 
}
