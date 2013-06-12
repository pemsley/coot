/* db-main/db-strand.hh
 * 
 * Copyright 2007 The University of York
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


#include "mini-mol.hh"

namespace coot {

   class db_strands {
      std::string ref_str_dir_str;
      std::vector<std::string> get_reference_pdb_list() const;
      CMMDBManager * get_mol(const std::string &filename) const;
      std::vector<coot::minimol::molecule>
      strand_analysis(CModel *model_p, CMMDBManager *mol,
		      const std::string &filename, int strand_length) const;
      std::pair<bool, clipper::RTop_orth>
      orient_strand_on_z(int SelHnd, CMMDBManager *mol) const;
      void apply_rtop_to_strand(int SelHnd, CMMDBManager *mol,
				const clipper::RTop_orth &rtop) const;
      std::vector<clipper::Coord_orth> z_control_points(int nres) const;
      void trim_to_mainchain(CMMDBManager *mol) const;
      
   public:
      db_strands();
      db_strands(const std::string &dir_ref_structs) {
	 ref_str_dir_str = dir_ref_structs;
      } 
      std::vector<minimol::molecule> get_reference_strands(int n_strands, int strand_length);
   };
} 
