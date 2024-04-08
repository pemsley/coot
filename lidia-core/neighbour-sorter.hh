/*
 * lidia-core/neighbour-sorter.hh
 *
 * Copyright 2015 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */



namespace coot {

   class chiral_neighbour_info_t {

   public:
      int idx_mmcif;
      int idx_atom_list;
      mmdb::Atom *at;
      chiral_neighbour_info_t(mmdb::Atom *at_in, int idx_mmcif_in, int idx_atom_list_in) {
	 at = at_in;
	 idx_mmcif = idx_mmcif_in;
	 idx_atom_list = idx_atom_list_in;
      }
      bool operator==(const chiral_neighbour_info_t &test_v) const {
	 if (at == test_v.at)
	    if (idx_mmcif == test_v.idx_mmcif)
	       if (idx_atom_list == test_v.idx_atom_list)
		  return true;
	 return false;
      } 
      static bool neighbour_sorter(const chiral_neighbour_info_t &v1,
				   const chiral_neighbour_info_t &v2);
   };
}
