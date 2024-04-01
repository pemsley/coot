/*
 * ideal/new-linked-residue-t.hh
 * 
 * Copyright 2020 by Medical Research Council
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


#ifndef NEW_LINKED_RESIDUE_T_HH
#define NEW_LINKED_RESIDUE_T_HH

#include <string>
#include <vector>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   class new_linked_residue_t {
   public:
      new_linked_residue_t(mmdb::Atom *at_1, mmdb::Atom *at_2,
			   mmdb::Residue *res_1_in, mmdb::Residue *res_2_in,
			   bool is_fixed_first_in, bool is_fixed_second_in,
			   const std::string &link_type,
			   bool order_switch_flag_in) : at_1(at_1), at_2(at_2), res_1(res_1_in), res_2(res_2_in),
                                                        link_type(link_type), order_switch_flag(order_switch_flag_in) {

	 is_fixed_first  = is_fixed_first_in;
	 is_fixed_second = is_fixed_second_in;

	 if (order_switch_flag) {
	    is_fixed_first = is_fixed_second_in;
	    is_fixed_second = is_fixed_first_in;
	    std::swap(res_1, res_2);
	    order_switch_flag = false;
	 }
      }
      mmdb::Atom *at_1, *at_2;
      mmdb::Residue *res_1, *res_2;
      bool is_fixed_first;
      bool is_fixed_second;
      std::string link_type;
      bool order_switch_flag;
      bool match_p(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 if (r1 == res_1 && r2 == res_2) {
	    return true;
	 } else {
	    return (r2 == res_1 && r1 == res_2);
	 }
      }
   };

   class new_linked_residue_list_t {
   public:
      new_linked_residue_list_t() {}
      std::vector<new_linked_residue_t> nlr_vec;
      void insert(mmdb::Atom *at_1, mmdb::Atom *at_2,
		  mmdb::Residue *res_1, mmdb::Residue *res_2,
		  bool is_fixed_first_in, bool is_fixed_second_in,
		  const std::string &link_type,
		  bool order_switch_flag) {
	 bool found = already_added_p(res_1, res_2);
	 if (! found) {
	    new_linked_residue_t nlr(at_1, at_2, res_1, res_2, is_fixed_first_in, is_fixed_second_in,
				     link_type, order_switch_flag);
	    nlr_vec.push_back(nlr);
	 }
      }
      bool already_added_p(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 bool found = false;
	 for (std::size_t i=0; i<nlr_vec.size(); i++) {
	    if (nlr_vec[i].match_p(r1, r2)) {
	       found = true;
	       break;
	    }
	 }
	 return found;
      }
      std::size_t size() const { return nlr_vec.size(); }
      bool empty() { return nlr_vec.empty(); }
      const new_linked_residue_t &operator[](std::size_t i) { return nlr_vec[i]; }
   };
}

#endif // NEW_LINED_RESIDUE_T_HH
