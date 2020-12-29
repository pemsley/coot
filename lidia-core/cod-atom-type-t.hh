/* lidia-core/cod-atom-type-t.hh
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

#ifndef COD_ATOM_TYPE_T_HH
#define COD_ATOM_TYPE_T_HH

#include <list>

#include "third-neighbour-info-t.hh"


namespace cod {

   int hybridization_to_int(RDKit::Atom::HybridizationType h);

   // return the number of rings and the ringinfo string.
   std::pair<int, std::string> make_ring_info_string(const RDKit::Atom *atom_p);

   class atom_level_2_type {
      class atom_level_2_component_type {
      public:
	 std::string element;
	 unsigned int number_of_rings;
	 std::string ring_info_string;
	 std::vector<int> neighb_degrees;
	 std::vector<int> neighb_extra_elect; // same size as above needed
	 std::string atom_name;
	 int n_extra_elect;
	 atom_level_2_component_type(const RDKit::Atom *at, const RDKit::ROMol &rdkm);
	 atom_level_2_component_type() {}
      };
      std::string str;
      std::string element;
      std::vector<atom_level_2_component_type> components;
      int n_extra_elect;
   public:

      atom_level_2_type() {}
      atom_level_2_type(const std::string &s) { str = s;} // read
      atom_level_2_type(const RDKit::Atom *p, const RDKit::ROMol &rdkm);

      std::string string() const {return str; }
      std::string extra_electron_type() const;
      int n_extra_electrons() const;
      static bool level_2_component_sorter(const atom_level_2_component_type &la,
					   const atom_level_2_component_type &lb);
      friend std::ostream &operator<<(std::ostream &s,
				      const atom_level_2_component_type &c);
   };
   std::ostream &operator<<(std::ostream &s,
			    const atom_level_2_type::atom_level_2_component_type &c);
   
   // a container for strings of COD atom types at various levels.
   // The first one was at the 4th level - most sophisticated and
   // with 3rd neighbour info
   //
   class atom_type_t {
      std::string neighb_degrees_str_; // string set on read from tables
      void set_neighb_degrees_string(); //
   public:
      enum { COLON_DEGREE_TYPE = 546 };

      atom_type_t() {}

      atom_type_t(const std::string &s1,
		  const std::string &s_neigh_degrees_colon,
		  const atom_level_2_type &l2,
		  const std::string &s3, const std::string &s4);

      atom_type_t(const std::string &l4) {
	 level_4 = l4;
	 level_3 = level_4_type_to_level_3_type(l4);
	 hash_value = -1;
      }
      atom_type_t(const std::string &s1, const std::string &s2) {
	 level_3 = s1;
	 level_4 = s2;
	 hash_value = -1;
      }
      atom_type_t(const std::vector<unsigned int> &degrees) {
	 neighb_degrees = degrees;
	 set_neighb_degrees_string();
      }
      std::string level_4;
      std::string level_3; // as 4 but without 3rd neighbour info
      atom_level_2_type level_2;
      std::vector<unsigned int> neighb_degrees; // for "colon type" - e.g. 3:3:2
      int hash_value; // can be zero (?), so use -1 for fail.
      std::list<third_neighbour_info_t> tnil;

      std::string neighb_degrees_str() {
	 if (neighb_degrees_str_.empty())
	    set_neighb_degrees_string(); // use neighb_degrees
	 return neighb_degrees_str_;
      }
      std::string neighb_degrees_str() const { return neighb_degrees_str_; }
      void extract_degree_info(const atom_type_t &at) {
	 neighb_degrees = at.neighb_degrees;
	 set_neighb_degrees_string();
      }

      // helper function
      static std::string level_4_type_to_level_3_type(const std::string &l4t);

      bool operator==(const atom_type_t &t_in) const {
	 return (t_in.level_4 == level_4);
      }

      bool l3_match(const atom_type_t &t_in) const {
	 return (t_in.level_3 == level_3);
      }

      bool operator<(const atom_type_t &t_in) const {
	 return (level_4 < t_in.level_4); // is this the right way round?
      }

      void set_hash_value(unsigned int hash_in) {
	 hash_value = hash_in;
      }
   };
}

#endif // COD_ATOM_TYPE_T_HH
