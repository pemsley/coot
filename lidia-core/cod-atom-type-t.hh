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

   class nb1nb2_type {
      class nb1nb2_component_type {
      public:
	 std::string element;
	 unsigned int number_of_rings;
	 std::string ring_info_string;
	 std::vector<int> neighb_degrees;
	 std::vector<int> neighb_extra_elect; // same size as above needed
	 std::string atom_name;
	 int n_extra_elect;
	 nb1nb2_component_type(const RDKit::Atom *at, const RDKit::ROMol &rdkm);
	 nb1nb2_component_type() {}
      };
      std::string str;
      std::string element;
      int n_extra_elect;
   public:

      nb1nb2_type() {}
      explicit nb1nb2_type(const std::string &s) : str(s) { n_extra_elect = 0; } // read
      nb1nb2_type(const RDKit::Atom *p, const RDKit::ROMol &rdkm);

      std::vector<nb1nb2_component_type> components;
      std::string string() const {return str; }
      std::string extra_electron_type() const;
      int n_extra_electrons() const;
      static bool nb1nb2_component_sorter(const nb1nb2_component_type &la,
					   const nb1nb2_component_type &lb);
      friend std::ostream &operator<<(std::ostream &s,
				      const nb1nb2_component_type &c);
   };
   std::ostream &operator<<(std::ostream &s,
			    const nb1nb2_type::nb1nb2_component_type &c);

   // a container for strings of COD atom types at various levels.
   // The first one was at the 4th level - most sophisticated and
   // with 3rd neighbour info
   //
   class atom_type_t {
      std::string nb2_extra_els_str_; // string set on read from tables
      void set_nb2_extra_els_string(); //
   public:
      enum { COLON_DEGREE_TYPE = 546 };

      atom_type_t() {}

      atom_type_t(const std::string &s1,
		  const std::string &s_neigh_degrees_colon,
		  const nb1nb2_type &l2,
		  const std::string &s3, const std::string &s4);

      atom_type_t(const std::string &l) {
	 full_type = l;
	 main_type = cod_type_to_main_type(l);
	 hash_value = -1;
      }
      //! called with main_type and full_type (not the other way around!)
      atom_type_t(const std::string &s1, const std::string &s2) {
	 main_type = s1;
	 full_type = s2;
	 hash_value = -1;
      }
      explicit atom_type_t(const std::vector<unsigned int> &degrees) : nb2_extra_els(degrees) {
	 set_nb2_extra_els_string();
      }
      std::string full_type;                   // acedrg FULL class, incl. trailing {...}
      std::string main_type;                  // as cod_type but without 3rd neighbour info
      nb1nb2_type nb1nb2;                      // acedrg NB1NB2
      std::string sp;      // acedrg sp string (SP1/SP2/SP3/SP-NON/SPDn), == libmol Atom*_sp
      std::string element; // acedrg Atom*_elem
      std::vector<unsigned int> nb2_extra_els; // for "colon type" - e.g. 3:3:2
      int hash_value; // can be zero (?), so use -1 for fail.
      std::list<third_neighbour_info_t> tnil;

      std::string nb2_extra_els_str() {
	 if (nb2_extra_els_str_.empty())
	    set_nb2_extra_els_string(); // use nb2_extra_els
	 return nb2_extra_els_str_;
      }
      std::string nb2_extra_els_str() const { return nb2_extra_els_str_; }
      void extract_degree_info(const atom_type_t &at) {
	 nb2_extra_els = at.nb2_extra_els;
	 set_nb2_extra_els_string();
      }

      // helper function
      static std::string cod_type_to_main_type(const std::string &l4t);

      bool operator==(const atom_type_t &t_in) const {
	 return (t_in.full_type == full_type);
      }

      bool main_type_match(const atom_type_t &t_in) const {
	 return (t_in.main_type == main_type);
      }

      bool operator<(const atom_type_t &t_in) const {
	 return (full_type < t_in.full_type); // is this the right way round?
      }

      void set_hash_value(unsigned int hash_in) {
	 hash_value = hash_in;
      }
   };
}

#endif // COD_ATOM_TYPE_T_HH
