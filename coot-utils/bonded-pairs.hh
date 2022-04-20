/* coot-utils/bonded-pairs.hh
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
 * Copyright 2015 by Medical Research Council
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

#ifndef HAVE_BONDED_PAIRS_HH
#define HAVE_BONDED_PAIRS_HH

#include <string>
#include <vector>

#include <mmdb2/mmdb_manager.h>

#include "geometry/protein-geometry.hh"

namespace coot {

   // The residues here are in order.  res_1 is comp_1 and res_2 is comp_2
   class bonded_pair_t {
      void delete_atom(mmdb::Residue *res, const std::string &atom_name);
   public:
      mmdb::Residue *res_1;
      mmdb::Residue *res_2;
      std::string link_type;
      bool is_fixed_first;
      bool is_fixed_second;
      bonded_pair_t(mmdb::Residue *r1, mmdb::Residue *r2, bool is_fixed_first_in, bool is_fixed_second_in,
		    const std::string &lt) {
	 res_1 = r1;
	 res_2 = r2;
	 link_type = lt;
	 is_fixed_first = is_fixed_first_in;
	 is_fixed_second = is_fixed_second_in;
      }
      bonded_pair_t() {
	 res_1=0;
	 res_2=0;
	 is_fixed_first  = false;
	 is_fixed_second = false;
      } 
      bool matches(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 if (r1 == res_1 && r2 == res_2) {
	    return true;
	 } else{
	    if (r1 == res_2 && r2 == res_1) {
	       return true;
	    } else {
	       return false;
	    } 
	 }
      }
      bonded_pair_t swapped() const {
	 bonded_pair_t bp;
	 bp.link_type = link_type;
	 bp.res_1 = res_2;
	 bp.res_2 = res_1;
	 bp.is_fixed_first  = is_fixed_second;
	 bp.is_fixed_second = is_fixed_first;
	 return bp;
      }
      void reorder_as_needed();
      // matches?, swap-is-needed-to-match?
      std::pair<bool, bool> matches_info(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 if (r1 == res_1 && r2 == res_2) {
	    return std::pair<bool, bool> (true, false);
	 } else{
	    if (r1 == res_2 && r2 == res_1) {
	       return std::pair<bool, bool> (true, true);
	    } else {
	       return std::pair<bool, bool> (false, false);
	    }
	 }
      }
      void apply_chem_mods(const protein_geometry &geom);
   };
   std::ostream &operator<<(std::ostream &s, bonded_pair_t bp);

   class bonded_pair_match_info_t {
   public:
      bool state;
      bool swap_needed;
      std::string link_type;
      bonded_pair_match_info_t(bool state_in, bool swap_needed_in, const std::string &link_type_in) :
         state(state_in), swap_needed(swap_needed_in), link_type(link_type_in) {}
   };
	 
   class bonded_pair_container_t {
      void reorder();
   public:
      std::vector<bonded_pair_t> bonded_residues;
      // the returned value is "was this pair already added?".  A return value of 0 means that th
      // pair was added.
      bool try_add(const bonded_pair_t &bp); // check for null residues too.
      unsigned int size() const { return bonded_residues.size(); }
      bonded_pair_t operator[](unsigned int i) { return bonded_residues[i]; }
      const bonded_pair_t &operator[](unsigned int i) const { return bonded_residues[i]; }
      bool linked_already_p(mmdb::Residue *r1, mmdb::Residue *r2) const;
      friend std::ostream& operator<<(std::ostream &s, bonded_pair_container_t bpc);
      // test order switch too.
      bool matches(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 bool r = false;
	 for (unsigned int i=0; i<bonded_residues.size(); i++) { 
	    if (bonded_residues[i].matches(r1, r2)) {
	       r = true;
	       break;
	    }
	 }
	 return r;
      }
      bonded_pair_match_info_t match_info(mmdb::Residue *r1, mmdb::Residue *r2) const {
	 bonded_pair_match_info_t mi(false, false, "");
	 for (unsigned int i=0; i<bonded_residues.size(); i++) {
	    std::pair<bool, bool> bb = bonded_residues[i].matches_info(r1, r2);
	    if (bb.first) {
	       mi.link_type = bonded_residues[i].link_type;
	       mi.state = true;
	       if (bb.second)
		  mi.swap_needed = true;
	       break;
	    }
	 }
	 return mi;
      }
      void filter(); // remove residue X 1-3 bonds if residue X 1-2 or 2-3 bonds exist.
      bool closer_exists_p(const bonded_pair_t &bp) const;
      void apply_chem_mods(const protein_geometry &geom);
   };
   std::ostream& operator<<(std::ostream &s, bonded_pair_container_t bpc);

}

#endif // HAVE_BONDED_PAIRS_HH
