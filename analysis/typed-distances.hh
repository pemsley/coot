/*
 * analysis/typed-distances.hh
 *
 * Copyright 2017 by Medical Research Council
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#ifndef TYPED_DISTANCES_HH
#define TYPED_DISTANCES_HH

// Implementation of ERRAT Colovos & Yeates (1993)

#include <map>
#include <vector>
#include <stdlib.h> // for non-C++-11, remove when this is necessary
#include <mmdb2/mmdb_manager.h>

namespace coot {

   class typed_distances {
   public:
      enum bin_type_t { UNSET= -1, CC, CO, CN, OO, ON, NN};
   private:
      enum atom_type_t { NONE, C, O, N};
      mmdb::Manager *mol;
      // old
      // std::map<mmdb::Residue *, std::vector<std::vector<unsigned int> > > bin_distances;

      // we can't iterate through the bin types bin_type_t.begin(), bin_type_t.end()
      // std::map<mmdb::Residue *, std::map<int, std::vector<float> > > residue_distances_map;
      // so let's index on an int
      //
      std::map<mmdb::Residue *, std::map<int, std::vector<float> > > residue_distances_map;

      // this vector should include self
      std::map<mmdb::Residue *, std::vector<mmdb::Residue *> > residues_within_window;
      void init();
      void setup_bin_distances();
      atom_type_t get_type(mmdb::Atom *at) const;
      // return -1 on no-bin
      int get_atom_pair_bin_id(const atom_type_t &t1, const atom_type_t &t2) const;
      unsigned int n_distance_bins;

      // this is nonsense
      unsigned int get_dist_bin_id(float &dist) const {
	 unsigned int id = 0;
	 if (dist > 4.75) { id = 8; } else {
	    if (dist > 4.5) { id = 7; } else {
	       if (dist > 4.25) { id = 6; } else {
		  if (dist > 4.0) { id = 5; } else {
		     if (dist > 3.75) { id = 4; } else {
			if (dist > 3.5) { id = 3; } else {
			   if (dist > 3.25) { id = 2; } else {
			      if (dist > 3.0) { id = 1; }}}}}}}}
	 return id;
      }
      bool in_self_or_bonded_residue(mmdb::Atom *central_at, mmdb::Atom *neighb_at) const {
	 bool bonded = false;
	 mmdb::Residue *r1 = central_at->GetResidue();
	 mmdb::Residue *r2 = neighb_at->GetResidue();
	 if (! r1) return true;
	 if (! r2) return true;
	 if (r1 == r2) {
	    bonded = true;
	 } else {
	    mmdb::Chain *c1 = r1->GetChain();
	    mmdb::Chain *c2 = r2->GetChain();
	    if (c1 == c2) {
	       int res_no_1 = r1->GetSeqNum();
	       int res_no_2 = r2->GetSeqNum();
#ifdef HAVE_CXX11
	       if (std::abs(res_no_2 - res_no_1) < 2)
		  bonded = true;
#else
	       if (abs(res_no_2 - res_no_1) < 2)
		  bonded = true;
#endif
	    }
	 }
	 
	 return bonded;
      }
      void find_residues_within_window(int wl); // fills residues_within_window
   public:
      typed_distances(mmdb::Manager *mol_in) : mol(mol_in) {
	 n_distance_bins = 9;
	 init();
      }
      void generate(int atom_selection_handle);
      void output() const;
   };

}

#endif // TYPED_DISTANCES_HH
