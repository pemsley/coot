/* ligand/monomer-utils.cc
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
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

#ifndef MONOMER_UTILS_HH
#define MONOMER_UTILS_HH

#ifndef HAVE_STRING
#include <string>
#endif
#ifndef HAVE_VECTOR
#include <vector>
#endif

#include <mmdb2/mmdb_manager.h>

// no dependency on coords files
// #include "coords/mmdb-extras.h"
// #include "coords/mmdb.h"

#include "clipper/core/coords.h"
#include "mini-mol/atom-quads.hh"

#include "coot-utils/contact-info.hh"

namespace coot {

   // caller needs to delete returned residue and its chain
   // 
   mmdb::Residue *deep_copy_residue(mmdb::Residue *residue);
	 
   // Atom names for a torsion bond, old/Stuart style
   // 
   class atom_name_pair {
   public:
      std::string atom1;
      std::string atom2;
      atom_name_pair(const std::string &atom_name_1,
		     const std::string &atom_name_2) {
	 atom1 = atom_name_1;
	 atom2 = atom_name_2;
      }
   };

   // For looking up atom names to atom indices in the vector of
   // coords that go to the mgtree.
   //
   class atom_index_pair {
   public:
      int index1;
      int index2;
      atom_index_pair() { // unassigned pair
	 index1 = -1;
	 index2 = -1;
      }
      atom_index_pair(int i1, int i2) {
	 index1 = i1;
	 index2 = i2;
      }
   };


   class monomer_utils {

      std::vector<atom_name_pair> atom_name_pair_list;
      std::vector<atom_name_quad> atom_name_quad_list;
      clipper::Coord_orth atom_to_co(mmdb::Atom *at) const;

   public:
      std::vector<atom_name_pair> AtomPairs() const {
	 return atom_name_pair_list;
      } 

      std::vector<atom_name_quad> AtomQuads() const {
	 return atom_name_quad_list;
      }

      void add_torsion_bond_by_name(const std::string &atom_name_1,
				    const std::string &atom_name_2);
      
      void add_torsion_bond_by_name(const std::string &atom_name_1,
				    const std::string &atom_name_2,
				    const std::string &atom_name_3,
				    const std::string &atom_name_4);
      
      std::vector<coot::atom_name_pair>
      atom_name_pairs(const std::string &res_type) const; 

      contact_info getcontacts(const atom_selection_container_t &asc) const; 

      std::vector<atom_index_pair> 
      get_atom_index_pairs(const std::vector<atom_name_pair> &atom_name_pairs,
			   const mmdb::PPAtom atoms, int nresatoms) const;

      std::vector<atom_index_quad> 
      get_atom_index_quads(const std::vector<atom_name_quad> &atom_name_pairs,
			   const mmdb::PPAtom atoms, int nresatoms) const;

      // files in ligand directory cannot depend on coords headers!
      //
      // static Cartesian coord_orth_to_cartesian(const clipper::Coord_orth &c);
      // static clipper::Coord_orth coord_orth_to_cart(const Cartesian &c);
      std::vector<atom_index_quad>
      get_quads(const std::vector<atom_name_quad> &atom_name_quads,
		mmdb::Residue *residue) const;
      std::vector<std::pair<int, float> > get_chi_angles(mmdb::Residue *res) const; // [1-indexed]
   };

   

}

#endif
