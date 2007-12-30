/* ligand/monomer-utils.cc
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb_manager.h"

#include "clipper/core/coords.h"

namespace coot {

   CResidue *deep_copy_residue(const CResidue *residue);
	 
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

   // 4 atom names for a torsion:
   class atom_name_quad {
   public:
      std::string atom1;
      std::string atom2;
      std::string atom3;
      std::string atom4;
      atom_name_quad(const std::string &atom_name_1,
		     const std::string &atom_name_2,
		     const std::string &atom_name_3,
		     const std::string &atom_name_4) {
	 atom1 = atom_name_1;
	 atom2 = atom_name_2;
	 atom3 = atom_name_3;
	 atom4 = atom_name_4;
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

   // For looking up atom names to atom indices in the vector of
   // coords that go to the mgtree.
   //
   class atom_index_quad {
   public:
      int index1;
      int index2;
      int index3;
      int index4;
      atom_index_quad() { // unassigned pair
	 index1 = -1;
	 index2 = -1;
	 index3 = -1;
	 index4 = -1;
      }
      atom_index_quad(int i1, int i2, int i3, int i4) {
	 index1 = i1;
	 index2 = i2;
	 index3 = i3;
	 index4 = i4;
      }
   };


   // Now in mmdb-extras
//    class contact_info {
//    public:
//       contact_info(PSContact con_in, int nc) {
// 	 pscontact = con_in;
// 	 n_contacts = nc;
//       }
//       PSContact pscontact;
//       int n_contacts;
//    };
 

   class monomer_utils {

      std::vector<atom_name_pair> atom_name_pair_list;
      std::vector<atom_name_quad> atom_name_quad_list;
      clipper::Coord_orth atom_to_co(CAtom *at) const;

   public:
      std::vector<atom_name_pair> AtomPairs() const {
	 return atom_name_pair_list;
      }

      void add_torsion_bond_by_name(const std::string &atom_name_1,
				    const std::string &atom_name_2);
      
      void add_torsion_bond_by_name(const std::string &atom_name_1,
				    const std::string &atom_name_2,
				    const std::string &atom_name_3,
				    const std::string &atom_name_4);
      
      std::vector<coot::atom_name_pair>
      atom_name_pairs(const std::string &res_type) const; 

      coot::contact_info getcontacts(const atom_selection_container_t &asc) const; 

      std::vector<coot::atom_index_pair> 
      get_atom_index_pairs(const std::vector<coot::atom_name_pair> &atom_name_pairs,
			   const PPCAtom atoms, int nresatoms) const;

      std::vector<coot::atom_index_quad> 
      get_atom_index_quads(const std::vector<coot::atom_name_quad> &atom_name_pairs,
			   const PPCAtom atoms, int nresatoms) const;

      static Cartesian coord_orth_to_cartesian(const clipper::Coord_orth &c);
      static clipper::Coord_orth coord_orth_to_cart(const Cartesian &c);
      std::vector<coot::atom_index_quad>
      get_quads(const std::vector<coot::atom_name_quad> &atom_name_quads,
		CResidue *residue) const;
      std::vector<std::pair<int, float> > get_chi_angles(CResidue *res) const; // [1-indexed]
   };

   

}

#endif
