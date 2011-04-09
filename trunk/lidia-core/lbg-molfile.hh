/* lbg/lbg-molefile.hh
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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

#ifndef LBG_MOLFILE_HH
#define LBG_MOLFILE_HH

#include "clipper/core/coords.h"

#include "mmdb_manager.h" // 20110408 for the construtor using restraints and CResidue
#include "protein-geometry.hh"

#include "lig-build.hh"


namespace lig_build {

   // throw an exception on unable to convert
   int string_to_int(const std::string &s);
   // throw an exception on unable to convert
   float string_to_float(const std::string &s);
   

   class molfile_atom_t {
   public:
      clipper::Coord_orth atom_position;
      std::string name;
      std::string element;
      bool aromatic;
      bool chiral;
      molfile_atom_t(const clipper::Coord_orth &pos_in,
		     const std::string &element_in,
		     const std::string &name_in) {
	 atom_position = pos_in;
	 element = element_in;
	 name = name_in;
      }
      molfile_atom_t(float x, float y, float z, std::string ele_in) {
	 atom_position = clipper::Coord_orth(x,y,z);
	 element = ele_in;
	 name = ele_in;
      }
      friend std::ostream &operator<<(std::ostream &s, const molfile_atom_t &a);
   };
   std::ostream &operator<<(std::ostream &s, const molfile_atom_t &a);

   class molfile_bond_t {
   public:
      int index_1;
      int index_2;
      bond_t::bond_type_t bond_type;
      molfile_bond_t(int i1, int i2, bond_t::bond_type_t bt) {
	 index_1 = i1;
	 index_2 = i2;
	 bond_type = bt;
      } 
   };

   class molfile_molecule_t {
      bond_t::bond_type_t get_bond_type(const std::string &restraints_bond_type) const;
   public:
      molfile_molecule_t() {}
      std::vector<molfile_atom_t> atoms;
      std::vector<molfile_bond_t> bonds;
      void read(const std::string &file_name);
      void add_atom(const molfile_atom_t &at) {
	 atoms.push_back(at);
      }
      void add_bond(const molfile_bond_t &b) {
	 bonds.push_back(b);
      }

      // Extra addtion, so that we can make a molfile_molecule_t from
      // an MMDB molecule and restraints - this is so that we can
      // bring topological filtering to the results of PRODRG.
      // 
      // We can (only) handle atoms in the residue that have altconf of "" -
      // others are ignored.
      //
      // This is a bit of bone-headedness on my part - this is not
      // needed for what I want - I do not need the atom position.
      // 
      molfile_molecule_t(CResidue *residue_p, const coot::dictionary_residue_restraints_t &restraints);

      // This is what I want - make fake atom positions and return a
      // list of chiral atoms (i.e. tetravalent and have
      // non-equivalent bonding atoms)
      // 
      molfile_molecule_t(const coot::dictionary_residue_restraints_t &restraints);
      
   };

   
}

#endif // LBG_MOLFILE_HH

