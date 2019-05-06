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

#include <algorithm>

#include <clipper/core/coords.h>

#include <mmdb2/mmdb_manager.h> // 20110408 for the constructor using restraints and mmdb::Residue
#include "geometry/protein-geometry.hh"

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
      int chiral; // encode direction information, CW CCW
      int formal_charge; 
      molfile_atom_t(const clipper::Coord_orth &pos_in,
		     const std::string &element_in,
		     const std::string &name_in) {
	 atom_position = pos_in;
	 element = element_in;
	 name = name_in;
	 formal_charge = 0;
	 aromatic = false;
	 chiral = 0;
      }
      molfile_atom_t(float x, float y, float z, std::string ele_in) {
	 atom_position = clipper::Coord_orth(x,y,z);
	 element = ele_in;
	 name = ele_in;
	 formal_charge = 0;
	 aromatic = false;
	 chiral = 0;
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
      friend std::ostream &operator<<(std::ostream &s, const molfile_bond_t &a);
   };
   std::ostream &operator<<(std::ostream &s, const molfile_bond_t &a);

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
      molfile_molecule_t(mmdb::Residue *residue_p,
			 const coot::dictionary_residue_restraints_t &restraints);

      // This is what I want - make fake atom positions and return a
      // list of chiral atoms (i.e. tetravalent and have
      // non-equivalent bonding atoms)
      // 
      molfile_molecule_t(const coot::dictionary_residue_restraints_t &restraints);

      std::pair<bool, double> get_scale_correction() const {

	 bool status = 0;
	 double scale = 1.0;
	 std::vector<double> bond_lengths; // process this to scale up the mol file molecule if
	 // needed.
   
	 for (unsigned int i=0; i<bonds.size(); i++) {
	    int index_1 = bonds[i].index_1;
	    int index_2 = bonds[i].index_2;
	    if (atoms[index_1].element != "H") { 
	       if (atoms[index_2].element != "H") { 
		  double l =
		     clipper::Coord_orth::length(atoms[index_1].atom_position,
						 atoms[index_2].atom_position);
		  bond_lengths.push_back(l);
	       }
	    }
	 }
	 if (bond_lengths.size() > 0) {
	    // for (unsigned int ibond=0; ibond<bond_lengths.size(); ibond++)
	    // std::cout << "  bond length " << ibond << " " << bond_lengths[ibond] << std::endl;
	    
	    status = 1;
	    std::sort(bond_lengths.begin(), bond_lengths.end());
	    int index = bond_lengths.size()/2;
	    double bll = bond_lengths[index];
	    scale = 1.0/bll;
	 }
	 return std::pair<bool, double> (status, scale);
      }

      bool delete_bond_between(int idx_1, int idx_2) {
	 bool status = false;
	 std::vector<molfile_bond_t>::iterator it;
	 for (it=bonds.begin(); it!=bonds.end(); it++) { 
	    int index_1 = it->index_1;
	    int index_2 = it->index_2;
	    if (idx_1 == index_1) { 
	       if (idx_2 == index_2) {
		  bonds.erase(it);
		  status = true;
		  break;
	       }
	    }
	 }
	 return status;
      } 

      void debug() const {
	 std::cout << "molfile_molecule_t: " << atoms.size() << " atoms" << std::endl;
	 std::cout << "molfile_molecule_t: " << bonds.size() << " bonds" << std::endl;
	 for (unsigned int iat=0; iat<atoms.size(); iat++) 
	    std::cout << "   " << iat << " " << atoms[iat] << std::endl;
	 for (unsigned int ib=0; ib<bonds.size(); ib++) 
	    std::cout << "   bond " << ib << ":  " << bonds[ib] << std::endl;

	 // bond length info
	 std::vector<double> bond_lengths;
	 for (unsigned int i=0; i<bonds.size(); i++) {
	    int index_1 = bonds[i].index_1;
	    int index_2 = bonds[i].index_2;
	    if (atoms[index_1].element != "H") { 
	       if (atoms[index_2].element != "H") { 
		  double l =
		     clipper::Coord_orth::length(atoms[index_1].atom_position,
						 atoms[index_2].atom_position);
		  bond_lengths.push_back(l);
	       }
	    }
	 }
	 for (unsigned int ibond=0; ibond<bond_lengths.size(); ibond++)
	    std::cout << "  bond length " << ibond << " " << bond_lengths[ibond] << std::endl;

      } 
   };
}

#endif // LBG_MOLFILE_HH

