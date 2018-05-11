/* lbg/lbg-molfile.cc
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

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

#include "lbg-molfile.hh"

std::ostream&
lig_build::operator<<(std::ostream &s, const molfile_atom_t &at) {

   s << "atom name :" << at.name << ":  ele :" << at.element << ": aromatic? "
     << at.aromatic << " chiral? " << at.chiral << " charge: " << at.formal_charge
     << " at " << at.atom_position.format();
   return s;
}

std::ostream&
lig_build::operator<<(std::ostream &s, const molfile_bond_t &b) {

   s << b.index_1 << " to " << b.index_2 << " type " << b.bond_type;
   return s;
} 



// throw an exception on unable to convert
int
lig_build::string_to_int(const std::string &s) {
   
   int i;
   std::istringstream myStream(s);
   
   if (myStream>>i) { 
      return i;
   } else {
      std::string mess = "Cannot convert \"";
      mess += s;
      mess += "\" to an integer";
      throw std::runtime_error(mess);
   }
} 

// throw an exception on unable to convert
float
lig_build::string_to_float(const std::string &s) {
   
   float f;
   std::istringstream myStream(s);
   
   if (myStream>>f) { 
      return f;
   } else {
      std::string mess = "Cannot convert \"";
      mess += s;
      mess += "\" to an integer";
      throw std::runtime_error(mess);
   }
} 


void 
lig_build::molfile_molecule_t::read(const std::string &file_name) {
   
   std::vector<std::string> header_lines(3);
   int n_chirals = 0;
   int n_atoms = 0;
   int n_bonds = 0;
   int extras_start = 0;
   int extras_end   = 0;
   
   std::vector<std::string> lines;
   std::ifstream in_file(file_name.c_str());
   if (!in_file) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {
      std::string line;
      while (std::getline(in_file, line)) { 
	 lines.push_back(line);
      }
   }

   if (false) // screen output debugging
      for (unsigned int i=0; i<lines.size(); i++)
	 std::cout << ":" << lines[i] << std::endl;

   std::cout << "File \"" << file_name << "\" contained " << lines.size()
	     << " lines " << std::endl;

   if (lines.size() > 3) {
      for (unsigned int i=0; i<3; i++)
	 header_lines[i] = lines[i];

      // count lines:
      //
      if (lines.size() > 3) {
	 int l = lines[3].length();
	 if (l > 15) {
	    try { 
	       std::string n_atoms_string =  lines[3].substr(0,3);
	       std::string n_bonds_string =  lines[3].substr(3,3);
	       n_atoms = lig_build::string_to_int(n_atoms_string);
	       n_bonds = lig_build::string_to_int(n_bonds_string);
	       std::cout << "n_atoms: " << n_atoms << " n_bonds: "
			 << n_bonds << std::endl;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << std::endl;
	    }
	 }

	 extras_start = 4 + n_atoms + n_bonds;
	 extras_end   = lines.size();

	 if (l > 15) {
	    std::string n_chirals_string = lines[3].substr(9,3);
	    try {
	       n_chirals = lig_build::string_to_int(n_chirals_string);
	    }
	    catch (const std::runtime_error &rte) {
	       // std::cout << rte.what() << std::endl; often happens that
	       // "   " cannot be converted to int.
	    }
	 }
	 // std::cout << "n_chirals: " << n_chirals << std::endl;
      }

      // atom block:
      // 
      unsigned int top_lim = n_atoms+4;
      for (unsigned int i=4; i<lines.size()  && i<top_lim; i++) {

	 int l = lines[i].length();
	 if (l > 31) {
	    std::string x_string =  lines[i].substr( 0,9);
	    std::string y_string =  lines[i].substr(10,9);
	    std::string z_string =  lines[i].substr(20,9);
	    std::string ele_str  = lines[i].substr(31,2);
	    try {
	       float x = lig_build::string_to_float(x_string);
	       float y = lig_build::string_to_float(y_string);
	       float z = lig_build::string_to_float(z_string);
	       std::string ele;
	       ele += ele_str[0];
	       if (ele_str[1] != ' ')
		  ele = ele_str;
	       lig_build::molfile_atom_t at(x,y,z,ele);
	       // std::cout << "ele :" << ele << ":" << std::endl;
	       atoms.push_back(at);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << std::endl; 
	    }
	 }
      }

      // bond block
      unsigned int n_bonds_and_n_atoms_lim = n_bonds+n_atoms+4;
      for (unsigned int i=top_lim; i<lines.size() && i<n_bonds_and_n_atoms_lim; i++) {
	 // std::cout << "parse this bond line line :" << lines[i] << std::endl;
	 int l = lines[i].length();
	 if (l > 8) {
	    std::string index_1_string =   lines[i].substr( 0,3);
	    std::string index_2_string =   lines[i].substr( 3,3);
	    std::string bond_type_string = lines[i].substr( 6,3);
	    try {

	       int stereo_bond = 0;
	       if (l> 11) {
		  std::string stereo_bond_string = lines[i].substr( 9,3);
		  stereo_bond = lig_build::string_to_int(stereo_bond_string);
	       }
	       
	       int index_1   = lig_build::string_to_int(index_1_string);
	       int index_2   = lig_build::string_to_int(index_2_string);
	       int bond_type = lig_build::string_to_int(bond_type_string);
	       lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
	       if (bond_type == 2)
		  bt =  lig_build::bond_t::DOUBLE_BOND;
	       if (bond_type == 3)
		  bt =  lig_build::bond_t::TRIPLE_BOND;
	       if (((index_1-1) >=0) && (index_1-1) < int(atoms.size())) { 
		  if (((index_2-1) >=0) && (index_2-1) < int(atoms.size())) { 
		     lig_build::molfile_bond_t bond(index_1-1, index_2-1, bt);
		     if (stereo_bond != 0) {
			// something interesting
			// std::cout << "atoms " << index_1 << " " << index_2
			// << " stereo_bond: " << stereo_bond << std::endl;
			if (stereo_bond == 6)
			   bond.bond_type = lig_build::bond_t::IN_BOND;
			if (stereo_bond == 1)
			   bond.bond_type = lig_build::bond_t::OUT_BOND;
		     }
		     if (0) 
			std::cout << "Added bond with type " << bond.bond_type << " "
				  << stereo_bond << std::endl;
		     bonds.push_back(bond);
		  }
	       }
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << std::endl; 
	    }
	 } 
      }

      // PROPERTY BLOCK
      // 
      // extras: charges, radicals, isotopes, ring-bond count and wotnot...
      //
      for (int i= extras_start; i<extras_end; i++) {
	 // std::cout << "Check for extras: " << lines[i] << std::endl;
	 if (lines[i].length() > 8) {
	    std::string extra_type = lines[i].substr(3,3);

	    // handle CHARGES
	    // 
	    if (extra_type == "CHG") {
	       // std::cout << "    parse charge on " << lines[i] << std::endl;
	       std::string n_charged_atoms_str = lines[i].substr(7,3);
	       try { 
		  int nca = lig_build::string_to_int(n_charged_atoms_str);
		  // std::cout << "trying to find " << nca << " charged atoms on line"
		  // << std::endl;
		  
		  for (int ic=0; ic<nca; ic++) {
		     if (int(lines[i].length()) > 10+ic*6+3) { 
			std::string atom_number_string = lines[i].substr(10+ic*8,3);
			std::string charge_string      = lines[i].substr(14+ic*8,3);
			// std::cout << "found atom_number_string :"
			// << atom_number_string << ":"
			//  << " charge_string :"
			// << charge_string
			// << ":" << std::endl;
			int atom_number = lig_build::string_to_int(atom_number_string);
			int charge      = lig_build::string_to_int(charge_string);
			if (atom_number < n_atoms) {
			   if (atom_number > 0) { 
			      atoms[atom_number-1].formal_charge = charge;
			   }
			} 
		     }
		  } 
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << rte.what() << std::endl;
	       } 
	    } 
	 } 
      } 
   }
}



// Extra addtion, so that we can make a molfile_molecule_t from
// an MMDB molecule and restraints - this is so that we can
// bring topological filtering to the results of PRODRG.
//
// We can only handle atoms in the residue that have altconf of "" - others are ignored.
//
lig_build::molfile_molecule_t::molfile_molecule_t(mmdb::Residue *residue_p,
						  const coot::dictionary_residue_restraints_t &restraints) {

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   molfile_atom_t blank_atom(0,0,0, "");
   atoms.push_back(blank_atom); // blank atom for 0-index.
   // make the atoms
   for (int iat=0; iat<n_residue_atoms; iat++) {
      clipper::Coord_orth pos(residue_atoms[iat]->x,
			      residue_atoms[iat]->y,
			      residue_atoms[iat]->z);
      molfile_atom_t atom(pos,
			  residue_atoms[iat]->element,
			  residue_atoms[iat]->name);
      atoms.push_back(atom);
   } 
   
   
   std::map<std::string, std::vector<mmdb::Atom *> > atom_map;
   std::map<std::string, std::vector<mmdb::Atom *> >::const_iterator it_1_atom_map;
   std::map<std::string, std::vector<mmdb::Atom *> >::const_iterator it_2_atom_map;
   std::map<mmdb::Atom *, int> atom_index_map;
   for (int iat=0; iat<n_residue_atoms; iat++) { 
      atom_map[residue_atoms[iat]->name].push_back(residue_atoms[iat]);
      atom_index_map[residue_atoms[iat]] = iat;
   }

   
   for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
      const coot::dict_bond_restraint_t &bond_restraint = restraints.bond_restraint[ibond];
      std::string atom_name_1 = bond_restraint.atom_id_1_4c();
      std::string atom_name_2 = bond_restraint.atom_id_2_4c();
      it_1_atom_map = atom_map.find(atom_name_1);
      it_2_atom_map = atom_map.find(atom_name_2);
      if (it_1_atom_map != atom_map.end()) { 
	 if (it_2_atom_map != atom_map.end()) {
	    const std::vector<mmdb::Atom *> &v_1 = it_1_atom_map->second;
	    const std::vector<mmdb::Atom *> &v_2 = it_2_atom_map->second;
	    for (unsigned int iat_1=0; iat_1<v_1.size(); iat_1++) {
	       std::string alt_conf_1 = v_1[iat_1]->altLoc;
	       if (alt_conf_1 == "") {
		  for (unsigned int iat_2=0; iat_2<v_2.size(); iat_2++) { 
		     std::string alt_conf_2 = v_2[iat_2]->altLoc;
		     if (alt_conf_2 == "") {

			// OK, so we have a bond, these are indices
			// into the atoms of the residue (they need to
			// be offset to match the atoms of the
			// molfile_molecule
			// 
			int idx_1 = atom_index_map[v_1[iat_1]];
			int idx_2 = atom_index_map[v_2[iat_2]];

			bond_t::bond_type_t bond_type = get_bond_type(bond_restraint.type());

			molfile_bond_t bond(idx_1+1, idx_2+1, bond_type);
			bonds.push_back(bond);
		     }
		  } 
	       }
	    }
	 }
      }
   }
}

// Extra addtion, so that we can make a molfile_molecule_t from
// an MMDB molecule and restraints - this is so that we can
// bring topological filtering to the results of PRODRG.
//
// We can only handle atoms in the residue that have altconf of "" - others are ignored.
//
lig_build::molfile_molecule_t::molfile_molecule_t(const coot::dictionary_residue_restraints_t &restraints) {

   std::map<std::string, int> atom_name_index;
   std::map<std::string, int>::const_iterator it_1_atom_name_index;
   std::map<std::string, int>::const_iterator it_2_atom_name_index;
   
   molfile_atom_t blank_atom(0,0,0, "");
   atoms.push_back(blank_atom); // blank atom for 0-index.
   for (unsigned int iat=0; iat<restraints.atom_info.size(); iat++) { 
      lig_build::molfile_atom_t atom(clipper::Coord_orth(0,0,0),
				     restraints.atom_info[iat].type_symbol,
				     restraints.atom_info[iat].atom_id);
      atoms.push_back(atom);
      atom_name_index[restraints.atom_info[iat].atom_id]=iat;
   }

   for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
      const coot::dict_bond_restraint_t &bond_restraint = restraints.bond_restraint[ibond];
      it_1_atom_name_index = atom_name_index.find(bond_restraint.atom_id_1_4c());
      it_2_atom_name_index = atom_name_index.find(bond_restraint.atom_id_2_4c());
      if (it_1_atom_name_index != atom_name_index.end()) {
	 if (it_2_atom_name_index != atom_name_index.end()) {
	    int idx_1 = it_1_atom_name_index->second;
	    int idx_2 = it_2_atom_name_index->second;
	    bond_t::bond_type_t bond_type = get_bond_type(bond_restraint.type());
	    molfile_bond_t bond(idx_1+1, idx_2+1, bond_type);
	    bonds.push_back(bond);
	 }
      }
   } 
} 

lig_build::bond_t::bond_type_t
lig_build::molfile_molecule_t::get_bond_type(const std::string &bond_restraint_type) const {

   bond_t::bond_type_t bond_type = lig_build::bond_t::BOND_UNDEFINED;
   if (bond_restraint_type == "single")
      bond_type = lig_build::bond_t::SINGLE_BOND;
   if (bond_restraint_type == "double")
      bond_type = lig_build::bond_t::DOUBLE_BOND;
   if (bond_restraint_type == "triple")
      bond_type = lig_build::bond_t::TRIPLE_BOND;
   if (bond_restraint_type == "triple")
      bond_type = lig_build::bond_t::TRIPLE_BOND;
   if (bond_restraint_type == "aromatic")
      bond_type = lig_build::bond_t::AROMATIC_BOND;

   return bond_type;
}
