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

   for (unsigned int i=0; i<lines.size(); i++) { 
      std::cout << ":" << lines[i] << std::endl;
   }

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
	       std::cout << "n_atoms: " << n_atoms << " n_bonds: " << n_bonds << std::endl;
	    }
	    catch (std::runtime_error rte) {
	       std::cout << rte.what() << std::endl;
	    }
	 }

	 if (l > 15) {
	    std::string n_chirals_string = lines[3].substr(9,3);
	    try {
	       n_chirals = lig_build::string_to_int(n_chirals_string);
	    }
	    catch (std::runtime_error rte) {
	       // std::cout << rte.what() << std::endl; often happens that
	       // "   " cannot be converted to int.
	    }
	 }
	 std::cout << "n_chirals: " << n_chirals << std::endl;
      }

      // atom block:
      //
      for (unsigned int i=4; i<lines.size()  && i<(n_atoms+4); i++) {
	 // std::cout << "parse this atom block :" << lines[i] << std::endl;
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
	       // std::cout << "read x y z " << x << " " << y << " " << z << std::endl;
	       std::string ele;
	       ele += ele_str[0];
	       if (ele_str[1] != ' ')
		  ele = ele_str;
	       lig_build::molfile_atom_t at(x,y,z,ele);
	       // std::cout << "ele :" << ele << ":" << std::endl;
	       atoms.push_back(at);
	    }
	    catch (std::runtime_error rte) {
	       std::cout << rte.what() << std::endl; 
	    }
	 }
      }

      // bond block
      for (unsigned int i=(4+n_atoms); i<lines.size()  && i<(n_bonds+n_atoms+4); i++) {
	 // std::cout << "parse this bond line line :" << lines[i] << std::endl;
	 int l = lines[i].length();
	 if (l > 8) {
	    std::string index_1_string =   lines[i].substr( 0,3);
	    std::string index_2_string =   lines[i].substr( 3,3);
	    std::string bond_type_string = lines[i].substr( 6,3);
	    try {
	       int index_1   = lig_build::string_to_int(index_1_string);
	       int index_2   = lig_build::string_to_int(index_2_string);
	       int bond_type = lig_build::string_to_int(bond_type_string);
	       lig_build::bond_t::bond_type_t bt = lig_build::bond_t::SINGLE_BOND;
	       if (bond_type == 2)
		  bt =  lig_build::bond_t::DOUBLE_BOND;
	       if (bond_type == 3)
		  bt =  lig_build::bond_t::TRIPLE_BOND;
	       if (((index_1-1) >=0) && (index_1-1) < atoms.size()) { 
		  if (((index_2-1) >=0) && (index_2-1) < atoms.size()) { 
		     lig_build::molfile_bond_t bond(index_1-1, index_2-1, bt);
		     bonds.push_back(bond);
		  }
	       }
	    }
	    catch (std::runtime_error rte) {
	       std::cout << rte.what() << std::endl; 
	    }
	 } 
      }
   }
}

	     
