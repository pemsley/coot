/* coords/phenix-geo.cc
 * 
 * Copyright 2014 by Medical Research Council
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
 * Lesser General Public License for more details.
 * 
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <optional>
#include "utils/coot-utils.hh"
#include "phenix-geo.hh"

std::ostream&
coot::operator<<(std::ostream &s, phenix_geo_bond gb) {

   std::cout << "[phenix geo: " << gb.atom_1 << " " << gb.atom_2 << " ]";
   return s;

}


coot::atom_spec_t
coot::phenix_geo_bonds::parse_line_for_atom_spec(const std::string &l) const {

   atom_spec_t atom_spec;
   std::string atom_name_1(l.substr(10,4));
   std::string chain_id(l.substr(19,1));
   unsigned int post_chain_id_char_idx = 20;

   // Fix multi-char chain-ids here with a while() {}

   unsigned int llen = l.length();
   unsigned int residue_number_str_idx = post_chain_id_char_idx + 1; // start
   std::string residue_number_str;
   std::string ins_code;
   while (residue_number_str_idx < llen) {
      char c = l[residue_number_str_idx];
      if (c >= '0' && c <= '9') {
	 residue_number_str += c;
      } else {
	 if (c != ' ')
	    ins_code = c;
	 if (residue_number_str.length())
	    break;
      }
      residue_number_str_idx++;
   }
   try {
      int res_no = util::string_to_int(residue_number_str);
      atom_spec = atom_spec_t(chain_id, res_no, ins_code, atom_name_1, "");
   }
   catch (const std::runtime_error &rte) {
      // parse fail. Heyho.
   }
   return atom_spec;
}

coot::phenix_geo_bonds::phenix_geo_bonds(const std::string &file_name) {

   auto parse_for_bond = [this] (const std::vector<std::string> &lines, unsigned int line_idx)
      -> std::optional<phenix_geo_bond> {

      const std::string &l = lines[line_idx];
      atom_spec_t atom_spec_1 = this->parse_line_for_atom_spec(l);
      if (! atom_spec_1.empty()) {
	 const std::string &l2 = lines[line_idx+1];
	 atom_spec_t atom_spec_2 = this->parse_line_for_atom_spec(l2);
	 if (! atom_spec_2.empty()) {
	    const std::string &lg = lines[line_idx+3];
	    std::vector<std::string> bits = util::split_string_no_blanks(lg);
	    if (bits.size() == 6) {
	       try {
		  std::vector<double> v(6,-1);
		  for (unsigned int i=0; i<6; i++) {
		     v[i] = util::string_to_float(bits[i]);
		  }
		  phenix_geo_bond gb(atom_spec_1, atom_spec_2);
		  // ideal model delta sigma weight residual
		  gb.set_geom(v[0], v[1], v[2], v[3], v[4], v[5]);
		  return gb;
	       }
	       catch(const std::runtime_error &rte) {
		  // number parse fail.
		  std::cout << "number parse fail " << lg << std::endl;
	       }
	    }
	 }
      }
      return std::nullopt;
   };

   bool debug = false;

   if (! file_exists(file_name)) {
      std::cout << "File not found: " << file_name << std::endl;
   } else {

      std::ifstream f(file_name.c_str());

      if (f) {
	 std::vector<std::string> lines;
	 std::string line;
	 while (std::getline(f, line)) {
	    lines.push_back(line);
	 }
	 if (lines.size() > 0) {
	    unsigned int line_idx = 0;
	    while (line_idx < lines.size()) {
	       const std::string &l = lines[line_idx];
	       if (l.length() > 22) {
		  if (l.substr(0,4) == "bond") {
		     auto bond = parse_for_bond(lines, line_idx);
		     if (bond) bonds.push_back(bond.value());
		  }
		  if (l.substr(0,5) == "angle") {
		  }
		  if (l.substr(0,8) == "dihedral") {
		  }
		  if (l.substr(0,9) == "chirality") {
		  }
		  if (l.substr(0,9) == "chirality") {
		  }
		  if (l.substr(0,9) == "planarity") {
		  }
		  if (l.substr(0,9) == "nonbonded") {
		  }
	       }
	       line_idx++; // next
	    }
	 }
      }
   }

   if (debug) {
      std::cout << "found " << bonds.size() << " bonds" << std::endl;

      for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
	 std::cout << "   " << bonds[ibond].atom_1 << " "<< bonds[ibond].atom_2 << std::endl;
      }
   }
}

//
