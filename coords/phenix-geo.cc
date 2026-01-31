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
coot::phenix_geo::operator<<(std::ostream &s, coot::phenix_geo::phenix_geo_bond gb) {

   std::cout << "[phenix geo: " << gb.atom_1 << " " << gb.atom_2 << " ]";
   return s;

}


coot::atom_spec_t
coot::phenix_geo::parse_line_for_atom_spec(const std::string &l) {

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

void
coot::phenix_geo::phenix_geometry::parse(const std::string &file_name) {

   auto parse_for_bond = [] (const std::vector<std::string> &lines, unsigned int line_idx)
      -> std::optional<phenix_geo_bond> {

      const std::string &l = lines[line_idx];
      atom_spec_t atom_spec_1 = parse_line_for_atom_spec(l);
      if (! atom_spec_1.empty()) {
	 const std::string &l2 = lines[line_idx+1];
	 atom_spec_t atom_spec_2 = parse_line_for_atom_spec(l2);
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

   auto parse_for_angle = [] (const std::vector<std::string> &lines, unsigned int line_idx)
      -> std::optional<phenix_geo_angle> {

      if (line_idx + 4 >= lines.size()) return std::nullopt;

      const std::string &l = lines[line_idx];
      atom_spec_t atom_spec_1 = parse_line_for_atom_spec(l);
      if (! atom_spec_1.empty()) {
         const std::string &l2 = lines[line_idx+1];
         atom_spec_t atom_spec_2 = parse_line_for_atom_spec(l2);
         if (! atom_spec_2.empty()) {
            const std::string &l3 = lines[line_idx+2];
            atom_spec_t atom_spec_3 = parse_line_for_atom_spec(l3);
            if (! atom_spec_3.empty()) {
               const std::string &lg = lines[line_idx+4];
               std::vector<std::string> bits = util::split_string_no_blanks(lg);
               if (bits.size() == 6) {
                  try {
                     std::vector<double> v(6,-1);
                     for (unsigned int i=0; i<6; i++) {
                        v[i] = util::string_to_float(bits[i]);
                     }
                     phenix_geo_angle ga(atom_spec_1, atom_spec_2, atom_spec_3);
                     ga.set_geom(v[0], v[1], v[2], v[3], v[4], v[5]);
                     return ga;
                  }
                  catch(const std::runtime_error &rte) {
                     std::cout << "angle number parse fail " << lg << std::endl;
                  }
               }
            }
         }
      }
      return std::nullopt;
   };

   auto parse_for_dihedral = [] (const std::vector<std::string> &lines, unsigned int line_idx)
      -> std::optional<phenix_geo_dihedral> {

      if (line_idx + 5 >= lines.size()) return std::nullopt;

      const std::string &l = lines[line_idx];
      atom_spec_t atom_spec_1 = parse_line_for_atom_spec(l);
      if (! atom_spec_1.empty()) {
         const std::string &l2 = lines[line_idx+1];
         atom_spec_t atom_spec_2 = parse_line_for_atom_spec(l2);
         if (! atom_spec_2.empty()) {
            const std::string &l3 = lines[line_idx+2];
            atom_spec_t atom_spec_3 = parse_line_for_atom_spec(l3);
            if (! atom_spec_3.empty()) {
               const std::string &l4 = lines[line_idx+3];
               atom_spec_t atom_spec_4 = parse_line_for_atom_spec(l4);
               if (! atom_spec_4.empty()) {
                  // Check if sinusoidal or harmonic from header line
                  const std::string &header_line = lines[line_idx+4];
                  bool is_sinusoidal = (header_line.find("sinusoidal") != std::string::npos);

                  const std::string &lg = lines[line_idx+5];
                  std::vector<std::string> bits = util::split_string_no_blanks(lg);
                  if (bits.size() == 7) {
                     try {
                        double ideal = util::string_to_float(bits[0]);
                        double model = util::string_to_float(bits[1]);
                        double delta = util::string_to_float(bits[2]);
                        int periodicity = util::string_to_int(bits[3]);
                        double sigma = util::string_to_float(bits[4]);
                        double weight = util::string_to_float(bits[5]);
                        double residual = util::string_to_float(bits[6]);
                        phenix_geo_dihedral gd(atom_spec_1, atom_spec_2, atom_spec_3, atom_spec_4);
                        gd.set_geom(ideal, model, delta, periodicity, sigma, weight, residual, is_sinusoidal);
                        return gd;
                     }
                     catch(const std::runtime_error &rte) {
                        std::cout << "dihedral number parse fail " << lg << std::endl;
                     }
                  }
               }
            }
         }
      }
      return std::nullopt;
   };

   auto parse_for_chirality = [] (const std::vector<std::string> &lines, unsigned int line_idx)
      -> std::optional<phenix_geo_chiral> {

      if (line_idx + 5 >= lines.size()) return std::nullopt;

      const std::string &l = lines[line_idx];
      atom_spec_t atom_spec_centre = parse_line_for_atom_spec(l);
      if (! atom_spec_centre.empty()) {
         const std::string &l2 = lines[line_idx+1];
         atom_spec_t atom_spec_1 = parse_line_for_atom_spec(l2);
         if (! atom_spec_1.empty()) {
            const std::string &l3 = lines[line_idx+2];
            atom_spec_t atom_spec_2 = parse_line_for_atom_spec(l3);
            if (! atom_spec_2.empty()) {
               const std::string &l4 = lines[line_idx+3];
               atom_spec_t atom_spec_3 = parse_line_for_atom_spec(l4);
               if (! atom_spec_3.empty()) {
                  const std::string &lg = lines[line_idx+5];
                  std::vector<std::string> bits = util::split_string_no_blanks(lg);
                  if (bits.size() == 7) {
                     try {
                        bool both_signs = (bits[0] == "True");
                        double ideal = util::string_to_float(bits[1]);
                        double model = util::string_to_float(bits[2]);
                        double delta = util::string_to_float(bits[3]);
                        double sigma = util::string_to_float(bits[4]);
                        double weight = util::string_to_float(bits[5]);
                        double residual = util::string_to_float(bits[6]);
                        phenix_geo_chiral gc(atom_spec_centre, atom_spec_1, atom_spec_2, atom_spec_3);
                        gc.set_geom(both_signs, ideal, model, delta, sigma, weight, residual);
                        return gc;
                     }
                     catch(const std::runtime_error &rte) {
                        std::cout << "chirality number parse fail " << lg << std::endl;
                     }
                  }
               }
            }
         }
      }
      return std::nullopt;
   };

   auto parse_for_nonbonded = [] (const std::vector<std::string> &lines, unsigned int line_idx)
      -> std::optional<phenix_geo_nonbonded> {

      if (line_idx + 3 >= lines.size()) return std::nullopt;

      const std::string &l = lines[line_idx];
      atom_spec_t atom_spec_1 = parse_line_for_atom_spec(l);
      if (! atom_spec_1.empty()) {
         const std::string &l2 = lines[line_idx+1];
         atom_spec_t atom_spec_2 = parse_line_for_atom_spec(l2);
         if (! atom_spec_2.empty()) {
            const std::string &lg = lines[line_idx+3];
            std::vector<std::string> bits = util::split_string_no_blanks(lg);
            if (bits.size() == 2) {
               try {
                  double model = util::string_to_float(bits[0]);
                  double vdw = util::string_to_float(bits[1]);
                  phenix_geo_nonbonded gnb(atom_spec_1, atom_spec_2);
                  gnb.set_geom(model, vdw);
                  return gnb;
               }
               catch(const std::runtime_error &rte) {
                  std::cout << "nonbonded number parse fail " << lg << std::endl;
               }
            }
         }
      }
      return std::nullopt;
   };


   // ----------------------------- main line ---------------------------
   bool debug = false;
   phenix_geo_bonds bonds;
   phenix_geo_angles angles;
   phenix_geo_dihedrals dihedrals;
   phenix_geo_chirals chirals;
   phenix_geo_nonbondeds nonbondeds;

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
		     if (bond) bonds.add_bond(bond.value());
		  }
		  if (l.substr(0,5) == "angle") {
		     auto angle = parse_for_angle(lines, line_idx);
		     if (angle) angles.add_angle(angle.value());
		  }
		  if (l.substr(0,8) == "dihedral") {
		     auto dihedral = parse_for_dihedral(lines, line_idx);
		     if (dihedral) dihedrals.add_dihedral(dihedral.value());
		  }
		  if (l.substr(0,9) == "chirality") {
		     auto chiral = parse_for_chirality(lines, line_idx);
		     if (chiral) chirals.add_chiral(chiral.value());
		  }
		  if (l.substr(0,9) == "nonbonded") {
		     auto nonbonded = parse_for_nonbonded(lines, line_idx);
		     if (nonbonded) nonbondeds.add_nonbonded(nonbonded.value());
		  }
	       }
	       line_idx++; // next
	    }
	 }
      }
   }

   // Store parsed results in class members
   geo_bonds = bonds;
   geo_angles = angles;
   geo_dihedrals = dihedrals;
   geo_chirals = chirals;
   geo_nonbondeds = nonbondeds;

   if (debug) {
      std::cout << "found " << bonds.size() << " bonds" << std::endl;
      std::cout << "found " << angles.size() << " angles" << std::endl;
      std::cout << "found " << dihedrals.size() << " dihedrals" << std::endl;
      std::cout << "found " << chirals.size() << " chirals" << std::endl;
      std::cout << "found " << nonbondeds.size() << " nonbondeds" << std::endl;

      for (unsigned int ibond=0; ibond<bonds.size(); ibond++) {
	 std::cout << "   " << bonds[ibond].atom_1 << " "<< bonds[ibond].atom_2 << std::endl;
      }
   }
}

//
