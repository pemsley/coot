/* geometry/make-shelx-restraints.cc
 * 
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

#include <algorithm>
#include <fstream>
#include <iomanip>
#include "utils/coot-utils.hh"
#include "geometry/protein-geometry.hh"

bool is_circular_permuations(const std::vector<std::string> &s_1,
			     std::vector<std::string> s_2) {

   bool v = false;

   if (s_1.size() != s_2.size()) return false;

   bool m_1 = std::equal(s_1.begin(), s_1.end(), s_2.begin());
   if (m_1) return true;

   for (unsigned int i=0; i<s_1.size()-1; i++) {
      std::rotate(s_2.begin(), s_2.begin()+1, s_2.end());
      bool m_r = std::equal(s_1.begin(), s_1.end(), s_2.begin());
      if (m_r) return true;
   }
   return v;
}

void output_flats(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   for (unsigned int i=0; i<rest.plane_restraint.size(); i++) { 
      const coot::dict_plane_restraint_t &pr= rest.plane_restraint[i];
      f << "FLAT_" << rest.residue_info.comp_id;
      for (int iat=0; iat<pr.n_atoms(); iat++)
	 f << " " << coot::util::remove_whitespace(pr.atom_id(iat));
      f << "\n";
   }
}

void output_chivs(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   if (rest.chiral_restraint.size() > 0) { 
      for (unsigned int i=0; i<rest.chiral_restraint.size(); i++) {
	 const coot::dict_chiral_restraint_t &r= rest.chiral_restraint[i];

	 // std::cout << "--------------- r: " << r << std::endl;

	 if (r.is_a_both_restraint())
	    continue;
	 
	 // SHELX alphabetically orders the neighbour atom names
	 // before using the chiral volume - CAs of amino acids have
	 // chiral volume target of about +2.5A^3
	 //
	 // So I need to test if the sorted list is a rotation
	 // permutation of the neighbour list
	 
	 std::vector<std::string> n;

	 if (! rest.is_hydrogen(r.atom_id_1_4c()))
	    n.push_back(coot::util::remove_whitespace(r.atom_id_1_4c()));
	 if (! rest.is_hydrogen(r.atom_id_2_4c()))
	    n.push_back(coot::util::remove_whitespace(r.atom_id_2_4c()));
	 if (! rest.is_hydrogen(r.atom_id_3_4c()))
	    n.push_back(coot::util::remove_whitespace(r.atom_id_3_4c()));


	 if (n.size() != 3) {

	    std::cout << "Oops, need 3-non-hydrogen neighbors for CHIV restraints "
		      << "for atom " << r.atom_id_c_4c() << std::endl;

	 } else {

	    std::vector<std::string> n_original = n;
	    double cv_mult = 1;
	    std::sort(n.begin(), n.end());
	    bool cp = is_circular_permuations(n_original, n);
	    if (! cp) cv_mult = -1;

	    if (r.has_unassigned_chiral_volume()) {
	       std::cout << "debug:: unassigned chiral volume " << r.atom_id_c_4c()
			 << std::endl;
	    } else {
	       double tv = r.target_volume();
	       f << "CHIV_" << rest.residue_info.comp_id;
	       f << " " << tv * cv_mult << " 0.1";
	       f << " " << coot::util::remove_whitespace(r.atom_id_c_4c());
	       for (unsigned int ii=0; ii<3; ii++)
		  f << " " << n[ii];
	       f << "\n";
	    }
	 }
      }
   }
}

void output_dfixs(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   for (unsigned int i=0; i<rest.bond_restraint.size(); i++) { 
      const coot::dict_bond_restraint_t &r = rest.bond_restraint[i];
      f << "DFIX_" << rest.residue_info.comp_id;
      f << " " << std::setprecision(4) << r.value_dist();
      f << " " << coot::util::remove_whitespace(r.atom_id_1());
      f << " " << coot::util::remove_whitespace(r.atom_id_2());
      f << "\n";
   }
}

void output_dangs(const coot::dictionary_residue_restraints_t &rest,
		  std::ofstream &f) {

   // we have to convert angles to 1-3 distances.  Need to find bonds
   // of the angle atoms.
   
   for (unsigned int i=0; i<rest.angle_restraint.size(); i++) { 
      const coot::dict_angle_restraint_t &r = rest.angle_restraint[i];
      if (! rest.is_hydrogen(r.atom_id_1())) { 
	 if (! rest.is_hydrogen(r.atom_id_3())) { 
	    coot::dict_bond_restraint_t b1(r.atom_id_1(), r.atom_id_2(), "", 0.0, 0.0, 0.0, 0.0, false);
	    coot::dict_bond_restraint_t b2(r.atom_id_2(), r.atom_id_3(), "", 0.0, 0.0, 0.0, 0.0, false);
	    for (unsigned int ib1=0; ib1<rest.bond_restraint.size(); ib1++) {
	       if (rest.bond_restraint[ib1].matches_names(b1)) {
		  for (unsigned int ib2=0; ib2<rest.bond_restraint.size(); ib2++) {
		     if (rest.bond_restraint[ib2].matches_names(b2)) {
			double d1 = rest.bond_restraint[ib1].value_dist();
			double d2 = rest.bond_restraint[ib2].value_dist();
			double angle = clipper::Util::d2rad(r.angle());
			double d_sqrd = d1*d1 + d2*d2 - 2*d1*d2*cos(angle);
			if (d_sqrd > 0) {
			   double d = sqrt(d_sqrd);
			   f << "DANG_" << rest.residue_info.comp_id;
			   f << " " << std::setprecision(4) << d;
			   f << " " << coot::util::remove_whitespace(r.atom_id_1());
			   f << " " << coot::util::remove_whitespace(r.atom_id_3());
			   f << "\n";
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}
   

int main(int argc, char **argv) {

   int status = 0;

   if (argc>1) {
      std::string mmcif_file_name(argv[1]);
      coot::protein_geometry geom;
      geom.set_verbose(false);
      int read_number = 0;
      geom.init_refmac_mon_lib(mmcif_file_name, read_number);
      std::vector<std::string> types = geom.monomer_types();

      if (types.size() > 0) {
	 
	 std::string file_name = "shelx-restraints.txt";
	 if (types.size() == 1)
	    file_name = "shelx-restraints-" + types[0] + ".txt";

	 for (unsigned int i=0; i<types.size(); i++) { 
	    std::pair<bool, coot::dictionary_residue_restraints_t> r =
	       geom.get_monomer_restraints(types[0], coot::protein_geometry::IMOL_ENC_ANY);
	    if (r.first) { // how can it not be?

	       if (! r.second.is_filled()) {
		  std::cout << "Input file " << mmcif_file_name
			    << " is not a full description" << std::endl;
	       } else {

		  std::ofstream f(file_name.c_str());
		  const coot::dictionary_residue_restraints_t &rest = r.second;
		  // FLAT CHIV DFIX DANG

		  f << "\n";

		  output_flats(rest, f);
		  output_chivs(rest, f);
		  output_dfixs(rest, f);
		  output_dangs(rest, f);

		  std::cout << "Output " << file_name << std::endl;
	       }
	    }
	 }
      }
   } else {
      std::cout << "Usage: coot-make-shelx-restraints cif-file-name"
		<< std::endl;
   } 
   return status;
}
