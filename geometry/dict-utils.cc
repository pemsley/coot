/* geometry/dict-utils.cc
 * 
 * Copyright 2014, 2015 by Medical Research Council
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

#include <fstream>
#include <map>
#include <algorithm>
#include <iomanip> // setw()
#include "utils/coot-utils.hh"
#include "protein-geometry.hh"
#include "dict-mismatches.hh"
#include "dict-utils.hh"

// quiet means don't tell me about matches
bool
coot::dictionary_residue_restraints_t::compare(const dictionary_residue_restraints_t &r,
					       double bond_length_tolerance,
					       double bond_esd_tolerance,
					       double angle_tolerance,
					       double angle_esd_tolerance,
					       bool compare_hydrogens,
					       bool output_energy_types,
					       bool quiet) const {

   std::cout << "Info:: using bond_length_tolerance " << bond_length_tolerance << " A" << std::endl;
   std::cout << "Info:: using angle_tolerance       " << angle_tolerance << " degrees " << std::endl;

   std::map<std::string, bool> hydrogen_status;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) {
      bool h_status = is_hydrogen(atom_info[iat].atom_id_4c);
      hydrogen_status[atom_info[iat].atom_id_4c] = h_status;
      hydrogen_status[atom_info[iat].atom_id   ] = h_status;
   }
   
   // residue info
   bool residue_info_matches = true;
   if (r.residue_info.comp_id != residue_info.comp_id) {
      // this can't happen, surely.
      std::cout << "Residue-Info:: mismatch comp_id \""
		<< residue_info.comp_id << "\" vs \""
		<< r.residue_info.comp_id << "\"" << std::endl;
      residue_info_matches = false;
   }
   std::string comp_id_s;
   if (output_energy_types)
      comp_id_s = std::string(" ") + residue_info.comp_id + std::string(" ");
   
   if (r.residue_info.three_letter_code != residue_info.three_letter_code) {
      std::cout << "Residue-Info:: " << comp_id_s << " mismatch three_letter_code "
		<< residue_info.three_letter_code << " vs "
		<< r.residue_info.three_letter_code
		<< std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.name != residue_info.name) {
      std::cout << "Residue-Info:: " << comp_id_s << " mismatch name \""
		<< residue_info.name << "\" vs \"" << r.residue_info.name << "\""
		<< std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.group != residue_info.group) {
      std::cout << "Residue-Info:: " << comp_id_s << " mismatch group "
		<< residue_info.group << " vs " << r.residue_info.group << std::endl;
      residue_info_matches = false;
   }
   if (compare_hydrogens) { 
      if (r.residue_info.number_atoms_all != residue_info.number_atoms_all) {
	 std::cout << "Residue-Info:: " << comp_id_s << " mismatch number_atoms_all "
		   << residue_info.number_atoms_all << " " << r.residue_info.number_atoms_all
		   << std::endl;
	 residue_info_matches = false;
      }
   }
   if (r.residue_info.number_atoms_nh != residue_info.number_atoms_nh) {
      std::cout << "Residue-Info:: " << comp_id_s << " mismatch number_atoms_nh "
		<< residue_info.number_atoms_nh << " " << r.residue_info.number_atoms_nh
		<< std::endl;
      residue_info_matches = false;
   }
   if (residue_info_matches)
      if (! quiet)
	 std::cout << "Residue-Info:: " << comp_id_s << "     all residue attributes match " << std::endl;

   // atom info
   if (atom_info.size() != r.atom_info.size()) {
      std::cout << "Atom-Info:: " << comp_id_s << " mismatch number of atoms " << atom_info.size() << " vs "
	 << r.atom_info.size() << std::endl;
   }


   // check for atoms that are in this atom_info, but not in r.
   //
   // check for atoms that are in r, but not in this atom info.
   //
   std::vector<std::string> missing_atoms;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) { 
      if (compare_hydrogens || ! is_hydrogen(atom_info[iat].atom_id_4c)) { 
	 const std::string &atom_id_refr = atom_info[iat].atom_id;
	 bool found = false;
	 for (unsigned int jat=0; jat<r.atom_info.size(); jat++) {
	    if (r.atom_info[jat].atom_id == atom_id_refr) {
	       found = true;
	       break;
	    } 
	 }
	 if (! found)
	    missing_atoms.push_back(atom_id_refr);
      }
   }
   if (missing_atoms.size()) {
      std::string s = "s";
      if (missing_atoms.size() == 1) s = "";
      std::cout << "Atom-Info:: atom" << s << " in dict-1, but not dict-2: ";
      for (unsigned int i=0; i<missing_atoms.size(); i++)
	 std::cout << missing_atoms[i] << " ";
      std::cout << std::endl;
   }
   missing_atoms.clear();
   
   for (unsigned int iat=0; iat<r.atom_info.size(); iat++) { 
      if (compare_hydrogens || ! r.is_hydrogen(r.atom_info[iat].atom_id_4c)) { 
	 const std::string &atom_id_refr = r.atom_info[iat].atom_id;
	 bool found = false;
	 for (unsigned int jat=0; jat<atom_info.size(); jat++) {
	    if (atom_info[jat].atom_id == atom_id_refr) {
	       found = true;
	       break;
	    } 
	 }
	 if (! found)
	    missing_atoms.push_back(atom_id_refr);
      }
   }
   if (missing_atoms.size()) {
      std::string s = "s";
      if (missing_atoms.size() == 1) s = "";
      std::cout << "Atom-Info:: atom" << s << " in dict-2, but not dict-1: ";
      for (unsigned int i=0; i<missing_atoms.size(); i++)
	 std::cout << missing_atoms[i] << " ";
      std::cout << std::endl;
   }

   
   for (unsigned int iat=0; iat<atom_info.size(); iat++) { 
      const std::string &atom_id_refr = atom_info[iat].atom_id;
      if (0) 
	 std::cout << "atom id ref: " << atom_id_refr << " hydrogen status "
		   << is_hydrogen(atom_info[iat].atom_id_4c) << std::endl;
      if (compare_hydrogens || ! is_hydrogen(atom_info[iat].atom_id_4c)) { 
	 for (unsigned int jat=0; jat<r.atom_info.size(); jat++) { 
	    const std::string &atom_id_test = r.atom_info[jat].atom_id;
	    if (atom_id_test == atom_id_refr) {
	       if (atom_info[iat].type_symbol != r.atom_info[jat].type_symbol) {
		  std::cout << "Atom-Info:: " << comp_id_s << " " << atom_id_refr << " type_symbol mismatch "
			    << atom_info[iat].type_symbol << " vs "
			    << r.atom_info[jat].type_symbol << " "
			    << std::endl;
	       }
	       if (atom_info[iat].type_energy != r.atom_info[jat].type_energy) {
		  std::cout << "Atom-Info:: " << comp_id_s << " " << atom_info[iat].atom_id_4c
			    << " type_energy mismatch "
			    << atom_info[iat].type_energy << " vs "
			    << r.atom_info[jat].type_energy << " "
			    << std::endl;
	       }
	    }
	 }
      }
   }

   // -------------------------  bond restraints -------------------
   //
   std::vector<bond_mismatch_t> bond_length_mismatches;
   std::vector<bond_mismatch_t> bond_esd_mismatches;
   if (bond_restraint.size() != r.bond_restraint.size())
      std::cout << "Bond-Restraint:: " << comp_id_s << " mismatch number of restraints "
		<< bond_restraint.size() << " " << r.bond_restraint.size()
		<< std::endl;
   bool bonds_match = true;
   for (unsigned int ib=0; ib<bond_restraint.size(); ib++) { 
      for (unsigned int jb=0; jb<r.bond_restraint.size(); jb++) {
	 if (bond_restraint[ib].matches_names(r.bond_restraint[jb])) {
	    bond_mismatch_t b(bond_restraint[ib].atom_id_1(),
			      bond_restraint[ib].atom_id_2(),
			      bond_restraint[ib].value_dist(),
			      r.bond_restraint[jb].value_dist());
	    bond_mismatch_t e(bond_restraint[ib].atom_id_1(),
			      bond_restraint[ib].atom_id_2(),
			      bond_restraint[ib].value_esd(),
			      r.bond_restraint[jb].value_esd());
	    bond_length_mismatches.push_back(b);
	    bond_esd_mismatches.push_back(e);

	    std::string type_1 = bond_restraint[ib].type().substr(0,4);
	    std::string type_2 = r.bond_restraint[jb].type().substr(0,4);
	    std::string up_1 = util::upcase(type_1);
	    std::string up_2 = util::upcase(type_2);
	    if (up_1 != up_2) {
	       std::cout << "Bond-Restraint:: " << comp_id_s << " mismatch bond order  between "
			 << bond_restraint[ib].atom_id_1() << " "
			 << bond_restraint[ib].atom_id_2() << " "
			 << bond_restraint[ib].type() << " vs "
			 << r.bond_restraint[jb].type()
			 << std::endl;
	    }
	 }
      }
   }
   std::sort(bond_length_mismatches.begin(), bond_length_mismatches.end());
   std::sort(bond_esd_mismatches.begin(), bond_esd_mismatches.end());
   std::cout.precision(3);
   std::cout << std::fixed;
   std::cout.width(6);
   std::cout << std::internal;
   
   for (unsigned int ib=0; ib<bond_length_mismatches.size(); ib++) {
      if (bond_length_mismatches[ib].abs_diff >= bond_length_tolerance) {
	 bool at_1_is_hydrogen = false;
	 bool at_2_is_hydrogen = false;
	 std::map<std::string, bool>::const_iterator it_1, it_2;
	 it_1 = hydrogen_status.find(bond_length_mismatches[ib].atom_id_1);
	 it_2 = hydrogen_status.find(bond_length_mismatches[ib].atom_id_2);
	 if (it_1 != hydrogen_status.end()) at_1_is_hydrogen = it_1->second;
	 if (it_2 != hydrogen_status.end()) at_2_is_hydrogen = it_2->second;
	 
	 if (compare_hydrogens || (!at_1_is_hydrogen && !at_2_is_hydrogen)) { 
	    std::cout << "Bond-Restraint:: " << comp_id_s << " mismatch bond length between     "
		      << bond_length_mismatches[ib].atom_id_1 << " "
		      << bond_length_mismatches[ib].atom_id_2 << "  "
		      << bond_length_mismatches[ib].dist_1 << " vs "
		      << bond_length_mismatches[ib].dist_2 << "  delta: ";
	    std::cout.width(6);
	    std::cout << std::right << bond_length_mismatches[ib].diff;
	    if (output_energy_types)
	       std::cout << " " << std::left << std::setw(4)
			 << r.type_energy(bond_length_mismatches[ib].atom_id_1)
			 << " " << std::setw(4)
			 << r.type_energy(bond_length_mismatches[ib].atom_id_2);
	    std::cout << "\n";
	    bonds_match = false;
	 }
      }
   }
   for (unsigned int ib=0; ib<bond_esd_mismatches.size(); ib++) {
      if (bond_esd_mismatches[ib].abs_diff > bond_esd_tolerance) { 
	 bool at_1_is_hydrogen = false;
	 bool at_2_is_hydrogen = false;
	 std::map<std::string, bool>::const_iterator it_1, it_2;
	 it_1 = hydrogen_status.find(bond_esd_mismatches[ib].atom_id_1);
	 it_2 = hydrogen_status.find(bond_esd_mismatches[ib].atom_id_2);
	 if (it_1 != hydrogen_status.end()) at_1_is_hydrogen = it_1->second;
	 if (it_2 != hydrogen_status.end()) at_2_is_hydrogen = it_2->second;
	 if (compare_hydrogens || (!at_1_is_hydrogen && !at_2_is_hydrogen)) { 
	    std::cout << "Bond-Restraint:: " << comp_id_s << " mismatch bond length esd between "
		      << bond_esd_mismatches[ib].atom_id_1 << " "
		      << bond_esd_mismatches[ib].atom_id_2 << "  "
		      << bond_esd_mismatches[ib].dist_1 << " vs "
		      << bond_esd_mismatches[ib].dist_2 << "  delta: ";
	    std::cout.width(6);
	    std::cout << bond_esd_mismatches[ib].diff  << "\n";
	    bonds_match = false;
	 }
      }
   }
   if (bonds_match)
      if (! quiet)
	 std::cout << "Bond-Restraint:: " << comp_id_s << "   all bonds  match within tolerance "
		   << std::endl;

   // -------------------------  angle restraints -------------------
   // 
   if (angle_restraint.size() != r.angle_restraint.size())
      std::cout << "Angle-Restraint:: " << comp_id_s << " mismatch number of restraints "
		<< angle_restraint.size() << " " << r.angle_restraint.size()
		<< std::endl;
   std::vector<angle_mismatch_t> angle_mismatches;
   std::vector<angle_mismatch_t> angle_esd_mismatches;
   bool angles_match = true;
   for (unsigned int ia=0; ia<angle_restraint.size(); ia++) { 
      for (unsigned int ja=0; ja<r.angle_restraint.size(); ja++) {
	 if (angle_restraint[ia].matches_names(r.angle_restraint[ja])) {
	    angle_mismatch_t a(angle_restraint[ia].atom_id_1(),
			       angle_restraint[ia].atom_id_2(),
			       angle_restraint[ia].atom_id_3(),
			       angle_restraint[ia].angle(),
			       r.angle_restraint[ja].angle());
	    angle_mismatch_t e(angle_restraint[ia].atom_id_1(),
			       angle_restraint[ia].atom_id_2(),
			       angle_restraint[ia].atom_id_3(),
			       angle_restraint[ia].esd(),
			       r.angle_restraint[ja].esd());
	    angle_mismatches.push_back(a);
	    angle_esd_mismatches.push_back(e);
	 }
      }
   }

   std::sort(angle_mismatches.begin(), angle_mismatches.end());
   std::sort(angle_esd_mismatches.begin(), angle_esd_mismatches.end());

   
   for (unsigned int ia=0; ia<angle_mismatches.size(); ia++) {
      if (angle_mismatches[ia].abs_diff >= angle_tolerance) {

	 bool at_1_is_hydrogen = false;
	 bool at_2_is_hydrogen = false;
	 bool at_3_is_hydrogen = false;
	 std::map<std::string, bool>::const_iterator it_1, it_2, it_3;
	 it_1 = hydrogen_status.find(angle_mismatches[ia].atom_id_1);
	 it_2 = hydrogen_status.find(angle_mismatches[ia].atom_id_2);
	 it_3 = hydrogen_status.find(angle_mismatches[ia].atom_id_3);
	 if (it_1 != hydrogen_status.end()) at_1_is_hydrogen = it_1->second;
	 if (it_2 != hydrogen_status.end()) at_2_is_hydrogen = it_2->second;
	 if (it_3 != hydrogen_status.end()) at_3_is_hydrogen = it_3->second;
	    
	 if (compare_hydrogens || (!at_1_is_hydrogen && !at_2_is_hydrogen && !at_3_is_hydrogen)) { 
	 
	    std::cout << "Angle-Restraint:: " << comp_id_s << " mismatch angle between     "
		      << angle_mismatches[ia].atom_id_1 << " "
		      << angle_mismatches[ia].atom_id_2 << "  "
		      << angle_mismatches[ia].atom_id_3 << "  "
		      << angle_mismatches[ia].angle_1 << " vs "
		      << angle_mismatches[ia].angle_2 << "  delta: ";
	    std::cout.width(6);
	    std::cout << std::right << angle_mismatches[ia].diff;
	    if (output_energy_types)
	       std::cout << " " << std::left << std::setw(4)
			 << r.type_energy(angle_mismatches[ia].atom_id_1)
			 << " " << std::setw(4)
			 << r.type_energy(angle_mismatches[ia].atom_id_2)
			 << " " << std::setw(4)
			 << r.type_energy(angle_mismatches[ia].atom_id_3);
	    std::cout << "\n";
	    angles_match = false;
	 }
      }
   }

   for (unsigned int ia=0; ia<angle_esd_mismatches.size(); ia++) {
      if (angle_esd_mismatches[ia].abs_diff > angle_esd_tolerance) { 
	 bool at_1_is_hydrogen = false;
	 bool at_2_is_hydrogen = false;
	 bool at_3_is_hydrogen = false;
	 std::map<std::string, bool>::const_iterator it_1, it_2, it_3;
	 it_1 = hydrogen_status.find(angle_esd_mismatches[ia].atom_id_1);
	 it_2 = hydrogen_status.find(angle_esd_mismatches[ia].atom_id_2);
	 it_3 = hydrogen_status.find(angle_esd_mismatches[ia].atom_id_3);
	 if (it_1 != hydrogen_status.end()) at_1_is_hydrogen = it_1->second;
	 if (it_2 != hydrogen_status.end()) at_2_is_hydrogen = it_2->second;
	 if (it_3 != hydrogen_status.end()) at_3_is_hydrogen = it_3->second;
	    
	 if (compare_hydrogens || (!at_1_is_hydrogen && !at_2_is_hydrogen && !at_3_is_hydrogen)) { 
	 
	    std::cout << "Angle-Restraint:: " << comp_id_s << " mismatch angle esd between "
		      << angle_esd_mismatches[ia].atom_id_1 << " "
		      << angle_esd_mismatches[ia].atom_id_2 << "  "
		      << angle_esd_mismatches[ia].atom_id_3 << "  "
		      << angle_esd_mismatches[ia].angle_1 << " vs "
		      << angle_esd_mismatches[ia].angle_2 << "  delta: ";
	    std::cout.width(6);
	    std::cout << angle_esd_mismatches[ia].diff << "\n";
	    angles_match = false;
	 }
      }
   }



   if (angles_match)
      if (! quiet)
	 std::cout << "Angle-Restraint:: " << comp_id_s << "  all angles match within tolerance "
		   << std::endl;

   // ------------------------- torsion restraints -------------------

   // ------------------------- chiral restraints -------------------
   // 
   if (chiral_restraint.size() != r.chiral_restraint.size())
      std::cout << "Chiral-Restraint:: " << comp_id_s << " mismatch number of restraints "
		<< chiral_restraint.size() << " " << r.chiral_restraint.size()
		<< std::endl;
   bool chirals_match = true;
   for (unsigned int i=0; i<chiral_restraint.size(); i++) { 
      for (unsigned int j=0; j<r.chiral_restraint.size(); j++) {
	 if (chiral_restraint[i].matches_names(r.chiral_restraint[j])) {
	    if (chiral_restraint[i].is_a_both_restraint() !=
		r.chiral_restraint[j].is_a_both_restraint()) {
	       std::cout << "Chiral-Restraint:: " << comp_id_s << " id "
			 << chiral_restraint[i].Chiral_Id()
			 << " mismatch 'both' type " 
			 << chiral_restraint[i].volume_sign << " vs "
			 << r.chiral_restraint[j].volume_sign
			 << std::endl;
	       chirals_match = false;
	    } else {
	       if (chiral_restraint[i].volume_sign !=
		   r.chiral_restraint[j].volume_sign) {
		  std::cout << "Chiral-Restraint:: " << comp_id_s << " id "
			    << chiral_restraint[i].Chiral_Id()
			    << " mismatch volume sign " 
			    << chiral_restraint[i].volume_sign << " vs "
			    << r.chiral_restraint[j].volume_sign
			    << std::endl;
		  chirals_match = false;
	       }
	    }
	 }
      }
   }
   if (chirals_match)
      if (! quiet)
	 std::cout << "Chiral-Restraint:: " << comp_id_s << " all chiral restraints match"
		   << std::endl;
   

   // -------------------------  plane restraints -------------------
   if (plane_restraint.size() != r.plane_restraint.size())
      std::cout << "Plane-Restraint:: " << comp_id_s << " mismatch number of restraints "
		<< plane_restraint.size() << " " << r.plane_restraint.size()
		<< std::endl;
   bool planes_match = false;
   unsigned int n_planes_matched = 0;
   for (unsigned int i=0; i<plane_restraint.size(); i++) {
      bool matched_this_plane = false;
      for (unsigned int j=0; j<r.plane_restraint.size(); j++) {
	 if (plane_restraint[i].matches_names(r.plane_restraint[j])) {
	    matched_this_plane = true;
	    break;
	 }
      }
      if (matched_this_plane) {
	 n_planes_matched++;
      } else { 
	 std::cout << "Plane-Restraint:: " << comp_id_s << " no match for plane restraints "
		   << plane_restraint[i].plane_id << std::endl;
      }
   }
   if (n_planes_matched == plane_restraint.size())
      if (n_planes_matched == r.plane_restraint.size())
	 planes_match = true;
   //
   if (planes_match) {
      if (! quiet)
	 std::cout << "Plane-Restraint:: " << comp_id_s << "  all plane restraints match" << std::endl;
   } else {
      std::cout << "Plane-Restraint:: " << comp_id_s << " plane restraints do not match" << std::endl;
   }


   // Statistics summary:
   //
   double sum_sqrd = 0.0;
   double sum = 0;
   double max_pos_diff = 0;
   double max_neg_diff = 0;
   double sum_fabs = 0;
   if (bond_length_mismatches.size() > 0) { 
      for (unsigned int ib=0; ib<bond_length_mismatches.size(); ib++) { 
	 sum      += bond_length_mismatches[ib].diff;
	 sum_sqrd += bond_length_mismatches[ib].diff * bond_length_mismatches[ib].diff;
	 sum_fabs += bond_length_mismatches[ib].abs_diff;
	 if (bond_length_mismatches[ib].diff > 0) { 
	    if (bond_length_mismatches[ib].diff > max_pos_diff)
	       max_pos_diff = bond_length_mismatches[ib].diff;
	 } else {
	    if (fabs(bond_length_mismatches[ib].diff) > fabs(max_neg_diff))
	       max_neg_diff = bond_length_mismatches[ib].diff;
	 }
      }
      double inv_N = 1.0/double(bond_length_mismatches.size());
      double average = sum * inv_N;
      double mean_fabs = sum_fabs * inv_N;
      double var     = sum_sqrd * inv_N - average * average;
      if (var < 0) var = 0;
      std::cout << "\nBond Length Difference Statistics Summary::\n";
      std::cout << "            Mean:           " << std::setw(6) << average      << " A\n";
      std::cout << "            Std. Dev:       " << std::setw(6) << sqrt(var)    << " A\n";
      std::cout << "            Mean Abs Values " << std::setw(6) << mean_fabs    << " A\n";
      std::cout << std::right;
      std::cout << "            Max Pos Diff:   ";
      std::cout.width(6); // affects the next thing written out
      std::cout << max_pos_diff << " A\n";
      std::cout << "            Max Neg Diff:   ";
      std::cout.width(6);
      std::cout << max_neg_diff << " A\n";
   }
   sum_sqrd = 0.0;
   sum = 0;
   max_pos_diff = 0;
   max_neg_diff = 0;
   if (bond_esd_mismatches.size() > 0) { 
      for (unsigned int ib=0; ib<bond_esd_mismatches.size(); ib++) { 
	 sum      += bond_esd_mismatches[ib].diff;
	 sum_sqrd += bond_esd_mismatches[ib].diff * bond_esd_mismatches[ib].diff;
	 if (bond_esd_mismatches[ib].diff > 0) { 
	    if (bond_esd_mismatches[ib].diff > max_pos_diff)
	       max_pos_diff = bond_esd_mismatches[ib].diff;
	 } else {
	    if (fabs(bond_esd_mismatches[ib].diff) > fabs(max_neg_diff))
	       max_neg_diff = bond_esd_mismatches[ib].diff;
	 }
      }
      double inv_N = 1.0/double(bond_esd_mismatches.size());
      double average = sum * inv_N;
      double var     = sum_sqrd * inv_N - average * average;
      if (var < 0) var = 0;
      std::cout << "\nBond Length ESD Difference Statistics Summary::\n";
      std::cout << "            Mean:         " << std::setw(6) << average      << " A\n";
      std::cout << "            Std. Dev:     " << std::setw(6) << sqrt(var)    << " A\n";
      std::cout << std::right;
      std::cout << "            Max Pos Diff: ";
      std::cout.width(6); // affects the next thing written out
      std::cout << max_pos_diff << " A\n";
      std::cout << "            Max Neg Diff: ";
      std::cout.width(6);
      std::cout << max_neg_diff << " A\n";
   }
   sum_sqrd = 0.0;
   sum = 0;
   max_pos_diff = 0;
   max_neg_diff = 0;
   if (angle_mismatches.size() > 0) { 
      for (unsigned int ia=0; ia<angle_mismatches.size(); ia++) { 
	 sum      += angle_mismatches[ia].diff;
	 sum_sqrd += angle_mismatches[ia].diff * angle_mismatches[ia].diff;
	 if (angle_mismatches[ia].diff > 0) { 
	    if (angle_mismatches[ia].diff > max_pos_diff)
	       max_pos_diff = angle_mismatches[ia].diff;
	 } else {
	    if (fabs(angle_mismatches[ia].diff) > fabs(max_neg_diff))
	       max_neg_diff = angle_mismatches[ia].diff;
	 }
      }
      double inv_N = 1.0/double(angle_mismatches.size());
      double average = sum * inv_N;
      double var     = sum_sqrd * inv_N - average * average;
      if (var < 0) var = 0;
      std::cout << "\nAngle Difference Statistics Summary::\n";
      std::cout << "            Mean:         " << std::setw(6) << average      << " degrees\n";
      std::cout << "            Std. Dev:     " << std::setw(6) << sqrt(var)    << " degrees\n";
      std::cout << std::right;
      std::cout << "            Max Pos Diff: ";
      std::cout.width(6); // affects the next thing written out
      std::cout << max_pos_diff << " degrees\n";
      std::cout << "            Max Neg Diff: ";
      std::cout.width(6);
      std::cout << max_neg_diff << " degrees\n";
   }
   sum_sqrd = 0.0;
   sum = 0;
   max_pos_diff = 0;
   max_neg_diff = 0;
   if (angle_esd_mismatches.size() > 0) { 
      for (unsigned int ia=0; ia<angle_esd_mismatches.size(); ia++) { 
	 sum      += angle_esd_mismatches[ia].diff;
	 sum_sqrd += angle_esd_mismatches[ia].diff * angle_esd_mismatches[ia].diff;
	 if (angle_esd_mismatches[ia].diff > 0) { 
	    if (angle_esd_mismatches[ia].diff > max_pos_diff)
	       max_pos_diff = angle_esd_mismatches[ia].diff;
	 } else {
	    if (fabs(angle_esd_mismatches[ia].diff) > fabs(max_neg_diff))
	       max_neg_diff = angle_esd_mismatches[ia].diff;
	 }
      }
      double inv_N = 1.0/double(angle_esd_mismatches.size());
      double average = sum * inv_N;
      double var     = sum_sqrd * inv_N - average * average;
      if (var < 0) var = 0;
      std::cout << "\nAngle ESD Difference Statistics Summary::\n";
      std::cout << "            Mean:         " << std::setw(6) << average      << " degrees\n";
      std::cout << "            Std. Dev:     " << std::setw(6) << sqrt(var)    << " degrees\n";
      std::cout << std::right;
      std::cout << "            Max Pos Diff: ";
      std::cout.width(6); // affects the next thing written out
      std::cout << max_pos_diff << " degrees\n";
      std::cout << "            Max Neg Diff: ";
      std::cout.width(6);
      std::cout << max_neg_diff << " degrees\n";
   }
   

   bool status = false;
   if (bonds_match == true)
      if (angles_match == true)
	 if (planes_match == true)
	    if (chirals_match == true)
	       status = true;
   
   return status;
}

// Return the number of matched atoms in first.
// 
// return a dictionary that is a copy of this dictionary, but
// trying to match the names of the atoms of ref.  Do graph
// matching to find the set of atom names that match/need to be
// changed.
// 
// std::pair<unsigned int, coot::dictionary_residue_restraints_t>

coot::dictionary_match_info_t
coot::dictionary_residue_restraints_t::match_to_reference(const coot::dictionary_residue_restraints_t &ref,
							  mmdb::Residue *residue_p,
							  const std::string &new_comp_id_in,
							  const std::string &new_compound_name) const {

   bool use_hydrogens = false; // pass this? Turning this on makes this function catatonic.
   
   dictionary_residue_restraints_t dict = *this;
   bool debug = false;
   typedef std::pair<std::string, std::string> SP;
   std::vector<SP> change_name;
   std::vector<std::string> same_names;


   std::string new_comp_id = new_comp_id_in;
   if (new_comp_id == "auto")
      new_comp_id = suggest_new_comp_id(residue_info.comp_id);

   int n_atoms = atom_info.size();
   if (use_hydrogens == false)
      n_atoms = number_of_non_hydrogen_atoms();

   mmdb::math::Graph *g_1 = make_graph(use_hydrogens);
   mmdb::math::Graph *g_2 = ref.make_graph(use_hydrogens);

   if (debug) {
      std::cout << "this-name ::::::::::::::::::: " << residue_info.comp_id << std::endl;
      std::cout << " ref-name ::::::::::::::::::: " << ref.residue_info.comp_id << std::endl;
      g_1->Print();
      g_2->Print();
   }
   
   g_1->SetName ("working-residue");
   g_1->MakeVertexIDs();
   
   g_2->SetName ("reference-residue");
   g_2->MakeVertexIDs();

   bool use_bond_order = false;
   use_bond_order = true;

   g_1->MakeSymmetryRelief(false);
   g_2->MakeSymmetryRelief(false);
   
   mmdb::math::GraphMatch match;
   int minMatch = ref.number_of_non_hydrogen_atoms() - 2;
   int n_top = int(0.75 * float(ref.number_of_non_hydrogen_atoms()));
   if (minMatch <  3) minMatch =  3;
   if (minMatch > 14) minMatch = 14;

   // 20160908-PE I removed this line.  I am not sure why it is needed and for matching
   // 3GP onto (ref) GTP if this line is in place, GetNofMatches() returns 0.
   // The reverse logic seems more sensible to me now.
   // 
   // if (minMatch < n_top) minMatch = n_top;
   // 
   if (minMatch > n_top) minMatch = n_top;

   bool vertext_type = true;
   std::string s;
   if (!use_hydrogens)
      s = " non-hydrogen";

   if (false)
      std::cout << "INFO:: Matching Graphs with minMatch " << minMatch << " with "
		<< n_atoms << s << " atoms in this and " << ref.number_of_non_hydrogen_atoms()
		<< " in ref" << std::endl;
   
   int build_result_1 = g_1->Build(use_bond_order);

   if (build_result_1 != 0) {
      std::cout << "Bad graph build result_1" << std::endl;
   } else {
      int build_result_2 = g_2->Build(use_bond_order);
      if (build_result_2 != 0) {
	 std::cout << "Bad graph build result_2" << std::endl;
      } else {
	 if (debug)
	    std::cout << "debug:: minMatch is " << minMatch << std::endl;

	 // match.SetMaxNofMatches(100, true); // only need find first 100 matches
	 match.SetTimeLimit(2); // seconds
	 match.MatchGraphs(g_1, g_2, minMatch, vertext_type);
	 int n_match = match.GetNofMatches();
	 if (debug) 
	    std::cout << "found " << n_match << " matches" << std::endl;
	 if (n_match > 0) {

	    // first find the best match:
	    //
	    int imatch_best = 0;
	    int n_best_match = 0;
	    for (int imatch=0; imatch<n_match; imatch++) {
	       int nv;
	       mmdb::realtype p1, p2;
	       mmdb::ivector FV1, FV2;
	       match.GetMatch(imatch, FV1, FV2, nv, p1, p2); // n p1 p2 set
	       if (nv > n_best_match) {
		  n_best_match = nv;
		  imatch_best = imatch;
	       }
	    }

	    if (n_best_match > 0) {

	       mmdb::realtype p1, p2;
	       mmdb::ivector FV1, FV2;
	       int nv;
	       match.GetMatch(imatch_best, FV1, FV2, nv, p1, p2); // n p1 p2 set

	       for (int ipair=1; ipair<=nv; ipair++) {
		  mmdb::math::Vertex *V1 = g_1->GetVertex ( FV1[ipair] );
		  mmdb::math::Vertex *V2 = g_2->GetVertex ( FV2[ipair] );
		  if ((!V1) || (!V2))  {
		     std::cout << "Can't get vertices for match " << ipair << std::endl;
		  } else {
		     const int &type_1 = V1->GetType();
		     const int &type_2 = V2->GetType();
		     if (type_1 == type_2) {
			std::string v1_name(V1->GetName());
			std::string v2_name(V2->GetName());
			if (debug)
			   std::cout << " atom " << V1->GetName()
				     << " in graph_1 (this-res) matches atom "
				     << V2->GetName() << " in graph_2 (ref comp-id)" << std::endl;
			if (v1_name != v2_name) {
			   change_name.push_back(SP(v1_name, v2_name));
			} else {
			   same_names.push_back(v1_name);
			}
		     } else {
			std::cout << "ERROR:: excluding match between "
				  << V1->GetName() << " and "
				  << V2->GetName() << " because types are " << type_1
				  << " and " << type_2 << std::endl;
		     }
		  }
	       }
	    }

	    if (debug) {
	       std::cout << "----- accumulated ------ " << change_name.size()
			 << " name changes " << "(and " << same_names.size()
			 << " matches of atoms with the same name)"
			 << " for " << n_atoms << " atoms"
			 << " --------- " << std::endl;
	       for (unsigned int i=0; i<change_name.size(); i++) { 
		  std::cout << i << "  " << change_name[i].first << " -> "
			    << change_name[i].second << std::endl;
	       }
	       std::cout << "--- these matched and were the same atom name" << std::endl;
	       for (unsigned int i=0; i<same_names.size(); i++) {
		  std::cout << i << "  " << same_names[i] << std::endl;
	       }
	       std::cout << "DEBUG:: in total: " << change_name.size() + same_names.size()
			 << " name matches" << std::endl;
	    }
	    
	    // also header info.
	    dict.residue_info.comp_id           = new_comp_id;
	    dict.residue_info.three_letter_code = new_comp_id;
	    dict.residue_info.name              = new_compound_name;
	    dict.residue_info.group             = ref.residue_info.group; // probably right.

	    // do any of the target (to) names exist in dict already?  If so,
	    // we will need to invent a new name for those already-existing
	    // atoms.
	    std::vector<SP> more_swaps_from_name_clashes = extra_name_swaps_from_name_clash(change_name);
	    for (unsigned int ii=0; ii<more_swaps_from_name_clashes.size(); ii++)
	       change_name.push_back(more_swaps_from_name_clashes[ii]);

	    // do the swap
	    dict.atom_id_swap(change_name);

	    // change the residue atom names too (if non-NULL).
	    if (residue_p) {
	       bool something_changed = change_names(residue_p, change_name, new_comp_id);
	       if (something_changed) {
		  mmdb::Manager *m = residue_p->chain->GetCoordHierarchy();
		  if (m)
		     m->FinishStructEdit();
	       }
	    } 
	 }
      }
   }

   delete g_1;
   delete g_2;

   dictionary_match_info_t dmi;
   dmi.n_matches = change_name.size() + same_names.size();
   dmi.dict = dict;
   dmi.name_swaps = change_name;
   dmi.same_names = same_names;
   dmi.new_comp_id = new_comp_id;
   
   return dmi;

}

bool
coot::dictionary_residue_restraints_t::change_names(mmdb::Residue *residue_p,
						    const std::vector<std::pair<std::string, std::string> > &change_name, const std::string &new_comp_id) const {

   bool changed_something = false;

   if (residue_p) {
      mmdb::PPAtom res_selection = NULL;
      int num_residue_atoms;
      residue_p->GetAtomTable(res_selection, num_residue_atoms);
      for (int iat=0; iat<num_residue_atoms; iat++) {
	 mmdb::Atom *at = res_selection[iat];
	 std::string atom_name = at->name;
	 for (unsigned int j=0; j<change_name.size(); j++) { 
	    if (change_name[j].first == atom_name) {
	       // 4 chars?
	       at->SetAtomName(change_name[j].second.c_str());
	       changed_something = true;
	       break;
	    }
	 }
      }
   }
   if (changed_something) {
      residue_p->SetResName(new_comp_id.c_str());
   } 
   return changed_something;
}


// do any of the target (to) names exist in dict already (and that to
// name is not assigned to be replaced)?  If so, we will need to
// invent a new name for those already-existing atoms.
// 
std::vector<std::pair<std::string, std::string> >
coot::dictionary_residue_restraints_t::extra_name_swaps_from_name_clash(const std::vector<std::pair<std::string, std::string> > &change_name) const {

   typedef std::pair<std::string, std::string> SP;
   std::vector<SP> r;

   std::vector<std::string> invented_names;
   for (unsigned int i=0; i<change_name.size(); i++) {
      for (unsigned int j=0; j<atom_info.size(); j++) {
	 const std::string &to_name = change_name[i].second;
	 // if it is an atom name that already exists in this residue...
	 if (to_name == atom_info[j].atom_id) {
	    // and if that atom name is not assigned to be changed... 
	    bool found = false;
	    for (unsigned int k=0; k<change_name.size(); k++) {
	       if (change_name[k].first == to_name) {
		  found = true;
		  break;
	       }
	    }
	    if (! found) {
	       // not assigned to be changed...
	       //
	       // so invent a new name.
	       std::string ele = "C";
	       for (unsigned int jj=0; jj<atom_info.size(); jj++) { 
		  if (atom_info[jj].atom_id == to_name) {
		     ele = atom_info[jj].type_symbol;
		     break;
		  } 
	       }
	       
	       std::string invented_name = invent_new_name(ele, invented_names);

	       if (0) 
		  std::cout << "extra_name_swaps_from_name_clash() " << i << " " << j
			    << " invented name: " << invented_name << std::endl;
	       invented_names.push_back(invented_name);
	       SP p(to_name, invented_name);
	       r.push_back(p);
	    } 
	 }
      }
   }

   return r;
}

std::string
coot::dictionary_residue_restraints_t::invent_new_name(const std::string &ele,
						       const std::vector<std::string> &other_invented_names) const {

   std::string new_name("XXX");
   std::string a("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
   bool found = false;

   std::vector<std::string> monomer_atom_names(atom_info.size());
   for (unsigned int iat=0; iat<atom_info.size(); iat++)
      monomer_atom_names[iat] = atom_info[iat].atom_id;

   std::vector<std::string> existing_atom_names = monomer_atom_names;
   for (unsigned int i=0; i<other_invented_names.size(); i++)
      existing_atom_names.push_back(other_invented_names[i]);
   
   for (unsigned int i=0; i<a.size(); i++) { 
      for (unsigned int j=0; j<a.size(); j++) {
	 std::string test_atom_name = "";
	 if (ele.length() == 1) { 
	    test_atom_name = " ";
	    test_atom_name += ele;
	 } else {
	    test_atom_name = ele;
	 }
	 test_atom_name += a[i];
	 test_atom_name += a[j];
	 if (std::find(existing_atom_names.begin(), existing_atom_names.end(), test_atom_name)
	     == existing_atom_names.end()) {
	    found = true;
	    new_name = test_atom_name;
	 }
	 if (found)
	    break;
      }
      if (found)
	 break;
   }
   return new_name;
} 



mmdb::math::Graph *
coot::dictionary_residue_restraints_t::make_graph(bool use_hydrogens) const { 

   std::map<std::string, unsigned int> name_map;

   
   // mol atom indexing -> graph vertex indexing (we need this because
   // (perhaps) not all atoms in the atom_info are atoms in the graph).
   //
   // This contains the atom indices of of atom_info (0-indexed).  An
   // mmdb graph is 1-indexed.
   // 
   int vertex_indexing[atom_info.size()];
   
   mmdb::math::Graph *graph = new mmdb::math::Graph;
   int i_atom = 0;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) { 
      std::string ele  = atom_info[iat].type_symbol;
      if (use_hydrogens || (ele != " H" && ele != "H")) { 
	 std::string name = atom_info[iat].atom_id_4c;
	 mmdb::math::Vertex *v = new mmdb::math::Vertex(ele.c_str(), name.c_str());
	 graph->AddVertex(v);
	 name_map[name] = i_atom;
	 i_atom++;
      }
   }

   for (unsigned int ib=0; ib<bond_restraint.size(); ib++) {
      const dict_bond_restraint_t &br = bond_restraint[ib];
      int mmdb_bond_type = br.mmdb_bond_type();
      // std::cout << "br: " << br << " mmdb_bond_type: " << mmdb_bond_type << std::endl;
      if (mmdb_bond_type != -1) {
	 std::map<std::string, unsigned int>::const_iterator it_1;
	 std::map<std::string, unsigned int>::const_iterator it_2;
	 it_1 = name_map.find(br.atom_id_1_4c());
	 it_2 = name_map.find(br.atom_id_2_4c());

	 if (it_1 == name_map.end()) {
	    if (use_hydrogens || !is_hydrogen(br.atom_id_1_4c()))
	       std::cout << "Not found in name map atom 1 :" << br.atom_id_1() << ":" << std::endl;
	 } else { 
	    if (it_2 == name_map.end()) {
	       if (use_hydrogens || !is_hydrogen(br.atom_id_2_4c()))
		  std::cout << "Not found in name map atom 2 :" << br.atom_id_2() << ":" << std::endl;
	    } else {

	       if (use_hydrogens || (!is_hydrogen(br.atom_id_1_4c()) &&
				     !is_hydrogen(br.atom_id_2_4c()))) {
		  mmdb::math::Edge *e = new mmdb::math::Edge(it_1->second + 1,  // 1-indexed
							     it_2->second + 1,
							     mmdb_bond_type);
		  graph->AddEdge(e);
	       }
	    }
	 }
      }
   }
   return graph;
}


bool
coot::dictionary_residue_restraints_t::comprised_of_organic_set() const { 

   bool r = true;
   if (atom_info.size() == 0) {
      r = false;
   } else {
      std::vector<std::string> e;
      // e.push_back(" H"); e.push_back(" C"); e.push_back(" N"); e.push_back(" O");
      // e.push_back(" I"); e.push_back("BR"); e.push_back("CL"); e.push_back(" F");
      // e.push_back(" S"); e.push_back(" P");
      e.push_back("H"); e.push_back("C"); e.push_back("N"); e.push_back("O");
      e.push_back("I"); e.push_back("BR"); e.push_back("CL"); e.push_back("F");
      e.push_back("S"); e.push_back("P");
      for (unsigned int i=0; i<atom_info.size(); i++) {
	 bool found_this = false;
	 for (unsigned int j=0; j<e.size(); j++) {
	    if (atom_info[i].type_symbol == e[j]) {
	       found_this = true;
	       break;
	    }
	 }
	 if (! found_this) {
	    std::cout << "INFO::organic_set_test: " << i << " " << atom_info[i]
		      << " \"" << atom_info[i].type_symbol << "\""
		      << " is not in the organic set" << std::endl;
	    r = false;
	    break;
	 }
      }
   }
   return r;
}


// are the number of atoms of each element the same ie. they have the same chemical formula?
//
bool coot::dictionary_residue_restraints_t::composition_matches(const dictionary_residue_restraints_t &refr) const {

   // Needs testing.

   bool r = true;

   std::map<std::string, int> this_atom_count;
   std::map<std::string, int> refr_atom_count;
   for (unsigned int i=0; i<atom_info.size(); i++)
      this_atom_count[atom_info[i].type_symbol]++;
   for (unsigned int i=0; i<refr.atom_info.size(); i++)
      refr_atom_count[refr.atom_info[i].type_symbol]++;

   std::map<std::string, int>::const_iterator it;
   std::map<std::string, int>::const_iterator it_refr;
   for (it=this_atom_count.begin(); it!=this_atom_count.end(); it++) {
      it_refr = refr_atom_count.find(it->first);
      if (it_refr == refr_atom_count.end()) {
	 r = false;
	 break;
      } else {
	 if (it->second != it_refr->second) {
	    r = false;
	    break;
	 }
      }
   }
   return r;
}


// for hydrogens
bool
coot::dictionary_residue_restraints_t::is_connected_to_donor(const std::string &H_at_name_4c,
							     const coot::energy_lib_t &energy_lib) const {

   bool result = false;
   for (unsigned int i=0; i<bond_restraint.size(); i++) {
      if (bond_restraint[i].atom_id_1_4c() == H_at_name_4c) {
	 // what is the energy type of bond[i].atom_id_2_4c()?
	 std::string energy_type = type_energy(bond_restraint[i].atom_id_2_4c());
	 std::map<std::string, energy_lib_atom>::const_iterator it = energy_lib.atom_map.find(energy_type);
	 if (it != energy_lib.atom_map.end()) {
	    if (it->second.hb_type == HB_DONOR || it->second.hb_type == HB_BOTH) {
	       result = true;
	       break;
	    }
	 }
      }
      if (bond_restraint[i].atom_id_2_4c() == H_at_name_4c) {
	 // what is the energy type of bond[i].atom_id_1_4c()?
	 std::string energy_type = type_energy(bond_restraint[i].atom_id_1_4c());
	 std::map<std::string, energy_lib_atom>::const_iterator it = energy_lib.atom_map.find(energy_type);
	 if (it != energy_lib.atom_map.end()) {
	    if (it->second.hb_type == HB_DONOR || it->second.hb_type == HB_BOTH) {
	       result = true;
	       break;
	    }
	 }
      }
   }
   return result;
}



// 1mzt

std::vector<std::string>
coot::comp_ids_in_dictionary_cif(const std::string &cif_dictionary_filename) {

   std::vector<std::string> v;
   coot::protein_geometry geom;
   geom.set_verbose(false);
   int read_number = 0;  // doesn't matter
   int imol_enc = protein_geometry::IMOL_ENC_ANY;
   geom.init_refmac_mon_lib(cif_dictionary_filename, read_number, imol_enc);
   v = geom.monomer_restraints_comp_ids();
   return v;
}



bool
coot::dict_plane_restraint_t::matches_names(const coot::dict_plane_restraint_t &r) const {

   bool status = true;
   unsigned int n_found = 0;
   if (atom_ids.size() != r.atom_ids.size())
      return false;
   if (atom_ids.size() > 0)
      status = false; // initial setting.
   
   for (unsigned int i=0; i<atom_ids.size(); i++) {
      const std::string &ref_atom = atom_ids[i].first;
      for (unsigned int j=0; j<r.atom_ids.size(); j++) { 
	 if (r.atom_ids[j].first == ref_atom) {
	    n_found++;
	    break;
	 }
      }
   }
   if (n_found == atom_ids.size())
      status = true;
   return status;
}


std::string
coot::atom_id_mmdb_expand(const std::string &atomname) { 
   std::string r;
   int ilen = atomname.length();
      
   if (ilen == 4) return atomname;
      
   if (ilen == 1) {
      r = " ";
      r += atomname;
      r += "  ";
   } else {
      if (ilen == 2) {

	 // 20180512-PE we have to be more clever here for metals.
	 // But what about CA! - argh! we shouldn't be using this function.
	 // We need to know the residue name to pad correctly.
	 //
	 bool done = false;
	 if (atomname == "MG" || atomname == "NA" || atomname == "LI" || atomname == "LI" || atomname == "AL" || atomname == "SI" ||
	     atomname == "CL" || atomname == "SC" || atomname == "TI" || atomname == "CR" || atomname == "MN" || atomname == "FE" ||
	     atomname == "CO" || atomname == "NI" || atomname == "CU" || atomname == "ZN" || atomname == "GA" || atomname == "AS" ||
	     atomname == "SE" || atomname == "BR" || atomname == "RB" || atomname == "SR" || atomname == "RE" || atomname == "OS" ||
	     atomname == "IR" || atomname == "PT" || atomname == "AU" || atomname == "HG" || atomname == "PB" || atomname == "BI") {
	    r += atomname;
	    r += "  ";
	 } else {
	    r = " ";
	    r += atomname;
	    r += " ";
	 }
      } else {
	 if (ilen == 3) {
	    r = " ";
	    r += atomname;
	 } else {
	    r = atomname;
	 }
      }
   }
   return r;
}

std::string
coot::atom_id_mmdb_expand(const std::string &atomname, const std::string &element) {

   std::string r = coot::atom_id_mmdb_expand(atomname);

   if (element.length() == 2 && element[0] != ' ') {
      if (atomname.length() == 1) { // unlikely
	 r = " ";
	 r += atomname;
	 r += "  ";
      } else {
	 if (atomname.length() == 2) {
	    r = atomname;
	    r += "  ";
	 } else {
	    if (atomname.length() == 3) {
	       r = atomname;
	       r += " ";
	    } else {
	       r = atomname;
	    }
	 }
      }
   }
   if (0)  // debug
      std::cout << "Given :" << atomname << ": and element :" <<
	 element << ": returning :" << r << ":" << std::endl;
   return r;
}


bool
coot::dict_torsion_restraint_t::is_pyranose_ring_torsion(const std::string &comp_id) const {

   // Needs fixup for PDBv3
   bool status = false;
   std::string ring_atoms[6] = { " C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 " };
   if (comp_id == "XYP")
      for (unsigned int i=0; i<6; i++)
	 ring_atoms[i][3] = 'B'; // danger on PDBv3 fixup.

   int n_matches = 0;
   for (unsigned int i=0; i<6; i++) { 
      if (atom_id_2_4c() == ring_atoms[i])
	 n_matches++;
      if (atom_id_3_4c() == ring_atoms[i])
	 n_matches++;
   }
   if (n_matches == 2)
      status = true;
   return status;
}

bool
coot::dict_link_torsion_restraint_t::is_pyranose_ring_torsion() const {

   // Needs fixup for PDBv3
   bool status = false;
   std::string ring_atoms[6] = { " C1 ", " C2 ", " C3 ", " C4 ", " C5 ", " O5 " };

   int n_matches = 0;
   for (unsigned int i=0; i<6; i++) { 
      if (atom_id_2_4c() == ring_atoms[i])
	 n_matches++;
      if (atom_id_3_4c() == ring_atoms[i])
	 n_matches++;
   }
   if (n_matches == 2)
      status = true;
   return status;
} 

bool
coot::dict_torsion_restraint_t::is_ring_torsion(const std::vector<std::vector<std::string> > &ring_atoms_sets) const {

   bool match = false; 
   std::vector<std::string> torsion_restraint_atom_names(2);
   torsion_restraint_atom_names[0] = atom_id_2_4c();
   torsion_restraint_atom_names[1] = atom_id_3_4c();
   
   for (unsigned int iring=0; iring<ring_atoms_sets.size(); iring++) { 
      const std::vector<std::string> &ring_atom_names = ring_atoms_sets[iring];

      int n_match = 0;
      for (unsigned int iname_1=0; iname_1<ring_atom_names.size(); iname_1++) {
	 for (unsigned int iname_2=0; iname_2<torsion_restraint_atom_names.size(); iname_2++) { 
	    if (ring_atom_names[iname_1] == torsion_restraint_atom_names[iname_2])
	       n_match++;
	 }
      }
      if (n_match == 2) {
	 match = true;
	 break;
      }
   }
   return match;
} 


bool
coot::dict_torsion_restraint_t::is_const() const {

   bool const_flag = 0;
   if (id_.length() > 5) {
      std::string bit = id_.substr(0,5);
      if (bit == "CONST")
	 const_flag = 1;
      if (bit == "const")
	 const_flag = 1;
   }
   return const_flag;
}



bool
coot::dict_atom::is_hydrogen() const {

   bool r = false;
   if (type_symbol == "H" ||
       type_symbol == " H" ||
       type_symbol == "D")
      r = true;
   return r;
}



void
coot::dict_atom::add_pos(int pos_type,
			 const std::pair<bool, clipper::Coord_orth> &model_pos) {

   if (pos_type == coot::dict_atom::IDEAL_MODEL_POS)
      pdbx_model_Cartn_ideal = model_pos;
   if (pos_type == coot::dict_atom::REAL_MODEL_POS) {
      model_Cartn = model_pos;
   }

}

// quote atom name as needed - i.e. CA -> CA, CA' -> "CA'"
std::string 
coot::dictionary_residue_restraints_t::quoted_atom_name(const std::string &an) const {

   std::string n = an;
   bool has_quotes = false;

   for (unsigned int i=0; i<an.size(); i++) {
      if (an[i] == '\'') {
	 has_quotes = true;
	 break;
      } 
   }
   if (has_quotes)
      n = "\"" + an + "\"";

   return n;
}


// make a connect file specifying the bonds to Hydrogens
bool
coot::protein_geometry::hydrogens_connect_file(const std::string &resname,
					       const std::string &filename) const {

   bool r = 0;
   std::pair<short int, dictionary_residue_restraints_t> p =
      get_monomer_restraints(resname, IMOL_ENC_ANY);

   if (p.first) {
      std::vector<dict_bond_restraint_t> bv = p.second.bond_restraint;
      if (bv.size() > 0) {
	 // try to open the file then:
	 std::ofstream connect_stream(filename.c_str());
	 if (connect_stream) {
	    int n_atoms = p.second.atom_info.size();
	    connect_stream << "# Generated by Coot" << std::endl;
	    connect_stream << "RESIDUE   " << resname << "   " << n_atoms << std::endl;
	    std::vector<std::pair<std::string, std::vector<std::string> > > assoc; 
	    for (unsigned int i=0; i<bv.size(); i++) {
	       std::string atom1 = bv[i].atom_id_1();
	       std::string atom2 = bv[i].atom_id_2();
	       // find atom1
	       bool found = 0;
	       int index_1 = -1;
	       int index_2 = -1;
	       for (unsigned int j=0; j<assoc.size(); j++) {
		  if (atom1 == assoc[j].first) {
		     found = 1;
		     index_1 = j;
		     break;
		  } 
	       }
	       if (found == 1) {
		  assoc[index_1].second.push_back(atom2);
	       } else {
		  // we need to add a new atom:
		  std::vector<std::string> vt;
		  vt.push_back(atom2);
		  std::pair<std::string, std::vector<std::string> > p(atom1, vt);
		  assoc.push_back(p);
	       }
	       // find atom2
	       found = 0;
	       for (unsigned int j=0; j<assoc.size(); j++) {
		  if (atom2 == assoc[j].first) {
		     found = 1;
		     index_2 = j;
		     break;
		  }
	       }
	       if (found == 1) {
		  assoc[index_2].second.push_back(atom1);
	       } else {
		  // we need to add a new atom:
		  std::vector<std::string> vt;
		  vt.push_back(atom1);
		  std::pair<std::string, std::vector<std::string> > p(atom2, vt);
		  assoc.push_back(p);
	       } 
	    }

	    r = 1;
	    // for each atom in assoc
	    for (unsigned int i=0; i<assoc.size(); i++) {
	       
	       connect_stream << "CONECT     " << assoc[i].first << "    "
			      << assoc[i].second.size();
	       for (unsigned int ii=0; ii<assoc[i].second.size(); ii++) {
		  connect_stream << assoc[i].second[ii] << " ";
	       }
	       connect_stream << std::endl;
	    }
	 }
      }
   }
   return r;
} 
      
