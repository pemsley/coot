
#include <map>
#include <algorithm>
#include <iomanip> // setw()
#include "utils/coot-utils.hh"
#include "protein-geometry.hh"
#include "dict-mismatches.hh"

// quiet means don't tell me about matches
bool
coot::dictionary_residue_restraints_t::compare(const dictionary_residue_restraints_t &r,
					       double bond_length_tolerance,
					       double bond_esd_tolerance,
					       double angle_tolerance,
					       double angle_esd_tolerance,
					       bool compare_hydrogens,
					       bool quiet) const {

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
   if (r.residue_info.three_letter_code != residue_info.three_letter_code) {
      std::cout << "Residue-Info:: mismatch three_letter_code "
		<< residue_info.three_letter_code << " vs "
		<< r.residue_info.three_letter_code
		<< std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.name != residue_info.name) {
      std::cout << "Residue-Info:: mismatch name \""
		<< residue_info.name << "\" vs \"" << r.residue_info.name << "\""
		<< std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.group != residue_info.group) {
      std::cout << "Residue-Info:: mismatch group "
		<< residue_info.group << " vs " << r.residue_info.group << std::endl;
      residue_info_matches = false;
   }
   if (compare_hydrogens) { 
      if (r.residue_info.number_atoms_all != residue_info.number_atoms_all) {
	 std::cout << "Residue-Info:: mismatch number_atoms_all "
		   << residue_info.number_atoms_all << " " << r.residue_info.number_atoms_all
		   << std::endl;
	 residue_info_matches = false;
      }
   }
   if (r.residue_info.number_atoms_nh != residue_info.number_atoms_nh) {
      std::cout << "Residue-Info:: mismatch number_atoms_nh "
		<< residue_info.number_atoms_nh << " " << r.residue_info.number_atoms_nh
		<< std::endl;
      residue_info_matches = false;
   }
   if (residue_info_matches)
      if (! quiet)
	 std::cout << "Residue-Info::     all residue attributes match " << std::endl;

   // atom info
   if (atom_info.size() != r.atom_info.size()) {
      std::cout << "Atom-Info:: mismatch number of atoms " << atom_info.size() << " vs "
	 << r.atom_info.size() << std::endl;
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
		  std::cout << "Atom-Info:: " << atom_id_refr << " type_symbol mismatch "
			    << atom_info[iat].type_symbol << " vs "
			    << r.atom_info[jat].type_symbol << " "
			    << std::endl;
	       }
	       if (atom_info[iat].type_energy != r.atom_info[jat].type_energy) {
		  std::cout << "Atom-Info:: " << atom_info[iat].atom_id_4c
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
      std::cout << "Bond-Restraint:: mismatch number of restraints "
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
	    std::string up_1 = coot::util::upcase(type_1);
	    std::string up_2 = coot::util::upcase(type_2);
	    if (up_1 != up_2) {
	       std::cout << "Bond-Restraint:: mismatch bond order  between "
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
      if (bond_length_mismatches[ib].abs_diff > bond_length_tolerance) {
	 bool at_1_is_hydrogen = false;
	 bool at_2_is_hydrogen = false;
	 std::map<std::string, bool>::const_iterator it_1, it_2;
	 it_1 = hydrogen_status.find(bond_length_mismatches[ib].atom_id_1);
	 it_2 = hydrogen_status.find(bond_length_mismatches[ib].atom_id_2);
	 if (it_1 != hydrogen_status.end()) at_1_is_hydrogen = it_1->second;
	 if (it_2 != hydrogen_status.end()) at_2_is_hydrogen = it_2->second;
	 
	 if (compare_hydrogens || (!at_1_is_hydrogen && !at_2_is_hydrogen)) { 
	    std::cout << "Bond-Restraint:: mismatch bond length between     "
		      << bond_length_mismatches[ib].atom_id_1 << " "
		      << bond_length_mismatches[ib].atom_id_2 << "  "
		      << bond_length_mismatches[ib].dist_1 << " vs "
		      << bond_length_mismatches[ib].dist_2 << "  delta: ";
	    std::cout.width(6);
	    std::cout << bond_length_mismatches[ib].diff << "\n";
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
	    std::cout << "Bond-Restraint:: mismatch bond length esd between "
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
	 std::cout << "Bond-Restraint::   all bonds  match within tolerance "
		   << std::endl;

   // -------------------------  angle restraints -------------------
   // 
   if (angle_restraint.size() != r.angle_restraint.size())
      std::cout << "Angle-Restraint:: mismatch number of restraints "
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
      if (angle_mismatches[ia].abs_diff > angle_tolerance) {

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
	 
	    std::cout << "Angle-Restraint:: mismatch angle between     "
		      << angle_mismatches[ia].atom_id_1 << " "
		      << angle_mismatches[ia].atom_id_2 << "  "
		      << angle_mismatches[ia].atom_id_3 << "  "
		      << angle_mismatches[ia].angle_1 << " vs "
		      << angle_mismatches[ia].angle_2 << "  delta: ";
	    std::cout.width(6);
	    std::cout << angle_mismatches[ia].diff << "\n";
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
	 
	    std::cout << "Angle-Restraint:: mismatch angle esd between "
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
	 std::cout << "Angle-Restraint::  all angles match within tolerance "
		   << std::endl;

   // ------------------------- torsion restraints -------------------

   // ------------------------- chiral restraints -------------------
   // 
   if (chiral_restraint.size() != r.chiral_restraint.size())
      std::cout << "Chiral-Restraint:: mismatch number of restraints "
		<< chiral_restraint.size() << " " << r.chiral_restraint.size()
		<< std::endl;
   bool chirals_match = true;
   for (unsigned int i=0; i<chiral_restraint.size(); i++) { 
      for (unsigned int j=0; j<r.chiral_restraint.size(); j++) {
	 if (chiral_restraint[i].matches_names(r.chiral_restraint[j])) {
	    if (chiral_restraint[i].is_a_both_restraint() !=
		r.chiral_restraint[j].is_a_both_restraint()) {
	       std::cout << "Chiral-Restraint:: id " << chiral_restraint[i].Chiral_Id()
			 << " mismatch 'both' type " 
			 << chiral_restraint[i].volume_sign << " vs "
			 << r.chiral_restraint[j].volume_sign
			 << std::endl;
	       chirals_match = false;
	    } else {
	       if (chiral_restraint[i].volume_sign !=
		   r.chiral_restraint[j].volume_sign) {
		  std::cout << "Chiral-Restraint:: id " << chiral_restraint[i].Chiral_Id()
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
	 std::cout << "Chiral-Restraint:: all chiral restraints match" << std::endl;
   

   // -------------------------  plane restraints -------------------
   if (plane_restraint.size() != r.plane_restraint.size())
      std::cout << "Plane-Restraint:: mismatch number of restraints "
		<< plane_restraint.size() << " " << r.plane_restraint.size()
		<< std::endl;
   bool planes_match = false;
   int n_planes_matched = 0;
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
	 std::cout << "Plane-Restraint:: no match for plane restraints "
		   << plane_restraint[i].plane_id << std::endl;
      }
   }
   if (n_planes_matched == plane_restraint.size())
      if (n_planes_matched == r.plane_restraint.size())
	 planes_match = true;
   //
   if (planes_match) {
      if (! quiet)
	 std::cout << "Plane-Restraint::  all plane restraints match" << std::endl;
   } else {
      std::cout << "Plane-Restraint:: plane restraints do not match" << std::endl;
   }


   // Statistics summary:
   //
   double sum_sqrd = 0.0;
   double sum = 0;
   double max_pos_diff = 0;
   double max_neg_diff = 0;
   if (bond_length_mismatches.size() > 0) { 
      for (unsigned int ib=0; ib<bond_length_mismatches.size(); ib++) { 
	 sum      += bond_length_mismatches[ib].diff;
	 sum_sqrd += bond_length_mismatches[ib].diff * bond_length_mismatches[ib].diff;
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
      double var     = sum_sqrd * inv_N - average * average;
      if (var < 0) var = 0;
      std::cout << "\nBond Length Difference Statistics Summary::\n";
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


// 1mzt
