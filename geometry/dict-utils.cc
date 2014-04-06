
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
	    std::string up_1 = coot::util::upcase(type_1);
	    std::string up_2 = coot::util::upcase(type_2);
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

// return a dictionary that is a copy of this dictionary, but
// trying to match the names of the atoms of ref.  Do graph
// matching to find the set of atom names that match/need to be
// changed.
// 
coot::dictionary_residue_restraints_t
coot::dictionary_residue_restraints_t::match_to_reference(const coot::dictionary_residue_restraints_t &ref,
							  CResidue *residue_p) {
   dictionary_residue_restraints_t dict = *this;
   typedef std::pair<std::string, std::string> SP;
   std::vector<SP> change_name;

   bool use_hydrogens = true;
   int n_atoms = atom_info.size();
   if (use_hydrogens == false)
      n_atoms = number_of_non_hydrogen_atoms();
   
   CGraph *g_1 = make_graph(use_hydrogens);
   CGraph *g_2 = ref.make_graph(use_hydrogens);

   if (0) { 
      std::cout << "this-name:::::::::::::::::::" << residue_info.comp_id << std::endl;
      std::cout << " ref-name:::::::::::::::::::" << ref.residue_info.comp_id << std::endl;
      g_1->Print();
      g_2->Print();
   }
   
   g_1->SetName (residue_info.name.c_str());
   g_2->MakeVertexIDs();

   Boolean use_bond_order = false;

   g_1->MakeSymmetryRelief(False);
   
   CGraphMatch match;
   int minMatch = static_cast<int>(0.9*float(n_atoms)) - 3;
   // minMatch = n_atoms;

   Boolean vertext_type = True;
   std::string s;
   if (!use_hydrogens)
      s = " non-hydrogen";
   std::cout << "MatchGraphs() with minMatch " << minMatch << " with "
	     << n_atoms << s << " atoms" << std::endl;
   int build_result_1 = g_1->Build(False);

   if (build_result_1 != 0) {
      std::cout << "Bad graph build result_1" << std::endl;
   } else {
      int build_result_2 = g_2->Build(use_bond_order);
      if (build_result_1 != 0) {
	 std::cout << "Bad graph build result_2" << std::endl;
      } else { 
	 match.MatchGraphs(g_1, g_2, minMatch, vertext_type);
	 int n_match = match.GetNofMatches();
	 if (1) 
	    std::cout << "found " << n_match << " matches" << std::endl;
	 if (n_match > 0) {

	    int imatch_best = 0;
	    for (int imatch=0; imatch<n_match; imatch++) {
	       int nv;
	       realtype p1, p2;
	       ivector FV1, FV2;
	       match.GetMatch(imatch, FV1, FV2, nv, p1, p2); // n p1 p2 set
	       // std::cout << "   imatch " << imatch << " " << nv << std::endl;
	       int n_type_match = 0;
	       for (int ipair=1; ipair<=nv; ipair++) {
		  CVertex *V1 = g_1->GetVertex ( FV1[ipair] );
		  CVertex *V2 = g_2->GetVertex ( FV2[ipair] );
		  if ((!V1) || (!V2))  {
		     std::cout << "Can't get vertices for match " << ipair << std::endl;
		  } else {
		     const int &type_1 = V1->GetType();
		     const int &type_2 = V2->GetType();
		     if (type_1 == type_2) {
			std::string v1_name(V1->GetName());
			std::string v2_name(V2->GetName());
			if (imatch == imatch_best)
			   if (0) 
			      std::cout << " imatch_best " << imatch
					<< " atom " << V1->GetName() << " in graph_1 matches atom "
					<< V2->GetName() << " in graph_2" << std::endl;
			if (imatch == imatch_best) { 
			   if (v1_name != v2_name) {
			      change_name.push_back(SP(v1_name, v2_name));
			   }
			}
		     } else {
			std::cout << "imatch " << imatch <<  "excluding match between "
				  << V1->GetName() << " and "
				  << V2->GetName() << " because types are " << type_1
				  << " and " << type_2 << std::endl;
		     } 
		  }
	       }
	    }

	    if (0) { 
	       std::cout << "----- accumulated ------ " << change_name.size() << " name changes "
			 << " --------- " << std::endl;
	       for (unsigned int i=0; i<change_name.size(); i++) { 
		  std::cout << i << "  " << change_name[i].first << " -> " << change_name[i].second
			    << std::endl;
	       }
	    }
	    // also header info.
	    dict.residue_info.comp_id           = ref.residue_info.comp_id;
	    dict.residue_info.three_letter_code = ref.residue_info.three_letter_code;
	    dict.residue_info.name              = ref.residue_info.name;
	    dict.residue_info.group             = ref.residue_info.group;

	    // change the residue atom names too (if non-NULL).
	    change_names(residue_p, change_name);

	    // do any of the target (to) names exist in dict already?  If so,
	    // we will need to invent a new name for those already-existing
	    // atoms.
	    std::vector<SP> more_swaps_from_name_clashes = extra_name_swaps_from_name_clash(change_name);
	    for (unsigned int ii=0; ii<more_swaps_from_name_clashes.size(); ii++)
	       change_name.push_back(more_swaps_from_name_clashes[ii]);

	    // do the swap
	    dict.atom_id_swap(change_name);
	 }
      }
   }

   delete g_1;
   delete g_2;
   
   return dict;
}

void
coot::dictionary_residue_restraints_t::change_names(CResidue *residue_p,
						    const std::vector<std::pair<std::string, std::string> > &change_name) const {

   if (residue_p) {
      PPCAtom res_selection = NULL;
      int num_residue_atoms;
      residue_p->GetAtomTable(res_selection, num_residue_atoms);
      for (unsigned int iat=0; iat<num_residue_atoms; iat++) {
	 CAtom *at = res_selection[iat];
	 std::string atom_name = at->name;
	 for (unsigned int j=0; j<change_name.size(); j++) { 
	    if (change_name[j].first == atom_name) {
	       // 4 chars?
	       at->SetAtomName(change_name[j].second.c_str());
	    } 
	 }
      }
   }
}


// do any of the target (to) names exist in dict already (and that to
// name is not assigned to be replaced)?  If so, we will need to
// invent a new name for those already-existing atoms.
// 
std::vector<std::pair<std::string, std::string> >
coot::dictionary_residue_restraints_t::extra_name_swaps_from_name_clash(const std::vector<std::pair<std::string, std::string> > &change_name) const {

   typedef std::pair<std::string, std::string> SP;
   std::vector<SP> r;
   for (unsigned int i=0; i<change_name.size(); i++) {
      for (unsigned int j=0; j<atom_info.size(); j++) {
	 const std::string &to_name = change_name[i].second;
	 // if it is an atom name that already exists in this residue...
	 if (to_name == atom_info[j].atom_id) {
	    // and if that atom name is not assigned to be changed... 
	    bool found = false;
	    for (unsigned int k=0; k<change_name.size(); k++) {
	       if (change_name[k].first == to_name) {
		  found == true;
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
	       
	       std::string invented_name = invent_new_name(ele);
	       SP p(to_name, invented_name);
	       r.push_back(p);
	    } 
	 }
      }
   }

   return r;
}

std::string
coot::dictionary_residue_restraints_t::invent_new_name(const std::string &ele) const {

   std::string new_name("XXX");
   std::string a("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
   bool found = false;

   std::vector<std::string> monomer_atom_names(atom_info.size());
   for (unsigned int iat=0; iat<atom_info.size(); iat++)
      monomer_atom_names[iat] = atom_info[iat].atom_id;
   
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
	 if (std::find(monomer_atom_names.begin(), monomer_atom_names.end(), test_atom_name)
	     == monomer_atom_names.end()) {
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



CGraph *
coot::dictionary_residue_restraints_t::make_graph(bool use_hydrogens) const { 

   std::map<std::string, unsigned int> name_map;

   
   // mol atom indexing -> graph vertex indexing (we need this because
   // (perhaps) not all atoms in the atom_info are atoms in the graph).
   //
   // This contains the atom indices of of atom_info (0-indexed).  An
   // mmdb graph is 1-indexed.
   // 
   int vertex_indexing[atom_info.size()];
   
   CGraph *graph = new CGraph;
   int i_atom = 0;
   for (unsigned int iat=0; iat<atom_info.size(); iat++) { 
      std::string ele  = atom_info[iat].type_symbol;
      if (use_hydrogens || (ele != " H" && ele != "H")) { 
	 std::string name = atom_info[iat].atom_id_4c;
	 CVertex *v = new CVertex(ele.c_str(), name.c_str());
	 graph->AddVertex(v);
	 name_map[name] = i_atom;
	 i_atom++;
      }
   }

   for (unsigned int ib=0; ib<bond_restraint.size(); ib++) {
      const dict_bond_restraint_t &br = bond_restraint[ib];
      int mmdb_bond_type = br.mmdb_bond_type();
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
		  CEdge *e = new CEdge(it_1->second + 1,  // 1-indexed
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




// 1mzt
