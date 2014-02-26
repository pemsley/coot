
#include "utils/coot-utils.hh"
#include "protein-geometry.hh"

// quite means don't tell me about matches
bool
coot::dictionary_residue_restraints_t::compare(const dictionary_residue_restraints_t &r,
					       bool quiet) const {

   bool check_hydrogens = false;
   double bond_length_tolerance = 0.01;
   double bond_esd_tolerance    = 0.005;
   double angle_tolerance       = 0.3;
   double angle_esd_tolerance   = 0.15;

   
   // residue info
   bool residue_info_matches = true;
   if (r.residue_info.comp_id != residue_info.comp_id) {
      std::cout << "Residue-Info:: mismatch comp_id "
		<< r.residue_info.comp_id << " " << residue_info.comp_id << std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.three_letter_code != residue_info.three_letter_code) {
      std::cout << "Residue-Info:: mismatch three_letter_code "
		<< r.residue_info.three_letter_code << " " << residue_info.three_letter_code
		<< std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.name != residue_info.name) {
      std::cout << "Residue-Info:: mismatch name "
		<< r.residue_info.name << " " << residue_info.name << std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.group != residue_info.group) {
      std::cout << "Residue-Info:: mismatch group "
		<< r.residue_info.group << " " << residue_info.group << std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.number_atoms_all != residue_info.number_atoms_all) {
      std::cout << "Residue-Info:: mismatch number_atoms_all "
		<< r.residue_info.number_atoms_all << " " << residue_info.number_atoms_all
		<< std::endl;
      residue_info_matches = false;
   }
   if (r.residue_info.number_atoms_nh != residue_info.number_atoms_nh) {
      std::cout << "Residue-Info:: mismatch number_atoms_nh "
		<< r.residue_info.number_atoms_nh << " " << residue_info.number_atoms_nh
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
	       std::cout << "Atom-Info:: " << atom_id_refr << " type_energy mismatch "
			 << atom_info[iat].type_energy << " vs "
			 << r.atom_info[jat].type_energy << " "
			 << std::endl;
	    }
	 }
      }
   }

   // -------------------------  bond restraints -------------------
   //
   if (bond_restraint.size() != r.bond_restraint.size())
      std::cout << "Bond-Restraint:: mismatch number of restraints "
		<< bond_restraint.size() << " " << r.bond_restraint.size()
		<< std::endl;
   bool bonds_match = true;
   for (unsigned int ib=0; ib<bond_restraint.size(); ib++) { 
      for (unsigned int jb=0; jb<r.bond_restraint.size(); jb++) {
	 if (bond_restraint[ib].matches_names(r.bond_restraint[jb])) {
	    double d_d = bond_restraint[ib].value_dist() -
	       r.bond_restraint[jb].value_dist();
	    double d_d_a = fabs(d_d);
	    if (d_d_a > bond_length_tolerance) {
	       std::cout << "Bond-Restraint:: mismatch bond length between "
			 << bond_restraint[ib].atom_id_1() << " "
			 << bond_restraint[ib].atom_id_2() << " "
			 << bond_restraint[ib].value_dist() << " vs "
			 << r.bond_restraint[jb].value_dist()
			 << std::endl;
	       bonds_match = false;
	    }
	    double esd_d = bond_restraint[ib].value_esd() -
	       r.bond_restraint[jb].value_esd();
	    double esd_d_a = fabs(esd_d);
	    if (esd_d_a > bond_esd_tolerance) {
	       std::cout << "Bond-Restraint:: mismatch bond length esd between "
			 << bond_restraint[ib].atom_id_1() << " "
			 << bond_restraint[ib].atom_id_2() << " "
			 << bond_restraint[ib].value_esd() << " vs "
			 << r.bond_restraint[jb].value_esd()
			 << std::endl;
	       bonds_match = false;
	    }
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
   bool angles_match = true;
   for (unsigned int ia=0; ia<angle_restraint.size(); ia++) { 
      for (unsigned int ja=0; ja<r.angle_restraint.size(); ja++) {
	 if (angle_restraint[ia].matches_names(r.angle_restraint[ja])) {
	    double a_d = angle_restraint[ia].angle() -
	       r.angle_restraint[ja].angle();
	    double a_d_a = fabs(a_d);
	    if (a_d_a > angle_tolerance) {
	       std::cout << "Angle-Restraint:: mismatch angle between "
			 << angle_restraint[ia].atom_id_1() << " "
			 << angle_restraint[ia].atom_id_2() << " "
			 << angle_restraint[ia].atom_id_3() << " "
			 << angle_restraint[ia].angle() << " vs "
			 << r.angle_restraint[ja].angle()
			 << std::endl;
	       angles_match = false;
	    }
	    double esd_d = angle_restraint[ia].esd() -
	       r.angle_restraint[ja].esd();
	    double esd_d_a = fabs(esd_d);
	    if (esd_d_a > angle_esd_tolerance) {
	       std::cout << "Angle-Restraint:: mismatch angle esd "
			 << angle_restraint[ia].atom_id_1() << " "
			 << angle_restraint[ia].atom_id_2() << " "
			 << angle_restraint[ia].atom_id_3() << " "
			 << angle_restraint[ia].esd() << " vs "
			 << r.angle_restraint[ja].esd()
			 << std::endl;
	       angles_match = false;
	    }
	 }
      }
   }
   if (angles_match)
      if (! quiet)
	 std::cout << "Angle-Restraint::  all angles match within tolerance "
		   << std::endl;

   // -------------------------  torsion restraints -------------------

   // -------------------------  chiral restraints -------------------
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
      planes_match = true;
   //
   if (planes_match) {
      if (! quiet)
	 std::cout << "Plane-Restraint::  all plane restraints match" << std::endl;
   } else {
      std::cout << "Plane-Restraint:: plane restraints do not match" << std::endl;
   }
   

   bool status = false;
   if (bonds_match == true)
      if (angles_match == true)
	 if (planes_match == true)
	    if (chirals_match == true)
	       status = true;
   
   return status;
} 
