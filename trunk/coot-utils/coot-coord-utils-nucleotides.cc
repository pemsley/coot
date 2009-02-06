/* coot-utils/coot-coord-utils-nucleotides.cc
 * 
 * Copyright 2008, 2009 by The University of Oxford
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
#include <stdexcept>

#include "coot-utils.hh"
#include "coot-coord-utils.hh"


// Throw an exception if it is not possible to generate pucker info
// 
coot::pucker_analysis_info_t::pucker_analysis_info_t(CResidue *res_p,
						     std::string altconf_in) {

   // The atoms are in the following order C1' C2' C3' C4' O4
   //
   altconf = altconf_in; // save for phosphate distance (if needed).
   
   std::vector<CAtom *> atoms(5);
   std::vector<coot::pucker_analysis_info_t::PUCKERED_ATOM_T> possible_puckers;
   possible_puckers.push_back(coot::pucker_analysis_info_t::C1_PRIME);
   possible_puckers.push_back(coot::pucker_analysis_info_t::C2_PRIME);
   possible_puckers.push_back(coot::pucker_analysis_info_t::C3_PRIME);
   possible_puckers.push_back(coot::pucker_analysis_info_t::C4_PRIME);
   possible_puckers.push_back(coot::pucker_analysis_info_t::O4_PRIME);

   PPCAtom residue_atoms = NULL;
   int n_residue_atoms;
   res_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atm_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (altconf == alt_name) { 
	 if (atm_name == " C1*") atoms[0] = residue_atoms[i];
	 if (atm_name == " C1'") atoms[0] = residue_atoms[i];
	 if (atm_name == " C2*") atoms[1] = residue_atoms[i];
	 if (atm_name == " C2'") atoms[1] = residue_atoms[i];
	 if (atm_name == " C3*") atoms[2] = residue_atoms[i];
	 if (atm_name == " C3'") atoms[2] = residue_atoms[i];
	 if (atm_name == " C4*") atoms[3] = residue_atoms[i];
	 if (atm_name == " C4'") atoms[3] = residue_atoms[i];
	 if (atm_name == " O4*") atoms[4] = residue_atoms[i];
	 if (atm_name == " O4'") atoms[4] = residue_atoms[i];
      }
   }
   if (! (atoms[0] && atoms[1] && atoms[2] && atoms[3] && atoms[4])) {
      std::string mess = "Not all atoms found in ribose.";
      throw std::runtime_error(mess);
   } else {
      for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
	 clipper::Coord_orth c(atoms[i_oop_atom]->x,
			       atoms[i_oop_atom]->y,
			       atoms[i_oop_atom]->z);
	 ribose_atoms_coords.push_back(c);
      }
      // oop: out of plane distance
      std::vector<std::pair<float, float> > pucker_distortion_and_oop_d(5);
      for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
	 std::vector<CAtom *> plane_atom;
	 std::vector<clipper::Coord_orth> plane_atom_coords;
	 for (int i=0; i<5; i++) {
	    if (i != i_oop_atom) {
	       clipper::Coord_orth c(atoms[i]->x, atoms[i]->y, atoms[i]->z);
	       plane_atom.push_back(atoms[i]);
	       plane_atom_coords.push_back(c);
	    }
	 }
	 // plane atom is now filled with 4 atoms from which the plane
	 // should be calculated.
	 clipper::Coord_orth pt(atoms[i_oop_atom]->x,
				atoms[i_oop_atom]->y,
				atoms[i_oop_atom]->z);
	 // lsq_plane_deviation returns pair(out-of-plane-dist, rms_deviation_plane);
	 std::pair<double, double> dev =
	    coot::lsq_plane_deviation(plane_atom_coords, pt);
	 pucker_distortion_and_oop_d[i_oop_atom] = dev;
      }

      // Find the biggest out-of-plane distance.  That is the pucker
      // of this ribose.
      puckered_atom_ = coot::pucker_analysis_info_t::NONE;
      std::pair<float, float> most_deviant(0,0);
      for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
// 	 std::cout << "   pucker_distortion_and_oop_d["
// 		   << i_oop_atom << "] "
// 		   << pucker_distortion_and_oop_d[i_oop_atom].first  << " " 
// 		   << pucker_distortion_and_oop_d[i_oop_atom].second  << std::endl;
	 if (fabs(pucker_distortion_and_oop_d[i_oop_atom].first) > fabs(most_deviant.first)) {
	    most_deviant = pucker_distortion_and_oop_d[i_oop_atom];
	    puckered_atom_ = possible_puckers[i_oop_atom];
	 }
      }
      out_of_plane_distance = most_deviant.first;
      plane_distortion = most_deviant.second;
   }
}

// Use the 3' phosphate of the following residue to calculate
// its out of plane distance.  Decide from that if this should
// have been 3' or 2'.  Check vs the actual puckering.
//
// Throw an exception if we can't do this.
// 
float 
coot::pucker_analysis_info_t::phosphate_distance(CResidue *following_res) {

   float oop = 0.0;
   PPCAtom residue_atoms = NULL;
   int n_residue_atoms;
   bool found = 0;
   following_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atm_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (atm_name == " P  ") { 
	 if (altconf == alt_name) {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    // lsq_plane_deviation returns pair(out-of-plane-dist, rms_deviation_plane);
	    std::pair<double, double> oop_plus_dev =
	       coot::lsq_plane_deviation(ribose_atoms_coords, pt);
	    oop = oop_plus_dev.first;
	    found = 1;
	    break;
	 }
      }
   }
   if (found == 0) {
      throw std::runtime_error("Failed to find following phosphate"); 
   } 
   return oop;
}

std::string
coot::pucker_analysis_info_t::puckered_atom() const {

   std::string s;
   if (puckered_atom_ == coot::pucker_analysis_info_t::C2_PRIME)
      s = " C2'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::C3_PRIME)
      s = " C3'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::NONE)
      s = "----";
   if (puckered_atom_ == coot::pucker_analysis_info_t::C1_PRIME)
      s = " C1'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::C4_PRIME)
      s = " C4'";
   if (puckered_atom_ == coot::pucker_analysis_info_t::O4_PRIME)
      s = " O4'";
   return s;
} 


// return "" on no canonical name found
std::string
coot::util::canonical_base_name(const std::string res_name_in, base_t rna_or_dna) {

   if (rna_or_dna == coot::RNA) {
      if (res_name_in == "C")
	 return "Cr";
      if (res_name_in == "A")
	 return "Ar";
      if (res_name_in == "G")
	 return "Gr";
      if (res_name_in == "T")
	 return "Tr";
      if (res_name_in == "U")
	 return "Ur";
      if (res_name_in == "Cr")
	 return "Cr";
      if (res_name_in == "Ar")
	 return "Ar";
      if (res_name_in == "Gr")
	 return "Gr";
      if (res_name_in == "Tr")
	 return "Tr";
      if (res_name_in == "Ur")
	 return "Ur";
      if (res_name_in == "Cd")
	 return "Cr";
      if (res_name_in == "Ad")
	 return "Ar";
      if (res_name_in == "Gd")
	 return "Gr";
      if (res_name_in == "Td")
	 return "Tr";
      if (res_name_in == "Ud")
	 return "Ur";
   }
   
   if (rna_or_dna == coot::RNA) {
      if (res_name_in == "C")
	 return "Cd";
      if (res_name_in == "A")
	 return "Ad";
      if (res_name_in == "G")
	 return "Gd";
      if (res_name_in == "T")
	 return "Td";
      if (res_name_in == "U")
	 return "Ud";
      if (res_name_in == "Cd")
	 return "Cd";
      if (res_name_in == "Ad")
	 return "Ad";
      if (res_name_in == "Gd")
	 return "Gd";
      if (res_name_in == "Td")
	 return "Td";
      if (res_name_in == "Ud")
	 return "Ud";
      if (res_name_in == "Cr")
	 return "Cd";
      if (res_name_in == "Ar")
	 return "Ad";
      if (res_name_in == "Gr")
	 return "Gd";
      if (res_name_in == "Tr")
	 return "Td";
      if (res_name_in == "Ur")
	 return "Ud";
   }

   return "";
}
