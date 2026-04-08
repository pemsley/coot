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

#include <optional>
#include <stdexcept>

#include "clipper/core/coords.h"
#include "geometry/residue-and-atom-specs.hh"
#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"

 
// Throw an exception if it is not possible to generate pucker info
// 
coot::pucker_analysis_info_t::pucker_analysis_info_t(mmdb::Residue *res_p,
                                                     std::string altconf_in) {

   auto get_base_lsq_plane = [] (const std::vector<clipper::Coord_orth> &coords) -> 
      std::optional<lsq_plane_info_t> {

      if (coords.size() < 3) return std::nullopt;
      lsq_plane_info_t lsq_plane(coords);
      return lsq_plane;
   };

   out_of_plane_distance = 0.0;
   plane_distortion = 0.0;
   C1_prime = nullptr;
   N1_or_9  = nullptr;
   // The atoms are in the following order C1' C2' C3' C4' O4
   //
   altconf = altconf_in; // save for phosphate distance (if needed).

   assign_base_atom_coords(res_p); // and C1_prime and N1_or_9 if possible
   std::optional<lsq_plane_info_t> lsq_plane = get_base_lsq_plane(base_atoms_coords);

   if (lsq_plane.has_value()) {

      // store the geometry
      markup_info.base_ring_centre = lsq_plane.value().centre();
      markup_info.base_ring_normal = lsq_plane.value().centre();

      std::vector<mmdb::Atom *> ribose_atoms(5, nullptr); // ribose atoms
      std::vector<coot::pucker_analysis_info_t::PUCKERED_ATOM_T> possible_puckers;
      possible_puckers.push_back(coot::pucker_analysis_info_t::C1_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::C2_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::C3_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::C4_PRIME);
      possible_puckers.push_back(coot::pucker_analysis_info_t::O4_PRIME);

      mmdb::PPAtom residue_atoms = NULL;
      int n_residue_atoms = 0;
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      // find the phosphorus atom
      for (int i=0; i<n_residue_atoms; i++) {
         mmdb::Atom *atm = residue_atoms[i];
         if (! atm->isTer()) {
            std::string atm_name(atm->name);
            std::string alt_name(atm->altLoc);
            if (altconf == alt_name) {
               if (atm_name == " P  ") { // PDBv3 FIXME
                  clipper::Coord_orth p(atm->x, atm->y, atm->z);
                  markup_info.phosphorus_position = p;
                  markup_info.projected_point = lsq_plane.value().projected_point(p);
               }
            }
         }
      }
      // find the ribose atoms
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         std::string atm_name(residue_atoms[i]->name);
         std::string alt_name(residue_atoms[i]->altLoc);
         if (altconf == alt_name) {
            if (atm_name == " C1*") ribose_atoms[0] = residue_atoms[i];
            if (atm_name == " C1'") ribose_atoms[0] = residue_atoms[i];
            if (atm_name == " C2*") ribose_atoms[1] = residue_atoms[i];
            if (atm_name == " C2'") ribose_atoms[1] = residue_atoms[i];
            if (atm_name == " C3*") ribose_atoms[2] = residue_atoms[i];
            if (atm_name == " C3'") ribose_atoms[2] = residue_atoms[i];
            if (atm_name == " C4*") ribose_atoms[3] = residue_atoms[i];
            if (atm_name == " C4'") ribose_atoms[3] = residue_atoms[i];
            if (atm_name == " O4*") ribose_atoms[4] = residue_atoms[i];
            if (atm_name == " O4'") ribose_atoms[4] = residue_atoms[i];
         }
      }
      if (! (ribose_atoms[0] && ribose_atoms[1] && ribose_atoms[2] && ribose_atoms[3] && ribose_atoms[4])) {
         std::string mess = "Not all atoms found in ribose.";
         throw std::runtime_error(mess);
      } else {
         for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
            clipper::Coord_orth c(ribose_atoms[i_oop_atom]->x,
                                  ribose_atoms[i_oop_atom]->y,
                                  ribose_atoms[i_oop_atom]->z);
            ribose_atoms_coords.push_back(c);
         }
         // oop: out of plane distance
         std::vector<std::pair<float, float> > pucker_distortion_and_oop_d(5);
         for (int i_oop_atom=0; i_oop_atom<5; i_oop_atom++) {
            std::vector<mmdb::Atom *> plane_atom;
            std::vector<clipper::Coord_orth> plane_atom_coords;
            for (int i=0; i<5; i++) {
               if (i != i_oop_atom) {
                  clipper::Coord_orth c(ribose_atoms[i]->x, ribose_atoms[i]->y, ribose_atoms[i]->z);
                  plane_atom.push_back(ribose_atoms[i]);
                  plane_atom_coords.push_back(c);
               }
            }
            // plane atom is now filled with 4 atoms from which the plane
            // should be calculated.
            clipper::Coord_orth pt(ribose_atoms[i_oop_atom]->x,
                                   ribose_atoms[i_oop_atom]->y,
                                   ribose_atoms[i_oop_atom]->z);
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
            //          std::cout << "   pucker_distortion_and_oop_d["
            //                    << i_oop_atom << "] "
            //                    << pucker_distortion_and_oop_d[i_oop_atom].first  << " "
            //                    << pucker_distortion_and_oop_d[i_oop_atom].second  << std::endl;
            if (fabs(pucker_distortion_and_oop_d[i_oop_atom].first) > fabs(most_deviant.first)) {
               most_deviant = pucker_distortion_and_oop_d[i_oop_atom];
               puckered_atom_ = possible_puckers[i_oop_atom];
            }
         }
         out_of_plane_distance = most_deviant.first;
         plane_distortion = most_deviant.second;
      }
   } else {
      // we throw on failure
      std::string mess = "base lsq plane has no value";
      throw std::runtime_error(mess);
   }
}

#include "json.hpp" // Assumes nlohmann/json.hpp is available
using json = nlohmann::json;

std::string
coot::pucker_analysis_info_t::to_json() const {

   json j;
   json j_altconf = altconf;
   json j_plane_distortion = plane_distortion;
   json j_out_of_plane_distance = out_of_plane_distance;
   json j_puckered_atom = puckered_atom();
   j["altconf"]               = j_altconf;
   j["plane_distortion"]      = j_plane_distortion;
   j["out_of_plane_distance"] = j_out_of_plane_distance;
   j["puckered_atom"]         = j_puckered_atom;

   // json["plane_distortion"] = plane_distortion;
   // json["out_of_plane_distance"] = out_of_plane_distance;
   // json["puckered_atom"] = puckered_atom();
   std::string s = j.dump(4);
   return s;
}



void
coot::pucker_analysis_info_t::assign_base_atom_coords(mmdb::Residue *residue_p) {

   std::vector<std::string> cytidine_base_names;
   std::vector<std::string> uracil_base_names;
   std::vector<std::string> adenine_base_names;
   std::vector<std::string> guanine_base_names;
   std::vector<std::string> thymine_base_names;

   cytidine_base_names.push_back(" N1 ");
   cytidine_base_names.push_back(" C2 ");
   cytidine_base_names.push_back(" N3 ");
   cytidine_base_names.push_back(" C4 ");
   cytidine_base_names.push_back(" C5 ");
   cytidine_base_names.push_back(" C6 ");
   cytidine_base_names.push_back(" O2 ");
   cytidine_base_names.push_back(" N4 ");

   uracil_base_names.push_back(" N1 ");
   uracil_base_names.push_back(" C2 ");
   uracil_base_names.push_back(" N3 ");
   uracil_base_names.push_back(" C4 ");
   uracil_base_names.push_back(" C5 ");
   uracil_base_names.push_back(" C6 ");
   uracil_base_names.push_back(" O2 ");
   uracil_base_names.push_back(" O4 ");

   adenine_base_names.push_back(" N9 ");
   adenine_base_names.push_back(" C8 ");
   adenine_base_names.push_back(" N7 ");
   adenine_base_names.push_back(" C5 ");
   adenine_base_names.push_back(" C4 ");
   adenine_base_names.push_back(" N1 ");
   adenine_base_names.push_back(" C2 ");
   adenine_base_names.push_back(" N3 ");
   adenine_base_names.push_back(" C6 ");
   adenine_base_names.push_back(" N6 ");

   guanine_base_names.push_back(" N9 ");
   guanine_base_names.push_back(" C8 ");
   guanine_base_names.push_back(" N7 ");
   guanine_base_names.push_back(" C5 ");
   guanine_base_names.push_back(" C4 ");
   guanine_base_names.push_back(" N1 ");
   guanine_base_names.push_back(" C2 ");
   guanine_base_names.push_back(" N3 ");
   guanine_base_names.push_back(" C6 ");
   guanine_base_names.push_back(" O6 ");
   guanine_base_names.push_back(" N2 ");

   thymine_base_names.push_back(" N1 ");
   thymine_base_names.push_back(" C2 ");
   thymine_base_names.push_back(" N3 ");
   thymine_base_names.push_back(" C4 ");
   thymine_base_names.push_back(" C5 ");
   thymine_base_names.push_back(" C6 ");
   thymine_base_names.push_back(" O2 ");
   thymine_base_names.push_back(" O4 ");
   thymine_base_names.push_back(" C5M");


   mmdb::PPAtom residue_atoms = NULL;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   // Assign N1_or_9 and C1_prime
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (alt_name == altconf) {
         if (atom_name == " N1 ")
            N1_or_9 = residue_atoms[i];
         if (atom_name == " N9 ")
            N1_or_9 = residue_atoms[i];
         if (atom_name == " C1*")
            C1_prime = residue_atoms[i];
         if (atom_name == " C1'")
            C1_prime = residue_atoms[i];
      }
   }

   // Fill base_names according to residue type/name.  If base_name is
   // empty after setting, just fall out (an exception is thrown in
   // the constructor if there are not enough base name atoms.

   std::vector<std::string> base_names;

   std::string residue_name(residue_p->GetResName());

   // current names
   if (residue_name == "C") base_names = cytidine_base_names;
   if (residue_name == "U") base_names = uracil_base_names;
   if (residue_name == "A") base_names = adenine_base_names;
   if (residue_name == "G") base_names = guanine_base_names;
   // old names
   if (residue_name == "Cr") base_names = cytidine_base_names;
   if (residue_name == "Ur") base_names = uracil_base_names;
   if (residue_name == "Ar") base_names = adenine_base_names;
   if (residue_name == "Gr") base_names = guanine_base_names;
   // modern (3.x) RNA base names
   if (residue_name == "CYT") base_names = cytidine_base_names;
   if (residue_name == "URA") base_names = uracil_base_names;
   if (residue_name == "ADE") base_names = adenine_base_names;
   if (residue_name == "GUA") base_names = guanine_base_names;

   if (base_names.size() > 0) {
      for (int i=0; i<n_residue_atoms; i++) {
         std::string atm_name(residue_atoms[i]->name);
         std::string alt_name(residue_atoms[i]->altLoc);
         for (unsigned int j=0; j<base_names.size(); j++) {
            if (base_names[j] == atm_name) {
               base_atoms_coords.push_back(clipper::Coord_orth(residue_atoms[i]->x,
                                                               residue_atoms[i]->y,
                                                               residue_atoms[i]->z));
            }
         }
      }
   }
} 


// Use the 3' phosphate of the following residue to calculate its out
// of plane distance (the plane being the base plane).  Decide from
// that if this should have been 3' or 2'.  Check vs the actual
// puckering.
//
// Throw an exception if we can't do this.
// 
float 
coot::pucker_analysis_info_t::phosphate_distance_to_base_plane(mmdb::Residue *following_res) {

   float oop = 0.0;
   mmdb::PPAtom residue_atoms = NULL;
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

            if (base_atoms_coords.size() < 4) {

               // construct an error message and throw an exception.
               // 
               std::string m = "Failed to find base atoms. Found ";
               m += coot::util::int_to_string(base_atoms_coords.size());
               m += " atoms. ";
               throw std::runtime_error(m);
            } else { 
               std::pair<double, double> oop_plus_dev =
                  coot::lsq_plane_deviation(base_atoms_coords, pt);
               oop = oop_plus_dev.first;
               found = 1;
               break;
            }
         }
      }
   }
   if (found == 0) {
      throw std::runtime_error("Failed to find following phosphate");
   }
   return oop;
}

// Throw an exception if the reference atoms are not found.
float 
coot::pucker_analysis_info_t::phosphate_distance(mmdb::Residue *following_res) {

   if (! C1_prime) { 
      std::string mess = "C1*/C1' not found in this residue";
      throw std::runtime_error(mess);
   }
   if (! N1_or_9) { 
      std::string mess = "N1/N9 not found in this residue";
      throw std::runtime_error(mess);
   }

   //                     X             
   //                    / \  90 degrees       ;
   //                 d / X \                  ;
   //                  /     \                 ;
   //                 /       \                       ;
   //                /         \               ; 
   //               /   pi-alpha\              ;
   //              --------------\C1'               ;
   //            P           alpha\            ; 
   //                              \           ;
   //                               \          ; 
   //                                \         ;
   //                                 \               ;
   //                                  \       ; 
   //                                   \ N1   ;
   
   
   float d = 0.0;
   mmdb::PPAtom residue_atoms = NULL;
   int n_residue_atoms;
   bool found = 0;
   following_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atm_name(residue_atoms[i]->name);
      std::string alt_name(residue_atoms[i]->altLoc);
      if (atm_name == " P  ") { 
         if (altconf == alt_name) {
            clipper::Coord_orth P_pt(residue_atoms[i]->x,
                                     residue_atoms[i]->y,
                                     residue_atoms[i]->z);
            clipper::Coord_orth N_pt( N1_or_9->x,  N1_or_9->y,  N1_or_9->z);
            clipper::Coord_orth C_pt(C1_prime->x, C1_prime->y, C1_prime->z);
            clipper::Coord_orth CN = N_pt - C_pt;
            clipper::Coord_orth CP = P_pt - C_pt;

            double CN_d = clipper::Coord_orth::length(N_pt, C_pt);
            double CP_d = clipper::Coord_orth::length(P_pt, C_pt);

            if (CN_d > 0.0) { 
               if (CP_d > 0.0) { 
                  found = 1;
                  double cos_alpha = clipper::Coord_orth::dot(CN, CP)/(CN_d*CP_d);
                  
                  double alpha = acos(cos_alpha);
                  double sin_pi_minus_alpha = sin(M_PI - alpha);
                  d = sin_pi_minus_alpha * CP_d;
               }
            }
         }
      }
   }

   if (! found) {
      std::string mess = "P not found in this residue";
      throw std::runtime_error(mess);
   }
   return d;
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
coot::util::canonical_base_name(const std::string &res_name_in, base_t rna_or_dna) {

   if (rna_or_dna == coot::RNA) {
      if (res_name_in == "C")
         return "C";
      if (res_name_in == "A")
         return "A";
      if (res_name_in == "G")
         return "G";
      if (res_name_in == "T")
         return "T";
      if (res_name_in == "U")
         return "U";
      if (res_name_in == "Cr")
         return "C";
      if (res_name_in == "Ar")
         return "A";
      if (res_name_in == "Gr")
         return "G";
      if (res_name_in == "Tr")
         return "T";
      if (res_name_in == "Ur")
         return "U";
      if (res_name_in == "Cd")
         return "C";
      if (res_name_in == "Ad")
         return "A";
      if (res_name_in == "Gd")
         return "G";
      if (res_name_in == "Td")
         return "T";
      if (res_name_in == "Ud")
         return "U";
   }
   
   if (rna_or_dna == coot::DNA) {
      if (res_name_in == "C")
         return "DC";
      if (res_name_in == "A")
         return "DA";
      if (res_name_in == "G")
         return "DG";
      if (res_name_in == "T")
         return "DT";
      if (res_name_in == "U")
         return "DU";
      if (res_name_in == "Cd")
         return "DC";
      if (res_name_in == "Ad")
         return "DA";
      if (res_name_in == "Gd")
         return "DG";
      if (res_name_in == "Td")
         return "DT";
      if (res_name_in == "Ud")
         return "DU";
      if (res_name_in == "Cr")
         return "DC";
      if (res_name_in == "Ar")
         return "DA";
      if (res_name_in == "Gr")
         return "DG";
      if (res_name_in == "Tr")
         return "DT";
      if (res_name_in == "Ur")
         return "DU";
   }

   return "";
}
