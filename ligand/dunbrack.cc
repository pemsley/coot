/* src/dunbrack.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2006 The University of York
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
 
#include <iostream>
#include <fstream>

#include <math.h>
#include <algorithm>
#include "utils/coot-utils.hh"
#include "dunbrack.hh"
#include "ccp4mg-utils/mgtree.h"


std::vector<float>
coot::dunbrack::probabilities() const {

   std::vector<float> p;
   std::string rt = Residue_Type();
   if (rt == "MSE")
      rt = "MET";
   std::vector<coot::simple_rotamer> rots = rotamers(rt, Probability_limit());

   for(unsigned int i=0; i<rots.size(); i++)
      p.push_back(rots[i].P_r1234());
   
   return p;
}

float
coot::dunbrack::d2rad(float degrees) const {

   return degrees*M_PI/180.0;
} 



void
coot::dunbrack_rotamer::add_simple_rotamer(coot::simple_rotamer const &rot) {
   rotamers.push_back(rot);
}


coot::dunbrack_rotamer::dunbrack_rotamer(const std::string &restype,
					 const simple_rotamer &rot) {
   residue_type = restype;
   rotamers.push_back(rot);
}


const float &
coot::simple_rotamer::operator[](int i) const {

   switch(i) {
   case 0:
      return chi1;
      break;
   case 1:
      return chi2;
      break;
   case 2:
      return chi3;
      break;
   case 3:
      return chi4;
      break;

   default:
      std::cout << "no such rotatable bond as " << i
		<< " in the dunbrack library\n";
      return minus_one;
   }
   return minus_one; // does not happen?
} 


void
coot::dunbrack::info() const {

   for (unsigned int i=0; i<typed_rotamers.size(); i++) {
      std::cout << i << "  " << typed_rotamers[i].Type()
		<< " " << typed_rotamers[i].n_rotamers() << std::endl;
   }
}



// now in monomer-utils.
// void
// coot::dunbrack_rotamer::add_torsion_bond_by_name(const std::string &atom_name_1,
// 						 const std::string &atom_name_2) {
//    atom_name_pair_list.push_back(coot::atom_name_pair(atom_name_1,
// 						      atom_name_2));
// }




std::vector<coot::simple_rotamer>
coot::dunbrack_rotamer::get_sorted_rotamers(float prob_cut) const {
   std::vector<coot::simple_rotamer> rots;
   for(unsigned int i=0; i< rotamers.size(); i++) {
      if (rotamers[i].P_r1234() > prob_cut) {
	 rots.push_back(rotamers[i]);
      }
   }

   std::sort(rots.begin(), rots.end(), coot::dunbrack_rotamer::compare_rotamers);

   return rots;
}


short int
coot::dunbrack_rotamer::compare_rotamers(const coot::simple_rotamer &a,
					 const coot::simple_rotamer &b) {

   return a.P_r1234() > b.P_r1234();
}


coot::contact_info
coot::dunbrack::getcontacts(const atom_selection_container_t &asc) const {

   mmdb::Contact *pscontact = NULL;
   int n_contacts;
   float min_dist = 0.1;
   float max_dist = 1.9; // CB->SG CYS 1.8A
   if (Residue_Type() == "MSE")
      max_dist = 2.0;
   long i_contact_group = 1;
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   asc.mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms,
			 asc.atom_selection, asc.n_selected_atoms,
			 min_dist, max_dist, // min, max distances
			 0,        // seqDist 0 -> in same res also
			 pscontact, n_contacts,
			 0, &my_matt, i_contact_group);

   return contact_info(pscontact, n_contacts);
} 




					    



void
coot::dunbrack::read_penultimate_library(const std::string &filename) {

   std::ifstream f(filename.c_str());
   int streamsize = 1024-1;
   char line[1024];

   if (f) {

      short int next_is_a_rotamer_flag = 0;
      std::string current_residue_type = "UNASSIGNED RESIDUE";
      while (!f.eof()) {
	 f.getline(line, streamsize);
	 std::vector<std::string> line_parts = coot::util::split_string(line, " ");
	 
	 if (end_of_a_rotamer_p(line_parts)) {
	    next_is_a_rotamer_flag = 0;
	    // std::cout << "ending a rotamer on:\n" << line << std::endl;
	 }

	 if (next_is_a_rotamer_flag) {
	    // std::cout << line << std::endl;
	    coot::simple_rotamer rotamer_info = parse_prl_rotamer_line(line, line_parts);
	    std::cout << current_residue_type << " " << rotamer_info.P_r1234() << std::endl;
	 }

	 // Like awking, this part must come after the test for next_is_a_rotamer_flag
	 // 
	 // std::cout << line_parts.size() << ":" << line << std::endl;
	 if (line_parts.size() > 0) {
	    if (is_a_residue_name(line_parts[0])) {
	       // std::cout << line_parts[0] << std::endl;
	       next_is_a_rotamer_flag = 1;
	       current_residue_type = convert_residue_name(line_parts[0]);
	       // std::cout << "starting a rotamer on:\n" << line << std::endl;
	    }
	 }
      }
   }
}

// static
short int
coot::dunbrack::is_a_residue_name(const std::string &line_part) {

   short int r = 0;
   std::vector<std::string> r_names;
   r_names.push_back("Arginine");
   r_names.push_back("Lysine");
   r_names.push_back("Methionine");
   r_names.push_back("Glutamate");
   r_names.push_back("Glutamine");
   r_names.push_back("Aspartate");
   r_names.push_back("Asparagine");
   r_names.push_back("Isoleucine");
   r_names.push_back("Leucine");
   r_names.push_back("Histidine");
   r_names.push_back("Tryptophan");
   r_names.push_back("Tyrosine");
   r_names.push_back("Phenylalanine");
   r_names.push_back("Proline");
   r_names.push_back("Threonine");
   r_names.push_back("Valine");
   r_names.push_back("Serine");
   r_names.push_back("Cysteine");
   r_names.push_back("Disulfideg");

   for (unsigned int i=0; i<r_names.size(); i++) {
      if (r_names[i] == line_part)
	 return 1;
   } 
   return r;
}

//static
std::string
coot::dunbrack::convert_residue_name(const std::string &name_in) {

   std::vector<std::pair<std::string, std::string> > r_names;
   r_names.push_back(std::pair<std::string, std::string>("Arginine", "ARG"));
   r_names.push_back(std::pair<std::string, std::string>("Lysine", "LYS"));
   r_names.push_back(std::pair<std::string, std::string>("Methionine", "MET"));
   r_names.push_back(std::pair<std::string, std::string>("Glutamate", "GLU"));
   r_names.push_back(std::pair<std::string, std::string>("Glutamine", "GLN"));
   r_names.push_back(std::pair<std::string, std::string>("Aspartate", "ASP"));
   r_names.push_back(std::pair<std::string, std::string>("Asparagine", "ASN"));
   r_names.push_back(std::pair<std::string, std::string>("Isoleucine", "ILE"));
   r_names.push_back(std::pair<std::string, std::string>("Leucine", "LEU"));
   r_names.push_back(std::pair<std::string, std::string>("Histidine", "HIS"));
   r_names.push_back(std::pair<std::string, std::string>("Tryptophan", "TRP"));
   r_names.push_back(std::pair<std::string, std::string>("Tyrosine", "TYR"));
   r_names.push_back(std::pair<std::string, std::string>("Phenylalanine", "PHE"));
   r_names.push_back(std::pair<std::string, std::string>("Proline", "PRO"));
   r_names.push_back(std::pair<std::string, std::string>("Threonine", "THR"));
   r_names.push_back(std::pair<std::string, std::string>("Valine", "VAL"));
   r_names.push_back(std::pair<std::string, std::string>("Serine", "SER"));
   r_names.push_back(std::pair<std::string, std::string>("Cysteine", "CYS"));
   //   r_names.push_back(std::pair<std::string, std::string>("Disulfideg");

   std::string r;
   for (unsigned int i=0; i<r_names.size(); i++) {
      if (name_in == r_names[i].first) {
	 r = r_names[i].second;
	 break;
      }
   }
   return r;
}

// static
short int
coot::dunbrack::end_of_a_rotamer_p(const std::vector<std::string> &parts) {

   short int r=0;
   if (parts.size() > 0) {

//       for (unsigned int i=0; i<parts.size(); i++)
// 	 std::cout << "Part " << i << "/" << parts.size() << " :"
// 		   << parts[i] << ":" << std::endl;

      // find the first non-blank element of parts
      int inb = 0;
      for (unsigned int i=0; i<parts.size(); i++)
	 if (parts[i] != "") {
	    inb = i;
	    break;
	 }
      
      std::vector<std::string> first_word_parts = coot::util::split_string(parts[inb], "%");
      if (first_word_parts.size() == 2) {
	 std::cout << "end: end of a rotamer on " << parts[0] << std::endl;
	 r = 1;
      }
   }

   return r;
}

coot::simple_rotamer
coot::dunbrack::parse_prl_rotamer_line(const std::string &line, const std::vector<std::string> &line_parts_in) {

   int rot1 = 0;
   int rot2 = 0;  
   int rot3 = 0;  
   int rot4 = 0;  
   int n_r1 = 0;
   int nr1234 = 0;
   float p_r1234 = 0.0; // critical value to assign
   float sig_p_r1234 = 0.0;
   float pr234_given_r1 = 0.0;
   float sig_pr234_given_r1 = 0.0;
   float chi1 = -999;
   float sig_chi1 = 0.0;
   float chi2 = -999;
   float sig_chi2 = 0.0;
   float chi3 = -999;
   float sig_chi3 = 0.0;
   float chi4 = -999;
   float sig_chi4 = 0.0;


   // remove blank elements of line_parts
   std::vector<std::string> line_parts;
   for (unsigned int i=0; i<line_parts_in.size(); i++)
      if (line_parts_in[i] != "")
	 line_parts.push_back(line_parts_in[i]);

   if (line_parts.size() > 6) {
      std::string percentage_str = line_parts[2];
      std::vector<std::string> percent_parts = coot::util::split_string(percentage_str, "%");
      int percent_number_parts_index = 0; // because we have to also parse percentage out of "<1%c"
      if (percent_parts.size() > 1) {
	 // float percentage = atof(percent_parts[percent_number_parts_index].c_str());
	 
      } else {
	 std::cout << "Ooops - can't find rotamer percentage" << std::endl;
      }
   }
				
   coot::simple_rotamer rotamer(rot1, rot2, rot3, rot4, n_r1, nr1234, p_r1234, sig_p_r1234,
				pr234_given_r1, sig_pr234_given_r1, chi1, sig_chi1,
				chi2, sig_chi2, chi3, sig_chi3, chi4, sig_chi4);

   return rotamer;
   
}
