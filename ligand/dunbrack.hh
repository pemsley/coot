/* src/dunbrack.hh
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

#ifndef DUNBRACK_HH
#define DUNBRACK_HH

#include "rotamer.hh"

namespace coot {


   class dunbrack : public rotamer {

      // void add_all_rotamers();  // an autogen function

      contact_info getcontacts(const atom_selection_container_t &asc) const;

      float d2rad(float degrees) const;
      std::vector<std::vector<std::string> >
      rotamer_atoms(const std::string &residue_name) const; 
      double probability_score(double chi_angle, int ichi, const coot::simple_rotamer &rot);


      // for penultime rotamer library
      static short int is_a_residue_name(const std::string &line_part);
      static short int end_of_a_rotamer_p(const std::vector<std::string> &parts);
      static std::string convert_residue_name(const std::string &name_in);
      simple_rotamer parse_prl_rotamer_line(const std::string &line, const std::vector<std::string> &line_parts);

   public:
      // We must be passed a deep copy of residue, the constructor
      // only copies the pointer.
      //
      // We only copy the mol pointer.  Why is it even needed you
      // might ask...
      // 
      // Well (sigh) it's needed in the calculation of the bonds which
      // uses mol->SeekContacts, even though we have a perfectly good
      // atom selection. There should be a version of SeekContacts
      // that does not need a CMMDBManager... Grumble grumble...
      //
      dunbrack(CResidue *residue,
	       const std::string &alt_conf_in,
	       CMMDBManager *mol,
	       float lowest_probability) :
	 rotamer(residue, alt_conf_in, 0) {
	 set_probability_limit(lowest_probability);
	 stored_mol = mol;
      	 // add_all_rotamers();
	 // setup_chi_atom_pairs(); // in dunbrack at least,
		 		 // setup_chi_atom_pairs should happen
				 // after add_all_rotamers().
      }

      dunbrack(CResidue *residue,
	       const std::string &alt_conf_in,
	       CMMDBManager *mol,
	       float lowest_probability,
	       short int add_extra_PHE_and_TYR_rotamers_flag) :
	 rotamer(residue, alt_conf_in, add_extra_PHE_and_TYR_rotamers_flag) {
	 set_probability_limit(lowest_probability);
	 stored_mol = mol;
      	 // add_all_rotamers();
	 // setup_chi_atom_pairs(); // in dunbrack at least,
		 		 // setup_chi_atom_pairs should happen
				 // after add_all_rotamers().
      }

            

      // For use with Z-score (which is analysis only: we don't move anything)
      // 
      dunbrack(CResidue *residue, const std::string &alt_conf_in) :
	 rotamer(residue, alt_conf_in, 0) {}
      
      // Return NULL if no residues available for this residue type
      // 
      std::vector<float> probabilities() const;

      void info() const;


      // maybe this will need to be a static, or a constructor that
      // gets passed to some setup function... not sure yet.  Or maybe
      // it should return some sort of internal data.
      void read_penultimate_library(const std::string &filename);
   };
}

#endif // DUNBRACK_HH
