/* coot-utils/coot-shelx-ins.hh
 * 
 * Copyright 2005, 2006, 2007 The University of York
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

#ifndef COOT_SHELX_HH
#define COOT_SHELX_HH


#ifndef HAVE_STRING
#include <string>
#define HAVE_STRING
#endif
#ifndef HAVE_VECTOR
#include <vector>
#define HAVE_VECTOR
#endif

#include <map>


#include <mmdb2/mmdb_manager.h>
#include "clipper/core/cell.h"

namespace coot {

   class shelx_card_info_t {
   public:
      shelx_card_info_t() {
	 spaced_start = 0;
      }
      std::string card;
      std::vector<std::string> words;
      short int spaced_start; 
      bool last_word_is_equal_symbol() const;
      void add_card(const shelx_card_info_t &other) {

	 // only add_card if the first 4 characters of the line are not spaces
	 short int add_it = 0;
	 if (other.card.length() > 3) {
	    if (other.card.substr(0,4) == "    ") {
	       add_it = 1; // it was a proper continuation line (according to GMS).
	    } else {
	       add_it = 0;
	    } 
	 } else {
	    add_it = 0;
	 }
	 if (add_it == 1) { 
	    card += other.card;
	    if (words.size() > 0) { 
	       if (words.back() == "=")
		  words.resize(words.size()-1);
	    }
	    for (unsigned int i=0; i<other.words.size(); i++) {
	       words.push_back(other.words[i]);
	    }
	 }
      }
      // return -1 on no bang, else return the index of the ! (bang).
      int bang_index() const;
      void strip_post_bang(); // removed commented lines
      void empty_yorself() {
	 card = "";
	 words.clear();
      } 
   };

   // This is a class for handling symmetry for/from shelx.
   //
   // Shelx has it's symmetry in the for 1, -2, 3 where if positive,
   // means centrosymmetric and negative means non-centrosymmetric.
   // 
   // We need to make symm cards for centering operations
   // (symm_cards_from_lat)
   // 
   class symm_card_composition_t {
      std::string fract_trans_to_str(int itrans_frac) const;
   public:
      int x_element[3], y_element[3], z_element[3];
      int frac_trans[3];
      symm_card_composition_t() {};
      // e.g. "X,Y,Z-1/2"
      explicit symm_card_composition_t(const std::string &symm_card);
      // e.g for F or I centering
      void add_centering_frac(int x_element_in,
			      int y_element_in,
			      int z_element_in);
      void invert();
      float trans_frac(unsigned int index) const {
	 return frac_trans[index]/12.0;
      } 
      std::vector<std::string> symm_cards_from_lat(int latt);
      std::string symm_card() const;
      friend std::ostream & operator<<(std::ostream &s, const symm_card_composition_t &sc);
   };
   std::ostream & operator<<(std::ostream &s, const symm_card_composition_t &sc);

   class shelx_rtab_chi_info_t {
   public:
      int chi;
      std::string resname;
      std::vector<std::string> chi_atoms;
   };

   class shelx_sump_info_t {
   public:
      int v1;
      float v2, v3;
      int v4;
      float v5;
      int v6;
      float v7;
      int v8;
      shelx_sump_info_t(int v1_in,
			float v2_in, float v3_in,
			int v4_in,
			float v5_in,
			int v6_in,
			float v7_in,
			int v8_in) {
	 v1 = v1_in;
	 v2 = v2_in;
	 v3 = v3_in;
	 v4 = v4_in;
	 v5 = v5_in;
	 v6 = v6_in;
	 v7 = v7_in;
	 v8 = v8_in;
      }
      shelx_sump_info_t() {
	 v1 = 1;
	 v2 = 0.01;
	 v3 = 1.0;
	 v4 = 14;
	 v5 = 1.0;
	 v6 = 15;
	 v7 = 1.0;
	 v8 = 16;
      }
   };

   class shelx_read_file_info_t {
   public:
      int udd_afix_handle;
      int status;
      bool is_protein_flag;
      mmdb::Manager *mol;
      shelx_read_file_info_t(int a, int b, mmdb::Manager *mol_in) {
	 status = a;
	 udd_afix_handle = b;
	 mol = mol_in;
	 is_protein_flag = 0;
      }
   };

   // to capture ISOR and CONN info (we need to know the ranges of HETATMs)
   // 
   class hetatom_range {
   public:
      mmdb::Atom *range_first;
      mmdb::Atom *range_last;
      int resno_offset;
      hetatom_range(mmdb::Atom *range_first_in, mmdb::Atom *range_last_in,
		    int resno_offset_in) {
	 range_first = range_first_in;
	 range_last  = range_last_in;
	 resno_offset = resno_offset_in;
      }
      hetatom_range() {
	 range_first = 0;
	 range_last  = 0;
	 resno_offset = 0;
      }
   };

   class ShelxIns {
      std::string title;
      short int filled_flag; 
      bool have_cell_flag;
      clipper::Cell cell;
      std::vector<std::string> sfac; // atom (types). RESI atom index into this vector
      std::vector<int> unit;
      std::vector<float> defs;
      std::vector<int> cgls;
      std::vector<std::string>  pre_atom_lines;
      std::vector<std::string> post_atom_lines;
      std::vector<float> fvars;
      int shel_1;
      float shel2;
      int fmap;
      int plan_1;
      float plan_2;
      int list;
      short int htab_flag;
      int wpdb;
      std::vector<shelx_rtab_chi_info_t> rtab_chi;
      shelx_sump_info_t sump;
      std::vector<std::string> symm_cards;
      std::vector<std::string> disp_cards;
      void init() { shel_1 = -10; shel2 = -1.0;
	 fmap = -1; plan_1 = -1; plan_2 = -1.0;
	 list = -1;
	 htab_flag = 0;
	 wpdb = -2;
	 have_cell_flag = 0;
	 udd_afix_handle = -1;
	 filled_flag = 0;
      }
      
      std::string make_atom_name(const std::string &atom_name_in,
				 const int &atomic_weight) const; // useless?
      std::string make_atom_name(const std::string &atom_name_in,
				 const std::string &element) const; // used now (20051017)
      std::string make_atom_element(const std::string &atom_name_in,
				    const int &atomic_weight) const;
      // return the success status of the assignment, 1 is OK, 0 is fail.
      bool try_assign_cell(mmdb::Manager *mol);

      bool is_unit_line(const std::string &s) const {
	 bool r = 0;
	 if (s.length() >= 4) {
	    if (s.substr(0, 4) == "UNIT") {
	       r = 1;
	    }
	 }
	 return r;
      }
      void write_sfac_line(std::ostream &f) const;
      void write_disp_lines(std::ostream &f) const;

      shelx_card_info_t read_line(std::ifstream &f);
      shelx_card_info_t read_card(std::ifstream &f);
      shelx_card_info_t read_card_extended(std::ifstream &f);
      int get_sfac_index(const std::string &element) const; // -1 on not found
      int altloc_to_part_no(const std::string &altloc) const;

      // As you add the atoms, check the residue name for being
      // standard, and use that to decide if the atom is Het or not.
      // And that means that the passed atom_vector needs to be a
      // reference (or a pointer, but we don't do that here) and no
      // longer const.
      mmdb::Residue *add_shelx_residue(std::vector<mmdb::Atom *> &atom_vector,
				  const std::string &current_res_name,
				  int &current_res_no) const;
      void save_fvars(const shelx_card_info_t &card);
      mmdb::Atom *make_atom(const coot::shelx_card_info_t &card, const std::string &altconf,
			    int udd_afix_handle,
			    int udd_riding_atom_flag_handle,
			    int udd_riding_atom_negative_u_value_handle,
			    bool have_udd_atoms, int current_afix,
			    clipper::Cell &cell, const std::vector<mmdb::Atom *> &atom_vector) const;

      // not const because we add atoms to sfac.
      std::pair<int, std::string> write_ins_file_internal(mmdb::Manager *mol_in,
							  const std::string &filename,
							  bool mol_is_from_shelx_ins=true);

      void write_orthodox_pre_atom_lines(std::ofstream &f) const;

      // not const because we add atoms to sfac.
      void write_synthetic_pre_atom_lines(mmdb::Manager *mol, std::ofstream &f);
      
      mmdb::Atom *previous_non_riding_atom(const std::vector<mmdb::Atom *> &atom_vector,
				      int udd_non_riding_atom_flag_handle) const {

	 if (atom_vector.size()==0) {
	    return NULL;
	 } else { 
	    for (int i=atom_vector.size()-1;  i>=0; i--) {
	       mmdb::Atom *at = atom_vector[i];
	       int ic = -1; // is riding atom flag
	       if (at->GetUDData(udd_non_riding_atom_flag_handle, ic) == mmdb::UDDATA_Ok) {
		  if (ic==1)
		     return at;
	       } else {
		  // this can happen (currently) when udd_non_riding_atom_flag is not set to 1 in
		  // make_atom() (as happens when it is a riding atom)
		  // std::cout << "ooops UDData get failed" << std::endl;
	       } 
	    }
	    // bad news if we get here.
	    return NULL;
	 }
      }
      std::string message_for_atom(const std::string &in_string, mmdb::Atom *at) const;
      std::map<std::string, unsigned int> get_atomic_contents(mmdb::Manager *mol) const;
      bool mol_needs_shelx_transfer(mmdb::Manager *mol) const;
      
   public:
      ShelxIns() {init(); }
      // pair: status (0: bad), udd_afix_handle (-1 bad)
      shelx_read_file_info_t read_file(const std::string &filename);
      explicit ShelxIns(const std::string &filename);
      // return status and message string
      std::pair<int, std::string> write_ins_file(mmdb::Manager *mol,
						 const std::string &filename,
						 bool mol_is_from_shelx_ins=true);
      // is this real shelx data or an empty holder?
      short int is_filled_p() const { return filled_flag; }
      int add_fvar(float f); // return the shelx index FVAR number for this fvar
      void set_fvar(int fvar_no, float f); // set FVAR number fvar_no to value f
      int udd_afix_handle;
      int udd_riding_atom_negative_u_value_handle;
      void debug() const;
      static int shelx_occ_to_fvar(float shelx_occ); // e.g. return 18 if shelx_occ 181.00,
                                                     // return -1 on a problem.
      void add_pre_atoms_line(const std::string &line) {
	 pre_atom_lines.push_back(line);
      }
      int new_chain_offset;
      void add_sfac(const std::string &ele); // do redundancy checking.
   };

   // Do lattice expansion and possible centro-symmetric expansion
   // and add in X,Y,Z which shelx ins file does not require.
   // 
   std::vector<std::string> clipper_symm_strings(const std::vector<std::string> &symm_vec,
						 int shelx_latt);

// return null on no conversion.
   mmdb::Manager *unshelx(mmdb::Manager *mol);
// return null on no conversion.
   mmdb::Manager *reshelx(mmdb::Manager *mol);

}

#endif // COOT_SHELX_HH
