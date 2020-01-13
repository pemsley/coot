/* src/chi-angles.hh
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

#ifndef CHI_ANGLES_HH
#define CHI_ANGLES_HH

#include "geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-utils.hh"

#include "monomer-utils.hh"

namespace coot {

   // Probabilities in percentages
   // 
   class simple_rotamer {
      short int rotamer_type; 
      std::string name;
      int rot1;  // rotamer indices
      int rot2;  // 
      int rot3;  //
      int rot4;  //
      int n_r1;
      int nr1234;
      float p_r1234;
      float sig_p_r1234;
      float pr234_given_r1;
      float sig_pr234_given_r1;
      float chi1;
      float sig_chi1;
      float chi2;
      float sig_chi2;
      float chi3;
      float sig_chi3;
      float chi4;
      float sig_chi4;
      float minus_one;
      std::vector<float> chi_vec;  // zero indexed, accessors must be offset.
      std::vector<float> sig_chi_vec; // ditto.

   public:
      simple_rotamer();
      simple_rotamer(int rot1,  
		     int rot2,  
		     int rot3,  
		     int rot4,  
		     int n_r1,
		     int nr1234,
		     float p_r1234,
		     float sig_p_r1234,
		     float pr234_given_r1,
		     float sig_pr234_given_r1,
		     float chi1,
		     float sig_chi1,
		     float chi2,
		     float sig_chi2,
		     float chi3,
		     float sig_chi3,
		     float chi4,
		     float sig_chi4);
      // constructor for richardson rotamer
      simple_rotamer(std::string rotamer_name,  
		     float percent_overall,
		     float percent_alpha,
		     float percent_beta,
		     float percent_other,
		     float chi_1_mode,
		     float chi_1_com,
		     float chi_2_mode,
		     float chi_2_com,
		     float chi_3_mode,
		     float chi_3_com,
		     float chi_4_mode,
		     float chi_4_com);

      enum rotamer_t { RICHARDSON_ROTAMER, DUNBRACK_ROTAMER};
      
      float P_r1234() const { return p_r1234; }
      float Probability_rich() const { return p_r1234; }
      const float & operator[](int i) const;
      float Chi1() const { return chi1; } 
      float Chi2() const { return chi2; } 
      float Chi3() const { return chi3; } 
      float Chi4() const { return chi4; }
      // start at ichi == 1 (not zero indexed)
      float get_chi(int ichi) const {
	 if (ichi == 1) return chi1;
	 else
	    if (ichi == 2) return chi2;
	    else
	       if (ichi == 3) return chi3;
	       else
		  if (ichi == 4) return chi4;
		  else return -999; }
	 

      int N_chi() const { return chi_vec.size(); }
      std::string rotamer_name() const {return name;} // richardson rotamer name (m, tt, p-90)
      
      short int has_chi2_p() const { return sig_chi2 > 0.0; };
      short int has_chi3_p() const { return sig_chi3 > 0.0; };
      short int has_chi4_p() const { return sig_chi4 > 0.0; };
      short int has_chi(int n) const {
	 int sig_chi_size = sig_chi_vec.size();
	 if (n < sig_chi_size)
	    if (sig_chi_vec[n-1] > 0.0)
	       return 1;
	    else
	       return 0;
	 else
	    return 0;
      }

      simple_rotamer rotate_chi2_180() const;
      // This works now with declaration coot::operator<<()
      friend std::ostream& operator<<(std::ostream &s, coot::simple_rotamer rot);
   };
   std::ostream& operator<<(std::ostream &s, coot::simple_rotamer rot);


   // a dunbrack_rotamer is a typed container of simple rotamers,
   // i.e. it is the set of all (e.g.) ARG rotamers.
   // 
   class dunbrack_rotamer : public monomer_utils {

      std::string residue_type;
      std::vector<simple_rotamer> rotamers;

   public:
      dunbrack_rotamer(const std::string &restype,
		       const simple_rotamer &rot);
      
      void add_simple_rotamer(const simple_rotamer &rot);
      std::string Type() const { return residue_type; }
      int n_rotamers() const { return rotamers.size(); }

      static short int compare_rotamers(const simple_rotamer &a,
					const simple_rotamer &b);
      std::vector<simple_rotamer> get_simple_rotamers() const { return rotamers; }
      std::vector<simple_rotamer> get_sorted_rotamers(float prob_cut) const;
   };


   
   class chi_angles {

   protected:

      mmdb::Residue *residue;
      std::string residue_type;
      // look it up in the appropriate dunbrack_rotamer:
      // 
      // On failure, return a vector of size one with pair with
      // pair.first = "empty"
      std::vector<atom_name_pair>
      atom_name_pair_list(const std::string &res_type) const;
      std::vector<coot::atom_index_pair> get_atom_index_pairs(const std::vector<coot::atom_name_pair> &atom_name_pairs,
							      const mmdb::PPAtom atoms, int nresatoms) const;

      std::vector<coot::atom_name_quad> atom_name_quad_list(const std::string &residue_type) const;
      std::vector<coot::atom_index_quad> get_atom_index_quads(const std::vector<coot::atom_name_quad> &atom_name_quads,
							      const mmdb::PPAtom atoms, int nresatoms) const;

      void add_all_rotamers();  // an autogen function
      void add_richardson_rotamers(); 
      void add_rotamer(std::string restype,
		       int rot1,  
		       int rot2,  
		       int rot3,  
		       int rot4,  
		       int n_r1,
		       int nr1234,
		       float p_r1234,
		       float sig_p_r1234,
		       float pr234_given_r1,
		       float sig_pr234_given_r1,
		       float chi1,
		       float sig_chi1,
		       float chi2,
		       float sig_chi2,
		       float chi3,
		       float sig_chi3,
		       float chi4,
		       float sig_chi4);

      // add_richardson_rotamer("ARG", "mmm-85", 22, 2, 2, 3, 3, -62, -62, 0, 0, 0, 0, 0, 0);
      void add_richardson_rotamer(std::string restype,
				  std::string rotamer_name,
				  float anumber,
				  float percent_overall,
				  float percent_alpha,
				  float percent_beta,
				  float percent_other,
				  float chi_1_mode,
				  float chi_1_com, 
				  float chi_2_mode,
				  float chi_2_com, 
				  float chi_3_mode,
				  float chi_3_com, 
				  float chi_4_mode,
				  float chi_4_com);
      void use_richardson_rotamers(); // uses above function and wipes
				      // out Dunbrack rotamers
	 

      // called by the change_bys, obviously.
      // return status and the new angle.
      std::pair<short int, float> change_by_internal(int ichi,
				   double diff,
				   const std::vector<coot::atom_name_pair> &atom_name_pairs,
				   const std::vector<std::vector<int> > &contact_indices,
				   mmdb::PPAtom residue_atoms,
				   int nResidueAtoms,
				   const coot::atom_spec_t &tree_base_atom);

      std::vector<coot::atom_name_pair>
      get_torsion_bonds_atom_pairs(const std::string &monomer_type,
				   int imol,
				   coot::protein_geometry *pg,
				   short int include_hydrogen_torsions_flag) const;

      void add_IUPAC_extras_PHE_and_TYR_rotamers();

   public:

      chi_angles(mmdb::Residue *residue_in, short int add_extra_PHE_and_TYR_rotamers_flag) {
	 residue = residue_in;
	 if (residue)
	    residue_type = residue->GetResName();
#ifdef USE_DUNBRACK_ROTAMERS			
	 add_all_rotamers();
#else
	 add_richardson_rotamers();
#endif // USE_DUNBRACK_ROTAMERS			
	 // debugging
// 	 if (add_extra_PHE_and_TYR_rotamers_flag)
// 	    std::cout << "+++ Adding in extra PHE/TYR\n";
// 	 else 
// 	    std::cout << "Not Adding in extra PHE/TYR\n";
	 if (add_extra_PHE_and_TYR_rotamers_flag)
	    add_IUPAC_extras_PHE_and_TYR_rotamers();
	 setup_chi_atom_quads();
      }

      mmdb::Residue *Residue() const { return residue; } 
      std::string Residue_Type() const { return residue_type; }
      void setup_chi_atom_quads();
      void add_chi_pair(const std::string &residue_type,
			const std::string &atom_name_1,
			const std::string &atom_name_2);
      void add_chi_quad(const std::string &residue_type,
			const std::string &atom_name_1,
			const std::string &atom_name_2,
			const std::string &atom_name_3,
			const std::string &atom_name_4);

      // std::vector<simple_rotamer> get_simple_rotamers(const std::string &res_type, float prob_cut) const;

#ifdef USE_DUNBRACK_ROTAMERS			
      std::vector<dunbrack_rotamer> typed_rotamers;
#else
      std::vector<dunbrack_rotamer> typed_rotamers; // may change in future?
#endif // USE_DUNBRACK_ROTAMERS			


      // Return success status, 
      // 0 means success
      // 
      // 1: we failed because we didn't find the residue type in
      // typed_rotamers
      // 
      // 2: ichi was wrong [ichi is not zero indexed].
      // 
      // 3: ... 
      // 
      // return status and the new angle.
      std::pair<short int, float> change_by(int ichi, double diff,
			  const std::vector<std::vector<int> > &contact_indices);
      std::pair<short int, float> change_by(int ichi, double diff,
					    coot::protein_geometry* geom_p);
      std::pair<short int, float> change_by(int imol,
					    int ichi, double diff,
					    const std::vector<std::vector<int> > &contact_indices,
					    coot::protein_geometry *pg,
					    const coot::atom_spec_t &tree_base_atom,
					    short int find_hydrogen_torsions);

      // so that we can highlight a bond, we need to know that atom
      // names of a given bond in the ligand:
      // return empty names on indexing failure
      std::pair<std::string, std::string> atom_names_of_bond(int i) const;

      // return the chi angle (e.g. 1, 2, 3) [not 0-based]
      std::vector<std::pair<int,float> > get_chi_angles() const;
   };
   
} // namespace coot

#endif // CHI_ANGLES_HH
