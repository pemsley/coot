/* db-main/db-main.hh
 * 
 * Copyright 2002, 2003, 2006 The University of York
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


#include <vector>
#include <string>

// clipper coords
#include "clipper/core/coords.h"
#include "mini-mol.hh"

// for directory reading
#include <sys/types.h>   
#include <dirent.h>     
// #include <unistd.h> not needed?

// #include <string.h> // for string matching // needed?

namespace coot { 

   bool matches_pdb_name(std::string file_str);
   
   class coord_fragment_wrapper_t {
      
   public:
      std::vector<clipper::Coord_orth> fragment; 
      std::string segment_id; 
      int size() { return fragment.size(); }
   };

   class main_fragment_t {
      int i_start_res_;
   public:
      int ilength; 
      int molecule_number; 
      std::vector<float> sqrt_eigen_values;
      float eigen_similarity_score; // how similar is this to target eigens?
      std::string segment_id;
      std::pair<bool,clipper::Coord_orth> middle_carbonyl_oxygen_position;
      main_fragment_t(int i_start_res, int mol_no,
		      std::vector<float> eigns_sqrt, std::string seg_id,
		      std::pair<bool, clipper::Coord_orth> mid_o_pt,
		      int ilen); 
      int i_start_res() const {return i_start_res_;}
   };

   class db_fitting_result {
   public:
      clipper::RTop_orth rtop;
      int db_frag_index;
      float deviance;
      int ilength;
      int istart_res_of_ca_target;
      db_fitting_result(const clipper::RTop_orth rtop_in, int db_frag_i, float devi, int ist_res, int ilen) { 
	 rtop = rtop_in; 
	 db_frag_index = db_frag_i; 
	 deviance = devi;
	 ilength = ilen;
	 istart_res_of_ca_target = ist_res; 
      }
   };

   class peptide_match_fragment_info_t {
   public:
      minimol::fragment fragment;
      float devi;
      peptide_match_fragment_info_t(const minimol::fragment &frag_in, float devi_in) {
	 fragment = frag_in;
	 devi = devi_in;
      } 
   }; 

   class weighted_residue : public minimol::residue {
   public:
      float weight_sum;
      float weight_sum_cb;
      short int have_cb_flag;
      int cb_index;

      weighted_residue() {
	 have_cb_flag = 0;
 	 weight_sum = 0.0;
 	 weight_sum_cb = 0.0;
      }

      // Adding residues with the weight_target_devi.
      // 
      // Given a residue with the weight_target_devi.  We apply
      // weight_pos_in_frag... 
      //
      // Note that we implicitly presume that all in_res has CA, C, N and O,
      // but only some have CB.  So we need to keep a separate track of the
      // sum of the CB weights (which is different from the sum of all
      // weights).
      //
      // Notice that when we start, the vector atoms is empty, so we do a
      // check that only fails for the first time, where we do an addatom,
      // rather than adding to the position of already esisting atoms.
      //
      // For the first time, we check that CA, C, N, O exist after we have
      // added in_res.
      // 
      void add_residue_pos(const minimol::residue &in_res,
			   const clipper::RTop_orth &rtop,
			   float weight);

      minimol::residue pull_residue() const; 
   };

   class db_main {
      std::vector<coot::main_fragment_t> mainchain_frag_db;
      std::vector<coot::minimol::molecule> molecule_list; 

      // For each set of 6 residues (Ca's thereof), we fit a set of
      // "similar by eigen" structures.  This is a vector of the sets
      // of (typically) 6 atoms (typically stepped by 3).
      //
      // We store also the starting residue number (of the target) of
      // the match in the db_fitting_result.
      //
      // Match target fragment fills this vector.
      // 
      std::vector <std::vector <db_fitting_result> > big_results;
      float max_devi;
      std::vector<weighted_residue> output_fragment;
      // Save these in match_target_fragment(), because they will be
      // needed in merge_fragments().
      int iresno_start;
      int iresno_end;
      std::string target_fragment_fragment_id; 

      // Given a molecule and a starting residue number, return a
      // vector of length ilength that contains the Ca positions of
      // that fragment.  We use (only) the first fragment (0) in
      // target.
      // 
      // You need to check the length of the returned vector to see if
      // this function was successful.
      // 
      std::vector<clipper::Coord_orth>
      get_target_ca_coords(int iresno, int ilength,
			   const coot::minimol::molecule &target) const;
      std::vector<db_fitting_result>
      fit_reference_structures(float max_devi,
			       std::vector<clipper::Coord_orth> target_ca,
			       int iresno_start,
			       int ilength) const; 

      short int similar_eigens(float tolerance_frac, 
			       const std::vector<float> &target, 
			       const std::vector<float> &frag) const; 

      std::vector<std::string> get_reference_pdb_list() const;
      
      coot::minimol::fragment get_fragment_ca_atoms(int istart, int ilength,
					   const coot::minimol::molecule &m) const;
      clipper::Matrix<float> make_cov_matrix(const std::vector<clipper::Coord_orth> &c) const;
      std::vector<clipper::Coord_orth>
      frag_to_coords(const coot::minimol::fragment &fragment) const;

      // This is an awkward function. We have to generate the fragment
      // from the molecule and then select the (ilength) Ca from the
      // right place (and segment id).
      // 
      std::vector<clipper::Coord_orth>
      mainchain_ca_coords_of_db_frag(int i, int ilength) const;
      float deviance(const std::vector<clipper::Coord_orth> &frag1, 
		     const std::vector<clipper::Coord_orth> &frag2, 
		     const clipper::RTop_orth &rtop) const; 

      // helper function for merge_fragments:
      float weight_pos_in_frag(int ipos, int ilength) const;
      // 
      minimol::residue pull_db_residue(const coot::db_fitting_result &res,
					     int ipos) const;

      minimol::fragment pull_db_fragment(const coot::main_fragment_t &res, int ilength);


      // Pepflip extras:
      std::pair<bool, clipper::Coord_orth> get_middle_ox_pos(const coot::minimol::fragment &fragment) const;
      // a vector of fragments that have been matches to a target 5 CA
      // peptide (that target 5 CA peptide is constructed, for the
      // moment, outside this class.
      std::vector<peptide_match_fragment_info_t> pepflip_fragments;
      void sort_mainchain_fragments_by_eigens(std::vector<float> target_eigens);
      static bool mainchain_fragment_sorter(const coot::main_fragment_t &fit_a,
					    const coot::main_fragment_t &fit_b);
      std::vector<float> tmp_target_eigens;
      void assign_eigen_similarity_scores(const std::vector<float> &target_eigens);

   public:
      db_main() { max_devi = 2.0; }

      // Fill the molecule_list with (minimol) molecules and
      // mainchain_frag_db with fragment info.
      // 
      int fill_with_fragments(int ilength);

      // So, we have a molecule, we want to chop it into fragments,
      // starting at residue istart.  We don't want to keep the actual
      // fragments, we will just store just the eigen values of the
      // fragment that we pull out (and the discard), along with
      // enough other information to get the fragment back if/when we
      // want to use it in recombination.
      //
      // 
      // Fill big results:
      // 
      void match_target_fragment(const coot::minimol::molecule &target_cas,
				 int iresno_target_start, // typically 1
				 int iresno_target_end, 
				 int ilength);

      // Molecular Graphics usage: We only want to load the graphics
      // once.  So let's do that and store it in a static.
      // 
      bool is_empty() const; 
      bool is_empty_of_pepflips() const; 

      // This is the big complex one:
      //
      void merge_fragments();

      // And we return the results in a minimol::fragment:
      // 
      minimol::fragment mainchain_fragment() const;

      // clear results
      void clear_results();

      // ----------------------- Pepflip extras --------------------
      void match_targets_for_pepflip(const minimol::fragment &target_ca_coords_5_res_frag);

      // Check that internal vector against the oxygen position.
      //
      // Return the fraction fo the best fitting peptides that have their
      // oxygen closer than d_crit
      float mid_oxt_outliers(const clipper::Coord_orth &my_peptide_oxt_pos, int resno_oxt, float cutoff);
      static bool pepflip_sorter(const coot::peptide_match_fragment_info_t &fit_a,
				 const coot::peptide_match_fragment_info_t &fit_b);

   };

} // namespace coot
