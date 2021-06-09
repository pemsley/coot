/* ligand/helix-placement.cc
 * 
 * Copyright 2005 The University of York
 * Author: Paul Emsley, Kevin Cowtan
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include "clipper/core/coords.h"
#include "clipper/core/xmap.h"
#include "mini-mol/mini-mol.hh"
#include "coot-utils/coot-map-utils.hh"

namespace coot {

   class scored_helix_info_t {
   public:
      minimol::molecule mol;
      float score;
      scored_helix_info_t(const minimol::molecule &mol_in, float score_in):mol(mol_in) {
	 score = score_in;
      }
      scored_helix_info_t() { score = -9999999.9;}
   };
   bool compare_scored_strands(const scored_helix_info_t &a, const scored_helix_info_t &b);

   class eigen_info_t {
   public:
      clipper::RTop_orth rtop;
      int best_eigen_value_index;
      std::vector<double> eigen_values;
      eigen_info_t(const clipper::RTop_orth &rtop_in,
		   int best_eigen_index_in,
		   const std::vector<double> &eigen_values_in) : rtop(rtop_in), eigen_values(eigen_values_in) {
	 best_eigen_value_index = best_eigen_index_in;
      }
   };

   class helix_placement_info_t {
   public:
      // the first mol [0] is the best fit.  If there is a second one
      // in [1] then it is the reverse fit.
      std::vector<minimol::molecule> mol;
      short int success;
      std::string failure_message;
      helix_placement_info_t(const minimol::molecule &mol_in,
			     short int success_in,
			     const std::string &f_message_in) {
	 mol.clear();
	 mol.push_back(mol_in);
	 success = success_in;
	 failure_message = f_message_in;
      }
      helix_placement_info_t(const std::vector<minimol::molecule> &mol_v,
			     short int success_in,
			     const std::string &f_message_in) :mol(mol_v), failure_message(f_message_in) {
	 success = success_in;
      }
   };
      
   class helix_placement {

      clipper::Xmap<float> xmap;
      clipper::Coord_orth
      move_helix_centre_point_guess(const clipper::Coord_orth &pt, float density_max) const;
      eigen_info_t helix_eigen_system(const clipper::Coord_orth pt, float search_radius) const;
      helix_placement_info_t get_20_residue_helix_standard_orientation(int n_residues,
								       float b_factor) const;
      helix_placement_info_t get_20_residue_helix(int nresidues) const;
      float score_helix_position(const minimol::molecule &m) const;
      // pair(other atoms, c-betas)
      std::vector<clipper::RTop_orth> 
      optimize_rotation_fit(const minimol::molecule &helix,
					       const clipper::RTop_orth &axis_ori,
					       const clipper::Coord_orth &helix_point) const;
      std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> >
      decompose_helix_by_cbeta(minimol::molecule &m) const;
      util::density_stats_info_t score_residue(const minimol::residue &residue) const;
      util::density_stats_info_t score_atoms(const std::vector<clipper::Coord_orth> &atom_pos) const;
      std::pair<int, int> trim_ends(minimol::fragment *m, float min_density_limit) const; // modify m by chopping off residues
      int trim_end(minimol::fragment *m, short int end_type, float min_density_limit) const;

      void build_on_N_end(minimol::fragment *m, float min_density_limit, float b_factor) const;
      void build_on_C_end(minimol::fragment *m, float min_density_limit, float b_factor) const;
      minimol::residue
      build_N_terminal_ALA(const clipper::Coord_orth &prev_n,
			   const clipper::Coord_orth &prev_ca,
			   const clipper::Coord_orth &prev_c,
			   int seqno,
			   float b_factor) const;
      minimol::residue
      build_C_terminal_ALA(const clipper::Coord_orth &prev_n,
			   const clipper::Coord_orth &prev_ca,
			   const clipper::Coord_orth &prev_c,
			   int seqno,
			   float b_factor) const;
      // tinker with m
      void trim_and_grow(minimol::molecule *m, float min_density_limit, float b_factor) const;

      // factoring out for strand tubes
      clipper::RTop_orth
      find_best_tube_orientation(clipper::Coord_orth ptc,
				 double cyl_len, double cyl_rad,
				 float density_level) const; // uses member data xmap

      scored_helix_info_t fit_strand(const minimol::molecule &mol,
				     const clipper::RTop_orth &rtop,
				     int imol,
				     float map_rmsd) const;
      
      std::vector<scored_helix_info_t>
      find_strand_candidates_by_shift_sampling(const minimol::molecule &mol,
					       const clipper::RTop_orth &rtop) const;

      
   public:
      explicit helix_placement(const clipper::Xmap<float> &xmap_in) : xmap(xmap_in) { }
      void discrimination_map() const; //debugging

      // Pass the initial testing helix length.  Typically start with
      // 20, try 12 if that fails.
      helix_placement_info_t place_alpha_helix_near(const clipper::Coord_orth &pt,
						    int n_helix_residues_start,
						    float density_level_for_trim,
						    float b_factor,
						    float map_rmsd) const;
      
      // Kevin's engine: do MR-like search on the surface of a
      // cylinder, not just the eigen vectors
      //
      // For triming and growing we grow into density that is at least
      // density_level_for_trim.
      //
      // For scoring (and to help not to build into "U" sites) we use
      // high_density_turning_point (that is we begin to score
      // progressively badly when the density is higher than
      // high_density_turning_point) i.e. high_density_turning_point
      // is the top of the "triangle".
      // 
      helix_placement_info_t place_alpha_helix_near_kc_version(const clipper::Coord_orth &pt,
							       int n_helix_residues_start,
							       float density_level_for_trim,
							       float high_density_turning_point,
							       float b_factor,
							       float map_rmds) const;

      // and now for strands, we use much of the same code, including
      // the perhaps mis-leading helper class names
      //
      // n_sample_strands is the number of sample strands that we should get from the database
      // 
      helix_placement_info_t place_strand(const clipper::Coord_orth &pt, int n_residues,
					  int n_sample_strands, float sigma_level);
      
   };
}
