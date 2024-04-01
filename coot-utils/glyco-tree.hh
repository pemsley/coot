/*
 * coot-utils/glyco-tree.hh
 *
 * Copyright 2024 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <string>

#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"

#include "coot-coord-extras.hh"
#include "tree.hh"

namespace coot {
   class glyco_tree_t {
   public:
      class residue_id_t {
      public:
	 enum prime_arm_flag_t { UNSET, PRIME, NON_PRIME };
	 std::string res_type; // this is tested as empty to see if this object is filled
	 std::string link_type;
	 std::string parent_res_type;
	 residue_spec_t parent_res_spec;
	 unsigned int level;
	 prime_arm_flag_t prime_arm_flag; // are we in the (4') arm?
	 residue_id_t() {}
	 residue_id_t(int level_in,
		      prime_arm_flag_t prime_flag_in,
		      const std::string &res_type_in,
		      const std::string &link_type_in,
		      const std::string &parent_res_type_in,
		      const residue_spec_t &parent_res_spec_in) : res_type(res_type_in),
								  link_type(link_type_in),
								  parent_res_type(parent_res_type_in),
								  parent_res_spec(parent_res_spec_in),
								  level(level_in),
								  prime_arm_flag(prime_flag_in) {}
      };
   private:
      protein_geometry *geom_p;
      std::vector<mmdb::Residue *> linked_residues;
      tree<linked_residue_t> glyco_tree; // 20170503 the constructor now stores its own tree

      bool is_pyranose(mmdb::Residue *r) const;
      tree<linked_residue_t> find_rooted_tree(mmdb::Residue *residue_root_p,
					      const std::vector<mmdb::Residue *> &residues) const;
      tree<linked_residue_t> find_ASN_rooted_tree(mmdb::Residue *residue_p,
						  const std::vector<mmdb::Residue *> &residues) const;
      tree<linked_residue_t> find_stand_alone_tree(const std::vector<mmdb::Residue *> &residues) const;
      void compare_vs_allowed_trees(const tree<linked_residue_t> &tr) const;
      bool compare_trees(const tree<linked_residue_t> &tree_for_testing,
			 const tree<linked_residue_t> &tree_reference) const;
      tree<linked_residue_t> oligomannose_tree() const;
      tree<linked_residue_t>      complex_tree() const;
      tree<linked_residue_t>       hybrid_tree() const;
#ifndef SWIG
      static bool residue_comparitor(mmdb::Residue *res1, mmdb::Residue *res2) {
	 return (residue_spec_t(res1) < residue_spec_t(res2));
      }
#endif
      void print(const tree<linked_residue_t> &glyco_tree) const;
      std::vector<mmdb::Residue *> residues(const tree<linked_residue_t> &glyco_tree) const;
      void output_internal_distances(mmdb::Residue *residue_p,
				     mmdb::Residue *parent_residue_p,
				     double dist_lim,
				     std::ofstream &f) const;
      // old, all-residue
      void output_internal_distances(mmdb::Residue *residue_p,
				     std::vector<mmdb::Residue *> residues,
				     double dist_lim,
				     std::ofstream &f) const;
      residue_id_t::prime_arm_flag_t get_prime(mmdb::Residue *residue_p) const;
      int get_level(mmdb::Residue *residue_p) const;

   public:
      glyco_tree_t(mmdb::Residue *residue_p, mmdb::Manager *mol, protein_geometry *geom_p_in);
      std::vector<mmdb::Residue *> residues(const coot::residue_spec_t &containing_res_spec) const;
      void internal_distances(double dist_lim, const std::string &file_name) const;
      residue_id_t get_id(mmdb::Residue *residue_p) const;
      // for tree comparison
      tree<linked_residue_t> get_glyco_tree() const { return glyco_tree; }
      bool compare_trees(const tree<linked_residue_t> &tree_in) const;
      std::vector<std::pair<coot::residue_spec_t, coot::residue_spec_t> > matched_pairs(const tree<linked_residue_t> &t_in) const;
   };

}