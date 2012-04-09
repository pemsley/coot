/* geometry/coot-atom-graph.cc
 * 
 * Copyright 2004  The University of York
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


#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#include "mini-mol.hh"
#include "coot-node-info.hh"
#include "sequence-assignment.hh"

namespace coot { 
   
   // a helper function for atom_graph:
   // 
   // Each node has a vector of these things (which get added in atom
   // assignment) as it may be possible that a single node has more
   // than one possible atom assignment.
   // 
   class graph_atom_info {
   public:
      int chain_number;
      int residue_number;
      short int is_water_flag;
      minimol::atom atom;
      float weight;
      std::vector <clipper::RTop_orth> rtops;
      graph_atom_info() {}
      graph_atom_info(int residue_number_in, int chain_number_in,
		      const minimol::atom &at_in,
		      std::vector <clipper::RTop_orth> rtops_in,
		      short int is_water_flag_in,
		      float weight_in) {
	 residue_number = residue_number_in;
	 chain_number = chain_number_in;
	 atom = at_in;
	 weight = weight_in;
	 is_water_flag = is_water_flag_in;
	 rtops = rtops_in;
      }
   };

   class chain_helper_info {
   public:
      int direction_score;
      int forward_score;
      int backward_score;
      int n_peptides;
      int chain_number;
      chain_helper_info() {}
      chain_helper_info(int direction_score_in,
			int forward_score_in,
			int backward_score_in,
			int n_peptides_in) {
	 direction_score = direction_score_in;
	 forward_score = forward_score_in;
	 backward_score = backward_score_in;
	 n_peptides = n_peptides_in;
      }
      // Return +1 for positive, -1 for negative and 0 for undetermined.
      // 
      int get_direction() const;
   };

   class atom_graph { 
      
      CMMDBManager *mol_internal_ptr_copy;
      std::vector<std::vector<node_info> > nodes;
      std::vector<coot::minimol::atom> atoms;
      std::vector<std::vector<graph_atom_info> > atom_info;
      std::vector<int> connectedness; // variable for recursion
      void sort_tips(std::vector <std::pair<int, int> > *tips) const;
      std::vector<std::vector<coot::node_info> > get_trace(int i_start_node) const;
      void assign_c_betas();
      void assign_waters();
      short int is_c_alpha_p(int index) const;

      void digraph_trace_along(int inode, 
			       int connection_number,
			       std::vector<std::vector<coot::node_info> > *trace,
			       std::vector<int> *con_local) const;
      chain_helper_info peptide_search(const std::vector<std::vector<coot::node_info> > &t,
			  int i_node_start,
			  const coot::chain_helper_info &hi_in,
			  short int make_assignments_flag);
      double squared(double a) const { return a*a;}

      clipper::Coord_orth tf (const coot::node_info &ni) const {
	 if (ni.symm_trans_needed_flag == 0) { 
	    return atoms[ni.index].pos;
	 }
      }

//       short int is_possible_ca_c_n_c(const std::vector<std::vector<coot::node_info> > &t,
// 				     int t_index_1, int t_index_2, int t_index_3, int t_index_4) const;
      short int is_possible_ca_c_n_c(int index_start,
				     const coot::node_info &t_index_2, 
				     const coot::node_info &t_index_3, 
				     const coot::node_info &t_index_4) const;

      clipper::Coord_orth get_transformed_atom(const clipper::Coord_orth &a,
					       const std::vector<clipper::RTop_orth> &transformations) const;

      // fill the atom_info vector
      //
      std::vector<clipper::RTop_orth>
      make_assignments(int direction,
		       double weight,
		       int res_no,
		       int chain_number,
		       int i_node_peptide_start,
		       const std::vector<clipper::RTop_orth> &running_ops,
		       const coot::node_info &atom_1,
		       const coot::node_info &atom_2,
		       const coot::node_info &atom_3);
      
      void write_molecule_from_atom_info(const std::string &file_name) const;
      std::string chain_id(int chain_number) const;

      std::vector<realtype> cell;

      // this will do for a space group for now:
      std::string spgr;

      // Takes into account potentially many symmetry operations
      // 
      double peptide_distortion_score(short int direction,
				      int index_start,
				      const coot::node_info &t_index_2, 
				      const coot::node_info &t_index_3, 
				      const coot::node_info &t_index_4, 
				      const coot::node_info &t_index_5, 
				      const coot::node_info &t_index_6, 
				      const coot::node_info &t_index_7) const;

      // Doesn't take into account symmetry operations
      // 
      double peptide_distortion_score(short int direction,
				      const clipper::Coord_orth &a1,
				      const clipper::Coord_orth &a2,
				      const clipper::Coord_orth &a3,
				      const clipper::Coord_orth &a4,
				      const clipper::Coord_orth &a5,
				      const clipper::Coord_orth &a6,
				      const clipper::Coord_orth &a7) const;

      // idx is coot::sequence_assignment::side_chain_name_index
      double distortion_score_side_chain(int idx,
					 const std::vector <clipper::Coord_orth> &c) const;
      // sc_nodes: a vector of sidechain nodes, starting at the C-beta.
      double distortion_score_side_chain_old(int i_node_ca,
					     const std::string &residue_type,
					     const std::vector<coot::node_info> &sc_nodes) const;

      // i_node_ca indexes the atom_info vector: 
      void score_all_side_chain_types(int i_node_ca,
				      const std::vector<coot::node_info> &sc_nodes,
				      coot::sequence_assignment::side_chain_score_t *scs) const;

      // The atom order of the vectors must match - and are sidechain
      // specific, obviously.
      // 
      double distortion_score_ser(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_ala(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_arg(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_asp(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_asn(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_gly(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_glu(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_gln(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_phe(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_tyr(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_trp(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_val(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_ile(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_leu(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_pro(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_cys(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_met(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_his(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_lys(const std::vector <clipper::Coord_orth> &c) const;
      double distortion_score_thr(const std::vector <clipper::Coord_orth> &c) const;

      // centre (chiral) atom first:
      double abs_chiral_vol(const clipper::Coord_orth &a1,
			    const clipper::Coord_orth &a2,
			    const clipper::Coord_orth &a3,
			    const clipper::Coord_orth &a4) const;

      void sidechains_search();
      std::vector<coot::node_info>
      get_side_chain_nodes(int i, coot::sequence_assignment::side_chain_name_index) const;
      std::vector<coot::node_info> get_side_chain_nodes(int i) const;

   public:

      // We need access to the molecule's cryst - so we pass the
      // CMMDBManager.  Don't bother to try and const it - it's just
      // too painful.  The danger is of course now that atom_graph is
      // used when mol has gone out of scope/deleted.  Don't do that,
      // then.  Yeah, I know - sorry (not really my fault - there's no
      // way to pass just the crystal-info.).
      // 
      atom_graph(CMMDBManager *mol,
		 const std::vector<std::vector<coot::node_info> > &connection_indices,
		 const std::vector<clipper::Coord_orth> &coords);
      void sort();
      void trace_along(int ic, int connection_number);
      static short int tip_compare(const std::pair<int, int> &a, 
				   const std::pair<int, int> &b);

      // Return a (traced) molecule.
      // 
      minimol::molecule traced() const;
   }; 
} 
