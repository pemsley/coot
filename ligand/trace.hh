/* ligand/trace.hh
 * 
 * Copyright 2015 by Medical Research Council
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


#ifndef LIGAND_TRACE_HH 
#define LIGAND_TRACE_HH

#include <map>
#include <clipper/core/xmap.h>

#include "mini-mol/mini-mol.hh"
#include "geometry/protein-geometry.hh"

#ifndef M_PI
#define M_PI       3.141592653589793238466
#endif

namespace coot {


   // This should be inside trace (it's too generic a name)
   
   // This is a node in the map, i.e. there is a vector of these for each atom_index
   // 
   // For reverse nodes, the atom_idx and the connecting node idx are swapped, so
   // that on construction of a peptide, this will need to be taken into account.
   // 
   class scored_node_t {
   public:
      unsigned int atom_idx;
      double spin_score;
      double alpha; // angle at which spin_score is recorded
      bool reversed_flag;
      bool udd_flag;
      scored_node_t(unsigned int idx_in, double score_in, double angle_in) {
	 atom_idx = idx_in;
	 spin_score = score_in;
	 alpha = angle_in;
	 reversed_flag = false;
	 udd_flag = false;
      }
      scored_node_t(unsigned int idx_in, double score_in, double angle_in,
		    bool reversed_flag_in) {
	 atom_idx = idx_in;
	 spin_score = score_in;
	 alpha = angle_in;
	 reversed_flag = reversed_flag_in;
	 udd_flag = false;
      }
      scored_node_t(const scored_node_t &sn, bool reversed_flag_in) {
	 atom_idx   = sn.atom_idx;
	 spin_score = sn.spin_score;
	 alpha      = sn.alpha;
	 reversed_flag = reversed_flag_in;
	 udd_flag = false;
      }
      scored_node_t() {
	 atom_idx = 999999;
	 spin_score = -9999;
	 alpha = -1;
	 reversed_flag = false;
	 udd_flag = false;
      }
      bool operator==(const scored_node_t &other) const
      { return (other.atom_idx == atom_idx); }
      static bool sort_scores(const scored_node_t &s1, const scored_node_t &s2) {
	 return (s2.spin_score < s1.spin_score);
      }
      static bool sort_pair_scores(const std::pair<unsigned int, scored_node_t> &s1, const std::pair<unsigned int, scored_node_t> &s2) {
	 return (s2.second.spin_score < s1.second.spin_score);
      }
      friend std::ostream& operator<<(std::ostream& s, scored_node_t sn);
   };
   std::ostream& operator<<(std::ostream &s, scored_node_t sn);

   class indexed_frag_t {
   public:
      scored_node_t node_1;
      scored_node_t node_2;
      minimol::fragment f1;
      minimol::fragment f2;
      double score;
      indexed_frag_t(const scored_node_t &n_1,
		     const scored_node_t &n_2,
		     const minimol::fragment &frag_1,
		     const minimol::fragment &frag_2) {
	 node_1 = n_1;
	 node_2 = n_2;
	 f1 = frag_1;
	 f2 = frag_2;
      }
      double get_score() const { return node_1.spin_score + node_2.spin_score; }
   };

   
   class trace {

      clipper::Xmap<float> xmap; // a copy
      float rmsd_cut_off;
      float flood_atom_mask_radius;
      float map_mean;
      float map_rmsd;

      float scale_CO;
      float scale_CO_low;
      float scale_CO_anti;
      float scale_perp;
      float scale_mid;
      float scale_non_line;
      float scale_N;
      float scale_N_accpt;
      
      coot::minimol::molecule get_flood_molecule() const;
      std::vector<std::pair<unsigned int, unsigned int> >
      atom_pairs_within_distance(const minimol::molecule &flood_mol,
				 double trans_dist,
				 double trans_dist_variation);
      // this sets members atom_selection and n_selected_atoms (if mol).
      std::vector<std::pair<unsigned int, unsigned int> >
      atom_pairs_within_distance(mmdb::Manager *mol,
				 double trans_dist,
				 double trans_dist_variation);
      
      std::vector<std::pair<unsigned int, coot::scored_node_t> > 
      spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd);
      
      std::pair<unsigned int, scored_node_t>
	 spin_score(unsigned int idx_1, unsigned int idx_2) const;

      // no longer use sas, because the indexing is strange.  We will
      // use pure mmdb mol and atom selection
      // 
      // minimol::molecule mol_for_sas;
      // std::vector<minimol::atom *> sas;

      mmdb::Manager *mol;
      mmdb::PAtom *atom_selection;
      int n_selected_atoms;
      int selhnd; // for atom_selection

      // presumes atom selection has been set
      clipper::Coord_orth index_to_pos(unsigned int idx) const {
	 return clipper::Coord_orth(atom_selection[idx]->x,
				    atom_selection[idx]->y,
				    atom_selection[idx]->z);
      }
      // presumes atom selection has been set
      std::string index_to_name(unsigned int idx) const {
	 return std::string(atom_selection[idx]->name);
      } 
      
      void trace_graph();
      
      // which atoms are connected to which other atoms
      // (forwards)
      std::map<unsigned int, std::vector<scored_node_t> > fwd_connection_map;
      // back
      std::map<unsigned int, std::vector<scored_node_t> > bck_connection_map;

      // fill the connection_map map with scores
      void make_connection_map(const std::vector<std::pair<unsigned int, scored_node_t> > &scores);

      double frag_score_crit;
      void set_frag_score_crit(const std::vector<std::pair<unsigned int, scored_node_t> > &scores);

      // 2 residues, one of which has two atom (CA, N)
      // (can't be const because it uses the map connection_map)
      minimol::fragment make_fragment(const std::pair<unsigned int, scored_node_t> &scored_node,
				      int res_no_base,
				      std::string chain_id);

      // minimol::residue
      minimol::fragment
      make_residue(const std::pair<unsigned int, scored_node_t> &scored_node_a,
		   const std::pair<unsigned int, scored_node_t> &scored_node_b,
		   int res_no_base,
		   std::string chain_id) const;

      void output_spin_score(const std::pair<unsigned int, scored_node_t> &score,
			     unsigned int atom_idx_1,
			     unsigned int atom_idx_2) const;
      


      enum dir_t { FORWARDS, BACKWARDS };
      
      void
      next_vertex(const std::vector<scored_node_t> &path,
		  unsigned int depth, scored_node_t this_scored_vertex);

      std::vector<scored_node_t>
      get_neighbours_of_vertex_excluding_path(unsigned int this_vertex,
					      const std::vector<scored_node_t> &path,
					      dir_t dir);
      
      void print_tree(const std::vector<unsigned int> &path) const;

      // accumlate interesting trees here
      std::vector<std::vector<scored_node_t> > interesting_trees;
      // 
      void add_tree_maybe(const std::vector<scored_node_t> &path);
      double path_candidate_angle(unsigned int candidate_vertex,
				  const std::vector<scored_node_t> &path) const;
      void print_interesting_trees() const;
      void sort_filter_interesting_trees(); // long trees to the top

      static bool sort_trees_by_length(const std::vector<scored_node_t> &v1,
				       const std::vector<scored_node_t> &v2) {
	 return (v1.size() > v2.size());
      }

      bool using_test_model;

      bool nice_fit(const minimol::residue &r1, const minimol::residue &r2) const;
      bool nice_fit(const minimol::fragment &f1) const;

      double get_fit_score(const minimol::residue &r1, const minimol::residue &r2) const;

      void frag_to_pdb(const minimol::fragment &frag, const std::string &fn) const;
      

      minimol::fragment merge_fragments(const minimol::fragment &f1,
					const minimol::fragment &f2) const;

      // we return an bool to allow a negative value to let the caller know that we didn't
      // find a good next node.
      // 
      std::pair<bool, scored_node_t>
      build_2_choose_1(unsigned int atom_idx, const std::vector<scored_node_t> &start_path,
		       int resno_base,
		       const std::string &chain_id,
		       dir_t dir);

      std::pair<std::vector<scored_node_t>, minimol::fragment>
      follow_fragment (unsigned int atom_idx, const std::vector<scored_node_t> &start_path,
		       int res_no_base,
		       const std::string &chain_id,
		       dir_t dir);
      void add_cbetas(minimol::fragment *frag); // modify frag

      // potentially modify frag_store
      // 
      void add_replace_reject(std::vector<std::pair<std::vector<scored_node_t>, minimol::fragment> > &frag_store, const std::pair<std::vector<scored_node_t>, minimol::fragment> &trial) const;

      void set_inital_density_scales() {
	 scale_CO       =  1.95;
	 scale_CO_low   = -0.6;
	 scale_CO_anti  = -0.15; // weak
	 scale_perp     = -0.8;
	 scale_mid      =  1.6;  // anything more than 1.5
	 scale_non_line =  0.4;
	 scale_N        =  1.2;
	 scale_N_accpt  =  0.0;
      }

      void set_scales(double s[8]) {
	 scale_CO       = s[0];
	 scale_CO_low   = s[1];
	 scale_CO_anti  = s[2];
	 scale_perp     = s[3];
	 scale_mid      = s[4];
	 scale_non_line = s[5];
	 scale_N        = s[6];
	 scale_N_accpt  = s[7];
      }

      double ks_test(const std::vector<std::pair<unsigned int, scored_node_t> > &scores); 

      // Rama terminal addition/refine trace.
      void multi_peptide(const std::vector<std::pair<std::vector<coot::scored_node_t>, coot::minimol::fragment> > &frag_store, const protein_geometry &geom, std::pair<float, float> &mv);

      void move_seeds_close_to_origin(std::vector<coot::minimol::fragment> *seeds) const;

   public:
      trace(const clipper::Xmap<float> &xmap_in);
      void set_atom_mask_radius(float r) { flood_atom_mask_radius = r; }
      void action();

      // testing
      void test_model(mmdb::Manager *mol);

      ~trace() {
	 // crash.  Why?
	 //  if (atom_selection)
	 // mol->DeleteSelection(selhnd);
      }

      double ks_test(); 
      void optimize_weights(mmdb::Manager *mol);
      std::vector<minimol::fragment> make_seeds();
      std::string frag_idx_to_chain_id(unsigned int idx) const;
      
   };

   class trace_path_eraser {
   public:
      std::vector<scored_node_t> lp;
      unsigned int crit_for_match;
      trace_path_eraser(const std::vector<scored_node_t> &long_path,
			unsigned int crit_for_match_in) {
	 lp = long_path;
	 crit_for_match = crit_for_match_in;
      }
      bool operator() (const std::vector<scored_node_t> &interesting_path) const {
	 bool r = false;
	 unsigned int n_match = 0;
	 for (unsigned int i=0; i<interesting_path.size(); i++) {
	    for (unsigned int j=0; j<lp.size(); j++) { 
	       if (lp[j] == interesting_path[i]) {
		  n_match++;
	       }
	    }
	    if (n_match >= crit_for_match) {
	       r = true;
	       break;
	    }
	    if (n_match >= (interesting_path.size() -1)) {
	       r = true;
	       break;
	    }
	 }
	 return r;
      }
   };

   class frag_store_eraser {
   public:
      unsigned int crit_for_match;
      std::pair<std::vector<scored_node_t>, minimol::fragment> trial;      
      const clipper::Xmap<float> *xmap;
      frag_store_eraser(const std::pair<std::vector<scored_node_t>, minimol::fragment> &trial_in,
			const clipper::Xmap<float> *xmap_in,
			unsigned int crit_for_match_in) {
	 xmap = xmap_in;
	 crit_for_match = crit_for_match_in;
	 trial = trial_in;
      }
      bool operator() (const std::pair<std::vector<scored_node_t>, minimol::fragment> &member) const {
	 bool r = false;
	 unsigned int n_match = 0;
	 for (unsigned int i=0; i<member.first.size(); i++) { 
	    for (unsigned int j=0; j<trial.first.size(); j++) { 
	       if (member.first[i] == trial.first[j]) {
		  n_match++;
	       } 
	    }
	 }
	 if (n_match >= crit_for_match)
	    r = true;
	 return r;
      } 

   };
}

#endif // LIGAND_TRACE_HH

