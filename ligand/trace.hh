

#ifndef LIGAND_TRACE_HH 
#define LIGAND_TRACE_HH

#include <map>
#include <clipper/core/xmap.h>

#include "../mini-mol/mini-mol.hh"


namespace coot {

   class trace_node {
   public:
      unsigned int atom_idx;
      // the score of the link to the parent
      double score_1;
      double score_2;
      trace_node() {}
      trace_node(unsigned int idx_in, double score_1_in, double score_2_in) {
	 atom_idx = idx_in;
	 score_1 = score_1_in;
	 score_2 = score_2_in;
      }
   };
   
   class trace {

      clipper::Xmap<float> xmap; // a copy
      float rmsd_cut_off;
      float flood_atom_mask_radius;
      coot::minimol::molecule get_flood_molecule() const;
      std::vector<std::pair<unsigned int, unsigned int> >
      atoms_pairs_within_distance(const minimol::molecule &flood_mol,
				  double trans_dist,
				  double trans_dist_variation);
      std::vector<std::pair<unsigned int, unsigned int> >
      atoms_pairs_within_distance(mmdb::Manager *mol,
				  double trans_dist,
				  double trans_dist_variation);
      
      // fill the tr map with scores
      void spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd);
      double spin_score(unsigned int idx_1, unsigned int idx_2) const;
      minimol::molecule mol_for_sas;
      std::vector<minimol::atom *> sas;
      void trace_graph();
      std::map<unsigned int, std::vector<unsigned int> > tr; // which atoms are connected to which other atoms
                                                             // backwards and forwards
      
      std::vector<unsigned int>
      next_vertex(const std::vector<unsigned int> &path,
		  unsigned int depth, unsigned int this_vertex);

      std::vector<unsigned int> get_neighbours_of_vertex_excluding_path(unsigned int this_vertex,
									const std::vector<unsigned int> &path);
      void print_tree(const std::vector<unsigned int> &path) const;

      // accumlate interesting trees here
      std::vector<std::vector<unsigned int> > interesting_trees;
      // 
      void add_tree_maybe(const std::vector<unsigned int> &path);
      double path_candidate_angle(const std::vector<unsigned int> &path,
				  unsigned int candidate_vertex) const;
      void print_interesting_trees() const;
      


   public:
      trace(const clipper::Xmap<float> &xmap_in);
      void set_atom_mask_radius(float r) { flood_atom_mask_radius = r; }
      void action();

      // testing
      void test_model(mmdb::Manager *mol);

   };
}

#endif // LIGAND_TRACE_HH

