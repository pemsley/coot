
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
      void spin_score_pairs(const std::vector<std::pair<unsigned int, unsigned int> > &apwd) const;
      double spin_score(unsigned int idx_1, unsigned int idx_2) const;
      std::vector<minimol::atom *> sas;
   public:
      trace(const clipper::Xmap<float> &xmap_in);
      void set_atom_mask_radius(float r) { flood_atom_mask_radius = r; }
      void action();

   };
}

