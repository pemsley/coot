
#include <clipper/core/xmap.h>
#include "mini-mol/mini-mol.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

   // Do I want to make this a class?

   minimol::fragment
   multi_build_N_terminal_ALA(mmdb::Residue *res_p,
			      mmdb::Residue *upstream_neighbour_p,
			      const std::string &chain_id,
			      float b_factor_in,
			      int n_trials,
			      const protein_geometry &geom,
			      const clipper::Xmap<float> &xmap_in,
			      std::pair<float, float> mean_and_variance,
			      bool debug_trials);

   minimol::fragment
   multi_build_C_terminal_ALA(mmdb::Residue *res_p,
			      mmdb::Residue *downstream_neighbour_p,
			      const std::string &chain_id,
			      float b_factor_in,
			      int n_trials,
			      const protein_geometry &geom,
			      const clipper::Xmap<float> &xmap_in,
			      std::pair<float, float> mean_and_variance,
			      bool debug_trials);

   minimol::fragment
   multi_build_terminal_ALA(int offset, // direction
			    mmdb::Residue *res_p,
			    mmdb::Residue *up_or_down_stream_neighbour_p, // depends on offset
			    const std::string &chain_id,
			    float b_factor_in,
			    int n_trials,
			    const protein_geometry &geom,
			    const clipper::Xmap<float> &xmap_in,
			    std::pair<float, float> mean_and_variance,
			    bool debug_trials);

   minimol::fragment
   multi_build_terminal_residue_addition_forwards_2018(mmdb::Residue *res_p,
						       mmdb::Residue *upstream_neighbour_p,
						       const std::string &chain_id,
						       float b_factor_in,
						       int n_trials,
						       const coot::protein_geometry &geom,
						       const clipper::Xmap<float> &xmap,
						       std::pair<float, float> mv,
						       bool debug_trials);

   void refine_end(coot::minimol::fragment *many_residues,
		   int seqnum, int offset,
		   const protein_geometry &geom,
		   const clipper::Xmap<float> &xmap_in);

   bool does_residue_fit(const coot::minimol::residue &res, const clipper::Xmap<float> &xmap,
			 std::pair<float, float> mv);

}
