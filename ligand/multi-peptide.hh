
#include <clipper/core/xmap.h>
#include "mini-mol/mini-mol.hh"
#include "geometry/protein-geometry.hh"

namespace coot {

   // Do I want to make this a class?

   minimol::fragment
   multi_build_N_terminal_ALA(mmdb::Residue *res_p,
			      const std::string &chain_id,
			      float b_factor_in,
			      int n_trials,
			      const protein_geometry &geom,
			      const clipper::Xmap<float> &xmap_in,
			      std::pair<float, float> mean_and_variance);

   minimol::fragment
   multi_build_C_terminal_ALA(mmdb::Residue *res_p,
			      const std::string &chain_id,
			      float b_factor_in,
			      int n_trials,
			      const protein_geometry &geom,
			      const clipper::Xmap<float> &xmap_in,
			      std::pair<float, float> mean_and_variance);

   minimol::fragment
   multi_build_terminal_ALA(int offset, // direction
			    mmdb::Residue *res_p,
			    const std::string &chain_id,
			    float b_factor_in,
			    int n_trials,
			    const protein_geometry &geom,
			    const clipper::Xmap<float> &xmap_in,
			    std::pair<float, float> mean_and_variance);

   void refine_end(coot::minimol::fragment *many_residues,
		   int seqnum, int offset,
		   const protein_geometry &geom,
		   const clipper::Xmap<float> &xmap_in);
   

   bool does_residue_fit(const coot::minimol::residue &res, const clipper::Xmap<float> &xmap,
			 std::pair<float, float> mv);

}
