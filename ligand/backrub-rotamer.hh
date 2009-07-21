
#ifndef BACKRUB_ROTAMER_HH
#define BACKRUB_ROTAMER_HH

#include <string>

#include <clipper/core/xmap.h>
#include "mini-mol.hh"

namespace coot {

   class backrub {
      CResidue *orig_this_residue;
      CResidue *orig_prev_residue;
      CResidue *orig_next_residue;
      std::string alt_conf;
      clipper::Coord_orth ca_prev;
      clipper::Coord_orth ca_next;
      // thow an exception on failure to find CA of prev or next.
      void setup_prev_next_ca_positions();
      minimol::fragment make_test_fragment(CResidue *r, double rotation_angle) const;
      std::string chain_id;
      float score_fragment(minimol::fragment &frag) const;
      clipper::Xmap<float> xmap;
      CMMDBManager *stored_mol;
      minimol::residue
      make_residue_include_only(CResidue *orig_prev_residue,
				const std::vector<std::string> &prev_res_atoms) const;
      
   public:

      // Throw an exception on failure to construct the backrub internals.
      // Throw an exception if this, previous or next residues are null
      backrub(const std::string &chain_id_in,
	      CResidue *this_r,
	      CResidue *prev_r,
	      CResidue *next_r,
	      const std::string &alt_conf_in,
	      CMMDBManager *mol_in,
	      const clipper::Xmap<float> &xmap_in) {
	 orig_this_residue = this_r;
	 orig_prev_residue = prev_r;
	 orig_next_residue = next_r;
	 setup_prev_next_ca_positions();
	 chain_id = chain_id_in;
	 xmap = xmap_in;
	 alt_conf = alt_conf_in;
	 stored_mol = mol_in;
      }

      // throw an exception on failure to get a good search result.
      std::pair<coot::minimol::molecule, float> search();
   };

}


#endif // BACKRUB_ROTAMER_HH
