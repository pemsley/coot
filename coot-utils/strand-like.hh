
#ifndef STRAND_LIKE_HH
#define STRAND_LIKE_HH

#include <mmdb2/mmdb_manager.h>
#include <vector>
#include <string>

namespace coot {

   // select the residues in a chain and call this for every chain
   // in the molecule.
   void like_a_strand(mmdb::Manager *mol, int selection_handle);

   class strand_results_t {
   public:
      bool is_beta_strand_like;
      bool is_beta_bulge_like;
      float sum_delta;
      helical_results_t() {
	 is_beta_bulge_like  = false;
	 is_beta_strand_like = false;
	 sum_delta = 0;
      }
   };

   // the atoms need to be in the correct order
   strand_results_t compare_to_strand(const std::vector<mmdb::Residue *> &test_residues,
                                      const std::vector<clipper::Coord_orth> &beta_strand_ref_positions);

   strand_results_t compare_to_strand(const std::vector<mmdb::Residue *> &strand_residues);

   std::vector<clipper::Coord_orth> beta_strand_reference_positions();
}

#endif
