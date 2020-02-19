
#include <mmdb2/mmdb_manager.h>
#include <vector>
#include <string>

namespace coot {

   // select the residues in a chain and call this for every chain
   // in the molecule.
   // return a vector of helical residues in that chain
   //
   std::vector<mmdb::Residue *> like_a_helix(mmdb::Manager *mol, int selection_handle);

   class helical_results_t {
   public:
      bool is_alpha_helix_like;
      bool is_pi_helix_like;
      bool is_3_10_helix_like;
      float sum_delta;
      helical_results_t() {
	 is_alpha_helix_like = false;
	 is_pi_helix_like    = false;
	 is_3_10_helix_like  = false;
	 sum_delta = 0;
      }
   };

   // test_helical_residues should be in the right order. 
   helical_results_t compare_to_helix(const std::vector<mmdb::Residue *> &test_helical_residues,
                                      const std::vector<clipper::Coord_orth> &alpha_helix_ref_positions);

   helical_results_t compare_to_helix(const std::vector<mmdb::Residue *> &helical_residues);

   std::vector<clipper::Coord_orth> alpha_helical_reference_positions();
}
