
#include "coot-utils/coot-coord-utils.hh"

namespace coot {
   
   class parallel_planes_t {
   public:
      std::vector<atom_spec_t> plane_1_atoms;
      std::vector<atom_spec_t> plane_2_atoms;
      parallel_planes_t(const std::vector<atom_spec_t> &ap1,
			const std::vector<atom_spec_t> &ap2) {
	 plane_1_atoms = ap1;
	 plane_2_atoms = ap2;
      } 
   };
}
