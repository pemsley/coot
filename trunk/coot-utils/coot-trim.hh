
#include "mmdb-extras.h"

namespace coot {

   namespace util { 

      enum {TRIM_BY_MAP_DELETE, TRIM_BY_MAP_ZERO_OCC};

      // return the number of trimmed atoms
      int
      trim_molecule_by_map(atom_selection_container_t asc,
			   const clipper::Xmap<float> &xmap,
			   float map_level,
			   short int remove_or_zero_occ_flag,
			   short int waters_only_flag);

   }

}   
