
#include "mmdb_manager.h"
#ifndef HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR
#include "protein-geometry.hh"

namespace coot {

   class nomenclature {
      CMMDBManager *mol_;
      int test_and_fix_PHE_TYR_nomenclature_errors(CResidue *residue_p); 

   public:
      nomenclature(CMMDBManager *mol) {
	 mol_ = mol;
      }
      // Here we rename atoms to fix nomeclature errors. Note ILEs are not fixed
      // by renaming atoms.
      // 
      std::vector<CResidue *> fix(protein_geometry *Geom_p); // adjust mol as needed.
                                                             // and Geom_p.
   };
}
