
#include <mmdb2/mmdb_manager.h>
#ifndef HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR
#include "geometry/protein-geometry.hh"

namespace coot {

   class nomenclature {
      mmdb::Manager *mol_;
      // return swap status.
      int test_and_fix_PHE_TYR_nomenclature_errors(mmdb::Residue *residue_p,
						   bool apply_swap_when_found);

      // 20110725 new function - don't use rotamers, simply use torsions -90 < tor < 90.
      // 
      int test_and_fix_ASP_GLU_nomenclature_errors(mmdb::Residue *residue_p,
						   bool apply_swap_when_found); 

      std::vector<mmdb::Residue *>  fix_and_swap_maybe(protein_geometry *Geom_p,
						  bool apply_swaps); // adjust mol as needed
                                                                     // if flagged.

   public:
      nomenclature(mmdb::Manager *mol) {
	 mol_ = mol;
      }
      // Here we rename atoms to fix nomeclature errors. Note ILEs are not fixed
      // by renaming atoms.
      // 
      std::vector<mmdb::Residue *>  fix(protein_geometry *Geom_p); // adjust mol as needed.

      std::vector<mmdb::Residue *> list(protein_geometry *Geom_p); // Don't touch the molecule,
                                                              // just analyse it.
   }; 
}
