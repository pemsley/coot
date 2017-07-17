
// #include <mmdb2/mmdb_manager.h>

#include "flev-annotations.hh"

namespace coot {

   // return something interesting at some stage
   std::vector<fle_ligand_bond_t>
   protein_ligand_interactions(mmdb::Residue *residue_p,
			       mmdb::Manager *mol,
			       protein_geometry *geom_p,
			       float h_bond_dist_max);

}
