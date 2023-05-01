
#include <string>
#include <vector>

#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"
#include "lidia-core/lbg-shared.hh"
#include "flev-annotations.hh"

//
std::ostream& coot::operator<<(std::ostream &s, coot::fle_ligand_bond_t flb) {

   s << "Ligand-H-bond: " << flb.bond_type << " lig-at: " << flb.ligand_atom_spec
     << " " << flb.interacting_residue_atom_spec << " length: " << flb.bond_length;
   if (flb.is_H_bond_to_water)
      s << " (water)";
   return s;
}
