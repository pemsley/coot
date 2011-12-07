
#include <string>
#include <vector>

#include <mmdb/mmdb_manager.h>
#include "protein-geometry.hh"
#include "lbg-shared.hh"
#include "coot-coord-utils.hh"

#include "flev-annotations.hh"

// 
std::ostream& coot::operator<<(std::ostream &s, coot::fle_ligand_bond_t flb) {

   s << "Ligand-H-bond: " << flb.bond_type << " lig-at: " << flb.ligand_atom_spec
     << " " << flb.interacting_residue_atom_spec << " length: " << flb.bond_length;
   
   return s;
}
