
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif
 
#include <iostream> // fixes undefined strchr, strchrr problems
#include <mmdb/mmdb_manager.h>

#include "geometry/protein-geometry.hh"
#include "coot-hydrogens.hh"


std::pair<bool, std::string>
coot::add_hydrogens(CResidue *residue_p,
		    const coot::dictionary_residue_restraints_t &restraints) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   return add_hydrogens_with_rdkit(residue_p, restraints);
#else    
   // return add_hydrogens_with_ccp4_tools(residue_p, restraints);
   return std::pair<bool, std::string> (0, "not implemented");
#endif   

}


