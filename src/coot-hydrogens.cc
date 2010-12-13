
#ifdef MAKE_ENTERPRISE_TOOLS
#include "rdkit-interface.hh"
#endif
 
#include <iostream> // fixes undefined strchr, strchrr problems
#include <mmdb_manager.h>

#include "protein-geometry.hh"
#include "coot-hydrogens.hh"


bool
coot::add_hydrogens(CResidue *residue_p,
		    const coot::dictionary_residue_restraints_t &restraints) {

#ifdef MAKE_ENTERPRISE_TOOLS
   return add_hydrogens_with_rdkit(residue_p, restraints);
#else    
   // return add_hydrogens_with_ccp4_tools(residue_p, restraints);
   return 0;
#endif   

}


