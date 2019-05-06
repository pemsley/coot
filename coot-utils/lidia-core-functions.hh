
#ifndef LIDIA_CORE_FUNCTIONS_HH
#define LIDIA_CORE_FUNCTIONS_HH

#include "atom-selection-container.hh"
// #include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"

namespace coot {
  // mdl mol file support
  atom_selection_container_t mdl_mol_to_asc(const lig_build::molfile_molecule_t &m);
  atom_selection_container_t mdl_mol_to_asc(const lig_build::molfile_molecule_t &m, float b_factor);
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
  atom_selection_container_t mol_to_asc_rdkit(const std::string &file_name);
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}


#endif // LIDIA_CORE_FUNCTIONS_HH

