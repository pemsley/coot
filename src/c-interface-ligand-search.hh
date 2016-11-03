
#ifndef C_INTERFACE_LIGAND_SEARCH_HH
#define C_INTERFACE_LIGAND_SEARCH_HH

#include "ligand/wligand.hh"
#include "c-interface-ligands-widgets.hh"

std::vector<int> execute_ligand_search_internal(coot::wligand *wlig_p);
ligand_wiggly_ligand_data_t ligand_search_install_wiggly_ligands();

#endif // C_INTERFACE_LIGAND_SEARCH_HH

