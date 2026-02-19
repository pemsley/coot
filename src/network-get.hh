
#include <string>

//! \brief download and read the model (and data) for the given PDB accession code from PDBe
//!
//! This is an asynchronous function and will return immediately and invoke a download subthread.
//! Check later so see if a molecule has loaded.
//!
//! @param pdb_accession_code is the PDB accession code to fetch
//! @param mode 1 means fetch the "mtz" (to make the map) and
//!             0 means fetch the coordinates model (.pdb or .cif file)
void network_get_accession_code_entity(const std::string &pdb_accession_code, int mode);
