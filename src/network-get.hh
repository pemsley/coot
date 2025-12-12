
#include <string>

//! download and read the model (and data) for the given PDB accession code from PDBe.
//!
//! @param pdb_accession_code is the PDB accession   code
//! @param mode 1 means "mtz" and 0 means the coordinates model (.pdb or .cif file)
void network_get_accession_code_entity(const std::string &pdb_accession_code, int mode);
