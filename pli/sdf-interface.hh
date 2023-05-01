/* src/sdf-interface.hh
 * 
 * Copyright 2012 by The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <string>
#include <mmdb2/mmdb_manager.h>

#include "geometry/protein-geometry.hh"

/*! \file
  \brief Coot Scripting Interface - SDF/Molfile interface
*/

// 20140226: now we change things so that the interface functions always get generated,
//           and what the function does depends on MAKE_ENHANCED_LIGAND_TOOLS
// 
// #ifdef MAKE_ENHANCED_LIGAND_TOOLS
// 
// was it OK or not?  (i.e. did we not catch an exception)
bool residue_to_sdf_file(int imol, mmdb::Residue *residue_p,  const char *sdf_file_name,
                         const coot::protein_geometry &geom, bool kekulize = true);
bool residue_to_mdl_file_for_mogul(int imol, mmdb::Residue *residue_p, const std::string &mdl_file_name,
                                   const coot::protein_geometry &geom);

// rdkit chemical features.
// bool show_feats(int imol, const char *chain_id, int resno, const char *ins_code); (make generic objects)
bool show_feats(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom);

// This is not an sdf function, perhaps it should be somewhere else - but
// if in rdkit-interface.hh, it would be a singleton.
//

//! \brief
//! import a molecule from a smiles string
//!
//! RDKit is used to interpret the SMILES string
//!
//! no dictionary is generated
//!
//! @return a molecule number or -1 on failure
int import_rdkit_mol_from_smiles(const std::string &smiles, const std::string &comp_id);

// #endif

