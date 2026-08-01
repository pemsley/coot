/*
 * lidia-core/pdbqt-ligand.hh
 *
 * Flexible-ligand PDBQT writer: Gasteiger partial charges and a torsion tree
 * (ROOT/BRANCH/ENDBRANCH/TORSDOF) built from acyclic, non-amide rotatable bonds.
 * This needs RDKit, hence it lives in lidia-core (guarded by
 * MAKE_ENHANCED_LIGAND_TOOLS). The shared, RDKit-free helpers are in
 * coot-utils/pdbqt.hh.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#ifndef LIDIA_CORE_PDBQT_LIGAND_HH
#define LIDIA_CORE_PDBQT_LIGAND_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"

namespace coot {

   namespace pdbqt {

      //! Write a flexible-ligand PDBQT file for a single residue. With a dictionary
      //! a torsion tree is written; without one (restraints_p null) the ligand is
      //! written rigidly. comp_id and cid are used only for the REMARK header.
      //!
      //! @return the number of atoms written (0 on failure)
      int write_ligand(mmdb::Residue *residue_p,
                       const dictionary_residue_restraints_t *restraints_p,
                       const std::string &comp_id,
                       const std::string &cid,
                       const std::string &file_name);

   }
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
#endif // LIDIA_CORE_PDBQT_LIGAND_HH
