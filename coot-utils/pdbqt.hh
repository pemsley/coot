/*
 * coot-utils/pdbqt.hh
 *
 * Write PDBQT (AutoDock/Vina input) files. PDBQT is PDB with two extra pieces of
 * per-atom data - a partial charge (columns 71-76) and an AutoDock atom type
 * (columns 78-79) - and, for flexible ligands, a torsion tree.
 *
 * This file holds the shared, RDKit-free machinery: the rigid receptor writer and
 * the per-atom formatting/typing helpers. The flexible-ligand writer (which needs
 * RDKit) lives in lidia-core/pdbqt-ligand.hh.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#ifndef COOT_UTILS_PDBQT_HH
#define COOT_UTILS_PDBQT_HH

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"
#include "lidia-core/pdbqt-common.hh"   // writable_atom_t, atom_line, ad_type_from_element

namespace coot {

   namespace pdbqt {

      //! Write a rigid, whole-molecule (receptor) PDBQT file. Native (no RDKit):
      //! AutoDock types come from the CCP4 energy-library type and hydrogen
      //! positions; charges from the reference table (standard residues) or a
      //! per-energy-type fallback. Non-polar hydrogens are merged; waters dropped.
      //!
      //! @return the number of atoms written (0 on failure)
      int write_receptor(mmdb::Manager *mol, protein_geometry &geom, int imol_enc,
                         const std::string &file_name);

   }
}

#endif // COOT_UTILS_PDBQT_HH
