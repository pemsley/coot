/*
 * lidia-core/pdbqt-common.hh
 *
 * Shared, RDKit-free pieces for writing PDBQT files: one writable atom, the
 * ATOM/HETATM line formatter, and element-only AutoDock typing. These live in
 * lidia-core (the lower library) so both the receptor writer (coot-utils/pdbqt)
 * and the ligand writer (lidia-core/pdbqt-ligand) can use them without a
 * backwards library dependency.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option) any
 * later version.
 */

#ifndef LIDIA_CORE_PDBQT_COMMON_HH
#define LIDIA_CORE_PDBQT_COMMON_HH

#include <string>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   namespace pdbqt {

      //! One atom to be written to a PDBQT file. The mmdb atom is authoritative for
      //! coordinates and naming; charge and AutoDock type are assigned by the writer.
      class writable_atom_t {
      public:
         int rdkit_idx;        // index into an RDKit mol (-1 if not RDKit-derived)
         mmdb::Atom *at;
         double charge;
         std::string ad_type;
         int serial;           // 1-based, assigned at write time
         writable_atom_t(int idx, mmdb::Atom *a, double q, const std::string &t) :
            rdkit_idx(idx), at(a), charge(q), ad_type(t), serial(-1) {}
      };

      //! Format a single PDBQT ATOM/HETATM record (columns 1-66 are standard PDB;
      //! 71-76 the partial charge, 78-79 the AutoDock type).
      std::string atom_line(const writable_atom_t &wa, bool is_het);

      //! Element-only AutoDock type (fallback when there is no dictionary topology).
      std::string ad_type_from_element(const std::string &element);

   }
}

#endif // LIDIA_CORE_PDBQT_COMMON_HH
