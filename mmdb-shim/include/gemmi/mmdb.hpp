// -*- mode: c++; -*-
//
// mmdb-shim: replacement for gemmi's own <gemmi/mmdb.hpp> bridge.
//
// The real gemmi/mmdb.hpp converts gemmi::Structure <-> mmdb::Manager by copying
// field-by-field through real-mmdb PUBLIC FIELDS (atom->x = ..., chain->seqRes,
// mmdb::newAtom(), Cryst::Z, ...). That is fundamentally incompatible with the
// accessor-based shim (x/y/z are methods, not fields).
//
// But it is also unnecessary: the shim's mmdb::Manager already OWNS a live
// gemmi::Structure (Manager::st), so the conversions are near-trivial. This header
// is picked up instead of gemmi's because the shim include dir precedes gemmi's on
// the -I path. Only the functions Coot actually calls are provided.
//
// Copyright 2026 by Medical Research Council Laboratory of Molecular Biology
#ifndef COOT_MMDB_SHIM_GEMMI_MMDB_HPP
#define COOT_MMDB_SHIM_GEMMI_MMDB_HPP

#include <mmdb2/mmdb_manager.h>   // shim mmdb::Manager (wraps gemmi::Structure)
#include <gemmi/model.hpp>        // gemmi::Structure

namespace gemmi {

// gemmi::Structure -> mmdb::Manager: the shim Manager IS a gemmi::Structure
// wrapper, so adopt the structure and (re)build the wrapper hierarchy.
inline void copy_to_mmdb(const Structure& st, mmdb::Manager* manager) {
  if (!manager) return;
  manager->st = st;
  manager->build_from_gemmi();
}

// mmdb::Manager -> gemmi::Structure: st is the source of truth (all wrappers
// resolve into it), so just return it.
inline Structure copy_from_mmdb(mmdb::Manager* manager) {
  return manager ? manager->st : Structure();
}

// SEQRES already travels inside Manager::st, so this is a no-op in the shim.
inline void transfer_seqres_to_mmdb(const Structure&, mmdb::Manager*) {}

}  // namespace gemmi

#endif  // COOT_MMDB_SHIM_GEMMI_MMDB_HPP
