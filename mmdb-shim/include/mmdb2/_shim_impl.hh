// mmdb-shim — architecture B implementation umbrella. Every public mmdb2/*.h
// wrapper resolves here (when COOT_USE_MMDB_SHIM is defined); this file just
// stitches the implementation layers together. See MMDB_SHIM_Recon_and_Plan.md.
//
// The shim presents MMDB's API but is backed by a live gemmi::Structure. The
// mmdb:: hierarchy classes are stable wrapper nodes over that structure (each holds
// Manager* + parent* + a cached sibling index, resolves to gemmi via g(), and
// exposes children as a vector<T*> that is both the identity cache and the
// PPAtom/PPResidue table). Field access is via accessors: numeric ->
// reference-returning x()/…; float occ/b_iso -> value get + set_*; char[] -> pstr
// getter + set_*. gemmi is used at the boundaries (I/O, symmetry, neighbour search)
// and the shim hand-implements only what gemmi lacks (UDData slots, the
// selection-handle engine, the stable-pointer wrapper tree, subgraph matching).
//
// The implementation is split into four dependency-ordered layers (each opens its
// own `namespace mmdb`); this file must include them in order:
//   1. _shim_types.hh     — typedefs, enums, UDStore, forward decls, record classes
//   2. _shim_hierarchy.hh — Atom / Residue / Chain / Model class definitions
//   3. _shim_manager.hh   — Manager (+ nested Selection / UDReg)
//   4. _shim_inline.hh    — g() resolvers, UDData helpers, out-of-line method bodies,
//                           the detail:: selection matchers, and free functions
// then the two sibling subsystems that build on the completed hierarchy.
#pragma once

#include "_shim_types.hh"
#include "_shim_hierarchy.hh"
#include "_shim_manager.hh"
#include "_shim_inline.hh"

// mmdb::mmcif::* (thin veneer over gemmi::cif) — re-opens mmdb{mmcif{...}}.
// pstr/cpstr/realtype are already in scope from the layers above.
#include "_mmcif_impl.hh"

// mmdb::math::{Vertex,Edge,Graph,GraphMatch} — molecular graph + subgraph match.
// Included last so Atom/Residue are complete (MakeGraph).
#include "_graph_impl.hh"
