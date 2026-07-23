// -*- mode: c++; -*-
//
// mmdb-shim: no-op replacement for <ssm/ssm_align.h>.
//
// The real SSM library (libssm) is a precompiled binary built against REAL MMDB
// 2.0.22, and its headers derive from mmdb::io::Stream / use mmdb binary
// serialization (DefineClass, io::RFile) that the gemmi-backed shim does not
// provide. Worse, even if the headers parsed, passing shim mmdb objects (thin
// gemmi wrappers) into precompiled ssm would be an ABI mismatch -> crash.
//
// Per the project decision to NOT build SSM code against the shim (SSM is being
// removed/replaced), this header shadows the real ssm/ssm_align.h (the shim
// include dir precedes ssm's on the -I path) with a self-contained, header-only,
// no-op ssm::Align. Everything compiles and links (no real libssm symbols are
// referenced); SSM structural superposition simply reports "no match" (identity
// transform), so callers take their no-alignment path.
//
// Copyright 2026 by Medical Research Council Laboratory of Molecular Biology
#ifndef COOT_MMDB_SHIM_SSM_ALIGN_H
#define COOT_MMDB_SHIM_SSM_ALIGN_H

#include <mmdb2/mmdb_manager.h>   // mmdb::mat44/realtype/ivector/rvector/PManager (shim)

namespace ssm {

// return codes (ssm_defs.h RETURN_CODE)
enum RETURN_CODE { RC_Ok, RC_NoHits, RC_NoSuperposition, RC_NoGraph,
                   RC_NoVertices, RC_NoGraph2, RC_NoVertices2, RC_TooFewMatches };
// precision levels (ssm_defs.h PRECISION)
enum PRECISION   { PREC_Highest, PREC_High, PREC_Normal, PREC_Low, PREC_Lowest };
// connectivity check modes (ssm_defs.h CONNECTIVITY)
enum CONNECTIVITY { CONNECT_None, CONNECT_Flexible, CONNECT_Strict };

// global tuning knobs (ssm_vxedge.h) — no-ops in the shim
inline void SetMatchPrecision(PRECISION) {}
inline void SetConnectivityCheck(CONNECTIVITY) {}

// SSM structure alignment result + driver. Public fields mirror ssm::Align so
// Coot's superpose code compiles; all methods are inert.
class Align {
 public:
  mmdb::mat44    TMatrix;                 // superposition matrix (identity => no move)
  mmdb::realtype rmsd = 0, Qscore = 0, ncombs = 0, seqIdentity = 0;
  int            cnCheck = 0;
  int            nres1 = 0, nres2 = 0;    // residues in each structure
  int            nsel1 = 0, nsel2 = 0;    // residues in aligned selections
  int            nalgn = 0, ngaps = 0, nmd = 0;
  int            selHndCa1 = 0, selHndCa2 = 0;
  mmdb::ivector  Ca1 = nullptr, Ca2 = nullptr;   // C-alpha correspondence vectors
  mmdb::rvector  dist1 = nullptr;                // optimised C-alpha distances

  Align() {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) TMatrix[i][j] = (i == j) ? 1.0 : 0.0;
  }
  ~Align() {}

  // no-op alignment: report "no hits" with an identity transform and zero counts,
  // so callers skip the Ca1/Ca2 correspondence loops (bounded by nsel*/nalgn == 0).
  int align(mmdb::PManager, mmdb::PManager, PRECISION, CONNECTIVITY,
            int /*selHnd1*/ = 0, int /*selHnd2*/ = 0) { return RC_NoHits; }
  int AlignSelectedMatch(mmdb::PManager, mmdb::PManager, PRECISION, CONNECTIVITY,
                         int /*selHnd1*/ = 0, int /*selHnd2*/ = 0, int /*nselect*/ = 0) {
    return RC_NoHits;
  }
  int GetNMatches() const { return 0; }
};
typedef Align *PAlign;

}  // namespace ssm

#endif  // COOT_MMDB_SHIM_SSM_ALIGN_H
