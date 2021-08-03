#ifndef RAMA_RSR_EXTEND_FRAGMENTS_HH
#define RAMA_RSR_EXTEND_FRAGMENTS_HH

#include <clipper/core/xmap.h>
#include <mmdb2/mmdb_manager.h>
#include "utils/ctpl.h"
#include "geometry/protein-geometry.hh"

// when this function updates mol, the update_count is updated.
// (this might need tricky locking of mol - but maybe not though - because the
//  chains are being extendended - not deleted here)
//
void rama_rsr_extend_fragments(mmdb::Manager *mol, const clipper::Xmap<float> &xmap, ctpl::thread_pool  *thread_pool_p, unsigned int n_threads,
                               float weight, unsigned int n_phi_psi_trials, const coot::protein_geometry &geom,
                               unsigned int *update_count);

#endif // RAMA_RSR_EXTEND_FRAGMENTS_HH
