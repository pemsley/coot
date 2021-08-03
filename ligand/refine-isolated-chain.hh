#ifndef REFINE_ISOLATED_CHAIN_HH
#define REFINE_ISOLATED_CHAIN_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "utils/ctpl.h"
#include "geometry/protein-geometry.hh"

namespace coot {
   void refine_isolated_chain(mmdb::Chain *chain_p, mmdb::Manager *mol_for_this_chain, const coot::protein_geometry &geom,
                              ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight,
                              const clipper::Xmap<float> &xmap);

}


#endif // REFINE_ISOLATED_CHAIN_HH
