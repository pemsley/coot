/*
 * ligand/refine-isolated-chain.hh
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
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
