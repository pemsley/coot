/*
 * ligand/rama-rsr-extend-fragments.hh
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
