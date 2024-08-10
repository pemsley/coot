/*
 * coot-utils/merge-C-and-N-terminii.hh
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
#ifndef MERGE_C_AND_N_TERMINII_HH
#define MERGE_C_AND_N_TERMINII_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

namespace coot {

   // not by finding overlapping fragments, try to merge using close N and C terminii and fitting
   // a possible missing residue or two between the N and C terminii and using symmetry
   //
   // use_symmetry, if set, will cause Coot to try to link symmetry-related fragment
   // but this is not coded up yet.
   //
   void merge_C_and_N_terminii(mmdb::Manager *mol,
                               const clipper::Xmap<float> &xmap,
                               bool use_symmetry=true, bool using_missing_loop_fit=true);

   void merge_C_and_N_terminii_0_gap(mmdb::Manager *mol); // don't use a map
   void merge_C_and_N_terminii_0_gap(mmdb::Manager *mol, const clipper::Xmap<float> &xmap);

}


#endif // MERGE_C_AND_N_TERMINII_HH
