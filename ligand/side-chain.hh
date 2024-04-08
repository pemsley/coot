/*
 * ligand/side-chain.hh
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
#ifndef SIDE_CHAIN_HH
#define SIDE_CHAIN_HH

#include "geometry/residue-and-atom-specs.hh"
#include "geometry/protein-geometry.hh"

namespace coot {
   void do_180_degree_side_chain_flip(const residue_spec_t &spec, const std::string &alt_conf,
                                      mmdb::Manager *mol, protein_geometry *pg);

}


#endif // SIDE_CHAIN_HH
