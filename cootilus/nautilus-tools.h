/*
 * cootilus/nautilus-tools.h
 *
 * Copyright 2011 by Kevin Cowtan
 * Author: Kevin Cowtan
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#ifndef NAUTILUS_TOOLS_H
#define NAUTILUS_TOOLS_H


#include <clipper/clipper-contrib.h>

#include "nucleicacid_db.h"


class NucleicAcidTools {
 public:
  NucleicAcidTools();
  static int base_index( char c )           { return bindex [int(c)]; }
  static int base_index_translate( char c ) { return bindext[int(c)]; }
  static clipper::MiniMol flag_chains( const clipper::MiniMol& mol );
  static clipper::RTop_orth symmetry_rtop( const std::vector<clipper::Coord_orth>& cowrk, clipper::Coord_orth& coref, const clipper::Spacegroup& spgr, const clipper::Cell& cell );
  static clipper::MiniMol chain_sort( const clipper::MiniMol& mol );
  static clipper::Coord_orth coord_adjust( const clipper::Coord_orth& co, const clipper::Coord_orth& cc3, const clipper::Coord_orth& cf3, const clipper::Coord_orth& cc4, const clipper::Coord_orth& cf4, double rad );
  static bool symm_match( clipper::MiniMol& molwrk, const clipper::MiniMol& molref );
 private:
  class MapFilterFn_g5 : public clipper::MapFilterFn_base { public:
    clipper::ftype operator() ( const clipper::ftype& radius ) const { return exp(-radius*radius/50.0); }
  };
  static int bindex[256], bindext[256];
};

#endif
