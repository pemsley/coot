/*
 * cootilus/nautilus-ss-find.h
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


#ifndef NAUTILUS_SS_FIND_H
#define NAUTILUS_SS_FIND_H


#include <clipper/clipper-contrib.h>


/*! Result class */
class SearchResult {
 public:
  float score; int rot; int trn;
  bool operator <( SearchResult other ) const { return score < other.score; }
};


//! class for fast secondary structure finding (alternative to fffear)
class SSfind {
 public:
  typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;

  void prep_xmap( const clipper::Xmap<float>& xmap, const double radius );
  void prep_search( const clipper::Xmap<float>& xmap );
  void prep_search( const clipper::Xmap<float>& xmap, const double rhocut, const double radcut, const clipper::Coord_orth centre );
  std::vector<SearchResult> search( const std::vector<Pair_coord>& target_cs, const std::vector<clipper::RTop_orth>& ops, const double rhocut, const double frccut = 0.0 ) const;

 private:
  std::vector<float> mapbox;
  std::vector<int> srctrn;
  clipper::Grid grid;
  clipper::Grid_range mxgr;
  clipper::Mat33<> grrot;
};


#endif
