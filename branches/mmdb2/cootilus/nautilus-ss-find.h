/*! \file nautiluss-find.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


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
