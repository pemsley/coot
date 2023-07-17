/*! \file buccaneer-prune.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_PRUNE_H
#define BUCCANEER_PRUNE_H
#pragma once

//! Class for pruning clashing Ca chains using density
class Ca_prune {
 public:
  Ca_prune( double rad = 3.0 ) : rad_(rad) {}
  static bool prune( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, double rad = 3.0 );
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const;
 private:
  static std::vector<float> score_positions( const clipper::MPolymer& mp, const std::vector<float>& scr );
  clipper::ftype rad_;
};

#endif
