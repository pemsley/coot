/*! \file buccaneer-filter.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_FILTER_H
#define BUCCANEER_FILTER_H
#pragma once

//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_filter {
 public:
  Ca_filter( double sig_cut = 3.0 ) : sigcut(sig_cut) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const;

  static bool filter( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, double sigcut, bool keep=true );
  static bool filter( clipper::MiniMol& mol, double sigcut );
 private:
  double sigcut;
};

#endif
