/*! \file buccaneer-build.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_BUILD_H
#define BUCCANEER_BUILD_H
#pragma once

//! Class for building Ca chains using density
class Ca_build {
 public:
 Ca_build( clipper::String newrestype="ALA", bool flexible=false ) : newrestype_(newrestype), flexible_(flexible) {}
  static bool build( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, clipper::String newrestype="ALA", bool flexible=false );
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const;
 private:
  struct Clash { int p1, m1, p2, m2; };
  static void build_rotate_rotamer( clipper::MMonomer& mm, int nr, int nc );
  static std::vector<std::pair<double,std::pair<int,int> > > score_rotamers( const clipper::MMonomer& mm, const clipper::Xmap<float>& xmap, const clipper::Map_stats& xstat, int nconf );
  static std::vector<Clash> find_clashes( const clipper::MiniMol& mol, const double& d );
  static void fix_clash  ( clipper::MMonomer& m1, clipper::MMonomer& m2, const clipper::Xmap<float>& xmap, const clipper::Map_stats& xstat, const double& d, clipper::String newrestype );
  static void fix_clashes( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const clipper::Map_stats& xstat, clipper::String newrestype );
  clipper::String newrestype_;
  bool flexible_;
};

#endif
