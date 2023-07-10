/*! \file buccaneer-ncsbuild.h buccaneer library */
/* (C) 2002-2007 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_NCSBUILD_H
#define BUCCANEER_NCSBUILD_H
#pragma once

//! Class for build NCS related chains using density
class Ca_ncsbuild {
 public:
  Ca_ncsbuild( double reliability = 0.5, double rmsd = 1.0, int nmin = 12 ) : reliability_(reliability), rmsd_(rmsd), nmin_(nmin) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const;
 private:
  double reliability_, rmsd_;
  int nmin_;
};

#endif
