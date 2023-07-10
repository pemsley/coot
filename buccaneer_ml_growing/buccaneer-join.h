/*! \file buccaneer-join.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */

#ifndef BUCCANEER_JOIN
#define BUCCANEER_JOIN

#include "buccaneer-prot.h"


//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_join {
 public:
  Ca_join( double rad_merge = 2.0, double rad_join = 2.0 ) : rmerg(rad_merge), rjoin(rad_join) {}
  bool operator() ( clipper::MiniMol& mol ) const;

  static bool join( clipper::MiniMol& mol, const double& rmerg, const double& rjoin, const clipper::Coord_orth& com );

  class Node { public: float score; std::vector<int> ptrs; };
  struct Tri_residue { int flag; double score; clipper::String type[3]; clipper::Coord_orth res[3][3];
  };
 private:
  static std::vector<int> best_chain( std::vector<Node>& nodes );
  double rmerg, rjoin;
};
#endif
