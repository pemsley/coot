/*! \file buccaneer-link.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_LINK_H
#define BUCCANEER_LINK_H
#pragma once

//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_link {
 public:
  Ca_link( clipper::ftype rad_link = 5.0, int torsion_sampling = 24 ) : torsion_sampling_(torsion_sampling), rlink(rad_link) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget );
  int num_linked() const;

 private:
  int torsion_sampling_;
  double rlink;
  int num_link;
};

#endif
