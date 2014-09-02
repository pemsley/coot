/*! \file nautilus-join.h nautilus library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#ifndef NAUTILUS_JOIN_H
#define NAUTILUS_JOIN_H


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


class NucleicAcidJoin {
 public:
  static clipper::MiniMol join( const clipper::MiniMol& mol );
 private:
  class Node { public: float score; std::vector<int> ptrs; };
  static std::vector<int> best_chain( std::vector<Node>& nodes );
};


#endif
