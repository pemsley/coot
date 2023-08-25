/*! \file buccaneer-correct.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_CORRECT_H
#define BUCCANEER_CORRECT_H
#pragma once

//! Class for correct Ca chains using density
class Ca_correct {
 public:
  Ca_correct( int torsion_sampling = 12 ) : torsion_sampling_(torsion_sampling) {}
  //! rebuild chain to fix insertions/deletions
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq );
  int num_corrected() const;
 private:
  //! score a residue
  static double score_chain_position( const clipper::MMonomer& mm, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget );
  //! score a chain
  static double score_chain_sequence( const clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget );
  //! rebuild a chain at a particular point and return best score
  clipper::MPolymer best_rebuild_sequence( const int& position, const clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget );

  int torsion_sampling_;
  int num_cor;
};

#endif
