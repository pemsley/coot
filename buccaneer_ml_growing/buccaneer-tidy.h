/*! \file buccaneer-tidy.h buccaneer library */
/* (C) 2002-2010 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_TIDY
#define BUCCANEER_TIDY

#include "buccaneer-prot.h"


//! Class for tidying model
class ModelTidy {
 public:
 ModelTidy( double rmsd = 1.0, int nmin = 12, clipper::String newrestype="ALA", bool verbose = false ) : rmsd_(rmsd), nmin_(nmin), newrestype_(newrestype), verbose_(verbose) {}
  bool tidy( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq ) const;

  static std::vector<int> chain_renumber( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq );
  static std::vector<int> chain_assign( const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const std::vector<int> seqnums, const double rmsd, const int nmin );
  static bool chain_move( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, std::vector<int> chnnums );
  static bool sequence_correct( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq, const std::vector<int> seqnums, clipper::String newrestype );

  static std::vector<int>      sequence_count( const clipper::MiniMol& mol );
  static clipper::Array2d<int> sequence_flags( const clipper::MiniMol& mol );

  static void trim( clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq );

 private:
  bool update_model( clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, std::vector<int> chnnums ) const;
  static int chain_renumber( clipper::MPolymer& mp, const clipper::MMoleculeSequence& seq );
  static std::vector<int> count_contacts( const clipper::MAtomNonBond& nb, const clipper::MiniMol& mol, const std::vector<int>& chnnums, const clipper::MPolymer mp, const double& rad );
  static bool move_chain( clipper::MPolymer& mp1, const clipper::MPolymer& mp0, const clipper::Spacegroup& spgr, const clipper::Cell& cell );
  static void best_closed_ncs_group( const clipper::Array2d<int>& super, const std::vector<int>& num_seq, const std::vector<int>& active, std::vector<int>& used_best, std::vector<int>& used );
  double rmsd_;
  int nmin_;
  clipper::String newrestype_;
  bool verbose_;
};


#endif
