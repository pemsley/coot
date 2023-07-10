/*! \file buccaneer-util.h buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */


#ifndef BUCCANEER_UTIL
#define BUCCANEER_UTIL

#include <clipper/clipper-minimol.h>


class BuccaneerUtil {
 public:
  static void set_reference( clipper::String& mtz, clipper::String& pdb );
  // coordinate utilities
  static void read_model( clipper::MiniMol& mol, clipper::String file, bool verbose );
};


class BuccaneerLog {
 public:
  BuccaneerLog( clipper::String& title ) : title_(title), currentcpu(0.0) { log(""); }
  void log( const clipper::String& id );
  void log( const clipper::String& id, const clipper::MiniMol& mol, bool view );
  clipper::String log( const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq );
  void xml( const clipper::String& xml ) const;
  void profile();
  std::vector<double>  evaluate(const clipper::MiniMol& mol, const clipper::MiniMol& mol_mr, const clipper::MMoleculeSequence& seq);
 private:
  struct cycdat { int nfrgs, nchns, nseq, nres, nmax, nunq; double cres, cchn; };
  std::vector<cycdat> data;
  std::vector<std::pair<std::string,double> > prof;
  clipper::String title_;
  double currentcpu;
};


#endif
