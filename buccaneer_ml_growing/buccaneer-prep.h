/*! \file buccaneer-prep.h buccaneer library */
/* (C) 2002-2008 Kevin Cowtan & University of York all rights reserved */

#ifndef BUCCANEER_PREP
#define BUCCANEER_PREP

#include "buccaneer-prot.h"


//! Class for merging overlapped Ca chains and grouping by symmetry
class Ca_prep {
 public:
  struct Rama_flt { double phi, psi, rad; };

  Ca_prep( double main_tgt_rad, double side_tgt_rad, Rama_flt rama_flt, bool correl, bool seqnc, bool debug=false ) : rama_flt_(rama_flt), main_tgt_rad_(main_tgt_rad), side_tgt_rad_(side_tgt_rad), correl_(correl), seqnc_(seqnc), debug_(debug) {}
  bool operator() ( LLK_map_target& llktgt, std::vector<LLK_map_target>& llkcls, const clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const;

  // ramachandran filter data
  static const Rama_flt rama_flt_all, rama_flt_helix, rama_flt_strand, rama_flt_nonhelix;

  static void set_cpus( int cpus ) { ncpu = cpus; }
 private:
  Rama_flt rama_flt_;
  double main_tgt_rad_, side_tgt_rad_;
  bool correl_, seqnc_, debug_;
  static int ncpu;
};


//! class for growing for Ca groups
class Prep_threaded : public clipper::Thread_base {
 public:
  Prep_threaded() {}
  Prep_threaded( std::vector<LLK_map_target>& targets, const clipper::Xmap<float>& xmap, const std::vector<std::vector<clipper::RTop_orth> >& rtops );
  void prep( const int& count );
  const std::vector<LLK_map_target> result() const { return targets_; }
  //! run single or multi-threaded
  bool operator() ( int nthread = 0 );
  //! merge results from multiple threads
  void merge( const Prep_threaded& other );
 private:
  void Run();        //!< the thread 'Run' method
  static int count;  //!< Thread control parameter
  // all data required for calculation is stored in the class
  std::vector<LLK_map_target> targets_;
  const clipper::Xmap<float>* xmap_;
  const std::vector<std::vector<clipper::RTop_orth> > rtops_;
  std::vector<bool> done;
};

#endif
