/*! \file buccaneer-sequence.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"

#ifndef BUCCANEER_SEQUENCE_H
#define BUCCANEER_SEQUENCE_H
#pragma once

//! Class for sequence Ca chains using density
class Ca_sequence {
 public:
  Ca_sequence( double reliability = 0.5 ) : reliability_(reliability) {}
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq );
  int num_sequenced() const;
  clipper::String format() const;

  static double phi_approx( double z );
  static void prepare_score( clipper::MMonomer& mm, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample );
  static void prepare_scores( clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample );
  static double sequence_overlap( const clipper::String& seq1, const clipper::String& seq2 );
  static double sequence_similarity( const clipper::String& seq1, const clipper::String& seq2 );
  static std::vector<bool> sequence_combine( const Score_list<clipper::String>& seq, const double& reliability );
  static std::pair<double,std::pair<int,int> > sequence_score( const std::vector<std::vector<double> >& scores, const clipper::String& subseq );
  static std::vector<clipper::String> sequence_align( const std::vector<std::vector<double> >& scores, const clipper::String& seq );
  static Score_list<clipper::String> sequence_match( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq );
  static Score_list<clipper::String> sequence_chain( clipper::MChain& chain, const clipper::MMoleculeSequence& seq );
  static void sequence_apply( clipper::MChain& chain, const Score_list<clipper::String>& seq, const std::vector<bool>& flags );
  static Score_list<clipper::String> sequence( clipper::MChain& chain, const clipper::MMoleculeSequence& seq, const double& reliability );

  static void set_semet( bool semet ) { semet_ = semet; }
  static void set_prior_model( const clipper::MiniMol& mol );

  static void set_cpus( int cpus ) { ncpu = cpus; }

  class Sequence_data {
   public:
    Sequence_data() {}
    Sequence_data( const Ca_group& c, const std::vector<double>& d ) :
      ca(c), data(d) {}
    Ca_group ca;
    std::vector<double> data;
  };
 private:
  double reliability_;
  int num_seq;
  std::vector<Score_list<clipper::String> > history;
  static int ncpu;
  static bool semet_;
  static clipper::MiniMol molprior;
};


//! class for sequenceing for Ca groups
class Sequence_score_threaded : public clipper::Thread_base {
 public:
  Sequence_score_threaded() {}
  Sequence_score_threaded( clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample );
  void sequence_score( const int& chn );
  const clipper::MPolymer& result() const { return mp_; }
  //! run single or multi-threaded
  bool operator() ( int nthread = 0 );
  //! merge results from multiple threads
  void merge( const Sequence_score_threaded& other );
 private:
  void Run();        //!< the thread 'Run' method
  static int count;  //!< Thread control parameter
  clipper::MPolymer mp_;
  const clipper::Xmap<float>* xmap_;
  const std::vector<LLK_map_target::Sampled>* llksample_;
  std::vector<bool> done;
};


//! class for sequenceing for Ca groups
class Sequence_threaded : public clipper::Thread_base {
 public:
  Sequence_threaded() {}
  Sequence_threaded( const clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq, const double& reliability );
  void sequence( const int& chn );
  const clipper::MiniMol& result() const { return mol_; }
  const std::vector<Score_list<clipper::String> >& history() const { return history_; }
  //! run single or multi-threaded
  bool operator() ( int nthread = 0 );
  //! merge results from multiple threads
  void merge( const Sequence_threaded& other );
 private:
  void Run();        //!< the thread 'Run' method
  static int count;  //!< Thread control parameter
  clipper::MiniMol mol_;
  clipper::MMoleculeSequence seq_;
  double reliability_;
  std::vector<Score_list<clipper::String> > history_;
  std::vector<bool> done;
};

#endif
