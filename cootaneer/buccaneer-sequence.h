/*! \file buccaneer-sequence.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"


//! Class for sequence Ca chains using density
class Ca_sequence {
 public:
  Ca_sequence( double reliability = 0.5 ) : reliability_(reliability) {}
  bool operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq );
  int num_sequenced() const;
  clipper::String format() const;

  class History {
  public:
    void append( const Ca_sequence& data );
    clipper::String format( const clipper::MiniMol& mol ) const;
  private:
    std::vector<std::pair<int,Score_list<clipper::String> > > history;
  };

  static double phi_approx( double z );
  static double sequence_overlap( const clipper::String& seq1, const clipper::String& seq2 );
  static double sequence_similarity( const clipper::String& seq1, const clipper::String& seq2 );
  Score_list<clipper::String> sequence_combine( const Score_list<clipper::String>& seq, const double& reliability );
  static std::pair<double,std::pair<int,int> > sequence_score( const std::vector<std::vector<double> >& scores, const clipper::String& subseq );
  static Score_list<clipper::String> sequence_match( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq );
  static Score_list<clipper::String> sequence_chain( const clipper::MChain& chain, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llktarget, const clipper::MMoleculeSequence& seq );
  static void sequence_apply( clipper::MChain& chain, const clipper::String& seq );

 private:
  double reliability_;
  int num_seq;
  std::vector<std::pair<int,Score_list<clipper::String> > > history;
  static int tag;
};
