/*! \file cootaneer-sequence.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-sequence.h"


class Coot_sequence {
 public:
  //! constructor: takes filename for resolution dependent LLK coefficients
  Coot_sequence( std::string filename );
  //! get sequence match and confidence for given chain
  std::pair<std::string,double> sequence_chain( const clipper::Xmap<float>& xmap, const std::vector<std::pair<std::string,std::string> >& sequence, CMMDBManager& mmdb, std::string chain_id );
  //! return the confidence of the the best sequence match
  double confidence() const { return prob; }
  //! return the best sequence match
  std::string best_sequence() const { return bestseq; }
  //! return the full sequence of the best sequence match
  std::string full_sequence() const { return fullseq; }
  //! return the chain number corresponding to the best match
  int chain_number() const { return bestchn; }
  //! return the chain sequence offset corresponding to the best match
  int chain_offset() const { return bestoff; }
  //! check if targets have been initialised
  bool is_null() const { return (llksmp.size() > 0); }

  static void write_targets( std::string name , const std::vector<LLK_map_target::Sampled> llksmp );
  static std::vector<LLK_map_target::Sampled> read_targets( std::string name );
 private:
  static std::pair<signed char, signed char> pack( double d );
  static double unpack( std::pair<signed char, signed char> pc );
  // llk data
  std::vector<LLK_map_target::Sampled> llksmp;
  // stored results:
  std::string bestseq, fullseq;
  int bestchn, bestoff;
  double prob;
};
