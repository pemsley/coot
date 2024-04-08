/*
 * cootaneer/buccaneer-sequence.h
 *
 * Copyright 2002-2006 by Kevin Cowtan
 * Author: Kevin Cowtan
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


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
