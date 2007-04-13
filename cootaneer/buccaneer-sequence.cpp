/*! \file buccaneer-sequence.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-sequence.h"

#include <clipper/clipper-contrib.h>

#include <set>

// tag counter
int Ca_sequence::tag = 0;


Score_list<clipper::String> Ca_sequence::sequence_match( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq )
{
  // prepare data for partial sequence match
  clipper::String nulseq = std::string( scores.size(), '?' );
  // find the best match of this chain against the sequence
  std::vector<double> score( scores.size() );
  std::vector<double> sccum( scores.size()+1 );
  Score_list<clipper::String> matches_tmp(100);
  for ( int seqchn = 0; seqchn < seq.size(); seqchn++ ) {
    clipper::String chnseq = nulseq + seq[seqchn].sequence() + nulseq;
    int maxoff = chnseq.length()-scores.size();
    for ( int seqoff = 0; seqoff <= maxoff; seqoff++ ) {
      // score the whole sequence
      clipper::String subseq = chnseq.substr( seqoff, scores.size() );
      for ( int r = 0; r < scores.size(); r++ ) {
	int t = ProteinTools::residue_index( subseq.substr( r, 1 ) );
	if ( t >= 0 ) score[r] = scores[r][t];  // add residue z-score
	else          score[r] = 0.10;          // or a penalty of +0.10
      }
      // calculate cumulative score
      sccum[0] = 0.0;
      for ( int r = 0; r < score.size(); r++ ) sccum[r+1] = sccum[r]+score[r];
      // now find grestest subsequence score
      int minseq = 1;
      int r1min = 0;
      int r2min = sccum.size()-1;
      double scmin = 0.0;
      for ( int r1 = 0; r1 < sccum.size()-minseq; r1++ )
	for ( int r2 = r1+minseq; r2 < sccum.size(); r2++ ) {
	  double l1 = ( r2 - r1 )/50.0;  // downweight very long sequences
	  double sc = ( sccum[r2]-sccum[r1] ) / pow( 1.0+l1*l1 , 0.25 );
	  if ( sc < scmin ) { r1min = r1; r2min = r2; scmin = sc; }
	}
      subseq = ( nulseq.substr(0,r1min) +
		 subseq.substr(r1min,r2min-r1min) +
		 nulseq.substr(r2min) );
      matches_tmp.add( scmin, subseq );
    }
  }

  // eliminate any duplicates
  Score_list<clipper::String> matches(10);
  for ( int i = 0; i < matches_tmp.size(); i++ ) {
    bool clash = false;
    for ( int j = 0; j < matches.size(); j++ )
      if ( matches[j] == matches_tmp[i] ) clash = true;
    if ( !clash )
      matches.add( matches_tmp.score(i), matches_tmp[i] );
  }

  return matches;
}


Score_list<clipper::String> Ca_sequence::sequence_chain( const clipper::MChain& chain, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llktarget, const clipper::MMoleculeSequence& seq, TYPE type )
{
  typedef clipper::MMonomer Mm;

  // for each chain, classify each residue against the density
  clipper::Coord_orth coord_n, coord_ca, coord_c;
  std::vector<std::vector<double> > scores;
  for ( int res = 0; res < chain.size(); res++ ) {
    std::vector<double> scores_type( llktarget.size(), 0.0 );
    // find ca, c, n
    int index_n  = chain[res].lookup( " N  ", clipper::MM::ANY );
    int index_ca = chain[res].lookup( " CA ", clipper::MM::ANY );
    int index_c  = chain[res].lookup( " C  ", clipper::MM::ANY );
    // if we have all three atoms, then add residue
    if ( index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
      coord_n  = chain[res][index_n].coord_orth();
      coord_ca = chain[res][index_ca].coord_orth();
      coord_c  = chain[res][index_c].coord_orth();
      Ca_group ca( coord_n, coord_ca, coord_c );
      if ( type == NORMAL )
	for ( int t = 0; t < llktarget.size(); t++ )
	  scores_type[t] = llktarget[t].llk( xmap, ca.rtop_beta_carbon() );
      else
	for ( int t = 0; t < llktarget.size(); t++ )
	  scores_type[t] = llktarget[t].correl( xmap, ca.rtop_beta_carbon() );
    }
    scores.push_back( scores_type );
  }

  // normalise down columns by mean and variance to get z-scores
  double s0, s1, s2;
  s0 = double( scores.size() );
  for ( int t = 0; t < llktarget.size(); t++ ) {
    s1 = s2 = 0.0;
    for ( int r = 0; r < scores.size(); r++ ) {
      s1 += scores[r][t];
      s2 += scores[r][t]*scores[r][t];
    }
    s1 /= s0;
    s2 /= s0;
    s2 = sqrt( s2 - s1*s1 );
    for ( int r = 0; r < scores.size(); r++ )
      scores[r][t] = ( scores[r][t] - s1 ) / s2;
  }

  // do sequence match
  return sequence_match( scores, seq );
}


bool Ca_sequence::operator() ( clipper::MiniMol& mol2, const clipper::MiniMol& mol1, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq )
{
  typedef clipper::MMonomer Mm;

  // extract the necessary bits of the likelihood targets
  std::vector<LLK_map_target::Sampled> llksample( llktarget.size() );
  for ( int t = 0; t < llktarget.size(); t++ )
    llksample[t] = llktarget[t].sampled();

  // split into separate chains
  ProteinTools::chain_tidy( mol2, mol1 );

  // set relibility criteria and docking count
  double score_cutoff = -20.0 * reliability_;  // min score
  double screl_cutoff =   0.3 * reliability_;  // min diff between 1st & 2nd
  num_seq = 0;

  // now loop over chains
  for ( int chn = 0; chn < mol2.size(); chn++ ) {
    // calculate possible matches to this chain
    Score_list<clipper::String> matches =
      sequence_chain( mol2[chn], xmap, llksample, seq );

    // test the reliability of the sequence match
    bool apply = false;
    if ( matches.size() >= 2 )
      if ( matches.score(0) < score_cutoff &&
	   (1.0-screl_cutoff)*matches.score(0) < matches.score(1) ) 
	apply = true;

    // apply sequence
    if ( apply ) {
      // copy the residue types from the best match into the model
      clipper::String bestseq = matches[0];
      for ( int res = 0; res < mol2[chn].size(); res++ ) {
	int t = ProteinTools::residue_index( bestseq.substr(res,1) );
	if ( t >= 0 ) {
	  mol2[chn][res].set_type( ProteinTools::residue_code_3( t ) );
	  num_seq++;
	}
      }
    }

    // flag tag for whether chain sequenced
    int tagx = tag;
    if ( !apply ) tagx += 10000;

    // tag the chains for later reference
    for ( int res = 0; res < mol2[chn].size(); res++ )
      mol2[chn][res].set_property( "TAG", clipper::Property<int>( tagx ) );

    // save some info for future output
    std::pair<int,Score_list<clipper::String> > histdat( tagx, matches );
    history.push_back( histdat );

    tag++;
  }

  return true;
}


int Ca_sequence::num_sequenced() const
{
  return num_seq;
}


clipper::String Ca_sequence::format() const
{
  clipper::String result = "";
  for ( int chn = 0; chn < history.size(); chn++ ) {
    result += "Chain number: " + clipper::String( chn, 4 ) + "    length: " + clipper::String( int(history[chn].second[0].length()) ) + "\n";
    for ( int res = 0; res < history[chn].second.size(); res++ ) {
      result += history[chn].second[res] + " \t" +
	clipper::String( history[chn].second.score(res), 10, 6 ) + "\n";
    }
  }
  return result;
}


// history class - for output

void Ca_sequence::History::append( const Ca_sequence& data )
{
  for ( int i = 0; i < data.history.size(); i++ )
    history.push_back( data.history[i] );
}


clipper::String Ca_sequence::History::format( const clipper::MiniMol& mol ) const
{
  clipper::String result = "";
  for ( int h = 0; h < history.size(); h++ ) {
    int tag = history[h].first;
    std::set<clipper::String> ids;
    clipper::String chn = "";
    // search for ids in the final model with this tag
    for ( int c = 0; c < mol.size(); c++ )
      for ( int r = 0; r < mol[c].size(); r++ ) {
	if ( mol[c][r].exists_property("TAG") ) {
	  int t = dynamic_cast<const clipper::Property<int>&>(mol[c][r].get_property("TAG")).value();
	  if ( t == tag ) ids.insert( mol[c].id() );
	}
      }
    if ( ids.size() > 0 ) {
      // get chain letters
      for ( std::set<clipper::String>::iterator i = ids.begin(); i != ids.end();
	    i++ ) chn += (*i) + " ";
      // now print chain infomation
      result += "Sequencing fragment of length: " + clipper::String( int(history[h].second[0].length()) ) + "  Final chain ID(s): " + chn + "\n";
      for ( int s = 0; s < clipper::Util::min(history[h].second.size(),5);
	    s++ ) {
	if ( s == 0 )
	  if ( tag < 10000 ) result += " Accepted: ";
	  else               result += " Rejected: ";
	else result += "           ";
	result += history[h].second[s] + " \t" +
	  clipper::String( history[h].second.score(s), 10, 6 ) + "\n";
      }
    }
  }
  return result;
}
