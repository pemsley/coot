/*! \file buccaneer-sequence.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-sequence.h"

#include <clipper/clipper-contrib.h>

#include <set>

// tag counter
int Ca_sequence::tag = 0;


// cumulative normal distribution function
//double phi(double z)
//{
//   return 0.5 * erfc( -0.70710678118654752440*z );
//}

// approximate cumulative normal distribution function
double Ca_sequence::phi_approx( double z )
{
  double p;
  if ( z < 0.0 )
    p = ( exp(-0.5*z*z) / (1.2533141373*(-z+sqrt(z*z+2.546479089470))) );
  else
    p = 1.0 - ( exp(-0.5*z*z) / (1.2533141373*(z+sqrt(z*z+2.546479089470))) );
  return p;
}

// return fraction of first sequence which overlaps the second
double Ca_sequence::sequence_overlap( const clipper::String& seq1, const clipper::String& seq2 )
{
  int lmin = (seq1.length()<seq2.length()) ? seq1.length() : seq2.length();
  int i1, i2;
  i1 = i2 = 0;
  for ( int i = 0; i < lmin; i++ ) {
    if ( isalpha(seq1[i]) ) i1++;
    if ( isalpha(seq1[i]) && isalpha(seq2[i]) ) i2++;
  }
  return double(i2)/double(i1);
}


// return fraction of sequenced residues which match
double Ca_sequence::sequence_similarity( const clipper::String& seq1, const clipper::String& seq2 )
{
  int lmin = (seq1.length()<seq2.length()) ? seq1.length() : seq2.length();
  int t1, t2, ns, nm;;
  ns = nm = 0;
  for ( int i = 0; i < lmin; i++ ) {
    t1 = ProteinTools::residue_index( seq1.substr( i, 1 ) );
    t2 = ProteinTools::residue_index( seq2.substr( i, 1 ) );
    if ( t1 >= 0 || t2 >= 0 ) {
      ns++;
      if ( t1 == t2 ) nm++;
    }
  }
  if ( ns == 0 ) return 0.0;
  return ( double(nm) / double(ns) );
}


/*! Combine multiple non-conflicting sequence alignments. */
Score_list<clipper::String> Ca_sequence::sequence_combine( const Score_list<clipper::String>& seq, const double& reliability )
{
  Score_list<clipper::String> result = seq;
  int len = seq[0].size();
  clipper::String totseq = std::string( len, '?' );
  clipper::String newseq;
  double totscr = 0.0;
  bool keep = false;
  for ( int i = 0; i < seq.size()-1; i++ ) {
    // ignore sequences which clash with current multisequence
    if ( sequence_overlap( seq[i], totseq ) < 0.40 ) {
      // now check whether this sequence meets the score criterion
      int j;
      for ( j = i+1; j < seq.size(); j++ )
	if ( sequence_overlap( seq[i], seq[j] ) > 0.20 ) break;
      double r = phi_approx( seq.score(i) - seq.score(j) );
      // if score difference good then add matched region to sequence
      if ( r < 1.0-reliability ) {
	clipper::String addseq = seq[i];
	newseq = "";
	for ( int k = 0; k < len; k++ )
	  newseq += ( totseq[k] == '?' ) ? addseq[k] : totseq[k];
	totseq = newseq;
	totscr = -999.0;
	keep = true;
      }
      //  mask the corresponding region
      int k1 = seq[i].find_first_not_of( "?" ) - 3;
      int k2 = seq[i].find_last_not_of( "?" ) + 3;
      newseq = "";
      for ( int k = 0; k < len; k++ )
	if ( k >= k1 && k <= k2 && totseq[k] == '?' ) newseq += 'x';
        else                                          newseq += totseq[k];
      totseq = newseq;
    }
  }
  // change 'x' back to '?'
  newseq = "";
  for ( int k = 0; k < len; k++ )
    newseq += ( totseq[k] == 'x' ) ? '?' : totseq[k];
  // add a new entry if a multisequence was found
  if ( keep ) result.add( totscr, newseq );
  return result;
}


/*! Returning highest scoring subsequence matching the given sequence
  to the supplied LLK scores. */
std::pair<double,std::pair<int,int> > Ca_sequence::sequence_score( const std::vector<std::vector<double> >& scores, const clipper::String& subseq )
{
  // accumulate scores
  std::vector<double> score( scores.size() );
  for ( int r = 0; r < scores.size(); r++ ) {
    int t = ProteinTools::residue_index( subseq.substr( r, 1 ) );
    if ( t >= 0 ) score[r] = scores[r][t];  // add residue z-score
    else          score[r] = 0.0;           // or a penalty of +0.0
  }
  // calculate cumulative score
  std::vector<double> sccum( scores.size()+1 );
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
      // double sc = (sccum[r2]-sccum[r1]) + 1.5*sqrt( double(r2-r1+1 ) );
      double sc = ( sccum[r2]-sccum[r1] ) / pow( 1.0+l1*l1 , 0.25 );
      if ( sc < scmin ) { r1min = r1; r2min = r2; scmin = sc; }
    }
  std::pair<double,std::pair<int,int> > result;
  result.first = scmin;
  result.second.first  = r1min;
  result.second.second = r2min;
  return result;
}


/*! Return a scored list of sequence matches between a set of LLK
  scores and available sequence ranges. */
Score_list<clipper::String> Ca_sequence::sequence_match( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq )
{
  // prepare data for partial sequence match
  clipper::String nulseq = std::string( scores.size(), '?' );
  int minoff = 5;
  // find the best match of this chain against the sequence
  std::pair<double,std::pair<int,int> > result, result_tmp;
  Score_list<std::pair<int,int> > matches_off(100);
  for ( int seqchn = 0; seqchn < seq.size(); seqchn++ ) {
    clipper::String chnseq = nulseq + seq[seqchn].sequence() + nulseq;
    int maxoff = chnseq.length()-scores.size();
    for ( int seqoff = minoff; seqoff <= maxoff-minoff; seqoff++ ) {
      // score the whole sequence
      clipper::String subseq = chnseq.substr( seqoff, scores.size() );
      result = sequence_score( scores, subseq );
      matches_off.add( result.first, std::pair<int,int>(seqchn,seqoff) );
    }
  }

  // refine the matches by insertion and deletion
  Score_list<clipper::String> matches_tmp(200);
  for ( int i = 0; i < matches_off.size(); i++ ) {
    int seqchn = matches_off[i].first;
    int seqoff = matches_off[i].second;
    clipper::String chnseq = nulseq + seq[seqchn].sequence() + nulseq;

    // try the unmutated sequence
    clipper::String subseq = chnseq.substr( seqoff, scores.size() );
    result_tmp = sequence_score( scores, subseq );
    int r1 = result_tmp.second.first;
    int r2 = result_tmp.second.second;
    double scrb = result_tmp.first;
    clipper::String seqb =
      nulseq.substr(0,r1)+subseq.substr(r1,r2-r1)+nulseq.substr(r2);

    // make a list of candidate chain mutations
    std::vector<clipper::String> seqmut;
    // insertions
    for ( int j = r1+minoff; j < r2-minoff; j++ ) {
      subseq = chnseq.substr(seqoff,j)+"+"+chnseq.substr(seqoff+j);
      for ( int k = 0; k <= 1; k++ )
	seqmut.push_back( subseq.substr( k, scores.size() ) );
    }
    // deletions
    for ( int j = r1+minoff; j < r2-minoff; j++ ) {
      subseq = chnseq.substr(seqoff-1,j)+"-"+chnseq.substr(seqoff+j+1);
      for ( int k = 0; k <= 1; k++ )
	seqmut.push_back( subseq.substr( k, scores.size() ) );
    }

    // and try them out
    for ( int j = 0; j < seqmut.size(); j++ ) {
      result_tmp = sequence_score( scores, seqmut[j] );
      double scr = result_tmp.first + 3.0;
      int s1 = result_tmp.second.first;
      int s2 = result_tmp.second.second;
      // if ( scr < scrb && (s2-s1) >= (r2-r1) && s1 < (2*r1+r2)/3 && s2 > (r1+2*r2)/3 ) {
      if ( scr < scrb && s1 <= r1 && s2 >= r2 ) {
	scrb = scr;
	seqb = nulseq.substr(0,s1)+seqmut[j].substr(s1,s2-s1)+nulseq.substr(s2);
      }
    }

    // add the best mutation
    matches_tmp.add( scrb, seqb );
  }

  // eliminate any duplicates
  Score_list<clipper::String> matches(50);
  for ( int i = 0; i < matches_tmp.size(); i++ ) {
    bool clash = false;
    for ( int j = 0; j < matches.size(); j++ )
      if ( sequence_similarity( matches[j], matches_tmp[i] ) > 0.25 )
	clash = true;
    if ( !clash )
      matches.add( matches_tmp.score(i), matches_tmp[i] );
  }

  return matches;
}


/*! Sequence a chain based on the map LLK target, and available
  sequence ranges. */
Score_list<clipper::String> Ca_sequence::sequence_chain( const clipper::MChain& chain, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llktarget, const clipper::MMoleculeSequence& seq )
{
  typedef clipper::MMonomer Mm;
  int nres = chain.size();
  int ntyp = llktarget.size();

  // for each chain, classify each residue against the density
  clipper::Coord_orth coord_n, coord_ca, coord_c;
  std::vector<std::vector<double> > scores;
  for ( int res = 0; res < nres; res++ ) {
    std::vector<double> scores_type( ntyp, 0.0 );
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
      for ( int t = 0; t < ntyp; t++ )
	scores_type[t] = llktarget[t].target( xmap, ca.rtop_beta_carbon() );
    }
    scores.push_back( scores_type );
  }

  // normalise across rows by mean with moving average
  std::vector<double> rscores( nres, 0.0 );
  for ( int r = 0; r < nres; r++ ) {
    for ( int t = 0; t < ntyp; t++ )
      rscores[r] += scores[r][t];
    rscores[r] /= double( ntyp );
  }
  // calc cummulative values
  std::vector<double> rcum( nres+1, 0.0 );
  for ( int r = 0; r < nres; r++ )
    rcum[r+1] = rscores[r] + rcum[r];
  // calc moving average
  int dr = 5;  // nres / 5 + 1;
  for ( int r = 0; r < nres; r++ ) {
    int r1 = clipper::Util::max( r - dr    ,    0 );
    int r2 = clipper::Util::min( r + dr + 1, nres );
    rscores[r] = ( rcum[r2] - rcum[r1] ) / double( r2 - r1 );
  }
  // and correct
  for ( int r = 0; r < nres; r++ )
    for ( int t = 0; t < ntyp; t++ )
      scores[r][t] = scores[r][t] - rscores[r];

  // normalise down columns by mean and variance to get z-scores
  double s0, s1, s2;
  s0 = double( nres );
  for ( int t = 0; t < ntyp; t++ ) {
    s1 = s2 = 0.0;
    for ( int r = 0; r < nres; r++ ) {
      s1 += scores[r][t];
      s2 += scores[r][t]*scores[r][t];
    }
    s1 /= s0;
    s2 /= s0;
    s2 = sqrt( s2 - s1*s1 );
    for ( int r = 0; r < nres; r++ )
      scores[r][t] = ( scores[r][t] - s1 ) / s2;
  }

  // do sequence match
  return sequence_match( scores, seq );
}


/*! Having found a sequence match, apply it to the chain, taking into
  account any existing sequence. */
void Ca_sequence::sequence_apply( clipper::MChain& chain, const clipper::String& seq )
{
  if ( chain.size() != seq.length() ) clipper::Message::message( clipper::Message_fatal( "Sequence: internal error - length mismatch" ) );

  // make old and new sequences
  int m, m1, m2;
  clipper::MChain oldseq, newseq;
  clipper::MMonomer mm;
  for ( m = 0; m < chain.size(); m++ ) {
    mm.set_type( chain[m].type() );
    oldseq.insert( mm );
    int t = ProteinTools::residue_index( seq.substr(m,1) );
    mm.set_type( "UNK" );
    if ( t >= 0 ) mm.set_type( ProteinTools::residue_code_3( t ) );
    if ( seq[m] == '+' ) mm.set_type( "+++" );
    if ( seq[m] == '-' ) mm.set_type( "---" );
    newseq.insert( mm );
  }

  // now find contiguous regions in oldseq trace
  std::vector<std::pair<int,int> > regions;
  m = 0;
  while ( m < oldseq.size() ) {
    while ( m < oldseq.size() ) {
      if ( oldseq[m].type() != "UNK" ) break;
      m++;
    }
    m1 = m;
    while ( m < oldseq.size() ) {
      if ( oldseq[m].type() == "UNK" ) break;
      m++;
    }
    m2 = m;
    if ( m2 > m1 ) regions.push_back( std::pair<int,int>( m1, m2 ) );
  }

  // check each region in turn for clashes
  for ( int i = 0; i < regions.size(); i++ ) {
    bool clash = false;
    m1 = regions[i].first;
    m2 = regions[i].second;
    for ( m = m1; m < m2; m++ )
      if ( newseq[m].type() != "UNK" && newseq[m].type() != oldseq[m].type() )
	clash = true;
    if ( m1 > 0 )
      if ( newseq[m1-1].type() != "UNK" && newseq[m1].type() == "UNK" )
	clash = true;
    if ( m2 < newseq.size() )
      if ( newseq[m2-1].type() == "UNK" && newseq[m2].type() != "UNK" )
	clash = true;
    if ( clash )
      for ( m = m1; m < m2; m++ ) oldseq[m].set_type( "UNK" );
  }

  // combine the remaining types
  for ( m = 0; m < newseq.size(); m++ )
    if ( newseq[m].type() != "UNK" )
      chain[m].set_type( newseq[m].type() );
    else
      chain[m].set_type( oldseq[m].type() );

  /*
  std::cout << "Applying sequence on chain " << chain.id() << " length " << chain.size() << "\n";
  for ( int i = 0; i < regions.size(); i++ ) std::cout << regions[i].first << " " << regions[i].second << "\n";
  for ( m = 0; m < chain.size(); m++ ) std::cout << m << "\t" << oldseq[m].type() << "\t" << newseq[m].type() << "\t" << chain[m].type() << "\n";
  */
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

    // now check for multiple matches
    matches = sequence_combine( matches, reliability_ );

    // test the reliability of the sequence match
    bool apply = false;
    if ( matches.size() >= 2 ) {
      double r = phi_approx( matches.score(0) - matches.score(1) );
      if ( r < 1.0-reliability_ ) apply = true;
    }

    // apply sequence
    if ( apply ) sequence_apply( mol2[chn], matches[0] );

    // flag tag for whether chain sequenced
    int tagx = tag;
    if ( !apply ) tagx += 10000;

    // tag the chains for later reference
    for ( int res = 0; res < mol2[chn].size(); res++ )
      mol2[chn][res].set_property( "TAG", clipper::Property<int>( tagx ) );

    // count sequenced residues
    for ( int res = 0; res < mol2[chn].size(); res++ )
      if ( mol2[chn][res].type() != "UNK" ) num_seq++;

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
