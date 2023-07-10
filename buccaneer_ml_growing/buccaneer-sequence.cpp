/*! \file buccaneer-sequence.cpp buccaneer library */
/* (C) 2006-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-sequence.h"
#include "buccaneer-tidy.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


int Ca_sequence::ncpu = 0;
bool Ca_sequence::semet_ = false;
clipper::MiniMol Ca_sequence::molprior;


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


// cache scores in residue properties
void Ca_sequence::prepare_score( clipper::MMonomer& mm, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample )
{
  bool cached = false;
  // check if a valid cached result is cached
  Ca_group ca( mm );
  if ( !ca.is_null() ) {
    if ( mm.exists_property( "SEQDAT" ) ) {
      const Sequence_data& sd = static_cast<const clipper::Property<Sequence_data>&>(mm.get_property( "SEQDAT" )).value();
      if ( ( ca.coord_n()  - sd.ca.coord_n()  ).lengthsq() < 1.0e-3 &&
           ( ca.coord_ca() - sd.ca.coord_ca() ).lengthsq() < 1.0e-3 &&
           ( ca.coord_c()  - sd.ca.coord_c()  ).lengthsq() < 1.0e-3 )
        cached = true;
    }
    // if not, calculate a result and cache
    if ( !cached ) {
      if ( mm.exists_property("SEQDAT") ) mm.delete_property("SEQDAT");
      const int ntyp = llksample.size();
      std::vector<double> scores( ntyp, 0.0 );
      for ( int t = 0; t < ntyp; t++ )
        scores[t] = llksample[t].target( xmap, ca.rtop_beta_carbon() );
      Sequence_data sd( ca, scores );
      mm.set_property( "SEQDAT", clipper::Property<Sequence_data>(sd) );
    }
  }
}


// cache scores in residue properties
void Ca_sequence::prepare_scores( clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample )
{
  /*
  for ( int m = 0; m < mp.size(); m++ )
    prepare_score( mp[m], xmap, llksample );
  */
  Sequence_score_threaded seqsc( mp, xmap, llksample );
  seqsc( ncpu );
  mp = seqsc.result();
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
    t1 = ProteinTools::residue_index_translate( seq1[i] );
    t2 = ProteinTools::residue_index_translate( seq2[i] );
    if ( t1 >= 0 || t2 >= 0 ) {
      ns++;
      if ( t1 == t2 ) nm++;
    }
  }
  if ( ns == 0 ) return 0.0;
  return ( double(nm) / double(ns) );
}


/*! Combine multiple non-conflicting sequence alignments. */
std::vector<bool> Ca_sequence::sequence_combine( const Score_list<clipper::String>& seq, const double& reliability )
{
  std::vector<bool> result( seq.size(), false );
  if ( seq.size() == 0 ) return result;
  int len = seq[0].size();
  clipper::String totseq = std::string( len, '?' );
  clipper::String newseq;
  for ( int i = 0; i < seq.size()-1; i++ ) {
    // ignore sequences which clash with current multisequence
    if ( sequence_overlap( seq[i], totseq ) < 0.40 ) {
      // now check whether this sequence meets the score criterion
      double r = 0.0;
      for ( int j = i+1; j < seq.size(); j++ )
        if ( sequence_overlap( seq[i], seq[j] ) > 0.20 ) {
          r = phi_approx( seq.score(i) - seq.score(j) );
          break;
        }
      // if score difference good then add matched region to sequence
      if ( r < 1.0-reliability ) {
        result[i] = true;
        clipper::String addseq = seq[i];
        newseq = "";
        for ( int k = 0; k < len; k++ )
          newseq += ( totseq[k] == '?' ) ? addseq[k] : totseq[k];
        totseq = newseq;
      }
      // mask the corresponding region
      int k1 = seq[i].find_first_not_of( "?" ) - 3;
      int k2 = seq[i].find_last_not_of( "?" ) + 3;
      newseq = "";
      for ( int k = 0; k < len; k++ )
        if ( k >= k1 && k <= k2 && totseq[k] == '?' ) newseq += 'x';
        else                                          newseq += totseq[k];
      totseq = newseq;
    }
  }
  return result;
}


/*! Returning highest scoring subsequence matching the given sequence
  to the supplied LLK scores. */
std::pair<double,std::pair<int,int> > Ca_sequence::sequence_score( const std::vector<std::vector<double> >& scores, const clipper::String& subseq )
{
  // accumulate scores
  std::vector<double> score( scores.size() );
  for ( int r = 0; r < scores.size(); r++ ) {
    score[r] = 0.0;    // UNK or unknown type
    if ( subseq[r] == '+' || subseq[r] == '-' ) {
      score[r] = 3.0;  // insertion/deletion penalty
    } else {
      int t = ProteinTools::residue_index_translate( subseq[r] );
      if ( t >= 0 ) score[r] = scores[r][t];  // add residue z-score
    }
  }
  // calculate cumulative score
  std::vector<double> sccum( scores.size()+1 ), scwt( scores.size()+1 );
  sccum[0] = 0.0;
  for ( int r = 0; r < score.size(); r++ ) sccum[r+1] = sccum[r]+score[r];
  for ( int r = 0; r < scwt.size(); r++ ) {
    double l1 = double(r)/50.0;
    scwt[r] = 1.0/pow(1.0+l1*l1,0.25);
  }
  // now find grestest subsequence score
  int minseq = 1;
  int r1min = 0;
  int r2min = sccum.size()-1;
  double scmin = 0.0;
  for ( int r1 = 0; r1 < sccum.size()-minseq; r1++ )
    for ( int r2 = r1+minseq; r2 < sccum.size(); r2++ ) {
      // double sc = (sccum[r2]-sccum[r1]) + 1.5*sqrt( double(r2-r1+1 ) );
      double sc = ( sccum[r2]-sccum[r1] ) * scwt[r2-r1];
      if ( sc < scmin ) { r1min = r1; r2min = r2; scmin = sc; }
    }
  std::pair<double,std::pair<int,int> > result;
  result.first = scmin;
  result.second.first  = r1min;
  result.second.second = r2min;
  return result;
}


/*! Perform sequence alignment between the given chain scores and the
  given sequence. */
std::vector<clipper::String> Ca_sequence::sequence_align( const std::vector<std::vector<double> >& scores, const clipper::String& seq )
{
  // set up data structures
  int n = scores.size();
  int m = seq.length();
  std::vector<clipper::String> result;
  if ( n < 3 || m < 3 ) return result;

  // construct a sequence alignment matrix from the sequencing data
  clipper::Matrix<double> s(m,n), c(m,n);
  for ( int i = 0; i < m; i++ ) {
    int t = ProteinTools::residue_index_translate( seq[i] );
    for ( int j = 0; j < n; j++ ) s(i,j) = scores[j][t];
  }

  // now make the cumulative matrix and the back pointers
  clipper::Matrix<int> b(m,n,0), e(m,n,0);
  // set up first row and column
  for ( int i = 0; i < m; i++ ) c(i,0) = s(i,0);
  for ( int j = 0; j < n; j++ ) c(0,j) = s(0,j);
  // set up second and third rows
  for ( int i = 1; i < m; i++ ) c(i,1) = c(i-1,0) + s(i,1);
  for ( int j = 1; j < n; j++ ) c(1,j) = c(0,j-1) + s(1,j);
  for ( int j = 1; j < n; j++ ) c(2,j) = c(1,j-1) + s(2,j);
  // now fill the rest of the matrix
  for ( int i = 3; i < m; i++ )
    for ( int j = 2; j < n; j++ ) {
      double s0 = c(i-1,j-1) + s(i,j);
      double s1 = c(i-1,j-2) + s(i,j) + 3.0;
      double s2 = c(i-3,j-2) + s(i,j) + 3.0;
      if ( s0 <= s1 && s0 <= s2 ) {  // sequence follows on
        c(i,j) = s0;
        b(i,j) = 0;
        e(i-1,j-1) = 1;
      } else if ( s1 <= s2 ) {       // skip a residue in chain
        c(i,j) = s1;
        b(i,j) = 1;
        e(i-1,j-2) = 1;
      } else {                       // skip a residue in sequence
        c(i,j) = s2;
        b(i,j) = 2;
        e(i-3,j-2) = 1;
      }
    }

  // make a list of candidate sequences
  for ( int i = 3; i < m; i++ )
    for ( int j = 2; j < n; j++ )
      if ( e(i,j) == 0 ) {  // if this is a sequence end
        std::vector<char> seqv(n,'?');
        int i1 = i;
        int j1 = j;
        while ( i1 >= 0 && j1 >= 0 ) {
          seqv[j1] = seq[i1];
          if ( b(i1,j1) == 1 ) {
            seqv[j1-1] = '+';
            i1 = i1 - 1;
            j1 = j1 - 2;
          } else if ( b(i1,j1) == 2 ) {
            seqv[j1-1] = '-';
            i1 = i1 - 3;
            j1 = j1 - 2;
          } else {
            i1 = i1 - 1;
            j1 = j1 - 1;
          }
        }
        std::string seqs( seqv.begin(), seqv.end() );
        result.push_back( seqs );
      }

  /*
  // diagnostics
  Score_list<clipper::String> scrs( 3 );
  for ( int i = 0; i < result.size(); i++ ) {
    std::pair<double,std::pair<int,int> > scr;
    scr = sequence_score( scores, result[i] );
    scrs.add( scr.first, result[i] );
  }
  std::cout << "DEBUG " << result.size() << " " << scrs.size() << std::endl;
  for ( int i = 0; i < std::min(int(scrs.size()),3); i++ )
    std::cout << i << "\t" << scrs.score(i) << "\t" << scrs[i] << std::endl;
  */

  // return the sequences for scoring
  return result;
}


/*! Return a scored list of sequence matches between a set of LLK
  scores and available sequence ranges. */
Score_list<clipper::String> Ca_sequence::sequence_match( const std::vector<std::vector<double> >& scores, const clipper::MMoleculeSequence& seq )
{
  // loop over chains and get possible alignments
  std::vector<std::pair<double,clipper::String> > matches_tmp;
  clipper::String nulseq = std::string( scores.size(), '?' );
  for ( int seqchn = 0; seqchn < seq.size(); seqchn++ ) {
    // get possible alignments
    std::vector<clipper::String> seqtmp = 
      sequence_align( scores, seq[seqchn].sequence() );
    // score alignments and add to list
    for ( int i = 0; i < seqtmp.size(); i++ ) {
      // get truncated alignment
      clipper::String subseq = seqtmp[i];
      std::pair<double,std::pair<int,int> > result_tmp =
        sequence_score( scores, subseq );
      // truncate
      int r1 = result_tmp.second.first;
      int r2 = result_tmp.second.second;
      double scrb = result_tmp.first;
      clipper::String seqb = nulseq.substr(0,r1)+subseq.substr(r1,r2-r1)+nulseq.substr(r2);
      // store
      std::pair<double,clipper::String> match( scrb, seqb );
      matches_tmp.push_back( match );
    }
  }
  std::sort( matches_tmp.begin(), matches_tmp.end() );

  // eliminate any duplicates
  Score_list<clipper::String> matches(50);
  for ( int i = 0; i < matches_tmp.size(); i++ )
    if ( matches.addable( matches_tmp[i].first ) ) {
      bool clash = false;
      for ( int j = 0; j < matches.size(); j++ )
        if ( sequence_similarity( matches[j], matches_tmp[i].second ) > 0.25 )
          clash = true;
      if ( !clash )
        matches.add( matches_tmp[i].first, matches_tmp[i].second );
    }

  // if first sequence is labelled, add the unmodified chain too.
  // (For use in sequins)
  if ( seq[0].id() == "TEST" && matches.size() > 0 ) {
    clipper::String s0 = seq[0].sequence();
    clipper::String s1 = matches[0];
    if ( s0.length() == s1.length() ) {
      int nmiss = 0;  // check that unmodified chain is different
      for ( int i = 0; i < s0.length(); i++ )
        if ( isupper( s1[i] ) && s1[i] != s0[i] ) nmiss++;
      if ( nmiss > 0 ) {
        std::pair<double,std::pair<int,int> > result_tmp =
          sequence_score( scores, s0 );
        matches.add( result_tmp.first, s0 );
      }
    }
  }

  // return result
  return matches;
}


/*! Sequence a chain based on the map LLK target, and available
  sequence ranges. */
Score_list<clipper::String> Ca_sequence::sequence_chain( clipper::MChain& chain, const clipper::MMoleculeSequence& seq )
{
  int nres = chain.size();

  // for each chain, classify each residue against the density
  std::vector<std::vector<double> > scores( nres );
  for ( int res = 0; res < nres; res++ ) {
    if ( chain[res].exists_property( "SEQDAT" ) ) {
      const Sequence_data& sd = static_cast<const clipper::Property<Sequence_data>&>(chain[res].get_property( "SEQDAT" )).value();
      scores[res] = sd.data;
    }
  }

  // check for valid types
  int ntyp = 0;
  for ( int res = 0; res < nres; res++ )
    ntyp = std::max( ntyp, int(scores[res].size()) );
  if ( ntyp == 0 ) return Score_list<clipper::String>();

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
    int r1 = std::max( r - dr    ,    0 );
    int r2 = std::min( r + dr + 1, nres );
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

  // adjust according to any prior sequence info (heavy atom or model)
  if ( !molprior.is_null() ) {
    std::vector<clipper::MAtomIndexSymmetry> atoms;
    const double nb_rad = 8.0;
    clipper::MAtomNonBond nb( molprior, nb_rad );
    for ( int r = 0; r < nres; r++ ) {
      Ca_group ca( chain[r] );
      if ( !ca.is_null() ) {
        // score against prior
        clipper::Coord_frac cf1, cf2;
        clipper::Coord_orth co = ca.coord_cb();
        atoms = nb.atoms_near( co, nb_rad );
        cf1 = co.coord_frac( molprior.cell() );
        for ( int i = 0; i < atoms.size(); i++ ) {
          const clipper::MAtom& atom = molprior.atom(atoms[i]);
          cf2 = atom.coord_orth().coord_frac( molprior.cell() );
          cf2 = molprior.spacegroup().symop( atoms[i].symmetry() ) * cf2;
          cf2 = cf2.lattice_copy_near( cf1 );
          double d2 = ( cf2 - cf1 ).lengthsq( molprior.cell() );
          if ( d2 <= nb_rad * nb_rad ) {
            int t = atoms[i].monomer();
            double rad = atom.u_iso();
            double w = atom.occupancy();
            double x = sqrt( d2 ) / rad;
            double f = 0.0;
            if      ( x < 1.0 ) f = 1.0 - 0.5 * clipper::Util::sqr( x );
            else if ( x < 2.0 ) f = 0.5 * clipper::Util::sqr( x - 2.0 );
            scores[r][t] = scores[r][t] - w * f;
          }
        }
      }
    }
  }

  // do sequence match
  return sequence_match( scores, seq );
}


/*! Having found a sequence match, apply it to the chain, taking into
  account any existing sequence. */
void Ca_sequence::sequence_apply( clipper::MChain& chain, const Score_list<clipper::String>& seq, const std::vector<bool>& flags )
{
  if ( seq.size() == 0 ) return;
  if ( seq.size() != flags.size() ) clipper::Message::message( clipper::Message_fatal( "Sequence: internal error - length mismatch" ) );
  if ( seq[0].size() != chain.size() ) clipper::Message::message( clipper::Message_fatal( "Sequence: internal error - length mismatch" ) );

  const clipper::String unktyp = "UNK";

  // get old and new sequences
  clipper::String curseq = ProteinTools::chain_sequence( chain );

  // merge sequences in turn
  //std::cout << "OLD: " << curseq << std::endl;
  for ( int i = seq.size()-1; i >= 0; i-- )
    if ( flags[i] ) {
      // check if the new sequence overlaps/contradicts existing sequence
      clipper::String newseq = seq[i];
      int match(0), mismatch(0);
      for ( int r = 0; r < newseq.size(); r++ ) {
        int t1 = ProteinTools::residue_index_translate( curseq[r] );
        int t2 = ProteinTools::residue_index_translate( newseq[r] );
        if ( t1 >= 0 && t2 >= 0 ) {
          if ( t1 == t2 ) match++;
          else mismatch++;
        }
      }
      // if there is a clash, mask the ends
      if ( match == 0 || mismatch > 0 ) {
        int r1 = newseq.find_first_not_of( "?" ) - 3;
        int r2 = newseq.find_last_not_of( "?" ) + 3;
        for ( int r = 0; r < curseq.length(); r++ )
          if ( r >= r1 && r <= r2 ) curseq[r] = '?';
      }
      // apply the sequence
      for ( int r = 0; r < curseq.size(); r++ )
        if ( newseq[r] != '?' ) curseq[r] = newseq[r];
    }
  //std::cout << "NEW: " << curseq << std::endl;

  // apply to chain
  for ( int r = 0; r < chain.size(); r++ ) {
    int t = ProteinTools::residue_index_translate( curseq[r] );
    clipper::String newtype = "UNK";
    if ( t >= 0 )
      newtype = ProteinTools::residue_code_3( t );
    else if ( curseq[r] == '+' )
      newtype = "+++";
    else if ( curseq[r] == '-' )
      newtype = "---";
    chain[r].set_type( newtype );
  }
}


Score_list<clipper::String> Ca_sequence::sequence( clipper::MChain& chain, const clipper::MMoleculeSequence& seq, const double& reliability )
{
  // do the sequencing
  Score_list<clipper::String> matches = sequence_chain( chain, seq );

  // now check for multiple matches
  std::vector<bool> flags = sequence_combine( matches, reliability );

  // apply sequence
  sequence_apply( chain, matches, flags );

  // translate MSE
  if ( semet_ )
    for ( int res = 0; res < chain.size(); res++ )
      if ( chain[res].type() == "MET" ) chain[res].set_type( "MSE" );

  return matches;
}


bool Ca_sequence::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq )
{
  // extract the necessary bits of the likelihood targets
  std::vector<LLK_map_target::Sampled> llksample( llktarget.size() );
  for ( int t = 0; t < llktarget.size(); t++ )
    llksample[t] = llktarget[t].sampled();

  // split into separate chains
  ProteinTools::split_chains_at_gap( mol );

  // score residues
  for ( int chn = 0; chn < mol.size(); chn++ )
    prepare_scores( mol[chn], xmap, llksample );

  // and sequence
  /*
  history = std::vector<Score_list<clipper::String> >( mol.size() );
  for ( int chn = 0; chn < mol.size(); chn++ )
    history[chn] = sequence( mol[chn], seq, reliability_ );  
  */
  Sequence_threaded seqnc( mol, seq, reliability_ );
  seqnc( ncpu );
  mol = seqnc.result();
  history = seqnc.history();

  // break chains where the sequence is broken
  ProteinTools::split_chains_at_unk( mol, xmap );

  // count sequenced residues
  num_seq = 0;
  for ( int chn = 0; chn < mol.size(); chn++ )
    if ( mol[chn].size() > 5 )
      for ( int res = 0; res < mol[chn].size(); res++ )
        if ( mol[chn][res].type() != "UNK" ) num_seq++;

  // trim trailing ends
  ModelTidy::trim( mol, seq );

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
    int nhist = std::min( history[chn].size(), 5 );
    if ( nhist > 0 ) {
      result += "Chain number: " + clipper::String( chn, 4 ) + "    length: " + clipper::String( int(history[chn][0].length()) ) + "\n";
      for ( int res = 0; res < nhist; res++ ) {
        result += history[chn][res] + " \t" +
          clipper::String( history[chn].score(res), 10, 6 ) + "\n";
      }
    }
  }
  return result;
}


/* set a prior model. Complete residues bias residue type by Cbeta
   position, and S/Se atoms weight in favour of CYS/MET/MSE. If the
   chain ID is '!', the radius and weight can be specified an B and
   occ. */
void Ca_sequence::set_prior_model( const clipper::MiniMol& mol )
{
  // make a prior model with 1 chain of 20 residues
  molprior.init( mol.spacegroup(), mol.cell() );
  clipper::MPolymer mp;
  for ( int t = 0; t < 20; t++ ) {  // one for each residue type
    clipper::MMonomer mm;
    mm.set_seqnum( t + 1 );
    mm.set_type( ProteinTools::residue_code_3( t ) );
    mp.insert( mm );
  }

  // set up prior dummy atoms
  clipper::MAtom atomcb, atomsu, atomse, atom;
  atomcb.set_id( " CB " ); atomcb.set_element( "C" );
  atomcb.set_u_iso( 4.0 ); atomcb.set_occupancy( 1.0 );
  atomsu.set_id( " S  " ); atomsu.set_element( "S" );
  atomsu.set_u_iso( 3.0 ); atomsu.set_occupancy( 4.0 );
  atomse.set_id( "SE  " ); atomse.set_element( "SE" );
  atomse.set_u_iso( 4.0 ); atomse.set_occupancy( 4.0 );
  const int t_cys = ProteinTools::residue_index_3( "CYS" );
  const int t_met = ProteinTools::residue_index_3( "MET" );

  // insert Cbeta atoms to bias residue probability
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ ) {
      Ca_group ca( mol[c][r] );
      if ( !ca.is_null() ) {
        int t = ProteinTools::residue_index_3( mol[c][r].type() );
        atom = atomcb;
        atom.set_coord_orth( ca.coord_cb() );
        mp[t].insert( atom );
      }
    }

  // insert Se atoms to bias MET/MSE probability
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      for ( int a = 0; a < mol[c][r].size(); a++ ) {
        clipper::String atm = mol[c][r][a].id().trim();
        clipper::String ele = mol[c][r][a].element();
        bool override = ( mol[c].id() == "!" );
        double u = clipper::Util::u2b( mol[c][r][a].u_iso() );
        double o = mol[c][r][a].occupancy();
        bool su = ( ele == "S"  || atm == "S"  || atm == "SG" || atm == "SD" );
        bool se = ( ele == "SE" || atm == "SE" );
        if ( su ) {
          atom = atomsu;
          atom.set_coord_orth( mol[c][r][a].coord_orth() );
          if ( override ) { atom.set_u_iso( u ); atom.set_occupancy( o ); }
          mp[ t_cys ].insert( atom );
        }
        if ( ( semet_ && se ) || ( !semet_ && su ) ) {
          atom = atomse;
          atom.set_coord_orth( mol[c][r][a].coord_orth() );
          if ( override ) { atom.set_u_iso( u ); atom.set_occupancy( o ); }
          mp[ t_met ].insert( atom );
        }
      }

  // Store the pseudo-chain
  molprior.insert( mp );
}


// thread methods

int Sequence_score_threaded::count = 0;

Sequence_score_threaded::Sequence_score_threaded( clipper::MPolymer& mp, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target::Sampled>& llksample ) : mp_(mp), xmap_(&xmap), llksample_(&llksample)
{
  // flag which chains were grown
  done = std::vector<bool>( mp_.size(), false );

  // init thread count
  count = 0;
}

void Sequence_score_threaded::sequence_score( const int& res )
{
  Ca_sequence::prepare_score( mp_[res], *xmap_, *llksample_ );
  done[res] = true;
}

bool Sequence_score_threaded::operator() ( int nthread )
{
  bool thread = ( nthread > 0 );
  // try running multi-threaded
  if ( thread ) {
    std::vector<Sequence_score_threaded> threads( nthread-1, (*this) );
    run();  for ( int i = 0; i < threads.size(); i++ ) threads[i].run();
    join(); for ( int i = 0; i < threads.size(); i++ ) threads[i].join();
    // check that it finished
    if ( count >= mp_.size() ) {
      for ( int i = 0; i < threads.size(); i++ ) merge( threads[i] );
    } else {
      thread = false;
    }
  }
  // else run in main thread
  if ( !thread ) {
    for ( int res = 0; res < mp_.size(); res++ ) sequence_score( res );
  }
  return true;
}

void Sequence_score_threaded::merge( const Sequence_score_threaded& other )
{
  for ( int res = 0; res < mp_.size(); res++ )
    if ( other.done[res] )
      mp_[res] = other.mp_[res];
}

void Sequence_score_threaded::Run()
{
  while (1) {
    lock();
    int res = count++;
    unlock();
    if ( res >= mp_.size() ) break;
    sequence_score( res );
  }
}


// thread methods

int Sequence_threaded::count = 0;

Sequence_threaded::Sequence_threaded( const clipper::MiniMol& mol, const clipper::MMoleculeSequence& seq, const double& reliability ) : mol_(mol), seq_(seq), reliability_( reliability )
{
  // flag which chains were grown
  done = std::vector<bool>( mol_.size(), false );
  history_ = std::vector<Score_list<clipper::String> >( mol_.size() );

  // init thread count
  count = 0;
}

void Sequence_threaded::sequence( const int& chn )
{
  history_[chn] = Ca_sequence::sequence( mol_[chn], seq_, reliability_ );
  done[chn] = true;
}

bool Sequence_threaded::operator() ( int nthread )
{
  bool thread = ( nthread > 0 );
  // try running multi-threaded
  if ( thread ) {
    std::vector<Sequence_threaded> threads( nthread-1, (*this) );
    run();  for ( int i = 0; i < threads.size(); i++ ) threads[i].run();
    join(); for ( int i = 0; i < threads.size(); i++ ) threads[i].join();
    // check that it finished
    if ( count >= mol_.size() ) {
      for ( int i = 0; i < threads.size(); i++ ) merge( threads[i] );
    } else {
      thread = false;
    }
  }
  // else run in main thread
  if ( !thread ) {
    for ( int chn = 0; chn < mol_.size(); chn++ ) sequence( chn );
  }
  return true;
}

void Sequence_threaded::merge( const Sequence_threaded& other )
{
  for ( int chn = 0; chn < mol_.size(); chn++ )
    if ( other.done[chn] ) {
      mol_[chn] = other.mol_[chn];
      history_[chn] = other.history_[chn];
    }
}

void Sequence_threaded::Run()
{
  while (1) {
    lock();
    int chn = count++;
    unlock();
    if ( chn >= mol_.size() ) break;
    sequence( chn );
  }
}
