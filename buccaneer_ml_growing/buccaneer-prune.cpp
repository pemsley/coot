/*! \file buccaneer-prune.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-prune.h"

#include <clipper/clipper-contrib.h>
#include "buccaneer-sequence.h"


std::vector<float> Ca_prune::score_positions( const clipper::MPolymer& mp, const std::vector<float>& scr )
{
  std::vector<float> scores( mp.size() );
  //int s0 = mp.size();
  float s0 = 0.0;
  float s1 = 0.0;
  for ( int r = 0; r < mp.size(); r++ ) {
    s0 += scr[r];
    if ( mp[r].type() != "UNK" ) s1 += 1.0;
  }
  for ( int r = 0; r < mp.size(); r++ ) {
    scores[r] = s0;
    if ( mp[r].type() != "UNK" ) scores[r] += 10000.0;
  }
  return scores;
}


bool Ca_prune::prune( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, double rad )
{
  double r2cut = rad*rad;

  clipper::Cell       cell = mol.cell();
  clipper::Spacegroup spgr = mol.spacegroup();

  // split into separate chains
  clipper::MiniMol moltmp = mol;
  ProteinTools::split_chains_at_gap( moltmp );

  // score chains by density
  std::vector<std::vector<float> > s;
  for ( int c = 0; c < moltmp.size(); c++ ) s.push_back( ProteinTools::main_chain_densities( moltmp[c], xmap ) );
  // calculate z-score
  double s0(0.0), s1(0.0), s2(0.0);
  for ( int c = 0; c < s.size(); c++ ) {
    for ( int r = 0; r < s[c].size(); r++ ) {
      s0 += 1.0;
      s1 += s[c][r];
      s2 += clipper::Util::sqr(s[c][r]);
    }
  }
  s0 = std::max( s0, 1.0 );
  s1 /= s0; s2 /= s0;
  s2 = sqrt( s2 - s1*s1 );
  for ( int c = 0; c < s.size(); c++ ) for ( int r = 0; r < s[c].size(); r++ ) s[c][r] = ( s[c][r] - s1 ) / s2;
  // calculate sigmoid values
  for ( int c = 0; c < s.size(); c++ ) for ( int r = 0; r < s[c].size(); r++ ) s[c][r] = Ca_sequence::phi_approx( s[c][r] );

  // now loop over chains
  clipper::Coord_frac cf1, cf2;
  for ( int chn1 = 0; chn1 < moltmp.size(); chn1++ ) {
    for ( int chn2 = chn1; chn2 < moltmp.size(); chn2++ ) {
      std::vector<float> score1 = score_positions( moltmp[chn1], s[chn1] );
      std::vector<float> score2 = score_positions( moltmp[chn2], s[chn2] );
      // find any clashing residues between these chains
      for ( int res1 = 0; res1 < moltmp[chn1].size(); res1++ ) {
        for ( int res2 = 0; res2 < moltmp[chn2].size(); res2++ ) {
          if ( ( chn1 != chn2 || res1 != res2 ) && 
               moltmp[chn1][res1].type() != "~~~" &&
               moltmp[chn2][res2].type() != "~~~" ) {
            int a1 = moltmp[chn1][res1].lookup( " CA ", clipper::MM::ANY );
            int a2 = moltmp[chn2][res2].lookup( " CA ", clipper::MM::ANY );
            if ( a1 >= 0 && a2 >= 0 ) {
              cf1 = moltmp[chn1][res1][a1].coord_orth().coord_frac(cell);
              cf2 = moltmp[chn2][res2][a2].coord_orth().coord_frac(cell);
              cf2 = cf2.symmetry_copy_near( spgr, cell, cf1 ) - cf1;
              if ( cf2.lengthsq(cell) < r2cut ) {
                // clash found: keep the one with the higher score.
                if ( score1[res1] > score2[res2] )
                  moltmp[chn2][res2].set_type( "~~~" );
                else
                  moltmp[chn1][res1].set_type( "~~~" );
              }
            }
          }
        }
      }
    }
  }

  // eliminate any sequences of less than 6 residues
  mol = clipper::MiniMol( spgr, cell );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < moltmp.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < moltmp[chn].size(); res++ ) {
      if ( moltmp[chn][res].type() != "~~~" ) {
        mp.insert( moltmp[chn][res] );
      } else {
        if ( mp.size() > 5 ) mol.insert( mp );
        mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol.insert( mp );
  }

  return true;
}


bool Ca_prune::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const
{
  return prune( mol, xmap, rad_ );
}
