/*! \file buccaneer-filter.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-filter.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


bool Ca_filter::filter( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, double sigcut, bool keep )
{
  clipper::MiniMol mol1 = mol;

  // determine sigma cutoff based on map
  //clipper::Map_stats stats( xmap );
  //double s1 = stats.mean();
  //double s2 = stats.std_dev();

  // now apply it a chain at a time
  for ( int chn = 0; chn < mol1.size(); chn++ ) {
    // score the residues
    std::vector<float> scores =
      ProteinTools::main_chain_densities( mol1[chn], xmap, 5 );
    double s0(0), s1(0), s2(0);
    for ( int res = 0; res < scores.size(); res++ ) {
      s0 += 1.0;
      s1 += scores[res];
      s2 += scores[res]*scores[res];
    }
    s2 = sqrt( s2*s0 - s1*s1 )/std::max(s0,1.0);
    s1 = s1/std::max(s0,1.0);
    for ( int res = 0; res < scores.size(); res++ )
      scores[res] = ( scores[res] - s1 ) / s2;
    
    // mark residues in poor density
    for ( int res = 0; res < mol1[chn].size(); res++ )
      if ( scores[res] <= -sigcut )
        if ( !keep || mol1[chn][res].type() == "UNK" )
          mol1[chn][res].set_type( "~~~" );
  }

  // eliminate any sequences of less than 6 residues
  clipper::MiniMol mol2( mol1.spacegroup(), mol1.cell() );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < mol1.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < mol1[chn].size(); res++ ) {
      if ( mol1[chn][res].type() != "~~~" ) {
        mp.insert( mol1[chn][res] );
      } else {
        if ( mp.size() > 5 ) mol2.insert( mp );
        mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol2.insert( mp );
  }

  mol = mol2;
  return true;
}


bool Ca_filter::filter( clipper::MiniMol& mol, double sigcut )
{
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    // filter step
    clipper::MPolymer mp0, mp1;
    mp0 = mol[chn];

    // score the residues
    std::vector<float> scores = ProteinTools::main_chain_u_values( mp0, 3 );
    double s0(0), s1(0), s2(0);
    for ( int res = 0; res < scores.size(); res++ ) {
      s0 += 1.0;
      s1 += scores[res];
      s2 += scores[res]*scores[res];
    }
    s2 = sqrt( s2*s0 - s1*s1 )/std::max(s0,1.0);
    s1 = s1/std::max(s0,1.0);
    for ( int res = 0; res < scores.size(); res++ )
      scores[res] = ( scores[res] - s1 ) / s2;

    // store residues with good B factors
    for ( int res = 0; res < mp0.size(); res++ )
      if ( scores[res] <= sigcut )
        mp1.insert( mp0[res] );

    // store
    mol[chn] = mp1;
  }

  return true;
}


bool Ca_filter::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const
{
  return filter( mol, xmap, sigcut );
}
