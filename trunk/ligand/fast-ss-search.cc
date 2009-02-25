/* ligand/fast-ss-search.cc
 * 
 * Copyright 2008 The University of York
 * Author: Kevin Cowtan
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

// Portability gubbins
#ifndef _MSC_VER
#include <unistd.h> // for getopt(3)
#else
#define stat _stat
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <algorithm> // for sort

#include <clipper/clipper.h>
#include "fast-ss-search.hh"
#include "mini-mol-utils.hh"


namespace coot {


SSfind::Target::Target( SSfind::SSTYPE type, int num_residues )
{
  const int sslen = num_residues;

  // set up targets
  float ta2[][2][3] = { { { 0.00, 0.00, 0.00}, { 2.50, 0.25, 0.00} },
			{ { 0.87, 0.00, 1.23}, { 1.00,-1.75,-0.25} },
			{ { 0.83, 0.00,-1.18}, { 0.25, 1.75, 0.50} } };
  float ta3[][2][3] = { { { 0.00, 0.00, 0.00}, { 3.00, 0.00, 0.00} },
			{ { 0.87, 0.00, 1.23}, {-1.00, 2.00, 0.25} } };
  float ta4[][2][3] = { { { 0.00, 0.00, 0.00}, { 0.75,-2.75, 0.00} },
			{ { 0.87, 0.00, 1.23}, { 1.25,-2.75, 0.00} } };
  float ta3s[][2][3]= { { { 0.00, 0.00, 0.00}, { 3.00, 0.00, 0.00} },
			{ { 0.87, 0.00, 1.23}, {-1.00, 2.00, 0.25} },
			{ { 0.85, 0.45,-0.70}, { 1.50, 2.00 ,1.50} } };
  float tb2[][2][3] = { { { 0.00, 0.00, 0.00}, {-1.00, 0.00,-1.75} },
			{ { 0.87, 0.00, 1.23}, { 2.00, 1.25, 0.50} },
			{ { 0.83, 0.00,-1.18}, { 1.75,-1.50,-0.25} } };
  float tb3[][2][3] = { { { 0.00, 0.00, 0.00}, { 2.00, 2.00,-0.25} },
			{ { 0.87, 0.00, 1.23}, { 3.00, 0.50,-0.25} } };
  float tb4[][2][3] = { { { 0.00, 0.00, 0.00}, { 3.25, 0.50,-0.25} },
			{ { 0.87, 0.00, 1.23}, { 3.25, 0.25, 0.00} } };
  float tb3s[][2][3]= { { { 0.00, 0.00, 0.00}, { 2.00, 2.00,-0.25} },
			{ { 0.87, 0.00, 1.23}, { 3.00, 0.50,-0.25} },
			{ { 0.35, 0.45,-1.10}, { 0.75,-2.75, 0.50} } };
  typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;
  std::vector<Pair_coord> rep_co, all_co;
  double phi0, psi0;
  if ( type == ALPHA2 ) {
    for ( int i = 0; i < sizeof(ta2)/sizeof(ta2[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta2[i][0][0],ta2[i][0][1],ta2[i][0][2]),
	clipper::Coord_orth(ta2[i][1][0],ta2[i][1][1],ta2[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == ALPHA3 ) {
    for ( int i = 0; i < sizeof(ta3)/sizeof(ta3[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta3[i][0][0],ta3[i][0][1],ta3[i][0][2]),
	clipper::Coord_orth(ta3[i][1][0],ta3[i][1][1],ta3[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == ALPHA4 ) {
    for ( int i = 0; i < sizeof(ta4)/sizeof(ta4[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta4[i][0][0],ta4[i][0][1],ta4[i][0][2]),
	clipper::Coord_orth(ta4[i][1][0],ta4[i][1][1],ta4[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == ALPHA3S ) {
    for ( int i = 0; i < sizeof(ta3s)/sizeof(ta3s[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta3s[i][0][0],ta3s[i][0][1],ta3s[i][0][2]),
	clipper::Coord_orth(ta3s[i][1][0],ta3s[i][1][1],ta3s[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  } else if ( type == BETA2 ) {
    for ( int i = 0; i < sizeof(tb2)/sizeof(tb2[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb2[i][0][0],tb2[i][0][1],tb2[i][0][2]),
	clipper::Coord_orth(tb2[i][1][0],tb2[i][1][1],tb2[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else if ( type == BETA3 ) {
    for ( int i = 0; i < sizeof(tb3)/sizeof(tb3[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb3[i][0][0],tb3[i][0][1],tb3[i][0][2]),
	clipper::Coord_orth(tb3[i][1][0],tb3[i][1][1],tb3[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else if ( type == BETA4 ) {
    for ( int i = 0; i < sizeof(tb4)/sizeof(tb4[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb4[i][0][0],tb4[i][0][1],tb4[i][0][2]),
	clipper::Coord_orth(tb4[i][1][0],tb4[i][1][1],tb4[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else if ( type == BETA3S ) {
    for ( int i = 0; i < sizeof(tb3s)/sizeof(tb3s[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(tb3s[i][0][0],tb3s[i][0][1],tb3s[i][0][2]),
	clipper::Coord_orth(tb3s[i][1][0],tb3s[i][1][1],tb3s[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-120.0); psi0 = clipper::Util::d2rad(120.0);
  } else {
    for ( int i = 0; i < sizeof(ta3)/sizeof(ta3[0]); i++ )  // repr coords
      rep_co.push_back( Pair_coord(
        clipper::Coord_orth(ta3[i][0][0],ta3[i][0][1],ta3[i][0][2]),
	clipper::Coord_orth(ta3[i][1][0],ta3[i][1][1],ta3[i][1][2]) ) );
    phi0 = clipper::Util::d2rad(-58.0); psi0 = clipper::Util::d2rad(-47.0);
  }

  // build residue
  clipper::Coord_orth coa( 0.00, 0.00, 0.00 );  //!< std C-a
  clipper::Coord_orth coc( 0.87, 0.00, 1.23 );  //!< std C
  clipper::Coord_orth con( 0.83, 0.00,-1.18 );  //!< std N
  std::vector<clipper::Coord_orth> mm;
  mm.push_back( con );
  mm.push_back( coa );
  mm.push_back( coc );

  // build secondary structure
  const double pi = clipper::Util::pi();
  std::vector<std::vector<clipper::Coord_orth> > mp;
  for ( int i = 0; i < sslen; i++ ) {
    mm[0] = clipper::Coord_orth( mm[0], mm[1], mm[2], 1.32, 1.99, psi0 );
    mm[1] = clipper::Coord_orth( mm[1], mm[2], mm[0], 1.47, 2.15, pi   );
    mm[2] = clipper::Coord_orth( mm[2], mm[0], mm[1], 1.53, 1.92, phi0 );
    mp.push_back( mm );
  }

  // get RTops
  const int ssmid = (sslen-1)/2;
  std::vector<clipper::RTop_orth> ssops( sslen );
  for ( int m = 0; m < mp.size(); m++ )
    ssops[m] = clipper::RTop_orth( mp[ssmid], mp[m] );

  // build whole ss repr coords
  target_cs.clear();
  for ( int i = 0; i < rep_co.size(); i++ )
    for ( int m = 0; m < ssops.size(); m++ )
      target_cs.push_back( Pair_coord( ssops[m] * rep_co[i].first,
				       ssops[m] * rep_co[i].second ) );

  // build ca coords
  calpha_cs.clear();
  for ( int m = 0; m < ssops.size(); m++ )
    calpha_cs.push_back( clipper::Coord_orth( ssops[m].trn() ) );
}


void SSfind::prep_xmap( const clipper::Xmap<float>& xmap, const double radius )
{
  // make a 1-d array of gridded density values covering ASU+border
  grid = xmap.grid_sampling();
  grrot = xmap.operator_orth_grid().rot();
  clipper::Grid_range gr0 = xmap.grid_asu();
  clipper::Grid_range gr1( xmap.cell(), xmap.grid_sampling(), radius );
  mxgr = clipper::Grid_range( gr0.min()+gr1.min(), gr0.max()+gr1.max() );
  mapbox = std::vector<float>( mxgr.size(), 0.0 );

  // make 1d list of densities
  clipper::Xmap<float>::Map_reference_index ix( xmap );
  for ( int i = 0; i < mapbox.size(); i++ ) {
    ix.set_coord( mxgr.deindex( i ) );
    mapbox[i] = xmap[ix];
  }
}


void SSfind::prep_search( const clipper::Xmap<float>& xmap )
{
  // make list of results
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  srctrn.clear();
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    srctrn.push_back( grid.index( ix.coord() ) );
}


void SSfind::prep_search( const clipper::Xmap<float>& xmap, const double rhocut, const double radcut, const clipper::Coord_orth centre )
{
  // make list of results
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  srctrn.clear();
  double r2cut = ( radcut > 0.0 ) ? radcut*radcut : 1.0e20;
  clipper::Coord_frac cf = centre.coord_frac( xmap.cell() );
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    if ( xmap[ix] > rhocut ) {
      clipper::Coord_frac df = ix.coord().coord_frac( xmap.grid_sampling() );
      df = df.symmetry_copy_near( xmap.spacegroup(), xmap.cell(), cf ) - cf;
      double r2 = df.lengthsq( xmap.cell() );
      if ( r2 < r2cut )
	srctrn.push_back( grid.index( ix.coord() ) );
    }
}


std::vector<SearchResult> SSfind::search( const std::vector<Pair_coord>& target_cs, const std::vector<clipper::RTop_orth>& ops, const double rhocut, const double frccut ) const
{
  // make a list of indexed, intergerized, rotated lists
  std::vector<std::vector<std::pair<int,int> > > index_lists;
  int i0 = mxgr.index( clipper::Coord_grid(0,0,0) );
  for ( int r = 0; r < ops.size(); r++ ) {
    clipper::RTop_orth op = ops[r];
    std::vector<std::pair<int,int> > tmp;
    for ( int i = 0; i < target_cs.size(); i++ ) {
      const clipper::Coord_map c1( grrot*(op*target_cs[i].first  ) );
      const clipper::Coord_map c2( grrot*(op*target_cs[i].second ) );
      tmp.push_back( std::pair<int,int>( mxgr.index(c1.coord_grid()) - i0,
					 mxgr.index(c2.coord_grid()) - i0 ) );
    }
    index_lists.push_back( tmp );
  }

  // make list of results
  SearchResult rsltnull = { 0.0, -1, -1 };
  std::vector<SearchResult> rslts( srctrn.size(), rsltnull );
  for ( int i = 0; i < rslts.size(); i++ ) rslts[i].trn = srctrn[i];

  // find ss elements
  float bestcut = 0.0;  // optimisation: abandon searches where score < bestcut
  const float bestscl( frccut ); 
  for ( int i = 0; i < rslts.size(); i++ ) {  // loop over map
    float bestscr = rslts[i].score;
    int   bestrot = rslts[i].rot;
    float bestlim = ( bestscr > bestcut ) ? bestscr : bestcut;
    clipper::Coord_grid cg = grid.deindex( rslts[i].trn );  // coord in grid
    const int index0 = mxgr.index( cg );                    // index in list
    if ( mapbox[index0] > rhocut ) {
      for ( int r = 0; r < index_lists.size(); r++ ) {      // loop over rotns
	const std::vector<std::pair<int,int> >& index_list( index_lists[r] );
	float hi = mapbox[index0+index_list[0].first ];
	float lo = mapbox[index0+index_list[0].second];
	int i = 1;
	while ( hi - lo > bestlim ) {                     // loop over points
	  hi = clipper::Util::min( hi, mapbox[index0+index_list[i].first ] );
	  lo = clipper::Util::max( lo, mapbox[index0+index_list[i].second] );
	  i++;
	  if ( !( i < index_list.size() ) ) break;
	}
	if ( hi - lo > bestlim ) {
	  bestlim = bestscr = hi - lo;
	  bestrot = r;
	}
      }
    }
    rslts[i].score = bestscr;  // store
    rslts[i].rot   = bestrot;
    bestcut = clipper::Util::max( bestscl*bestscr, bestcut );  // optimisation
  }

  // eliminate any results which would have been eliminated by the cutoff
  for ( int i = 0; i < rslts.size(); i++ )
    if ( rslts[i].score < bestcut ) rslts[i] = rsltnull;

  return rslts;
}


// the wrapper class for coot

int fast_secondary_structure_search::join_offset( const std::vector<clipper::Coord_orth>& frag1, const std::vector<clipper::Coord_orth>& frag2 )
{
  double d2min = 1.0e20;
  int a1min = 0;
  int a2min = 0;
  for ( int a1 = 0; a1 < frag1.size(); a1++ )
    for ( int a2 = 0; a2 < frag2.size(); a2++ ) {
      double d2 = ( frag1[a1] - frag2[a2] ).lengthsq();
      if ( d2 < d2min ) {
	d2min = d2;
	a1min = a1;
	a2min = a2;
      }
    }
  return a2min - a1min;
}


double fast_secondary_structure_search::join_score( const std::vector<clipper::Coord_orth>& frag1, const std::vector<clipper::Coord_orth>& frag2 )
{
  int offset = join_offset( frag1, frag2 );
  double score = 0.0;
  for ( int a1 = 0; a1 < frag1.size(); a1++ ) {
    int a2 = a1 + offset;
    if ( a2 >= 0 && a2 < frag2.size() ) {
      double d2 = ( frag1[a1] - frag2[a2] ).lengthsq() / 1.9;
      if ( d2 < 1.0 ) score += ( 1.0 - d2*d2 );
    }
  }
  return score;
}


std::vector<std::vector<clipper::Coord_orth> > fast_secondary_structure_search::join( std::vector<std::vector<clipper::Coord_orth> > frags, int extn )
{
  // first extend the fragments to increase overlap
  // can we do it?
  for ( int f = 0; f < frags.size(); f++ ) if ( frags[f].size() < 4 ) extn = 0;
  // extend fragments
  for ( int e = 0; e < extn; e++ )
    for ( int f = 0; f < frags.size(); f++ ) {
      using clipper::Coord_orth;
      const std::vector<Coord_orth> fr = frags[f];
      int l = fr.size();
      Coord_orth c0( fr[2], fr[1], fr[0],
		     Coord_orth::length (fr[1],fr[0]),
		     Coord_orth::angle  (fr[2],fr[1],fr[0]),
		     Coord_orth::torsion(fr[3],fr[2],fr[1],fr[0]) );
      Coord_orth c1( fr[l-3], fr[l-2], fr[l-1],
		     Coord_orth::length (fr[l-2],fr[l-1]),
		     Coord_orth::angle  (fr[l-3],fr[l-2],fr[l-1]),
		     Coord_orth::torsion(fr[l-4],fr[l-3],fr[l-2],fr[l-1]) );
      frags[f].resize( l+2 );
      frags[f][ 0 ] = c0;
      for ( int a = 0; a < l; a++ ) frags[f][a+1] = fr[a];
      frags[f][l+1] = c1;
    }

  // create weight lists
  std::vector<std::vector<double> > wghts( frags.size() );
  for ( int f = 0; f < wghts.size(); f++ ) {
    wghts[f].resize( frags[f].size() );
    double da = 0.5*(wghts[f].size()-1);
    double sa = 0.5*(wghts[f].size()+1-2*extn);
    for ( int a = 0; a < wghts[f].size(); a++ )
      wghts[f][a] = clipper::Util::max( 1.0-fabs((double(a)-da)/sa), 1.0e-6 );
  }

  // find best pair to merge
  clipper::Matrix<double> scores(frags.size(),frags.size(),0.0);
  for ( int f1 = 0; f1 < frags.size(); f1++ )
    for ( int f2 = f1+1; f2 < frags.size(); f2++ )
      scores(f1,f2) = join_score( frags[f1], frags[f2] );

  // do iterative merge
  for ( int i = 0; i < frags.size()*frags.size(); i++ ) {
    // find best score
    int j1(0), j2(0);
    double js = 0.0;
    for ( int f1 = 0; f1 < frags.size(); f1++ )
      for ( int f2 = f1+1; f2 < frags.size(); f2++ )
	if ( scores(f1,f2) > js ) { js = scores(f1,f2); j1=f1; j2=f2; }
    if ( js < 1.0 ) break;
    // merge them
    std::vector<clipper::Coord_orth> newfrag;
    std::vector<double> newwght;
    int offset = join_offset( frags[j1], frags[j2] );
    int i1 = clipper::Util::min(0,-offset);
    int i2 = clipper::Util::max(frags[j1].size(),frags[j2].size()-offset);
    for ( int a1 = i1; a1 < i2; a1++ ) {
      int a2 = a1 + offset;
      clipper::Coord_orth co( 0.0, 0.0, 0.0 );
      double wt = 0.0;
      if ( a1 >= 0 && a1 < frags[j1].size() ) {
	co += wghts[j1][a1] * frags[j1][a1];
	wt += wghts[j1][a1];
      }
      if ( a2 >= 0 && a2 < frags[j2].size() ) {
	co += wghts[j2][a2] * frags[j2][a2];
	wt += wghts[j2][a2];
      }
      newfrag.push_back( (1.0/wt) * co );
      newwght.push_back( wt );
    }
    frags[j1] = newfrag; frags[j2].clear();
    wghts[j1] = newwght; wghts[j2].clear();
    // update scores
    for ( int f2 = j1+1; f2 < frags.size(); f2++ )
      scores(j1,f2) = join_score( frags[j1], frags[f2] );
    for ( int f2 = j2+1; f2 < frags.size(); f2++ )
      scores(j2,f2) = 0.0;
  }

  // remove fragment extensions
  std::vector<std::vector<clipper::Coord_orth> > frags_new;  
  for ( int f = 0; f < frags.size(); f++ )
    if ( frags[f].size() > 2*extn ) {
      std::vector<clipper::Coord_orth> fr( frags[f].size() - 2*extn );
      for ( int a = 0; a < fr.size(); a++ ) fr[a] = frags[f][a+extn];
      frags_new.push_back( fr );
    }

  return frags_new;
}


void fast_secondary_structure_search::operator()( const clipper::Xmap<float>& xmap, const clipper::Coord_orth& centre, double radius, std::vector<SSfind::Target> targets )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup& spgr    = xmap.spacegroup();
  const clipper::Cell& cell          = xmap.cell();
  const clipper::Grid_sampling& grid = xmap.grid_sampling();

  // make a list of rotations
  std::vector<clipper::RTop_orth> rots;
  // make a list of rotation ops to try
  const float step = 15.0;
  float glim = 360.0;  // gamma
  float blim = 180.0;  // beta
  float alim = 360.0;  // alpha
  // do a uniformly sampled search of orientation space
  float anglim = clipper::Util::min( alim, glim );
  for ( float bdeg=step/2; bdeg < 180.0; bdeg += step ) {
    float beta = clipper::Util::d2rad(bdeg);
    float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
    float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
    for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
      for ( float thmi=smi/2; thmi < 360.0; thmi += smi ) {
	float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
	float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);
	if ( adeg <= alim && bdeg <= blim && gdeg <= glim ) {
	  float alpha = clipper::Util::d2rad(adeg);
	  float gamma = clipper::Util::d2rad(gdeg);
	  clipper::Euler_ccp4 euler( alpha, beta, gamma );
	  rots.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
	}
      }
  }

  // get model radius from targets
  double r2( 0.0 ), d2( 0.0 );
  for ( int t = 0; t < targets.size(); t++ ) {
    const std::vector<Pair_coord>& target_cs = targets[t].target_coords();
    for ( int i = 0; i < target_cs.size(); i++ ) {
      d2 = target_cs[i].first.lengthsq();
      if ( d2 > r2 ) r2 = d2;
      d2 = target_cs[i].second.lengthsq();
      if ( d2 > r2 ) r2 = d2;
    }
  }
  double rad = sqrt( r2 ) + 1.0;

  // get cutoff (for optimisation)
  clipper::Map_stats stats( xmap );
  double sigcut = stats.mean() + 1.0 * stats.std_dev();

  // make a 1-d array of gridded density values covering ASU+border
  prep_xmap( xmap, rad );

  // make a list of results
  prep_search( xmap, sigcut, radius, centre );

  // loop over targets and find:
  std::vector<std::vector<clipper::Coord_orth> > sscoord;
  for ( int t = 0; t < targets.size(); t++ ) {
    const std::vector<Pair_coord>& target_cs = targets[t].target_coords();
    const std::vector<clipper::Coord_orth>& calpha_cs = targets[t].calpha_coords();

    // find ss elements
    std::vector<SearchResult> rslts = search( target_cs, rots, sigcut, 0.4 );

    // filter
    std::vector<SearchResult> rsltf;
    for ( int i = 0; i < rslts.size(); i++ )
      if ( rslts[i].rot>=0 && rslts[i].trn>=0 ) rsltf.push_back( rslts[i] );

    // sort
    std::sort( rsltf.begin(), rsltf.end() );
    std::reverse( rsltf.begin(), rsltf.end() );

    // build result list
    std::vector<clipper::RTop_orth> rtops;
    for ( int i = 0; i < rsltf.size(); i++ ) {
      int ir = rsltf[i].rot;
      int it = rsltf[i].trn;
      clipper::RTop_orth rtop( rots[ir].rot(),
			       xmap.coord_orth(grid.deindex(it).coord_map()) );
      rtops.push_back( rtop );
    }

    // now do some magic to get the chains near the given centre
    clipper::Coord_frac cf = centre.coord_frac( cell );
    for ( int i = 0; i < rtops.size(); i++ ) {
      clipper::RTop_frac rtof = rtops[i].rtop_frac( cell );
      int smin = 0;
      double d2min = 1.0e12;
      for ( int s = 0; s < spgr.num_symops(); s++ ) {
	clipper::Coord_frac cfs( spgr.symop(s) * rtof.trn() );
	cfs = cfs.lattice_copy_near( cf ) - cf;
	double d2 = cfs.lengthsq( cell );
	if ( d2 < d2min ) { smin  = s; d2min = d2; }
      }
      rtof = clipper::RTop_frac( spgr.symop(smin) * rtof );
      clipper::Coord_frac df( rtof.trn() );
      df = df.lattice_copy_near( cf ) - df;
      rtof = clipper::RTop_frac( rtof.rot(), rtof.trn() + df );
      rtops[i] = rtof.rtop_orth( cell );
    }

    // build the results
    std::vector<clipper::Coord_orth> ss( calpha_cs.size() );
    std::vector<std::vector<clipper::Coord_orth> > result( rtops.size() );
    for ( int i = 0; i < rtops.size(); i++ ) {
      for ( int r = 0; r < ss.size(); r++ )
	ss[r] = rtops[i] * calpha_cs[r];
      result[i] = ss;
    }

    // join fragments
    std::vector<std::vector<clipper::Coord_orth> > ssj = join( result, 1 );

    // and add to final model
    for ( int i = 0; i < ssj.size(); i++ ) sscoord.push_back( ssj[i] );
  }

  // NOW PREPARE THE MODEL FOR COOT
  float acell[6];
  acell[0] = xmap.cell().descr().a();
  acell[1] = xmap.cell().descr().b();
  acell[2] = xmap.cell().descr().c();
  acell[3] = clipper::Util::rad2d(xmap.cell().descr().alpha());
  acell[4] = clipper::Util::rad2d(xmap.cell().descr().beta());
  acell[5] = clipper::Util::rad2d(xmap.cell().descr().gamma());
  std::string spacegroup_str_hm = xmap.spacegroup().symbol_hm();

  mol.set_cell(acell);
  mol.set_spacegroup(spacegroup_str_hm);

  success = ( sscoord.size() > 0 );

  std::string cnames = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  for ( int c = 0; c < sscoord.size(); c++ ) {
    std::string cname = "";
    int c1 = c / cnames.length() - 1;
    int c2 = c % cnames.length();
    if ( c1 >= 0 && c1 < cnames.length() ) cname += cnames.substr(c1,1);
    if ( c2 >= 0 && c2 < cnames.length() ) cname += cnames.substr(c2,1);
    minimol::fragment mf( cname );
    for ( int r = 0; r < sscoord[c].size(); r++ ) {
      minimol::residue mr( r+1, "UNK" );
      minimol::atom ma( " CA ", "C", sscoord[c][r], "", 30.0);
      mr.atoms.push_back( ma );
      mf.residues.push_back( mr );
    }
    mol.fragments.push_back( mf );
  }
}


} //namespace coot
