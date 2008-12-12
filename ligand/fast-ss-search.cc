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

#include <clipper/clipper.h>
#include "fast-ss-search.hh"
#include "mini-mol-utils.hh"


namespace coot {


SSfind::SSfind( const double step )
{
  if ( step < 1.0 )
    step_ = clipper::Util::rad2d( step );
  else
    step_ = step;
}


std::vector<std::vector<clipper::Coord_orth> > SSfind::operator() ( const clipper::Xmap<float>& xmap, int num_residues, SSfind::SSTYPE type )
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
  float tb2[][2][3] = { { { 0.00, 0.00, 0.00}, {-1.00, 0.00,-1.75} },
			{ { 0.87, 0.00, 1.23}, { 2.00, 1.25, 0.50} },
			{ { 0.83, 0.00,-1.18}, { 1.75,-1.50,-0.25} } };
  float tb3[][2][3] = { { { 0.00, 0.00, 0.00}, { 2.00, 2.00,-0.25} },
			{ { 0.87, 0.00, 1.23}, { 3.00, 0.50,-0.25} } };
  float tb4[][2][3] = { { { 0.00, 0.00, 0.00}, { 3.25, 0.50,-0.25} },
			{ { 0.87, 0.00, 1.23}, { 3.25, 0.25, 0.00} } };
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
  for ( int i = 0; i < rep_co.size(); i++ )
    for ( int m = 0; m < ssops.size(); m++ )
      all_co.push_back( Pair_coord( ssops[m] * rep_co[i].first,
				    ssops[m] * rep_co[i].second ) );

  // do the search
  std::vector<clipper::RTop_orth> rtops = (*this)( xmap, all_co );

  // build the results
  std::vector<clipper::Coord_orth> ss( ssops.size() );
  std::vector<std::vector<clipper::Coord_orth> > result( rtops.size() );
  for ( int i = 0; i < rtops.size(); i++ ) {
    for ( int r = 0; r < ssops.size(); r++ )
      ss[r] = rtops[i] * clipper::Coord_orth( ssops[r].trn() );
    result[i] = ss;
  }

  return result;
}


std::vector<clipper::RTop_orth> SSfind::operator() ( const clipper::Xmap<float>& xmap, std::vector<std::pair<clipper::Coord_orth,clipper::Coord_orth> >& all_co )
{
  // make a list of rotations
  std::vector<clipper::RTop_orth> rots;
  // make a list of rotation ops to try
  float glim = 360.0;  // gamma
  float blim = 180.0;  // beta
  float alim = 360.0;  // alpha
  // do a uniformly sampled search of orientation space
  float anglim = clipper::Util::min( alim, glim );
  for ( float bdeg=step_/2; bdeg < 180.0; bdeg += step_ ) {
    float beta = clipper::Util::d2rad(bdeg);
    float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step_+1);
    float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step_+1);
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

  // get cutoff (for optimisation)
  clipper::Map_stats stats( xmap );
  double sigcut = stats.mean() + 1.0 * stats.std_dev();

  // get map radius
  double r2( 0.0 ), d2( 0.0 );
  for ( int i = 0; i < all_co.size(); i++ ) {
    d2 = all_co[i].first.lengthsq();
    if ( d2 > r2 ) r2 = d2;
    d2 = all_co[i].second.lengthsq();
    if ( d2 > r2 ) r2 = d2;
  }
  double rad = sqrt( r2 ) + 1.0;

  // make a 1-d array of gridded density values covering ASU+border
  clipper::Grid_range gr0 = xmap.grid_asu();
  clipper::Grid_range gr1( xmap.cell(), xmap.grid_sampling(), rad );
  clipper::Grid_range mxgr( gr0.min()+gr1.min(), gr0.max()+gr1.max() );
  std::vector<float> mapbox( mxgr.size(), 0.0 );
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  MRI ix( xmap );
  for ( int i = 0; i < mapbox.size(); i++ ) {
    ix.set_coord( mxgr.deindex( i ) );
    mapbox[i] = xmap[ix];
  }

  // make a list of indexed, intergerized, rotated lists
  std::vector<std::vector<std::pair<int,int> > > index_lists;
  clipper::Mat33<> grrot( xmap.operator_orth_grid().rot() );
  int i0 = mxgr.index( clipper::Coord_grid(0,0,0) );
  for ( int r = 0; r < rots.size(); r++ ) {
    std::vector<std::pair<int,int> > tmp;
    for ( int i = 0; i < all_co.size(); i++ ) {
      const clipper::Coord_map c1( grrot*(rots[r]*all_co[i].first  ) );
      const clipper::Coord_map c2( grrot*(rots[r]*all_co[i].second ) );
      tmp.push_back( std::pair<int,int>( mxgr.index( c1.coord_grid() ) - i0,
					 mxgr.index( c2.coord_grid() ) - i0 ) );
    }
    index_lists.push_back( tmp );
  }

  // find ss elements
  std::vector<std::pair<float,std::pair<int,int> > > results;
  float bestcut = 0.0;  // abandon any searches where score < bestcut
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {  // loop over map
    clipper::Coord_grid cg = ix.coord();                  // coord in grid
    const int index0 = mxgr.index( cg );                  // index in list
    if ( mapbox[index0] > sigcut ) {
      int   bestrot = -1;
      float bestscr = bestcut;
      for ( int r = 0; r < index_lists.size(); r++ ) {    // loop over rotns
	const std::vector<std::pair<int,int> >& index_list( index_lists[r] );
	float hi = mapbox[index0+index_list[0].first ];
	float lo = mapbox[index0+index_list[0].second];
	int i = 1;
	while ( hi - lo > bestscr ) {                     // loop over points
	  hi = clipper::Util::min( hi, mapbox[index0+index_list[i].first ] );
	  lo = clipper::Util::max( lo, mapbox[index0+index_list[i].second] );
	  i++;
	  if ( !( i < index_list.size() ) ) break;
	}
	if ( hi - lo > bestscr ) {
	  bestscr = hi - lo;
	  bestrot = r;
	}
      }
      if ( bestrot >= 0 ) {
	std::pair<int,int> rt( bestrot, index0 );
	std::pair<float,std::pair<int,int> > scr( bestscr, rt );
	results.push_back( scr );
	// increase search cutoff if required (for speed)
	float tempcut = 0.4 * bestscr;
	if ( tempcut > bestcut ) bestcut = tempcut;
      }
    }
  }

  // sort
  std::sort( results.begin(), results.end() );
  std::reverse( results.begin(), results.end() );

  // build result list
  scores_.clear();
  std::vector<clipper::RTop_orth> rtops;
  for ( int i = 0; i < results.size(); i++ )
    if ( results[i].first > bestcut ) {
      int ir = results[i].second.first ;
      int it = results[i].second.second;
      clipper::RTop_orth
	rtop( rots[ir].rot(), xmap.coord_orth( mxgr.deindex(it).coord_map() ) );
      scores_.push_back( results[i].first );
      rtops.push_back( rtop );
    }

  return rtops;
}



fast_secondary_structure_search::fast_secondary_structure_search( const clipper::Xmap<float> &xmap, int num_residues, SSTYPE type )
{
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

   std::vector<std::vector<clipper::Coord_orth> > sscoord;
   sscoord = (*this)( xmap, num_residues, type );
   success = ( sscoord.size() > 0 );

   std::string cnames = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
   for ( int c = 0; c < sscoord.size(); c++ ) {
     std::string cname = "";
     if ( c < cnames.length() ) cname = cnames.substr(c,1);
     minimol::fragment mf( cname );
     for ( int r = 0; r < sscoord[c].size(); r++ ) {
       minimol::residue mr( r+1, "UNK" );
       minimol::atom ma( " CA ", "C", sscoord[c][r], "" );
       mr.atoms.push_back( ma );
       mf.residues.push_back( mr );
     }
     mol.fragments.push_back( mf );
   }
}


} //namespace coot
