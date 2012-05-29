/*! \file ss-find.cpp buccaneer library */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */

#include "nautilus-ss-find.h"


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
	  hi = std::min( hi, mapbox[index0+index_list[i].first ] );
	  lo = std::max( lo, mapbox[index0+index_list[i].second] );
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
    bestcut = std::max( bestscl*bestscr, bestcut );  // optimisation
  }

  // eliminate any results which would have been eliminated by the cutoff
  for ( int i = 0; i < rslts.size(); i++ )
    if ( rslts[i].score < bestcut ) rslts[i] = rsltnull;

  return rslts;
}
