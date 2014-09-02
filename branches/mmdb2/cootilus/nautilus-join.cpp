/*! Nautilus join */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#include "nautilus-join.h"

#include "nautilus-tools.h"


#include <set>
#include <algorithm>


std::vector<int> NucleicAcidJoin::best_chain( std::vector<Node>& nodes )
{
  // first make a list of loop starts
  std::vector<float> node_score( nodes.size(), 0.0 );
  // loop starts are now marks as 0
  // declare a list of 'dirty' nodes
  std::set<int> dirty;
  for ( int i = 0; i < node_score.size(); i++ )
    if ( node_score[i] == 0.0 ) dirty.insert( i );
  // now propogate the values
  std::vector<int> bck_ptrs( nodes.size(), -1 );
  while ( dirty.size() > 0 ) {
    // get a node from the dirty list and remove it from the list
    std::set<int>::iterator iter = dirty.begin();
    int node = *iter;
    dirty.erase( iter );
    // now check the children
    for ( int j = 0; j < nodes[node].ptrs.size(); j++ ) {
      int next_node = nodes[node].ptrs[j];
      // test whether we've found a longer route
      float next_score = node_score[node] + nodes[next_node].score;
      if ( node_score[next_node] < next_score ) {
	// check here for loops and broken paths
	int back_node = node;
	int back_node_next;
	while ( bck_ptrs[back_node] >= 0 ) {
	  back_node_next = bck_ptrs[back_node];
	  if ( back_node_next == next_node ) break;
	  if ( node_score[back_node_next] >= node_score[back_node] ) break;
	  back_node = back_node_next;
	}
	// if the path to this node is clean, we can update the node
	if ( bck_ptrs[back_node] < 0 ) {
	  // if this is a longer non-looped route, store it
	  node_score[next_node] = next_score;
	  bck_ptrs[next_node] = node;
	  dirty.insert( next_node );
	}
      }
    }
  }
  // we've found all the long routes, now find the longest and back-trace it
  int node_max = 0;
  for ( int i = 1; i < node_score.size(); i++ )
    if ( node_score[i] > node_score[node_max] )
      node_max = i;
  // and back-trace
  std::vector<int> result;
  int node = node_max;
  result.push_back( node );
  while ( bck_ptrs[node] >= 0 ) {
    int next = bck_ptrs[node];
    bck_ptrs[node] = -1;
    node = next;
    result.push_back( node );
  }
  // reverse the list
  std::reverse( result.begin(), result.end() );
  return result;
}


clipper::MiniMol NucleicAcidJoin::join( const clipper::MiniMol& mol )
{
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();

  clipper::MiniMol mol_wrk( spgr, cell ), mol_new( spgr, cell );
  for ( int c = 0; c < mol.size(); c++ )
    if ( mol[c].exists_property( "NON-NA" ) ) mol_new.insert( mol[c] );
    else                                      mol_wrk.insert( mol[c] );

  // sort by size
  clipper::MiniMol mols = NucleicAcidTools::chain_sort( mol_wrk );

  // dissassemble into lists
  std::vector<clipper::MMonomer> nas;
  for ( int c = 0; c < mols.size(); c++ )
    for ( int r = 0; r < mols[c].size(); r++ ) {
      mols[c][r].set_seqnum( r+1 );
      int i1 = mols[c][r].lookup( " C1'", clipper::MM::ANY );
      int i3 = mols[c][r].lookup( " C3'", clipper::MM::ANY );
      int i4 = mols[c][r].lookup( " C4'", clipper::MM::ANY );
      if ( i1 >= 0 && i3 >= 0 && i4 >= 0 ) nas.push_back( mols[c][r] );
    }

  // make nnb model
  clipper::MiniMol molnb( spgr, cell );
  clipper::MPolymer chnnb;
  for ( int r = 0; r < nas.size(); r++ ) {
    int a = nas[r].lookup( " C1'", clipper::MM::ANY );
    clipper::MMonomer mm;
    mm.insert( nas[r][a] );
    chnnb.insert( mm );
  }
  molnb.insert( chnnb );
  clipper::MAtomNonBond nb( molnb, 4.0 );

  // make list of equivalents
  std::vector<int> equivalent( nas.size() );
  for ( int r = 0; r < nas.size(); r++ ) equivalent[r] = r;

  // find equivalents
  double d2(1.0*1.0), j2(8.0*8.0);
  for ( int r1 = 0; r1 < chnnb.size()-1; r1++ ) { // for each NA
    const int a1 = nas[r1].lookup( " C1'", clipper::MM::ANY );
    const int a3 = nas[r1].lookup( " C3'", clipper::MM::ANY );
    const int a4 = nas[r1].lookup( " C4'", clipper::MM::ANY );
    clipper::Coord_frac ca1 = nas[r1][a1].coord_orth().coord_frac(cell);
    clipper::Coord_frac ca3 = nas[r1][a3].coord_orth().coord_frac(cell);
    clipper::Coord_frac ca4 = nas[r1][a4].coord_orth().coord_frac(cell);
    std::vector<clipper::MAtomIndexSymmetry> atoms =
      nb( chnnb[r1][0].coord_orth(), 2.0 );
    //std::cout << " NB " << r1 << ": ";
    for ( int i = 0; i < atoms.size(); i++ ) { // find other equivalent NAs
      int r2 = atoms[i].monomer();
      //std::cout << r2 << " ";
      // if there is a match AND the residues hasn't already been matched
      if ( equivalent[r2] > equivalent[r1] ) {
	const int b1 = nas[r2].lookup( " C1'", clipper::MM::ANY );
	const int b3 = nas[r2].lookup( " C3'", clipper::MM::ANY );
	const int b4 = nas[r2].lookup( " C4'", clipper::MM::ANY );
	clipper::Coord_frac cb1 = nas[r2][b1].coord_orth().coord_frac(cell);
	clipper::Coord_frac cb3 = nas[r2][b3].coord_orth().coord_frac(cell);
	clipper::Coord_frac cb4 = nas[r2][b4].coord_orth().coord_frac(cell);
	cb1 = cb1.symmetry_copy_near( spgr, cell, ca1 );
	cb3 = cb3.symmetry_copy_near( spgr, cell, ca3 );
	cb4 = cb4.symmetry_copy_near( spgr, cell, ca4 );
	//std::cout << "(" << ( cb1 - ca1 ).lengthsq( cell ) << "," << ( cb3 - ca3 ).lengthsq( cell ) << "," << ( cb4 - ca4 ).lengthsq( cell ) << ") "; 
	if ( ( cb1 - ca1 ).lengthsq( cell ) < d2 &&
	     ( cb3 - ca3 ).lengthsq( cell ) < d2 &&
	     ( cb4 - ca4 ).lengthsq( cell ) < d2 )
	  equivalent[r2] = equivalent[r1];
      }
    }
    //std::cout << std::endl;
  }

  // find links
  Node nodenull; nodenull.score = 1.0;
  std::vector<Node> joins( nas.size(), nodenull );
  for ( int r1 = 0; r1 < nas.size()-1; r1++ ) {
    int r2 = r1 + 1;
    if ( nas[r2].seqnum() == nas[r1].seqnum()+1 ) {
      int a4 = nas[r1].lookup( " C4'", clipper::MM::ANY );
      int b3 = nas[r2].lookup( " C3'", clipper::MM::ANY );
      clipper::Coord_frac ca4 = nas[r1][a4].coord_orth().coord_frac(cell);
      clipper::Coord_frac cb3 = nas[r2][b3].coord_orth().coord_frac(cell);
      cb3 = cb3.symmetry_copy_near( spgr, cell, ca4 );
      if ( ( cb3 - ca4 ).lengthsq( cell ) < j2 ) {
	int e1 = equivalent[r1];
	int e2 = equivalent[r2];
	bool found = false;
	for ( int i = 0; i < joins[e1].ptrs.size(); i++ )
	  if ( joins[e1].ptrs[i] == e2 ) found = true;
	if ( !found ) joins[e1].ptrs.push_back( e2 );
      }
    }
  }

  /*
  for ( int r = 0; r < joins.size(); r++ ) {
    std::cout << r << ":\t" << equivalent[r] << "\t: ";
    for ( int i = 0; i < joins[r].ptrs.size(); i++ ) {
      std::cout << joins[r].ptrs[i] << " ";
    }
    std::cout << std::endl;
  }
  */

  // use threading to extract successive longest chains
  std::vector<int> flags( nas.size(), 1 );
  std::vector<std::vector<int> > chns;
  while (1) {
    // get longest remaining chain
    std::vector<int> chn = best_chain( joins );
    if ( chn.size() < 3 ) break;
    // add longest chain to list
    chns.push_back( chn );
    // remove used fragments
    for ( int r = 0; r < chn.size(); r++ )
      flags[chn[r]] = 0;
    // remove links from used fragments
    for ( int f = 0; f < joins.size(); f++ )
      if ( flags[f] == 0 )
	joins[f].ptrs.clear();
    // and links to used fragments
    for ( int f = 0; f < joins.size(); f++ )
      for ( int j = joins[f].ptrs.size()-1; j >= 0; j-- )
	if ( flags[joins[f].ptrs[j]] == 0 )
	  joins[f].ptrs.erase( joins[f].ptrs.begin() + j );
  }

  /*
  for ( int c = 0; c < chns.size(); c++ ) {
    std::cout << c << ":\t";
    for ( int r = 0; r < chns[c].size(); r++ ) std::cout << chns[c][r] << " ";
    std::cout << std::endl;
  }
  */

  // build chains from successive NAs, with symmetry shift
  for ( int c = 0; c < chns.size(); c++ ) {
    clipper::MPolymer mp;
    // set a reference coord to build near
    clipper::Coord_orth cref( clipper::Coord_orth::null() );
    for ( int r = 0; r < chns[c].size(); r++ ) {
      clipper::MMonomer mm = nas[chns[c][r]];
      if ( !cref.is_null() ) {
	// get nearest symmetry copy
	std::vector<clipper::Coord_orth> cwrk;
	for ( int a = 0; a < mm.size(); a++ )
	  cwrk.push_back( mm[a].coord_orth() );
	mm.transform(NucleicAcidTools::symmetry_rtop(cwrk,cref,spgr,cell));
      }
      mp.insert( mm );
      int a = mm.lookup( " O3'", clipper::MM::ANY );
      if ( a < 0 ) a = mm.lookup( " C3'", clipper::MM::ANY );
      if ( a < 0 ) a = mm.lookup( " C4'", clipper::MM::ANY );
      if ( a >= 0 ) cref = mm[a].coord_orth();
    }
    mol_new.insert( mp );
  }

  // Build final molecule
  return mol_new;
}
