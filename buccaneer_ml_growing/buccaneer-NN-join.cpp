/*! \file buccaneer-join.cpp buccaneer library */
/* (C) 2006-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-join.h"

#include <clipper/clipper-contrib.h>
#include "buccaneer-NN-Select.h"
#include <algorithm>
#include <set>

/*
struct Tri_residue {
  int flag;
  double score;
  clipper::String type[3];
  clipper::Coord_orth res[3][3];
};
*/

// find best single chain ignoring loops
std::vector<int> Ca_join::best_chain( std::vector<Node>& nodes )
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


bool Ca_join::join( clipper::MiniMol& mol, const double& rmerg, const double& rjoin, const clipper::Coord_orth& com )
{
  typedef clipper::MMonomer Mm;
  double r2merg = rmerg*rmerg;
  double r2join = rjoin*rjoin;

  const clipper::MiniMol& mol1 = mol;
  const clipper::Cell       cell = mol1.cell();
  const clipper::Spacegroup spgr = mol1.spacegroup();

  clipper::Coord_frac comf = com.coord_frac(cell);

  // create 3-residue segments
  std::vector<Tri_residue> fragments;
  Tri_residue fragment;
  for ( int chn = 0; chn < mol1.size(); chn++ ) {
    for ( int res = 1; res < mol1[chn].size()-1; res++ ) {
      if ( Mm::protein_peptide_bond( mol1[chn][res-1], mol1[chn][res] ) &&
           Mm::protein_peptide_bond( mol1[chn][res], mol1[chn][res+1] ) ) {
        int n0 = mol1[chn][res-1].lookup( " N  ", clipper::MM::ANY );
        int a0 = mol1[chn][res-1].lookup( " CA ", clipper::MM::ANY );
        int c0 = mol1[chn][res-1].lookup( " C  ", clipper::MM::ANY );
        int n1 = mol1[chn][res  ].lookup( " N  ", clipper::MM::ANY );
        int a1 = mol1[chn][res  ].lookup( " CA ", clipper::MM::ANY );
        int c1 = mol1[chn][res  ].lookup( " C  ", clipper::MM::ANY );
        int n2 = mol1[chn][res+1].lookup( " N  ", clipper::MM::ANY );
        int a2 = mol1[chn][res+1].lookup( " CA ", clipper::MM::ANY );
        int c2 = mol1[chn][res+1].lookup( " C  ", clipper::MM::ANY );
        if ( n0 >= 0 && a0 >= 0 && c0 >= 0 &&
             n1 >= 0 && a1 >= 0 && c1 >= 0 &&
             n2 >= 0 && a2 >= 0 && c2 >= 0 ) {
          fragment.type[0] = mol1[chn][res-1].type();
          fragment.type[1] = mol1[chn][res  ].type();
          fragment.type[2] = mol1[chn][res+1].type();
          fragment.res[0][0] = mol1[chn][res-1][n0].coord_orth();
          fragment.res[0][1] = mol1[chn][res-1][a0].coord_orth();
          fragment.res[0][2] = mol1[chn][res-1][c0].coord_orth();
          fragment.res[1][0] = mol1[chn][res  ][n1].coord_orth();
          fragment.res[1][1] = mol1[chn][res  ][a1].coord_orth();
          fragment.res[1][2] = mol1[chn][res  ][c1].coord_orth();
          fragment.res[2][0] = mol1[chn][res+1][n2].coord_orth();
          fragment.res[2][1] = mol1[chn][res+1][a2].coord_orth();
          fragment.res[2][2] = mol1[chn][res+1][c2].coord_orth();
          fragment.flag = 1;
          fragment.score = 1.0;
          // upweight sequenced fragments, and flag core seqeunced regions.
          if ( mol1[chn][res].type() != "UNK" ) {
            int r1, r2;
            for ( r1 = res-1; r1 >= 0; r1-- )
              if ( mol1[chn][r1].type() == "UNK" ) break;
            for ( r2 = res+1; r2 < mol1[chn].size(); r2++ )
              if ( mol1[chn][r2].type() == "UNK" ) break;
            int d = std::min( res-r1, r2-res );
            fragment.score += 0.1 * double(d);  // upweight sequenced
            if ( d > 15 ) fragment.flag = 2;    // flag core sequenced
          }
          fragments.push_back( fragment );
        } // if mainchain atoms present
      } // if connected residues
    } // loop over residues
  } // loop over chains


  // now merge equivalent fragments
  for ( int f1 = 0; f1 < fragments.size(); f1++ ) {  // NOT -1 (sign compare)
    clipper::Coord_frac cx0 = fragments[f1].res[0][1].coord_frac(cell);
    clipper::Coord_frac cx1 = fragments[f1].res[1][1].coord_frac(cell);
    clipper::Coord_frac cx2 = fragments[f1].res[2][1].coord_frac(cell);
    for ( int f2 = f1+1; f2 < fragments.size(); f2++ ) {
      if ( fragments[f1].flag == 1 && fragments[f2].flag == 1 ) {
        clipper::Coord_frac cy0 = fragments[f2].res[0][1].coord_frac(cell);
        clipper::Coord_frac cy1 = fragments[f2].res[1][1].coord_frac(cell);
        clipper::Coord_frac cy2 = fragments[f2].res[2][1].coord_frac(cell);
        cy0 = cy0.symmetry_copy_near( spgr, cell, cx1 );
        cy1 = cy1.symmetry_copy_near( spgr, cell, cx1 );
        cy2 = cy2.symmetry_copy_near( spgr, cell, cx1 );
        if ( ( cy0 - cx0 ).lengthsq(cell) < r2merg &&
             ( cy1 - cx1 ).lengthsq(cell) < r2merg &&
             ( cy2 - cx2 ).lengthsq(cell) < r2merg ) {
          double s1 = fragments[f1].score / (fragments[f1].score+1.0);
          double s2 =                 1.0 / (fragments[f1].score+1.0);
          for ( int r = 0; r < 3; r++ )
            for ( int a = 0; a < 3; a++ ) {
              cy1 = fragments[f2].res[r][a].coord_frac(cell);
              cy1 = cy1.symmetry_copy_near( spgr, cell, cx1 );
              fragments[f1].res[r][a] =
                s1 * fragments[f1].res[r][a] + s2 * ( cy1.coord_orth(cell) );
            }
          fragments[f1].score += fragments[f2].score;
          fragments[f2].score = 0.0;
          fragments[f2].flag = 0;

        }
      }
    }
  }





  //NN
  SelectNN::predict(fragments,mol,r2join); // if nnselect keyword not used, nothing will happen here




  // make a list of joins
  std::vector<Node> joins( fragments.size() );
  for ( int f1 = 0; f1 < fragments.size(); f1++ )
    if ( fragments[f1].flag != 0 ) {
      joins[f1].score = 1.0;
      clipper::Coord_frac cx1 = fragments[f1].res[1][1].coord_frac(cell);
      clipper::Coord_frac cx2 = fragments[f1].res[2][1].coord_frac(cell);
      for ( int f2 = 0; f2 < fragments.size(); f2++ )
        if ( fragments[f2].flag != 0 ) {
          if ( f1 != f2 ) {
            clipper::Coord_frac cy0 = fragments[f2].res[0][1].coord_frac(cell);
            clipper::Coord_frac cy1 = fragments[f2].res[1][1].coord_frac(cell);
            clipper::Coord_frac cy2 = fragments[f2].res[2][1].coord_frac(cell);
            cy0 = cy0.symmetry_copy_near( spgr, cell, cx1 );
            cy1 = cy1.symmetry_copy_near( spgr, cell, cx1 );
            cy2 = cy2.symmetry_copy_near( spgr, cell, cx1 );
            if ( (cx1-cy0).lengthsq(cell) < r2join &&
                 (cx2-cy1).lengthsq(cell) < r2join ) {
              if ( fragments[f1].flag == 1 && fragments[f2].flag == 1 )
                joins[f1].ptrs.push_back( f2 );
              else
                if ( f2 == f1+1 )
                  joins[f1].ptrs.push_back( f2 );
            }
          }
        }
    }

  /*
  for ( int j1 = 0; j1 < joins.size(); j1++ ) {
    std::cout << j1 << "(" << fragments[j1].flag << "):\t";
    for ( int j2 = 0; j2 < joins[j1].size(); j2++ ) {
      std::cout << joins[j1][j2] << "(" << fragments[joins[j1][j2]].flag << ")\t";
    }
    std::cout << std::endl;
  }
  */

  // use threading to extract successive longest chains
  std::vector<std::vector<int> > chns;
  while (1) {
    // get longest remaining chain
    std::vector<int> chn = best_chain( joins );
    if ( chn.size() < 6 ) break;
    // add longest chain to list
    chns.push_back( chn );
    // remove used fragments
    for ( int r = 0; r < chn.size(); r++ )
      fragments[chn[r]].flag = 0;
    // remove links from used fragments
    for ( int f = 0; f < joins.size(); f++ )
      if ( fragments[f].flag == 0 )
        joins[f].ptrs.clear();
    // and links to used fragments
    for ( int f = 0; f < joins.size(); f++ )
      for ( int j = joins[f].ptrs.size()-1; j >= 0; j-- )
        if ( fragments[joins[f].ptrs[j]].flag == 0 )
          joins[f].ptrs.erase( joins[f].ptrs.begin() + j );
  }

  /*
  for ( int c = 0; c < chns.size(); c++ ) {
    std::cout << c << ":\t";
    for ( int r = 0; r < chns[c].size(); r++ ) {
      std::cout << chns[c][r] << " ";
    }
    std::cout << std::endl;
  }
  */

  // now join the fragments
  char atomid[3][5] = { " N  ", " CA ", " C  " };
  char atomel[3][2] = { "N", "C", "C" };

  clipper::MiniMol mol2( spgr, cell );
  for ( int c = 0; c < chns.size(); c++ ) {
    // chain and atom info
    clipper::MPolymer chain;
    int ires = 1;
    clipper::MAtom atom = clipper::Atom::null();
    atom.set_occupancy(1.0);
    atom.set_u_iso( clipper::Util::nan() );

    const std::vector<int>& chn = chns[c];

    clipper::Coord_frac cx, cy;
    cx = fragments[chn[0]].res[0][1].coord_frac(cell);  // reference coord
    cx = cx.symmetry_copy_near( spgr, cell, comf );

    for ( int f = -1; f < int(chn.size())+1; f++ ) {
      // residue info
      clipper::MMonomer residue;
      residue.set_type("UNK");

      // add this residue
      for ( int a = 0; a < 3; a++ ) {
        clipper::Coord_orth co( 0.0, 0.0, 0.0 );
        double s = 0.0;
        for ( int r = -1; r <= 1; r++ )
          if ( f+r >= 0 && f+r < chn.size() ) {
            s += 1.0;
            cy = fragments[chn[f+r]].res[1-r][a].coord_frac(cell);
            cy = cy.symmetry_copy_near( spgr, cell, cx );
            co += cy.coord_orth(cell);
          }
        co = (1.0/s) * co;
        atom.set_element( atomel[a] );
        atom.set_id( atomid[a] );
        atom.set_coord_orth( co );
        residue.insert( atom );
      }
      residue.set_seqnum( ires++ );
      chain.insert( residue );

      cx = residue[1].coord_orth().coord_frac(cell);  // update reference coord
    }
    mol2.insert( chain );
  }

  // tidy up the peptide bonds
  for ( int chn = 0; chn < mol2.size(); chn++ )
    for ( int res = 0; res < mol2[chn].size()-1; res++ )
      ProteinTools::tidy_peptide_bond( mol2[chn][res], mol2[chn][res+1] );

  // globularise
  ProteinTools::globularise( mol2, comf );

  // restore the residue types, if any
  ProteinTools::copy_residue_types( mol2, mol1 );

  mol = mol2;
  return true;
}


// build chains by merging and joining tri-residue fragments
bool Ca_join::operator() ( clipper::MiniMol& mol ) const
{
  const clipper::Cell       cell = mol.cell();
  const clipper::Spacegroup spgr = mol.spacegroup();

  // first calculate a convenient ASU centre for output of results
  // (cosmetic only)
  // calc mask
  clipper::Resolution reso( 2.0 );
  clipper::Grid_sampling grid( spgr, cell, reso, 1.0 );
  clipper::Xmap<float> xmap( spgr, cell, grid ), xflt( spgr, cell, grid );
  clipper::EDcalc_mask<float> maskcalc( 2.0 );
  maskcalc( xmap, mol.atom_list() );
  // calc smoothing radius
  double rad = 0.5 * pow( cell.volume()/spgr.num_symops(), 0.333 );
  clipper::MapFilterFn_linear fn( rad );
  clipper::MapFilter_fft<float>
    fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( xflt, xmap );
  // find peak
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  MRI iy = xflt.first();
  for ( MRI ix = xflt.first(); !ix.last(); ix.next() )
    if ( xflt[ix] > xflt[iy] ) iy = ix;
  clipper::Coord_orth com = iy.coord_orth();

  return join( mol, rmerg, rjoin, com );
}
