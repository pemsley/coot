/*! \file buccaneer-link.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-link.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


bool Ca_link::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget )
{
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();

  // establish map statistics to determine stopping value for building
  double cutoff = llktarget.llk_distribution( 0.005 );

  // now do some rebuilding
  clipper::Coord_frac cf1, cf2;
  num_link = 0;

  // tag the terminal residues
  for ( int c = 0; c < mol.size(); c++ ) {
    int m = mol[c].size()-1;
    clipper::Property<int> p(c);
    mol[c][0].set_property( "CHNID", p );
    mol[c][m].set_property( "CHNID", p );
  }

  // search over possible chains and find short links:
  std::vector<std::pair<double, std::pair<int,int> > > links;
  for ( int chnn = 0; chnn < mol.size(); chnn++ )
    for ( int chnc = 0; chnc < mol.size(); chnc++ ) 
      if ( chnn != chnc && mol[chnn].size() > 3 && mol[chnc].size() > 3 ) {
        int m = mol[chnc].size()-1;
        int in = mol[chnn][0].lookup( " CA ", clipper::MM::ANY );
        int ic = mol[chnc][m].lookup( " CA ", clipper::MM::ANY );
        cf1 = mol[chnn][0][in].coord_orth().coord_frac(cell);
        cf2 = mol[chnc][m][ic].coord_orth().coord_frac(cell);
        cf1 = cf1.symmetry_copy_near( spgr, cell, cf2 );
        double d2 = (cf1-cf2).lengthsq(cell);
        if ( d2 < rlink*rlink ) links.push_back( std::pair<double,std::pair<int,int> >( d2, std::pair<int,int>( chnn, chnc ) ) );
      }

  // sort the links by length
  std::sort( links.begin(), links.end() );

  // now try each link in turn and see if it can be rebuilt
  for ( int i = 0; i < links.size(); i++ ) {
    int chnidn = links[i].second.first;
    int chnidc = links[i].second.second;
    int chnn, chnc;
    chnn = chnc = -1;
    for ( int c = 0; c < mol.size(); c++ ) {
      int m = mol[c].size()-1;
      if ( mol[c][0].exists_property("CHNID") )
        if ( static_cast<const clipper::Property<int>&>(mol[c][0].get_property("CHNID")).value() == chnidn ) chnn = c;
      if ( mol[c][m].exists_property("CHNID") )
        if ( static_cast<const clipper::Property<int>&>(mol[c][m].get_property("CHNID")).value() == chnidc ) chnc = c;
    }
    if ( chnn >= 0 && chnc >= 0 && chnn != chnc ) {
      // transform second chain to match first
      int m = mol[chnc].size()-1;
      int in = mol[chnn][0].lookup( " CA ", clipper::MM::ANY );
      int ic = mol[chnc][m].lookup( " CA ", clipper::MM::ANY );
      cf1 = mol[chnn][0][in].coord_orth().coord_frac(cell);
      cf2 = mol[chnc][m][ic].coord_orth().coord_frac(cell);
      double d2min = 1.0e9;
      clipper::RTop_orth rto( clipper::RTop_orth::identity() );
      for ( int s = 0; s < spgr.num_symops(); s++ ) {
        clipper::Coord_frac cf1s = spgr.symop(s) * cf1;
        clipper::Coord_frac cf1o = cf1s.lattice_copy_near( cf2 );
        double d2 = (cf1o-cf2).lengthsq(cell);
        if ( d2 < d2min ) {
          d2min = d2;
          clipper::RTop_frac rtf( spgr.symop(s).rot(),
                                  spgr.symop(s).trn() + cf1o - cf1s );
          rto = rtf.rtop_orth(cell);
        }
      }
      clipper::MPolymer        mpn = mol[chnn];
      const clipper::MPolymer& mpc = mol[chnc];
      mpn.transform( rto );
      // try to rebuild
      int resnbest = -1;
      int rescbest = -1;
      double scrbest = 1.0e12;
      ProteinLoop::CoordList<8> r8best;
      for ( int resc = mpc.size()-2; resc < mpc.size(); resc++ )
        for ( int resn = 0; resn < 2; resn++ ) 
          if ( resc > 1 && resn < mol[resn].size() - 1 ) {
            int index_cc0 = mpc[resc-1].lookup( " C  ", clipper::MM::ANY );
            int index_cn1 = mpc[resc  ].lookup( " N  ", clipper::MM::ANY );
            int index_ca1 = mpc[resc  ].lookup( " CA ", clipper::MM::ANY );
            int index_ca4 = mpn[resn  ].lookup( " CA ", clipper::MM::ANY );
            int index_cc4 = mpn[resn  ].lookup( " C  ", clipper::MM::ANY );
            int index_cn5 = mpn[resn+1].lookup( " N  ", clipper::MM::ANY );
            if ( index_cc0 >= 0 && index_cn1 >= 0 && index_ca1 >= 0 &&
                 index_ca4 >= 0 && index_cc4 >= 0 && index_cn5 >= 0 ) {
              // rebuild loop
              ProteinLoop pl( torsion_sampling_ );
              std::vector<ProteinLoop::CoordList<8> > r8;
              r8 = pl.rebuild8atoms( mpc[resc-1][index_cc0].coord_orth(),
                                     mpc[resc  ][index_cn1].coord_orth(),
                                     mpc[resc  ][index_ca1].coord_orth(),
                                     mpn[resn  ][index_ca4].coord_orth(),
                                     mpn[resn  ][index_cc4].coord_orth(),
                                     mpn[resn+1][index_cn5].coord_orth() );
              for ( int i = 0; i < r8.size(); i++ ) {
                Ca_group ca1( r8[i][1], r8[i][2], r8[i][3] );
                double s1 = llktarget.llk( xmap, ca1.rtop_from_std_ori() );
                if ( s1 < scrbest ) {
                  Ca_group ca2( r8[i][4], r8[i][5], r8[i][6] );
                  double s2 = llktarget.llk( xmap, ca2.rtop_from_std_ori() );
                  if ( s2 < scrbest ) {
                    scrbest = std::max( s1, s2 );
                    r8best = r8[i];
                    resnbest = resn;
                    rescbest = resc;
                  }
                }
              }
            }
          }
      // now test whether the link is good enough
      if ( scrbest < cutoff ) {
        clipper::MPolymer mp;
        clipper::MAtom ca( clipper::MAtom::null() ), cn, cc;
        ca.set_occupancy( 1.0 );
        ca.set_u_iso( clipper::Util::nan() );
        cn = cc = ca;
        cn.set_element( "N" ); cn.set_id( "N"  );
        ca.set_element( "C" ); ca.set_id( "CA" );
        cc.set_element( "C" ); cc.set_id( "C"  );
        // add first chain
        for ( int r = 0         ; r < rescbest  ; r++ ) mp.insert( mpc[r] );
        // and final residue
        clipper::MMonomer mm0 = mpc[rescbest];
        int i0 = mm0.lookup( " C  ", clipper::MM::ANY );
        if ( i0 >= 0 ) mm0[i0].set_coord_orth( r8best[0] );
        // first interpolated residue
        cn.set_coord_orth( r8best[1] );
        ca.set_coord_orth( r8best[2] );
        cc.set_coord_orth( r8best[3] );
        clipper::MMonomer mm1; mm1.set_type( "UNK" );
        mm1.insert( cn ); mm1.insert( ca ); mm1.insert( cc );
        // second interpolated residue
        cn.set_coord_orth( r8best[4] );
        ca.set_coord_orth( r8best[5] );
        cc.set_coord_orth( r8best[6] );
        clipper::MMonomer mm2; mm2.set_type( "UNK" );
        mm2.insert( cn ); mm2.insert( ca ); mm2.insert( cc );
        // first residue of next chain
        clipper::MMonomer mm3 = mpn[resnbest];
        int i3 = mm3.lookup( " N  ", clipper::MM::ANY );
        if ( i3 >= 0 ) mm3[i3].set_coord_orth( r8best[7] );
        // add the joining residues
        mp.insert( mm0 );
        mp.insert( mm1 );
        mp.insert( mm2 );
        mp.insert( mm3 );
        // and the remaining residues
        for ( int r = resnbest+1; r < mpn.size(); r++ ) mp.insert( mpn[r] );
        // store the first chain and delete the second
        int chnlo = ( chnn < chnc ) ? chnn : chnc ;
        int chnhi = ( chnn < chnc ) ? chnc : chnn ;
        clipper::MiniMol moltmp( mol.spacegroup(), mol.cell() );
        for ( int c = 0; c < mol.size(); c++ )
          if      ( c == chnlo ) moltmp.insert( mp );
          else if ( c != chnhi ) moltmp.insert( mol[c] );
        mol = moltmp;
        // and go back for another chain
        num_link++;
      } // if we build a new join
    } // if we found the chains to join 
  } // loop over possible joins

  // untag the terminal residues
  for ( int c = 0; c < mol.size(); c++ )
    for ( int r = 0; r < mol[c].size(); r++ )
      if ( mol[c][r].exists_property( "CHNID" ) )
        mol[c][r].delete_property( "CHNID" );

  return true;
}


int Ca_link::num_linked() const
{
  return num_link;
}
