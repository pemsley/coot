/*! \file buccaneer-ncsbuild.cpp buccaneer library */
/* (C) 2006-2020 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-ncsbuild.h"


#include "buccaneer-join.h"
#include "buccaneer-sequence.h"
#include "buccaneer-filter.h"


bool Ca_ncsbuild::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const
{
  //clipper::Cell       cell = xmap.cell();
  clipper::Spacegroup spgr = xmap.spacegroup();

  // get center of mass
  clipper::Coord_orth com( 0.0, 0.0, 0.0 );
  clipper::Atom_list atoms = mol.atom_list();
  for ( int i = 0; i < atoms.size(); i++ ) com = com + atoms[i].coord_orth();
  com = ( 1.0 / double( atoms.size() ) ) * com;

  // extract the necessary bits of the likelihood targets
  std::vector<LLK_map_target::Sampled> llksample( llktarget.size() );
  for ( int t = 0; t < llktarget.size(); t++ )
    llksample[t] = llktarget[t].sampled();

  // split into separate chains
  ProteinTools::split_chains_at_gap( mol );

  // now loop over chains
  for ( int chn1 = 0; chn1 < mol.size(); chn1++ ) {
    // assemble combined model for this chain
    clipper::MPolymer mp1, mp2;
    mp1 = mol[chn1];
    clipper::MiniMol mol_wrk;
    mol_wrk.init( mol.spacegroup(), mol.cell() );
    mol_wrk.insert( mp1 );
    // add any other matching chains
    for ( int chn2 = 0; chn2 < mol.size(); chn2++ ) {
      if ( chn2 != chn1 ) {
        clipper::RTop_orth rtop =
          ProteinTools::superpose( mol[chn2], mol[chn1], rmsd_, nmin_, nmin_ );
        if ( !rtop.is_null() ) {
          clipper::MPolymer mp = mol[chn2];
          mp.transform( rtop );
          mol_wrk.insert( mp );
        }
      }
    }
    // were any matches found?
    if ( mol_wrk.size() > 1 ) {
      // remove sequence
      for ( int c = 0; c < mol_wrk.size(); c++ )
        for ( int r = 0; r < mol_wrk[c].size(); r++ )
          mol_wrk[c][r].set_type( "UNK" );
      // join
      Ca_join::join( mol_wrk, 2.0, 2.0, com );
      if ( mol_wrk.size() > 0 ) {  // trap empty models
        if ( mol_wrk[0].size() > mp1.size() ) {  // optimisation
          // sequence
          Ca_sequence::prepare_scores( mol_wrk[0], xmap, llksample );
          Ca_sequence::sequence( mol_wrk[0], seq, reliability_ );
          // filter
          Ca_filter::filter( mol_wrk, xmap, 1.0 );
          // tidy
          ProteinTools::split_chains_at_gap( mol_wrk );
          // make new chain
          if ( mol_wrk.size() > 0 ) mp2 = mol_wrk[0];
          mp2.copy( mp1, clipper::MM::COPY_MP );
          // test if the chain is improved
          int l0, l1, s0, s1;
          l0 = mp1.size();
          l1 = mp2.size();
          s0 = s1 = 0;
          for ( int r = 0; r < mp1.size(); r++ )
            if ( ProteinTools::residue_index_3( mp1[r].type() ) >= 0 ) s0++;
          for ( int r = 0; r < mp2.size(); r++ )
            if ( ProteinTools::residue_index_3( mp2[r].type() ) >= 0 ) s1++;
          // if new chain is better, keep it
          if ( l1 > l0 && s1 > s0 ) {
            mol[chn1] = mp2;
          }
        }
      }
    }
  }

  return true;
}
