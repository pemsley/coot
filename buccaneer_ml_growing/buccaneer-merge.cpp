/*! \file buccaneer-merge.cpp buccaneer library */
/* (C) 2009 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-merge.h"


#include "buccaneer-join.h"
#include "buccaneer-sequence.h"
#include "buccaneer-prune.h"
#include "buccaneer-build.h"


bool Ca_merge::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const std::vector<LLK_map_target>& llktarget, const clipper::MMoleculeSequence& seq ) const
{
  // extract the necessary bits of the likelihood targets
  std::vector<LLK_map_target::Sampled> llksample( llktarget.size() );
  for ( int t = 0; t < llktarget.size(); t++ )
    llksample[t] = llktarget[t].sampled();

  // constants
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();
  const clipper::Coord_orth com( 0.0, 0.0, 0.0 );
  const double rad = 2.0;

  // update sequencing
  for ( int c1 = 0; c1 < mol.size(); c1++ ) {
    Ca_sequence::prepare_scores( mol[c1], xmap, llksample );
    Ca_sequence::sequence( mol[c1], seq, reliability_ );
  }
  ProteinTools::split_chains_at_unk( mol, xmap );
  ProteinTools::split_chains_at_gap( mol );

  // now try to augment each chain in turn to increasing the sequencing
  clipper::Coord_frac f1, f2;
  for ( int c1 = 0; c1 < mol.size()-1; c1++ ) {
    std::vector<clipper::Coord_frac> atoms;
    for ( int r1 = 0; r1 < mol[c1].size(); r1++ ) {
        int a1 = mol[c1][r1].lookup( "CA", clipper::MM::ANY );
        if ( a1 >= 0 )
          atoms.push_back( mol[c1][r1][a1].coord_orth().coord_frac(cell) );
    }
    for ( int c2 = c1+1; c2 < mol.size(); c2++ ) {
      std::cout << c1 << " " << c2 << std::endl;

      // test if this chain overlaps the target chain
      bool overlap = false;
      for ( int r2 = 0; r2 < mol[c2].size(); r2++ ) {
        int a2 = mol[c2][r2].lookup( "CA", clipper::MM::ANY );
        if ( a2 >= 0 ) {
          f2 = mol[c2][r2][a2].coord_orth().coord_frac(cell);
          for ( int i1 = 0; i1 < atoms.size(); i1++ ) {
            f1 = atoms[i1];
            f2 = f2.symmetry_copy_near( spgr, cell, f1 );
            double d2 = ( f2 - f1 ).lengthsq( cell );
            if ( d2 < rad*rad ) { overlap = true; break; }
          }
        }
        if ( overlap ) break;
      }

      // chains overlap - try combining them
      if ( overlap ) {
        clipper::MiniMol mol_wrk( spgr, cell );
        mol_wrk.insert( mol[c1] );
        mol_wrk.insert( mol[c2] );
        Ca_join::join( mol_wrk, 2.0, 2.0, com );
        if ( mol_wrk.size() > 0 ) {
          Ca_sequence::prepare_scores( mol_wrk[0], xmap, llksample );
          Ca_sequence::sequence( mol_wrk[0], seq, reliability_ );
          // test if the chain is improved
          const clipper::MPolymer& mp1 = mol[c1];
          const clipper::MPolymer& mp2 = mol_wrk[0];
          int s0, s1;
          s0 = s1 = 0;
          for ( int r = 0; r < mp1.size(); r++ )
            if ( ProteinTools::residue_index_3( mp1[r].type() ) >= 0 ) s0++;
          for ( int r = 0; r < mp2.size(); r++ )
            if ( ProteinTools::residue_index_3( mp2[r].type() ) >= 0 ) s1++;
          // if new chain is better, keep it
          if ( s1 > s0 ) mol[c1] = mol_wrk[0];
        }
      }
    }
  }

  ProteinTools::split_chains_at_gap( mol );  // split chains
  Ca_prune::prune( mol, xmap );
  Ca_build::build( mol, xmap );
  ProteinTools::split_chains_at_gap( mol );  // rename chains
  return true;

  /*
    clipper::MiniMol mol_nb( spgr, cell );
    mol_nb.insert( mol[c1].select("* CA") );  !!!!Correct
    clipper::MAtomNonBond nb( mol_nb, 2.0*rad );
          f1 = mol[c2][r2][a2].coord_orth().coord_frac(cell);
          atoms = nb.atoms_near( mol[c2][r2][a2].coord_orth(), rad );
          for ( int i = 0; i < atoms.size(); i++ ) {
            const clipper::MAtom& atom = mol_nb.atom(atoms[i]);
            f2 = atom.coord_orth().coord_frac(cell);
            f2 = spgr.symop(atoms[i].symmetry()) * f2;
            f2 = f2.lattice_copy_near( f1 );
            double d2 = ( f2 - f1 ).lengthsq( cell );
            if ( d2 < rad*rad ) overlap = true;
          }
  */
}


std::vector<int> Ca_merge::merge_mr( clipper::MiniMol& mol, clipper::MiniMol& mol_mr, double sigcut, int nseed, bool mr_filter, bool mr_seed )
{
  clipper::MPolymer mp;
  const double rad = 3.0;
  clipper::MAtomNonBond nb( mol, rad );

  // result
  std::vector<int> nres(3,0);

  // go through a chain at a time
  for ( int chn = 0; chn < mol_mr.size(); chn++ ) {
    // filter step
    clipper::MPolymer mp0, mp1, mp2;
    mp0 = mol_mr[chn];

    // score the residues
    if ( mr_filter ) {
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
    } else {
      mp1 = mp0;
    }
    
    // now convert to seeds if required
    if ( mr_seed ) {
      for ( int res = 0; res < mp1.size(); res+=nseed ) mp2.insert( mp1[res] );
    } else {
      mp2 = mp1;
    }
    
    // accumulate non-clashing residues
    for ( int res = 0; res < mp2.size(); res++ ) {
      const int cn = mp2[res].lookup( " N  ", clipper::MM::ANY );
      const int ca = mp2[res].lookup( " CA ", clipper::MM::ANY );
      const int cc = mp2[res].lookup( " C  ", clipper::MM::ANY );
      if ( cn >= 0 && ca >= 0 && cc >= 0 ) {
        const int n =
          nb( mp2[res][cn].coord_orth(), rad ).size() +
          nb( mp2[res][ca].coord_orth(), rad ).size() +
          nb( mp2[res][cc].coord_orth(), rad ).size();
        if ( n == 0 ) mp.insert( mp2[res] );
      }
    }

    nres[0] += mp0.size();
    nres[1] += mp1.size();
    nres[2] += mp2.size();
    //std::cout << chn << " " << mp1.size() << " " << mp2.size() << " " << mp.size() << std::endl;
  }

  // convert residues to UNK
  for ( int res = 0; res < mp.size(); res++ ) {
    mp[res].set_type( "UNK" );
    mp[res].protein_sidechain_build_rotamer( 0 );
  }

  // add to input molecule
  if ( mp.size() > 0 ) {
    mol.insert( mp );
    ProteinTools::split_chains_at_gap( mol );
  }

  return nres;
}
