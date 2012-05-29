/*! \file cootilus-build.cpp nautilus library */
/* (C) 2012 Kevin Cowtan & University of York all rights reserved */

#include "cootilus-build.h"

#include "nautilus-target.h"
#include "nautilus-join.h"


Coot_nucleic_acid_build::Coot_nucleic_acid_build( std::string filename )
{
  filename_ = filename;
}


bool Coot_nucleic_acid_build::build( CMMDBManager* mmdb, const clipper::Xmap<float>& xmap, const clipper::Coord_orth& centre, double radius ) const
{
  clipper::MMDBfile* mmdbfile = static_cast<clipper::MMDBfile*>( mmdb );

  // Get reference model
  NucleicAcidTargets natools; NucleicAcidJoin na_join;
  natools.add_pdb( filename_ );
  natools.set_search_sphere( centre, radius );
  natools.init_stats( xmap, 20 );

  // read existing model
  clipper::MiniMol mol_wrk;
  mmdbfile->import_minimol( mol_wrk );
  mmdbfile->DeleteAllModels();

  // force spacegroup/cell to match map
  mol_wrk.init( xmap.spacegroup(), xmap.cell() );

  int n1(0), n2(0);

  // count monomers before
  for ( int c = 0; c < mol_wrk.size(); c++ ) n1 += mol_wrk[c].size();
  
  // find chains
  mol_wrk = natools.find( xmap, mol_wrk, 5, 5, 18.0 );
  // grow chains
  mol_wrk = natools.grow( xmap, mol_wrk, 25, 0.001 );
  // join chains
  mol_wrk = na_join.join( mol_wrk );
  // prune
  mol_wrk = natools.prune( mol_wrk );
  // rebuild chains
  mol_wrk = natools.rebuild_chain( xmap, mol_wrk );
  // rebuild bases
  mol_wrk = natools.rebuild_bases( xmap, mol_wrk );

  const clipper::String chainid = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for ( int i = 0; i < std::min(mol_wrk.size(),51); i++ ) {
    mol_wrk[i].set_id( chainid.substr( i, 1 ) );
    for ( int j = 0; j < mol_wrk[i].size(); j++ ) {
      mol_wrk[i][j].set_seqnum( j+1 );
      mol_wrk[i][j].set_type( "  U" );
    }
  }

  // count monomers after
  for ( int c = 0; c < mol_wrk.size(); c++ ) n2 += mol_wrk[c].size();

  // update model
  mmdbfile->export_minimol( mol_wrk );

  return ( n2 > n1 );
}
