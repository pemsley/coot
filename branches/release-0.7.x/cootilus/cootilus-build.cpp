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

  // move all chains to view centre
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();
  clipper::Coord_frac cf0, cf1, cf2, df;
  cf0 = centre.coord_frac(cell);
  for ( int c = 0; c < mol_wrk.size(); c++ ) {
    double d2min = 1.0e30;
    int smin = -1;
    for ( int r = 0; r < mol_wrk[c].size(); r++ ) {
      for ( int a = 0; a < mol_wrk[c][r].size(); a++ ) {
	cf1 = mol_wrk[c][r][a].coord_orth().coord_frac(cell);
	for ( int s = 0; s < spgr.num_symops(); s++ ) {
	  cf2 = spgr.symop(s) * cf1;
	  double d2 = (cf2.lattice_copy_near(cf0)-cf0).lengthsq(cell);
	  if ( d2 < d2min ) {
	    d2min = d2;
	    smin = s;
	    df = cf2.lattice_copy_near(cf0) - cf2;
	  }
	}
      }
    }
    if ( smin >= 0 ) {
      clipper::RTop_frac rtf( spgr.symop(smin).rot(),
			      spgr.symop(smin).trn()+df );
      clipper::RTop_orth rto = rtf.rtop_orth( cell );
      mol_wrk[c].transform( rto );
    }
  }


  // label chains
  const clipper::String chainid = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for ( int i = 0; i < std::min(mol_wrk.size(),52); i++ ) {
    mol_wrk[i].set_id( chainid.substr( i, 1 ) );
    for ( int j = 0; j < mol_wrk[i].size(); j++ ) {
      mol_wrk[i][j].set_seqnum( j+1 );
      mol_wrk[i][j].set_type( "U" );
    }
  }

  // count monomers after
  for ( int c = 0; c < mol_wrk.size(); c++ ) n2 += mol_wrk[c].size();

  // update model
  mmdbfile->DeleteAllModels();
  mmdbfile->export_minimol( mol_wrk );

  return ( n2 > n1 );
}
