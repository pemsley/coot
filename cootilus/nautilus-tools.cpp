/*! Nautilus tools main chain database */
/* (C) 2011 Kevin Cowtan & University of York all rights reserved */


#include "nautilus-tools.h"

#include <algorithm>


int NucleicAcidTools::bindex[256], NucleicAcidTools::bindext[256];


NucleicAcidTools::NucleicAcidTools()
{
  for ( int i = 0; i < 256; i++ ) {
    char c = char(i);
    int t = -1;
    if ( c == 'A' ) t = 0;
    if ( c == 'C' ) t = 1;
    if ( c == 'G' ) t = 2;
    if ( c == 'T' ) t = 3;
    if ( c == 'U' ) t = 4;
    bindex[i]  = t;
    bindext[i] = std::min(t,3);
  }
}


// static functions

clipper::MiniMol NucleicAcidTools::flag_chains( const clipper::MiniMol& mol )
{
  // flag any chain containing at least one non-NA
  clipper::MiniMol mol_new = mol;
  for ( int c = 0; c < mol_new.size(); c++ ) {
    bool flg = false;
    for ( int r = 0; r < mol_new[c].size(); r++ ) {
      NucleicAcidDB::NucleicAcid na( mol_new[c][r] );
      if ( na.flag() == NucleicAcidDB::NucleicAcid::NONE ) {
	flg = true;
	break;
      }
    }
    if ( flg )
      mol_new[c].set_property( "NON-NA", clipper::Property<bool>( true ) );
  }
  // for the remaining chains, adjust names to a single character
  for ( int c = 0; c < mol_new.size(); c++ ) {
    if ( !mol_new[c].exists_property( "NON-NA" ) ) {
      std::vector<int> flag( mol_new[c].size(), 1 );  // flag NAs as good
      for ( int r = 0; r < mol_new[c].size(); r++ ) {
	clipper::String type = mol_new[c][r].type().trim();
	while ( type.length() > 1 ) type = type.substr(1);
	mol_new[c][r].set_type( type );
	if ( type == "U" ) {  // detect unknown NAs
	  if ( mol_new[c][r].lookup( " O4 ", clipper::MM::ANY ) < 0 ) {
	    if ( mol_new[c][r].lookup( " C4 ", clipper::MM::ANY ) < 0 )
	      flag[r] =  0;  // not base atoms: maybe unknown
	    else
	      flag[r] = -1;  // C4 but no O4: definitely unknown
	  }
	}
      }
      if ( flag[0] == 0 )             flag[0] = -1;
      if ( flag[flag.size()-1] == 0 ) flag[flag.size()-1] = -1;
      for ( int r = 1; r < flag.size()-1; r++ )
        if ( flag[r] == 0 && flag[r-1] == -1 ) flag[r] = -1;
      for ( int r = flag.size()-2; r > 0; r-- )
        if ( flag[r] == 0 && flag[r+1] == -1 ) flag[r] = -1;
      for ( int r = 0; r < mol_new[c].size(); r++ )
	if ( flag[r] == -1 ) mol_new[c][r].set_type( "?" );
    }
  }
  return mol_new;
}


clipper::RTop_orth NucleicAcidTools::symmetry_rtop( const std::vector<clipper::Coord_orth>& cowrk, clipper::Coord_orth& coref, const clipper::Spacegroup& spgr, const clipper::Cell& cell )
{
  std::vector<clipper::Coord_frac> cwrk( cowrk.size() );
  for ( int a = 0; a < cowrk.size(); a++ )
    cwrk[a] = cowrk[a].coord_frac(cell);
  clipper::Coord_frac cref = coref.coord_frac(cell);
  clipper::Coord_frac c1, c2;
  double d2, d2min(1.0e12);
  int smin(0);
  clipper::Coord_frac dmin(0.0,0.0,0.0);
  for ( int s = 0; s < spgr.num_symops(); s++ )
    for ( int a = 0; a < cwrk.size(); a++ ) {
      c1 = ( spgr.symop(s) * cwrk[a] );
      c2 = c1.lattice_copy_near( cref );
      d2 = ( c2 - cref ).lengthsq( cell );
      if ( d2 < d2min ) {
	d2min = d2;
	smin = s;
	dmin = c2 - c1;
      }
    }
  clipper::RTop_frac rf( spgr.symop(smin).rot(), spgr.symop(smin).trn()+dmin );
  return rf.rtop_orth( cell );
}


clipper::MiniMol NucleicAcidTools::chain_sort( const clipper::MiniMol& mol )
{
  std::vector<std::pair<int,int> > chnsiz( mol.size() ); 
  for ( int chn = 0; chn < mol.size(); chn++ )
    chnsiz[chn] = std::pair<int,int>( -mol[chn].size(), chn );
  std::sort( chnsiz.begin(), chnsiz.end() );
  clipper::MiniMol molnew( mol.spacegroup(), mol.cell() );
  for ( int chn = 0; chn < mol.size(); chn++ )
    molnew.insert( mol[chnsiz[chn].second] );
  return molnew;
}


clipper::Coord_orth NucleicAcidTools::coord_adjust( const clipper::Coord_orth& co, const clipper::Coord_orth& cc3, const clipper::Coord_orth& cf3, const clipper::Coord_orth& cc4, const clipper::Coord_orth& cf4, double rad )
{
  if ( co.is_null() ) return co;
  clipper::Coord_orth result = co;
  double w3 = 1.0 - sqrt( ( co - cf3 ).lengthsq() ) / rad;
  double w4 = 1.0 - sqrt( ( co - cf4 ).lengthsq() ) / rad;
  if ( w3 > 0.0 ) result += w3 * ( cc3 - cf3 );
  if ( w4 > 0.0 ) result += w4 * ( cc4 - cf4 );
  return result;
}


bool NucleicAcidTools::symm_match( clipper::MiniMol& molwrk, const clipper::MiniMol& molref )
{
  clipper::Spacegroup spg1 = clipper::Spacegroup(clipper::Spacegroup::P1);
  clipper::Spacegroup spgr = molwrk.spacegroup();
  clipper::Cell       cell = molwrk.cell();

  // calculate extent of model
  clipper::Atom_list atomr = molref.atom_list();
  clipper::Range<clipper::ftype> urange, vrange, wrange;
  clipper::Coord_frac cfr( 0.0, 0.0, 0.0 );
  for ( int i = 0; i < atomr.size(); i++ ) {
    clipper::Coord_frac cf = atomr[i].coord_orth().coord_frac( cell );
    cfr += cf;
    urange.include( cf.u() );
    vrange.include( cf.v() );
    wrange.include( cf.w() );
  }
  clipper::Coord_frac cf0( urange.min(), vrange.min(), wrange.min() );
  clipper::Coord_frac cf1( urange.max(), vrange.max(), wrange.max() );
  cfr = (1.0/double(atomr.size())) * cfr;

  // calculate mask using wrk cell and ref atoms
  clipper::Resolution reso( 5.0 );
  clipper::Grid_sampling grid( spg1, cell, reso );
  clipper::Grid_range    grng( grid,  cf0,  cf1 );
  grng.add_border(4);
  clipper::NXmap<float> nxmap( cell, grid, grng ), nxflt( cell, grid, grng );
  clipper::EDcalc_mask<float> maskcalc( 2.0 );
  nxmap = 0.0;
  maskcalc( nxmap, atomr );
  MapFilterFn_g5 fn;
  clipper::MapFilter_fft<float>
    fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( nxflt, nxmap );

  // now score each chain, symmetry and offset in turn
  for ( int c = 0; c < molwrk.size(); c++ ) {
    double              bestscr = 0.0;
    int                 bestsym = 0;
    clipper::Coord_frac bestoff( 0.0, 0.0, 0.0 );
    const clipper::Coord_frac cfh( 0.5, 0.5, 0.5 );
    for ( int sym = 0; sym < spgr.num_symops(); sym++ ) {
      clipper::Atom_list atomw = molwrk[c].atom_list();
      clipper::RTop_orth rtop = spgr.symop(sym).rtop_orth( cell );
      clipper::Coord_orth cow( 0.0, 0.0, 0.0 );
      for ( int a = 0; a < atomw.size(); a++ ) {
	atomw[a].transform( rtop );
	cow += atomw[a].coord_orth();
      }
      if ( atomw.size() > 0 ) cow = (1.0/double(atomw.size())) * cow;
      clipper::Coord_frac cfw = cow.coord_frac( cell );
      clipper::Coord_frac cfwt = cfw.lattice_copy_near( cfr - cfh );
      clipper::Coord_frac off0 = cfwt - cfw;

      // try offsets
      for ( double du = 0.0; du <= 1.01; du += 1.0 )
	for ( double dv = 0.0; dv < 1.01; dv += 1.0 )
	  for ( double dw = 0.0; dw < 1.01; dw += 1.0 ) {
	    clipper::Coord_frac off( rint( off0.u() ) + du,
				     rint( off0.v() ) + dv,
				     rint( off0.w() ) + dw );
	    clipper::Coord_orth ofo = off.coord_orth( cell );
	    double scr = 0.0;
	    for ( int a = 0; a < atomw.size(); a++ ) {
	      clipper::Coord_orth coa = atomw[a].coord_orth() + ofo;
	      clipper::Coord_grid cga = nxflt.coord_map( coa ).coord_grid();
	      if ( nxflt.in_map( cga ) ) scr += nxflt.get_data( cga );
	    }
	    if ( scr > bestscr ) {
	      bestscr = scr;
	      bestsym = sym;
	      bestoff = off;
	    }
	  }
    }
    // now transform using the best operator
    clipper::Coord_orth cot = bestoff.coord_orth( cell );
    clipper::RTop_orth rtop = spgr.symop(bestsym).rtop_orth( cell );
    rtop = clipper::RTop_orth( rtop.rot(), rtop.trn()+cot );
    molwrk[c].transform( rtop );
  }

  return true;
}
