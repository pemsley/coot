/*! \file buccaneer-known.cpp buccaneer library */
/* (C) 2010 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-known.h"
#include "buccaneer-prot.h"


KnownStructure::KnownStructure( const clipper::MiniMol& mol, const std::vector<std::pair<clipper::String,double> >& ids, double nprad )
{
  // make null radius entries
  clipper::MiniMol tmp = mol;
  double nulrad = -1.0;
  for ( int p = 0; p < tmp.size(); p++ )
    for ( int m = 0; m < tmp[p].size(); m++ )
      for ( int a = 0; a < tmp[p][m].size(); a++ )
        tmp[p][m][a].set_property( "RADI", PROP( nulrad ) );

  // fill in nonprotein atoms
  if ( nprad >= 0.0 ) {
    for ( int p = 0; p < tmp.size(); p++ ) {
      for ( int m = 0; m < tmp[p].size(); m++ ) {
        if ( ! ProteinTools::is_protein( tmp[p][m] ) ) {
          bool sel = false;
          if ( tmp[p][m].size() > 1 ) sel = true;
          if ( tmp[p][m].size() == 1 && tmp[p][m][0].element().trim() != "O" ) sel = true;
          if ( sel ) {
            for ( int a = 0; a < tmp[p][m].size(); a++ ) {
              tmp[p][m][a].delete_property( "RADI" );
              tmp[p][m][a].set_property( "RADI", PROP( nprad ) );
            }
          }
        }
      }
    }
  }

  // now fill in selected values
  radius_max = nprad > 0.0 ? nprad : 0.0;
  for ( int id = 0; id < ids.size(); id++ ) {
    clipper::String sel = ids[id].first;
    double radius = ids[id].second;
    if ( radius > radius_max ) radius_max = radius;
    std::vector<clipper::MAtomIndex> atoms = tmp.select_atom_index( sel, clipper::MM::ANY );
    for ( int a = 0; a < atoms.size(); a++ ) {
      clipper::MAtom& atom = tmp.atom( atoms[a] );
      double r = dynamic_cast<const PROP&>(atom.get_property("RADI")).value();
      if ( radius > r ) {
        atom.delete_property( "RADI" );
        atom.set_property( "RADI", PROP( radius ) );
      }
    }
  }

  // at this point:
  //  all ignored atoms have RADI = -1,
  //  all pass-through atoms will have RADI = 0,
  //  all protect atoms will have RADI > 0
  // now make a new molecule with just the atoms of interest
  //  known contains RADI > 0 atoms
  //  knownall contains RADI >= 0 atoms
  known.init( tmp.spacegroup(), tmp.cell() );
  knownall.init( tmp.spacegroup(), tmp.cell() );
  for ( int p = 0; p < tmp.size(); p++ ) {
    clipper::MPolymer mp, mpall;
    mp.copy( tmp[p], clipper::MM::COPY_M );
    mpall.copy( tmp[p], clipper::MM::COPY_M );
    for ( int m = 0; m < tmp[p].size(); m++ ) {
      clipper::MMonomer mm, mmall;
      mm.copy( tmp[p][m], clipper::MM::COPY_M );
      mmall.copy( tmp[p][m], clipper::MM::COPY_M );
      for ( int a = 0; a < tmp[p][m].size(); a++ ) {
        const clipper::MAtom& atom = tmp[p][m][a];
        double r = dynamic_cast<const PROP&>(atom.get_property("RADI")).value();
        if ( r > 0.0 ) mm.insert( atom );
        if ( r >= 0.0 ) mmall.insert( atom );
      }
      if ( mm.size() > 0 ) mp.insert( mm );
      if ( mmall.size() > 0 ) mpall.insert( mmall );
    }
    if ( mp.size() > 0 ) known.insert( mp );
    if ( mpall.size() > 0 ) knownall.insert( mpall );
  }

  // make non-bond object
  if ( known.size() > 0 ) knownnb = clipper::MAtomNonBond( known, 4.0 );
}


bool KnownStructure::copy_to( clipper::MiniMol& mol, bool includeAll ) const
{
  clipper::MiniMol prior = includeAll ? knownall : known;

  // check for id clash
  bool idclash = false;
  for ( int c1 = 0; c1 < mol.size(); c1++ )
    for ( int c2 = 0; c2 < prior.size(); c2++ )
      if ( mol[c1].id() == prior[c2].id() ) idclash = true;
  if ( idclash )
    for ( int c1 = 0; c1 < mol.size(); c1++ ) mol[c1].set_id( "" );

  // combine prior and built chains
  clipper::MiniMol built = mol;
  ProteinTools::symm_match( prior, built );
  mol = clipper::MiniMol( built.spacegroup(), built.cell() );
  for ( int c = 0; c < prior.size(); c++ ) mol.insert( prior[c] );
  for ( int c = 0; c < built.size(); c++ ) mol.insert( built[c] );  

  // return true if copy was without clash
  return !idclash;
}


bool KnownStructure::clash( const clipper::Coord_orth& coord ) const
{
  // fast path
  if ( known.size() == 0 ) return false;

  // use nnb calc to detect clashes with known structure
  clipper::Coord_frac f1, f2;
  const clipper::Spacegroup& spgr = known.spacegroup();
  const clipper::Cell&       cell = known.cell();

  const std::vector<clipper::MAtomIndexSymmetry> atoms =
    knownnb.atoms_near( coord, radius_max );
  f1 = coord.coord_frac( cell );
  for ( int i = 0; i < atoms.size(); i++ ) {
    const clipper::MAtom& atom = known.atom(atoms[i]);
    f2 = atom.coord_orth().coord_frac(cell);
    f2 = spgr.symop(atoms[i].symmetry()) * f2;
    f2 = f2.lattice_copy_near( f1 );
    const double d2 = ( f2 - f1 ).lengthsq( cell );
    const double r2 = clipper::Util::sqr(
      dynamic_cast<const PROP&>(atom.get_property("RADI")).value());
    if ( d2 < r2 ) return true;
  }

  return false;
}


bool KnownStructure::prune( clipper::MiniMol& mol ) const
{
  // fast path
  if ( known.size() == 0 ) return false;

  clipper::Cell       cell = mol.cell();
  clipper::Spacegroup spgr = mol.spacegroup();

  // loop over chains and remove clashing residues
  clipper::MiniMol moltmp = mol;
  clipper::Coord_frac cf1, cf2;
  for ( int chn = 0; chn < moltmp.size(); chn++ ) {
    for ( int res = 0; res < moltmp[chn].size(); res++ ) {
      for ( int atm = 0; atm < moltmp[chn][res].size(); atm++ ) {
        const clipper::MAtom& atom = moltmp[chn][res][atm];
        if ( atom.id() == " N  " || atom.id() == " CA " ||
             atom.id() == " C  " || atom.id() == " O  " ) {
          if ( clash( atom.coord_orth() ) ) {
            moltmp[chn][res].set_type( "~~~" );
            break;
          }
        }
      }
    }
  }

  // eliminate any sequences of less than 6 residues
  mol = clipper::MiniMol( spgr, cell );
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < moltmp.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < moltmp[chn].size(); res++ ) {
      if ( moltmp[chn][res].type() != "~~~" ) {
        mp.insert( moltmp[chn][res] );
      } else {
        if ( mp.size() > 5 ) mol.insert( mp );
        mp = mpnull;
      }
    }
    if ( mp.size() > 5 ) mol.insert( mp );
  }

  return true;
}


std::pair<clipper::String,double> KnownStructure::parse( clipper::String arg )
{
  std::vector<clipper::String> argsplit = arg.split( ":" );
  if ( argsplit.size() == 1 ) argsplit.push_back( 0.0 );
  return std::pair<clipper::String,double>( argsplit[0], argsplit[1].f() );
}


void KnownStructure::debug() const
{
  for ( int chn = 0; chn < knownall.size(); chn++ ) {
    for ( int res = 0; res < knownall[chn].size(); res++ ) {
      for ( int atm = 0; atm < knownall[chn][res].size(); atm++ ) {
        const double r = dynamic_cast<const PROP&>(knownall[chn][res][atm].get_property("RADI")).value();
        std::cout << "/" << knownall[chn].id() << "/" << knownall[chn][res].id() << "/" << knownall[chn][res][atm].id() << " : " << r << std::endl;
      }
    }
  }
}
