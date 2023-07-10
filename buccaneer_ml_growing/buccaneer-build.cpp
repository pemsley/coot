/*! \file buccaneer-build.cpp buccaneer library */
/* (C) 2006-2008 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-build.h"

#include <clipper/clipper-contrib.h>

#include <algorithm>


void Ca_build::build_rotate_rotamer( clipper::MMonomer& mm, int nr, int nc )
{
  // if no valid rotamer, truncate
  if ( nr < 0 || nr >= mm.protein_sidechain_number_of_rotamers() ) {
    clipper::MMonomer mm1;
    mm1.copy( mm, clipper::MM::COPY_MP );
    for ( int atm = 0; atm < mm.size(); atm++ ) {
      if ( mm[atm].name() == " N  " || mm[atm].name() == " CA " ||
	   mm[atm].name() == " C  " || mm[atm].name() == " O  " ||
	   ( mm.type() != "GLY" && mm[atm].name() == " CB " ) )
	mm1.insert( mm[atm] );
    }
    mm = mm1;
    return;
  }

  // build rotamer
  mm.protein_sidechain_build_rotamer( nr );

  // rotate rotamer
  if ( nc == 0 ) return;
  int a = mm.lookup( " CA ", clipper::MM::ANY );
  int b = mm.lookup( " CB ", clipper::MM::ANY );
  int g = mm.lookup( " CG ", clipper::MM::ANY );
  if ( a < 0 || b < 0 || g < 0 ) return;
  const double dchis[25][2] =
    { { 0.0, 0.0},{-1.0, 0.0},{ 1.0, 0.0},{-2.0, 0.0},{ 2.0, 0.0},
      { 0.0,-1.0},{-1.0,-1.0},{ 1.0,-1.0},{-2.0,-1.0},{ 2.0,-1.0},
      { 0.0, 1.0},{-1.0, 1.0},{ 1.0, 1.0},{-2.0, 1.0},{ 2.0, 1.0},
      { 0.0,-2.0},{-1.0,-2.0},{ 1.0,-2.0},{-2.0,-2.0},{ 2.0,-2.0},
      { 0.0, 2.0},{-1.0, 2.0},{ 1.0, 2.0},{-2.0, 2.0},{ 2.0, 2.0} };
  const double dchi = clipper::Util::d2rad(9.0);
  const double c1 = cos(dchi*dchis[nc][0]);
  const double s1 = sin(dchi*dchis[nc][0]);
  const double c2 = cos(dchi*dchis[nc][1]);
  const double s2 = sin(dchi*dchis[nc][1]);
  const clipper::Coord_orth ca = mm[a].coord_orth();
  const clipper::Coord_orth cb = mm[b].coord_orth();
  const clipper::Coord_orth cg = mm[g].coord_orth();
  const clipper::Vec3<> axis1 = (cb-ca).unit();
  const clipper::Vec3<> axis2 = clipper::Vec3<>::cross(cg-cb,axis1).unit();
  const clipper::Rotation q1( c1, s1*axis1[0], s1*axis1[1], s1*axis1[2] );
  const clipper::Mat33<> mat1 = q1.matrix();
  const clipper::Rotation q2( c2, s2*axis2[0], s2*axis2[1], s2*axis2[2] );
  const clipper::Mat33<> mat2 = q2.matrix();
  const clipper::Mat33<> mat = mat1 * mat2;
  const clipper::RTop_orth rtop( mat, cb - mat*cb );
  for ( int atm = 0; atm < mm.size(); atm++ )
    if ( mm[atm].name() != " N  " && mm[atm].name() != " CA " &&
	 mm[atm].name() != " C  " && mm[atm].name() != " O  " &&
	 mm[atm].name() != " CB " )
      mm[atm].set_coord_orth( rtop * mm[atm].coord_orth() );
}


std::vector<std::pair<double,std::pair<int,int> > > Ca_build::score_rotamers( const clipper::MMonomer& mm, const clipper::Xmap<float>& xmap, const clipper::Map_stats& xstat, int nconf )
{
  std::vector<std::pair<double,std::pair<int,int> > > result;
  int nr = mm.protein_sidechain_number_of_rotamers();
  if ( nr <= 0 ) {
    result.push_back( std::pair<double,std::pair<int,int> >
		      ( 0.0, std::pair<int,int>( 0, 0 ) ) );
    return result;
  }

  // make a list of atom linkage distances
  clipper::MMonomer mr = mm;
  int index = ProteinTools::residue_index_3( mr.type() );
  mr.set_type( ProteinTools::residue_code_3( index ) );
  mr.protein_sidechain_build_rotamer( 0 );
  std::vector<int> atyp( mr.size(), -1 );
  int maxtype = -1;
  for ( int atm = 0; atm < mr.size(); atm++ ) {
    const char& c = mr[atm].id()[2];
    if      ( c == 'G' ) atyp[atm] = 0;
    else if ( c == 'D' ) atyp[atm] = 1;
    else if ( c == 'E' ) atyp[atm] = 2;
    else if ( c == 'Z' ) atyp[atm] = 3;
    else if ( c == 'H' ) atyp[atm] = 4;
    if ( atyp[atm] > maxtype ) maxtype = atyp[atm];
  }
  int nc = 1;
  if      ( maxtype == 1 ) nc = 1;
  else if ( maxtype == 2 ) nc = 5;
  else if ( maxtype >= 3 ) nc = 15;
  if ( nc > nconf ) nc = nconf;

  // calculate density score, discarding Eta atoms
  double zwt = 2.0;  // EXPECTED Z-DIFF BETWEEN CORRECT AND RANDOM SCORES
  double wgts[] = { 1.0, 1.0, 1.0, 1.0, 0.0 };
  for ( int r = 0; r < nr; r++ ) {
    double p = mr.protein_sidechain_build_rotamer( 0 );
    for ( int c = 0; c < nc; c++ ) {
      build_rotate_rotamer( mr, r, c );
      double s = 0.0;
      double w = 0.0;
      for ( int atm = 0; atm < mr.size(); atm++ ) {
	int a = atyp[atm];
	if ( a >= 0 ) {
	  double z = ( xmap.interp<clipper::Interp_cubic>( mr[atm].coord_orth().coord_frac( xmap.cell() ) ) - xstat.mean() ) / xstat.std_dev();
	  s += wgts[a] * z;
	  w += wgts[a];
	}
      }
      s /= w;
      double z0 = -s;
      double z1 = -log( p );
      s = zwt*z0 + z1;
      std::pair<int,int> rc( r, c );
      result.push_back( std::pair<double,std::pair<int,int> >( s, rc ) );
    }
  }

  return result;
}


bool Ca_build::build( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, clipper::String newrestype, bool flexible )
{
  typedef clipper::MMonomer Mm;
  std::vector<std::pair<double,std::pair<int,int> > > scrs;

  // how much can we flex the rotamer?
  int nconf = flexible ? 15 : 5;

  // get map stats
  clipper::Map_stats xstat( xmap );

  // grow side chains
  for ( int chn = 0; chn < mol.size(); chn++ ) {
    // build sidechains
    for ( int res = 0; res < mol[chn].size(); res++ ) {
      // preserve residue name, but build "UNK" as newrestype
      clipper::String oldrestype = mol[chn][res].type();
      if ( ProteinTools::residue_index_3( oldrestype ) < 0 )
	mol[chn][res].set_type( newrestype );
      // get rotamer scores
      scrs = score_rotamers( mol[chn][res], xmap, xstat, nconf );
      std::sort( scrs.begin(), scrs.end() );
      // build best
      build_rotate_rotamer( mol[chn][res], scrs[0].second.first,
                                           scrs[0].second.second );
      // restore name
      mol[chn][res].set_type( oldrestype );
    }

    // build oxygens
    for ( int res = 0; res < mol[chn].size() - 1; res++ )
      if ( Mm::protein_peptide_bond( mol[chn][res], mol[chn][res+1] ) )
        mol[chn][res].protein_mainchain_build_carbonyl_oxygen(mol[chn][res+1]);
      else
	mol[chn][res].protein_mainchain_build_carbonyl_oxygen();
    mol[chn][mol[chn].size()-1].protein_mainchain_build_carbonyl_oxygen();
  }

  // fix clashes
  fix_clashes( mol, xmap, xstat, newrestype );

  return true;
}


bool Ca_build::operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap ) const
{
  return build( mol, xmap, newrestype_, flexible_ );
}


std::vector<Ca_build::Clash> Ca_build::find_clashes( const clipper::MiniMol& mol, const double& d )
{
  // now search for any clashes
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();
  Clash clash;
  std::vector<Clash> clashes;
  clipper::MAtomNonBond nb( mol, 8.0 );
  std::vector<clipper::MAtomIndexSymmetry> atoms;
  std::vector<clipper::Coord_orth> catoms;
  clipper::Coord_orth o1, o2;
  clipper::Coord_frac f1, f2;
  for ( int p = 0; p < mol.size(); p++ ) {
    for ( int m = 0; m < mol[p].size(); m++ ) {
      int a = mol[p][m].lookup( " CA ", clipper::MM::ANY );
      if ( a >= 0 ) {
	o1 = mol[p][m][a].coord_orth();
	f1 = o1.coord_frac( cell );
	atoms = nb( o1, 8.0 );
	catoms.resize( atoms.size() );
	for ( int i = 0; i < atoms.size(); i++ ) {
	  o2 = mol.atom(atoms[i]).coord_orth();
	  f2 = o2.coord_frac( cell );
	  f2 = spgr.symop(atoms[i].symmetry()) * f2;
	  f2 = f2.lattice_copy_near( f1 );
	  catoms[i] = f2.coord_orth(cell);
	}
	// search over atoms in this residue
	double d2min = 1.0e9;
	int i2min = 0;
	for ( int i1 = 0; i1 < atoms.size(); i1++ )
	  if ( atoms[i1].polymer() == p && atoms[i1].monomer() == m ) {
	    // is native atom is movable?
	    clipper::String id = mol[atoms[i1].polymer()]
	      [atoms[i1].monomer()][atoms[i1].atom()].name();
	    if ( id != " CA " && id != " N  " &&
		 id != " C  " && id != " O  " && id != " CB " ) {
	      // and atoms from elsewhere
	      for ( int i2 = 0; i2 < atoms.size(); i2++ )
		if ( atoms[i2].polymer() != p || atoms[i2].monomer() != m ) {
		  // check for a clash
		  double d2 = (catoms[i1]-catoms[i2]).lengthsq();
		  if ( d2 < d2min ) {
		    d2min = d2;
		    i2min = i2;
		  }
		}
	    }
	  }
	// now check the closest clash
	if ( d2min < d*d ) {
	  clash.p1 = p;
	  clash.m1 = m;
	  clash.p2 = atoms[i2min].polymer();
	  clash.m2 = atoms[i2min].monomer();
	  clashes.push_back( clash );
	}
      }
    }
  }
  return clashes;
}


void Ca_build::fix_clash  ( clipper::MMonomer& m1, clipper::MMonomer& m2, const clipper::Xmap<float>& xmap, const clipper::Map_stats& xstat, const double& d, clipper::String newrestype )
{
  // useful data
  const clipper::Spacegroup& spgr = xmap.spacegroup();
  const clipper::Cell&       cell = xmap.cell();

  // deal with unknow residues
  clipper::String oldrestype1 = m1.type();
  clipper::String oldrestype2 = m2.type();
  if ( ProteinTools::residue_index_3(oldrestype1) < 0 ) m1.set_type(newrestype);
  if ( ProteinTools::residue_index_3(oldrestype2) < 0 ) m2.set_type(newrestype);

  // set up two monomers
  clipper::MMonomer mm1(m1), mm2(m2);
  // move second to be close to first
  clipper::Coord_orth co( 0.0, 0.0, 0.0 );
  for ( int i1 = 0; i1 < mm1.size(); i1++ ) co += mm1[i1].coord_orth();
  co = (1.0/double(mm1.size())) * co;
  clipper::Coord_frac cf1 = co.coord_frac( cell );
  for ( int i2 = 0; i2 < mm2.size(); i2++ ) {
    clipper::Coord_frac cf2 = mm2[i2].coord_orth().coord_frac( cell );
    cf2 = cf2.symmetry_copy_near( spgr, cell, cf1 );
    mm2[i2].set_coord_orth( cf2.coord_orth( cell ) );
  }
  // fetch the rotamer scores
  std::vector<std::pair<double,std::pair<int,int> > > score1, score2;
  score1 = score_rotamers( mm1, xmap, xstat, 5 );
  score2 = score_rotamers( mm2, xmap, xstat, 5 );
  // now search over orientations
  double smin = 1.0e9;
  int r1min, r2min;
  r1min = r2min = -1;
  for ( int r1 = 0; r1 < score1.size(); r1++ ) {
    for ( int r2 = 0; r2 < score2.size(); r2++ ) {
      // if this combination gives a viable score, check for clashes...
      double s = score1[r1].first + score2[r2].first;
      if ( s < smin ) {
	build_rotate_rotamer( mm1, score1[r1].second.first,
			           score1[r1].second.second );
	build_rotate_rotamer( mm2, score2[r2].second.first,
			           score2[r2].second.second );
	double d2min = 1.0e9;
	for ( int a1 = 0; a1 < mm1.size(); a1++ )
	  for ( int a2 = 0; a2 < mm2.size(); a2++ ) {
	    double d2 = (mm1[a1].coord_orth()-mm2[a2].coord_orth()).lengthsq();
	    if ( d2 < d2min ) d2min = d2;
	  }
	if ( d2min < d*d ) s = s + 10.0;
	// keep the best
	if ( s < smin ) {
	  smin = s;
	  r1min = r1;
	  r2min = r2;
	}
      }
    }
  }

  // rebuild
  if ( r1min >= 0 && r2min >= 0 ) {
    build_rotate_rotamer( m1, score1[r1min].second.first,
			      score1[r1min].second.second );
    build_rotate_rotamer( m2, score2[r2min].second.first,
			      score2[r2min].second.second );
  }
  m1.set_type( oldrestype1 );
  m2.set_type( oldrestype2 );
}


void Ca_build::fix_clashes( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const clipper::Map_stats& xstat, clipper::String newrestype )
{
  std::vector<Clash> clashes;
  const double d = 1.25;

  // check for clashes
  clashes = find_clashes( mol, d );

  // and try and fix them
  for ( int i = 0; i < clashes.size(); i++ )
    fix_clash( mol[clashes[i].p1][clashes[i].m1],
	       mol[clashes[i].p2][clashes[i].m2], xmap, xstat, d, newrestype );

  // check for any remaining clashes
  clashes = find_clashes( mol, d );

  // and delete the offending atoms
  for ( int i = 0; i < clashes.size(); i++ ) {
    build_rotate_rotamer( mol[clashes[i].p1][clashes[i].m1], -1, -1 );
    build_rotate_rotamer( mol[clashes[i].p2][clashes[i].m2], -1, -1 );
  }
}
