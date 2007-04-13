/*! \file buccaneer-prot.cpp buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */

#include "buccaneer-prot.h"

#include <clipper/clipper-contrib.h>


clipper::RTop_orth Ca_group::rtop_from_std_ori() const
{
  const clipper::Coord_orth c1 = coord_c() - coord_ca();
  const clipper::Coord_orth c2 = coord_n() - coord_ca();
  const clipper::Vec3<> v1( (c1.unit()+c2.unit()).unit() );
  const clipper::Vec3<> v2( clipper::Vec3<>::cross(c1,c2).unit() );
  const clipper::Vec3<> v3( clipper::Vec3<>::cross(v1,v2).unit() );
  return
    clipper::RTop_orth( clipper::Mat33<>( v1[0], v2[0], v3[0],
					  v1[1], v2[1], v3[1],
					  v1[2], v2[2], v3[2] ), coord_ca() );
}

clipper::RTop_orth Ca_group::rtop_beta_carbon() const
{
  clipper::RTop_orth rtop = rtop_from_std_ori();
  return clipper::RTop_orth( rtop.rot(), rtop.trn() + coord_cb() - coord_ca() );
}

Ca_group Ca_group::next_ca_group(const clipper::ftype& psi, const clipper::ftype& phi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth n ( coord_n(),  coord_ca(), coord_c(), 1.32, 1.99, psi );
  clipper::Coord_orth ca( coord_ca(), coord_c(),  n,         1.47, 2.15, pi  );
  clipper::Coord_orth c ( coord_c(),  n,          ca,        1.53, 1.92, phi );
  return Ca_group( n, ca, c );
}

Ca_group Ca_group::prev_ca_group(const clipper::ftype& phi, const clipper::ftype& psi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth c ( coord_c(),  coord_ca(), coord_n(), 1.32, 2.15, phi );
  clipper::Coord_orth ca( coord_ca(), coord_n(),  c,         1.53, 1.99, pi  );
  clipper::Coord_orth n ( coord_n(),  c,          ca,        1.47, 1.92, psi );
  return Ca_group( n, ca, c );
}

clipper::ftype Ca_chain::ramachandran_phi( const int& resno ) const
{
  const Ca_chain& chain = (*this);
  if ( resno-1 >= 0 && resno < chain.size() )
    return clipper::Coord_orth::torsion( chain[resno-1].coord_c(), chain[resno].coord_n(), chain[resno].coord_ca(), chain[resno].coord_c() );
  else
    return clipper::Util::nan();
}

clipper::ftype Ca_chain::ramachandran_psi( const int& resno ) const
{
  const Ca_chain& chain = (*this);
  if ( resno >= 0 && resno+1 < chain.size() )
    return clipper::Coord_orth::torsion( chain[resno].coord_n(), chain[resno].coord_ca(), chain[resno].coord_c(), chain[resno+1].coord_n() );
  else
    return clipper::Util::nan();
}


Pr_group::Pr_group( const clipper::Coord_orth& ca, const clipper::Coord_orth& c, const clipper::Coord_orth& other, const TYPE& type )
{
  coord_ca_ = ca;
  coord_c_  = c;
  if ( type == CaCN ) {
    coord_n_ = other;
  } else {
    const clipper::Vec3<> v0 = ( ca - c ).unit();
    const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, other - c ) ) ).unit();
    // length 1.32 angle 1.99
    coord_n_ = c + clipper::Coord_orth( -0.537*v0 + 1.206*v1 );
  }
}

clipper::Coord_orth Pr_group::coord_o() const
{
  const clipper::Vec3<> v0 = ( coord_ca() - coord_c() ).unit();
  const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, coord_n_next() - coord_c() ) ) ).unit();
  // length 1.24 angle 2.11
  return coord_c() + clipper::Coord_orth( -0.637*v0 + 1.064*v1 );
}

clipper::Coord_orth Pr_group::coord_ca_next() const
{
  const clipper::Vec3<> v0 = ( coord_c() - coord_n_next() ).unit();
  const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, coord_ca() - coord_c() ) ) ).unit();
  // length 1.47 angle 2.15
  return coord_n_next() + clipper::Coord_orth( -0.805*v0 + 1.230*v1 );
}

clipper::RTop_orth Pr_group::rtop_from_std_ori() const
{
  clipper::Coord_orth ca = coord_ca_next();
  const clipper::Coord_orth cen = 0.5*( coord_ca() + ca );
  const clipper::Coord_orth c1 = ca - cen;
  const clipper::Coord_orth c2 = coord_c() - cen;
  const clipper::Vec3<> v1( c1.unit() );
  const clipper::Vec3<> v2( clipper::Vec3<>::cross(c1,c2).unit() );
  const clipper::Vec3<> v3( clipper::Vec3<>::cross(v1,v2).unit() );
  return
    clipper::RTop_orth( clipper::Mat33<>( v1[0], v2[0], v3[0],
					  v1[1], v2[1], v3[1],
					  v1[2], v2[2], v3[2] ), cen );
}

Pr_group Pr_group::next_pr_group( const clipper::ftype& phi, const clipper::ftype& psi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth ca(coord_ca(), coord_c(), coord_n_next(), 1.47,2.15,pi);
  clipper::Coord_orth c (coord_c(), coord_n_next(), ca,         1.53,1.92,phi);
  clipper::Coord_orth n (coord_n_next(), ca, c,                 1.32,1.99,psi);
  return Pr_group( ca, c, n, CaCN );
}

Pr_group Pr_group::prev_pr_group( const clipper::ftype& psi, const clipper::ftype& phi ) const
{
  const clipper::ftype pi = clipper::Util::pi();
  clipper::Coord_orth n (coord_n_next(), coord_c(), coord_ca(), 1.47,1.92,psi);
  clipper::Coord_orth c (coord_c(), coord_ca(), n,              1.32,2.15,phi);
  clipper::Coord_orth ca(coord_ca(), n, c,                      1.53,1.99,pi);
  return Pr_group( ca, c, n, CaCN );
}


bool ProteinTools::chain_tidy( clipper::MiniMol& target, const clipper::MiniMol& source )
{
  typedef clipper::MMonomer Mm;
  // create new minimol
  target = clipper::MiniMol( source.spacegroup(), source.cell() );
  // now separate unlinked fragments into separate chains
  clipper::MPolymer mp, mpnull;
  for ( int chn = 0; chn < source.size(); chn++ ) {
    mp = mpnull;
    for ( int res = 0; res < source[chn].size(); res++ ) {
      mp.insert( source[chn][res] );
      if ( res < source[chn].size()-1 )
	if ( !Mm::protein_peptide_bond(source[chn][res],source[chn][res+1]) ) {
	  target.insert( mp );
	  mp = mpnull;
	}
    }
    target.insert( mp );
  }
  // now relabel the chains and residues
  clipper::String labels =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for ( int chn = 0; chn < target.size(); chn++ ) {
    int c1 = chn % labels.length();
    int c2 = chn / labels.length();
    target[chn].set_id( labels.substr( c1, 1 ) );
    for ( int res = 0; res < target[chn].size(); res++ )
      target[chn][res].set_seqnum( 1000*c2 + res + 1 );
  }
  return true;
}

bool ProteinTools::copy_residue_types( clipper::MiniMol& target, const clipper::MiniMol& source )
{
  clipper::Coord_frac ca1, ca2;
  // loop over source model
  for ( int c1 = 0; c1 < source.size(); c1++ )
    for ( int r1 = 0; r1 < source[c1].size(); r1++ ) {
      // residues with Ca and known type
      int a1 = source[c1][r1].lookup( " CA ", clipper::MM::ANY );
      if ( a1 >= 0 && source[c1][r1].type() != "UNK" )
	for ( int c2 = 0; c2 < target.size(); c2++ )
	  for ( int r2 = 0; r2 < target[c2].size(); r2++ ) {
	    // residues with Ca and unknown type
	    int a2 = target[c2][r2].lookup( " CA ", clipper::MM::ANY );
	    if ( a2 >= 0 && target[c2][r2].type() == "UNK" ) {
	      ca1 = source[c1][r1][a1].coord_orth().coord_frac( target.cell() );
	      ca2 = target[c2][r2][a2].coord_orth().coord_frac( target.cell() );
	      ca2 = ca2.symmetry_copy_near( target.spacegroup(), target.cell(),
					    ca1 );
	      // if target is near source, set its residue type to match
	      if ( (ca2-ca1).lengthsq( target.cell() ) < 1.0 )
		target[c2][r2].set_type( source[c1][r1].type() );
	    }
	  }
    }
  return true;
}


bool ProteinTools::globularise( clipper::MiniMol& mol, const clipper::Coord_frac cent )
{
  typedef clipper::MMonomer Mm;
  const clipper::Spacegroup& spgr = mol.spacegroup();
  const clipper::Cell&       cell = mol.cell();
  clipper::Coord_frac cf1, cf2, cf3, cfmin;

  for ( int chn = 0; chn < mol.size(); chn++ ) {
    int res0 = 0;
    for ( int res = 0; res < mol[chn].size(); res++ ) {
      bool chnbrk;
      if ( res == mol[chn].size()-1 )
	chnbrk = true;
      else
	chnbrk = !Mm::protein_peptide_bond( mol[chn][res], mol[chn][res+1] );
      if ( chnbrk ) {
	// process res0...res
	// get COM of current chain
	clipper::Coord_orth co( 0.0, 0.0, 0.0 );
	double no = 0.0;
	for ( int r1 = res0; r1 <= res; r1++ )
	  for ( int a1 = 0; a1 < mol[chn][r1].size(); a1++ ) {
	    co += mol[chn][r1][a1].coord_orth();
	    no += 1.0;
	  }
	co = (1.0/no) * co;
	cf1 = co.coord_frac( cell );
	// find symop and lattice shift which brings it close to cent
	double r2min = 1.0e9;
	int symin = 0;
	cfmin = clipper::Coord_frac( 0.0, 0.0, 0.0 );
	for ( int s = 0; s < spgr.num_symops(); s++ ) {
	  cf2 = spgr.symop(s) * cf1;
	  cf3 = cf2.lattice_copy_near( cent );
	  double r2 = (cf3-cent).lengthsq(cell);
	  if ( r2 < r2min ) {
	    r2min = r2;
	    symin = s;
	    cfmin = cf3 - cf2;
	  }
	}
	// apply the shifts
	for ( int r1 = res0; r1 <= res; r1++ )
	  for ( int a1 = 0; a1 < mol[chn][r1].size(); a1++ ) {
	    cf1 = mol[chn][r1][a1].coord_orth().coord_frac( cell );
	    cf1 = spgr.symop(symin) * cf1 + cfmin;
	    mol[chn][r1][a1].set_coord_orth( cf1.coord_orth( cell ) );
	  }
	// done processing
	res0 = res+1;
      }
    }
  }
  return true;
}

bool ProteinTools::globularise( clipper::MiniMol& mol )
{
  // iteratively globularise the model
  for ( int i = 0; i < 3; i++ ) {
    // find centre of mass
    clipper::Coord_orth co( 0.0, 0.0, 0.0 );
    double no = 0.0;
    for ( int chn = 0; chn < mol.size(); chn++ )
      for ( int res = 0; res < mol[chn].size(); res++ )
	for ( int atm = 0; atm < mol[chn][res].size(); atm++ ) {
	  co += mol[chn][res][atm].coord_orth();
	  no += 1.0;
	}
    co = (1.0/no) * co;
    clipper::Coord_frac cf = co.coord_frac( mol.cell() );
    // globularise about it
    globularise( mol, cf );
  }
  return true;
}

const int ProteinTools::ntype = 21;
const char ProteinTools::rtype1[21]    = {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',  'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',  'M'};
const char ProteinTools::rtype3[21][4] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE"};

int ProteinTools::residue_index( clipper::String code, bool translate )
{
  int r = -1;
  if ( code.length() == 3 )
    for ( r = 0; r < ntype; r++ )
      if ( strncmp( code.c_str(), rtype3[r], 3 ) == 0 ) break;
  if ( code.length() == 1 )
    for ( r = 0; r < ntype; r++ )
      if ( code[0] == rtype1[r] ) break;
  if ( r == ntype ) return -1;
  if ( translate ) {
    if ( r == 20 ) r = 12;
  }
  return r;
}

clipper::String ProteinTools::residue_code_1( int index )
{
  if ( index < 0 ) return "";
  return clipper::String( rtype1[index] );
}

clipper::String ProteinTools::residue_code_3( int index )
{
  if ( index < 0 ) return "";
  return clipper::String( rtype3[index] );
}

clipper::String ProteinTools::residue_code( clipper::String code, bool translate )
{
  int r = residue_index( code, translate );
  if ( code.length() == 1 )
    return residue_code_3(r);
  if ( code.length() == 3 )
    return residue_code_1(r);
  return "";
}
