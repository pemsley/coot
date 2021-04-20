/*! \file buccaneer-prot.cpp buccaneer library */
/* (C) 2002-2006 University of York all rights reserved */

#include <algorithm>

#include <string.h> // fro strncmp
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
  if ( resno-1 >= 0 && resno < int(chain.size()) )
    return clipper::Coord_orth::torsion( chain[resno-1].coord_c(), chain[resno].coord_n(), chain[resno].coord_ca(), chain[resno].coord_c() );
  else
    return clipper::Util::nan();
}

clipper::ftype Ca_chain::ramachandran_psi( const int& resno ) const
{
  const Ca_chain& chain = (*this);
  if ( resno >= 0 && resno+1 < int(chain.size()) )
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


// Protein loops

ProteinLoop::ProteinLoop( int torsion_sampling )
{
  rama = clipper::Ramachandran( clipper::Ramachandran::All );
  ntor = torsion_sampling;
}


clipper::Coord_orth ProteinLoop::Coord_O( const clipper::Coord_orth ca0, const clipper::Coord_orth c0, const clipper::Coord_orth n1 ) const
{
  clipper::Vec3<> v = ( (c0-n1).unit() + (c0-ca0).unit() ).unit();
  return clipper::Coord_orth( c0 + 1.23*v );
}


clipper::Coord_orth ProteinLoop::Coord_Cb( const clipper::Coord_orth n0, const clipper::Coord_orth ca0, const clipper::Coord_orth c0 ) const
{
  clipper::Vec3<> v1 = ( (ca0-n0).unit() + (ca0-c0).unit() ).unit();
  clipper::Vec3<> v2 = clipper::Coord_orth::cross( (ca0-n0), (ca0-c0) );
  return clipper::Coord_orth( ca0 + 1.04*v1 + 0.53*v2 );
}


/* torsion-rebuild 5 atoms */
std::vector<ProteinLoop::CoordList<5> > ProteinLoop::rebuild5atoms( const clipper::Coord_orth c0, const clipper::Coord_orth n1, const clipper::Coord_orth ca1, const clipper::Coord_orth ca3, const clipper::Coord_orth c3, const clipper::Coord_orth n4 ) const
{
  using clipper::Coord_orth;
  const double dtor = clipper::Util::twopi() / double( ntor );
  std::vector<CoordList<5> > result;
  CoordList<5> r;
  for ( int iphi1 = 0; iphi1 < ntor; iphi1++ ) {
    // build c1
    double phi1 = dtor * double( iphi1 );
    r[0] = Coord_orth( c0, n1, ca1, 1.53, 1.92, phi1 );
    // find 2 possible positions for ca2
    std::vector<Coord_orth> ca2s = constrained_coords( ca1, r[0]-ca1, 3.8, 0.360, ca3, 3.8 );
    for ( unsigned int ipsi1 = 0; ipsi1 < ca2s.size(); ipsi1++ ) {
      r[2] = ca2s[ipsi1];
      double psi1 = Coord_orth::torsion( n1, ca1, r[0], r[2] );
      if ( rama.allowed( phi1, psi1 ) ) {
	// build n2
	r[1] = Coord_orth( n1, ca1, r[0], 1.33, 1.99, psi1 );
	// find 2 possible positions for c2
	std::vector<Coord_orth> cc2s = constrained_coords( r[2], r[2]-r[1], 1.53, 1.22, ca3, 2.43 );
	for ( unsigned int iphi2 = 0; iphi2 < cc2s.size(); iphi2++ ) {
	  // build c2, n3
	  r[3] = cc2s[iphi2];
	  double phi2 = Coord_orth::torsion( r[0], r[1], r[2], r[3] );
	  double psi2 = Coord_orth::torsion( r[1], r[2], r[3], ca3 );
	  if ( rama.allowed( phi2, psi2 ) ) {
	    r[4] = Coord_orth( r[1], r[2], r[3], 1.33, 1.99, psi2 );
	    double ang = Coord_orth::angle( r[4], ca3, c3 );
	    if ( ang > 1.75 && ang < 2.10 ) result.push_back( r );
	  }
	}
      }
    }
  }
  return result;
}


/* torsion-rebuild 8 atoms */
std::vector<ProteinLoop::CoordList<8> > ProteinLoop::rebuild8atoms( const clipper::Coord_orth c0, const clipper::Coord_orth n1, const clipper::Coord_orth ca1, const clipper::Coord_orth ca4, const clipper::Coord_orth c4, const clipper::Coord_orth n5 ) const
{
  using clipper::Coord_orth;
  const double dtor = clipper::Util::twopi() / double( ntor );
  const double pi = clipper::Util::pi();
  std::vector<CoordList<8> > result;
  CoordList<8> r;
  std::vector<CoordList<5> > r5;
  for ( int iphi1 = 0; iphi1 < ntor; iphi1++ ) {
    // build c1
    double phi1 = dtor * double( iphi1 );
    for ( int ipsi1 = 0; ipsi1 < ntor; ipsi1++ ) {
      double psi1 = dtor * double( ipsi1 );
      if ( rama.allowed( phi1, psi1 ) ) {
	r[0] = Coord_orth(  c0,   n1,  ca1, 1.53, 1.92, phi1 );
	r[1] = Coord_orth(  n1,  ca1, r[0], 1.33, 1.99, psi1 );
	r[2] = Coord_orth( ca1, r[0], r[1], 1.47, 2.15, pi   );
	r5 = rebuild5atoms( r[0], r[1], r[2], ca4, c4, n5 );
	for ( unsigned int i = 0; i < r5.size(); i++ ) {
	  for ( int j = 0; j < 5; j++ ) r[j+3] = r5[i][j];
	  result.push_back(r);
	}
      }
    }
  }
  return result;
}


/* find up to 2 possible positions obtained by rotating an arm of
   given base, length and angle to the axis about a rotation axis, to
   give a point a given distance from a target point */
std::vector<clipper::Coord_orth> ProteinLoop::constrained_coords( const clipper::Coord_orth& srcpos, const clipper::Coord_orth& rtnvec, const double& length, const double& angle, const clipper::Coord_orth& tgtpos, const double& tgtdst ) const
{
  using clipper::Coord_orth;
  std::vector<Coord_orth> result;
  Coord_orth v0( rtnvec.unit() );
  Coord_orth v1( Coord_orth::cross( v0, tgtpos-srcpos ).unit() );
  Coord_orth v2( Coord_orth::cross( v1, v0 ).unit() );
  double lcos = length * cos(angle);
  double lsin = length * sin(angle);
  Coord_orth p0 = lcos*v0 + srcpos;
  double z = fabs( (tgtpos-p0) * v0 );
  double x = fabs( (tgtpos-p0) * v2 );
  if ( tgtdst <= z ) return result;
  double d2 = tgtdst*tgtdst - z*z;
  double cost = (lsin*lsin + x*x - d2)/(2.0*lsin*x);
  if ( cost*cost <= 0.995 ) {
    double sint = sqrt(1.0-cost*cost);
    result.push_back( p0 + (lsin*cost)*v2 - (lsin*sint)*v1 );
    result.push_back( p0 + (lsin*cost)*v2 + (lsin*sint)*v1 );
  } else if ( cost*cost <= 1.1 ) {
    result.push_back( p0 + (lsin*cost)*v2 );
  }
  return result;
}


// Protein tools

clipper::String ProteinTools::chain_sequence( const clipper::MPolymer& mp )
{
  clipper::String seq = "";
  for ( int res = 0; res < mp.size(); res++ ) {
    if        ( mp[res].type() == "+++" ) {
      seq += "+";
    } else if ( mp[res].type() == "---" ) {
      seq += "-";
    } else {
      int t = ProteinTools::residue_index( mp[res].type() );
      if ( t >= 0 ) seq += ProteinTools::residue_code_1( t );
      else          seq += "?";
    }
  }
  return seq;
}


std::pair<int,int> ProteinTools::chain_sequence_match( const clipper::String& chnseq, const clipper::MMoleculeSequence& seq )
{
  // convert sequences to unique strings
  std::vector<clipper::String> seqs( seq.size() );
  for ( int chn = 0; chn < seq.size(); chn++ ) {
    clipper::String s = "";
    for ( unsigned int res = 0; res < seq[chn].sequence().length(); res++ )
      s += residue_code_1( residue_index( seq[chn].sequence().substr(res,1) ) );
    seqs[chn] = s;
  }
  // set minimum match threshold
  int minscr = 0;
  for ( unsigned int i = 0; i < chnseq.size(); i++ )
    if ( isupper(chnseq[i]) ) minscr++;
  minscr = minscr/3+4;
  // now find best match
  int bestchn = -1;
  int bestoff = -1;
  int bestscr = minscr;
  int lenc = chnseq.length();
  for ( unsigned int chn = 0; chn < seqs.size(); chn++ ) {
    int lens = seqs[chn].length();
    for ( int off = -lenc+bestscr; off < lens-bestscr; off++ ) {
      int scr = 0;
      int sl = seqs[chn].length();
      for ( int i = 0; i < sl; i++ )
	if ( i-off >= 0 && i-off < int(chnseq.length()) )
	  if ( seqs[chn][i] == chnseq[i-off] )
	    if ( isupper(chnseq[i-off]) )
	      scr++;
      if ( scr > bestscr ) {
	bestchn = chn;
	bestoff = off;
	bestscr = scr;
      }
    }
  }
  // if ( bestchn >= 0 ) { clipper::String tmp = "???????????????????????????????????????????????????????????????????????????????????????????????????"+seqs[bestchn]+"???????????????????????????????????????????????????????????????????????????????????????????????????"; std::cout << tmp.substr( bestoff+99, chnseq.length() ) << "\n"; }
  // return best match
  return std::pair<int,int>( bestchn, bestoff );
}


bool ProteinTools::chain_renumber( clipper::MPolymer& pol, const clipper::MMoleculeSequence& seq )
{
  std::pair<int,int> match = chain_sequence_match( chain_sequence(pol), seq );
  int chn = match.first;
  int off = match.second;
  if ( chn >= 0 ) {
    for ( int res = 0; res < pol.size(); res++ )
      pol[res].set_seqnum( res + off + 1 );
    return true;
  } else {
    return false;
  }
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
    if ( mp.size() > 0 ) target.insert( mp );
  }

  // sort the chains by size
  std::vector<std::pair<int,int> > chnsiz( target.size() ); 
  for ( int chn = 0; chn < target.size(); chn++ )
    chnsiz[chn] = std::pair<int,int>( -target[chn].size(), chn );
  std::sort( chnsiz.begin(), chnsiz.end() );
  clipper::MiniMol temp = target;
  for ( int chn = 0; chn < target.size(); chn++ )
    target[chn] = temp[chnsiz[chn].second];

  // now relabel the chains and residues
  clipper::String labels1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  clipper::String labels2 = "abcdefghijklmnopqrstuvwxyz";
  for ( int chn = 0; chn < target.size(); chn++ ) {
    if ( chn < int(labels1.length()) ) {
      target[chn].set_id( labels1.substr( chn, 1 ) );
      for ( int res = 0; res < target[chn].size(); res++ )
	target[chn][res].set_seqnum( res + 1 );
    } else {
      int c = chn - labels1.length();
      int c1 = c % labels2.length();
      int c2 = c / labels2.length();
      target[chn].set_id( labels2.substr( c1, 1 ) );
      for ( int res = 0; res < target[chn].size(); res++ )
	target[chn][res].set_seqnum( res + 1 + 1000*c2 );
    }
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


bool ProteinTools::globularise( clipper::MiniMol& mol, const clipper::Coord_frac cent, bool nucleotides )
{

   auto check_for_chain_break = [] (const clipper::MMonomer& m1,
                                    const clipper::MMonomer& m2,
                                    const std::string &at_name_1,
                                    const std::string &at_name_2,
                                    double dist_crit = (1.32+0.2)) {
                                   int c1 = m1.lookup(at_name_1, clipper::MM::ANY);
                                   int n2 = m1.lookup(at_name_2, clipper::MM::ANY);
                                   if (c1 >= 0) {
                                      if (n2 >=0) {
                                         double dd = (m1[c1].coord_orth() - m2[n2].coord_orth()).lengthsq();
                                         if (dd > dist_crit * dist_crit) {
                                            return true;
                                         }
                                      }
                                   }
                                   return false; // I don't know what we are looking at
                                };

   const clipper::Spacegroup& spgr = mol.spacegroup();
   const clipper::Cell&       cell = mol.cell();
   if (cell.is_null()) {
      std::cout << "WARNING:: no cell in globularise() " << std::endl;
      return false;
   }

   clipper::Coord_frac cf1, cf2, cf3, cfmin;
   double dist_crit = 1.32 * 0.2;
   std::string at_name_1 = " C  ";
   std::string at_name_2 = " N  ";
   if (nucleotides) {
      dist_crit = 1.57 + 0.33;
      at_name_1 = " O3'";
      at_name_2 = " P  ";
   }

   for ( int chn = 0; chn < mol.size(); chn++ ) {
      int res0 = 0;
      for ( int res = 0; res < mol[chn].size(); res++ ) {
         bool chnbrk;
         if ( res == mol[chn].size()-1 )
            chnbrk = true;
         else
            chnbrk = check_for_chain_break(mol[chn][res], mol[chn][res+1], at_name_1, at_name_2, dist_crit);
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
            std::cout << "apply the shifts with res " << res << std::endl;
            for ( int r1 = res0; r1 <= res; r1++ ) {
               std::cout << "r1 size " << mol[chn][r1].size() << std::endl;
               for ( int a1 = 0; a1 < mol[chn][r1].size(); a1++ ) {
                  clipper::MAtom mat = mol[chn][r1][a1];
                  std::cout << "   " << mat.id_tidy("a") << std::endl;
                  clipper::Coord_orth co_at = mat.coord_orth();
                  cf1 = co_at.coord_frac( cell );
                  cf1 = spgr.symop(symin) * cf1 + cfmin;
                  mol[chn][r1][a1].set_coord_orth( cf1.coord_orth( cell ) );
               }
               // done processing
               res0 = res+1;
            }
         }
      }
   }
   return true;
}

bool ProteinTools::globularise( clipper::MiniMol& mol, bool nucleotides)
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
    globularise( mol, cf, nucleotides);
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
  clipper::String result;
  if ( index >= 0 && index < ntype ) result += rtype1[index];
  return result;
}

clipper::String ProteinTools::residue_code_3( int index )
{
  clipper::String result;
  if ( index >= 0 && index < ntype ) result = rtype3[index];
  return result;
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
