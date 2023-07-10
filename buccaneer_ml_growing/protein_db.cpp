/*! ProteinDB Top 500 main chain database */
/* (C) 2008-2009 Kevin Cowtan & University of York all rights reserved */


#include "protein_db.h"
#include <algorithm>
#include <fstream>


namespace ProteinDB {


const int Residue::TypeMask::msks[] =
  { 0x00000000, 0x00000001, 0x00000000, 0x00000010,   // @ABC
    0x00000008, 0x00000040, 0x00002000, 0x00000080,   // DEFG
    0x00000100, 0x00000200, 0x00000000, 0x00000800,   // HIJK
    0x00000400, 0x00001000, 0x00000004, 0x00000000,   // LMNO
    0x00004000, 0x00000020, 0x00000002, 0x00008000,   // PQRS
    0x00010000, 0x00000000, 0x00080000, 0x00020000,   // TUVW
    0x00000000, 0x00040000, 0x00000000, 0x00000000,   // XYZ[
    0x00000000, 0x00000000, 0x00000000, 0x000fffff }; // \]^?
const char Residue::rtype1[Residue::ntype] =
  {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
     'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
     'M',  '?'};
const char Residue::rtype3[Residue::ntype][4] =
  {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
   "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
   "MSE","UNK"};


Residue::Residue( clipper::Coord_orth& ca, const clipper::String& type )
{
  typ = residue_type( type );
  cax = float( ca.x() );
  cay = float( ca.y() );
  caz = float( ca.z() );
  set_flag(CALPHA);
}


Residue::Residue( clipper::Coord_orth& cn, clipper::Coord_orth& ca, clipper::Coord_orth& cc, const clipper::String& type )
{
  typ = residue_type( type );
  nnx = float( cn.x() );
  nny = float( cn.y() );
  nnz = float( cn.z() );
  cax = float( ca.x() );
  cay = float( ca.y() );
  caz = float( ca.z() );
  ccx = float( cc.x() );
  ccy = float( cc.y() );
  ccz = float( cc.z() );
  set_flag(NORMAL);
}


Residue::Residue( const clipper::MMonomer& mm ) {
  typ = residue_type( mm.type() );
  int cn = mm.lookup( " N  ", clipper::MM::ANY );
  int ca = mm.lookup( " CA ", clipper::MM::ANY );
  int cc = mm.lookup( " C  ", clipper::MM::ANY );
  if ( ca >= 0 && cn >= 0 && cc >= 0 ) {
    nnx = float( mm[cn].coord_orth().x() );
    nny = float( mm[cn].coord_orth().y() );
    nnz = float( mm[cn].coord_orth().z() );
    cax = float( mm[ca].coord_orth().x() );
    cay = float( mm[ca].coord_orth().y() );
    caz = float( mm[ca].coord_orth().z() );
    ccx = float( mm[cc].coord_orth().x() );
    ccy = float( mm[cc].coord_orth().y() );
    ccz = float( mm[cc].coord_orth().z() );
    set_flag(NORMAL);
  } else if ( ca >= 0 ) {
    cax = float( mm[ca].coord_orth().x() );
    cay = float( mm[ca].coord_orth().y() );
    caz = float( mm[ca].coord_orth().z() );
    set_flag(CALPHA);
  } else {
    set_flag(NONE);
  }
}


clipper::MMonomer Residue::mmonomer() const
{
  clipper::MMonomer mm;
  clipper::MAtom ma_n, ma_a, ma_c;
  ma_n = ma_a = ma_c = clipper::MAtom::null();
  ma_n.set_u_iso ( 0.25 ); ma_n.set_occupancy( 1.0 );
  ma_n.set_id( " N  " ); ma_n.set_element( "N" );
  ma_a.set_u_iso ( 0.25 ); ma_a.set_occupancy( 1.0 );
  ma_a.set_id( " CA " ); ma_a.set_element( "C" );
  ma_c.set_u_iso ( 0.25 ); ma_c.set_occupancy( 1.0 );
  ma_c.set_id( " C  " ); ma_c.set_element( "C" );
  ma_n.set_coord_orth( coord_n () );
  ma_a.set_coord_orth( coord_ca() );
  ma_c.set_coord_orth( coord_c () );
  mm.insert( ma_n );
  mm.insert( ma_a );
  mm.insert( ma_c );
  return mm;
}


clipper::Coord_orth Residue::coord_n () const
{ return clipper::Coord_orth( double(nnx), double(nny), double(nnz) ); }


clipper::Coord_orth Residue::coord_ca() const
{ return clipper::Coord_orth( double(cax), double(cay), double(caz) ); }


clipper::Coord_orth Residue::coord_c () const
{ return clipper::Coord_orth( double(ccx), double(ccy), double(ccz) ); }


void Residue::transform( const clipper::RTop_orth& rtop )
{
  clipper::Coord_orth cn = rtop * coord_n ();
  nnx = float( cn.x() );
  nny = float( cn.y() );
  nnz = float( cn.z() );
  clipper::Coord_orth ca = rtop * coord_ca();
  cax = float( ca.x() );
  cay = float( ca.y() );
  caz = float( ca.z() );
  clipper::Coord_orth cc = rtop * coord_c ();
  ccx = float( cc.x() );
  ccy = float( cc.y() );
  ccz = float( cc.z() );
}


bool Residue::merge( const Residue& other, const double wn, const double wa, const double wc )
{
  if ( flag() != NORMAL || other.flag() != NORMAL ) return false;
  nnx = (1.0-wn)*nnx + wn*other.nnx;
  nny = (1.0-wn)*nny + wn*other.nny;
  nnz = (1.0-wn)*nnz + wn*other.nnz;
  cax = (1.0-wa)*cax + wa*other.cax;
  cay = (1.0-wa)*cay + wa*other.cay;
  caz = (1.0-wa)*caz + wa*other.caz;
  ccx = (1.0-wc)*ccx + wc*other.ccx;
  ccy = (1.0-wc)*ccy + wc*other.ccy;
  ccz = (1.0-wc)*ccz + wc*other.ccz;
  return true;
}


void Residue::data_import( const char* d )
{
  unpack_float( d+ 0, nnx );
  unpack_float( d+ 2, nny );
  unpack_float( d+ 4, nnz );
  unpack_float( d+ 6, cax );
  unpack_float( d+ 8, cay );
  unpack_float( d+10, caz );
  unpack_float( d+12, ccx );
  unpack_float( d+14, ccy );
  unpack_float( d+16, ccz );
  typ = d[18];
  flg = d[19];
}


void Residue::data_export( char* d ) const
{
  pack_float( d+ 0, nnx );
  pack_float( d+ 2, nny );
  pack_float( d+ 4, nnz );
  pack_float( d+ 6, cax );
  pack_float( d+ 8, cay );
  pack_float( d+10, caz );
  pack_float( d+12, ccx );
  pack_float( d+14, ccy );
  pack_float( d+16, ccz );
  d[18] = typ;
  d[19] = flg;
}


char Residue::residue_type( const clipper::String& type )
{
  if ( type.length() == 3 ) {
    for ( int t = 0; t < ntype; t++ )
      if ( type[0] == rtype3[t][0] &&
	   type[1] == rtype3[t][1] &&
	   type[2] == rtype3[t][2] ) return rtype1[t];
  } else if ( type.length() == 1 ) {
    for ( int t = 0; t < ntype; t++ )
      if ( type[0] == rtype1[t] ) return rtype1[t];
  }
  return ' ';
}


bool Chain::add_pdb( const clipper::String file )
{
  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  clipper::MMDBfile mfile;
  clipper::MiniMol mol;
  mfile.SetFlag( mmdbflags );
  mfile.read_file( file );
  mfile.import_minimol( mol );
  if ( mol.size() == 0 ) return false;
  for ( int c = 0; c < mol.size(); c++ ) {
    clipper::MPolymer mp;
    // select residues by occupancy
    for ( int r = 0; r < mol[c].size(); r++ ) {
      int a = mol[c][r].lookup( " CA ", clipper::MM::ANY );
      if ( a >= 0 ) {
	if ( mol[c][r][a].occupancy() > 0.01 &&
	     mol[c][r][a].u_iso() < clipper::Util::b2u(40.0) )
	  mp.insert( mol[c][r] );
      }
    }
    // shift centre-of-mass of chain to the origin
    clipper::Coord_orth cm( 0.0, 0.0, 0.0 );
    double              sm = 0.0;
    for ( int r = 0; r < mp.size(); r++ ) {
      int a = mp[r].lookup( " CA ", clipper::MM::ANY );
      cm += mp[r][a].coord_orth();
      sm += 1.0;
    }
    cm = (-1.0/sm) * cm;
    clipper::RTop_orth rt( clipper::Mat33<>::identity(), cm );
    mp.transform( rt );
    // now add the chain to the db
    for ( int r = 0; r < mp.size(); r++ ) {
      Residue rp( mp[r] );
      if ( rp.flag() == Residue::NORMAL && rp.type() != ' ' )
	add_residue( rp );
    }
  }
  return true;
}


bool Chain::save_db( const clipper::String file ) const
{
 std::ofstream fs( file.c_str(), std::ios::out | std::ios::binary );
  char d[20];
  for ( int r = 0; r < dbresidues.size(); r++ ) {
    dbresidues[r].data_export( d );
    fs.write( d, 20 );
  }
  fs.close();
  return true;
}


bool Chain::load_db( const clipper::String file )
{
  dbresidues.clear();
  // read whole file (for speed)
  std::ifstream fs( file.c_str(), std::ios::in | std::ios::binary );
  if ( !fs ) return false;
  fs.seekg( 0, std::ios::end );
  int i2 = fs.tellg();
  fs.seekg( 0, std::ios::beg );
  int i1 = fs.tellg();
  int l = i2 - i1;
  std::vector<char> d(l);
  fs.read( d.data(), l );
  fs.close();
  if ( l%20 != 0 ) return false;
  // import file data
  dbresidues.resize( l/20 );
  for ( int r = 0; r < dbresidues.size(); r++ ) {
    dbresidues[r].data_import( d.data() + 20*r );
  }
  return true;
}


bool Chain::merge( const Chain& other, const std::vector<double>& wgt )
{
  if ( other.size() != size() ) return false;
  if ( wgt.size() != 3*size() ) return false;
  for ( int r = 0; r < dbresidues.size(); r++ )
    dbresidues[r].merge( other.dbresidues[r],
			 wgt[3*r], wgt[3*r+1], wgt[3*r+2] );
  return true;
}


Chain Chain::extract( int offset, int len ) const
{
  Chain dbc;
  for ( int i = 0; i < len; i++ ) dbc.add_residue( dbresidues[offset+i] );
  return dbc;
}


bool Chain::is_continuous() const
{
  // go through and find elements where there is a chain break
  const double dmin = 4.0;
  std::vector<bool> cterm( dbresidues.size(), false );
  for ( int i = 0; i < dbresidues.size()-1; i++ ) {
    int j = i + 1;
    if ( !dbresidues[i].is_null() && !dbresidues[j].is_null() ) {
      clipper::Coord_orth co1 = dbresidues[i].coord_ca();
      clipper::Coord_orth co2 = dbresidues[j].coord_ca();
      const double d2 = ( co1 - co2 ).lengthsq();
      if ( d2 > dmin*dmin ) return false;
    }
  }
  return true;
}


void Chain::transform( const clipper::RTop_orth& rtop )
{
  for ( int r = 0; r < dbresidues.size(); r++ )
    dbresidues[r].transform( rtop );
}


void Chain::lsq_superpose( const Chain& frag )
{
  std::vector<clipper::Coord_orth> c1, c2;
  for ( int i = 0; i < frag.size(); i++ ) {
    if ( !dbresidues[i].is_null() && !frag.dbresidues[i].is_null() ) {
      c1.push_back( dbresidues[i].coord_ca() );
      c2.push_back( frag.dbresidues[i].coord_ca() );
    }
  }
  transform( clipper::RTop_orth( c1, c2 ) );
}


void Chain::lsq_superpose( const Chain& frag, const std::vector<double>& wgts )
{
  std::vector<clipper::Coord_orth> c1, c2;
  std::vector<double> w;
  for ( int i = 0; i < frag.size(); i++ ) {
    if ( !dbresidues[i].is_null() && !frag.dbresidues[i].is_null() ) {
      c1.push_back( dbresidues[i].coord_ca() );
      c2.push_back( frag.dbresidues[i].coord_ca() );
      w.push_back( wgts[i] );
    }
  }
  transform( clipper::RTop_orth( c1, c2, w ) );
}


double Chain::rmsd( const Chain& other ) const
{
  double s0(0.0), s1(0.0);
  for ( int i = 0; i < dbresidues.size(); i++ ) {
    if ( !dbresidues[i].is_null() && !other.dbresidues[i].is_null() ) {
      s0 += 1.0;
      s1 +=
	(dbresidues[i].coord_ca()-other.dbresidues[i].coord_ca()).lengthsq();
    }
  }
  return sqrt( s1 / s0 );
}


double Chain::rmsd( const Chain& other, const std::vector<double>& wgts ) const
{
  double s0(0.0), s1(0.0);
  for ( int i = 0; i < dbresidues.size(); i++ ) {
    if ( !dbresidues[i].is_null() && !other.dbresidues[i].is_null() ) {
      s0 += wgts[i];
      s1 += wgts[i] *
	(dbresidues[i].coord_ca()-other.dbresidues[i].coord_ca()).lengthsq();
    }
  }
  return sqrt( s1 / s0 );
} 


void Chain::debug() const
{
  double x1(0.0),y1(0.0),z1(0.0),x2(0.0),y2(0.0),z2(0.0);
  for ( int r = 0; r < dbresidues.size(); r++ ) {
    x1 = clipper::Util::min( x1, dbresidues[r].coord_ca().x() );
    y1 = clipper::Util::min( y1, dbresidues[r].coord_ca().y() );
    z1 = clipper::Util::min( z1, dbresidues[r].coord_ca().z() );
    x2 = clipper::Util::max( x2, dbresidues[r].coord_ca().x() );
    y2 = clipper::Util::max( y2, dbresidues[r].coord_ca().y() );
    z2 = clipper::Util::max( z2, dbresidues[r].coord_ca().z() );
  }
  std::cout << "DEBUG Nres: " << dbresidues.size() << std::endl;
  std::cout << "DEBUG Cmin: " << x1 << " " << y1 << " " << z1 << std::endl;
  std::cout << "DEBUG Cmax: " << x2 << " " << y2 << " " << z2 << std::endl;
}


// ChainDB methods

void ChainDB::init( const clipper::String file )
{
  load_db( file );
  calc_distances();
}


void ChainDB::calc_distances()
{
  if ( size() == 0 ) return;

  // go through and find elements where there is a chain break
  const double dmin = 4.0;
  std::vector<bool> cterm( dbresidues.size(), false );
  for ( int i = 0; i < dbresidues.size()-1; i++ ) {
    int j = i + 1;
    if ( !dbresidues[i].is_null() && !dbresidues[j].is_null() ) {
      clipper::Coord_orth co1 = dbresidues[i].coord_ca();
      clipper::Coord_orth co2 = dbresidues[j].coord_ca();
      const double d2 = ( co1 - co2 ).lengthsq();
      if ( d2 > dmin*dmin ) cterm[i] = true;
    }
  }
  cterm.back() = true;

  // fill out the distance matrix table
  dbdistvecs.resize( dbresidues.size() );
  for ( int i = 0; i < dbdistvecs.size(); i++ ) {
    const clipper::Coord_orth co1 = dbresidues[i].coord_ca();
    for ( int d = 0; d < ndist; d++ ) dbdistvecs[i].data[d] = -1.0;
    for ( int d = 0; d < ndist; d++ ) {
      int j = i + d + 1;
      //if ( j == dbresidues.size() ) break;
      if ( cterm[j-1] ) break;
      if ( !dbresidues[i].is_null() && !dbresidues[j].is_null() ) {
	const clipper::Coord_orth co2 = dbresidues[j].coord_ca();
	double dist = sqrt( (co2-co1).lengthsq() );
	dbdistvecs[i].data[d] = float( dist );
      }
    }
  }
}


/*!
  Return a fast distance score of a fragment against the nth fragment
  in the DB. If the nth fragment in the DB is broken across a chain
  boundary, then the return value is -1.0;

  The score is roughly comparable to a sum of squared differences
  between the search fragment and the DB fragment over all atom
  distances present in the fragment. The number of distances is equal
  to the number of pairs of non-null C-alphas in the searh fragment.

  \param frag The data for the search fragment
  \param offset The fragment in the DB against which to score
  \return The fragment score, or -1.0 if the DB fragment is invalid.
*/
double ChainDB::score_distance( const ChainDB& frag, int offset ) const
{
  double score = 0.0;
  const int len = frag.dbdistvecs.size() - 1;
  for ( int i = 0; i < len; i++ ) {
    for ( int j = 0; j < len - i; j++ ) {
      const float& dbdist = dbdistvecs[offset+i].data[j];
      if ( dbdist <= 0.0 ) return -1.0;
      const float& frdist = frag.dbdistvecs[i].data[j];
      if ( frdist >  0.0 ) {
	const double d = double( dbdist - frdist );
	score += d * d;
      }
    }
  }
  return score;
}


/*!
  Return a fast distance score of a fragment against the nth fragment
  in the DB. If the nth fragment in the DB is broken across a chain
  boundary, then the return value is -1.0;

  The score is roughly comparable to a sum of squared differences
  between the search fragment and the DB fragment over all atom
  distances present in the fragment. The number of distances is equal
  to the number of pairs of non-null C-alphas in the searh fragment.

  This version of the scoring function contains an optimisation: it
  will terminate early and return -1 if the score exceeds the
  specified cutoff value. The cutoff value probably needs to be
  determined empirically: it probably varies with the square of the
  desired RMSD cutoff and the number of non-null C-alpha distances in
  the search fragment.

  \param frag The data for the search fragment
  \param offset The fragment in the DB against which to score
  \param scut Score cutoff value (optimisation)
  \return The fragment score, or -1.0 if the DB fragment is invalid.
*/
double ChainDB::score_distance( const ChainDB& frag, int offset, double scut ) const
{
  double score = 0.0;
  const int len = frag.dbdistvecs.size() - 1;
  for ( int i = 0; i < len; i++ ) {
    for ( int j = 0; j < len - i; j++ ) {
      const float& dbdist = dbdistvecs[offset+i].data[j];
      if ( dbdist <= 0.0 ) return -1.0;
      const float& frdist = frag.dbdistvecs[i].data[j];
      if ( frdist >  0.0 ) {
	const double d = double( dbdist - frdist );
	score += d * d;
	if ( score > scut ) return -1.0;
      }
    }
  }
  return score;
}


/*!
  Return a fast distance score of a fragment against the nth fragment
  in the DB. If the nth fragment in the DB is broken across a chain
  boundary, then the return value is -1.0;

  The score is roughly comparable to a sum of squared differences
  between the search fragment and the DB fragment over all atom
  distances present in the fragment. The number of distances is equal
  to the number of pairs of non-null C-alphas in the searh fragment.

  This version of the scoring function contains an optimisation: it
  will terminate early and return -1 if the score exceeds the
  specified cutoff value. The cutoff value probably needs to be
  determined empirically: it probably varies with the square of the
  desired RMSD cutoff and the number of non-null C-alpha distances in
  the search fragment.

  This version of the scoring function also takes a list of residue
  type masks, which act as a pre-filter on the DB fragment. If the DB
  fragment does not match the type mask at each position on the
  fragment, the return value is -1.0.

  \param frag The data for the search fragment
  \param types List of residue tpye masks.
  \param offset The fragment in the DB against which to score
  \param scut Score cutoff value (optimisation)
  \return The fragment score, or -1.0 if the DB fragment is invalid.
*/
double ChainDB::score_distance( const ChainDB& frag, const std::vector<Residue::TypeMask>& types, int offset, double scut ) const
{
  for ( int i = 0; i < types.size(); i++ )
    if ( !( Residue::TypeMask(dbresidues[offset+i].type()).mask() &
	    types[i].mask() ) ) return -1.0;
  return score_distance( frag, offset, scut );
}


/*!
  Return a list of fragment offsets from the DB which MAY match the
  given fragment.

  \param frag The data for the search fragment
  \param nhit The maximum number of intermediate matches to use.
*/
std::vector<int> ChainDB::match_fragment_preliminary( const ChainDB& fragdb, int nhit ) const
{
  const std::vector<Residue::TypeMask> types;
  return match_fragment_preliminary( fragdb, types, nhit );
}


/*!
  Return a list of fragment offsets from the DB which MAY match the
  given fragment.

  \param frag The data for the search fragment
  \param nhit The maximum number of intermediate matches to use.
  \param types Allowed types for each residue in the search fragment
*/
std::vector<int> ChainDB::match_fragment_preliminary( const ChainDB& fragdb, const std::vector<Residue::TypeMask>& types, int nhit ) const
{
  // find preliminary matching db fragments and sort
  std::vector<std::pair<double,int> > scores_dist;
  double scut = 1.0e20;
  for ( int i = 0; i < dbresidues.size(); i++ ) {
    double scr = score_distance( fragdb, types, i, scut );
    if ( scr >= 0.0 )
      scores_dist.push_back( std::pair<double,int>( scr, i ) );
    if ( scores_dist.size() >= 2*nhit ) {  // optimisation
      std::sort( scores_dist.begin(), scores_dist.end() );
      scores_dist.resize( nhit );
      scut = scores_dist.back().first;
    }
  }
  // sort the list of hits
  std::sort( scores_dist.begin(), scores_dist.end() );

  // construct result list
  std::vector<int> result;
  for ( int i = 0; i < nhit; i++ )
    if ( i < scores_dist.size() )
      result.push_back( scores_dist[i].second );
  return result;
}


/*!
  Return a list of fragments from the DB which match the given fragment.

  \param frag The data for the search fragment
  \param nlsq The maximum number of matches to return, sorted by RMSD
  \param nhit The maximum number of intermediate matches to use to get nlsq
  \return A list of up to nlsq matching fragments.

  If nhit is zero or omitted it is set to 10*nlsq. Occasionally a few
  of the best lsq fits may be missed, this can be improved by
  increasing nhit at a cost in performance.
*/
std::vector<Chain> ChainDB::match_fragment( const ChainDB& fragdb, int nlsq, int nhit ) const
{
  const std::vector<Residue::TypeMask> types;
  return match_fragment( fragdb, types, nlsq, nhit );
}


/*!
  Return a list of fragments from the DB which match the given fragment.

  \param frag The data for the search fragment
  \param types Allowed types for each residue in the search fragment
  \param nlsq The maximum number of matches to return, sorted by RMSD
  \param nhit The maximum number of intermediate matches to use to get nlsq
  \return A list of up to nlsq matching fragments.

  If nhit is zero or omitted it is set to 10*nlsq. Occasionally a few
  of the best lsq fits may be missed, this can be improved by
  increasing nhit at a cost in performance.
*/
std::vector<Chain> ChainDB::match_fragment( const ChainDB& fragdb, const std::vector<Residue::TypeMask>& types, int nlsq, int nhit ) const
{
  if ( nhit == 0 ) nhit = 10 * nlsq;

  // find preliminary matching db fragments and sort
  std::vector<int> scores_dist = match_fragment_preliminary( fragdb, types, nhit );

  // score fragments
  std::vector<std::pair<double,int> > scores_lsq;
  // for best matches, calc lsq supersposition
  for ( int i = 0; i < scores_dist.size(); i++ ) {
    const int offset = scores_dist[i];
    Chain fragnew = extract( offset, fragdb.size() );
    fragnew.lsq_superpose( fragdb );
    const double d = fragdb.rmsd( fragnew );
    scores_lsq.push_back( std::pair<double,int>( d, offset ) );
  }
  std::sort( scores_lsq.begin(), scores_lsq.end() );
  if ( scores_lsq.size() > nlsq ) scores_lsq.resize( nlsq );

  // build results
  std::vector<Chain> result;
  for ( int i = 0; i < scores_lsq.size(); i++ ) {
    const int offset = scores_lsq[i].second;
    Chain fragnew = extract( offset, fragdb.size() );
    fragnew.lsq_superpose( fragdb );
    result.push_back( fragnew );
  }

  return result;
}


} // namespace ProteinDB
