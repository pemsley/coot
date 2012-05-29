/*! NucleicAcidDB Top 500 main chain database */
/* (C) 2008-2009 Kevin Cowtan & University of York all rights reserved */


#include "nucleicacid_db.h"
#include <fstream>


namespace NucleicAcidDB {


NucleicAcid::NucleicAcid( const clipper::Coord_orth& cp, const clipper::Coord_orth& co5, const clipper::Coord_orth& cc5, const clipper::Coord_orth& cc4, const clipper::Coord_orth& co4, const clipper::Coord_orth& cc3, const clipper::Coord_orth& co3, const clipper::Coord_orth& cc2, const clipper::Coord_orth& cc1, const clipper::Coord_orth& cn, const clipper::String& type )
{
  clipper::String t = type + "?";
  typ = t.trim()[0];
  clipper::Util::set_null( p_x );
  clipper::Util::set_null( o5x );
  clipper::Util::set_null( c5x );
  clipper::Util::set_null( c4x );
  clipper::Util::set_null( o4x );
  clipper::Util::set_null( c3x );
  clipper::Util::set_null( o3x );
  clipper::Util::set_null( c2x );
  clipper::Util::set_null( c1x );
  clipper::Util::set_null( n_x );
  if ( !cp.is_null()  ) { p_x = cp.x(); p_y = cp.y(); p_z = cp.z(); }
  if ( !co5.is_null() ) { o5x = co5.x(); o5y = co5.y(); o5z = co5.z(); }
  if ( !cc5.is_null() ) { c5x = cc5.x(); c5y = cc5.y(); c5z = cc5.z(); }
  if ( !cc4.is_null() ) { c4x = cc4.x(); c4y = cc4.y(); c4z = cc4.z(); }
  if ( !co4.is_null() ) { o4x = co4.x(); o4y = co4.y(); o4z = co4.z(); }
  if ( !cc3.is_null() ) { c3x = cc3.x(); c3y = cc3.y(); c3z = cc3.z(); }
  if ( !co3.is_null() ) { o3x = co3.x(); o3y = co3.y(); o3z = co3.z(); }
  if ( !cc2.is_null() ) { c2x = cc2.x(); c2y = cc2.y(); c2z = cc2.z(); }
  if ( !cc1.is_null() ) { c1x = cc1.x(); c1y = cc1.y(); c1z = cc1.z(); }
  if ( !cn.is_null()  ) { n_x = cn.x(); n_y = cn.y(); n_z = cn.z(); }
  set_flag();
}

NucleicAcid::NucleicAcid( const clipper::MMonomer& mm )
{
  clipper::String t = mm.type() + "?";
  typ = t.trim()[0];
  clipper::Util::set_null( p_x );
  clipper::Util::set_null( o5x );
  clipper::Util::set_null( c5x );
  clipper::Util::set_null( c4x );
  clipper::Util::set_null( o4x );
  clipper::Util::set_null( c3x );
  clipper::Util::set_null( o3x );
  clipper::Util::set_null( c2x );
  clipper::Util::set_null( c1x );
  clipper::Util::set_null( n_x );
  int ip  = mm.lookup( " P  ", clipper::MM::ANY );
  int io5 = mm.lookup( " O5'", clipper::MM::ANY );
  int ic5 = mm.lookup( " C5'", clipper::MM::ANY );
  int ic4 = mm.lookup( " C4'", clipper::MM::ANY );
  int io4 = mm.lookup( " O4'", clipper::MM::ANY );
  int ic3 = mm.lookup( " C3'", clipper::MM::ANY );
  int io3 = mm.lookup( " O3'", clipper::MM::ANY );
  int ic2 = mm.lookup( " C2'", clipper::MM::ANY );
  int ic1 = mm.lookup( " C1'", clipper::MM::ANY );
  int in  = mm.lookup( " N9 ", clipper::MM::ANY );
  if ( in < 0 ) in = mm.lookup( " N1 ", clipper::MM::ANY );
  if ( ip  >= 0 ) {
    p_x = mm[ip].coord_orth().x();
    p_y = mm[ip].coord_orth().y();
    p_z = mm[ip].coord_orth().z();
  }
  if ( io5 >= 0 ) {
    o5x = mm[io5].coord_orth().x();
    o5y = mm[io5].coord_orth().y();
    o5z = mm[io5].coord_orth().z();
  }
  if ( ic5 >= 0 ) {
    c5x = mm[ic5].coord_orth().x();
    c5y = mm[ic5].coord_orth().y();
    c5z = mm[ic5].coord_orth().z();
  }
  if ( ic4 >= 0 ) {
    c4x = mm[ic4].coord_orth().x();
    c4y = mm[ic4].coord_orth().y();
    c4z = mm[ic4].coord_orth().z();
  }
  if ( io4 >= 0 ) {
    o4x = mm[io4].coord_orth().x();
    o4y = mm[io4].coord_orth().y();
    o4z = mm[io4].coord_orth().z();
  }
  if ( ic3 >= 0 ) {
    c3x = mm[ic3].coord_orth().x();
    c3y = mm[ic3].coord_orth().y();
    c3z = mm[ic3].coord_orth().z();
  }
  if ( io3 >= 0 ) {
    o3x = mm[io3].coord_orth().x();
    o3y = mm[io3].coord_orth().y();
    o3z = mm[io3].coord_orth().z();
  }
  if ( ic2 >= 0 ) {
    c2x = mm[ic2].coord_orth().x();
    c2y = mm[ic2].coord_orth().y();
    c2z = mm[ic2].coord_orth().z();
  }
  if ( ic1 >= 0 ) {
    c1x = mm[ic1].coord_orth().x();
    c1y = mm[ic1].coord_orth().y();
    c1z = mm[ic1].coord_orth().z();
  }
  if ( in  >= 0 ) {
    n_x = mm[in].coord_orth().x();
    n_y = mm[in].coord_orth().y();
    n_z = mm[in].coord_orth().z();
  }
  set_flag();
}

clipper::Coord_orth NucleicAcid::coord_p () const
{
  if ( clipper::Util::is_null( p_x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( p_x, p_y, p_z );
}

clipper::Coord_orth NucleicAcid::coord_o5() const
{
  if ( clipper::Util::is_null( o5x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( o5x, o5y, o5z );
}

clipper::Coord_orth NucleicAcid::coord_c5() const
{
  if ( clipper::Util::is_null( c5x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c5x, c5y, c5z );
}

clipper::Coord_orth NucleicAcid::coord_c4() const
{
  if ( clipper::Util::is_null( c4x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c4x, c4y, c4z );
}

clipper::Coord_orth NucleicAcid::coord_c3() const
{
  if ( clipper::Util::is_null( c3x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c3x, c3y, c3z );
}

clipper::Coord_orth NucleicAcid::coord_o3() const
{
  if ( clipper::Util::is_null( o3x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( o3x, o3y, o3z );
}

clipper::Coord_orth NucleicAcid::coord_c2() const
{
  if ( clipper::Util::is_null( c2x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c2x, c2y, c2z );
}

clipper::Coord_orth NucleicAcid::coord_c1() const
{
  if ( clipper::Util::is_null( c1x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( c1x, c1y, c1z );
}

clipper::Coord_orth NucleicAcid::coord_o4() const
{
  if ( clipper::Util::is_null( o4x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( o4x, o4y, o4z );
}

clipper::Coord_orth NucleicAcid::coord_n () const
{
  if ( clipper::Util::is_null( n_x ) )
    return clipper::Coord_orth( clipper::Coord_orth::null() );
  else
    return clipper::Coord_orth( n_x, n_y, n_z );
}

clipper::MMonomer NucleicAcid::mmonomer() const
{
  clipper::MMonomer mm;
  clipper::MAtom ma_p, ma_n, ma_c, ma_o;
  clipper::MAtom mao5, mac5, mac4, mac3, mao3, mac2, mac1, mao4;
  ma_p = ma_n = ma_c = ma_o = clipper::MAtom::null();
  ma_p.set_u_iso ( 0.25 ); ma_p.set_occupancy( 1.0 );
  ma_p.set_id( " P  " ); ma_p.set_element( "P" );
  ma_n.set_u_iso ( 0.25 ); ma_n.set_occupancy( 1.0 );
  ma_n.set_id( " N1 " ); ma_n.set_element( "N" );
  ma_c.set_u_iso ( 0.25 ); ma_c.set_occupancy( 1.0 );
  ma_c.set_id( " C  " ); ma_c.set_element( "C" );
  ma_o.set_u_iso ( 0.25 ); ma_o.set_occupancy( 1.0 );
  ma_o.set_id( " O  " ); ma_o.set_element( "O" );
  mao5 = mao4 = mao3 = ma_o;
  mac5 = mac4 = mac3 = mac2 = mac1 = ma_c;
  mao5.set_id( " O5'" ); mao4.set_id( " O4'" ); mao3.set_id( " O3'" ); 
  mac5.set_id( " C5'" ); mac4.set_id( " C4'" ); mac3.set_id( " C3'" ); 
  mac2.set_id( " C2'" ); mac1.set_id( " C1'" );
  ma_p.set_coord_orth( coord_p() );
  ma_n.set_coord_orth( coord_n() );
  mao5.set_coord_orth( coord_o5() );
  mao4.set_coord_orth( coord_o4() );
  mao3.set_coord_orth( coord_o3() );
  mac5.set_coord_orth( coord_c5() );
  mac4.set_coord_orth( coord_c4() );
  mac3.set_coord_orth( coord_c3() );
  mac2.set_coord_orth( coord_c2() );
  mac1.set_coord_orth( coord_c1() );
  if ( !ma_p.coord_orth().is_null() ) mm.insert( ma_p );
  if ( !mao5.coord_orth().is_null() ) mm.insert( mao5 );
  if ( !mac5.coord_orth().is_null() ) mm.insert( mac5 );
  if ( !mac4.coord_orth().is_null() ) mm.insert( mac4 );
  if ( !mao4.coord_orth().is_null() ) mm.insert( mao4 );
  if ( !mac3.coord_orth().is_null() ) mm.insert( mac3 );
  if ( !mao3.coord_orth().is_null() ) mm.insert( mao3 );
  if ( !mac2.coord_orth().is_null() ) mm.insert( mac2 );
  if ( !mac1.coord_orth().is_null() ) mm.insert( mac1 );
  if ( !ma_n.coord_orth().is_null() ) mm.insert( ma_n );
  mm.set_type( std::string( 1, typ ) );
  return mm;
}


void NucleicAcid::transform( const clipper::RTop_orth& rtop )
{
  if ( !clipper::Util::is_null( p_x ) ) {
    clipper::Coord_orth c = rtop * coord_p();
    p_x = float( c.x() );
    p_y = float( c.y() );
    p_z = float( c.z() );
  }
  if ( !clipper::Util::is_null( o5x ) ) {
    clipper::Coord_orth c = rtop * coord_o5();
    o5x = float( c.x() );
    o5y = float( c.y() );
    o5z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c5x ) ) {
    clipper::Coord_orth c = rtop * coord_c5();
    c5x = float( c.x() );
    c5y = float( c.y() );
    c5z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c4x ) ) {
    clipper::Coord_orth c = rtop * coord_c4();
    c4x = float( c.x() );
    c4y = float( c.y() );
    c4z = float( c.z() );
  }
  if ( !clipper::Util::is_null( o4x ) ) {
    clipper::Coord_orth c = rtop * coord_o4();
    o4x = float( c.x() );
    o4y = float( c.y() );
    o4z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c3x ) ) {
    clipper::Coord_orth c = rtop * coord_c3();
    c3x = float( c.x() );
    c3y = float( c.y() );
    c3z = float( c.z() );
  }
  if ( !clipper::Util::is_null( o3x ) ) {
    clipper::Coord_orth c = rtop * coord_o3();
    o3x = float( c.x() );
    o3y = float( c.y() );
    o3z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c2x ) ) {
    clipper::Coord_orth c = rtop * coord_c2();
    c2x = float( c.x() );
    c2y = float( c.y() );
    c2z = float( c.z() );
  }
  if ( !clipper::Util::is_null( c1x ) ) {
    clipper::Coord_orth c = rtop * coord_c1();
    c1x = float( c.x() );
    c1y = float( c.y() );
    c1z = float( c.z() );
  }
  if ( !clipper::Util::is_null( n_x ) ) {
    clipper::Coord_orth c = rtop * coord_n();
    n_x = float( c.x() );
    n_y = float( c.y() );
    n_z = float( c.z() );
  }
}


void NucleicAcid::set_flag()
{
  if ( !clipper::Util::is_null( c1x ) &&
       !clipper::Util::is_null( c3x ) &&
       !clipper::Util::is_null( c4x ) ) {
    if ( !clipper::Util::is_null( n_x ) &&
	 !clipper::Util::is_null( p_x ) &&
	 !clipper::Util::is_null( c2x ) &&
	 !clipper::Util::is_null( c5x ) &&
	 !clipper::Util::is_null( o3x ) &&
	 !clipper::Util::is_null( o4x ) &&
	 !clipper::Util::is_null( o5x ) ) {
      flg = COMPLETE;
    } else {
      flg = INCOMPLETE;
    }
  }
  else {
    flg = NONE;
  }
}


// Chain classes

bool Chain::add_pdb( const clipper::String file )
{
  const int mmdbflags = MMDBF_IgnoreBlankLines | MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreNonCoorPDBErrors | MMDBF_IgnoreRemarks;
  clipper::MMDBfile mfile;
  clipper::MiniMol mol;
  mfile.SetFlag( mmdbflags );
  mfile.read_file( file );
  mfile.import_minimol( mol );
  if ( mol.size() == 0 ) return false;
  for ( int c = 0; c < mol.size(); c++ ) {
    clipper::MPolymer mp;
    // select monomers by occupancy
    for ( int r = 0; r < mol[c].size(); r++ ) {
      if ( mol[c][r].lookup( " C1'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " C2'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " C3'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " C4'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " C5'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " O3'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " O4'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " O5'", clipper::MM::ANY ) >= 0 &&
	   mol[c][r].lookup( " P  ", clipper::MM::ANY ) >= 0 ) {
	int a = mol[c][r].lookup( " C4'", clipper::MM::ANY );
	if ( mol[c][r][a].occupancy() > 0.01 &&
	     mol[c][r][a].u_iso() < clipper::Util::b2u(100.0) )
	  mp.insert( mol[c][r] );
      }
    }
    // shift centre-of-mass of chain to the origin
    clipper::Coord_orth cm( 0.0, 0.0, 0.0 );
    double              sm = 0.0;
    for ( int r = 0; r < mp.size(); r++ ) {
      int a = mp[r].lookup( " C4'", clipper::MM::ANY );
      cm += mp[r][a].coord_orth();
      sm += 1.0;
    }
    cm = (-1.0/sm) * cm;
    clipper::RTop_orth rt( clipper::Mat33<>::identity(), cm );
    mp.transform( rt );
    // now add the chain to the db
    for ( int r = 0; r < mp.size(); r++ ) {
      NucleicAcid rp( mp[r] );
      if ( rp.flag() == NucleicAcid::COMPLETE && rp.type() != ' ' )
	add_monomer( rp );
    }
  }
  return true;
}


bool Chain::save_db( const clipper::String file ) const
{
  /*
  std::ofstream fs( file.c_str(), std::ios::out | std::ios::binary );
  char d[20];
  for ( int r = 0; r < dbmonomers.size(); r++ ) {
    dbmonomers[r].data_export( d );
    fs.write( d, 20 );
  }
  fs.close();
  */
  return true;
}


bool Chain::load_db( const clipper::String file )
{
  /*
  dbmonomers.clear();
  // read whole file (for speed)
  std::ifstream fs( file.c_str(), std::ios::in | std::ios::binary );
  if ( !fs ) return false;
  fs.seekg( 0, std::ios::end );
  int i2 = fs.tellg();
  fs.seekg( 0, std::ios::beg );
  int i1 = fs.tellg();
  int l = i2 - i1;
  char d[l];
  fs.read( d, l );
  fs.close();
  if ( l%20 != 0 ) return false;
  // import file data
  dbmonomers.resize( l/20 );
  for ( int r = 0; r < dbmonomers.size(); r++ ) {
    dbmonomers[r].data_import( d + 20*r );
  }
  */
  return true;
}


bool Chain::merge( const Chain& other, const std::vector<double>& wgt )
{
  /*
  if ( other.size() != size() ) return false;
  if ( wgt.size() != 3*size() ) return false;
  for ( int r = 0; r < dbmonomers.size(); r++ )
    dbmonomers[r].merge( other.dbmonomers[r],
			 wgt[3*r], wgt[3*r+1], wgt[3*r+2] );
  */
  return true;
}


Chain Chain::extract( int offset, int len ) const
{
  Chain dbc;
  for ( int i = 0; i < len; i++ ) dbc.add_monomer( dbmonomers[offset+i] );
  return dbc;
}


bool Chain::is_continuous() const
{
  // go through and find elements where there is a chain break
  const double dmin = 2.0;
  std::vector<bool> cterm( dbmonomers.size(), false );
  for ( int i = 0; i < dbmonomers.size()-1; i++ ) {
    int j = i + 1;
    const clipper::Coord_orth co1 = dbmonomers[i].coord_o3();
    const clipper::Coord_orth co2 = dbmonomers[j].coord_p();
    if ( co1.is_null() || co2.is_null() ) return false;
    const double d2 = ( co1 - co2 ).lengthsq();
    if ( d2 > dmin*dmin ) return false;
  }
  return true;
}


void Chain::transform( const clipper::RTop_orth& rtop )
{
  for ( int r = 0; r < dbmonomers.size(); r++ )
    dbmonomers[r].transform( rtop );
}


}
