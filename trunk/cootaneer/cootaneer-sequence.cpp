/*! \file cootaneer-sequence.cpp buccaneer library */
/* (C) 2006 Kevin Cowtan & University of York all rights reserved */

#include "cootaneer-sequence.h"

#include <stdlib.h>  // for labs
#include <iostream>
#include <fstream>


Coot_sequence::Coot_sequence( std::string filename )
{
  llksmp = read_targets( filename );
}


std::pair<std::string,double> Coot_sequence::sequence_chain( const clipper::Xmap<float>& xmap, const std::vector<std::pair<std::string,std::string> >& sequence, CMMDBManager& mmdb, std::string chain_id )
{
  // convert the sequence to minimol form
  clipper::MMoleculeSequence seq;
  for ( int c = 0; c < sequence.size(); c++ ) {
    clipper::MPolymerSequence pol;
    pol.set_id( sequence[c].first );
    pol.set_sequence( sequence[c].second );
    seq.insert( pol );
  }
  // convert the molecule to minimol form
  clipper::MiniMol mmol;
  static_cast<clipper::MMDBfile&>(mmdb).import_minimol( mmol );
  // select the chain
  clipper::MPolymer mchn = mmol.find( chain_id );
  // sequence it
  Score_list<clipper::String> scores = 
    Ca_sequence::sequence_chain( mchn, xmap, llksmp, seq );

  // store resutls
  bestseq = scores[0];
  prob = Ca_sequence::phi_approx(scores.score(1) - scores.score(0));

  // get renumbering and chain info
  fullseq = "";
  std::pair<int,int> info = ProteinTools::chain_sequence_match( bestseq, seq );
  bestchn = info.first;
  bestoff = info.second;
  if ( bestchn >= 0 ) {
    std::string nullseq( mchn.size(), '?' );
    fullseq = nullseq + seq[bestchn].sequence() + nullseq;
    fullseq = fullseq.substr( mchn.size() + bestoff, mchn.size() );
  }

  // assemble result
  std::pair<std::string,double> result;
  result.first = bestseq;
  result.second = prob;
  return result;
}


std::pair<signed char, signed char> Coot_sequence::pack( double d )
{
  std::pair<signed char, signed char> pc;
  int m, e;
  m = e = 0;
  while ( fabs(d) <= 0.5 ) {
    d *= 2.0;
    e--;
    if ( labs(e) >= 127 ) break;
  }
  while ( fabs(d) >= 1.0 ) {
    d *= 0.5;
    e++;
    if ( labs(e) >= 127 ) break;
  }
  if ( labs(e) >= 127 ) m = e = 0;
  else                 m = int( 128.0 * d );
  if ( m >=  128 ) m =  127;
  if ( m <= -128 ) m = -127;
  pc.first  = m;
  pc.second = e;
  return pc;
}

double Coot_sequence::unpack( std::pair<signed char, signed char> pc )
{
  int m = int( pc.first );
  int e = int( pc.second );
  double d = ( m / 128.0 );
  while ( e < 0 ) {
    d *= 0.5;
    e++;
  }
  while ( e > 0 ) {
    d *= 2.0;
    e--;
  }
  return d;
}


void Coot_sequence::write_targets( std::string name , const std::vector<LLK_map_target::Sampled> llksmp )
{
  int nt = 20;
  std::ofstream file ( name.c_str(), std::ios::binary|std::ios::out|std::ios::trunc );
  if ( file.is_open() ) {
    // extract
    typedef std::pair<signed char, signed char> PC;
    int nt = 20;
    int ls = 2*(3+2*nt);
    int ns = llksmp[0].size();
    int size = ns*ls;
    char *memblock = new char[size];
    for ( int s = 0; s < ns; s++ ) {
      PC x, y, z, tgt, wgt;
      x = pack( llksmp[0].coord_orth(s).x() );
      y = pack( llksmp[0].coord_orth(s).y() );
      z = pack( llksmp[0].coord_orth(s).z() );
      memblock[ls*s+0] = x.first;
      memblock[ls*s+1] = x.second;
      memblock[ls*s+2] = y.first;
      memblock[ls*s+3] = y.second;
      memblock[ls*s+4] = z.first;
      memblock[ls*s+5] = z.second;
      for ( int t = 0; t < 20; t++ ) {
	tgt = pack( llksmp[t].target(s) );
	wgt = pack( llksmp[t].weight(s) );
	memblock[ls*s+6+4*t+0] = tgt.first;
	memblock[ls*s+6+4*t+1] = tgt.second;
	memblock[ls*s+6+4*t+2] = wgt.first;
	memblock[ls*s+6+4*t+3] = wgt.second;
      }
    }
    // write data
    file.write( memblock, size );
    file.close();
    // clean up
    delete[] memblock;    
  }
}

std::vector<LLK_map_target::Sampled> Coot_sequence::read_targets( std::string name )
{
  std::vector<LLK_map_target::Sampled> result;
  std::ifstream file ( name.c_str(), std::ios::binary|std::ios::in );
  std::ifstream::pos_type size;
  if ( file.is_open() ) {
    // read data
    file.seekg( 0, std::ios::end );
    size = file.tellg();
    char *memblock = new char[size];
    file.seekg( 0, std::ios::beg );
    file.read( memblock, size );
    file.close();
    // extract
    typedef std::pair<signed char, signed char> PC;
    int nt = 20;
    int ls = 2*(3+2*nt);
    int ns = int(size)/ls;
    result.resize( nt );
    for ( int t = 0; t < nt; t++ )
      result[t].set_type( LLK_map_target::CORREL );
    for ( int s = 0; s < ns; s++ ) {
      double x, y, z, tgt, wgt;
      x = unpack( PC( memblock[ls*s+0], memblock[ls*s+1] ) );
      y = unpack( PC( memblock[ls*s+2], memblock[ls*s+3] ) );
      z = unpack( PC( memblock[ls*s+4], memblock[ls*s+5] ) );
      clipper::Coord_orth co( x, y, z );
      for ( int t = 0; t < nt; t++ ) {
	tgt = unpack( PC( memblock[ls*s+6+4*t+0], memblock[ls*s+6+4*t+1] ) );
	wgt = unpack( PC( memblock[ls*s+6+4*t+2], memblock[ls*s+6+4*t+3] ) );
	result[t].insert( co, tgt, wgt );
      }
    }
    // clean up
    delete[] memblock;
  }
  return result;
}
