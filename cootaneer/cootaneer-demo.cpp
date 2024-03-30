/*
 * cootaneer/cootaneer-demo.cpp
 *
 * Copyright 2003-2006 by Kevin Cowtan
 * Author: Kevin Cowtan
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */



#include <stdlib.h>  // needed for exit()

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "cootaneer-sequence.h"


int main( int argc, char** argv )
{
   if ( argc < 5 ) {
      std::cout << "cootaneer-demo <llk_data_file> <xmap_file> <mmdb_file> <chain_id> "
		<<  "<sequence>" << std::endl;
      exit(1);
   }

  std::string llkdfile( argv[1] );
  std::string xmapfile( argv[2] );
  std::string mmdbfile( argv[3] );
  std::string chain_id( argv[4] );
  std::string sequence( argv[5] );

  // read xmap
  clipper::Xmap<float> xmap;
  clipper::CCP4MAPfile mapfile;
  mapfile.open_read( xmapfile );
  mapfile.import_xmap( xmap );
  mapfile.close_read();

  // read mmdb
  mmdb::Manager mmdb;
  mmdb.SetFlag( mmdb::MMDBF_AutoSerials | mmdb::MMDBF_IgnoreDuplSeqNum );
  mmdb.ReadPDBASCII( mmdbfile.c_str() );

  // make sequence object
  std::vector<std::pair<std::string,std::string> > seq;
  seq.push_back( std::pair<std::string,std::string>( "A", sequence ) );

  // create sequencer
  Coot_sequence sequencer( llkdfile );

  // and apply
  sequencer.sequence_chain( xmap, seq, mmdb, chain_id );
  std::string bestseq = sequencer.best_sequence();
  std::string fullseq = sequencer.full_sequence();
  double conf = sequencer.confidence();
  int chnnum = sequencer.chain_number();
  int chnoff = sequencer.chain_offset();

  // write results
  std::cout << "\nSequence: " << bestseq << "\nCondidence: " << conf << "\n";
  if ( chnnum >= 0 ) std::cout << "\nFrom    : " << fullseq << "\nChain id: " << chnnum << "\tOffset: " << chnoff+1 << "\n";
}
