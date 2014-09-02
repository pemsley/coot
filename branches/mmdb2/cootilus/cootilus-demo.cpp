// Clipper cootilus-demo
/* Copyright 2012 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

#include "cootilus-build.h"
#include <stdlib.h>


int main( int argc, char** argv )
{
   if ( argc < 4 ) {
      std::cout << "cootaneer-demo <lib_file> <xmap_file> <mmdb_file> <radius>"
		<< std::endl;
      return 1;
   }

  std::string iplib( argv[1] );
  std::string ipmap( argv[2] );
  std::string ippdb( argv[3] );
  double radius = atof( argv[4] );

  // get map
  clipper::Xmap<float> xmap;
  clipper::CCP4MAPfile mapin;
  mapin.open_read( ipmap );
  mapin.import_xmap( xmap );
  mapin.close_read();

  // get initial coord
  clipper::MMDBfile pdbin;
  clipper::MiniMol mol_in;
  pdbin.read_file( ippdb );
  pdbin.import_minimol( mol_in );

  clipper::MMDBfile mmdb;
  
  Coot_nucleic_acid_build nabuild( iplib );
  nabuild.build( &mmdb, xmap, mol_in[0][0][0].coord_orth(), radius );

  mmdb.write_file( "cootilus.pdb" );
}
