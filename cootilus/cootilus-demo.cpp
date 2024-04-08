/*
 * cootilus/cootilus-demo.cpp
 *
 * Copyright 2012 by Kevin Cowtan
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
