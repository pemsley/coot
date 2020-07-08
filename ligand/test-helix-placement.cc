/* ligand/test-helix-placement.cc
 * 
 * Copyright 2005 by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 * 02110-1301, USA.
 */

// Portability (getopt) gubbins
#include <unistd.h> // for getopt(3)
#if !defined WINDOWS_MINGW
// conflicts in Windows and doesnt seem to be needed anyway
// actually not sure if needed at all...
#include <stdlib.h> // for atof
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <iostream>
#include <string>

#include "coot-map-utils.hh"
#include "helix-placement.hh"

int main(int argc, char **argv) {

   clipper::Xmap<float> xmap;
   int istatus = 0;

   if (argc < 7) {
      istatus = 1;
      std::cout << "Usage: " << argv[0]
	 // 		<< " --pdbin pdb-in-filename "
		<< " --hklin mtz-filename "
		<< " --f f_col_label"
		<< " --phi phi_col_label"
		<< std::endl;
   } else {
      std::string pdb_file_name("");
      std::string  mtz_filename("");
      std::string         f_col("");
      std::string       phi_col("");
      short int is_diff_map = 0;
      bool have_x_flag = 0; 
      bool have_y_flag = 0; 
      bool have_z_flag = 0;
      float x = 0.0; 
      float y = 0.0; 
      float z = 0.0; 


      // There is something broken with passing -x 4 -y 5 -z 6.  Understand why and fix later.
      // 
      const char *optstr = "i:h:f:p:s:x:y:z";
      struct option long_options[] = {
	 {"pdbin",  1, 0, 0},
	 {"hklin",  1, 0, 0},
	 {"f",      1, 0, 0},
	 {"phi",    1, 0, 0},
	 {"x",    1, 0, 0},
	 {"y",    1, 0, 0},
	 {"z",    1, 0, 0},
	 {"strand", 0, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      bool do_strand = 0;
      while ( -1 != 
	      (ch = getopt_long(argc, argv, optstr, long_options, &option_index))) { 

// 	 if (coot_optarg)
// 	    std::cout << "DEBUG:: index " << option_index << " " << coot_optarg << std::endl;
// 	 else 
// 	    std::cout << "DEBUG:: index " << option_index << " null coot_optarg" << std::endl;

	 switch(ch) { 

	    // long arguments, no ch
	 case 0:
	    // std::cout << "processing... " << long_options[option_index].name << std::endl;
	    if (coot_optarg) {
	       
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "pdbin") { 
		  pdb_file_name = coot_optarg;
	       } 
	       if (arg_str == "hklin") { 
		  mtz_filename = coot_optarg;
	       } 
	       if (arg_str == "f") { 
		  f_col = coot_optarg;
	       } 
	       if (arg_str == "phi") {
		  phi_col = coot_optarg;
	       }
	       if (arg_str == "x") {
		  have_x_flag = 1;
		  x = atof(coot_optarg);
	       }
	       if (arg_str == "y") {
		  have_y_flag = 1;
		  y = atof(coot_optarg);
	       }
	       if (arg_str == "z") {
		  have_z_flag = 1;
		  z = atof(coot_optarg);
	       }
	       
	    } else { 
	       // long arg, no coot_optarg
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "strand") {
		  do_strand = 1;
	       }
	    }
	    break;

	    // 1-char args:
	 case 'i':
	    pdb_file_name = coot_optarg;
	    break;
	    
	 case 'h':
	    mtz_filename = coot_optarg;
	    break;
	    
	 case 'f':
	    f_col = coot_optarg;
	    break;
	    
	 case 'p':
	    phi_col = coot_optarg;
	    break;

	 case 'x':
	    if (coot_optarg) { 
	       //std::cout << "processing... x " << std::endl;
	       have_x_flag = 1;
	       x = atof(coot_optarg);
	       //std::cout << "set value of x " <<  x << std::endl;
	    } else {
	       //std::cout << "no coot_optarg for x!" << std::endl;
	    } 
	    
	 case 'y':
	    if (coot_optarg) { 
	       //std::cout << "processing... y " << std::endl;
	       have_y_flag = 1;
	       y = atof(coot_optarg);
	       // std::cout << "set value of y " <<  y << std::endl;
	    } else { 
	       //std::cout << "no coot_optarg for y!" << std::endl;
	    } 
	    
	 case 'z':
	    if (coot_optarg) { 
	       //std::cout << "processing... z " << std::endl;
	       have_z_flag = 1;
	       //std::cout << "seting value of z from " <<  coot_optarg << std::endl;
	       z = atof(coot_optarg);
	       //std::cout << "set value of z " <<  z << std::endl;
	    } else {
	       //std::cout << "no coot_optarg for z!" << std::endl;
	    } 
	    
	 default:
	    std::cout << "default coot_optarg: " << coot_optarg << std::endl;
	    break;
	 }
      }

      if (mtz_filename.length() > 0) {

	 if (f_col.length() > 0) { 
	    if (phi_col.length() > 0) { 
   
	       short int pdb_file_name_was_set = 0;
	       if (pdb_file_name.length() > 0)
		  pdb_file_name_was_set = 1;

	       std::string weight_col("");
	       short int use_weights = 0;
	       coot::util::map_fill_from_mtz(&xmap, mtz_filename,
					     f_col, phi_col,
					     weight_col, use_weights,
					     is_diff_map);

	       coot::helix_placement p(xmap);
	       // clipper::Coord_orth pt(23, 34, 45);
	       // clipper::Coord_orth pt(34.4170570, 54.75373, 52.24235);

	       // for 1cos-A:
	       // clipper::Coord_orth pt(8, -11, 14); tricky
	       clipper::Coord_orth pt(10, -7, 9);

	       // did the user specify a position?
	       if (have_x_flag && have_y_flag && have_z_flag) {
		  pt = clipper::Coord_orth(x, y, z);
	       } 
	       if (do_strand == 0) { 
		  p.place_alpha_helix_near(pt, 20, 0.2, 30);
	       } else { 
		  coot::helix_placement_info_t si = p.place_strand(pt, 6, 5, 0.2);
		  if (si.success) {
		     si.mol[0].write_file("strand.pdb", 30.0);
		  } else {
		     std::cout << "Strand placement failed." << std::endl;
		  } 
	       } 

	       // p.discrimination_map();
	       
	    } else {
	       std::cout << "must specify --phi column label" << std::endl;
	    } 
	 } else {
	    std::cout << "must specify --f column label" << std::endl;
	 }
      } else {
	 std::cout << "must specify --hklin file-name" << std::endl;
      }
   }
   return istatus; 
}

