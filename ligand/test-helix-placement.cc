// Portability (getopt) gubbins
#include <unistd.h> // for getopt(3)

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

      char *optstr = "i:h:f:p";
      struct option long_options[] = {
	 {"pdbin", 1, 0, 0},
	 {"hklin", 1, 0, 0},
	 {"f", 1, 0, 0},
	 {"phi", 1, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = getopt_long(argc, argv, optstr, long_options, &option_index))) { 

	 std::cout << "DEBUG:: " << option_index << " " << optarg << std::endl;

	 switch(ch) { 
	    
	 case 0:
	    if (optarg) { 
	       std::string arg_str = long_options[option_index].name;

	       if (arg_str == "pdbin") { 
		  pdb_file_name = optarg;
	       } 
	       if (arg_str == "hklin") { 
		  mtz_filename = optarg;
	       } 
	       if (arg_str == "f") { 
		  f_col = optarg;
	       } 
	       if (arg_str == "phi") {
		  phi_col = optarg;
	       } 
	       
	    } else { 
	       std::cout << "Malformed option: "
			 << long_options[option_index].name << std::endl;
	    }
	    break;

	 case 'i':
	    pdb_file_name = optarg;
	    break;
	    
	 case 'h':
	    mtz_filename = optarg;
	    break;
	    
	 case 'f':
	    f_col = optarg;
	    break;
	    
	 case 'p':
	    phi_col = optarg;
	    break;
	    
	 default:
	    std::cout << "default optarg: " << optarg << std::endl;
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
	       p.place_alpha_helix_near(pt, 20);
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

