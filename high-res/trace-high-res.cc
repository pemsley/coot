/* high-res/trace-high-res.cc
 * 
 * Copyright 2002, 2003, 2004 The University of York
 * Author Paul Emsley
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
 * 02110-1301, USA
 */

// Portability gubbins
#ifndef _MSC_VER
#include <unistd.h> // for getopt(3)
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <iostream>

#include "high-res.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/peak-search.hh"

int
main(int argc, char **argv) { 

   if (argc < 6) { 
      std::cout << "Usage: " << argv[0] 
		<< " --hklin mtz-filename"
		<< " --f f_col_label"
		<< " --phi phi_col_label"
		<< " --weight weight_col_label"
		<< " [--pdbout waters-filename]"
		<< " [ --sigma sigma-level]"
		<< "\n";

   } else { 

      short int use_weights = 0;
      std::string  mtz_filename;
      std::string         f_col;
      std::string       phi_col;
      std::string    weight_col;
      std::string    output_pdb;
      std::string     sigma_str;

      const char *optstr = "h:f:p:w:o:s";
      struct option long_options[] = {
	 {"hklin",  1, 0, 0},
	 {"f",      1, 0, 0},
	 {"phi",    1, 0, 0},
	 {"weight", 1, 0, 0},
	 {"pdbout", 1, 0, 0},
	 {"sigma",  1, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) { 

	 switch(ch) { 
	    
	 case 0:
	    if (optarg) { 
	       std::string arg_str = long_options[option_index].name;

	       if (arg_str == "pdbout") { 
		  output_pdb = optarg;
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
	       if (arg_str == "weight") {
		  weight_col = optarg;
	       } 
	       if (arg_str == "sigma") {
		  sigma_str = optarg;
	       } 
	       
	    } else { 
	       std::cout << "Malformed option: "
			 << long_options[option_index].name << std::endl;
	    }
	    break;

	 case 'o':
	    output_pdb = optarg;
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
	    
	 case 'w':
	    weight_col = optarg;
	    break;
	    
	 case 's':
	    sigma_str = optarg;
	    break;

	 default:
	    break;
	 }
      }

      short int do_it = 0;
      if (output_pdb.length() == 0) {
	 output_pdb = "a.pdb";
	 std::cout << "using default output filename: " << output_pdb << "\n";
      }

      if (weight_col.length() == 0) {
	 use_weights = 0;
      } else {
	 use_weights = 1;
      }

      if (mtz_filename.length() == 0) { 
	 std::cout << "Missing MTZ file\n";
	 exit(1);
      } else { 
	 if (f_col.length() == 0) { 
	    std::cout << "Missing F column name\n";
	    exit(1);
	 } else {
	    if (phi_col.length() == 0) { 
	       std::cout << "Missing PHI column name\n";
	       exit(1);
	    } else { 
	       if (sigma_str.length() == 0) {
		  sigma_str = "1.0";
	       }
	       do_it = 1;
	    }
	 }
      }

      if (do_it) { 
	 float sigma_level = atof(sigma_str.c_str());
	 short int is_diff_map = 0;
	 clipper::Xmap<float> xmap;

	 coot::util::map_fill_from_mtz(&xmap, mtz_filename, f_col, phi_col,
				       weight_col, use_weights, is_diff_map);
	 coot::peak_search ps(xmap);
	 std::vector<clipper::Coord_orth> ps_peaks = ps.get_peaks(xmap, 4.8);
	 ps.mask_map(&xmap, ps_peaks); // modify xmap
	 std::vector<clipper::Coord_orth> smaller_peaks = ps.get_peaks(xmap, sigma_level);

	 std::cout << "INFO:: Found " << ps_peaks.size() << " main peaks and "
		   << smaller_peaks.size() << " subsiduary peaks\n";

	 ps.add_peak_vectors(&ps_peaks, smaller_peaks);
	 coot::minimol::molecule ps_mol(ps_peaks, "HOH", " OW1", "W");
	 std::string spg(xmap.spacegroup().descr().symbol_hm());
	 ps_mol.set_spacegroup(spg);
	 ps_mol.set_cell(xmap.cell());
	 ps_mol.write_file("trace-peaksearch-peaks.pdb", 20.0);
	 coot::high_res hr(ps_mol);
	 hr.output_pdb("trace-coords.pdb");
      } else { 
	 std::cout << "not doing it\n";
      } 

      return 0;
   }
}
