/* high-res/globularize.hh
 * 
 * Copyright 2004  The University of York
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

#ifdef __GNU_LIBRARY__
#include "coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <stdlib.h> // for atof
#include "mmdb-extras.h"
#include "mmdb.h"
#include "high-res.hh"

class glob_command_line_options {
public:
   std::string input_pdb_file;
   std::string output_pdb_file;
   short int have_centres; // 3 is the success number
   clipper::Coord_orth centre;
   glob_command_line_options() {
      have_centres = 0;
      input_pdb_file = "";
      output_pdb_file = "";
   }
   void parse_command_line_options(int argc, char **argv);
};


void 
glob_command_line_options::parse_command_line_options(int argc, char **argv) {

#ifdef _GETOPT_H // remove this?

   const char *optstr = "c:"; 
   static struct option long_options[] = {
      {"centre-x" , 1, 0, 0},
      {"centre-y" , 1, 0, 0},
      {"centre-z" , 1, 0, 0},
      {"pdb-in" , 1, 0, 0},
      {"pdb-out" , 1, 0, 0},
      {0, 0, 0, 0}	       // must have blanks at end
   };

   int option_index = 0; 
   int ch;
   short int have_x_flag = 0;
   short int have_y_flag = 0;
   short int have_z_flag = 0;
   double x = 0;
   double y = 0;
   double z = 0;

   while( -1 != 
	  (ch = getopt_long(argc, argv,optstr, long_options, &option_index) )) {

      switch(ch) {
	 
	 // a long option
	 
      case 0:
	 if (coot_optarg) {
	    
	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "centre-x") {
	       x = atof(coot_optarg);
	       have_x_flag = 1;
	    }
	    if (arg_str == "centre-y") {
	       y = atof(coot_optarg);
	       have_y_flag = 1;
	    }
	    if (arg_str == "centre-z") {
	       z = atof(coot_optarg);
	       have_z_flag = 1;
	    }
	    if (arg_str == "pdb-in") {
	       input_pdb_file = coot_optarg;
	    }
	    if (arg_str == "pdb-out") {
	       output_pdb_file = coot_optarg;
	    }
	 }
      }
   }
   if (have_x_flag)
      have_centres++;
   if (have_y_flag)
      have_centres++;
   if (have_z_flag)
      have_centres++;
   if (have_centres == 3) { 
      centre = clipper::Coord_orth(x, y, z);
   }

#endif // _GETOPT_H
   // pass back the number of command line options read:
}

int
main(int argc, char **argv) { 


   short int command_line_failure = 0;
   int inarg = 1;
   int outarg = 2;
      
   glob_command_line_options clo;
   if (argc > 2) {

      clo.parse_command_line_options(argc, argv);

      if (clo.input_pdb_file == "")
	 command_line_failure = 1;
      if (clo.output_pdb_file == "")
	 command_line_failure = 1;

      if (!((clo.have_centres == 3) || (clo.have_centres == 0)))
	 command_line_failure = 1;
      
   } else {
      command_line_failure = 1;
   }

   if (command_line_failure) { 
      std::cout << "Usage: " << argv[0]
		<< " [--centre-x x --centre-y y --centre-z z] "
		<< "--pdb-in pdb-in-filename "
		<< "--pdb-out pdb-out-filename\n"
		<< "   Only 0 or 3 --centre-? args is permissible"
		<< std::endl;
   } else {

      // Test here that SYMINFO is set and the file exists.

      
      std::string pdb_filename(clo.input_pdb_file);
      std::string out_pdb_filename(clo.output_pdb_file);
      atom_selection_container_t asc = get_atom_selection(pdb_filename, 1);
      if (asc.read_success) { 
	 coot::minimol::molecule mmol(asc.mol);
	 if (clo.have_centres == 3) {
	    std::cout << "Targetting protein at centre: "
		      << clo.centre.format() << std::endl;
	    coot::high_res hr(mmol, clo.centre);
	    hr.output_pdb(out_pdb_filename);
	 } else {
	    coot::high_res hr(mmol, 0);
	    hr.buccafilter();
	    hr.add_cbetas();
	    hr.add_os();
	    hr.output_pdb(out_pdb_filename);
	 }
      } else {
	 std::cout << " Failed to read input pdb file" << std::endl;
      } 
   }

   return 0;
}
