/* ligand/make-ligands-db.cc
 * 
 * Copyright 2016 by Medical Research Council
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
 * 02110-1301, USA
 */

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#ifdef USE_SQLITE3
#include <sqlite3.h>
#include "ligands-db.hh"
#endif
#include <iostream>


void print_usage(const std::string &argv_0) {

   std::cout << "Usage: " << argv_0 << "\n"
	     << "     --create  <create new database from ligand metrics file)\n"
	     << "     --ligand-metrics <ligand-metrics-table (default: ligand-metrics.tab)>\n"
	     << "     --database <db-file-name (default: ligands.db)>\n"
	     << "     --update-resolutions <resolutions-table-file-name>\n"
	     << "     --update-headers <headers-table-file-name>\n"
	     << "     --update-edstats <edstats-results-table-file-name>\n"
	     << "     --help   this output\n";

}

int main(int argc, char **argv) {

   int status = 0;
#ifdef USE_SQLITE3

   std::string resolutions_table;
   std::string headers_table;
   std::string database_name = "ligands.db";
   std::string ligand_metrics_tab_file_name = "ligand-metrics.tab";
   std::string edstats_wwpdb_results_tab_file_name;
   
   bool do_create = false;
   const char *optstr = "c:u:h";
   struct option long_options[] = {
      {"create",             0, 0, 0},
      {"ligand-metrics",     1, 0, 0},
      {"update-resolutions", 1, 0, 0},
      {"update-headers",     1, 0, 0},
      {"update-edstats",     1, 0, 0},  // from edstats table (file-name)
      {"database",           1, 0, 0},  // database file-name
      {"help",               0, 0, 0},
      {0, 0, 0, 0}
   };
      
   int ch;
   int option_index = 0;
   while ( -1 != 
	   (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

      switch(ch) {

      case 0:

	 if (! coot_optarg) { // no args for these options
	    std::string arg_str = long_options[option_index].name;
	    if (arg_str == "help") {
	       print_usage(argv[0]);
	    }
	    if (arg_str == "create") {
	       do_create = true;
	    }
	 } else {
	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "ligand-metrics") {
	       ligand_metrics_tab_file_name = coot_optarg; // Hmm.
	    }
	    if (arg_str == "database") {
	       database_name = coot_optarg;
	    }
	    if (arg_str == "update-resolutions") {
	       resolutions_table = coot_optarg;
	    }
	    if (arg_str == "update-headers") {
	       headers_table = coot_optarg;
	    }
	    if (arg_str == "update-edstats") {
	       edstats_wwpdb_results_tab_file_name = coot_optarg;
	    }
	 }
      }
   }

   if (argc == 1) {
      print_usage(argv[0]);
   } else { 
      if (do_create) {
	 coot::ligand_metrics lm;
	 std::cout << "calling parse_core_metrics(); " << ligand_metrics_tab_file_name
		   << " " << database_name << std::endl;
	 lm.parse_core_metrics(ligand_metrics_tab_file_name, database_name);
      }
      if (! resolutions_table.empty()) {
	 coot::ligand_metrics lm(database_name);
	 lm.update_resolutions(resolutions_table);
      }
      if (! headers_table.empty()) {
	 coot::ligand_metrics lm(database_name);
	 lm.update_headers(headers_table);
      }
      if (! edstats_wwpdb_results_tab_file_name.empty()) {
	 coot::ligand_metrics lm(database_name);
	 lm.update_edstats_results(edstats_wwpdb_results_tab_file_name);
      }
   }
   return status;
#else
   return 0;
#endif    

} 
