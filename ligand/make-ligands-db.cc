
#ifdef USE_SQLITE3
#include <sqlite3.h>
#include "ligands-db.hh"
#endif
#include <iostream>

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

void print_usage(const std::string &argv_0) {

   std::cout << "Usage: " << argv_0 << "\n"
	     << "     --create <ligand-metrics-table>\n"
	     << "     --database <db-file-name>\n"
	     << "     --update-resolutions <resolutions-table-file-name>\n"
	     << "     --update-headers <headers-table-file-name>\n"
	     << "     --help   this output\n";

}

int main(int argc, char **argv) {

   int status = 0;
#ifdef USE_SQLITE3

   std::string resolutions_table;
   std::string headers_table;
   std::string database_name = "ligands-2015.db";
   std::string ligand_metrics_tab_file_name = "ligand-metrics.tab";
   bool do_create = false;
   const char *optstr = "c:u:h";
   struct option long_options[] = {
      {"create",             1, 0, 0},
      {"update-resolutions", 1, 0, 0},
      {"update-headers", 1, 0, 0},
      {"database",           1, 0, 0},
      {"help",               0, 0, 0},
      {0, 0, 0, 0}
   };
      
   int ch;
   int option_index = 0;
   while ( -1 != 
	   (ch = getopt_long(argc, argv, optstr, long_options, &option_index))) {

      switch(ch) {

      case 0:

	 if (optarg) { 

	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "create") {
	       ligand_metrics_tab_file_name = optarg; // Hmm.
	       do_create = true;
	    }
	    if (arg_str == "database") {
	       database_name = optarg;
	    }
	    if (arg_str == "update-resolutions") {
	       resolutions_table = optarg;
	    }
	    if (arg_str == "update-headers") {
	       headers_table = optarg;
	    }
	 } else {
	    std::string arg_str = long_options[option_index].name;
	    if (arg_str == "help") {
	       print_usage(argv[0]);
	    }
	 }
      }
   }

   if (argc == 1) {
      print_usage(argv[0]);
   } else { 
      if (do_create) {
	 coot::ligand_metrics lm;
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
   }
   return status;
#else
   return 0;
#endif    

} 
