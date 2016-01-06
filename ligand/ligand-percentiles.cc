
#ifdef USE_SQLITE3
#include <stdexcept>
#include <sqlite3.h>
#include "ligands-db.hh"
#include "utils/coot-utils.hh"
#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif
#endif

#include <iostream>

int main(int argc, char **argv) {

   int status = 0;

#ifdef USE_SQLITE3

   if (argc > 1) {
      std::string metric_str = "density_correlation";
      std::string value_str;
      std::string database_name = "new-ligands.db";
      bool do_help = false;

      const char *optstr = "m:d:v:h";
      struct option long_options[] = {
	 {"metric",        1, 0, 0},
	 {"database",      1, 0, 0},
	 {"value",         1, 0, 0},
	 {"update-resolutions", 1, 0, 0},
	 {"help",          0, 0, 0},
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

	       if (arg_str == "metric") {
		  metric_str = optarg;
	       } 
	       if (arg_str == "database") {
		  database_name = optarg;
	       } 
	       if (arg_str == "value") {
		  value_str = optarg;
	       }

	    } else {

	       // no options to this command line argument
	       std::string arg_str(long_options[option_index].name);
	       
	       if (arg_str == "help") { 
		  do_help = true;
	       } 
	    }
	    break;

	 case 'm':
	    metric_str = optarg;
	    break;
	    
	 case 'd':
	    database_name = optarg;
	    break;
	    
	 case 'v':
	    value_str = optarg;
	    break;
	    
	 case 'h':
	    do_help = true;
	    break;
	    
	 default:
	    std::cout << "default optarg: " << optarg << std::endl;
	    break;
	 }
      }

      if (do_help) {

	 std::cout << "Usage: " << argv[0] << "\n"
		   << "     --database <db-file-name>\n"
		   << "     --metric <metric-name>\n"
		   << "               e.g density_correlation, bumps_1, mogul_z_worst, diff_map_stat_1\n"
		   << "     --value <value>\n"
		   << "     --update-database <resolutions-table-file-name>\n"
		   << "     --help   this output\n";
	 
      } else {

	 try {
	    if (! value_str.empty()) { 
	       double dc = coot::util::string_to_float(value_str);
	       
	       coot::ligand_metrics lm(database_name);
	       std::pair<int, int> idx_pair = lm.get_index(dc, metric_str);
	       std::cout << "idx_pair: " << idx_pair.first << " " << idx_pair.second << std::endl;
	       if (idx_pair.second != 0) {
		  double ratio = double(idx_pair.first)/double(idx_pair.second);
		  std::cout << metric_str << " percentile for " << dc << " is " << 100 * ratio << "%"
			    << std::endl;
	       }
	    }
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << " error " << rte.what() << std::endl;
	 }
      }
   }
#endif // USE_SQLITE3 
   return status;
}
