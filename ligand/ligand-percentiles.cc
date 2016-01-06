
#ifdef USE_SQLITE3
#include <stdexcept>
#include <map>
#include <algorithm>
#include <math.h>
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


#ifdef USE_SQLITE3

std::pair<int, int> get_index(double val, const std::vector<double> &v, bool low_is_good) {

   int len = 0;
   int idx = 0;

   if (low_is_good) { 
      for (unsigned int i=0; i<v.size(); i++) { 
	 if (val <= v[i]) {
	    len = v.size();
	    idx = len -i;
	    break;
	 }
      }
   } else {
      for (unsigned int i=0; i<v.size(); i++) { 
	 if (val < v[i]) {
	    idx = i;
	    len = v.size();
	    break;
	 }
      }
   }
   return std::pair<int, int>(idx, len);
}

bool string_double_sorter(const std::pair<std::string, double> &p1,
			  const std::pair<std::string, double> &p2) {

   bool r = false;
   return (p1.second<p2.second);
} 

void ranks(const std::string &database_name) { 

   coot::ligand_metrics lm(database_name);
   std::vector<std::string> metric_names(4);
   metric_names[0] = "density_correlation";
   metric_names[1] = "coot_diff_map_correlation";
   metric_names[2] = "bumps_1";
   metric_names[3] = "mogul_z_worst";
   // std::string metric_5 = "bfactor_stat_1"; // not in database yet!
   // needs ligand-metrics import update

   std::vector<std::string> primary_keys = lm.get_primary_keys();

   std::cout << "got " << primary_keys.size() << " primary keys" << std::endl;

   std::vector<std::pair<std::string, double> > tri_ranks;
   std::map<std::string, std::vector<double> > values_store;

   for (unsigned int i=0; i<metric_names.size(); i++) {
      const std::string &metric_name = metric_names[i];
      std::vector<double> v = lm.get_values(metric_name);
      if (metric_name == "coot_diff_map_correlation")
	 for (unsigned int i=0; i<v.size(); i++)
	    v[i] = fabs(v[i]);
      std::sort(v.begin(), v.end());
      values_store[metric_name] = v;
   }

   for (unsigned int i=0; i<primary_keys.size(); i++) {
      std::vector<int> rank_indices;
      const std::string &code = primary_keys[i]; // accession_code
      for (unsigned int j=0; j<metric_names.size(); j++) {
	 const std::string &metric_name = metric_names[j];

	 bool from_int = false;
	 if (metric_name == "bumps_1") from_int = true;
	 double value = lm.get_value(code, metric_name, from_int);
	 
	 bool reverse_order = false;
	 if (metric_name == "coot_diff_map_correlation") reverse_order = true;
	 if (metric_name == "bumps_1") reverse_order = true;
	 if (metric_name == "mogul_z_worst") reverse_order = true;
	 std::pair<int, int> idx_pair = get_index(value, values_store[metric_name], reverse_order);

	 std::cout << code << " " << metric_name << " " << value << "   "
		   << idx_pair.first << " / " << idx_pair.second << std::endl;

	 if (idx_pair.second != 0) {
	    rank_indices.push_back(idx_pair.first);
	 } else {
	    std::cout << "bad idx " << code << " " << metric_name << std::endl;
	 } 
      }

      if (rank_indices.size() == metric_names.size()) {
 	 double vi = 0;
 	 for (unsigned int i=0; i<metric_names.size(); i++) {
 	    vi += rank_indices[i] * rank_indices[i];
 	 }
 	 tri_ranks.push_back(std::pair<std::string, double> (code, vi));
      }
   }

   std::sort(tri_ranks.begin(), tri_ranks.end(), coot::util::sd_compare);
      
   std::vector<std::pair<std::string, double> >::const_iterator it;
   for (it = tri_ranks.begin(); it != tri_ranks.end(); it++) {
      std::cout << it->first << " " << it->second << std::endl;
   }
   
}

#endif // USE_SQLITE3


int main(int argc, char **argv) {

   int status = 0;

#ifdef USE_SQLITE3

   if (argc > 1) {
      std::string metric_str = "density_correlation";
      std::string value_str;
      std::string database_name = std::string(PKGDATADIR) + "/data/ligands-2016.db";
      bool do_help = false;
      bool do_rank = false;

      const char *optstr = "m:d:v:h";
      struct option long_options[] = {
	 {"metric",        1, 0, 0},
	 {"database",      1, 0, 0},
	 {"value",         1, 0, 0},
	 {"update-resolutions", 1, 0, 0},
	 {"rank",          0, 0, 0},
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

	       if (arg_str == "rank") {
		  do_rank = true;
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

	 if (do_rank) {

	    ranks(database_name);
	    
	 } else {

	    // normal usage

	    try {
	       if (! value_str.empty()) { 
		  double dc = coot::util::string_to_float(value_str);
	       
		  coot::ligand_metrics lm(database_name);
		  bool reverse_order = false;
		  std::pair<int, int> idx_pair = lm.get_index(dc, metric_str, reverse_order);
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
   }
#endif // USE_SQLITE3 
   return status;
}
