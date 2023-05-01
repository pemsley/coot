/* ligand/ligand-percentiles.cc
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

#ifdef USE_SQLITE3
#include <stdexcept>
#include <map>
#include <algorithm>
#include <iomanip>
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

bool string_double_sorter(const std::pair<std::string, double> &p1,
			  const std::pair<std::string, double> &p2) {

   bool r = false;
   return (p1.second<p2.second);
}

class ligand_index_score_t {

public:
   ligand_index_score_t(int idx_1_in, int idx_2_in, int idx_3_in, int idx_4_in, double score_in) {
      idx_1 = idx_1_in;
      idx_2 = idx_2_in;
      idx_3 = idx_3_in;
      idx_4 = idx_4_in;
      score = score_in;
   }
   int idx_1, idx_2, idx_3, idx_4;
   double score;
   static bool sort_compare_fn(const std::pair<std::string, ligand_index_score_t> & lis_1,
			       const std::pair<std::string, ligand_index_score_t> & lis_2) { 
      return (lis_1.second.score < lis_2.second.score);
   }
};


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

   std::cout << "# Found " << primary_keys.size() << " primary keys" << std::endl;

   std::vector<std::pair<std::string, ligand_index_score_t> > tri_ranks;
   std::map<std::string, std::vector<double> > values_store;

   for (unsigned int i=0; i<metric_names.size(); i++) {
      const std::string &metric_name = metric_names[i];
      std::vector<double> v = lm.get_values(metric_name);
      if (metric_name == "coot_diff_map_correlation")
	 for (unsigned int ii=0; ii<v.size(); ii++)
	    v[ii] = fabs(v[ii]);
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
	 if (metric_name == "coot_diff_map_correlation")
	    value = fabs(value);
	     
	 
	 bool reverse_order = false;
	 if (metric_name == "coot_diff_map_correlation") reverse_order = true;
	 if (metric_name == "bumps_1") reverse_order = true;
	 if (metric_name == "mogul_z_worst") reverse_order = true;
	 std::pair<int, int> idx_pair = lm.get_index(value, values_store[metric_name], reverse_order);

	 if (false)
	    std::cout << "DEBUG:: ranks() " << code << " " << metric_name << " " << value << "   "
		      << idx_pair.first << " / " << idx_pair.second << std::endl;

	 if (idx_pair.second != 0) {
	    rank_indices.push_back(idx_pair.first);
	 } else {
	    std::cout << "bad idx " << code << " " << metric_name << " has data size() "
		      << values_store[metric_name].size() << std::endl;
	 } 
      }

      if (rank_indices.size() == metric_names.size()) {
 	 double vi = 0;
 	 for (unsigned int ii=0; ii<metric_names.size(); ii++) {
 	    vi += rank_indices[ii] * rank_indices[ii];
 	 }
	 ligand_index_score_t lis(rank_indices[0], rank_indices[1], rank_indices[2], rank_indices[3], vi);
 	 tri_ranks.push_back(std::pair<std::string, ligand_index_score_t> (code, lis));
      }
   }

   std::sort(tri_ranks.begin(), tri_ranks.end(), ligand_index_score_t::sort_compare_fn);

   std::cout << "code rank direct-map-correl diff-map-correl bumps mogul_z_worst ";
   std::cout << "comp_id chain_id res_no\n";
   int n = primary_keys.size();
   double m = metric_names.size() * n * n;
   double rm = 1.0/m;
   std::vector<std::pair<std::string, ligand_index_score_t> >::const_iterator it;
   for (it = tri_ranks.begin(); it != tri_ranks.end(); it++) {
      std::pair<coot::residue_spec_t, std::string> residue_spec_and_type = lm.get_spec_and_type(it->first);
      std::cout << it->first << "  "
		<< std::setw(12) << rm * it->second.score << " "
		<< std::setw(5) << it->second.idx_1 << " "
		<< std::setw(5) << it->second.idx_2 << " "
		<< std::setw(5) << it->second.idx_3 << " "
		<< std::setw(5) << it->second.idx_4 << " "
		<< residue_spec_and_type.second << " "
		<< residue_spec_and_type.first.chain_id << " "
		<< residue_spec_and_type.first.res_no << "\n";
   }
   
}

#endif // USE_SQLITE3


int main(int argc, char **argv) {

   int status = 0;

#ifdef USE_SQLITE3

   if (argc > 1) {
      std::string metric_str = "density_correlation";
      std::string value_str;
      std::string database_name = coot::package_data_dir() + "/data/ligands-2016.db";
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
	      (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

	 switch(ch) {
	    
	 case 0:

	    if (coot_optarg) { 

	       std::string arg_str = long_options[option_index].name;

	       if (arg_str == "metric") {
		  metric_str = coot_optarg;
	       } 
	       if (arg_str == "database") {
		  database_name = coot_optarg;
	       } 
	       if (arg_str == "value") {
		  value_str = coot_optarg;
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
	    metric_str = coot_optarg;
	    break;
	    
	 case 'd':
	    database_name = coot_optarg;
	    break;
	    
	 case 'v':
	    value_str = coot_optarg;
	    break;
	    
	 case 'h':
	    do_help = true;
	    break;
	    
	 default:
	    std::cout << "default coot_optarg: " << coot_optarg << std::endl;
	    break;
	 }
      }

      if (do_help) {

	 std::cout << "Usage: " << argv[0] << "\n"
		   << "     --database <db-file-name>\n"
		   << "     --metric <metric-name>\n"
		   << "               e.g density_correlation, bumps_1, mogul_z_worst, diff_map_stat_1\n"
		   << "     --value <metric-name>\n"
		   << "     --update-database <resolutions-table-file-name>\n"
		   << "     --rank   rank PDB ligands\n"
		   << "     --help   this output\n";
	 
	 // e.g. ligand-percentiles --database ligands-2016.db --value 0.6
	 //                         --metric density_correlation

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
