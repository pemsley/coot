/* ligand/ligands-db.cc
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
#include <math.h> 
#include <sqlite3.h>
#include "ligands-db.hh"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "utils/coot-utils.hh"
#include "ligands-db.hh"

coot::ligand_metrics::ligand_metrics(const std::string &db_file_name) {

   init();
   if (file_exists(db_file_name)) {
      int rc = sqlite3_open(db_file_name.c_str(), &db_);
   } else {
      std::cout << "WARNING:: File not found " << db_file_name << std::endl;
   } 
}

void
coot::ligand_metrics::init() {

   table_name = "LIGANDS";
   db_ = NULL;
}

std::vector<double>
coot::ligand_metrics::get_values(const std::string &col_name) const {

   std::vector<double> v;

   if (db_) { 
      std::string cmd = "select " + col_name + " from " + table_name + " ;";
      char *zErrMsg = 0;

      void *data_pointer = static_cast<void *>(&v);

      int rc = sqlite3_exec(db_, cmd.c_str(), db_select_callback, data_pointer, &zErrMsg);
      if (rc !=  SQLITE_OK) {
	 if (zErrMsg) { 
	    std::cout << "ERROR: processing command '" << cmd << "' "
		      << zErrMsg << std::endl;
	 } else { 
	    std::cout << "ERROR: processing command " << cmd << std::endl;
	    sqlite3_free(zErrMsg);
	 }
      }
   } else {
      std::cout << "invalid database" << std::endl;
   }

   return v;
}

// return high idx for good ligands
//
// return second of 0 on failure
// 
std::pair<int, int>
coot::ligand_metrics::get_index(double val, const std::string &col_name, bool low_is_good) const {

   int len = 0;
   int idx = 0;
   std::vector<double> v = get_values(col_name);

   if (v.size()) { 

      if (col_name == "coot_diff_map_correlation")
	 for (unsigned int i=0; i<v.size(); i++)
	    v[i] = fabs(v[i]);

      std::sort(v.begin(), v.end()); // sorted low to high

      if (low_is_good) { 
	 for (unsigned int i=0; i<v.size(); i++) { 
	    if (val <= v[i]) {
	       len = v.size();
	       idx = len -i;
	       break;
	    }
	 }
      } else {
	 bool found = false;
	 for (unsigned int i=0; i<v.size(); i++) { 
	    if (val < v[i]) {
	       idx = i;
	       len = v.size();
	       break;
	    }
	 }

	 if (! found) { 
	    // OK! was it the top ranked item?
	    //
	    if (val == v.back()) {
	       len = v.size();
	       idx = len;
	    } 
	 }
      }
   } else {
      std::cout << "No data to index" << std::endl;
   } 
   return std::pair<int, int>(idx, len);
}

std::pair<int, int>
coot::ligand_metrics::get_index(double val, const std::vector<double> &v, bool low_is_good) const {

   int len = 0;
   int idx = 0;

   if (v.size()) {
      if (low_is_good) { 
	 for (unsigned int i=0; i<v.size(); i++) { 
	    if (val <= v[i]) {
	       len = v.size();
	       idx = len -i;
	       break;
	    }
	 }
      } else {
	 bool found = false;
	 for (unsigned int i=0; i<v.size(); i++) { 
	    if (val < v[i]) {
	       idx = i;
	       len = v.size();
	       found = true;
	       break;
	    }
	 }

	 if (! found) { 
	    // OK! was it the top ranked item?
	    //
	    if (val == v.back()) {
	       len = v.size();
	       idx = len;
	    } 
	 }
      }
   } else {
      std::cout << "No data to index " << std::endl;
   }
   return std::pair<int, int>(idx, len);
}


std::pair<coot::residue_spec_t, std::string>
coot::ligand_metrics::get_spec_and_type(const std::string &accession_code) const {

   std::string table_name = "LIGANDS";
   std::string cmd = "SELECT chain_id,res_no,comp_id from " + table_name +
      " where accession_code = '" + accession_code + "' ;" ;

   std::pair<residue_spec_t, std::string> p;
   void *data_pointer = static_cast<void *>(&p);
   char *zErrMsg = 0;
   int rc = sqlite3_exec(db_, cmd.c_str(), db_select_spec_callback, data_pointer, &zErrMsg);
   if (rc !=  SQLITE_OK) {
      if (zErrMsg) { 
	 std::cout << "ERROR: processing command: " << cmd << " " << zErrMsg << std::endl;
      } else { 
	 std::cout << "ERROR when processing command: " << cmd << std::endl;
	 sqlite3_free(zErrMsg);
      }
   }

   return p;

}



void
coot::ligand_metrics::process_ligand_metrics_tab_line(const std::string &line,
						      sqlite3 *db) {

   std::string table_name = "LIGANDS";
   bool debug = false;
   std::vector<std::string> bits = coot::util::split_string_no_blanks(line, " ");

   
   if (false) {
      std::cout << bits.size() << " fields:" << std::endl;
      for (unsigned int i=0; i<bits.size(); i++) { 
	 std::cout << " (" << i << " '" << bits[i] << "')\n";
      }
      std::cout << std::endl;
   }

   // This is super-fragile!
   
   if (bits.size() == 34) {
      std::string quoted_chain_id      = bits[ 2];
      std::string resno_str            = bits[ 3];
      std::string quoted_ins_code      = bits[ 4]; 
      std::string quoted_comp_id       = bits[ 5];
      std::string correlation_str      = bits[ 7];
      std::string mogul_z_worst_str    = bits[ 9];
      std::string bumps_1_str          = bits[11];
      std::string bumps_2_str          = bits[12];
      std::string bumps_3_str          = bits[13];
      std::string bumps_4_str          = bits[14];
      std::string bumps_5_str          = bits[15];
      std::string diff_map_stat_1_str  = bits[17];
      std::string diff_map_stat_2_str  = bits[18];
      std::string diff_map_stat_3_str  = bits[19];
      std::string diff_map_stat_4_str  = bits[20];
      std::string diff_map_stat_5_str  = bits[21];
      std::string diff_map_stat_6_str  = bits[22];
      std::string diff_map_stat_7_str  = bits[23];
      std::string diff_map_stat_8_str  = bits[24];
      std::string diff_map_stat_9_str  = bits[25];
      std::string diff_map_stat_10_str = bits[26];
      std::string diff_map_stat_11_str = bits[27];
      std::string diff_map_stat_12_str = bits[28];
      std::string bfactor_stat_1_str   = bits[30];
      std::string bfactor_stat_2_str   = bits[31];
      std::string bfactor_stat_3_str   = bits[32];
      std::string bfactor_stat_4_str   = bits[33];

      std::string pdb_file_name = bits[1];
      std::vector<std::string> pdb_bits = util::split_string_no_blanks(pdb_file_name, "/");
      std::string pdb_tail = pdb_bits.back();
      std::string accession_code = util::name_sans_extension(pdb_tail);

      std::string chain_id = quoted_chain_id.substr(1, 1);
      int lc = quoted_comp_id.length();
      std::string comp_id = quoted_comp_id.substr(1, lc-2);

      try {

	 double corr = util::string_to_float(correlation_str);
	 if (corr <= 1.0) {
	    if (corr >= -1.0) {

	       // diff map stats:
	       // 1: correlation
	       // 2: var_x
	       // 3: var_y
	       // 4: n
	       // 5: sum_x
	       // 6: sum_y
	       // 7: D
	       // 8: D2 (based on mean/sd of the map at the ligand) 
	       // 9: map_mean
	       // 10: map_mean_at_ligand
	       // 11: map_sd
	       // 12: map_sd_at_ligand
	       // 13: "b-factor-metrics:"
	       // 14: median-ratio
	       // 15: median-ligand
	       // 16: median-env
	       // 17: ks-test-D

	       std::string cmd = "INSERT INTO " + table_name +
		  std::string("(accession_code, chain_id, res_no, comp_id, density_correlation, ") +
		  std::string("mogul_z_worst, bumps_1, coot_diff_map_correlation, coot_diff_map_n_grid, ") +
		  std::string("coot_diff_map_KS_1, coot_diff_map_KS_2, ") + 
		  std::string("coot_diff_map_mean, coot_diff_map_mean_at_ligand, ") + 
		  std::string("coot_diff_map_sd, coot_diff_map_sd_at_ligand,") +
		  std::string("b_factor_median_ratio, b_factor_median_ligand, b_factor_median_env, b_factor_KS_D") +
		  std::string(") ") +
		  std::string(" VALUES ") + 
		  std::string(" ( ") +
		  util::single_quote(accession_code, "'")    + ", " +
		  util::single_quote(chain_id, "'")          + ", " +
		  util::single_quote(resno_str, "'")         + ", " +
		  util::single_quote(comp_id, "'")           + ", " +
		  util::single_quote(correlation_str, "'")   + ", " +
		  util::single_quote(mogul_z_worst_str, "'") + ", " +
		  util::single_quote(bumps_1_str, "'")       + ", " +
		  //util::single_quote(diff_map_stat_1_str, "'") +
		  diff_map_stat_1_str  + ", " + 
		  diff_map_stat_4_str  + ", " + 
		  diff_map_stat_7_str  + ", " + 
		  diff_map_stat_8_str  + ", " + 
		  diff_map_stat_9_str  + ", " + 
		  diff_map_stat_10_str + ", " + 
		  diff_map_stat_11_str + ", " + 
		  diff_map_stat_12_str + ", " +
		  bfactor_stat_1_str   + ", " +
		  bfactor_stat_2_str   + ", " +
		  bfactor_stat_3_str   + ", " +
		  bfactor_stat_4_str   + 
		  std::string(");");
	 
	       char *zErrMsg = 0;
	       int rc = sqlite3_exec(db, cmd.c_str(), db_callback, 0, &zErrMsg);
	       if (true) { 
		  if (zErrMsg) { 
		     std::cout << "rc for " << cmd << " " << rc << " " << zErrMsg << std::endl;
		  } else {
		     // std::cout << "rc for " << cmd << " " << rc << " " << std::endl;
		  }
	       }
	    } else {
	       std::cout << "bad correl: " << line << std::endl;
	    } 
	 } else {
	    std::cout << "bad correl: " << line << std::endl;
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Error parsing " << line << " " << rte.what() << std::endl;
      }
   }
}

void
coot::ligand_metrics::parse_core_metrics(const std::string &input_file_name,
					 const std::string &db_file_name) {

   if (! file_exists(db_file_name)) {

      std::ifstream f(input_file_name.c_str());
      if (f) {
	 std::vector<std::string> lines;
	 std::string line;
	 while (std::getline(f, line)) {
	    lines.push_back(line);
	 }

	 sqlite3 *db = NULL;
	 if (! file_exists(db_file_name))
	    db = make_db(db_file_name); // database empty but with ligands table

	 if (db) {
	    char *zErrMsg = 0;
	    sqlite3_exec(db, "BEGIN", db_callback, 0, &zErrMsg);
	    for (unsigned int i=0; i<lines.size(); i++)
	       process_ligand_metrics_tab_line(lines[i], db);
	    sqlite3_exec(db, "END", db_callback, 0, &zErrMsg);
	 }
	 std::cout << "INFO:: database " << db_file_name << " created." << std::endl;
      } else {
	 std::cout << "WARNING:: metrics file " << input_file_name << " not found."
		   << std::endl;
      }
   } else {
      std::cout << "WARNING:: database " << db_file_name << " already exists - skipping action"
		<< std::endl;
   } 
} 

// return status, return pointer on success, null on fail
// 
sqlite3 *
coot::ligand_metrics::make_db(const std::string &db_file_name) const {

   sqlite3 *db = NULL;
   if (! file_exists(db_file_name)) {
      int rc = sqlite3_open(db_file_name.c_str(), &db);
      char *zErrMsg = 0;
      std::string command;

      if (rc == 0) {
	 // happy path
	 std::vector<std::string> commands;
	 std::string create_cmd =
	    std::string("CREATE TABLE LIGANDS (") +
	    std::string("accession_code TEXT PRIMARY KEY, ") +
	    std::string("date TEXT, year TEXT, ") +
	    std::string("header TEXT, ") +
	    std::string("nominal_resolution REAL, ") +
	    std::string("effective_resolution REAL, ") +
	    std::string("chain_id TEXT, ") +
	    std::string("res_no INT, ") +
	    std::string("comp_id TEXT, ") + 
            std::string("density_correlation REAL, ") +
            std::string("mogul_z_average REAL, " )+ 
            std::string("mogul_z_worst   REAL, ") +
            std::string("bumps_1 INT, ") +
            std::string("bumps_2 INT, ") +
            std::string("bumps_3 INT, ") +
            std::string("bumps_4 INT, ") +
            std::string("bumps_5 INT, ") +
            std::string("coot_diff_map_correlation  REAL, ") +
	    std::string("coot_diff_map_n_grid  REAL, ") +
	    std::string("coot_diff_map_KS_1  REAL, ") +
	    std::string("coot_diff_map_KS_2  REAL, ") +
	    std::string("coot_diff_map_mean  REAL, ") +
	    std::string("coot_diff_map_mean_at_ligand  REAL, ") +
	    std::string("coot_diff_map_sd  REAL, ") +
	    std::string("coot_diff_map_sd_at_ligand  REAL, ") +
	    std::string("b_factor_median_ratio REAL,") +
	    std::string("b_factor_median_ligand REAL,") +
	    std::string("b_factor_median_env REAL,") +
	    std::string("b_factor_KS_D REAL,") +
	    edstats_columns() + 
	    ");";

	 rc = sqlite3_exec(db, create_cmd.c_str(), db_callback, 0, &zErrMsg);
	 if (rc) { 
	    if (zErrMsg) { 
	       std::cout << "rc for " << create_cmd << " " << rc << " " << zErrMsg << std::endl;
	    } else {
	       std::cout << "rc for " << create_cmd << " " << rc << " " << std::endl;
	    }
	    // exit(1); no exit from libraries
	    db = NULL; // something bad happend
	 }

      }
   }
   return db;
}

std::string
coot::ligand_metrics::edstats_columns() const {

   std::string r;
   std::vector<std::pair<std::string, std::string> > v;

   // these column names have to be consistent with SQL naming convention.
   // - -> _minus_
   // _ -> _plus_
   //
   
   // v.push_back(std::pair<std::string, std::string> ("RT",   "REAL"));    
   // v.push_back(std::pair<std::string, std::string> ("CI",   "REAL"));    
   // v.push_back(std::pair<std::string, std::string> ("RN",   "REAL"));
   
   if (false) {
      // these columns don't make sense for ligands
      v.push_back(std::pair<std::string, std::string> ("BAm",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("NPm",  "INT"));   
      v.push_back(std::pair<std::string, std::string> ("Rm",   "REAL"));    
      v.push_back(std::pair<std::string, std::string> ("RGm",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("SRGm", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("CCSm", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("CCPm", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("ZCCm", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("ZOm",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("ZDm",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("ZD_minus_m", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("ZD_plus_m", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("BAs",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("NPs",  "INT"));   
      v.push_back(std::pair<std::string, std::string> ("Rs",   "REAL"));    
      v.push_back(std::pair<std::string, std::string> ("RGs",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("SRGs", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("CCSs", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("CCPs", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("ZCCs", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("ZOs",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("ZDs",  "REAL"));   
      v.push_back(std::pair<std::string, std::string> ("ZD_minus_s", "REAL"));  
      v.push_back(std::pair<std::string, std::string> ("ZD_plus_s", "REAL"));
   }

   v.push_back(std::pair<std::string, std::string> ("edstats_BAa",  "REAL"));   
   v.push_back(std::pair<std::string, std::string> ("edstats_NPa",  "INT"));   
   v.push_back(std::pair<std::string, std::string> ("edstats_Ra",   "REAL"));    
   v.push_back(std::pair<std::string, std::string> ("edstats_RGa",  "REAL"));   
   v.push_back(std::pair<std::string, std::string> ("edstats_SRGa", "REAL"));  
   v.push_back(std::pair<std::string, std::string> ("edstats_CCSa", "REAL"));  
   v.push_back(std::pair<std::string, std::string> ("edstats_CCPa", "REAL"));  
   v.push_back(std::pair<std::string, std::string> ("edstats_ZCCa", "REAL"));  
   v.push_back(std::pair<std::string, std::string> ("edstats_ZOa",  "REAL"));   
   v.push_back(std::pair<std::string, std::string> ("edstats_ZDa",  "REAL"));   
   v.push_back(std::pair<std::string, std::string> ("edstats_ZD_minus_a", "REAL"));  
   v.push_back(std::pair<std::string, std::string> ("edstats_ZD_plus_a", "REAL"));

   for (unsigned int i=0; i<v.size(); i++) { 
      r += v[i].first;
      r += " ";
      r += v[i].second;
      if (i < (v.size()-1))
	 r += ", ";
      else 
	 r += " ";
   }

   return r;
}


static int coot::db_callback(void *NotUsed, int argc, char **argv, char **azColName) {
   int status = 0;
   // std::cout << "in callback argc is " << argc << std::endl;
   for(int i=0; i<argc; i++){
      std::cout << " db_callback(): " << i << " " << argv[1];
   }
   std::cout << std::endl;
   return status;
}


// fill data_store with as a std::vector<double>
// 
static int coot::db_select_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0;
   std::vector<double> *v = (std::vector<double> *)(data_store);
   for(int i=0; i<argc; i++){
      if (argv[i] != NULL)
	 v->push_back(util::string_to_float(argv[i]));
      else
	 std::cout << "null argv " << i << "!" << std::endl;
   }
   return status;
}

// fill data_store with as a std::vector<double>
// 
static int coot::db_select_primary_key_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0;
   std::vector<std::string> *v = static_cast<std::vector<std::string> *>(data_store);
   for(int i=0; i<argc; i++){
      if (argv[i] != NULL)
	 v->push_back(argv[i]);
      else
	 std::cout << "null argv " << i << "!" << std::endl;
   }
   return status;
}

// fill data_store with as a std::vector<double>
// 
static int coot::db_select_spec_callback(void *data_store, int argc, char **argv, char **azColName) {

   int status = 0;
   std::pair<coot::residue_spec_t, std::string> *p =
      static_cast<std::pair<coot::residue_spec_t, std::string> *>(data_store);
   std::string comp_id;
   residue_spec_t spec;

   if (argc == 3) { 
      if (argv[0] != NULL) {
	 std::string chain_id = argv[0];
	 if (argv[1] != NULL) {
	    // an int?
	    std::string res_no_str = argv[1];
	    if (argv[2] != NULL) {
	       std::string comp_id = argv[2];
	       try {
		  int res_no = coot::util::string_to_int(res_no_str);
		  spec = residue_spec_t(chain_id, res_no, "");
		  *p = std::pair<residue_spec_t, std::string> (spec, comp_id);
	       }
	       catch (const std::exception &e) {
		  std::cout << "ERROR:: " << e.what() << std::endl;
	       } 
	    }
	 }
      }
   }
   return status;
}



bool
coot::ligand_metrics::update_resolutions_by_line(const std::string &line) {

   bool status = false;
   std::vector<std::string> bits = coot::util::split_string_no_blanks(line, " ");
   if (bits.size() == 3) {
      try {
	 const std::string &code = bits[0];
	 double res_nominal   = util::string_to_float(bits[1]);
	 double res_effective = util::string_to_float(bits[2]);

	 if (res_nominal > 0) { 
	    if (res_effective > 0) { 

	       std::string cmd = "UPDATE " + table_name + " ";
	       cmd += " set ";
	       cmd += "nominal_resolution = " + bits[1] + ", ";
	       cmd += "effective_resolution = " + bits[2] + " ";
	       cmd += " WHERE accession_code = ";
	       cmd += util::single_quote(code, "'");
	       cmd += ";";
	 
	       char *zErrMsg = 0;
	       int rc = sqlite3_exec(db_, cmd.c_str(), db_callback, 0, &zErrMsg);
	       if (rc !=  SQLITE_OK) {
		  if (zErrMsg) { 
		     std::cout << "ERROR: processing command: " << cmd << " " << zErrMsg << std::endl;
		  } else { 
		     std::cout << "ERROR when processing command: " << cmd << std::endl;
		     sqlite3_free(zErrMsg);
		  }
	       } else {
		  // was OK!
		  status = true;
	       }
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Error parsing " << line << " " << rte.what() << std::endl;
      }
   }
   return status;
}

void
coot::ligand_metrics::update_resolutions(const std::string &resolutions_table_file_name) {

   if (db_) {

      if (file_exists(resolutions_table_file_name)) {

	 std::ifstream f(resolutions_table_file_name.c_str());
	 if (f) {
	    std::vector<std::string> lines;
	    std::string line;
	    while (std::getline(f, line)) {
	       lines.push_back(line);
	    }
	    int n_success = 0;
	    char *zErrMsg = 0;
	    sqlite3_exec(db_, "BEGIN", db_callback, 0, &zErrMsg);
	    for (unsigned int i=0; i<lines.size(); i++) { 
	       bool status = update_resolutions_by_line(lines[i]);
	       if (status)
		  n_success++;
	    } 
	    sqlite3_exec(db_, "END", db_callback, 0, &zErrMsg);
	    std::cout << "INFO:: " << n_success << " records updated successfully"
		      << std::endl;
	 }
      }
   } 
}

void
coot::ligand_metrics::update_headers(const std::string &headers_file_name) {

   // this is the same as update_resolutions() execpt for the by_line function.
   // Hmm...

   if (db_) {

      if (file_exists(headers_file_name)) {

	 std::ifstream f(headers_file_name.c_str());
	 if (f) {
	    std::vector<std::string> lines;
	    std::string line;
	    while (std::getline(f, line)) {
	       lines.push_back(line);
	    }
	    int n_success = 0;
	    char *zErrMsg = 0;
	    sqlite3_exec(db_, "BEGIN", db_callback, 0, &zErrMsg);
	    for (unsigned int i=0; i<lines.size(); i++) { 
	       bool status = update_headers_by_line(lines[i]);
	       if (status)
		  n_success++;
	    } 
	    sqlite3_exec(db_, "END", db_callback, 0, &zErrMsg);
	    std::cout << "INFO:: " << n_success << " records updated successfully"
		      << " out of " << lines.size() << std::endl;
	 }
      } else {
	 std::cout << "WARNING:: input file \"" << headers_file_name << "\" not found"
		   << std::endl;
      } 
   } 
} 

bool
coot::ligand_metrics::update_headers_by_line(const std::string &line) {

   bool status = false;
   std::vector<std::string> bits = coot::util::split_string_no_blanks(line, " ");

   if (bits.size() >= 3) {
      try {
	 std::string up_code = bits.back();
	 std::string code = util::downcase(up_code);
	 std::string date = bits.end()[-2];
	 std::vector<std::string> year_bits = coot::util::split_string_no_blanks(date, "-");

	 if (false) { 
	    std::cout << "up_code: " << up_code << std::endl;
	    std::cout << "code: " << code << std::endl;
	    std::cout << "date: " << date << std::endl;
	    std::cout << "year_bits ";
	    for (unsigned int i=0; i<year_bits.size(); i++)
	       std::cout << "   " << year_bits[i];
	    std::cout << std::endl;
	 }

	 std::string molecule_type;
	 for (unsigned int i=1; i<(bits.size()-2); i++) {
	    if (i != 1)
	       molecule_type += " ";
	    std::string part = bits[i];
	    // Here remove single quotes e.g. as in " 3'-TERMINAL" from part
	    molecule_type += part;
	 }
	 
	 if (year_bits.size() == 3) { 
	    const std::string &year_str = year_bits.back();

	    int year = util::string_to_int(year_str);
	    if (year < 50)
	       year += 2000;
	    else
	       year += 1900;

	    std::string new_year_str = util::int_to_string(year);

	    std::string cmd = "UPDATE " + table_name + " ";
	    cmd += " set ";
	    cmd += "year = "   + new_year_str + ", ";
	    cmd += "date = "   + util::single_quote(date, "'") + ", ";
	    cmd += "header = " + util::single_quote(molecule_type, "'") + " ";
	    cmd += " WHERE accession_code = ";
	    cmd += util::single_quote(code, "'");
	    cmd += ";";

	    // std::cout << "command: " << cmd << std::endl;
	 
	    char *zErrMsg = 0;
	    int rc = sqlite3_exec(db_, cmd.c_str(), db_callback, 0, &zErrMsg);
	    if (rc !=  SQLITE_OK) {
	       if (zErrMsg) { 
		  std::cout << "ERROR: processing command: " << cmd << " " << zErrMsg << std::endl;
	       } else { 
		  std::cout << "ERROR when processing command: " << cmd << std::endl;
		  sqlite3_free(zErrMsg);
	       }
	    } else {
	       // was OK!
	       status = true;
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Error parsing " << line << " " << rte.what() << std::endl;
      }
   }
   return status;
}

void
coot::ligand_metrics::update_edstats_results(const std::string &edstats_file_name) {
   
   // this is the same as update_resolutions() execpt for the by_line function.

   if (db_) {

      if (file_exists(edstats_file_name)) {

	 std::ifstream f(edstats_file_name.c_str());
	 if (f) {
	    std::vector<std::string> lines;
	    std::string line;
	    while (std::getline(f, line)) {
	       lines.push_back(line);
	    }
	    int n_success = 0;
	    char *zErrMsg = 0;
	    sqlite3_exec(db_, "BEGIN", db_callback, 0, &zErrMsg);
	    for (unsigned int i=0; i<lines.size(); i++) { 
	       bool status = update_edstats_results_by_line(lines[i]);
	       if (status)
		  n_success++;
	    } 
	    sqlite3_exec(db_, "END", db_callback, 0, &zErrMsg);
	    std::cout << "INFO:: " << n_success << " records updated successfully"
		      << " out of " << lines.size() << std::endl;
	 }
      } else {
	 std::cout << "WARNING:: input file \"" << edstats_file_name << "\" not found"
		   << std::endl;
      }
   }
}
   
bool
coot::ligand_metrics::update_edstats_results_by_line(const std::string &line) {

   bool status = false;
   std::vector<std::string> bits = coot::util::split_string_no_blanks(line, " ");
   if (bits.size() >= 41) {
      try {

	 // BAa   all-atom: weighted average B-iso                                 
	 // NPa   all-atom: Real-space Zdiff score (RSZD) i.e. max(-RSZD-,RSZD+)   
	 // Ra    all-atom: Real space R-factor                                    
	 // RGa   all-atom: Generalized Real space R-factor?                       
	 // SRGa  all-atom: Standard uncertainty of the Generalized Real space R-fa
	 // CCSa  all-atom: Real Space sample correlation coefficient (RSCC)       
	 // CCPa  all-atom: Real Space population correlation coefficient (RSCC)   
	 // ZCCa  all-atom: Z-score of real-space sample correlation coefficient   
	 // ZOa   all-atom: Real-space Zobs score (RSZO)                           
	 // ZDa   all-atom: Real-space Zdiff score (RSZD) i.e. max(-RSZD-,RSZD+)   
	 // ZD-a  all-atom: Real-space Zdiff score for negative differences (RSZD-)
	 // ZD+a  all-atom: Real-space Zdiff score for positive differences (RSZD+)
	 
	 std::string code = bits[1];

	 if (code.length() == 4) {
	    std::string BAa  = bits[29];
	    std::string NPa  = bits[30];
	    std::string Ra   = bits[31];
	    std::string RGa  = bits[32];
	    std::string SRGa = bits[33];
	    std::string CCSa = bits[34];
	    std::string CCPa = bits[35];
	    std::string ZCCa = bits[36];
	    std::string ZOa  = bits[37];
	    std::string ZDa  = bits[38];
	    std::string ZD_minus_a  = bits[39];
	    std::string ZD_plus_a   = bits[40];

	    if (false)
	       std::cout << "debug:: parse: " << code << " " << BAa << " " << NPa << " "
			 << ZD_minus_a << " " << ZD_plus_a << std::endl;

	    std::string cmd = "UPDATE " + table_name + " ";
	    cmd += " set ";
	    cmd += "edstats_BAa = "   + BAa + ", ";
	    cmd += "edstats_NPa = "   + NPa + ", ";
	    cmd += "edstats_Ra = "    + Ra + ", ";
	    cmd += "edstats_RGa = "   + RGa + ", ";
	    cmd += "edstats_SRGa = "  + SRGa + ", ";
	    cmd += "edstats_CCSa = "  + CCSa + ", ";
	    cmd += "edstats_CCPa = "  + CCPa + ", ";
	    cmd += "edstats_ZCCa = "  + ZCCa + ", ";
	    cmd += "edstats_ZOa = "   + ZOa  + ", ";
	    cmd += "edstats_ZDa = "   + ZDa + ", ";
	    cmd += "edstats_ZD_minus_a = " + ZD_minus_a + ", ";
	    cmd += "edstats_ZD_plus_a = "  + ZD_plus_a + " ";
	    cmd += " WHERE accession_code = ";
	    cmd += util::single_quote(code, "'");
	    cmd += ";";

	    char *zErrMsg = 0;
	    int rc = sqlite3_exec(db_, cmd.c_str(), db_callback, 0, &zErrMsg);
	    if (rc !=  SQLITE_OK) {
	       if (zErrMsg) { 
		  std::cout << "ERROR: in processing\n"
			    << " line:    \"" << line << "\"\n"
			    << " command: " << cmd << "\n"
			    << " error:   " << zErrMsg << std::endl;
	       } else { 
		  std::cout << "ERROR when processing command: " << cmd << std::endl;
		  sqlite3_free(zErrMsg);
	       }
	    } else {
	       // was OK!
	       status = true;
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Error parsing " << line << " " << rte.what() << std::endl;
      }
   }
   return status;
}


double
coot::ligand_metrics::get_value(const std::string &accession_code,
				const std::string &metric_name,
				bool from_int) const {

   double val = -1;
   std::vector<double> v;

   if (db_) { 
      std::string cmd = "select " + metric_name + " from " + table_name +
	 " where accession_code = '" + accession_code + "' ;";
      char *zErrMsg = 0;

      void *data_pointer = static_cast<void *>(&v);

      // std::cout << "processing " << cmd << std::endl;

      int rc = sqlite3_exec(db_, cmd.c_str(), db_select_callback, data_pointer, &zErrMsg);
      if (rc !=  SQLITE_OK) {
	 if (zErrMsg) {
	    std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
	 } else { 
	    std::cout << "ERROR: processing command " << cmd << std::endl;
	    sqlite3_free(zErrMsg);
	 }
      }
   } else {
      std::cout << "invalid database" << std::endl;
   }

   if (v.size() == 1)
      val = v[0];

   return val;
}

std::vector<std::string>
coot::ligand_metrics::get_primary_keys() const {

   std::vector<std::string> v;

   if (db_) { 
      std::string cmd = "select accession_code from " + table_name + " ;";
      char *zErrMsg = 0;

      void *data_pointer = static_cast<void *>(&v);

      int rc = sqlite3_exec(db_, cmd.c_str(), db_select_primary_key_callback, data_pointer, &zErrMsg);
      if (rc !=  SQLITE_OK) {
	 if (zErrMsg) { 
	    std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
	 } else { 
	    std::cout << "ERROR: processing command " << cmd << std::endl;
	    sqlite3_free(zErrMsg);
	 }
      }
   } else {
      std::cout << "invalid database" << std::endl;
   }
   

   return v;
} 


// this file is only compiled if we have configured with sqlite3
#endif // USE_SQLITE3
