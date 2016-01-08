
#ifndef LIGANDS_DB_HH
#define LIGANDS_DB_HH

#ifdef USE_SQLITE3
#include <sqlite3.h>

#include <string>
#include <vector>
#include "coot-utils/residue-and-atom-specs.hh"

namespace coot {

   static int db_callback(void *NotUsed, int argc, char **argv, char **azColName);
   static int db_select_callback(void *NotUsed, int argc, char **argv, char **azColName);
   static int db_select_primary_key_callback(void *NotUsed, int argc, char **argv, char **azColName);
   static int db_select_spec_callback(void *NotUsed, int argc, char **argv, char **azColName);

   class ligand_metrics {
      void process_ligand_metrics_tab_line(const std::string &line, sqlite3 *db);
      sqlite3 *db_;
      std::string table_name; // "LIGANDS"

      // return status, return pointer on success, null on fail
      // 
      sqlite3 *make_db(const std::string &db_file_name) const;
      bool update_resolutions_by_line(const std::string &line);
      bool update_headers_by_line(const std::string &line);
      bool update_edstats_results_by_line(const std::string &edstats_file_name);
      std::string edstats_columns() const;
      void init();
      
   public:
      ligand_metrics() { init(); }
      ligand_metrics(const std::string &file_name); // set db_
      // make the database output file if it doesn't exist, otherwise add/replace
      void parse_core_metrics(const std::string &input_file_name,
			      const std::string &output_db_file_name);
      double get_value(const std::string &accession_code, const std::string &metric_name, bool from_int) const;
      std::vector<double> get_values(const std::string &col_name) const;
      // return 0,0 on failure (e.g. wrong col_name)
      // return high idx for good ligands (use this for singletons)
      std::pair<int, int> get_index(double val, const std::string &col_name, bool low_is_good) const;
      // pre-fetch the data and use this get_index() when we we want multiple calls.
      std::pair<int, int> get_index(double val, const std::vector<double> &v, bool low_is_good) const;
      void update_resolutions(const std::string &resolutions_table_file_name);
      void update_headers(const std::string &headers_file_name);
      void update_edstats_results(const std::string &edstats_file_name);
      std::vector<std::string> get_primary_keys() const;
      std::pair<residue_spec_t, std::string> get_spec_and_type(const std::string &accession_code) const;
   };
} 

#endif // USE_SQLITE3
#endif // LIGANDS_DB_HH
