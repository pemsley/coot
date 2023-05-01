
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef USE_SQLITE3

#include <sqlite3.h>

#include <GraphMol/MolOps.h> // needed?

#include "cod-atom-types.hh"
#include "utils/coot-utils.hh"
// #include "coot-utils/coot-coord-utils.hh" out of order now
#include "bond-table-record-t.hh"

class hybridization_info_t {
public:
   std::string string;
   bool order_switch_needed;
   hybridization_info_t(const std::string &s, bool order_switch_in) :
      string(s), order_switch_needed(order_switch_in) {}
};


// fill data_store with as a std::vector<double>
// 
int db_select_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0;
   std::vector<double> *v = static_cast<std::vector<double> * > (data_store);
   std::cout << "argc " << argc << std::endl;
   for (int i=0; i<argc; i++) {
      if (argv[i] != NULL) {
         std::cout << " .... " << argv[i] << std::endl;
         v->push_back(coot::util::string_to_float(argv[i]));
      } else {
         std::cout << "null argv " << i << "!" << std::endl;
      }
   }
   return status;
}

int db_select_index_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0;
   int *i_p = static_cast<int *> (data_store);
   if (argc == 1) {
      *i_p = coot::util::string_to_float(argv[0]);
   }
   return status;
}
int db_select_type_3s_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0; // return something other than 1 if you want the select to stop
   std::vector<int> *iv_p = static_cast<std::vector<int> *> (data_store);
   if (argc == 1) {
      iv_p->push_back(coot::util::string_to_int(argv[0]));
   }
   return status;
}

int db_select_bond_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0;
   cod::bond_table_record_t *btr_p = static_cast<cod::bond_table_record_t *> (data_store);
   // std::cout << "    in db_select_bond_callback() argc: " << argc << std::endl;
   if (false)
      for (int i=0; i<argc; i++)
         if (argv[i] != NULL)
            std::cout << " .... " << argv[i] << std::endl;

   if (argc == 3) {
      cod::atom_type_t dum_1;
      cod::atom_type_t dum_2;
      if (false)
         std::cout << "    in db_select_bond_callback() filling with "
                   << argv[0] << " " << argv[1] << " " << argv[2]
                   << std::endl;
      cod::bond_table_record_t new_btr(dum_1, dum_2,
                                       coot::util::string_to_double(argv[0]),
                                       coot::util::string_to_double(argv[1]),
                                       coot::util::string_to_int(argv[2]));
      new_btr.file_name = "from-db";
      if (btr_p->count > 0) {
         std::cout << "POSSIBLE ERROR:: bond was already set in db_select_bond_callback()\n"
                   << "        This doesn't happen often. Was " << *btr_p
                   << ",  new: " << new_btr << std::endl;
      }
      if (btr_p->count < new_btr.count)
         *btr_p = new_btr;
   } else {
      std::cout << "    in db_select_bond_callback() fail to set bond" << std::endl;
   }
   return status;
}

// type of data: std::vector<cod::bond_table_record_t> accumulate_bonds;
//

int db_multi_select_bond_callback(void *data_store, int argc, char **argv, char **azColName) {
   int status = 0;
   std::vector<cod::bond_table_record_t> *btrv_p = static_cast<std::vector<cod::bond_table_record_t> *> (data_store);

   if (argc == 3) {
      cod::atom_type_t dum_1;
      cod::atom_type_t dum_2;
      cod::bond_table_record_t new_btr(dum_1, dum_2,
                                       coot::util::string_to_double(argv[0]),
                                       coot::util::string_to_double(argv[1]),
                                       coot::util::string_to_int(argv[2]));
      new_btr.file_name = "from-db";
      btrv_p->push_back(new_btr);
   } else {
      std::cout << "    in db_multi_select_bond_callback() fail to set bond" << std::endl;
   }
   return status;
}


int db_select_multi_bond_callback(void *data_store, int argc, char **argv, char **azColName) {
   /// hmm
   return 0;
}


// return -1 on failure
int
get_atom_index(const cod::atom_type_t &cod_1, sqlite3 *db) {

   int idx = -1;

   if (db) {
      std::string cmd = "SELECT atom_index from COD_TYPE_4_INDICES WHERE level_4_atom_type = " +
         coot::util::single_quote(cod_1.level_4, "'");
      int index = -1; // unset
      void *data_pointer = static_cast<void *>(&index);
      char *zErrMsg = 0;
      int rc = sqlite3_exec(db, cmd.c_str(), db_select_index_callback, data_pointer, &zErrMsg);
      if (rc !=  SQLITE_OK) {
         if (zErrMsg) {
            std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
         } else { 
            std::cout << "ERROR: processing command " << cmd << std::endl;
            sqlite3_free(zErrMsg);
         }
      } else {
         std::cout << "    " << cod_1.level_4 << " " << index << std::endl;
         idx = index;
      }
   }
   return idx;
}

// needed for below function
std::string
make_get_bond_record_query_string_from_type_2_indices(const std::vector<int> &v1,
                                                      const std::vector<int> &v2) {

   std::string s;
   if ((v1.size() > 0) && (v2.size() > 0)) {
      s += "SELECT mean, sd, count from COD_TYPE_4_BONDS WHERE ";
      s += " ( ";
      unsigned int v_size = v1.size();
      if (v_size > 495) {
         std::cout << "    info:: capping index_1 atoms at 495" << std::endl;
         v_size = 495;
      }
      unsigned int last_i_index = v_size - 1;
      for (unsigned int i=0; i<v_size; i++) {
         s += " atom_index_1 = ";
         s += coot::util::int_to_string(v1[i]);
         if (i != last_i_index)
            s += " OR ";
      }
      s += " ) ";
      s += " AND ";
      s += " ( ";
      v_size = v2.size();
      if (v_size > 495) {
         std::cout << "info:: capping index_2 atoms at 495" << std::endl;
         v_size = 495;
      }
      last_i_index = v_size - 1;
      for (unsigned int i=0; i<v_size; i++) {
         s += " atom_index_2 = ";
         s += coot::util::int_to_string(v2[i]);
         if (i != last_i_index)
            s += " OR ";
      }
      s += " ) ";
   }
   return s;
}

cod::bond_table_record_t
get_bond_record_by_combination(const std::vector<cod::bond_table_record_t> &accumulated_bonds,
                               cod::bond_table_record_t::approximation_level_t al) {

   cod::bond_table_record_t btr;
   double sum = 0;
   double sum_var = 0;
   unsigned int n = 0;

   if (accumulated_bonds.size() > 0) {
      for (unsigned int i=0; i<accumulated_bonds.size(); i++) {
         const cod::bond_table_record_t &b = accumulated_bonds[i];
         // std::cout << "adding bond " << b.mean << " " << b.std_dev << " " << b.count << std::endl;
         sum += b.mean * b.count;
         sum_var += b.std_dev * b.std_dev * b.count;
         n += b.count;
      }
      double mean = sum/double(n);
      double std_dev = sqrt(sum_var);
      cod::atom_type_t dum_1;
      cod::atom_type_t dum_2;
      btr = cod::bond_table_record_t(dum_1, dum_1, mean, std_dev, n, al);
   } else {
      std::cout << "oops accumulated_bonds size was 0" << std::endl;
   }
   return btr;
}


unsigned int
sum_n_bonds(const std::vector<cod::bond_table_record_t> &accumulate_bonds) {
   unsigned int n = 0;
   for (unsigned int i=0; i<accumulate_bonds.size(); i++) {
      n += accumulate_bonds[i].count;
   }
   return n;
}


cod::bond_table_record_t
get_bond_record_try_type_2(const cod::atom_type_t &cod_1,
                           const cod::atom_type_t &cod_2,
                           int idx_1, int idx_2, unsigned int n_count_min, sqlite3 *db) {

   std::cout << "    try type 2" << std::endl;
   cod::atom_type_t dum_1;
   cod::atom_type_t dum_2;
   cod::bond_table_record_t btr(dum_1, dum_2, 0.0, 0.0, 0);
   if (db) {
      // Like the level-3, find the atom indices that have the save level_2 as idx_2

      std::string cmd_1 = "SELECT atom_index from COD_TYPE_4_INDICES WHERE level_2_atom_type = ";
      std::string cmd_2 = "SELECT atom_index from COD_TYPE_4_INDICES WHERE level_2_atom_type = ";
      cmd_1 += coot::util::single_quote(cod_1.level_2.extra_electron_type(), "'");
      cmd_2 += coot::util::single_quote(cod_2.level_2.extra_electron_type(), "'");
      cmd_1 += " AND hash_code = " + coot::util::int_to_string(cod_1.hash_value);
      cmd_2 += " AND hash_code = " + coot::util::int_to_string(cod_2.hash_value);
      // cmd_1 += coot::util::int_to_string(cod_1.hash_value);
      // cmd_2 += coot::util::int_to_string(cod_2.hash_value);
      // this doesn't select anything at the moment because atom_type_2 in the tables are not
      // what they used to be - now they are electron-based types
      // std::cout << "debug:: cmd_1: " << cmd_1 << std::endl;
      std::vector<int> cod_1_type_2_atom_indices;
      std::vector<int> cod_2_type_2_atom_indices;
      void *data_pointer_1 = static_cast<void *> (&cod_1_type_2_atom_indices);
      void *data_pointer_2 = static_cast<void *> (&cod_2_type_2_atom_indices);
      char *zErrMsg_1 = 0;
      char *zErrMsg_2 = 0;
      int rc_1 = sqlite3_exec(db, cmd_1.c_str(), db_select_type_3s_callback, data_pointer_1, &zErrMsg_1);
      int rc_2 = sqlite3_exec(db, cmd_2.c_str(), db_select_type_3s_callback, data_pointer_2, &zErrMsg_2);
      if (rc_1 !=  SQLITE_OK) {
         if (zErrMsg_1) {
            std::cout << "ERROR: processing command " << cmd_1 << " " << zErrMsg_1 << std::endl;
         } else { 
            std::cout << "ERROR: processing command " << cmd_1 << std::endl;
            sqlite3_free(zErrMsg_1);
         }
      } else {
         // not so much rigour with the second select
         if (rc_2 ==  SQLITE_OK) {
            // ran without error
            std::cout << "    cod_1 type_2_atom_indices size() is "
                      << cod_1_type_2_atom_indices.size() << std::endl;
            std::cout << "    cod_2 type_2_atom_indices size() is "
                      << cod_2_type_2_atom_indices.size() << std::endl;
            std::string cmd = make_get_bond_record_query_string_from_type_2_indices(cod_1_type_2_atom_indices,
                                                                                    cod_2_type_2_atom_indices);
            // std::cout << "qs: " << cmd << std::endl; // debugging
            std::vector<cod::bond_table_record_t> accumulate_bonds;
            void *data_pointer = static_cast<void *> (&accumulate_bonds);
            char *zErrMsg = 0;
            int rc = sqlite3_exec(db, cmd.c_str(), db_multi_select_bond_callback, data_pointer, &zErrMsg);
            if (rc !=  SQLITE_OK) {
               if (zErrMsg) {
                  std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
               } else {
                  std::cout << "ERROR: processing command " << cmd << std::endl;
                  sqlite3_free(zErrMsg);
               }
            } else {
               std::cout << "    get_bond_record_try_type_2(), here with accumulate_bonds records size: "
                         << accumulate_bonds.size() << std::endl;
               if (accumulate_bonds.size() > 0) {
                  unsigned int n = sum_n_bonds(accumulate_bonds);
                  if (n >= 4) {
                     cod::bond_table_record_t::approximation_level_t al = cod::bond_table_record_t::EXTRA_ELECTRON;
                     btr = get_bond_record_by_combination(accumulate_bonds, al);
                  }
               }
            }
         }
      }
   } else {
      std::cout << "No database" << std::endl;
   }
   return btr;
}

// we pass a vector and an int - and they can be either way round.
//
std::vector<cod::bond_table_record_t>
get_type_3_bonds_inner(const std::vector<int> &cod_1_type_3_atom_indices,
                       int idx_2, int one_or_two_is_vector,
                       hybridization_info_t hybridization_info,
                       std::string ring_bool_str,
                       sqlite3 *db) {

   std::vector<cod::bond_table_record_t> accumulate_bonds;
   for (unsigned int i=0; i<cod_1_type_3_atom_indices.size(); i++) {
      std::string cmd;
      if (one_or_two_is_vector == 1)
         cmd = "SELECT mean, sd, count from COD_TYPE_4_BONDS WHERE atom_index_1 = " +
            coot::util::int_to_string(cod_1_type_3_atom_indices[i]) + " AND atom_index_2 = " +
            coot::util::int_to_string(idx_2) + 
            " AND hybridization = " + coot::util::single_quote(hybridization_info.string, "'") +
            " AND ring = " + coot::util::single_quote(ring_bool_str, "'");
      if (one_or_two_is_vector == 2)
         cmd = "SELECT mean, sd, count from COD_TYPE_4_BONDS WHERE atom_index_1 = " +
            coot::util::int_to_string(idx_2) + " AND atom_index_2 = " +
            coot::util::int_to_string(cod_1_type_3_atom_indices[i]) +
            " AND hybridization = " + coot::util::single_quote(hybridization_info.string, "'") +
            " AND ring = " + coot::util::single_quote(ring_bool_str, "'");

      cod::atom_type_t dum_1;
      cod::atom_type_t dum_2;
      cod::bond_table_record_t btr(dum_1, dum_2, 0.0, 0.0, 0);
      void *data_pointer = static_cast<void *> (&btr);
      char *zErrMsg = 0;
      // this should definately not give a "bond was already set" message
      int rc = sqlite3_exec(db, cmd.c_str(), db_select_bond_callback, data_pointer, &zErrMsg);
      if (rc !=  SQLITE_OK) {
         if (zErrMsg) {
            std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
         } else {
            std::cout << "ERROR: processing command " << cmd << std::endl;
            sqlite3_free(zErrMsg);
         }
      } else {
         if (btr.count > 0) {
            accumulate_bonds.push_back(btr);
         }
      }
   }
   return accumulate_bonds;
}


cod::bond_table_record_t
get_bond_record_try_colon_types(const cod::atom_type_t &cod_1,
                                const cod::atom_type_t &cod_2,
                                hybridization_info_t hybridization_info,
                                std::string ring_bool_str,
                                sqlite3 *db) {
   cod::bond_table_record_t btr;

   if (db) {
      std::string cmd_1 = "SELECT atom_index from COD_TYPE_4_INDICES WHERE colon_degree_atom_type = ";
      std::string cmd_2 = "SELECT atom_index from COD_TYPE_4_INDICES WHERE colon_degree_atom_type = ";
      cmd_1 += coot::util::single_quote(cod_1.neighb_degrees_str(), "'");
      cmd_2 += coot::util::single_quote(cod_2.neighb_degrees_str(), "'");
      cmd_1 += " AND hash_code = " + coot::util::int_to_string(cod_1.hash_value);
      cmd_2 += " AND hash_code = " + coot::util::int_to_string(cod_2.hash_value);
      std::vector<int> cod_1_atom_indices;
      std::vector<int> cod_2_atom_indices;
      void *data_pointer_1 = static_cast<void *> (&cod_1_atom_indices);
      void *data_pointer_2 = static_cast<void *> (&cod_2_atom_indices);
      char *zErrMsg_1 = 0;
      char *zErrMsg_2 = 0;
      int rc_1 = sqlite3_exec(db, cmd_1.c_str(), db_select_type_3s_callback, data_pointer_1, &zErrMsg_1);
      int rc_2 = sqlite3_exec(db, cmd_2.c_str(), db_select_type_3s_callback, data_pointer_2, &zErrMsg_2);
      if (rc_1 ==  SQLITE_OK) {
         if (rc_2 ==  SQLITE_OK) {
            std::cout << "    Here with cod-1 size " << cod_1_atom_indices.size() << std::endl;
            std::cout << "    Here with cod-2 size " << cod_2_atom_indices.size() << std::endl;
            unsigned int v_size_1 = cod_1_atom_indices.size();
            unsigned int v_size_2 = cod_2_atom_indices.size();
            if (v_size_1 > 0) {
               if (v_size_2 > 0) {
                  if (v_size_1 > 495) {
                     std::cout << "    info:: capping index_1 atoms at 495" << std::endl;
                     v_size_1 = 495;
                  }
                  if (v_size_2 > 495) {
                     std::cout << "    info:: capping index_1 atoms at 495" << std::endl;
                     v_size_2 = 495;
                  }
                  unsigned int last_i_index_1 = v_size_1 - 1;
                  unsigned int last_i_index_2 = v_size_2 - 1;
                  std::string s;
                  s += "SELECT mean, sd, count from COD_TYPE_4_BONDS WHERE ";
                  s += " ( ";
                  for (unsigned int i=0; i<v_size_1; i++) {
                     s += " atom_index_1 = ";
                     s += coot::util::int_to_string(cod_1_atom_indices[i]);
                     if (i != last_i_index_1)
                        s += " OR ";
                  }
                  s += " ) ";
                  s += " AND ";
                  s += " ( ";
                  for (unsigned int i=0; i<v_size_2; i++) {
                     s += " atom_index_2 = ";
                     s += coot::util::int_to_string(cod_2_atom_indices[i]);
                     if (i != last_i_index_2)
                        s += " OR ";
                  }
                  s += " ) ";
                  std::vector<cod::bond_table_record_t> accumulate_bonds;
                  void *data_pointer = static_cast<void *> (&accumulate_bonds);
                  char *zErrMsg = 0;
                  int rc = sqlite3_exec(db, s.c_str(), db_multi_select_bond_callback, data_pointer, &zErrMsg);
                  if (rc !=  SQLITE_OK) {
                     if (zErrMsg) {
                        std::cout << "ERROR: processing command " << s << " " << zErrMsg << std::endl;
                     } else {
                        std::cout << "ERROR: processing command " << s << std::endl;
                        sqlite3_free(zErrMsg);
                     }
                  } else {
                     std::cout << "    get_bond_record_try_colon_types(), here with accumulate_bonds records size: "
                               << accumulate_bonds.size() << std::endl;
                     if (accumulate_bonds.size() > 0) {
                        unsigned int n = sum_n_bonds(accumulate_bonds);
                        if (n >= 4) {
                           cod::bond_table_record_t::approximation_level_t al = cod::bond_table_record_t::COLON_DEGREE;
                           btr = get_bond_record_by_combination(accumulate_bonds, al);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return btr;
}


cod::bond_table_record_t
get_bond_record_try_type_3(const cod::atom_type_t &cod_1,
                           const cod::atom_type_t &cod_2,
                           int idx_1, int idx_2,
                           hybridization_info_t hybridization_info,
                           std::string ring_bool_str,
                           unsigned int n_count_min, sqlite3 *db) {

   cod::bond_table_record_t btr;
   //
   bool done_2 = false;
   if (cod_1.level_4 == cod_1.level_3) {
      if (cod_2.level_4 == cod_2.level_3) {
         btr = get_bond_record_try_type_2(cod_1, cod_2, idx_1, idx_2, n_count_min, db);
         done_2 = true;
      }
   }

   if (! done_2) {
      // perhaps there is some hope that there are enough type_3s?
      //
      // how do we find them though?
      //
      // which rows in the COD_TYPE_4_INDICES have the same type_3 as cod_1?
      // (and similar thinking for cod_2)
      //
      std::string cmd_1 = "SELECT atom_index from COD_TYPE_4_INDICES WHERE level_3_atom_type = ";
      std::string cmd_2 = "SELECT atom_index from COD_TYPE_4_INDICES WHERE level_3_atom_type = ";
      cmd_1 += coot::util::single_quote(cod_1.level_3, "'");
      cmd_2 += coot::util::single_quote(cod_2.level_3, "'");
      std::vector<int> cod_1_type_3_atom_indices;
      std::vector<int> cod_2_type_3_atom_indices;
      void *data_pointer_1 = static_cast<void *> (&cod_1_type_3_atom_indices);
      void *data_pointer_2 = static_cast<void *> (&cod_2_type_3_atom_indices);
      char *zErrMsg_1 = 0;
      char *zErrMsg_2 = 0;
      int rc_1 = sqlite3_exec(db, cmd_1.c_str(), db_select_type_3s_callback, data_pointer_1, &zErrMsg_1);
      int rc_2 = sqlite3_exec(db, cmd_2.c_str(), db_select_type_3s_callback, data_pointer_2, &zErrMsg_2);
      if (rc_1 !=  SQLITE_OK) {
         if (zErrMsg_1) {
            std::cout << "ERROR: processing command " << cmd_1 << " " << zErrMsg_1 << std::endl;
         } else { 
            std::cout << "ERROR: processing command " << cmd_1 << std::endl;
            sqlite3_free(zErrMsg_1);
         }
      } else {
         // not so much rigour with the second select
         if (rc_2 ==  SQLITE_OK) {
            // ran without error
            std::cout << "    now cod_1 type_3_atom_indices size() is "
                      << cod_1_type_3_atom_indices.size() << std::endl;
            std::cout << "    now cod_2 type_3_atom_indices size() is "
                      << cod_2_type_3_atom_indices.size() << std::endl;


            std::vector<cod::bond_table_record_t> accumulate_bonds_all;

            // 1 or 2 marks which atom has the vector of indicies
            //
            if (cod_1_type_3_atom_indices.size() > 0) {
               std::vector<cod::bond_table_record_t> accumulate_bonds_1 =
                  get_type_3_bonds_inner(cod_1_type_3_atom_indices, idx_2, 1, hybridization_info, ring_bool_str, db);

               std::vector<cod::bond_table_record_t> accumulate_bonds_2 =
                  get_type_3_bonds_inner(cod_2_type_3_atom_indices, idx_1, 2, hybridization_info, ring_bool_str, db);

               std::cout << "    here with " << accumulate_bonds_1.size()
                         << " and " << accumulate_bonds_2.size()
                         << " accumulated bonds" << std::endl;

               accumulate_bonds_all = accumulate_bonds_1;
               accumulate_bonds_all.insert(accumulate_bonds_all.end(),
                                           accumulate_bonds_2.begin(),
                                           accumulate_bonds_2.end());
               unsigned int n_bond_count = 0;
               for (unsigned int i=0; i<accumulate_bonds_all.size(); i++)
                  n_bond_count += accumulate_bonds_all[i].count;
               if (n_bond_count < n_count_min) {
                  std::cout << "    try level 2 generalization " << std::endl;
                  btr = get_bond_record_try_type_2(cod_1, cod_2, idx_1, idx_2, n_count_min, db);
               } else {
                  cod::bond_table_record_t::approximation_level_t al = cod::bond_table_record_t::MAIN_SECTION;
                  btr = get_bond_record_by_combination(accumulate_bonds_all, al);
               }
            }
         }
      }
   }

   return btr;
}


cod::bond_table_record_t
get_bond_record(const cod::atom_type_t &cod_1,
                const cod::atom_type_t &cod_2,
                int idx_1, int idx_2,
                hybridization_info_t hybridization_info,
                std::string ring_bool_str,
                unsigned int n_count_min, sqlite3 *db) {

   std::cout << "    in get_bond_record() " << idx_1 << " " << idx_2
             << " with hybridization_info: " << hybridization_info.string << std::endl;

   cod::atom_type_t dum_1;
   cod::atom_type_t dum_2;
   cod::bond_table_record_t btr(dum_1, dum_2, 0.0, 0.0, 0);
   if (db) {
      std::string cmd = "SELECT mean, sd, count from COD_TYPE_4_BONDS WHERE atom_index_1 = " +
         coot::util::int_to_string(idx_1) + " AND atom_index_2 = " + coot::util::int_to_string(idx_2) +
         " AND hybridization = " + coot::util::single_quote(hybridization_info.string, "'") +
         " AND ring = " + coot::util::single_quote(ring_bool_str, "'");

      void *data_pointer = static_cast<void *> (&btr);
      char *zErrMsg = 0;
      // this should not give a "bond was already set" message - btr is fresh
      int rc = sqlite3_exec(db, cmd.c_str(), db_select_bond_callback, data_pointer, &zErrMsg);
      if (rc !=  SQLITE_OK) {
         if (zErrMsg) {
            std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
         } else {
            std::cout << "ERROR: processing command " << cmd << std::endl;
            sqlite3_free(zErrMsg);
         }
      } else {
         // std::cout << "    in get_bond_record(): btr: " << btr << std::endl;
         if (btr.count < n_count_min) {
            std::cout << "    needs generalization - found (only) counts: " << btr.count << std::endl;

            // try this again with any hybridization
            cmd = "SELECT mean, sd, count from COD_TYPE_4_BONDS WHERE atom_index_1 = " +
               coot::util::int_to_string(idx_1) + " AND atom_index_2 = " + coot::util::int_to_string(idx_2) +
               " AND ring = " + coot::util::single_quote(ring_bool_str, "'");
            void *data_pointer = static_cast<void *> (&btr);
            char *zErrMsg = 0;
            // this might give a "bond was already set" message?
            // unlikely? btr can be set by any atom types pair
            int rc = sqlite3_exec(db, cmd.c_str(), db_select_bond_callback, data_pointer, &zErrMsg);
            if (rc !=  SQLITE_OK) {
               if (zErrMsg) {
                  std::cout << "ERROR: processing command " << cmd << " " << zErrMsg << std::endl;
               } else {
                  std::cout << "ERROR: processing command " << cmd << std::endl;
                  sqlite3_free(zErrMsg);
               }
            } else {
               if (btr.count < n_count_min) {
                  // this also will try type 2 if type 3 fail to find enough
                  btr = get_bond_record_try_type_3(cod_1, cod_2, idx_1, idx_2, hybridization_info, ring_bool_str,
                                                   n_count_min, db);
                  if (btr.count < n_count_min) {
                     std::cout << "    Found " << btr.count << " bonds - try the colon degree types" << std::endl;
                     btr = get_bond_record_try_colon_types(cod_1, cod_2, hybridization_info, ring_bool_str, db);
                  }
               }
            }
         } else {
            // golden
         }
      }
   } else {
      std::cout << "No database" << std::endl;
   }
   return btr;
}

// should be part of rdkit-interface or maybe cod-atom/bond.
//
hybridization_info_t
get_hybridization_info_string(RDKit::RWMol &rdkm, unsigned int mol_atom_idx_1, unsigned int mol_atom_idx_2) {
   std::string hif;

   // Acedrg tables contains: SP-NON_SP1 SP1_SP1 SP-NON_SP-NON SP1_SP2 SP-NON_SP2 SP-NON_SP3 SP2_SP3
   //                         SP3_SP3 SP2_SP2
   // which means that the atoms in the bond are order dependent
   //
   // SP > SP-NON > SP1 > SP2 > SP3 # perverse (or lexographical?)

   const RDKit::Atom *at_1 = rdkm[mol_atom_idx_1].get();
   const RDKit::Atom *at_2 = rdkm[mol_atom_idx_2].get();

   RDKit::Atom::HybridizationType h_1 = at_1->getHybridization();
   RDKit::Atom::HybridizationType h_2 = at_2->getHybridization();

   std::string part_1 = "SP-NON";
   std::string part_2 = "SP-NON";
   if (h_1 == RDKit::Atom::SP)    part_1 = "SP";
   if (h_1 == RDKit::Atom::SP2)   part_1 = "SP2";
   if (h_1 == RDKit::Atom::SP3)   part_1 = "SP3";
   if (h_1 == RDKit::Atom::SP3D)  part_1 = "SP3"; // what's this
   if (h_1 == RDKit::Atom::SP3D2) part_1 = "SP3"; // or this?
   if (h_2 == RDKit::Atom::SP)    part_2 = "SP";
   if (h_2 == RDKit::Atom::SP2)   part_2 = "SP2";
   if (h_2 == RDKit::Atom::SP3)   part_2 = "SP3";
   if (h_2 == RDKit::Atom::SP3D)  part_2 = "SP3";
   if (h_2 == RDKit::Atom::SP3D2) part_2 = "SP3";

   // The following is needed to match Acedrg bond lookup for PA-O* and PB-O*
   // in CCD ADP.
   //
   // Acedrg rules: if this is O-P and the O is charged, make the
   // hybridization on the O SP2 (rather than SP3)
   // Note that the O in P-O-P is SP3 - we shouldn't adjust that
   // (let's filter that out using the number of bonds)
   //
   unsigned int n_bonds_at_1 = 0;
   unsigned int n_bonds_at_2 = 0;
   RDKit::ROMol::OEDGE_ITER current, end;
   boost::tie(current, end) = rdkm.getAtomBonds(at_1);
   while (current != end) {
      RDKit::BOND_SPTR bond= rdkm[*current];
      current++;
      n_bonds_at_1++;
   }
   boost::tie(current, end) = rdkm.getAtomBonds(at_2);
   while (current != end) {
      RDKit::BOND_SPTR bond= rdkm[*current];
      current++;
      n_bonds_at_2++;
   }
   
   if (at_1->getAtomicNum() == 8) {
      if (at_2->getAtomicNum() == 15) {
         if (n_bonds_at_1 == 1)
            part_1 = "SP2";
      }
   }
   if (at_2->getAtomicNum() == 8) {
      if (at_1->getAtomicNum() == 15) {
         if (n_bonds_at_2 == 1)
            part_2 = "SP2";
      }
   }

   bool order_switch_needed = false;
   if (part_1 > part_2)
      order_switch_needed = true;
   hif = part_1 + "_" + part_2;
   if (order_switch_needed)
      hif = part_2 + "_" + part_1;
   hybridization_info_t hi(hif, order_switch_needed);
   return hi;
}

// ring atom info are in ints, not unsigned ints.
std::string
get_ring_info_bool_str(RDKit::RingInfo* ring_info_p, int idx_1, int idx_2) {

   std::string ri = "N";
   if (ring_info_p) {
      unsigned int n_rings = ring_info_p->numRings();
      std::vector<std::vector<int> > atom_rings = ring_info_p->atomRings();
      for (unsigned int i_ring=0; i_ring<n_rings; i_ring++) {
         unsigned int n_found = 0;
         std::vector<int> ring_atom_indices = atom_rings[i_ring];
         for (unsigned int j=0; j<ring_atom_indices.size(); j++) {
            if (ring_atom_indices[j] == idx_1) n_found++;
            if (ring_atom_indices[j] == idx_2) n_found++;
         }
         if (n_found == 2) {
            ri = "Y";
            break;
         }
      }
   } else {
      std::cout << "ERROR:: No ring info!" << std::endl;
   }
   return ri;
}

void debug_atom_types(const RDKit::RWMol &rdkm,
                      const std::vector<cod::atom_type_t> &v,
                      const coot::dictionary_residue_restraints_t &dict) {
   
   if (v.size() != dict.number_of_atoms()) {
      std::cout << "debug_atom_types(): mismatch types " << v.size() << " "
                << dict.number_of_atoms() << std::endl;
   } else {
      for (unsigned int i=0; i<v.size(); i++) {
         std::string name;
         RDKit::Atom *at = rdkm[i].get();
         at->getProp("name", name);
         std::string ee_type = v[i].level_2.extra_electron_type();
         std::cout << " types for " << name
                   << " type-4 " << v[i].level_4
                   << " nee " << v[i].level_2.n_extra_electrons()
                   << " std-level_2: " << std::setw(20) << v[i].level_2.string()
                   << " ee: " << ee_type
                   << " colon-hybrid: " << v[i].neighb_degrees_str()
                   << " hash: " << v[i].hash_value
                   << std::endl;
      }
   }
}

void validate_bonds(mmdb::Residue *residue_p,
                    const coot::dictionary_residue_restraints_t &dict,
                    RDKit::RWMol &rdkm,
                    int n_added,
                    sqlite3 *db) {

   unsigned int n_count_min = 5;
   cod::atom_types_t t;
   std::vector<cod::atom_type_t> v = t.get_cod_atom_types(rdkm); // for every atom
   debug_atom_types(rdkm, v, dict);

   if (v.size() != dict.number_of_atoms()) {
      std::cout << "ERROR:: mismatched size of types " << v.size() << " " << dict.number_of_atoms() << std::endl;
   } else {

      RDKit::RingInfo* ring_info_p = rdkm.getRingInfo();

      // is this the correct type for the data? c.f. bond_record_container_t?
      // std::map<cod::atom_type_t, cod::bond_table_record_t> bond_map;
      std::map<cod::atom_type_t, int> atom_indices;

      std::string atom_name_1;
      std::string atom_name_2;
      for (unsigned int i=0; i<dict.bond_restraint.size(); i++) {
         const coot::dict_bond_restraint_t &rest = dict.bond_restraint[i];
         if (! dict.is_hydrogen(rest.atom_id_1_4c())) {
            if (! dict.is_hydrogen(rest.atom_id_2_4c())) {
               mmdb::Atom **residue_atoms = 0;
               int n_residue_atoms;
               std::pair<unsigned int, mmdb::Atom *> atom_1(0, 0);
               std::pair<unsigned int, mmdb::Atom *> atom_2(0, 0);
               residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  mmdb::Atom *at = residue_atoms[iat];
                  std::string atom_name = at->name;
                  if (atom_name == rest.atom_id_1_4c()) {
                     atom_1.first = int(iat);
                     atom_1.second = at;
                     atom_name_1 = atom_name;
                  }
                  if (atom_name == rest.atom_id_2_4c()) {
                     atom_2.first = int(iat);
                     atom_2.second = at;
                     atom_name_2 = atom_name;
                  }
                  if (atom_1.second && atom_2.second) break;
               }
               if (atom_1.second && atom_2.second) {
                  double l_sqrd = (coot::co(atom_1.second)-coot::co(atom_2.second)).lengthsq();
                  double bl = sqrt(l_sqrd);
                  // std::cout << rest << " bl: " << bl << std::endl;
                  cod::atom_type_t cod_1 = v[atom_1.first];
                  cod::atom_type_t cod_2 = v[atom_2.first];

                  // make a std::map<std::pair<cod::atom_type_t, cod::atom_type_t>, cod::bond_table_record_t>
                  // so that we don't do the db search for bonds that we
                  // have looked up already.
                  // cod::bond_table_record_t btr = get_bond_table_record(cod_1, cod_2, db);

                  std::cout << "##### " << atom_name_1 << " " << atom_name_2 << " types:  l4: "
                            << cod_1.level_4 << "  l2: " << cod_1.level_2.extra_electron_type() << " and l4: "
                            << cod_2.level_4 << "  l2: " << cod_2.level_2.extra_electron_type() << std::endl;

                  hybridization_info_t hybridization_info =
                     get_hybridization_info_string(rdkm, atom_1.first, atom_2.first);

                  // check that these are not already in the index cache/map first
                  // need residue_atom index -> cod db atom index map
                  // std::map<int, int> cod_db_atom_index_map;
                  //
                  int idx_1 = get_atom_index(cod_1, db); // db atom type index, I mean.
                  int idx_2 = get_atom_index(cod_2, db);

                  if (idx_1 != -1) {
                     if (idx_2 != -1) {
                        std::string ring_bool_str =
                           get_ring_info_bool_str(ring_info_p, atom_1.first, atom_2.first);
                        if (hybridization_info.order_switch_needed) {
                           std::swap(cod_1, cod_2);
                           std::swap(idx_1, idx_2);
                           std::swap(atom_name_1, atom_name_2);
                        }
                        cod::bond_table_record_t btr = get_bond_record(cod_1, cod_2,
                                                                       idx_1, idx_2,
                                                                       hybridization_info,
                                                                       ring_bool_str,
                                                                       n_count_min, db);
                        std::cout << "    validate bond " << atom_name_1 << " " << atom_name_2 << " "
                                  << " with btr: " << btr << std::endl;
                        if (btr.count >= n_count_min) {
                           double delta = bl - btr.mean;
                           std::cout.setf(std::ios::fixed, std::ios::floatfield); // reset
                           std::cout << "    delta: "
                                     << atom_name_1 << " " << atom_name_2 << " "
                                     << std::setw(12) << delta << " A "
                                     << std::right << std::setprecision(4) << std::fixed
                                     << std::setw(12) << delta * 100.0 << " pm " << std::endl;
                           std::cout.setf(std::ios::fixed, std::ios::floatfield); // reset

                        }
                     } else {
                        std::cout << "No db atom type index for atom type " << cod_2.level_4 << std::endl;
                        // something like:
                        // cod::bond_table_record_t btr = get_bond_record_try_type_3(cod_1, cod_2, -1, -1, n_count_min, db);
                     }
                  } else {
                     std::cout << "No db atom type index for atom type " << cod_1.level_4 << std::endl;
                     // so... what do we do now?
                  }

               } else {
                  std::cout << "failed to find both atoms in " << rest << std::endl;
               }
            }
         }
      }
   }
}

void validate_bonds_for_residue_dict_db_file_name(mmdb::Residue *residue_p,
                                                  coot::dictionary_residue_restraints_t dict,
                                                  bool repronate_flag,
                                                  const std::string &db_file_name) {

   if (residue_p) {
      try {
         RDKit::RWMol rdkm = rdkit_mol(residue_p, dict);

         int n_added = 0;
         if (repronate_flag) {
            bool deloc_remaining_bonds = true;
            n_added += coot::remove_phosphate_hydrogens(&rdkm, deloc_remaining_bonds);
            n_added += coot::remove_sulphate_hydrogens(&rdkm, deloc_remaining_bonds);
            n_added += coot::remove_carboxylate_hydrogens(&rdkm, deloc_remaining_bonds);
            coot::charge_guanidinos(&rdkm);
            // restore ring info:
            std::vector<std::vector<int> > ring_info; // fill this
            RDKit::MolOps::findSSSR(rdkm, ring_info);
            // fix the dictionary to match
            dict.remove_phosphate_hydrogens(); // and charge the P
            dict.remove_sulphate_hydrogens();
            dict.remove_carboxylate_hydrogens();
         }

         if (coot::file_exists(db_file_name)) {
            sqlite3 *db;
            int rc = sqlite3_open(db_file_name.c_str(), &db);
            validate_bonds(residue_p, dict, rdkm, n_added, db);
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   } else {
      std::cout << "WARNING:: null residue " << std::endl;
   }
}

int main(int argc, char **argv) {

   int status = 0;
   bool reprotonate_flag = true; // make user-settable

   if (argc > 4) {
      std::string pdb_file_name = argv[1];
      std::string chain_id      = argv[2];
      std::string res_no_str    = argv[3];
      std::string db_file_name  = argv[4];
      try {
         int res_no = coot::util::string_to_int(res_no_str);
         if (! coot::file_exists(pdb_file_name)) {
            std::cout << "File not found " << pdb_file_name << std::endl;
         } else {
            mmdb::Manager *mol = new mmdb::Manager;
            mmdb::ERROR_CODE err = mol->ReadCoorFile(pdb_file_name.c_str());
            if (err) {
               std::cout << "ERROR on reading " << pdb_file_name << std::endl;
            } else {
               coot::residue_spec_t spec(chain_id, res_no, "");
               mmdb::Residue *residue_p = coot::util::get_residue(spec, mol);
               if (residue_p) {
                  // for now we will get a dictionary based on the residue
                  // for future we can read in a dictionary cif
                  //
                  // Recall that we need a fully hydrogenated ligand
                  // (otherwise the atom types are not correct).
                  coot::dictionary_residue_restraints_t dict(residue_p);

                  validate_bonds_for_residue_dict_db_file_name(residue_p, dict, 
                                                               reprotonate_flag, db_file_name);

               } else {
                  std::cout << "No residue found." << std::endl;
               }
            }
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "FAIL to get res_no from " << res_no_str << std::endl;
      }
   } else {

      if (argc == 3) {

         // the argument is a dictionary file name, get the residue from the dictionary
         std::string dict_file_name = argv[1];
         std::string   db_file_name = argv[2];
         if (! coot::file_exists(dict_file_name)) {
            std::cout << "WARNING:: file not found " << dict_file_name << std::endl;
         } else {
            coot::protein_geometry geom;
            coot::read_refmac_mon_lib_info_t rmit = geom.init_refmac_mon_lib(dict_file_name, 0);
            if (rmit.n_bonds == 0) {
               std::cout << "WARNING:: failed to find bonds in dictionary " << dict_file_name
                         << std::endl;
            } else {
               if (rmit.monomer_idx >= 0) {
                  const coot::dictionary_residue_restraints_t &dict =
                     geom.get_monomer_restraints(rmit.monomer_idx);
                  int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
                  mmdb::Residue *residue_p = dict.GetResidue(false, 20);
                  validate_bonds_for_residue_dict_db_file_name(residue_p, dict, 
                                                               reprotonate_flag, db_file_name);
               } else {
                  std::cout << "WARNING:: bad monomer index " << std::endl;
               }
            }
         }
      }
   }
   return status;
}

#else
#include <iostream>
int main(int argc, char **argv) {
   std::cout << "Not compiled with SQLITE3" << std::endl;
   return 1;
}

#endif // USE_SQLITE3

#else

#include <iostream>
int main(int argc, char **argv) {
   std::cout << "Not compiled for enhanced-ligands" << std::endl;
   return 1;
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
