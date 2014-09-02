
#include "rotamer.hh"
#include "coot-utils.hh"
#include <sqlite3.h>

void store_data_1d(sqlite3 *db, const coot::a_rotamer_table &t);


static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   std::cout << "in callback argc is " << argc << std::endl;
   for(int i=0; i<argc; i++){
      // printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      std::string col_name;
   }
   printf("\n");
   return 0;
}

double callback_probability(void *NotUsed, int argc, char **argv, char **azColName){

   double probability = -1;
   for(int i=0; i<argc; i++){
      // printf("probability %s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      double probability = coot::util::string_to_float(argv[i]);
      std::cout << "prob " << probability << std::endl;
   }
   printf("\n");
   return probability;
}

void test_read(const std::string &db_file_name) {

   sqlite3 *db;
   int rc = sqlite3_open(db_file_name.c_str(), &db);
   std::cout << "rc: open read " << rc << std::endl;
   char *zErrMsg = 0;
   std::string command;
   std::cout << "rc: read " << rc << std::endl;
//    for (unsigned int ich=0; ich<10; ich++) { 
//       double chi = ich + 0.5;
//       std::string command = "select probability from val where chi_1 = ";
//       command += coot::util::float_to_string(chi).substr(0,4);
//       command += ";";
//       std::cout << "   " << command << std::endl;
//       rc = sqlite3_exec(db, command.c_str(), callback_probability, 0, &zErrMsg);
//       // std::cout << "rc: select " << rc << std::endl;
//    }
   for (unsigned int ich=0; ich<10; ich++) { 
      double chi = ich + 0.5;
      std::string command = "select probability from val where chi_1 = ";
      command += coot::util::float_to_string(chi).substr(0,4);
      command += ";";
      std::cout << "   " << command << std::endl;
      sqlite3_stmt *stmt;
      const char **pzTail = NULL;
      sqlite3_prepare_v2(db, command.c_str(), command.length()+1, &stmt, pzTail);
      int r = sqlite3_step(stmt);
      if (r == SQLITE_ROW) {
	 double p = sqlite3_column_double(stmt, 0);
	 std::cout << "   p: " << p << std::endl;
      } 
   }
   
}

void store_data(sqlite3 *db, const coot::a_rotamer_table &t) {

   char *zErrMsg = 0;
   std::cout << " " << t.residue_name << std::endl;
   if (t.residue_name == "VAL") {
      store_data_1d(db, t);
   } 
}

void store_data_1d(sqlite3 *db, const coot::a_rotamer_table &t) {

      for (unsigned int ich=0; ich<t.pr_chi_1.size(); ich++) { 
	 double chi = ich + 0.5;
	 std::string s = "insert into ";
	 s += t.residue_name; // table name
	 s += " (chi_1, probability) values ";
	 s += "(";
	 s += coot::util::float_to_string(chi);
	 s += ",";
	 s += coot::util::float_to_string_using_dec_pl(t.pr_chi_1[ich], 8);
	 s += ");";
	 // std::cout << "   " << s << std::endl;
	 int rc = sqlite3_exec(db, s.c_str(), callback, 0, &zErrMsg);
	 // std::cout << "rc: " << rc << std::endl;
      }
}


int main(int argc, char **argv) {

   if (argc > 1) {
      std::string db_file_name = argv[1];
      std::string dir = "../../coot/trunk/rama-data";
      coot::rotamer_probability_tables tables;
      tables.fill_tables(dir);
      sqlite3 *db;
      int rc = sqlite3_open(db_file_name.c_str(), &db);
      std::cout << "rc: " << rc << std::endl;
      char *zErrMsg = 0;

      std::vector<std::string> commands;
      commands.push_back("CREATE TABLE VAL (chi_1 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE SER (chi_1 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE THR (chi_1 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE CYS (chi_1 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE PRO (chi_1 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE ASN (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE ASP (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE PHE (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE TYR (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE TRP (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE HIS (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE ILE (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE LEU (chi_1 DECIMAL(5.1), chi_2 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE MET (chi_1 INT, chi_2 INT, chi_3 INT, probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE MSE (chi_1 INT, chi_2 INT, chi_3 INT, probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE GLU (chi_1 INT, chi_2 INT, chi_3 DECIMAL(5.1), probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE GLN (chi_1 INT, chi_2 INT, chi_3 INT, probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE ARG (chi_1 INT, chi_2 INT, chi_3 INT, chi_4 INT, probability DECIMAL (10.8) );");
      commands.push_back("CREATE TABLE LYS (chi_1 INT, chi_2 INT, chi_3 INT, chi_4 INT, probability DECIMAL (10.8) );");
      
      for (unsigned int ic=0; ic<commands.size(); ic++)
	 rc = sqlite3_exec(db, commands[ic].c_str(), callback, 0, &zErrMsg);

      for (unsigned int irt=0; irt<tables.n_tables(); irt++) {
	 const coot::a_rotamer_table &t = tables[irt];
	 store_data(db, t);
      }
      sqlite3_close(db);
      test_read(db_file_name);
   }
   return 0;
} 
