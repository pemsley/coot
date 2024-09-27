
#include "mmdb2/mmdb_defs.h"
#include "mmdb2/mmdb_model.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <mmdb2/mmdb_manager.h>

int n_tests = 0;
static std::vector<std::pair<std::string, int> > test_results;

namespace mmcif_tests {

   void write_test_name(const std::string &test_name);
   int run_test(int (*test_func) (), const std::string &test_name);
   int run_tests(bool last_test_only);
   int read_pdb();
   enum status_t { BAD, GOOD };
}

void
mmcif_tests::write_test_name(const std::string &test_name) {

   std::ofstream f(".current-test");
   f << "\"" << test_name << "\"" << "\n";
   f.close();
}

int
mmcif_tests::run_test(int (*test_func) (), const std::string &test_name) {

   n_tests++;
   write_test_name(test_name);
   int status = test_func();
   std::string status_string = "FAIL: ";
   std::string uncol = "[m";
   std::string col = "[31m";
   if (status == 1) {
      status_string = "PASS: ";
      col = "[32m";
   }
   std::cout << status_string << std::setw(40) << std::left << test_name << col << " ⬤ " << uncol << std::endl;
   test_results.push_back(std::make_pair(test_name, status));

   return status;
}

int mmcif_tests::read_pdb() {

   int status = 0;
   mmdb::Manager *mol = new mmdb::Manager;
   // 6DGD for DNA bases
   // https://files.rcsb.org/download/6GDG.cif
   mmdb::ERROR_CODE read_status = mol->ReadCoorFile("6gdg.cif");
   std::cout << "TEST read_pdb() with read_status " << read_status << std::endl;
   if (read_status == mmdb::Error_NoError) {
      for (int imod = 1; imod <= mol->GetNumberOfModels(); imod++) {
        mmdb::Model *model_p = mol->GetModel(imod);
        if (model_p) {
          int n_links = model_p->GetNumberOfLinks();
          std::cout << "n_links: " << n_links << std::endl;
          for (int i_link = 0; i_link < n_links; i_link++) {
            mmdb::Link *link_p = model_p->GetLink(i_link);
            std::cout << "Link " << i_link << " " << link_p << std::endl;
          }
        }
      }
   }
   return status;
};

// return "all pass" status (1 for true)
int
mmcif_tests::run_tests(bool last_test_only) {

   if (! last_test_only) {
      run_test(read_pdb, "read_pdb");
   }
   int status = 1;
   for (const auto &t : test_results) {
      if (t.second != 1)
         status = 0;
   }
   return status;
}

