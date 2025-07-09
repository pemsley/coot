
#include <string>
#include <filesystem>
#include <clipper/ccp4/ccp4_mtz_io.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include "born.hh"

int main(int argc, char **argv) {

   int status = 0;
   if (argc > 3) {
      std::string mtz_file_name = argv[1];
      std::filesystem::path mtz_file_path(mtz_file_name);
      if (std::filesystem::exists(mtz_file_path)) {

         try {
            std::string f_col   = argv[2];
            std::string phi_col = argv[3];
            coot::born(mtz_file_name, f_col, phi_col, "born.map");
         }

         catch (const std::runtime_error &rte) {
            std::cout << "ERROR" << rte.what() << std::endl;
         }
      }
   }
   return status;
}
