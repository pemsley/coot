
#include <iostream>
#include "srs-interface.hh"
#include <stdlib.h>

#include "utils/coot-utils.hh"

// use environment variables or ccp4/prefix-dir fall-back
std::string
coot::get_srs_dir() {

   std::string dir;
   const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
   const char *d2 = getenv("CCP4");

   if (d1) {
      if (file_exists(d1))
	 dir = d1;
   } else {
      if (d2) {
	 std::string dir_a = util::append_dir_dir(d2, "share");
	 std::string dir_b = util::append_dir_dir(dir_a, "ccp4srs");
	 if (file_exists(dir_b))
	    dir = dir_b;
      }
   }
   
   if (dir.length()) {
      std::cout << "INFO:: CCP4SRS::loadIndex from dir: " << dir << std::endl;
   }

   return dir;
}
