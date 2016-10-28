
#ifndef SRS_INTERFACE_HH
#define SRS_INTERFACE_HH

#include <string>

#define MONOMER_DIR_STR "COOT_CCP4SRS_DIR"

namespace coot { 
   int get_min_match(const int &n1, const float similarity);
   std::string get_srs_dir(); // use environment variables or ccp4/prefix-dir fall-back
}

#endif // SRS_INTERFACE_HH
