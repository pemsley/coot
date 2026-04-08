#ifndef COOT_UTILS_CREMER_POPLE_HH

#define COOT_UTILS_CREMER_POPLE_HH
#include <mmdb2/mmdb_manager.h>

namespace coot {

   class cremer_pople_t {
   public:
      cremer_pople_t(mmdb::Residue *residue_p);
      double theta;
   };

}

#endif // COOT_UTILS_CREMER_POPLE_HH

