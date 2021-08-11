
#ifndef SFCALC_GENMAP_HH
#define SFCALC_GENMAP_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/hkl_datatypes.h>
#include <clipper/core/xmap.h>

namespace coot {
   namespace util {
      void sfcalc_genmap(mmdb::Manager *mol,
                         const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                         const clipper::HKL_data<clipper::data32::Flag> &free,
                         clipper::Xmap<float> *xmap_p);

   }
}

#endif // SFCALC_GENMAP_HH

