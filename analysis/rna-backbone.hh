#ifndef COOT_UTILS_RNA_BACKBONE_HH
#define COOT_UTILS_RNA_BACKBONE_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "analysis/stats.hh"

namespace coot {

   class rna_backbone_t {
   public:
      stats::single base_samples;
      stats::single backbone_samples;

      rna_backbone_t(mmdb::Manager *mol, const clipper::Xmap<float> &xmap);
      void scan_all();
   };

}


#endif // COOT_UTILS_RNA_BACKBONE_HH
