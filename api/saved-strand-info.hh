#ifndef SAVED_STRAND_INFO_HH
#define SAVED_STRAND_INFO_HH

#include "geometry/residue-and-atom-specs.hh"

namespace coot {
   class saved_strand_info_t {
   public:
      coot::residue_spec_t start;
      coot::residue_spec_t end;
      int strand_idx;
      saved_strand_info_t(const coot::residue_spec_t &s, const coot::residue_spec_t &e, int idx) : start(s), end(e) {
         strand_idx = idx;
      }
   };

}


#endif // SAVED_STRAND_INFO_HH
