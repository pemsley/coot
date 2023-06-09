#ifndef COOOT_API_SYMMETRY_INFO_HH
#define COOOT_API_SYMMETRY_INFO_HH

#include <vector>

#include "coords/mmdb-crystal.h"
#include "cell.hh"

namespace coot {

   class symmetry_info_t {
   public:
      std::vector<std::pair<symm_trans_t, Cell_Translation> > symm_trans;
      Cell cell;
      symmetry_info_t(const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans_in,
                      const Cell &cell_in) : symm_trans(symm_trans_in), cell(cell_in) {}
      symmetry_info_t() {}
   };

}


#endif // COOOT_API_SYMMETRY_INFO_HH
