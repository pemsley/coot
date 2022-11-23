#ifndef PHI_PSI_PROB_HH
#define PHI_PSI_PROB_HH

#include "coot-utils/coot-rama.hh"
#include "coords/Cartesian.h"
#include "coords/ramachandran-container.hh"

namespace coot {
   class phi_psi_prob_t {
      public:
         phi_psi_prob_t() {}
      phi_psi_prob_t(const util::phi_psi_t &pp, const Cartesian &pos, const ramachandrans_container_t &rama_container);
      util::phi_psi_t phi_psi;
      Cartesian position;
      double probability;
      bool is_allowed_flag;
      bool is_allowed() const { return is_allowed_flag; }
      std::string residue_name() const { return phi_psi.residue_name(); }
   };
}


#endif // PHI_PSI_PROB_HH

