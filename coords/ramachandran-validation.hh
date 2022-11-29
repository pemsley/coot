#ifndef RAMACHANDRAN_VALIDATION_HH
#define RAMACHANDRAN_VALIDATION_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "rama-plot-phi-psi.hh"
#include "phi-psi-prob.hh"

namespace coot {

   std::vector<phi_psi_prob_t> ramachandran_validation(mmdb::Manager *mol, const ramachandrans_container_t &rc);

}


#endif // RAMACHANDRAN_VALIDATION_HH
