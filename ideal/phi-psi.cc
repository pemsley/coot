
#include "phi-psi.hh"


std::ostream &
operator<<(std::ostream &s, const coot::phi_psi_t &pp) {
   s << pp.phi << " " << pp.psi;
   return s;
}

