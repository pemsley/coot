
#ifndef PHI_PSI_T
#define PHI_PSI_T

#include <iosfwd>

namespace coot {

   // consolidate with phi_psi_pair (i.e. make that disappear)
   class phi_psi_t {
      // and now with tau!
   public:
      phi_psi_t(float a, float b) {
	 phi = a;
	 psi = b;
	 tau = -20; // unset value
      }
      phi_psi_t(float a, float b, float c) {
	 phi = a;
	 psi = b;
	 tau = c;
      }
      float phi;
      float psi;
      float tau;
      friend std::ostream &operator<<(std::ostream &s, const phi_psi_t &pp);
   };
   std::ostream &operator<<(std::ostream &s, const phi_psi_t &pp);

}

#endif // PHI_PSI_T
