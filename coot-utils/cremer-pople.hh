#ifndef COOT_UTILS_CREMER_POPLE_HH
#define COOT_UTILS_CREMER_POPLE_HH

#include <vector>
#include <clipper/core/coords.h>

namespace coot {

   // Cremer-Pople puckering parameters for a 5- or 6-membered ring.
   //
   // The caller supplies the ring atom positions ALREADY IN RING ORDER. That
   // ordering is not a labelling convenience: rotating the start atom by one
   // position flips the sign of q3 and so swaps 4C1 and 1C4. Choosing the
   // ordering is the caller's job (in this project, cod_rings.py derives it
   // from the acedrg atom types).
   //
   // up_reference, when non-null, selects which ring face is "up": the mean-plane
   // normal is flipped if necessary so that (up_reference - centroid) . normal > 0.
   // Pass the position of a ring substituent to get a chirality-aware face. When
   // null, the normal follows the right-hand rule implied by the given ordering.
   class cremer_pople_t {
   public:
      cremer_pople_t(const std::vector<clipper::Coord_orth> &ring_positions,
                     const clipper::Coord_orth *up_reference = nullptr);

      bool filled;             // false if the ring size was not 5 or 6
      unsigned int n_ring;     // 5 or 6
      double Q;                // total puckering amplitude (A)
      double theta;            // radians, [0, pi]; NaN when n_ring == 5
      double phi;              // radians, [0, 2pi)
      double q2;
      double q3;               // 0.0 when n_ring == 5 (no m=3 term exists)
      std::vector<double> z;   // per-atom displacement from the CP mean plane
      clipper::Coord_orth normal;
   };

}

#endif // COOT_UTILS_CREMER_POPLE_HH
