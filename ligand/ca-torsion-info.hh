
#ifndef CA_TORSION_INFO_HH
#define CA_TORSION_INFO_HH

#include "angles/AngleInfo.h"
#include "mini-mol/mini-mol.hh"

namespace coot {

   class CA_torsion_info_t {

      // this class works in radians

   public:
      CA_torsion_info_t(const AngleInfo &ai_in) : ai(ai_in) {
	 status = false;
      }
      const AngleInfo &ai;
      bool status;
      double angle;
      double torsion;
      // seqnum is the residue number of the first residue to be added (N)
      // if offset = 1, building forwards on the C-terminus:
      // Let's score residues N-2, N-1, N, N+1.
      // We could also score N-3, N-2, N-1, N - which is perhaps more interesting?
      std::pair<bool, double> score_fragment(const minimol::fragment &frag,
					     mmdb::Residue *res_p,
					     mmdb::Residue *res_prev_p,
					     int seqnum, int offset);
   };

}



#endif // CA_TORSION_INFO_HH
