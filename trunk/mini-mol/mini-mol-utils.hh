
#include "mini-mol.hh"

namespace coot { 

   // C-beta position return a pair, the first of which indicates if
   // the cbeta was properly positioned.
   std::pair<short int, clipper::Coord_orth> cbeta_position(const minimol::residue &res);

   // Carbonyl O position return a pair, the first of which indicates
   // if the O was properly positioned.  There are distance sanity
   // checks applied also.
   std::pair<short int, clipper::Coord_orth> o_position(const minimol::residue &res_with_CA_C, const minimol::residue &res_with_N);

}
