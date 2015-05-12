
#include "coot-utils/coot-least-squares.hh"

class crunch_model_t {

public:
   crunch_model_t(const coot::least_squares_fit &lsq_in,
		  const clipper::Coord_orth &co) {
      lsq = lsq_in;
      centre = co;
   }
   coot::least_squares_fit lsq;
   clipper::Coord_orth centre;
};
