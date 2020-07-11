
#ifndef ROTAMER_CONTAINER_HH
#include "utils/colour-holder.hh"

class rotamer_markup_container_t {

public:
   clipper::Coord_orth pos;
   coot::colour_holder col;

   rotamer_markup_container_t(const clipper::Coord_orth &pos_in,
			      const coot::colour_holder &col_in) {
      pos = pos_in;
      col = col_in;
   }

   rotamer_markup_container_t() {}
};

#endif // ROTAMER_CONTAINER_HH

