
#ifndef ROTAMER_CONTAINER_HH
#define ROTAMER_CONTAINER_HH

#include <clipper/core/coords.h>
#include "utils/colour-holder.hh"

#include "ligand/rotamer.hh"

class rotamer_markup_container_t {

public:
   coot::residue_spec_t spec;
   clipper::Coord_orth pos;
   coot::colour_holder col;
   coot::rotamer_probability_info_t rpi;

   rotamer_markup_container_t(const coot::residue_spec_t &spec_in,
                              const clipper::Coord_orth &pos_in,
			      const coot::colour_holder &col_in,
                              const coot::rotamer_probability_info_t &prob_in) :
      spec(spec_in), pos(pos_in), col(col_in), rpi(prob_in) {}
   rotamer_markup_container_t() {}
};

#endif // ROTAMER_CONTAINER_HH

