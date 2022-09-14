#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "geometry/residue-and-atom-specs.hh"

class molecules_container_t {

public:
   atom_selection_container_t atom_sel;
   clipper::Xmap<float> map;
   int flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf);

};

#endif // MOLECULES_CONTAINER_HH
