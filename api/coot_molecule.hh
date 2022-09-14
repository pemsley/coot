#ifndef COOT_MOLECULE_HH
#define COOT_MOLECULE_HH

#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "geometry/residue-and-atom-specs.hh"

class coot_molecule_t {

   atom_selection_container_t atom_sel;
   clipper::Xmap<float> map;

public:
   coot_molecule_t() {}

   bool is_valid_model_molecule();
   bool is_valid_map_molecule();

   int flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf);

};

#endif // COOT_MOLECULE_HH
