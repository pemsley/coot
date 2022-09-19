#ifndef COOT_MOLECULE_HH
#define COOT_MOLECULE_HH

#include <clipper/core/xmap.h>
#include "coot-utils/atom-selection-container.hh"
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class molecule_t {

   public:
      atom_selection_container_t atom_sel;
      molecule_t() {}
      explicit molecule_t(atom_selection_container_t asc) : atom_sel(asc) {}
      clipper::Xmap<float> xmap; // public because the filling function needs access

      // utils

      bool is_valid_model_molecule();
      bool is_valid_map_molecule();
      std::pair<bool, coot::residue_spec_t> cid_to_residue_spec(const std::string &cid);

      // functions

      int flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf);
   };
}


#endif // COOT_MOLECULE_HH
