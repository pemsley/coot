#ifndef API_RIGID_BODY_FIT_HH
#define API_RIGID_BODY_FIT_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>
#include "mini-mol/mini-mol.hh"

namespace coot {
   namespace api {
      void rigid_body_fit(mmdb::Manager *mol, int udd_atom_selection_fitting_atoms, const clipper::Xmap<float> &xmap);

      // above function calls this
      minimol::molecule
      rigid_body_fit_inner(const minimol::molecule &mol_without_moving_atoms,
                           const minimol::molecule &mol_for_moving_atoms,
                           const clipper::Xmap<float> &xmap);
   }
}

#endif // API_RIGID_BODY_REFINE_HH
