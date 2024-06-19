
#ifndef PLI_SOLVENT_EXPOSURE_DIFFERENCE_HH
#define PLI_SOLVENT_EXPOSURE_DIFFERENCE_HH

#include "geometry/residue-and-atom-specs.hh"

namespace pli {

   // a trivial class to hold the residue and the solvent exposure,
   // including and not including the ligand.
   class solvent_exposure_difference_helper_t {
   public:
      coot::residue_spec_t res_spec;
      double exposure_fraction_holo;
      double exposure_fraction_apo;
      solvent_exposure_difference_helper_t(coot::residue_spec_t res_spec_in, double h, double a) :
         res_spec(res_spec_in), exposure_fraction_holo(h), exposure_fraction_apo(a) {}
   };
}

#endif // PLI_SOLVENT_EXPOSURE_DIFFERENCE_HH
