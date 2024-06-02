#ifndef PLI_DOTS_REPRESENTATION_INFO_HH
#define PLI_DOTS_REPRESENTATION_INFO_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"
#include "solvent-exposure-difference-helper.hh"

namespace pli {

   // Is this for the 3D representation of the ligand? I think so.

   class dots_representation_info_t {

   public:
      dots_representation_info_t() {}
      std::vector<std::pair<coot::atom_spec_t, float> >
      solvent_accessibilities(mmdb::Residue *res_ref,
                              const std::vector<mmdb::Residue *> &filtered_residues);
      std::vector<solvent_exposure_difference_helper_t>
      solvent_exposure_differences(mmdb::Residue *res_ref, const std::vector<mmdb::Residue *> &residues) const;
   };

}


#endif // DOTS_REPRESENTATION_INFO_HH
