#ifndef RESIDUE_VALIDATION_INFORMATION_HH
#define RESIDUE_VALIDATION_INFORMATION_HH

#include <string>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {
   class residue_validation_information_t {
   public:
      residue_validation_information_t(const coot::residue_spec_t &rs,
                                       const coot::atom_spec_t &atom_spec_in,
                                       double distortion_in, const std::string &l);
      residue_spec_t residue_spec;
      atom_spec_t atom_spec;
      double distortion;
      
      // what is this for? The color is computed in the widget at draw-time. 
      std::string block_colour;
      std::string label;
   };
}

#endif // RESIDUE_VALIDATION_INFORMATION_HH
