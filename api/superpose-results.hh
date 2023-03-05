
#include <string>
#include "validation-information.hh"

//! superposition results
class superpose_results_t {
   public:
      //! the JSON string with the superposition statistics
      std::string suppose_info;
      //! the superposition alignment
      std::pair<std::string, std::string> alignment;
      //! post-alignment distances between aligned residues
      //! encoded as validation information.
      coot::validation_information_t alignment_info;
};
