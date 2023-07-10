
#include <string>
#include "validation-information.hh"

//! superposition results
class superpose_results_t {
   public:
      //! the JSON string with the superposition statistics
      std::string superpose_info;
      //! the superposition alignment
      std::pair<std::string, std::string> alignment;
      //! post-alignment distances between aligned residues
      //! encoded as validation information for the reference sequence
      std::vector<coot::validation_information_t> alignment_info_vec;
      std::vector<std::pair<coot::residue_validation_information_t, coot::residue_validation_information_t> > aligned_pairs;
};
