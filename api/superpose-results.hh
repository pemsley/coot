
#include <string>
#include "validation-information.hh"

class superpose_results_t {
   public:
      std::string suppose_info;
      std::pair<std::string, std::string> alignment;
      coot::validation_information_t alignment_info;
};
