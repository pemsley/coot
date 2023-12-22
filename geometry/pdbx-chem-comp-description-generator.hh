#ifndef PDBX_CHEM_COMP_DESCRIPTOR_HH
#define PDBX_CHEM_COMP_DESCRIPTOR_HH

#include <string>

namespace coot {

   class pdbx_chem_comp_description_generator_t {
   public:
      pdbx_chem_comp_description_generator_t() {}
      pdbx_chem_comp_description_generator_t(const std::string &pn, const std::string &pv, const std::string &d) :
         program_name(pn), program_version(pv), descriptor(d) {}
      std::string program_name;
      std::string program_version;
      std::string descriptor;
   };

}

#endif // PDBX_CHEM_COMP_DESCRIPTOR_HH
