#ifndef VALIDATION_INFORMATION_HH
#define VALIDATION_INFORMATION_HH

#include <vector>
#include "residue-validation-information.hh"

namespace coot {

   class chain_validation_information_t {
   public:
      std::string chain_id;
      std::string name;
      std::string type;
      std::vector<residue_validation_information_t> rviv;
      explicit chain_validation_information_t(const std::string &chain_id_in);
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi);
   };

   class validation_information_t {
   public:
      std::string name;
      std::string type;
      std::vector<chain_validation_information_t> cviv;
      unsigned int get_index_for_chain(const std::string &chain_id);
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi, const std::string &chain_id);
   };

}

#endif // VALIDATION_INFORMATION_HH
