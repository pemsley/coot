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
      explicit chain_validation_information_t(const std::string &chain_id_in) : chain_id(chain_id_in) {}
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi) {
         rviv.push_back(rvi);
      }
   };

   class validation_information_t {
   public:
      std::string name;
      std::string type;
      std::vector<chain_validation_information_t> cviv;
      unsigned int get_index_for_chain(const std::string &chain_id) {
         for (unsigned int i=0; i<cviv.size(); i++) {
            if (chain_id == cviv[i].chain_id)
               return i;
         }
         chain_validation_information_t cvi(chain_id);
         cviv.push_back(cvi);
         return cviv.size() -1;
      }
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi, const std::string &chain_id) {
         unsigned int idx = get_index_for_chain(chain_id);
         cviv[idx].add_residue_valiation_informtion(rvi);
      }
   };

}

#endif // VALIDATION_INFORMATION_HH
