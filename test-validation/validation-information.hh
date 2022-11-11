#ifndef VALIDATION_INFORMATION_HH
#define VALIDATION_INFORMATION_HH

#include <vector>
#include "residue-validation-information.hh"

namespace coot {

   class chain_validation_information_t {
   public:
      std::string chain_id;
      std::vector<residue_validation_information_t> rviv;
      explicit chain_validation_information_t(const std::string &chain_id_in) : chain_id(chain_id_in) {}
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi) {
         rviv.push_back(rvi);
      }
   };

   class validation_information_min_max_t {
   public:
      bool is_set;
      double min;
      double max;
      validation_information_min_max_t() : is_set(false), min(0), max(0) {}
      validation_information_min_max_t(const double &min_in, const double &max_in) : is_set(true), min(min_in), max(max_in) {}
   };

   enum graph_data_type { UNSET, DISTORTION, ENERGY, PROBABILITY, CORRELATION, LOG_PROBABILITY };

   class validation_information_t {
   public:
      std::string name;
      // std::string type;
      enum graph_data_type type;
      validation_information_min_max_t min_max;
      std::vector<chain_validation_information_t> cviv;

      validation_information_t() : type(UNSET), min_max(validation_information_min_max_t()) {}

      validation_information_t(graph_data_type gdt, const validation_information_min_max_t &min_max_in) : type(gdt), min_max(min_max_in) {}

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
