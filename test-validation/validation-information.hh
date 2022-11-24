#ifndef VALIDATION_INFORMATION_HH
#define VALIDATION_INFORMATION_HH

#include <vector>
#include "residue-validation-information.hh"

namespace coot {

   enum class graph_data_type {
      /// Between 0 and 100
      Distortion, 
      /// Just normal auto scale
      Energy, 
      /// Values from 0.0 - 1.0
      Probability, 
      /// Between 1.0 and -0.2 (let's say)
      /// Hanging down
      Correlation, 
      /// Values from 0.0 - 1.0, (negative) log scale
      LogProbability
   };
   class chain_validation_information_t {
   public:
      std::string chain_id;
      /// what is this used for?
      std::string name;
      /// if this corresponds to `validation_information_t::type`, then I guess it should be removed
      /// (in order to have a single source of truth)
      std::string type;
      std::vector<residue_validation_information_t> rviv;
      explicit chain_validation_information_t(const std::string &chain_id_in);
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi);
   };

   class validation_information_t {
   public:
      std::string name;
      graph_data_type type;
      std::vector<chain_validation_information_t> cviv;
      unsigned int get_index_for_chain(const std::string &chain_id);
      void add_residue_valiation_informtion(const residue_validation_information_t &rvi, const std::string &chain_id);
   };

}

#endif // VALIDATION_INFORMATION_HH
