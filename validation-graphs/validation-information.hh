#ifndef VALIDATION_INFORMATION_HH
#define VALIDATION_INFORMATION_HH

#include <vector>
#include "residue-validation-information.hh"

namespace coot {

   enum class graph_data_type {
      // not yet set
      UNSET,
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
      LogProbability,
      /// score is like probablity, but no cap
      Score
   };
   
   bool should_hang_down(graph_data_type type);
   bool is_probability_plot(graph_data_type type);

   class chain_validation_information_t {
   public:
      std::string chain_id;
      /// what is this used for?
      std::string name;
      /// if this corresponds to `validation_information_t::type`, then I guess it should be removed
      /// (in order to have a single source of truth)
      std::string type;
      std::vector<residue_validation_information_t> rviv;
      explicit chain_validation_information_t(const std::string &chain_id_in) : chain_id(chain_id_in) {}
      void add_residue_validation_information(const residue_validation_information_t &rvi) {
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

   class validation_information_t {
   public:
      std::string name;
      std::vector<chain_validation_information_t> cviv;
      unsigned int get_index_for_chain(const std::string &chain_id);
      /// This is only used with EMSCRIPTEN. The validation graph widget doesn't use it.
      validation_information_min_max_t min_max;

#ifdef EMSCRIPTEN
      std::string type;
      validation_information_t() : min_max(validation_information_min_max_t()), type("UNSET") {}
      validation_information_t(const std::string  &gdt, const validation_information_min_max_t &min_max_in) : min_max(min_max_in), type(gdt) {}
#else
      enum graph_data_type type;
      validation_information_t() : min_max(validation_information_min_max_t()) {type = graph_data_type::UNSET;}
      validation_information_t(graph_data_type gdt, const validation_information_min_max_t &min_max_in) : min_max(min_max_in), type(gdt) {}
#endif

      void add_residue_validation_information(const residue_validation_information_t &rvi, const std::string &chain_id) {
         unsigned int idx = get_index_for_chain(chain_id);
         cviv[idx].add_residue_validation_information(rvi);
      }
      bool empty() const { return cviv.empty(); }
      void set_min_max() {
         unsigned int n = 0;
         double min =  9999999999999;
         double max = -9999999999999;
         for (const auto &chain : cviv) {
            for (const auto &res : chain.rviv) {
               n++;
               if (res.function_value < min) min = res.function_value;
               if (res.function_value > max) max = res.function_value;
            }
         }
         if (n > 0) {
            min_max.min = min;
         } else {
            min_max.max = max;
         }
      }
   };

}

#endif // VALIDATION_INFORMATION_HH
