/*
 * validation-graphs/validation-information.hh
 *
 * Copyright 2023 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef VALIDATION_INFORMATION_HH
#define VALIDATION_INFORMATION_HH

#include <vector>
#include "residue-validation-information.hh"

namespace coot {
   /// Validation graph type
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
   /// Probablity plots use logarithmic scale
   bool is_probability_plot(graph_data_type type);

   /// Represents graph data for a single chain in validation graphs
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

   /// Represents a set minimum-maximum boundaries for residue values used for drawing validation graphs
   /// 
   /// This is only used for EMSCRIPTEN. The validation graph widget doesn't use it
   class validation_information_min_max_t {
      public:
         bool is_set;
         double min;
         double max;
         validation_information_min_max_t() : is_set(false), min(0), max(0) {}
         validation_information_min_max_t(const double &min_in, const double &max_in) : is_set(true), min(min_in), max(max_in) {}
   };

   /// Represents graph data a validation graph
   class validation_information_t {
   public:
      std::string name;
      enum graph_data_type type;
      std::vector<chain_validation_information_t> cviv;
      unsigned int get_index_for_chain(const std::string &chain_id);
      /// This is only used with EMSCRIPTEN. The validation graph widget doesn't use it.
      validation_information_min_max_t min_max;

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

#ifdef EMSCRIPTEN
      
      validation_information_t() : min_max(validation_information_min_max_t()), type("UNSET") {}
      validation_information_t(graph_data_type gdt, const validation_information_min_max_t &min_max_in) : min_max(min_max_in), type(gdt) {}

      
#else // No EMSCRIPTEN
      
      validation_information_t() : type(graph_data_type::UNSET) {}
      validation_information_t(graph_data_type gdt) : type(gdt) {}
#endif // EMSCRIPTEN

      void add_residue_validation_information(const residue_validation_information_t &rvi, const std::string &chain_id) {
         unsigned int idx = get_index_for_chain(chain_id);
         cviv[idx].add_residue_validation_information(rvi);
      }
      bool empty() const { return cviv.empty(); }
   };

}

#endif // VALIDATION_INFORMATION_HH
