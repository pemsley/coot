/*
 * MoleculesToTriangles/CXXClasses/ColorRule-gemmi.hh
 *
 * gemmi-native twin of ColorRule.h (mmdb->gemmi migration). colorForAtom takes a
 * gemmi::CRA (chain/residue/atom). Type is coot::m2t::ColorRule. Lives alongside
 * the original.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef ColorRule_gemmi_hh
#define ColorRule_gemmi_hh

#include <string>
#include <memory>
#include <gemmi/model.hpp>   // gemmi::CRA
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include "MoleculesToTriangles/CXXClasses/gemmi-selection.hh"

namespace coot {
   namespace m2t {

      class ColorRule {
      protected:
         std::shared_ptr<compound_selection_t> compoundSelection;
         float rank;
         int type;
      public:
         ColorRule() : compoundSelection(nullptr), rank(-1) {}
         virtual ~ColorRule() = default;
         enum ColorRuleType { Solid, AtomPropertyRamp };
         float getRank() const { return rank; }
         void setRank(const float &_rank) { rank = _rank; }
         std::shared_ptr<compound_selection_t> getCompoundSelection() { return compoundSelection; }
         void setCompoundSelection(std::shared_ptr<compound_selection_t> _cs) { compoundSelection = _cs; }
         int getType() { return type; }
         virtual FCXXCoord colorForAtom(const gemmi::CRA &cra) = 0;
         static bool compareRank(const std::shared_ptr<ColorRule> rule1,
                                 const std::shared_ptr<ColorRule> rule2) {
            return rule1->getRank() < rule2->getRank();
         }
      };
   }
}

#endif // ColorRule_gemmi_hh
