/*
 * MoleculesToTriangles/CXXClasses/SolidColorRule-gemmi.hh
 *
 * gemmi-native twin of SolidColorRule.h. coot::m2t::SolidColorRule. Colour-name
 * resolution reuses the original SolidColorRule's static (mmdb-free) colour table.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef SolidColorRule_gemmi_hh
#define SolidColorRule_gemmi_hh

#include "ColorRule-gemmi.hh"
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include <string>
#include <memory>

namespace coot {
   namespace m2t {

      class SolidColorRule : public ColorRule {
      private:
         FCXXCoord color;
      public:
         SolidColorRule() { type = Solid; }
         FCXXCoord getColor() const { return color; }
         void setColor(const FCXXCoord _color) { color = _color; }
         virtual FCXXCoord colorForAtom(const gemmi::CRA &) override { return color; }

         static std::shared_ptr<SolidColorRule>
            colorRuleForSelectionStringAndName(std::string selectionString, std::string _name);
         static std::shared_ptr<SolidColorRule>
            colorRuleForSelectionAndName(std::shared_ptr<compound_selection_t> selection, std::string _name);
         static std::shared_ptr<SolidColorRule>
            colorRuleForSelectionAndColor(std::shared_ptr<compound_selection_t> selection, FCXXCoord colorCoord);
         static std::shared_ptr<SolidColorRule>
            colorRuleForSelectionStringAndColor(std::string selectionString, FCXXCoord color);
      };
   }
}

#endif // SolidColorRule_gemmi_hh
