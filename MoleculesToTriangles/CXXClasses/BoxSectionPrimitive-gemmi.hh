/*
 * MoleculesToTriangles/CXXClasses/BoxSectionPrimitive-gemmi.hh
 *
 * gemmi-native twin of BoxSectionPrimitive.h. Subclass of coot::m2t::Cylinders
 * Primitive (shares its atom-free CylinderPoint list via the friend declaration);
 * renders a rectangular box cross-section (beta-strand ribbon). Type is
 * coot::m2t::BoxSectionPrimitive.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef BoxSectionPrimitive_gemmi_hh
#define BoxSectionPrimitive_gemmi_hh

#include <memory>
#include <utility>
#include "CylindersPrimitive-gemmi.hh"

class Renderer;

namespace coot {
   namespace m2t {

      class BoxSectionPrimitive : public CylindersPrimitive {
      public:
         BoxSectionPrimitive();
         virtual ~BoxSectionPrimitive() { emptyArrays(); }
         virtual void generateArrays() override;
         virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer) override;
         std::pair<bool, FCXXCoord> boxEdgeColor;
         void setBoxEdgeColor(const FCXXCoord &col_in) { boxEdgeColor.first = true; boxEdgeColor.second = col_in; }
         std::pair<bool, FCXXCoord> getBoxEdgeColor() { return boxEdgeColor; }
      };
   }
}

#endif // BoxSectionPrimitive_gemmi_hh
