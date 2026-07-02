/*
 * MoleculesToTriangles/CXXClasses/FlatFanPrimitive-gemmi.hh
 *
 * gemmi-native twin of FlatFanPrimitive.h. Stores ring vertex POSITIONS (the
 * original stored mmdb::Atom* and used only their coords) + a colour. Type is
 * coot::m2t::FlatFanPrimitive.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef FlatFanPrimitive_gemmi_hh
#define FlatFanPrimitive_gemmi_hh

#include <vector>
#include "VertexColorNormalPrimitive.h"

namespace coot {
   namespace m2t {

      class FlatFanPrimitive : public VertexColorNormalPrimitive {
      private:
         std::vector<FCXXCoord> positions;
         FCXXCoord color;
      public:
         FlatFanPrimitive(const std::vector<FCXXCoord> &positions_in, FCXXCoord &color_in)
            : positions(positions_in), color(color_in) {}
         virtual ~FlatFanPrimitive() { emptyArrays(); }
         virtual void generateArrays();
      };
   }
}

#endif // FlatFanPrimitive_gemmi_hh
