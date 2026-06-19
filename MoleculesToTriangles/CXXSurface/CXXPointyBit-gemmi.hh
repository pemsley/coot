/*
 * MoleculesToTriangles/CXXSurface/CXXPointyBit-gemmi.hh
 *
 * gemmi-native twin of CXXPointyBit.h. Atom token is const gemmi::Atom*.
 * Type is coot::m2t::PointyBit.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXPointyBit_gemmi_included
#define CXXPointyBit_gemmi_included

#include "CXXCoord.h"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class PointyBit {
      public:
         const gemmi::Atom* atomI;
         const gemmi::Atom* atomJ;
         CXXCoord<CXXCoord_ftype> coord;
         int isNull;
         PointyBit() : atomI(0), atomJ(0), isNull(1) {}
         PointyBit(const gemmi::Atom* _atomI, const gemmi::Atom* _atomJ,
                   const CXXCoord<CXXCoord_ftype> &_coord, int _isNull) :
            atomI(_atomI), atomJ(_atomJ), coord(_coord), isNull(_isNull) {}
      };
   }
}

#endif
