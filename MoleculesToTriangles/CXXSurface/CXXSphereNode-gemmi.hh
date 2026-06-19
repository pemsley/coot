/*
 * MoleculesToTriangles/CXXSurface/CXXSphereNode-gemmi.hh
 *
 * gemmi-native twin of CXXSphereNode.h. Atom token const gemmi::Atom*;
 * CXXCircle used only via pointer (forward-declared coot::m2t::CXXCircle).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXSphereNode_gemmi_included
#define CXXSphereNode_gemmi_included
#include <vector>
#include "CXXCoord.h"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXCircle;  // pointer only

      class CXXSphereNode {
      private:
         CXXCoord<CXXCoord_ftype> theVertex;
         const CXXCircle *theIntersector;
         int shouldBeDrawn;
         const gemmi::Atom* theAtom;
      public:
         CXXSphereNode();
         CXXSphereNode(const CXXCoord<CXXCoord_ftype>&aCoord);
         const CXXCoord<CXXCoord_ftype>& vertex() const;
         int setDoDraw(const int yesNo);
         int doDraw() const;
         int setVertex(const CXXCoord<CXXCoord_ftype>&);
         void setIntersector(const CXXCircle *aCircle);
         const CXXCircle *getIntersector() const;
         int operator < (const CXXSphereNode &comparator) const {
            return (theVertex[2] < comparator.vertex()[2]);
         }
         void setAtom(const gemmi::Atom* anAtom) { theAtom = anAtom; }
         const gemmi::Atom* getAtom() const { return theAtom; }
      };
   }
}

#endif
