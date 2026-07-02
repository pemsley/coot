/*
 * MoleculesToTriangles/CXXSurface/CXXSphereFlatTriangle-gemmi.hh
 *
 * gemmi-native twin of CXXSphereFlatTriangle.h (derives from coot::m2t::CXXTriangle).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXSphereFlatTriangle_gemmi_included
#define CXXSphereFlatTriangle_gemmi_included

#include "CXXTriangle-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"

namespace coot {
   namespace m2t {

      class CXXCircle;

      class CXXSphereFlatTriangle : public CXXTriangle {
      private:
         const CXXCircle *edgeCircles[3];
         CXXCircleNode circleNodes[3];
      public:
         CXXSphereFlatTriangle() : CXXTriangle() {
            for (int i=0; i<3; i++) edgeCircles[i] = 0;
         }
         CXXSphereFlatTriangle(size_t i, size_t j, size_t k, size_t l) : CXXTriangle(i, j, k, l) {
            for (size_t ii=0; ii<3; ii++) edgeCircles[ii] = 0;
         }
         CXXSphereFlatTriangle(int i, int j, int k) : CXXTriangle(i, j, k) {
            for (int ii=0; ii<3; ii++) edgeCircles[ii] = 0;
         }
         void setEdgeCircle(const int i, const CXXCircle *edgeCircle) { edgeCircles[i] = edgeCircle; }
         const CXXCircle *getEdgeCircle(const int i) const { return edgeCircles[i]; }
         void setCircleNode(const int i, const CXXCircleNode &aNode) { circleNodes[i] = aNode; }
         const CXXCircleNode &getCircleNode(const int i) const { return circleNodes[i]; }
      };
   }
}

#endif
