/*
 * MoleculesToTriangles/CXXSurface/CXXSphereTriangleEdge-gemmi.hh
 *
 * gemmi-native twin of CXXSphereTriangleEdge.h. No atom token; stores pointers to
 * the coot::m2t CXXSphereElement / CXXCircle.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXSphereTriangleEdge_gemmi_included
#define CXXSphereTriangleEdge_gemmi_included
#include <vector>
#include "CXXCoord.h"

namespace coot {
   namespace m2t {

      class CXXSphereNode;
      class CXXSphereElement;
      class CXXCircle;

      class CXXSphereTriangleEdge {
      private:
         double edgeLength;
         CXXCoord<CXXCoord_ftype> edgeNormal;
         size_t edgeVertices[2];
         CXXCoord<CXXCoord_ftype> theEdgeCentre;
         CXXCoord<CXXCoord_ftype> theSphereCentre;
         double calculateLength();
         double edgeRadius;
         CXXSphereElement *theSphereElement;
         CXXCircle *theCircle;
      public:
         CXXSphereTriangleEdge();
         ~CXXSphereTriangleEdge();
         CXXSphereTriangleEdge(const CXXCoord<CXXCoord_ftype>&aNormal, size_t iV1, size_t iV2,
                               const CXXCoord<CXXCoord_ftype>&anEdgeCentre, const CXXCoord<CXXCoord_ftype>&aSphereCentre,
                               double anEdgeRadius, CXXSphereElement *aSphereElement);
         double length() const { return edgeLength; }
         const CXXCoord<CXXCoord_ftype>& normal() const;
         size_t vertex(size_t) const;
         const CXXCoord<CXXCoord_ftype>& edgeCentre() const;
         const CXXCoord<CXXCoord_ftype>& sphereCentre() const;
         CXXCoord<CXXCoord_ftype> midpoint() const;
         double radius() const;
         CXXSphereElement *sphereElement() const;
         int setSphereElement(CXXSphereElement *se);
         int setSphereCentre(const CXXCoord<CXXCoord_ftype>&crd);
         int setEdgeCentre(const CXXCoord<CXXCoord_ftype>&crd);
         void setLength(double aLength) { edgeLength = aLength; }
         void setCircle(CXXCircle *aCircle) { theCircle = aCircle; }
         CXXCircle *getCircle() const { return theCircle; }
         void setVertex(size_t index, size_t value) { edgeVertices[index] = value; }
      };
   }
}

#endif
