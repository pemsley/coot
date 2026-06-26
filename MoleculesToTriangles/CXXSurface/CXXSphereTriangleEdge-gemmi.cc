/*
 * MoleculesToTriangles/CXXSurface/CXXSphereTriangleEdge-gemmi.cc
 *
 * gemmi-native twin of CXXSphereTriangleEdge.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */
#define _USE_MATH_DEFINES
#include <math.h>
#include "CXXSphereTriangleEdge-gemmi.hh"
#include "CXXSphereElement-gemmi.hh"
#include "CXXSphereNode-gemmi.hh"
#include "CXXCoord.h"

using namespace coot::m2t;

CXXSphereTriangleEdge::CXXSphereTriangleEdge() {}
CXXSphereTriangleEdge::~CXXSphereTriangleEdge() {}

CXXSphereTriangleEdge::CXXSphereTriangleEdge(const CXXCoord<CXXCoord_ftype>&aNormal, size_t iV1, size_t iV2,
                                             const CXXCoord<CXXCoord_ftype>&anEdgeCentre, const CXXCoord<CXXCoord_ftype>&aSphereCentre,
                                             double anEdgeRadius, CXXSphereElement *aSphereElement) :
   edgeNormal(aNormal), theEdgeCentre(anEdgeCentre), theSphereCentre(aSphereCentre),
   edgeRadius(anEdgeRadius), theSphereElement(aSphereElement)
{
   edgeVertices[0] = iV1;
   edgeVertices[1] = iV2;
   calculateLength();
}

const CXXCoord<CXXCoord_ftype>& CXXSphereTriangleEdge::normal() const { return edgeNormal; }
size_t CXXSphereTriangleEdge::vertex(size_t i) const { return edgeVertices[i]; }
const CXXCoord<CXXCoord_ftype>& CXXSphereTriangleEdge::edgeCentre() const { return theEdgeCentre; }
const CXXCoord<CXXCoord_ftype>& CXXSphereTriangleEdge::sphereCentre() const { return theSphereCentre; }

double CXXSphereTriangleEdge::calculateLength() {
   CXXCoord<CXXCoord_ftype>v1 = theSphereElement->vertex(edgeVertices[0]).vertex() - theEdgeCentre;
   CXXCoord<CXXCoord_ftype>v2 = theSphereElement->vertex(edgeVertices[1]).vertex() - theEdgeCentre;
   v1.normalise();
   v2.normalise();
   edgeLength = edgeNormal.angleBetween(v1, v2);
   return edgeLength;
}

CXXCoord<CXXCoord_ftype> CXXSphereTriangleEdge::midpoint() const {
   CXXCoord<CXXCoord_ftype>result;
   CXXCoord<CXXCoord_ftype>v1 = theSphereElement->vertex(edgeVertices[0]).vertex() - theEdgeCentre;
   CXXCoord<CXXCoord_ftype>v2 = theSphereElement->vertex(edgeVertices[1]).vertex() - theEdgeCentre;
   if (fabs(edgeLength-0.)> 1e-8 && fabs(edgeLength-M_PI)>1e-8) {
      CXXCoord<CXXCoord_ftype>v3 = v1 + v2;
      v3.normalise();
      v3.scale(edgeRadius);
      if (edgeLength>M_PI) v3.scale(-1.);
      result = theEdgeCentre + v3;
   }
   else if (fabs(edgeLength-0.)<= 1e-8) {
      result = theSphereElement->vertex(edgeVertices[0]).vertex();
   }
   else {
      result = edgeNormal^v1;
      result = result + theEdgeCentre;
   }
   return CXXCoord<CXXCoord_ftype>(result);
}

double CXXSphereTriangleEdge::radius() const { return edgeRadius; }
CXXSphereElement *CXXSphereTriangleEdge::sphereElement() const { return theSphereElement; }
int CXXSphereTriangleEdge::setSphereElement(CXXSphereElement *se) { theSphereElement = se; return 0; }
int CXXSphereTriangleEdge::setSphereCentre(const CXXCoord<CXXCoord_ftype>&crd) { theSphereCentre = crd; return 0; }
int CXXSphereTriangleEdge::setEdgeCentre(const CXXCoord<CXXCoord_ftype>&crd) { theEdgeCentre = crd; return 0; }
