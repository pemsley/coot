/*
 * MoleculesToTriangles/CXXSurface/CXXSphereTriangle-gemmi.hh
 *
 * gemmi-native twin of CXXSphereTriangle.h. Atom token const gemmi::Atom*;
 * CXXSphereElement used via pointer (forward-declared).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXSphereTriangle_gemmi_included
#define CXXSphereTriangle_gemmi_included

#include "CXXCoord.h"
#include "CXXSphereTriangleEdge-gemmi.hh"
#include "CXXSphereNode-gemmi.hh"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXSphereElement;  // pointer only

      class CXXSphereTriangle {
      private:
         size_t triangleVertices[3];
         size_t triangleEdges[3];
         double theRadius;
         CXXCoord<CXXCoord_ftype> theCentre;
         CXXSphereElement *theSphereElement;
         const gemmi::Atom* theAtom;
      public:
         CXXSphereTriangle();
         CXXSphereTriangle(CXXSphereElement *se, size_t *vertices, size_t *edges,
                           double aRadius, CXXCoord<CXXCoord_ftype>&aCentre);
         CXXSphereTriangle(CXXSphereElement *se, size_t *vertices, size_t *edges,
                           double aRadius, CXXCoord<CXXCoord_ftype>&aCentre, const gemmi::Atom* anAtom);

         size_t vertex(size_t) const;
         size_t edge(size_t) const;
         double radius() const;
         const CXXCoord<CXXCoord_ftype>& centre() const;
         CXXSphereElement *sphereElement() const;
         int bisect(double);
         const gemmi::Atom* getAtom() const;

         int setAtom(const gemmi::Atom* anAtom);
         int setEdge(size_t i, size_t);
         int setVertex(size_t, size_t);
         int setSphereElement(CXXSphereElement *se);
         int setCentre(const CXXCoord<CXXCoord_ftype>&);
         int setRadius(const double radius);
      };
   }
}

#endif
