/*
 * MoleculesToTriangles/CXXSurface/CXXSphereElement-gemmi.hh
 *
 * gemmi-native twin of CXXSphereElement.h. Atom token const gemmi::Atom*.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXSphereElement_gemmi_included
#define CXXSphereElement_gemmi_included
#include <vector>
#include <list>
#include <map>
#include "CXXCoord.h"
#include "CXXSphereTriangleEdge-gemmi.hh"
#include "CXXSphereNode-gemmi.hh"
#include "CXXSphereTriangle-gemmi.hh"
#include "CXXSphereFlatTriangle-gemmi.hh"
#include "CXXCircle-gemmi.hh"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXTriangle;
      class CXXSurface;
      class CXXCircleNode;
      class CXXTorusElement;
      class CXXBall;

      class TriangleEdgePair {
      public:
         int triangle;
         int edge;
         TriangleEdgePair(int i, int j) : triangle(i), edge(j) {}
      };

      class CXXSphereElement {
      private:
         const gemmi::Atom* theAtom;
         CXXCoord<CXXCoord_ftype> theCentre;
         std::vector<CXXSphereNode> theVertices;
         std::vector<CXXSphereTriangle> theTriangles;
         std::vector<CXXSphereTriangleEdge> theEdges;
         std::list<CXXSphereFlatTriangle> flatTriangles;
         std::vector<std::vector<CXXCircle> > theCircles;
         std::map<CXXSphereFlatTriangle *, int> edgeTriangles;
         double theRadius;
         double deltaRadians;
         void init();
      public:
         CXXSphereElement();
         CXXSphereElement(const CXXCoord<CXXCoord_ftype>&position, double radius, double del);
         CXXSphereElement(const gemmi::Atom* anAtom, double delta);

         void initWith(const CXXCircleNode &aNode, double del, double radius_in, bool *includeAtoms, int UseOrGenerate);
         CXXSphereElement(const CXXSphereElement &oldOne);
         size_t addVertex(const CXXSphereNode &vert);
         size_t addTriangle(const CXXSphereTriangle &);
         size_t addEdge(const CXXSphereTriangleEdge &);
         int addFlatTriangle(const CXXSphereFlatTriangle &);
         void flattenLastTriangle(void);

         const std::vector<CXXSphereNode> &getVertices() const { return theVertices; }
         const std::vector<CXXSphereTriangle> &getTriangles() const { return theTriangles; }
         const std::vector<CXXSphereTriangleEdge> &getEdges() const { return theEdges; }
         const std::list<CXXSphereFlatTriangle> &getFlatTriangles() const { return flatTriangles; }
         const std::vector<std::vector<CXXCircle> > &getCircles() const { return theCircles; }
         const CXXCoord<CXXCoord_ftype>& centre() const;
         const size_t nVertices() const;
         const CXXSphereNode &vertex(const size_t iVertex) const;
         void moveVertex(const int iVertex, const CXXCoord<CXXCoord_ftype>&position);
         const CXXSphereTriangle &triangle(const size_t iTriangle) const;
         const size_t nFlatTriangles() const;
         const double radius() const;
         const size_t nEdges() const;
         const CXXSphereTriangleEdge &edge(const size_t iEdge) const;
         const gemmi::Atom* getAtom() const;
         const double delta() const;

         int calculate();
         int trimBy(const CXXCircle &aCircle, int carefully);
         int upLoad(CXXSurface *aSurface, double probeRadius, const int sense) const;
         int translateBy(const CXXCoord<CXXCoord_ftype>&crd);
         int scaleBy(const double factor);
         int setAtom(const gemmi::Atom* anAtom);

         int addTriangularPatch(const CXXCoord<CXXCoord_ftype>&u1,
                                const CXXCoord<CXXCoord_ftype>&u2,
                                const CXXCoord<CXXCoord_ftype>&u3, const gemmi::Atom* ,
                                std::vector<CXXCircle> &circles,
                                int UseOrGenerate);

         static const int Reentrant = 1;
         static const int Contact = 2;
         static const int VDW = 3;
         static const int Accessible = 4;
         static const int GenerateCircles = 0;
         static const int UseCircles = 1;

         CXXCoord<CXXCoord_ftype> voronoiPoint(const CXXCoord<CXXCoord_ftype>&a, const CXXCoord<CXXCoord_ftype>&b, const CXXCoord<CXXCoord_ftype>&c) const;

         void identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXCircleNode*, std::vector<CXXCoord<CXXCoord_ftype> > >&raggedEdges);
         void identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXBall*, std::vector<CXXCoord<CXXCoord_ftype> > >&raggedEdges);

         size_t flagCutTriangles(const CXXCircle &aCircle);
         void addTorusVertices(const CXXTorusElement &aTorus);
         void addCircleVertices(const CXXCircle &theCircle, int iEdge, double delta);
         void clearCutFlags() { edgeTriangles.clear(); }

         int nDrawnTriangles;
         int getNDrawnTriangles() const { return nDrawnTriangles; }
         int countDrawnTriangles() const;
         int addVertex(const CXXCircleNode &aCircle);

         void initWith(const CXXCoord<CXXCoord_ftype>&aCentre, const gemmi::Atom* atomI, const gemmi::Atom* atomJ, const gemmi::Atom* atomK,
                       double delta, double radius_in, const bool *includeAtoms);
      };
   }
}

#endif
