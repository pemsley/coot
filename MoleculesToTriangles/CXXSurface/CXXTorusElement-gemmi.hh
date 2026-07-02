/*
 * MoleculesToTriangles/CXXSurface/CXXTorusElement-gemmi.hh
 *
 * gemmi-native twin of CXXTorusElement.h (holds coot::m2t::CXXTorusNode/CXXTriangle).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */
#ifndef CXXTorusElement_gemmi_included
#define CXXTorusElement_gemmi_included
#include "CXXTorusNode-gemmi.hh"
#include "CXXTriangle-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"
#include "CXXCoord.h"
#include <vector>
#include <list>

namespace coot {
   namespace m2t {

      class CXXSurface;
      class CXXCircle;

      class CXXTorusElement {
      private:
         static CXXCircle nullCircle;
         const CXXCircle &theCircle;
         std::vector<CXXTorusNode> nodes;
         std::list<CXXTriangle> flatTriangles;
         double omega1;
         double omega2;
         int nOmegaSteps;
         double deltaOmega;
         int nThetaSteps;
         double deltaTheta;
         double absoluteStartOmega;
         double theta1;
         double theta2;
         std::list<CXXTriangle *> edgeTriangles;
         CXXCoord<CXXCoord_ftype> v1unit;
         CXXCoord<CXXCoord_ftype> v2unit;
         CXXCoord<CXXCoord_ftype> n1unit;
         CXXCoord<CXXCoord_ftype> torusAxis;
         CXXCoord<CXXCoord_ftype> torusCentre;
         double rProbe, rTraj;
         void init();
         int debug;
      public:
         CXXTorusElement();
         ~CXXTorusElement();
         CXXTorusElement(const CXXCircle &aCircle, int iEdge, double delta, double probeRadius);
         void deleteLastTriangle(void);
         size_t addNode(CXXTorusNode &aNode);
         const size_t numberOfTorusNodes(void);
         const CXXCoord<CXXCoord_ftype> coordFromThetaOmega(double theta, double omega) const;
         const CXXCoord<CXXCoord_ftype> normalToProbeAtTheta(CXXCoord<CXXCoord_ftype>&p, double theta) const;
         const CXXCoord<CXXCoord_ftype> probeAtOmega(double omega) const;
         const CXXTorusNode &getNode(const int i) const;
         int upload(CXXSurface *aSurface);
         void addEdgeVertex(CXXCircleNode &aNode);
         int getNOmegaSteps() const { return nOmegaSteps; }
         double getDeltaOmega() const { return deltaOmega; }
         double getAbsoluteStartOmega() const { return absoluteStartOmega; }
         const CXXCircle &getCircle() const { return theCircle; }
         double getTheta2() const { return theta2; }
         const CXXTorusNode &node(size_t i) const { return nodes[i]; }
         size_t nTorusNodes() const { return nodes.size(); }
         size_t nFlatTriangles() const { return flatTriangles.size(); }
         std::list<CXXTriangle>::const_iterator firstTriangle() const { return flatTriangles.begin(); }
         std::list<CXXTriangle>::const_iterator endOfTriangles() const { return flatTriangles.end(); }
      };
   }
}

#endif
