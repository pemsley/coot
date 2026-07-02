/*
 * MoleculesToTriangles/CXXSurface/CXXNewHood-gemmi.hh
 *
 * gemmi-native twin of CXXNewHood.h. Atom token const gemmi::Atom*; selection
 * membership is a std::set<const gemmi::Atom*> (replaces the mmdb selHnd).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */
#ifndef CXXNewHood_gemmi_included
#define CXXNewHood_gemmi_included
#include <vector>
#include <list>
#include <map>
#include <set>
#include "CXXCoord.h"
#include "CXXCircle-gemmi.hh"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXCircleNode;
      class CXXSurface;
      class CXXSphereElement;
      class CXXBall;

      class CXXNewHood {
      private:
         const gemmi::Atom* theAtomI;
         const CXXBall *theBall;
         double theRadius;
         double theProbeRadius;
         CXXCoord<CXXCoord_ftype> theCentre;
         std::list<CXXCircle> theCircles;
         void init();
      public:
         CXXNewHood();
         CXXNewHood(const gemmi::Atom* centralAtom, double radiusOfAtom1, double probeRadius);
         CXXNewHood(const CXXCircleNode &aNode, double probeRadius);

         int addAtom(const gemmi::Atom* candidate, double radiusOfAtom2);
         int findSegments();

         void initWith(const gemmi::Atom* atomI, double radiusOfAtom1, double probeRadius) {
            theAtomI = atomI;
            theRadius = radiusOfAtom1 + probeRadius;
            theProbeRadius = probeRadius;
            theCentre = CXXCoord<CXXCoord_ftype>(atomI->pos.x, atomI->pos.y, atomI->pos.z);
         }
         void initWith(const CXXCircleNode &aNode, double probeRadius);
         void initWith(const CXXBall *aBall);

         const gemmi::Atom* getAtomI() const;
         double getRadius() const;
         double getProbeRadius() const;
         int getNCircles() const;
         const CXXCoord<CXXCoord_ftype>& getCentre() const;

         int addNode(const CXXCircleNode &aNode);
         int addBall(const CXXBall &aBall);
         size_t nCircles() const;
         int nNodes() const;
         const CXXCircleNode &getNode(const int iNode) const;
         static bool containsDrawable(const CXXNewHood &aHood);

         // selHnd -> selection set
         void identifyUniqueNodes(std::vector<CXXCircleNode> &circleNodes,
                                  const std::set<const gemmi::Atom*> &selSet) const;

         void triangulateAsRegularHoodInto(CXXSurface &aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const;
         const std::list<CXXCircle> &getCircles() const { return theCircles; }
         std::list<CXXCircle> &getCircles() { return theCircles; }
         void triangulateAsBallHoodInto(CXXSurface &aSurface, double delta,
                                        std::map<const CXXBall*, std::vector<CXXCoord<CXXCoord_ftype> > > &raggedEdges,
                                        bool useEdges, int insideOrOutside, const CXXSphereElement &unitCellAtOriginForDelta);
      };
   }
}

#endif
