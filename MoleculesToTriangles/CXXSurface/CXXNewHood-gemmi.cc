/*
 * MoleculesToTriangles/CXXSurface/CXXNewHood-gemmi.cc
 *
 * gemmi-native twin of CXXNewHood.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <iostream>
#include <sstream>
#include "CXXSurfaceVertex.h"
#include "CXXSurface-gemmi.hh"
#include "CXXNewHood-gemmi.hh"
#include "CXXCircle-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"
#include "CXXTorusElement-gemmi.hh"
#include "CXXBall-gemmi.hh"
#include "CXXSphereElement-gemmi.hh"

using namespace std;
using namespace coot::m2t;

void CXXNewHood::init() {
   theAtomI = 0;
   theRadius = 0;
   theProbeRadius = 0;
   theCentre = CXXCoord<CXXCoord_ftype>(0.,0.,0.);
   theCircles.resize(0);
}

CXXNewHood::CXXNewHood() { init(); }

CXXNewHood::CXXNewHood(const gemmi::Atom* centralAtom, double radiusOfAtom1, double probeRadius) :
   theAtomI(centralAtom), theRadius(radiusOfAtom1+probeRadius), theProbeRadius(probeRadius) {
   theCentre = CXXCoord<CXXCoord_ftype>(theAtomI->pos.x, theAtomI->pos.y, theAtomI->pos.z);
}

void CXXNewHood::initWith(const CXXCircleNode &aNode, double probeRadius) {
   theAtomI = 0;
   theRadius = probeRadius;
   theProbeRadius = 0;
   theCentre = aNode.getCoord();
}

void CXXNewHood::initWith(const CXXBall *aBall) {
   theBall = aBall;
   theRadius = aBall->getRadius();
   theProbeRadius = 0.;
   theCentre = aBall->getCoord();
}

int CXXNewHood::addAtom(const gemmi::Atom* anAtomJ, double radiusOfAtom2) {
   if (anAtomJ->serial == theAtomI->serial) {
      return 0;
   }
   double radiusOfAtomJ = radiusOfAtom2 + theProbeRadius;
   CXXCoord<CXXCoord_ftype>centreOfAtomJ = CXXCoord<CXXCoord_ftype>(anAtomJ->pos.x, anAtomJ->pos.y, anAtomJ->pos.z);
   CXXCoord<CXXCoord_ftype>centreOfAtomI = theCentre;
   int result = 0;
   if (!centreOfAtomI.isNearly(centreOfAtomJ, 0.0001)) {
      CXXCoord<CXXCoord_ftype>separ = centreOfAtomJ - centreOfAtomI;
      double sumrad = radiusOfAtomJ + theRadius;
      if (separ.get3DLengthSq()<(sumrad*sumrad)) {
         CXXCircle aCircle(this, anAtomJ, radiusOfAtom2, theProbeRadius);
         theCircles.push_back(aCircle);
         result = 1;
      }
   }
   return result;
}

int CXXNewHood::addBall(const CXXBall &aBall) {
   if (!theCentre.isNearly(aBall.getCoord(), 0.0001)) {
      CXXCoord<CXXCoord_ftype>separ = aBall.getCoord() - theCentre;
      double sumrad = theRadius + aBall.getRadius();
      if (separ.get3DLengthSq()<(sumrad*sumrad)) {
         theCircles.push_back(CXXCircle(this, aBall));
         return 1;
      }
   }
   return 0;
}

int CXXNewHood::findSegments() {
   vector<CXXCoord<CXXCoord_ftype> > theIntersections(2);
   CXXCircleNode circleNode;
   std::list<CXXCircle>::iterator circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::iterator circle1Iter = theCircles.begin(); circle1Iter != circlesEnd; ++circle1Iter) {
      CXXCircle &circle1(*circle1Iter);
      if (!circle1.getEaten()) {
         for (std::list<CXXCircle>::iterator circle2Iter = circle1Iter; circle2Iter != circlesEnd; ++circle2Iter) {
            CXXCircle &circle2(*circle2Iter);
            if (&circle1 != &circle2) {
               int eatenBy = circle1.meetsCircle(circle2, theIntersections);
               if (eatenBy == 2) circle1.setEaten(1);
               else if (eatenBy == 1) circle2.setEaten(1);
               else if (eatenBy==0) {
                  for (int i=0; i<2; i++) {
                     bool addNodes = true;
                     for (std::list<CXXCircle>::iterator circle3Iter = theCircles.begin();
                          circle3Iter != circlesEnd && addNodes; ++circle3Iter) {
                        if (circle3Iter != circle1Iter && circle3Iter != circle2Iter) {
                           CXXCircle &circle3(*circle3Iter);
                           if (!circle3.getEaten()) {
                              if (!circle3.accIsBehind(theIntersections[i])) addNodes = false;
                           }
                        }
                     }
                     if (addNodes) {
                        circleNode.setParent(&circle1);
                        circleNode.setOtherCircle(&circle2);
                        circleNode.setCoord(theIntersections[i]);
                        circleNode.setFlag((i==0?1:2));
                        circle1.addNode(circleNode);
                        circleNode.setParent(&circle2);
                        circleNode.setOtherCircle(&circle1);
                        circleNode.setCoord(theIntersections[i]);
                        circleNode.setFlag((i==0?2:1));
                        circle2.addNode(circleNode);
                     } else {
                        circle1.setContainsEatenNodes(1);
                        circle2.setContainsEatenNodes(1);
                     }
                  }
               }
            }
         }
      }
   }
   circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::iterator circleIter = theCircles.begin(); circleIter != circlesEnd; ++circleIter) {
      CXXCircle &circle1(*circleIter);
      bool completeOrbit = (circle1.getNNodes() == 0 && circle1.getContainsEatenNodes()==0);
      if (!circle1.getEaten()) {
         if (completeOrbit) circle1.setArbitraryReference();
         circle1.sortNodes();
         circle1.newIdentifyArcs();
      }
   }
   return 0;
}

const CXXCoord<CXXCoord_ftype>& CXXNewHood::getCentre() const { return theCentre; }
double CXXNewHood::getRadius() const { return theRadius; }
size_t CXXNewHood::nCircles() const { return theCircles.size(); }
const gemmi::Atom* CXXNewHood::getAtomI() const { return theAtomI; }
double CXXNewHood::getProbeRadius() const { return theProbeRadius; }

bool CXXNewHood::containsDrawable(const CXXNewHood &aHood) {
   const std::list<CXXCircle>& theCircles(aHood.getCircles());
   bool result = (theCircles.size() == 0);
   std::list<CXXCircle>::const_iterator circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::const_iterator circleIter = theCircles.begin();
        circleIter != circlesEnd && !result; ++circleIter) {
      const CXXCircle &theCircle(*circleIter);
      if (theCircle.nSegments()!=0) result = true;
   }
   return result;
}

void CXXNewHood::triangulateAsRegularHoodInto(CXXSurface &aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const {
   const CXXNewHood &theNewHood(*this);
   CXXSphereElement vdwSphere(*unitSphereAtOrigin);
   vdwSphere.scaleBy(theNewHood.getRadius());
   vdwSphere.translateBy(theNewHood.getCentre());
   vdwSphere.setAtom(theNewHood.getAtomI());

   std::list<CXXCircle>::const_iterator circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::const_iterator circleIter = theCircles.begin();
        circleIter != circlesEnd && vdwSphere.getNDrawnTriangles(); ++circleIter) {
      const CXXCircle &theCircle(*circleIter);
      if (!theCircle.getEaten()) vdwSphere.trimBy(theCircle, 1);
   }

   map<const CXXCircle *, vector<CXXCircleNode> > rawEdges;
   if (vdwSphere.nFlatTriangles()) {
      for (unsigned iVertex = 0; iVertex < vdwSphere.nVertices(); iVertex++) {
         if (vdwSphere.vertex(iVertex).doDraw()) {
            const CXXCircle *correspondingCircle = vdwSphere.vertex(iVertex).getIntersector();
            if (correspondingCircle != 0) {
               CXXCircleNode extraNode(correspondingCircle, 0, vdwSphere.vertex(iVertex).vertex(), iVertex);
               extraNode.setReference(correspondingCircle->getReferenceUnitRadius());
               rawEdges[correspondingCircle].push_back(extraNode);
            }
         }
      }
   }

   circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::const_iterator circleIter = theCircles.begin(); circleIter != circlesEnd; ++circleIter) {
      const CXXCircle &theCircle(*circleIter);
      if ((!theCircle.getEaten()) && theCircle.getNNodes() != 0) {
         vdwSphere.flagCutTriangles(theCircle);
         for (unsigned iTorus = 0; iTorus < theCircle.nSegments(); iTorus++) {
            CXXTorusElement theTorus(theCircle, iTorus, delta, getProbeRadius());
            for (unsigned iRawEdge=0; iRawEdge < rawEdges[&theCircle].size(); iRawEdge++) {
               theTorus.addEdgeVertex(rawEdges[&theCircle][iRawEdge]);
            }
            aSurface.uploadTorus(theTorus);
            vdwSphere.addTorusVertices(theTorus);
         }
      }
   }
   aSurface.upLoadSphere(vdwSphere, getProbeRadius(), CXXSphereElement::Contact);
}

void CXXNewHood::identifyUniqueNodes(vector<CXXCircleNode> &circleNodes, const std::set<const gemmi::Atom*> &selSet) const {
   std::list<CXXCircle>::const_iterator circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::const_iterator circleIter = theCircles.begin(); circleIter != circlesEnd; ++circleIter) {
      const CXXCircle &theCircle(*circleIter);
      if (!theCircle.getEaten()) {
         for (unsigned iSegment = 0; iSegment<theCircle.nSegments(); iSegment++) {
            const CXXCircleNode *ends[2];
            ends[0] = theCircle.start(iSegment);
            ends[1] = theCircle.stop(iSegment);
            for (int iTerminus = 0; iTerminus<2; iTerminus++) {
               const CXXCircleNode &aNode(*ends[iTerminus]);
               if (aNode.getOtherCircle()) {
                  bool includeNode = false;
                  bool jIn = selSet.count(aNode.getAtomJ()) > 0;
                  bool kIn = selSet.count(aNode.getAtomK()) > 0;
                  if (jIn && kIn) {
                     if (aNode.getAtomI() < aNode.getAtomJ() && aNode.getAtomJ() < aNode.getAtomK()) includeNode = true;
                  }
                  else if (!jIn && !kIn) {
                     if (aNode.getAtomJ() < aNode.getAtomK()) includeNode = true;
                  }
                  else if (jIn && !kIn) {
                     if (aNode.getAtomI() < aNode.getAtomJ()) includeNode = true;
                  }
                  else if (kIn) {
                  }
                  if (includeNode) {
#pragma omp critical (circleNodes)
                     circleNodes.push_back(aNode);
                  }
               }
            }
         }
      }
   }
}

void CXXNewHood::triangulateAsBallHoodInto(CXXSurface &aSurface, double delta,
                                           std::map<const CXXBall*, std::vector<CXXCoord<CXXCoord_ftype> > > &raggedEdges,
                                           bool useEdges, int insideOrOutside,
                                           const CXXSphereElement &unitCellAtOriginForDelta) {
   CXXSphereElement ballSphere;
   theBall->initSphereElement(ballSphere, delta, unitCellAtOriginForDelta);
   int wasTrimmed = 0;
   std::list<CXXCircle>::const_iterator circlesEnd = theCircles.end();
   for (std::list<CXXCircle>::const_iterator circleIter = theCircles.begin();
        circleIter != circlesEnd && ballSphere.getNDrawnTriangles(); ++circleIter) {
      const CXXCircle &theCircle(*circleIter);
      if (theCircle.getBallJ()!=0) {
         if (ballSphere.trimBy(theCircle, 1)==1) wasTrimmed = 1;
      }
   }
   if (!useEdges && !wasTrimmed) {
      aSurface.upLoadSphere(ballSphere, theBall->getRadius(), insideOrOutside);
   }
   else if (!useEdges && wasTrimmed && ballSphere.getNDrawnTriangles()>0) {
      std::list<CXXCircle>::iterator circlesEnd2 = theCircles.end();
      for (std::list<CXXCircle>::iterator circleIter = theCircles.begin(); circleIter != circlesEnd2; ++circleIter) {
         CXXCircle &theCircle(*circleIter);
         if (theCircle.getBallJ()!=0) ballSphere.identifyRaggedEdges(theCircle, raggedEdges);
      }
   }
   else if (useEdges && wasTrimmed && ballSphere.getNDrawnTriangles()>0) {
      std::list<CXXCircle>::iterator circlesEnd2 = theCircles.end();
      for (std::list<CXXCircle>::iterator circleIter = theCircles.begin();
           circleIter != circlesEnd2 && ballSphere.getNDrawnTriangles()>0; ++circleIter) {
         CXXCircle &theCircle(*circleIter);
         const CXXBall *ballJ(theCircle.getBallJ());
         ballSphere.flagCutTriangles(theCircle);
         std::vector<CXXCoord<CXXCoord_ftype> > &ballsEdges(raggedEdges[ballJ]);
         std::vector<CXXCoord<CXXCoord_ftype> >::iterator ballsEdgesEnd(ballsEdges.end());
         for (std::vector<CXXCoord<CXXCoord_ftype> >::iterator ballsEdge = ballsEdges.begin();
              ballsEdge!= ballsEdgesEnd; ballsEdge++) {
            ballSphere.addVertex(CXXCircleNode(&theCircle, 0, *ballsEdge, 0));
         }
      }
#pragma omp critical (mainTriangles)
      aSurface.upLoadSphere(ballSphere, theBall->getRadius(), insideOrOutside);
   }
}
