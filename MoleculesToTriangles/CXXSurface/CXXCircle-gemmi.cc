/*
 * MoleculesToTriangles/CXXSurface/CXXCircle-gemmi.cc
 *
 * gemmi-native twin of CXXCircle.cpp.
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
#include <iostream>
#include "CXXSurfaceVertex.h"
#include "CXXCircle-gemmi.hh"
#include "CXXNewHood-gemmi.hh"
#include "CXXBall-gemmi.hh"

using namespace std;
using namespace coot::m2t;

CXXCircle::CXXCircle() :
   theAtomJ(0), theBallJ(0), theParent(0),
   centreOfSecondSphere(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   theNormal(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   radiusOfSecondSphere(0.), radiusOfSphere(0.),
   centreOfCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   centreToCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   referenceUnitRadius(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   radiusOfCircle(0.), nIntersectingCircles(0),
   completelyEaten(0), nodeNumber(0), containsEatenNodes(0) {}

CXXCircle::CXXCircle(CXXNewHood *aHood, const gemmi::Atom* atom2, double radiusOfAtom2, double probeRadius) :
   theAtomJ(atom2), theBallJ(0), theParent(aHood),
   centreOfSecondSphere(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   theNormal(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   radiusOfSecondSphere(0.), radiusOfSphere(aHood->getRadius()),
   centreOfCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   centreToCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   referenceUnitRadius(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   radiusOfCircle(0.), nIntersectingCircles(0),
   completelyEaten(0), nodeNumber(0), containsEatenNodes(0)
{
   centreOfSecondSphere = CXXCoord<CXXCoord_ftype>(theAtomJ->pos.x, theAtomJ->pos.y, theAtomJ->pos.z);
   theNormal = centreOfSecondSphere - getCentreOfSphere();
   radiusOfSecondSphere = radiusOfAtom2 + probeRadius;
   performPrecalculations();
}

CXXCircle::CXXCircle(CXXNewHood *aHood, const CXXBall &aBall) :
   theAtomJ(aBall.getAtomI()), theBallJ(&aBall), theParent(aHood),
   centreOfSecondSphere(aBall.getCoord()),
   theNormal(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   radiusOfSecondSphere(aBall.getRadius()), radiusOfSphere(aHood->getRadius()),
   centreOfCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   centreToCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   referenceUnitRadius(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
   radiusOfCircle(0.), nIntersectingCircles(0),
   completelyEaten(0), nodeNumber(0), containsEatenNodes(0)
{
   theNormal = centreOfSecondSphere - getCentreOfSphere();
   performPrecalculations();
}

void CXXCircle::performPrecalculations() {
   theNormal.normalise();
   centreToCircle = centreOfSecondSphere - getCentreOfSphere();
   double Rs, Rc, Rcs;
   Rs = radiusOfSecondSphere;
   Rc = radiusOfSphere;
   Rcs = centreToCircle.get3DLength();
   centreToCircle /= Rcs;
   centreToCircle *= (Rc*Rc + Rcs*Rcs - Rs*Rs)/(2*Rcs);
   radiusOfCircle = sqrt(radiusOfSphere*radiusOfSphere - centreToCircle.get3DLengthSq());
   centreOfCircle = getCentreOfSphere() + centreToCircle;
   CXXCoord<CXXCoord_ftype>unitCentreToCircle = centreToCircle;
   unitCentreToCircle.normalise();
}

int CXXCircle::meetsCircle(const CXXCircle &otherCircle, vector<CXXCoord<CXXCoord_ftype> > &nodeList) const {
   CXXCoord<CXXCoord_ftype>intersectA, intersectB;
   CXXCoord<CXXCoord_ftype>rj(getCentreToCircle());
   CXXCoord<CXXCoord_ftype>rk(otherCircle.getCentreToCircle());
   CXXCoord<CXXCoord_ftype>rjxrk;
   double rjrk, rjSquared, rkSquared;
   double a, b, distSq, radiusOfSphereSq, x;
   rjSquared = rj*rj;
   rkSquared = rk*rk;
   rjrk = rj*rk;
   a = ((rjSquared - rjrk)*rkSquared)/(rjSquared*rkSquared - rjrk*rjrk);
   b = ((rkSquared - rjrk)*rjSquared)/(rjSquared*rkSquared - rjrk*rjrk);
   CXXCoord<CXXCoord_ftype>planeIntersect((rj*a) + (rk*b));
   distSq = planeIntersect*planeIntersect;
   radiusOfSphereSq = radiusOfSphere*radiusOfSphere;
   if (distSq < radiusOfSphereSq) {
      rjxrk = rj ^ rk;
      rjxrk.normalise();
      x = sqrt(radiusOfSphereSq - distSq);
      rjxrk.scale(x);
      CXXCoord<CXXCoord_ftype>unitCentreToCircle(centreToCircle);
      unitCentreToCircle.normalise();
      int forwardOne = (theNormal*unitCentreToCircle > 0.?1:-1);
      CXXCoord<CXXCoord_ftype>otherCircleUnitCentreToCircle(otherCircle.getCentreToCircle());
      otherCircleUnitCentreToCircle.normalise();
      int forwardTwo = (otherCircle.getNormal()*otherCircleUnitCentreToCircle > 0.?1:-1);
      if (forwardOne*forwardTwo > 0) {
         intersectA = planeIntersect+rjxrk;
         intersectB = planeIntersect-rjxrk;
      } else {
         intersectB = planeIntersect+rjxrk;
         intersectA = planeIntersect-rjxrk;
      }
      nodeList[0] = intersectA + getCentreOfSphere();
      nodeList[1] = intersectB + getCentreOfSphere();
      return 0;
   } else {
      if (otherCircle.accIsBehind(getCentreOfCircle()) == 0) return 2;
      if (accIsBehind(otherCircle.getCentreOfCircle()) == 0) return 1;
      return 3;
   }
}

int CXXCircle::isSomewhereInsideSphere(const CXXCoord<CXXCoord_ftype>&centre, const double radius) const {
   CXXCoord<CXXCoord_ftype>xAxis(1.0, 0.0, 0.0);
   CXXCoord<CXXCoord_ftype>yAxis(0.0, 1.0, 0.0);
   CXXCoord<CXXCoord_ftype>radiusVector(0.,0.,0.);
   if (fabs(theNormal*xAxis)<0.9999999) radiusVector = xAxis ^ theNormal;
   else radiusVector = yAxis ^ theNormal;
   radiusVector.normalise();
   radiusVector *= radiusOfCircle;
   CXXCoord<CXXCoord_ftype>pointOnCircle = getCentreOfCircle() + radiusVector;
   CXXCoord<CXXCoord_ftype>vectorToSphere = pointOnCircle - centre;
   double distanceToSphereSq = vectorToSphere.get3DLengthSq();
   if (distanceToSphereSq < radius*radius) return 1;
   else return 0;
}

const CXXCoord<CXXCoord_ftype>& CXXCircle::getCentreOfSecondSphere() const { return centreOfSecondSphere; }
const CXXCoord<CXXCoord_ftype>& CXXCircle::getCentreOfSphere() const { return theParent->getCentre(); }
const CXXCoord<CXXCoord_ftype>& CXXCircle::getCentreToCircle() const { return centreToCircle; }
double CXXCircle::getRadiusOfSphere() const { return radiusOfSphere; }
double CXXCircle::getRadiusOfSecondSphere() const { return radiusOfSecondSphere; }
const gemmi::Atom* CXXCircle::getAtomJ() const { return theAtomJ; }
const CXXCoord<CXXCoord_ftype>& CXXCircle::getNormal() const { return theNormal; }

int CXXCircle::trimNodesBy(const CXXCircle &otherCircle) {
   int deletedCount = 0;
   list<CXXCircleNode>::iterator lastNode = theNodes.end();
   for (list<CXXCircleNode>::iterator nodeIter = theNodes.begin(); nodeIter != lastNode; ++nodeIter) {
      CXXCircleNode &node(*nodeIter);
      if (node.getOtherCircle()!=this) {
         if (!node.isDeleted() && node.getFlag()!=-1) {
            if (!otherCircle.accIsBehind(node.getCoord())) {
               node.setDeleted(1);
               deletedCount ++;
            }
         }
      }
   }
   return deletedCount;
}

bool CXXCircle::abBracketsC(const CXXCircleNode &nodea, const CXXCircleNode &nodeb, const CXXCircleNode &nodec) const {
   double aCrossBdotN = (nodea.getUnitRadius()^nodeb.getUnitRadius())*getNormal();
   double aCrossCdotN = (nodea.getUnitRadius()^nodec.getUnitRadius())*getNormal();
   double bCrossCdotN = (nodeb.getUnitRadius()^nodec.getUnitRadius())*getNormal();
   if (nodea.getUnitRadius().isNearly(nodec.getUnitRadius(), 0.000000000001) ||
       nodeb.getUnitRadius().isNearly(nodec.getUnitRadius(), 0.000000000001)) {
      return false;
   }
   bool liesBetween = false;
   if (aCrossBdotN>0.) liesBetween = aCrossCdotN>0. && bCrossCdotN<0.;
   else liesBetween = aCrossCdotN>0. || bCrossCdotN<0.;
   return liesBetween;
}

bool CXXCircle::smallabBracketsC(const CXXCircleNode &nodea, const CXXCircleNode &nodeb, const CXXCircleNode &nodec) const {
   double aCrossBdotN = (nodea.getUnitRadius()^nodeb.getUnitRadius())*getNormal();
   double aCrossCdotN = (nodea.getUnitRadius()^nodec.getUnitRadius())*getNormal();
   double bCrossCdotN = (nodeb.getUnitRadius()^nodec.getUnitRadius())*getNormal();
   bool liesBetween = false;
   if (aCrossBdotN>0.) liesBetween = aCrossCdotN>0. && bCrossCdotN<0.;
   else liesBetween = aCrossCdotN<0. && bCrossCdotN>0.;
   return liesBetween;
}

void CXXCircle::trimOwnNodes() {
   list<CXXCircle> &otherCircles = theParent->getCircles();
   int nodesToDraw = countDrawnNodes();
   list<CXXCircle>::iterator otherCirclesEnd = otherCircles.end();
   for (list<CXXCircle>::iterator otherCircleIter = otherCircles.begin();
        otherCircleIter != otherCirclesEnd && nodesToDraw > 0; ++otherCircleIter) {
      CXXCircle &otherCircle(*otherCircleIter);
      if (!otherCircle.getEaten() && &otherCircle != this) {
         list<CXXCircleNode>::iterator nodesEnd = theNodes.end();
         for (list<CXXCircleNode>::iterator nodeIter = theNodes.begin();
              nodeIter != nodesEnd && nodesToDraw>0; ++nodeIter) {
            CXXCircleNode &nodec(*nodeIter);
            if (nodec.getOtherCircle() != &otherCircle) {
               if (nodec.isDeleted() == 0) {
                  if (!otherCircle.accIsBehind(nodec.getCoord())) {
                     nodec.setDeleted(1);
                     setContainsEatenNodes(1);
                     nodesToDraw--;
                  }
               }
            }
         }
      }
   }
}

int CXXCircle::sortNodes() {
   if (theNodes.empty()) return 1;
   size_t initialNNodes = theNodes.size();
   if (initialNNodes == 2 && theNodes.front().getFlag() == -1 && theNodes.back().getFlag() == -1) {
      referenceUnitRadius = theNodes.begin()->getUnitRadius();
      theNodes.front().setAngle(0);
      theNodes.back().setFlag(2);
      theNodes.front().setAngle(2.*M_PI);
      theNodes.back().setFlag(1);
      return 0;
   }
   theNodes.remove_if(CXXCircleNode::shouldDelete);
   if (theNodes.size()%2) std::cout << "Seem to have non-even number of nodes (after trimming)\n";
   CXXCircleNode *startNode = 0;
   list<CXXCircleNode>::iterator nodesEnd = theNodes.end();
   for (list<CXXCircleNode>::iterator nodeIter = theNodes.begin(); nodeIter != nodesEnd; ++nodeIter) {
      CXXCircleNode &theNode(*nodeIter);
      if (theNode.getFlag()==2) startNode = &theNode;
   }
   if (startNode == 0) {
      if (theNodes.size()>0) std::cout << theNodes.size() << "nodes but no startpoint\n";
      theNodes.resize(0);
      return 0;
   }
   referenceUnitRadius = startNode->getUnitRadius();
   nodesEnd = theNodes.end();
   for (list<CXXCircleNode>::iterator nodeIter = theNodes.begin(); nodeIter != nodesEnd; ++nodeIter) {
      CXXCircleNode &theNode(*nodeIter);
      if (&theNode == startNode) theNode.setAngle(0.);
      else theNode.setReference(referenceUnitRadius);
   }
   theNodes.sort(CXXCircleNode::angleLessThan);
   return 0;
}

int CXXCircle::newIdentifyArcs() {
   size_t nNodes = theNodes.size();
   if (getEaten() || nNodes==0) {
      theStarts.resize(0);
      theStops.resize(0);
      return 1;
   }
   const gemmi::Atom* atomK(theNodes.begin()->getAtomK());
   if (!atomK && theNodes.size()>1) {
      theStarts.push_back(&(*theNodes.begin()));
      theStops.push_back(&theNodes.back());
      return 0;
   }
   list<CXXCircleNode>::iterator nodesEnd = theNodes.end();
   for (list<CXXCircleNode>::iterator nodeIter = theNodes.begin(); nodeIter != nodesEnd; ++nodeIter) {
      CXXCircleNode &theNode(*nodeIter);
      if (theNode.isDeleted() == 0) {
         if (theNode.getFlag() == 2) theStarts.push_back(&theNode);
         else theStops.push_back(&theNode);
      }
   }
   if (theStarts.size() != theStops.size()) {
      std::cout << "uneven count of starts and stops\n";
      theStarts.resize(0);
   }
   return 0;
}

CXXNewHood *CXXCircle::getParent() const { return theParent; }
size_t CXXCircle::getNNodes() const { return theNodes.size(); }
int CXXCircle::addNode(const CXXCircleNode &aNode) { theNodes.push_back(aNode); return 0; }
int CXXCircle::getEaten() const { return completelyEaten; }
void CXXCircle::setEaten(int flag) { completelyEaten = flag; if (flag) theNodes.clear(); }
size_t CXXCircle::nSegments() const { return theStarts.size(); }
CXXCircleNode* CXXCircle::start(const int i) const { return theStarts[i]; }
CXXCircleNode* CXXCircle::stop(const int i) const { return theStops[i]; }

double CXXCircle::getRadiusOfVdWCircle() const {
   return getRadiusOfCircle() * (getRadiusOfSphere()-theParent->getProbeRadius()) / getRadiusOfSphere();
}

const CXXCoord<CXXCoord_ftype> CXXCircle::getCentreOfVdWCircle() const {
   CXXCoord<CXXCoord_ftype>diff = getCentreToCircle();
   double frac = (getRadiusOfSphere()-theParent->getProbeRadius()) / getRadiusOfSphere();
   diff *= frac;
   return getCentreOfSphere()+diff;
}

void CXXCircle::dumpVdw() const {
   CXXCoord<CXXCoord_ftype>aRadius;
   if (theNormal*CXXCoord<CXXCoord_ftype>(1.,0.,0.) < 0.9999999) aRadius = theNormal^CXXCoord<CXXCoord_ftype>(1.,0.,0.);
   else aRadius = theNormal^CXXCoord<CXXCoord_ftype>(0.,1.,0.);
   aRadius.normalise();
   CXXCoord<CXXCoord_ftype>bRadius = theNormal^aRadius;
   CXXCoord<CXXCoord_ftype>offset = theNormal;
   offset *= 0.001;
   CXXCoord<CXXCoord_ftype>cCent = getCentreOfVdWCircle()+offset;
   aRadius *= getRadiusOfVdWCircle();
   bRadius *= getRadiusOfVdWCircle();
   for (int i=0; i<12; i++) {
      cout << "add triangle ";
      for (int j=0; j<3; j++) cout << cCent[j] << " ";
      for (int j=0; j<3; j++) cout << cCent[j] + cos(30.*M_PI/180.*double(i+1))*aRadius[j] + sin(30.*M_PI/180.*double(i+1))*bRadius[j] << " ";
      for (int j=0; j<3; j++) cout << cCent[j] + cos(30.*M_PI/180.*double(i))*aRadius[j] + sin(30.*M_PI/180.*double(i))*bRadius[j] << " ";
      cout << endl;
   }
}

int CXXCircle::vdwIsBehind(const CXXCoord<CXXCoord_ftype>&crd) const {
   CXXCoord<CXXCoord_ftype>diff = crd - getCentreOfVdWCircle();
   return (diff*getNormal() < 0.);
}

const CXXCoord<CXXCoord_ftype> CXXCircle::vdwPlaneIntersect(const CXXCoord<CXXCoord_ftype>&A, const CXXCoord<CXXCoord_ftype>&B) const {
   CXXCoord<CXXCoord_ftype>AO = getCentreOfVdWCircle() - A;
   CXXCoord<CXXCoord_ftype>AB = B - A;
   double frac = (AO*getNormal()) / (AB*getNormal());
   AB *= frac;
   return A + AB;
}

const CXXCoord<CXXCoord_ftype> CXXCircle::accPlaneIntersect(const CXXCoord<CXXCoord_ftype>&A, const CXXCoord<CXXCoord_ftype>&B) const {
   CXXCoord<CXXCoord_ftype>AO(getCentreOfCircle() - A);
   CXXCoord<CXXCoord_ftype>AB(B - A);
   double frac = (AO*getNormal()) / (AB*getNormal());
   AB *= frac;
   CXXCoord<CXXCoord_ftype>result(A + AB);
   CXXCoord<CXXCoord_ftype>OC = result - getCentreOfCircle();
   double radialLength = OC.get3DLength();
   OC *= getRadiusOfCircle()/radialLength;
   return getCentreOfCircle() + OC;
}

const CXXCoord<CXXCoord_ftype>& CXXCircle::getReferenceUnitRadius() const { return referenceUnitRadius; }

void CXXCircle::setArbitraryReference() {
   CXXCoord<CXXCoord_ftype>v1unit;
   CXXCoord<CXXCoord_ftype>xAxis(1.,0.,0.,0.);
   CXXCoord<CXXCoord_ftype>zAxis(0.,0.,1.,0.);
   if (fabs(getNormal()*zAxis) < 0.9999) { v1unit = (getNormal() ^ zAxis); v1unit.normalise(); }
   else { v1unit = (getNormal() ^ xAxis); v1unit.normalise(); }
   referenceUnitRadius = v1unit;
   v1unit *= getRadiusOfCircle();
   CXXCoord<CXXCoord_ftype>arbitraryEdgePoint(getCentreOfCircle() + v1unit);
   addNode(CXXCircleNode(this, 0, arbitraryEdgePoint, -1));
   addNode(CXXCircleNode(this, 0, arbitraryEdgePoint, -1));
}

int CXXCircle::countDrawnNodes() const {
   int answer = 0;
   list<CXXCircleNode>::const_iterator nodesEnd = theNodes.end();
   for (list<CXXCircleNode>::const_iterator nodeIter = theNodes.begin(); nodeIter != nodesEnd; ++nodeIter) {
      const CXXCircleNode &theNode(*nodeIter);
      if (theNode.isDeleted()==0) answer++;
   }
   return answer;
}
