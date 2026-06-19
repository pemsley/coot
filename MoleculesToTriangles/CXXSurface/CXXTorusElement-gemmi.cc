/*
 * MoleculesToTriangles/CXXSurface/CXXTorusElement-gemmi.cc
 *
 * gemmi-native twin of CXXTorusElement.cpp.
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
#include "CXXTorusElement-gemmi.hh"
#include "CXXTorusNode-gemmi.hh"
#include "CXXCoord.h"
#include "CXXSurface-gemmi.hh"
#include "CXXTriangle-gemmi.hh"
#include "CXXCircle-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"
#include "CXXNewHood-gemmi.hh"

using namespace std;
using namespace coot::m2t;

CXXCircle CXXTorusElement::nullCircle = CXXCircle();

CXXTorusElement::CXXTorusElement() : theCircle(nullCircle), debug(0) {
   init();
}

CXXTorusElement::~CXXTorusElement() {}

void CXXTorusElement::init() {
   flatTriangles.resize(0);
}

CXXTorusElement::CXXTorusElement(const CXXCircle &aCircle, int iEdge, double delta, double probeRadius) :
   theCircle(aCircle), rProbe(probeRadius), debug(0)
{
   init();
   const CXXCoord<CXXCoord_ftype>xAxis(1.,0.,0.,0.);
   const CXXCoord<CXXCoord_ftype>yAxis(0.,1.,0.,0.);
   const CXXCoord<CXXCoord_ftype>zAxis(0.,0.,1.,0.);

   torusCentre = theCircle.getCentreOfCircle();
   torusAxis = theCircle.getNormal();
   rTraj = theCircle.getRadiusOfCircle();

   const CXXCircleNode &startNode(*theCircle.start(iEdge));
   const CXXCircleNode &endNode(*theCircle.stop(iEdge));
   v1unit = startNode.getUnitRadius();
   v2unit = endNode.getUnitRadius();
   CXXCoord<CXXCoord_ftype>temp1 = startNode.getCoord();

   n1unit = v1unit ^ torusAxis;

   CXXCoord<CXXCoord_ftype>diff = theCircle.getCentreOfSecondSphere()-temp1;
   diff.normalise();
   theta1 = acos(torusAxis*diff);
   diff = theCircle.getCentreOfSphere()-temp1;
   diff.normalise();
   theta2 = acos(torusAxis*diff);
   double halfWay1= (theta2+theta1)/2.;
   double halfWay2 = halfWay1;

   absoluteStartOmega = startNode.getAngle();

   omega1 = 0;
   if (startNode.getOtherCircle()) {
      omega2 = endNode.getAngle() - startNode.getAngle();
      while (omega2 < 0.) omega2 += 2.*M_PI;
   } else {
      omega2 = 2.*M_PI;
   }

   deltaOmega=omega2-omega1;
   deltaTheta=theta2-halfWay2;

   nOmegaSteps = 1;
   while (fabs(deltaOmega)>delta) {
      nOmegaSteps*=2;
      deltaOmega /= 2.;
   }

   nThetaSteps = 1;
   while (fabs(deltaTheta)>delta) {
      nThetaSteps*=2;
      deltaTheta/=2.;
   }

   if (rTraj <= probeRadius) {
      halfWay1 = asin(rTraj/probeRadius);
      halfWay2 = M_PI - halfWay1;
      nThetaSteps = int((fabs(theta2-halfWay2)) / deltaTheta) + 1;
   }

   nodes.resize((nThetaSteps+1)*(nOmegaSteps+1));

   double omega = omega1;
   for (int iOmega=0; iOmega<=nOmegaSteps; iOmega++) {
      double theta = theta2;
      for (int iTheta=0; iTheta<=nThetaSteps; iTheta++) {
         CXXTorusNode aNode(theta,omega);
         CXXCoord<CXXCoord_ftype>xyz(coordFromThetaOmega(theta, omega));
         aNode.setCoord(xyz);
         nodes[iTheta+(iOmega*(nThetaSteps+1))] = aNode;
         theta-= deltaTheta;
         if (deltaTheta > 0.) theta = (theta > halfWay2 ? theta :halfWay2);
         else theta = (theta < halfWay2 ? theta :halfWay2);
      }
      omega+=deltaOmega;
   }

   for (int iOmega=0; iOmega<nOmegaSteps; iOmega++) {
      for (int iTheta=0; iTheta<nThetaSteps; iTheta++) {
         size_t iTriangle = flatTriangles.size();
         flatTriangles.push_back(CXXTriangle(((iOmega+1)*(nThetaSteps+1))+(iTheta+0),
                                             ((iOmega+0)*(nThetaSteps+1))+(iTheta+0),
                                             ((iOmega+0)*(nThetaSteps+1))+(iTheta+1),
                                             iTriangle));
         if (iTheta == 0) edgeTriangles.push_back(&flatTriangles.back());
         iTriangle = flatTriangles.size();
         flatTriangles.push_back(CXXTriangle(((iOmega+0)*(nThetaSteps+1))+(iTheta+1),
                                             ((iOmega+1)*(nThetaSteps+1))+(iTheta+1),
                                             ((iOmega+1)*(nThetaSteps+1))+(iTheta+0),
                                             iTriangle));
      }
   }
   for (unsigned int iNode = 0; iNode<nodes.size(); iNode++) {
      nodes[iNode].setAtom(theCircle.getParent()->getAtomI());
   }
}

size_t CXXTorusElement::addNode(CXXTorusNode &aNode) {
   CXXTorusNode newNode(aNode);
   CXXCoord<CXXCoord_ftype>xyz = coordFromThetaOmega(newNode.getTheta(), newNode.getOmega());
   newNode.setCoord(xyz);
   nodes.push_back(newNode);
   return (nodes.size() - 1);
}

const size_t CXXTorusElement::numberOfTorusNodes(void) {
   return nodes.size();
}

const CXXCoord<CXXCoord_ftype> CXXTorusElement::probeAtOmega(double omega) const {
   CXXCoord<CXXCoord_ftype>temp1 = v1unit;
   temp1.scale(cos(omega));
   CXXCoord<CXXCoord_ftype>temp2 = n1unit;
   temp2.scale(sin(omega));
   CXXCoord<CXXCoord_ftype>temp3 = temp1 + temp2;
   temp3.scale(rTraj);
   return torusCentre + temp3;
}

const CXXCoord<CXXCoord_ftype> CXXTorusElement::normalToProbeAtTheta(CXXCoord<CXXCoord_ftype>&p, double theta) const {
   CXXCoord<CXXCoord_ftype>temp1 = torusAxis;
   temp1.scale(cos(theta));
   CXXCoord<CXXCoord_ftype>temp2 = torusCentre - p;
   temp2.normalise();
   temp2.scale(sin(theta));
   return temp1 + temp2;
}

const CXXCoord<CXXCoord_ftype> CXXTorusElement::coordFromThetaOmega(double theta, double omega) const {
   CXXCoord<CXXCoord_ftype>p(probeAtOmega(omega));
   CXXCoord<CXXCoord_ftype>normal(normalToProbeAtTheta(p, theta));
   normal.scale(rProbe);
   return  p + normal;
}

int CXXTorusElement::upload(CXXSurface *aSurface) {
   size_t oldVertexCount;
   {
      std::vector<double> verticesBuffer(nodes.size()*3);
      std::vector<double> accessiblesBuffer(nodes.size()*3);
      std::vector<double> normalsBuffer(nodes.size()*3);
      for (unsigned int i=0; i< nodes.size(); i++) {
         for (int j=0; j<3; j++) verticesBuffer[3*i+j] = nodes[i].coord().element(j);
         CXXCoord<CXXCoord_ftype>accessible = probeAtOmega(nodes[i].getOmega());
         for (int j=0; j<3; j++) accessiblesBuffer[3*i+j] = accessible.element(j);
         CXXCoord<CXXCoord_ftype>normal = normalToProbeAtTheta(accessible, nodes[i].getTheta());
         normal.scale(-1.);
         for (int j=0; j<3; j++) normalsBuffer[3*i+j] = normal.element(j);
      }
      oldVertexCount = aSurface->numberOfVertices();
      aSurface->updateWithVectorData(nodes.size(), "vertices", oldVertexCount, verticesBuffer.data());
      aSurface->updateWithVectorData(nodes.size(), "accessibles", oldVertexCount, accessiblesBuffer.data());
      aSurface->updateWithVectorData(nodes.size(), "normals",  oldVertexCount, normalsBuffer.data());
   }
   {
      std::vector<void *>atomBuffer(nodes.size());
      for (unsigned int i=0; i< nodes.size(); i++) {
         atomBuffer[i] = (void *)nodes[i].getAtom();
      }
      aSurface->updateWithPointerData(nodes.size(), "atom", oldVertexCount, atomBuffer.data());
   }
   {
      std::vector<int> triangleBuffer(flatTriangles.size()*3);
      int nToDraw = 0;
      list<CXXTriangle>::iterator trianglesEnd(flatTriangles.end());
      for (list<CXXTriangle>::iterator triangle=flatTriangles.begin(); triangle != trianglesEnd; ++triangle) {
         CXXTriangle &flatTriangle(*triangle);
         if (flatTriangle.doDraw()) {
            for (unsigned int j=0; j<3; j++) {
               triangleBuffer[(3*nToDraw)+j] = int(flatTriangle[j] + oldVertexCount);
            }
            nToDraw++;
         }
      }
      aSurface->extendTriangles(triangleBuffer.data(), nToDraw);
   }
   return 0;
}

const CXXTorusNode &CXXTorusElement::getNode(const int i) const {
   return nodes[i];
}

void CXXTorusElement::addEdgeVertex(CXXCircleNode &aNode) {
   double omega = aNode.getAngle() - absoluteStartOmega;
   while (omega < 0.) omega += 2.*M_PI;
   if (omega < omega2) {
      int triangleFound = 0;
      list<CXXTriangle *>::iterator matchingTriangle;
      list<CXXTriangle *>::iterator edgeTrianglesEnd = edgeTriangles.end();
      for (list<CXXTriangle *>::iterator triangle = edgeTriangles.begin();
           triangle != edgeTrianglesEnd && !triangleFound; ++triangle) {
         CXXTriangle &theFlatTriangle(**triangle);
         double omegaStart = nodes[theFlatTriangle[1]].getOmega();
         double omegaEnd   = nodes[theFlatTriangle[0]].getOmega();
         if (omega >= omegaStart && omega <= omegaEnd) {
            triangleFound = 1;
            matchingTriangle = triangle;
         }
      }

      if (triangleFound) {
         CXXTorusNode aNode2(theta2,omega);
         CXXCoord<CXXCoord_ftype>xyz(coordFromThetaOmega(theta2, omega));
         aNode2.setCoord(xyz);
         aNode2.setAtom(theCircle.getParent()->getAtomI());
         nodes.push_back(aNode2);

         CXXTriangle &theFT(**matchingTriangle);
         theFT.setDoDraw(0);
         edgeTriangles.erase(matchingTriangle);

         flatTriangles.push_back(CXXTriangle(theFT[0], nodes.size()-1, theFT[2]));
         edgeTriangles.push_back(&flatTriangles.back());
         flatTriangles.push_back(CXXTriangle(nodes.size()-1, theFT[1], theFT[2]));
         edgeTriangles.push_back(&flatTriangles.back());
      }
   }
}
