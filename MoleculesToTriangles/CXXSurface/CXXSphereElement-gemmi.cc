/*
 * MoleculesToTriangles/CXXSurface/CXXSphereElement-gemmi.cc
 *
 * gemmi-native twin of CXXSphereElement.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#define _USE_MATH_DEFINES
#include "CXXSurfaceVertex.h"
#include "CXXSphereElement-gemmi.hh"
#include "CXXCoord.h"
#include "CXXSphereTriangle-gemmi.hh"
#include "CXXSphereTriangleEdge-gemmi.hh"
#include "CXXSphereNode-gemmi.hh"
#include "CXXSphereFlatTriangle-gemmi.hh"
#include "CXXSurface-gemmi.hh"
#include "CXXCircleNode-gemmi.hh"
#include "CXXCircle-gemmi.hh"
#include "CXXNewHood-gemmi.hh"
#include "CXXTorusElement-gemmi.hh"
#include <math.h>

using namespace std;
using namespace coot::m2t;

void CXXSphereElement::init() {
   nDrawnTriangles=0;
}
CXXSphereElement::CXXSphereElement() : theAtom(0) {
   init();
}

CXXSphereElement::CXXSphereElement(const gemmi::Atom* anAtom, double del) :
   theAtom(anAtom), deltaRadians(del) {
   init();
   theCentre = CXXCoord<CXXCoord_ftype>(anAtom->pos.x, anAtom->pos.y, anAtom->pos.z);
   calculate();
}

CXXSphereElement::CXXSphereElement(const CXXCoord<CXXCoord_ftype>&position, double radius, double del) :
   theCentre(position), theRadius(radius), deltaRadians(del) {
   init();
   calculate();
   std::list<CXXSphereFlatTriangle>::iterator trianglesEnd = flatTriangles.end();
   for (std::list<CXXSphereFlatTriangle>::iterator triangle = flatTriangles.begin(); triangle!=trianglesEnd; ++triangle) {
      CXXSphereFlatTriangle &flatTriangle(*triangle);
      for (int j=0; j<3; j++) flatTriangle.setEdgeCircle(j,0);
   }
}

int CXXSphereElement::calculate() {
   theVertices.resize(0);
   theTriangles.resize(0);
   flatTriangles.resize(0);

   CXXCoord<CXXCoord_ftype>xAxis(1.,0.,0.,0.);
   CXXCoord<CXXCoord_ftype>yAxis(0.,1.,0.,0.);
   CXXCoord<CXXCoord_ftype>zAxis(0.,0.,1.,0.);
   CXXCoord<CXXCoord_ftype>normal1(zAxis);
   CXXCoord<CXXCoord_ftype>normal2(xAxis);
   CXXCoord<CXXCoord_ftype>normal3(yAxis);
   xAxis.scale(theRadius);
   yAxis.scale(theRadius);
   zAxis.scale(theRadius);
   CXXSphereNode v1(theCentre + xAxis);
   size_t iV1 = addVertex(v1);
   CXXSphereNode v2(theCentre + yAxis);
   size_t iV2 = addVertex(v2);
   CXXSphereNode v3(theCentre + zAxis);
   size_t iV3 = addVertex(v3);

   CXXSphereTriangleEdge edge1(normal1, iV1, iV2, theCentre, theCentre, theRadius, this);
   size_t iEdge1 = addEdge(edge1);
   CXXSphereTriangleEdge edge2(normal2, iV2, iV3, theCentre, theCentre, theRadius, this);
   size_t iEdge2 = addEdge(edge2);
   CXXSphereTriangleEdge edge3(normal3, iV3, iV1, theCentre, theCentre, theRadius, this);
   size_t iEdge3 = addEdge(edge3);
   size_t vertices[] = {iV1, iV2, iV3};
   size_t edges[] = {iEdge1, iEdge2, iEdge3};
   CXXSphereTriangle aTriangle(this, vertices, edges, theRadius, theCentre);
   addTriangle(aTriangle);

   while (theTriangles.size()) {
      if (!theTriangles.back().bisect(delta())) flattenLastTriangle();
   }
   normal1.scale(-1.);
   normal2.scale(-1.);
   normal3.scale(-1.);

   edge1 = CXXSphereTriangleEdge(normal1, iV2, iV1, theCentre, theCentre, theRadius, this);
   iEdge1 = addEdge(edge1);
   edge2 = CXXSphereTriangleEdge(normal2, iV3, iV2, theCentre, theCentre, theRadius, this);
   iEdge2 = addEdge(edge2);
   edge3 = CXXSphereTriangleEdge(normal3, iV1, iV3, theCentre, theCentre, theRadius, this);
   iEdge3 = addEdge(edge3);
   size_t vertices1[] = {iV1, iV3, iV2};
   size_t edges1[] = {iEdge3, iEdge2, iEdge1};
   aTriangle = CXXSphereTriangle(this, vertices1, edges1, theRadius, theCentre);
   addTriangle(aTriangle);

   while (theTriangles.size()) {
      if (!theTriangles.back().bisect(delta())) flattenLastTriangle();
   }
   nDrawnTriangles = countDrawnTriangles();
   return 0;
}

int CXXSphereElement::addTriangularPatch(const CXXCoord<CXXCoord_ftype>&u1_in, const CXXCoord<CXXCoord_ftype>&u2_in,
                                         const CXXCoord<CXXCoord_ftype>&u3_in, const gemmi::Atom* anAtom, vector<CXXCircle> &circles,
                                         int UseOrGenerate)
{
   theAtom = anAtom;
   CXXCoord<CXXCoord_ftype>u1(u1_in);
   CXXCoord<CXXCoord_ftype>u2(u2_in);
   CXXCoord<CXXCoord_ftype>u3(u3_in);

   CXXCoord<CXXCoord_ftype>normal1(u1^u2);
   normal1.normalise();
   if (normal1*u3<0.) {
      CXXCoord<CXXCoord_ftype>temp(u1);
      u1 = u2;
      u2 = temp;
      normal1 *= -1.;
   }
   CXXCoord<CXXCoord_ftype>normal2(u2^u3);
   normal2.normalise();
   CXXCoord<CXXCoord_ftype>normal3(u3^u1);
   normal3.normalise();

   CXXCoord<CXXCoord_ftype>v1Coord(u1);
   CXXCoord<CXXCoord_ftype>v2Coord(u2);
   CXXCoord<CXXCoord_ftype>v3Coord(u3);
   v1Coord *= theRadius;
   v2Coord *= theRadius;
   v3Coord *= theRadius;
   CXXSphereNode v1(theCentre + v1Coord);
   v1.setAtom(anAtom);
   size_t iV1 = addVertex(v1);
   CXXSphereNode v2(theCentre + v2Coord);
   v2.setAtom(anAtom);
   size_t iV2 = addVertex(v2);
   CXXSphereNode v3(theCentre + v3Coord);
   v3.setAtom(anAtom);
   size_t iV3 = addVertex(v3);

   circles.resize(3);
   CXXCircle &edgeCircle1(circles[0]);
   CXXCircle &edgeCircle2(circles[1]);
   CXXCircle &edgeCircle3(circles[2]);

   edgeCircle1.setNormal(normal1);
   edgeCircle1.setCentreToCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.,0.));
   edgeCircle1.setCentreOfCircle(theCentre);
   edgeCircle1.setRadiusOfCircle(theRadius);

   edgeCircle2.setNormal(normal2);
   edgeCircle2.setCentreToCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.,0.));
   edgeCircle2.setCentreOfCircle(theCentre);
   edgeCircle2.setRadiusOfCircle(theRadius);

   edgeCircle3.setNormal(normal3);
   edgeCircle3.setCentreToCircle(CXXCoord<CXXCoord_ftype>(0.,0.,0.,0.));
   edgeCircle3.setCentreOfCircle(theCentre);
   edgeCircle3.setRadiusOfCircle(theRadius);

   CXXSphereTriangleEdge edge1(normal1, iV1, iV2, theCentre, theCentre, theRadius, this);
   edge1.setCircle(&edgeCircle1);
   size_t iEdge1 = addEdge(edge1);

   CXXSphereTriangleEdge edge2(normal2, iV2, iV3, theCentre, theCentre, theRadius, this);
   edge2.setCircle(&edgeCircle2);
   size_t iEdge2 = addEdge(edge2);
   CXXSphereTriangleEdge edge3(normal3, iV3, iV1, theCentre, theCentre, theRadius, this);
   edge3.setCircle(&edgeCircle3);
   size_t iEdge3 = addEdge(edge3);

   size_t vertices[] = {iV1, iV2, iV3};
   size_t edges[] = {iEdge1, iEdge2, iEdge3};

   CXXSphereTriangle aTriangle(this, vertices, edges, theRadius, theCentre);
   aTriangle.setAtom(anAtom);
   addTriangle(aTriangle);

   while (theTriangles.size()) {
      if (!theTriangles.back().bisect(delta())) flattenLastTriangle();
   }
   clearCutFlags();
   nDrawnTriangles = countDrawnTriangles();
   return 0;
}

CXXSphereElement::CXXSphereElement(const CXXSphereElement &oldOne) :
   theAtom(oldOne.getAtom()),
   theCentre(oldOne.centre()),
   theVertices(oldOne.getVertices()),
   theTriangles(oldOne.getTriangles()),
   flatTriangles(oldOne.getFlatTriangles()),
   theCircles(oldOne.getCircles()),
   theRadius(oldOne.radius()),
   deltaRadians(oldOne.delta()),
   nDrawnTriangles(oldOne.getNDrawnTriangles())
{
   std::vector<CXXSphereTriangle>::iterator trianglesEnd = theTriangles.end();
   for (std::vector<CXXSphereTriangle>::iterator triangle = theTriangles.begin(); triangle!= trianglesEnd; ++triangle) {
      triangle->setSphereElement(this);
   }
   theEdges.resize(oldOne.nEdges());
   for (int i=0; i<oldOne.nEdges(); i++) {
      theEdges[i]=CXXSphereTriangleEdge(oldOne.edge(i));
      theEdges[i].setSphereElement(this);
   }
}

void CXXSphereElement::initWith(const CXXCircleNode &aNode, double delta,
                                double radius_in, bool *includeAtoms, int UseOrGenerate) {
   theCentre=aNode.getCoord();
   theRadius=radius_in;
   deltaRadians=delta;
   init();

   const gemmi::Atom* atomK=aNode.getAtomK();
   const gemmi::Atom* atomJ=aNode.getAtomJ();
   const gemmi::Atom* atomI=aNode.getAtomI();

   CXXCoord<CXXCoord_ftype>u1(atomK->pos.x, atomK->pos.y, atomK->pos.z);
   CXXCoord<CXXCoord_ftype>u2(atomJ->pos.x, atomJ->pos.y, atomJ->pos.z);
   CXXCoord<CXXCoord_ftype>u3(atomI->pos.x, atomI->pos.y, atomI->pos.z);

   u1 = u1 - theCentre;
   u1.normalise();
   u2 = u2-theCentre;
   u2.normalise();
   u3 = u3-theCentre;
   u3.normalise();

   CXXCoord<CXXCoord_ftype>u12;
   u12 = u1 + u2;
   u12.normalise();
   CXXCoord<CXXCoord_ftype>u23;
   u23 = u2 + u3;
   u23.normalise();
   CXXCoord<CXXCoord_ftype>u31;
   u31 = u3 + u1;
   u31.normalise();
   CXXCoord<CXXCoord_ftype>u123(u1 + u2 + u3);
   u123.normalise();

   int nIncludeAtoms = 0;
   for (int i=0; i<3; i++) if (includeAtoms[i]) nIncludeAtoms++;
   theCircles.resize(2*nIncludeAtoms);

   int iTriangularPatch = 0;
   if (includeAtoms[0]) {
      addTriangularPatch(u123, u3,  u23, atomI, theCircles[iTriangularPatch++], UseOrGenerate);
      addTriangularPatch(u123, u31, u3,  atomI, theCircles[iTriangularPatch++], UseOrGenerate);
   }
   if (includeAtoms[1]) {
      addTriangularPatch(u123, u2,  u12, atomJ, theCircles[iTriangularPatch++], UseOrGenerate);
      addTriangularPatch(u123, u23, u2,  atomJ, theCircles[iTriangularPatch++], UseOrGenerate);
   }
   if (includeAtoms[2]) {
      addTriangularPatch(u123, u1,  u31, atomK, theCircles[iTriangularPatch++], UseOrGenerate);
      addTriangularPatch(u123, u12, u1,  atomK, theCircles[iTriangularPatch++], UseOrGenerate);
   }
}

void CXXSphereElement::initWith(const CXXCoord<CXXCoord_ftype>&aCentre, const gemmi::Atom* atomI, const gemmi::Atom* atomJ, const gemmi::Atom* atomK,
                                double delta, double radius_in, const bool *includeAtoms) {
   int UseOrGenerate = CXXSphereElement::GenerateCircles;
   theCentre=aCentre;
   theRadius=radius_in;
   deltaRadians=delta;
   init();

   CXXCoord<CXXCoord_ftype>u1(atomK->pos.x, atomK->pos.y, atomK->pos.z);
   CXXCoord<CXXCoord_ftype>u2(atomJ->pos.x, atomJ->pos.y, atomJ->pos.z);
   CXXCoord<CXXCoord_ftype>u3(atomI->pos.x, atomI->pos.y, atomI->pos.z);

   u1 = u1 - theCentre;
   u1.normalise();
   u2 = u2-theCentre;
   u2.normalise();
   u3 = u3-theCentre;
   u3.normalise();

   CXXCoord<CXXCoord_ftype>u12;
   u12 = u1 + u2;
   u12.normalise();
   CXXCoord<CXXCoord_ftype>u23;
   u23 = u2 + u3;
   u23.normalise();
   CXXCoord<CXXCoord_ftype>u31;
   u31 = u3 + u1;
   u31.normalise();
   CXXCoord<CXXCoord_ftype>u123(u1 + u2 + u3);
   u123.normalise();

   int nIncludeAtoms = 0;
   for (int i=0; i<3; i++) if (includeAtoms[i]) nIncludeAtoms++;
   theCircles.resize(2*nIncludeAtoms);

   int iTriangularPatch = 0;
   if (includeAtoms[0]) {
      addTriangularPatch(u123, u3,  u23, atomI, theCircles[iTriangularPatch++], UseOrGenerate);
      addTriangularPatch(u123, u31, u3,  atomI, theCircles[iTriangularPatch++], UseOrGenerate);
   }
   if (includeAtoms[1]) {
      addTriangularPatch(u123, u2,  u12, atomJ, theCircles[iTriangularPatch++], UseOrGenerate);
      addTriangularPatch(u123, u23, u2,  atomJ, theCircles[iTriangularPatch++], UseOrGenerate);
   }
   if (includeAtoms[2]) {
      addTriangularPatch(u123, u1,  u31, atomK, theCircles[iTriangularPatch++], UseOrGenerate);
      addTriangularPatch(u123, u12, u1,  atomK, theCircles[iTriangularPatch++], UseOrGenerate);
   }
}

size_t CXXSphereElement::addVertex(const CXXSphereNode &vert) {
   size_t iVertex = theVertices.size();
   theVertices.push_back(vert);
   return iVertex;
}

size_t CXXSphereElement::addTriangle(const CXXSphereTriangle &newTri) {
   size_t iTriangle = theTriangles.size();
   theTriangles.push_back(newTri);
   return iTriangle;
}

size_t CXXSphereElement::addEdge(const CXXSphereTriangleEdge &newEdge) {
   size_t iEdge = theEdges.size();
   theEdges.push_back(newEdge);
   return iEdge;
}

void CXXSphereElement::flattenLastTriangle(void) {
   size_t i = theTriangles.back().vertex(0);
   size_t j = theTriangles.back().vertex(1);
   size_t k = theTriangles.back().vertex(2);
   CXXSphereFlatTriangle newTriangle(i, k, j, flatTriangles.size());
   newTriangle.setAtom(theTriangles.back().getAtom());
   newTriangle.setEdgeCircle(0, theEdges[theTriangles.back().edge(2)].getCircle());
   newTriangle.setEdgeCircle(1, theEdges[theTriangles.back().edge(1)].getCircle());
   newTriangle.setEdgeCircle(2, theEdges[theTriangles.back().edge(0)].getCircle());
   flatTriangles.push_back(newTriangle);
   theTriangles.pop_back();
}

const CXXCoord<CXXCoord_ftype>& CXXSphereElement::centre() const { return theCentre; }
const double CXXSphereElement::radius() const { return theRadius; }
const size_t CXXSphereElement::nVertices() const { return theVertices.size(); }
const size_t CXXSphereElement::nEdges() const { return theEdges.size(); }
const size_t CXXSphereElement::nFlatTriangles() const { return flatTriangles.size(); }
const CXXSphereNode &CXXSphereElement::vertex(const size_t iVertex) const { return theVertices[iVertex]; }
void CXXSphereElement::moveVertex(const int iVertex, const CXXCoord<CXXCoord_ftype>&position) { theVertices[iVertex].setVertex(position); }
const CXXSphereTriangle &CXXSphereElement::triangle(const size_t iTriangle) const { return theTriangles[iTriangle]; }
const CXXSphereTriangleEdge &CXXSphereElement::edge(const size_t iEdge) const { return theEdges[iEdge]; }

int CXXSphereElement::trimBy(const CXXCircle &aCircle, int carefully) {
   CXXSphereFlatTriangle newTriangle;
   int wasTrimmed = 0;
   vector<CXXCoord<CXXCoord_ftype> > theIntersections(2);

   list<CXXSphereFlatTriangle>newFlatTriangles;

   list<CXXSphereFlatTriangle>::iterator trianglesEnd = flatTriangles.end();
   for (list<CXXSphereFlatTriangle>::iterator triangle = flatTriangles.begin(); triangle != trianglesEnd; ++triangle) {
      CXXSphereFlatTriangle &flatTriangle(*triangle);
      if (flatTriangle.doDraw()) {
         size_t iNodes[] = {flatTriangle[0], flatTriangle[1], flatTriangle[2]};
         CXXSphereNode &vertex0(theVertices[iNodes[0]]);
         CXXSphereNode &vertex1(theVertices[iNodes[1]]);
         CXXSphereNode &vertex2(theVertices[iNodes[2]]);

         int isInFront[] = {
            aCircle.accIsBehind(vertex0.vertex()),
            aCircle.accIsBehind(vertex1.vertex()),
            aCircle.accIsBehind(vertex2.vertex())
         };

         int isInFrontMask = isInFront[0] + (isInFront[1]<<1) + (isInFront[2] << 2);
         int nInFront = isInFront[0] + isInFront[1] + isInFront[2];
         if (nInFront != 3) wasTrimmed=1;
         if (nInFront == 3) {
         }
         else if (nInFront == 0) {
            flatTriangle.setDoDraw(0);
            nDrawnTriangles--;
            for (int j=0; j<3; j++) theVertices[iNodes[j]].setDoDraw(0);
         }
         else if (nInFront == 1) {
            int i=0;
            if ((isInFrontMask & 1) == 1) i = 0;
            else if ((isInFrontMask & 2) == 2) i = 1;
            else if ((isInFrontMask & 4) == 4) i = 2;

            flatTriangle.setDoDraw(0);
            nDrawnTriangles--;

            theVertices[iNodes[(i+1)%3]].setDoDraw(0);
            theVertices[iNodes[(i+2)%3]].setDoDraw(0);

            CXXCoord<CXXCoord_ftype>newVertex1 = aCircle.accPlaneIntersect(theVertices[iNodes[i]].vertex(), theVertices[iNodes[(i+1)%3]].vertex());

            const CXXCircle *oldCircle = flatTriangle.getEdgeCircle(i);
            if (oldCircle &&carefully) {
               if (!aCircle.meetsCircle(*oldCircle, theIntersections)) {
                  CXXCoord<CXXCoord_ftype>diff1(theIntersections[0]- newVertex1);
                  CXXCoord<CXXCoord_ftype>diff2(theIntersections[1]- newVertex1);
                  newVertex1 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
               }
            }

            CXXSphereNode newNode1(newVertex1);
            newNode1.setAtom(theVertices[iNodes[i]].getAtom());
            iNodes[(i+1)%3] = addVertex(newNode1);
            theVertices[iNodes[(i+1)%3]].setIntersector(&aCircle);

            CXXCoord<CXXCoord_ftype>newVertex2 = aCircle.accPlaneIntersect(theVertices[iNodes[i]].vertex(), theVertices[iNodes[(i+2)%3]].vertex());

            const CXXCircle *oldCircle1 = flatTriangle.getEdgeCircle((i+2)%3);
            if (oldCircle1&&carefully) {
               if (!aCircle.meetsCircle(*oldCircle1, theIntersections)) {
                  CXXCoord<CXXCoord_ftype>diff1(theIntersections[0]- newVertex2);
                  CXXCoord<CXXCoord_ftype>diff2(theIntersections[1]- newVertex2);
                  newVertex2 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
               }
            }

            CXXSphereNode newNode2(newVertex2);
            newNode2.setAtom(theVertices[iNodes[i]].getAtom());
            iNodes[(i+2)%3] = addVertex(newNode2);
            theVertices[iNodes[(i+2)%3]].setIntersector(&aCircle);

            newTriangle[0] = iNodes[i];
            newTriangle[1] = iNodes[(i+1)%3];
            newTriangle[2] = iNodes[(i+2)%3];
            newTriangle.setAtom(flatTriangle.getAtom());
            newTriangle.setEdgeCircle(0, oldCircle);
            newTriangle.setEdgeCircle(1, &aCircle);
            newTriangle.setEdgeCircle(2, oldCircle1);

            newFlatTriangles.push_back(newTriangle);
            nDrawnTriangles++;
         }
         else {
            int i=0;
            if ((isInFrontMask & 1) != 1) i = 0;
            else if ((isInFrontMask & 2) != 2) i = 1;
            else if ((isInFrontMask & 4) != 4) i = 2;

            flatTriangle.setDoDraw(0);
            nDrawnTriangles--;

            theVertices[iNodes[i]].setDoDraw(0);

            CXXCoord<CXXCoord_ftype>newVertex1 = aCircle.accPlaneIntersect(theVertices[iNodes[(i+1)%3]].vertex(), theVertices[iNodes[i]].vertex());

            const CXXCircle *oldCircle = flatTriangle.getEdgeCircle(i);
            if (oldCircle&&carefully) {
               if (!aCircle.meetsCircle(*oldCircle, theIntersections)) {
                  CXXCoord<CXXCoord_ftype>diff1(theIntersections[0]- newVertex1);
                  CXXCoord<CXXCoord_ftype>diff2(theIntersections[1]- newVertex1);
                  newVertex1 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
               }
            }

            CXXSphereNode newNode1(newVertex1);
            newNode1.setAtom(theVertices[iNodes[i]].getAtom());
            size_t iNewNode1 = addVertex(newNode1);
            theVertices[iNewNode1].setIntersector(&aCircle);

            newTriangle[0] = iNewNode1;
            newTriangle[1] = iNodes[(i+1)%3];
            newTriangle[2] = iNodes[(i+2)%3];
            newTriangle.setAtom(flatTriangle.getAtom());
            newTriangle.setEdgeCircle(0, flatTriangle.getEdgeCircle(i));
            newTriangle.setEdgeCircle(1, flatTriangle.getEdgeCircle((i+1)%3));
            newTriangle.setEdgeCircle(2, 0);
            newFlatTriangles.push_back(newTriangle);
            nDrawnTriangles++;

            CXXCoord<CXXCoord_ftype>newVertex2 = aCircle.accPlaneIntersect(theVertices[iNodes[(i+2)%3]].vertex(), theVertices[iNodes[i]].vertex());

            const CXXCircle *oldCircle1 = flatTriangle.getEdgeCircle((i+2)%3);
            if (oldCircle1&&carefully) {
               if (!aCircle.meetsCircle(*oldCircle1, theIntersections)) {
                  CXXCoord<CXXCoord_ftype>diff1(theIntersections[0]- newVertex2);
                  CXXCoord<CXXCoord_ftype>diff2(theIntersections[1]- newVertex2);
                  newVertex2 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
               }
            }

            CXXSphereNode newNode2(newVertex2);
            newNode2.setAtom(theVertices[iNodes[i]].getAtom());
            size_t iNewNode2 = addVertex(newNode2);
            theVertices[iNewNode2].setIntersector(&aCircle);

            newTriangle[0] = iNewNode2;
            newTriangle[1] = iNewNode1;
            newTriangle[2] = iNodes[(i+2)%3];
            newTriangle.setAtom(flatTriangle.getAtom());
            newTriangle.setEdgeCircle(0, &aCircle);
            newTriangle.setEdgeCircle(1, 0);
            newTriangle.setEdgeCircle(2, flatTriangle.getEdgeCircle((i+2)%3));
            newFlatTriangles.push_back(newTriangle);
            nDrawnTriangles++;
         }
      }
   }
   flatTriangles.insert(flatTriangles.end(), newFlatTriangles.begin(), newFlatTriangles.end());
   flatTriangles.remove_if(CXXTriangle::doNotDraw);

   return wasTrimmed;
}

const double CXXSphereElement::delta() const { return deltaRadians; }

int CXXSphereElement::translateBy(const CXXCoord<CXXCoord_ftype>&crd) {
   theCentre = theCentre + crd;
   for (unsigned i=0; i<nVertices(); i++) {
      const CXXCoord<CXXCoord_ftype>result = theVertices[i].vertex() + crd;
      theVertices[i].setVertex(result);
   }
   std::vector<CXXSphereTriangle>::iterator trianglesEnd = theTriangles.end();
   for (std::vector<CXXSphereTriangle>::iterator triangle = theTriangles.begin(); triangle != trianglesEnd; ++triangle) {
      CXXSphereTriangle &theTriangle(*triangle);
      CXXCoord<CXXCoord_ftype>result = theTriangle.centre() + crd;
      theTriangle.setCentre(result);
   }
   for (int i=0; i<nEdges(); i++) {
      CXXCoord<CXXCoord_ftype>result = theEdges[i].edgeCentre() + crd;
      theEdges[i].setEdgeCentre(result);
      result = theEdges[i].sphereCentre() + crd;
      theEdges[i].setSphereCentre(result);
   }
   return 0;
}

int CXXSphereElement::scaleBy(const double factor) {
   double oldRadius = theRadius;
   theRadius = theRadius * factor;
   for (unsigned i=0; i<nVertices(); i++) {
      CXXCoord<CXXCoord_ftype>result = theVertices[i].vertex() - theCentre;
      result.scale(theRadius / oldRadius);
      result = result + theCentre;
      theVertices[i].setVertex(result);
   }
   return 0;
}

int CXXSphereElement::setAtom(const gemmi::Atom* anAtom) { theAtom = anAtom; return 0; }
const gemmi::Atom* CXXSphereElement::getAtom() const { return theAtom; }

int CXXSphereElement::addFlatTriangle(const CXXSphereFlatTriangle &theTriangle) {
   flatTriangles.push_back(theTriangle);
   return 0;
}

CXXCoord<CXXCoord_ftype> CXXSphereElement::voronoiPoint(const CXXCoord<CXXCoord_ftype>&a, const CXXCoord<CXXCoord_ftype>&b, const CXXCoord<CXXCoord_ftype>&c) const {
   CXXCoord<CXXCoord_ftype>edge1(b-a);
   CXXCoord<CXXCoord_ftype>edge2(c-a);
   CXXCoord<CXXCoord_ftype>planeNormal(edge1^edge2);
   CXXCoord<CXXCoord_ftype>trivialResult(a+b+c);
   trivialResult /= 3.;
   if (planeNormal.get3DLengthSq()<1e-8) return trivialResult;
   planeNormal.normalise();
   CXXCoord<CXXCoord_ftype>edgeNormal1(planeNormal^edge1);
   CXXCoord<CXXCoord_ftype>edgeNormal2(planeNormal^edge2);

   CXXCoord<CXXCoord_ftype>midPoint1(edge1);
   midPoint1 /= 2.;
   midPoint1 = midPoint1 + a;

   CXXCoord<CXXCoord_ftype>midPoint2(edge2);
   midPoint2 /= 2.;
   midPoint2 = midPoint2 + a;

   double k2 = (edgeNormal1[1]*(midPoint2[0]-midPoint1[0]) - edgeNormal1[0]*(midPoint2[1]-midPoint1[1])) /
      (edgeNormal1[0]*edgeNormal2[1]-edgeNormal2[0]*edgeNormal1[1]);

   CXXCoord<CXXCoord_ftype>result(edgeNormal2);
   result *= k2;
   result = result + midPoint2;
   CXXCoord<CXXCoord_ftype>diff(result - a);
   double c1 = diff*edge1;
   double c2 = diff*edge2;
   if (c1> 0. && c1 < 1. && c2 > 0. && c2 < 1. && c1 + c2 < 1.) return result;
   else return trivialResult;
}

void CXXSphereElement::identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXBall*, vector<CXXCoord<CXXCoord_ftype> > >&raggedEdges) {
   bool noneForThisCircle = true;
   vector<CXXCoord<CXXCoord_ftype> >*raggedEdgesOfBall = 0;
   std::list<CXXSphereFlatTriangle>::iterator trianglesEnd = flatTriangles.end();
   for (std::list<CXXSphereFlatTriangle>::iterator triangle = flatTriangles.begin(); triangle!=trianglesEnd; ++triangle) {
      CXXSphereFlatTriangle &flatTriangle(*triangle);
      if (flatTriangle.doDraw()) {
         int spun = 0;
         for (int j=0; j<3 && !spun; j++) {
            if (flatTriangle.getEdgeCircle(j) == &theCircle) {
               spun = 1;
               if (noneForThisCircle) raggedEdgesOfBall = &raggedEdges[theCircle.getBallJ()];
               const CXXCoord<CXXCoord_ftype>&v1(theVertices[flatTriangle[ j     ]].vertex());
               raggedEdgesOfBall->push_back(v1);
               const CXXCoord<CXXCoord_ftype>&v2(theVertices[flatTriangle[(j+1)%3]].vertex());
               raggedEdgesOfBall->push_back(v2);
            }
         }
      }
   }
}

size_t CXXSphereElement::flagCutTriangles(const CXXCircle &theCircle) {
   clearCutFlags();
   std::list<CXXSphereFlatTriangle>::iterator trianglesEnd = flatTriangles.end();
   for (std::list<CXXSphereFlatTriangle>::iterator triangle = flatTriangles.begin(); triangle!=trianglesEnd; ++triangle) {
      CXXSphereFlatTriangle &flatTriangle(*triangle);
      if (flatTriangle.doDraw()) {
         int spun = 0;
         for (int j=0; j<3 && !spun; j++) {
            if (flatTriangle.getEdgeCircle(j) == &theCircle) {
               spun = 1;
               edgeTriangles[&flatTriangle] = j;
               const CXXCoord<CXXCoord_ftype>&v1(theVertices[flatTriangle[ j     ]].vertex());
               CXXCircleNode cn1(&theCircle, 0, v1, 0);
               flatTriangle.setCircleNode( j     , cn1);
               const CXXCoord<CXXCoord_ftype>&v2(theVertices[flatTriangle[(j+1)%3]].vertex());
               CXXCircleNode cn2(&theCircle, 0, v2, 0);
               flatTriangle.setCircleNode((j+1)%3, cn2);
            }
         }
      }
   }
   return edgeTriangles.size();
}

void CXXSphereElement::addCircleVertices(const CXXCircle &theCircle, int iEdge, double delta) {
   const CXXCircleNode &node1(*theCircle.start(iEdge));
   const CXXCircleNode &node2(*theCircle.stop(iEdge));
   double startOmega = node1.getAngle();
   double stopOmega = node2.getAngle();
   double deltaOmega = stopOmega - startOmega;
   int nSteps = 1;
   while (deltaOmega > delta/2) {
      deltaOmega /= 2.;
      nSteps *= 2;
   }
   CXXCoord<CXXCoord_ftype>normalUnitRadius = theCircle.getReferenceUnitRadius()^theCircle.getNormal();
   double omega = startOmega;
   int nPlaced = 0;
   for (int iStep = 0; iStep <= nSteps+1; iStep++) {
      CXXCoord<CXXCoord_ftype>newNode = (theCircle.getReferenceUnitRadius()*cos(omega) + normalUnitRadius*sin(omega)) * theCircle.getRadiusOfCircle();
      CXXCoord<CXXCoord_ftype>translatedNewNode = newNode + theCircle.getCentreOfCircle();
      CXXCircleNode extraNode(&theCircle, 0, translatedNewNode, 0);
      extraNode.setAngle(omega);
      if (!addVertex(extraNode)) nPlaced++;
      omega += deltaOmega;
   }
   if (nPlaced != nSteps) std::cout << "Placed "<< nPlaced << " of " << nSteps+1<<std::endl;
}

void CXXSphereElement::addTorusVertices(const CXXTorusElement &aTorus) {
   for (int iStep=0; iStep <= aTorus.getNOmegaSteps(); iStep++) {
      double omega = aTorus.getAbsoluteStartOmega() + double(iStep) * aTorus.getDeltaOmega();
      while (omega < 0.) omega += 2.*M_PI;
      while (omega > 2.*M_PI) omega -= 2.*M_PI;
      const CXXCoord<CXXCoord_ftype>newNode(aTorus.probeAtOmega(omega-aTorus.getAbsoluteStartOmega()));
      CXXCircleNode extraNode(&aTorus.getCircle(), 0, newNode, 0);
      extraNode.setAngle(omega);
      addVertex(extraNode);
   }
}

int CXXSphereElement::addVertex(const CXXCircleNode &extraNode) {
   const CXXCoord<CXXCoord_ftype>&newNode = extraNode.getCoord();
   const CXXCircle *aCircle = extraNode.getParent();
   map<CXXSphereFlatTriangle *, int>::iterator cutTriangle;
   map<CXXSphereFlatTriangle *, int>::iterator edgeTrianglesEnd = edgeTriangles.end();
   map<CXXSphereFlatTriangle *, int>::iterator removeTriangle;
   map<CXXSphereFlatTriangle *, int> extraEdgeTriangles;
   int cutFound = 0;
   for (cutTriangle=edgeTriangles.begin(); cutTriangle!=edgeTrianglesEnd && !cutFound; ++cutTriangle) {
      CXXSphereFlatTriangle &flatTriangle(*cutTriangle->first);
      int iEdge = cutTriangle->second;
      if (flatTriangle.doDraw()) {
         if (flatTriangle.getEdgeCircle(iEdge) == aCircle) {
            const CXXCircleNode &cn1 = flatTriangle.getCircleNode( iEdge     );
            const CXXCircleNode &cn2 = flatTriangle.getCircleNode((iEdge+1)%3);
            const CXXSphereNode &vert1(vertex(flatTriangle[iEdge]));
            if (aCircle->smallabBracketsC(cn1, cn2, extraNode)) {
               cutFound = 1;
               removeTriangle = cutTriangle;
               CXXSphereNode newSphereNode(newNode);
               newSphereNode.setAtom(vert1.getAtom());
               size_t iv3 = addVertex(newSphereNode);

               flatTriangles.push_back(CXXSphereFlatTriangle(flatTriangle));
               flatTriangles.back().setElement   ((iEdge+1)%3, iv3);
               flatTriangles.back().setCircleNode((iEdge+1)%3, extraNode);
               nDrawnTriangles++;
               extraEdgeTriangles[&flatTriangles.back()] = iEdge;

               flatTriangles.push_back(CXXSphereFlatTriangle(flatTriangle));
               flatTriangles.back().setElement   ( iEdge,      iv3);
               flatTriangles.back().setCircleNode( iEdge,      extraNode);
               nDrawnTriangles++;
               extraEdgeTriangles[&flatTriangles.back()] = iEdge;

               flatTriangle.setDoDraw(0);
               nDrawnTriangles--;
            }
         }
      }
   }
   if (cutFound) {
      edgeTriangles.erase(removeTriangle);
      edgeTriangles.insert(extraEdgeTriangles.begin(), extraEdgeTriangles.end());
      return 0;
   }
   else return 1;
}

int CXXSphereElement::countDrawnTriangles() const {
   int nToDraw = 0;
   std::list<CXXSphereFlatTriangle>::const_iterator trianglesEnd = flatTriangles.end();
   for (std::list<CXXSphereFlatTriangle>::const_iterator triangle = flatTriangles.begin(); triangle != trianglesEnd; ++triangle) {
      if (triangle->doDraw()) nToDraw++;
   }
   return nToDraw;
}
