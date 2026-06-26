/*
 * MoleculesToTriangles/CXXSurface/CXXSurface-gemmi.cc
 *
 * gemmi-native twin of CXXSurface.cpp. The grasp-file I/O and mmdb-UDD colour
 * helpers (readGraspFile/writeAsGrasp/colorByColourArray/getIntegerUDDataOfAtom)
 * are not ported (electrostatics/Phase 6 / file I/O not needed for the surface
 * mesh path). The per-vertex atom rides the opaque void* "atom" pointer property.
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
#include <math.h>
#include <set>
#include "string.h"
#include "CXXSurface-gemmi.hh"
#include "CXXNewHood-gemmi.hh"
#include "CXXTorusElement-gemmi.hh"
#include "CXXTorusNode-gemmi.hh"
#include "CXXSphereElement-gemmi.hh"
#include "CXXTriangle-gemmi.hh"
#include "CXXSphereFlatTriangle-gemmi.hh"
#include "CXXSurfaceVertex.h"
#include <limits>
#include <stdint.h>

using namespace std;
using namespace coot::m2t;

SurfaceParameters CXXSurface::measuredProperties()
{
   double area = 0.;
   int i;
   CXXCoord<CXXCoord_ftype>AB, AC, ABxAC;
   SurfaceParameters result;

   result.nVertices = vertices.size();
   result.nTriangles = triangles.size();

   if (vectors.find("vertices") != vectors.end()) {
      for (i=0, area = 0.; i<(int)nTriangles; i++) {
         const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["vertices"], i, 0);
         const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["vertices"], i, 1);
         const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["vertices"], i, 2);
         AB = coordA - coordB;
         AC = coordC - coordB;
         ABxAC = AB ^ AC;
         area += 0.5 * fabs(ABxAC.get3DLength());
      }
      result.MSA = area;
   }

   if (vectors.find("accessibles") != vectors.end()) {
      for (i=0, area = 0.; i<(int)nTriangles; i++) {
         const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["accessibles"], i, 0);
         const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["accessibles"], i, 1);
         const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["accessibles"], i, 2);
         AB = coordA - coordB;
         AC = coordC - coordB;
         ABxAC = AB ^ AC;
         area += 0.5 * fabs(ABxAC.get3DLength());
      }
      result.ASA = area;
   }

   map<string, size_t>::iterator scalar;
   int j;
   for (scalar = scalars.begin(), j=0; scalar != scalars.end(); ++scalar, j++) {
      double potential, pMin, pMax, pMean;
      for (i=0, pMin = 1e30, pMax = -1e30, pMean = 0.; i<int(vertices.size()); i++) {
         potential = vertices[i].scalar(scalars[scalar->first]);
         pMin = (potential<pMin?potential:pMin);
         pMax = (potential>pMax?potential:pMax);
         pMean += potential;
      }
      pMean /= vertices.size();
      result.pMins[scalar->first] = pMin;
      result.pMins[scalar->first] = pMax;
      result.pMins[scalar->first] = pMean;
   }
   return result;
}

std::string CXXSurface::report() {
   double area = 0.;
   double pMin, pMax, pMean;
   int i;
   CXXCoord<CXXCoord_ftype>AB, AC, ABxAC;
   std::ostringstream output;

   output << "This surface has " << vertices.size() << " vertices and " << nTriangles << " triangles\n";

   if (vectors.find("vertices") != vectors.end()) {
      for (i=0, area = 0.; i<(int)nTriangles; i++) {
         const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["vertices"], i, 0);
         const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["vertices"], i, 1);
         const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["vertices"], i, 2);
         AB = coordA - coordB;
         AC = coordC - coordB;
         ABxAC = AB ^ AC;
         area += 0.5 * fabs(ABxAC.get3DLength());
      }
      output << "The Molecular surface area is " << area << endl;
   }

   if (vectors.find("accessibles") != vectors.end()) {
      for (i=0, area = 0.; i<(int)nTriangles; i++) {
         const CXXCoord<CXXCoord_ftype>&coordA = coordRef(vectors["accessibles"], i, 0);
         const CXXCoord<CXXCoord_ftype>&coordB = coordRef(vectors["accessibles"], i, 1);
         const CXXCoord<CXXCoord_ftype>&coordC = coordRef(vectors["accessibles"], i, 2);
         AB = coordA - coordB;
         AC = coordC - coordB;
         ABxAC = AB ^ AC;
         area += 0.5 * fabs(ABxAC.get3DLength());
      }
      output << "The Solvent Accessible surface area is " << area << endl;
   }

   map<string, size_t>::iterator scalar;
   int j;
   for (scalar = scalars.begin(), j=0; scalar != scalars.end(); ++scalar, j++) {
      double potential;
      for (i=0, pMin = 1e30, pMax = -1e30, pMean = 0.; i<int(vertices.size()); i++) {
         potential = vertices[i].scalar(scalars[scalar->first]);
         pMin = (potential<pMin?potential:pMin);
         pMax = (potential>pMax?potential:pMax);
         pMean += potential;
      }
      pMean /= vertices.size();
      output << "Property Number " << scalar->second << " Range " << pMin << " to " << pMax << " mean " << pMean << endl;
   }
   return output.str();
}

int CXXSurface::assignAtom(gemmi::Structure *structure, const std::set<const gemmi::Atom*> &selSet) {
   (void) structure;
   std::vector<void*> pointerBuffer(vertices.size());
   for (int i=0; i< int(vertices.size()); i++) {
      const CXXCoord<CXXCoord_ftype>&vertex = coordRef(vectors["vertices"], i);
      double minDistSq = 1e30;
      const gemmi::Atom* best = 0;
      for (const gemmi::Atom* a : selSet) {
         CXXCoord<CXXCoord_ftype>atom(a->pos.x, a->pos.y, a->pos.z);
         CXXCoord<CXXCoord_ftype>diff = atom - vertex;
         double dxsq = diff.x() * diff.x();
         if (dxsq<minDistSq) {
            double dysq = diff.y() * diff.y();
            if (dysq<minDistSq) {
               double dzsq = diff.z() * diff.z();
               if (dzsq<minDistSq) {
                  if (dxsq + dysq + dzsq < minDistSq) {
                     best = a;
                     minDistSq = dxsq + dysq + dzsq;
                  }
               }
            }
         }
      }
      pointerBuffer[i] = (void *) best;
   }
   addPerVertexPointer("atom", pointerBuffer.data());
   return 0;
}

int CXXSurface::colorByAssignedAtom() {
   double redColour[] = {1.,0.,0.};
   double greenColour[] = {0.,1.,0.};
   double blueColour[] = {0.,0.,1.};
   double yellowColour[] = {1.,1.,0.};
   double greyColour[] = {1.,0.,0.};

   double *colourBuffer = new double[3*vertices.size()];
   for (int i=0; i< int(vertices.size()); i++) {
      for (int j=0; j<3; j++) colourBuffer[3*i + j] = greyColour[j];
   }
   updateWithVectorData(vertices.size(), "colour", 0, colourBuffer);
   delete [] colourBuffer;

   for (int i=0; i< int(vertices.size()); i++) {
      const gemmi::Atom* theAtom = (const gemmi::Atom*) vertices[i].pointer(pointers["atom"]);
      if (theAtom != 0) {
         char e = theAtom->element.name()[0];
         switch (e) {
            case 'C': vertices[i].setXyz(vectors["colour"], greenColour);  break;
            case 'N': vertices[i].setXyz(vectors["colour"], blueColour);   break;
            case 'O': vertices[i].setXyz(vectors["colour"], redColour);    break;
            case 'S': vertices[i].setXyz(vectors["colour"], yellowColour); break;
         }
      }
      else {
         vertices[i].setXyz(vectors["colour"], greyColour);
      }
   }
   return 0;
}

size_t CXXSurface::extendWithVectorData(size_t count, const string name, double *vectorBuffer) {
   size_t start = vertices.size();
   updateWithVectorData(count, name, start, vectorBuffer);
   return vertices.size();
}

size_t CXXSurface::updateWithVectorData(size_t count, const string name, size_t start, double *vectorBuffer) {
   size_t iVector  = getVectorHandle(name);
   if (vertices.size() < start + count) vertices.resize(start+count);
   for (unsigned int i=0; i<(unsigned int) count; i++) {
      vertices[i+start].setXyz(iVector, &(vectorBuffer[3*i]));
   }
   return vertices.size();
}

size_t CXXSurface::updateWithPointerData(size_t count, const string name, size_t start, void* *pointerBuffer) {
   size_t iPointer = getPointerHandle(name);
   if (vertices.size() < start + count) vertices.resize(start+count);
   for (size_t i=0; i< count; i++) {
      vertices[i+start].setPointer(iPointer, pointerBuffer[i]);
   }
   return vertices.size();
}

size_t CXXSurface::extendTriangles(int *triangleBuffer, int count) {
   triangles.resize(nTriangles+count);
   for (int i=0; i<count; i++) {
      int l = triangleBuffer[3*i];
      int m = triangleBuffer[3*i+1];
      int n = triangleBuffer[3*i+2];
      triangles[i+nTriangles] = CXXTriangle(l, m, n, i+nTriangles);
   }
   nTriangles = triangles.size();
   return nTriangles;
}

size_t CXXSurface::numberOfTriangles() const { return triangles.size(); }
size_t CXXSurface::numberOfVertices() const { return vertices.size(); }
size_t CXXSurface::vertex(size_t iTriangle, size_t iCorner) const { return triangles[iTriangle][iCorner]; }

int CXXSurface::upLoadSphere(CXXSphereElement &theSphere, double probeRadius, const int sense) {
   size_t oldVertexCount;
   CXXCoord<CXXCoord_ftype>theCentre = theSphere.centre();

   vector<int> equivalence(theSphere.nVertices());
   vector<int> uniqueAndDrawn(theSphere.nVertices());

   int nDrawn = 0;
   for (unsigned  i=0; i< theSphere.nVertices(); i++) {
      uniqueAndDrawn[i] = 0;
      if (theSphere.vertex(i).doDraw()) {
         uniqueAndDrawn[i] = 1;
         if (uniqueAndDrawn[i]) equivalence[i] = nDrawn++;
      }
   }
   static const std::string vertexName("vertices");
   static const std::string accessiblesName("accessibles");
   static const std::string normalsName("normals");
   {
      oldVertexCount = numberOfVertices();
      vertices.resize(oldVertexCount+nDrawn);
      size_t verticesHandle = getVectorHandle(vertexName);
      size_t accessiblesHandle = getVectorHandle(accessiblesName);
      size_t normalsHandle = getVectorHandle(normalsName);
      int iDraw = 0;
      for (unsigned int i=0; i< theSphere.nVertices(); i++) {
         if (uniqueAndDrawn[i]) {
            CXXCoord<CXXCoord_ftype>vertexCoord = theSphere.vertex(i).vertex();
            if (sense == CXXSphereElement::Contact) {
               vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord);
               CXXCoord<CXXCoord_ftype>normal = vertexCoord - theCentre;
               CXXCoord<CXXCoord_ftype>diff(normal);
               diff *= (theSphere.radius() - probeRadius) / theSphere.radius();
               normal.normalise();
               vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
               CXXCoord<CXXCoord_ftype>vertex = theCentre + diff;
               vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertex);
            }
            else if (sense == CXXSphereElement::Reentrant) {
               vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
               CXXCoord<CXXCoord_ftype>normal = theCentre - vertexCoord;
               normal.normalise();
               vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
               vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, theCentre);
            }
            else if (sense == CXXSphereElement::VDW) {
               vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
               CXXCoord<CXXCoord_ftype>normal = vertexCoord - theCentre;
               normal.normalise();
               vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
               vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord+normal);
            }
            else if (sense == CXXSphereElement::Accessible) {
               vertices[oldVertexCount+iDraw].setCoord(verticesHandle, vertexCoord);
               CXXCoord<CXXCoord_ftype>normal = vertexCoord - theCentre;
               normal.normalise();
               vertices[oldVertexCount+iDraw].setCoord(normalsHandle, normal);
               vertices[oldVertexCount+iDraw].setCoord(accessiblesHandle, vertexCoord);
            }
            iDraw++;
         }
      }
   }

   {
      std::vector<void*> atomBuffer(nDrawn);
      int iDraw = 0;
      for (unsigned int i=0; i< theSphere.nVertices(); i++) {
         if (uniqueAndDrawn[i]) {
            const gemmi::Atom* anAtom;
            if ((anAtom = theSphere.vertex(i).getAtom())!=0) {
               atomBuffer[iDraw] = (void*) anAtom;
            }
            else atomBuffer[iDraw] = (void*) theSphere.getAtom();
            iDraw++;
         }
      }
      updateWithPointerData(nDrawn, "atom", oldVertexCount, (void **)&(atomBuffer[0]));
   }

   {
      std::vector<int> triangleBuffer(theSphere.nFlatTriangles()*3);
      int drawCount = 0;
      std::list<CXXSphereFlatTriangle>::const_iterator trianglesEnd = theSphere.getFlatTriangles().end();
      for (std::list<CXXSphereFlatTriangle>::const_iterator triangle = theSphere.getFlatTriangles().begin();
           triangle != trianglesEnd; ++triangle) {
         const CXXSphereFlatTriangle &theTriangle(*triangle);
         if (theTriangle.doDraw()) {
            if (sense == CXXSphereElement::Contact || sense == CXXSphereElement::VDW || sense == CXXSphereElement::Accessible) {
               for (unsigned int j=0; j<3; j++) {
                  int index = equivalence[theTriangle[2-j]];
                  triangleBuffer[3*drawCount+j] = int(index + oldVertexCount);
               }
            }
            else {
               for (unsigned int j=0; j<3; j++) {
                  int index = equivalence[theTriangle[j]];
                  triangleBuffer[3*drawCount+j] = int(index + oldVertexCount);
               }
            }
            drawCount++;
         }
      }
      extendTriangles((int *) &(triangleBuffer[0]), drawCount);
   }
   return 0;
}

int CXXSurface::uploadTorus(CXXTorusElement &theTorus) {
   size_t oldVertexCount;
   {
      std::vector<double> verticesBuffer(theTorus.nTorusNodes()*3);
      std::vector<double> accessiblesBuffer(theTorus.nTorusNodes()*3);
      std::vector<double> normalsBuffer(theTorus.nTorusNodes()*3);
      for (unsigned int i=0; i< theTorus.nTorusNodes(); i++) {
         for (int j=0; j<3; j++) verticesBuffer[3*i+j] = theTorus.node(i).coord().element(j);
         CXXCoord<CXXCoord_ftype>accessible = theTorus.probeAtOmega(theTorus.node(i).getOmega());
         for (int j=0; j<3; j++) accessiblesBuffer[3*i+j] = accessible.element(j);
         CXXCoord<CXXCoord_ftype>normal = theTorus.normalToProbeAtTheta(accessible, theTorus.node(i).getTheta());
         normal.scale(-1.);
         for (int j=0; j<3; j++) normalsBuffer[3*i+j] = normal.element(j);
      }
      oldVertexCount = numberOfVertices();
      updateWithVectorData(theTorus.nTorusNodes(), "vertices", oldVertexCount, verticesBuffer.data());
      updateWithVectorData(theTorus.nTorusNodes(), "accessibles", oldVertexCount, accessiblesBuffer.data());
      updateWithVectorData(theTorus.nTorusNodes(), "normals",  oldVertexCount, normalsBuffer.data());
   }
   {
      std::vector<void *>atomBuffer(theTorus.nTorusNodes());
      for (unsigned int i=0; i< theTorus.nTorusNodes(); i++) {
         atomBuffer[i] = (void *)theTorus.node(i).getAtom();
      }
      updateWithPointerData(theTorus.nTorusNodes(), "atom", oldVertexCount, atomBuffer.data());
   }
   {
      std::vector<int> triangleBuffer(theTorus.nFlatTriangles()*3);
      int nToDraw = 0;
      for (list<CXXTriangle>::const_iterator triangle = theTorus.firstTriangle(); triangle != theTorus.endOfTriangles(); ++triangle) {
         const CXXTriangle &flatTriangle(*triangle);
         if (flatTriangle.doDraw()) {
            for (unsigned int j=0; j<3; j++) {
               triangleBuffer[(3*nToDraw)+j] = int(flatTriangle[j] + oldVertexCount);
            }
            nToDraw++;
         }
      }
      extendTriangles(triangleBuffer.data(), nToDraw);
   }
   return 0;
}

gemmi::Structure* CXXSurface::getStructure() const { return allAtomsManager; }

size_t CXXSurface::setCoord(const string &name, size_t iVertex, const CXXCoord<CXXCoord_ftype>&crd) {
   size_t iVector = vectors[name];
   if (vertices.size() <= iVertex) vertices.resize(iVertex+1);
   vertices[iVertex].setCoord(iVector, crd);
   return vertices.size();
}

size_t CXXSurface::getVectorHandle(const string name) {
   if (vectors.find(name) == vectors.end()) {
      size_t oldSize = vectors.size();
      vectors[name] = oldSize + 1;
   }
   return vectors[name];
}

size_t CXXSurface::getReadVectorHandle(const string name) {
   if (vectors.find(name) == vectors.end()) return SIZE_MAX;
   return vectors[name];
}

size_t CXXSurface::getScalarHandle(const string name) {
   if (scalars.find(name) == scalars.end()) {
      size_t oldSize = scalars.size();
      scalars[name] = oldSize + 1;
   }
   return scalars[name];
}

void CXXSurface::setScalar(size_t scalarHandle, size_t iVertex, double &value) {
   vertices[iVertex].setScalar(scalarHandle, value);
}

void CXXSurface::setScalar(const string name, size_t iVertex, double &value) {
   size_t scalarHandle=getScalarHandle(name);
   vertices[iVertex].setScalar(scalarHandle, value);
}

size_t CXXSurface::getReadScalarHandle(const string name) {
   if (scalars.find(name) == scalars.end()) return SIZE_MAX;
   return scalars[name];
}

size_t CXXSurface::getPointerHandle(const string name) {
   if (pointers.find(name) == pointers.end()) {
      size_t oldSize = pointers.size();
      pointers[name] = oldSize + 1;
   }
   return pointers[name];
}

const CXXCoord<CXXCoord_ftype>& CXXSurface::coordRef(size_t coordType, size_t iTriangle, size_t iCorner) const {
   const CXXTriangle &theTriangle(triangles[iTriangle]);
   size_t iVertex = theTriangle[iCorner];
   return vertices[iVertex].coordRef(coordType);
}

const CXXCoord<CXXCoord_ftype>& CXXSurface::coordRef(size_t coordType, size_t iVertex) const {
   return vertices[iVertex].coordRef(coordType);
}

int CXXSurface::addPerVertexVector(const string name, double *coordBuffer) {
   size_t vectorHandle = getVectorHandle(name);
   for (int i=0; i<int(vertices.size()); i++) vertices[i].setXyz(vectorHandle, &(coordBuffer[3*i]));
   return 0;
}

int CXXSurface::addPerVertexScalar(const string name, double *scalarBuffer) {
   size_t scalarHandle = getScalarHandle(name);
   for (int i=0; i<int(vertices.size()); i++) vertices[i].setScalar(scalarHandle, scalarBuffer[i]);
   return 0;
}

int CXXSurface::addPerVertexPointer(const string name, void **pointerBuffer) {
   size_t pointerHandle = getPointerHandle(name);
   for (int i=0; i<int(vertices.size()); i++) vertices[i].setPointer(pointerHandle, pointerBuffer[i]);
   return 0;
}

int CXXSurface::getCoord(const string &type, const size_t iTriangle, const size_t corner, double *buffer) {
   size_t relevantVector = getReadVectorHandle(type);
   if (relevantVector==SIZE_MAX) return 1;
   size_t iVertex = triangles[iTriangle][corner];
   return getCoord(relevantVector, iVertex, buffer);
}

int CXXSurface::getCoord(const string &type, const size_t iVertex, double *buffer) {
   size_t relevantVector = getReadVectorHandle(type);
   if (relevantVector==SIZE_MAX) return 1;
   return getCoord(relevantVector, iVertex, buffer);
}

int CXXSurface::getCoord(const size_t type, const size_t iVertex, double *buffer) {
   if (type <= vectors.size()) {
      const CXXCoord<CXXCoord_ftype>&theCoord = coordRef(type, iVertex);
      for (int i=0; i<4; i++) buffer[i] = theCoord.element(i);
      return 0;
   }
   else return 1;
}

int CXXSurface::getPointer(const string &type, size_t iVertex, void **return_p) {
   if (pointers.find(type) != pointers.end()) {
      *return_p = vertices[iVertex].pointer(pointers[type]);
      return 0;
   }
   *return_p = 0;
   return 1;
}

int CXXSurface::getScalar(int handle, int iVertex, double &result) {
   result = vertices[iVertex].scalar(handle);
   return 0;
}

int CXXSurface::addTriangle(const CXXTriangle &aTriangle) {
   triangles.push_back(aTriangle);
   return 0;
}

void CXXSurface::appendSurface(const CXXSurface &otherSurface) {
   size_t oldNVertices = vertices.size();
   size_t oldNTriangles = triangles.size();

   if (vectors.size() == 0) vectors = otherSurface.vectorNames();
   if (scalars.size() == 0) scalars = otherSurface.scalarNames();
   if (pointers.size() == 0) pointers = otherSurface.pointerNames();

   vertices.insert(vertices.end(),otherSurface.getVertices().begin(), otherSurface.getVertices().end());
   triangles.insert(triangles.end(),otherSurface.getTriangles().begin(), otherSurface.getTriangles().end());
   nTriangles = triangles.size();

   vector<CXXTriangle>::iterator triangle;
   vector<CXXTriangle>::iterator triangleEnd = triangles.end();
   for (triangle = (triangles.begin() + oldNTriangles); triangle!=triangleEnd; ++triangle) {
      for (int i=0; i<3; i++) (*triangle)[i] = (*triangle)[i] + oldNVertices;
   }
}

void CXXSurface::compress(double tolerance) {
   vector<CXXSurfaceVertex> compressedVertices;
   compressedVertices.reserve(vertices.size());
   vector<CXXTriangle>compressedTriangles;
   compressedTriangles.reserve(triangles.size());

   size_t vertexHandle = getVectorHandle("vertices");
   size_t normalHandle = getVectorHandle("normals");
   size_t atomHandle = getVectorHandle("atom");

   vector<size_t>equivalences(vertices.size());

   for (int i=0; i<(int)vertices.size(); i++) {
      bool uniqueAndDrawn = true;
      for (int j=0; j<(int)compressedVertices.size() && uniqueAndDrawn; j++) {
         if ((vertices[i].coordRef(vertexHandle).isNearly(compressedVertices[j].coordRef(vertexHandle), tolerance) &&
              vertices[i].coordRef(normalHandle).isNearly(compressedVertices[j].coordRef(normalHandle), tolerance)) &&
             (compressedVertices[j].pointer(atomHandle) == vertices[i].pointer(atomHandle)) ) {
            uniqueAndDrawn = false;
            equivalences[i] = j;
         }
      }
      if (uniqueAndDrawn) {
         compressedVertices.push_back(vertices[i]);
         equivalences[i] = compressedVertices.size()-1;
      }
   }
   size_t equivalent[3];
   for (int i=0; i<(int)triangles.size(); i++) {
      for (int j=0; j<3; j++) equivalent[j] = equivalences[vertex(i,j)];
      if ((equivalent[0] != equivalent[1] && equivalent[0] != equivalent[2] && equivalent[1] != equivalent[2])) {
         triangles[i][0] = equivalent[0];
         triangles[i][1] = equivalent[1];
         triangles[i][2] = equivalent[2];
         compressedTriangles.push_back(triangles[i]);
      }
   }
   vertices=compressedVertices;
   triangles=compressedTriangles;
   nTriangles = triangles.size();
}
