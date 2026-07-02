/*
 * MoleculesToTriangles/CXXClasses/CylindersPrimitive-gemmi.cc
 *
 * gemmi-native twin of CylindersPrimitive.cpp - identical tube/dome geometry,
 * fed from coordinates (no mmdb::Atom*; the picking atomArray is not populated).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <cmath>
#include "CylindersPrimitive-gemmi.hh"
#include "Renderer.h"

void coot::m2t::CylindersPrimitive::addHalfAtomBondWithCoords(FCXXCoord &coord1, FCXXCoord &atom1Color,
                                                              FCXXCoord &coord2, FCXXCoord &atom2Color,
                                                              float cylinderRadius)
{
   static FCXXCoord xAxis(1., 0., 0., 0.);
   static FCXXCoord yAxis(0., 1., 0., 0.);
   static FCXXCoord zAxis(0., 0., 1., 0.);

   FCXXCoord vecAB = coord2 - coord1;
   vecAB.normalise();
   float dotX = vecAB*xAxis;
   float dotY = vecAB*yAxis;
   float dotZ = vecAB*zAxis;
   FCXXCoord normal1;
   if (dotX <= dotY && dotX <= dotZ)      normal1 = vecAB^xAxis;
   else if (dotY <= dotZ && dotY <= dotX) normal1 = vecAB^yAxis;
   else                                   normal1 = vecAB^zAxis;
   normal1.normalise();
   FCXXCoord normal2 = vecAB^normal1;

   FCXXCoord halfWay = (coord1 + coord2) / 2.;

   CylinderPoint cylinderPoint0(coord1, atom1Color, normal1, normal2, cylinderRadius, cylinderRadius, CylinderPoint::CylinderPointTypeStart);
   addPoint(cylinderPoint0);
   CylinderPoint cylinderPoint1(halfWay, atom1Color, normal1, normal2, cylinderRadius, cylinderRadius);
   addPoint(cylinderPoint1);
   CylinderPoint cylinderPoint2(halfWay, atom2Color, normal1, normal2, cylinderRadius, cylinderRadius, CylinderPoint::CylinderPointTypeStart);
   addPoint(cylinderPoint2);
   CylinderPoint cylinderPoint3(coord2, atom2Color, normal1, normal2, cylinderRadius, cylinderRadius);
   addPoint(cylinderPoint3);
}

void coot::m2t::CylindersPrimitive::generateArrays()
{
   float angularStep = (2.*M_PI) / (float)angularSampling;

   int nCaps = 0;
   if (points.size() >= 2) {
      nCaps++; // First segment start
      for (size_t i = 1; i < points.size(); i++) {
         if (points[i].type == CylinderPoint::CylinderPointTypeStart)
            nCaps += 2;
      }
      nCaps++; // Last segment end
   }

   int nCapRings = 3;
   int capVertices = nCaps * (nCapRings * angularSampling + 1);
   int capTriangles = nCaps * (2 * nCapRings * angularSampling + angularSampling);

   int totalVertices = angularSampling * points.size() + capVertices;
   unsigned long totalIndices = 6 * angularSampling * points.size() + 3 * capTriangles;

   vertexColorNormalArray = new VertexColorNormal[totalVertices];
   indexArray = new GLIndexType[totalIndices];

   int iGLVertex = 0;
   for (int i = 0; i < (int)points.size(); i++) {
      for (int j = 0; j < angularSampling; j++) {
         float angle = (float)j * angularStep;
         FCXXCoord vertex = points[i].vertex + points[i].normalOne*(points[i].radiusOne*sinf(angle)) + points[i].normalTwo*(points[i].radiusTwo*cosf(angle));
         FCXXCoord normal = points[i].normalOne*sinf(angle)/points[i].radiusOne +
                            points[i].normalTwo*cosf(angle)/points[i].radiusTwo;
         normal.normalise();
         for (int k = 0; k < 4; k++) {
            vertexColorNormalArray[iGLVertex].vertex[k] = vertex[k];
            float floatComponentValue = points[i].color[k] * 255.;
            int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255.? 255 : floatComponentValue));
            vertexColorNormalArray[iGLVertex].color[k] = uintComponentValue;
            vertexColorNormalArray[iGLVertex].normal[k] = normal[k];
         }
         iGLVertex++;
      }
   }

   _nVertices = iGLVertex;
   _nTriangles = 0;
   int iIndex = 0;

   if (points.size() > 0) {
      for (int i = 0; i < (int)points.size()-1; i++) {
         if (points[i+1].type != CylinderPoint::CylinderPointTypeStart) {
            for (int j = 0; j < angularSampling; j++) {
               indexArray[iIndex++] = (j%angularSampling) + (angularSampling * i);
               indexArray[iIndex++] = ((j+1)%angularSampling) + (angularSampling * i);
               indexArray[iIndex++] = (j%angularSampling) + (angularSampling * (i+1));
               _nTriangles++;
               indexArray[iIndex++] = (j%angularSampling) + (angularSampling * (i+1));
               indexArray[iIndex++] = ((j+1)%angularSampling) + (angularSampling * i);
               indexArray[iIndex++] = ((j+1)%angularSampling) + (angularSampling * (i+1));
               _nTriangles++;
            }
         }
      }
   }

   if (points.size() >= 2) {
      auto addDomeCap = [&](int tubeRingIndex, FCXXCoord capDirection, bool isStartCap) {
         const CylinderPoint& pt = points[tubeRingIndex];
         int tubeRingVertexBase = tubeRingIndex * angularSampling;

         float maxPhi = M_PI / 6.0f;
         float capDepth = (pt.radiusOne + pt.radiusTwo) * 0.5f * 0.6f;

         FCXXCoord apex = pt.vertex + capDirection * capDepth;
         int apexIndex = iGLVertex;
         for (int k = 0; k < 4; k++) {
            vertexColorNormalArray[iGLVertex].vertex[k] = apex[k];
            float floatComponentValue = pt.color[k] * 255.;
            int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255. ? 255 : floatComponentValue));
            vertexColorNormalArray[iGLVertex].color[k] = uintComponentValue;
            vertexColorNormalArray[iGLVertex].normal[k] = capDirection[k];
         }
         iGLVertex++;

         int ringBases[3];
         for (int ring = 0; ring < 3; ring++) {
            float phi = maxPhi * (float)(ring + 1) / 4.0f;
            float cosPhi = cosf(phi);
            float sinPhi = sinf(phi);
            FCXXCoord ringCenter = pt.vertex + capDirection * (capDepth * sinPhi);
            float ringRadiusOne = pt.radiusOne * cosPhi;
            float ringRadiusTwo = pt.radiusTwo * cosPhi;
            ringBases[ring] = iGLVertex;
            for (int j = 0; j < angularSampling; j++) {
               float angle = (float)j * angularStep;
               FCXXCoord vertex = ringCenter + pt.normalOne * (ringRadiusOne * sinf(angle)) + pt.normalTwo * (ringRadiusTwo * cosf(angle));
               FCXXCoord radialNormal = pt.normalOne * sinf(angle) / powf(pt.radiusOne, 0.5f) +
                                        pt.normalTwo * cosf(angle) / powf(pt.radiusTwo, 0.5f);
               radialNormal.normalise();
               FCXXCoord normal = radialNormal * cosPhi + capDirection * sinPhi;
               normal.normalise();
               for (int k = 0; k < 4; k++) {
                  vertexColorNormalArray[iGLVertex].vertex[k] = vertex[k];
                  float floatComponentValue = pt.color[k] * 255.;
                  int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255. ? 255 : floatComponentValue));
                  vertexColorNormalArray[iGLVertex].color[k] = uintComponentValue;
                  vertexColorNormalArray[iGLVertex].normal[k] = normal[k];
               }
               iGLVertex++;
            }
         }

         for (int j = 0; j < angularSampling; j++) {
            int j1 = (j + 1) % angularSampling;
            int outerA = tubeRingVertexBase + j;
            int outerB = tubeRingVertexBase + j1;
            int innerA = ringBases[0] + j;
            int innerB = ringBases[0] + j1;
            if (isStartCap) {
               indexArray[iIndex++] = outerA; indexArray[iIndex++] = innerA; indexArray[iIndex++] = outerB; _nTriangles++;
               indexArray[iIndex++] = outerB; indexArray[iIndex++] = innerA; indexArray[iIndex++] = innerB; _nTriangles++;
            } else {
               indexArray[iIndex++] = outerA; indexArray[iIndex++] = outerB; indexArray[iIndex++] = innerA; _nTriangles++;
               indexArray[iIndex++] = outerB; indexArray[iIndex++] = innerB; indexArray[iIndex++] = innerA; _nTriangles++;
            }
         }
         for (int ring = 0; ring < 2; ring++) {
            for (int j = 0; j < angularSampling; j++) {
               int j1 = (j + 1) % angularSampling;
               int outerA = ringBases[ring] + j;
               int outerB = ringBases[ring] + j1;
               int innerA = ringBases[ring + 1] + j;
               int innerB = ringBases[ring + 1] + j1;
               if (isStartCap) {
                  indexArray[iIndex++] = outerA; indexArray[iIndex++] = innerA; indexArray[iIndex++] = outerB; _nTriangles++;
                  indexArray[iIndex++] = outerB; indexArray[iIndex++] = innerA; indexArray[iIndex++] = innerB; _nTriangles++;
               } else {
                  indexArray[iIndex++] = outerA; indexArray[iIndex++] = outerB; indexArray[iIndex++] = innerA; _nTriangles++;
                  indexArray[iIndex++] = outerB; indexArray[iIndex++] = innerB; indexArray[iIndex++] = innerA; _nTriangles++;
               }
            }
         }
         for (int j = 0; j < angularSampling; j++) {
            int j1 = (j + 1) % angularSampling;
            int innerA = ringBases[2] + j;
            int innerB = ringBases[2] + j1;
            if (isStartCap) {
               indexArray[iIndex++] = innerA; indexArray[iIndex++] = apexIndex; indexArray[iIndex++] = innerB;
            } else {
               indexArray[iIndex++] = innerA; indexArray[iIndex++] = innerB; indexArray[iIndex++] = apexIndex;
            }
            _nTriangles++;
         }
      };

      int segmentStart = 0;
      for (int i = 0; i <= (int)points.size(); i++) {
         bool isSegmentEnd = (i == (int)points.size()) ||
                            (i > 0 && points[i].type == CylinderPoint::CylinderPointTypeStart);
         if (isSegmentEnd && i > segmentStart) {
            int segmentEnd = i - 1;
            if (segmentEnd > segmentStart) {
               FCXXCoord startDir = points[segmentStart].vertex - points[segmentStart + 1].vertex;
               startDir.normalise();
               addDomeCap(segmentStart, startDir, true);
               FCXXCoord endDir = points[segmentEnd].vertex - points[segmentEnd - 1].vertex;
               endDir.normalise();
               addDomeCap(segmentEnd, endDir, false);
            }
            segmentStart = i;
         }
      }
   }

   _nVertices = iGLVertex;
}

coot::m2t::CylindersPrimitive::CylindersPrimitive() : VertexColorNormalPrimitive() {
   angularSampling = 12;
   vertexColorNormalArray = 0;
   indexArray = 0;
   emptyArrays();
   drawModeGL = DrawAsTriangleStrip;
   enableColorGL = true;
   primitiveType = DisplayPrimitive::PrimitiveType::CylinderPrimitive;
}

void coot::m2t::CylindersPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
   if (vertexColorNormalArray == 0)
      generateArrays();
   renderer->renderVertexColorNormalPrimitive(this);
}
