/*
 * MoleculesToTriangles/CXXClasses/CylindersPrimitive.cpp
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "CylindersPrimitive.h"
#include "Renderer.h"
#include "mmdb2/mmdb_manager.h"

void CylindersPrimitive::addHalfAtomBond(mmdb::Atom* atom1, FCXXCoord &atom1Color, mmdb::Atom* atom2, FCXXCoord &atom2Color, float cylinderRadius)
{
    FCXXCoord coord1(atom1->x, atom1->y, atom1->z);
    FCXXCoord coord2(atom2->x, atom2->y, atom2->z);
    
    addHalfAtomBondWithCoords(coord1, atom1, atom1Color, coord2,atom2, atom2Color, cylinderRadius);
}

void CylindersPrimitive::addHalfAtomBondWithCoords(FCXXCoord &coord1, mmdb::Atom* atom1, FCXXCoord &atom1Color, FCXXCoord &coord2,
                                                   mmdb::Atom* atom2, FCXXCoord &atom2Color, float cylinderRadius)
{
    static FCXXCoord xAxis (1., 0., 0., 0.);
    static FCXXCoord yAxis (0., 1., 0., 0.);
    static FCXXCoord zAxis (0., 0., 1., 0.);
    
    FCXXCoord vecAB = coord2-coord1;
    vecAB.normalise();
    float dotX = vecAB*xAxis;
    float dotY = vecAB*yAxis;
    float dotZ = vecAB*zAxis;
    FCXXCoord normal1;
    if (dotX<=dotY && dotX<=dotZ){
        normal1 = vecAB^xAxis;
    }
    else if (dotY<=dotZ && dotY<=dotX){
        normal1 = vecAB^yAxis;
    }
    else {
        normal1 = vecAB^zAxis;
    }
    normal1.normalise();
    FCXXCoord normal2 = vecAB^normal1;
    
    FCXXCoord halfWay = (coord1 + coord2) / 2.;
    
    CylinderPoint cylinderPoint0(coord1, atom1Color, normal1, normal2, cylinderRadius, cylinderRadius, CylinderPoint::CylinderPointTypeStart, atom1);
    addPoint(cylinderPoint0);
    CylinderPoint cylinderPoint1(halfWay, atom1Color, normal1, normal2, cylinderRadius, cylinderRadius, atom1);
    addPoint(cylinderPoint1);
    
    CylinderPoint cylinderPoint2(halfWay, atom2Color, normal1, normal2, cylinderRadius, cylinderRadius, CylinderPoint::CylinderPointTypeStart, atom2);
    addPoint (cylinderPoint2);
    CylinderPoint cylinderPoint3(coord2, atom2Color, normal1, normal2, cylinderRadius, cylinderRadius, atom2);
    addPoint (cylinderPoint3);
}

void CylindersPrimitive::generateArrays()
{
    //std::cout << "In cylinder generateArrays "<<points.size()<<" "<<angularSampling;
    float angularStep = (2.*M_PI) / (float)angularSampling;

    // Count how many separate tube segments we have (separated by CylinderPointTypeStart)
    int nCaps = 0;
    if (points.size() >= 2) {
        // Count start caps (each segment needs a start cap)
        nCaps++; // First segment start
        for (size_t i = 1; i < points.size(); i++) {
            if (points[i].type == CylinderPoint::CylinderPointTypeStart) {
                nCaps += 2; // End cap for previous segment, start cap for new segment
            }
        }
        nCaps++; // Last segment end
    }

    // Cap geometry: 3 intermediate rings + 1 apex per cap
    // Extra vertices per cap: 3*angularSampling (intermediate rings) + 1 (apex)
    // Extra triangles per cap: 2*angularSampling * 3 (between rings) + angularSampling (last ring to apex)
    int nCapRings = 3;
    int capVertices = nCaps * (nCapRings * angularSampling + 1);
    int capTriangles = nCaps * (2 * nCapRings * angularSampling + angularSampling);

    int totalVertices = angularSampling * points.size() + capVertices;
    unsigned long totalIndices = 6 * angularSampling * points.size() + 3 * capTriangles;

    vertexColorNormalArray = new VertexColorNormal[totalVertices];
    atomArray = new const mmdb::Atom*[totalVertices];
    indexArray = new GLIndexType[totalIndices];

    int iGLVertex = 0;
    for (int i=0; i< points.size(); i++){
        for (int j=0; j<angularSampling; j++){
            //First work out coords and vertices and copy these into relevant arrays
            float angle = (float)j * angularStep;
            FCXXCoord vertex = points[i].vertex + points[i].normalOne*(points[i].radiusOne*sinf(angle)) + points[i].normalTwo*(points[i].radiusTwo*cosf(angle));
            FCXXCoord normal = points[i].normalOne*sinf(angle)/pow(points[i].radiusOne,0.5) +
            points[i].normalTwo*cosf(angle)/pow(points[i].radiusTwo,0.5);
            normal.normalise();
            for (int k=0; k<4; k++){
                vertexColorNormalArray[iGLVertex].vertex[k] = vertex[k];

                float floatComponentValue = points[i].color[k] * 255.;
                int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255.? 255 : floatComponentValue));

                vertexColorNormalArray[iGLVertex].color[k] = uintComponentValue;
                vertexColorNormalArray[iGLVertex].normal[k] = normal[k];
            }
            atomArray[iGLVertex] = points[i].atom;
            iGLVertex++;
        }
    }

    int tubeVertexCount = iGLVertex;
    _nVertices = iGLVertex;

    _nTriangles = 0;
    int iIndex = 0;

    int highestIndex = -1;
    if (points.size() > 0) {
       for (int i=0; i< points.size()-1; i++){
	  if (points[i+1].type != CylinderPoint::CylinderPointTypeStart){
	     for (int j=0; j<angularSampling; j++){

                int i0, i1, i2, i3, i4, i5;
                //Now work out indices
                indexArray[iIndex++] = i0 = (j%angularSampling) + (angularSampling * i);
                indexArray[iIndex++] = i1 = ((j+1)%angularSampling) + (angularSampling * i);
                indexArray[iIndex++] = i2 = (j%angularSampling) + (angularSampling * (i+1));
                _nTriangles++;

                indexArray[iIndex++] = i3 = (j%angularSampling) + (angularSampling * (i+1));
                indexArray[iIndex++] = i4 = ((j+1)%angularSampling) + (angularSampling * i) ;
                indexArray[iIndex++] = i5 = ((j+1)%angularSampling) + (angularSampling * (i+1));
                _nTriangles++;

                highestIndex = (highestIndex > i0?highestIndex:i0);
                highestIndex = (highestIndex > i1?highestIndex:i1);
                highestIndex = (highestIndex > i2?highestIndex:i2);
                highestIndex = (highestIndex > i3?highestIndex:i3);
                highestIndex = (highestIndex > i4?highestIndex:i4);
                highestIndex = (highestIndex > i5?highestIndex:i5);
	     }
	  }
       }
    }

    // Generate dome end caps
    if (points.size() >= 2) {
        // Lambda to add a dome cap
        // tubeRingIndex: index of the tube ring (in points array) to cap
        // capDirection: unit vector pointing outward from tube end
        // isStartCap: true if this is the start of a segment (affects winding order)
        auto addDomeCap = [&](int tubeRingIndex, FCXXCoord capDirection, bool isStartCap) {
            const CylinderPoint& pt = points[tubeRingIndex];
            int tubeRingVertexBase = tubeRingIndex * angularSampling;

            // Cap parameters - 3 intermediate rings for smooth curvature
            float maxPhi = M_PI / 6.0f;  // 30 degrees total - shallow dome
            float capDepth = (pt.radiusOne + pt.radiusTwo) * 0.5f * 0.6f;

            // Apex
            FCXXCoord apex = pt.vertex + capDirection * capDepth;
            int apexIndex = iGLVertex;

            // Add apex vertex
            for (int k = 0; k < 4; k++) {
                vertexColorNormalArray[iGLVertex].vertex[k] = apex[k];
                float floatComponentValue = pt.color[k] * 255.;
                int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255. ? 255 : floatComponentValue));
                vertexColorNormalArray[iGLVertex].color[k] = uintComponentValue;
                vertexColorNormalArray[iGLVertex].normal[k] = capDirection[k];
            }
            atomArray[iGLVertex] = pt.atom;
            iGLVertex++;

            // Add 3 intermediate ring vertices
            int ringBases[3];
            for (int ring = 0; ring < 3; ring++) {
                // Latitude angles: distribute evenly from near-equator to near-apex
                // ring 0: closest to tube edge, ring 2: closest to apex
                float phi = maxPhi * (float)(ring + 1) / 4.0f;  // 7.5°, 15°, 22.5°
                float cosPhi = cosf(phi);
                float sinPhi = sinf(phi);

                FCXXCoord ringCenter = pt.vertex + capDirection * (capDepth * sinPhi);
                float ringRadiusOne = pt.radiusOne * cosPhi;
                float ringRadiusTwo = pt.radiusTwo * cosPhi;

                ringBases[ring] = iGLVertex;
                for (int j = 0; j < angularSampling; j++) {
                    float angle = (float)j * angularStep;
                    FCXXCoord vertex = ringCenter + pt.normalOne * (ringRadiusOne * sinf(angle)) + pt.normalTwo * (ringRadiusTwo * cosf(angle));

                    // Normal blends between radial and axial
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
                    atomArray[iGLVertex] = pt.atom;
                    iGLVertex++;
                }
            }

            // Triangles: tube ring to first intermediate ring
            for (int j = 0; j < angularSampling; j++) {
                int j1 = (j + 1) % angularSampling;
                int outerA = tubeRingVertexBase + j;
                int outerB = tubeRingVertexBase + j1;
                int innerA = ringBases[0] + j;
                int innerB = ringBases[0] + j1;

                if (isStartCap) {
                    indexArray[iIndex++] = outerA;
                    indexArray[iIndex++] = innerA;
                    indexArray[iIndex++] = outerB;
                    _nTriangles++;
                    indexArray[iIndex++] = outerB;
                    indexArray[iIndex++] = innerA;
                    indexArray[iIndex++] = innerB;
                    _nTriangles++;
                } else {
                    indexArray[iIndex++] = outerA;
                    indexArray[iIndex++] = outerB;
                    indexArray[iIndex++] = innerA;
                    _nTriangles++;
                    indexArray[iIndex++] = outerB;
                    indexArray[iIndex++] = innerB;
                    indexArray[iIndex++] = innerA;
                    _nTriangles++;
                }
            }

            // Triangles: between intermediate rings
            for (int ring = 0; ring < 2; ring++) {
                for (int j = 0; j < angularSampling; j++) {
                    int j1 = (j + 1) % angularSampling;
                    int outerA = ringBases[ring] + j;
                    int outerB = ringBases[ring] + j1;
                    int innerA = ringBases[ring + 1] + j;
                    int innerB = ringBases[ring + 1] + j1;

                    if (isStartCap) {
                        indexArray[iIndex++] = outerA;
                        indexArray[iIndex++] = innerA;
                        indexArray[iIndex++] = outerB;
                        _nTriangles++;
                        indexArray[iIndex++] = outerB;
                        indexArray[iIndex++] = innerA;
                        indexArray[iIndex++] = innerB;
                        _nTriangles++;
                    } else {
                        indexArray[iIndex++] = outerA;
                        indexArray[iIndex++] = outerB;
                        indexArray[iIndex++] = innerA;
                        _nTriangles++;
                        indexArray[iIndex++] = outerB;
                        indexArray[iIndex++] = innerB;
                        indexArray[iIndex++] = innerA;
                        _nTriangles++;
                    }
                }
            }

            // Triangles: innermost ring to apex (fan)
            for (int j = 0; j < angularSampling; j++) {
                int j1 = (j + 1) % angularSampling;
                int innerA = ringBases[2] + j;
                int innerB = ringBases[2] + j1;

                if (isStartCap) {
                    indexArray[iIndex++] = innerA;
                    indexArray[iIndex++] = apexIndex;
                    indexArray[iIndex++] = innerB;
                } else {
                    indexArray[iIndex++] = innerA;
                    indexArray[iIndex++] = innerB;
                    indexArray[iIndex++] = apexIndex;
                }
                _nTriangles++;
            }
        };

        // Find segment boundaries and add caps
        int segmentStart = 0;
        for (int i = 0; i <= (int)points.size(); i++) {
            bool isSegmentEnd = (i == (int)points.size()) ||
                               (i > 0 && points[i].type == CylinderPoint::CylinderPointTypeStart);

            if (isSegmentEnd && i > segmentStart) {
                int segmentEnd = i - 1;

                // Only add caps if segment has at least 2 points
                if (segmentEnd > segmentStart) {
                    // Start cap
                    FCXXCoord startDir = points[segmentStart].vertex - points[segmentStart + 1].vertex;
                    startDir.normalise();
                    addDomeCap(segmentStart, startDir, true);

                    // End cap
                    FCXXCoord endDir = points[segmentEnd].vertex - points[segmentEnd - 1].vertex;
                    endDir.normalise();
                    addDomeCap(segmentEnd, endDir, false);
                }

                segmentStart = i;
            }
        }
    }

    _nVertices = iGLVertex;
    //std::cout << "Highest index is " <<  highestIndex << std::endl;
}

CylindersPrimitive::CylindersPrimitive() : VertexColorNormalPrimitive () {
    vertexColorNormalArray = 0;
    indexArray = 0;
    emptyArrays();
    drawModeGL = DrawAsTriangleStrip;
    enableColorGL = true;
    primitiveType = DisplayPrimitive::PrimitiveType::CylinderPrimitive;
}


void CylindersPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
    //std::cout << "in Cylinder render\n";
    if (vertexColorNormalArray == 0){
        generateArrays();
    }
    renderer->renderVertexColorNormalPrimitive(this);
}

