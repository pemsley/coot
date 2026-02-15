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
    vertexColorNormalArray = new VertexColorNormal[angularSampling*points.size()];
    atomArray = new const mmdb::Atom*[angularSampling*points.size()];
    unsigned long nIndices = 6*angularSampling*points.size();
    indexArray = new GLIndexType[nIndices];

    int iGLVertex = 0;
    for (int i=0; i< points.size(); i++){
        for (int j=0; j<angularSampling; j++){
            //First work out coords and vertices and copy these into relevant arrays
            float angle = (float)j * angularStep;
            FCXXCoord vertex = points[i].vertex + points[i].normalOne*(points[i].radiusOne*sinf(angle)) + points[i].normalTwo*(points[i].radiusTwo*cosf(angle));
            FCXXCoord normal = points[i].normalOne*sinf(angle)/points[i].radiusOne +
            points[i].normalTwo*cosf(angle)/points[i].radiusTwo;
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

