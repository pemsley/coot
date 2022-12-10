//
//  BoxSectionPrimitive.cpp
//  AesopCD_macos
//
//  Created by Martin Noble on 11/08/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#include <iostream>

// #include <epoxy/gl.h>

#include "BoxSectionPrimitive.h"
#include "Renderer.h"

void BoxSectionPrimitive::generateArrays()
{   
    unsigned long nPoints = 8*points.size();
    vertexColorNormalArray = new VertexColorNormal[nPoints];
    atomArray = new const mmdb::Atom*[nPoints];
    unsigned long nIndices = 24*points.size();
    indexArray = new GLIndexType[nIndices];    
 
    float vertexCoefficients[8][2] = {
        -1,1, 
        1,1, 
        
        1,1, 
        1,-1,  
        
        1,-1,
        -1,-1,
        
        -1,-1,
        -1,1
    };
    
    float normalCoefficients[8][2] = {
        0,1,
        0,1,
        
        1,0,
        1,0,
        
        0,-1,
        0,-1,
        
        -1,0,
        -1,0
    };

    int iGLVertex = 0;
    for (int i=0; i< points.size(); i++){
        for (int j=0; j<8; j++){
            //First work out coords and vertices and copy these into relevant arrays
            FCXXCoord vertex = points[i].vertex + points[i].normalOne*(points[i].radiusOne*vertexCoefficients[j][0]) + 
                                                 points[i].normalTwo*(points[i].radiusTwo*vertexCoefficients[j][1]);
            FCXXCoord normal = points[i].normalOne*normalCoefficients[j][0] + points[i].normalTwo*normalCoefficients[j][1];
            for (int k=0; k<4; k++){
                vertexColorNormalArray[iGLVertex].vertex[k] = vertex[k];

                float floatComponentValue = points[i].color[k] * 255.;
                int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255.? 255 : floatComponentValue));

                if (boxEdgeColor.first) {
                   if (! (j&2)) {
                      vertexColorNormalArray[iGLVertex].color[k] = uintComponentValue;
                   } else {
                      // should use the same protection as above
                      vertexColorNormalArray[iGLVertex].color[k] = 255.0 * boxEdgeColor.second[k];
                   }
                } else {
                   vertexColorNormalArray[iGLVertex].color[k] = 255.0 * uintComponentValue;
                }
                vertexColorNormalArray[iGLVertex].normal[k] = normal[k];
            }
            atomArray[iGLVertex] = points[i].atom;
            iGLVertex++;
        }
    }
    _nVertices = iGLVertex;
    
    _nTriangles = 0;
    int iIndex = 0;
    
    for (int i=0; i< points.size()-1; i++){
        int majorOffset = 8*i;
        int nextMajorOffset = 8*(i+1);
        if (points[i+1].type != CylinderPoint::CylinderPointTypeStart){
            for (int j=0; j<4; j++){
                int minorOffset = 2*j;
                indexArray[iIndex++] = (0+minorOffset)%8+majorOffset;
                indexArray[iIndex++] = (1+minorOffset)%8+majorOffset;
                indexArray[iIndex++] = (0+minorOffset)%8+nextMajorOffset;
                _nTriangles++;

                indexArray[iIndex++] = (1+minorOffset)%8+majorOffset;
                indexArray[iIndex++] = (1+minorOffset)%8+nextMajorOffset;
                indexArray[iIndex++] = (0+minorOffset)%8+nextMajorOffset;
                _nTriangles++;
            }
        }
    }
}


void BoxSectionPrimitive::renderWithRenderer(std::shared_ptr<Renderer> renderer)
{
    if (vertexColorNormalArray == 0){
        generateArrays();
    }
    renderer->renderVertexColorNormalPrimitive(this);
}

BoxSectionPrimitive::BoxSectionPrimitive() : CylindersPrimitive () {
    primitiveType = DisplayPrimitive::PrimitiveType::BoxSectionPrimitive;
    //std::cout << "Creating boxSectionPrimitive" << std::endl;
}
