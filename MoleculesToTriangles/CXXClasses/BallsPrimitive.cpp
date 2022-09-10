//
//  BallsPrimitive.cpp
//  AesopCD_ios
//
//  Created by Martin Noble on 09/02/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#include "BallsPrimitive.h"
#include "Renderer.h"
#include "oglPolyhedron.h"
#include "CXXMatrix.h"

void BallsPrimitive::generateArrays()
{
    oglPolyhedron *poly = oglPolyhedron::dodecaHedron4();
    
    unsigned long requestedVertices = balls.size() * ((Polyhedron *)poly)->nVertices();
    unsigned long requestedIndices = 3 * balls.size() * ((Polyhedron *)poly)->nFaces();
    vertexColorNormalArray = new VertexColorNormal[requestedVertices];
    indexArray = new GLIndexType[requestedIndices];

    std::vector<Ball>::iterator ball = balls.begin();
    std::vector<Ball>::iterator ballEnd = balls.end();
    _nVertices = 0;
    _nTriangles = 0;
    int iIndex = 0;
    unsigned long highestIndex = 0;
    
    static FCXXCoord xAxis(1.,0.,0.,0.);
    static FCXXCoord yAxis(0.,1.,0.,0.);
    static FCXXCoord zAxis(0.,0.,1.,0.);
    
    
    for (; ball != ballEnd; ++ball){
        float dotX = ball->normal*xAxis;
        float dotY = ball->normal*yAxis;
        float dotZ = ball->normal*zAxis;
        FCXXCoord normal1;
        if (dotX<=dotY && dotX<=dotZ){
            normal1 = ball->normal^xAxis;
        }
        else if (dotY<=dotZ && dotY<=dotX){
            normal1 = ball->normal^yAxis;
        }
        else {
            normal1 = ball->normal^zAxis;
        }
        normal1.normalise();
        FCXXCoord normal2 = ball->normal^normal1;
        
        FCXXCoord scaledNormal1 = normal1*ball->radius;
        FCXXCoord scaledNormal2 = normal2*ball->radius;
        FCXXCoord scaledNormal = ball->normal * ball->radiusAlongNormal;
        
        float rotationAndScaleAsFloats[] = {
            scaledNormal1.x(), scaledNormal1.y(), scaledNormal1.z(), 0.,
            scaledNormal2.x(), scaledNormal2.y(), scaledNormal2.z(), 0.,
            scaledNormal.x(), scaledNormal.y(), scaledNormal.z(), 0.,
            0.,0.,0.,1.};
        CXXMatrix rotationAndScale = CXXMatrix(rotationAndScaleAsFloats);

        FCXXCoord inverseScaledNormal1 = normal1/ball->radius;
        FCXXCoord inverseScaledNormal2 = normal2/ball->radius;
        FCXXCoord inverseScaledNormal = ball->normal / ball->radiusAlongNormal;
        float rotationAndInverseScaleAsFloats[] = {
            inverseScaledNormal1.x(), inverseScaledNormal1.y(), inverseScaledNormal1.z(), 0.,
            inverseScaledNormal2.x(), inverseScaledNormal2.y(), inverseScaledNormal2.z(), 0.,
            inverseScaledNormal.x(), inverseScaledNormal.y(), inverseScaledNormal.z(), 0.,
            0.,0.,0.,1.};
        CXXMatrix rotationAndInverseScale = CXXMatrix(rotationAndInverseScaleAsFloats);

        
        unsigned long vertexOffset = _nVertices;
        
        for (int iVertex = 0; iVertex < ((Polyhedron *)poly)->nVertices(); iVertex++){
            FCXXCoord vertex = ball->centre + (rotationAndScale * poly->vertex(iVertex));
            FCXXCoord normal = rotationAndInverseScale * poly->vertex(iVertex);
            normal.normalise();

            for (int k=0; k<4; k++){
                vertexColorNormalArray[_nVertices].vertex[k] = vertex[k];
                float floatComponentValue = ball->color[k] * 255.;
                int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255.? 255 : floatComponentValue));
                vertexColorNormalArray[_nVertices].color[k] = uintComponentValue;
                vertexColorNormalArray[_nVertices].normal[k] = normal[k];
            }
            _nVertices++;
        }
        for (int iFace = 0; iFace < ((Polyhedron *)poly)->nFaces(); iFace++){
            PolyhedronFace &face = ((Polyhedron *)poly)->face(iFace);
            for (int i=0; i<3; i++){
                unsigned long iVertex = face[i] + vertexOffset;
                highestIndex = (highestIndex > iVertex ? highestIndex : iVertex);
                indexArray[iIndex++] = GLIndexType(iVertex);
            }
            _nTriangles++;
        }
    }
    //std::cout << "Highest balls index is " <<  highestIndex << std::endl;
}

