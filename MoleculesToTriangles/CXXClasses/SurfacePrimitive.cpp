/*
 *  SurfacePrimitive.cpp
 *  AesopCoreData
 *
 *  Created by Martin Noble on 23/03/2010.
 *  Copyright 2010 LMB, Oxford University. All rights reserved.
 *
 */

#include "SurfacePrimitive.h"
#include "Renderer.h"
#include "MoleculesToTriangles/CXXSurface/CXXSurface.h"
#include "MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.h"

SurfacePrimitive::SurfacePrimitive() : VertexColorNormalPrimitive (){
    drawModeGL = DrawAsTriangles;
    enableColorGL = true;
    cxxSurfaceMaker = 0;
    primitiveType = DisplayPrimitive::PrimitiveType::SurfacePrimitive;
}

void SurfacePrimitive::generateArrays()
{
    CXXSurfaceMaker &containerSurface(*cxxSurfaceMaker);

    double coords[4];

    _nVertices = 0;
    _nTriangles = 0;
    for (vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
         surfIter != containerSurface.getChildSurfaces().end();
         ++surfIter){
        _nVertices += surfIter->numberOfVertices();
        _nTriangles += surfIter->numberOfTriangles();
    }
    
    vertexColorNormalArray = new VertexColorNormal[_nVertices];
    atomArray = new const mmdb::Atom*[_nVertices];
    
    size_t iGlobalVert=0;
    for (vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
         surfIter != containerSurface.getChildSurfaces().end();
         ++surfIter){
        CXXSurface &childSurface = *surfIter;
        for (size_t iLocalVert=0; iLocalVert< childSurface.numberOfVertices(); iLocalVert++){
            VertexColorNormal &vcn(vertexColorNormalArray[iGlobalVert]);
            //Copy vertex into vertices array
            if (!childSurface.getCoord("vertices", iLocalVert, coords)){
                for (size_t k=0; k<3; k++) vcn.vertex[k] = coords[k];
            }
            vcn.vertex[3] = 1.;
            if (!childSurface.getCoord("normals", iLocalVert, coords)){
                for (size_t k=0; k<3; k++) vcn.normal[k] = coords[k];
            }
            vcn.normal[3] = 1.;
            
            
            //Copy color into vertices array
            if (!childSurface.getCoord("colour", iLocalVert, coords)){
                for (int k=0; k<4; k++) {
                    
                    float floatComponentValue = coords[k]*255.;
                    int uintComponentValue = (floatComponentValue < 0. ? 0 : (floatComponentValue > 255.? 255 : floatComponentValue));
                    
                    
                    vcn.color[k] = uintComponentValue;
                }
            }
            else for (int k=0; k<4; k++) vcn.color[k] = 0.5;
            
            //Copy atom pointer into atom pointers array
            mmdb::Atom *atomPointer;
            int result = childSurface.getPointer("atom", iLocalVert, (void**)&(atomPointer));
            if (result) atomArray[iGlobalVert] = atomPointer;
            iGlobalVert++;
        }
    }
    indexArray = new GLIndexType[3*_nTriangles];
    size_t iOffset = 0;
    int idx=0;
    for (vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
         surfIter != containerSurface.getChildSurfaces().end();
         ++surfIter){
        CXXSurface &childSurface = *surfIter;
        for (int i=0; i< childSurface.numberOfTriangles(); i++){
            for (int j=0; j<3; j++){
                indexArray[idx++] = GLIndexType(childSurface.vertex(i,j) + iOffset);
            }
        }
        iOffset += childSurface.numberOfVertices();
    }
    
    //Here an idea to liberate memory 
    delete cxxSurfaceMaker;
    cxxSurfaceMaker = 0;
}

SurfacePrimitive::SurfacePrimitive(mmdb::Manager *mmdb, int chunkHndl, int selHnd, std::shared_ptr<ColorScheme> _colorScheme, enum SurfaceType type, float probeRadius, float radiusMultiplier) {
    primitiveType = DisplayPrimitive::PrimitiveType::SurfacePrimitive;
    cxxSurfaceMaker = 0;
    vertexColorNormalArray = 0;
    indexArray = 0;
    colorScheme = _colorScheme;
    cxxSurfaceMaker = new CXXSurfaceMaker();
    try {
        if (type == AccessibleSurface){
            cxxSurfaceMaker->calculateAccessibleFromAtoms(mmdb, chunkHndl, selHnd, probeRadius, 30.*(M_PI/180.), radiusMultiplier, false);
        }
        else if (type == VdWSurface){
            cxxSurfaceMaker->calculateVDWFromAtoms(mmdb, chunkHndl, selHnd, probeRadius, 30.*(M_PI/180.), radiusMultiplier, false);
        }
        else if (type == MolecularSurface){
            cxxSurfaceMaker->calculateFromAtoms(mmdb, chunkHndl, selHnd, probeRadius, 30.*(M_PI/180.), radiusMultiplier, false);
        }
    }
    catch (exception &e) {
        std::cout << "exception caught" << std::endl;
        delete cxxSurfaceMaker;
        cxxSurfaceMaker = 0;
    }
    
    if (cxxSurfaceMaker){
        std::map<std::shared_ptr<ColorRule>,int> handles = colorScheme->prepareForMMDB(mmdb);
        //Assign colors based on ColorScheme
        CXXSurfaceMaker &mySurfaceMaker = *cxxSurfaceMaker;
        for (vector<CXXSurface>::iterator surfIter = mySurfaceMaker.getChildSurfaces().begin();
             surfIter != mySurfaceMaker.getChildSurfaces().end();
             ++surfIter){
            CXXSurface &mySurface = *surfIter;
            for (int iVertex = 0; iVertex<mySurface.numberOfVertices(); iVertex++){
                mmdb::Atom* theAtom;
                int result = mySurface.getPointer("atom", iVertex, (void **) &theAtom);
                if (result == 0){
                    FCXXCoord color = colorScheme->colorForAtom(theAtom, handles);
                    mySurface.setCoord("colour", iVertex, CXXCoord<CXXCoord_ftype>(color[0],color[1],color[2],color[3]));
                }
                else {
                    std::cout << "Anable to assign atom to scheme" << theAtom;
                }
            }
        }
        for (vector<CXXSurface>::iterator surfIter = mySurfaceMaker.getChildSurfaces().begin();
             surfIter != mySurfaceMaker.getChildSurfaces().end();
             ++surfIter){
            CXXSurface &mySurface = *surfIter;
            for (int iVertex = 0; iVertex<mySurface.numberOfVertices(); iVertex++){
                mmdb::Atom* theAtom;
                int result = mySurface.getPointer("atom", iVertex, (void **) &theAtom);
                if (result == 0){
                    FCXXCoord color = colorScheme->colorForAtom(theAtom, handles);
                    mySurface.setCoord("colour", iVertex, CXXCoord<CXXCoord_ftype>(color[0],color[1],color[2],color[3]));
                }
                else {
                    std::cout << "Anable to assign atom to scheme" << theAtom;
                }
            }
        }
        std::cout << cxxSurfaceMaker->report();
    }
};


