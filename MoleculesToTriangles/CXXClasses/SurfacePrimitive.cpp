/*
 * MoleculesToTriangles/CXXClasses/SurfacePrimitive.cpp
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

#include "SurfacePrimitive.h"
#include "Renderer.h"
#include "MoleculesToTriangles/CXXSurface/CXXSurface.h"
#include "MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.h"

SurfacePrimitive::SurfacePrimitive() : VertexColorNormalPrimitive (){
    drawModeGL = DrawAsTriangles;
    enableColorGL = true;
    cxxSurfaceMaker = 0;
    atom2Array = 0;
    atomWeightArray = 0;
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
    // Value-initialised, so a vertex whose atom cannot be found reads back as null rather than
    // as whatever was on the heap.
    atomArray = new const mmdb::Atom*[_nVertices]();
    // The second atom of a saddle, null for a convex cap where there is only one. The weight
    // then stays at 1 and the vertex belongs to its own atom outright.
    atom2Array = new const mmdb::Atom*[_nVertices]();
    atomWeightArray = new float[_nVertices];
    for (size_t i=0; i<_nVertices; i++) atomWeightArray[i] = 1.0f;

    size_t iGlobalVert=0;
    for (vector<CXXSurface>::iterator surfIter = containerSurface.getChildSurfaces().begin();
         surfIter != containerSurface.getChildSurfaces().end();
         ++surfIter){
        CXXSurface &childSurface = *surfIter;
        // Looked up once per child surface rather than per vertex, and SIZE_MAX when this
        // surface has no saddles at all - an isolated atom, say.
        const size_t weightHandle = childSurface.getReadScalarHandle("atomWeight");
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
            // CXXSurface::getPointer returns 0 on success, as the two loops over "atom" in the
            // constructor below both assume. The test here was the other way round, so the
            // pointer was stored only when the lookup had failed - and every vertex whose atom
            // WAS found kept uninitialised heap memory instead. Nothing read atomArray, so it
            // went unnoticed.
            mmdb::Atom *atomPointer;
            int result = childSurface.getPointer("atom", iLocalVert, (void**)&(atomPointer));
            if (result == 0) atomArray[iGlobalVert] = atomPointer;

            // ...and the other atom of the saddle, where this vertex is on one. Only the
            // toroidal patches set these, so a cap leaves the null and the weight of 1 that
            // were put there above. getScalar is asked for the weight only when there really
            // is a second atom, because an unset scalar reads back as zero rather than as one.
            mmdb::Atom *atom2Pointer = 0;
            if (childSurface.getPointer("atom2", iLocalVert, (void**)&(atom2Pointer)) == 0
                && atom2Pointer != 0) {
                atom2Array[iGlobalVert] = atom2Pointer;
                if (weightHandle != SIZE_MAX) {
                    double weight = 1.;
                    if (childSurface.getScalar(int(weightHandle), int(iLocalVert), weight) == 0)
                        atomWeightArray[iGlobalVert] = float(weight);
                }
            }
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
        for (std::size_t i=0; i< childSurface.numberOfTriangles(); i++){
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
    atom2Array = 0;
    atomWeightArray = 0;
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
            for (std::size_t iVertex = 0; iVertex<mySurface.numberOfVertices(); iVertex++){
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
            for (std::size_t iVertex = 0; iVertex<mySurface.numberOfVertices(); iVertex++){
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


