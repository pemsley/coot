/*
 * MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.h
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
#ifndef VertexColorNormalPrimitive_h
#define VertexColorNormalPrimitive_h

#include "DisplayPrimitive.h"
#include "GLIndexType.h"
#include "mmdb2/mmdb_manager.h"

class Renderer;

class VertexColorNormalPrimitive : public DisplayPrimitive {
public:
    typedef struct VertexColorNormal_ {
        float vertex[4];
        float color[4];
        float normal[4];
    } VertexColorNormal;
    typedef struct VertexNormal_ {
        float vertex[4];
        float normal[4];
    } VertexNormal;
    typedef struct VertexColor_ {
        float vertex[4];
        float color[4];
    } VertexColor;
protected:
    unsigned long _nTriangles;
    GLIndexType *indexArray;
    unsigned long _nVertices;
    unsigned long _nLines;
    VertexColorNormal *vertexColorNormalArray;
    VertexNormal *vertexNormalArray;
    VertexColor *vertexColorArray;
    const mmdb::Atom **atomArray;
    
    //Here store whether this is to be drawn with glDrawElements using mode GL_TRIANGLES or GL_TRIANGLE_STRIP
    int drawModeGL;
    bool enableColorGL;
public:
    enum GL_DRAW_MODES {
        DrawAsTriangles, DrawAsTriangleStrip, DrawAsLines
    };
    VertexColorNormalPrimitive (){
        vertexColorNormalArray = 0;
        vertexNormalArray = 0;
        vertexColorArray = 0;
        indexArray = 0;
        atomArray = 0;
    }
    virtual ~VertexColorNormalPrimitive(){
        emptyArrays();
    };
    void emptyArrays() {
        delete [] vertexColorNormalArray;
        vertexColorNormalArray = 0;
        delete [] vertexNormalArray;
        vertexNormalArray = 0;
        delete [] vertexColorArray;
        vertexColorArray = 0;
        delete [] indexArray;
        indexArray = 0;
        delete [] atomArray;
        atomArray = 0;
    };
    
    //This class implements the render method of the base DisplayPrimitiveClass
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
    
    const VertexColorNormal *getVertexColorNormalArray() const {
        return vertexColorNormalArray;
    };
    
    VertexColorNormal *getVertexColorNormalArray() {
        return vertexColorNormalArray;
    };
    
    const VertexNormal *getVertexNormalArray() const {
        return vertexNormalArray;
    };
    
    const VertexColor *getVertexColorArray() const {
        return vertexColorArray;
    };
    
    const GLIndexType *getIndexArray() const {
        return indexArray;
    };
    const mmdb::Atom **getAtomArray() const {
        return atomArray;
    };
    unsigned long nTriangles(){
        return _nTriangles;
    };
    unsigned long nVertices(){
        return _nVertices;
    };
    unsigned long nLines(){
        return _nLines;
    };
};

#endif
