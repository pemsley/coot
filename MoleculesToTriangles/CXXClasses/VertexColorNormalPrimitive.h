/*
 *  VertexColorNormalPrimitive.h
 *  AesopCoreData
 *
 *  Created by Martin Noble on 24/03/2010.
 *  Copyright 2010 Apple Inc. All rights reserved.
 *
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
