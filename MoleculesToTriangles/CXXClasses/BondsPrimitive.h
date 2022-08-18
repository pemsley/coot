/*
 *  BondsPrimitive.h
 *  AesopCD2.0
 *
 *  Created by Martin Noble on 09/07/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef BondsPrimitive_H
#define BondsPrimitive_H

#include "DisplayPrimitive.h"
#include <vector>
#include <map>
#include <memory>
#include "GLIndexType.h"
#include "mmdb2/mmdb_manager.h"

class ColorScheme;
class ColorRule;
class Renderer;

class BondsPrimitive : public DisplayPrimitive {
public:
    typedef struct _VertexColor {
        float vertex[4];
        float color[4];
    } VertexColor;
private:
    VertexColor *vertexColorArray;
    GLIndexType *indexArray;
    std::map<mmdb::Atom *, std::vector<mmdb::Atom *> > bonds;
    int nBonds;
    std::shared_ptr<ColorScheme> colorScheme;
public:
    BondsPrimitive () : vertexColorArray(0), indexArray(0), nBonds(0){};
    ~BondsPrimitive () {
        invalidateGLPrimitives();
    };
    void invalidateGLPrimitives();
    void evaluateGLPrimitives(std::map<std::shared_ptr<ColorRule>,int> &handles);
    void setColorScheme(std::shared_ptr<ColorScheme> _colorScheme){
        colorScheme = _colorScheme;
    };
    void addPair(mmdb::Atom *atom1, mmdb::Atom *atom2) {
        std::vector<mmdb::Atom *> &atom1List(bonds[atom1]);
        atom1List.push_back(atom2);
        nBonds++;
    };
    int getNBonds(){
        return nBonds;
    };
    VertexColor *getVertexColorArray(){
        return vertexColorArray;
    };
    GLIndexType *getIndexArray(){
        return indexArray;
    };
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
    
    //Children of this class will have to implement the following
    virtual void generateArrays(){};
    

};

#endif
