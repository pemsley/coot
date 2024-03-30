/*
 * MoleculesToTriangles/CXXClasses/BondsPrimitive.h
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
