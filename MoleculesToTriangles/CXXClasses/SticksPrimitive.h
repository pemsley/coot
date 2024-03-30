/*
 * MoleculesToTriangles/CXXClasses/SticksPrimitive.h
 *
 * Copyright 2011 by Martin Noble, University of Oxford
 * Author: Martin Noble
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

//
//  SticksPrimitive.h
//  AesopCD_ios
//

#ifndef SticksPrimitive_H
#define SticksPrimitive_H

#include "VertexColorNormalPrimitive.h"
#include <vector>
#include <memory>
#include <map>
#include "GLIndexType.h"
#include "mmdb2/mmdb_manager.h"

class Renderer;
class ColorScheme;

class SticksPrimitive : public VertexColorNormalPrimitive {
public:
private:
    std::map<mmdb::Atom *, std::vector<mmdb::Atom *> > sticks;
    std::shared_ptr<ColorScheme> colorScheme;
    mmdb::Manager *mmdb;
public:
    void addPair(mmdb::Atom *atom1, mmdb::Atom *atom2) {
        std::vector<mmdb::Atom *> &atom1List(sticks[atom1]);
        atom1List.push_back(atom2);
        _nLines++;
    };
    void setColorScheme(std::shared_ptr<ColorScheme> aScheme){
        colorScheme = aScheme;
    };
    void setMmdb(mmdb::Manager *aMmdb){
        mmdb = aMmdb;
    };
    SticksPrimitive() {
        _nLines = 0;
        vertexColorNormalArray = 0;
        indexArray = 0;
        emptyArrays();
        drawModeGL = DrawAsLines;
        enableColorGL = true;
    };
    void emptyArrays(){
        _nLines = 0;
        delete [] vertexColorArray;
        vertexColorArray = 0;
        delete [] indexArray;
        indexArray = 0;
    };
    virtual ~SticksPrimitive(){
        emptyArrays();
        //std::cout << "In Sticks destructor" << std::endl;
    };
    virtual void generateArrays();
    
    virtual void renderWithRenderer(std::shared_ptr<Renderer> renderer);
    
};

#endif
