//
//  SticksPrimitive.h
//  AesopCD_ios
//
//  Created by Martin Noble on 14/02/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
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
