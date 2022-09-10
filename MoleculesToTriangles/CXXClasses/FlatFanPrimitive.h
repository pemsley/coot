//
//  FlatFanPrimitive.h
//  MoleculesToTriangles
//
//  Created by Martin Noble on 16/03/2017.
//  Copyright © 2017 MartinNobleSoftware. All rights reserved.
//

#ifndef FlatFanPrimitive_h
#define FlatFanPrimitive_h

#include <stdio.h>
#include "VertexColorNormalPrimitive.h"
#include "mmdb2/mmdb_manager.h"

class FlatFanPrimitive : public VertexColorNormalPrimitive
{
private:
    std::vector<mmdb::Atom*> atoms;
    FCXXCoord color;
public:
    FlatFanPrimitive(const std::vector<mmdb::Atom*> &atoms_in, FCXXCoord &color_in) {
        atoms = atoms_in;
        color = color_in;
    };
    virtual ~FlatFanPrimitive(){
        emptyArrays();
        //std::cout << "In FlatFanPrimitive destructor" << std::endl;
    };
    virtual void generateArrays();
    
};


#endif /* FlatFanPrimitive_h */
