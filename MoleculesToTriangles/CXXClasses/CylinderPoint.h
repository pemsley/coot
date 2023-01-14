//
//  CylinderPoint.h
//  AesopCD_ios
//
//  Created by Martin Noble on 09/02/2011.
//  Copyright 2011 Dept. of Biochemistry, Oxford University. All rights reserved.
//

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

#ifndef CylinderPoint_H
#define CylinderPoint_H
#include "mmdb2/mmdb_manager.h"

class CylinderPoint {
public:
    enum CylinderPointType {
        CylinderPointTypeStart,
        CylinderPointTypeContinue
    };
    FCXXCoord vertex;
    FCXXCoord color;
    FCXXCoord normalOne;
    FCXXCoord normalTwo;
    float radiusOne;
    float radiusTwo;
    int type;
    const mmdb::Atom *atom;
    
    CylinderPoint(FCXXCoord &_vertex, FCXXCoord &_color, FCXXCoord &_normalOne, FCXXCoord &_normalTwo, float _radiusOne, float _radiusTwo, const mmdb::Atom *_atom) :
    vertex(_vertex), color(_color), normalOne(_normalOne), normalTwo(_normalTwo), radiusOne(_radiusOne), radiusTwo(_radiusTwo), type(CylinderPointTypeContinue), atom(_atom) {
    };
    CylinderPoint(FCXXCoord &_vertex, FCXXCoord &_color, FCXXCoord &_normalOne, FCXXCoord &_normalTwo, float _radiusOne, float _radiusTwo, int _type, const mmdb::Atom *_atom) :
    vertex(_vertex), color(_color), normalOne(_normalOne), normalTwo(_normalTwo), radiusOne(_radiusOne), radiusTwo(_radiusTwo), type(_type), atom(_atom) {
    };
};

#endif
