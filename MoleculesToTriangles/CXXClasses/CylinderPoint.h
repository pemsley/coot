/*
 * MoleculesToTriangles/CXXClasses/CylinderPoint.h
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
//  CylinderPoint.h
//  AesopCD_ios
//
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
