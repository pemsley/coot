/*
 *  PositionedEntity.h
 *  Aesop
 *
 *  Created by Martin Noble on 14/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PositionedEntity_H
#define PositionedEntity_H

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class PositionedEntity {
private:
    FCXXCoord position;
public:
    float getPositionX() const{
        return position[0];
    };
    float getPositionY() const{
        return position[1];
    };
    float getPositionZ() const{
        return position[2];
    };
    void setPositionX(float _x) {
        position[0] = _x;
    };
    void setPositionY(float _y) {
        position[1] = _y;
    };
    void setPositionZ(float _z) {
        position[2] = _z;
    };
    void setPosition(const FCXXCoord &_position){
        position = _position;
    };
    const FCXXCoord &getPosition() const{
        return position;
    };
    FCXXCoord &getPosition() {
        return position;
    };
};
    
#endif
