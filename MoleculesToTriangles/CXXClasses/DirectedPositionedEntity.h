/*
 *  DirectedPositionedEntity.h
 *  Aesop
 *
 *  Created by Martin Noble on 14/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DirectedPositionedEntity_H
#define DirectedPositionedEntity_H

#include "PositionedEntity.h"

class DirectedPositionedEntity : public PositionedEntity {
protected:
    FCXXCoord direction;
public:
    const FCXXCoord &getDirection() const{
        return direction;
    };
    FCXXCoord &getDirection() {
        return direction;
    };
    float getDirectionX() const {
        return direction[0];
    };
    float getDirectionY() const {
        return direction[1];
    };
    float getDirectionZ() const {
        return direction[2];
    };
    void setDirectionX(const float _x){
        direction[0] = _x;
    };
    void setDirectionY(const float _y){
        direction[1] = _y;
    };
    void setDirectionZ(const float _z){
        direction[2] = _z;
    };
    void setDirection(const FCXXCoord &_direction){
        direction = _direction;
    };
};

#endif
