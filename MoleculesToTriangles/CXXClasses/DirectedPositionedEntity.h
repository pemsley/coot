/*
 * MoleculesToTriangles/CXXClasses/DirectedPositionedEntity.h
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
