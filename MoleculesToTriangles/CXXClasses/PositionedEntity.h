/*
 * MoleculesToTriangles/CXXClasses/PositionedEntity.h
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
