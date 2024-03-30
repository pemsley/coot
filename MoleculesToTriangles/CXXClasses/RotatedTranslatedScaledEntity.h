/*
 * MoleculesToTriangles/CXXClasses/RotatedTranslatedScaledEntity.h
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

#ifndef RotatedTranslatedScaledEntity_H
#define RotatedTranslatedScaledEntity_H

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class RotatedTranslatedScaledEntity {
protected:
    FCXXCoord  rotation;
    FCXXCoord  translation;
    float scale;
public:
    RotatedTranslatedScaledEntity () {
        rotation = FCXXCoord  (0., 1., 0., 0.);
        translation = FCXXCoord  (0., 0., 0., 0.);
        scale = 1.;
    };
    FCXXCoord  getRotation() {
        return rotation;
    };
    FCXXCoord  getTranslation() {
        return translation;
    };
    float getScale() const {
        return scale;
    };
    void setRotation(const FCXXCoord  &_rotation){
        rotation = _rotation;
    };
    void setTranslation(const FCXXCoord  &_translation){
        translation = _translation;
    };
    void setScale(const float &_scale){
        scale = _scale;
    };
    
    float getRotationV() const {
        return rotation[0];
    };
    float getRotationX() const {
        return rotation[1];
    };
    float getRotationY() const {
        return rotation[2];
    };
    float getRotationZ() const {
        return rotation[3];
    };
    float getTranslationX() const {
        return translation[0];
    };
    float getTranslationY() const {
        return translation[1];
    };
    float getTranslationZ() const {
        return translation[2];
    };
    
    void setRotationV(const float &v){
        rotation[0] = v;
    };
    void setRotationX(const float &x){
        rotation[1] = x;
    };
    void setRotationY(const float &y){
        rotation[2] = y;
    };
    void setRotationZ(const float &z){
        rotation[3] = z;
    };
    
    void setTranslationX(const float &x){
        translation[0] = x;
    };
    void setTranslationY(const float &y){
        translation[1] = y;
    };
    void setTranslationZ(const float &z){
        translation[2] = z;
    };
    void translateBy(const FCXXCoord  &vector){
        translation += vector;
    };
    void rotateBy(const FCXXCoord  &vector);
};
    
#endif
