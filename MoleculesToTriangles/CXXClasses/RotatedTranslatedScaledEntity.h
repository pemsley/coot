/*
 *  RotatedTranslatedScaledEntity.h
 *  Aesop
 *
 *  Created by Martin Noble on 14/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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
