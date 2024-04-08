/*
 * MoleculesToTriangles/CXXClasses/Light.h
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

#ifndef Light_H
#define Light_H
#include <memory>

#include "RotatedTranslatedScaledEntity.h"
class Renderer;
class SceneSetup;

class Light : public RotatedTranslatedScaledEntity {
protected:
    FCXXCoord  ambient;
    FCXXCoord  diffuse;
    FCXXCoord  specular;
    float intensity;
    float exponent;
    SceneSetup *sceneSetup;
    int lightType;
    bool drawLight;
public:
    ~Light(){
#ifdef DEBUG_MINE
        std::cout<<"In light destructor";
#endif
    };
    enum LightType {
        Directional, Positional, SpotLight
    };
    
    FCXXCoord  getAmbient() const{
        return ambient;
    };
    FCXXCoord  getDiffuse() const {
        return diffuse;
    };
    FCXXCoord  getSpecular() const {
        return specular;
    };
    float getIntensity() const {
        return intensity;
    };
    float getExponent() const {
        return exponent;
    };
    SceneSetup *getSceneSetup() {
        return sceneSetup;
    };
    void setSceneSetup(SceneSetup *_sceneSetup){
        sceneSetup = _sceneSetup;
    };
    int getLightType() const {
        return lightType;
    };
    void setAmbient(const FCXXCoord &_ambient){
        ambient = _ambient;
    };
    void setDiffuse(const FCXXCoord &_diffuse){
        diffuse = _diffuse;
    };
    void setSpecular(const FCXXCoord &_specular){
        specular = _specular;
    };
    void setIntensity(const float &_intensity){
        intensity = _intensity;
    };
    void setExponent(const float &_exponent){
        exponent = _exponent;
    };
    void setLightType(int _lightType){
        lightType = _lightType;
    };
    void setupInRendererAsIndexFromViewpointWithScale (std::shared_ptr<Renderer> renderer, int lightEnum, FCXXCoord  fromViewpoint, float scale);
    void setDrawLight(bool yesOrNo){
        drawLight = yesOrNo;
    };
    bool getDrawLight() const {
        return drawLight;
    };
    static std::shared_ptr<Light> defaultLight();
};
#endif
