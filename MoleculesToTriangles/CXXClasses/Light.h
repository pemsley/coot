/*
 *  Light.h
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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
