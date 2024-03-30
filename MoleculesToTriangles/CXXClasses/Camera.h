/*
 *  Camera.h
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CAMERA_H
#define CAMERA_H
#include <memory>

#include "RotatedTranslatedScaledEntity.h"

class SceneSetup;
class Renderer;

class Camera : public RotatedTranslatedScaledEntity { 
protected:
    bool active;
    std::shared_ptr<SceneSetup> sceneSetup;
    float zClipFront;
    float zClipDepth;
    float fogFront;
    float fogDepthRange;
    float fovy;
public:
    bool getActive() const {
        return active;
    };
    void setActive(const bool &_active) {
        active = _active;
    };
    void setSceneSetup(std::shared_ptr<SceneSetup> _sceneSetup){
        sceneSetup = _sceneSetup;
    };
    std::shared_ptr<SceneSetup> getSceneSetup(){
        return sceneSetup;
    };
    float getFovy() const {
        return fovy;
    };
    float getZClipFront() const {
        return zClipFront;
    };
    void setFovy(const float &_fovy){
        fovy = _fovy;
    };
    void setZClipFront(const float &_zClipFront){
        zClipFront = _zClipFront;
    };
    float getZClipDepth() const {
        return zClipDepth;
    };
    void setZClipDepth(const float &_zClipDepth){
        zClipDepth = _zClipDepth;
    };
    float getFogFront() const {
        return fogFront;
    };    
    void setFogFront(const float _fogFront){
        fogFront = _fogFront;
    };
    float getFogDepthRange() const {
        return fogDepthRange;
    };
    void setFogDepthRange(const float _fogDepthRange){
        fogDepthRange = _fogDepthRange;
    };
    static std::shared_ptr<Camera> defaultCamera() {
        std::shared_ptr<Camera>result(new Camera());
        result->setRotation(FCXXCoord (0., 1., 0., 0.));
        result->setTranslation(FCXXCoord (0., 0., 1000.));
        result->setSceneSetup(0);
        result->setZClipFront(100.);
        result->setZClipDepth(200.);
        result->setFogDepthRange(40);
        result->setFogFront(-20);
        result->setFovy(11);
        return result;
    };
    void renderWithRenderer(std::shared_ptr<Renderer> renderer);
};
#endif
