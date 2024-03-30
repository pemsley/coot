/*
 * MoleculesToTriangles/CXXClasses/Camera.h
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
