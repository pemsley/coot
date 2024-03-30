/*
 *  SceneSetup.h
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SceneSetup_H
#define SceneSetup_H

#include "RotatedTranslatedScaledEntity.h"
#include <list>
#include <set>
#include <algorithm>
#include <memory>

class Light;
class RepresentationInstance;
class MolecularRepresentationInstance;
class MolecularRepresentation;
class Renderer;
class Camera;

class SceneSetup : public RotatedTranslatedScaledEntity {
private:
    FCXXCoord backgroundColor;
    std::list<std::shared_ptr<Light> >lights;
    std::list<std::shared_ptr<RepresentationInstance> >representationInstances;
    std::list<std::shared_ptr<Camera> >cameras;
public:
    const FCXXCoord &getBackgroundColor() const {
        return backgroundColor;
    };
    void setBackgroundColor(const FCXXCoord &color){
        backgroundColor = color;
    };
    void renderWithRendererFromViewpoint(std::shared_ptr<Renderer> renderer, const FCXXCoord &viewpoint);
    void renderWithRendererFromViewpointDontChangeMVM(std::shared_ptr<Renderer> renderer, const FCXXCoord &viewpoint);
    void addLight(std::shared_ptr<Light> light){
        //Note that to avoid accumulating multiple references to the same
        //object, I screen here
        if (std::find(lights.begin(), lights.end(), light) == lights.end()){
            lights.push_back(light);
        }
    };
    void removeLight(std::shared_ptr<Light> light){
        lights.remove(light);
    };
    void clearRepresentatationInstances(){
        representationInstances.clear();
    };
    void addCamera(std::shared_ptr<Camera> camera);
    void removeCamera(std::shared_ptr<Camera> camera) {
        cameras.remove(camera);
    };
    void addRepresentationInstance(std::shared_ptr<RepresentationInstance> representationInstance) {
        if (std::find(representationInstances.begin(), representationInstances.end(), representationInstance) ==representationInstances.end()){
            representationInstances.push_back(representationInstance);
        }
    };
    void removeRepresentationInstance(std::shared_ptr<RepresentationInstance> representationInstance) {
        if (std::find(representationInstances.begin(), representationInstances.end(), representationInstance) !=representationInstances.end()){
            representationInstances.remove(representationInstance);
        }
    };
    std::list<std::shared_ptr<Light> >::iterator lightsBegin() {
        return lights.begin();
    };
    std::list<std::shared_ptr<Light> >::iterator lightsEnd() {
        return lights.end();
    };
    std::shared_ptr<Light> getLight(int iLight) {
        auto returnPntr = lightsBegin();
        for (int i=0; i<iLight && returnPntr != lightsEnd(); i++) returnPntr++;
        if (returnPntr != lightsEnd()) return *returnPntr;
        return 0;
    };
    std::shared_ptr<RepresentationInstance> getRepresentationInstance(int iRepresentationInstance) {
        auto returnPntr = representationInstancesBegin();
        for (int i=0; i<iRepresentationInstance && returnPntr != representationInstancesEnd(); i++) returnPntr++;
        if (returnPntr != representationInstancesEnd()) return *returnPntr;
        return 0;
    };
    std::list<std::shared_ptr<RepresentationInstance> >::iterator representationInstancesBegin() {
        return representationInstances.begin();
    };
    std::list<std::shared_ptr<RepresentationInstance> >::iterator representationInstancesEnd() {
        return representationInstances.end();
    };
    
    std::set<std::shared_ptr<MolecularRepresentation> > molecularRepresentations();
    
    static std::shared_ptr<SceneSetup> defaultSceneSetup();
    
};

#endif
