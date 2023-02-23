/*
 *  SceneSetup.cpp
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "SceneSetup.h"

#include "CXXMatrix.h"
#include "Light.h"
#include "MolecularRepresentation.h"
#include "RepresentationInstance.h"
#include "Renderer.h"
#include "Camera.h"


// which header do I want for this?
// #ifndef APIENTRY
// #define APIENTRY
// #endif
// #ifndef GLAPI
// #define GLAPI extern
// #endif
// typedef unsigned int GLenum;
// GLAPI GLenum APIENTRY glGetError (void);

// #include <GL/gl.h>

// #include <epoxy/gl.h>



void SceneSetup::renderWithRendererFromViewpointDontChangeMVM(std::shared_ptr<Renderer> renderer, const FCXXCoord &viewpoint){
#ifdef DEBUG_MINE
    std::cout << "In scenesetup render\n";
#endif
    renderWithRendererFromViewpoint(renderer, viewpoint);
}

void SceneSetup::renderWithRendererFromViewpoint(std::shared_ptr<Renderer> renderer, const FCXXCoord &viewpoint){
#ifdef DEBUG_MINE
    std::cout << "In scenesetup render\n";
#endif
    renderer->setupScene(this);
    GLenum err = glGetError(); std::cout << "Here in renderWithRendererFromViewpoint() 2 err " << err << std::endl;
    auto lightsPntr = lights.begin();
    auto lightsEnd = lights.end();
    int iLight = 0;
    for (; lightsPntr != lightsEnd; ++lightsPntr){
        (*lightsPntr)->setupInRendererAsIndexFromViewpointWithScale(renderer, iLight++, viewpoint, scale);
    }
      
    // Representation instances
    auto representationInstancesPntr = representationInstances.begin();
    auto representationInstancesEnd = representationInstances.end();
    for (; representationInstancesPntr != representationInstancesEnd; ++representationInstancesPntr){
        err = glGetError(); std::cout << "Here in renderWithRendererFromViewpoint() 4 err " << err << std::endl;
        (*representationInstancesPntr)->renderInRenderer(renderer);
        err = glGetError(); std::cout << "Here in renderWithRendererFromViewpoint() 5 err " << err << std::endl;
    }
}

void SceneSetup::addCamera(std::shared_ptr<Camera> camera) {
    //Note that to avoid accumulating multiple references to the same
    //object, I screen here
    if (std::find(cameras.begin(), cameras.end(), camera) == cameras.end()){
        cameras.push_back(camera);
    }
};

std::shared_ptr<SceneSetup> SceneSetup::defaultSceneSetup() {
    std::shared_ptr<SceneSetup> result(new SceneSetup());
    result->setBackgroundColor(FCXXCoord (0., 0., 0., 1.0));
    result->setRotation(FCXXCoord (0.,1.,0.,0.));
    result->setTranslation(FCXXCoord (0.,0.,0.));
    result->setScale(1.0);
    return result;
}

std::set<std::shared_ptr<MolecularRepresentation> > SceneSetup::molecularRepresentations()
{
    std::set<std::shared_ptr<MolecularRepresentation> > result;
    auto inst = representationInstancesBegin();
    auto lastInst = representationInstancesEnd();
    for (; inst != lastInst; ++inst){
        std::shared_ptr<MolecularRepresentation> molRep = std::dynamic_pointer_cast<MolecularRepresentation>((*inst)->getRepresentation());
        if (molRep) result.insert(molRep);
    }
    return result;
}

