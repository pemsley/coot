/*
 * MoleculesToTriangles/CXXClasses/SceneSetup.cpp
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

