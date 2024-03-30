/*
 *  RepresentationInstance.mm
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "RepresentationInstance.h"

#include "Renderer.h"
#include "Representation.h"

void RepresentationInstance::renderInRenderer(std::shared_ptr<Renderer> renderer){
#ifdef DEBUG_MINE
    std::cout << "In representation instance render "<<getDoDraw()<<"\n";
#endif
    if (getDoDraw()) {
        renderer->setupRepresentationInstance(this);
        representation->renderWithRenderer(renderer);
        renderer->restoreModelviewMatrix();
    }
}

void RepresentationInstance::reset()
{
    setRotation(FCXXCoord (0.,1.,0.,0.));
    setTranslation(FCXXCoord (0.,0.,0.,0.));
    setScale(1.);
}

RepresentationInstance *RepresentationInstance::defaultRepresentationInstance(std::shared_ptr<Representation>ofRepresentation){
    RepresentationInstance *inst = new RepresentationInstance();
    inst->setRepresentation(ofRepresentation);
    return inst;
}

void RepresentationInstance::setDoDraw(const bool yesOrNo){
    doDraw = yesOrNo;
    if (yesOrNo && representation) representation->setDoDraw(yesOrNo);
};
