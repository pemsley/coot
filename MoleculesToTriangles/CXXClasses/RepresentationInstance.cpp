/*
 * MoleculesToTriangles/CXXClasses/RepresentationInstance.cpp
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
