/*
 * MoleculesToTriangles/CXXClasses/RepresentationInstance.h
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

#ifndef RepresentationInstance_H
#define RepresentationInstance_H
#include <memory>

// #define GLM_ENABLE_EXPERIMENTAL // # for norm things (needed?)
// #include <glm/ext.hpp>

#include "RotatedTranslatedScaledEntity.h"

class MolecularRepresentation;
class Renderer;
class Representation;

class RepresentationInstance : public RotatedTranslatedScaledEntity {
    friend class MolecularRepresentationInstance;
    bool doDraw;
    std::shared_ptr<Representation> representation;
    std::shared_ptr<RotatedTranslatedScaledEntity> postTransformation;
public:
    RepresentationInstance() : doDraw(true), representation(0), postTransformation(0) {
        rotation = FCXXCoord (0., 1., 0., 0.);
        translation = FCXXCoord (0., 0., 0., 0.);
        scale = 1.0;
    };
    bool getDoDraw() const {
        return doDraw;
    };
    void setDoDraw(const bool yesOrNo);
    void setRepresentation(std::shared_ptr<Representation>holder){
        representation = holder;
    };
    void reset();
    std::shared_ptr<Representation> getRepresentation(){
        return representation;
    };
    std::shared_ptr<RotatedTranslatedScaledEntity> getPostTransformation() {
        return postTransformation;
    };
    void setPostTransformation(std::shared_ptr<RotatedTranslatedScaledEntity> value){
        postTransformation = value;
    };
    void setPostTranslation(const FCXXCoord  &_translation){
        if (!postTransformation) postTransformation=std::shared_ptr<RotatedTranslatedScaledEntity>(new RotatedTranslatedScaledEntity());
        postTransformation->setTranslation(_translation);
    };
    void setPostRotation(const FCXXCoord  &_rotation){
        if (!postTransformation) postTransformation=std::shared_ptr<RotatedTranslatedScaledEntity>(new RotatedTranslatedScaledEntity());
        postTransformation->setRotation(_rotation);
    };

    void renderInRenderer(std::shared_ptr<Renderer> renderer);
    static RepresentationInstance *defaultRepresentationInstance(std::shared_ptr<Representation> ofRepresentation);
};

#endif
