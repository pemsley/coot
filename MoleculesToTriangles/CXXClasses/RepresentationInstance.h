/*
 *  RepresentationInstance.h
 *  Aesop
 *
 *  Created by Martin Noble on 02/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RepresentationInstance_H
#define RepresentationInstance_H
#include <memory>

#define GLM_ENABLE_EXPERIMENTAL // # for norm things (needed?)
#include <glm/ext.hpp>

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
